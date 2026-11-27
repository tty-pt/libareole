/* ========================================================================= */
/* Chunk Population, Parent/Child Vertex Reuse & Async Staging               */
/* ========================================================================= */

#ifdef TTYPT_WQ_H
static void ar_async_populate_worker(void *ctx)
{
	ar_node_meta_t *meta = (ar_node_meta_t *)ctx;
	ar_system_t *sys = meta ? meta->async_sys : NULL;
	if (!sys || !meta->staging_samples || !meta->staging_values || !meta->staging_rgba)
		return;

	int32_t gdim = (int32_t)sys->grid_dim;
	if (gdim < 1) gdim = 1;
	int32_t vdim = gdim + 1;
	uint16_t world = meta->async_world;
	uint8_t level = meta->async_level;
	int32_t ox = meta->async_ox;
	int32_t oy = meta->async_oy;
	int32_t w = meta->async_w;
	int32_t h = meta->async_h;

	int32_t min_val = INT32_MAX;
	int32_t max_val = INT32_MIN;
	uint32_t flags = 0;
	if (sys->sparsity) {
		if (sys->sparsity(sys->user_ctx, world, level, ox, oy, w, h))
			flags |= AR_NODE_FLAG_SPARSE;
	}

	size_t idx = 0;
	if (sys->row_sampler) {
		for (int r = 0; r <= gdim; r++) {
			int32_t wy = oy + (int32_t)(((int64_t)r * h) / gdim);
			sys->row_sampler(sys->user_ctx, world, level, wy, ox, w, (uint32_t)vdim, &meta->staging_samples[idx]);
			for (int c = 0; c <= gdim; c++) {
				ar_sample_t *s = &meta->staging_samples[idx];
				if (s->val < min_val) min_val = s->val;
				if (s->val > max_val) max_val = s->val;
				flags |= (s->flags & AR_NODE_USER_FLAGS_MASK);
				meta->staging_values[idx] = s->val;
				meta->staging_rgba[idx] = s->rgba ? s->rgba : sys->default_rgba;
				idx++;
			}
		}
	} else {
		for (int r = 0; r <= gdim; r++) {
			int32_t wy = oy + (int32_t)(((int64_t)r * h) / gdim);
			for (int c = 0; c <= gdim; c++) {
				int32_t wx = ox + (int32_t)(((int64_t)c * w) / gdim);
				ar_sample_t *s = &meta->staging_samples[idx];
				sys->sampler(sys->user_ctx, world, level, wx, wy, s);
				if (s->val < min_val) min_val = s->val;
				if (s->val > max_val) max_val = s->val;
				flags |= (s->flags & AR_NODE_USER_FLAGS_MASK);
				meta->staging_values[idx] = s->val;
				meta->staging_rgba[idx] = s->rgba ? s->rgba : sys->default_rgba;
				idx++;
			}
		}
	}

	meta->async_min_val = min_val;
	meta->async_max_val = max_val;
	meta->async_flags = flags;
}
#endif

static ar_node_t *ar_get_or_populate(ar_system_t *sys, uint16_t world,
                                     uint8_t level, int32_t ox, int32_t oy)
{
	uint64_t key = morton_quad_key(world, level, ox, oy, sys->base_w, sys->base_h);

	/* 1. Fast O(1) Islet spatial map lookup */
	uint32_t hit_slot = 0;
	int slot = -1;
	bool upgrading = false;
	if (ar_map_get(sys, key, &hit_slot)) {
		if (hit_slot < sys->pool_capacity && sys->pool_meta[hit_slot].state == 1) {
			ar_node_t *node = &sys->pool[hit_slot];
			sys->pool_meta[hit_slot].last_used_tick = sys->current_tick;
			sys->pool_meta[hit_slot].last_touch_frame = sys->current_tick;
			ar_lru_touch(sys, (int32_t)hit_slot);
			sys->stats.cache_hits_total++;
			sys->metrics.frame_hits++;
#ifdef TTYPT_WQ_H
			if ((node->flags & AR_NODE_FLAG_DEFERRED) &&
			    sys->wq && sys->staging_arena) {
				ar_node_meta_t *hmeta = &sys->pool_meta[hit_slot];
				uint8_t st = atomic_load_explicit(&hmeta->async_task.state, memory_order_acquire);
				if (st == STL_WQ_QUEUED || st == STL_WQ_RUNNING || st == STL_WQ_DONE)
					return node;
				if (sys->max_populates_per_frame > 0 && sys->frame_populates_count >= sys->max_populates_per_frame) {
					hmeta->async_sys = sys;
					hmeta->async_slot = hit_slot;
					hmeta->async_key = key;
					hmeta->async_world = world;
					hmeta->async_level = level;
					hmeta->async_ox = node->ox;
					hmeta->async_oy = node->oy;
					hmeta->async_w = node->w;
					hmeta->async_h = node->h;
					(void)ar_wq_submit(sys->wq, &hmeta->async_task, ar_async_populate_worker, hmeta);
					return node;
				}
			}
#endif
			if ((node->flags & AR_NODE_FLAG_DEFERRED) &&
			    (sys->max_populates_per_frame == 0 || sys->frame_populates_count < sys->max_populates_per_frame)) {
				slot = (int)hit_slot;
				upgrading = true;
			} else {
				return node;
			}
		}
	}

	uint64_t t_pop = ar_now_us(sys->timing_enabled);
	sys->stats.cache_misses_total++;
	sys->metrics.frame_misses++;

	/* 2. Allocate free slot or evict O(1) from LRU tail */
	if (!upgrading) {
		slot = ar_pool_alloc_or_evict(sys);
		if (slot < 0)
			return NULL;
		ar_lru_insert_head(sys, (int32_t)slot);
	}

	ar_node_t *node = &sys->pool[slot];
	ar_node_meta_t *meta = &sys->pool_meta[slot];
#ifdef TTYPT_WQ_H
	if (sys->wq && stl_wq_busy(&meta->async_task)) {
		ar_wq_wait(sys->wq, &meta->async_task);
	}
	meta->async_task.state = STL_WQ_IDLE;
#endif
	node->world = world;
	node->level = level;
	node->ox = ox;
	node->oy = oy;
	node->w = sys->base_w << level;
	node->h = sys->base_h << level;
	node->flags = 0;
	node->min_val = INT32_MAX;
	node->max_val = INT32_MIN;
	meta->last_used_tick = sys->current_tick;
	meta->last_touch_frame = sys->current_tick;
	if (!upgrading)
		meta->last_active_tick = 0;
	meta->active_view_mask = 0;
	meta->prev_active_view_mask = 0;
	meta->state = 1;
	meta->morton_key = key;

	if (!upgrading) {
		islet_quad_map_put(&sys->map, key, (uint32_t)slot);
	} else {
		sys->metrics.frame_deferred_upgrades++;
	}

	ar_node_bind_samples(sys, slot, level);

	/* 3. Query sparsity if provider supplied */
	if (sys->sparsity) {
		if (sys->sparsity(sys->user_ctx, world, level, ox, oy, node->w, node->h)) {
			node->flags |= AR_NODE_FLAG_SPARSE;
		}
	}

	/* 4. Check cached parent and child nodes for vertex reuse */
	int32_t gdim = (int32_t)sys->grid_dim;
	if (gdim < 1) gdim = 1;
	int32_t vdim = gdim + 1;
	int32_t half_gdim = gdim >> 1;
	bool can_reuse_hierarchy = ((gdim & 1) == 0 && half_gdim > 0);

	const ar_node_t *parent_node = NULL;
	int parent_r_off = 0, parent_c_off = 0;
	const ar_node_t *child_nodes[4] = { NULL, NULL, NULL, NULL };

	if (can_reuse_hierarchy) {
		uint64_t pkey = morton_quad_parent(key);
		if (pkey != key) {
			uint32_t pslot = 0;
			if (ar_map_get(sys, pkey, &pslot) &&
			    pslot < sys->pool_capacity && sys->pool_meta[pslot].state == 1 &&
			    sys->pool[pslot].world == world && sys->pool[pslot].level == level + 1 &&
			    (sys->pool[pslot].flags & AR_NODE_FLAG_DEFERRED) == 0) {
				parent_node = &sys->pool[pslot];
				uint8_t q = morton_quad_quadrant(key);
				parent_c_off = (q & 1) ? half_gdim : 0;
				parent_r_off = (q & 2) ? half_gdim : 0;
			}
		}
		if (level > 0) {
			for (uint8_t q = 0; q < 4; q++) {
				uint64_t ckey = morton_quad_child(key, q);
				uint32_t cslot = 0;
				if (ar_map_get(sys, ckey, &cslot) &&
				    cslot < sys->pool_capacity && sys->pool_meta[cslot].state == 1 &&
				    sys->pool[cslot].world == world && sys->pool[cslot].level == level - 1 &&
				    (sys->pool[cslot].flags & AR_NODE_FLAG_DEFERRED) == 0) {
					child_nodes[q] = &sys->pool[cslot];
				}
			}
		}
	}

	int32_t     *mut_values  = (int32_t *)node->values;
	uint32_t    *mut_rgba    = (uint32_t *)node->rgba;
	ar_sample_t *mut_samples = (ar_sample_t *)node->samples;

	/* Deferred bilinear upsampling fallback if per-frame populate budget exceeded */
	bool budget_exceeded = (sys->max_populates_per_frame > 0 &&
	                        sys->frame_populates_count >= sys->max_populates_per_frame);
#ifdef TTYPT_WQ_H
	bool has_async_wq = (sys->wq != NULL && sys->staging_arena != NULL);
#else
	bool has_async_wq = false;
#endif
	if (budget_exceeded && parent_node && (!mut_samples || parent_node->samples != NULL || has_async_wq)) {
		size_t idx = 0;
		int vdim_p = (int32_t)parent_node->vert_dim;
		const ar_sample_t *ps = parent_node->samples;
		const int32_t *pv = parent_node->values;
		const uint32_t *pc_arr = parent_node->rgba;
		for (int r = 0; r <= gdim; r++) {
			int r_even = ((r & 1) == 0);
			int pr0 = parent_r_off + (r >> 1);
			int pr1 = r_even ? pr0 : (pr0 + 1);
			if (pr1 > gdim) pr1 = gdim;
			for (int c = 0; c <= gdim; c++) {
				int c_even = ((c & 1) == 0);
				int pc0 = parent_c_off + (c >> 1);
				int pc1 = c_even ? pc0 : (pc0 + 1);
				if (pc1 > gdim) pc1 = gdim;

				int32_t v_out;
				uint32_t rgba_out;
				if (mut_samples && ps) {
					ar_sample_t *s = &mut_samples[idx];
					if (r_even && c_even) {
						*s = ps[pr0 * vdim_p + pc0];
					} else if (r_even) {
						const ar_sample_t *s0 = &ps[pr0 * vdim_p + pc0];
						const ar_sample_t *s1 = &ps[pr0 * vdim_p + pc1];
						s->val = (int32_t)(((int64_t)s0->val + (int64_t)s1->val) >> 1);
						s->flags = s0->flags | s1->flags;
						s->rgba = (((s0->rgba & 0xFEFEFEFEu) >> 1) + ((s1->rgba & 0xFEFEFEFEu) >> 1));
						for (uint8_t a = 0; a < AR_SAMPLE_ATTR_CAP; a++) {
							s->attr[a] = (uint8_t)(((uint16_t)s0->attr[a] + (uint16_t)s1->attr[a]) >> 1);
						}
					} else if (c_even) {
						const ar_sample_t *s0 = &ps[pr0 * vdim_p + pc0];
						const ar_sample_t *s1 = &ps[pr1 * vdim_p + pc0];
						s->val = (int32_t)(((int64_t)s0->val + (int64_t)s1->val) >> 1);
						s->flags = s0->flags | s1->flags;
						s->rgba = (((s0->rgba & 0xFEFEFEFEu) >> 1) + ((s1->rgba & 0xFEFEFEFEu) >> 1));
						for (uint8_t a = 0; a < AR_SAMPLE_ATTR_CAP; a++) {
							s->attr[a] = (uint8_t)(((uint16_t)s0->attr[a] + (uint16_t)s1->attr[a]) >> 1);
						}
					} else {
						const ar_sample_t *s00 = &ps[pr0 * vdim_p + pc0];
						const ar_sample_t *s10 = &ps[pr0 * vdim_p + pc1];
						const ar_sample_t *s01 = &ps[pr1 * vdim_p + pc0];
						const ar_sample_t *s11 = &ps[pr1 * vdim_p + pc1];
						s->val = (int32_t)(((int64_t)s00->val + (int64_t)s10->val + (int64_t)s01->val + (int64_t)s11->val) >> 2);
						s->flags = s00->flags | s10->flags | s01->flags | s11->flags;
						s->rgba = (((s00->rgba & 0xFCFCFCFCu) >> 2) + ((s10->rgba & 0xFCFCFCFCu) >> 2) +
						           ((s01->rgba & 0xFCFCFCFCu) >> 2) + ((s11->rgba & 0xFCFCFCFCu) >> 2));
						for (uint8_t a = 0; a < AR_SAMPLE_ATTR_CAP; a++) {
							s->attr[a] = (uint8_t)(((uint16_t)s00->attr[a] + (uint16_t)s10->attr[a] +
							                        (uint16_t)s01->attr[a] + (uint16_t)s11->attr[a]) >> 2);
						}
					}
					v_out = s->val;
					rgba_out = s->rgba ? s->rgba : sys->default_rgba;
					node->flags |= (s->flags & AR_NODE_USER_FLAGS_MASK);
				} else {
					if (r_even && c_even) {
						v_out = pv[pr0 * vdim_p + pc0];
						rgba_out = pc_arr[pr0 * vdim_p + pc0];
					} else if (r_even) {
						v_out = (int32_t)(((int64_t)pv[pr0 * vdim_p + pc0] + (int64_t)pv[pr0 * vdim_p + pc1]) >> 1);
						rgba_out = (((pc_arr[pr0 * vdim_p + pc0] & 0xFEFEFEFEu) >> 1) +
						            ((pc_arr[pr0 * vdim_p + pc1] & 0xFEFEFEFEu) >> 1));
					} else if (c_even) {
						v_out = (int32_t)(((int64_t)pv[pr0 * vdim_p + pc0] + (int64_t)pv[pr1 * vdim_p + pc0]) >> 1);
						rgba_out = (((pc_arr[pr0 * vdim_p + pc0] & 0xFEFEFEFEu) >> 1) +
						            ((pc_arr[pr1 * vdim_p + pc0] & 0xFEFEFEFEu) >> 1));
					} else {
						v_out = (int32_t)(((int64_t)pv[pr0 * vdim_p + pc0] + (int64_t)pv[pr0 * vdim_p + pc1] +
						                   (int64_t)pv[pr1 * vdim_p + pc0] + (int64_t)pv[pr1 * vdim_p + pc1]) >> 2);
						rgba_out = (((pc_arr[pr0 * vdim_p + pc0] & 0xFCFCFCFCu) >> 2) +
						            ((pc_arr[pr0 * vdim_p + pc1] & 0xFCFCFCFCu) >> 2) +
						            ((pc_arr[pr1 * vdim_p + pc0] & 0xFCFCFCFCu) >> 2) +
						            ((pc_arr[pr1 * vdim_p + pc1] & 0xFCFCFCFCu) >> 2));
					}
					node->flags |= (parent_node->flags & AR_NODE_USER_FLAGS_MASK);
					if (mut_samples) {
						ar_sample_t *s = &mut_samples[idx];
						memset(s, 0, sizeof(*s));
						s->val = v_out;
						s->rgba = rgba_out;
						s->flags = (uint8_t)(parent_node->flags & AR_NODE_USER_FLAGS_MASK);
					}
				}
				if (v_out < node->min_val) node->min_val = v_out;
				if (v_out > node->max_val) node->max_val = v_out;

				mut_values[idx] = v_out;
				mut_rgba[idx] = rgba_out;
				idx++;
			}
		}
		node->flags |= AR_NODE_FLAG_DEFERRED;
		node->occ_grid_k = 0;
		node->rev++;
#ifdef TTYPT_WQ_H
		if (has_async_wq && !stl_wq_busy(&meta->async_task)) {
			meta->async_sys = sys;
			meta->async_slot = (uint32_t)slot;
			meta->async_key = key;
			meta->async_world = world;
			meta->async_level = level;
			meta->async_ox = ox;
			meta->async_oy = oy;
			meta->async_w = node->w;
			meta->async_h = node->h;
			(void)ar_wq_submit(sys->wq, &meta->async_task, ar_async_populate_worker, meta);
		}
#endif
		if (sys->timing_enabled)
			sys->metrics.us_populate += (uint32_t)(ar_now_us(true) - t_pop);
		return node;
	}

#ifdef TTYPT_WQ_H
	/* Fast 5x5 coarse preview + async full bake when budget exceeded and no cached parent exists */
	if (budget_exceeded && has_async_wq && gdim >= 4) {
		ar_sample_t cgrid[25];
		for (int cr = 0; cr <= 4; cr++) {
			int32_t wy = oy + (int32_t)(((int64_t)cr * node->h) >> 2);
			if (sys->row_sampler) {
				sys->row_sampler(sys->user_ctx, world, level, wy, ox, node->w, 5u, &cgrid[cr * 5]);
			} else {
				for (int cc = 0; cc <= 4; cc++) {
					int32_t wx = ox + (int32_t)(((int64_t)cc * node->w) >> 2);
					sys->sampler(sys->user_ctx, world, level, wx, wy, &cgrid[cr * 5 + cc]);
				}
			}
		}
		for (int ci = 0; ci < 25; ci++) {
			if (!cgrid[ci].rgba)
				cgrid[ci].rgba = sys->default_rgba;
			node->flags |= (cgrid[ci].flags & AR_NODE_USER_FLAGS_MASK);
		}
		size_t idx = 0;
		for (int r = 0; r <= gdim; r++) {
			int32_t ry_q8 = (r * 4 * 256) / gdim;
			int cr0 = ry_q8 >> 8;
			int cr1 = (cr0 < 4) ? (cr0 + 1) : 4;
			int fy = ry_q8 & 255;
			int ify = 256 - fy;
			for (int c = 0; c <= gdim; c++) {
				int32_t cx_q8 = (c * 4 * 256) / gdim;
				int cc0 = cx_q8 >> 8;
				int cc1 = (cc0 < 4) ? (cc0 + 1) : 4;
				int fx = cx_q8 & 255;
				int ifx = 256 - fx;

				const ar_sample_t *s00 = &cgrid[cr0 * 5 + cc0];
				const ar_sample_t *s10 = &cgrid[cr0 * 5 + cc1];
				const ar_sample_t *s01 = &cgrid[cr1 * 5 + cc0];
				const ar_sample_t *s11 = &cgrid[cr1 * 5 + cc1];

				int64_t v0 = ((int64_t)s00->val * ifx + (int64_t)s10->val * fx) >> 8;
				int64_t v1 = ((int64_t)s01->val * ifx + (int64_t)s11->val * fx) >> 8;
				int32_t v_out = (int32_t)((v0 * ify + v1 * fy) >> 8);

				uint32_t c00 = s00->rgba, c10 = s10->rgba, c01 = s01->rgba, c11 = s11->rgba;
				uint32_t r0 = ((c00 & 0xFFu) * (uint32_t)ifx + (c10 & 0xFFu) * (uint32_t)fx) >> 8;
				uint32_t r1 = ((c01 & 0xFFu) * (uint32_t)ifx + (c11 & 0xFFu) * (uint32_t)fx) >> 8;
				uint32_t g0 = (((c00 >> 8) & 0xFFu) * (uint32_t)ifx + ((c10 >> 8) & 0xFFu) * (uint32_t)fx) >> 8;
				uint32_t g1 = (((c01 >> 8) & 0xFFu) * (uint32_t)ifx + ((c11 >> 8) & 0xFFu) * (uint32_t)fx) >> 8;
				uint32_t b0 = (((c00 >> 16) & 0xFFu) * (uint32_t)ifx + ((c10 >> 16) & 0xFFu) * (uint32_t)fx) >> 8;
				uint32_t b1 = (((c01 >> 16) & 0xFFu) * (uint32_t)ifx + ((c11 >> 16) & 0xFFu) * (uint32_t)fx) >> 8;
				const ar_sample_t *nearest = (fy < 128) ? ((fx < 128) ? s00 : s10)
				                                        : ((fx < 128) ? s01 : s11);
				uint32_t rgba_out = (((r0 * (uint32_t)ify + r1 * (uint32_t)fy) >> 8) & 0xFFu) |
				                    ((((g0 * (uint32_t)ify + g1 * (uint32_t)fy) >> 8) & 0xFFu) << 8) |
				                    ((((b0 * (uint32_t)ify + b1 * (uint32_t)fy) >> 8) & 0xFFu) << 16) |
				                    (nearest->rgba & 0xFF000000u);

				if (v_out < node->min_val) node->min_val = v_out;
				if (v_out > node->max_val) node->max_val = v_out;
				mut_values[idx] = v_out;
				mut_rgba[idx] = rgba_out;
				if (mut_samples) {
					ar_sample_t *ds = &mut_samples[idx];
					*ds = *nearest;
					ds->val = v_out;
					ds->rgba = rgba_out;
				}
				idx++;
			}
		}
		node->flags |= AR_NODE_FLAG_DEFERRED;
		node->occ_grid_k = 0;
		node->rev++;
		if (!stl_wq_busy(&meta->async_task)) {
			meta->async_sys = sys;
			meta->async_slot = (uint32_t)slot;
			meta->async_key = key;
			meta->async_world = world;
			meta->async_level = level;
			meta->async_ox = ox;
			meta->async_oy = oy;
			meta->async_w = node->w;
			meta->async_h = node->h;
			(void)ar_wq_submit(sys->wq, &meta->async_task, ar_async_populate_worker, meta);
		}
		if (sys->timing_enabled)
			sys->metrics.us_populate += (uint32_t)(ar_now_us(true) - t_pop);
		return node;
	}
#endif

	sys->frame_populates_count++;
	sys->metrics.frame_populates++;
	if (level < AR_METRICS_MAX_LEVELS)
		sys->metrics.level_populates[level]++;
	node->flags &= ~AR_NODE_FLAG_DEFERRED;

	bool can_reuse_parent = (parent_node != NULL && (!mut_samples || parent_node->samples != NULL));
	bool can_reuse_children = false;
	for (int q = 0; q < 4; q++) {
		if (child_nodes[q] && (!mut_samples || child_nodes[q]->samples != NULL)) {
			can_reuse_children = true;
			break;
		}
	}
	bool has_reusable_nodes = can_reuse_parent || can_reuse_children;

	size_t idx = 0;
	if (sys->row_sampler && !has_reusable_nodes) {
		for (int r = 0; r <= gdim; r++) {
			int32_t wy = oy + (int32_t)(((int64_t)r * node->h) / gdim);
			ar_sample_t *row_dst = mut_samples ? &mut_samples[idx] : sys->row_scratch;
			sys->row_sampler(sys->user_ctx, world, level, wy, ox, node->w, (uint32_t)vdim, row_dst);
			sys->stats.sampler_invocations_total += (uint32_t)vdim;
			sys->metrics.frame_row_sampler_calls++;
			sys->metrics.frame_verts_sampled += (uint32_t)vdim;
			for (int c = 0; c <= gdim; c++) {
				ar_sample_t *s = &row_dst[c];
				if (s->val < node->min_val) node->min_val = s->val;
				if (s->val > node->max_val) node->max_val = s->val;
				node->flags |= (s->flags & AR_NODE_USER_FLAGS_MASK);

				mut_values[idx] = s->val;
				mut_rgba[idx] = s->rgba ? s->rgba : sys->default_rgba;
				idx++;
			}
		}
	} else {
		for (int r = 0; r <= gdim; r++) {
			int32_t wy = oy + (int32_t)(((int64_t)r * node->h) / gdim);

			/* Fast path: if neither child covering row r exists and parent cannot
			 * supply any vertex on row r, batch the entire row via row_sampler! */
			bool row_has_child = false;
			if (can_reuse_children) {
				int qy = (r > half_gdim) ? 1 : 0;
				const ar_node_t *c0 = child_nodes[0 | (qy << 1)];
				const ar_node_t *c1 = child_nodes[1 | (qy << 1)];
				if ((c0 && (!mut_samples || c0->samples != NULL)) ||
				    (c1 && (!mut_samples || c1->samples != NULL)))
					row_has_child = true;
			}
			bool row_has_parent = (can_reuse_parent && ((r & 1) == 0));

			if (sys->row_sampler && !row_has_child && !row_has_parent) {
				ar_sample_t *row_dst = mut_samples ? &mut_samples[idx] : sys->row_scratch;
				sys->row_sampler(sys->user_ctx, world, level, wy, ox, node->w, (uint32_t)vdim, row_dst);
				sys->stats.sampler_invocations_total += (uint32_t)vdim;
				sys->metrics.frame_row_sampler_calls++;
				sys->metrics.frame_verts_sampled += (uint32_t)vdim;
				for (int c = 0; c <= gdim; c++) {
					ar_sample_t *s = &row_dst[c];
					if (s->val < node->min_val) node->min_val = s->val;
					if (s->val > node->max_val) node->max_val = s->val;
					node->flags |= (s->flags & AR_NODE_USER_FLAGS_MASK);

					mut_values[idx] = s->val;
					mut_rgba[idx] = s->rgba ? s->rgba : sys->default_rgba;
					idx++;
				}
				continue;
			}

			int32_t miss_wx[4];
			size_t miss_idx[4];
			int nmiss = 0;

			for (int c = 0; c <= gdim; c++) {
				int32_t wx = ox + (int32_t)(((int64_t)c * node->w) / gdim);
				ar_sample_t tmp_s;
				ar_sample_t *s = mut_samples ? &mut_samples[idx] : &tmp_s;
				bool reused = false;

				if (row_has_child) {
					int qx = (c > half_gdim) ? 1 : 0;
					int qy = (r > half_gdim) ? 1 : 0;
					uint8_t q = (uint8_t)(qx | (qy << 1));
					const ar_node_t *cn = child_nodes[q];
					if (cn && (!mut_samples || cn->samples != NULL)) {
						int cc = (c - (qx ? half_gdim : 0)) * 2;
						int cr = (r - (qy ? half_gdim : 0)) * 2;
						int32_t c_wx = cn->ox + (int32_t)(((int64_t)cc * cn->w) / gdim);
						int32_t c_wy = cn->oy + (int32_t)(((int64_t)cr * cn->h) / gdim);
						if (c_wx == wx && c_wy == wy) {
							if (cn->samples) {
								*s = cn->samples[cr * vdim + cc];
							} else {
								s->val = cn->values[cr * vdim + cc];
								s->rgba = cn->rgba[cr * vdim + cc];
								s->flags = cn->flags & AR_NODE_USER_FLAGS_MASK;
							}
							reused = true;
							sys->metrics.frame_verts_reused_child++;
						}
					}
				}
				if (!reused && row_has_parent && ((c & 1) == 0)) {
					int pr = parent_r_off + (r >> 1);
					int pc = parent_c_off + (c >> 1);
					int32_t p_wx = parent_node->ox + (int32_t)(((int64_t)pc * parent_node->w) / gdim);
					int32_t p_wy = parent_node->oy + (int32_t)(((int64_t)pr * parent_node->h) / gdim);
					if (p_wx == wx && p_wy == wy) {
						if (parent_node->samples) {
							*s = parent_node->samples[pr * vdim + pc];
						} else {
							s->val = parent_node->values[pr * vdim + pc];
							s->rgba = parent_node->rgba[pr * vdim + pc];
							s->flags = parent_node->flags & AR_NODE_USER_FLAGS_MASK;
						}
						reused = true;
						sys->metrics.frame_verts_reused_parent++;
					}
				}

				if (reused) {
					if (s->val < node->min_val) node->min_val = s->val;
					if (s->val > node->max_val) node->max_val = s->val;
					node->flags |= (s->flags & AR_NODE_USER_FLAGS_MASK);

					mut_values[idx] = s->val;
					mut_rgba[idx] = s->rgba ? s->rgba : sys->default_rgba;
					idx++;
				} else if (sys->eval_row4) {
					miss_wx[nmiss] = wx;
					miss_idx[nmiss] = idx;
					nmiss++;
					idx++;
					if (nmiss == 4) {
						ar_sample_t batch[4];
						memset(batch, 0, sizeof(batch));
						sys->eval_row4(sys->user_ctx, world, level, miss_wx, wy, batch);
						sys->stats.sampler_invocations_total += 4u;
						sys->metrics.frame_row_sampler_calls++;
						sys->metrics.frame_verts_sampled += 4u;
						for (int k = 0; k < 4; k++) {
							size_t mi = miss_idx[k];
							if (mut_samples)
								mut_samples[mi] = batch[k];
							if (batch[k].val < node->min_val) node->min_val = batch[k].val;
							if (batch[k].val > node->max_val) node->max_val = batch[k].val;
							node->flags |= (batch[k].flags & AR_NODE_USER_FLAGS_MASK);
							mut_values[mi] = batch[k].val;
							mut_rgba[mi] = batch[k].rgba ? batch[k].rgba : sys->default_rgba;
						}
						nmiss = 0;
					}
				} else {
					sys->sampler(sys->user_ctx, world, level, wx, wy, s);
					sys->stats.sampler_invocations_total++;
					sys->metrics.frame_scalar_sampler_calls++;
					sys->metrics.frame_verts_sampled++;

					if (s->val < node->min_val) node->min_val = s->val;
					if (s->val > node->max_val) node->max_val = s->val;
					node->flags |= (s->flags & AR_NODE_USER_FLAGS_MASK);

					mut_values[idx] = s->val;
					mut_rgba[idx] = s->rgba ? s->rgba : sys->default_rgba;
					idx++;
				}
			}

			for (int k = 0; k < nmiss; k++) {
				size_t mi = miss_idx[k];
				ar_sample_t tmp_s;
				ar_sample_t *s = mut_samples ? &mut_samples[mi] : &tmp_s;
				sys->sampler(sys->user_ctx, world, level, miss_wx[k], wy, s);
				sys->stats.sampler_invocations_total++;
				sys->metrics.frame_scalar_sampler_calls++;
				sys->metrics.frame_verts_sampled++;
				if (s->val < node->min_val) node->min_val = s->val;
				if (s->val > node->max_val) node->max_val = s->val;
				node->flags |= (s->flags & AR_NODE_USER_FLAGS_MASK);
				mut_values[mi] = s->val;
				mut_rgba[mi] = s->rgba ? s->rgba : sys->default_rgba;
			}
		}
	}
	ar_node_compute_occ_cell_min(node);
	node->rev++;

	if (sys->timing_enabled)
		sys->metrics.us_populate += (uint32_t)(ar_now_us(true) - t_pop);
	return node;
}
