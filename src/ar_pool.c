/* ========================================================================= */
/* System Lifecycle, Contiguous Arena & O(1) Intrusive LRU Pool              */
/* ========================================================================= */

static void ar_pool_init_lists(ar_system_t *sys)
{
	sys->lru_head = -1;
	sys->lru_tail = -1;
	sys->pool_in_use = 0;
	if (sys->pool_capacity == 0) {
		sys->free_head = -1;
		return;
	}
	sys->free_head = 0;
	for (uint16_t i = 0; i < sys->pool_capacity; i++) {
		sys->pool_meta[i].lru_prev = -1;
		sys->pool_meta[i].lru_next = (i + 1u < sys->pool_capacity) ? (int32_t)(i + 1u) : -1;
	}
}

static void ar_lru_unlink(ar_system_t *sys, int32_t slot)
{
	ar_node_meta_t *m = &sys->pool_meta[slot];
	int32_t prev = m->lru_prev;
	int32_t next = m->lru_next;
	if (prev >= 0) sys->pool_meta[prev].lru_next = next;
	else           sys->lru_head = next;
	if (next >= 0) sys->pool_meta[next].lru_prev = prev;
	else           sys->lru_tail = prev;
	m->lru_prev = -1;
	m->lru_next = -1;
}

static void ar_lru_insert_head(ar_system_t *sys, int32_t slot)
{
	ar_node_meta_t *m = &sys->pool_meta[slot];
	int32_t old_head = sys->lru_head;
	m->lru_prev = -1;
	m->lru_next = old_head;
	if (old_head >= 0) sys->pool_meta[old_head].lru_prev = slot;
	else               sys->lru_tail = slot;
	sys->lru_head = slot;
}

static void ar_lru_touch(ar_system_t *sys, int32_t slot)
{
	if (slot < 0 || slot == sys->lru_head)
		return;
	ar_lru_unlink(sys, slot);
	ar_lru_insert_head(sys, slot);
}

static int ar_pool_alloc_or_evict(ar_system_t *sys)
{
	if (sys->free_head >= 0) {
		int32_t slot = sys->free_head;
		sys->free_head = sys->pool_meta[slot].lru_next;
		sys->pool_meta[slot].lru_prev = -1;
		sys->pool_meta[slot].lru_next = -1;
		sys->pool_in_use++;
		return (int)slot;
	}

	/* Pool full: walk backwards from lru_tail (oldest first).
	 * Common case terminates on the very first node in O(1). */
	int32_t fallback_slot = -1;
	for (int32_t cur = sys->lru_tail; cur >= 0; cur = sys->pool_meta[cur].lru_prev) {
#ifdef TTYPT_WQ_H
		if (stl_wq_busy(&sys->pool_meta[cur].async_task))
			continue;
#endif
		if (fallback_slot < 0)
			fallback_slot = cur;
		if (sys->pool_meta[cur].active_view_mask == 0u) {
			ar_lru_unlink(sys, cur);
			islet_quad_map_del(&sys->map, sys->pool_meta[cur].morton_key);
			sys->metrics.frame_evictions++;
			return (int)cur;
		}
	}

	if (fallback_slot >= 0) {
		ar_lru_unlink(sys, fallback_slot);
		islet_quad_map_del(&sys->map, sys->pool_meta[fallback_slot].morton_key);
		sys->metrics.frame_evictions++;
		sys->metrics.frame_active_evictions++;
		return (int)fallback_slot;
	}
	return -1;
}

ar_system_t *ar_system_create_ex(const ar_config_t *cfg)
{
	if (!cfg || cfg->base_w <= 0 || cfg->base_h <= 0 ||
	    cfg->grid_dim == 0 || cfg->pool_capacity == 0 ||
	    cfg->max_active == 0 || cfg->max_levels == 0 || !cfg->sampler)
		return NULL;

	uint16_t grid_dim = cfg->grid_dim;
	uint16_t vert_dim = (uint16_t)(grid_dim + 1u);
	uint32_t vert_count = (uint32_t)vert_dim * (uint32_t)vert_dim;
	uint16_t pool_cap = cfg->pool_capacity;
	bool samples_bounded = (cfg->samples_pool_capacity > 0 && cfg->samples_pool_capacity < pool_cap);
	uint16_t samples_cap = samples_bounded ? cfg->samples_pool_capacity : pool_cap;
	uint16_t max_act = cfg->max_active;
	uint8_t max_lvls = cfg->max_levels;
	uint32_t def_rgba = cfg->default_rgba;

	size_t sys_bytes = ar_align8(sizeof(ar_system_t));
	size_t pool_bytes = ar_align8((size_t)pool_cap * sizeof(ar_node_t));
	size_t pool_meta_bytes = ar_align8((size_t)pool_cap * sizeof(ar_node_meta_t));
	size_t active_bytes = ar_align8((size_t)max_act * sizeof(ar_node_t *));
	size_t cands_bytes = ar_align8((size_t)(max_act + 8u) * sizeof(ar_candidate_t));
	size_t active_cands_bytes = ar_align8((size_t)max_act * sizeof(ar_candidate_t));
	size_t heap_bytes = ar_align8((size_t)(max_act + 8u) * sizeof(uint32_t));
	size_t free_bytes = ar_align8((size_t)(max_act + 8u) * sizeof(uint32_t));
	size_t sample_owner_bytes = samples_bounded ? ar_align8((size_t)samples_cap * sizeof(int32_t)) : 0u;
	size_t row_scratch_bytes = samples_bounded ? ar_align8((size_t)vert_dim * sizeof(ar_sample_t)) : 0u;

	size_t per_node_vals = ar_align8((size_t)vert_count * sizeof(int16_t));
	size_t per_node_rgba = ar_align8((size_t)vert_count * sizeof(uint32_t));
	size_t per_node_samples = ar_align8((size_t)vert_count * sizeof(ar_sample_t));
	size_t per_node_gpu_arena = per_node_vals + per_node_rgba;

	size_t total_bytes = sys_bytes + pool_bytes + pool_meta_bytes + active_bytes +
	                     cands_bytes + active_cands_bytes + heap_bytes + free_bytes +
	                     sample_owner_bytes + row_scratch_bytes +
	                     (size_t)pool_cap * per_node_gpu_arena +
	                     (size_t)samples_cap * per_node_samples;

	uint8_t *arena = (uint8_t *)calloc(1, total_bytes);
	if (!arena)
		return NULL;

	ar_system_t *sys = (ar_system_t *)arena;
	uint8_t *cursor = arena + sys_bytes;

	sys->pool = (ar_node_t *)cursor;
	cursor += pool_bytes;

	sys->pool_meta = (ar_node_meta_t *)cursor;
	cursor += pool_meta_bytes;

	sys->active_chunks = (const ar_node_t **)cursor;
	cursor += active_bytes;

	sys->cands_scratch = (ar_candidate_t *)cursor;
	cursor += cands_bytes;

	sys->active_cands_scratch = (ar_candidate_t *)cursor;
	cursor += active_cands_bytes;

	sys->heap_scratch = (uint32_t *)cursor;
	cursor += heap_bytes;

	sys->free_scratch = (uint32_t *)cursor;
	cursor += free_bytes;

	if (samples_bounded) {
		sys->sample_owner = (int32_t *)cursor;
		cursor += sample_owner_bytes;
		for (uint16_t k = 0; k < samples_cap; k++)
			sys->sample_owner[k] = -1;
		sys->row_scratch = (ar_sample_t *)cursor;
		cursor += row_scratch_bytes;
	} else {
		sys->sample_owner = NULL;
		sys->row_scratch = NULL;
	}

	for (uint16_t i = 0; i < pool_cap; i++) {
		sys->pool[i].grid_dim = grid_dim;
		sys->pool[i].vert_dim = vert_dim;
		sys->pool[i].vert_count = vert_count;
		sys->pool[i].values = (const int16_t *)cursor;
		cursor += per_node_vals;
		sys->pool[i].rgba = (const uint32_t *)cursor;
		cursor += per_node_rgba;
		sys->pool[i].samples = NULL;
		sys->pool_meta[i].sample_slot = -1;
	}

	sys->samples_arena = cursor;
	sys->per_node_samples_stride = per_node_samples;
	if (!samples_bounded) {
		for (uint16_t i = 0; i < pool_cap; i++) {
			sys->pool[i].samples = (const ar_sample_t *)(sys->samples_arena + (size_t)i * per_node_samples);
			sys->pool_meta[i].sample_slot = (int32_t)i;
		}
	}

	float val_scale = (cfg->value_scale > 0.0f) ? cfg->value_scale : AR_DEFAULT_VALUE_SCALE;

	sys->base_w = cfg->base_w;
	sys->base_h = cfg->base_h;
	sys->grid_dim = grid_dim;
	sys->vert_dim = vert_dim;
	sys->vert_count = vert_count;
	sys->pool_capacity = pool_cap;
	sys->samples_pool_capacity = samples_cap;
	sys->samples_max_level = cfg->samples_max_level;
	sys->samples_bounded = samples_bounded;
	sys->max_active = max_act;
	sys->max_levels = max_lvls;
	sys->default_rgba = def_rgba;
	sys->value_scale = val_scale;
	sys->max_root_cells = (cfg->max_root_cells > 0) ? cfg->max_root_cells : AR_DEFAULT_MAX_ROOT_CELLS;
	sys->probe_ancestor_cap = (cfg->probe_ancestor_cap > 0) ? cfg->probe_ancestor_cap : AR_DEFAULT_PROBE_ANCESTOR_CAP;
	sys->period_x = cfg->period_x;
	sys->period_y = cfg->period_y;
	sys->z_min = (cfg->z_min != 0.0f || cfg->z_max != 0.0f) ? cfg->z_min : -((float)cfg->base_h * val_scale);
	sys->z_max = (cfg->z_min != 0.0f || cfg->z_max != 0.0f) ? cfg->z_max : ((float)cfg->base_h * val_scale);
	sys->target_quad_px = (cfg->target_quad_px > 0.0f) ? cfg->target_quad_px : AR_DEFAULT_TARGET_QUAD_PX;
	sys->margin = (cfg->margin >= 0.0f) ? cfg->margin : 0.0f;
	sys->sampler = cfg->sampler;
	sys->row_sampler = cfg->row_sampler;
	sys->eval_row4 = cfg->eval_row4;
	sys->sparsity = cfg->sparsity;
	sys->user_ctx = cfg->user_ctx;
	sys->max_populates_per_frame = cfg->max_populates_per_frame;
	sys->frame_populates_count = 0;
	sys->map_capacity = (cfg->map_capacity > 0) ? cfg->map_capacity : ((uint32_t)pool_cap * 2u);
	sys->interp_mode = cfg->interp_mode;
	sys->interp_fn = cfg->interp_fn;
	sys->default_pitch = (cfg->default_pitch > 0.001f) ? cfg->default_pitch : AR_DEFAULT_PITCH;
	sys->default_cam_dist = (cfg->default_cam_dist > 1.0f) ? cfg->default_cam_dist : (float)(cfg->base_w * AR_DEFAULT_CAM_DIST_MULT);
	sys->default_topdown_dist = (cfg->default_topdown_dist > 1.0f) ? cfg->default_topdown_dist : (float)(cfg->base_h * AR_DEFAULT_TOPDOWN_DIST_MULT);
	sys->default_viewport_w = (cfg->default_viewport_w > 0.0f) ? cfg->default_viewport_w : (float)(cfg->base_w * AR_DEFAULT_VIEWPORT_W_MULT);
	sys->default_viewport_h = (cfg->default_viewport_h > 0.0f) ? cfg->default_viewport_h : (float)(cfg->base_h * AR_DEFAULT_VIEWPORT_H_MULT);
	sys->sag_radius = (cfg->sag_radius > 0.0f) ? cfg->sag_radius : 0.0f;
	sys->sag_tangent_slope = (cfg->sag_tangent_slope > 0.0f) ? cfg->sag_tangent_slope : AR_DEFAULT_SAG_TANGENT_SLOPE;
	sys->split_hysteresis = (cfg->split_hysteresis > 0.0f) ? cfg->split_hysteresis : 0.0f;
	sys->bounds_fn = cfg->bounds_fn;
	sys->bounds_3d_fn = cfg->bounds_3d_fn;
	sys->cull_3d_fn = cfg->cull_3d_fn;
	sys->bounds_ctx = cfg->bounds_ctx;
	sys->attr_count = (cfg->attr_count == 0) ? AR_SAMPLE_ATTR_CAP
	                  : ((cfg->attr_count <= AR_SAMPLE_ATTR_CAP) ? cfg->attr_count : 0u);
	sys->attr_interp_mode = cfg->attr_interp_mode;
	sys->interpolate_attrs = cfg->interpolate_attrs;
	sys->current_tick = 1;
	sys->begin_frame_cleared_mask = 0u;
	sys->wq = cfg->wq;
	sys->staging_arena = NULL;
	sys->staging_bytes = 0;

	ar_pool_init_lists(sys);

	if (!islet_quad_map_init(&sys->map, (size_t)sys->map_capacity)) {
		free(arena);
		return NULL;
	}

	if (cfg->wq != NULL) {
		ar_system_set_wq(sys, cfg->wq);
	}

	sys->stats.pool_capacity = pool_cap;
	sys->stats.min_level_observed = 255;
	sys->stats.max_level_observed = 0;
	sys->stats.arena_bytes = total_bytes;

	return sys;
}

ar_system_t *ar_system_create(int32_t base_w, int32_t base_h,
                              uint16_t grid_dim,
                              uint16_t pool_capacity,
                              uint16_t max_active,
                              uint8_t max_levels,
                              ar_sampler_fn sampler,
                              ar_sparsity_fn sparsity,
                              void *user_ctx)
{
	ar_config_t cfg;
	memset(&cfg, 0, sizeof(cfg));
	cfg.base_w = base_w;
	cfg.base_h = base_h;
	cfg.grid_dim = grid_dim;
	cfg.pool_capacity = pool_capacity;
	cfg.max_active = max_active;
	cfg.max_levels = max_levels;
	cfg.sampler = sampler;
	cfg.sparsity = sparsity;
	cfg.user_ctx = user_ctx;
	return ar_system_create_ex(&cfg);
}

#ifdef TTYPT_WQ_H
#if defined(__GNUC__) && !defined(__wasm__)
extern int          stl_wq_submit(stl_wq_t *wq, stl_wq_task_t *task, stl_wq_fn fn, void *arg) __attribute__((weak));
extern unsigned int stl_wq_thread_count(const stl_wq_t *wq) __attribute__((weak));
extern unsigned int stl_wq_step(stl_wq_t *wq, unsigned int max_tasks) __attribute__((weak));
extern void         stl_wq_wait(stl_wq_t *wq, stl_wq_task_t *task) __attribute__((weak));
extern void         stl_wq_drain(stl_wq_t *wq) __attribute__((weak));
#endif
static inline int ar_wq_submit(stl_wq_t *wq, stl_wq_task_t *task, stl_wq_fn fn, void *arg)
{
#if defined(__GNUC__) && !defined(__wasm__)
	if (stl_wq_submit)
		return stl_wq_submit(wq, task, fn, arg);
	if (!task || !fn) return -1;
	task->fn = fn;
	task->arg = arg;
	atomic_store_explicit(&task->state, (uint8_t)STL_WQ_RUNNING, memory_order_relaxed);
	fn(arg);
	atomic_store_explicit(&task->state, (uint8_t)STL_WQ_DONE, memory_order_release);
	return 0;
#else
	return stl_wq_submit(wq, task, fn, arg);
#endif
}
static inline unsigned int ar_wq_thread_count(const stl_wq_t *wq)
{
#if defined(__GNUC__) && !defined(__wasm__)
	return stl_wq_thread_count ? stl_wq_thread_count(wq) : 0u;
#else
	return stl_wq_thread_count(wq);
#endif
}
static inline unsigned int ar_wq_step(stl_wq_t *wq, unsigned int max_tasks)
{
#if defined(__GNUC__) && !defined(__wasm__)
	return stl_wq_step ? stl_wq_step(wq, max_tasks) : 0u;
#else
	return stl_wq_step(wq, max_tasks);
#endif
}
static inline void ar_wq_wait(stl_wq_t *wq, stl_wq_task_t *task)
{
#if defined(__GNUC__) && !defined(__wasm__)
	if (stl_wq_wait)
		stl_wq_wait(wq, task);
	else if (task)
		atomic_store_explicit(&task->state, (uint8_t)STL_WQ_DONE, memory_order_release);
#else
	stl_wq_wait(wq, task);
#endif
}
static inline void ar_wq_drain(stl_wq_t *wq)
{
#if defined(__GNUC__) && !defined(__wasm__)
	if (stl_wq_drain)
		stl_wq_drain(wq);
#else
	stl_wq_drain(wq);
#endif
}
#endif

void ar_system_set_wq(ar_system_t *sys, struct stl_wq *wq)
{
	if (!sys)
		return;
#ifdef TTYPT_WQ_H
	if (sys->wq && sys->wq != wq) {
		ar_wq_drain(sys->wq);
	}
	sys->wq = wq;
	if (wq && !sys->staging_arena) {
		size_t per_node_vals = ar_align8((size_t)sys->vert_count * sizeof(int16_t));
		size_t per_node_rgba = ar_align8((size_t)sys->vert_count * sizeof(uint32_t));
		size_t per_node_samples = ar_align8((size_t)sys->vert_count * sizeof(ar_sample_t));
		size_t per_node_arena = per_node_vals + per_node_rgba + per_node_samples;
		uint8_t *s_arena = (uint8_t *)calloc((size_t)sys->pool_capacity, per_node_arena);
		if (s_arena) {
			sys->staging_arena = s_arena;
			sys->staging_bytes = (size_t)sys->pool_capacity * per_node_arena;
			uint8_t *cur = s_arena;
			for (uint16_t i = 0; i < sys->pool_capacity; i++) {
				sys->pool_meta[i].staging_values = (int16_t *)cur;
				cur += per_node_vals;
				sys->pool_meta[i].staging_rgba = (uint32_t *)cur;
				cur += per_node_rgba;
				sys->pool_meta[i].staging_samples = (ar_sample_t *)cur;
				cur += per_node_samples;
			}
		}
	}
#else
	(void)wq;
#endif
}

void ar_system_destroy(ar_system_t *sys)
{
	if (!sys)
		return;
#ifdef TTYPT_WQ_H
	if (sys->wq) {
		for (size_t i = 0; i < sys->pool_capacity; i++) {
			if (stl_wq_busy(&sys->pool_meta[i].async_task))
				ar_wq_wait(sys->wq, &sys->pool_meta[i].async_task);
		}
	}
#endif
	free(sys->staging_arena);
	islet_quad_map_destroy(&sys->map);
	free(sys);
}

static void ar_node_release_samples(ar_system_t *sys, int slot)
{
	if (!sys || !sys->samples_bounded || slot < 0 || slot >= (int)sys->pool_capacity)
		return;
	int32_t sslot = sys->pool_meta[slot].sample_slot;
	if (sslot >= 0 && sslot < (int32_t)sys->samples_pool_capacity) {
		if (sys->sample_owner[sslot] == slot)
			sys->sample_owner[sslot] = -1;
	}
	sys->pool_meta[slot].sample_slot = -1;
	sys->pool[slot].samples = NULL;
}

static void ar_node_bind_samples(ar_system_t *sys, int slot, uint8_t level)
{
	if (!sys || !sys->samples_bounded || slot < 0 || slot >= (int)sys->pool_capacity)
		return;
	if (level > sys->samples_max_level) {
		ar_node_release_samples(sys, slot);
		return;
	}
	if (sys->pool_meta[slot].sample_slot >= 0) {
		sys->pool[slot].samples = (const ar_sample_t *)(sys->samples_arena +
			(size_t)sys->pool_meta[slot].sample_slot * sys->per_node_samples_stride);
		return;
	}

	int32_t chosen = -1;
	for (uint16_t k = 0; k < sys->samples_pool_capacity; k++) {
		if (sys->sample_owner[k] < 0) {
			chosen = (int32_t)k;
			break;
		}
	}
	if (chosen < 0) {
		uint32_t oldest_inactive = 0xFFFFFFFFu;
		int32_t best_inactive = -1;
		uint32_t oldest_any = 0xFFFFFFFFu;
		int32_t best_any = -1;
		for (uint16_t k = 0; k < sys->samples_pool_capacity; k++) {
			int32_t owner = sys->sample_owner[k];
			if (owner < 0 || owner >= (int32_t)sys->pool_capacity || owner == slot)
				continue;
#ifdef TTYPT_WQ_H
			if (stl_wq_busy(&sys->pool_meta[owner].async_task))
				continue;
#endif
			uint32_t tick = sys->pool_meta[owner].last_used_tick;
			if (sys->pool_meta[owner].active_view_mask == 0 && tick < oldest_inactive) {
				oldest_inactive = tick;
				best_inactive = (int32_t)k;
			}
			if (tick < oldest_any) {
				oldest_any = tick;
				best_any = (int32_t)k;
			}
		}
		chosen = (best_inactive >= 0) ? best_inactive : best_any;
		if (chosen >= 0) {
			int32_t prev_owner = sys->sample_owner[chosen];
			if (prev_owner >= 0 && prev_owner < (int32_t)sys->pool_capacity) {
				sys->pool_meta[prev_owner].sample_slot = -1;
				sys->pool[prev_owner].samples = NULL;
				sys->metrics.frame_samples_evictions++;
			}
		}
	}
	if (chosen >= 0) {
		sys->sample_owner[chosen] = slot;
		sys->pool_meta[slot].sample_slot = chosen;
		sys->pool[slot].samples = (const ar_sample_t *)(sys->samples_arena +
			(size_t)chosen * sys->per_node_samples_stride);
	} else {
		sys->pool_meta[slot].sample_slot = -1;
		sys->pool[slot].samples = NULL;
	}
}

void ar_system_clear(ar_system_t *sys)
{
	if (!sys)
		return;
#ifdef TTYPT_WQ_H
	if (sys->wq) {
		for (size_t i = 0; i < sys->pool_capacity; i++) {
			if (stl_wq_busy(&sys->pool_meta[i].async_task))
				ar_wq_wait(sys->wq, &sys->pool_meta[i].async_task);
			sys->pool_meta[i].async_task.state = STL_WQ_IDLE;
		}
	}
#endif
	for (size_t i = 0; i < sys->pool_capacity; i++) {
		sys->pool_meta[i].state = 0;
		sys->pool_meta[i].last_used_tick = 0;
		sys->pool_meta[i].last_active_tick = 0;
		sys->pool_meta[i].active_view_mask = 0;
		sys->pool_meta[i].prev_active_view_mask = 0;
		if (sys->samples_bounded) {
			sys->pool_meta[i].sample_slot = -1;
			sys->pool[i].samples = NULL;
		}
	}
	if (sys->samples_bounded && sys->sample_owner) {
		for (uint16_t k = 0; k < sys->samples_pool_capacity; k++)
			sys->sample_owner[k] = -1;
	}
	ar_pool_init_lists(sys);
	islet_quad_map_clear(&sys->map);
	sys->active_count = 0;
	sys->begin_frame_cleared_mask = 0u;
	sys->stats.pool_in_use = 0;
	sys->stats.active_chunks_last_frame = 0;
	ar_system_reset_frame_metrics(sys);
}

void ar_system_begin_frame(ar_system_t *sys)
{
	if (!sys)
		return;
	sys->current_tick++;
	sys->frame_populates_count = 0;
	sys->begin_frame_cleared_mask = ~0u;
	ar_system_reset_frame_metrics(sys);

#ifdef TTYPT_WQ_H
	if (sys->wq && sys->staging_arena) {
		if (ar_wq_thread_count(sys->wq) == 0u)
			ar_wq_step(sys->wq, 4u);
		size_t vals_bytes = (size_t)sys->vert_count * sizeof(int16_t);
		size_t rgba_bytes = (size_t)sys->vert_count * sizeof(uint32_t);
		size_t samp_bytes = (size_t)sys->vert_count * sizeof(ar_sample_t);
		for (int32_t cur = sys->lru_head; cur >= 0; cur = sys->pool_meta[cur].lru_next) {
			ar_node_meta_t *meta = &sys->pool_meta[cur];
			meta->prev_active_view_mask = meta->active_view_mask;
			meta->active_view_mask = 0u;
			if (atomic_load_explicit(&meta->async_task.state, memory_order_acquire) == STL_WQ_DONE) {
				ar_node_t *node = &sys->pool[cur];
				if (meta->state == 1 && meta->morton_key == meta->async_key &&
				    (node->flags & AR_NODE_FLAG_DEFERRED)) {
					memcpy((int16_t *)node->values, meta->staging_values, vals_bytes);
					memcpy((uint32_t *)node->rgba, meta->staging_rgba, rgba_bytes);
					if (node->samples && meta->staging_samples)
						memcpy((ar_sample_t *)node->samples, meta->staging_samples, samp_bytes);
					node->min_val = meta->async_min_val;
					node->max_val = meta->async_max_val;
					node->flags = meta->async_flags;
					ar_node_compute_occ_cell_min(node);
					node->rev++;
					sys->stats.sampler_invocations_total += sys->vert_count;
					sys->metrics.frame_deferred_upgrades++;
				}
				atomic_store_explicit(&meta->async_task.state, (uint8_t)STL_WQ_IDLE, memory_order_relaxed);
			}
		}
		return;
	}
#endif
	for (int32_t cur = sys->lru_head; cur >= 0; cur = sys->pool_meta[cur].lru_next) {
		ar_node_meta_t *meta = &sys->pool_meta[cur];
		meta->prev_active_view_mask = meta->active_view_mask;
		meta->active_view_mask = 0u;
	}
}
