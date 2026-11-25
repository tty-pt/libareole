/* ========================================================================= */
/* Bounds Envelopes, Binary Max-Heap LOD Solver & Multi-View Queries         */
/* ========================================================================= */

static inline void ar_compute_bounds_envelope(const ar_system_t *sys,
                                              const ar_cull_opts_t *opts,
                                              float ox, float oy, float w, float h,
                                              float *out_z_min, float *out_z_max)
{
	float z_min = (opts && (opts->z_min != 0.0f || opts->z_max != 0.0f)) ? opts->z_min : sys->z_min;
	float z_max = (opts && (opts->z_min != 0.0f || opts->z_max != 0.0f)) ? opts->z_max : sys->z_max;
	float margin = (opts && opts->margin > 0.0f) ? opts->margin : sys->margin;
	ar_bounds_fn b_fn = (opts && opts->bounds_fn) ? opts->bounds_fn : sys->bounds_fn;
	void *b_ctx = (opts && opts->bounds_fn) ? opts->bounds_ctx : sys->bounds_ctx;
	float sag_r = (opts && opts->sag_radius > 0.0f) ? opts->sag_radius : sys->sag_radius;
	float sag_slope = (opts && opts->sag_tangent_slope > 0.0f) ? opts->sag_tangent_slope : sys->sag_tangent_slope;
	if (sag_slope <= 0.0f) sag_slope = AR_DEFAULT_SAG_TANGENT_SLOPE;

	if (b_fn) {
		b_fn(b_ctx, ox, oy, w, h, &z_min, &z_max);
	} else if (sag_r > 0.0f) {
		float fx = opts ? opts->focal_x : 0.0f;
		float fy = opts ? opts->focal_y : 0.0f;

		float corners[4][2] = {
			{ ox, oy },
			{ ox + w, oy },
			{ ox, oy + h },
			{ ox + w, oy + h }
		};
		float r2_max = 0.0f;
		for (int k = 0; k < 4; k++) {
			float dx = corners[k][0] - fx;
			float dy = corners[k][1] - fy;
			float r2 = dx * dx + dy * dy;
			if (r2 > r2_max) r2_max = r2;
		}
		float r_max = sqrtf(r2_max);
		float delta_sag;
		if (r_max < sag_r) {
			float diff = sag_r * sag_r - r2_max;
			if (diff < 1.0f) diff = 1.0f;
			delta_sag = sag_r - sqrtf(diff);
		} else {
			delta_sag = sag_r + (r_max - sag_r) * sag_slope;
		}
		z_min -= delta_sag;
	}
	*out_z_min = z_min - margin;
	*out_z_max = z_max + margin;
}

static inline void ar_apply_cached_z_bounds(ar_system_t *sys,
                                            const ar_cull_opts_t *opts,
                                            uint16_t world, uint8_t level,
                                            int32_t ox, int32_t oy,
                                            float *inout_z_min, float *inout_z_max)
{
	if (!opts || !opts->use_cached_bounds || level > opts->cached_bounds_max_level)
		return;
	uint64_t key = morton_quad_key(world, level, ox, oy, sys->base_w, sys->base_h);
	uint32_t slot = 0;
	if (!ar_map_get(sys, key, &slot) || slot >= sys->pool_capacity)
		return;
	if (sys->pool_meta[slot].state != 1)
		return;
	const ar_node_t *node = &sys->pool[slot];
	if ((node->flags & AR_NODE_FLAG_DEFERRED) != 0u)
		return;
	float inv_s = 1.0f / sys->value_scale;
	float cz_min = (float)node->min_val * inv_s - opts->cached_bounds_pad;
	float cz_max = (float)node->max_val * inv_s + opts->cached_bounds_pad;
	if (cz_min > *inout_z_min) *inout_z_min = cz_min;
	if (cz_max < *inout_z_max) *inout_z_max = cz_max;
	if (*inout_z_min > *inout_z_max) {
		float mid = 0.5f * (*inout_z_min + *inout_z_max);
		*inout_z_min = mid;
		*inout_z_max = mid;
	}
}

static inline void ar_compute_node_3d_box(ar_system_t *sys,
                                          const ar_cull_opts_t *opts,
                                          uint16_t world, uint8_t level,
                                          int32_t iox, int32_t ioy,
                                          int32_t iw, int32_t ih,
                                          const float eye_xyz[3], bool has_eye,
                                          float out_min[3], float out_max[3],
                                          float out_center[3],
                                          float *out_span_w, float *out_span_h,
                                          float *out_dist_sq,
                                          float *out_z_min, float *out_z_max)
{
	float ox = (float)iox, oy = (float)ioy;
	float w  = (float)iw,  h  = (float)ih;
	ar_bounds_3d_fn b3d = (opts && opts->bounds_3d_fn) ? opts->bounds_3d_fn : sys->bounds_3d_fn;
	void *b_ctx = (opts && (opts->bounds_3d_fn || opts->bounds_fn)) ? opts->bounds_ctx : sys->bounds_ctx;

	if (b3d) {
		float z_min = (opts && (opts->z_min != 0.0f || opts->z_max != 0.0f)) ? opts->z_min : sys->z_min;
		float z_max = (opts && (opts->z_min != 0.0f || opts->z_max != 0.0f)) ? opts->z_max : sys->z_max;
		float margin = (opts && opts->margin > 0.0f) ? opts->margin : sys->margin;
		ar_apply_cached_z_bounds(sys, opts, world, level, iox, ioy, &z_min, &z_max);
		float eff_z_min = z_min - margin;
		float eff_z_max = z_max + margin;
		float sw = w, sh = h;
		b3d(b_ctx, ox, oy, w, h, eff_z_min, eff_z_max,
		    out_min, out_max, out_center, &sw, &sh);
		*out_span_w = (sw > 0.0f) ? sw : w;
		*out_span_h = (sh > 0.0f) ? sh : h;
		if (out_z_min) *out_z_min = eff_z_min;
		if (out_z_max) *out_z_max = eff_z_max;

		float ex = has_eye ? eye_xyz[0] : 0.0f;
		float ey = has_eye ? eye_xyz[1] : 0.0f;
		float ez = has_eye ? eye_xyz[2] : 0.0f;
		float dx = 0.0f, dy = 0.0f, dz = 0.0f;
		if (ex < out_min[0]) dx = out_min[0] - ex; else if (ex > out_max[0]) dx = ex - out_max[0];
		if (ey < out_min[1]) dy = out_min[1] - ey; else if (ey > out_max[1]) dy = ey - out_max[1];
		if (ez < out_min[2]) dz = out_min[2] - ez; else if (ez > out_max[2]) dz = ez - out_max[2];
		*out_dist_sq = dx * dx + dy * dy + dz * dz;
		return;
	}

	float node_z_min, node_z_max;
	ar_compute_bounds_envelope(sys, opts, ox, oy, w, h, &node_z_min, &node_z_max);
	if (opts && opts->use_cached_bounds && level <= opts->cached_bounds_max_level) {
		float margin = (opts->margin > 0.0f) ? opts->margin : sys->margin;
		float z0 = node_z_min + margin;
		float z1 = node_z_max - margin;
		ar_apply_cached_z_bounds(sys, opts, world, level, iox, ioy, &z0, &z1);
		node_z_min = z0 - margin;
		node_z_max = z1 + margin;
	}
	out_min[0] = ox;
	out_min[1] = oy;
	out_min[2] = node_z_min;
	out_max[0] = ox + w;
	out_max[1] = oy + h;
	out_max[2] = node_z_max;
	float cz = opts ? opts->focal_z : (node_z_min + node_z_max) * 0.5f;
	if (cz < node_z_min) cz = node_z_min;
	if (cz > node_z_max) cz = node_z_max;
	out_center[0] = ox + w * 0.5f;
	out_center[1] = oy + h * 0.5f;
	out_center[2] = cz;
	*out_span_w = w;
	*out_span_h = h;
	if (out_z_min) *out_z_min = node_z_min;
	if (out_z_max) *out_z_max = node_z_max;

	if (has_eye) {
		float dx = 0.0f, dy = 0.0f;
		if (eye_xyz[0] < out_min[0]) dx = out_min[0] - eye_xyz[0];
		else if (eye_xyz[0] > out_max[0]) dx = eye_xyz[0] - out_max[0];
		if (eye_xyz[1] < out_min[1]) dy = out_min[1] - eye_xyz[1];
		else if (eye_xyz[1] > out_max[1]) dy = eye_xyz[1] - out_max[1];
		float dz = eye_xyz[2] - cz;
		*out_dist_sq = dx * dx + dy * dy + dz * dz;
	} else {
		float fx = opts ? opts->focal_x : 0.0f;
		float fy = opts ? opts->focal_y : 0.0f;
		float dx = 0.0f, dy = 0.0f;
		if (fx < ox) dx = ox - fx;
		else if (fx > ox + w) dx = fx - (ox + w);
		if (fy < oy) dy = oy - fy;
		else if (fy > oy + h) dy = fy - (oy + h);
		*out_dist_sq = dx * dx + dy * dy;
	}
}

static inline float ar_cand_split_ratio(ar_system_t *sys, uint16_t world,
                                        const ar_candidate_t *cand,
                                        float split_threshold, float hyst,
                                        uint32_t view_mask)
{
	float node_thresh = split_threshold;
	if (hyst > 0.0f) {
		/* Avoid Morton lookups when candidate is clearly outside the hysteresis band */
		float band_lo = split_threshold * (1.0f - hyst);
		float band_hi = split_threshold * (1.0f + hyst);
		if (cand->proj_w >= band_lo && cand->proj_w <= band_hi) {
			uint64_t pkey = morton_quad_key(world, cand->level, cand->ox, cand->oy,
			                                sys->base_w, sys->base_h);
			bool use_prev = (sys->begin_frame_cleared_mask & view_mask) != 0u;
			bool kids_were_active = false;
			for (uint8_t q = 0; q < 4; q++) {
				uint64_t ckey = morton_quad_child(pkey, q);
				uint32_t cslot = 0;
				if (ar_map_get(sys, ckey, &cslot) &&
				    cslot < sys->pool_capacity && sys->pool_meta[cslot].state == 1) {
					uint32_t mask = use_prev ? sys->pool_meta[cslot].prev_active_view_mask
					                         : sys->pool_meta[cslot].active_view_mask;
					if ((mask & view_mask) != 0u) {
						kids_were_active = true;
						break;
					}
				}
			}
			node_thresh *= kids_were_active ? (1.0f - hyst) : (1.0f + hyst);
		} else {
			/* Outside hysteresis deadband: apply same directional scaling so heap priority matches */
			node_thresh *= (cand->proj_w > band_hi) ? (1.0f + hyst) : (1.0f - hyst);
		}
	}
	float ratio = (node_thresh > 1e-5f) ? (cand->proj_w / node_thresh) : 0.0f;
	if (hyst > 0.0f && split_threshold > 1e-5f) {
		bool raw_split = (cand->proj_w > split_threshold);
		bool hyst_split = (ratio > 1.0f);
		if (raw_split != hyst_split)
			sys->metrics.frame_hysteresis_holds++;
	}
	return ratio;
}

static inline bool ar_heap_better(const ar_candidate_t *cands, uint32_t a, uint32_t b)
{
	float ra = cands[a].split_ratio;
	float rb = cands[b].split_ratio;
	if (ra != rb)
		return ra > rb;
	float da = cands[a].dist_sq;
	float db = cands[b].dist_sq;
	if (da != db)
		return da < db;
	return a < b;
}

static inline void ar_heap_push(const ar_candidate_t *cands, uint32_t *heap,
                                size_t *heap_len, uint32_t idx)
{
	size_t i = (*heap_len)++;
	heap[i] = idx;
	while (i > 0) {
		size_t p = (i - 1u) >> 1;
		if (!ar_heap_better(cands, heap[i], heap[p]))
			break;
		uint32_t tmp = heap[p];
		heap[p] = heap[i];
		heap[i] = tmp;
		i = p;
	}
}

static inline uint32_t ar_heap_pop(const ar_candidate_t *cands, uint32_t *heap,
                                   size_t *heap_len)
{
	uint32_t top = heap[0];
	size_t n = --(*heap_len);
	if (n == 0)
		return top;
	heap[0] = heap[n];
	size_t i = 0;
	for (;;) {
		size_t l = (i << 1) + 1u;
		if (l >= n)
			break;
		size_t r = l + 1u;
		size_t best = (r < n && ar_heap_better(cands, heap[r], heap[l])) ? r : l;
		if (!ar_heap_better(cands, heap[best], heap[i]))
			break;
		uint32_t tmp = heap[i];
		heap[i] = heap[best];
		heap[best] = tmp;
		i = best;
	}
	return top;
}

size_t ar_system_query_matrix(ar_system_t *sys, uint16_t world,
                              const float mvp[16],
                              const ar_cull_opts_t *opts,
                              const ar_node_t **out_nodes,
                              size_t max_out)
{
	if (!sys || !mvp || !out_nodes || max_out == 0)
		return 0;

	uint64_t t_query = ar_now_us(sys->timing_enabled);
	sys->current_tick++;
	size_t out_count = 0;
	size_t max_candidates = (max_out < sys->max_active) ? max_out : sys->max_active;

	ar_frustum_t frustum;
	ar_frustum_from_matrix(&frustum, mvp);

	float inv_mvp[16];
	bool has_inv = ar_mat4_invert(inv_mvp, mvp);

	float fx = opts ? opts->focal_x : 0.0f;
	float fy = opts ? opts->focal_y : 0.0f;
	float fz = opts ? opts->focal_z : 0.0f;
	float sag_r = (opts && opts->sag_radius > 0.0f) ? opts->sag_radius : sys->sag_radius;
	float sag_slope = (opts && opts->sag_tangent_slope > 0.0f) ? opts->sag_tangent_slope : sys->sag_tangent_slope;
	if (sag_slope <= 0.0f) sag_slope = AR_DEFAULT_SAG_TANGENT_SLOPE;
	float margin = (opts && opts->margin > 0.0f) ? opts->margin : sys->margin;
	float vp_w = (opts && opts->viewport_w > 0.0f) ? opts->viewport_w : sys->default_viewport_w;
	float vp_h = (opts && opts->viewport_h > 0.0f) ? opts->viewport_h : sys->default_viewport_h;
	int64_t max_root_cells = (opts && opts->max_root_cells > 0)
	                         ? (int64_t)opts->max_root_cells
	                         : (int64_t)((sys->max_root_cells > 0) ? sys->max_root_cells : AR_DEFAULT_MAX_ROOT_CELLS);

	float ray_span_mult = (opts && opts->max_ray_span_mult > 0.0f) ? opts->max_ray_span_mult : AR_DEFAULT_MAX_RAY_SPAN_MULT;
	float coord_bound = (opts && opts->max_coord_bound > 0.0f) ? opts->max_coord_bound : AR_DEFAULT_MAX_COORD_BOUND;

	float z_min = (opts && (opts->z_min != 0.0f || opts->z_max != 0.0f)) ? opts->z_min : sys->z_min;
	float z_max = (opts && (opts->z_min != 0.0f || opts->z_max != 0.0f)) ? opts->z_max : sys->z_max;
	bool is_perspective = (fabsf(mvp[3]) > 1e-5f || fabsf(mvp[7]) > 1e-5f || fabsf(mvp[11]) > 1e-5f);

	float eye_xyz[3] = { fx, fy, fz };
	bool has_eye = false;
	if (is_perspective && has_inv && fabsf(inv_mvp[11]) > 1e-12f) {
		float inv_ew = 1.0f / inv_mvp[11];
		float ex = inv_mvp[8]  * inv_ew;
		float ey = inv_mvp[9]  * inv_ew;
		float ez = inv_mvp[10] * inv_ew;
		if (isfinite(ex) && isfinite(ey) && isfinite(ez)) {
			eye_xyz[0] = ex;
			eye_xyz[1] = ey;
			eye_xyz[2] = ez;
			has_eye = true;
		}
	}

	float focal_pt[4] = { fx, fy, fz, 1.0f };
	float focal_clip[4];
	for (int i = 0; i < 4; i++) {
		focal_clip[i] = mvp[0 * 4 + i] * focal_pt[0] +
		               mvp[1 * 4 + i] * focal_pt[1] +
		               mvp[2 * 4 + i] * focal_pt[2] +
		               mvp[3 * 4 + i] * focal_pt[3];
	}
	float dist = (is_perspective && focal_clip[3] > 0.01f) ? fabsf(focal_clip[3]) : 1.0f;
	float z_span = fabsf(z_max - z_min);
	float max_ray_span = is_perspective ? ((dist * ray_span_mult > z_span) ? (dist * ray_span_mult) : z_span) : coord_bound;

	float len_row0 = sqrtf(mvp[0]*mvp[0] + mvp[4]*mvp[4] + mvp[8]*mvp[8]);
	float len_row1 = sqrtf(mvp[1]*mvp[1] + mvp[5]*mvp[5] + mvp[9]*mvp[9]);
	if (len_row0 < 1e-7f) len_row0 = 0.001f;
	if (len_row1 < 1e-7f) len_row1 = 0.001f;

	ar_bounds_3d_fn b3d = (opts && opts->bounds_3d_fn) ? opts->bounds_3d_fn : sys->bounds_3d_fn;
	ar_bounds_fn b_fn = (opts && opts->bounds_fn) ? opts->bounds_fn : sys->bounds_fn;

	float bb_min_x = 1e30f, bb_max_x = -1e30f;
	float bb_min_y = 1e30f, bb_max_y = -1e30f;
	size_t unproj_pts = 0;
	bool frustum_misses_slab = false;

	float fc[8][3];
	float faabb_min[3] = {  1e30f,  1e30f,  1e30f };
	float faabb_max[3] = { -1e30f, -1e30f, -1e30f };
	bool has_faabb = false;
	if (has_inv) {
		int ok_cnt = 0;
		for (int k = 0; k < 8; k++) {
			float nx = (k & 1) ? 1.0f : -1.0f;
			float ny = (k & 2) ? 1.0f : -1.0f;
			float nz = (k & 4) ? 1.0f : -1.0f;
			if (ar_unproject_ndc(inv_mvp, nx, ny, nz, fc[k]) &&
			    isfinite(fc[k][0]) && isfinite(fc[k][1]) && isfinite(fc[k][2])) {
				for (int a = 0; a < 3; a++) {
					if (fc[k][a] < faabb_min[a]) faabb_min[a] = fc[k][a];
					if (fc[k][a] > faabb_max[a]) faabb_max[a] = fc[k][a];
				}
				ok_cnt++;
			}
		}
		if (ok_cnt == 8) {
			for (int a = 0; a < 3; a++) {
				float eps_lo = 1e-5f * fabsf(faabb_min[a]) + 1e-2f;
				float eps_hi = 1e-5f * fabsf(faabb_max[a]) + 1e-2f;
				faabb_min[a] -= eps_lo;
				faabb_max[a] += eps_hi;
			}
			has_faabb = true;
		}
	}

	if (b3d) {
		uint8_t top_lvl = (opts && opts->max_level > 0 && opts->max_level < sys->max_levels)
			? opts->max_level : (uint8_t)(sys->max_levels - 1);
		float half_gw = (sys->period_x > 0)
			? ((float)sys->period_x * 0.5f)
			: ((float)(sys->base_w << top_lvl) * 2.0f);
		float half_gh = (sys->period_y > 0)
			? ((float)sys->period_y * 0.5f)
			: ((float)(sys->base_h << top_lvl) * 2.0f);
		bb_min_x = fx - half_gw;
		bb_max_x = fx + half_gw;
		bb_min_y = fy - half_gh;
		bb_max_y = fy + half_gh;
	} else {
		bool all_8_ok = has_faabb;
		if (has_inv) {
			if (!opts) {
				float c_near[3], c_far[3];
				if (ar_unproject_ndc(inv_mvp, 0.0f, 0.0f, -1.0f, c_near) &&
				    (ar_unproject_ndc(inv_mvp, 0.0f, 0.0f,  1.0f, c_far) ||
				     ar_unproject_ndc(inv_mvp, 0.0f, 0.0f,  0.0f, c_far))) {
					float c_dz = c_far[2] - c_near[2];
					if (fabsf(c_dz) > 1e-6f) {
						float tc = -c_near[2] / c_dz;
						if (tc >= 0.0f) {
							fx = c_near[0] + tc * (c_far[0] - c_near[0]);
							fy = c_near[1] + tc * (c_far[1] - c_near[1]);
						}
					}
				}
			}

			if (all_8_ok) {
				float slab_z_min = z_min - margin;
				float slab_z_max = z_max + margin;

				if (sag_r > 0.0f) {
					float max_r2 = 0.0f;
					for (int k = 0; k < 8; k++) {
						float dx = fc[k][0] - fx;
						float dy = fc[k][1] - fy;
						float r2 = dx * dx + dy * dy;
						if (r2 > max_r2) max_r2 = r2;
					}
					float r_far = sqrtf(max_r2);
					if (r_far > coord_bound) r_far = coord_bound;
					float delta_sag;
					if (r_far < sag_r) {
						float diff = sag_r * sag_r - r_far * r_far;
						if (diff < 1.0f) diff = 1.0f;
						delta_sag = sag_r - sqrtf(diff);
					} else {
						delta_sag = sag_r + (r_far - sag_r) * sag_slope;
					}
					slab_z_min -= delta_sag;
				}

				for (int k = 0; k < 8; k++) {
					if (fc[k][2] >= slab_z_min && fc[k][2] <= slab_z_max) {
						float ix = fc[k][0];
						float iy = fc[k][1];
						if (ix < -coord_bound) ix = -coord_bound;
						if (ix >  coord_bound) ix =  coord_bound;
						if (iy < -coord_bound) iy = -coord_bound;
						if (iy >  coord_bound) iy =  coord_bound;
						if (ix < bb_min_x) bb_min_x = ix;
						if (ix > bb_max_x) bb_max_x = ix;
						if (iy < bb_min_y) bb_min_y = iy;
						if (iy > bb_max_y) bb_max_y = iy;
						unproj_pts++;
					}
				}

				static const uint8_t edges[12][2] = {
					{0, 1}, {1, 3}, {3, 2}, {2, 0},
					{4, 5}, {5, 7}, {7, 6}, {6, 4},
					{0, 4}, {1, 5}, {2, 6}, {3, 7}
				};
				const float slab_planes[2] = { slab_z_min, slab_z_max };

				for (int e = 0; e < 12; e++) {
					const float *pa = fc[edges[e][0]];
					const float *pb = fc[edges[e][1]];
					float za = pa[2], zb = pb[2];
					float dz = zb - za;
					if (fabsf(dz) <= 1e-6f) continue;
					for (int p = 0; p < 2; p++) {
						float pz = slab_planes[p];
						if ((za - pz) * (zb - pz) <= 0.0f) {
							float t = (pz - za) / dz;
							if (t < 0.0f) t = 0.0f;
							if (t > 1.0f) t = 1.0f;
							float ix = pa[0] + t * (pb[0] - pa[0]);
							float iy = pa[1] + t * (pb[1] - pa[1]);
							if (ix < -coord_bound) ix = -coord_bound;
							if (ix >  coord_bound) ix =  coord_bound;
							if (iy < -coord_bound) iy = -coord_bound;
							if (iy >  coord_bound) iy =  coord_bound;
							if (ix < bb_min_x) bb_min_x = ix;
							if (ix > bb_max_x) bb_max_x = ix;
							if (iy < bb_min_y) bb_min_y = iy;
							if (iy > bb_max_y) bb_max_y = iy;
							unproj_pts++;
						}
					}
				}
			}
		}

		if (unproj_pts >= 2 && bb_max_x >= bb_min_x && bb_max_y >= bb_min_y) {
			bb_min_x -= margin;
			bb_max_x += margin;
			bb_min_y -= margin;
			bb_max_y += margin;
			if (!opts) {
				if (fx < bb_min_x || fx > bb_max_x || fy < bb_min_y || fy > bb_max_y) {
					fx = (bb_min_x + bb_max_x) * 0.5f;
					fy = (bb_min_y + bb_max_y) * 0.5f;
				}
			}
		} else if (all_8_ok && unproj_pts == 0 && !b_fn) {
			frustum_misses_slab = true;
		} else {
			float persp_aspect = (opts && opts->fallback_persp_aspect > 0.0f)
			                     ? opts->fallback_persp_aspect : AR_DEFAULT_FALLBACK_PERSP_ASPECT;
			float ground_w = (2.0f * dist) / len_row0;
			float ground_h = (2.0f * dist * (is_perspective ? persp_aspect : 1.0f)) / len_row1;

			if (sag_r > 0.0f && is_perspective) {
				float max_sag_span = sag_r * 2.0f;
				if (max_ray_span * 2.0f < max_sag_span)
					max_sag_span = max_ray_span * 2.0f;
				if (ground_w < max_sag_span) ground_w = max_sag_span;
				if (ground_h < max_sag_span) ground_h = max_sag_span;
			}

			float half_gw = ground_w * 0.5f + margin;
			float half_gh = ground_h * 0.5f + margin;

			bb_min_x = fx - half_gw;
			bb_max_x = fx + half_gw;
			bb_min_y = fy - half_gh;
			bb_max_y = fy + half_gh;
		}
	}

	ar_occluder_t *occ = (opts && opts->occluder && opts->occluder->valid) ? opts->occluder : NULL;
	float occ_box_pad = opts ? opts->occluder_box_z_pad : 0.0f;
	bool occ_pad_all = (b3d != NULL);

	ar_candidate_t *cands = sys->cands_scratch;
	size_t cand_count = 0;

	ar_cull_3d_fn c3d = (opts && opts->cull_3d_fn) ? opts->cull_3d_fn : sys->cull_3d_fn;
	void *c3d_ctx = (opts && (opts->cull_3d_fn || opts->bounds_3d_fn || opts->bounds_fn)) ? opts->bounds_ctx : sys->bounds_ctx;

	uint8_t root_level = 0;
	if (!frustum_misses_slab) {
		while (root_level + 1 < sys->max_levels) {
			int32_t rw = sys->base_w << root_level;
			int32_t rh = sys->base_h << root_level;
			int32_t gi0 = morton_fdiv(ar_floor_i32(bb_min_x), rw);
			int32_t gi1 = morton_fdiv(ar_floor_i32(bb_max_x), rw);
			int32_t gj0 = morton_fdiv(ar_floor_i32(bb_min_y), rh);
			int32_t gj1 = morton_fdiv(ar_floor_i32(bb_max_y), rh);
			int64_t cols = (int64_t)gi1 - gi0 + 1;
			int64_t rows = (int64_t)gj1 - gj0 + 1;
			if (cols > 0 && rows > 0 && cols * rows <= max_root_cells)
				break;
			root_level++;
		}

		if (opts && opts->min_level > 0 && root_level < opts->min_level)
			root_level = opts->min_level;
		if (opts && opts->max_level > 0 && root_level > opts->max_level)
			root_level = opts->max_level;

		int32_t root_w = sys->base_w << root_level;
		int32_t root_h = sys->base_h << root_level;

		int32_t i0 = morton_fdiv(ar_floor_i32(bb_min_x), root_w);
		int32_t i1 = morton_fdiv(ar_floor_i32(bb_max_x), root_w);
		int32_t j0 = morton_fdiv(ar_floor_i32(bb_min_y), root_h);
		int32_t j1 = morton_fdiv(ar_floor_i32(bb_max_y), root_h);

		uint8_t eff_max_lvl = (opts && opts->max_level > 0 && opts->max_level < sys->max_levels)
			? opts->max_level : (uint8_t)(sys->max_levels - 1);
		while (((int64_t)i1 - i0 + 1) * ((int64_t)j1 - j0 + 1) > max_root_cells &&
		       root_level < eff_max_lvl) {
			root_level++;
			root_w = sys->base_w << root_level;
			root_h = sys->base_h << root_level;
			i0 = morton_fdiv(ar_floor_i32(bb_min_x), root_w);
			i1 = morton_fdiv(ar_floor_i32(bb_max_x), root_w);
			j0 = morton_fdiv(ar_floor_i32(bb_min_y), root_h);
			j1 = morton_fdiv(ar_floor_i32(bb_max_y), root_h);
		}

		if (sys->period_x > 0 && root_w > 0) {
			int64_t max_cells_x = (sys->period_x + root_w - 1) / root_w;
			if (max_cells_x < 1) max_cells_x = 1;
			if ((int64_t)i1 - i0 + 1 > max_cells_x) {
				float mid_x = (bb_min_x + bb_max_x) * 0.5f;
				int32_t ic = morton_fdiv(ar_floor_i32(mid_x), root_w);
				i0 = ic - (int32_t)(max_cells_x / 2);
				i1 = i0 + (int32_t)max_cells_x - 1;
			}
		}

		if (sys->period_y > 0 && root_h > 0) {
			int64_t max_cells_y = (sys->period_y + root_h - 1) / root_h;
			if (max_cells_y < 1) max_cells_y = 1;
			if ((int64_t)j1 - j0 + 1 > max_cells_y) {
				float mid_y = (bb_min_y + bb_max_y) * 0.5f;
				int32_t jc = morton_fdiv(ar_floor_i32(mid_y), root_h);
				j0 = jc - (int32_t)(max_cells_y / 2);
				j1 = j0 + (int32_t)max_cells_y - 1;
			}
		}

		for (int32_t j = j0; j <= j1; j++) {
			for (int32_t i = i0; i <= i1; i++) {
				int32_t rox = i * root_w;
				int32_t roy = j * root_h;
				float b_min[3], b_max[3], center_xyz[3], span_w, span_h, dist_sq;
				float nz_min, nz_max;
				sys->metrics.frame_cands_evaluated++;
				ar_compute_node_3d_box(sys, opts, world, root_level,
				                       rox, roy, root_w, root_h,
				                       eye_xyz, has_eye,
				                       b_min, b_max, center_xyz, &span_w, &span_h, &dist_sq,
				                       &nz_min, &nz_max);

				if (!ar_frustum_test_box_sat(&frustum, has_faabb, faabb_min, faabb_max,
				                             b_min, b_max)) {
					sys->metrics.frame_frustum_culled++;
					continue;
				}
				if (c3d && !c3d(c3d_ctx, (float)rox, (float)roy, (float)root_w, (float)root_h,
				                nz_min, nz_max, center_xyz, span_w, span_h)) {
					sys->metrics.frame_cull3d_culled++;
					continue;
				}
				if (ar_occ_test_padded(sys, occ, b_min, b_max, occ_box_pad, occ_pad_all))
					continue;

				float pw = ar_projected_screen_width_box_3d(mvp, opts, center_xyz,
				                                            eye_xyz, has_eye,
				                                            span_w, span_h, vp_w, vp_h);
				size_t dst_idx = cand_count;
				if (cand_count < max_candidates) {
					cand_count++;
				} else {
					size_t far_idx = 0;
					float far_d2 = cands[0].dist_sq;
					for (size_t k = 1; k < cand_count; k++) {
						if (cands[k].dist_sq > far_d2) {
							far_d2 = cands[k].dist_sq;
							far_idx = k;
						}
					}
					if (dist_sq >= far_d2)
						continue;
					dst_idx = far_idx;
				}
				cands[dst_idx].ox = rox;
				cands[dst_idx].oy = roy;
				cands[dst_idx].level = root_level;
				cands[dst_idx].proj_w = pw;
				cands[dst_idx].dist_sq = dist_sq;
				memcpy(cands[dst_idx].b_min, b_min, sizeof(b_min));
				memcpy(cands[dst_idx].b_max, b_max, sizeof(b_max));
				cands[dst_idx].active = true;
				cands[dst_idx].split_blocked = false;
			}
		}
	}

	float target_quad_px = (opts && opts->target_quad_px > 0.0f) ? opts->target_quad_px : sys->target_quad_px;
	float split_threshold = target_quad_px * (float)sys->grid_dim;
	float hyst = (opts && opts->split_hysteresis > 0.0f) ? opts->split_hysteresis : sys->split_hysteresis;
	if (hyst > AR_DEFAULT_MAX_HYSTERESIS) hyst = AR_DEFAULT_MAX_HYSTERESIS;
	uint32_t view_mask = (opts && opts->view_id < 32) ? (1u << opts->view_id) : 1u;
	uint8_t min_lvl = opts ? opts->min_level : 0;
	uint32_t *heap = sys->heap_scratch;
	size_t heap_len = 0;
	uint32_t *free_slots = sys->free_scratch;
	size_t free_count = 0;
	size_t live_cands = cand_count;

	for (size_t k = 0; k < cand_count; k++) {
		if (!cands[k].active || cands[k].level <= min_lvl) {
			cands[k].split_ratio = 0.0f;
			continue;
		}
		float ratio = ar_cand_split_ratio(sys, world, &cands[k], split_threshold, hyst, view_mask);
		cands[k].split_ratio = ratio;
		if (ratio > 1.0f) {
			ar_heap_push(cands, heap, &heap_len, (uint32_t)k);
		}
	}

	while (heap_len > 0) {
		uint32_t best_idx = ar_heap_pop(cands, heap, &heap_len);
		if (!cands[best_idx].active || cands[best_idx].split_blocked || cands[best_idx].level <= min_lvl)
			continue;

		uint8_t child_lvl = (uint8_t)(cands[best_idx].level - 1);
		int32_t cw = sys->base_w << child_lvl;
		int32_t ch = sys->base_h << child_lvl;
		int32_t pox = cands[best_idx].ox;
		int32_t poy = cands[best_idx].oy;

		struct {
			int32_t ox, oy;
			float   proj_w;
			float   dist_sq;
			float   b_min[3];
			float   b_max[3];
			bool    visible;
		} kids[4];
		size_t vis_kids = 0;

		for (int c = 0; c < 4; c++) {
			kids[c].ox = pox + ((c & 1) ? cw : 0);
			kids[c].oy = poy + ((c & 2) ? ch : 0);
			kids[c].visible = false;
			float center_xyz[3], span_w, span_h, dist_sq, nz_min, nz_max;
			sys->metrics.frame_cands_evaluated++;
			ar_compute_node_3d_box(sys, opts, world, child_lvl,
			                       kids[c].ox, kids[c].oy, cw, ch,
			                       eye_xyz, has_eye,
			                       kids[c].b_min, kids[c].b_max, center_xyz,
			                       &span_w, &span_h, &dist_sq,
			                       &nz_min, &nz_max);
			if (!ar_frustum_test_box_sat(&frustum, has_faabb, faabb_min, faabb_max,
			                             kids[c].b_min, kids[c].b_max)) {
				sys->metrics.frame_frustum_culled++;
				continue;
			}
			if (c3d && !c3d(c3d_ctx, (float)kids[c].ox, (float)kids[c].oy,
			                (float)cw, (float)ch,
			                nz_min, nz_max,
			                center_xyz, span_w, span_h)) {
				sys->metrics.frame_cull3d_culled++;
				continue;
			}
			if (ar_occ_test_padded(sys, occ, kids[c].b_min, kids[c].b_max,
			                       occ_box_pad, occ_pad_all))
				continue;

			kids[c].visible = true;
			kids[c].proj_w = ar_projected_screen_width_box_3d(mvp, opts,
			                                                  center_xyz,
			                                                  eye_xyz, has_eye,
			                                                  span_w, span_h, vp_w, vp_h);
			kids[c].dist_sq = dist_sq;
			vis_kids++;
		}

		if (vis_kids == 0) {
			cands[best_idx].active = false;
			free_slots[free_count++] = best_idx;
			live_cands--;
			continue;
		}

		if (live_cands + vis_kids - 1 > max_candidates) {
			cands[best_idx].split_blocked = true;
			sys->metrics.frame_split_blocks_cap++;
			continue;
		}

		sys->metrics.frame_splits++;

		bool replaced = false;
		for (int c = 0; c < 4; c++) {
			if (!kids[c].visible) continue;
			uint32_t dst_idx;
			if (!replaced) {
				dst_idx = best_idx;
				replaced = true;
			} else {
				dst_idx = (free_count > 0) ? free_slots[--free_count] : (uint32_t)(cand_count++);
				live_cands++;
			}
			cands[dst_idx].ox = kids[c].ox;
			cands[dst_idx].oy = kids[c].oy;
			cands[dst_idx].level = child_lvl;
			cands[dst_idx].proj_w = kids[c].proj_w;
			cands[dst_idx].dist_sq = kids[c].dist_sq;
			memcpy(cands[dst_idx].b_min, kids[c].b_min, sizeof(kids[c].b_min));
			memcpy(cands[dst_idx].b_max, kids[c].b_max, sizeof(kids[c].b_max));
			cands[dst_idx].active = true;
			cands[dst_idx].split_blocked = false;
			cands[dst_idx].cached = false;
			if (child_lvl > min_lvl) {
				float ratio = ar_cand_split_ratio(sys, world, &cands[dst_idx],
				                                  split_threshold, hyst, view_mask);
				cands[dst_idx].split_ratio = ratio;
				if (ratio > 1.0f) {
					ar_heap_push(cands, heap, &heap_len, dst_idx);
				}
			} else {
				cands[dst_idx].split_ratio = 0.0f;
			}
		}
	}

	size_t active_count = 0;
	ar_candidate_t *active_cands = sys->active_cands_scratch;
	for (size_t i = 0; i < cand_count; i++) {
		if (cands[i].active && active_count < max_candidates) {
			active_cands[active_count++] = cands[i];
		}
	}

	for (size_t a = 1; a < active_count; a++) {
		ar_candidate_t key = active_cands[a];
		size_t b = a;
		while (b > 0 && (active_cands[b - 1].dist_sq > key.dist_sq ||
		                 (active_cands[b - 1].dist_sq == key.dist_sq && active_cands[b - 1].level < key.level))) {
			active_cands[b] = active_cands[b - 1];
			b--;
		}
		active_cands[b] = key;
	}

	/* Post-split progressive populate budget:
	 * Preserve 100% of already-cached leaves and the closest max_pops uncached leaves.
	 * Only when uncached leaves exceed max_pops, coarsen the farthest uncached groups
	 * bottom-up (preferring already-cached ancestors that maintain forward progress,
	 * else coalescing uncached siblings) while maintaining exact non-overlapping coverage. */
	uint16_t max_pops = opts ? opts->max_populates_per_query : 0;
	if (max_pops > 0 && active_count > 0) {
		uint32_t uncached_total = 0;
		for (size_t i = 0; i < active_count; i++) {
			uint64_t lkey = morton_quad_key(world, active_cands[i].level,
			                                active_cands[i].ox, active_cands[i].oy,
			                                sys->base_w, sys->base_h);
			uint32_t lslot = 0;
			bool is_cached = (ar_map_get(sys, lkey, &lslot) &&
			                  lslot < sys->pool_capacity && sys->pool_meta[lslot].state == 1);
			active_cands[i].cached = is_cached;
			if (!is_cached)
				uncached_total++;
		}

		while (uncached_total > max_pops) {
			bool merged = false;
			/* Coarsen breadth-first by lowest level (farthest first within each level)
			 * so the quadtree degrades uniformly by 1 level instead of collapsing one quadrant repeatedly. */
			uint8_t min_u_lvl = root_level;
			for (size_t i = 0; i < active_count; i++) {
				if (!active_cands[i].cached && active_cands[i].level < min_u_lvl)
					min_u_lvl = active_cands[i].level;
			}
			for (uint8_t try_lvl = min_u_lvl; try_lvl < root_level && !merged; try_lvl++) {
			for (size_t rev = active_count; rev > 0 && !merged; rev--) {
				size_t idx = rev - 1;
				if (active_cands[idx].cached || active_cands[idx].level != try_lvl)
					continue;

				uint8_t best_plvl = 0;
				int32_t best_pox = 0, best_poy = 0;
				bool best_pcached = false;
				uint32_t best_g_uncached = 0;
				int best_tier = 99;

				for (uint8_t p_lvl = (uint8_t)(active_cands[idx].level + 1); p_lvl <= root_level; p_lvl++) {
					int32_t pw = sys->base_w << p_lvl;
					int32_t ph = sys->base_h << p_lvl;
					int32_t pox = morton_fdiv(active_cands[idx].ox, pw) * pw;
					int32_t poy = morton_fdiv(active_cands[idx].oy, ph) * ph;
					int64_t px1 = (int64_t)pox + (int64_t)pw;
					int64_t py1 = (int64_t)poy + (int64_t)ph;

					uint64_t pkey = morton_quad_key(world, p_lvl, pox, poy,
					                                sys->base_w, sys->base_h);
					uint32_t pslot = 0;
					bool p_cached = (ar_map_get(sys, pkey, &pslot) &&
					                 pslot < sys->pool_capacity && sys->pool_meta[pslot].state == 1);

					uint32_t g_uncached = 0;
					for (size_t j = 0; j < active_count; j++) {
						if (active_cands[j].ox >= pox && (int64_t)active_cands[j].ox < px1 &&
						    active_cands[j].oy >= poy && (int64_t)active_cands[j].oy < py1) {
							if (!active_cands[j].cached) g_uncached++;
						}
					}

					/* Never let a distant cached ancestor (e.g. orbital node) swallow dozens of finer leaves */
					bool valid = false;
					if (p_cached && g_uncached >= 1 && p_lvl <= active_cands[idx].level + 2 &&
					    uncached_total > g_uncached) {
						valid = true;
					} else if (!p_cached && g_uncached >= 2) {
						valid = true;
					} else if (p_cached && g_uncached >= 1 && p_lvl == active_cands[idx].level + 1) {
						valid = true;
					}

					if (valid) {
						best_tier = 0;
						best_plvl = p_lvl;
						best_pox = pox;
						best_poy = poy;
						best_pcached = p_cached;
						best_g_uncached = g_uncached;
						break;
					}
				}

				if (best_tier < 99) {
					int32_t pw = sys->base_w << best_plvl;
					int32_t ph = sys->base_h << best_plvl;
					int64_t px1 = (int64_t)best_pox + (int64_t)pw;
					int64_t py1 = (int64_t)best_poy + (int64_t)ph;

					size_t write_i = 0;
					for (size_t j = 0; j < active_count; j++) {
						if (active_cands[j].ox >= best_pox && (int64_t)active_cands[j].ox < px1 &&
						    active_cands[j].oy >= best_poy && (int64_t)active_cands[j].oy < py1) {
							continue;
						}
						active_cands[write_i++] = active_cands[j];
					}

					ar_candidate_t pcand;
					memset(&pcand, 0, sizeof(pcand));
					float center_xyz[3], span_w, span_h, nz_min, nz_max;
					ar_compute_node_3d_box(sys, opts, world, best_plvl,
					                       best_pox, best_poy, pw, ph,
					                       eye_xyz, has_eye,
					                       pcand.b_min, pcand.b_max, center_xyz,
					                       &span_w, &span_h, &pcand.dist_sq,
					                       &nz_min, &nz_max);
					pcand.ox = best_pox;
					pcand.oy = best_poy;
					pcand.level = best_plvl;
					pcand.active = true;
					pcand.cached = best_pcached;

					/* Insert pcand maintaining front-to-back sort by dist_sq */
					size_t ins = write_i;
					while (ins > 0 && (active_cands[ins - 1].dist_sq > pcand.dist_sq ||
					                   (active_cands[ins - 1].dist_sq == pcand.dist_sq &&
					                    active_cands[ins - 1].level < pcand.level))) {
						active_cands[ins] = active_cands[ins - 1];
						ins--;
					}
					active_cands[ins] = pcand;
					active_count = write_i + 1;

					uint32_t reduction = best_pcached ? best_g_uncached : (best_g_uncached - 1u);
					uncached_total = (uncached_total > reduction) ? (uncached_total - reduction) : 0u;
					sys->metrics.frame_split_blocks_budget++;
					merged = true;
				}
			}
			}
			if (!merged)
				break;
		}
	}

	/* Clear previous active bit for this view (skip if begin_frame already cleared it) */
	if ((sys->begin_frame_cleared_mask & view_mask) == 0u) {
		for (int32_t cur = sys->lru_head; cur >= 0; cur = sys->pool_meta[cur].lru_next) {
			if (sys->pool[cur].world == world) {
				sys->pool_meta[cur].active_view_mask &= ~view_mask;
			}
		}
	}
	sys->begin_frame_cleared_mask &= ~view_mask;

	ar_point_3d_fn pt3d = opts ? opts->point_3d_fn : NULL;
	void *pt3d_ctx = (opts && opts->bounds_ctx) ? opts->bounds_ctx : sys->bounds_ctx;

	for (size_t i = 0; i < active_count && out_count < max_candidates; i++) {
		if (occ && ar_occ_test_padded(sys, occ, active_cands[i].b_min, active_cands[i].b_max,
		                              occ_box_pad, occ_pad_all)) {
			continue;
		}

		ar_node_t *node = ar_get_or_populate(sys, world,
		                                     active_cands[i].level,
		                                     active_cands[i].ox,
		                                     active_cands[i].oy);
		if (node) {
			size_t slot_idx = (size_t)(node - sys->pool);
			sys->pool_meta[slot_idx].active_view_mask |= view_mask;
			sys->pool_meta[slot_idx].last_touch_frame = sys->current_tick;
			out_nodes[out_count] = node;
			sys->active_chunks[out_count] = node;
			out_count++;
			if (node->level > sys->stats.max_level_observed)
				sys->stats.max_level_observed = node->level;
			if (node->level < sys->stats.min_level_observed)
				sys->stats.min_level_observed = node->level;
			if (node->level < AR_METRICS_MAX_LEVELS)
				sys->metrics.level_active[node->level]++;
			if ((node->flags & AR_NODE_FLAG_DEFERRED) != 0u)
				sys->metrics.frame_deferred_active++;

			if (occ && !(opts && opts->occluder_no_write)) {
				uint64_t t_occ = ar_now_us(sys->timing_enabled);
				float texel_span = (float)((node->w > node->h) ? node->w : node->h) /
				                   (float)((node->grid_dim > 0) ? node->grid_dim : 1u);
				float coarse_k = (node->grid_dim < AR_OCC_CELL_K) ? (float)node->grid_dim : (float)AR_OCC_CELL_K;
				float cell_span = (float)((node->w > node->h) ? node->w : node->h) / coarse_k;
				float node_z_bias = opts->occluder_z_bias +
				                    texel_span * opts->occluder_texel_bias +
				                    (cell_span * cell_span * 0.125f) * opts->occluder_concave_inv_r;
				ar_occluder_rasterize_node(occ, node, sys->value_scale, node_z_bias,
				                           pt3d, pt3d_ctx);
				if ((node->flags & (AR_NODE_FLAG_SPARSE | AR_NODE_FLAG_DEFERRED)) == 0u) {
					uint32_t k = (node->grid_dim < AR_OCC_CELL_K) ? node->grid_dim : AR_OCC_CELL_K;
					sys->metrics.frame_occ_nodes_rasterized++;
					sys->metrics.frame_occ_tris_rasterized += 2u * k * k + 8u * k;
				}
				if (sys->timing_enabled)
					sys->metrics.us_occ_raster += (uint32_t)(ar_now_us(true) - t_occ);
			}
		}
	}

	sys->active_count = out_count;
	sys->stats.active_chunks_last_frame = out_count;
	sys->stats.pool_in_use = sys->pool_in_use;

	if (sys->timing_enabled)
		sys->metrics.us_query_total += (uint32_t)(ar_now_us(true) - t_query);
	return out_count;
}

size_t ar_system_query_modelview_projection(ar_system_t *sys, uint16_t world,
                                            const float modelview[16],
                                            const float projection[16],
                                            const ar_cull_opts_t *opts,
                                            const ar_node_t **out_nodes,
                                            size_t max_out)
{
	float mvp[16];
	ar_mat4_multiply(mvp, projection, modelview);
	return ar_system_query_matrix(sys, world, mvp, opts, out_nodes, max_out);
}

static size_t ar_system_update_view(ar_system_t *sys, uint8_t view_id, uint16_t world,
                                    const ar_view_t *view,
                                    const ar_node_t **out_nodes,
                                    size_t max_out)
{
	if (!sys || !view || !out_nodes || max_out == 0)
		return 0;

	float zoom = (view->zoom > 0.0f) ? view->zoom : 1.0f;
	bool has_custom_dist = (view->eff_cam_dist > 0.0f);
	float eff_cam_dist = has_custom_dist ? view->eff_cam_dist : sys->default_cam_dist;
	float eye_dist = eff_cam_dist / zoom;

	float sin_y = sinf(view->yaw);
	float cos_y = cosf(view->yaw);

	float eye[3], target[3], up[3];
	target[0] = view->cam_x;
	target[1] = view->cam_y;
	target[2] = view->cam_z;

	bool is_top_down = (view->pitch < 0.0f) || (view->fov <= 0.001f && view->pitch <= 0.001f);
	if (is_top_down) {
		float topdown_dist = has_custom_dist ? eye_dist : (sys->default_topdown_dist / zoom);
		eye[0] = view->cam_x;
		eye[1] = view->cam_y;
		eye[2] = view->cam_z + topdown_dist;
		up[0] = -sin_y;
		up[1] = -cos_y;
		up[2] = 0.0f;
	} else {
		float pitch = (view->pitch > 0.001f) ? view->pitch : sys->default_pitch;
		float cos_p = cosf(pitch);
		float sin_p = sinf(pitch);
		float d_horiz = eye_dist * cos_p;
		float d_vert  = eye_dist * sin_p;
		eye[0] = view->cam_x + d_horiz * sin_y;
		eye[1] = view->cam_y + d_horiz * cos_y;
		eye[2] = view->cam_z + d_vert;
		up[0] = -sin_p * sin_y;
		up[1] = -sin_p * cos_y;
		up[2] = cos_p;
	}

	float modelview[16], proj[16], mvp[16];
	ar_mat4_lookat(modelview, eye, target, up);

	float vp_w = (view->viewport_w > 0.0f) ? view->viewport_w : sys->default_viewport_w;
	float vp_h = (view->viewport_h > 0.0f) ? view->viewport_h : sys->default_viewport_h;

	if (view->fov <= 0.001f) {
		float hw = (vp_w * 0.5f) / zoom;
		float hh = (vp_h * 0.5f) / zoom;
		float z_bound = fabsf(sys->z_max - sys->z_min);
		float near_z = (view->near_z != 0.0f || view->far_z != 0.0f)
		               ? view->near_z
		               : ((eye_dist * 2.0f > z_bound) ? (-eye_dist * 2.0f) : -z_bound);
		float far_z = (view->near_z != 0.0f || view->far_z != 0.0f)
		              ? view->far_z
		              : ((eye_dist * 2.0f > z_bound) ? (eye_dist * 2.0f) : z_bound);
		ar_mat4_ortho(proj, -hw, hw, -hh, hh, near_z, far_z);
	} else {
		float aspect = vp_w / vp_h;
		float near_z = (view->near_z > 0.0f)
		               ? view->near_z
		               : ((eye_dist * AR_DEFAULT_PERSP_NEAR_RATIO > 1.0f) ? (eye_dist * AR_DEFAULT_PERSP_NEAR_RATIO) : 1.0f);
		float far_z = (view->far_z > near_z)
		              ? view->far_z
		              : ((eye_dist * AR_DEFAULT_PERSP_FAR_DIST_MULT > fabsf(sys->z_max - sys->z_min) * AR_DEFAULT_PERSP_FAR_ZSPAN_MULT)
		                 ? (eye_dist * AR_DEFAULT_PERSP_FAR_DIST_MULT)
		                 : (fabsf(sys->z_max - sys->z_min) * AR_DEFAULT_PERSP_FAR_ZSPAN_MULT));
		ar_mat4_perspective(proj, view->fov, aspect, near_z, far_z);
	}
	ar_mat4_multiply(mvp, proj, modelview);

	ar_cull_opts_t opts;
	memset(&opts, 0, sizeof(opts));
	opts.z_min = sys->z_min;
	opts.z_max = sys->z_max;
	opts.sag_radius = (view->sag_radius > 0.0f) ? view->sag_radius : sys->sag_radius;
	opts.sag_tangent_slope = sys->sag_tangent_slope;
	opts.focal_x = view->cam_x;
	opts.focal_y = view->cam_y;
	opts.focal_z = view->cam_z;
	opts.margin = (view->margin > 0.0f) ? view->margin : ((view->fov <= 0.001f) ? 0.0f : sys->margin);
	opts.target_quad_px = (view->target_quad_px > 0.0f) ? view->target_quad_px : sys->target_quad_px;
	opts.viewport_w = vp_w;
	opts.viewport_h = vp_h;
	opts.min_level = view->min_level;
	opts.max_level = view->max_level;
	opts.view_id = (uint8_t)(view_id & 31u);
	opts.split_hysteresis = (view->split_hysteresis > 0.0f) ? view->split_hysteresis : sys->split_hysteresis;
	opts.bounds_fn = sys->bounds_fn;
	opts.bounds_3d_fn = sys->bounds_3d_fn;
	opts.bounds_ctx = sys->bounds_ctx;

	return ar_system_query_matrix(sys, world, mvp, &opts, out_nodes, max_out);
}

size_t ar_system_update(ar_system_t *sys, uint16_t world,
                        const ar_view_t *view,
                        const ar_node_t **out_nodes,
                        size_t max_out)
{
	return ar_system_update_view(sys, 0, world, view, out_nodes, max_out);
}

const ar_node_t *ar_system_find(ar_system_t *sys, uint16_t world,
                                uint8_t level, int32_t ox, int32_t oy)
{
	if (!sys) return NULL;
	uint64_t key = morton_quad_key(world, level, ox, oy, sys->base_w, sys->base_h);
	return ar_system_find_by_key(sys, key);
}

const ar_node_t *ar_system_find_by_key(ar_system_t *sys, uint64_t morton_key)
{
	if (!sys) return NULL;
	uint32_t slot = 0;
	if (ar_map_get(sys, morton_key, &slot)) {
		if (slot < sys->pool_capacity && sys->pool_meta[slot].state == 1)
			return &sys->pool[slot];
	}
	return NULL;
}

const ar_node_t *ar_system_ensure(ar_system_t *sys, uint16_t world,
                                  uint8_t level, int32_t ox, int32_t oy)
{
	if (!sys) return NULL;
	return ar_get_or_populate(sys, world, level, ox, oy);
}

size_t ar_system_get_active_count(const ar_system_t *sys)
{
	return sys ? sys->active_count : 0;
}

const ar_node_t *ar_system_get_active_node(const ar_system_t *sys, size_t index)
{
	if (!sys || index >= sys->active_count)
		return NULL;
	return sys->active_chunks[index];
}

static inline ptrdiff_t ar_node_slot_index(const ar_system_t *sys, const ar_node_t *node)
{
	if (!sys || !node || !sys->pool)
		return -1;
	ptrdiff_t idx = node - sys->pool;
	if (idx < 0 || (size_t)idx >= sys->pool_capacity)
		return -1;
	if (sys->pool_meta[idx].state != 1)
		return -1;
	return idx;
}

uint64_t ar_system_node_morton_key(const ar_system_t *sys, const ar_node_t *node)
{
	ptrdiff_t idx = ar_node_slot_index(sys, node);
	return (idx >= 0) ? sys->pool_meta[idx].morton_key : 0u;
}

uint32_t ar_system_node_active_views(const ar_system_t *sys, const ar_node_t *node)
{
	ptrdiff_t idx = ar_node_slot_index(sys, node);
	return (idx >= 0) ? sys->pool_meta[idx].active_view_mask : 0u;
}

bool ar_system_node_is_active_in_view(const ar_system_t *sys, const ar_node_t *node, uint8_t view_id)
{
	if (view_id >= 32u)
		return false;
	return (ar_system_node_active_views(sys, node) & (1u << view_id)) != 0u;
}

size_t ar_system_eval_view(ar_system_t *sys, uint8_t view_id, uint16_t world,
                           const ar_view_t *view,
                           ar_draw_item_t *out_items,
                           size_t max_items)
{
	if (!sys || !view || !out_items || max_items == 0)
		return 0;

	size_t cap = max_items < sys->max_active ? max_items : sys->max_active;
	size_t count = ar_system_update_view(sys, view_id, world, view, sys->active_chunks, cap);

	for (size_t i = 0; i < count; i++) {
		const ar_node_t *node = sys->active_chunks[i];
		out_items[i].node = node;
		out_items[i].world_bounds[0] = (float)node->ox;
		out_items[i].world_bounds[1] = (float)node->oy;
		out_items[i].world_bounds[2] = (float)(node->ox + node->w);
		out_items[i].world_bounds[3] = (float)(node->oy + node->h);
		out_items[i].uv_offset[0] = 0.0f;
		out_items[i].uv_offset[1] = 0.0f;
		out_items[i].uv_scale = 1.0f;
		out_items[i].view_lod = node->level;
	}
	return count;
}

size_t ar_system_eval_ortho_ex(ar_system_t *sys, uint8_t view_id, uint16_t world,
                               const ar_ortho_opts_t *opts)
{
	if (!sys || !opts || opts->world_w <= 0.0f || opts->world_h <= 0.0f)
		return 0;

	sys->current_tick++;
	uint32_t view_mask = 1u << (view_id & 31u);
	if ((sys->begin_frame_cleared_mask & view_mask) == 0u) {
		for (int32_t cur = sys->lru_head; cur >= 0; cur = sys->pool_meta[cur].lru_next) {
			if (sys->pool[cur].world == world) {
				sys->pool_meta[cur].active_view_mask &= ~view_mask;
			}
		}
	}
	sys->begin_frame_cleared_mask &= ~view_mask;

	float world_w = opts->world_w;
	float world_h = opts->world_h;
	float vp_w = (opts->viewport_w > 1.0f) ? opts->viewport_w : world_w;
	float vp_h = (opts->viewport_h > 1.0f) ? opts->viewport_h : world_h;
	float units_per_px_x = world_w / vp_w;
	float units_per_px_y = world_h / vp_h;
	float target_sample_step = (units_per_px_x > units_per_px_y) ? units_per_px_x : units_per_px_y;
	float quad_thresh = (opts->target_quad_px > 0.0f)
	                    ? opts->target_quad_px
	                    : ((sys->target_quad_px > 0.0f) ? sys->target_quad_px : AR_DEFAULT_ORTHO_QUAD_PX);
	int64_t max_root_cells = (opts->max_root_cells > 0)
	                         ? (int64_t)opts->max_root_cells
	                         : ((sys->max_root_cells > 0) ? (int64_t)sys->max_root_cells : AR_DEFAULT_MAX_ROOT_CELLS);

	float min_x = opts->center_x - world_w * 0.5f;
	float max_x = opts->center_x + world_w * 0.5f;
	float min_y = opts->center_y - world_h * 0.5f;
	float max_y = opts->center_y + world_h * 0.5f;

	uint8_t sys_max_lvl = (uint8_t)(sys->max_levels - 1u);
	uint8_t max_lvl = (opts->max_level > 0 && opts->max_level < sys_max_lvl) ? opts->max_level : sys_max_lvl;
	uint8_t level = (opts->min_level <= max_lvl) ? opts->min_level : max_lvl;
	while (level < max_lvl) {
		int32_t cw = sys->base_w << level;
		int32_t ch = sys->base_h << level;
		int32_t i0 = morton_fdiv((int32_t)min_x, cw);
		int32_t i1 = morton_fdiv((int32_t)max_x, cw);
		int32_t j0 = morton_fdiv((int32_t)min_y, ch);
		int32_t j1 = morton_fdiv((int32_t)max_y, ch);
		int64_t cols = (int64_t)i1 - i0 + 1;
		int64_t rows = (int64_t)j1 - j0 + 1;

		int32_t next_cw = sys->base_w << (level + 1);
		float next_step = (float)next_cw / (float)sys->grid_dim;

		if (cols > 0 && rows > 0 && cols * rows <= max_root_cells && next_step > target_sample_step * quad_thresh)
			break;
		level++;
	}

	int32_t cw = sys->base_w << level;
	int32_t ch = sys->base_h << level;
	int32_t i0 = morton_fdiv((int32_t)min_x, cw);
	int32_t i1 = morton_fdiv((int32_t)max_x, cw);
	int32_t j0 = morton_fdiv((int32_t)min_y, ch);
	int32_t j1 = morton_fdiv((int32_t)max_y, ch);

	size_t count = 0;
	for (int32_t j = j0; j <= j1; j++) {
		for (int32_t i = i0; i <= i1; i++) {
			ar_node_t *node = ar_get_or_populate(sys, world, level, i * cw, j * ch);
			if (node) {
				size_t slot_idx = (size_t)(node - sys->pool);
				sys->pool_meta[slot_idx].active_view_mask |= view_mask;
				sys->pool_meta[slot_idx].last_touch_frame = sys->current_tick;
				if (count < sys->max_active)
					sys->active_chunks[count] = node;
				count++;
			}
		}
	}
	sys->active_count = (count < sys->max_active) ? count : sys->max_active;
	sys->stats.pool_in_use = sys->pool_in_use;
	return count;
}

size_t ar_system_eval_ortho(ar_system_t *sys, uint8_t view_id, uint16_t world,
                            float center_x, float center_y,
                            float world_w, float world_h,
                            float viewport_w, float viewport_h)
{
	ar_ortho_opts_t opts;
	memset(&opts, 0, sizeof(opts));
	opts.center_x = center_x;
	opts.center_y = center_y;
	opts.world_w = world_w;
	opts.world_h = world_h;
	opts.viewport_w = viewport_w;
	opts.viewport_h = viewport_h;
	return ar_system_eval_ortho_ex(sys, view_id, world, &opts);
}
