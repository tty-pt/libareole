/* ========================================================================= */
/* Single-Pass Point Probing, Interpolation & Telemetry Snapshot             */
/* ========================================================================= */

static inline float ar_eval_mode_weight(ar_interp_mode_t mode, float t)
{
	switch (mode) {
	case AR_INTERP_LINEAR:
		return t;
	case AR_INTERP_QUINTIC:
		return t * t * t * (t * (t * 6.0f - 15.0f) + 10.0f);
	case AR_INTERP_NEAREST:
		return (t >= 0.5f) ? 1.0f : 0.0f;
	case AR_INTERP_SMOOTHSTEP:
	default:
		return t * t * (3.0f - 2.0f * t);
	}
}

static inline float ar_eval_interp_weight(const ar_system_t *sys, float t)
{
	if (sys->interp_fn) {
		float w = sys->interp_fn(t);
		if (w < 0.0f) return 0.0f;
		if (w > 1.0f) return 1.0f;
		return w;
	}
	return ar_eval_mode_weight(sys->interp_mode, t);
}

/* Single-pass hierarchical ancestor lookup: checks active view bit and records
 * finest cached ancestor in the same walk (halving worst-case hash lookups). */
static const ar_node_t *ar_probe_find_node(ar_system_t *sys, uint32_t view_mask, uint16_t world,
                                           int32_t wx, int32_t wy, uint8_t max_ancestor_override)
{
	uint64_t key = morton_quad_key(world, 0, wx, wy, sys->base_w, sys->base_h);
	uint8_t max_lvl = (sys->max_levels > 0) ? (uint8_t)(sys->max_levels - 1) : 16;
	uint8_t def_cap = (max_ancestor_override > 0) ? max_ancestor_override : sys->probe_ancestor_cap;
	if (def_cap == 0) def_cap = AR_DEFAULT_PROBE_ANCESTOR_CAP;
	uint8_t cap_lvl = (view_mask != 0) ? max_lvl : (max_lvl < def_cap ? max_lvl : def_cap);
	uint8_t walk_max = (view_mask != 0) ? max_lvl : cap_lvl;

	const ar_node_t *first_cached = NULL;
	for (uint8_t lvl = 0; lvl <= walk_max; lvl++) {
		uint32_t slot = 0;
		sys->metrics.frame_probe_hash_lookups++;
		if (ar_map_get(sys, key, &slot) &&
		    slot < sys->pool_capacity && sys->pool_meta[slot].state == 1) {
			const ar_node_t *cand = &sys->pool[slot];
			if (view_mask == 0u)
				return cand;
			if ((sys->pool_meta[slot].active_view_mask & view_mask) != 0u)
				return cand;
			if (!first_cached && lvl <= cap_lvl)
				first_cached = cand;
		}
		if (lvl >= walk_max || lvl >= 16u)
			break;
		key = morton_quad_parent(key);
	}
	return first_cached;
}

static bool ar_probe_val_internal(ar_system_t *sys, uint32_t view_mask, uint16_t world,
                                  int32_t wx, int32_t wy, uint8_t max_ancestor_override,
                                  int16_t *out_val)
{
	if (!sys || !out_val)
		return false;

	uint64_t t0 = ar_now_us(sys->timing_enabled);
	sys->metrics.frame_probes++;

	const ar_node_t *best_node = ar_probe_find_node(sys, view_mask, world, wx, wy, max_ancestor_override);
	if (best_node) {
		int gdim = (int)sys->grid_dim;
		int vdim = (int)sys->vert_dim;
		float fx = ((float)(wx - best_node->ox) / (float)best_node->w) * (float)gdim;
		float fy = ((float)(wy - best_node->oy) / (float)best_node->h) * (float)gdim;

		int col0 = (int)fx;
		int row0 = (int)fy;
		if (col0 < 0) col0 = 0;
		if (col0 >= gdim) col0 = gdim - 1;
		if (row0 < 0) row0 = 0;
		if (row0 >= gdim) row0 = gdim - 1;

		float u = fx - (float)col0;
		float v = fy - (float)row0;
		if (u < 0.0f) u = 0.0f;
		if (u > 1.0f) u = 1.0f;
		if (v < 0.0f) v = 0.0f;
		if (v > 1.0f) v = 1.0f;
		float wu = ar_eval_interp_weight(sys, u);
		float wv = ar_eval_interp_weight(sys, v);

		const int16_t *r0 = &best_node->values[row0 * vdim + col0];
		const int16_t *r1 = r0 + vdim;
		float h0 = (float)r0[0] * (1.0f - wu) + (float)r0[1] * wu;
		float h1 = (float)r1[0] * (1.0f - wu) + (float)r1[1] * wu;
		float h  = h0 * (1.0f - wv) + h1 * wv;
		*out_val = (int16_t)(h >= 0.0f ? (h + 0.5f) : (h - 0.5f));
		if (sys->timing_enabled)
			sys->metrics.us_probe += (uint32_t)(ar_now_us(true) - t0);
		return true;
	}

	ar_sample_t tmp;
	sys->sampler(sys->user_ctx, world, 0, wx, wy, &tmp);
	sys->metrics.frame_probe_sampler_fallbacks++;
	*out_val = tmp.val;
	if (sys->timing_enabled)
		sys->metrics.us_probe += (uint32_t)(ar_now_us(true) - t0);
	return true;
}

static bool ar_probe_internal(ar_system_t *sys, uint32_t view_mask, uint16_t world,
                              int32_t wx, int32_t wy, uint8_t max_ancestor_override,
                              ar_sample_t *out)
{
	if (!sys || !out)
		return false;

	uint64_t t0 = ar_now_us(sys->timing_enabled);
	sys->metrics.frame_probes++;

	const ar_node_t *best_node = ar_probe_find_node(sys, view_mask, world, wx, wy, max_ancestor_override);

	if (best_node) {
		int gdim = (int)sys->grid_dim;
		int vdim = (int)sys->vert_dim;
		float fx = ((float)(wx - best_node->ox) / (float)best_node->w) * (float)gdim;
		float fy = ((float)(wy - best_node->oy) / (float)best_node->h) * (float)gdim;

		int col0 = (int)fx;
		int row0 = (int)fy;
		if (col0 < 0) col0 = 0;
		if (col0 >= gdim) col0 = gdim - 1;
		if (row0 < 0) row0 = 0;
		if (row0 >= gdim) row0 = gdim - 1;

		int col1 = col0 + 1;
		int row1 = row0 + 1;

		float u = fx - (float)col0;
		float v = fy - (float)row0;
		if (u < 0.0f) u = 0.0f;
		if (u > 1.0f) u = 1.0f;
		if (v < 0.0f) v = 0.0f;
		if (v > 1.0f) v = 1.0f;
		float wu = ar_eval_interp_weight(sys, u);
		float wv = ar_eval_interp_weight(sys, v);

		int idx00 = row0 * vdim + col0;
		int idx10 = row0 * vdim + col1;
		int idx01 = row1 * vdim + col0;
		int idx11 = row1 * vdim + col1;

		float h0 = (float)best_node->values[idx00] * (1.0f - wu) + (float)best_node->values[idx10] * wu;
		float h1 = (float)best_node->values[idx01] * (1.0f - wu) + (float)best_node->values[idx11] * wu;
		float h  = h0 * (1.0f - wv) + h1 * wv;
		out->val = (int16_t)(h >= 0.0f ? (h + 0.5f) : (h - 0.5f));

		uint8_t n_attr = sys->attr_count;
		uint32_t c00, c10, c01, c11;
		if (best_node->samples) {
			const ar_sample_t *s00 = &best_node->samples[idx00];
			const ar_sample_t *s10 = &best_node->samples[idx10];
			const ar_sample_t *s01 = &best_node->samples[idx01];
			const ar_sample_t *s11 = &best_node->samples[idx11];
			const ar_sample_t *s_near = (u < 0.5f)
			                            ? ((v < 0.5f) ? s00 : s01)
			                            : ((v < 0.5f) ? s10 : s11);
			bool interp_attr = sys->interpolate_attrs ||
			                   (sys->attr_interp_mode == AR_INTERP_LINEAR) ||
			                   (sys->attr_interp_mode == AR_INTERP_QUINTIC);
			float awu = wu, awv = wv;
			if (interp_attr && sys->attr_interp_mode != 0) {
				awu = ar_eval_mode_weight(sys->attr_interp_mode, u);
				awv = ar_eval_mode_weight(sys->attr_interp_mode, v);
			}
			for (uint8_t k = 0; k < AR_SAMPLE_ATTR_CAP; k++) {
				if (k >= n_attr) {
					out->attr[k] = 0;
				} else if (interp_attr) {
					float a0 = (float)s00->attr[k] * (1.0f - awu) + (float)s10->attr[k] * awu;
					float a1 = (float)s01->attr[k] * (1.0f - awu) + (float)s11->attr[k] * awu;
					int aval = (int)(a0 * (1.0f - awv) + a1 * awv + 0.5f);
					if (aval < 0) aval = 0;
					if (aval > 255) aval = 255;
					out->attr[k] = (uint8_t)aval;
				} else {
					out->attr[k] = s_near->attr[k];
				}
			}
			c00 = s00->rgba; c10 = s10->rgba; c01 = s01->rgba; c11 = s11->rgba;
			out->flags = s_near->flags;
		} else {
			ar_sample_t tmp_s;
			memset(&tmp_s, 0, sizeof(tmp_s));
			sys->sampler(sys->user_ctx, world, best_node->level, wx, wy, &tmp_s);
			sys->metrics.frame_probe_sampler_fallbacks++;
			for (uint8_t k = 0; k < AR_SAMPLE_ATTR_CAP; k++) {
				out->attr[k] = (k < n_attr) ? tmp_s.attr[k] : 0;
			}
			c00 = best_node->rgba[idx00];
			c10 = best_node->rgba[idx10];
			c01 = best_node->rgba[idx01];
			c11 = best_node->rgba[idx11];
			out->flags = tmp_s.flags;
		}

		uint32_t r_b = 0, g_b = 0, b_b = 0, a_b = 0;
		for (int ch = 0; ch < 4; ch++) {
			int shift = ch * 8;
			float v00 = (float)((c00 >> shift) & 0xFF);
			float v10 = (float)((c10 >> shift) & 0xFF);
			float v01 = (float)((c01 >> shift) & 0xFF);
			float v11 = (float)((c11 >> shift) & 0xFF);
			float top = v00 * (1.0f - wu) + v10 * wu;
			float bot = v01 * (1.0f - wu) + v11 * wu;
			float val = top * (1.0f - wv) + bot * wv;
			int ival = (int)(val + 0.5f);
			if (ival < 0) ival = 0;
			if (ival > 255) ival = 255;
			if (ch == 0) r_b = (uint32_t)ival;
			else if (ch == 1) g_b = (uint32_t)ival << 8;
			else if (ch == 2) b_b = (uint32_t)ival << 16;
			else a_b = (uint32_t)ival << 24;
		}
		out->rgba = r_b | g_b | b_b | a_b;
		if (sys->timing_enabled)
			sys->metrics.us_probe += (uint32_t)(ar_now_us(true) - t0);
		return true;
	}

	sys->sampler(sys->user_ctx, world, 0, wx, wy, out);
	sys->metrics.frame_probe_sampler_fallbacks++;
	for (uint8_t k = sys->attr_count; k < AR_SAMPLE_ATTR_CAP; k++)
		out->attr[k] = 0;
	if (sys->timing_enabled)
		sys->metrics.us_probe += (uint32_t)(ar_now_us(true) - t0);
	return true;
}

bool ar_system_probe_view(ar_system_t *sys, uint8_t view_id, uint16_t world,
                          int32_t wx, int32_t wy,
                          ar_sample_t *out)
{
	uint32_t view_mask = 1u << (view_id & 31u);
	return ar_probe_internal(sys, view_mask, world, wx, wy, 0, out);
}

bool ar_system_probe_ex(ar_system_t *sys, uint16_t world,
                        int32_t wx, int32_t wy,
                        uint8_t max_ancestor_levels,
                        ar_sample_t *out)
{
	return ar_probe_internal(sys, 0, world, wx, wy, max_ancestor_levels, out);
}

bool ar_system_probe(ar_system_t *sys, uint16_t world,
                     int32_t wx, int32_t wy,
                     ar_sample_t *out)
{
	return ar_probe_internal(sys, 0, world, wx, wy, 0, out);
}

bool ar_system_probe_val(ar_system_t *sys, uint16_t world,
                         int32_t wx, int32_t wy,
                         int16_t *out_val)
{
	return ar_probe_val_internal(sys, 0, world, wx, wy, 0, out_val);
}

bool ar_system_probe_view_val(ar_system_t *sys, uint8_t view_id, uint16_t world,
                              int32_t wx, int32_t wy,
                              int16_t *out_val)
{
	uint32_t view_mask = 1u << (view_id & 31u);
	return ar_probe_val_internal(sys, view_mask, world, wx, wy, 0, out_val);
}

float ar_system_get_value_scale(const ar_system_t *sys)
{
	return (sys && sys->value_scale > 0.0f) ? sys->value_scale : AR_DEFAULT_VALUE_SCALE;
}

float ar_system_value_to_float(const ar_system_t *sys, int16_t val)
{
	return ar_value_to_float_ex(val, ar_system_get_value_scale(sys));
}

int16_t ar_system_value_from_float(const ar_system_t *sys, float val)
{
	return ar_value_from_float_ex(val, ar_system_get_value_scale(sys));
}

float ar_system_level_texel_span(const ar_system_t *sys, uint8_t level)
{
	if (!sys)
		return 0.0f;
	return ar_level_texel_span(sys->base_w, sys->grid_dim, level);
}

void ar_system_get_stats(const ar_system_t *sys, ar_stats_t *out_stats)
{
	if (!sys || !out_stats)
		return;
	*out_stats = sys->stats;
}

void ar_system_reset_frame_metrics(ar_system_t *sys)
{
	if (!sys)
		return;
	memset(&sys->metrics, 0, sizeof(sys->metrics));
}

void ar_system_set_timing_enabled(ar_system_t *sys, bool enabled)
{
	if (!sys)
		return;
	sys->timing_enabled = enabled;
}

void ar_system_get_metrics(const ar_system_t *sys, ar_metrics_t *out_metrics)
{
	if (!sys || !out_metrics)
		return;
	*out_metrics = sys->metrics;

	islet_quad_map_stats_t ms;
	islet_quad_map_inspect(&sys->map, &ms);
	out_metrics->map_capacity = (uint32_t)ms.cap;
	out_metrics->map_live_entries = (uint32_t)ms.live_entries;
	out_metrics->map_tombstones = (uint32_t)ms.tombstones;
	out_metrics->map_max_cluster_len = ms.max_cluster_len;
	out_metrics->arena_bytes = sys->stats.arena_bytes;
	out_metrics->map_bytes = ms.bytes_allocated;
	out_metrics->staging_bytes = sys->staging_bytes;
}
