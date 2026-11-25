/* ========================================================================= */
/* Conservative 2D Software Hi-Z Occluder (Linear Clip-Space W Depth)        */
/* ========================================================================= */

static inline void ar_occ_mark_tile_dirty(ar_occluder_t *occ, int ti)
{
	occ->tile_dirty[ti >> 5] |= (1u << (ti & 31));
}

static inline void ar_occ_refresh_tile(ar_occluder_t *occ, int tx, int ty)
{
	int ti = ty * AR_OCC_TW + tx;
	if ((occ->tile_dirty[ti >> 5] & (1u << (ti & 31))) == 0u)
		return;
	int x0 = tx * AR_OCC_TILE;
	int y0 = ty * AR_OCC_TILE;
	int x1 = x0 + AR_OCC_TILE;
	int y1 = y0 + AR_OCC_TILE;
	if (x1 > AR_OCC_W) x1 = AR_OCC_W;
	if (y1 > AR_OCC_H) y1 = AR_OCC_H;
	float mx = 0.0f;
	for (int y = y0; y < y1; y++) {
		const float *row = &occ->depth[y * AR_OCC_W];
		for (int x = x0; x < x1; x++) {
			if (row[x] > mx) mx = row[x];
		}
	}
	occ->tile_max[ti] = mx;
	occ->tile_dirty[ti >> 5] &= ~(1u << (ti & 31));
}

void ar_occluder_reset(ar_occluder_t *occ, const float mvp[16],
                       float near_w, float depth_bias)
{
	if (!occ) return;
	if (!mvp) {
		occ->valid = false;
		return;
	}
	memcpy(occ->mvp, mvp, sizeof(occ->mvp));
	for (int i = 0; i < AR_OCC_W * AR_OCC_H; i++)
		occ->depth[i] = FLT_MAX;
	for (int i = 0; i < AR_OCC_TILES; i++)
		occ->tile_max[i] = FLT_MAX;
	occ->tile_dirty[0] = 0u;
	occ->tile_dirty[1] = 0u;
	occ->tile_dirty[2] = 0u;
	occ->near_w = (near_w > 0.0f) ? near_w : 0.5f;
	occ->depth_bias = (depth_bias > 0.0f) ? depth_bias : 0.002f;
	occ->culled_count = 0u;
	occ->tested_count = 0u;
	occ->valid = true;
}

void ar_occluder_rasterize_tri(ar_occluder_t *occ,
                               const float p0[3], const float p1[3], const float p2[3])
{
	if (!occ || !occ->valid || !p0 || !p1 || !p2) return;
	const float *m = occ->mvp;
	const float *pts[3] = { p0, p1, p2 };
	float sx[3], sy[3], cw[3];

	for (int k = 0; k < 3; k++) {
		float x = pts[k][0], y = pts[k][1], z = pts[k][2];
		float w = m[3] * x + m[7] * y + m[11] * z + m[15];
		if (w <= occ->near_w)
			return; /* Skip near/behind-camera triangle (strictly conservative) */
		float cx = m[0] * x + m[4] * y + m[8]  * z + m[12];
		float cy = m[1] * x + m[5] * y + m[9]  * z + m[13];
		float inv_w = 1.0f / w;
		sx[k] = (cx * inv_w * 0.5f + 0.5f) * (float)AR_OCC_W;
		sy[k] = (cy * inv_w * 0.5f + 0.5f) * (float)AR_OCC_H;
		cw[k] = w;
	}

	float area = (sx[1] - sx[0]) * (sy[2] - sy[0]) - (sy[1] - sy[0]) * (sx[2] - sx[0]);
	if (fabsf(area) < 1e-4f)
		return;
	if (area < 0.0f) {
		float tx = sx[1]; sx[1] = sx[2]; sx[2] = tx;
		float ty = sy[1]; sy[1] = sy[2]; sy[2] = ty;
		float tw = cw[1]; cw[1] = cw[2]; cw[2] = tw;
	}

	float min_x = sx[0]; if (sx[1] < min_x) min_x = sx[1]; if (sx[2] < min_x) min_x = sx[2];
	float max_x = sx[0]; if (sx[1] > max_x) max_x = sx[1]; if (sx[2] > max_x) max_x = sx[2];
	float min_y = sy[0]; if (sy[1] < min_y) min_y = sy[1]; if (sy[2] < min_y) min_y = sy[2];
	float max_y = sy[0]; if (sy[1] > max_y) max_y = sy[1]; if (sy[2] > max_y) max_y = sy[2];

	if (max_x < 0.0f || min_x >= (float)AR_OCC_W || max_y < 0.0f || min_y >= (float)AR_OCC_H)
		return;

	int ix0 = (int)floorf(min_x);
	int ix1 = (int)floorf(max_x);
	int iy0 = (int)floorf(min_y);
	int iy1 = (int)floorf(max_y);
	if (ix0 < 0) ix0 = 0;
	if (iy0 < 0) iy0 = 0;
	if (ix1 >= AR_OCC_W) ix1 = AR_OCC_W - 1;
	if (iy1 >= AR_OCC_H) iy1 = AR_OCC_H - 1;
	if (ix0 > ix1 || iy0 > iy1)
		return;

	float tri_w = cw[0];
	if (cw[1] > tri_w) tri_w = cw[1];
	if (cw[2] > tri_w) tri_w = cw[2];
	float w_write = tri_w * (1.0f + occ->depth_bias) + 1e-4f;

	float dx01 = sx[1] - sx[0], dy01 = sy[1] - sy[0];
	float dx12 = sx[2] - sx[1], dy12 = sy[2] - sy[1];
	float dx20 = sx[0] - sx[2], dy20 = sy[0] - sy[2];

	for (int py = iy0; py <= iy1; py++) {
		float cy_p = (float)py + 0.5f;
		float cx0 = (float)ix0 + 0.5f;
		float e0 = dx01 * (cy_p - sy[0]) - dy01 * (cx0 - sx[0]);
		float e1 = dx12 * (cy_p - sy[1]) - dy12 * (cx0 - sx[1]);
		float e2 = dx20 * (cy_p - sy[2]) - dy20 * (cx0 - sx[2]);
		float *row = &occ->depth[py * AR_OCC_W];
		uint32_t row_tx_dirty = 0u;

		for (int px = ix0; px <= ix1; px++) {
			if (e0 >= 0.0f && e1 >= 0.0f && e2 >= 0.0f) {
				if (w_write < row[px]) {
					row[px] = w_write;
					row_tx_dirty |= (1u << (px >> 3));
				}
			}
			e0 -= dy01;
			e1 -= dy12;
			e2 -= dy20;
		}
		if (row_tx_dirty) {
			int tile_row_base = (py >> 3) * AR_OCC_TW;
			while (row_tx_dirty) {
				int tx = __builtin_ctz(row_tx_dirty);
				row_tx_dirty &= (row_tx_dirty - 1u);
				ar_occ_mark_tile_dirty(occ, tile_row_base + tx);
			}
		}
	}
}

static void ar_node_compute_occ_cell_min(ar_node_t *node)
{
	if (!node || !node->values || node->grid_dim == 0) {
		if (node) node->occ_grid_k = 0;
		return;
	}
	uint32_t gd = node->grid_dim;
	uint32_t vd = node->vert_dim;
	uint32_t K = (gd < AR_OCC_CELL_K) ? gd : AR_OCC_CELL_K;
	uint32_t s[AR_OCC_CELL_K + 1u];
	for (uint32_t c = 0; c <= K; c++)
		s[c] = (c * gd) / K;
	for (uint32_t gr = 0; gr < K; gr++) {
		for (uint32_t gc = 0; gc < K; gc++) {
			int16_t mn = 32767;
			for (uint32_t vy = s[gr]; vy <= s[gr + 1]; vy++) {
				const int16_t *vrow = &node->values[vy * vd];
				for (uint32_t vx = s[gc]; vx <= s[gc + 1]; vx++) {
					if (vrow[vx] < mn) mn = vrow[vx];
				}
			}
			node->occ_cell_min[gr * K + gc] = mn;
		}
	}
	node->occ_grid_k = (uint8_t)K;
}

void ar_occluder_rasterize_node(ar_occluder_t *occ, const ar_node_t *node,
                                float value_scale, float z_bias,
                                ar_point_3d_fn point_fn, void *point_ctx)
{
	if (!occ || !occ->valid || !node || !node->values || node->grid_dim == 0)
		return;
	if ((node->flags & (AR_NODE_FLAG_SPARSE | AR_NODE_FLAG_DEFERRED)) != 0u)
		return;

	float inv_scale = 1.0f / ((value_scale > 0.0f) ? value_scale : AR_DEFAULT_VALUE_SCALE);
	uint32_t gd = node->grid_dim;
	uint32_t vd = node->vert_dim;
	uint32_t K = (gd < AR_OCC_CELL_K) ? gd : AR_OCC_CELL_K;

	uint32_t s[AR_OCC_CELL_K + 1u];
	for (uint32_t c = 0; c <= K; c++)
		s[c] = (c * gd) / K;

	int16_t cell_min_fallback[AR_OCC_CELL_K * AR_OCC_CELL_K];
	const int16_t *cell_min;
	if (node->occ_grid_k == K) {
		cell_min = node->occ_cell_min;
	} else {
		for (uint32_t gr = 0; gr < K; gr++) {
			for (uint32_t gc = 0; gc < K; gc++) {
				int16_t mn = 32767;
				for (uint32_t vy = s[gr]; vy <= s[gr + 1]; vy++) {
					const int16_t *vrow = &node->values[vy * vd];
					for (uint32_t vx = s[gc]; vx <= s[gc + 1]; vx++) {
						if (vrow[vx] < mn) mn = vrow[vx];
					}
				}
				cell_min_fallback[gr * K + gc] = mn;
			}
		}
		cell_min = cell_min_fallback;
	}

	float verts[(AR_OCC_CELL_K + 1u) * (AR_OCC_CELL_K + 1u)][3];
	float vert_wz[(AR_OCC_CELL_K + 1u) * (AR_OCC_CELL_K + 1u)];
	float inv_gd = 1.0f / (float)gd;
	for (uint32_t vr = 0; vr <= K; vr++) {
		uint32_t r0 = (vr > 0u) ? (vr - 1u) : 0u;
		uint32_t r1 = (vr < K)  ? vr : (K - 1u);
		float wy = (float)node->oy + ((float)s[vr] * inv_gd) * (float)node->h;
		for (uint32_t vc = 0; vc <= K; vc++) {
			uint32_t c0 = (vc > 0u) ? (vc - 1u) : 0u;
			uint32_t c1 = (vc < K)  ? vc : (K - 1u);
			int16_t mn = 32767;
			for (uint32_t gr = r0; gr <= r1; gr++) {
				for (uint32_t gc = c0; gc <= c1; gc++) {
					if (cell_min[gr * K + gc] < mn)
						mn = cell_min[gr * K + gc];
				}
			}
			float wx = (float)node->ox + ((float)s[vc] * inv_gd) * (float)node->w;
			float wz = (float)mn * inv_scale - z_bias;
			vert_wz[vr * (K + 1u) + vc] = wz;
			float *dst = verts[vr * (K + 1u) + vc];
			if (point_fn) {
				point_fn(point_ctx, wx, wy, wz, dst);
			} else {
				dst[0] = wx;
				dst[1] = wy;
				dst[2] = wz;
			}
		}
	}

	for (uint32_t gr = 0; gr < K; gr++) {
		for (uint32_t gc = 0; gc < K; gc++) {
			const float *v00 = verts[gr * (K + 1u) + gc];
			const float *v10 = verts[gr * (K + 1u) + gc + 1u];
			const float *v01 = verts[(gr + 1u) * (K + 1u) + gc];
			const float *v11 = verts[(gr + 1u) * (K + 1u) + gc + 1u];
			ar_occluder_rasterize_tri(occ, v00, v10, v11);
			ar_occluder_rasterize_tri(occ, v00, v11, v01);
		}
	}

	float node_span_z = (float)(node->max_val - node->min_val) * inv_scale;
	float max_dim = (node->w > node->h) ? (float)node->w : (float)node->h;
	float skirt_drop = (node_span_z > 32.0f) ? node_span_z : 32.0f;
	if (0.05f * max_dim > skirt_drop)
		skirt_drop = 0.05f * max_dim;
	if (skirt_drop > 512.0f)
		skirt_drop = 512.0f;

	float sk_a[3], sk_b[3];
	for (int edge = 0; edge < 2; edge++) {
		uint32_t vr = (edge == 0) ? 0u : K;
		float wy = (float)node->oy + ((float)s[vr] * inv_gd) * (float)node->h;
		for (uint32_t gc = 0; gc < K; gc++) {
			uint32_t i0 = vr * (K + 1u) + gc;
			uint32_t i1 = i0 + 1u;
			float wx0 = (float)node->ox + ((float)s[gc] * inv_gd) * (float)node->w;
			float wx1 = (float)node->ox + ((float)s[gc + 1u] * inv_gd) * (float)node->w;
			float wz0 = vert_wz[i0] - skirt_drop;
			float wz1 = vert_wz[i1] - skirt_drop;
			if (point_fn) {
				point_fn(point_ctx, wx0, wy, wz0, sk_a);
				point_fn(point_ctx, wx1, wy, wz1, sk_b);
			} else {
				sk_a[0] = wx0; sk_a[1] = wy; sk_a[2] = wz0;
				sk_b[0] = wx1; sk_b[1] = wy; sk_b[2] = wz1;
			}
			ar_occluder_rasterize_tri(occ, verts[i0], verts[i1], sk_b);
			ar_occluder_rasterize_tri(occ, verts[i0], sk_b, sk_a);
		}
	}
	for (int edge = 0; edge < 2; edge++) {
		uint32_t vc = (edge == 0) ? 0u : K;
		float wx = (float)node->ox + ((float)s[vc] * inv_gd) * (float)node->w;
		for (uint32_t gr = 0; gr < K; gr++) {
			uint32_t i0 = gr * (K + 1u) + vc;
			uint32_t i1 = (gr + 1u) * (K + 1u) + vc;
			float wy0 = (float)node->oy + ((float)s[gr] * inv_gd) * (float)node->h;
			float wy1 = (float)node->oy + ((float)s[gr + 1u] * inv_gd) * (float)node->h;
			float wz0 = vert_wz[i0] - skirt_drop;
			float wz1 = vert_wz[i1] - skirt_drop;
			if (point_fn) {
				point_fn(point_ctx, wx, wy0, wz0, sk_a);
				point_fn(point_ctx, wx, wy1, wz1, sk_b);
			} else {
				sk_a[0] = wx; sk_a[1] = wy0; sk_a[2] = wz0;
				sk_b[0] = wx; sk_b[1] = wy1; sk_b[2] = wz1;
			}
			ar_occluder_rasterize_tri(occ, verts[i0], verts[i1], sk_b);
			ar_occluder_rasterize_tri(occ, verts[i0], sk_b, sk_a);
		}
	}
}

bool ar_occluder_test_box(ar_occluder_t *occ,
                          const float b_min[3], const float b_max[3])
{
	if (!occ || !occ->valid || !b_min || !b_max)
		return false;

	occ->tested_count++;
	const float *m = occ->mvp;
	float min_nx =  1e30f, max_nx = -1e30f;
	float min_ny =  1e30f, max_ny = -1e30f;
	float min_w  =  FLT_MAX;

	for (int c = 0; c < 8; c++) {
		float x = (c & 1) ? b_max[0] : b_min[0];
		float y = (c & 2) ? b_max[1] : b_min[1];
		float z = (c & 4) ? b_max[2] : b_min[2];
		float w = m[3] * x + m[7] * y + m[11] * z + m[15];
		if (w <= occ->near_w)
			return false;
		float cx = m[0] * x + m[4] * y + m[8]  * z + m[12];
		float cy = m[1] * x + m[5] * y + m[9]  * z + m[13];
		float inv_w = 1.0f / w;
		float nx = cx * inv_w;
		float ny = cy * inv_w;
		if (nx < min_nx) min_nx = nx;
		if (nx > max_nx) max_nx = nx;
		if (ny < min_ny) min_ny = ny;
		if (ny > max_ny) max_ny = ny;
		if (w  < min_w)  min_w  = w;
	}

	float min_px = (min_nx * 0.5f + 0.5f) * (float)AR_OCC_W;
	float max_px = (max_nx * 0.5f + 0.5f) * (float)AR_OCC_W;
	float min_py = (min_ny * 0.5f + 0.5f) * (float)AR_OCC_H;
	float max_py = (max_ny * 0.5f + 0.5f) * (float)AR_OCC_H;

	int x0 = (int)floorf(min_px) - 1;
	int x1 = (int)ceilf(max_px) + 1;
	int y0 = (int)floorf(min_py) - 1;
	int y1 = (int)ceilf(max_py) + 1;
	if (x0 < 0) x0 = 0;
	if (y0 < 0) y0 = 0;
	if (x1 >= AR_OCC_W) x1 = AR_OCC_W - 1;
	if (y1 >= AR_OCC_H) y1 = AR_OCC_H - 1;
	if (x0 > x1 || y0 > y1)
		return false;

	int tx0 = x0 / AR_OCC_TILE;
	int tx1 = x1 / AR_OCC_TILE;
	int ty0 = y0 / AR_OCC_TILE;
	int ty1 = y1 / AR_OCC_TILE;

	for (int ty = ty0; ty <= ty1; ty++) {
		int tpy0 = ty * AR_OCC_TILE;
		int tpy1 = tpy0 + AR_OCC_TILE - 1;
		if (tpy1 >= AR_OCC_H) tpy1 = AR_OCC_H - 1;
		for (int tx = tx0; tx <= tx1; tx++) {
			ar_occ_refresh_tile(occ, tx, ty);
			float tmax = occ->tile_max[ty * AR_OCC_TW + tx];
			if (tmax < min_w)
				continue;

			int tpx0 = tx * AR_OCC_TILE;
			int tpx1 = tpx0 + AR_OCC_TILE - 1;
			if (tpx1 >= AR_OCC_W) tpx1 = AR_OCC_W - 1;

			if (tpx0 >= x0 && tpx1 <= x1 && tpy0 >= y0 && tpy1 <= y1)
				return false;

			int sx0 = (x0 > tpx0) ? x0 : tpx0;
			int sx1 = (x1 < tpx1) ? x1 : tpx1;
			int sy0 = (y0 > tpy0) ? y0 : tpy0;
			int sy1 = (y1 < tpy1) ? y1 : tpy1;
			for (int py = sy0; py <= sy1; py++) {
				const float *row = &occ->depth[py * AR_OCC_W];
				for (int px = sx0; px <= sx1; px++) {
					if (row[px] >= min_w)
						return false;
				}
			}
		}
	}

	occ->culled_count++;
	return true;
}

static inline bool ar_occ_test_padded(ar_system_t *sys, ar_occluder_t *occ,
                                      const float b_min[3], const float b_max[3],
                                      float pad, bool pad_all_axes)
{
	if (!occ || !occ->valid)
		return false;
	uint64_t t0 = ar_now_us(sys && sys->timing_enabled);
	bool culled;
	if (pad <= 0.0f) {
		culled = ar_occluder_test_box(occ, b_min, b_max);
	} else {
		float pmin[3] = { b_min[0] - (pad_all_axes ? pad : 0.0f),
		                  b_min[1] - (pad_all_axes ? pad : 0.0f),
		                  b_min[2] - pad };
		float pmax[3] = { b_max[0] + (pad_all_axes ? pad : 0.0f),
		                  b_max[1] + (pad_all_axes ? pad : 0.0f),
		                  b_max[2] + pad };
		culled = ar_occluder_test_box(occ, pmin, pmax);
	}
	if (sys) {
		sys->metrics.frame_occ_tested++;
		if (culled) sys->metrics.frame_occ_culled++;
		if (sys->timing_enabled)
			sys->metrics.us_occ_test += (uint32_t)(ar_now_us(true) - t0);
	}
	return culled;
}
