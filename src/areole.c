#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "ttypt/morton.h"
#include "ttypt/islet_quad.h"
#include "ttypt/areole.h"

/* ========================================================================= */
/* Matrix Mathematics (Zero-Dependency Pure C, Column-Major 4x4)             */
/* ========================================================================= */

void ar_mat4_identity(float out[16])
{
	memset(out, 0, 16 * sizeof(float));
	out[0] = 1.0f;
	out[5] = 1.0f;
	out[10] = 1.0f;
	out[15] = 1.0f;
}

void ar_mat4_multiply(float out[16], const float a[16], const float b[16])
{
	float tmp[16];
	for (int c = 0; c < 4; c++) {
		for (int r = 0; r < 4; r++) {
			tmp[c * 4 + r] = a[0 * 4 + r] * b[c * 4 + 0] +
			                 a[1 * 4 + r] * b[c * 4 + 1] +
			                 a[2 * 4 + r] * b[c * 4 + 2] +
			                 a[3 * 4 + r] * b[c * 4 + 3];
		}
	}
	memcpy(out, tmp, 16 * sizeof(float));
}

void ar_mat4_ortho(float out[16], float left, float right, float bottom, float top, float near_z, float far_z)
{
	memset(out, 0, 16 * sizeof(float));
	float rl = right - left;
	float tb = top - bottom;
	float fn = far_z - near_z;
	if (fabsf(rl) < 1e-7f || fabsf(tb) < 1e-7f || fabsf(fn) < 1e-7f) {
		ar_mat4_identity(out);
		return;
	}
	out[0]  = 2.0f / rl;
	out[5]  = 2.0f / tb;
	out[10] = -2.0f / fn;
	out[12] = -(right + left) / rl;
	out[13] = -(top + bottom) / tb;
	out[14] = -(far_z + near_z) / fn;
	out[15] = 1.0f;
}

void ar_mat4_perspective(float out[16], float fov_rad, float aspect, float near_z, float far_z)
{
	memset(out, 0, 16 * sizeof(float));
	float half_fov = fov_rad * 0.5f;
	float sin_fov = sinf(half_fov);
	float cos_fov = cosf(half_fov);
	if (fabsf(sin_fov) < 1e-7f || fabsf(aspect) < 1e-7f || fabsf(far_z - near_z) < 1e-7f) {
		ar_mat4_identity(out);
		return;
	}
	float f = cos_fov / sin_fov;
	float nf = 1.0f / (near_z - far_z);

	out[0]  = f / aspect;
	out[5]  = f;
	out[10] = (far_z + near_z) * nf;
	out[11] = -1.0f;
	out[14] = 2.0f * far_z * near_z * nf;
}

/* Standard right-handed OpenGL/WebGL LookAt:
 *   f = normalize(target - eye)
 *   s = normalize(f x up)  == normalize(up x z) where z = -f
 *   u = s x f              == z x s
 */
void ar_mat4_lookat(float out[16], const float eye[3], const float target[3], const float up[3])
{
	float z0 = eye[0] - target[0];
	float z1 = eye[1] - target[1];
	float z2 = eye[2] - target[2];
	float len_z = sqrtf(z0 * z0 + z1 * z1 + z2 * z2);
	if (len_z < 1e-7f) {
		z0 = 0.0f; z1 = 0.0f; z2 = 1.0f;
	} else {
		float inv_z = 1.0f / len_z;
		z0 *= inv_z; z1 *= inv_z; z2 *= inv_z;
	}

	/* Right-handed camera right vector: x = normalize(up x z) */
	float upx = up[0], upy = up[1], upz = up[2];
	float x0 = upy * z2 - upz * z1;
	float x1 = upz * z0 - upx * z2;
	float x2 = upx * z1 - upy * z0;
	float len_x = sqrtf(x0 * x0 + x1 * x1 + x2 * x2);
	if (len_x < 1e-7f) {
		x0 = 1.0f; x1 = 0.0f; x2 = 0.0f;
	} else {
		float inv_x = 1.0f / len_x;
		x0 *= inv_x; x1 *= inv_x; x2 *= inv_x;
	}

	/* Right-handed camera up vector: y = z x x */
	float y0 = z1 * x2 - z2 * x1;
	float y1 = z2 * x0 - z0 * x2;
	float y2 = z0 * x1 - z1 * x0;

	out[0] = x0; out[1] = y0; out[2] = z0; out[3] = 0.0f;
	out[4] = x1; out[5] = y1; out[6] = z1; out[7] = 0.0f;
	out[8] = x2; out[9] = y2; out[10] = z2; out[11] = 0.0f;
	out[12] = -(x0 * eye[0] + x1 * eye[1] + x2 * eye[2]);
	out[13] = -(y0 * eye[0] + y1 * eye[1] + y2 * eye[2]);
	out[14] = -(z0 * eye[0] + z1 * eye[1] + z2 * eye[2]);
	out[15] = 1.0f;
}

/* Left-handed / East-flipped LookAt convention: x = -(up x z), y = x x z */
void ar_mat4_lookat_lh(float out[16], const float eye[3], const float target[3], const float up[3])
{
	float z0 = eye[0] - target[0];
	float z1 = eye[1] - target[1];
	float z2 = eye[2] - target[2];
	float len_z = sqrtf(z0 * z0 + z1 * z1 + z2 * z2);
	if (len_z < 1e-7f) {
		z0 = 0.0f; z1 = 0.0f; z2 = 1.0f;
	} else {
		float inv_z = 1.0f / len_z;
		z0 *= inv_z; z1 *= inv_z; z2 *= inv_z;
	}

	float upx = up[0], upy = up[1], upz = up[2];
	float x0 = -(upy * z2 - upz * z1);
	float x1 = -(upz * z0 - upx * z2);
	float x2 = -(upx * z1 - upy * z0);
	float len_x = sqrtf(x0 * x0 + x1 * x1 + x2 * x2);
	if (len_x < 1e-7f) {
		x0 = 1.0f; x1 = 0.0f; x2 = 0.0f;
	} else {
		float inv_x = 1.0f / len_x;
		x0 *= inv_x; x1 *= inv_x; x2 *= inv_x;
	}

	float y0 = x1 * z2 - x2 * z1;
	float y1 = x2 * z0 - x0 * z2;
	float y2 = x0 * z1 - x1 * z0;

	out[0] = x0; out[1] = y0; out[2] = z0; out[3] = 0.0f;
	out[4] = x1; out[5] = y1; out[6] = z1; out[7] = 0.0f;
	out[8] = x2; out[9] = y2; out[10] = z2; out[11] = 0.0f;
	out[12] = -(x0 * eye[0] + x1 * eye[1] + x2 * eye[2]);
	out[13] = -(y0 * eye[0] + y1 * eye[1] + y2 * eye[2]);
	out[14] = -(z0 * eye[0] + z1 * eye[1] + z2 * eye[2]);
	out[15] = 1.0f;
}

/* 4x4 Matrix Inversion (double-precision cofactors for high near/far ratios) */
bool ar_mat4_invert(float out[16], const float mf[16])
{
	double m[16];
	for (int i = 0; i < 16; i++)
		m[i] = (double)mf[i];

	double inv[16];
	inv[0] = m[5]  * m[10] * m[15] -
	         m[5]  * m[11] * m[14] -
	         m[9]  * m[6]  * m[15] +
	         m[9]  * m[7]  * m[14] +
	         m[13] * m[6]  * m[11] -
	         m[13] * m[7]  * m[10];

	inv[4] = -m[4]  * m[10] * m[15] +
	          m[4]  * m[11] * m[14] +
	          m[8]  * m[6]  * m[15] -
	          m[8]  * m[7]  * m[14] -
	          m[12] * m[6]  * m[11] +
	          m[12] * m[7]  * m[10];

	inv[8] = m[4]  * m[9] * m[15] -
	         m[4]  * m[11] * m[13] -
	         m[8]  * m[5] * m[15] +
	         m[8]  * m[7] * m[13] +
	         m[12] * m[5] * m[11] -
	         m[12] * m[7] * m[9];

	inv[12] = -m[4]  * m[9] * m[14] +
	           m[4]  * m[10] * m[13] +
	           m[8]  * m[5] * m[14] -
	           m[8]  * m[6] * m[13] -
	           m[12] * m[5] * m[10] +
	           m[12] * m[6] * m[9];

	inv[1] = -m[1]  * m[10] * m[15] +
	          m[1]  * m[11] * m[14] +
	          m[9]  * m[2] * m[15] -
	          m[9]  * m[3] * m[14] -
	          m[13] * m[2] * m[11] +
	          m[13] * m[3] * m[10];

	inv[5] = m[0]  * m[10] * m[15] -
	         m[0]  * m[11] * m[14] -
	         m[8]  * m[2] * m[15] +
	         m[8]  * m[3] * m[14] +
	         m[12] * m[2] * m[11] -
	         m[12] * m[3] * m[10];

	inv[9] = -m[0]  * m[9] * m[15] +
	          m[0]  * m[11] * m[13] +
	          m[8]  * m[1] * m[15] -
	          m[8]  * m[3] * m[13] -
	          m[12] * m[1] * m[11] +
	          m[12] * m[3] * m[9];

	inv[13] = m[0]  * m[9] * m[14] -
	          m[0]  * m[10] * m[13] -
	          m[8]  * m[1] * m[14] +
	          m[8]  * m[2] * m[13] +
	          m[12] * m[1] * m[10] -
	          m[12] * m[2] * m[9];

	inv[2] = m[1]  * m[6] * m[15] -
	         m[1]  * m[7] * m[14] -
	         m[5]  * m[2] * m[15] +
	         m[5]  * m[3] * m[14] +
	         m[13] * m[2] * m[7] -
	         m[13] * m[3] * m[6];

	inv[6] = -m[0]  * m[6] * m[15] +
	          m[0]  * m[7] * m[14] +
	          m[4]  * m[2] * m[15] -
	          m[4]  * m[3] * m[14] -
	          m[12] * m[2] * m[7] +
	          m[12] * m[3] * m[6];

	inv[10] = m[0]  * m[5] * m[15] -
	          m[0]  * m[7] * m[13] -
	          m[4]  * m[1] * m[15] +
	          m[4]  * m[3] * m[13] +
	          m[12] * m[1] * m[7] -
	          m[12] * m[3] * m[5];

	inv[14] = -m[0]  * m[5] * m[14] +
	           m[0]  * m[6] * m[13] +
	           m[4]  * m[1] * m[14] -
	           m[4]  * m[2] * m[13] -
	           m[12] * m[1] * m[6] +
	           m[12] * m[2] * m[5];

	inv[3] = -m[1] * m[6] * m[11] +
	          m[1] * m[7] * m[10] +
	          m[5] * m[2] * m[11] -
	          m[5] * m[3] * m[10] -
	          m[9] * m[2] * m[7] +
	          m[9] * m[3] * m[6];

	inv[7] = m[0] * m[6] * m[11] -
	         m[0] * m[7] * m[10] -
	         m[4] * m[2] * m[11] +
	         m[4] * m[3] * m[10] +
	         m[8] * m[2] * m[7] -
	         m[8] * m[3] * m[6];

	inv[11] = -m[0] * m[5] * m[11] +
	           m[0] * m[7] * m[9] +
	           m[4] * m[1] * m[11] -
	           m[4] * m[3] * m[9] -
	           m[8] * m[1] * m[7] +
	           m[8] * m[3] * m[5];

	inv[15] = m[0] * m[5] * m[10] -
	          m[0] * m[6] * m[9] -
	          m[4] * m[1] * m[10] +
	          m[4] * m[2] * m[9] +
	          m[8] * m[1] * m[6] -
	          m[8] * m[2] * m[5];

	double det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];
	if (fabs(det) < 1e-35 || !isfinite(det)) {
		ar_mat4_identity(out);
		return false;
	}

	double inv_det = 1.0 / det;
	for (int i = 0; i < 16; i++) {
		out[i] = (float)(inv[i] * inv_det);
	}
	return true;
}

/* ========================================================================= */
/* Frustum Extraction & AABB Culling (Gribb-Hartmann 6-Plane Extraction)     */
/* ========================================================================= */

void ar_frustum_from_matrix(ar_frustum_t *out_frustum, const float m[16])
{
	if (!out_frustum || !m) return;

	/* Left: row3 + row0 */
	out_frustum->planes[0][0] = m[3]  + m[0];
	out_frustum->planes[0][1] = m[7]  + m[4];
	out_frustum->planes[0][2] = m[11] + m[8];
	out_frustum->planes[0][3] = m[15] + m[12];

	/* Right: row3 - row0 */
	out_frustum->planes[1][0] = m[3]  - m[0];
	out_frustum->planes[1][1] = m[7]  - m[4];
	out_frustum->planes[1][2] = m[11] - m[8];
	out_frustum->planes[1][3] = m[15] - m[12];

	/* Bottom: row3 + row1 */
	out_frustum->planes[2][0] = m[3]  + m[1];
	out_frustum->planes[2][1] = m[7]  + m[5];
	out_frustum->planes[2][2] = m[11] + m[9];
	out_frustum->planes[2][3] = m[15] + m[13];

	/* Top: row3 - row1 */
	out_frustum->planes[3][0] = m[3]  - m[1];
	out_frustum->planes[3][1] = m[7]  - m[5];
	out_frustum->planes[3][2] = m[11] - m[9];
	out_frustum->planes[3][3] = m[15] - m[13];

	/* Near: row3 + row2 (OpenGL NDC z in [-1, 1]) */
	out_frustum->planes[4][0] = m[3]  + m[2];
	out_frustum->planes[4][1] = m[7]  + m[6];
	out_frustum->planes[4][2] = m[11] + m[10];
	out_frustum->planes[4][3] = m[15] + m[14];

	/* Far: row3 - row2 */
	out_frustum->planes[5][0] = m[3]  - m[2];
	out_frustum->planes[5][1] = m[7]  - m[6];
	out_frustum->planes[5][2] = m[11] - m[10];
	out_frustum->planes[5][3] = m[15] - m[14];

	for (int i = 0; i < 6; i++) {
		float a = out_frustum->planes[i][0];
		float b = out_frustum->planes[i][1];
		float c = out_frustum->planes[i][2];
		float len = sqrtf(a * a + b * b + c * c);
		if (len > 1e-7f) {
			float inv_len = 1.0f / len;
			out_frustum->planes[i][0] *= inv_len;
			out_frustum->planes[i][1] *= inv_len;
			out_frustum->planes[i][2] *= inv_len;
			out_frustum->planes[i][3] *= inv_len;
		}
	}
}

bool ar_frustum_test_box(const ar_frustum_t *frustum,
                         float min_x, float min_y, float min_z,
                         float max_x, float max_y, float max_z)
{
	if (!frustum) return true;

	for (int i = 0; i < 6; i++) {
		float a = frustum->planes[i][0];
		float b = frustum->planes[i][1];
		float c = frustum->planes[i][2];
		float d = frustum->planes[i][3];

		float px = (a > 0.0f) ? max_x : min_x;
		float py = (b > 0.0f) ? max_y : min_y;
		float pz = (c > 0.0f) ? max_z : min_z;

		if (a * px + b * py + c * pz + d < 0.0f) {
			return false;
		}
	}
	return true;
}

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

	/* Max vertex clip-w is an exact upper bound for perspective-interpolated w on the triangle */
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
		int ty = py / AR_OCC_TILE;

		for (int px = ix0; px <= ix1; px++) {
			if (e0 >= 0.0f && e1 >= 0.0f && e2 >= 0.0f) {
				if (w_write < row[px]) {
					row[px] = w_write;
					ar_occ_mark_tile_dirty(occ, ty * AR_OCC_TW + (px / AR_OCC_TILE));
				}
			}
			e0 -= dy01;
			e1 -= dy12;
			e2 -= dy20;
		}
	}
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
	uint32_t K = (gd < 8u) ? gd : 8u;

	uint32_t s[9];
	for (uint32_t c = 0; c <= K; c++)
		s[c] = (c * gd) / K;

	/* Step 1: exact minimum sample inside each coarse cell [s[gc]..s[gc+1]] x [s[gr]..s[gr+1]] */
	int16_t cell_min[64];
	for (uint32_t gr = 0; gr < K; gr++) {
		for (uint32_t gc = 0; gc < K; gc++) {
			int16_t mn = 32767;
			for (uint32_t vy = s[gr]; vy <= s[gr + 1]; vy++) {
				const int16_t *vrow = &node->values[vy * vd];
				for (uint32_t vx = s[gc]; vx <= s[gc + 1]; vx++) {
					if (vrow[vx] < mn) mn = vrow[vx];
				}
			}
			cell_min[gr * K + gc] = mn;
		}
	}

	/* Step 2: each coarse vertex takes min of all incident coarse cells,
	 * guaranteeing the entire piecewise-linear watertight coarse mesh lies <= H(x, y) everywhere */
	float verts[81][3];
	float vert_wz[81];
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

	/* Step 3: rasterize watertight K x K lowered mesh */
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

	/* Step 4: conservative downward skirts along all 4 outer borders.
	 * Because every boundary vertex already satisfies wz <= H(x, y) - z_bias, dropping
	 * straight down along the seam lies strictly inside solid terrain volume {(x,y,z): z <= H(x,y)}
	 * and seals inter-chunk elevation steps and cross-LOD T-junction cracks without false occlusion. */
	float node_span_z = (float)(node->max_val - node->min_val) * inv_scale;
	float max_dim = (node->w > node->h) ? (float)node->w : (float)node->h;
	float skirt_drop = (node_span_z > 32.0f) ? node_span_z : 32.0f;
	if (0.05f * max_dim > skirt_drop)
		skirt_drop = 0.05f * max_dim;
	if (skirt_drop > 512.0f)
		skirt_drop = 512.0f;

	float sk_a[3], sk_b[3];
	/* Top (vr = 0) and Bottom (vr = K) edges */
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
	/* Left (vc = 0) and Right (vc = K) edges */
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
			return false; /* Straddling or behind near plane -> never occluded */
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

	/* Outer-conservative 1-pixel expansion on all 4 sides (zero silhouette leaks) */
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
				continue; /* Entire 8x8 tile is strictly nearer than box */

			int tpx0 = tx * AR_OCC_TILE;
			int tpx1 = tpx0 + AR_OCC_TILE - 1;
			if (tpx1 >= AR_OCC_W) tpx1 = AR_OCC_W - 1;

			/* If query rect covers the entire tile, tmax >= min_w proves at least one pixel is visible */
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

/* ========================================================================= */
/* Vertical Bounds Envelope & Screen Projection Helpers                      */
/* ========================================================================= */

static inline void ar_compute_bounds_envelope(const ar_system_t *sys,
                                              const ar_cull_opts_t *opts,
                                              float ox, float oy, float w, float h,
                                              float *out_z_min, float *out_z_max);

static inline float ar_projected_screen_width_box_3d(const float mvp[16],
                                                     const ar_cull_opts_t *opts,
                                                     const float center_xyz[3],
                                                     const float eye_xyz[3], bool has_eye,
                                                     float span_w, float span_h,
                                                     float vp_w, float vp_h)
{
	bool is_perspective = (fabsf(mvp[3]) > 1e-5f || fabsf(mvp[7]) > 1e-5f || fabsf(mvp[11]) > 1e-5f);

	float len_row0 = sqrtf(mvp[0]*mvp[0] + mvp[4]*mvp[4] + mvp[8]*mvp[8]);
	float len_row1 = sqrtf(mvp[1]*mvp[1] + mvp[5]*mvp[5] + mvp[9]*mvp[9]);
	if (len_row0 < 1e-7f) len_row0 = 0.001f;
	if (len_row1 < 1e-7f) len_row1 = 0.001f;

	float w_eff = 1.0f;
	if (is_perspective) {
		float cw = mvp[3] * center_xyz[0] +
		           mvp[7] * center_xyz[1] +
		           mvp[11] * center_xyz[2] +
		           mvp[15];
		float clip_w_thresh = (opts && opts->behind_cam_clip_w > 0.0f)
		                      ? opts->behind_cam_clip_w : AR_DEFAULT_BEHIND_CAM_CLIP_W;
		if (has_eye) {
			float dx = center_xyz[0] - eye_xyz[0];
			float dy = center_xyz[1] - eye_xyz[1];
			float dz = center_xyz[2] - eye_xyz[2];
			float len_row3 = sqrtf(mvp[3]*mvp[3] + mvp[7]*mvp[7] + mvp[11]*mvp[11]);
			float d_center = sqrtf(dx * dx + dy * dy + dz * dz) * len_row3;
			w_eff = (d_center > clip_w_thresh) ? d_center : clip_w_thresh;
		} else {
			w_eff = (cw > clip_w_thresh) ? cw : clip_w_thresh;
		}
	}

	float px_w = (span_w * len_row0 * vp_w * 0.5f) / w_eff;
	float px_h = (span_h * len_row1 * vp_h * 0.5f) / w_eff;

	return (px_w > px_h) ? px_w : px_h;
}

/* ========================================================================= */
/* System Lifecycle & Single-Arena Buffer Management                         */
/* ========================================================================= */

typedef struct {
	int32_t ox, oy;
	uint8_t level;
	float   proj_w;
	float   dist_sq;
	float   b_min[3];
	float   b_max[3];
	bool    active;
	bool    split_blocked;
} ar_candidate_t;

typedef struct {
	uint64_t       morton_key;
	uint32_t       last_used_tick;
	uint32_t       last_touch_frame;
	uint32_t       active_view_mask;
	uint8_t        state;
#ifdef TTYPT_WQ_H
	stl_wq_task_t  async_task;
	struct ar_system *async_sys;
	uint32_t       async_slot;
	uint64_t       async_key;
	uint16_t       async_world;
	uint8_t        async_level;
	int32_t        async_ox;
	int32_t        async_oy;
	int32_t        async_w;
	int32_t        async_h;
	uint32_t       async_flags;
	int16_t        async_min_val;
	int16_t        async_max_val;
	int16_t       *staging_values;
	uint32_t      *staging_rgba;
	ar_sample_t   *staging_samples;
#endif
} ar_node_meta_t;

struct ar_system {
	int32_t           base_w;
	int32_t           base_h;
	uint16_t          grid_dim;
	uint16_t          vert_dim;
	uint32_t          vert_count;
	uint16_t          pool_capacity;
	uint16_t          max_active;
	uint8_t           max_levels;
	uint32_t          default_rgba;
	float             z_min;
	float             z_max;
	float             target_quad_px;
	float             margin;
	float             value_scale;
	uint16_t          max_root_cells;
	uint8_t           probe_ancestor_cap;
	int32_t           period_x;
	int32_t           period_y;
	ar_sampler_fn     sampler;
	ar_row_sampler_fn row_sampler;
	ar_sparsity_fn    sparsity;
	void             *user_ctx;
	uint32_t          max_populates_per_frame;
	uint32_t          frame_populates_count;
	uint32_t          map_capacity;
	ar_interp_mode_t  interp_mode;
	ar_interp_fn      interp_fn;
	float             default_pitch;
	float             default_cam_dist;
	float             default_topdown_dist;
	float             default_viewport_w;
	float             default_viewport_h;
	float             sag_radius;
	float             sag_tangent_slope;
	float             split_hysteresis;
	ar_bounds_fn      bounds_fn;
	ar_bounds_3d_fn   bounds_3d_fn;
	ar_cull_3d_fn     cull_3d_fn;
	void             *bounds_ctx;
	uint8_t           attr_count;
	ar_interp_mode_t  attr_interp_mode;
	bool              interpolate_attrs;
	uint32_t          current_tick;
	struct stl_wq    *wq;

	/* Pointers into single contiguous allocation */
	ar_node_t        *pool;
	ar_node_meta_t   *pool_meta;
	const ar_node_t **active_chunks;
	ar_candidate_t   *cands_scratch;
	ar_candidate_t   *active_cands_scratch;
	size_t            active_count;
	uint8_t          *staging_arena;

	/* Backing libislet in-memory quadtree spatial database */
	islet_quad_map_t  map;

	/* Telemetry & diagnostics */
	ar_stats_t        stats;
};

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

static inline void ar_apply_cached_z_bounds(const ar_system_t *sys,
                                            const ar_cull_opts_t *opts,
                                            uint16_t world, uint8_t level,
                                            int32_t ox, int32_t oy,
                                            float *inout_z_min, float *inout_z_max)
{
	if (!opts || !opts->use_cached_bounds || level > opts->cached_bounds_max_level)
		return;
	uint64_t key = morton_quad_key(world, level, ox, oy, sys->base_w, sys->base_h);
	uint32_t slot = 0;
	if (!islet_quad_map_get(&sys->map, key, &slot) || slot >= sys->pool_capacity)
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

static inline void ar_compute_node_3d_box(const ar_system_t *sys,
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

		/* Distance squared from true camera eye (or origin if unavailable) to 3D AABB */
		float ex = has_eye ? eye_xyz[0] : 0.0f;
		float ey = has_eye ? eye_xyz[1] : 0.0f;
		float ez = has_eye ? eye_xyz[2] : 0.0f;
		float d2 = 0.0f;
		float dx = 0.0f, dy = 0.0f, dz = 0.0f;
		if (ex < out_min[0]) dx = out_min[0] - ex; else if (ex > out_max[0]) dx = ex - out_max[0];
		if (ey < out_min[1]) dy = out_min[1] - ey; else if (ey > out_max[1]) dy = ey - out_max[1];
		if (ez < out_min[2]) dz = out_min[2] - ez; else if (ez > out_max[2]) dz = ez - out_max[2];
		d2 = dx * dx + dy * dy + dz * dz;
		*out_dist_sq = d2;
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

static inline size_t ar_align8(size_t n)
{
	return (n + 7u) & ~(size_t)7u;
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
	uint16_t max_act = cfg->max_active;
	uint8_t max_lvls = cfg->max_levels;
	uint32_t def_rgba = cfg->default_rgba;

	size_t sys_bytes = ar_align8(sizeof(ar_system_t));
	size_t pool_bytes = ar_align8((size_t)pool_cap * sizeof(ar_node_t));
	size_t pool_meta_bytes = ar_align8((size_t)pool_cap * sizeof(ar_node_meta_t));
	size_t active_bytes = ar_align8((size_t)max_act * sizeof(ar_node_t *));
	size_t cands_bytes = ar_align8((size_t)(max_act + 8u) * sizeof(ar_candidate_t));
	size_t active_cands_bytes = ar_align8((size_t)max_act * sizeof(ar_candidate_t));

	size_t per_node_vals = ar_align8((size_t)vert_count * sizeof(int16_t));
	size_t per_node_rgba = ar_align8((size_t)vert_count * sizeof(uint32_t));
	size_t per_node_samples = ar_align8((size_t)vert_count * sizeof(ar_sample_t));
	size_t per_node_arena = per_node_vals + per_node_rgba + per_node_samples;

	size_t total_bytes = sys_bytes + pool_bytes + pool_meta_bytes + active_bytes +
	                     cands_bytes + active_cands_bytes + (size_t)pool_cap * per_node_arena;

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

	for (uint16_t i = 0; i < pool_cap; i++) {
		sys->pool[i].grid_dim = grid_dim;
		sys->pool[i].vert_dim = vert_dim;
		sys->pool[i].vert_count = vert_count;
		sys->pool[i].values = (const int16_t *)cursor;
		cursor += per_node_vals;
		sys->pool[i].rgba = (const uint32_t *)cursor;
		cursor += per_node_rgba;
		sys->pool[i].samples = (const ar_sample_t *)cursor;
		cursor += per_node_samples;
	}

	float val_scale = (cfg->value_scale > 0.0f) ? cfg->value_scale : AR_DEFAULT_VALUE_SCALE;

	sys->base_w = cfg->base_w;
	sys->base_h = cfg->base_h;
	sys->grid_dim = grid_dim;
	sys->vert_dim = vert_dim;
	sys->vert_count = vert_count;
	sys->pool_capacity = pool_cap;
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
	sys->wq = cfg->wq;
	sys->staging_arena = NULL;

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
__attribute__((weak)) uint32_t stl_cpu_count(void) { return 1u; }
__attribute__((weak)) stl_wq_t *stl_wq_create(const char *name, uint32_t nthreads, uint32_t capacity)
{
	(void)name; (void)nthreads; (void)capacity;
	return NULL;
}
__attribute__((weak)) int stl_wq_submit(stl_wq_t *wq, stl_wq_task_t *task, stl_wq_fn fn, void *arg)
{
	(void)wq;
	if (!task || !fn) return -1;
	task->fn = fn;
	task->arg = arg;
	task->state = STL_WQ_RUNNING;
	fn(arg);
	task->state = STL_WQ_DONE;
	return 0;
}
__attribute__((weak)) uint32_t stl_wq_step(stl_wq_t *wq, uint32_t max_tasks)
{
	(void)wq; (void)max_tasks;
	return 0u;
}
__attribute__((weak)) void stl_wq_wait(stl_wq_t *wq, stl_wq_task_t *task)
{
	(void)wq;
	if (task) task->state = STL_WQ_DONE;
}
__attribute__((weak)) void stl_wq_drain(stl_wq_t *wq)
{
	(void)wq;
}
__attribute__((weak)) void stl_wq_destroy(stl_wq_t *wq)
{
	(void)wq;
}
#endif

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

	int16_t min_val = 32767;
	int16_t max_val = -32768;
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

void ar_system_set_wq(ar_system_t *sys, struct stl_wq *wq)
{
	if (!sys)
		return;
#ifdef TTYPT_WQ_H
	if (sys->wq && sys->wq != wq) {
		stl_wq_drain(sys->wq);
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
				stl_wq_wait(sys->wq, &sys->pool_meta[i].async_task);
		}
	}
#endif
	free(sys->staging_arena);
	islet_quad_map_destroy(&sys->map);
	free(sys);
}

void ar_system_clear(ar_system_t *sys)
{
	if (!sys)
		return;
#ifdef TTYPT_WQ_H
	if (sys->wq) {
		for (size_t i = 0; i < sys->pool_capacity; i++) {
			if (stl_wq_busy(&sys->pool_meta[i].async_task))
				stl_wq_wait(sys->wq, &sys->pool_meta[i].async_task);
			sys->pool_meta[i].async_task.state = STL_WQ_IDLE;
		}
	}
#endif
	for (size_t i = 0; i < sys->pool_capacity; i++) {
		sys->pool_meta[i].state = 0;
		sys->pool_meta[i].last_used_tick = 0;
		sys->pool_meta[i].active_view_mask = 0;
	}
	islet_quad_map_clear(&sys->map);
	sys->active_count = 0;
	sys->stats.pool_in_use = 0;
	sys->stats.active_chunks_last_frame = 0;
}

void ar_system_begin_frame(ar_system_t *sys)
{
	if (!sys)
		return;
	sys->current_tick++;
	sys->frame_populates_count = 0;
#ifdef TTYPT_WQ_H
	if (sys->wq && sys->staging_arena) {
		/* Cooperatively step up to 16 background chunks when running without worker threads */
		stl_wq_step(sys->wq, 16u);
		size_t vals_bytes = (size_t)sys->vert_count * sizeof(int16_t);
		size_t rgba_bytes = (size_t)sys->vert_count * sizeof(uint32_t);
		size_t samp_bytes = (size_t)sys->vert_count * sizeof(ar_sample_t);
		for (size_t i = 0; i < sys->pool_capacity; i++) {
			ar_node_meta_t *meta = &sys->pool_meta[i];
			meta->active_view_mask = 0;
			if (stl_wq_done(&meta->async_task)) {
				ar_node_t *node = &sys->pool[i];
				if (meta->state == 1 && meta->morton_key == meta->async_key &&
				    (node->flags & AR_NODE_FLAG_DEFERRED)) {
					memcpy((int16_t *)node->values, meta->staging_values, vals_bytes);
					memcpy((uint32_t *)node->rgba, meta->staging_rgba, rgba_bytes);
					memcpy((ar_sample_t *)node->samples, meta->staging_samples, samp_bytes);
					node->min_val = meta->async_min_val;
					node->max_val = meta->async_max_val;
					node->flags = meta->async_flags;
					node->rev++;
					sys->stats.sampler_invocations_total += sys->vert_count;
				}
				meta->async_task.state = STL_WQ_IDLE;
			}
		}
		return;
	}
#endif
	for (size_t i = 0; i < sys->pool_capacity; i++) {
		sys->pool_meta[i].active_view_mask = 0;
	}
}

static ar_node_t *ar_get_or_populate(ar_system_t *sys, uint16_t world,
                                     uint8_t level, int32_t ox, int32_t oy)
{
	uint64_t key = morton_quad_key(world, level, ox, oy, sys->base_w, sys->base_h);

	/* 1. Fast O(1) Islet spatial map lookup */
	uint32_t hit_slot = 0;
	int slot = -1;
	bool upgrading = false;
	if (islet_quad_map_get(&sys->map, key, &hit_slot)) {
		if (hit_slot < sys->pool_capacity && sys->pool_meta[hit_slot].state == 1) {
			ar_node_t *node = &sys->pool[hit_slot];
			sys->pool_meta[hit_slot].last_used_tick = sys->current_tick;
			sys->pool_meta[hit_slot].last_touch_frame = sys->current_tick;
			sys->stats.cache_hits_total++;
#ifdef TTYPT_WQ_H
			if ((node->flags & AR_NODE_FLAG_DEFERRED) &&
			    sys->wq && sys->staging_arena &&
			    stl_wq_busy(&sys->pool_meta[hit_slot].async_task)) {
				return node;
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

	sys->stats.cache_misses_total++;

	/* 2. Find an empty slot or evict oldest LRU slot not active this frame */
	if (!upgrading) {
		for (size_t i = 0; i < sys->pool_capacity; i++) {
			if (sys->pool_meta[i].state == 0) {
				slot = (int)i;
				break;
			}
		}

		if (slot < 0) {
			uint32_t oldest_tick = 0xFFFFFFFFu;
			for (size_t i = 0; i < sys->pool_capacity; i++) {
#ifdef TTYPT_WQ_H
				if (stl_wq_busy(&sys->pool_meta[i].async_task))
					continue;
#endif
				if (sys->pool_meta[i].active_view_mask == 0) {
					if (sys->pool_meta[i].last_used_tick < oldest_tick) {
						oldest_tick = sys->pool_meta[i].last_used_tick;
						slot = (int)i;
					}
				}
			}
			if (slot < 0) {
				oldest_tick = 0xFFFFFFFFu;
				for (size_t i = 0; i < sys->pool_capacity; i++) {
#ifdef TTYPT_WQ_H
					if (stl_wq_busy(&sys->pool_meta[i].async_task))
						continue;
#endif
					if (sys->pool_meta[i].last_used_tick < oldest_tick) {
						oldest_tick = sys->pool_meta[i].last_used_tick;
						slot = (int)i;
					}
				}
			}
			if (slot >= 0 && sys->pool_meta[slot].state == 1) {
				islet_quad_map_del(&sys->map, sys->pool_meta[slot].morton_key);
			}
		}
	}

	if (slot < 0)
		return NULL;

	ar_node_t *node = &sys->pool[slot];
	ar_node_meta_t *meta = &sys->pool_meta[slot];
#ifdef TTYPT_WQ_H
	if (sys->wq && stl_wq_busy(&meta->async_task)) {
		stl_wq_wait(sys->wq, &meta->async_task);
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
	node->min_val = 32767;
	node->max_val = -32768;
	meta->last_used_tick = sys->current_tick;
	meta->last_touch_frame = sys->current_tick;
	meta->active_view_mask = 0;
	meta->state = 1;
	meta->morton_key = key;

	if (!upgrading) {
		islet_quad_map_put(&sys->map, key, (uint32_t)slot);
	}

	/* 3. Query sparsity if provider supplied */
	if (sys->sparsity) {
		if (sys->sparsity(sys->user_ctx, world, level, ox, oy, node->w, node->h)) {
			node->flags |= AR_NODE_FLAG_SPARSE;
		}
	}

	/* 4. Evaluate sample grid: (grid_dim + 1) x (grid_dim + 1) vertices.
	 * Reuse already-cached parent or child nodes at even grid points to separate
	 * field computation caching from per-view rendering. */
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
			if (islet_quad_map_get(&sys->map, pkey, &pslot) &&
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
				if (islet_quad_map_get(&sys->map, ckey, &cslot) &&
				    cslot < sys->pool_capacity && sys->pool_meta[cslot].state == 1 &&
				    sys->pool[cslot].world == world && sys->pool[cslot].level == level - 1 &&
				    (sys->pool[cslot].flags & AR_NODE_FLAG_DEFERRED) == 0) {
					child_nodes[q] = &sys->pool[cslot];
				}
			}
		}
	}

	int16_t     *mut_values  = (int16_t *)node->values;
	uint32_t    *mut_rgba    = (uint32_t *)node->rgba;
	ar_sample_t *mut_samples = (ar_sample_t *)node->samples;

	/* If per-frame populate budget is exceeded and a populated parent quadtree node exists,
	 * synthesize node geometry/texture via O(1) integer bilinear upsampling without stalling */
	bool budget_exceeded = (sys->max_populates_per_frame > 0 &&
	                        sys->frame_populates_count >= sys->max_populates_per_frame);
	if (budget_exceeded && parent_node) {
		size_t idx = 0;
		int vdim_p = (int32_t)parent_node->vert_dim;
		const ar_sample_t *ps = parent_node->samples;
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

				ar_sample_t *s = &mut_samples[idx];
				if (r_even && c_even) {
					*s = ps[pr0 * vdim_p + pc0];
				} else if (r_even) {
					const ar_sample_t *s0 = &ps[pr0 * vdim_p + pc0];
					const ar_sample_t *s1 = &ps[pr0 * vdim_p + pc1];
					s->val = (int16_t)(((int32_t)s0->val + (int32_t)s1->val) >> 1);
					s->flags = s0->flags | s1->flags;
					s->rgba = (((s0->rgba & 0xFEFEFEFEu) >> 1) + ((s1->rgba & 0xFEFEFEFEu) >> 1));
					for (uint8_t a = 0; a < AR_SAMPLE_ATTR_CAP; a++) {
						s->attr[a] = (uint8_t)(((uint16_t)s0->attr[a] + (uint16_t)s1->attr[a]) >> 1);
					}
				} else if (c_even) {
					const ar_sample_t *s0 = &ps[pr0 * vdim_p + pc0];
					const ar_sample_t *s1 = &ps[pr1 * vdim_p + pc0];
					s->val = (int16_t)(((int32_t)s0->val + (int32_t)s1->val) >> 1);
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
					s->val = (int16_t)(((int32_t)s00->val + (int32_t)s10->val + (int32_t)s01->val + (int32_t)s11->val) >> 2);
					s->flags = s00->flags | s10->flags | s01->flags | s11->flags;
					s->rgba = (((s00->rgba & 0xFCFCFCFCu) >> 2) + ((s10->rgba & 0xFCFCFCFCu) >> 2) +
					           ((s01->rgba & 0xFCFCFCFCu) >> 2) + ((s11->rgba & 0xFCFCFCFCu) >> 2));
					for (uint8_t a = 0; a < AR_SAMPLE_ATTR_CAP; a++) {
						s->attr[a] = (uint8_t)(((uint16_t)s00->attr[a] + (uint16_t)s10->attr[a] +
						                        (uint16_t)s01->attr[a] + (uint16_t)s11->attr[a]) >> 2);
					}
				}
				if (s->val < node->min_val) node->min_val = s->val;
				if (s->val > node->max_val) node->max_val = s->val;
				node->flags |= (s->flags & AR_NODE_USER_FLAGS_MASK);

				mut_values[idx] = s->val;
				mut_rgba[idx] = s->rgba ? s->rgba : sys->default_rgba;
				idx++;
			}
		}
		node->flags |= AR_NODE_FLAG_DEFERRED;
		node->rev++;
#ifdef TTYPT_WQ_H
		if (sys->wq && sys->staging_arena && !stl_wq_busy(&meta->async_task)) {
			meta->async_sys = sys;
			meta->async_slot = (uint32_t)slot;
			meta->async_key = key;
			meta->async_world = world;
			meta->async_level = level;
			meta->async_ox = ox;
			meta->async_oy = oy;
			meta->async_w = node->w;
			meta->async_h = node->h;
			(void)stl_wq_submit(sys->wq, &meta->async_task, ar_async_populate_worker, meta);
		}
#endif
		return node;
	}

	sys->frame_populates_count++;
	node->flags &= ~AR_NODE_FLAG_DEFERRED;

	bool has_reusable_nodes = (parent_node != NULL ||
	                           child_nodes[0] != NULL || child_nodes[1] != NULL ||
	                           child_nodes[2] != NULL || child_nodes[3] != NULL);

	size_t idx = 0;
	if (sys->row_sampler && !has_reusable_nodes) {
		for (int r = 0; r <= gdim; r++) {
			int32_t wy = oy + (int32_t)(((int64_t)r * node->h) / gdim);
			sys->row_sampler(sys->user_ctx, world, level, wy, ox, node->w, (uint32_t)vdim, &mut_samples[idx]);
			sys->stats.sampler_invocations_total += (uint32_t)vdim;
			for (int c = 0; c <= gdim; c++) {
				ar_sample_t *s = &mut_samples[idx];
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
			for (int c = 0; c <= gdim; c++) {
				int32_t wx = ox + (int32_t)(((int64_t)c * node->w) / gdim);
				ar_sample_t *s = &mut_samples[idx];
				bool reused = false;

				/* Check cached finer child first (every grid point (c, r) in parent maps to
				 * even grid point (2*(c - qx*half), 2*(r - qy*half)) in child quadrant q) */
				if (can_reuse_hierarchy) {
					int qx = (c > half_gdim) ? 1 : 0;
					int qy = (r > half_gdim) ? 1 : 0;
					uint8_t q = (uint8_t)(qx | (qy << 1));
					const ar_node_t *cn = child_nodes[q];
					if (cn) {
						int cc = (c - (qx ? half_gdim : 0)) * 2;
						int cr = (r - (qy ? half_gdim : 0)) * 2;
						int32_t c_wx = cn->ox + (int32_t)(((int64_t)cc * cn->w) / gdim);
						int32_t c_wy = cn->oy + (int32_t)(((int64_t)cr * cn->h) / gdim);
						if (c_wx == wx && c_wy == wy) {
							*s = cn->samples[cr * vdim + cc];
							reused = true;
						}
					}
					/* Otherwise check cached parent at even grid coordinates */
					if (!reused && parent_node && ((r & 1) == 0) && ((c & 1) == 0)) {
						int pr = parent_r_off + (r >> 1);
						int pc = parent_c_off + (c >> 1);
						int32_t p_wx = parent_node->ox + (int32_t)(((int64_t)pc * parent_node->w) / gdim);
						int32_t p_wy = parent_node->oy + (int32_t)(((int64_t)pr * parent_node->h) / gdim);
						if (p_wx == wx && p_wy == wy) {
							*s = parent_node->samples[pr * vdim + pc];
							reused = true;
						}
					}
				}

				if (!reused) {
					sys->sampler(sys->user_ctx, world, level, wx, wy, s);
					sys->stats.sampler_invocations_total++;
				}

				if (s->val < node->min_val) node->min_val = s->val;
				if (s->val > node->max_val) node->max_val = s->val;
				node->flags |= (s->flags & AR_NODE_USER_FLAGS_MASK);

				mut_values[idx] = s->val;
				mut_rgba[idx] = s->rgba ? s->rgba : sys->default_rgba;
				idx++;
			}
		}
	}
	node->rev++;

	return node;
}

/* ========================================================================= */
/* Matrix Traversal & Frustum Quadtree Traversal                             */
/* ========================================================================= */

static inline float ar_dist_to_box_sq(float fx, float fy, float ox, float oy, float w, float h)
{
	float dx = 0.0f;
	if (fx < ox) dx = ox - fx;
	else if (fx > ox + w) dx = fx - (ox + w);

	float dy = 0.0f;
	if (fy < oy) dy = oy - fy;
	else if (fy > oy + h) dy = fy - (oy + h);

	return dx * dx + dy * dy;
}

static inline bool ar_unproject_ndc(const float inv_mvp[16], float ndc_x, float ndc_y, float ndc_z, float out_world[3])
{
	float v[4] = { ndc_x, ndc_y, ndc_z, 1.0f };
	float p[4];
	for (int i = 0; i < 4; i++) {
		p[i] = inv_mvp[0 * 4 + i] * v[0] +
		       inv_mvp[1 * 4 + i] * v[1] +
		       inv_mvp[2 * 4 + i] * v[2] +
		       inv_mvp[3 * 4 + i] * v[3];
	}
	if (fabsf(p[3]) < 1e-9f) return false;
	float inv_w = 1.0f / p[3];
	out_world[0] = p[0] * inv_w;
	out_world[1] = p[1] * inv_w;
	out_world[2] = p[2] * inv_w;
	return true;
}

static inline int32_t ar_floor_i32(float v)
{
	if (v < -1e9f) return -1000000000;
	if (v >  1e9f) return  1000000000;
	return (int32_t)floorf(v);
}

static inline bool ar_occ_test_padded(ar_occluder_t *occ,
                                      const float b_min[3], const float b_max[3],
                                      float pad, bool pad_all_axes)
{
	if (!occ || !occ->valid)
		return false;
	if (pad <= 0.0f)
		return ar_occluder_test_box(occ, b_min, b_max);
	float pmin[3] = { b_min[0] - (pad_all_axes ? pad : 0.0f),
	                  b_min[1] - (pad_all_axes ? pad : 0.0f),
	                  b_min[2] - pad };
	float pmax[3] = { b_max[0] + (pad_all_axes ? pad : 0.0f),
	                  b_max[1] + (pad_all_axes ? pad : 0.0f),
	                  b_max[2] + pad };
	return ar_occluder_test_box(occ, pmin, pmax);
}

static inline bool ar_frustum_test_box_sat(const ar_frustum_t *f,
                                           bool has_faabb,
                                           const float f_min[3], const float f_max[3],
                                           const float b_min[3], const float b_max[3])
{
	if (has_faabb) {
		if (b_max[0] < f_min[0] || b_min[0] > f_max[0] ||
		    b_max[1] < f_min[1] || b_min[1] > f_max[1] ||
		    b_max[2] < f_min[2] || b_min[2] > f_max[2])
			return false;
	}
	return ar_frustum_test_box(f, b_min[0], b_min[1], b_min[2],
	                           b_max[0], b_max[1], b_max[2]);
}

size_t ar_system_query_matrix(ar_system_t *sys, uint16_t world,
                              const float mvp[16],
                              const ar_cull_opts_t *opts,
                              const ar_node_t **out_nodes,
                              size_t max_out)
{
	if (!sys || !mvp || !out_nodes || max_out == 0)
		return 0;

	sys->current_tick++;
	size_t out_count = 0;
	size_t max_candidates = (max_out < sys->max_active) ? max_out : sys->max_active;

	/* 1. Extract 6 frustum planes */
	ar_frustum_t frustum;
	ar_frustum_from_matrix(&frustum, mvp);

	/* 2. Analytical Inversion for Viewport Ray Unprojection & Camera Eye Extraction */
	float inv_mvp[16];
	bool has_inv = ar_mat4_invert(inv_mvp, mvp);

	/* 3. Determine camera eye & root footprint bounds */
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

	/* Extract true camera eye position in MVP world space:
	 * For any perspective P*V, P*(0,0,0,1)^T = (0,0,P[14],0)^T -> eye = inv_mvp.col2 / inv_mvp[11]. */
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

	/* Unproject all 8 NDC frustum corners once to build exact frustum AABB (SAT culling)
	 * and exact Frustum intersect Slab footprint for flat worlds. */
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

			/* Exact convex Frustum intersect Slab [slab_z_min, slab_z_max]:
			 * Clip the 12 frustum edges against the slab planes. */
			if (all_8_ok) {
				float slab_z_min = z_min - margin;
				float slab_z_max = z_max + margin;

				if (sag_r > 0.0f) {
					/* Compute conservative sag lower bound across actual frustum reach */
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
					{0, 1}, {1, 3}, {3, 2}, {2, 0}, /* Near face */
					{4, 5}, {5, 7}, {7, 6}, {6, 4}, /* Far face */
					{0, 4}, {1, 5}, {2, 6}, {3, 7}  /* Side edges */
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
			/* Frustum lies completely outside the terrain Z slab (e.g. looking straight into sky) */
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

	if (!frustum_misses_slab) {
		/* 4. Determine root level that naturally covers the bounding box (cols * rows <= max_root_cells) */
		uint8_t root_level = 0;
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

		/* Cap grid span to at most one full world period so periodic/toroidal worlds
		 * never spawn duplicate wrapped chunks overlapping in 3D. */
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
				ar_compute_node_3d_box(sys, opts, world, root_level,
				                       rox, roy, root_w, root_h,
				                       eye_xyz, has_eye,
				                       b_min, b_max, center_xyz, &span_w, &span_h, &dist_sq,
				                       &nz_min, &nz_max);

				if (ar_frustum_test_box_sat(&frustum, has_faabb, faabb_min, faabb_max,
				                            b_min, b_max) &&
				    (!c3d || c3d(c3d_ctx, (float)rox, (float)roy, (float)root_w, (float)root_h,
				                 nz_min, nz_max, center_xyz, span_w, span_h)) &&
				    !ar_occ_test_padded(occ, b_min, b_max, occ_box_pad, occ_pad_all)) {
					float pw = ar_projected_screen_width_box_3d(mvp, opts, center_xyz,
					                                            eye_xyz, has_eye,
					                                            span_w, span_h, vp_w, vp_h);
					size_t dst_idx = cand_count;
					if (cand_count < max_candidates) {
						cand_count++;
					} else {
						/* Keep closest roots if max_candidates is smaller than visible root count */
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
	}

	/* 5. Priority Refinement Loop (with optional split/merge hysteresis and per-query populate budget) */
	float target_quad_px = (opts && opts->target_quad_px > 0.0f) ? opts->target_quad_px : sys->target_quad_px;
	float split_threshold = target_quad_px * (float)sys->grid_dim;
	float hyst = (opts && opts->split_hysteresis > 0.0f) ? opts->split_hysteresis : sys->split_hysteresis;
	if (hyst > AR_DEFAULT_MAX_HYSTERESIS) hyst = AR_DEFAULT_MAX_HYSTERESIS;
	uint32_t view_mask = (opts && opts->view_id < 32) ? (1u << opts->view_id) : 1u;
	uint8_t min_lvl = opts ? opts->min_level : 0;
	uint16_t max_pops = opts ? opts->max_populates_per_query : 0;
	uint32_t pops_used = 0;

	if (max_pops > 0) {
		for (size_t k = 0; k < cand_count; k++) {
			if (!cands[k].active) continue;
			uint64_t rkey = morton_quad_key(world, cands[k].level, cands[k].ox, cands[k].oy,
			                                sys->base_w, sys->base_h);
			uint32_t rslot = 0;
			bool cached = (islet_quad_map_get(&sys->map, rkey, &rslot) &&
			               rslot < sys->pool_capacity && sys->pool_meta[rslot].state == 1);
			if (!cached)
				pops_used++;
		}
	}

	while (cand_count > 0) {
		int best_idx = -1;
		float best_ratio = 1.0f;

		for (size_t k = 0; k < cand_count; k++) {
			if (!cands[k].active || cands[k].split_blocked || cands[k].level <= min_lvl)
				continue;

			float node_thresh = split_threshold;
			if (hyst > 0.0f) {
				uint64_t pkey = morton_quad_key(world, cands[k].level, cands[k].ox, cands[k].oy,
				                                sys->base_w, sys->base_h);
				bool kids_were_active = false;
				for (uint8_t q = 0; q < 4; q++) {
					uint64_t ckey = morton_quad_child(pkey, q);
					uint32_t cslot = 0;
					if (islet_quad_map_get(&sys->map, ckey, &cslot) &&
					    cslot < sys->pool_capacity && sys->pool_meta[cslot].state == 1 &&
					    (sys->pool_meta[cslot].active_view_mask & view_mask) != 0) {
						kids_were_active = true;
						break;
					}
				}
				/* If children were already active in this view, keep them split down to (1 - hyst);
				 * if this node was unsplit, require (1 + hyst) before splitting. */
				node_thresh *= kids_were_active ? (1.0f - hyst) : (1.0f + hyst);
			}

			float ratio = (node_thresh > 1e-5f) ? (cands[k].proj_w / node_thresh) : 0.0f;
			if (ratio > best_ratio ||
			    (ratio == best_ratio && best_idx >= 0 && cands[k].dist_sq < cands[best_idx].dist_sq)) {
				best_ratio = ratio;
				best_idx = (int)k;
			}
		}

		if (best_idx < 0)
			break;

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
			bool    cached;
		} kids[4];
		size_t vis_kids = 0;
		uint32_t uncached_vis_kids = 0;

		for (int c = 0; c < 4; c++) {
			kids[c].ox = pox + ((c & 1) ? cw : 0);
			kids[c].oy = poy + ((c & 2) ? ch : 0);
			kids[c].cached = false;
			float center_xyz[3], span_w, span_h, dist_sq, nz_min, nz_max;
			ar_compute_node_3d_box(sys, opts, world, child_lvl,
			                       kids[c].ox, kids[c].oy, cw, ch,
			                       eye_xyz, has_eye,
			                       kids[c].b_min, kids[c].b_max, center_xyz,
			                       &span_w, &span_h, &dist_sq,
			                       &nz_min, &nz_max);
			kids[c].visible = ar_frustum_test_box_sat(&frustum, has_faabb, faabb_min, faabb_max,
			                                          kids[c].b_min, kids[c].b_max) &&
			                  (!c3d || c3d(c3d_ctx, (float)kids[c].ox, (float)kids[c].oy,
			                               (float)cw, (float)ch,
			                               nz_min, nz_max,
			                               center_xyz, span_w, span_h)) &&
			                  !ar_occ_test_padded(occ, kids[c].b_min, kids[c].b_max,
			                                      occ_box_pad, occ_pad_all);
			if (kids[c].visible) {
				kids[c].proj_w = ar_projected_screen_width_box_3d(mvp, opts,
				                                                  center_xyz,
				                                                  eye_xyz, has_eye,
				                                                  span_w, span_h, vp_w, vp_h);
				kids[c].dist_sq = dist_sq;
				if (max_pops > 0) {
					uint64_t ckey = morton_quad_key(world, child_lvl, kids[c].ox, kids[c].oy,
					                                sys->base_w, sys->base_h);
					uint32_t cslot = 0;
					kids[c].cached = (islet_quad_map_get(&sys->map, ckey, &cslot) &&
					                  cslot < sys->pool_capacity && sys->pool_meta[cslot].state == 1);
					if (!kids[c].cached)
						uncached_vis_kids++;
				}
				vis_kids++;
			}
		}

		if (vis_kids == 0) {
			/* Swap-remove dead candidate immediately so cand_count budget is never starved */
			cands[best_idx] = cands[--cand_count];
			continue;
		}

		uint32_t parent_uncached = 0;
		if (max_pops > 0) {
			uint64_t pkey = morton_quad_key(world, cands[best_idx].level, pox, poy,
			                                sys->base_w, sys->base_h);
			uint32_t pslot = 0;
			bool p_cached = (islet_quad_map_get(&sys->map, pkey, &pslot) &&
			                 pslot < sys->pool_capacity && sys->pool_meta[pslot].state == 1);
			if (!p_cached && pops_used > 0)
				parent_uncached = 1;
		}
		uint32_t net_new_pops = (uncached_vis_kids > parent_uncached) ? (uncached_vis_kids - parent_uncached) : 0;

		if (max_pops > 0 && uncached_vis_kids > 0 && pops_used + net_new_pops > max_pops) {
			cands[best_idx].split_blocked = true;
			continue;
		}

		if (cand_count + vis_kids - 1 > max_candidates) {
			cands[best_idx].split_blocked = true;
			continue;
		}

		if (max_pops > 0) {
			pops_used += net_new_pops;
		}

		bool replaced = false;
		for (int c = 0; c < 4; c++) {
			if (!kids[c].visible) continue;
			size_t dst_idx = replaced ? cand_count++ : (size_t)best_idx;
			replaced = true;
			cands[dst_idx].ox = kids[c].ox;
			cands[dst_idx].oy = kids[c].oy;
			cands[dst_idx].level = child_lvl;
			cands[dst_idx].proj_w = kids[c].proj_w;
			cands[dst_idx].dist_sq = kids[c].dist_sq;
			memcpy(cands[dst_idx].b_min, kids[c].b_min, sizeof(kids[c].b_min));
			memcpy(cands[dst_idx].b_max, kids[c].b_max, sizeof(kids[c].b_max));
			cands[dst_idx].active = true;
			cands[dst_idx].split_blocked = false;
		}
	}

	/* 6. Populate Active Output Leaf Nodes (Front-to-Back for Conservative Occlusion) */
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

	for (size_t i = 0; i < sys->pool_capacity; i++) {
		if (sys->pool_meta[i].state == 1 && sys->pool[i].world == world) {
			sys->pool_meta[i].active_view_mask &= ~view_mask;
		}
	}

	ar_point_3d_fn pt3d = opts ? opts->point_3d_fn : NULL;
	void *pt3d_ctx = (opts && opts->bounds_ctx) ? opts->bounds_ctx : sys->bounds_ctx;

	for (size_t i = 0; i < active_count && out_count < max_candidates; i++) {
		if (occ && ar_occ_test_padded(occ, active_cands[i].b_min, active_cands[i].b_max,
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

			if (occ && !(opts && opts->occluder_no_write)) {
				float texel_span = (float)((node->w > node->h) ? node->w : node->h) /
				                   (float)((node->grid_dim > 0) ? node->grid_dim : 1u);
				float coarse_k = (node->grid_dim < 8u) ? (float)node->grid_dim : 8.0f;
				float cell_span = (float)((node->w > node->h) ? node->w : node->h) / coarse_k;
				float node_z_bias = opts->occluder_z_bias +
				                    texel_span * opts->occluder_texel_bias +
				                    (cell_span * cell_span * 0.125f) * opts->occluder_concave_inv_r;
				ar_occluder_rasterize_node(occ, node, sys->value_scale, node_z_bias,
				                           pt3d, pt3d_ctx);
			}
		}
	}

	sys->active_count = out_count;
	sys->stats.active_chunks_last_frame = out_count;

	size_t in_use = 0;
	for (size_t i = 0; i < sys->pool_capacity; i++) {
		if (sys->pool_meta[i].state == 1)
			in_use++;
	}
	sys->stats.pool_in_use = in_use;

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
	if (islet_quad_map_get(&sys->map, morton_key, &slot)) {
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
	for (size_t i = 0; i < sys->pool_capacity; i++) {
		if (sys->pool_meta[i].state == 1 && sys->pool[i].world == world) {
			sys->pool_meta[i].active_view_mask &= ~view_mask;
		}
	}

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

static bool ar_probe_internal(ar_system_t *sys, uint32_t view_mask, uint16_t world,
                              int32_t wx, int32_t wy, uint8_t max_ancestor_override,
                              ar_sample_t *out)
{
	if (!sys || !out)
		return false;

	uint64_t key0 = morton_quad_key(world, 0, wx, wy, sys->base_w, sys->base_h);
	uint8_t max_lvl = (sys->max_levels > 0) ? (uint8_t)(sys->max_levels - 1) : 16;
	const ar_node_t *best_node = NULL;

	/* 1. If a specific view_mask is requested, prefer the node active in that view */
	if (view_mask != 0) {
		uint64_t key = key0;
		for (uint8_t lvl = 0; lvl <= max_lvl; lvl++) {
			uint32_t slot = 0;
			if (islet_quad_map_get(&sys->map, key, &slot)) {
				if (slot < sys->pool_capacity && sys->pool_meta[slot].state == 1 &&
				    (sys->pool_meta[slot].active_view_mask & view_mask) != 0) {
					best_node = &sys->pool[slot];
					break;
				}
			}
			key = morton_quad_parent(key);
		}
	}

	/* 2. Fallback to finest cached ancestor (up to max_lvl for view probes, or configured ancestor cap for ground probes) */
	if (!best_node) {
		uint32_t slot = 0;
		uint8_t found_lvl = 0;
		uint8_t def_cap = (max_ancestor_override > 0) ? max_ancestor_override : sys->probe_ancestor_cap;
		if (def_cap == 0) def_cap = AR_DEFAULT_PROBE_ANCESTOR_CAP;
		uint8_t cap_lvl = (view_mask != 0) ? max_lvl : (max_lvl < def_cap ? max_lvl : def_cap);
		if (islet_quad_map_find_ancestor(&sys->map, key0, cap_lvl, &slot, &found_lvl)) {
			if (slot < sys->pool_capacity && sys->pool_meta[slot].state == 1) {
				best_node = &sys->pool[slot];
			}
		}
	}

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

		const ar_sample_t *s00 = &best_node->samples[row0 * vdim + col0];
		const ar_sample_t *s10 = &best_node->samples[row0 * vdim + col1];
		const ar_sample_t *s01 = &best_node->samples[row1 * vdim + col0];
		const ar_sample_t *s11 = &best_node->samples[row1 * vdim + col1];

		float h0 = (float)s00->val * (1.0f - wu) + (float)s10->val * wu;
		float h1 = (float)s01->val * (1.0f - wu) + (float)s11->val * wu;
		float h  = h0 * (1.0f - wv) + h1 * wv;
		out->val = (int16_t)(h >= 0.0f ? (h + 0.5f) : (h - 0.5f));

		const ar_sample_t *s_near = (u < 0.5f)
		                            ? ((v < 0.5f) ? s00 : s01)
		                            : ((v < 0.5f) ? s10 : s11);
		uint8_t n_attr = sys->attr_count;
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

		uint32_t c00 = s00->rgba, c10 = s10->rgba, c01 = s01->rgba, c11 = s11->rgba;
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
		out->flags = s_near->flags;
		return true;
	}

	sys->sampler(sys->user_ctx, world, 0, wx, wy, out);
	for (uint8_t k = sys->attr_count; k < AR_SAMPLE_ATTR_CAP; k++)
		out->attr[k] = 0;
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
