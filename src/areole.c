#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ttypt/morton.h"
#include "ttypt/areole.h"
#include "ttypt/chunk.h"
#include "ttypt/islet_quad.h"

/* Floor division for integer coordinate lattices */
static inline int32_t ar_fdiv(int32_t a, int32_t b)
{
	int32_t q = a / b;
	int32_t r = a % b;
	if (r != 0 && ((a < 0) ^ (b < 0)))
		q--;
	return q;
}

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

/* 4x4 Matrix Inversion */
bool ar_mat4_invert(float out[16], const float m[16])
{
	float inv[16];
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

	float det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];
	if (fabsf(det) < 1e-35f || !isfinite(det)) {
		ar_mat4_identity(out);
		return false;
	}

	float inv_det = 1.0f / det;
	for (int i = 0; i < 16; i++) {
		out[i] = inv[i] * inv_det;
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
/* Vertical Bounds Envelope & Screen Projection Helpers                      */
/* ========================================================================= */

static inline void ar_compute_bounds_envelope(const ar_cull_opts_t *opts,
                                              float ox, float oy, float w, float h,
                                              float *out_z_min, float *out_z_max)
{
	float z_min = opts ? opts->z_min : -2000.0f;
	float z_max = opts ? opts->z_max : 2000.0f;
	float margin = opts ? opts->margin : 0.0f;

	if (opts && opts->bounds_fn) {
		opts->bounds_fn(opts->bounds_ctx, ox, oy, w, h, &z_min, &z_max);
	} else if (opts && opts->sag_radius > 0.0f) {
		float sag_r = opts->sag_radius;
		float fx = opts->focal_x;
		float fy = opts->focal_y;

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
			delta_sag = sag_r + (r_max - sag_r) * 2.0f;
		}
		z_min -= delta_sag;
	}
	*out_z_min = z_min - margin;
	*out_z_max = z_max + margin;
}

static inline float ar_projected_screen_width(const float mvp[16],
                                              float ox, float oy, float w, float h,
                                              float z_mid,
                                              float vp_w, float vp_h)
{
	float cx = ox + w * 0.5f;
	float cy = oy + h * 0.5f;

	float cw = mvp[0 * 4 + 3] * cx +
	           mvp[1 * 4 + 3] * cy +
	           mvp[2 * 4 + 3] * z_mid +
	           mvp[3 * 4 + 3] * 1.0f;

	bool is_perspective = (fabsf(mvp[3]) > 1e-5f || fabsf(mvp[7]) > 1e-5f || fabsf(mvp[11]) > 1e-5f);

	float len_row0 = sqrtf(mvp[0]*mvp[0] + mvp[4]*mvp[4] + mvp[8]*mvp[8]);
	float len_row1 = sqrtf(mvp[1]*mvp[1] + mvp[5]*mvp[5] + mvp[9]*mvp[9]);
	if (len_row0 < 1e-7f) len_row0 = 0.001f;
	if (len_row1 < 1e-7f) len_row1 = 0.001f;

	float w_eff = 1.0f;
	if (is_perspective) {
		if (cw <= 1.0f) return (vp_w > vp_h ? vp_w : vp_h) * 2.0f;
		w_eff = cw;
	}

	float px_w = (w * len_row0 * (vp_w > 0.0f ? vp_w : 1280.0f) * 0.5f) / w_eff;
	float px_h = (h * len_row1 * (vp_h > 0.0f ? vp_h : 720.0f) * 0.5f) / w_eff;

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
	bool    active;
} ar_candidate_t;

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
	ar_sampler_fn     sampler;
	ar_sparsity_fn    sparsity;
	void             *user_ctx;
	uint32_t          current_tick;

	/* Pointers into single contiguous allocation */
	ar_node_t        *pool;
	const ar_node_t **active_chunks;
	ar_candidate_t   *cands_scratch;
	ar_candidate_t   *active_cands_scratch;
	size_t            active_count;

	/* Backing libislet in-memory quadtree spatial database */
	islet_quad_map_t  map;

	/* Telemetry & diagnostics */
	ar_stats_t        stats;
};

static inline size_t ar_align8(size_t n)
{
	return (n + 7u) & ~(size_t)7u;
}

ar_system_t *ar_system_create_ex(const ar_config_t *cfg)
{
	if (!cfg || cfg->base_w <= 0 || cfg->base_h <= 0 || !cfg->sampler)
		return NULL;

	uint16_t grid_dim = cfg->grid_dim ? cfg->grid_dim : AR_DEFAULT_GRID_DIM;
	uint16_t vert_dim = (uint16_t)(grid_dim + 1u);
	uint32_t vert_count = (uint32_t)vert_dim * (uint32_t)vert_dim;
	uint16_t pool_cap = cfg->pool_capacity ? cfg->pool_capacity : AR_DEFAULT_POOL_CAPACITY;
	uint16_t max_act = cfg->max_active ? cfg->max_active : AR_DEFAULT_MAX_ACTIVE;
	uint8_t max_lvls = cfg->max_levels ? cfg->max_levels : AR_DEFAULT_MAX_LEVELS;
	uint32_t def_rgba = cfg->default_rgba ? cfg->default_rgba : AR_DEFAULT_FALLBACK_RGBA;

	size_t sys_bytes = ar_align8(sizeof(ar_system_t));
	size_t pool_bytes = ar_align8((size_t)pool_cap * sizeof(ar_node_t));
	size_t active_bytes = ar_align8((size_t)max_act * sizeof(ar_node_t *));
	size_t cands_bytes = ar_align8((size_t)(max_act + 8u) * sizeof(ar_candidate_t));
	size_t active_cands_bytes = ar_align8((size_t)max_act * sizeof(ar_candidate_t));

	size_t per_node_vals = ar_align8((size_t)vert_count * sizeof(int16_t));
	size_t per_node_rgba = ar_align8((size_t)vert_count * sizeof(uint32_t));
	size_t per_node_samples = ar_align8((size_t)vert_count * sizeof(ar_sample_t));
	size_t per_node_arena = per_node_vals + per_node_rgba + per_node_samples;

	size_t total_bytes = sys_bytes + pool_bytes + active_bytes + cands_bytes +
	                     active_cands_bytes + (size_t)pool_cap * per_node_arena;

	uint8_t *arena = (uint8_t *)calloc(1, total_bytes);
	if (!arena)
		return NULL;

	ar_system_t *sys = (ar_system_t *)arena;
	uint8_t *cursor = arena + sys_bytes;

	sys->pool = (ar_node_t *)cursor;
	cursor += pool_bytes;

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
		sys->pool[i].values = (int16_t *)cursor;
		cursor += per_node_vals;
		sys->pool[i].rgba = (uint32_t *)cursor;
		cursor += per_node_rgba;
		sys->pool[i].samples = (ar_sample_t *)cursor;
		cursor += per_node_samples;
	}

	sys->base_w = cfg->base_w;
	sys->base_h = cfg->base_h;
	sys->grid_dim = grid_dim;
	sys->vert_dim = vert_dim;
	sys->vert_count = vert_count;
	sys->pool_capacity = pool_cap;
	sys->max_active = max_act;
	sys->max_levels = max_lvls;
	sys->default_rgba = def_rgba;
	sys->sampler = cfg->sampler;
	sys->sparsity = cfg->sparsity;
	sys->user_ctx = cfg->user_ctx;
	sys->current_tick = 1;

	if (!islet_quad_map_init(&sys->map, (size_t)pool_cap * 2u)) {
		free(arena);
		return NULL;
	}

	sys->stats.pool_capacity = pool_cap;
	sys->stats.min_level_observed = 255;
	sys->stats.max_level_observed = 0;

	return sys;
}

ar_system_t *ar_system_create(int32_t base_w, int32_t base_h,
                              ar_sampler_fn sampler,
                              ar_sparsity_fn sparsity,
                              void *user_ctx)
{
	ar_config_t cfg;
	memset(&cfg, 0, sizeof(cfg));
	cfg.base_w = base_w;
	cfg.base_h = base_h;
	cfg.sampler = sampler;
	cfg.sparsity = sparsity;
	cfg.user_ctx = user_ctx;
	return ar_system_create_ex(&cfg);
}

void ar_system_destroy(ar_system_t *sys)
{
	if (!sys)
		return;
	islet_quad_map_destroy(&sys->map);
	free(sys);
}

void ar_system_clear(ar_system_t *sys)
{
	if (!sys)
		return;
	for (size_t i = 0; i < sys->pool_capacity; i++) {
		sys->pool[i].state = 0;
		sys->pool[i].last_used_tick = 0;
		sys->pool[i].active_view_mask = 0;
	}
	islet_quad_map_clear(&sys->map);
	sys->active_count = 0;
}

void ar_system_begin_frame(ar_system_t *sys)
{
	if (!sys)
		return;
	sys->current_tick++;
	for (size_t i = 0; i < sys->pool_capacity; i++) {
		sys->pool[i].active_view_mask = 0;
	}
}

static ar_node_t *ar_get_or_populate(ar_system_t *sys, uint16_t world,
                                     uint8_t level, int32_t ox, int32_t oy)
{
	uint64_t key = morton_quad_key(world, level, ox, oy, sys->base_w, sys->base_h);

	/* 1. Fast O(1) Islet spatial map lookup */
	uint32_t hit_slot = 0;
	if (islet_quad_map_get(&sys->map, key, &hit_slot)) {
		if (hit_slot < sys->pool_capacity && sys->pool[hit_slot].state == 1) {
			ar_node_t *node = &sys->pool[hit_slot];
			node->last_used_tick = sys->current_tick;
			node->last_touch_frame = sys->current_tick;
			sys->stats.cache_hits_total++;
			return node;
		}
	}

	sys->stats.cache_misses_total++;

	/* 2. Find an empty slot or evict oldest LRU slot not active this frame */
	int slot = -1;
	for (size_t i = 0; i < sys->pool_capacity; i++) {
		if (sys->pool[i].state == 0) {
			slot = (int)i;
			break;
		}
	}

	if (slot < 0) {
		uint32_t oldest_tick = sys->current_tick;
		for (size_t i = 0; i < sys->pool_capacity; i++) {
			if (sys->pool[i].active_view_mask == 0) {
				if (sys->pool[i].last_used_tick < oldest_tick) {
					oldest_tick = sys->pool[i].last_used_tick;
					slot = (int)i;
				}
			}
		}
		if (slot < 0) {
			oldest_tick = sys->current_tick;
			for (size_t i = 0; i < sys->pool_capacity; i++) {
				if (sys->pool[i].last_used_tick < oldest_tick) {
					oldest_tick = sys->pool[i].last_used_tick;
					slot = (int)i;
				}
			}
		}
		if (slot < 0) {
			oldest_tick = 0xFFFFFFFFu;
			for (size_t i = 0; i < sys->pool_capacity; i++) {
				if (sys->pool[i].last_used_tick < oldest_tick) {
					oldest_tick = sys->pool[i].last_used_tick;
					slot = (int)i;
				}
			}
		}
		if (slot >= 0 && sys->pool[slot].state == 1) {
			islet_quad_map_del(&sys->map, sys->pool[slot].morton_key);
		}
	}

	if (slot < 0)
		return NULL;

	ar_node_t *node = &sys->pool[slot];
	node->world = world;
	node->level = level;
	node->ox = ox;
	node->oy = oy;
	node->w = sys->base_w << level;
	node->h = sys->base_h << level;
	node->last_used_tick = sys->current_tick;
	node->last_touch_frame = sys->current_tick;
	node->active_view_mask = 0;
	node->state = 1;
	node->flags = 0;
	node->min_val = 32767;
	node->max_val = -32768;
	node->morton_key = key;

	islet_quad_map_put(&sys->map, key, (uint32_t)slot);

	/* 3. Query sparsity if provider supplied */
	if (sys->sparsity) {
		if (sys->sparsity(sys->user_ctx, world, level, ox, oy, node->w, node->h)) {
			node->flags |= (1u << 31);
		}
	}

	/* 4. Evaluate sample grid: (grid_dim + 1) x (grid_dim + 1) vertices */
	int32_t gdim = (int32_t)sys->grid_dim;
	int32_t step_x = node->w / gdim;
	int32_t step_y = node->h / gdim;
	if (step_x < 1) step_x = 1;
	if (step_y < 1) step_y = 1;

	size_t idx = 0;
	for (int r = 0; r <= gdim; r++) {
		int32_t wy = oy + r * step_y;
		for (int c = 0; c <= gdim; c++) {
			int32_t wx = ox + c * step_x;
			ar_sample_t *s = &node->samples[idx];
			sys->sampler(sys->user_ctx, world, level, wx, wy, s);
			sys->stats.sampler_invocations_total++;

			if (s->val < node->min_val) node->min_val = s->val;
			if (s->val > node->max_val) node->max_val = s->val;
			node->flags |= s->flags;

			node->values[idx] = s->val;
			node->rgba[idx] = s->rgba ? s->rgba : sys->default_rgba;
			idx++;
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

	/* 2. Analytical Inversion for Viewport Ray Unprojection */
	float inv_mvp[16];
	bool has_inv = ar_mat4_invert(inv_mvp, mvp);

	/* 3. Determine root footprint bounds */
	float fx = opts ? opts->focal_x : 0.0f;
	float fy = opts ? opts->focal_y : 0.0f;
	float sag_r = opts ? opts->sag_radius : 0.0f;
	float margin = (opts && opts->margin > 0.0f) ? opts->margin : 150.0f;
	float vp_w = (opts && opts->viewport_w > 0.0f) ? opts->viewport_w : 1280.0f;
	float vp_h = (opts && opts->viewport_h > 0.0f) ? opts->viewport_h : 720.0f;

	float z_min = opts ? opts->z_min : -2000.0f;
	float z_max = opts ? opts->z_max : 2000.0f;
	if (z_min == 0.0f && z_max == 0.0f) {
		z_min = -2000.0f;
		z_max = 2000.0f;
	}
	float eff_z_min = z_min;
	if (sag_r > 0.0f) eff_z_min -= sag_r;
	float eff_z_max = z_max;

	bool is_perspective = (fabsf(mvp[3]) > 1e-5f || fabsf(mvp[7]) > 1e-5f || fabsf(mvp[11]) > 1e-5f);
	float focal_pt[4] = { fx, fy, 0.0f, 1.0f };
	float focal_clip[4];
	for (int i = 0; i < 4; i++) {
		focal_clip[i] = mvp[0 * 4 + i] * focal_pt[0] +
		               mvp[1 * 4 + i] * focal_pt[1] +
		               mvp[2 * 4 + i] * focal_pt[2] +
		               mvp[3 * 4 + i] * focal_pt[3];
	}
	float dist = (is_perspective && focal_clip[3] > 0.01f) ? fabsf(focal_clip[3]) : 1.0f;

	float len_row0 = sqrtf(mvp[0]*mvp[0] + mvp[4]*mvp[4] + mvp[8]*mvp[8]);
	float len_row1 = sqrtf(mvp[1]*mvp[1] + mvp[5]*mvp[5] + mvp[9]*mvp[9]);
	if (len_row0 < 1e-7f) len_row0 = 0.001f;
	if (len_row1 < 1e-7f) len_row1 = 0.001f;

	float bb_min_x = 1e9f, bb_max_x = -1e9f;
	float bb_min_y = 1e9f, bb_max_y = -1e9f;
	size_t unproj_pts = 0;

	if (has_inv) {
		if (!opts) {
			float c_near[3], c_far[3];
			if (ar_unproject_ndc(inv_mvp, 0.0f, 0.0f, -1.0f, c_near) &&
			    ar_unproject_ndc(inv_mvp, 0.0f, 0.0f,  1.0f, c_far)) {
				float c_dz = c_far[2] - c_near[2];
				if (fabsf(c_dz) > 1e-6f) {
					float tc = -c_near[2] / c_dz;
					if (tc >= -0.5f && tc <= 1.5f) {
						fx = c_near[0] + tc * (c_far[0] - c_near[0]);
						fy = c_near[1] + tc * (c_far[1] - c_near[1]);
					}
				}
			}
		}

		const float ndc_corners[4][2] = {
			{ -1.0f, -1.0f },
			{  1.0f, -1.0f },
			{  1.0f,  1.0f },
			{ -1.0f,  1.0f }
		};

		for (int k = 0; k < 4; k++) {
			float p_near[3], p_far[3];
			if (!ar_unproject_ndc(inv_mvp, ndc_corners[k][0], ndc_corners[k][1], -1.0f, p_near) ||
			    !ar_unproject_ndc(inv_mvp, ndc_corners[k][0], ndc_corners[k][1],  1.0f, p_far)) {
				continue;
			}

			float dx = p_far[0] - p_near[0];
			float dy = p_far[1] - p_near[1];
			float dz = p_far[2] - p_near[2];

			if (fabsf(dz) > 1e-6f) {
				float t0 = -p_near[2] / dz;
				if (t0 >= -0.2f && t0 <= 1.2f) {
					float ix = p_near[0] + t0 * dx;
					float iy = p_near[1] + t0 * dy;
					if (ix >= -100000000.0f && ix <= 100000000.0f &&
					    iy >= -100000000.0f && iy <= 100000000.0f) {
						if (ix < bb_min_x) bb_min_x = ix;
						if (ix > bb_max_x) bb_max_x = ix;
						if (iy < bb_min_y) bb_min_y = iy;
						if (iy > bb_max_y) bb_max_y = iy;
						unproj_pts++;
					}
				}

				if (!is_perspective || sag_r > 0.0f || (opts && opts->bounds_fn)) {
					const float test_z[2] = { eff_z_min, eff_z_max };
					for (int p = 0; p < 2; p++) {
						float t = (test_z[p] - p_near[2]) / dz;
						if (t >= -0.2f && t <= 1.2f) {
							float ix = p_near[0] + t * dx;
							float iy = p_near[1] + t * dy;
							if (ix >= -100000000.0f && ix <= 100000000.0f &&
							    iy >= -100000000.0f && iy <= 100000000.0f) {
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
		}
	}

	if (unproj_pts >= 2 && bb_max_x > bb_min_x && bb_max_y > bb_min_y) {
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
	} else {
		float ground_w = (2.0f * dist) / len_row0;
		float ground_h = (2.0f * dist * (is_perspective ? 1.5f : 1.0f)) / len_row1;

		if (sag_r > 0.0f && is_perspective) {
			if (ground_w < sag_r * 2.0f) ground_w = sag_r * 2.0f;
			if (ground_h < sag_r * 2.0f) ground_h = sag_r * 2.0f;
		}

		float half_gw = ground_w * 0.5f + margin;
		float half_gh = ground_h * 0.5f + margin;

		bb_min_x = fx - half_gw;
		bb_max_x = fx + half_gw;
		bb_min_y = fy - half_gh;
		bb_max_y = fy + half_gh;
	}

	/* 4. Determine root level that naturally covers the bounding box (cols * rows <= 16) */
	uint8_t root_level = 0;
	while (root_level + 1 < sys->max_levels) {
		int32_t rw = sys->base_w << root_level;
		int32_t rh = sys->base_h << root_level;
		int32_t gi0 = ar_fdiv((int32_t)bb_min_x, rw);
		int32_t gi1 = ar_fdiv((int32_t)bb_max_x, rw);
		int32_t gj0 = ar_fdiv((int32_t)bb_min_y, rh);
		int32_t gj1 = ar_fdiv((int32_t)bb_max_y, rh);
		int64_t cols = (int64_t)gi1 - gi0 + 1;
		int64_t rows = (int64_t)gj1 - gj0 + 1;
		if (cols > 0 && rows > 0 && cols * rows <= 16)
			break;
		root_level++;
	}

	if (opts && opts->min_level > 0 && root_level < opts->min_level)
		root_level = opts->min_level;
	if (opts && opts->max_level > 0 && root_level > opts->max_level)
		root_level = opts->max_level;

	int32_t root_w = sys->base_w << root_level;
	int32_t root_h = sys->base_h << root_level;

	int32_t i0 = ar_fdiv((int32_t)bb_min_x, root_w);
	int32_t i1 = ar_fdiv((int32_t)bb_max_x, root_w);
	int32_t j0 = ar_fdiv((int32_t)bb_min_y, root_h);
	int32_t j1 = ar_fdiv((int32_t)bb_max_y, root_h);

	while (((int64_t)i1 - i0 + 1) * ((int64_t)j1 - j0 + 1) > 16 && root_level + 1 < sys->max_levels) {
		root_level++;
		root_w = sys->base_w << root_level;
		root_h = sys->base_h << root_level;
		i0 = ar_fdiv((int32_t)bb_min_x, root_w);
		i1 = ar_fdiv((int32_t)bb_max_x, root_w);
		j0 = ar_fdiv((int32_t)bb_min_y, root_h);
		j1 = ar_fdiv((int32_t)bb_max_y, root_h);
	}

	ar_candidate_t *cands = sys->cands_scratch;
	size_t cand_count = 0;

	for (int32_t j = j0; j <= j1; j++) {
		for (int32_t i = i0; i <= i1; i++) {
			int32_t rox = i * root_w;
			int32_t roy = j * root_h;
			float node_z_min, node_z_max;
			ar_compute_bounds_envelope(opts, (float)rox, (float)roy, (float)root_w, (float)root_h,
			                           &node_z_min, &node_z_max);

			if (ar_frustum_test_box(&frustum, (float)rox, (float)roy, node_z_min,
			                        (float)(rox + root_w), (float)(roy + root_h), node_z_max)) {
				if (cand_count < max_candidates) {
					float z_mid = (node_z_min + node_z_max) * 0.5f;
					float pw = ar_projected_screen_width(mvp, (float)rox, (float)roy,
					                                     (float)root_w, (float)root_h,
					                                     z_mid, vp_w, vp_h);
					cands[cand_count].ox = rox;
					cands[cand_count].oy = roy;
					cands[cand_count].level = root_level;
					cands[cand_count].proj_w = pw;
					cands[cand_count].dist_sq = ar_dist_to_box_sq(fx, fy, (float)rox, (float)roy,
					                                              (float)root_w, (float)root_h);
					cands[cand_count].active = true;
					cand_count++;
				}
			}
		}
	}

	/* 5. Priority Refinement Loop */
	float target_quad_px = (opts && opts->target_quad_px > 0.0f) ? opts->target_quad_px : 2.0f;
	float split_threshold = target_quad_px * (float)sys->grid_dim;
	uint8_t min_lvl = opts ? opts->min_level : 0;

	while (cand_count + 3 <= max_candidates) {
		int best_idx = -1;
		float max_err = split_threshold;

		for (size_t k = 0; k < cand_count; k++) {
			if (!cands[k].active || cands[k].level <= min_lvl)
				continue;
			if (cands[k].proj_w > max_err ||
			    (cands[k].proj_w == max_err && best_idx >= 0 && cands[k].dist_sq < cands[best_idx].dist_sq)) {
				max_err = cands[k].proj_w;
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
			bool    visible;
		} kids[4];
		size_t vis_kids = 0;

		for (int c = 0; c < 4; c++) {
			kids[c].ox = pox + ((c & 1) ? cw : 0);
			kids[c].oy = poy + ((c & 2) ? ch : 0);
			float kid_z_min, kid_z_max;
			ar_compute_bounds_envelope(opts, (float)kids[c].ox, (float)kids[c].oy,
			                           (float)cw, (float)ch, &kid_z_min, &kid_z_max);
			kids[c].visible = ar_frustum_test_box(&frustum,
			                                      (float)kids[c].ox, (float)kids[c].oy, kid_z_min,
			                                      (float)(kids[c].ox + cw), (float)(kids[c].oy + ch), kid_z_max);
			if (kids[c].visible) {
				float z_mid = (kid_z_min + kid_z_max) * 0.5f;
				kids[c].proj_w = ar_projected_screen_width(mvp, (float)kids[c].ox, (float)kids[c].oy,
				                                           (float)cw, (float)ch, z_mid, vp_w, vp_h);
				kids[c].dist_sq = ar_dist_to_box_sq(fx, fy, (float)kids[c].ox, (float)kids[c].oy,
				                                     (float)cw, (float)ch);
				vis_kids++;
			}
		}

		if (vis_kids == 0) {
			cands[best_idx].active = false;
			continue;
		}

		if (cand_count + vis_kids - 1 > max_candidates)
			break;

		bool replaced = false;
		for (int c = 0; c < 4; c++) {
			if (!kids[c].visible) continue;
			if (!replaced) {
				cands[best_idx].ox = kids[c].ox;
				cands[best_idx].oy = kids[c].oy;
				cands[best_idx].level = child_lvl;
				cands[best_idx].proj_w = kids[c].proj_w;
				cands[best_idx].dist_sq = kids[c].dist_sq;
				cands[best_idx].active = true;
				replaced = true;
			} else {
				cands[cand_count].ox = kids[c].ox;
				cands[cand_count].oy = kids[c].oy;
				cands[cand_count].level = child_lvl;
				cands[cand_count].proj_w = kids[c].proj_w;
				cands[cand_count].dist_sq = kids[c].dist_sq;
				cands[cand_count].active = true;
				cand_count++;
			}
		}
	}

	/* 6. Populate Active Output Leaf Nodes */
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

	for (size_t i = 0; i < active_count && out_count < max_candidates; i++) {
		ar_node_t *node = ar_get_or_populate(sys, world,
		                                     active_cands[i].level,
		                                     active_cands[i].ox,
		                                     active_cands[i].oy);
		if (node) {
			out_nodes[out_count++] = node;
			if (node->level > sys->stats.max_level_observed)
				sys->stats.max_level_observed = node->level;
			if (node->level < sys->stats.min_level_observed)
				sys->stats.min_level_observed = node->level;
		}
	}

	sys->active_count = out_count;
	sys->stats.active_chunks_last_frame = out_count;

	size_t in_use = 0;
	for (size_t i = 0; i < sys->pool_capacity; i++) {
		if (sys->pool[i].state == 1)
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

size_t ar_system_update(ar_system_t *sys, uint16_t world,
                        const ar_view_t *view,
                        const ar_node_t **out_nodes,
                        size_t max_out)
{
	if (!sys || !view || !out_nodes || max_out == 0)
		return 0;

	float zoom = (view->zoom > 0.0f) ? view->zoom : 1.0f;
	float eff_cam_dist = (view->eff_cam_dist > 1.0f) ? view->eff_cam_dist : 280.0f;
	float eye_dist = eff_cam_dist / zoom;

	float eye[3], target[3], up[3];
	target[0] = view->cam_x;
	target[1] = view->cam_y;
	target[2] = view->cam_z;

	if (view->fov <= 0.001f && view->pitch <= 0.001f) {
		eye[0] = view->cam_x;
		eye[1] = view->cam_y;
		eye[2] = view->cam_z + (eye_dist > 100.0f ? eye_dist : 500.0f);
		up[0] = 0.0f;
		up[1] = -1.0f;
		up[2] = 0.0f;
	} else {
		float pitch = (view->pitch > 0.001f) ? view->pitch : (55.0f * 3.14159265f / 180.0f);
		float cos_p = cosf(pitch);
		float sin_p = sinf(pitch);
		eye[0] = view->cam_x;
		eye[1] = view->cam_y + eye_dist * cos_p;
		eye[2] = view->cam_z + eye_dist * sin_p;
		up[0] = 0.0f;
		up[1] = -sin_p;
		up[2] = cos_p;
	}

	float modelview[16], proj[16], mvp[16];
	ar_mat4_lookat(modelview, eye, target, up);

	float vp_w = (view->viewport_w > 0.0f) ? view->viewport_w : 1280.0f;
	float vp_h = (view->viewport_h > 0.0f) ? view->viewport_h : 720.0f;

	if (view->fov <= 0.001f) {
		float hw = (vp_w * 0.5f) / zoom;
		float hh = (vp_h * 0.5f) / zoom;
		float near_z = (eye_dist * 2.0f > 2000.0f) ? (-eye_dist * 2.0f) : -2000.0f;
		float far_z = (eye_dist * 2.0f > 2000.0f) ? (eye_dist * 2.0f) : 2000.0f;
		ar_mat4_ortho(proj, -hw, hw, -hh, hh, near_z, far_z);
	} else {
		float aspect = vp_w / vp_h;
		float near_z = (eye_dist * 0.001f > 1.0f) ? (eye_dist * 0.001f) : 1.0f;
		float far_z = (eye_dist * 10.0f > 50000.0f) ? (eye_dist * 10.0f) : 50000.0f;
		ar_mat4_perspective(proj, view->fov, aspect, near_z, far_z);
	}
	ar_mat4_multiply(mvp, proj, modelview);

	ar_cull_opts_t opts;
	memset(&opts, 0, sizeof(opts));
	opts.z_min = -2000.0f;
	opts.z_max = 2000.0f;
	opts.sag_radius = 0.0f;
	opts.focal_x = view->cam_x;
	opts.focal_y = view->cam_y;
	opts.margin = (view->fov <= 0.001f) ? 0.0f : 150.0f;
	opts.target_quad_px = 2.0f;
	opts.viewport_w = vp_w;
	opts.viewport_h = vp_h;
	opts.min_level = view->min_level;
	opts.max_level = view->max_level;

	return ar_system_query_matrix(sys, world, mvp, &opts, out_nodes, max_out);
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
		if (slot < sys->pool_capacity && sys->pool[slot].state == 1)
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

size_t ar_system_eval_view(ar_system_t *sys, uint8_t view_id, uint16_t world,
                           const ar_view_t *view,
                           ar_draw_item_t *out_items,
                           size_t max_items)
{
	if (!sys || !view || !out_items || max_items == 0)
		return 0;

	size_t cap = max_items < sys->max_active ? max_items : sys->max_active;
	size_t count = ar_system_update(sys, world, view, sys->active_chunks, cap);

	for (size_t i = 0; i < count; i++) {
		const ar_node_t *node = sys->active_chunks[i];
		if (node) {
			((ar_node_t *)node)->active_view_mask |= (1 << (view_id & 7));
			((ar_node_t *)node)->last_touch_frame = sys->current_tick;
		}

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

bool ar_system_probe(ar_system_t *sys, uint16_t world,
                     int32_t wx, int32_t wy,
                     ar_sample_t *out)
{
	if (!sys || !out)
		return false;

	uint64_t key0 = morton_quad_key(world, 0, wx, wy, sys->base_w, sys->base_h);
	uint32_t slot = 0;
	uint8_t found_lvl = 0;
	const ar_node_t *best_node = NULL;

	if (islet_quad_map_find_ancestor(&sys->map, key0, 4, &slot, &found_lvl)) {
		if (slot < sys->pool_capacity && sys->pool[slot].state == 1) {
			best_node = &sys->pool[slot];
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
		float wu = u * u * (3.0f - 2.0f * u);
		float wv = v * v * (3.0f - 2.0f * v);

		const ar_sample_t *s00 = &best_node->samples[row0 * vdim + col0];
		const ar_sample_t *s10 = &best_node->samples[row0 * vdim + col1];
		const ar_sample_t *s01 = &best_node->samples[row1 * vdim + col0];
		const ar_sample_t *s11 = &best_node->samples[row1 * vdim + col1];

		float h0 = (float)s00->val * (1.0f - wu) + (float)s10->val * wu;
		float h1 = (float)s01->val * (1.0f - wu) + (float)s11->val * wu;
		float h  = h0 * (1.0f - wv) + h1 * wv;
		out->val = (int16_t)(h >= 0.0f ? (h + 0.5f) : (h - 0.5f));

		for (int k = 0; k < 8; k++) {
			float a0 = (float)s00->attr[k] * (1.0f - wu) + (float)s10->attr[k] * wu;
			float a1 = (float)s01->attr[k] * (1.0f - wu) + (float)s11->attr[k] * wu;
			float a  = a0 * (1.0f - wv) + a1 * wv;
			int ival = (int)(a + 0.5f);
			if (ival < 0) ival = 0;
			if (ival > 255) ival = 255;
			out->attr[k] = (uint8_t)ival;
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
		out->flags = s00->flags | s10->flags | s01->flags | s11->flags;
		return true;
	}

	sys->sampler(sys->user_ctx, world, 0, wx, wy, out);
	return true;
}

void ar_system_get_stats(const ar_system_t *sys, ar_stats_t *out_stats)
{
	if (!sys || !out_stats)
		return;
	*out_stats = sys->stats;
}

/* ========================================================================= */
/* Legacy chunk_* ABI Symbol Exports (for shared library compatibility)      */
/* ========================================================================= */

#undef chunk_system_create_ex
#undef chunk_system_create
#undef chunk_system_destroy
#undef chunk_system_clear
#undef chunk_system_begin_frame
#undef chunk_system_eval_view
#undef chunk_mat4_identity
#undef chunk_mat4_multiply
#undef chunk_mat4_ortho
#undef chunk_mat4_perspective
#undef chunk_mat4_lookat
#undef chunk_mat4_lookat_lh
#undef chunk_mat4_invert
#undef chunk_frustum_from_matrix
#undef chunk_frustum_test_box
#undef chunk_system_query_matrix
#undef chunk_system_query_modelview_projection
#undef chunk_system_update
#undef chunk_system_probe
#undef chunk_system_find
#undef chunk_system_find_by_key
#undef chunk_system_ensure
#undef chunk_system_get_stats

ar_system_t *chunk_system_create(int32_t base_w, int32_t base_h,
                                 ar_sampler_fn sampler,
                                 ar_sparsity_fn sparsity,
                                 void *user_ctx)
{
	return ar_system_create(base_w, base_h, sampler, sparsity, user_ctx);
}

void chunk_system_destroy(ar_system_t *cs) { ar_system_destroy(cs); }
void chunk_system_clear(ar_system_t *cs) { ar_system_clear(cs); }
void chunk_system_begin_frame(ar_system_t *cs) { ar_system_begin_frame(cs); }

size_t chunk_system_eval_view(ar_system_t *cs, uint8_t view_id, uint16_t world,
                              const ar_view_t *view,
                              ar_draw_item_t *out_items,
                              size_t max_items)
{
	return ar_system_eval_view(cs, view_id, world, view, out_items, max_items);
}

void chunk_mat4_identity(float out[16]) { ar_mat4_identity(out); }
void chunk_mat4_multiply(float out[16], const float a[16], const float b[16]) { ar_mat4_multiply(out, a, b); }
void chunk_mat4_ortho(float out[16], float l, float r, float b, float t, float n, float f) { ar_mat4_ortho(out, l, r, b, t, n, f); }
void chunk_mat4_perspective(float out[16], float fov, float asp, float n, float f) { ar_mat4_perspective(out, fov, asp, n, f); }
void chunk_mat4_lookat(float out[16], const float eye[3], const float target[3], const float up[3]) { ar_mat4_lookat(out, eye, target, up); }
bool chunk_mat4_invert(float out[16], const float m[16]) { return ar_mat4_invert(out, m); }
void chunk_frustum_from_matrix(ar_frustum_t *f, const float mvp[16]) { ar_frustum_from_matrix(f, mvp); }
bool chunk_frustum_test_box(const ar_frustum_t *f, float x0, float y0, float z0, float x1, float y1, float z1) { return ar_frustum_test_box(f, x0, y0, z0, x1, y1, z1); }

size_t chunk_system_query_matrix(ar_system_t *cs, uint16_t world,
                                 const float mvp[16],
                                 const ar_cull_opts_t *opts,
                                 const ar_node_t **out_chunks,
                                 size_t max_out)
{
	return ar_system_query_matrix(cs, world, mvp, opts, out_chunks, max_out);
}

size_t chunk_system_query_modelview_projection(ar_system_t *cs, uint16_t world,
                                               const float modelview[16],
                                               const float projection[16],
                                               const ar_cull_opts_t *opts,
                                               const ar_node_t **out_chunks,
                                               size_t max_out)
{
	return ar_system_query_modelview_projection(cs, world, modelview, projection, opts, out_chunks, max_out);
}

size_t chunk_system_update(ar_system_t *cs, uint16_t world,
                           const ar_view_t *view,
                           const ar_node_t **out_chunks,
                           size_t max_out)
{
	return ar_system_update(cs, world, view, out_chunks, max_out);
}

bool chunk_system_probe(ar_system_t *cs, uint16_t world,
                        int32_t wx, int32_t wy,
                        ar_sample_t *out)
{
	return ar_system_probe(cs, world, wx, wy, out);
}

const ar_node_t *chunk_system_find(ar_system_t *cs, uint16_t world,
                                   uint8_t level, int32_t ox, int32_t oy)
{
	return ar_system_find(cs, world, level, ox, oy);
}

const ar_node_t *chunk_system_find_by_key(ar_system_t *cs, uint64_t morton_key)
{
	return ar_system_find_by_key(cs, morton_key);
}

const ar_node_t *chunk_system_ensure(ar_system_t *cs, uint16_t world,
                                     uint8_t level, int32_t ox, int32_t oy)
{
	return ar_system_ensure(cs, world, level, ox, oy);
}

void chunk_system_get_stats(const ar_system_t *cs, ar_stats_t *out_stats)
{
	ar_system_get_stats(cs, out_stats);
}
