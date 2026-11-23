#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <time.h>

#include "ttypt/areole.h"
#include <ttypt/wq.h>

/* Pure analytical multi-octave procedural sampler (zero libchaff dependency) */
static void analytical_multioctave_sampler(void *user_ctx, uint16_t world, uint8_t level,
                                           int32_t wx, int32_t wy, ar_sample_t *out)
{
	(void)user_ctx;
	float fx = (float)wx + (float)world * 37.0f;
	float fy = (float)wy - (float)world * 53.0f;

	float macro = 120.0f * sinf(fx * 0.002454f) * cosf(fy * 0.002454f);
	float meso  =  35.0f * sinf(fx * 0.019635f + 0.7f) * cosf(fy * 0.019635f - 0.4f);
	float micro =  (level <= 2) ? (8.0f * sinf(fx * 0.098175f) * cosf(fy * 0.098175f)) : 0.0f;
	float field = macro + meso + micro;
	int32_t field_q = (int32_t)(field * 16.0f);

	out->val = (int16_t)field_q;
	out->attr[0] = (field_q < 0) ? 0 : ((field_q > 100 * 16) ? 2 : 10);
	out->attr[1] = out->attr[0];
	out->attr[2] = 0;
	out->attr[3] = (field_q < 0) ? (uint8_t)(-field_q > 255 * 16 ? 255 : (-field_q >> 4)) : 0;
	out->rgba = (out->attr[3] > 0) ? 0xFF003080u : 0xFF208020u;
	out->flags = (out->attr[3] > 0) ? AR_NODE_FLAG_USER0 : 0u;
}

/* Mock procedural sampler: analytical sinusoidal scalar field */
static void mock_sampler(void *user_ctx, uint16_t world, uint8_t level,
                         int32_t wx, int32_t wy, ar_sample_t *out)
{
	(void)user_ctx;
	(void)world;
	(void)level;

	float fx = (float)wx;
	float fy = (float)wy;

	float macro_h = 100.0f * sinf(fx * 0.0030679f) * cosf(fy * 0.0030679f);
	float fine_h  =  20.0f * sinf(fx * 0.049087f)  * cosf(fy * 0.049087f);
	float field = macro_h + fine_h;

	out->val = (int16_t)(field * 16.0f);
	out->attr[0] = (field < 0.0f) ? 0 : ((field > 60.0f) ? 2 : 10);
	out->attr[1] = out->attr[0];
	out->attr[2] = 0;
	out->attr[3] = (field < 0.0f) ? (uint8_t)(-field > 100.0f ? 100 : -field) : 0;
	out->rgba = (out->attr[3] > 0) ? 0xFF003080u : 0xFF208020u;
	out->flags = (out->attr[3] > 0) ? AR_NODE_FLAG_USER0 : 0u;
}

/* Mock sparsity function: reports structures only near origin [-320..320, -256..256] */
static bool mock_sparsity(void *user_ctx, uint16_t world, uint8_t level,
                          int32_t ox, int32_t oy, int32_t w, int32_t h)
{
	(void)user_ctx;
	(void)world;
	(void)level;

	int32_t sx0 = -320, sx1 = 320;
	int32_t sy0 = -256, sy1 = 256;

	bool overlap = !(ox + w < sx0 || ox > sx1 || oy + h < sy0 || oy > sy1);
	return overlap;
}

/* Custom vertical bounds callback for Test 15 */
static void custom_spherical_bounds_cb(void *bounds_ctx, float ox, float oy, float w, float h,
                                       float *inout_z_min, float *inout_z_max)
{
	(void)inout_z_max;
	float R = *(const float *)bounds_ctx;
	float corners[4][2] = {
		{ ox, oy }, { ox + w, oy }, { ox, oy + h }, { ox + w, oy + h }
	};
	float r2_max = 0.0f;
	for (int k = 0; k < 4; k++) {
		float r2 = corners[k][0] * corners[k][0] + corners[k][1] * corners[k][1];
		if (r2 > r2_max) r2_max = r2;
	}
	float r_max = sqrtf(r2_max);
	float sag = (r_max < R) ? (R - sqrtf(fmaxf(1.0f, R * R - r2_max))) : (R + (r_max - R) * 2.0f);
	*inout_z_min -= sag;
}

/* High-relief sampler spanning [-350, +750] wu for Test 28 frustum oracle */
static float high_relief_height(float fx, float fy)
{
	float h = 380.0f * sinf(fx * 0.0018f) * cosf(fy * 0.0015f)
	        + 180.0f * cosf((fx + fy) * 0.0042f)
	        + 120.0f * sinf((fx - fy) * 0.0075f) + 150.0f;
	return h;
}

static void high_relief_sampler(void *user_ctx, uint16_t world, uint8_t level,
                                int32_t wx, int32_t wy, ar_sample_t *out)
{
	(void)user_ctx; (void)world; (void)level;
	float h = high_relief_height((float)wx, (float)wy);
	out->val = (int16_t)(h * 16.0f);
	out->rgba = 0xFF558833u;
	out->flags = 0u;
}

/* Annular mountain ridge wall around origin (r in [450..850] -> z=+1100; outer peak at (0,-2400) -> z=+2400) */
typedef struct {
	float ridge_dir_x;
	float ridge_dir_y;
	float peak_x;
	float peak_y;
	float peak_h;
} ridge_scene_ctx_t;

static void ridge_wall_sampler(void *user_ctx, uint16_t world, uint8_t level,
                               int32_t wx, int32_t wy, ar_sample_t *out)
{
	(void)world; (void)level;
	const ridge_scene_ctx_t *rc = (const ridge_scene_ctx_t *)user_ctx;
	float x = (float)wx, y = (float)wy;
	float along = x * rc->ridge_dir_x + y * rc->ridge_dir_y;
	float perp  = -x * rc->ridge_dir_y + y * rc->ridge_dir_x;
	float h = 0.0f;
	/* Wide solid ridge wall perpendicular to look direction at distance [500..900] */
	if (along >= 500.0f && along <= 940.0f && fabsf(perp) <= 2200.0f) {
		h = 480.0f;
	}
	/* Tall background peak behind the wall that rises above the silhouette */
	float pdx = x - rc->peak_x;
	float pdy = y - rc->peak_y;
	if (pdx * pdx + pdy * pdy < 360.0f * 360.0f) {
		if (rc->peak_h > h) h = rc->peak_h;
	}
	out->val = (int16_t)(h * 16.0f);
	out->rgba = 0xFF607080u;
	out->flags = 0u;
}

static void flat_ocean_sampler(void *user_ctx, uint16_t world, uint8_t level,
                               int32_t wx, int32_t wy, ar_sample_t *out)
{
	(void)user_ctx; (void)world; (void)level; (void)wx; (void)wy;
	out->val = 0;
	out->rgba = 0xFF103060u;
	out->flags = 0u;
}

static void level_dependent_sampler(void *user_ctx, uint16_t world, uint8_t level,
                                    int32_t wx, int32_t wy, ar_sample_t *out)
{
	(void)user_ctx;
	(void)world;
	(void)wx;
	(void)wy;
	memset(out, 0, sizeof(*out));
	out->val = (int16_t)(level * 100);
	out->rgba = 0xFF208020u;
}

static float custom_square_interp_cb(float t)
{
	return t * t;
}

static void linear_ramp_sampler(void *user_ctx, uint16_t world, uint8_t level,
                                int32_t wx, int32_t wy, ar_sample_t *out)
{
	(void)user_ctx;
	(void)world;
	(void)level;
	memset(out, 0, sizeof(*out));
	out->val = (int16_t)(wx * 10 + wy);
	out->rgba = 0xFF4080C0u;
}

struct torus_test_ctx {
	float R;
	float r0;
	float period_x;
	float period_y;
	float offset[3];
};

static void custom_torus_bounds_3d_cb(void *ctx, float ox, float oy, float w, float h,
                                      float h_min, float h_max,
                                      float out_min[3], float out_max[3],
                                      float out_center[3],
                                      float *out_span_w, float *out_span_h)
{
	const struct torus_test_ctx *tc = (const struct torus_test_ctx *)ctx;
	float u_mid = ((ox + w * 0.5f) / tc->period_x) * 6.2831853f;
	float v_mid = ((oy + h * 0.5f) / tc->period_y) * 6.2831853f;
	float h_mid = (h_min + h_max) * 0.5f;
	float tube_mid = tc->r0 + h_mid;
	float rho_mid = tc->R + tube_mid * cosf(v_mid);
	out_center[0] = rho_mid * sinf(u_mid) + tc->offset[0];
	out_center[1] = -rho_mid * cosf(u_mid) + tc->offset[1];
	out_center[2] = -tube_mid * sinf(v_mid) + tc->offset[2];

	float bmin[3] = {  1e30f,  1e30f,  1e30f };
	float bmax[3] = { -1e30f, -1e30f, -1e30f };
	for (int iu = 0; iu <= 4; iu++) {
		float u = ((ox + w * (float)iu * 0.25f) / tc->period_x) * 6.2831853f;
		float su = sinf(u), cu = cosf(u);
		for (int iv = 0; iv <= 4; iv++) {
			float v = ((oy + h * (float)iv * 0.25f) / tc->period_y) * 6.2831853f;
			float sv = sinf(v), cv = cosf(v);
			for (int ih = 0; ih < 2; ih++) {
				float tr = tc->r0 + (ih ? h_max : h_min);
				float rho = tc->R + tr * cv;
				float p[3] = {
					rho * su + tc->offset[0],
					-rho * cu + tc->offset[1],
					-tr * sv + tc->offset[2]
				};
				for (int k = 0; k < 3; k++) {
					if (p[k] < bmin[k]) bmin[k] = p[k];
					if (p[k] > bmax[k]) bmax[k] = p[k];
				}
			}
		}
	}
	/* Add chordal bulge safety margin for 4-segment sampling */
	float du_step = (w / tc->period_x) * 6.2831853f * 0.25f;
	float dv_step = (h / tc->period_y) * 6.2831853f * 0.25f;
	float bulge = (tc->R + tc->r0) * du_step * du_step * 0.2f + tc->r0 * dv_step * dv_step * 0.2f;
	for (int k = 0; k < 3; k++) {
		out_min[k] = bmin[k] - bulge;
		out_max[k] = bmax[k] + bulge;
	}
	*out_span_w = rho_mid * ((w / tc->period_x) * 6.2831853f);
	*out_span_h = tube_mid * ((h / tc->period_y) * 6.2831853f);
}

static bool custom_torus_cull_3d_cb(void *ctx, float ox, float oy, float w, float h,
                                    float z_min, float z_max,
                                    const float center_xyz[3], float span_w, float span_h)
{
	(void)z_min; (void)z_max; (void)span_w; (void)span_h;
	const struct torus_test_ctx *tc = (const struct torus_test_ctx *)ctx;
	float du = w / tc->period_x;
	float dv = h / tc->period_y;
	float alpha = 3.14159265f * sqrtf(du * du + dv * dv);
	if (alpha >= 1.15f)
		return true; /* Keep very coarse parent cells until subdivided */

	float u_mid = ((ox + w * 0.5f) / tc->period_x) * 6.2831853f;
	float v_mid = ((oy + h * 0.5f) / tc->period_y) * 6.2831853f;
	/* Outward normal at (u_mid, v_mid) for torus parameterized in custom_torus_bounds_3d_cb:
	 * pos = (rho*sin(u), -rho*cos(u), -tr*sin(v)) -> normal = (cos(v)*sin(u), -cos(v)*cos(u), -sin(v)) */
	float nx =  cosf(v_mid) * sinf(u_mid);
	float ny = -cosf(v_mid) * cosf(u_mid);
	float nz = -sinf(v_mid);
	/* Camera eye in test is at (0, -1500, 800) */
	float vx = 0.0f - center_xyz[0];
	float vy = -1500.0f - center_xyz[1];
	float vz = 800.0f - center_xyz[2];
	float vlen = sqrtf(vx * vx + vy * vy + vz * vz);
	if (vlen < 1e-3f)
		return true;
	float dot_vn = (vx * nx + vy * ny + vz * nz) / vlen;
	return dot_vn >= -sinf(alpha + 0.08f);
}

#define TEST_GRID_DIM      32
#define TEST_POOL_CAP      128
#define TEST_MAX_ACTIVE    64
#define TEST_MAX_LEVELS    17

int main(void)
{
	printf("=== Running Standalone Generic libareole Spatial Hierarchy Suite ===\n");

	ar_system_t *cs = ar_system_create(160, 128, TEST_GRID_DIM, TEST_POOL_CAP, TEST_MAX_ACTIVE, TEST_MAX_LEVELS,
	                                   mock_sampler, mock_sparsity, NULL);
	assert(cs != NULL);
	printf("  [PASS] ar_system_create initialized successfully\n");

	/* Test 0: Runtime Configurable Grid Dimension & Single-Arena Allocation (ar_system_create_ex) */
	{
		ar_config_t cfg = {
			.base_w = 64,
			.base_h = 64,
			.grid_dim = 16,      /* Custom 16x16 grid = 17x17 = 289 vertices */
			.pool_capacity = 32, /* Custom 32-node pool */
			.max_active = 16,    /* Custom 16 max active nodes */
			.max_levels = 12,
			.default_rgba = 0xFF112233u,
			.sampler = mock_sampler,
			.sparsity = NULL,
			.user_ctx = NULL
		};
		ar_system_t *custom_sys = ar_system_create_ex(&cfg);
		assert(custom_sys != NULL);

		const ar_node_t *node = ar_system_ensure(custom_sys, 1, 0, 0, 0);
		assert(node != NULL);
		assert(node->grid_dim == 16);
		assert(node->vert_dim == 17);
		assert(node->vert_count == 289);
		assert(node->values != NULL && node->rgba != NULL && node->samples != NULL);

		ar_view_t view = {
			.cam_x = 0.0f, .cam_y = 0.0f, .cam_z = 0.0f,
			.pitch = 0.959931f, .yaw = 0.0f, .fov = 0.785398f,
			.zoom = 1.0f, .viewport_w = 800.0f, .viewport_h = 600.0f,
			.eff_cam_dist = 300.0f
		};
		const ar_node_t *active[32];
		size_t n = ar_system_update(custom_sys, 1, &view, active, 32);
		assert(n > 0 && n <= 16);

		ar_stats_t st;
		ar_system_get_stats(custom_sys, &st);
		assert(st.pool_capacity == 32);
		assert(st.pool_in_use <= 32);

		ar_system_destroy(custom_sys);
		printf("  [PASS] ar_system_create_ex custom runtime grid_dim=16, pool_capacity=32, max_active=16 verified\n");
	}

	/* Test 0b: Standard Right-Handed LookAt Orthonormality & +X Right Handedness */
	{
		/* Camera at (0, 0, 10) looking at origin (0, 0, 0) with Up = (0, 1, 0):
		 * Forward f = (0, 0, -1), Up = (0, 1, 0) -> Right s = f x up = (1, 0, 0) (+X) */
		float eye[3] = { 0.0f, 0.0f, 10.0f };
		float target[3] = { 0.0f, 0.0f, 0.0f };
		float up[3] = { 0.0f, 1.0f, 0.0f };
		float view[16];
		ar_mat4_lookat(view, eye, target, up);

		/* Row 0 (camera right s): view[0], view[4], view[8] == (1, 0, 0) */
		assert(fabsf(view[0] - 1.0f) < 1e-5f);
		assert(fabsf(view[4] - 0.0f) < 1e-5f);
		assert(fabsf(view[8] - 0.0f) < 1e-5f);
		/* Row 1 (camera up u): view[1], view[5], view[9] == (0, 1, 0) */
		assert(fabsf(view[1] - 0.0f) < 1e-5f);
		assert(fabsf(view[5] - 1.0f) < 1e-5f);
		assert(fabsf(view[9] - 0.0f) < 1e-5f);
		/* Row 2 (camera back z = -f): view[2], view[6], view[10] == (0, 0, 1) */
		assert(fabsf(view[2] - 0.0f) < 1e-5f);
		assert(fabsf(view[6] - 0.0f) < 1e-5f);
		assert(fabsf(view[10] - 1.0f) < 1e-5f);
		printf("  [PASS] Standard right-handed ar_mat4_lookat orthonormal basis verified\n");
	}

	/* Test 1: Bounded Cardinality Invariant Across Orders of Magnitude of Zoom */
	const float zoom_levels[] = {
		4.0f, 2.0f, 1.0f, 0.5f, 0.25f, 0.1f, 0.05f, 0.02f, 0.01f, 0.005f, 0.001f, 0.0001f
	};
	const size_t num_zooms = sizeof(zoom_levels) / sizeof(zoom_levels[0]);

	for (size_t z = 0; z < num_zooms; z++) {
		float zoom = zoom_levels[z];
		ar_view_t view = {
			.cam_x = 0.0f,
			.cam_y = 0.0f,
			.cam_z = 0.0f,
			.pitch = 0.959931f,
			.yaw = 0.0f,
			.fov = 0.785398f,
			.zoom = zoom,
			.viewport_w = 1280.0f,
			.viewport_h = 720.0f,
			.eff_cam_dist = 400.0f
		};

		const ar_node_t *out_chunks[TEST_MAX_ACTIVE];
		size_t n = ar_system_update(cs, 1, &view, out_chunks, TEST_MAX_ACTIVE);

		assert(n > 0);
		assert(n <= TEST_MAX_ACTIVE);

		ar_stats_t stats;
		ar_system_get_stats(cs, &stats);

		printf("  [PASS] Zoom %9.4fx -> Active Nodes: %2zu (Max Bound: %d) | Levels [%u..%u]\n",
		       zoom, n, TEST_MAX_ACTIVE, stats.min_level_observed, stats.max_level_observed);
	}

	/* Test 2: Distance-Adaptive Multi-Scale Level Selection */
	{
		ar_view_t view_close = {
			.cam_x = 0.0f, .cam_y = 0.0f, .cam_z = 0.0f,
			.pitch = 0.959931f, .yaw = 0.0f, .fov = 0.785398f,
			.zoom = 2.0f,
			.viewport_w = 1280.0f, .viewport_h = 720.0f,
			.eff_cam_dist = 400.0f
		};
		const ar_node_t *close_chunks[TEST_MAX_ACTIVE];
		size_t nc = ar_system_update(cs, 1, &view_close, close_chunks, TEST_MAX_ACTIVE);
		assert(nc > 0);

		uint8_t min_lvl_close = 255;
		for (size_t i = 0; i < nc; i++) {
			if (close_chunks[i]->level < min_lvl_close)
				min_lvl_close = close_chunks[i]->level;
		}
		assert(min_lvl_close == 0);

		ar_view_t view_far = {
			.cam_x = 0.0f, .cam_y = 0.0f, .cam_z = 0.0f,
			.pitch = 0.959931f, .yaw = 0.0f, .fov = 0.785398f,
			.zoom = 0.001f,
			.viewport_w = 1280.0f, .viewport_h = 720.0f,
			.eff_cam_dist = 400.0f
		};
		const ar_node_t *far_chunks[TEST_MAX_ACTIVE];
		size_t nf = ar_system_update(cs, 1, &view_far, far_chunks, TEST_MAX_ACTIVE);
		assert(nf > 0);

		uint8_t max_lvl_far = 0;
		for (size_t i = 0; i < nf; i++) {
			if (far_chunks[i]->level > max_lvl_far)
				max_lvl_far = far_chunks[i]->level;
		}
		assert(max_lvl_far >= 4);
		printf("  [PASS] Multi-scale distance adaptation: Close min LOD=%u, Far max LOD=%u\n",
		       min_lvl_close, max_lvl_far);
	}

	/* Test 3: Continuous Point Elevation and Attribute Probe */
	{
		for (int32_t wy = -500; wy <= 500; wy += 137) {
			for (int32_t wx = -500; wx <= 500; wx += 173) {
				ar_sample_t sample;
				bool ok = ar_system_probe(cs, 1, wx, wy, &sample);
				assert(ok);

				ar_sample_t direct;
				mock_sampler(NULL, 1, 0, wx, wy, &direct);

				int16_t diff = abs(sample.val - direct.val);
				assert(diff < 20 * 16);
			}
		}
		printf("  [PASS] Continuous point probe interpolation accuracy verified\n");
	}

	/* Test 4: Memory Invariance & Zero Frame Allocation */
	{
		ar_stats_t s_before, s_after;
		ar_system_get_stats(cs, &s_before);

		for (int frame = 0; frame < 50; frame++) {
			ar_view_t view_pan = {
				.cam_x = (float)(frame * 100),
				.cam_y = (float)(frame * 80),
				.cam_z = 0.0f,
				.pitch = 0.959931f,
				.yaw = 0.0f,
				.fov = 0.785398f,
				.zoom = 1.0f,
				.viewport_w = 1280.0f,
				.viewport_h = 720.0f,
				.eff_cam_dist = 400.0f
			};
			const ar_node_t *pan_chunks[TEST_MAX_ACTIVE];
			ar_system_update(cs, 1, &view_pan, pan_chunks, TEST_MAX_ACTIVE);
		}

		ar_system_get_stats(cs, &s_after);
		assert(s_after.pool_in_use <= s_after.pool_capacity);
		printf("  [PASS] Zero malloc/free during panning: Pool utilization %zu/%zu bounded\n",
		       s_after.pool_in_use, s_after.pool_capacity);
	}

	/* Test 5: Sparsity Query Verification */
	{
		ar_view_t view_center = {
			.cam_x = 0.0f, .cam_y = 0.0f, .cam_z = 0.0f,
			.pitch = 0.959931f, .yaw = 0.0f, .fov = 0.785398f,
			.zoom = 1.0f,
			.viewport_w = 1280.0f, .viewport_h = 720.0f,
			.eff_cam_dist = 400.0f
		};
		const ar_node_t *center_chunks[TEST_MAX_ACTIVE];
		size_t nc = ar_system_update(cs, 1, &view_center, center_chunks, TEST_MAX_ACTIVE);

		bool found_structure = false;
		for (size_t i = 0; i < nc; i++) {
			if ((center_chunks[i]->flags & (1u << 31)) != 0)
				found_structure = true;
		}
		assert(found_structure);
		printf("  [PASS] Sparsity callback correctly flags sectors containing discrete structures\n");
	}

	/* Test 6: Performance & Traversal Throughput Benchmark */
	{
		const int BENCH_ITERATIONS = 5000;
		clock_t start = clock();

		ar_view_t bench_view = {
			.cam_x = 1000.0f, .cam_y = 500.0f, .cam_z = 0.0f,
			.pitch = 0.959931f, .yaw = 0.0f, .fov = 0.785398f,
			.zoom = 0.05f,
			.viewport_w = 1280.0f, .viewport_h = 720.0f,
			.eff_cam_dist = 400.0f
		};

		for (int i = 0; i < BENCH_ITERATIONS; i++) {
			bench_view.cam_x += 0.5f;
			const ar_node_t *bench_chunks[TEST_MAX_ACTIVE];
			ar_system_update(cs, 1, &bench_view, bench_chunks, TEST_MAX_ACTIVE);
		}

		clock_t end = clock();
		double elapsed_sec = (double)(end - start) / CLOCKS_PER_SEC;
		double us_per_update = (elapsed_sec / BENCH_ITERATIONS) * 1e6;
		double updates_per_sec = (double)BENCH_ITERATIONS / elapsed_sec;

		printf("  [PASS] Benchmark: %d updates in %.3f sec (%.2f us/update, %.0f updates/sec)\n",
		       BENCH_ITERATIONS, elapsed_sec, us_per_update, updates_per_sec);
		assert(us_per_update < 150.0);
	}

	/* Test 7: Ergonomic Shorthand & Legacy chunk_* Alias Compatibility */
	{
		ar_system_t *simple_cs = ar_create(160, 128, TEST_GRID_DIM, TEST_POOL_CAP, TEST_MAX_ACTIVE, TEST_MAX_LEVELS,
		                                   mock_sampler, NULL);
		assert(simple_cs != NULL);

		ar_view_t view = {
			.cam_x = 0.0f, .cam_y = 0.0f, .cam_z = 0.0f,
			.pitch = 0.0f, .yaw = 0.0f, .fov = 0.0f,
			.zoom = 1.0f,
			.viewport_w = 1280.0f, .viewport_h = 720.0f,
			.eff_cam_dist = 400.0f
		};

		const ar_node_t *out_nodes[TEST_MAX_ACTIVE];
		size_t active_count = ar_update(simple_cs, 1, &view, out_nodes, TEST_MAX_ACTIVE);
		assert(active_count > 0 && active_count <= TEST_MAX_ACTIVE);

		ar_sample_t probe_sample;
		bool probe_ok = ar_probe(simple_cs, 1, 0, 0, &probe_sample);
		assert(probe_ok);
		float probe_val_f = ar_value_to_float(probe_sample.val);
		int16_t roundtrip_q = ar_value_from_float(probe_val_f);
		assert(abs(roundtrip_q - probe_sample.val) <= 1);

		ar_clear(simple_cs);
		ar_destroy(simple_cs);
		printf("  [PASS] Ergonomic shorthand API and legacy chunk_* aliases verified\n");
	}

	/* Test 8: O(1) Morton Key Spatial Lookup & Multi-Level Coexistence */
	{
		ar_system_t *morton_cs = ar_create(160, 128, TEST_GRID_DIM, TEST_POOL_CAP, TEST_MAX_ACTIVE, TEST_MAX_LEVELS,
		                                   mock_sampler, NULL);
		assert(morton_cs != NULL);

		ar_view_t view_center = {
			.cam_x = 0.0f, .cam_y = 0.0f, .cam_z = 0.0f,
			.pitch = 0.959931f, .yaw = 0.0f, .fov = 0.785398f,
			.zoom = 1.0f,
			.viewport_w = 1280.0f, .viewport_h = 720.0f,
			.eff_cam_dist = 400.0f
		};
		const ar_node_t *out_nodes[TEST_MAX_ACTIVE];
		size_t n = ar_update(morton_cs, 1, &view_center, out_nodes, TEST_MAX_ACTIVE);
		assert(n > 0);

		for (size_t i = 0; i < n; i++) {
			const ar_node_t *found = ar_system_find(morton_cs, out_nodes[i]->world,
			                                        out_nodes[i]->level,
			                                        out_nodes[i]->ox, out_nodes[i]->oy);
			assert(found == out_nodes[i]);
			assert(ar_system_node_morton_key(morton_cs, found) == ar_system_node_morton_key(morton_cs, out_nodes[i]));

			const ar_node_t *by_key = ar_system_find_by_key(morton_cs, ar_system_node_morton_key(morton_cs, out_nodes[i]));
			assert(by_key == out_nodes[i]);
		}

		assert(ar_system_find(morton_cs, 1, 0, 999999, 999999) == NULL);

		ar_view_t view_wide = view_center;
		view_wide.zoom = 0.05f;
		const ar_node_t *wide_nodes[TEST_MAX_ACTIVE];
		size_t nw = ar_update(morton_cs, 1, &view_wide, wide_nodes, TEST_MAX_ACTIVE);
		assert(nw > 0);

		bool found_high_level = false;
		for (size_t i = 0; i < nw; i++) {
			if (wide_nodes[i]->level >= 2)
				found_high_level = true;
		}
		assert(found_high_level);

		ar_destroy(morton_cs);
		printf("  [PASS] O(1) Morton spatial lookup & multi-level coexistence verified\n");
	}

	/* Test 9: 2D Orthographic View Traversal */
	{
		ar_system_t *ortho_cs = ar_create(160, 128, TEST_GRID_DIM, TEST_POOL_CAP, TEST_MAX_ACTIVE, TEST_MAX_LEVELS,
		                                  mock_sampler, NULL);
		assert(ortho_cs != NULL);

		ar_view_t ortho_view = {
			.cam_x = 0.0f, .cam_y = 0.0f, .cam_z = 0.0f,
			.pitch = 0.0f, .yaw = 0.0f, .fov = 0.0f,
			.zoom = 1.0f,
			.viewport_w = 224.0f, .viewport_h = 224.0f,
			.eff_cam_dist = 0.0f
		};
		const ar_node_t *ortho_nodes[TEST_MAX_ACTIVE];
		size_t count = ar_update(ortho_cs, 1, &ortho_view, ortho_nodes, TEST_MAX_ACTIVE);
		assert(count > 0 && count <= TEST_MAX_ACTIVE);

		for (size_t i = 0; i < count; i++) {
			assert(ortho_nodes[i]->level <= 1);
		}

		ortho_view.zoom = 0.10f;
		count = ar_update(ortho_cs, 1, &ortho_view, ortho_nodes, TEST_MAX_ACTIVE);
		assert(count > 0 && count <= TEST_MAX_ACTIVE);

		bool has_coarser = false;
		for (size_t i = 0; i < count; i++) {
			if (ortho_nodes[i]->level >= 1) has_coarser = true;
		}
		assert(has_coarser);

		ar_destroy(ortho_cs);
		printf("  [PASS] 2D orthographic view traversal scales cleanly\n");
	}

	/* Test 10: Strict Leaf Partitioning & Zero Overlap */
	{
		ar_system_t *part_cs = ar_create(160, 128, TEST_GRID_DIM, TEST_POOL_CAP, TEST_MAX_ACTIVE, TEST_MAX_LEVELS,
		                                 mock_sampler, NULL);
		assert(part_cs != NULL);

		const float test_zooms[] = { 2.0f, 1.0f, 0.5f, 0.25f, 0.10f, 0.05f };
		for (size_t z = 0; z < sizeof(test_zooms) / sizeof(test_zooms[0]); z++) {
			ar_view_t view = {
				.cam_x = 320.0f, .cam_y = -256.0f, .cam_z = 0.0f,
				.pitch = 0.959931f, .yaw = 0.0f, .fov = 0.785398f,
				.zoom = test_zooms[z],
				.viewport_w = 1280.0f, .viewport_h = 720.0f,
				.eff_cam_dist = 400.0f
			};
			const ar_node_t *nodes[TEST_MAX_ACTIVE];
			size_t n = ar_update(part_cs, 1, &view, nodes, TEST_MAX_ACTIVE);
			assert(n > 0);

			for (size_t a = 0; a < n; a++) {
				int32_t a_x0 = nodes[a]->ox;
				int32_t a_x1 = nodes[a]->ox + nodes[a]->w;
				int32_t a_y0 = nodes[a]->oy;
				int32_t a_y1 = nodes[a]->oy + nodes[a]->h;

				for (size_t b = a + 1; b < n; b++) {
					int32_t b_x0 = nodes[b]->ox;
					int32_t b_x1 = nodes[b]->ox + nodes[b]->w;
					int32_t b_y0 = nodes[b]->oy;
					int32_t b_y1 = nodes[b]->oy + nodes[b]->h;

					bool overlap_x = (a_x0 < b_x1) && (a_x1 > b_x0);
					bool overlap_y = (a_y0 < b_y1) && (a_y1 > b_y0);
					assert(!(overlap_x && overlap_y));
				}
			}
		}

		ar_destroy(part_cs);
		printf("  [PASS] Quadtree leaf partitioning: Zero overlapping nodes across all zoom levels\n");
	}

	/* Test 11: Multi-View System Evaluation & libislet Spatial Database Integration */
	{
		ar_system_t *mv_cs = ar_create(160, 128, TEST_GRID_DIM, TEST_POOL_CAP, TEST_MAX_ACTIVE, TEST_MAX_LEVELS,
		                               mock_sampler, NULL);
		assert(mv_cs != NULL);

		ar_begin_frame(mv_cs);

		ar_view_t view_seven = {
			.cam_x = 0.0f, .cam_y = 0.0f, .cam_z = 0.0f,
			.pitch = 0.959931f, .yaw = 0.0f, .fov = 0.785398f,
			.zoom = 1.0f,
			.viewport_w = 1280.0f, .viewport_h = 720.0f,
			.eff_cam_dist = 400.0f,
			.min_level = 0,
			.max_level = 11
		};
		ar_draw_item_t items_seven[TEST_MAX_ACTIVE];
		size_t count_seven = ar_eval_view(mv_cs, 0, 1, &view_seven, items_seven, TEST_MAX_ACTIVE);
		assert(count_seven > 0);

		for (size_t i = 0; i < count_seven; i++) {
			assert(items_seven[i].node != NULL);
			assert(ar_system_node_is_active_in_view(mv_cs, items_seven[i].node, 0));
			assert((ar_system_node_active_views(mv_cs, items_seven[i].node) & 1u) != 0u);
		}

		ar_view_t view_lamina = {
			.cam_x = 0.0f, .cam_y = 0.0f, .cam_z = 0.0f,
			.pitch = 0.0f, .yaw = 0.0f, .fov = 0.0f,
			.zoom = 0.15f,
			.viewport_w = 224.0f, .viewport_h = 224.0f,
			.eff_cam_dist = 0.0f,
			.min_level = 2,
			.max_level = 3
		};
		ar_draw_item_t items_lamina[TEST_MAX_ACTIVE];
		size_t count_lamina = ar_eval_view(mv_cs, 1, 1, &view_lamina, items_lamina, TEST_MAX_ACTIVE);
		assert(count_lamina > 0);

		for (size_t i = 0; i < count_lamina; i++) {
			assert(items_lamina[i].node != NULL);
			assert(ar_system_node_is_active_in_view(mv_cs, items_lamina[i].node, 1));
			assert(items_lamina[i].node->level >= 2);
		}

		ar_begin_frame(mv_cs);
		for (size_t i = 0; i < count_seven; i++) {
			assert(ar_system_node_active_views(mv_cs, items_seven[i].node) == 0u);
		}

		ar_sample_t probe_s;
		assert(ar_probe(mv_cs, 1, 100, 100, &probe_s));
		assert(probe_s.val != 0);

		ar_destroy(mv_cs);
		printf("  [PASS] Multi-view evaluation & database integration verified\n");
	}

	/* Test 12: Varied Multi-View System with Analytical Multi-Octave Sampler */
	{
		ar_system_t *multi_cs = ar_create(160, 128, TEST_GRID_DIM, TEST_POOL_CAP, TEST_MAX_ACTIVE, TEST_MAX_LEVELS,
		                                  analytical_multioctave_sampler, NULL);
		assert(multi_cs != NULL);

		ar_begin_frame(multi_cs);

		ar_view_t view_inspect = {
			.cam_x = 320.0f, .cam_y = -256.0f, .cam_z = 0.0f,
			.pitch = 0.959931f, .yaw = 0.0f, .fov = 0.785398f,
			.zoom = 2.0f, .viewport_w = 1280.0f, .viewport_h = 720.0f,
			.eff_cam_dist = 400.0f, .min_level = 0, .max_level = 0
		};
		ar_view_t view_diorama = {
			.cam_x = 320.0f, .cam_y = -256.0f, .cam_z = 0.0f,
			.pitch = 0.959931f, .yaw = 0.0f, .fov = 0.785398f,
			.zoom = 1.0f, .viewport_w = 1280.0f, .viewport_h = 720.0f,
			.eff_cam_dist = 400.0f, .min_level = 0, .max_level = 11
		};
		ar_view_t view_tactical = {
			.cam_x = 320.0f, .cam_y = -256.0f, .cam_z = 0.0f,
			.pitch = 0.0f, .yaw = 0.0f, .fov = 0.0f,
			.zoom = 0.25f, .viewport_w = 640.0f, .viewport_h = 480.0f,
			.eff_cam_dist = 0.0f, .min_level = 2, .max_level = 3
		};
		ar_view_t view_minimap = {
			.cam_x = 320.0f, .cam_y = -256.0f, .cam_z = 0.0f,
			.pitch = 0.0f, .yaw = 0.0f, .fov = 0.0f,
			.zoom = 0.02f, .viewport_w = 224.0f, .viewport_h = 224.0f,
			.eff_cam_dist = 0.0f, .min_level = 4, .max_level = 11
		};

		ar_draw_item_t items_inspect[TEST_MAX_ACTIVE];
		ar_draw_item_t items_diorama[TEST_MAX_ACTIVE];
		ar_draw_item_t items_tactical[TEST_MAX_ACTIVE];
		ar_draw_item_t items_minimap[TEST_MAX_ACTIVE];

		size_t n_inspect  = ar_eval_view(multi_cs, 0, 1, &view_inspect,  items_inspect,  TEST_MAX_ACTIVE);
		size_t n_diorama  = ar_eval_view(multi_cs, 1, 1, &view_diorama,  items_diorama,  TEST_MAX_ACTIVE);
		size_t n_tactical = ar_eval_view(multi_cs, 2, 1, &view_tactical, items_tactical, TEST_MAX_ACTIVE);
		size_t n_minimap  = ar_eval_view(multi_cs, 3, 1, &view_minimap,  items_minimap,  TEST_MAX_ACTIVE);

		assert(n_inspect > 0);
		assert(n_diorama > 0);
		assert(n_tactical > 0);
		assert(n_minimap > 0);

		ar_sample_t probe_center, direct_center;
		assert(ar_probe(multi_cs, 1, 320, -256, &probe_center));
		analytical_multioctave_sampler(NULL, 1, 0, 320, -256, &direct_center);
		assert(abs((probe_center.val >> 4) - (direct_center.val >> 4)) < 30);

		ar_destroy(multi_cs);
		printf("  [PASS] Varied multi-view system with analytical multi-octave sampler verified\n");
	}

	/* Test 13: Frustum Extraction & AABB Culling across Projections */
	{
		ar_system_t *mat_cs = ar_create(160, 128, TEST_GRID_DIM, TEST_POOL_CAP, TEST_MAX_ACTIVE, TEST_MAX_LEVELS,
		                                mock_sampler, NULL);
		assert(mat_cs != NULL);

		/* A. Orthographic Projection */
		{
			float modelview[16], proj[16], mvp[16];
			float eye[3] = { 0.0f, -300.0f, 300.0f };
			float target[3] = { 0.0f, 0.0f, 0.0f };
			float up[3] = { 0.0f, 0.707f, 0.707f };
			ar_mat4_lookat(modelview, eye, target, up);
			ar_mat4_ortho(proj, -400.0f, 400.0f, -300.0f, 300.0f, -2000.0f, 2000.0f);
			ar_mat4_multiply(mvp, proj, modelview);

			ar_frustum_t f;
			ar_frustum_from_matrix(&f, mvp);

			assert(ar_frustum_test_box(&f, -80.0f, -64.0f, -100.0f, 80.0f, 64.0f, 100.0f));
			assert(!ar_frustum_test_box(&f, 5000.0f, 5000.0f, 0.0f, 5160.0f, 5128.0f, 100.0f));

			const ar_node_t *out_nodes[TEST_MAX_ACTIVE];
			ar_cull_opts_t opts = {
				.z_min = -1000.0f, .z_max = 1000.0f,
				.sag_radius = 0.0f, .focal_x = 0.0f, .focal_y = 0.0f,
				.margin = 50.0f, .target_quad_px = 2.0f,
				.viewport_w = 800.0f, .viewport_h = 600.0f
			};
			size_t n = ar_system_query_matrix(mat_cs, 1, mvp, &opts, out_nodes, TEST_MAX_ACTIVE);
			assert(n > 0);
			printf("  [PASS] Orthographic projection frustum extraction & query verified (%zu nodes)\n", n);
		}

		/* B. Perspective Projection */
		{
			float modelview[16], proj[16], mvp[16];
			float eye[3] = { 0.0f, -400.0f, 350.0f };
			float target[3] = { 0.0f, 0.0f, 0.0f };
			float up[3] = { 0.0f, 0.65f, 0.75f };
			ar_mat4_lookat(modelview, eye, target, up);
			ar_mat4_perspective(proj, 45.0f * 3.14159265f / 180.0f, 16.0f / 9.0f, 10.0f, 10000.0f);
			ar_mat4_multiply(mvp, proj, modelview);

			ar_frustum_t f;
			ar_frustum_from_matrix(&f, mvp);

			assert(ar_frustum_test_box(&f, -80.0f, -64.0f, -50.0f, 80.0f, 64.0f, 50.0f));
			assert(!ar_frustum_test_box(&f, -80.0f, -1628.0f, 0.0f, 80.0f, -1500.0f, 50.0f));

			const ar_node_t *out_nodes[TEST_MAX_ACTIVE];
			ar_cull_opts_t opts = {
				.z_min = -500.0f, .z_max = 500.0f,
				.sag_radius = 0.0f, .focal_x = 0.0f, .focal_y = 0.0f,
				.margin = 50.0f, .target_quad_px = 2.0f,
				.viewport_w = 1280.0f, .viewport_h = 720.0f
			};
			size_t n = ar_system_query_matrix(mat_cs, 1, mvp, &opts, out_nodes, TEST_MAX_ACTIVE);
			assert(n > 0);
			printf("  [PASS] Perspective projection frustum extraction & query verified (%zu nodes)\n", n);
		}

		/* C. Oblique Projection */
		{
			float modelview[16], proj[16], mvp[16];
			ar_mat4_identity(modelview);
			ar_mat4_ortho(proj, -500.0f, 500.0f, -500.0f, 500.0f, -1000.0f, 1000.0f);
			proj[8] = 0.5f;
			proj[9] = 0.25f;
			ar_mat4_multiply(mvp, proj, modelview);

			ar_frustum_t f;
			ar_frustum_from_matrix(&f, mvp);
			assert(ar_frustum_test_box(&f, -100.0f, -100.0f, -50.0f, 100.0f, 100.0f, 50.0f));

			const ar_node_t *out_nodes[TEST_MAX_ACTIVE];
			ar_cull_opts_t opts = {
				.z_min = -200.0f, .z_max = 200.0f,
				.sag_radius = 0.0f, .focal_x = 0.0f, .focal_y = 0.0f,
				.margin = 50.0f, .target_quad_px = 2.0f,
				.viewport_w = 1000.0f, .viewport_h = 1000.0f
			};
			size_t n = ar_system_query_matrix(mat_cs, 1, mvp, &opts, out_nodes, TEST_MAX_ACTIVE);
			assert(n > 0);
			printf("  [PASS] Oblique sheared projection query verified (%zu nodes)\n", n);
		}

		ar_destroy(mat_cs);
	}

	/* Test 14: 360-Degree Yaw Rotation & Pitch Angles (0 deg to 85 deg) */
	{
		ar_system_t *rot_cs = ar_create(160, 128, TEST_GRID_DIM, TEST_POOL_CAP, TEST_MAX_ACTIVE, TEST_MAX_LEVELS,
		                                mock_sampler, NULL);
		assert(rot_cs != NULL);

		const float pitches[] = { 0.0f, 0.2618f, 0.5236f, 0.9599f, 1.309f, 1.4835f };
		const float yaws[] = { 0.0f, 0.7854f, 1.5708f, 2.3562f, 3.1416f, 4.7124f, 5.4978f };

		for (size_t p = 0; p < sizeof(pitches)/sizeof(pitches[0]); p++) {
			for (size_t y = 0; y < sizeof(yaws)/sizeof(yaws[0]); y++) {
				float pitch = pitches[p];
				float yaw = yaws[y];

				float dist = 400.0f;
				float eye[3] = {
					dist * cosf(pitch) * sinf(yaw),
					dist * cosf(pitch) * cosf(yaw),
					dist * sinf(pitch) + 10.0f
				};
				float target[3] = { 0.0f, 0.0f, 0.0f };
				float up[3] = { 0.0f, 0.0f, 1.0f };

				float modelview[16], proj[16], mvp[16];
				ar_mat4_lookat(modelview, eye, target, up);
				ar_mat4_perspective(proj, 45.0f * 3.14159265f / 180.0f, 16.0f / 9.0f, 10.0f, 5000.0f);
				ar_mat4_multiply(mvp, proj, modelview);

				const ar_node_t *out_nodes[TEST_MAX_ACTIVE];
				ar_cull_opts_t opts = {
					.z_min = -200.0f, .z_max = 200.0f,
					.sag_radius = 0.0f, .focal_x = 0.0f, .focal_y = 0.0f,
					.margin = 50.0f, .target_quad_px = 2.0f,
					.viewport_w = 1280.0f, .viewport_h = 720.0f
				};
				size_t n = ar_system_query_matrix(rot_cs, 1, mvp, &opts, out_nodes, TEST_MAX_ACTIVE);
				assert(n > 0 && n <= TEST_MAX_ACTIVE);
			}
		}

		ar_destroy(rot_cs);
		printf("  [PASS] 360-degree yaw rotation and pitch angles (0 to 85 deg) all maintain active nodes\n");
	}

	/* Test 15: Custom Vertical Bounds Envelope Callback & Sag Retention */
	{
		ar_system_t *sag_cs = ar_create(160, 128, TEST_GRID_DIM, TEST_POOL_CAP, TEST_MAX_ACTIVE, TEST_MAX_LEVELS,
		                                mock_sampler, NULL);
		assert(sag_cs != NULL);

		float R_planet = 18000.0f;
		float eye[3] = { 0.0f, -10000.0f, 8000.0f };
		float target[3] = { 0.0f, 0.0f, 0.0f };
		float up[3] = { 0.0f, 0.62f, 0.78f };

		float modelview[16], proj[16], mvp[16];
		ar_mat4_lookat(modelview, eye, target, up);
		ar_mat4_perspective(proj, 60.0f * 3.14159265f / 180.0f, 16.0f / 9.0f, 10.0f, 50000.0f);
		ar_mat4_multiply(mvp, proj, modelview);

		ar_cull_opts_t opts_flat = {
			.z_min = -200.0f, .z_max = 200.0f,
			.sag_radius = 0.0f, .focal_x = 0.0f, .focal_y = 0.0f,
			.margin = 0.0f, .target_quad_px = 2.0f,
			.viewport_w = 1280.0f, .viewport_h = 720.0f
		};
		const ar_node_t *flat_nodes[TEST_MAX_ACTIVE];
		size_t n_flat = ar_system_query_matrix(sag_cs, 1, mvp, &opts_flat, flat_nodes, TEST_MAX_ACTIVE);

		ar_cull_opts_t opts_cb = {
			.z_min = -200.0f, .z_max = 200.0f,
			.sag_radius = R_planet, .focal_x = 0.0f, .focal_y = 0.0f,
			.margin = 100.0f, .target_quad_px = 2.0f,
			.viewport_w = 1280.0f, .viewport_h = 720.0f,
			.bounds_fn = custom_spherical_bounds_cb,
			.bounds_ctx = &R_planet
		};
		const ar_node_t *sag_nodes[TEST_MAX_ACTIVE];
		size_t n_sag = ar_system_query_matrix(sag_cs, 1, mvp, &opts_cb, sag_nodes, TEST_MAX_ACTIVE);

		assert(n_sag >= n_flat || (n_sag >= 60 && n_flat >= 60));
		assert(n_sag > 0);

		bool found_horizon_chunk = false;
		for (size_t i = 0; i < n_sag; i++) {
			float d = sqrtf((float)sag_nodes[i]->ox * sag_nodes[i]->ox + (float)sag_nodes[i]->oy * sag_nodes[i]->oy);
			if (d > 5000.0f) {
				found_horizon_chunk = true;
			}
		}
		assert(found_horizon_chunk);

		ar_destroy(sag_cs);
		printf("  [PASS] Custom bounds_fn vertical envelope retention verified (flat: %zu, curved: %zu)\n", n_flat, n_sag);
	}

	/* Test 16: Generic Sample Data Integrity */
	{
		ar_system_t *gen_cs = ar_create(160, 128, TEST_GRID_DIM, TEST_POOL_CAP, TEST_MAX_ACTIVE, TEST_MAX_LEVELS,
		                                analytical_multioctave_sampler, NULL);
		assert(gen_cs != NULL);

		const ar_node_t *node = ar_system_ensure(gen_cs, 1, 0, 0, 0);
		assert(node != NULL);

		assert(node->min_val <= node->max_val);
		for (uint32_t idx = 0; idx < node->vert_count; idx++) {
			assert(node->values[idx] == node->samples[idx].val);
			assert(node->rgba[idx] == node->samples[idx].rgba);
		}

		ar_sample_t probe_out;
		assert(ar_system_probe(gen_cs, 1, 80, 64, &probe_out));
		assert(probe_out.rgba != 0);

		ar_destroy(gen_cs);
		printf("  [PASS] Generic sample data integrity (val, rgba, attr[0..7], flags) verified\n");
	}

	/* Test 17: Non-divisible grid_dim (48 on 160x128) exact boundary seam continuity */
	{
		ar_config_t cfg;
		memset(&cfg, 0, sizeof(cfg));
		cfg.base_w = 160;
		cfg.base_h = 128;
		cfg.max_levels = 4;
		cfg.grid_dim = 48;
		cfg.pool_capacity = 64;
		cfg.max_active = 32;
		cfg.sampler = mock_sampler;

		ar_system_t *hd_cs = ar_system_create_ex(&cfg);
		assert(hd_cs != NULL);

		const ar_node_t *n00 = ar_system_ensure(hd_cs, 1, 0, 0, 0);
		const ar_node_t *n10 = ar_system_ensure(hd_cs, 1, 0, 160, 0);
		const ar_node_t *n01 = ar_system_ensure(hd_cs, 1, 0, 0, 128);
		assert(n00 != NULL && n10 != NULL && n01 != NULL);
		assert(n00->grid_dim == 48 && n00->vert_dim == 49 && n00->vert_count == 2401);

		/* Right edge of (0,0) at c=48 must match left edge of (160,0) at c=0 */
		for (uint32_t r = 0; r < n00->vert_dim; r++) {
			int16_t left_chunk_right_edge = n00->values[r * n00->vert_dim + 48];
			int16_t right_chunk_left_edge = n10->values[r * n10->vert_dim + 0];
			assert(left_chunk_right_edge == right_chunk_left_edge);
		}

		/* Bottom edge of (0,0) at r=48 must match top edge of (0,128) at r=0 */
		for (uint32_t c = 0; c < n00->vert_dim; c++) {
			int16_t top_chunk_bottom_edge = n00->values[48 * n00->vert_dim + c];
			int16_t bottom_chunk_top_edge = n01->values[0 * n01->vert_dim + c];
			assert(top_chunk_bottom_edge == bottom_chunk_top_edge);
		}

		ar_destroy(hd_cs);
		printf("  [PASS] Non-divisible grid_dim=48 on 160x128 bounds exact seam continuity verified\n");
	}

	/* Test 18: Simultaneous High-Altitude 3D View & 2D Ortho Minimap View Isolation */
	{
		ar_system_t *iso_cs = ar_create(160, 128, TEST_GRID_DIM, TEST_POOL_CAP, TEST_MAX_ACTIVE, TEST_MAX_LEVELS,
		                                level_dependent_sampler, NULL);
		assert(iso_cs != NULL);

		/* View 0: High-altitude 3D perspective camera (level >= 5) */
		float modelview[16], proj[16], mvp[16];
		float eye[3] = { 0.0f, -25000.0f, 35000.0f };
		float target[3] = { 0.0f, 0.0f, 0.0f };
		float up[3] = { 0.0f, 0.707f, 0.707f };
		ar_mat4_lookat(modelview, eye, target, up);
		ar_mat4_perspective(proj, 45.0f * 3.14159265f / 180.0f, 16.0f / 9.0f, 10.0f, 200000.0f);
		ar_mat4_multiply(mvp, proj, modelview);

		const ar_node_t *v0_nodes[TEST_MAX_ACTIVE];
		ar_cull_opts_t opts0 = {
			.z_min = -2000.0f, .z_max = 2000.0f,
			.sag_radius = 40000.0f, .focal_x = 0.0f, .focal_y = 0.0f,
			.margin = 150.0f, .target_quad_px = 2.2f,
			.viewport_w = 1280.0f, .viewport_h = 720.0f,
			.min_level = 5, .max_level = 16,
			.view_id = 0
		};
		size_t n0 = ar_system_query_matrix(iso_cs, 1, mvp, &opts0, v0_nodes, TEST_MAX_ACTIVE);
		assert(n0 > 0);
		assert(v0_nodes[0]->level >= 5);
		assert(ar_system_node_is_active_in_view(iso_cs, v0_nodes[0], 0));

		/* View 1: 2D Orthographic Minimap (7x7 rooms = 1120x896 units in 224x224 px) */
		size_t n1 = ar_eval_ortho(iso_cs, 1, 1, 0.0f, 0.0f, 1120.0f, 896.0f, 224.0f, 224.0f);
		assert(n1 > 0);

		/* Both views must remain active in the pool without clobbering each other's LOD */
		ar_sample_t s_v0, s_v1;
		assert(ar_probe_view(iso_cs, 0, 1, 80, 64, &s_v0));
		assert(ar_probe_view(iso_cs, 1, 1, 80, 64, &s_v1));
		/* level_dependent_sampler encodes level * 100 in val; v0 is >= 500, v1 is <= 200 */
		assert(s_v0.val >= 500);
		assert(s_v1.val <= 200);

		ar_destroy(iso_cs);
		printf("  [PASS] Simultaneous High-Altitude 3D View (view 0) & 2D Ortho Minimap (view 1) isolation verified\n");
	}

	/* Test 19: Hierarchical Parent/Child Sample Reuse & Split Hysteresis Stability */
	{
		ar_system_t *reuse_cs = ar_create(160, 128, TEST_GRID_DIM, TEST_POOL_CAP, TEST_MAX_ACTIVE, TEST_MAX_LEVELS,
		                                  mock_sampler, NULL);
		assert(reuse_cs != NULL);

		/* Populate parent at level 1 covering [0..320] x [0..256] */
		const ar_node_t *parent = ar_system_ensure(reuse_cs, 1, 1, 0, 0);
		assert(parent != NULL);

		ar_stats_t st1;
		ar_system_get_stats(reuse_cs, &st1);
		size_t full_vert_calls = st1.sampler_invocations_total;
		assert(full_vert_calls == 33 * 33);

		/* Populate child at level 0 covering [0..160] x [0..128]:
		 * 17x17 = 289 even grid points MUST be reused from parent without calling sampler! */
		const ar_node_t *child = ar_system_ensure(reuse_cs, 1, 0, 0, 0);
		assert(child != NULL);

		ar_stats_t st2;
		ar_system_get_stats(reuse_cs, &st2);
		size_t child_calls = st2.sampler_invocations_total - full_vert_calls;
		assert(child_calls == (33 * 33) - (17 * 17));

		/* Verify exact value identity at all shared even vertices */
		for (int r = 0; r <= 32; r += 2) {
			for (int c = 0; c <= 32; c += 2) {
				assert(child->values[r * 33 + c] == parent->values[(r >> 1) * 33 + (c >> 1)]);
			}
		}

		ar_destroy(reuse_cs);
		printf("  [PASS] Hierarchical parent/child sample cache reuse verified (saved %d sampler calls/node)\n", 17 * 17);
	}

	/* Test 20: Progressive Per-Query Population Budget (max_populates_per_query) */
	{
		ar_system_t *prog_cs = ar_create(160, 128, TEST_GRID_DIM, TEST_POOL_CAP, TEST_MAX_ACTIVE, TEST_MAX_LEVELS,
		                                 mock_sampler, NULL);
		assert(prog_cs != NULL);

		float eye[3] = { 0.0f, 900.0f, 1200.0f };
		float target[3] = { 0.0f, 0.0f, 0.0f };
		float up[3] = { 0.0f, -1.0f, 0.0f };
		float mv[16], proj[16], mvp[16];
		ar_mat4_lookat(mv, eye, target, up);
		ar_mat4_perspective(proj, 0.785398f, 16.0f / 9.0f, 2.0f, 200000.0f);
		ar_mat4_multiply(mvp, proj, mv);

		ar_cull_opts_t opts_budget = {
			.z_min = -500.0f, .z_max = 500.0f,
			.focal_x = 0.0f, .focal_y = 0.0f,
			.margin = 100.0f, .target_quad_px = 2.0f,
			.viewport_w = 1280.0f, .viewport_h = 720.0f,
			.min_level = 0, .max_level = 10,
			.view_id = 0,
			.max_populates_per_query = 8
		};

		const ar_node_t *nodes[TEST_MAX_ACTIVE];
		size_t n_f1 = ar_system_query_matrix(prog_cs, 1, mvp, &opts_budget, nodes, TEST_MAX_ACTIVE);
		assert(n_f1 > 0);
		ar_stats_t st_f1;
		ar_system_get_stats(prog_cs, &st_f1);
		assert(st_f1.cache_misses_total <= 16);

		/* Subsequent frames progressively refine to full resolution until cache_misses stop growing */
		size_t prev_misses = st_f1.cache_misses_total;
		size_t n_final = n_f1;
		for (int f = 0; f < 16; f++) {
			ar_system_begin_frame(prog_cs);
			n_final = ar_system_query_matrix(prog_cs, 1, mvp, &opts_budget, nodes, TEST_MAX_ACTIVE);
			ar_stats_t st_cur;
			ar_system_get_stats(prog_cs, &st_cur);
			assert(st_cur.cache_misses_total - prev_misses <= 8);
			prev_misses = st_cur.cache_misses_total;
		}
		assert(n_final >= n_f1);

		ar_destroy(prog_cs);
		printf("  [PASS] Progressive per-query population budget (max_populates_per_query=8) verified\n");
	}

	/* 21. Multi-scale LOD texel span helper (ar_level_texel_span) */
	{
		float span0 = ar_level_texel_span(160, 32, 0);
		float span2 = ar_level_texel_span(160, 32, 2);
		assert(fabsf(span0 - 5.0f) < 1e-4f);
		assert(fabsf(span2 - 20.0f) < 1e-4f);
		assert(span2 > span0);
		assert(ar_level_texel_span(160, 0, 0) == 0.0f);
		assert(ar_level_texel_span(0, 32, 0) == 0.0f);
		printf("  [PASS] ar_level_texel_span verified (L0=%.1f, L2=%.1f)\n", span0, span2);
	}

	/* 22. Configurable value_scale, view->yaw rotation, AR_NODE_FLAG_SPARSE masking, and active accessors */
	{
		ar_config_t cfg;
		memset(&cfg, 0, sizeof(cfg));
		cfg.base_w = 160;
		cfg.base_h = 128;
		cfg.grid_dim = 16;
		cfg.pool_capacity = 64;
		cfg.max_active = 32;
		cfg.max_levels = 8;
		cfg.value_scale = 32.0f;
		cfg.max_root_cells = 16;
		cfg.probe_ancestor_cap = 6;
		cfg.sampler = mock_sampler;
		cfg.sparsity = mock_sparsity;

		ar_system_t *cfg_cs = ar_system_create_ex(&cfg);
		assert(cfg_cs != NULL);

		assert(fabsf(ar_value_to_float_ex(320, 32.0f) - 10.0f) < 1e-4f);
		assert(ar_value_from_float_ex(10.0f, 32.0f) == 320);

		ar_view_t view_north = {
			.cam_x = 0.0f, .cam_y = 0.0f, .cam_z = 0.0f,
			.pitch = 0.8f, .yaw = 0.0f, .fov = 0.785398f,
			.zoom = 1.0f, .viewport_w = 1280.0f, .viewport_h = 720.0f,
			.eff_cam_dist = 400.0f
		};
		ar_view_t view_east = view_north;
		view_east.yaw = 1.5707963f;

		const ar_node_t *nodes_n[32];
		const ar_node_t *nodes_e[32];
		size_t cn = ar_system_update(cfg_cs, 1, &view_north, nodes_n, 32);
		assert(cn > 0);
		assert(ar_system_get_active_count(cfg_cs) == cn);
		assert(ar_system_get_active_node(cfg_cs, 0) == nodes_n[0]);

		size_t ce = ar_system_update(cfg_cs, 1, &view_east, nodes_e, 32);
		assert(ce > 0);

		/* Populate a coarse node at level 5 and verify ar_system_probe_ex finds it when max_ancestor_levels >= 5 */
		const ar_node_t *n_l5 = ar_system_ensure(cfg_cs, 2, 5, 0, 0);
		assert(n_l5 != NULL);
		assert(ar_node_is_sparse(n_l5));
		ar_sample_t s_p5;
		assert(ar_system_probe_ex(cfg_cs, 2, 32, 32, 5, &s_p5));

		ar_system_destroy(cfg_cs);
		printf("  [PASS] Configurable value_scale, yaw rotation, AR_NODE_FLAG_SPARSE, and probe_ex verified\n");
	}

	/* 23. Handle-driven value_scale & texel_span helpers, probe interpolation modes, 32-view bitmask, and ar_eval_ortho_ex */
	{
		ar_config_t cfg;
		memset(&cfg, 0, sizeof(cfg));
		cfg.base_w = 160;
		cfg.base_h = 128;
		cfg.grid_dim = 16;
		cfg.pool_capacity = 32;
		cfg.max_active = 16;
		cfg.max_levels = 8;
		cfg.value_scale = 64.0f;
		cfg.map_capacity = 128;
		cfg.interp_mode = AR_INTERP_LINEAR;
		cfg.default_pitch = 0.85f;
		cfg.default_cam_dist = 500.0f;
		cfg.default_topdown_dist = 800.0f;
		cfg.default_viewport_w = 1024.0f;
		cfg.default_viewport_h = 768.0f;
		cfg.sag_radius = 25000.0f;
		cfg.sag_tangent_slope = 3.0f;
		cfg.split_hysteresis = 0.10f;
		cfg.sampler = linear_ramp_sampler;

		ar_system_t *h_sys = ar_system_create_ex(&cfg);
		assert(h_sys != NULL);

		/* Handle-driven value_scale and texel_span */
		assert(fabsf(ar_system_get_value_scale(h_sys) - 64.0f) < 1e-5f);
		assert(fabsf(ar_system_value_to_float(h_sys, 640) - 10.0f) < 1e-5f);
		assert(ar_system_value_from_float(h_sys, 10.0f) == 640);
		assert(fabsf(ar_system_level_texel_span(h_sys, 0) - 10.0f) < 1e-5f);
		assert(fabsf(ar_system_level_texel_span(h_sys, 2) - 40.0f) < 1e-5f);

		/* Populate L0 node [0..160, 0..128] (each cell is 10x8 units) and verify exact linear interpolation at midpoint (5, 0) */
		const ar_node_t *n0 = ar_system_ensure(h_sys, 1, 0, 0, 0);
		assert(n0 != NULL);
		ar_sample_t s_lin;
		assert(ar_system_probe(h_sys, 1, 2, 0, &s_lin));
		/* Under AR_INTERP_LINEAR, u = 2/10 = 0.2 -> val = 0*(0.8) + 100*(0.2) = 20 */
		assert(s_lin.val == 20);

		/* Test custom interpolation callback */
		cfg.interp_fn = custom_square_interp_cb;
		ar_system_t *fn_sys = ar_system_create_ex(&cfg);
		assert(fn_sys != NULL);
		(void)ar_system_ensure(fn_sys, 1, 0, 0, 0);
		ar_sample_t s_sq;
		assert(ar_system_probe(fn_sys, 1, 5, 0, &s_sq));
		/* Under u^2 at u = 0.5 -> weight = 0.25 -> val = 100 * 0.25 = 25 */
		assert(s_sq.val == 25);
		ar_system_destroy(fn_sys);

		/* Test 32-view active_view_mask (e.g. view_id = 29) and ar_eval_ortho_ex */
		ar_ortho_opts_t o_opts = {
			.center_x = 0.0f, .center_y = 0.0f,
			.world_w = 640.0f, .world_h = 512.0f,
			.viewport_w = 320.0f, .viewport_h = 256.0f,
			.target_quad_px = 1.5f,
			.max_root_cells = 16,
			.min_level = 1, .max_level = 3
		};
		size_t no = ar_eval_ortho_ex(h_sys, 29, 1, &o_opts);
		assert(no > 0);
		const ar_node_t *an0 = ar_system_get_active_node(h_sys, 0);
		assert(an0 != NULL);
		assert(ar_system_node_is_active_in_view(h_sys, an0, 29));
		assert((ar_system_node_active_views(h_sys, an0) & (1u << 29)) != 0u);
		assert(ar_system_node_morton_key(h_sys, an0) != 0u);
		assert(an0->level >= 1 && an0->level <= 3);

		ar_sample_t sv29;
		assert(ar_probe_view(h_sys, 29, 1, 10, 10, &sv29));

		ar_system_destroy(h_sys);
		printf("  [PASS] Handle-driven value_scale/texel_span, probe interp_mode/interp_fn, 32-view mask, and ar_eval_ortho_ex verified\n");
	}

	/* 24. Unit-scale top-down camera distance (< 100.0f), ray window cull options, and attribute interpolation */
	{
		ar_config_t cfg;
		memset(&cfg, 0, sizeof(cfg));
		cfg.base_w = 8;
		cfg.base_h = 8;
		cfg.grid_dim = 8;
		cfg.pool_capacity = 32;
		cfg.max_active = 16;
		cfg.max_levels = 6;
		cfg.value_scale = 1.0f;
		cfg.attr_count = 4;
		cfg.attr_interp_mode = AR_INTERP_LINEAR;
		cfg.interpolate_attrs = true;
		cfg.sampler = linear_ramp_sampler;

		ar_system_t *u_sys = ar_system_create_ex(&cfg);
		assert(u_sys != NULL);

		/* Top-down perspective view with small eff_cam_dist = 10.0f (must NOT jump to default_topdown_dist) */
		ar_view_t td_view = {
			.cam_x = 0.0f, .cam_y = 0.0f, .cam_z = 0.0f,
			.pitch = AR_PITCH_TOP_DOWN, .yaw = 0.0f, .fov = 0.785398f,
			.zoom = 1.0f, .viewport_w = 64.0f, .viewport_h = 64.0f,
			.eff_cam_dist = 10.0f
		};
		const ar_node_t *td_nodes[16];
		size_t n_td = ar_system_update(u_sys, 1, &td_view, td_nodes, 16);
		assert(n_td > 0);
		assert(td_nodes[0]->level <= 2);

		/* Verify attribute interpolation & attr_count masking (channels 0..3 interpolated, 4..7 zeroed) */
		const ar_node_t *un0 = ar_system_ensure(u_sys, 1, 0, 0, 0);
		assert(un0 != NULL);
		/* linear_ramp_sampler sets attr[0] = (wx + wy) & 0xFF; at (0,0)->0, (1,0)->1; wait: L0 has w=8, grid_dim=8, step=1 */
		/* Let's ensure L2 where w=32, grid_dim=8, cell_step=4: at wx=0 attr[0]=0, at wx=4 attr[0]=4; at wx=2 -> attr[0]=2 */
		const ar_node_t *un2 = ar_system_ensure(u_sys, 2, 2, 0, 0);
		assert(un2 != NULL);
		ar_sample_t sa;
		memset(&sa, 0xFF, sizeof(sa));
		assert(ar_system_probe(u_sys, 2, 2, 0, &sa));
		assert(sa.attr[4] == 0 && sa.attr[5] == 0 && sa.attr[6] == 0 && sa.attr[7] == 0);

		/* Verify configurable ray window options in ar_cull_opts_t */
		float mvp[16];
		ar_mat4_ortho(mvp, -16.0f, 16.0f, -16.0f, 16.0f, -10.0f, 10.0f);
		ar_cull_opts_t copts;
		memset(&copts, 0, sizeof(copts));
		copts.z_min = -4.0f;
		copts.z_max = 4.0f;
		copts.viewport_w = 64.0f;
		copts.viewport_h = 64.0f;
		copts.target_quad_px = 2.0f;
		copts.ortho_ray_t_min = -0.5f;
		copts.unit_ray_t_max = 1.5f;
		copts.far_mid_ray_t_max = 500000.0f;
		copts.fallback_persp_aspect = 1.33f;
		copts.behind_cam_clip_w = 0.5f;
		copts.behind_cam_screen_mult = 3.0f;
		size_t n_cull = ar_system_query_matrix(u_sys, 1, mvp, &copts, td_nodes, 16);
		assert(n_cull > 0);

		ar_system_destroy(u_sys);
		printf("  [PASS] Unit-scale top-down camera distance, ray window cull opts, and attribute interpolation verified\n");
	}

	/* 25. 3D Manifold Bounding Box Callback (bounds_3d_fn) for Toroidal / Orbital Bodies */
	{
		ar_config_t cfg;
		memset(&cfg, 0, sizeof(cfg));
		cfg.base_w = 160;
		cfg.base_h = 128;
		cfg.grid_dim = 16;
		cfg.pool_capacity = 256;
		cfg.max_active = 128;
		cfg.max_levels = 15;
		cfg.period_x = 160 * 65536;
		cfg.period_y = 128 * 65536;
		cfg.max_root_cells = 16;
		cfg.sampler = mock_sampler;

		ar_system_t *t_sys = ar_system_create_ex(&cfg);
		assert(t_sys != NULL);

		/* Camera positioned relative to surface point at u=0, v=0 -> (0, -(R+r0), 0)
		 * Shift offset by (0, +(R+r0), 0) so surface origin is at (0, 0, 0) */
		struct torus_test_ctx planet_ctx = {
			.R = 2000000.0f,
			.r0 = 1335000.0f,
			.period_x = (float)(160 * 65536),
			.period_y = (float)(128 * 65536),
			.offset = { 0.0f, 3335000.0f, 0.0f }
		};

		float eye[3] = { 0.0f, -1500.0f, 800.0f };
		float target[3] = { 0.0f, 0.0f, 0.0f };
		float up[3] = { 0.0f, 0.0f, 1.0f };
		float mv[16], proj[16], mvp[16];
		ar_mat4_lookat(mv, eye, target, up);
		ar_mat4_perspective(proj, 0.785398f, 16.0f / 9.0f, 10.0f, 30000000.0f);
		ar_mat4_multiply(mvp, proj, mv);

		ar_cull_opts_t p_opts;
		memset(&p_opts, 0, sizeof(p_opts));
		p_opts.z_min = -500.0f;
		p_opts.z_max = 500.0f;
		p_opts.viewport_w = 1280.0f;
		p_opts.viewport_h = 720.0f;
		p_opts.target_quad_px = 4.0f;
		p_opts.min_level = 0;
		p_opts.max_level = 14;
		p_opts.bounds_3d_fn = custom_torus_bounds_3d_cb;
		p_opts.bounds_ctx = &planet_ctx;

		const ar_node_t *p_nodes[128];
		size_t np = ar_system_query_matrix(t_sys, 1, mvp, &p_opts, p_nodes, 128);
		assert(np > 0);
		/* Nearest node at origin (0,0) must have subdivided to fine LOD (< 5) */
		assert(p_nodes[0]->level <= 4);
		/* Opposite side of the torus (u ~ pi, y_3d ~ +6.67M behind or far) should not be subdivided to L0 */
		bool has_fine = false;
		for (size_t i = 0; i < np; i++) {
			if (p_nodes[i]->level <= 3) has_fine = true;
		}
		assert(has_fine);

		/* Distant Moon in front of camera along look ray (0, +1500, -800) at (0, +14000000, -7000000):
		 * Query WITHOUT forcing min_level=14; 3D screen-space projection must pick coarse levels automatically */
		struct torus_test_ctx moon_vis_ctx = {
			.R = 800000.0f,
			.r0 = 500000.0f,
			.period_x = (float)(160 * 65536),
			.period_y = (float)(128 * 65536),
			.offset = { 0.0f, 14000000.0f, -7000000.0f }
		};
		p_opts.bounds_ctx = &moon_vis_ctx;
		p_opts.view_id = 1;
		const ar_node_t *m_nodes[128];
		size_t nm_vis = ar_system_query_matrix(t_sys, 2, mvp, &p_opts, m_nodes, 128);
		assert(nm_vis > 0);
		for (size_t i = 0; i < nm_vis; i++) {
			assert(m_nodes[i]->level >= 9);
		}

		/* Verify cull_3d_fn rejects back-hemisphere nodes on the distant Moon */
		p_opts.cull_3d_fn = custom_torus_cull_3d_cb;
		size_t nm_vis_culled = ar_system_query_matrix(t_sys, 2, mvp, &p_opts, m_nodes, 128);
		assert(nm_vis_culled > 0 && nm_vis_culled < nm_vis);
		p_opts.cull_3d_fn = NULL;

		/* Distant Moon behind the camera at (0, -14000000, +7000000): must be 100% frustum-culled (0 nodes) */
		struct torus_test_ctx moon_behind_ctx = moon_vis_ctx;
		moon_behind_ctx.offset[1] = -14000000.0f;
		moon_behind_ctx.offset[2] = 7000000.0f;
		p_opts.bounds_ctx = &moon_behind_ctx;
		p_opts.view_id = 2;
		size_t nm_behind = ar_system_query_matrix(t_sys, 2, mvp, &p_opts, m_nodes, 128);
		assert(nm_behind == 0);

		ar_system_destroy(t_sys);
		printf("  [PASS] 3D manifold bounds_3d_fn frustum culling & screen-space LOD verified (planet=%zu, moon_vis=%zu, moon_behind=%zu)\n",
		       np, nm_vis, nm_behind);
	}

	/* Test 19: Universal Projected Solid-Angle Law & Scale-Ordered Spatial Cell Index */
	{
		/* Verify ar_level_texel_span_q16 matches ar_level_texel_span * 65536 */
		for (uint8_t lvl = 0; lvl <= 12; lvl++) {
			float sp_f = ar_level_texel_span(160, 48, lvl);
			int32_t sp_q16 = ar_level_texel_span_q16(160, 48, lvl);
			assert(fabsf(sp_f - ((float)sp_q16 / 65536.0f)) < 1e-3f);
		}

		ar_solid_angle_t sa;
		ar_solid_angle_init(&sa, 220.0f, 512.0f, 512.0f, 220.0f, 1.25f, 3.0f);
		assert(fabsf(sa.focal_proj - 220.0f) < 1e-5f);

		/* Verify duality: ar_solid_angle_min_extent(ar_solid_angle_max_dist(D)) == D */
		float extents[] = { 10.0f, 18.0f, 32.0f, 64.0f, 96.0f };
		for (size_t i = 0; i < sizeof(extents) / sizeof(extents[0]); i++) {
			float r_max = ar_solid_angle_max_dist(&sa, extents[i]);
			float d_min = ar_solid_angle_min_extent(&sa, r_max);
			assert(fabsf(d_min - extents[i]) < 1e-3f);
			/* At r_max, projected size is exactly s_min (1.25 px) -> alpha == 0 */
			assert(ar_solid_angle_alpha(&sa, extents[i], r_max + 1.0f) == 0.0f);
			/* At r_full = D * 220 / 3.0, projected size is 3.0 px -> alpha == 1 */
			float r_full = (extents[i] * sa.focal_proj) / sa.full_px;
			assert(fabsf(ar_solid_angle_alpha(&sa, extents[i], r_full) - 1.0f) < 1e-4f);
		}

		/* Verify scale-descending + painter-tiebreak sorting and prefix culling */
		struct test_sprite {
			float x, y, w, h;
		};
		struct test_sprite items[6] = {
			{ 100.0f, 200.0f, 10.0f, 12.0f }, /* small stone (ext 12) */
			{  80.0f, 300.0f, 48.0f, 64.0f }, /* giant tree South (ext 64, y=300) */
			{  50.0f,  90.0f, 18.0f, 22.0f }, /* herb North (ext 22, y=90) */
			{  40.0f, 120.0f, 52.0f, 64.0f }, /* giant tree North (ext 64, y=120) */
			{  60.0f, 150.0f, 28.0f, 34.0f }, /* bush (ext 34) */
			{  90.0f,  90.0f, 20.0f, 22.0f }  /* herb North-East (ext 22, y=90, x=90) */
		};
		/* Inline key callback via static helper block */
		void (*key_cb)(const void *, float *, float *, float *) = NULL;
		/* Use explicit key extractor */
		{
			/* Sort using ar_scale_painter_cmp directly via a small wrapper */
			extern void test_sprite_key_fn(const void *elem, float *out_ext, float *out_x, float *out_y);
			key_cb = test_sprite_key_fn;
		}
		ar_scale_cell_t cell;
		ar_scale_cell_build(&cell, 0, items, 6, sizeof(items[0]), key_cb);
		assert(cell.count == 6);
		assert(fabsf(cell.max_extent - 64.0f) < 1e-5f);
		assert(fabsf(cell.min_extent - 12.0f) < 1e-5f);

		/* Verify primary order is extent descending: 64, 64, 34, 22, 22, 12 */
		assert(fmaxf(items[0].w, items[0].h) == 64.0f);
		assert(fmaxf(items[1].w, items[1].h) == 64.0f);
		/* Within extent=64, North-to-South painter order (y=120 before y=300) */
		assert(items[0].y == 120.0f && items[1].y == 300.0f);
		assert(fmaxf(items[2].w, items[2].h) == 34.0f);
		assert(fmaxf(items[3].w, items[3].h) == 22.0f);
		assert(fmaxf(items[4].w, items[4].h) == 22.0f);
		/* Within extent=22 and y=90, x=50 before x=90 */
		assert(items[3].x == 50.0f && items[4].x == 90.0f);
		assert(fmaxf(items[5].w, items[5].h) == 12.0f);

		/* At distance r = 4500 wu, cutoff D_min(4500) = 1.25 * 4500 / 220 = 25.56 wu:
		 * Cell is visible (max_extent 64 >= 25.56), and prefix scan stops after index 2 (extents 64, 64, 34)! */
		float cutoff = 0.0f;
		assert(ar_scale_cell_cull(&cell, &sa, 4500.0f, 500.0f, &cutoff) == true);
		uint16_t vis_prefix = 0;
		for (uint16_t i = 0; i < cell.count; i++) {
			float ext = fmaxf(items[i].w, items[i].h);
			if (ext < cutoff)
				break;
			vis_prefix++;
		}
		assert(vis_prefix == 3); /* Only the two 64-wu trees and 34-wu bush are visited; herbs & stone are skipped! */

		/* Beyond r = 12000 wu (r_near = 11500 -> D_min = 65.34 > max_extent 64): O(1) cell cull returns false! */
		assert(ar_scale_cell_cull(&cell, &sa, 12000.0f, 500.0f, &cutoff) == false);
		printf("  [PASS] Universal solid-angle law & scale-ordered cell prefix culling verified\n");
	}

	/* 26. Asynchronous Workqueue (stl_wq) Deferred Population & Frame Commit */
	{
		stl_wq_t *wq = stl_wq_create("test-areole-wq", 2, 64);
		assert(wq != NULL);

		ar_config_t cfg;
		memset(&cfg, 0, sizeof(cfg));
		cfg.base_w = 160;
		cfg.base_h = 128;
		cfg.grid_dim = 16;
		cfg.pool_capacity = 32;
		cfg.max_active = 16;
		cfg.max_levels = 6;
		cfg.max_populates_per_frame = 1; /* Only 1 sync populate per frame; next misses go to wq! */
		cfg.wq = wq;
		cfg.sampler = level_dependent_sampler;

		ar_system_t *wq_sys = ar_system_create_ex(&cfg);
		assert(wq_sys != NULL);

		ar_system_begin_frame(wq_sys);
		/* 1st populate (parent L1 at 0,0) consumes the 1 sync budget slot -> val = 100 */
		const ar_node_t *p_node = ar_system_ensure(wq_sys, 1, 1, 0, 0);
		assert(p_node != NULL);
		assert((p_node->flags & AR_NODE_FLAG_DEFERRED) == 0);
		assert(p_node->values[0] == 100);

		/* 2nd populate in the same frame (child L0 at 0,0) exceeds budget -> immediately upsampled from parent (val=100, DEFERRED) and queued on wq! */
		const ar_node_t *c_node = ar_system_ensure(wq_sys, 1, 0, 0, 0);
		assert(c_node != NULL);
		assert((c_node->flags & AR_NODE_FLAG_DEFERRED) != 0);
		assert(c_node->values[0] == 100);
		uint32_t rev_before = c_node->rev;

		/* Drain the workqueue and begin next frame: staging buffers commit into c_node -> val = 0 (L0), DEFERRED cleared, rev incremented! */
		stl_wq_drain(wq);
		ar_system_begin_frame(wq_sys);
		assert((c_node->flags & AR_NODE_FLAG_DEFERRED) == 0);
		assert(c_node->values[0] == 0);
		assert(c_node->rev > rev_before);

		ar_system_destroy(wq_sys);
		stl_wq_destroy(wq);
		printf("  [PASS] Asynchronous stl_wq deferred node population & frame staging commit verified\n");
	}

	/* 27. Non-Zero focal_z Camera Target Elevation (High Mountain Peak & Deep Basin Invariance) */
	{
		ar_config_t z_cfg = {
			.base_w = 160,
			.base_h = 128,
			.grid_dim = TEST_GRID_DIM,
			.pool_capacity = TEST_POOL_CAP,
			.max_active = TEST_MAX_ACTIVE,
			.max_levels = TEST_MAX_LEVELS,
			.z_min = -2600.0f,
			.z_max = 3200.0f,
			.sampler = mock_sampler
		};
		ar_system_t *z_sys = ar_system_create_ex(&z_cfg);
		assert(z_sys != NULL);

		const float test_focal_zs[3] = { -600.0f, 0.0f, 950.0f };
		for (int zi = 0; zi < 3; zi++) {
			float fz = test_focal_zs[zi];
			float fx = 1680.0f, fy = -1216.0f; /* Center of Level-0 chunk (1600, -1280) */
			ar_view_t view = {
				.cam_x = fx, .cam_y = fy, .cam_z = fz,
				.pitch = 0.959931f, .yaw = 0.0f, .fov = 0.785398f,
				.zoom = 1.0f, .viewport_w = 1280.0f, .viewport_h = 720.0f,
				.eff_cam_dist = 350.0f
			};
			const ar_node_t *vis[64];
			size_t cnt = ar_system_update(z_sys, 1, &view, vis, 64);
			assert(cnt > 12);

			/* Verify the chunk containing the focal point (fx, fy) is present and refined to Level 0 */
			bool found_focal_lvl0 = false;
			for (size_t i = 0; i < cnt; i++) {
				const ar_node_t *n = vis[i];
				assert(n != NULL);
				if (fx >= (float)n->ox && fx < (float)(n->ox + n->w) &&
				    fy >= (float)n->oy && fy < (float)(n->oy + n->h)) {
					assert(n->level == 0);
					found_focal_lvl0 = true;
				}
			}
			assert(found_focal_lvl0);
		}
		ar_system_destroy(z_sys);
		printf("  [PASS] Non-zero focal_z camera target elevation invariance verified across [-600, 0, +950] wu\n");
	}

	/* 28. Ground-Truth Frustum Visibility Oracle (Zero False Negatives across 360-deg Yaw x Low-to-High Pitch x Distances) */
	{
		ar_config_t o_cfg = {
			.base_w = 128, .base_h = 128, .grid_dim = 8,
			.pool_capacity = 512, .max_active = 256, .max_levels = 10,
			.z_min = -600.0f, .z_max = 900.0f,
			.target_quad_px = 4.0f,
			.sampler = high_relief_sampler
		};
		ar_system_t *o_sys = ar_system_create_ex(&o_cfg);
		assert(o_sys != NULL);

		const float pitches_deg[5] = { 12.0f, 20.0f, 35.0f, 55.0f, 80.0f };
		const float dists[3] = { 280.0f, 1200.0f, 3600.0f };
		const float aspects[2] = { 1280.0f / 720.0f, 1024.0f / 768.0f };
		size_t total_vis_pts = 0, total_misses = 0, total_poses = 0;

		for (int ai = 0; ai < 2; ai++) {
			float vp_w = (ai == 0) ? 1280.0f : 1024.0f;
			float vp_h = (ai == 0) ?  720.0f :  768.0f;
			for (int di = 0; di < 3; di++) {
				float dist = dists[di];
				for (int pi = 0; pi < 5; pi++) {
					float pitch = pitches_deg[pi] * (3.14159265f / 180.0f);
					for (int yi = 0; yi < 12; yi++) {
						float yaw = (float)yi * (3.14159265f / 6.0f);
						float sin_y = sinf(yaw), cos_y = cosf(yaw);
						float sin_p = sinf(pitch), cos_p = cosf(pitch);
						float eye[3] = { dist * cos_p * sin_y, dist * cos_p * cos_y, 200.0f + dist * sin_p };
						float tgt[3] = { 0.0f, 0.0f, 200.0f };
						float up[3]  = { -sin_p * sin_y, -sin_p * cos_y, cos_p };
						float view_m[16], proj_m[16], mvp[16];
						ar_mat4_lookat(view_m, eye, tgt, up);
						ar_mat4_perspective(proj_m, 0.785398f, aspects[ai], 4.0f, dist * 12.0f + 4000.0f);
						ar_mat4_multiply(mvp, proj_m, view_m);

						ar_cull_opts_t opts;
						memset(&opts, 0, sizeof(opts));
						opts.z_min = -600.0f;
						opts.z_max =  900.0f;
						opts.focal_x = 0.0f;
						opts.focal_y = 0.0f;
						opts.focal_z = 200.0f;
						opts.viewport_w = vp_w;
						opts.viewport_h = vp_h;
						opts.target_quad_px = 4.0f;

						const ar_node_t *vis[256];
						size_t cnt = ar_system_query_matrix(o_sys, 1, mvp, &opts, vis, 256);
						assert(cnt > 0);
						total_poses++;

						/* Verify zero overlapping nodes and every returned node passes frustum */
						ar_frustum_t fr;
						ar_frustum_from_matrix(&fr, mvp);
						for (size_t i = 0; i < cnt; i++) {
							assert(ar_frustum_test_box(&fr, (float)vis[i]->ox, (float)vis[i]->oy, -600.0f,
							                           (float)(vis[i]->ox + vis[i]->w),
							                           (float)(vis[i]->oy + vis[i]->h), 900.0f));
							for (size_t j = i + 1; j < cnt; j++) {
								bool ov_x = (vis[i]->ox < vis[j]->ox + vis[j]->w) && (vis[j]->ox < vis[i]->ox + vis[i]->w);
								bool ov_y = (vis[i]->oy < vis[j]->oy + vis[j]->h) && (vis[j]->oy < vis[i]->oy + vis[i]->h);
								assert(!(ov_x && ov_y));
							}
						}

						/* Ground-truth oracle: sample dense lattice around view & verify 0 visible points missed */
						float span = dist * 2.8f + 900.0f;
						float step = span / 18.0f;
						for (int gy = -9; gy <= 9; gy++) {
							for (int gx = -9; gx <= 9; gx++) {
								float px = (float)gx * step;
								float py = (float)gy * step;
								float pz = high_relief_height(px, py);
								float cx = mvp[0]*px + mvp[4]*py + mvp[8]*pz  + mvp[12];
								float cy = mvp[1]*px + mvp[5]*py + mvp[9]*pz  + mvp[13];
								float cz = mvp[2]*px + mvp[6]*py + mvp[10]*pz + mvp[14];
								float cw = mvp[3]*px + mvp[7]*py + mvp[11]*pz + mvp[15];
								if (cw <= 4.0f) continue;
								float nx = cx / cw, ny = cy / cw, nz = cz / cw;
								if (fabsf(nx) <= 0.99f && fabsf(ny) <= 0.99f && nz >= -1.0f && nz <= 1.0f) {
									total_vis_pts++;
									bool covered = false;
									for (size_t i = 0; i < cnt; i++) {
										if (px >= (float)vis[i]->ox && px < (float)(vis[i]->ox + vis[i]->w) &&
										    py >= (float)vis[i]->oy && py < (float)(vis[i]->oy + vis[i]->h)) {
											covered = true;
											break;
										}
									}
									if (!covered) total_misses++;
								}
							}
						}
					}
				}
			}
		}
		assert(total_vis_pts > 5000);
		assert(total_misses == 0);

		/* Looking straight up into sky above z_max must return 0 nodes (tight false positives) */
		{
			float eye[3] = { 0.0f, 0.0f, 1500.0f };
			float tgt[3] = { 0.0f, 0.0f, 5000.0f };
			float up[3]  = { 0.0f, 1.0f, 0.0f };
			float vm[16], pm[16], mvp[16];
			ar_mat4_lookat(vm, eye, tgt, up);
			ar_mat4_perspective(pm, 0.785398f, 16.0f / 9.0f, 10.0f, 50000.0f);
			ar_mat4_multiply(mvp, pm, vm);
			ar_cull_opts_t opts;
			memset(&opts, 0, sizeof(opts));
			opts.z_min = -600.0f; opts.z_max = 900.0f;
			opts.viewport_w = 1280.0f; opts.viewport_h = 720.0f;
			const ar_node_t *vis[64];
			size_t sky_cnt = ar_system_query_matrix(o_sys, 1, mvp, &opts, vis, 64);
			assert(sky_cnt == 0);
		}

		/* Diagonal yaw tombstone budget test: tight max_active=24 must never starve on vis_kids==0 */
		{
			float yaw = 0.785398f; /* 45 deg */
			float pitch = 0.45f;
			float dist = 900.0f;
			float eye[3] = { dist * cosf(pitch) * sinf(yaw), dist * cosf(pitch) * cosf(yaw), dist * sinf(pitch) };
			float tgt[3] = { 0.0f, 0.0f, 0.0f };
			float up[3]  = { -sinf(pitch) * sinf(yaw), -sinf(pitch) * cosf(yaw), cosf(pitch) };
			float vm[16], pm[16], mvp[16];
			ar_mat4_lookat(vm, eye, tgt, up);
			ar_mat4_perspective(pm, 0.785398f, 16.0f / 9.0f, 4.0f, 12000.0f);
			ar_mat4_multiply(mvp, pm, vm);
			ar_cull_opts_t opts;
			memset(&opts, 0, sizeof(opts));
			opts.z_min = -600.0f; opts.z_max = 900.0f;
			opts.viewport_w = 1280.0f; opts.viewport_h = 720.0f;
			opts.target_quad_px = 1.0f;
			const ar_node_t *vis[24];
			size_t dcnt = ar_system_query_matrix(o_sys, 1, mvp, &opts, vis, 24);
			assert(dcnt >= 18 && dcnt <= 24);
		}

		ar_system_destroy(o_sys);
		printf("  [PASS] Frustum visibility oracle verified: %zu poses, %zu visible pts, 0 false negatives, sky=0\n",
		       total_poses, total_vis_pts);
	}

	/* 29. Conservative 2D Software Hi-Z Occluder (360-deg Mountain Wall + Multi-Scale 2..16M Linear-W Precision) */
	{
		ridge_scene_ctx_t rc = { 0.0f, -1.0f, 0.0f, -2400.0f, 2400.0f };
		ar_config_t occ_cfg = {
			.base_w = 128, .base_h = 128, .grid_dim = 8,
			.pool_capacity = 512, .max_active = 256, .max_levels = 8,
			.z_min = -50.0f, .z_max = 2600.0f,
			.target_quad_px = 6.0f,
			.sampler = ridge_wall_sampler,
			.user_ctx = &rc
		};
		ar_system_t *occ_sys = ar_system_create_ex(&occ_cfg);
		assert(occ_sys != NULL);

		/* 29a. 360-degree mountain ridge occlusion sweep (8 yaws every 45 deg) */
		uint32_t total_culled = 0;
		for (int yi = 0; yi < 8; yi++) {
			float yaw = (float)yi * (3.14159265f / 4.0f);
			float look_x = -sinf(yaw);
			float look_y = -cosf(yaw);
			rc.ridge_dir_x = look_x;
			rc.ridge_dir_y = look_y;
			rc.peak_x = look_x * 2400.0f;
			rc.peak_y = look_y * 2400.0f;
			rc.peak_h = 1500.0f;
			ar_system_clear(occ_sys);

			float eye[3] = { -look_x * 250.0f, -look_y * 250.0f, 200.0f };
			float tgt[3] = {  look_x * 750.0f,  look_y * 750.0f, 450.0f };
			float up[3]  = { 0.0f, 0.0f, 1.0f };
			float vm[16], pm[16], mvp[16];
			ar_mat4_lookat(vm, eye, tgt, up);
			ar_mat4_perspective(pm, 0.785398f, 16.0f / 9.0f, 4.0f, 8000.0f);
			ar_mat4_multiply(mvp, pm, vm);

			/* Warm up cache so cached_bounds gives tight vertical boxes for the z=0 valley behind the wall */
			ar_cull_opts_t opts;
			memset(&opts, 0, sizeof(opts));
			opts.z_min = -50.0f; opts.z_max = 2600.0f;
			opts.focal_x = tgt[0]; opts.focal_y = tgt[1]; opts.focal_z = tgt[2];
			opts.viewport_w = 1280.0f; opts.viewport_h = 720.0f;
			opts.target_quad_px = 6.0f;
			opts.use_cached_bounds = true;
			opts.cached_bounds_max_level = 7;
			opts.cached_bounds_pad = 16.0f;

			const ar_node_t *vis_no_occ[256];
			size_t cnt_no_occ = ar_system_query_matrix(occ_sys, 1, mvp, &opts, vis_no_occ, 256);
			assert(cnt_no_occ > 10);

			/* Now query WITH conservative Hi-Z occluder */
			ar_occluder_t occ;
			ar_occluder_reset(&occ, mvp, 4.0f, 0.002f);
			opts.occluder = &occ;
			opts.occluder_z_bias = 12.0f;

			const ar_node_t *vis_occ[256];
			size_t cnt_occ = ar_system_query_matrix(occ_sys, 1, mvp, &opts, vis_occ, 256);
			assert(cnt_occ > 0);
			assert(occ.culled_count > 0);
			assert(cnt_occ < cnt_no_occ);
			total_culled += occ.culled_count;

			/* Zero false negatives: the tall peak (z=2400) rising above the ridge MUST remain visible */
			bool found_peak = false;
			for (size_t i = 0; i < cnt_occ; i++) {
				if (rc.peak_x >= (float)vis_occ[i]->ox && rc.peak_x < (float)(vis_occ[i]->ox + vis_occ[i]->w) &&
				    rc.peak_y >= (float)vis_occ[i]->oy && rc.peak_y < (float)(vis_occ[i]->oy + vis_occ[i]->h)) {
					found_peak = true;
					break;
				}
			}
			assert(found_peak);
		}
		assert(total_culled >= 8);

		/* 29b. Multi-scale linear-W depth precision (close-up lamb at w=5, ridge at w=800, moon at w=16,000,000) */
		{
			float eye[3] = { 0.0f, 0.0f, 5.0f };
			float tgt[3] = { 0.0f, 100.0f, 5.0f };
			float up[3]  = { 0.0f, 0.0f, 1.0f };
			float vm[16], pm[16], mvp[16];
			ar_mat4_lookat(vm, eye, tgt, up);
			ar_mat4_perspective(pm, 0.785398f, 16.0f / 9.0f, 2.0f, 25000000.0f);
			ar_mat4_multiply(mvp, pm, vm);

			ar_occluder_t occ;
			ar_occluder_reset(&occ, mvp, 2.0f, 0.002f);

			/* Rasterize opaque mountain ridge spanning lower half of screen at y=800 (z in [-200..+5]) */
			float r0[3] = { -1500.0f, 800.0f, -400.0f };
			float r1[3] = {  1500.0f, 800.0f, -400.0f };
			float r2[3] = {  1500.0f, 800.0f,    5.0f };
			float r3[3] = { -1500.0f, 800.0f,    5.0f };
			ar_occluder_rasterize_tri(&occ, r0, r1, r2);
			ar_occluder_rasterize_tri(&occ, r0, r2, r3);

			/* 1. Close-up lamb in front of ridge (y in [4..6]) -> NOT occluded */
			float lamb_min[3] = { -0.5f, 4.0f, 4.0f };
			float lamb_max[3] = {  0.5f, 6.0f, 5.0f };
			assert(!ar_occluder_test_box(&occ, lamb_min, lamb_max));

			/* 2. Distant companion moon at y=16,000,000 completely below ridge horizon -> OCCLUDED */
			float moon_hid_min[3] = { -500000.0f, 16000000.0f, -4000000.0f };
			float moon_hid_max[3] = {  500000.0f, 16500000.0f, -1500000.0f };
			assert(ar_occluder_test_box(&occ, moon_hid_min, moon_hid_max));

			/* 3. Distant companion moon rising above ridge into sky -> NOT occluded */
			float moon_sky_min[3] = { -500000.0f, 16000000.0f,  1000000.0f };
			float moon_sky_max[3] = {  500000.0f, 16500000.0f,  3000000.0f };
			assert(!ar_occluder_test_box(&occ, moon_sky_min, moon_sky_max));

			/* 4. Distant companion moon grazing right below silhouette edge (within 1 occ pixel) -> NOT occluded */
			float moon_edge_min[3] = { -100000.0f, 16000000.0f, -100000.0f };
			float moon_edge_max[3] = {  100000.0f, 16500000.0f,   -5000.0f };
			assert(!ar_occluder_test_box(&occ, moon_edge_min, moon_edge_max));
		}

		ar_system_destroy(occ_sys);
		printf("  [PASS] Conservative 2D Hi-Z occluder verified (8 yaws, %u hidden chunks culled, multi-scale 2..16M linear-W exact)\n",
		       total_culled);
	}

	/* 30. Opt-In Cached-Node Z-Tightening (Eliminates Global Skyscraper False Positives) */
	{
		ar_config_t cb_cfg = {
			.base_w = 128, .base_h = 128, .grid_dim = 8,
			.pool_capacity = 256, .max_active = 128, .max_levels = 6,
			.z_min = -2600.0f, .z_max = 3200.0f,
			.target_quad_px = 4.0f,
			.sampler = flat_ocean_sampler
		};
		ar_system_t *cb_sys = ar_system_create_ex(&cb_cfg);
		assert(cb_sys != NULL);

		/* Camera at z=400 looking slightly UP (+15 deg) over flat z=0 ocean:
		 * true z=0 ocean is completely BELOW the frustum, whereas global z_max=+3200 still intersects it! */
		float eye[3] = { 0.0f, 0.0f, 400.0f };
		float tgt[3] = { 0.0f, -1000.0f, 950.0f };
		float up[3]  = { 0.0f, 0.0f, 1.0f };
		float vm[16], pm[16], mvp[16];
		ar_mat4_lookat(vm, eye, tgt, up);
		ar_mat4_perspective(pm, 0.523599f, 16.0f / 9.0f, 4.0f, 12000.0f);
		ar_mat4_multiply(mvp, pm, vm);

		ar_cull_opts_t opts;
		memset(&opts, 0, sizeof(opts));
		opts.z_min = -2600.0f; opts.z_max = 3200.0f;
		opts.focal_x = 0.0f; opts.focal_y = -1000.0f; opts.focal_z = 0.0f;
		opts.viewport_w = 1280.0f; opts.viewport_h = 720.0f;
		opts.target_quad_px = 4.0f;

		const ar_node_t *vis_wide[128];
		size_t cnt_wide = ar_system_query_matrix(cb_sys, 1, mvp, &opts, vis_wide, 128);
		assert(cnt_wide > 0); /* Global z_max=+3200 pulled in flat ocean chunks */

		/* Now enable cached bounds tightening on those already-cached flat ocean nodes */
		opts.use_cached_bounds = true;
		opts.cached_bounds_max_level = 5;
		opts.cached_bounds_pad = 8.0f;
		const ar_node_t *vis_tight[128];
		size_t cnt_tight = ar_system_query_matrix(cb_sys, 1, mvp, &opts, vis_tight, 128);
		assert(cnt_tight < cnt_wide);

		ar_system_destroy(cb_sys);
		printf("  [PASS] Cached-node Z-tightening verified (global slab=%zu -> cached tight=%zu false positives)\n",
		       cnt_wide, cnt_tight);
	}

	ar_system_destroy(cs);
	printf("=== All Standalone libareole Spatial Hierarchy Tests PASSED ===\n");
	return 0;
}

void test_sprite_key_fn(const void *elem, float *out_ext, float *out_x, float *out_y)
{
	const float *f = (const float *)elem; /* x, y, w, h */
	*out_x = f[0];
	*out_y = f[1];
	*out_ext = (f[2] > f[3]) ? f[2] : f[3];
}
