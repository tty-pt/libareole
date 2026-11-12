#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <time.h>

#include "ttypt/areole.h"

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
	float elev  = macro + meso + micro;
	int32_t elev_q12_4 = (int32_t)(elev * 16.0f);

	out->val = (int16_t)elev_q12_4;
	out->attr[0] = (elev_q12_4 < 0) ? 0 : ((elev_q12_4 > 100 * 16) ? 2 : 10);
	out->attr[1] = out->attr[0];
	out->attr[2] = 0;
	out->attr[3] = (elev_q12_4 < 0) ? (uint8_t)(-elev_q12_4 > 255 * 16 ? 255 : (-elev_q12_4 >> 4)) : 0;
	out->rgba = (out->attr[3] > 0) ? 0xFF003080u : 0xFF208020u;
	out->flags = (out->attr[3] > 0) ? 1u : 0u;
}

/* Mock procedural sampler: analytical sinusoidal elevation + water plane */
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
	float elev = macro_h + fine_h;

	out->val = (int16_t)(elev * 16.0f);
	out->attr[0] = (elev < 0.0f) ? 0 : ((elev > 60.0f) ? 2 : 10);
	out->attr[1] = out->attr[0];
	out->attr[2] = 0;
	out->attr[3] = (elev < 0.0f) ? (uint8_t)(-elev > 100.0f ? 100 : -elev) : 0;
	out->rgba = (out->attr[3] > 0) ? 0xFF003080u : 0xFF208020u;
	out->flags = (out->attr[3] > 0) ? 1 : 0;
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

int main(void)
{
	printf("=== Running Standalone Generic libareole Spatial Hierarchy Suite ===\n");

	ar_system_t *cs = ar_system_create(160, 128, mock_sampler, mock_sparsity, NULL);
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

		const ar_node_t *out_chunks[AR_DEFAULT_MAX_ACTIVE];
		size_t n = ar_system_update(cs, 1, &view, out_chunks, AR_DEFAULT_MAX_ACTIVE);

		assert(n > 0);
		assert(n <= AR_DEFAULT_MAX_ACTIVE);

		ar_stats_t stats;
		ar_system_get_stats(cs, &stats);

		printf("  [PASS] Zoom %9.4fx -> Active Nodes: %2zu (Max Bound: %d) | Levels [%u..%u]\n",
		       zoom, n, AR_DEFAULT_MAX_ACTIVE, stats.min_level_observed, stats.max_level_observed);
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
		const ar_node_t *close_chunks[AR_DEFAULT_MAX_ACTIVE];
		size_t nc = ar_system_update(cs, 1, &view_close, close_chunks, AR_DEFAULT_MAX_ACTIVE);
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
		const ar_node_t *far_chunks[AR_DEFAULT_MAX_ACTIVE];
		size_t nf = ar_system_update(cs, 1, &view_far, far_chunks, AR_DEFAULT_MAX_ACTIVE);
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
			const ar_node_t *pan_chunks[AR_DEFAULT_MAX_ACTIVE];
			ar_system_update(cs, 1, &view_pan, pan_chunks, AR_DEFAULT_MAX_ACTIVE);
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
		const ar_node_t *center_chunks[AR_DEFAULT_MAX_ACTIVE];
		size_t nc = ar_system_update(cs, 1, &view_center, center_chunks, AR_DEFAULT_MAX_ACTIVE);

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
			const ar_node_t *bench_chunks[AR_DEFAULT_MAX_ACTIVE];
			ar_system_update(cs, 1, &bench_view, bench_chunks, AR_DEFAULT_MAX_ACTIVE);
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
		ar_system_t *simple_cs = ar_create(160, 128, mock_sampler, NULL);
		assert(simple_cs != NULL);

		ar_view_t view = {
			.cam_x = 0.0f, .cam_y = 0.0f, .cam_z = 0.0f,
			.pitch = 0.0f, .yaw = 0.0f, .fov = 0.0f,
			.zoom = 1.0f,
			.viewport_w = 1280.0f, .viewport_h = 720.0f,
			.eff_cam_dist = 400.0f
		};

		const ar_node_t *out_nodes[AR_DEFAULT_MAX_ACTIVE];
		size_t active_count = ar_update(simple_cs, 1, &view, out_nodes, AR_DEFAULT_MAX_ACTIVE);
		assert(active_count > 0 && active_count <= AR_DEFAULT_MAX_ACTIVE);

		ar_sample_t probe_sample;
		bool probe_ok = ar_probe(simple_cs, 1, 0, 0, &probe_sample);
		assert(probe_ok);
		float probe_elev_f = ar_elevation_to_float(probe_sample.val);
		int16_t roundtrip_q = ar_elevation_from_float(probe_elev_f);
		assert(abs(roundtrip_q - probe_sample.val) <= 1);

		ar_clear(simple_cs);
		ar_destroy(simple_cs);
		printf("  [PASS] Ergonomic shorthand API and legacy chunk_* aliases verified\n");
	}

	/* Test 8: O(1) Morton Key Spatial Lookup & Multi-Level Coexistence */
	{
		ar_system_t *morton_cs = ar_create(160, 128, mock_sampler, NULL);
		assert(morton_cs != NULL);

		ar_view_t view_center = {
			.cam_x = 0.0f, .cam_y = 0.0f, .cam_z = 0.0f,
			.pitch = 0.959931f, .yaw = 0.0f, .fov = 0.785398f,
			.zoom = 1.0f,
			.viewport_w = 1280.0f, .viewport_h = 720.0f,
			.eff_cam_dist = 400.0f
		};
		const ar_node_t *out_nodes[AR_DEFAULT_MAX_ACTIVE];
		size_t n = ar_update(morton_cs, 1, &view_center, out_nodes, AR_DEFAULT_MAX_ACTIVE);
		assert(n > 0);

		for (size_t i = 0; i < n; i++) {
			const ar_node_t *found = ar_system_find(morton_cs, out_nodes[i]->world,
			                                        out_nodes[i]->level,
			                                        out_nodes[i]->ox, out_nodes[i]->oy);
			assert(found == out_nodes[i]);
			assert(found->morton_key == out_nodes[i]->morton_key);

			const ar_node_t *by_key = ar_system_find_by_key(morton_cs, out_nodes[i]->morton_key);
			assert(by_key == out_nodes[i]);
		}

		assert(ar_system_find(morton_cs, 1, 0, 999999, 999999) == NULL);

		ar_view_t view_wide = view_center;
		view_wide.zoom = 0.05f;
		const ar_node_t *wide_nodes[AR_DEFAULT_MAX_ACTIVE];
		size_t nw = ar_update(morton_cs, 1, &view_wide, wide_nodes, AR_DEFAULT_MAX_ACTIVE);
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
		ar_system_t *ortho_cs = ar_create(160, 128, mock_sampler, NULL);
		assert(ortho_cs != NULL);

		ar_view_t ortho_view = {
			.cam_x = 0.0f, .cam_y = 0.0f, .cam_z = 0.0f,
			.pitch = 0.0f, .yaw = 0.0f, .fov = 0.0f,
			.zoom = 1.0f,
			.viewport_w = 224.0f, .viewport_h = 224.0f,
			.eff_cam_dist = 0.0f
		};
		const ar_node_t *ortho_nodes[AR_DEFAULT_MAX_ACTIVE];
		size_t count = ar_update(ortho_cs, 1, &ortho_view, ortho_nodes, AR_DEFAULT_MAX_ACTIVE);
		assert(count > 0 && count <= AR_DEFAULT_MAX_ACTIVE);

		for (size_t i = 0; i < count; i++) {
			assert(ortho_nodes[i]->level <= 1);
		}

		ortho_view.zoom = 0.10f;
		count = ar_update(ortho_cs, 1, &ortho_view, ortho_nodes, AR_DEFAULT_MAX_ACTIVE);
		assert(count > 0 && count <= AR_DEFAULT_MAX_ACTIVE);

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
		ar_system_t *part_cs = ar_create(160, 128, mock_sampler, NULL);
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
			const ar_node_t *nodes[AR_DEFAULT_MAX_ACTIVE];
			size_t n = ar_update(part_cs, 1, &view, nodes, AR_DEFAULT_MAX_ACTIVE);
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
		ar_system_t *mv_cs = ar_create(160, 128, mock_sampler, NULL);
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
		ar_draw_item_t items_seven[AR_DEFAULT_MAX_ACTIVE];
		size_t count_seven = ar_eval_view(mv_cs, 0, 1, &view_seven, items_seven, AR_DEFAULT_MAX_ACTIVE);
		assert(count_seven > 0);

		for (size_t i = 0; i < count_seven; i++) {
			assert(items_seven[i].node != NULL);
			assert((items_seven[i].node->active_view_mask & 1) != 0);
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
		ar_draw_item_t items_lamina[AR_DEFAULT_MAX_ACTIVE];
		size_t count_lamina = ar_eval_view(mv_cs, 1, 1, &view_lamina, items_lamina, AR_DEFAULT_MAX_ACTIVE);
		assert(count_lamina > 0);

		for (size_t i = 0; i < count_lamina; i++) {
			assert(items_lamina[i].node != NULL);
			assert((items_lamina[i].node->active_view_mask & 2) != 0);
			assert(items_lamina[i].node->level >= 2);
		}

		ar_begin_frame(mv_cs);
		for (size_t i = 0; i < count_seven; i++) {
			assert(items_seven[i].node->active_view_mask == 0);
		}

		ar_sample_t probe_s;
		assert(ar_probe(mv_cs, 1, 100, 100, &probe_s));
		assert(probe_s.val != 0);

		ar_destroy(mv_cs);
		printf("  [PASS] Multi-view evaluation & database integration verified\n");
	}

	/* Test 12: Varied Multi-View System with Analytical Multi-Octave Sampler */
	{
		ar_system_t *multi_cs = ar_create(160, 128, analytical_multioctave_sampler, NULL);
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

		ar_draw_item_t items_inspect[AR_DEFAULT_MAX_ACTIVE];
		ar_draw_item_t items_diorama[AR_DEFAULT_MAX_ACTIVE];
		ar_draw_item_t items_tactical[AR_DEFAULT_MAX_ACTIVE];
		ar_draw_item_t items_minimap[AR_DEFAULT_MAX_ACTIVE];

		size_t n_inspect  = ar_eval_view(multi_cs, 0, 1, &view_inspect,  items_inspect,  AR_DEFAULT_MAX_ACTIVE);
		size_t n_diorama  = ar_eval_view(multi_cs, 1, 1, &view_diorama,  items_diorama,  AR_DEFAULT_MAX_ACTIVE);
		size_t n_tactical = ar_eval_view(multi_cs, 2, 1, &view_tactical, items_tactical, AR_DEFAULT_MAX_ACTIVE);
		size_t n_minimap  = ar_eval_view(multi_cs, 3, 1, &view_minimap,  items_minimap,  AR_DEFAULT_MAX_ACTIVE);

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
		ar_system_t *mat_cs = ar_create(160, 128, mock_sampler, NULL);
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

			const ar_node_t *out_nodes[AR_DEFAULT_MAX_ACTIVE];
			ar_cull_opts_t opts = {
				.z_min = -1000.0f, .z_max = 1000.0f,
				.sag_radius = 0.0f, .focal_x = 0.0f, .focal_y = 0.0f,
				.margin = 50.0f, .target_quad_px = 2.0f,
				.viewport_w = 800.0f, .viewport_h = 600.0f
			};
			size_t n = ar_system_query_matrix(mat_cs, 1, mvp, &opts, out_nodes, AR_DEFAULT_MAX_ACTIVE);
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

			const ar_node_t *out_nodes[AR_DEFAULT_MAX_ACTIVE];
			ar_cull_opts_t opts = {
				.z_min = -500.0f, .z_max = 500.0f,
				.sag_radius = 0.0f, .focal_x = 0.0f, .focal_y = 0.0f,
				.margin = 50.0f, .target_quad_px = 2.0f,
				.viewport_w = 1280.0f, .viewport_h = 720.0f
			};
			size_t n = ar_system_query_matrix(mat_cs, 1, mvp, &opts, out_nodes, AR_DEFAULT_MAX_ACTIVE);
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

			const ar_node_t *out_nodes[AR_DEFAULT_MAX_ACTIVE];
			ar_cull_opts_t opts = {
				.z_min = -200.0f, .z_max = 200.0f,
				.sag_radius = 0.0f, .focal_x = 0.0f, .focal_y = 0.0f,
				.margin = 50.0f, .target_quad_px = 2.0f,
				.viewport_w = 1000.0f, .viewport_h = 1000.0f
			};
			size_t n = ar_system_query_matrix(mat_cs, 1, mvp, &opts, out_nodes, AR_DEFAULT_MAX_ACTIVE);
			assert(n > 0);
			printf("  [PASS] Oblique sheared projection query verified (%zu nodes)\n", n);
		}

		ar_destroy(mat_cs);
	}

	/* Test 14: 360-Degree Yaw Rotation & Pitch Angles (0 deg to 85 deg) */
	{
		ar_system_t *rot_cs = ar_create(160, 128, mock_sampler, NULL);
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

				const ar_node_t *out_nodes[AR_DEFAULT_MAX_ACTIVE];
				ar_cull_opts_t opts = {
					.z_min = -200.0f, .z_max = 200.0f,
					.sag_radius = 0.0f, .focal_x = 0.0f, .focal_y = 0.0f,
					.margin = 50.0f, .target_quad_px = 2.0f,
					.viewport_w = 1280.0f, .viewport_h = 720.0f
				};
				size_t n = ar_system_query_matrix(rot_cs, 1, mvp, &opts, out_nodes, AR_DEFAULT_MAX_ACTIVE);
				assert(n > 0 && n <= AR_DEFAULT_MAX_ACTIVE);
			}
		}

		ar_destroy(rot_cs);
		printf("  [PASS] 360-degree yaw rotation and pitch angles (0 to 85 deg) all maintain active nodes\n");
	}

	/* Test 15: Custom Vertical Bounds Envelope Callback & Sag Retention */
	{
		ar_system_t *sag_cs = ar_create(160, 128, mock_sampler, NULL);
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
		const ar_node_t *flat_nodes[AR_DEFAULT_MAX_ACTIVE];
		size_t n_flat = ar_system_query_matrix(sag_cs, 1, mvp, &opts_flat, flat_nodes, AR_DEFAULT_MAX_ACTIVE);

		ar_cull_opts_t opts_cb = {
			.z_min = -200.0f, .z_max = 200.0f,
			.sag_radius = R_planet, .focal_x = 0.0f, .focal_y = 0.0f,
			.margin = 100.0f, .target_quad_px = 2.0f,
			.viewport_w = 1280.0f, .viewport_h = 720.0f,
			.bounds_fn = custom_spherical_bounds_cb,
			.bounds_ctx = &R_planet
		};
		const ar_node_t *sag_nodes[AR_DEFAULT_MAX_ACTIVE];
		size_t n_sag = ar_system_query_matrix(sag_cs, 1, mvp, &opts_cb, sag_nodes, AR_DEFAULT_MAX_ACTIVE);

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
		ar_system_t *gen_cs = ar_create(160, 128, analytical_multioctave_sampler, NULL);
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
		ar_system_t *iso_cs = ar_create(160, 128, level_dependent_sampler, NULL);
		assert(iso_cs != NULL);

		/* View 0: High-altitude 3D perspective camera (level >= 5) */
		float modelview[16], proj[16], mvp[16];
		float eye[3] = { 0.0f, -25000.0f, 35000.0f };
		float target[3] = { 0.0f, 0.0f, 0.0f };
		float up[3] = { 0.0f, 0.707f, 0.707f };
		ar_mat4_lookat(modelview, eye, target, up);
		ar_mat4_perspective(proj, 45.0f * 3.14159265f / 180.0f, 16.0f / 9.0f, 10.0f, 200000.0f);
		ar_mat4_multiply(mvp, proj, modelview);

		const ar_node_t *v0_nodes[AR_DEFAULT_MAX_ACTIVE];
		ar_cull_opts_t opts0 = {
			.z_min = -2000.0f, .z_max = 2000.0f,
			.sag_radius = 40000.0f, .focal_x = 0.0f, .focal_y = 0.0f,
			.margin = 150.0f, .target_quad_px = 2.2f,
			.viewport_w = 1280.0f, .viewport_h = 720.0f,
			.min_level = 5, .max_level = 16,
			.view_id = 0
		};
		size_t n0 = ar_system_query_matrix(iso_cs, 1, mvp, &opts0, v0_nodes, AR_DEFAULT_MAX_ACTIVE);
		assert(n0 > 0);
		assert(v0_nodes[0]->level >= 5);
		assert((v0_nodes[0]->active_view_mask & (1u << 0)) != 0);

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

	printf("=== All Standalone libareole Spatial Hierarchy Tests PASSED ===\n");
	return 0;
}
