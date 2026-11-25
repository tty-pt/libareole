#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ttypt/areole.h>

static void bench_sampler(void *ctx, uint16_t world, uint8_t level,
                          int32_t wx, int32_t wy, ar_sample_t *out)
{
	(void)ctx; (void)world; (void)level;
	int32_t v = ((wx * 13) ^ (wy * 29)) & 511;
	out->val = (int16_t)(v - 256);
	out->rgba = 0xFF408040u;
	out->flags = 0;
	for (uint32_t i = 0; i < AR_SAMPLE_ATTR_CAP; i++)
		out->attr[i] = (uint8_t)((v + (int32_t)i * 17) & 0xFF);
}

static void bench_row_sampler(void *ctx, uint16_t world, uint8_t level,
                              int32_t wy, int32_t ox, int32_t w,
                              uint32_t vert_dim, ar_sample_t *out_row)
{
	uint32_t gdim = (vert_dim > 1u) ? (vert_dim - 1u) : 1u;
	for (uint32_t c = 0; c < vert_dim; c++) {
		int32_t wx = ox + (int32_t)(((int64_t)c * w) / (int64_t)gdim);
		bench_sampler(ctx, world, level, wx, wy, &out_row[c]);
	}
}

static void build_mvp(float cam_x, float cam_y, float dist, float mvp[16])
{
	float eye[3]    = { cam_x, cam_y + dist * 0.8f, dist * 0.6f };
	float target[3] = { cam_x, cam_y, 0.0f };
	float up[3]     = { 0.0f, -0.6f, 0.8f };
	float view[16], proj[16];
	ar_mat4_lookat(view, eye, target, up);
	ar_mat4_perspective(proj, 1.0471975f, 16.0f / 9.0f, dist * 0.02f, dist * 40.0f);
	ar_mat4_multiply(mvp, proj, view);
}

static void print_metrics(const char *label, const ar_metrics_t *m)
{
	double avg_probe = m->frame_map_lookups > 0
		? (double)m->frame_map_probe_steps / (double)m->frame_map_lookups
		: 0.0;
	printf("  %-22s q=%4uus pop=%3u ev=%3u actEv=%u hits=%3u rowC=%4u scC=%5u reuseP=%5u tombs=%u avgProbe=%.2f hyst=%u\n",
	       label, m->us_query_total, m->frame_populates, m->frame_evictions,
	       m->frame_active_evictions, m->frame_hits,
	       m->frame_row_sampler_calls, m->frame_scalar_sampler_calls,
	       m->frame_verts_reused_parent, m->map_tombstones, avg_probe,
	       m->frame_hysteresis_holds);
}

int main(void)
{
	ar_config_t cfg;
	memset(&cfg, 0, sizeof(cfg));
	cfg.base_w = 240;
	cfg.base_h = 240;
	cfg.grid_dim = 48;
	cfg.pool_capacity = 512;
	cfg.samples_pool_capacity = 192;
	cfg.samples_max_level = 2;
	cfg.max_active = 160;
	cfg.max_levels = 16;
	cfg.sampler = bench_sampler;
	cfg.row_sampler = bench_row_sampler;
	cfg.target_quad_px = 10.0f;
	cfg.split_hysteresis = 0.15f;
	cfg.z_min = -100.0f;
	cfg.z_max = 100.0f;

	ar_system_t *sys = ar_system_create_ex(&cfg);
	assert(sys != NULL);
	ar_system_set_timing_enabled(sys, true);

	const ar_node_t *nodes[160];
	ar_cull_opts_t opts;
	memset(&opts, 0, sizeof(opts));
	opts.viewport_w = 1280.0f;
	opts.viewport_h = 720.0f;
	opts.target_quad_px = 10.0f;
	opts.split_hysteresis = 0.15f;
	opts.view_id = 0;

	printf("=== libareole Telemetry & Performance Regression Benchmark ===\n");

	/* 1. Cold start + steady stationary frame */
	float mvp[16];
	build_mvp(0.0f, 0.0f, 180.0f, mvp);
	opts.focal_x = 0.0f; opts.focal_y = 0.0f;
	ar_system_begin_frame(sys);
	size_t n0 = ar_system_query_matrix(sys, 0, mvp, &opts, nodes, 160);
	assert(n0 > 0);
	ar_metrics_t m_cold;
	ar_system_get_metrics(sys, &m_cold);
	print_metrics("1. cold_start", &m_cold);
	assert(m_cold.frame_row_sampler_calls > 0);
	assert(m_cold.map_tombstones == 0);

	ar_system_begin_frame(sys);
	size_t n1 = ar_system_query_matrix(sys, 0, mvp, &opts, nodes, 160);
	assert(n1 == n0);
	ar_metrics_t m_steady;
	ar_system_get_metrics(sys, &m_steady);
	print_metrics("2. steady_stationary", &m_steady);
	assert(m_steady.frame_populates == 0);
	assert(m_steady.frame_evictions == 0);
	assert(m_steady.map_tombstones == 0);

	/* 2. Hysteresis verification with begin_frame() across LOD split boundaries */
	uint32_t total_hyst_holds = 0;
	for (int f = 0; f < 60; f++) {
		float d = (f < 30) ? (220.0f + (float)f * 35.0f) : (220.0f + (float)(59 - f) * 35.0f);
		build_mvp(0.0f, 0.0f, d, mvp);
		ar_system_begin_frame(sys);
		ar_system_query_matrix(sys, 0, mvp, &opts, nodes, 160);
		ar_metrics_t mh;
		ar_system_get_metrics(sys, &mh);
		total_hyst_holds += mh.frame_hysteresis_holds;
	}
	printf("  %-22s total_hysteresis_holds=%u across 60 zoom-sweep frames\n",
	       "3. hysteresis_sweep", total_hyst_holds);
	assert(total_hyst_holds > 0);

	/* 3. Zoom-in from coarse to fine: proves simultaneous parent reuse + odd-row row_sampler batching */
	ar_system_clear(sys);
	build_mvp(5000.0f, 5000.0f, 800.0f, mvp);
	opts.focal_x = 5000.0f; opts.focal_y = 5000.0f;
	ar_system_begin_frame(sys);
	ar_system_query_matrix(sys, 0, mvp, &opts, nodes, 160);

	build_mvp(5000.0f, 5000.0f, 180.0f, mvp);
	ar_system_begin_frame(sys);
	ar_system_query_matrix(sys, 0, mvp, &opts, nodes, 160);
	ar_metrics_t m_zoom;
	ar_system_get_metrics(sys, &m_zoom);
	print_metrics("4. zoom_in_reuse", &m_zoom);
	assert(m_zoom.frame_verts_reused_parent > 0);
	assert(m_zoom.frame_row_sampler_calls > 0);

	/* 4. Sustained 400-frame pan churn (forces thousands of LRU evictions) */
	uint64_t total_ev = 0, total_act_ev = 0, total_lookups = 0, total_steps = 0, total_us = 0;
	for (int f = 0; f < 400; f++) {
		float cx = (float)f * 260.0f;
		float cy = sinf((float)f * 0.1f) * 600.0f;
		build_mvp(cx, cy, 220.0f, mvp);
		opts.focal_x = cx; opts.focal_y = cy;
		ar_system_begin_frame(sys);
		ar_system_query_matrix(sys, 0, mvp, &opts, nodes, 160);
		ar_metrics_t mp;
		ar_system_get_metrics(sys, &mp);
		assert(mp.map_tombstones == 0);
		total_ev += mp.frame_evictions;
		total_act_ev += mp.frame_active_evictions;
		total_lookups += mp.frame_map_lookups;
		total_steps += mp.frame_map_probe_steps;
		total_us += mp.us_query_total;
	}
	double pan_avg_probe = total_lookups > 0 ? (double)total_steps / (double)total_lookups : 0.0;
	printf("  %-22s 400 frames: avg_q=%.1fus evictions=%llu actEv=%llu tombs=0 avgProbe=%.2f\n",
	       "5. pan_churn_400f", (double)total_us / 400.0,
	       (unsigned long long)total_ev, (unsigned long long)total_act_ev, pan_avg_probe);
	assert(total_ev > 1000);
	assert(total_act_ev == 0);
	assert(pan_avg_probe < 2.2);

	/* 5. Single-pass probe hash lookup budget check */
	ar_system_reset_frame_metrics(sys);
	for (int i = 0; i < 1000; i++) {
		int16_t v = 0;
		ar_system_probe_view_val(sys, 0, 0, 100000 + i * 37, -100000 + i * 19, &v);
	}
	ar_metrics_t m_probe;
	ar_system_get_metrics(sys, &m_probe);
	double lookups_per_probe = (double)m_probe.frame_probe_hash_lookups / (double)m_probe.frame_probes;
	printf("  %-22s probes=%u lookups/probe=%.2f fallbacks=%u us=%u\n",
	       "6. miss_probes_1k", m_probe.frame_probes, lookups_per_probe,
	       m_probe.frame_probe_sampler_fallbacks, m_probe.us_probe);
	assert(lookups_per_probe <= 16.0);

	ar_system_destroy(sys);
	printf("PASS: bench_areole\n");
	return 0;
}
