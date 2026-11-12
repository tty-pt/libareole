#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <ttypt/areole.h>
#include <ttypt/iso.h>

/* Built-in procedural continuous sampler for WASM diorama/reactor */
static void default_wasm_sampler(void *user_ctx, uint16_t world, uint8_t level,
                                 int32_t wx, int32_t wy, ar_sample_t *out)
{
	(void)user_ctx;
	(void)world;
	(void)level;

	float fx = (float)wx;
	float fy = (float)wy;

	float macro_h = 120.0f * sinf(fx * 0.0030679f) * cosf(fy * 0.0030679f);
	float fine_h  =  25.0f * sinf(fx * 0.049087f)  * cosf(fy * 0.049087f);
	float detail  =   6.0f * sinf(fx * 0.19635f)   * cosf(fy * 0.19635f);
	float elev = macro_h + fine_h + detail;

	out->val = (int16_t)(elev * 16.0f);
	out->attr[0] = (elev < 0.0f) ? 0 : ((elev > 70.0f) ? 2 : ((elev > 30.0f) ? 1 : 10));
	out->attr[1] = out->attr[0];
	out->attr[2] = 0;
	out->attr[3] = (elev < 0.0f) ? (uint8_t)(-elev > 100.0f ? 100 : -elev) : 0;
	out->rgba = (out->attr[3] > 0) ? 0xFF003080u : 0xFF208020u;
	out->flags = (out->attr[3] > 0) ? 1 : 0;
}

static const ar_node_t **s_active_nodes = NULL;
static size_t s_active_cap = 0;
static size_t s_active_count = 0;

ISO_EXPORT(ar_system_t *, ar_wasm_create, int32_t base_w, int32_t base_h,
           uint32_t grid_dim, uint32_t pool_capacity, uint32_t max_active, uint32_t max_levels) {
	if (max_active > s_active_cap) {
		const ar_node_t ** next = (const ar_node_t **)realloc((void *)s_active_nodes, max_active * sizeof(ar_node_t *));
		if (!next) return NULL;
		s_active_nodes = next;
		s_active_cap = max_active;
	}
	return ar_system_create(base_w, base_h,
	                        (uint16_t)grid_dim,
	                        (uint16_t)pool_capacity,
	                        (uint16_t)max_active,
	                        (uint8_t)max_levels,
	                        default_wasm_sampler, NULL, NULL);
}

ISO_EXPORT(void, ar_wasm_destroy, ar_system_t *sys) {
	if (sys) ar_system_destroy(sys);
}

ISO_EXPORT(void, ar_wasm_clear, ar_system_t *sys) {
	if (sys) ar_system_clear(sys);
}

ISO_EXPORT(uint32_t, ar_wasm_update, ar_system_t *sys, uint32_t world,
                          float cam_x, float cam_y, float cam_z,
                          float pitch, float yaw, float fov,
                          float zoom, float vp_w, float vp_h, float eff_dist) {
	if (!sys || !s_active_nodes || s_active_cap == 0) return 0;
	ar_view_t view = {
		.cam_x = cam_x, .cam_y = cam_y, .cam_z = cam_z,
		.pitch = pitch, .yaw = yaw, .fov = fov,
		.zoom = zoom,
		.viewport_w = vp_w, .viewport_h = vp_h,
		.eff_cam_dist = eff_dist
	};
	s_active_count = ar_system_update(sys, (uint16_t)world, &view, s_active_nodes, s_active_cap);
	return (uint32_t)s_active_count;
}

ISO_EXPORT(uint32_t, ar_wasm_get_node_count, void) {
	return (uint32_t)s_active_count;
}

ISO_EXPORT(const ar_node_t *, ar_wasm_get_node_ptr, uint32_t index) {
	if (index >= s_active_count || !s_active_nodes) return NULL;
	return s_active_nodes[index];
}

ISO_EXPORT(int32_t, ar_wasm_probe, ar_system_t *sys, uint32_t world, int32_t wx, int32_t wy, ar_sample_t *out) {
	if (!sys || !out) return 0;
	return ar_system_probe(sys, (uint16_t)world, wx, wy, out) ? 1 : 0;
}
