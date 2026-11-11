#ifndef TTYPT_AREOLE_GLSL_H
#define TTYPT_AREOLE_GLSL_H

/*
 * ttypt/areole_glsl.h — Canonical GPU Shader Contract & GLSL / C Parity for libareole.
 *
 * LICENSE: BSD-2-Clause
 * Copyright (c) 2026, tty-pt
 *
 * Defines uniform names/layouts, GLSL ES 3.0 functions, and bit/float parity C reference
 * functions for hardware-accelerated spatial and terrain rendering.
 *
 * Zero OpenGL/EGL compile-time or runtime dependencies on server.
 */

#include <stdint.h>
#include <stdbool.h>
#include <math.h>

/* Standard Uniform Names */
#define AR_UNIFORM_ORIGIN        "uAreoleOrigin"
#define AR_UNIFORM_SPAN          "uAreoleSpan"
#define AR_UNIFORM_MVP           "uMVP"
#define AR_UNIFORM_FOCAL_CENTER  "uFocalCenter"
#define AR_UNIFORM_CURVE_RADIUS  "uCurveRadius"

/* Uniform buffer representation in C */
typedef struct {
	float origin[2];       /* uAreoleOrigin: (ox, oy) in world units */
	float span[2];         /* uAreoleSpan: (w, h) in world units */
	float mvp[16];         /* uMVP: Column-major 4x4 matrix */
	float focal_center[2]; /* uFocalCenter: Camera focal center (tx, ty) */
	float curve_radius;    /* uCurveRadius: Spherical curvature radius (0 = flat) */
	float padding;         /* Alignment */
} ar_glsl_uniforms_t;

/* GLSL ES 3.0 Shader Source Functions */
#define AR_GLSL_FUNCTIONS \
	"vec2 ar_uv_to_world(vec2 uv, vec2 origin, vec2 span) {\n" \
	"    return origin + uv * span;\n" \
	"}\n" \
	"vec3 ar_apply_spherical_sag(vec2 world_pos, float base_z, vec2 focal_center, float curve_radius) {\n" \
	"    if (curve_radius <= 0.0) return vec3(world_pos.x, world_pos.y, base_z);\n" \
	"    vec2 dFocal = world_pos - focal_center;\n" \
	"    float r2 = dot(dFocal, dFocal);\n" \
	"    float r = sqrt(r2);\n" \
	"    float curveSag = (r < curve_radius) ? (curve_radius - sqrt(max(1.0, curve_radius * curve_radius - r2))) : (curve_radius + (r - curve_radius) * 2.0);\n" \
	"    return vec3(world_pos.x, world_pos.y, base_z - curveSag);\n" \
	"}\n" \
	"vec4 ar_project(vec3 curved_pos, mat4 mvp) {\n" \
	"    return mvp * vec4(curved_pos, 1.0);\n" \
	"}\n"

static inline const char *ar_glsl_source(void)
{
	return AR_GLSL_FUNCTIONS;
}

/* Matching C Reference Functions for Bit/Float Parity Verification */
static inline void ar_glsl_ref_uv_to_world(const float uv[2], const float origin[2], const float span[2], float out_world[2])
{
	out_world[0] = origin[0] + uv[0] * span[0];
	out_world[1] = origin[1] + uv[1] * span[1];
}

static inline void ar_glsl_ref_apply_spherical_sag(const float world_pos[2], float base_z,
                                                   const float focal_center[2], float curve_radius,
                                                   float out_curved[3])
{
	if (curve_radius <= 0.0f) {
		out_curved[0] = world_pos[0];
		out_curved[1] = world_pos[1];
		out_curved[2] = base_z;
		return;
	}
	float dx = world_pos[0] - focal_center[0];
	float dy = world_pos[1] - focal_center[1];
	float r2 = dx * dx + dy * dy;
	float r = sqrtf(r2);
	float curveSag;
	if (r < curve_radius) {
		float diff = curve_radius * curve_radius - r2;
		if (diff < 1.0f) diff = 1.0f;
		curveSag = curve_radius - sqrtf(diff);
	} else {
		curveSag = curve_radius + (r - curve_radius) * 2.0f;
	}
	out_curved[0] = world_pos[0];
	out_curved[1] = world_pos[1];
	out_curved[2] = base_z - curveSag;
}

static inline void ar_glsl_ref_project(const float curved_pos[3], const float mvp[16], float out_clip[4])
{
	float x = curved_pos[0], y = curved_pos[1], z = curved_pos[2], w = 1.0f;
	out_clip[0] = mvp[0] * x + mvp[4] * y + mvp[8]  * z + mvp[12] * w;
	out_clip[1] = mvp[1] * x + mvp[5] * y + mvp[9]  * z + mvp[13] * w;
	out_clip[2] = mvp[2] * x + mvp[6] * y + mvp[10] * z + mvp[14] * w;
	out_clip[3] = mvp[3] * x + mvp[7] * y + mvp[11] * z + mvp[15] * w;
}

#endif /* TTYPT_AREOLE_GLSL_H */
