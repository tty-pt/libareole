#ifndef TTYPT_CHUNK_GLSL_H
#define TTYPT_CHUNK_GLSL_H

/*
 * ttypt/chunk_glsl.h — Backward-compatibility alias header for <ttypt/areole_glsl.h>.
 *
 * LICENSE: BSD-2-Clause
 * Copyright (c) 2026, tty-pt
 */

#include <ttypt/areole_glsl.h>

#define CHUNK_UNIFORM_ORIGIN        "uChunkOrigin"
#define CHUNK_UNIFORM_SPAN          "uChunkSpan"
#define CHUNK_UNIFORM_MVP           AR_UNIFORM_MVP
#define CHUNK_UNIFORM_FOCAL_CENTER  AR_UNIFORM_FOCAL_CENTER
#define CHUNK_UNIFORM_PLANET_RADIUS "uPlanetRadius"

typedef struct {
	float origin[2];
	float span[2];
	float mvp[16];
	float focal_center[2];
	float planet_radius;
	float padding;
} chunk_glsl_uniforms_t;

#define CHUNK_GLSL_FUNCTIONS \
	AR_GLSL_FUNCTIONS \
	"vec2 chunk_uv_to_world(vec2 uv, vec2 origin, vec2 span) { return ar_uv_to_world(uv, origin, span); }\n" \
	"vec3 chunk_apply_spherical_sag(vec2 world_pos, float base_z, vec2 focal_center, float planet_radius) { return ar_apply_spherical_sag(world_pos, base_z, focal_center, planet_radius); }\n" \
	"vec4 chunk_project(vec3 curved_pos, mat4 mvp) { return ar_project(curved_pos, mvp); }\n"

static inline const char *chunk_glsl_source(void)
{
	return CHUNK_GLSL_FUNCTIONS;
}

#define chunk_glsl_ref_uv_to_world         ar_glsl_ref_uv_to_world
#define chunk_glsl_ref_apply_spherical_sag ar_glsl_ref_apply_spherical_sag
#define chunk_glsl_ref_project             ar_glsl_ref_project

#endif /* TTYPT_CHUNK_GLSL_H */
