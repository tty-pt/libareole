#ifndef TTYPT_CHUNK_H
#define TTYPT_CHUNK_H

/*
 * ttypt/chunk.h — Backward-compatibility alias header for <ttypt/areole.h>.
 *
 * LICENSE: BSD-2-Clause
 * Copyright (c) 2026, tty-pt
 */

#include <ttypt/areole.h>

#define CHUNK_GRID_DIM       AR_DEFAULT_GRID_DIM
#define CHUNK_VERT_DIM       AR_DEFAULT_VERT_DIM
#define CHUNK_VERT_COUNT     AR_DEFAULT_VERT_COUNT
#define CHUNK_MAX_LEVELS     AR_DEFAULT_MAX_LEVELS
#define CHUNK_MAX_ACTIVE     AR_DEFAULT_MAX_ACTIVE
#define CHUNK_POOL_CAPACITY  AR_DEFAULT_POOL_CAPACITY

typedef ar_sample_t          chunk_sample_t;
typedef ar_sampler_fn        chunk_sampler_fn;
typedef ar_sparsity_fn       chunk_sparsity_fn;
typedef ar_bounds_fn         chunk_bounds_fn;
typedef ar_config_t          chunk_config_t;
typedef ar_node_t            chunk_node_t;
typedef ar_frustum_t         chunk_frustum_t;
typedef ar_cull_opts_t       chunk_cull_opts_t;
typedef ar_view_t            chunk_view_t;
typedef ar_draw_item_t       chunk_draw_item_t;
typedef ar_system_t          chunk_system_t;
typedef ar_system_t          chunk_t;
typedef ar_stats_t           chunk_stats_t;

#define chunk_system_create_ex              ar_system_create_ex
#define chunk_system_create                 ar_system_create
#define chunk_system_destroy                ar_system_destroy
#define chunk_system_clear                  ar_system_clear
#define chunk_system_begin_frame            ar_system_begin_frame
#define chunk_system_eval_view              ar_system_eval_view
#define chunk_mat4_identity                 ar_mat4_identity
#define chunk_mat4_multiply                 ar_mat4_multiply
#define chunk_mat4_ortho                    ar_mat4_ortho
#define chunk_mat4_perspective              ar_mat4_perspective
#define chunk_mat4_lookat                   ar_mat4_lookat
#define chunk_mat4_lookat_lh                ar_mat4_lookat_lh
#define chunk_mat4_invert                   ar_mat4_invert
#define chunk_frustum_from_matrix           ar_frustum_from_matrix
#define chunk_frustum_test_box              ar_frustum_test_box
#define chunk_system_query_matrix           ar_system_query_matrix
#define chunk_system_query_modelview_projection ar_system_query_modelview_projection
#define chunk_system_update                 ar_system_update
#define chunk_system_probe                  ar_system_probe
#define chunk_system_find                   ar_system_find
#define chunk_system_find_by_key            ar_system_find_by_key
#define chunk_system_ensure                 ar_system_ensure
#define chunk_system_get_stats              ar_system_get_stats

#define chunk_node_find                     ar_node_find
#define chunk_create                        ar_create
#define chunk_destroy                       ar_destroy
#define chunk_clear                         ar_clear
#define chunk_begin_frame                   ar_begin_frame
#define chunk_eval_view                     ar_eval_view
#define chunk_update                        ar_update
#define chunk_probe                         ar_probe
#define chunk_elevation_to_float            ar_elevation_to_float
#define chunk_elevation_from_float          ar_elevation_from_float

#endif /* TTYPT_CHUNK_H */
