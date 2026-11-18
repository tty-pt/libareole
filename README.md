<p align="center">
  <img src="assets/logo.png" alt="libareole logo" width="220" />
</p>

# libareole

[![npm](https://img.shields.io/badge/npm-@tty--pt/areole-CB3837?logo=npm)](npm/README.md)
[![C99](https://img.shields.io/badge/C-C99-555?logo=c)](#)
[![BSD-2-Clause](https://img.shields.io/badge/License-BSD--2--Clause-blue)](#)
[![Spatial LOD](https://img.shields.io/badge/spatial-quadtree%20LOD-16A34A)](#)
[![Frustum Culling](https://img.shields.io/badge/culling-6--plane%20Gribb--Hartmann-blueviolet)](#)

> Generic multi-scale spatial hierarchy, 6-plane frustum culling, and quadtree LOD engine in pure C. Canonical headers are `<ttypt/areole.h>` and `<ttypt/areole_glsl.h>`.

Pure C (C99/C11), zero-dependency, standalone spatial hierarchy engine. Designed as the generic multi-scale spatial index and frustum culling engine for real-time 3D, 2.5D, and 2D procedural terrain, planetary rendering, and continuous simulation.

---

## 1. Quickstart

Initialize a single-arena spatial hierarchy and query visible LOD chunks against a $4 \times 4$ camera ModelView-Projection (MVP) matrix:

```c
#include <ttypt/areole.h>

/* 1. Create system at startup */
ar_config_t cfg = {
    .base_w = 160,
    .base_h = 128,
    .grid_dim = 32,
    .pool_capacity = 128,
    .max_active = 64,
    .max_levels = 16,
    .sampler = my_sampler,
    .user_ctx = NULL
};
ar_system_t *sys = ar_system_create_ex(&cfg);

/* 2. Query using exact camera MVP matrix each frame */
float modelview[16], proj[16], mvp[16];
ar_mat4_lookat(modelview, eye, target, up);
ar_mat4_perspective(proj, fov_rad, aspect, near_z, far_z);
ar_mat4_multiply(mvp, proj, modelview);

ar_cull_opts_t opts = {
    .z_min = -2000.0f,
    .z_max = 2000.0f,
    .sag_radius = 24000.0f,
    .focal_x = cam_x,
    .focal_y = cam_y,
    .margin = 150.0f,
    .target_quad_px = 2.0f,
    .viewport_w = 1920.0f,
    .viewport_h = 1080.0f
};

const ar_node_t *visible[64];
size_t count = ar_system_query_matrix(sys, world_id, mvp, &opts, visible, 64);

/* 3. Clean up on shutdown */
ar_system_destroy(sys);
```

### Build & Link
```sh
cc app.c $(pkg-config --cflags --libs areole) -lm -o app
```

---

## 2. Architectural Highlights

1. **Domain-Agnostic Spatial Primitives (`ar_` API)**:
   - Zero game-specific concepts (elevation, biomes, water, structures) in the core engine.
   - Generic sample points provide a 16-bit scalar displacement (`val`), packed 32-bit color/attributes (`rgba`), 32-bit bitmask (`flags`), and an 8-byte application payload (`attr[8]`).
2. **Runtime-Configurable Single-Arena Allocation (`ar_config_t`)**:
   - Configurable grid resolution (`grid_dim`), pool capacity (`pool_capacity`), active query budget (`max_active`), hierarchy depth (`max_levels`), value scale (`value_scale`), and toroidal periods (`period_x`, `period_y`).
   - Entire `ar_system_t` instance, node pool, contiguous vertex buffers, and spatial hash table are allocated in a single contiguous memory arena at initialization with zero runtime heap allocations during camera motion.
3. **Rigorous Gribb-Hartmann 6-Plane Frustum Culling**:
   - Extracts exact normalized frustum bounding planes directly from any column-major $4 \times 4$ ModelView-Projection (MVP) matrix ($row_3 \pm row_0, row_3 \pm row_1, row_3 \pm row_2$).
   - Seamlessly supports orthographic, perspective, oblique, and rotated cameras across full $360^\circ$ yaw and $0^\circ$ to $85^\circ$ pitch angles.
4. **Pluggable Vertical, 3D Manifold Envelope & Horizon Culling Callbacks (`ar_bounds_fn`, `ar_bounds_3d_fn`, `ar_cull_3d_fn`)**:
   - Supports arbitrary 2.5D vertical envelopes (`bounds_fn`), full 3D non-Euclidean manifold bounding boxes (`bounds_3d_fn`), and manifold normal-cone/horizon visibility culling (`cull_3d_fn`) for curved or planetary surfaces.
5. **Coarse-First Priority Refinement, Scanline Batching & Async Workqueue Population ($O(1)$ Active Chunks)**:
   - Traversal maintains $\le \text{max\_active}$ chunks across all zoom levels with guaranteed 4-quadrant frustum coverage, optional vectorized scanline sampling (`ar_row_sampler_fn`), per-frame populate caps (`max_populates_per_frame`), and non-blocking background worker population via `<ttypt/wq.h>` (`cfg.wq`, `AR_NODE_FLAG_DEFERRED`).
6. **Configurable Point & Attribute Probe**:
   - Interpolates scalar `val`, blends `rgba` color channels, and interpolates application `attr[0..7]` using configurable curves (`AR_INTERP_SMOOTHSTEP`, `AR_INTERP_LINEAR`, `AR_INTERP_QUINTIC`, `AR_INTERP_NEAREST`, or custom `ar_interp_fn`) with $O(1)$ hierarchical ancestor lookup.
7. **GLSL ES 3.0 Uniform & Shader Parity (`<ttypt/areole_glsl.h>`)**:
   - Standardized uniform layouts and GLSL functions (`ar_apply_spherical_sag_ex`, `ar_glsl_ref_eval_vertex`) with bit/float parity tested C reference implementations.

---

## 3. Public C API (`<ttypt/areole.h>`)

### Core Types

```c
#include <ttypt/areole.h>

/* Generic sample evaluated at each grid vertex */
typedef struct {
    int16_t  val;        /* Generic scalar displacement / field value */
    uint32_t rgba;       /* Packed 32-bit color or RGBA attributes */
    uint32_t flags;      /* Generic bitmask flags (bits 0..29 user-defined; bit 30 AR_NODE_FLAG_DEFERRED; bit 31 AR_NODE_FLAG_SPARSE) */
    uint8_t  attr[8];    /* Generic 8-byte application payload (materials, blend) */
} ar_sample_t;

/* Runtime configuration for ar_system_create_ex(&cfg) */
typedef struct {
    int32_t           base_w, base_h;   /* Base level-0 spatial extent in world units */
    uint16_t          grid_dim;         /* Quads per chunk edge (required > 0) */
    uint16_t          pool_capacity;    /* Total cached chunk slots (required > 0) */
    uint16_t          max_active;       /* Max active chunks per query (required > 0) */
    uint8_t           max_levels;       /* Max hierarchy levels (required > 0) */
    uint32_t          default_rgba;     /* Fallback RGBA when sample->rgba == 0 */
    float             z_min, z_max;     /* Default vertical bounds envelope */
    float             target_quad_px;   /* Default screen-space quad subdivision threshold */
    float             margin;           /* Default world-space culling margin */
    float             value_scale;      /* Fixed-point scalar divisor (0 = 16.0f) */
    int32_t           period_x, period_y; /* Optional toroidal world periods */
    ar_sampler_fn     sampler;          /* Required vertex sampler callback */
    ar_row_sampler_fn row_sampler;      /* Optional vectorized scanline/row sampler */
    ar_sparsity_fn    sparsity;         /* Optional discrete sparsity callback */
    void             *user_ctx;         /* Opaque context pointer passed to callbacks */
    uint32_t          max_populates_per_frame; /* Optional per-frame new-chunk populate budget */
    ar_bounds_fn      bounds_fn;        /* Optional vertical envelope callback */
    ar_bounds_3d_fn   bounds_3d_fn;     /* Optional 3D manifold AABB callback */
    ar_cull_3d_fn     cull_3d_fn;       /* Optional 3D manifold normal-cone/horizon culling callback */
    void             *bounds_ctx;       /* Opaque context passed to bounds/cull callbacks */
    struct stl_wq    *wq;               /* Optional stele workqueue for async background chunk population */
} ar_config_t;
```

---

## 4. Installation & Build

### Prebuilt Packages
Prebuilt packages for Linux (APT / Alpine / Arch / Fedora-RHEL), macOS (Homebrew), Windows (winget / MSYS2), and OpenBSD are distributed via [tty.pt packages](https://github.com/tty-pt/ci/blob/main/docs/install.md). Use **libareole** as the package name.

### Build from Source
```sh
make              # builds lib/libareole.so and lib/libchunk.so
make wasm         # builds WebAssembly reactor (npm/dist/areole.wasm)
make test         # runs standalone C unit and GLSL parity test suite
make npm-test     # runs Node.js WebAssembly integration test suite
sudo make install # install headers and libraries to system
```

---

## 5. JavaScript & WebAssembly Bindings (`@tty-pt/areole`)

For browser, Node.js, and WebGL2 engines, `libareole` ships an isomorphic WebAssembly spatial hierarchy and culling package on npm:

```sh
npm install @tty-pt/areole
```

See [npm/README.md](npm/README.md) for the full TypeScript and WebAssembly API reference.

---

## License

BSD 2-Clause License — Copyright (c) 2026, tty-pt. See `LICENSE`.
