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
    .sample_fn = my_sampler,
    .user_ctx = NULL
};
ar_system_t *sys = ar_init(&cfg);

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
ar_free(sys);
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
   - Configurable grid resolution (`grid_dim`), pool capacity (`pool_capacity`), active query budget (`max_active`), default background color (`default_rgba`), and custom allocator callbacks (`alloc_fn`, `free_fn`).
   - Entire `ar_system_t` instance, node pool, contiguous vertex buffers, and spatial hash table are allocated in a single contiguous memory arena at initialization with zero runtime heap allocations during camera motion.
3. **Rigorous Gribb-Hartmann 6-Plane Frustum Culling**:
   - Extracts exact normalized frustum bounding planes directly from any column-major $4 \times 4$ ModelView-Projection (MVP) matrix ($row_3 \pm row_0, row_3 \pm row_1, row_3 \pm row_2$).
   - Seamlessly supports orthographic, perspective, oblique, and rotated cameras across full $360^\circ$ yaw and $0^\circ$ to $85^\circ$ pitch angles.
4. **Pluggable Vertical Envelope Callbacks (`ar_bounds_fn`)**:
   - Supports arbitrary surface curvature or displacement envelopes via `bounds_fn`, plus a built-in spherical curvature sag helper (`ar_spherical_sag_bounds_cb`).
5. **Coarse-First Priority Refinement ($O(1)$ Active Chunks)**:
   - Traversal maintains $\le \text{max\_active}$ chunks across all zoom levels with guaranteed 4-quadrant frustum coverage.
6. **Bilinear Point & Attribute Probe**:
   - Bilinearly interpolates scalar `val`, blends `rgba` color channels, interpolates application `attr[0..7]`, and bitwise-ORs `flags` with $O(1)$ hierarchical ancestor lookup.
7. **GLSL ES 3.0 Uniform & Shader Parity (`<ttypt/areole_glsl.h>`)**:
   - Standardized uniform layouts and GLSL functions with bit/float parity tested C reference implementations.

---

## 3. Public C API (`<ttypt/areole.h>`)

### Core Types

```c
#include <ttypt/areole.h>

/* Generic sample evaluated at each grid vertex */
typedef struct {
    int16_t  val;        /* Generic scalar displacement / field value */
    uint32_t rgba;       /* Packed 32-bit color or RGBA attributes */
    uint32_t flags;      /* Generic bitmask flags (e.g. discrete, surface type) */
    uint8_t  attr[8];    /* Generic 8-byte application payload (materials, blend) */
} ar_sample_t;

/* Runtime configuration for ar_init(&cfg) */
typedef struct {
    int32_t       base_w, base_h;   /* Base level-0 spatial extent in world units */
    uint16_t      grid_dim;         /* Quads per chunk edge (0 = default 32) */
    uint16_t      pool_capacity;    /* Total cached chunk slots (0 = default 128) */
    uint16_t      max_active;       /* Max active chunks per query (0 = default 64) */
    uint32_t      default_rgba;     /* Fallback RGBA for empty probes (0 = 0xFF000000u) */
    ar_sample_fn  sample_fn;        /* Point-sampling callback */
    ar_fill_fn    fill_fn;          /* Optional bulk-fill callback */
    void         *user_ctx;         /* Opaque context pointer passed to callbacks */
    void        *(*alloc_fn)(size_t);
    void         (*free_fn)(void *);
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
