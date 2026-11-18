<p align="center">
  <img src="https://raw.githubusercontent.com/tty-pt/libareole/main/assets/logo.png" alt="@tty-pt/areole logo" width="220" />
</p>

# @tty-pt/areole

[![npm version](https://img.shields.io/npm/v/@tty-pt/areole.svg)](https://www.npmjs.com/package/@tty-pt/areole)
[![License: BSD-2-Clause](https://img.shields.io/badge/License-BSD--2--Clause-blue)](https://github.com/tty-pt/libareole/blob/main/LICENSE)
[![TypeScript Ready](https://img.shields.io/badge/TypeScript-Ready-3178C6?logo=typescript)](#)
[![Node.js Native](https://img.shields.io/badge/Node.js-N--API-339933?logo=node.js)](#)
[![WebAssembly](https://img.shields.io/badge/WebAssembly-Reactor-654FF0?logo=webassembly)](#)
[![Spatial LOD](https://img.shields.io/badge/spatial-quadtree%20LOD-16A34A)](#)

> Universal Node.js and Browser bindings for **`libareole`**, the distance-adaptive multi-scale spatial quadtree, 6-plane Gribb-Hartmann frustum culling, and LOD hierarchy engine. Canonical C headers: `<ttypt/areole.h>` and `<ttypt/areole_glsl.h>`.

📖 **Full C & Architectural Documentation**: See the main [**libareole GitHub README**](https://github.com/tty-pt/libareole#readme) for single-arena memory layouts (`ar_config_t`), 3D manifold AABB callbacks, and GLSL ES 3.0 vertex parity.

---

## 1. Overview & Domain Concepts

`libareole` (`@tty-pt/areole`) is a domain-agnostic spatial hierarchy and level-of-detail (LOD) culling engine for real-time 3D, 2.5D, and 2D procedural worlds:

1. **Single-Arena Zero-Allocation Traversal**:
   - Allocates the entire quadtree node pool, vertex sample buffers, and spatial hash table in a single contiguous arena at startup (`createAreoleSystem`), performing zero heap allocations during per-frame camera motion.
2. **Rigorous 6-Plane Gribb-Hartmann & 3D Manifold Horizon Culling**:
   - Extracts exact frustum planes from the active camera view (supporting full $360^\circ$ yaw, $0^\circ\text{–}85^\circ$ pitch, perspective/orthographic projection, planetary horizon curvature sag, and pluggable 3D manifold AABB / normal-cone horizon culling callbacks `ar_bounds_3d_fn` and `ar_cull_3d_fn`).
3. **Coarse-First Priority LOD Selection, Scanline Batching & Async Workqueue Population ($O(1)$ Active Budget)**:
   - Guarantees complete 4-quadrant horizon coverage within a bounded `maxActive` chunk budget while subdividing near-camera quads to target screen-space pixel density, with vectorized scanline callbacks (`ar_row_sampler_fn`), watertight cross-LOD seam welding, and background worker population (`<ttypt/wq.h>`).
4. **Hierarchical Point & Attribute Probing (`.probe(wx, wy)`)**:
   - Evaluates continuous interpolated scalar displacement (`value`), packed 32-bit `rgba`, `flags`, and 8-byte material/biome payload (`attr[0..7]`) at any world coordinate.

---

## 2. Installation

```bash
npm install @tty-pt/areole
```

**Native OS Packages (`tty.pt`)**: Prebuilt native C packages (`libareole`) for Linux (`apt`, `apk`, `pacman`, `dnf`), macOS (`brew`), Windows (`winget`, MSYS2), and OpenBSD (`pkg_add`) are available via the [**tty.pt package repositories**](https://github.com/tty-pt/ci/blob/main/docs/install.md).

---

## 3. Quickstart

```typescript
import { createAreoleSystem } from "@tty-pt/areole";

// 1. Allocate a spatial quadtree hierarchy (baseW = 160, baseH = 128)
const sys = createAreoleSystem(160, 128);

// 2. Query visible LOD chunks for the current 3D camera frustum
const nodes = sys.update({
  camX: 0,
  camY: 0,
  camZ: 120,
  pitch: 0.6,
  yaw: 0.0,
  fov: 1.0,
  zoom: 1.0,
  viewportW: 1920,
  viewportH: 1080,
  effCamDist: 400
}, 1 /* worldId */);

console.log("Visible LOD chunks:", nodes.length, nodes[0]);

// 3. Probe exact interpolated surface attributes at world coordinate (wx, wy)
const sample = sys.probe(32, 48, 1);
console.log("Probed sample:", sample);

// 4. Release arena allocation when finished
sys.destroy();
```

---

## 4. JavaScript / TypeScript API Reference

### `createAreoleSystem(baseW?, baseH?, gridDim?, poolCap?, maxActive?, maxLevels?): AreoleSystem`
Creates an `ar_system_t` spatial hierarchy instance and returns an `AreoleSystem` controller:
- `baseW`, `baseH`: Base level-0 chunk dimensions in world units (defaults: `160`, `128`).
- `gridDim`: Quads per chunk edge (default: `32`).
- `poolCap`: Total cached chunk slots in the arena (default: `128`).
- `maxActive`: Maximum active chunks returned per frame query (default: `64`).
- `maxLevels`: Maximum quadtree subdivision depth (default: `16`).

### `AreoleSystem` Methods
- `.update(view: AreoleView, world?: number): AreoleNodeSummary[]`: Runs 6-plane frustum culling and LOD selection, returning `{ ox, oy, w, h, level, flags, minValue, maxValue }` for each active chunk.
- `.probe(wx: number, wy: number, world?: number): AreoleSample | null`: Probes interpolated `{ value, rgba, flags, attr }` at world coordinates `(wx, wy)`.
- `.clear(): void`: Evicts cached quadtree nodes without freeing the arena.
- `.destroy(): void`: Frees the underlying `ar_system_t` arena.

### `getAreoleWasm(): WebAssembly.Instance`
Returns the initialized `areole.wasm` reactor instance on `@tty-pt/stele`'s shared linear memory substrate.

---

## 5. C Development Headers & Native Integration

```makefile
CFLAGS += -I./node_modules/@tty-pt/areole/include
LDLIBS += -lareole -lm
```

---

## License

BSD-2-Clause — Copyright (c) 2026, Paulo Quirino <quirinpa@gmail.com>. See [LICENSE](https://github.com/tty-pt/libareole/blob/main/LICENSE).
