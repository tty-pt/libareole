<p align="center">
  <img src="../assets/logo.png" alt="@tty-pt/areole logo" width="220" />
</p>

# @tty-pt/areole

[![npm version](https://img.shields.io/badge/npm-v1.0.0-CB3837?logo=npm)](#)
[![License: BSD-2-Clause](https://img.shields.io/badge/License-BSD--2--Clause-blue)](../LICENSE)
[![TypeScript Ready](https://img.shields.io/badge/TypeScript-Ready-3178C6?logo=typescript)](#)
[![WebAssembly](https://img.shields.io/badge/WebAssembly-Reactor-654FF0?logo=webassembly)](#)
[![Spatial LOD](https://img.shields.io/badge/spatial-quadtree%20LOD-16A34A)](#)

> Universal Node.js and Browser WebAssembly bindings for **`libareole`**, the generic distance-adaptive multi-scale spatial quadtree, 6-plane Gribb-Hartmann frustum culling, and LOD hierarchy engine.

---

## 1. Quickstart

```bash
npm install @tty-pt/areole @tty-pt/stele
```

```javascript
import { createAreoleSystem } from "@tty-pt/areole";

// 1. Create a spatial hierarchy system (baseW = 160, baseH = 128)
const sys = createAreoleSystem(160, 128);

// 2. Query active LOD nodes for the current camera frustum
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
  effCamDist: 400,
}, 1 /* worldId */);

console.log("Visible LOD chunks:", nodes.length);

// 3. Probe terrain elevation & biome at world coordinates (wx, wy)
const sample = sys.probe(32, 48, 1);
console.log("Sample:", sample);

// 4. Release WASM arena allocation when finished
sys.destroy();
```

---

## 2. JavaScript / TypeScript API Reference

### `createAreoleSystem(baseW?: number, baseH?: number)`
Allocates an `ar_system_t` instance in WebAssembly linear memory and returns a handle with:
- `.update(view, world?: number)`: Runs 6-plane frustum culling and quadtree LOD selection, returning an array of visible chunk summaries (`{ ox, oy, w, h, level, hasWater, minElevation, maxElevation }`).
- `.probe(wx: number, wy: number, world?: number)`: Probes exact terrain attributes (`{ elevation, biome, biome2, blend, waterDepth, flags }`) at world coordinates `(wx, wy)`.
- `.clear()`: Clears cached quadtree nodes without destroying the arena.
- `.destroy()`: Frees the underlying `ar_system_t` allocation in WASM linear memory.

### `getAreoleWasm(): WebAssembly.Instance`
Returns the initialized `areole.wasm` reactor instance on the `@tty-pt/stele` shared substrate.

---

## 3. Building & Testing

```bash
npm run build
npm test
```

## License

BSD-2-Clause — Author: Paulo Quirino <q@qnixsoft.com>
