import test from "node:test";
import assert from "node:assert/strict";
import {
  getAreoleWasm,
  createAreoleSystem
} from "../dist/index.js";

test("@tty-pt/areole: WebAssembly exports are available and functional", () => {
  const wasm = getAreoleWasm();
  assert.ok(wasm, "WASM instance not loaded");
  assert.equal(typeof wasm.exports.ar_wasm_create, "function");
  assert.equal(typeof wasm.exports.ar_wasm_update, "function");
  assert.equal(typeof wasm.exports.ar_wasm_probe, "function");
  assert.equal(typeof wasm.exports.ar_wasm_destroy, "function");
});

test("@tty-pt/areole: createAreoleSystem lifecycle and bounded active nodes", () => {
  const cs = createAreoleSystem(160, 128);
  assert.ok(cs, "createAreoleSystem returned null");

  const zooms = [4.0, 2.0, 1.0, 0.5, 0.25, 0.1, 0.01, 0.001, 0.0001];
  for (const zoom of zooms) {
    const active = cs.update({
      camX: 0,
      camY: 0,
      zoom,
      viewportW: 1280,
      viewportH: 720
    });
    assert.ok(active.length > 0, `No active nodes at zoom ${zoom}`);
    assert.ok(active.length <= 64, `Active node count ${active.length} exceeded 64 bound at zoom ${zoom}`);
  }

  cs.destroy();
});

test("@tty-pt/areole: hierarchical level scaling with distance / zoom", () => {
  const cs = createAreoleSystem(160, 128);

  const fg = cs.update({ camX: 0, camY: 0, zoom: 4.0, viewportW: 1280, viewportH: 720 });
  const maxFgLevel = Math.max(...fg.map(c => c.level));
  assert.equal(maxFgLevel, 0, "Foreground view should only select Level 0");

  const orbital = cs.update({ camX: 0, camY: 0, zoom: 0.001, viewportW: 1280, viewportH: 720 });
  const maxOrbitalLevel = Math.max(...orbital.map(c => c.level));
  assert.ok(maxOrbitalLevel >= 4, `Orbital view should select high hierarchy levels (got ${maxOrbitalLevel})`);

  cs.destroy();
});

test("@tty-pt/areole: point elevation and attribute probe", () => {
  const cs = createAreoleSystem(160, 128);

  const sample = cs.probe(0, 0);
  assert.ok(sample !== null, "Probe returned null at (0, 0)");
  assert.equal(typeof sample.elevation, "number");
  assert.ok(sample.biome >= 0 && sample.biome <= 255);
  assert.ok(sample.waterDepth >= 0);

  cs.destroy();
});
