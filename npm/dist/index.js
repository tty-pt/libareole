import { loadWasmModule, getSubstrate } from "@tty-pt/stele";

let wasmInstance = null;

export function getAreoleWasm() {
  if (wasmInstance) return wasmInstance;
  wasmInstance = loadWasmModule(new URL("./areole.wasm", import.meta.url));
  return wasmInstance;
}

export const getChunkWasm = getAreoleWasm;

export function createAreoleSystem(baseW = 160, baseH = 128) {
  const wasm = getAreoleWasm();
  const substrate = getSubstrate();
  const handle = wasm.exports.ar_wasm_create(baseW, baseH);
  if (!handle) {
    throw new Error("Failed to create ar_system in WASM linear memory");
  }

  return {
    update(view, world = 1) {
      const camX = view.camX || 0;
      const camY = view.camY || 0;
      const camZ = view.camZ || 0;
      const pitch = view.pitch || 0;
      const yaw = view.yaw || 0;
      const fov = view.fov || 0;
      const zoom = view.zoom || 1.0;
      const vpW = view.viewportW || 1280;
      const vpH = view.viewportH || 720;
      const effDist = view.effCamDist || 400;

      const activeCount = wasm.exports.ar_wasm_update(
        handle,
        world,
        camX,
        camY,
        camZ,
        pitch,
        yaw,
        fov,
        zoom,
        vpW,
        vpH,
        effDist
      );

      const summaries = [];
      const memBuffer = (wasm.exports.memory || substrate.memory).buffer;
      for (let i = 0; i < activeCount; i++) {
        const nodePtr = wasm.exports.ar_wasm_get_node_ptr(i);
        if (!nodePtr) continue;
        const u32 = new Uint32Array(memBuffer, nodePtr, 16);
        const i32 = new Int32Array(memBuffer, nodePtr, 8);
        const u8 = new Uint8Array(memBuffer, nodePtr, 32);
        const i16 = new Int16Array(memBuffer, nodePtr, 32);

        summaries.push({
          ox: i32[0],
          oy: i32[1],
          w: i32[2],
          h: i32[3],
          level: u8[16],
          hasWater: (u32[7] & 1) !== 0,
          minElevation: i16[16] / 16.0,
          maxElevation: i16[17] / 16.0
        });
      }
      return summaries;
    },

    probe(wx, wy, world = 1) {
      const outPtr = wasm.exports.malloc(20);
      const ok = wasm.exports.ar_wasm_probe(handle, world, wx, wy, outPtr);
      if (!ok) {
        wasm.exports.free(outPtr);
        return null;
      }
      const mem = (wasm.exports.memory || substrate.memory).buffer;
      const i16 = new Int16Array(mem, outPtr, 1);
      const u8 = new Uint8Array(mem, outPtr, 20);
      const u32 = new Uint32Array(mem, outPtr, 5);

      const res = {
        elevation: i16[0] / 16.0,
        biome: u8[12],
        biome2: u8[13],
        blend: u8[14],
        waterDepth: u8[15],
        flags: u32[2]
      };
      wasm.exports.free(outPtr);
      return res;
    },

    clear() {
      wasm.exports.ar_wasm_clear(handle);
    },

    destroy() {
      wasm.exports.ar_wasm_destroy(handle);
    }
  };
}

export const createChunkSystem = createAreoleSystem;
