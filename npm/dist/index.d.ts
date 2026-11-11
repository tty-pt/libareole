export interface ChunkSample {
  elevation: number;
  biome: number;
  biome2: number;
  blend: number;
  waterDepth: number;
  flags: number;
}

export interface ChunkNodeSummary {
  ox: number;
  oy: number;
  w: number;
  h: number;
  level: number;
  hasWater: boolean;
  minElevation: number;
  maxElevation: number;
}

export interface ChunkView {
  camX: number;
  camY: number;
  camZ?: number;
  pitch?: number;
  yaw?: number;
  fov?: number;
  zoom: number;
  viewportW: number;
  viewportH: number;
  effCamDist?: number;
}

export interface ChunkSystem {
  update(view: ChunkView, world?: number): ChunkNodeSummary[];
  probe(wx: number, wy: number, world?: number): ChunkSample | null;
  clear(): void;
  destroy(): void;
}

export function getChunkWasm(): any;
export function createChunkSystem(baseW?: number, baseH?: number): ChunkSystem;
