export interface AreoleSample {
  elevation: number;
  biome: number;
  biome2: number;
  blend: number;
  waterDepth: number;
  flags: number;
}

export interface AreoleNodeSummary {
  ox: number;
  oy: number;
  w: number;
  h: number;
  level: number;
  hasWater: boolean;
  minElevation: number;
  maxElevation: number;
}

export interface AreoleView {
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

export interface AreoleSystem {
  update(view: AreoleView, world?: number): AreoleNodeSummary[];
  probe(wx: number, wy: number, world?: number): AreoleSample | null;
  clear(): void;
  destroy(): void;
}

export function getAreoleWasm(): any;
export function createAreoleSystem(baseW?: number, baseH?: number): AreoleSystem;
