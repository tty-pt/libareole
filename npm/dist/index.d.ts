export interface AreoleSample {
  value: number;
  rgba: number;
  flags: number;
  attr: number[];
}

export interface AreoleNodeSummary {
  ox: number;
  oy: number;
  w: number;
  h: number;
  level: number;
  flags: number;
  minValue: number;
  maxValue: number;
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
export function createAreoleSystem(
  baseW?: number,
  baseH?: number,
  gridDim?: number,
  poolCap?: number,
  maxActive?: number,
  maxLevels?: number
): AreoleSystem;
