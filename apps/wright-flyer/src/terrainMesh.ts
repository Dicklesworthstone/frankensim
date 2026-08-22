// E2.3 terrain core (bead wf-root-guzez.3.3): PURE heightfield sampling +
// vertex/color generation from the committed E1.3 grid — headless-tested;
// the three.js BufferGeometry assembly is a thin consumer. World frame:
// x east, z SOUTH (three.js right-handed y-up: north = -z), y up.

export interface TerrainGridJson {
  grid_n: number;
  spacing_m: number;
  rows_south_to_north: number[][];
}

/** Bilinear height at (east_m, north_m) from the SW corner; clamps to the
 * tile edge (presentation may clamp; PHYSICS refusal lives in Rust). */
export function heightAt(g: TerrainGridJson, eastM: number, northM: number): number {
  const n = g.grid_n;
  const fx = Math.min(Math.max(eastM / g.spacing_m, 0), n - 1);
  const fy = Math.min(Math.max(northM / g.spacing_m, 0), n - 1);
  const c0 = Math.min(Math.floor(fx), n - 2);
  const r0 = Math.min(Math.floor(fy), n - 2);
  const tx = fx - c0;
  const ty = fy - r0;
  const z00 = g.rows_south_to_north[r0]![c0]!;
  const z01 = g.rows_south_to_north[r0]![c0 + 1]!;
  const z10 = g.rows_south_to_north[r0 + 1]![c0]!;
  const z11 = g.rows_south_to_north[r0 + 1]![c0 + 1]!;
  return z00 * (1 - tx) * (1 - ty) + z01 * tx * (1 - ty) + z10 * (1 - tx) * ty + z11 * tx * ty;
}

/** Material class per the E1.3 derivation (battery cross-checks). */
export function materialClass(elevM: number): "water" | "sand" | "dune" {
  return elevM <= 0.05 ? "water" : elevM > 6 ? "dune" : "sand";
}

const COLORS: Record<string, [number, number, number]> = {
  water: [0.23, 0.36, 0.45],
  sand: [0.72, 0.66, 0.51],
  dune: [0.78, 0.7, 0.52],
};

/** Dense (res+1)² vertex grid over the tile: positions are three.js
 * world coords CENTERED on the tile (y-up, north = -z); colors follow the
 * material classes. Launch-flat origin: the scene places the rail at the
 * E1.3 launch region's center, which this function reports. */
export function buildTerrainArrays(
  g: TerrainGridJson,
  res: number,
): {
  positions: Float32Array;
  colors: Float32Array;
  uvs: Float32Array;
  indices: Uint32Array;
  launch: [number, number, number];
} {
  const extent = (g.grid_n - 1) * g.spacing_m; // 2000 m
  const half = extent / 2;
  const nv = res + 1;
  const positions = new Float32Array(nv * nv * 3);
  const colors = new Float32Array(nv * nv * 3);
  // Planar UVs over the tile: a textured material samples NOTHING
  // without them (WebGL zero-fills missing attributes, collapsing the
  // whole map to one texel).
  const uvs = new Float32Array(nv * nv * 2);
  for (let j = 0; j <= res; j++) {
    for (let i = 0; i <= res; i++) {
      const east = (i / res) * extent;
      const north = (j / res) * extent;
      const h = heightAt(g, east, north);
      const k = (j * nv + i) * 3;
      positions[k] = east - half;
      positions[k + 1] = h;
      positions[k + 2] = -(north - half);
      const c = COLORS[materialClass(h)]!;
      colors[k] = c[0];
      colors[k + 1] = c[1];
      colors[k + 2] = c[2];
      const ku = (j * nv + i) * 2;
      uvs[ku] = i / res;
      uvs[ku + 1] = j / res;
    }
  }
  const indices = new Uint32Array(res * res * 6);
  let t = 0;
  for (let j = 0; j < res; j++) {
    for (let i = 0; i < res; i++) {
      const a = j * nv + i;
      const b = a + 1;
      const c = a + nv;
      const d = c + 1;
      indices[t++] = a; indices[t++] = c; indices[t++] = b;
      indices[t++] = b; indices[t++] = c; indices[t++] = d;
    }
  }
  // Launch region center (E1.3: rows 10-16, cols 3-13 of the 17-grid):
  // row 13, col 8 → north 13*125, east 8*125.
  const launchEast = 8 * g.spacing_m;
  const launchNorth = 13 * g.spacing_m;
  const launch: [number, number, number] = [
    launchEast - half,
    heightAt(g, launchEast, launchNorth),
    -(launchNorth - half),
  ];
  return { positions, colors, uvs, indices, launch };
}

/** The arrival-shot camera path (pure in t): a slow descending dolly from
 * over the sound toward the launch flat, settling beside the rail. */
export function arrivalCamera(
  t: number,
  launch: [number, number, number],
): { pos: [number, number, number]; look: [number, number, number] } {
  const s = Math.min(t / 14, 1); // 14-second arrival
  const ease = s * s * (3 - 2 * s);
  const pos: [number, number, number] = [
    launch[0] - 380 + 355 * ease,
    launch[1] + 60 - 55.5 * ease, // settles 4.5 m above the flat
    launch[2] + 240 - 226 * ease,
  ];
  return { pos, look: [launch[0], launch[1] + 1.5, launch[2]] };
}
