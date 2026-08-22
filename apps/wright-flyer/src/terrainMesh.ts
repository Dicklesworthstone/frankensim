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

/* ---------- presentation-only dune relief (bead guzez.13+) ----------
 * The surveyed 17x17 grid is honest terrain but nearly flat at 125 m
 * spacing (max 18 m over the whole tile), so the rendered world reads
 * as a table-top. This layer adds DETERMINISTIC, BOUNDED dune relief
 * for the RENDER MESH ONLY — `heightAt` (the physics/presentation
 * height oracle) is untouched, so no visual can move physics. */

/** Deterministic 2D hash in [0, 1) — integer lattice avalanche, no
 * state, same value every call (replays render the same dunes). */
function hash2(ix: number, iy: number): number {
  let h = (Math.imul(ix | 0, 374761393) + Math.imul(iy | 0, 668265263)) | 0;
  h = Math.imul(h ^ (h >>> 13), 1274126177);
  return ((h ^ (h >>> 16)) >>> 0) / 4294967296;
}

/** Smooth (quintic-edge) value noise in [0, 1). */
export function valueNoise2(x: number, y: number): number {
  const ix = Math.floor(x);
  const iy = Math.floor(y);
  const tx = x - ix;
  const ty = y - iy;
  const sx = tx * tx * (3 - 2 * tx);
  const sy = ty * ty * (3 - 2 * ty);
  const a = hash2(ix, iy);
  const b = hash2(ix + 1, iy);
  const c = hash2(ix, iy + 1);
  const d = hash2(ix + 1, iy + 1);
  return a * (1 - sx) * (1 - sy) + b * sx * (1 - sy) + c * (1 - sx) * sy + d * sx * sy;
}

export interface DuneDetailOptions {
  /** Radius around the launch point kept at survey height [m] — the
   * rail corridor, the camp, and the landing flat must stay true. */
  flatRadiusM?: number;
  /** Feather distance beyond that radius where dunes reach full
   * amplitude [m] (smoothstep ramp). */
  featherM?: number;
  /** Full-amplitude crest height [m]. */
  amplitudeM?: number;
}

const DUNE_DEFAULTS: Required<DuneDetailOptions> = {
  flatRadiusM: 55,
  featherM: 130,
  amplitudeM: 3.6,
};

/** Launch-relative dune displacement [m] for the render mesh. Ridged
 * fBm (crest lines stretched along the flight line — wind scours
 * across them) plus one broad swell octave; masked to ZERO inside the
 * launch/camp flat so rails, props, and figures sit exactly where the
 * survey says. Pure: same (x, z) -> same sand, every run. */
export function duneDetail(
  xRel: number,
  zRel: number,
  opts: DuneDetailOptions = {},
): number {
  if (!Number.isFinite(xRel) || !Number.isFinite(zRel)) {
    throw new RangeError(`dune coords must be finite, got ${xRel}, ${zRel}`);
  }
  const o = { ...DUNE_DEFAULTS, ...opts };
  const dist = Math.hypot(xRel, zRel);
  if (dist <= o.flatRadiusM) {
    return 0;
  }
  const t = Math.min(1, (dist - o.flatRadiusM) / o.featherM);
  const mask = t * t * (3 - 2 * t);
  // Crests elongated east-west (x compressed less than z).
  const nx = xRel * 0.0135;
  const nz = zRel * 0.021;
  let ridge = (1 - Math.abs(2 * valueNoise2(nx, nz) - 1)) * 0.68;
  ridge += (1 - Math.abs(2 * valueNoise2(nx * 2.3 + 7.3, nz * 2.3 + 3.1) - 1)) * 0.32;
  const swell = valueNoise2(xRel * 0.006 + 11.7, zRel * 0.006 + 5.9) - 0.5;
  return mask * (o.amplitudeM * ridge + o.amplitudeM * 0.9 * swell);
}

/** Slope/crest shade multiplier in [0.72, 1.14] for vertex tinting:
 * lee faces darken, crests catch light. Samples the SAME detail fn —
 * colors and geometry can never disagree. */
export function duneShade(
  detailAt: (x: number, z: number) => number,
  x: number,
  z: number,
): number {
  const e = 4;
  const dx = detailAt(x + e, z) - detailAt(x - e, z);
  const dz = detailAt(x, z + e) - detailAt(x, z - e);
  const slope = Math.hypot(dx, dz) / (2 * e); // ~tan of the local grade
  const self = detailAt(x, z);
  const crest = Math.max(0, self) * 0.05; // brighten standing crests a touch
  return Math.max(0.72, Math.min(1.14, 1 - Math.min(slope * 0.55, 0.28) + crest));
}

const COLORS: Record<string, [number, number, number]> = {
  water: [0.23, 0.36, 0.45],
  sand: [0.72, 0.66, 0.51],
  dune: [0.78, 0.7, 0.52],
};

/** Dense (res+1)² vertex grid over the tile: positions are three.js
 * world coords CENTERED on the tile (y-up, north = -z); colors follow the
 * material classes. Launch-flat origin: the scene places the rail at the
 * E1.3 launch region's center, which this function reports.
 *
 * `detail` is the OPTIONAL presentation-only relief layer (duneDetail):
 * when given, the render mesh is displaced by it (launch-relative coords)
 * and vertex colors are slope/crest-shaded to match. Omit it and the
 * arrays are exactly the survey surface. */
export function buildTerrainArrays(
  g: TerrainGridJson,
  res: number,
  detail?: (xRel: number, zRel: number) => number,
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
  // Launch region center (E1.3: rows 10-16, cols 3-13 of the 17-grid):
  // row 13, col 8 → north 13*125, east 8*125. Needed BEFORE the loop:
  // the detail layer samples launch-relative coords.
  const launchEast = 8 * g.spacing_m;
  const launchNorth = 13 * g.spacing_m;
  const launch: [number, number, number] = [
    launchEast - half,
    heightAt(g, launchEast, launchNorth),
    -(launchNorth - half),
  ];
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
      const wx = east - half;
      const wz = -(north - half);
      const xRel = wx - launch[0];
      const zRel = wz - launch[2];
      const hSurvey = heightAt(g, east, north);
      const d = detail === undefined ? 0 : detail(xRel, zRel);
      const h = hSurvey + d;
      const k = (j * nv + i) * 3;
      positions[k] = wx;
      positions[k + 1] = h;
      positions[k + 2] = wz;
      let c = COLORS[materialClass(hSurvey)]!;
      if (detail !== undefined && d !== 0) {
        const shade = duneShade(detail, xRel, zRel);
        c = [c[0] * shade, c[1] * shade, c[2] * shade];
      }
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
