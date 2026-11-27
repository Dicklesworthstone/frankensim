/**
 * Kept hot until lunch: yesterday's chili, reheated and packed hot in insulated carriers, then carried to a lunch at
 * 11:30. Four carriers side by side: six pans in one and a single pan in another, each either primed first with a pan
 * of near-boiling water for half an hour or loaded cold from the shelf.
 *
 * Each pan of chili conducts its heat out cell by cell (a quarter of the pan, by symmetry, in a three-dimensional grid)
 * to its faces. Its top gives heat to its lid across the air under it, and the steam condenses on the lid's underside.
 * The lid, the pan's floor and its walls each pass heat across their air layer to whatever faces them: another pan in
 * the next slot, or the carrier's own inner wall. That wall is a moulded polyethylene skin, foamed between it and the
 * outer skin, and the heat crosses it into the room. Each pan also hangs by its steel rim on the moulded rails, and
 * heat runs up its wall into them. Nothing is fitted: the chili's properties are Choi and Okos's, the air layers'
 * the textbook correlations, the shell's its materials'. US patent 7,134,552 reports that a carrier of this kind
 * holds food heated to 185 °F above 141 °F "for at least about six hours" and "as much as fifteen hours". Here a lone
 * pan in a primed carrier holds its middle above 141 °F for about six; six pans together hold their coolest middle
 * above it for about twenty, since the middle pans only lose heat sideways, through 15 cm of chili, and no air leaks
 * past the door's gasket here (sim/holdHot.test.ts).
 */
import { CHILI, EPS_STEEL, across, fToC, foodProps } from "./coolDown";

const SIGMA = 5.670374419e-8;
const INCH = 0.0254;

/** The carrier as drawn (three/models/support.ts, NSF-CARRIER-INSULATED): the cavity and the shell outside it (m). */
export const CAVITY = { w: 0.36, h: 0.525, d: 0.54, floorY: 0.06 } as const;
export const OUTER = { w: 0.43, h: 0.605, d: 0.635 } as const;
/** Each slot's pan rim, m up from the carrier's base, bottom slot first. */
export const RIM_Y: readonly number[] = [0, 1, 2, 3, 4, 5].map((i) => 0.135 + 0.085 * i);
/** A GN 1/1 hotel pan 65 mm deep with 2 in of chili, 7 litres (as sim/coolDown.ts's PAN_2), its flat lid on. */
export const PAN = { panM: 0.065, foodM: 2 * INCH, w: 0.29, l: 0.475, lidW: 0.329, lidL: 0.534 } as const;
/** The pan's floor below its rim; its lid's top above it; the air between the chili's top and the lid. */
const FLOOR_DROP = 0.0655;
const LID_UP = 0.002;
const HEADSPACE = 0.064 - PAN.foodM;
/** Stainless 0.8 mm: its heat capacity per m² (J/m²·K). */
const STEEL_C = 7900 * 500 * 0.0008;

/** The shell: a polyethylene skin each side of rigid polyurethane foam. */
const SKIN = { t: 0.0035, k: 0.4, rc: 940 * 2300 } as const;
const FOAM = { t: 0.03, k: 0.03, rc: 40 * 1400, n: 6 } as const;
export const WALL_T = 2 * SKIN.t + FOAM.t;
/** Each of the wall's nodes' depth from the cavity (m): the inner skin, the foam's cells, the outer skin. */
export const WALL_POS: readonly number[] = [SKIN.t / 2, ...Array.from({ length: FOAM.n }, (_, i) => SKIN.t + (FOAM.t / FOAM.n) * (i + 0.5)), WALL_T - SKIN.t / 2];
const A_IN = 2 * (CAVITY.w * CAVITY.h + CAVITY.w * CAVITY.d + CAVITY.h * CAVITY.d);
const A_OUT = 2 * (OUTER.w * OUTER.h + OUTER.w * OUTER.d + OUTER.h * OUTER.d);
const EPS_PE = 0.9;

/** The room the carriers stand in, and the van: 72 °F. */
export const ROOM_C = fToC(72);
/** Priming: a 2-inch pan of water near the boil in the bottom slot, door shut, this long before the chili goes in. */
export const PRIME_MIN = 30;
const PRIME_C = fToC(190);
/** What a pan passes (W/K) up its long walls to its rim and into the rails it rests on: the 0.8 mm steel up 20 mm of
 *  wall along both sides (0.64 W/K), the rim's 12 mm on the rails at about 200 W/m²·K (2.4), and the rail itself into
 *  the skin (about 1.5), in series. */
export const RIM_G = 0.38;

export const HOT_LINE_C = fToC(135);
export const SAMPLE_MIN = 5;
/** Steps (s). The chili's cells are 6.4 mm deep: its explicit step is stable below about 135 s, the lids and the shell
 *  are stepped implicitly. A finer grid at 10 s gives the same hours (sim/holdHot.test.ts). */
const STEP_S = 60;
/** Where the edge is read: an inch in from a pan's side and its end, halfway down. */
const EDGE_IN = INCH;
/** The pans' cut faces are kept this long after packing (the stage shows up to 7 h from 6:00 to 13:00). */
export const SLICE_MIN = 8 * 60;

/** The grid over a quarter of the chili: across half its width, along half its length, through its depth. */
export const NX = 6;
export const NZ = 6;
export const NY = 8;

export type HoldKey = "fullPrimed" | "fullCold" | "onePrimed" | "oneCold";
export interface HoldWay {
  key: HoldKey;
  pans: number;
  primed: boolean;
  name: string;
  short: string;
}
export const HOLD_WAYS: readonly HoldWay[] = [
  { key: "fullPrimed", pans: 6, primed: true, name: "Six pans, primed", short: "Full, primed" },
  { key: "fullCold", pans: 6, primed: false, name: "Six pans, cold carrier", short: "Full, cold" },
  { key: "onePrimed", pans: 1, primed: true, name: "One pan, primed", short: "1 pan, primed" },
  { key: "oneCold", pans: 1, primed: false, name: "One pan, cold carrier", short: "1 pan, cold" },
];

export interface HoldRun {
  way: HoldWay;
  startC: number;
  /** The cut's cells: across from the middle out, and up. */
  grid: [number, number];
  /** Samples every SAMPLE_MIN from −PRIME_MIN (priming, or the carrier standing empty) to `endMin`. */
  endMin: number;
  /** Each pan's middle (°C), bottom slot first. */
  centre: Float32Array[];
  /** The edge: where a probe an inch in from a pan's corner, halfway down, reads; the coldest pan's. The chili's
   *  coldest spot is at its corners, in a cell's width of them, so a fixed point a probe can find is what is judged. */
  coldest: Float32Array;
  /** The wall from inside out: inner skin, the foam's cells, outer skin. */
  wall: Float32Array;
  /** The priming water (°C), NaN when not priming. */
  water: Float32Array;
  /** Each pan's cut through its middle, `grid` per pan (across from the middle out, bottom up), to SLICE_MIN. */
  slice: Float32Array;
  /** Minutes after packing when the edge falls below 135 °F, and when the coolest pan's middle falls below 141 °F;
   *  null if not within the run (0 if packed under the line). */
  hotMin: number | null;
  mid141: number | null;
  /** Energy (J): stored at packing and at the end, and what left the carrier in between. */
  energy: { at0: number; atEnd: number; lost: number };
}
export const WALL_NODES = FOAM.n + 2;
export const samplesTo = (min: number) => Math.round((min + PRIME_MIN) / SAMPLE_MIN) + 1;
/** The sample index at minute `m` after packing (fractional). */
export const sampleAt = (m: number) => (m + PRIME_MIN) / SAMPLE_MIN;

/* Air and its layers. */
const airK = (k: number) => 0.02587 * (k / 293.15) ** 0.85;
const airNu = (k: number) => 1.516e-5 * (k / 293.15) ** 1.75;
const radH = (a: number, b: number, eps: number) => {
  const ka = a + 273.15;
  const kb = b + 273.15;
  return eps * SIGMA * (ka * ka + kb * kb) * (ka + kb);
};
const epsPair = (a: number, b: number) => 1 / (1 / a + 1 / b - 1);
/** Heat (W/m²·K) across a horizontal air layer `gap` deep, `lo` the lower face, `hi` the upper: Hollands's
 *  correlation when heated from below, conduction alone when the warmer face is on top; and radiation. */
function layerH(lo: number, hi: number, gap: number, eps: number): number {
  const fm = (lo + hi) / 2 + 273.15;
  const nu = airNu(fm);
  const ra = lo > hi ? ((9.80665 / fm) * (lo - hi) * gap ** 3) / ((nu * nu) / 0.71) : 0;
  const nus = 1 + 1.44 * Math.max(0, 1 - 1708 / Math.max(ra, 1)) + Math.max(0, Math.cbrt(ra / 5830) - 1);
  return (nus * airK(fm)) / gap + radH(lo, hi, eps);
}
/** Heat (W/m²·K) from a vertical face `h` m tall to a wall `gap` away across still air: the face's film and the
 *  wall's in series (each at half the difference), or conduction across the gap if that is more; and radiation. */
function sideH(t: number, wall: number, h: number, gap: number, eps: number): number {
  const fm = (t + wall) / 2 + 273.15;
  const film = 1.42 * (Math.abs(t - wall) / 2 / h) ** 0.25;
  return Math.max(film / 2, airK(fm) / gap) + radH(t, wall, eps);
}
/** The carrier's outside to the room: free convection on its walls and radiation. */
const outerH = (t: number) => 1.42 * (Math.abs(t - ROOM_C) / OUTER.h) ** 0.25 + radH(t, ROOM_C, EPS_PE);

/* The chili's properties by temperature, every quarter degree from 0 to 100 °C. */
const TAB_N = 401;
const TAB_K = new Float64Array(TAB_N);
const TAB_RC = new Float64Array(TAB_N);
for (let i = 0; i < TAB_N; i++) {
  const p = foodProps(CHILI, i * 0.25);
  TAB_K[i] = p.k;
  TAB_RC[i] = p.rho * p.cp;
}
const tab = (a: Float64Array, c: number) => {
  const x = Math.min(TAB_N - 1.000001, Math.max(0, c / 0.25));
  const i = Math.floor(x);
  return a[i]! + (x - i) * (a[i + 1]! - a[i]!);
};
/** The chili's heat content per m³ from 0 °C (J/m³), by the table. */
function enthalpy(c: number): number {
  let h = 0;
  const n = Math.max(1, Math.ceil(Math.abs(c) / 0.5));
  const d = c / n;
  for (let i = 0; i < n; i++) h += tab(TAB_RC, (i + 0.5) * d) * d;
  return h;
}

/** Where each loaded pan stands: its slot, and what its floor and lid face across their air layers. */
interface Slot {
  slot: number;
  /** The pan below's index (its lid faces this pan's floor), or −1: the cavity's floor. */
  below: number;
  above: number;
  gapBelow: number;
  gapAbove: number;
}
function slotsOf(pans: number): Slot[] {
  const slots = Array.from({ length: pans }, (_, i) => i);
  return slots.map((s, i) => {
    const floorY = RIM_Y[s]! - FLOOR_DROP;
    const lidY = RIM_Y[s]! + LID_UP;
    const below = i > 0 ? i - 1 : -1;
    const above = i < pans - 1 ? i + 1 : -1;
    return {
      slot: s,
      below,
      above,
      gapBelow: below >= 0 ? floorY - (RIM_Y[slots[below]!]! + LID_UP) : floorY - CAVITY.floorY,
      gapAbove: above >= 0 ? RIM_Y[slots[above]!]! - FLOOR_DROP - lidY : CAVITY.floorY + CAVITY.h - lidY,
    };
  });
}

export interface HoldOpts {
  /** Run at least this long after packing (min); on until the edge is under 135 °F, but no longer than
   *  `maxMin`. */
  minMin?: number;
  maxMin?: number;
  /** The grid (cells across, through, along a quarter pan) and the step (s): for checking they are fine enough. */
  grid?: readonly [number, number, number];
  dt?: number;
  /** Stop once the middles are under 141 °F as well (for fitting). */
  to141?: boolean;
}

export function holdRun(way: HoldWay, startC: number, o: HoldOpts = {}): HoldRun {
  const [nx, ny, nz] = o.grid ?? [NX, NY, NZ];
  const dt = o.dt ?? STEP_S;
  const dx = PAN.w / 2 / nx;
  const dy = PAN.foodM / ny;
  const dz = PAN.l / 2 / nz;
  const minMin = o.minMin ?? SLICE_MIN;
  const maxMin = o.maxMin ?? 40 * 60;
  const P = way.pans;
  const slots = slotsOf(P);
  const NC = nx * ny * nz;
  const at = (i: number, j: number, k: number) => (k * ny + j) * nx + i;
  // Each pan's quarter: temperatures, and the lid.
  const T = Array.from({ length: P }, () => new Float64Array(NC).fill(startC));
  const dT = new Float64Array(NC);
  const lid = new Float64Array(P).fill(startC);
  const A_FOOT = PAN.w * PAN.l;
  /** The pan's floor, under the chili's footprint: what faces the lid below across their air layer. */
  const A_FLOOR = (PAN.lidW - 0.03) * (PAN.lidL - 0.03);
  const A_LID = PAN.lidW * PAN.lidL;
  const C_LID = A_LID * STEEL_C;
  // The wall's nodes: inner skin, foam cells, outer skin; their places (from the inside) and areas there.
  const NW = WALL_NODES;
  const wPos = new Float64Array(NW);
  const wT = new Float64Array(NW).fill(ROOM_C);
  const wC = new Float64Array(NW);
  const wG = new Float64Array(NW - 1);
  const areaAt = (s: number) => A_IN + ((A_OUT - A_IN) * s) / WALL_T;
  wPos[0] = SKIN.t / 2;
  for (let i = 0; i < FOAM.n; i++) wPos[1 + i] = SKIN.t + (FOAM.t / FOAM.n) * (i + 0.5);
  wPos[NW - 1] = WALL_T - SKIN.t / 2;
  wC[0] = SKIN.rc * SKIN.t * areaAt(wPos[0]!);
  for (let i = 1; i < NW - 1; i++) wC[i] = FOAM.rc * (FOAM.t / FOAM.n) * areaAt(wPos[i]!);
  wC[NW - 1] = SKIN.rc * SKIN.t * areaAt(wPos[NW - 1]!);
  const half = (i: number) => (i === 0 || i === NW - 1 ? [SKIN.t / 2, SKIN.k] : [FOAM.t / FOAM.n / 2, FOAM.k]) as [number, number];
  for (let i = 0; i < NW - 1; i++) {
    const s = i === 0 ? SKIN.t : i === NW - 2 ? WALL_T - SKIN.t : wPos[i]! + FOAM.t / FOAM.n / 2;
    const [ta, ka] = half(i);
    const [tb, kb] = half(i + 1);
    wG[i] = areaAt(s) / (ta / ka + tb / kb);
  }
  const triA = new Float64Array(NW);
  const triB = new Float64Array(NW);
  const triC = new Float64Array(NW);
  const triD = new Float64Array(NW);
  let lost = 0;
  /** One implicit step of the wall, `q` (W) into its inner skin from the cavity. */
  function wallStep(q: number): void {
    // The outer skin's loss to the room, at the step's start: the same in the solve and in what is counted lost.
    const gOut = outerH(wT[NW - 1]!) * A_OUT;
    for (let i = 0; i < NW; i++) {
      const gl = i > 0 ? wG[i - 1]! : 0;
      const gr = i < NW - 1 ? wG[i]! : gOut;
      triA[i] = -gl;
      triC[i] = i < NW - 1 ? -gr : 0;
      triB[i] = wC[i]! / dt + gl + gr;
      triD[i] = (wC[i]! / dt) * wT[i]! + (i === 0 ? q : 0) + (i === NW - 1 ? gr * ROOM_C : 0);
    }
    for (let i = 1; i < NW; i++) {
      const m = triA[i]! / triB[i - 1]!;
      triB[i] = triB[i]! - m * triC[i - 1]!;
      triD[i] = triD[i]! - m * triD[i - 1]!;
    }
    wT[NW - 1] = triD[NW - 1]! / triB[NW - 1]!;
    for (let i = NW - 2; i >= 0; i--) wT[i] = (triD[i]! - triC[i]! * wT[i + 1]!) / triB[i]!;
    lost += gOut * (wT[NW - 1]! - ROOM_C) * dt;
  }

  const n = samplesTo(maxMin);
  const centre = Array.from({ length: P }, () => new Float32Array(n).fill(NaN));
  const coldest = new Float32Array(n).fill(NaN);
  const wall = new Float32Array(n * NW).fill(NaN);
  const water = new Float32Array(n).fill(NaN);
  const nSlice = samplesTo(SLICE_MIN);
  const slice = new Float32Array(nSlice * P * nx * ny).fill(NaN);
  let k0 = 0;

  // Priming: the water pan in the bottom slot, open, the door shut. An empty carrier stands at the room's warmth.
  const steps = Math.round((PRIME_MIN * 60) / dt);
  const perSample = Math.round((SAMPLE_MIN * 60) / dt);
  let tw = PRIME_C;
  const C_W = 7.0 * 4186 + (A_FOOT + 2 * (PAN.w + PAN.l) * PAN.panM) * STEEL_C;
  const s0 = slotsOf(1)[0]!;
  const record = (k: number, mAfter: number) => {
    if (k >= n) return;
    for (let i = 0; i < NW; i++) wall[k * NW + i] = wT[i]!;
    if (mAfter < 0) {
      if (way.primed) water[k] = tw;
      return;
    }
    let cold = Infinity;
    // The edge's place in the grid's cells (cell centres at (i + ½)·dx): an inch in from the side and the end.
    const fx = Math.min(nx - 1, Math.max(0, (PAN.w / 2 - EDGE_IN) / dx - 0.5));
    const fz = Math.min(nz - 1, Math.max(0, (PAN.l / 2 - EDGE_IN) / dz - 0.5));
    const fy = ny / 2 - 0.5;
    const [i0, j0, k1] = [Math.floor(fx), Math.floor(fy), Math.floor(fz)];
    const [i1, j1, k2] = [Math.min(nx - 1, i0 + 1), Math.min(ny - 1, j0 + 1), Math.min(nz - 1, k1 + 1)];
    const [ax, ay, az] = [fx - i0, fy - j0, fz - k1];
    for (let p = 0; p < P; p++) {
      const t = T[p]!;
      centre[p]![k] = (t[at(0, ny / 2 - 1, 0)]! + t[at(0, ny / 2, 0)]!) / 2;
      const lerp = (a: number, b: number, f: number) => a + (b - a) * f;
      const plane = (k: number) => lerp(lerp(t[at(i0, j0, k)]!, t[at(i1, j0, k)]!, ax), lerp(t[at(i0, j1, k)]!, t[at(i1, j1, k)]!, ax), ay);
      cold = Math.min(cold, lerp(plane(k1), plane(k2), az));
      if (k < nSlice) for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) slice[((k * P + p) * ny + j) * nx + i] = t[at(i, j, 0)]!;
    }
    coldest[k] = cold;
  };
  for (let s = 0; s < steps; s++) {
    if (s % perSample === 0) record(k0++, -PRIME_MIN + (s * dt) / 60);
    if (!way.primed) continue;
    const ts = wT[0]!;
    const sides = 2 * (PAN.w + PAN.l) * PAN.foodM;
    const q =
      across(tw, ts, CAVITY.floorY + CAVITY.h - (RIM_Y[0]! + LID_UP) + HEADSPACE) * A_FOOT +
      layerH(ts, tw, s0.gapBelow, epsPair(EPS_STEEL, EPS_PE)) * A_FOOT * (tw - ts) +
      sideH(tw, ts, PAN.panM, (CAVITY.w - PAN.w) / 2, epsPair(EPS_STEEL, EPS_PE)) * sides * (tw - ts);
    tw -= (q * dt) / C_W;
    wallStep(q);
  }
  // The chili goes in; the water comes out.
  const wallStored = () => {
    let e = 0;
    for (let i = 0; i < NW; i++) e += wC[i]! * wT[i]!;
    return e;
  };
  const cellV = dx * dy * dz;
  const stored = () => {
    let e = wallStored();
    for (let p = 0; p < P; p++) {
      const t = T[p]!;
      for (let k = 0; k < nz; k++)
        for (let j = 0; j < ny; j++)
          for (let i = 0; i < nx; i++) {
            const c = t[at(i, j, k)]!;
            e += 4 * (enthalpy(c) * cellV + steelOf(i, j, k) * c);
          }
      e += C_LID * lid[p]!;
    }
    return e;
  };
  /** The steel (J/K) on a cell's outer faces: the pan's floor and walls (a quarter-pan's cell). */
  function steelOf(i: number, j: number, k: number): number {
    let a = 0;
    if (j === 0) a += dx * dz;
    if (i === nx - 1) a += dy * dz;
    if (k === nz - 1) a += dx * dy;
    return a * STEEL_C;
  }
  lost = 0;
  const at0 = stored();

  const epsSS = epsPair(EPS_STEEL, EPS_STEEL);
  const epsSP = epsPair(EPS_STEEL, EPS_PE);
  const bottomMean = new Float64Array(P);
  const topMean = new Float64Array(P);
  const sideMean = new Float64Array(P);
  const endMean = new Float64Array(P);
  let hotMin: number | null = null;
  let mid141: number | null = null;
  let prevCold = startC;
  let prevMid = startC;
  let step = 0;
  const kOf = Array.from({ length: P }, () => new Float64Array(NC));
  const capOf = Array.from({ length: P }, () => new Float64Array(NC));
  const gBot = Array.from({ length: P }, () => new Float64Array(nx * nz));
  const gBotSum = new Float64Array(P);
  const gBotT = new Float64Array(P);
  const hTopOf = new Float64Array(P);
  const lidNew = new Float64Array(P);
  for (;;) {
    const m = (step * dt) / 60;
    if (step % perSample === 0) {
      record(k0++, m);
      const kk = k0 - 1;
      const cold = coldest[kk]!;
      let mid = Infinity;
      for (let p = 0; p < P; p++) mid = Math.min(mid, centre[p]![kk]!);
      // Where the line was crossed between this sample and the last; at packing, if packed under it.
      const cross = (line: number, prev: number, now: number) => (m === 0 || !(prev > now) ? m : m - (SAMPLE_MIN * (line - now)) / (prev - now));
      if (hotMin === null && cold < HOT_LINE_C) hotMin = cross(HOT_LINE_C, prevCold, cold);
      if (mid141 === null && mid < fToC(141)) mid141 = cross(fToC(141), prevMid, mid);
      prevCold = cold;
      prevMid = mid;
      const done = m >= minMin && hotMin !== null && (!o.to141 || mid141 !== null);
      if (done || m >= maxMin || k0 >= n) break;
    }
    const ts = wT[0]!;
    // Each pan's faces' mean temperatures.
    for (let p = 0; p < P; p++) {
      const t = T[p]!;
      let b = 0;
      let u = 0;
      let sd = 0;
      let en = 0;
      for (let k = 0; k < nz; k++)
        for (let i = 0; i < nx; i++) {
          b += t[at(i, 0, k)]!;
          u += t[at(i, ny - 1, k)]!;
        }
      for (let k = 0; k < nz; k++) for (let j = 0; j < ny; j++) sd += t[at(nx - 1, j, k)]!;
      for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) en += t[at(i, j, nz - 1)]!;
      bottomMean[p] = b / (nx * nz);
      topMean[p] = u / (nx * nz);
      sideMean[p] = sd / (ny * nz);
      endMean[p] = en / (nx * ny);
    }
    let toWall = 0;
    // The conductances this step, from last step's temperatures. A pan's floor faces the lid below it (or the cavity's
    // floor) across their air layer: one conductance per floor cell, used from both sides, so what one gives the other
    // takes.
    for (let p = 0; p < P; p++) {
      const t = T[p]!;
      const sl = slots[p]!;
      const kp = kOf[p]!;
      const cp = capOf[p]!;
      for (let c = 0; c < NC; c++) {
        kp[c] = tab(TAB_K, t[c]!);
        cp[c] = tab(TAB_RC, t[c]!) * cellV;
      }
      const tBelow = sl.below >= 0 ? lid[sl.below]! : ts;
      const hB = layerH(tBelow, bottomMean[p]!, sl.gapBelow, sl.below >= 0 ? epsSS : epsSP) * (A_FLOOR / A_FOOT);
      let g = 0;
      let gt = 0;
      const gb = gBot[p]!;
      for (let k = 0; k < nz; k++)
        for (let i = 0; i < nx; i++) {
          const c = at(i, 0, k);
          const u = (dx * dz) / (1 / hB + dy / 2 / kp[c]!);
          gb[k * nx + i] = u;
          g += 4 * u;
          gt += 4 * u * t[c]!;
        }
      gBotSum[p] = g;
      gBotT[p] = gt;
    }
    // The lids, each between its own chili's top and what is above it: the next pan's floor or the cavity's roof.
    for (let p = 0; p < P; p++) {
      const t = T[p]!;
      const sl = slots[p]!;
      const kp = kOf[p]!;
      const dTop = topMean[p]! - lid[p]!;
      const hTop = Math.max(0.05, Math.abs(dTop) > 1e-6 ? across(topMean[p]!, lid[p]!, HEADSPACE) / dTop : across(lid[p]! + 0.01, lid[p]!, HEADSPACE) / 0.01);
      hTopOf[p] = hTop;
      let gIn = 0;
      let gInT = 0;
      for (let k = 0; k < nz; k++)
        for (let i = 0; i < nx; i++) {
          const c = at(i, ny - 1, k);
          const u = (4 * dx * dz) / (1 / hTop + dy / 2 / kp[c]!);
          gIn += u;
          gInT += u * t[c]!;
        }
      let gUp: number;
      let gUpT: number;
      if (sl.above >= 0) {
        gUp = gBotSum[sl.above]!;
        gUpT = gBotT[sl.above]!;
      } else {
        gUp = layerH(lid[p]!, ts, sl.gapAbove, epsSP) * A_LID;
        gUpT = gUp * ts;
      }
      lidNew[p] = ((C_LID / dt) * lid[p]! + gInT + gUpT) / (C_LID / dt + gIn + gUp);
      if (sl.above < 0) toWall += gUp * (lidNew[p]! - ts);
    }
    // Conduction inside each pan, and its faces.
    for (let p = 0; p < P; p++) {
      const t = T[p]!;
      const sl = slots[p]!;
      const kp = kOf[p]!;
      const cp = capOf[p]!;
      const gb = gBot[p]!;
      const tBelow = sl.below >= 0 ? lidNew[sl.below]! : ts;
      // The long walls: across the air to the side walls, and up the steel to the rim on its rail.
      const hS = sideH(sideMean[p]!, ts, PAN.panM, (CAVITY.w - PAN.w) / 2, epsSP) + RIM_G / (2 * PAN.l * PAN.foodM);
      const hE = sideH(endMean[p]!, ts, PAN.panM, (CAVITY.d - PAN.l) / 2, epsSP);
      const hTop = hTopOf[p]!;
      const lp = lidNew[p]!;
      dT.fill(0);
      for (let k = 0; k < nz; k++)
        for (let j = 0; j < ny; j++)
          for (let i = 0; i < nx; i++) {
            const c = at(i, j, k);
            const tc = t[c]!;
            let q = 0;
            if (i + 1 < nx) {
              const d = at(i + 1, j, k);
              const g = (((kp[c]! + kp[d]!) / 2) * dy * dz) / dx;
              q -= g * (tc - t[d]!);
              dT[d] = dT[d]! + g * (tc - t[d]!);
            } else {
              const g = (dy * dz) / (1 / hS + dx / 2 / kp[c]!);
              q -= g * (tc - ts);
              toWall += 4 * g * (tc - ts);
            }
            if (j + 1 < ny) {
              const d = at(i, j + 1, k);
              const g = (((kp[c]! + kp[d]!) / 2) * dx * dz) / dy;
              q -= g * (tc - t[d]!);
              dT[d] = dT[d]! + g * (tc - t[d]!);
            } else {
              const g = (dx * dz) / (1 / hTop + dy / 2 / kp[c]!);
              q -= g * (tc - lp);
            }
            if (j === 0) {
              const g = gb[k * nx + i]!;
              q -= g * (tc - tBelow);
              if (sl.below < 0) toWall += 4 * g * (tc - ts);
            }
            if (k + 1 < nz) {
              const d = at(i, j, k + 1);
              const g = (((kp[c]! + kp[d]!) / 2) * dx * dy) / dz;
              q -= g * (tc - t[d]!);
              dT[d] = dT[d]! + g * (tc - t[d]!);
            } else {
              const g = (dx * dy) / (1 / hE + dz / 2 / kp[c]!);
              q -= g * (tc - ts);
              toWall += 4 * g * (tc - ts);
            }
            dT[c] = dT[c]! + q;
          }
      for (let k = 0; k < nz; k++)
        for (let j = 0; j < ny; j++)
          for (let i = 0; i < nx; i++) {
            const c = at(i, j, k);
            t[c] = t[c]! + (dT[c]! * dt) / (cp[c]! + steelOf(i, j, k));
          }
    }
    lid.set(lidNew);
    wallStep(toWall);
    step++;
  }
  return { way, startC, grid: [nx, ny], endMin: (step * dt) / 60, centre, coldest, wall, water, slice, hotMin, mid141, energy: { at0, atEnd: stored(), lost } };
}

/** Every way's run for chili packed at `startC`. */
export function holdAll(startC: number, o: HoldOpts = {}): HoldRun[] {
  return HOLD_WAYS.map((w) => holdRun(w, startC, o));
}

/** A run's value at minute `m` after packing, straight between samples. */
export function sampleOf(a: Float32Array, m: number): number {
  const x = Math.min(a.length - 1, Math.max(0, sampleAt(m)));
  const i = Math.min(a.length - 2, Math.floor(x));
  const t = x - i;
  const v0 = a[i]!;
  const v1 = a[i + 1]!;
  if (Number.isNaN(v1)) return v0;
  if (Number.isNaN(v0)) return v1;
  return v0 + (v1 - v0) * t;
}
