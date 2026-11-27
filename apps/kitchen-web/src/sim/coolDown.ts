/**
 * Cooked the day before, cooled in time (chapter 06, second figure): a batch of chili off the stove, cooled four ways
 * in the walk-in and the walk-in freezer, against the Food Code's two limits (135 °F to 70 °F within 2 hours, and to
 * 41 °F within 6 hours in all, FDA Food Code 3-501.14).
 *
 * The food conducts its heat out from its middle: a slab across its depth in a pan (a full-size pan is 12 by 20 in;
 * at 2 or 3 in deep its middle is as far from the top and bottom as anything, the sides 6 in away count for nothing
 * there), and a cylinder, across and up, in a pot. Its conductivity, density and heat capacity at each temperature
 * come from its make-up (USDA's chili con carne with beans) by Choi and Okos's equations. Out of each face, every
 * minute:
 *   - the air's draw, h (T_s − T_air), one coefficient for the whole box: the only number fitted, to Roberts et al.'s
 *     four runs in a school's walk-in cooler and freezer loaded as they are on a working day;
 *   - radiation to the box's walls at the air's temperature (open food 0.94, stainless 0.25, film 0.9);
 *   - off open food, steam: the water vapour's pull from the wet surface into the box air, by the same draw (the
 *     Chilton-Colburn analogy), and the heat it takes;
 *   - under a cover, the air between the food and the cover (Hollands's layer heated from below), radiation across
 *     it, and steam carried up to condense on the cover, which then loses it all outside.
 * Ice forms on the freezer's pans only in their last minutes there and only in the outer millimetres: left out, it
 * changes their time to 41 °F by under 1%.
 *
 * Against Roberts's four runs and six more of Olds's (an empty walk-in, its own coefficient fitted the same way), the
 * times to 70 °F and to 41 °F land within a fifth of all twenty measurements (sim/coolDown.test.ts).
 */
import { POT_40QT } from "./potBoil";

/** Mass fractions of a food: water, protein, fat, carbohydrate other than fibre, fibre, ash. */
export interface Makeup {
  water: number;
  protein: number;
  fat: number;
  carb: number;
  fiber: number;
  ash: number;
}
/** USDA FoodData Central, SR Legacy 172097, chili con carne with beans: 75.7 g water, 5.8 protein, 3.47 fat, 13.1
 *  carbohydrate (3.3 of it fibre) and 1.94 ash in 100 g. */
export const CHILI: Makeup = { water: 0.757, protein: 0.058, fat: 0.0347, carb: 0.131 - 0.033, fiber: 0.033, ash: 0.0194 };

/** Conductivity (W/m·K), density (kg/m³) and specific heat (J/kg·K) of a food at `c` °C: Choi and Okos (1986), each
 *  part's own by temperature, heat capacity and density by mass, conductivity by volume (ASHRAE Handbook,
 *  Refrigeration, "Thermal Properties of Foods"). */
export function foodProps(f: Makeup, c: number): { k: number; rho: number; cp: number } {
  const t = c;
  const t2 = t * t;
  const parts: [number, number, number, number][] = [
    [f.water, 0.57109 + 1.7625e-3 * t - 6.7036e-6 * t2, 997.18 + 3.1439e-3 * t - 3.7574e-3 * t2, 4176.2 - 9.0864e-2 * t + 5.4731e-3 * t2],
    [f.protein, 0.17881 + 1.1958e-3 * t - 2.7178e-6 * t2, 1329.9 - 0.5184 * t, 2008.2 + 1.2089 * t - 1.3129e-3 * t2],
    [f.fat, 0.18071 - 2.7604e-4 * t - 1.7749e-7 * t2, 925.59 - 0.41757 * t, 1984.2 + 1.4733 * t - 4.8008e-3 * t2],
    [f.carb, 0.20141 + 1.3874e-3 * t - 4.3312e-6 * t2, 1599.1 - 0.31046 * t, 1548.8 + 1.9625 * t - 5.9399e-3 * t2],
    [f.fiber, 0.18331 + 1.2497e-3 * t - 3.1683e-6 * t2, 1311.5 - 0.36589 * t, 1845.9 + 1.8306 * t - 4.6509e-3 * t2],
    [f.ash, 0.32962 + 1.4011e-3 * t - 2.9069e-6 * t2, 2423.8 - 0.28063 * t, 1092.6 + 1.8896 * t - 3.6817e-3 * t2],
  ];
  let vol = 0;
  let cp = 0;
  for (const [x, , r, c0] of parts) {
    vol += x / r;
    cp += x * c0;
  }
  let k = 0;
  for (const [x, kk, r] of parts) k += (x / r / vol) * kk;
  return { k, rho: 1 / vol, cp };
}

const SIGMA = 5.670374419e-8;
/** °F to °C. */
export const fToC = (f: number) => (f - 32) / 1.8;
export const cToF = (c: number) => c * 1.8 + 32;
/** The Food Code's lines: hot holding, the end of the first stage, cold holding (°C). */
export const HOT_C = fToC(135);
export const WARM_C = fToC(70);
export const COLD_C = fToC(41);
/** Minutes allowed from 135 °F: to 70 °F, and to 41 °F in all. */
export const STAGE1_MIN = 120;
export const TOTAL_MIN = 360;
/** Chili goes into its pans at 190 °F, as it came off the stove in Olds and Sneed's runs. */
export const START_C = fToC(190);

/** Water vapour's saturation pressure (Pa) at `c` °C (Magnus, Alduchov and Eskridge), its density at relative
 *  humidity `rh` (kg/m³), and the heat it takes to evaporate (J/kg). */
export const psat = (c: number) => 610.94 * Math.exp((17.625 * c) / (c + 243.04));
const vapour = (c: number, rh: number) => (rh * psat(c)) / (461.5 * (c + 273.15));
export const latent = (c: number) => (2501 - 2.37 * c) * 1000;
/** Water activity at the chili's surface. */
const AW = 0.98;

/** A box's air: its temperature, humidity, and its draw on what is cooling in it (W/m²·K). */
export interface Box {
  c: number;
  rh: number;
  h: number;
}
/** Roberts et al. (2013), Food Protection Trends 33(1): a school's walk-in cooler at 39 °F and freezer at 0 °F, each
 *  80% loaded as the schools they surveyed load them. The draw is fitted to their four runs. */
export const H_FULL = 6.5;
export const WALK_IN: Box = { c: 3.9, rh: 0.85, h: H_FULL };
export const FREEZER: Box = { c: -17.8, rh: 0.85, h: H_FULL };
/** Olds (2019), J. Foodservice Mgmt & Educ. 13(1): a university walk-in at 34.6 °F, empty but for the test pans,
 *  spaced apart on wire shelving. Its draw, fitted to its six runs. */
export const H_EMPTY = 9.5;
export const EMPTY_WALK_IN: Box = { c: fToC(34.6), rh: 0.85, h: H_EMPTY };

/** What a face of the food looks out through. */
export type Face =
  | { kind: "open" }
  | { kind: "wall"; eps: number }
  /** A cover `gap` above the food: film (0.9 both sides), a stainless lid (0.25 outside) or foil (0.05 outside); its
   *  inside, beaded with condensate, 0.9. */
  | { kind: "cover"; gap: number; epsOut: number };
export const EPS_FOOD = 0.94;
export const EPS_STEEL = 0.25;
export const EPS_ALUMINIUM = 0.1;

/** Steam's pull (kg/m²·s) off open food at `ts` into the box. */
function steamFlux(ts: number, box: Box): number {
  const hm = box.h / (1.2 * 1006 * 0.85 ** (2 / 3));
  return Math.max(0, hm * (AW * psat(ts) / (461.5 * (ts + 273.15)) - vapour(box.c, box.rh)));
}

/** Heat (W/m²) out of a face at surface temperature `ts`, and how much of it went as steam off open food; under a
 *  cover, what crosses the air layer to the cover, at `cover` °C. */
function faceFlux(ts: number, f: Face, box: Box, cover: number): { q: number; steam: number } {
  const tk = ts + 273.15;
  const ak = box.c + 273.15;
  if (f.kind === "open") {
    const steam = steamFlux(ts, box) * latent(ts);
    return { q: box.h * (ts - box.c) + EPS_FOOD * SIGMA * (tk ** 4 - ak ** 4) + steam, steam };
  }
  if (f.kind === "wall") return { q: box.h * (ts - box.c) + f.eps * SIGMA * (tk ** 4 - ak ** 4), steam: 0 };
  return { q: across(ts, cover, f.gap), steam: 0 };
}
/** Heat (W/m²) across an air layer `gap` m deep from food at `ts` up to a cover at `tc`: the layer's convection
 *  (Hollands, heated from below), radiation, and steam carried up that condenses on the cover. */
export function across(ts: number, tc: number, gap: number): number {
  const fm = (ts + tc) / 2 + 273.15;
  const nu = 1.516e-5 * (fm / 293.15) ** 1.75;
  const ka = 0.02587 * (fm / 293.15) ** 0.85;
  const dv = 2.5e-5 * (fm / 293.15) ** 1.75;
  const ra = Math.max(0, ((9.80665 / fm) * (ts - tc) * gap ** 3) / ((nu * nu) / 0.71));
  const nus = 1 + 1.44 * Math.max(0, 1 - 1708 / Math.max(ra, 1)) + Math.max(0, Math.cbrt(ra / 5830) - 1);
  const rad = (SIGMA * ((ts + 273.15) ** 4 - (tc + 273.15) ** 4)) / (1 / EPS_FOOD + 1 / 0.9 - 1);
  const cond = Math.max(0, ((nus * dv) / gap) * ((AW * psat(ts)) / (461.5 * (ts + 273.15)) - vapour(tc, 1))) * latent(ts);
  return ((nus * ka) / gap) * (ts - tc) + rad + cond;
}
/** The cover's temperature (one piece: a lid or a film): what reaches it from the food's top, surface temperatures
 *  `ts` over areas `a`, equals what it loses outside. */
function coverC(ts: ArrayLike<number>, a: ArrayLike<number>, f: Face & { kind: "cover" }, box: Box): number {
  const ak = box.c + 273.15;
  let area = 0;
  let hi = box.c;
  for (let i = 0; i < ts.length; i++) {
    area += a[i]!;
    hi = Math.max(hi, ts[i]!);
  }
  let lo = box.c;
  for (let it = 0; it < 40; it++) {
    const m = (lo + hi) / 2;
    let inW = 0;
    for (let i = 0; i < ts.length; i++) inW += across(ts[i]!, m, f.gap) * a[i]!;
    const outW = (box.h * (m - box.c) + f.epsOut * SIGMA * ((m + 273.15) ** 4 - ak ** 4)) * area;
    if (inW > outW) lo = m;
    else hi = m;
  }
  return (lo + hi) / 2;
}

/** The surface temperature of a face whose nearest cell, at `tc`, lies `dist` inside it (conductivity `k`): heat
 *  conducted to the face equals heat out of it. Newton's method from `guess` (last step's). */
function surfaceC(tc: number, dist: number, k: number, f: Face, box: Box, cover: number, guess: number): number {
  let ts = guess;
  for (let i = 0; i < 12; i++) {
    const g = (k * (tc - ts)) / dist - faceFlux(ts, f, box, cover).q;
    const d = 1e-3;
    const g2 = (k * (tc - ts - d)) / dist - faceFlux(ts + d, f, box, cover).q;
    const step = (g * d) / (g2 - g);
    ts -= step;
    if (Math.abs(step) < 1e-5) break;
  }
  return ts;
}

/** The food's properties by temperature, every quarter degree from −30 to 100 °C. */
const TAB0 = -30;
const TAB_STEP = 0.25;
interface Table {
  k: Float64Array;
  rc: Float64Array;
  rho: Float64Array;
}
const tables = new Map<Makeup, Table>();
function table(f: Makeup): Table {
  const hit = tables.get(f);
  if (hit) return hit;
  const n = Math.round(130 / TAB_STEP) + 1;
  const t: Table = { k: new Float64Array(n), rc: new Float64Array(n), rho: new Float64Array(n) };
  for (let i = 0; i < n; i++) {
    const p = foodProps(f, TAB0 + i * TAB_STEP);
    t.k[i] = p.k;
    t.rc[i] = p.rho * p.cp;
    t.rho[i] = p.rho;
  }
  tables.set(f, t);
  return t;
}
const lookup = (a: Float64Array, c: number) => {
  const x = Math.min(a.length - 1.000001, Math.max(0, (c - TAB0) / TAB_STEP));
  const i = Math.floor(x);
  return a[i]! + (x - i) * (a[i + 1]! - a[i]!);
};

/** What a run is cooled in, from which minute: the freezer, then the walk-in once its middle reaches 41 °F. */
export interface Plan {
  box: Box;
  /** Where it goes once its middle is at 41 °F, if anywhere; or at minute `at`, if given (a pan on a rack moved when
   *  the rack's other pans are cold). */
  then?: Box;
  at?: number;
}
/** Every minute of a run that its samples and the stage need. */
export const SAMPLE_MIN = 2;
export const SNAP_MIN = 10;

export interface CoolRun {
  /** The middle's temperature every SAMPLE_MIN from going in (°C). */
  centre: Float32Array;
  /** The steam off its open food every SAMPLE_MIN, W/m² (0 under a cover). */
  steam: Float32Array;
  /** The whole field every SNAP_MIN: a pan's cells bottom to top; a pot's rows bottom to top, each from the axis
   *  out. */
  snaps: Float32Array[];
  shape: { kind: "slab"; n: number; depth: number } | { kind: "pot"; nr: number; nz: number; r: number; fill: number };
  /** Minutes from going in to its middle reaching 135 °F, 70 °F and 41 °F (null: not within the run). */
  t135: number | null;
  t70: number | null;
  t41: number | null;
  /** The minute it moved from the freezer to the walk-in, if it did. */
  movedMin: number | null;
  /** The run's length (minutes). */
  endMin: number;
  /** Heat out of the food (J per m² of top, slab; J, pot), and what that leaves the food's heat content at; and the
   *  water it lost to steam (kg per m² of top, slab; kg, pot). For checking, and for the words. */
  heatOut: number;
  heatDrop: number;
  /** Of heatOut, the heat the steam took (same units). */
  steamJ: number;
  steamKg: number;
}

export interface RunOpts {
  food?: Makeup;
  startC?: number;
  /** Run until the middle is at 41 °F, this many minutes have passed, and pastHotMin more since it fell to 135 °F;
   *  or for maxMin at most. */
  minMin?: number;
  pastHotMin?: number;
  maxMin?: number;
  /** Test hook: hold the properties at those of this temperature. */
  constAt?: number;
}

/** A pan's food: `depth` m deep, its top open or under a cover, its bottom on the pan's floor. */
export function coolSlab(depth: number, top: Face, bottom: Face, plan: Plan, o: RunOpts = {}): CoolRun {
  const food = o.food ?? CHILI;
  const tab = table(food);
  const n = Math.max(12, Math.round(depth / 0.002));
  const dx = depth / n;
  const T = new Float64Array(n).fill(o.startC ?? START_C);
  const kAt = (c: number) => lookup(tab.k, o.constAt ?? c);
  const rcAt = (c: number) => lookup(tab.rc, o.constAt ?? c);
  let alphaMax = 0;
  for (let i = 0; i < tab.k.length; i++) alphaMax = Math.max(alphaMax, tab.k[i]! / tab.rc[i]!);
  // Explicit steps, well inside the stability limit, a whole number to each sample.
  const per = Math.ceil((SAMPLE_MIN * 60) / Math.min(20, (0.4 * dx * dx) / alphaMax));
  const dt = (SAMPLE_MIN * 60) / per;
  const mid = () => (n % 2 ? T[(n - 1) / 2]! : (T[n / 2 - 1]! + T[n / 2]!) / 2);
  const q = new Float64Array(n + 1);
  const run = newRun({ kind: "slab", n, depth });
  const heat0 = slabHeat(T, dx, tab, o.constAt);
  let box = plan.box;
  let steamNow = 0;
  const record = (min: number) => {
    run.centre.push(mid());
    run.steam.push(steamNow);
    if (min % SNAP_MIN === 0) run.snaps.push(Float32Array.from(T));
  };
  record(0);
  const maxMin = o.maxMin ?? 72 * 60;
  let tb = T[0]!;
  let tt = T[n - 1]!;
  let cover = (T[n - 1]! + box.c) / 2;
  const one = [0];
  const unit = [1];
  for (let min = 0; min < maxMin; ) {
    for (let s = 0; s < per; s++) {
      tb = surfaceC(T[0]!, dx / 2, kAt(T[0]!), bottom, box, cover, tb);
      tt = surfaceC(T[n - 1]!, dx / 2, kAt(T[n - 1]!), top, box, cover, tt);
      const fb = faceFlux(tb, bottom, box, cover);
      const ft = faceFlux(tt, top, box, cover);
      if (top.kind === "cover") {
        one[0] = tt;
        cover = coverC(one, unit, top, box);
      }
      q[0] = -fb.q;
      q[n] = ft.q;
      for (let i = 1; i < n; i++) q[i] = ((kAt(T[i - 1]!) + kAt(T[i]!)) / 2) * ((T[i - 1]! - T[i]!) / dx);
      for (let i = 0; i < n; i++) T[i] = T[i]! + ((q[i]! - q[i + 1]!) * dt) / (rcAt(T[i]!) * dx);
      run.heatOut += (fb.q + ft.q) * dt;
      run.steamKg += ((ft.steam + fb.steam) / latent(tt)) * dt;
      run.steamJ += (ft.steam + fb.steam) * dt;
      steamNow = ft.steam;
    }
    min += SAMPLE_MIN;
    record(min);
    const m = mid();
    marks(run, m, min);
    if (movesNow(run, plan, min)) {
      run.movedMin = min;
      box = plan.then!;
    }
    if (done(run, min, o)) break;
  }
  run.endMin = (run.centre.length - 1) * SAMPLE_MIN;
  run.heatDrop = heat0 - slabHeat(T, dx, tab, o.constAt);
  return finish(run);
}

/** A pot's food, `r` m across its middle and `fill` m deep, under its lid or open, its floor and wall stainless. */
export function coolPot(r: number, fill: number, top: Face, bottom: Face, side: Face, plan: Plan, o: RunOpts = {}): CoolRun {
  const food = o.food ?? CHILI;
  const tab = table(food);
  const nr = 20;
  const dr = r / nr;
  const nz = Math.max(8, Math.round(fill / dr));
  const dz = fill / nz;
  const T = new Float64Array(nr * nz).fill(o.startC ?? START_C);
  const kAt = (c: number) => lookup(tab.k, o.constAt ?? c);
  const rcAt = (c: number) => lookup(tab.rc, o.constAt ?? c);
  let alphaMax = 0;
  for (let i = 0; i < tab.k.length; i++) alphaMax = Math.max(alphaMax, tab.k[i]! / tab.rc[i]!);
  const per = Math.ceil((SAMPLE_MIN * 60) / Math.min(60, (0.2 * Math.min(dr, dz) ** 2) / alphaMax));
  const dt = (SAMPLE_MIN * 60) / per;
  const at = (i: number, j: number) => j * nr + i;
  // The middle: on the axis, halfway up.
  const mid = () => (nz % 2 ? T[at(0, (nz - 1) / 2)]! : (T[at(0, nz / 2 - 1)]! + T[at(0, nz / 2)]!) / 2);
  const ring = (i: number) => Math.PI * (((i + 1) * dr) ** 2 - (i * dr) ** 2);
  const dT = new Float64Array(nr * nz);
  const run = newRun({ kind: "pot", nr, nz, r, fill });
  const potHeat = () => {
    let e = 0;
    for (let j = 0; j < nz; j++) for (let i = 0; i < nr; i++) e += cellHeat(T[at(i, j)]!, tab, o.constAt) * ring(i) * dz;
    return e;
  };
  const heat0 = potHeat();
  let box = plan.box;
  let steamNow = 0;
  const record = (min: number) => {
    run.centre.push(mid());
    run.steam.push(steamNow);
    if (min % SNAP_MIN === 0) run.snaps.push(Float32Array.from(T));
  };
  record(0);
  const maxMin = o.maxMin ?? 72 * 60;
  // Last step's surface temperatures: the top and the floor across, the wall up; and the lid's.
  const tTop = Float64Array.from({ length: nr }, (_, i) => T[at(i, nz - 1)]!);
  const tBot = Float64Array.from({ length: nr }, (_, i) => T[at(i, 0)]!);
  const tSide = Float64Array.from({ length: nz }, (_, j) => T[at(nr - 1, j)]!);
  const rings = Float64Array.from({ length: nr }, (_, i) => ring(i));
  let cover = (START_C + box.c) / 2;
  for (let min = 0; min < maxMin; ) {
    for (let s = 0; s < per; s++) {
      dT.fill(0);
      let out = 0;
      let steamW = 0;
      let steamKg = 0;
      for (let j = 0; j < nz; j++)
        for (let i = 0; i < nr; i++) {
          const c = at(i, j);
          const tc = T[c]!;
          const kc = kAt(tc);
          let qin = 0;
          // Across, out from the axis.
          if (i < nr - 1) {
            const k = (kc + kAt(T[at(i + 1, j)]!)) / 2;
            const flow = ((k * (T[at(i + 1, j)]! - tc)) / dr) * 2 * Math.PI * (i + 1) * dr * dz;
            qin += flow;
            dT[at(i + 1, j)] = dT[at(i + 1, j)]! - flow;
          } else {
            tSide[j] = surfaceC(tc, dr / 2, kc, side, box, cover, tSide[j]!);
            const f = faceFlux(tSide[j]!, side, box, cover);
            const a = 2 * Math.PI * r * dz;
            qin -= f.q * a;
            out += f.q * a;
          }
          // Up.
          if (j < nz - 1) {
            const k = (kc + kAt(T[at(i, j + 1)]!)) / 2;
            const flow = ((k * (T[at(i, j + 1)]! - tc)) / dz) * ring(i);
            qin += flow;
            dT[at(i, j + 1)] = dT[at(i, j + 1)]! - flow;
          } else {
            const ts = (tTop[i] = surfaceC(tc, dz / 2, kc, top, box, cover, tTop[i]!));
            const f = faceFlux(ts, top, box, cover);
            qin -= f.q * ring(i);
            out += f.q * ring(i);
            steamW += f.steam * ring(i);
            steamKg += (f.steam / latent(ts)) * ring(i);
          }
          if (j === 0) {
            tBot[i] = surfaceC(tc, dz / 2, kc, bottom, box, cover, tBot[i]!);
            const f = faceFlux(tBot[i]!, bottom, box, cover);
            qin -= f.q * ring(i);
            out += f.q * ring(i);
          }
          dT[c] = dT[c]! + qin;
        }
      for (let j = 0; j < nz; j++) for (let i = 0; i < nr; i++) T[at(i, j)] = T[at(i, j)]! + (dT[at(i, j)]! * dt) / (rcAt(T[at(i, j)]!) * rings[i]! * dz);
      if (top.kind === "cover") cover = coverC(tTop, rings, top, box);
      run.heatOut += out * dt;
      run.steamKg += steamKg * dt;
      run.steamJ += steamW * dt;
      steamNow = steamW / (Math.PI * r * r);
    }
    min += SAMPLE_MIN;
    record(min);
    marks(run, mid(), min);
    if (movesNow(run, plan, min)) {
      run.movedMin = min;
      box = plan.then!;
    }
    if (done(run, min, o)) break;
  }
  run.endMin = (run.centre.length - 1) * SAMPLE_MIN;
  run.heatDrop = heat0 - potHeat();
  return finish(run);
}

/** A cell's heat content per m³ above 0 °C: ∫ ρc dT. */
function cellHeat(c: number, tab: Table, constAt?: number): number {
  if (constAt !== undefined) return lookup(tab.rc, constAt) * c;
  // Trapezoids over the table, from 0 °C.
  const steps = Math.max(1, Math.ceil(Math.abs(c) / 0.5));
  const h = c / steps;
  let e = 0;
  for (let s = 0; s < steps; s++) e += ((lookup(tab.rc, s * h) + lookup(tab.rc, (s + 1) * h)) / 2) * h;
  return e;
}
function slabHeat(T: Float64Array, dx: number, tab: Table, constAt?: number): number {
  let e = 0;
  for (const c of T) e += cellHeat(c, tab, constAt) * dx;
  return e;
}

interface Building extends Omit<CoolRun, "centre" | "steam"> {
  centre: number[];
  steam: number[];
}
function newRun(shape: CoolRun["shape"]): Building {
  return { centre: [], steam: [], snaps: [], shape, t135: null, t70: null, t41: null, movedMin: null, endMin: 0, heatOut: 0, heatDrop: 0, steamKg: 0, steamJ: 0 };
}
/** Whether a run has gone far enough (RunOpts). */
function done(run: Building, min: number, o: RunOpts): boolean {
  return run.t41 !== null && min >= (o.minMin ?? 0) && (o.pastHotMin === undefined || (run.t135 !== null && min >= run.t135 + o.pastHotMin));
}
/** Whether a run moves out of its first box at this minute. */
function movesNow(run: Building, plan: Plan, min: number): boolean {
  return !!plan.then && run.movedMin === null && (plan.at !== undefined ? min >= plan.at : run.t41 !== null);
}
/** The minute (to a fraction, straight between samples) the middle crossed each line. */
function marks(run: Building, m: number, min: number): void {
  const prev = run.centre[run.centre.length - 2] ?? m;
  const cross = (line: number) => min - SAMPLE_MIN * ((line - m) / (prev - m || 1));
  if (run.t135 === null && m <= HOT_C) run.t135 = prev > HOT_C ? cross(HOT_C) : min - SAMPLE_MIN;
  if (run.t70 === null && m <= WARM_C) run.t70 = cross(WARM_C);
  if (run.t41 === null && m <= COLD_C) run.t41 = cross(COLD_C);
}
function finish(b: Building): CoolRun {
  return { ...b, centre: Float32Array.from(b.centre), steam: Float32Array.from(b.steam) };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* The four ways on the page.                                                                                     */
/* ------------------------------------------------------------------------------------------------------------ */

const IN = 0.0254;
/** Full-size hotel pans: their depth, and the litres of food at the depth they are filled to (the drawn pan's
 *  tapered inside, models/smallwares hotelPan, measured at that depth). */
export const PAN_2 = { panM: 0.065, foodM: 2 * IN, litres: 7.0 };
export const PAN_3 = { panM: 0.1, foodM: 3 * IN, litres: 10.2 };

export type WayKey = "pot" | "deep" | "shallow" | "freezer";
export interface Way {
  key: WayKey;
  name: string;
  /** For the stage's labels: few words. */
  short: string;
}
export const WAYS: readonly Way[] = [
  { key: "pot", name: "Left in its pot", short: "In its pot" },
  { key: "deep", name: "3-inch pans", short: "3 in deep" },
  { key: "shallow", name: "2-inch pans", short: "2 in deep" },
  { key: "freezer", name: "2-inch pans, freezer first", short: "Freezer first" },
];
/** The batches offered (litres), as the pot holds them. */
export const MIN_LITRES = 8;
export const MAX_LITRES = 32;

export interface WayRun {
  way: Way;
  run: CoolRun;
  /** Pans it takes (0 for the pot). */
  pans: number;
  /** The last pan's run when it holds less than the rest (coolRest), else null. */
  rest: CoolRun | null;
  /** Minutes from 135 °F to 70 °F and to 41 °F, and whether each is in time. */
  stage1: number | null;
  total: number | null;
  ok1: boolean;
  ok2: boolean;
}

/** The pot's food: its depth for `litres` in the 40-quart pot, and the room over it under the lid. */
export function potFill(litres: number): { r: number; fill: number; gap: number } {
  const r = POT_40QT.d / 2;
  const fill = litres / 1000 / (Math.PI * r * r);
  return { r, fill, gap: Math.max(0.01, POT_40QT.depth - fill) };
}

/** The pans a way takes for `litres` (none for the pot): each filled to its depth, the last with what is left. */
export function pansFor(key: WayKey, litres: number): number {
  return key === "pot" ? 0 : Math.ceil(litres / (key === "deep" ? PAN_3 : PAN_2).litres - 1e-9);
}
/** The litres in a way's last pan when it holds less than the full pans before it (0 when every pan is full, when
 *  there is only the one pan, and for the pot). */
export function restLitres(key: WayKey, litres: number): number {
  const n = pansFor(key, litres);
  if (key === "pot" || n < 2) return 0;
  const cap = (key === "deep" ? PAN_3 : PAN_2).litres;
  const r = litres - (n - 1) * cap;
  return r < cap - 1e-6 ? r : 0;
}
/** How deep a pan way's first pan is filled (m): to the pan's depth, or less when the whole batch is less than a pan. */
export function panDepth(key: Exclude<WayKey, "pot">, litres: number): number {
  const pan = key === "deep" ? PAN_3 : PAN_2;
  return pan.foodM * Math.min(1, litres / pan.litres);
}
/** The key a batch's runs go by. */
export const coolKey = (litres: number, covered: boolean) => `${litres}|${covered ? "c" : "o"}`;

const STEEL: Face = { kind: "wall", eps: EPS_STEEL };
/** A pan way's run for `depthM` of chili in its pan, open or under film; freezer first, it moves at `movedAt` if given
 *  (with its rack), else when its own middle is at 41 °F. */
function panRun(key: Exclude<WayKey, "pot">, depthM: number, covered: boolean, o: RunOpts, movedAt?: number): CoolRun {
  const pan = key === "deep" ? PAN_3 : PAN_2;
  const top: Face = covered ? { kind: "cover", gap: pan.panM - depthM, epsOut: 0.9 } : { kind: "open" };
  return coolSlab(depthM, top, STEEL, key === "freezer" ? { box: FREEZER, then: WALK_IN, at: movedAt } : { box: WALK_IN }, o);
}
/**
 * The last pan's run when it holds less than the rest (restLitres): shallower, so it cools faster. It stands on the
 * same rack, so freezer first it moves when the full pans do (`movedMin`, theirs; never, if null). Null when there is
 * no such pan. Only the stage's colours use it: the full pans are the slowest, and the verdicts are theirs.
 */
export function coolRest(key: WayKey, litres: number, covered: boolean, movedMin: number | null, o: RunOpts = {}): CoolRun | null {
  const rest = restLitres(key, litres);
  if (key === "pot" || rest <= 0) return null;
  const pan = key === "deep" ? PAN_3 : PAN_2;
  return panRun(key, Math.max(0.003, (pan.foodM * rest) / pan.litres), covered, o, movedMin ?? Infinity);
}

/** One way's run, for `litres` of chili, open or covered. */
export function coolWay(key: WayKey, litres: number, covered: boolean, o: RunOpts = {}): WayRun {
  const way = WAYS.find((w) => w.key === key)!;
  let run: CoolRun;
  let pans = 0;
  if (key === "pot") {
    const p = potFill(litres);
    const top: Face = covered ? { kind: "cover", gap: p.gap, epsOut: EPS_STEEL } : { kind: "open" };
    run = coolPot(p.r, p.fill, top, STEEL, STEEL, { box: WALK_IN }, o);
  } else {
    pans = pansFor(key, litres);
    run = panRun(key, panDepth(key, litres), covered, o);
  }
  return wayOf(way, run, pans, key === "pot" ? null : coolRest(key, litres, covered, run.movedMin, o));
}
/** A way's verdicts from its run. */
export function wayOf(way: Way, run: CoolRun, pans: number, rest: CoolRun | null): WayRun {
  const stage1 = run.t135 !== null && run.t70 !== null ? run.t70 - run.t135 : null;
  const total = run.t135 !== null && run.t41 !== null ? run.t41 - run.t135 : null;
  return { way, run, pans, rest, stage1, total, ok1: stage1 !== null && stage1 <= STAGE1_MIN, ok2: total !== null && total <= TOTAL_MIN };
}

/** The stage's day, and the chart's span from 135 °F: every way runs for at least both, and on to 41 °F if that
 *  comes within MAX_MIN of going in. */
export const DAY_MIN = 24 * 60;
export const CHART_MIN = 24 * 60;
export const MAX_MIN = 72 * 60;
export const PAGE_RUN: RunOpts = { minMin: DAY_MIN, pastHotMin: CHART_MIN, maxMin: MAX_MIN };
export function coolAll(litres: number, covered: boolean): WayRun[] {
  return WAYS.map((w) => coolWay(w.key, litres, covered, PAGE_RUN));
}

/** The temperature of a run's middle at minute `t` (straight between samples; the last past its end). */
export function centreAt(run: CoolRun, t: number): number {
  const x = Math.max(0, t / SAMPLE_MIN);
  const i = Math.min(run.centre.length - 1, Math.floor(x));
  const j = Math.min(run.centre.length - 1, i + 1);
  return run.centre[i]! + (x - i) * (run.centre[j]! - run.centre[i]!);
}

/* ------------------------------------------------------------------------------------------------------------ */
/* The measured runs it is held to.                                                                               */
/* ------------------------------------------------------------------------------------------------------------ */

export interface Measured {
  study: "Roberts 2013" | "Olds 2019";
  name: string;
  /** Hours:minutes, 135 °F to 70 °F and 135 °F to 41 °F, the studies' means. */
  to70: string;
  to41: string;
  run: (h: number) => CoolRun;
}
const box = (b: Box, h: number): Box => ({ ...b, h });
const wall = (eps: number): Face => ({ kind: "wall", eps });
/** Olds's 20-quart aluminium stockpot: 12 in across, 10½ in high, 3 gallons of chili in it. */
const OLDS_POT = { r: 6 * IN, fill: (3 * 3.785411784) / 1000 / (Math.PI * (6 * IN) ** 2), h: 10.5 * IN };
export const MEASURED: readonly Measured[] = [
  { study: "Roberts 2013", name: "2 in, walk-in", to70: "3:09", to41: "9:42", run: (h) => coolSlab(PAN_2.foodM, { kind: "open" }, wall(EPS_STEEL), { box: box(WALK_IN, h) }) },
  { study: "Roberts 2013", name: "3 in, walk-in", to70: "5:31", to41: "16:30", run: (h) => coolSlab(PAN_3.foodM, { kind: "open" }, wall(EPS_STEEL), { box: box(WALK_IN, h) }) },
  { study: "Roberts 2013", name: "2 in, freezer", to70: "1:47", to41: "3:10", run: (h) => coolSlab(PAN_2.foodM, { kind: "open" }, wall(EPS_STEEL), { box: box(FREEZER, h) }) },
  { study: "Roberts 2013", name: "3 in, freezer", to70: "2:38", to41: "4:50", run: (h) => coolSlab(PAN_3.foodM, { kind: "open" }, wall(EPS_STEEL), { box: box(FREEZER, h) }) },
  { study: "Olds 2019", name: "2 in, open", to70: "2:21", to41: "6:23", run: (h) => coolSlab(PAN_2.foodM, { kind: "open" }, wall(EPS_STEEL), { box: box(EMPTY_WALK_IN, h) }) },
  { study: "Olds 2019", name: "2 in, film", to70: "3:54", to41: "10:42", run: (h) => coolSlab(PAN_2.foodM, { kind: "cover", gap: 0.5 * IN, epsOut: 0.9 }, wall(EPS_STEEL), { box: box(EMPTY_WALK_IN, h) }) },
  { study: "Olds 2019", name: "3 in, open", to70: "3:47", to41: "10:15", run: (h) => coolSlab(PAN_3.foodM, { kind: "open" }, wall(EPS_STEEL), { box: box(EMPTY_WALK_IN, h) }) },
  { study: "Olds 2019", name: "3 in, film", to70: "4:54", to41: "13:40", run: (h) => coolSlab(PAN_3.foodM, { kind: "cover", gap: 1 * IN, epsOut: 0.9 }, wall(EPS_STEEL), { box: box(EMPTY_WALK_IN, h) }) },
  { study: "Olds 2019", name: "3 gal pot, open", to70: "5:53", to41: "16:12", run: (h) => coolPot(OLDS_POT.r, OLDS_POT.fill, { kind: "open" }, wall(EPS_ALUMINIUM), wall(EPS_ALUMINIUM), { box: box(EMPTY_WALK_IN, h) }) },
  { study: "Olds 2019", name: "3 gal pot, foil", to70: "7:21", to41: "20:17", run: (h) => coolPot(OLDS_POT.r, OLDS_POT.fill, { kind: "cover", gap: OLDS_POT.h - OLDS_POT.fill, epsOut: 0.05 }, wall(EPS_ALUMINIUM), wall(EPS_ALUMINIUM), { box: box(EMPTY_WALK_IN, h) }) },
];
export const hoursOf = (s: string) => {
  const [h, m] = s.split(":").map(Number);
  return h! + m! / 60;
};
