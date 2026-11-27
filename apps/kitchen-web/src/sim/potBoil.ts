/**
 * Water brought to the boil on a range (chapter 02, second figure): how long a pot takes, and what a lid and a second
 * pot on a second burner change.
 *
 * The water is one lump at one temperature (Berick measured it even to a tenth or two of a degree, top to bottom and
 * side to side, all the way to the boil). Heat in: the burner's input times the share of it that reaches the pot (an
 * open gas burner's 42%, as the stove kernel has it; the Food Service Technology Center found that share changes
 * little with the flame). Heat out, every second:
 *   - steam off the open water, at Berick's measured rate per square metre of water from 50 to 98 °C (a 25 cm
 *     stainless pot, 5 L, a room at 10 °C), straight between his points on a log scale and on at his last points'
 *     rate past them; less for a wider pot by the quarter power of its width (natural convection over a warm
 *     surface, laminar at these sizes: Sh ∝ Ra^1/4). A lid stops 92% of it (Berick, at 98 °C);
 *   - radiation from the open water (ε 0.94) or the lid (0.25), and from the pot's side as high as the water (0.25,
 *     shiny stainless): Berick's measured emissivities;
 *   - warm air rising off the side (Churchill and Chu's vertical plate) and off the top (a hot plate facing up).
 * The pot's stainless shell and its aluminium base heat with the water. The water boils at 100 °C; what it loses to
 * steam on the way leaves the pot.
 *
 * Every number here is computed from those: no time to boil is put in.
 */

/** Berick (Aprovecho Research Center), "Heat Losses in a Cook Pot at Constant Temperature", Figure 5: water lost to
 *  evaporation, grams an hour, from 0.0491 m² of open water in a 25 cm pot held at each temperature. */
export const BERICK_C = [50, 60, 70, 80, 90, 92, 94, 96, 98] as const;
export const BERICK_G_PER_H = [74, 146, 240, 410, 648, 706, 824, 954, 1128] as const;
const BERICK_AREA_M2 = 0.0491;
const BERICK_D_M = 0.25;
/** The share of the open water's evaporation a lid lets out (Berick: the lid was 92% effective at 98 °C). */
export const LID_LETS_OUT = 0.08;
export const EPS_WATER = 0.94;
export const EPS_SHINY_STEEL = 0.25;
const SIGMA = 5.670374419e-8;
const G = 9.80665;
const C_WATER = 4186;
const STEEL_J_PER_M3K = 7900 * 500;
const ALU_J_PER_M3K = 2700 * 900;

/** An open gas burner's input, 30 000 BTU/h (a commercial six-burner range's), and the share that reaches the pot. */
export const BURNER_W = 30000 * 0.29307107;
export const GAS_INTO_POT = 0.42;
/** Water from the cold tap, the kitchen's air, and where water boils (at sea level). */
export const TAP_C = 18;
export const ROOM_C = 25;
export const BOIL_C = 100;

/** Latent heat of vaporisation of water at `c` °C, J/kg (within a third of a percent from 0 to 100 °C). */
export const latentJPerKg = (c: number) => (2501 - 2.37 * c) * 1000;

/** Grams an hour off a square metre of open water at `c` °C, in Berick's 25 cm pot. */
export function evaporationGPerHM2(c: number): number {
  const n = BERICK_C.length;
  const i = c <= BERICK_C[0] ? 0 : c >= BERICK_C[n - 1]! ? n - 2 : BERICK_C.findIndex((x, k) => k + 1 < n && x <= c && c < BERICK_C[k + 1]!);
  const [a, b] = [BERICK_C[i]!, BERICK_C[i + 1]!];
  const [ga, gb] = [BERICK_G_PER_H[i]!, BERICK_G_PER_H[i + 1]!];
  return (ga * Math.exp(((c - a) / (b - a)) * Math.log(gb / ga))) / BERICK_AREA_M2;
}

/** Air at film temperature `c` °C: kinematic viscosity (m²/s), conductivity (W/m·K), Prandtl number. Power laws in
 *  absolute temperature through the table values, within 2% from 0 to 100 °C. */
function air(c: number): { nu: number; k: number; pr: number } {
  const r = (c + 273.15) / 293.15;
  return { nu: 1.516e-5 * r ** 1.75, k: 0.02587 * r ** 0.85, pr: 0.71 };
}
/** Rayleigh number for a surface at `ts` in air at `ta`, over length `l`. */
function rayleigh(ts: number, ta: number, l: number): { ra: number; k: number; pr: number } {
  const f = (ts + ta) / 2;
  const { nu, k, pr } = air(f);
  const beta = 1 / (f + 273.15);
  return { ra: (G * beta * Math.max(ts - ta, 1e-6) * l ** 3) / (nu * (nu / pr)), k, pr };
}
/** Natural convection off a vertical wall `h` high: Churchill and Chu. W/m²·K. */
export function hSide(ts: number, ta: number, h: number): number {
  if (h <= 0) return 0;
  const { ra, k, pr } = rayleigh(ts, ta, h);
  const nu = (0.825 + (0.387 * ra ** (1 / 6)) / (1 + (0.492 / pr) ** (9 / 16)) ** (8 / 27)) ** 2;
  return (nu * k) / h;
}
/** Natural convection off a hot disc `d` wide facing up (length area over perimeter, d/4): laminar 0.54 Ra^1/4 to
 *  1e7, turbulent 0.15 Ra^1/3 past it. W/m²·K. */
export function hTop(ts: number, ta: number, d: number): number {
  const l = d / 4;
  const { ra, k } = rayleigh(ts, ta, l);
  return ((ra < 1e7 ? 0.54 * ra ** 0.25 : 0.15 * ra ** (1 / 3)) * k) / l;
}

/** A pot: inside diameter and depth (m), and its shell's heat capacity (J/K), lid included when it has one on. */
export interface PotSpec {
  key: string;
  name: string;
  /** Inside diameter and the depth it holds water to (m). */
  d: number;
  depth: number;
  /** Litres it holds to the brim. */
  litres: number;
  /** The stainless shell and its aluminium base, and the lid (J/K). */
  shellJK: number;
  lidJK: number;
}
/**
 * A stockpot as the page draws it (models/smallwares `stockpot`, the archetype's outer width, depth and height):
 * a 1.4 mm drawn 304 shell, an 8.5 mm impact-bonded base (6 mm of aluminium under 2.5 mm of stainless), a 1 mm lid.
 */
function stockpot(key: string, name: string, w: number, h: number): PotSpec {
  const R = w / 2 - 0.0005;
  const t = 0.0014;
  const body = h - 0.012 - 0.024;
  const d = 2 * (R - t);
  const depth = body - 0.004 - (0.0085 + t);
  const base = Math.PI * R * R;
  return {
    key,
    name,
    d,
    depth,
    litres: Math.PI * (d / 2) ** 2 * depth * 1000,
    shellJK: 2 * Math.PI * R * body * t * STEEL_J_PER_M3K + base * (0.006 * ALU_J_PER_M3K + 0.0025 * STEEL_J_PER_M3K),
    lidJK: base * 0.001 * STEEL_J_PER_M3K,
  };
}
export const POT_20QT = stockpot("20qt", "20-quart", 0.31, 0.295);
export const POT_40QT = stockpot("40qt", "40-quart", 0.38, 0.375);

/** Where each loss went, watts, at one moment. */
export interface Losses {
  steam: number;
  radiation: number;
  air: number;
}
/** The heat a pot of water at `c` °C loses (W): `d` inside diameter, `wet` the water's depth, `lid` on or off,
 *  `film` an oil film on the water (Berick's control: no evaporation, ε 0.85), in a room at `room` °C. */
export function lossesW(c: number, d: number, wet: number, lid: boolean, room = ROOM_C, film = false): Losses {
  const top = (Math.PI * d * d) / 4;
  const side = Math.PI * d * wet;
  const steam = film ? 0 : (evaporationGPerHM2(c) / 3.6e6) * (BERICK_D_M / d) ** 0.25 * top * latentJPerKg(c) * (lid ? LID_LETS_OUT : 1);
  const k4 = (c + 273.15) ** 4 - (room + 273.15) ** 4;
  const eTop = lid ? EPS_SHINY_STEEL : film ? 0.85 : EPS_WATER;
  const radiation = SIGMA * (eTop * top + EPS_SHINY_STEEL * side) * k4;
  const air = (hSide(c, room, wet) * side + hTop(c, room, d) * top) * (c - room);
  return { steam, radiation, air };
}

/** One pot heating: what it starts as and what heats it. */
export interface HeatUp {
  /** Inside diameter (m), litres of water, the shell's heat capacity (J/K). */
  d: number;
  litres: number;
  shellJK: number;
  lid: boolean;
  /** Watts reaching the pot. */
  powerW: number;
  waterC?: number;
  roomC?: number;
  boilC?: number;
  film?: boolean;
  /** Stop here if it has not boiled (minutes). */
  maxMin?: number;
}
/** A pot's run: its temperature and the litres left, every `STEP_MIN`; when it reached the boil, if it did. */
export interface HeatRun {
  t: number[];
  c: number[];
  litres: number[];
  /** Minutes to the boil, or null if it was still short of it at `maxMin`. */
  boilMin: number | null;
  /** Litres gone to steam by the boil (or by the end). */
  steamL: number;
  /** Where the heat went, by the boil (or the end): J into the pot, to steam, radiated, to the air. */
  heat: { into: number; steam: number; radiation: number; air: number };
}
export const STEP_MIN = 0.25;
const DT_S = 1;

export function heatUp(p: HeatUp): HeatRun {
  const room = p.roomC ?? ROOM_C;
  const boil = p.boilC ?? BOIL_C;
  const maxS = (p.maxMin ?? 240) * 60;
  const top = (Math.PI * p.d * p.d) / 4;
  let kg = p.litres * 0.999;
  const kg0 = kg;
  let c = p.waterC ?? TAP_C;
  const run: HeatRun = { t: [0], c: [c], litres: [p.litres], boilMin: null, steamL: 0, heat: { into: 0, steam: 0, radiation: 0, air: 0 } };
  let next = STEP_MIN * 60;
  for (let s = 0; s < maxS; s += DT_S) {
    const wet = kg / 998 / top;
    const l = lossesW(c, p.d, wet, p.lid, room, p.film);
    const q = p.powerW - l.steam - l.radiation - l.air;
    const c1 = c + (q * DT_S) / (kg * C_WATER + p.shellJK);
    kg -= (l.steam * DT_S) / latentJPerKg(c);
    run.heat.into += p.powerW * DT_S;
    run.heat.steam += l.steam * DT_S;
    run.heat.radiation += l.radiation * DT_S;
    run.heat.air += l.air * DT_S;
    if (c1 >= boil) {
      // The second it crossed: in proportion to the rise left.
      const f = (boil - c) / (c1 - c);
      run.boilMin = (s + f * DT_S) / 60;
      c = boil;
      run.t.push(run.boilMin);
      run.c.push(boil);
      run.litres.push(kg / 0.999);
      break;
    }
    c = c1;
    if (s + DT_S >= next - 1e-9) {
      run.t.push(next / 60);
      run.c.push(c);
      run.litres.push(kg / 0.999);
      next += STEP_MIN * 60;
    }
  }
  run.steamL = (kg0 - kg) / 0.999;
  return run;
}

/** The four ways the page puts the same water on: one big pot or two smaller ones on two burners, lid on or off. */
export interface Setup {
  key: "one-open" | "one-lid" | "two-open" | "two-lid";
  name: string;
  pots: 1 | 2;
  pot: PotSpec;
  lid: boolean;
}
export const SETUPS: readonly Setup[] = [
  { key: "one-open", name: "One pot, no lid", pots: 1, pot: POT_40QT, lid: false },
  { key: "one-lid", name: "One pot, lid on", pots: 1, pot: POT_40QT, lid: true },
  { key: "two-open", name: "Two pots, no lids", pots: 2, pot: POT_20QT, lid: false },
  { key: "two-lid", name: "Two pots, lids on", pots: 2, pot: POT_20QT, lid: true },
];
/** The most water the page offers: two 20-quart pots filled to nine tenths. */
export const MAX_LITRES = Math.floor(2 * 0.9 * POT_20QT.litres);

export interface SetupRun extends HeatRun {
  setup: Setup;
  /** Litres in each pot. */
  each: number;
  /** Burner-minutes of gas to the boil (each burner lit until its pot boils), or to the end if it never did. */
  burnerMin: number;
}
/** Each set-up's run for `litres` of water in all, the burners at `flame` (0..1 of full), cut off at `maxMin`. */
export function boilAll(litres: number, flame: number, maxMin = 180, waterC = TAP_C): SetupRun[] {
  return SETUPS.map((setup) => {
    const each = litres / setup.pots;
    const run = heatUp({ d: setup.pot.d, litres: each, shellJK: setup.pot.shellJK + (setup.lid ? setup.pot.lidJK : 0), lid: setup.lid, powerW: BURNER_W * GAS_INTO_POT * flame, waterC, maxMin });
    return { ...run, setup, each, burnerMin: setup.pots * (run.boilMin ?? maxMin) };
  });
}
