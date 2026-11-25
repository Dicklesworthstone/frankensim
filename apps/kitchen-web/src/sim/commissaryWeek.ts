/**
 * A commissary's week: one crew, one menu, and a forecast that is never quite right.
 *
 * Six production days, Monday to Saturday. Each day's orders are known that morning; the days after it are only
 * forecast. Every meal takes some of each dish on the menu, and every dish costs the crew hands-on minutes per meal.
 * A dish that keeps (cook-chilled in the blast chiller, stored in the walk-in, reheated and held hot on the day) can
 * be cooked on an earlier, quieter day; a dish that does not keep is cooked the day it is served.
 *
 *   "day"    Every dish is cooked the day it is served. Slow days leave the crew idle; busy days run into overtime,
 *            and past the overtime cap, meals go out late.
 *   "ahead"  Each morning, once the day's own work is counted, the spare crew time cooks for the coming days that
 *            the forecast says will overflow: as much of each keeping dish as levels them (plus a margin that grows
 *            with the forecast error), the dishes that move the most minutes per pan of rack space first, the
 *            batch for the nearest busy day first, within the chiller's pans a day and the walk-in's racks.
 *            Chilled stock is used oldest first; what passes its keeping time is thrown out.
 *
 * Demand: each day's meals are the forecast times a lognormal error with spread `sigma`, correlated across the week
 * (a slow week tends to be slow every day). `simulateWeek` runs one week; `runCommissary` runs the forecast week
 * itself and a Monte Carlo of many weeks on a fixed seed, so the same settings always give the same answer.
 */

export const DAYS = ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat"] as const;
export type Station = "kettle" | "skillet" | "combi" | "prep";

export interface Dish {
  key: string;
  name: string;
  station: Station;
  /** Hands-on minutes per meal that includes it (the batch's work spread over its portions). */
  laborMin: number;
  /** Days a cook-chilled batch keeps in the walk-in: 0 means it is served the day it is cooked. */
  keeps: number;
  /** Share of the day's meals that include it. */
  share: number;
  /** Portions in one 4 in full hotel pan. */
  perPan: number;
  /** Hotel-pan food the 3D view draws (smallwares food index). */
  food: number;
}

/** The menu as written. Three dishes have a variant that keeps (see `SWAPS`). */
export const MENU: readonly Dish[] = [
  { key: "chicken", name: "Grilled chicken", station: "skillet", laborMin: 0.6, keeps: 0, share: 0.45, perPan: 24, food: 0 },
  { key: "shortrib", name: "Short rib braise", station: "combi", laborMin: 0.5, keeps: 4, share: 0.25, perPan: 20, food: 0 },
  { key: "pasta", name: "Pasta to order", station: "kettle", laborMin: 0.35, keeps: 0, share: 0.35, perPan: 24, food: 1 },
  { key: "rice", name: "Rice pilaf", station: "combi", laborMin: 0.15, keeps: 1, share: 0.45, perPan: 40, food: 1 },
  { key: "soup", name: "Tomato soup", station: "kettle", laborMin: 0.2, keeps: 3, share: 0.4, perPan: 36, food: 3 },
  { key: "salad", name: "Green salad", station: "prep", laborMin: 0.45, keeps: 0, share: 0.45, perPan: 30, food: 2 },
];

/** Dishes that do the same job on the menu and keep: braised rather than grilled, baked rather than boiled to order. */
export const SWAPS: Readonly<Record<string, Dish>> = {
  chicken: { key: "chicken", name: "Braised chicken thighs", station: "skillet", laborMin: 0.5, keeps: 3, share: 0.45, perPan: 24, food: 0 },
  pasta: { key: "pasta", name: "Baked rigatoni", station: "combi", laborMin: 0.4, keeps: 2, share: 0.35, perPan: 24, food: 1 },
  salad: { key: "salad", name: "Farro and herb salad", station: "prep", laborMin: 0.4, keeps: 2, share: 0.45, perPan: 30, food: 2 },
};

/** The crew, the room and the week's forecast. */
export interface Kitchen {
  cooks: number;
  shiftMin: number;
  /** Overtime per cook that still gets the day's meals out on time (minutes); work past it leaves meals late. */
  otCapMin: number;
  /** Setup and cleaning per cook per day (minutes): paid, worked, and the same every day. */
  fixedMin: number;
  /** Portioning and packing, minutes per meal, always on the day. */
  packMin: number;
  /** Reheating a chilled portion and holding it hot, minutes per portion, on the day it is served. */
  rethermMin: number;
  /** Blast chiller: pans per cycle and minutes per cycle (a cycle runs only inside the shift). */
  chillPans: number;
  chillCycleMin: number;
  /** Walk-in: speed racks and pans per rack. */
  racks: number;
  rackPans: number;
  /** Meals forecast for each production day. */
  forecast: readonly number[];
}

export const KITCHEN: Kitchen = {
  cooks: 3,
  shiftMin: 480,
  otCapMin: 120,
  fixedMin: 45,
  packMin: 0.55,
  rethermMin: 0.05,
  chillPans: 15,
  chillCycleMin: 90,
  racks: 3,
  rackPans: 20,
  forecast: [600, 360, 1250, 760, 1120, 540],
};

export type Strategy = "day" | "ahead";

export interface Settings {
  strategy: Strategy;
  /** Keys of the dishes swapped for the variant that keeps. */
  swaps: readonly string[];
  /** Forecast error: the spread (lognormal sigma) of actual meals around the forecast, 0..0.5. */
  sigma: number;
  /** Speed racks in the walk-in, if not the kitchen's own. */
  racks?: number;
}

/** Cooking done ahead on one day for a later day. */
export interface Batch {
  dish: string;
  portions: number;
  pans: number;
  minutes: number;
  /** The day it was cooked for. */
  forDay: number;
}

export interface DayResult {
  meals: number;
  /** Minutes of work today's meals took on the day (cooking what was not in stock, reheating what was, packing,
   *  setup and cleaning). */
  dayOfMin: number;
  /** Minutes spent cooking ahead for later days. */
  aheadMin: number;
  /** The batches cooked ahead today. */
  batches: Batch[];
  /** Overtime worked, crew total (minutes). */
  otMin: number;
  /** Paid crew time nobody had work for (minutes). */
  idleMin: number;
  /** Overtime past what the delivery window allows (minutes), and the meals that leaves late. */
  lateMin: number;
  lateMeals: number;
  /** Portions served from chilled stock, and cooked on the day. */
  fromStock: number;
  cookedToday: number;
  /** Portions thrown out at the end of the day (past their keeping time). */
  wasted: number;
  /** Pans in the walk-in at the start and the end of the day, and pans through the chiller today. */
  pansStart: number;
  pansEnd: number;
  pansChilled: number;
  /** Crew minutes on the clock today (shift plus overtime) and of them worked. */
  paidMin: number;
  workedMin: number;
  /** Minutes of each station's cooking today (on the day and ahead). */
  stationMin: Record<Station, number>;
}

export interface WeekResult {
  days: DayResult[];
  otMin: number;
  idleMin: number;
  lateMeals: number;
  wasted: number;
  /** Portions still good in the walk-in at the end of Saturday (they serve next week, not counted as waste). */
  carried: number;
  /** Worked minutes over paid minutes. */
  utilization: number;
  meals: number;
}

/** The menu as it stands with these swaps. */
export function menuFor(swaps: readonly string[]): Dish[] {
  return MENU.map((d) => (swaps.includes(d.key) ? SWAPS[d.key]! : d));
}

interface Lot {
  dish: number;
  portions: number;
  /** Last day it may be served. */
  until: number;
  forDay: number;
}

const pansOf = (portions: number, perPan: number) => (portions <= 1e-9 ? 0 : Math.ceil(portions / perPan - 1e-9));

/**
 * One week. `actual[d]` is day d's meals as they turn out; the planner sees `k.forecast` for the days after today.
 * `margin` is the share of extra portions planned above the forecast for days cooked ahead (from the error).
 */
export function simulateWeek(k: Kitchen, menu: readonly Dish[], strategy: Strategy, actual: readonly number[], margin: number): WeekResult {
  const n = k.forecast.length;
  const capMin = k.cooks * k.shiftMin;
  const fixed = k.cooks * k.fixedMin;
  const rackMax = k.racks * k.rackPans;
  const chillMaxPans = Math.floor(k.shiftMin / k.chillCycleMin) * k.chillPans;
  const lots: Lot[] = [];
  const days: DayResult[] = [];
  let wasted = 0;

  const stockPans = () => {
    let p = 0;
    for (const l of lots) p += pansOf(l.portions, menu[l.dish]!.perPan);
    return p;
  };
  /** Minutes a day's meals would take on the day, given the portions of each dish already cooked for it. */
  const dayOfWork = (meals: number, covered: readonly number[]) => {
    let m = fixed + meals * k.packMin;
    menu.forEach((dish, i) => {
      const need = meals * dish.share;
      const have = Math.min(need, covered[i] ?? 0);
      m += (need - have) * dish.laborMin + have * k.rethermMin;
    });
    return m;
  };

  for (let d = 0; d < n; d++) {
    const meals = actual[d]!;
    const pansStart = stockPans();
    const stationMin: Record<Station, number> = { kettle: 0, skillet: 0, combi: 0, prep: 0 };

    // Serve today: chilled stock first, oldest (soonest to expire) first; the rest is cooked now.
    lots.sort((a, b) => a.until - b.until || a.forDay - b.forDay);
    let fromStock = 0;
    let cookedToday = 0;
    let work = fixed + meals * k.packMin;
    stationMin.prep += meals * k.packMin;
    menu.forEach((dish, i) => {
      let need = meals * dish.share;
      for (const l of lots) {
        if (need <= 1e-9) break;
        if (l.dish !== i || l.portions <= 1e-9 || l.until < d) continue;
        const take = Math.min(need, l.portions);
        l.portions -= take;
        need -= take;
        fromStock += take;
        work += take * k.rethermMin;
        stationMin.combi += take * k.rethermMin;
      }
      cookedToday += need;
      work += need * dish.laborMin;
      stationMin[dish.station] += need * dish.laborMin;
    });
    for (let i = lots.length - 1; i >= 0; i--) if (lots[i]!.portions <= 1e-9) lots.splice(i, 1);

    let otMin = 0;
    let lateMin = 0;
    let idleMin = 0;
    let aheadMin = 0;
    let pansChilled = 0;
    const batches: Batch[] = [];
    if (work > capMin) {
      // The crew stays until the work is done. Past the overtime the delivery window allows, meals leave late.
      otMin = work - capMin;
      lateMin = Math.max(0, otMin - k.cooks * k.otCapMin);
    } else {
      let slack = capMin - work;
      if (strategy === "ahead") {
        // The coming days' overflow at the forecast (plus the margin), less what is already cooked for each.
        const coveredFor = (D: number) => {
          const c = menu.map(() => 0);
          for (const l of lots) if (l.until >= D && l.forDay === D) c[l.dish]! += l.portions;
          return c;
        };
        // Levelled to 92% of the shift: a little room left for the unforeseen.
        const level = 0.92 * capMin;
        const order = menu
          .map((dish, i) => ({ dish, i, perPan: dish.laborMin * dish.perPan }))
          .filter((x) => x.dish.keeps > 0 && x.dish.laborMin > k.rethermMin)
          .sort((a, b) => b.perPan - a.perPan || a.i - b.i);
        // Each day's spare time: today's for real, the coming days' at the forecast (negative where they overflow).
        const plan = (j: number) => k.forecast[j]! * (1 + margin);
        const spare = k.forecast.map((_, j) => (j < d ? 0 : j === d ? slack : level - dayOfWork(plan(j), coveredFor(j))));
        // Every overflowing day in turn, nearest first, pulls its keeping dishes back into the latest spare day the
        // dish keeps for: the food is as fresh as it can be, and the day before a rush is used before the one before
        // that. Only what lands on today is cooked now; the rest is booked, and planned again tomorrow morning.
        for (let D = d + 1; D < n; D++) {
          let over = -spare[D]!;
          if (over <= 1e-6) continue;
          const covered = coveredFor(D);
          for (const { dish, i } of order) {
            const saving = dish.laborMin - k.rethermMin;
            let want = Math.max(0, plan(D) * dish.share - covered[i]!);
            for (let j = D - 1; j >= Math.max(d, D - dish.keeps) && want > 1e-9 && over > 1e-6; j--) {
              if (spare[j]! <= 1e-6) continue;
              let portions = Math.min(want, over / saving, spare[j]! / dish.laborMin);
              if (j === d) {
                // Cooked today in whole pans: as many as the time, the walk-in's racks and the chiller allow.
                const pans = Math.min(
                  Math.round(portions / dish.perPan),
                  rackMax - stockPans(),
                  chillMaxPans - pansChilled,
                  Math.floor(spare[j]! / (dish.laborMin * dish.perPan) + 1e-9),
                );
                if (pans <= 0) continue;
                portions = pans * dish.perPan;
                const minutes = portions * dish.laborMin;
                lots.push({ dish: i, portions, until: d + dish.keeps, forDay: D });
                batches.push({ dish: dish.key, portions, pans, minutes, forDay: D });
                pansChilled += pans;
                stationMin[dish.station] += minutes;
                aheadMin += minutes;
              }
              spare[j]! -= portions * dish.laborMin;
              want -= portions;
              over -= portions * saving;
            }
          }
        }
        slack = spare[d]!;
      }
      idleMin = slack;
    }

    // End of the day: what has passed its keeping time goes.
    let dayWaste = 0;
    for (let i = lots.length - 1; i >= 0; i--) {
      if (lots[i]!.until <= d) {
        dayWaste += lots[i]!.portions;
        lots.splice(i, 1);
      }
    }
    wasted += dayWaste;

    const perMeal = meals > 0 ? (work - fixed) / meals : 0;
    const workedMin = work + aheadMin;
    days.push({
      meals,
      dayOfMin: work,
      aheadMin,
      batches,
      otMin,
      idleMin,
      lateMin,
      lateMeals: perMeal > 0 ? lateMin / perMeal : 0,
      fromStock,
      cookedToday,
      wasted: dayWaste,
      pansStart,
      pansEnd: stockPans(),
      pansChilled,
      paidMin: capMin + otMin,
      workedMin,
      stationMin,
    });
  }

  const sum = (f: (r: DayResult) => number) => days.reduce((s, r) => s + f(r), 0);
  const paid = sum((r) => r.paidMin);
  return {
    days,
    otMin: sum((r) => r.otMin),
    idleMin: sum((r) => r.idleMin),
    lateMeals: sum((r) => r.lateMeals),
    wasted,
    carried: lots.reduce((s, l) => s + l.portions, 0),
    utilization: paid > 0 ? sum((r) => r.workedMin) / paid : 0,
    meals: sum((r) => r.meals),
  };
}

/** Deterministic uniform and normal draws (mulberry32, Box-Muller). */
export function rng(seed: number): { u(): number; z(): number } {
  let s = seed >>> 0;
  const u = () => {
    s = (s + 0x6d2b79f5) >>> 0;
    let t = s;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
  let spare: number | null = null;
  const z = () => {
    if (spare !== null) {
      const v = spare;
      spare = null;
      return v;
    }
    const a = Math.max(1e-12, u());
    const b = u();
    const r = Math.sqrt(-2 * Math.log(a));
    spare = r * Math.sin(2 * Math.PI * b);
    return r * Math.cos(2 * Math.PI * b);
  };
  return { u, z };
}

export interface Spread {
  p10: number;
  p50: number;
  p90: number;
  mean: number;
}

export function spread(xs: number[]): Spread {
  const s = [...xs].sort((a, b) => a - b);
  const q = (p: number) => {
    const x = p * (s.length - 1);
    const i = Math.floor(x);
    return s[i]! + (s[Math.min(s.length - 1, i + 1)]! - s[i]!) * (x - i);
  };
  return { p10: q(0.1), p50: q(0.5), p90: q(0.9), mean: s.reduce((a, b) => a + b, 0) / Math.max(1, s.length) };
}

export interface CommissaryRun {
  settings: Settings;
  menu: Dish[];
  /** The forecast week, as planned (actual meals equal to the forecast). */
  plan: WeekResult;
  /** Across the simulated weeks. */
  weeks: number;
  ot: Spread;
  idle: Spread;
  late: Spread;
  wasted: Spread;
  utilization: Spread;
  /** Each day's work demanded on the day plus cooked ahead (crew minutes). */
  dayWork: Spread[];
}

/** Correlation of the forecast error between the days of one week. */
export const WEEK_RHO = 0.35;
/** Extra planned above the forecast for a day cooked ahead, per unit of forecast error. */
export const MARGIN_PER_SIGMA = 0.8;

export function runCommissary(settings: Settings, kitchen: Kitchen = KITCHEN, weeks = 400, seed = 20261124): CommissaryRun {
  const k = settings.racks === undefined ? kitchen : { ...kitchen, racks: settings.racks };
  const menu = menuFor(settings.swaps);
  const margin = MARGIN_PER_SIGMA * settings.sigma;
  const plan = simulateWeek(k, menu, settings.strategy, k.forecast, margin);
  const r = rng(seed);
  const sig = settings.sigma;
  const ot: number[] = [];
  const idle: number[] = [];
  const late: number[] = [];
  const wasted: number[] = [];
  const util: number[] = [];
  const dayWork: number[][] = k.forecast.map(() => []);
  for (let w = 0; w < weeks; w++) {
    const common = r.z();
    const actual = k.forecast.map((f) => {
      const z = Math.sqrt(WEEK_RHO) * common + Math.sqrt(1 - WEEK_RHO) * r.z();
      return Math.round(f * Math.exp(sig * z - (sig * sig) / 2));
    });
    const wk = simulateWeek(k, menu, settings.strategy, actual, margin);
    ot.push(wk.otMin / 60);
    idle.push(wk.idleMin / 60);
    late.push(wk.lateMeals);
    wasted.push(wk.wasted);
    util.push(wk.utilization);
    wk.days.forEach((dr, i) => dayWork[i]!.push(dr.dayOfMin + dr.aheadMin));
  }
  return {
    settings,
    menu,
    plan,
    weeks,
    ot: spread(ot),
    idle: spread(idle),
    late: spread(late),
    wasted: spread(wasted),
    utilization: spread(util),
    dayWork: dayWork.map(spread),
  };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* Every plan the week figure can show: the crew, when to cook, the walk-in's racks, and which dishes are swapped. */
/* ------------------------------------------------------------------------------------------------------------ */

export const PLAN_CREWS: readonly number[] = [2, 3, 4];
export const PLAN_RACKS = 5;
const SWAP_KEYS = Object.keys(SWAPS);
export const WEEK_PLANS = PLAN_CREWS.length * 2 * PLAN_RACKS * (1 << SWAP_KEYS.length);

export interface WeekPlan {
  cooks: number;
  strategy: Strategy;
  racks: number;
  swaps: string[];
}

/** Plan `i` of WEEK_PLANS: swaps vary fastest, then racks, then when to cook, then the crew. */
export function weekPlan(i: number): WeekPlan {
  const nm = 1 << SWAP_KEYS.length;
  const m = i % nm;
  const racks = 1 + (Math.floor(i / nm) % PLAN_RACKS);
  const strategy: Strategy = Math.floor(i / (nm * PLAN_RACKS)) % 2 === 0 ? "day" : "ahead";
  const cooks = PLAN_CREWS[Math.floor(i / (nm * PLAN_RACKS * 2))]!;
  return { cooks, strategy, racks, swaps: SWAP_KEYS.filter((_, b) => m & (1 << b)) };
}

/** The number of a plan (see weekPlan). */
export function weekPlanIndex(p: WeekPlan): number {
  const nm = 1 << SWAP_KEYS.length;
  const m = SWAP_KEYS.reduce((a, k, b) => (p.swaps.includes(k) ? a | (1 << b) : a), 0);
  return ((PLAN_CREWS.indexOf(p.cooks) * 2 + (p.strategy === "ahead" ? 1 : 0)) * PLAN_RACKS + (p.racks - 1)) * nm + m;
}

/** The kitchen with this many cooks. */
export const withCooks = (cooks: number, kitchen: Kitchen = KITCHEN): Kitchen => (cooks === kitchen.cooks ? kitchen : { ...kitchen, cooks });

/** Crew-hours paid in a week: every cook's shift on every day, and the overtime of a median week. */
export const paidHours = (run: CommissaryRun, cooks: number, kitchen: Kitchen = KITCHEN): number => (cooks * kitchen.shiftMin * kitchen.forecast.length) / 60 + run.ot.p50;
