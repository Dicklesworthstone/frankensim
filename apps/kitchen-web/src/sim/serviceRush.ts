/**
 * The lunch rush: one service on a fast-casual line, stepped every 15 seconds.
 *
 * Orders come in from 11:00 to 14:00 round a peak near 12:20. The forecast knows the day's size and the peak's timing
 * only roughly. Three cooks each run a station (the fryer, the griddle, the combi) and plate their own dishes. Each
 * dish is cooked one of two ways:
 *
 *   "order"  A batch is started for the orders waiting as soon as its station has a basket, a zone or a shelf free
 *            and its cook is free. The guest waits for the queue, the cooking and the plating.
 *   "hold"   The dish is cooked ahead and boxed into the heated holding cabinet, enough to cover the forecast's
 *            orders over the next stretch of the service (a share of its holding time), from 10:00 on, while the
 *            cabinet has room. The boxing is done with the batch, so a guest who finds it on the shelf waits only for
 *            it to be bagged. Portions past their holding time, and all held food left at close, are thrown out.
 *
 * A free cook first starts a batch that guests are waiting for, then plates the oldest order whose food is ready,
 * then tops up the dishes held. `simulateService` runs one day; `runRush` runs many on a fixed seed, so the same
 * settings always give the same answer.
 */
import { rng, spread, type Spread } from "./commissaryWeek";

export type StationKey = "fryer" | "griddle" | "combi";
export type Mode = "order" | "hold";

export interface RushDish {
  key: string;
  name: string;
  station: StationKey;
  /** Share of the orders. */
  share: number;
  /** Portions in one batch: a basket, a griddle load, a hotel pan. */
  batch: number;
  /** Minutes a batch cooks, unattended once it is in. */
  cookMin: number;
  /** Hands-on minutes per batch (setting up, pulling and panning it when done) and per portion in it. */
  loadBase: number;
  loadPer: number;
  /** Minutes a portion keeps in the holding cabinet. */
  holdMin: number;
  /** Portions in one half pan on a cabinet shelf. */
  perPan: number;
  /** Hotel-pan food the 3D view draws (smallwares food index). */
  food: 0 | 1 | 2 | 3;
}

export const RUSH_MENU: readonly RushDish[] = [
  { key: "tenders", name: "Chicken tenders", station: "fryer", share: 0.28, batch: 6, cookMin: 5, loadBase: 0.3, loadPer: 0.05, holdMin: 20, perPan: 12, food: 0 },
  { key: "burger", name: "Smash burger", station: "griddle", share: 0.24, batch: 8, cookMin: 4, loadBase: 0.2, loadPer: 0.25, holdMin: 8, perPan: 8, food: 0 },
  { key: "thighs", name: "Grilled chicken thighs", station: "griddle", share: 0.16, batch: 16, cookMin: 8, loadBase: 0.3, loadPer: 0.06, holdMin: 45, perPan: 12, food: 2 },
  { key: "mac", name: "Mac and cheese", station: "combi", share: 0.14, batch: 24, cookMin: 9, loadBase: 0.4, loadPer: 0.05, holdMin: 60, perPan: 12, food: 1 },
  { key: "beef", name: "Braised beef, reheated", station: "combi", share: 0.18, batch: 24, cookMin: 7, loadBase: 0.4, loadPer: 0.05, holdMin: 120, perPan: 12, food: 0 },
];

/**
 * The menu re-engineered for space: the thighs roast on sheet pans in the combi instead of grilling, taking the
 * griddle's second dish to the station with room in the rush.
 */
export const RUSH_SWAP: RushDish = { key: "thighs", name: "Roast chicken thighs", station: "combi", share: 0.16, batch: 24, cookMin: 12, loadBase: 0.4, loadPer: 0.04, holdMin: 45, perPan: 12, food: 2 };
export type RushMenu = "written" | "balanced";
export const rushMenu = (m: RushMenu): RushDish[] => RUSH_MENU.map((d) => (m === "balanced" && d.key === RUSH_SWAP.key ? RUSH_SWAP : d));

/** Hands-on minutes for a batch of n portions. */
export const loadOf = (d: RushDish, n: number) => d.loadBase + d.loadPer * n;

export const STATIONS: Readonly<Record<StationKey, { name: string; slots: number }>> = {
  fryer: { name: "Fryer", slots: 4 },
  griddle: { name: "Griddle", slots: 4 },
  combi: { name: "Combi", slots: 5 },
};
export const STATION_KEYS: readonly StationKey[] = ["fryer", "griddle", "combi"];

/** The service's clock, in minutes from 10:00. */
export const SERVICE = {
  /** Orders from 11:00 until 14:00; the kitchen runs from 10:00 to 14:30, and on past it if guests are still waiting. */
  open: 60,
  lastOrder: 240,
  close: 270,
  /** The peak of the forecast's order curve (12:20) and its spread. */
  peak: 140,
  peakSd: 24,
  /** The curve's floor: the share of its height that runs all through the service. */
  floor: 0.2,
  /** Minutes to plate an order at the station: the sides, the sauce, the box, the ticket. */
  plateMin: 0.6,
  /** A held dish is boxed as it comes out of the oven (minutes a portion, done with the batch), and at the order
   *  only bagged and handed over. */
  packMin: 0.5,
  handMin: 0.15,
  /** Half-pan places in the holding cabinet: 2.5-inch pans, two to a rail, on the 15 rails above its water pan. */
  cabinetPans: 30,
  /** A held dish is cooked to cover the forecast's orders for this share of its holding time ahead, once the
   *  shortfall reaches this share of a batch. */
  cover: 0.7,
  topUp: 0.25,
  /** The simulation step (minutes), and a hard stop for a kitchen still serving. */
  dt: 0.25,
  end: 330,
} as const;

/** Meals forecast for a quiet day and a busy one. */
export const RUSH_DAYS = { quiet: 220, busy: 420 } as const;
export type RushDay = keyof typeof RUSH_DAYS;

export interface RushSettings {
  /** How each dish is cooked, by key (missing: to order). */
  modes: Readonly<Record<string, Mode>>;
  day: RushDay;
  /** Forecast error: the day's size (lognormal spread) and the peak's timing (50 minutes per unit). */
  sigma: number;
  /** The menu as written, or with the thighs moved from the griddle to the combi. */
  menu?: RushMenu;
  /** Whether a cook with nothing of their own to do plates for the other stations. */
  pool?: boolean;
}

/** The dishes that keep at least half an hour: what "hold what keeps" holds. */
export const KEEPERS: readonly string[] = RUSH_MENU.filter((d) => d.holdMin >= 30).map((d) => d.key);

/** The forecast's order curve: the expected share of the day's orders in minute `m` (with the peak moved by `shift`). */
const CURVE = (() => {
  const raw = (m: number, shift: number) => {
    if (m < SERVICE.open || m >= SERVICE.lastOrder) return 0;
    const z = (m + 0.5 - SERVICE.peak - shift) / SERVICE.peakSd;
    return SERVICE.floor + Math.exp(-0.5 * z * z);
  };
  return (shift: number) => {
    const w = new Float64Array(SERVICE.end);
    let sum = 0;
    for (let m = 0; m < SERVICE.end; m++) sum += w[m] = raw(m, shift);
    for (let m = 0; m < SERVICE.end; m++) w[m]! /= sum;
    return w;
  };
})();
const FORECAST = CURVE(0);

/** Expected orders of a dish between minute a and b on the forecast (fractional minutes allowed). */
function forecastBetween(meals: number, share: number, a: number, b: number): number {
  let s = 0;
  const m0 = Math.max(0, Math.floor(a));
  const m1 = Math.min(SERVICE.end, Math.ceil(b));
  for (let m = m0; m < m1; m++) s += FORECAST[m]! * (Math.min(b, m + 1) - Math.max(a, m));
  return meals * share * s;
}

export interface Order {
  dish: number;
  arrive: number;
  /** When its food was there for it: on the shelf when it came, or out of the batch it waited for. */
  food: number;
  plated: number;
}

export interface DayRun {
  meals: number;
  orders: Order[];
  /** Per dish: portions cooked, and thrown out (past their time, or left at close). */
  cooked: number[];
  wasted: number[];
  /** Per cook (one per station), busy minutes in each minute of the clock, and the part of them spent plating. */
  busy: Float32Array[];
  plating: Float32Array[];
  /** Per held dish, portions in the cabinet at each minute; per station, batches cooking at each minute. */
  stock: Float32Array[];
  cooking: Float32Array[];
  /** Ticket times (minutes), sorted. */
  waits: Float64Array;
  waitP50: number;
  waitP90: number;
  /** Idle crew time between 10:00 and 14:30, as a share of it. */
  idle: number;
  /** The crew's busy share before the doors open (10:00 to 11:00), and in its busiest hour. */
  prepLoad: number;
  rushLoad: number;
  /** Minutes past 14:30 the last order was plated (0 if before). */
  over: number;
  wasteTotal: number;
}

interface Batch {
  dish: number;
  portions: number;
  from: number;
  done: number;
}

const quantile = (s: ArrayLike<number>, p: number) => {
  if (!s.length) return 0;
  const x = p * (s.length - 1);
  const i = Math.floor(x);
  return s[i]! + ((s[Math.min(s.length - 1, i + 1)] ?? s[i]!) - s[i]!) * (x - i);
};

/** One service. `meals` is the forecast; the day's own size and peak come from the seed and `sigma`. */
export function simulateService(settings: RushSettings, seed: number): DayRun {
  const menu = rushMenu(settings.menu ?? "written");
  const S = SERVICE;
  const forecast = RUSH_DAYS[settings.day];
  const r = rng(seed);
  const sig = settings.sigma;
  const meals = forecast * Math.exp(sig * r.z() - (sig * sig) / 2);
  const shift = 50 * sig * r.z();
  const curve = CURVE(shift);
  const held = menu.map((d) => settings.modes[d.key] === "hold");
  const stations = STATION_KEYS.map((k) => ({ key: k, slots: new Float64Array(STATIONS[k].slots) }));
  const stIdx = menu.map((d) => STATION_KEYS.indexOf(d.station));
  const cum: number[] = [];
  menu.reduce((a, d, i) => (cum[i] = a + d.share), 0);
  const shareSum = cum[cum.length - 1]!;

  const orders: Order[] = [];
  /** Per dish: orders waiting for food (indices into `orders`), and held portions (the time each was cooked). */
  const queue: number[][] = menu.map(() => []);
  const shelf: number[][] = menu.map(() => []);
  /** Orders with their food, waiting to be plated, per station (oldest first by arrival). */
  const ready: number[][] = STATION_KEYS.map(() => []);
  const cooking: Batch[] = [];
  const cooked = menu.map(() => 0);
  const wasted = menu.map(() => 0);
  const minutes = S.end;
  const busy = STATION_KEYS.map(() => new Float32Array(minutes));
  const plating = STATION_KEYS.map(() => new Float32Array(minutes));
  const stock = menu.map(() => new Float32Array(minutes));
  const cookingAt = STATION_KEYS.map(() => new Float32Array(minutes));
  const free = STATION_KEYS.map(() => 0);
  let plated = 0;

  const addBusy = (c: number, a: number, b: number, into = busy) => {
    for (let m = Math.floor(a); m < b && m < minutes; m++) into[c]![m]! += Math.min(b, m + 1) - Math.max(a, m);
  };
  /** Portions of a dish on the way: cooking for it, not yet promised to a guest waiting. */
  const inOven = (d: number) => cooking.reduce((s, b) => s + (b.dish === d ? b.portions : 0), 0);
  const pansHeld = () => menu.reduce((s, d, i) => s + (held[i] ? Math.ceil((shelf[i]!.length + inOven(i)) / d.perPan) : 0), 0);
  const giveFood = (oi: number, t: number) => {
    const o = orders[oi]!;
    o.food = t;
    const st = ready[stIdx[o.dish]!]!;
    // Kept in arrival order: the oldest guest is plated first.
    let k = st.length;
    while (k > 0 && orders[st[k - 1]!]!.arrive > o.arrive) k--;
    st.splice(k, 0, oi);
  };
  const freeSlot = (s: number, t: number) => {
    const sl = stations[s]!.slots;
    for (let i = 0; i < sl.length; i++) if (sl[i]! <= t + 1e-9) return i;
    return -1;
  };
  const start = (d: number, portions: number, s: number, slot: number, t: number, c: number) => {
    const dish = menu[d]!;
    // A held batch is boxed as it is panned: that work moves from the orders to the batch.
    const load = loadOf(dish, portions) + (held[d] ? S.packMin * portions : 0);
    const done = t + load + dish.cookMin;
    cooking.push({ dish: d, portions, from: t, done });
    stations[s]!.slots[slot] = done;
    cooked[d]! += portions;
    free[c] = t + load;
    addBusy(c, t, t + load);
  };
  /**
   * How far a held dish is below its cover at time t (portions): the forecast's orders from now until a stretch of
   * its holding time past when a batch started now would be done, less what is on the shelf or cooking for it and
   * not promised to a guest already waiting.
   */
  const holdShort = (d: number, t: number) => {
    const dish = menu[d]!;
    const due = t + loadOf(dish, dish.batch) + S.packMin * dish.batch + dish.cookMin;
    const want = forecastBetween(forecast, dish.share / shareSum, t, Math.min(S.lastOrder, due + dish.holdMin * S.cover));
    return want - (shelf[d]!.length + inOven(d) - queue[d]!.length);
  };
  /** A held dish's batch: what covers the shortfall, up to a full one. */
  const holdBatch = (d: number, short: number) => Math.max(1, Math.min(menu[d]!.batch, Math.ceil(short)));

  /** What a free cook does next at time t; false when there is nothing to do. */
  const act = (c: number, t: number): boolean => {
    // 1. A batch guests are waiting for: the dish with the longest-waiting guest not yet covered by a batch.
    let best = -1;
    let bestAge = -1;
    let bestN = 0;
    for (let d = 0; d < menu.length; d++) {
      if (stIdx[d] !== c) continue;
      const q = queue[d]!;
      const uncovered = q.length - inOven(d);
      if (uncovered <= 0) continue;
      const oldest = orders[q[q.length - uncovered]!]!.arrive;
      if (t - oldest > bestAge) {
        bestAge = t - oldest;
        best = d;
        bestN = uncovered;
      }
    }
    if (best >= 0) {
      const slot = freeSlot(c, t);
      if (slot >= 0) {
        const dish = menu[best]!;
        // One cooked to order: what is waiting. One held and run out: that, and its cover again, onto the shelf.
        start(best, held[best] ? Math.max(Math.min(dish.batch, bestN), holdBatch(best, holdShort(best, t))) : Math.min(dish.batch, bestN), c, slot, t, c);
        return true;
      }
    }
    // 2. Plate the oldest order whose food is ready: this station's, or with pooling, any station's.
    let from = ready[c]!.length ? c : -1;
    if (from < 0 && settings.pool) {
      let oldest = Infinity;
      ready.forEach((rq, k) => {
        if (rq.length && orders[rq[0]!]!.arrive < oldest) {
          oldest = orders[rq[0]!]!.arrive;
          from = k;
        }
      });
    }
    if (from >= 0) {
      const o = orders[ready[from]!.shift()!]!;
      o.plated = t + (held[o.dish] ? S.handMin : S.plateMin);
      plated++;
      free[c] = o.plated;
      addBusy(c, t, o.plated);
      addBusy(c, t, o.plated, plating);
      return true;
    }
    // 3. Top up a held dish below its cover, the one furthest below first.
    let need = 0;
    let pick = -1;
    let pickShort = 0;
    for (let d = 0; d < menu.length; d++) {
      if (!held[d] || stIdx[d] !== c) continue;
      const short = holdShort(d, t);
      // Topped up once a quarter of a batch is short: smaller batches are not worth the load.
      if (short >= Math.max(1, menu[d]!.batch * S.topUp) && short / menu[d]!.batch > need) {
        need = short / menu[d]!.batch;
        pick = d;
        pickShort = short;
      }
    }
    if (pick >= 0) {
      const slot = freeSlot(c, t);
      const dish = menu[pick]!;
      const n = holdBatch(pick, pickShort);
      if (slot >= 0 && pansHeld() + Math.ceil(n / dish.perPan) <= S.cabinetPans) {
        start(pick, n, c, slot, t, c);
        return true;
      }
    }
    return false;
  };

  let t = 0;
  let lastPlated = 0;
  for (; t < S.end; t += S.dt) {
    // Batches done: their portions go to the guests waiting, the rest to the shelf.
    for (let i = cooking.length - 1; i >= 0; i--) {
      const b = cooking[i]!;
      if (b.done > t + 1e-9) continue;
      cooking.splice(i, 1);
      let n = b.portions;
      const q = queue[b.dish]!;
      while (n > 0 && q.length) {
        giveFood(q.shift()!, b.done);
        n--;
      }
      for (; n > 0; n--) shelf[b.dish]!.push(b.done);
    }
    // Held food past its time is thrown out.
    menu.forEach((dish, d) => {
      const sh = shelf[d]!;
      while (sh.length && t - sh[0]! > dish.holdMin) {
        sh.shift();
        wasted[d]!++;
      }
    });
    // Guests arrive.
    const m = Math.floor(t);
    if (t >= S.open && t < S.lastOrder) {
      const lam = meals * curve[m]! * S.dt;
      // Poisson by inversion (lambda is well under 10 a step).
      let n = 0;
      let p = Math.exp(-lam);
      let F = p;
      const u = r.u();
      while (u > F && n < 50) {
        n++;
        p *= lam / n;
        F += p;
      }
      for (let k = 0; k < n; k++) {
        const v = r.u() * shareSum;
        let d = 0;
        while (d < cum.length - 1 && v > cum[d]!) d++;
        // On the step (15 s): a cook acting within it never starts on an order before it came in.
        const oi = orders.push({ dish: d, arrive: t, food: NaN, plated: NaN }) - 1;
        if (shelf[d]!.length) {
          shelf[d]!.shift();
          giveFood(oi, t);
        } else queue[d]!.push(oi);
      }
    }
    // Each cook works through the step.
    for (let c = 0; c < STATION_KEYS.length; c++) {
      let guard = 0;
      while (free[c]! < t + S.dt && guard++ < 40) {
        const at = Math.max(free[c]!, t);
        if (!act(c, at)) break;
      }
    }
    if (m < minutes) {
      menu.forEach((_, d) => (stock[d]![m] = shelf[d]!.length));
      // Batch-minutes a minute: each step adds its share, so a minute sums to the batches cooking through it.
      const mid = t + S.dt / 2;
      for (const b of cooking) if (b.from <= mid && b.done > mid) cookingAt[stIdx[b.dish]!]![m]! += S.dt;
    }
    if (t >= S.close && plated === orders.length) break;
  }
  // Close: what is left on the shelf is thrown out; any guest never plated waited until the stop.
  menu.forEach((_, d) => {
    wasted[d]! += shelf[d]!.length;
    shelf[d]!.length = 0;
  });
  for (const o of orders) {
    if (Number.isNaN(o.plated)) o.plated = S.end;
    if (Number.isNaN(o.food)) o.food = o.plated;
    lastPlated = Math.max(lastPlated, o.plated);
  }
  const waits = Float64Array.from(orders.map((o) => o.plated - o.arrive)).sort();
  let worked = 0;
  for (const b of busy) for (let m = 0; m < S.close; m++) worked += b[m]!;
  const crew = new Float64Array(minutes + 1);
  for (let m = 0; m < minutes; m++) crew[m + 1] = crew[m]! + busy.reduce((a, b) => a + b[m]!, 0);
  let rushLoad = 0;
  for (let m = 0; m + 60 <= minutes; m++) rushLoad = Math.max(rushLoad, crew[m + 60]! - crew[m]!);
  rushLoad /= 60 * STATION_KEYS.length;
  const prepLoad = crew[S.open]! / (S.open * STATION_KEYS.length);
  return {
    meals: orders.length,
    orders,
    cooked,
    wasted,
    busy,
    plating,
    stock,
    cooking: cookingAt,
    waits,
    waitP50: quantile(waits, 0.5),
    waitP90: quantile(waits, 0.9),
    idle: 1 - worked / (STATION_KEYS.length * S.close),
    prepLoad,
    rushLoad,
    over: Math.max(0, lastPlated - S.close),
    wasteTotal: wasted.reduce((a, b) => a + b, 0),
  };
}

export interface RushRun {
  settings: RushSettings;
  days: number;
  waitP50: Spread;
  waitP90: Spread;
  idle: Spread;
  prepLoad: Spread;
  rushLoad: Spread;
  waste: Spread;
  /** The day the 3D view and the chart show: the median day by orders. Its orders do not depend on the plan, the
   *  menu or the crew, so every setting is shown on the same day. */
  shown: DayRun;
}

export function runRush(settings: RushSettings, days = 200, seed = 20261124): RushRun {
  const runs: DayRun[] = [];
  for (let i = 0; i < days; i++) runs.push(simulateService(settings, seed + i * 7919));
  const p90s = runs.map((d) => d.waitP90);
  const byOrders = runs.map((d, i) => [d.meals, i] as const).sort((a, b) => a[0] - b[0] || a[1] - b[1]);
  const shown = runs[byOrders[Math.floor(days / 2)]![1]]!;
  return {
    settings,
    days,
    waitP50: spread(runs.map((d) => d.waitP50)),
    waitP90: spread(p90s),
    idle: spread(runs.map((d) => d.idle)),
    prepLoad: spread(runs.map((d) => d.prepLoad)),
    rushLoad: spread(runs.map((d) => d.rushLoad)),
    waste: spread(runs.map((d) => d.wasteTotal)),
    shown,
  };
}

export const modesFor = (hold: readonly string[]): Record<string, Mode> => Object.fromEntries(RUSH_MENU.map((d) => [d.key, hold.includes(d.key) ? "hold" : "order"]));
