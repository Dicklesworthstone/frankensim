/**
 * A week of meat and fish through a commissary's walk-in: ordered off the distributor's trucks, kept as long as the
 * cold storage chart allows, pulled for each day's catered lunches, thrown out once past its time.
 *
 * How long each keeps is the federal cold storage chart's (foodsafety.gov, as the Virginia Department of Health
 * publishes it, last reviewed 2021), refrigerator at 40 °F or below and freezer at 0 °F, counting the day it comes in
 * as day 1 the way the FDA Food Code counts date marking (§3-501.17): fresh chicken pieces and ground beef 1 to 2 days,
 * fatty fish such as salmon 1 to 3, beef roasts 3 to 5; frozen, 2 to 12 months, longer than any week here. Frozen
 * meat thaws in the walk-in: USDA's guidance is a full day even for a pound of ground meat or boneless chicken, and for
 * a large frozen item a day for every 4 to 5 lb. So loose salmon portions and chicken thighs, frozen piece by piece, thaw in
 * a day, and the 10 lb packs of ground beef and short ribs in two. Thawed, ground meat, poultry and seafood keep 1 to
 * 2 more days and red meat cuts 3 to 5. Each takes the longest of its range.
 *
 * What the kitchen knows when: the client confirms each lunch's final count by noon the day before. A truck's order
 * goes in by 2 pm the day before it comes (6 am), so it knows the final count for the truck's own day and only the
 * plan for later days; frozen meat is pulled to thaw at 8 am, its thawing time ahead of its lunch, before the count
 * is in, so that too goes by the plan. Ordering to the plan adds the cushion. Both come in whole cases, and frozen is
 * pulled to thaw by the case. Whatever is short at the morning's prep is bought locally on the day ("short"); whatever
 * is still there at the night's check on its last day is thrown out.
 *
 * Space: the meat and fish have three shelves of four case places in the walk-in, and as many in the freezer, two
 * cases of a delivery stacked in a place (a case is never stacked on another delivery's, so the older is never buried).
 * On the shelves raw meat stands by the temperature it is cooked to, the hottest lowest, so nothing drips onto a food
 * cooked less than itself (the Food Code's annex: chicken below beef), and fish above beef cooked to the same. A case
 * with no shelf place that keeps that order stands on the floor, which the Food Code does not allow for food
 * (§3-305.11: at least 6 inches above it).
 */

export const DAY = ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"] as const;
export const DAY3 = ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"] as const;

export type ItemKey = "salmon" | "ribs" | "beef" | "thighs";
export interface Item {
  key: ItemKey;
  name: string;
  dish: string;
  /** Raw weight in a lunch (lb). */
  lbPer: number;
  /** Fresh comes in cases of this many lb. */
  caseLb: number;
  /** Days fresh keeps in the walk-in, the day it comes in counted as day 1. */
  freshDays: number;
  /** Days it keeps once thawed, the day it is ready counted as day 1. */
  thawedDays: number;
  /** Days it takes to thaw in the walk-in: pulled at 8 am this many days before its lunch. */
  thawDays: number;
  /** Its quality frozen at 0 °F, months (the chart's longest). */
  frozenMonths: number;
  /** The Food Code's cooking temperature for it, °F (§3-401.11). */
  cookF: number;
}
/** In the order they stand on the shelves, top to bottom: the higher a food's cooking temperature, the lower it is
 *  kept, so nothing drips onto a food cooked less (the Food Code's annex: chicken below beef). */
export const ITEMS: readonly Item[] = [
  { key: "salmon", name: "Salmon fillets", dish: "Roast salmon", lbPer: 0.375, caseLb: 20, freshDays: 3, thawedDays: 2, thawDays: 1, frozenMonths: 3, cookF: 145 },
  { key: "ribs", name: "Beef short ribs", dish: "Short rib braise", lbPer: 0.75, caseLb: 40, freshDays: 5, thawedDays: 5, thawDays: 2, frozenMonths: 12, cookF: 145 },
  { key: "beef", name: "Ground beef", dish: "Chili", lbPer: 0.25, caseLb: 40, freshDays: 2, thawedDays: 2, thawDays: 2, frozenMonths: 4, cookF: 155 },
  { key: "thighs", name: "Chicken thighs", dish: "Roast chicken thighs", lbPer: 0.375, caseLb: 40, freshDays: 2, thawedDays: 2, thawDays: 1, frozenMonths: 9, cookF: 165 },
];
export const ITEM = Object.fromEntries(ITEMS.map((i) => [i.key, i])) as Record<ItemKey, Item>;

/** The week's catered lunches, as planned: a slow Tuesday before the big Wednesday. Nothing at the weekend. */
export const MENU: readonly { day: number; item: ItemKey; plan: number }[] = [
  { day: 0, item: "thighs", plan: 320 },
  { day: 1, item: "salmon", plan: 220 },
  { day: 2, item: "thighs", plan: 700 },
  { day: 3, item: "ribs", plan: 450 },
  { day: 4, item: "beef", plan: 380 },
];
/** This week's final counts, against the plan: what the figure plays. */
export const THIS_WEEK: readonly number[] = [334, 205, 741, 431, 396];
/** How far a final count moves from the plan, one standard deviation, for the many weeks. */
export const COUNT_SD = 0.1;

export type Trucks = 2 | 3 | 5;
/** The days the distributor's truck comes, for 2, 3 or 5 trucks a week. */
export const TRUCK_DAYS: Record<Trucks, readonly number[]> = { 2: [0, 3], 3: [0, 2, 4], 5: [0, 1, 2, 3, 4] };

/** The case places on the shelves the meat and fish have, in the walk-in and in the freezer: three shelves, two
 *  units of two cases each along them, a case stacked on each. */
export const SHELVES = 3;
export const COLUMNS = 4;
export const LEVELS = 2;
export const WALKIN_PLACES = SHELVES * COLUMNS * LEVELS;
export const FREEZER_PLACES = SHELVES * COLUMNS * LEVELS;

export interface Policy {
  trucks: Trucks;
  frozen: boolean;
  /** Ordered over the plan, as a fraction (0.1 = 10%). */
  cushion: number;
}

export type Where = "walkin" | "freezer";
export type State = "fresh" | "frozen" | "thawing" | "thawed";
export interface Lot {
  id: number;
  item: ItemKey;
  lb: number;
  where: Where;
  state: State;
  /** Absolute day it came in (or was pulled to thaw). */
  inDay: number;
  /** Absolute day it can first be used (a thawing bag the day after it is pulled). */
  readyDay: number;
  /** Absolute last day it can be used; Infinity frozen. */
  lastDay: number;
}
/** Cases a lot fills on a shelf: whole cases, an opened one counted whole. */
export const casesOf = (l: { item: ItemKey; lb: number }): number => Math.ceil(l.lb / ITEM[l.item].caseLb - 1e-9);

export type EventKind = "truck" | "thaw" | "prep" | "check";
/** What happened at one moment of the week, and the stock just after it. Hours from Monday 00:00 of the week played. */
export interface WeekEvent {
  hour: number;
  day: number;
  kind: EventKind;
  /** Lots that came in, were pulled to thaw (the new thawing lot), were used from, or were thrown out. */
  lots: { id: number; item: ItemKey; lb: number }[];
  /** Prep: the lunches of the day, the lb used, lunches short and the lb bought locally. */
  lunches?: number;
  item?: ItemKey;
  shortLunches?: number;
  stock: Lot[];
}

export interface WeekResult {
  /** The stock on Monday at midnight, carried in from the week before. */
  start: Lot[];
  policy: Policy;
  /** Thrown out over the week, lb, in all and by item. */
  thrownLb: number;
  thrownBy: Record<ItemKey, number>;
  /** Lunches whose meat was not there at prep. */
  shortLunches: number;
  shortBy: Record<ItemKey, number>;
  /** The fullest the walk-in's and freezer's places got over the week, in cases (more than the places: on the floor). */
  walkinPeak: number;
  freezerPeak: number;
  /** Bought over the week, lb. */
  boughtLb: number;
  trucks: number;
  events: WeekEvent[];
}

const zero = (): Record<ItemKey, number> => ({ salmon: 0, ribs: 0, beef: 0, thighs: 0 });
const menuOn = (d: number) => MENU.find((m) => m.day === ((d % 7) + 7) % 7) ?? null;

/**
 * Runs three weeks with these final counts (lunches by absolute day, 0 at the weekend) and reports the middle one, its
 * stock carried in from the first. `counts` has 21 entries.
 */
export function runWeeks(p: Policy, counts: readonly number[]): WeekResult {
  const W0 = 7;
  const W1 = 14;
  const truckOn = (d: number) => TRUCK_DAYS[p.trucks].includes(((d % 7) + 7) % 7);
  const nextTruck = (d: number) => {
    for (let k = d + 1; k < d + 15; k++) if (truckOn(k)) return k;
    return d + 7;
  };
  const planLb = (u: number) => {
    const m = menuOn(u);
    return m ? m.plan * ITEM[m.item].lbPer * (1 + p.cushion) : 0;
  };
  let lots: Lot[] = [];
  let nextId = 1;
  const r: WeekResult = { start: [], policy: p, thrownLb: 0, thrownBy: zero(), shortLunches: 0, shortBy: zero(), walkinPeak: 0, freezerPeak: 0, boughtLb: 0, trucks: 0, events: [] };
  const inWeek = (d: number) => d >= W0 && d < W1;
  const snap = () => lots.map((l) => ({ ...l }));
  const log = (d: number, hour: number, kind: EventKind, ev: Omit<WeekEvent, "hour" | "day" | "kind" | "stock">) => {
    if (!inWeek(d)) return;
    r.events.push({ hour: (d - W0) * 24 + hour, day: d - W0, kind, ...ev, stock: snap() });
  };
  const peaks = (d: number) => {
    if (!inWeek(d)) return;
    let wi = 0;
    let fz = 0;
    for (const l of lots) {
      if (l.lb <= 1e-9) continue;
      if (l.where === "walkin") wi += casesOf(l);
      else fz += casesOf(l);
    }
    r.walkinPeak = Math.max(r.walkinPeak, wi);
    r.freezerPeak = Math.max(r.freezerPeak, fz);
  };
  /** lb of `item` still unmet for these days' needs (ascending), once the lots on hand are given to them first-expiring
   *  first, each only to days it is ready for and still good on. */
  const unmet = (item: ItemKey, need: { u: number; lb: number }[], pool: Lot[]) => {
    const left = pool.filter((l) => l.item === item && l.lb > 1e-9).map((l) => ({ lb: l.lb, ready: l.readyDay, last: l.lastDay })).sort((a, b) => a.last - b.last);
    let short = 0;
    for (const n of need) {
      let want = n.lb;
      for (const l of left) {
        if (want <= 1e-9) break;
        if (l.lb <= 1e-9 || l.ready > n.u || l.last < n.u) continue;
        const take = Math.min(l.lb, want);
        l.lb -= take;
        want -= take;
      }
      short += want;
    }
    return short;
  };

  for (let d = 0; d < 21; d++) {
    if (d === W0) r.start = snap();
    // 06:00 the truck, its order placed by 2 pm yesterday: today's final count, the plan beyond.
    if (truckOn(d)) {
      if (inWeek(d)) r.trucks++;
      const next = nextTruck(d);
      const came: WeekEvent["lots"] = [];
      for (const it of ITEMS) {
        if (!p.frozen) {
          // Fresh: the lunches until the next truck that it is still good for.
          const need: { u: number; lb: number }[] = [];
          for (let u = d; u < next && u <= d + it.freshDays - 1; u++) {
            const m = menuOn(u);
            if (m?.item !== it.key) continue;
            need.push({ u, lb: u === d ? counts[u]! * it.lbPer : planLb(u) });
          }
          if (!need.length) continue;
          const lb = unmet(it.key, need, lots.filter((l) => l.where === "walkin"));
          if (lb <= 1e-6) continue;
          const cases = Math.ceil(lb / it.caseLb - 1e-9);
          const lot: Lot = { id: nextId++, item: it.key, lb: cases * it.caseLb, where: "walkin", state: "fresh", inDay: d, readyDay: d, lastDay: d + it.freshDays - 1 };
          lots.push(lot);
          came.push({ id: lot.id, item: lot.item, lb: lot.lb });
        } else {
          // Frozen: what the pulls from this morning until the next truck's will take (each its thawing time ahead of
          // its lunch), by the plan, less what is in the freezer.
          let want = 0;
          for (let x = d; x < next; x++) if (menuOn(x + it.thawDays)?.item === it.key) want += planLb(x + it.thawDays);
          const have = lots.filter((l) => l.item === it.key && l.state === "frozen").reduce((s, l) => s + l.lb, 0);
          const lb = want - have;
          if (lb <= 1e-6) continue;
          const cases = Math.ceil(lb / it.caseLb - 1e-9);
          const lot: Lot = { id: nextId++, item: it.key, lb: cases * it.caseLb, where: "freezer", state: "frozen", inDay: d, readyDay: d, lastDay: Infinity };
          lots.push(lot);
          came.push({ id: lot.id, item: lot.item, lb: lot.lb });
        }
      }
      for (const c of came) if (inWeek(d)) r.boughtLb += c.lb;
      log(d, 6, "truck", { lots: came });
      peaks(d);
    }
    // 08:00 frozen: each lunch its thawing time ahead, by the plan, pulled to thaw in the walk-in, less what is already
    // thawing or thawed and still good that day; by the case, first in first out.
    if (p.frozen) {
      const pulled: WeekEvent["lots"] = [];
      for (const it of ITEMS) {
        const u = d + it.thawDays;
        if (menuOn(u)?.item !== it.key) continue;
        const lb = unmet(it.key, [{ u, lb: planLb(u) }], lots.filter((l) => l.state === "thawing" || l.state === "thawed"));
        let pull = lb > 1e-6 ? Math.ceil(lb / it.caseLb - 1e-9) * it.caseLb : 0;
        if (pull <= 0) continue;
        let got = 0;
        for (const l of lots.filter((x) => x.item === it.key && x.state === "frozen").sort((a, b) => a.inDay - b.inDay)) {
          const take = Math.min(l.lb, pull - got);
          l.lb -= take;
          got += take;
          if (got >= pull - 1e-9) break;
        }
        pull = got;
        if (pull <= 1e-9) continue;
        const lot: Lot = { id: nextId++, item: it.key, lb: pull, where: "walkin", state: "thawing", inDay: d, readyDay: u, lastDay: u + it.thawedDays - 1 };
        lots.push(lot);
        pulled.push({ id: lot.id, item: lot.item, lb: lot.lb });
      }
      if (pulled.length) {
        lots = lots.filter((l) => l.lb > 1e-9);
        log(d, 8, "thaw", { lots: pulled });
        peaks(d);
      }
    }
    // Thawed by morning.
    for (const l of lots) if (l.state === "thawing" && l.readyDay <= d) l.state = "thawed";
    // 09:00 prep: today's final count, first-expiring first; what is not there is bought on the day.
    const m = menuOn(d);
    if (m && counts[d]! > 0) {
      const it = ITEM[m.item];
      let want = counts[d]! * it.lbPer;
      const used: WeekEvent["lots"] = [];
      for (const l of lots.filter((x) => x.item === it.key && x.where === "walkin" && x.readyDay <= d && x.lastDay >= d).sort((a, b) => a.lastDay - b.lastDay || a.inDay - b.inDay)) {
        if (want <= 1e-9) break;
        const take = Math.min(l.lb, want);
        l.lb -= take;
        want -= take;
        used.push({ id: l.id, item: l.item, lb: take });
      }
      const short = want > 1e-6 ? Math.ceil(want / it.lbPer - 1e-9) : 0;
      lots = lots.filter((l) => l.lb > 1e-9);
      if (inWeek(d)) {
        r.shortLunches += short;
        r.shortBy[it.key] += short;
      }
      log(d, 9, "prep", { lots: used, lunches: counts[d]!, item: it.key, shortLunches: short });
    }
    // 22:00 the night's check: whatever is still there on its last day goes.
    const out = lots.filter((l) => l.lastDay <= d && l.where === "walkin" && l.lb > 1e-9);
    if (out.length) {
      if (inWeek(d))
        for (const l of out) {
          r.thrownLb += l.lb;
          r.thrownBy[l.item] += l.lb;
        }
      lots = lots.filter((l) => !out.includes(l));
      log(d, 22, "check", { lots: out.map((l) => ({ id: l.id, item: l.item, lb: l.lb })) });
    }
  }
  return r;
}

/** Final counts for three weeks: this week's in the middle, the weeks either side as planned. */
export function thisWeekCounts(): number[] {
  const c: number[] = [];
  for (let d = 0; d < 21; d++) {
    const m = menuOn(d);
    c.push(m ? (d >= 7 && d < 14 ? THIS_WEEK[m.day]! : m.plan) : 0);
  }
  return c;
}

/** A deterministic generator (mulberry32) and a standard normal from it. */
function rng(seed: number): () => number {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
function normal(r: () => number): number {
  const u = Math.max(1e-12, r());
  return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * r());
}

export interface ManyWeeks {
  /** Medians over the weeks, and the 9th decile (1 week in 10 is worse). */
  thrownLb: { mid: number; p90: number };
  shortLunches: { mid: number; p90: number };
  walkinPeak: { mid: number; p90: number };
  freezerPeak: { mid: number; p90: number };
  /** Share of weeks with any lunch short, and with any case on the floor (the walk-in's, the freezer's). */
  anyShort: number;
  anyFloor: { walkin: number; freezer: number };
  boughtLb: number;
  thrownBy: Record<ItemKey, number>;
  shortBy: Record<ItemKey, number>;
}
const quant = (xs: number[], q: number) => {
  const s = [...xs].sort((a, b) => a - b);
  const i = (s.length - 1) * q;
  const lo = Math.floor(i);
  return s[lo]! + (s[Math.min(s.length - 1, lo + 1)]! - s[lo]!) * (i - lo);
};

/** The same policy over `n` weeks of final counts drawn around the plan (the same weeks for every policy). */
export function manyWeeks(p: Policy, n = 400, seed = 7): ManyWeeks {
  const run = manyWeeksRun(p, n, seed);
  run.step(n);
  return run.result();
}

/** manyWeeks a slice at a time: `step(k)` runs up to k more of the weeks and says whether all n are done; `result()`
 *  then gives what manyWeeks would. */
export function manyWeeksRun(p: Policy, n = 400, seed = 7): { step: (k: number) => boolean; result: () => ManyWeeks } {
  const g = rng(seed);
  const rows: WeekResult[] = [];
  const floors: Layout["floorMost"][] = [];
  const step = (k: number) => {
    for (let i = 0; i < k && rows.length < n; i++) {
      const counts: number[] = [];
      for (let d = 0; d < 21; d++) {
        const m = menuOn(d);
        counts.push(m ? Math.max(0, Math.round(m.plan * (1 + COUNT_SD * normal(g)))) : 0);
      }
      const w = runWeeks(p, counts);
      rows.push(w);
      floors.push(layoutWeek(w).floorMost);
    }
    return rows.length >= n;
  };
  return { step, result: () => summarise(rows, floors) };
}

function summarise(rows: WeekResult[], floors: Layout["floorMost"][]): ManyWeeks {
  const n = rows.length;
  const col = (f: (w: WeekResult) => number) => rows.map(f);
  const both = (xs: number[]) => ({ mid: quant(xs, 0.5), p90: quant(xs, 0.9) });
  const mean = (f: (w: WeekResult) => number) => rows.reduce((s, w) => s + f(w), 0) / rows.length;
  const by = (f: (w: WeekResult) => Record<ItemKey, number>) => Object.fromEntries(ITEMS.map((i) => [i.key, mean((w) => f(w)[i.key])])) as Record<ItemKey, number>;
  return {
    thrownLb: both(col((w) => w.thrownLb)),
    shortLunches: both(col((w) => w.shortLunches)),
    walkinPeak: both(col((w) => w.walkinPeak)),
    freezerPeak: both(col((w) => w.freezerPeak)),
    anyShort: mean((w) => (w.shortLunches > 0 ? 1 : 0)),
    anyFloor: { walkin: floors.filter((f) => f.walkin > 0).length / n, freezer: floors.filter((f) => f.freezer > 0).length / n },
    boughtLb: mean((w) => w.boughtLb),
    thrownBy: by((w) => w.thrownBy),
    shortBy: by((w) => w.shortBy),
  };
}

/* ------------------------------------------------------------------------------------------------------------ */
/* The cases, for the figure: where each one stands, and how each is carried, event by event.                    */
/* ------------------------------------------------------------------------------------------------------------ */

/** A case's place. On a shelf: its shelf (tier 0 the lowest), its column (0 nearest the door), and level 0 on the
 *  shelf or 1 on a case of its own lot there. On the floor: its row place in `col`, and its level. */
export interface Spot {
  box: Where;
  floor: boolean;
  tier: number;
  col: number;
  level: number;
}
/** Floor places in each box, clear of the way in, each a stack of one delivery's cases up to FLOOR_LEVELS high. The
 *  figure's week never needs more; past them a case takes the next place number all the same. Frozen cases stand
 *  4 high, 76 cm; the walk-in's, 57 cm. */
export const FLOOR_ROWS: Record<Where, number> = { walkin: 3, freezer: 2 };
export const FLOOR_LEVELS: Record<Where, number> = { walkin: 3, freezer: 4 };
/** Where a carried case comes from or goes: a place, the truck at the dock, or the kitchen (the day's prep, the bin). */
export type End = Spot | "dock" | "kitchen";
export interface Move {
  caseId: number;
  lot: number;
  item: ItemKey;
  from: End;
  to: End;
}
export interface CaseAt {
  caseId: number;
  lot: number;
  item: ItemKey;
  at: Spot;
}
export interface Staged {
  ev: WeekEvent;
  /** The cases carried, in the order they go. */
  moves: Move[];
  /** Every case on hand once they have. */
  after: CaseAt[];
}
export interface Layout {
  start: CaseAt[];
  events: Staged[];
  /** Cases put on the floor though a shelf place was free: none there kept the hottest-cooked lowest, even with cases
   *  moved down to make room. */
  shelvedOut: number;
  /** The most cases on each floor at once. */
  floorMost: Record<Where, number>;
}
/** Each item's place in the shelves' order, 0 the highest (ITEMS is top to bottom). */
const RANK = Object.fromEntries(ITEMS.map((i, n) => [i.key, n])) as Record<ItemKey, number>;
const spotKey = (s: Spot) => `${s.box}|${s.floor ? 1 : 0}|${s.tier}|${s.col}|${s.level}`;

/**
 * Where each case goes. A case goes on top of a case of its own lot if one has room, else on a free shelf place
 * where the shelves' order holds in its column, top to bottom: salmon, short ribs, ground beef, chicken (ITEMS'
 * order), each as near the door as it can. When no free place keeps that order, the stacks just over a free place are
 * moved down a shelf each, or the stacks just under it up a shelf each (moves of their own, shelf to shelf), and the
 * case goes where the last of them stood. With no shelf place left that keeps the order, it stands on the floor, on
 * its own lot's stack there if one has room. Cases leave a lot from the floor first, then from half stacks, the top of
 * a stack first. A case pulled from the freezer to thaw is the same case, put in the walk-in.
 */
export function layoutWeek(w: WeekResult): Layout {
  const at = new Map<number, CaseAt>();
  const where = new Map<string, number>();
  const byLot = new Map<number, number[]>();
  let nextCase = 1;
  let shelvedOut = 0;
  const floorMost: Record<Where, number> = { walkin: 0, freezer: 0 };
  const snapshot = () => [...at.values()].map((c) => ({ ...c, at: { ...c.at } }));
  const who = (s: Spot) => {
    const id = where.get(spotKey(s));
    return id === undefined ? undefined : at.get(id);
  };
  const put = (c: CaseAt, s: Spot) => {
    where.delete(spotKey(c.at));
    c.at = { ...s };
    where.set(spotKey(s), c.caseId);
  };
  /** The moves made to make room, waiting to go in the event's list before the case they made room for. */
  let restacks: Move[] = [];
  function place(box: Where, item: ItemKey, lot: number): Spot {
    const k = RANK[item];
    const mine = (byLot.get(lot) ?? []).map((id) => at.get(id)!).filter((c) => c.at.box === box);
    // On top of its own lot's case: off the floor first.
    const onTop = mine
      .map((c) => ({ ...c.at, level: c.at.level + 1 }))
      .filter((s) => s.level < (s.floor ? FLOOR_LEVELS[box] : LEVELS) && !who(s))
      .sort((a, b) => Number(a.floor) - Number(b.floor) || a.col - b.col || a.tier - b.tier || a.level - b.level);
    if (onTop.length) return onTop[0]!;
    const here = [...at.values()].filter((c) => c.at.box === box && !c.at.floor);
    /** Whether the case may stand at `s`, the cases `under` and `over` (being moved) to end up below and above it. */
    const fits = (s: Spot, under: CaseAt[] = [], over: CaseAt[] = []) =>
      here.every((c) => {
        if (c.at.col !== s.col) return true;
        if (under.includes(c)) return RANK[c.item] >= k;
        if (over.includes(c)) return RANK[c.item] <= k;
        return c.at.tier < s.tier ? RANK[c.item] >= k : c.at.tier > s.tier ? RANK[c.item] <= k : true;
      });
    const free: Spot[] = [];
    for (let tier = 0; tier < SHELVES; tier++) for (let col = 0; col < COLUMNS; col++) if (!who({ box, floor: false, tier, col, level: 0 })) free.push({ box, floor: false, tier, col, level: 0 });
    const good = free.filter((s) => fits(s));
    // Ground beef and chicken as low as they go, salmon and short ribs as high.
    const low = ITEM[item].cookF >= 155;
    const order = (a: Spot, b: Spot) => (low ? a.tier - b.tier : b.tier - a.tier) || a.col - b.col;
    if (good.length) return good.sort(order)[0]!;
    // Make room in a column with a free place: move the stacks just over it down a shelf each (the case to go above
    // them), or the stacks just under it up a shelf each (the case to go below them), the nearest first, and put the
    // case where the last of them stood, if it fits there. The fewest stacks that do, down before up.
    let best: { at: Spot; block: CaseAt[][]; step: -1 | 1 } | null = null;
    for (const step of [-1, 1] as const)
      for (const s of free) {
        const block: CaseAt[][] = [];
        for (let t = s.tier - step; t >= 0 && t < SHELVES; t -= step) {
          const a = who({ ...s, tier: t, level: 0 });
          if (!a) break;
          const b = who({ ...s, tier: t, level: 1 });
          block.push(b ? [a, b] : [a]);
          if (best && best.block.length <= block.length) break;
          const spot: Spot = { ...s, tier: s.tier - step * block.length };
          if (step < 0 ? fits(spot, block.flat()) : fits(spot, [], block.flat())) {
            best = { at: spot, block: [...block], step };
            break;
          }
        }
      }
    if (best) {
      for (const [a, b] of best.block) {
        const to: Spot = { ...a!.at, tier: a!.at.tier + best.step, level: 0 };
        // The top case first, onto the next shelf; the one under it then onto that.
        const first = b ?? a!;
        const from1 = { ...first.at };
        put(first, to);
        restacks.push({ caseId: first.caseId, lot: first.lot, item: first.item, from: from1, to: { ...to } });
        if (b) {
          const from2 = { ...a!.at };
          const on: Spot = { ...to, level: 1 };
          put(a!, on);
          restacks.push({ caseId: a!.caseId, lot: a!.lot, item: a!.item, from: from2, to: { ...on } });
        }
      }
      return best.at;
    }
    // The floor (a case of its lot there is stacked on, above, when it has room): no shelf place takes it.
    if (free.length) shelvedOut++;
    let i = 0;
    while (who({ box, floor: true, tier: 0, col: i, level: 0 })) i++;
    const s: Spot = { box, floor: true, tier: 0, col: i, level: 0 };
    return s;
  }
  /** `n` of a lot's cases: off the floor first, from the top of a stack first, then nearest the door. */
  function take(lot: number, n: number): number[] {
    const ids = byLot.get(lot) ?? [];
    // Place by place, the top case first: the floor's first, then half stacks, so what is left stands in whole stacks.
    const placeOf = (s: Spot) => spotKey({ ...s, level: 0 });
    const inPlace = new Map<string, number>();
    for (const id of ids) inPlace.set(placeOf(at.get(id)!.at), (inPlace.get(placeOf(at.get(id)!.at)) ?? 0) + 1);
    const order = [...ids].sort((a, b) => {
      const A = at.get(a)!.at;
      const B = at.get(b)!.at;
      return Number(B.floor) - Number(A.floor) || inPlace.get(placeOf(A))! - inPlace.get(placeOf(B))! || A.col - B.col || B.tier - A.tier || B.level - A.level;
    });
    const out = order.slice(0, Math.max(0, n));
    byLot.set(lot, ids.filter((i) => !out.includes(i)));
    for (const id of out) where.delete(spotKey(at.get(id)!.at));
    return out;
  }
  const add = (caseId: number, lot: Lot, spot: Spot) => {
    const c: CaseAt = { caseId, lot: lot.id, item: lot.item, at: { ...spot } };
    at.set(caseId, c);
    where.set(spotKey(spot), caseId);
    if (!byLot.has(lot.id)) byLot.set(lot.id, []);
    byLot.get(lot.id)!.push(caseId);
    if (spot.floor) floorMost[spot.box] = Math.max(floorMost[spot.box], [...at.values()].filter((x) => x.at.box === spot.box && x.at.floor).length);
  };
  for (const l of w.start) for (let k = 0; k < casesOf(l); k++) add(nextCase++, l, place(l.where, l.item, l.id));
  // Moves made before the week starts are not shown.
  restacks = [];
  const start = snapshot();
  let prev = new Map(w.start.map((l) => [l.id, l]));
  const events: Staged[] = [];
  for (const ev of w.events) {
    const now = new Map(ev.stock.map((l) => [l.id, l]));
    const moves: Move[] = [];
    const left = (id: number) => {
      const before = prev.get(id);
      const after = now.get(id);
      return (before ? casesOf(before) : 0) - (after ? casesOf(after) : 0);
    };
    const arrive = (id: number, lot: Lot, from: End) => {
      const to = place(lot.where, lot.item, lot.id);
      moves.push(...restacks);
      restacks = [];
      add(id, lot, to);
      moves.push({ caseId: id, lot: lot.id, item: lot.item, from, to });
    };
    if (ev.kind === "truck") {
      // Chicken first, on the lowest shelves, the rest above it as it comes off.
      for (const x of [...ev.lots].sort((a, b) => ITEM[b.item].cookF - ITEM[a.item].cookF)) {
        const lot = now.get(x.id)!;
        for (let k = 0; k < casesOf(lot); k++) arrive(nextCase++, lot, "dock");
      }
    } else if (ev.kind === "thaw") {
      for (const x of ev.lots) {
        const lot = now.get(x.id)!;
        for (const [id, l] of prev) {
          if (l.where !== "freezer" || l.item !== lot.item) continue;
          for (const cid of take(id, left(id))) {
            const from = { ...at.get(cid)!.at };
            at.delete(cid);
            arrive(cid, lot, from);
          }
        }
      }
    } else {
      // Prep takes the cases it empties; the night's check, every case of a lot past its day.
      for (const x of ev.lots) {
        for (const cid of take(x.id, left(x.id))) {
          const c = at.get(cid)!;
          at.delete(cid);
          moves.push({ caseId: cid, lot: x.id, item: c.item, from: c.at, to: "kitchen" });
        }
      }
    }
    prev = now;
    events.push({ ev, moves, after: snapshot() });
  }
  return { start, events, shelvedOut, floorMost };
}
