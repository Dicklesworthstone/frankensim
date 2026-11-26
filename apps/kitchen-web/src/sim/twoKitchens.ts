/**
 * One order from two kitchens: the hub, and a second kitchen across town drafted for the morning.
 *
 * Both kitchens are the morning's commissary (cateringMorning.ts): a six-burner range, two combis, a worktable. The
 * hub's van leaves at 10:30 as before; the second kitchen is farther from the client, so its van leaves at 10:00, and
 * both arrive for lunch at 11:30. Each place on the menu (the sides too) is cooked all at the hub, half and half, or
 * all across town, and each kitchen plans its share of the morning on its own (planMorning, with its van's time).
 *
 * A dish cooked across town is held half an hour longer before lunch, so only a dish that keeps can go: one that has
 * to be cooked in the last minutes before the van has no time at all over there. Lunch is when the later van arrives;
 * a portion is past its best if it was cooked longer before lunch than its dish keeps.
 *
 * `bestPair` tries every split and keeps the best; `pairLimit` finds the most meals some split gets to lunch on time
 * with every portion at its best. Deterministic.
 */
import { DOOR, SERVE, madeIt, menuDishes, morningLimit, planMorning, portionsFor, type Crew, type Morning } from "./cateringMorning";

/** When the van from across town leaves: half an hour farther to drive. */
export const FAR_DOOR = DOOR - 30;
/** How much of a place's portions the kitchen across town cooks. */
export const SPLITS = [0, 0.5, 1] as const;

export interface PairSettings {
  meals: number;
  menu: number;
  /** Cooks at the hub, and across town. */
  hubCooks: Crew;
  farCooks: Crew;
}

export interface PairPlan {
  settings: PairSettings;
  /** Per place (the sides last), the share cooked across town. */
  parts: number[];
  /** Per place, the portions at the hub and across town. */
  here: number[];
  there: number[];
  hub: Morning;
  /** Null when nothing is cooked across town. */
  far: Morning | null;
  /** Both vans away on time and every portion at its best at lunch. */
  made: boolean;
  /** How late lunch is: the later van's lateness, minutes. */
  lunchLate: number;
  /** Portions past their best at that lunch, per place and in all. */
  pastBestBySlot: number[];
  pastBest: number;
  /** Crew-minutes of hands-on work at each kitchen. */
  hubHands: number;
  farHands: number;
}

/** What a plan and a limit are asked and answered by. */
export const planKey = (s: PairSettings): string => `${s.meals}|${s.menu}|${s.hubCooks}|${s.farCooks}`;
export const limitKey = (menu: number, hubCooks: Crew, farCooks: Crew): string => `${menu}|${hubCooks}|${farCooks}`;

/** A dish can go across town only if it keeps longer than the time from that van to lunch. */
export const travels = (holdMin: number): boolean => holdMin > SERVE - FAR_DOOR;

/** Each place's portions at the hub and across town for these shares (the hub takes the odd portion). */
export function splitPortions(meals: number, parts: readonly number[]): { here: number[]; there: number[] } {
  const all = portionsFor(meals);
  const there = all.map((p, k) => Math.floor(p * (parts[k] ?? 0)));
  return { here: all.map((p, k) => p - there[k]!), there };
}

/** Every split worth planning for a menu, the hub alone first: a dish that cannot travel stays at the hub. */
export function splitsFor(menu: number): number[][] {
  const can = menuDishes(menu).map((d) => travels(d.holdMin));
  let out: number[][] = [[]];
  for (let k = 0; k < can.length; k++) out = out.flatMap((c) => (can[k] ? SPLITS : [0]).map((f) => [...c, f]));
  return out;
}

/** Plans for one kitchen, keyed by what it cooks: splits share them. */
export type PlanCache = Map<string, Morning>;

function kitchenPlan(cache: PlanCache | null, meals: number, menu: number, cooks: Crew, portions: number[], door: number): Morning {
  const key = `${menu}|${cooks}|${door}|${portions.join(",")}`;
  let m = cache?.get(key);
  if (!m) {
    m = planMorning({ meals, cooks, menu, ahead: false, portions, door });
    cache?.set(key, m);
  }
  return m;
}

export function planPair(s: PairSettings, parts: readonly number[], cache: PlanCache | null = null): PairPlan {
  const { here, there } = splitPortions(s.meals, parts);
  const hub = kitchenPlan(cache, s.meals, s.menu, s.hubCooks, here, DOOR);
  const far = there.some((p) => p > 0) ? kitchenPlan(cache, s.meals, s.menu, s.farCooks, there, FAR_DOOR) : null;
  const lunchLate = Math.max(hub.late, far?.late ?? 0);
  const lunch = SERVE + lunchLate;
  const dishes = menuDishes(s.menu);
  const pastBestBySlot = dishes.map(() => 0);
  for (const m of far ? [hub, far] : [hub])
    for (const b of m.batches) if (b.cooked >= 0 && lunch - b.cooked > dishes[b.slot]!.holdMin + 1e-6) pastBestBySlot[b.slot]! += b.portions;
  const pastBest = pastBestBySlot.reduce((a, b) => a + b, 0);
  return {
    settings: { ...s },
    parts: [...parts],
    here,
    there,
    hub,
    far,
    made: madeIt(hub) && (!far || madeIt(far)) && lunchLate <= 0.5 && pastBest === 0,
    lunchLate,
    pastBestBySlot,
    pastBest,
    hubHands: hub.handsMin,
    farHands: far?.handsMin ?? 0,
  };
}

/** How close a kitchen came to missing its van: its last pan's time less the van's (negative: time to spare). */
const tightness = (m: Morning, door: number) => (m.finished ? Math.max(...m.lastPacked) - door : Infinity);

/**
 * The better of two plans. One that makes lunch at its best first, and of two that do, the one that carries less
 * across town (the kitchen across town is drafted only for what the hub cannot do itself). Of two that miss: the sooner
 * lunch, then fewer portions past their best. Then the more time to spare in the kitchen nearest to missing its van;
 * then less carried.
 */
export function better(a: PairPlan, b: PairPlan): boolean {
  const carried = (p: PairPlan) => p.there.reduce((x, y) => x + y, 0);
  if (a.made !== b.made) return a.made;
  if (a.made && carried(a) !== carried(b)) return carried(a) < carried(b);
  if (Math.abs(a.lunchLate - b.lunchLate) > 0.5) return a.lunchLate < b.lunchLate;
  if (a.pastBest !== b.pastBest) return a.pastBest < b.pastBest;
  const ta = Math.max(tightness(a.hub, DOOR), a.far ? tightness(a.far, FAR_DOOR) : -Infinity);
  const tb = Math.max(tightness(b.hub, DOOR), b.far ? tightness(b.far, FAR_DOOR) : -Infinity);
  if (Math.abs(ta - tb) > 0.5) return ta < tb;
  return carried(a) < carried(b);
}

/** The best of every split for these settings; with `alone`, the hub on its own. */
export function bestPair(s: PairSettings, alone = false, cache: PlanCache | null = null): PairPlan {
  const splits = alone ? [splitsFor(s.menu)[0]!] : splitsFor(s.menu);
  let best: PairPlan | null = null;
  for (const parts of splits) {
    const p = planPair(s, parts, cache);
    if (!best || better(p, best)) best = p;
  }
  return best!;
}

/**
 * The most meals, in tens, that some split gets to lunch at 11:30 with every portion at its best, that split, and the
 * hub alone's most (morningLimit, where the search starts). Each split is tried at ten more than the best found so far
 * and, when it makes that, raised by halving.
 */
export function pairLimit(s: Omit<PairSettings, "meals">, hi = 4000, cache: PlanCache | null = null): PairLimit {
  const steps = pairLimitSteps(s, hi, cache);
  for (;;) {
    const r = steps.next();
    if (r.done) return r.value;
  }
}

export interface PairLimit {
  meals: number;
  parts: number[];
  alone: number;
}

/** pairLimit, a split at a time: it stops after each, for a worker to answer something newer in between. */
export function* pairLimitSteps(s: Omit<PairSettings, "meals">, hi = 4000, cache: PlanCache | null = null): Generator<void, PairLimit, void> {
  const alone = morningLimit({ cooks: s.hubCooks, menu: s.menu, ahead: false });
  let best = { meals: alone, parts: splitsFor(s.menu)[0]! };
  const made = (meals: number, parts: number[]) => planPair({ ...s, meals }, parts, cache).made;
  for (const parts of splitsFor(s.menu).slice(1)) {
    yield;
    if (best.meals + 10 > hi || !made(best.meals + 10, parts)) continue;
    let lo = (best.meals + 10) / 10;
    let top = hi / 10 + 1;
    while (top - lo > 1) {
      const m = Math.floor((lo + top) / 2);
      if (made(m * 10, parts)) lo = m;
      else top = m;
    }
    best = { meals: lo * 10, parts };
  }
  return { ...best, alone };
}
