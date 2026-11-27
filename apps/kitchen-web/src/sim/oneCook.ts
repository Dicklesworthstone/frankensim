/**
 * One cook, two dishes: what camera glasses on a cook let the twin keep for each dish, and how it plans two of them
 * for one pair of hands.
 *
 * A dish's card is its stages as the morning's model has them (sim/cateringMorning.ts): where each happens, how long
 * it takes and how much of that time it needs the cook's hands. A stage at the table needs them all the time; a pot
 * stirred throughout needs them all the time, a pan only minded half (one cook minds two); a combi needs them only to
 * load and unload it, and nobody while its program runs.
 *
 * `oneCook` plans the two dishes for one cook as soon as the work allows (`soonest`: the morning's cooks, stations and
 * rules, with no holding time and no van), and each dish alone, the same way; one after the other is the first alone,
 * then the second alone from the minute the first is done.
 */
import { DISHES, soonest, type Job, type Morning, type MorningDish, type Station } from "./cateringMorning";

/** The dishes a cook can be given, in the menu card's order. */
export const COOK_DISHES: readonly MorningDish[] = [DISHES.piccata, DISHES.thighs, DISHES.risotto, DISHES.mac, DISHES.tacos, DISHES.braise];

/** One row of a dish's card. */
export interface CardStep {
  verb: string;
  at: Station;
  /** Minutes: per portion at the table, per pot or pan on the range, per program in a combi. */
  min: number;
  per: "portion" | "batch" | "load";
  /** The share of the stage's minutes it needs the cook's hands, for one batch of the dish. */
  hands: number;
  /** Oven: hands-on minutes per pan, to load it, and again to unload it. */
  perPan?: number;
}

/** A dish's card: its stages for one batch. */
export function cardOf(d: MorningDish): CardStep[] {
  const pans = Math.ceil(d.batch / d.pan);
  return d.stages.map((s) => {
    if (s.at === "table") return { verb: s.verb, at: s.at, min: s.min, per: "portion", hands: 1 };
    if (s.at === "range") return { verb: s.verb, at: s.at, min: s.min, per: "batch", hands: s.hands ?? 1 };
    const io = 2 * pans * (s.perPan ?? 0);
    return { verb: s.verb, at: s.at, min: s.min, per: "load", hands: io / (s.min + io), perPan: s.perPan ?? 0 };
  });
}

/** When the last of a morning's work is done, minutes from its start. */
export const endOf = (m: Morning): number => m.batches.reduce((a, b) => Math.max(a, b.packed), 0);

/** One cook's plan for some dishes, as soon as the work allows. */
const plan = (dishes: readonly MorningDish[], portions: number): Morning =>
  soonest({ meals: portions * dishes.length, cooks: 1, menu: 0, ahead: false, dishes, portions: dishes.map(() => portions) });

/**
 * Two plans one after the other: `b`'s, from the minute `a`'s work is done, as the second dish (place 1) of the same
 * morning, with its batches numbered after `a`'s.
 */
export function oneAfterOther(a: Morning, b: Morning, dishes: readonly MorningDish[]): Morning {
  const at = endOf(a);
  const nb = a.batches.length;
  const job = (j: Job): Job => ({ ...j, batch: j.batch + nb, slot: 1, t0: j.t0 + at, t1: j.t1 + at });
  /** A time in `b`'s plan, moved on (-1, never, stays). */
  const on = (x: number) => (x < 0 ? x : x + at);
  return {
    ...a,
    settings: { ...a.settings, meals: a.settings.meals + b.settings.meals, dishes, portions: [...(a.settings.portions ?? []), ...(b.settings.portions ?? [])] },
    jobs: [...a.jobs, ...b.jobs.map(job)],
    burns: [...a.burns, ...b.burns.map((x) => ({ ...x, batch: x.batch + nb, slot: 1, t0: x.t0 + at, t1: x.t1 + at }))],
    loads: [...a.loads, ...b.loads.map((x) => ({ ...x, batches: x.batches.map((i) => i + nb), slot: 1, t0: x.t0 + at, t1: x.t1 + at, t2: x.t2 + at, tu: x.tu + at, t3: x.t3 + at }))],
    waits: [...a.waits, ...b.waits.map((w) => ({ ...w, t0: w.t0 + at, t1: w.t1 + at }))],
    batches: [...a.batches, ...b.batches.map((x) => ({ ...x, slot: 1, release: x.release + at, started: on(x.started), cooked: on(x.cooked), packed: on(x.packed) }))],
    depart: Math.max(a.depart, b.depart + at),
    late: Math.max(a.late, b.late + at),
    finished: a.finished && b.finished,
    pastBest: 0,
    pastBestBySlot: [0, 0],
    lastPacked: [a.lastPacked[0] ?? 0, (b.lastPacked[0] ?? 0) + at],
    onTimePastBest: null,
    pans: a.pans + b.pans,
    handsMin: a.handsMin + b.handsMin,
    aheadMin: 0,
  };
}

export interface OneCook {
  dishes: readonly [MorningDish, MorningDish];
  portions: number;
  /** Both dishes for one cook, planned together; and one after the other. */
  together: Morning;
  apart: Morning;
  /** When each is done: together, each alone, and one after the other. */
  end: number;
  endA: number;
  endB: number;
  endApart: number;
  /** Minutes the cook's hands are not needed, together (a pot only minded leaves half of them free). */
  free: number;
}

/** Two dishes, `portions` of each, for one cook: planned together, and one after the other. */
export function oneCook(a: MorningDish, b: MorningDish, portions: number): OneCook {
  const dishes = [a, b] as const;
  const together = plan(dishes, portions);
  const aloneA = plan([a], portions);
  const aloneB = plan([b], portions);
  const apart = oneAfterOther(aloneA, aloneB, dishes);
  const end = endOf(together);
  const endA = endOf(aloneA);
  const endB = endOf(aloneB);
  return { dishes, portions, together, apart, end, endA, endB, endApart: endA + endB, free: end - together.handsMin };
}
