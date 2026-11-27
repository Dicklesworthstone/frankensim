/**
 * A combi trips (chapter 07, fifth figure): the morning above with three cooks, planned as usual, then the same
 * morning with one of its two combis out of use for a while and the kitchen keeping to the plan (carryOn). For every
 * trip on a grid (when it trips, every ten minutes from 6:00 to 10:20, and how long it is out, ten to ninety minutes)
 * and for each combi: how much later the van leaves, and how many more portions are past their best at lunch.
 */
import { HORIZON, OVENS, carryOn, type Morning, type MorningSettings } from "./cateringMorning";

/** When a trip can start, minutes from 6:00: 6:00 to 10:20, every ten minutes. */
export const TRIP_AT = { from: 0, to: 260, step: 10 } as const;
/** How long a combi can be out, minutes. */
export const TRIP_FOR = { from: 10, to: 90, step: 10 } as const;
export const AT_N = (TRIP_AT.to - TRIP_AT.from) / TRIP_AT.step + 1;
export const FOR_N = (TRIP_FOR.to - TRIP_FOR.from) / TRIP_FOR.step + 1;
/** The crew for every trip: the morning's three cooks. */
export const TRIP_COOKS = 3;

export interface TripSettings {
  menu: number;
  meals: number;
}
export const tripKey = (s: TripSettings): string => `${s.menu}/${s.meals}`;
/** The plan's settings: three cooks, nothing cooked the day before, the kitchen as built. */
export const tripBase = (s: TripSettings): MorningSettings => ({ meals: s.meals, cooks: TRIP_COOKS, menu: s.menu, ahead: false });

/** A trip: which combi, when it trips and for how long, minutes. */
export interface Trip {
  oven: number;
  at: number;
  mins: number;
}

/** The morning `plan` with this trip, the kitchen keeping to the plan. */
export function tripMorning(plan: Morning, trip: Trip): Morning {
  return carryOn({ ...plan.settings, out: [{ oven: trip.oven, t0: trip.at, t1: trip.at + trip.mins }] }, plan);
}

/** What a trip costs against the plan: minutes later the van leaves (Infinity: not by 4 pm), portions more past their
 *  best at lunch. A list schedule can come out ahead of its plan by a little; that counts as costing nothing. */
export interface TripCost {
  late: number;
  pastBest: number;
}
export function costOf(m: Morning, plan: Morning): TripCost {
  if (!m.finished) return { late: Infinity, pastBest: 0 };
  return { late: Math.max(0, m.depart - Math.min(plan.depart, HORIZON)), pastBest: Math.max(0, m.pastBest - plan.pastBest) };
}
/** Whether a trip costs nothing: the van no later (to the half minute), and nothing more past its best. */
export const free = (c: TripCost): boolean => c.late <= 0.5 && c.pastBest === 0;

/** Each combi's grid of costs: index [row * AT_N + col], row the time out (TRIP_FOR), col when it trips (TRIP_AT). */
export interface TripMap {
  late: Float32Array[];
  pastBest: Int32Array[];
}
export const cellOf = (at: number, dur: number): number => Math.round((dur - TRIP_FOR.from) / TRIP_FOR.step) * AT_N + Math.round((at - TRIP_AT.from) / TRIP_AT.step);

export function tripMap(plan: Morning): TripMap {
  const late: Float32Array[] = [];
  const pastBest: Int32Array[] = [];
  for (let oven = 0; oven < OVENS; oven++) {
    const L = new Float32Array(AT_N * FOR_N);
    const P = new Int32Array(AT_N * FOR_N);
    for (let r = 0; r < FOR_N; r++)
      for (let c = 0; c < AT_N; c++) {
        const k = r * AT_N + c;
        const cost = costOf(tripMorning(plan, { oven, at: TRIP_AT.from + c * TRIP_AT.step, mins: TRIP_FOR.from + r * TRIP_FOR.step }), plan);
        L[k] = cost.late;
        P[k] = cost.pastBest;
      }
    late.push(L);
    pastBest.push(P);
  }
  return { late, pastBest };
}
