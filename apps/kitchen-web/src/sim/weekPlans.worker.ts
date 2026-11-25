/**
 * Every plan for the week, off the page's main thread: two, three or four cooks; cooking on the day or ahead; one to
 * five racks in the walk-in; any of the three dishes swapped for the variant that keeps. Each runs over the same
 * simulated weeks as the figure (runCommissary), on the figure's week (reshapeWeek), and its scores are posted back
 * as they come.
 */
import { WEEK_PLANS, paidHours, reshapeWeek, runCommissary, weekPlan, withCooks } from "./commissaryWeek";

export interface WeekPlanAsk {
  key: string;
  sigma: number;
  /** The swing between quiet and busy days (see reshapeWeek). */
  swing: number;
}
export interface WeekPlanScore {
  key: string;
  /** The plan's number (see weekPlan). */
  plan: number;
  /** Crew-hours paid a week (median week). */
  paid: number;
  /** Meals out late: the median week, and 1 week in 10. */
  late: number;
  lateP90: number;
  /** Portions thrown out, median week. */
  waste: number;
}

const scope = self as unknown as { onmessage: ((e: MessageEvent<WeekPlanAsk>) => void) | null; postMessage(m: WeekPlanScore): void };
scope.onmessage = (e) => {
  const { key, sigma, swing } = e.data;
  const week = reshapeWeek(swing);
  for (let plan = 0; plan < WEEK_PLANS; plan++) {
    const p = weekPlan(plan);
    const r = runCommissary({ strategy: p.strategy, swaps: p.swaps, sigma, racks: p.racks }, withCooks(p.cooks, week));
    scope.postMessage({ key, plan, paid: paidHours(r, p.cooks, week), late: r.late.p50, lateP90: r.late.p90, waste: r.wasted.p50 });
  }
};
