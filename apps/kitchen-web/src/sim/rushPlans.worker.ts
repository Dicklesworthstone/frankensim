/**
 * Every plan for one lunch service, off the page's main thread: which dishes are held (any of the five), whether a
 * free cook plates for the others, and where the thighs cook. Each is run over the same simulated days as the figure
 * (runRush), and its medians are posted back as they come.
 */
import { PLANS, modesFor, planSettings, runRush, type RushDay } from "./serviceRush";

export interface PlanAsk {
  key: string;
  day: RushDay;
  sigma: number;
  days: number;
}
export interface PlanScore {
  key: string;
  /** The plan's number (see planSettings). */
  plan: number;
  waitP90: number;
  waste: number;
  prepLoad: number;
}

const scope = self as unknown as { onmessage: ((e: MessageEvent<PlanAsk>) => void) | null; postMessage(m: PlanScore): void };
scope.onmessage = (e) => {
  const { key, day, sigma, days } = e.data;
  for (let plan = 0; plan < PLANS; plan++) {
    const p = planSettings(plan);
    const r = runRush({ modes: modesFor(p.hold), day, sigma, pool: p.pool, menu: p.menu }, days);
    scope.postMessage({ key, plan, waitP90: r.waitP90.p50, waste: r.waste.p50, prepLoad: r.prepLoad.p50 });
  }
};
