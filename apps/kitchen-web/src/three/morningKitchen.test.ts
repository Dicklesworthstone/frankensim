import { describe, expect, test } from "bun:test";
import { MENUS, planMorning } from "../sim/cateringMorning";
import { BUYS, buyKey, morningOf } from "../sim/buyOne";
import { END_HOLD_S, SEC_PER_MIN, planActs, planRoutines, type CookRig } from "./morningKitchen";
import { planningBody, type Cohort } from "./models/people";

/** Closer than this, metres between where two cooks stand, and their bodies meet: two at one range stand 0.46 apart. */
const BODY = 0.44;

describe("morningKitchen", () => {
  test("no cook stands in another: never two at one place, nor closer than a body's width, idle or at work", () => {
    for (const b of BUYS)
      for (let menu = 0; menu < MENUS; menu++)
        for (const meals of [200, 700, 1400]) {
          const m = planMorning(morningOf(meals, menu, b));
          const acts = planActs(m, m.settings.cooks, m.depart + 6);
          for (let c = 0; c < acts.length; c++)
            for (let d = c + 1; d < acts.length; d++)
              for (const x of acts[c]!)
                for (const y of acts[d]!) {
                  if (!(x.t0 < y.t1 - 1e-9 && y.t0 < x.t1 - 1e-9)) continue;
                  const gap = Math.hypot(x.spot.x - y.spot.x, x.spot.z - y.spot.z);
                  if (gap < BODY) throw new Error(`${buyKey(b)} menu ${menu} ${meals} meals: cooks ${c} and ${d} at ${x.spot.key} / ${y.spot.key}, ${gap.toFixed(2)} m apart, ${x.t0.toFixed(1)}-${x.t1.toFixed(1)} s / ${y.t0.toFixed(1)}-${y.t1.toFixed(1)} s`);
                }
        }
    expect(true).toBe(true);
  }, 60_000);

  test("the morning figure's crews too: two, three or four cooks, cooking ahead or not", () => {
    for (const cooks of [2, 3, 4] as const)
      for (const ahead of [false, true])
        for (let menu = 0; menu < MENUS; menu++)
          for (const meals of [200, 700, 1400]) {
            const m = planMorning({ meals, menu, cooks, ahead });
            const acts = planActs(m, cooks, m.depart + 6);
            for (let c = 0; c < acts.length; c++)
              for (let d = c + 1; d < acts.length; d++)
                for (const x of acts[c]!)
                  for (const y of acts[d]!)
                    if (x.t0 < y.t1 - 1e-9 && y.t0 < x.t1 - 1e-9 && Math.hypot(x.spot.x - y.spot.x, x.spot.z - y.spot.z) < BODY)
                      throw new Error(`${cooks} cooks${ahead ? ", ahead" : ""}, menu ${menu}, ${meals} meals: cooks ${c} and ${d} at ${x.spot.key} / ${y.spot.key}`);
          }
    expect(true).toBe(true);
  }, 60_000);

  test("on the move too: whatever is bought, whatever the menu, cooks' bodies touch for under 8% of the loop", () => {
    const cohorts: Cohort[] = ["p50", "p95", "p05", "p50"];
    // A body seen from above: 0.46 m across the shoulders, 0.30 m front to back.
    const reach = (yaw: number, ux: number, uz: number) => {
      const c = ux * Math.cos(yaw) - uz * Math.sin(yaw);
      const s = ux * Math.sin(yaw) + uz * Math.cos(yaw);
      return (0.23 * 0.15) / Math.hypot(0.15 * c, 0.23 * s);
    };
    const ms = BUYS.flatMap((b) => [0, 6, MENUS - 1].map((menu) => planMorning(morningOf(1400, menu, b))));
    for (const cooks of [2, 3, 4]) for (const menu of [0, MENUS - 1]) ms.push(planMorning({ meals: 700, menu, cooks, ahead: false }));
    for (const m of ms) {
      const crew = m.settings.cooks;
      const endMin = m.depart + 14;
      const rigs: CookRig[] = cohorts.slice(0, crew).map((c) => ({ st: planningBody(c), routine: null, walks: new Map(), tools: new Map(), held: null }));
      const loopS = planRoutines(rigs, planActs(m, crew, endMin), endMin * SEC_PER_MIN + END_HOLD_S);
      let touching = 0;
      for (let t = 0; t < loopS; t += 0.05) {
        const p = rigs.map((k) => k.routine!.rootAt(t));
        for (let i = 0; i < crew; i++)
          for (let j = i + 1; j < crew; j++) {
            const d = Math.hypot(p[j]!.x - p[i]!.x, p[j]!.z - p[i]!.z);
            const ux = (p[j]!.x - p[i]!.x) / Math.max(1e-6, d);
            const uz = (p[j]!.z - p[i]!.z) / Math.max(1e-6, d);
            if (d < reach(p[i]!.yaw, ux, uz) + reach(p[j]!.yaw, -ux, -uz) - 0.02) touching += 0.05;
          }
      }
      expect(touching / loopS).toBeLessThan(0.08);
    }
  });

  test("each cook's acts run end to end over the morning, in order, each stretch the work of one of the plan's cooks", () => {
    for (const b of BUYS)
      for (const menu of [0, 6, MENUS - 1]) {
        const m = planMorning(morningOf(1400, menu, b));
        const end = m.depart + 6;
        const acts = planActs(m, m.settings.cooks, end);
        for (const xs of acts) {
          expect(xs[0]!.t0).toBe(0);
          expect(xs[xs.length - 1]!.t1).toBeCloseTo(end * SEC_PER_MIN, 6);
          for (let k = 1; k < xs.length; k++) {
            expect(xs[k]!.t0).toBeCloseTo(xs[k - 1]!.t1, 6);
            expect(xs[k]!.t1).toBeGreaterThan(xs[k]!.t0);
          }
        }
        // Handing work over (swapTails) keeps every one of the plan's cooks on screen exactly once at every moment, so
        // the pause badge over a cook's head is the right cook's.
        for (let t = 0.05; t < end * SEC_PER_MIN; t += 0.5) {
          const who = acts.map((xs) => xs.find((x) => x.t0 <= t && x.t1 > t)!.who).sort();
          expect(who).toEqual(acts.map((_, c) => c));
        }
      }
  });
});
