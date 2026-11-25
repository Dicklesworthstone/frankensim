import { describe, expect, test } from "bun:test";
import {
  DAYS,
  KITCHEN,
  MENU,
  SWAPS,
  WEEK_PLANS,
  menuFor,
  paidHours,
  reshapeWeek,
  rng,
  runCommissary,
  simulateWeek,
  weekPlan,
  weekPlanIndex,
  withCooks,
  type Kitchen,
  type Strategy,
} from "./commissaryWeek";

const ALL_SWAPS = Object.keys(SWAPS);
const share = (swaps: readonly string[]) => menuFor(swaps).reduce((s, d) => s + d.share, 0);

/** Weeks of random demand around the forecast, on a fixed seed. */
function weeks(count: number, sigma: number, seed = 7): number[][] {
  const r = rng(seed);
  return Array.from({ length: count }, () => KITCHEN.forecast.map((f) => Math.round(f * Math.exp(sigma * r.z()))));
}

const CASES: { strategy: Strategy; swaps: string[]; racks: number }[] = [];
for (const strategy of ["day", "ahead"] as const) for (const swaps of [[], ALL_SWAPS, ["chicken"]]) for (const racks of [1, 3, 5]) CASES.push({ strategy, swaps, racks });

describe("commissary week", () => {
  test("every portion of every meal is served, from stock or cooked on the day", () => {
    for (const c of CASES) {
      const k: Kitchen = { ...KITCHEN, racks: c.racks };
      for (const actual of weeks(40, 0.3)) {
        const wk = simulateWeek(k, menuFor(c.swaps), c.strategy, actual, 0.2);
        wk.days.forEach((d, i) => expect(d.fromStock + d.cookedToday).toBeCloseTo(actual[i]! * share(c.swaps), 6));
      }
    }
  });

  test("what is cooked ahead is served, thrown out, or carried into next week: nothing appears or vanishes", () => {
    for (const c of CASES) {
      const k: Kitchen = { ...KITCHEN, racks: c.racks };
      for (const actual of weeks(40, 0.3)) {
        const wk = simulateWeek(k, menuFor(c.swaps), c.strategy, actual, 0.2);
        const made = wk.days.reduce((s, d) => s + d.batches.reduce((t, b) => t + b.portions, 0), 0);
        const served = wk.days.reduce((s, d) => s + d.fromStock, 0);
        expect(served + wk.wasted + wk.carried).toBeCloseTo(made, 6);
      }
    }
  });

  test("the walk-in's racks and the blast chiller are never over-filled; batches are whole pans", () => {
    const chillMax = Math.floor(KITCHEN.shiftMin / KITCHEN.chillCycleMin) * KITCHEN.chillPans;
    for (const c of CASES) {
      const k: Kitchen = { ...KITCHEN, racks: c.racks };
      for (const actual of weeks(40, 0.3)) {
        const wk = simulateWeek(k, menuFor(c.swaps), c.strategy, actual, 0.2);
        for (const d of wk.days) {
          expect(d.pansEnd).toBeLessThanOrEqual(c.racks * KITCHEN.rackPans);
          expect(d.pansChilled).toBeLessThanOrEqual(chillMax);
          for (const b of d.batches) {
            const dish = menuFor(c.swaps).find((m) => m.key === b.dish)!;
            expect(b.portions).toBeCloseTo(b.pans * dish.perPan, 9);
          }
        }
      }
    }
  });

  test("a batch is cooked only for a day it keeps to, and only by a dish that keeps", () => {
    for (const c of CASES) {
      const menu = menuFor(c.swaps);
      for (const actual of weeks(20, 0.3)) {
        const wk = simulateWeek({ ...KITCHEN, racks: c.racks }, menu, c.strategy, actual, 0.2);
        wk.days.forEach((d, i) => {
          for (const b of d.batches) {
            const dish = menu.find((m) => m.key === b.dish)!;
            expect(b.forDay).toBeGreaterThan(i);
            expect(b.forDay - i).toBeLessThanOrEqual(dish.keeps);
          }
        });
      }
    }
  });

  test("the crew's time adds up: paid = shift + overtime = worked + idle, and nobody is both idle and on overtime", () => {
    const cap = KITCHEN.cooks * KITCHEN.shiftMin;
    for (const c of CASES) {
      for (const actual of weeks(20, 0.3)) {
        const wk = simulateWeek({ ...KITCHEN, racks: c.racks }, menuFor(c.swaps), c.strategy, actual, 0.2);
        for (const d of wk.days) {
          expect(d.paidMin).toBeCloseTo(cap + d.otMin, 6);
          expect(d.workedMin + d.idleMin).toBeCloseTo(d.paidMin, 6);
          expect(d.otMin > 0 && d.idleMin > 1e-9).toBe(false);
          expect(d.lateMin).toBeLessThanOrEqual(d.otMin + 1e-9);
        }
      }
    }
  });

  test("cooking on the day holds no stock and throws nothing out", () => {
    const run = runCommissary({ strategy: "day", swaps: ALL_SWAPS, sigma: 0.3 });
    for (const d of run.plan.days) expect(d.pansEnd + d.aheadMin + d.wasted).toBe(0);
    expect(run.wasted.p90).toBe(0);
  });

  test("the forecast week: cooking ahead moves Wednesday's overflow into Tuesday's idle time", () => {
    const day = runCommissary({ strategy: "day", swaps: [], sigma: 0 }).plan;
    const ahead = runCommissary({ strategy: "ahead", swaps: [], sigma: 0 }).plan;
    const tue = DAYS.indexOf("Tue");
    const wed = DAYS.indexOf("Wed");
    expect(day.days[wed]!.otMin).toBeGreaterThan(0);
    expect(day.days[tue]!.idleMin).toBeGreaterThan(0);
    expect(ahead.days[tue]!.aheadMin).toBeGreaterThan(0);
    expect(ahead.days[wed]!.otMin).toBeLessThan(day.days[wed]!.otMin);
    expect(ahead.days[tue]!.idleMin).toBeLessThan(day.days[tue]!.idleMin);
    expect(ahead.otMin).toBeLessThan(day.otMin);
    expect(ahead.utilization).toBeGreaterThan(day.utilization);
  });

  test("menu engineering: dishes that keep let the week level further than the menu as written", () => {
    const asWritten = runCommissary({ strategy: "ahead", swaps: [], sigma: 0 }).plan;
    const swapped = runCommissary({ strategy: "ahead", swaps: ALL_SWAPS, sigma: 0 }).plan;
    expect(asWritten.otMin).toBeGreaterThan(0);
    expect(swapped.otMin).toBeLessThan(asWritten.otMin);
    expect(swapped.lateMeals).toBe(0);
  });

  test("a bigger walk-in never makes the forecast week's overtime worse", () => {
    for (const swaps of [[], ALL_SWAPS]) {
      let prev = Infinity;
      for (let racks = 1; racks <= 6; racks++) {
        const ot = runCommissary({ strategy: "ahead", swaps, sigma: 0, racks }).plan.otMin;
        expect(ot).toBeLessThanOrEqual(prev + 1e-9);
        prev = ot;
      }
    }
  });

  test("the same settings give the same answer, and the swaps only change the dishes they name", () => {
    const a = runCommissary({ strategy: "ahead", swaps: ["salad"], sigma: 0.25 });
    const b = runCommissary({ strategy: "ahead", swaps: ["salad"], sigma: 0.25 });
    expect(JSON.stringify(a)).toBe(JSON.stringify(b));
    const menu = menuFor(["salad"]);
    menu.forEach((d, i) => expect(d === MENU[i] || d.key === "salad").toBe(true));
    for (const [key, alt] of Object.entries(SWAPS)) {
      const base = MENU.find((m) => m.key === key)!;
      expect(alt.share).toBe(base.share);
      expect(alt.keeps).toBeGreaterThan(base.keeps);
    }
  });
});

test("a fourth cook on the day and three cooking ahead both get the week's meals out, the fourth cook mostly idle", () => {
  const plan = (strategy: Strategy, cooks: number, swaps: readonly string[] = [], racks = 3) => runCommissary({ strategy, swaps, sigma: 0.1, racks }, { ...KITCHEN, cooks });
  const threeDay = plan("day", 3);
  const fourDay = plan("day", 4);
  const threeAhead = plan("ahead", 3);
  expect(threeDay.late.p50).toBeGreaterThan(20);
  expect(fourDay.late.p50).toBeLessThan(1);
  expect(threeAhead.late.p50).toBeLessThan(1);
  expect(fourDay.idle.p50).toBeGreaterThan(threeAhead.idle.p50 + 30);
  // Two cooks cannot carry the week, whatever is cooked ahead or swapped.
  expect(plan("ahead", 2, ALL_SWAPS, 5).late.p50).toBeGreaterThan(100);
});

test("the plan map's plans: every crew, timing, rack count and set of swaps once, numbered both ways", () => {
  const seen = new Set<string>();
  for (let i = 0; i < WEEK_PLANS; i++) {
    const p = weekPlan(i);
    expect(weekPlanIndex(p)).toBe(i);
    seen.add(`${p.cooks}|${p.strategy}|${p.racks}|${[...p.swaps].sort().join("+")}`);
  }
  expect(WEEK_PLANS).toBe(3 * 2 * 5 * 8);
  expect(seen.size).toBe(WEEK_PLANS);
  // The order the swaps are listed in does not change the plan.
  expect(weekPlanIndex({ cooks: 3, strategy: "ahead", racks: 2, swaps: ["salad", "chicken"] })).toBe(weekPlanIndex({ cooks: 3, strategy: "ahead", racks: 2, swaps: ["chicken", "salad"] }));
});

test("reshapeWeek keeps the week's meals and moves each day toward the average or away from it", () => {
  const total = (k: Kitchen) => k.forecast.reduce((a, b) => a + b, 0);
  expect(reshapeWeek(1)).toBe(KITCHEN);
  for (const swing of [0, 0.25, 0.5, 0.75, 1.25, 1.5]) {
    const k = reshapeWeek(swing);
    expect(total(k)).toBe(total(KITCHEN));
    for (const n of k.forecast) expect(Number.isInteger(n) && n > 0).toBe(true);
  }
  const flat = reshapeWeek(0).forecast;
  expect(Math.max(...flat) - Math.min(...flat)).toBeLessThanOrEqual(1);
  // The order of the days stays: Wednesday the busiest, Tuesday the quietest.
  const sharp = reshapeWeek(1.5).forecast;
  expect(sharp.indexOf(Math.max(...sharp))).toBe(2);
  expect(sharp.indexOf(Math.min(...sharp))).toBe(1);
  expect(sharp[1]!).toBeLessThan(KITCHEN.forecast[1]!);
});

test("the same meals, spread flat or sharp: cooking ahead buys nothing on a flat week, not enough on a sharp one", () => {
  const plan = (week: Kitchen, strategy: Strategy, cooks: number, swaps: readonly string[] = [], racks = 3) =>
    runCommissary({ strategy, swaps, sigma: 0.1, racks }, withCooks(cooks, week));
  const flat = reshapeWeek(0);
  const sharp = reshapeWeek(1.5);
  // Flat: three cooks on the day get every meal out with no overtime.
  const flatDay = plan(flat, "day", 3);
  expect(flatDay.late.p90).toBeLessThan(0.5);
  expect(flatDay.ot.p50).toBeLessThan(0.5);
  // Sharp: three cooks fall behind even cooking ahead; a fourth cook cooking ahead does not.
  expect(plan(sharp, "ahead", 3, ["chicken"], 2).late.p50).toBeGreaterThan(20);
  expect(plan(sharp, "ahead", 4, ["chicken"], 2).late.p90).toBeLessThan(0.5);
});

test("paid crew-hours: every cook's shift on every day, plus a median week's overtime", () => {
  for (const cooks of [2, 3, 4]) {
    const r = runCommissary({ strategy: "day", swaps: [], sigma: 0.1 }, withCooks(cooks));
    expect(paidHours(r, cooks)).toBeCloseTo(cooks * 8 * DAYS.length + r.ot.p50, 9);
  }
  expect(withCooks(KITCHEN.cooks)).toBe(KITCHEN);
});
