import { describe, expect, test } from "bun:test";
import { DAYS, KITCHEN, MENU, SWAPS, menuFor, rng, runCommissary, simulateWeek, type Kitchen, type Strategy } from "./commissaryWeek";

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
