import { describe, expect, test } from "bun:test";
import { BURNERS, DOOR, MENUS, OVENS, TABLE_SPOTS, madeIt, menuDishes, morningLimit, planMorning, portionsFor, type Crew, type Morning, type MorningSettings } from "./cateringMorning";

const ALL: MorningSettings[] = [];
for (const meals of [150, 700, 1100]) for (const cooks of [2, 3, 4] as Crew[]) for (let menu = 0; menu < MENUS; menu++) for (const ahead of [false, true]) ALL.push({ meals, cooks, menu, ahead });

/** The most of anything at once: intervals with a weight, swept by their ends (an end before a start at a tie, times
 * within a microsecond being a tie: the scheduler finishes what ends within that of the moment it is at). */
function peak(xs: { t0: number; t1: number; w: number }[]): number {
  const us = (t: number) => Math.round(t * 1e6);
  const ev = xs.flatMap((x) => [
    [us(x.t0), x.w],
    [us(x.t1), -x.w],
  ]);
  ev.sort((a, b) => a[0]! - b[0]! || a[1]! - b[1]!);
  let cur = 0;
  let max = 0;
  for (const [, w] of ev) max = Math.max(max, (cur += w!));
  return max;
}

const runs = new Map<MorningSettings, Morning>();
const run = (s: MorningSettings) => {
  let m = runs.get(s);
  if (!m) runs.set(s, (m = planMorning(s)));
  return m;
};

describe("cateringMorning", () => {
  test("every meal is cooked and packed once, each place's dish in its share", () => {
    for (const s of ALL) {
      const m = run(s);
      expect(m.finished).toBe(true);
      const want = portionsFor(s.meals);
      expect(want.slice(0, 3).reduce((a, b) => a + b, 0)).toBe(s.meals);
      const got = [0, 0, 0, 0];
      for (const b of m.batches) {
        got[b.slot]! += b.portions;
        expect(b.packed).toBeGreaterThan(0);
        expect(b.packed).toBeLessThanOrEqual(m.depart + 1e-6);
      }
      expect(got).toEqual(want);
      const dishes = menuDishes(s.menu);
      expect(m.pans).toBe(want.reduce((a, p, k) => a + Math.ceil(p / dishes[k]!.pan - 1e-6), 0));
      expect(m.batches.reduce((a, b) => a + b.pans, 0)).toBe(m.pans);
    }
  });

  test("no station or cook is ever asked for more than it has", () => {
    for (const s of ALL) {
      const m = run(s);
      expect(peak(m.burns.map((b) => ({ t0: b.t0, t1: b.t1, w: 1 })))).toBeLessThanOrEqual(BURNERS);
      for (let k = 0; k < BURNERS; k++) expect(peak(m.burns.filter((b) => b.burner === k).map((b) => ({ t0: b.t0, t1: b.t1, w: 1 })))).toBeLessThanOrEqual(1);
      for (let o = 0; o < OVENS; o++) expect(peak(m.loads.filter((l) => l.oven === o).map((l) => ({ t0: l.t0, t1: l.t3, w: 1 })))).toBeLessThanOrEqual(1);
      for (const l of m.loads) {
        expect(l.pans).toBeLessThanOrEqual(10);
        expect(l.t0 <= l.t1 && l.t1 <= l.t2 && l.t2 <= l.tu && l.tu <= l.t3).toBe(true);
        // Unloading starts when a cook does: the load's unload job runs from tu to t3.
        expect(m.jobs.some((j) => j.kind === "unload" && j.spot === l.oven && Math.abs(j.t0 - l.tu) < 1e-6 && Math.abs(j.t1 - l.t3) < 1e-6)).toBe(true);
      }
      const table = m.jobs.filter((j) => j.at === "table");
      expect(peak(table.map((j) => ({ t0: j.t0, t1: j.t1, w: 1 })))).toBeLessThanOrEqual(TABLE_SPOTS);
      for (let k = 0; k < TABLE_SPOTS; k++) expect(peak(table.filter((j) => j.spot === k).map((j) => ({ t0: j.t0, t1: j.t1, w: 1 })))).toBeLessThanOrEqual(1);
      for (let c = 0; c < s.cooks; c++) expect(peak(m.jobs.filter((j) => j.cook === c).map((j) => ({ t0: j.t0, t1: j.t1, w: j.hands })))).toBeLessThanOrEqual(1 + 1e-9);
      expect(m.jobs.every((j) => j.cook >= 0 && j.cook < s.cooks)).toBe(true);
    }
  });

  test("each batch's stages run in their order, one after another", () => {
    for (const s of ALL.filter((x) => x.meals === 700)) {
      const m = run(s);
      const dishes = menuDishes(s.menu);
      const byBatch = new Map<number, { t0: number; t1: number; kind: string }[]>();
      const add = (b: number, x: { t0: number; t1: number; kind: string }) => byBatch.set(b, [...(byBatch.get(b) ?? []), x]);
      for (const j of m.jobs) if (j.at !== "oven") add(j.batch, j);
      for (const l of m.loads) for (const b of l.batches) add(b, { t0: l.t0, t1: l.t3, kind: "oven" });
      m.batches.forEach((b, i) => {
        const xs = (byBatch.get(i) ?? []).sort((a, c) => a.t0 - c.t0);
        const d = dishes[b.slot]!;
        const stages = s.ahead && d.reheat ? d.reheat : s.ahead ? d.stages.filter((st) => !st.keeps) : d.stages;
        expect(xs.length).toBe(stages.length);
        for (let k = 1; k < xs.length; k++) expect(xs[k]!.t0).toBeGreaterThanOrEqual(xs[k - 1]!.t1 - 1e-6);
        expect(xs[xs.length - 1]!.kind).toBe("pack");
      });
    }
  });

  test("a plan that makes the van has nothing past its best", () => {
    for (const s of ALL) {
      const m = run(s);
      if (m.late <= 0.5) expect(m.pastBest).toBe(0);
      expect(m.depart).toBeGreaterThanOrEqual(DOOR);
      // A late plan's cost of making the van is never nothing: a way to make it with nothing past its best is the plan.
      if (m.onTimePastBest !== null && m.late > 0.5) {
        expect(m.onTimePastBest).toBeGreaterThan(0);
        expect(m.onTimePastBest).toBeLessThanOrEqual(s.meals);
      }
    }
  });

  test("a morning that only makes the van started at 6:00 still says what that costs", () => {
    // Two cooks, the written menu: moving starts earlier by the lateness never gets there, starting at 6:00 does.
    const m = planMorning({ meals: 500, cooks: 2, menu: 0, ahead: false });
    expect(m.late).toBeGreaterThan(30);
    expect(m.onTimePastBest).not.toBeNull();
    expect(m.onTimePastBest!).toBeGreaterThan(0);
  });

  test("the same settings give the same morning", () => {
    const s: MorningSettings = { meals: 700, cooks: 3, menu: 0, ahead: false };
    expect(JSON.stringify(planMorning({ ...s }))).toBe(JSON.stringify(planMorning({ ...s })));
  });

  test("at 700 meals with three cooks, the menu as written misses the van and the oven menu makes it", () => {
    const written = run(ALL.find((x) => x.meals === 700 && x.cooks === 3 && x.menu === 0 && !x.ahead)!);
    const swapped = run(ALL.find((x) => x.meals === 700 && x.cooks === 3 && x.menu === MENUS - 1 && !x.ahead)!);
    expect(written.late).toBeGreaterThan(30);
    expect(swapped.late).toBe(0);
    // The written menu keeps its cooks at the range: more hands-on time than the oven menu.
    expect(written.handsMin).toBeGreaterThan(1.5 * swapped.handsMin);
  });

  test("cooking the day before only moves work off the morning", () => {
    for (const s of ALL.filter((x) => x.ahead)) {
      const today = run(ALL.find((x) => x.meals === s.meals && x.cooks === s.cooks && x.menu === s.menu && !x.ahead)!);
      const m = run(s);
      expect(m.aheadMin).toBeGreaterThan(0);
      expect(m.handsMin).toBeLessThan(today.handsMin);
    }
  });

  test("the limit is the most meals in tens that make the van", () => {
    for (const menu of [0, MENUS - 1]) {
      const s = { cooks: 3 as Crew, menu, ahead: false };
      const lim = morningLimit(s);
      expect(lim).toBeGreaterThan(0);
      expect(madeIt(planMorning({ ...s, meals: lim }))).toBe(true);
      expect(madeIt(planMorning({ ...s, meals: lim + 10 }))).toBe(false);
    }
    expect(morningLimit({ cooks: 3, menu: MENUS - 1, ahead: false })).toBeGreaterThan(morningLimit({ cooks: 3, menu: 0, ahead: false }));
  });
});
