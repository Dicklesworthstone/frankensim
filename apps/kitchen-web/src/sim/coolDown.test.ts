import { describe, expect, test } from "bun:test";
import {
  CHART_MIN,
  CHILI,
  DAY_MIN,
  EMPTY_WALK_IN,
  H_EMPTY,
  H_FULL,
  MEASURED,
  PAGE_RUN,
  SAMPLE_MIN,
  SNAP_MIN,
  STAGE1_MIN,
  TOTAL_MIN,
  WALK_IN,
  coolAll,
  coolPot,
  coolSlab,
  coolWay,
  foodProps,
  hoursOf,
  restLitres,
  type Face,
} from "./coolDown";

/** Bessel J0 and J1: Abramowitz and Stegun 9.4.1-9.4.6 (to about 1e-8). */
function j0(x: number): number {
  const a = Math.abs(x);
  if (a <= 3) {
    const y = (x / 3) ** 2;
    return 1 + y * (-2.2499997 + y * (1.2656208 + y * (-0.3163866 + y * (0.0444479 + y * (-0.0039444 + y * 0.00021)))));
  }
  const y = 3 / a;
  const f = 0.79788456 + y * (-0.00000077 + y * (-0.0055274 + y * (-0.00009512 + y * (0.00137237 + y * (-0.00072805 + y * 0.00014476)))));
  const th = a - 0.78539816 + y * (-0.04166397 + y * (-0.00003954 + y * (0.00262573 + y * (-0.00054125 + y * (-0.00029333 + y * 0.00013558)))));
  return (f * Math.cos(th)) / Math.sqrt(a);
}
function j1(x: number): number {
  const a = Math.abs(x);
  if (a <= 3) {
    const y = (x / 3) ** 2;
    return x * (0.5 + y * (-0.56249985 + y * (0.21093573 + y * (-0.03954289 + y * (0.00443319 + y * (-0.00031761 + y * 0.00001109))))));
  }
  const y = 3 / a;
  const f = 0.79788456 + y * (0.00000156 + y * (0.01659667 + y * (0.00017105 + y * (-0.00249511 + y * (0.00113653 + y * -0.00020033)))));
  const th = a - 2.35619449 + y * (0.12499612 + y * (0.0000565 + y * (-0.00637879 + y * (0.00074348 + y * (0.00079824 + y * -0.00029166)))));
  return (Math.sign(x) * f * Math.cos(th)) / Math.sqrt(a);
}
/** The roots of f on (0, top), by scanning for sign changes and halving. */
function roots(f: (x: number) => number, top: number, count: number): number[] {
  const out: number[] = [];
  let a = 1e-9;
  let fa = f(a);
  for (let b = 0.01; b < top && out.length < count; b += 0.01) {
    const fb = f(b);
    if (fa * fb < 0) {
      let lo = a;
      let hi = b;
      for (let i = 0; i < 60; i++) {
        const m = (lo + hi) / 2;
        if (f(lo) * f(m) <= 0) hi = m;
        else lo = m;
      }
      out.push((lo + hi) / 2);
    }
    a = b;
    fa = fb;
  }
  return out;
}
/** The middle of a slab 2L thick, both faces cooled by h, as a share of its start over the air: the series solution. */
function slabMiddle(bi: number, fo: number): number {
  let s = 0;
  for (const l of roots((x) => x * Math.sin(x) - bi * Math.cos(x), 60, 12)) s += ((4 * Math.sin(l)) / (2 * l + Math.sin(2 * l))) * Math.exp(-l * l * fo);
  return s;
}
/** The axis of an infinite cylinder, its side cooled by h, likewise. */
function cylinderAxis(bi: number, fo: number): number {
  let s = 0;
  for (const z of roots((x) => x * j1(x) - bi * j0(x), 60, 12)) s += ((2 * j1(z)) / (z * (j0(z) ** 2 + j1(z) ** 2))) * Math.exp(-z * z * fo);
  return s;
}

const plain = (): Face => ({ kind: "wall", eps: 0 });

describe("coolDown", () => {
  test("the food's properties: water's own from its parts, and chili's in the range of a wet, beany stew", () => {
    const w = foodProps({ water: 1, protein: 0, fat: 0, carb: 0, fiber: 0, ash: 0 }, 20);
    expect(Math.abs(w.k / 0.598 - 1)).toBeLessThan(0.015);
    expect(Math.abs(w.cp / 4182 - 1)).toBeLessThan(0.005);
    expect(Math.abs(w.rho / 998.2 - 1)).toBeLessThan(0.003);
    const c = foodProps(CHILI, 20);
    expect(c.k).toBeGreaterThan(0.5);
    expect(c.k).toBeLessThan(0.56);
    expect(c.cp).toBeGreaterThan(3500);
    expect(c.cp).toBeLessThan(3700);
    expect(c.rho).toBeGreaterThan(1050);
    expect(c.rho).toBeLessThan(1080);
  });

  test("the Bessel functions the next test leans on", () => {
    expect(j0(1)).toBeCloseTo(0.7651977, 6);
    expect(j1(1)).toBeCloseTo(0.4400506, 6);
    expect(j0(5)).toBeCloseTo(-0.1775968, 6);
    expect(j1(5)).toBeCloseTo(-0.3275791, 6);
    expect(Math.abs(j0(2.4048256))).toBeLessThan(1e-6);
    expect(Math.abs(j1(3.8317060))).toBeLessThan(1e-6);
  });

  test("with fixed properties and plain convection, the solvers follow the exact solutions", () => {
    const box = { c: 0, rh: 0.85, h: 12 };
    const p = foodProps(CHILI, 20);
    const alpha = p.k / (p.rho * p.cp);
    // A slab 5 cm thick from 50 °C.
    const depth = 0.05;
    const slab = coolSlab(depth, plain(), plain(), { box }, { startC: 50, constAt: 20, minMin: 300, maxMin: 300 });
    const bis = (box.h * depth) / 2 / p.k;
    for (const min of [30, 60, 120, 240]) {
      const exact = 50 * slabMiddle(bis, (alpha * min * 60) / (depth / 2) ** 2);
      expect(Math.abs(slab.centre[min / SAMPLE_MIN]! - exact)).toBeLessThan(0.02);
    }
    // A pot of it 30 cm across and 20 cm deep: a slab across its depth times a cylinder across its width.
    const r = 0.15;
    const fill = 0.2;
    const pot = coolPot(r, fill, plain(), plain(), plain(), { box }, { startC: 50, constAt: 20, minMin: 600, maxMin: 600 });
    const biz = (box.h * fill) / 2 / p.k;
    const bir = (box.h * r) / p.k;
    for (const min of [60, 180, 360, 600]) {
      const t = min * 60;
      const exact = 50 * slabMiddle(biz, (alpha * t) / (fill / 2) ** 2) * cylinderAxis(bir, (alpha * t) / r ** 2);
      expect(Math.abs(pot.centre[min / SAMPLE_MIN]! - exact)).toBeLessThan(0.1);
    }
  });

  test("the heat out of every face adds up to what the food lost, open and covered, pan and pot", () => {
    for (const covered of [false, true])
      for (const key of ["pot", "shallow", "freezer"] as const) {
        const w = coolWay(key, 12, covered, { maxMin: 8 * 60, minMin: 8 * 60 });
        expect(Math.abs(w.run.heatOut / w.run.heatDrop - 1)).toBeLessThan(1e-3);
        if (covered) expect(w.run.steamKg).toBe(0);
        else expect(w.run.steamKg).toBeGreaterThan(0);
        // The steam's heat: its water at the latent heat somewhere between 100 °C and 0 °C.
        expect(w.run.steamJ).toBeGreaterThanOrEqual(2.25e6 * w.run.steamKg);
        expect(w.run.steamJ).toBeLessThanOrEqual(2.51e6 * w.run.steamKg);
        expect(w.run.steamJ).toBeLessThan(w.run.heatOut);
      }
  });

  test("against the measured runs: each time to 70 °F and to 41 °F within a fifth, each study's draw the best fit", () => {
    const errs: number[] = [];
    const score = (study: string, h: number) => {
      let s = 0;
      for (const m of MEASURED.filter((x) => x.study === study)) {
        const r = m.run(h);
        const a = (r.t70! - r.t135!) / 60 / hoursOf(m.to70);
        const b = (r.t41! - r.t135!) / 60 / hoursOf(m.to41);
        s += Math.log(a) ** 2 + Math.log(b) ** 2;
        if (h === H_FULL || h === H_EMPTY) errs.push(a - 1, b - 1);
      }
      return s;
    };
    for (const [study, h] of [
      ["Roberts 2013", H_FULL],
      ["Olds 2019", H_EMPTY],
    ] as const) {
      const best = score(study, h);
      expect(score(study, h - 0.5)).toBeGreaterThan(best);
      expect(score(study, h + 0.5)).toBeGreaterThan(best);
    }
    expect(errs.length).toBe(MEASURED.length * 2);
    for (const e of errs) expect(Math.abs(e)).toBeLessThan(0.2);
    expect(errs.reduce((s, e) => s + Math.abs(e), 0) / errs.length).toBeLessThan(0.12);
  });

  test("the measured verdicts: in Roberts's walk-ins only 2-inch pans in the freezer met both limits, as here", () => {
    for (const m of MEASURED.filter((x) => x.study === "Roberts 2013")) {
      const r = m.run(H_FULL);
      const ok = r.t70! - r.t135! <= STAGE1_MIN && r.t41! - r.t135! <= TOTAL_MIN;
      const measured = hoursOf(m.to70) <= 2 && hoursOf(m.to41) <= 6;
      expect(ok).toBe(measured);
    }
    // And in Olds's empty walk-in none met both.
    for (const m of MEASURED.filter((x) => x.study === "Olds 2019")) {
      const r = m.run(H_EMPTY);
      expect(r.t70! - r.t135! <= STAGE1_MIN && r.t41! - r.t135! <= TOTAL_MIN).toBe(false);
    }
    expect(EMPTY_WALK_IN.h).toBe(H_EMPTY);
    expect(WALK_IN.h).toBe(H_FULL);
  });

  test("the four ways: slower the deeper, slower covered, only the pot slower for a bigger batch", () => {
    const at = (l: number, c: boolean) => Object.fromEntries(coolAll(l, c).map((w) => [w.way.key, w]));
    for (const covered of [false, true]) {
      const w = at(30, covered);
      expect(w.pot!.total!).toBeGreaterThan(w.deep!.total!);
      expect(w.deep!.total!).toBeGreaterThan(w.shallow!.total!);
      expect(w.shallow!.total!).toBeGreaterThan(w.freezer!.total!);
      const small = at(12, covered);
      expect(small.pot!.total!).toBeLessThan(w.pot!.total!);
      for (const k of ["deep", "shallow", "freezer"]) expect(small[k]!.total).toBe(w[k]!.total);
      expect(small.shallow!.pans).toBeLessThan(w.shallow!.pans);
    }
    const open = at(30, false);
    const film = at(30, true);
    for (const k of ["pot", "deep", "shallow", "freezer"]) expect(film[k]!.total!).toBeGreaterThan(open[k]!.total!);
  });

  test("a part-filled last pan: shallower, faster, on the same rack; every page run reaches the chart's end", () => {
    expect(restLitres("shallow", 22)).toBeCloseTo(1, 9);
    expect(restLitres("shallow", 21)).toBe(0);
    expect(restLitres("deep", 30)).toBeCloseTo(30 - 2 * 10.2, 9);
    expect(restLitres("pot", 30)).toBe(0);
    for (const key of ["shallow", "freezer"] as const) {
      const w = coolWay(key, 22, false, PAGE_RUN);
      expect(w.pans).toBe(4);
      expect(w.rest).not.toBeNull();
      expect(w.rest!.shape.kind === "slab" && w.rest!.shape.depth).toBeCloseTo((2 * 0.0254) / 7, 9);
      expect(w.rest!.t41!).toBeLessThan(w.run.t41!);
      // It moves with its rack, when the full pans are cold.
      expect(w.rest!.movedMin).toBe(key === "freezer" ? w.run.movedMin : null);
      for (const r of [w.run, w.rest!]) expect(r.endMin).toBeGreaterThanOrEqual(r.t135! + CHART_MIN);
    }
    expect(coolWay("shallow", 21, false, PAGE_RUN).rest).toBeNull();
    // Less than a pan: the one pan, filled shallower, is the way's run.
    const one = coolWay("deep", 8, false, PAGE_RUN);
    expect(one.pans).toBe(1);
    expect(one.rest).toBeNull();
    expect(one.run.shape.kind === "slab" && one.run.shape.depth).toBeCloseTo((3 * 0.0254 * 8) / 10.2, 9);
    expect(one.total!).toBeLessThan(coolWay("deep", 12, false, PAGE_RUN).total!);
    expect(coolWay("pot", 22, false, PAGE_RUN).rest).toBeNull();
  });

  test("the page's batch: only the freezer gets it there in time; every run covers the day, sampled and snapped", () => {
    const ways = coolAll(30, false);
    for (const w of ways) {
      expect(w.ok1 && w.ok2).toBe(w.way.key === "freezer");
      expect(w.run.endMin).toBeGreaterThanOrEqual(DAY_MIN);
      expect(w.run.centre.length).toBe(w.run.endMin / SAMPLE_MIN + 1);
      expect(w.run.snaps.length).toBe(Math.floor(w.run.endMin / SNAP_MIN) + 1);
      expect(w.run.t41).not.toBeNull();
    }
    const fz = ways.find((w) => w.way.key === "freezer")!;
    expect(fz.run.movedMin).not.toBeNull();
    expect(fz.run.movedMin! - fz.run.t41!).toBeLessThan(SAMPLE_MIN);
    // In the walk-in after, it settles to the walk-in's air and does not freeze.
    expect(fz.run.centre[fz.run.centre.length - 1]!).toBeGreaterThan(WALK_IN.c - 0.5);
    expect(fz.run.centre[fz.run.centre.length - 1]!).toBeLessThan(WALK_IN.c + 0.5);
  });
});
