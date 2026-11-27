import { describe, expect, test } from "bun:test";
import { cToF, fToC } from "./coolDown";
import { HOLD_WAYS, PRIME_MIN, ROOM_C, SAMPLE_MIN, SLICE_MIN, WALL_NODES, holdAll, holdRun, sampleOf, samplesTo, type HoldKey, type HoldRun } from "./holdHot";

const way = (k: HoldKey) => HOLD_WAYS.find((w) => w.key === k)!;
const hot = new Map<number, HoldRun[]>();
const runsAt = (f: number) => {
  let r = hot.get(f);
  if (!r) hot.set(f, (r = holdAll(fToC(f), { to141: true })));
  return r;
};
const of = (rs: HoldRun[], k: HoldKey) => rs.find((r) => r.way.key === k)!;

describe("chili kept hot in a carrier", () => {
  test("what the chili and the shell lose is what left through the shell", () => {
    for (const r of runsAt(185)) {
      const { at0, atEnd, lost } = r.energy;
      expect(at0 - atEnd).toBeGreaterThan(0);
      expect(Math.abs(at0 - atEnd - lost) / (at0 - atEnd)).toBeLessThan(1e-3);
    }
  });

  test("chili at the room's warmth in a carrier at the room's warmth stays as it is", () => {
    const r = holdRun(way("fullCold"), ROOM_C, { minMin: 60, maxMin: 60 });
    for (let k = samplesTo(0) - 1; k < samplesTo(60); k++) expect(Math.abs(r.coldest[k]! - ROOM_C)).toBeLessThan(1e-5);
  });

  test("a finer grid and a shorter step give the same hours", () => {
    for (const k of ["oneCold", "fullPrimed"] as const) {
      const a = holdRun(way(k), fToC(185), { to141: true });
      const b = holdRun(way(k), fToC(185), { to141: true, grid: [10, 14, 10], dt: 10 });
      expect(Math.abs(b.hotMin! - a.hotMin!) / b.hotMin!).toBeLessThan(0.03);
      expect(Math.abs(b.mid141! - a.mid141!) / b.mid141!).toBeLessThan(0.02);
    }
  }, 120000);

  test("full beats one pan, primed beats cold, hotter packed holds longer", () => {
    for (const f of [185, 165]) {
      const rs = runsAt(f);
      const h = (k: HoldKey) => of(rs, k).hotMin!;
      expect(h("fullPrimed")).toBeGreaterThan(h("fullCold"));
      expect(h("fullCold")).toBeGreaterThan(h("onePrimed"));
      expect(h("onePrimed")).toBeGreaterThan(h("oneCold"));
    }
    for (const w of HOLD_WAYS) expect(of(runsAt(185), w.key).hotMin!).toBeGreaterThan(of(runsAt(165), w.key).hotMin!);
    // Priming matters far more for a lone pan, which has a whole carrier's walls to warm, than for six.
    const gain = (a: HoldKey, b: HoldKey) => of(runsAt(185), a).hotMin! / of(runsAt(185), b).hotMin!;
    expect(gain("onePrimed", "oneCold")).toBeGreaterThan(1.3);
    expect(gain("fullPrimed", "fullCold")).toBeLessThan(1.15);
  });

  test("against the carrier patent: a lone primed pan holds about six hours above 141 °F, six pans longer than fifteen", () => {
    const rs = runsAt(185);
    expect(of(rs, "onePrimed").mid141! / 60).toBeGreaterThan(5);
    expect(of(rs, "onePrimed").mid141! / 60).toBeLessThan(7.5);
    expect(of(rs, "fullPrimed").mid141! / 60).toBeGreaterThan(15);
  });

  test("priming warms the shell; a cold carrier stands at the room's warmth; the edge is always cooler than the middle", () => {
    const rs = runsAt(185);
    const k0 = samplesTo(0) - 1;
    for (const r of rs) {
      const inner = cToF(r.wall[k0 * WALL_NODES]!);
      const outer = cToF(r.wall[(k0 + 1) * WALL_NODES - 1]!);
      if (r.way.primed) {
        expect(inner).toBeGreaterThan(115);
        expect(inner).toBeLessThan(160);
        expect(outer).toBeLessThan(85);
        expect(r.water[k0 - 1]!).toBeLessThan(r.water[0]!);
      } else {
        expect(Math.abs(inner - 72)).toBeLessThan(1e-3);
        expect(Number.isNaN(r.water[0]!)).toBe(true);
      }
      for (let m = SAMPLE_MIN; m <= 6 * 60; m += SAMPLE_MIN) for (const c of r.centre) expect(sampleOf(r.coldest, m)).toBeLessThan(sampleOf(c, m));
    }
  });

  test("a page's run covers the stage's hours, the cut faces filled to the end of them", () => {
    const r = holdRun(way("oneCold"), fToC(165));
    expect(r.endMin).toBeGreaterThanOrEqual(SLICE_MIN);
    const [nx, ny] = r.grid;
    const per = r.way.pans * nx * ny;
    for (let k = samplesTo(0) - 1; k < samplesTo(SLICE_MIN); k++) expect(Number.isNaN(r.slice[k * per]!)).toBe(false);
    expect(Number.isNaN(r.slice[(samplesTo(0) - 2) * per]!)).toBe(true);
    expect(PRIME_MIN).toBe(30);
  });
});
