import { describe, expect, test } from "bun:test";
import {
  BERICK_C,
  BERICK_G_PER_H,
  BOIL_C,
  BURNER_W,
  GAS_INTO_POT,
  MAX_LITRES,
  POT_20QT,
  POT_40QT,
  SETUPS,
  TAP_C,
  boilAll,
  evaporationGPerHM2,
  heatUp,
  latentJPerKg,
  lossesW,
} from "./potBoil";

const byKey = (litres: number, flame: number) => Object.fromEntries(boilAll(litres, flame).map((r) => [r.setup.key, r]));

describe("potBoil: the measured inputs", () => {
  test("the evaporation rate passes through Berick's points and rises with temperature", () => {
    BERICK_C.forEach((c, i) => expect(evaporationGPerHM2(c) * 0.0491).toBeCloseTo(BERICK_G_PER_H[i]!, 6));
    let last = 0;
    for (let c = 20; c <= 100; c += 0.5) {
      const g = evaporationGPerHM2(c);
      expect(g).toBeGreaterThan(last);
      last = g;
    }
  });

  test("latent heat is the steam tables' to within a third of a percent", () => {
    const table: [number, number][] = [[0, 2.501e6], [20, 2.454e6], [40, 2.406e6], [60, 2.358e6], [80, 2.309e6], [100, 2.257e6]];
    for (const [c, h] of table) expect(Math.abs(latentJPerKg(c) / h - 1)).toBeLessThan(1 / 300);
  });

  test("Berick's pot at 98 °C: the lid takes nearly all the steam, and steam is most of an open pot's loss", () => {
    const wet = 5 / 1000 / ((Math.PI * 0.25 ** 2) / 4);
    const open = lossesW(98, 0.25, wet, false, 10);
    const lid = lossesW(98, 0.25, wet, true, 10);
    // 1128 g/h at 98 °C is 707 W of steam.
    expect(open.steam).toBeGreaterThan(690);
    expect(open.steam).toBeLessThan(720);
    expect(lid.steam / open.steam).toBeCloseTo(0.08, 6);
    expect(open.steam).toBeGreaterThan(2 * (open.radiation + open.air));
    expect(lid.radiation).toBeLessThan(open.radiation);
  });

  test("the pots are the drawn 20- and 40-quart stockpots", () => {
    expect(POT_20QT.litres).toBeGreaterThan(17);
    expect(POT_20QT.litres).toBeLessThan(19.5);
    expect(POT_40QT.litres).toBeGreaterThan(34);
    expect(POT_40QT.litres).toBeLessThan(39);
    expect(MAX_LITRES).toBeLessThanOrEqual(2 * 0.9 * POT_20QT.litres);
  });
});

describe("potBoil: heating to the boil", () => {
  test("Berick's Figure 7: fitted to his oil-film run, the open pot boils when he measured it (27 min)", () => {
    const pot = { d: 0.25, litres: 5, shellJK: 600, waterC: 15, roomC: 10, boilC: 99.6 };
    let lo = 300;
    let hi = 4000;
    let P = 0;
    for (let i = 0; i < 40; i++) {
      P = (lo + hi) / 2;
      const r = heatUp({ ...pot, lid: false, film: true, powerW: P });
      if ((r.boilMin ?? Infinity) > 23) lo = P;
      else hi = P;
    }
    const open = heatUp({ ...pot, lid: false, powerW: P });
    expect(open.boilMin).not.toBeNull();
    expect(open.boilMin!).toBeGreaterThan(26);
    expect(open.boilMin!).toBeLessThan(28.5);
  });

  test("with power to spare the time is the heat needed over the power", () => {
    const P = 2e6;
    const r = heatUp({ d: 0.3, litres: 10, shellJK: 3000, lid: true, powerW: P });
    const want = ((10 * 0.999 * 4186 + 3000) * (BOIL_C - TAP_C)) / P / 60;
    expect(r.boilMin!).toBeGreaterThan(want * 0.995);
    expect(r.boilMin!).toBeLessThan(want * 1.01 + 1 / 60);
  });

  test("every joule is accounted for: in = warming + steam + radiation + air", () => {
    for (const r of boilAll(24, 0.75)) {
      const kg0 = r.litres[0]! * 0.999;
      const kg1 = r.litres[r.litres.length - 1]! * 0.999;
      const shell = r.setup.pot.shellJK + (r.setup.lid ? r.setup.pot.lidJK : 0);
      const warm = ((kg0 + kg1) / 2) * 4186 * (BOIL_C - TAP_C) + shell * (BOIL_C - TAP_C);
      const out = warm + r.heat.steam + r.heat.radiation + r.heat.air;
      expect(Math.abs(out / r.heat.into - 1)).toBeLessThan(0.01);
    }
  });

  test("the samples are a quarter minute apart, rising, and end at the boil", () => {
    for (const r of boilAll(20, 1)) {
      expect(r.t[0]).toBe(0);
      expect(r.c[0]).toBe(TAP_C);
      for (let i = 1; i < r.t.length - 1; i++) {
        expect(r.t[i]! - r.t[i - 1]!).toBeCloseTo(0.25, 9);
        expect(r.c[i]!).toBeGreaterThan(r.c[i - 1]!);
      }
      expect(r.t[r.t.length - 1]).toBe(r.boilMin!);
      expect(r.c[r.c.length - 1]).toBe(BOIL_C);
    }
  });
});

describe("potBoil: the four ways", () => {
  test("a lid always helps, and two pots on two burners beat one on one", () => {
    for (const flame of [1, 0.75])
      for (const litres of [10, 20, MAX_LITRES]) {
        const k = byKey(litres, flame);
        expect(k["one-lid"]!.boilMin!).toBeLessThan(k["one-open"]!.boilMin!);
        expect(k["two-lid"]!.boilMin!).toBeLessThan(k["two-open"]!.boilMin!);
        expect(k["two-open"]!.boilMin!).toBeLessThan(k["one-open"]!.boilMin!);
        expect(k["two-lid"]!.boilMin!).toBeLessThan(k["one-lid"]!.boilMin!);
        // Twice the heat into the same water: about half the time, for about the same gas.
        const x = k["one-lid"]!.boilMin! / k["two-lid"]!.boilMin!;
        expect(x).toBeGreaterThan(1.8);
        expect(x).toBeLessThan(2.2);
        expect(Math.abs(k["two-lid"]!.burnerMin / k["one-lid"]!.burnerMin - 1)).toBeLessThan(0.1);
        // An open pot loses water on the way; a lid keeps nearly all of it.
        expect(k["one-open"]!.steamL).toBeGreaterThan(5 * k["one-lid"]!.steamL);
      }
  });

  test("more water, or less flame, takes longer", () => {
    for (const s of SETUPS) {
      const t = (litres: number, flame: number) => boilAll(litres, flame).find((r) => r.setup.key === s.key)!.boilMin!;
      expect(t(20, 1)).toBeGreaterThan(t(12, 1));
      expect(t(20, 0.75)).toBeGreaterThan(t(20, 1));
    }
  });

  test("on half a flame the big open pot does not boil in 3 h: it levels off short of 100 °C, steam carrying most of the heat", () => {
    const k = byKey(30, 0.5);
    const oo = k["one-open"]!;
    expect(oo.boilMin).toBeNull();
    const peak = oo.c[oo.c.length - 1]!;
    expect(peak).toBeGreaterThan(97);
    expect(peak).toBeLessThan(BOIL_C);
    // It has all but settled: the last ten minutes add under a quarter of a degree, and at that temperature it
    // loses all but a hair of what the burner gives it.
    expect(peak - oo.c[oo.c.length - 1 - 40]!).toBeLessThan(0.25);
    const l = lossesW(peak, POT_40QT.d, oo.litres[oo.litres.length - 1]! / 1000 / ((Math.PI * POT_40QT.d ** 2) / 4), false);
    const share = (l.steam + l.radiation + l.air) / (BURNER_W * GAS_INTO_POT * 0.5);
    expect(share).toBeGreaterThan(0.97);
    expect(share).toBeLessThan(1);
    // At every amount of water the page offers: levelled off, short of the boil, most of the heat going as steam.
    for (let litres = 8; litres <= MAX_LITRES; litres++) {
      const r = boilAll(litres, 0.5).find((x) => x.setup.key === "one-open")!;
      const n = r.c.length;
      expect(r.boilMin).toBeNull();
      expect(r.c[n - 1]!).toBeLessThan(BOIL_C);
      expect(r.c[n - 1]! - r.c[n - 41]!).toBeLessThan(0.3);
      const at = lossesW(r.c[n - 1]!, POT_40QT.d, r.litres[n - 1]! / 1000 / ((Math.PI * POT_40QT.d ** 2) / 4), false);
      expect(at.steam / (BURNER_W * GAS_INTO_POT * 0.5)).toBeGreaterThan(0.75);
    }
    // With lids, or split in two, the same water does boil.
    expect(k["one-lid"]!.boilMin).not.toBeNull();
    expect(k["two-lid"]!.boilMin).not.toBeNull();
  });
});
