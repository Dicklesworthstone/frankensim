import { describe, expect, test } from "bun:test";
import { CC_CFM, CODE_CFM, LPS_PER_CFM, SETS, SETUPS, catches, caughtCount, needCfm, perFoot, warmBtuh } from "./hoodCapture";

describe("broilers under a hood that stays put", () => {
  test("the lab's rates: less air as the broilers stand further back, less again with the gap closed", () => {
    for (const key of ["open", "sealed"] as const) for (let i = 1; i < 3; i++) expect(CC_CFM[key][i]!).toBeLessThan(CC_CFM[key][i - 1]!);
    for (let i = 0; i < 3; i++) expect(CC_CFM.sealed[i]!).toBeLessThan(CC_CFM.open[i]!);
    // The per-foot figures the lab gave with them, and the handbook's (L/s per metre) for the two open ones nearest the edge.
    expect([perFoot(CC_CFM.open[0]), perFoot(CC_CFM.open[1]), perFoot(CC_CFM.open[2]), perFoot(CC_CFM.sealed[2])]).toEqual([510, 440, 400, 280]);
    expect((CC_CFM.open[0] * LPS_PER_CFM) / 3.048).toBeCloseTo(790, -1);
    expect((CC_CFM.open[1] * LPS_PER_CFM) / 3.048).toBeCloseTo(680, -1);
  });

  test("where the broilers stand: the hood's depth is their depth, the gap and the push-back at every setting", () => {
    const depths = SETS.map((s) => s.gapIn + s.backIn);
    for (const d of depths) expect(d).toBeCloseTo(depths[0]!, 6);
    for (const s of SETS) expect(s.cookIn - s.backIn).toBe(6);
  });

  test("at the code's 4,000 cfm: open, only the furthest back is caught; sealed, all three are", () => {
    expect(CODE_CFM).toBe(4000);
    const got = SETUPS.map((s) => `${s.at}${s.sealed ? "s" : "o"}:${catches(CODE_CFM, s)}`);
    expect(got).toEqual(["0o:false", "0s:true", "1o:false", "1s:true", "2o:true", "2s:true"]);
    expect(caughtCount(CODE_CFM)).toBe(4);
    expect(caughtCount(2799)).toBe(0);
    expect(caughtCount(5100)).toBe(6);
    // Closing the gap does more than a foot of push-back.
    expect(needCfm({ at: 0, sealed: true })).toBeLessThan(needCfm({ at: 2, sealed: false }));
  });

  test("warming the fan's air: the standard-air rule, and nothing to warm on a warm day", () => {
    expect(warmBtuh(4000, 40)).toBeCloseTo(1.08 * 4000 * 30, 6);
    expect(warmBtuh(4000, 70)).toBe(0);
    expect(warmBtuh(4000, 95)).toBe(0);
    // Linear in the air: half the fan, half the heat.
    expect(warmBtuh(2000, 40) * 2).toBeCloseTo(warmBtuh(4000, 40), 6);
  });
});
