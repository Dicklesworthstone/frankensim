import { describe, expect, test } from "bun:test";
import { planningBody } from "../three/models/people";
import { COHORT_SPEED, CREW_COHORTS, DWELL, FROM_S, KITCHENS, MAX_CREW, MORNINGS, RUN_S, footprintsOf, obstacles, runCrowd, runMornings, stanceOf, type CrowdRun, type CrowdSet } from "./crewCrowd";

const cache = new Map<string, CrowdRun>();
const run = (k: string, n: number, ghosts = false) => {
  const key = `${k}/${n}/${ghosts}`;
  let r = cache.get(key);
  if (!r) cache.set(key, (r = runCrowd(KITCHENS[k]!, n, RUN_S, FROM_S, ghosts)));
  return r;
};
const sets = new Map<string, CrowdSet>();
const pooled = (k: string, n: number, ghosts = false) => {
  const key = `${k}/${n}/${ghosts}`;
  let r = sets.get(key);
  if (!r) sets.set(key, (r = runMornings(KITCHENS[k]!, n, ghosts)));
  return r;
};
/** One pan more or less in the counted span, an hour. */
const ONE_PAN = 3600 / (RUN_S - FROM_S);

describe("the kitchens", () => {
  test("every place a cook stands is on the floor and clear of every piece by a walking body's width", () => {
    for (const k of Object.values(KITCHENS)) {
      const rects = obstacles(k);
      for (const cohort of CREW_COHORTS) {
        const { rWalk } = footprintsOf(planningBody(cohort));
        for (const s of k.spots) {
          const p = stanceOf(k, s, rWalk);
          expect(p.x - rWalk).toBeGreaterThan(k.bounds.minX);
          expect(p.x + rWalk).toBeLessThan(k.bounds.maxX);
          expect(p.z - rWalk).toBeGreaterThan(k.bounds.minZ);
          expect(p.z + rWalk).toBeLessThan(k.bounds.maxZ);
          for (const r of rects) {
            const dx = Math.max(r.minX - p.x, 0, p.x - r.maxX);
            const dz = Math.max(r.minZ - p.z, 0, p.z - r.maxZ);
            expect(Math.hypot(dx, dz)).toBeGreaterThanOrEqual(rWalk);
          }
        }
      }
    }
  });
  test("the island pulled back is the same island 18 in further from the line, the rest unmoved", () => {
    const a = KITCHENS.built!.pieces;
    const b = KITCHENS.moved!.pieces;
    expect(b.length).toBe(a.length);
    a.forEach((p, i) => {
      const q = b[i]!;
      expect(q.id).toBe(p.id);
      expect(q.fx).toBeCloseTo(p.fx, 9);
      if (p.yaw === Math.PI || !p.id) expect(q.fz - p.fz).toBeCloseTo(0.457, 3);
      else expect(q.fz).toBeCloseTo(p.fz, 9);
    });
  });
});

describe("one cook", () => {
  test("nobody is ever in the way: bodies and ghosts turn out the same pans, and the cook never waits on anyone", () => {
    for (const k of Object.keys(KITCHENS)) {
      const b = run(k, 1);
      const g = run(k, 1, true);
      expect(b.pans).toBe(g.pans);
      expect(b.cooks[0]!.waitBody).toBe(0);
      expect(b.waits.length).toBe(0);
      expect(b.pans).toBeGreaterThan(0);
    }
  });
  test("a pan takes the loop's work and the walks between: no faster than the work alone", () => {
    const b = run("built", 1);
    const perPan = DWELL.cold + DWELL.prep + DWELL.hot + DWELL.pack;
    expect(b.perHour).toBeLessThanOrEqual(3600 / perPan + ONE_PAN);
    expect(b.perHour).toBeGreaterThan(3600 / (perPan + 60));
  });
});

describe("a crew", () => {
  test("is the same every time it runs", () => {
    const a = runCrowd(KITCHENS.built!, 6);
    const b = runCrowd(KITCHENS.built!, 6);
    expect(JSON.stringify(b)).toBe(JSON.stringify(a));
  });

  for (const k of Object.keys(KITCHENS))
    for (const n of [2, 4, 6, 8])
      test(`${k}, ${n} cooks: nobody walks into anybody, nobody gridlocks, and bodies never out-cook ghosts`, () => {
        const b = run(k, n);
        const g = run(k, n, true);
        expect(b.gridlock).toBe(false);
        expect(b.minGap).not.toBeNull();
        expect(b.minGap!).toBeGreaterThanOrEqual(0);
        expect(b.perHour).toBeLessThanOrEqual(g.perHour + ONE_PAN);
        expect(b.pans).toBeGreaterThan(0);
      });

  test("every cook's day holds together: legs in order, nobody leaves before their work is done, walks at their pace", () => {
    for (const k of Object.keys(KITCHENS))
      for (const ghosts of [false, true]) {
        const r = run(k, MAX_CREW, ghosts);
        let pans = 0;
        r.cooks.forEach((day, i) => {
          expect(day.cohort).toBe(CREW_COHORTS[i % CREW_COHORTS.length]!);
          let free = 0;
          for (const l of day.legs) {
            expect(l.ready).toBeGreaterThanOrEqual(free - 1e-6);
            expect(l.go).toBeGreaterThanOrEqual(l.ready - 1e-6);
            expect(l.arrive).toBeGreaterThanOrEqual(l.go);
            expect(l.speed).toBeLessThanOrEqual(COHORT_SPEED[day.cohort] + 1e-9);
            // Ghosts wait only for a free place, never on a body.
            if (ghosts) expect(l.blocked).toBe(false);
            free = l.arrive + l.dwell;
          }
          if (ghosts) {
            expect(day.waitBody).toBe(0);
            expect(day.aside).toBe(0);
          }
          pans += day.legs.filter((l) => l.kind === "pack" && l.arrive + l.dwell > r.from && l.arrive + l.dwell <= r.seconds + 1e-6).length;
        });
        expect(pans).toBe(r.pans);
        expect(r.perHour).toBeCloseTo((r.pans * 3600) / (r.seconds - r.from), 9);
      }
  });

  test("eight cooks, three mornings pooled, lose more to each other as built than with the island pulled back", () => {
    const lost = (k: string) => 1 - pooled(k, 8).perHour / pooled(k, 8, true).perHour;
    expect(lost("built")).toBeGreaterThan(0.08);
    expect(lost("moved")).toBeLessThan(lost("built") - 0.08);
  });

  test("the pooled mornings add up: pans over all the counted spans, the morning shown the middle one", () => {
    const s = pooled("built", 8);
    const runs = MORNINGS.map((st) => runCrowd(KITCHENS.built!, 8, RUN_S, FROM_S, false, undefined, st));
    const pans = runs.reduce((a, r) => a + r.pans, 0);
    expect(s.perHour).toBeCloseTo((pans * 3600) / (MORNINGS.length * (RUN_S - FROM_S)), 9);
    expect(s.waits.length).toBe(runs.reduce((a, r) => a + r.waits.length, 0));
    expect(JSON.stringify(s.shown)).toBe(JSON.stringify(run("built", 8)));
    expect(s.minGap).toBe(Math.min(...runs.map((r) => r.minGap!)));
  });

  test("waits on someone are counted where they happened, and only from the counted span on", () => {
    const r = run("built", 8);
    expect(r.waits.length).toBeGreaterThan(0);
    for (const [x, z, s] of r.waits) {
      expect(s).toBeGreaterThan(0);
      expect(x).toBeGreaterThan(KITCHENS.built!.bounds.minX);
      expect(x).toBeLessThan(KITCHENS.built!.bounds.maxX);
      expect(z).toBeGreaterThan(KITCHENS.built!.bounds.minZ);
      expect(z).toBeLessThan(KITCHENS.built!.bounds.maxZ);
    }
    const total = r.cooks.reduce((a, c) => a + c.waitBody, 0);
    expect(r.waits.reduce((a, w) => a + w[2], 0)).toBeLessThanOrEqual(total + 1e-6);
  });
});
