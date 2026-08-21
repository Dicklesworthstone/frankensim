// E7.1b battery (bead wf-root-guzez.8.2): tracer laws. Deterministic
// checkpoint reconstruction (bitwise continuation equality); the
// 10-minute soak (72000 ticks) holds the retention cap; pathline /
// streakline / streamline VERIFIED DISTINCT on an unsteady flow and
// coincident on a steady one; snapshot-PAIR binding (bad pairs
// refuse); cancellation; caps at cap AND cap+1; no field data on the
// checkpoint (dense-history rejection witness).
// Repro: node --test test/tracers.test.ts

import assert from "node:assert/strict";
import { test } from "node:test";

import { integrateStreamlines } from "../src/fieldViz.ts";
import {
  MAX_POINTS_PER_TRACER,
  MAX_TRACERS,
  TracerService,
  type SnapshotPair,
} from "../src/tracers.ts";

function jlog(kase: string, payload: string): void {
  console.log(`{"suite":"wf-app-tracers","case":"${kase}",${payload}}`);
}

const DT = 1 / 120;

/** Unsteady fixture: u = (1, k·t, 0) — the classic distinctness flow. */
function unsteadyPair(tick: number): SnapshotPair {
  const u = (t: number) => (_p: readonly [number, number, number]) =>
    [1, 0.4 * t, 0] as const;
  return {
    tickA: tick,
    tickB: tick + 1,
    idA: `snap-${tick}`,
    idB: `snap-${tick + 1}`,
    samplerA: u(tick * DT),
    samplerB: u((tick + 1) * DT),
  };
}

function steadyPair(tick: number): SnapshotPair {
  const u = (_p: readonly [number, number, number]) => [1, 0.25, 0] as const;
  return {
    tickA: tick,
    tickB: tick + 1,
    idA: `snap-${tick}`,
    idB: `snap-${tick + 1}`,
    samplerA: u,
    samplerB: u,
  };
}

test("checkpoint reconstruction is deterministic (bitwise continuation)", () => {
  const a = new TracerService();
  assert.ok(a.release(0, [0, 0, 1]).ok);
  assert.ok(a.release(0, [0.5, 0, 1]).ok);
  for (let t = 0; t < 50; t += 1) assert.ok(a.advance(unsteadyPair(t), DT).ok);
  const cp = a.checkpoint();
  // Continue the original 30 more ticks…
  for (let t = 50; t < 80; t += 1) assert.ok(a.advance(unsteadyPair(t), DT).ok);
  // …and the restored twin the same 30 ticks.
  const b = TracerService.restore(cp);
  for (let t = 50; t < 80; t += 1) assert.ok(b.advance(unsteadyPair(t), DT).ok);
  const pa = a.pathline(0);
  const pb = b.pathline(0);
  assert.ok(pa.ok && pb.ok);
  if (pa.ok && pb.ok) assert.deepEqual(pa.value, pb.value, "bitwise reconstruction");
  // The checkpoint holds NO field data: its size depends on particle
  // count and trail cap, never on how many snapshots passed through.
  assert.ok(!JSON.stringify(cp).includes("sampler"));
  jlog("checkpoint", `"tracers":${cp.tracers.length}`);
});

test("10-minute soak holds the retention cap", () => {
  const s = new TracerService();
  for (let i = 0; i < 8; i += 1) assert.ok(s.release(0, [i, 0, 1]).ok);
  for (let t = 0; t < 72_000; t += 1) {
    const r = s.advance(steadyPair(t), DT);
    assert.ok(r.ok);
  }
  const points = s.retainedPoints();
  assert.ok(
    points <= 8 * (MAX_POINTS_PER_TRACER + 1),
    `soak retained ${points} points beyond the cap`
  );
  assert.ok(points > 8 * 32, "thinning must still RETAIN a usable trail");
  jlog("soak", `"ticks":72000,"retained_points":${points}`);
});

test("pathline, streakline, streamline are DISTINCT on unsteady flow, coincident on steady", () => {
  // Unsteady: release from the SAME point at successive ticks.
  const s = new TracerService();
  const origin = [0, 0, 1] as const;
  for (let t = 0; t < 40; t += 1) {
    assert.ok(s.release(t, origin).ok);
    assert.ok(s.advance(unsteadyPair(t), DT).ok);
  }
  const path = s.pathline(0);
  assert.ok(path.ok);
  const streak = s.streakline(origin);
  // Streamline at the FINAL instant (fieldViz — instantaneous field).
  const finalT = 40 * DT;
  const stream = integrateStreamlines(
    () => [1, 0.4 * finalT, 0] as const,
    [[origin[0], origin[1], origin[2]]],
    DT,
    39,
  );
  assert.ok(stream.ok);
  if (!path.ok || !stream.ok) return;
  // Compare the three curves' terminal y/x slopes: they must differ
  // pairwise on the unsteady flow (the classic textbook fact).
  const slopeOf = (pts: Float64Array | readonly number[]) => {
    const n = pts.length;
    return (pts[n - 2] ?? 0) / (pts[n - 3] ?? 1);
  };
  const sp = slopeOf(path.value);
  const sk = slopeOf(streak);
  const sl = slopeOf(stream.value[0]?.points ?? new Float64Array());
  assert.ok(Math.abs(sp - sk) > 1e-3, `pathline vs streakline: ${sp} vs ${sk}`);
  assert.ok(Math.abs(sp - sl) > 1e-3, `pathline vs streamline: ${sp} vs ${sl}`);
  assert.ok(Math.abs(sk - sl) > 1e-3, `streakline vs streamline: ${sk} vs ${sl}`);
  // Steady flow: the three coincide (up to integrator class).
  const st = new TracerService();
  for (let t = 0; t < 40; t += 1) {
    assert.ok(st.release(t, origin).ok);
    assert.ok(st.advance(steadyPair(t), DT).ok);
  }
  const pathS = st.pathline(0);
  const streakS = st.streakline(origin);
  assert.ok(pathS.ok);
  if (pathS.ok) {
    const a = slopeOf(pathS.value);
    const b = slopeOf(streakS);
    assert.ok(Math.abs(a - b) < 1e-9, `steady coincidence: ${a} vs ${b}`);
    assert.ok(Math.abs(a - 0.25) < 1e-9, "steady slope is the field's");
  }
  jlog("distinctness", `"pathline":${sp},"streakline":${sk},"streamline":${sl}`);
});

test("pair binding, cancellation, caps", () => {
  const s = new TracerService();
  assert.ok(s.release(0, [0, 0, 1]).ok);
  // Non-adjacent pair refuses; empty snapshot id refuses.
  const bad = s.advance({ ...steadyPair(0), tickB: 5 }, DT);
  assert.ok(!bad.ok && bad.refusal.code === "tracer-pair-invalid");
  const noId = s.advance({ ...steadyPair(0), idA: "" }, DT);
  assert.ok(!noId.ok && noId.refusal.code === "tracer-pair-invalid");
  const badDt = s.advance(steadyPair(0), 0);
  assert.ok(!badDt.ok && badDt.refusal.code === "tracer-dt-invalid");
  // Retroactive release refuses.
  assert.ok(s.advance(steadyPair(0), DT).ok);
  const retro = s.release(0, [0, 0, 1]);
  assert.ok(!retro.ok && retro.refusal.code === "tracer-release-invalid");
  // Tracer cap AND cap+1.
  const big = new TracerService();
  for (let i = 0; i < MAX_TRACERS; i += 1) {
    assert.ok(big.release(0, [i, 0, 1]).ok, `release ${i}`);
  }
  const over = big.release(0, [0, 0, 1]);
  assert.ok(!over.ok && over.refusal.code === "tracer-count-exceeded");
  // Cancellation: every later mutation refuses.
  big.cancel();
  const afterCancel = big.advance(steadyPair(0), DT);
  assert.ok(!afterCancel.ok && afterCancel.refusal.code === "tracer-cancelled");
  const releaseAfter = big.release(0, [0, 0, 1]);
  assert.ok(!releaseAfter.ok && releaseAfter.refusal.code === "tracer-cancelled");
  jlog("caps", `"max_tracers":${MAX_TRACERS}`);
});
