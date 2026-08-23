// Timing/input-semantics battery (bead wf-root-guzez.1.8, E0.8). Receipts:
// clock-offset recovery, deterministic quantization, requested/applied tick
// protocol with ApplyNextEligibleTickAndFlag, canonical-vs-acquisition
// trace separation, backlog bound + re-anchor, and metric separation
// (service vs lateness vs backlog).

import assert from "node:assert/strict";
import { test } from "node:test";
import {
  InputScheduler,
  estimateClockOffsetMs,
  quantizeControl,
} from "../src/transport/inputClock.ts";
import { TickScheduler, percentile } from "../src/transport/schedule.ts";

const jlog = (obj: object): void =>
  console.log(JSON.stringify({ suite: "wf-timing", ...obj }));

test("clock sync recovers a known offset via the min-RTT midpoint", () => {
  // Remote clock = local + 250 ms. Asymmetric noise on most samples; one
  // clean low-RTT sample must dominate.
  const samples = [
    { localSentMs: 100, remoteMs: 100 + 250 + 9, localReceivedMs: 130 },
    { localSentMs: 200, remoteMs: 202 + 250, localReceivedMs: 204 },
    { localSentMs: 300, remoteMs: 300 + 250 + 22, localReceivedMs: 360 },
  ];
  const offset = estimateClockOffsetMs(samples);
  assert.ok(Math.abs(offset - 250) < 1.0, `offset ${offset} != 250±1`);
  assert.throws(() => estimateClockOffsetMs([]), /at least one/);
  jlog({ case: "clock-sync", offset });
});

test("control quantization is deterministic, clamped, and finite-safe", () => {
  assert.equal(quantizeControl(0.5), quantizeControl(0.5));
  assert.equal(quantizeControl(2.0), 1);
  assert.equal(quantizeControl(-2.0), -1);
  assert.equal(quantizeControl(Number.NaN), 0);
  assert.equal(quantizeControl(1 / 3), Math.round((1 / 3) * 4096) / 4096);
});

test("on-time inputs hit their requested tick; late inputs flag ApplyNextEligible", () => {
  const tickMs = 1000 / 120;
  const sched = new InputScheduler(tickMs, /*epoch*/ 0, /*lead*/ 1);
  // Sampled during tick-window 10 → requested = 11 + wait: raw=ceil(t/tickMs), lead 1.
  const sampleMs = 10.2 * tickMs;
  const requested = sched.requestedTickFor(sampleMs); // ceil(10.2)+1 = 12
  assert.equal(requested, 12);
  // Sim is at tick 5 → next eligible 6 < requested → applies at requested, on time.
  const onTime = sched.admit(
    { channel: 0, deviceSampleMs: sampleMs, quantizedValue: 0.25, sequence: 1 },
    5,
  );
  assert.equal(onTime.appliedTick, 12);
  // Sim already at tick 20 → the input is LATE → next eligible 21, flagged.
  const late = sched.admit(
    { channel: 0, deviceSampleMs: sampleMs, quantizedValue: 0.5, sequence: 2 },
    20,
  );
  assert.equal(late.appliedTick, 21);
  const acq = sched.acquisitionTrace();
  assert.equal(acq[0]!.lateByTicks, 0);
  assert.equal(acq[1]!.lateByTicks, 21 - 12);
  // Canonical trace carries NO acquisition-clock fields (Round-5 boundary).
  const canonical = sched.appliedTrace();
  assert.deepEqual(Object.keys(canonical[0]!).sort(), [
    "appliedTick",
    "channel",
    "ordinalWithinTick",
    "quantizedValue",
  ]);
  jlog({ case: "late-input", lateBy: acq[1]!.lateByTicks });
});

test("two inputs on one tick receive distinct ordinals (canonical determinism)", () => {
  const sched = new InputScheduler(1, 0, 1);
  const a = sched.admit({ channel: 0, deviceSampleMs: 100, quantizedValue: 1, sequence: 1 }, 200);
  const b = sched.admit({ channel: 1, deviceSampleMs: 100, quantizedValue: -1, sequence: 2 }, 200);
  assert.equal(a.appliedTick, b.appliedTick);
  assert.equal(a.ordinalWithinTick, 0);
  assert.equal(b.ordinalWithinTick, 1);
});

test("scheduler: exact tick count, bounded burst, re-anchor on long stall", () => {
  const tickMs = 10;
  const sched = new TickScheduler(tickMs, 0, /*maxCatchup*/ 3);
  const ran: number[] = [];
  const fakeNow = (): number => 0;
  // Advance 35 ms in small steps: ticks due at 0,10,20,30 → 4 ticks total.
  sched.pump(0, (t) => ran.push(t), fakeNow);
  sched.pump(15, (t) => ran.push(t), fakeNow);
  sched.pump(35, (t) => ran.push(t), fakeNow);
  assert.deepEqual(ran, [1, 2, 3, 4]);
  assert.equal(sched.metrics.reanchors, 0);
  // A 500 ms stall (50 ticks of backlog) must re-anchor, not burst.
  const before = sched.metrics.ticksRun;
  sched.pump(535, (t) => ran.push(t), fakeNow);
  const burst = sched.metrics.ticksRun - before;
  assert.ok(burst <= 3 + 1, `burst ${burst} exceeded the bound`);
  assert.equal(sched.metrics.reanchors, 1, "long stall must re-anchor");
  assert.ok(sched.metrics.maxBacklogObserved >= 50);
  jlog({
    case: "scheduler-bound",
    burst,
    reanchors: sched.metrics.reanchors,
    maxBacklog: sched.metrics.maxBacklogObserved,
  });
});

test("metric separation: service, lateness, and backlog are independent series", () => {
  const sched = new TickScheduler(10, 0, 3);
  let clock = 0;
  const nowFn = (): number => {
    clock += 2; // every service call "costs" 2 ms
    return clock;
  };
  sched.pump(12, () => true, nowFn); // one tick 12 ms late? due 0 and 10 → 2 ticks
  assert.equal(sched.metrics.latenessMs.length, sched.metrics.serviceMs.length);
  assert.ok(sched.metrics.latenessMs[0]! >= 12 - 0.001, "start lateness measured");
  assert.ok(sched.metrics.serviceMs[0]! > 0, "service time measured");
  assert.equal(percentile([1, 2, 3, 4], 0.5), 3);
  assert.equal(percentile([], 0.99), 0);
});
