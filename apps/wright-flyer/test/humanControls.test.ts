// E5.3a battery (bead wf-root-guzez.6.4): physical mapping oracles at
// cap AND past-cap (engine authority is the full scale), the
// ApplyNextEligibleTickAndFlag hold law (wait-before-first, ZOH,
// late-by receipts, supersede), the latency ledger lifecycle incl. the
// bounded-drop path, and the quantization-grid identity (a value off
// the 1/4096 grid would silently fork replay identity — falsified).
// Repro: node --test test/humanControls.test.ts

import assert from "node:assert/strict";
import { test } from "node:test";

import {
  ControlHold,
  LatencyLedger,
  MAX_LEVER_FORCE_N,
  MAX_WARP_RAD,
  latencyLine,
  toPhysical,
} from "../src/sim/humanControls.ts";
import { quantizeControl } from "../src/transport/inputClock.ts";
import { NEUTRAL, stepCommand } from "../src/input.ts";

test("toPhysical: exact scale at full command; quantization is the identity grid", () => {
  const full = toPhysical({ canard: 1, warp: -1, mode: "keyboard-rate" });
  assert.equal(full.leverForceN, MAX_LEVER_FORCE_N, "full pull = engine cap exactly");
  assert.equal(full.warpCmdRad, -MAX_WARP_RAD, "full warp = ±8.5° cap exactly");
  const zero = toPhysical(NEUTRAL);
  assert.equal(zero.leverForceN, 0);
  assert.equal(zero.warpCmdRad, 0);
  // Off-grid input is re-quantized (the falsifier: a raw 0.3 is NOT on
  // the 1/4096 grid, so skipping the re-quantization would produce a
  // different force than the trace records).
  const off = toPhysical({ canard: 0.3, warp: 0, mode: "keyboard-rate" });
  assert.equal(off.leverForceN, quantizeControl(0.3) * MAX_LEVER_FORCE_N);
  assert.notEqual(off.leverForceN, 0.3 * MAX_LEVER_FORCE_N, "0.3 is off-grid");
  // Past cap clamps via the quantizer (cap AND cap+1 class).
  const over = toPhysical({ canard: 1.5, warp: 2, mode: "keyboard-rate" });
  assert.equal(over.leverForceN, MAX_LEVER_FORCE_N);
  assert.equal(over.warpCmdRad, MAX_WARP_RAD);
});

test("keyboard-rate transducer composes: slew then map stays in authority", () => {
  let cmd = NEUTRAL;
  for (let i = 0; i < 600; i += 1) {
    cmd = stepCommand(
      cmd,
      { canardUp: true, canardDown: false, warpLeft: false, warpRight: true, recenter: false },
      1 / 60,
    );
  }
  // 10 s of held key saturates the command; force saturates at the cap.
  const phys = toPhysical(cmd);
  assert.equal(phys.leverForceN, MAX_LEVER_FORCE_N);
  assert.equal(phys.warpCmdRad, MAX_WARP_RAD);
});

test("ControlHold: waits before first input, ZOH after, receipts exact", () => {
  const hold = new ControlHold();
  // The sim WAITS: no value at any tick before the first admission.
  assert.equal(hold.valueAt(1), null);
  assert.equal(hold.valueAt(50), null);
  // Admit at currentTick 10, requested 8 (device was early → late by 3:
  // applied = max(8, 11) = 11).
  const r1 = hold.admit(1, { leverForceN: 40, warpCmdRad: 0.01 }, 8, 10);
  assert.equal(r1.appliedTick, 11);
  assert.equal(r1.lateByTicks, 3);
  // Not yet applied at tick 10; applies from 11; holds through 200.
  assert.equal(hold.valueAt(10), null, "never applies before its tick");
  assert.equal(hold.valueAt(11)?.leverForceN, 40);
  assert.equal(hold.valueAt(200)?.leverForceN, 40, "zero-order hold");
  // On-time input: requested 202 > next eligible 201 → late by 0.
  const r2 = hold.admit(2, { leverForceN: -60, warpCmdRad: 0 }, 202, 200);
  assert.equal(r2.appliedTick, 202);
  assert.equal(r2.lateByTicks, 0);
  assert.equal(hold.valueAt(201)?.leverForceN, 40, "old value until applied tick");
  assert.equal(hold.valueAt(202)?.leverForceN, -60, "supersedes at its tick");
  // Non-finite refuses.
  assert.throws(
    () => hold.admit(3, { leverForceN: Number.NaN, warpCmdRad: 0 }, 300, 299),
    RangeError,
  );
  assert.equal(hold.receiptLog().length, 2);
});

test("latency ledger: full lifecycle emits one line with every stage", () => {
  const lines: string[] = [];
  const ledger = new LatencyLedger((l) => lines.push(l));
  ledger.sent(7, 1000, 1002);
  ledger.acked(7, 480, 1, 1010);
  ledger.published(479, 1015); // earlier tick: must NOT complete seq 7
  assert.equal(lines.length, 0);
  ledger.presented(1016);
  assert.equal(lines.length, 0, "not presented before published");
  ledger.published(480, 1020);
  ledger.presented(1024);
  assert.equal(lines.length, 1);
  const rec = JSON.parse(lines[0]!) as Record<string, unknown>;
  assert.equal(rec.suite, "wf-input-latency");
  assert.equal(rec.seq, 7);
  assert.equal(rec.device_to_sent_ms, 2);
  assert.equal(rec.sent_to_ack_ms, 8);
  assert.equal(rec.ack_to_published_ms, 10);
  assert.equal(rec.published_to_present_ms, 4);
  assert.equal(rec.device_to_present_ms, 24);
  assert.equal(rec.applied_tick, 480);
  assert.equal(rec.late_by_ticks, 1);
  assert.equal(ledger.inflightCount(), 0);
});

test("latency ledger: unmeasured stages stay null, never fake zeros", () => {
  const r = {
    sequence: 1,
    deviceMs: 100,
    sentMs: 101,
    ackMs: null,
    appliedTick: null,
    lateByTicks: null,
    publishedMs: null,
    presentedMs: null,
  };
  const line = JSON.parse(latencyLine(r)) as Record<string, unknown>;
  assert.equal(line.device_to_sent_ms, 1);
  assert.equal(line.sent_to_ack_ms, null);
  assert.equal(line.ack_to_published_ms, null);
  assert.equal(line.device_to_present_ms, null);
});

test("latency ledger: bounded — overflow drops the OLDEST, loudly", () => {
  const lines: string[] = [];
  const ledger = new LatencyLedger((l) => lines.push(l), 4);
  for (let seq = 1; seq <= 5; seq += 1) {
    ledger.sent(seq, seq * 10, seq * 10 + 1);
  }
  assert.equal(ledger.inflightCount(), 4, "cap held");
  const dropped = lines.map((l) => JSON.parse(l) as Record<string, unknown>);
  assert.equal(dropped.length, 1);
  assert.equal(dropped[0]!.dropped_seq, 1, "oldest dropped, and it is REPORTED");
});
