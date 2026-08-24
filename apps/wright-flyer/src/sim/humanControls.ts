// E5.3a authentic mechanical controls (bead wf-root-guzez.6.4): PURE
// human-input plumbing between the transducer (input.ts, keyboard-rate
// — devices cannot claim to measure force, plan law) and the engine's
// physical units, plus the ApplyNextEligibleTickAndFlag hold law and
// the latency-decomposition record. No DOM, no Worker: headless-tested.

import { quantizeControl, type AppliedEvent } from "../transport/inputClock.ts";
import type { PilotCommand } from "../input.ts";

/** Engine lever-force cap [N] (mirror of simloop's ±220 clamp — the
 * command maps to the pilot's PUSH/PULL on the canard lever, and full
 * scale is exactly the engine's admitted authority). */
export const MAX_LEVER_FORCE_N = 220;
/** Engine warp-command cap [rad] (mirror of simloop's ±0.148 clamp =
 * the ±8.5° verified warp limit). */
export const MAX_WARP_RAD = 0.148;

/** Physical control values for one engine tick. */
export interface PhysicalControl {
  readonly leverForceN: number;
  readonly warpCmdRad: number;
}

/** Map the quantized pilot command to engine units. The command is
 * re-quantized defensively (the 1/4096 grid is the trace identity —
 * a value off the grid here would silently fork replay identity). */
export function toPhysical(cmd: PilotCommand): PhysicalControl {
  return {
    leverForceN: quantizeControl(cmd.canard) * MAX_LEVER_FORCE_N,
    warpCmdRad: quantizeControl(cmd.warp) * MAX_WARP_RAD,
  };
}

/**
 * The worker-side hold law for Human mode. The ENGINE requires input
 * every tick (control-input-missing is its law); the pilot's device
 * samples arrive slower than 120 Hz. Resolution (declared, mirrors
 * LateInputPolicy::ApplyNextEligibleTickAndFlag):
 *   - before the FIRST admitted control, the sim WAITS (no step, no
 *     refusal spam — the run starts when the pilot touches the stick);
 *   - each admitted control applies from the next eligible tick and
 *     HOLDS (zero-order) until superseded;
 *   - every admission is receipted with requested/applied/late-by.
 */
export class ControlHold {
  private held: PhysicalControl | null = null;
  private pendingValue: PhysicalControl | null = null;
  private pendingFromTick = 0;
  private readonly receipts: {
    sequence: number;
    appliedTick: number;
    lateByTicks: number;
  }[] = [];

  /** Admit a control for application at the NEXT ELIGIBLE tick
   * (`currentTick + 1`). `requestedTick` is retained only as the
   * wall-clock target for the lateness receipt — never as the
   * application tick. */
  admit(
    sequence: number,
    value: PhysicalControl,
    requestedTick: number,
    currentTick: number,
  ): { appliedTick: number; lateByTicks: number } {
    if (!Number.isFinite(value.leverForceN) || !Number.isFinite(value.warpCmdRad)) {
      throw new RangeError("non-finite control never enters the hold");
    }
    const nextEligible = currentTick + 1;
    // ApplyNextEligibleTickAndFlag — the declared law of this class.
    // Applying at the future wall-clock target instead DEADLOCKS the
    // paused engine: ticks advance only by stepping (H-2), so a first
    // touch requested even one tick past `nextEligible` could never come
    // due and every run would wait for its first step forever. Found
    // live by the E6.4c e2e row (bead frankensim-nty3a).
    const appliedTick = nextEligible;
    const lateByTicks = appliedTick - requestedTick;
    this.pendingValue = value;
    this.pendingFromTick = appliedTick;
    const receipt = { sequence, appliedTick, lateByTicks };
    this.receipts.push(receipt);
    return receipt;
  }

  /** The control the engine gets at `tick`, or null while waiting for
   * the first input (the sim does not step). */
  valueAt(tick: number): PhysicalControl | null {
    if (this.pendingValue !== null && tick >= this.pendingFromTick) {
      this.held = this.pendingValue;
      this.pendingValue = null;
    }
    return this.held;
  }

  receiptLog(): readonly { sequence: number; appliedTick: number; lateByTicks: number }[] {
    return this.receipts;
  }
}

/** The §6 latency decomposition stages, one record per control sample.
 * All times are main-thread monotonic ms; absent stages stay null so a
 * gap in measurement is visible, never a fake zero. */
export interface LatencyRecord {
  readonly sequence: number;
  /** Device event timestamp. */
  readonly deviceMs: number;
  /** Quantized physical command carried by this sample. The pump
   * heartbeats NEUTRAL commands by design, so these fields separate
   * real pilot actuation from the idle stream (e2e efficacy gate). */
  readonly leverN: number;
  readonly warpRad: number;
  /** postMessage to the worker. */
  readonly sentMs: number;
  /** Ack received (worker admitted the control). */
  ackMs: number | null;
  /** Tick the control applied at (from the ack). */
  appliedTick: number | null;
  lateByTicks: number | null;
  /** First snapshot with tick >= appliedTick arrived. */
  publishedMs: number | null;
  /** First rAF present after that snapshot. */
  presentedMs: number | null;
}

/** Format the decomposition line (JSONL, e2e logging contract). */
export function latencyLine(r: LatencyRecord): string {
  const seg = (a: number | null, b: number | null): number | null =>
    a !== null && b !== null ? Number((b - a).toFixed(2)) : null;
  return JSON.stringify({
    suite: "wf-input-latency",
    seq: r.sequence,
    device_to_sent_ms: seg(r.deviceMs, r.sentMs),
    sent_to_ack_ms: seg(r.sentMs, r.ackMs),
    ack_to_published_ms: seg(r.ackMs, r.publishedMs),
    published_to_present_ms: seg(r.publishedMs, r.presentedMs),
    device_to_present_ms: seg(r.deviceMs, r.presentedMs),
    applied_tick: r.appliedTick,
    late_by_ticks: r.lateByTicks,
    lever_n: r.leverN,
    warp_rad: r.warpRad,
  });
}

/** Bounded latency ledger: tracks in-flight records, completes them as
 * acks/snapshots/frames arrive, emits finished lines via `emit`. */
export class LatencyLedger {
  private readonly inflight = new Map<number, LatencyRecord>();
  private readonly emit: (line: string) => void;
  private readonly cap: number;

  constructor(emit: (line: string) => void, cap = 256) {
    this.emit = emit;
    this.cap = cap;
  }

  sent(sequence: number, deviceMs: number, sentMs: number, leverN: number, warpRad: number): void {
    if (this.inflight.size >= this.cap) {
      // Bounded (logging contract): drop the OLDEST, loudly.
      const oldest = this.inflight.keys().next().value;
      if (oldest !== undefined) {
        this.inflight.delete(oldest);
        this.emit(JSON.stringify({ suite: "wf-input-latency", dropped_seq: oldest }));
      }
    }
    this.inflight.set(sequence, {
      sequence,
      deviceMs,
      sentMs,
      leverN,
      warpRad,
      ackMs: null,
      appliedTick: null,
      lateByTicks: null,
      publishedMs: null,
      presentedMs: null,
    });
  }

  acked(sequence: number, appliedTick: number, lateByTicks: number, nowMs: number): void {
    const r = this.inflight.get(sequence);
    if (r !== undefined) {
      r.ackMs = nowMs;
      r.appliedTick = appliedTick;
      r.lateByTicks = lateByTicks;
    }
  }

  /** A snapshot at `tick` arrived; then the next present completes any
   * record whose applied tick it covers. */
  published(tick: number, nowMs: number): void {
    for (const r of this.inflight.values()) {
      if (r.publishedMs === null && r.appliedTick !== null && tick >= r.appliedTick) {
        r.publishedMs = nowMs;
      }
    }
  }

  presented(nowMs: number): void {
    for (const [seq, r] of [...this.inflight.entries()]) {
      if (r.publishedMs !== null && r.presentedMs === null) {
        r.presentedMs = nowMs;
        this.emit(latencyLine(r));
        this.inflight.delete(seq);
      }
    }
  }

  inflightCount(): number {
    return this.inflight.size;
  }
}

// Re-export for worker convenience (the applied trace type is the
// inputClock's — one canonical identity, per the Round-5 boundary).
export type { AppliedEvent };
