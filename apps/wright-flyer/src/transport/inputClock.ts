// Input-clock semantics (plan §6 item 2, Round-2/Round-4/Round-5): how a
// main-thread device event becomes a simulation tick. Bead
// wf-root-guzez.1.8 (E0.8).
//
// - Main/worker monotonic clocks are synchronized by a ping exchange
//   (midpoint offset estimate) with periodic drift checks.
// - Every device sample is deterministically QUANTIZED, then targeted at a
//   requested tick with a minimum lead.
// - LateInputPolicy::ApplyNextEligibleTickAndFlag — NO interactive
//   rollback: a late input applies at the next eligible tick and the trace
//   records requested_tick, applied_tick, and late_by_ticks.
// - The CANONICAL applied trace (channel, applied_tick, ordinal, value) is
//   kept separate from the acquisition attachment (device times, requested
//   ticks, lateness) — one InputTraceId preimage, per the Round-5 boundary.

export interface ClockSyncSample {
  /** Local monotonic time when the ping was sent [ms]. */
  readonly localSentMs: number;
  /** Remote monotonic time stamped inside the pong [ms]. */
  readonly remoteMs: number;
  /** Local monotonic time when the pong arrived [ms]. */
  readonly localReceivedMs: number;
}

/**
 * Midpoint offset estimate: remote ≈ local + offset. Uses the
 * minimum-round-trip sample (least queueing distortion) per NTP practice.
 */
export function estimateClockOffsetMs(samples: readonly ClockSyncSample[]): number {
  if (samples.length === 0) {
    throw new Error("clock sync requires at least one sample");
  }
  let best = samples[0]!;
  let bestRtt = best.localReceivedMs - best.localSentMs;
  for (const s of samples) {
    const rtt = s.localReceivedMs - s.localSentMs;
    if (rtt < bestRtt) {
      best = s;
      bestRtt = rtt;
    }
  }
  const midpoint = (best.localSentMs + best.localReceivedMs) / 2;
  return best.remoteMs - midpoint;
}

/** Deterministic control quantization: fixed-point 1/4096 steps in [-1, 1]. */
export function quantizeControl(value: number): number {
  if (!Number.isFinite(value)) {
    return 0;
  }
  const clamped = Math.min(1, Math.max(-1, value));
  return Math.round(clamped * 4096) / 4096;
}

export interface InputPacket {
  readonly channel: number;
  readonly deviceSampleMs: number; // main-thread monotonic
  readonly quantizedValue: number;
  readonly sequence: number;
}

/** Canonical applied event (the InputTraceV1 domain). */
export interface AppliedEvent {
  readonly channel: number;
  readonly appliedTick: number;
  readonly ordinalWithinTick: number;
  readonly quantizedValue: number;
}

/** Acquisition attachment row (NEVER part of the canonical trace hash). */
export interface AcquisitionRow {
  readonly sequence: number;
  readonly deviceSampleMs: number;
  readonly requestedTick: number;
  readonly appliedTick: number;
  readonly lateByTicks: number;
}

export class InputScheduler {
  private readonly tickMs: number;
  private readonly epochMs: number;
  private readonly minLeadTicks: number;
  private readonly applied: AppliedEvent[] = [];
  private readonly acquisition: AcquisitionRow[] = [];
  private ordinalTick = -1;
  private ordinal = 0;

  constructor(tickMs: number, simEpochMs: number, minLeadTicks = 1) {
    this.tickMs = tickMs;
    this.epochMs = simEpochMs;
    this.minLeadTicks = minLeadTicks;
  }

  /** Target tick for a device sample under the minimum-lead rule. */
  requestedTickFor(deviceSampleMs: number): number {
    const raw = Math.ceil((deviceSampleMs - this.epochMs) / this.tickMs);
    return raw + this.minLeadTicks;
  }

  /**
   * Admit a packet at simulation time `currentTick` (the next tick to run
   * is `currentTick + 1`). Applies ApplyNextEligibleTickAndFlag.
   */
  admit(packet: InputPacket, currentTick: number): AppliedEvent {
    const requested = this.requestedTickFor(packet.deviceSampleMs);
    const nextEligible = currentTick + 1;
    const applied = Math.max(requested, nextEligible);
    const lateBy = applied - requested;
    if (applied !== this.ordinalTick) {
      this.ordinalTick = applied;
      this.ordinal = 0;
    }
    const event: AppliedEvent = {
      channel: packet.channel,
      appliedTick: applied,
      ordinalWithinTick: this.ordinal,
      quantizedValue: packet.quantizedValue,
    };
    this.ordinal += 1;
    this.applied.push(event);
    this.acquisition.push({
      sequence: packet.sequence,
      deviceSampleMs: packet.deviceSampleMs,
      requestedTick: requested,
      appliedTick: applied,
      lateByTicks: lateBy,
    });
    return event;
  }

  /** The canonical applied trace (InputTraceV1 events). */
  appliedTrace(): readonly AppliedEvent[] {
    return this.applied;
  }

  /** The acquisition attachment (replay evidence, not identity). */
  acquisitionTrace(): readonly AcquisitionRow[] {
    return this.acquisition;
  }
}
