// E5.2c record→replay core (bead wf-root-guzez.6.3.3): PURE flight
// recording + replay-identity checking + ghost lookup. The identity
// claim is the ENGINE's chained digest (org.frankensim.wf.sim-digest.v1)
// — equal scenario + equal engine ⇒ equal digest, and any divergence
// is a typed verdict, never a silent overlay of two different flights.
// (The full scrubber/ABComparisonReceiptV1 machinery is E6.1; this is
// the milestone ghost slice.)

import { PAYLOAD_F64S, PAYLOAD_F64S_V1, type ScenarioInit } from "./protocol.ts";
import { decodeSnapshot, type SimSnapshot } from "./snapshotView.ts";

export const RECORDING_SCHEMA_V1 = "org.frankensim.wf.flight-recording.v1";
export const RECORDING_SCHEMA = "org.frankensim.wf.flight-recording.v2";

/** Scenario identity fields (bigints carried as decimal strings). */
export interface RecordedScenario {
  readonly seed: string;
  readonly rhoKgM3: number;
  readonly headwindMps: number;
  readonly mode: number;
  readonly member: number;
  readonly railLengthM: number;
  readonly maxTicks: string;
  readonly assist: boolean;
  readonly catapult: boolean;
}

export interface FlightRecording {
  readonly schema: typeof RECORDING_SCHEMA | typeof RECORDING_SCHEMA_V1;
  readonly scenario: RecordedScenario;
  readonly runIntentId: string;
  readonly tick0Digest: string;
  /** Terminal phase word + engine digest at the end of the run. */
  readonly terminalPhase: string;
  readonly finalDigest: string;
  /** Per-snapshot ticks (parallel to `frames`). */
  readonly ticks: readonly number[];
  /** Flat payload transcript: v2 uses 14 words; legacy v1 uses 12. */
  readonly frames: readonly number[];
}

/** Payload stride carried by a recording schema. Stored v1 flights
 * retain their 12-word frames; v2 appends roll and heading. */
export function recordingPayloadWords(
  recording: Pick<FlightRecording, "schema">,
): typeof PAYLOAD_F64S | typeof PAYLOAD_F64S_V1 {
  return recording.schema === RECORDING_SCHEMA_V1 ? PAYLOAD_F64S_V1 : PAYLOAD_F64S;
}

export function scenarioToRecorded(s: ScenarioInit): RecordedScenario {
  return {
    seed: s.seed.toString(),
    rhoKgM3: s.rhoKgM3,
    headwindMps: s.headwindMps,
    mode: s.mode,
    member: s.member,
    railLengthM: s.railLengthM,
    maxTicks: s.maxTicks.toString(),
    assist: s.assist,
    catapult: s.catapult,
  };
}

export function recordedToScenario(r: RecordedScenario): ScenarioInit {
  return {
    seed: BigInt(r.seed),
    rhoKgM3: r.rhoKgM3,
    headwindMps: r.headwindMps,
    mode: r.mode,
    member: r.member,
    railLengthM: r.railLengthM,
    maxTicks: BigInt(r.maxTicks),
    assist: r.assist,
    catapult: r.catapult,
  };
}

/** Streaming recorder: append monotone-tick payloads, then seal. */
export class FlightRecorder {
  private readonly ticks: number[] = [];
  private readonly frames: number[] = [];

  append(tick: number, payload: Float64Array): void {
    if (payload.length !== PAYLOAD_F64S) {
      throw new RangeError(`payload ${payload.length} != v2 ${PAYLOAD_F64S}`);
    }
    const last = this.ticks[this.ticks.length - 1];
    if (last !== undefined && tick <= last) {
      throw new RangeError(`non-monotone tick ${tick} after ${last}`);
    }
    this.ticks.push(tick);
    for (let i = 0; i < PAYLOAD_F64S; i += 1) {
      this.frames.push(payload[i]!);
    }
  }

  frameCount(): number {
    return this.ticks.length;
  }

  seal(meta: {
    scenario: ScenarioInit;
    runIntentId: string;
    tick0Digest: string;
    terminalPhase: string;
    finalDigest: string;
  }): FlightRecording {
    if (this.ticks.length === 0) {
      throw new RangeError("empty recording cannot be sealed");
    }
    return {
      schema: RECORDING_SCHEMA,
      scenario: scenarioToRecorded(meta.scenario),
      runIntentId: meta.runIntentId,
      tick0Digest: meta.tick0Digest,
      terminalPhase: meta.terminalPhase,
      finalDigest: meta.finalDigest,
      ticks: [...this.ticks],
      frames: [...this.frames],
    };
  }
}

export type ReplayVerdict =
  | { readonly kind: "identical"; readonly digest: string }
  | {
      readonly kind: "diverged";
      readonly expectedDigest: string;
      readonly observedDigest: string;
    };

/** Compare a recording's sealed digest against a re-run's digest. */
export function replayVerdict(recording: FlightRecording, observedDigest: string): ReplayVerdict {
  return recording.finalDigest === observedDigest
    ? { kind: "identical", digest: observedDigest }
    : {
        kind: "diverged",
        expectedDigest: recording.finalDigest,
        observedDigest,
      };
}

/** Fail-closed parse of a serialized recording (hostile-twin hardened). */
export function parseRecording(json: string): FlightRecording | { readonly error: string } {
  let value: unknown;
  try {
    value = JSON.parse(json);
  } catch (e) {
    return { error: `not JSON: ${String(e)}` };
  }
  if (typeof value !== "object" || value === null) {
    return { error: "not an object" };
  }
  const r = value as Record<string, unknown>;
  if (r.schema !== RECORDING_SCHEMA && r.schema !== RECORDING_SCHEMA_V1) {
    return { error: `unsupported schema ${String(r.schema)}` };
  }
  if (
    typeof r.runIntentId !== "string" ||
    typeof r.tick0Digest !== "string" ||
    typeof r.terminalPhase !== "string" ||
    typeof r.finalDigest !== "string" ||
    !/^[0-9a-f]{64}$/.test(r.finalDigest)
  ) {
    return { error: "identity fields missing or malformed" };
  }
  if (!Array.isArray(r.ticks) || !Array.isArray(r.frames)) {
    return { error: "transcript arrays missing" };
  }
  const payloadWords = r.schema === RECORDING_SCHEMA_V1 ? PAYLOAD_F64S_V1 : PAYLOAD_F64S;
  if (r.frames.length !== r.ticks.length * payloadWords) {
    return { error: `frames ${r.frames.length} != ticks ${r.ticks.length} * ${payloadWords}` };
  }
  if (r.ticks.length === 0) {
    return { error: "empty transcript" };
  }
  for (let i = 1; i < r.ticks.length; i += 1) {
    if (!(Number(r.ticks[i]) > Number(r.ticks[i - 1]))) {
      return { error: `non-monotone ticks at index ${i}` };
    }
  }
  for (const v of r.frames as unknown[]) {
    if (typeof v !== "number" || !Number.isFinite(v)) {
      return { error: "non-finite frame value" };
    }
  }
  return value as FlightRecording;
}

/**
 * Ghost lookup: the recorded snapshot at the LATEST recorded tick ≤ the
 * live tick (hold-last semantics past the end — the ghost freezes at
 * its terminal state; it never extrapolates a flight that ended).
 */
export function ghostAt(recording: FlightRecording, liveTick: number): SimSnapshot | null {
  const ticks = recording.ticks;
  if (liveTick < ticks[0]!) {
    return null;
  }
  // Binary search: greatest i with ticks[i] <= liveTick.
  let lo = 0;
  let hi = ticks.length - 1;
  while (lo < hi) {
    const mid = (lo + hi + 1) >> 1;
    if (ticks[mid]! <= liveTick) {
      lo = mid;
    } else {
      hi = mid - 1;
    }
  }
  const payloadWords = recordingPayloadWords(recording);
  const payload = new Float64Array(payloadWords);
  for (let i = 0; i < payloadWords; i += 1) {
    payload[i] = recording.frames[lo * payloadWords + i]!;
  }
  return decodeSnapshot(ticks[lo]!, payload);
}
