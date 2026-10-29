// E5.2c record→replay core (bead wf-root-guzez.6.3.3): PURE flight
// recording + replay-identity checking + ghost lookup. The identity
// claim is the ENGINE's chained digest (org.frankensim.wf.sim-digest.v1)
// — equal scenario + engine + applied input ⇒ equal digest, and any divergence
// is a typed verdict, never a silent overlay of two different flights.
// (The full scrubber/ABComparisonReceiptV1 machinery is E6.1; this is
// the milestone ghost slice.)
// v5 closes the human-input path. v1–v4 human flights remain ghost-only;
// their device commands cannot be reconstructed from render snapshots.

import {
  PAYLOAD_F64S, PAYLOAD_F64S_V1, type ScenarioInit, type SpatialWindInit, spatialWindError,
  type MaterialAirInit, materialAirError, MODE_HUMAN, type ReplayInitialization,
} from "./protocol.ts";
import { decodeSnapshot, type SimSnapshot } from "./snapshotView.ts";
import { appliedControlTraceError, copyAppliedControls, MAX_CONTROL_TICKS, type AppliedControlTrace } from "./appliedControls.ts";

export const RECORDING_SCHEMA_V1 = "org.frankensim.wf.flight-recording.v1";
export const RECORDING_SCHEMA = "org.frankensim.wf.flight-recording.v2";
/** Spatial recordings require a reader that preserves their force model. */
export const RECORDING_SCHEMA_SPATIAL = "org.frankensim.wf.flight-recording.v3";
/** Material recordings require a reader that retains the gas-law inputs. */
export const RECORDING_SCHEMA_MATERIAL = "org.frankensim.wf.flight-recording.v4";
/** Human replay requires the controls actually consumed by the native engine. */
export const RECORDING_SCHEMA_HUMAN = "org.frankensim.wf.flight-recording.v5";

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
  readonly spatialWind?: SpatialWindInit;
  readonly materialAir?: MaterialAirInit;
}

export interface FlightRecording {
  readonly schema: typeof RECORDING_SCHEMA | typeof RECORDING_SCHEMA_V1 |
    typeof RECORDING_SCHEMA_SPATIAL | typeof RECORDING_SCHEMA_MATERIAL | typeof RECORDING_SCHEMA_HUMAN;
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
  readonly appliedControls?: AppliedControlTrace;
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
    ...(s.spatialWind === undefined ? {} : {
      spatialWind: { ...s.spatialWind, correlationTimesS: [...s.spatialWind.correlationTimesS] },
    }),
    ...(s.materialAir === undefined ? {} : { materialAir: { ...s.materialAir } }),
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
    ...(r.spatialWind === undefined ? {} : {
      spatialWind: { ...r.spatialWind, correlationTimesS: [...r.spatialWind.correlationTimesS] },
    }),
    ...(r.materialAir === undefined ? {} : { materialAir: { ...r.materialAir } }),
  };
}

/** Streaming recorder: append monotone-tick payloads, then seal. */
export class FlightRecorder {
  private readonly ticks: number[] = [];
  private readonly frames: number[] = [];

  append(tick: number, payload: Float64Array): void {
    if (!Number.isSafeInteger(tick) || tick < 1 || tick > MAX_CONTROL_TICKS ||
        payload.some((value) => !Number.isFinite(value))) {
      throw new RangeError("snapshot requires a bounded integer tick and finite payload");
    }
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

  /** The SAB reader may already have sampled this exact terminal tick. */
  appendTerminal(tick: number, payload: Float64Array): void {
    if (this.ticks[this.ticks.length - 1] !== tick) {
      this.append(tick, payload);
      return;
    }
    const start = this.frames.length - PAYLOAD_F64S;
    if (payload.length !== PAYLOAD_F64S || payload.some((v, i) => !Object.is(v, this.frames[start + i]))) {
      throw new RangeError("terminal payload differs from the recorded native snapshot");
    }
  }

  seal(meta: {
    scenario: ScenarioInit;
    runIntentId: string;
    tick0Digest: string;
    terminalPhase: string;
    finalDigest: string;
    appliedControls?: AppliedControlTrace;
  }): FlightRecording {
    if (this.ticks.length === 0) {
      throw new RangeError("empty recording cannot be sealed");
    }
    const appliedControls = meta.appliedControls === undefined ? undefined : copyAppliedControls(meta.appliedControls);
    if (meta.scenario.mode === MODE_HUMAN) {
      if (appliedControls === undefined || appliedControls.endTick !== this.ticks[this.ticks.length - 1]) {
        throw new RangeError("human recording requires the complete applied controls through its terminal tick");
      }
    } else if (appliedControls !== undefined) {
      throw new RangeError("applied human controls cannot label a nonhuman recording");
    }
    return {
      schema: meta.scenario.mode === MODE_HUMAN ? RECORDING_SCHEMA_HUMAN :
        meta.scenario.materialAir !== undefined ? RECORDING_SCHEMA_MATERIAL :
        meta.scenario.spatialWind === undefined ? RECORDING_SCHEMA : RECORDING_SCHEMA_SPATIAL,
      scenario: scenarioToRecorded(meta.scenario),
      runIntentId: meta.runIntentId,
      tick0Digest: meta.tick0Digest,
      terminalPhase: meta.terminalPhase,
      finalDigest: meta.finalDigest,
      ticks: [...this.ticks],
      frames: [...this.frames],
      ...(appliedControls === undefined ? {} : { appliedControls }),
    };
  }
}

/** Legacy human recordings remain usable as ghosts, but cannot recreate input. */
export function replayInitialization(recording: FlightRecording): ReplayInitialization | { readonly error: string } {
  if (![recording.runIntentId, recording.tick0Digest].every((id) => /^[0-9a-f]{64}$/.test(id))) {
    return { error: "replay requires native run and tick-zero identities" };
  }
  if (recording.scenario.mode === MODE_HUMAN) {
    const trace = recording.appliedControls;
    if (recording.schema !== RECORDING_SCHEMA_HUMAN || trace === undefined) {
      return { error: "legacy human recording has no applied controls; start a new flight to record replayable input" };
    }
    const error = appliedControlTraceError(trace);
    if (error !== null) return { error };
    if (trace.endTick !== recording.ticks[recording.ticks.length - 1]) return { error: "applied trace does not reach the terminal tick" };
    return { runIntentId: recording.runIntentId, tick0Digest: recording.tick0Digest, appliedControls: copyAppliedControls(trace) };
  }
  if (recording.appliedControls !== undefined) return { error: "nonhuman recording contains human controls" };
  return { runIntentId: recording.runIntentId, tick0Digest: recording.tick0Digest };
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
  if (r.schema !== RECORDING_SCHEMA && r.schema !== RECORDING_SCHEMA_V1 &&
      r.schema !== RECORDING_SCHEMA_SPATIAL && r.schema !== RECORDING_SCHEMA_MATERIAL &&
      r.schema !== RECORDING_SCHEMA_HUMAN) {
    return { error: `unsupported schema ${String(r.schema)}` };
  }
  if (typeof r.scenario !== "object" || r.scenario === null || Array.isArray(r.scenario)) {
    return { error: "scenario missing or malformed" };
  }
  const scenario = r.scenario as Record<string, unknown>;
  for (const key of ["rhoKgM3", "headwindMps", "railLengthM"]) {
    if (typeof scenario[key] !== "number" || !Number.isFinite(scenario[key])) return { error: `scenario ${key} must be finite` };
  }
  if (typeof scenario.seed !== "string" || !/^(0|[1-9][0-9]{0,19})$/.test(scenario.seed) ||
      BigInt(scenario.seed) > 18446744073709551615n ||
      typeof scenario.maxTicks !== "string" || !/^[1-9][0-9]{0,4}$/.test(scenario.maxTicks) ||
      Number(scenario.maxTicks) > MAX_CONTROL_TICKS ||
      ![0, 1, MODE_HUMAN].includes(scenario.mode as number) ||
      !Number.isSafeInteger(scenario.member) || Number(scenario.member) < 0 ||
      typeof scenario.assist !== "boolean" || typeof scenario.catapult !== "boolean") {
    return { error: "scenario seed, duration, mode, member or flags malformed" };
  }
  const humanSchema = r.schema === RECORDING_SCHEMA_HUMAN;
  const hasSpatialModel = r.schema === RECORDING_SCHEMA_SPATIAL || r.schema === RECORDING_SCHEMA_MATERIAL;
  if (!humanSchema && hasSpatialModel !== (scenario.spatialWind !== undefined)) {
    return { error: "spatial recordings require a spatial schema and their complete wind configuration" };
  }
  if (scenario.spatialWind !== undefined) {
    const error = spatialWindError(scenario.spatialWind);
    if (error !== null) return { error };
  }
  if ((!humanSchema && (r.schema === RECORDING_SCHEMA_MATERIAL) !== (scenario.materialAir !== undefined)) ||
      (scenario.materialAir !== undefined && scenario.spatialWind === undefined)) {
    return { error: "material recordings require the v4 schema and their complete air configuration" };
  }
  if (scenario.materialAir !== undefined) {
    const error = materialAirError(scenario.materialAir);
    if (error !== null) return { error };
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
  for (let i = 0; i < r.ticks.length; i += 1) {
    if (!Number.isSafeInteger(r.ticks[i]) || r.ticks[i] < 1 || r.ticks[i] > Number(scenario.maxTicks) ||
        (i > 0 && !(r.ticks[i] > r.ticks[i - 1]))) {
      return { error: `non-monotone ticks at index ${i}` };
    }
  }
  for (const v of r.frames as unknown[]) {
    if (typeof v !== "number" || !Number.isFinite(v)) {
      return { error: "non-finite frame value" };
    }
  }
  if (humanSchema) {
    if (scenario.mode !== MODE_HUMAN) return { error: "v5 records human applied input only" };
    const replay = replayInitialization(value as FlightRecording);
    if ("error" in replay) return replay;
  } else if (r.appliedControls !== undefined) {
    return { error: "applied controls require the v5 human schema" };
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
