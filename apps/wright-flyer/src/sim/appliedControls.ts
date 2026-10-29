// The applied command at the native tick boundary is the replay authority.
// Device acknowledgements may be superseded before a tick and are not a trace.
import { ControlHold, type PhysicalControl } from "./humanControls.ts";
import { parseStepEnvelope } from "./engineFacade.ts";

export const MAX_CONTROL_TICKS = 72_000;
/** JSON numbers erase negative zero; native command digests need its bits. */
type RecordedNumber = number | "-0";
export interface AppliedControlTrace {
  readonly endTick: number;
  readonly changes: readonly {
    readonly tick: number;
    readonly leverForceN: RecordedNumber;
    readonly warpCmdRad: RecordedNumber;
  }[];
}

export function appliedControlTraceError(value: unknown): string | null {
  if (typeof value !== "object" || value === null) return "applied controls missing";
  const trace = value as AppliedControlTrace;
  if (!Number.isSafeInteger(trace.endTick) || trace.endTick < 1 || trace.endTick > MAX_CONTROL_TICKS ||
      !Array.isArray(trace.changes) || trace.changes.length === 0 || trace.changes.length > trace.endTick) {
    return "applied controls require a bounded nonempty tick extent";
  }
  let previous = 0;
  for (const change of trace.changes) {
    if (typeof change !== "object" || change === null || !Number.isSafeInteger(change.tick) ||
        change.tick <= previous || change.tick > trace.endTick ||
        (previous === 0 && change.tick !== 1)) return "applied controls must start at tick 1 and increase";
    for (const v of [change.leverForceN, change.warpCmdRad]) {
      if (v !== "-0" && (typeof v !== "number" || !Number.isFinite(v))) return "applied command must be finite";
    }
    previous = change.tick;
  }
  return null;
}

export function copyAppliedControls(trace: AppliedControlTrace): AppliedControlTrace {
  const error = appliedControlTraceError(trace);
  if (error !== null) throw new RangeError(error);
  return { endTick: trace.endTick, changes: trace.changes.map((change) => ({ ...change })) };
}

const encode = (value: number): RecordedNumber => Object.is(value, -0) ? "-0" : value;
const decode = (value: RecordedNumber): number => value === "-0" ? -0 : value;

/** Live acquisition and playback share exactly one successful-step boundary. */
export class AppliedControls {
  private readonly hold = new ControlHold();
  private readonly playback: AppliedControlTrace | undefined;
  private cursor = 0;
  private endTick = 0;
  private failed = false;
  private readonly changes: AppliedControlTrace["changes"][number][] = [];

  constructor(playback?: AppliedControlTrace) {
    this.playback = playback === undefined ? undefined : copyAppliedControls(playback);
  }

  admit(sequence: number, value: PhysicalControl, requestedTick: number, currentTick: number):
    { appliedTick: number; lateByTicks: number } | null {
    if (this.playback !== undefined || this.failed) return null;
    return this.hold.admit(sequence, value, requestedTick, currentTick);
  }

  step(tick: number, nativeStep: (value: PhysicalControl) => string):
    ReturnType<typeof parseStepEnvelope> | { readonly kind: "waiting" } {
    if (this.failed || tick !== this.endTick + 1 || tick > MAX_CONTROL_TICKS ||
        (this.playback !== undefined && tick > this.playback.endTick)) {
      this.failed = true;
      return { kind: "malformed", detail: "applied control tick outside the consecutive recorded schedule" };
    }
    let value: PhysicalControl | null;
    if (this.playback === undefined) {
      value = this.hold.valueAt(tick);
    } else {
      const changes = this.playback.changes;
      while (this.cursor + 1 < changes.length && changes[this.cursor + 1]!.tick <= tick) this.cursor += 1;
      const change = changes[this.cursor]!;
      value = { leverForceN: decode(change.leverForceN), warpCmdRad: decode(change.warpCmdRad) };
    }
    if (value === null) return { kind: "waiting" };
    let step: ReturnType<typeof parseStepEnvelope>;
    try {
      step = parseStepEnvelope(nativeStep(value));
    } catch (error) {
      this.failed = true;
      return { kind: "malformed", detail: `native step threw: ${String(error)}` };
    }
    if (step.kind !== "ok" || step.tick !== tick) {
      this.failed = true;
      return step.kind === "ok" ? { kind: "malformed", detail: "native tick differs from applied control tick" } : step;
    }
    const last = this.changes[this.changes.length - 1];
    const leverForceN = encode(value.leverForceN);
    const warpCmdRad = encode(value.warpCmdRad);
    if (last === undefined || !Object.is(last.leverForceN, leverForceN) || !Object.is(last.warpCmdRad, warpCmdRad)) {
      this.changes.push({ tick, leverForceN, warpCmdRad });
    }
    this.endTick = tick;
    return step;
  }

  seal(): AppliedControlTrace {
    if (this.failed) throw new RangeError("a failed native step cannot seal an applied trace");
    return copyAppliedControls({ endTick: this.endTick, changes: this.changes });
  }
}
