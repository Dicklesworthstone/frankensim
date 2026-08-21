// Persistent Lagrangian tracer service (bead wf-root-guzez.8.2,
// E7.1b). Tick-addressed release, temporal interpolation over the
// ordered snapshot PAIR (Round-4 S-02: tracers bind two snapshot
// ids, never a raw tick), compact trajectory retention (deterministic
// stride-doubling thinning — a dense field-history ring is REJECTED:
// this service never stores a field, only its own particle states),
// checkpoint/replay reconstruction, cancellation, memory caps.
//
// Terminology, verified distinct by the battery on an unsteady flow:
//   - PATHLINE: one particle's trajectory through time;
//   - STREAKLINE: the locus of ALL particles released from one point
//     at successive ticks (what smoke shows);
//   - STREAMLINE: an instantaneous field line — that is
//     fieldViz.integrateStreamlines' job, deliberately NOT here.

export interface TracerRefusal {
  readonly code: string;
  readonly message: string;
  readonly rankedRepairs: readonly string[];
}

export type TracerResult<T> = { ok: true; value: T } | { ok: false; refusal: TracerRefusal };

function refuse<T>(code: string, message: string, repair: string): TracerResult<T> {
  return { ok: false, refusal: { code, message, rankedRepairs: [repair] } };
}

/** Velocity sampler bound to ONE immutable snapshot. */
export type SnapshotSampler = (
  p: readonly [number, number, number],
) => readonly [number, number, number] | null;

/** The ordered snapshot pair a tracer step interpolates across. */
export interface SnapshotPair {
  readonly tickA: number;
  readonly tickB: number;
  readonly idA: string;
  readonly idB: string;
  readonly samplerA: SnapshotSampler;
  readonly samplerB: SnapshotSampler;
}

/** Tracer budget. */
export const MAX_TRACERS = 512;
/** Retained points per tracer (the compact-retention cap). */
export const MAX_POINTS_PER_TRACER = 256;

interface TracerState {
  id: number;
  releaseTick: number;
  releasePoint: [number, number, number];
  position: [number, number, number];
  alive: boolean;
  /** retained trajectory, thinned; stride doubles when full. */
  trail: number[];
  stride: number;
  sinceKept: number;
}

/** A serializable tracer checkpoint (NO field data — ever). */
export interface TracerCheckpoint {
  readonly tick: number;
  readonly tracers: ReadonlyArray<{
    readonly id: number;
    readonly releaseTick: number;
    readonly releasePoint: readonly [number, number, number];
    readonly position: readonly [number, number, number];
    readonly alive: boolean;
    readonly trail: readonly number[];
    readonly stride: number;
    readonly sinceKept: number;
  }>;
}

/** The service. */
export class TracerService {
  private tracers: TracerState[] = [];
  private tick = 0;
  private cancelled = false;
  private nextId = 0;

  /** Tick-addressed release of one tracer. */
  release(tick: number, p: readonly [number, number, number]): TracerResult<number> {
    if (this.cancelled) {
      return refuse("tracer-cancelled", "service cancelled", "create a new service");
    }
    if (this.tracers.length >= MAX_TRACERS) {
      return refuse(
        "tracer-count-exceeded",
        `${this.tracers.length} tracers at the cap`,
        "retire tracers before releasing more",
      );
    }
    if (!p.every(Number.isFinite) || !Number.isFinite(tick) || tick < this.tick) {
      return refuse(
        "tracer-release-invalid",
        `tick ${tick} at service tick ${this.tick}`,
        "releases are tick-addressed, never retroactive",
      );
    }
    const id = this.nextId;
    this.nextId += 1;
    this.tracers.push({
      id,
      releaseTick: tick,
      releasePoint: [p[0], p[1], p[2]],
      position: [p[0], p[1], p[2]],
      alive: true,
      trail: [p[0], p[1], p[2]],
      stride: 1,
      sinceKept: 0,
    });
    return { ok: true, value: id };
  }

  /**
   * Advance one tick across the snapshot PAIR: velocity is the
   * temporal interpolation uA·(1−s) + uB·s at s = midpoint of the
   * tick interval (RK2 midpoint in space). The pair is used and
   * DROPPED — the service retains no field history by construction.
   */
  advance(pair: SnapshotPair, dtS: number): TracerResult<number> {
    if (this.cancelled) {
      return refuse("tracer-cancelled", "service cancelled", "create a new service");
    }
    if (!(dtS > 0) || !Number.isFinite(dtS)) {
      return refuse("tracer-dt-invalid", `dt ${dtS}`, "positive finite dt");
    }
    if (pair.tickB !== pair.tickA + 1 || pair.idA === "" || pair.idB === "") {
      return refuse(
        "tracer-pair-invalid",
        `pair (${pair.tickA}, ${pair.tickB})`,
        "tracers bind the ORDERED pair of adjacent snapshot ids",
      );
    }
    const sample = (
      p: readonly [number, number, number],
      s: number,
    ): readonly [number, number, number] | null => {
      const a = pair.samplerA(p);
      const b = pair.samplerB(p);
      if (a === null || b === null) return null;
      return [
        a[0] * (1 - s) + b[0] * s,
        a[1] * (1 - s) + b[1] * s,
        a[2] * (1 - s) + b[2] * s,
      ];
    };
    let advanced = 0;
    for (const t of this.tracers) {
      if (!t.alive || t.releaseTick > this.tick) continue;
      const u0 = sample(t.position, 0);
      if (u0 === null) {
        t.alive = false;
        continue;
      }
      const mid: [number, number, number] = [
        t.position[0] + 0.5 * dtS * u0[0],
        t.position[1] + 0.5 * dtS * u0[1],
        t.position[2] + 0.5 * dtS * u0[2],
      ];
      const um = sample(mid, 0.5);
      if (um === null) {
        t.alive = false;
        continue;
      }
      t.position = [
        t.position[0] + dtS * um[0],
        t.position[1] + dtS * um[1],
        t.position[2] + dtS * um[2],
      ];
      advanced += 1;
      // Compact retention: keep every stride-th point; when the
      // trail is full, thin it 2:1 and double the stride —
      // deterministic, bounded, never a dense history.
      t.sinceKept += 1;
      if (t.sinceKept >= t.stride) {
        t.sinceKept = 0;
        t.trail.push(t.position[0], t.position[1], t.position[2]);
        if (t.trail.length / 3 > MAX_POINTS_PER_TRACER) {
          const kept: number[] = [];
          for (let i = 0; i < t.trail.length / 3; i += 2) {
            kept.push(t.trail[3 * i] ?? 0, t.trail[3 * i + 1] ?? 0, t.trail[3 * i + 2] ?? 0);
          }
          t.trail = kept;
          t.stride *= 2;
        }
      }
    }
    this.tick += 1;
    return { ok: true, value: advanced };
  }

  /** Cancel: every later mutation refuses. */
  cancel(): void {
    this.cancelled = true;
  }

  /** Total retained points (the memory-cap witness). */
  retainedPoints(): number {
    return this.tracers.reduce((a, t) => a + t.trail.length / 3, 0);
  }

  /** One tracer's PATHLINE (its own retained trajectory). */
  pathline(id: number): TracerResult<Float64Array> {
    const t = this.tracers.find((x) => x.id === id);
    if (t === undefined) {
      return refuse("tracer-unknown", `id ${id}`, "release() returns the id");
    }
    return { ok: true, value: Float64Array.from(t.trail) };
  }

  /**
   * The STREAKLINE through a release point: current positions of
   * every particle released from that exact point, in release order.
   */
  streakline(releasePoint: readonly [number, number, number]): Float64Array {
    const pts: number[] = [];
    for (const t of this.tracers) {
      if (
        t.alive &&
        t.releasePoint[0] === releasePoint[0] &&
        t.releasePoint[1] === releasePoint[1] &&
        t.releasePoint[2] === releasePoint[2] &&
        t.releaseTick <= this.tick
      ) {
        pts.push(t.position[0], t.position[1], t.position[2]);
      }
    }
    return Float64Array.from(pts);
  }

  /** Serialize (particle states only — a field in here is a bug). */
  checkpoint(): TracerCheckpoint {
    return {
      tick: this.tick,
      tracers: this.tracers.map((t) => ({
        id: t.id,
        releaseTick: t.releaseTick,
        releasePoint: [...t.releasePoint] as const,
        position: [...t.position] as const,
        alive: t.alive,
        trail: [...t.trail],
        stride: t.stride,
        sinceKept: t.sinceKept,
      })),
    };
  }

  /** Reconstruct from a checkpoint (deterministic). */
  static restore(cp: TracerCheckpoint): TracerService {
    const s = new TracerService();
    s.tick = cp.tick;
    s.tracers = cp.tracers.map((t) => ({
      id: t.id,
      releaseTick: t.releaseTick,
      releasePoint: [...t.releasePoint] as [number, number, number],
      position: [...t.position] as [number, number, number],
      alive: t.alive,
      trail: [...t.trail],
      stride: t.stride,
      sinceKept: t.sinceKept,
    }));
    s.nextId = cp.tracers.reduce((a, t) => Math.max(a, t.id + 1), 0);
    return s;
  }
}
