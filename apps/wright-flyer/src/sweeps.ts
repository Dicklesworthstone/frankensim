// Worker-pool sweep engine (bead wf-root-guzez.9.1, E8.1). 1-D/2-D
// config grids, common-random-number ensembles, progress streaming,
// RunSpecId-keyed caching, CSV export — QoS-SUBORDINATED: a sweep
// step dispatches ONLY when the QoS gate reports headroom, so sweeps
// can never increase sim deadline misses (the V-14 clause is a
// structural property of the dispatcher, verified by the battery
// with a closed gate).
//
// CRN law: ensemble member k uses the SAME realization seed at EVERY
// design point — differences between design points are then design
// effects, never realization luck.

export interface SweepRefusal {
  readonly code: string;
  readonly message: string;
  readonly rankedRepairs: readonly string[];
}

export type SweepResult<T> = { ok: true; value: T } | { ok: false; refusal: SweepRefusal };

function refuse<T>(code: string, message: string, repair: string): SweepResult<T> {
  return { ok: false, refusal: { code, message, rankedRepairs: [repair] } };
}

export interface SweepAxis {
  readonly name: string;
  readonly values: readonly number[];
}

/** Per-axis point cap. */
export const MAX_AXIS_POINTS = 64;
/** Total design-point cap. */
export const MAX_DESIGN_POINTS = 1_024;
/** Ensemble cap. */
export const MAX_ENSEMBLE = 32;

export interface DesignPoint {
  readonly index: number;
  readonly config: Readonly<Record<string, number>>;
}

/** Build the (1-D or 2-D) design grid, row-major, deterministic. */
export function makeSweepGrid(
  axis1: SweepAxis,
  axis2?: SweepAxis,
): SweepResult<DesignPoint[]> {
  for (const axis of axis2 === undefined ? [axis1] : [axis1, axis2]) {
    if (axis.values.length === 0 || axis.values.length > MAX_AXIS_POINTS) {
      return refuse(
        "sweep-axis-invalid",
        `${axis.name}: ${axis.values.length} points outside [1, ${MAX_AXIS_POINTS}]`,
        "sweeps are bounded studies",
      );
    }
    if (axis.values.some((v) => !Number.isFinite(v)) || axis.name.trim() === "") {
      return refuse("sweep-axis-invalid", `${axis.name}: non-finite or unnamed`, "finite named axes");
    }
  }
  const total = axis1.values.length * (axis2?.values.length ?? 1);
  if (total > MAX_DESIGN_POINTS) {
    return refuse(
      "sweep-grid-too-large",
      `${total} design points > ${MAX_DESIGN_POINTS}`,
      "coarsen an axis",
    );
  }
  const points: DesignPoint[] = [];
  for (const v1 of axis1.values) {
    if (axis2 === undefined) {
      points.push({ index: points.length, config: { [axis1.name]: v1 } });
    } else {
      for (const v2 of axis2.values) {
        points.push({
          index: points.length,
          config: { [axis1.name]: v1, [axis2.name]: v2 },
        });
      }
    }
  }
  return { ok: true, value: points };
}

/** The RunSpecId: exact key over point config + seed + model version. */
export function runSpecId(
  point: DesignPoint,
  realizationSeed: number,
  modelVersion: string,
): string {
  const cfg = Object.entries(point.config)
    .sort(([a], [b]) => (a < b ? -1 : 1))
    .map(([k, v]) => `${k}=${v}`)
    .join(",");
  return `${modelVersion}|${cfg}|seed=${realizationSeed}`;
}

export interface SweepRunRecord {
  readonly pointIndex: number;
  readonly member: number;
  readonly seed: number;
  readonly value: number;
  readonly fromCache: boolean;
}

export interface SweepProgress {
  readonly completed: number;
  readonly total: number;
}

/** Runner: executes one (design point, realization seed) unit. */
export type SweepRunner = (point: DesignPoint, seed: number) => number;

/** QoS gate: true = the sim loop has headroom for one sweep unit. */
export type QosGate = () => boolean;

/** The stepwise sweep engine (worker-pool dispatch seam). */
export class SweepEngine {
  private readonly points: DesignPoint[];
  private readonly seeds: number[];
  private readonly modelVersion: string;
  private readonly cache: Map<string, number>;
  private cursor = 0;
  readonly records: SweepRunRecord[] = [];

  constructor(
    points: DesignPoint[],
    ensembleSeeds: number[],
    modelVersion: string,
    cache: Map<string, number>,
  ) {
    this.points = points;
    this.seeds = ensembleSeeds;
    this.modelVersion = modelVersion;
    this.cache = cache;
  }

  /** Total work units. */
  total(): number {
    return this.points.length * this.seeds.length;
  }

  /** Progress snapshot (streamed to the UI). */
  progress(): SweepProgress {
    return { completed: this.cursor, total: this.total() };
  }

  /**
   * Dispatch ONE unit if the QoS gate has headroom. Returns whether a
   * unit ran — with the gate closed, NO work happens and NO runner is
   * called (the deadline-miss protection law).
   */
  step(gate: QosGate, runner: SweepRunner): boolean {
    if (this.cursor >= this.total() || !gate()) {
      return false;
    }
    const pointIndex = Math.floor(this.cursor / this.seeds.length);
    const member = this.cursor % this.seeds.length;
    const point = this.points[pointIndex];
    const seed = this.seeds[member];
    if (point === undefined || seed === undefined) return false;
    const key = runSpecId(point, seed, this.modelVersion);
    const cached = this.cache.get(key);
    const value = cached ?? runner(point, seed);
    if (cached === undefined) this.cache.set(key, value);
    this.records.push({
      pointIndex,
      member,
      seed,
      value,
      fromCache: cached !== undefined,
    });
    this.cursor += 1;
    return true;
  }
}

/** Build a sweep engine (validates the ensemble). */
export function makeSweepEngine(
  points: DesignPoint[],
  ensembleSeeds: readonly number[],
  modelVersion: string,
  cache: Map<string, number>,
): SweepResult<SweepEngine> {
  if (ensembleSeeds.length === 0 || ensembleSeeds.length > MAX_ENSEMBLE) {
    return refuse(
      "sweep-ensemble-invalid",
      `${ensembleSeeds.length} members outside [1, ${MAX_ENSEMBLE}]`,
      "a bounded CRN ensemble",
    );
  }
  if (new Set(ensembleSeeds).size !== ensembleSeeds.length) {
    return refuse(
      "sweep-ensemble-invalid",
      "duplicate realization seeds",
      "distinct seeds per member (CRN pairs across POINTS, not members)",
    );
  }
  if (modelVersion.trim() === "") {
    return refuse("sweep-model-version-missing", "unnamed model", "pin the model version");
  }
  return {
    ok: true,
    value: new SweepEngine(points, [...ensembleSeeds], modelVersion, cache),
  };
}

/** Deterministic CSV export (grid order; stable header). */
export function exportCsv(engine: SweepEngine, points: DesignPoint[]): string {
  const axisNames = Object.keys(points[0]?.config ?? {}).sort();
  const header = ["point", ...axisNames, "member", "seed", "value", "from_cache"].join(",");
  const rows = engine.records.map((r) => {
    const p = points[r.pointIndex];
    const cfg = axisNames.map((n) => p?.config[n] ?? Number.NaN);
    return [r.pointIndex, ...cfg, r.member, r.seed, r.value, r.fromCache].join(",");
  });
  return [header, ...rows].join("\n");
}
