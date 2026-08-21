// Field-visualization logic (bead wf-root-guzez.8.3, E7.2, plan
// §5.5 consumers): glyph instancing, streamline integration,
// divergence overlay duals, wake age-fade, probe gizmos, legends.
// Pure logic — three.js consumes the typed arrays this module emits,
// so every law here is testable under node --test.
//
// Doctrine carried from the field service:
//   - the divergence overlay shows BOTH absolute |div u| and the
//     normalized value, masking normalized where the gradient norm is
//     under the floor or inside singularity cores (the teaching
//     toggle switches between the analytic and FD duals — same mask);
//   - legends are ALWAYS visible and never label a sum "total flow"
//     while a force-coupled component is omitted (the omissions are
//     named instead);
//   - wake filaments with age-fade are PRESENTATION-ONLY and say so.

export interface VizRefusal {
  readonly code: string;
  readonly message: string;
  readonly rankedRepairs: readonly string[];
}

export type VizResult<T> = { ok: true; value: T } | { ok: false; refusal: VizRefusal };

function refuse<T>(code: string, message: string, repair: string): VizResult<T> {
  return { ok: false, refusal: { code, message, rankedRepairs: [repair] } };
}

/** The §5.5 arrays a snapshot-lease sample delivers to the renderer. */
export interface FieldArrays {
  readonly n: number;
  /** xyz per point. */
  readonly points: Float64Array;
  /** velocity xyz per point [m/s]. */
  readonly u: Float64Array;
  readonly divAnalytic: Float64Array;
  readonly divFd: Float64Array;
  /** Frobenius norm of grad_u per point [1/s]. */
  readonly gradNorm: Float64Array;
  readonly validity: Uint8Array;
  readonly singularityCore: Uint8Array;
  /** Names of components absent from the sum (meta passthrough). */
  readonly omittedComponents: readonly string[];
  /** Force-coupled component names the source model supports. */
  readonly forceCoupledSupported: readonly string[];
}

/** Instanced-glyph budget (plan: <= 30k). */
export const MAX_GLYPHS = 30_000;

/** Mirror of the native gradient-norm floor. */
export const GRAD_NORM_FLOOR = 1e-6;

export interface GlyphInstances {
  readonly count: number;
  /** xyz per glyph. */
  readonly positions: Float64Array;
  /** unit direction xyz per glyph. */
  readonly directions: Float64Array;
  /** |u| per glyph [m/s]. */
  readonly magnitudes: Float64Array;
}

/**
 * Build glyph instances from VALID, non-core points only. Refuses
 * TYPED beyond the instancing budget (AT the cap admits) — silent
 * downsampling would misrepresent the field.
 */
export function buildGlyphInstances(field: FieldArrays): VizResult<GlyphInstances> {
  const idx: number[] = [];
  for (let i = 0; i < field.n; i += 1) {
    if (field.validity[i] === 1 && field.singularityCore[i] === 0) idx.push(i);
  }
  if (idx.length > MAX_GLYPHS) {
    return refuse(
      "glyph-count-exceeded",
      `${idx.length} candidate glyphs > ${MAX_GLYPHS}`,
      "sample a coarser grid; silent downsampling misrepresents the field",
    );
  }
  const positions = new Float64Array(idx.length * 3);
  const directions = new Float64Array(idx.length * 3);
  const magnitudes = new Float64Array(idx.length);
  idx.forEach((p, k) => {
    const ux = field.u[3 * p] ?? 0;
    const uy = field.u[3 * p + 1] ?? 0;
    const uz = field.u[3 * p + 2] ?? 0;
    const mag = Math.sqrt(ux * ux + uy * uy + uz * uz);
    magnitudes[k] = mag;
    positions[3 * k] = field.points[3 * p] ?? 0;
    positions[3 * k + 1] = field.points[3 * p + 1] ?? 0;
    positions[3 * k + 2] = field.points[3 * p + 2] ?? 0;
    if (mag > 1e-12) {
      directions[3 * k] = ux / mag;
      directions[3 * k + 1] = uy / mag;
      directions[3 * k + 2] = uz / mag;
    }
  });
  return { ok: true, value: { count: idx.length, positions, directions, magnitudes } };
}

export type DivergenceDual = "analytic" | "finite-difference";

export interface DivergenceOverlay {
  readonly dual: DivergenceDual;
  /** |div u| per point (always shown). */
  readonly absolute: Float64Array;
  /** normalized eps per point; NaN where masked (floor/core/invalid). */
  readonly normalized: Float64Array;
  /** true where normalized is masked. */
  readonly masked: Uint8Array;
}

/**
 * The divergence overlay: BOTH absolute and normalized, from the
 * SELECTED dual (the teaching toggle). The mask law is identical for
 * both duals — switching duals never unmasks a point.
 */
export function divergenceOverlay(field: FieldArrays, dual: DivergenceDual): DivergenceOverlay {
  const src = dual === "analytic" ? field.divAnalytic : field.divFd;
  const absolute = new Float64Array(field.n);
  const normalized = new Float64Array(field.n);
  const masked = new Uint8Array(field.n);
  for (let i = 0; i < field.n; i += 1) {
    absolute[i] = Math.abs(src[i] ?? 0);
    const mask =
      field.validity[i] !== 1 ||
      field.singularityCore[i] === 1 ||
      (field.gradNorm[i] ?? 0) < GRAD_NORM_FLOOR;
    if (mask) {
      masked[i] = 1;
      normalized[i] = Number.NaN;
    } else {
      normalized[i] = (absolute[i] ?? 0) / (field.gradNorm[i] ?? 1);
    }
  }
  return { dual, absolute, normalized, masked };
}

/** Streamline seed/step budgets. */
export const MAX_SEEDS = 256;
export const MAX_STEPS = 4_096;

export interface Streamline {
  /** xyz polyline (3 * nPoints). */
  readonly points: Float64Array;
  /** why integration stopped. */
  readonly ended: "steps-exhausted" | "left-domain";
}

/**
 * Fixed-step RK4 streamline integration over a sampler bound to the
 * LATEST snapshot (the sampler returns null outside its domain).
 * Deterministic: fixed step, fixed order, no adaptivity.
 */
export function integrateStreamlines(
  sampler: (p: readonly [number, number, number]) => readonly [number, number, number] | null,
  seeds: ReadonlyArray<readonly [number, number, number]>,
  stepS: number,
  nSteps: number,
): VizResult<Streamline[]> {
  if (!(Number.isFinite(stepS) && stepS > 0)) {
    return refuse("streamline-step-invalid", `step ${stepS}`, "a finite positive step");
  }
  if (seeds.length === 0 || seeds.length > MAX_SEEDS) {
    return refuse(
      "streamline-seeds-invalid",
      `${seeds.length} seeds outside [1, ${MAX_SEEDS}]`,
      "seed from the probe gizmos or the glyph grid",
    );
  }
  if (nSteps < 1 || nSteps > MAX_STEPS) {
    return refuse(
      "streamline-steps-invalid",
      `${nSteps} steps outside [1, ${MAX_STEPS}]`,
      "budget the integration to the frame",
    );
  }
  const lines: Streamline[] = [];
  for (const seed of seeds) {
    const pts: number[] = [seed[0], seed[1], seed[2]];
    let p: [number, number, number] = [seed[0], seed[1], seed[2]];
    let ended: Streamline["ended"] = "steps-exhausted";
    for (let s = 0; s < nSteps; s += 1) {
      const k1 = sampler(p);
      if (k1 === null) {
        ended = "left-domain";
        break;
      }
      const at = (
        q: [number, number, number],
      ): readonly [number, number, number] | null => sampler(q);
      const mid1: [number, number, number] = [
        p[0] + 0.5 * stepS * k1[0],
        p[1] + 0.5 * stepS * k1[1],
        p[2] + 0.5 * stepS * k1[2],
      ];
      const k2 = at(mid1);
      if (k2 === null) {
        ended = "left-domain";
        break;
      }
      const mid2: [number, number, number] = [
        p[0] + 0.5 * stepS * k2[0],
        p[1] + 0.5 * stepS * k2[1],
        p[2] + 0.5 * stepS * k2[2],
      ];
      const k3 = at(mid2);
      if (k3 === null) {
        ended = "left-domain";
        break;
      }
      const end: [number, number, number] = [
        p[0] + stepS * k3[0],
        p[1] + stepS * k3[1],
        p[2] + stepS * k3[2],
      ];
      const k4 = at(end);
      if (k4 === null) {
        ended = "left-domain";
        break;
      }
      p = [
        p[0] + (stepS / 6) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]),
        p[1] + (stepS / 6) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]),
        p[2] + (stepS / 6) * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]),
      ];
      pts.push(p[0], p[1], p[2]);
    }
    lines.push({ points: Float64Array.from(pts), ended });
  }
  return { ok: true, value: lines };
}

/** Presentation-only marker for the wake-filament rendering path. */
export const WAKE_FADE_PRESENTATION_ONLY =
  "wake filament age-fade is presentation-only; it never alters the physical replay";

/**
 * Age-fade opacity for a wake filament row: 1 at age 0, monotone
 * down to 0 at maxAge, clamped outside.
 */
export function wakeAgeFade(ageTicks: number, maxAgeTicks: number): number {
  if (!(maxAgeTicks > 0) || !Number.isFinite(ageTicks)) return 0;
  const t = ageTicks / maxAgeTicks;
  if (t <= 0) return 1;
  if (t >= 1) return 0;
  // Smoothstep-down: continuous, monotone.
  return 1 - t * t * (3 - 2 * t);
}

export interface LegendConfig {
  readonly label: string;
  readonly units: string;
  readonly transfer: string;
  readonly alwaysVisible: true;
  /** Named omissions the legend must display alongside the label. */
  readonly omittedComponents: readonly string[];
}

/**
 * Legend for a rendered quantity. The HONESTY law lives here: the
 * velocity-sum legend is titled "total flow" ONLY when no supported
 * force-coupled component is omitted; otherwise it is "selected
 * components" with the omissions named.
 */
export function legendConfig(
  quantity: "velocity" | "vorticity" | "divergence-absolute" | "divergence-normalized",
  field: FieldArrays,
): LegendConfig {
  const omittedForceCoupled = field.omittedComponents.filter((c) =>
    field.forceCoupledSupported.includes(c),
  );
  const base = {
    transfer: "viridis-class",
    alwaysVisible: true as const,
    omittedComponents: field.omittedComponents,
  };
  switch (quantity) {
    case "velocity":
      return {
        ...base,
        label: omittedForceCoupled.length === 0 ? "total flow" : "selected components",
        units: "m/s",
      };
    case "vorticity":
      return { ...base, label: "vorticity", units: "1/s" };
    case "divergence-absolute":
      return { ...base, label: "|div u| (absolute)", units: "1/s" };
    case "divergence-normalized":
      return { ...base, label: "eps_div (normalized, masked)", units: "1" };
    default: {
      const never: never = quantity;
      throw new Error(`unreachable legend quantity ${String(never)}`);
    }
  }
}

/** A probe gizmo: position + the sampled row it points at. */
export interface ProbeGizmo {
  readonly position: readonly [number, number, number];
  readonly pointIndex: number;
  readonly valid: boolean;
}

/** Probe cap (HUD budget). */
export const MAX_PROBES = 16;

/**
 * Bind probe gizmos to their nearest field points (deterministic
 * first-minimum tie-break by index).
 */
export function bindProbes(
  field: FieldArrays,
  probes: ReadonlyArray<readonly [number, number, number]>,
): VizResult<ProbeGizmo[]> {
  if (probes.length > MAX_PROBES) {
    return refuse(
      "probe-count-exceeded",
      `${probes.length} probes > ${MAX_PROBES}`,
      "the HUD shows at most 16 probes",
    );
  }
  const out: ProbeGizmo[] = [];
  for (const probe of probes) {
    let best = 0;
    let bestD = Number.POSITIVE_INFINITY;
    for (let i = 0; i < field.n; i += 1) {
      const dx = (field.points[3 * i] ?? 0) - probe[0];
      const dy = (field.points[3 * i + 1] ?? 0) - probe[1];
      const dz = (field.points[3 * i + 2] ?? 0) - probe[2];
      const d = dx * dx + dy * dy + dz * dz;
      if (d < bestD) {
        bestD = d;
        best = i;
      }
    }
    out.push({ position: probe, pointIndex: best, valid: field.validity[best] === 1 });
  }
  return { ok: true, value: out };
}
