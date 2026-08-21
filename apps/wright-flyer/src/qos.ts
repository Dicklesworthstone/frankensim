// E5.6 runtime QoS governor (bead wf-root-guzez.6.9): PURE hysteretic
// Normal → Constrained → Critical machine over render frame times.
// LAW: the governor owns PRESENTATION quality only — profiles are
// validated against an allowlist and a profile that names a physics
// knob REFUSES (tier immutability is a gate, not a convention). The
// sim tick, solver schedules, and every digest-bearing quantity are
// untouchable from here.

export type QosState = "normal" | "constrained" | "critical";

/** Atomic presentation profile (applied as a unit, never piecemeal). */
export interface PresentationProfile {
  readonly pixelRatioCap: number;
  readonly terrainDetail: "full" | "reduced";
  readonly ghostVisible: boolean;
  /** Field/sweep renderer throttle level (0 = none, 2 = heaviest). */
  readonly fieldThrottle: 0 | 1 | 2;
  /** The honest badge (null in Normal). */
  readonly badge: string | null;
}

const BADGE = "visual analysis reduced; physics unchanged";

export const PROFILES: Readonly<Record<QosState, PresentationProfile>> = Object.freeze({
  normal: Object.freeze({
    pixelRatioCap: 2,
    terrainDetail: "full" as const,
    ghostVisible: true,
    fieldThrottle: 0 as const,
    badge: null,
  }),
  constrained: Object.freeze({
    pixelRatioCap: 1.5,
    terrainDetail: "full" as const,
    ghostVisible: true,
    fieldThrottle: 1 as const,
    badge: BADGE,
  }),
  critical: Object.freeze({
    pixelRatioCap: 1,
    terrainDetail: "reduced" as const,
    ghostVisible: false,
    fieldThrottle: 2 as const,
    badge: BADGE,
  }),
});

/** The ONLY keys a presentation profile may carry. */
const ALLOWED_KEYS = new Set([
  "pixelRatioCap",
  "terrainDetail",
  "ghostVisible",
  "fieldThrottle",
  "badge",
]);

/** Tier-immutability gate: any unknown key — a physics knob, a solver
 * setting, anything — refuses. Executed by the battery's hostile twin. */
export function validatePresentationProfile(profile: Record<string, unknown>): void {
  for (const key of Object.keys(profile)) {
    if (!ALLOWED_KEYS.has(key)) {
      throw new RangeError(
        `presentation profile may not carry '${key}' — physics-tier immutability (plan law)`,
      );
    }
  }
}

export interface QosSpec {
  /** Escalate Normal→Constrained above this frame time [ms]. */
  readonly enterConstrainedMs: number;
  /** De-escalate Constrained→Normal below this [ms] (must sit BELOW
   * the enter threshold — that gap is the hysteresis). */
  readonly exitConstrainedMs: number;
  readonly enterCriticalMs: number;
  readonly exitCriticalMs: number;
  /** Consecutive qualifying frames before any transition (chatter
   * resistance). */
  readonly dwellFrames: number;
  /** Persistent-Critical budget: after this many consecutive Critical
   * frames the governor emits ONE typed performance refusal. */
  readonly refusalAfterCriticalFrames: number;
}

export const QOS_V1: QosSpec = {
  enterConstrainedMs: 22, // ~45 fps
  exitConstrainedMs: 15, // ~66 fps
  enterCriticalMs: 33, // ~30 fps
  exitCriticalMs: 22,
  dwellFrames: 45, // ~0.75 s at 60 fps
  refusalAfterCriticalFrames: 1800, // ~30 s persistently critical
};

export interface QosRefusal {
  readonly code: "performance-budget-missed";
  readonly message: string;
  readonly ranked_repairs: readonly string[];
}

export interface QosSample {
  readonly state: QosState;
  /** True exactly when the state changed on this sample. */
  readonly changed: boolean;
  readonly profile: PresentationProfile;
  /** Present exactly once, on the persistent-Critical budget miss. */
  readonly refusal?: QosRefusal;
}

export class QosGovernor {
  private readonly spec: QosSpec;
  private state: QosState = "normal";
  private escalateStreak = 0;
  private deescalateStreak = 0;
  private criticalStreak = 0;
  private refused = false;

  constructor(spec: QosSpec = QOS_V1) {
    if (
      !(spec.exitConstrainedMs < spec.enterConstrainedMs) ||
      !(spec.exitCriticalMs < spec.enterCriticalMs) ||
      spec.dwellFrames < 1 ||
      spec.refusalAfterCriticalFrames < 1
    ) {
      throw new RangeError("QoS spec: exit thresholds must sit below enter thresholds");
    }
    this.spec = spec;
  }

  current(): QosState {
    return this.state;
  }

  /** Feed one frame time; returns the (possibly new) state atomically. */
  sample(frameMs: number): QosSample {
    if (!Number.isFinite(frameMs) || frameMs < 0) {
      throw new RangeError(`frame time out of domain: ${frameMs}`);
    }
    const s = this.spec;
    let changed = false;
    // Escalation pressure (one level at a time; dwell-gated).
    const escalateThreshold =
      this.state === "normal" ? s.enterConstrainedMs : s.enterCriticalMs;
    const deescalateThreshold =
      this.state === "critical" ? s.exitCriticalMs : s.exitConstrainedMs;
    if (this.state !== "critical" && frameMs > escalateThreshold) {
      this.escalateStreak += 1;
      this.deescalateStreak = 0;
      if (this.escalateStreak >= s.dwellFrames) {
        this.state = this.state === "normal" ? "constrained" : "critical";
        this.escalateStreak = 0;
        changed = true;
      }
    } else if (this.state !== "normal" && frameMs < deescalateThreshold) {
      this.deescalateStreak += 1;
      this.escalateStreak = 0;
      if (this.deescalateStreak >= s.dwellFrames) {
        this.state = this.state === "critical" ? "constrained" : "normal";
        this.deescalateStreak = 0;
        changed = true;
      }
    } else {
      // Between thresholds: hysteresis holds the state, streaks decay.
      this.escalateStreak = 0;
      this.deescalateStreak = 0;
    }
    let refusal: QosRefusal | undefined;
    if (this.state === "critical") {
      this.criticalStreak += 1;
      if (this.criticalStreak >= s.refusalAfterCriticalFrames && !this.refused) {
        this.refused = true;
        refusal = {
          code: "performance-budget-missed",
          message: `critical presentation quality for ${this.criticalStreak} consecutive frames`,
          ranked_repairs: [
            "close other tabs / plug in the device (the sim tick is UNCHANGED)",
            "lower the window size",
            "the physics tier is immutable — no silent degradation happened",
          ],
        };
      }
    } else {
      this.criticalStreak = 0;
      this.refused = false;
    }
    const profile = PROFILES[this.state];
    validatePresentationProfile(profile as unknown as Record<string, unknown>);
    return refusal !== undefined
      ? { state: this.state, changed, profile, refusal }
      : { state: this.state, changed, profile };
  }
}
