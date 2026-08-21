// Lesson scaffolding + pilot-perception view (bead
// wf-root-guzez.8.7, E7.5). Curated overlay scripts — each lesson is
// a deterministic sequence of steps that name the overlays they
// require and the claims they may voice — plus the perception view,
// which shows the pilot model's CUE STATE beside the true state.
//
// Laws:
//   - a lesson step may only voice claims from its declared claim
//     list (curated scripts never improvise physics claims);
//   - the perception view consumes the DETERMINISTIC cue state the
//     sim plane publishes — never renderer output (structural: the
//     API takes cue/true rows, there is no render-tap input);
//   - the perceived-vs-true gap is DISPLAYED per cue, never hidden.

export interface LessonRefusal {
  readonly code: string;
  readonly message: string;
  readonly rankedRepairs: readonly string[];
}

export type LessonResult<T> = { ok: true; value: T } | { ok: false; refusal: LessonRefusal };

function refuse<T>(code: string, message: string, repair: string): LessonResult<T> {
  return { ok: false, refusal: { code, message, rankedRepairs: [repair] } };
}

/** Overlay ids a lesson step may require (the E7 render stack). */
export type OverlayId =
  | "glyphs"
  | "streamlines"
  | "divergence"
  | "force-overlay"
  | "eigenmodes"
  | "ground-images"
  | "perception";

export interface LessonStep {
  readonly title: string;
  readonly overlays: readonly OverlayId[];
  /** claims this step voices — every one must be in the lesson's declared list. */
  readonly voicedClaims: readonly string[];
}

export interface Lesson {
  readonly id: string;
  readonly title: string;
  /** the curated claim list (with provenance strings). */
  readonly declaredClaims: readonly string[];
  readonly steps: readonly LessonStep[];
}

/** Step budget per lesson. */
export const MAX_STEPS = 32;

/** Validate a lesson: steps bounded, every voiced claim declared. */
export function validateLesson(lesson: Lesson): LessonResult<Lesson> {
  if (lesson.steps.length === 0 || lesson.steps.length > MAX_STEPS) {
    return refuse(
      "lesson-steps-invalid",
      `${lesson.steps.length} steps outside [1, ${MAX_STEPS}]`,
      "lessons are short and curated",
    );
  }
  for (const [i, step] of lesson.steps.entries()) {
    for (const claim of step.voicedClaims) {
      if (!lesson.declaredClaims.includes(claim)) {
        return refuse(
          "lesson-undeclared-claim",
          `step ${i} ('${step.title}') voices an undeclared claim: '${claim}'`,
          "curated scripts never improvise physics claims; declare it with provenance",
        );
      }
    }
    if (step.overlays.length === 0) {
      return refuse(
        "lesson-step-blind",
        `step ${i} shows nothing`,
        "every step names the overlays it teaches with",
      );
    }
  }
  return { ok: true, value: lesson };
}

/** Run state: a validated lesson plus a cursor (deterministic). */
export interface LessonRun {
  readonly lesson: Lesson;
  readonly stepIndex: number;
  readonly done: boolean;
}

/** Start a lesson (validates first). */
export function startLesson(lesson: Lesson): LessonResult<LessonRun> {
  const v = validateLesson(lesson);
  if (!v.ok) return v;
  return { ok: true, value: { lesson: v.value, stepIndex: 0, done: false } };
}

/** Advance to the next step; sets done past the last. */
export function advanceLesson(run: LessonRun): LessonResult<LessonRun> {
  if (run.done) {
    return refuse("lesson-already-done", run.lesson.id, "restart to run again");
  }
  const next = run.stepIndex + 1;
  return {
    ok: true,
    value: {
      lesson: run.lesson,
      stepIndex: Math.min(next, run.lesson.steps.length - 1),
      done: next >= run.lesson.steps.length,
    },
  };
}

/** The three curated v1 lessons (runnable end-to-end). */
export function curatedLessons(): Lesson[] {
  const groundClaims = [
    "ground effect raises lift at fixed alpha near the certified flat (V-06a receipt)",
    "the image system is a boundary device, not physical vorticity",
    "the 1903 flights flew IN ground effect for their whole length",
  ];
  const anhedralClaims = [
    "the 1903 anhedral was intentional: it damps the upset a gust starts (historical flag)",
    "dihedral stabilizes spirally but couples gusts into roll",
  ];
  const stabilityClaims = [
    "fixed-stick poles differ from free-control poles (V-02b sign/tendency evidence class)",
    "the canard's free-control branch is the porpoising ingredient",
    "time-to-double is ln 2 over the unstable pole's real part",
  ];
  return [
    {
      id: "lesson-ground-effect",
      title: "Ground effect at Kill Devil Hills",
      declaredClaims: groundClaims,
      steps: [
        { title: "Free-air baseline", overlays: ["glyphs", "force-overlay"], voicedClaims: [] },
        {
          title: "Descend toward the certified flat",
          overlays: ["ground-images", "force-overlay"],
          voicedClaims: [groundClaims[0] ?? ""],
        },
        {
          title: "The mirror airplane",
          overlays: ["ground-images", "streamlines"],
          voicedClaims: [groundClaims[1] ?? "", groundClaims[2] ?? ""],
        },
      ],
    },
    {
      id: "lesson-anhedral",
      title: "Why the wings droop",
      declaredClaims: anhedralClaims,
      steps: [
        { title: "A gust hits", overlays: ["glyphs", "perception"], voicedClaims: [] },
        {
          title: "Anhedral vs dihedral",
          overlays: ["force-overlay", "eigenmodes"],
          voicedClaims: [...anhedralClaims],
        },
      ],
    },
    {
      id: "lesson-fixed-vs-free",
      title: "Fixed stick, free stick",
      declaredClaims: stabilityClaims,
      steps: [
        { title: "Hold the canard", overlays: ["eigenmodes"], voicedClaims: [stabilityClaims[0] ?? ""] },
        {
          title: "Let go",
          overlays: ["eigenmodes", "perception"],
          voicedClaims: [stabilityClaims[1] ?? ""],
        },
        {
          title: "Time to double",
          overlays: ["eigenmodes"],
          voicedClaims: [stabilityClaims[2] ?? ""],
        },
      ],
    },
  ];
}

/** One cue row: the modeled percept beside the true state value. */
export interface CueRow {
  readonly cue: string;
  /** the pilot model's PUBLISHED perceived value. */
  readonly perceived: number;
  /** the true state value, same units. */
  readonly truth: number;
  readonly units: string;
}

export interface PerceptionViewRow {
  readonly cue: string;
  readonly perceived: number;
  readonly truth: number;
  /** the DISPLAYED gap (perceived − truth) — never hidden. */
  readonly gap: number;
  readonly units: string;
}

/** Cue budget. */
export const MAX_CUES = 16;

/**
 * Build the perception view rows: verbatim cue state beside truth
 * with the gap displayed. There is no renderer input to this
 * function BY CONSTRUCTION — it types only the published cue rows.
 */
export function perceptionView(rows: readonly CueRow[]): LessonResult<PerceptionViewRow[]> {
  if (rows.length === 0 || rows.length > MAX_CUES) {
    return refuse(
      "perception-cues-invalid",
      `${rows.length} cues outside [1, ${MAX_CUES}]`,
      "the perception model publishes a bounded cue set",
    );
  }
  const out: PerceptionViewRow[] = [];
  for (const r of rows) {
    if (!Number.isFinite(r.perceived) || !Number.isFinite(r.truth) || r.cue.trim() === "") {
      return refuse(
        "perception-cue-invalid",
        `cue '${r.cue}'`,
        "finite perceived/truth, named cue",
      );
    }
    out.push({
      cue: r.cue,
      perceived: r.perceived,
      truth: r.truth,
      gap: r.perceived - r.truth,
      units: r.units,
    });
  }
  return { ok: true, value: out };
}
