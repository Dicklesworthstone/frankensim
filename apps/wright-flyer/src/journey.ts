// Guided journey orchestration (plan §2.1 five-minute arc; the E9.2
// onboarding scope narrowed to its URL-addressable spine). PURE stage
// math: a stage number maps to the URL the app already honors plus the
// honest caption that stage owes the player — no third config path,
// no outcome promises in any copy (plan Round-2 rule).
//
//   Stage 1 — watch a MODELED pilot fly the historical hypothesis.
//   Stage 2 — take the controls with the bounded Training Assist.
//   Stage 3 — authentic controls: raw mechanical instability.
//
// The journey advances by navigation: each terminal screen links the
// next stage's URL, so every stage stays bookmarkable and replayable.

export interface JourneyStage {
  readonly index: 1 | 2 | 3;
  /** The URL that runs this stage (existing params only). */
  readonly url: string;
  /** Caption shown while the stage runs (copy law compliant). */
  readonly caption: string;
  /** Prompt shown when the stage's run ends. */
  readonly prompt: string;
}

export const JOURNEY_STAGES: readonly JourneyStage[] = [
  {
    index: 1,
    url: "?sim=1&mode=historical&journey=1",
    caption:
      "GUIDED JOURNEY 1/3 — WATCH A MODELED PILOT: a computed hypothesis of December 17, drawn wind included. Not footage.",
    prompt: "Now fly it yourself with training assist.",
  },
  {
    index: 2,
    url: "?sim=1&mode=human&assist=1&journey=2",
    caption:
      "GUIDED JOURNEY 2/3 — TRAINING ASSIST: a bounded aid (30% canard authority). Drag = hip cradle · WASD/arrows · gamepad stick.",
    prompt: "Ready for the real thing?",
  },
  {
    index: 3,
    url: "?sim=1&mode=human&journey=3",
    caption:
      "GUIDED JOURNEY 3/3 — AUTHENTIC CONTROLS: the machine as built — unstable in pitch by design. Expect to work.",
    prompt: "The journey ends here; the engineer's views remain yours.",
  },
];

/** Parse a `journey` URL param into its stage, or null. */
export function journeyStage(param: string | null): JourneyStage | null {
  if (param === null) {
    return null;
  }
  const n = Number(param);
  if (!Number.isInteger(n)) {
    return null;
  }
  return JOURNEY_STAGES.find((s) => s.index === n) ?? null;
}

/** The URL of the NEXT stage after `n`, or null when the journey ends. */
export function journeyNextUrl(n: number): string | null {
  return journeyStage(String(n + 1))?.url ?? null;
}
