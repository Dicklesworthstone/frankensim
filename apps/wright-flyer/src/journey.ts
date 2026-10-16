// Guided journey orchestration (plan §2.1 five-minute arc; the E9.2
// onboarding scope narrowed to its URL-addressable spine). PURE stage
// math: a stage number maps to the URL the app already honors plus the
// honest caption that stage owes the player — no third config path,
// no outcome promises in any copy (plan Round-2 rule).
//
//   Stage 1 — watch a MODELED pilot fly the historical hypothesis.
//   Stage 2 — take the controls with the bounded Training Assist.
//   Stage 3 — authentic controls: raw mechanical instability.
//   Stage 4 — open the existing "Why it porpoises" instrument.
//
// The journey advances by navigation: each terminal screen links the
// next stage's URL, so every stage stays bookmarkable and replayable.

export interface JourneyStage {
  readonly index: 1 | 2 | 3 | 4;
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
      "GUIDED JOURNEY 1/4 — WATCH A MODELED PILOT: a computed hypothesis of December 17, drawn wind included. Not footage.",
    prompt: "Now fly it yourself with training assist.",
  },
  {
    index: 2,
    url: "?sim=1&mode=human&assist=1&journey=2",
    caption:
      "GUIDED JOURNEY 2/4 — TRAINING ASSIST: a bounded aid (30% canard authority). Drag = hip cradle · WASD/arrows · gamepad stick.",
    prompt: "Ready for the real thing?",
  },
  {
    index: 3,
    url: "?sim=1&mode=human&journey=3",
    caption:
      "GUIDED JOURNEY 3/4 — AUTHENTIC CONTROLS: the machine as built — unstable in pitch by design. Expect to work.",
    prompt: "Now inspect why it porpoises.",
  },
  {
    index: 4,
    url: "?sim=1&mode=human&journey=4&inst=porpoises",
    caption:
      "GUIDED JOURNEY 4/4 — WHY IT PORPOISES: an illustrative live view of command, response, lag, and growth cues — not released engine eigenmode authority.",
    prompt: "Journey complete. Keep exploring the engineer's views.",
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
