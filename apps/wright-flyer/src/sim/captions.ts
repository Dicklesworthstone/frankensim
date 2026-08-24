// E5.3b-ii honest hypothesis captioning (bead wf-root-guzez.6.5.2):
// PURE caption stream for the historical ride-along (§2.1). Every
// caption carries an EVIDENCE LABEL — the formatter REFUSES an
// unlabeled claim (the evidence UX is a gate, not a style). Event
// detection runs on the live snapshot stream: rail run, liftoff,
// undulation crests/troughs, terminal. The undulation caption links
// toward the WHY-IT-PORPOISES view (E7.4 anchor, stubbed until it
// lands — the link target is declared, never faked).

import type { SimSnapshot } from "./snapshotView.ts";

/** The evidence taxonomy (plan evidence UX). */
export type EvidenceLabel = "Verified" | "Estimated" | "Hypothesis";

export interface Caption {
  readonly atTick: number;
  readonly label: EvidenceLabel;
  readonly text: string;
  /** Optional view anchor (e.g. the porpoises view). */
  readonly link?: "why-it-porpoises";
}

/** Format one caption line. Throws on an empty text or a text that
 * already smuggles a bracket tag (double-labeling hides the gate). */
export function formatCaption(c: Caption): string {
  if (c.text.length === 0) {
    throw new RangeError("caption text must be non-empty");
  }
  if (/^\[/.test(c.text)) {
    throw new RangeError("caption text must not carry its own label bracket");
  }
  const link = c.link !== undefined ? ` → ${c.link}` : "";
  return `[${c.label}] ${c.text}${link}`;
}

/** Streaming event detector: feed every snapshot in tick order. */
export class CaptionStream {
  private readonly captions: Caption[] = [];
  private started = false;
  private airborne = false;
  private lastQSign = 0;
  private undulationCount = 0;
  private endedTick: number | null = null;

  /** All captions so far (append-only). */
  all(): readonly Caption[] {
    return this.captions;
  }

  /** Captions at or before `tick` (the HUD shows the latest few). */
  upTo(tick: number): readonly Caption[] {
    return this.captions.filter((c) => c.atTick <= tick);
  }

  feed(s: SimSnapshot): void {
    if (this.endedTick !== null) {
      return;
    }
    if (!this.started) {
      this.started = true;
      this.captions.push({
        atTick: s.tick,
        label: "Verified",
        text: "Rail run: 60 ft launch rail, headwind start (Dec 17 procedure).",
      });
    }
    if (!this.airborne && s.phase === "airborne") {
      this.airborne = true;
      this.captions.push({
        atTick: s.tick,
        label: "Estimated",
        text: `Liftoff at ${(s.tick / 120).toFixed(1)} s — lift exceeds weight under the coupled build-up.`,
      });
    }
    if (this.airborne && !s.ended) {
      const sign = s.qRadS > 1e-3 ? 1 : s.qRadS < -1e-3 ? -1 : 0;
      if (sign !== 0 && this.lastQSign !== 0 && sign !== this.lastQSign) {
        this.undulationCount += 1;
        if (this.undulationCount % 2 === 0) {
          this.captions.push({
            atTick: s.tick,
            label: "Hypothesis",
            text: `Undulation ${this.undulationCount / 2}: the unstable canard pitch mode, barely held by the pilot model.`,
            link: "why-it-porpoises",
          });
        }
      }
      if (sign !== 0) {
        this.lastQSign = sign;
      }
    }
    if (s.ended && this.endedTick === null) {
      this.endedTick = s.tick;
      const text =
        s.phase === "ended:ground-contact"
          ? `Ground contact at ${(s.tick / 120).toFixed(1)} s, ${s.xM.toFixed(0)} m — the historical ending class.`
          : s.phase === "ended:rail-end-without-lift"
            ? "Ran off the rail without lift — several December attempts ended this way."
            : s.phase === "ended:envelope-exceeded"
              ? `Final plunge left the certified aero envelope at ${(s.tick / 120).toFixed(1)} s — the run ends with a receipt, not a guess.`
              : s.phase === "ended:damage-model-unavailable"
                ? "A swept-feature strike ended the physical run — cinematic continuation is presentation only."
              : "Tick budget reached — still flying.";
      this.captions.push({
        atTick: s.tick,
        label: s.phase === "ended:envelope-exceeded" ? "Verified" : "Estimated",
        text,
      });
    }
  }

  undulations(): number {
    return Math.floor(this.undulationCount / 2);
  }
}
