// Shared General MIDI percussion mapping for score playback and live input.
// Positions are fractions of the playing surface radius; g scales physical
// stick speed, never PCM gain. Instrument substitutions retain their labels.

export type PieceId = "kick" | "snare" | "snare-off" | "tom-12" | "tom-14" | "hihat" | "crash" | "ride";
export const PIECE_NAMES: Record<PieceId, string> = {
  kick: "Kick", snare: "Snare", "snare-off": "Snare, wires off", "tom-12": "Rack tom", "tom-14": "Floor tom",
  hihat: "Hi-hat", crash: "Crash", ride: "Ride",
};

/** MIDI velocity to a factor on a piece's standard stroke: 100 → 1, 64 → 0.68, 127 → 1.24 (clamped 0.25..1.35). */
export function stickFactor(velocity: number): number {
  return Math.max(0.25, Math.min(1.35, (Math.max(1, velocity) / 100) ** 0.88));
}

export interface DrumRule { piece: PieceId; u: number; v: number; g: number; label?: string }
export const DRUM_MAP: Record<number, DrumRule> = {
  35: { piece: "kick", u: 0, v: 0.02, g: 1 },
  36: { piece: "kick", u: 0, v: 0.02, g: 1 },
  37: { piece: "snare", u: 0.86, v: 0.05, g: 0.75, label: "Side stick → snare, near the rim" },
  38: { piece: "snare", u: 0.35, v: 0.1, g: 1 },
  39: { piece: "snare", u: 0.55, v: -0.2, g: 0.8, label: "Hand clap → snare" },
  40: { piece: "snare", u: 0.3, v: 0.05, g: 1.05 },
  41: { piece: "tom-14", u: 0.35, v: 0.1, g: 1 },
  43: { piece: "tom-14", u: 0.3, v: 0.12, g: 1 },
  42: { piece: "hihat", u: 0.85, v: 0, g: 1 },
  44: { piece: "hihat", u: 0, v: 0, g: 0 },
  46: { piece: "hihat", u: 0.85, v: 0, g: 1 },
  45: { piece: "tom-12", u: 0.4, v: 0.1, g: 1 },
  47: { piece: "tom-12", u: 0.35, v: 0.1, g: 1 },
  48: { piece: "tom-12", u: 0.3, v: 0.12, g: 1 },
  50: { piece: "tom-12", u: 0.25, v: 0.1, g: 1.05 },
  49: { piece: "crash", u: 0.85, v: 0, g: 1 },
  57: { piece: "crash", u: 0.8, v: 0.1, g: 1 },
  52: { piece: "crash", u: 0.88, v: -0.1, g: 1, label: "China cymbal → crash" },
  55: { piece: "crash", u: 0.9, v: 0.05, g: 0.8, label: "Splash cymbal → crash" },
  51: { piece: "ride", u: 0.6, v: 0, g: 1 },
  59: { piece: "ride", u: 0.62, v: 0.06, g: 1 },
  53: { piece: "ride", u: 0.12, v: 0, g: 1.1, label: "Ride bell → ride, on the bell" },
  60: { piece: "tom-12", u: 0.7, v: 0, g: 0.8, label: "Bongo → rack tom, off centre" },
  61: { piece: "tom-12", u: 0.75, v: 0.1, g: 0.75, label: "Bongo → rack tom, off centre" },
  62: { piece: "tom-14", u: 0.7, v: 0, g: 0.6, label: "Conga → floor tom" },
  63: { piece: "tom-14", u: 0.55, v: 0.05, g: 0.8, label: "Conga → floor tom" },
  64: { piece: "tom-14", u: 0.5, v: 0.1, g: 0.85, label: "Conga → floor tom" },
  65: { piece: "snare-off", u: 0.4, v: 0.1, g: 1, label: "Timbale → snare with the wires off" },
  66: { piece: "snare-off", u: 0.5, v: 0.1, g: 0.9, label: "Timbale → snare with the wires off" },
};
