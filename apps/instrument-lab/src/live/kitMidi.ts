import { DRUM_MAP, stickFactor, type PieceId } from "../midi/drums";
import { HIHAT_CLOSED_FORCE_N, HIHAT_TRANSITION_SECONDS } from "./hiHat";
import type { MidiChannelMessage } from "./midiInput";

/** The seven physical pieces present in the playable kit. */
export type LiveKitPieceId = Exclude<PieceId, "snare-off">;

export interface KitMidiActions {
  /** u/v are fractions of radius; g multiplies the piece's standard stick speed. */
  strike(piece: LiveKitPieceId, u: number, v: number, g: number, frame: number): void;
  pedal(forceN: number, transitionSeconds: number, frame: number): void;
  unsupported?(note: number): void;
}

const byte = (value: number): boolean => Number.isInteger(value) && value >= 0 && value <= 127;
const channelNumber = (value: number): boolean => Number.isInteger(value) && value >= 0 && value <= 15;

/**
 * A hardware performance, without the file arranger's thinning or randomization.
 * CC4 prescribes 0..80 N on the existing hi-hat rod. Once that channel supplies
 * CC4, note articulations cannot overwrite its half-open pedal position.
 */
export class LiveKitMidi {
  private channel: number | null = 9;
  private readonly continuousPedal = new Set<number>();
  private readonly actions: KitMidiActions;

  constructor(actions: KitMidiActions) { this.actions = actions; }

  /** Zero-based channel, or null for all channels. Changing selection forgets controllers. */
  setChannel(channel: number | null): void {
    if (channel !== null && !channelNumber(channel)) throw new RangeError("MIDI channel must be 0..15 or null for all channels");
    if (channel !== this.channel) {
      this.channel = channel;
      this.reset();
    }
  }

  /** Forget input state on disconnect/stop; the physical instrument keeps ringing as it is. */
  reset(_frame?: number): void { this.continuousPedal.clear(); }

  /** Both actions of an articulated strike use the caller's one audio-frame timestamp. */
  receive(message: MidiChannelMessage, frame: number): void {
    if (!channelNumber(message.channel) || (this.channel !== null && message.channel !== this.channel)) return;
    if (!Number.isFinite(frame) || frame < 0) throw new RangeError("MIDI performance frame must be finite and nonnegative");
    if (message.kind === "cc") {
      if (!byte(message.controller) || !byte(message.value)) return;
      if (message.controller === 4) {
        this.continuousPedal.add(message.channel);
        this.actions.pedal((message.value / 127) * HIHAT_CLOSED_FORCE_N, HIHAT_TRANSITION_SECONDS, frame);
      } else if (message.controller === 120 || message.controller === 121 || message.controller === 123) {
        // Input reset, not a fake acoustic choke or a new foot/stick impulse.
        this.continuousPedal.delete(message.channel);
      }
      return;
    }
    if (message.kind !== "noteon" || !byte(message.key) || !byte(message.velocity) || message.velocity === 0) return;
    const rule = DRUM_MAP[message.key];
    if (!rule || rule.piece === "snare-off") {
      this.actions.unsupported?.(message.key);
      return;
    }
    if (rule.piece === "hihat" && !this.continuousPedal.has(message.channel)) {
      this.actions.pedal(message.key === 46 ? 0 : HIHAT_CLOSED_FORCE_N, HIHAT_TRANSITION_SECONDS, frame);
    }
    // A foot note never launches a stick, including when CC4 already moved it.
    if (message.key === 44) return;
    this.actions.strike(rule.piece, rule.u, rule.v, rule.g * stickFactor(message.velocity), frame);
  }
}
