import type { MidiChannelMessage } from './midiInput';

export type PianoMidiChannel = 'pitched' | number;

export interface PianoMidiPedals {
  sustain: number;
  sostenuto: boolean;
  unaCorda: boolean;
}

export interface PianoMidiTarget {
  /** Every note-on is a strike; acquireHold is true only for the first MIDI owner of this key. */
  noteOn: (key: number, velocity: number, acquireHold: boolean) => void;
  /** Called only after the last MIDI channel holding this key has released it. */
  noteOff: (key: number) => void;
  pedals: (pedals: PianoMidiPedals) => void;
}

interface ChannelState extends PianoMidiPedals {
  keys: Set<number>;
}

const pedalsUp = (): PianoMidiPedals => ({ sustain: 0, sostenuto: false, unaCorda: false });

/** One physical piano has one set of pedals: every input may keep a pedal down. */
export function combinePianoPedals(a: PianoMidiPedals, b: PianoMidiPedals): PianoMidiPedals {
  return {
    sustain: Math.max(a.sustain, b.sustain),
    sostenuto: a.sostenuto || b.sostenuto,
    unaCorda: a.unaCorda || b.unaCorda,
  };
}

/**
 * MIDI owns at most one hold per channel/key, and one section hold per key.
 * Repeated note-ons retrigger without accumulating unmatched holds. Channels
 * share the existing piano's keys and pedals, rather than independent voices.
 */
export class PianoMidi {
  private channel: PianoMidiChannel = 'pitched';
  private readonly channels = new Map<number, ChannelState>();
  private pedals = pedalsUp();

  constructor(private readonly target: PianoMidiTarget) {}

  selectChannel(channel: PianoMidiChannel): void {
    if (channel !== 'pitched' && (!Number.isInteger(channel) || channel < 0 || channel > 15)) {
      throw new RangeError('MIDI channel must be 0–15, or all pitched channels.');
    }
    if (channel === this.channel) return;
    this.reset();
    this.channel = channel;
  }

  receive(message: MidiChannelMessage): void {
    if (this.channel === 'pitched' ? message.channel === 9 : message.channel !== this.channel) return;
    if (message.kind !== 'cc') {
      if (message.key < 21 || message.key > 108) return;
      if (message.kind === 'noteoff' || message.velocity === 0) {
        const state = this.channels.get(message.channel);
        if (state?.keys.delete(message.key) && !this.held(message.key)) this.target.noteOff(message.key);
      } else {
        const acquireHold = !this.held(message.key);
        this.state(message.channel).keys.add(message.key);
        this.target.noteOn(message.key, message.velocity, acquireHold);
      }
      return;
    }
    const { channel, controller, value } = message;
    if (controller === 123 || controller === 120) {
      this.releaseChannel(channel);
      // All Sound Off releases this channel's pedal claims as well. Acoustic
      // motion and other inputs' holds remain in the one shared physical piano.
      if (controller === 120) this.resetControllers(channel);
    } else if (controller === 121) {
      this.resetControllers(channel);
    } else if (controller === 64 || controller === 66 || controller === 67) {
      const state = this.state(channel);
      if (controller === 64) state.sustain = value / 127;
      else if (controller === 66) state.sostenuto = value >= 64;
      else state.unaCorda = value >= 64;
      this.publishPedals();
    }
  }

  /** Disconnect, source/channel change, or the section's existing all-input silence. */
  reset(): void {
    const keys = new Set<number>();
    for (const state of this.channels.values()) for (const key of state.keys) keys.add(key);
    this.channels.clear();
    for (const key of keys) this.target.noteOff(key);
    this.publishPedals();
  }

  private state(channel: number): ChannelState {
    let state = this.channels.get(channel);
    if (!state) {
      state = { keys: new Set(), ...pedalsUp() };
      this.channels.set(channel, state);
    }
    return state;
  }

  private held(key: number): boolean {
    for (const state of this.channels.values()) if (state.keys.has(key)) return true;
    return false;
  }

  private releaseChannel(channel: number): void {
    const state = this.channels.get(channel);
    if (!state) return;
    const keys = [...state.keys];
    state.keys.clear();
    for (const key of keys) if (!this.held(key)) this.target.noteOff(key);
  }

  private resetControllers(channel: number): void {
    const state = this.channels.get(channel);
    if (!state) return;
    Object.assign(state, pedalsUp());
    this.publishPedals();
  }

  private publishPedals(): void {
    let pedals = pedalsUp();
    for (const state of this.channels.values()) pedals = combinePianoPedals(pedals, state);
    if (pedals.sustain === this.pedals.sustain && pedals.sostenuto === this.pedals.sostenuto
      && pedals.unaCorda === this.pedals.unaCorda) return;
    this.pedals = pedals;
    this.target.pedals({ ...pedals });
  }
}
