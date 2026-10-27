// Playback of two already-rendered native results, never a second synthesizer.
import {compareNativePressure, nativePcmSamples} from './native-wave.mjs';

export class NativeComparisonPlayback {
  constructor(createContext = () => new AudioContext({sampleRate:48000})) {
    this.createContext = createContext;
    this.context = null;
    this.data = null;
    this.buffers = null;
    this.nodes = null;
    this.epoch = 0;
    this.side = 'a';
    this.volume = 0.25;
    this.duration = 0;
    this.startedAt = 0;
    this.offsetSeconds = 0;
    this.mix = {from:0, to:0, start:0, end:0};
  }

  // Admission and decoding finish before replacing a playable pair. Copies of
  // the PCM samples protect the reference from later editor/attachment changes.
  setPair(aBytes, aReceipt, bBytes, bReceipt) {
    const comparison = compareNativePressure(aBytes, aReceipt, bBytes, bReceipt);
    const data = [nativePcmSamples(aBytes), nativePcmSamples(bBytes)]
      .map(samples => samples.slice(0, comparison.samples));
    this.stop();
    this.data = data;
    this.buffers = null;
    this.duration = comparison.durationSeconds;
    return comparison;
  }

  clear() {
    this.stop(); this.data = null; this.buffers = null; this.duration = 0;
  }

  get playing() { return this.nodes !== null; }
  get positionSeconds() {
    return this.playing ? Math.max(this.offsetSeconds, Math.min(this.duration, this.context.currentTime - this.startedAt)) : 0;
  }

  async play(offsetSeconds = 0) {
    if (!this.data || !Number.isFinite(offsetSeconds) || offsetSeconds < 0 || offsetSeconds >= this.duration) {
      throw new Error('Choose a position inside a calibrated comparison.');
    }
    this.stop();
    const epoch = this.epoch;
    const context = this.context ??= this.createContext();
    await context.resume();
    // Stop, a new result, or leaving the page can interrupt a pending resume.
    if (epoch !== this.epoch) return false;
    if (!this.buffers) {
      this.buffers = this.data.map(samples => {
        const buffer = context.createBuffer(1, samples.length, 48000);
        buffer.copyToChannel(samples, 0); return buffer;
      });
    }
    const master = context.createGain();
    master.gain.value = this.volume; master.connect(context.destination);
    const nodes = {master, sources:[], gains:[]};
    this.nodes = nodes;
    const start = context.currentTime + 0.02;
    this.startedAt = start - offsetSeconds;
    this.offsetSeconds = offsetSeconds;
    const mix = this.side === 'b' ? 1 : 0;
    this.mix = {from:mix, to:mix, start, end:start};
    try {
      for (let index = 0; index < 2; index++) {
        const source = context.createBufferSource(), gain = context.createGain();
        nodes.sources.push(source); nodes.gains.push(gain);
        source.buffer = this.buffers[index];
        gain.gain.value = index === 0 ? 1 - mix : mix;
        source.connect(gain); gain.connect(master);
        // Both sources share one absolute start, offset and finite window.
        // Switching the audible side never restarts either physical timeline.
        source.start(start, offsetSeconds, this.duration - offsetSeconds);
      }
      nodes.sources[0].onended = () => { if (this.nodes === nodes) this.stop(); };
      return true;
    } catch (error) { this.stop(); throw error; }
  }

  select(side) {
    if (side !== 'a' && side !== 'b') throw new Error('Comparison side must be A or B.');
    this.side = side;
    if (!this.nodes) return;
    const now = this.context.currentTime, {from, to, start, end} = this.mix;
    const amount = end > start ? Math.max(0, Math.min(1, (now - start) / (end - start))) : 1;
    const present = from + amount * (to - from), next = side === 'b' ? 1 : 0;
    // A complementary 5 ms crossfade prevents switch clicks. It is playback
    // only: original WAVs, pressure measurements and native clocks are untouched.
    this.mix = {from:present, to:next, start:now, end:now + 0.005};
    this.nodes.gains.forEach((node, index) => {
      node.gain.cancelScheduledValues(now);
      node.gain.setValueAtTime(index === 0 ? 1 - present : present, now);
      node.gain.linearRampToValueAtTime(index === 0 ? 1 - next : next, now + 0.005);
    });
  }

  setVolume(value) {
    if (!Number.isFinite(value) || value < 0 || value > 1) throw new Error('Shared playback volume must be between zero and one.');
    this.volume = value;
    if (this.nodes) this.nodes.master.gain.setValueAtTime(value, this.context.currentTime);
  }

  stop() {
    this.epoch++;
    const nodes = this.nodes; this.nodes = null;
    if (!nodes) return;
    for (const source of nodes.sources) {
      source.onended = null;
      try { source.stop(); } catch { /* A completed source can already be stopped. */ }
      source.disconnect();
    }
    for (const gain of nodes.gains) gain.disconnect();
    nodes.master.disconnect();
  }

  async close() {
    this.stop();
    const context = this.context; this.context = null; this.buffers = null;
    if (context && context.state !== 'closed') await context.close();
  }
}
