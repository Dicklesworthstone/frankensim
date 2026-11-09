// AudioWorklet sink for the virtual condenser microphone.
// The main thread renders mono 48 kHz blocks from the wasm engine (PocketWatchSim::render_audio_block)
// and posts them here; this processor plays them from a ring buffer and reports its fill level so the
// producer can keep latency near a target instead of guessing from wall-clock time.

class MicRingProcessor extends AudioWorkletProcessor {
  constructor() {
    super();
    this.capacity = 48000; // 1 s ring
    this.ring = new Float32Array(this.capacity);
    this.readIdx = 0;
    this.writeIdx = 0;
    this.available = 0;
    this.underruns = 0;
    this.quanta = 0;
    this.lastSample = 0;
    this.port.onmessage = (ev) => {
      const msg = ev.data;
      if (msg && msg.type === "chunk" && msg.samples) {
        this.push(msg.samples);
      } else if (msg && msg.type === "flush") {
        this.readIdx = 0;
        this.writeIdx = 0;
        this.available = 0;
      }
    };
  }

  push(samples) {
    const n = samples.length;
    for (let i = 0; i < n; i++) {
      this.ring[this.writeIdx] = samples[i];
      this.writeIdx = (this.writeIdx + 1) % this.capacity;
    }
    this.available += n;
    if (this.available > this.capacity) {
      // Overflow: drop the oldest audio, keep the freshest.
      const drop = this.available - this.capacity;
      this.readIdx = (this.readIdx + drop) % this.capacity;
      this.available = this.capacity;
    }
  }

  process(_inputs, outputs) {
    const out = outputs[0];
    if (!out || out.length === 0) return true;
    const ch0 = out[0];
    const n = ch0.length;
    let underrun = false;
    for (let i = 0; i < n; i++) {
      let s;
      if (this.available > 0) {
        s = this.ring[this.readIdx];
        this.readIdx = (this.readIdx + 1) % this.capacity;
        this.available--;
        this.lastSample = s;
      } else {
        // Decay toward zero instead of clicking when the producer falls behind.
        this.lastSample *= 0.97;
        s = this.lastSample;
        underrun = true;
      }
      ch0[i] = s;
    }
    for (let c = 1; c < out.length; c++) out[c].set(ch0);
    if (underrun) this.underruns++;
    if (++this.quanta % 4 === 0) {
      this.port.postMessage({ type: "level", available: this.available, underruns: this.underruns });
    }
    return true;
  }
}

registerProcessor("fs-mic-ring", MicRingProcessor);
