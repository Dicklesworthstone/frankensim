// fs-live-mix: sums frame-stamped audio from several live engines in lock-step with the
// AudioContext clock, so a piano and a drum kit computed in different Workers stay
// sample-aligned with each other and with anything scheduled on the same context.
//
// Like fs-live-sink it adds no gain, limiter, fade or resampling: a lane's samples are summed
// exactly as its engine produced them. A sample that arrives after its frame has played is
// dropped and counted (late); a lane that sent nothing for a frame contributes silence.
//
// Protocol
//   main -> node.port   {type:'lane', id, port}   adopt a producer MessagePort as lane `id`
//   main -> node.port   {type:'drop', id}         forget a lane
//   lane -> mix         {type:'audio', frame, channels:[Float32Array, ...]}
//                       `frame` = context frame of the first sample (transfer the buffers).
//                       1 channel is copied to every output channel; N channels map 1:1.
//   mix  -> lane        {type:'clock', frame, sampleRate}   every 2 render quanta:
//                       the context frame at which the NEXT quantum starts.
//   mix  -> node.port   {type:'stats', frame, lanes:{id:{queuedFrames, lateFrames, gapFrames}}}
//                       about every 0.1 s. gapFrames: silence inside a lane that was sending.

class FsLiveMix extends AudioWorkletProcessor {
  constructor() {
    super();
    this.lanes = new Map(); // id -> {port, queue:[{frame, channels, len}], lateFrames, gapFrames, streaming}
    this.quantum = 0;
    this.port.onmessage = (event) => this.control(event.data);
  }

  control(message) {
    if (!message || typeof message !== 'object') return;
    if (message.type === 'lane' && message.port) {
      const lane = { port: message.port, queue: [], lateFrames: 0, gapFrames: 0, streaming: false };
      message.port.onmessage = (event) => this.receive(lane, event.data);
      this.lanes.set(message.id, lane);
      message.port.postMessage({ type: 'clock', frame: currentFrame, sampleRate });
    } else if (message.type === 'drop') {
      const lane = this.lanes.get(message.id);
      if (lane) {
        lane.port.onmessage = null;
        this.lanes.delete(message.id);
      }
    }
  }

  receive(lane, message) {
    if (!message || message.type !== 'audio' || !Array.isArray(message.channels) || !message.channels.length) return;
    const len = message.channels[0].length;
    const frame = Math.round(message.frame);
    if (!len || !Number.isFinite(frame)) return;
    if (frame + len <= currentFrame) {
      lane.lateFrames += len;
      return;
    }
    const q = lane.queue;
    // producers send in order; keep the queue sorted even if one does not
    let i = q.length;
    while (i > 0 && q[i - 1].frame > frame) i--;
    q.splice(i, 0, { frame, channels: message.channels, len });
  }

  process(_inputs, outputs) {
    const out = outputs[0];
    const n = out[0].length;
    const f0 = currentFrame;
    const f1 = f0 + n;
    for (const lane of this.lanes.values()) {
      const q = lane.queue;
      let covered = 0;
      for (let k = 0; k < q.length; k++) {
        const c = q[k];
        if (c.frame >= f1) break;
        const a = Math.max(f0, c.frame);
        const b = Math.min(f1, c.frame + c.len);
        if (b <= a) continue;
        for (let ch = 0; ch < out.length; ch++) {
          const src = c.channels[c.channels.length === 1 ? 0 : Math.min(ch, c.channels.length - 1)];
          const dst = out[ch];
          for (let i = a; i < b; i++) dst[i - f0] += src[i - c.frame];
        }
        covered += b - a;
      }
      // a hole: frames missing while the lane was streaming and still has audio queued after them
      // (a lane that falls silent on purpose simply stops sending and leaves its queue empty)
      if (covered < n && (covered > 0 || lane.streaming) && q.some((c) => c.frame + c.len > f1)) lane.gapFrames += n - covered;
      lane.streaming = covered === n;
      while (q.length && q[0].frame + q[0].len <= f1) q.shift();
    }
    this.quantum += 1;
    if (this.quantum % 2 === 0) {
      for (const lane of this.lanes.values()) lane.port.postMessage({ type: 'clock', frame: f1, sampleRate });
    }
    if (this.quantum % 38 === 0) {
      const lanes = {};
      for (const [id, lane] of this.lanes) {
        let queued = 0;
        for (const c of lane.queue) queued += c.len;
        lanes[id] = { queuedFrames: queued, lateFrames: lane.lateFrames, gapFrames: lane.gapFrames };
      }
      this.port.postMessage({ type: 'stats', frame: f1, lanes });
    }
    return true;
  }
}

registerProcessor('fs-live-mix', FsLiveMix);
