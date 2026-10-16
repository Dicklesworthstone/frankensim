// fs-live-mix: sums frame-stamped audio from several live engines in lock-step with the
// AudioContext clock, so a piano and a drum kit computed in different Workers stay
// sample-aligned with each other and with anything scheduled on the same context.
//
// Like fs-live-sink it adds no gain, limiter, fade or resampling: a lane's samples are summed
// exactly as its engine produced them. A sample that arrives after its frame has played is
// dropped and counted (late); a lane that sent nothing for a frame contributes silence.
//
// Protocol
//   main -> node.port   {type:'lane', id, port}   adopt a producer MessagePort as lane `id` (a lane
//                       already under `id` is dropped first)
//   main -> node.port   {type:'drop', id}         forget a lane and close its port
//   main -> node.port   {type:'dispose'}          drop every lane; process() then returns false, so
//                       the disconnected node can be collected
//   main -> node.port   {type:'flush', id?, from?, epoch?}  discard queued audio of lane `id`
//                       (every lane when omitted) at context frames >= `from` (default: now).
//                       With `epoch`, only audio stamped with an older epoch is discarded, and
//                       older audio still in flight is refused when it arrives. Used to stop a
//                       fixed score that the lanes rendered far ahead; nothing is faded here.
//   lane -> mix         {type:'audio', frame, channels:[Float32Array, ...], epoch?, rest?}
//                       `frame` = context frame of the first sample (transfer the buffers).
//                       1 channel is copied to every output channel; N channels map 1:1.
//                       `rest`: true when the lane sent nothing on purpose since its previous
//                       message (its engine was silent and not stepped, or it starts here);
//                       false when this message continues that stream, so any frames between the
//                       two are a hole. A lane that never sends `rest` keeps the older gap rule
//                       (see process), so an older lane.js still gets the counts it got before.
//   mix  -> lane        {type:'clock', frame, sampleRate, wallMs}   every 2 render quanta:
//                       the context frame at which the NEXT quantum starts, and Date.now() when
//                       the tick was sent (Date is the one clock this scope shares with Workers;
//                       it has no `performance`), so a lane can tell how long a tick queued.
//   mix  -> node.port   {type:'stats', frame, lanes:{id:{queuedFrames, lateFrames, gapFrames}}}
//                       about every 0.1 s.
//                       lateFrames: frames that arrived after their frame had played, so were
//                       dropped: whole messages, and the played-out head of a partly late one.
//                       gapFrames: frames a streaming lane never sent (it fell behind and stepped
//                       ahead), counted when the message after the hole arrives. A lane that
//                       rests and wakes later is not a gap, and after a flush the next message
//                       starts afresh. A frame is never both, so late + gap is the audio lost.

class FsLiveMix extends AudioWorkletProcessor {
  constructor() {
    super();
    // id -> {port, queue:[{frame, channels, len}], lateFrames, gapFrames,
    //        end: frame after the lane's last message, marked: it sends `rest`, streaming: older rule}
    this.lanes = new Map();
    this.quantum = 0;
    this.alive = true;
    this.port.onmessage = (event) => this.control(event.data);
  }

  // forget lane `id`: its queued audio and its port (a port left open would keep receiving, and its
  // audio would pile up in a queue nothing plays)
  forget(id) {
    const lane = this.lanes.get(id);
    if (!lane) return;
    lane.port.onmessage = null;
    lane.port.close();
    lane.queue = [];
    this.lanes.delete(id);
  }

  control(message) {
    if (!message || typeof message !== 'object') return;
    if (message.type === 'lane' && message.port) {
      this.forget(message.id);
      const lane = { port: message.port, queue: [], lateFrames: 0, gapFrames: 0, end: undefined, marked: false, streaming: false };
      message.port.onmessage = (event) => this.receive(lane, event.data);
      this.lanes.set(message.id, lane);
      message.port.postMessage({ type: 'clock', frame: currentFrame, sampleRate, wallMs: Date.now() });
    } else if (message.type === 'drop') {
      this.forget(message.id);
    } else if (message.type === 'dispose') {
      for (const id of [...this.lanes.keys()]) this.forget(id);
      this.alive = false;
      this.port.onmessage = null;
    } else if (message.type === 'flush') {
      const from = Number.isFinite(message.from) ? Math.round(message.from) : currentFrame;
      const epoch = Number.isFinite(message.epoch) ? message.epoch : undefined;
      for (const [id, lane] of this.lanes) {
        if (message.id !== undefined && message.id !== id) continue;
        const kept = [];
        for (const c of lane.queue) {
          if (c.frame + c.len <= from || (epoch !== undefined && c.epoch >= epoch)) kept.push(c);
          else if (c.frame < from) {
            const len = from - c.frame;
            kept.push({ frame: c.frame, len, epoch: c.epoch, channels: c.channels.map((x) => x.subarray(0, len)) });
          }
        }
        lane.queue = kept;
        lane.end = undefined; // what follows a flush is not measured against what came before
        lane.streaming = false;
        if (epoch !== undefined) lane.epoch = epoch;
      }
    }
  }

  receive(lane, message) {
    if (!message || message.type !== 'audio' || !Array.isArray(message.channels) || !message.channels.length) return;
    const len = message.channels[0].length;
    const frame = Math.round(message.frame);
    if (!len || !Number.isFinite(frame)) return;
    const epoch = Number.isFinite(message.epoch) ? message.epoch : 0;
    if (lane.epoch !== undefined && epoch < lane.epoch) return; // rendered before a flush: stale
    if (typeof message.rest === 'boolean') {
      lane.marked = true;
      if (!message.rest && lane.end !== undefined && frame > lane.end) lane.gapFrames += frame - lane.end;
    }
    lane.end = frame + len;
    // currentFrame is the first frame of the next quantum, so frames before it have played out
    // without this message: drop and count them, a whole message or its head
    const late = currentFrame - frame;
    if (late > 0) {
      lane.lateFrames += late < len ? late : len;
      if (late >= len) return;
    }
    const q = lane.queue;
    // producers send in order; keep the queue sorted even if one does not
    let i = q.length;
    while (i > 0 && q[i - 1].frame > frame) i--;
    q.splice(i, 0, { frame, channels: message.channels, len, epoch });
  }

  process(_inputs, outputs) {
    if (!this.alive) return false;
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
      if (!lane.marked) {
        // older rule, for a lane that does not say when it rests: frames missing from a quantum
        // while the lane was streaming and still has audio queued after them
        if (covered < n && (covered > 0 || lane.streaming) && q.some((c) => c.frame + c.len > f1)) lane.gapFrames += n - covered;
        lane.streaming = covered === n;
      }
      while (q.length && q[0].frame + q[0].len <= f1) q.shift();
    }
    this.quantum += 1;
    if (this.quantum % 2 === 0) {
      const wallMs = Date.now();
      for (const lane of this.lanes.values()) lane.port.postMessage({ type: 'clock', frame: f1, sampleRate, wallMs });
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
