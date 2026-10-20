// fs-live-sink: a dumb, honest audio sink for live FrankenSim engines.
//
// It plays Float32 chunks exactly as received. No gain, no limiter, no fade,
// no resampling, no concealment: when the FIFO runs dry it outputs exact
// zeros and counts the underrun.
//
// Protocol (identical on the node's own port and on any port handed over with
// {type:'connect', port}; a Worker normally owns such a port):
//
//   in   {type:'audio', channels:[Float32Array, ...]}   transfer the buffers.
//        1 channel is copied to every output channel; N channels map 1:1.
//   in   {type:'flush', epoch?: number, primeFrames?: number}
//        drop everything queued and start a new epoch (used on rebuild: the
//        old instrument stops, the new one starts silent). Underrun counters
//        restart with the epoch. Playback of the new epoch begins once
//        `primeFrames` frames are queued (default 0); the wait is silence and
//        is not an underrun.
//   in   {type:'connect', port: MessagePort}            adopt another producer.
//   in   {type:'config', levelEveryQuanta?: number}     level cadence (default 2).
//   in   {type:'dispose'}   close producer ports and end processing permanently.
//   out  {type:'level', epoch, queuedFrames, receivedFrames, playedFrames,
//         underruns, underrunFrames, currentFrame, currentTime, sampleRate}
//        receivedFrames is cumulative over the node's life; playedFrames,
//        underruns and underrunFrames count since the last flush (the epoch).
//
// `receivedFrames` lets a producer compute what is still in flight:
//   estimatedQueue = queuedFrames + (framesSent - receivedFrames).

class FsLiveSink extends AudioWorkletProcessor {
  constructor() {
    super();
    this.queue = [];          // [{channels:[Float32Array], offset:number}]
    this.queuedFrames = 0;
    this.receivedFrames = 0;
    this.playedFrames = 0;
    this.underruns = 0;       // quanta that had to be padded with zeros
    this.underrunFrames = 0;  // zero frames emitted because the FIFO was empty
    this.started = false;     // no underrun is counted before playback starts
    this.primeFrames = 0;     // frames that must be queued before playback starts
    this.epoch = 0;
    this.levelEvery = 2;
    this.quantum = 0;
    this.ports = [this.port];
    this.alive = true;
    this.port.onmessage = (event) => this.receive(event.data);
  }

  receive(message) {
    if (!this.alive || !message || typeof message !== 'object') return;
    switch (message.type) {
      case 'audio': {
        const channels = message.channels;
        if (!Array.isArray(channels) || channels.length === 0) return;
        const frames = channels[0].length;
        if (frames === 0) return;
        this.queue.push({ channels, offset: 0 });
        this.queuedFrames += frames;
        this.receivedFrames += frames;
        if (this.queuedFrames >= this.primeFrames) this.started = true;
        break;
      }
      case 'flush':
        this.queue.length = 0;
        this.queuedFrames = 0;
        this.playedFrames = 0;
        this.underruns = 0;
        this.underrunFrames = 0;
        this.started = false;
        this.primeFrames = Number.isFinite(message.primeFrames) ? Math.max(0, message.primeFrames) : 0;
        this.epoch = Number.isFinite(message.epoch) ? message.epoch : this.epoch + 1;
        this.postLevel();
        break;
      case 'connect':
        if (message.port) {
          message.port.onmessage = (event) => this.receive(event.data);
          this.ports.push(message.port);
          this.postLevel();
        }
        break;
      case 'config':
        if (Number.isFinite(message.levelEveryQuanta) && message.levelEveryQuanta >= 1) {
          this.levelEvery = Math.floor(message.levelEveryQuanta);
        }
        break;
      case 'dispose':
        this.alive = false;
        this.queue.length = 0;
        this.queuedFrames = 0;
        for (const port of this.ports) {
          port.onmessage = null;
          port.close();
        }
        this.ports = [];
        break;
      default:
        break;
    }
  }

  postLevel() {
    const level = {
      type: 'level',
      epoch: this.epoch,
      queuedFrames: this.queuedFrames,
      receivedFrames: this.receivedFrames,
      playedFrames: this.playedFrames,
      underruns: this.underruns,
      underrunFrames: this.underrunFrames,
      currentFrame,
      currentTime,
      sampleRate,
    };
    for (const port of this.ports) port.postMessage(level);
  }

  process(_inputs, outputs) {
    if (!this.alive) return false;
    const output = outputs[0];
    const frames = output[0].length;
    let written = 0;
    while (this.started && written < frames && this.queue.length > 0) {
      const head = this.queue[0];
      const available = head.channels[0].length - head.offset;
      const take = Math.min(available, frames - written);
      for (let c = 0; c < output.length; c++) {
        const source = head.channels[head.channels.length === 1 ? 0 : Math.min(c, head.channels.length - 1)];
        output[c].set(source.subarray(head.offset, head.offset + take), written);
      }
      head.offset += take;
      written += take;
      if (head.offset >= head.channels[0].length) this.queue.shift();
    }
    this.queuedFrames -= written;
    this.playedFrames += written;
    if (written < frames) {
      // The AudioWorklet hands us zero-filled buffers: the remainder is silence.
      if (this.started) {
        this.underruns += 1;
        this.underrunFrames += frames - written;
      }
    }
    this.quantum += 1;
    if (this.quantum % this.levelEvery === 0) this.postLevel();
    return true;
  }
}

registerProcessor('fs-live-sink', FsLiveSink);
