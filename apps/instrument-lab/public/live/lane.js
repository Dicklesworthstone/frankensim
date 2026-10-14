// Pacing for a live engine Worker that feeds one lane of `fs-live-mix` (mix-worklet.js).
//
// The Worker renders a fixed distance ahead of the AudioContext clock and stamps every block
// with the context frame its first sample plays at, so lanes computed in different Workers
// stay sample-aligned. Events (a stroke, a key) carry the frame they act on; an event without
// one acts on the next rendered sample. A sequence repeats a list of events every
// `loopFrames`, starting at `startFrame`.
//
// While the engine reports itself idle (a silent drum is not stepped) nothing is rendered and
// the lane jumps ahead with the clock, so the next stroke lands `aheadFrames` after it arrives.

export function createLane(port, engine, { aheadFrames = 1024, chunk = 128 } = {}) {
  let clock = { frame: 0, at: 0, rate: 48000, seen: false };
  let next = -1;                 // context frame of the next sample to render
  const queue = [];              // pending events, sorted by frame (undefined frame = now)
  let seq = null;                // {startFrame, loopFrames, events, loop}
  let pumping = false;
  const stats = { renderedFrames: 0, skippedFrames: 0, lateEvents: 0, wallMs: 0, audioMs: 0 };

  port.onmessage = (event) => {
    const m = event.data;
    if (m && m.type === 'clock') {
      // ticks that queued up while this Worker was busy are stale: keep the extrapolated clock
      // unless the tick is current, or the estimate is old enough to need re-anchoring anyway
      const at = performance.now();
      if (!clock.seen || m.frame > now() - 1024 || at - clock.at > 50) clock = { frame: m.frame, at, rate: m.sampleRate || 48000, seen: true };
      pump();
    }
  };

  const now = () => clock.frame + ((performance.now() - clock.at) * clock.rate) / 1000;

  function insert(e) {
    let i = queue.length;
    while (i > 0 && (queue[i - 1].frame ?? -Infinity) > (e.frame ?? -Infinity)) i--;
    queue.splice(i, 0, e);
  }

  // queue every sequenced event whose frame falls before `until`
  function extend(until) {
    if (!seq) return;
    while (seq.startFrame + seq.loop * seq.loopFrames < until) {
      const base = seq.startFrame + seq.loop * seq.loopFrames;
      for (const e of seq.events) insert({ ...e, frame: base + e.offset });
      seq.loop += 1;
    }
  }

  function pump() {
    if (pumping || !clock.seen) return;
    pumping = true;
    try {
      const t = now();
      const horizon = t + aheadFrames;
      extend(horizon + chunk);
      const due = queue.length ? queue[0].frame : undefined;
      if (engine.idle() && (due === undefined || due > horizon) && !queue.some((e) => e.frame === undefined)) {
        // silent and nothing due within reach: keep pace with the clock without stepping
        // (never past a scheduled event, which must still land on its own frame)
        if (next < horizon) next = Math.round(due !== undefined ? Math.min(horizon, due) : horizon);
        return;
      }
      if (next < t + chunk) {
        // never rendered, woken from idle, or fallen behind: resume a safe distance ahead
        const resume = Math.round(t + aheadFrames / 2);
        if (next >= 0 && !engine.idle()) stats.skippedFrames += Math.max(0, resume - next);
        next = resume;
      }
      if (engine.idle() && due !== undefined && due > next) next = Math.round(due); // skip silence up to the event
      while (next < horizon) {
        // act on everything due at this sample, then render up to the next event or a full chunk
        while (queue.length && (queue[0].frame === undefined || queue[0].frame <= next)) {
          const e = queue.shift();
          if (e.frame !== undefined && e.frame < next) stats.lateEvents += 1;
          engine.apply(e);
        }
        let n = chunk;
        if (queue.length && queue[0].frame !== undefined) n = Math.max(1, Math.min(n, Math.round(queue[0].frame) - next));
        const t0 = performance.now();
        const samples = engine.render(n);
        stats.wallMs += performance.now() - t0;
        stats.audioMs += (n * 1000) / clock.rate;
        stats.renderedFrames += n;
        port.postMessage({ type: 'audio', frame: next, channels: [samples] }, [samples.buffer]);
        next += n;
        if (engine.idle() && !queue.length && !seq) break;
      }
    } finally {
      pumping = false;
    }
  }

  // the clock arrives every ~5 ms while the context runs; this keeps a lane alive between ticks
  setInterval(pump, 10);

  return {
    stats,
    /** Next frame to be rendered (the context frame the current engine state belongs to). */
    get frame() { return next; },
    now,
    event(e) { insert(e); pump(); },
    sequence(startFrame, loopFrames, events) {
      for (let i = queue.length - 1; i >= 0; i--) if (queue[i].seq) queue.splice(i, 1);
      const sorted = events.map((e) => ({ ...e, seq: true })).sort((a, b) => a.offset - b.offset);
      seq = { startFrame: Math.round(startFrame), loopFrames: Math.round(loopFrames), events: sorted, loop: 0 };
      pump();
    },
    stopSequence() {
      seq = null;
      for (let i = queue.length - 1; i >= 0; i--) if (queue[i].seq) queue.splice(i, 1);
    },
    pump,
  };
}
