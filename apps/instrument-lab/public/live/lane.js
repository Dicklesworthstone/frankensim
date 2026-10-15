// Pacing for a live engine Worker that feeds one lane of `fs-live-mix` (mix-worklet.js).
//
// The Worker renders a fixed distance ahead of the AudioContext clock and stamps every block
// with the context frame its first sample plays at, so lanes computed in different Workers
// stay sample-aligned. The clock comes from the mixer's ticks, each advanced by the time it
// waited in this Worker's queue, so a busy Worker still knows when it is late (it then steps
// ahead and counts skippedFrames and lateEvents). Events (a stroke, a key) carry the frame they
// act on; an event without one acts on the next rendered sample. A sequence repeats a list of
// events every `loopFrames`, starting at `startFrame`.
//
// While the engine reports itself idle (a silent drum is not stepped) nothing is rendered and
// the lane jumps ahead with the clock, so the next stroke lands `aheadFrames` after it arrives.
// Each audio message says whether the silence before it was such a rest (`rest`), so the mixer
// counts as gaps only the frames a lane skipped after falling behind (= skippedFrames).
//
// A score is a fixed, non-looping list of events with absolute frames (a whole MIDI file): it is
// sorted once and fed into the queue only as far as the render horizon, so a 10,000-event score
// costs what a short loop does. `resync()` forgets how far ahead the lane has rendered, so after
// the mixer drops a lane's queued audio the next samples are rendered for the clock as it is now.
//
// Rendering happens in turns, one task each. A clock tick only moves the clock and asks for a
// turn, so ticks that queued while the Worker was busy are read in microseconds instead of
// costing a render each. A turn renders for at most TURN_MS of wall time and then yields through
// a message the lane posts to itself (no timer clamping), so messages from the page (seek, pause,
// a new score) and fresh ticks are handled between turns even while the lane is behind.
//
// The engine is called for BLOCK frames at a time (fewer up to the next event). `chunk` is the most
// frames one audio message carries: with a larger chunk the blocks are gathered into one message
// before it is sent, which changes nothing but the number of messages. A lane that renders a score
// a second ahead can send 512 frames at a time, a quarter of the tasks for the mixer's realtime
// thread; a partly filled message waits for its next blocks only while its first frame is at least
// a quarter of the look-ahead (and 4 chunks) ahead of the clock, so it never arrives late for it.

const TURN_MS = 16;
const BLOCK = 128;
// The clock is extrapolated from the last tick for at most this long. Ticks arrive every ~5 ms while the
// context renders, and a busy lane reads them between turns, so a longer silence means the context is
// not rendering (suspended, interrupted, its device gone): the lane then stops rendering ahead instead of
// filling the mixer's queue with audio for a clock that is not moving.
const EXTRAPOLATE_MS = 1000;

export function createLane(port, engine, { aheadFrames = 1024, chunk = BLOCK } = {}) {
  chunk = Math.max(BLOCK, Math.round(chunk) || BLOCK);
  const keepAhead = Math.max(4 * chunk, aheadFrames / 4);
  let clock = { frame: 0, at: 0, rate: 48000, seen: false };
  let next = -1;                 // context frame of the next sample to render
  const queue = [];              // pending events, sorted by frame (undefined frame = now)
  let seq = null;                // {startFrame, loopFrames, events, loop}
  let score = null;              // {events (sorted by frame), at}
  let epoch = 0;                 // stamped on every block; raised by resync (see mix-worklet flush)
  let pumping = false;
  let open = null;               // the audio message being gathered: {frame, epoch, rest, buf, len}
  let sentEnd = -1;              // context frame after the last block delivered
  let hole = false;              // the next block follows frames this lane skipped (it fell behind)
  const stats = { renderedFrames: 0, skippedFrames: 0, lateEvents: 0, wallMs: 0, audioMs: 0 };

  port.onmessage = (event) => {
    const m = event.data;
    if (m && m.type === 'clock') {
      // A tick that queued while this Worker was busy is stale by the time it is read. A tick
      // stamped with the wall clock (`wallMs`, Date.now() in the mixer) says by how much: advance
      // its frame by that age, less 1 ms for Date.now()'s whole-millisecond steps so the clock is
      // never put ahead, and "now" stays right however far behind the Worker falls (so the
      // fell-behind step below fires and counts). An unstamped tick is taken as current: keep the
      // extrapolated clock unless it is, or the estimate is old enough to need re-anchoring anyway.
      const at = performance.now();
      const rate = m.sampleRate || 48000;
      let frame = m.frame;
      if (Number.isFinite(m.wallMs)) {
        frame += (Math.max(0, Date.now() - m.wallMs - 1) * rate) / 1000;
        // Date.now() is wall time: a step of the system clock (NTP, waking from sleep) between the tick's
        // send and its read would date it minutes or hours ahead, and the lane would stop rendering until
        // the real clock got there. The context never runs ahead of this lane's own (monotonic) estimate by
        // more than a burst of quanta, so a tick dated further ahead than that keeps its own frame.
        if (clock.seen && frame > now() + rate / 4) frame = m.frame;
      }
      if (!clock.seen || frame > now() - 1024 || at - clock.at > 50) clock = { frame, at, rate, seen: true };
      schedulePump();
    }
  };

  const now = () => clock.frame + (Math.min(EXTRAPOLATE_MS, performance.now() - clock.at) * clock.rate) / 1000;

  function insert(e) {
    let i = queue.length;
    while (i > 0 && (queue[i - 1].frame ?? -Infinity) > (e.frame ?? -Infinity)) i--;
    queue.splice(i, 0, e);
  }

  // queue every sequenced or scored event whose frame falls before `until`
  function extend(until) {
    if (score) {
      const ev = score.events;
      while (score.at < ev.length && ev[score.at].frame < until) insert(ev[score.at++]);
      if (score.at >= ev.length) score = null;
    }
    if (!seq) return;
    while (seq.startFrame + seq.loop * seq.loopFrames < until) {
      const base = seq.startFrame + seq.loop * seq.loopFrames;
      for (const e of seq.events) insert({ ...e, frame: base + e.offset });
      seq.loop += 1;
    }
  }

  // at most one turn pending: queued behind whatever the page and the mixer have already posted
  let scheduled = false;
  const wake = new MessageChannel();
  wake.port1.onmessage = () => {
    scheduled = false;
    pump();
  };
  function schedulePump() {
    if (scheduled) return;
    scheduled = true;
    wake.port2.postMessage(0);
  }

  function send() {
    const s = open.len === open.buf.length ? open.buf : open.buf.subarray(0, open.len);
    port.postMessage({ type: 'audio', frame: open.frame, channels: [s], epoch: open.epoch, rest: open.rest }, [s.buffer]);
    open = null;
  }

  // a rendered block, to the mixer at once or gathered into the open message. `rest` tells the
  // mixer that the silence before this block is on purpose (the lane rested, or starts here); a
  // block after a skip is not a rest, so the mixer counts the skipped frames as a gap.
  function deliver(samples, frame) {
    const rest = frame !== sentEnd && !hole;
    hole = false;
    sentEnd = frame + samples.length;
    if (open && (open.frame + open.len !== frame || open.epoch !== epoch)) send();
    if (chunk === BLOCK) {
      port.postMessage({ type: 'audio', frame, channels: [samples], epoch, rest }, [samples.buffer]);
      return;
    }
    if (!open) open = { frame, epoch, rest, buf: new Float32Array(chunk), len: 0 };
    open.buf.set(samples, open.len);
    open.len += samples.length;
    if (open.len + BLOCK > chunk) send(); // the next block may not fit
  }

  function pump() {
    if (pumping || !clock.seen) return;
    pumping = true;
    try {
      const t = now();
      const lead = seq ? Math.min(seq.startFrame, t + aheadFrames) : t;
      const horizon = Math.max(t, lead) + aheadFrames;
      extend(horizon + BLOCK);
      const due = queue.length ? queue[0].frame : undefined;
      if (engine.idle() && (due === undefined || due > horizon) && !queue.some((e) => e.frame === undefined)) {
        // silent and nothing due within reach: keep pace with the clock without stepping
        // (never past a scheduled event, which must still land on its own frame). Nothing has been
        // sent past sentEnd, so this mark also comes back when the clock does (a context that stalled
        // while the estimate ran on): the next stroke still lands `aheadFrames` after it arrives
        next = Math.max(sentEnd, Math.round(due !== undefined ? Math.min(horizon, due) : horizon));
        return;
      }
      if (next < t + BLOCK) {
        if (next < 0 || engine.idle()) {
          // never rendered or woken from idle: start at the next scheduled event or a short lead-in
          next = due !== undefined && due >= t + BLOCK
            ? Math.round(Math.min(due, horizon))
            : Math.round(t + Math.min(aheadFrames / 2, 4096));
        } else {
          // momentarily fell behind while streaming: step ahead by a small ~42 ms cushion (never a half-second dropout)
          const resume = Math.round(t + Math.min(aheadFrames / 4, 2048));
          stats.skippedFrames += Math.max(0, resume - next);
          hole = resume > next;
          while (queue.length && queue[0].frame !== undefined && queue[0].frame < resume) {
            stats.lateEvents += 1;
            engine.apply(queue.shift());
          }
          next = resume;
        }
      }
      // skip silence up to the event; an idle lane has sent nothing past sentEnd, so it also steps back to an
      // event scheduled closer than the look-ahead it kept pace at (a sequence or score starting soon)
      if (engine.idle() && due !== undefined && due >= Math.max(sentEnd, t + BLOCK)) next = Math.round(due);
      const turnLimit = next + Math.max(BLOCK * 16, 2048);
      const stopAt = Math.min(horizon, turnLimit);
      const yieldAt = performance.now() + TURN_MS;
      while (next < stopAt) {
        // act on everything due at this sample, then render up to the next event or a full chunk
        while (queue.length && (queue[0].frame === undefined || queue[0].frame <= next)) {
          const e = queue.shift();
          if (e.frame !== undefined && e.frame < next) stats.lateEvents += 1;
          engine.apply(e);
        }
        let n = BLOCK;
        if (queue.length && queue[0].frame !== undefined) n = Math.max(1, Math.min(n, Math.round(queue[0].frame) - next));
        const t0 = performance.now();
        const samples = engine.render(n);
        const t1 = performance.now();
        stats.wallMs += t1 - t0;
        stats.audioMs += (n * 1000) / clock.rate;
        stats.renderedFrames += n;
        deliver(samples, next);
        next += n;
        if (engine.idle() && !queue.length && !seq) {
          if (open) send();
          break;
        }
        if (t1 >= yieldAt) break; // the turn has had its time; the rest is rendered after a yield
      }
      const nextDue = queue.length ? queue[0].frame : undefined;
      if (next < horizon && (!engine.idle() || (nextDue !== undefined && nextDue <= horizon))) {
        schedulePump();
      }
    } finally {
      pumping = false;
      if (open && open.frame < now() + keepAhead) send(); // close to its time: send what there is
    }
  }

  // the clock arrives every ~5 ms while the context runs; this keeps a lane alive between ticks
  setInterval(schedulePump, 10);

  return {
    stats,
    /** Next frame to be rendered (the context frame the current engine state belongs to). */
    get frame() { return next; },
    /** The epoch stamped on blocks rendered now (see resync). */
    get epoch() { return epoch; },
    now,
    event(e) { insert(e); pump(); },
    sequence(startFrame, loopFrames, events) {
      for (let i = queue.length - 1; i >= 0; i--) if (queue[i].seq) queue.splice(i, 1);
      const sorted = events.map((e) => ({ ...e, seq: true })).sort((a, b) => a.offset - b.offset);
      const start = Math.round(startFrame);
      seq = { startFrame: start, loopFrames: Math.round(loopFrames), events: sorted, loop: 0 };
      if (engine.idle() && next < start) next = start;
      pump();
    },
    stopSequence() {
      seq = null;
      for (let i = queue.length - 1; i >= 0; i--) if (queue[i].seq) queue.splice(i, 1);
    },
    /** Replace the score: events [{frame, ...}] with absolute context frames, played once. */
    score(events) {
      for (let i = queue.length - 1; i >= 0; i--) if (queue[i].score) queue.splice(i, 1);
      const sorted = [];
      for (const e of events) if (Number.isFinite(e.frame)) { e.frame = Math.round(e.frame); e.score = true; sorted.push(e); }
      sorted.sort((a, b) => a.frame - b.frame);
      score = sorted.length ? { events: sorted, at: 0 } : null;
      pump();
    },
    stopScore() {
      score = null;
      for (let i = queue.length - 1; i >= 0; i--) if (queue[i].score) queue.splice(i, 1);
    },
    /**
     * Render the next samples for the clock as it is now, stamped with `newEpoch`: the mixer was
     * told to flush this lane with the same epoch, so blocks rendered before this call are refused.
     */
    resync(newEpoch) {
      epoch = Number.isFinite(newEpoch) ? newEpoch : epoch + 1;
      open = null; // rendered for the old epoch: the mixer refuses it
      next = -1;
      sentEnd = -1;
      hole = false;
      pump();
    },
    pump,
  };
}
