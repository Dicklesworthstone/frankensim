// Playback: one Web Audio player per comparison. Both sides of an A/B pair play in lock-step through their own gain,
// so switching between the simulation and the recording keeps the playhead and only crossfades (20 ms).
import { DATA } from './app.js';

let ctx = null;
let current = null;
const decoded = new Map();

export function audioContext() {
  if (!ctx) ctx = new (window.AudioContext || window.webkitAudioContext)();
  return ctx;
}

/** Decode a site audio file once (FLAC, lossless). */
export function loadBuffer(rel) {
  if (!decoded.has(rel)) {
    decoded.set(rel, fetch(DATA + rel)
      .then(r => {
        if (!r.ok) throw new Error(`${rel}: HTTP ${r.status}`);
        return r.arrayBuffer();
      })
      .then(bytes => audioContext().decodeAudioData(bytes)));
  }
  return decoded.get(rel);
}

export function stopAll() {
  if (current) current.pause();
  current = null;
}

export class Player {
  /** `sources`: [{rel, label}], one (a single clip) or two (simulation, recording). */
  constructor(sources) {
    this.sources = sources;
    this.buffers = [];
    this.active = 0;
    this.offset = 0;
    this.playing = false;
    this.loop = false;
    this.nodes = [];
    this.listeners = new Set();
    this.ready = Promise.all(sources.map(src => loadBuffer(src.rel))).then(b => {
      this.buffers = b;
      const ac = audioContext();
      this.gains = b.map((_, i) => {
        const g = ac.createGain();
        g.gain.value = i === this.active ? 1 : 0;
        g.connect(ac.destination);
        return g;
      });
      this.emit();
      return this;
    });
  }

  get duration() {
    return this.buffers.reduce((m, b) => Math.max(m, b.duration), 0);
  }

  get position() {
    if (!this.playing) return this.offset;
    return Math.min(Math.max(audioContext().currentTime - this.t0, 0), this.duration);
  }

  on(fn) {
    this.listeners.add(fn);
    return () => this.listeners.delete(fn);
  }

  emit() {
    for (const fn of this.listeners) fn(this);
  }

  async play() {
    await this.ready;
    if (this.playing) return;
    if (current && current !== this) current.pause();
    current = this;
    const ac = audioContext();
    await ac.resume();
    if (this.offset >= this.duration - 0.01) this.offset = 0;
    const start = ac.currentTime + 0.03;
    this.t0 = start - this.offset;
    this.nodes = this.buffers.map((b, i) => {
      const src = ac.createBufferSource();
      src.buffer = b;
      src.connect(this.gains[i]);
      src.start(start, Math.min(this.offset, b.duration));
      return src;
    });
    const longest = this.nodes[this.buffers.findIndex(b => b.duration === this.duration)];
    const mine = this.nodes;
    longest.onended = () => {
      if (this.nodes !== mine || !this.playing) return; // stopped or restarted on purpose
      this.playing = false;
      this.offset = 0;
      if (this.loop) this.play();
      else this.emit();
    };
    this.playing = true;
    this.tick();
    this.emit();
  }

  pause() {
    if (!this.playing) return;
    this.offset = this.position;
    this.playing = false;
    for (const n of this.nodes) {
      try { n.stop(); } catch (e) { /* already stopped */ }
    }
    this.nodes = [];
    this.emit();
  }

  toggle() {
    return this.playing ? this.pause() : this.play();
  }

  seek(t) {
    const wasPlaying = this.playing;
    this.pause();
    this.offset = Math.min(Math.max(t, 0), this.duration);
    if (wasPlaying) this.play();
    else this.emit();
  }

  setActive(i) {
    if (i === this.active || i >= this.sources.length) return;
    this.active = i;
    if (this.gains) {
      const t = audioContext().currentTime;
      this.gains.forEach((g, k) => {
        g.gain.cancelScheduledValues(t);
        g.gain.setValueAtTime(g.gain.value, t);
        g.gain.linearRampToValueAtTime(k === i ? 1 : 0, t + 0.02);
      });
    }
    this.emit();
  }

  tick() {
    if (!this.playing) return;
    this.emit();
    requestAnimationFrame(() => this.tick());
  }
}

/** Transport controls for a player: play/pause, the A/B switch (when there are two sources), time and loop. */
export function controls(player, { names } = {}) {
  const el = (tag, cls, text) => {
    const e = document.createElement(tag);
    if (cls) e.className = cls;
    if (text) e.textContent = text;
    return e;
  };
  const bar = el('div', 'transport');
  const play = el('button', 'play', '▶');
  play.title = 'Play or pause (space)';
  play.addEventListener('click', () => player.toggle());
  bar.append(play);
  let switches = [];
  if (player.sources.length > 1) {
    const seg = el('div', 'segmented');
    switches = player.sources.map((src, i) => {
      const b = el('button', '', (names || [])[i] || src.label);
      b.title = `Listen to ${src.label} (key ${i + 1})`;
      b.addEventListener('click', () => player.setActive(i));
      seg.append(b);
      return b;
    });
    bar.append(seg);
  }
  const time = el('span', 'time', '0:00.0');
  bar.append(time);
  const loop = el('label', 'loop');
  const box = el('input');
  box.type = 'checkbox';
  box.addEventListener('change', () => { player.loop = box.checked; });
  loop.append(box, ' loop');
  bar.append(loop);
  const status = el('span', 'status', 'decoding…');
  bar.append(status);
  player.ready.then(() => { status.textContent = ''; }, err => { status.textContent = `audio unavailable: ${err.message}`; });
  player.on(p => {
    play.textContent = p.playing ? '❚❚' : '▶';
    time.textContent = `${fmtTime(p.position)} / ${fmtTime(p.duration)}`;
    switches.forEach((b, i) => b.classList.toggle('on', i === p.active));
  });
  switches.forEach((b, i) => b.classList.toggle('on', i === player.active));
  return bar;
}

function fmtTime(t) {
  if (!Number.isFinite(t)) return '0:00.0';
  const m = Math.floor(t / 60);
  return `${m}:${(t - 60 * m).toFixed(1).padStart(4, '0')}`;
}

// Keyboard: space plays or pauses the current player, 1 and 2 pick its side.
document.addEventListener('keydown', e => {
  if (!current || e.target.closest('input, textarea, select')) return;
  if (e.key === ' ') {
    e.preventDefault();
    current.toggle();
  } else if (e.key === '1' || e.key === '2') {
    current.setActive(Number(e.key) - 1);
  }
});

/** Make `p` the player the keyboard drives even before it has played. */
export function focusPlayer(p) {
  if (current && current !== p) current.pause();
  current = p;
}
