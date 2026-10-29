// Playback: one Web Audio player per clip or comparison. Both sides of an A/B pair play in lock-step through their own
// gain, so switching between the simulation and the recording keeps the playhead and only crossfades (20 ms).
// Audio is fetched and decoded on first play (or an explicit load()), so a page of thirty phrases costs nothing until
// one is played. One player sounds at a time; the "now playing" bar follows it.
import { DATA } from './app.js';

let ctx = null;
let current = null;
const decoded = new Map();
const watchers = new Set();

export function audioContext() {
  if (!ctx) {
    // iOS mutes Web Audio under the ringer switch unless the page's audio session is declared as playback.
    try {
      if (navigator.audioSession) navigator.audioSession.type = 'playback';
    } catch (e) { /* older browsers */ }
    ctx = new (window.AudioContext || window.webkitAudioContext)();
  }
  return ctx;
}

/** Decode a site audio file once. */
export function loadBuffer(rel) {
  if (!decoded.has(rel)) {
    decoded.set(rel, fetch(DATA + rel)
      .then(r => {
        if (!r.ok) throw new Error(`${rel}: HTTP ${r.status}`);
        return r.arrayBuffer();
      })
      .then(bytes => new Promise((ok, fail) => audioContext().decodeAudioData(bytes, ok, fail)))
      .catch(err => {
        decoded.delete(rel);
        throw err;
      }));
  }
  return decoded.get(rel);
}

// Bumped on every page change. A player from a page that has gone cannot start: a seek-and-play whose audio finished
// decoding after the reader left would otherwise begin sounding on the next page.
let generation = 0;

export function stopAll() {
  generation++;
  if (current) current.pause();
  setCurrent(null);
}

function setCurrent(p) {
  current = p;
  for (const fn of watchers) fn(p);
}

export function onCurrent(fn) {
  watchers.add(fn);
  return () => watchers.delete(fn);
}

export class Player {
  /** `sources`: [{rel, label, side?: 'sim'|'real'}], one (a single clip) or two (simulation, recording).
   *  `meta`: {title, sub} for the now-playing bar and the OS media controls. */
  constructor(sources, meta = {}) {
    this.sources = sources;
    this.meta = meta;
    this.buffers = [];
    this.active = 0;
    this.offset = 0;
    this.playing = false;
    this.loading = false;
    this.loop = false;
    this.nodes = [];
    this.listeners = new Set();
    this.knownDuration = meta.duration || 0;
    this._ready = null;
    this.generation = generation;
  }

  load() {
    if (!this._ready) {
      this.loading = true;
      this.emit();
      this._ready = Promise.all(this.sources.map(src => loadBuffer(src.rel))).then(b => {
        this.buffers = b;
        const ac = audioContext();
        this.gains = b.map((_, i) => {
          const g = ac.createGain();
          g.gain.value = i === this.active ? 1 : 0;
          g.connect(ac.destination);
          return g;
        });
        this.loading = false;
        this.emit();
        return this;
      }, err => {
        this.loading = false;
        this.error = err;
        this._ready = null;
        this.emit();
        throw err;
      });
    }
    return this._ready;
  }

  get ready() {
    return this.load();
  }

  get duration() {
    return this.buffers.reduce((m, b) => Math.max(m, b.duration), 0) || this.knownDuration;
  }

  get position() {
    if (!this.playing) return this.offset;
    return Math.min(Math.max(audioContext().currentTime - this.t0, 0), this.duration);
  }

  /** Call `fn(player)` on every change. With `owner` (a DOM node) the listener goes away once the node has been on the
   *  page and left it, so views rebuilt over a cached player do not pile up listeners. */
  on(fn, owner) {
    const entry = { fn, owner, seen: false };
    this.listeners.add(entry);
    return () => this.listeners.delete(entry);
  }

  emit() {
    for (const l of this.listeners) {
      if (l.owner) {
        if (l.owner.isConnected) l.seen = true;
        else if (l.seen) {
          this.listeners.delete(l);
          continue;
        }
      }
      l.fn(this);
    }
  }

  /** Start playing (after loading, the first time). Never rejects: a load failure shows as `error`. */
  async play() {
    if (this.generation !== generation) return;
    if (current && current !== this) current.pause();
    setCurrent(this);
    const want = {};
    this.wanted = want; // a pause or another play while this one loads cancels it
    const ac = audioContext();
    const resumed = ac.resume();
    try {
      await this.load();
      await resumed;
    } catch (err) {
      if (this.wanted === want) this.wanted = null;
      return;
    }
    if (this.playing || current !== this || this.wanted !== want) return;
    this.wanted = null;
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
    if ('mediaSession' in navigator && window.MediaMetadata) {
      navigator.mediaSession.metadata = new MediaMetadata({ title: this.meta.title || this.sources[this.active].label, artist: this.meta.sub || 'FrankenSim instrument tests' });
    }
    this.tick();
    this.emit();
  }

  pause() {
    if (this.wanted) {
      this.wanted = null;
      this.emit();
    }
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
    return this.playing || this.wanted ? this.pause() : this.play();
  }

  seek(t) {
    const wasPlaying = this.playing;
    const pending = this.wanted; // a seek during loading keeps the play that is waiting for it
    this.pause();
    this.wanted = pending;
    this.offset = Math.min(Math.max(t, 0), this.duration || t);
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

const ICON = {
  play: '<svg viewBox="0 0 24 24" aria-hidden="true"><path d="M8 5.5v13a1 1 0 0 0 1.5.86l10.6-6.5a1 1 0 0 0 0-1.72L9.5 4.64A1 1 0 0 0 8 5.5Z" fill="currentColor"/></svg>',
  pause: '<svg viewBox="0 0 24 24" aria-hidden="true"><rect x="6" y="5" width="4.2" height="14" rx="1.2" fill="currentColor"/><rect x="13.8" y="5" width="4.2" height="14" rx="1.2" fill="currentColor"/></svg>',
  close: '<svg viewBox="0 0 24 24" aria-hidden="true"><path d="M6 6l12 12M18 6 6 18" stroke="currentColor" stroke-width="2.2" stroke-linecap="round"/></svg>',
};

/** A round play/pause button bound to `player`. */
export function playButton(player, { small = false, label } = {}) {
  const b = document.createElement('button');
  b.type = 'button';
  b.className = `play${small ? ' small' : ''}`;
  const name = label || player.meta.title || 'clip';
  const draw = p => {
    b.innerHTML = p.playing ? ICON.pause : ICON.play;
    b.setAttribute('aria-label', `${p.playing ? 'Pause' : 'Play'} ${name}`);
    b.classList.toggle('busy', p.loading);
  };
  b.addEventListener('click', e => {
    e.stopPropagation();
    player.toggle();
  });
  player.on(draw, b);
  draw(player);
  return b;
}

/** Transport controls: play/pause, the A/B switch (two sources), a seek bar, the time and loop. */
export function controls(player, { names, keys = true } = {}) {
  const el = (tag, cls, text) => {
    const e = document.createElement(tag);
    if (cls) e.className = cls;
    if (text) e.textContent = text;
    return e;
  };
  const bar = el('div', 'transport');
  bar._player = player;
  if (watchTransports) watchTransports.observe(bar);
  bar.append(playButton(player));
  let switches = [];
  if (player.sources.length > 1) {
    const seg = el('div', 'segmented');
    seg.setAttribute('role', 'group');
    seg.setAttribute('aria-label', 'Which to hear');
    switches = player.sources.map((src, i) => {
      const b = el('button', src.side || '', (names || [])[i] || src.label);
      b.type = 'button';
      b.title = `Hear ${src.label} (key ${i + 1})`;
      b.addEventListener('click', () => player.setActive(i));
      seg.append(b);
      return b;
    });
    bar.append(seg);
  }
  const scrub = el('div', 'scrub');
  const fill = el('span');
  scrub.append(fill);
  scrub.tabIndex = 0;
  scrub.setAttribute('role', 'slider');
  scrub.setAttribute('aria-label', 'Seek');
  scrub.setAttribute('aria-valuemin', '0');
  scrub.addEventListener('click', e => {
    const r = scrub.getBoundingClientRect();
    player.load().then(() => player.seek(((e.clientX - r.left) / r.width) * player.duration), () => {});
  });
  scrub.addEventListener('keydown', e => {
    const step = e.shiftKey ? 5 : 1;
    if (e.key === 'ArrowRight') { player.seek(player.position + step); e.preventDefault(); }
    if (e.key === 'ArrowLeft') { player.seek(player.position - step); e.preventDefault(); }
  });
  bar.append(scrub);
  const time = el('span', 'time', '0:00.0');
  bar.append(time);
  const loop = el('label', 'loop');
  const box = el('input');
  box.type = 'checkbox';
  box.addEventListener('change', () => { player.loop = box.checked; });
  loop.append(box, 'loop');
  bar.append(loop);
  const status = el('span', 'status');
  status.setAttribute('aria-live', 'polite');
  bar.append(status);
  if (keys && player.sources.length > 1) {
    const k = el('span', 'keys');
    k.innerHTML = '<kbd>space</kbd> play · <kbd>1</kbd> <kbd>2</kbd> switch · <kbd>←</kbd> <kbd>→</kbd> seek';
    bar.append(k);
  }
  const draw = p => {
    const d = p.duration || 0;
    time.textContent = `${fmtTime(p.position)} / ${fmtTime(d)}`;
    fill.style.width = d ? `${(100 * p.position) / d}%` : '0';
    scrub.className = `scrub${p.sources[p.active] && p.sources[p.active].side === 'real' ? ' real' : ''}`;
    scrub.setAttribute('aria-valuemax', d.toFixed(1));
    scrub.setAttribute('aria-valuenow', p.position.toFixed(1));
    switches.forEach((b, i) => {
      b.classList.toggle('on', i === p.active);
      b.setAttribute('aria-pressed', String(i === p.active));
    });
    status.textContent = p.error ? `audio unavailable: ${p.error.message}` : p.loading ? 'loading…' : '';
  };
  player.on(draw, bar);
  draw(player);
  return bar;
}

function fmtTime(t) {
  if (!Number.isFinite(t) || t < 0) return '0:00.0';
  const tenths = Math.round(t * 10); // round first, so 59.97 s reads 1:00.0, not 0:60.0
  const m = Math.floor(tenths / 600);
  return `${m}:${((tenths - 600 * m) / 10).toFixed(1).padStart(4, '0')}`;
}

// Keyboard: space plays or pauses the current player, 1 and 2 pick its side, arrows seek.
document.addEventListener('keydown', e => {
  if (!current || e.target.closest('input:not([type=checkbox]), textarea, select') || e.metaKey || e.ctrlKey || e.altKey) return;
  // 1 and 2 work from anywhere (right after pressing Play the focus is still on the button); space and the arrows
  // keep their own meaning on controls that use them.
  if (e.key !== '1' && e.key !== '2' && e.target.closest('input, button, a, [role="slider"], [tabindex="0"]')) return;
  if (e.key === ' ') {
    e.preventDefault();
    current.toggle();
  } else if (e.key === '1' || e.key === '2') {
    current.setActive(Number(e.key) - 1);
  } else if (e.key === 'ArrowRight' || e.key === 'ArrowLeft') {
    e.preventDefault();
    current.seek(current.position + (e.key === 'ArrowRight' ? 1 : -1) * (e.shiftKey ? 5 : 1));
  }
});

/** Make `p` the player the keyboard drives even before it has played. */
export function focusPlayer(p) {
  if (current && current !== p) current.pause();
  setCurrent(p);
}

// The transports on screen. The floating bar stays away while the playing sound's own controls are in view (a phrase
// page's sticky player, say), where it would only repeat them and cover the page.
const transportsInView = new Set();
let refreshNowPlaying = () => {};
const watchTransports = typeof IntersectionObserver === 'function' ? new IntersectionObserver(entries => {
  for (const e of entries) {
    if (e.isIntersecting) transportsInView.add(e.target);
    else transportsInView.delete(e.target);
  }
  refreshNowPlaying();
}) : null;
const ownControlsInView = p => [...transportsInView].some(b => b.isConnected && b._player === p);

/** The floating "now playing" bar: what sounds, which side, pause and stop. */
export function mountNowPlaying(el) {
  let off = null;
  const title = document.createElement('div');
  title.className = 'np-title';
  const sub = document.createElement('div');
  sub.className = 'np-sub';
  const text = document.createElement('div');
  text.className = 'np-text';
  text.append(title, sub);
  const side = document.createElement('span');
  const barFill = document.createElement('span');
  const bar = document.createElement('div');
  bar.className = 'np-bar';
  bar.append(barFill);
  const close = document.createElement('button');
  close.className = 'iconbtn';
  close.type = 'button';
  close.innerHTML = ICON.close;
  close.setAttribute('aria-label', 'Stop and close');
  close.addEventListener('click', () => { if (current) { current.pause(); current.seek(0); } el.hidden = true; });
  let shown = null;
  onCurrent(p => {
    if (off) off();
    off = null;
    if (!p) { el.hidden = true; return; }
    const draw = q => {
      if (!q.playing && !q.position) return;
      el.hidden = ownControlsInView(q);
      if (shown !== q) {
        shown = q;
        el.replaceChildren(playButton(q, { small: true }), text, side, close, bar);
      }
      title.textContent = q.meta.title || q.sources[q.active].label;
      sub.textContent = q.meta.sub || '';
      const s = q.sources[q.active];
      side.className = `np-side ${s.side || 'sim'}`;
      side.textContent = q.sources.length > 1 ? s.label : (s.side === 'real' ? 'recording' : 'simulation');
      barFill.style.width = q.duration ? `${(100 * q.position) / q.duration}%` : '0';
      barFill.style.background = s.side === 'real' ? 'var(--real)' : 'var(--sim)';
    };
    off = p.on(draw);
    refreshNowPlaying = () => draw(p);
    draw(p);
  });
}
