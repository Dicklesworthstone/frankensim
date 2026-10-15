// The site's LIVE PIANO section: the reduced real-time concert grand running
// in WebAssembly, played from an 88-key strip, a computer keyboard or MIDI,
// with the strings and the soundboard drawn from the engine's own state, the
// live energy ledger, and an honest account of what the live model omits.
//
// Nothing is fetched and no audio starts until the visitor presses
// "Start the live piano".

import './live-piano.css';
import { resumeContext } from '../audio';
import { LivePianoClient, midiVelocityToMS, type PianoInfo, type PianoSnapshot, type PianoStats } from './pianoClient';
import { LiveRefusal } from './liveAudio';
import { PianoStage, noteName, type PianoView, type StageScales } from './pianoStage';

const FIRST = 21;
const LAST = 108;
const BLACK = new Set([1, 3, 6, 8, 10]);
/** Black-key offsets from the white-key boundary, in black-key widths (as on a real keyboard). */
const BLACK_OFFSET: Record<number, number> = { 1: -0.12, 3: 0.12, 6: -0.16, 8: 0, 10: 0.16 };
/** Computer keys -> semitones above the current C (A = C, W = C sharp, ...). */
const KEY_OFFSETS: Record<string, number> = {
  KeyA: 0, KeyW: 1, KeyS: 2, KeyE: 3, KeyD: 4, KeyF: 5, KeyT: 6, KeyG: 7, KeyY: 8, KeyH: 9,
  KeyU: 10, KeyJ: 11, KeyK: 12, KeyO: 13, KeyL: 14, KeyP: 15, Semicolon: 16,
};
const SPEED_MIN = 0.5;
const SPEED_MAX = 4.5;

interface DemoStep { at: number; run: () => void }
interface Demo { label: string; text: string; build: (s: LivePianoSection) => DemoStep[] }

const DEMOS: Record<string, Demo> = {
  dynamics: {
    label: 'pp → ff on one key',
    text: 'A3 three times, at 0.6, 2.0 and 4.5 m/s. The wool felt is a stiffening power law, so a faster hammer drives it into its stiff range: the contact is shorter and stronger and excites more of the upper partials. Watch the felt force on the hammer, listen for the brighter attack.',
    build: (s) => [
      { at: 0, run: () => s.press(57, 0.6, 'demo') },
      { at: 1300, run: () => s.release(57, 'demo') },
      { at: 1500, run: () => s.press(57, 2.0, 'demo') },
      { at: 2800, run: () => s.release(57, 'demo') },
      { at: 3000, run: () => s.press(57, 4.5, 'demo') },
      { at: 4700, run: () => s.release(57, 'demo') },
    ],
  },
  octave: {
    label: 'Octave with pedal',
    text: 'Pedal down, then C3 and C4. The keys come up at once, but with the dampers lifted both notes keep ringing, coupled through the one bridge and soundboard they share. After five seconds the pedal comes up.',
    build: (s) => [
      { at: 0, run: () => s.setPedal(1, 'demo') },
      { at: 60, run: () => s.press(48, 2.2, 'demo') },
      { at: 420, run: () => s.press(60, 2.2, 'demo') },
      { at: 700, run: () => { s.release(48, 'demo'); s.release(60, 'demo'); } },
      { at: 5600, run: () => s.setPedal(0, 'demo') },
    ],
  },
  chord: {
    label: 'A chord, then lift the pedal',
    text: 'Pedal down and a C major chord across two octaves, released immediately: the pedal holds it. After three seconds the pedal lifts and the dampers fall on every string at once. Watch the damper share of the dissipated energy in the ledger.',
    build: (s) => {
      const chord = [48, 52, 55, 60, 64, 67];
      return [
        { at: 0, run: () => s.setPedal(1, 'demo') },
        { at: 80, run: () => chord.forEach((m) => s.press(m, 2.4, 'demo')) },
        { at: 420, run: () => chord.forEach((m) => s.release(m, 'demo')) },
        { at: 3400, run: () => s.setPedal(0, 'demo') },
      ];
    },
  },
  lowest: {
    label: 'Lowest string (A0)',
    text: 'A0: one wound string, 2.01 m long. Its partials follow from its tension, bending stiffness and mass, so they stretch progressively above whole multiples of 27.5 Hz. This live model keeps its first 200 of them.',
    build: (s) => [
      { at: 0, run: () => s.press(21, 2.8, 'demo') },
      { at: 5200, run: () => s.release(21, 'demo') },
    ],
  },
};

let serial = 0;

function el<K extends keyof HTMLElementTagNameMap>(tag: K, className?: string, text?: string): HTMLElementTagNameMap[K] {
  const node = document.createElement(tag);
  if (className) node.className = className;
  if (text !== undefined) node.textContent = text;
  return node;
}

function q<T extends Element>(root: ParentNode, selector: string): T {
  const found = root.querySelector<T>(selector);
  if (!found) throw new Error(`live piano: missing ${selector}`);
  return found;
}

function setText(node: Element, text: string): void {
  if (node.textContent !== text) node.textContent = text;
}

function millijoules(value: number): string {
  const mj = value * 1e3;
  if (mj !== 0 && Math.abs(mj) < 0.01) return `${(value * 1e6).toFixed(Math.abs(value) < 1e-7 ? 3 : 1)} µJ`;
  if (Math.abs(mj) >= 100) return `${mj.toFixed(1)} mJ`;
  if (Math.abs(mj) >= 1) return `${mj.toFixed(2)} mJ`;
  return `${mj.toFixed(3)} mJ`;
}

function scientific(value: number): string {
  if (value === 0) return '0 J';
  const [mantissa, exponent] = value.toExponential(1).split('e');
  return `${mantissa} × 10${(exponent ?? '0').replace('+', '').replace('-', '⁻').replace(/\d/g, (d) => '⁰¹²³⁴⁵⁶⁷⁸⁹'[Number(d)] ?? d)} J`;
}

function isTyping(target: EventTarget | null): boolean {
  if (!(target instanceof HTMLElement)) return false;
  if (target.isContentEditable) return true;
  const tag = target.tagName;
  if (tag === 'TEXTAREA' || tag === 'SELECT') return true;
  if (tag === 'INPUT') {
    const type = (target as HTMLInputElement).type;
    return !['range', 'checkbox', 'radio', 'button'].includes(type);
  }
  return false;
}

export class LivePianoSection {
  private readonly root: HTMLElement;
  private readonly id = `lp${++serial}`;
  private readonly stage: PianoStage;
  private readonly stageEl: HTMLElement;
  private readonly kbd: HTMLElement;
  private readonly scroll: HTMLElement;
  private readonly keyEls = new Map<number, HTMLElement>();
  private client: LivePianoClient | null = null;
  private starting = false;
  private disposed = false;
  private visible = false;
  private raf = 0;
  private lastFrame = 0;
  private queue: Array<{ due: number; snapshot: PianoSnapshot }> = [];
  private lastShown: PianoSnapshot | null = null;
  private readonly down = new Map<number, number>();
  private readonly pointerKeys = new Map<number, number>();
  private readonly codeKeys = new Map<string, number>();
  private readonly demoKeys = new Set<number>();
  private touchMode: 'position' | 'fixed' = 'position';
  private fixedSpeed = 2.5;
  private octave = 4;
  private pedal = 0;
  private pedalBeforeSpace = 0;
  private demoTimers: number[] = [];
  private demoName: string | null = null;
  private lastStats: PianoStats | null = null;
  private sounding = new Set<number>();
  private lastAnnounce = 0;
  private midiAccess: MIDIAccess | null = null;
  private readonly reduced = window.matchMedia('(prefers-reduced-motion: reduce)');
  private readonly cleanups: Array<() => void> = [];

  constructor(root: HTMLElement) {
    this.root = root;
    root.classList.add('lp', 'live-scope');
    root.innerHTML = this.skeleton();
    this.stageEl = q<HTMLElement>(root, '.lp-stage');
    this.kbd = q<HTMLElement>(root, '.lp-kbd');
    this.scroll = q<HTMLElement>(root, '.lp-kbd-scroll');
    this.stage = new PianoStage(this.stageEl);
    this.stage.reducedMotion = this.reduced.matches;
    this.stage.onScales = (scales) => this.paintScales(scales);
    this.buildKeyboard();
    this.buildDemos();
    this.buildHonesty(null);
    this.wire();
    this.layoutKeyboard(true);
    this.stage.draw(0);
    this.syncView();
  }

  // ------------------------------------------------------------------ markup
  private skeleton(): string {
    const id = this.id;
    return `
<div class="lp-instrument">
  <div class="stage lp-stage" role="group" aria-label="Live piano: sounding strings and soundboard" aria-describedby="${id}-stage-help">
    <div class="stage-top">
      <span class="tag live">live</span>
      <div class="seg live lp-seg-sm lp-view" role="group" aria-label="View">
        <button type="button" data-view="strings" aria-pressed="true">Strings</button>
        <button type="button" data-view="board" aria-pressed="false">Soundboard</button>
      </div>
    </div>
    <div class="stage-hud"><span class="lp-hud-strings"></span><span class="lp-hud-board"></span></div>
    <div class="stage-note lp-note">
      <div class="lp-note-body">
        <span class="tag live">Live engine</span>
        <p class="lp-note-title">A concert grand solved in your browser while you play it.</p>
        <p class="small lp-note-text">Starting downloads FrankenSim's reduced piano engine as WebAssembly (about 700&nbsp;kB
          compressed, the Model&nbsp;D soundboard included) and assembles all 88 keys from the published string scale.
          The piano starts silent and at rest.</p>
        <button type="button" class="btn live lp-start">Start the live piano</button>
        <p class="small lp-note-status" role="status"></p>
      </div>
    </div>
  </div>
  <p id="${id}-stage-help" class="sr-only">The picture shows, for every sounding key, the real transverse shape of its first
    string between the agraffe and the bridge, its hammer and felt force, and the soundboard's signed displacement.</p>

  <div class="lp-deck">
    <div class="lp-controls">
      <div class="field lp-touch">
        <span class="label"><span>Hammer speed</span><output class="lp-speed-out num">from key position</output></span>
        <div class="lp-row">
          <div class="seg live lp-seg-sm lp-touch-mode" role="group" aria-label="Hammer speed source">
            <button type="button" data-touch="position" aria-pressed="true">Where you press</button>
            <button type="button" data-touch="fixed" aria-pressed="false">Fixed</button>
          </div>
          <input type="range" class="lp-speed" min="0.3" max="4.5" step="0.1" value="2.5" aria-label="Fixed hammer speed in metres per second">
        </div>
      </div>
      <div class="field lp-pedal">
        <span class="label"><span>Sustain pedal</span><output class="lp-pedal-out num">up</output></span>
        <div class="lp-row">
          <button type="button" class="btn sm lp-pedal-btn" aria-pressed="false">Pedal</button>
          <input type="range" class="lp-pedal-range" min="0" max="1" step="0.01" value="0" aria-label="Sustain pedal travel, 0 to 1">
        </div>
      </div>
      <div class="field lp-midi" hidden>
        <span class="label"><span>MIDI</span><output class="lp-midi-out">not connected</output></span>
        <div class="lp-row"><button type="button" class="btn sm lp-midi-btn">Connect MIDI</button></div>
      </div>
    </div>
    <div class="lp-scrub" aria-hidden="true"><div class="lp-scrub-keys"></div><div class="lp-scrub-window"></div></div>
    <div class="kbd-scroll lp-kbd-scroll">
      <div class="kbd live tall lp-kbd" role="group" aria-label="Piano keyboard, A0 to C8" aria-describedby="${id}-legend"></div>
    </div>
    <p id="${id}-legend" class="small muted lp-legend">
      <span class="lp-legend-keys"><kbd>A</kbd><kbd>W</kbd><kbd>S</kbd><kbd>E</kbd><kbd>D</kbd><kbd>F</kbd><kbd>T</kbd><kbd>G</kbd><kbd>Y</kbd><kbd>H</kbd><kbd>U</kbd><kbd>J</kbd><kbd>K</kbd><kbd>O</kbd><kbd>L</kbd><kbd>P</kbd><kbd>;</kbd>
        play <b class="lp-range">C4–E5</b></span>
      <span><kbd>Z</kbd> <kbd>X</kbd> octave</span>
      <span><kbd>Space</kbd> pedal (hold)</span>
      <span class="lp-touch-hint">Press near the front edge of a key for a faster hammer; slide for a glissando.</span>
    </p>
  </div>
</div>

<div class="grid wide-left lp-grid">
  <div class="stack">
    <div class="panel lp-try" aria-labelledby="${id}-try">
      <div class="panel-head"><h3 id="${id}-try">Try this</h3><button type="button" class="btn sm ghost lp-demo-stop" hidden>Stop</button></div>
      <div class="lp-demos"></div>
      <p class="small lp-demo-text" role="status">Each one only presses keys and the pedal on a timeline; everything you hear is computed.</p>
    </div>
    <div class="panel lp-ledger-panel" aria-labelledby="${id}-ledger">
      <div class="panel-head"><h3 id="${id}-ledger">Energy ledger</h3><span class="tag live">live</span></div>
      <div class="ledger lp-ledger">
        <div class="ledger-bar" role="img" aria-label="Share of the supplied energy now stored, dissipated and retired">
          <i class="s"></i><i class="d"></i><i class="r"></i>
        </div>
        <div class="ledger-keys">
          <span><i style="background:var(--supplied)"></i>supplied <b class="num" data-k="input">0 mJ</b></span>
          <span><i style="background:var(--stored)"></i>stored <b class="num" data-k="stored">0 mJ</b></span>
          <span><i style="background:var(--lost)"></i>dissipated <b class="num" data-k="dissipated">0 mJ</b></span>
          <span><i style="background:var(--supplied);opacity:.55"></i>retired <b class="num" data-k="retired">0 mJ</b></span>
        </div>
        <p class="lp-breakdown small num muted" data-k="breakdown"></p>
        <p class="lp-residual"><span class="muted">balance residual</span> <b class="num" data-k="residual">0 J</b>
          <span class="muted small">supplied − stored − dissipated − retired</span></p>
      </div>
      <p class="small muted lp-ledger-note">Supplied is the work the hammers bring in. It is either still stored (strings,
        board, hammers, felt), dissipated (string and board losses, dampers, felt and shank, the hammer caught on its
        rest), or retired: relative string energy discarded when a quiet key returns to its resting image, or when the
        polyphony budget forces the quietest key out. The residual is what the arithmetic fails to close.</p>
    </div>
  </div>
  <div class="stack lp-side">
    <dl class="stats lp-stats">
      <div class="stat"><dt>Sounding keys</dt><dd data-s="voices">–</dd></div>
      <div class="stat"><dt>Strings · modes</dt><dd data-s="modes">–</dd></div>
      <div class="stat"><dt>Real time</dt><dd data-s="rtf">–</dd></div>
      <div class="stat"><dt>Latency</dt><dd data-s="latency">–</dd></div>
      <div class="stat"><dt>Underruns</dt><dd data-s="underruns">–</dd></div>
      <div class="stat"><dt>Peak pressure</dt><dd data-s="peak">–</dd></div>
      <div class="stat"><dt>Clipped</dt><dd data-s="clipped">–</dd></div>
      <div class="stat"><dt>Forced out</dt><dd data-s="forced">–</dd></div>
      <div class="stat lp-stat-wide"><dt>Quality granted</dt><dd data-s="quality">–</dd></div>
    </dl>
    <p class="small muted lp-stats-note">Latency is the audio queued ahead of the playhead, plus the sound's travel to
      the microphone 1&nbsp;m above the board, plus your output device. The hammer's 2&nbsp;mm flight to the string adds
      0.4 to 4&nbsp;ms. Peak pressure is at that microphone; 10&nbsp;Pa is digital full scale.</p>
    <div class="refusal lp-refusal" hidden><p></p><pre></pre></div>
  </div>
</div>
<div class="panel lp-honesty"></div>
<div class="sr-only" aria-live="polite" data-announce></div>`;
  }

  private buildKeyboard(): void {
    const scrubKeys = q<HTMLElement>(this.root, '.lp-scrub-keys');
    for (let midi = FIRST; midi <= LAST; midi++) {
      const black = BLACK.has(midi % 12);
      const key = el('div', black ? 'k' : 'w');
      key.dataset.midi = String(midi);
      if (!black && midi % 12 === 0) key.append(el('span', 'lp-oct', noteName(midi)));
      this.kbd.append(key);
      this.keyEls.set(midi, key);
      if (!black) scrubKeys.append(el('i', midi % 12 === 0 ? 'c' : undefined));
    }
  }

  private layoutKeyboard(centre: boolean): void {
    const avail = this.scroll.clientWidth || 320;
    const coarse = window.matchMedia('(pointer: coarse)').matches;
    const ww = Math.max(coarse ? 22 : 15, avail / 52);
    const bw = ww * 0.6;
    let white = 0;
    for (let midi = FIRST; midi <= LAST; midi++) {
      const key = this.keyEls.get(midi);
      if (!key) continue;
      const pitch = midi % 12;
      if (BLACK.has(pitch)) {
        key.style.left = `${(white * ww - bw / 2 + (BLACK_OFFSET[pitch] ?? 0) * bw).toFixed(2)}px`;
        key.style.width = `${bw.toFixed(2)}px`;
      } else {
        key.style.left = `${(white * ww).toFixed(2)}px`;
        key.style.width = `${(ww + 0.5).toFixed(2)}px`;
        white++;
      }
    }
    this.kbd.style.width = `${(52 * ww).toFixed(1)}px`;
    const overflow = 52 * ww > avail + 1;
    q<HTMLElement>(this.root, '.lp-scrub').hidden = !overflow;
    this.root.classList.toggle('lp-narrow', overflow);
    if (centre && overflow) {
      const c4 = this.keyEls.get(60);
      if (c4) this.scroll.scrollLeft = c4.offsetLeft + ww / 2 - avail / 2;
    }
    this.paintScrub();
  }

  private paintScrub(): void {
    const total = this.kbd.offsetWidth || 1;
    const avail = this.scroll.clientWidth;
    const win = q<HTMLElement>(this.root, '.lp-scrub-window');
    win.style.left = `${((this.scroll.scrollLeft / total) * 100).toFixed(2)}%`;
    win.style.width = `${Math.min(100, (avail / total) * 100).toFixed(2)}%`;
  }

  private buildDemos(): void {
    const host = q<HTMLElement>(this.root, '.lp-demos');
    for (const [name, demo] of Object.entries(DEMOS)) {
      const button = el('button', 'btn sm lp-demo', demo.label);
      button.type = 'button';
      button.dataset.demo = name;
      button.disabled = true;
      host.append(button);
    }
  }

  private buildHonesty(info: PianoInfo | null): void {
    const host = q<HTMLElement>(this.root, '.lp-honesty');
    host.replaceChildren();
    const details = el('details', 'evidence');
    details.append(el('summary', undefined, 'What this live model is — and is not'));
    const body = el('div', 'evidence-body lp-evidence');
    body.append(el('p', undefined, 'The reduced real-time model. It steps the same equations, with the same time discretisation, as FrankenSim\'s offline grand-piano reference, compiled from the reference\'s own source files: stiff strings whose partials follow from tension, bending stiffness, mass and length (never assigned pitches), a moving bridge coupled reciprocally to the soundboard, a felt hammer with crush memory and relaxation solved implicitly against the strings, and sound pressure from the Rayleigh integral over the moving board. What it reduces is listed here, with the measured cost of each reduction.'));
    const lists = el('div', 'lp-incl');
    const column = (title: string, items: string[], cls: string): HTMLElement => {
      const box = el('div');
      box.append(el('h4', undefined, title));
      const list = el('ul', `lp-list ${cls}`);
      for (const item of items) list.append(el('li', undefined, item));
      box.append(list);
      return box;
    };
    lists.append(
      column('In the live model', [
        'All 88 keys of the published scale, 234 speaking strings on the bridge. Only the sounding keys are dynamic: every unison string and duplex segment, up to 200 string modes per key (low trichords keep partials to about 5–7 kHz).',
        'Every other string stays attached as a quasi-static chord (its mass and tension load the bridge); a key becomes dynamic at its strike and returns to rest when quiet.',
        'One shared soundboard: 32 coordinates of the Model D board (16 exact modes to 639 Hz plus 16 bridge-response vectors from the 0–2400 Hz slice), with its full damping matrix.',
        'The published hammer shank (rigid rotation plus bending), the per-key wool-felt card with crush memory and Prony relaxation, contact solved implicitly every step.',
        'The estimated damper pad on every string up to E6 (four quadrature points), and a continuous sustain pedal.',
        'Pressure from the time-domain Rayleigh integral at a microphone 1 m above the board, at a fixed 10 Pa full scale. No samples, no normalisation, no EQ, no reverb.',
      ], 'in'),
      column('Reduced or absent', [
        'One mechanical step per sample at 48 kHz; the offline reference uses four.',
        '32 board coordinates instead of 128: measured on A3, about 10 dB less energy above 2 kHz than the offline render.',
        'Strings of silent keys neither resonate nor absorb: offline, their dampers soak up about 10 % of an A3 strike, so the live notes decay more slowly early on (≈ −9.6 against −15.2 dB/s).',
        `A polyphony budget of ${info ? info.maxDynamicModes.toLocaleString('en-US') : '2,400'} string modes and ${info ? info.maxSoundingKeys : 24} keys. Beyond it the quietest key is returned to rest early; this is counted ("Forced out") and its energy appears as retired in the ledger. A load governor may lower the budget for new keys on a slow device, and says so.`,
        'No una corda or sostenuto, no key or action beyond the shank, no room, lid or case, no air loading of the board.',
      ], 'out'),
    );
    body.append(lists);
    body.append(el('h4', undefined, 'Checked against the offline reference'));
    body.append(el('p', undefined, 'With nothing but one key on the bridge, the live kernel and the reference engine (the same source files) agree on bridge velocity to within 2.3 × 10⁻⁹ relative over one second, for eight keys at three hammer speeds, and their energy books agree to about ten digits. For an isolated C4 on the same 32-coordinate board, the live model\'s pressure, compiled natively, equals the offline renderer\'s file in all but 137 of 48,000 float32 samples, which differ by at most 1.5 × 10⁻⁸ Pa; the WebAssembly build in this tab sums its string modes in two lanes, which moves its output by less than 10⁻¹⁰ (rms, relative). Against the full offline render, the first 15 partials of A3 sit within 2.5 cents and the fitted inharmonicity is 2.42 × 10⁻⁴ (offline 2.47 × 10⁻⁴).'));
    const chips = (title: string, cls: string, items: string[]): HTMLElement => {
      const row = el('div', 'lp-chiprow');
      row.append(el('h4', undefined, title));
      const list = el('div', 'chips');
      for (const item of items) list.append(el('span', `chip ${cls}`, item));
      row.append(list);
      return row;
    };
    body.append(
      chips('Published', 'pub', ['RT-0425 string scale, 84 keys', 'RT-0425 hammer mass, exponent, stiffness', 'Model D dimensions', 'Model D board drawing (Boutillon et al.)', 'Hammer shank geometry (Chabassier & Durufle)']),
      chips('Estimated', 'est', ['4 end keys of the scale', 'Felt crush and relaxation', 'Unison detune 0.8 cents', 'String loss law', 'Damper pads, 0.4 N s/m', 'Board damping 1.5 %', 'Bridge stations, rib sections']),
      chips('Not modelled live', 'not', ['Silent-string resonance', 'Silent-string damper loss', 'Board above 32 coordinates', 'Una corda, sostenuto', 'Action beyond the shank', 'Room, lid, case']),
    );
    if (info) {
      const prov = el('p', 'small muted');
      prov.textContent = `Engine statement: ${info.boardCoordinates} board coordinates (${info.boardExactLowModes} exact); ${info.speakingStrings} strings on the bridge; ${info.tablePartials.toLocaleString('en-US')} partials in the string tables; microphone at ${info.microphoneM.join(', ')} m. Board: ${info.boardProvenance}.`;
      body.append(prov);
    }
    details.append(body);
    host.append(details);
  }

  // ------------------------------------------------------------------ wiring
  private listen<K extends keyof WindowEventMap>(target: Window, type: K, fn: (e: WindowEventMap[K]) => void, options?: AddEventListenerOptions): void;
  private listen<K extends keyof DocumentEventMap>(target: Document, type: K, fn: (e: DocumentEventMap[K]) => void, options?: AddEventListenerOptions): void;
  private listen(target: Window | Document, type: string, fn: (e: Event) => void, options?: AddEventListenerOptions): void {
    target.addEventListener(type, fn, options);
    this.cleanups.push(() => target.removeEventListener(type, fn, options));
  }

  private wire(): void {
    q<HTMLButtonElement>(this.root, '.lp-start').addEventListener('click', () => void this.start());

    for (const button of this.root.querySelectorAll<HTMLButtonElement>('.lp-view button')) {
      button.addEventListener('click', () => {
        this.stage.view = (button.dataset.view as PianoView) ?? 'strings';
        this.syncView();
        this.stage.draw(0);
      });
    }

    // hammer speed source
    const speed = q<HTMLInputElement>(this.root, '.lp-speed');
    const paintSpeed = (): void => {
      this.fixedSpeed = Number(speed.value);
      speed.style.setProperty('--fill', `${(((this.fixedSpeed - 0.3) / 4.2) * 100).toFixed(1)}%`);
      if (this.touchMode === 'fixed') setText(q(this.root, '.lp-speed-out'), `${this.fixedSpeed.toFixed(1)} m/s`);
    };
    speed.addEventListener('input', paintSpeed);
    paintSpeed();
    for (const button of this.root.querySelectorAll<HTMLButtonElement>('.lp-touch-mode button')) {
      button.addEventListener('click', () => {
        this.touchMode = button.dataset.touch === 'fixed' ? 'fixed' : 'position';
        for (const b of this.root.querySelectorAll<HTMLButtonElement>('.lp-touch-mode button')) {
          b.setAttribute('aria-pressed', String(b === button));
        }
        setText(q(this.root, '.lp-speed-out'), this.touchMode === 'fixed' ? `${this.fixedSpeed.toFixed(1)} m/s` : 'from key position');
        this.root.classList.toggle('lp-fixed', this.touchMode === 'fixed');
      });
    }

    // pedal
    const pedalRange = q<HTMLInputElement>(this.root, '.lp-pedal-range');
    pedalRange.addEventListener('input', () => this.setPedal(Number(pedalRange.value), 'ui'));
    q<HTMLButtonElement>(this.root, '.lp-pedal-btn').addEventListener('click', () => this.setPedal(this.pedal > 0 ? 0 : 1, 'ui'));
    this.paintPedal();

    // MIDI
    if (typeof navigator.requestMIDIAccess === 'function') {
      q<HTMLElement>(this.root, '.lp-midi').hidden = false;
      q<HTMLButtonElement>(this.root, '.lp-midi-btn').addEventListener('click', () => void this.connectMidi());
    }

    // keyboard strip: pointer capture per pointer, glissando by hit-testing
    this.kbd.addEventListener('pointerdown', (e) => this.onPointerDown(e));
    this.kbd.addEventListener('pointermove', (e) => this.onPointerMove(e));
    const end = (e: PointerEvent): void => this.onPointerEnd(e);
    this.kbd.addEventListener('pointerup', end);
    this.kbd.addEventListener('pointercancel', end);
    this.kbd.addEventListener('lostpointercapture', end);
    this.kbd.addEventListener('contextmenu', (e) => e.preventDefault());
    this.scroll.addEventListener('scroll', () => this.paintScrub(), { passive: true });
    this.wireScrub();

    // demos
    for (const button of this.root.querySelectorAll<HTMLButtonElement>('.lp-demo')) {
      button.addEventListener('click', () => this.runDemo(button.dataset.demo ?? ''));
    }
    q<HTMLButtonElement>(this.root, '.lp-demo-stop').addEventListener('click', () => this.cancelDemo(true));

    // computer keyboard
    this.listen(window, 'keydown', (e) => this.onKeyDown(e));
    this.listen(window, 'keyup', (e) => this.onKeyUp(e));

    // lifecycle
    this.listen(window, 'blur', () => this.silence());
    this.listen(window, 'pagehide', () => this.dispose());
    this.listen(document, 'visibilitychange', () => {
      if (document.hidden) {
        this.silence();
      } else if (this.client) {
        void resumeContext().catch(() => undefined);
      }
      this.kick();
    });
    const onMotion = (): void => {
      this.stage.reducedMotion = this.reduced.matches;
      this.stage.draw(0);
    };
    this.reduced.addEventListener('change', onMotion);
    this.cleanups.push(() => this.reduced.removeEventListener('change', onMotion));

    const resize = new ResizeObserver(() => {
      this.stage.resize();
      this.syncView();
      this.layoutKeyboard(false);
      this.stage.draw(0);
    });
    resize.observe(this.stageEl);
    resize.observe(this.scroll);
    this.cleanups.push(() => resize.disconnect());
    const io = new IntersectionObserver((entries) => {
      this.visible = entries.some((entry) => entry.isIntersecting);
      this.kick();
    }, { rootMargin: '80px 0px' });
    io.observe(this.root);
    this.cleanups.push(() => io.disconnect());
  }

  private wireScrub(): void {
    const scrub = q<HTMLElement>(this.root, '.lp-scrub');
    let active = -1;
    const seek = (e: PointerEvent): void => {
      const rect = scrub.getBoundingClientRect();
      const t = Math.max(0, Math.min(1, (e.clientX - rect.left) / rect.width));
      this.scroll.scrollLeft = t * this.kbd.offsetWidth - this.scroll.clientWidth / 2;
    };
    scrub.addEventListener('pointerdown', (e) => {
      active = e.pointerId;
      scrub.setPointerCapture(e.pointerId);
      seek(e);
    });
    scrub.addEventListener('pointermove', (e) => { if (e.pointerId === active) seek(e); });
    const stop = (e: PointerEvent): void => { if (e.pointerId === active) active = -1; };
    scrub.addEventListener('pointerup', stop);
    scrub.addEventListener('pointercancel', stop);
  }

  private syncView(): void {
    const split = this.stage.split;
    q<HTMLElement>(this.root, '.lp-view').hidden = split;
    for (const b of this.root.querySelectorAll<HTMLButtonElement>('.lp-view button')) {
      b.setAttribute('aria-pressed', String(b.dataset.view === this.stage.view));
    }
  }

  // ------------------------------------------------------------------ start
  private async start(): Promise<void> {
    if (this.starting || this.client || this.disposed) return;
    this.starting = true;
    const button = q<HTMLButtonElement>(this.root, '.lp-start');
    const status = q<HTMLElement>(this.root, '.lp-note-status');
    button.disabled = true;
    status.textContent = 'Starting the audio clock…';
    try {
      const ctx = await resumeContext();
      status.textContent = 'Loading the engine and assembling 88 keys…';
      const client = await LivePianoClient.create(ctx, {});
      if (this.disposed) {
        client.dispose();
        return;
      }
      this.client = client;
      this.stage.setBoard(client.board);
      this.stage.setGeometry(client.keyGeometry);
      client.onSnapshot((s) => this.onSnapshot(s));
      client.onStats((s) => this.onStats(s));
      client.onError((e) => this.onRefusal(e));
      q<HTMLElement>(this.root, '.lp-note').hidden = true;
      for (const b of this.root.querySelectorAll<HTMLButtonElement>('.lp-demo')) b.disabled = false;
      this.buildHonesty(client.info);
      this.root.classList.add('lp-running');
      this.setPedal(this.pedal, 'ui');
      this.announce(`Live piano ready in ${Math.round(client.initMs)} milliseconds. Play the keys, or use the computer keyboard.`);
      this.kick();
    } catch (error) {
      const message = error instanceof Error ? error.message : String(error);
      status.textContent = `The engine did not start: ${message}`;
      button.disabled = false;
      button.textContent = 'Try again';
    } finally {
      this.starting = false;
    }
  }

  private nudgeStart(): void {
    const status = q<HTMLElement>(this.root, '.lp-note-status');
    if (!this.client && !this.starting) status.textContent = 'Press “Start the live piano” first.';
  }

  // ------------------------------------------------------------------ notes
  /** Strike (counted per key: a key already down is not struck again). */
  press(midi: number, speed: number, source: 'pointer' | 'kbd' | 'midi' | 'demo'): void {
    if (!this.client || midi < FIRST || midi > LAST) return;
    const count = this.down.get(midi) ?? 0;
    this.down.set(midi, count + 1);
    if (source === 'demo') this.demoKeys.add(midi);
    if (count > 0) return;
    const v = Math.max(0.05, Math.min(8, speed));
    this.client.noteOn(midi, v);
    this.keyEls.get(midi)?.classList.add('on');
    setText(q(this.root, '.lp-speed-out'), `${noteName(midi)} at ${v.toFixed(2)} m/s`);
  }

  release(midi: number, source: 'pointer' | 'kbd' | 'midi' | 'demo'): void {
    const count = this.down.get(midi) ?? 0;
    if (count <= 0) return;
    if (source === 'demo') this.demoKeys.delete(midi);
    if (count > 1) {
      this.down.set(midi, count - 1);
      return;
    }
    this.down.delete(midi);
    this.client?.noteOff(midi);
    this.keyEls.get(midi)?.classList.remove('on');
  }

  setPedal(value: number, _source: 'ui' | 'space' | 'midi' | 'demo'): void {
    this.pedal = Math.max(0, Math.min(1, value));
    this.client?.sustain(this.pedal);
    this.paintPedal();
  }

  private paintPedal(): void {
    const range = q<HTMLInputElement>(this.root, '.lp-pedal-range');
    range.value = String(this.pedal);
    range.style.setProperty('--fill', `${(this.pedal * 100).toFixed(0)}%`);
    q<HTMLButtonElement>(this.root, '.lp-pedal-btn').setAttribute('aria-pressed', String(this.pedal > 0));
    setText(q(this.root, '.lp-pedal-out'), this.pedal === 0 ? 'up' : this.pedal === 1 ? 'down' : `${Math.round(this.pedal * 100)} % travel`);
  }

  /** Every key up, the pedal up, the engine told so (blur, hidden tab). */
  private silence(): void {
    this.cancelDemo(false);
    for (const [pointer] of this.pointerKeys) {
      try { this.kbd.releasePointerCapture(pointer); } catch { /* not captured */ }
    }
    this.pointerKeys.clear();
    this.codeKeys.clear();
    this.down.clear();
    for (const key of this.keyEls.values()) key.classList.remove('on');
    this.pedal = 0;
    this.paintPedal();
    this.client?.allNotesOff();
  }

  private keyFromPoint(x: number, y: number): number | null {
    const hit = document.elementFromPoint(x, y)?.closest<HTMLElement>('[data-midi]');
    if (!hit || !this.kbd.contains(hit)) return null;
    return Number(hit.dataset.midi);
  }

  private speedAt(midi: number, clientY: number): number {
    if (this.touchMode === 'fixed') return this.fixedSpeed;
    const key = this.keyEls.get(midi);
    if (!key) return this.fixedSpeed;
    const rect = key.getBoundingClientRect();
    const t = Math.max(0, Math.min(1, (clientY - rect.top) / Math.max(1, rect.height)));
    return SPEED_MIN + (SPEED_MAX - SPEED_MIN) * t;
  }

  private onPointerDown(e: PointerEvent): void {
    if (e.button !== 0 && e.pointerType === 'mouse') return;
    if (!this.client) {
      this.nudgeStart();
      return;
    }
    e.preventDefault();
    try { this.kbd.setPointerCapture(e.pointerId); } catch { /* synthetic pointer */ }
    const midi = this.keyFromPoint(e.clientX, e.clientY);
    if (midi === null) return;
    this.pointerKeys.set(e.pointerId, midi);
    this.press(midi, this.speedAt(midi, e.clientY), 'pointer');
  }

  private onPointerMove(e: PointerEvent): void {
    const current = this.pointerKeys.get(e.pointerId);
    if (current === undefined) return;
    const midi = this.keyFromPoint(e.clientX, e.clientY);
    if (midi === null || midi === current) return;
    this.release(current, 'pointer');
    this.pointerKeys.set(e.pointerId, midi);
    this.press(midi, this.speedAt(midi, e.clientY), 'pointer');
  }

  private onPointerEnd(e: PointerEvent): void {
    const current = this.pointerKeys.get(e.pointerId);
    if (current === undefined) return;
    this.pointerKeys.delete(e.pointerId);
    this.release(current, 'pointer');
  }

  private onKeyDown(e: KeyboardEvent): void {
    if (!this.client || !this.visible || e.ctrlKey || e.metaKey || e.altKey || isTyping(e.target)) return;
    if (e.code === 'Space') {
      e.preventDefault();
      if (!e.repeat) {
        this.pedalBeforeSpace = this.pedal;
        this.setPedal(1, 'space');
      }
      return;
    }
    if (e.code === 'KeyZ' || e.code === 'KeyX') {
      e.preventDefault();
      if (!e.repeat) {
        this.octave = Math.max(1, Math.min(7, this.octave + (e.code === 'KeyZ' ? -1 : 1)));
        const low = 12 * (this.octave + 1);
        setText(q(this.root, '.lp-range'), `${noteName(low)}–${noteName(Math.min(LAST, low + 16))}`);
      }
      return;
    }
    const offset = KEY_OFFSETS[e.code];
    if (offset === undefined) return;
    e.preventDefault();
    if (e.repeat || this.codeKeys.has(e.code)) return;
    const midi = 12 * (this.octave + 1) + offset;
    if (midi < FIRST || midi > LAST) return;
    this.codeKeys.set(e.code, midi);
    this.press(midi, this.fixedSpeed, 'kbd');
  }

  private onKeyUp(e: KeyboardEvent): void {
    if (e.code === 'Space') {
      if (this.client && this.visible && !isTyping(e.target)) {
        e.preventDefault();
        this.setPedal(this.pedalBeforeSpace, 'space');
      }
      return;
    }
    const midi = this.codeKeys.get(e.code);
    if (midi === undefined) return;
    this.codeKeys.delete(e.code);
    this.release(midi, 'kbd');
  }

  private async connectMidi(): Promise<void> {
    const out = q<HTMLElement>(this.root, '.lp-midi-out');
    if (!this.client) {
      out.textContent = 'start the piano first';
      return;
    }
    try {
      const access = await navigator.requestMIDIAccess();
      this.midiAccess = access;
      const attach = (): void => {
        const names: string[] = [];
        access.inputs.forEach((input) => {
          input.onmidimessage = (event) => this.onMidi(event);
          names.push(input.name ?? 'input');
        });
        out.textContent = names.length ? names.join(', ') : 'no input found';
      };
      access.onstatechange = attach;
      attach();
    } catch (error) {
      out.textContent = error instanceof Error ? error.message : 'access refused';
    }
  }

  private onMidi(event: MIDIMessageEvent): void {
    const data = event.data;
    if (!data || data.length < 2) return;
    const kind = (data[0] ?? 0) & 0xf0;
    const a = data[1] ?? 0;
    const b = data[2] ?? 0;
    if (kind === 0x90 && b > 0) this.press(a, midiVelocityToMS(b), 'midi');
    else if (kind === 0x80 || (kind === 0x90 && b === 0)) this.release(a, 'midi');
    else if (kind === 0xb0 && a === 64) this.setPedal(b / 127, 'midi');
  }

  // ------------------------------------------------------------------ demos
  private runDemo(name: string): void {
    const demo = DEMOS[name];
    if (!demo || !this.client) return;
    this.cancelDemo(true);
    this.demoName = name;
    setText(q(this.root, '.lp-demo-text'), demo.text);
    q<HTMLElement>(this.root, '.lp-demo-stop').hidden = false;
    for (const b of this.root.querySelectorAll<HTMLButtonElement>('.lp-demo')) {
      b.classList.toggle('is-running', b.dataset.demo === name);
    }
    const steps = demo.build(this);
    let last = 0;
    for (const step of steps) {
      last = Math.max(last, step.at);
      this.demoTimers.push(window.setTimeout(step.run, step.at));
    }
    this.demoTimers.push(window.setTimeout(() => this.finishDemo(), last + 400));
  }

  private finishDemo(): void {
    this.demoTimers = [];
    this.demoName = null;
    q<HTMLElement>(this.root, '.lp-demo-stop').hidden = true;
    for (const b of this.root.querySelectorAll<HTMLButtonElement>('.lp-demo')) b.classList.remove('is-running');
  }

  private cancelDemo(releasePedal: boolean): void {
    for (const timer of this.demoTimers) window.clearTimeout(timer);
    const wasRunning = this.demoName !== null;
    for (const midi of [...this.demoKeys]) this.release(midi, 'demo');
    this.demoKeys.clear();
    if (wasRunning && releasePedal) this.setPedal(0, 'demo');
    this.finishDemo();
  }

  // ------------------------------------------------------------------ engine reports
  private onSnapshot(snapshot: PianoSnapshot): void {
    const latency = this.client?.latencyEstimateMs();
    // Show each state when it is heard: after the queued audio, the
    // receiver's travel and filter delay, and the output device.
    const delay = snapshot.aheadMs + (latency ? latency.receiverMs + latency.outputMs : 0);
    this.queue.push({ due: performance.now() + delay, snapshot });
    if (this.queue.length > 24) this.queue.splice(0, this.queue.length - 24);
    this.kick();
  }

  private onStats(stats: PianoStats): void {
    const previous = this.lastStats;
    this.lastStats = stats;
    const set = (key: string, html: string): void => {
      const node = this.root.querySelector(`[data-s="${key}"]`);
      if (node && node.innerHTML !== html) node.innerHTML = html;
    };
    const latency = this.client?.latencyEstimateMs();
    set('voices', String(stats.voices));
    set('modes', `${stats.dynamicStrings}<small>·</small>${stats.dynamicModes.toLocaleString('en-US')}`);
    set('rtf', stats.realtimeFactor ? `${stats.realtimeFactor.toFixed(2)}<small>×</small>` : '–');
    set('latency', latency ? `${latency.totalMs.toFixed(0)}<small>ms</small>` : '–');
    set('underruns', String(stats.underruns));
    set('peak', `${stats.peakPa.toFixed(2)}<small>Pa</small>`);
    set('clipped', String(stats.clippedSamples));
    set('forced', String(stats.forcedRetirements));
    set('quality', `${stats.quality.maxDynamicModes.toLocaleString('en-US')}<small>modes</small> ${stats.quality.courseModeBudget}<small>per key</small>`);
    if (previous) {
      const now = performance.now();
      if (stats.governorSteps > previous.governorSteps) {
        this.announce(`Rendering ran close to real time on this device, so newly struck keys now get at most ${stats.quality.courseModeBudget} string modes (${stats.quality.maxDynamicModes} in total).`, true);
      } else if (stats.forcedRetirements > previous.forcedRetirements && now - this.lastAnnounce > 4000) {
        this.announce(`The mode budget was full, so the quietest sounding key was returned to rest (${stats.forcedRetirements} so far). Its remaining energy is booked as retired.`, true);
      } else if (stats.contactFailures > previous.contactFailures) {
        this.announce('A strike could not be resolved at this clock and was withdrawn; its energy is booked in the ledger.', true);
      }
    }
  }

  private onRefusal(error: LiveRefusal): void {
    if (error.code === 'refused') {
      this.announce(error.message, true);
      return;
    }
    const box = q<HTMLElement>(this.root, '.lp-refusal');
    box.hidden = false;
    setText(q(box, 'p'), error.code === 'engine' ? 'The engine declined a step and said why:' : 'The engine reported a problem:');
    setText(q(box, 'pre'), error.message);
  }

  private announce(text: string, visible = false): void {
    this.lastAnnounce = performance.now();
    const region = q<HTMLElement>(this.root, '[data-announce]');
    region.textContent = text;
    if (visible) setText(q(this.root, '.lp-demo-text'), text);
  }

  // ------------------------------------------------------------------ drawing
  private running(): boolean {
    return this.visible && !document.hidden && !this.disposed;
  }

  private kick(): void {
    if (!this.raf && this.running()) this.raf = requestAnimationFrame((t) => this.frame(t));
  }

  private frame(now: number): void {
    this.raf = 0;
    if (!this.running()) return;
    const dt = this.lastFrame ? Math.min(100, now - this.lastFrame) : 16;
    this.lastFrame = now;
    let shown: PianoSnapshot | null = null;
    while (this.queue.length && (this.queue[0]?.due ?? Infinity) <= now) shown = this.queue.shift()?.snapshot ?? shown;
    if (shown) {
      this.lastShown = shown;
      this.stage.show(shown);
      this.paintLedger(shown);
      this.paintSounding(shown);
    }
    this.stage.draw(dt);
    if (this.client) this.raf = requestAnimationFrame((t) => this.frame(t));
  }

  private paintScales(scales: StageScales): void {
    const strings = scales.laneM > 0
      ? `strings: lane <b>±${scales.laneM >= 1e-3 ? `${(scales.laneM * 1e3).toFixed(2)} mm` : `${(scales.laneM * 1e6).toFixed(0)} µm`}</b> · height <b>×${Math.round(scales.exaggeration).toLocaleString('en-US')}</b>`
      : 'strings: at rest';
    const board = scales.boardM > 0 && this.lastShown && this.lastShown.keys.length
      ? `board: colour <b>±${scales.boardM >= 1e-6 ? `${(scales.boardM * 1e6).toFixed(2)} µm` : `${(scales.boardM * 1e9).toFixed(0)} nm`}</b> <span class="lp-up">up</span> / <span class="lp-down">down</span>`
      : 'board: at rest';
    const hudS = q<HTMLElement>(this.root, '.lp-hud-strings');
    const hudB = q<HTMLElement>(this.root, '.lp-hud-board');
    const showS = scales.split || this.stage.view === 'strings';
    const showB = scales.split || this.stage.view === 'board';
    if (hudS.innerHTML !== strings) hudS.innerHTML = strings;
    if (hudB.innerHTML !== board) hudB.innerHTML = board;
    hudS.hidden = !showS;
    hudB.hidden = !showB;
  }

  private paintLedger(s: PianoSnapshot): void {
    const l = s.ledger;
    const basis = Math.max(l.inputJ, l.storedJ + l.dissipatedJ + l.retiredJ, 1e-12);
    const bar = this.root.querySelectorAll<HTMLElement>('.lp-ledger .ledger-bar i');
    const parts = [l.storedJ, l.dissipatedJ, l.retiredJ];
    bar.forEach((node, i) => { node.style.flexBasis = `${((Math.max(0, parts[i] ?? 0) / basis) * 100).toFixed(2)}%`; });
    const set = (k: string, text: string): void => setText(q(this.root, `[data-k="${k}"]`), text);
    set('input', millijoules(l.inputJ));
    set('stored', millijoules(l.storedJ));
    set('dissipated', millijoules(l.dissipatedJ));
    set('retired', millijoules(l.retiredJ));
    set('residual', scientific(l.residualJ));
    set('breakdown', `of which strings ${millijoules(l.stringLossJ)} · board ${millijoules(l.boardLossJ)} · dampers ${millijoules(l.damperLossJ)} · felt and shank ${millijoules(l.hammerLossJ)} · hammer caught ${millijoules(l.catchLossJ)}`);
  }

  /** Keys that sound with no finger on them (pedal or still ringing) get the held outline. */
  private paintSounding(s: PianoSnapshot): void {
    const now = new Set<number>();
    for (const k of s.keys) if (k.relativeEnergyJ > 1e-7) now.add(k.midi);
    for (const midi of this.sounding) if (!now.has(midi)) this.keyEls.get(midi)?.classList.remove('held');
    for (const midi of now) this.keyEls.get(midi)?.classList.toggle('held', !this.down.has(midi));
    this.sounding = now;
  }

  // ------------------------------------------------------------------ teardown
  dispose(): void {
    if (this.disposed) return;
    this.silence();
    this.disposed = true;
    if (this.raf) cancelAnimationFrame(this.raf);
    this.raf = 0;
    for (const cleanup of this.cleanups.splice(0)) cleanup();
    if (this.midiAccess) {
      this.midiAccess.inputs.forEach((input) => { input.onmidimessage = null; });
      this.midiAccess.onstatechange = null;
    }
    this.client?.dispose();
    this.client = null;
  }
}

/** Build the live piano inside `root` (called once, when the section scrolls into view). */
export function mountLivePiano(root: HTMLElement): void {
  new LivePianoSection(root);
}
