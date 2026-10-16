// The site's LIVE DRUM section: a three.js stage over FrankenSim's own drum
// mechanics running in WebAssembly, with rebuild-from-geometry controls, the
// live energy ledger and an honest account of what the reduced model omits.
//
// Nothing is fetched and no audio starts until the visitor presses
// "Start the live drum".

import './live-drum.css';
import { getContext, resumeContext } from '../audio';
import {
  DRUM_PRESETS,
  LiveDrumClient,
  type DrumCalibration,
  type DrumInfo,
  type DrumMonitor,
  type DrumParamPatch,
  type DrumParams,
  type DrumReady,
  type DrumSnapshot,
  type DrumStats,
} from './drumClient';
import { LiveRefusal } from './liveAudio';
import { DrumStage, type StageMode, type StageView } from './drumStage';

type PresetName = keyof typeof DRUM_PRESETS;
type ParamKey = 'batter_tension_n_m' | 'resonant_tension_n_m' | 'radius_m' | 'depth_m' | 'batter_thickness_m' | 'damping_ratio';

interface SliderSpec {
  key: ParamKey | 'speed';
  label: string;
  min: number;
  max: number;
  step: number;
  initial: number;
  /** Slider runs on log10 of the value. */
  log?: boolean;
  format(value: number): string;
}

const SPEED: SliderSpec = {
  key: 'speed', label: 'Stick speed', min: 0.2, max: 8, step: 0.05, initial: 2,
  format: (v) => `${v.toFixed(2)} m/s`,
};

const DRUM_SLIDERS: SliderSpec[] = [
  { key: 'batter_tension_n_m', label: 'Batter tension', min: 800, max: 6000, step: 25, initial: 3000, format: (v) => `${v.toFixed(0)} N/m` },
  { key: 'resonant_tension_n_m', label: 'Resonant tension', min: 400, max: 4000, step: 25, initial: 1500, format: (v) => `${v.toFixed(0)} N/m` },
  {
    key: 'radius_m', label: 'Head radius', min: 0.1, max: 0.25, step: 0.0005, initial: 0.1703,
    format: (v) => `${(v * 1000).toFixed(1)} mm · ${((2 * v) / 0.0254).toFixed(1)} in dia.`,
  },
  {
    key: 'depth_m', label: 'Shell depth', min: 0.06, max: 0.4, step: 0.001, initial: 0.1651,
    format: (v) => `${(v * 1000).toFixed(0)} mm · ${(v / 0.0254).toFixed(1)} in`,
  },
  {
    key: 'batter_thickness_m', label: 'Batter film', min: 0.0001, max: 0.0004, step: 0.000002, initial: 0.000254,
    format: (v) => `${(v * 1000).toFixed(3)} mm · ${(v / 0.0000254).toFixed(1)} mil`,
  },
  { key: 'damping_ratio', label: 'Damping ratio', min: 0.0003, max: 0.01, step: 0.01, initial: 0.001, log: true, format: (v) => `ζ ${v.toPrecision(2)}` },
];

const PRESET_TEXT: Record<PresetName, { title: string; band: string; mesh: string }> = {
  reference: { title: 'Reference', band: '80–500 Hz', mesh: '5 × 32 mesh' },
  wide: { title: 'Wide', band: '60–800 Hz', mesh: '10 × 48 mesh' },
  wider: { title: 'Wider', band: '60–1000 Hz', mesh: '12 × 64 mesh' },
};

const RATES = [16, 8, 4, 2, 1] as const;

/** Measured natively: relative RMS deviation of batter-microphone pressure from the 768 kHz image. */
const DEVIATION: Array<[string, string]> = [
  ['384 kHz', '0.27 %'],
  ['192 kHz', '0.85 %'],
  ['96 kHz', '2.1 %'],
  ['48 kHz', '5.4 %'],
];

const ENERGY_REF_J = 0.1;
const ENERGY_RANGE_DB = 80;
const AXIS_LO_HZ = 50;
const AXIS_HI_HZ = 1200;

let serial = 0;

function el<K extends keyof HTMLElementTagNameMap>(tag: K, className?: string, text?: string): HTMLElementTagNameMap[K] {
  const node = document.createElement(tag);
  if (className) node.className = className;
  if (text !== undefined) node.textContent = text;
  return node;
}

function q<T extends Element>(root: ParentNode, selector: string): T {
  const found = root.querySelector<T>(selector);
  if (!found) throw new Error(`live drum: missing ${selector}`);
  return found;
}

function joules(value: number): string {
  if (value === 0) return '0 J';
  const v = Math.abs(value);
  if (v >= 1) return `${value.toFixed(3)} J`;
  if (v >= 1e-3) return `${(value * 1e3).toFixed(2)} mJ`;
  if (v >= 1e-6) return `${(value * 1e6).toFixed(1)} µJ`;
  return `${value.toExponential(1)} J`;
}

function setText(node: Element, text: string): void {
  if (node.textContent !== text) node.textContent = text;
}

function presetOf(params: DrumParams): PresetName | null {
  for (const name of Object.keys(DRUM_PRESETS) as PresetName[]) {
    const p = DRUM_PRESETS[name];
    if (p.radial_intervals === params.radial_intervals && p.azimuths === params.azimuths
      && p.band_low_hz === params.band_low_hz && p.band_high_hz === params.band_high_hz) return name;
  }
  return null;
}

function paramValue(params: DrumParams, key: SliderSpec['key']): number {
  switch (key) {
    case 'speed':
      return Number.NaN;
    case 'damping_ratio':
      return params.batter_damping_ratio;
    default:
      return params[key];
  }
}

function axisPercent(frequencyHz: number): number {
  const t = Math.log(frequencyHz / AXIS_LO_HZ) / Math.log(AXIS_HI_HZ / AXIS_LO_HZ);
  return Math.max(0, Math.min(1, t)) * 100;
}

class LiveDrumSection {
  private readonly root: HTMLElement;
  private readonly id = `ld${++serial}`;
  private readonly stage: DrumStage;
  private readonly stageEl: HTMLElement;
  private client: LiveDrumClient | null = null;
  private starting = false;
  private info: DrumInfo | null = null;
  private calibration: DrumCalibration | null = null;
  private speed = SPEED.initial;
  private aim: [number, number] = [0.06, 0.01];
  private aimVisible = false;
  private monitor: DrumMonitor = 'batter';
  private pending: DrumParamPatch = {};
  private debounce = 0;
  private queue: Promise<unknown> = Promise.resolve();
  private strokePeak = 0;
  private lastStrikes = -1;
  private lastDomUpdate = 0;
  private lastAnnounce = 0;
  private announceTimer = 0;
  private modeBars: HTMLElement[] = [];
  private railTicks: HTMLElement[] = [];
  private modeFrequencies: number[] = [];
  private readonly sliders = new Map<string, { spec: SliderSpec; input: HTMLInputElement; output: HTMLOutputElement }>();

  constructor(root: HTMLElement) {
    this.root = root;
    root.classList.add('ld', 'live-scope');
    root.innerHTML = this.skeleton();
    this.stageEl = q<HTMLElement>(root, '.ld-stage');
    this.stage = new DrumStage(this.stageEl, {
      strike: (x, y, source) => this.strike(x, y, source),
      scale: (factor) => {
        setText(q(this.root, '.ld-scale b'), factor > 0 ? `×${factor.toLocaleString('en-US')}` : 'at rest');
        setText(q(this.root, '.ld-scale .muted'), factor > 0 ? 'auto' : '');
      },
      hint: (text) => this.flashHint(text),
    });
    this.buildSliders();
    this.buildModels();
    this.buildHonesty(null);
    this.wire();
    this.updateHint();
  }

  // ------------------------------------------------------------ markup
  private skeleton(): string {
    const id = this.id;
    return `
<div class="grid wide-left ld-grid">
  <div class="stack ld-left">
    <div class="stage ld-stage pan-y" tabindex="0" role="group" aria-roledescription="drum"
         aria-label="Live drum" aria-describedby="${id}-help">
      <div class="stage-top">
        <div class="seg live ld-seg-sm" role="group" aria-label="View">
          <button type="button" data-view="batter" aria-pressed="true" aria-label="Batter head, seen from above">Top</button>
          <button type="button" data-view="resonant" aria-pressed="false" aria-label="Resonant head, seen from below">Below</button>
        </div>
        <div class="seg live ld-seg-sm" role="group" aria-label="Pointer action">
          <button type="button" data-mode="strike" aria-pressed="true">Strike</button>
          <button type="button" data-mode="orbit" aria-pressed="false">Orbit</button>
        </div>
      </div>
      <div class="stage-hud">
        <span class="ld-scale">displacement <b>at rest</b> <span class="muted"></span></span>
        <span class="ld-hint"></span>
      </div>
      <div class="stage-note ld-note">
        <div class="ld-note-body">
          <span class="tag live">Live engine</span>
          <p class="ld-note-title">A drum solved in your browser, while you play it.</p>
          <p class="small ld-note-text">Starting downloads FrankenSim's drum engine as WebAssembly (about 180&nbsp;kB),
            meshes and eigen-solves both heads, and measures this device to choose the fastest mechanical rate it can
            sustain. The drum starts silent and at rest.</p>
          <button type="button" class="btn live ld-start">Start the live drum</button>
          <p class="small ld-note-status" role="status"></p>
        </div>
      </div>
    </div>
    <p id="${id}-help" class="sr-only">When this area has focus, Space or Enter strikes the batter head at the marked point and
      the arrow keys move that point. With a pointer, press on the head to strike exactly there; drag beside the drum to
      turn it.</p>
    <div class="ld-strokebar">
      <div class="field ld-speed"></div>
      <button type="button" class="btn sm ld-strike-btn" disabled>Strike</button>
    </div>
    <p class="small muted ld-caption">The head is the engine's actual finite-element mesh, coloured by signed
      displacement (<span class="ld-key-up">ivory = up</span>, <span class="ld-key-down">ember = down</span>) and
      drawn from the accepted modal state about 60 times a second, in step with what you hear. The real head moves
      200 to 500 times a second, so the picture is a sampled view.</p>

    <section class="panel ld-ledger-panel" aria-labelledby="${id}-ledger">
      <div class="panel-head">
        <h3 id="${id}-ledger">Energy ledger</h3>
        <span class="tag live">live</span>
      </div>
      <div class="ledger ld-ledger">
        <div class="ledger-bar" role="img" aria-label="Energy ledger">
          <i class="s"></i><i class="d"></i><i class="r"></i>
        </div>
        <div class="ledger-keys">
          <span><i style="background:var(--supplied)"></i>supplied <b class="num" data-k="supplied">0 J</b></span>
          <span><i style="background:var(--stored)"></i>stored <b class="num" data-k="stored">0 J</b></span>
          <span><i style="background:var(--lost)"></i>dissipated <b class="num" data-k="dissipated">0 J</b></span>
          <span><i style="background:var(--supplied);opacity:.55"></i>recovered by the hand <b class="num" data-k="recovered">0 J</b></span>
        </div>
        <p class="ld-residual"><span class="muted">balance residual</span>
          <b class="num" data-k="residual">0 J</b>
          <span class="muted small" data-k="relative"></span></p>
      </div>
      <p class="small muted ld-ledger-note">Every joule the stick brings in is accounted for: still stored in the heads,
        the air spring and the contact; dissipated by damping and contact loss; or carried off by the rebounding stick
        when the hand catches it. The residual is what the arithmetic fails to close. Nothing is corrected.</p>
    </section>

    <dl class="stats ld-stats">
      <div class="stat"><dt>Mechanical rate</dt><dd data-s="rate">–</dd></div>
      <div class="stat"><dt>Real time</dt><dd data-s="rtf">–</dd></div>
      <div class="stat"><dt>Latency</dt><dd data-s="latency">–</dd></div>
      <div class="stat"><dt>Underruns</dt><dd data-s="underruns">–</dd></div>
      <div class="stat"><dt>FMA</dt><dd data-s="fma">–</dd></div>
      <div class="stat"><dt>Modes</dt><dd data-s="modes">–</dd></div>
      <div class="stat"><dt>Last rebuild</dt><dd data-s="rebuild">–</dd></div>
      <div class="stat"><dt>Stroke peak</dt><dd data-s="peak">–</dd></div>
    </dl>
    <p class="small ld-warn" hidden></p>
  </div>

  <div class="stack ld-right">
    <section class="panel" aria-labelledby="${id}-drum">
      <div class="panel-head"><h3 id="${id}-drum">The drum</h3></div>
      <div class="ld-fields"></div>
      <p class="small ld-status" role="status">Changes rebuild the drum from geometry: new mesh, new eigen-solve.</p>
      <div class="refusal ld-refusal" hidden><p></p><pre></pre>
        <button type="button" class="btn sm ld-reset" hidden>Reset the voice</button></div>
    </section>

    <section class="panel" aria-labelledby="${id}-model">
      <div class="panel-head"><h3 id="${id}-model">Model and clock</h3></div>
      <div class="stack ld-model">
        <div class="field"><span class="label"><span>Head model</span><output class="ld-preset-now"></output></span>
          <div class="seg live ld-presets" role="group" aria-label="Head model"></div></div>
        <div class="field"><span class="label"><span>Mechanical rate</span><output class="ld-rate-now"></output></span>
          <div class="seg live ld-rates" role="group" aria-label="Mechanical rate"></div>
          <div class="ld-rate-row"><button type="button" class="btn sm ghost ld-calibrate" disabled>Auto-calibrate</button>
            <span class="small muted ld-rate-note">Auto picks the highest rate that runs at least 1.6× faster than
              real time on this device during a 20 Hz roll.</span></div></div>
        <div class="field"><span class="label"><span>Microphone</span></span>
          <div class="seg live ld-monitor" role="group" aria-label="Microphone">
            <button type="button" data-monitor="batter" aria-pressed="true">Batter side<small>0.35 m above</small></button>
            <button type="button" data-monitor="resonant" aria-pressed="false">Resonant side<small>0.35 m below</small></button>
            <button type="button" data-monitor="split" aria-pressed="false">Both<small>L top · R bottom</small></button>
          </div></div>
      </div>
    </section>

    <section class="panel" aria-labelledby="${id}-modes">
      <div class="panel-head"><h3 id="${id}-modes">Modes</h3><span class="small muted ld-modes-count"></span></div>
      <div class="ld-rail" aria-hidden="true"><div class="ld-window"></div><div class="ld-rail-ticks"></div></div>
      <div class="ld-axis small num" aria-hidden="true"><span>50</span><span>100</span><span>200</span><span>500</span><span>1 kHz</span></div>
      <div class="ld-mode-groups"></div>
      <p class="small muted ld-modes-note">Frequencies come only from the eigen-solve of the heads you set; nothing is
        tuned. Bars show each mode's energy, (v² + ω²q²)/2 of the mass-normalised state, on an 80 dB scale below
        100 mJ. The shaded band is the retained eigen-window.</p>
    </section>
  </div>
</div>
<section class="panel ld-honesty"></section>
<div class="sr-only" aria-live="polite" data-announce></div>`;
  }

  private buildSliders(): void {
    const make = (spec: SliderSpec, host: HTMLElement): void => {
      const inputId = `${this.id}-${spec.key}`;
      const field = host.classList.contains('field') ? host : host.appendChild(el('div', 'field'));
      const label = el('label');
      label.htmlFor = inputId;
      label.append(el('span', undefined, spec.label));
      const output = el('output');
      output.htmlFor.add(inputId);
      label.append(output);
      const input = el('input');
      input.type = 'range';
      input.id = inputId;
      input.min = String(spec.log ? Math.log10(spec.min) : spec.min);
      input.max = String(spec.log ? Math.log10(spec.max) : spec.max);
      input.step = String(spec.step);
      field.append(label, input);
      this.sliders.set(spec.key, { spec, input, output });
      this.setSlider(spec.key, spec.initial);
      input.addEventListener('input', () => {
        const value = this.sliderValue(spec.key);
        this.paintSlider(spec.key);
        if (spec.key === 'speed') {
          this.speed = value;
          return;
        }
        if (spec.key === 'damping_ratio') this.pending.damping_ratio = value;
        else this.pending[spec.key] = value;
        this.scheduleRebuild();
      });
    };
    make(SPEED, q<HTMLElement>(this.root, '.ld-speed'));
    const fields = q<HTMLElement>(this.root, '.ld-fields');
    for (const spec of DRUM_SLIDERS) make(spec, fields);
  }

  private sliderValue(key: string): number {
    const slider = this.sliders.get(key);
    if (!slider) return Number.NaN;
    const raw = Number(slider.input.value);
    return slider.spec.log ? 10 ** raw : raw;
  }

  private setSlider(key: string, value: number): void {
    const slider = this.sliders.get(key);
    if (!slider || !Number.isFinite(value)) return;
    slider.input.value = String(slider.spec.log ? Math.log10(value) : value);
    this.paintSlider(key);
  }

  private paintSlider(key: string): void {
    const slider = this.sliders.get(key);
    if (!slider) return;
    const { spec, input, output } = slider;
    const min = Number(input.min);
    const max = Number(input.max);
    const raw = Number(input.value);
    input.style.setProperty('--fill', `${(((raw - min) / (max - min)) * 100).toFixed(1)}%`);
    output.value = spec.format(spec.log ? 10 ** raw : raw);
  }

  private buildModels(): void {
    const presets = q<HTMLElement>(this.root, '.ld-presets');
    for (const name of Object.keys(DRUM_PRESETS) as PresetName[]) {
      const button = el('button');
      button.type = 'button';
      button.dataset.preset = name;
      button.setAttribute('aria-pressed', String(name === 'reference'));
      button.append(PRESET_TEXT[name].title, el('small', undefined, PRESET_TEXT[name].band), el('small', undefined, PRESET_TEXT[name].mesh));
      presets.append(button);
    }
    const rates = q<HTMLElement>(this.root, '.ld-rates');
    for (const substeps of RATES) {
      const button = el('button');
      button.type = 'button';
      button.dataset.substeps = String(substeps);
      button.setAttribute('aria-pressed', 'false');
      button.disabled = true;
      button.append(`${(48 * substeps).toLocaleString('en-US')} kHz`, el('small', undefined, '–'));
      rates.append(button);
    }
  }

  private buildHonesty(info: DrumInfo | null): void {
    const host = q<HTMLElement>(this.root, '.ld-honesty');
    host.replaceChildren();
    const details = el('details', 'evidence');
    details.append(el('summary', undefined, 'What this live model is — and is not'));
    const body = el('div', 'evidence-body ld-evidence');
    body.append(
      el('p', undefined, info?.provenance ?? 'The reduced real-time model: FrankenSim\'s own linear drum mechanics compiled to WebAssembly. Press Start to load it; the engine then states its own provenance here.'),
    );
    if (info) body.append(el('p', 'muted', info.no_claim));

    const lists = el('div', 'ld-incl');
    const column = (title: string, items: string[], cls: string): HTMLElement => {
      const box = el('div');
      box.append(el('h4', undefined, title));
      const list = el('ul', `ld-list ${cls}`);
      for (const item of items) list.append(el('li', undefined, item));
      box.append(list);
      return box;
    };
    lists.append(
      column('In the live model', [
        'Both heads as finite-element films: bending and membrane prestress on the actual mesh.',
        'Their eigenmodes inside the declared window, solved in the browser on every rebuild.',
        'The sealed air between the heads as one uniform-pressure spring coupling them.',
        'A rigid stick with the tip mass of a 5A profile, launched at the speed you set, with Hertz contact.',
        'Exact zero-order-hold modal stepping with the contact force solved every mechanical step.',
        'Sound pressure from the time-domain Rayleigh integral over each moving head, at a fixed microphone.',
      ], 'in'),
      column('Left out (the offline reference has most of these)', [
        'Boundary-element radiation of the whole closed drum. Here each head radiates into its own infinite baffle, which overstates low-frequency level.',
        'Radiation loading on the heads, and air modes inside the shell.',
        'Head stretching at large amplitude (pitch glide), snare wires and vents.',
        'Shell, hoop and lug flexibility, stick flexure and the hand.',
        'Every mode above the declared window: absent, not approximated.',
      ], 'out'),
    );
    body.append(lists);

    const deviation = el('div', 'table-scroll');
    const table = el('table');
    const head = el('tr');
    head.append(el('th', undefined, 'Mechanical rate'), el('th', undefined, 'Deviation from the 768 kHz image'));
    table.append(head);
    for (const [rate, value] of DEVIATION) {
      const row = el('tr');
      row.append(el('td', 'num', rate), el('td', 'num', value));
      table.append(row);
    }
    deviation.append(table);
    body.append(
      el('h4', undefined, 'How much the mechanical rate matters'),
      deviation,
      el('p', 'small muted', 'Relative RMS difference of the batter-side microphone pressure over 0.5 s after a 2 m/s stroke (reference model), measured against the same engine at 768 kHz, the offline reference rate. The browser build is bit-identical to the native build at the same rate, and its mechanics reproduce the offline drum-modal binary byte for byte.'),
    );

    const chips = (title: string, cls: string, items: string[]): HTMLElement => {
      const box = el('div', 'ld-chiprow');
      box.append(el('h4', undefined, title));
      const row = el('div', 'chips');
      for (const item of items) row.append(el('span', `chip ${cls}`, item));
      box.append(row);
      return box;
    };
    body.append(
      chips('Published dimensions', 'pub', [
        '14 × 6.5 in nominal shell', 'films 0.254 / 0.0762 mm (10 / 3 mil)', '5A stick 406.4 mm × Ø 14.2 mm',
      ]),
      chips('Estimated inputs', 'est', [
        'clear radius 170.3 mm', 'tensions 3000 / 1500 N/m', 'film E 4 GPa, ν 0.38, 1390 kg/m³', 'damping ratio 0.001',
        'stick tip mass 7.24 g', 'Hertz tip 0.8 GPa, R 3 mm', 'air 1.2 kg/m³, 343 m/s', 'microphones 0.35 m',
      ]),
      chips('Not modelled', 'not', [
        'BEM radiation', 'radiation loading', 'head stretching', 'snare wires', 'shell flex', 'air modes', 'vents', 'stick flex',
      ]),
    );
    details.append(body);
    host.append(details);
  }

  // ------------------------------------------------------------ wiring
  private wire(): void {
    q<HTMLButtonElement>(this.root, '.ld-start').addEventListener('click', () => void this.start());
    q<HTMLButtonElement>(this.root, '.ld-strike-btn').addEventListener('click', () => this.strike(this.aim[0], this.aim[1], 'keyboard'));
    q<HTMLButtonElement>(this.root, '.ld-calibrate').addEventListener('click', () => {
      this.enqueue('Measuring this device…', async (client) => {
        const result = await client.calibrate();
        this.calibration = result;
        this.renderRates();
        this.status(result.sustained
          ? `Calibrated: ${result.chosen ? 48 * result.chosen : '?'} kHz runs ${this.factorOf(result.chosen)} real time here.`
          : `No rate reached ${result.minFactor}× real time here; running the cheapest, ${result.chosen ? 48 * result.chosen : '?'} kHz. Expect dropouts.`);
      });
    });
    q<HTMLButtonElement>(this.root, '.ld-reset').addEventListener('click', () => {
      this.enqueue('Resetting…', async (client) => {
        await client.reset();
        this.hideRefusal();
        this.stage.rest();
        this.status('The voice was reset: same drum, at rest.');
      });
    });
    for (const button of this.root.querySelectorAll<HTMLButtonElement>('[data-view]')) {
      button.addEventListener('click', () => {
        const view = (button.dataset.view as StageView | undefined) ?? 'batter';
        this.pressOnly('[data-view]', button);
        this.stage.setView(view);
        this.updateHint();
      });
    }
    for (const button of this.root.querySelectorAll<HTMLButtonElement>('[data-mode]')) {
      button.addEventListener('click', () => {
        const mode = (button.dataset.mode as StageMode | undefined) ?? 'strike';
        this.pressOnly('[data-mode]', button);
        this.stage.setMode(mode);
        this.updateHint();
      });
    }
    for (const button of this.root.querySelectorAll<HTMLButtonElement>('[data-monitor]')) {
      button.addEventListener('click', () => {
        const monitor = (button.dataset.monitor as DrumMonitor | undefined) ?? 'batter';
        this.monitor = monitor;
        this.pressOnly('[data-monitor]', button);
        if (this.client) void this.client.setMonitor(monitor).catch((error: unknown) => this.refuse(error));
      });
    }
    for (const button of this.root.querySelectorAll<HTMLButtonElement>('[data-preset]')) {
      button.addEventListener('click', () => {
        const name = (button.dataset.preset as PresetName | undefined) ?? 'reference';
        this.pressOnly('[data-preset]', button);
        if (!this.client) {
          Object.assign(this.pending, DRUM_PRESETS[name]);
          return;
        }
        this.enqueue(`Rebuilding with the ${PRESET_TEXT[name].title.toLowerCase()} model, then measuring this device…`, async (client) => {
          const ready = await client.configure({ ...DRUM_PRESETS[name] });
          const result = await client.calibrate();
          this.calibration = result;
          this.renderRates();
          const modes = ready.info.heads[0].modes + ready.info.heads[1].modes;
          this.status(`Rebuilt in ${ready.buildMs.toFixed(0)} ms (eigen-solve re-ran): ${modes} modes, calibrated to ${result.chosen ? 48 * result.chosen : '?'} kHz.`);
        });
      });
    }
    for (const button of this.root.querySelectorAll<HTMLButtonElement>('[data-substeps]')) {
      button.addEventListener('click', () => {
        const substeps = Number(button.dataset.substeps) as 1 | 2 | 4 | 8 | 16;
        this.enqueue('Changing the mechanical rate…', async (client) => {
          const ready = await client.setSubsteps(substeps);
          const factor = this.measuredFactor(substeps);
          const slow = factor !== null && factor < 1;
          this.status(`Now ${this.rateText()} (rebuilt in ${ready.buildMs.toFixed(0)} ms, eigen-solve reused).${slow ? ' This rate measured slower than real time here: expect dropouts.' : ''}`);
        });
      });
    }

    // Keyboard on the stage.
    this.stageEl.addEventListener('keydown', (event) => {
      if (!this.client) return;
      const step = event.shiftKey ? 0.0025 : 0.01;
      let moved = false;
      if (event.key === ' ' || event.key === 'Enter') {
        event.preventDefault();
        this.strike(this.aim[0], this.aim[1], 'keyboard');
        return;
      }
      if (event.key === 'ArrowLeft') { this.aim[0] -= step; moved = true; }
      if (event.key === 'ArrowRight') { this.aim[0] += step; moved = true; }
      if (event.key === 'ArrowUp') { this.aim[1] += step; moved = true; }
      if (event.key === 'ArrowDown') { this.aim[1] -= step; moved = true; }
      if (!moved) return;
      event.preventDefault();
      this.clampAim();
      this.stage.setAim(this.aim[0], this.aim[1], true);
      this.announce(`Strike point ${this.describe(this.aim)}`);
    });
    this.stageEl.addEventListener('focus', () => {
      this.aimVisible = true;
      this.stage.setAim(this.aim[0], this.aim[1], !!this.client);
    });
    this.stageEl.addEventListener('blur', () => {
      this.aimVisible = false;
      this.stage.setAim(this.aim[0], this.aim[1], false);
    });

    // Audio interruptions and page lifecycle.
    const context = (): AudioContext | null => (this.client ? getContext() : null);
    document.addEventListener('visibilitychange', () => {
      const ctx = context();
      if (!document.hidden && ctx && ctx.state !== 'running') void resumeContext().catch(() => this.showResume());
    });
    window.addEventListener('pagehide', () => this.stop('The page was hidden, so the engine was stopped. Start it again to play.'));
  }

  private pressOnly(selector: string, pressed: HTMLElement): void {
    for (const button of this.root.querySelectorAll<HTMLElement>(selector)) button.setAttribute('aria-pressed', String(button === pressed));
  }

  // ------------------------------------------------------------ engine
  private async start(): Promise<void> {
    if (this.client || this.starting) return;
    this.starting = true;
    const button = q<HTMLButtonElement>(this.root, '.ld-start');
    const note = q<HTMLElement>(this.root, '.ld-note-status');
    button.disabled = true;
    button.textContent = 'Starting…';
    note.textContent = 'Loading the engine, meshing and eigen-solving both heads, measuring this device…';
    // Resume inside the gesture, before anything awaits.
    const resumed = resumeContext();
    try {
      const context = await resumed;
      const client = await LiveDrumClient.create(context, {
        params: { ...this.pending },
        autoCalibrate: true,
        snapshotHz: 60,
        monitor: this.monitor,
        syncSnapshotsToAudio: true,
      });
      this.pending = {};
      this.client = client;
      this.calibration = client.ready.calibration ?? null;
      client.onReady((ready) => this.applyReady(ready));
      client.onSnapshot((snapshot) => this.onSnapshot(snapshot));
      client.onStats((stats) => this.onStats(stats));
      client.onCalibration((result) => {
        this.calibration = result;
        this.renderRates();
      });
      client.onError((error) => this.refuse(error, true));
      context.addEventListener('statechange', this.onContextState);
      this.applyReady(client.ready);
      this.buildHonesty(client.ready.info);
      this.stageEl.classList.remove('pan-y');
      q<HTMLElement>(this.root, '.ld-note').hidden = true;
      for (const b of this.root.querySelectorAll<HTMLButtonElement>('.ld-strike-btn, .ld-calibrate, [data-substeps]')) b.disabled = false;
      this.stage.setStarted(true);
      this.stage.setAim(this.aim[0], this.aim[1], false);
      this.stageEl.focus({ preventScroll: true });
      this.status(`Built in ${client.ready.buildMs.toFixed(0)} ms: both heads meshed and eigen-solved; running at ${this.rateText()}.`);
      this.announce(`Live drum ready at ${this.rateText()}. Press the head, or Space, to strike.`, true);
    } catch (error) {
      const refusal = error instanceof LiveRefusal ? error : null;
      note.textContent = refusal ? `Refused (${refusal.code}): ${refusal.message}` : `Could not start: ${String((error as Error)?.message ?? error)}`;
      button.disabled = false;
      button.textContent = 'Try again';
    } finally {
      this.starting = false;
    }
  }

  private stop(message: string): void {
    if (!this.client) return;
    getContext().removeEventListener('statechange', this.onContextState);
    this.client.dispose();
    this.client = null;
    this.stage.setStarted(false);
    this.stage.rest();
    this.stageEl.classList.add('pan-y');
    const note = q<HTMLElement>(this.root, '.ld-note');
    note.hidden = false;
    const button = q<HTMLButtonElement>(this.root, '.ld-start');
    button.disabled = false;
    button.textContent = 'Start the live drum';
    q<HTMLElement>(this.root, '.ld-note-status').textContent = message;
    for (const b of this.root.querySelectorAll<HTMLButtonElement>('.ld-strike-btn, .ld-calibrate, [data-substeps]')) b.disabled = true;
  }

  private readonly onContextState = (): void => {
    if (!this.client) return;
    if (getContext().state === 'running') this.hideResume();
    else this.showResume();
  };

  private showResume(): void {
    if (!this.client || this.root.querySelector('.ld-resume')) return;
    const overlay = el('div', 'stage-note ld-resume');
    const body = el('div', 'ld-note-body');
    const button = el('button', 'btn live', 'Resume audio');
    button.type = 'button';
    button.addEventListener('click', () => void resumeContext().then(() => this.hideResume()).catch(() => undefined));
    body.append(el('p', undefined, 'The browser paused audio. The drum is waiting where it was.'), button);
    overlay.append(body);
    this.stageEl.append(overlay);
  }

  private hideResume(): void {
    this.root.querySelector('.ld-resume')?.remove();
  }

  private enqueue(message: string, task: (client: LiveDrumClient) => Promise<void>): void {
    this.queue = this.queue.then(async () => {
      const client = this.client;
      if (!client) return;
      this.status(message);
      this.hideRefusal();
      try {
        await task(client);
      } catch (error) {
        this.refuse(error);
      }
    });
  }

  private scheduleRebuild(): void {
    if (!this.client) {
      this.status('These values will be used when the drum starts.');
      return;
    }
    window.clearTimeout(this.debounce);
    this.debounce = window.setTimeout(() => {
      const patch = this.pending;
      this.pending = {};
      if (Object.keys(patch).length === 0) return;
      this.enqueue('Rebuilding the drum from geometry…', async (client) => {
        const ready = await client.configure(patch);
        this.status(`Rebuilt in ${ready.buildMs.toFixed(0)} ms — eigen-solve ${ready.eigenSolved ? 're-ran' : 'reused'}. The new drum starts at rest.`);
        this.announce(`Rebuilt. Lowest batter mode ${ready.info.heads[0].frequencies_hz[0]?.toFixed(0) ?? '?'} hertz.`);
      });
    }, 250);
  }

  private applyReady(ready: DrumReady): void {
    const info = ready.info;
    this.info = info;
    const params = info.params;
    this.stage.setGeometry(params.radius_m, params.depth_m);
    this.stage.setMesh(ready.mesh.nodes, ready.mesh.triangles);
    this.stage.rest();
    this.clampAim();
    this.stage.setAim(this.aim[0], this.aim[1], this.aimVisible);
    // Controls show the instrument that actually exists.
    for (const spec of DRUM_SLIDERS) {
      if (!(spec.key in this.pending)) this.setSlider(spec.key, paramValue(params, spec.key));
    }
    const preset = presetOf(params);
    for (const button of this.root.querySelectorAll<HTMLElement>('[data-preset]')) {
      button.setAttribute('aria-pressed', String(button.dataset.preset === preset));
    }
    this.buildModes(info);
    this.renderRates();
    const modes = `${info.heads[0].modes} + ${info.heads[1].modes}`;
    this.stat('modes', modes);
    this.stat('rate', `${info.mechanical_rate_hz / 1000}`, `${info.substeps} mechanical steps per 48 kHz audio frame`, 'kHz');
    this.stat('fma', info.fma.startsWith('hardware') ? 'hardware' : info.fma.startsWith('software') ? 'software' : info.fma,
      `Fused multiply-add path: ${info.fma}. Results are bit-identical either way; only speed differs.`);
    if (ready.buildMs > 0 || ready.eigenSolved) {
      this.stat('rebuild', ready.buildMs.toFixed(0), ready.eigenSolved ? 'meshing and eigen-solve of both heads' : 'eigen-solve reused', 'ms');
    }
    setText(q(this.root, '.ld-preset-now'), `${modes} modes`);
    this.strokePeak = 0;
    this.lastStrikes = -1;
    this.stat('peak', '–');
    q<HTMLElement>(this.root, '.ld-warn').hidden = true;
  }

  private buildModes(info: DrumInfo): void {
    const groups = q<HTMLElement>(this.root, '.ld-mode-groups');
    const ticks = q<HTMLElement>(this.root, '.ld-rail-ticks');
    const previous = new Map<string, number>(this.railTicks.map((tick) => [tick.dataset.key ?? '', Number(tick.dataset.f)]));
    groups.replaceChildren();
    this.modeBars = [];
    this.modeFrequencies = [];
    const fresh: HTMLElement[] = [];
    for (const [h, head] of info.heads.entries()) {
      const group = el('div', 'ld-mode-group');
      const title = el('h4', undefined, h === 0 ? 'Batter head' : 'Resonant head');
      title.append(el('span', 'muted', ` · ${head.modes} modes`));
      const list = el('ol', 'ld-mode-list');
      for (const [i, frequency] of head.frequencies_hz.entries()) {
        const item = el('li');
        item.append(el('span', 'num', frequency.toFixed(1)));
        const bar = el('span', 'ld-bar');
        const fill = el('i', h === 0 ? 'b' : 'r');
        bar.append(fill);
        item.append(bar);
        list.append(item);
        this.modeBars.push(fill);
        this.modeFrequencies.push(frequency);
        const key = `${h}-${i}`;
        const tick = this.railTicks.find((t) => t.dataset.key === key) ?? el('i', h === 0 ? 'b' : 'r');
        tick.dataset.key = key;
        tick.dataset.f = String(frequency);
        if (!previous.has(key)) tick.style.left = `${axisPercent(frequency)}%`;
        fresh.push(tick);
      }
      group.append(title, list);
      groups.append(group);
    }
    ticks.replaceChildren(...fresh);
    // Commit the old positions, then let the CSS transition glide to the new ones.
    void ticks.offsetWidth;
    for (const tick of fresh) tick.style.left = `${axisPercent(Number(tick.dataset.f))}%`;
    this.railTicks = fresh;
    const band = q<HTMLElement>(this.root, '.ld-window');
    const lo = axisPercent(Math.max(AXIS_LO_HZ, info.params.band_low_hz));
    const hi = axisPercent(info.params.band_high_hz);
    band.style.left = `${lo}%`;
    band.style.width = `${Math.max(0, hi - lo)}%`;
    const coupled = info.coupled_frequencies_hz;
    setText(q(this.root, '.ld-modes-count'), `air-coupled lowest: ${coupled[0]?.toFixed(1) ?? '–'} Hz`);
  }

  private renderRates(): void {
    const info = this.info;
    for (const button of this.root.querySelectorAll<HTMLButtonElement>('[data-substeps]')) {
      const substeps = Number(button.dataset.substeps);
      button.setAttribute('aria-pressed', String(info?.substeps === substeps));
      const small = button.querySelector('small');
      if (small) setText(small, this.factorOf(substeps));
    }
    setText(q(this.root, '.ld-rate-now'), info ? this.rateText() : '');
  }

  private measuredFactor(substeps: number | null): number | null {
    const result = this.calibration?.results.find((r) => r.substeps === substeps);
    return result?.measured && result.realtimeFactor !== undefined ? result.realtimeFactor : null;
  }

  private factorOf(substeps: number | null): string {
    const result = this.calibration?.results.find((r) => r.substeps === substeps);
    if (!result) return '–';
    if (result.measured && result.realtimeFactor !== undefined) return `${result.realtimeFactor.toFixed(1)}×`;
    if (result.predictedRealtimeFactor !== undefined) return `<${Math.max(0.1, result.predictedRealtimeFactor * 1.5).toFixed(1)}×`;
    return 'refused';
  }

  private rateText(): string {
    return this.info ? `${this.info.mechanical_rate_hz / 1000} kHz` : '–';
  }

  // ----------------------------------------------------------- playing
  private strike(x: number, y: number, source: 'pointer' | 'keyboard'): void {
    if (!this.client) return;
    this.aim = [x, y];
    this.clampAim();
    this.stage.showMarker(this.aim[0], this.aim[1]);
    this.stage.setAim(this.aim[0], this.aim[1], this.aimVisible);
    const [ax, ay] = this.aim;
    this.client.strike(ax, ay, this.speed).catch((error: unknown) => {
      if (error instanceof LiveRefusal && error.code === 'strike_outside_head') this.flashHint('That point is off the meshed head: strike further in');
      else this.refuse(error);
    });
    if (source === 'keyboard') this.announce(`Struck at ${this.describe(this.aim)}, ${this.speed.toFixed(1)} metres per second`);
  }

  private clampAim(): void {
    const radius = (this.info?.params.radius_m ?? 0.17) * 0.9;
    const r = Math.hypot(this.aim[0], this.aim[1]);
    if (r > radius) this.aim = [(this.aim[0] * radius) / r, (this.aim[1] * radius) / r];
  }

  private describe([x, y]: [number, number]): string {
    const cm = (v: number): string => `${Math.abs(v * 100).toFixed(0)} centimetres`;
    const horizontal = Math.abs(x) < 0.005 ? 'centre' : `${cm(x)} ${x > 0 ? 'right' : 'left'}`;
    const vertical = Math.abs(y) < 0.005 ? 'centre' : `${cm(y)} ${y > 0 ? 'back' : 'front'}`;
    return `${horizontal}, ${vertical}`;
  }

  private onSnapshot(s: DrumSnapshot): void {
    this.stage.setDisplacements(s.batterDisplacementM, s.resonantDisplacementM);
    if (s.strikesLaunched !== this.lastStrikes) {
      if (s.strikesLaunched > this.lastStrikes && this.lastStrikes >= 0) this.strokePeak = 0;
      this.lastStrikes = s.strikesLaunched;
    }
    const mic = this.monitor === 'resonant' ? s.peakPa[1] : this.monitor === 'split' ? Math.max(s.peakPa[0], s.peakPa[1]) : s.peakPa[0];
    if (mic > this.strokePeak) this.strokePeak = mic;
    const now = performance.now();
    if (now - this.lastDomUpdate < 50) return;
    this.lastDomUpdate = now;

    const ledger = s.ledger;
    const supplied = ledger.suppliedJ;
    const part = (value: number): string => (supplied > 0 ? `${Math.max(0, Math.min(100, (value / supplied) * 100)).toFixed(2)}%` : '0%');
    const bar = this.root.querySelectorAll<HTMLElement>('.ld-ledger .ledger-bar i');
    bar[0]?.style.setProperty('flex-basis', part(ledger.storedJ));
    bar[1]?.style.setProperty('flex-basis', part(ledger.dissipatedJ));
    bar[2]?.style.setProperty('flex-basis', part(ledger.recoveredJ));
    setText(q(this.root, '[data-k="supplied"]'), joules(supplied));
    setText(q(this.root, '[data-k="stored"]'), joules(ledger.storedJ));
    setText(q(this.root, '[data-k="dissipated"]'), joules(ledger.dissipatedJ));
    setText(q(this.root, '[data-k="recovered"]'), joules(ledger.recoveredJ));
    setText(q(this.root, '[data-k="residual"]'), `${ledger.balanceResidualJ.toExponential(2)} J`);
    setText(q(this.root, '[data-k="relative"]'), supplied > 0 ? `${(Math.abs(ledger.balanceResidualJ) / supplied).toExponential(1)} of supplied` : '');

    const modes = this.modeFrequencies.length;
    for (let i = 0; i < modes; i++) {
      const omega = 2 * Math.PI * (this.modeFrequencies[i] ?? 0);
      const qv = s.q[i] ?? 0;
      const vv = s.v[i] ?? 0;
      const energy = 0.5 * (vv * vv + omega * omega * qv * qv);
      const level = energy > 0 ? (10 * Math.log10(energy / ENERGY_REF_J) + ENERGY_RANGE_DB) / ENERGY_RANGE_DB : 0;
      const fraction = Math.max(0, Math.min(1, level));
      this.modeBars[i]?.style.setProperty('transform', `scaleX(${fraction.toFixed(3)})`);
      this.railTicks[i]?.style.setProperty('--e', fraction.toFixed(3));
    }

    if (this.strokePeak > 0) this.stat('peak', `${this.strokePeak.toFixed(this.strokePeak < 10 ? 1 : 0)}`, undefined, 'Pa');
    const warn = q<HTMLElement>(this.root, '.ld-warn');
    const fullScale = this.info?.params.full_scale_pa ?? 400;
    if (s.clippedFrames > 0) {
      warn.hidden = false;
      setText(warn, `${s.clippedFrames} audio frames went past the fixed full scale of ${fullScale} Pa and were clipped by your output. A softer stroke stays clean; nothing here limits or rescales the sound.`);
    }
  }

  private onStats(stats: DrumStats): void {
    this.stat('rtf', stats.realtimeFactor ? `${stats.realtimeFactor.toFixed(1)}×` : '–');
    this.stat('latency', `${stats.latencyEstimateMs.toFixed(0)}`, 'queue + device output + anti-alias filter; acoustic travel to the microphone (about 1 ms) is physics, not latency', 'ms');
    this.stat('underruns', String(stats.underruns), stats.targetRaises > 0 ? `queue raised to ${stats.targetMs.toFixed(0)} ms after dropouts` : undefined);
    if (stats.mechanicalRateHz) this.stat('rate', `${stats.mechanicalRateHz / 1000}`, undefined, 'kHz');
  }

  private stat(key: string, value: string, title?: string, unit?: string): void {
    const node = this.root.querySelector<HTMLElement>(`[data-s="${key}"]`);
    if (!node) return;
    const text = unit ? `${value} ${unit}` : value;
    if (node.dataset.text !== text) {
      node.dataset.text = text;
      node.textContent = value;
      if (unit) node.append(el('small', undefined, unit));
    }
    if (title !== undefined) node.title = title;
  }

  // ------------------------------------------------------------ voice
  private status(text: string): void {
    setText(q(this.root, '.ld-status'), text);
  }

  private refuse(error: unknown, halted = false): void {
    const box = q<HTMLElement>(this.root, '.ld-refusal');
    const refusal = error instanceof LiveRefusal ? error : null;
    box.hidden = false;
    setText(q(box, 'p'), halted
      ? 'The engine stopped on a refusal and is silent. Reset the voice to keep playing this drum.'
      : 'Refused. The previous drum is still the one playing.');
    setText(q(box, 'pre'), refusal ? `${refusal.code}: ${refusal.message}` : String((error as Error)?.message ?? error));
    q<HTMLElement>(box, '.ld-reset').hidden = !halted;
    this.status('');
    // Show the instrument that exists, not the request that was refused.
    if (this.client) {
      this.pending = {};
      const params = this.client.info.params;
      for (const spec of DRUM_SLIDERS) this.setSlider(spec.key, paramValue(params, spec.key));
      const preset = presetOf(params);
      for (const button of this.root.querySelectorAll<HTMLElement>('[data-preset]')) button.setAttribute('aria-pressed', String(button.dataset.preset === preset));
    }
    this.announce(refusal ? `Refused: ${refusal.message}` : 'Refused.', true);
  }

  private hideRefusal(): void {
    q<HTMLElement>(this.root, '.ld-refusal').hidden = true;
  }

  private updateHint(): void {
    const touch = window.matchMedia('(pointer: coarse)').matches;
    const mode = this.root.querySelector('[data-mode][aria-pressed="true"]')?.getAttribute('data-mode') ?? 'strike';
    const below = this.root.querySelector('[data-view="resonant"]')?.getAttribute('aria-pressed') === 'true';
    let hint = mode === 'strike'
      ? (touch ? 'Tap the head to strike · drag outside it, or with two fingers, to turn' : 'Click the head to strike · drag beside the drum to turn it')
      : `Drag to turn the drum · ${touch ? 'tap' : 'click'} the head to strike`;
    if (below) hint = 'Resonant head from below · the stick still lands on the batter head above the point you touch';
    this.hint = hint;
    setText(q(this.root, '.ld-hint'), hint);
  }

  private hint = '';
  private hintTimer = 0;
  private flashHint(text: string): void {
    const node = q(this.root, '.ld-hint');
    setText(node, text);
    window.clearTimeout(this.hintTimer);
    this.hintTimer = window.setTimeout(() => setText(node, this.hint), 1800);
  }

  private announce(text: string, force = false): void {
    const node = q<HTMLElement>(this.root, '[data-announce]');
    const now = performance.now();
    window.clearTimeout(this.announceTimer);
    const wait = force ? 0 : Math.max(0, 1200 - (now - this.lastAnnounce));
    this.announceTimer = window.setTimeout(() => {
      this.lastAnnounce = performance.now();
      node.textContent = text;
    }, wait);
  }
}

/**
 * Build the live drum section inside `root`. Call once, when the section
 * scrolls into view. Nothing is fetched and no audio starts until the visitor
 * presses "Start the live drum".
 */
export function mountLiveDrum(root: HTMLElement): void {
  if (root.dataset.liveDrum === 'mounted') return;
  root.dataset.liveDrum = 'mounted';
  new LiveDrumSection(root);
}
