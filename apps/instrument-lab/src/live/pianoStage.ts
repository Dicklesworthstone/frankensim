// Canvas stage of the live piano: the sounding strings and the soundboard,
// drawn from the engine's snapshots (no animation of its own: every pixel is
// a value the mechanics produced). One canvas; the section owns the rAF loop.

import type { PianoBoard, PianoKeySnapshot, PianoSnapshot } from './pianoClient';

export type PianoView = 'strings' | 'board';

const NAMES = ['C', 'C♯', 'D', 'D♯', 'E', 'F', 'F♯', 'G', 'G♯', 'A', 'A♯', 'B'];
export function noteName(midi: number): string {
  return `${NAMES[((midi % 12) + 12) % 12] ?? '?'}${Math.floor(midi / 12) - 1}`;
}

const VERDIGRIS: [number, number, number] = [99, 210, 192];
const EMBER: [number, number, number] = [232, 106, 60];
const INSIDE: [number, number, number] = [22, 26, 33];
const EDGE: [number, number, number] = [92, 90, 84];
const REST_M = 0.02; // hammer rest below the string plane (engine CATCH_DISTANCE)

interface Rect { x: number; y: number; w: number; h: number }

interface Raster {
  cols: number;
  rows: number;
  /** Nearest surface point per pixel, -1 outside the board. */
  nearest: Int32Array;
  edge: Uint8Array;
  minX: number;
  minY: number;
  spanX: number;
  spanY: number;
  canvas: HTMLCanvasElement;
  image: ImageData;
}

export interface StageScales {
  /** Displacement at a lane edge [m], 0 when nothing sounds. */
  laneM: number;
  /** Vertical exaggeration relative to the longest sounding string drawn across the width. */
  exaggeration: number;
  /** Board colour full scale [m]. */
  boardM: number;
  split: boolean;
}

export class PianoStage {
  private readonly canvas: HTMLCanvasElement;
  private readonly ctx: CanvasRenderingContext2D;
  private cssW = 0;
  private cssH = 0;
  private dpr = 1;
  private board: PianoBoard | null = null;
  private modes = 0;
  private points = 0;
  private raster: Raster | null = null;
  private field = new Float32Array(0);
  private geometry = new Map<number, { length: number; strike: number; unison: number; f1: number }>();
  private snapshot: PianoSnapshot | null = null;
  private fieldDirty = true;
  private laneScale = 0; // px per metre
  private boardScale = 0; // metres at full colour
  private readonly envelopes = new Map<number, Float32Array>();
  private readonly fade = new Map<number, number>();
  view: PianoView = 'strings';
  reducedMotion = false;
  onScales: (scales: StageScales) => void = () => {};

  constructor(private readonly host: HTMLElement) {
    this.canvas = document.createElement('canvas');
    this.canvas.setAttribute('aria-hidden', 'true');
    host.prepend(this.canvas);
    const ctx = this.canvas.getContext('2d', { alpha: false });
    if (!ctx) throw new Error('2D canvas is unavailable');
    this.ctx = ctx;
    this.resize();
  }

  get split(): boolean {
    return this.cssW >= 720;
  }

  resize(): void {
    const rect = this.host.getBoundingClientRect();
    this.dpr = Math.min(2, window.devicePixelRatio || 1);
    this.cssW = Math.max(1, Math.round(rect.width));
    this.cssH = Math.max(1, Math.round(rect.height));
    this.canvas.width = Math.round(this.cssW * this.dpr);
    this.canvas.height = Math.round(this.cssH * this.dpr);
    this.canvas.style.width = `${this.cssW}px`;
    this.canvas.style.height = `${this.cssH}px`;
  }

  setGeometry(values: Float64Array): void {
    this.geometry.clear();
    for (let i = 0; i + 4 < values.length; i += 5) {
      this.geometry.set(values[i] ?? 0, {
        length: values[i + 1] ?? 1, strike: values[i + 2] ?? 0.12, unison: values[i + 3] ?? 1, f1: values[i + 4] ?? 0,
      });
    }
  }

  setBoard(board: PianoBoard): void {
    this.board = board;
    this.points = Math.floor(board.xy.length / 2);
    this.modes = this.points > 0 ? Math.floor(board.shapes.length / this.points) : 0;
    this.field = new Float32Array(this.points);
    this.raster = this.buildRaster(board);
    this.fieldDirty = true;
  }

  show(snapshot: PianoSnapshot): void {
    this.snapshot = snapshot;
    this.fieldDirty = true;
  }

  clear(): void {
    this.snapshot = null;
    this.envelopes.clear();
    this.fade.clear();
    this.fieldDirty = true;
  }

  /** Nearest-point raster of the FE surface quadrature: the board's own outline. */
  private buildRaster(board: PianoBoard): Raster {
    let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
    for (let i = 0; i < this.points; i++) {
      const x = board.xy[2 * i] ?? 0, y = board.xy[2 * i + 1] ?? 0;
      minX = Math.min(minX, x); maxX = Math.max(maxX, x); minY = Math.min(minY, y); maxY = Math.max(maxY, y);
    }
    const spacing = Math.sqrt(((maxX - minX) * (maxY - minY)) / Math.max(1, this.points));
    const pad = spacing * 0.9;
    minX -= pad; minY -= pad; maxX += pad; maxY += pad;
    const pixel = spacing * 0.4;
    const cols = Math.max(8, Math.ceil((maxX - minX) / pixel));
    const rows = Math.max(8, Math.ceil((maxY - minY) / pixel));
    const cell = spacing * 1.6;
    const gx = Math.ceil((maxX - minX) / cell) + 1;
    const gy = Math.ceil((maxY - minY) / cell) + 1;
    const buckets: number[][] = Array.from({ length: gx * gy }, () => []);
    for (let i = 0; i < this.points; i++) {
      const cx = Math.floor(((board.xy[2 * i] ?? 0) - minX) / cell);
      const cy = Math.floor(((board.xy[2 * i + 1] ?? 0) - minY) / cell);
      buckets[cy * gx + cx]?.push(i);
    }
    const nearest = new Int32Array(cols * rows).fill(-1);
    const near = new Uint8Array(cols * rows);
    const limit2 = (spacing * 1.3) ** 2;
    for (let r = 0; r < rows; r++) {
      for (let c = 0; c < cols; c++) {
        const x = minX + (c + 0.5) * ((maxX - minX) / cols);
        const y = minY + (r + 0.5) * ((maxY - minY) / rows);
        const cx = Math.floor((x - minX) / cell), cy = Math.floor((y - minY) / cell);
        let best = -1, bestD = Infinity;
        for (let dy = -2; dy <= 2; dy++) {
          for (let dx = -2; dx <= 2; dx++) {
            const bucket = buckets[(cy + dy) * gx + (cx + dx)];
            if (!bucket || cx + dx < 0 || cx + dx >= gx) continue;
            for (const i of bucket) {
              const ddx = (board.xy[2 * i] ?? 0) - x, ddy = (board.xy[2 * i + 1] ?? 0) - y;
              const d = ddx * ddx + ddy * ddy;
              if (d < bestD) { bestD = d; best = i; }
            }
          }
        }
        nearest[r * cols + c] = best;
        if (bestD <= limit2) near[r * cols + c] = 1;
      }
    }
    // Outside = pixels far from every point AND connected to the raster's
    // border; enclosed gaps between sparse quadrature points stay inside.
    const outside = new Uint8Array(cols * rows);
    const stack: number[] = [];
    const seed = (i: number): void => { if (!near[i] && !outside[i]) { outside[i] = 1; stack.push(i); } };
    for (let c = 0; c < cols; c++) { seed(c); seed((rows - 1) * cols + c); }
    for (let r = 0; r < rows; r++) { seed(r * cols); seed(r * cols + cols - 1); }
    while (stack.length) {
      const i = stack.pop() ?? 0;
      const r = Math.floor(i / cols), c = i % cols;
      if (r > 0) seed(i - cols);
      if (r < rows - 1) seed(i + cols);
      if (c > 0) seed(i - 1);
      if (c < cols - 1) seed(i + 1);
    }
    const edge = new Uint8Array(cols * rows);
    for (let r = 0; r < rows; r++) {
      for (let c = 0; c < cols; c++) {
        const i = r * cols + c;
        if (outside[i]) { nearest[i] = -1; continue; }
        const out = (rr: number, cc: number): boolean =>
          rr < 0 || cc < 0 || rr >= rows || cc >= cols || outside[rr * cols + cc] === 1;
        if (out(r - 1, c) || out(r + 1, c) || out(r, c - 1) || out(r, c + 1)) edge[i] = 1;
      }
    }
    const canvas = document.createElement('canvas');
    canvas.width = cols;
    canvas.height = rows;
    const image = new ImageData(cols, rows);
    return { cols, rows, nearest, edge, minX, minY, spanX: maxX - minX, spanY: maxY - minY, canvas, image };
  }

  // ------------------------------------------------------------------ draw
  draw(dt: number): void {
    const ctx = this.ctx;
    ctx.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    ctx.fillStyle = '#0a0b0f';
    ctx.fillRect(0, 0, this.cssW, this.cssH);
    const top = 66, bottom = 34;
    const area: Rect = { x: 0, y: top, w: this.cssW, h: Math.max(40, this.cssH - top - bottom) };
    let laneM = 0, exaggeration = 0;
    if (this.split) {
      const sw = Math.round(this.cssW * 0.6);
      const strings: Rect = { x: area.x + 14, y: area.y, w: sw - 24, h: area.h };
      const board: Rect = { x: sw + 10, y: area.y, w: this.cssW - sw - 24, h: area.h };
      ctx.strokeStyle = 'rgba(236,226,205,0.1)';
      ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(sw + 0.5, top - 24); ctx.lineTo(sw + 0.5, this.cssH - 16); ctx.stroke();
      [laneM, exaggeration] = this.drawStrings(strings, dt);
      this.drawBoard(board, dt);
    } else if (this.view === 'strings') {
      [laneM, exaggeration] = this.drawStrings({ x: 12, y: area.y, w: this.cssW - 24, h: area.h }, dt);
    } else {
      this.drawBoard({ x: 12, y: area.y, w: this.cssW - 24, h: area.h }, dt);
    }
    this.onScales({ laneM, exaggeration, boardM: this.boardScale, split: this.split });
  }

  private label(text: string, x: number, y: number, color = '#8f897d', align: CanvasTextAlign = 'left'): void {
    const ctx = this.ctx;
    ctx.font = '500 10.5px "JetBrains Mono", ui-monospace, Menlo, monospace';
    ctx.textAlign = align;
    ctx.textBaseline = 'middle';
    ctx.fillStyle = color;
    ctx.fillText(text, x, y);
  }

  private drawStrings(rect: Rect, dt: number): [number, number] {
    const ctx = this.ctx;
    const snap = this.snapshot;
    const keys: PianoKeySnapshot[] = snap ? [...snap.keys].sort((a, b) => b.midi - a.midi) : [];
    this.label(this.split ? 'STRINGS  ·  shape y(x) along x/L' : 'shape y(x) along x/L', rect.x, rect.y - 12);
    const labelW = 46;
    const x0 = rect.x + labelW, x1 = rect.x + rect.w - 8;
    const w = Math.max(10, x1 - x0);
    if (keys.length === 0) {
      ctx.strokeStyle = 'rgba(236,226,205,0.12)';
      ctx.beginPath(); ctx.moveTo(x0, rect.y + rect.h / 2 + 0.5); ctx.lineTo(x1, rect.y + rect.h / 2 + 0.5); ctx.stroke();
      this.label('no key is sounding', x0 + w / 2, rect.y + rect.h / 2 - 16, '#5f5b53', 'center');
      this.label('agraffe', x0, rect.y + rect.h / 2 + 14, '#5f5b53');
      this.label('bridge', x1, rect.y + rect.h / 2 + 14, '#5f5b53', 'right');
      this.envelopes.clear();
      return [0, 0];
    }
    const lane = Math.min(150, (rect.h - 14) / keys.length);
    const used = lane * keys.length;
    const start = rect.y + 14 + (rect.h - 14 - used) / 2;
    // Common vertical scale: the largest displacement fills 40 % of a lane.
    let peak = 0;
    let longest = 0;
    for (const k of keys) {
      for (const y of k.stringShapeM) peak = Math.max(peak, Math.abs(y));
      peak = Math.max(peak, Math.abs(k.bridgeDisplacementM));
      longest = Math.max(longest, this.geometry.get(k.midi)?.length ?? 0);
    }
    if (peak > 0) {
      const target = (lane * 0.4) / peak;
      if (this.laneScale === 0 || target < this.laneScale || this.reducedMotion) this.laneScale = target;
      else this.laneScale += (target - this.laneScale) * Math.min(1, dt / 900);
    }
    const scale = this.laneScale;
    const alive = new Set<number>();
    keys.forEach((k, index) => {
      alive.add(k.midi);
      const base = Math.round(start + (index + 0.5) * lane) + 0.5;
      const energy = Math.max(1e-14, k.relativeEnergyJ);
      const targetAlpha = Math.max(0.16, Math.min(1, (Math.log10(energy) + 9.5) / 7));
      const previous = this.fade.get(k.midi) ?? targetAlpha;
      const alpha = previous + (targetAlpha - previous) * Math.min(1, dt / 250);
      this.fade.set(k.midi, alpha);
      // lane rule
      ctx.strokeStyle = 'rgba(236,226,205,0.08)';
      ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(x0, base); ctx.lineTo(x1, base); ctx.stroke();
      ctx.fillStyle = 'rgba(236,226,205,0.25)';
      ctx.fillRect(x0 - 1, base - 5, 2, 10); // agraffe
      ctx.fillRect(x1 - 2, base - 6, 3, 12); // bridge pin
      // note label
      this.label(noteName(k.midi), rect.x, base - (lane > 34 ? 6 : 0), k.held ? '#f1eadb' : '#c9c2b3');
      if (lane > 34) this.label(k.damped ? 'damped' : k.held ? 'key down' : 'pedal', rect.x, base + 8, k.damped ? '#b5694f' : '#5f5b53');
      // envelope (peak hold) and instantaneous shape
      const n = k.stringShapeM.length;
      let env = this.envelopes.get(k.midi);
      if (!env || env.length !== n) { env = new Float32Array(n); this.envelopes.set(k.midi, env); }
      const decay = Math.exp(-dt / 140);
      for (let j = 0; j < n; j++) env[j] = Math.max(Math.abs(k.stringShapeM[j] ?? 0), (env[j] ?? 0) * decay);
      const xs = (j: number): number => x0 + ((j + 0.5) / n) * w;
      const clampY = (v: number): number => Math.max(-lane * 0.48, Math.min(lane * 0.48, v * scale));
      ctx.beginPath();
      ctx.moveTo(x0, base);
      for (let j = 0; j < n; j++) ctx.lineTo(xs(j), base - clampY(env[j] ?? 0));
      ctx.lineTo(x1, base);
      for (let j = n - 1; j >= 0; j--) ctx.lineTo(xs(j), base + clampY(env[j] ?? 0));
      ctx.closePath();
      ctx.fillStyle = `rgba(99,210,192,${(0.13 * alpha).toFixed(3)})`;
      ctx.fill();
      if (!this.reducedMotion) {
        ctx.beginPath();
        ctx.moveTo(x0, base);
        for (let j = 0; j < n; j++) ctx.lineTo(xs(j), base - clampY(k.stringShapeM[j] ?? 0));
        ctx.lineTo(x1, base - clampY(k.bridgeDisplacementM));
        ctx.strokeStyle = `rgba(99,210,192,${alpha.toFixed(3)})`;
        ctx.lineWidth = 1.6;
        ctx.lineJoin = 'round';
        ctx.stroke();
      }
      // hammer at its strike station; travel 0..20 mm drawn to the lane, not exaggerated
      const g = this.geometry.get(k.midi);
      if (g) {
        const hx = x0 + g.strike * w;
        const travel = Math.max(0, Math.min(1.1, -k.hammerPositionM / REST_M));
        const room = Math.max(6, lane * 0.42 - 4);
        const hy = base + 3 + travel * room;
        const contact = k.hammerForceN > 0.01;
        ctx.fillStyle = contact ? '#e86a3c' : `rgba(182,73,63,${Math.max(0.35, alpha).toFixed(3)})`;
        ctx.beginPath();
        ctx.roundRect(hx - 6, hy, 12, Math.max(5, Math.min(10, lane * 0.12)), 3);
        ctx.fill();
        if (contact) {
          const bar = Math.min(lane * 0.4, 4 + Math.log10(1 + k.hammerForceN) * 12);
          ctx.fillStyle = '#e86a3c';
          ctx.fillRect(hx + 8, base - bar, 2, bar);
          this.label(`${k.hammerForceN.toFixed(k.hammerForceN < 10 ? 1 : 0)} N`, hx + 13, base - bar + 5, '#f2a184');
        }
      }
    });
    for (const midi of [...this.envelopes.keys()]) if (!alive.has(midi)) { this.envelopes.delete(midi); this.fade.delete(midi); }
    this.label('agraffe', x0, rect.y + 4, '#5f5b53');
    this.label('bridge', x1, rect.y + 4, '#5f5b53', 'right');
    const laneM = scale > 0 ? (lane * 0.5) / scale : 0;
    const exaggeration = longest > 0 && scale > 0 ? scale / (w / longest) : 0;
    return [laneM, exaggeration];
  }

  private drawBoard(rect: Rect, dt: number): void {
    const ctx = this.ctx;
    this.label(this.split ? 'SOUNDBOARD  ·  Model D, signed displacement' : 'Model D, signed displacement', rect.x, rect.y - 12);
    const raster = this.raster;
    const board = this.board;
    if (!raster || !board) {
      ctx.strokeStyle = 'rgba(236,226,205,0.12)';
      ctx.strokeRect(rect.x + 0.5, rect.y + 0.5, rect.w - 1, rect.h - 1);
      return;
    }
    const snap = this.snapshot;
    if (this.fieldDirty) {
      this.fieldDirty = false;
      const coords = snap?.boardCoordinates;
      let peak = 0;
      for (let p = 0; p < this.points; p++) {
        let w = 0;
        if (coords) {
          const row = p * this.modes;
          for (let j = 0; j < this.modes; j++) w += (board.shapes[row + j] ?? 0) * (coords[j] ?? 0);
        }
        // Reduced motion: a slowly released magnitude instead of the
        // instantaneous (aliased, flickering) signed field.
        this.field[p] = this.reducedMotion ? Math.max(Math.abs(w), (this.field[p] ?? 0) * 0.86) : w;
        peak = Math.max(peak, Math.abs(this.field[p] ?? 0));
      }
      if (peak > 0) {
        if (this.boardScale === 0 || peak > this.boardScale || this.reducedMotion) this.boardScale = peak;
        else this.boardScale += (peak - this.boardScale) * Math.min(1, dt / 600);
      }
      this.paintField(raster);
    }
    // Fit, rotating the long axis along the wider side of the rectangle.
    const rotate = (rect.w > rect.h) !== (raster.spanX > raster.spanY);
    const bw = rotate ? raster.spanY : raster.spanX;
    const bh = rotate ? raster.spanX : raster.spanY;
    const s = Math.min(rect.w / bw, rect.h / bh) * 0.96;
    const cx = rect.x + rect.w / 2, cy = rect.y + rect.h / 2;
    ctx.save();
    ctx.translate(cx, cy);
    if (rotate) ctx.rotate(-Math.PI / 2);
    // Board +y points up on screen (the drawing's convention), x to the right.
    ctx.scale(1, -1);
    ctx.imageSmoothingEnabled = true;
    ctx.imageSmoothingQuality = 'high';
    const dw = raster.spanX * s, dh = raster.spanY * s;
    ctx.drawImage(raster.canvas, -dw / 2, -dh / 2, dw, dh);
    ctx.restore();
  }

  private paintField(raster: Raster): void {
    const data = raster.image.data;
    const scale = this.boardScale > 0 ? this.boardScale : 1;
    const { cols, rows, nearest, edge } = raster;
    for (let r = 0; r < rows; r++) {
      // ImageData rows run top-down; the raster's row 0 is the smallest y.
      for (let c = 0; c < cols; c++) {
        const i = r * cols + c;
        const o = i * 4;
        const p = nearest[i] ?? -1;
        if (p < 0) { data[o + 3] = 0; continue; }
        if (edge[i]) {
          data[o] = EDGE[0]; data[o + 1] = EDGE[1]; data[o + 2] = EDGE[2]; data[o + 3] = 255;
          continue;
        }
        const t = Math.max(-1, Math.min(1, (this.field[p] ?? 0) / scale));
        const a = Math.abs(t) ** 0.7;
        const hue = this.reducedMotion ? VERDIGRIS : t >= 0 ? VERDIGRIS : EMBER;
        data[o] = INSIDE[0] + (hue[0] - INSIDE[0]) * a;
        data[o + 1] = INSIDE[1] + (hue[1] - INSIDE[1]) * a;
        data[o + 2] = INSIDE[2] + (hue[2] - INSIDE[2]) * a;
        data[o + 3] = 255;
      }
    }
    raster.canvas.getContext('2d')?.putImageData(raster.image, 0, 0);
  }
}
