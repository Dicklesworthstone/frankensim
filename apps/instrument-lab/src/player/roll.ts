// The player's face: a perforated paper roll moving down over a brass tracker bar, the 88 keys it
// plays, and the percussion lanes of an orchestrion beside them. One canvas; nothing animates on
// its own: every hole is a note of the file, every key goes down when the engines are told to
// strike it, and the paper's position is the audio clock.

import type { Arrangement, PieceId } from "../midi/arrange";
import { PIECE_NAMES } from "../midi/arrange";

const BLACK = [false, true, false, true, false, false, true, false, true, false, true, false];
const isBlack = (m: number): boolean => BLACK[m % 12]!;
const PIECE_ABBR: Record<PieceId, string> = { kick: "K", snare: "S", "snare-off": "S·", "tom-12": "T", "tom-14": "F", "hihat-closed": "H", "hihat-open": "O", crash: "C", ride: "R" };
const CYMBAL = new Set<PieceId>(["hihat-closed", "hihat-open", "crash", "ride"]);

interface Lane { x: number; w: number; black: boolean }

export class RollView {
  readonly canvas: HTMLCanvasElement;
  private g: CanvasRenderingContext2D;
  private W = 0;
  private H = 0;
  private dpr = 1;
  private arr: Arrangement | null = null;
  private lanes = new Map<number, Lane>();
  private keyLo = 21;
  private keyHi = 108;
  private drumX = new Map<PieceId, number>();
  private drumW = 10;
  private paperL = 0;
  private paperR = 0;
  private keysL = 0;
  private keysR = 0;
  private trackerY = 0;
  private keysTop = 0;
  private keysH = 0;
  private textX = 0;
  private pps = 100;
  private texture: HTMLCanvasElement | null = null;
  private maxDur = 1;
  private hitFlash = new Map<PieceId, number>();
  /** score seconds shown above the tracker bar */
  visibleSeconds = 4;
  compact = false;
  onScrub: ((t: number, done: boolean) => void) | null = null;
  scrubbable = true;

  constructor(host: HTMLElement) {
    this.canvas = document.createElement("canvas");
    this.canvas.className = "pp-roll";
    this.canvas.setAttribute("role", "img");
    this.canvas.setAttribute("aria-label", "The paper roll and the keyboard it plays");
    host.append(this.canvas);
    this.g = this.canvas.getContext("2d", { alpha: false })!;
    new ResizeObserver(() => this.resize()).observe(this.canvas);
    this.bindScrub();
    this.resize();
  }

  setArrangement(a: Arrangement | null): void {
    this.arr = a;
    this.maxDur = a ? a.notes.reduce((m, n) => Math.max(m, n.off - n.t), 0.5) : 1;
    this.layout();
  }

  private resize(): void {
    const r = this.canvas.getBoundingClientRect();
    this.dpr = Math.min(2, window.devicePixelRatio || 1);
    this.W = Math.max(10, Math.round(r.width));
    this.H = Math.max(10, Math.round(r.height));
    this.canvas.width = Math.round(this.W * this.dpr);
    this.canvas.height = Math.round(this.H * this.dpr);
    this.compact = this.W < 640;
    this.layout();
    this.paperTexture();
  }

  /** Keys across the width: all 88 on a wide screen, the file's own range (whole octaves) on a narrow one. */
  private layout(): void {
    const a = this.arr;
    if (this.compact && a && a.notes.length) {
      let lo = 108, hi = 21;
      for (const n of a.notes) { lo = Math.min(lo, n.key); hi = Math.max(hi, n.key); }
      lo = Math.max(21, lo - (lo % 12)); // down to a C
      hi = Math.min(108, hi + (11 - (hi % 12))); // up to a B
      while (hi - lo < 35) { if (lo > 21) lo = Math.max(21, lo - 12); if (hi < 108) hi = Math.min(108, hi + 12); if (lo === 21 && hi === 108) break; }
      this.keyLo = lo; this.keyHi = hi;
    } else { this.keyLo = 21; this.keyHi = 108; }
    const pieces = a?.pieces ?? [];
    const pad = this.compact ? 6 : 14;
    const pedalW = this.compact ? 9 : 14;
    const textW = this.compact ? 30 : 54;
    this.drumW = pieces.length ? (this.compact ? 11 : 15) : 0;
    const drumsW = pieces.length ? pieces.length * this.drumW + (this.compact ? 6 : 12) : 0;
    this.paperL = pad;
    this.paperR = this.W - pad;
    this.textX = this.paperL + pedalW + 4;
    this.keysL = this.paperL + pedalW + textW;
    this.keysR = this.paperR - drumsW - (this.compact ? 3 : 8);
    let whites = 0;
    for (let m = this.keyLo; m <= this.keyHi; m++) if (!isBlack(m)) whites++;
    const ww = (this.keysR - this.keysL) / whites;
    this.lanes.clear();
    let wi = 0;
    for (let m = this.keyLo; m <= this.keyHi; m++) {
      if (isBlack(m)) this.lanes.set(m, { x: this.keysL + wi * ww - ww * 0.32, w: ww * 0.64, black: true });
      else { this.lanes.set(m, { x: this.keysL + wi * ww, w: ww, black: false }); wi++; }
    }
    this.drumX.clear();
    const d0 = this.keysR + (this.compact ? 6 : 12);
    pieces.forEach((p, i) => this.drumX.set(p, d0 + i * this.drumW));
    this.keysH = Math.round(Math.min(132, Math.max(this.compact ? 54 : 70, ww * 5.4)));
    this.keysTop = this.H - this.keysH - 4;
    this.trackerY = this.keysTop - (this.compact ? 16 : 22);
    this.pps = (this.trackerY - 18) / this.visibleSeconds;
  }

  /** Paper fibres, drawn once and moved with the paper. */
  private paperTexture(): void {
    const c = document.createElement("canvas");
    c.width = 256; c.height = 512;
    const g = c.getContext("2d")!;
    g.fillStyle = "#ece2cb";
    g.fillRect(0, 0, 256, 512);
    let s = 12345;
    const rnd = (): number => { s = (s * 16807) % 2147483647; return s / 2147483647; };
    for (let i = 0; i < 2600; i++) {
      const x = rnd() * 256, y = rnd() * 512, l = 2 + rnd() * 9;
      g.strokeStyle = rnd() < 0.5 ? `rgba(120, 96, 60, ${0.03 + rnd() * 0.05})` : `rgba(255, 252, 240, ${0.05 + rnd() * 0.08})`;
      g.lineWidth = 0.6;
      g.beginPath();
      g.moveTo(x, y);
      g.lineTo(x + (rnd() - 0.5) * 2, y + l);
      g.stroke();
    }
    this.texture = c;
  }

  flash(piece: PieceId, at: number): void {
    this.hitFlash.set(piece, at);
  }

  // ------------------------------------------------------------------- draw
  /** `t`: score seconds at the tracker bar. `wall`: performance.now() for flashes. */
  draw(t: number, wall: number, playing: boolean): void {
    const g = this.g;
    const { W, H } = this;
    g.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    // the cabinet: dark walnut
    const wood = g.createLinearGradient(0, 0, W, 0);
    wood.addColorStop(0, "#1b120d"); wood.addColorStop(0.5, "#24180f"); wood.addColorStop(1, "#1a110c");
    g.fillStyle = wood;
    g.fillRect(0, 0, W, H);
    this.drawPaper(t, playing);
    this.drawTracker(t);
    this.drawKeys(t, wall);
  }

  private yOf(sec: number, now: number): number {
    return this.trackerY - (sec - now) * this.pps;
  }

  private drawPaper(now: number, playing: boolean): void {
    const g = this.g;
    const a = this.arr;
    const top = 0, bottom = this.trackerY + 8;
    const L = this.paperL, R = this.paperR;
    // paper with fibres that travel with it
    g.save();
    g.beginPath();
    g.rect(L, top, R - L, bottom - top);
    g.clip();
    if (this.texture) {
      const pat = g.createPattern(this.texture, "repeat")!;
      const off = ((now * this.pps) % 512 + 512) % 512;
      pat.setTransform(new DOMMatrix([1, 0, 0, 1, L, off]));
      g.fillStyle = pat;
    } else g.fillStyle = "#ece2cb";
    g.fillRect(L, top, R - L, bottom - top);
    // curvature: the paper leaves the supply spool at the top and darkens toward the edges
    const shade = g.createLinearGradient(L, 0, R, 0);
    shade.addColorStop(0, "rgba(60, 40, 20, 0.28)"); shade.addColorStop(0.035, "rgba(60, 40, 20, 0)");
    shade.addColorStop(0.965, "rgba(60, 40, 20, 0)"); shade.addColorStop(1, "rgba(60, 40, 20, 0.3)");
    g.fillStyle = shade;
    g.fillRect(L, top, R - L, bottom - top);
    if (a) {
      const t0 = now - (bottom - this.trackerY) / this.pps - 0.05;
      const t1 = now + this.trackerY / this.pps + 0.05;
      this.drawPrint(a, now, t0, t1);
      this.drawPedal(a, now, t0, t1);
      this.drawHoles(a, now, t0, t1);
      this.drawDrumHoles(a, now, t0, t1);
    } else {
      g.fillStyle = "rgba(70, 52, 30, 0.55)";
      g.font = `500 ${this.compact ? 12 : 14}px "Fraunces Variable", Georgia, serif`;
      g.textAlign = "center";
      g.fillText("Choose a roll from the shelf, or drop a MIDI file", (L + R) / 2, this.trackerY * 0.45);
    }
    // the supply spool's shadow at the top
    const spool = g.createLinearGradient(0, 0, 0, 34);
    spool.addColorStop(0, "rgba(10, 6, 3, 0.75)"); spool.addColorStop(1, "rgba(10, 6, 3, 0)");
    g.fillStyle = spool;
    g.fillRect(L, 0, R - L, 34);
    g.restore();
    void playing;
  }

  private drawPrint(a: Arrangement, now: number, t0: number, t1: number): void {
    const g = this.g;
    const L = this.keysL, R = this.keysR;
    // bar lines, printed faintly; every fourth bar numbered in the margin
    g.lineWidth = 1;
    for (const b of a.bars) {
      if (b.t < t0 || b.t > t1) continue;
      const y = Math.round(this.yOf(b.t, now)) + 0.5;
      g.strokeStyle = b.n % 4 === 0 ? "rgba(110, 82, 48, 0.32)" : "rgba(110, 82, 48, 0.14)";
      g.beginPath(); g.moveTo(L, y); g.lineTo(R, y); g.stroke();
      if (b.n % 4 === 0 && !this.compact) {
        g.fillStyle = "rgba(110, 82, 48, 0.55)";
        g.font = '500 9px "JetBrains Mono", monospace';
        g.textAlign = "right";
        g.fillText(String(b.n + 1), L - 4, y - 3);
      }
    }
    // sections: a double rule across and the name in small capitals
    for (const s of a.sections) {
      if (s.t < t0 || s.t > t1) continue;
      const y = Math.round(this.yOf(s.t, now)) + 0.5;
      g.strokeStyle = "rgba(122, 42, 32, 0.55)";
      g.beginPath(); g.moveTo(this.paperL + 4, y); g.lineTo(this.paperR - 4, y); g.moveTo(this.paperL + 4, y + 3); g.lineTo(this.paperR - 4, y + 3); g.stroke();
      g.fillStyle = "rgba(122, 42, 32, 0.85)";
      g.font = `600 ${this.compact ? 9 : 10.5}px "JetBrains Mono", monospace`;
      g.textAlign = "left";
      g.fillText(s.text.toUpperCase(), this.keysL + 6, y - 5);
    }
    // chord symbols printed down the margin, like the words on an old song roll
    g.textAlign = "left";
    for (const c of a.chords) {
      if (c.t < t0 - 0.5 || c.t > t1) continue;
      const y = this.yOf(c.t, now);
      const past = c.t <= now;
      g.fillStyle = past ? "rgba(70, 52, 30, 0.45)" : "rgba(52, 36, 18, 0.92)";
      g.font = `600 italic ${this.compact ? 11 : 14}px "Fraunces Variable", Georgia, serif`;
      g.save();
      g.translate(this.textX, y);
      g.fillText(c.label, 0, -2, this.keysL - this.textX - 4);
      g.restore();
    }
  }

  private drawPedal(a: Arrangement, now: number, t0: number, t1: number): void {
    const g = this.g;
    const x = this.paperL + 3, w = (this.compact ? 9 : 14) - 4;
    let down: number | null = null;
    let level = 0;
    const segs: [number, number, number][] = [];
    for (const p of a.pedal) {
      if (p.value >= 0.5 && down === null) { down = p.t; level = p.value; }
      else if (p.value < 0.5 && down !== null) { segs.push([down, p.t, level]); down = null; }
    }
    if (down !== null) segs.push([down, a.durationS, level]);
    for (const [s, e] of segs) {
      if (e < t0 || s > t1) continue;
      const y0 = this.yOf(e, now), y1 = this.yOf(s, now);
      const lit = s <= now && e > now;
      this.hole(x, y0, w, y1 - y0, lit ? 1 : 0, lit ? "#e0b45c" : null);
    }
  }

  private hole(x: number, y: number, w: number, h: number, glow: number, color: string | null): void {
    const g = this.g;
    const r = Math.min(w / 2, 3);
    g.beginPath();
    g.roundRect(x, y, w, Math.max(w * 0.9, h), r);
    if (glow > 0 && color) {
      g.fillStyle = color;
      g.shadowColor = color;
      g.shadowBlur = 10 * glow;
      g.fill();
      g.shadowBlur = 0;
    } else {
      g.fillStyle = "#1d150d";
      g.fill();
      // the cut edge catches light on the lower side
      g.strokeStyle = "rgba(255, 248, 230, 0.35)";
      g.lineWidth = 0.6;
      g.stroke();
    }
  }

  private drawHoles(a: Arrangement, now: number, t0: number, t1: number): void {
    const notes = a.notes;
    // first note that could still be visible
    let lo = 0, hi = notes.length;
    const from = t0 - this.maxDur;
    while (lo < hi) { const m = (lo + hi) >> 1; if (notes[m]!.t < from) lo = m + 1; else hi = m; }
    for (let i = lo; i < notes.length; i++) {
      const n = notes[i]!;
      if (n.t > t1) break;
      if (n.off < t0) continue;
      const lane = this.lanes.get(n.key);
      if (!lane) continue;
      const w = Math.max(2, lane.w * (lane.black ? 0.62 : 0.5));
      const x = lane.x + (lane.w - w) / 2;
      const y0 = this.yOf(n.off, now), y1 = this.yOf(n.t, now);
      const lit = n.t <= now && n.off > now;
      const v = Math.min(1, n.vel / 110);
      this.hole(x, y0, w, y1 - y0 - 1, lit ? 0.5 + v * 0.5 : 0, lit ? `rgb(${Math.round(70 + 60 * v)}, ${Math.round(190 + 40 * v)}, ${Math.round(172 + 30 * v)})` : null);
    }
  }

  private drawDrumHoles(a: Arrangement, now: number, t0: number, t1: number): void {
    const hits = a.hits;
    let lo = 0, hi = hits.length;
    while (lo < hi) { const m = (lo + hi) >> 1; if (hits[m]!.t < t0) lo = m + 1; else hi = m; }
    const g = this.g;
    // lane guides
    for (const [, x] of this.drumX) {
      g.fillStyle = "rgba(110, 82, 48, 0.08)";
      g.fillRect(x + 1, 0, this.drumW - 2, this.trackerY + 8);
    }
    for (let i = lo; i < hits.length; i++) {
      const h = hits[i]!;
      if (h.t > t1) break;
      const x = this.drumX.get(h.piece);
      if (x === undefined) continue;
      const y = this.yOf(h.t, now);
      const d = Math.max(4, Math.min(this.drumW - 3, 3 + h.g * (this.drumW - 6)));
      const lit = now - h.t >= 0 && now - h.t < 0.09;
      g.beginPath();
      g.arc(x + this.drumW / 2, y - d / 2, d / 2, 0, Math.PI * 2);
      if (lit) {
        const c = CYMBAL.has(h.piece) ? "#f0c674" : "#ef8a5a";
        g.fillStyle = c; g.shadowColor = c; g.shadowBlur = 10; g.fill(); g.shadowBlur = 0;
      } else {
        g.fillStyle = "#1d150d"; g.fill();
        g.strokeStyle = "rgba(255, 248, 230, 0.35)"; g.lineWidth = 0.6; g.stroke();
      }
    }
  }

  private drawTracker(now: number): void {
    const g = this.g;
    const y = this.trackerY;
    const L = this.paperL - 4, R = this.paperR + 4;
    const h = this.compact ? 9 : 12;
    const brass = g.createLinearGradient(0, y - h / 2, 0, y + h / 2);
    brass.addColorStop(0, "#f3d68f"); brass.addColorStop(0.35, "#c79a45"); brass.addColorStop(0.7, "#8a6424"); brass.addColorStop(1, "#5b3f14");
    g.fillStyle = brass;
    g.beginPath();
    g.roundRect(L, y - h / 2, R - L, h, 3);
    g.fill();
    g.fillStyle = "rgba(255, 245, 210, 0.5)";
    g.fillRect(L + 2, y - h / 2 + 1, R - L - 4, 1);
    // the holes of the bar: one per key, lit while the key's hole passes over it
    const a = this.arr;
    const active = new Map<number, number>();
    if (a) {
      const notes = a.notes;
      let lo = 0, hi = notes.length;
      while (lo < hi) { const m = (lo + hi) >> 1; if (notes[m]!.t < now - this.maxDur) lo = m + 1; else hi = m; }
      for (let i = lo; i < notes.length && notes[i]!.t <= now; i++) { const n = notes[i]!; if (n.off > now) active.set(n.key, Math.max(active.get(n.key) ?? 0, n.vel)); }
    }
    for (const [m, lane] of this.lanes) {
      const cx = lane.x + lane.w / 2;
      const v = active.get(m);
      g.beginPath();
      g.arc(cx, y, lane.black ? 1.3 : 1.6, 0, Math.PI * 2);
      if (v !== undefined) { g.fillStyle = "#9ff0e2"; g.shadowColor = "#63d2c0"; g.shadowBlur = 8; g.fill(); g.shadowBlur = 0; }
      else { g.fillStyle = "#2a1c09"; g.fill(); }
    }
    this.active = active;
  }
  private active = new Map<number, number>();

  private drawKeys(now: number, wall: number): void {
    const g = this.g;
    const top = this.keysTop, h = this.keysH;
    const L = this.paperL, R = this.paperR;
    // the key slip and a strip of red felt under the fallboard
    g.fillStyle = "#0b0806";
    g.fillRect(L - 4, top - 6, R - L + 8, h + 10);
    g.fillStyle = "#6e1f1a";
    g.fillRect(this.keysL, top - 3, this.keysR - this.keysL, 3);
    const lit = this.active;
    for (const pass of [false, true]) {
      for (const [m, lane] of this.lanes) {
        if (lane.black !== pass) continue;
        const v = lit.get(m);
        const down = v !== undefined;
        const k = down ? Math.min(1, v / 110) : 0;
        if (!lane.black) {
          const y = top + (down ? 2 : 0);
          const grad = g.createLinearGradient(0, y, 0, y + h);
          if (down) { grad.addColorStop(0, "#cfc5b1"); grad.addColorStop(1, "#b9e9df"); }
          else { grad.addColorStop(0, "#f7f0e2"); grad.addColorStop(0.85, "#e6dccb"); grad.addColorStop(1, "#cfc4b0"); }
          g.fillStyle = grad;
          g.fillRect(lane.x + 0.5, y, lane.w - 1, h - (down ? 2 : 0));
          if (down) { g.fillStyle = `rgba(99, 210, 192, ${0.25 + 0.5 * k})`; g.fillRect(lane.x + 0.5, y + h * 0.55, lane.w - 1, h * 0.45 - 2); }
        } else {
          const bh = h * 0.62;
          const y = top + (down ? 1.5 : 0);
          const grad = g.createLinearGradient(0, y, 0, y + bh);
          if (down) { grad.addColorStop(0, "#1c3a36"); grad.addColorStop(1, "#2f8f82"); }
          else { grad.addColorStop(0, "#2b2d33"); grad.addColorStop(0.8, "#121317"); grad.addColorStop(1, "#08090b"); }
          g.fillStyle = grad;
          g.fillRect(lane.x, y, lane.w, bh);
          g.fillStyle = down ? "rgba(160, 240, 226, 0.35)" : "rgba(255, 255, 255, 0.16)";
          g.fillRect(lane.x + 1, y + bh - (down ? 5 : 7), lane.w - 2, 2);
        }
      }
    }
    // the percussion pads under their lanes
    const a = this.arr;
    if (a) {
      for (const [p, x] of this.drumX) {
        const at = this.hitFlash.get(p) ?? -1e9;
        const f = Math.max(0, 1 - (wall - at) / 180);
        const cx = x + this.drumW / 2, cy = top + h * 0.38, r = this.drumW * 0.44;
        g.beginPath();
        g.arc(cx, cy, r, 0, Math.PI * 2);
        const cym = CYMBAL.has(p);
        g.fillStyle = f > 0 ? (cym ? `rgba(240, 198, 116, ${0.35 + 0.65 * f})` : `rgba(239, 138, 90, ${0.35 + 0.65 * f})`) : cym ? "#5a4320" : "#3b2a22";
        g.fill();
        g.strokeStyle = cym ? "#c79a45" : "#a8a29a";
        g.lineWidth = 1;
        g.stroke();
        g.fillStyle = f > 0.3 ? "#1b1407" : "rgba(241, 234, 219, 0.75)";
        g.font = `600 ${this.compact ? 7.5 : 9}px "JetBrains Mono", monospace`;
        g.textAlign = "center";
        g.textBaseline = "middle";
        g.fillText(PIECE_ABBR[p], cx, cy + 0.5);
        g.textBaseline = "alphabetic";
      }
      if (!this.compact && this.drumX.size) {
        g.save();
        g.fillStyle = "rgba(201, 194, 179, 0.6)";
        g.font = '500 9px "JetBrains Mono", monospace';
        g.textAlign = "center";
        const xs = [...this.drumX.values()];
        g.fillText("DRUMS", (xs[0]! + xs[xs.length - 1]! + this.drumW) / 2, top + h * 0.8);
        g.restore();
      }
    }
    void now;
  }

  /** The piece under a lane at canvas x (for the legend tooltip). */
  pieceAt(x: number): PieceId | null {
    for (const [p, px] of this.drumX) if (x >= px && x < px + this.drumW) return p;
    return null;
  }
  pieceName(p: PieceId): string { return PIECE_NAMES[p]; }

  // ------------------------------------------------------------------ scrub
  private bindScrub(): void {
    let drag: { y: number; t: number } | null = null;
    this.canvas.addEventListener("pointerdown", (e) => {
      // a finger scrolls the page; the paper is dragged with a mouse or a pen (touch uses the slider)
      if (e.pointerType === "touch" || !this.scrubbable || !this.arr || !this.onScrub) return;
      drag = { y: e.clientY, t: this.scrubFrom() };
      this.canvas.setPointerCapture(e.pointerId);
      this.canvas.classList.add("dragging");
    });
    this.canvas.addEventListener("pointermove", (e) => {
      if (!drag || !this.onScrub) return;
      const t = Math.max(0, Math.min(this.arr?.durationS ?? 0, drag.t + (e.clientY - drag.y) / this.pps));
      this.onScrub(t, false);
    });
    const end = (e: PointerEvent): void => {
      if (!drag || !this.onScrub) return;
      const t = Math.max(0, Math.min(this.arr?.durationS ?? 0, drag.t + (e.clientY - drag.y) / this.pps));
      drag = null;
      this.canvas.classList.remove("dragging");
      this.onScrub(t, true);
    };
    this.canvas.addEventListener("pointerup", end);
    this.canvas.addEventListener("pointercancel", end);
  }
  /** Set by the owner: the position a drag starts from. */
  scrubFrom: () => number = () => 0;
}
