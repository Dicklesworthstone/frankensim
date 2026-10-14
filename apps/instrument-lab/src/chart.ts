import { dpr, el, fmt, hidpi } from "./util";

export interface Series {
  label: string;
  color: string;
  x: ArrayLike<number>;
  y: ArrayLike<number>;
  axis?: "left" | "right";
  width?: number;
  dash?: number[];
  scatter?: boolean;
}
export interface ChartOptions {
  xLabel: string;
  yLabel: string;
  y2Label?: string;
  xLog?: boolean;
  y2Log?: boolean;
  aspect?: number;
  xRange?: [number, number];
  yRange?: [number, number];
  y2Range?: [number, number];
  xFormat?: (v: number) => string;
  yFormat?: (v: number) => string;
  y2Format?: (v: number) => string;
  zeroLine?: boolean;
  bands?: { from: number; to: number; label?: string }[];
}

function niceTicks(lo: number, hi: number, count = 5): number[] {
  const span = hi - lo;
  if (!(span > 0)) return [lo];
  const step0 = span / count;
  const mag = 10 ** Math.floor(Math.log10(step0));
  const step = [1, 2, 2.5, 5, 10].map((m) => m * mag).find((s) => span / s <= count + 1) ?? mag * 10;
  const out: number[] = [];
  for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-9; v += step) out.push(Math.abs(v) < step * 1e-9 ? 0 : v);
  return out;
}

function extent(series: Series[], key: "x" | "y", axis?: "left" | "right"): [number, number] {
  let lo = Infinity;
  let hi = -Infinity;
  for (const s of series) {
    if (axis && (s.axis ?? "left") !== axis) continue;
    const a = s[key];
    for (let i = 0; i < a.length; i++) {
      const v = a[i]!;
      if (!Number.isFinite(v)) continue;
      if (v < lo) lo = v;
      if (v > hi) hi = v;
    }
  }
  if (lo === hi) {
    lo -= 1;
    hi += 1;
  }
  return [lo, hi];
}

/** A small, honest line chart: real axes, units, hover readout. */
export class LineChart {
  private series: Series[] = [];
  private hoverX: number | null = null;
  private readonly tip: HTMLElement;

  constructor(private readonly canvas: HTMLCanvasElement, private readonly opts: ChartOptions, legend?: HTMLElement | null) {
    this.tip = el("div", { class: "chart-tip" });
    canvas.parentElement?.append(this.tip);
    new ResizeObserver(() => this.draw()).observe(canvas);
    canvas.addEventListener("pointermove", (e) => {
      const r = canvas.getBoundingClientRect();
      this.hoverX = e.clientX - r.left;
      this.draw();
    });
    canvas.addEventListener("pointerleave", () => {
      this.hoverX = null;
      this.tip.style.opacity = "0";
      this.draw();
    });
    this.legend = legend ?? null;
  }
  private legend: HTMLElement | null;

  set(series: Series[]): void {
    this.series = series;
    if (this.legend) {
      this.legend.replaceChildren(...series.map((s) => {
        const i = el("i");
        i.style.background = s.color;
        if (s.dash) i.style.background = `repeating-linear-gradient(to right, ${s.color} 0 4px, transparent 4px 7px)`;
        return el("span", {}, i, s.label + (s.axis === "right" ? " (right axis)" : ""));
      }));
    }
    this.draw();
  }

  draw(): void {
    const c = this.canvas;
    hidpi(c, this.opts.aspect ?? 0.56, 200);
    const g = c.getContext("2d");
    if (!g || this.series.length === 0) return;
    const d = dpr();
    const W = c.width;
    const H = c.height;
    const hasRight = this.series.some((s) => s.axis === "right");
    const pad = { l: 58 * d, r: (hasRight ? 58 : 14) * d, t: 14 * d, b: 40 * d };
    const pw = W - pad.l - pad.r;
    const ph = H - pad.t - pad.b;
    g.fillStyle = "#07080b";
    g.fillRect(0, 0, W, H);
    const [x0, x1] = this.opts.xRange ?? extent(this.series, "x");
    let [y0, y1] = this.opts.yRange ?? extent(this.series, "y", "left");
    if (!this.opts.yRange) {
      const m = (y1 - y0) * 0.06;
      y0 -= m;
      y1 += m;
    }
    let [r0, r1] = this.opts.y2Range ?? (hasRight ? extent(this.series, "y", "right") : [0, 1]);
    const xl = this.opts.xLog;
    const X = (v: number): number => pad.l + (xl ? Math.log(v / x0) / Math.log(x1 / x0) : (v - x0) / (x1 - x0)) * pw;
    const Y = (v: number): number => pad.t + (1 - (v - y0) / (y1 - y0)) * ph;
    const y2l = this.opts.y2Log;
    if (y2l) {
      r0 = Math.max(r0, 1e-30);
    }
    const Y2 = (v: number): number => pad.t + (1 - (y2l ? Math.log(v / r0) / Math.log(r1 / r0) : (v - r0) / (r1 - r0))) * ph;

    g.font = `${10.5 * d}px "JetBrains Mono", monospace`;
    g.textBaseline = "middle";
    for (const b of this.opts.bands ?? []) {
      g.fillStyle = "rgba(217,169,79,0.07)";
      g.fillRect(X(b.from), pad.t, X(b.to) - X(b.from), ph);
    }
    g.strokeStyle = "rgba(236,226,205,0.08)";
    g.fillStyle = "rgba(201,194,179,0.75)";
    g.lineWidth = d;
    const yf = this.opts.yFormat ?? ((v: number) => fmt(v));
    for (const v of niceTicks(y0, y1)) {
      const y = Y(v);
      g.beginPath();
      g.moveTo(pad.l, y);
      g.lineTo(pad.l + pw, y);
      g.stroke();
      g.textAlign = "right";
      g.fillText(yf(v), pad.l - 6 * d, y);
    }
    if (this.opts.zeroLine && y0 < 0 && y1 > 0) {
      g.strokeStyle = "rgba(236,226,205,0.22)";
      g.beginPath();
      g.moveTo(pad.l, Y(0));
      g.lineTo(pad.l + pw, Y(0));
      g.stroke();
      g.strokeStyle = "rgba(236,226,205,0.08)";
    }
    if (hasRight) {
      const f2 = this.opts.y2Format ?? ((v: number) => fmt(v));
      const ticks = y2l ? (() => { const out: number[] = []; for (let e = Math.ceil(Math.log10(r0)); e <= Math.floor(Math.log10(r1)); e++) out.push(10 ** e); return out; })() : niceTicks(r0, r1);
      g.textAlign = "left";
      for (const v of ticks) g.fillText(f2(v), pad.l + pw + 6 * d, Y2(v));
    }
    const xf = this.opts.xFormat ?? ((v: number) => fmt(v));
    const xticks = xl ? [10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000].filter((v) => v >= x0 && v <= x1) : niceTicks(x0, x1, 6);
    g.textAlign = "center";
    g.textBaseline = "top";
    for (const v of xticks) {
      const x = X(v);
      g.beginPath();
      g.moveTo(x, pad.t);
      g.lineTo(x, pad.t + ph);
      g.stroke();
      g.fillText(xf(v), x, pad.t + ph + 6 * d);
    }
    g.fillStyle = "rgba(143,137,125,0.95)";
    g.fillText(this.opts.xLabel, pad.l + pw / 2, pad.t + ph + 22 * d);
    g.save();
    g.translate(13 * d, pad.t + ph / 2);
    g.rotate(-Math.PI / 2);
    g.textBaseline = "middle";
    g.fillText(this.opts.yLabel, 0, 0);
    g.restore();
    if (hasRight && this.opts.y2Label) {
      g.save();
      g.translate(W - 10 * d, pad.t + ph / 2);
      g.rotate(Math.PI / 2);
      g.textBaseline = "middle";
      g.fillText(this.opts.y2Label, 0, 0);
      g.restore();
    }
    g.save();
    g.beginPath();
    g.rect(pad.l, pad.t, pw, ph);
    g.clip();
    for (const s of this.series) {
      const right = s.axis === "right";
      g.strokeStyle = s.color;
      g.fillStyle = s.color;
      g.lineWidth = (s.width ?? 1.6) * d;
      g.setLineDash((s.dash ?? []).map((v) => v * d));
      g.beginPath();
      let first = true;
      for (let i = 0; i < s.x.length; i++) {
        const xv = s.x[i]!;
        const yv = s.y[i]!;
        if (!Number.isFinite(xv) || !Number.isFinite(yv) || (xl && xv <= 0)) continue;
        const px = X(xv);
        const py = right ? Y2(Math.max(yv, y2l ? r0 : -Infinity)) : Y(yv);
        if (s.scatter) {
          g.moveTo(px + 2 * d, py);
          g.arc(px, py, 2 * d, 0, Math.PI * 2);
        } else if (first) {
          g.moveTo(px, py);
          first = false;
        } else g.lineTo(px, py);
      }
      if (s.scatter) g.fill();
      else g.stroke();
    }
    g.restore();
    g.setLineDash([]);
    if (this.hoverX !== null) {
      const hx = this.hoverX * d;
      if (hx >= pad.l && hx <= pad.l + pw) {
        g.strokeStyle = "rgba(241,234,219,0.4)";
        g.beginPath();
        g.moveTo(hx, pad.t);
        g.lineTo(hx, pad.t + ph);
        g.stroke();
        const frac = (hx - pad.l) / pw;
        const xv = xl ? x0 * Math.pow(x1 / x0, frac) : x0 + frac * (x1 - x0);
        const lines = [`${this.opts.xLabel.split(" [")[0]}: ${xf(xv)}`];
        for (const s of this.series) {
          let best = -1;
          let bd = Infinity;
          for (let i = 0; i < s.x.length; i++) {
            const dd = Math.abs(s.x[i]! - xv);
            if (dd < bd) {
              bd = dd;
              best = i;
            }
          }
          if (best >= 0) lines.push(`${s.label}: ${(s.axis === "right" ? this.opts.y2Format ?? fmt : yf)(s.y[best]!)}`);
        }
        this.tip.textContent = lines.join("\n");
        this.tip.style.whiteSpace = "pre";
        this.tip.style.opacity = "1";
        const left = this.hoverX + 14;
        const cw = this.canvas.clientWidth;
        this.tip.style.left = `${left > cw - 180 ? this.hoverX - 190 : left}px`;
        this.tip.style.top = "10px";
      }
    }
  }
}
