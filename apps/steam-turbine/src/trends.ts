// Trend recorder, event log and CSV export for the turbine page.
//
// Every fixed engine step (20 ms of simulated time) is recorded into a fine ring (the last 60 s)
// and, every 0.5 s of simulated time, into a coarse ring (the last 2 h). Values are the decoded
// snapshot's fields (snapshotLayout.ts holds the one slot map); nothing here computes physics.
// The event log records the status, alarm-bit and breaker changes the engine reports, plus the
// operator's actions.

import { ALARM_FLAG_TEXT, STATUS_LABELS } from "./snapshotLayout.ts";
import type { TurbineSnapshot } from "./snapshotLayout.ts";

/** One recordable telemetry channel: a named field of the decoded snapshot, optionally scaled. */
export interface TrendChannel {
  id: string;
  label: string;
  unit: string;
  digits: number;
  key: NumericSnapshotKey;
  scale?: number;
}

type NumericSnapshotKey = { [K in keyof TurbineSnapshot]: TurbineSnapshot[K] extends number ? K : never }[keyof TurbineSnapshot];

const ch = (id: string, label: string, unit: string, digits: number, key: NumericSnapshotKey, scale?: number): TrendChannel =>
  scale === undefined ? { id, label, unit, digits, key } : { id, label, unit, digits, key, scale };

export const TREND_CHANNELS: readonly TrendChannel[] = [
  ch("mw", "Electric power", "MW", 2, "elecMw"),
  ch("rpm", "Shaft speed", "rpm", 1, "rpm"),
  ch("hp_mw", "HP section power", "MW", 2, "hpMw"),
  ch("iplp_mw", "IP/LP section power", "MW", 2, "iplpMw"),
  ch("q_mvar", "Reactive power", "MVAr", 2, "reactiveMvar"),
  ch("delta", "Load angle", "°", 2, "loadAngleDeg"),
  ch("live_flow", "Live steam flow", "kg/s", 2, "liveFlowKgS"),
  ch("chest_p", "Steam chest pressure", "bar", 2, "chestBar"),
  ch("chest_t", "Steam chest temperature", "°C", 1, "chestC"),
  ch("ext_p", "Extraction header pressure", "bar", 2, "extHeaderBar"),
  ch("ext_flow", "Extraction flow", "kg/s", 2, "extFlowKgS"),
  ch("lp_flow", "LP exhaust flow", "kg/s", 2, "lpFlowKgS"),
  ch("cond_p", "Condenser pressure", "kPa", 2, "condKpa"),
  ch("cond_duty", "Condenser duty", "MW", 2, "condDutyMw"),
  ch("lp_x", "LP exhaust dryness", "–", 4, "lpQualityX"),
  ch("cv", "Control valve lift", "%", 1, "cvLift", 100),
  ch("eff_e", "Electric efficiency", "%", 2, "elecEffPct"),
  ch("eff_chp", "CHP utilisation", "%", 2, "chpEffPct"),
  ch("t_hp_rotor", "HP rotor temperature", "°C", 2, "tHpRotorC"),
  ch("t_hp_casing", "HP casing temperature", "°C", 2, "tHpCasingC"),
  ch("t_ip_rotor", "IP/LP rotor temperature", "°C", 2, "tIpRotorC"),
  ch("t_ip_casing", "IP casing temperature", "°C", 2, "tIpCasingC"),
  ch("t_lp_hood", "LP hood temperature", "°C", 2, "tLpHoodC"),
  ch("diff_exp", "Differential expansion", "mm", 3, "diffExpansionMm"),
  ch("hp_stress", "HP thermal stress", "MPa", 1, "hpThermalStressMpa"),
  ch("vib_front", "Front vibration S(p-p)", "µm", 1, "frontVibUmPp"),
  ch("vib_rear", "Rear vibration S(p-p)", "µm", 1, "rearVibUmPp"),
  ch("re_lambda", "Rotor stability max Re(λ)", "1/s", 2, "spectralAbscissaPerS"),
  ch("t_brg_front", "Front Babbitt temperature", "°C", 2, "tFrontBrgC"),
  ch("t_brg_rear", "Rear Babbitt temperature", "°C", 2, "tRearBrgC"),
  ch("t_oil", "Lube oil temperature", "°C", 2, "tLubeOilC"),
  ch("film_front", "Front min oil film", "µm", 1, "frontFilmUm"),
  ch("thrust", "Net axial thrust", "kN", 1, "axialThrustKn"),
];

export const DEFAULT_LANES: readonly string[] = ["mw", "rpm", "t_hp_rotor", "vib_front"];

/** A fixed-capacity ring of (time, all-channel values) rows. */
export class TrendRing {
  readonly cap: number;
  readonly t: Float64Array;
  readonly v: Float32Array;
  private head = 0;
  n = 0;

  constructor(cap: number) {
    this.cap = cap;
    this.t = new Float64Array(cap);
    this.v = new Float32Array(cap * TREND_CHANNELS.length);
  }

  push(snap: TurbineSnapshot): void {
    const i = this.head;
    this.t[i] = snap.timeS;
    const base = i * TREND_CHANNELS.length;
    for (let c = 0; c < TREND_CHANNELS.length; c++) {
      const chn = TREND_CHANNELS[c]!;
      this.v[base + c] = snap[chn.key] * (chn.scale ?? 1);
    }
    this.head = (i + 1) % this.cap;
    this.n = Math.min(this.n + 1, this.cap);
  }

  clear(): void {
    this.head = 0;
    this.n = 0;
  }

  /** Ring index of the k-th oldest row (k in 0..n). */
  index(k: number): number {
    return (this.head - this.n + k + this.cap) % this.cap;
  }

  time(k: number): number {
    return this.t[this.index(k)]!;
  }

  value(k: number, channel: number): number {
    return this.v[this.index(k) * TREND_CHANNELS.length + channel]!;
  }
}

/** Seconds of simulated time between coarse samples. */
export const COARSE_DT_S = 0.5;

export class TrendRecorder {
  /** Every engine step (20 ms): the last 60 s. */
  readonly fine = new TrendRing(3000);
  /** Every 0.5 s of simulated time: the last 2 h. */
  readonly coarse = new TrendRing(14400);
  private nextCoarseT = 0;

  record(snap: TurbineSnapshot): void {
    const t = snap.timeS;
    // The engine clock went back (a reset the caller did not announce): start a new record.
    if (this.fine.n > 0 && t < this.fine.time(this.fine.n - 1) - 1e-9) this.clear();
    this.fine.push(snap);
    if (t + 1e-9 >= this.nextCoarseT) {
      this.coarse.push(snap);
      this.nextCoarseT = (Math.floor(t / COARSE_DT_S + 1e-9) + 1) * COARSE_DT_S;
    }
  }

  clear(): void {
    this.fine.clear();
    this.coarse.clear();
    this.nextCoarseT = 0;
  }

  /**
   * The ring that holds the last `windowS` seconds at the finest resolution: the fine ring for
   * windows up to its 60 s, and for longer windows while it has not yet wrapped (it then still
   * holds the whole record); otherwise the coarse ring.
   */
  ringFor(windowS: number): TrendRing {
    const fineHoldsAll = this.fine.n < this.fine.cap;
    return windowS <= 60 || fineHoldsAll || this.coarse.n < 2 ? this.fine : this.coarse;
  }
}

export type EventKind = "operator" | "alarm" | "trip" | "info";

export interface TurbineEvent {
  t: number;
  kind: EventKind;
  text: string;
}

const statusText = (code: number): string => STATUS_LABELS[code]?.text ?? `status ${code}`;
const maxOf = (a: number, b: number) => Math.max(a, b);

/**
 * What each alarm bit (snapshotLayout's ALARM_FLAG_TEXT order) logs, and the value it quotes.
 * `prev` is the previous step's snapshot: the governor's vibration trip acts on the S(p-p) the
 * rotor model reported at the end of the previous step, so that is the reading the trip quotes.
 */
const FLAG_EVENTS: ReadonlyArray<{
  bit: number;
  kind: EventKind;
  detail: (s: TurbineSnapshot, prev: TurbineSnapshot) => string;
}> = [
  { bit: 0, kind: "alarm", detail: (s) => `S(p-p) ${maxOf(s.frontVibUmPp, s.rearVibUmPp).toFixed(0)} µm` },
  { bit: 1, kind: "alarm", detail: (s) => `Babbitt ${maxOf(s.tFrontBrgC, s.tRearBrgC).toFixed(1)} °C` },
  { bit: 2, kind: "alarm", detail: (s) => `max Re(λ) ${s.spectralAbscissaPerS.toFixed(2)} 1/s` },
  { bit: 3, kind: "trip", detail: (s) => `${s.rpm.toFixed(0)} rpm` },
  { bit: 4, kind: "trip", detail: (s) => `Babbitt ${maxOf(s.tFrontBrgC, s.tRearBrgC).toFixed(1)} °C` },
  { bit: 5, kind: "trip", detail: (_s, prev) => `S(p-p) ${maxOf(prev.frontVibUmPp, prev.rearVibUmPp).toFixed(0)} µm` },
  { bit: 6, kind: "trip", detail: () => "" },
  { bit: 7, kind: "trip", detail: (s) => `load angle ${s.loadAngleDeg.toFixed(0)}°` },
  { bit: 8, kind: "alarm", detail: (s) => `${(s.deratingFrac * 100).toFixed(0)} % output allowed` },
];

/** Status, alarm and breaker changes reported by the engine, plus operator actions. */
export class EventLog {
  readonly events: TurbineEvent[] = [];
  private lastFlags = 0;
  private lastStatus = -1;
  private lastBreaker = -1;
  private lastT = 0;
  private prev: TurbineSnapshot | null = null;
  private readonly cap = 2000;
  onChange: (() => void) | null = null;

  /** Compare a new snapshot with the previous one and log what changed. */
  observe(s: TurbineSnapshot): void {
    const t = s.timeS;
    if (t < this.lastT - 1e-9) {
      this.restart("Engine clock restarted", s);
      return;
    }
    this.lastT = t;
    const flags = s.alarmFlags | 0;
    const status = s.statusCode | 0;
    const breaker = s.breakerClosed ? 1 : 0;
    const changed = flags ^ this.lastFlags;
    for (const f of FLAG_EVENTS) {
      if (changed & (1 << f.bit)) {
        const on = (flags & (1 << f.bit)) !== 0;
        const detail = on ? f.detail(s, this.prev ?? s) : "";
        const name = ALARM_FLAG_TEXT[f.bit] ?? `bit ${f.bit}`;
        this.add(t, on ? f.kind : "info", `${on ? "Raised" : "Cleared"}: ${name}${detail ? ` (${detail})` : ""}`);
      }
    }
    if (this.lastStatus >= 0 && status !== this.lastStatus) {
      this.add(t, "info", `Status ${statusText(this.lastStatus)} → ${statusText(status)}`);
    }
    if (this.lastBreaker >= 0 && breaker !== this.lastBreaker) {
      this.add(t, "info", `Generator breaker ${breaker ? "closed" : "opened"}`);
    }
    this.lastFlags = flags;
    this.lastStatus = status;
    this.lastBreaker = breaker;
    this.prev = s;
  }

  /** Log an operator action at the current simulated time. */
  operator(text: string): void {
    this.add(this.lastT, "operator", text);
  }

  /** Start a new run: the engine's clock restarted. */
  restart(text: string, s: TurbineSnapshot): void {
    this.events.length = 0;
    this.lastFlags = s.alarmFlags | 0;
    this.lastStatus = s.statusCode | 0;
    this.lastBreaker = s.breakerClosed ? 1 : 0;
    this.lastT = s.timeS;
    this.prev = s;
    this.add(this.lastT, "operator", text);
  }

  private add(t: number, kind: EventKind, text: string): void {
    this.events.push({ t, kind, text });
    if (this.events.length > this.cap) this.events.splice(0, this.events.length - this.cap);
    this.onChange?.();
  }
}

/** `m:ss.s` (or `h:mm:ss`) for a simulated time in seconds. */
export function formatSimTime(t: number): string {
  if (!Number.isFinite(t)) return "–";
  const sign = t < 0 ? "-" : "";
  const a = Math.abs(t);
  const h = Math.floor(a / 3600);
  const m = Math.floor((a % 3600) / 60);
  const sec = a % 60;
  if (h > 0) return `${sign}${h}:${String(m).padStart(2, "0")}:${sec.toFixed(0).padStart(2, "0")}`;
  return `${sign}${m}:${sec.toFixed(1).padStart(4, "0")}`;
}

// ---------------------------------------------------------------------------------------------
// CSV export
// ---------------------------------------------------------------------------------------------

function csvCell(v: number, digits: number): string {
  return Number.isFinite(v) ? v.toFixed(digits) : "";
}

/** All channels of one ring as CSV, oldest row first. */
export function trendCsv(ring: TrendRing): string {
  const head = ["t_s", ...TREND_CHANNELS.map((c) => `${c.id} [${c.unit}]`)].join(",");
  const rows: string[] = [head];
  for (let k = 0; k < ring.n; k++) {
    const cells = [ring.time(k).toFixed(3)];
    for (let c = 0; c < TREND_CHANNELS.length; c++) {
      cells.push(csvCell(ring.value(k, c), TREND_CHANNELS[c]!.digits + 1));
    }
    rows.push(cells.join(","));
  }
  return rows.join("\n") + "\n";
}

export function eventsCsv(log: EventLog): string {
  const esc = (s: string) => `"${s.replace(/"/g, '""')}"`;
  return ["t_s,kind,event", ...log.events.map((e) => `${e.t.toFixed(3)},${e.kind},${esc(e.text)}`)].join("\n") + "\n";
}

export function downloadText(filename: string, text: string, mime = "text/csv"): void {
  const blob = new Blob([text], { type: `${mime};charset=utf-8` });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

// ---------------------------------------------------------------------------------------------
// Strip chart
// ---------------------------------------------------------------------------------------------

export const LANE_COLORS: readonly string[] = ["#7dd3fc", "#fbbf24", "#f472b6", "#86efac"];
/** Height of each lane's label row (name and value left, plotted range right). */
const LABEL_BAND = 14;

/**
 * Draw up to four stacked lanes over the last `windowS` seconds of simulated time (or the whole
 * record when `windowS` is Infinity), sharing one time axis. `hoverFrac` (0..1 across the plot)
 * draws a cursor and reads every lane at that time.
 */
export function drawTrendLanes(
  canvas: HTMLCanvasElement,
  rec: TrendRecorder,
  laneIds: readonly string[],
  windowS: number,
  hoverFrac: number | null,
): void {
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  const cssW = Math.max(160, canvas.clientWidth);
  const cssH = Math.max(80, canvas.clientHeight);
  if (canvas.width !== Math.round(cssW * dpr) || canvas.height !== Math.round(cssH * dpr)) {
    canvas.width = Math.round(cssW * dpr);
    canvas.height = Math.round(cssH * dpr);
  }
  const ctx = canvas.getContext("2d");
  if (!ctx) return;
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, cssW, cssH);

  const ring = rec.ringFor(windowS);
  const padL = 8;
  const padR = 8;
  const axisH = 16;
  const plotW = cssW - padL - padR;
  const lanes = laneIds
    .map((id) => TREND_CHANNELS.findIndex((c) => c.id === id))
    .filter((i) => i >= 0)
    .slice(0, 4);
  ctx.font = '10px "JetBrains Mono", ui-monospace, monospace';
  if (ring.n < 2 || lanes.length === 0) {
    ctx.fillStyle = "rgba(164, 183, 201, 0.8)";
    ctx.fillText(lanes.length === 0 ? "Pick a channel below." : "Recording…", padL + 4, 18);
    return;
  }
  const tEnd = ring.time(ring.n - 1);
  const tStartData = ring.time(0);
  const tStart = Number.isFinite(windowS) ? Math.max(tStartData, tEnd - windowS) : tStartData;
  const span = Math.max(1e-6, tEnd - tStart);
  // First row inside the window (binary search on the monotone ring times).
  let lo = 0;
  let hi = ring.n - 1;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (ring.time(mid) < tStart) lo = mid + 1;
    else hi = mid;
  }
  const k0 = Math.max(0, lo - 1);
  const xOf = (t: number) => padL + ((t - tStart) / span) * plotW;
  const laneGap = 6;
  const laneH = (cssH - axisH - laneGap * (lanes.length - 1)) / lanes.length;
  const tHover = hoverFrac === null ? null : tStart + Math.min(1, Math.max(0, hoverFrac)) * span;
  // The row nearest the cursor (binary search on the monotone times).
  let kHover = ring.n - 1;
  if (tHover !== null) {
    let a = k0;
    let b = ring.n - 1;
    while (a < b) {
      const mid = (a + b) >> 1;
      if (ring.time(mid) < tHover) a = mid + 1;
      else b = mid;
    }
    kHover = a > k0 && Math.abs(ring.time(a - 1) - tHover) < Math.abs(ring.time(a) - tHover) ? a - 1 : a;
  }

  lanes.forEach((ch, li) => {
    const top = li * (laneH + laneGap);
    const color = LANE_COLORS[li % LANE_COLORS.length]!;
    const meta = TREND_CHANNELS[ch]!;
    let vMin = Infinity;
    let vMax = -Infinity;
    for (let k = k0; k < ring.n; k++) {
      const v = ring.value(k, ch);
      if (Number.isFinite(v)) {
        vMin = Math.min(vMin, v);
        vMax = Math.max(vMax, v);
      }
    }
    ctx.fillStyle = "rgba(8, 12, 18, 0.55)";
    ctx.fillRect(padL, top + LABEL_BAND, plotW, laneH - LABEL_BAND);
    ctx.strokeStyle = "rgba(126, 166, 198, 0.16)";
    ctx.strokeRect(padL + 0.5, top + LABEL_BAND + 0.5, plotW - 1, laneH - LABEL_BAND - 1);
    if (!Number.isFinite(vMin)) {
      ctx.fillStyle = "rgba(164, 183, 201, 0.8)";
      ctx.fillText(`${meta.label}: no data`, padL + 1, top + 10);
      return;
    }
    if (vMax - vMin < 1e-9 * Math.max(1, Math.abs(vMax))) {
      const pad = Math.max(Math.abs(vMax) * 0.01, 10 ** -meta.digits);
      vMin -= pad;
      vMax += pad;
    } else {
      const pad = (vMax - vMin) * 0.08;
      vMin -= pad;
      vMax += pad;
    }
    const plotTop = top + LABEL_BAND;
    const plotH = laneH - LABEL_BAND;
    const yOf = (v: number) => plotTop + plotH - 2 - ((v - vMin) / (vMax - vMin)) * (plotH - 4);

    ctx.save();
    ctx.beginPath();
    ctx.rect(padL, plotTop, plotW, plotH);
    ctx.clip();
    ctx.strokeStyle = color;
    ctx.lineWidth = 1.4;
    ctx.beginPath();
    let pen = false;
    // More rows than pixel columns: draw each column's minimum and maximum, so a spike or an
    // oscillation keeps its envelope instead of aliasing away.
    const step = Math.max(1, Math.floor((ring.n - k0) / plotW));
    for (let k = k0; k < ring.n; k += step) {
      let a = Number.NaN;
      let b = Number.NaN;
      for (let j = k; j < Math.min(ring.n, k + step); j++) {
        const w = ring.value(j, ch);
        if (Number.isFinite(w)) {
          a = Number.isFinite(a) ? Math.min(a, w) : w;
          b = Number.isFinite(b) ? Math.max(b, w) : w;
        }
      }
      if (!Number.isFinite(a)) {
        pen = false;
        continue;
      }
      const x = xOf(ring.time(k));
      if (pen) ctx.lineTo(x, yOf(a));
      else ctx.moveTo(x, yOf(a));
      if (b !== a) ctx.lineTo(x, yOf(b));
      pen = true;
    }
    ctx.stroke();
    ctx.restore();

    // Value read-out: the latest value, or the value under the cursor.
    const shownV = ring.value(kHover, ch);
    const shown = Number.isFinite(shownV) ? shownV.toFixed(meta.digits) : "—";
    ctx.fillStyle = color;
    ctx.fillText(`${meta.label} ${shown} ${meta.unit}`, padL + 1, top + 10);
    ctx.fillStyle = "rgba(164, 183, 201, 0.75)";
    const rangeTxt = `${vMin.toFixed(meta.digits)} … ${vMax.toFixed(meta.digits)}`;
    ctx.fillText(rangeTxt, padL + plotW - ctx.measureText(rangeTxt).width, top + 10);
  });

  // Shared time axis: about five ticks in simulated time.
  const axisY = cssH - axisH + 11;
  ctx.fillStyle = "rgba(164, 183, 201, 0.75)";
  const nice = [0.1, 0.2, 0.5, 1, 2, 5, 10, 15, 30, 60, 120, 300, 600, 900, 1800, 3600];
  const tickDt = nice.find((d) => span / d <= 6) ?? 3600;
  for (let tt = Math.ceil(tStart / tickDt) * tickDt; tt <= tEnd + 1e-9; tt += tickDt) {
    const x = xOf(tt);
    ctx.strokeStyle = "rgba(126, 166, 198, 0.12)";
    ctx.beginPath();
    ctx.moveTo(x + 0.5, 0);
    ctx.lineTo(x + 0.5, cssH - axisH);
    ctx.stroke();
    const label = formatSimTime(tt);
    ctx.fillText(label, Math.min(cssW - padR - ctx.measureText(label).width, Math.max(padL, x - 14)), axisY);
  }
  if (tHover !== null) {
    const x = xOf(ring.time(kHover));
    const label = formatSimTime(ring.time(kHover));
    const lw = ctx.measureText(label).width + 8;
    const lx = Math.min(cssW - padR - lw, Math.max(padL, x - lw / 2));
    ctx.fillStyle = "rgba(16, 24, 34, 0.96)";
    ctx.fillRect(lx, cssH - axisH + 1, lw, axisH - 2);
    ctx.fillStyle = "rgba(237, 243, 248, 0.95)";
    ctx.fillText(label, lx + 4, axisY);
    ctx.strokeStyle = "rgba(237, 243, 248, 0.55)";
    ctx.setLineDash([3, 3]);
    ctx.beginPath();
    for (let li = 0; li < lanes.length; li++) {
      const plotTop = li * (laneH + laneGap) + LABEL_BAND;
      ctx.moveTo(x + 0.5, plotTop);
      ctx.lineTo(x + 0.5, plotTop + laneH - LABEL_BAND);
    }
    ctx.stroke();
    ctx.setLineDash([]);
  }
}
