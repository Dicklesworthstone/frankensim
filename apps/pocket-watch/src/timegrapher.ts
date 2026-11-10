// Timegrapher (Witschi-style beat-deviation chart) and the 12-node material/phase ledger.
//
// Beat sources, best first:
//   1. "impacts"   – engine lock/drop impacts from `timegrapher_drain()` (lab time, side, energy)
//   2. "crossings" – balance passing the line of centres, found by cubic Hermite interpolation of the
//                    snapshot (theta, omega) between frames; this is the oscillator's own timing
//   3. "readout"   – no beat stream yet: show the engine's rate / beat-error readouts only
//
// Each beat's deviation from the nominal 0.2 s grid (18,000 vph) is plotted against time. A least
// squares fit dev = a + b t + c s (s = -1/2 tick, +1/2 tock) gives rate = -86400 b and
// beat error = |c|, exactly how an acoustic timegrapher reads a movement.

import { WATCH_NODES_METADATA, type WatchSnapshot } from "./snapshotLayout";

export const NOMINAL_BEAT_S = 0.2;
const WINDOW_S = 30;
const MAX_POINTS = 360;
const BOUNCE_GATE_S = 0.06;
const GAP_RESET_S = 1.6;

export type BeatSource = "impacts" | "crossings" | "readout";

interface BeatPoint {
  readonly t: number;
  readonly devS: number;
  readonly side: 0 | 1;
  readonly energy: number;
}

export interface TimegrapherReadings {
  readonly source: BeatSource;
  readonly rateSDay: number | null;
  readonly beatErrorMs: number | null;
  readonly beats: number;
  readonly spanS: number;
  /** No beat for more than a second of lab time: the balance has stopped. */
  readonly stopped: boolean;
}

export class HorologicalTimegrapher {
  private points: BeatPoint[] = [];
  private t0 = 0;
  private k = 0;
  private lastBeatT = -Infinity;
  private hasRef = false;
  private sourceKind: BeatSource = "readout";
  private prevT = Number.NaN;
  private prevTheta = 0;
  private prevOmega = 0;
  private fit: { a: number; b: number; c: number; tMean: number } | null = null;
  private fitDirty = true;
  private ctx2d: CanvasRenderingContext2D | null = null;
  private cssW = 0;
  private cssH = 0;
  private dpr = 0;
  private lastLabT = 0;

  public reset(): void {
    this.points = [];
    this.hasRef = false;
    this.k = 0;
    this.lastBeatT = -Infinity;
    this.prevT = Number.NaN;
    this.fit = null;
    this.fitDirty = true;
  }

  public source(): BeatSource {
    return this.sourceKind;
  }

  /** Flat triples [lab_time_s, side (0 tick / 1 tock), impact_energy_nJ] from the engine. */
  public ingestImpacts(flat: Float64Array): void {
    if (this.sourceKind !== "impacts") {
      this.sourceKind = "impacts";
      this.reset();
    }
    for (let i = 0; i + 2 < flat.length; i += 3) {
      const side = (flat[i + 1] ?? 0) >= 0.5 ? 1 : 0;
      this.addBeat(flat[i] ?? 0, side, flat[i + 2] ?? 0);
    }
  }

  /**
   * Fallback beat detector: the balance crossing theta = 0 (line of centres, where the impulse pin
   * meets the fork). Between two frames the motion is interpolated with a cubic Hermite in
   * (theta, omega), which is accurate to well under a millisecond at 2.5 Hz.
   */
  public ingestBalance(labT: number, theta: number, omega: number): void {
    if (this.sourceKind === "impacts") return;
    const dt = labT - this.prevT;
    if (!(dt > 0)) {
      if (Number.isNaN(this.prevT) || dt < 0) {
        this.prevT = labT;
        this.prevTheta = theta;
        this.prevOmega = omega;
      }
      return;
    }
    if (dt <= 0.06 && Math.sign(theta) !== Math.sign(this.prevTheta) && this.prevTheta !== 0) {
      const th0 = this.prevTheta;
      const th1 = theta;
      const m0 = this.prevOmega * dt;
      const m1 = omega * dt;
      const p = (s: number): number => {
        const s2 = s * s;
        const s3 = s2 * s;
        return (2 * s3 - 3 * s2 + 1) * th0 + (s3 - 2 * s2 + s) * m0 + (-2 * s3 + 3 * s2) * th1 + (s3 - s2) * m1;
      };
      let lo = 0;
      let hi = 1;
      const sLo = Math.sign(th0);
      for (let it = 0; it < 30; it++) {
        const mid = 0.5 * (lo + hi);
        if (Math.sign(p(mid)) === sLo) lo = mid;
        else hi = mid;
      }
      const s = 0.5 * (lo + hi);
      const tCross = this.prevT + s * dt;
      // Swinging positive = tick, negative = tock.
      const side: 0 | 1 = th1 > th0 ? 0 : 1;
      if (this.sourceKind !== "crossings") {
        this.sourceKind = "crossings";
        this.points = [];
        this.hasRef = false;
      }
      this.addBeat(tCross, side, 0);
    }
    this.prevT = labT;
    this.prevTheta = theta;
    this.prevOmega = omega;
  }

  private addBeat(t: number, side: 0 | 1, energy: number): void {
    if (t - this.lastBeatT < BOUNCE_GATE_S) return; // pallet bounce / double strike
    if (!this.hasRef || t - this.lastBeatT > GAP_RESET_S) {
      this.points = [];
      this.t0 = t;
      this.k = 0;
      this.hasRef = true;
    } else {
      this.k += Math.max(1, Math.round((t - this.lastBeatT) / NOMINAL_BEAT_S));
    }
    this.lastBeatT = t;
    this.points.push({ t, devS: t - this.t0 - this.k * NOMINAL_BEAT_S, side, energy });
    if (this.points.length > MAX_POINTS) this.points.splice(0, this.points.length - MAX_POINTS);
    this.fitDirty = true;
  }

  private windowPoints(): BeatPoint[] {
    const n = this.points.length;
    if (n === 0) return this.points;
    const tEnd = this.points[n - 1]!.t;
    let i = n - 1;
    while (i > 0 && tEnd - this.points[i - 1]!.t <= WINDOW_S) i--;
    return this.points.slice(i);
  }

  private solveFit(): void {
    this.fitDirty = false;
    const pts = this.windowPoints();
    if (pts.length < 12) {
      this.fit = null;
      return;
    }
    let tMean = 0;
    for (const p of pts) tMean += p.t;
    tMean /= pts.length;
    // Normal equations for [1, t', s] with t' = t - tMean.
    let s11 = 0, s1t = 0, s1s = 0, stt = 0, sts = 0, sss = 0;
    let y1 = 0, yt = 0, ys = 0;
    for (const p of pts) {
      const tt = p.t - tMean;
      const ss = p.side === 1 ? 0.5 : -0.5;
      s11 += 1; s1t += tt; s1s += ss; stt += tt * tt; sts += tt * ss; sss += ss * ss;
      y1 += p.devS; yt += tt * p.devS; ys += ss * p.devS;
    }
    const M = [
      [s11, s1t, s1s],
      [s1t, stt, sts],
      [s1s, sts, sss],
    ];
    const det =
      M[0]![0]! * (M[1]![1]! * M[2]![2]! - M[1]![2]! * M[2]![1]!) -
      M[0]![1]! * (M[1]![0]! * M[2]![2]! - M[1]![2]! * M[2]![0]!) +
      M[0]![2]! * (M[1]![0]! * M[2]![1]! - M[1]![1]! * M[2]![0]!);
    if (Math.abs(det) < 1e-18) {
      this.fit = null;
      return;
    }
    const solve = (col: number): number => {
      const R = M.map((row) => row.slice());
      R[0]![col] = y1;
      R[1]![col] = yt;
      R[2]![col] = ys;
      return (
        (R[0]![0]! * (R[1]![1]! * R[2]![2]! - R[1]![2]! * R[2]![1]!) -
          R[0]![1]! * (R[1]![0]! * R[2]![2]! - R[1]![2]! * R[2]![0]!) +
          R[0]![2]! * (R[1]![0]! * R[2]![1]! - R[1]![1]! * R[2]![0]!)) /
        det
      );
    };
    this.fit = { a: solve(0), b: solve(1), c: solve(2), tMean };
  }

  public readings(labT: number): TimegrapherReadings {
    if (this.fitDirty) this.solveFit();
    const pts = this.windowPoints();
    const span = pts.length > 1 ? pts[pts.length - 1]!.t - pts[0]!.t : 0;
    const stopped = this.hasRef && labT - this.lastBeatT > 1.0;
    if (!this.fit || stopped) {
      return {
        source: this.points.length ? this.sourceKind : "readout",
        rateSDay: null,
        beatErrorMs: null,
        beats: pts.length,
        spanS: span,
        stopped,
      };
    }
    return {
      source: this.sourceKind,
      rateSDay: -this.fit.b * 86_400,
      beatErrorMs: Math.abs(this.fit.c) * 1000,
      beats: pts.length,
      spanS: span,
      stopped: false,
    };
  }

  private ensureCanvas(canvas: HTMLCanvasElement): CanvasRenderingContext2D | null {
    const dpr = Math.min(2, window.devicePixelRatio || 1);
    const w = canvas.clientWidth;
    const h = canvas.clientHeight;
    if (w < 2 || h < 2) return null;
    if (!this.ctx2d || w !== this.cssW || h !== this.cssH || dpr !== this.dpr) {
      canvas.width = Math.round(w * dpr);
      canvas.height = Math.round(h * dpr);
      this.cssW = w;
      this.cssH = h;
      this.dpr = dpr;
      this.ctx2d = canvas.getContext("2d");
      this.ctx2d?.setTransform(dpr, 0, 0, dpr, 0, 0);
    }
    return this.ctx2d;
  }

  public draw(canvas: HTMLCanvasElement, snap: WatchSnapshot): void {
    const ctx = this.ensureCanvas(canvas);
    if (!ctx) return;
    if (this.fitDirty) this.solveFit();
    this.lastLabT = snap.labTimeS;
    const W = this.cssW;
    const H = this.cssH;
    const padL = 40;
    const padR = 8;
    const padT = 8;
    const padB = 18;
    const plotW = W - padL - padR;
    const plotH = H - padT - padB;
    ctx.clearRect(0, 0, W, H);
    ctx.fillStyle = "#070a10";
    ctx.fillRect(0, 0, W, H);

    const pts = this.windowPoints();
    const tEnd = Math.max(pts.length ? pts[pts.length - 1]!.t : 0, this.lastLabT);
    const tStart = tEnd - WINDOW_S;
    // Centre on the window mean so a drifting trace stays on screen; the slope is the rate.
    let mean = 0;
    for (const p of pts) mean += p.devS;
    mean = pts.length ? mean / pts.length : 0;
    let maxAbs = 0;
    for (const p of pts) maxAbs = Math.max(maxAbs, Math.abs(p.devS - mean));
    const steps = [1, 2, 5, 10, 20, 50, 100, 200];
    let scaleMs = 0;
    for (const s of steps) {
      if (maxAbs * 1000 * 1.15 <= s) {
        scaleMs = s;
        break;
      }
    }
    // Beyond ±200 ms in the window the trace wraps modulo one beat, as on a paper-strip
    // timegrapher whose drum turns once per nominal beat: a steep rate becomes parallel bands
    // instead of a line pinned to the chart edge. The reference is the run's first beat, so the
    // bands only scroll sideways.
    const wrap = scaleMs === 0;
    const beatMs = NOMINAL_BEAT_S * 1000;
    if (wrap) scaleMs = beatMs / 2;
    const centre = wrap ? 0 : mean;
    const devMs = (devS: number): number => {
      const d = (devS - centre) * 1000;
      if (!wrap) return d;
      return ((((d + scaleMs) % beatMs) + beatMs) % beatMs) - scaleMs;
    };
    const yOfMs = (dMs: number): number => padT + plotH / 2 - (dMs / scaleMs) * (plotH / 2);
    const yOf = (devS: number): number => yOfMs(devMs(devS));
    const xOf = (t: number): number => padL + ((t - tStart) / WINDOW_S) * plotW;

    // Grid: centre line, ±scale/2, ±scale; 5 s verticals.
    ctx.lineWidth = 1;
    ctx.strokeStyle = "rgba(231, 196, 122, 0.08)";
    ctx.beginPath();
    for (let q = -2; q <= 2; q++) {
      const y = Math.round(padT + plotH / 2 - (q / 2) * (plotH / 2)) + 0.5;
      ctx.moveTo(padL, y);
      ctx.lineTo(W - padR, y);
    }
    for (let s = 0; s <= WINDOW_S; s += 5) {
      const x = Math.round(padL + (s / WINDOW_S) * plotW) + 0.5;
      ctx.moveTo(x, padT);
      ctx.lineTo(x, padT + plotH);
    }
    ctx.stroke();
    ctx.strokeStyle = "rgba(231, 196, 122, 0.22)";
    ctx.beginPath();
    ctx.moveTo(padL, Math.round(padT + plotH / 2) + 0.5);
    ctx.lineTo(W - padR, Math.round(padT + plotH / 2) + 0.5);
    ctx.stroke();

    ctx.fillStyle = "rgba(200, 196, 186, 0.55)";
    ctx.font = '500 9.5px "JetBrains Mono", ui-monospace, monospace';
    ctx.textAlign = "right";
    ctx.textBaseline = "middle";
    ctx.fillText(`+${scaleMs}`, padL - 6, padT + 5);
    ctx.fillText("0 ms", padL - 6, padT + plotH / 2);
    ctx.fillText(`−${scaleMs}`, padL - 6, padT + plotH - 5);
    ctx.textAlign = "left";
    ctx.textBaseline = "alphabetic";

    // Fitted tick / tock lines (their separation is the beat error, their slope the rate).
    if (this.fit && pts.length > 1) {
      const f = this.fit;
      for (const side of [0, 1] as const) {
        const s = side === 1 ? 0.5 : -0.5;
        const ta = Math.max(tStart, pts[0]!.t);
        const tb = tEnd;
        ctx.strokeStyle = side === 0 ? "rgba(232, 199, 125, 0.45)" : "rgba(111, 240, 192, 0.45)";
        ctx.lineWidth = 1;
        ctx.beginPath();
        const devAt = (t: number): number => f.a + f.b * (t - f.tMean) + f.c * s;
        if (!wrap) {
          ctx.moveTo(xOf(ta), yOf(devAt(ta)));
          ctx.lineTo(xOf(tb), yOf(devAt(tb)));
        } else {
          // Split the straight fit at every wrap so each band is drawn edge to edge.
          const slopeMs = f.b * 1000;
          let t = ta;
          let guard = 0;
          while (t < tb && guard++ < 400) {
            const d = devMs(devAt(t));
            const toEdge = slopeMs > 0 ? (scaleMs - d) / slopeMs : slopeMs < 0 ? (-scaleMs - d) / slopeMs : Infinity;
            const tNext = Math.min(tb, t + Math.max(toEdge, 1e-6));
            ctx.moveTo(xOf(t), yOfMs(d));
            ctx.lineTo(xOf(tNext), yOfMs(d + slopeMs * (tNext - t)));
            t = tNext + 1e-6;
          }
        }
        ctx.stroke();
      }
    }

    // Beat dots: tick = filled gold disc, tock = mint ring (shape differs, not colour alone).
    for (const p of pts) {
      const x = xOf(p.t);
      if (x < padL - 2) continue;
      const y = Math.max(padT + 2, Math.min(padT + plotH - 2, yOf(p.devS)));
      ctx.beginPath();
      ctx.arc(x, y, 2.1, 0, Math.PI * 2);
      if (p.side === 0) {
        ctx.fillStyle = "#e8c77d";
        ctx.fill();
      } else {
        ctx.strokeStyle = "#6ff0c0";
        ctx.lineWidth = 1.3;
        ctx.stroke();
      }
    }

    ctx.fillStyle = "rgba(200, 196, 186, 0.6)";
    ctx.font = '500 9.5px "JetBrains Mono", ui-monospace, monospace';
    const silent = pts.length === 0 || tEnd - (pts[pts.length - 1]?.t ?? -1e9) > 1.0;
    let note: string;
    if (snap.statusCode >= 2 && silent) note = "NO BEAT — ESCAPEMENT STOPPED";
    else if (this.sourceKind === "impacts") note = `● tick  ○ tock · lock impacts · 30 s${wrap ? " · wrapped" : ""}`;
    else if (this.sourceKind === "crossings") note = `● tick  ○ tock · line of centres · 30 s${wrap ? " · wrapped" : ""}`;
    else note = "waiting for beats…";
    ctx.fillText(note, padL, H - 5);
  }
}

type CellState = "hot" | "warm" | "cold" | "nominal" | "muted" | "molten" | "frozen" | "gold" | "ice";

function setCell(cell: HTMLTableCellElement, text: string, state: CellState): void {
  if (cell.textContent !== text) cell.textContent = text;
  if (cell.dataset.state !== state) cell.dataset.state = state;
}

export function renderNodeLedgerTable(tbody: HTMLElement, snap: WatchSnapshot): void {
  // Build rows once, then update cells in place.
  if (tbody.children.length !== WATCH_NODES_METADATA.length) {
    tbody.replaceChildren();
    for (const meta of WATCH_NODES_METADATA) {
      const tr = document.createElement("tr");
      const th = document.createElement("th");
      th.scope = "row";
      th.title = `${meta.fullName} · ${meta.materialLabel} · ${meta.matdbPack}`;
      th.textContent = meta.shortName;
      tr.append(th, document.createElement("td"), document.createElement("td"), document.createElement("td"));
      tbody.appendChild(tr);
    }
  }

  for (let i = 0; i < WATCH_NODES_METADATA.length; i++) {
    const row = tbody.children[i] as HTMLTableRowElement | undefined;
    const node = snap.nodes[i];
    if (!row || !node) continue;
    const cTemp = row.children[1] as HTMLTableCellElement;
    const cMelt = row.children[2] as HTMLTableCellElement;
    const cFluid = row.children[3] as HTMLTableCellElement;

    const tempC = node.temperatureK - 273.15;
    setCell(
      cTemp,
      `${node.temperatureK.toFixed(0)} K · ${tempC >= 0 ? "+" : "−"}${Math.abs(tempC).toFixed(0)}°C`,
      node.temperatureK > 950 ? "hot" : node.temperatureK < 210 ? "cold" : node.temperatureK > 330 ? "warm" : "nominal",
    );

    if (i <= 3) {
      if (node.meltFraction > 0.005 || node.plasticSlump > 0.01) {
        setCell(
          cMelt,
          `φ ${(node.meltFraction * 100).toFixed(0)}% · slump ${(node.plasticSlump * 100).toFixed(0)}%`,
          node.meltFraction > 0.01 ? "molten" : "gold",
        );
      } else {
        setCell(cMelt, "Solid Au", "muted");
      }
    } else if (node.moltenGoldG > 0.001) {
      setCell(cMelt, `Liquid Au ${node.moltenGoldG.toFixed(2)} g`, "molten");
    } else if (node.frozenGoldG > 0.001) {
      setCell(cMelt, `Frozen Au ${node.frozenGoldG.toFixed(2)} g`, "frozen");
    } else {
      setCell(cMelt, "Clean", "muted");
    }

    const parts: string[] = [];
    if (node.oilFilmUm > 0.05) parts.push(`${node.oilFilmUm.toFixed(1)} µm oil`);
    if (node.frostThicknessUm > 0.5) parts.push(`${node.frostThicknessUm.toFixed(0)} µm ice`);
    if (node.ln2FilmG > 0.01) parts.push("LN₂ film");
    setCell(
      cFluid,
      parts.length > 0 ? parts.join(" · ") : "—",
      node.frostThicknessUm > 1.0 ? "ice" : node.oilFilmUm > 6.0 ? "gold" : "muted",
    );
  }
}
