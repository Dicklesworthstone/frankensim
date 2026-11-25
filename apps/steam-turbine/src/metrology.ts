// High-DPI Canvas2D scientific metrology plots for the steam turbine model:
// 1. Exact IAPWS-IF97 Mollier (h-s) steam expansion diagram with saturation dome & moisture lines
// 2. 17-Stage Axial Flowpath Profile (Pressure & Temp, Mach & Wetness x, Blade σ_vM & Stage MW)
// 3. 4-DOF gyroscopic Campbell diagram (fs-mbd::rotordynamics) with 1x / 0.46x rays & criticals
// 4. Real-time front-bearing shaft whirl orbit (x-y) with ISO 7919 & labyrinth rub clearances
// 5. L4 ASCENT fs-dfo CMA-ES + Nelder-Mead optimization convergence sparkline

import {
  CAMPBELL_STRIDE,
  MOLLIER_PAIR_STRIDE,
  type EngineLimits,
  type StageTelemetry,
  type TurbineSnapshot,
} from "./snapshotLayout";

export type StageChartMode = "pt" | "mach_wet" | "stress_mw";

function prepareCanvas(canvas: HTMLCanvasElement, fallbackW: number, fallbackH: number): {
  ctx: CanvasRenderingContext2D | null;
  w: number;
  h: number;
} {
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  const w = canvas.clientWidth || fallbackW;
  const h = canvas.clientHeight || fallbackH;
  const pw = Math.round(w * dpr);
  const ph = Math.round(h * dpr);
  if (canvas.width !== pw || canvas.height !== ph) {
    canvas.width = pw;
    canvas.height = ph;
  }
  const ctx = canvas.getContext("2d");
  if (ctx) ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  return { ctx, w, h };
}

/**
 * Draw the live IAPWS-IF97 Mollier (h-s) diagram.
 * `dome` = flat `[s_f0, h_f0, s_g0, h_g0, ...]` in kJ/(kg·K) & kJ/kg.
 * `path` = flat 5-station pairs `[s0, h0, ..., s4, h4]` (Boiler, HP Chest, 12-bar Ext, Exhaust, Exhaust_is).
 */
export function drawMollierDiagram(
  canvas: HTMLCanvasElement,
  dome: Float64Array,
  path: Float64Array,
  snap: TurbineSnapshot,
): void {
  const { ctx, w, h } = prepareCanvas(canvas, 320, 200);
  if (!ctx) return;

  ctx.fillStyle = "#080c12";
  ctx.fillRect(0, 0, w, h);

  const padL = 38;
  const padR = 12;
  const padT = 14;
  const padB = 24;
  const pw = w - padL - padR;
  const ph = h - padT - padB;

  const sMin = 5.4;
  const sMax = 8.8;
  const hMin = 1950;
  const hMax = 3650;

  const xOfS = (s: number) => padL + ((s - sMin) / (sMax - sMin)) * pw;
  const yOfH = (hv: number) => padT + (1 - (hv - hMin) / (hMax - hMin)) * ph;

  ctx.strokeStyle = "rgba(116, 158, 190, 0.14)";
  ctx.fillStyle = "rgba(156, 180, 202, 0.72)";
  ctx.font = '9px "JetBrains Mono", monospace';
  ctx.lineWidth = 1;
  for (let hv = 2000; hv <= 3600; hv += 400) {
    const y = yOfH(hv);
    ctx.beginPath();
    ctx.moveTo(padL, y);
    ctx.lineTo(padL + pw, y);
    ctx.stroke();
    ctx.fillText(`${hv}`, 4, y + 3);
  }
  for (let sv = 5.5; sv <= 8.5; sv += 0.5) {
    const x = xOfS(sv);
    ctx.beginPath();
    ctx.moveTo(x, padT);
    ctx.lineTo(x, padT + ph);
    ctx.stroke();
    ctx.fillText(sv.toFixed(1), x - 8, h - 8);
  }
  ctx.fillStyle = "rgba(176, 198, 218, 0.7)";
  ctx.fillText("h [kJ/kg]", 4, 10);
  ctx.fillText("s [kJ/(kg·K)]", w - 78, h - 1);

  ctx.save();
  ctx.beginPath();
  ctx.rect(padL, padT, pw, ph);
  ctx.clip();

  const nDome = Math.floor(dome.length / 4);
  if (nDome >= 4) {
    ctx.beginPath();
    for (let i = 0; i < nDome; i++) {
      const sg = dome[i * 4 + 2] ?? 0;
      const hg = dome[i * 4 + 3] ?? 0;
      const x = xOfS(sg);
      const y = yOfH(hg);
      if (i === 0) ctx.moveTo(x, y);
      else ctx.lineTo(x, y);
    }
    ctx.lineTo(padL + pw, padT + ph);
    ctx.lineTo(padL, padT + ph);
    ctx.closePath();
    ctx.fillStyle = "rgba(76, 134, 184, 0.14)";
    ctx.fill();

    ctx.strokeStyle = "#69a8d6";
    ctx.lineWidth = 1.6;
    ctx.beginPath();
    for (let i = 0; i < nDome; i++) {
      const sg = dome[i * 4 + 2] ?? 0;
      const hg = dome[i * 4 + 3] ?? 0;
      const x = xOfS(sg);
      const y = yOfH(hg);
      if (i === 0) ctx.moveTo(x, y);
      else ctx.lineTo(x, y);
    }
    ctx.stroke();

    for (const qVal of [0.90, 0.80]) {
      ctx.strokeStyle = "rgba(105, 168, 214, 0.34)";
      ctx.setLineDash([3, 3]);
      ctx.lineWidth = 1;
      ctx.beginPath();
      for (let i = 0; i < nDome; i++) {
        const sf = dome[i * 4] ?? 0;
        const hf = dome[i * 4 + 1] ?? 0;
        const sg = dome[i * 4 + 2] ?? 0;
        const hg = dome[i * 4 + 3] ?? 0;
        const sx = sf + qVal * (sg - sf);
        const hx = hf + qVal * (hg - hf);
        const x = xOfS(sx);
        const y = yOfH(hx);
        if (i === 0) ctx.moveTo(x, y);
        else ctx.lineTo(x, y);
      }
      ctx.stroke();
      ctx.setLineDash([]);
    }

    ctx.fillStyle = "rgba(128, 186, 228, 0.78)";
    ctx.fillText("x=1.0 sat", xOfS(5.65), yOfH(2825));
    ctx.fillText("x=0.90", xOfS(5.65), yOfH(2525));
    ctx.fillText("x=0.80", xOfS(5.65), yOfH(2255));
  }

  const nPts = Math.floor(path.length / MOLLIER_PAIR_STRIDE);
  const stationLabels = ["Boiler", "HP Chest", "12-bar Ext", "LP Exh"];
  const stationColors = ["#ff7849", "#f5a54a", "#eab663", "#5ad2f4"];

  if (nPts >= 5) {
    const sChest = path[1 * MOLLIER_PAIR_STRIDE] ?? 6.5;
    const hChest = path[1 * MOLLIER_PAIR_STRIDE + 1] ?? 3400;
    const sIs = path[4 * MOLLIER_PAIR_STRIDE] ?? sChest;
    const hIs = path[4 * MOLLIER_PAIR_STRIDE + 1] ?? 2200;
    ctx.strokeStyle = "rgba(255, 255, 255, 0.36)";
    ctx.setLineDash([4, 3]);
    ctx.lineWidth = 1.1;
    ctx.beginPath();
    ctx.moveTo(xOfS(sChest), yOfH(hChest));
    ctx.lineTo(xOfS(sIs), yOfH(hIs));
    ctx.stroke();
    ctx.setLineDash([]);

    ctx.strokeStyle = "#f7c46c";
    ctx.lineWidth = 2.2;
    ctx.beginPath();
    for (let i = 0; i < 4; i++) {
      const sv = path[i * MOLLIER_PAIR_STRIDE] ?? 6.5;
      const hv = path[i * MOLLIER_PAIR_STRIDE + 1] ?? 3000;
      const x = xOfS(sv);
      const y = yOfH(hv);
      if (i === 0) ctx.moveTo(x, y);
      else ctx.lineTo(x, y);
    }
    ctx.stroke();

    for (let i = 0; i < 4; i++) {
      const sv = path[i * MOLLIER_PAIR_STRIDE] ?? 6.5;
      const hv = path[i * MOLLIER_PAIR_STRIDE + 1] ?? 3000;
      const x = xOfS(sv);
      const y = yOfH(hv);
      ctx.fillStyle = stationColors[i] ?? "#fff";
      ctx.beginPath();
      ctx.arc(x, y, 4, 0, Math.PI * 2);
      ctx.fill();
      ctx.strokeStyle = "#080c12";
      ctx.lineWidth = 1;
      ctx.stroke();

      ctx.fillStyle = "rgba(232, 240, 248, 0.9)";
      const label = stationLabels[i] ?? "";
      const dx = i === 0 ? -36 : 6;
      const dy = i === 1 ? -6 : 3;
      ctx.fillText(label, x + dx, y + dy);
    }
  }

  ctx.restore();

  const moistPct = Math.max(0, (1 - snap.lpQualityX) * 100);
  ctx.fillStyle = moistPct > 12 ? "#ff7070" : "#8bd5f7";
  ctx.font = '10px "JetBrains Mono", monospace';
  ctx.fillText(`LP wetness: ${moistPct.toFixed(1)}%`, w - 118, 12);
}

/**
 * Draw the 17-Stage Row-by-Row Axial Flowpath Profile from `stage_meanline.rs`
 * (Curtis Stage 0 + HP 1..7 + IP 8..13 + LP 14..16).
 */
export function drawStageProfileChart(
  canvas: HTMLCanvasElement,
  stages: StageTelemetry[],
  mode: StageChartMode,
): void {
  const { ctx, w, h } = prepareCanvas(canvas, 320, 152);
  if (!ctx) return;

  ctx.fillStyle = "#080c12";
  ctx.fillRect(0, 0, w, h);
  if (stages.length === 0) return;

  const padL = 32;
  const padR = 32;
  const padT = 18;
  const padB = 20;
  const pw = w - padL - padR;
  const ph = h - padT - padB;
  const n = stages.length;
  const barStep = pw / n;
  const barW = Math.max(4, barStep * 0.62);

  // Background section bands: HP (0..7), IP (8..13), LP (14..16)
  ctx.fillStyle = "rgba(249, 115, 22, 0.06)";
  ctx.fillRect(padL, padT, 8 * barStep, ph);
  ctx.fillStyle = "rgba(234, 179, 8, 0.05)";
  ctx.fillRect(padL + 8 * barStep, padT, 6 * barStep, ph);
  ctx.fillStyle = "rgba(56, 189, 248, 0.07)";
  ctx.fillRect(padL + 14 * barStep, padT, 3 * barStep, ph);

  // Horizontal grid
  ctx.strokeStyle = "rgba(116, 158, 190, 0.14)";
  ctx.lineWidth = 1;
  for (let g = 0; g <= 4; g++) {
    const y = padT + (g / 4) * ph;
    ctx.beginPath();
    ctx.moveTo(padL, y);
    ctx.lineTo(padL + pw, y);
    ctx.stroke();
  }

  // Extraction station marker between stage 7 & 8
  const extX = padL + 8 * barStep;
  ctx.strokeStyle = "rgba(250, 204, 21, 0.35)";
  ctx.setLineDash([2, 2]);
  ctx.beginPath();
  ctx.moveTo(extX, padT);
  ctx.lineTo(extX, padT + ph);
  ctx.stroke();
  ctx.setLineDash([]);

  ctx.font = '8.5px "JetBrains Mono", monospace';

  if (mode === "pt") {
    // Bars: Temperature [0..560 C], Line: log10 Pressure [0.05..160 bar]
    const tMax = 560;
    const lpMin = Math.log10(0.05);
    const lpMax = Math.log10(160);

    for (let i = 0; i < n; i++) {
      const s = stages[i]!;
      const x = padL + i * barStep + (barStep - barW) * 0.5;
      const frac = Math.max(0, Math.min(1, s.tOutC / tMax));
      const bh = frac * ph;
      ctx.fillStyle = i < 8 ? "rgba(249, 115, 22, 0.62)" : i < 14 ? "rgba(234, 179, 8, 0.58)" : "rgba(56, 189, 248, 0.62)";
      ctx.fillRect(x, padT + ph - bh, barW, bh);
    }

    ctx.strokeStyle = "#7dd3fc";
    ctx.lineWidth = 1.9;
    ctx.beginPath();
    for (let i = 0; i < n; i++) {
      const s = stages[i]!;
      const cx = padL + (i + 0.5) * barStep;
      const lp = Math.log10(Math.max(0.05, s.pOutBar));
      const y = padT + (1 - Math.max(0, Math.min(1, (lp - lpMin) / (lpMax - lpMin)))) * ph;
      if (i === 0) ctx.moveTo(cx, y);
      else ctx.lineTo(cx, y);
    }
    ctx.stroke();

    ctx.fillStyle = "#fb923c";
    ctx.fillText("560°C", 2, padT + 6);
    ctx.fillText("0°C", 6, padT + ph);
    ctx.fillStyle = "#7dd3fc";
    ctx.fillText("160b", w - 30, padT + 6);
    ctx.fillText("0.05", w - 30, padT + ph);
    ctx.fillStyle = "rgba(220, 232, 244, 0.85)";
    ctx.fillText("Bars: T [°C] · Line: log P [bar]", padL + 2, 11);
  } else if (mode === "mach_wet") {
    // Bars: Stage Max Mach [0..1.4], Line: Steam Quality x [0.80..1.00]
    const mMax = 1.4;
    for (let i = 0; i < n; i++) {
      const s = stages[i]!;
      const x = padL + i * barStep + (barStep - barW) * 0.5;
      const frac = Math.max(0, Math.min(1, s.machMax / mMax));
      const bh = frac * ph;
      ctx.fillStyle = s.machMax > 1.0 ? "rgba(244, 63, 94, 0.76)" : "rgba(45, 212, 191, 0.64)";
      ctx.fillRect(x, padT + ph - bh, barW, bh);
    }

    // Ma = 1.0 sonic reference line
    const ySonic = padT + (1 - 1.0 / mMax) * ph;
    ctx.strokeStyle = "rgba(244, 63, 94, 0.45)";
    ctx.setLineDash([3, 2]);
    ctx.beginPath();
    ctx.moveTo(padL, ySonic);
    ctx.lineTo(padL + pw, ySonic);
    ctx.stroke();
    ctx.setLineDash([]);

    ctx.strokeStyle = "#38bdf8";
    ctx.lineWidth = 1.9;
    ctx.beginPath();
    for (let i = 0; i < n; i++) {
      const s = stages[i]!;
      const cx = padL + (i + 0.5) * barStep;
      const qFrac = Math.max(0, Math.min(1, (s.qualityX - 0.80) / 0.20));
      const y = padT + (1 - qFrac) * ph;
      if (i === 0) ctx.moveTo(cx, y);
      else ctx.lineTo(cx, y);
    }
    ctx.stroke();

    ctx.fillStyle = "#2dd4bf";
    ctx.fillText("1.4M", 2, padT + 6);
    ctx.fillText("0.0", 6, padT + ph);
    ctx.fillStyle = "#38bdf8";
    ctx.fillText("x=1", w - 28, padT + 6);
    ctx.fillText("0.80", w - 30, padT + ph);
    ctx.fillStyle = "rgba(220, 232, 244, 0.85)";
    ctx.fillText("Bars: Mach · Line: Dryness x", padL + 2, 11);
  } else {
    // Bars: Blade Root Von Mises Stress [0..600 MPa], Line: Stage Power [0..10 MW]
    const sMax = 600;
    const pMax = Math.max(6.0, ...stages.map((st) => st.powerMw)) * 1.15;
    for (let i = 0; i < n; i++) {
      const s = stages[i]!;
      const x = padL + i * barStep + (barStep - barW) * 0.5;
      const frac = Math.max(0, Math.min(1, s.rootVmMpa / sMax));
      const bh = frac * ph;
      ctx.fillStyle =
        s.rootVmMpa > 420
          ? "rgba(239, 68, 68, 0.78)"
          : s.rootVmMpa > 260
            ? "rgba(250, 204, 21, 0.72)"
            : "rgba(52, 211, 153, 0.64)";
      ctx.fillRect(x, padT + ph - bh, barW, bh);
    }

    ctx.strokeStyle = "#facc15";
    ctx.lineWidth = 1.9;
    ctx.beginPath();
    for (let i = 0; i < n; i++) {
      const s = stages[i]!;
      const cx = padL + (i + 0.5) * barStep;
      const y = padT + (1 - Math.max(0, Math.min(1, s.powerMw / pMax))) * ph;
      if (i === 0) ctx.moveTo(cx, y);
      else ctx.lineTo(cx, y);
    }
    ctx.stroke();

    ctx.fillStyle = "#f87171";
    ctx.fillText("600", 4, padT + 6);
    ctx.fillText("0MPa", 2, padT + ph);
    ctx.fillStyle = "#facc15";
    ctx.fillText(`${pMax.toFixed(0)}MW`, w - 31, padT + 6);
    ctx.fillText("0MW", w - 28, padT + ph);
    ctx.fillStyle = "rgba(220, 232, 244, 0.85)";
    ctx.fillText("Bars: Root σ_vM · Line: Stage MW", padL + 2, 11);
  }

  // Stage axis labels
  ctx.fillStyle = "rgba(156, 180, 202, 0.78)";
  for (const idx of [0, 3, 7, 11, 14, 16]) {
    const cx = padL + (idx + 0.5) * barStep;
    ctx.fillText(`S${idx + 1}`, cx - 6, h - 5);
  }
}

/**
 * Draw the 4-DOF gyroscopic Campbell diagram (`fs_mbd::rotordynamics`).
 */
export function drawCampbellDiagram(
  canvas: HTMLCanvasElement,
  table: Float64Array,
  snap: TurbineSnapshot,
): void {
  const { ctx, w, h } = prepareCanvas(canvas, 320, 175);
  if (!ctx) return;

  ctx.fillStyle = "#080c12";
  ctx.fillRect(0, 0, w, h);

  const padL = 34;
  const padR = 10;
  const padT = 14;
  const padB = 22;
  const pw = w - padL - padR;
  const ph = h - padT - padB;

  const rpmMax = 4500;
  const fMax = 110;
  const xOfRpm = (r: number) => padL + (Math.max(0, Math.min(rpmMax, r)) / rpmMax) * pw;
  const yOfHz = (f: number) => padT + (1 - Math.max(0, Math.min(fMax, f)) / fMax) * ph;

  ctx.strokeStyle = "rgba(116, 158, 190, 0.14)";
  ctx.fillStyle = "rgba(156, 180, 202, 0.72)";
  ctx.font = '9px "JetBrains Mono", monospace';
  ctx.lineWidth = 1;
  for (let f = 20; f <= 100; f += 20) {
    const y = yOfHz(f);
    ctx.beginPath();
    ctx.moveTo(padL, y);
    ctx.lineTo(padL + pw, y);
    ctx.stroke();
    ctx.fillText(`${f}Hz`, 3, y + 3);
  }
  for (let r = 1000; r <= 4000; r += 1000) {
    const x = xOfRpm(r);
    ctx.beginPath();
    ctx.moveTo(x, padT);
    ctx.lineTo(x, padT + ph);
    ctx.stroke();
    ctx.fillText(`${r}`, x - 12, h - 6);
  }

  ctx.strokeStyle = "rgba(90, 210, 244, 0.45)";
  ctx.setLineDash([4, 3]);
  ctx.beginPath();
  ctx.moveTo(xOfRpm(0), yOfHz(0));
  ctx.lineTo(xOfRpm(rpmMax), yOfHz(rpmMax / 60));
  ctx.stroke();

  ctx.strokeStyle = "rgba(255, 112, 112, 0.45)";
  ctx.beginPath();
  ctx.moveTo(xOfRpm(0), yOfHz(0));
  ctx.lineTo(xOfRpm(rpmMax), yOfHz((rpmMax / 60) * 0.46));
  ctx.stroke();
  ctx.setLineDash([]);

  ctx.fillStyle = "rgba(90, 210, 244, 0.75)";
  ctx.fillText("1×", xOfRpm(4200), yOfHz(73));
  ctx.fillStyle = "rgba(255, 128, 128, 0.75)";
  ctx.fillText("0.46×", xOfRpm(3950), yOfHz(34));

  const n = Math.floor(table.length / CAMPBELL_STRIDE);
  const modeColors = ["#60a5fa", "#34d399", "#f59e0b", "#c084fc"];
  if (n >= 2) {
    for (let m = 0; m < 4; m++) {
      ctx.strokeStyle = modeColors[m] ?? "#fff";
      ctx.lineWidth = m < 2 ? 1.8 : 1.25;
      ctx.beginPath();
      for (let i = 0; i < n; i++) {
        const rpm = table[i * CAMPBELL_STRIDE] ?? 0;
        const hz = table[i * CAMPBELL_STRIDE + 1 + m] ?? 0;
        const x = xOfRpm(rpm);
        const y = yOfHz(hz);
        if (i === 0) ctx.moveTo(x, y);
        else ctx.lineTo(x, y);
      }
      ctx.stroke();
    }
  }

  if (snap.crit1Rpm > 300 && snap.crit1Rpm < rpmMax) {
    const cx = xOfRpm(snap.crit1Rpm);
    const cy = yOfHz(snap.crit1Rpm / 60);
    ctx.fillStyle = "#facc15";
    ctx.beginPath();
    ctx.moveTo(cx, cy - 5);
    ctx.lineTo(cx + 5, cy);
    ctx.lineTo(cx, cy + 5);
    ctx.lineTo(cx - 5, cy);
    ctx.closePath();
    ctx.fill();
    ctx.fillText(`Nc1 ${Math.round(snap.crit1Rpm)}`, Math.max(padL + 2, cx - 26), cy - 8);
  }

  const opX = xOfRpm(snap.rpm);
  ctx.strokeStyle = snap.spectralAbscissaPerS > 0 ? "#ff5252" : "#ffffff";
  ctx.lineWidth = 1.6;
  ctx.beginPath();
  ctx.moveTo(opX, padT);
  ctx.lineTo(opX, padT + ph);
  ctx.stroke();

  ctx.fillStyle = snap.spectralAbscissaPerS > 0 ? "#ff6b6b" : "#e2e8f0";
  ctx.fillText(`${Math.round(snap.rpm)} rpm`, Math.min(w - 58, opX + 4), padT + 10);
}

/**
 * Draw the front-journal orbit (the engine samples it inside its 0.5 ms rotor sub-steps, up to
 * 48 points per revolution) with the ISO 7919 zone circles and the bearing clearance. A circle of
 * radius r is an orbit whose peak-to-peak is 2r, so each ISO peak-to-peak limit is drawn at half.
 */
export function drawOrbitPlot(
  canvas: HTMLCanvasElement,
  orbitFlat: Float32Array,
  snap: TurbineSnapshot,
  limits: EngineLimits,
  bearingClearanceUm = 220,
): void {
  const { ctx, w, h } = prepareCanvas(canvas, 320, 160);
  if (!ctx) return;

  ctx.fillStyle = "#080c12";
  ctx.fillRect(0, 0, w, h);

  const cx = w * 0.36;
  const cy = h * 0.52;
  const radiusPx = Math.min(w * 0.30, h * 0.43);

  const brgClear = Math.max(100, bearingClearanceUm);
  const peakAmp = Math.max(snap.frontVibUmPp * 0.5, snap.rearVibUmPp * 0.5, snap.sealExcursionUm);
  const maxUm = peakAmp > 180 ? Math.max(420, Math.ceil(brgClear * 1.25)) : peakAmp < 45 ? 60 : Math.max(240, Math.ceil(brgClear + 20));
  const scale = radiusPx / maxUm;

  ctx.strokeStyle = "rgba(116, 158, 190, 0.18)";
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(cx - radiusPx, cy);
  ctx.lineTo(cx + radiusPx, cy);
  ctx.moveTo(cx, cy - radiusPx);
  ctx.lineTo(cx, cy + radiusPx);
  ctx.stroke();

  const rings: Array<{ rUm: number; color: string; dash: number[] }> = [
    { rUm: limits.vibZoneAbUm / 2, color: "rgba(52, 211, 153, 0.38)", dash: [2, 2] },
    { rUm: limits.vibAlarmUm / 2, color: "rgba(250, 204, 21, 0.45)", dash: [3, 2] },
    { rUm: limits.vibTripUm / 2, color: "rgba(255, 82, 82, 0.72)", dash: [4, 2] },
    { rUm: brgClear, color: "rgba(96, 165, 250, 0.55)", dash: [] },
  ];
  for (const ring of rings) {
    if (ring.rUm > maxUm * 1.02) continue;
    ctx.strokeStyle = ring.color;
    ctx.setLineDash(ring.dash);
    ctx.beginPath();
    ctx.arc(cx, cy, ring.rUm * scale, 0, Math.PI * 2);
    ctx.stroke();
  }
  ctx.setLineDash([]);

  const n = Math.floor(orbitFlat.length / 2);
  if (n >= 2) {
    for (let i = 1; i < n; i++) {
      const alpha = 0.12 + 0.88 * (i / n);
      const x0 = cx + (orbitFlat[(i - 1) * 2] ?? 0) * scale;
      const y0 = cy - (orbitFlat[(i - 1) * 2 + 1] ?? 0) * scale;
      const x1 = cx + (orbitFlat[i * 2] ?? 0) * scale;
      const y1 = cy - (orbitFlat[i * 2 + 1] ?? 0) * scale;
      ctx.strokeStyle =
        snap.rubForceKn > 0.1
          ? `rgba(255, 92, 92, ${alpha.toFixed(2)})`
          : snap.spectralAbscissaPerS > 0
            ? `rgba(251, 191, 36, ${alpha.toFixed(2)})`
            : `rgba(90, 210, 244, ${alpha.toFixed(2)})`;
      ctx.lineWidth = i > n - 12 ? 2.0 : 1.25;
      ctx.beginPath();
      ctx.moveTo(x0, y0);
      ctx.lineTo(x1, y1);
      ctx.stroke();
    }
    const lx = cx + (orbitFlat[(n - 1) * 2] ?? 0) * scale;
    const ly = cy - (orbitFlat[(n - 1) * 2 + 1] ?? 0) * scale;
    ctx.fillStyle = "#ffffff";
    ctx.beginPath();
    ctx.arc(lx, ly, 3.2, 0, Math.PI * 2);
    ctx.fill();
  }

  const lx = w * 0.67;
  ctx.font = '9.5px "JetBrains Mono", monospace';
  ctx.fillStyle = "rgba(180, 200, 218, 0.85)";
  ctx.fillText(`Scale: ±${maxUm} µm`, lx, 22);
  ctx.fillStyle = "#34d399";
  ctx.fillText(`${limits.vibZoneAbUm.toFixed(0)} µm pp ISO A/B`, lx, 42);
  ctx.fillStyle = "#facc15";
  ctx.fillText(`${limits.vibAlarmUm.toFixed(0)} µm pp Alarm`, lx, 60);
  ctx.fillStyle = "#ff5c66";
  ctx.fillText(`${limits.vibTripUm.toFixed(0)} µm pp Trip`, lx, 78);
  ctx.fillStyle = "#60a5fa";
  ctx.fillText(`${Math.round(brgClear)} µm Clearance`, lx, 96);
  ctx.fillStyle = snap.rubForceKn > 0.1 ? "#ff6b6b" : "#94a3b8";
  ctx.fillText(`Seal: ${snap.sealExcursionUm.toFixed(0)} µm`, lx, 122);
  ctx.fillText(`Rub:  ${snap.rubForceKn.toFixed(1)} kN`, lx, 140);
}

/**
 * Draw the L4 ASCENT `fs-dfo` CMA-ES + Nelder-Mead optimization convergence trace.
 */
export function drawOptConvergenceChart(
  canvas: HTMLCanvasElement,
  historyObj: number[],
): void {
  const { ctx, w, h } = prepareCanvas(canvas, 286, 86);
  if (!ctx) return;

  ctx.fillStyle = "#070b11";
  ctx.fillRect(0, 0, w, h);

  if (historyObj.length < 2) {
    ctx.fillStyle = "rgba(148, 163, 184, 0.65)";
    ctx.font = '9.5px "JetBrains Mono", monospace';
    ctx.fillText("Run fs-dfo CMA-ES + Nelder-Mead to view convergence", 10, h * 0.54);
    return;
  }

  const padL = 8;
  const padR = 8;
  const padT = 14;
  const padB = 14;
  const pw = w - padL - padR;
  const ph = h - padT - padB;

  const minV = Math.min(...historyObj);
  const maxV = Math.max(...historyObj);
  const span = Math.max(1e-6, maxV - minV);

  ctx.strokeStyle = "rgba(52, 211, 153, 0.22)";
  ctx.lineWidth = 1;
  ctx.strokeRect(padL, padT, pw, ph);

  ctx.strokeStyle = "#34d399";
  ctx.lineWidth = 1.9;
  ctx.beginPath();
  for (let i = 0; i < historyObj.length; i++) {
    const x = padL + (i / (historyObj.length - 1)) * pw;
    // Higher (worse) objective at top -> lower (better minV) descends toward bottom!
    const y = padT + (1 - (historyObj[i]! - minV) / span) * ph;
    if (i === 0) ctx.moveTo(x, y);
    else ctx.lineTo(x, y);
  }
  ctx.stroke();

  ctx.fillStyle = "rgba(186, 230, 253, 0.85)";
  ctx.font = '8.5px "JetBrains Mono", monospace';
  ctx.fillText(`Objective f(x) ↓  (${historyObj.length} trace pts)`, padL + 2, 10);
  ctx.fillStyle = "#34d399";
  ctx.fillText(`best=${minV.toFixed(3)}`, w - 82, 10);
}
