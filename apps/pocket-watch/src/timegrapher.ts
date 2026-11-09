// Swiss Horological Timegrapher (Witschi-style acoustic beat strip chart)
// and 12-node Material & Enthalpy Phase Ledger renderer.

import { WATCH_NODES_METADATA, type WatchSnapshot } from "./snapshotLayout";

interface BeatDot {
  readonly beatIndex: number;
  readonly isTick: boolean;
  readonly yOffsetPx: number;
  readonly rateSPerDay: number;
  readonly amplitudeDeg: number;
  readonly jammed: boolean;
}

export class HorologicalTimegrapher {
  private readonly dots: BeatDot[] = [];
  private lastBeatRecorded = 0;
  private cumulativeY = 0;

  public reset(): void {
    this.dots.length = 0;
    this.lastBeatRecorded = 0;
    this.cumulativeY = 0;
  }

  public update(snap: WatchSnapshot, canvas: HTMLCanvasElement): void {
    const beats = Math.floor(snap.totalBeats);
    if (beats > this.lastBeatRecorded) {
      const deltaBeats = Math.min(beats - this.lastBeatRecorded, 8);
      for (let b = 0; b < deltaBeats; b++) {
        const idx = this.lastBeatRecorded + b + 1;
        const isTick = idx % 2 === 1;
        // Effective rate combines measured period drift + instantaneous thermoelastic/magnetic rate
        const effRate = 0.4 * snap.measuredRateSPerDay + 0.6 * snap.theoreticalRateSPerDay;
        // Slope on the strip chart is proportional to rate (s/day)
        this.cumulativeY += (effRate / 300.0) * 0.48;
        if (this.cumulativeY > 42) this.cumulativeY -= 84;
        if (this.cumulativeY < -42) this.cumulativeY += 84;

        // Beat error separates the Tick and Tock traces vertically
        const beatSplit = (isTick ? -1 : 1) * Math.min(28, snap.measuredBeatErrorMs * 2.4 + Math.abs(snap.bFieldEscapementMt) * 0.45);
        this.dots.push({
          beatIndex: idx,
          isTick,
          yOffsetPx: this.cumulativeY + beatSplit,
          rateSPerDay: effRate,
          amplitudeDeg: snap.balanceAmplitudeDeg,
          jammed: snap.statusCode >= 2,
        });
        if (this.dots.length > 180) {
          this.dots.shift();
        }
      }
      this.lastBeatRecorded = beats;
    }

    this.drawChart(canvas, snap);
  }

  private drawChart(canvas: HTMLCanvasElement, snap: WatchSnapshot): void {
    const ctx = canvas.getContext("2d");
    if (!ctx) return;
    const W = canvas.width;
    const H = canvas.height;
    const midY = H / 2;

    ctx.fillStyle = "#06080c";
    ctx.fillRect(0, 0, W, H);

    // Horizontal reference zero-rate line and +/- drift guide rails
    ctx.strokeStyle = "rgba(217, 177, 95, 0.16)";
    ctx.lineWidth = 1;
    ctx.setLineDash([4, 4]);
    ctx.beginPath();
    ctx.moveTo(0, midY);
    ctx.lineTo(W, midY);
    ctx.moveTo(0, midY - 28);
    ctx.lineTo(W, midY - 28);
    ctx.moveTo(0, midY + 28);
    ctx.lineTo(W, midY + 28);
    ctx.stroke();
    ctx.setLineDash([]);

    // Plot the dual acoustic Tick (gold) and Tock (emerald) trace dots
    const n = this.dots.length;
    for (let i = 0; i < n; i++) {
      const d = this.dots[i]!;
      const x = ((i + 1) / 180) * (W - 16) + 8;
      const y = midY - d.yOffsetPx;
      ctx.fillStyle = d.isTick ? "#f4d284" : "#5ef2b8";
      ctx.beginPath();
      ctx.arc(x, Math.max(8, Math.min(H - 8, y)), 2.1, 0, Math.PI * 2);
      ctx.fill();
    }

    ctx.fillStyle = "rgba(152, 148, 138, 0.85)";
    ctx.font = '9.5px "JetBrains Mono", monospace';
    const statusNote =
      snap.statusCode >= 2
        ? "NO ACOUSTIC BEAT — ESCAPEMENT HALTED"
        : "● TICK (ENTRY RUBY)   ● TOCK (EXIT RUBY)";
    ctx.fillText(statusNote, 8, H - 6);
  }
}

export function renderNodeLedgerTable(tbody: HTMLElement, snap: WatchSnapshot): void {
  // Lazily build rows once, then update cells in place for zero DOM thrash
  if (tbody.children.length !== WATCH_NODES_METADATA.length) {
    tbody.innerHTML = "";
    for (const meta of WATCH_NODES_METADATA) {
      const tr = document.createElement("tr");
      tr.innerHTML = `<td title="${meta.fullName} (${meta.matdbPack})">${meta.shortName}</td><td></td><td></td><td></td>`;
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
    cTemp.textContent = `${node.temperatureK.toFixed(0)}K (${tempC >= 0 ? "+" : ""}${tempC.toFixed(0)}°C)`;
    if (node.temperatureK > 950) {
      cTemp.style.color = "#ff7847";
    } else if (node.temperatureK < 210) {
      cTemp.style.color = "#56d8ff";
    } else if (node.temperatureK > 330) {
      cTemp.style.color = "#ffc168";
    } else {
      cTemp.style.color = "#ece8df";
    }

    if (i <= 3) {
      // Gold Case node: show liquid fraction and plastic slump %
      if (node.meltFraction > 0.005 || node.plasticSlump > 0.01) {
        cMelt.textContent = `φ=${(node.meltFraction * 100).toFixed(0)}% slp=${(node.plasticSlump * 100).toFixed(0)}%`;
        cMelt.style.color = node.meltFraction > 0.01 ? "#ff6b2b" : "#f4d284";
      } else {
        cMelt.textContent = "Solid Au";
        cMelt.style.color = "#98948a";
      }
    } else {
      // Inner mechanism node: show liquid / frozen gold fouling mass
      if (node.moltenGoldG > 0.001) {
        cMelt.textContent = `Liq Au ${node.moltenGoldG.toFixed(2)}g`;
        cMelt.style.color = "#ff6b2b";
      } else if (node.frozenGoldG > 0.001) {
        cMelt.textContent = `Solid Au ${node.frozenGoldG.toFixed(2)}g`;
        cMelt.style.color = "#ff7a8e";
      } else {
        cMelt.textContent = "Clean";
        cMelt.style.color = "#98948a";
      }
    }

    const parts: string[] = [];
    if (node.oilFilmUm > 0.05) {
      parts.push(`${node.oilFilmUm.toFixed(1)}µm oil`);
    }
    if (node.frostThicknessUm > 0.5) {
      parts.push(`${node.frostThicknessUm.toFixed(0)}µm ice`);
    }
    if (node.ln2FilmG > 0.01) {
      parts.push(`LN₂`);
    }
    cFluid.textContent = parts.length > 0 ? parts.join(" · ") : "—";
    cFluid.style.color =
      node.frostThicknessUm > 1.0
        ? "#56d8ff"
        : node.oilFilmUm > 6.0
          ? "#f4d284"
          : "#98948a";
  }
}
