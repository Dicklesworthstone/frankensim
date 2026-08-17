// Tick scheduler with the §7.2.1 metric separation (Round-4 S-15): compute
// SERVICE time, scheduling LATENESS, and BACKLOG are distinct quantities —
// a 3 ms tick scheduled 6 ms late still misses. Bounded catch-up with
// re-anchoring implements the visibility-pause/no-unbounded-catch-up rule.
// Bead wf-root-guzez.1.8 (E0.8).

export interface TickMetrics {
  ticksRun: number;
  reanchors: number;
  maxBacklogObserved: number;
  latenessMs: number[]; // start lateness per tick (caller may drain)
  serviceMs: number[]; // compute service time per tick (caller may drain)
}

export class TickScheduler {
  private readonly tickMs: number;
  private readonly maxCatchupTicks: number;
  private nextDueMs: number;
  private tick = 0;
  readonly metrics: TickMetrics = {
    ticksRun: 0,
    reanchors: 0,
    maxBacklogObserved: 0,
    latenessMs: [],
    serviceMs: [],
  };

  constructor(tickMs: number, startMs: number, maxCatchupTicks = 3) {
    this.tickMs = tickMs;
    this.nextDueMs = startMs;
    this.maxCatchupTicks = maxCatchupTicks;
  }

  currentTick(): number {
    return this.tick;
  }

  /**
   * Run due ticks at monotonic time `nowMs`. The backlog (due-but-unrun
   * ticks) is measured BEFORE running; if it exceeds the bounded-catch-up
   * limit the schedule re-anchors (pause semantics: simulation time simply
   * does not advance for the gap — no burst, no wall-clock catch-up).
   */
  pump(nowMs: number, runTick: (tick: number) => void, nowFn: () => number): number {
    let backlog = Math.floor((nowMs - this.nextDueMs) / this.tickMs) +
      (nowMs >= this.nextDueMs ? 1 : 0);
    if (backlog < 0) {
      backlog = 0;
    }
    if (backlog > this.metrics.maxBacklogObserved) {
      this.metrics.maxBacklogObserved = backlog;
    }
    if (backlog > this.maxCatchupTicks) {
      this.metrics.reanchors += 1;
      this.nextDueMs = nowMs; // re-anchor: run exactly one tick now
    }
    let ran = 0;
    while (this.nextDueMs <= nowMs) {
      this.metrics.latenessMs.push(nowMs - this.nextDueMs);
      this.tick += 1;
      const t0 = nowFn();
      runTick(this.tick);
      this.metrics.serviceMs.push(nowFn() - t0);
      this.metrics.ticksRun += 1;
      this.nextDueMs += this.tickMs;
      ran += 1;
      if (ran > this.maxCatchupTicks) {
        break; // hard safety: never a burst beyond the bound
      }
    }
    return ran;
  }
}

/** Percentile over a sample array (nearest-rank; empty → 0). */
export function percentile(samples: readonly number[], q: number): number {
  if (samples.length === 0) {
    return 0;
  }
  const sorted = [...samples].sort((a, b) => a - b);
  const index = Math.min(sorted.length - 1, Math.floor(q * sorted.length));
  return sorted[index]!;
}
