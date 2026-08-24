// Microbench harness (bead wf-root-guzez.1.6.1, E0.6a). Grounds the §7.2
// budget rows: warmup discard, batched timing, nearest-rank percentiles,
// JSONL rows. Numbers from this harness are MEASUREMENTS OF A HOST — the
// plan's acceptance contract (§7.2.1) is satisfied only by qualified-device
// rows, never by any single machine.

export interface BenchResult {
  readonly name: string;
  readonly batchSize: number;
  readonly samples: number;
  readonly p50_us: number;
  readonly p95_us: number;
  readonly p99_us: number;
  /** Throughput over the timed window, or null when the host timer
   * quantum resolved EVERY sample to zero elapsed time (coarsened
   * non-isolated origins): totals are unresolvable, not zero-cost. */
  readonly opsPerSec: number | null;
}

export function percentileOf(sorted: readonly number[], q: number): number {
  if (sorted.length === 0) {
    return 0;
  }
  return sorted[Math.min(sorted.length - 1, Math.floor(q * sorted.length))]!;
}

/**
 * Time `fn` (which runs `batchSize` operations internally) `samples` times
 * after `warmup` discarded runs. Per-OPERATION microseconds reported.
 */
export function bench(
  name: string,
  batchSize: number,
  fn: () => void,
  samples = 60,
  warmup = 10,
): BenchResult {
  for (let i = 0; i < warmup; i += 1) {
    fn();
  }
  const perOpUs: number[] = [];
  let totalMs = 0;
  for (let i = 0; i < samples; i += 1) {
    const t0 = performance.now();
    fn();
    const dt = performance.now() - t0;
    totalMs += dt;
    perOpUs.push((dt * 1000) / batchSize);
  }
  perOpUs.sort((a, b) => a - b);
  return {
    name,
    batchSize,
    samples,
    p50_us: Number(percentileOf(perOpUs, 0.5).toFixed(4)),
    p95_us: Number(percentileOf(perOpUs, 0.95).toFixed(4)),
    p99_us: Number(percentileOf(perOpUs, 0.99).toFixed(4)),
    opsPerSec: totalMs > 0 ? Math.round((samples * batchSize * 1000) / totalMs) : null,
  };
}
