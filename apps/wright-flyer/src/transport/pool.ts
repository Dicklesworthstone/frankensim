// Transferable buffer pool — the non-SAB fallback transport (plan §7.3,
// Round-2/Round-3: a buffer backed by wasm linear memory cannot be
// transferred away while the sim keeps using that memory, so every
// publication performs ONE explicit pack/copy into a pooled standalone
// ArrayBuffer, ownership transfers to the consumer, and the buffer returns
// via an acknowledgement channel). Bead wf-root-guzez.1.7 (E0.7).
//
// Copy bytes/time, pool starvation, and acknowledgement latency are
// MEASURED — they are §7.2 acceptance-budget lines, not noise.

export interface PoolMetrics {
  publications: number;
  drops: number;
  copiedBytes: number;
  copyTimeMs: number;
  starvationEvents: number;
  outstanding: number;
}

export class TransferablePool {
  private readonly free: ArrayBuffer[] = [];
  private readonly bytesPerBuffer: number;
  private readonly capacity: number;
  private readonly metrics: PoolMetrics = {
    publications: 0,
    drops: 0,
    copiedBytes: 0,
    copyTimeMs: 0,
    starvationEvents: 0,
    outstanding: 0,
  };

  constructor(capacity: number, bytesPerBuffer: number) {
    this.capacity = capacity;
    this.bytesPerBuffer = bytesPerBuffer;
    for (let i = 0; i < capacity; i += 1) {
      this.free.push(new ArrayBuffer(bytesPerBuffer));
    }
  }

  /**
   * Pack `source` into a pooled buffer for transfer. Returns the buffer
   * (ownership moves to the caller/consumer) or null on starvation — the
   * caller drops the publication rather than blocking or allocating.
   */
  pack(source: Float64Array): ArrayBuffer | null {
    const buffer = this.free.pop();
    if (buffer === undefined) {
      this.metrics.drops += 1;
      this.metrics.starvationEvents += 1;
      return null;
    }
    const t0 = performance.now();
    new Float64Array(buffer).set(source);
    this.metrics.copyTimeMs += performance.now() - t0;
    this.metrics.copiedBytes += source.byteLength;
    this.metrics.publications += 1;
    this.metrics.outstanding += 1;
    return buffer;
  }

  /** Acknowledgement path: the consumer returns the (possibly detached-and-
   * replaced) buffer to the pool. A detached buffer is replaced by a fresh
   * allocation of the same size so capacity is conserved. */
  acknowledge(buffer: ArrayBuffer): void {
    this.metrics.outstanding -= 1;
    if (this.free.length < this.capacity) {
      this.free.push(
        buffer.byteLength === this.bytesPerBuffer
          ? buffer
          : new ArrayBuffer(this.bytesPerBuffer),
      );
    }
  }

  snapshotMetrics(): PoolMetrics {
    return { ...this.metrics };
  }
}
