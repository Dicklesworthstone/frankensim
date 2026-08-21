// Seqlock state ring — the sim→render snapshot transport (plan §7.3,
// Round-4 P-08 header). Bead wf-root-guzez.1.7 (E0.7 prototype; E5.0 ships
// the production ABI on this exact protocol).
//
// Protocol: the writer marks a slot's sequence ODD, writes the payload,
// then publishes an EVEN sequence; readers snapshot the sequence, copy the
// payload, and retry if the sequence changed or was odd. All header traffic
// uses Atomics. The header additionally carries a RUN EPOCH and layout
// identity so a reader can never consume a complete-but-stale slot after a
// worker restart or ABI change (the Round-4 ABA hazard).

/** Int32 header slots (per ring, then per-slot sequence words). */
export const H_ABI_VERSION = 0;
export const H_LAYOUT_HASH = 1;
export const H_PAYLOAD_F64S = 2;
export const H_RUN_EPOCH = 3;
export const H_ANCHOR_PREFIX = 4;
export const H_PUBLISHED_SLOT = 5;
export const H_TICK_LO = 6;
export const H_TICK_HI = 7;
export const HEADER_I32S = 8;

export const ABI_VERSION = 1;

export interface SeqlockLayout {
  readonly slots: number;
  readonly payloadF64s: number;
}

/** Bytes for the shared buffer: header + per-slot sequence + payloads. */
export function seqlockBytes(layout: SeqlockLayout): number {
  const headerBytesAligned = Math.ceil((4 * (HEADER_I32S + layout.slots)) / 8) * 8;
  return headerBytesAligned + 8 * layout.slots * layout.payloadF64s;
}

function seqIndex(slot: number): number {
  return HEADER_I32S + slot;
}

function payloadOffsetF64(layout: SeqlockLayout, slot: number): number {
  // Payload region starts after the i32 header+sequence block, 8-aligned.
  const i32s = HEADER_I32S + layout.slots;
  const alignedBytes = Math.ceil((4 * i32s) / 8) * 8;
  return alignedBytes / 8 + slot * layout.payloadF64s;
}

// MEMORY-MODEL LAW (measured, not theoretical): the ES spec gives
// UNORDERED (plain) SAB accesses no ordering against the seq-cst
// Atomics — on Apple-Silicon ARM a writer's plain payload stores were
// OBSERVED becoming visible before its odd-mark store, and the
// transport battery consumed a torn payload (2026-08-20, wf-transport
// seqlock-stress, 2/3 repro in isolation). Every payload word
// therefore moves through Atomics on an i32 view (f64 split into two
// i32 halves via a thread-local scratch); the f64 view exists only for
// the scratch conversion. Do NOT "optimize" this back to a bulk copy.
export class SeqlockWriter {
  private readonly i32: Int32Array;
  private readonly layout: SeqlockLayout;
  private readonly scratchF64: Float64Array;
  private readonly scratchI32: Int32Array;
  private next = 0;

  constructor(
    sab: SharedArrayBuffer,
    layout: SeqlockLayout,
    runEpoch: number,
    layoutHash: number,
    anchorPrefix: number,
  ) {
    this.i32 = new Int32Array(sab);
    this.layout = layout;
    this.scratchF64 = new Float64Array(layout.payloadF64s);
    this.scratchI32 = new Int32Array(this.scratchF64.buffer);
    Atomics.store(this.i32, H_ABI_VERSION, ABI_VERSION);
    Atomics.store(this.i32, H_LAYOUT_HASH, layoutHash | 0);
    Atomics.store(this.i32, H_PAYLOAD_F64S, layout.payloadF64s);
    Atomics.store(this.i32, H_RUN_EPOCH, runEpoch | 0);
    Atomics.store(this.i32, H_ANCHOR_PREFIX, anchorPrefix | 0);
    Atomics.store(this.i32, H_PUBLISHED_SLOT, -1);
    for (let s = 0; s < layout.slots; s += 1) {
      Atomics.store(this.i32, seqIndex(s), 0);
    }
  }

  /** Publish one snapshot; never blocks. `fill` writes into a scratch
   * buffer; the transfer to shared memory is fully atomic (see the
   * memory-model law above). */
  publish(tick: number, fill: (payload: Float64Array) => void): void {
    const slot = this.next;
    this.next = (this.next + 1) % this.layout.slots;
    const si = seqIndex(slot);
    const seq = Atomics.load(this.i32, si);
    Atomics.store(this.i32, si, seq + 1); // odd: write in progress
    fill(this.scratchF64);
    const baseI32 = payloadOffsetF64(this.layout, slot) * 2;
    for (let k = 0; k < this.scratchI32.length; k += 1) {
      Atomics.store(this.i32, baseI32 + k, this.scratchI32[k]!);
    }
    Atomics.store(this.i32, H_TICK_LO, tick & 0x7fffffff);
    Atomics.store(this.i32, H_TICK_HI, Math.floor(tick / 0x80000000));
    Atomics.store(this.i32, si, seq + 2); // even: published
    Atomics.store(this.i32, H_PUBLISHED_SLOT, slot);
  }
}

export type SeqlockRejection =
  | "abi-version"
  | "layout-hash"
  | "payload-size"
  | "run-epoch"
  | "anchor-prefix"
  | "nothing-published"
  | "torn";

export class SeqlockReader {
  private readonly i32: Int32Array;
  private readonly layout: SeqlockLayout;
  private readonly scratchF64: Float64Array;
  private readonly scratchI32: Int32Array;
  private readonly expect: {
    runEpoch: number;
    layoutHash: number;
    anchorPrefix: number;
  };

  constructor(
    sab: SharedArrayBuffer,
    layout: SeqlockLayout,
    expect: { runEpoch: number; layoutHash: number; anchorPrefix: number },
  ) {
    this.i32 = new Int32Array(sab);
    this.layout = layout;
    this.scratchF64 = new Float64Array(layout.payloadF64s);
    this.scratchI32 = new Int32Array(this.scratchF64.buffer);
    this.expect = expect;
  }

  /**
   * Try to copy the latest published snapshot into `out`. Returns the tick
   * on success or a typed rejection. Header identity is validated BEFORE
   * any payload access (Round-4 P-08: a seqlock alone cannot stop a reader
   * consuming a complete stale slot).
   */
  read(out: Float64Array, retries = 4): number | SeqlockRejection {
    if (Atomics.load(this.i32, H_ABI_VERSION) !== ABI_VERSION) {
      return "abi-version";
    }
    if ((Atomics.load(this.i32, H_LAYOUT_HASH) | 0) !== (this.expect.layoutHash | 0)) {
      return "layout-hash";
    }
    if (Atomics.load(this.i32, H_PAYLOAD_F64S) !== this.layout.payloadF64s) {
      return "payload-size";
    }
    if ((Atomics.load(this.i32, H_RUN_EPOCH) | 0) !== (this.expect.runEpoch | 0)) {
      return "run-epoch";
    }
    if ((Atomics.load(this.i32, H_ANCHOR_PREFIX) | 0) !== (this.expect.anchorPrefix | 0)) {
      return "anchor-prefix";
    }
    for (let attempt = 0; attempt < retries; attempt += 1) {
      const slot = Atomics.load(this.i32, H_PUBLISHED_SLOT);
      if (slot < 0) {
        return "nothing-published";
      }
      const si = seqIndex(slot);
      const before = Atomics.load(this.i32, si);
      if ((before & 1) === 1) {
        continue; // write in progress
      }
      // Atomic payload copy (memory-model law above) into the scratch,
      // then out — `out` is only touched on a validated read.
      const baseI32 = payloadOffsetF64(this.layout, slot) * 2;
      for (let k = 0; k < this.scratchI32.length; k += 1) {
        this.scratchI32[k] = Atomics.load(this.i32, baseI32 + k);
      }
      const tickLo = Atomics.load(this.i32, H_TICK_LO);
      const tickHi = Atomics.load(this.i32, H_TICK_HI);
      const after = Atomics.load(this.i32, si);
      if (after === before) {
        out.set(this.scratchF64);
        return tickHi * 0x80000000 + tickLo;
      }
    }
    return "torn";
  }
}
