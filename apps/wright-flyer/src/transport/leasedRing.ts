// Leased immutable snapshot ring — the sim→field FieldSourceSnapshot
// transport (plan §7.3, Round-2/Round-4). Bead wf-root-guzez.1.7 (E0.7).
//
// States per slot: FREE → WRITING → PUBLISHED → LEASED → FREE, all
// transitions via Atomics.compareExchange. The WRITER NEVER BLOCKS: with no
// FREE slot it skips the publication and increments a drop counter (the
// field-age UI surfaces staleness). A reader leases the latest PUBLISHED
// slot; payload ownership is immutable while LEASED. Minimum three slots —
// a double buffer is unsafe while a reader holds a lease.

export const SLOT_FREE = 0;
export const SLOT_WRITING = 1;
export const SLOT_PUBLISHED = 2;
export const SLOT_LEASED = 3;

/** Control words per ring: [dropCount, latestPublishedSlot, runEpoch]. */
const C_DROPS = 0;
const C_LATEST = 1;
const C_EPOCH = 2;
const CONTROL_I32S = 3;

export interface LeasedRingLayout {
  readonly slots: number; // must be >= 3
  readonly payloadF64s: number;
}

export function leasedRingBytes(layout: LeasedRingLayout): number {
  const controlBytesAligned =
    Math.ceil((4 * (CONTROL_I32S + 2 * layout.slots)) / 8) * 8;
  return controlBytesAligned + 8 * layout.slots * layout.payloadF64s;
}

function stateIndex(slot: number): number {
  return CONTROL_I32S + slot;
}

function tickIndex(layout: LeasedRingLayout, slot: number): number {
  return CONTROL_I32S + layout.slots + slot;
}

function payloadOffsetF64(layout: LeasedRingLayout, slot: number): number {
  const i32s = CONTROL_I32S + 2 * layout.slots;
  const alignedBytes = Math.ceil((4 * i32s) / 8) * 8;
  return alignedBytes / 8 + slot * layout.payloadF64s;
}

export class LeasedRingWriter {
  private readonly i32: Int32Array;
  private readonly f64: Float64Array;
  private readonly layout: LeasedRingLayout;

  constructor(sab: SharedArrayBuffer, layout: LeasedRingLayout, runEpoch: number) {
    if (layout.slots < 3) {
      throw new Error("leased ring requires >= 3 slots (double buffer is unsafe)");
    }
    this.i32 = new Int32Array(sab);
    this.f64 = new Float64Array(sab);
    this.layout = layout;
    Atomics.store(this.i32, C_DROPS, 0);
    Atomics.store(this.i32, C_LATEST, -1);
    Atomics.store(this.i32, C_EPOCH, runEpoch | 0);
    for (let s = 0; s < layout.slots; s += 1) {
      Atomics.store(this.i32, stateIndex(s), SLOT_FREE);
    }
  }

  /**
   * Publish a snapshot into a FREE slot, or drop (never blocks, never
   * touches a LEASED slot). A superseded PUBLISHED slot is reclaimed —
   * only the latest publication matters to the field service.
   */
  publish(tick: number, fill: (payload: Float64Array) => void): boolean {
    let slot = -1;
    for (let s = 0; s < this.layout.slots; s += 1) {
      if (
        Atomics.compareExchange(this.i32, stateIndex(s), SLOT_FREE, SLOT_WRITING) ===
        SLOT_FREE
      ) {
        slot = s;
        break;
      }
    }
    if (slot < 0) {
      // Reclaim a stale PUBLISHED slot (not the latest) if any; else drop.
      const latest = Atomics.load(this.i32, C_LATEST);
      for (let s = 0; s < this.layout.slots; s += 1) {
        if (s === latest) {
          continue;
        }
        if (
          Atomics.compareExchange(
            this.i32,
            stateIndex(s),
            SLOT_PUBLISHED,
            SLOT_WRITING,
          ) === SLOT_PUBLISHED
        ) {
          slot = s;
          break;
        }
      }
    }
    if (slot < 0) {
      Atomics.add(this.i32, C_DROPS, 1);
      return false;
    }
    const base = payloadOffsetF64(this.layout, slot);
    fill(this.f64.subarray(base, base + this.layout.payloadF64s));
    Atomics.store(this.i32, tickIndex(this.layout, slot), tick | 0);
    Atomics.store(this.i32, stateIndex(slot), SLOT_PUBLISHED);
    Atomics.store(this.i32, C_LATEST, slot);
    return true;
  }

  dropCount(): number {
    return Atomics.load(this.i32, C_DROPS);
  }
}

export interface Lease {
  readonly slot: number;
  readonly tick: number;
  readonly payload: Float64Array;
}

export class LeasedRingReader {
  private readonly i32: Int32Array;
  private readonly f64: Float64Array;
  private readonly layout: LeasedRingLayout;
  private readonly runEpoch: number;

  constructor(sab: SharedArrayBuffer, layout: LeasedRingLayout, runEpoch: number) {
    this.i32 = new Int32Array(sab);
    this.f64 = new Float64Array(sab);
    this.layout = layout;
    this.runEpoch = runEpoch | 0;
  }

  /** Lease the latest PUBLISHED slot (null if none / epoch mismatch). */
  lease(): Lease | null {
    if (Atomics.load(this.i32, C_EPOCH) !== this.runEpoch) {
      return null; // restart happened; a stale ring is never consumed
    }
    for (let attempt = 0; attempt < this.layout.slots; attempt += 1) {
      const slot = Atomics.load(this.i32, C_LATEST);
      if (slot < 0) {
        return null;
      }
      if (
        Atomics.compareExchange(
          this.i32,
          stateIndex(slot),
          SLOT_PUBLISHED,
          SLOT_LEASED,
        ) === SLOT_PUBLISHED
      ) {
        const base = payloadOffsetF64(this.layout, slot);
        return {
          slot,
          tick: Atomics.load(this.i32, tickIndex(this.layout, slot)),
          payload: this.f64.subarray(base, base + this.layout.payloadF64s),
        };
      }
      // The latest moved or is being rewritten; retry with the new latest.
    }
    return null;
  }

  release(lease: Lease): void {
    Atomics.store(this.i32, stateIndex(lease.slot), SLOT_FREE);
  }
}
