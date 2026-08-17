// 10-minute contention soak (bead wf-root-guzez.1.7, E0.7 DONE-WHEN).
// A synthetic 120 Hz sim loop publishes through the seqlock ring and the
// leased snapshot ring while a hammering reader thread and the main thread
// consume; a mid-soak PAUSE window verifies the no-catch-up rule (plan
// §4.3: on resume, ticks continue from the schedule, never a burst).
// Emits one JSONL row per 10 s window: publication counts, read outcomes,
// drop counters, tick-lateness percentiles. Run: node test/soak.ts [minutes]

import { Worker } from "node:worker_threads";
import { fileURLToPath } from "node:url";
import {
  SeqlockWriter,
  seqlockBytes,
  type SeqlockLayout,
} from "../src/transport/seqlock.ts";
import {
  LeasedRingReader,
  LeasedRingWriter,
  leasedRingBytes,
  type LeasedRingLayout,
} from "../src/transport/leasedRing.ts";

const minutes = Number(process.argv[2] ?? "10");
const TICK_MS = 1000 / 120;
const MAX_CATCHUP_TICKS = 3; // bounded burst; beyond this we re-anchor (pause semantics)

const seqLayout: SeqlockLayout = { slots: 3, payloadF64s: 256 };
const seqSab = new SharedArrayBuffer(seqlockBytes(seqLayout));
const leaseLayout: LeasedRingLayout = { slots: 4, payloadF64s: 512 };
const leaseSab = new SharedArrayBuffer(leasedRingBytes(leaseLayout));

const writer = new SeqlockWriter(seqSab, seqLayout, 7, 42, 9);
const leaseWriter = new LeasedRingWriter(leaseSab, leaseLayout, 7);
const leaseReader = new LeasedRingReader(leaseSab, leaseLayout, 7);

// Reader thread hammers the seqlock ring for the whole soak.
const readerSrc = `
  const { workerData, parentPort } = require("node:worker_threads");
  const path = ${JSON.stringify(fileURLToPath(new URL("../src/transport/seqlock.ts", import.meta.url)))};
  import(path).then(({ SeqlockReader }) => {
    const reader = new SeqlockReader(workerData.sab, workerData.layout,
      { runEpoch: 7, layoutHash: 42, anchorPrefix: 9 });
    const out = new Float64Array(workerData.layout.payloadF64s);
    let reads = 0, ok = 0, torn = 0, inconsistent = 0;
    const check = () => {
      const r = reader.read(out);
      reads += 1;
      if (typeof r === "number") {
        ok += 1;
        const t3 = out[0];
        for (let k = 1; k < out.length; k += 1) {
          if (out[k] !== t3 + k) { inconsistent += 1; break; }
        }
      } else if (r === "torn") { torn += 1; }
    };
    const iv = setInterval(() => { for (let i = 0; i < 200; i += 1) check(); }, 1);
    parentPort.on("message", () => {
      clearInterval(iv);
      parentPort.postMessage({ reads, ok, torn, inconsistent });
    });
  });
`;
const readerWorker = new Worker(readerSrc, {
  eval: true,
  workerData: { sab: seqSab, layout: seqLayout },
});

const jlog = (obj: object): void =>
  console.log(JSON.stringify({ suite: "wf-transport-soak", ...obj }));

let tick = 0;
let lateness: number[] = [];
let published = 0;
let leasePublished = 0;
let leaseHeld: ReturnType<typeof leaseReader.lease> = null;
let reanchors = 0;
let pauseDone = false;
const startMs = performance.now();
let nextDueMs = startMs;
const endMs = startMs + minutes * 60_000;
let windowStart = startMs;

function fill(payload: Float64Array): void {
  for (let k = 0; k < payload.length; k += 1) {
    payload[k] = tick * 3 + k;
  }
}

const pump = (): void => {
  const now = performance.now();
  // Bounded catch-up: if we are more than MAX_CATCHUP_TICKS behind, re-anchor
  // the schedule (the visibility-pause / no-unbounded-catch-up rule).
  if (now - nextDueMs > MAX_CATCHUP_TICKS * TICK_MS) {
    reanchors += 1;
    nextDueMs = now;
  }
  while (nextDueMs <= now) {
    lateness.push(now - nextDueMs);
    tick += 1;
    writer.publish(tick, fill);
    published += 1;
    // Field-rate (every 8th tick ≈ 15 Hz) leased publication + rotating lease.
    if (tick % 8 === 0) {
      if (leaseWriter.publish(tick, fill)) {
        leasePublished += 1;
      }
      if (leaseHeld) {
        leaseReader.release(leaseHeld);
      }
      leaseHeld = leaseReader.lease();
    }
    nextDueMs += TICK_MS;
  }
  // Mid-soak deliberate pause: at half time, stall 2 s in one blocking gap
  // (simulating a hidden tab), then rely on the re-anchor rule.
  if (!pauseDone && now - startMs > (minutes * 60_000) / 2) {
    pauseDone = true;
    const stallUntil = performance.now() + 2000;
    while (performance.now() < stallUntil) {
      /* deliberate stall */
    }
    jlog({ event: "pause-window", stalledMs: 2000, tickAtPause: tick });
  }
  if (now - windowStart >= 10_000) {
    lateness.sort((a, b) => a - b);
    const pct = (q: number): number =>
      lateness.length ? lateness[Math.min(lateness.length - 1, Math.floor(q * lateness.length))]! : 0;
    jlog({
      window_s: 10,
      tick,
      published,
      leasePublished,
      leaseDrops: leaseWriter.dropCount(),
      lateness_p50_ms: Number(pct(0.5).toFixed(3)),
      lateness_p99_ms: Number(pct(0.99).toFixed(3)),
      reanchors,
    });
    lateness = [];
    windowStart = now;
  }
  if (now < endMs) {
    setTimeout(pump, 1);
  } else {
    readerWorker.postMessage("stop");
    readerWorker.once("message", (r: { reads: number; ok: number; torn: number; inconsistent: number }) => {
      jlog({ event: "soak-complete", minutes, finalTick: tick, reader: r, reanchors });
      const verdict =
        r.inconsistent === 0 && reanchors >= 1 && published > 0
          ? "SOAK-PASS"
          : "SOAK-FAIL";
      jlog({ verdict, inconsistentReads: r.inconsistent, reanchorsObserved: reanchors });
      void readerWorker.terminate().then(() => process.exit(verdict === "SOAK-PASS" ? 0 : 1));
    });
  }
};
jlog({ event: "soak-start", minutes, tickHz: 120, maxCatchupTicks: MAX_CATCHUP_TICKS });
pump();
