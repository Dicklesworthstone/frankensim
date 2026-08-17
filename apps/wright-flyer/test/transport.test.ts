// Transport protocol battery (bead wf-root-guzez.1.7, E0.7). Runs under
// plain `node --test` with worker_threads + SharedArrayBuffer — the same
// Atomics semantics the browser provides. Every assertion below is a
// DONE-WHEN receipt: torn-read impossibility under stress, header-identity
// rejection before payload access (the restart/slot-reuse ABA twin), the
// never-blocks drop counter, lease immutability, and pool starvation/ack
// accounting. Structured JSONL result lines go to stdout for the E2E logs.

import assert from "node:assert/strict";
import { test } from "node:test";
import { Worker } from "node:worker_threads";
import { fileURLToPath } from "node:url";
import {
  SeqlockReader,
  SeqlockWriter,
  seqlockBytes,
  type SeqlockLayout,
} from "../src/transport/seqlock.ts";
import {
  LeasedRingReader,
  LeasedRingWriter,
  SLOT_LEASED,
  leasedRingBytes,
  type LeasedRingLayout,
} from "../src/transport/leasedRing.ts";
import { TransferablePool } from "../src/transport/pool.ts";

const jlog = (obj: object): void =>
  console.log(JSON.stringify({ suite: "wf-transport", ...obj }));

// A payload pattern where every word equals tick*3+k lets a reader prove a
// copy is internally consistent — any interleaving of two writes breaks it.
function fillPattern(payload: Float64Array, tick: number): void {
  for (let k = 0; k < payload.length; k += 1) {
    payload[k] = tick * 3 + k;
  }
}

function patternConsistent(payload: Float64Array): boolean {
  const tick3 = payload[0]!;
  for (let k = 1; k < payload.length; k += 1) {
    if (payload[k] !== tick3 + k) {
      return false;
    }
  }
  return true;
}

test("seqlock: no torn reads under a hammering cross-thread writer", async () => {
  const layout: SeqlockLayout = { slots: 3, payloadF64s: 256 };
  const sab = new SharedArrayBuffer(seqlockBytes(layout));
  // Writer thread hammers publications as fast as it can for ~2 s.
  const writerSrc = `
    const { workerData, parentPort } = require("node:worker_threads");
    const path = ${JSON.stringify(fileURLToPath(new URL("../src/transport/seqlock.ts", import.meta.url)))};
    import(path).then(({ SeqlockWriter }) => {
      const writer = new SeqlockWriter(workerData.sab, workerData.layout, 7, 42, 9);
      let tick = 0;
      const end = Date.now() + 2000;
      while (Date.now() < end) {
        tick += 1;
        writer.publish(tick, (p) => {
          for (let k = 0; k < p.length; k += 1) p[k] = tick * 3 + k;
        });
      }
      parentPort.postMessage(tick);
    });
  `;
  const worker = new Worker(writerSrc, {
    eval: true,
    workerData: { sab, layout },
  });
  const writerTicks: Promise<number> = new Promise((resolve, reject) => {
    worker.once("message", resolve);
    worker.once("error", reject);
  });
  // Reader hammers concurrently on the main thread.
  const reader = new SeqlockReader(sab, layout, {
    runEpoch: 7,
    layoutHash: 42,
    anchorPrefix: 9,
  });
  const out = new Float64Array(layout.payloadF64s);
  let reads = 0;
  let successes = 0;
  let torn = 0;
  let inconsistent = 0;
  const end = Date.now() + 2000;
  while (Date.now() < end) {
    const result = reader.read(out);
    reads += 1;
    if (typeof result === "number") {
      successes += 1;
      if (!patternConsistent(out)) {
        inconsistent += 1;
      }
    } else if (result === "torn") {
      torn += 1; // retries exhausted is legal; CONSUMING torn data is not
    }
  }
  const ticks = await writerTicks;
  await worker.terminate();
  jlog({ case: "seqlock-stress", reads, successes, torn, inconsistent, writerTicks: ticks });
  assert.ok(successes > 1000, `expected heavy successful traffic, got ${successes}`);
  assert.equal(inconsistent, 0, "a torn payload was CONSUMED — protocol violation");
});

test("seqlock: header identity rejects stale/foreign rings before payload (ABA twin)", () => {
  const layout: SeqlockLayout = { slots: 3, payloadF64s: 8 };
  const sab = new SharedArrayBuffer(seqlockBytes(layout));
  const writer = new SeqlockWriter(sab, layout, /*epoch*/ 1, /*hash*/ 42, /*anchor*/ 9);
  writer.publish(5, (p) => fillPattern(p, 5));
  const out = new Float64Array(layout.payloadF64s);
  // Correct identity reads fine.
  const good = new SeqlockReader(sab, layout, { runEpoch: 1, layoutHash: 42, anchorPrefix: 9 });
  assert.equal(good.read(out), 5);
  // Restart twin: a reader expecting the NEXT epoch must reject the intact,
  // fully published, internally consistent old ring.
  const afterRestart = new SeqlockReader(sab, layout, { runEpoch: 2, layoutHash: 42, anchorPrefix: 9 });
  assert.equal(afterRestart.read(out), "run-epoch");
  // Layout / anchor / size mismatches likewise reject before payload.
  const wrongHash = new SeqlockReader(sab, layout, { runEpoch: 1, layoutHash: 43, anchorPrefix: 9 });
  assert.equal(wrongHash.read(out), "layout-hash");
  const wrongAnchor = new SeqlockReader(sab, layout, { runEpoch: 1, layoutHash: 42, anchorPrefix: 8 });
  assert.equal(wrongAnchor.read(out), "anchor-prefix");
  const wrongSize = new SeqlockReader(sab, { slots: 3, payloadF64s: 16 }, { runEpoch: 1, layoutHash: 42, anchorPrefix: 9 });
  assert.equal(wrongSize.read(new Float64Array(16)), "payload-size");
  jlog({ case: "seqlock-aba-twin", verdict: "all stale identities rejected" });
});

test("leased ring: writer never blocks, drops count, leases stay immutable", () => {
  const layout: LeasedRingLayout = { slots: 3, payloadF64s: 32 };
  const sab = new SharedArrayBuffer(leasedRingBytes(layout));
  const writer = new LeasedRingWriter(sab, layout, 1);
  const reader = new LeasedRingReader(sab, layout, 1);

  assert.ok(writer.publish(1, (p) => fillPattern(p, 1)));
  const leaseA = reader.lease();
  assert.ok(leaseA && leaseA.tick === 1);
  const frozen = Array.from(leaseA!.payload);

  // Fill the remaining capacity while A stays leased…
  assert.ok(writer.publish(2, (p) => fillPattern(p, 2)));
  assert.ok(writer.publish(3, (p) => fillPattern(p, 3)));
  // …then keep publishing: slot A is LEASED, slot(latest=3) is protected,
  // so the writer reclaims the superseded published slot (tick 2) and, once
  // options run out, DROPS rather than blocking or touching the lease.
  assert.ok(writer.publish(4, (p) => fillPattern(p, 4)), "reclaims superseded slot");
  const leaseB = reader.lease(); // leases latest (tick 4)
  assert.ok(leaseB && leaseB.tick === 4);
  const dropped = writer.publish(5, (p) => fillPattern(p, 5));
  const dropped2 = writer.publish(6, (p) => fillPattern(p, 6));
  // With two leases held and one latest-published slot, at least one of
  // these publications must drop (never block, never corrupt).
  assert.ok(!dropped || !dropped2, "with 2/3 slots leased the writer must drop");
  assert.ok(writer.dropCount() >= 1, "drop counter must increment");
  // Lease A's payload never changed underneath the reader.
  assert.deepEqual(Array.from(leaseA!.payload), frozen, "leased payload mutated!");
  assert.ok(patternConsistent(leaseA!.payload));
  reader.release(leaseA!);
  reader.release(leaseB!);
  jlog({ case: "leased-ring", drops: writer.dropCount(), verdict: "never-blocks + immutability hold" });
});

test("leased ring: epoch mismatch refuses leases after restart", () => {
  const layout: LeasedRingLayout = { slots: 3, payloadF64s: 8 };
  const sab = new SharedArrayBuffer(leasedRingBytes(layout));
  const writer = new LeasedRingWriter(sab, layout, 1);
  writer.publish(1, (p) => fillPattern(p, 1));
  const staleReader = new LeasedRingReader(sab, layout, 2);
  assert.equal(staleReader.lease(), null, "stale-epoch lease must refuse");
  jlog({ case: "leased-ring-epoch-twin", verdict: "stale epoch refused" });
});

test("leased ring: fewer than three slots is a construction refusal", () => {
  const layout: LeasedRingLayout = { slots: 2, payloadF64s: 8 };
  const sab = new SharedArrayBuffer(leasedRingBytes(layout));
  assert.throws(() => new LeasedRingWriter(sab, layout, 1), /3 slots/);
});

test("transferable pool: starvation drops, ack restores, metrics account bytes", () => {
  const pool = new TransferablePool(2, 8 * 64);
  const source = new Float64Array(64).fill(3.25);
  const a = pool.pack(source);
  const b = pool.pack(source);
  assert.ok(a && b);
  assert.equal(pool.pack(source), null, "exhausted pool must drop, not allocate");
  let m = pool.snapshotMetrics();
  assert.equal(m.publications, 2);
  assert.equal(m.drops, 1);
  assert.equal(m.starvationEvents, 1);
  assert.equal(m.outstanding, 2);
  assert.equal(m.copiedBytes, 2 * 8 * 64);
  assert.ok(new Float64Array(a!)[0] === 3.25, "pack copies content");
  pool.acknowledge(a!);
  assert.ok(pool.pack(source), "acknowledged buffer is reusable");
  m = pool.snapshotMetrics();
  assert.equal(m.publications, 3);
  jlog({ case: "pool", metrics: m });
});
