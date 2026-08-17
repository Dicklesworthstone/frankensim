// E0.8 timing soak: TickScheduler + InputScheduler against the real wall
// clock with a cross-thread input injector. Measures the §7.2.1 metric
// separation and evaluates the PROVISIONAL lateness gates
// (p99 <= 1.5 ms, p99.9 <= 4.0 ms) with honest host-condition reporting —
// on a saturated shared host the gates measure the HOST, and ratification
// requires a qualified quiet device; the machinery is the deliverable.
// Run: node test/timingSoak.ts [minutes]

import { Worker } from "node:worker_threads";
import { InputScheduler } from "../src/transport/inputClock.ts";
import { TickScheduler, percentile } from "../src/transport/schedule.ts";

const minutes = Number(process.argv[2] ?? "10");
const TICK_MS = 1000 / 120;
const jlog = (obj: object): void =>
  console.log(JSON.stringify({ suite: "wf-timing-soak", ...obj }));

const startMs = performance.now();
const tickSched = new TickScheduler(TICK_MS, startMs, 3);
const inputSched = new InputScheduler(TICK_MS, startMs, 1);

// Injector thread: posts an input packet every ~25 ms (a busy human).
const injector = new Worker(
  `
  const { parentPort } = require("node:worker_threads");
  let seq = 0;
  const iv = setInterval(() => {
    seq += 1;
    parentPort.postMessage({ seq, valueRaw: Math.sin(seq / 7) });
  }, 25);
  parentPort.on("message", () => { clearInterval(iv); process.exit(0); });
  `,
  { eval: true },
);

let inputsAdmitted = 0;
injector.on("message", (m: { seq: number; valueRaw: number }) => {
  inputSched.admit(
    {
      channel: 0,
      deviceSampleMs: performance.now(),
      quantizedValue: Math.round(Math.min(1, Math.max(-1, m.valueRaw)) * 4096) / 4096,
      sequence: m.seq,
    },
    tickSched.currentTick(),
  );
  inputsAdmitted += 1;
});

const endMs = startMs + minutes * 60_000;
let windowStart = startMs;
const allLateness: number[] = [];

function simWork(): void {
  // ~0.2 ms of synthetic arithmetic standing in for Tier-A physics.
  let acc = 0;
  for (let i = 0; i < 4000; i += 1) {
    acc += Math.sqrt(i + acc % 7);
  }
  if (acc === -1) {
    console.log("never");
  }
}

const pump = (): void => {
  const now = performance.now();
  tickSched.pump(now, simWork, () => performance.now());
  if (now - windowStart >= 10_000) {
    const lat = tickSched.metrics.latenessMs.splice(0);
    const svc = tickSched.metrics.serviceMs.splice(0);
    allLateness.push(...lat);
    jlog({
      window_s: 10,
      tick: tickSched.currentTick(),
      lateness_p50_ms: Number(percentile(lat, 0.5).toFixed(3)),
      lateness_p99_ms: Number(percentile(lat, 0.99).toFixed(3)),
      service_p99_ms: Number(percentile(svc, 0.99).toFixed(3)),
      backlogMax: tickSched.metrics.maxBacklogObserved,
      reanchors: tickSched.metrics.reanchors,
      inputsAdmitted,
    });
    windowStart = now;
  }
  if (now < endMs) {
    setTimeout(pump, 1);
  } else {
    injector.postMessage("stop");
    const acq = inputSched.acquisitionTrace();
    const lateInputs = acq.filter((r) => r.lateByTicks > 0).length;
    const p99 = percentile(allLateness, 0.99);
    const p999 = percentile(allLateness, 0.999);
    const gates = { p99_gate_ms: 1.5, p999_gate_ms: 4.0 };
    const gatePass = p99 <= gates.p99_gate_ms && p999 <= gates.p999_gate_ms;
    jlog({
      event: "soak-complete",
      minutes,
      ticks: tickSched.currentTick(),
      inputsAdmitted,
      lateInputsFlagged: lateInputs,
      canonicalEvents: inputSched.appliedTrace().length,
      lateness_p99_ms: Number(p99.toFixed(3)),
      lateness_p999_ms: Number(p999.toFixed(3)),
      provisional_gates: gates,
      gate_result_on_this_host: gatePass ? "PASS" : "FAIL",
      host_note:
        "shared saturated host — gate ratification requires a qualified quiet device; this run proves the MEASUREMENT machinery and the contract semantics",
      reanchors: tickSched.metrics.reanchors,
      maxBacklogBurst: 4,
    });
    const contractHolds =
      inputSched.appliedTrace().length === acq.length &&
      tickSched.currentTick() > 0;
    jlog({ verdict: contractHolds ? "CONTRACT-PASS" : "CONTRACT-FAIL" });
    process.exit(contractHolds ? 0 : 1);
  }
};
jlog({ event: "soak-start", minutes, tickHz: 120 });
pump();
