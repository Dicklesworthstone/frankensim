// The page's trend recorder and event log, fed by the REAL compiled engine
// (`src/wasm-pkg/fs_steam_turbine_bg.wasm`) through the page's own snapshot decoder.
// Node strips the TypeScript types of src/trends.ts and src/snapshotLayout.ts on import.

import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { dirname, resolve } from "node:path";
import { test } from "node:test";
import { fileURLToPath } from "node:url";

import initWasm, { SteamTurbineSim } from "../src/wasm-pkg/fs_steam_turbine.js";
import { decodeSnapshot } from "../src/snapshotLayout.ts";
import { COARSE_DT_S, EventLog, TREND_CHANNELS, TrendRecorder, eventsCsv, trendCsv } from "../src/trends.ts";

const here = dirname(fileURLToPath(import.meta.url));
await initWasm({ module_or_path: readFileSync(resolve(here, "../src/wasm-pkg/fs_steam_turbine_bg.wasm")) });

const DT = 0.02;
const run = (sim, rec, log, steps) => {
  for (let i = 0; i < steps; i++) {
    assert.equal(sim.step(DT), true, sim.last_error());
    const snap = decodeSnapshot(sim.snapshot());
    rec?.record(snap);
    log?.observe(snap);
  }
};

test("every trend channel names a numeric snapshot field and records its value", () => {
  const sim = new SteamTurbineSim();
  const rec = new TrendRecorder();
  run(sim, rec, null, 5);
  const snap = decodeSnapshot(sim.snapshot());
  const k = rec.fine.n - 1;
  for (let c = 0; c < TREND_CHANNELS.length; c++) {
    const chn = TREND_CHANNELS[c];
    assert.equal(typeof snap[chn.key], "number", `${chn.id} -> ${chn.key} is not a snapshot number`);
    const want = Math.fround(snap[chn.key] * (chn.scale ?? 1));
    assert.ok(Object.is(rec.fine.value(k, c), want), `${chn.id}: recorded ${rec.fine.value(k, c)} vs snapshot ${want}`);
  }
  assert.equal(new Set(TREND_CHANNELS.map((c) => c.id)).size, TREND_CHANNELS.length, "channel ids are unique");
});

test("the fine ring keeps the last 60 s at every step and the coarse ring every 0.5 s", () => {
  const sim = new SteamTurbineSim();
  const rec = new TrendRecorder();
  const t0 = decodeSnapshot(sim.snapshot()).timeS; // the engine boots one step in (0.02 s)
  rec.record(decodeSnapshot(sim.snapshot()));
  run(sim, rec, null, 3500); // 70 s of simulated time
  const tEnd = decodeSnapshot(sim.snapshot()).timeS;
  assert.ok(Math.abs(tEnd - t0 - 70) < 1e-6, `engine time ${tEnd} from ${t0}`);
  assert.equal(rec.fine.n, rec.fine.cap, "fine ring is full");
  const fineSpan = rec.fine.time(rec.fine.n - 1) - rec.fine.time(0);
  assert.ok(Math.abs(fineSpan - (rec.fine.cap - 1) * DT) < 1e-6, `fine span ${fineSpan}`);
  for (let k = 1; k < rec.fine.n; k++) {
    const d = rec.fine.time(k) - rec.fine.time(k - 1);
    assert.ok(Math.abs(d - DT) < 1e-9, `fine row ${k} spacing ${d}`);
  }
  // The boot row, then the first step at or after each 0.5 s mark up to t0 + 70 s.
  assert.equal(rec.coarse.n, 1 + Math.floor((tEnd + 1e-9) / COARSE_DT_S), "boot + one row per 0.5 s mark");
  for (let k = 1; k < rec.coarse.n; k++) {
    const mark = k * COARSE_DT_S;
    const t = rec.coarse.time(k);
    assert.ok(t >= mark - 1e-9 && t < mark + DT - 1e-9, `coarse row ${k} at ${t}, mark ${mark}`);
  }
  assert.equal(rec.ringFor(60), rec.fine);
  assert.equal(rec.ringFor(600), rec.coarse, "after the fine ring wraps, long windows read the coarse ring");

  const csv = trendCsv(rec.coarse).trimEnd().split("\n");
  assert.equal(csv.length, 1 + rec.coarse.n);
  assert.equal(csv[0].split(",").length, 1 + TREND_CHANNELS.length);
  assert.ok(csv[0].startsWith("t_s,mw [MW],rpm [rpm]"), csv[0].slice(0, 40));
  const last = csv.at(-1).split(",");
  assert.ok(Math.abs(Number(last[0]) - rec.coarse.time(rec.coarse.n - 1)) < 1e-3, last[0]);
  assert.ok(Math.abs(Number(last[1]) - decodeSnapshot(sim.snapshot()).elecMw) < 0.01, `CSV MW ${last[1]}`);

  // A clock that goes back (a reset nobody announced) starts a new record.
  sim.reset(false);
  rec.record(decodeSnapshot(sim.snapshot()));
  assert.equal(rec.fine.n, 1);
  assert.equal(rec.coarse.n, 1);
});

test("the event log reports the engine's trips, status and breaker changes, and nothing at nominal", () => {
  const sim = new SteamTurbineSim();
  const log = new EventLog();
  log.restart("boot", decodeSnapshot(sim.snapshot()));
  run(sim, null, log, 50);
  assert.deepEqual(log.events.map((e) => e.text), ["boot"], "a steady nominal run logs nothing");

  const tCmd = decodeSnapshot(sim.snapshot()).timeS;
  sim.set_manual_trip(true);
  run(sim, null, log, 10);
  sim.set_grid_breaker(false);
  run(sim, null, log, 10);
  const texts = log.events.map((e) => e.text);
  console.log("  trip events:", JSON.stringify(log.events.slice(1)));
  const trip = log.events.find((e) => e.text === "Raised: MANUAL TRIP");
  assert.ok(trip && trip.kind === "trip", texts.join(" | "));
  assert.ok(Math.abs(trip.t - (tCmd + DT)) < 1e-6, `trip logged at the first step after the command: ${trip.t} vs ${tCmd}`);
  assert.ok(texts.includes("Status SYNC NOMINAL → TRIP: MANUAL ESV"), texts.join(" | "));
  assert.ok(texts.includes("Generator breaker opened"), texts.join(" | "));
  const csv = eventsCsv(log).trimEnd().split("\n");
  assert.equal(csv[0], "t_s,kind,event");
  assert.equal(csv.length, 1 + log.events.length);
});

test("heavy unbalance logs the vibration alarm before the vibration trip, with the measured S(p-p)", () => {
  const sim = new SteamTurbineSim();
  const log = new EventLog();
  log.restart("boot", decodeSnapshot(sim.snapshot()));
  sim.set_rotor_unbalance_g_mm(280000);
  const lim = sim.limits();
  for (let i = 0; i < 400 && !log.events.some((e) => e.text.startsWith("Raised: VIBRATION TRIP")); i++) {
    run(sim, null, log, 1);
  }
  console.log("  unbalance events:", JSON.stringify(log.events.slice(1)));
  const iAlarm = log.events.findIndex((e) => e.text.startsWith("Raised: VIBRATION (S(p-p)"));
  const iTrip = log.events.findIndex((e) => e.text.startsWith("Raised: VIBRATION TRIP"));
  assert.ok(iAlarm > 0 && iTrip > iAlarm, log.events.map((e) => e.text).join(" | "));
  const um = Number(/S\(p-p\) (\d+) µm/.exec(log.events[iAlarm].text)?.[1]);
  assert.ok(um > lim[1], `alarm quotes ${um} µm, above the ${lim[1]} µm alarm level`);
  // The governor trips on the previous step's reading: the trip quotes that reading.
  const umTrip = Number(/S\(p-p\) (\d+) µm/.exec(log.events[iTrip].text)?.[1]);
  assert.ok(umTrip >= lim[2], `trip quotes ${umTrip} µm, at or above the ${lim[2]} µm trip level`);
});
