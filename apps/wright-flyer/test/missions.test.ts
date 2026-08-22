// Mission battery (PurpleCliff slice, wf-root-guez): flight table
// integrity vs the plan §3.2 anchors, seed separation + determinism,
// and the honesty law of missionOutcome (distribution language, never
// a "beat this" target; every verdict names modeled-vs-accounts).
// Repro: node --test test/missions.test.ts

import assert from "node:assert/strict";
import { test } from "node:test";

import {
  DEC17_FLIGHTS,
  flightByIndex,
  flightSeed,
  missionOutcome,
} from "../src/missions/flights.ts";

test("the four Dec-17 anchors are present with pilots and bands", () => {
  assert.equal(DEC17_FLIGHTS.length, 4);
  assert.deepEqual(
    DEC17_FLIGHTS.map((f) => [f.id, f.pilot]),
    [
      [1, "Orville"],
      [2, "Wilbur"],
      [3, "Orville"],
      [4, "Wilbur"],
    ],
  );
  for (const f of DEC17_FLIGHTS) {
    assert.ok(f.lowM > 0 && f.highM >= f.lowM, `band sane for flight ${f.id}`);
    assert.ok(f.durationS > 0 && f.note.length > 10);
    // Metre bands must agree with the feet figures in the notes.
    if (f.id === 4) {
      assert.equal(f.precision, "measured");
      assert.ok(f.highM - f.lowM < 2, "a surveyed flight is a narrow tolerance band");
      assert.ok(Math.abs((f.lowM + f.highM) / 2 - 259.7) < 0.2, "852 ft = 259.7 m midpoint");
    } else {
      assert.equal(f.precision, "accounts");
      assert.ok(f.highM > f.lowM, "account variance is a real band");
    }
  }
});

test("flight lookup refuses non-integer and out-of-range ids", () => {
  for (let id = 1; id <= 4; id++) {
    assert.equal(flightByIndex(id)?.id, id);
  }
  assert.equal(flightByIndex(0), null);
  assert.equal(flightByIndex(5), null);
  assert.equal(flightByIndex(2.5), null);
});

test("per-flight seeds differ, stay deterministic, and keep the base case", () => {
  const seeds = DEC17_FLIGHTS.map((f) => flightSeed(f.id));
  assert.equal(new Set(seeds).size, 4, "four distinct ensembles");
  for (const s of seeds) {
    assert.ok(s !== 1903n, "missions never reuse the default seed");
  }
  assert.deepEqual(seeds, DEC17_FLIGHTS.map((f) => flightSeed(f.id)), "deterministic");
});

test("outcome wording obeys the honesty law in all three regions", () => {
  const f4 = flightByIndex(4)!;
  // Within band.
  const within = missionOutcome(f4, 259.8, 52.0);
  // Short of band — explains why, never scolds.
  const short = missionOutcome(f4, 100.0, 20.0);
  assert.match(short.verdict, /SHORT/);
  assert.match(short.verdict, /conditions and your inputs differ/);
  // Beyond band — explicitly NOT a record claim.
  const beyond = missionOutcome(f4, 400.0, 70.0);
  assert.match(beyond.verdict, /BEYOND/);
  assert.match(beyond.verdict, /not a record claim/);
  // Every outcome carries the not-a-scoreboard footer + run numbers.
  for (const o of [within, short, beyond]) {
    assert.equal(o.lines.length, 5);
    assert.match(o.lines[4]!, /not a scoreboard/);
    assert.match(o.lines[1]!, /your run: \d+\.\d m downrange/);
    assert.match(o.lines[0]!, /MISSION: flight 4 · Wilbur/);
  }
  // Accounts-precision flights print a wide BAND, measured prints its
  // narrow surveyed tolerance band.
  const f1 = flightByIndex(1)!;
  assert.match(missionOutcome(f1, 33.0, 12.0).lines[2]!, /historical band: 30\.5–36\.6 m/);
  assert.match(missionOutcome(f4, 260.0, 59.0).lines[2]!, /historical band: 258\.8–260\.6 m \(surveyed\)/);
});
