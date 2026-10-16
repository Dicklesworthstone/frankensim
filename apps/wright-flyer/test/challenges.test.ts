// E9.2c challenge battery: exact preset IDs/queries, one shared route
// into ScenarioInit, fail-closed lookup, and mutation-sensitive config
// identities. No local descriptor is mislabeled PhysicalScenarioId.
// Repro: node --test test/challenges.test.ts

import assert from "node:assert/strict";
import { test } from "node:test";

import {
  CHALLENGE_PRESETS,
  challengeById,
  challengeQuery,
  challengeScenario,
  challengeScenarioIdentity,
  scenarioConfigIdentity,
} from "../src/challenges.ts";
import { scenarioFromQuery } from "../src/menu.ts";
import { MODE_FIXED, MODE_HISTORICAL, MODE_HUMAN } from "../src/sim/protocol.ts";

test("challenge rail is an exact, unique, honest preset catalog", () => {
  assert.deepEqual(
    CHALLENGE_PRESETS.map((preset) => preset.id),
    [
      "assisted-first-flight",
      "authentic-fourth-flight",
      "modeled-fourth-flight",
      "huffman-catapult",
    ],
  );
  assert.equal(new Set(CHALLENGE_PRESETS.map((preset) => preset.id)).size, 4);
  for (const preset of CHALLENGE_PRESETS) {
    assert.ok(preset.title.length > 0 && preset.description.length > 0);
    assert.ok(!/will (complete|reach|fly)/i.test(preset.description), preset.id);
  }
});

test("challenge queries use only the existing scenario route", () => {
  assert.deepEqual(
    CHALLENGE_PRESETS.map((preset) => challengeQuery(preset)),
    [
      "?sim=1&mode=human&assist=1&flight=1",
      "?sim=1&mode=human&flight=4",
      "?sim=1&mode=historical&flight=4",
      "?sim=1&site=huffman",
    ],
  );
  const allowed = new Set(["sim", "mode", "assist", "flight", "site"]);
  for (const preset of CHALLENGE_PRESETS) {
    const params = new URLSearchParams(challengeQuery(preset).slice(1));
    for (const key of params.keys()) {
      assert.ok(allowed.has(key), `${preset.id} invented query field ${key}`);
    }
    assert.deepEqual(scenarioFromQuery(params), challengeScenario(preset));
  }
});

test("each preset loads a distinct deterministic ScenarioInit identity", () => {
  const scenarios = CHALLENGE_PRESETS.map(challengeScenario);
  assert.deepEqual(
    scenarios.map((scenario) => scenario.mode),
    [MODE_HUMAN, MODE_HUMAN, MODE_HISTORICAL, MODE_FIXED],
  );
  assert.deepEqual(
    scenarios.map((scenario) => scenario.assist),
    [true, false, false, false],
  );
  assert.deepEqual(
    scenarios.map((scenario) => scenario.catapult),
    [false, false, false, true],
  );
  const identities = CHALLENGE_PRESETS.map(challengeScenarioIdentity);
  assert.equal(new Set(identities).size, CHALLENGE_PRESETS.length);
  assert.deepEqual(identities, CHALLENGE_PRESETS.map(challengeScenarioIdentity));
  for (const identity of identities) {
    const parsed = JSON.parse(identity);
    assert.equal(parsed.schema, "wf-scenario-init-v1");
    assert.equal(typeof parsed.seed, "string");
    assert.equal(typeof parsed.maxTicks, "string");
  }
});

test("identity covers every field and changes when a field changes", () => {
  const base = challengeScenario(CHALLENGE_PRESETS[0]!);
  const original = scenarioConfigIdentity(base);
  const variants = [
    { ...base, seed: base.seed + 1n },
    { ...base, rhoKgM3: base.rhoKgM3 + 0.01 },
    { ...base, headwindMps: base.headwindMps + 0.01 },
    { ...base, mode: base.mode + 1 },
    { ...base, member: base.member + 1 },
    { ...base, railLengthM: base.railLengthM + 0.01 },
    { ...base, maxTicks: base.maxTicks + 1n },
    { ...base, assist: !base.assist },
    { ...base, catapult: !base.catapult },
  ];
  for (const variant of variants) {
    assert.notEqual(scenarioConfigIdentity(variant), original);
  }
});

test("lookup fails closed for unknown challenge ids", () => {
  for (const preset of CHALLENGE_PRESETS) {
    assert.equal(challengeById(preset.id), preset);
  }
  assert.equal(challengeById(""), null);
  assert.equal(challengeById("assisted-first-flight "), null);
  assert.equal(challengeById("unknown"), null);
});

test("hostile URL combinations cannot leak assist or missions into unsupported modes/sites", () => {
  const fixedAssist = scenarioFromQuery(new URLSearchParams("mode=fixed&assist=1&flight=1"));
  assert.equal(fixedAssist.assist, false);
  const huffman = scenarioFromQuery(
    new URLSearchParams("mode=human&site=huffman&assist=1&flight=4"),
  );
  assert.equal(huffman.assist, false);
  assert.equal(huffman.catapult, true);
  assert.equal(huffman.seed, 1903n, "Dec-17 mission identity cannot leak to Huffman");
  const invalidFlight = scenarioFromQuery(new URLSearchParams("mode=human&flight=99"));
  assert.equal(invalidFlight.seed, 1903n);
});
