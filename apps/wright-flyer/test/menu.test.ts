// Landing-menu battery (UI overhaul, root guzez): the selection->URL
// map emits exactly the params the app honors — no third config path,
// assist never leaks into modes that cannot engage it, and the key
// card mirrors input.ts bindings.
// Repro: node --test test/menu.test.ts

import assert from "node:assert/strict";
import { test } from "node:test";

import {
  DEFAULT_SELECTION,
  KEY_LINES,
  MODE_CARDS,
  assistAvailable,
  menuQuery,
} from "../src/menu.ts";
import { keysFrom } from "../src/input.ts";
import { journeyNextUrl, journeyStage, JOURNEY_STAGES } from "../src/journey.ts";

function jlog(kase: string, payload: string): void {
  console.log(`{"suite":"wf-app-menu","case":"${kase}",${payload}}`);
}

test("selection to query, exact strings", () => {
  assert.equal(menuQuery({ mode: "human", site: "kdh", assist: false }), "?sim=1&mode=human");
  assert.equal(
    menuQuery({ mode: "human", site: "kdh", assist: true }),
    "?sim=1&mode=human&assist=1",
  );
  assert.equal(
    menuQuery({ mode: "human", site: "huffman", assist: true }),
    "?sim=1&mode=human&site=huffman",
    "Huffman pins assist off (protocol.ts), so the URL drops it",
  );
  assert.equal(menuQuery({ mode: "historical", site: "kdh", assist: false }), "?sim=1&mode=historical");
  assert.equal(menuQuery({ mode: "fixed", site: "kdh", assist: false }), "?sim=1");
  assert.equal(menuQuery({ mode: "fixed", site: "huffman", assist: false }), "?sim=1&site=huffman");
  jlog("exact-queries", `"count":5`);
});

test("assist never leaks where it cannot engage", () => {
  assert.equal(
    menuQuery({ mode: "historical", site: "kdh", assist: true }),
    "?sim=1&mode=historical",
    "historical ignores assist",
  );
  assert.equal(menuQuery({ mode: "fixed", site: "huffman", assist: true }), "?sim=1&site=huffman");
  // huffmanScenario pins assist OFF (protocol.ts), so even HUMAN mode
  // at Huffman must not emit the param — the URL would claim an aid
  // the scenario silently drops.
  assert.equal(
    menuQuery({ mode: "human", site: "huffman", assist: true }),
    "?sim=1&mode=human&site=huffman",
    "human at Huffman drops assist (scenario cannot engage it)",
  );
  assert.equal(assistAvailable({ mode: "human", site: "kdh" }), true);
  assert.equal(assistAvailable({ mode: "human", site: "huffman" }), false);
  assert.equal(assistAvailable({ mode: "historical", site: "kdh" }), false);
  jlog("assist-gating", `"leak":false`);
});

test("mode cards cover all three modes exactly once", () => {
  const modes = MODE_CARDS.map((c) => c.mode).sort();
  assert.deepEqual(modes, ["fixed", "historical", "human"]);
  for (const card of MODE_CARDS) {
    assert.ok(card.title.length > 0 && card.blurb.length > 0);
  }
  jlog("mode-cards", `"count":${MODE_CARDS.length}`);
});

test("default selection is a flyable human start", () => {
  assert.equal(DEFAULT_SELECTION.mode, "human");
  assert.equal(menuQuery(DEFAULT_SELECTION), "?sim=1&mode=human");
  jlog("default", `"query":"?sim=1&mode=human"`);
});

test("key card mirrors input.ts bindings (the ONLY binding authority)", () => {
  // The card claims S/↓ pulls the canard: input.ts must agree.
  const pull = keysFrom(new Set(["KeyS"]));
  assert.equal(pull.canardUp, true, "S pulls (nose up)");
  const pullArrow = keysFrom(new Set(["ArrowDown"]));
  assert.equal(pullArrow.canardUp, true);
  const push = keysFrom(new Set(["KeyW"]));
  assert.equal(push.canardDown, true);
  const left = keysFrom(new Set(["ArrowLeft"]));
  assert.equal(left.warpLeft, true);
  const recenter = keysFrom(new Set(["Space"]));
  assert.equal(recenter.recenter, true);
  // Every claim in the card names a real binding word.
  const joined = KEY_LINES.join(" ");
  for (const word of ["canard", "warp", "recenter", "replay", "telemetry"]) {
    assert.ok(joined.includes(word), `card mentions ${word}`);
  }
  jlog("key-card", `"lines":${KEY_LINES.length}`);
});

test("flight param emits only for valid missions at Kill Devil Hills", () => {
  assert.equal(menuQuery({ mode: "human", site: "kdh", assist: false, flight: 4 }), "?sim=1&mode=human&flight=4");
  assert.equal(menuQuery({ mode: "human", site: "kdh", assist: true, flight: 1 }), "?sim=1&mode=human&assist=1&flight=1");
  // Huffman drops the mission (the URL never claims what will not run).
  assert.equal(
    menuQuery({ mode: "historical", site: "huffman", assist: false, flight: 2 }),
    "?sim=1&mode=historical&site=huffman",
  );
  // Out-of-range ids emit nothing (no third config path, no junk params).
  assert.equal(menuQuery({ mode: "fixed", site: "kdh", assist: false, flight: 0 }), "?sim=1");
  assert.equal(menuQuery({ mode: "fixed", site: "kdh", assist: false, flight: 9 }), "?sim=1");
  jlog("flight-gating", `"cases":5`);
});

test("journey stages chain watch -> assist -> authentic with honest copy", () => {
  assert.equal(JOURNEY_STAGES.length, 3);
  assert.equal(journeyStage(null), null);
  assert.equal(journeyStage("0"), null);
  assert.equal(journeyStage("4"), null);
  assert.equal(journeyStage("abc"), null);
  const s2 = journeyStage("2")!;
  assert.match(s2.url, /mode=human&assist=1&journey=2/);
  // Copy law: stage 1 names itself a hypothesis; no stage promises an outcome.
  assert.match(journeyStage("1")!.caption, /hypothesis/);
  for (const s of JOURNEY_STAGES) {
    assert.ok(!/will (fly|reach|travel)/i.test(s.caption), `no promises in stage ${s.index}`);
    for (const other of JOURNEY_STAGES) {
      if (other.index === s.index + 1) {
        assert.equal(journeyNextUrl(s.index), other.url, `stage ${s.index} chains to ${other.index}`);
      }
    }
  }
  assert.equal(journeyNextUrl(3), null, "stage 3 ends the journey");
  jlog("journey-chain", `"stages":${JOURNEY_STAGES.length}`);
});
