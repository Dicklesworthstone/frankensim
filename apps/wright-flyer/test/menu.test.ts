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
