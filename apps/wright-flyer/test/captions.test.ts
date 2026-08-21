// E5.3b-ii caption battery (bead wf-root-guzez.6.5.2): the label gate
// (unlabeled or double-labeled claims REFUSE), event-detection oracles
// on synthetic snapshot sequences (liftoff tick exact, undulation
// count, every terminal's caption + label), the porpoises-view link on
// undulation captions, and idempotence past the terminal.
// Repro: node --test test/captions.test.ts

import assert from "node:assert/strict";
import { test } from "node:test";

import { CaptionStream, formatCaption } from "../src/sim/captions.ts";
import type { SimSnapshot } from "../src/sim/snapshotView.ts";

function snap(overrides: Partial<SimSnapshot>): SimSnapshot {
  return {
    tick: 1,
    phase: "on-rail",
    ended: false,
    xM: 0,
    hM: 2.4,
    uMps: 11,
    wMps: 0,
    qRadS: 0,
    thetaRad: 0,
    dcRad: 0.14,
    warpRad: 0,
    omegaPropRadS: 48,
    gustWMps: 0,
    assistActive: false,
    ...overrides,
  };
}

test("format: label bracket present; empty and double-labeled refuse", () => {
  const line = formatCaption({ atTick: 1, label: "Hypothesis", text: "the pitch mode." });
  assert.equal(line, "[Hypothesis] the pitch mode.");
  const linked = formatCaption({
    atTick: 1,
    label: "Hypothesis",
    text: "undulation.",
    link: "why-it-porpoises",
  });
  assert.match(linked, /→ why-it-porpoises$/);
  assert.throws(() => formatCaption({ atTick: 1, label: "Verified", text: "" }), RangeError);
  assert.throws(
    () => formatCaption({ atTick: 1, label: "Verified", text: "[Verified] smuggled" }),
    RangeError,
  );
});

test("stream: rail, liftoff tick exact, undulations counted, ground contact", () => {
  const cs = new CaptionStream();
  cs.feed(snap({ tick: 10 }));
  cs.feed(snap({ tick: 11 }));
  cs.feed(snap({ tick: 626, phase: "airborne", qRadS: 0.4 }));
  // Two q sign flips = ONE undulation (captioned on the even flip).
  cs.feed(snap({ tick: 700, phase: "airborne", qRadS: -0.4 }));
  cs.feed(snap({ tick: 800, phase: "airborne", qRadS: 0.4 }));
  cs.feed(
    snap({ tick: 1450, phase: "ended:ground-contact", ended: true, xM: 46.7 }),
  );
  const all = cs.all();
  assert.equal(all.length, 4, JSON.stringify(all));
  assert.equal(all[0]!.atTick, 10, "rail caption at first feed");
  assert.equal(all[0]!.label, "Verified");
  assert.equal(all[1]!.atTick, 626, "liftoff at the exact transition tick");
  assert.equal(all[1]!.label, "Estimated");
  assert.equal(all[2]!.label, "Hypothesis", "undulation attribution is a hypothesis");
  assert.equal(all[2]!.link, "why-it-porpoises");
  assert.equal(cs.undulations(), 1);
  assert.match(all[3]!.text, /Ground contact at 12\.1 s, 47 m/);
  // Past-terminal feeds are ignored (idempotent tail).
  cs.feed(snap({ tick: 1451, phase: "ended:ground-contact", ended: true }));
  assert.equal(cs.all().length, 4);
});

test("stream: every terminal phase gets a labeled caption", () => {
  const terminals: [SimSnapshot["phase"], RegExp][] = [
    ["ended:rail-end-without-lift", /Ran off the rail/],
    ["ended:envelope-exceeded", /certified aero envelope/],
    ["ended:max-ticks", /Tick budget/],
  ];
  for (const [phase, re] of terminals) {
    const cs = new CaptionStream();
    cs.feed(snap({ tick: 5 }));
    cs.feed(snap({ tick: 6, phase, ended: true }));
    const last = cs.all().at(-1)!;
    assert.match(last.text, re);
    assert.ok(["Verified", "Estimated", "Hypothesis"].includes(last.label));
    // The formatter accepts every emitted caption (the gate holds).
    assert.doesNotThrow(() => formatCaption(last));
  }
});

test("upTo filters by tick for the HUD tail", () => {
  const cs = new CaptionStream();
  cs.feed(snap({ tick: 10 }));
  cs.feed(snap({ tick: 626, phase: "airborne", qRadS: 0.3 }));
  assert.equal(cs.upTo(20).length, 1);
  assert.equal(cs.upTo(626).length, 2);
});
