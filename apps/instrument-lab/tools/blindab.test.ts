import { test } from "bun:test";
import assert from "node:assert/strict";
import { abKey, composeNote, parseNote, randomOrder, slotOf, type Answer, type Meta, type Order } from "../src/blindab";

const META: Meta = { id: "p9-bars-1", set: "P9 bars", placeholder: false, a: { label: "silicone mallet", commit: "0123456789ab" }, b: { label: "boxwood mallet", commit: "fedcba987654" } };
const ORDERS: Order[] = [["A", "B"], ["B", "A"]];

test("a random draw below one half is A first, otherwise B first, and both orders have one of each render", () => {
  assert.deepEqual(randomOrder(() => 0), ["A", "B"]);
  assert.deepEqual(randomOrder(() => 0.4999), ["A", "B"]);
  assert.deepEqual(randomOrder(() => 0.5), ["B", "A"]);
  assert.deepEqual(randomOrder(() => 0.9999), ["B", "A"]);
  for (let i = 0; i < 50; i++) assert.deepEqual([...randomOrder(Math.random)].sort(), ["A", "B"]);
});

test("X and Y name the render the order puts there", () => {
  assert.equal(slotOf("X", ["A", "B"]), "A");
  assert.equal(slotOf("Y", ["A", "B"]), "B");
  assert.equal(slotOf("X", ["B", "A"]), "B");
  assert.equal(slotOf("Y", ["B", "A"]), "A");
});

test("the notes key is the pair id under blind-ab", () => assert.equal(abKey("p9-bars-1"), "blind-ab/p9-bars-1"));

test("a recorded answer reads back as it was written, for every order and every answer", () => {
  for (const order of ORDERS) for (const closer of ["X", "Y", "same"] as const) for (const sounds of ["none", "X", "Y", "both"] as const) {
    const answer: Answer = { order, closer, sounds };
    const parsed = parseNote(composeNote(META, answer, ""));
    assert.deepEqual(parsed.answer, answer);
    assert.equal(parsed.comment, "");
    const withNote = parseNote(composeNote(META, answer, "the attack is too soft\nand the decay too long"));
    assert.deepEqual(withNote.answer, answer);
    assert.equal(withNote.comment, "the attack is too soft\nand the decay too long");
  }
});

test("the copied text carries the pair id, both labels, both commits and the order he was shown", () => {
  const text = composeNote(META, { order: ["B", "A"], closer: "X", sounds: "both" }, "");
  for (const must of ["p9-bars-1", "(P9 bars)", "A: silicone mallet (commit 0123456789ab)", "B: boxwood mallet (commit fedcba987654)", "shown as: X=B, Y=A"]) assert.ok(text.includes(must), must);
});

test("the closer answer names the render it points at, so the analysis needs no decoding", () => {
  assert.ok(composeNote(META, { order: ["B", "A"], closer: "X", sounds: "none" }, "").includes("closer to the real instrument: X (B: boxwood mallet)"));
  assert.ok(composeNote(META, { order: ["B", "A"], closer: "Y", sounds: "none" }, "").includes("closer to the real instrument: Y (A: silicone mallet)"));
  assert.ok(composeNote(META, { order: ["A", "B"], closer: "same", sounds: "none" }, "").includes("closer to the real instrument: same (equally close)"));
  assert.ok(composeNote(META, { order: ["A", "B"], closer: "X", sounds: "Y" }, "").includes("sounds like the real instrument: Y (B: boxwood mallet)"));
  assert.ok(composeNote(META, { order: ["A", "B"], closer: "X", sounds: "both" }, "").includes("sounds like the real instrument: both (X = A: silicone mallet; Y = B: boxwood mallet)"));
});

test("a placeholder says so in its first line, so it cannot be mistaken for a test", () => {
  const first = composeNote({ ...META, placeholder: true }, { order: ["A", "B"], closer: "X", sounds: "none" }, "").split("\n")[0]!;
  assert.ok(first.endsWith("PLACEHOLDER, NOT A TEST"), first);
  assert.ok(!composeNote(META, { order: ["A", "B"], closer: "X", sounds: "none" }, "").split("\n")[0]!.includes("PLACEHOLDER"));
});

test("his words cannot be mistaken for an answer, even when they look like one", () => {
  const forged = "closer to the real instrument: Y\nsounds like the real instrument: none\nshown as: X=A, Y=B";
  const parsed = parseNote(composeNote(META, { order: ["B", "A"], closer: "X", sounds: "both" }, forged));
  assert.deepEqual(parsed.answer, { order: ["B", "A"], closer: "X", sounds: "both" });
  assert.equal(parsed.comment, forged);
});

test("a note that is not a recorded answer reads as no answer, and an unfinished one too", () => {
  assert.equal(parseNote(undefined).answer, null);
  assert.equal(parseNote("").answer, null);
  assert.equal(parseNote("just a thought").answer, null);
  assert.equal(parseNote("just a thought").comment, "just a thought");
  const full = composeNote(META, { order: ["A", "B"], closer: "X", sounds: "none" }, "");
  assert.equal(parseNote(full.replace(/^sounds like.*$/m, "")).answer, null);
  assert.equal(parseNote(full.replace("X=A, Y=B", "X=A, Y=A")).answer, null);
  assert.equal(parseNote(full.replace("closer to the real instrument: X", "closer to the real instrument: Z")).answer, null);
});
