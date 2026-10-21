import { test } from "bun:test";
import assert from "node:assert/strict";
import { candidateKey, commentOf, compose, tally, verdictOf, withComment, withVerdict } from "../src/verdict";

test("a note carries its verdict on the first line and nothing else is read as one", () => {
  assert.equal(verdictOf("verdict: VALID"), "VALID");
  assert.equal(verdictOf("verdict: INVALID\nthe room is too wet"), "INVALID");
  assert.equal(verdictOf("verdict: VALID\r\nwindows line ending"), "VALID");
  assert.equal(verdictOf(undefined), null);
  assert.equal(verdictOf(""), null);
  assert.equal(verdictOf("I think this is verdict: VALID"), null); // a sentence of his own, not a mark
  assert.equal(verdictOf("a line first\nverdict: VALID"), null);
  assert.equal(verdictOf("verdict: valid"), null);
  assert.equal(verdictOf("verdict: VALIDATED"), null);
});

test("the comment is everything but the verdict line", () => {
  assert.equal(commentOf("verdict: VALID\nclose-miked, dry"), "close-miked, dry");
  assert.equal(commentOf("verdict: INVALID"), "");
  assert.equal(commentOf("no verdict, just words"), "no verdict, just words");
  assert.equal(commentOf(undefined), "");
});

test("setting a verdict keeps the comment, and a comment keeps the verdict", () => {
  const note = "reverb tail after each note";
  assert.equal(withVerdict(note, "INVALID"), "verdict: INVALID\nreverb tail after each note");
  assert.equal(withVerdict("verdict: INVALID\nreverb tail after each note", "VALID"), "verdict: VALID\nreverb tail after each note");
  assert.equal(withVerdict("verdict: VALID\nx", null), "x");
  assert.equal(withComment("verdict: VALID", "dry and close"), "verdict: VALID\ndry and close");
  assert.equal(withComment("verdict: VALID\nold", ""), "verdict: VALID");
  assert.equal(withComment(undefined, "words"), "words");
  assert.equal(withComment("old words", ""), "");
});

test("clearing a verdict and a comment leaves an empty note the page can drop", () => {
  assert.equal(compose(null, ""), "");
  assert.equal(compose("VALID", ""), "verdict: VALID");
  assert.equal(compose("VALID", "   \n"), "verdict: VALID");
});

test("a comment that itself starts with a verdict word is stored as a comment, under the real verdict", () => {
  const stored = withComment("verdict: VALID", "verdict: INVALID is what I wrote first");
  assert.equal(stored.split("\n")[0], "verdict: VALID");
  assert.equal(verdictOf(stored), "VALID");
});

test("the tally counts only the candidates asked about", () => {
  const notes: Record<string, string> = {
    [candidateKey("violin-yt-a")]: "verdict: VALID",
    [candidateKey("violin-yt-b")]: "verdict: INVALID\nreverb",
    [candidateKey("violin-yt-c")]: "just a remark",
    [candidateKey("cello-yt-z")]: "verdict: VALID",
    "violin/G3": "verdict: VALID", // a card note that happens to look the same: not a candidate
  };
  assert.deepEqual(tally(["violin-yt-a", "violin-yt-b", "violin-yt-c", "violin-yt-d"], notes), { valid: 1, invalid: 1, undecided: 2 });
  assert.deepEqual(tally([], notes), { valid: 0, invalid: 0, undecided: 0 });
});
