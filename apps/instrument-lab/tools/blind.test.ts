import { test } from "bun:test";
import assert from "node:assert/strict";
import { blindLine, blindOrder, letterOf, type Order } from "../src/blind";

const KEYS = ["clarinet/A4", "clarinet/D3", "violin/G3", "violin/E5", "cello/C2", "voice/C3-a", "voice/C3-e", "trombone/Bb2", "guitar/B3", "viola/D4"];

test("a note opens in the same order every time, until it is shuffled", () => {
  for (const k of KEYS) assert.deepEqual(blindOrder(k, 0), blindOrder(k, 0));
});

test("every order has exactly one real clip and one render", () => {
  for (const k of KEYS) for (let n = 0; n < 4; n++) assert.deepEqual([...blindOrder(k, n)].sort(), ["ours", "real"]);
});

test("across notes and across shuffles both orders occur, and about half each", () => {
  const first = (k: string, n: number): string => blindOrder(k, n)[0];
  const orders = Array.from({ length: 400 }, (_, i) => first(`family${i % 7}/note${i}`, i % 5));
  const realFirst = orders.filter((o) => o === "real").length;
  assert.ok(realFirst > 150 && realFirst < 250, `real was X ${realFirst} times in 400`);
  for (const k of KEYS) assert.ok(new Set(Array.from({ length: 24 }, (_, n) => first(k, n))).size === 2, `${k} never changes under Shuffle again`);
});

test("the stored line names the pick and where the real clip was", () => {
  const realIsX: Order = ["real", "ours"];
  const realIsY: Order = ["ours", "real"];
  assert.equal(blindLine(realIsX, 0), "blind: picked X, real was X");
  assert.equal(blindLine(realIsX, 1), "blind: picked Y, real was X");
  assert.equal(blindLine(realIsY, 0), "blind: picked X, real was Y");
  assert.equal(blindLine(realIsY, 1), "blind: picked Y, real was Y");
  assert.equal(letterOf("real", realIsY), "Y");
  assert.equal(letterOf("ours", realIsY), "X");
});
