import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { test } from 'node:test';

// Exercise the production scheduler, with controlled clock/render boundaries.
// Sound and real Worker behavior are checked separately in the browser.
const source = readFileSync(new URL('../public/live/piano/worker.js', import.meta.url), 'utf8');
const body = source.slice(source.indexOf('function estimateQueue(now) {'), source.indexOf('function renderChunk() {'));
type Level = { at: number; queuedFrames: number; receivedFrames: number };
const estimate = new Function('now', 'level', 'sentTotal', 'sampleRate', `const CHUNK = 128; ${body}; return estimateQueue(now);`) as
  (now: number, level: Level, sent: number, rate: number) => number;

test('queue prediction retains normal latency and in-flight accounting', () => {
  assert.equal(estimate(102, { at: 100, queuedFrames: 896, receivedFrames: 1024 }, 1152, 48000), 928);
});

test('refill accounts for consumption while rendering without enlarging the queue target', () => {
  const pump = source.slice(source.indexOf('function pump() {'), source.indexOf('function onSinkMessage(event) {'));
  const run = new Function(`
    const CHUNK = 128, sampleRate = 48000, targetFrames = 864;
    let now = 100, sentTotal = 512, chunks = 0;
    const level = { at: 100, queuedFrames: 512, receivedFrames: 512 };
    const performance = { now: () => now };
    const engine = { lastError: () => null }, sinkPort = {};
    const snapshotIntervalMs = Infinity, lastSnapshotAt = 0, lastStatsAt = Infinity;
    const govern = () => {};
    function renderChunk() { now += 3; sentTotal += CHUNK; chunks++; }
    ${body}
    ${pump}
    pump();
    return { chunks, queued: estimateQueue(now) };
  `) as () => { chunks: number; queued: number };
  // The former loop stopped after three chunks, leaving the queue short.
  assert.deepEqual(run(), { chunks: 5, queued: 896 });
});

test('a stalled sink cannot accumulate seconds of stale audio; acknowledgement resumes it', () => {
  for (const rate of [44100, 48000, 96000, 192000]) {
    const target = Math.ceil(rate * .018 / 128) * 128;
    let sent = target;
    const level = { at: 1, queuedFrames: target, receivedFrames: target };
    for (let now = 9; now < 30000; now += 8) {
      while (estimate(now, level, sent, rate) < target) sent += 128;
    }
    assert.equal(sent, target + 256);
    const resumed = { at: 30001, queuedFrames: 0, receivedFrames: sent };
    assert.equal(estimate(30001, resumed, sent, rate), 0);
    assert.ok(estimate(30002, resumed, sent + target, rate) > 0);
  }
});

test('drawing yields under audio load, stays bounded, and resumes without changing audio', () => {
  const pump = source.slice(source.indexOf('function pump() {'), source.indexOf('function onSinkMessage(event) {'));
  const run = new Function(`
    const CHUNK = 128, sampleRate = 48000, targetFrames = 864;
    let now = 100, sentTotal = 2048, loadAverage = 0.9;
    let lastSnapshotAt = 0, snapshotIntervalMs = 1000 / 30;
    const lastStatsAt = Infinity, rendered = 2048;
    const level = { at: 1, queuedFrames: 2048, receivedFrames: 2048 };
    const performance = { now: () => now }, sinkPort = {};
    const updates = [];
    const post = message => updates.push(message.type);
    const engine = { lastError: () => null, snapshot: () => new Float32Array(1) };
    const govern = () => {};
    const renderChunk = () => { throw new Error('already full; must not render more audio'); };
    ${body}
    ${pump}
    pump(); const deferred = updates.length;
    now = 250; pump(); const bounded = updates.length;
    now = 300; loadAverage = 0.2; pump(); const resumed = updates.length;
    now = 1000; snapshotIntervalMs = 0; pump();
    return { deferred, bounded, resumed, final: updates.length, sentTotal, rendered };
  `);
  assert.deepEqual(run(), { deferred: 0, bounded: 1, resumed: 2, final: 2, sentTotal: 2048, rendered: 2048 });
});

const governSource = source.slice(source.indexOf('function govern(now) {'), source.indexOf('function pump() {'));
type Quality = { courseModeBudget: number; partialBandHz: number; maxDynamicModes: number };
const governorRun = new Function('asked', 'initial', 'loads', 'admit', `
  let granted = { ...initial }, governor = true, lastGovernorAt = 0, governorSteps = 0, loadAverage = 0;
  const calls = [];
  const engine = { setQuality: (...args) => { calls.push(args); return admit; } };
  ${governSource}
  let now = 400;
  for (const load of loads) { loadAverage = load; govern(now); now += 400; }
  return { granted, governorSteps, calls };
`) as (asked: Quality, initial: Quality, loads: number[], admit: boolean) => {
  granted: Quality; governorSteps: number; calls: number[][];
};

test('governor cannot raise a user-selected low course budget under overload', () => {
  const requested = { courseModeBudget: 32, partialBandHz: 12000, maxDynamicModes: 900 };
  const result = governorRun(requested, requested, [1, 1, 1, 1, 1, 1], true);
  assert.equal(result.granted.courseModeBudget, 32);
  assert.equal(result.granted.maxDynamicModes, 400);
  assert.ok(result.calls.every(([modes, , total]) => modes <= 32 && total <= 900));
});

test('governor commits no state or reduction counter when the engine refuses', () => {
  const requested = { courseModeBudget: 200, partialBandHz: 12000, maxDynamicModes: 2400 };
  const reduced = { ...requested, courseModeBudget: 100, maxDynamicModes: 900 };
  for (const [state, load] of [[requested, 1], [reduced, 0.2]] as const) {
    const result = governorRun(requested, state, [load], false);
    assert.deepEqual(result.granted, state);
    assert.equal(result.governorSteps, 0);
    assert.equal(result.calls.length, 1);
  }
});

test('governor preserves ordinary reduction and recovery decisions', () => {
  const requested = { courseModeBudget: 200, partialBandHz: 12000, maxDynamicModes: 2400 };
  const result = governorRun(requested, requested, [1, 0.2, 0.2, 0.2], true);
  assert.deepEqual(result.calls, [[170, 12000, 1920], [184, 12000, 2112], [199, 12000, 2324], [200, 12000, 2400]]);
  assert.deepEqual(result.granted, requested);
  assert.equal(result.governorSteps, 1);
});
