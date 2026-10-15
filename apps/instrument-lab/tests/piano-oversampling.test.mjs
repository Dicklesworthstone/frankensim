// Real factory/worker selection tests. The recording Engine is NOT piano physics.
import test from 'node:test';
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { createPianoEngine } from '../public/live/piano/engine-factory.js';

function recordingEngine() {
  const calls = [];
  class Engine {
    constructor(...args) { calls.push(['new', ...args]); }
    static withBoard(...args) { calls.push(['board', ...args]); return { kind: 'board' }; }
    static withOversampling(...args) { calls.push(['oversampling', ...args]); return { kind: 'oversampling' }; }
    static withBoardOversampling(...args) { calls.push(['board-oversampling', ...args]); return { kind: 'board-oversampling' }; }
  }
  return { Engine, calls };
}

test('omitted or explicit 1x preserves both original constructor signatures', () => {
  const { Engine, calls } = recordingEngine();
  const bytes = new Uint8Array([1, 2, 3]);
  createPianoEngine(Engine, 48000);
  createPianoEngine(Engine, 48000, { substeps: 1, boardCoordinates: 16, courseModeBudget: 80, maxDynamicModes: 900, fullScalePa: 2 });
  createPianoEngine(Engine, 48000, { substeps: 1 }, bytes);
  assert.deepEqual(calls, [['new', 48000, 0, 0, 0, 0], ['new', 48000, 16, 80, 900, 2], ['board', 48000, bytes, 0, 0, 0]]);
});

test('both explicit higher clocks select their actual physical constructor', () => {
  const { Engine, calls } = recordingEngine();
  const bytes = new Uint8Array([4, 5]);
  for (const substeps of [2, 4]) {
    const options = Object.freeze({ substeps, boardCoordinates: 64, courseModeBudget: 200, maxDynamicModes: 2400, fullScalePa: 10 });
    assert.equal(createPianoEngine(Engine, 48000, options).kind, 'oversampling');
    assert.equal(createPianoEngine(Engine, 48000, options, bytes).kind, 'board-oversampling');
    assert.deepEqual(calls.at(-2), ['oversampling', 48000, 0, 200, 2400, 10, substeps]);
    assert.deepEqual(calls.at(-1), ['board-oversampling', 48000, bytes, 200, 2400, 10, substeps]);
    assert.equal(calls.at(-1)[2], bytes);
  }
});

test('invalid or excessive mechanical clocks refuse before engine creation', () => {
  const { Engine, calls } = recordingEngine();
  for (const substeps of [null, 0, -1, 1.5, 3, 8, NaN, Infinity, 'bad']) {
    assert.throws(() => createPianoEngine(Engine, 48000, { substeps }), RangeError);
  }
  for (const [rate, substeps] of [[0, 1], [7999, 1], [48000.5, 1], [48001, 4], [96001, 2], [192001, 1]]) {
    assert.throws(() => createPianoEngine(Engine, rate, { substeps }), RangeError);
  }
  assert.equal(calls.length, 0);
  createPianoEngine(Engine, 96000, { substeps: 2 });
  createPianoEngine(Engine, 44100, { substeps: 4 });
  assert.equal(calls.length, 2);
});

test('old deployed WASM is compatible at 1x but cannot silently drop oversampling', () => {
  const { Engine, calls } = recordingEngine();
  Engine.withOversampling = undefined;
  Engine.withBoardOversampling = undefined;
  createPianoEngine(Engine, 48000);
  for (const bytes of [null, new Uint8Array([1])]) {
    assert.throws(() => createPianoEngine(Engine, 48000, { substeps: 4 }, bytes), /rebuilding.*WASM/);
  }
  assert.equal(calls.length, 1);
});

test('both actual worker factories and their probes forward the selected clock', async () => {
  // Execute the small factory closures FROM the real workers; do not import
  // their WASM or create a background Worker in this Node selection test.
  const worker = await readFile(new URL('../public/live/piano/worker.js', import.meta.url), 'utf8');
  const lane = await readFile(new URL('../public/live/piano/lane-worker.js', import.meta.url), 'utf8');
  const { Engine, calls } = recordingEngine();
  const makeBody = worker.match(/function makeEngine\(options, bytes = boardBytes\) \{([\s\S]*?)\n\}/)?.[1];
  assert.ok(makeBody, 'standalone factory remains inspectable');
  const make = new Function('createPianoEngine', 'PianoEngine', 'sampleRate', 'options', 'bytes', makeBody);
  make(createPianoEngine, Engine, 48000, { substeps: 4 }, null);
  assert.deepEqual(calls.at(-1), ['oversampling', 48000, 0, 0, 0, 0, 4]);
  const laneBody = lane.match(/const make = \(bytes, coords, cap\) => \{([\s\S]*?)\n      \};/)?.[1];
  assert.ok(laneBody, 'scheduled-lane factory remains inspectable');
  const makeLane = new Function('createPianoEngine', 'PianoEngine', 'o', 'full', 'atCap', 'bytes', 'coords', 'cap', laneBody);
  const board = new Uint8Array([9]);
  for (const substeps of [2, 4]) {
    for (const bytes of [null, board]) {
      const caps = [];
      makeLane(createPianoEngine, Engine, { substeps }, {}, (e, cap) => { caps.push(cap); return e; }, bytes, 16, 900);
      assert.deepEqual(calls.at(-1), [bytes ? 'board-oversampling' : 'oversampling', 48000, bytes ?? 16, 0, 0, 0, substeps]);
      assert.deepEqual(caps, [900]);
    }
  }
  assert.match(worker, /makeEngine\(\{ \.\.\.engineOptions, maxDynamicModes: 1e9 \}\)/);
  assert.match(lane, /chooseQuality\(\{ make:/);
});
