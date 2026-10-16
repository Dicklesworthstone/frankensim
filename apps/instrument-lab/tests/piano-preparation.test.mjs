// Exercise the production constructor selector. The recording Engine does NOT
// parse Rust scales, run WASM, or simulate piano physics.
import test from 'node:test';
import assert from 'node:assert/strict';
import { createPianoEngine } from '../public/live/piano/engine-factory.js';

function recording() {
  const calls = [];
  class Engine {
    constructor(...args) { calls.push(['legacy', ...args]); }
    static withBoard(...args) { calls.push(['board', ...args]); return { kind: 'board' }; }
    static withOversampling(...args) { calls.push(['rate', ...args]); return { kind: 'rate' }; }
    static withBoardOversampling(...args) { calls.push(['board-rate', ...args]); return { kind: 'board-rate' }; }
    static withPreparation(...args) { calls.push(['prepared', ...args]); return { kind: 'prepared' }; }
  }
  return { Engine, calls };
}
const scale = '# supplied SI scale text is passed to Rust, never interpreted here\n';

test('import alone forwards exact text and an empty tuning packet', () => {
  const { Engine, calls } = recording();
  const options = Object.freeze({ scaleCsv: scale, boardCoordinates: 16 });
  assert.equal(createPianoEngine(Engine, 48000, options).kind, 'prepared');
  assert.deepEqual(calls[0], ['prepared', 48000, new Uint8Array(), 16, 0, 0, 0, 1, scale, new Float64Array()]);
  assert.equal(options.scaleCsv, scale);
});

test('complete preparation reaches both board paths at all admitted clocks', () => {
  const { Engine, calls } = recording();
  const board = new Uint8Array([4, 5, 6]);
  const bass = Object.freeze([0, 1, 1, 0]), treble = Object.freeze([1, 0.2, 0, 0]);
  const stretchTuning = Object.freeze({ concertAHz: 442, bassOctaveWeights: bass,
    trebleOctaveWeights: treble, maxPartialHz: 20000 });
  for (const substeps of [1, 2, 4]) {
    for (const bytes of [null, board]) {
      const options = Object.freeze({ scaleCsv: scale, stretchTuning, substeps,
        boardCoordinates: 64, courseModeBudget: 180, maxDynamicModes: 2100, fullScalePa: 7 });
      createPianoEngine(Engine, 48000, options, bytes);
      assert.deepEqual(calls.at(-1), ['prepared', 48000, bytes ?? new Uint8Array(),
        0, 180, 2100, 7, substeps, scale,
        new Float64Array([442, 0,1,1,0, 1,0.2,0,0, 20000])]);
      if (bytes) assert.equal(calls.at(-1)[2], board);
    }
  }
  createPianoEngine(Engine, 48000, { stretchTuning: { concertAHz: 440 } });
  assert.deepEqual(calls.at(-1).at(-1), new Float64Array([440, 1,0,0,0, 1,0,0,0, 0]));
  assert.equal(calls.at(-1).at(-2), ''); // embedded scale, not a fabricated CSV
  assert.deepEqual(bass, [0,1,1,0]); assert.deepEqual(treble, [1,0.2,0,0]);
});

test('absent preparation preserves old constructors and old-WASM behavior', () => {
  const { Engine, calls } = recording();
  const board = new Uint8Array([9]);
  Engine.withPreparation = undefined;
  createPianoEngine(Engine, 48000);
  createPianoEngine(Engine, 48000, { substeps: 1, boardCoordinates: 16 }, board);
  createPianoEngine(Engine, 48000, { substeps: 4 });
  createPianoEngine(Engine, 48000, { substeps: 2 }, board);
  assert.deepEqual(calls, [['legacy',48000,0,0,0,0], ['board',48000,board,0,0,0],
    ['rate',48000,0,0,0,0,4], ['board-rate',48000,board,0,0,0,2]]);
  for (const bytes of [null, board]) {
    for (const options of [{ scaleCsv: scale }, { stretchTuning: { concertAHz: 440 } }]) {
      assert.throws(() => createPianoEngine(Engine, 48000, options, bytes), /rebuilding.*WASM/);
    }
  }
  assert.equal(calls.length, 4); // never substitutes the original ET scale
});

test('malformed physical controls and clocks refuse before engine allocation', () => {
  const { Engine, calls } = recording();
  for (const stretchTuning of [null, [], '', {}, { concertAHz: 429 }, { concertAHz: 451 },
    { concertAHz: NaN }, { concertAHz: '440' }, { concertAHz: 440, concertPitch: 442 },
    { concertAHz: 440, bassOctaveWeights: [0,0,0,0] },
    { concertAHz: 440, bassOctaveWeights: [-1,1,0,0] },
    { concertAHz: 440, trebleOctaveWeights: [1,0,0] },
    { concertAHz: 440, trebleOctaveWeights: [Infinity,0,0,0] },
    { concertAHz: 440, trebleOctaveWeights: new Array(4) },
    { concertAHz: 440, trebleOctaveWeights: [Number.MAX_VALUE,Number.MIN_VALUE,0,0] },
    { concertAHz: 440, maxPartialHz: 0 }, { concertAHz: 440, maxPartialHz: NaN }]) {
    assert.throws(() => createPianoEngine(Engine, 48000, { stretchTuning }));
  }
  for (const scaleCsv of [null, '', '  ', 42, new Uint8Array([1])]) {
    assert.throws(() => createPianoEngine(Engine, 48000, { scaleCsv }));
  }
  for (const substeps of [0,3,8]) {
    assert.throws(() => createPianoEngine(Engine, 48000, { scaleCsv: scale, substeps }));
  }
  assert.throws(() => createPianoEngine(Engine, 96000, { scaleCsv: scale, substeps: 4 }));
  assert.equal(calls.length, 0);
});

test('a physical-owner refusal is preserved rather than retried with defaults', () => {
  const { Engine, calls } = recording();
  const refusal = new Error('key 81 requires supplied octave neighbour 69');
  Engine.withPreparation = () => { throw refusal; };
  assert.throws(() => createPianoEngine(Engine, 48000,
    { scaleCsv: scale, stretchTuning: { concertAHz: 442 } }), (error) => error === refusal);
  assert.equal(calls.length, 0);
});

test('text allowance is UTF-8 bytes and typed weights are transported numerically', () => {
  const { Engine, calls } = recording();
  const exact = '#'+ 'é'.repeat(131071) + '\n'; // exactly 256 KiB UTF-8
  assert.equal(new TextEncoder().encode(exact).length, 256 * 1024);
  createPianoEngine(Engine, 48000, { scaleCsv: exact,
    stretchTuning: { concertAHz: 440, bassOctaveWeights: new Float64Array([1,2,3,4]) } });
  assert.equal(calls[0].at(-2), exact);
  assert.deepEqual(calls[0].at(-1), new Float64Array([440,1,2,3,4,1,0,0,0,0]));
  assert.throws(() => createPianoEngine(Engine, 48000, { scaleCsv: exact + 'x' }), /256 KiB/);
  assert.equal(calls.length, 1);
});
