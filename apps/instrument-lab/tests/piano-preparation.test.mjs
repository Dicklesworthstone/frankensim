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
    static withAtmosphere(...args) { calls.push(['atmosphere', ...args]); return { kind: 'atmosphere' }; }
    static withStringMaterials(...args) { calls.push(['materials', ...args]); return { kind: 'materials' }; }
  }
  return { Engine, calls };
}
const scale = '# supplied SI scale text is passed to Rust, never interpreted here\n';
const geometry = 'frankensim-piano-materials-v2\ntension,target-pitch\nstrength,ultimate_tensile_strength\n'
  + 'key,layer,outer_radius_m,wire_diameter_m,pitch_m,material_index\n69,core,0.0004,,,0\n69,core,0.0006,,,1\n';
const materialInput = () => ({ packs: [new Uint8Array([1, 2, 3]), new Uint8Array([4, 5]).buffer],
  geometry, temperatureAxis: 'T', temperatureK: 293.15, viscosityProperty: 'solid_viscosity',
  frequencyAxis: 'frequency', frequencyHz: 1000 });

test('canonical material framing preserves byte views and identities through every board and clock', () => {
  const { Engine, calls } = recording();
  const original = new Uint8Array([90, 1, 2, 3, 91]);
  const packs = Object.freeze([original.subarray(1, 4), new Uint8Array([4, 5]).buffer]);
  const stringMaterials = Object.freeze({ ...materialInput(), packs });
  const board = new Uint8Array([7, 8, 9]);
  const gas = 'mixture 300 90000 mole 2 N2 0.75 O2 0.25';
  for (const substeps of [1, 2, 4]) {
    for (const bytes of [null, board]) {
      for (const boardCoordinates of [0, 16, 64]) {
        const options = Object.freeze({ stringMaterials, scaleCsv: scale,
          stretchTuning: Object.freeze({ concertAHz: 442, maxPartialHz: 18000 }), gas, substeps,
          boardCoordinates, courseModeBudget: 180, maxDynamicModes: 2100, fullScalePa: 7 });
        assert.equal(createPianoEngine(Engine, 48000, options, bytes).kind, 'materials');
        assert.deepEqual(calls.at(-1), ['materials', 48000, bytes ?? new Uint8Array(),
          bytes || boardCoordinates > 16 ? 0 : boardCoordinates, 180, 2100, 7, substeps, scale,
          new Float64Array([442,1,0,0,0,1,0,0,0,18000]), gas,
          new Uint8Array([1,2,3,4,5]), new Uint32Array([3,2]), geometry,
          'T', 293.15, 'solid_viscosity', 'frequency', 1000]);
        if (bytes) assert.equal(calls.at(-1)[2], board);
      }
    }
  }
  calls.at(-1)[11][0] = 255;
  assert.deepEqual(original, new Uint8Array([90,1,2,3,91]));
  assert.deepEqual(new Uint8Array(packs[1]), new Uint8Array([4,5]));
});

test('material-only construction keeps the embedded scale and an explicit standard atmosphere', () => {
  const { Engine, calls } = recording();
  Engine.withPreparation = Engine.withAtmosphere = undefined;
  createPianoEngine(Engine, 48000, { stringMaterials: materialInput() });
  assert.equal(calls[0][8], '');
  assert.deepEqual(calls[0][9], new Float64Array());
  assert.equal(calls[0][10], '293.15 101325 0');
});

test('only the Rust owner interprets geometry syntax, source bytes and section coverage', () => {
  const { Engine, calls } = recording();
  for (const text of [`# original source receipt\n\n ${geometry.replaceAll('\n', '\r\n')}\n`,
    'nonempty text whose geometry validity belongs to Rust']) {
    createPianoEngine(Engine, 48000, { stringMaterials: { ...materialInput(), geometry: text } });
    assert.equal(calls.at(-1)[13], text);
  }
});

test('requested materials never fall back on old WASM or swallow a canonical-owner refusal', () => {
  const { Engine, calls } = recording();
  const options = { stringMaterials: materialInput(), scaleCsv: scale, gas: '293.15 50000 0' };
  Engine.withStringMaterials = undefined;
  assert.throws(() => createPianoEngine(Engine, 48000, options), /rebuilding.*WASM.*withStringMaterials/);
  const refusal = new Error('selected material geometry does not cover every scale key');
  Engine.withStringMaterials = () => { throw refusal; };
  assert.throws(() => createPianoEngine(Engine, 48000, options), (error) => error === refusal);
  assert.equal(calls.length, 0);
  // An old engine still supports every previously exposed route when materials were not requested.
  createPianoEngine(Engine, 48000);
  createPianoEngine(Engine, 48000, { scaleCsv: scale });
  createPianoEngine(Engine, 48000, { gas: options.gas });
  assert.deepEqual(calls.map((call) => call[0]), ['legacy', 'prepared', 'atmosphere']);
});

test('malformed or incomplete material controls refuse before any engine allocation', () => {
  const { Engine, calls } = recording();
  const valid = materialInput();
  const invalid = [null, [], '', {}, { ...valid, temperature: 300 }];
  for (const key of Object.keys(valid)) {
    invalid.push(Object.fromEntries(Object.entries(valid).filter(([name]) => name !== key)));
  }
  for (const packs of [null, [], new Uint8Array([1]), [new Uint8Array()], [new ArrayBuffer(0)],
    [new Uint16Array([1])], [new DataView(new ArrayBuffer(1))], [[1,2]], new Array(2),
    Array.from({ length: 129 }, () => new Uint8Array([1]))]) invalid.push({ ...valid, packs });
  for (const geometry of [null, '', ' ', 42]) invalid.push({ ...valid, geometry });
  for (const key of ['temperatureAxis', 'frequencyAxis', 'viscosityProperty']) {
    for (const value of ['', ' ', ' T', 'T ', null, 123, 'é'.repeat(65)]) invalid.push({ ...valid, [key]: value });
  }
  invalid.push({ ...valid, frequencyAxis: valid.temperatureAxis });
  for (const key of ['temperatureK', 'frequencyHz']) {
    for (const value of [0, -1, Infinity, NaN, '300', null]) invalid.push({ ...valid, [key]: value });
  }
  for (const stringMaterials of invalid) {
    assert.throws(() => createPianoEngine(Engine, 48000, { stringMaterials }));
  }
  for (const gas of ['', null, 'é'.repeat(4097)]) {
    assert.throws(() => createPianoEngine(Engine, 48000, { stringMaterials: valid, gas }));
  }
  assert.equal(calls.length, 0);
});

test('material count and UTF-8/pack byte budgets admit exact limits and refuse overflow', () => {
  const { Engine, calls } = recording();
  const valid = materialInput();
  const prefix = geometry + '#';
  const remaining = 256 * 1024 - new TextEncoder().encode(prefix).length;
  const exactGeometry = prefix + 'é'.repeat(Math.floor(remaining / 2)) + (remaining % 2 ? 'x' : '');
  const largePack = new Uint8Array(16 * 1024 * 1024);
  largePack[0] = 8; largePack[largePack.length - 1] = 9;
  createPianoEngine(Engine, 48000, { stringMaterials: { ...valid, geometry: exactGeometry,
    packs: [largePack], temperatureAxis: 'é'.repeat(64) } });
  assert.equal(calls[0][11].length, largePack.length);
  assert.equal(calls[0][11][0], 8); assert.equal(calls[0][11].at(-1), 9);
  assert.deepEqual(calls[0][12], new Uint32Array([largePack.length]));
  assert.equal(calls[0][13], exactGeometry);
  assert.throws(() => createPianoEngine(Engine, 48000, { stringMaterials: { ...valid,
    packs: [largePack, new Uint8Array([1])] } }), /16 MiB/);
  assert.throws(() => createPianoEngine(Engine, 48000, { stringMaterials: { ...valid,
    geometry: exactGeometry + 'x' } }), /256 KiB/);
  createPianoEngine(Engine, 48000, { stringMaterials: { ...valid,
    packs: Array.from({ length: 128 }, () => new Uint8Array([1])) } });
  assert.equal(calls[1][12].length, 128);
  assert.equal(calls.length, 2);
});

test('one atmosphere reaches both board paths, physical preparation and every audio clock', () => {
  const { Engine, calls } = recording();
  const gas = 'mixture 293.15 100000 mole 3 N2 0.7 O2 0.2 CO2 0.1';
  const board = new Uint8Array([4, 5, 6]);
  for (const substeps of [1, 2, 4]) {
    for (const bytes of [null, board]) {
      const options = Object.freeze({ gas, substeps, boardCoordinates: 16, scaleCsv: scale,
        stretchTuning: { concertAHz: 442 } });
      assert.equal(createPianoEngine(Engine, 48000, options, bytes).kind, 'atmosphere');
      assert.deepEqual(calls.at(-1), ['atmosphere', 48000, bytes ?? new Uint8Array(),
        bytes ? 0 : 16, 0, 0, 0, substeps, scale,
        new Float64Array([442, 1,0,0,0, 1,0,0,0, 0]), gas]);
    }
  }
  createPianoEngine(Engine, 48000, { gas });
  assert.equal(calls.at(-1).at(-3), '');
  assert.deepEqual(calls.at(-1).at(-2), new Float64Array());
  assert.equal(calls.at(-1).at(-1), gas);
});

test('requested atmosphere never falls back to room air on old WASM or a physical refusal', () => {
  const { Engine, calls } = recording();
  const gas = '293.15 50000 0';
  Engine.withAtmosphere = undefined;
  assert.throws(() => createPianoEngine(Engine, 48000, { gas }), /rebuilding.*WASM/);
  const refusal = new Error('material gas outside its domain');
  Engine.withAtmosphere = () => { throw refusal; };
  assert.throws(() => createPianoEngine(Engine, 48000, { gas }), (error) => error === refusal);
  assert.equal(calls.length, 0);
});

test('malformed atmosphere input is refused before allocating an engine', () => {
  const { Engine, calls } = recording();
  for (const gas of [null, '', ' ', 123, {}, 'é'.repeat(4097)]) {
    assert.throws(() => createPianoEngine(Engine, 48000, { gas }));
  }
  assert.equal(calls.length, 0);
});

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
