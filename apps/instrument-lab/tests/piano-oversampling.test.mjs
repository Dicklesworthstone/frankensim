// Real factory/worker selection tests. The recording Engine is NOT piano physics.
import test from 'node:test';
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { createPianoEngine } from '../public/live/piano/engine-factory.js';
import { chooseBoard, chooseQuality, probeRealtime, probeSustained } from '../public/live/piano/boards.js';

// Host-boundary recording only: no string dynamics or acoustic timing oracle.
function probeKeyboard(keys, report = true) {
  const events = [], absent = [];
  let frame = 0;
  const engine = {
    keyGeometry: report ? () => Float64Array.from(keys.flatMap(k => [k, 1, 0.1, 1, 440])) : undefined,
    noteOn(key, speed) { events.push(['on', frame, key, speed]); if (!keys.includes(key)) absent.push(key); return keys.includes(key); },
    noteOff(key) { events.push(['off', frame, key]); if (!keys.includes(key)) absent.push(key); },
    render(n) { frame += n; },
  };
  return { engine, events, absent, frames: () => frame };
}

test('the built-in keyboard preserves both existing startup probe performances', () => {
  const keys = Array.from({ length: 88 }, (_, i) => 21 + i);
  for (const [probe, frames] of [[probeRealtime, 6656], [probeSustained, 16384]]) {
    const legacy = probeKeyboard(keys, false), reported = probeKeyboard(keys);
    probe(legacy.engine); probe(reported.engine);
    assert.deepEqual(reported.events, legacy.events, 'reported full geometry changes no note, speed or event frame');
    assert.equal(reported.frames(), frames);
    assert.deepEqual(reported.absent, []);
  }
});

test('startup probes exercise actual imported keys instead of timing absent-note refusals', () => {
  for (const keys of [[60], [86, 90, 97], [21, 28, 33, 35, 42, 46, 53, 58, 60, 68, 75, 84, 99, 108]]) {
    for (const [probe, frames] of [[probeRealtime, 6656], [probeSustained, 16384]]) {
      const f = probeKeyboard(keys);
      probe(f.engine);
      assert.deepEqual(f.absent, [], 'every strike and release must name a physical course in this engine');
      const opening = f.events.filter(e => e[0] === 'on').slice(0, Math.min(12, keys.length));
      assert.equal(new Set(opening.map(e => e[2])).size, Math.min(12, keys.length));
      assert.equal(f.frames(), frames, 'the measured audio span is unchanged');
    }
  }
});

test('empty reported key geometry refuses before timing an idle engine', () => {
  for (const probe of [probeRealtime, probeSustained]) {
    const f = probeKeyboard([]);
    assert.throws(() => probe(f.engine), /key geometry/);
    assert.deepEqual(f.events, []);
    assert.equal(f.frames(), 0);
  }
});

for (const [name, choose] of [['keyboard', chooseBoard], ['scheduled', chooseQuality]]) {
  test(`${name} timing refusal frees the probe and its unplayed copy`, async () => {
    const engines = [], refusal = new Error('probe refused');
    const make = () => {
      const engine = { freed: 0, copy: () => make(),
        info: () => JSON.stringify({ maxDynamicModes: 2400 }),
        free() { this.freed++; } };
      engines.push(engine);
      return engine;
    };
    await assert.rejects(choose({ make, measure: () => { throw refusal; }, target: 1.3, sizes: [] }), e => e === refusal);
    assert.equal(engines.length, 2);
    assert.ok(engines.every(e => e.freed === 1), 'neither the failed probe nor its pristine copy can play');
  });
}

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

test('startup rejects unsigned mode-budget coercion before constructing any engine', () => {
  for (const field of ['courseModeBudget', 'maxDynamicModes']) {
    for (const value of [-1, 1.5, NaN, Infinity, 2 ** 32, 2 ** 32 + 1]) {
      for (const substeps of [1, 2, 4]) {
        const { Engine, calls } = recordingEngine();
        assert.throws(() => createPianoEngine(Engine, 48000, { [field]: value, substeps }), /nonnegative 32-bit integer/);
        assert.deepEqual(calls, [], 'refusal must precede WASM allocation or budget conversion');
      }
    }
  }
  const { Engine, calls } = recordingEngine();
  createPianoEngine(Engine, 48000, { courseModeBudget: 0, maxDynamicModes: 0 });
  assert.deepEqual(calls, [['new', 48000, 0, 0, 0, 0]], 'explicit zero retains startup defaults');
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
