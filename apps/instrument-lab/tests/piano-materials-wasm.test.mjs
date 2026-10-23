// The production constructor selector driving the REAL piano WASM with the committed string constructions
// and the compiled seed packs the engine-piano tests use (bead pexfp.369, gap 5). public/live/piano/engine.js
// and engine_bg.wasm are build products of engine-piano/build-wasm.sh (not source), so these tests run only
// where that build exists and report a skip, not a pass, elsewhere. piano-preparation.test.mjs covers the
// selector everywhere with a recording engine.
import test from 'node:test';
import assert from 'node:assert/strict';
import { existsSync, readFileSync } from 'node:fs';
import { createPianoEngine } from '../public/live/piano/engine-factory.js';

const at = (path) => new URL(path, import.meta.url);
const wasm = at('../public/live/piano/engine_bg.wasm');
const skip = existsSync(wasm) && existsSync(at('../public/live/piano/engine.js')) ? false
  : 'the piano WASM is not built here (engine-piano/build-wasm.sh installs it)';
const pack = (name) => new Uint8Array(readFileSync(at(`../engine-piano/tests/fixtures/packs/${name}.fsmatpk`)));
const construction = (name) => readFileSync(at(`../engine-piano/data/constructions/${name}.csv`), 'utf8');
// The operating point the native tests use: room temperature, the steel pack's in-band excess property.
const point = { temperatureAxis: 'T', temperatureK: 293.15, viscosityProperty: 'bending_loss_tangent',
  frequencyAxis: 'frequency', frequencyHz: 1000 };
const steel = () => [pack('steel-music-wire-string-heptonstall-ligo-p1400229'), pack('copper-pure-nbs-c447-string-wrap'),
  pack('steel-music-wire-min-strength-qq-w-470b')];

async function engine() {
  const module = await import(at('../public/live/piano/engine.js').href);
  module.initSync({ module: readFileSync(wasm) });
  return module.PianoEngine;
}
/** The message of a refusal the WASM throws. */
function refusal(run) {
  try { run(); } catch (error) { return String(error?.message ?? error); }
  assert.fail('expected the Rust owner to refuse');
}

test('real WASM: steel strings from packs play; gold at concert pitch refuses all 88 keys with each highest pitch',
  { skip }, async () => {
    const PianoEngine = await engine();
    const piano = createPianoEngine(PianoEngine, 48000,
      { stringMaterials: { ...point, packs: steel(), geometry: construction('steel-copper-wound-rt0425') } });
    for (const key of [21, 60, 96]) assert.equal(piano.noteOn(key, 2.0), true, `key ${key}`);
    const pcm = piano.render(4800);
    assert.equal(pcm.length, 4800);
    assert.ok(pcm.every(Number.isFinite) && pcm.some((x) => x !== 0));
    assert.equal(piano.isPoisoned(), false);

    const gold = pack('gold-pure-nbs-c447');
    const message = refusal(() => createPianoEngine(PianoEngine, 48000,
      { stringMaterials: { ...point, packs: [gold, gold], geometry: construction('gold-concert-pitch-rt0425') } }));
    const json = message.slice(message.indexOf('{"refusal"'));
    assert.match(json, /^\{"refusal":"overstressed-strings","limitProperty":"ultimate_tensile_strength_cold_rolled_wire"/);
    assert.match(json, /"limitCaveats":\["strength extrapolated; 'ultimate_tensile_strength_cold_rolled_wire' from gold-pure-nbs-c447 holds only for 1\.2700 mm wire/);
    assert.match(json, /"geometryCaveat":"core and wrap diameters of keys 21-40 are back-solved/);
    assert.match(json, /"keys":88,/);
    const keys = [...json.matchAll(/\{"key":(\d+),[^}]*"fMaxHz":([0-9.e+-]+)/g)];
    assert.deepEqual(keys.map((m) => Number(m[1])), Array.from({ length: 88 }, (_, i) => 21 + i));
    assert.ok(keys.every((m) => Number(m[2]) > 0), 'every key names its highest attainable pitch');
  });

test('real WASM: the retired v1 format, a diameter-free limit setting yield-fraction tension and an undeclared one refuse',
  { skip }, async () => {
    const PianoEngine = await engine();
    const v1 = 'frankensim-piano-materials-v1\nkey,outer_radius_m,material_index\n69,0.0004,0\n';
    assert.match(refusal(() => createPianoEngine(PianoEngine, 48000,
      { stringMaterials: { ...point, packs: steel(), geometry: v1 } })), /frankensim-piano-materials-v1 declares no tension policy/);
    const silver = [pack('silver-pure-nbs-c447'), pack('silver-pure-nbs-c447'), pack('silver-drawn-hard-uts-nbs-c101')];
    const concert = construction('silver-concert-pitch-rt0425');
    const half = concert.replace('\ntension,target-pitch\n', '\ntension,yield-fraction,0.5\n');
    assert.notEqual(half, concert);
    assert.match(refusal(() => createPianoEngine(PianoEngine, 48000,
      { stringMaterials: { ...point, packs: silver, geometry: half } })),
      /yield-fraction.*silver-drawn-hard-uts-nbs-c101.*nothing back-solves a yield/s);
    const unflagged = concert.replace(',diameter-free,', ',');
    assert.notEqual(unflagged, concert);
    assert.match(refusal(() => createPianoEngine(PianoEngine, 48000,
      { stringMaterials: { ...point, packs: silver, geometry: unflagged } })), /state no wire diameter.*diameter-free/s);
  });
