// Execute the actual worker selection/rebuild path with recording WASM and
// browser boundaries. These tests check delivery of physical inputs, not sound.
import test from 'node:test';
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { runInNewContext } from 'node:vm';
import { gasSourcePacket } from '../public/live/gas-sources.js';

const url = new URL('../public/live/kit/worker.js', import.meta.url);
const source = (await readFile(url, 'utf8'))
  .replace(/^import .*;\n/gm, '')
  .replaceAll('import.meta.url', JSON.stringify(url.href));

async function run(gasRecord, { copies = true, slowFirst = false, refuseFine = false,
  missingWithGas = false, refuseGas = false, scales = [3, 5], materials,
  missingWithMaterials = false, refuseMaterials = false,
  missingWithBendingLoss = false, refuseBendingLoss = false, gasSources,
  missingWithGasSources = false, refuseGasSources = false, missingWithStatedUnknown = false } = {}) {
  const calls = [], messages = [], kept = [];
  let id = 0, active = 0, clock = 0, attached = 0, correctionCalls = 0;
  function piece(name, scale, record, materials) {
    const instance = {
      id: ++id, name, scale, record, materials,
      copy() {
        if (!copies) throw new Error('old engine cannot copy');
        const copy = { ...this }; kept.push(copy); return copy;
      },
      strike() {}, radius() { return 0.2; }, free() {}, asleep() { return true; },
      render(n) { active = this.id; return new Float32Array(n); },
      info() { return JSON.stringify({ piece: name, gasRecord: record, materials }); },
      surfaceXY() { return new Float32Array([0, 0]); }, surfaceH() { return 0.01; },
    };
    return instance;
  }
  class KitPiece {
    constructor(name, scale) { calls.push(['default', name, scale]); return piece(name, scale); }
    static withGas(name, scale, record) {
      calls.push(['gas', name, scale, record]);
      if (refuseGas) throw new Error('outside admitted gas domain');
      if (refuseFine && scale < 5) throw new Error('fine grid needs another mechanical clock');
      return piece(name, scale, record);
    }
    static withMaterials(name, scale, record, ...inputs) {
      calls.push(['materials', name, scale, record, ...inputs]);
      if (refuseMaterials) throw new Error('material card has no admitted isotropic elasticity');
      if (refuseFine && scale < 5) throw new Error('fine grid needs another mechanical clock');
      return piece(name, scale, record, { importedPacks: 3 });
    }
    static withMaterialBendingLoss(name, scale, record, ...inputs) {
      calls.push(['bending', name, scale, record, ...inputs]);
      if (refuseBendingLoss) throw new Error('bending-viscosity source must cover internal Nyquist');
      if (refuseFine && scale < 5) throw new Error('fine grid needs another mechanical clock');
      return piece(name, scale, record, { importedPacks: 3, bendingLoss: { model: 'kelvin-voigt-plane-stress' } });
    }
    static withMaterialsStatedUnknownFilmLoss(name, scale, record, ...inputs) {
      calls.push(['stated-unknown', name, scale, record, ...inputs]);
      if (refuseFine && scale < 5) throw new Error('fine grid needs another mechanical clock');
      return piece(name, scale, record, { importedPacks: 3, filmLossStatedUnknown: 'receipt' });
    }
    static withMaterialsStatedUnknownFilmLossAndGasSources(name, scale, record, ...inputs) {
      calls.push(['source-stated-unknown', name, scale, record, ...inputs]);
      if (refuseGasSources) throw new Error('canonical gas source domain refused');
      if (refuseFine && scale < 5) throw new Error('fine grid needs another mechanical clock');
      return piece(name, scale, record, { importedPacks: 3, filmLossStatedUnknown: 'receipt' });
    }
    static withGasSources(name, scale, record, ...inputs) {
      calls.push(['source-gas', name, scale, record, ...inputs]);
      if (refuseGasSources) throw new Error('canonical gas source domain refused');
      if (refuseFine && scale < 5) throw new Error('fine grid needs another mechanical clock');
      return piece(name, scale, record);
    }
    static withMaterialsAndGasSources(name, scale, record, ...inputs) {
      calls.push(['source-materials', name, scale, record, ...inputs]);
      if (refuseGasSources) throw new Error('canonical gas source domain refused');
      if (refuseFine && scale < 5) throw new Error('fine grid needs another mechanical clock');
      return piece(name, scale, record, { importedPacks: 3 });
    }
    static withMaterialBendingLossAndGasSources(name, scale, record, ...inputs) {
      calls.push(['source-bending', name, scale, record, ...inputs]);
      if (refuseGasSources) throw new Error('canonical gas source domain refused');
      if (refuseFine && scale < 5) throw new Error('fine grid needs another mechanical clock');
      return piece(name, scale, record, { importedPacks: 3, bendingLoss: { model: 'kelvin-voigt-plane-stress' } });
    }
  }
  if (missingWithGas) KitPiece.withGas = undefined;
  if (missingWithMaterials) KitPiece.withMaterials = undefined;
  if (missingWithBendingLoss) KitPiece.withMaterialBendingLoss = undefined;
  if (missingWithStatedUnknown) {
    KitPiece.withMaterialsStatedUnknownFilmLoss = undefined;
    KitPiece.withMaterialsStatedUnknownFilmLossAndGasSources = undefined;
  }
  if (missingWithGasSources) {
    KitPiece.withGasSources = undefined;
    KitPiece.withMaterialsAndGasSources = undefined;
    KitPiece.withMaterialBendingLossAndGasSources = undefined;
  }
  const self = { postMessage(message) { messages.push(message); } };
  runInNewContext(source, { self, KitPiece, init: async () => {}, gasSourcePacket,
    createLane: () => { attached++; return { frame: 0 }; }, levelFix: () => { correctionCalls++; return 2; },
    Float32Array, Uint8Array, Uint32Array, URL,
    performance: { now() { clock += slowFirst && active === 1 ? 10 : 1; return clock; } },
    fetch: async () => ({ ok: true, clone() { return this; }, arrayBuffer: async () => new ArrayBuffer(0),
      json: async () => ({}) }),
  });
  await self.onmessage({ data: { type: 'init', piece: 'snare', gasRecord, gasSources, materials, port: {}, scales, minRealtime: 1.7 } });
  return { calls, messages, kept, attached, correctionCalls };
}

const materialInput = Object.freeze({
  surfacePack: new Uint8Array([1, 2, 3]), secondaryPack: new Uint8Array([4, 5]),
  strikerPack: new Uint8Array([6, 7, 8]), temperatureAxis: 'T', temperatureK: 293.15,
  strikerVolumeM3: 9e-6, strikerTipRadiusM: 0.005, impactLossSM: 0.05,
});

const sourceRecord = 'thermochemical 300 100000 mole 2 custom-alpha 0.7 custom-beta 0.3';
const sourceInput = () => ({ phase: 'frozen', waterSpecies: null, components: [
  { speciesPack: new Uint8Array([99, 1, 2, 99]).subarray(1, 3), nasa9Pack: new Uint8Array([3]), transportPack: new Uint8Array([4, 5]) },
  { speciesPack: new Uint8Array([6]), nasa9Pack: new Uint8Array([7, 8]), transportPack: new Uint8Array([9]) },
] });

test('canonical gas bytes compose with both solid routes through every grid, retry and copy fallback', async () => {
  const bendingLoss = { model: 'kelvin-voigt-plane-stress', frequencyAxis: 'frequency',
    frequencyHz: 1000, primaryAngleRad: 0.3, secondaryAngleRad: -0.7 };
  for (const materials of [undefined, materialInput, { ...materialInput, bendingLoss }]) {
    const expected = materials ? (materials.bendingLoss ? 'source-bending' : 'source-materials') : 'source-gas';
    for (const options of [{}, { copies: false, slowFirst: true }, { refuseFine: true, scales: [2, 3] }]) {
      const result = await run(sourceRecord, { ...options, materials, gasSources: sourceInput() });
      assert.equal(result.messages.at(-1).type, 'ready');
      assert.equal(result.messages.at(-1).gasRecord, sourceRecord);
      assert.ok(result.calls.length >= 2);
      for (const call of result.calls) {
        assert.equal(call[0], expected);
        assert.equal(call[3], sourceRecord);
        assert.deepEqual([...call.at(-4)], [1, 2, 3, 4, 5, 6, 7, 8, 9]);
        assert.deepEqual([...call.at(-3)], [2, 1, 2, 1, 2, 1]);
        assert.equal(call.at(-2), '');
        assert.equal(call.at(-1), 'frozen');
      }
      assert.equal(result.correctionCalls, 0, 'source atmosphere retains physical pressure gain');
      assert.equal(result.messages.at(-1).levels, 'none');
      if (materials) assert.equal(result.messages.at(-1).info.materials.importedPacks, 3);
      if (materials?.bendingLoss) assert.equal(result.messages.at(-1).info.materials.bendingLoss.model, bendingLoss.model);
    }
  }
});

test('source gas refusal or an older WASM cannot fall back to ordinary gas or discard solid selection', async () => {
  for (const materials of [undefined, materialInput, { ...materialInput,
    bendingLoss: { model: 'kelvin-voigt-plane-stress', frequencyAxis: 'f', frequencyHz: 1000,
      primaryAngleRad: 0, secondaryAngleRad: 0 } }]) {
    for (const options of [{ missingWithGasSources: true }, { refuseGasSources: true }]) {
      const result = await run(sourceRecord, { ...options, materials, gasSources: sourceInput() });
      assert.equal(result.messages.at(-1).type, 'error');
      assert.equal(result.attached, 0);
      assert.ok(result.calls.every(call => call[0].startsWith('source-')));
    }
  }
});

test('source associations and phase declarations are explicit and bounded before building grids', async () => {
  const valid = sourceInput();
  for (const gasSources of [null, {}, { ...valid, components: [] }, { ...valid, components: Array(129).fill(valid.components[0]) },
    { ...valid, waterSpecies: undefined }, { ...valid, waterSpecies: '' }, { ...valid, phase: 'equilibrium' },
    { ...valid, phase: 'liquid-water' }, { ...valid, components: [{ ...valid.components[0], nasa9Pack: new Uint8Array() }] },
    { ...valid, components: [{ ...valid.components[0], speciesPack: [1, 2] }] }]) {
    const result = await run(sourceRecord, { gasSources });
    assert.equal(result.messages.at(-1).type, 'error');
    assert.equal(result.calls.length, 0);
  }
  for (const record of [undefined, '300 100000 0', 'mixture 300 100000 mole 1 N2 1']) {
    const result = await run(record, { gasSources: valid });
    assert.equal(result.messages.at(-1).type, 'error');
    assert.equal(result.calls.length, 0);
  }
  const wet = await run(sourceRecord, { gasSources: { ...valid, phase: 'liquid-water', waterSpecies: 'custom-water' } });
  assert.equal(wet.messages.at(-1).type, 'ready');
  assert.ok(wet.calls.every(call => call.at(-2) === 'custom-water' && call.at(-1) === 'liquid-water'));
});

test('source viscosity selection and orientations reach every actual constructor path', async () => {
  const bendingLoss = { model: 'kelvin-voigt-plane-stress', frequencyAxis: 'frequency',
    frequencyHz: 1000, primaryAngleRad: 0.3, secondaryAngleRad: -0.7 };
  const materials = { ...materialInput, bendingLoss };
  for (const options of [{}, { copies: false, slowFirst: true }, { refuseFine: true, scales: [2, 3] },
    { missingWithMaterials: true }]) {
    const result = await run(undefined, { ...options, materials });
    assert.ok(result.calls.length >= 2);
    for (const call of result.calls) {
      assert.equal(call[0], 'bending');
      assert.deepEqual(call.slice(3), ['293.15 101325 0', materialInput.surfacePack,
        materialInput.secondaryPack, materialInput.strikerPack, 'T', 293.15,
        9e-6, 0.005, 0.05, 'frequency', 1000, 0.3, -0.7]);
    }
    assert.equal(result.messages.at(-1).info.materials.bendingLoss.model, bendingLoss.model);
    assert.equal(result.correctionCalls, 0);
  }
});

test('selected source loss cannot silently use old WASM, authored loss, or an unsupported source', async () => {
  const bendingLoss = { model: 'kelvin-voigt-plane-stress', frequencyAxis: 'frequency',
    frequencyHz: 1000, primaryAngleRad: 0, secondaryAngleRad: 0 };
  for (const options of [{ missingWithBendingLoss: true }, { refuseBendingLoss: true }]) {
    const result = await run(undefined, { ...options, materials: { ...materialInput, bendingLoss } });
    assert.equal(result.attached, 0);
    assert.equal(result.messages.at(-1).type, 'error');
    assert.match(result.messages.at(-1).message, /source bending-loss support|internal Nyquist/);
    assert.ok(result.calls.every(call => call[0] === 'bending'));
  }
  for (const loss of [null, {}, { ...bendingLoss, model: 'infer-from-density' },
    { ...bendingLoss, frequencyHz: NaN }, { ...bendingLoss, secondaryAngleRad: Infinity }]) {
    const result = await run(undefined, { materials: { ...materialInput, bendingLoss: loss } });
    assert.equal(result.attached, 0);
    assert.equal(result.calls.length, 0);
    assert.equal(result.messages.at(-1).type, 'error');
  }
});

test('a stated-unknown head film loss reaches only its own constructors, with or without source gas', async () => {
  const materials = { ...materialInput, filmLoss: 'stated-unknown' };
  for (const [record, gas, expected] of [[undefined, undefined, 'stated-unknown'],
    [sourceRecord, true, 'source-stated-unknown']]) {
    for (const options of [{}, { copies: false, slowFirst: true }, { refuseFine: true, scales: [2, 3] },
      { missingWithMaterials: true }]) {
      const result = await run(record, { ...options, materials, gasSources: gas && sourceInput() });
      assert.equal(result.messages.at(-1).type, 'ready');
      assert.ok(result.calls.length >= 2);
      for (const call of result.calls) {
        assert.equal(call[0], expected);
        assert.deepEqual(call.slice(4, 12), [materialInput.surfacePack, materialInput.secondaryPack,
          materialInput.strikerPack, 'T', 293.15, 9e-6, 0.005, 0.05]);
      }
      assert.equal(result.messages.at(-1).info.materials.filmLossStatedUnknown, 'receipt');
    }
  }
});

test('a film-loss declaration never reaches a preset-decay constructor or an older WASM', async () => {
  const declared = { ...materialInput, filmLoss: 'stated-unknown' };
  for (const gasSources of [undefined, sourceInput()]) {
    const result = await run(gasSources ? sourceRecord : undefined,
      { missingWithStatedUnknown: true, materials: declared, gasSources });
    assert.equal(result.attached, 0);
    assert.equal(result.messages.at(-1).type, 'error');
    assert.match(result.messages.at(-1).message, /stated-unknown film-loss route/);
    assert.equal(result.calls.length, 0);
  }
  const bendingLoss = { model: 'kelvin-voigt-plane-stress', frequencyAxis: 'frequency',
    frequencyHz: 1000, primaryAngleRad: 0, secondaryAngleRad: 0 };
  for (const materials of [{ ...materialInput, filmLoss: 'zero' }, { ...materialInput, filmLoss: null },
    { ...declared, bendingLoss }]) {
    const result = await run(undefined, { materials });
    assert.equal(result.attached, 0);
    assert.equal(result.calls.length, 0);
    assert.equal(result.messages.at(-1).type, 'error');
  }
});

test('material-card bytes and impactor geometry survive every grid, retry and clone fallback', async () => {
  const gas = 'mixture 300 90000 mole 2 N2 0.8 CO2 0.2';
  for (const options of [{}, { copies: false, slowFirst: true }, { refuseFine: true, scales: [2, 3] }]) {
    const result = await run(gas, { ...options, materials: materialInput });
    assert.ok(result.calls.length >= 2);
    for (const call of result.calls) {
      assert.equal(call[0], 'materials');
      assert.deepEqual(call.slice(3), [gas, materialInput.surfacePack, materialInput.secondaryPack,
        materialInput.strikerPack, 'T', 293.15, 9e-6, 0.005, 0.05]);
    }
    assert.equal(result.messages.at(-1).info.materials.importedPacks, 3);
    assert.equal(result.messages.at(-1).levelFix, 1);
    assert.equal(result.messages.at(-1).levels, 'none');
    assert.equal(result.correctionCalls, 0, 'default-material level calibration must not normalize a material swap');
  }
});

test('material import defaults only the declared reusable secondary and standard atmosphere', async () => {
  const { secondaryPack, ...materials } = materialInput;
  const result = await run(undefined, { materials });
  for (const call of result.calls) {
    assert.equal(call[3], '293.15 101325 0');
    assert.deepEqual(call[5], new Uint8Array());
  }
  assert.equal(result.attached, 1);
});

test('unsupported or refused material data never falls back to preset bronze or Mylar', async () => {
  for (const options of [{ missingWithMaterials: true }, { refuseMaterials: true }]) {
    const result = await run('293.15 101325 0', { ...options, materials: materialInput });
    assert.equal(result.attached, 0);
    assert.equal(result.messages.at(-1).type, 'error');
    assert.match(result.messages.at(-1).message, /material-card support|no admitted isotropic elasticity/);
    assert.ok(result.calls.every(call => call[0] === 'materials'));
  }
  for (const materials of [null, {}, { ...materialInput, surfacePack: [] },
    { ...materialInput, strikerPack: new Uint8Array() }, { ...materialInput, temperatureK: NaN }]) {
    const result = await run(undefined, { materials });
    assert.equal(result.attached, 0);
    assert.equal(result.calls.length, 0);
    assert.equal(result.messages.at(-1).type, 'error');
  }
});

test('the chosen clone and every automatic grid trial retain the exact mixture', async () => {
  const record = 'mixture 300 80000 mole 3 N2 0.7 O2 0.2 CO2 0.1';
  const result = await run(record);
  assert.deepEqual(result.calls, [['gas', 'snare', 5, record], ['gas', 'snare', 3, record]]);
  assert.ok(result.kept.every(copy => copy.record === record));
  const ready = result.messages.find(m => m.type === 'ready');
  assert.equal(ready.gasRecord, record);
  assert.equal(ready.info.gasRecord, record);
  assert.equal(ready.scale, 3);
  assert.equal(result.attached, 1);
});

test('warmup retries and old-WASM clone fallbacks use the same gas constructor', async () => {
  const record = '293.15 65000 0.5';
  const result = await run(record, { copies: false, slowFirst: true });
  assert.deepEqual(result.calls, [['gas', 'snare', 5, record], ['gas', 'snare', 5, record],
    ['gas', 'snare', 3, record], ['gas', 'snare', 3, record]]);
  assert.equal(result.messages.at(-1).gasRecord, record);
});

test('fallback after every requested grid refuses still retains mass-fraction composition', async () => {
  const record = 'mixture 295 100000 mass 2 N2 0.8 CO2 0.2';
  const result = await run(record, { scales: [2, 3], copies: false, refuseFine: true });
  assert.deepEqual(result.calls, [['gas', 'snare', 3, record], ['gas', 'snare', 2, record],
    ['gas', 'snare', 5, record], ['gas', 'snare', 5, record]]);
  assert.equal(result.messages.at(-1).gasRecord, record);
});

test('missing atmosphere support or rejected material data never silently creates default air', async () => {
  for (const options of [{ missingWithGas: true }, { refuseGas: true }]) {
    const result = await run('mixture 300 80000 mole 1 CO2 1', options);
    assert.equal(result.attached, 0);
    assert.equal(result.messages.at(-1).type, 'error');
    assert.match(result.messages.at(-1).message, /rebuilding with atmosphere support|outside admitted gas domain/);
    assert.ok(result.calls.every(call => call[0] === 'gas'));
  }
  for (const record of ['', null, 7]) {
    const result = await run(record);
    assert.equal(result.attached, 0);
    assert.match(result.messages.at(-1).message, /nonempty atmosphere record/);
    assert.equal(result.calls.length, 0);
  }
});

test('a restarted worker adopts its new atmosphere while omitted records keep the standard path', async () => {
  const air = await run('293.15 101325 0');
  const mixture = await run('mixture 303.15 95000 mole 1 Ar 1');
  assert.equal(air.messages.at(-1).info.gasRecord, '293.15 101325 0');
  assert.equal(mixture.messages.at(-1).info.gasRecord, 'mixture 303.15 95000 mole 1 Ar 1');
  const legacy = await run(undefined);
  assert.ok(legacy.calls.every(call => call[0] === 'default'));
  assert.equal(legacy.messages.at(-1).type, 'ready');
});
