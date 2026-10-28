// Source-edit and original-packet regressions; native physics owns gas evaluation.
import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import {nativeGasRegions, editNativeGasRegion} from '../src/native-gas-records.mjs';
import {serializeNativeExperiment, validateInput} from '../src/native-experiment.mjs';

const fixture = name => readFile(new URL(`../../../crates/fs-couple/examples/${name}`, import.meta.url), 'utf8');

test('one data editor reaches canonical flute, reed, chimney and string gases', async () => {
  const cases = [
    ['flue-dg-baffled.performance', ['gas', 'dg_exterior']],
    ['flue-dg-molecular.performance', ['gas', 'dg_exterior']],
    ['resolved-dg.performance', ['gas', 'dg_exterior']],
    ['reed-two-medium-tone-hole.performance', ['ambient', 'exterior', 'cylinder']],
    ['material_modal_valve/paired-dg-chimney.performance', ['ambient', 'dg_chimney_gas', 'dg_exterior']],
    ['material_strings/rayleigh-observed.material-strings', ['gas']],
  ];
  for (const [name, keys] of cases) {
    const source = await fixture(name), regions = nativeGasRegions(source);
    assert.deepEqual(regions.map(region => region.key), keys);
    for (const region of regions) assert.equal(editNativeGasRegion(source, region, region.values), source);
  }
});

test('adding helium changes only the selected exterior data and keeps clocks and excitation exact', async () => {
  const source = await fixture('flue-dg-baffled.performance');
  const exterior = nativeGasRegions(source).find(region => region.key === 'dg_exterior');
  const changed = editNativeGasRegion(source, exterior, {...exterior.values, temperature:'310',
    components:[{species:'N2', fraction:'0.6'}, {species:'CO2', fraction:'0.2'}, {species:'He', fraction:'0.2'}]});
  const expected = source.replace('dg_exterior mixture 293.15 101325 mole 2 N2 0.8 CO2 0.2',
    'dg_exterior mixture 310 101325 mole 3 N2 0.6 CO2 0.2 He 0.2');
  assert.equal(changed, expected);
  const selected = nativeGasRegions(changed).find(region => region.key === 'dg_exterior');
  const withoutCO2 = editNativeGasRegion(changed, selected, {...selected.values,
    components:[{species:'N2', fraction:'0.6'}, {species:'He', fraction:'0.4'}]});
  assert.match(withoutCO2, /^dg_exterior mixture 310 101325 mole 2 N2 0.6 He 0.4$/m);
  const packet = validateInput(new TextEncoder().encode(serializeNativeExperiment(changed, [], {renderer:'wind', decimate:true})), true);
  assert.equal(new TextDecoder().decode(packet.source), expected);
  assert.equal(packet.decimate, true);
});

test('common absolute pressure reaches explicit regions without changing inherited pressure or aliases', () => {
  const source = 'frankensim-reed-performance-v1\r\n'
    + 'ambient\t293.15\t101325\t0.1\r\n'
    + 'exterior same\r\n'
    + 'cylinder 0.01 0.2 gas 300 0.4\r\n'
    + 'cone 0.01 0.02 0.3 gas mixture 310 101325 mass 2 N2 0.7 O2 0.3\r\n'
    + 'duct_section 0 1 0.3 0.01 0 gas thermochemical 290 101325 mole 1 oxygen-source 1 loss viscothermal\r\n';
  const regions = nativeGasRegions(source);
  assert.equal(regions.length, 4);
  assert.equal(regions[1].pressureSpan, null);
  const changed = editNativeGasRegion(source, regions[0], {...regions[0].values, pressure:'80000'}, true);
  assert.equal(changed, source.replace('ambient\t293.15\t101325\t0.1', 'ambient\t293.15 80000 0.1').replaceAll(' 101325 ', ' 80000 '));
  assert.match(changed, /exterior same\r\n/);
  assert.match(changed, /cylinder 0.01 0.2 gas 300 0.4\r\n/);
  assert.match(changed, /oxygen-source 1 loss viscothermal\r\n/);
});

test('source-backed species edits keep original NASA-9 and transport associations in a saved experiment', () => {
  const source = 'frankensim-material-strings-v3\n'
    + 'gas thermochemical 298.15 101325 mass 2 source-O2 0.3 source-N2 0.7\n'
    + 'observer_gas stokes-kirchhoff 256 1000000 50000000\n';
  const components = ['source-O2', 'source-N2', 'source-He'].map(id => ({
    speciesPath:`gas/${id}.species`, nasa9Path:`gas/${id}.nasa9`, transportPath:`gas/${id}.transport`}));
  const options = {renderer:'strings', gasSources:{phase:'frozen', waterSpecies:null, components}};
  const files = components.flatMap(component => Object.values(component)).map((path, i) =>
    ({path, bytes:new Uint8Array([0, 255, i, 13, 10])}));
  const region = nativeGasRegions(source)[0];
  const changed = editNativeGasRegion(source, region, {...region.values, pressure:'90000',
    components:[{species:'source-O2', fraction:'0.2'}, {species:'source-N2', fraction:'0.6'}, {species:'source-He', fraction:'0.2'}]});
  const packet = validateInput(new TextEncoder().encode(serializeNativeExperiment(changed, files, options)), true);
  assert.deepEqual(packet.attachments, files);
  assert.deepEqual(packet.gasSources, options.gasSources);
  assert.match(new TextDecoder().decode(packet.source), /gas thermochemical 298.15 90000 mass 3 source-O2 0.2 source-N2 0.6 source-He 0.2/);
});

test('unchanged records preserve original whitespace and unknown source forms stay outside the editor', () => {
  const source = '\ufefffrankensim-flue-performance-v1\r\n# µ original atmosphere\r\n'
    + 'gas\tmixture\t293.15  101325 mole 2 N2 0.8 O2 0.2\r\n'
    + 'exterior same\r\ngas pack original-source\r\nobserver_gas stokes-kirchhoff 1 2 3\r\n'
    + 'gas mixture 300 101325 mole 3 N2 0.8 O2 0.2\r\n'
    + 'track\tgas\tbow\t0\r\n';
  const regions = nativeGasRegions(source);
  assert.equal(regions.length, 1);
  assert.equal(regions[0].line, 3);
  assert.equal(editNativeGasRegion(source, regions[0], regions[0].values), source);
  assert.throws(() => editNativeGasRegion('# changed\n' + source, regions[0], regions[0].values), /performance changed/);
});

test('invalid fractions are refused without normalization and changing the gas model is explicit', () => {
  const source = 'gas 293.15 101325 0\n';
  const region = nativeGasRegions(source)[0];
  const mixture = {kind:'mixture', temperature:'293.15', pressure:'101325', basis:'mole',
    components:[{species:'N2', fraction:'0.8'}, {species:'He', fraction:'0.2'}]};
  assert.equal(editNativeGasRegion(source, region, mixture), 'gas mixture 293.15 101325 mole 2 N2 0.8 He 0.2\n');
  for (const components of [[{species:'N2', fraction:'0.8'}], [{species:'N2', fraction:'1.1'}],
    [{species:'N2', fraction:'0.5'}, {species:'N2', fraction:'0.5'}], [{species:'N2\ngas', fraction:'1'}]]) {
    assert.throws(() => editNativeGasRegion(source, region, {...mixture, components}));
  }
  for (const pressure of ['0', '-10', 'Infinity', 'NaN', '0xff']) {
    assert.throws(() => editNativeGasRegion(source, region, {...mixture, pressure}), /pressure/);
  }
  assert.equal(source, 'gas 293.15 101325 0\n');
});

test('resolved dry-air temperature and common pressure reach the original primary and exterior records', async () => {
  const source = await fixture('resolved-dg.performance');
  const [primary, exterior] = nativeGasRegions(source);
  assert.equal(primary.values.kind, 'dry-air');
  assert.deepEqual(primary.supportedKinds, ['dry-air', 'thermochemical']);
  assert.deepEqual(exterior.supportedKinds, ['moist-air', 'mixture', 'thermochemical']);
  const changed = editNativeGasRegion(source, primary,
    {...primary.values, temperature:'315', pressure:'90000'}, true);
  assert.equal(changed, source.replace('gas dry-air 300 100000', 'gas dry-air 315 90000')
    .replace('dg_exterior 300 100000 0', 'dg_exterior 300 90000 0'));
  const exteriorChanged = editNativeGasRegion(source, exterior,
    {...exterior.values, temperature:'305', pressure:'85000'}, true);
  assert.equal(exteriorChanged, source.replace('gas dry-air 300 100000', 'gas dry-air 300 85000')
    .replace('dg_exterior 300 100000 0', 'dg_exterior 305 85000 0'));
});

test('resolved primary gas refuses snapshot models while exterior models retain their original choices', async () => {
  const source = await fixture('resolved-dg.performance');
  const [primary, exterior] = nativeGasRegions(source);
  for (const values of [
    {kind:'moist-air', temperature:'300', pressure:'100000', humidity:'0', components:[]},
    {kind:'mixture', temperature:'300', pressure:'100000', basis:'mole', components:[{species:'N2', fraction:'1'}]},
  ]) {
    assert.throws(() => editNativeGasRegion(source, primary, values), /supports dry-air or thermochemical/);
    assert.notEqual(editNativeGasRegion(source, exterior, {...values, temperature:'310'}), source);
  }
  assert.throws(() => editNativeGasRegion(source, exterior, primary.values), /supports moist-air or mixture or thermochemical/);
  const changed = editNativeGasRegion(source, exterior, {kind:'thermochemical', temperature:'300',
    pressure:'100000', basis:'mass', components:[{species:'original-N2', fraction:'1'}]});
  assert.equal(changed, source.replace('dg_exterior 300 100000 0',
    'dg_exterior thermochemical 300 100000 mass 1 original-N2 1'));
});

test('resolved gas model changes require the user to author or remove energy-solve controls', async () => {
  const source = (await fixture('resolved-dg.performance')).replaceAll('\n', '\r\n');
  const primary = nativeGasRegions(source)[0];
  assert.equal(editNativeGasRegion(source, primary, primary.values), source);
  const thermochemical = {kind:'thermochemical', temperature:'300', pressure:'100000',
    basis:'mass', components:[{species:'original-N2', fraction:'1'}]};
  assert.throws(() => editNativeGasRegion(source, primary, thermochemical), /author an energy_solve.*performance editor/);
  const authored = source.replace('bulk_viscosity', 'energy_solve 0.01 24\r\nbulk_viscosity');
  const changed = editNativeGasRegion(authored, nativeGasRegions(authored)[0], thermochemical);
  assert.equal(changed, authored.replace('gas dry-air 300 100000',
    'gas thermochemical 300 100000 mass 1 original-N2 1'));
  assert.throws(() => editNativeGasRegion(changed, nativeGasRegions(changed)[0], primary.values),
    /remove the thermochemical energy_solve.*performance editor/);
  const removed = changed.replace('energy_solve 0.01 24\r\n', '');
  assert.equal(editNativeGasRegion(removed, nativeGasRegions(removed)[0], primary.values), source);
});

test('resolved composition edits retain original geometry, clocks, source attachments and solver controls', async () => {
  const source = (await fixture('resolved-dg.performance'))
    .replace('gas dry-air 300 100000', 'gas thermochemical 300 100000 mole 2 source-N2 0.8 source-O2 0.2')
    .replace('bulk_viscosity', 'energy_solve 0.01 24\nbulk_viscosity');
  const components = ['source-N2', 'source-O2', 'source-He'].map(id => ({
    speciesPath:`gas/${id}.species`, nasa9Path:`gas/${id}.nasa9`, transportPath:`gas/${id}.transport`}));
  const options = {renderer:'resolved', decimate:true, gasSources:{phase:'frozen', waterSpecies:null, components}};
  const files = components.flatMap(component => Object.values(component)).map((path, i) =>
    ({path, bytes:new Uint8Array([0, 255, i, 13, 10])}));
  const primary = nativeGasRegions(source)[0];
  const changed = editNativeGasRegion(source, primary, {...primary.values, temperature:'310',
    components:[{species:'source-N2', fraction:'0.6'}, {species:'source-O2', fraction:'0.2'}, {species:'source-He', fraction:'0.2'}]});
  const expected = source.replace('gas thermochemical 300 100000 mole 2 source-N2 0.8 source-O2 0.2',
    'gas thermochemical 310 100000 mole 3 source-N2 0.6 source-O2 0.2 source-He 0.2');
  assert.equal(changed, expected);
  const packet = validateInput(new TextEncoder().encode(serializeNativeExperiment(changed, files, options)), true);
  assert.equal(new TextDecoder().decode(packet.source), expected);
  assert.deepEqual(packet.attachments, files);
  assert.deepEqual(packet.gasSources, options.gasSources);
  assert.equal(packet.renderer, 'resolved');
  assert.equal(packet.decimate, true);
});
