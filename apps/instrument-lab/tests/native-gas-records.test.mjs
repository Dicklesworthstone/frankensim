// Source-edit and original-packet regressions; native physics owns gas evaluation.
import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import {nativeGasRegions, editNativeGasRegion} from '../src/native-gas-records.mjs';
import {serializeNativeExperiment, validateInput} from '../src/native-experiment.mjs';

const fixture = name => readFile(new URL(`../../../crates/fs-couple/examples/${name}`, import.meta.url), 'utf8');

const fieldBore = [
  'dg_bore 2 open 2 2 8 50 2000 0.01 1e-11 1024 4 gas_field 3 mole 2 Z-N2 A-CO2',
  'dg_segment cylinder 0.001 0.02',
  '# original cone and independently declared station columns',
  'dg_segment cone 0.001 0.002 0.03',
  'gas_station\t0\t300\t1 0',
  'gas_station\t0.02\t310\t0.75 0.25',
  'gas_station\t0.05\t320\t0.5 0.5',
].join('\r\n') + '\r\n';
const upstreamField = [
  'upstream_tract',
  'dg_bore 1 closed 1 2 8 50 2000 0.01 1e-11 1024 4 gas_field 2 mole 2 A-CO2 Z-N2',
  'dg_segment cylinder 0.003 0.04',
  'gas_station 0 310 0.1 0.9',
  'gas_station 0.04 300 0.2 0.8',
].join('\r\n') + '\r\n';
const fieldSource = '\ufefffrankensim-material-modal-dg-valve-v1\r\n'
  + 'ambient thermochemical 300 101325 mole 2 Z-N2 0.8 A-CO2 0.2\r\n'
  + fieldBore + 'dg_molecular declared 2 Z-N2 0 A-CO2 0\r\n' + upstreamField;

test('complete main and upstream DG fields expose original station columns and preserve no-op CRLF', () => {
  const regions = nativeGasRegions(fieldSource), stations = regions.filter(region => region.field);
  assert.equal(regions.length, 6);
  assert.equal(stations.length, 5);
  assert.deepEqual(stations.map(region => region.field.position), ['0', '0.02', '0.05', '0', '0.04']);
  assert.deepEqual(stations.map(region => region.field.station), [1, 2, 3, 1, 2]);
  assert.deepEqual(stations[0].field.columns, ['Z-N2', 'A-CO2']);
  assert.deepEqual(stations[3].field.columns, ['A-CO2', 'Z-N2']);
  assert.notEqual(stations[0].field.headerLine, stations[3].field.headerLine);
  for (const region of regions) {
    assert.equal(editNativeGasRegion(fieldSource, region, region.values), fieldSource);
    if (region.field) {
      assert.equal(region.pressureSpan, null);
      assert.equal(region.values.pressure, '101325');
      assert.deepEqual(region.supportedKinds, ['thermochemical']);
    }
  }
});

test('a field station edit changes only its temperature and fractions, preserving header, geometry and upstream gas', () => {
  const selected = nativeGasRegions(fieldSource).filter(region => region.field)[1];
  const changed = editNativeGasRegion(fieldSource, selected, {...selected.values, temperature:'315',
    components:[{species:'Z-N2', fraction:'0.6'}, {species:'A-CO2', fraction:'0.4'}]});
  assert.equal(changed, fieldSource.replace('gas_station\t0.02\t310\t0.75 0.25',
    'gas_station\t0.02\t315 0.6 0.4'));
  const upstream = nativeGasRegions(changed).filter(region => region.field)[3];
  const twice = editNativeGasRegion(changed, upstream, {...upstream.values, temperature:'305'});
  assert.equal(twice, changed.replace('gas_station 0 310 0.1 0.9', 'gas_station 0 305 0.1 0.9'));
});

test('field station edits refuse changes to inherited pressure, basis and original species columns', () => {
  const selected = nativeGasRegions(fieldSource).find(region => region.field);
  assert.throws(() => editNativeGasRegion(fieldSource, selected, selected.values, true), /inherits static pressure/);
  assert.throws(() => editNativeGasRegion(fieldSource, selected, {...selected.values, pressure:'90000'}), /inherits static pressure/);
  for (const changed of [
    {basis:'mass'},
    {components:[...selected.values.components].reverse()},
    {components:[{species:'He', fraction:'1'}, {species:'A-CO2', fraction:'0'}]},
    {components:[{species:'Z-N2', fraction:'1'}]},
    {components:[...selected.values.components, {species:'He', fraction:'0'}]},
  ]) {
    assert.throws(() => editNativeGasRegion(fieldSource, selected, {...selected.values, ...changed}), /declared species columns/);
  }
  assert.throws(() => editNativeGasRegion(fieldSource, selected, {...selected.values, temperature:'0'}), /Temperature/);
  assert.throws(() => editNativeGasRegion(fieldSource, selected, {...selected.values,
    components:[{species:'Z-N2', fraction:'0.6'}, {species:'A-CO2', fraction:'0.5'}]}), /sum/);
});

test('common pressure changes the ambient owner while every original station row stays exact', () => {
  const ambient = nativeGasRegions(fieldSource)[0];
  const changed = editNativeGasRegion(fieldSource, ambient, {...ambient.values, pressure:'90000'}, true);
  assert.equal(changed, fieldSource.replace('300 101325 mole', '300 90000 mole'));
  const stations = nativeGasRegions(changed).filter(region => region.field);
  assert.equal(stations.length, 5);
  assert.ok(stations.every(region => region.values.pressure === '90000'));
});

test('malformed, oversized and unknown field blocks stay outside the station editor', () => {
  const ambient = 'ambient 300 101325 0\r\n';
  const malformed = [
    fieldBore.replace('gas_field 3 mole', 'gas_field 4 mole'),
    fieldBore.replace('gas_field 3 mole', 'gas_field 4097 mole'),
    fieldBore.replace('gas_field 3 mole 2', 'gas_field 3 mole 129'),
    fieldBore.replace('gas_field 3 mole', 'gas_field 3 mass'),
    fieldBore.replace('Z-N2 A-CO2', 'Z-N2 Z-N2'),
    fieldBore.replace('dg_segment cone', 'dg_segment profile'),
    fieldBore.replace('0.02\t310', '0\t310'),
    fieldBore.replace('0.05\t320', '0.06\t320'),
    fieldBore.replace('0.75 0.25', '0.75'),
    fieldBore.replace('0.75 0.25', '0.75 NaN'),
    fieldBore.replace('0.75 0.25', '0.75 0.35'),
    fieldBore + 'gas_station 0.05 300 0.5 0.5\r\n',
    'gas_field continuous 16 2 mole 2 Z-N2 A-CO2\r\ngas_station 0 300 1 0\r\ngas_station 0.05 320 0.5 0.5\r\n',
  ];
  for (const block of malformed) {
    const source = ambient + block, regions = nativeGasRegions(source);
    assert.equal(regions.filter(region => region.field).length, 0);
    assert.equal(editNativeGasRegion(source, regions[0], regions[0].values), source);
  }
  // A malformed first bore cannot donate its columns to a later valid tract.
  const later = nativeGasRegions(ambient + malformed[0] + upstreamField).filter(region => region.field);
  assert.equal(later.length, 2);
  assert.deepEqual(later[0].field.columns, ['A-CO2', 'Z-N2']);
});

test('field station experiments preserve raw source, both fields and every original source triple', () => {
  // The transport admits an exact first-line schema; BOM preservation is
  // covered separately by the source editor's no-op test above.
  const source = fieldSource.slice(1);
  const components = ['Z-N2', 'A-CO2'].map(id => ({
    speciesPath:'gas/' + id + '.species', nasa9Path:'gas/' + id + '.nasa9', transportPath:'gas/' + id + '.transport'}));
  const options = {renderer:'wind', decimate:true, gasSources:{phase:'frozen', waterSpecies:null, components}};
  const files = components.flatMap(component => Object.values(component)).map((path, i) =>
    ({path, bytes:new Uint8Array([0, 255, i, 13, 10])}));
  const selected = nativeGasRegions(source).filter(region => region.field)[4];
  const changed = editNativeGasRegion(source, selected, {...selected.values, temperature:'299'});
  const packet = validateInput(new TextEncoder().encode(serializeNativeExperiment(changed, files, options)), true);
  assert.equal(new TextDecoder().decode(packet.source), changed);
  assert.deepEqual(packet.attachments, files);
  assert.deepEqual(packet.gasSources, options.gasSources);
  assert.equal(packet.decimate, true);
});

test('one data editor reaches canonical flute, reed, chimney and string gases', async () => {
  const cases = [
    ['flue-dg-baffled.performance', ['gas', 'dg_exterior']],
    ['flue-dg-molecular.performance', ['gas', 'dg_exterior']],
    ['flue-dg-conical-molecular.performance', ['gas', 'dg_exterior']],
    ['flue-dg-field-molecular.performance', ['gas', 'gas_station', 'gas_station', 'gas_station', 'dg_exterior']],
    ['resolved-dg.performance', ['gas', 'dg_exterior']],
    ['reed-two-medium-tone-hole.performance', ['ambient', 'exterior', 'cylinder']],
    ['material_modal_valve/paired-dg-chimney.performance', ['ambient', 'dg_chimney_gas', 'dg_exterior']],
    ['material_strings/rayleigh-observed.material-strings', ['gas']],
    ['material_strings/regional-cavity.material-strings', ['gas', 'cavity_gas']],
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

test('conical molecular atmosphere edits preserve geometry, constitutive laws and saved source exactly', async () => {
  const source = await fixture('flue-dg-conical-molecular.performance');
  assert.match(source, /^dg_bore /m);
  assert.match(source, /^dg_segment cone /m);
  const [internal, exterior] = nativeGasRegions(source);
  assert.deepEqual([internal.key, exterior.key], ['gas', 'dg_exterior']);
  const changed = editNativeGasRegion(source, internal, {...internal.values, pressure:'85000'}, true);
  assert.equal(changed, source.replaceAll(' 101325 ', ' 85000 '));
  const packet = validateInput(new TextEncoder().encode(serializeNativeExperiment(changed, [], {
    renderer:'wind', decimate:true,
  })), true);
  assert.equal(new TextDecoder().decode(packet.source), changed);
  assert.equal(packet.decimate, true);
  assert.deepEqual(packet.attachments, []);
});

test('string cavity gas edits preserve the exterior, body, forces and saved performer', async () => {
  const source = await fixture('material_strings/regional-cavity.material-strings');
  const [exterior, cavity] = nativeGasRegions(source);
  assert.deepEqual([exterior.key, cavity.key], ['gas', 'cavity_gas']);
  const changed = editNativeGasRegion(source, cavity, {kind:'mixture', temperature:'303.15',
    pressure:'101325', basis:'mole', components:[{species:'He', fraction:'0.5'}, {species:'N2', fraction:'0.5'}]});
  assert.equal(changed, source.replace(/^cavity_gas .+$/m,
    'cavity_gas mixture 303.15 101325 mole 2 He 0.5 N2 0.5'));
  assert.deepEqual(nativeGasRegions(changed)[0].values, exterior.values);
  const packet = validateInput(new TextEncoder().encode(serializeNativeExperiment(changed, [], {
    renderer:'strings', decimate:false,
  })), true);
  assert.equal(new TextDecoder().decode(packet.source), changed);
  assert.equal(packet.renderer, 'strings');
  assert.equal(packet.decimate, false);
  const common = editNativeGasRegion(source, exterior, {...exterior.values, pressure:'85000'}, true);
  assert.equal(common, source.replaceAll(' 101325 ', ' 85000 '));
  assert.deepEqual(nativeGasRegions(common).map(region => region.values.pressure), ['85000', '85000']);
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
