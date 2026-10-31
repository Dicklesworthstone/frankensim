// Source-edit and original-packet regressions; native physics owns gas evaluation.
import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import {nativeGasRegions, editNativeGasRegion} from '../src/native-gas-records.mjs';
import {editableAttachmentText, editedAttachmentBytes, serializeNativeExperiment, validateInput} from '../src/native-experiment.mjs';

const fixture = name => readFile(new URL(`../../../crates/fs-couple/examples/${name}`, import.meta.url), 'utf8');


test('attached finite gas records expose their original state without inventing a performance owner', async () => {
  const source = await fixture('finite_gas_string/nitrogen-reservoir.gas');
  const [region] = nativeGasRegions(source);
  assert.equal(nativeGasRegions(source).length, 1);
  assert.deepEqual(region.supportedKinds, ['thermochemical']);
  assert.equal(region.values.temperature, '300');
  assert.equal(region.values.pressure, '95000');
  assert.deepEqual(region.values.components, [{species:'N2', fraction:'1'}]);
  assert.equal(editNativeGasRegion(source, region, region.values), source);
  assert.equal(editNativeGasRegion(source, region, {...region.values,
    temperature:'310', pressure:'90000',
    components:[{species:'N2', fraction:'0.8'}, {species:'CO2', fraction:'0.2'}]}),
    'thermochemical 310 90000 mole 2 N2 0.8 CO2 0.2\n');
});

test('bare gas edits retain comments, BOM, line endings and the explicitly selected fraction basis', () => {
  const source = '\ufeff# original reservoir, not a new material card\r\n\tthermochemical\t300 95000 mass 1 N2 1\r\n# original source IDs\r\n';
  const [region] = nativeGasRegions(source);
  assert.equal(region.line, 2);
  assert.equal(region.values.basis, 'mass');
  assert.equal(editNativeGasRegion(source, region, region.values), source);
  assert.equal(editNativeGasRegion(source, region, {...region.values, pressure:'88000'}),
    '\ufeff# original reservoir, not a new material card\r\n\tthermochemical 300 88000 mass 1 N2 1\r\n# original source IDs\r\n');
  const withBom = '\ufeffthermochemical 300 95000 mole 1 N2 1\r\n';
  const [first] = nativeGasRegions(withBom);
  assert.equal(editNativeGasRegion(withBom, first, {...first.values, temperature:'305'}),
    '\ufeffthermochemical 305 95000 mole 1 N2 1\r\n');
});

test('finite gas sources keep their native thermochemical law and reject ambiguous complete records', () => {
  const source = 'thermochemical 300 95000 mole 1 N2 1\n';
  const [region] = nativeGasRegions(source);
  assert.throws(() => editNativeGasRegion(source, region, {...region.values, kind:'mixture'}), /supports thermochemical/);
  assert.throws(() => editNativeGasRegion(source, region, {...region.values, kind:'dry-air'}), /supports thermochemical/);
  assert.throws(() => editNativeGasRegion(source, region, {...region.values, pressure:'0'}), /Absolute pressure/);
  for (const invalid of [
    source + source, source + 'unclaimed row\n', source.trim() + ' trailing\n',
    'thermochemical 300 95000 mole 2 N2 1\n',
    'mixture 300 95000 mole 1 N2 1\n',
  ]) assert.deepEqual(nativeGasRegions(invalid), []);
  assert.throws(() => editNativeGasRegion(source + '# changed\n', region, region.values), /changed/);
});


test('edited finite atmosphere reaches the saved native experiment with independent reservoir bytes and source bindings', async () => {
  const encode = text => new TextEncoder().encode(text);
  const source = await fixture('finite_gas_string/browser-feedback-flow.string');
  const initial = encode('thermochemical 300 101325 mole 1 N2 1\r\n');
  const reservoir = encode(await fixture('finite_gas_string/nitrogen-reservoir.gas'));
  const original = editableAttachmentText(initial);
  const [region] = nativeGasRegions(original);
  const changed = editNativeGasRegion(original, region, {...region.values, temperature:'315', pressure:'88000'});
  const edited = editedAttachmentBytes(initial, changed);
  const components = [{speciesPath:'species.pack', nasa9Path:'nasa9.pack', transportPath:'transport.pack'}];
  const files = [
    {path:'string.pack', bytes:new Uint8Array([0, 1, 128, 255])},
    {path:'nitrogen.gas', bytes:edited},
    {path:'nitrogen-reservoir.gas', bytes:reservoir},
    ...Object.values(components[0]).map(path => ({path, bytes:new Uint8Array([255, 0, 127])})),
  ];
  const options = {renderer:'finite-gas-string', finiteGasString:{solidPath:'string.pack', gasRecordPath:'nitrogen.gas'},
    gasSources:{components, phase:'frozen', waterSpecies:null}};
  const packet = validateInput(encode(serializeNativeExperiment(source, files, options)), true);
  assert.deepEqual(packet.source, encode(source));
  assert.deepEqual(packet.finiteGasString, options.finiteGasString);
  assert.deepEqual(packet.gasSources, options.gasSources);
  assert.deepEqual(packet.attachments, files);
  assert.equal(editableAttachmentText(packet.attachments[1].bytes),
    'thermochemical 315 88000 mole 1 N2 1\r\n');
  assert.deepEqual(packet.attachments[2].bytes, reservoir);
  assert.equal(editedAttachmentBytes(initial, original.replace(/\r\n/g, '\n')), initial);
});

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
  assert.deepEqual(primary.supportedKinds, ['dry-air', 'mixture', 'thermochemical']);
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

test('resolved primary gas refuses moist-air snapshots while exterior models retain their original choices', async () => {
  const source = await fixture('resolved-dg.performance');
  const [primary, exterior] = nativeGasRegions(source);
  const moist = {kind:'moist-air', temperature:'300', pressure:'100000', humidity:'0', components:[]};
  assert.throws(() => editNativeGasRegion(source, primary, moist), /supports dry-air or mixture or thermochemical/);
  assert.notEqual(editNativeGasRegion(source, exterior, {...moist, temperature:'310'}), source);
  assert.throws(() => editNativeGasRegion(source, exterior, primary.values), /supports moist-air or mixture or thermochemical/);
  const changed = editNativeGasRegion(source, exterior, {kind:'thermochemical', temperature:'300',
    pressure:'100000', basis:'mass', components:[{species:'original-N2', fraction:'1'}]});
  assert.equal(changed, source.replace('dg_exterior 300 100000 0',
    'dg_exterior thermochemical 300 100000 mass 1 original-N2 1'));
});

test('resolved named mixtures reach both original atmosphere records and portable experiments without rewriting physics', async () => {
  const original = await fixture('resolved-dg-radiating-mouth.performance');
  // The existing numerical window belongs to the author, so a model switch
  // first identifies its incompatibility instead of quietly changing limits.
  const mixture = {kind:'mixture', temperature:'310', pressure:'90000', basis:'mole',
    components:[{species:'N2', fraction:'0.7'}, {species:'O2', fraction:'0.2'}, {species:'He', fraction:'0.1'}]};
  assert.throws(() => editNativeGasRegion(original, nativeGasRegions(original)[0], mixture),
    /fluid_limits.*MINIMUM_TEMPERATURE_K.*MAXIMUM_TEMPERATURE_K.*Discard changes/);
  const source = original.replace('fluid_limits 144 250 600', 'fluid_limits 144 280 320');
  assert.notEqual(source, original);
  for (const basis of ['mole', 'mass']) {
    const values = {...mixture, basis};
    const primary = nativeGasRegions(source).find(region => region.key === 'gas');
    const inside = editNativeGasRegion(source, primary, values, true);
    const exterior = nativeGasRegions(inside).find(region => region.key === 'dg_exterior');
    const changed = editNativeGasRegion(inside, exterior, values);
    const declaration = `mixture 310 90000 ${basis} 3 N2 0.7 O2 0.2 He 0.1`;
    assert.equal(changed, source.replace('gas dry-air 300 100000', `gas ${declaration}`)
      .replace('dg_exterior 300 100000 0', `dg_exterior ${declaration}`));
    assert.deepEqual(nativeGasRegions(changed).map(region => region.values), [values, values]);
    assert.equal(editNativeGasRegion(changed, nativeGasRegions(changed)[0], values), changed);
    assert.match(changed, /^reservoir_aperture 1 high reference 1$/m);
    assert.match(changed, /^reservoir_feedback -20 20 16 16$/m);
    const packet = validateInput(new TextEncoder().encode(serializeNativeExperiment(changed, [],
      {renderer:'resolved', decimate:true})), true);
    assert.equal(new TextDecoder().decode(packet.source), changed);
    assert.equal(packet.renderer, 'resolved');
    assert.equal(packet.decimate, true);
    assert.deepEqual(packet.attachments, []);
  }
});

test('resolved named mixtures retain temperature-window and thermochemical energy-solve ownership', async () => {
  const original = await fixture('resolved-dg.performance');
  const source = original.replace('fluid_limits 144 250 600', 'fluid_limits 144 280 320');
  const values = {kind:'mixture', temperature:'300', pressure:'100000', basis:'mass',
    components:[{species:'N2', fraction:'0.8'}, {species:'He', fraction:'0.2'}]};
  for (const row of ['', 'fluid_limits 144 280', 'fluid_limits 144 NaN 320 0.2 1 17 1',
    'fluid_limits 144 250 320 0.2 1 17 1', 'fluid_limits 144 280 600 0.2 1 17 1',
    'fluid_limits 144 320 280 0.2 1 17 1', 'fluid_limits 144 300 300 0.2 1 17 1',
    'fluid_limits 144 280 320 0.2 1 17 1\nfluid_limits 144 280 320 0.2 1 17 1']) {
    const invalid = source.replace('fluid_limits 144 280 320 0.2 1 17 1', row);
    assert.throws(() => editNativeGasRegion(invalid, nativeGasRegions(invalid)[0], values), /fluid_limits/);
  }
  const authored = source.replace('gas dry-air 300 100000',
    'gas thermochemical 300 100000 mole 1 source-N2 1')
    .replace('bulk_viscosity', 'energy_solve 0.01 24\nbulk_viscosity');
  assert.throws(() => editNativeGasRegion(authored, nativeGasRegions(authored)[0], values),
    /remove the thermochemical energy_solve.*performance editor/);
  const removed = authored.replace('energy_solve 0.01 24\n', '');
  assert.equal(editNativeGasRegion(removed, nativeGasRegions(removed)[0], values),
    source.replace('gas dry-air 300 100000', 'gas mixture 300 100000 mass 2 N2 0.8 He 0.2'));
});

test('source-bound reservoir fixtures retain their binding through common-pressure edits and native packets', async () => {
  for (const name of ['resolved-dg-radiating-mouth.performance', 'resolved-dg-side-aperture.performance',
    'resolved-dg-slit-edge.performance', 'resolved-dg-si-slit-edge.performance', 'resolved-dg-transverse-drive.performance']) {
    const source = await fixture(name);
    const binding = source.split('\n').find(line => /^reservoir_aperture(?:_m)? /.test(line));
    assert.match(binding, /^reservoir_aperture(?:_m)? [012] (?:low|high) reference /);
    const primary = nativeGasRegions(source).find(region => region.key === 'gas');
    const changed = editNativeGasRegion(source, primary, {...primary.values, pressure:'90000'}, true);
    const regions = nativeGasRegions(changed);
    assert(regions.every(region => region.values.pressure === '90000'));
    assert.equal(changed, source.replace('gas dry-air 300 100000', 'gas dry-air 300 90000')
      .replace('dg_exterior 300 100000 0', 'dg_exterior 300 90000 0'));
    const bytes = new TextEncoder().encode(serializeNativeExperiment(changed, [], {renderer:'resolved'}));
    const packet = validateInput(bytes, true);
    assert.equal(new TextDecoder().decode(packet.source), changed);
    assert.equal(packet.renderer, 'resolved');
    // The reference state is evaluated by the Rust material owner. The browser
    // does not calculate a replacement density or rewrite physical controls.
    assert(changed.split('\n').includes(binding));
  }
});

test('source-bound radiating mouth keeps the same aperture through thermochemical atmosphere substitutions', async () => {
  const source = (await fixture('resolved-dg-radiating-mouth.performance'))
    .replace('bulk_viscosity', 'energy_solve 0.01 24\nbulk_viscosity');
  const values = {kind:'thermochemical', temperature:'310', pressure:'90000', basis:'mole',
    components:[{species:'source-N2', fraction:'0.7'}, {species:'source-O2', fraction:'0.2'},
      {species:'source-He', fraction:'0.1'}]};
  const primary = nativeGasRegions(source).find(region => region.key === 'gas');
  const inside = editNativeGasRegion(source, primary, values, true);
  const exterior = nativeGasRegions(inside).find(region => region.key === 'dg_exterior');
  const changed = editNativeGasRegion(inside, exterior, values);
  const declaration = 'thermochemical 310 90000 mole 3 source-N2 0.7 source-O2 0.2 source-He 0.1';
  assert.equal(changed, source.replace('gas dry-air 300 100000', `gas ${declaration}`)
    .replace('dg_exterior 300 100000 0', `dg_exterior ${declaration}`));
  assert.match(changed, /^reservoir_aperture 1 high reference 1$/m);
  assert.match(changed, /^reservoir_feedback -20 20 16 16$/m);
  assert.deepEqual(nativeGasRegions(changed).map(region => region.values.components),
    [values.components, values.components]);
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
