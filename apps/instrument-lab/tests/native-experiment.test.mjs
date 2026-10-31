// Exact source-packet round trips, not a surrogate physical model.
import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import {posix} from 'node:path';
import {editableAttachmentText, editedAttachmentBytes, prepareNativeTextExperiment, serializeNativeExperiment, validateInput,
  MAX_SOURCE_BYTES, MAX_ATTACHMENTS, MAX_ATTACHMENT_FILE_BYTES, MAX_ATTACHMENT_BYTES, MAX_REQUEST_BYTES} from '../src/native-experiment.mjs';
import {renderNativePerformance} from '../src/native-wind-client.mjs';
import {validateInput as serviceInput} from '../tools/native-wind-protocol.mjs';

const encode = text => new TextEncoder().encode(text);
const decode = bytes => new TextDecoder('utf-8', {fatal:true, ignoreBOM:true}).decode(bytes);
const root = new URL('../../../', import.meta.url);
const readOriginal = async path => decode(await readFile(new URL(path, root)));

test('complete canonical string and reed starters reach the native request with original physical sources', async () => {
  const strings = JSON.parse(await readFile(new URL('../src/native-string-starters.json', import.meta.url), 'utf8'));
  const reeds = JSON.parse(await readFile(new URL('../src/native-reed-starters.json', import.meta.url), 'utf8'));
  assert.equal(Object.keys(reeds).length, 7);
  assert.ok(reeds['material_modal_valve/paired-spatial.performance']);
  const starters = {...strings, ...reeds};
  assert.equal(Object.keys(starters).length, Object.keys(strings).length + Object.keys(reeds).length);
  const status = {available:true, renderers:{wind:true, strings:true}, token:'test-session',
    maximumSourceBytes:MAX_SOURCE_BYTES, maximumOutputBytes:16 * 1024 * 1024,
    maximumAttachments:MAX_ATTACHMENTS, maximumAttachmentFileBytes:MAX_ATTACHMENT_FILE_BYTES,
    maximumAttachmentBytes:MAX_ATTACHMENT_BYTES, maximumRequestBytes:MAX_REQUEST_BYTES,
    supportsMaterialManifests:true};
  for (const [name, inputs] of Object.entries(starters)) {
    const source = await readOriginal('crates/fs-couple/examples/' + name);
    const reed = Object.hasOwn(reeds, name);
    const decimate = reed && !['material_modal_valve/paired-walls.performance',
      'material_modal_valve/paired-spatial.performance'].includes(name);
    const options = {renderer:reed ? 'wind' : 'strings', decimate, materialManifests:inputs.materialManifests,
      stringProfiles:inputs.stringProfiles};
    const prepared = await prepareNativeTextExperiment(source, inputs.sources, options, readOriginal);
    const files = new Map(prepared.attachments.map(file => [file.path, decode(file.bytes)]));
    assert.equal(decode(prepared.source), source, name);
    for (const [path, original] of Object.entries(inputs.sources))
      assert.equal(files.get(path), await readOriginal(original), `${name}: ${path}`);
    for (const binding of prepared.materialManifests) {
      // Check the original manifest's relative source references; no source
      // parsing or material coefficients are added to browser production code.
      const manifest = files.get(binding.manifestPath);
      if (reed) assert.deepEqual(binding.materialState, {chemistry:'synthetic', phase:'solid',
        process:binding.packPath === 'wall.matpack' ? 'manufactured-membrane' : 'manufactured', revision:0},
        `${name}: full material card identity must reach the compiler`);
      for (const row of manifest.split(/\r?\n/).filter(line => line.startsWith('source\t'))) {
        const sourcePath = posix.join(posix.dirname(binding.manifestPath), row.split('\t')[2]);
        assert.ok(files.has(sourcePath), `${name}: compiler source ${sourcePath}`);
      }
    }
    const packs = new Set(prepared.materialManifests.map(binding => binding.packPath));
    for (const text of [source, ...Object.values(inputs.stringProfiles ?? {}).map(path => files.get(path))])
      for (const pack of text.match(/[^\s]+\.(?:fsmatpk|matpack)\b/g) ?? [])
        assert.ok(packs.has(pack), `${name}: physical source ${pack}`);
    let submitted;
    await assert.rejects(renderNativePerformance(source, prepared.attachments, status, prepared, undefined,
      async (url, request) => {
        assert.equal(url, '/api/native-wind/render');
        assert.equal(request.method, 'POST');
        submitted = serviceInput(encode(request.body), true);
        // Deliberate native-boundary refusal: this test proves source delivery,
        // not execution of the physical solver or an artificial sound.
        return new Response(JSON.stringify({error:'native-boundary witness'}), {
          status:400, headers:{'Content-Type':'application/json'},
        });
      }), /native-boundary witness/);
    assert.deepEqual(submitted, prepared, name);
    assert.deepEqual(validateInput(encode(serializeNativeExperiment(source, prepared.attachments, prepared)), true), prepared);
  }
});

test('reed and wall source edits survive saving without replacing the coupled performance or other materials', async () => {
  const starters = JSON.parse(await readFile(new URL('../src/native-reed-starters.json', import.meta.url), 'utf8'));
  const name = 'material_modal_valve/paired-walls.performance', inputs = starters[name];
  const source = await readOriginal('crates/fs-couple/examples/' + name);
  const before = await prepareNativeTextExperiment(source, inputs.sources,
    {renderer:'wind', decimate:false, materialManifests:inputs.materialManifests}, readOriginal);
  const originals = new Map(before.attachments.map(file => [file.path, file.bytes]));
  const originalBlade = decode(originals.get('blade.tsv'));
  const editedBlade = originalBlade.replace('young_modulus\tsynthetic\tyoung_modulus\t4000000000\tPa',
    'young_modulus\tsynthetic\tyoung_modulus\t3200000000\tPa');
  assert.notEqual(editedBlade, originalBlade);
  const files = before.attachments.map(file => file.path === 'blade.tsv'
    ? {...file, bytes:editedAttachmentBytes(file.bytes, editedBlade)} : file);
  const stored = serializeNativeExperiment(source, files, before);
  const restored = serviceInput(encode(stored), true);
  assert.equal(decode(restored.source), source);
  assert.equal(restored.decimate, false);
  assert.deepEqual(restored.materialManifests, before.materialManifests);
  for (const file of restored.attachments) {
    if (file.path === 'blade.tsv') assert.equal(decode(file.bytes), editedBlade);
    else assert.deepEqual(file.bytes, originals.get(file.path), file.path);
  }
  assert.deepEqual(validateInput(encode(serializeNativeExperiment(source, restored.attachments, restored)), true), restored);
});

test('starter source failures cannot return a partial experiment or discard required profile bindings', async () => {
  const source = 'frankensim-material-strings-v3\n';
  const options = {renderer:'strings', stringProfiles:{bodyPath:'body.profile'}};
  const sources = {'wire.tsv':'original-wire', 'body.profile':'original-body'};
  await assert.rejects(prepareNativeTextExperiment(source, sources, options, async path => {
    if (path === 'original-body') throw new Error('original source unavailable');
    return 'original wire';
  }), /original source unavailable/);
  await assert.rejects(prepareNativeTextExperiment(source, sources, options, async () => undefined), /Could not load original source/);
  await assert.rejects(prepareNativeTextExperiment(source, {'wire.tsv':'original-wire'}, options,
    async () => 'original wire'), /String profiles/);
  const before = await prepareNativeTextExperiment(source, sources, options, async path => path + '\r\n');
  const changed = await prepareNativeTextExperiment(source, sources, options,
    async path => (path === 'original-wire' ? 'changed material table' : path) + '\r\n');
  assert.equal(decode(before.attachments[0].bytes), 'original-wire\r\n');
  assert.equal(decode(changed.attachments[0].bytes), 'changed material table\r\n');
  assert.deepEqual(changed.source, before.source);
  assert.deepEqual(changed.stringProfiles, before.stringProfiles);
});

test('complete string experiments retain original moving contact, thermal, material and gas data', () => {
  const source = 'frankensim-material-strings-v3\r\n# µ original specimen\r\naudio 144000 432000 5 128\r\n'
    + 'state 307 1000\r\ngas thermochemical 298.15 82000 mole 2 oxygen-source 0.3 inert-source 0.7\r\n'
    + '# bow-stations and all following physical records are opaque transport bytes\r\n';
  const options = {renderer:'strings', decimate:true,
    stringProfiles:{bodyPath:'specimen/body.profile', stopPath:'contacts/bow.stops', contactThermalPath:'contacts/heat.profile'},
    coatingStrainGuard:{maximumBytes:536870912, maximumWork:1000000000000},
    gasSources:{components:[{speciesPath:'gas/oxygen.species', nasa9Path:'gas/oxygen.nasa9', transportPath:'gas/oxygen.transport'},
      {speciesPath:'gas/inert.species', nasa9Path:'gas/inert.nasa9', transportPath:'gas/inert.transport'}],
      waterSpecies:null, phase:'frozen'}};
  const paths = ['materials/wire.fsmatpk', 'materials/receiver.fsmatpk',
    'contacts/normal.fsintpk', 'contacts/normal.fsmodelpk', ...Object.values(options.stringProfiles),
    ...options.gasSources.components.flatMap(component => Object.values(component))];
  const files = paths.map((path, i) => ({path, bytes:new Uint8Array([0, 255, i, 13, 10, 128])}));
  const stored = serializeNativeExperiment(source, files, options);
  const restored = validateInput(encode(stored), true);
  assert.equal(serviceInput, validateInput, 'the browser and server must admit the same envelope');
  assert.equal(decode(restored.source), source);
  assert.deepEqual(restored.attachments, files);
  assert.equal(restored.renderer, 'strings');
  assert.equal(restored.decimate, true);
  assert.deepEqual(restored.stringProfiles, options.stringProfiles);
  assert.deepEqual(restored.coatingStrainGuard, options.coatingStrainGuard);
  assert.deepEqual(restored.gasSources, options.gasSources);
  assert.deepEqual(validateInput(encode(serializeNativeExperiment(decode(restored.source), restored.attachments, restored)), true), restored);
  // A later data-only material substitution cannot change a saved experiment.
  files[0].bytes[0] = 17;
  const variant = validateInput(encode(serializeNativeExperiment(source, files, options)), true);
  assert.equal(variant.attachments[0].bytes[0], 17);
  assert.equal(restored.attachments[0].bytes[0], 0);
  assert.deepEqual(variant.source, restored.source);
  assert.deepEqual(variant.stringProfiles, restored.stringProfiles);
  assert.deepEqual(variant.coatingStrainGuard, restored.coatingStrainGuard);
});

test('coating strain budgets require explicit exact limits and a material-string body profile', () => {
  const guard = {maximumBytes:536870912, maximumWork:1000000000000};
  const files = [{path:'body.profile', bytes:encode('original coated body')}];
  const packet = JSON.parse(serializeNativeExperiment('frankensim-material-strings-v3\n', files,
    {renderer:'strings', stringProfiles:{bodyPath:'body.profile'}}));
  for (const value of [null, [], {}, {...guard, extra:1},
    ...[0, -1, 1.5, '2', Number.MAX_SAFE_INTEGER + 1].flatMap(invalid =>
      [{...guard, maximumBytes:invalid}, {...guard, maximumWork:invalid}])]) {
    assert.throws(() => validateInput(encode(JSON.stringify({...packet, coatingStrainGuard:value})), true),
      /Coating strain checks/);
  }
  const {stringProfiles, ...withoutBody} = packet;
  assert.throws(() => validateInput(encode(JSON.stringify({...withoutBody, coatingStrainGuard:guard})), true),
    /Coating strain checks/);
  const wind = serializeNativeExperiment('frankensim-reed-performance-v1\n', [],
    {renderer:'wind', coatingStrainGuard:guard});
  assert.throws(() => validateInput(encode(wind), true), /Coating strain checks/);
  assert.equal(validateInput(encode(JSON.stringify(packet)), true).coatingStrainGuard, undefined);
});

test('wind experiments preserve independently declared internal, regional and external atmospheres', () => {
  const source = 'frankensim-material-modal-dg-valve-v1\n'
    + '# The native parser alone owns the following regional and external physics.\n'
    + 'ambient mixture 310 90000 mole 2 O2 0.25 N2 0.75\n'
    + 'exterior mixture 290 90000 mole 2 N2 0.9 He 0.1\n';
  const packet = serializeNativeExperiment(source, [], {renderer:'wind', decimate:true});
  const restored = validateInput(encode(packet), true);
  assert.equal(decode(restored.source), source);
  assert.equal(restored.decimate, true);
  assert.equal(restored.renderer ?? 'wind', 'wind');
  assert.deepEqual(restored.attachments, []);
  assert.equal(restored.gasSources, undefined, 'built-in mixture declarations do not acquire guessed source associations');
  assert.deepEqual(validateInput(encode(serializeNativeExperiment(decode(restored.source), restored.attachments, restored)), true), restored);
});

test('an incomplete experiment cannot drop material profiles or gas source bindings while opening', () => {
  const packet = JSON.parse(serializeNativeExperiment('frankensim-material-strings-v3\n',
    [{path:'wire.fsmatpk', bytes:new Uint8Array([1, 2, 3])}], {renderer:'strings'}));
  assert.throws(() => validateInput(encode(JSON.stringify({...packet, stringProfiles:{stopPath:'missing.stops'}})), true), /String profiles/);
  assert.throws(() => validateInput(encode(JSON.stringify({...packet, gasSources:{phase:'frozen', waterSpecies:null,
    components:[{speciesPath:'wire.fsmatpk', nasa9Path:'missing.model', transportPath:'wire.fsmatpk'}]}})), true), /attached source files/);
  assert.throws(() => validateInput(encode(JSON.stringify({...packet, attachments:[{path:'../wire.fsmatpk', bytes:'AQID'}]})), true), /relative paths/);
  assert.throws(() => validateInput(encode(JSON.stringify({...packet, renderer:'wind'})), true), /native reed/);
});

test('opening or applying an unchanged text source preserves every original byte', () => {
  const original = encode('\ufefffrankensim-string-stops-v1\r\n# µ specimen\r\n');
  const snapshot = original.slice();
  const text = editableAttachmentText(original);
  assert.equal(text, '\ufefffrankensim-string-stops-v1\r\n# µ specimen\r\n');
  assert.deepEqual(original, snapshot, 'viewing and cancelling edits cannot mutate the attachment');
  assert.equal(editedAttachmentBytes(original, text), original);
  assert.equal(editedAttachmentBytes(original, text.replace(/\r\n?/g, '\n')), original,
    'a textarea normalizing CRLF must not replace untouched source bytes');
});

test('an edited material-stop profile round-trips with its original bulk interface and model packs', () => {
  const path = 'contacts/fingers.profile';
  const source = 'frankensim-material-strings-v3\naudio 144000 432000 5 128\n';
  const text = 'frankensim-string-stops-v1\n'
    + 'material-stop finger 0 0.25 positive fixed -1e-9 materials/receiver.fsmatpk contacts/normal.fsintpk contacts/normal.fsmodelpk 5000 5000 0.1 0.00005 0.00005 0.01 0.02 0.03 10 1e-8\n';
  const files = [{path, bytes:encode(text)},
    ...['materials/wire.fsmatpk', 'materials/receiver.fsmatpk',
      'contacts/normal.fsintpk', 'contacts/normal.fsmodelpk'].map((path, i) =>
      ({path, bytes:new Uint8Array([0, 255, i, 13, 10, 128])}))];
  const options = {renderer:'strings', decimate:true, stringProfiles:{stopPath:path}};
  const before = validateInput(encode(serializeNativeExperiment(source, files, options)), true);
  const changed = editedAttachmentBytes(files[0].bytes, text.replace('0.25 positive', '0.35 positive'));
  const after = validateInput(encode(serializeNativeExperiment(source,
    [{path, bytes:changed}, ...files.slice(1)], options)), true);
  assert.equal(decode(after.attachments[0].bytes), text.replace('0.25 positive', '0.35 positive'));
  assert.equal(decode(before.attachments[0].bytes), text);
  assert.equal(decode(files[0].bytes), text, 'the edit creates a candidate without changing the old attachment');
  assert.deepEqual(after.attachments.slice(1), before.attachments.slice(1));
  assert.deepEqual(after.source, before.source);
  assert.deepEqual(after.stringProfiles, before.stringProfiles);
  assert.equal(after.decimate, true);
});

test('binary source packs and invalid edited text cannot enter the text editor', () => {
  for (const bytes of [new Uint8Array([255]), encode('FSMODELPK\u0000binary'), new Uint8Array()]) {
    assert.throws(() => editableAttachmentText(bytes), /UTF-8|binary|nonempty/);
  }
  const original = encode('frankensim-contact-thermal-v1\n');
  for (const text of ['', 'profile\u0000data']) assert.throws(() => editedAttachmentBytes(original, text), /nonempty|binary/);
  assert.equal(decode(original), 'frankensim-contact-thermal-v1\n');
});

test('resolved experiments retain their performer, original gas sources and physical field files', () => {
  const source = 'frankensim-resolved-dg-performance-v1\r\n# original geometry and fields µ\r\n';
  const gasSources = {components:[{speciesPath:'gas/species', nasa9Path:'gas/nasa9',
    transportPath:'gas/transport'}], waterSpecies:null, phase:'frozen'};
  const attachments = ['fields/initial.txt', 'wall.fsmatpk', ...Object.values(gasSources.components[0])]
    .map((path, index) => ({path, bytes:new Uint8Array([0, 255, index, 13, 10])}));
  const options = {renderer:'resolved', decimate:true, gasSources};
  const restored = validateInput(encode(serializeNativeExperiment(source, attachments, options)), true);
  assert.equal(restored.renderer, 'resolved');
  assert.equal(decode(restored.source), source);
  assert.deepEqual(restored.attachments, attachments);
  assert.deepEqual(restored.gasSources, gasSources);
  assert.equal(restored.decimate, true);
  assert.deepEqual(validateInput(encode(serializeNativeExperiment(decode(restored.source), restored.attachments, restored)), true), restored);
});

test('finite-gas string experiments retain exact solid, atmosphere, reservoir and feedback inputs', () => {
  const source = 'frankensim-finite-gas-string-v1\r\naudio 48000 2400 0.01\r\n'
    + 'nozzle gas/reservoir.gas 0.000002 0.0000000001 1 1e-14 1e-9 0.05\r\n'
    + 'drive current 0.1\r\nelectrical_feedback 0.000005 1e-12 64\r\n';
  const options = {renderer:'finite-gas-string',
    finiteGasString:{solidPath:'materials/wire.pack', gasRecordPath:'gas/local.gas'},
    gasSources:{components:[{speciesPath:'gas/species', nasa9Path:'gas/nasa9', transportPath:'gas/transport'}],
      waterSpecies:null, phase:'frozen'}};
  const files = [...Object.values(options.finiteGasString), 'gas/reservoir.gas',
    ...Object.values(options.gasSources.components[0])]
    .map((path, i) => ({path, bytes:new Uint8Array([0, 255, i, 13, 10])}));
  const restored = validateInput(encode(serializeNativeExperiment(source, files, options)), true);
  assert.equal(restored.renderer, 'finite-gas-string');
  assert.equal(decode(restored.source), source);
  assert.equal(restored.decimate, false);
  assert.deepEqual(restored.finiteGasString, options.finiteGasString);
  assert.deepEqual(restored.gasSources, options.gasSources);
  assert.deepEqual(restored.attachments, files);
  assert.deepEqual(validateInput(encode(serializeNativeExperiment(decode(restored.source), restored.attachments, restored)), true), restored);
});

test('finite-gas string packets require complete frozen sources and preserve native output names', () => {
  const source = 'frankensim-finite-gas-string-v1\n';
  const files = ['solid', 'local.gas', 'species', 'nasa9', 'transport']
    .map(path => ({path, bytes:new Uint8Array([1])}));
  const options = {renderer:'finite-gas-string', finiteGasString:{solidPath:'solid', gasRecordPath:'local.gas'},
    gasSources:{components:[{speciesPath:'species', nasa9Path:'nasa9', transportPath:'transport'}], waterSpecies:null, phase:'frozen'}};
  const packet = JSON.parse(serializeNativeExperiment(source, files, options));
  const admit = value => validateInput(encode(JSON.stringify(value)), true);
  for (const changes of [{finiteGasString:undefined}, {gasSources:undefined}, {decimate:true},
    {gasSources:{...options.gasSources, phase:'liquid-water', waterSpecies:'H2O'}},
    {finiteGasString:{solidPath:'missing', gasRecordPath:'local.gas'}},
    {finiteGasString:{solidPath:'solid', gasRecordPath:'../local.gas'}},
    {finiteGasString:{...options.finiteGasString, gain:2}}, {stringProfiles:{}}]) {
    assert.throws(() => admit({...packet, ...changes}));
  }
  for (const path of ['output.history.csv', 'output.history.csv/child'])
    assert.throws(() => admit({...packet, attachments:[...packet.attachments, {path, bytes:'AQ=='}]}), /collide/);
  assert.throws(() => admit({...packet, renderer:'wind',
    source:Buffer.from('frankensim-reed-performance-v1\n').toString('base64')}));
});

test('material manifest experiments retain original source data, explicit identities and generated gas bindings', () => {
  const source = 'frankensim-material-strings-v3\r\n# unchanged measured setup and temperature_known 0\r\n';
  const files = [
    {path:'sources/manifest.tsv', bytes:encode('original manifest\r\n')},
    {path:'sources/properties.tsv', bytes:encode('original source table · µ\r\n')},
  ];
  const materialState = {chemistry:'authored mixture', phase:'solid', process:'measured-state', revision:0};
  const materialManifests = [
    {packPath:'materials/string.pack', manifestPath:files[0].path},
    {packPath:'gas/transport.pack', manifestPath:files[0].path, materialState},
    {packPath:'gas/species.pack', manifestPath:files[0].path},
    {packPath:'gas/nasa9.pack', manifestPath:files[0].path},
  ];
  const options = {renderer:'strings', materialManifests,
    gasSources:{components:[{speciesPath:'gas/species.pack', nasa9Path:'gas/nasa9.pack',
      transportPath:'gas/transport.pack'}], phase:'frozen', waterSpecies:null}};
  const saved = serializeNativeExperiment(source, files, options);
  const restored = validateInput(encode(saved), true);
  assert.deepEqual(restored.materialManifests, materialManifests);
  assert.deepEqual(restored.attachments, files);
  assert.deepEqual(restored.gasSources, options.gasSources);
  assert.equal(decode(restored.source), source);
  assert.deepEqual(validateInput(encode(serializeNativeExperiment(source, restored.attachments, restored)), true), restored);
  files[1].bytes = encode('authored material-table variant\n');
  const changed = validateInput(encode(serializeNativeExperiment(source, files, options)), true);
  assert.notDeepEqual(changed.attachments[1].bytes, restored.attachments[1].bytes);
  assert.deepEqual(changed.source, restored.source);
  assert.deepEqual(changed.materialManifests, restored.materialManifests);
});

test('material compilation requires original manifest inputs, new confined targets and explicit complete identities', () => {
  const source = 'frankensim-material-modal-valve-v1\n';
  const files = [{path:'source/manifest.tsv', bytes:encode('original source')}];
  const binding = {packPath:'reed.pack', manifestPath:files[0].path};
  const check = value => validateInput(encode(serializeNativeExperiment(source, files,
    {materialManifests:value})), true);
  for (const value of [[], null, {}, [null], [{...binding, manifestPath:'absent.tsv'}],
    [{...binding, packPath:files[0].path}], [binding, binding],
    [{...binding, packPath:'output.wav'}], [{...binding, packPath:'../reed.pack'}],
    [{...binding, packPath:'source'}], [{...binding, packPath:'reed.pack/part'}, binding],
    [binding, {packPath:'other.pack', manifestPath:'reed.pack'}],
    [binding, ...Array.from({length:31}, (_, i) => ({...binding, packPath:'pack-' + i}))],
    ...[{}, {chemistry:'c', phase:'s', process:'p'}, {chemistry:'c', phase:'s', process:'p', revision:-1},
      {chemistry:'c', phase:'s', process:'p', revision:4294967296},
      {chemistry:'c', phase:'s', process:'p', revision:'0'}].map(materialState => [{...binding, materialState}])]) {
    assert.throws(() => check(value));
  }
  assert.equal(check([binding]).materialManifests[0].materialState, undefined,
    'plain source claims, models, interfaces and species acquire no guessed material identity');
  const old = validateInput(encode(serializeNativeExperiment(source, files)), true);
  assert.equal(old.materialManifests, undefined);
});
