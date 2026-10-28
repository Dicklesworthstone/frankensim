// Exact source-packet round trips, not a surrogate physical model.
import {test} from 'node:test';
import assert from 'node:assert/strict';
import {editableAttachmentText, editedAttachmentBytes, serializeNativeExperiment, validateInput} from '../src/native-experiment.mjs';
import {validateInput as serviceInput} from '../tools/native-wind-protocol.mjs';

const encode = text => new TextEncoder().encode(text);
const decode = bytes => new TextDecoder('utf-8', {fatal:true, ignoreBOM:true}).decode(bytes);

test('complete string experiments retain original moving contact, thermal, material and gas data', () => {
  const source = 'frankensim-material-strings-v3\r\n# µ original specimen\r\naudio 144000 432000 5 128\r\n'
    + 'state 307 1000\r\ngas thermochemical 298.15 82000 mole 2 oxygen-source 0.3 inert-source 0.7\r\n'
    + '# bow-stations and all following physical records are opaque transport bytes\r\n';
  const options = {renderer:'strings', decimate:true,
    stringProfiles:{bodyPath:'specimen/body.profile', stopPath:'contacts/bow.stops', contactThermalPath:'contacts/heat.profile'},
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
  assert.deepEqual(restored.gasSources, options.gasSources);
  assert.deepEqual(validateInput(encode(serializeNativeExperiment(decode(restored.source), restored.attachments, restored)), true), restored);
  // A later data-only material substitution cannot change a saved experiment.
  files[0].bytes[0] = 17;
  const variant = validateInput(encode(serializeNativeExperiment(source, files, options)), true);
  assert.equal(variant.attachments[0].bytes[0], 17);
  assert.equal(restored.attachments[0].bytes[0], 0);
  assert.deepEqual(variant.source, restored.source);
  assert.deepEqual(variant.stringProfiles, restored.stringProfiles);
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
