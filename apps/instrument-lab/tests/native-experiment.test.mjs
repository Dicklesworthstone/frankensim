// Exact source-packet round trips, not a surrogate physical model.
import {test} from 'node:test';
import assert from 'node:assert/strict';
import {serializeNativeExperiment, validateInput} from '../src/native-experiment.mjs';
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
  const paths = ['materials/wire.fsmatpk', ...Object.values(options.stringProfiles),
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
