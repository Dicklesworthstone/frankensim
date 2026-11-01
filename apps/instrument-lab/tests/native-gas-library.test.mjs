// Original material data reaches the existing native protocol; no surrogate
// coefficients or synthetic instrument acoustics are evaluated in the browser.
import {test} from 'node:test';
import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import {posix} from 'node:path';
import {NATIVE_GAS_LIBRARY, loadNativeGasSource, mergeNativeGasSource} from '../src/native-gas-library.mjs';
import {prepareNativeTextExperiment, serializeNativeExperiment, validateInput,
  MAX_SOURCE_BYTES, MAX_ATTACHMENTS, MAX_ATTACHMENT_FILE_BYTES, MAX_ATTACHMENT_BYTES, MAX_REQUEST_BYTES} from '../src/native-experiment.mjs';
import {nativeGasRegions, editNativeGasRegion} from '../src/native-gas-records.mjs';
import {renderNativePerformance} from '../src/native-wind-client.mjs';
import {validateInput as serviceInput} from '../tools/native-wind-protocol.mjs';

const encode = text => new TextEncoder().encode(text);
const decode = bytes => new TextDecoder('utf-8', {fatal:true, ignoreBOM:true}).decode(bytes);
const root = new URL('../../../', import.meta.url);
const readOriginal = path => readFile(new URL(path, root), 'utf8');

test('all library gases retain complete canonical manifests, source tables and exact species IDs', async () => {
  let total = 0;
  assert.deepEqual(NATIVE_GAS_LIBRARY.map(entry => entry.species), ['N2', 'O2', 'CO2', 'He', 'Ar', 'Ne', 'SF6', 'H2O']);
  for (const entry of NATIVE_GAS_LIBRARY) {
    const originalFiles = new Map();
    const added = await loadNativeGasSource(entry.id, async path => {
      const text = await readOriginal(path); originalFiles.set(path, text); return text;
    });
    const files = new Map(added.attachments.map(file => [file.path, decode(file.bytes)]));
    assert.equal(files.size, 6);
    assert.equal(originalFiles.size, 6);
    assert.deepEqual([...files.values()].sort(), [...originalFiles.values()].sort());
    for (const binding of added.materialManifests) {
      assert.equal(binding.materialState, undefined, 'retain source claims/model/species envelopes without guessing identity');
      const manifest = files.get(binding.manifestPath);
      assert.ok(manifest.startsWith('frankensim.matdb-manifest.'));
      for (const row of manifest.split(/\r?\n/).filter(line => line.startsWith('source\t'))) {
        const path = posix.join(posix.dirname(binding.manifestPath), row.split('\t')[2]);
        assert.ok(files.has(path), `${entry.id}: original compiler input ${path}`);
      }
    }
    const species = [...files.entries()].find(([path]) => path.endsWith('/species.tsv'))[1];
    assert.equal(species.split(/\r?\n/).find(row => row.startsWith('species\t')).split('\t')[1], entry.species);
    assert.deepEqual(Object.values(added.component), added.materialManifests.map(binding => binding.packPath));
    total += added.attachments.reduce((sum, file) => sum + file.bytes.length, 0);
  }
  assert.ok(total < MAX_ATTACHMENT_BYTES, 'all eight original gas sources fit the existing aggregate byte limit');
});

test('multigas material-string experiments preserve their mechanics and reach the native compiler/render request', async () => {
  const starters = JSON.parse(await readOriginal('apps/instrument-lab/src/native-string-starters.json'));
  const input = starters['material_strings/pluck-and-bow.material-strings'];
  const original = await readOriginal('crates/fs-couple/examples/material_strings/pluck-and-bow.material-strings');
  const prepared = await prepareNativeTextExperiment(original, input.sources,
    {renderer:'strings', materialManifests:input.materialManifests}, readOriginal);
  let merged = {attachments:prepared.attachments, materialManifests:prepared.materialManifests};
  const components = [];
  for (const id of ['nitrogen', 'oxygen', 'helium', 'sulfur-hexafluoride']) {
    const added = await loadNativeGasSource(id, readOriginal);
    merged = mergeNativeGasSource(merged.attachments, merged.materialManifests, added);
    components.push(added.component);
  }
  assert.ok(merged.attachments.length + merged.materialManifests.length > 32);
  assert.ok(merged.attachments.length + merged.materialManifests.length <= MAX_ATTACHMENTS);
  const [region] = nativeGasRegions(original);
  const source = editNativeGasRegion(original, region, {kind:'thermochemical', temperature:region.values.temperature,
    pressure:region.values.pressure, basis:'mole', components:[{species:'N2', fraction:'0.6'},
      {species:'O2', fraction:'0.2'}, {species:'He', fraction:'0.1'}, {species:'SF6', fraction:'0.1'}]});
  assert.equal(source.slice(source.indexOf('\nbody ')), original.slice(original.indexOf('\nbody ')),
    'the specimen, bow, tuning and pressure scale are not changed by a gas source choice');
  const gasSources = {components, waterSpecies:null, phase:'frozen'};
  const options = {renderer:'strings', gasSources, materialManifests:merged.materialManifests};
  const saved = serializeNativeExperiment(source, merged.attachments, options);
  const restored = validateInput(encode(saved), true);
  assert.deepEqual(restored.gasSources, gasSources);
  assert.deepEqual(restored.materialManifests, merged.materialManifests);
  const status = {available:true, renderers:{wind:true, strings:true}, token:'test-session',
    maximumSourceBytes:MAX_SOURCE_BYTES, maximumOutputBytes:16 * 1024 * 1024,
    maximumAttachments:MAX_ATTACHMENTS, maximumAttachmentFileBytes:MAX_ATTACHMENT_FILE_BYTES,
    maximumAttachmentBytes:MAX_ATTACHMENT_BYTES, maximumRequestBytes:MAX_REQUEST_BYTES,
    supportsMaterialManifests:true, supportsGasSources:true};
  let submitted;
  await assert.rejects(renderNativePerformance(source, merged.attachments, status, options, undefined,
    async (url, request) => {
      assert.equal(url, '/api/native-wind/render');
      submitted = serviceInput(encode(request.body), true);
      return new Response(JSON.stringify({error:'native-boundary witness'}), {
        status:400, headers:{'Content-Type':'application/json'},
      });
    }), /native-boundary witness/);
  assert.deepEqual(submitted, restored);
});

test('gas source additions retain edited/uploaded sources and refuse conflicting paths before publishing', async () => {
  const added = await loadNativeGasSource('nitrogen', readOriginal);
  const empty = {attachments:[], materialManifests:[]};
  const first = mergeNativeGasSource(empty.attachments, empty.materialManifests, added);
  assert.deepEqual(empty, {attachments:[], materialManifests:[]});
  const repeated = mergeNativeGasSource(first.attachments, first.materialManifests, added);
  assert.deepEqual(repeated, first, 'reselecting the identical original source is idempotent');
  const edited = first.attachments.map((file, index) => index === 0
    ? {...file, bytes:encode(decode(file.bytes) + '# user source edit\n')} : file);
  const before = structuredClone(edited);
  assert.throws(() => mergeNativeGasSource(edited, first.materialManifests, added), /different data/);
  assert.deepEqual(edited, before);
  for (const attachments of [
    [{path:added.component.speciesPath, bytes:encode('original user pack')}],
    [{path:'gas-sources', bytes:encode('original file, not a directory')}],
    [{path:'../source', bytes:encode('invalid source')}],
  ]) assert.throws(() => mergeNativeGasSource(attachments, [], added), /already supplied|collide|traversal/);
  const named = [{...first.materialManifests[0], materialState:{chemistry:'other', phase:'gas', process:'specified', revision:0}}];
  assert.throws(() => mergeNativeGasSource([], named, added), /another source/);
  assert.throws(() => mergeNativeGasSource([], [], added, 8), /8 file allowance/);
  await assert.rejects(loadNativeGasSource('unsupported', readOriginal), /complete/);
  await assert.rejects(loadNativeGasSource('nitrogen', async path => {
    if (path.endsWith('/properties.tsv')) throw new Error('source chunk unavailable');
    return readOriginal(path);
  }), /source chunk unavailable/);
});
