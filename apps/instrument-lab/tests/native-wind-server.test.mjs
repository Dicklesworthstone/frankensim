// Actual HTTP and subprocess tests with an explicitly synthetic executable.
// This proves the adapter, not Rust execution, acoustics, or instrument realism.
import {test} from 'node:test';
import assert from 'node:assert/strict';
import {mkdtemp, mkdir, readFile, readdir, rm, symlink, writeFile} from 'node:fs/promises';
import {tmpdir} from 'node:os';
import {join, isAbsolute} from 'node:path';
import http from 'node:http';
import {nativeGasRegions, editNativeGasRegion} from '../src/native-gas-records.mjs';
import {serializeNativeExperiment} from '../src/native-experiment.mjs';
import {createNativeWindServer} from '../tools/serve-native-wind.mjs';
import {MAX_SOURCE_BYTES, MAX_ATTACHMENTS, MAX_ATTACHMENT_FILE_BYTES, validateInput, validateSource} from '../tools/native-wind-protocol.mjs';

const source = (mode = 'ok') => Buffer.from(`frankensim-reed-performance-v1\r\n${mode}\r\n# UTF-8: 音; $() '--help'\r\n`);
const bundle = (bytes = source(), files = []) => Buffer.from(JSON.stringify({
  format:'frankensim-native-wind-input-v1', source:bytes.toString('base64'),
  attachments:files.map(([path, bytes]) => ({path, bytes:bytes.toString('base64')})),
}));
const pause = ms => new Promise(resolve => setTimeout(resolve, ms));
async function until(fn) { for (let i = 0; i < 400; i++) { if (await fn()) return; await pause(10); } throw new Error('Timed out waiting for test condition.'); }
async function setup(options = {}) {
  const {withStrings = false, withFiniteGasString = false, withMaterialCompiler = false, ...serverOptions} = options;
  const dir = await mkdtemp(join(tmpdir(), 'native-wind-test-')), dist = join(dir, 'dist'), tempRoot = join(dir, 'jobs');
  await mkdir(dist); await mkdir(tempRoot); await writeFile(join(dist, 'index.html'), '<title>Synthetic adapter test</title>');
  const receipt = join(dir, 'receipt.json'), engine = join(dir, 'fixture.mjs');
  const fixture = `#!${process.execPath}
// Synthetic process fixture: never presented as native physics.
import {readFileSync, writeFileSync, statSync} from 'node:fs';
const args = process.argv.slice(2), offset = ['wind', 'resolved'].includes(args[0]) ? 1 : 0;
const finite = args.includes('--history-output');
const inputPath = args[offset], outputPath = args[finite ? 2 : offset + 1], input = readFileSync(inputPath);
const text = input.toString(), requested = /^# read-files (.+)$/m.exec(text);
const files = requested ? JSON.parse(requested[1]).map(path => ({
 path, bytes:readFileSync(path).toString('base64'), mode:statSync(path).mode & 0o777,
})) : [];
writeFileSync(${JSON.stringify(receipt)}, JSON.stringify({args, cwd:process.cwd(), bytes:input.toString('base64'), files, pid:process.pid,
 directoryMode:statSync(process.cwd()).mode & 0o777, sourceMode:statSync(inputPath).mode & 0o777}));
if (text.includes('test-refuse')) { console.error('synthetic parser refusal'); process.exit(2); }
if (text.includes('test-log')) { process.stderr.write(Buffer.alloc(300000, 120)); await new Promise(() => setInterval(() => {}, 1000)); }
if (text.includes('test-hang')) { process.on('SIGTERM', () => {}); await new Promise(() => setInterval(() => {}, 1000)); }
if (text.includes('test-invalid')) { writeFileSync(outputPath, 'not a wave'); process.exit(0); }
if (text.includes('test-large')) { writeFileSync(outputPath, Buffer.alloc(4097)); process.exit(0); }
if (text.includes('test-delay')) await new Promise(resolve => setTimeout(resolve, 200));
const wav = Buffer.alloc(48); wav.write('RIFF'); wav.writeUInt32LE(40, 4); wav.write('WAVEfmt ', 8);
wav.writeUInt32LE(16, 16); wav.writeUInt16LE(1, 20); wav.writeUInt16LE(1, 22); wav.writeUInt32LE(48000, 24);
wav.writeUInt32LE(96000, 28); wav.writeUInt16LE(2, 32); wav.writeUInt16LE(16, 34); wav.write('data', 36); wav.writeUInt32LE(4, 40);
if (text.includes('test-wrong-rate')) { wav.writeUInt32LE(16000, 24); wav.writeUInt32LE(32000, 28); }
writeFileSync(outputPath, wav);
if (finite) {
 const history = text.includes('test-history-large') ? Buffer.alloc(257, 49)
  : text.includes('test-history-invalid') ? Buffer.from([255])
  : Buffer.from('time_s,pressure_pa\\r\\n0.000020833333333333333,0\\r\\n');
 if (!text.includes('test-history-missing')) writeFileSync(args[args.indexOf('--history-output') + 1], history);
 if (!text.includes('test-no-calibration')) console.error('samples=2; rate_Hz=48000; clips=0; full_scale_Pa=1e-2');
}
`;
  await writeFile(engine, fixture, {mode:0o700});
  const compiler = join(dir, 'material-compiler.mjs'), compilerReceipt = join(dir, 'compiler.json');
  if (withMaterialCompiler) await writeFile(compiler, `#!${process.execPath}
// Synthetic compiler transport witness; no material normalization or physics.
import {readFileSync, writeFileSync, truncateSync} from 'node:fs';
import {dirname, join} from 'node:path';
const args = process.argv.slice(2);
if (args[0] !== 'matdb-pack') process.exit(3);
const manifest = args[args.indexOf('--manifest') + 1], output = args[args.indexOf('--out') + 1];
const source = readFileSync(manifest), table = readFileSync(join(dirname(manifest), 'properties.tsv'));
writeFileSync(${JSON.stringify(compilerReceipt)}, JSON.stringify({args, pid:process.pid,
 manifest:source.toString('base64'), table:table.toString('base64')}));
if (source.toString().includes('compile-refuse')) { console.error('synthetic material refusal'); process.exit(2); }
if (source.toString().includes('compile-hang')) {
 process.on('SIGTERM', () => {}); await new Promise(() => setInterval(() => {}, 1000));
}
writeFileSync(output, Buffer.concat([Buffer.from('compiled transport witness\\0'), source, table]), {flag:'wx', mode:0o600});
if (source.toString().includes('compile-large')) truncateSync(output, ${MAX_ATTACHMENT_FILE_BYTES + 1});
`, {mode:0o700});
  const app = createNativeWindServer({engine, materialCompiler:withMaterialCompiler ? compiler : null, stringsEngine:withStrings ? engine : null,
    finiteGasStringEngine:withFiniteGasString ? engine : null, dist, tempRoot, ...serverOptions});
  await new Promise((accept, reject) => {app.server.once('error', reject); app.server.listen(0, '127.0.0.1', accept);});
  const url = `http://127.0.0.1:${app.server.address().port}`;
  const status = () => fetch(url + '/api/native-wind/status').then(r => r.json());
  const initial = await status();
  const run = (bytes = source(), extra = {}) => fetch(url + '/api/native-wind/render', {method:'POST',
    headers:{'Content-Type':'text/plain;charset=utf-8', 'X-Native-Wind-Token':initial.token}, body:bytes, ...extra});
  const runBundle = bytes => run(bytes, {headers:{'Content-Type':'application/json', 'X-Native-Wind-Token':initial.token}});
  return {...app, url, dir, receipt, compilerReceipt, tempRoot, status, initial, run, runBundle,
    async cleanup() { await app.shutdown(); await rm(dir, {recursive:true, force:true}); }};
}
const native = {skip:process.platform === 'win32'};

test('material strings reuse the exact source packet, gas associations and explicit profile paths', native, async () => {
  const app = await setup({withStrings:true, engine:null});
  try {
    assert.deepEqual(app.initial.renderers, {wind:false, strings:true, resolved:false, 'finite-gas-string':false});
    assert.equal(app.initial.supportsCoatingStrainGuard, true);
    assert.equal(app.initial.available, true);
    assert.ok(app.initial.supportedSchemas.includes('frankensim-material-strings-v3'));
    const names = ['wire.pack', 'body.pack', 'shape.body', 'finger.stops', 'contact.thermal', 'gas/s', 'gas/n', 'gas/t'];
    const bytes = Buffer.from('frankensim-material-strings-v3\r\n# read-files ' + JSON.stringify(names) + '\r\n# original 音\r\n');
    const files = names.map((name, i) => [name, Buffer.from([0, 255, i, 13, 10])]);
    const packet = JSON.parse(bundle(bytes, files));
    packet.renderer = 'strings'; packet.decimate = true;
    packet.stringProfiles = {bodyPath:'shape.body', stopPath:'finger.stops', contactThermalPath:'contact.thermal'};
    packet.coatingStrainGuard = {maximumBytes:536870912, maximumWork:1000000000000};
    packet.gasSources = {components:[{speciesPath:'gas/s', nasa9Path:'gas/n', transportPath:'gas/t'}], waterSpecies:null, phase:'frozen'};
    const response = await app.runBundle(JSON.stringify(packet));
    assert.equal(response.status, 200); assert.equal((await response.arrayBuffer()).byteLength, 48);
    const receipt = JSON.parse(await readFile(app.receipt, 'utf8'));
    assert.equal(receipt.bytes, bytes.toString('base64'));
    assert.deepEqual(receipt.args, [join(receipt.cwd, 'input.performance'), join(receipt.cwd, 'output.wav'),
      '--self-contained', '--decimate', ...names.flatMap(name => ['--attached-source', name]),
      '--gas-water', 'none', '--gas-phase', 'frozen', '--gas-component', 'gas/s', 'gas/n', 'gas/t',
      '--body-profile', 'shape.body', '--stop-profile', 'finger.stops', '--contact-thermal', 'contact.thermal',
      '--coating-strain-guard', '536870912', '1000000000000']);
    assert.deepEqual(receipt.files.map(file => [file.path, file.bytes]), files.map(([path, bytes]) => [path, bytes.toString('base64')]));
    await until(async () => !(await app.status()).busy);
    assert.equal((await app.run()).status, 503, 'a configured string performer must not substitute for a missing wind performer');
  } finally { await app.cleanup(); }
});

test('resolved gas fields and material sources run through their own native command unchanged', native, async () => {
  const app = await setup({timeoutMs:120000});
  try {
    assert.equal(app.initial.renderers.resolved, true);
    assert.equal(app.initial.timeoutMs, 120000);
    assert.ok(app.initial.supportedSchemas.includes('frankensim-resolved-dg-performance-v1'));
    const names = ['fields/initial.txt', 'wall.fsmatpk', 'gas/species', 'gas/nasa9', 'gas/transport'];
    const bytes = Buffer.from('frankensim-resolved-dg-performance-v1\r\n# read-files '
      + JSON.stringify(names) + '\r\n# original physical states µ\r\n');
    const files = names.map((name, i) => [name, Buffer.from([0, 255, i, 13, 10])]);
    const packet = JSON.parse(bundle(bytes, files));
    packet.renderer = 'resolved'; packet.decimate = true;
    packet.gasSources = {components:[{speciesPath:'gas/species', nasa9Path:'gas/nasa9',
      transportPath:'gas/transport'}], waterSpecies:null, phase:'frozen'};
    const response = await app.runBundle(JSON.stringify(packet));
    assert.equal(response.status, 200); assert.equal((await response.arrayBuffer()).byteLength, 48);
    const record = JSON.parse(await readFile(app.receipt, 'utf8'));
    assert.equal(record.bytes, bytes.toString('base64'));
    assert.deepEqual(record.args, ['resolved', join(record.cwd, 'input.performance'), join(record.cwd, 'output.wav'),
      '--self-contained', '--decimate', ...names.flatMap(name => ['--attached-source', name]),
      '--gas-water', 'none', '--gas-phase', 'frozen', '--gas-component', 'gas/species', 'gas/nasa9', 'gas/transport']);
    assert.deepEqual(record.files.map(file => [file.path, file.bytes]), files.map(([path, bytes]) => [path, bytes.toString('base64')]));
  } finally { await app.cleanup(); }
});


test('edited resolved named mixtures reach the native service as the exact portable source', native, async () => {
  const app = await setup();
  try {
    const original = await readFile(new URL('../../../crates/fs-couple/examples/resolved-dg-radiating-mouth.performance', import.meta.url), 'utf8');
    const source = original.replace('fluid_limits 144 250 600', 'fluid_limits 144 280 320');
    const primary = nativeGasRegions(source).find(region => region.key === 'gas');
    const values = {kind:'mixture', temperature:'310', pressure:'90000', basis:'mass',
      components:[{species:'N2', fraction:'0.7'}, {species:'O2', fraction:'0.2'}, {species:'He', fraction:'0.1'}]};
    const changed = editNativeGasRegion(source, primary, values, true);
    assert.equal(changed, source.replace('gas dry-air 300 100000',
      'gas mixture 310 90000 mass 3 N2 0.7 O2 0.2 He 0.1')
      .replace('dg_exterior 300 100000 0', 'dg_exterior 300 90000 0'));
    const response = await app.runBundle(serializeNativeExperiment(changed, [], {renderer:'resolved', decimate:true}));
    assert.equal(response.status, 200);
    assert.equal((await response.arrayBuffer()).byteLength, 48);
    const record = JSON.parse(await readFile(app.receipt, 'utf8'));
    assert.equal(record.bytes, Buffer.from(changed).toString('base64'));
    assert.deepEqual(record.args, ['resolved', join(record.cwd, 'input.performance'),
      join(record.cwd, 'output.wav'), '--self-contained', '--decimate']);
    assert.deepEqual(record.files, []);
  } finally { await app.cleanup(); }
});

test('renderer and profile declarations refuse mismatched schemas, unlisted paths and arbitrary flags', () => {
  const bytes = Buffer.from('frankensim-material-strings-v1\nopaque');
  for (const renderer of ['strings', 'wind', 'resolved', 'shell', null]) {
    const packet = {...JSON.parse(bundle(bytes, [['body.profile', Buffer.from('original')]])), renderer};
    if (renderer === 'strings') {
      assert.equal(validateInput(Buffer.from(JSON.stringify(packet)), true).renderer, 'strings');
      for (const stringProfiles of [{bodyPath:'missing'}, {stopPath:'../body.profile'}, {observerStore:'body.profile'}, {bodyPath:null}])
        assert.throws(() => validateInput(Buffer.from(JSON.stringify({...packet, stringProfiles})), true), /String profiles/);
    } else assert.throws(() => validateInput(Buffer.from(JSON.stringify(packet)), true));
  }
  const wind = {...JSON.parse(bundle(source(), [['body.profile', Buffer.from('original')]])), stringProfiles:{bodyPath:'body.profile'}};
  assert.throws(() => validateInput(Buffer.from(JSON.stringify(wind)), true), /String profiles/);
  const resolved = {...JSON.parse(bundle(Buffer.from('frankensim-resolved-dg-performance-v1\n'))), renderer:'resolved'};
  assert.equal(validateInput(Buffer.from(JSON.stringify(resolved)), true).renderer, 'resolved');
  for (const renderer of ['wind', 'strings'])
    assert.throws(() => validateInput(Buffer.from(JSON.stringify({...resolved, renderer})), true));
  assert.throws(() => validateInput(Buffer.from(JSON.stringify({...resolved, stringProfiles:{}})), true), /String profiles/);
});

test('string source packets cannot fall back to the wind executable', native, async () => {
  const app = await setup();
  try {
    const packet = {...JSON.parse(bundle(Buffer.from('frankensim-material-strings-v3\nopaque'))), renderer:'strings'};
    const response = await app.runBundle(JSON.stringify(packet));
    assert.equal(response.status, 503);
    assert.match((await response.json()).error, /material-string executable/);
    await assert.rejects(readFile(app.receipt), {code:'ENOENT'});
  } finally { await app.cleanup(); }
});

test('transport retains exact UTF-8 and admits the original native wind source families', () => {
  const bytes = source(); assert.equal(validateSource(bytes), bytes);
  const flute = Buffer.from('frankensim-flue-dg-performance-v1\nopaque'); assert.equal(validateSource(flute), flute);
  const nolle = Buffer.from('frankensim-flue-performance-v1\nopaque'); assert.equal(validateSource(nolle), nolle);
  for (const bad of [Buffer.alloc(0), Buffer.from([0xff]), Buffer.from('frankensim-material-strings-v1\nunsupported'), Buffer.alloc(MAX_SOURCE_BYTES + 1)])
    assert.throws(() => validateSource(bad));
  assert.throws(() => createNativeWindServer({engine:'relative-executable'}), /absolute/);
});

test('loopback status/static service, token, origin and confined assets', native, async () => {
  const app = await setup({engine:'/nonexistent/music_render'});
  try {
    assert.equal(app.initial.available, false); assert.equal(app.initial.busy, false);
    assert.equal(app.initial.maximumSourceBytes, MAX_SOURCE_BYTES);
    assert.equal((await app.run()).status, 503);
    assert.match(await (await fetch(app.url)).text(), /Synthetic adapter test/);
    await writeFile(join(app.dir, 'private.js'), 'private'); await symlink(join(app.dir, 'private.js'), join(app.dir, 'dist', 'leak.js'));
    assert.equal((await fetch(app.url + '/leak.js')).status, 404);
    assert.equal((await fetch(app.url + '/api/native-wind/status', {headers:{Origin:'https://example.com'}})).status, 403);
    const forbiddenHost = await new Promise(resolve => http.get(app.url + '/api/native-wind/status',
      {headers:{Host:'attacker.example'}}, res => {res.resume(); resolve(res.statusCode);}));
    assert.equal(forbiddenHost, 403);
    assert.equal((await app.run(source(), {headers:{'Content-Type':'text/plain'}})).status, 403);
    assert.equal((await app.run(source(), {headers:{'Content-Type':'application/octet-stream', 'X-Native-Wind-Token':app.initial.token}})).status, 415);
    assert.equal((await fetch(app.url + '/api/execute')).status, 404);
  } finally { await app.cleanup(); }
});

test('built MIDI and opaque binary assets retain their bytes without opening paths outside dist', native, async () => {
  const app = await setup({engine:'/nonexistent/music_render'});
  try {
    const bytes = Buffer.from([77, 84, 104, 100, 0, 255, 128, 1]);
    for (const [name, type] of [['phrase.mid', 'audio/midi'], ['phrase.midi', 'audio/midi'], ['body.fsb', 'application/octet-stream']]) {
      await writeFile(join(app.dir, 'dist', name), bytes);
      const response = await fetch(app.url + '/' + name);
      assert.equal(response.status, 200); assert.equal(response.headers.get('content-type'), type);
      assert.equal(response.headers.get('content-length'), String(bytes.length)); assert.equal(response.headers.get('accept-ranges'), 'bytes');
      assert.deepEqual(Buffer.from(await response.arrayBuffer()), bytes);
    }
    const head = await fetch(app.url + '/phrase.mid', {method:'HEAD', headers:{Range:'bytes=1-2'}});
    assert.equal(head.status, 200); assert.equal(head.headers.get('content-type'), 'audio/midi');
    assert.equal(head.headers.get('content-length'), String(bytes.length)); assert.equal((await head.arrayBuffer()).byteLength, 0);
    await writeFile(join(app.dir, 'private.fsb'), bytes);
    await symlink(join(app.dir, 'private.fsb'), join(app.dir, 'dist', 'leak.fsb'));
    await writeFile(join(app.dir, 'dist', '.hidden.fsb'), bytes);
    for (const path of ['/leak.fsb', '/.hidden.fsb', '/%2e%2e%2fprivate.fsb'])
      assert.equal((await fetch(app.url + path)).status, 404);
  } finally { await app.cleanup(); }
});

test('static audio serves single byte ranges, suffixes, full fallbacks and unsatisfiable 416 responses', native, async () => {
  const app = await setup({engine:'/nonexistent/music_render'});
  try {
    const bytes = Buffer.from('0123456789abcdef');
    await writeFile(join(app.dir, 'dist', 'clip.wav'), bytes);
    for (const [range, start, end] of [['bytes=2-5', 2, 5], ['bytes=10-', 10, 15],
      ['bytes=-3', 13, 15], ['bytes=14-99', 14, 15], ['bytes=-99', 0, 15]]) {
      const response = await fetch(app.url + '/clip.wav', {headers:{Range:range}});
      assert.equal(response.status, 206); assert.equal(response.headers.get('content-type'), 'audio/wav');
      assert.equal(response.headers.get('content-range'), `bytes ${start}-${end}/${bytes.length}`);
      assert.equal(response.headers.get('content-length'), String(end - start + 1));
      assert.deepEqual(Buffer.from(await response.arrayBuffer()), bytes.subarray(start, end + 1));
    }
    for (const range of ['bytes=16-', 'bytes=-0', 'bytes=9-2']) {
      const response = await fetch(app.url + '/clip.wav', {headers:{Range:range}});
      assert.equal(response.status, 416); assert.equal(response.headers.get('content-range'), 'bytes */16');
      assert.equal(response.headers.get('content-length'), '0'); assert.equal((await response.arrayBuffer()).byteLength, 0);
    }
    for (const headers of [{Range:'bytes=0-1,4-5'}, {Range:'items=1-3'}, {Range:'bytes=0-1', 'If-Range':'"unknown"'}]) {
      const response = await fetch(app.url + '/clip.wav', {headers}); assert.equal(response.status, 200);
      assert.equal(response.headers.get('content-length'), '16'); assert.deepEqual(Buffer.from(await response.arrayBuffer()), bytes);
    }
    await writeFile(join(app.dir, 'dist', 'empty.bin'), Buffer.alloc(0));
    const empty = await fetch(app.url + '/empty.bin'); assert.equal(empty.status, 200); assert.equal(empty.headers.get('content-length'), '0');
    const range = await fetch(app.url + '/empty.bin', {headers:{Range:'bytes=0-'}});
    assert.equal(range.status, 416); assert.equal(range.headers.get('content-range'), 'bytes */0');
  } finally { await app.cleanup(); }
});

test('exact native arguments and bytes, private cwd, real returned WAV and cleanup', native, async () => {
  const app = await setup();
  try {
    const bytes = source(), response = await app.run(bytes);
    assert.equal(response.status, 200); assert.match(response.headers.get('content-type'), /^audio\/wav/);
    const wave = Buffer.from(await response.arrayBuffer()); assert.equal(wave.toString('ascii', 0, 4), 'RIFF'); assert.equal(wave.length, 48);
    const receipt = JSON.parse(await readFile(app.receipt, 'utf8'));
    assert.equal(receipt.bytes, bytes.toString('base64')); assert.equal(receipt.args.length, 4); assert.equal(receipt.args[0], 'wind');
    assert.equal(receipt.args[1], join(receipt.cwd, 'input.performance')); assert.equal(receipt.args[2], join(receipt.cwd, 'output.wav'));
    assert.equal(receipt.args[3], '--self-contained');
    assert.ok(isAbsolute(receipt.cwd)); assert.ok(receipt.cwd.startsWith(app.tempRoot + '/'));
    assert.equal(receipt.directoryMode, 0o700); assert.equal(receipt.sourceMode, 0o600);
    await until(async () => !(await app.status()).busy); assert.deepEqual(await readdir(app.tempRoot), []);
  } finally { await app.cleanup(); }
});

test('one active render; native refusal, malformed WAV and output/log limits never publish audio', native, async () => {
  const app = await setup({maximumOutputBytes:1024});
  try {
    const first = app.run(source('test-delay'));
    await until(async () => (await app.status()).busy);
    assert.equal((await app.run()).status, 409); assert.equal((await first).status, 200);
    await until(async () => !(await app.status()).busy);
    for (const [mode, status, message] of [['test-refuse', 400, /synthetic parser refusal/],
      ['test-invalid', 502, /WAV/], ['test-large', 502, /output limit/], ['test-log', 502, /diagnostics/]]) {
      const response = await app.run(source(mode)); assert.equal(response.status, status);
      assert.match(response.headers.get('content-type'), /application\/json/); assert.match((await response.json()).error, message);
      assert.equal((await app.status()).busy, false); assert.deepEqual(await readdir(app.tempRoot), []);
    }
  } finally { await app.cleanup(); }
});

test('bounded source uploads including chunked bodies', native, async () => {
  const app = await setup();
  try {
    assert.equal((await app.run(Buffer.alloc(MAX_SOURCE_BYTES + 1))).status, 413);
    const result = await new Promise((accept, reject) => {
      const req = http.request(app.url + '/api/native-wind/render', {method:'POST', headers:{
        'Content-Type':'text/plain', 'X-Native-Wind-Token':app.initial.token, 'Transfer-Encoding':'chunked'}}, res => {
          res.resume(); res.once('end', () => accept(res.statusCode)); });
      req.on('error', reject); req.write(Buffer.alloc(MAX_SOURCE_BYTES)); req.end(Buffer.from('x'));
    });
    assert.equal(result, 413); await until(async () => !(await app.status()).busy);
    assert.deepEqual(await readdir(app.tempRoot), []);
  } finally { await app.cleanup(); }
});

test('timeout, browser abort and shutdown kill the child and remove its private files', native, async () => {
  const timed = await setup({timeoutMs:150});
  try { const result = await timed.run(source('test-hang')); assert.equal(result.status, 504);
    assert.match((await result.json()).error, /execution budget/); assert.equal((await timed.status()).busy, false);
    const record = JSON.parse(await readFile(timed.receipt, 'utf8')); assert.throws(() => process.kill(record.pid, 0));
    assert.deepEqual(await readdir(timed.tempRoot), []);
  } finally { await timed.cleanup(); }
  for (const shutdown of [false, true]) {
    const app = await setup(), controller = new AbortController();
    try {
      const pending = app.run(source('test-hang'), {signal:controller.signal}).catch(e => e);
      await until(async () => { try { await readFile(app.receipt); return true; } catch { return false; } });
      const record = JSON.parse(await readFile(app.receipt, 'utf8'));
      if (shutdown) await app.shutdown(); else controller.abort();
      await pending;
      if (!shutdown) await until(async () => !(await app.status()).busy);
      assert.throws(() => process.kill(record.pid, 0)); assert.deepEqual(await readdir(app.tempRoot), []);
    } finally { await app.cleanup(); }
  }
});

test('material valve source and nested binary packs reach the native process without translation', native, async () => {
  const app = await setup();
  try {
    const paths = ['materials/reed.fsmcdpk', 'materials/wall.fsmcdpk'];
    const bytes = Buffer.from('frankensim-material-modal-dg-valve-v1\n# read-files ' + JSON.stringify(paths) + '\n# µ source\n');
    const first = Buffer.from([0, 255, 128, 13, 10, 39, 36, 40, 41]);
    const second = Buffer.from('independent source material bytes');
    assert.equal(app.initial.maximumAttachments, MAX_ATTACHMENTS);
    assert.equal(app.initial.maximumAttachmentFileBytes, MAX_ATTACHMENT_FILE_BYTES);
    const response = await app.runBundle(bundle(bytes, [[paths[0], first], [paths[1], second]]));
    assert.equal(response.status, 200);
    assert.equal((await response.arrayBuffer()).byteLength, 48);
    const record = JSON.parse(await readFile(app.receipt, 'utf8'));
    assert.equal(record.bytes, bytes.toString('base64'));
    assert.deepEqual(record.files, [
      {path:paths[0], bytes:first.toString('base64'), mode:0o600},
      {path:paths[1], bytes:second.toString('base64'), mode:0o600},
    ]);
    assert.equal(record.args[0], 'wind');
    assert.deepEqual(record.args.slice(3), ['--self-contained', '--attached-source', paths[0], '--attached-source', paths[1]]);
    await until(async () => !(await app.status()).busy);
    assert.deepEqual(await readdir(app.tempRoot), []);
    // A material swap is only a source-byte change; no CLI flags or PCM gain change.
    const replacement = Buffer.from([91, 4, 3, 2, 1, 0]);
    assert.equal((await app.runBundle(bundle(bytes, [[paths[0], replacement], [paths[1], second]]))).status, 200);
    const swapped = JSON.parse(await readFile(app.receipt, 'utf8'));
    assert.equal(swapped.files[0].bytes, replacement.toString('base64'));
    assert.equal(swapped.bytes, record.bytes); assert.equal(swapped.args.length, 8);
  } finally { await app.cleanup(); }
});

test('attachment envelopes refuse traversal, collisions, malformed bytes and excess files before execution', native, async () => {
  const app = await setup();
  try {
    const one = Buffer.from([1,2,3]);
    for (const path of ['../outside', '/absolute', 'a/../b', 'a//b', 'a\\b', 'a:b',
      '.hidden', 'a/.hidden', 'a b', '-flag', 'input.performance', 'output.wav', 'output.provenance.json', 'input.performance/child']) {
      const response = await app.runBundle(bundle(source(), [[path, one]]));
      assert.equal(response.status, 400, path);
    }
    for (const paths of [['same', 'same'], ['parent', 'parent/child'], ['parent/child', 'parent']]) {
      assert.equal((await app.runBundle(bundle(source(), paths.map(path => [path, one])))).status, 400);
    }
    const malformed = JSON.parse(bundle().toString());
    for (const bytes of ['AA=A', 'AB==', 'abc', '!!!!', '']) {
      malformed.attachments = [{path:'wall.fsmcdpk', bytes}];
      assert.equal((await app.runBundle(Buffer.from(JSON.stringify(malformed)))).status, 400);
    }
    assert.equal((await app.runBundle(bundle(source(), Array.from({length:MAX_ATTACHMENTS+1}, (_, i) => [String(i), one])))).status, 400);
    await assert.rejects(readFile(app.receipt), {code:'ENOENT'});
    assert.deepEqual(await readdir(app.tempRoot), []);
  } finally { await app.cleanup(); }
});

test('native missing-file refusal is visible and failed source jobs keep no attachments', native, async () => {
  const app = await setup();
  try {
    const bytes = Buffer.from('frankensim-material-modal-dg-valve-v1\n# read-files ["missing.fsmcdpk"]\n');
    const response = await app.runBundle(bundle(bytes, [['different.fsmcdpk', Buffer.from([1])]]));
    assert.equal(response.status, 400); assert.match((await response.json()).error, /Native renderer refused/);
    assert.deepEqual(await readdir(app.tempRoot), []);
  } finally { await app.cleanup(); }
});

test('source and per-file byte budgets apply to the decoded attachment envelope', () => {
  const packet = bundle(source(), [['wall.fsmcdpk', Buffer.alloc(MAX_ATTACHMENT_FILE_BYTES + 1, 7)]]);
  assert.throws(() => validateInput(packet, true), /byte limit|base64/);
  assert.throws(() => validateInput(bundle(Buffer.alloc(MAX_SOURCE_BYTES + 1, 97)), true), /byte limit|base64/);
  assert.throws(() => validateInput(Buffer.from([0xff]), true), /UTF-8 JSON/);
});

test('native decimation is explicit while original source clocks and bytes stay unchanged', native, async () => {
  const app = await setup();
  try {
    const bytes = Buffer.from('frankensim-flue-dg-performance-v1\nconfig 96000 96 37 0.2 1e-5\n');
    const packet = JSON.parse(bundle(bytes).toString());
    packet.decimate = true;
    const response = await app.runBundle(Buffer.from(JSON.stringify(packet)));
    assert.equal(response.status, 200);
    assert.equal((await response.arrayBuffer()).byteLength, 48);
    const record = JSON.parse(await readFile(app.receipt, 'utf8'));
    assert.equal(record.bytes, bytes.toString('base64'));
    assert.deepEqual(record.args.slice(3), ['--self-contained', '--decimate']);
    for (const value of [1, 'true', null]) {
      packet.decimate = value;
      assert.throws(() => validateInput(Buffer.from(JSON.stringify(packet)), true), /decimation/);
    }
    delete packet.decimate; packet.gain = 2;
    assert.throws(() => validateInput(Buffer.from(JSON.stringify(packet)), true), /unsupported fields/);
  } finally { await app.cleanup(); }
});

test('gas declarations bind only explicitly attached original pack bytes and fixed native arguments', native, async () => {
  const app = await setup();
  try {
    assert.equal(app.initial.supportsGasSources, true);
    const paths = ['gas/custom.species', 'gas/custom.nasa9', 'gas/custom.transport'];
    const bytes = Buffer.from('frankensim-flue-performance-v1\n# read-files ' + JSON.stringify(paths) + '\n');
    const originals = paths.map((path, index) => [path, Buffer.from([0, 255, index, 13, 10])]);
    const packet = JSON.parse(bundle(bytes, originals).toString());
    packet.gasSources = {components:[{speciesPath:paths[0], nasa9Path:paths[1], transportPath:paths[2]}],
      waterSpecies:null, phase:'frozen'};
    const response = await app.runBundle(Buffer.from(JSON.stringify(packet)));
    assert.equal(response.status, 200); assert.equal((await response.arrayBuffer()).byteLength, 48);
    const record = JSON.parse(await readFile(app.receipt, 'utf8'));
    assert.equal(record.bytes, bytes.toString('base64'));
    assert.deepEqual(record.files, originals.map(([path, bytes]) => ({path, bytes:bytes.toString('base64'), mode:0o600})));
    assert.deepEqual(record.args.slice(3), ['--self-contained',
      ...paths.flatMap(path => ['--attached-source', path]),
      '--gas-water', 'none', '--gas-phase', 'frozen', '--gas-component', ...paths]);
    await until(async () => !(await app.status()).busy);
    assert.deepEqual(await readdir(app.tempRoot), []);
  } finally { await app.cleanup(); }
});

test('gas associations refuse missing source paths and absent phase or water policy before execution', () => {
  const paths = ['gas/s', 'gas/n', 'gas/t'];
  const base = JSON.parse(bundle(source(), paths.map(path => [path, Buffer.from([1])])).toString());
  const component = {speciesPath:paths[0], nasa9Path:paths[1], transportPath:paths[2]};
  const valid = {components:[component], waterSpecies:null, phase:'frozen'};
  const packet = gasSources => Buffer.from(JSON.stringify({...base, gasSources}));
  assert.deepEqual(validateInput(packet(valid), true).gasSources, valid);
  assert.equal(validateInput(packet({...valid, phase:'liquid-water', waterSpecies:'authored-water'}), true).gasSources.waterSpecies, 'authored-water');
  for (const bad of [null, [], {}, {...valid, phase:undefined}, {...valid, waterSpecies:undefined},
    {...valid, waterSpecies:''}, {...valid, waterSpecies:'none'}, {...valid, waterSpecies:'two words'},
    {...valid, phase:'liquid-water'}, {...valid, phase:'automatic'}, {...valid, components:[]},
    {...valid, components:Array.from({length:129}, () => component)}, {...valid, arbitraryFlag:'--gain'},
    {...valid, components:[{...component, speciesPath:'../outside'}]},
    {...valid, components:[{...component, nasa9Path:'unlisted.nasa9'}]},
    {...valid, components:[{...component, transportPath:undefined}]},
    {...valid, components:[{...component, modelCoefficients:[1, 2, 3]}]}]) {
    assert.throws(() => validateInput(packet(bad), true), /[Gg]as|[Ee]very gas|[Cc]hoose an explicit|liquid-water/);
  }
});

function finiteGasPacket(mode = '') {
  const paths = ['materials/wire.pack', 'gas/local.gas', 'gas/reservoir.gas', 'gas/species', 'gas/nasa9', 'gas/transport'];
  const bytes = Buffer.from('frankensim-finite-gas-string-v1\r\n# read-files ' + JSON.stringify(paths)
    + '\r\n# original temperature and feedback source µ\r\n' + mode + '\r\n');
  const files = paths.map((path, i) => [path, Buffer.from([0, 255, i, 13, 10])]);
  const packet = {...JSON.parse(bundle(bytes, files)), renderer:'finite-gas-string', decimate:false,
    finiteGasString:{solidPath:paths[0], gasRecordPath:paths[1]},
    gasSources:{components:[{speciesPath:paths[3], nasa9Path:paths[4], transportPath:paths[5]}],
      waterSpecies:null, phase:'frozen'}};
  return {packet, paths, files, bytes};
}

test('finite-gas strings run the existing executable with exact physical inputs and return calibrated WAV plus CSV', native, async () => {
  const app = await setup({withFiniteGasString:true, engine:null});
  try {
    assert.equal(app.initial.available, true);
    assert.equal(app.initial.renderers['finite-gas-string'], true);
    assert.equal(app.initial.renderers.wind, false);
    assert.equal(app.initial.maximumHistoryBytes, 16 * 1024 * 1024);
    assert.ok(app.initial.supportedSchemas.includes('frankensim-finite-gas-string-v1'));
    const f = finiteGasPacket();
    const response = await app.runBundle(JSON.stringify(f.packet));
    assert.equal(response.status, 200); assert.match(response.headers.get('content-type'), /application\/json/);
    const result = await response.json();
    assert.equal(result.format, 'frankensim-finite-gas-string-result-v1');
    const wave = Buffer.from(result.wave, 'base64');
    assert.equal(wave.length, 48); assert.equal(wave.readUInt32LE(24), 48000);
    assert.equal(Buffer.from(result.history, 'base64').toString(), 'time_s,pressure_pa\r\n0.000020833333333333333,0\r\n');
    assert.equal(result.receipt.renderer, 'finite-gas-string');
    assert.equal(result.receipt.samples, 2); assert.equal(result.receipt.fullScalePa, 0.01);
    const record = JSON.parse(await readFile(app.receipt, 'utf8'));
    assert.equal(record.bytes, f.bytes.toString('base64'));
    assert.deepEqual(record.args, [join(record.cwd, 'input.performance'), f.paths[0], join(record.cwd, 'output.wav'),
      '--gas-record', f.paths[1], '--self-contained', '--history-output', join(record.cwd, 'output.history.csv'),
      '--require-sample-rate', '48000', ...f.paths.flatMap(path => ['--attached-source', path]),
      '--gas-water', 'none', '--gas-phase', 'frozen', '--gas-component', ...f.paths.slice(3)]);
    assert.deepEqual(record.files.map(file => [file.path, file.bytes]), f.files.map(([path, bytes]) => [path, bytes.toString('base64')]));
    await until(async () => !(await app.status()).busy);
    assert.deepEqual(await readdir(app.tempRoot), []);
    assert.equal((await app.run()).status, 503, 'finite-gas rendering cannot substitute for a missing wind executable');
  } finally { await app.cleanup(); }
});

test('an unconfigured finite-gas executable refuses without substituting another performer', native, async () => {
  const app = await setup({withStrings:true});
  try {
    assert.equal(app.initial.renderers['finite-gas-string'], false);
    const response = await app.runBundle(JSON.stringify(finiteGasPacket().packet));
    assert.equal(response.status, 503); assert.match((await response.json()).error, /finite.gas/i);
    await assert.rejects(readFile(app.receipt), {code:'ENOENT'});
  } finally { await app.cleanup(); }
});

test('finite-gas history, sample-rate and calibration failures never return a partial experiment', native, async () => {
  const app = await setup({withFiniteGasString:true, maximumHistoryBytes:64});
  try {
    assert.equal(app.initial.maximumHistoryBytes, 64);
    for (const mode of ['test-history-large', 'test-history-invalid', 'test-history-missing', 'test-no-calibration', 'test-wrong-rate']) {
      const response = await app.runBundle(JSON.stringify(finiteGasPacket(mode).packet));
      assert.equal(response.status, 502, mode);
      const result = await response.json();
      assert.equal(typeof result.error, 'string'); assert.equal(result.wave, undefined); assert.equal(result.history, undefined);
      await until(async () => !(await app.status()).busy);
      assert.deepEqual(await readdir(app.tempRoot), []);
    }
  } finally { await app.cleanup(); }
});

test('attached material manifests compile through xtask arguments and feed the actual native source paths', native, async () => {
  const app = await setup({withMaterialCompiler:true});
  try {
    assert.equal(app.initial.supportsMaterialManifests, true);
    const names = ['sources/manifest.tsv', 'sources/properties.tsv', 'reed.pack'];
    const input = Buffer.from('frankensim-material-modal-valve-v1\r\n# read-files ' + JSON.stringify(names) + '\r\n');
    const manifest = Buffer.from('original manifest µ\r\n');
    const packet = JSON.parse(bundle(input, [[names[0], manifest], [names[1], Buffer.from('original table\r\n')]]));
    packet.materialManifests = [{packPath:'reed.pack', manifestPath:names[0],
      materialState:{chemistry:'authored composite', phase:'solid', process:'measured-state', revision:0}}];
    for (const table of [Buffer.from('original table\r\n'), Buffer.from('authored data-only material variant\r\n')]) {
      packet.attachments[1].bytes = table.toString('base64');
      const response = await app.runBundle(JSON.stringify(packet));
      assert.equal(response.status, 200); assert.equal((await response.arrayBuffer()).byteLength, 48);
      const compiled = JSON.parse(await readFile(app.compilerReceipt, 'utf8'));
      const record = JSON.parse(await readFile(app.receipt, 'utf8'));
      assert.deepEqual(compiled.args, ['matdb-pack', '--manifest', join(record.cwd, names[0]),
        '--out', join(record.cwd, 'reed.pack'), '--material-state',
        'authored composite', 'solid', 'measured-state', '0']);
      assert.equal(compiled.manifest, manifest.toString('base64'));
      assert.equal(compiled.table, table.toString('base64'));
      const expectedPack = Buffer.concat([Buffer.from('compiled transport witness\0'), manifest, table]);
      assert.deepEqual(record.files.map(file => [file.path, file.bytes]),
        [[names[0], manifest.toString('base64')], [names[1], table.toString('base64')],
          ['reed.pack', expectedPack.toString('base64')]]);
      assert.equal(record.bytes, input.toString('base64'));
      assert.deepEqual(record.args, ['wind', join(record.cwd, 'input.performance'), join(record.cwd, 'output.wav'),
        '--self-contained', ...names.flatMap(name => ['--attached-source', name])]);
      await until(async () => !(await app.status()).busy);
    }
  } finally { await app.cleanup(); }
});

test('unavailable, refusing or oversized material compilation never starts the native performer', native, async () => {
  for (const [mode, expected] of [[null, 503], ['compile-refuse', 400], ['compile-large', 502]]) {
    const app = await setup({withMaterialCompiler:mode !== null});
    try {
      const packet = JSON.parse(bundle(source(), [
        ['source/manifest.tsv', Buffer.from(mode ?? 'original manifest')],
        ['source/properties.tsv', Buffer.from('original source table')]]));
      packet.materialManifests = [{packPath:'reed.pack', manifestPath:'source/manifest.tsv'}];
      assert.equal(app.initial.supportsMaterialManifests, mode !== null);
      const response = await app.runBundle(JSON.stringify(packet));
      assert.equal(response.status, expected);
      assert.match((await response.json()).error, /compiler|material pack/);
      await until(async () => !(await app.status()).busy);
      await assert.rejects(readFile(app.receipt), {code:'ENOENT'});
    } finally { await app.cleanup(); }
  }
});

test('browser cancellation includes a running material compiler and prevents subsequent rendering', native, async () => {
  const app = await setup({withMaterialCompiler:true});
  try {
    const packet = JSON.parse(bundle(source(), [
      ['source/manifest.tsv', Buffer.from('compile-hang')],
      ['source/properties.tsv', Buffer.from('original source table')]]));
    packet.materialManifests = [{packPath:'reed.pack', manifestPath:'source/manifest.tsv'}];
    const controller = new AbortController();
    const pending = app.run(JSON.stringify(packet), {signal:controller.signal,
      headers:{'Content-Type':'application/json', 'X-Native-Wind-Token':app.initial.token}}).catch(error => error);
    await until(async () => { try { await readFile(app.compilerReceipt); return true; } catch { return false; } });
    const record = JSON.parse(await readFile(app.compilerReceipt, 'utf8'));
    controller.abort(); await pending;
    await until(async () => !(await app.status()).busy);
    assert.throws(() => process.kill(record.pid, 0));
    await assert.rejects(readFile(app.receipt), {code:'ENOENT'});
  } finally { await app.cleanup(); }
});
