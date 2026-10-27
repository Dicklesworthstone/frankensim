// Actual HTTP and subprocess tests with an explicitly synthetic executable.
// This proves the adapter, not Rust execution, acoustics, or instrument realism.
import {test} from 'node:test';
import assert from 'node:assert/strict';
import {mkdtemp, mkdir, readFile, readdir, rm, symlink, writeFile} from 'node:fs/promises';
import {tmpdir} from 'node:os';
import {join, isAbsolute} from 'node:path';
import http from 'node:http';
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
  const {withStrings = false, ...serverOptions} = options;
  const dir = await mkdtemp(join(tmpdir(), 'native-wind-test-')), dist = join(dir, 'dist'), tempRoot = join(dir, 'jobs');
  await mkdir(dist); await mkdir(tempRoot); await writeFile(join(dist, 'index.html'), '<title>Synthetic adapter test</title>');
  const receipt = join(dir, 'receipt.json'), engine = join(dir, 'fixture.mjs');
  const fixture = `#!${process.execPath}
// Synthetic process fixture: never presented as native physics.
import {readFileSync, writeFileSync, statSync} from 'node:fs';
const args = process.argv.slice(2), offset = args[0] === 'wind' ? 1 : 0;
const inputPath = args[offset], outputPath = args[offset + 1], input = readFileSync(inputPath);
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
writeFileSync(outputPath, wav);
`;
  await writeFile(engine, fixture, {mode:0o700});
  const app = createNativeWindServer({engine, stringsEngine:withStrings ? engine : null, dist, tempRoot, ...serverOptions});
  await new Promise((accept, reject) => {app.server.once('error', reject); app.server.listen(0, '127.0.0.1', accept);});
  const url = `http://127.0.0.1:${app.server.address().port}`;
  const status = () => fetch(url + '/api/native-wind/status').then(r => r.json());
  const initial = await status();
  const run = (bytes = source(), extra = {}) => fetch(url + '/api/native-wind/render', {method:'POST',
    headers:{'Content-Type':'text/plain;charset=utf-8', 'X-Native-Wind-Token':initial.token}, body:bytes, ...extra});
  const runBundle = bytes => run(bytes, {headers:{'Content-Type':'application/json', 'X-Native-Wind-Token':initial.token}});
  return {...app, url, dir, receipt, tempRoot, status, initial, run, runBundle,
    async cleanup() { await app.shutdown(); await rm(dir, {recursive:true, force:true}); }};
}
const native = {skip:process.platform === 'win32'};

test('material strings reuse the exact source packet, gas associations and explicit profile paths', native, async () => {
  const app = await setup({withStrings:true, engine:null});
  try {
    assert.deepEqual(app.initial.renderers, {wind:false, strings:true});
    assert.equal(app.initial.available, true);
    assert.ok(app.initial.supportedSchemas.includes('frankensim-material-strings-v3'));
    const names = ['wire.pack', 'body.pack', 'shape.body', 'finger.stops', 'gas/s', 'gas/n', 'gas/t'];
    const bytes = Buffer.from('frankensim-material-strings-v3\r\n# read-files ' + JSON.stringify(names) + '\r\n# original 音\r\n');
    const files = names.map((name, i) => [name, Buffer.from([0, 255, i, 13, 10])]);
    const packet = JSON.parse(bundle(bytes, files));
    packet.renderer = 'strings'; packet.decimate = true;
    packet.stringProfiles = {bodyPath:'shape.body', stopPath:'finger.stops'};
    packet.gasSources = {components:[{speciesPath:'gas/s', nasa9Path:'gas/n', transportPath:'gas/t'}], waterSpecies:null, phase:'frozen'};
    const response = await app.runBundle(JSON.stringify(packet));
    assert.equal(response.status, 200); assert.equal((await response.arrayBuffer()).byteLength, 48);
    const receipt = JSON.parse(await readFile(app.receipt, 'utf8'));
    assert.equal(receipt.bytes, bytes.toString('base64'));
    assert.deepEqual(receipt.args, [join(receipt.cwd, 'input.performance'), join(receipt.cwd, 'output.wav'),
      '--self-contained', '--decimate', ...names.flatMap(name => ['--attached-source', name]),
      '--gas-water', 'none', '--gas-phase', 'frozen', '--gas-component', 'gas/s', 'gas/n', 'gas/t',
      '--body-profile', 'shape.body', '--stop-profile', 'finger.stops']);
    assert.deepEqual(receipt.files.map(file => [file.path, file.bytes]), files.map(([path, bytes]) => [path, bytes.toString('base64')]));
    await until(async () => !(await app.status()).busy);
    assert.equal((await app.run()).status, 503, 'a configured string performer must not substitute for a missing wind performer');
  } finally { await app.cleanup(); }
});

test('renderer and profile declarations refuse mismatched schemas, unlisted paths and arbitrary flags', () => {
  const bytes = Buffer.from('frankensim-material-strings-v1\nopaque');
  for (const renderer of ['strings', 'wind', 'shell', null]) {
    const packet = {...JSON.parse(bundle(bytes, [['body.profile', Buffer.from('original')]])), renderer};
    if (renderer === 'strings') {
      assert.equal(validateInput(Buffer.from(JSON.stringify(packet)), true).renderer, 'strings');
      for (const stringProfiles of [{bodyPath:'missing'}, {stopPath:'../body.profile'}, {observerStore:'body.profile'}, {bodyPath:null}])
        assert.throws(() => validateInput(Buffer.from(JSON.stringify({...packet, stringProfiles})), true), /String profiles/);
    } else assert.throws(() => validateInput(Buffer.from(JSON.stringify(packet)), true));
  }
  const wind = {...JSON.parse(bundle(source(), [['body.profile', Buffer.from('original')]])), stringProfiles:{bodyPath:'body.profile'}};
  assert.throws(() => validateInput(Buffer.from(JSON.stringify(wind)), true), /String profiles/);
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
