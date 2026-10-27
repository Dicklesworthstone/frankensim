// Actual HTTP and subprocess tests with an explicitly synthetic executable.
// This proves the adapter, not Rust execution, acoustics, or instrument realism.
import {test} from 'node:test';
import assert from 'node:assert/strict';
import {mkdtemp, mkdir, readFile, readdir, rm, symlink, writeFile} from 'node:fs/promises';
import {tmpdir} from 'node:os';
import {join, isAbsolute} from 'node:path';
import http from 'node:http';
import {createNativeWindServer} from '../tools/serve-native-wind.mjs';
import {MAX_SOURCE_BYTES, validateSource} from '../tools/native-wind-protocol.mjs';

const source = (mode = 'ok') => Buffer.from(`frankensim-reed-performance-v1\r\n${mode}\r\n# UTF-8: 音; $() '--help'\r\n`);
const pause = ms => new Promise(resolve => setTimeout(resolve, ms));
async function until(fn) { for (let i = 0; i < 400; i++) { if (await fn()) return; await pause(10); } throw new Error('Timed out waiting for test condition.'); }
async function setup(options = {}) {
  const dir = await mkdtemp(join(tmpdir(), 'native-wind-test-')), dist = join(dir, 'dist'), tempRoot = join(dir, 'jobs');
  await mkdir(dist); await mkdir(tempRoot); await writeFile(join(dist, 'index.html'), '<title>Synthetic adapter test</title>');
  const receipt = join(dir, 'receipt.json'), engine = join(dir, 'fixture.mjs');
  const fixture = `#!${process.execPath}
// Synthetic process fixture: never presented as native physics.
import {readFileSync, writeFileSync, statSync} from 'node:fs';
const args = process.argv.slice(2), input = readFileSync(args[1]);
writeFileSync(${JSON.stringify(receipt)}, JSON.stringify({args, cwd:process.cwd(), bytes:input.toString('base64'), pid:process.pid,
 directoryMode:statSync(process.cwd()).mode & 0o777, sourceMode:statSync(args[1]).mode & 0o777}));
const text = input.toString();
if (text.includes('test-refuse')) { console.error('synthetic parser refusal'); process.exit(2); }
if (text.includes('test-log')) { process.stderr.write(Buffer.alloc(300000, 120)); await new Promise(() => setInterval(() => {}, 1000)); }
if (text.includes('test-hang')) { process.on('SIGTERM', () => {}); await new Promise(() => setInterval(() => {}, 1000)); }
if (text.includes('test-invalid')) { writeFileSync(args[2], 'not a wave'); process.exit(0); }
if (text.includes('test-large')) { writeFileSync(args[2], Buffer.alloc(4097)); process.exit(0); }
if (text.includes('test-delay')) await new Promise(resolve => setTimeout(resolve, 200));
const wav = Buffer.alloc(48); wav.write('RIFF'); wav.writeUInt32LE(40, 4); wav.write('WAVEfmt ', 8);
wav.writeUInt32LE(16, 16); wav.writeUInt16LE(1, 20); wav.writeUInt16LE(1, 22); wav.writeUInt32LE(48000, 24);
wav.writeUInt32LE(96000, 28); wav.writeUInt16LE(2, 32); wav.writeUInt16LE(16, 34); wav.write('data', 36); wav.writeUInt32LE(4, 40);
writeFileSync(args[2], wav);
`;
  await writeFile(engine, fixture, {mode:0o700});
  const app = createNativeWindServer({engine, dist, tempRoot, ...options});
  await new Promise((accept, reject) => {app.server.once('error', reject); app.server.listen(0, '127.0.0.1', accept);});
  const url = `http://127.0.0.1:${app.server.address().port}`;
  const status = () => fetch(url + '/api/native-wind/status').then(r => r.json());
  const initial = await status();
  const run = (bytes = source(), extra = {}) => fetch(url + '/api/native-wind/render', {method:'POST',
    headers:{'Content-Type':'text/plain;charset=utf-8', 'X-Native-Wind-Token':initial.token}, body:bytes, ...extra});
  return {...app, url, dir, receipt, tempRoot, status, initial, run,
    async cleanup() { await app.shutdown(); await rm(dir, {recursive:true, force:true}); }};
}
const native = {skip:process.platform === 'win32'};

test('transport retains exact UTF-8 and admits only self-contained native source families', () => {
  const bytes = source(); assert.equal(validateSource(bytes), bytes);
  const flute = Buffer.from('frankensim-flue-dg-performance-v1\nopaque'); assert.equal(validateSource(flute), flute);
  for (const bad of [Buffer.alloc(0), Buffer.from([0xff]), Buffer.from('frankensim-material-modal-valve-v1\nsection /etc/passwd'), Buffer.alloc(MAX_SOURCE_BYTES + 1)])
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
    assert.equal((await app.run(source(), {headers:{'Content-Type':'application/json', 'X-Native-Wind-Token':app.initial.token}})).status, 415);
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
    assert.equal(receipt.bytes, bytes.toString('base64')); assert.equal(receipt.args.length, 3); assert.equal(receipt.args[0], 'wind');
    assert.equal(receipt.args[1], join(receipt.cwd, 'input.performance')); assert.equal(receipt.args[2], join(receipt.cwd, 'output.wav'));
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
