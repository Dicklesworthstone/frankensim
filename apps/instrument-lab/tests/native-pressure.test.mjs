// Manufactured PCM and explicit subprocess fixtures. No native physics claim.
import {test} from 'node:test';
import assert from 'node:assert/strict';
import {createHash} from 'node:crypto';
import {mkdtemp, writeFile, mkdir, open} from 'node:fs/promises';
import {join} from 'node:path';
import {tmpdir} from 'node:os';
import {inspectNativeWave, validatePressureReceipt, verifyPressureReceipt, nativePcmSamples,
  compareNativePressure, PRESSURE_RECEIPT_FORMAT, PRESSURE_RECEIPT_HEADER} from '../src/native-wave.mjs';
import {nativePressureReceipt, MAX_PRESSURE_SIDECAR_BYTES} from '../tools/native-pressure-receipt.mjs';
import {createNativeWindServer} from '../tools/serve-native-wind.mjs';
import {validateWave} from '../tools/native-wind-protocol.mjs';
import {nativeWindStatus, renderNativePerformance, renderNativeWindBundle} from '../src/native-wind-client.mjs';

function wave(samples = [8000, -8000, 4000, -4000], junk = false) {
  const extra = junk ? 10 : 0, out = Buffer.alloc(44 + extra + samples.length * 2);
  out.write('RIFF'); out.writeUInt32LE(out.length - 8, 4); out.write('WAVE', 8);
  let pos = 12;
  if (junk) { out.write('JUNK', pos); out.writeUInt32LE(1, pos + 4); out[pos + 8] = 255; pos += 10; }
  out.write('fmt ', pos); out.writeUInt32LE(16, pos + 4); out.writeUInt16LE(1, pos + 8);
  out.writeUInt16LE(1, pos + 10); out.writeUInt32LE(48000, pos + 12);
  out.writeUInt32LE(96000, pos + 16); out.writeUInt16LE(2, pos + 20); out.writeUInt16LE(16, pos + 22);
  pos += 24; out.write('data', pos); out.writeUInt32LE(samples.length * 2, pos + 4);
  samples.forEach((sample, i) => out.writeInt16LE(sample, pos + 8 + 2 * i));
  return out;
}
function receipt(bytes, changes = {}) {
  return {format:PRESSURE_RECEIPT_FORMAT, renderer:'wind', sampleRateHz:48000,
    samples:inspectNativeWave(bytes).samples, fullScalePa:200, clippedSamples:0,
    wavSha256:createHash('sha256').update(bytes).digest('hex'), ...changes};
}
const close = (a, b) => assert.ok(Math.abs(a - b) <= 1e-12 * Math.max(1, Math.abs(b)), `${a} != ${b}`);
const scratch = () => mkdtemp(join(tmpdir(), 'frankensim-pressure-test-'));

test('one WAV reader handles padded chunks and offset byte views without altering native PCM', () => {
  for (const junk of [false, true]) {
    const bytes = wave([32767, -32768, 0, 12], junk), padded = Buffer.concat([Buffer.from([1, 2, 3]), bytes]);
    const view = padded.subarray(3), info = inspectNativeWave(view);
    assert.equal(info.samples, 4); assert.equal(info.dataOffset, junk ? 54 : 44);
    assert.equal(validateWave(view), view);
    assert.deepEqual([...nativePcmSamples(view)], [32767 / 32768, -1, 0, 12 / 32768]);
    assert.deepEqual(view, bytes);
  }
});

test('service and browser share refusal of incomplete or non-native PCM', () => {
  const good = wave(), variants = [Buffer.alloc(0), good.subarray(0, 43), Buffer.concat([good, Buffer.from([0])]), wave([])];
  for (const [offset, value, bits] of [[20, 3, 16], [22, 2, 16], [24, 44100, 32],
    [28, 12345, 32], [32, 4, 16], [34, 24, 16], [40, 9999999, 32]]) {
    const bad = Buffer.from(good); bad[bits === 16 ? 'writeUInt16LE' : 'writeUInt32LE'](value, offset); variants.push(bad);
  }
  const duplicate = Buffer.concat([good, good.subarray(36)]); duplicate.writeUInt32LE(duplicate.length - 8, 4); variants.push(duplicate);
  for (const bytes of variants) {
    assert.throws(() => inspectNativeWave(bytes), /WAV/);
    assert.throws(() => validateWave(bytes), /WAV/);
  }
  assert.throws(() => inspectNativeWave(good, good.length - 1), /WAV/);
});

test('pressure receipt is bound to native output dimensions, performer and exact WAV bytes', async () => {
  const bytes = wave(), source = receipt(bytes), original = Buffer.from(bytes);
  assert.deepEqual(await verifyPressureReceipt(source, bytes, 'wind'), source);
  assert.ok(Object.isFrozen(validatePressureReceipt(source, inspectNativeWave(bytes))));
  for (const changes of [{format:'other'}, {renderer:'piano'}, {samples:5}, {sampleRateHz:96000},
    {fullScalePa:0}, {fullScalePa:Infinity}, {fullScalePa:NaN}, {clippedSamples:-1},
    {clippedSamples:5}, {clippedSamples:0.5}, {wavSha256:'missing'}]) {
    assert.throws(() => validatePressureReceipt({...source, ...changes}, inspectNativeWave(bytes)), /calibration/);
  }
  await assert.rejects(verifyPressureReceipt(source, bytes, 'strings'), /calibration/);
  const changed = Buffer.from(bytes); changed[44] ^= 1;
  await assert.rejects(verifyPressureReceipt(source, changed), /different WAV/);
  assert.deepEqual(bytes, original);
});

test('comparison retains a fourfold amplitude change at one native pressure scale', () => {
  const a = wave(), b = wave([2000, -2000, 1000, -1000]);
  const result = compareNativePressure(a, receipt(a), b, receipt(b));
  close(result.rmsB / result.rmsA, 0.25); close(result.peakB / result.peakA, 0.25);
  close(result.rmsDifference / result.rmsA, 0.75);
  close(result.levelDifferenceDb, 20 * Math.log10(0.25));
  close(result.peakA, 8000 / 32767 * 200);
  assert.equal(result.identical, false); assert.equal(result.truncated, false);
  assert.equal(compareNativePressure(a, receipt(a), a, receipt(a)).identical, true);
});

test('comparison keeps original time zero and explicitly bounds unequal durations', () => {
  const a = wave([100, 0, 0, 0]), b = wave([0, 100, 0]);
  const result = compareNativePressure(a, receipt(a), b, receipt(b));
  assert.equal(result.samples, 3); assert.equal(result.truncated, true);
  assert.equal(result.durationSeconds, 3 / 48000); assert.ok(result.rmsDifference > 0);
  close(result.levelDifferenceDb, 0);
  const silent = wave([0, 0, 0]);
  assert.equal(compareNativePressure(silent, receipt(silent), b, receipt(b)).levelDifferenceDb, null);
});

test('no hidden calibration conversion and no calibrated comparison of clipped signals', () => {
  const a = wave(), b = wave();
  assert.throws(() => compareNativePressure(a, receipt(a), b, receipt(b, {fullScalePa:100})), /same native full-scale/);
  for (const side of ['a', 'b']) {
    assert.throws(() => compareNativePressure(a, receipt(a, {clippedSamples:side === 'a' ? 1 : 0}),
      b, receipt(b, {clippedSamples:side === 'b' ? 1 : 0})), /clipped/);
  }
});

test('wind calibration comes from the native sidecar and strings from the output-clock summary', async () => {
  const dir = await scratch(), output = join(dir, 'result.wav'), bytes = wave();
  assert.equal(await nativePressureReceipt('wind', output, {stderr:''}, bytes), null);
  const native = {schema:'frankensim-music-render-provenance-v1', sample_rate_hz:48000,
    samples:4, full_scale_pa:200, clipped_samples:0};
  await writeFile(join(dir, 'result.provenance.json'), JSON.stringify(native));
  assert.deepEqual(await nativePressureReceipt('wind', output, {stderr:'irrelevant'}, bytes), receipt(bytes));
  assert.deepEqual(await nativePressureReceipt('resolved', output, {stderr:'irrelevant'}, bytes),
    receipt(bytes, {renderer:'resolved'}));
  assert.deepEqual(await nativePressureReceipt('strings', output,
    {stderr:'unrelated\nsamples=4; rate_Hz=48000; clips=0; full_scale_Pa=200\nobstacle_resolves=0\n'}, bytes),
    receipt(bytes, {renderer:'strings'}));
  assert.equal(await nativePressureReceipt('strings', output, {stderr:'legacy output'}, bytes), null);
  await writeFile(join(dir, 'result.provenance.json'), JSON.stringify({...native, samples:8}));
  await assert.rejects(nativePressureReceipt('wind', output, {stderr:''}, bytes), /calibration/);
});

test('malformed emitted calibration refuses instead of silently becoming uncalibrated', async () => {
  const dir = await scratch(), output = join(dir, 'result.wav'), path = join(dir, 'result.provenance.json'), bytes = wave();
  for (const text of ['not JSON', '{}', '{"schema":"other"}', '\xff']) {
    await writeFile(path, text); await assert.rejects(nativePressureReceipt('wind', output, {stderr:''}, bytes), /calibration/);
  }
  const file = await open(path, 'w'); await file.truncate(MAX_PRESSURE_SIDECAR_BYTES + 1); await file.close();
  await assert.rejects(nativePressureReceipt('wind', output, {stderr:''}, bytes), /calibration/);
  for (const stderr of ['samples=4; rate_Hz=96000; clips=0; full_scale_Pa=200',
    'samples=4; rate_Hz=48000; clips=0; full_scale_Pa=NaN', 'samples=4; incomplete',
    'samples=4; rate_Hz=48000; clips=0; full_scale_Pa=200\nsamples=4; rate_Hz=48000; clips=0; full_scale_Pa=200']) {
    await assert.rejects(nativePressureReceipt('strings', output, {stderr}, bytes), /calibration/);
  }
});

async function server(t) {
  const dir = await scratch(), dist = join(dir, 'dist'), engine = join(dir, 'fixture.mjs');
  await mkdir(dist);
  const bytes = wave();
  await writeFile(engine, `#!${process.execPath}
// Explicit synthetic native-output fixture, not a physical simulation.
import {readFileSync, writeFileSync} from 'node:fs';
import {dirname, join} from 'node:path';
const args = process.argv.slice(2), strings = !['wind', 'resolved'].includes(args[0]), offset = strings ? 0 : 1;
const text = readFileSync(args[offset], 'utf8'), output = args[offset + 1];
writeFileSync(output, Buffer.from('${bytes.toString('base64')}', 'base64'));
if (!text.includes('no-receipt')) {
  const count = text.includes('wrong-receipt') ? 99 : 4;
  const clips = text.includes('clipped-receipt') ? 1 : 0;
  if (strings) {
    process.stderr.write('samples=' + count + '; rate_Hz=48000; ');
    await new Promise(resolve => setTimeout(resolve, 5));
    process.stdout.write('separate stdout must not corrupt the stderr summary\\n');
    process.stderr.write('clips=' + clips + '; full_scale_Pa=200\\n');
  } else writeFileSync(join(dirname(output), 'output.provenance.json'), JSON.stringify({
    schema:'frankensim-music-render-provenance-v1', sample_rate_hz:48000,
    samples:count, clipped_samples:clips, full_scale_pa:200}));
}
`, {mode:0o700});
  const app = createNativeWindServer({engine, stringsEngine:engine, dist, tempRoot:dir});
  t.after(() => app.shutdown());
  await new Promise((accept, reject) => { app.server.once('error', reject); app.server.listen(0, '127.0.0.1', accept); });
  const url = `http://127.0.0.1:${app.server.address().port}`;
  const request = (path, init) => fetch(url + path, init);
  const status = await nativeWindStatus(undefined, request);
  const ready = async () => {
    for (let i = 0; i < 200; i++) {
      if (!(await nativeWindStatus(undefined, request)).busy) return;
      await new Promise(resolve => setTimeout(resolve, 5));
    }
    throw new Error('fixture cleanup did not finish');
  };
  return {request, status, bytes, ready};
}

for (const renderer of ['wind', 'strings', 'resolved']) {
  test(`${renderer}: native service and client preserve WAV and expose matching pressure calibration`, {skip:process.platform === 'win32'}, async t => {
    const app = await server(t), header = {wind:'frankensim-reed-performance-v1',
      strings:'frankensim-material-strings-v3', resolved:'frankensim-resolved-dg-performance-v1'}[renderer];
    assert.equal(app.status.supportsPressureReceipts, true);
    const result = await renderNativePerformance(header + '\n', [], app.status, {renderer}, undefined, app.request);
    assert.deepEqual(Buffer.from(await result.wave.arrayBuffer()), app.bytes);
    assert.deepEqual(result.receipt, receipt(app.bytes, {renderer}));
    await app.ready();
    const legacy = await renderNativePerformance(header + '\nno-receipt\n', [], app.status, {renderer}, undefined, app.request);
    assert.equal(legacy.receipt, null); assert.deepEqual(Buffer.from(await legacy.wave.arrayBuffer()), app.bytes);
    await app.ready();
    const blob = await renderNativeWindBundle(header + '\n', [], app.status, {renderer}, undefined, app.request);
    assert.ok(blob instanceof Blob); assert.deepEqual(Buffer.from(await blob.arrayBuffer()), app.bytes);
    await app.ready();
    await assert.rejects(renderNativePerformance(header + '\nwrong-receipt\n', [], app.status,
      {renderer}, undefined, app.request), /calibration/);
    await app.ready();
    const clipped = await renderNativePerformance(header + '\nclipped-receipt\n', [], app.status,
      {renderer}, undefined, app.request);
    assert.equal(clipped.receipt.clippedSamples, 1);
  });
}

test('client refuses a corrupt or mismatched calibration header even on otherwise valid WAV', async () => {
  const bytes = wave(), status = {available:true, token:'test', maximumSourceBytes:1000, maximumOutputBytes:1000};
  for (const header of ['not JSON', 'x'.repeat(2049), JSON.stringify(receipt(bytes, {wavSha256:'0'.repeat(64)})),
    JSON.stringify(receipt(bytes, {renderer:'strings'}))]) {
    const request = async () => new Response(bytes, {headers:{'Content-Type':'audio/wav', [PRESSURE_RECEIPT_HEADER]:header}});
    await assert.rejects(renderNativePerformance('source', [], status, {}, undefined, request), /calibration/);
  }
});
