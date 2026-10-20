// Execute the actual worker selection/rebuild path with recording WASM and
// browser boundaries. These tests check delivery of physical inputs, not sound.
import test from 'node:test';
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { runInNewContext } from 'node:vm';

const url = new URL('../public/live/kit/worker.js', import.meta.url);
const source = (await readFile(url, 'utf8'))
  .replace(/^import .*;\n/gm, '')
  .replaceAll('import.meta.url', JSON.stringify(url.href));

async function run(gasRecord, { copies = true, slowFirst = false, refuseFine = false,
  missingWithGas = false, refuseGas = false, scales = [3, 5] } = {}) {
  const calls = [], messages = [], kept = [];
  let id = 0, active = 0, clock = 0, attached = 0;
  function piece(name, scale, record) {
    const instance = {
      id: ++id, name, scale, record,
      copy() {
        if (!copies) throw new Error('old engine cannot copy');
        const copy = { ...this }; kept.push(copy); return copy;
      },
      strike() {}, radius() { return 0.2; }, free() {}, asleep() { return true; },
      render(n) { active = this.id; return new Float32Array(n); },
      info() { return JSON.stringify({ piece: name, gasRecord: record }); },
      surfaceXY() { return new Float32Array([0, 0]); }, surfaceH() { return 0.01; },
    };
    return instance;
  }
  class KitPiece {
    constructor(name, scale) { calls.push(['default', name, scale]); return piece(name, scale); }
    static withGas(name, scale, record) {
      calls.push(['gas', name, scale, record]);
      if (refuseGas) throw new Error('outside admitted gas domain');
      if (refuseFine && scale < 5) throw new Error('fine grid needs another mechanical clock');
      return piece(name, scale, record);
    }
  }
  if (missingWithGas) KitPiece.withGas = undefined;
  const self = { postMessage(message) { messages.push(message); } };
  runInNewContext(source, { self, KitPiece, init: async () => {},
    createLane: () => { attached++; return { frame: 0 }; }, levelFix: () => 1,
    Float32Array, Uint8Array, URL,
    performance: { now() { clock += slowFirst && active === 1 ? 10 : 1; return clock; } },
    fetch: async () => ({ ok: true, clone() { return this; }, arrayBuffer: async () => new ArrayBuffer(0),
      json: async () => ({}) }),
  });
  await self.onmessage({ data: { type: 'init', piece: 'snare', gasRecord, port: {}, scales, minRealtime: 1.7 } });
  return { calls, messages, kept, attached };
}

test('the chosen clone and every automatic grid trial retain the exact mixture', async () => {
  const record = 'mixture 300 80000 mole 3 N2 0.7 O2 0.2 CO2 0.1';
  const result = await run(record);
  assert.deepEqual(result.calls, [['gas', 'snare', 5, record], ['gas', 'snare', 3, record]]);
  assert.ok(result.kept.every(copy => copy.record === record));
  const ready = result.messages.find(m => m.type === 'ready');
  assert.equal(ready.gasRecord, record);
  assert.equal(ready.info.gasRecord, record);
  assert.equal(ready.scale, 3);
  assert.equal(result.attached, 1);
});

test('warmup retries and old-WASM clone fallbacks use the same gas constructor', async () => {
  const record = '293.15 65000 0.5';
  const result = await run(record, { copies: false, slowFirst: true });
  assert.deepEqual(result.calls, [['gas', 'snare', 5, record], ['gas', 'snare', 5, record],
    ['gas', 'snare', 3, record], ['gas', 'snare', 3, record]]);
  assert.equal(result.messages.at(-1).gasRecord, record);
});

test('fallback after every requested grid refuses still retains mass-fraction composition', async () => {
  const record = 'mixture 295 100000 mass 2 N2 0.8 CO2 0.2';
  const result = await run(record, { scales: [2, 3], copies: false, refuseFine: true });
  assert.deepEqual(result.calls, [['gas', 'snare', 3, record], ['gas', 'snare', 2, record],
    ['gas', 'snare', 5, record], ['gas', 'snare', 5, record]]);
  assert.equal(result.messages.at(-1).gasRecord, record);
});

test('missing atmosphere support or rejected material data never silently creates default air', async () => {
  for (const options of [{ missingWithGas: true }, { refuseGas: true }]) {
    const result = await run('mixture 300 80000 mole 1 CO2 1', options);
    assert.equal(result.attached, 0);
    assert.equal(result.messages.at(-1).type, 'error');
    assert.match(result.messages.at(-1).message, /rebuilding with atmosphere support|outside admitted gas domain/);
    assert.ok(result.calls.every(call => call[0] === 'gas'));
  }
  for (const record of ['', null, 7]) {
    const result = await run(record);
    assert.equal(result.attached, 0);
    assert.match(result.messages.at(-1).message, /nonempty atmosphere record/);
    assert.equal(result.calls.length, 0);
  }
});

test('a restarted worker adopts its new atmosphere while omitted records keep the standard path', async () => {
  const air = await run('293.15 101325 0');
  const mixture = await run('mixture 303.15 95000 mole 1 Ar 1');
  assert.equal(air.messages.at(-1).info.gasRecord, '293.15 101325 0');
  assert.equal(mixture.messages.at(-1).info.gasRecord, 'mixture 303.15 95000 mole 1 Ar 1');
  const legacy = await run(undefined);
  assert.ok(legacy.calls.every(call => call[0] === 'default'));
  assert.equal(legacy.messages.at(-1).type, 'ready');
});
