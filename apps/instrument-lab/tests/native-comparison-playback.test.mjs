// Exercises the actual shared player against an observable AudioContext boundary.
import {test} from 'node:test';
import assert from 'node:assert/strict';
import {createHash} from 'node:crypto';
import {NativeComparisonPlayback} from '../src/native-comparison-playback.mjs';
import {PRESSURE_RECEIPT_FORMAT} from '../src/native-wave.mjs';

function wave(value, count = 4800) {
  const bytes = Buffer.alloc(44 + 2 * count);
  bytes.write('RIFF'); bytes.writeUInt32LE(bytes.length - 8, 4); bytes.write('WAVEfmt ', 8);
  bytes.writeUInt32LE(16, 16); bytes.writeUInt16LE(1, 20); bytes.writeUInt16LE(1, 22);
  bytes.writeUInt32LE(48000, 24); bytes.writeUInt32LE(96000, 28);
  bytes.writeUInt16LE(2, 32); bytes.writeUInt16LE(16, 34); bytes.write('data', 36);
  bytes.writeUInt32LE(2 * count, 40);
  for (let i = 0; i < count; i++) bytes.writeInt16LE(value, 44 + 2 * i);
  return bytes;
}
const receipt = bytes => ({format:PRESSURE_RECEIPT_FORMAT, renderer:'strings', samples:(bytes.length - 44) / 2,
  sampleRateHz:48000, fullScalePa:200, clippedSamples:0, wavSha256:createHash('sha256').update(bytes).digest('hex')});
function context() {
  const sourceNodes = [], gainNodes = [];
  const node = () => ({connections:[], disconnects:0,
    connect(other) { this.connections.push(other); }, disconnect() { this.disconnects++; }});
  return {currentTime:10, state:'suspended', destination:{}, sourceNodes, gainNodes,
    async resume() { this.state = 'running'; }, async close() { this.state = 'closed'; },
    createBuffer(channels, count, rate) {
      assert.equal(channels, 1); assert.equal(rate, 48000);
      return {data:new Float32Array(count), copyToChannel(samples, channel) { assert.equal(channel, 0); this.data.set(samples); }};
    },
    createBufferSource() {
      const source = {...node(), starts:[], stops:0, onended:null,
        start(...args) { this.starts.push(args); }, stop() { this.stops++; }};
      sourceNodes.push(source); return source;
    },
    createGain() {
      const gain = {...node(), gain:{value:1, events:[],
        cancelScheduledValues(time) { this.events.push(['cancel', time]); },
        setValueAtTime(value, time) { this.events.push(['set', value, time]); },
        linearRampToValueAtTime(value, time) { this.events.push(['ramp', value, time]); }}};
      gainNodes.push(gain); return gain;
    },
  };
}
function setup(ctx = context()) {
  const player = new NativeComparisonPlayback(() => ctx), a = wave(8000), b = wave(2000);
  const stats = player.setPair(a, receipt(a), b, receipt(b));
  return {player, ctx, a, b, stats};
}
const near = (a, b) => assert.ok(Math.abs(a - b) < 1e-12, `${a} != ${b}`);

test('both native streams use one clock, offset, window and shared volume without normalization', async () => {
  const {player, ctx, a, stats} = setup();
  near(stats.rmsB / stats.rmsA, 0.25);
  a.fill(0, 44); // Changing caller memory cannot mutate the retained reference.
  assert.equal(await player.play(0.025), true);
  const [first, second] = ctx.sourceNodes;
  assert.deepEqual(first.starts, second.starts);
  assert.deepEqual(first.starts[0], [10.02, 0.025, 0.07500000000000001]);
  assert.equal(first.buffer.data[0], 8000 / 32768);
  assert.equal(second.buffer.data[0], 2000 / 32768);
  assert.equal(ctx.gainNodes[0].gain.value, 0.25);
  assert.equal(ctx.gainNodes[1].gain.value, 1); assert.equal(ctx.gainNodes[2].gain.value, 0);
  near(player.positionSeconds, 0.025); // Scheduled start has not happened yet.
  ctx.currentTime = 10.045; near(player.positionSeconds, 0.05);
  player.select('b'); player.setVolume(0.4);
  assert.equal(ctx.sourceNodes.length, 2); assert.equal(first.starts.length, 1);
  assert.deepEqual(ctx.gainNodes[0].gain.events.at(-1), ['set', 0.4, 10.045]);
  assert.deepEqual(ctx.gainNodes[1].gain.events.at(-1), ['ramp', 0, 10.05]);
  assert.deepEqual(ctx.gainNodes[2].gain.events.at(-1), ['ramp', 1, 10.05]);
  first.onended(); assert.equal(player.playing, false);
  assert.ok(ctx.sourceNodes.every(node => node.disconnects === 1));
});

test('rapid switches continue complementary crossfades rather than restarting either note', async () => {
  const {player, ctx} = setup(); await player.play();
  ctx.currentTime = 10.03; player.select('b');
  ctx.currentTime = 10.032; player.select('a');
  const a = ctx.gainNodes[1].gain.events.at(-2), b = ctx.gainNodes[2].gain.events.at(-2);
  assert.equal(a[0], 'set'); assert.equal(b[0], 'set'); near(a[1] + b[1], 1);
  near(a[1], 0.6); near(b[1], 0.4);
  assert.equal(ctx.sourceNodes.length, 2);
  player.stop(); assert.equal(player.playing, false);
});

test('new results, explicit stop and close invalidate an asynchronous audio resume', async () => {
  for (const action of ['stop', 'clear', 'close', 'replace']) {
    const ctx = context(); let resume;
    ctx.resume = () => new Promise(resolve => { resume = resolve; });
    const {player, a, b} = setup(ctx), pending = player.play();
    if (action === 'replace') player.setPair(a, receipt(a), b, receipt(b));
    else await player[action]();
    resume(); assert.equal(await pending, false); assert.equal(ctx.sourceNodes.length, 0);
  }
});

test('unequal windows stop at the common horizon and invalid replacements preserve the accepted pair', async () => {
  const {player, ctx, a, b} = setup(), short = wave(1000, 2400);
  assert.throws(() => player.setPair(a, receipt(a), b, {...receipt(b), fullScalePa:100}), /same native/);
  near(player.duration, 0.1);
  const stats = player.setPair(a, receipt(a), short, receipt(short));
  assert.equal(stats.truncated, true); near(player.duration, 0.05);
  await player.play();
  assert.ok(ctx.sourceNodes.every(source => source.buffer.data.length === 2400));
  assert.deepEqual(ctx.sourceNodes[0].starts[0], [10.02, 0, 0.05]);
  await assert.rejects(player.play(0.05), /position/);
  assert.throws(() => player.setVolume(1.01), /volume/);
  assert.throws(() => player.select('other'), /side/);
  await player.close(); assert.equal(ctx.state, 'closed');
});
