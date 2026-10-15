import {test} from 'node:test';
import assert from 'node:assert/strict';
import {spectrum, envelopeIndices, nearest, format} from '../visuals.mjs';
import {parseCSV, toCSV, validateSnapshot} from '../protocol.mjs';
test('Hann FFT preserves physical amplitude and frequency; silence stays zero', () => {
  const count = 1024, rate = 16000, amplitude = 0.003, bin = 32;
  const samples = Array.from({length:count}, (_, i) => amplitude * Math.sin(2*Math.PI*bin*i/count));
  const result = spectrum(samples, rate);
  const peak = result.y.indexOf(Math.max(...result.y));
  assert.equal(result.x[peak], 500);
  assert.ok(Math.abs(result.y[peak] - amplitude) < 1e-8);
  assert.equal(Math.max(...spectrum(samples.map(()=>0),rate).y),0);
  assert.equal(spectrum([1,2,3],rate),null);
});
test('FFT uses actual record length for coherent gain, including non-power-of-two records', () => {
  const result = spectrum(Array(800).fill(0.025),16000);
  assert.equal(result.x.length,513);
  assert.ok(Math.abs(result.y[0]-.025)<1e-14);
  assert.equal(result.x.at(-1),8000);
});
test('pixel envelope retains narrow positive and negative extrema in temporal order', () => {
  const x = Array.from({length:1000},(_,i)=>i), y = x.map(()=>0); y[7]=10; y[9]=-11; y[996]=20;
  const indices = envelopeIndices(x,y,0,999,20);
  assert.ok(indices.includes(7)); assert.ok(indices.includes(9)); assert.ok(indices.includes(996));
  assert.ok(indices.every((v,i)=>i===0||v>indices[i-1]));
  assert.ok(indices.length<=42);
  assert.ok(envelopeIndices(x,y,100,200,20).every(i=>i>=100&&i<=200));
});
test('cursor selection clamps safely to actual sample indices', () => {
  const x=[0,.001,.007,.01];
  assert.equal(nearest(x,-1),0); assert.equal(nearest(x,1),3); assert.equal(nearest(x,.005),2);
  assert.equal(nearest(x,.001),1); assert.equal(format(undefined),'—');
});
test('CSV exports disambiguate duplicate display labels without dropping or rounding samples', () => {
  const s=parseCSV('time_s,a_pa,b_pa\n0,1e-8,-2\n1,3,-4');
  s.traces.series.forEach(c=>c.label='identical');
  const result=parseCSV(toCSV(validateSnapshot(s)));
  assert.deepEqual(result.traces.series.map(c=>c.values),[[1e-8,3],[-2,-4]]);
  assert.deepEqual(result.traces.series.map(c=>c.unit),['Pa','Pa']);
  s.traces.series[0].label='x'.repeat(100);
  assert.deepEqual(parseCSV(toCSV(s)).traces.series[0].values,[1e-8,3]);
});
