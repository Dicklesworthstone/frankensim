import {test} from 'node:test';
import assert from 'node:assert/strict';
import {SCHEMA, parseCSV, validateSnapshot, validateInputs, parseNative, toCSV, audioRate} from '../protocol.mjs';
export const inputs = {voltage: 0.2, duration: 0.1, resistivity: 2, support: 'tension', strain: 0.0001, limit: 0.01};
export const csv = 'time_s,cold_pressure_pa,heated_pressure_pa\n' + Array.from({length:800}, (_, i) => `${i/16000},${i%2*0.01},${i%3*0.02}`).join('\n');
export const log = 'support=None; tension_cold_N=10; tension_heated_N=10; heated_layer_stresses_Pa=None\nT_layers_K=[301, 302]; EA_cold_N=150000; EA_heated_N=149000; electrical_J=1; stored_heat_J=0.99999999; thermal_residual_J=1e-8\n';
test('native adapter consumes the current CLI layout without inventing history', () => {
  const s = parseNative(csv, log, inputs);
  assert.equal(s.schema, SCHEMA); assert.equal(audioRate(s),16000);
  assert.deepEqual(s.section.hot_temperature_k, [301,302]);
  assert.equal(s.metrics.hot_tension_n,10); assert.equal(s.traces.series[1].values[2],0.04);
  assert.equal(s.section.hot_stress_pa, undefined);
});
test('native changed output and missing evidence fail closed', () => {
  assert.throws(()=>parseNative(csv, log.replace('electrical_J','other'), inputs), /electrical_J/);
  assert.throws(()=>parseNative(csv, log, {...inputs, support:'strain'}), /stresses/);
  assert.throws(()=>parseNative(csv.replace('cold_pressure_pa','unknown_pa'),log,inputs), /layout/);
});
test('native constrained stresses remain signed', () => {
  const s = parseNative(csv,log.replace('None\n','Some([-2e6, 3e6])\n'), {...inputs, support:'strain'});
  assert.deepEqual(s.section.hot_stress_pa,[-2e6,3e6]);
});
test('CSV imports quoted names, BOM and CRLF and round trips exact values', () => {
  const s = parseCSV('\uFEFFtime_s,"sensor, A_pa"\r\n0,1e-6\r\n0.01,-2.0\r\n');
  assert.equal(s.traces.series[0].unit,'Pa');
  assert.deepEqual(parseCSV(toCSV(s)).traces.series[0].values,[1e-6,-2]);
});
for (const [name, value] of Object.entries({blank:'time_s,x\n0,\n1,2',duplicate:'time_s,x\n0,1\n0,2',
  infinity:'time_s,x\n0,1e999\n1,2',decreasing:'time_s,x\n1,1\n0,2',columns:'time_s,x\n0,1,2\n1,2',
  quotes:'time_s,"x\n0,1',double_header:'time_s,x,x\n0,1,2\n1,2,3',hex:'time_s,x\n0,0xff\n1,2'})) {
  test(`reject malformed CSV: ${name}`,()=>assert.throws(()=>parseCSV(value)));
}
test('finite sample budgets, units and source notes are mandatory',()=> {
  const s = parseCSV('time_s,sensor\n0,1\n1,2');
  assert.equal(s.traces.series[0].unit,'unspecified'); assert.equal(audioRate(s),null);
  assert.throws(()=>validateSnapshot({...s, source:{...s.source,notes:null}}));
  s.traces.series[0].values[1]=NaN; assert.throws(()=>validateSnapshot(s));
});
test('irregular pressure data is not silently resampled for playback',()=> {
  assert.equal(audioRate(parseCSV('time_s,x_pa\n0,1\n0.0001,2\n0.00021,3')),null);
});
test('request parameters are finite, bounded, and never shell arguments',()=> {
  assert.deepEqual(validateInputs(inputs),inputs);
  for (const change of [{voltage:'; touch x'}, {duration:Infinity}, {support:'none'}, {strain:0.1},
    {limit:1e-9}, {resistivity:0}, {command:'ls'}]) assert.throws(()=>validateInputs({...inputs,...change}));
});
test('geometry and source text survive as data, not markup',()=> {
  const s = parseNative(csv,log,inputs); s.title='<img src=x onerror=alert(1)>';
  assert.equal(validateSnapshot(s).title,s.title);
  s.section.radii_m=[0.5,0.1]; assert.throws(()=>validateSnapshot(s));
});
