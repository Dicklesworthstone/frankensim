// Shared, dependency-free boundary for the local process adapter and browser.
// A visualization snapshot is not an authenticated Design Ledger certificate.
export const SCHEMA = 'frankensim.observatory/v1';
export const MAX_BYTES = 8 * 1024 * 1024;
export const MAX_POINTS = 32768;
const finite = (x, label) => {
  if (typeof x !== 'number' || !Number.isFinite(x) || Math.abs(x) > 1e100)
    throw new Error(`${label} must be a finite number (magnitude <= 1e100).`);
  return x;
};
const text = (x, label, max = 200) => {
  if (typeof x !== 'string' || !x.trim() || x.length > max) throw new Error(`Invalid ${label}.`);
  return x;
};
const vector = (x, label, length) => {
  if (!Array.isArray(x) || x.length !== length) throw new Error(`${label}: inconsistent array length.`);
  return x.map((v) => finite(v, label));
};
export function validateInputs(raw) {
  if (!raw || typeof raw !== 'object' || Array.isArray(raw)) throw new Error('Expected experiment parameters.');
  const allowed = ['voltage', 'duration', 'resistivity', 'support', 'strain', 'limit'];
  if (Object.keys(raw).some(k => !allowed.includes(k))) throw new Error('Unknown experiment parameter.');
  const result = {};
  for (const [key, lo, hi] of [['voltage', 0, 1], ['duration', 0.001, 10],
    ['resistivity', 0.1, 10], ['strain', 1e-9, 0.01], ['limit', 1e-9, 0.1]]) {
    const v = finite(raw[key], key);
    if (v < lo || v > hi) throw new Error(`${key} must be between ${lo} and ${hi}.`);
    result[key] = v;
  }
  if (!['tension', 'strain'].includes(raw.support)) throw new Error('Unknown support prescription.');
  if (result.strain > result.limit) throw new Error('Strain exceeds the declared strain limit.');
  result.support = raw.support;
  return result;
}
export function validateSnapshot(raw) {
  if (!raw || raw.schema !== SCHEMA) throw new Error(`Expected schema ${SCHEMA}.`);
  const t = raw.traces?.time_s;
  if (!Array.isArray(t) || t.length < 2 || t.length > MAX_POINTS) throw new Error(`Expected 2–${MAX_POINTS} time samples.`);
  const time = t.map(v => finite(v, 'time_s'));
  if (time[0] < 0 || time.some((v, i) => i > 0 && v <= time[i - 1])) throw new Error('Times must be nonnegative and strictly increasing.');
  const list = raw.traces.series;
  if (!Array.isArray(list) || !list.length || list.length > 12) throw new Error('Expected 1–12 signal channels.');
  const series = list.map((s, i) => ({id: `s${i}`, label: text(s.label, 'channel label', 100),
    unit: text(s.unit, 'channel unit', 24), values: vector(s.values, 'channel', time.length)}));
  const source = raw.source;
  if (!source || !Array.isArray(source.notes) || source.notes.length > 12) throw new Error('Source notes are required.');
  const result = {schema: SCHEMA, title: text(raw.title, 'title'), source: {
    model: text(source.model, 'model'), material: text(source.material, 'material basis'),
    notes: source.notes.map(n => text(n, 'source note', 1500))}, traces: {time_s: time, series}};
  if (raw.inputs) result.inputs = validateInputs(raw.inputs);
  if (raw.section) {
    const s = raw.section;
    if (s.kind !== 'concentric-section') throw new Error('Unsupported geometry.');
    const length = finite(s.length_m, 'length_m');
    const radii = vector(s.radii_m, 'radii_m', 2);
    if (length <= 0 || radii[0] <= 0 || radii[1] <= radii[0]) throw new Error('Invalid section geometry.');
    const cold = vector(s.cold_temperature_k, 'cold temperature', 2);
    const hot = vector(s.hot_temperature_k, 'hot temperature', 2);
    if ([...cold, ...hot].some(x => x <= 0)) throw new Error('Absolute temperatures must be positive.');
    result.section = {kind: s.kind, length_m: length, radii_m: radii,
      cold_temperature_k: cold, hot_temperature_k: hot};
    if (s.hot_stress_pa != null) result.section.hot_stress_pa = vector(s.hot_stress_pa, 'layer stress', 2);
  }
  if (raw.metrics) {
    result.metrics = {};
    for (const key of ['cold_tension_n', 'hot_tension_n', 'cold_axial_n', 'hot_axial_n',
      'electrical_energy_j', 'stored_heat_j', 'thermal_residual_j']) {
      if (raw.metrics[key] !== undefined) result.metrics[key] = finite(raw.metrics[key], key);
    }
  }
  return result;
}
// Numeric CSV, including quoted headers/CRLF/BOM. No guessing of missing samples.
function rows(csv) {
  const out = []; let row = [], cell = '', quoted = false, closed = false;
  for (let i = 0; i <= csv.length; i++) {
    const c = i === csv.length ? '\n' : csv[i];
    if (quoted) {
      if (c === '"' && csv[i + 1] === '"') { cell += '"'; i++; }
      else if (c === '"') { quoted = false; closed = true; }
      else if (i === csv.length) throw new Error('Unclosed CSV quote.');
      else cell += c;
    } else if (c === '"' && !cell && !closed) quoted = true;
    else if (c === ',' || c === '\n' || c === '\r') {
      row.push(cell); cell = ''; closed = false;
      if (row.length > 13) throw new Error('At most 12 signal channels are supported.');
      if (c !== ',') {
        if (row.some(x => x !== '')) out.push(row);
        row = [];
        if (c === '\r' && csv[i + 1] === '\n') i++;
        if (out.length > MAX_POINTS + 1) throw new Error('CSV exceeds the sample budget.');
      }
    } else {
      if (closed || c === '"') throw new Error('Malformed CSV quoting.');
      cell += c;
    }
  }
  return out;
}
export function parseCSV(csv, title = 'Imported signal') {
  if (typeof csv !== 'string' || new TextEncoder().encode(csv).length > MAX_BYTES) throw new Error('CSV exceeds 8 MiB.');
  const data = rows(csv.replace(/^\uFEFF/, ''));
  const header = data.shift()?.map(x => x.trim());
  if (!header || header[0] !== 'time_s' || header.length < 2 || new Set(header).size !== header.length)
    throw new Error('CSV needs a unique header: time_s,channel_name,...');
  const numbers = data.map((row, i) => {
    if (row.length !== header.length) throw new Error(`CSV row ${i + 2}: missing or extra columns.`);
    return row.map(v => {
      if (!/^[+-]?(?:\d+\.?\d*|\.\d+)(?:e[+-]?\d+)?$/i.test(v.trim())) throw new Error(`CSV row ${i + 2}: expected a number.`);
      return finite(Number(v), 'CSV value');
    });
  });
  const units = {pa: 'Pa', k: 'K', n: 'N', j: 'J', m: 'm', w: 'W', hz: 'Hz'};
  return validateSnapshot({schema: SCHEMA, title, source: {model: 'Imported time series',
    material: 'Not provided by CSV', notes: ['Imported data; provenance and solver validity are not authenticated.',
      'Units are read from column suffixes; unspecified units remain unspecified. No section geometry is inferred.']},
    traces: {time_s: numbers.map(r => r[0]), series: header.slice(1).map((name, i) => ({
      label: name.replace(/_/g, ' '), unit: units[name.split('_').at(-1).toLowerCase()] || 'unspecified',
      values: numbers.map(r => r[i + 1])}))}});
}
export function toCSV(snapshot) {
  const escape = s => `"${(/^[=+@-]/.test(s) ? "'" + s : s).replaceAll('"', '""')}"`;
  return ['time_s,' + snapshot.traces.series.map(s => escape(`${s.label.replace(/ /g, '_')}_${s.unit.toLowerCase()}`)).join(','),
    ...snapshot.traces.time_s.map((t, i) => [t, ...snapshot.traces.series.map(s => s.values[i])].join(','))].join('\n') + '\n';
}
export function audioRate(snapshot, series = snapshot.traces.series[0]) {
  const t = snapshot.traces.time_s, dt = t[1] - t[0], rate = Math.round(1 / dt);
  if (series.unit !== 'Pa' || rate < 3000 || rate > 192000 || Math.abs(rate * dt - 1) > 1e-6) return null;
  return t.every((x, i) => Math.abs(x - t[0] - i / rate) <= Math.max(1e-9, dt * 1e-5)) ? rate : null;
}
export function parseNative(csv, log, inputs, packs = false) {
  const snapshot = parseCSV(csv, 'Heated concentric section');
  if (snapshot.traces.series.length !== 2 || snapshot.traces.series[0].label !== 'cold pressure pa'
    || snapshot.traces.series[1].label !== 'heated pressure pa' || snapshot.traces.time_s.length !== 800
    || audioRate(snapshot) !== 16000 || snapshot.traces.time_s[0] !== 0)
    throw new Error('Native output layout changed; update the heated_section_sound adapter.');
  const number = key => {
    const match = log.match(new RegExp(`(?:^|[;\\s])${key}=([^;\\s]+)`, 'm'));
    if (!match) throw new Error(`Native output is missing ${key}; rebuild the current example.`);
    return finite(Number(match[1]), key);
  };
  const temps = log.match(/T_layers_K=\[([^\]]+)\]/);
  if (!temps) throw new Error('Native output is missing layer temperatures.');
  const stresses = log.match(/heated_layer_stresses_Pa=Some\(\[([^\]]+)\]\)/);
  if (inputs.support === 'strain' && !stresses) throw new Error('Native output is missing fixed-support stresses.');
  snapshot.source = {model: 'fs-couple / heated_section_sound',
    material: packs ? 'User-supplied canonical material packs' : 'Synthetic demonstration cards · not measured materials',
    notes: ['Two frozen-state acoustic observations of the same specimen and pluck, not sound during heating.',
      'Reduced two-layer thermal network with uniform temperature per layer, not a resolved radial temperature field.',
      'One-way thermal-to-quasistatic mechanics on fixed geometry. No mechanical-work feedback, melting or buckling.',
      'Reported residual is thermal accounting, not a discretization-error bound or authenticated certificate.']};
  snapshot.inputs = inputs;
  snapshot.section = {kind: 'concentric-section', length_m: 0.65, radii_m: [0.00025, 0.0005],
    cold_temperature_k: [300, 300], hot_temperature_k: temps[1].split(',').map(Number),
    hot_stress_pa: stresses ? stresses[1].split(',').map(Number) : null};
  snapshot.metrics = Object.fromEntries(Object.entries({cold_tension_n: 'tension_cold_N', hot_tension_n: 'tension_heated_N',
    cold_axial_n: 'EA_cold_N', hot_axial_n: 'EA_heated_N', electrical_energy_j: 'electrical_J',
    stored_heat_j: 'stored_heat_J', thermal_residual_j: 'thermal_residual_J'}).map(([k, v]) => [k, number(v)]));
  snapshot.traces.series[0].label = 'Cold pressure'; snapshot.traces.series[1].label = 'Heated pressure';
  return validateSnapshot(snapshot);
}
