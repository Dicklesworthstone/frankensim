import {MAX_BYTES, validateSnapshot as validateSignalSnapshot, validateInputs, parseCSV, toCSV, audioRate} from './protocol.mjs';
import {SectionView, drawPlot, format, palette, nearest, spectrum} from './visuals.mjs';
import {MESH_SCHEMA, MESH_LIMITS, validateMeshSnapshot} from './mesh-data.mjs';
import {parseVTU} from './import-vtu.mjs';
import {MeshWorkspace} from './mesh-ui.mjs';
const validateSnapshot = raw => raw?.schema === MESH_SCHEMA ? validateMeshSnapshot(raw) : validateSignalSnapshot(raw);
const $ = id => document.getElementById(id);
const element = (tag, text, className) => {const n = document.createElement(tag); if (text != null) n.textContent = text; if (className) n.className = className; return n;};
const BUILD = 'cargo build --release -p fs-couple --example heated_section_sound';
const STORE = 'frankensim.observatory.runs.v1';
let engine = {ready: false, material: 'synthetic'}, job = null, runs = [], current = null, baseline = null;
let state = 'hot', layer = 0, group = 'Pa', hidden = new Set(), cursorIndex = 0, plotLayout = null;
let audio = null, audioSource = null, audioGeneration = 0, noticeTimer, rendering = false;
const fftCache = new WeakMap();
const scene = new SectionView($('scene'), selected => {layer = selected; updateView();});
const meshWorkspace = new MeshWorkspace({notify, onBack: () => {current = null; updateView();}});
const activeView = () => meshWorkspace.active ? meshWorkspace.view : scene;
$('file-input').accept += ',.vtu';
$('drop-zone').querySelector('span').textContent = 'VTU or mesh JSON · up to 32 MiB / signal CSV or JSON · up to 8 MiB';
function notify(message, error = false) {
  clearTimeout(noticeTimer); $('notification').textContent = message; $('notification').hidden = false;
  $('notification').classList.toggle('error', error); $('activity').textContent = message;
  noticeTimer = setTimeout(() => {$('notification').hidden = true;}, error ? 14000 : 5000);
}
function download(data, type, name) {
  const url = URL.createObjectURL(data instanceof Blob ? data : new Blob([data], {type}));
  const a = element('a'); a.href = url; a.download = name; document.body.append(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 30000);
}
function nameFor(extension) {return (current?.snapshot.title || 'frankensim').replace(/[^a-z0-9_-]+/gi, '-').slice(0, 70) + '.' + extension;}
function values() {
  return validateInputs(Object.fromEntries(['voltage', 'duration', 'resistivity', 'support', 'strain', 'limit']
    .map(k => [k, k === 'support' ? $(k).value : Number($(k).value)])));
}
function draft() {
  $('strain-fields').hidden = $('support').value !== 'strain';
  let changed = false;
  try {changed = !current?.snapshot.inputs || JSON.stringify(values()) !== JSON.stringify(validateInputs(current.snapshot.inputs));} catch {changed = true;}
  $('draft-label').textContent = changed && current ? 'Edited parameters · not yet simulated' : current?.snapshot.inputs ? 'Controls match this run' : 'Ready to define an experiment';
  $('draft-label').classList.toggle('changed', changed && !!current);
}
function setInputs(v) {
  if (!v) return;
  for (const k of ['voltage', 'duration', 'resistivity', 'support', 'strain', 'limit']) {
    $(k).value = v[k]; if ($(k + '-range')) $(k + '-range').value = v[k];
  }
  if (engine.material === 'packs') {$('resistivity').value = '1'; $('resistivity-range').value = '1';}
  draft();
}
async function status(announce = false) {
  try {
    const res = await fetch('/api/status'); if (!res.ok) throw new Error('Local companion is unavailable.');
    engine = await res.json();
    $('engine-status').classList.toggle('ready', engine.ready);
    $('engine-label').textContent = engine.ready ? 'Native engine ready' : 'Local viewer mode';
    $('setup').hidden = engine.ready; $('run').disabled = !engine.ready || !!job;
    const packs = engine.material === 'packs';
    $('material-source').textContent = packs ? 'Your canonical material packs' : 'Synthetic demonstration cards';
    $('resistivity').disabled = packs; $('resistivity-range').disabled = packs;
    if (packs) {$('resistivity').value = '1'; $('resistivity-range').value = '1';}
    if (announce) notify(engine.ready ? 'Native engine found. Ready to run.' : 'Viewer is ready. Build the native example to enable local experiments.');
  } catch (e) {
    engine.ready = false; $('engine-label').textContent = 'Companion disconnected'; $('engine-status').classList.remove('ready');
    $('setup').hidden = false; $('run').disabled = true; if (announce) notify(e.message, true);
  }
  draft();
}
function persist() {
  try {
    const saved = []; let bytes = 2;
    for (const r of runs) {
      const json = JSON.stringify({snapshot: r.snapshot, added: r.added});
      if (bytes + json.length + 1 > 3 * 1024 * 1024) continue;
      saved.push(json); bytes += json.length + 1;
    }
    localStorage.setItem(STORE, '[' + saved.join(',') + ']');
    if (saved.length < runs.length) notify('Browser storage is bounded. Export large or older runs for durable retention.');
  } catch {notify('Browser storage is unavailable or full. Your run is still open; export it to keep a copy.');}
}
function restore() {
  try {
    const raw = localStorage.getItem(STORE); if (!raw || raw.length > 4 * 1024 * 1024) return;
    const saved = JSON.parse(raw); if (!Array.isArray(saved)) return;
    for (const r of saved.slice(0, 8)) {
      try {runs.push({id: crypto.randomUUID(), snapshot: validateSnapshot(r.snapshot), origin: 'Restored',
        added: typeof r.added === 'string' ? r.added.slice(0, 50) : ''});} catch { /* A bad entry cannot poison valid runs. */ }
    }
    if (runs.length) selectRun(runs[0]);
  } catch { /* Storage is optional; native runs and imports remain available. */ }
}
function addRun(snapshot, origin) {
  const r = {id: crypto.randomUUID(), snapshot: validateSnapshot(snapshot), origin, added: new Date().toLocaleTimeString([], {hour:'2-digit', minute:'2-digit'})};
  runs.unshift(r); runs = runs.slice(0, 8); selectRun(r); persist();
}
function selectRun(run) {
  if (job) return notify('Cancel the running experiment before switching runs.');
  stopAudio(); current = run; hidden = new Set(); cursorIndex = 0; $('cursor').value = 0; $('zoom').value = '1';
  group = run.snapshot.traces ? keyFor(run.snapshot.traces.series[0]) : 'mesh'; setInputs(run.snapshot.inputs); updateView();
}
function library() {
  const list = $('run-list'); if (!runs.length) return;
  list.replaceChildren();
  for (const r of runs) {
    const b = element('button', null, 'run-card' + (r === current ? ' active' : ''));
    b.setAttribute('aria-pressed', String(r === current)); b.disabled = !!job;
    b.append(element('strong', r.snapshot.title), element('span', `${r.origin.toUpperCase()} · ${r.added}`),
      element('span', r.snapshot.inputs ? `${r.snapshot.inputs.voltage} V · ${r.snapshot.inputs.duration} s · ${r.snapshot.inputs.support === 'strain' ? 'fixed ends' : '10 N'}` : r.snapshot.mesh ? `${r.snapshot.mesh.cell_count} cells · ${r.snapshot.mesh.time_s.length} recorded frames` : `${r.snapshot.traces.series.length} channels · imported values`, 'run-kind'));
    b.addEventListener('click', () => selectRun(r)); list.append(b);
  }
  $('baseline-label').textContent = baseline ? baseline.snapshot.title : 'No baseline selected';
  $('pin-baseline').disabled = !current || current === baseline; $('clear-baseline').hidden = !baseline;
}
function metric(id, value, unit) {
  const n = $(id); n.replaceChildren(document.createTextNode(format(value, 5))); if (value != null) n.append(element('em', unit));
}
function updateView() {
  const s = current?.snapshot, sec = s?.section, m = s?.metrics;
  if (s?.mesh) {meshWorkspace.enter(current, baseline); library(); return;}
  meshWorkspace.exit();
  $('cold').classList.toggle('selected', state === 'cold'); $('cold').setAttribute('aria-pressed', String(state === 'cold'));
  $('hot').classList.toggle('selected', state === 'hot'); $('hot').setAttribute('aria-pressed', String(state === 'hot'));
  $('core').setAttribute('aria-pressed', String(layer === 0)); $('sleeve').setAttribute('aria-pressed', String(layer === 1));
  const temperatures = sec?.[state + '_temperature_k'], stresses = state === 'hot' ? sec?.hot_stress_pa : null;
  const fieldValues = $('field').value === 'stress' ? stresses : temperatures, unit = $('field').value === 'stress' ? 'Pa' : 'K';
  $('core-value').textContent = fieldValues ? `${format(fieldValues[0], 5)} ${unit}` : '—';
  $('sleeve-value').textContent = fieldValues ? `${format(fieldValues[1], 5)} ${unit}` : '—';
  $('scene-state').textContent = sec ? `${current.origin.toUpperCase()} · ${state.toUpperCase()} PHYSICAL STATE` : s ? 'SIGNAL-ONLY IMPORT · NO GEOMETRY' : 'GEOMETRY PREVIEW · NO RUN DATA';
  $('scene-heading').textContent = sec || !s ? 'Concentric section' : 'No section geometry supplied';
  $('scene-caption').textContent = sec ? 'UNIFORM FIELD PER LAYER · LIGHTING IS VISUAL ONLY' : 'ONE SPECIMEN. TWO MATERIAL LAYERS.';
  $('scene-subtitle').textContent = sec ? ($('field').value === 'stress' && !stresses ? 'Stress was not reported for this state. No values are inferred.' : 'Drag to orbit · click a layer to inspect · scroll to zoom')
    : s ? 'The reference diagram is not a reconstruction of this imported signal.' : 'Run a native experiment or import a snapshot to reveal its state.';
  scene.update({snapshot: s, state, field: $('field').value, selected: layer, origin: current?.origin, baseline: baseline?.snapshot});
  metric('metric-temperature', temperatures?.[0], 'K');
  metric('metric-tension', m?.[state + '_tension_n'], 'N');
  metric('metric-heat', m?.stored_heat_j, 'J'); metric('metric-residual', m?.thermal_residual_j, 'J');
  $('metric-temperature-note').textContent = temperatures ? `${format(temperatures[0] - sec.cold_temperature_k[0], 4)} K relative to cold state` : 'Not supplied by this data';
  $('metric-heat-note').textContent = current?.origin === 'Local engine' ? 'Reported by the local native process' : 'Source-reported · not authenticated';
  $('metric-tension-note').textContent = m?.hot_tension_n != null && m?.cold_tension_n != null ? `${format(m.hot_tension_n - m.cold_tension_n, 4)} N cold → heated change` : 'Actual support resultant';
  $('selected-layer').textContent = layer ? 'Sleeve' : 'Core';
  const dl = $('layer-details'); dl.replaceChildren();
  const rows = sec ? [['Outer radius', `${format(sec.radii_m[layer] * 1000, 5)} mm`], ['Temperature', `${format(temperatures?.[layer], 5)} K`],
    ['Layer stress', stresses ? `${format(stresses[layer], 5)} Pa` : 'Not reported for this state'],
    ['Section stiffness (EA)', m?.[state + '_axial_n'] != null ? `${format(m[state + '_axial_n'], 5)} N` : 'Not reported'], ['Material basis', s.source.material]]
    : [['Outer radius', s ? 'Not provided' : `${layer ? '0.500' : '0.250'} mm`], ['Material basis', 'Reference diagram only']];
  for (const [key, value] of rows) dl.append(element('dt', key), element('dd', value));
  for (const id of ['export-json', 'export-csv', 'pin-baseline']) $(id).disabled = !s;
  $('empty-chart').hidden = !!s; $('sample-count').textContent = `${s ? s.traces.time_s.length.toLocaleString() : 0} samples`;
  const rate = s && audioRate(s);
  $('listen-cold').disabled = !rate; $('listen-hot').disabled = !s || !s.traces.series[1] || !audioRate(s, s.traces.series[1]);
  $('listen-cold').textContent = '▷ ' + (s?.traces.series[0]?.label || 'Cold').replace(/ pressure$/i, '');
  $('listen-hot').textContent = '▷ ' + (s?.traces.series[1]?.label || 'Heated').replace(/ pressure$/i, '');
  $('audio-note').textContent = rate ? `${format(s.traces.time_s.length / rate * 1000, 4)} ms of original samples at ${rate.toLocaleString()} Hz. Played once, at the original rate.`
    : 'Audio requires a Pa channel sampled uniformly at an audio rate. No resampling or inferred sound.';
  tabs(); library(); draft(); scheduleCharts();
}
const keyFor = s => s.unit === 'unspecified' ? `unknown:${s.id}` : `unit:${s.unit}`;
function tabs() {
  const box = $('plot-tabs'); box.replaceChildren(); if (!current) return;
  const groups = new Map(current.snapshot.traces.series.map(s => [keyFor(s), s.unit === 'unspecified' ? s.label + ' · unspecified units' : s.unit]));
  if (!groups.has(group)) group = groups.keys().next().value;
  for (const [key, label] of groups) {
    const b = element('button', label, key === group ? 'selected' : ''); b.setAttribute('aria-pressed', String(key === group));
    b.addEventListener('click', () => {group = key; tabs(); scheduleCharts();}); box.append(b);
  }
}
function scheduleCharts() {
  if (rendering || document.hidden) return;
  rendering = true; requestAnimationFrame(() => {rendering = false; charts();});
}
function spectra(s) {
  if (!fftCache.has(s)) fftCache.set(s, s.traces.series.map(channel => {
    const rate = audioRate(s, channel); return rate ? spectrum(channel.values, rate) : null;
  }));
  return fftCache.get(s);
}
function charts() {
  if (meshWorkspace.active) {meshWorkspace.drawProbe(); return;}
  const s = current?.snapshot, colors = palette(), box = $('channels'); box.replaceChildren();
  if (!s) {drawPlot($('trace'), []); drawPlot($('spectrum'), [], {xunit:'Hz', positive:true}); return;}
  cursorIndex = Math.min(cursorIndex, s.traces.time_s.length - 1);
  const t = s.traces.time_s, span = (t.at(-1) - t[0]) / Number($('zoom').value), at = t[cursorIndex];
  const begin = Math.max(t[0], Math.min(at - span / 2, t.at(-1) - span)), range = [begin, begin + span];
  const lines = [], frequencyLines = [];
  let yunit = '';
  s.traces.series.forEach((channel, i) => {
    if (keyFor(channel) !== group) return;
    yunit = channel.unit;
    const b = element('button'); b.setAttribute('aria-pressed', String(!hidden.has(i)));
    const dot = element('i', null, 'dot channel-' + i % 6);
    b.append(dot, document.createTextNode(channel.label + ' '), element('output', `${format(channel.values[cursorIndex], 5)} ${channel.unit}`));
    b.addEventListener('click', () => {if (hidden.has(i)) hidden.delete(i); else hidden.add(i); scheduleCharts();}); box.append(b);
    if (hidden.has(i)) return;
    lines.push({x: t, y: channel.values, color: colors[i % colors.length]});
    const sp = spectra(s)[i]; if (sp) frequencyLines.push({...sp, color: colors[i % colors.length]});
  });
  if (baseline?.snapshot.traces && baseline !== current && group !== 'unspecified') {
    const bs = baseline.snapshot, index = bs.traces.series.findLastIndex(ch => ch.unit === yunit && yunit !== 'unspecified');
    if (index >= 0) {
      const ch = bs.traces.series[index]; lines.push({x: bs.traces.time_s, y: ch.values, color: colors[2], dashed:true});
      const sp = spectra(bs)[index]; if (sp) frequencyLines.push({...sp, color:colors[2], dashed:true});
      const covered = at >= bs.traces.time_s[0] && at <= bs.traces.time_s.at(-1);
      box.append(element('span', `┄ Baseline: ${ch.label} · ${covered ? format(ch.values[nearest(bs.traces.time_s, at)], 5) + ' ' + ch.unit : 'outside recorded time range'}`, 'small muted'));
    }
  }
  plotLayout = drawPlot($('trace'), lines, {range, cursor: at, xunit: t.at(-1) < 1 ? 'ms' : 's', yunit});
  $('cursor-value').textContent = `${format(at, 5)} s`;
  const maxFrequency = frequencyLines.length ? Math.max(...frequencyLines.map(sp => sp.x.at(-1))) : 1;
  drawPlot($('spectrum'), frequencyLines, {range:[0, maxFrequency], xunit:'Hz', yunit:frequencyLines.length ? yunit : '', positive:true});
  $('spectrum-note').textContent = frequencyLines.length ? 'Full-record, single-sided amplitude · Hann coherent-gain correction · no normalization or inferred fundamental.'
    : 'Spectrum unavailable: this selection is not uniformly sampled audio-rate pressure data.';
}
async function runExperiment(e) {
  e?.preventDefault(); if (job) return;
  if (meshWorkspace.active) return notify('Choose “New native experiment” to return to the solver controls.');
  if (!$('experiment-form').reportValidity()) return;
  let input; try {input = values();} catch (error) {return notify(error.message, true);}
  await status(); if (!engine.ready) return notify('Build the native example, then click the engine indicator to refresh.', true);
  if (job) return;
  stopAudio(); const controller = new AbortController(); job = controller;
  $('run').disabled = true; $('cancel').hidden = false; $('parameters').disabled = true; document.body.classList.add('busy'); library();
  const start = performance.now(); const ticker = setInterval(() => {$('activity').textContent = `Native solver running · ${((performance.now() - start) / 1000).toFixed(1)} s elapsed · no partial result published`;}, 250);
  const timeout = setTimeout(() => controller.abort('client-timeout'), 65000);
  try {
    const res = await fetch('/api/run', {method:'POST', headers:{'Content-Type':'application/json', 'x-observatory-token':engine.token}, body:JSON.stringify(input), signal:controller.signal});
    const data = await res.json(); if (!res.ok) throw new Error(data.error || 'Native execution failed.');
    const snapshot = validateSnapshot(data.snapshot); job = null; addRun(snapshot, 'Local engine');
    notify(`Accepted native result · ${((performance.now() - start) / 1000).toFixed(2)} s wall time. Previous runs are retained.`);
  } catch (error) {
    notify(controller.signal.aborted ? (controller.signal.reason === 'client-timeout' ? 'Execution timed out. Previous result retained.' : 'Experiment cancelled. Previous result retained.') : error.message, true);
  } finally {
    clearTimeout(timeout); clearInterval(ticker); job = null; $('cancel').hidden = true; $('parameters').disabled = false;
    document.body.classList.remove('busy'); await status(); library();
  }
}
function stopAudio() {
  audioGeneration++;
  if (audioSource) {try {audioSource.stop();} catch {} audioSource.disconnect(); audioSource = null;}
  $('stop-audio').disabled = true;
}
async function listen(index) {
  const s = current?.snapshot, channel = s?.traces?.series[index], rate = s && channel && audioRate(s, channel);
  if (!rate) return;
  const scale = Number($('full-scale').value);
  if (!Number.isFinite(scale) || scale < 1e-12 || scale > 1e12) return notify('Shared full scale must be between 1e-12 and 1e12 Pa.', true);
  try {
    stopAudio(); const generation = audioGeneration; audio ||= new AudioContext(); await audio.resume();
    if (generation !== audioGeneration || current?.snapshot !== s) return;
    const buffer = audio.createBuffer(1, channel.values.length, rate), out = buffer.getChannelData(0); let clipped = 0;
    for (let i = 0; i < out.length; i++) {const x = channel.values[i] / scale; if (Math.abs(x) > 1) clipped++; out[i] = Math.max(-1, Math.min(1, x));}
    const source = audio.createBufferSource(); source.buffer = buffer; source.connect(audio.destination); audioSource = source;
    source.onended = () => {source.disconnect(); if (audioSource === source) {audioSource = null; $('stop-audio').disabled = true;}};
    $('stop-audio').disabled = false; source.start();
    $('audio-level').textContent = `Shared full scale ${format(scale)} Pa · ${clipped} clipped samples · no normalization`;
    if (clipped) notify(`${clipped} samples clipped at the selected listening gain. Increase shared full scale to avoid clipping.`, true);
  } catch (error) {stopAudio(); notify(`Audio unavailable: ${error.message}`, true);}
}
function modal(title) {
  $('dialog-title').textContent = title; $('dialog-content').replaceChildren();
  if (!$('dialog').open) $('dialog').showModal(); return $('dialog-content');
}
function help() {
  const box = modal('Your local physics workspace');
  box.append(element('p', 'Run the real heated-section example, or import existing results. Everything runs on this machine; no account, remote service or dependency download is needed to open the viewer.'));
  box.append(element('h3', 'Enable the native experiment'), element('code', BUILD), element('p', 'Use the normal repository setup with pinned sibling crates. Then refresh the engine indicator. Parameters marked as edited are drafts, not new simulation results.'));
  box.append(element('h3', 'Bring your own data'), element('p', 'Drop a snapshot JSON or numeric CSV up to 8 MiB. CSV starts with time_s followed by 1–12 unique channel names; suffixes such as _pa, _k and _n declare units. Times must strictly increase. Missing values and non-finite data are rejected, never filled in.'));
  box.append(element('h3', 'Explore mesh fields'), element('p', 'Import ASCII VTU from fs-viz or a recorded mesh JSON up to 32 MiB. VTU import requires you to declare coordinate units. Pick points or cells, inspect vector components, and replay recorded frames. Color limits span the whole record; no intermediate physics is synthesized. Use Shift-drag to pan the mesh and the frame controls to inspect time.'));
  box.append(element('h3', 'Move quickly'));
  const table = element('table', null, 'evidence-table');
  for (const [key, action] of [['R', 'Run experiment'], ['I', 'Import a run'], ['⌘ / Ctrl K', 'Command palette'], ['?', 'Open this guide'], ['Esc', 'Close a dialog or cancel a running experiment'], ['Arrow keys on model', 'Rotate camera; + / − zoom; 0 reset'], ['Arrow keys on trace', 'Inspect adjacent samples']]) {
    const row = element('tr'); row.append(element('th', key), element('td', action)); table.append(row);
  }
  box.append(table, element('p', 'The viewport shows two uniform-temperature layers, not a continuous spatial solution. Cold and heated audio are separate, short frozen-state recordings. All lighting, the ground grid and the exaggerated transverse scale are visual aids.'));
}
function evidence() {
  const box = modal('Evidence, not decoration'), s = current?.snapshot;
  if (!s) {box.append(element('p', 'No accepted result is loaded. The model is a geometry preview; the dashboard deliberately leaves numerical readouts empty.')); return;}
  box.append(element('p', `${current.origin}: ${s.source.model}`), element('p', s.source.material));
  if (current.origin !== 'Local engine') box.append(element('p', 'This data was imported or restored. Its stated source is not authenticated, even when it names a native solver.'));
  if (s.inputs) {
    box.append(element('h3', 'Inputs attached to this result')); const t = element('table', null, 'evidence-table');
    for (const [k, v] of Object.entries(s.inputs)) {const r = element('tr'); r.append(element('th', k), element('td', String(v))); t.append(r);} box.append(t);
  }
  box.append(element('h3', 'Model boundaries'));
  s.source.notes.forEach(note => box.append(element('p', note)));
  box.append(element('p', s.mesh ? `${s.mesh.positions_m.length} source points, ${s.mesh.cell_count} source cells, ${s.mesh.time_s.length} recorded frames. Mesh exports retain all source values. Probe CSV contains the selected entity and component. Display interpolation never changes source samples.` : `${s.traces.time_s.length} original samples per channel. Chart envelopes retain pixel-bin extrema; CSV/JSON exports retain every sample. Spectra use the complete record, not the zoomed view.`));
}
function commands() {
  const box = modal('Jump to an action'), input = element('input'), list = element('div');
  input.placeholder = 'Search actions…'; input.setAttribute('aria-label', 'Search workspace actions'); box.append(input, list);
  const actions = [['Run experiment', 'R', () => $('experiment-form').requestSubmit()], ['Import snapshot or CSV', 'I', () => $('file-input').click()],
    ['Inspect evidence', '', evidence], ['Export current snapshot', '', exportJSON], ['Pin current run as baseline', '', pin],
    ['Reset camera', '0', () => activeView().reset()], ['Switch color theme', '', theme], ['Help and shortcuts', '?', help]];
  function populate() {
    list.replaceChildren();
    for (const [label, key, action] of actions.filter(a => a[0].toLowerCase().includes(input.value.toLowerCase()))) {
      const b = element('button', null, 'command-item'); b.append(element('span', label), element('kbd', key || '↵'));
      b.addEventListener('click', () => {$('dialog').close(); action();}); list.append(b);
    }
    if (!list.children.length) list.append(element('p', 'No matching actions.'));
  }
  input.addEventListener('input', populate); input.addEventListener('keydown', e => {
    if (e.key === 'Enter') {e.preventDefault(); list.querySelector('button')?.click();}
    if (e.key === 'ArrowDown') {e.preventDefault(); list.querySelector('button')?.focus();}
  }); populate(); input.focus();
}
function theme() {
  const value = document.documentElement.dataset.theme === 'dark' ? 'light' : 'dark'; document.documentElement.dataset.theme = value;
  try {localStorage.setItem('frankensim.observatory.theme', value);} catch {}
  activeView().invalidate(); scheduleCharts();
}
function exportJSON() {
  if (!current) return notify('Load a run before exporting.');
  download(JSON.stringify(current.snapshot, null, 2) + '\n', 'application/json', nameFor('json'));
}
function pin() {if (!current) return notify('Load a run before pinning a baseline.'); baseline = current; updateView(); notify('Baseline pinned. Select or run another experiment to overlay its matching-unit signal.');}
async function importFile(file) {
  if (!file) return;
  if (job) return notify('Cancel the experiment before importing another run.');
  const csv = file.name.toLowerCase().endsWith('.csv'), vtu = file.name.toLowerCase().endsWith('.vtu');
  if (file.size > (csv ? MAX_BYTES : MESH_LIMITS.bytes)) return notify('File exceeds the import budget (8 MiB CSV; 32 MiB mesh).', true);
  try {
    const text = await file.text();
    let snapshot;
    if (vtu) {
      const lengthUnit = await coordinateUnit(); if (lengthUnit == null) return;
      snapshot = parseVTU(text, {lengthUnit, title: file.name.slice(0, 190)});
    } else if (csv) snapshot = parseCSV(text, file.name.slice(0, 190));
    else {
      const raw = JSON.parse(text);
      if (raw.schema !== MESH_SCHEMA && file.size > MAX_BYTES) throw new Error('Signal snapshots remain limited to 8 MiB.');
      snapshot = validateSnapshot(raw);
    }
    if (job) throw new Error('An experiment started during import. Cancel it before importing.');
    addRun(snapshot, 'Imported'); notify(`Imported ${file.name}. Source claims remain unverified.`);
  } catch (e) {notify(`Import refused: ${e.message} Previous result retained.`, true);}
}
function coordinateUnit() {
  return new Promise(resolve => {
    const box = modal('Declare the mesh coordinate unit');
    box.append(element('p', 'VTU coordinates do not necessarily declare units. Choose the unit used by the exporter. Field Unit attributes are retained unchanged. No material properties are inferred.'));
    const label = element('label', 'Coordinate unit'); label.htmlFor = 'import-length-unit';
    const select = element('select'); select.id = 'import-length-unit';
    for (const [value, text] of [['', 'Choose explicitly…'], ['m', 'Metres (m)'], ['mm', 'Millimetres (mm)'], ['cm', 'Centimetres (cm)']]) {
      const option = element('option', text); option.value = value; select.append(option);
    }
    const button = element('button', 'Import mesh', 'button primary'); button.id = 'confirm-mesh-import'; button.disabled = true;
    select.addEventListener('change', () => {button.disabled = !select.value;});
    const close = () => {resolve(null);}; $('dialog').addEventListener('close', close, {once: true});
    button.addEventListener('click', () => {$('dialog').removeEventListener('close', close); resolve(select.value); $('dialog').close();});
    box.append(label, select, button); select.focus();
  });
}
for (const id of ['voltage', 'duration', 'resistivity']) {
  $(id + '-range').addEventListener('input', () => {$(id).value = $(id + '-range').value; draft();});
  $(id).addEventListener('input', () => {$(id + '-range').value = $(id).value; draft();});
}
for (const id of ['support', 'strain', 'limit']) $(id).addEventListener('input', draft);
$('experiment-form').addEventListener('submit', runExperiment);
$('cancel').addEventListener('click', () => job?.abort());
$('engine-status').addEventListener('click', () => status(true));
for (const id of ['import', 'import-nav']) $(id).addEventListener('click', () => $('file-input').click());
$('file-input').addEventListener('change', e => {importFile(e.target.files[0]); e.target.value = '';});
for (const id of ['help', 'footer-help']) $(id).addEventListener('click', help);
for (const id of ['show-evidence', 'evidence-nav']) $(id).addEventListener('click', evidence);
$('command').addEventListener('click', commands); $('theme').addEventListener('click', theme);
$('close-dialog').addEventListener('click', () => $('dialog').close());
$('dialog').addEventListener('click', e => {if (e.target === $('dialog')) {const r = $('dialog').getBoundingClientRect(); if (e.clientX < r.left || e.clientX > r.right || e.clientY < r.top || e.clientY > r.bottom) $('dialog').close();}});
$('copy-build').addEventListener('click', async () => {try {await navigator.clipboard.writeText(BUILD); notify('Build command copied.');} catch {help();}});
$('library-nav').addEventListener('click', () => {$('runs-heading').scrollIntoView({block:'center', behavior:matchMedia('(prefers-reduced-motion: reduce)').matches ? 'instant' : 'smooth'}); $('run-list').querySelector('button')?.focus({preventScroll:true});});
$('studio-nav').addEventListener('click', () => {$('workspace').scrollIntoView({block:'start'}); $('workspace').focus({preventScroll:true});});
$('field').addEventListener('change', updateView);
for (const id of ['cold', 'hot']) $(id).addEventListener('click', () => {state = id; updateView();});
for (const [id, index] of [['core', 0], ['sleeve', 1]]) $(id).addEventListener('click', () => {layer = index; updateView();});
$('zoom-in').addEventListener('click', () => activeView().zoomBy(1.15)); $('zoom-out').addEventListener('click', () => activeView().zoomBy(1 / 1.15));
$('reset-view').addEventListener('click', () => activeView().reset());
$('snapshot-image').addEventListener('click', () => {if (meshWorkspace.active) return meshWorkspace.screenshot(download, nameFor('png')); scene.draw(); $('scene').toBlob(blob => {if (blob) download(blob, 'image/png', nameFor('png'));});});
$('export-json').addEventListener('click', exportJSON); $('export-csv').addEventListener('click', () => {if (current) download(toCSV(meshWorkspace.active ? {traces: meshWorkspace.trace()} : current.snapshot), 'text/csv', nameFor('csv'));});
$('pin-baseline').addEventListener('click', pin); $('clear-baseline').addEventListener('click', () => {baseline = null; updateView();});
$('listen-cold').addEventListener('click', () => listen(0)); $('listen-hot').addEventListener('click', () => listen(1)); $('stop-audio').addEventListener('click', stopAudio);
$('zoom').addEventListener('change', scheduleCharts);
function cursorTo(time) {
  if (!current) return;
  if (meshWorkspace.active) return meshWorkspace.cursor(time);
  const t = current.snapshot.traces.time_s; cursorIndex = nearest(t, time);
  $('cursor').value = String((t[cursorIndex] - t[0]) / (t.at(-1) - t[0]) * 1000); scheduleCharts();
}
$('cursor').addEventListener('input', () => {if (current && !meshWorkspace.active) {const t = current.snapshot.traces.time_s; cursorTo(t[0] + Number($('cursor').value) / 1000 * (t.at(-1) - t[0]));}});
$('trace').addEventListener('pointerdown', e => {const plot = meshWorkspace.active ? meshWorkspace.plot : plotLayout; if (plot) {const x = e.clientX - $('trace').getBoundingClientRect().left; cursorTo(plot.range[0] + (x - plot.start) / plot.width * (plot.range[1] - plot.range[0]));}});
$('trace').addEventListener('keydown', e => {if (current && ['ArrowLeft', 'ArrowRight'].includes(e.key)) {e.preventDefault(); if (meshWorkspace.active) return meshWorkspace.step(e.key === 'ArrowLeft' ? -1 : 1); cursorIndex = Math.max(0, Math.min(current.snapshot.traces.time_s.length - 1, cursorIndex + (e.key === 'ArrowLeft' ? -1 : 1))); cursorTo(current.snapshot.traces.time_s[cursorIndex]);}});
new ResizeObserver(scheduleCharts).observe($('trace')); new ResizeObserver(scheduleCharts).observe($('spectrum'));
document.addEventListener('visibilitychange', () => {if (!document.hidden) {activeView().invalidate(); scheduleCharts();}});
document.addEventListener('keydown', e => {
  if ((e.metaKey || e.ctrlKey) && e.key.toLowerCase() === 'k') {e.preventDefault(); if ($('dialog').open) $('dialog').close(); else commands(); return;}
  if (!$('dialog').open && e.key === 'Escape' && job) {job.abort(); return;}
  if ($('dialog').open || e.metaKey || e.ctrlKey || e.altKey || ['INPUT', 'SELECT', 'TEXTAREA'].includes(e.target.tagName)) return;
  if (e.key === 'Escape' && job) {job.abort(); return;}
  if (e.key.toLowerCase() === 'r') {e.preventDefault(); $('experiment-form').requestSubmit();}
  if (e.key.toLowerCase() === 'i') {e.preventDefault(); $('file-input').click();}
  if (e.key === '?') {e.preventDefault(); help();}
});
let dragDepth = 0;
document.addEventListener('dragenter', e => {if (e.dataTransfer.types.includes('Files')) {dragDepth++; $('drop-zone').hidden = false;}});
document.addEventListener('dragleave', () => {if (--dragDepth <= 0) {dragDepth = 0; $('drop-zone').hidden = true;}});
document.addEventListener('dragover', e => {if (e.dataTransfer.types.includes('Files')) {e.preventDefault(); e.dataTransfer.dropEffect = 'copy';}});
document.addEventListener('drop', e => {e.preventDefault(); dragDepth = 0; $('drop-zone').hidden = true; if (e.dataTransfer.files.length !== 1) notify('Import one run at a time.', true); else importFile(e.dataTransfer.files[0]);});
window.addEventListener('pagehide', () => {job?.abort(); meshWorkspace.pause(); stopAudio(); audio?.close();});
try {const saved = localStorage.getItem('frankensim.observatory.theme'); if (['dark', 'light'].includes(saved)) document.documentElement.dataset.theme = saved;} catch {}
restore(); updateView(); status();
