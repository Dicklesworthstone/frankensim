// One family: listen (simulation against recording), grades, grading history, qualification evidence, timeline and
// the screened reference recordings.
import { h, getJSON, fmt, badge, splitBar, link, chip, empty, DATA } from '../app.js';
import { Player, controls, stopAll, focusPlayer, loadBuffer } from '../audio.js';
import { mapFrame, computeMap } from '../spectro.js';
import { rulerBar, historyChart, bodyModes, bowedSweep, bandDecays } from '../charts.js';

const TABS = [
  ['listen', 'Listen'],
  ['grades', 'Grades'],
  ['history', 'History'],
  ['qualification', 'Qualification'],
  ['timeline', 'Timeline'],
  ['references', 'References'],
];

export async function render(root, r, sameFamily) {
  const doc = await getJSON(`family/${r.id}.json`);
  let page = root.querySelector(`.familypage[data-id="${CSS.escape(r.id)}"]`);
  if (!sameFamily || !page) {
    page = h('div.familypage', { 'data-id': r.id }, header(doc), h('nav.tabs'), h('div.tabbody'));
    root.replaceChildren(page);
  }
  const counts = {
    listen: doc.pairs.filter(p => p.kind !== 'reference-only').length + doc.extras.length + doc.fixtures.length,
    grades: doc.pairs.filter(p => p.kind === 'graded').length,
    history: doc.pairs.filter(p => p.history.length > 1).length,
    qualification: doc.milestones.length + doc.evidence.length,
    timeline: doc.timeline_total + doc.closes.length,
    references: doc.references.length,
  };
  page.querySelector('nav.tabs').replaceChildren(...TABS.map(([id, label]) =>
    h(`a${id === r.tab ? '.active' : ''}`, { href: link('f', doc.id, id) }, label, h('span.count', null, String(counts[id])))));
  stopAll();
  const body = page.querySelector('.tabbody');
  const view = { listen, grades, history, qualification, timeline, references }[r.tab] || listen;
  body.replaceChildren(view(doc, r.item));
}

function header(doc) {
  const st = doc.status || {};
  const v = doc.verdict;
  const physics = st.physics_present || st.physics_missing
    ? h('details.physics', null, h('summary', null, 'What the simulation includes, and what it still lacks'),
      h('div.cols', null,
        h('div', null, h('h4', null, 'Included'), h('ul', null, (st.physics_present || []).map(x => h('li', null, x)))),
        h('div', null, h('h4', null, 'Missing'), h('ul', null, (st.physics_missing || []).map(x => h('li', null, x))))))
    : null;
  return h('header.familyhead', null,
    h('a.back', { href: '#/' }, '← All instruments'),
    h('div.titlerow', null, h('h1', null, doc.title), badge(v.badge, v.label), h('span.section', null, doc.section)),
    st.reference_instrument ? h('p.small.muted', null, `Reference instrument: ${st.reference_instrument}`) : null,
    h(doc.summary_derived ? 'p.summary.muted' : 'p.summary', null, doc.summary),
    v.blocked ? h('p.blocked', null, `Blocked: ${typeof v.blocked === 'string' ? v.blocked : JSON.stringify(v.blocked)}`) : null,
    physics);
}

// ---------------------------------------------------------------- Listen

function listenItems(doc) {
  const items = [];
  const add = (group, key, title, sub, kind, data) => items.push({ group, key, title, sub, kind, data });
  for (const p of doc.pairs) {
    const sub = p.kind === 'graded' && p.tally
      ? `${p.tally.loose.inside}/${p.tally.loose.inside + p.tally.loose.outside} inside`
      : p.kind === 'refused' ? 'refused' : p.kind === 'listening' ? 'no recording' : 'recording only';
    const group = { graded: 'Graded against a recording', listening: 'Rendered, no recording to grade against', refused: 'Refused by the grader', 'reference-only': 'Recordings the grader holds, not rendered yet' }[p.kind];
    add(group, p.id, [p.note, p.dynamic].filter(Boolean).join(' '), sub, p.kind, p);
  }
  for (const x of doc.extras) add('Other renders (not graded)', `x-${x.id}`, x.title, fmt.date(x.date), 'extra', x);
  for (const x of doc.fixtures) add('August listening fixtures (before the 2026-10-21 pivot)', `fx-${x.id}`, x.title, fmt.date(x.date), 'fixture', x);
  return items;
}

function listen(doc, itemKey) {
  const items = listenItems(doc);
  if (!items.length) return empty('Nothing has been rendered for this family yet. Its screened reference recordings are under References.');
  const sel = items.find(i => i.key === itemKey) || items.find(i => i.kind === 'graded') || items[0];
  const groups = [...new Set(items.map(i => i.group))];
  const list = h('aside.itemlist', null, groups.map(g => h('div.itemgroup', null, h('h4', null, g),
    items.filter(i => i.group === g).map(i => h(`a.item${i === sel ? '.on' : ''}`, { href: link('f', doc.id, 'listen', i.key) },
      h('span.t', null, i.title), h(`span.s.k-${i.kind}`, null, i.sub))))));
  const detail = h('div.detail');
  const parts = sel.kind === 'extra' || sel.kind === 'fixture' ? clipDetail(sel.data, sel.kind) : pairDetail(doc, sel.data);
  detail.append(...parts.filter(Boolean));
  return h('div.listen', null, list, detail);
}

function pairDetail(doc, p) {
  const out = [];
  out.push(h('div.detailhead', null, h('h2', null, [p.note, p.dynamic].filter(Boolean).join(' ')), chip(kindLabel(p.kind), `k-${p.kind}`),
    p.duration_s ? h('span.small.muted', null, `${p.duration_s} s analysed`) : null));
  if (p.kind === 'graded' && p.tally) {
    const t = p.tally;
    out.push(h('p.verdictline', null, splitBar(t.loose.inside, t.loose.outside),
      ` ${t.loose.inside} of ${t.loose.inside + t.loose.outside} measures inside the real spread`,
      t.strict.inside + t.strict.outside ? ` (same-note ruler: ${t.strict.inside} of ${t.strict.inside + t.strict.outside})` : '',
      ' · ', h('a', { href: `${link('f', doc.id, 'grades')}` }, 'see the measures')));
  }
  if (p.kind === 'refused' && p.refused) {
    out.push(h('div.note.warn', null, h('b', null, 'The grader refused this render. '), String(p.refused.reason || ''),
      p.refused.render_sha ? h('span.small.muted', null, ` (render ${p.refused.render_sha})`) : null));
  }
  if (p.kind === 'listening' && p.no_reference) out.push(h('div.note', null, h('b', null, 'No recording to grade against. '), p.no_reference));
  if (p.kind === 'reference-only') out.push(h('div.note', null, 'The grader holds this recording, but the family has no render of the note yet.'));

  const sources = [];
  if (p.audio.ours) sources.push({ rel: p.audio.ours, label: 'Simulation', side: 'ours' });
  if (p.audio.real) sources.push({ rel: p.audio.real, label: 'Recording', side: 'real' });
  if (sources.length) {
    const player = new Player(sources);
    focusPlayer(player);
    const maps = h('div.maps');
    const figKinds = Object.keys(p.figures);
    let figKind = figKinds.includes('spectrogram') ? 'spectrogram' : figKinds[0];
    const drawMaps = () => {
      const fig = p.figures[figKind];
      maps.replaceChildren();
      if (!fig) return;
      const sides = [['ours', 'Simulation'], ['real', 'Recording'], ['diff', 'Difference: simulation minus recording']];
      for (const [side, label] of sides) {
        const rel = fig.img[side];
        if (!rel) continue;
        const img = h('img', { src: DATA + rel, alt: `${figKind} of the ${label.toLowerCase()}`, loading: 'lazy' });
        maps.append(mapFrame(img, { fLo: fig.f_lo_hz, fHi: fig.f_hi_hz, duration: p.duration_s, player, label }));
      }
      if (fig.img.diff) {
        maps.append(h('div.diffkey', null, h('span', null, `simulation ${fig.diff_range_db} dB quieter`), h('span.grad'), h('span', null, `${fig.diff_range_db} dB louder`)));
      }
      maps.append(h('p.small.muted', null, figKind === 'spectrogram'
        ? `Short-time Fourier spectrogram, 10 ms frames, log frequency ${fig.f_lo_hz}–${fig.f_hi_hz} Hz, ${fig.span_db} dB below the recording's peak. Click a map to play from that point.`
        : `Morlet scalogram (wavelet analysis: sharper in time at high frequency), log frequency ${fig.f_lo_hz}–${fig.f_hi_hz} Hz, ${fig.span_db} dB below the recording's peak.`));
    };
    const figTabs = figKinds.length > 1
      ? h('div.segmented.small', null, figKinds.map(k => {
        const b = h(`button${k === figKind ? '.on' : ''}`, null, k);
        b.addEventListener('click', () => {
          figKind = k;
          for (const x of figTabs.children) x.classList.toggle('on', x === b);
          drawMaps();
        });
        return b;
      }))
      : null;
    drawMaps();
    out.push(h('div.playerbox', null, controls(player),
      sources.length > 1 ? h('p.small.muted', null, 'Both play in step: switching keeps the playhead. Keys: space plays or pauses, 1 and 2 switch. Loudness is matched.') : null,
      sources.length > 1 ? blindTest(doc, p, sources) : null));
    out.push(figTabs, maps);
  } else {
    out.push(empty('No audio for this entry.'));
  }
  out.push(provenance(p));
  return out;
}

function kindLabel(kind) {
  return { graded: 'graded', listening: 'not graded', refused: 'refused', 'reference-only': 'recording only', extra: 'not graded', fixture: 'pre-pivot fixture' }[kind] || kind;
}

function provenance(p) {
  const rows = [];
  if (p.ours) {
    rows.push(['Render', [p.ours.observer, p.ours.git_sha ? ` (commit ${p.ours.git_sha})` : '']]);
    if (p.ours.recipe) rows.push(['Recipe', p.ours.recipe]);
  }
  if (p.real) {
    rows.push(['Recording', [p.real.file || p.real.reference_id || 'unnamed', p.real.segment_s ? ` (segment ${p.real.segment_s.map(x => fmt.num(x, 2)).join('–')} s)` : '']]);
    if (p.real.credit) rows.push(['Credit', [p.real.credit, p.real.licence ? ` · ${p.real.licence}` : '']]);
  }
  if (p.f0_hz) rows.push(['Fundamental', `simulation ${fmt.num(p.f0_hz.ours, 2)} Hz · recording ${fmt.num(p.f0_hz.real, 2)} Hz`]);
  if (p.level_match) rows.push(['Loudness', p.level_match]);
  if (!rows.length) return null;
  return h('dl.prov', null, rows.flatMap(([k, v]) => [h('dt', null, k), h('dd', null, v)]));
}

/** Blind A/B: the two sides in a random order; guess which is the recording. Scores stay in this browser. */
function blindTest(doc, p, sources) {
  const box = h('div.blind');
  const start = h('button.ghost', null, 'Blind test: which is the recording?');
  box.append(start);
  start.addEventListener('click', () => {
    const order = Math.random() < 0.5 ? [0, 1] : [1, 0];
    const blind = new Player(order.map((i, k) => ({ rel: sources[i].rel, label: k ? 'B' : 'A' })));
    focusPlayer(blind);
    const realIs = order.indexOf(sources.findIndex(x => x.side === 'real'));
    const result = h('p.small');
    const guess = k => {
      const right = k === realIs;
      const score = blindScore(doc.id, right);
      result.replaceChildren(h(`b.${right ? 'good' : 'bad'}`, null, right ? 'Right.' : 'Wrong.'),
        ` ${'AB'[realIs]} was the recording. On this browser you have picked the recording ${score.right} of ${score.total} times for ${doc.title}.`);
      for (const b of picks.querySelectorAll('button')) b.disabled = true;
    };
    const picks = h('div.picks', null,
      h('button', { onclick: () => guess(0) }, 'A is the recording'),
      h('button', { onclick: () => guess(1) }, 'B is the recording'));
    box.replaceChildren(h('p.small.muted', null, 'The maps are hidden below while you listen. Switch with 1 and 2.'), controls(blind, { names: ['A', 'B'] }), picks, result,
      h('button.ghost', { onclick: () => { blind.pause(); box.closest('.detail').classList.remove('blinding'); box.replaceChildren(start); } }, 'Done'));
    box.closest('.detail').classList.add('blinding');
  });
  return box;
}

function blindScore(family, right) {
  let all = {};
  try { all = JSON.parse(localStorage.getItem('ie-blind') || '{}'); } catch (e) { all = {}; }
  const s = all[family] || { right: 0, total: 0 };
  s.total += 1;
  if (right) s.right += 1;
  all[family] = s;
  try { localStorage.setItem('ie-blind', JSON.stringify(all)); } catch (e) { /* private mode */ }
  return s;
}

/** A single clip with no grader images: play it and compute its map here. */
function clipDetail(x, kind) {
  const out = [h('div.detailhead', null, h('h2', null, x.title), chip(kindLabel(kind), `k-${kind}`), h('span.small.muted', null, fmt.date(x.date)))];
  if (x.era) out.push(h('div.note', null, x.era));
  if (x.question) out.push(h('p', null, h('b', null, 'Question it was made for: '), x.question));
  if (x.verdict) out.push(h('p', null, h('b', null, 'Verdict recorded then: '), x.verdict));
  if (x.detail) out.push(h('p.small', null, x.detail));
  if (!x.audio) return [...out, empty('No audio.')];
  const player = new Player([{ rel: x.audio, label: x.title }]);
  focusPlayer(player);
  const mapBox = h('div.maps', null, h('p.small.muted', null, 'Computing the spectrogram…'));
  player.ready.then(() => loadBuffer(x.audio)).then(buf => {
    const m = computeMap(buf);
    mapBox.replaceChildren(mapFrame(m.canvas, { fLo: m.fLo, fHi: m.fHi, duration: m.duration, player, label: 'Simulation' }),
      h('p.small.muted', null, `Spectrogram computed in the browser: log frequency ${m.fLo}–${Math.round(m.fHi)} Hz, 80 dB below the clip's own peak.`));
  }).catch(err => mapBox.replaceChildren(empty(`Could not decode: ${err.message}`)));
  out.push(h('div.playerbox', null, controls(player)), mapBox);
  if (x.provenance && Object.keys(x.provenance).length) {
    out.push(h('dl.prov', null, Object.entries(x.provenance).flatMap(([k, v]) => [h('dt', null, k), h('dd', null, String(v))])));
  }
  return out;
}

// ---------------------------------------------------------------- Grades

function ratioClass(m) {
  if (!Number.isFinite(m.gap)) return 'na';
  if (!Number.isFinite(m.ruler_max) || m.ruler_max <= 0) return m.gap === 0 ? 'in' : 'na';
  const q = Math.abs(m.gap) / m.ruler_max;
  return q <= 1 ? 'in' : q <= 2 ? 'near' : 'far';
}

function grades(doc) {
  const graded = doc.pairs.filter(p => p.kind === 'graded');
  if (!graded.length) return empty('No note of this family has been graded against a recording yet.');
  const out = [];
  out.push(h('div.note', null,
    h('p', null, 'Each measure compares the simulation with a real recording of the same note. Its ruler is how far real recordings differ from each other on the same measure: other dynamics, other takes, neighbouring notes. A measure is ', h('b', null, 'inside'), ' when the simulation is no further from the recording than the largest of those real-vs-real distances. The dashed line is the stricter same-note ruler.'),
    doc.baseline ? h('p', null, `How demanding is the ruler? ${doc.baseline.recordings} real recordings from ${doc.baseline.source} (${doc.baseline.compared_with}) scored on the same rulers fall inside ${fmt.pct(doc.baseline.real.inside_fraction)} of the time; the simulation's renders of the same notes, ${fmt.pct(doc.baseline.ours.inside_fraction)}.`) : null));
  const ids = graded[0].metrics.map(m => m.id);
  const labels = Object.fromEntries(graded[0].metrics.map(m => [m.id, m.label]));
  out.push(h('div.heatwrap', null, h('table.heat', null,
    h('thead', null, h('tr', null, h('th', null, 'Note'), ids.map(id => h('th', { title: labels[id] }, h('span', null, labels[id]))))),
    h('tbody', null, graded.map(p => h('tr', null,
      h('th', null, h('a', { href: `#grade-${p.id}` }, [p.note, p.dynamic].filter(Boolean).join(' '))),
      ids.map(id => {
        const m = p.metrics.find(x => x.id === id);
        if (!m) return h('td.na');
        const q = Number.isFinite(m.ruler_max) && m.ruler_max > 0 ? Math.abs(m.gap) / m.ruler_max : null;
        return h(`td.${ratioClass(m)}`, { title: `${m.label}: ours ${fmt.num(m.gap, 3)} ${m.unit}, real spread up to ${fmt.num(m.ruler_max, 3)}` }, q === null ? '·' : q.toFixed(1));
      })))))),
  h('p.small.muted', null, 'Cells are the simulation\'s distance divided by the ruler: 1.0 or less is inside (green), up to 2 is near (amber), beyond 2 is far (red).'));
  for (const p of graded) {
    out.push(h('section.gradecard', { id: `grade-${p.id}` },
      h('div.detailhead', null, h('h3', null, [p.note, p.dynamic].filter(Boolean).join(' ')),
        splitBar(p.tally.loose.inside, p.tally.loose.outside),
        h('span.small', null, `${p.tally.loose.inside} inside, ${p.tally.loose.outside} outside`),
        h('a.small', { href: link('f', doc.id, 'listen', p.id) }, 'listen')),
      h('table.metrics', null,
        h('thead', null, h('tr', null, ['Measure', 'Simulation vs recording', 'Largest real vs real', '', 'Verdict'].map(t => h('th', null, t)))),
        h('tbody', null, p.metrics.map(m => h('tr', null,
          h('td', { title: m.note }, m.label),
          h('td.num', null, `${fmt.num(m.gap, 2)} ${m.unit}`),
          h('td.num', null, Number.isFinite(m.ruler_max) ? fmt.num(m.ruler_max, 2) : 'no ruler'),
          h('td', null, rulerBar(m)),
          h('td', null, m.inside === null ? chip('no ruler', 'na') : m.inside ? chip('inside', 'in') : chip('outside', 'out'),
            m.inside_strict === false && m.inside ? chip('outside same-note', 'near') : null)))))));
  }
  return h('div.grades', null, out);
}

// ---------------------------------------------------------------- History

function history(doc) {
  const notes = doc.pairs.filter(p => p.history.length > 1 || (p.history.length && p.kind === 'graded'));
  if (!notes.length) return empty('No grading history: no note of this family has been graded more than once.');
  const series = notes.map(p => ({
    name: [p.note, p.dynamic].filter(Boolean).join(' '),
    points: p.history.filter(s => s.tally && s.tally.loose.inside + s.tally.loose.outside > 0).map(s => ({
      t: new Date(s.when.replace(/Z?$/, 'Z')),
      y: s.tally.loose.inside / (s.tally.loose.inside + s.tally.loose.outside),
      label: `${fmt.datetime(s.when)}: ${s.tally.loose.inside} of ${s.tally.loose.inside + s.tally.loose.outside} inside (render ${s.git_sha || '?'})`,
    })),
  })).filter(x => x.points.length);
  const out = [h('p.note', null, 'Each time a note is re-rendered and re-graded, the previous grading is kept. The chart shows the share of measures inside the real spread at each grading; refused gradings have no point.')];
  const chart = historyChart(series);
  if (chart) out.push(chart);
  for (const p of notes) {
    const snaps = p.history;
    const moves = movers(snaps);
    out.push(h('section.histcard', null, h('h3', null, [p.note, p.dynamic].filter(Boolean).join(' ')),
      h('table.snaps', null, h('thead', null, h('tr', null, ['Graded', 'Render commit', 'Result'].map(t => h('th', null, t)))),
        h('tbody', null, snaps.slice().reverse().map(s => h(`tr${s.current ? '.current' : ''}`, null,
          h('td', null, fmt.datetime(s.when), s.current ? chip('current') : null),
          h('td.mono', null, s.git_sha || '–'),
          h('td', null, s.tally && s.tally.loose.inside + s.tally.loose.outside
            ? [splitBar(s.tally.loose.inside, s.tally.loose.outside), ` ${s.tally.loose.inside} of ${s.tally.loose.inside + s.tally.loose.outside} inside`]
            : s.refused ? h('span.small', null, `refused: ${s.refused}`) : s.kind))))),
      moves.length ? h('div.movers', null, h('h4', null, 'What moved at the last re-grading'), h('ul', null, moves.map(m => h('li', null, m)))) : null));
  }
  return h('div.history', null, out);
}

/** Measures that crossed the ruler, or moved by more than a fifth of it, between the last two graded snapshots. */
function movers(snaps) {
  const graded = snaps.filter(s => s.tally && Object.keys(s.ratios || {}).length);
  if (graded.length < 2) return [];
  const [a, b] = graded.slice(-2);
  const out = [];
  for (const [id, rb] of Object.entries(b.ratios)) {
    const ra = a.ratios[id];
    if (!Number.isFinite(ra) || !Number.isFinite(rb)) continue;
    const crossed = (ra <= 1) !== (rb <= 1);
    if (crossed || Math.abs(rb - ra) > 0.2) {
      out.push(`${id}: ${fmt.num(ra, 2)} → ${fmt.num(rb, 2)} of the ruler${crossed ? (rb <= 1 ? ' (now inside)' : ' (now outside)') : ''}`);
    }
  }
  return out;
}

// ---------------------------------------------------------------- Qualification

const GATE_DESC = {
  0: 'First audible baseline',
  QP1: 'Passive chain built from data, and the instrument speaks',
  Q2: 'Recognisable as the instrument',
  Q3: 'Inside the real-recording ruler',
  Q4: 'Passes blind listening tests',
  Q5: 'Plays musical phrases',
};

function qualification(doc) {
  const out = [];
  if (doc.milestones.length) {
    out.push(h('h3', null, 'Qualification ladder'),
      h('p.small.muted', null, 'Each gate is a milestone bead. The bar is the state of the work it depends on (closed, in progress, open), which is progress, not qualification: a gate is qualified only when its milestone closes with evidence.'),
      h('table.milestones', null, h('tbody', null, doc.milestones.map(m => {
        const pr = m.prereq;
        return h('tr', null,
          h('td.mono', null, m.code),
          h('td', null, h('b', null, GATE_DESC[m.gate] || m.gate_text), h('div.small.muted', null, m.title)),
          h('td', null, chip(m.status.replace('_', ' '), `st-${m.status}`)),
          h('td.progress', null, h('div.pbar', { title: `${pr.closed} closed, ${pr.in_progress} in progress, ${pr.open} open of ${pr.total}` },
            h('span.c', { style: { flex: String(pr.closed) } }), h('span.p', { style: { flex: String(pr.in_progress) } }), h('span.o', { style: { flex: String(pr.open) } })),
          h('span.small', null, `${pr.closed}/${pr.total} prerequisites closed`)));
      }))));
  }
  for (const e of doc.evidence) {
    let chart = null;
    if (e.kind === 'body-modes') {
      chart = h('div', null, bodyModes(e), h('div.legend', null,
        h('span', null, h('i.shape.circle'), 'A0 (air)'), h('span', null, h('i.shape.square'), 'CBR'), h('span', null, h('i.shape.diamond'), 'B1-'), h('span', null, h('i.shape.triangle'), 'B1+'),
        h('span', null, h('i.hitkey'), 'inside its band'), h('span', null, h('i.misskey'), 'outside')));
    } else if (e.kind === 'bowed-sweep') chart = bowedSweep(e);
    else if (e.kind === 'band-decays') chart = h('div', null, bandDecays(e), h('ul.small', null, e.pieces.map(p => h('li', null, h('b', null, `${p.piece}: `), p.what))));
    else if (e.kind === 'before-after') chart = beforeAfter(e);
    out.push(h('section.evidence', null, h('h3', null, e.title), h('p.small', null, e.caption), chart, h('p.small.muted.mono', null, e.source)));
  }
  if (doc.claims.length) {
    out.push(h('h3', null, 'Physics claims register'),
      h('p.small.muted', null, 'Rows of instrument-claims.json for the physics this family uses. Green means a gate test exists and passes; ungated means the claim has no gate yet.'),
      h('table.claims', null, h('thead', null, h('tr', null, ['Physics', 'Image', 'Quantity', 'Gate', 'Live default', 'Evidence'].map(t => h('th', null, t)))),
        h('tbody', null, doc.claims.map(c => h('tr', null, h('td', null, c.filling), h('td', null, c.image), h('td', null, c.qoi),
          h('td', null, chip(c.gate, c.gate === 'green' ? 'in' : c.gate === 'refused' ? 'out' : 'na')), h('td', null, c.live_default),
          h('td.small.mono', null, (c.evidence || []).map(x => h('div', null, x))))))));
  }
  const bs = doc.beads;
  out.push(h('p.small.muted', null, `Beads naming this family: ${Object.entries(bs).map(([k, v]) => `${v} ${k.replace('_', ' ')}`).join(', ') || 'none'}.`));
  return out.length ? h('div.qual', null, out) : empty('No qualification evidence for this family yet.');
}

function beforeAfter(e) {
  const [a, b] = e.rows;
  if (!a || !b) return null;
  const keys = Object.keys(a.distance).filter(k => k in b.distance);
  return h('div', null, h('table.metrics', null,
    h('thead', null, h('tr', null, ['', a.label, b.label, ''].map(t => h('th', null, t)))),
    h('tbody', null,
      h('tr', null, h('td', null, 'fundamental (recording ' + fmt.num(a.real_f0, 2) + ' Hz)'), h('td.num', null, `${fmt.num(a.ours_f0, 2)} Hz`), h('td.num', null, `${fmt.num(b.ours_f0, 2)} Hz`),
        h('td', null, Math.abs(b.ours_f0 - b.real_f0) < Math.abs(a.ours_f0 - a.real_f0) ? chip('closer', 'in') : chip('further', 'out'))),
      keys.map(k => h('tr', null, h('td', null, k.replaceAll('_', ' ')), h('td.num', null, fmt.num(a.distance[k], 2)), h('td.num', null, fmt.num(b.distance[k], 2)),
        h('td', null, Math.abs(b.distance[k]) < Math.abs(a.distance[k]) ? chip('closer', 'in') : chip('further', 'out')))))));
}

// ---------------------------------------------------------------- Timeline

function timeline(doc) {
  const events = [
    ...doc.timeline.map(c => ({ t: c.when, kind: 'commit', c })),
    ...doc.closes.map(b => ({ t: (b.closed_at || '').slice(0, 16).replace('T', ' '), kind: 'close', b })),
  ].sort((x, y) => (x.t < y.t ? 1 : -1));
  if (!events.length) return empty('No commits or bead closes name this family yet.');
  const list = h('ol.timeline');
  let filter = 'all';
  let q = '';
  const draw = () => {
    const rx = q ? new RegExp(q.replace(/[.*+?^${}()|[\]\\]/g, '\\$&'), 'i') : null;
    const shown = events.filter(e => (filter === 'all' || e.kind === filter) && (!rx || rx.test(e.kind === 'commit' ? e.c.subject : `${e.b.title} ${e.b.reason}`)));
    list.replaceChildren(...shown.slice(0, 300).map(e => e.kind === 'commit'
      ? h('li.commit', null, h('span.when', null, fmt.datetime(e.t)), chip(e.c.type, `t-${e.c.type}`), h('span.subj', null, e.c.subject), h('span.mono.small.muted', null, ` ${e.c.sha}`))
      : h('li.close', null, h('span.when', null, fmt.datetime(e.t)), chip('closed', 'in'), h('span.subj', null, e.b.title), e.b.reason ? h('div.small.muted', null, e.b.reason) : null)));
    if (!shown.length) list.append(empty('Nothing matches.'));
  };
  const seg = h('div.segmented.small', null, [['all', 'All'], ['commit', 'Commits'], ['close', 'Bead closes']].map(([k, label]) => {
    const b = h(`button${k === filter ? '.on' : ''}`, null, label);
    b.addEventListener('click', () => {
      filter = k;
      for (const x of seg.children) x.classList.toggle('on', x === b);
      draw();
    });
    return b;
  }));
  const search = h('input.search', { type: 'search', placeholder: 'Filter…' });
  search.addEventListener('input', () => { q = search.value.trim(); draw(); });
  draw();
  return h('div', null, h('div.toolbar', null, seg, search,
    h('span.small.muted', null, `${doc.timeline_total} commits since ${fmt.date('2026-10-15')} mention this family${doc.timeline_total > doc.timeline.length ? ` (newest ${doc.timeline.length} shown)` : ''}.`)), list);
}

// ---------------------------------------------------------------- References

function references(doc) {
  const refs = doc.references;
  const cov = doc.reference_coverage || {};
  if (!refs.length) return empty('No screened reference recording for this family yet.');
  const panel = h('div.refplayer', null, h('p.small.muted', null, 'Pick a recording below to hear it and see its spectrogram.'));
  const credits = [...new Map(refs.map(r => [r.licence_id, r])).values()];
  const tbody = h('tbody');
  let q = '';
  const open = r => {
    const player = new Player([{ rel: r.audio, label: r.id }]);
    focusPlayer(player);
    const mapBox = h('div.maps', null, h('p.small.muted', null, 'Decoding…'));
    panel.replaceChildren(h('div.detailhead', null, h('h3', null, `${r.note || ''} ${r.dynamic || ''}`), h('span.small.muted', null, `${r.instrument || ''} · ${r.kind || ''}${r.piece ? ` · ${r.piece}` : ''}`)),
      controls(player), mapBox,
      h('p.small', null, r.author ? `${r.author} · ` : '', r.url ? h('a', { href: r.url, target: '_blank', rel: 'noopener' }, 'source') : '', ` · ${r.licence}`, r.credit ? ` · ${r.credit}` : ''));
    player.ready.then(() => loadBuffer(r.audio)).then(buf => {
      const m = computeMap(buf);
      mapBox.replaceChildren(mapFrame(m.canvas, { fLo: m.fLo, fHi: m.fHi, duration: m.duration, player, label: 'Recording' }));
      player.play();
    }).catch(err => mapBox.replaceChildren(empty(`Could not decode: ${err.message}`)));
  };
  const draw = () => {
    const rx = q ? new RegExp(q.replace(/[.*+?^${}()|[\]\\]/g, '\\$&'), 'i') : null;
    tbody.replaceChildren(...refs.filter(r => !rx || rx.test([r.note, r.dynamic, r.kind, r.piece, r.instrument, r.author, r.licence].join(' '))).map(r => {
      const tr = h('tr.clickable', null, h('td', null, r.note || '–'), h('td', null, r.dynamic || '–'), h('td', null, r.kind || '–'),
        h('td', null, r.piece || ''), h('td.num', null, r.duration_s ? `${fmt.num(r.duration_s, 1)} s` : ''), h('td.small', null, r.instrument || ''), h('td.small', null, r.licence));
      tr.addEventListener('click', () => {
        for (const x of tbody.children) x.classList.remove('on');
        tr.classList.add('on');
        open(r);
      });
      return tr;
    }));
  };
  const search = h('input.search', { type: 'search', placeholder: 'Filter by note, dynamic, source…' });
  search.addEventListener('input', () => { q = search.value.trim(); draw(); });
  draw();
  return h('div.refs', null,
    h('p.note', null, `${refs.length} recordings passed the screening (${cov.valid ?? refs.length} of ${cov.rows ?? '?'} candidate clips; the others failed checks such as signal-to-noise, a single clean note, decay or pitch). They stay on this machine: several are not licensed for redistribution.`),
    panel,
    h('div.toolbar', null, search),
    h('div.tablewrap', null, h('table.reftable', null, h('thead', null, h('tr', null, ['Note', 'Dynamic', 'Kind', 'Piece', 'Length', 'Source', 'Licence'].map(t => h('th', null, t)))), tbody)),
    h('div.credits', null, h('h4', null, 'Credits'), h('ul.small', null, credits.map(r => h('li', null, `${r.licence}: `, r.credit || r.instrument || r.source, r.licence_notes ? h('span.muted', null, ` (${r.licence_notes})`) : null)))));
}
