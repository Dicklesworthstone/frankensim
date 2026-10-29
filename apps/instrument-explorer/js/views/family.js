// One family: listen (simulation against recording), grades, grading history, qualification evidence, timeline and
// the screened reference recordings.
import { h, getJSON, fmt, badge, splitBar, link, chip, empty, DATA, scrollMotion, stack } from '../app.js';
import { Player, controls, focusPlayer, loadBuffer } from '../audio.js';
import { mapFrame, computeMap } from '../spectro.js';
import { rulerBar, historyChart, bodyModes, bowedSweep, bandDecays } from '../charts.js';
import { phraseCard } from './phrase.js';

let lastTitle = null;
export function title() {
  return lastTitle;
}

const TABS = [
  ['hear', 'Hear it'],
  ['listen', 'Compare'],
  ['grades', 'Grades'],
  ['history', 'History'],
  ['qualification', 'Qualification'],
  ['timeline', 'Timeline'],
  ['references', 'References'],
];

const VIEWS = { hear, listen, grades, history, qualification, timeline, references };

/** One page per family: every section stacked, a section bar that follows the scroll, and each section built as it
 *  comes near the screen. A section link (#/f/<family>/<section>[/<item>]) scrolls to it and, with an item, rebuilds
 *  only that section around the item; the rest of the page and anything playing stay put. */
export async function render(root, r, sameFamily) {
  const doc = await getJSON(`family/${r.id}.json`);
  doc.phrases = doc.phrases || [];
  lastTitle = doc.title;
  let page = root.querySelector(`.familypage[data-id="${CSS.escape(r.id)}"]`);
  const counts = {
    hear: doc.phrases.length,
    listen: doc.pairs.filter(p => p.kind !== 'reference-only').length + doc.extras.length + doc.fixtures.length,
    grades: doc.pairs.filter(p => p.kind === 'graded').length,
    history: doc.pairs.filter(p => p.history.length > 1).length,
    qualification: doc.milestones.length + doc.evidence.length,
    timeline: doc.timeline_total + doc.closes.length,
    references: doc.references.length,
  };
  const shown = TABS.filter(([id]) => id !== 'hear' || counts.hear);
  if (!sameFamily || !page) {
    const nav = h('nav.tabs.sectionbar', { 'aria-label': `${doc.title} sections` }, shown.map(([id, label]) =>
      h('a', { href: link('f', doc.id, id), 'data-sec': id }, label, h('span.count', null, String(counts[id])))));
    const sections = shown.map(([id, label]) => {
      const body = h('div.fbody', { 'data-built': '' });
      return h(`section.fsec`, { id: `sec-${id}`, 'data-sec': id, 'aria-labelledby': `h-${id}` },
        h('h2.fsechead', { id: `h-${id}` }, label, h('span.count', null, String(counts[id]))), body);
    });
    if (nav.firstChild) nav.firstChild.classList.add('active'); // until the scroll reaches a section
    page = h('div.familypage', { 'data-id': r.id }, header(doc), nav, h('div.fsections', null, sections));
    root.replaceChildren(page);
    // A section rebuilt around a newly picked item keeps its list's scroll position and shows the pick inside it.
    const build = (sec, item) => {
      const body = sec.querySelector('.fbody');
      const kept = [...body.querySelectorAll('.itemlist, .plist')].map(l => l.scrollTop);
      body.dataset.built = item || '-';
      body.replaceChildren(VIEWS[sec.dataset.sec](doc, item));
      body.querySelectorAll('.itemlist, .plist').forEach((l, k) => {
        if (kept[k]) l.scrollTop = kept[k];
        const on = l.querySelector('.on');
        if (!on || !l.isConnected) return;
        const d = on.getBoundingClientRect().top - l.getBoundingClientRect().top;
        if (d < 0 || d + on.offsetHeight > l.clientHeight) l.scrollTop += d - l.clientHeight / 3;
      });
    };
    // Every section is built now: all of them together take a few milliseconds, and a page whose height never changes
    // under the reader keeps tab jumps, the scrollbar and back/forward positions exact.
    for (const sec of sections) build(sec, sec.dataset.sec === r.tab ? r.item : null);
    // The bar marks the last section whose top has passed under it. It slides sideways to keep that tab visible by
    // setting its own scrollLeft: scrollIntoView would also scroll the page, which cancels a smooth jump halfway.
    let current = null;
    let wanted = null; // a short last section can never reach the bar; at the page's end the tab that was asked for wins
    const pick = () => {
      if (!page.isConnected) return removeEventListener('scroll', onScroll);
      const edge = nav.getBoundingClientRect().bottom + 24;
      let cur = sections[0];
      for (const s of sections) if (s.getBoundingClientRect().top <= edge) cur = s;
      if (wanted && innerHeight + scrollY >= document.documentElement.scrollHeight - 2) cur = wanted;
      if (cur === current) return;
      current = cur;
      for (const a of nav.querySelectorAll('a')) {
        const on = a.dataset.sec === cur.dataset.sec;
        a.classList.toggle('active', on);
        if (on && (a.offsetLeft < nav.scrollLeft || a.offsetLeft + a.offsetWidth > nav.scrollLeft + nav.clientWidth)) {
          nav.scrollLeft = a.offsetLeft - (nav.clientWidth - a.offsetWidth) / 2;
        }
      }
    };
    let queued = false;
    const onScroll = () => {
      if (queued) return;
      queued = true;
      requestAnimationFrame(() => { queued = false; pick(); });
    };
    addEventListener('scroll', onScroll, { passive: true });
    page._want = sec => { wanted = sec; };
    page._build = build;
  }
  if (r.tab) {
    const sec = page.querySelector(`.fsec[data-sec="${CSS.escape(r.tab)}"]`);
    if (sec) {
      const body = sec.querySelector('.fbody');
      // A section tab names no item: it goes to the section and keeps what was picked there (and anything playing).
      // An item link rebuilds the section around its item, and so does back/forward to the section's bare address.
      const keep = !r.item && r.restoreY === undefined;
      if (!keep && body.dataset.built !== (r.item || '-')) page._build(sec, r.item);
      page._want(sec);
      // A section link goes to the section; an item picked inside it goes to its detail, and only if that is off
      // screen. Back and forward return to the remembered position instead (r.restoreY, set by the router).
      const target = (r.item && sec.querySelector('.detail')) || sec;
      if (r.restoreY === undefined) {
        requestAnimationFrame(() => {
          const top = target.getBoundingClientRect().top;
          if (!r.item || top < 0 || top > innerHeight * 0.6) aim(target);
        });
      }
    }
  }
}

/** Scrolls a section (or an item's detail) to just under the bars once the web fonts have settled the layout, and
 *  corrects the landing where a smooth scroll ends if something above still moved (a small correction only: a larger
 *  difference means the reader scrolled away during the jump). */
function aim(target) {
  const go = behavior => target.scrollIntoView({ behavior, block: 'start' });
  const off = () => Math.abs(target.getBoundingClientRect().top
    - (parseFloat(getComputedStyle(document.documentElement).scrollPaddingTop) || 0)
    - (parseFloat(getComputedStyle(target).scrollMarginTop) || 0));
  (document.fonts ? document.fonts.ready : Promise.resolve()).then(() => {
    if (off() <= 4) return; // already there: no scroll, so no scrollend either
    const smooth = scrollMotion(off()) === 'smooth';
    go(smooth ? 'smooth' : 'auto');
    if (!smooth || !('onscrollend' in window)) return;
    addEventListener('scrollend', () => {
      const d = off();
      if (d > 4 && d < 300) go('auto');
    }, { once: true });
  });
}

// ---------------------------------------------------------------- Phrases

function hear(doc) {
  if (!doc.phrases.length) return empty('No phrases have been rendered for this family yet. Single notes against recordings are under Compare.');
  return h('div.pgrid.infamily', null, doc.phrases.map(p => phraseCard(p, doc.title)));
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
  const summary = h(`p.summary${doc.summary_derived ? '.muted' : ''}${(doc.summary || '').length > 320 ? '.clamped' : ''}`, null, doc.summary);
  const more = summary.classList.contains('clamped')
    ? h('button.more', { type: 'button', 'aria-expanded': 'false', onclick: e => {
      const open = summary.classList.toggle('clamped');
      e.target.textContent = open ? 'Read more' : 'Show less';
      e.target.setAttribute('aria-expanded', String(!open));
    } }, 'Read more')
    : null;
  return h('header.familyhead', null,
    h('nav.crumbs', { 'aria-label': 'Breadcrumb' }, h('a', { href: '#/instruments' }, 'Instruments'), h('span', { 'aria-hidden': 'true' }, '/'), h('span', null, doc.section)),
    h('div.titlerow', null, h('h1', null, doc.title), badge(v.badge, v.label)),
    st.reference_instrument ? h('p.small.muted', null, `Graded against: ${st.reference_instrument}`) : null,
    summary, more,
    st.draft ? h('p.small.muted', { role: 'note' }, st.draft) : null,
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
      : p.kind === 'refused' ? 'refused' : p.kind === 'listening' ? 'no recording' : p.kind === 'compared' ? 'with recording' : 'recording only';
    const group = { graded: 'Graded against a recording', listening: 'Rendered, no recording to grade against', refused: 'Refused by the grader', compared: 'Beside a recording (not graded)', 'reference-only': 'Recordings the grader holds, not rendered yet' }[p.kind];
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
  if (p.kind === 'compared' && p.compare) out.push(compareTable(p));

  const sources = [];
  if (p.audio.ours) sources.push({ rel: p.audio.ours, label: 'Simulation', side: 'sim' });
  if (p.audio.real) sources.push({ rel: p.audio.real, label: 'Recording', side: 'real' });
  if (p.real && p.real.withheld) {
    out.push(h('div.note', null, h('b', null, 'The recording is not playable here. '), p.real.withheld));
  }
  if (sources.length) {
    const player = new Player(sources, { title: `${doc.title} ${[p.note, p.dynamic].filter(Boolean).join(' ')}`, sub: sources.length > 1 ? 'simulation and recording' : 'simulation', duration: p.duration_s });
    player.load().catch(() => {});
    focusPlayer(player);
    const maps = h('div.maps');
    const figKinds = Object.keys(p.figures);
    let figKind = figKinds.includes('spectrogram') ? 'spectrogram' : figKinds[0];
    const drawMaps = () => {
      const fig = p.figures[figKind];
      maps.replaceChildren();
      if (!fig) return;
      const sides = [['ours', 'Simulation', 'sim'], ['real', 'Recording', 'real'], ['diff', 'Difference: simulation minus recording', null]];
      for (const [side, label, cls] of sides) {
        const rel = fig.img[side];
        if (!rel) continue;
        const img = h('img', { src: DATA + rel, alt: `${figKind} of the ${label.toLowerCase()}`, loading: 'lazy' });
        maps.append(mapFrame(img, { fLo: fig.f_lo_hz, fHi: fig.f_hi_hz, duration: p.duration_s, player, label, side: cls }));
      }
      if (fig.img.diff) {
        maps.append(h('div.diffkey', null, h('span', null, `simulation ${fig.diff_range_db} dB quieter`),
          h('span.grad', fig.diff_key ? { style: { background: fig.diff_key } } : null), h('span', null, `${fig.diff_range_db} dB louder`)));
      }
      maps.append(h('p.small.muted', null, figKind === 'spectrogram'
        ? `Short-time Fourier spectrogram, 10 ms frames, log frequency ${fig.f_lo_hz}–${fig.f_hi_hz} Hz, ${fig.span_db} dB below the recording's peak. Click a map to play from that point.`
        : `Morlet scalogram (wavelet analysis: sharper in time at high frequency), log frequency ${fig.f_lo_hz}–${fig.f_hi_hz} Hz, ${fig.span_db} dB below the recording's peak.`));
    };
    const figTabs = figKinds.length > 1
      ? h('div.segmented.small', null, figKinds.map(k => {
        const b = h(`button${k === figKind ? '.on' : ''}`, { type: 'button' }, k);
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
      sources.length > 1 ? h('p.small.muted', null, 'Both play in step, so switching keeps your place. Their loudness is matched.') : null,
      sources.length > 1 ? blindTest(doc, p, sources) : null));
    out.push(figTabs, maps);
  } else {
    out.push(empty('No audio for this entry.'));
  }
  out.push(provenance(p));
  return out;
}

function kindLabel(kind) {
  return { graded: 'graded', listening: 'not graded', refused: 'refused', compared: 'not graded', 'reference-only': 'recording only', extra: 'not graded', fixture: 'pre-pivot fixture' }[kind] || kind;
}

/** A kit piece beside its recording: the kit grader's descriptors, and its octave-balance gap against real-vs-real distances. */
function compareTable(p) {
  const c = p.compare;
  const r = c.ruler || {};
  const label = { centroid_hz: ['Spectral centroid', 'Hz'], above_1k_pct: ['Energy above 1 kHz', '%'], above_4k_pct: ['Energy above 4 kHz', '%'] };
  // Two significant figures, never exponent notation: these span 0.001 % to 11 kHz, and below 0.001 a band is all but
  // empty (the kit11 kick's energy above 4 kHz is 2.5e-8 %, which toString() writes in exponent form)
  const sig = v => (v === null || v === undefined ? '–' : v >= 100 ? v.toFixed(0) : v === 0 ? '0' : v < 0.001 ? '< 0.001'
    : Number(v.toPrecision(2)).toString());
  const rows = Object.entries(c.descriptors).map(([k, v]) => {
    const [name, unit] = label[k] || [k.replaceAll('_', ' '), ''];
    return h('tr', null, h('td', null, name), h('td.num', null, `${sig(v.ours)} ${unit}`), h('td.num', null, `${sig(v.real)} ${unit}`));
  });
  return h('div', null,
    h('p.small', null, `Octave-band balance over the first 250 ms: the render is ${fmt.num(c.balance_gap_db, 1)} dB (RMS over bands) from the recording. `,
      `Two strokes of the real drum at the same dynamic differ by ${fmt.num(r.real_same_stroke_db, 1)} dB, and one dynamic apart by ${fmt.num(r.real_other_dynamic_db, 1)} dB.`),
    h('div.tablewrap', null, h('table.metrics', null,
      h('thead', null, h('tr', null, ['First 250 ms after the stroke', 'Simulation', 'Recording'].map(t => h('th', null, t)))),
      h('tbody', null, rows))),
    c.caveat ? h('div.note.warn', null, h('p', null, c.caveat)) : null);
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
  if (p.f0_hz) {
    rows.push(['Fundamental', [Number.isFinite(p.f0_hz.ours) ? `simulation ${fmt.num(p.f0_hz.ours, 2)} Hz` : null,
      Number.isFinite(p.f0_hz.real) ? `recording ${fmt.num(p.f0_hz.real, 2)} Hz` : null].filter(Boolean).join(' · ') || '–']);
  }
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
  const player = new Player([{ rel: x.audio, label: x.title, side: 'sim' }], { title: x.title, sub: 'simulation, not graded' });
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
  if (m.inside === true) return 'in';  // the build's verdict on unrounded values; the rounded ratio can sit on 1.0
  if (!Number.isFinite(m.ruler_max) || m.ruler_max <= 0) return m.gap === 0 ? 'in' : 'na';
  const q = Math.abs(m.gap) / m.ruler_max;
  return q <= 1 && m.inside !== false ? 'in' : q <= 2 ? 'near' : 'far';
}

function grades(doc) {
  const graded = doc.pairs.filter(p => p.kind === 'graded');
  if (!graded.length) return empty('No note of this family has been graded against a recording yet.');
  const out = [];
  out.push(h('div.note', null,
    h('p', null, 'Each measure compares the simulation with a real recording of the same note. Its ruler is how far real recordings differ from each other on the same measure: other dynamics, other takes, neighbouring notes. A measure is ', h('b', null, 'inside'), ' when the simulation is no further from the recording than the largest of those real-vs-real distances. The dashed line is the stricter same-note ruler.'),
    doc.baseline ? h('p', null, `How demanding is the ruler? ${doc.baseline.recordings} real recordings from ${doc.baseline.source}, measured against ${doc.baseline.compared_with}, fall inside ${fmt.pct(doc.baseline.real.inside_fraction)} of the time; the simulation's renders of the same notes, ${fmt.pct(doc.baseline.ours.inside_fraction)}.`) : null));
  const ids = graded[0].metrics.map(m => m.id);
  const labels = Object.fromEntries(graded[0].metrics.map(m => [m.id, m.label]));
  out.push(h('div.heatwrap', null, h('table.heat', null,
    h('thead', null, h('tr', null, h('th', null, 'Note'), ids.map(id => h('th', { title: labels[id] }, h('span', null, labels[id]))))),
    h('tbody', null, graded.map(p => h('tr', null,
      // a fragment link would be read as a route (the home page), so the note scrolls to its card itself
      h('th', null, h('a', { href: link('f', doc.id, 'grades'), onclick: e => {
        e.preventDefault();
        const card = document.getElementById(`grade-${p.id}`);
        if (card) {
          if (card.tagName === 'DETAILS') card.open = true;
          aim(card);
        }
      } }, [p.note, p.dynamic].filter(Boolean).join(' '))),
      ids.map(id => {
        const m = p.metrics.find(x => x.id === id);
        if (!m) return h('td.na');
        const q = Number.isFinite(m.ruler_max) && m.ruler_max > 0 ? Math.abs(m.gap) / m.ruler_max : null;
        return h(`td.${ratioClass(m)}`, { title: `${m.label}: ours ${fmt.num(m.gap, 3)} ${m.unit}, real spread up to ${fmt.num(m.ruler_max, 3)}` }, q === null ? '·' : q.toFixed(1));
      })))))),
  h('p.small.muted', null, 'Cells are the simulation\'s distance divided by the ruler: 1.0 or less is inside (green), up to 2 is near (amber), beyond 2 is far (red).'));
  // With more than a few notes the heat map above is the overview and each note's card opens on demand: 21 open cards
  // made the clarinet's Grades section 36,000 px tall on a phone.
  const fold = graded.length > 3;
  for (const p of graded) {
    const head = h('div.detailhead', null, h('h3', null, [p.note, p.dynamic].filter(Boolean).join(' ')),
      splitBar(p.tally.loose.inside, p.tally.loose.outside),
      h('span.small', null, `${p.tally.loose.inside} inside, ${p.tally.loose.outside} outside`),
      h('a.small', { href: link('f', doc.id, 'listen', p.id) }, 'listen'));
    out.push(h(fold ? 'details.gradecard' : 'section.gradecard', { id: `grade-${p.id}` },
      fold ? h('summary', null, head) : head,
      stack(h('table.metrics.grade', null,
        h('thead', null, h('tr', null, ['Measure', 'Simulation vs recording', 'Largest real vs real', '', 'Verdict'].map(t => h('th', null, t)))),
        h('tbody', null, p.metrics.map(m => h('tr', null,
          h('td', { title: m.note }, m.label),
          h('td.num', null, `${fmt.num(m.gap, 2)} ${m.unit}`),
          h('td.num', null, Number.isFinite(m.ruler_max) ? fmt.num(m.ruler_max, 2) : 'no ruler'),
          h('td', null, rulerBar(m)),
          h('td', null, m.inside === null ? chip('no ruler', 'na') : m.inside ? chip('inside', 'in') : chip('outside', 'out'),
            m.inside_strict === false && m.inside ? chip('outside same-note', 'near') : null))))),
      ['', 'simulation', 'real vs real up to', '', ''])));
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
      tally: s.tally.loose,
    })),
  })).filter(x => x.points.length);
  // Too many notes to tell apart by colour: each is a thin grey line (its points keep their titles) under one line for
  // the family, the share of all its measures inside, counting each note's latest grading at that moment. A grading run
  // grades its notes seconds apart, so gradings less than 30 minutes apart count as one moment (the run's last).
  const many = series.length > 6;
  if (many) {
    const runs = series.map(x => x.points.map(p => ({ t: +p.t, tally: p.tally })).sort((a, b) => a.t - b.t));
    const times = [];
    for (const t of [...new Set(runs.flat().map(g => g.t))].sort((a, b) => a - b)) {
      if (times.length && t - times[times.length - 1] < 18e5) times[times.length - 1] = t;
      else times.push(t);
    }
    const points = times.map(t => {
      let inside = 0;
      let n = 0;
      let count = 0;
      for (const r of runs) {
        const last = r.filter(g => g.t <= t).pop();
        if (last) [inside, n, count] = [inside + last.tally.inside, n + last.tally.inside + last.tally.outside, count + 1];
      }
      return { t: new Date(t), y: inside / n, label: `${fmt.datetime(new Date(t).toISOString())}: ${inside} of ${n} measures inside, over ${fmt.plural(count, 'note')}` };
    });
    for (const x of series) x.muted = true;
    series.push({ name: 'All notes', color: 'var(--sim)', points });
  }
  const out = [h('p.note', null, `Each time a note is re-rendered and re-graded, the previous grading is kept. The chart shows the share of measures inside the real spread at each grading; refused gradings have no point.${many ? ' The grey lines are single notes; the bold line is all of them together, counting each note at its latest grading.' : ''}`)];
  const chart = historyChart(series);
  if (chart) out.push(chart);
  // the chart is the overview; with more than a few notes each note's table opens on demand, as on Grades
  const fold = notes.length > 3;
  for (const p of notes) {
    const snaps = p.history;
    const moves = movers(snaps);
    const name = h('h3', null, [p.note, p.dynamic].filter(Boolean).join(' '));
    const last = snaps.filter(s => s.tally && s.tally.loose.inside + s.tally.loose.outside).pop();
    const head = fold ? h('summary', null, h('div.detailhead', null, name,
      h('span.small.muted', null, `${fmt.plural(snaps.length, 'grading')}${last ? `; latest ${last.tally.loose.inside} of ${last.tally.loose.inside + last.tally.loose.outside} inside` : ''}`))) : name;
    out.push(h(fold ? 'details.histcard' : 'section.histcard', null, head,
      stack(h('table.snaps', null, h('thead', null, h('tr', null, ['Graded', 'Render commit', 'Result'].map(t => h('th', null, t)))),
        h('tbody', null, snaps.slice().reverse().map(s => h(`tr${s.current ? '.current' : ''}`, null,
          h('td', null, fmt.datetime(s.when), s.current ? chip('current') : null),
          h('td.mono', null, s.git_sha || '–'),
          h('td', null, s.tally && s.tally.loose.inside + s.tally.loose.outside
            ? [splitBar(s.tally.loose.inside, s.tally.loose.outside), ` ${s.tally.loose.inside} of ${s.tally.loose.inside + s.tally.loose.outside} inside`]
            : s.refused ? h('span.small', null, `refused: ${s.refused}`) : s.kind))))), ['', 'render', '']),
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
      h('div.tablewrap', null, stack(h('table.milestones', null, h('tbody', null, doc.milestones.map(m => {
        const pr = m.prereq;
        return h('tr', null,
          h('td.mono', null, m.code),
          h('td', null, h('b', null, GATE_DESC[m.gate] || m.gate_text), h('div.small.muted', null, m.title)),
          h('td', null, chip(m.status.replace('_', ' '), `st-${m.status}`)),
          h('td.progress', null, h('div.pbar', { title: `${pr.closed} closed, ${pr.in_progress} in progress, ${pr.open} open of ${pr.total}` },
            h('span.c', { style: { flex: String(pr.closed) } }), h('span.p', { style: { flex: String(pr.in_progress) } }), h('span.o', { style: { flex: String(pr.open) } })),
          h('span.small', null, `${pr.closed}/${pr.total} prerequisites closed`)));
      }))))));
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
      h('div.tablewrap', null, stack(h('table.claims', null, h('thead', null, h('tr', null, ['Physics', 'Image', 'Quantity', 'Gate', 'Live default', 'Evidence'].map(t => h('th', null, t)))),
        h('tbody', null, doc.claims.map(c => h('tr', null, h('td', null, c.filling), h('td', null, c.image), h('td', null, c.qoi),
          h('td', null, chip(c.gate, c.gate === 'green' ? 'in' : c.gate === 'refused' ? 'out' : 'na')), h('td', null, c.live_default),
          h('td.small.mono', null, evidenceList(c.evidence || [])))))))));
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

// A claim's evidence names tests as file::test. They are grouped by file, and a list of three or more folds behind its
// count: one trumpet claim listed six full paths, which ran to 30 lines of monospace on a phone.
function evidenceList(list) {
  const files = new Map();
  const other = [];
  for (const x of list) {
    const m = /^([^\s:]+)::(.+)$/.exec(x);
    if (!m) { other.push(x); continue; }
    if (!files.has(m[1])) files.set(m[1], []);
    files.get(m[1]).push(m[2]);
  }
  const rows = [...other.map(x => h('div', null, x)),
    ...[...files].map(([f, tests]) => h('div.evfile', null, h('span.muted', null, `${f} `), tests.map(t => h('span.evtest', null, t))))];
  if (list.length < 3) return rows;
  const tests = list.length - other.length;
  const label = [tests && `${tests} test${tests > 1 ? 's' : ''}`, other.length && `${other.length} record${other.length > 1 ? 's' : ''}`].filter(Boolean).join(' and ');
  return h('details.evlist', null, h('summary', null, label), rows);
}

// ---------------------------------------------------------------- Timeline

function timeline(doc) {
  const events = [
    ...doc.timeline.map(c => ({ t: c.when, kind: 'commit', c })),
    ...doc.closes.map(b => ({ t: (b.closed_at || '').slice(0, 16).replace('T', ' '), kind: 'close', b })),
  ].sort((x, y) => (x.t < y.t ? 1 : -1));
  if (!events.length) return empty('No commits or bead closes name this family yet.');
  const list = h('ol.timeline');
  // The newest 25 first, more on request: 300 wrapped commit subjects made the piano's timeline 11,000 px tall on a phone.
  // The filter and the search still look through every event.
  const STEP = 25;
  let filter = 'all';
  let q = '';
  let limit = STEP;
  const draw = () => {
    const rx = q ? new RegExp(q.replace(/[.*+?^${}()|[\]\\]/g, '\\$&'), 'i') : null;
    const shown = events.filter(e => (filter === 'all' || e.kind === filter) && (!rx || rx.test(e.kind === 'commit' ? e.c.subject : `${e.b.title} ${e.b.reason}`)));
    list.replaceChildren(...shown.slice(0, limit).map(e => e.kind === 'commit'
      ? h('li.commit', null, h('span.when', null, fmt.datetime(e.t)), chip(e.c.type, `t-${e.c.type}`), h('span.subj', null, e.c.subject), h('span.mono.small.muted', null, ` ${e.c.sha}`))
      : h('li.close', null, h('span.when', null, fmt.datetime(e.t)), chip('closed', 'in'), h('span.subj', null, e.b.title), e.b.reason ? h('div.small.muted', null, e.b.reason) : null)));
    if (!shown.length) list.append(empty('Nothing matches.'));
    if (shown.length > limit) {
      const more = h('button.btn.small', { type: 'button' }, `Show ${Math.min(STEP * 2, shown.length - limit)} more (${shown.length - limit} not shown)`);
      more.addEventListener('click', () => {
        const first = limit;
        limit += STEP * 2;
        draw();
        // Keyboard focus lands on the first newly shown event instead of falling back to the page.
        const subj = list.children[first]?.querySelector('.subj');
        if (subj) { subj.setAttribute('tabindex', '-1'); subj.focus({ preventScroll: true }); }
      });
      list.append(h('li.more', null, more));
    }
  };
  const seg = h('div.segmented.small', null, [['all', 'All'], ['commit', 'Commits'], ['close', 'Bead closes']].map(([k, label]) => {
    const b = h(`button${k === filter ? '.on' : ''}`, null, label);
    b.addEventListener('click', () => {
      filter = k;
      limit = STEP;
      for (const x of seg.children) x.classList.toggle('on', x === b);
      draw();
    });
    return b;
  }));
  const search = h('input.search', { type: 'search', placeholder: 'Filter…' });
  search.addEventListener('input', () => { q = search.value.trim(); limit = STEP; draw(); });
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
  // One credit line per licence and author: attribution licences need every author named, not one per licence.
  const credits = [...new Map(refs.map(r => [`${r.licence_id}|${r.author || ''}`, r])).values()];
  const tbody = h('tbody');
  let q = '';
  const open = r => {
    const parts = [r.author, r.url ? h('a', { href: r.url, target: '_blank', rel: 'noopener' }, 'source') : null, r.licence, r.credit].filter(Boolean);
    const credit = h('p.small', null, parts.flatMap((x, i) => (i ? [' · ', x] : [x])));
    const head = h('div.detailhead', null, h('h3', null, `${r.note || ''} ${r.dynamic || ''}`), h('span.small.muted', null, `${r.instrument || ''} · ${r.kind || ''}${r.piece ? ` · ${r.piece}` : ''}`));
    panel.classList.add('picked'); // the empty prompt does not stick over the table; a picked recording does
    if (!r.audio) {
      panel.replaceChildren(head, h('div.note', null, r.withheld || 'This recording\'s licence does not allow it to be played here; follow the source link to hear it.'), credit);
      return;
    }
    const player = new Player([{ rel: r.audio, label: 'Recording', side: 'real' }], { title: `${doc.title} ${r.note || ''} ${r.dynamic || ''}`.trim(), sub: `recording · ${r.licence}` });
    focusPlayer(player);
    const mapBox = h('div.maps', null, h('p.small.muted', null, 'Decoding…'));
    panel.replaceChildren(head, controls(player), mapBox, credit);
    player.ready.then(() => loadBuffer(r.audio)).then(buf => {
      if (!mapBox.isConnected) return; // another recording was picked while this one decoded
      const m = computeMap(buf);
      mapBox.replaceChildren(mapFrame(m.canvas, { fLo: m.fLo, fHi: m.fHi, duration: m.duration, player, label: 'Recording' }));
      player.play();
    }).catch(err => mapBox.replaceChildren(empty(`Could not decode: ${err.message}`)));
  };
  // The first 30 rows, the rest on request (the piano lists 129; a scroll box inside the page would fight the page's own
  // scrolling on a phone). The filter still looks through every recording.
  const STEP = 30;
  let limit = STEP;
  const moreBox = h('div.morebox');
  const draw = () => {
    const rx = q ? new RegExp(q.replace(/[.*+?^${}()|[\]\\]/g, '\\$&'), 'i') : null;
    const shown = refs.filter(r => !rx || rx.test([r.note, r.dynamic, r.kind, r.piece, r.instrument, r.author, r.licence].join(' ')));
    moreBox.replaceChildren();
    if (shown.length > limit) {
      const more = h('button.btn.small', { type: 'button' }, `Show all ${shown.length} recordings (${shown.length - limit} more)`);
      more.addEventListener('click', () => {
        const first = limit;
        limit = Infinity;
        draw();
        tbody.children[first]?.focus({ preventScroll: true });
      });
      moreBox.append(more);
    }
    tbody.replaceChildren(...shown.slice(0, limit).map(r => {
      // Source and Licence give way on a phone (the credits below list both); "link only" then shows by the note
      const tr = h('tr.clickable', null, h('td', null, r.note || '–', r.audio ? null : h('span.muted.phoneonly', null, ' · link only')),
        h('td', null, r.dynamic || '–'), h('td', null, r.kind || '–'),
        h('td', null, r.piece || ''), h('td.num', null, r.duration_s ? `${fmt.num(r.duration_s, 1)} s` : ''), h('td.small.wideonly', null, r.instrument || ''),
        h('td.small.wideonly', null, r.licence, r.audio ? null : h('span.muted', null, ' · link only')));
      const pick = () => {
        for (const x of tbody.children) x.classList.remove('on');
        tr.classList.add('on');
        open(r);
      };
      tr.tabIndex = 0;
      tr.addEventListener('click', pick);
      tr.addEventListener('keydown', e => {
        if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); pick(); }
      });
      return tr;
    }));
  };
  const search = h('input.search', { type: 'search', placeholder: 'Filter by note, dynamic, source…' });
  search.addEventListener('input', () => { q = search.value.trim(); limit = STEP; draw(); });
  draw();
  return h('div.refs', null,
    h('p.note', null, `${refs.length} recordings passed the screening (${cov.valid ?? refs.length} of ${cov.rows ?? '?'} candidate clips; the others failed checks such as signal-to-noise, a single clean note, decay or pitch). Each is credited with its licence below; recordings whose licence does not allow redistribution are listed with a link to their source but not played here.`),
    panel,
    h('div.toolbar', null, search),
    h('div.tablewrap', null, h('table.reftable', null, h('thead', null, h('tr', null, ['Note', 'Dynamic', 'Kind', 'Piece', 'Length', 'Source', 'Licence'].map((t, i) => h(i > 4 ? 'th.wideonly' : 'th', null, t)))), tbody)),
    moreBox,
    h('div.credits', null, h('h4', null, 'Credits'), h('ul.small', null, credits.map(r => h('li', null, `${r.licence}: `, [r.author, r.credit || r.instrument || r.source].filter(Boolean).join(' · '),
      r.licence_notes ? h('span.muted', null, ` (${r.licence_notes})`) : null)))));
}
