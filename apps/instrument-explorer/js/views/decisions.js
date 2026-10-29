// Decisions: the choices that steer the instrument program (owner decisions, rulings, findings, method changes), newest
// first, and the owner-decision beads as the tracker has them.
import { h, getJSON, fmt, chip, stack } from '../app.js';

const KINDS = {
  owner: ['Owner decision', 'k-owner'],
  ruling: ['Ruling', 'k-ruling'],
  finding: ['Finding', 'k-finding'],
  method: ['Grading method', 'k-method'],
};

export function title() {
  return 'Decisions';
}

export async function render(root) {
  const doc = await getJSON('decisions.json');
  const entries = (doc.entries || []).slice().reverse();
  const shown = new Set(Object.keys(KINDS));
  const list = h('ol.decisions');
  const draw = () => list.replaceChildren(...entries.filter(e => shown.has(e.kind)).map(e => {
    const [label, cls] = KINDS[e.kind] || [e.kind, ''];
    return h(`li.decision.${cls}`, null,
      h('div.dhead', null, h('span.when', null, fmt.date(e.date)), chip(label, cls), h('h3', null, e.title)),
      h('dl.dbody', null,
        e.decided ? [h('dt', null, 'Decided'), h('dd', null, e.decided)] : null,
        e.why ? [h('dt', null, 'Why'), h('dd', null, e.why)] : null,
        e.effect ? [h('dt', null, 'What it changed'), h('dd', null, e.effect)] : null,
        e.refs && e.refs.length ? [h('dt', null, 'Record'), h('dd.mono.small', null, e.refs.join(' · '))] : null));
  }));
  const filters = h('div.segmented.small', null, Object.entries(KINDS).map(([k, [label]]) => {
    const b = h('button.on', { type: 'button' }, label);
    b.addEventListener('click', () => {
      if (shown.has(k)) shown.delete(k); else shown.add(k);
      b.classList.toggle('on', shown.has(k));
      draw();
    });
    return b;
  }));
  draw();
  const beads = doc.owner_beads || [];
  root.replaceChildren(h('article.decisions-page', null,
    h('section.intro', null, h('h1', null, 'Decisions'),
      h('p.lede', null, 'What was decided, why, and what it changed, newest first. This is the reasoning behind the renders: why some families are labelled diagnostic baselines, why the public site waits, why a measure was retired.'),
      h('p.small.muted', null, doc.about || '')),
    h('section', null, h('div.toolbar', null, h('span.small.muted', null, 'Show:'), filters), list),
    beads.length ? h('section', null, h('h2', null, 'Questions waiting on the owner'),
      h('p.small.muted', null, 'Owner-decision beads in the tracker. Open ones still apply their plan default until answered.'),
      h('div.tablewrap', null, stack(h('table.metrics', null,
        h('thead', null, h('tr', null, ['Bead', 'Question', 'Status', 'Answer'].map(t => h('th', null, t)))),
        h('tbody', null, beads.map(b => h('tr', null,
          h('td.mono.small', null, b.id.replace(/^frankensim-irp-v5-instrument-realism-/, '')),
          h('td', null, b.title),
          h('td', null, b.status === 'closed' ? chip(`decided ${fmt.date(b.closed)}`, 'in') : chip(b.status.replace('_', ' '), 'near')),
          h('td.small', null, b.reason || '–'))))), ['', '', '', 'Answer']))) : null));
}
