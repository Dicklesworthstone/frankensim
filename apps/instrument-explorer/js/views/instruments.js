// Instruments: where every family stands, grouped by section, with the numbers behind its badge.
import { h, getJSON, fmt, badge, splitBar, link, BADGES } from '../app.js';

const GATES = ['QP1', 'Q2', 'Q3', 'Q4', 'Q5'];
const GATE_TEXT = { QP1: 'speaks', Q2: 'alive', Q3: 'inside the ruler', Q4: 'blind test', Q5: 'phrases' };

export function title() {
  return 'Instruments';
}

export function ladder(milestones) {
  if (!milestones.length) return h('div.ladder.none', null, 'no qualification milestones defined');
  return h('div.ladder', { role: 'img', 'aria-label': `Qualification gates: ${milestones.map(m => `${m.gate} ${m.status.replace('_', ' ')}`).join(', ')}` }, GATES.map(g => {
    const m = milestones.find(x => x.gate === g);
    if (!m) return h('span.rung.absent', { title: `${g}: no milestone` }, g);
    const frac = m.total ? m.closed / m.total : 0;
    return h(`span.rung.st-${m.status}`, { title: `${m.code} (${GATE_TEXT[g]}): ${m.status.replace('_', ' ')}; ${m.closed} of ${m.total} prerequisite beads closed` },
      h('span.fill', { style: { width: `${100 * frac}%` } }), h('span.lbl', null, g));
  }));
}

export function familyCard(f) {
  const c = f.counts;
  const facts = [];
  if (c.graded) facts.push(fmt.plural(c.graded, 'graded note'));
  if (c.listening) facts.push(`${c.listening} ungraded`);
  if (c.refused) facts.push(`${c.refused} refused`);
  if (f.renders && !c.graded && !c.listening) facts.push(fmt.plural(f.renders, 'render'));
  facts.push(fmt.plural(f.references, 'recording'));
  return h('a.card.family', { href: link('f', f.id) },
    h('div.cardhead', null, h('h3', null, f.title), badge(f.badge, f.label)),
    f.phrases ? h('span.hears', null, '♪ ', fmt.plural(f.phrases, 'phrase'), ' to hear') : null,
    h('p.facts', null, facts.join(' · ')),
    c.inside + c.outside
      ? h('div.measure', null, splitBar(c.inside, c.outside), h('span.small', null, `${c.inside} of ${c.inside + c.outside} measures inside`))
      : null,
    f.summary && !f.summary_derived ? h('p.summary', null, f.summary) : null,
    ladder(f.milestones));
}

export async function render(root) {
  const idx = await getJSON('index.json');
  const t = idx.totals;
  const measured = t.measures_inside + t.measures_outside;
  const intro = h('section.intro', null,
    h('p.eyebrow', null, `${idx.families.length} instrument families`),
    h('h1', null, 'Where every instrument stands'),
    h('p.lede', null,
      `A badge sums up how a family's graded notes compare with real recordings of the same notes. A measure is `, h('b', null, 'inside'),
      ` when the simulation differs from the recording by no more than real recordings of the instrument differ from each other. `,
      `${t.graded_notes} notes are graded so far; ${t.measures_inside} of their ${measured} measures (${fmt.pct(t.measures_inside / measured)}) are inside.`),
    h('div.badgekey', null, Object.keys(BADGES).map(b => badge(b))),
    h('p.small.muted', null, 'The five boxes on each card are the qualification gates, QP1 to Q5. Each fills as the work it depends on is closed; green means the gate itself has passed. ',
      h('a', { href: '#/method' }, 'What the gates mean')));
  const sections = idx.sections.map(sec => {
    const fams = idx.families.filter(f => f.section === sec);
    if (!fams.length) return null;
    return h('section.group', null, h('h2', null, sec), h('div.grid', null, fams.map(familyCard)));
  });
  root.replaceChildren(intro, h('div.groups', null, sections));
}
