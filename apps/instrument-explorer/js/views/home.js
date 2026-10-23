// Home: what the program has produced so far, family by family.
import { h, getJSON, fmt, badge, splitBar, link, BADGES } from '../app.js';

const GATES = ['QP1', 'Q2', 'Q3', 'Q4', 'Q5'];
const GATE_TEXT = { QP1: 'speaks', Q2: 'recognisable', Q3: 'inside the ruler', Q4: 'blind test', Q5: 'phrases' };

function ladder(milestones) {
  if (!milestones.length) return h('div.ladder.none', null, 'no qualification milestones defined');
  return h('div.ladder', null, GATES.map(g => {
    const m = milestones.find(x => x.gate === g);
    if (!m) return h('span.rung.absent', { title: `${g}: no milestone` }, g);
    const frac = m.total ? m.closed / m.total : 0;
    return h(`span.rung.st-${m.status}`, { title: `${m.code} (${GATE_TEXT[g]}): ${m.status}; ${m.closed} of ${m.total} prerequisite beads closed` },
      h('span.fill', { style: { width: `${100 * frac}%` } }), h('span.lbl', null, g));
  }));
}

function card(f) {
  const c = f.counts;
  const facts = [];
  if (c.graded) facts.push(fmt.plural(c.graded, 'graded note'));
  if (c.listening) facts.push(`${c.listening} ungraded`);
  if (c.refused) facts.push(`${c.refused} refused`);
  if (f.renders && !c.graded && !c.listening) facts.push(fmt.plural(f.renders, 'render'));
  facts.push(fmt.plural(f.references, 'recording'));
  return h('a.card.family', { href: link('f', f.id) },
    h('div.cardhead', null, h('h3', null, f.title), badge(f.badge, f.label)),
    h('p.facts', null, facts.join(' · ')),
    c.inside + c.outside
      ? h('div.measure', null, splitBar(c.inside, c.outside), h('span.small', null, `${c.inside} of ${c.inside + c.outside} measures inside`))
      : null,
    f.summary && !f.summary_derived ? h('p.summary', null, f.summary) : null,
    ladder(f.milestones),
    h('p.small.muted', null, f.last_activity ? `Last activity ${fmt.date(f.last_activity)} · ${f.commits} commits` : `${f.commits} commits`));
}

export async function render(root) {
  const idx = await getJSON('index.json');
  const t = idx.totals;
  const measured = t.measures_inside + t.measures_outside;
  const anyClose = idx.families.some(f => f.badge === 'close');
  const intro = h('section.intro', null,
    h('h1', null, 'Simulated instruments, heard next to real ones'),
    h('p.lede', null,
      `Every sound on these pages is computed from physics: geometry, materials and air, with no recorded samples and nothing fitted to a recording. `,
      `${t.graded_notes} notes from ${idx.families.filter(f => f.counts.graded).length} families have been rendered and graded against real recordings of the same note. `,
      `Of their ${measured} measured differences, ${t.measures_inside} (${fmt.pct(t.measures_inside / measured)}) fall inside the spread between real recordings of the instrument; `,
      anyClose ? 'some families are inside on every measure.' : 'no family is inside on every measure yet.'),
    h('p.lede.muted', null,
      `${t.milestones_closed} of ${t.milestones_total} qualification milestones are closed. Many bodies are still simplified, and several renders are labelled diagnostic baselines. `,
      'Listen with that in mind: this is a working record, not a showcase.'),
    h('div.tiles', null,
      tile(t.graded_notes, 'notes graded against recordings'),
      tile(fmt.pct(t.measures_inside / measured), 'of measures inside the real spread'),
      tile(t.references, 'screened reference recordings'),
      tile(`${t.beads_closed} / ${t.beads_total}`, 'program beads closed'),
      tile(t.commits, `commits since ${fmt.date('2026-10-15')}`)),
    h('div.badgekey', null, Object.keys(BADGES).map(b => badge(b))));
  const sections = idx.sections.map(sec => {
    const fams = idx.families.filter(f => f.section === sec);
    if (!fams.length) return null;
    return h('section.group', null, h('h2', null, sec), h('div.grid', null, fams.map(card)));
  });
  root.replaceChildren(intro, ...sections);
}

function tile(value, label) {
  return h('div.tile', null, h('b', null, String(value)), h('span', null, label));
}
