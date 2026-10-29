// Process: how the work has gone, program-wide. Commits per day by instrument section, bead closes, the milestone
// matrix and the newest findings.
import { h, getJSON, fmt, chip, link } from '../app.js';
import { stackedBars, legend, SECTION_COLORS } from '../charts.js';

const GATES = ['QP1', 'Q2', 'Q3', 'Q4', 'Q5'];

export function title() {
  return 'Activity';
}

export async function render(root) {
  const [proc, idx] = await Promise.all([getJSON('process.json'), getJSON('index.json')]);
  const keys = [...proc.sections, 'Several families', 'Shared physics and tools'];
  const commitRows = proc.commits_by_day.map(d => ({ label: d.date, by: d.by }));
  const closeRows = proc.closes_by_day.map(d => ({ label: d.date, by: { closed: d.n } }));
  const b = proc.beads;
  const statusOrder = ['closed', 'in_progress', 'open', 'blocked', 'deferred'];
  const statusTotal = statusOrder.reduce((a, k) => a + (b[k] || 0), 0);
  root.replaceChildren(
    h('section.intro', null, h('h1', null, 'How the work has gone'),
      h('p.lede', null, `Since ${fmt.date(proc.since)}: ${proc.commits_total} commits to main (${proc.commits_irp} tagged to the instrument-realism program) and ${proc.beads_total} program beads, ${b.closed || 0} of them closed. `,
        'A commit subject in this repository states one finding, so the feed at the bottom reads as a lab notebook.')),
    h('section', null, h('h2', null, 'Commits per day, by the instrument section they name'),
      stackedBars(commitRows, keys, SECTION_COLORS, { yLabel: 'commits' }), legend(keys, SECTION_COLORS),
      h('p.small.muted', null, 'A commit counts toward a section when its subject or bead names a family of that section; "Shared physics and tools" names none (solvers, materials, meshing, infrastructure).')),
    h('section', null, h('h2', null, 'Beads closed per day'),
      stackedBars(closeRows, ['closed'], { closed: '#2f9e44' }, { height: 160, yLabel: 'closes' }),
      // the bar carries counts only where a segment has room for them; the caption names every status
      h('div.statusbar', null, statusOrder.filter(k => b[k]).map(k => h(`span.st-${k}`, { style: { flex: String(b[k]) }, title: `${b[k]} ${k.replace('_', ' ')}` },
        b[k] / statusTotal >= 0.06 ? String(b[k]) : ''))),
      h('p.small.muted', null, `${statusTotal} program beads by status: ${statusOrder.filter(k => b[k]).map(k => `${b[k]} ${k.replace('_', ' ')}`).join(', ')}. Labels: ${Object.entries(proc.labels).map(([k, v]) => `${v} ${k}`).join(', ')}. A bead is closed only after its evidence is verified on the build host.`)),
    h('section', null, h('h2', null, 'Qualification milestones'),
      h('p.small.muted', null, 'Rows are families, columns the gates from "speaks" to "plays phrases". Each cell shows how much of the work the milestone depends on is closed. No milestone is qualified until it closes.'),
      h('div.tablewrap', null, h('table.matrix', null,
        h('thead', null, h('tr', null, h('th', null, 'Family'), GATES.map(g => h('th', null, g)))),
        h('tbody', null, idx.families.filter(f => f.milestones.length).map(f => h('tr', null,
          h('th', null, h('a', { href: link('f', f.id, 'qualification') }, f.title)),
          GATES.map(g => {
            const m = f.milestones.find(x => x.gate === g);
            if (!m) return h('td.absent', null, '');
            const frac = m.total ? m.closed / m.total : 0;
            return h(`td.st-${m.status}`, { title: `${m.code}: ${m.status}, ${m.closed}/${m.total} prerequisites closed` },
              // a narrow phone may break "37/122" after the slash rather than push Q5 off the screen
              h('div.cellbar', null, h('span', { style: { width: `${100 * frac}%` } })), h('span.small', null, `${m.closed}/`, h('wbr'), `${m.total}`));
          }))))))),
    h('section', null, h('h2', null, 'Newest findings'),
      h('ol.timeline', null, proc.feed.map(c => h('li.commit', null, h('span.when', null, fmt.datetime(c.when)), chip(c.type, `t-${c.type}`),
        h('span.subj', null, c.subject), h('span.mono.small.muted', null, ` ${c.sha}`))))));
}
