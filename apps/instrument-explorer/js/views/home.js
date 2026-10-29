// Listen: what the simulated instruments can play today. A featured player up top, then every phrase (scales, double
// stops, melodies, grooves, rows of notes) as cards you can play in place, then the way into the evidence.
import { h, getJSON, fmt, link, scrollMotion } from '../app.js';
import { phraseCard, phraseStage, FILTERS, pitchTally } from './phrase.js';

const FEATURED = [
  ['piano__melody', 'Piano'],
  ['kit__groove', 'Drum kit'],
  ['violin__melody', 'Violin'],
  ['trumpet__fanfare', 'Trumpet'],
];

// The gallery leads with what already sounds like the instrument and ends with the rows that show what does not yet.
const FAMILY_ORDER = ['piano', 'kit', 'violin', 'viola', 'cello', 'double-bass', 'trumpet', 'voice', 'tin-whistle', 'pan-pipes',
  'trombone', 'clarinet', 'flute', 'recorder'];
const KIND_ORDER = ['melody', 'groove', 'chords', 'scale', 'arpeggio', 'slur', 'vowels', 'row'];
const rank = (list, x) => (list.includes(x) ? list.indexOf(x) : list.length);

export function title() {
  return 'Listen';
}

let filter = 'all';

export async function render(root, r) {
  const [idx, ph] = await Promise.all([getJSON('index.json'), getJSON('phrases.json')]);
  const phrases = ph.phrases.slice().sort((a, b) => rank(FAMILY_ORDER, a.family) - rank(FAMILY_ORDER, b.family)
    || rank(KIND_ORDER, a.kind) - rank(KIND_ORDER, b.kind));
  const titles = ph.titles;
  const t = idx.totals;
  const measured = t.measures_inside + t.measures_outside;
  const byId = Object.fromEntries(phrases.map(p => [p.id, p]));
  const featured = FEATURED.filter(([id]) => byId[id]);

  // ---- hero
  const stageBox = h('div');
  const pick = h('div.segmented.small', { role: 'group', 'aria-label': 'Featured phrase' });
  const show = id => {
    for (const b of pick.children) {
      b.classList.toggle('on', b.dataset.id === id);
      b.setAttribute('aria-pressed', String(b.dataset.id === id));
    }
    const p = byId[id];
    stageBox.replaceChildren(featureCard(p, titles[p.family]));
  };
  for (const [id, label] of featured) {
    const b = h('button', { type: 'button', 'data-id': id }, label);
    b.addEventListener('click', () => show(id));
    pick.append(b);
  }
  const instrumentsHeard = new Set(phrases.map(p => p.family)).size;
  const hero = h('section.hero', null,
    h('div', null,
      h('p.eyebrow', null, 'FrankenSim · instrument tests'),
      h('h1', null, 'Instruments computed from ', h('em', null, 'physics'), ', played next to the real thing'),
      h('p.lede', null, 'Strings, bores, reeds, lips, plates and air, described from published measurements and solved in time. No recorded samples, nothing fitted to a recording. ',
        'Every clip below comes with its spectrogram, the notes it was asked to play, and how far it landed from them.'),
      h('p.lede.muted', null, `Some of it already sounds like the instrument; some of it plainly does not. Both are shown. ${t.graded_notes} single notes are graded against real recordings: ${fmt.pct(t.measures_inside / measured)} of their measures fall inside the spread between real instruments.`),
      h('div.actions', null,
        h('button.btn.primary', { type: 'button', onclick: () => {
          const to = document.getElementById('phrases');
          to.scrollIntoView({ behavior: scrollMotion(to.getBoundingClientRect().top) });
        } }, `Hear all ${phrases.length} phrases`),
        h('a.btn', { href: '#/instruments' }, 'Where each instrument stands'),
        h('a.btn', { href: '#/method' }, 'How we grade'))),
    h('div.feature', null, h('div.fhead', null, h('b', null, 'Featured'), pick), stageBox));
  if (featured.length) show(featured[0][0]);

  // ---- stats
  const inTune = phrases.reduce((a, p) => { const x = pitchTally(p); a.n += x.in + x.near + x.out; a.in += x.in; return a; }, { n: 0, in: 0 });
  const stats = h('div.stats', null,
    stat(phrases.length, `phrases from ${instrumentsHeard} instruments`),
    stat(`${inTune.in} / ${inTune.n}`, 'phrase notes within 10 cents of the note asked for'),
    stat(t.graded_notes, 'single notes graded against recordings'),
    stat(fmt.pct(t.measures_inside / measured), 'of graded measures inside the real spread'),
    stat(t.references, 'screened real recordings'));

  // ---- gallery
  const grid = h('div.pgrid');
  const filters = h('div.filters', { role: 'group', 'aria-label': 'Filter phrases' });
  const draw = () => {
    const kinds = (FILTERS.find(f => f[0] === filter) || FILTERS[0])[2];
    const shown = phrases.filter(p => !kinds || kinds.includes(p.kind));
    grid.replaceChildren(...shown.map(p => phraseCard(p, titles[p.family])));
    for (const b of filters.children) b.setAttribute('aria-pressed', String(b.dataset.f === filter));
  };
  for (const [id, label, kinds] of FILTERS) {
    const n = kinds ? phrases.filter(p => kinds.includes(p.kind)).length : phrases.length;
    if (!n) continue;
    const b = h('button', { type: 'button', 'data-f': id }, `${label} `, h('span.muted', null, String(n)));
    b.addEventListener('click', () => { filter = id; draw(); });
    filters.append(b);
  }
  draw();
  const gallery = h('section', { id: 'phrases' },
    h('div.sechead', null,
      h('div', null, h('h2', null, 'Hear them play'),
        h('p', null, 'Press play on any card. The coloured bars under a card are its notes: green within 10 cents of the note asked for, amber within 50, red further off. Open a phrase for its full spectrogram, note by note.')),
      filters),
    grid);

  // ---- way in
  const more = h('section', null,
    h('div.sechead', null, h('div', null, h('h2', null, 'The evidence behind it'),
      h('p', null, 'Every verdict on this site comes from a measurement you can open.'))),
    h('div.starthere', null,
      h('a.card', { href: link('f', 'violin', 'listen') }, h('b', null, 'Simulation and recording, note by note'),
        h('span', null, 'Violin G3 against a real violin: playback in lock-step, both spectrograms and their difference, a blind test, every measure against its ruler.')),
      h('a.card', { href: '#/method' }, h('b', null, 'How we grade'),
        h('span', null, 'The qualification gates, the 16 measures and how often each is passed, the rulers checked against a second library, and the score over time.')),
      h('a.card', { href: '#/decisions' }, h('b', null, 'Decisions'),
        h('span', null, 'What was decided and why: the first-principles pivot, materials and air as data, the grading changes, the findings that changed the plan.')),
      h('a.card', { href: '#/instruments' }, h('b', null, 'All 27 families'),
        h('span', null, 'Which instruments speak, which are graded, which are still blocked on geometry or data, and the qualification ladder of each.'))));
  const unknown = r && r.unknown
    ? h('div.note.warn', { role: 'status' }, h('p', null, `There is no page at ${r.unknown}, so this is the start page.`))
    : null;
  root.replaceChildren(...[unknown, hero, stats, gallery, more].filter(Boolean));
}

function featureCard(p, familyTitle) {
  const stage = phraseStage(p, familyTitle);
  stage.classList.add('compact');
  const tally = pitchTally(p);
  const total = tally.in + tally.near + tally.out;
  return h('div', null, stage,
    h('p.verdict', null,
      h('a.open', { href: link('p', p.id) }, h('span', null, p.real || p.ref ? 'Compare with the real one' : 'Full analysis'), h('span.arrow', { 'aria-hidden': 'true' }, '→')),
      h('span', null, total ? `${tally.in} of ${total} notes within 10 cents` : p.real ? 'switch to the recording while it plays' : 'one continuous render'),
      h('a', { href: link('f', p.family) }, `All of ${familyTitle}`)));
}

function stat(value, label) {
  return h('div.stat', null, h('b', null, String(value)), h('span', null, label));
}
