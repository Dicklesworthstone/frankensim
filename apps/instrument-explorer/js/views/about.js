// About: what the badges and measures mean, where every number comes from, and how to rebuild.
import { h, badge, BADGES } from '../app.js';

export function title() {
  return 'About';
}

export async function render(root) {
  root.replaceChildren(h('article.about', null,
    h('h1', null, 'About this explorer'),
    h('p', null, 'FrankenSim renders musical instruments from first principles: the geometry of the instrument, the measured properties of its materials and of the air, and the physics that couples them. Nothing is sampled, and no parameter is tuned to make a render resemble a recording. This explorer puts those renders next to real recordings of the same notes, with the grading, the qualification work and the history behind them.'),
    h('h2', null, 'What the badges mean'),
    h('ul.badgelist', null, Object.keys(BADGES).map(b => h('li', null, badge(b), ' ', {
      close: 'every graded measure of every graded note is inside the spread between real recordings.',
      mixed: 'at least half of the graded measures are inside, on both the loose and the same-note ruler.',
      off: 'fewer than half of the graded measures are inside.',
      listening: 'renders exist, but there is no recording of the same note to grade them against (recorder, tin whistle, pan pipes).',
      ungraded: 'renders exist but none is graded: the grader refused them (wrong note, no signal) or they are work in progress.',
      none: 'nothing has been rendered yet; the recordings the grader holds are listed.',
    }[b]))),
    h('h2', null, 'How a note is graded'),
    h('p', null, 'The grader (instrgrade decompose) aligns the render with a University of Iowa recording of the same note and dynamic, matches their loudness, and measures 15 or 16 differences: spectral and envelope distances, spectrogram and scalogram differences, partial frequencies and levels, inharmonicity, odd/even balance, noise, attack, brightness and its motion, and decay for struck and plucked strings. The ruler for each measure is how far real recordings of the instrument differ from each other on that measure: another dynamic, another take, the neighbouring notes. A measure is inside when the render is no further from the recording than the largest of those real-vs-real distances; the strict ruler uses only the same note. As a check of the rulers, recordings from a second library (TinySOL) are scored the same way; most of them fall inside.'),
    h('p', null, 'Inside the ruler is necessary, not sufficient: a render can pass every measure and still sound synthetic. Blind listening (Q4) is the gate that answers that, and the blind test on each Listen page is a private version of it.'),
    h('h2', null, 'Where the data comes from'),
    h('ul', null,
      h('li', null, 'Renders, recordings at the grading segment, spectrogram and scalogram images and grades: the grading tree on the analysis host (frankensim.instrument-pair.v1), mirrored by tools/sync.sh. Earlier gradings of a note are kept and shown under History.'),
      h('li', null, 'Family descriptions, the reference screening, the physics claims register, the violin body-mode tables and the bow sweeps: files committed to the repository, read at origin/main (never the shared working tree).'),
      h('li', null, 'Milestones, bead closes and progress: the beads tracker export. Timelines: git history of origin/main.'),
      h('li', null, 'Work in progress (drum renders, the assembled violin): renders by lanes that have not been committed or graded yet. They are labelled as such.')),
    h('h2', null, 'Licences'),
    h('p', null, 'No recording is ever committed to the repository. The University of Iowa Musical Instrument Samples are free to use without restriction. Recordings under CC BY are played with their credit, shown in each instrument\'s References section. Recordings whose licence does not allow redistribution (YouTube videos under the standard licence) were used only to grade: they are listed with their source link and are not played here.'),
    h('h2', null, 'Rebuilding'),
    h('pre', null, 'apps/instrument-explorer/tools/sync.sh           # mirror renders, grades and recordings\npython3 apps/instrument-explorer/tools/build.py  # build the site data\npython3 apps/instrument-explorer/tools/serve.py  # http://localhost:8742')));
}
