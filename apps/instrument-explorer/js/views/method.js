// How we grade: the target, the gate ladder, how one note is graded, every measure with how often the simulations pass it,
// the rulers' own check, the score over time, the pitch of every playable note, and what the grader cannot see.
import { h, getJSON, fmt, chip, link, splitBar, DATA, stack } from '../app.js';
import { rulerBar, historyChart } from '../charts.js';

// What each grader measure compares, and the kind of defect it catches.
const MEASURES = {
  spec_db: ['The average level difference between the two notes\' long-term spectra.', 'the overall tone colour: too bright, too dull, a missing band'],
  env_db: ['Level differences in third-octave bands, roughly how the ear groups frequencies.', 'body and bore resonance regions that are too loud or too quiet'],
  harm_db: ['The level of each harmonic relative to the fundamental, harmonic by harmonic.', 'the wrong excitation spectrum (bow, reed, lips, hammer) or radiation'],
  stft_lsd_db: ['The spectrograms compared frame by frame (10 ms frames, log frequency).', 'anything that differs in time and frequency at once: onsets, decays, wobble'],
  cwt_lsd_db: ['The same comparison with a wavelet transform, sharper in time at high frequency.', 'attack transients and low-frequency detail the spectrogram smears'],
  cepstral_env_db: ['The smooth outline of the spectrum, independent of the pitch.', 'the shape of the body or bore resonances (the "formants" of the instrument)'],
  partial_traj_db: ['How each partial\'s level moves over the note.', 'partials that swell, decay or beat differently from the real instrument'],
  partial_freq_cents: ['The frequency of each partial against the recording\'s, in cents.', 'wrong stiffness, mistuned resonances, a partial series that does not line up'],
  f0_cents: ['The pitch of the fundamental, in cents.', 'the wrong note or the wrong tuning'],
  inharmonicity_B: ['The stiffness coefficient B fitted to the partial series f_n = n f0 sqrt(1 + B n^2).', 'a string or bar that is too stiff or too floppy (strings and bars only)'],
  odd_even_db: ['The energy in odd harmonics against even ones.', 'a wrong bore or pipe character, a wrong bowing or plucking point'],
  hnr_db: ['The energy in the harmonics against the noise between them.', 'a tone that is too clean (no breath or bow noise) or too noisy'],
  attack_log2: ['The time the note takes to rise from 10 % to 90 %, compared in octaves of time (twice as long = 1).', 'onsets that are too abrupt or too slow'],
  decay_db_s: ['How fast each partial dies away, in dB per second.', 'a struck or plucked note that rings too long or dies too fast'],
  centroid_st: ['The brightness (spectral centroid) over the note, in semitones.', 'brightness that is wrong or does not move the way a real note opens up'],
  flux_db: ['How much the spectrum changes from frame to frame.', 'a frozen, synthetic steady state; this is also the Q2 "alive" measure'],
};

const GATES = [
  ['QP1', 'Passive chain and speaks', 'The instrument\'s passive behaviour matches measurements (input impedance per fingering, bridge admittance, modal ratios); it sounds the intended note within 10 cents, with thresholds inside the published range.'],
  ['Q2', 'Alive', 'A sustained note\'s spectral flux after the attack lies inside the spread of real recordings of its instrument, on at least 80 % of graded notes. Real recordings from a second library must pass and five lifeless null signals must fail in the same run.'],
  ['Q3', 'Inside the ruler', 'Every bound measure is inside the real-vs-real spread, and the way the sound changes from pp to ff and along a scale matches.'],
  ['Q4', 'Owner blind test', 'A balanced two-choice listening test: the owner cannot tell the simulation from the recording better than chance (both confidence bounds inside 40-60 %). A spectrogram-pair test runs beside it.'],
  ['Q5', 'Phrases', 'Scales, slurs, articulation and dynamics changes pass Q3 and Q4.'],
];

const BLIND_SPOTS = [
  'Analysis windows of 0.30-1.10 s (sustained) and 0.02-0.52 s (decaying) leave out the attack and the late decay.',
  'One straight-line decay fit cannot see a double decay or beating.',
  'Harmonics-to-noise counts stiff partials (B·k³ above about 0.06) as noise, and cannot tell pitch-synchronous noise from steady noise.',
  'A 30-50 Hz high-pass removes the fundamentals of the tuba and the double bass.',
  'The University of Iowa recordings are not level-calibrated, so absolute loudness is matched, never verified.',
  'Inside the ruler is necessary, not sufficient: a render can pass every measure and still sound synthetic. Only the blind test (Q4) answers that.',
];

function pitchChip(n) {
  if (n.cents_off === null || n.cents_off === undefined) return chip(n.sounds || 'silent', 'na');
  const a = Math.abs(n.cents_off);
  if (a <= 10) return chip('in tune', 'in');
  if (a <= 50) return chip(`${n.cents_off > 0 ? '+' : ''}${Math.round(n.cents_off)} c`, 'near');
  return chip('wrong note', 'out');
}

function exampleFigure(p) {
  const fig = p && p.figures && p.figures.spectrogram;
  if (!fig) return null;
  const frame = (side, label) => fig.img[side]
    ? h('figure.mfig', null, h('img', { src: DATA + fig.img[side], alt: label, loading: 'lazy' }), h('figcaption', null, label))
    : null;
  return h('div.mfigs', null, frame('ours', 'Simulation'), frame('real', 'Recording'), frame('diff', 'Difference (red: simulation louder, blue: quieter)'));
}

export function title() {
  return 'How we grade';
}

export async function render(root) {
  const [proc, idx] = await Promise.all([getJSON('process.json'), getJSON('index.json')]);
  let example = null;
  try {
    const violin = await getJSON('family/violin.json');
    example = violin.pairs.find(p => p.kind === 'graded') || null;
  } catch (err) {
    example = null;
  }
  const exIn = example && example.metrics.find(m => m.inside === true && m.ruler.length > 2);
  const exOut = example && example.metrics.find(m => m.inside === false && m.ruler.length > 2);

  const fams = idx.families.filter(f => f.milestones.length);
  const gateRows = GATES.map(([code, name, test]) => {
    const ms = fams.map(f => f.milestones.find(m => m.gate === code)).filter(Boolean);
    const closed = ms.filter(m => m.status === 'closed').length;
    const pre = ms.reduce((a, m) => a + m.closed, 0);
    const preT = ms.reduce((a, m) => a + m.total, 0);
    return h('tr', null, h('th', null, code), h('td', null, h('b', null, name)), h('td', null, test),
      h('td.num', null, `${closed} of ${ms.length} families`, h('div.small.muted', null, `${pre} of ${preT} prerequisite beads closed`)));
  });

  const ps = proc.program_score || [];
  const last = ps[ps.length - 1];
  const series = [{ name: 'All families', points: ps.map(p => ({ t: new Date(p.when.replace(/Z?$/, 'Z')), y: p.inside / p.measured,
    label: `${fmt.datetime(p.when)}: ${p.inside} of ${p.measured} inside over ${p.notes} notes` })) }];
  const famSeries = Object.entries(proc.scores || {}).map(([fid, pts]) => ({ name: proc.titles[fid] || fid,
    points: pts.map(p => ({ t: new Date(p.when.replace(/Z?$/, 'Z')), y: p.inside / p.measured,
      label: `${fmt.datetime(p.when)}: ${p.inside} of ${p.measured} inside over ${p.notes} notes` })) }));

  const pitchRows = (proc.pitch || []).map(f => {
    const ok = f.notes.filter(n => n.cents_off !== null && Math.abs(n.cents_off) <= 10).length;
    const scored = f.notes.filter(n => n.target).length;
    return h('section.pitchfam', null,
      h('h3', null, h('a', { href: link('f', f.id) }, f.title), h('span.small.muted', null, scored ? ` ${ok} of ${scored} in tune` : ' no target pitch: the lip finds its own partial'),
        f.sha ? h('span.small.muted.mono', null, ` · commit ${f.sha.slice(0, 9)}`) : null),
      // a note with no target is a render id (hold-t5e-1), so it is numbered instead
      h('div.pitchnotes', null, f.notes.map((n, i) => h('div.pn', { title: n.f0_hz ? `${n.f0_hz} Hz` : 'no pitch' },
        h('b', null, n.target ? n.label : `Held note ${i + 1}`),
        h('span', null, n.sounds ? `sounds ${n.sounds}${n.sounds_alt ? ` or ${n.sounds_alt}` : ''}` : 'silent'), pitchChip(n)))));
  });

  root.replaceChildren(h('article.method', null,
    h('section.intro', null,
      h('h1', null, 'How we grade'),
      h('p.lede', null, 'The target is the owner\'s: "a spectrogram of our purely synthetic/simulated instrument is basically indistinguishable to a well-recorded sample of a real instrument." Every render is put next to a real recording of the same note and measured. The tolerance on each measure comes from real recordings themselves, not from a number we chose.'),
      h('div.note.warn', null,
        h('p', null, h('b', null, 'The grader checks the work; it does not steer it. '),
          'Until 2026-10-21 the program did hill-climb the grader: find the measure that is furthest off, add a mechanism, re-grade. The owner stopped that as "hill climbing with citations", because a mechanism added to close a gap can raise the score without making the instrument more real. Since then the work is driven by what the real instrument has that the model lacks (geometry, materials, the player). The grader is the final check. ',
          h('a', { href: '#/decisions' }, 'See the decisions.')))),

    h('section', null, h('h2', null, 'The ladder every family climbs'),
      h('p.small.muted', null, 'A family qualifies at a gate only when the gate\'s milestone bead closes, and that needs every bead it depends on closed and verified on the build host.'),
      h('div.tablewrap', null, stack(h('table.metrics.gates', null,
        h('thead', null, h('tr', null, ['Gate', '', 'What it takes', 'Where the program is'].map(t => h('th', null, t)))),
        h('tbody', null, gateRows)), ['', '', '', '']))),

    h('section', null, h('h2', null, 'How one note is graded'),
      h('ol.steps', null,
        h('li', null, h('b', null, 'Pair. '), 'The render is paired with a University of Iowa recording of the same note and dynamic (another library where Iowa has none).'),
        h('li', null, h('b', null, 'Align and match loudness. '), 'Both are cut to the same steady segment and given one shared gain, so neither is louder.'),
        h('li', null, h('b', null, 'Measure. '), `Up to ${Object.keys(MEASURES).length} differences between them, listed below.`),
        h('li', null, h('b', null, 'Build the ruler. '), 'The same differences are measured between real recordings of the instrument: other dynamics, other takes, the neighbouring notes. Their largest value is the ruler for that measure. The strict ruler uses the same note only.'),
        h('li', null, h('b', null, 'Judge. '), 'A measure is inside when the simulation is no further from the recording than real recordings are from each other.')),
      example ? h('div.example', null,
        h('p', null, `A real example: the violin's ${example.note} ${example.dynamic || ''}, graded ${fmt.datetime((example.history.slice(-1)[0] || {}).when)}. `,
          h('a', { href: link('f', 'violin', 'listen', example.id) }, 'Listen to it.')),
        exampleFigure(example),
        h('table.metrics', null, h('tbody', null, [exIn, exOut].filter(Boolean).map(m => h('tr', null,
          h('td', null, m.label), h('td.num', null, `${fmt.num(m.gap, 2)} ${m.unit}`),
          h('td', null, rulerBar(m)), h('td', null, m.inside ? chip('inside', 'in') : chip('outside', 'out')))))),
        h('p.small.muted', null, 'Green band: the spread between real recordings, from zero to the largest real-vs-real distance. Grey ticks: each real-vs-real distance. Dot: the simulation\'s distance from the recording. Dashed line: the same-note ruler.')) : null),

    h('section', null, h('h2', null, 'The measures, and how often the simulations pass them'),
      h('p.small.muted', null, `Every current graded note of every family, hardest measure first. "Typical miss" is the median of the simulation's distance divided by the ruler: 1.0 is the edge of the real spread.`),
      h('div.tablewrap', null, stack(h('table.metrics.measures', null,
        h('thead', null, h('tr', null, ['Measure', 'What it compares', 'What it catches', 'Inside', 'Typical miss', 'Missed most by'].map(t => h('th', null, t)))),
        h('tbody', null, (proc.measures || []).map(m => {
          const d = MEASURES[m.id] || ['', ''];
          const top = Object.entries(m.missed_by).slice(0, 3).map(([k, v]) => `${k} ${v}`).join(', ');
          return h('tr', null, h('td', null, h('b', null, m.label)), h('td.small', null, d[0]), h('td.small', null, d[1]),
            h('td.num', null, splitBar(m.inside, m.graded - m.inside), h('span.small', null, `${m.inside} of ${m.graded}`)),
            h('td.num', null, fmt.num(m.median_ratio, 2)), h('td.small.muted', null, top || 'none'));
        }))), ['', 'Compares', 'Catches', 'Inside', 'Typical miss', 'Missed most by']))),

    h('section', null, h('h2', null, 'Is the ruler fair?'),
      h('p.small.muted', null, 'Real recordings from a second library (TinySOL, IRCAM: another room, microphone and player) are scored on exactly the same rulers. If the rulers were too strict, real instruments would fail them too.'),
      h('div.tablewrap', null, stack(h('table.metrics.calib', null,
        h('thead', null, h('tr', null, ['Family', 'Real recordings inside', 'Our renders inside', 'Recordings checked'].map(t => h('th', null, t)))),
        h('tbody', null, (proc.calibration || []).map(c => h('tr', null,
          h('td', null, h('a', { href: link('f', c.id, 'grades') }, c.title)),
          h('td.num', null, h('div.cbar', null, h('span.real', { style: { width: `${100 * (c.real || 0)}%` } })), fmt.pct(c.real)),
          h('td.num', null, h('div.cbar', null, h('span.ours', { style: { width: `${100 * (c.ours || 0)}%` } })), fmt.pct(c.ours)),
          h('td.num', null, String(c.recordings || '–')))))), ['', 'Real recordings inside', 'Our renders inside', 'Recordings checked']))),

    h('section', null, h('h2', null, 'The score over time'),
      last ? h('p.small.muted', null, `Share of all graded measures inside the real spread after each grading, every note counted at its latest grading. Now: ${last.inside} of ${last.measured} (${fmt.pct(last.inside / last.measured)}) over ${last.notes} notes. The set of graded notes changes as families are added, so a step can come from coverage as well as quality. The newest grading is from ${fmt.datetime(last.when)}: renders made after it (the assembled violin, the playable notes, the kit) are not graded yet.`) : null,
      ps.length ? historyChart(series, { height: 220 }) : h('p.empty', null, 'No grading history yet.'),
      famSeries.length ? h('details.famscores', null, h('summary', null, 'By family'), historyChart(famSeries, { height: 320 })) : null),

    h('section', null, h('h2', null, 'Does it play the right note?'),
      h('p.small.muted', null, `Every note each family can play today (each family's renders name the commit they were made at), with the pitch it actually sounds, measured by ${proc.pitch_method || 'autocorrelation'}. Q1 asks for the intended note within 10 cents. In tune: within 10 cents. Amber: the right note, out of tune. Red: a different note.`),
      pitchRows.length ? pitchRows : h('p.empty', null, 'No pitch check in this build.')),

    h('section', null, h('h2', null, 'What the grader cannot see'),
      h('ul', null, BLIND_SPOTS.map(t => h('li', null, t))))));
}
