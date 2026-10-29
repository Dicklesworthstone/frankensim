// A phrase (scale, double stops, melody, groove, a row of notes): its player, its spectrogram with the notes it asks
// for drawn on it, the measured pitch of every note, and how it was made. Shared by the home gallery and family pages.
import { h, fmt, link, DATA } from '../app.js';
import { Player, controls, playButton } from '../audio.js';
import { mapFrame, pitchClass, unpitched } from '../spectro.js';

export const KIND_LABEL = {
  scale: 'scale', arpeggio: 'arpeggio', slur: 'lip slur', chords: 'chords', melody: 'melody', groove: 'groove',
  vowels: 'vowels', row: 'every note',
};

export const FILTERS = [
  ['all', 'All'],
  ['scales', 'Scales and arpeggios', ['scale', 'arpeggio', 'slur']],
  ['chords', 'Chords', ['chords']],
  ['melody', 'Melodies and grooves', ['melody', 'groove']],
  ['rows', 'Rows of notes', ['row', 'vowels']],
];

const players = new Map();

/** One player per phrase, shared by every view that shows it (so a card and its page stay in step). */
export function phrasePlayer(p, familyTitle) {
  if (!players.has(p.id)) {
    const sources = [{ rel: p.audio, label: 'Simulation', side: 'sim' }];
    if (p.real && p.real.audio) sources.push({ rel: p.real.audio, label: 'Recording', side: 'real' });
    else if (p.ref && p.ref.audio) sources.push({ rel: p.ref.audio, label: 'Real notes', side: 'real' });
    players.set(p.id, new Player(sources, { title: p.title, sub: `${familyTitle || p.family} · ${KIND_LABEL[p.kind] || p.kind}`, duration: p.duration_s }));
  }
  return players.get(p.id);
}

/** Pitch verdict counts for a phrase's notes that have a target. */
export function pitchTally(p) {
  const t = { in: 0, near: 0, out: 0, refused: 0 };
  // A double stop has no single pitch to check, so a refused chord is silent, not a wrong note.
  const checked = p.kind !== 'chords';
  for (const n of p.notes || []) {
    const c = pitchClass(n);
    if (c === 'refused') {  // asked for and not played: a miss, counted among the wrong notes and named apart
      if (checked) t.out += 1;
      t.refused += 1;
    } else if (c !== 'none') t[c] += 1;
  }
  return t;
}

function centsText(n) {
  if (n.refused) return 'refused';
  if (unpitched(n)) return 'no steady pitch';
  if (n.cents === null || n.cents === undefined) return n.f0_hz ? `${fmt.num(n.f0_hz, 0)} Hz` : '–';
  const c = Math.round(n.cents);
  if (Math.abs(c) >= 600) {
    const semis = Math.round(c / 100);
    return `${semis > 0 ? '+' : '−'}${Math.abs(semis)} st`;
  }
  return `${c > 0 ? '+' : c < 0 ? '−' : '±'}${Math.abs(c)}¢`;
}

/** What a card says when no note of the phrase has a pitch to check (chords, a single render). */
function untimedText(p) {
  const n = (p.notes || []).length;
  if (p.native || !n) return 'one continuous render';
  const refused = p.notes.filter(x => x.refused).length;
  return fmt.plural(n, p.kind === 'chords' ? 'chord' : 'part') + (refused ? `, ${refused} refused by the solver` : '');
}

/** A gallery card: the phrase's map as a thumbnail you can play in place, and a link to its page. */
export function phraseCard(p, familyTitle) {
  const player = phrasePlayer(p, familyTitle);
  const img = h('img', { src: DATA + p.map, alt: `Spectrogram of ${p.title}`, loading: 'lazy', decoding: 'async' });
  const head = h('div.playhead');
  const thumb = h('div.thumb', null, img, head, playButton(player, { small: true, label: p.title }),
    h('span.dur', null, `${fmt.num(p.duration_s, 0)} s`));
  thumb.addEventListener('click', () => player.toggle());
  const t = pitchTally(p);
  const total = t.in + t.near + t.out;
  const page = link('p', p.id);
  const card = h('article.pcard', null, thumb,
    h('div.body', null,
      h('div.fam', null, h('b', null, familyTitle || p.family), h('span.chip.kind', null, KIND_LABEL[p.kind] || p.kind)),
      h('h3', null, h('a', { href: page }, p.title)),
      total
        ? h('div.pitchline', null, h('span.pitchdots', { 'aria-hidden': 'true' }, (p.notes || []).filter(n => pitchClass(n) !== 'none').slice(0, 24).map(n => h(`i.${pitchClass(n)}`))),
          h('span', null, t.in === total ? `all ${total} notes in tune` : `${t.in} of ${total} in tune${t.out - t.refused ? `, ${t.out - t.refused} wrong` : ''}${t.refused ? `, ${t.refused} refused` : ''}`))
        : h('div.pitchline', null, h('span', null, untimedText(p))),
      h('a.open', { href: page, 'aria-label': `Full analysis of ${p.title}` },
        h('span', null, p.real || p.ref ? 'Compare with the real one' : 'Full analysis'), h('span.arrow', { 'aria-hidden': 'true' }, '→'))));
  player.on(q => {
    card.classList.toggle('playing', q.playing);
    head.style.left = `${(100 * Math.min(q.position, p.duration_s)) / p.duration_s}%`;
    head.style.opacity = q.playing || q.position > 0 ? '1' : '0';
  }, card);
  return card;
}

/** The full phrase view: transport, map(s) with note overlays, the note strip, the pitch summary and the method. */
export function phraseStage(p, familyTitle) {
  const player = phrasePlayer(p, familyTitle);
  player.load().catch(() => {});
  const fig = { fLo: p.f_lo_hz, fHi: p.f_hi_hz, duration: p.duration_s, player };
  const maps = h('div.maps', null,
    mapFrame(h('img', { src: DATA + p.map, alt: `Spectrogram of the simulation: ${p.title}` }), { ...fig, notes: p.notes, label: 'Simulation', side: 'sim', tall: true }));
  if (p.real && p.real.map) {
    maps.append(mapFrame(h('img', { src: DATA + p.real.map, alt: `Spectrogram of the recording: ${p.title}` }),
      { ...fig, duration: p.real.duration_s || p.duration_s, label: 'Recording', side: 'real', tall: true }));
  }
  const notes = p.notes || [];
  const strip = notes.length ? h('div.notestrip', { role: 'group', 'aria-label': 'Notes and their measured pitch; each plays from its start' }, notes.map(n => {
    const c = pitchClass(n);
    const chip = h(`button.nchip.${c}`, { type: 'button', title: noteTitle(n), 'aria-label': `${noteTitle(n)}: ${centsText(n)}` },
      h('b', null, n.note), h('span', null, centsText(n)));
    chip.addEventListener('click', () => player.load().then(() => { player.seek(n.at); player.play(); }, () => {}));
    return chip;
  })) : null;
  if (strip) {
    const chips = [...strip.children];
    let last = -1;
    player.on(q => {
      const i = q.playing || q.position ? notes.findIndex(n => q.position >= n.at && q.position < n.at + n.dur) : -1;
      if (i === last) return;
      if (last >= 0 && chips[last]) chips[last].classList.remove('now');
      if (i >= 0 && chips[i]) {
        chips[i].classList.add('now');
        if (q.playing) {
          // keep the sounding note in view by scrolling the strip itself, never the page under the reader
          const c = chips[i].getBoundingClientRect();
          const s = strip.getBoundingClientRect();
          if (c.left < s.left || c.right > s.right) strip.scrollLeft += c.left - s.left - (s.width - c.width) / 2;
        }
      }
      last = i;
    }, strip);
  }
  const t = pitchTally(p);
  const total = t.in + t.near + t.out;
  const legend = h('div.legendrow', null,
    notes.some(n => n.target_hz) ? [
      h('span', null, h('i.c-in'), 'note asked for, played within 10¢'),
      h('span', null, h('i.c-near'), 'within 50¢'),
      h('span', null, h('i.c-out'), 'further off'),
      notes.some(n => pitchClass(n) === 'out') ? h('span', null, h('i.solid.c-out'), 'where it actually sounds') : null,
      notes.some(n => n.refused) ? h('span', null, h('i.c-refused'), 'refused by the solver: silent') : null,
    ] : null,
    p.real ? [h('span', null, h('span.sw.sim'), 'simulation'), h('span', null, h('span.sw.real'), 'recording')] : null);
  const summary = total ? h('p.pitchsum', null,
    h('b', null, t.in === total ? `Every note in tune (within 10 cents)` : `${t.in} of ${total} notes within 10 cents`),
    t.near ? h('span', null, `${t.near} within 50`) : null, t.out - t.refused ? h('span', null, `${t.out - t.refused} further off`) : null,
    t.refused ? h('span', null, `${t.refused} refused by the solver`) : null,
    h('a', { href: '#/method' }, 'how pitch is measured')) : null;
  const how = h('div.howbox', null,
    h('div', null, h('h4', null, 'How it was made'), h('p', null, p.how || '')),
    h('div', null, h('h4', null, 'What to listen for'), listenFor(p)));
  return h('div.stage', null,
    h('div.detailhead', null, h('h2', null, p.title), h('span.chip.kind', null, KIND_LABEL[p.kind] || p.kind)),
    summary,
    controls(player, { keys: Boolean(p.real) }),
    maps, legend, strip, how);
}

function noteTitle(n) {
  const parts = [n.note];
  if (n.string) parts.push(`on the ${n.string} string${n.semitones ? `, stopped ${n.semitones} semitone${n.semitones > 1 ? 's' : ''} up` : ' (open)'}`);
  if (n.target_hz) parts.push(`asked for ${fmt.num(n.target_hz, 1)} Hz`);
  if (n.refused) parts.push(`not played: the solver refused it (${n.refused})`);
  else if (unpitched(n)) parts.push('no steady pitch: the waveform never settles into one period');
  else if (n.f0_hz) parts.push(`sounds ${fmt.num(n.f0_hz, 1)} Hz`);
  if (n.octave_check) parts.push(`(${n.octave_check})`);
  return parts.join(' · ');
}

function listenFor(p) {
  const notes = p.notes || [];
  const t = pitchTally(p);
  const total = t.in + t.near + t.out;
  const out = [];
  if (p.real) out.push(h('p', null, 'Switch between the simulation and the recording while it plays (keys 1 and 2): the playhead stays put, so you hear the same moment of both.'));
  if (total && t.in === total) out.push(h('p', null, 'The pitch is right throughout; listen instead for the tone: how each note starts, how bright it is, and what the body adds.'));
  const worst = notes.filter(n => pitchClass(n) === 'out');
  if (worst.length) {
    out.push(h('p', null, `${worst.length} note${worst.length > 1 ? 's' : ''} land${worst.length > 1 ? '' : 's'} more than 50 cents from the note asked for (${worst.slice(0, 5).map(n => `${n.note} ${centsText(n)}`).join(', ')}${worst.length > 5 ? ', …' : ''}). On the map, the solid line sits where it actually sounds.`));
  }
  const refused = notes.filter(n => n.refused);
  if (refused.length) {
    out.push(h('p', null, `${refused.map(n => n.note).join(', ')}: silent. The solver refused ${refused.length > 1 ? 'these notes' : 'this note'} part way through ${refused.length > 1 ? 'their renders' : 'its render'} rather than return a wrong answer: `,
      refused.map((n, i) => [i ? '; ' : '', h('code', null, n.refused)]), '.'));
  }
  if ((p.family === 'double-bass' || p.family === 'cello') && t.out) {
    out.push(h('p', null, 'Short notes on the lowest strings do not settle into a steady tone in this model: the body\'s single resonance near 280 Hz takes over. Longer notes do settle.'));
  }
  if (!out.length) out.push(h('p', null, p.native ? 'Every change happens inside one simulation: nothing is spliced.' : 'Each part is its own render; the map shows where one ends and the next begins.'));
  return out;
}
