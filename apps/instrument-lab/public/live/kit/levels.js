// Grid-independent output level for the live kit.
//
// Every page sets one gain per piece, measured on the 5x grid (the live kit's gains against the
// Salamander overheads, Strandbar's against its studio stems). The Worker then picks the grid per
// device and piece (2-10x), and the same stroke's loudness moves with the grid, by piece and not
// monotonically: the crash's standard stroke is 5 dB louder at 4x than at 5x and the closed hi-hat
// 5 dB quieter, while a drum moves by well under 1 dB. engine_levels.json, written next to the wasm
// by tools/kit_levels.mjs (engine-kit/build-wasm.sh runs it), holds each piece's standard stroke per
// grid: its loudness (the largest 400 ms RMS, the momentary-loudness window, so a cymbal's wash
// counts and not only its first milliseconds) and its peak. levelFix() is the factor that gives a
// stroke on `scale` the loudness it has on the reference grid.
//
// The standard stroke stays at or under PEAK_CEILING of full scale whenever its measured peak
// and page gain are valid. This can require an additional cut even at the reference grid.

export const PEAK_CEILING = 0.55;

/** The table's row for `scale`: exact, else interpolated in log level between the neighbours. */
export function levelAt(rows, scale) {
  const keys = Object.keys(rows).map(Number).filter((s) => Number.isFinite(s) && rows[s]).sort((a, b) => a - b);
  if (!keys.length || !Number.isFinite(scale)) return null;
  if (rows[scale]) return rows[scale];
  if (scale <= keys[0]) return rows[keys[0]];
  if (scale >= keys[keys.length - 1]) return rows[keys[keys.length - 1]];
  let i = 1;
  while (keys[i] < scale) i++;
  const lo = rows[keys[i - 1]], hi = rows[keys[i]], t = (scale - keys[i - 1]) / (keys[i] - keys[i - 1]);
  const mix = (a, b) => Math.exp((1 - t) * Math.log(a) + t * Math.log(b));
  return { level: mix(lo.level, hi.level), peak: mix(lo.peak, hi.peak) };
}

/** The gain factor for `piece` on grid `scale` given the page's (reference-grid) gain; 1 without a table. */
export function levelFix(table, piece, scale, gain) {
  const rows = table && table.pieces && table.pieces[piece];
  if (!rows) return 1;
  const ref = rows[table.reference], at = levelAt(rows, scale);
  if (!ref || !at || !Number.isFinite(ref.level) || !(ref.level > 0)
    || !Number.isFinite(at.level) || !(at.level > 0)) return 1;
  const fix = ref.level / at.level;
  const room = Number.isFinite(gain) && gain > 0 && Number.isFinite(at.peak) && at.peak > 0
    ? PEAK_CEILING / (gain * at.peak) : Infinity;
  return Math.min(fix, room);
}
