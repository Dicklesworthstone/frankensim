import test from 'node:test';
import assert from 'node:assert/strict';
import { levelFix, PEAK_CEILING } from '../public/live/kit/levels.js';

const table = rows => ({ reference: 5, pieces: { hihat: rows } });

test('the reference grid still cuts a measured standard stroke above the peak ceiling', () => {
  const measured = table({ 5: { level: 1, peak: 10 } });
  const gain = .1, fix = levelFix(measured, 'hihat', 5, gain);
  assert.equal(fix, .55);
  assert.equal(gain * fix * measured.pieces.hihat[5].peak, PEAK_CEILING);
});

test('a loudness cut also respects the actual peak instead of stopping at unity gain', () => {
  const measured = table({ 4: { level: 2, peak: 30 }, 5: { level: 1, peak: 10 } });
  const gain = .1, fix = levelFix(measured, 'hihat', 4, gain);
  assert.equal(fix, PEAK_CEILING / 3);
  assert(fix < .5, 'the relative loudness cut alone still exceeds the ceiling');
  assert(Math.abs(gain * fix * measured.pieces.hihat[4].peak - PEAK_CEILING) <= Number.EPSILON);
});

test('safe boosts and unavailable peak measurements retain their declared loudness behavior', () => {
  const measured = table({ 4: { level: .5, peak: 1 }, 5: { level: 1, peak: 10 } });
  assert.equal(levelFix(measured, 'hihat', 4, .1), 2);
  for (const peak of [undefined, 0, -1, NaN, Infinity]) {
    const noPeak = table({ 4: { level: .5, peak }, 5: { level: 1, peak: 10 } });
    assert.equal(levelFix(noPeak, 'hihat', 4, .1), 2);
  }
  for (const gain of [undefined, 0, -1, NaN, Infinity]) {
    assert.equal(levelFix(measured, 'hihat', 4, gain), 2);
  }
  for (const level of [0, -1, NaN, Infinity]) {
    assert.equal(levelFix(table({ 5: { level, peak: 10 } }), 'hihat', 5, .1), 1);
  }
  assert.equal(levelFix(null, 'hihat', 5, .1), 1);
  assert.equal(levelFix(measured, 'missing', 5, .1), 1);
});
