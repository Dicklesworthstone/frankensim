// Read the native writer's result, not a parallel performance/physics parser.
import {open} from 'node:fs/promises';
import {constants} from 'node:fs';
import {createHash} from 'node:crypto';
import {basename, dirname, extname, join} from 'node:path';
import {inspectNativeWave, validatePressureReceipt, PRESSURE_RECEIPT_FORMAT} from '../src/native-wave.mjs';

export const MAX_PRESSURE_SIDECAR_BYTES = 4 * 1024 * 1024;
const invalid = () => new Error('Native pressure calibration is invalid or does not match the WAV.');

async function windSummary(output) {
  let file;
  try {
    file = await open(join(dirname(output), basename(output, extname(output)) + '.provenance.json'),
      constants.O_RDONLY | (constants.O_NOFOLLOW || 0) | (constants.O_NONBLOCK || 0));
    const info = await file.stat();
    if (!info.isFile() || info.size > MAX_PRESSURE_SIDECAR_BYTES) throw invalid();
    // Read a bounded prefix even if a faulty producer changes its file after stat.
    const bytes = Buffer.alloc(Math.min(info.size + 1, MAX_PRESSURE_SIDECAR_BYTES + 1));
    let length = 0;
    while (length < bytes.length) {
      const read = await file.read(bytes, length, bytes.length - length, length);
      if (!read.bytesRead) break;
      length += read.bytesRead;
    }
    if (length !== info.size) throw invalid();
    const source = JSON.parse(new TextDecoder('utf-8', {fatal:true}).decode(bytes.subarray(0, length)));
    if (source.schema !== 'frankensim-music-render-provenance-v1') throw invalid();
    return {sampleRateHz:source.sample_rate_hz, samples:source.samples,
      fullScalePa:source.full_scale_pa, clippedSamples:source.clipped_samples};
  } catch (error) {
    // Old executables may provide WAV only. Such audio is playable but is NOT
    // eligible for calibrated comparisons; missing calibration never means 1 Pa.
    if (error.code === 'ENOENT') return null;
    throw invalid();
  } finally { await file?.close(); }
}

function stringSummary(stderr) {
  const lines = stderr.split(/\r?\n/).filter(line => line.startsWith('samples='));
  if (!lines.length) return null;
  if (lines.length !== 1) throw invalid();
  const match = /^samples=(\d+); rate_Hz=(\d+); clips=(\d+); full_scale_Pa=(\S+)$/.exec(lines[0]);
  if (!match) throw invalid();
  return {samples:Number(match[1]), sampleRateHz:Number(match[2]),
    clippedSamples:Number(match[3]), fullScalePa:Number(match[4])};
}

/** Missing legacy receipts return null; malformed emitted receipts always refuse. */
export async function nativePressureReceipt(renderer, output, diagnostics, wave) {
  const source = renderer === 'strings' ? stringSummary(diagnostics.stderr) : await windSummary(output);
  if (!source) return null;
  return validatePressureReceipt({format:PRESSURE_RECEIPT_FORMAT, renderer, ...source,
    wavSha256:createHash('sha256').update(wave).digest('hex')}, inspectNativeWave(wave), renderer);
}
