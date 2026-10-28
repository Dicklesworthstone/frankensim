// One reader for the native service and calibrated browser comparisons.
// This observes the native PCM result; it never normalizes or synthesizes sound.
export const NATIVE_WAVE_LIMIT = 16 * 1024 * 1024;
export const PRESSURE_RECEIPT_HEADER = 'X-Native-Pressure-Receipt';
export const PRESSURE_RECEIPT_FORMAT = 'frankensim-native-pressure-v1';

export function inspectNativeWave(bytes, limit = NATIVE_WAVE_LIMIT) {
  const bad = () => { throw new Error('Native renderer did not produce a complete 48 kHz mono PCM16 WAV.'); };
  if (!(bytes instanceof Uint8Array) || bytes.byteLength < 44 || bytes.byteLength > limit) bad();
  const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  const tag = at => view.getUint32(at, false);
  if (tag(0) !== 0x52494646 || tag(8) !== 0x57415645 || view.getUint32(4, true) + 8 !== bytes.byteLength) bad();
  let position = 12, format = false, dataOffset = -1, dataBytes = 0;
  while (position + 8 <= bytes.byteLength) {
    const name = tag(position), size = view.getUint32(position + 4, true), start = position + 8;
    const end = start + size;
    if (end > bytes.byteLength) bad();
    if (name === 0x666d7420) {
      if (format || size < 16 || view.getUint16(start, true) !== 1
          || view.getUint16(start + 2, true) !== 1 || view.getUint32(start + 4, true) !== 48000
          || view.getUint32(start + 8, true) !== 96000 || view.getUint16(start + 12, true) !== 2
          || view.getUint16(start + 14, true) !== 16) bad();
      format = true;
    } else if (name === 0x64617461) {
      if (dataOffset !== -1 || size === 0 || size % 2 !== 0) bad();
      dataOffset = start; dataBytes = size;
    }
    position = end + (size & 1);
  }
  if (position !== bytes.byteLength || !format || dataOffset === -1) bad();
  return Object.freeze({sampleRateHz:48000, samples:dataBytes / 2, dataOffset, dataBytes});
}

/** Calibration is emitted by the native writer, never guessed from input text. */
export function validatePressureReceipt(value, info, renderer) {
  if (!value || value.format !== PRESSURE_RECEIPT_FORMAT
      || !['wind', 'strings', 'resolved'].includes(value.renderer)
      || (renderer !== undefined && value.renderer !== renderer)
      || value.sampleRateHz !== info.sampleRateHz || value.samples !== info.samples
      || !Number.isFinite(value.fullScalePa) || value.fullScalePa <= 0
      || !Number.isSafeInteger(value.clippedSamples) || value.clippedSamples < 0
      || value.clippedSamples > info.samples
      || typeof value.wavSha256 !== 'string' || !/^[0-9a-f]{64}$/.test(value.wavSha256)) {
    throw new Error('Native pressure calibration does not match the returned WAV.');
  }
  return Object.freeze({format:PRESSURE_RECEIPT_FORMAT, renderer:value.renderer,
    sampleRateHz:info.sampleRateHz, samples:info.samples, fullScalePa:value.fullScalePa,
    clippedSamples:value.clippedSamples, wavSha256:value.wavSha256});
}

export async function verifyPressureReceipt(value, bytes, renderer) {
  const receipt = validatePressureReceipt(value, inspectNativeWave(bytes), renderer);
  const digest = await globalThis.crypto.subtle.digest('SHA-256', bytes);
  const actual = Array.from(new Uint8Array(digest), byte => byte.toString(16).padStart(2, '0')).join('');
  if (actual !== receipt.wavSha256) throw new Error('Native pressure calibration belongs to different WAV bytes.');
  return receipt;
}

/** Decode original integers for playback, without browser decodeAudioData resampling. */
export function nativePcmSamples(bytes) {
  const info = inspectNativeWave(bytes);
  const view = new DataView(bytes.buffer, bytes.byteOffset + info.dataOffset, info.dataBytes);
  const samples = new Float32Array(info.samples);
  for (let i = 0; i < samples.length; i++) samples[i] = view.getInt16(2 * i, true) / 32768;
  return samples;
}

/** Pa reconstructed from the native encoder's i16::MAX scale, after quantization. */
export function compareNativePressure(aBytes, aReceipt, bBytes, bReceipt) {
  const aInfo = inspectNativeWave(aBytes), bInfo = inspectNativeWave(bBytes);
  const a = validatePressureReceipt(aReceipt, aInfo), b = validatePressureReceipt(bReceipt, bInfo);
  if (a.fullScalePa !== b.fullScalePa) {
    throw new Error('Use the same native full-scale pressure for A and B, then render again. No automatic level matching is applied.');
  }
  if (a.clippedSamples || b.clippedSamples) {
    throw new Error('A or B clipped in the native encoder. Render both at the same larger full-scale pressure; playback gain cannot restore lost pressure.');
  }
  const samples = Math.min(a.samples, b.samples);
  const av = new DataView(aBytes.buffer, aBytes.byteOffset + aInfo.dataOffset, aInfo.dataBytes);
  const bv = new DataView(bBytes.buffer, bBytes.byteOffset + bInfo.dataOffset, bInfo.dataBytes);
  let aSquared = 0, bSquared = 0, differenceSquared = 0, aPeak = 0, bPeak = 0;
  for (let i = 0; i < samples; i++) {
    const x = av.getInt16(i * 2, true), y = bv.getInt16(i * 2, true);
    aSquared += x * x; bSquared += y * y; differenceSquared += (y - x) * (y - x);
    aPeak = Math.max(aPeak, Math.abs(x)); bPeak = Math.max(bPeak, Math.abs(y));
  }
  const toPa = value => value / 32767 * a.fullScalePa;
  const rms = sum => toPa(Math.sqrt(sum / samples));
  const rmsA = rms(aSquared), rmsB = rms(bSquared), rmsDifference = rms(differenceSquared);
  const peakA = toPa(aPeak), peakB = toPa(bPeak);
  if (![rmsA, rmsB, rmsDifference, peakA, peakB].every(Number.isFinite)) {
    throw new Error('PCM pressure comparison is not representable at this native scale.');
  }
  return Object.freeze({samples, sampleRateHz:a.sampleRateHz, durationSeconds:samples / a.sampleRateHz,
    truncated:a.samples !== b.samples, fullScalePa:a.fullScalePa,
    rmsA, rmsB, rmsDifference, peakA, peakB,
    // Null, rather than Infinity/NaN, keeps silent references honest in the UI.
    levelDifferenceDb:aSquared > 0 && bSquared > 0 ? 10 * Math.log10(bSquared / aSquared) : null,
    identical:differenceSquared === 0});
}
