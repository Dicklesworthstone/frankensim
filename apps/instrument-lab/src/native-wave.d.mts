export const NATIVE_WAVE_LIMIT: number;
export const PRESSURE_RECEIPT_HEADER: string;
export const PRESSURE_RECEIPT_FORMAT: string;
export interface NativeWaveInfo {
  readonly sampleRateHz: number;
  readonly samples: number;
  readonly dataOffset: number;
  readonly dataBytes: number;
}
export interface NativePressureReceipt {
  readonly format: string;
  readonly renderer: 'wind' | 'strings';
  readonly sampleRateHz: number;
  readonly samples: number;
  readonly fullScalePa: number;
  readonly clippedSamples: number;
  readonly wavSha256: string;
}
export interface NativePressureComparison {
  readonly samples: number;
  readonly sampleRateHz: number;
  readonly durationSeconds: number;
  readonly truncated: boolean;
  readonly fullScalePa: number;
  readonly rmsA: number;
  readonly rmsB: number;
  readonly rmsDifference: number;
  readonly peakA: number;
  readonly peakB: number;
  readonly levelDifferenceDb: number | null;
  readonly identical: boolean;
}
export function inspectNativeWave(bytes: Uint8Array, limit?: number): NativeWaveInfo;
export function validatePressureReceipt(value: unknown, info: NativeWaveInfo, renderer?: string): NativePressureReceipt;
export function verifyPressureReceipt(value: unknown, bytes: Uint8Array, renderer?: string): Promise<NativePressureReceipt>;
export function nativePcmSamples(bytes: Uint8Array): Float32Array;
export function compareNativePressure(aBytes: Uint8Array, aReceipt: NativePressureReceipt,
  bBytes: Uint8Array, bReceipt: NativePressureReceipt): NativePressureComparison;
