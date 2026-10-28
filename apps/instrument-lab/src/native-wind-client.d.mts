/** Transport contract of tools/serve-native-wind.mjs; no instrument or material defaults. */
import type {NativePressureReceipt} from './native-wave.mjs';
export interface NativeWindStatus {
  available: boolean;
  renderers?: { wind: boolean; strings: boolean; resolved?: boolean };
  busy: boolean;
  token?: string;
  maximumSourceBytes: number;
  maximumOutputBytes: number;
  timeoutMs: number;
  maximumAttachments?: number;
  maximumAttachmentFileBytes?: number;
  maximumAttachmentBytes?: number;
  maximumRequestBytes?: number;
  supportsGasSources?: boolean;
  supportsPressureReceipts?: boolean;
}
export type NativeWindRequest = typeof globalThis.fetch;
export function nativeWindStatus(signal?: AbortSignal, request?: NativeWindRequest): Promise<NativeWindStatus>;
export function renderNativeWind(source: string, status: NativeWindStatus, signal?: AbortSignal, request?: NativeWindRequest): Promise<Blob>;

export interface NativeWindAttachment { path: string; bytes: Uint8Array; }
export interface NativeWindGasComponent {
  speciesPath: string;
  nasa9Path: string;
  transportPath: string;
}
export interface NativeWindGasSources {
  components: readonly NativeWindGasComponent[];
  waterSpecies: string | null;
  phase: "frozen" | "liquid-water";
}
export interface NativeStringProfiles { bodyPath?: string; stopPath?: string; contactThermalPath?: string; }
export interface NativeWindOptions {
  renderer?: "wind" | "strings" | "resolved";
  decimate?: boolean;
  gasSources?: NativeWindGasSources;
  stringProfiles?: NativeStringProfiles;
}
export function renderNativeWindBundle(source: string, attachments: readonly NativeWindAttachment[], status: NativeWindStatus, options?: NativeWindOptions, signal?: AbortSignal, request?: NativeWindRequest): Promise<Blob>;
export interface NativePerformanceResult { wave: Blob; receipt: NativePressureReceipt | null; }
export function renderNativePerformance(source: string, attachments: readonly NativeWindAttachment[], status: NativeWindStatus, options?: NativeWindOptions, signal?: AbortSignal, request?: NativeWindRequest): Promise<NativePerformanceResult>;
