/** Transport contract of tools/serve-native-wind.mjs; no instrument or material defaults. */
import type {NativePressureReceipt} from './native-wave.mjs';
export interface NativeWindStatus {
  available: boolean;
  renderers?: { wind: boolean; strings: boolean; resolved?: boolean; 'finite-gas-string'?: boolean };
  busy: boolean;
  token?: string;
  maximumSourceBytes: number;
  maximumOutputBytes: number;
  maximumHistoryBytes?: number;
  timeoutMs: number;
  maximumAttachments?: number;
  maximumAttachmentFileBytes?: number;
  maximumAttachmentBytes?: number;
  maximumRequestBytes?: number;
  supportsGasSources?: boolean;
  supportsPressureReceipts?: boolean;
  supportsCoatingStrainGuard?: boolean;
  supportsMaterialManifests?: boolean;
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
export interface NativeCoatingStrainGuard { maximumBytes: number; maximumWork: number; }
export interface NativeMaterialManifest {
  packPath: string;
  manifestPath: string;
  materialState?: {chemistry: string; phase: string; process: string; revision: number};
}
export interface NativeFiniteGasStringSources { solidPath: string; gasRecordPath: string; }
export type NativeRenderer = 'wind' | 'strings' | 'resolved' | 'finite-gas-string';
export interface NativeWindOptions {
  renderer?: NativeRenderer;
  decimate?: boolean;
  gasSources?: NativeWindGasSources;
  stringProfiles?: NativeStringProfiles;
  coatingStrainGuard?: NativeCoatingStrainGuard;
  finiteGasString?: NativeFiniteGasStringSources;
  materialManifests?: readonly NativeMaterialManifest[];
}
export function renderNativeWindBundle(source: string, attachments: readonly NativeWindAttachment[], status: NativeWindStatus, options?: NativeWindOptions, signal?: AbortSignal, request?: NativeWindRequest): Promise<Blob>;
export interface NativePerformanceResult { wave: Blob; receipt: NativePressureReceipt | null; history?: Blob; }
export function renderNativePerformance(source: string, attachments: readonly NativeWindAttachment[], status: NativeWindStatus, options?: NativeWindOptions, signal?: AbortSignal, request?: NativeWindRequest): Promise<NativePerformanceResult>;
