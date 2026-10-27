/** Transport contract of tools/serve-native-wind.mjs; no instrument or material defaults. */
export interface NativeWindStatus {
  available: boolean;
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
export interface NativeWindOptions { decimate?: boolean; gasSources?: NativeWindGasSources; }
export function renderNativeWindBundle(source: string, attachments: readonly NativeWindAttachment[], status: NativeWindStatus, options?: NativeWindOptions, signal?: AbortSignal, request?: NativeWindRequest): Promise<Blob>;
