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
}
export type NativeWindRequest = typeof globalThis.fetch;
export function nativeWindStatus(signal?: AbortSignal, request?: NativeWindRequest): Promise<NativeWindStatus>;
export function renderNativeWind(source: string, status: NativeWindStatus, signal?: AbortSignal, request?: NativeWindRequest): Promise<Blob>;

export interface NativeWindAttachment { path: string; bytes: Uint8Array; }
export interface NativeWindOptions { decimate?: boolean; }
export function renderNativeWindBundle(source: string, attachments: readonly NativeWindAttachment[], status: NativeWindStatus, options?: NativeWindOptions, signal?: AbortSignal, request?: NativeWindRequest): Promise<Blob>;
