/** Transport contract of tools/serve-native-wind.mjs; no instrument or material defaults. */
export interface NativeWindStatus {
  available: boolean;
  busy: boolean;
  token?: string;
  maximumSourceBytes: number;
  maximumOutputBytes: number;
  timeoutMs: number;
}
export type NativeWindRequest = typeof globalThis.fetch;
export function nativeWindStatus(signal?: AbortSignal, request?: NativeWindRequest): Promise<NativeWindStatus>;
export function renderNativeWind(source: string, status: NativeWindStatus, signal?: AbortSignal, request?: NativeWindRequest): Promise<Blob>;
