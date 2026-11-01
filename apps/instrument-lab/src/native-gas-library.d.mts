import type {NativeWindAttachment, NativeWindGasComponent, NativeMaterialManifest} from './native-wind-client.mjs';

export const NATIVE_GAS_LIBRARY: readonly {
  id: string; label: string; species: string;
  speciesSource: string; nasa9Source: string; transportSource: string;
}[];
export interface NativeGasSourceAddition {
  attachments: NativeWindAttachment[];
  materialManifests: NativeMaterialManifest[];
  component: NativeWindGasComponent;
}
export function loadNativeGasSource(id: string, readSource: (path: string) => Promise<string>): Promise<NativeGasSourceAddition>;
export function mergeNativeGasSource(attachments: readonly NativeWindAttachment[],
  materialManifests: readonly NativeMaterialManifest[] | undefined, added: NativeGasSourceAddition,
  maximumAttachments?: number): {attachments: NativeWindAttachment[]; materialManifests: NativeMaterialManifest[]};
