import type {NativeWindAttachment, NativeWindGasSources, NativeWindOptions, NativeStringProfiles, NativeCoatingStrainGuard, NativeFiniteGasStringSources, NativeRenderer} from './native-wind-client.mjs';

export const MAX_SOURCE_BYTES: number;
export const MAX_ATTACHMENTS: number;
export const MAX_ATTACHMENT_FILE_BYTES: number;
export const MAX_ATTACHMENT_BYTES: number;
export const MAX_REQUEST_BYTES: number;
export const MAX_HISTORY_BYTES: number;
export const FINITE_GAS_STRING_RESULT_FORMAT: string;
export const SUPPORTED_SCHEMAS: readonly string[];
export const STRING_SCHEMAS: readonly string[];
export const RESOLVED_SCHEMAS: readonly string[];
export const FINITE_GAS_STRING_SCHEMAS: readonly string[];
export class NativeWindError extends Error { status: number; constructor(status: number, message: string); }
export interface NativeExperiment {
  source: Uint8Array;
  attachments: NativeWindAttachment[];
  renderer?: Exclude<NativeRenderer, 'wind'>;
  decimate: boolean;
  gasSources?: NativeWindGasSources;
  stringProfiles?: NativeStringProfiles;
  coatingStrainGuard?: NativeCoatingStrainGuard;
  finiteGasString?: NativeFiniteGasStringSources;
}
export function editableAttachmentText(bytes: Uint8Array): string;
export function editedAttachmentBytes(original: Uint8Array, text: string): Uint8Array;
export function validateSource(bytes: Uint8Array, renderer?: NativeRenderer): Uint8Array;
export function validateInput(bytes: Uint8Array, bundled?: boolean): NativeExperiment;
export function serializeNativeExperiment(source: string, attachments: readonly NativeWindAttachment[], options?: NativeWindOptions): string;
