import type {NativeWindAttachment, NativeWindGasSources, NativeWindOptions, NativeStringProfiles} from './native-wind-client.mjs';

export const MAX_SOURCE_BYTES: number;
export const MAX_ATTACHMENTS: number;
export const MAX_ATTACHMENT_FILE_BYTES: number;
export const MAX_ATTACHMENT_BYTES: number;
export const MAX_REQUEST_BYTES: number;
export const SUPPORTED_SCHEMAS: readonly string[];
export const STRING_SCHEMAS: readonly string[];
export class NativeWindError extends Error { status: number; constructor(status: number, message: string); }
export interface NativeExperiment {
  source: Uint8Array;
  attachments: NativeWindAttachment[];
  renderer?: 'strings';
  decimate: boolean;
  gasSources?: NativeWindGasSources;
  stringProfiles?: NativeStringProfiles;
}
export function editableAttachmentText(bytes: Uint8Array): string;
export function editedAttachmentBytes(original: Uint8Array, text: string): Uint8Array;
export function validateSource(bytes: Uint8Array, renderer?: 'wind' | 'strings'): Uint8Array;
export function validateInput(bytes: Uint8Array, bundled?: boolean): NativeExperiment;
export function serializeNativeExperiment(source: string, attachments: readonly NativeWindAttachment[], options?: NativeWindOptions): string;
