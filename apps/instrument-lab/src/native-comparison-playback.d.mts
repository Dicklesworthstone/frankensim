import type {NativePressureReceipt, NativePressureComparison} from './native-wave.mjs';
export class NativeComparisonPlayback {
  constructor(createContext?: () => AudioContext);
  readonly playing: boolean;
  readonly positionSeconds: number;
  readonly duration: number;
  readonly side: 'a' | 'b';
  setPair(aBytes: Uint8Array, aReceipt: NativePressureReceipt, bBytes: Uint8Array,
    bReceipt: NativePressureReceipt): NativePressureComparison;
  clear(): void;
  play(offsetSeconds?: number): Promise<boolean>;
  select(side: 'a' | 'b'): void;
  setVolume(value: number): void;
  stop(): void;
  close(): Promise<void>;
}
