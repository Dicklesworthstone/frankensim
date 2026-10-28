export interface NativeGasValues {
  kind: 'moist-air' | 'mixture' | 'thermochemical';
  temperature: string;
  pressure: string;
  humidity?: string;
  basis?: 'mole' | 'mass';
  components: {species:string; fraction:string}[];
}
export interface NativeGasRegion {
  source: string;
  start: number;
  end: number;
  line: number;
  key: string;
  inherited: boolean;
  label: string;
  values: NativeGasValues;
  pressureSpan: [number, number] | null;
}
export function nativeGasRegions(source:string): NativeGasRegion[];
export function editNativeGasRegion(source:string, region:NativeGasRegion, values:NativeGasValues, commonPressure?:boolean): string;
