export interface NativeGasValues {
  kind: 'dry-air' | 'moist-air' | 'mixture' | 'thermochemical';
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
  supportedKinds: NativeGasValues['kind'][];
  label: string;
  values: NativeGasValues;
  pressureSpan: [number, number] | null;
  /** Station metadata is inherited from one complete native dg_bore field. */
  field?: {
    columns: string[];
    headerLine: number;
    position: string;
    station: number;
    stations: number;
  };
}
export function nativeGasRegions(source:string): NativeGasRegion[];
export function editNativeGasRegion(source:string, region:NativeGasRegion, values:NativeGasValues, commonPressure?:boolean): string;
