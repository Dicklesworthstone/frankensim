// The wasm pkg is a DERIVED artifact (npm run wasm builds it from
// crates/fs-flyer-wasm via wasm-pack --target web; gitignored). This
// declaration is the typed boundary the worker compiles against.
declare module "*/wasm-pkg/fs_flyer_wasm.js" {
  export default function init(input?: unknown): Promise<unknown>;
  export function flyer_engine_init(
    seed: bigint,
    rho_kg_m3: number,
    headwind_mps: number,
    mode: number,
    member: number,
    rail_length_m: number,
    max_ticks: bigint,
    assist: boolean,
    catapult: boolean,
  ): string;
  export function flyer_engine_step(
    has_input: boolean,
    lever_force_n: number,
    warp_cmd_rad: number,
  ): string;
  export function flyer_engine_digest(): string;
}
