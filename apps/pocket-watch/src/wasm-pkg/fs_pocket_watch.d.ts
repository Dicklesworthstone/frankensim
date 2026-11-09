/* tslint:disable */
/* eslint-disable */

/**
 * Top-level coupled simulation engine for the FrankenSim Grand Complication Pocket Watch.
 */
export class PocketWatchSim {
    free(): void;
    [Symbol.dispose](): void;
    /**
     * Stokes-Stuart oscillating-cylinder aerodynamic drag on the balance wheel ($\text{s}^{-1}$) (`fs-vkplate::air_loss`).
     */
    balance_aero_drag_s_inv(): number;
    /**
     * Return the 64x8 (512-node) 3D Thermo-Elasto-Viscoplastic + Surface-Tension Case Shell snapshot
     * (`512 * 6 = 3072` `f32`s): `[x_mm, y_mm, z_mm, thickness_mm, temperature_k, melt_fraction]`
     * for each 3D Lagrangian node `(it, iz)` on the Bassine 18k gold case.
     */
    case_shell_snapshot(): Float32Array;
    /**
     * Ultrasonic clean excess flooded oil and restore factory $1.6\text{ }\mu\text{m}$ jewel lubrication.
     */
    clean_excess_oil(): void;
    /**
     * Return the dynamic recoil click lever (*cliquet*) deflection angle ($\text{rad}$).
     */
    click_lever_angle_rad(): number;
    /**
     * Anytime-valid COSC chronometer compliance e-process value $E_t$ (`fs-eproc`).
     */
    cosc_e_value(): number;
    /**
     * Return the keyless works crown wheel rotation angle ($\text{rad}$).
     */
    crown_wheel_angle_rad(): number;
    /**
     * Lagrange cost-optimal GD&T manufacturing tolerance allocation cost (`fs-toleralloc`).
     */
    gdt_allocated_cost(): number;
    /**
     * Return the 18 discretized elastica mainspring coil radii ($\text{mm}$) inside the barrel.
     */
    mainspring_coils_mm(): Float32Array;
    /**
     * ISO 532-1 loudness level in phon at the 3D virtual microphone (`fs-psycho`).
     */
    mic_loudness_phon(): number;
    /**
     * ISO 532-1 Zwicker stationary loudness in sone at the 3D virtual microphone (`fs-psycho`).
     */
    mic_loudness_sone(): number;
    /**
     * DIN 45692 psychoacoustic sharpness in acum at the 3D virtual microphone (`fs-psycho`).
     */
    mic_sharpness_acum(): number;
    /**
     * Sound Pressure Level ($\text{dB SPL re } 20\text{ }\mu\text{Pa}$) at the 3D virtual microphone (`fs-psycho`).
     */
    mic_spl_db(): number;
    /**
     * Multi-crate spectral, topological, regime, UQ, and port-Hamiltonian coupling audit summary
     * exercising `fs-couple`, `fs-regime`, `fs-cheb`, `fs-fft`, `fs-orbit`, `fs-asbuilt`,
     * `fs-rep-sdf`, `fs-rep-mesh`, `fs-topo`, `fs-vfit`, `fs-uq`, and `fs-viz`.
     */
    multiphysics_crate_audit_summary(): Float64Array;
    /**
     * Construct a new, chronometer-calibrated Patek-style Tourbillon pocket watch at $293.15\text{ K}$.
     */
    constructor();
    /**
     * Upper bound of certified interval rate enclosure ($\text{s/day}$) (`fs-ivl`).
     */
    rate_interval_hi_s_day(): number;
    /**
     * Lower bound of certified interval rate enclosure ($\text{s/day}$) (`fs-ivl`).
     */
    rate_interval_lo_s_day(): number;
    /**
     * Exact forward-mode AD thermoelastic rate sensitivity $\partial r / \partial T$ ($\text{s/day/K}$) (`fs-ad`).
     */
    rate_sensitivity_s_day_per_k(): number;
    /**
     * Render `frames` mono audio samples at $48\text{ kHz}$ from the 3D virtual microphone
     * and return them as a `Vec<f32>` (converted to `Float32Array` across `wasm-bindgen`).
     */
    render_audio_block(frames: number): Float32Array;
    /**
     * Reset the entire watch and all materials back to pristine $20.0^\circ\text{C}$ chronometer state.
     */
    reset(): void;
    /**
     * Configure the virtual oxy-butane blowtorch ($T_{\text{flame}} = 2250\text{ K}$) in 3D watch coordinates ($\text{mm}$).
     */
    set_blowtorch(power: number, x_mm: number, y_mm: number, z_mm: number): void;
    /**
     * Toggle between traditional High-Carbon Blued Steel (`false`) and Compensated Invar/Fe-36Ni (`true`) hairspring.
     */
    set_compensated_hairspring(compensated: boolean): void;
    /**
     * Configure the cryogenic liquid nitrogen ($\text{LN}_2$ at $77.36\text{ K}$) pour stream in 3D watch coordinates ($\text{mm}$).
     */
    set_liquid_nitrogen(pour_rate: number, x_mm: number, y_mm: number, z_mm: number): void;
    /**
     * Configure the external sintered $\text{NdFeB N42}$ super-magnet ($B_r = 1.32\text{ T}$) in 3D watch coordinates ($\text{mm}$).
     */
    set_magnet(enabled: boolean, strength: number, x_mm: number, y_mm: number, z_mm: number): void;
    /**
     * Position the 3D studio measurement microphone in watch coordinates ($\text{mm}$) and set gain.
     */
    set_microphone(x_mm: number, y_mm: number, z_mm: number, gain: number): void;
    /**
     * Configure the watch oil dropper / flood stream (`pennzane-shf-x-2000-bearing-oil`) in 3D watch coordinates ($\text{mm}$).
     */
    set_oil_pour(pour_rate: number, x_mm: number, y_mm: number, z_mm: number): void;
    /**
     * Return the current 176-element `f64` state snapshot for 3D rendering and telemetry.
     */
    snapshot(): Float64Array;
    /**
     * Return the 24x24 (576-cell) 2D Finite-Volume Spatial PDE grid snapshot (`576 * 7 = 4032` `f32`s):
     * `[temperature_k, melt_fraction, case_slump_frac, liquid_gold_mm, frozen_gold_mm, frost_um, oil_um]`
     * for each cell `(ix, iy)` in row-major order.
     */
    spatial_grid_snapshot(): Float32Array;
    /**
     * Deterministic BLAKE3 digest of the current 176-float state snapshot (`fs-blake3`).
     */
    state_digest_hex(): string;
    /**
     * Human-readable status banner for the current horological & multiphysics state.
     */
    status_label(): string;
    /**
     * Advance the coupled multiphysics and horological simulation by `dt_s` seconds.
     */
    step(dt_s: number): void;
    /**
     * Fire the alternating-field degaussing coil to wipe out residual ferromagnetic remanence $M_{\text{rem}}$.
     */
    trigger_degauss(): void;
    /**
     * Wind the mainspring crown by `delta_turns` (where `8.5` turns = full 42-hour wind).
     */
    wind_crown(delta_turns: number): void;
}

export type InitInput = RequestInfo | URL | Response | BufferSource | WebAssembly.Module;

export interface InitOutput {
    readonly memory: WebAssembly.Memory;
    readonly __wbg_pocketwatchsim_free: (a: number, b: number) => void;
    readonly pocketwatchsim_balance_aero_drag_s_inv: (a: number) => number;
    readonly pocketwatchsim_case_shell_snapshot: (a: number) => [number, number];
    readonly pocketwatchsim_clean_excess_oil: (a: number) => void;
    readonly pocketwatchsim_click_lever_angle_rad: (a: number) => number;
    readonly pocketwatchsim_cosc_e_value: (a: number) => number;
    readonly pocketwatchsim_crown_wheel_angle_rad: (a: number) => number;
    readonly pocketwatchsim_gdt_allocated_cost: (a: number) => number;
    readonly pocketwatchsim_mainspring_coils_mm: (a: number) => [number, number];
    readonly pocketwatchsim_mic_loudness_phon: (a: number) => number;
    readonly pocketwatchsim_mic_loudness_sone: (a: number) => number;
    readonly pocketwatchsim_mic_sharpness_acum: (a: number) => number;
    readonly pocketwatchsim_mic_spl_db: (a: number) => number;
    readonly pocketwatchsim_multiphysics_crate_audit_summary: (a: number) => [number, number];
    readonly pocketwatchsim_new: () => number;
    readonly pocketwatchsim_rate_interval_hi_s_day: (a: number) => number;
    readonly pocketwatchsim_rate_interval_lo_s_day: (a: number) => number;
    readonly pocketwatchsim_rate_sensitivity_s_day_per_k: (a: number) => number;
    readonly pocketwatchsim_render_audio_block: (a: number, b: number) => [number, number];
    readonly pocketwatchsim_reset: (a: number) => void;
    readonly pocketwatchsim_set_blowtorch: (a: number, b: number, c: number, d: number, e: number) => void;
    readonly pocketwatchsim_set_compensated_hairspring: (a: number, b: number) => void;
    readonly pocketwatchsim_set_liquid_nitrogen: (a: number, b: number, c: number, d: number, e: number) => void;
    readonly pocketwatchsim_set_magnet: (a: number, b: number, c: number, d: number, e: number, f: number) => void;
    readonly pocketwatchsim_set_microphone: (a: number, b: number, c: number, d: number, e: number) => void;
    readonly pocketwatchsim_set_oil_pour: (a: number, b: number, c: number, d: number, e: number) => void;
    readonly pocketwatchsim_snapshot: (a: number) => [number, number];
    readonly pocketwatchsim_spatial_grid_snapshot: (a: number) => [number, number];
    readonly pocketwatchsim_state_digest_hex: (a: number) => [number, number];
    readonly pocketwatchsim_status_label: (a: number) => [number, number];
    readonly pocketwatchsim_step: (a: number, b: number) => void;
    readonly pocketwatchsim_trigger_degauss: (a: number) => void;
    readonly pocketwatchsim_wind_crown: (a: number, b: number) => void;
    readonly __wbindgen_externrefs: WebAssembly.Table;
    readonly __wbindgen_free: (a: number, b: number, c: number) => void;
    readonly __wbindgen_start: () => void;
}

export type SyncInitInput = BufferSource | WebAssembly.Module;

/**
 * Instantiates the given `module`, which can either be bytes or
 * a precompiled `WebAssembly.Module`.
 *
 * @param {{ module: SyncInitInput }} module - Passing `SyncInitInput` directly is deprecated.
 *
 * @returns {InitOutput}
 */
export function initSync(module: { module: SyncInitInput } | SyncInitInput): InitOutput;

/**
 * If `module_or_path` is {RequestInfo} or {URL}, makes a request and
 * for everything else, calls `WebAssembly.instantiate` directly.
 *
 * @param {{ module_or_path: InitInput | Promise<InitInput> }} module_or_path - Passing `InitInput` directly is deprecated.
 *
 * @returns {Promise<InitOutput>}
 */
export default function __wbg_init (module_or_path?: { module_or_path: InitInput | Promise<InitInput> } | InitInput | Promise<InitInput>): Promise<InitOutput>;
