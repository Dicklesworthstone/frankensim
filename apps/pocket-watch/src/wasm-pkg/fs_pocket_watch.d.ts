/* tslint:disable */
/* eslint-disable */

/**
 * Top-level coupled simulation of the pocket watch.
 */
export class PocketWatchSim {
    free(): void;
    [Symbol.dispose](): void;
    /**
     * Air damping of the balance as a decay rate [1/s].
     */
    balance_aero_drag_s_inv(): number;
    /**
     * The case-band shell, 64 × 8 nodes × `[x, y, z, thickness, T, melt]`
     * (mm, K, 0–1).
     */
    case_shell_snapshot(): Float32Array;
    /**
     * Ultrasonic clean: strip flooded oil, frost and LN₂ from the movement
     * and re-oil the jewels.
     */
    clean_excess_oil(): void;
    /**
     * Click (ratchet pawl) angle [rad].
     */
    click_lever_angle_rad(): number;
    /**
     * Anytime-valid e-value that the watch is outside the COSC −4/+6 s/day
     * band (`fs-eproc` betting e-process).
     */
    cosc_e_value(): number;
    /**
     * Crown-wheel angle [rad].
     */
    crown_wheel_angle_rad(): number;
    /**
     * The 16-value escapement state: θ, θ̇, φ, φ̇, ψ, ψ̇ (balance, lever and
     * escape wheel in the cage frame), the active-contact mask (bit 0 entry
     * stone, 1 exit stone, 2 pin/fork, 3 banking, 4 guard), the last impulse
     * energy [µJ], escapement efficiency, measured lift [deg], balance Q,
     * escape-wheel drive [µN m], barrel torque [mN m], power reserve [h],
     * balance pivot friction [nN m] and the air share of the losses.
     */
    escapement_state(): Float64Array;
    /**
     * Half the side of the plan grid [mm].
     */
    grid_half_span_mm(): number;
    /**
     * Cells per side of the plan grid behind `spatial_grid_snapshot`.
     */
    grid_n(): number;
    /**
     * Mainspring coil radii inside the barrel [mm], innermost first.
     */
    mainspring_coils_mm(): Float32Array;
    /**
     * ISO 532-1 loudness level at the microphone [phon] (`fs-psycho`).
     */
    mic_loudness_phon(): number;
    /**
     * ISO 532-1 loudness at the microphone [sone] (`fs-psycho`).
     */
    mic_loudness_sone(): number;
    /**
     * DIN 45692 sharpness at the microphone [acum] (`fs-psycho`).
     */
    mic_sharpness_acum(): number;
    /**
     * Sound pressure level at the microphone [dB re 20 µPa].
     */
    mic_spl_db(): number;
    /**
     * A regulated, wound tourbillon pocket watch at 20 °C.
     */
    constructor();
    /**
     * Upper end of the ±2σ interval on the measured rate [s/day].
     */
    rate_interval_hi_s_day(): number;
    /**
     * Lower end of the ±2σ interval on the measured rate [s/day].
     */
    rate_interval_lo_s_day(): number;
    /**
     * Rate change per kelvin of the balance and hairspring [s/day/K], from
     * the alloy's thermoelastic coefficient and the balance's expansion.
     */
    rate_sensitivity_s_day_per_k(): number;
    /**
     * Render `frames` (1–8192) mono samples at 48 kHz from the microphone.
     */
    render_audio_block(frames: number): Float32Array;
    /**
     * Back to the pristine watch at 20 °C (the hairspring choice is kept).
     */
    reset(): void;
    /**
     * Compressed-air nozzle: power 0–1 and position [mm].
     */
    set_air_jet(power: number, x_mm: number, y_mm: number, z_mm: number): void;
    /**
     * Bell-jar pressure [atm], 0 (vacuum) to 5.
     */
    set_ambient_pressure_atm(pressure_atm: number): void;
    /**
     * Air–propane torch: power 0–1 and nozzle position in watch coordinates [mm].
     */
    set_blowtorch(power: number, x_mm: number, y_mm: number, z_mm: number): void;
    /**
     * Hairspring alloy: carbon steel (`false`) or a Nivarox-type
     * thermocompensating alloy (`true`). Each is regulated at the factory.
     */
    set_compensated_hairspring(compensated: boolean): void;
    /**
     * Gravity / centrifuge acceleration in g, watch frame (x to 3 o'clock,
     * y to 12 o'clock, z out of the dial).
     */
    set_gravity_vector(gx: number, gy: number, gz: number): void;
    /**
     * Liquid-nitrogen pour: rate 0–1 and spout position [mm].
     */
    set_liquid_nitrogen(pour_rate: number, x_mm: number, y_mm: number, z_mm: number): void;
    /**
     * N42 NdFeB magnet: enabled, strength (0–2, scales the remanence) and
     * position [mm].
     */
    set_magnet(enabled: boolean, strength: number, x_mm: number, y_mm: number, z_mm: number): void;
    /**
     * Measurement microphone position [mm] and gain (0–4).
     */
    set_microphone(x_mm: number, y_mm: number, z_mm: number, gain: number): void;
    /**
     * Watch-oil dropper: rate 0–1 and nozzle position [mm].
     */
    set_oil_pour(pour_rate: number, x_mm: number, y_mm: number, z_mm: number): void;
    /**
     * The 176-value state snapshot (32 header values, then 12 per node).
     */
    snapshot(): Float64Array;
    /**
     * The plan grid, `grid_n()²` cells row-major from (−span, −span), each
     * `[T K, melt 0–1, case slump 0–1, liquid gold mm, frozen gold mm,
     * frost µm, oil µm]`.
     */
    spatial_grid_snapshot(): Float32Array;
    /**
     * BLAKE3 digest of the snapshot (`fs-blake3`), for replay checks.
     */
    state_digest_hex(): string;
    /**
     * Status banner: which mechanism dominates the movement's behaviour.
     */
    status_label(): string;
    /**
     * Advance the coupled simulation by `dt_s` seconds (clamped to 0.1 ms –
     * 100 ms), in slices of at most 10 ms: the fields first (placing the
     * lever and cage for the magnet), then the movement in the environment
     * the fields leave, then the sound sources it produced.
     */
    step(dt_s: number): void;
    /**
     * Lock drops since the last call, as `[lab_time_s, side, energy_nJ]`
     * triples (side 0 = entry, 1 = exit), for the timegrapher.
     */
    timegrapher_drain(): Float64Array;
    /**
     * Fire the 50 Hz degaussing coil (a 3 s geometrically decaying alternating field).
     */
    trigger_degauss(): void;
    /**
     * A drop or knock: a half-sine acceleration pulse of `g_peak` (10–5000 g)
     * applied to the movement and the liquids on it.
     */
    trigger_mechanical_shock(g_peak: number): void;
    /**
     * Turn the crown by `delta_turns` (the crown-to-barrel gearing and the
     * click are simulated; the bridle slips past full wind).
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
    readonly pocketwatchsim_escapement_state: (a: number) => [number, number];
    readonly pocketwatchsim_grid_half_span_mm: (a: number) => number;
    readonly pocketwatchsim_grid_n: (a: number) => number;
    readonly pocketwatchsim_mainspring_coils_mm: (a: number) => [number, number];
    readonly pocketwatchsim_mic_loudness_phon: (a: number) => number;
    readonly pocketwatchsim_mic_loudness_sone: (a: number) => number;
    readonly pocketwatchsim_mic_sharpness_acum: (a: number) => number;
    readonly pocketwatchsim_mic_spl_db: (a: number) => number;
    readonly pocketwatchsim_new: () => number;
    readonly pocketwatchsim_rate_interval_hi_s_day: (a: number) => number;
    readonly pocketwatchsim_rate_interval_lo_s_day: (a: number) => number;
    readonly pocketwatchsim_rate_sensitivity_s_day_per_k: (a: number) => number;
    readonly pocketwatchsim_render_audio_block: (a: number, b: number) => [number, number];
    readonly pocketwatchsim_reset: (a: number) => void;
    readonly pocketwatchsim_set_air_jet: (a: number, b: number, c: number, d: number, e: number) => void;
    readonly pocketwatchsim_set_ambient_pressure_atm: (a: number, b: number) => void;
    readonly pocketwatchsim_set_blowtorch: (a: number, b: number, c: number, d: number, e: number) => void;
    readonly pocketwatchsim_set_compensated_hairspring: (a: number, b: number) => void;
    readonly pocketwatchsim_set_gravity_vector: (a: number, b: number, c: number, d: number) => void;
    readonly pocketwatchsim_set_liquid_nitrogen: (a: number, b: number, c: number, d: number, e: number) => void;
    readonly pocketwatchsim_set_magnet: (a: number, b: number, c: number, d: number, e: number, f: number) => void;
    readonly pocketwatchsim_set_microphone: (a: number, b: number, c: number, d: number, e: number) => void;
    readonly pocketwatchsim_set_oil_pour: (a: number, b: number, c: number, d: number, e: number) => void;
    readonly pocketwatchsim_snapshot: (a: number) => [number, number];
    readonly pocketwatchsim_spatial_grid_snapshot: (a: number) => [number, number];
    readonly pocketwatchsim_state_digest_hex: (a: number) => [number, number];
    readonly pocketwatchsim_status_label: (a: number) => [number, number];
    readonly pocketwatchsim_step: (a: number, b: number) => void;
    readonly pocketwatchsim_timegrapher_drain: (a: number) => [number, number];
    readonly pocketwatchsim_trigger_degauss: (a: number) => void;
    readonly pocketwatchsim_trigger_mechanical_shock: (a: number, b: number) => void;
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
