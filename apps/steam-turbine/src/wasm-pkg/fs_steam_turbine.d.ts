/* tslint:disable */
/* eslint-disable */

/**
 * Top-level interactive steam turbine model exported to native Rust & WebAssembly.
 */
export class SteamTurbineSim {
    free(): void;
    [Symbol.dispose](): void;
    /**
     * Apply new physical design parameters, rebuilding the coupled 12-component `fs-system` plant,
     * 17-stage meanline, 2D axisymmetric continuum, and 4-DOF Campbell diagram while keeping live operator controls!
     */
    apply_design_params(tip_clearance_mm: number, hp_reaction: number, lp_last_blade_h_m: number, bearing_clearance_um: number, bearing_preload: number, casing_wall_mm: number): void;
    /**
     * Export all 34 blade/vane rows as flat `[stage_idx, family_id, is_rotor, z_m, r_hub_m, r_tip_m, blade_count, axial_chord_m, has_shroud, ...]`,
     * dynamically reflecting the live `lp_last_blade_h_m` engineering design parameter!
     */
    blade_rows_flat(): Float64Array;
    /**
     * Campbell diagram data `[rpm, f1_hz, f2_hz, f3_hz, f4_hz, max_re_per_s, ...]` (24 speeds × 6).
     */
    campbell_data(): Float64Array;
    /**
     * Packed 2D axisymmetric `(64 × 24 × 8)` continuum field tensor for GPU false-color DataTextures:
     * `[T_K, sigma_vM_MPa, P_bar, Mach, q_r_kWm2, q_z_kWm2, wetness, creep_log]`.
     */
    continuum_field_flat(): Float32Array;
    /**
     * Axisymmetric 2D continuum grid bounds `[nz, nr, z_min_m, z_max_m, r_max_m]`.
     */
    continuum_grid_dims(): Float64Array;
    /**
     * Current engineer-editable design parameters `[tip_mm, hp_react, lp_h_m, brg_c_um, brg_preload, wall_mm]`.
     */
    design_params_flat(): Float64Array;
    /**
     * Certified FEM exhaust-hood structural shell eigenmodes from `fs-plate::shell::modes_shell`
     * and Levine-Schwinger diffuser radiation resistance from `fs-duct`:
     * `[f0_hz, f1_hz, f2_hz, stiffness_trace_n_m, rad_resistance, certified_from_fem]`.
     */
    fem_shell_certificate_flat(): Float64Array;
    /**
     * Export FrankenSim geometry provenance & exact Green's-theorem rotor inertia metrics:
     * `[shaft_mass_kg, blades_mass_kg, total_mass_kg, com_z_m, polar_inertia_kg_m2,
     *   transverse_inertia_kg_m2, shaft_surface_area_m2, rotor_vertices, rotor_triangles,
     *   rotor_hausdorff_m, rotor_euler_chi, rotor_watertight, total_base_tris, total_scene_tris]`.
     */
    geometry_receipt_flat(): Float64Array;
    /**
     * Last solver error string (empty when healthy).
     */
    last_error(): string;
    /**
     * The protection limits the engine applies, for the page to draw and colour against:
     * `[vib_zone_ab_um, vib_alarm_um, vib_trip_um, bearing_derate_c, bearing_alarm_c,
     *   bearing_trip_c, overspeed_trip_rpm]`.
     */
    limits(): Float64Array;
    mesh_cap_start(mesh_id: number): number;
    /**
     * Total number of exported base meshes (`24`: 17 revolved solids/cutaway shells + 7 NURBS airfoils).
     */
    mesh_count(): number;
    mesh_indices(mesh_id: number): Uint32Array;
    mesh_normals(mesh_id: number): Float32Array;
    mesh_positions(mesh_id: number): Float32Array;
    mesh_uvs(mesh_id: number): Float32Array;
    /**
     * Microphone position [m] and gain, as the engine holds them: `[x, y, z, gain]`.
     */
    mic_settings(): Float64Array;
    /**
     * Exact IAPWS-IF97 Region 4 saturation dome `[s_f, h_f, s_g, h_g, ...]` (`kJ/kg-K`, `kJ/kg`).
     */
    mollier_dome(): Float64Array;
    /**
     * Live 5-station steam expansion path on the Mollier (h-s) diagram:
     * `[s_live, h_live, s_chest, h_chest, s_ext, h_ext, s_exhaust, h_exhaust, s_exhaust_is, h_exhaust_is]`
     * in `[kJ/(kg·K), kJ/kg]`.
     */
    mollier_path(): Float64Array;
    /**
     * Create a warm, synchronized 50 MW steam turbine running at 3,000 rpm.
     */
    constructor();
    /**
     * Recent 128-point shaft whirl orbit `[x0_um, y0_um, ...]` at the front journal bearing.
     */
    orbit_history(): Float32Array;
    /**
     * Live engineering probe at any axisymmetric `(z_m, r_m)` coordinate on the 3D CAD mesh:
     * returns `[temp_c, sigma_vm_mpa, p_bar, mach, q_mag_kw_m2, wetness_pct, creep_log, region_id]`.
     */
    probe_continuum(z_m: number, r_m: number): Float64Array;
    /**
     * Force an immediate Campbell diagram recomputation using live bearing & Alford matrices.
     */
    refresh_campbell(): void;
    /**
     * Render `frames` mono 48 kHz physical audio samples in `[-1.0, 1.0]`.
     */
    render_audio_block(frames: number): Float32Array;
    /**
     * Restart the operating state: nominal warm 50 MW synchronized operation (`cold_start = false`)
     * or a cold roll-up from turning gear (`cold_start = true`). The machine's hardware (the
     * active design) and the microphone placement are kept; only the operating state restarts.
     */
    reset(cold_start: boolean): void;
    /**
     * Restore everything the page can change: the default design, the default microphone and
     * warm nominal operation.
     */
    reset_to_defaults(): void;
    reset_trips(): void;
    /**
     * Run FrankenSim's L4 ASCENT `fs-dfo::cma::cmaes` + `fs-dfo::nelder_mead` multi-variable
     * design optimizer over the live operating envelope!
     *
     * `objective_id`: `0` = Max Net Electric Efficiency, `1` = Max CHP Utilization, `2` = Max Stability & Creep Life.
     */
    run_design_optimization(objective_id: number, max_evals: number): Float64Array;
    set_condenser_air_fraction(w_air: number): void;
    set_cw_flow_frac(frac: number): void;
    set_cw_inlet_temp_c(t_c: number): void;
    set_excitation_pu(pu: number): void;
    set_extraction_demand(frac: number): void;
    set_grid_breaker(closed: boolean): void;
    set_live_steam_bar(p_bar: number): void;
    set_live_steam_temp_c(t_c: number): void;
    set_load_setpoint(frac: number): void;
    set_lube_oil_pressure_frac(frac: number): void;
    set_manual_trip(tripped: boolean): void;
    set_mic_gain(gain: number): void;
    set_mic_pos(x: number, y: number, z: number): void;
    set_oil_whip_destabilize(level: number): void;
    set_rotor_unbalance_g_mm(g_mm: number): void;
    set_speed_ref_rpm(rpm: number): void;
    /**
     * Flat 81-element telemetry snapshot array (see `snapshotLayout.ts` for field mapping).
     */
    snapshot(): Float64Array;
    /**
     * Row-by-row 17-stage meanline profile (`17 × 16 = 272` floats):
     * `[z, r_hub, r_tip, P_bar, T_C, x, Mach, eta_pct, MW, c_a, c_u1, c_u2, R, sigma_cf, sigma_vm, thrust_kN]`.
     */
    stage_profile_flat(): Float64Array;
    /**
     * Deterministic BLAKE3 content digest of the current accepted multiphysics state.
     */
    state_digest_hex(): string;
    /**
     * Advance the coupled `fs_system::System` + 4-DOF rotordynamics + 17-stage meanline + 2D continuum by `dt_s`.
     */
    step(dt_s: number): boolean;
}

export type InitInput = RequestInfo | URL | Response | BufferSource | WebAssembly.Module;

export interface InitOutput {
    readonly memory: WebAssembly.Memory;
    readonly __wbg_steamturbinesim_free: (a: number, b: number) => void;
    readonly steamturbinesim_apply_design_params: (a: number, b: number, c: number, d: number, e: number, f: number, g: number) => void;
    readonly steamturbinesim_blade_rows_flat: (a: number) => [number, number];
    readonly steamturbinesim_campbell_data: (a: number) => [number, number];
    readonly steamturbinesim_continuum_field_flat: (a: number) => [number, number];
    readonly steamturbinesim_continuum_grid_dims: (a: number) => [number, number];
    readonly steamturbinesim_design_params_flat: (a: number) => [number, number];
    readonly steamturbinesim_fem_shell_certificate_flat: (a: number) => [number, number];
    readonly steamturbinesim_geometry_receipt_flat: (a: number) => [number, number];
    readonly steamturbinesim_last_error: (a: number) => [number, number];
    readonly steamturbinesim_limits: (a: number) => [number, number];
    readonly steamturbinesim_mesh_cap_start: (a: number, b: number) => number;
    readonly steamturbinesim_mesh_count: (a: number) => number;
    readonly steamturbinesim_mesh_indices: (a: number, b: number) => [number, number];
    readonly steamturbinesim_mesh_normals: (a: number, b: number) => [number, number];
    readonly steamturbinesim_mesh_positions: (a: number, b: number) => [number, number];
    readonly steamturbinesim_mesh_uvs: (a: number, b: number) => [number, number];
    readonly steamturbinesim_mic_settings: (a: number) => [number, number];
    readonly steamturbinesim_mollier_dome: (a: number) => [number, number];
    readonly steamturbinesim_mollier_path: (a: number) => [number, number];
    readonly steamturbinesim_new: () => number;
    readonly steamturbinesim_orbit_history: (a: number) => [number, number];
    readonly steamturbinesim_probe_continuum: (a: number, b: number, c: number) => [number, number];
    readonly steamturbinesim_refresh_campbell: (a: number) => void;
    readonly steamturbinesim_render_audio_block: (a: number, b: number) => [number, number];
    readonly steamturbinesim_reset: (a: number, b: number) => void;
    readonly steamturbinesim_reset_to_defaults: (a: number) => void;
    readonly steamturbinesim_reset_trips: (a: number) => void;
    readonly steamturbinesim_run_design_optimization: (a: number, b: number, c: number) => [number, number];
    readonly steamturbinesim_set_condenser_air_fraction: (a: number, b: number) => void;
    readonly steamturbinesim_set_cw_flow_frac: (a: number, b: number) => void;
    readonly steamturbinesim_set_cw_inlet_temp_c: (a: number, b: number) => void;
    readonly steamturbinesim_set_excitation_pu: (a: number, b: number) => void;
    readonly steamturbinesim_set_extraction_demand: (a: number, b: number) => void;
    readonly steamturbinesim_set_grid_breaker: (a: number, b: number) => void;
    readonly steamturbinesim_set_live_steam_bar: (a: number, b: number) => void;
    readonly steamturbinesim_set_live_steam_temp_c: (a: number, b: number) => void;
    readonly steamturbinesim_set_load_setpoint: (a: number, b: number) => void;
    readonly steamturbinesim_set_lube_oil_pressure_frac: (a: number, b: number) => void;
    readonly steamturbinesim_set_manual_trip: (a: number, b: number) => void;
    readonly steamturbinesim_set_mic_gain: (a: number, b: number) => void;
    readonly steamturbinesim_set_mic_pos: (a: number, b: number, c: number, d: number) => void;
    readonly steamturbinesim_set_oil_whip_destabilize: (a: number, b: number) => void;
    readonly steamturbinesim_set_rotor_unbalance_g_mm: (a: number, b: number) => void;
    readonly steamturbinesim_set_speed_ref_rpm: (a: number, b: number) => void;
    readonly steamturbinesim_snapshot: (a: number) => [number, number];
    readonly steamturbinesim_stage_profile_flat: (a: number) => [number, number];
    readonly steamturbinesim_state_digest_hex: (a: number) => [number, number];
    readonly steamturbinesim_step: (a: number, b: number) => number;
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
