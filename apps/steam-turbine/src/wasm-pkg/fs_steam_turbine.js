/* @ts-self-types="./fs_steam_turbine.d.ts" */

/**
 * Top-level interactive Siemens SST-600 digital twin exported to native Rust & WebAssembly.
 */
export class SteamTurbineSim {
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        SteamTurbineSimFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_steamturbinesim_free(ptr, 0);
    }
    /**
     * Apply new physical design parameters, rebuilding the coupled 12-component `fs-system` plant,
     * 17-stage meanline, 2D axisymmetric continuum, and 4-DOF Campbell diagram while keeping live operator controls!
     * @param {number} tip_clearance_mm
     * @param {number} hp_reaction
     * @param {number} lp_last_blade_h_m
     * @param {number} bearing_clearance_um
     * @param {number} bearing_preload
     * @param {number} casing_wall_mm
     */
    apply_design_params(tip_clearance_mm, hp_reaction, lp_last_blade_h_m, bearing_clearance_um, bearing_preload, casing_wall_mm) {
        wasm.steamturbinesim_apply_design_params(this.__wbg_ptr, tip_clearance_mm, hp_reaction, lp_last_blade_h_m, bearing_clearance_um, bearing_preload, casing_wall_mm);
    }
    /**
     * Export all 34 blade/vane rows as flat `[stage_idx, family_id, is_rotor, z_m, r_hub_m, r_tip_m, blade_count, axial_chord_m, has_shroud, ...]`,
     * dynamically reflecting the live `lp_last_blade_h_m` engineering design parameter!
     * @returns {Float64Array}
     */
    blade_rows_flat() {
        const ret = wasm.steamturbinesim_blade_rows_flat(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * Campbell diagram data `[rpm, f1_hz, f2_hz, f3_hz, f4_hz, max_re_per_s, ...]` (24 speeds × 6).
     * @returns {Float64Array}
     */
    campbell_data() {
        const ret = wasm.steamturbinesim_campbell_data(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * Packed 2D axisymmetric `(64 × 24 × 8)` continuum field tensor for GPU false-color DataTextures:
     * `[T_K, sigma_vM_MPa, P_bar, Mach, q_r_kWm2, q_z_kWm2, wetness, creep_log]`.
     * @returns {Float32Array}
     */
    continuum_field_flat() {
        const ret = wasm.steamturbinesim_continuum_field_flat(this.__wbg_ptr);
        var v1 = getArrayF32FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
        return v1;
    }
    /**
     * Axisymmetric 2D continuum grid bounds `[nz, nr, z_min_m, z_max_m, r_max_m]`.
     * @returns {Float64Array}
     */
    continuum_grid_dims() {
        const ret = wasm.steamturbinesim_continuum_grid_dims(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * Current engineer-editable design parameters `[tip_mm, hp_react, lp_h_m, brg_c_um, brg_preload, wall_mm]`.
     * @returns {Float64Array}
     */
    design_params_flat() {
        const ret = wasm.steamturbinesim_design_params_flat(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * Certified FEM exhaust-hood structural shell eigenmodes from `fs-plate::shell::modes_shell`
     * and Levine-Schwinger diffuser radiation resistance from `fs-duct`:
     * `[f0_hz, f1_hz, f2_hz, stiffness_trace_n_m, rad_resistance, certified_from_fem]`.
     * @returns {Float64Array}
     */
    fem_shell_certificate_flat() {
        const ret = wasm.steamturbinesim_fem_shell_certificate_flat(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * Export FrankenSim geometry provenance & exact Green's-theorem rotor inertia metrics:
     * `[shaft_mass_kg, blades_mass_kg, total_mass_kg, com_z_m, polar_inertia_kg_m2,
     *   transverse_inertia_kg_m2, shaft_surface_area_m2, rotor_vertices, rotor_triangles,
     *   rotor_hausdorff_m, rotor_euler_chi, rotor_watertight, total_base_tris, total_scene_tris]`.
     * @returns {Float64Array}
     */
    geometry_receipt_flat() {
        const ret = wasm.steamturbinesim_geometry_receipt_flat(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * Last solver error string (empty when healthy).
     * @returns {string}
     */
    last_error() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.steamturbinesim_last_error(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * @param {number} mesh_id
     * @returns {number}
     */
    mesh_cap_start(mesh_id) {
        const ret = wasm.steamturbinesim_mesh_cap_start(this.__wbg_ptr, mesh_id);
        return ret >>> 0;
    }
    /**
     * Total number of exported base meshes (`24`: 17 revolved solids/cutaway shells + 7 NURBS airfoils).
     * @returns {number}
     */
    mesh_count() {
        const ret = wasm.steamturbinesim_mesh_count(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * @param {number} mesh_id
     * @returns {Uint32Array}
     */
    mesh_indices(mesh_id) {
        const ret = wasm.steamturbinesim_mesh_indices(this.__wbg_ptr, mesh_id);
        var v1 = getArrayU32FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
        return v1;
    }
    /**
     * @param {number} mesh_id
     * @returns {Float32Array}
     */
    mesh_normals(mesh_id) {
        const ret = wasm.steamturbinesim_mesh_normals(this.__wbg_ptr, mesh_id);
        var v1 = getArrayF32FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
        return v1;
    }
    /**
     * @param {number} mesh_id
     * @returns {Float32Array}
     */
    mesh_positions(mesh_id) {
        const ret = wasm.steamturbinesim_mesh_positions(this.__wbg_ptr, mesh_id);
        var v1 = getArrayF32FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
        return v1;
    }
    /**
     * @param {number} mesh_id
     * @returns {Float32Array}
     */
    mesh_uvs(mesh_id) {
        const ret = wasm.steamturbinesim_mesh_uvs(this.__wbg_ptr, mesh_id);
        var v1 = getArrayF32FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
        return v1;
    }
    /**
     * Exact IAPWS-IF97 Region 4 saturation dome `[s_f, h_f, s_g, h_g, ...]` (`kJ/kg-K`, `kJ/kg`).
     * @returns {Float64Array}
     */
    mollier_dome() {
        const ret = wasm.steamturbinesim_mollier_dome(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * Live 5-station steam expansion path on the Mollier (h-s) diagram:
     * `[s_live, h_live, s_chest, h_chest, s_ext, h_ext, s_exhaust, h_exhaust, s_exhaust_is, h_exhaust_is]`
     * in `[kJ/(kg·K), kJ/kg]`.
     * @returns {Float64Array}
     */
    mollier_path() {
        const ret = wasm.steamturbinesim_mollier_path(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * Create a warm, synchronized 50 MW Siemens SST-600 steam turbine running at 3,000 rpm.
     */
    constructor() {
        const ret = wasm.steamturbinesim_new();
        this.__wbg_ptr = ret;
        SteamTurbineSimFinalization.register(this, this.__wbg_ptr, this);
        return this;
    }
    /**
     * Recent 128-point shaft whirl orbit `[x0_um, y0_um, ...]` at the front journal bearing.
     * @returns {Float32Array}
     */
    orbit_history() {
        const ret = wasm.steamturbinesim_orbit_history(this.__wbg_ptr);
        var v1 = getArrayF32FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
        return v1;
    }
    /**
     * Live engineering probe at any axisymmetric `(z_m, r_m)` coordinate on the 3D CAD mesh:
     * returns `[temp_c, sigma_vm_mpa, p_bar, mach, q_mag_kw_m2, wetness_pct, creep_log, region_id]`.
     * @param {number} z_m
     * @param {number} r_m
     * @returns {Float64Array}
     */
    probe_continuum(z_m, r_m) {
        const ret = wasm.steamturbinesim_probe_continuum(this.__wbg_ptr, z_m, r_m);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * Force an immediate Campbell diagram recomputation using live bearing & Alford matrices.
     */
    refresh_campbell() {
        wasm.steamturbinesim_refresh_campbell(this.__wbg_ptr);
    }
    /**
     * Render `frames` mono 48 kHz physical audio samples in `[-1.0, 1.0]`.
     * @param {number} frames
     * @returns {Float32Array}
     */
    render_audio_block(frames) {
        const ret = wasm.steamturbinesim_render_audio_block(this.__wbg_ptr, frames);
        var v1 = getArrayF32FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
        return v1;
    }
    /**
     * Reset to nominal warm 50 MW synchronized operation (`cold_start = false`) or
     * cold start roll-up from turning gear (`cold_start = true`).
     * @param {boolean} cold_start
     */
    reset(cold_start) {
        wasm.steamturbinesim_reset(this.__wbg_ptr, cold_start);
    }
    reset_trips() {
        wasm.steamturbinesim_reset_trips(this.__wbg_ptr);
    }
    /**
     * Run FrankenSim's L4 ASCENT `fs-dfo::cma::cmaes` + `fs-dfo::nelder_mead` multi-variable
     * design optimizer over the live operating envelope!
     *
     * `objective_id`: `0` = Max Net Electric Efficiency, `1` = Max CHP Utilization, `2` = Max Stability & Creep Life.
     * @param {number} objective_id
     * @param {number} max_evals
     * @returns {Float64Array}
     */
    run_design_optimization(objective_id, max_evals) {
        const ret = wasm.steamturbinesim_run_design_optimization(this.__wbg_ptr, objective_id, max_evals);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * @param {number} w_air
     */
    set_condenser_air_fraction(w_air) {
        wasm.steamturbinesim_set_condenser_air_fraction(this.__wbg_ptr, w_air);
    }
    /**
     * @param {number} frac
     */
    set_cw_flow_frac(frac) {
        wasm.steamturbinesim_set_cw_flow_frac(this.__wbg_ptr, frac);
    }
    /**
     * @param {number} t_c
     */
    set_cw_inlet_temp_c(t_c) {
        wasm.steamturbinesim_set_cw_inlet_temp_c(this.__wbg_ptr, t_c);
    }
    /**
     * @param {number} pu
     */
    set_excitation_pu(pu) {
        wasm.steamturbinesim_set_excitation_pu(this.__wbg_ptr, pu);
    }
    /**
     * @param {number} frac
     */
    set_extraction_demand(frac) {
        wasm.steamturbinesim_set_extraction_demand(this.__wbg_ptr, frac);
    }
    /**
     * @param {boolean} closed
     */
    set_grid_breaker(closed) {
        wasm.steamturbinesim_set_grid_breaker(this.__wbg_ptr, closed);
    }
    /**
     * @param {number} p_bar
     */
    set_live_steam_bar(p_bar) {
        wasm.steamturbinesim_set_live_steam_bar(this.__wbg_ptr, p_bar);
    }
    /**
     * @param {number} t_c
     */
    set_live_steam_temp_c(t_c) {
        wasm.steamturbinesim_set_live_steam_temp_c(this.__wbg_ptr, t_c);
    }
    /**
     * @param {number} frac
     */
    set_load_setpoint(frac) {
        wasm.steamturbinesim_set_load_setpoint(this.__wbg_ptr, frac);
    }
    /**
     * @param {number} frac
     */
    set_lube_oil_pressure_frac(frac) {
        wasm.steamturbinesim_set_lube_oil_pressure_frac(this.__wbg_ptr, frac);
    }
    /**
     * @param {boolean} tripped
     */
    set_manual_trip(tripped) {
        wasm.steamturbinesim_set_manual_trip(this.__wbg_ptr, tripped);
    }
    /**
     * @param {number} gain
     */
    set_mic_gain(gain) {
        wasm.steamturbinesim_set_mic_gain(this.__wbg_ptr, gain);
    }
    /**
     * @param {number} x
     * @param {number} y
     * @param {number} z
     */
    set_mic_pos(x, y, z) {
        wasm.steamturbinesim_set_mic_pos(this.__wbg_ptr, x, y, z);
    }
    /**
     * @param {number} level
     */
    set_oil_whip_destabilize(level) {
        wasm.steamturbinesim_set_oil_whip_destabilize(this.__wbg_ptr, level);
    }
    /**
     * @param {number} g_mm
     */
    set_rotor_unbalance_g_mm(g_mm) {
        wasm.steamturbinesim_set_rotor_unbalance_g_mm(this.__wbg_ptr, g_mm);
    }
    /**
     * @param {number} rpm
     */
    set_speed_ref_rpm(rpm) {
        wasm.steamturbinesim_set_speed_ref_rpm(this.__wbg_ptr, rpm);
    }
    /**
     * Flat 80-element telemetry snapshot array (see `snapshotLayout.ts` for field mapping).
     * @returns {Float64Array}
     */
    snapshot() {
        const ret = wasm.steamturbinesim_snapshot(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * Row-by-row 17-stage meanline profile (`17 × 16 = 272` floats):
     * `[z, r_hub, r_tip, P_bar, T_C, x, Mach, eta_pct, MW, c_a, c_u1, c_u2, R, sigma_cf, sigma_vm, thrust_kN]`.
     * @returns {Float64Array}
     */
    stage_profile_flat() {
        const ret = wasm.steamturbinesim_stage_profile_flat(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * Deterministic BLAKE3 content digest of the current accepted multiphysics state.
     * @returns {string}
     */
    state_digest_hex() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.steamturbinesim_state_digest_hex(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * Advance the coupled `fs_system::System` + 4-DOF rotordynamics + 17-stage meanline + 2D continuum by `dt_s`.
     * @param {number} dt_s
     * @returns {boolean}
     */
    step(dt_s) {
        const ret = wasm.steamturbinesim_step(this.__wbg_ptr, dt_s);
        return ret !== 0;
    }
}
if (Symbol.dispose) SteamTurbineSim.prototype[Symbol.dispose] = SteamTurbineSim.prototype.free;
function __wbg_get_imports() {
    const import0 = {
        __proto__: null,
        __wbg___wbindgen_throw_5d9e815e6fdf150f: function(arg0, arg1) {
            throw new Error(getStringFromWasm0(arg0, arg1));
        },
        __wbindgen_init_externref_table: function() {
            const table = wasm.__wbindgen_externrefs;
            const offset = table.grow(4);
            table.set(0, undefined);
            table.set(offset + 0, undefined);
            table.set(offset + 1, null);
            table.set(offset + 2, true);
            table.set(offset + 3, false);
        },
    };
    return {
        __proto__: null,
        "./fs_steam_turbine_bg.js": import0,
    };
}

const SteamTurbineSimFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_steamturbinesim_free(ptr, 1));

function getArrayF32FromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return getFloat32ArrayMemory0().subarray(ptr / 4, ptr / 4 + len);
}

function getArrayF64FromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return getFloat64ArrayMemory0().subarray(ptr / 8, ptr / 8 + len);
}

function getArrayU32FromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return getUint32ArrayMemory0().subarray(ptr / 4, ptr / 4 + len);
}

let cachedFloat32ArrayMemory0 = null;
function getFloat32ArrayMemory0() {
    if (cachedFloat32ArrayMemory0 === null || cachedFloat32ArrayMemory0.byteLength === 0) {
        cachedFloat32ArrayMemory0 = new Float32Array(wasm.memory.buffer);
    }
    return cachedFloat32ArrayMemory0;
}

let cachedFloat64ArrayMemory0 = null;
function getFloat64ArrayMemory0() {
    if (cachedFloat64ArrayMemory0 === null || cachedFloat64ArrayMemory0.byteLength === 0) {
        cachedFloat64ArrayMemory0 = new Float64Array(wasm.memory.buffer);
    }
    return cachedFloat64ArrayMemory0;
}

function getStringFromWasm0(ptr, len) {
    return decodeText(ptr >>> 0, len);
}

let cachedUint32ArrayMemory0 = null;
function getUint32ArrayMemory0() {
    if (cachedUint32ArrayMemory0 === null || cachedUint32ArrayMemory0.byteLength === 0) {
        cachedUint32ArrayMemory0 = new Uint32Array(wasm.memory.buffer);
    }
    return cachedUint32ArrayMemory0;
}

let cachedUint8ArrayMemory0 = null;
function getUint8ArrayMemory0() {
    if (cachedUint8ArrayMemory0 === null || cachedUint8ArrayMemory0.byteLength === 0) {
        cachedUint8ArrayMemory0 = new Uint8Array(wasm.memory.buffer);
    }
    return cachedUint8ArrayMemory0;
}

let cachedTextDecoder = new TextDecoder('utf-8', { ignoreBOM: true, fatal: true });
cachedTextDecoder.decode();
const MAX_SAFARI_DECODE_BYTES = 2146435072;
let numBytesDecoded = 0;
function decodeText(ptr, len) {
    numBytesDecoded += len;
    if (numBytesDecoded >= MAX_SAFARI_DECODE_BYTES) {
        cachedTextDecoder = new TextDecoder('utf-8', { ignoreBOM: true, fatal: true });
        cachedTextDecoder.decode();
        numBytesDecoded = len;
    }
    return cachedTextDecoder.decode(getUint8ArrayMemory0().subarray(ptr, ptr + len));
}

let wasmModule, wasmInstance, wasm;
function __wbg_finalize_init(instance, module) {
    wasmInstance = instance;
    wasm = instance.exports;
    wasmModule = module;
    cachedFloat32ArrayMemory0 = null;
    cachedFloat64ArrayMemory0 = null;
    cachedUint32ArrayMemory0 = null;
    cachedUint8ArrayMemory0 = null;
    wasm.__wbindgen_start();
    return wasm;
}

async function __wbg_load(module, imports) {
    if (typeof Response === 'function' && module instanceof Response) {
        if (!module.ok) {
            throw new Error(`failed to fetch Wasm: ${module.status} ${module.statusText} fetching '${module.url}'`);
        }

        if (typeof WebAssembly.instantiateStreaming === 'function') {
            try {
                return await WebAssembly.instantiateStreaming(module, imports);
            } catch (e) {
                const validResponse = expectedResponseType(module.type);

                if (validResponse && module.headers.get('Content-Type') !== 'application/wasm') {
                    console.warn("`WebAssembly.instantiateStreaming` failed because your server does not serve Wasm with `application/wasm` MIME type. Falling back to `WebAssembly.instantiate` which is slower. Original error:\n", e);

                } else { throw e; }
            }
        }

        const bytes = await module.arrayBuffer();
        return await WebAssembly.instantiate(bytes, imports);
    } else {
        const instance = await WebAssembly.instantiate(module, imports);

        if (instance instanceof WebAssembly.Instance) {
            return { instance, module };
        } else {
            return instance;
        }
    }

    function expectedResponseType(type) {
        switch (type) {
            case 'basic': case 'cors': case 'default': return true;
        }
        return false;
    }
}

function initSync(module) {
    if (wasm !== undefined) return wasm;


    if (module !== undefined) {
        if (Object.getPrototypeOf(module) === Object.prototype) {
            ({module} = module)
        } else {
            console.warn('using deprecated parameters for `initSync()`; pass a single object instead')
        }
    }

    const imports = __wbg_get_imports();
    if (!(module instanceof WebAssembly.Module)) {
        module = new WebAssembly.Module(module);
    }
    const instance = new WebAssembly.Instance(module, imports);
    return __wbg_finalize_init(instance, module);
}

async function __wbg_init(module_or_path) {
    if (wasm !== undefined) return wasm;


    if (module_or_path !== undefined) {
        if (Object.getPrototypeOf(module_or_path) === Object.prototype) {
            ({module_or_path} = module_or_path)
        } else {
            console.warn('using deprecated parameters for the initialization function; pass a single object instead')
        }
    }

    if (module_or_path === undefined) {
        module_or_path = new URL('fs_steam_turbine_bg.wasm', import.meta.url);
    }
    const imports = __wbg_get_imports();

    if (typeof module_or_path === 'string' || (typeof Request === 'function' && module_or_path instanceof Request) || (typeof URL === 'function' && module_or_path instanceof URL)) {
        module_or_path = fetch(module_or_path);
    }

    const { instance, module } = await __wbg_load(await module_or_path, imports);

    return __wbg_finalize_init(instance, module);
}

export { initSync, __wbg_init as default };
