/* @ts-self-types="./fs_pocket_watch.d.ts" */

/**
 * Top-level coupled simulation engine for the FrankenSim Grand Complication Pocket Watch.
 */
export class PocketWatchSim {
    __destroy_into_raw() {
        const ptr = this.__wbg_ptr;
        this.__wbg_ptr = 0;
        PocketWatchSimFinalization.unregister(this);
        return ptr;
    }
    free() {
        const ptr = this.__destroy_into_raw();
        wasm.__wbg_pocketwatchsim_free(ptr, 0);
    }
    /**
     * Stokes-Stuart oscillating-cylinder aerodynamic drag on the balance wheel ($\text{s}^{-1}$) (`fs-vkplate::air_loss`).
     * @returns {number}
     */
    balance_aero_drag_s_inv() {
        const ret = wasm.pocketwatchsim_balance_aero_drag_s_inv(this.__wbg_ptr);
        return ret;
    }
    /**
     * Return the 64x8 (512-node) 3D Thermo-Elasto-Viscoplastic + Surface-Tension Case Shell snapshot
     * (`512 * 6 = 3072` `f32`s): `[x_mm, y_mm, z_mm, thickness_mm, temperature_k, melt_fraction]`
     * for each 3D Lagrangian node `(it, iz)` on the Bassine 18k gold case.
     * @returns {Float32Array}
     */
    case_shell_snapshot() {
        const ret = wasm.pocketwatchsim_case_shell_snapshot(this.__wbg_ptr);
        var v1 = getArrayF32FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
        return v1;
    }
    /**
     * Ultrasonic clean excess flooded oil and restore factory $1.6\text{ }\mu\text{m}$ jewel lubrication.
     */
    clean_excess_oil() {
        wasm.pocketwatchsim_clean_excess_oil(this.__wbg_ptr);
    }
    /**
     * Return the dynamic recoil click lever (*cliquet*) deflection angle ($\text{rad}$).
     * @returns {number}
     */
    click_lever_angle_rad() {
        const ret = wasm.pocketwatchsim_click_lever_angle_rad(this.__wbg_ptr);
        return ret;
    }
    /**
     * Anytime-valid COSC chronometer compliance e-process value $E_t$ (`fs-eproc`).
     * @returns {number}
     */
    cosc_e_value() {
        const ret = wasm.pocketwatchsim_cosc_e_value(this.__wbg_ptr);
        return ret;
    }
    /**
     * Return the keyless works crown wheel rotation angle ($\text{rad}$).
     * @returns {number}
     */
    crown_wheel_angle_rad() {
        const ret = wasm.pocketwatchsim_crown_wheel_angle_rad(this.__wbg_ptr);
        return ret;
    }
    /**
     * Lagrange cost-optimal GD&T manufacturing tolerance allocation cost (`fs-toleralloc`).
     * @returns {number}
     */
    gdt_allocated_cost() {
        const ret = wasm.pocketwatchsim_gdt_allocated_cost(this.__wbg_ptr);
        return ret;
    }
    /**
     * Return the 18 discretized elastica mainspring coil radii ($\text{mm}$) inside the barrel.
     * @returns {Float32Array}
     */
    mainspring_coils_mm() {
        const ret = wasm.pocketwatchsim_mainspring_coils_mm(this.__wbg_ptr);
        var v1 = getArrayF32FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
        return v1;
    }
    /**
     * ISO 532-1 loudness level in phon at the 3D virtual microphone (`fs-psycho`).
     * @returns {number}
     */
    mic_loudness_phon() {
        const ret = wasm.pocketwatchsim_mic_loudness_phon(this.__wbg_ptr);
        return ret;
    }
    /**
     * ISO 532-1 Zwicker stationary loudness in sone at the 3D virtual microphone (`fs-psycho`).
     * @returns {number}
     */
    mic_loudness_sone() {
        const ret = wasm.pocketwatchsim_mic_loudness_sone(this.__wbg_ptr);
        return ret;
    }
    /**
     * DIN 45692 psychoacoustic sharpness in acum at the 3D virtual microphone (`fs-psycho`).
     * @returns {number}
     */
    mic_sharpness_acum() {
        const ret = wasm.pocketwatchsim_mic_sharpness_acum(this.__wbg_ptr);
        return ret;
    }
    /**
     * Sound Pressure Level ($\text{dB SPL re } 20\text{ }\mu\text{Pa}$) at the 3D virtual microphone (`fs-psycho`).
     * @returns {number}
     */
    mic_spl_db() {
        const ret = wasm.pocketwatchsim_mic_spl_db(this.__wbg_ptr);
        return ret;
    }
    /**
     * Multi-crate spectral, topological, regime, UQ, and port-Hamiltonian coupling audit summary
     * exercising `fs-couple`, `fs-regime`, `fs-cheb`, `fs-fft`, `fs-orbit`, `fs-asbuilt`,
     * `fs-rep-sdf`, `fs-rep-mesh`, `fs-topo`, `fs-vfit`, `fs-uq`, and `fs-viz`.
     * @returns {Float64Array}
     */
    multiphysics_crate_audit_summary() {
        const ret = wasm.pocketwatchsim_multiphysics_crate_audit_summary(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * Construct a new, chronometer-calibrated Patek-style Tourbillon pocket watch at $293.15\text{ K}$.
     */
    constructor() {
        const ret = wasm.pocketwatchsim_new();
        this.__wbg_ptr = ret;
        PocketWatchSimFinalization.register(this, this.__wbg_ptr, this);
        return this;
    }
    /**
     * Upper bound of certified interval rate enclosure ($\text{s/day}$) (`fs-ivl`).
     * @returns {number}
     */
    rate_interval_hi_s_day() {
        const ret = wasm.pocketwatchsim_rate_interval_hi_s_day(this.__wbg_ptr);
        return ret;
    }
    /**
     * Lower bound of certified interval rate enclosure ($\text{s/day}$) (`fs-ivl`).
     * @returns {number}
     */
    rate_interval_lo_s_day() {
        const ret = wasm.pocketwatchsim_rate_interval_lo_s_day(this.__wbg_ptr);
        return ret;
    }
    /**
     * Exact forward-mode AD thermoelastic rate sensitivity $\partial r / \partial T$ ($\text{s/day/K}$) (`fs-ad`).
     * @returns {number}
     */
    rate_sensitivity_s_day_per_k() {
        const ret = wasm.pocketwatchsim_rate_sensitivity_s_day_per_k(this.__wbg_ptr);
        return ret;
    }
    /**
     * Render `frames` mono audio samples at $48\text{ kHz}$ from the 3D virtual microphone
     * and return them as a `Vec<f32>` (converted to `Float32Array` across `wasm-bindgen`).
     * @param {number} frames
     * @returns {Float32Array}
     */
    render_audio_block(frames) {
        const ret = wasm.pocketwatchsim_render_audio_block(this.__wbg_ptr, frames);
        var v1 = getArrayF32FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
        return v1;
    }
    /**
     * Reset the entire watch and all materials back to pristine $20.0^\circ\text{C}$ chronometer state.
     */
    reset() {
        wasm.pocketwatchsim_reset(this.__wbg_ptr);
    }
    /**
     * Configure the high-velocity compressed air / argon nozzle ($300\text{ K}, M \approx 0.35$) in 3D watch coordinates ($\text{mm}$).
     * Applies forced convective cooling and aerodynamic surface shear that blows molten gold and watch oil across the movement.
     * @param {number} power
     * @param {number} x_mm
     * @param {number} y_mm
     * @param {number} z_mm
     */
    set_air_jet(power, x_mm, y_mm, z_mm) {
        wasm.pocketwatchsim_set_air_jet(this.__wbg_ptr, power, x_mm, y_mm, z_mm);
    }
    /**
     * Set the ambient bell-jar chamber pressure in atmospheres ($0.0\text{ atm}$ vacuum $\to 1.0\text{ atm} \to 5.0\text{ atm}$ hyperbaric).
     * @param {number} pressure_atm
     */
    set_ambient_pressure_atm(pressure_atm) {
        wasm.pocketwatchsim_set_ambient_pressure_atm(this.__wbg_ptr, pressure_atm);
    }
    /**
     * Configure the virtual oxy-butane blowtorch ($T_{\text{flame}} = 2250\text{ K}$) in 3D watch coordinates ($\text{mm}$).
     * @param {number} power
     * @param {number} x_mm
     * @param {number} y_mm
     * @param {number} z_mm
     */
    set_blowtorch(power, x_mm, y_mm, z_mm) {
        wasm.pocketwatchsim_set_blowtorch(this.__wbg_ptr, power, x_mm, y_mm, z_mm);
    }
    /**
     * Toggle between traditional High-Carbon Blued Steel (`false`) and Compensated Invar/Fe-36Ni (`true`) hairspring.
     * @param {boolean} compensated
     */
    set_compensated_hairspring(compensated) {
        wasm.pocketwatchsim_set_compensated_hairspring(this.__wbg_ptr, compensated);
    }
    /**
     * Set the 3D gravitational / centrifuge acceleration vector $[g_x, g_y, g_z]$ in Earth-$g$ units ($9.80665\text{ m/s}^2$).
     * @param {number} gx
     * @param {number} gy
     * @param {number} gz
     */
    set_gravity_vector(gx, gy, gz) {
        wasm.pocketwatchsim_set_gravity_vector(this.__wbg_ptr, gx, gy, gz);
    }
    /**
     * Configure the cryogenic liquid nitrogen ($\text{LN}_2$ at $77.36\text{ K}$) pour stream in 3D watch coordinates ($\text{mm}$).
     * @param {number} pour_rate
     * @param {number} x_mm
     * @param {number} y_mm
     * @param {number} z_mm
     */
    set_liquid_nitrogen(pour_rate, x_mm, y_mm, z_mm) {
        wasm.pocketwatchsim_set_liquid_nitrogen(this.__wbg_ptr, pour_rate, x_mm, y_mm, z_mm);
    }
    /**
     * Configure the external sintered $\text{NdFeB N42}$ super-magnet ($B_r = 1.32\text{ T}$) in 3D watch coordinates ($\text{mm}$).
     * @param {boolean} enabled
     * @param {number} strength
     * @param {number} x_mm
     * @param {number} y_mm
     * @param {number} z_mm
     */
    set_magnet(enabled, strength, x_mm, y_mm, z_mm) {
        wasm.pocketwatchsim_set_magnet(this.__wbg_ptr, enabled, strength, x_mm, y_mm, z_mm);
    }
    /**
     * Position the 3D studio measurement microphone in watch coordinates ($\text{mm}$) and set gain.
     * @param {number} x_mm
     * @param {number} y_mm
     * @param {number} z_mm
     * @param {number} gain
     */
    set_microphone(x_mm, y_mm, z_mm, gain) {
        wasm.pocketwatchsim_set_microphone(this.__wbg_ptr, x_mm, y_mm, z_mm, gain);
    }
    /**
     * Configure the watch oil dropper / flood stream (`pennzane-shf-x-2000-bearing-oil`) in 3D watch coordinates ($\text{mm}$).
     * @param {number} pour_rate
     * @param {number} x_mm
     * @param {number} y_mm
     * @param {number} z_mm
     */
    set_oil_pour(pour_rate, x_mm, y_mm, z_mm) {
        wasm.pocketwatchsim_set_oil_pour(this.__wbg_ptr, pour_rate, x_mm, y_mm, z_mm);
    }
    /**
     * Return the current 176-element `f64` state snapshot for 3D rendering and telemetry.
     * @returns {Float64Array}
     */
    snapshot() {
        const ret = wasm.pocketwatchsim_snapshot(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * Return the 24x24 (576-cell) 2D Finite-Volume Spatial PDE grid snapshot (`576 * 7 = 4032` `f32`s):
     * `[temperature_k, melt_fraction, case_slump_frac, liquid_gold_mm, frozen_gold_mm, frost_um, oil_um]`
     * for each cell `(ix, iy)` in row-major order.
     * @returns {Float32Array}
     */
    spatial_grid_snapshot() {
        const ret = wasm.pocketwatchsim_spatial_grid_snapshot(this.__wbg_ptr);
        var v1 = getArrayF32FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
        return v1;
    }
    /**
     * Deterministic BLAKE3 digest of the current 176-float state snapshot (`fs-blake3`).
     * @returns {string}
     */
    state_digest_hex() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.pocketwatchsim_state_digest_hex(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * Human-readable status banner for the current horological & multiphysics state.
     * @returns {string}
     */
    status_label() {
        let deferred1_0;
        let deferred1_1;
        try {
            const ret = wasm.pocketwatchsim_status_label(this.__wbg_ptr);
            deferred1_0 = ret[0];
            deferred1_1 = ret[1];
            return getStringFromWasm0(ret[0], ret[1]);
        } finally {
            wasm.__wbindgen_free(deferred1_0, deferred1_1, 1);
        }
    }
    /**
     * Advance the coupled multiphysics and horological simulation by `dt_s` seconds.
     * @param {number} dt_s
     */
    step(dt_s) {
        wasm.pocketwatchsim_step(this.__wbg_ptr, dt_s);
    }
    /**
     * Fire the alternating-field degaussing coil to wipe out residual ferromagnetic remanence $M_{\text{rem}}$.
     */
    trigger_degauss() {
        wasm.pocketwatchsim_trigger_degauss(this.__wbg_ptr);
    }
    /**
     * Deliver a mechanical bench-drop / G-shock impulse ($G_{\text{peak}}$ in $\text{g}$-units, e.g. $500\text{ G}$)
     * to the Incabloc ruby jewel suspension, balance wheel, 18k gold case shell, and molten gold pool.
     * @param {number} g_peak
     */
    trigger_mechanical_shock(g_peak) {
        wasm.pocketwatchsim_trigger_mechanical_shock(this.__wbg_ptr, g_peak);
    }
    /**
     * Wind the mainspring crown by `delta_turns` (where `8.5` turns = full 42-hour wind).
     * @param {number} delta_turns
     */
    wind_crown(delta_turns) {
        wasm.pocketwatchsim_wind_crown(this.__wbg_ptr, delta_turns);
    }
}
if (Symbol.dispose) PocketWatchSim.prototype[Symbol.dispose] = PocketWatchSim.prototype.free;
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
        "./fs_pocket_watch_bg.js": import0,
    };
}

const PocketWatchSimFinalization = (typeof FinalizationRegistry === 'undefined')
    ? { register: () => {}, unregister: () => {} }
    : new FinalizationRegistry(ptr => wasm.__wbg_pocketwatchsim_free(ptr, 1));

function getArrayF32FromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return getFloat32ArrayMemory0().subarray(ptr / 4, ptr / 4 + len);
}

function getArrayF64FromWasm0(ptr, len) {
    ptr = ptr >>> 0;
    return getFloat64ArrayMemory0().subarray(ptr / 8, ptr / 8 + len);
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
        module_or_path = new URL('fs_pocket_watch_bg.wasm', import.meta.url);
    }
    const imports = __wbg_get_imports();

    if (typeof module_or_path === 'string' || (typeof Request === 'function' && module_or_path instanceof Request) || (typeof URL === 'function' && module_or_path instanceof URL)) {
        module_or_path = fetch(module_or_path);
    }

    const { instance, module } = await __wbg_load(await module_or_path, imports);

    return __wbg_finalize_init(instance, module);
}

export { initSync, __wbg_init as default };
