/* @ts-self-types="./fs_pocket_watch.d.ts" */

/**
 * Top-level coupled simulation of the pocket watch.
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
     * Air damping of the balance as a decay rate [1/s].
     * @returns {number}
     */
    balance_aero_drag_s_inv() {
        const ret = wasm.pocketwatchsim_balance_aero_drag_s_inv(this.__wbg_ptr);
        return ret;
    }
    /**
     * The case-band shell, 64 × 8 nodes × `[x, y, z, thickness, T, melt]`
     * (mm, K, 0–1).
     * @returns {Float32Array}
     */
    case_shell_snapshot() {
        const ret = wasm.pocketwatchsim_case_shell_snapshot(this.__wbg_ptr);
        var v1 = getArrayF32FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
        return v1;
    }
    /**
     * Ultrasonic clean: strip flooded oil, frost and LN₂ from the movement
     * and re-oil the jewels.
     */
    clean_excess_oil() {
        wasm.pocketwatchsim_clean_excess_oil(this.__wbg_ptr);
    }
    /**
     * Click (ratchet pawl) angle [rad].
     * @returns {number}
     */
    click_lever_angle_rad() {
        const ret = wasm.pocketwatchsim_click_lever_angle_rad(this.__wbg_ptr);
        return ret;
    }
    /**
     * Anytime-valid e-value that the watch is outside the COSC −4/+6 s/day
     * band (`fs-eproc` betting e-process).
     * @returns {number}
     */
    cosc_e_value() {
        const ret = wasm.pocketwatchsim_cosc_e_value(this.__wbg_ptr);
        return ret;
    }
    /**
     * Crown-wheel angle [rad].
     * @returns {number}
     */
    crown_wheel_angle_rad() {
        const ret = wasm.pocketwatchsim_crown_wheel_angle_rad(this.__wbg_ptr);
        return ret;
    }
    /**
     * The 16-value escapement state: θ, θ̇, φ, φ̇, ψ, ψ̇ (balance, lever and
     * escape wheel in the cage frame), the active-contact mask (bit 0 entry
     * stone, 1 exit stone, 2 pin/fork, 3 banking, 4 guard), the last impulse
     * energy [µJ], escapement efficiency, measured lift [deg], balance Q,
     * escape-wheel drive [µN m], barrel torque [mN m], power reserve [h],
     * balance pivot friction [nN m] and the air share of the losses.
     * @returns {Float64Array}
     */
    escapement_state() {
        const ret = wasm.pocketwatchsim_escapement_state(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * Half the side of the plan grid [mm].
     * @returns {number}
     */
    grid_half_span_mm() {
        const ret = wasm.pocketwatchsim_grid_half_span_mm(this.__wbg_ptr);
        return ret;
    }
    /**
     * Cells per side of the plan grid behind `spatial_grid_snapshot`.
     * @returns {number}
     */
    grid_n() {
        const ret = wasm.pocketwatchsim_grid_n(this.__wbg_ptr);
        return ret >>> 0;
    }
    /**
     * Mainspring coil radii inside the barrel [mm], innermost first.
     * @returns {Float32Array}
     */
    mainspring_coils_mm() {
        const ret = wasm.pocketwatchsim_mainspring_coils_mm(this.__wbg_ptr);
        var v1 = getArrayF32FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
        return v1;
    }
    /**
     * ISO 532-1 loudness level at the microphone [phon] (`fs-psycho`).
     * @returns {number}
     */
    mic_loudness_phon() {
        const ret = wasm.pocketwatchsim_mic_loudness_phon(this.__wbg_ptr);
        return ret;
    }
    /**
     * ISO 532-1 loudness at the microphone [sone] (`fs-psycho`).
     * @returns {number}
     */
    mic_loudness_sone() {
        const ret = wasm.pocketwatchsim_mic_loudness_sone(this.__wbg_ptr);
        return ret;
    }
    /**
     * DIN 45692 sharpness at the microphone [acum] (`fs-psycho`).
     * @returns {number}
     */
    mic_sharpness_acum() {
        const ret = wasm.pocketwatchsim_mic_sharpness_acum(this.__wbg_ptr);
        return ret;
    }
    /**
     * Sound pressure level at the microphone [dB re 20 µPa].
     * @returns {number}
     */
    mic_spl_db() {
        const ret = wasm.pocketwatchsim_mic_spl_db(this.__wbg_ptr);
        return ret;
    }
    /**
     * A regulated, wound tourbillon pocket watch at 20 °C.
     */
    constructor() {
        const ret = wasm.pocketwatchsim_new();
        this.__wbg_ptr = ret;
        PocketWatchSimFinalization.register(this, this.__wbg_ptr, this);
        return this;
    }
    /**
     * Upper end of the ±2σ interval on the measured rate [s/day].
     * @returns {number}
     */
    rate_interval_hi_s_day() {
        const ret = wasm.pocketwatchsim_rate_interval_hi_s_day(this.__wbg_ptr);
        return ret;
    }
    /**
     * Lower end of the ±2σ interval on the measured rate [s/day].
     * @returns {number}
     */
    rate_interval_lo_s_day() {
        const ret = wasm.pocketwatchsim_rate_interval_lo_s_day(this.__wbg_ptr);
        return ret;
    }
    /**
     * Rate change per kelvin of the balance and hairspring [s/day/K], from
     * the alloy's thermoelastic coefficient and the balance's expansion.
     * @returns {number}
     */
    rate_sensitivity_s_day_per_k() {
        const ret = wasm.pocketwatchsim_rate_sensitivity_s_day_per_k(this.__wbg_ptr);
        return ret;
    }
    /**
     * Render `frames` (1–8192) mono samples at 48 kHz from the microphone.
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
     * Back to the pristine watch at 20 °C (the hairspring choice is kept).
     */
    reset() {
        wasm.pocketwatchsim_reset(this.__wbg_ptr);
    }
    /**
     * Compressed-air nozzle: power 0–1 and position [mm].
     * @param {number} power
     * @param {number} x_mm
     * @param {number} y_mm
     * @param {number} z_mm
     */
    set_air_jet(power, x_mm, y_mm, z_mm) {
        wasm.pocketwatchsim_set_air_jet(this.__wbg_ptr, power, x_mm, y_mm, z_mm);
    }
    /**
     * Bell-jar pressure [atm], 0 (vacuum) to 5.
     * @param {number} pressure_atm
     */
    set_ambient_pressure_atm(pressure_atm) {
        wasm.pocketwatchsim_set_ambient_pressure_atm(this.__wbg_ptr, pressure_atm);
    }
    /**
     * Air–propane torch: power 0–1 and nozzle position in watch coordinates [mm].
     * @param {number} power
     * @param {number} x_mm
     * @param {number} y_mm
     * @param {number} z_mm
     */
    set_blowtorch(power, x_mm, y_mm, z_mm) {
        wasm.pocketwatchsim_set_blowtorch(this.__wbg_ptr, power, x_mm, y_mm, z_mm);
    }
    /**
     * Hairspring alloy: carbon steel (`false`) or a Nivarox-type
     * thermocompensating alloy (`true`). Each is regulated at the factory.
     * @param {boolean} compensated
     */
    set_compensated_hairspring(compensated) {
        wasm.pocketwatchsim_set_compensated_hairspring(this.__wbg_ptr, compensated);
    }
    /**
     * Gravity / centrifuge acceleration in g, watch frame (x to 3 o'clock,
     * y to 12 o'clock, z out of the dial).
     * @param {number} gx
     * @param {number} gy
     * @param {number} gz
     */
    set_gravity_vector(gx, gy, gz) {
        wasm.pocketwatchsim_set_gravity_vector(this.__wbg_ptr, gx, gy, gz);
    }
    /**
     * Liquid-nitrogen pour: rate 0–1 and spout position [mm].
     * @param {number} pour_rate
     * @param {number} x_mm
     * @param {number} y_mm
     * @param {number} z_mm
     */
    set_liquid_nitrogen(pour_rate, x_mm, y_mm, z_mm) {
        wasm.pocketwatchsim_set_liquid_nitrogen(this.__wbg_ptr, pour_rate, x_mm, y_mm, z_mm);
    }
    /**
     * N42 NdFeB magnet: enabled, strength (0–2, scales the remanence) and
     * position [mm].
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
     * Measurement microphone position [mm] and gain (0–4).
     * @param {number} x_mm
     * @param {number} y_mm
     * @param {number} z_mm
     * @param {number} gain
     */
    set_microphone(x_mm, y_mm, z_mm, gain) {
        wasm.pocketwatchsim_set_microphone(this.__wbg_ptr, x_mm, y_mm, z_mm, gain);
    }
    /**
     * Watch-oil dropper: rate 0–1 and nozzle position [mm].
     * @param {number} pour_rate
     * @param {number} x_mm
     * @param {number} y_mm
     * @param {number} z_mm
     */
    set_oil_pour(pour_rate, x_mm, y_mm, z_mm) {
        wasm.pocketwatchsim_set_oil_pour(this.__wbg_ptr, pour_rate, x_mm, y_mm, z_mm);
    }
    /**
     * The 176-value state snapshot (32 header values, then 12 per node).
     * @returns {Float64Array}
     */
    snapshot() {
        const ret = wasm.pocketwatchsim_snapshot(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * The plan grid, `grid_n()²` cells row-major from (−span, −span), each
     * `[T K, melt 0–1, case slump 0–1, liquid gold mm, frozen gold mm,
     * frost µm, oil µm]`.
     * @returns {Float32Array}
     */
    spatial_grid_snapshot() {
        const ret = wasm.pocketwatchsim_spatial_grid_snapshot(this.__wbg_ptr);
        var v1 = getArrayF32FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 4, 4);
        return v1;
    }
    /**
     * BLAKE3 digest of the snapshot (`fs-blake3`), for replay checks.
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
     * Status banner: which mechanism dominates the movement's behaviour.
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
     * Advance the coupled simulation by `dt_s` seconds (clamped to 0.1 ms –
     * 100 ms), in slices of at most 10 ms: the fields first (placing the
     * lever and cage for the magnet), then the movement in the environment
     * the fields leave, then the sound sources it produced.
     * @param {number} dt_s
     */
    step(dt_s) {
        wasm.pocketwatchsim_step(this.__wbg_ptr, dt_s);
    }
    /**
     * Lock drops since the last call, as `[lab_time_s, side, energy_nJ]`
     * triples (side 0 = entry, 1 = exit), for the timegrapher.
     * @returns {Float64Array}
     */
    timegrapher_drain() {
        const ret = wasm.pocketwatchsim_timegrapher_drain(this.__wbg_ptr);
        var v1 = getArrayF64FromWasm0(ret[0], ret[1]).slice();
        wasm.__wbindgen_free(ret[0], ret[1] * 8, 8);
        return v1;
    }
    /**
     * Fire the 50 Hz degaussing coil (a 3 s geometrically decaying alternating field).
     */
    trigger_degauss() {
        wasm.pocketwatchsim_trigger_degauss(this.__wbg_ptr);
    }
    /**
     * A drop or knock: a half-sine acceleration pulse of `g_peak` (10–5000 g)
     * applied to the movement and the liquids on it.
     * @param {number} g_peak
     */
    trigger_mechanical_shock(g_peak) {
        wasm.pocketwatchsim_trigger_mechanical_shock(this.__wbg_ptr, g_peak);
    }
    /**
     * Turn the crown by `delta_turns` (the crown-to-barrel gearing and the
     * click are simulated; the bridle slips past full wind).
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
