use fs_flux::lc::{TeslaCoilParams, step_tesla_coil};

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

#[cfg_attr(target_arch = "wasm32", wasm_bindgen)]
pub fn tesla_coil_step(
    resonant_freq_khz: f64,
    input_kv: f64,
    spark_gap_mm: f64,
    q_factor: f64,
) -> String {
    let params = TeslaCoilParams {
        resonant_freq_khz,
        input_kv,
        spark_gap_mm,
        q_factor,
    };

    // House canonical emission: the Franken-only dependency policy
    // (Decalogue P1) forbids a serde edge here, and one flat record needs
    // nothing smarter than explicit field formatting.
    let result = step_tesla_coil(&params);
    format!(
        "{{\"ok\":{{\"resonant_freq_khz\":{},\"secondary_potential_mv\":{},\
         \"streamer_length_inches\":{},\"streamer_length_meters\":{}}}}}",
        result.resonant_freq_khz,
        result.secondary_potential_mv,
        result.streamer_length_inches,
        result.streamer_length_meters,
    )
}
