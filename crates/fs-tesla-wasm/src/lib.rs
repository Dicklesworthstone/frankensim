use fs_flux::lc::{step_tesla_coil, TeslaCoilParams};
use serde_json::json;

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
    
    let result = step_tesla_coil(&params);
    json!({ "ok": result }).to_string()
}
