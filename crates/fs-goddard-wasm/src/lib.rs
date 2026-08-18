use fs_mbd::goddard::{step_goddard_rocket, GoddardParams};
use serde_json::json;

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

#[cfg_attr(target_arch = "wasm32", wasm_bindgen)]
pub fn goddard_rocket_step(
    chamber_pressure_psi: f64,
    fuel_flow_kg_per_sec: f64,
    throat_area_cm2: f64,
    expansion_ratio: f64,
) -> String {
    let params = GoddardParams {
        chamber_pressure_psi,
        fuel_flow_kg_per_sec,
        throat_area_cm2,
        expansion_ratio,
    };
    
    let result = step_goddard_rocket(&params);
    json!({ "ok": result }).to_string()
}
