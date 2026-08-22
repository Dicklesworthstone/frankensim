use fs_mbd::goddard::{GoddardParams, step_goddard_rocket};

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

    // House canonical emission: the Franken-only dependency policy
    // (Decalogue P1) forbids a serde edge here, and one flat record needs
    // nothing smarter than explicit field formatting.
    let result = step_goddard_rocket(&params);
    format!(
        "{{\"ok\":{{\"chamber_pressure_psi\":{},\"chamber_pressure_pa\":{},\
         \"exhaust_velocity_mps\":{},\"thrust_newtons\":{},\
         \"specific_impulse_sec\":{},\"mach_exit\":{}}}}}",
        result.chamber_pressure_psi,
        result.chamber_pressure_pa,
        result.exhaust_velocity_mps,
        result.thrust_newtons,
        result.specific_impulse_sec,
        result.mach_exit,
    )
}
