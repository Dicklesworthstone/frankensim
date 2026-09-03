//! Browser boundary for a source-bounded US 5,121,329 physical screen.
//!
//! Generic laws remain owned by `fs-flux` and `fs-conduction`. This L6 crate
//! admits browser scalars and serializes one deterministic receipt.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use core::fmt::Write as _;
use fs_conduction::{FirstModeSlabInput, step_first_mode_slab_cooling};
use fs_flux::{CircularCapillaryInput, step_newtonian_circular_capillary};

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

fn refusal_json(code: &str, message: &str, repair: &str) -> String {
    format!(
        "{{\"refusal\":{{\"code\":\"{code}\",\"message\":\"{message}\",\"repairs\":[\"{repair}\"]}}}}"
    )
}

/// Compose generic capillary-flow and reduced slab-cooling owners.
///
/// Every argument is SI: Pa·s, m, m³/s, m²/s, and K. The response explicitly
/// carries reduced-model boundaries; it is not a historic performance claim.
#[cfg_attr(target_arch = "wasm32", wasm_bindgen)]
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn crump_fdm_step(
    dynamic_viscosity_pa_s: f64,
    capillary_length_m: f64,
    capillary_radius_m: f64,
    volumetric_flow_m3_s: f64,
    layer_thickness_m: f64,
    thermal_diffusivity_m2_s: f64,
    initial_temperature_k: f64,
    boundary_temperature_k: f64,
    threshold_temperature_k: f64,
) -> String {
    let capillary = match step_newtonian_circular_capillary(CircularCapillaryInput {
        dynamic_viscosity_pa_s,
        length_m: capillary_length_m,
        radius_m: capillary_radius_m,
        volumetric_flow_m3_s,
    }) {
        Ok(step) => step,
        Err(_) => {
            return refusal_json(
                "capillary-input-outside-domain",
                "the generic circular-capillary owner refused a non-finite, non-positive, reverse-flow, or overflowing input",
                "Use finite positive viscosity and geometry with a finite non-negative SI flow rate",
            );
        }
    };

    let thermal_input = FirstModeSlabInput {
        thickness_m: layer_thickness_m,
        thermal_diffusivity_m2_s,
        initial_temperature_k,
        boundary_temperature_k,
        threshold_temperature_k,
    };
    let thermal = match step_first_mode_slab_cooling(thermal_input) {
        Ok(step) => step,
        Err(_) => {
            return refusal_json(
                "thermal-screen-input-outside-domain",
                "the generic first-mode slab owner requires finite positive SI inputs with initial temperature above threshold above boundary",
                "Supply a monotone cooling interval and positive layer thickness and diffusivity",
            );
        }
    };
    let threshold_temperature_check_k = match thermal
        .temperature_at_s(thermal_input, thermal.time_to_threshold_s)
    {
        Ok(value) => value,
        Err(_) => {
            return refusal_json(
                "thermal-screen-evaluation-failed",
                "the admitted thermal-screen result could not be evaluated at its threshold time",
                "Inspect the fs-conduction reduced-slab contract",
            );
        }
    };

    let mut output = String::with_capacity(720);
    let _ = write!(
        output,
        "{{\"ok\":{{\
         \"capillary_owner\":\"fs-flux::capillary::step_newtonian_circular_capillary\",\
         \"thermal_owner\":\"fs-conduction::reduced_slab::step_first_mode_slab_cooling\",\
         \"pressure_drop_pa\":{},\
         \"wall_shear_rate_per_s\":{},\
         \"hydraulic_power_w\":{},\
         \"cooling_time_constant_s\":{},\
         \"time_to_threshold_s\":{},\
         \"threshold_temperature_check_k\":{},\
         \"capillary_boundary\":\"newtonian-incompressible-fully-developed-laminar-no-slip-circular-land\",\
         \"thermal_boundary\":\"one-dimensional-fixed-boundary-first-mode-screen-no-phase-change\"\
        }}}}",
        capillary.pressure_drop_pa,
        capillary.wall_shear_rate_per_s,
        capillary.hydraulic_power_w,
        thermal.time_constant_s,
        thermal.time_to_threshold_s,
        threshold_temperature_check_k,
    );
    output
}

#[cfg(test)]
mod tests {
    use super::*;

    fn admitted_step(flow_m3_s: f64) -> String {
        crump_fdm_step(
            280.0, 1.6e-3, 0.2e-3, flow_m3_s, 0.2e-3, 0.082e-6, 498.15, 298.15, 378.15,
        )
    }

    #[test]
    fn serializes_both_generic_owners_and_si_outputs() {
        let output = admitted_step(4.05e-9);
        assert!(output.contains("\"capillary_owner\":\"fs-flux::capillary"));
        assert!(output.contains("\"thermal_owner\":\"fs-conduction::reduced_slab"));
        assert!(output.contains("\"pressure_drop_pa\":2887707.287459"));
        assert!(output.contains("\"threshold_temperature_check_k\":378.15"));
        assert!(!output.contains("\"refusal\""));
    }

    #[test]
    fn admits_a_stationary_zero_flow_screen() {
        let output = admitted_step(0.0);
        assert!(output.contains("\"pressure_drop_pa\":0"));
        assert!(output.contains("\"hydraulic_power_w\":0"));
    }

    #[test]
    fn emits_typed_refusals_instead_of_partial_results() {
        let bad_radius = crump_fdm_step(
            280.0, 1.6e-3, 0.0, 4.05e-9, 0.2e-3, 0.082e-6, 498.15, 298.15, 378.15,
        );
        assert!(bad_radius.contains("\"code\":\"capillary-input-outside-domain\""));
        assert!(!bad_radius.contains("\"ok\""));

        let bad_threshold = crump_fdm_step(
            280.0, 1.6e-3, 0.2e-3, 4.05e-9, 0.2e-3, 0.082e-6, 498.15, 298.15, 520.0,
        );
        assert!(bad_threshold.contains("\"code\":\"thermal-screen-input-outside-domain\""));
        assert!(!bad_threshold.contains("\"ok\""));
    }
}
