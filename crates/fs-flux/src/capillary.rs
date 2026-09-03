//! Reduced, deterministic circular-capillary flow screening.
//!
//! This module owns the Hagen–Poiseuille relation for a declared Newtonian,
//! incompressible fluid in a straight, rigid, circular capillary with
//! fully-developed laminar flow and no slip. It deliberately does not turn
//! that screening law into a polymer-rheology claim: entrance contraction,
//! shear thinning, viscoelasticity, wall slip, heating, and free-surface flow
//! are outside this rung.

use core::f64::consts::PI;

/// SI inputs admitted by the circular-capillary screening law.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CircularCapillaryInput {
    /// Dynamic viscosity in pascal-seconds.
    pub dynamic_viscosity_pa_s: f64,
    /// Straight capillary land length in metres.
    pub length_m: f64,
    /// Circular capillary radius in metres.
    pub radius_m: f64,
    /// Volumetric flow rate in cubic metres per second.
    pub volumetric_flow_m3_s: f64,
}

/// Deterministic output of the admitted reduced law.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CircularCapillaryStep {
    /// Required pressure drop in pascals.
    pub pressure_drop_pa: f64,
    /// Newtonian wall shear rate, `4Q/(πR³)`, in reciprocal seconds.
    pub wall_shear_rate_per_s: f64,
    /// Hydraulic power `ΔP Q` in watts.
    pub hydraulic_power_w: f64,
}

/// Typed refusal at the reduced-law boundary.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CircularCapillaryError {
    /// A named input was NaN or infinite.
    NonFinite {
        /// Stable input field name.
        field: &'static str,
    },
    /// A named input that must be strictly positive was zero or negative.
    NonPositive {
        /// Stable input field name.
        field: &'static str,
    },
    /// Reverse flow is outside this scalar browser-facing screening rung.
    NegativeFlow,
    /// The admitted calculation overflowed or otherwise became non-finite.
    NonFiniteResult,
}

/// Evaluate Hagen–Poiseuille flow within the declared reduced-law boundary.
///
/// Zero flow is admitted and returns zero pressure, shear rate, and power.
pub fn step_newtonian_circular_capillary(
    input: CircularCapillaryInput,
) -> Result<CircularCapillaryStep, CircularCapillaryError> {
    for (field, value) in [
        ("dynamic_viscosity_pa_s", input.dynamic_viscosity_pa_s),
        ("length_m", input.length_m),
        ("radius_m", input.radius_m),
        ("volumetric_flow_m3_s", input.volumetric_flow_m3_s),
    ] {
        if !value.is_finite() {
            return Err(CircularCapillaryError::NonFinite { field });
        }
    }
    for (field, value) in [
        ("dynamic_viscosity_pa_s", input.dynamic_viscosity_pa_s),
        ("length_m", input.length_m),
        ("radius_m", input.radius_m),
    ] {
        if value <= 0.0 {
            return Err(CircularCapillaryError::NonPositive { field });
        }
    }
    if input.volumetric_flow_m3_s < 0.0 {
        return Err(CircularCapillaryError::NegativeFlow);
    }

    let radius_cubed = input.radius_m.powi(3);
    let pressure_drop_pa =
        8.0 * input.dynamic_viscosity_pa_s * input.length_m * input.volumetric_flow_m3_s
            / (PI * input.radius_m.powi(4));
    let wall_shear_rate_per_s = 4.0 * input.volumetric_flow_m3_s / (PI * radius_cubed);
    let hydraulic_power_w = pressure_drop_pa * input.volumetric_flow_m3_s;

    if !pressure_drop_pa.is_finite()
        || !wall_shear_rate_per_s.is_finite()
        || !hydraulic_power_w.is_finite()
    {
        return Err(CircularCapillaryError::NonFiniteResult);
    }

    Ok(CircularCapillaryStep {
        pressure_drop_pa,
        wall_shear_rate_per_s,
        hydraulic_power_w,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_the_closed_form_and_quartic_radius_scaling() {
        let base = step_newtonian_circular_capillary(CircularCapillaryInput {
            dynamic_viscosity_pa_s: 280.0,
            length_m: 1.6e-3,
            radius_m: 0.2e-3,
            volumetric_flow_m3_s: 4.05e-9,
        })
        .expect("admitted capillary");
        let wide = step_newtonian_circular_capillary(CircularCapillaryInput {
            radius_m: 0.4e-3,
            ..CircularCapillaryInput {
                dynamic_viscosity_pa_s: 280.0,
                length_m: 1.6e-3,
                radius_m: 0.2e-3,
                volumetric_flow_m3_s: 4.05e-9,
            }
        })
        .expect("wider admitted capillary");

        assert!((base.pressure_drop_pa - 2_887_707.287_459_349).abs() < 1.0e-6);
        assert!((base.wall_shear_rate_per_s - 644.577_519_522_176).abs() < 1.0e-9);
        assert!((base.pressure_drop_pa / wide.pressure_drop_pa - 16.0).abs() < 1.0e-12);
        assert!((base.hydraulic_power_w - base.pressure_drop_pa * 4.05e-9).abs() < 1.0e-15);
    }

    #[test]
    fn zero_flow_is_a_valid_stationary_state() {
        let result = step_newtonian_circular_capillary(CircularCapillaryInput {
            dynamic_viscosity_pa_s: 1.0,
            length_m: 1.0,
            radius_m: 0.1,
            volumetric_flow_m3_s: 0.0,
        })
        .expect("zero-flow state");
        assert_eq!(result.pressure_drop_pa, 0.0);
        assert_eq!(result.wall_shear_rate_per_s, 0.0);
        assert_eq!(result.hydraulic_power_w, 0.0);
    }

    #[test]
    fn refuses_invalid_physical_inputs() {
        let input = CircularCapillaryInput {
            dynamic_viscosity_pa_s: 1.0,
            length_m: 1.0,
            radius_m: 0.1,
            volumetric_flow_m3_s: 1.0,
        };
        assert_eq!(
            step_newtonian_circular_capillary(CircularCapillaryInput {
                radius_m: 0.0,
                ..input
            }),
            Err(CircularCapillaryError::NonPositive { field: "radius_m" })
        );
        assert_eq!(
            step_newtonian_circular_capillary(CircularCapillaryInput {
                volumetric_flow_m3_s: -1.0,
                ..input
            }),
            Err(CircularCapillaryError::NegativeFlow)
        );
        assert_eq!(
            step_newtonian_circular_capillary(CircularCapillaryInput {
                length_m: f64::NAN,
                ..input
            }),
            Err(CircularCapillaryError::NonFinite { field: "length_m" })
        );
    }
}
