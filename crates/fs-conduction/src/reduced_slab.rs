//! First-mode one-dimensional slab cooling screen.
//!
//! This deliberately small reduced model evaluates
//! `T(t) = T_boundary + (T_initial - T_boundary) exp(-t/τ)` with
//! `τ = thickness²/(π² α)`. It admits only monotone cooling to a fixed,
//! uniform boundary temperature. It is not the transient finite-element rung,
//! does not infer convection coefficients or contact resistance, and does not
//! call a glass-transition crossing "solidification".

use core::f64::consts::PI;

/// SI inputs to the declared first-mode cooling screen.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FirstModeSlabInput {
    /// Slab thickness in metres.
    pub thickness_m: f64,
    /// Thermal diffusivity in square metres per second.
    pub thermal_diffusivity_m2_s: f64,
    /// Initial uniform temperature in kelvin.
    pub initial_temperature_k: f64,
    /// Fixed boundary temperature in kelvin.
    pub boundary_temperature_k: f64,
    /// Screening threshold in kelvin.
    pub threshold_temperature_k: f64,
}

/// Admitted first-mode cooling result.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FirstModeSlabStep {
    /// First-mode time constant in seconds.
    pub time_constant_s: f64,
    /// Time at which the exponential screen crosses the declared threshold.
    pub time_to_threshold_s: f64,
}

impl FirstModeSlabStep {
    /// Evaluate the declared exponential screen at non-negative time.
    pub fn temperature_at_s(
        self,
        input: FirstModeSlabInput,
        elapsed_s: f64,
    ) -> Result<f64, FirstModeSlabError> {
        if !elapsed_s.is_finite() {
            return Err(FirstModeSlabError::NonFinite { field: "elapsed_s" });
        }
        if elapsed_s < 0.0 {
            return Err(FirstModeSlabError::NegativeTime);
        }
        Ok(input.boundary_temperature_k
            + (input.initial_temperature_k - input.boundary_temperature_k)
                * (-elapsed_s / self.time_constant_s).exp())
    }
}

/// Typed refusal at the reduced thermal boundary.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FirstModeSlabError {
    /// A named input was NaN or infinite.
    NonFinite {
        /// Stable input field name.
        field: &'static str,
    },
    /// A named geometric or material quantity was zero or negative.
    NonPositive {
        /// Stable input field name.
        field: &'static str,
    },
    /// Temperatures do not describe monotone cooling through the threshold.
    ThresholdOutsideCoolingInterval,
    /// A requested evaluation time was negative.
    NegativeTime,
    /// The admitted calculation became non-finite.
    NonFiniteResult,
}

/// Evaluate the declared fixed-boundary first-mode cooling screen.
pub fn step_first_mode_slab_cooling(
    input: FirstModeSlabInput,
) -> Result<FirstModeSlabStep, FirstModeSlabError> {
    for (field, value) in [
        ("thickness_m", input.thickness_m),
        ("thermal_diffusivity_m2_s", input.thermal_diffusivity_m2_s),
        ("initial_temperature_k", input.initial_temperature_k),
        ("boundary_temperature_k", input.boundary_temperature_k),
        ("threshold_temperature_k", input.threshold_temperature_k),
    ] {
        if !value.is_finite() {
            return Err(FirstModeSlabError::NonFinite { field });
        }
    }
    for (field, value) in [
        ("thickness_m", input.thickness_m),
        ("thermal_diffusivity_m2_s", input.thermal_diffusivity_m2_s),
        ("initial_temperature_k", input.initial_temperature_k),
        ("boundary_temperature_k", input.boundary_temperature_k),
        ("threshold_temperature_k", input.threshold_temperature_k),
    ] {
        if value <= 0.0 {
            return Err(FirstModeSlabError::NonPositive { field });
        }
    }
    if !(input.initial_temperature_k > input.threshold_temperature_k
        && input.threshold_temperature_k > input.boundary_temperature_k)
    {
        return Err(FirstModeSlabError::ThresholdOutsideCoolingInterval);
    }

    let time_constant_s = input.thickness_m.powi(2) / (PI.powi(2) * input.thermal_diffusivity_m2_s);
    let temperature_ratio = (input.initial_temperature_k - input.boundary_temperature_k)
        / (input.threshold_temperature_k - input.boundary_temperature_k);
    let time_to_threshold_s = time_constant_s * temperature_ratio.ln();
    if !time_constant_s.is_finite() || !time_to_threshold_s.is_finite() {
        return Err(FirstModeSlabError::NonFiniteResult);
    }
    Ok(FirstModeSlabStep {
        time_constant_s,
        time_to_threshold_s,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const INPUT: FirstModeSlabInput = FirstModeSlabInput {
        thickness_m: 0.2e-3,
        thermal_diffusivity_m2_s: 0.082e-6,
        initial_temperature_k: 498.15,
        boundary_temperature_k: 298.15,
        threshold_temperature_k: 378.15,
    };

    #[test]
    fn crosses_the_threshold_at_the_reported_time() {
        let step = step_first_mode_slab_cooling(INPUT).expect("admitted cooling screen");
        let temperature = step
            .temperature_at_s(INPUT, step.time_to_threshold_s)
            .expect("non-negative crossing time");
        assert!((temperature - INPUT.threshold_temperature_k).abs() < 1.0e-10);
        assert!(step.time_to_threshold_s < step.time_constant_s);
    }

    #[test]
    fn thickness_has_the_expected_quadratic_scaling() {
        let base = step_first_mode_slab_cooling(INPUT).expect("base screen");
        let thick = step_first_mode_slab_cooling(FirstModeSlabInput {
            thickness_m: INPUT.thickness_m * 2.0,
            ..INPUT
        })
        .expect("thick screen");
        assert!((thick.time_constant_s / base.time_constant_s - 4.0).abs() < 1.0e-12);
        assert!((thick.time_to_threshold_s / base.time_to_threshold_s - 4.0).abs() < 1.0e-12);
    }

    #[test]
    fn refuses_unordered_thresholds_and_invalid_time() {
        assert_eq!(
            step_first_mode_slab_cooling(FirstModeSlabInput {
                threshold_temperature_k: INPUT.initial_temperature_k,
                ..INPUT
            }),
            Err(FirstModeSlabError::ThresholdOutsideCoolingInterval)
        );
        let step = step_first_mode_slab_cooling(INPUT).expect("base screen");
        assert_eq!(
            step.temperature_at_s(INPUT, -0.1),
            Err(FirstModeSlabError::NegativeTime)
        );
    }
}
