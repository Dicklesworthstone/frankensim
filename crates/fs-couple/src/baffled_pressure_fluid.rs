//! Original Rayleigh surface paths with generic-fluid classical bulk absorption.
//! Every signed quadrature path retains its own distance, attenuation and delay.
use super::gas::PathFilter;
use super::{BaffledPressure, RayleighMedium, SurfaceSample};
use crate::air_path::fluid::{FluidObserverFilter, FluidPropagationReceipt};
use fs_material::fluid_acoustics::thermoviscous::absorption::ClassicalBulkAbsorption;

/// The existing material-path work limits, shared without a second budget law.
pub use super::gas::GasPropagationLimits as FluidPropagationLimits;

/// Original material coefficients and independent numerical/work controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FluidPropagation {
    /// Original fluid snapshot and explicit constant bulk-viscosity law.
    pub absorption: ClassicalBulkAbsorption,
    /// Relative impulse-response tail energy ceiling for each physical path.
    pub maximum_truncated_energy: f64,
    /// Cold preparation bounds shared with the established gas-path owner.
    pub limits: FluidPropagationLimits,
}

/// Per-path material transfers, before the shared sinc and causal padding.
#[derive(Clone, Debug)]
pub struct FluidPropagationReport {
    /// Index of each original surface sample's bit-exact distance response.
    pub surface_filter_indices: Vec<usize>,
    /// One original response report per distinct exact path distance.
    pub filters: Vec<FluidPropagationReceipt>,
    /// Common whole-sample padding, preserving all relative physical delays.
    pub additional_latency_samples: usize,
    /// Actual cold convolution and signed-projection scalar products.
    pub kernel_terms: usize,
}

impl BaffledPressure {
    /// Observe a flat vibrating surface through its original stationary fluid.
    /// Classical absorption is applied on each exact surface-to-receiver path
    /// before summing signed source contributions. It is not a filter on the
    /// already mixed pressure, and it introduces no second propagation delay.
    ///
    /// Reuses the complete gas-path geometry, cold convolution, cancellation,
    /// work limits, sinc interpolation and allocation-free stepping owner.
    /// The actual fluid supplies `rho*dA/(2*pi*r)` and `r/c`; no gas EOS or
    /// atmospheric coefficients are substituted. Preparation is atomic.
    ///
    /// This remains prescribed linear motion in an infinite rigid baffle,
    /// without backreaction, free-surface scattering or phase motion. The
    /// absorption is a weak-loss, in-band minimum-phase approximation; its
    /// model/source limits remain on each original filter receipt.
    #[allow(clippy::too_many_arguments)]
    pub fn new_multirate_fluid_differentiating(
        surface: &[SurfaceSample],
        inputs: usize,
        mechanics_rate: u32,
        output_rate: u32,
        position_m: [f64; 3],
        propagation: &FluidPropagation,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Self, FluidPropagationReport), String> {
        let acoustic = propagation.absorption.state().acoustic();
        let medium = RayleighMedium {
            density: acoustic.density_kg_m3(),
            sound_speed: acoustic.sound_speed_m_s(),
        };
        let (receiver, paths) = Self::new_multirate_material_differentiating(
            surface,
            inputs,
            mechanics_rate,
            output_rate,
            position_m,
            medium,
            propagation.limits,
            cancelled,
            |distance| {
                let filter = FluidObserverFilter::new(
                    propagation.absorption,
                    distance,
                    1.0 / f64::from(output_rate),
                    propagation.maximum_truncated_energy,
                )
                .map_err(|error| format!("Rayleigh fluid path at {distance} m: {error}"))?;
                Ok(PathFilter {
                    response: filter.taps().to_vec(),
                    receipt: *filter.receipt(),
                    receiver_delay_s: filter.receipt().receiver_delay_s,
                    low_frequency_delay_s: filter.receipt().low_frequency_delay_s,
                })
            },
        )?;
        Ok((
            receiver,
            FluidPropagationReport {
                surface_filter_indices: paths.surface_filter_indices,
                filters: paths.filters,
                additional_latency_samples: paths.additional_latency_samples,
                kernel_terms: paths.kernel_terms,
            },
        ))
    }
}

#[cfg(test)]
#[path = "baffled_pressure_fluid_tests.rs"]
mod tests;
