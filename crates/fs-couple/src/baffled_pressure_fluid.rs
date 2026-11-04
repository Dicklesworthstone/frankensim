//! Original Rayleigh surface paths with selected generic-fluid bulk propagation.
//! Every signed quadrature path retains its own distance, attenuation and delay.
use super::gas::PathFilter;
use super::{BaffledPressure, RayleighMedium, SurfaceSample};
use crate::air_path::fluid::{
    FluidObserverFilter, FluidPropagationModel, FluidPropagationReceipt,
};
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

impl FluidPropagation {
    pub(super) fn rayleigh_medium(&self) -> RayleighMedium {
        let acoustic = self.absorption.state().acoustic();
        RayleighMedium {
            density: acoustic.density_kg_m3(),
            sound_speed: acoustic.sound_speed_m_s(),
        }
    }
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
    /// Align independently prepared fluid receivers before accepting any trace.
    /// Use the greatest `additional_latency_samples` in their original reports
    /// to sum pressures on a common numerical clock. Physical travel, source
    /// attenuation, relative phase, quadrature and filter receipts stay intact.
    /// The added padding appears once in [`Self::filter_delay_samples`].
    ///
    /// Receivers must use the same output clock and derivative reconstruction.
    /// Supply this receiver's own unchanged report. A stale padding value,
    /// decreasing latency, exhausted history budget, or any previously accepted
    /// trace refuses before changing the receiver or report. Even a silent
    /// accepted trace advances the clock and closes this cold operation.
    pub fn align_fluid_latency(
        &mut self,
        report: &mut FluidPropagationReport,
        latency: usize,
    ) -> Result<(), String> {
        self.align_material_latency(&mut report.additional_latency_samples, latency)
    }

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
        Self::new_multirate_fluid_differentiating_with_model(
            surface,
            inputs,
            mechanics_rate,
            output_rate,
            position_m,
            propagation,
            FluidPropagationModel::ClassicalMinimumPhase,
            cancelled,
        )
    }

    /// Select classical minimum phase or the original fluid's thermoviscous
    /// acoustic dispersion on every physical surface-to-receiver path.
    /// Complex phase, exact-distance caching, signed projection and one shared
    /// numerical padding retain the same preparation/stepping owner.
    ///
    /// The selected path law remains within the observer's supported weak-loss
    /// band. It does not replace the Rayleigh source normalization or impose
    /// thermoviscous boundary conditions on the vibrating surface. The `r/c0`
    /// reference and reported numerical latency do not assert a finite front
    /// for the parabolic Navier–Stokes–Fourier model.
    #[allow(clippy::too_many_arguments)]
    pub fn new_multirate_fluid_differentiating_with_model(
        surface: &[SurfaceSample],
        inputs: usize,
        mechanics_rate: u32,
        output_rate: u32,
        position_m: [f64; 3],
        propagation: &FluidPropagation,
        model: FluidPropagationModel,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Self, FluidPropagationReport), String> {
        let medium = propagation.rayleigh_medium();
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
                let filter = FluidObserverFilter::with_model(
                    propagation.absorption,
                    distance,
                    1.0 / f64::from(output_rate),
                    propagation.maximum_truncated_energy,
                    model,
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
