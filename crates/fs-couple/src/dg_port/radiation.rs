//! Bind a physical baffled exterior to the original DG endpoint.
//!
//! This reuses the existing geometry/gas BEM sampling, rational bake and common
//! load admission. A circular outlet at the origin of an infinite rigid baffle
//! is explicit source physics, not a free-field bell/body approximation.

use super::{DgBoundary, DgDuct, DgPortLoad, invalid};
use crate::bernoulli_aperture::radiation::environment::multiport::input::MimoRadiationInputError;
use crate::bernoulli_aperture::radiation::environment::multiport::source::ComputedBaffledInput;
use crate::pcm_wav::baffled::openings::CircularOpening;
use crate::pcm_wav::baffled::{BaffledPressure, CircularOutletReceiver};
use crate::acoustic_realize::AcousticRealizeError;
use fs_exec::CancelGate;
use fs_material::gas_species::input::GasRecord;
use fs_math::det;
use fs_vfit::mimo::load::{MimoLoadLimits, PreparedMimoLoad};

impl DgPortLoad {
    /// Compute a causal radiation load from the actual outlet area and explicit
    /// exterior material, preserving the complete source record. This is cold
    /// preparation; a material substitution recomputes source and runtime once.
    ///
    /// The bore and exterior have one equilibrium static pressure. Different
    /// compositions/temperatures are allowed, but a pressure jump would require
    /// a mean-flow/interface model and is refused by this quiescent adapter.
    ///
    /// The original source bake must resolve the full mechanical Nyquist band.
    /// Its numerical fit and GridOnly checks retain their actual report classes.
    /// Selecting a microphone later does not recompute or alter this reaction.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_computed_baffled(
        duct: DgDuct,
        mechanics_rate_hz: u32,
        exponential_tolerance: f64,
        maximum_dofs: usize,
        exterior: &GasRecord,
        source: ComputedBaffledInput,
        gate: &CancelGate,
    ) -> Result<Self, MimoRadiationInputError> {
        if gate.is_requested() { return Err(MimoRadiationInputError::Cancelled); }
        let dt = 1.0 / f64::from(mechanics_rate_hz);
        if mechanics_rate_hz == 0
            || !matches!(duct.inlet(), DgBoundary::Port(_))
            || !matches!(duct.outlet(), DgBoundary::Port(_))
            || duct.total_dofs() > maximum_dofs
            || duct.viscothermal_band_rad_s().is_some_and(|band|
                band[1] * dt > core::f64::consts::PI)
            || exterior.state().pressure.to_bits() != duct.outlet_gas().pressure.to_bits()
        {
            return Err(MimoRadiationInputError::Physics(invalid(
                "computed DG radiation requires owned ports, a resolved clock/state budget and equal equilibrium pressures",
            )));
        }
        let area = duct.outlet_area_m2();
        let radius_m = det::sqrt(area / core::f64::consts::PI);
        if !radius_m.is_finite() || radius_m <= 0.0 {
            return Err(MimoRadiationInputError::Physics(invalid(
                "computed DG radiation requires a finite positive original outlet",
            )));
        }
        let (model, record) = source.build_apertures(
            exterior, &[CircularOpening { center_m: [0.0; 2], radius_m }],
            mechanics_rate_hz, gate,
        )?;
        let limits = source.limits();
        // The band endpoint is the same physical Nyquist as the source's Hz
        // declaration, evaluated in the runtime's own clock arithmetic.
        let terminal = PreparedMimoLoad::new(&model.impedance, dt, MimoLoadLimits {
            band_rad_s: [
                core::f64::consts::TAU * model.band_hz[0],
                core::f64::consts::PI / dt,
            ],
            sampled_passivity_points: limits.sampled_passivity_points,
            max_states: limits.max_states,
            max_memory_bytes: limits.max_memory_bytes,
        }).map_err(|error| MimoRadiationInputError::Physics(invalid(error.0)))?;
        if gate.is_requested() { return Err(MimoRadiationInputError::Cancelled); }
        let mut load = Self::new_with_terminal(
            duct, dt, exponential_tolerance, maximum_dofs, terminal,
        ).map_err(MimoRadiationInputError::Physics)?;
        load.radiation_record = Some(record);
        if gate.is_requested() { return Err(MimoRadiationInputError::Cancelled); }
        Ok(load)
    }
}

impl DgPortLoad {
    /// Prepare an observer of this computed exterior's accepted terminal flow.
    /// Source geometry, full gas state, clock and radiation/wall bands remain
    /// owned by this load; a microphone never substitutes or re-bakes them.
    /// The original Rayleigh owner admits quadrature, delay/history and its
    /// same-clock differentiating image. This adds no acoustic reaction.
    pub fn baffled_receiver(
        &self,
        mechanics_rate_hz: u32,
        receiver: CircularOutletReceiver,
    ) -> Result<BaffledPressure, AcousticRealizeError> {
        if mechanics_rate_hz == 0 || self.state().step_count() != 0
            || self.propagator().dt().to_bits()
                != (1.0 / f64::from(mechanics_rate_hz)).to_bits()
        {
            return Err(invalid("DG receiver requires a cold source and its original physical clock"));
        }
        let source = self.radiation_record().ok_or_else(|| invalid(
            "baffled DG observation requires the original computed radiation source",
        ))?;
        let [aperture] = source.apertures() else {
            return Err(invalid("baffled DG observation requires one bound circular outlet"));
        };
        let area = core::f64::consts::PI * aperture.radius_m * aperture.radius_m;
        let actual_area = self.duct().outlet_area_m2();
        if aperture.center_m != [0.0; 2]
            || !area.is_finite()
            || (area - actual_area).abs() > 64.0 * f64::EPSILON * actual_area
            || source.band_hz()[1].to_bits()
                != (0.5 * f64::from(mechanics_rate_hz)).to_bits()
        {
            return Err(invalid(
                "baffled DG observation must use its original outlet geometry and physical clock",
            ));
        }
        // Compare in the wall owner's native angular-frequency arithmetic.
        // Converting its endpoint back to Hz can reject the identical endpoint
        // solely because multiplication and division round differently.
        if receiver.maximum_frequency_hz > source.band_hz()[1]
            || self.duct().viscothermal_band_rad_s().is_some_and(|band|
                receiver.maximum_frequency_hz * core::f64::consts::TAU > band[1])
        {
            return Err(invalid("baffled DG receiver exceeds its radiation or wall-loss band"));
        }
        BaffledPressure::circular_openings_multirate_differentiating(
            source.apertures(), mechanics_rate_hz, mechanics_rate_hz,
            receiver, source.gas_record().state().into(),
        ).map_err(AcousticRealizeError::Nonlinear)
    }
}
