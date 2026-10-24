//! Bind a physical baffled exterior to the original DG endpoint.
//!
//! This reuses the existing geometry/gas BEM sampling, rational bake and common
//! load admission. A circular outlet at the origin of an infinite rigid baffle
//! is explicit source physics, not a free-field bell/body approximation.

use super::{DgBoundary, DgDuct, DgPortLoad, invalid};
use crate::bernoulli_aperture::radiation::environment::multiport::input::MimoRadiationInputError;
use crate::bernoulli_aperture::radiation::environment::multiport::source::ComputedBaffledInput;
use crate::pcm_wav::baffled::openings::CircularOpening;
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
