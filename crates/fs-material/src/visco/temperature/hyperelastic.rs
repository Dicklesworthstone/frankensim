//! Temperature-dependent kinetics in the original finite-strain return map.
//!
//! Thermorheological simplicity changes branch relaxation times, not the
//! hyperelastic energy, internal viscous metric or pressure constraint. The
//! original full and mixed/remainder solvers therefore keep their work and
//! heat accounting and their consistent algorithmic tangents.

use super::{ShiftedViscoelastic, TemperatureShift, scaled};
use crate::MaterialError;
use crate::hyper::Hyperelastic;
use crate::visco::ViscoError;
use crate::visco_hyper::ViscoHyperelastic;

/// Original finite-strain law at one declared temperature.
pub type ShiftedViscoHyperelastic = ShiftedViscoelastic<ViscoHyperelastic>;

/// Preserve the original calibration and finite-strain admission errors.
#[derive(Clone, Debug, PartialEq)]
pub enum HyperelasticTemperatureError {
    /// Temperature calibration or representability of a shifted time.
    Shift(ViscoError),
    /// Original equilibrium or branch parameters are inadmissible.
    Material(MaterialError),
}

impl core::fmt::Display for HyperelasticTemperatureError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Shift(error) => write!(f, "finite-strain temperature shift: {error}"),
            Self::Material(error) => write!(f, "finite-strain reference material: {error}"),
        }
    }
}

impl std::error::Error for HyperelasticTemperatureError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Shift(error) => Some(error),
            Self::Material(error) => Some(error),
        }
    }
}

impl TemperatureShift {
    /// Prepare the existing finite-strain law with `tau_j(T) = a_T tau_j(ref)`.
    ///
    /// The equilibrium energy, each Neo-Hookean/Ogden branch energy, branch
    /// order, deformation domain and original Newton controls are preserved.
    /// No eigensolve or return map runs during preparation. The returned law
    /// drives the same `step_response`, `freeze`, `freeze_remainder` and
    /// mixed finite-strain consumers.
    ///
    /// Reuse the original accepted viscous metrics when temperature changes;
    /// do not birth a fresh history. Always shift the ORIGINAL reference law,
    /// rather than shifting a previously shifted result. With unchanged
    /// energies, changing temperature alone creates no elastic work or heat.
    /// For one frozen temperature the original update is equivalent to the
    /// reference update over reduced duration `dt / a_T`.
    ///
    /// This is a constitutive kinetic shift, not a thermal expansion, caloric
    /// or entropy law. A varying-temperature trajectory requires explicit
    /// temporal refinement of its frozen-temperature intervals. Source support
    /// and qualification remain with the caller; this law has no source card.
    ///
    /// # Errors
    /// Refuses out-of-calibration temperatures, invalid original public
    /// parameters, or overflow/underflow of a positive relaxation time.
    /// The original law and every accepted history remain unchanged.
    pub fn shift_hyperelastic(
        &self,
        reference: &ViscoHyperelastic,
        temperature_k: f64,
    ) -> Result<ShiftedViscoHyperelastic, HyperelasticTemperatureError> {
        let state = self
            .evaluate(temperature_k)
            .map_err(HyperelasticTemperatureError::Shift)?;
        // Public reference fields can have changed since construction.
        // Reuse the original owners' parameter checks before scaling.
        let equilibrium = Hyperelastic::new(
            reference.equilibrium.model.clone(),
            reference.equilibrium.stretch_limit,
        )
        .map_err(HyperelasticTemperatureError::Material)?;
        let mut law = ViscoHyperelastic::new(equilibrium, reference.branches.clone())
            .and_then(|law| {
                law.with_local_newton(
                    reference.local_newton_tolerance(),
                    reference.local_newton_max_iterations(),
                )
            })
            .map_err(HyperelasticTemperatureError::Material)?;
        for branch in &mut law.branches {
            branch.tau_s = scaled(branch.tau_s, state.factor())
                .map_err(HyperelasticTemperatureError::Shift)?;
        }
        Ok(ShiftedViscoelastic { law, state })
    }
}
