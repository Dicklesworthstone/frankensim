//! Spatial transport coefficients on the existing source-coupled thermal slabs.
//!
//! These are numerical relaxation times, not another material law. Scalars
//! per cell describe isotropic transport with uniform thermal storage. The
//! existing BGK, low-Mach and boundary approximations remain in force.

use super::{
    ThermalLbm, ThermalLbm3, ThermalSourceBalance, ThermalSourceError,
    ThermalSourceLimits, invalid, source_step_transport_mapped,
};

/// Complete, step-local momentum and thermal relaxation fields.
///
/// Supply one value per cell, including both wall planes. Values must be finite
/// and representably greater than 1/2. No floor or cap changes the material law.
/// Flow uses nu=(tau-1/2)/3 in both dimensions. Heat uses alpha=(tau-1/2)/3
/// for D2Q5, but alpha=(tau-1/2)/4 for D3Q7; equal physical diffusivity does NOT
/// generally mean equal thermal relaxation times.
///
/// The fields are held fixed over this one step. They do not reset populations
/// or persist into subsequent calls. Existing scalar diagnostics such as
/// `diffusivity`, `rayleigh` and `nusselt` describe the background setup, not
/// the spatial field; retain the supplied coefficients for local diagnostics.
#[derive(Debug, Clone, Copy)]
pub struct RelaxationField<'a> {
    /// Local momentum relaxation time (dimensionless lattice steps).
    pub momentum: &'a [f64],
    /// Local thermal relaxation time (dimensionless lattice steps).
    pub thermal: &'a [f64],
}
impl RelaxationField<'_> {
    pub(super) fn admit(self, cells: usize) -> Result<(), ThermalSourceError> {
        if self.momentum.len() != cells || self.thermal.len() != cells {
            return Err(invalid("transport requires one momentum and thermal relaxation time per cell"));
        }
        if self.momentum.iter().chain(self.thermal).any(|t| !t.is_finite() || *t <= 0.5) {
            return Err(invalid("local transport relaxation must be finite and greater than one half"));
        }
        Ok(())
    }
}

impl ThermalLbm {
    /// Advance the original D2Q9/D2Q5 solver under local isotropic transport,
    /// prescribed heat increments and external lattice force densities.
    /// The same independent wall/source balance and transactional checks as
    /// `step_with_sources` apply. Coefficients are borrowed for this step only;
    /// no viscosity remains behind when a later call omits the thermal field.
    /// This is not variable-density/heat-capacity storage or a stability proof.
    pub fn step_with_sources_and_transport(
        &mut self, heat_per_step: &[f64], forces: &[[f64; 2]],
        transport: RelaxationField<'_>, limits: ThermalSourceLimits,
    ) -> Result<ThermalSourceBalance, ThermalSourceError> {
        source_step_transport_mapped::<Self, 2, 5, _>(
            self, heat_per_step, forces, limits, Some(transport), |_, balance| Ok(balance),
        )
    }
}
impl ThermalLbm3 {
    /// Three-dimensional counterpart, using the existing D3Q19/D3Q7 kernels.
    /// Local rates, source/force semantics, balance and rollback are shared.
    pub fn step_with_sources_and_transport(
        &mut self, heat_per_step: &[f64], forces: &[[f64; 3]],
        transport: RelaxationField<'_>, limits: ThermalSourceLimits,
    ) -> Result<ThermalSourceBalance, ThermalSourceError> {
        source_step_transport_mapped::<Self, 3, 7, _>(
            self, heat_per_step, forces, limits, Some(transport), |_, balance| Ok(balance),
        )
    }
}
