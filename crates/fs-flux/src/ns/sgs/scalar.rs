//! Scalar heat/species mixing from the SAME accepted planar eddy viscosity.

use super::{PlanarSgs, PlanarSgsError, PlanarSgsField};
use crate::ns::transport::scalar::{
    PreparedScalarTransport, ScalarBoundary, ScalarTransportControls, ScalarTransportError,
};

#[derive(Clone, Debug, PartialEq)]
pub enum PlanarSgsScalarError {
    Subgrid(PlanarSgsError),
    Scalar(ScalarTransportError),
    InvalidCell { cell: usize },
}
impl core::fmt::Display for PlanarSgsScalarError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "planar SGS scalar: {self:?}")
    }
}
impl std::error::Error for PlanarSgsScalarError {}

/// One velocity endpoint supplies BOTH original momentum closure diagnostics
/// and scalar advection/diffusion. No second velocity reconstruction is used.
pub struct PlanarSgsScalarTransport<'m> {
    pub transport: PreparedScalarTransport<'m>,
    pub sgs: PlanarSgsField,
    pub turbulent_diffusivity_m2_s: Vec<f64>,
}

impl<'s, 'm> PlanarSgs<'s, 'm> {
    /// Prepare actual passive thermal/species evolution with
    /// `alpha_eff = alpha_molecular + nu_t / ratio`. Supply positive finite
    /// turbulent Prandtl numbers for temperature or Schmidt numbers for a
    /// passive species; neither coefficient nor molecular diffusivity is guessed.
    /// All three original SGS models retain their authored widths and values.
    ///
    /// The accepted velocity may come from any native steady/Picard/Newton/BDF
    /// flow path. Scalar BDF1/BDF2 then uses this frozen endpoint velocity and
    /// diffusivity on the same mesh. Reprepare at each new flow endpoint; no
    /// feedback into molecular viscosity, buoyancy or momentum is implied.
    /// Molecular and modeled momentum dissipation are NOT added as heat.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_scalar_transport(
        &self,
        velocity: &[f64],
        molecular_diffusivity_m2_s: &[f64],
        turbulent_ratio: &[f64],
        boundaries: &[ScalarBoundary<'_>],
        controls: ScalarTransportControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PlanarSgsScalarTransport<'m>, PlanarSgsScalarError> {
        if cancelled() {
            return Err(PlanarSgsScalarError::Scalar(
                ScalarTransportError::Cancelled,
            ));
        }
        let n = self.system.mesh.tris.len();
        if molecular_diffusivity_m2_s.len() != n || turbulent_ratio.len() != n {
            return Err(PlanarSgsScalarError::Scalar(ScalarTransportError::Invalid(
                "one molecular diffusivity and turbulent ratio per original cell required",
            )));
        }
        for cell in 0..n {
            if !molecular_diffusivity_m2_s[cell].is_finite()
                || molecular_diffusivity_m2_s[cell] < 0.0
                || !turbulent_ratio[cell].is_finite()
                || turbulent_ratio[cell] <= 0.0
            {
                return Err(PlanarSgsScalarError::InvalidCell { cell });
            }
        }
        let sgs = self
            .field(velocity)
            .map_err(PlanarSgsScalarError::Subgrid)?;
        let mut turbulent = Vec::with_capacity(n);
        let mut effective = Vec::with_capacity(n);
        for (cell, diagnostic) in sgs.cells.iter().enumerate() {
            if cancelled() {
                return Err(PlanarSgsScalarError::Scalar(
                    ScalarTransportError::Cancelled,
                ));
            }
            let modeled = diagnostic.eddy_viscosity_m2_s / turbulent_ratio[cell];
            let total = molecular_diffusivity_m2_s[cell] + modeled;
            if !modeled.is_finite()
                || !total.is_finite()
                || (diagnostic.eddy_viscosity_m2_s > 0.0 && modeled == 0.0)
            {
                return Err(PlanarSgsScalarError::InvalidCell { cell });
            }
            turbulent.push(modeled);
            effective.push(total);
        }
        let transport = self
            .system
            .prepare_scalar_transport(velocity, &effective, boundaries, controls, cancelled)
            .map_err(PlanarSgsScalarError::Scalar)?;
        Ok(PlanarSgsScalarTransport {
            transport,
            sgs,
            turbulent_diffusivity_m2_s: turbulent,
        })
    }
}
