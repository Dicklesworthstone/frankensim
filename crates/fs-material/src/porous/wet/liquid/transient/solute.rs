//! Passive neutral solutes carried by the original porous liquid interval.
//!
//! This is dilute Fick transport in the stationary solid frame:
//! J_k=q*c_k-(theta*D_k/tau)*grad(c_k), theta=phi*S_l.
//! Each solute has its explicit aqueous diffusivity; tau is the transport divisor
//! in this formula. Geometric or acoustic tortuosity conventions require an
//! explicit source-supported conversion. The same caller
//! ordering is required in every constituent array. The caller qualifies D/tau
//! and the passive dilute approximation at the actual aqueous state. Solutes do not alter water
//! pressure, saturation, viscosity, density or retention. Activity/cross-friction,
//! ionic migration, reaction, adsorption and hydrodynamic dispersion need their
//! own constitutive owners and are not inferred from these passive coefficients.

use super::{LiquidSaturationCandidate, LiquidSaturationColumn};

mod transport;
pub use transport::{
    DiluteSoluteBoundary, DiluteSoluteControls, DiluteSoluteError, DiluteSoluteReport,
    DiluteSoluteStep,
};

impl LiquidSaturationCandidate {
    /// Advance passive dissolved inventories over this SAME accepted water step.
    ///
    /// The candidate binds its exact original saturation/time and native cells/
    /// calibrations. Both original and final water volume must be strictly
    /// positive in every cell; deposition at a drying endpoint requires another
    /// phase/storage owner. Input amounts are mol per bulk transverse m², whereas
    /// boundary and returned concentrations are mol per liquid m³.
    ///
    /// Backward Euler uses the accepted hydraulic face integrals and the actual
    /// final water capacity. This preserves uniform concentration under changing
    /// water storage up to the separately admitted hydraulic residual. Endpoint
    /// theta*D/tau half-cell resistances add in series. An exact zero D seals
    /// molecular exchange only; the original water may still advect solute.
    ///
    /// Donor-cell hydraulic transport and diffusion share one tridiagonal
    /// M-matrix per solute. Native adjacent-pivot LU costs O(cells*species);
    /// positive storage and nonnegative sources give a positive implicit map.
    /// No concentration clipping, limiter, invented water flux or diffusion CFL
    /// is used. Physical cell equations and global constituent balances are
    /// independently checked after solving. Failure/cancellation changes neither
    /// this hydraulic candidate, the original column nor the initial solutes.
    ///
    /// This is one first-order endpoint approximation. Refine physical water
    /// intervals for temporal accuracy; an adaptive water window's cumulative
    /// face report cannot substitute for its individual accepted intervals.
    #[allow(clippy::too_many_arguments)]
    pub fn transport_dilute_solutes<const N: usize>(
        &self,
        original: &LiquidSaturationColumn,
        initial_amount_mol_m2: &[[f64; N]],
        aqueous_diffusivity_m2_s: &[[f64; N]],
        transport_tortuosity: &[f64],
        boundaries: [DiluteSoluteBoundary<N>; 2],
        controls: DiluteSoluteControls<N>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DiluteSoluteStep<N>, DiluteSoluteError> {
        if cancelled() {
            return Err(DiluteSoluteError::Cancelled);
        }
        let n = original.cells.len();
        transport::admit_extent(n, controls)?;
        if self.original_elapsed_s != original.elapsed_s
            || self.original_saturations != original.saturations
            || self.state.cells != original.cells
            || self.state.elapsed_s != original.elapsed_s + self.step.duration_s
        {
            return Err(DiluteSoluteError::Invalid(
                "solute transport requires the exact original hydraulic interval and material cells",
            ));
        }
        let mut before = Vec::with_capacity(n);
        let mut after = Vec::with_capacity(n);
        let mut widths = Vec::with_capacity(n);
        for i in 0..n {
            if i % 128 == 0 && cancelled() {
                return Err(DiluteSoluteError::Cancelled);
            }
            let cell = &original.cells[i];
            before.push(cell.capacity_m() * original.saturations[i]);
            after.push(cell.capacity_m() * self.state.saturations[i]);
            widths.push(cell.thickness_m);
        }
        let mut faces = Vec::with_capacity(n + 1);
        for (i, &flux) in self.step.face_flux_m_s.iter().enumerate() {
            if i % 128 == 0 && cancelled() {
                return Err(DiluteSoluteError::Cancelled);
            }
            let volume = self.step.duration_s * flux;
            if !volume.is_finite() || (flux != 0.0 && volume == 0.0) {
                return Err(DiluteSoluteError::Invalid(
                    "original integrated water flux is unrepresentable",
                ));
            }
            faces.push(volume);
        }
        transport::transport(
            [original.elapsed_s, self.state.elapsed_s],
            self.step.duration_s,
            &before,
            &after,
            &widths,
            &faces,
            initial_amount_mol_m2,
            aqueous_diffusivity_m2_s,
            transport_tortuosity,
            boundaries,
            controls,
            cancelled,
        )
    }
}
