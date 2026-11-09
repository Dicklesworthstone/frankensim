//! Coherent candidate material coefficients and cell-indexed physical loads.

use super::*;

impl<'m> FluxSystem<'m> {
    // The material owner prepares every cell's viscosity and force together
    // before quadrature. The original mixed owner re-evaluates this closure at
    // its final accepted velocity, so its nonlinear residual includes BOTH
    // source fields rather than an earlier buoyancy or rheology iterate.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn bdf1_navier_stokes_mixed_prepared_cell_forces_checked(
        &self,
        params: FluxParams,
        prepare: &dyn Fn(&[f64]) -> Result<Vec<f64>, StokesBoundaryError>,
        force: &dyn Fn(usize, [f64; 2]) -> [f64; 2],
        boundaries: &[NavierStokesBoundary<'_>],
        previous_velocity: &[f64],
        timestep_s: f64,
        controls: NavierStokesControls,
        source_is_valid: &dyn Fn() -> bool,
    ) -> Result<NavierStokesStep, StokesBoundaryError> {
        self.solve_navier_stokes_mixed_bdf1(
            params,
            boundaries,
            previous_velocity,
            timestep_s,
            controls,
            0.0,
            &|admitted, velocity| {
                let viscosities = prepare(velocity)?;
                if viscosities.len() != self.mesh.tris.len() {
                    return Err(StokesBoundaryError::InvalidInput(
                        "one prepared material viscosity and force per triangle required",
                    ));
                }
                self.assemble_stokes_admitted_conservative(
                    params,
                    &|cell, _| viscosities[cell],
                    force,
                    admitted,
                    None,
                    timestep_s.recip(),
                    Some(previous_velocity),
                    Some(&velocity[..self.n_u]),
                    None,
                    None,
                    None,
                )
            },
            source_is_valid,
        )
    }
}
