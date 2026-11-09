//! Original source viscosity and body force at the SAME nonlinear endpoint.

use super::*;
use std::cell::RefCell;

impl<F: FluidTransport + ?Sized> MaterialFluxSystem<'_, '_, '_, F> {
    /// BDF1 convection and rheology with a source-owned constant force per cell.
    ///
    /// At each actual objective shear query, `resolve` returns the original
    /// material owner and its additional physical force density [N/m³]. Both
    /// fields are prepared together before the native volume quadrature. The
    /// separate `force_density` supplies any held spatial load and is added
    /// once. Neither closure may advance physical history during iteration.
    ///
    /// Density, consistent inertia, traction, slip, pressure and material heat
    /// retain the same reference normalization as the ordinary rheology entry.
    /// Source forces use original triangle identities, including disconnected
    /// and nonconvex domains; there is no point-to-cell search or interpolation.
    /// The original final nonlinear reassembly includes the endpoint forces,
    /// and native constrained reactions retain their complete volume load.
    ///
    /// Returns the accepted flow/material owners and exact accepted cell forces
    /// in triangle order. Force work is not automatically deposited as heat.
    /// This adds neither an EOS nor density/thermal evolution.
    ///
    /// # Errors
    /// Retains the source's original refusal, and rejects nonfinite forces,
    /// changed density and all native history/boundary/nonlinear failures.
    /// No partial endpoint or replacement history is published on refusal.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn bdf1_navier_stokes_mixed_rheology_with_cell_forces<G: FluidTransport>(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialNavierStokesBoundary<'_>],
        previous_velocity: &[f64],
        timestep_s: f64,
        controls: NavierStokesControls,
        mut resolve: impl FnMut(
            MaterialFluxRheologyPoint,
            &FluidTransportState,
        ) -> Result<(G, [f64; 2]), MaterialError>,
    ) -> Result<(MaterialNavierStokesRheologyStep<G>, Vec<[f64; 2]>), MaterialFluxError> {
        let density = self.transport.density_kg_m3();
        let forces = RefCell::new(vec![None; self.system.mesh.tris.len()]);
        let field = RheologyField::new(
            &self.system.mesh.centroids,
            &self.transport,
            |point, reference| {
                let (material, force) = resolve(point, reference)?;
                if force.iter().any(|value| !value.is_finite()) {
                    return Err(MaterialError::State {
                        what: "finite source-owned cell force density required".into(),
                    });
                }
                forces.borrow_mut()[point.cell] = Some(force);
                Ok(material)
            },
        );
        let endpoint = with_kinematic_boundaries(density, boundaries, |native| {
            self.system.bdf1_navier_stokes_mixed_prepared_cell_forces_checked(
                self.params,
                &|velocity| {
                    let rates = self.system.shear_rates(velocity);
                    if rates.iter().any(|rate| !rate.is_finite()) {
                        return Err(crate::ns::boundary::StokesBoundaryError::InvalidInput(
                            "nonfinite objective shear rate",
                        ));
                    }
                    Ok(rates.iter().enumerate().map(|(cell, &rate)| {
                        field.viscosity(cell, rate)
                    }).collect())
                },
                &|cell, point| {
                    let held = force_density(point);
                    // A refused source only finishes private assembly; its
                    // exact error is returned before a solve or publication.
                    let local = forces.borrow()[cell].unwrap_or([f64::NAN; 2]);
                    core::array::from_fn(|axis| held[axis] / density + local[axis] / density)
                },
                native,
                previous_velocity,
                timestep_s,
                controls,
                &|| field.is_valid(),
            )
        });
        let (endpoint, samples) = field.finish(endpoint.map_err(MaterialFluxError::Boundary))?;
        let forces = forces.into_inner().into_iter().map(|force| {
            force.ok_or(MaterialFluxError::Invalid("accepted cell force is missing"))
        }).collect::<Result<Vec<_>, _>>()?;
        Ok((self.material_navier_stokes_rheology_step(endpoint, samples)?, forces))
    }
}

#[cfg(test)]
#[path = "forcing/tests.rs"]
mod tests;
