//! Shared-source SI pressure and inertia for nonlinear mixed-boundary flow.

use super::{MaterialFluxError, MaterialFluxSystem, MaterialStokesSolution};
use crate::ns::boundary::{NavierStokesBoundary, NavierStokesControls, OpenBoundaryPolicy};
use fs_material::fluid_transport::FluidTransport;

/// Physical endpoint conditions on the original FEEC boundary edges.
#[derive(Clone, Copy)]
pub enum MaterialNavierStokesBoundary<'a> {
    /// Prescribed physical velocity.
    Velocity {
        /// Original mesh boundary edge.
        edge: usize,
        /// Velocity [m/s] at position [m].
        value: &'a dyn Fn([f64; 2]) -> [f64; 2],
    },
    /// Open physical stress and the separately selected backflow treatment.
    Open {
        /// Original mesh boundary edge.
        edge: usize,
        /// Total traction [Pa], evaluated at position [m] and outward normal.
        traction: &'a dyn Fn([f64; 2], [f64; 2]) -> [f64; 2],
        /// NaturalTraction imposes the supplied stress literally. The
        /// dissipative alternative adds `rho*min(u.n,0)*u` to that stress,
        /// using this SAME immutable source density, rather than changing p.
        backflow: OpenBoundaryPolicy,
    },
}

/// Accepted nonlinear endpoint and its original shared-material snapshot.
pub struct MaterialNavierStokesStep {
    /// Physical pressure, velocity history, kinetic energy and source state.
    pub solution: MaterialStokesSolution,
    /// Completed frozen-advection solves.
    pub iterations: usize,
    /// Final velocity correction before relaxation.
    pub relative_update: f64,
    /// Independently reassembled nonlinear endpoint defect.
    pub nonlinear_residual: f64,
    /// Actual open-boundary policies, retaining the meaning of supplied loads.
    pub open_boundary_policies: Vec<(usize, OpenBoundaryPolicy)>,
}

impl<F: FluidTransport + ?Sized> MaterialFluxSystem<'_, '_, '_, F> {
    /// Advance fully implicit BDF1 Navier–Stokes with endpoint force density
    /// [N/m³], velocities [m/s] and physical total tractions [Pa]. The same
    /// retained material density scales inertia, convection, viscous stress,
    /// force, traction, pressure recovery, backflow resistance and energy.
    ///
    /// This source snapshot has fixed density and Newtonian viscosity during
    /// the step. Mechanical pressure is not an EOS state. The underlying
    /// checked Picard iteration requires the original and row-scaled nonlinear
    /// residuals and the unrelaxed correction, with no changed history or
    /// partially accepted source state on refusal. Every callback must be pure.
    ///
    /// # Errors
    /// Preserves native boundary, source-independent load, nonlinear and linear
    /// refusal. Nonfinite recovered physical pressure or energy also refuses.
    pub fn bdf1_navier_stokes_mixed_step(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialNavierStokesBoundary<'_>],
        previous_velocity: &[f64],
        timestep_s: f64,
        controls: NavierStokesControls,
    ) -> Result<MaterialNavierStokesStep, MaterialFluxError> {
        let density = self.transport.density_kg_m3();
        // One homogeneous closure allocation; no boxed callback per edge.
        let scaled: Vec<_> = boundaries
            .iter()
            .map(|boundary| {
                move |point, normal| match boundary {
                    MaterialNavierStokesBoundary::Open { traction, .. } => {
                        let t = traction(point, normal);
                        [t[0] / density, t[1] / density]
                    }
                    MaterialNavierStokesBoundary::Velocity { .. } => {
                        unreachable!("only open edges request traction samples")
                    }
                }
            })
            .collect();
        let native: Vec<_> = boundaries
            .iter()
            .zip(&scaled)
            .map(|(boundary, traction)| match *boundary {
                MaterialNavierStokesBoundary::Velocity { edge, value } => {
                    NavierStokesBoundary::Velocity { edge, value }
                }
                MaterialNavierStokesBoundary::Open { edge, backflow, .. } => {
                    NavierStokesBoundary::Open {
                        edge,
                        traction,
                        backflow,
                    }
                }
            })
            .collect();
        let endpoint = self
            .system
            .bdf1_navier_stokes_mixed_step(
                self.params,
                &|point| {
                    let f = force_density(point);
                    [f[0] / density, f[1] / density]
                },
                &native,
                previous_velocity,
                timestep_s,
                controls,
            )
            .map_err(MaterialFluxError::Boundary)?;
        Ok(MaterialNavierStokesStep {
            solution: self.material_stokes_solution(endpoint.solution)?,
            iterations: endpoint.iterations,
            relative_update: endpoint.relative_update,
            nonlinear_residual: endpoint.nonlinear_residual,
            open_boundary_policies: endpoint.open_boundary_policies,
        })
    }
}
