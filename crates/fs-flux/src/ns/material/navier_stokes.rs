//! Shared-source SI pressure and inertia for nonlinear mixed-boundary flow.

use super::rheology::RheologyField;
use super::{
    MaterialFluxError, MaterialFluxRheologyPoint, MaterialFluxRheologySample, MaterialFluxSystem,
    MaterialStokesSolution,
};
use crate::ns::boundary::{
    NavierStokesBoundary, NavierStokesControls, NavierStokesStep, OpenBoundaryPolicy,
};
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};

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
    /// Literal physical stress and the velocity carried by incoming fluid.
    /// The original material density scales inertia and momentum transport;
    /// the reservoir velocity does not impose the boundary's volume flow rate.
    OpenReservoir {
        /// Original mesh boundary edge.
        edge: usize,
        /// Total traction [Pa] at position [m] and outward normal.
        traction: &'a dyn Fn([f64; 2], [f64; 2]) -> [f64; 2],
        /// Incoming velocity [m/s], in the same frame as the flow history.
        /// Its original two edge-Gauss samples define the retained affine trace.
        velocity: &'a dyn Fn([f64; 2]) -> [f64; 2],
    },
}

/// Accepted nonlinear endpoint and its original shared-material snapshot.
pub struct MaterialNavierStokesStep {
    /// Physical pressure, velocity history, kinetic energy and source state.
    pub solution: MaterialStokesSolution,
    /// Completed frozen-advection/constitutive solves.
    pub iterations: usize,
    /// Final velocity correction before relaxation.
    pub relative_update: f64,
    /// Independently reassembled nonlinear endpoint defect.
    pub nonlinear_residual: f64,
    /// Actual open-boundary policies, retaining the meaning of supplied loads.
    pub open_boundary_policies: Vec<(usize, OpenBoundaryPolicy)>,
}

/// Accepted pressure-driven transient and its exact endpoint material owners.
pub struct MaterialNavierStokesRheologyStep<G> {
    /// Physical flow, convergence evidence and the declared open-boundary laws.
    /// The retained transport snapshot supplies the constant reference density;
    /// endpoint viscosities are in [`Self::samples`], not that initial snapshot.
    pub step: MaterialNavierStokesStep,
    samples: Vec<MaterialFluxRheologySample<G>>,
}

impl<G> MaterialNavierStokesRheologyStep<G> {
    /// One original constitutive owner per triangle in mesh order, evaluated
    /// at the returned velocity and used in its accepted nonlinear residual.
    /// Reading these retained responses does not query the source again.
    #[must_use]
    pub fn samples(&self) -> &[MaterialFluxRheologySample<G>] {
        &self.samples
    }
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
        self.bdf1_navier_stokes_mixed_step_with_velocity_tolerance(
            force_density,
            boundaries,
            previous_velocity,
            timestep_s,
            controls,
            0.0,
        )
    }

    /// The same source-owned nonlinear endpoint with an absolute unrelaxed
    /// affine normal-velocity trace correction bound [m/s], in addition to the
    /// original relative criterion. Both nonlinear equation residual gates
    /// remain mandatory. Zero retains the existing relative-only behavior;
    /// neither mass history nor physical material transport is altered.
    ///
    /// # Errors
    /// Native refusals and invalid finite/nonnegative velocity tolerances.
    pub fn bdf1_navier_stokes_mixed_step_with_velocity_tolerance(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialNavierStokesBoundary<'_>],
        previous_velocity: &[f64],
        timestep_s: f64,
        controls: NavierStokesControls,
        absolute_velocity_tolerance_m_s: f64,
    ) -> Result<MaterialNavierStokesStep, MaterialFluxError> {
        let density = self.transport.density_kg_m3();
        let endpoint = with_kinematic_boundaries(density, boundaries, |native| {
            self.system
                .bdf1_navier_stokes_mixed_step_with_velocity_tolerance(
                    self.params,
                    &|point| {
                        let f = force_density(point);
                        [f[0] / density, f[1] / density]
                    },
                    native,
                    previous_velocity,
                    timestep_s,
                    controls,
                    absolute_velocity_tolerance_m_s,
                )
        })
        .map_err(MaterialFluxError::Boundary)?;
        self.material_navier_stokes_step(endpoint)
    }

    /// Advance pressure-driven generalized-Newtonian BDF1 with the original
    /// material law resolved at each candidate's actual cell strain [1/s].
    /// The returned velocity must satisfy BOTH implicit convection and the
    /// shear-dependent viscous stress in its independently rebuilt residual.
    ///
    /// Each source response must retain the reference density. That same rho
    /// scales dynamic viscosity, inertia, force density [N/m³], total traction
    /// [Pa], explicit backflow resistance, recovered pressure and energy.
    /// Mechanical traction pressure is never used as thermodynamic pressure.
    /// resolve receives the cell centroid [m] and strain invariant; captured
    /// temperature/composition remain fixed during this step. Unchanged
    /// cell/rate queries are cached, and the exact final owners are returned
    /// by [`MaterialNavierStokesRheologyStep::samples`].
    ///
    /// The initial material snapshot remains the reference, while returned
    /// samples supply the endpoint viscosities. No constitutive fit, floor,
    /// thermal update, viscoelastic memory or density evolution is inferred.
    /// Boundary callbacks must be pure and describe the new time level.
    ///
    /// # Errors
    /// Preserves the original source refusal, rejects changed density and all
    /// native boundary/coefficient/nonlinear failures, and never publishes a
    /// fallback source or partial endpoint. History is unchanged on refusal.
    #[allow(clippy::too_many_arguments)]
    pub fn bdf1_navier_stokes_mixed_rheology_step<G: FluidTransport>(
        &self,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[MaterialNavierStokesBoundary<'_>],
        previous_velocity: &[f64],
        timestep_s: f64,
        controls: NavierStokesControls,
        resolve: impl FnMut(MaterialFluxRheologyPoint, &FluidTransportState) -> Result<G, MaterialError>,
    ) -> Result<MaterialNavierStokesRheologyStep<G>, MaterialFluxError> {
        let density = self.transport.density_kg_m3();
        let field = RheologyField::new(&self.system.mesh.centroids, &self.transport, resolve);
        let result = with_kinematic_boundaries(density, boundaries, |native| {
            self.system.bdf1_navier_stokes_mixed_rheology_step_checked(
                self.params,
                &|cell, _, rate| field.viscosity(cell, rate),
                &|point| {
                    let f = force_density(point);
                    [f[0] / density, f[1] / density]
                },
                native,
                previous_velocity,
                timestep_s,
                controls,
                &|| field.is_valid(),
            )
        })
        .map_err(MaterialFluxError::Boundary);
        let (endpoint, samples) = field.finish(result)?;
        Ok(MaterialNavierStokesRheologyStep {
            step: self.material_navier_stokes_step(endpoint)?,
            samples,
        })
    }

    fn material_navier_stokes_step(
        &self,
        endpoint: NavierStokesStep,
    ) -> Result<MaterialNavierStokesStep, MaterialFluxError> {
        Ok(MaterialNavierStokesStep {
            solution: self.material_stokes_solution(endpoint.solution)?,
            iterations: endpoint.iterations,
            relative_update: endpoint.relative_update,
            nonlinear_residual: endpoint.nonlinear_residual,
            open_boundary_policies: endpoint.open_boundary_policies,
        })
    }
}

fn with_kinematic_boundaries<T>(
    density: f64,
    boundaries: &[MaterialNavierStokesBoundary<'_>],
    solve: impl FnOnce(&[NavierStokesBoundary<'_>]) -> T,
) -> T {
    // One homogeneous closure allocation; no boxed callback per edge.
    let scaled: Vec<_> = boundaries
        .iter()
        .map(|boundary| {
            move |point, normal| match boundary {
                MaterialNavierStokesBoundary::Open { traction, .. }
                | MaterialNavierStokesBoundary::OpenReservoir { traction, .. } => {
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
            MaterialNavierStokesBoundary::OpenReservoir { edge, velocity, .. } => {
                NavierStokesBoundary::OpenReservoir {
                    edge,
                    traction,
                    velocity,
                }
            }
        })
        .collect();
    solve(&native)
}
