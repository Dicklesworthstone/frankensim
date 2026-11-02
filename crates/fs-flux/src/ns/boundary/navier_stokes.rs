//! Fully implicit pressure-driven flow on the original BDM1/P0 operators.

use super::{
    AdmittedBoundary, StokesAssembly, StokesBoundary, StokesBoundaryError, StokesRheologyControls,
    StokesSolution,
};
use crate::ns::{FluxParams, FluxSystem};

/// Explicit treatment of convective momentum at an open traction boundary.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OpenBoundaryPolicy {
    /// Impose the supplied physical total traction exactly, including during
    /// backflow. The conservative flux uses the interior velocity on either
    /// sign. Inward flow can bring kinetic energy into the domain; no backflow
    /// stability or externally supplied inflow-velocity model is inferred.
    NaturalTraction,
    /// Add the declared resistance `rho * min(u.n,0) * u` to the supplied
    /// physical traction, using the endpoint velocity. In kinematic units,
    /// `(2 nu D(u)-p I)n = t + min(u.n,0)u`. This is NOT the same boundary
    /// load as NaturalTraction during inflow, and is tied to the fixed frame.
    /// The resulting conservative boundary flux is `max(u.n,0)u`.
    /// With divergence-free advection its boundary contribution to kinetic
    /// energy dissipation is `|u.n| |u|²/2`, excluding the supplied load work.
    /// It does not establish whole-solver stability for arbitrary penalties,
    /// data, nonlinear iterations or timesteps.
    DissipativeBackflow,
}

/// One velocity or explicitly selected open-traction owner per boundary edge.
#[derive(Clone, Copy)]
pub enum NavierStokesBoundary<'a> {
    /// Prescribed physical velocity, strong in normal BDM moments and weak in
    /// tangential SIP stress. Incoming advected velocity uses these same data.
    Velocity {
        /// Boundary edge in the original mesh.
        edge: usize,
        /// Endpoint velocity [m/s] at physical position [m].
        value: &'a dyn Fn([f64; 2]) -> [f64; 2],
    },
    /// Open total stress with an explicit treatment of reverse flow.
    Open {
        /// Boundary edge in the original mesh.
        edge: usize,
        /// Endpoint kinematic traction [m²/s²] at position and outward normal.
        /// A pressure-only load is `-p*n`, where p is divided by density.
        traction: &'a dyn Fn([f64; 2], [f64; 2]) -> [f64; 2],
        /// Whether the supplied traction is literal or has added resistance.
        backflow: OpenBoundaryPolicy,
    },
}

/// Bounded Picard controls shared with the existing mixed constitutive solver.
pub type NavierStokesControls = StokesRheologyControls;

/// A nonlinear endpoint accepted against its own convection operator.
pub struct NavierStokesStep {
    /// Velocity and pressure with component-local datums and actual endpoint
    /// residuals, rather than only the last frozen Oseen solve's residuals.
    pub solution: StokesSolution,
    /// Number of completed frozen-advection solves.
    pub iterations: usize,
    /// Relative velocity correction before relaxation.
    pub relative_update: f64,
    /// Maximum original and row-scaled nonlinear endpoint residual.
    pub nonlinear_residual: f64,
    /// Retained (edge, policy) pairs for the physical interpretation of loads.
    pub open_boundary_policies: Vec<(usize, OpenBoundaryPolicy)>,
}

impl<'m> FluxSystem<'m> {
    /// Assemble one BDF1 Oseen operator with an explicitly frozen advecting
    /// velocity. The nonlinear endpoint is supplied by
    /// [`Self::bdf1_navier_stokes_mixed_step`]. All callbacks belong to the new
    /// time level and must be pure. Pressure/force are divided by one constant
    /// reference density, velocity is m/s, mesh coordinates metres, time seconds.
    ///
    /// Reuses the original mass, symmetric stress, conservative cell convection
    /// and normal-continuous upwind face assembly. Reversing affine face fluxes
    /// use the original exact split quadrature. The matrix remains available
    /// for its existing linear transpose solve; it is not a nonlinear Jacobian.
    ///
    /// # Errors
    /// Refuses invalid timestep/history/advection, viscosity, ownership,
    /// component-local flux, loads and unrepresentable assembly.
    #[allow(clippy::too_many_arguments)]
    pub fn assemble_navier_stokes_mixed_bdf1<'a>(
        &'a self,
        params: FluxParams,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[NavierStokesBoundary<'_>],
        advecting_velocity: &[f64],
        previous_velocity: &[f64],
        timestep_s: f64,
    ) -> Result<StokesAssembly<'a, 'm>, StokesBoundaryError> {
        if !params.nu.is_finite() || params.nu <= 0.0 {
            return Err(StokesBoundaryError::InvalidInput(
                "invalid Navier–Stokes viscosity",
            ));
        }
        self.assemble_navier_stokes_mixed_bdf1_variable_viscosity(
            params,
            &|_, _| params.nu,
            force,
            boundaries,
            advecting_velocity,
            previous_velocity,
            timestep_s,
        )
    }

    /// The same Oseen assembly with positive endpoint `nu(cell, point)` [m²/s].
    /// This field stays frozen within the step; `params.nu` is unused. No
    /// variable density or implicit constitutive evolution is inferred.
    ///
    /// # Errors
    /// Inherits checked mixed-boundary assembly and coefficient refusals.
    #[allow(clippy::too_many_arguments)]
    pub fn assemble_navier_stokes_mixed_bdf1_variable_viscosity<'a>(
        &'a self,
        params: FluxParams,
        viscosity: &dyn Fn(usize, [f64; 2]) -> f64,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[NavierStokesBoundary<'_>],
        advecting_velocity: &[f64],
        previous_velocity: &[f64],
        timestep_s: f64,
    ) -> Result<StokesAssembly<'a, 'm>, StokesBoundaryError> {
        self.check_mixed_bdf1_history(previous_velocity, timestep_s)?;
        if advecting_velocity.len() != self.n_u
            || advecting_velocity.iter().any(|value| !value.is_finite())
        {
            return Err(StokesBoundaryError::InvalidInput(
                "finite advecting velocity moments required",
            ));
        }
        let admitted = self.admit_navier_stokes_boundaries(params, boundaries)?;
        self.assemble_stokes_admitted_operator(
            params,
            viscosity,
            force,
            &admitted,
            None,
            timestep_s.recip(),
            Some(previous_velocity),
            Some(advecting_velocity),
        )
    }

    /// Solve the fully implicit incompressible BDF1 equation
    /// `(u-u_old)/dt + div(u tensor u) - div(2 nu D(u)) + grad(p) = f`.
    ///
    /// The original bounded mixed Picard owner reassembles convection at the
    /// endpoint before accepting both original and row-scaled residuals and
    /// the unrelaxed correction. A frozen-advection solve alone is insufficient.
    /// Boundary data are sampled once, and all iterates keep the same mass
    /// history and source viscosity. Physical inertia admits all-open domains;
    /// sealed components retain independent pressure gauges and flux checks.
    /// No history, pressure, or partial candidate is published on refusal.
    ///
    /// This is fixed-mesh, constant-density flow. Natural backflow can supply
    /// energy; choosing resistance changes the actual boundary stress as
    /// documented in [`OpenBoundaryPolicy`]. Picard convergence, turbulence,
    /// nonlinear adjoints and scalable-throughput claims are not supplied.
    ///
    /// # Errors
    /// Refuses invalid inputs, failed linear solves or exhausted nonlinear work.
    #[allow(clippy::too_many_arguments)]
    pub fn bdf1_navier_stokes_mixed_step(
        &self,
        params: FluxParams,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[NavierStokesBoundary<'_>],
        previous_velocity: &[f64],
        timestep_s: f64,
        controls: NavierStokesControls,
    ) -> Result<NavierStokesStep, StokesBoundaryError> {
        if !params.nu.is_finite() || params.nu <= 0.0 {
            return Err(StokesBoundaryError::InvalidInput(
                "invalid Navier–Stokes viscosity",
            ));
        }
        self.bdf1_navier_stokes_mixed_variable_viscosity_step(
            params,
            &|_, _| params.nu,
            force,
            boundaries,
            previous_velocity,
            timestep_s,
            controls,
        )
    }

    /// Fully implicit convection with a frozen positive spatial viscosity
    /// field. The original harmonic one-sided SIP stress weighting is retained.
    ///
    /// # Errors
    /// Inherits [`Self::bdf1_navier_stokes_mixed_step`] acceptance and refusals.
    #[allow(clippy::too_many_arguments)]
    pub fn bdf1_navier_stokes_mixed_variable_viscosity_step(
        &self,
        params: FluxParams,
        viscosity: &dyn Fn(usize, [f64; 2]) -> f64,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[NavierStokesBoundary<'_>],
        previous_velocity: &[f64],
        timestep_s: f64,
        controls: NavierStokesControls,
    ) -> Result<NavierStokesStep, StokesBoundaryError> {
        self.check_mixed_bdf1_history(previous_velocity, timestep_s)?;
        if controls.maximum_iterations == 0
            || !controls.relative_tolerance.is_finite()
            || controls.relative_tolerance <= 0.0
            || !controls.relaxation.is_finite()
            || controls.relaxation <= 0.0
            || controls.relaxation > 1.0
        {
            return Err(StokesBoundaryError::InvalidInput(
                "positive nonlinear budget/tolerance and relaxation in (0,1] required",
            ));
        }
        let admitted = self.admit_navier_stokes_boundaries(params, boundaries)?;
        let solved = self.solve_mixed_picard(
            previous_velocity,
            controls,
            &|velocity| {
                self.assemble_stokes_admitted_operator(
                    params,
                    viscosity,
                    force,
                    &admitted,
                    None,
                    timestep_s.recip(),
                    Some(previous_velocity),
                    Some(&velocity[..self.n_u]),
                )
            },
            &|| true,
        )?;
        Ok(NavierStokesStep {
            solution: solved.solution,
            iterations: solved.iterations,
            relative_update: solved.relative_update,
            nonlinear_residual: solved.nonlinear_residual,
            open_boundary_policies: admitted
                .open_policies
                .as_ref()
                .expect("open policies")
                .iter()
                .enumerate()
                .filter_map(|(edge, policy)| policy.map(|p| (edge, p)))
                .collect(),
        })
    }

    fn check_mixed_bdf1_history(
        &self,
        previous: &[f64],
        dt: f64,
    ) -> Result<(), StokesBoundaryError> {
        if previous.len() != self.n_u
            || previous.iter().any(|value| !value.is_finite())
            || !dt.is_finite()
            || dt <= 0.0
            || !dt.recip().is_finite()
            || dt.recip() <= 0.0
        {
            return Err(StokesBoundaryError::InvalidInput(
                "finite history and representable positive timestep required",
            ));
        }
        Ok(())
    }

    fn admit_navier_stokes_boundaries(
        &self,
        params: FluxParams,
        boundaries: &[NavierStokesBoundary<'_>],
    ) -> Result<AdmittedBoundary, StokesBoundaryError> {
        let stresses: Vec<_> = boundaries
            .iter()
            .map(|boundary| match *boundary {
                NavierStokesBoundary::Velocity { edge, value } => {
                    StokesBoundary::Velocity { edge, value }
                }
                NavierStokesBoundary::Open { edge, traction, .. } => {
                    StokesBoundary::TotalTraction {
                        edge,
                        value: traction,
                    }
                }
            })
            .collect();
        let mut admitted = self.admit_stokes_problem_with_velocity_anchors(
            params,
            &stresses,
            &vec![true; self.pressure_roots.len()],
        )?;
        let mut policies = vec![None; self.mesh.edges.len()];
        for boundary in boundaries {
            if let NavierStokesBoundary::Open { edge, backflow, .. } = *boundary {
                policies[edge] = Some(backflow);
            }
        }
        admitted.open_policies = Some(policies);
        Ok(admitted)
    }
}
