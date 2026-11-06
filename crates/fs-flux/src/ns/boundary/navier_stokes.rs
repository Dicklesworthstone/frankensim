//! Steady and fully implicit pressure-driven flow on the original BDM1/P0 operators.

use super::{
    AdmittedBoundary, StokesAssembly, StokesBoundary, StokesBoundaryError, StokesRheologyControls,
    StokesSolution,
};
use crate::ns::{FluxParams, FluxSystem};

/// Explicit treatment of convective momentum at an open traction boundary.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OpenBoundaryPolicy {
    /// Impose the supplied physical total traction exactly, including during
    /// backflow. On an `Open` boundary the conservative flux uses the interior
    /// velocity on either sign. `OpenReservoir` retains this literal stress but
    /// supplies its separately admitted incoming velocity. Inward flow can
    /// bring kinetic energy; no whole-scheme backflow stability is inferred.
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
    /// Strong normal wall velocity with a physical tangential Robin law.
    /// Uses the same stress, mass, convection and reaction owners as Velocity.
    /// See StokesBoundary::NavierSlip for its precise traction convention.
    /// Advected momentum uses the interior trace on both signs: the wall law
    /// does not author an incoming reservoir velocity. An impermeable wall
    /// requires zero normal speed relative to this mesh. Nonzero normal data
    /// instead describe a prescribed-normal transport boundary, with literal
    /// conservative inflow and no additional backflow resistance.
    NavierSlip {
        /// Original mesh boundary edge.
        edge: usize,
        /// Endpoint wall velocity [m/s], in the flow history's frame.
        value: &'a dyn Fn([f64; 2]) -> [f64; 2],
        /// Constant edge resistance beta/rho [m/s], finite and nonnegative.
        /// Zero is perfect slip; exact no-slip remains Velocity.
        tangential_resistance_m_s: f64,
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
    /// Literal total traction with independently supplied incoming momentum.
    /// The normal flow rate is solved, not constrained to the reservoir value.
    /// For the retained advector w, the flux is (w.n)u on outflow and
    /// (w.n)u_reservoir on inflow. No backflow resistance is added.
    OpenReservoir {
        /// Boundary edge in the original mesh.
        edge: usize,
        /// Endpoint kinematic total traction [m²/s²], as on `Open`.
        traction: &'a dyn Fn([f64; 2], [f64; 2]) -> [f64; 2],
        /// Incoming physical velocity [m/s] in this mesh's fixed frame.
        /// Two original edge-Gauss samples define one affine trace, retained
        /// through every Picard iterate and the endpoint reaction. Supply all
        /// components, including the reservoir's own normal momentum.
        velocity: &'a dyn Fn([f64; 2]) -> [f64; 2],
    },
}

/// Bounded Picard controls shared with the existing mixed constitutive solver.
pub type NavierStokesControls = StokesRheologyControls;

/// A steady or implicit endpoint accepted against its own convection and viscosity.
pub struct NavierStokesStep {
    /// Velocity and pressure with component-local datums and actual endpoint
    /// residuals, rather than only the last frozen Oseen solve's residuals.
    pub solution: StokesSolution,
    /// Number of completed frozen-advection/constitutive solves.
    pub iterations: usize,
    /// Relative velocity correction before relaxation.
    pub relative_update: f64,
    /// Maximum original and row-scaled nonlinear endpoint residual.
    pub nonlinear_residual: f64,
    /// Retained (edge, stress policy) pairs. Reservoir edges use literal
    /// NaturalTraction and are identified by their boundary reaction owner.
    pub open_boundary_policies: Vec<(usize, OpenBoundaryPolicy)>,
}

impl<'m> FluxSystem<'m> {
    /// Solve steady incompressible flow with mixed velocity, slip and open
    /// total-traction boundaries:
    /// `div(u tensor u) - div(2 nu D(u)) + grad(p) = f`.
    ///
    /// The initial velocity is only a nonlinear guess; there is no mass term,
    /// history load or artificial timestep. Each connected component needs the
    /// original steady velocity/slip support; all-open components are refused.
    /// The same admitted boundary samples, conservative convection and endpoint
    /// reactions are retained. Both nonlinear residuals and the unrelaxed
    /// correction must satisfy the existing controls before publication.
    ///
    /// Force and pressure are divided by the same constant reference density,
    /// velocity is m/s and mesh coordinates metres. Callbacks must remain pure
    /// and fixed. Open-boundary policies retain their literal stress meaning;
    /// this solve does not imply convergence, uniqueness or stability for
    /// arbitrary flow data.
    ///
    /// # Errors
    /// Refuses invalid guesses/controls, unanchored components, inadmissible
    /// boundaries/coefficients, failed linear solves and exhausted iterations.
    pub fn solve_navier_stokes_mixed(
        &self,
        params: FluxParams,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[NavierStokesBoundary<'_>],
        initial_velocity: &[f64],
        controls: NavierStokesControls,
    ) -> Result<NavierStokesStep, StokesBoundaryError> {
        if !params.nu.is_finite() || params.nu <= 0.0 {
            return Err(StokesBoundaryError::InvalidInput(
                "invalid Navier–Stokes viscosity",
            ));
        }
        self.solve_navier_stokes_mixed_variable_viscosity(
            params,
            &|_, _| params.nu,
            force,
            boundaries,
            initial_velocity,
            controls,
        )
    }

    /// Steady convection with a frozen positive spatial `nu(cell, point)`
    /// [m²/s]. The original symmetric stress and harmonic one-sided SIP weights
    /// are retained. `params.nu` is unused; no strain queries or material
    /// evolution are added to this frozen-viscosity path.
    ///
    /// # Errors
    /// Inherits [`Self::solve_navier_stokes_mixed`] acceptance and refusals.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_navier_stokes_mixed_variable_viscosity(
        &self,
        params: FluxParams,
        viscosity: &dyn Fn(usize, [f64; 2]) -> f64,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[NavierStokesBoundary<'_>],
        initial_velocity: &[f64],
        controls: NavierStokesControls,
    ) -> Result<NavierStokesStep, StokesBoundaryError> {
        self.solve_navier_stokes_mixed_endpoint(
            params,
            boundaries,
            initial_velocity,
            controls,
            0.0,
            &[],
            &|admitted, velocity| {
                self.assemble_stokes_admitted_operator(
                    params,
                    viscosity,
                    force,
                    admitted,
                    None,
                    0.0,
                    None,
                    Some(&velocity[..self.n_u]),
                )
            },
            &|| true,
        )
    }

    /// Steady generalized-Newtonian flow with viscosity and convection at the
    /// same accepted velocity. `law(cell, point, sqrt(2 D:D))` returns positive
    /// finite kinematic viscosity [m²/s]; rate is in 1/s. Its temperature and
    /// composition context stays fixed. `params.nu` is unused.
    ///
    /// The initial guess supplies the first objective strain. Every candidate
    /// uses the full symmetric-stress operator, including one-sided face
    /// tractions, and acceptance reassembles that law at the returned endpoint.
    /// No pseudo-time mass, thermal evolution or viscoelastic memory is added.
    ///
    /// # Errors
    /// Inherits [`Self::solve_navier_stokes_mixed`] acceptance and refusals.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_navier_stokes_mixed_rheology(
        &self,
        params: FluxParams,
        law: &dyn Fn(usize, [f64; 2], f64) -> f64,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[NavierStokesBoundary<'_>],
        initial_velocity: &[f64],
        controls: NavierStokesControls,
    ) -> Result<NavierStokesStep, StokesBoundaryError> {
        self.solve_navier_stokes_mixed_rheology_checked(
            params,
            law,
            force,
            boundaries,
            initial_velocity,
            controls,
            &|| true,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn solve_navier_stokes_mixed_rheology_checked(
        &self,
        params: FluxParams,
        law: &dyn Fn(usize, [f64; 2], f64) -> f64,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[NavierStokesBoundary<'_>],
        initial_velocity: &[f64],
        controls: NavierStokesControls,
        source_is_valid: &dyn Fn() -> bool,
    ) -> Result<NavierStokesStep, StokesBoundaryError> {
        self.solve_navier_stokes_mixed_endpoint(
            params,
            boundaries,
            initial_velocity,
            controls,
            0.0,
            &[],
            &|admitted, velocity| {
                let rates = self.shear_rates(velocity);
                if rates.iter().any(|rate| !rate.is_finite()) {
                    return Err(StokesBoundaryError::InvalidInput("nonfinite shear rate"));
                }
                self.assemble_stokes_admitted_operator(
                    params,
                    &|cell, point| law(cell, point, rates[cell]),
                    force,
                    admitted,
                    None,
                    0.0,
                    None,
                    Some(&velocity[..self.n_u]),
                )
            },
            source_is_valid,
        )
    }

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
        self.bdf1_navier_stokes_mixed_step_with_velocity_tolerance(
            params,
            force,
            boundaries,
            previous_velocity,
            timestep_s,
            controls,
            0.0,
        )
    }

    /// The same nonlinear BDF1 owner with an explicitly authored absolute
    /// velocity-correction tolerance [m/s]. The maximum unrelaxed change of
    /// the ORIGINAL affine edge-normal velocity traces must meet this bound
    /// OR the original relative-L2 correction bound. The absolute norm is
    /// max_edges(|delta_mean| + 6*|delta_linear_moment|), exactly the largest
    /// trace-endpoint change. No pressure coefficient participates; increasing
    /// the mesh size does not loosen this maximum.
    /// Both freshly assembled nonlinear equation residuals must still meet
    /// `controls.relative_tolerance`. The returned `relative_update` remains
    /// the actual ratio and may exceed it when the absolute criterion passes.
    /// Zero selects exactly the original relative-only correction criterion.
    ///
    /// # Errors
    /// Original refusals remain; the absolute bound must be finite/nonnegative.
    #[allow(clippy::too_many_arguments)]
    pub fn bdf1_navier_stokes_mixed_step_with_velocity_tolerance(
        &self,
        params: FluxParams,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[NavierStokesBoundary<'_>],
        previous_velocity: &[f64],
        timestep_s: f64,
        controls: NavierStokesControls,
        absolute_velocity_tolerance_m_s: f64,
    ) -> Result<NavierStokesStep, StokesBoundaryError> {
        if !params.nu.is_finite() || params.nu <= 0.0 {
            return Err(StokesBoundaryError::InvalidInput(
                "invalid Navier–Stokes viscosity",
            ));
        }
        self.bdf1_navier_stokes_mixed_variable_viscosity_step_with_velocity_tolerance(
            params,
            &|_, _| params.nu,
            force,
            boundaries,
            previous_velocity,
            timestep_s,
            controls,
            absolute_velocity_tolerance_m_s,
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
        self.bdf1_navier_stokes_mixed_variable_viscosity_step_with_velocity_tolerance(
            params,
            viscosity,
            force,
            boundaries,
            previous_velocity,
            timestep_s,
            controls,
            0.0,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn bdf1_navier_stokes_mixed_variable_viscosity_step_with_velocity_tolerance(
        &self,
        params: FluxParams,
        viscosity: &dyn Fn(usize, [f64; 2]) -> f64,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[NavierStokesBoundary<'_>],
        previous_velocity: &[f64],
        timestep_s: f64,
        controls: NavierStokesControls,
        absolute_velocity_tolerance_m_s: f64,
    ) -> Result<NavierStokesStep, StokesBoundaryError> {
        self.solve_navier_stokes_mixed_bdf1(
            params,
            boundaries,
            previous_velocity,
            timestep_s,
            controls,
            absolute_velocity_tolerance_m_s,
            &|admitted, velocity| {
                self.assemble_stokes_admitted_operator(
                    params,
                    viscosity,
                    force,
                    admitted,
                    None,
                    timestep_s.recip(),
                    Some(previous_velocity),
                    Some(&velocity[..self.n_u]),
                )
            },
            &|| true,
        )
    }

    // Rigid ALE transports ABSOLUTE velocity in current mesh axes with the
    // relative advector u-w_mesh. Q rotates the original history field into
    // those axes under its original mass quadrature. Both material laws use
    // this same nonlinear endpoint/reaction/transport assembly.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn bdf1_navier_stokes_mixed_rigid_ale_step(
        &self,
        params: FluxParams,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[NavierStokesBoundary<'_>],
        previous_velocity: &[f64],
        history_rotation: [f64; 2],
        mesh_velocity: &[f64],
        timestep_s: f64,
        controls: NavierStokesControls,
        absolute_velocity_tolerance_m_s: f64,
    ) -> Result<NavierStokesStep, StokesBoundaryError> {
        self.solve_navier_stokes_mixed_rigid_ale(
            params,
            None,
            force,
            boundaries,
            previous_velocity,
            history_rotation,
            mesh_velocity,
            timestep_s,
            controls,
            absolute_velocity_tolerance_m_s,
            &|| true,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn bdf1_navier_stokes_mixed_rigid_ale_rheology_step_checked(
        &self,
        params: FluxParams,
        law: &dyn Fn(usize, [f64; 2], f64) -> f64,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[NavierStokesBoundary<'_>],
        previous_velocity: &[f64],
        history_rotation: [f64; 2],
        mesh_velocity: &[f64],
        timestep_s: f64,
        controls: NavierStokesControls,
        absolute_velocity_tolerance_m_s: f64,
        source_is_valid: &dyn Fn() -> bool,
    ) -> Result<NavierStokesStep, StokesBoundaryError> {
        self.solve_navier_stokes_mixed_rigid_ale(
            params,
            Some(law),
            force,
            boundaries,
            previous_velocity,
            history_rotation,
            mesh_velocity,
            timestep_s,
            controls,
            absolute_velocity_tolerance_m_s,
            source_is_valid,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn solve_navier_stokes_mixed_rigid_ale(
        &self,
        params: FluxParams,
        law: Option<&dyn Fn(usize, [f64; 2], f64) -> f64>,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[NavierStokesBoundary<'_>],
        previous_velocity: &[f64],
        history_rotation: [f64; 2],
        mesh_velocity: &[f64],
        timestep_s: f64,
        controls: NavierStokesControls,
        absolute_velocity_tolerance_m_s: f64,
        source_is_valid: &dyn Fn() -> bool,
    ) -> Result<NavierStokesStep, StokesBoundaryError> {
        if (law.is_none() && (!params.nu.is_finite() || params.nu <= 0.0))
            || mesh_velocity.len() != self.n_u
            || mesh_velocity.iter().any(|v| !v.is_finite())
            || history_rotation.iter().any(|v| !v.is_finite())
            || (history_rotation[0].hypot(history_rotation[1]) - 1.0).abs() > 32.0 * f64::EPSILON
        {
            return Err(StokesBoundaryError::InvalidInput(
                "finite rigid ALE advection, proper history rotation and positive viscosity required",
            ));
        }
        // Reuse one normal-continuous advection workspace across all frozen
        // solves and the independent nonlinear endpoint reassembly.
        let relative = std::cell::RefCell::new(vec![0.0; self.n_u]);
        self.solve_navier_stokes_mixed_bdf1(
            params,
            boundaries,
            previous_velocity,
            timestep_s,
            controls,
            absolute_velocity_tolerance_m_s,
            &|admitted, velocity| {
                let mut advection = relative.borrow_mut();
                for ((a, u), w) in advection.iter_mut().zip(velocity).zip(mesh_velocity) {
                    *a = *u - *w;
                }
                if advection.iter().any(|v| !v.is_finite()) {
                    return Err(StokesBoundaryError::InvalidInput(
                        "rigid ALE relative advection overflow",
                    ));
                }
                // D(u) is objective and unchanged by subtraction of a rigid
                // mesh field. Use the ORIGINAL absolute candidate so the law
                // never mistakes frame angular speed for physical shear.
                // The Newtonian path adds no strain evaluation/allocation.
                let rates = law.map(|_| self.shear_rates(velocity));
                if rates
                    .as_ref()
                    .is_some_and(|rates| rates.iter().any(|rate| !rate.is_finite()))
                {
                    return Err(StokesBoundaryError::InvalidInput("nonfinite shear rate"));
                }
                self.assemble_stokes_admitted_rotated_history(
                    params,
                    &|cell, point| match (law, rates.as_ref()) {
                        (Some(law), Some(rates)) => law(cell, point, rates[cell]),
                        _ => params.nu,
                    },
                    force,
                    admitted,
                    None,
                    timestep_s.recip(),
                    Some(previous_velocity),
                    Some(advection.as_slice()),
                    Some(history_rotation),
                )
            },
            source_is_valid,
        )
    }

    /// Advance pressure-driven generalized-Newtonian Navier–Stokes with both
    /// convection and `nu(cell, point, gamma)` implicit at the new time level.
    /// Here `gamma = sqrt(2 D(u):D(u))` [1/s] is computed from each cell's
    /// affine BDM1 velocity, and the law returns positive kinematic viscosity
    /// [m²/s]. The strain and advecting velocity always come from the SAME
    /// candidate, including the independently reassembled endpoint residual.
    ///
    /// The actual previous velocity initializes the law; there is no assumed
    /// zero-shear probe or viscosity floor. Captured temperature/composition
    /// and the density stay fixed, while viscosity may vary with position and
    /// strain. `params.nu` is unused. One-sided harmonic SIP stress, admitted
    /// boundary samples, mass history and the explicit backflow policy retain
    /// their original owners. This adds no thermal or viscoelastic evolution.
    ///
    /// # Errors
    /// Inherits mixed nonlinear admission and convergence refusals, and rejects
    /// nonfinite strain or a nonpositive/nonfinite constitutive response.
    #[allow(clippy::too_many_arguments)]
    pub fn bdf1_navier_stokes_mixed_rheology_step(
        &self,
        params: FluxParams,
        law: &dyn Fn(usize, [f64; 2], f64) -> f64,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[NavierStokesBoundary<'_>],
        previous_velocity: &[f64],
        timestep_s: f64,
        controls: NavierStokesControls,
    ) -> Result<NavierStokesStep, StokesBoundaryError> {
        self.bdf1_navier_stokes_mixed_rheology_step_checked(
            params,
            law,
            force,
            boundaries,
            previous_velocity,
            timestep_s,
            controls,
            &|| true,
        )
    }

    /// Advance equal-step BDF2 generalized-Newtonian Navier–Stokes with the
    /// SAME endpoint strain, viscosity and convection as the BDF1 owner.
    ///
    /// `previous_velocity` is u_n and `older_velocity=Some(u_(n-1))` uses
    /// `(3 u_new - 4 u_n + u_(n-1))/(2 dt)`. Both history intervals must equal
    /// `timestep_s`; `None` uses the original BDF1 startup without fabricated
    /// history. Boundary data, force and any captured temperature/composition
    /// describe the new time level. Density and mesh remain fixed.
    ///
    /// Only the mass load uses `(4 u_n-u_(n-1))/3`. The actual u_n initializes
    /// the constitutive law and nonlinear iteration: the mass-history vector
    /// is not a physical strain state. Constructing the mass load requires no
    /// constitutive query or per-iteration allocation. Original endpoint
    /// residual, total-traction, backflow and source-refusal checks remain.
    ///
    /// BDF2 is second-order in time for a smooth resolved solution with
    /// converged nonlinear solves. Its two-history energy is different from
    /// endpoint kinetic energy; endpoint viscous heating is an instantaneous
    /// rate, not an integrated heat or whole-step energy balance. This adds
    /// neither variable-step integration, ALE history nor a general nonlinear
    /// stability guarantee.
    ///
    /// # Errors
    /// Refuses invalid or unrepresentable histories/timestep before callbacks,
    /// and preserves the original coefficient, boundary and solve refusals.
    #[allow(clippy::too_many_arguments)]
    pub fn bdf2_navier_stokes_mixed_rheology_step(
        &self,
        params: FluxParams,
        law: &dyn Fn(usize, [f64; 2], f64) -> f64,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[NavierStokesBoundary<'_>],
        previous_velocity: &[f64],
        older_velocity: Option<&[f64]>,
        timestep_s: f64,
        controls: NavierStokesControls,
    ) -> Result<NavierStokesStep, StokesBoundaryError> {
        self.bdf2_navier_stokes_mixed_rheology_step_checked(
            params,
            law,
            force,
            boundaries,
            previous_velocity,
            older_velocity,
            timestep_s,
            controls,
            &|| true,
        )
    }

    // A material source retains its typed refusal outside the scalar law.
    // The common Picard owner checks it before any solve or final publication.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn bdf1_navier_stokes_mixed_rheology_step_checked(
        &self,
        params: FluxParams,
        law: &dyn Fn(usize, [f64; 2], f64) -> f64,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[NavierStokesBoundary<'_>],
        previous_velocity: &[f64],
        timestep_s: f64,
        controls: NavierStokesControls,
        source_is_valid: &dyn Fn() -> bool,
    ) -> Result<NavierStokesStep, StokesBoundaryError> {
        self.bdf2_navier_stokes_mixed_rheology_step_checked(
            params,
            law,
            force,
            boundaries,
            previous_velocity,
            None,
            timestep_s,
            controls,
            source_is_valid,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn bdf2_navier_stokes_mixed_rheology_step_checked(
        &self,
        params: FluxParams,
        law: &dyn Fn(usize, [f64; 2], f64) -> f64,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[NavierStokesBoundary<'_>],
        previous_velocity: &[f64],
        older_velocity: Option<&[f64]>,
        timestep_s: f64,
        controls: NavierStokesControls,
        source_is_valid: &dyn Fn() -> bool,
    ) -> Result<NavierStokesStep, StokesBoundaryError> {
        // Reuse the original mass/history quadrature and its reaction rows.
        // Construct this load once; keep the real previous velocity as the
        // initial iterate so no material law sees an extrapolated strain.
        let history = if let Some(older) = older_velocity {
            self.check_mixed_bdf1_history(previous_velocity, timestep_s)?;
            if older.len() != self.n_u || older.iter().any(|value| !value.is_finite()) {
                return Err(StokesBoundaryError::InvalidInput(
                    "BDF2 needs two finite original velocity histories",
                ));
            }
            Some(
                previous_velocity
                    .iter()
                    .zip(older)
                    .map(|(&recent, &old)| {
                        // Preserve constant history exactly. Scale first only
                        // when the otherwise finite difference would overflow.
                        let difference = recent - old;
                        let correction = if difference.is_finite() {
                            difference / 3.0
                        } else {
                            recent / 3.0 - old / 3.0
                        };
                        recent + correction
                    })
                    .collect::<Vec<_>>(),
            )
        } else {
            None
        };
        let effective_timestep = if history.is_some() {
            (2.0 / 3.0) * timestep_s
        } else {
            timestep_s
        };
        let mass_history = history.as_deref().unwrap_or(previous_velocity);
        if history.is_some() {
            self.check_mixed_bdf1_history(mass_history, effective_timestep)?;
        }
        self.solve_navier_stokes_mixed_bdf1(
            params,
            boundaries,
            previous_velocity,
            effective_timestep,
            controls,
            0.0,
            &|admitted, velocity| {
                let rates = self.shear_rates(velocity);
                if rates.iter().any(|rate| !rate.is_finite()) {
                    return Err(StokesBoundaryError::InvalidInput("nonfinite shear rate"));
                }
                self.assemble_stokes_admitted_operator(
                    params,
                    &|cell, point| law(cell, point, rates[cell]),
                    force,
                    admitted,
                    None,
                    effective_timestep.recip(),
                    Some(mass_history),
                    Some(&velocity[..self.n_u]),
                )
            },
            source_is_valid,
        )
    }

    // Preserve one admission/reporting path without adding strain work to the
    // frozen-viscosity entry. Only the candidate-dependent operator differs.
    #[allow(clippy::too_many_arguments)]
    fn solve_navier_stokes_mixed_bdf1<'a>(
        &'a self,
        params: FluxParams,
        boundaries: &[NavierStokesBoundary<'_>],
        previous_velocity: &[f64],
        timestep_s: f64,
        controls: NavierStokesControls,
        absolute_velocity_tolerance_m_s: f64,
        assemble: &dyn Fn(
            &AdmittedBoundary,
            &[f64],
        ) -> Result<StokesAssembly<'a, 'm>, StokesBoundaryError>,
        source_is_valid: &dyn Fn() -> bool,
    ) -> Result<NavierStokesStep, StokesBoundaryError> {
        self.check_mixed_bdf1_history(previous_velocity, timestep_s)?;
        self.solve_navier_stokes_mixed_endpoint(
            params,
            boundaries,
            previous_velocity,
            controls,
            absolute_velocity_tolerance_m_s,
            &vec![true; self.pressure_roots.len()],
            assemble,
            source_is_valid,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn solve_navier_stokes_mixed_endpoint<'a>(
        &'a self,
        params: FluxParams,
        boundaries: &[NavierStokesBoundary<'_>],
        initial_velocity: &[f64],
        controls: NavierStokesControls,
        absolute_velocity_tolerance_m_s: f64,
        anchors_velocity: &[bool],
        assemble: &dyn Fn(
            &AdmittedBoundary,
            &[f64],
        ) -> Result<StokesAssembly<'a, 'm>, StokesBoundaryError>,
        source_is_valid: &dyn Fn() -> bool,
    ) -> Result<NavierStokesStep, StokesBoundaryError> {
        if initial_velocity.len() != self.n_u
            || initial_velocity.iter().any(|value| !value.is_finite())
            || controls.maximum_iterations == 0
            || !absolute_velocity_tolerance_m_s.is_finite()
            || absolute_velocity_tolerance_m_s < 0.0
            || !controls.relative_tolerance.is_finite()
            || controls.relative_tolerance <= 0.0
            || !controls.relaxation.is_finite()
            || controls.relaxation <= 0.0
            || controls.relaxation > 1.0
        {
            return Err(StokesBoundaryError::InvalidInput(
                "finite velocity dofs, positive nonlinear budget/tolerance and relaxation in (0,1] required",
            ));
        }
        let admitted = self.admit_navier_stokes_boundaries_with_velocity_anchors(
            params,
            boundaries,
            anchors_velocity,
        )?;
        let solved = self.solve_mixed_picard_with_velocity_tolerance(
            initial_velocity,
            controls,
            absolute_velocity_tolerance_m_s,
            &|velocity| assemble(&admitted, velocity),
            source_is_valid,
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

    pub(super) fn admit_navier_stokes_boundaries(
        &self,
        params: FluxParams,
        boundaries: &[NavierStokesBoundary<'_>],
    ) -> Result<AdmittedBoundary, StokesBoundaryError> {
        self.admit_navier_stokes_boundaries_with_velocity_anchors(
            params,
            boundaries,
            &vec![true; self.pressure_roots.len()],
        )
    }

    fn admit_navier_stokes_boundaries_with_velocity_anchors(
        &self,
        params: FluxParams,
        boundaries: &[NavierStokesBoundary<'_>],
        anchors_velocity: &[bool],
    ) -> Result<AdmittedBoundary, StokesBoundaryError> {
        let stresses: Vec<_> = boundaries
            .iter()
            .map(|boundary| match *boundary {
                NavierStokesBoundary::Velocity { edge, value } => {
                    StokesBoundary::Velocity { edge, value }
                }
                NavierStokesBoundary::NavierSlip {
                    edge,
                    value,
                    tangential_resistance_m_s,
                } => StokesBoundary::NavierSlip {
                    edge,
                    value,
                    tangential_resistance_m_s,
                },
                NavierStokesBoundary::Open { edge, traction, .. }
                | NavierStokesBoundary::OpenReservoir { edge, traction, .. } => {
                    StokesBoundary::TotalTraction {
                        edge,
                        value: traction,
                    }
                }
            })
            .collect();
        let mut admitted =
            self.admit_stokes_problem_with_velocity_anchors(params, &stresses, anchors_velocity)?;
        let mut policies = vec![None; self.mesh.edges.len()];
        let mut reservoirs = boundaries
            .iter()
            .any(|boundary| matches!(boundary, NavierStokesBoundary::OpenReservoir { .. }))
            .then(|| vec![None; self.mesh.edges.len()]);
        // Edge ownership, geometry and stress admission above precede every
        // reservoir query. The same original samples own assembly and reporting;
        // a changing upwind split never queries the source at new positions.
        for boundary in boundaries {
            match *boundary {
                NavierStokesBoundary::Open { edge, backflow, .. } => {
                    policies[edge] = Some(backflow);
                }
                NavierStokesBoundary::OpenReservoir { edge, velocity, .. } => {
                    policies[edge] = Some(OpenBoundaryPolicy::NaturalTraction);
                    let geometry = &self.mesh.edges[edge];
                    let points = crate::bdm::edge_gauss_pub(
                        self.mesh.verts[geometry.verts.0],
                        self.mesh.verts[geometry.verts.1],
                    );
                    let mut samples = [[0.0; 2]; 2];
                    for (q, (point, _)) in points.into_iter().enumerate() {
                        samples[q] = velocity(point);
                        if samples[q].iter().any(|value| !value.is_finite()) {
                            return Err(StokesBoundaryError::NonFiniteBoundary { edge });
                        }
                    }
                    reservoirs.as_mut().expect("reservoir storage")[edge] = Some(samples);
                }
                NavierStokesBoundary::Velocity { .. } | NavierStokesBoundary::NavierSlip { .. } => {
                }
            }
        }
        admitted.open_policies = Some(policies);
        admitted.reservoir_velocities = reservoirs;
        Ok(admitted)
    }
}
