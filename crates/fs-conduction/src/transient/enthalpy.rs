//! Stationary reference-mass enthalpy transport on the existing P1 mesh.
//!
//! With specific enthalpy `h` [J/kg], frozen reference density `rho0` and
//! `q = -k(T) grad(T)`, the weak balance is
//! `integral(rho0 N_i dh/dt) + integral(grad(N_i) k grad(T)) = load_i`.
//! Neumann flux is positive OUTWARD; Robin contributes `htc*(T-T_ref)`
//! outward. Existing source, boundary and contact assemblers own all spatial
//! quadrature. Energetic internal variables and geometry are frozen here.
//!
//! Nodal constitutive evaluation gives `T_i=T(h_i)` and the existing P1
//! temperature interpolant. Row-sum storage gives invariant reference masses
//! `m_i=sum_e rho0 V_e/4`. Backward Euler solves, in joules,
//! `R(h)=m*(h-h_old)+dt*(A(T(h))*T(h)-b)=0`. Its Jacobian action is
//! `J_h v=m*v+dt*J_T(T)*(T'(h)*v)`: the derivative scales COLUMNS, not rows.
//! The chart's exact latent plateau has `T'=0`; positive mass keeps these
//! columns nonsingular. No apparent heat capacity or extra latent account is
//! introduced. Conductivity, including its exact `k'(T)` contribution, remains
//! owned by the existing conduction assembly.
//!
//! Uniform and explicitly assigned heterogeneous reference materials share this
//! transport kernel, natural flux/Robin boundaries and finite-resistance
//! contacts. Distinct materials require distinct interface vertices; no nodal
//! mixture law is inferred. Every Dirichlet row is
//! refused: temperature on a latent plateau does not determine enthalpy.
//! Equilibrium density from the chart never changes the stored reference mass.
//! This does not model expansion, material motion, pressure work or remapping.

/// Implicit endpoint derivatives through the accepted physical balance.
pub mod adjoint;
/// Explicit heterogeneous phase charts and frozen reference densities.
pub mod heterogeneous;
/// Implicit ambient radiation with immutable enthalpy history and joule gates.
pub mod radiation;
pub use radiation::adjoint::{EnthalpyRadiationStepGradient, EnthalpyRadiationStepLinearization};

use std::{cell::RefCell, fmt};

use fs_exec::Cx;
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, PhaseStateError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, NewtonStallDiagnosis,
    NonlinearProblem, SolverCallback, SolverRunProgress,
};
use fs_sparse::Csr;

use crate::{
    ConductionError, ConductionMesh, ConductionProblem, InterfaceFlux, RobinFlux, ThermalBc,
    ThermalInterfaces,
    assemble::{
        ASSEMBLY_TILE, AssembledSystem, DofMap, assemble_jacobian_with_optional_interfaces,
        assemble_operator_scaled_with_interfaces,
    },
    solve::energy_balance,
};

/// Spatial work and allocation admission for one prepared reference mass.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EnthalpyBudget {
    /// Maximum nodal enthalpy unknowns, including boundary nodes.
    pub max_vertices: usize,
    /// Maximum tetrahedra used by storage and transport assembly.
    pub max_elements: usize,
}

/// Explicit nonlinear/Krylov work and absolute energy-closure policy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EnthalpyStepConfig {
    /// Shared Newton controls; absolute residual tolerance is in joules.
    /// Restart and cycle bounds apply to each bounded Newton attempt.
    pub newton: NewtonKrylovConfig,
    /// Maximum attempts for this endpoint; no unconverged endpoint is returned.
    pub max_newton_iterations: usize,
    /// Maximum absolute storage-minus-external-input defect, in joules.
    pub energy_tolerance_j: f64,
}

/// A refused construction or step returns no new physical state.
#[derive(Debug, Clone, PartialEq)]
pub enum EnthalpyError {
    /// Invalid shape, nonpositive physical input, or unrepresentable size/work.
    InvalidInput(&'static str),
    /// The caller's spatial admission cap is insufficient.
    Budget {
        /// Which admitted resource exceeded its cap.
        resource: &'static str,
        /// Actual mesh requirement.
        required: usize,
        /// Caller-provided ceiling.
        limit: usize,
    },
    /// Fixed temperature does not select a unique state on a latent plateau.
    UnsupportedDirichlet,
    /// Existing spatial/material/interface refusal, including cancellation.
    Conduction(ConductionError),
    /// An exact nodal chart or derivative query was refused.
    Phase {
        /// Nodal state that could not be resolved.
        vertex: usize,
        /// Original constitutive refusal.
        source: PhaseStateError,
    },
    /// Shared Newton configuration or finite-residual admission failed.
    Newton(NewtonError),
    /// Shared Newton/Krylov work or globalization failed to converge.
    NotConverged(NewtonReport),
    /// A nonlinear callback unwound; the incomplete Newton attempt was discarded.
    SolverCallbackPanicked(SolverCallback),
    /// Independently integrated input did not balance the published storage.
    EnergyBalance {
        /// Actual storage minus net external input, joules.
        residual_j: f64,
        /// Declared absolute allowance, joules.
        tolerance_j: f64,
    },
    /// Finite physical inputs produced an unrepresentable numerical quantity.
    NonFiniteArithmetic,
}

impl fmt::Display for EnthalpyError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidInput(reason) => write!(formatter, "enthalpy input: {reason}"),
            Self::Budget {
                resource,
                required,
                limit,
            } => {
                write!(
                    formatter,
                    "enthalpy {resource} requires {required}, cap {limit}"
                )
            }
            Self::UnsupportedDirichlet => write!(
                formatter,
                "enthalpy transport requires Neumann/Robin boundaries; a prescribed temperature does not select latent enthalpy"
            ),
            Self::Conduction(error) => error.fmt(formatter),
            Self::Phase { vertex, source } => {
                write!(formatter, "enthalpy vertex {vertex}: {source}")
            }
            Self::Newton(error) => error.fmt(formatter),
            Self::NotConverged(report) => write!(
                formatter,
                "enthalpy Newton refused after {} attempts: {:?}; residual {} J",
                report.iterations, report.diagnosis, report.residual_norm
            ),
            Self::SolverCallbackPanicked(callback) => {
                write!(
                    formatter,
                    "enthalpy Newton {callback:?} callback panicked; step discarded"
                )
            }
            Self::EnergyBalance {
                residual_j,
                tolerance_j,
            } => write!(
                formatter,
                "enthalpy energy defect {residual_j} J exceeds {tolerance_j} J"
            ),
            Self::NonFiniteArithmetic => {
                write!(formatter, "enthalpy arithmetic is not representable")
            }
        }
    }
}

impl std::error::Error for EnthalpyError {}

impl From<ConductionError> for EnthalpyError {
    fn from(error: ConductionError) -> Self {
        Self::Conduction(error)
    }
}

/// Accepted nodal enthalpy, temperature and phase with a direct energy audit.
#[derive(Debug, Clone)]
pub struct EnthalpyStepSolution {
    /// Published specific enthalpy at each vertex, J/kg; use as the next history.
    pub specific_enthalpy_j_kg: Vec<f64>,
    /// Chart-resolved nodal absolute temperature, K.
    pub temperature: Vec<f64>,
    /// Chart-resolved liquid mass fractions in [0,1].
    pub liquid_mass_fraction: Vec<f64>,
    /// Endpoint volumetric generation integrated by the existing source rule, W.
    pub source_w: f64,
    /// Endpoint outward Neumann heat, W; negative means imposed heating.
    pub neumann_out_w: f64,
    /// Endpoint outward Robin heat, W.
    pub robin_out_w: f64,
    /// `sum_i m_i*(h_new_i-h_old_i)` of the actual returned state, J.
    pub stored_energy_change_j: f64,
    /// Storage minus dt times net external input, J.
    pub energy_residual_j: f64,
    /// Existing per-region endpoint convection reports.
    pub robin_fluxes: Vec<RobinFlux>,
    /// Internal contact exchanges, excluded from net external input.
    pub contact_fluxes: Vec<InterfaceFlux>,
    /// Shared Newton's true residual and bounded inner-solve telemetry.
    pub newton: NewtonReport,
}

/// Fixed reference nodal masses and one immutable equilibrium enthalpy chart.
#[derive(Debug)]
pub struct EnthalpyBackwardEuler<'m, 'c> {
    mesh: &'m ConductionMesh,
    curve: &'c EquilibriumEnthalpyPhaseCurve,
    // Only the heterogeneous wrapper constructs this private override. The
    // public uniform API, including phase_curve(), keeps its original meaning.
    nodal_curves: Option<Vec<&'c EquilibriumEnthalpyPhaseCurve>>,
    masses: Vec<f64>,
}

impl<'m, 'c> EnthalpyBackwardEuler<'m, 'c> {
    /// Prepare invariant P1 nodal masses from declared reference density [kg/m3].
    /// This is mass quadrature, not a heat-capacity declaration. The curve's
    /// state-dependent density is never substituted into these masses.
    pub fn uniform(
        cx: &Cx<'_>,
        mesh: &'m ConductionMesh,
        curve: &'c EquilibriumEnthalpyPhaseCurve,
        reference_density_kg_m3: f64,
        budget: EnthalpyBudget,
    ) -> Result<Self, EnthalpyError> {
        poll(cx, 0)?;
        let n = mesh.vertex_count();
        for (resource, required, limit) in [
            ("vertices", n, budget.max_vertices),
            ("elements", mesh.element_count(), budget.max_elements),
        ] {
            if required > limit {
                return Err(EnthalpyError::Budget {
                    resource,
                    required,
                    limit,
                });
            }
        }
        if n == 0
            || mesh.element_count() == 0
            || !reference_density_kg_m3.is_finite()
            || reference_density_kg_m3 <= 0.0
        {
            return Err(EnthalpyError::InvalidInput(
                "nonempty mesh and positive finite reference density required",
            ));
        }
        // Mesh/contact assembly is at most dense in its admitted vertex set.
        // Check both dense-index and element-staging products before allocation.
        n.checked_mul(n)
            .and_then(|v| v.checked_mul(24))
            .and_then(|_| mesh.element_count().checked_mul(32 * 24))
            .ok_or(EnthalpyError::InvalidInput(
                "spatial allocation size overflow",
            ))?;
        let mut masses = vec![0.0; n];
        for (element, vertices) in mesh.complex().tets.iter().enumerate() {
            if element % ASSEMBLY_TILE == 0 {
                poll(cx, element)?;
            }
            let contribution =
                finite(reference_density_kg_m3 * mesh.element_volume(element) / 4.0)?;
            if contribution <= 0.0 {
                return Err(EnthalpyError::InvalidInput(
                    "positive reference mass is not representable",
                ));
            }
            for &vertex in vertices {
                masses[vertex as usize] = finite(masses[vertex as usize] + contribution)?;
            }
        }
        for (vertex, &mass) in masses.iter().enumerate() {
            if vertex % ASSEMBLY_TILE == 0 {
                poll(cx, vertex)?;
            }
            if mass <= 0.0 {
                return Err(EnthalpyError::InvalidInput(
                    "every vertex requires positive reference mass",
                ));
            }
        }
        poll(cx, n)?;
        Ok(Self {
            mesh,
            curve,
            nodal_curves: None,
            masses,
        })
    }

    /// Invariant reference masses in nodal mesh order, kg.
    #[must_use]
    pub fn reference_nodal_masses_kg(&self) -> &[f64] {
        &self.masses
    }

    /// Immutable constitutive chart retained by this spatial storage model.
    #[must_use]
    pub const fn phase_curve(&self) -> &EquilibriumEnthalpyPhaseCurve {
        self.curve
    }

    // Callers have already admitted the nodal shape. Keep constitutive lookup
    // common to residuals, Newton tangents, phase reporting and adjoint guards.
    fn curve_for_vertex(&self, vertex: usize) -> &EquilibriumEnthalpyPhaseCurve {
        self.nodal_curves
            .as_ref()
            .map_or(self.curve, |curves| curves[vertex])
    }

    /// Advance from immutable specific-enthalpy history, with current endpoint
    /// source and natural boundary data. Conductivity is evaluated at every
    /// actual nonlinear trial. No returned field escapes before both Newton's
    /// residual gate and the independent absolute joule-balance gate pass.
    ///
    /// The mesh caps bound spatial buffers. Before Newton allocation, checked
    /// products bound the restart basis/Hessenberg sizes and the total maximum
    /// Krylov columns `max_newton_iterations*restart*cycles`. Each shared Newton
    /// attempt uses its declared globalization bounds; Cx is polled between
    /// attempts and inside constitutive, assembly and sparse-action tiles.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_specific_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
    ) -> Result<EnthalpyStepSolution, EnthalpyError> {
        self.admit_step(cx, problem, old_specific_h, dt_s, config)?;
        let context = StepContext {
            storage: self,
            cx,
            problem,
            interfaces,
            old: old_specific_h,
            dt: dt_s,
        };
        // Resolve history and the initial tangent before shared Newton admits it.
        // A chart/domain refusal here retains its original typed error.
        let mut stage = context.stage(old_specific_h)?;
        let mut newton =
            match NewtonKrylovState::new(&stage, old_specific_h.to_vec(), config.newton) {
                Ok(state) => state,
                Err(error) => {
                    return Err(stage.take_failure().unwrap_or(EnthalpyError::Newton(error)));
                }
            };
        poll(cx, 0)?;
        for attempt in 0..config.max_newton_iterations {
            poll(cx, attempt)?;
            let progress = newton.run_cancellable(&stage, 1, cx);
            match progress.progress {
                SolverRunProgress::Complete => {}
                SolverRunProgress::Paused => {
                    return Err(ConductionError::Cancelled {
                        stage: "spatial-enthalpy-solver",
                        at: attempt,
                    }
                    .into());
                }
                SolverRunProgress::CallbackPanicked(callback) => {
                    return Err(EnthalpyError::SolverCallbackPanicked(callback));
                }
                SolverRunProgress::DimensionMismatch { expected, actual } => {
                    return Err(EnthalpyError::Newton(NewtonError::Dimension {
                        problem: actual,
                        state: expected,
                    }));
                }
            }
            if let Some(error) = stage.fatal.take() {
                return Err(error);
            }
            poll(cx, attempt)?;
            let report = progress.report;
            if report.converged {
                return context.finish(newton.x, report, config.energy_tolerance_j);
            }
            if report.diagnosis != Some(NewtonStallDiagnosis::BudgetExhausted)
                || attempt + 1 == config.max_newton_iterations
            {
                return Err(stage
                    .take_failure()
                    .unwrap_or(EnthalpyError::NotConverged(report)));
            }
            // One shared Newton attempt uses one exact tangent point. Rebuild
            // after acceptance rather than freezing k(T) across Newton updates.
            stage = context.stage(&newton.x)?;
        }
        unreachable!("positive Newton attempt budget admitted above")
    }

    fn admit_step(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        old_specific_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
    ) -> Result<(), EnthalpyError> {
        poll(cx, 0)?;
        if !std::ptr::eq(self.mesh, problem.mesh) || old_specific_h.len() != self.masses.len() {
            return Err(EnthalpyError::InvalidInput(
                "reference mass and history must match the exact transport mesh",
            ));
        }
        if problem
            .boundary
            .conditions()
            .iter()
            .any(|bc| matches!(bc, ThermalBc::Dirichlet { .. }))
        {
            return Err(EnthalpyError::UnsupportedDirichlet);
        }
        validate_work(self.masses.len(), dt_s, config)?;
        for condition in problem.boundary.conditions() {
            condition.validate(self.masses.len())?;
        }
        Ok(())
    }
}

struct StepContext<'a, 'm, 'c, 'cx> {
    storage: &'a EnthalpyBackwardEuler<'m, 'c>,
    cx: &'a Cx<'cx>,
    problem: ConductionProblem<'a>,
    interfaces: Option<&'a ThermalInterfaces>,
    old: &'a [f64],
    dt: f64,
}

impl StepContext<'_, '_, '_, '_> {
    fn temperatures(&self, h: &[f64]) -> Result<Vec<f64>, EnthalpyError> {
        let mut temperature = Vec::with_capacity(h.len());
        for (vertex, &value) in h.iter().enumerate() {
            if vertex % ASSEMBLY_TILE == 0 {
                poll(self.cx, vertex)?;
            }
            let state = self
                .storage
                .curve_for_vertex(vertex)
                .state_at_specific_enthalpy(value)
                .map_err(|source| EnthalpyError::Phase { vertex, source })?;
            temperature.push(state.temperature_k());
        }
        Ok(temperature)
    }

    fn assemble(&self, temperature: &[f64]) -> Result<AssembledSystem, EnthalpyError> {
        Ok(assemble_operator_scaled_with_interfaces(
            self.cx,
            self.storage.mesh,
            self.problem.boundary,
            self.problem.material,
            self.problem.source,
            temperature,
            None,
            self.interfaces,
            self.problem.element_materials,
        )?)
    }

    fn stage(&self, h: &[f64]) -> Result<NewtonStage<'_, '_, '_, '_, '_>, EnthalpyError> {
        let temperature = self.temperatures(h)?;
        let picard = self.assemble(&temperature)?;
        let tangent = assemble_jacobian_with_optional_interfaces(
            self.cx,
            self.storage.mesh,
            self.problem.boundary,
            self.problem.material,
            &temperature,
            self.interfaces,
            self.problem.element_materials,
        )?;
        let mut slope = Vec::with_capacity(h.len());
        let mut inverse = Vec::with_capacity(h.len());
        for (vertex, &value) in h.iter().enumerate() {
            if vertex % ASSEMBLY_TILE == 0 {
                poll(self.cx, vertex)?;
            }
            let derivative = self
                .storage
                .curve_for_vertex(vertex)
                .temperature_derivative_at_specific_enthalpy(value)
                .map_err(|source| EnthalpyError::Phase { vertex, source })?;
            let diagonal = finite(
                self.storage.masses[vertex]
                    + self.dt * picard.operator.get(vertex, vertex) * derivative,
            )?;
            if diagonal <= 0.0 {
                return Err(EnthalpyError::InvalidInput(
                    "positive enthalpy Picard diagonal required",
                ));
            }
            slope.push(derivative);
            inverse.push(finite(1.0 / diagonal)?);
        }
        poll(self.cx, h.len())?;
        Ok(NewtonStage {
            context: self,
            tangent,
            slope,
            inverse,
            fatal: RefCell::new(None),
            trial_refusal: RefCell::new(None),
        })
    }

    fn residual(&self, h: &[f64], out: &mut [f64]) -> Result<(), EnthalpyError> {
        let temperature = self.temperatures(h)?;
        let system = self.assemble(&temperature)?;
        for (vertex, value) in out.iter_mut().enumerate() {
            if vertex % ASSEMBLY_TILE == 0 {
                poll(self.cx, vertex)?;
            }
            let (columns, entries) = system.operator.row(vertex);
            let transport: f64 = columns
                .iter()
                .zip(entries)
                .map(|(&column, &entry)| entry * temperature[column])
                .sum();
            *value = finite(
                self.storage.masses[vertex] * (h[vertex] - self.old[vertex])
                    + self.dt * (transport - system.load[vertex]),
            )?;
        }
        Ok(())
    }

    fn finish(
        &self,
        h: Vec<f64>,
        newton: NewtonReport,
        tolerance_j: f64,
    ) -> Result<EnthalpyStepSolution, EnthalpyError> {
        let temperature = self.temperatures(&h)?;
        let system = self.assemble(&temperature)?;
        let mut liquid_mass_fraction = Vec::with_capacity(h.len());
        let mut storage = 0.0;
        for (vertex, &value) in h.iter().enumerate() {
            if vertex % ASSEMBLY_TILE == 0 {
                poll(self.cx, vertex)?;
            }
            storage = finite(storage + self.storage.masses[vertex] * (value - self.old[vertex]))?;
            liquid_mass_fraction.push(
                self.storage
                    .curve_for_vertex(vertex)
                    .state_at_specific_enthalpy(value)
                    .map_err(|source| EnthalpyError::Phase { vertex, source })?
                    .liquid_mass_fraction(),
            );
        }
        // Independent physical source/face integration, not the sum of Newton
        // residuals or an inferred heat input from storage. This shared owner
        // pass and its contact reporting are bounded by the admitted mesh.
        let dofs = DofMap::new(self.problem.boundary, h.len())?;
        let (energy, robin_fluxes) = energy_balance(
            self.storage.mesh,
            self.problem.boundary,
            self.problem.source,
            &system,
            &dofs,
            &temperature,
        );
        poll(self.cx, h.len())?;
        let net = finite(energy.source_w - energy.neumann_out_w - energy.robin_out_w)?;
        let energy_residual_j = finite(storage - self.dt * net)?;
        if energy_residual_j.abs() > tolerance_j {
            return Err(EnthalpyError::EnergyBalance {
                residual_j: energy_residual_j,
                tolerance_j,
            });
        }
        let contact_fluxes = match self.interfaces {
            Some(interfaces) => interfaces.fluxes(&temperature)?,
            None => Vec::new(),
        };
        poll(self.cx, h.len())?;
        Ok(EnthalpyStepSolution {
            specific_enthalpy_j_kg: h,
            temperature,
            liquid_mass_fraction,
            source_w: energy.source_w,
            neumann_out_w: energy.neumann_out_w,
            robin_out_w: energy.robin_out_w,
            stored_energy_change_j: storage,
            energy_residual_j,
            robin_fluxes,
            contact_fluxes,
            newton,
        })
    }
}

// Private adapter valid for precisely one shared Newton attempt. All line-
// search residuals evaluate the full nonlinear model; its Jacobian and Jacobi
// are the exact current Newton point's operators, reused across Krylov columns.
struct NewtonStage<'s, 'a, 'm, 'c, 'cx> {
    context: &'s StepContext<'a, 'm, 'c, 'cx>,
    tangent: Csr,
    slope: Vec<f64>,
    inverse: Vec<f64>,
    fatal: RefCell<Option<EnthalpyError>>,
    trial_refusal: RefCell<Option<EnthalpyError>>,
}

impl NewtonStage<'_, '_, '_, '_, '_> {
    fn take_failure(&self) -> Option<EnthalpyError> {
        self.fatal.take().or_else(|| self.trial_refusal.take())
    }

    fn poll_action(&self, at: usize) -> bool {
        if self.fatal.borrow().is_some() {
            return false;
        }
        if let Err(error) = poll(self.context.cx, at) {
            self.fatal.replace(Some(error));
            return false;
        }
        true
    }
}

impl NonlinearProblem for NewtonStage<'_, '_, '_, '_, '_> {
    fn dimension(&self) -> usize {
        self.slope.len()
    }

    fn residual(&self, x: &[f64], output: &mut [f64]) {
        output.fill(f64::NAN);
        if !self.poll_action(0) {
            return;
        }
        match self.context.residual(x, output) {
            Ok(()) => {
                self.trial_refusal.take();
            }
            Err(error) => {
                // Out-of-domain trial states are rejected by the shared line
                // search. Constitutive/assembly cancellation is always fatal.
                let recoverable = matches!(
                    error,
                    EnthalpyError::Phase {
                        source: PhaseStateError::OutsideEnthalpyDomain { .. },
                        ..
                    } | EnthalpyError::Conduction(ConductionError::OutsideTemperatureSpan { .. })
                );
                if recoverable {
                    self.trial_refusal.replace(Some(error));
                } else {
                    self.fatal.replace(Some(error));
                }
                output.fill(f64::NAN);
            }
        }
    }

    fn jacobian_apply(&self, _x: &[f64], direction: &[f64], output: &mut [f64]) {
        for (row, value) in output.iter_mut().enumerate() {
            if row % ASSEMBLY_TILE == 0 && !self.poll_action(row) {
                output.fill(f64::NAN);
                return;
            }
            let (columns, entries) = self.tangent.row(row);
            let applied: f64 = columns
                .iter()
                .zip(entries)
                .map(|(&column, &entry)| entry * (self.slope[column] * direction[column]))
                .sum();
            *value = self.context.storage.masses[row] * direction[row] + self.context.dt * applied;
        }
    }

    fn preconditioner_apply(
        &self,
        _x: &[f64],
        _outer: usize,
        _inner: usize,
        residual: &[f64],
        output: &mut [f64],
    ) {
        for (row, value) in output.iter_mut().enumerate() {
            if row % ASSEMBLY_TILE == 0 && !self.poll_action(row) {
                output.fill(f64::NAN);
                return;
            }
            *value = residual[row] * self.inverse[row];
        }
    }
}

fn validate_work(n: usize, dt: f64, config: EnthalpyStepConfig) -> Result<(), EnthalpyError> {
    if !dt.is_finite()
        || dt <= 0.0
        || config.max_newton_iterations == 0
        || !config.energy_tolerance_j.is_finite()
        || config.energy_tolerance_j <= 0.0
    {
        return Err(EnthalpyError::InvalidInput(
            "positive finite step, energy tolerance and Newton budget required",
        ));
    }
    let restart = config.newton.linear_restart;
    // Conservative numerical-buffer size: Arnoldi and flexible bases,
    // Hessenberg, rotations, nonlinear trials and physical nodal scratch.
    restart
        .checked_add(1)
        .and_then(|r| r.checked_mul(restart))
        .and_then(|h| {
            n.checked_mul(restart.checked_mul(2)?.checked_add(32)?)?
                .checked_add(h)
        })
        .and_then(|v| v.checked_add(restart.checked_mul(8)?))
        .and_then(|v| v.checked_mul(std::mem::size_of::<f64>()))
        .filter(|&bytes| isize::try_from(bytes).is_ok())
        .ok_or(EnthalpyError::InvalidInput(
            "Newton workspace size overflow",
        ))?;
    config
        .max_newton_iterations
        .checked_mul(restart)
        .and_then(|v| v.checked_mul(config.newton.max_linear_cycles))
        .ok_or(EnthalpyError::InvalidInput("Newton work budget overflow"))?;
    Ok(())
}

fn finite(value: f64) -> Result<f64, EnthalpyError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(EnthalpyError::NonFiniteArithmetic)
    }
}

fn poll(cx: &Cx<'_>, at: usize) -> Result<(), EnthalpyError> {
    cx.checkpoint().map_err(|_| {
        EnthalpyError::Conduction(ConductionError::Cancelled {
            stage: "spatial-enthalpy",
            at,
        })
    })
}
