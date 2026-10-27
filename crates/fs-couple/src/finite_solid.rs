//! Source-bound finite-deformation solid motion on the native solid/time owners.
//!
//! The same immutable material card supplies a finite-strain constitutive model
//! and reference density at one exact typed external context. The reference
//! mass stays fixed as the body deforms. Native hyper3 supplies Piola forces,
//! consistent tangents and stored energy; native generalized-alpha at rho_inf=1
//! supplies implicit-midpoint motion. No constitutive law or integrator is
//! duplicated here.
//!
//! This bounded path has fixed supports, held referential nodal/body loads,
//! compressible isothermal hyperelasticity and conforming P1 tetrahedra.
//! Nonlinear midpoint has a signed energy defect: it is reported and gated,
//! never converted into physical heat. No thermal eigen-deformation, plasticity,
//! contact, moving supports, follower pressure, mixed incompressible formulation,
//! or continuous-in-time nonintersection claim is made.

use std::cell::RefCell;
use std::sync::Arc;

use fs_exec::Cx;
use fs_matdb::{MaterialCard, QueryPoint};
use fs_material::hyper::Hyperelastic;
use fs_material::hyper::cards::{HyperelasticMaterial, ResolvedHyperelastic};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_qty::Density;
use fs_solid::hyper3::{HyperTetError, HyperTetProblem, HyperelasticLawField};
use fs_solver::{FlexiblePreconditioner, SolveReport};
use fs_sparse::{Coo, Csr};
use fs_time::galpha::initialization::{
    InitialSolveConfig, InitialSolveError, second_order_acceleration,
};
use fs_time::galpha::{
    ImplicitSolveConfig, ImplicitStepTelemetry, OperatorGeneralizedAlpha,
    SecondOrderOperatorWeights, SecondOrderProblem, SecondOrderState, TimeSolveError,
};

/// One immutable, same-card binding of constitutive law and reference density.
#[derive(Debug, Clone, PartialEq)]
pub struct FiniteSolidMaterial {
    law: ResolvedHyperelastic,
    density: ResolvedMaterialStatePoint,
    model_index: usize,
}

impl FiniteSolidMaterial {
    /// Resolve an explicit model and density from the SAME material card.
    ///
    /// Temperature, pressure and composition may describe held preparation
    /// conditions. Stress, strain, rate and history source coordinates require
    /// evolving owners; neither model nor density silently freezes those axes.
    pub fn resolve(
        card: &MaterialCard,
        model_index: usize,
        context: &QueryPoint,
        selection: MaterialPropertySelection,
    ) -> Result<Self, FiniteSolidError> {
        let model = card.models().get(model_index)
            .ok_or(FiniteSolidError::Invalid("model index is outside the material card"))?;
        let law = ResolvedHyperelastic::try_from_model_card_at(model, context)
            .map_err(|error| FiniteSolidError::Source(error.to_string()))?;
        let requirement = ScalarPropertyRequirement::try_new(
            DENSITY_PROPERTY, Density::DIMS, ScalarAdmissibility::StrictlyPositive,
        ).map_err(|error| FiniteSolidError::Source(error.to_string()))?;
        let density = resolve_material_state_point(
            card, law.context(), &[requirement], selection,
        ).map_err(|error| FiniteSolidError::Source(error.to_string()))?;
        if !density.properties()[0].is_static_at_fixed_external_conditions() {
            return Err(FiniteSolidError::Source(
                "reference density requires validity-wide support at fixed external conditions; evolving mechanical axes are not owned".into(),
            ));
        }
        Ok(Self { law, density, model_index })
    }

    /// Original admitted constitutive model and exact external context.
    pub const fn law(&self) -> &ResolvedHyperelastic { &self.law }

    /// Source density bundle including card identity, typed query and receipt.
    pub const fn density(&self) -> &ResolvedMaterialStatePoint { &self.density }

    /// Explicit model position in the immutable material card.
    pub const fn model_index(&self) -> usize { self.model_index }

    /// Density per reference volume; fixed for this material motion.
    pub fn reference_density_kg_m3(&self) -> f64 { self.density.properties()[0].value_si() }
}

impl HyperelasticMaterial for FiniteSolidMaterial {
    fn hyperelastic(&self) -> &Hyperelastic { self.law.law() }
}

/// Budget for the existing native consistent-initial-acceleration solve.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FiniteSolidInitialization {
    /// Bounded FGMRES controls used by fs-time.
    pub solve: InitialSolveConfig,
    /// Scalar workspace ceiling passed to the native mass solve.
    pub maximum_workspace_components: usize,
}

impl Default for FiniteSolidInitialization {
    fn default() -> Self {
        Self {
            solve: InitialSolveConfig::default(),
            maximum_workspace_components: 1_000_000,
        }
    }
}

/// Controls for one candidate mechanical interval.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FiniteSolidControls {
    /// Existing Newton/Krylov controls; no second time solver is introduced.
    pub solve: ImplicitSolveConfig,
    /// Absolute admission gate on signed kinetic + strain - load-work defect.
    /// This is an integration check, not a physical dissipation or accuracy bound.
    pub energy_tolerance_j: f64,
    /// Total accepted-step cap, including the native retained history.
    pub maximum_steps: usize,
}

impl Default for FiniteSolidControls {
    fn default() -> Self {
        Self {
            solve: ImplicitSolveConfig::default(),
            energy_tolerance_j: 1.0e-8,
            maximum_steps: 100_000,
        }
    }
}

/// Structured refusal; a failed candidate never changes an accepted body.
#[derive(Debug, Clone)]
pub enum FiniteSolidError {
    /// Invalid dimensions, controls, boundary state or representability.
    Invalid(&'static str),
    /// Original material source refused its fixed-condition binding.
    Source(String),
    /// Original native solid geometry/material/cancellation refusal.
    Solid(HyperTetError),
    /// Original native consistent-initialization refusal.
    Initialization(InitialSolveError),
    /// Original native nonlinear time-step refusal.
    Time(TimeSolveError),
    /// Signed mechanical integration defect exceeded its requested gate.
    Energy { defect_j: f64, allowance_j: f64 },
}

impl std::fmt::Display for FiniteSolidError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Invalid(message) => write!(f, "FS-COUPLE-FINITE-SOLID: {message}"),
            Self::Source(message) => write!(f, "FS-COUPLE-FINITE-SOLID-SOURCE: {message}"),
            Self::Solid(error) => std::fmt::Display::fmt(error, f),
            Self::Initialization(error) => std::fmt::Display::fmt(error, f),
            Self::Time(error) => std::fmt::Display::fmt(error, f),
            Self::Energy { defect_j, allowance_j } => write!(
                f, "FS-COUPLE-FINITE-SOLID-ENERGY: signed defect {defect_j:e} J exceeds {allowance_j:e} J"
            ),
        }
    }
}
impl std::error::Error for FiniteSolidError {}
impl From<HyperTetError> for FiniteSolidError {
    fn from(error: HyperTetError) -> Self { Self::Solid(error) }
}
impl From<InitialSolveError> for FiniteSolidError {
    fn from(error: InitialSolveError) -> Self { Self::Initialization(error) }
}
impl From<TimeSolveError> for FiniteSolidError {
    fn from(error: TimeSolveError) -> Self { Self::Time(error) }
}

/// Physical work and storage for one accepted candidate.
#[derive(Debug, Clone)]
pub struct FiniteSolidStep {
    /// Original accepted Newton/Krylov/time receipt.
    pub time: ImplicitStepTelemetry,
    pub kinetic_before_j: f64,
    pub kinetic_after_j: f64,
    pub strain_before_j: f64,
    pub strain_after_j: f64,
    /// Held referential load dotted with actual nodal displacement increment.
    pub source_work_j: f64,
    /// Delta kinetic + delta strain - source work. Signed numerical defect;
    /// it is never exported as thermal energy.
    pub integration_defect_j: f64,
    /// Midpoint support forces M delta-v/dt + internal - external, in original
    /// prescribed-DOF order. Supports are stationary and perform zero work.
    pub support_reactions_n: Vec<(usize, f64)>,
}

/// Accepted finite-strain motion with immutable source and reference operators.
#[derive(Debug, Clone)]
pub struct FiniteSolidBody<'a> {
    problem: HyperTetProblem<'a, FiniteSolidMaterial>,
    mass: Arc<Csr>,
    free_mass: Arc<Csr>,
    free: Arc<[usize]>,
    fixed: Arc<[f64]>,
    reference_mass_kg: f64,
    initialization: SolveReport,
    state: SecondOrderState,
}

impl<'a> FiniteSolidBody<'a> {
    /// Bind native geometry and same-card materials, then solve the initial
    /// physical acceleration with fs-time's existing bounded mass solver.
    ///
    /// Prescribed displacements are already attained at the initial time.
    /// Their velocities must be exactly zero; moving supports need another rung.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        problem: HyperTetProblem<'a, FiniteSolidMaterial>,
        displacement_m: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        time_s: f64,
        controls: FiniteSolidInitialization,
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteSolidError> {
        checkpoint(cx)?;
        if !time_s.is_finite()
            || displacement_m.len() != problem.nodes_m.len()
            || velocity_m_s.len() != problem.nodes_m.len()
            || displacement_m.iter().flatten().chain(velocity_m_s.iter().flatten())
                .any(|value| !value.is_finite())
            || problem.nodes_m.len() > problem.budget.maximum_nodes
            || problem.tetrahedra.len() > problem.budget.maximum_tetrahedra
        {
            return Err(FiniteSolidError::Invalid("invalid initial state or mesh budget"));
        }
        let dofs = problem.nodes_m.len().checked_mul(3)
            .ok_or(FiniteSolidError::Invalid("DOF count overflow"))?;
        let densities: Vec<f64> = match problem.material {
            HyperelasticLawField::Uniform(material) =>
                vec![material.reference_density_kg_m3(); problem.tetrahedra.len()],
            HyperelasticLawField::PerElement(materials) => materials.iter()
                .map(FiniteSolidMaterial::reference_density_kg_m3).collect(),
        };
        let mass = problem.reference_mass_matrix(&densities, cx)?;
        let mut fixed = vec![0.0; dofs];
        let mut is_fixed = vec![false; dofs];
        let u: Vec<f64> = displacement_m.iter().flatten().copied().collect();
        let v: Vec<f64> = velocity_m_s.iter().flatten().copied().collect();
        // Native geometry validation above also checks support indices/duplicates.
        for &(dof, value) in problem.prescribed_m {
            if u[dof] != value || v[dof] != 0.0 {
                return Err(FiniteSolidError::Invalid(
                    "initial state must attain fixed displacements with zero supported velocity",
                ));
            }
            fixed[dof] = value;
            is_fixed[dof] = true;
        }
        let free: Vec<usize> = (0..dofs).filter(|&dof| !is_fixed[dof]).collect();
        if free.is_empty() {
            return Err(FiniteSolidError::Invalid("dynamic body requires a free DOF"));
        }
        let free_mass = reduce_matrix(&mass, &free, cx)?;
        if (0..free.len()).any(|row| {
            let value = free_mass.get(row, row);
            !value.is_finite() || value <= 0.0
        }) {
            return Err(FiniteSolidError::Invalid(
                "each free DOF needs finite positive reference mass; unused free nodes are refused",
            ));
        }
        let mut reference_mass_kg = 0.0;
        for row in 0..dofs {
            if row % 4096 == 0 { checkpoint(cx)?; }
            for value in mass.row(row).1 { reference_mass_kg += value / 3.0; }
        }
        if !reference_mass_kg.is_finite() || reference_mass_kg <= 0.0 {
            return Err(FiniteSolidError::Invalid("total reference mass is not representable"));
        }
        problem.stored_energy_j(&u, cx)?;
        let force = problem.external_force_n(cx)?;
        let q_free: Vec<f64> = free.iter().map(|&dof| u[dof]).collect();
        let v_free: Vec<f64> = free.iter().map(|&dof| v[dof]).collect();
        let forcing: Vec<f64> = free.iter().map(|&dof| force[dof]).collect();
        let fixed: Arc<[f64]> = Arc::from(fixed);
        let free: Arc<[usize]> = Arc::from(free);
        let free_mass = Arc::new(free_mass);
        let operator = DynamicOperator::new(&problem, &free_mass, &free, &fixed, cx);
        let preconditioner = MassDiagonal { mass: &free_mass };
        let result = second_order_acceleration(
            &operator, &q_free, &v_free, &forcing, &preconditioner, controls.solve,
            controls.maximum_workspace_components, &mut || cx.checkpoint().is_err(),
        );
        if let Some(error) = operator.error.take() { return Err(error.into()); }
        let initial = result?;
        drop(operator);
        let state = SecondOrderState::new(time_s, &q_free, &v_free, &initial.value);
        let body = Self {
            problem, mass: Arc::new(mass), free_mass, free, fixed,
            reference_mass_kg, initialization: initial.report, state,
        };
        body.kinetic_energy_j()?;
        checkpoint(cx)?;
        Ok(body)
    }

    /// Original geometry, sources and stationary supports. Load fields here
    /// describe initialization; each candidate explicitly supplies held loads.
    pub const fn problem(&self) -> &HyperTetProblem<'a, FiniteSolidMaterial> { &self.problem }
    pub const fn native_state(&self) -> &SecondOrderState { &self.state }
    pub const fn initialization_report(&self) -> &SolveReport { &self.initialization }
    pub fn free_dofs(&self) -> &[usize] { &self.free }
    pub const fn reference_mass_kg(&self) -> f64 { self.reference_mass_kg }

    /// Full node-major displacement, including fixed supports.
    pub fn displacement_m(&self) -> Vec<[f64; 3]> {
        inflate(&self.state.q, &self.free, &self.fixed).as_chunks::<3>().0.to_vec()
    }
    /// Full node-major velocity; fixed supports remain stationary.
    pub fn velocity_m_s(&self) -> Vec<[f64; 3]> {
        let mut velocity = vec![[0.0; 3]; self.problem.nodes_m.len()];
        for (&dof, value) in self.free.iter().zip(&self.state.v) {
            velocity[dof / 3][dof % 3] = *value;
        }
        velocity
    }
    /// Actual kinetic storage on the fixed reference mass.
    pub fn kinetic_energy_j(&self) -> Result<f64, FiniteSolidError> {
        let mut momentum = vec![0.0; self.free.len()];
        self.free_mass.spmv(&self.state.v, &mut momentum);
        let energy = 0.5 * dot(&self.state.v, &momentum);
        if !energy.is_finite() || energy < 0.0 {
            return Err(FiniteSolidError::Invalid("kinetic energy is not representable"));
        }
        Ok(energy)
    }

    /// Advance a candidate using held dead loads over a positive interval.
    ///
    /// The original body remains unchanged on every refusal, including a
    /// converged nonlinear solve that subsequently fails the physical energy
    /// gate. Accept the returned body explicitly to publish its q/v/a/history.
    pub fn step_candidate(
        &self,
        duration_s: f64,
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        controls: FiniteSolidControls,
        cx: &Cx<'_>,
    ) -> Result<(Self, FiniteSolidStep), FiniteSolidError> {
        checkpoint(cx)?;
        if !duration_s.is_finite() || duration_s <= 0.0
            || !controls.energy_tolerance_j.is_finite() || controls.energy_tolerance_j <= 0.0
            || controls.solve.max_newton_iterations == 0
            || self.state.steps >= controls.maximum_steps
        {
            return Err(FiniteSolidError::Invalid("invalid interval, controls or step budget"));
        }
        // Reborrow into this interval: transient loads never acquire the source
        // material's lifetime and are not retained in the accepted body.
        let loading = HyperTetProblem {
            nodes_m: self.problem.nodes_m,
            tetrahedra: self.problem.tetrahedra,
            material: self.problem.material,
            prescribed_m: self.problem.prescribed_m,
            nodal_forces_n,
            body_force_n_m3,
            budget: self.problem.budget,
            settings: self.problem.settings,
        };
        let force = loading.external_force_n(cx)?;
        let forcing: Vec<f64> = self.free.iter().map(|&dof| force[dof]).collect();
        let before_u = inflate(&self.state.q, &self.free, &self.fixed);
        let strain_before_j = self.problem.stored_energy_j(&before_u, cx)?;
        let kinetic_before_j = self.kinetic_energy_j()?;
        let operator = DynamicOperator::new(
            &self.problem, &self.free_mass, &self.free, &self.fixed, cx,
        );
        let integrator = OperatorGeneralizedAlpha::new(
            self.free.len(), duration_s, 1.0, controls.solve,
        );
        // Immutable matrices are shared; only the physical trajectory is copied.
        let mut candidate = self.clone();
        let result = integrator.step_controlled(
            &mut candidate.state, &operator, &forcing, &mut || cx.checkpoint().is_err(),
        );
        if let Some(error) = operator.error.take() { return Err(error.into()); }
        let time = result?;
        let after_u = inflate(&candidate.state.q, &self.free, &self.fixed);
        let strain_after_j = self.problem.stored_energy_j(&after_u, cx)?;
        let kinetic_after_j = candidate.kinetic_energy_j()?;
        let mut source_work_j = 0.0;
        for (dof, ((before, after), load)) in
            before_u.iter().zip(&after_u).zip(&force).enumerate()
        {
            if dof % 4096 == 0 { checkpoint(cx)?; }
            source_work_j = load.mul_add(after - before, source_work_j);
        }
        let integration_defect_j =
            (kinetic_after_j - kinetic_before_j) + (strain_after_j - strain_before_j)
            - source_work_j;
        if !source_work_j.is_finite() || !integration_defect_j.is_finite() {
            return Err(FiniteSolidError::Invalid("mechanical work is not representable"));
        }
        if integration_defect_j.abs() > controls.energy_tolerance_j {
            return Err(FiniteSolidError::Energy {
                defect_j: integration_defect_j,
                allowance_j: controls.energy_tolerance_j,
            });
        }
        let midpoint: Vec<f64> = before_u.iter().zip(&after_u)
            .map(|(before, after)| 0.5 * before + 0.5 * after).collect();
        let internal = self.problem.internal_force(&midpoint, cx)?;
        let mut acceleration = vec![0.0; before_u.len()];
        for ((&dof, before), after) in
            self.free.iter().zip(&self.state.v).zip(&candidate.state.v)
        {
            acceleration[dof] = (after - before) / duration_s;
        }
        let mut inertial = vec![0.0; before_u.len()];
        self.mass.spmv(&acceleration, &mut inertial);
        if internal.iter().chain(&inertial).any(|value| !value.is_finite()) {
            return Err(FiniteSolidError::Invalid("midpoint support force is not representable"));
        }
        let support_reactions_n: Vec<(usize, f64)> = self.problem.prescribed_m.iter()
            .map(|&(dof, _)| (dof, inertial[dof] + internal[dof] - force[dof])).collect();
        if support_reactions_n.iter().any(|(_, value)| !value.is_finite()) {
            return Err(FiniteSolidError::Invalid("support reaction accumulation overflow"));
        }
        checkpoint(cx)?;
        Ok((candidate, FiniteSolidStep {
            time, kinetic_before_j, kinetic_after_j, strain_before_j, strain_after_j,
            source_work_j, integration_defect_j, support_reactions_n,
        }))
    }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), HyperTetError> {
    cx.checkpoint().map_err(|_| HyperTetError::Cancelled)
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).fold(0.0, |sum, (x, y)| x.mul_add(*y, sum))
}

fn inflate(reduced: &[f64], free: &[usize], fixed: &[f64]) -> Vec<f64> {
    let mut full = fixed.to_vec();
    for (&dof, &value) in free.iter().zip(reduced) { full[dof] = value; }
    full
}

/// Extract the original native matrix's free/free submatrix in canonical order.
fn reduce_matrix(full: &Csr, free: &[usize], cx: &Cx<'_>) -> Result<Csr, HyperTetError> {
    let mut row_of = vec![usize::MAX; full.nrows()];
    for (row, &dof) in free.iter().enumerate() { row_of[dof] = row; }
    let mut coo = Coo::new(free.len(), free.len());
    for (row, &dof) in free.iter().enumerate() {
        if row % 4096 == 0 { checkpoint(cx)?; }
        let (columns, values) = full.row(dof);
        for (&column, &value) in columns.iter().zip(values) {
            if !value.is_finite() {
                return Err(HyperTetError::InvalidInput {
                    what: "non-finite native matrix entry".into(),
                });
            }
            if row_of[column] != usize::MAX { coo.push(row, row_of[column], value); }
        }
    }
    checkpoint(cx)?;
    Ok(coo.assemble())
}

struct MassDiagonal<'a> { mass: &'a Csr }
impl FlexiblePreconditioner for MassDiagonal<'_> {
    fn apply(&self, _iteration: usize, residual: &[f64], output: &mut [f64]) {
        for (row, (out, value)) in output.iter_mut().zip(residual).enumerate() {
            *out = value / self.mass.get(row, row);
        }
    }
}

/// Borrowed numerical adapter. All nonlinear mechanics remain in hyper3.
struct DynamicOperator<'a, 'b, 'cx> {
    problem: &'b HyperTetProblem<'a, FiniteSolidMaterial>,
    mass: &'b Csr,
    free: &'b [usize],
    fixed: &'b [f64],
    cx: &'b Cx<'cx>,
    tangent: RefCell<Option<(Vec<f64>, Csr)>>,
    // Void native operator callbacks communicate a persistent error here.
    // Speculative validate_state refusals remain recoverable by line search.
    error: RefCell<Option<HyperTetError>>,
}

impl<'a, 'b, 'cx> DynamicOperator<'a, 'b, 'cx> {
    fn new(
        problem: &'b HyperTetProblem<'a, FiniteSolidMaterial>,
        mass: &'b Csr,
        free: &'b [usize],
        fixed: &'b [f64],
        cx: &'b Cx<'cx>,
    ) -> Self {
        Self {
            problem, mass, free, fixed, cx,
            tangent: RefCell::new(None), error: RefCell::new(None),
        }
    }

    fn callback_ready(&self, output: &mut [f64]) -> bool {
        if self.error.borrow().is_some() {
            output.fill(f64::NAN);
            return false;
        }
        if let Err(error) = checkpoint(self.cx) {
            self.refuse(error, output);
            return false;
        }
        true
    }

    fn refuse(&self, error: HyperTetError, output: &mut [f64]) {
        if self.error.borrow().is_none() { self.error.replace(Some(error)); }
        output.fill(f64::NAN);
    }

    fn ensure_tangent(&self, q: &[f64]) -> Result<(), HyperTetError> {
        let rebuild = self.tangent.borrow().as_ref()
            .is_none_or(|(at, _)| at.as_slice() != q);
        if rebuild {
            let full = self.problem.tangent_matrix(&inflate(q, self.free, self.fixed), self.cx)?;
            let reduced = reduce_matrix(&full, self.free, self.cx)?;
            self.tangent.replace(Some((q.to_vec(), reduced)));
        }
        Ok(())
    }

    fn finite_output(&self, output: &mut [f64]) {
        if output.iter().any(|value| !value.is_finite()) {
            self.refuse(HyperTetError::InvalidInput {
                what: "non-finite native dynamic operator output".into(),
            }, output);
        } else if let Err(error) = checkpoint(self.cx) {
            self.refuse(error, output);
        }
    }
}

impl SecondOrderProblem for DynamicOperator<'_, '_, '_> {
    fn dimension(&self) -> usize { self.free.len() }

    fn validate_state(&self, q: &[f64], v: &[f64]) -> Result<(), String> {
        if q.len() != self.free.len() || v.len() != self.free.len()
            || q.iter().chain(v).any(|value| !value.is_finite())
        {
            return Err("finite dimension-matched displacement/velocity required".into());
        }
        self.problem.stored_energy_j(&inflate(q, self.free, self.fixed), self.cx)
            .map(|_| ()).map_err(|error| error.to_string())
    }

    fn mass_apply(&self, input: &[f64], output: &mut [f64]) {
        if !self.callback_ready(output) { return; }
        self.mass.spmv(input, output);
        self.finite_output(output);
    }

    fn damping_apply(&self, _input: &[f64], output: &mut [f64]) {
        if !self.callback_ready(output) { return; }
        output.fill(0.0);
    }

    fn internal_force(&self, q: &[f64], output: &mut [f64]) {
        if !self.callback_ready(output) { return; }
        match self.problem.internal_force(&inflate(q, self.free, self.fixed), self.cx) {
            Ok(full) => {
                for (value, &dof) in output.iter_mut().zip(self.free) { *value = full[dof]; }
                self.finite_output(output);
            }
            Err(error) => self.refuse(error, output),
        }
    }

    fn tangent_apply(&self, q: &[f64], direction: &[f64], output: &mut [f64]) {
        if !self.callback_ready(output) { return; }
        if let Err(error) = self.ensure_tangent(q) {
            self.refuse(error, output);
            return;
        }
        if let Some((_, matrix)) = self.tangent.borrow().as_ref() {
            matrix.spmv(direction, output);
        } else {
            output.fill(f64::NAN);
        }
        self.finite_output(output);
    }

    fn preconditioner_apply(
        &self,
        q: &[f64],
        weights: SecondOrderOperatorWeights,
        _outer_iteration: usize,
        _inner_iteration: usize,
        residual: &[f64],
        output: &mut [f64],
    ) {
        if !self.callback_ready(output) { return; }
        if let Err(error) = self.ensure_tangent(q) {
            self.refuse(error, output);
            return;
        }
        let tangent = self.tangent.borrow();
        let Some((_, matrix)) = tangent.as_ref() else {
            output.fill(f64::NAN);
            return;
        };
        for (row, (out, value)) in output.iter_mut().zip(residual).enumerate() {
            let diagonal = weights.mass.mul_add(
                self.mass.get(row, row), weights.tangent * matrix.get(row, row),
            ).abs();
            *out = if diagonal > 0.0 { value / diagonal } else { *value };
        }
        self.finite_output(output);
    }
}
