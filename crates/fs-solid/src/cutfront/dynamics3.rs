//! Material-backed CutFEM motion through the original generalized-alpha owner.
//!
//! The retained native stiffness and consistent mass provide M a + K q = f.
//! A held eigenstrain changes resistance to K q - b_free using its original
//! native equivalent force; caller-supplied forcing remains external only.
//! Optional frozen initial stress adds G(q-q_0), prepared from that same held
//! eigenstrain and an admitted reference displacement. Its signed tangent is
//! K+G; the geometric force vanishes at q_0 without changing physical inertia.
//! Strong homogeneous supports are removed once from every action and load;
//! adaptive hanging constraints remain inside the original native operators.
//! No mass lumping, stabilization inertia, constitutive fit or new time scheme
//! is introduced. Geometry, occupancy and stiffness scales remain fixed. Each
//! immutable model holds one source and free strain; thermal_midpoint_candidate
//! pairs accepted and candidate source endpoints through the original midpoint
//! owner. Evolving prestress, nonzero/moving supports, follower loads and
//! physical damping are outside this elastic adapter.
//!
//! Native mass and strain callbacks poll cancellation within their kernels.
//! The original stiffness action remains indivisible; cancellation is checked
//! before and after that action on the admitted, bounded native grid.

use std::cell::RefCell;
use std::ops::ControlFlow;

use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
use fs_cutfem::elastic3::eigenstrain::PreparedEigenstrain3;
use fs_cutfem::elastic3::mass::MaterialMass3;
use fs_cutfem::elastic3::{CutElasticity3, ElasticityError3};
use fs_exec::CancelGate;
use fs_solver::LinearOp;
use fs_time::galpha::second_order_adjoint::SecondOrderVjp;
use fs_time::galpha::{
    ImplicitStepTelemetry, OperatorGeneralizedAlpha, SecondOrderOperatorWeights,
    SecondOrderProblem, SecondOrderState, TimeSolveError,
};

mod prestress;
mod thermal;
use prestress::HeldPrestress3;
pub use thermal::{CutThermalMidpointSettings3, CutThermalMidpointStep3};

/// The original fixed-geometry elasticity owner; this view creates no mesh.
#[derive(Clone, Copy)]
pub enum CutElasticityRef3<'a> {
    Cartesian(&'a CutElasticity3),
    Adaptive(&'a AdaptiveElasticity3),
}

impl<'a> CutElasticityRef3<'a> {
    pub(super) fn operator(self) -> &'a dyn LinearOp {
        match self {
            Self::Cartesian(op) => op,
            Self::Adaptive(op) => op,
        }
    }

    pub(super) fn fixed(self) -> &'a [bool] {
        match self {
            Self::Cartesian(op) => op.fixed(),
            Self::Adaptive(op) => op.fixed(),
        }
    }

    pub(super) fn validate_mass(self, mass: &MaterialMass3) -> Result<(), ElasticityError3> {
        match self {
            Self::Cartesian(op) => op.validate_material_mass(mass),
            Self::Adaptive(op) => op.validate_material_mass(mass),
        }
    }

    fn admit(self, q: &[f64], gate: &CancelGate) -> Result<(), ElasticityError3> {
        match self {
            Self::Cartesian(op) => op.admit_material_domain(q, || checkpoint(gate)),
            Self::Adaptive(op) => op.admit_material_domain(q, || checkpoint(gate)),
        }
    }

    fn validate_eigenstrain(self, prepared: &PreparedEigenstrain3) -> Result<(), ElasticityError3> {
        match self {
            Self::Cartesian(op) => op.validate_eigenstrain(prepared),
            Self::Adaptive(op) => op.validate_eigenstrain(prepared),
        }
    }

    fn eigenstrain_load(
        self,
        prepared: &PreparedEigenstrain3,
        gate: &CancelGate,
    ) -> Result<Vec<f64>, ElasticityError3> {
        match self {
            Self::Cartesian(op) => op.eigenstrain_load(prepared, || checkpoint(gate)),
            Self::Adaptive(op) => op.eigenstrain_load(prepared, || checkpoint(gate)),
        }
    }

    fn admit_eigenstrain(
        self,
        held: &HeldEigenstrain3<'_>,
        q: &[f64],
        gate: &CancelGate,
    ) -> Result<(), ElasticityError3> {
        match self {
            Self::Cartesian(op) => op.admit_eigenstrain(
                held.prepared,
                q,
                held.maximum_kinematic_strain,
                held.maximum_points,
                || checkpoint(gate),
            ),
            Self::Adaptive(op) => op.admit_eigenstrain(
                held.prepared,
                q,
                held.maximum_kinematic_strain,
                held.maximum_points,
                || checkpoint(gate),
            ),
        }
    }
}

fn checkpoint(gate: &CancelGate) -> ControlFlow<()> {
    if gate.is_requested() {
        ControlFlow::Break(())
    } else {
        ControlFlow::Continue(())
    }
}

fn poll(gate: &CancelGate) -> Result<(), ElasticityError3> {
    if checkpoint(gate).is_break() {
        Err(ElasticityError3::Cancelled)
    } else {
        Ok(())
    }
}

fn invalid(message: &'static str) -> ElasticityError3 {
    ElasticityError3::Invalid(message)
}

struct Scratch {
    input: Vec<f64>,
    output: Vec<f64>,
    /// Allocated only when a frozen reference prestress is selected.
    geometric: Vec<f64>,
}

#[derive(Clone, Copy)]
enum Action {
    Mass,
    Tangent,
    TangentTranspose,
    Resistance,
}

struct HeldEigenstrain3<'a> {
    prepared: &'a PreparedEigenstrain3,
    /// Already in the same free coordinates as M and K, assembled only once.
    force: Vec<f64>,
    maximum_kinematic_strain: f64,
    maximum_points: usize,
}

/// Fixed material pencil in free native coordinates, usable by structural
/// stepping, consistent acceleration initialization and discrete adjoints.
///
/// Source identity is checked by the native mass owner, not by shape equality.
/// Physical displacement admission uses the original material strain ceilings
/// at every generalized-alpha stage and endpoint. A held eigenstrain uses
/// elastic strain for those ceilings and a separate total-kinematic budget,
/// evaluated at the native owner's retained quadrature points. All actions retain
/// O(native DOFs) scratch; no dense pencil or copied material field is prepared.
/// Positive mass diagonals are necessary admission, not a nonsingularity proof.
/// Native solver residuals and spatial/model validity remain separate claims.
pub struct CutSecondOrderProblem3<'a> {
    source: CutElasticityRef3<'a>,
    mass: &'a MaterialMass3,
    occupancy: &'a [f64],
    gate: &'a CancelGate,
    free: Vec<usize>,
    mass_diagonal: Vec<f64>,
    eigenstrain: Option<HeldEigenstrain3<'a>>,
    prestress: Option<HeldPrestress3>,
    scratch: RefCell<Scratch>,
}

impl<'a> CutSecondOrderProblem3<'a> {
    /// Bind one exact native pencil and its independent physical cell fractions.
    /// `maximum_native_dofs` bounds the original coordinate scratch before
    /// allocation. Fractions are finite in [0,1]; stiffness scales are never
    /// reused as mass fractions. Zero-inertia free rows and empty spaces refuse.
    pub fn new(
        source: CutElasticityRef3<'a>,
        mass: &'a MaterialMass3,
        occupancy: &'a [f64],
        maximum_native_dofs: usize,
        gate: &'a CancelGate,
    ) -> Result<Self, ElasticityError3> {
        poll(gate)?;
        source.validate_mass(mass)?;
        let n = source.operator().n();
        if n == 0 || n > maximum_native_dofs {
            return Err(invalid("CutFEM dynamics exceeds the native-DOF allowance"));
        }
        let diagonal = mass.diagonal(occupancy, || checkpoint(gate))?;
        let mut free = Vec::new();
        let mut mass_diagonal = Vec::new();
        for (dof, &value) in diagonal.iter().enumerate() {
            if dof % 192 == 0 {
                poll(gate)?;
            }
            if !source.fixed()[dof / 3] {
                if value <= 0.0 || !value.is_finite() {
                    return Err(invalid(
                        "CutFEM dynamics requires inertia on every free coordinate",
                    ));
                }
                free.push(dof);
                mass_diagonal.push(value);
            }
        }
        if free.is_empty() {
            return Err(invalid("CutFEM dynamics requires a free coordinate"));
        }
        poll(gate)?;
        Ok(Self {
            source,
            mass,
            occupancy,
            gate,
            free,
            mass_diagonal,
            eigenstrain: None,
            prestress: None,
            scratch: RefCell::new(Scratch {
                input: vec![0.0; n],
                output: vec![0.0; n],
                geometric: Vec::new(),
            }),
        })
    }

    /// Bind a prescribed thermal or transformation strain for the entire
    /// trajectory. Native preparation must belong to this exact geometry,
    /// material assignment, supports and displacement space; equal dimensions
    /// are insufficient. The immutable source borrow also freezes its scales.
    /// Its equivalent force is assembled at those scales and retained once in
    /// free coordinates, independently of the physical mass fractions.
    ///
    /// The equation is `M a + (K q - b_free) = f_external`. Supply EXTERNAL
    /// forcing only to [`Self::advance`], consistent acceleration initialization
    /// and the existing discrete adjoints; adding `b_free` again doubles it.
    /// Tangents and state/force VJPs hold free strain fixed. This constructor
    /// neither advances temperature nor supplies a thermal-parameter gradient.
    ///
    /// `maximum_kinematic_strain` is the finite total-strain Frobenius budget
    /// in (0,1); the original material separately admits actual and B-bar
    /// elastic strains. Both checks use the native retained bulk points, not
    /// certified field extrema. `maximum_eigenstrain_points` bounds this work
    /// at every physical stage and endpoint without allocating stress outputs.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_eigenstrain(
        source: CutElasticityRef3<'a>,
        mass: &'a MaterialMass3,
        occupancy: &'a [f64],
        eigenstrain: &'a PreparedEigenstrain3,
        maximum_kinematic_strain: f64,
        maximum_eigenstrain_points: usize,
        maximum_native_dofs: usize,
        gate: &'a CancelGate,
    ) -> Result<Self, ElasticityError3> {
        poll(gate)?;
        source.validate_eigenstrain(eigenstrain)?;
        if !maximum_kinematic_strain.is_finite()
            || maximum_kinematic_strain <= 0.0
            || maximum_kinematic_strain >= 1.0
        {
            return Err(invalid("invalid CutFEM dynamics total-strain budget"));
        }
        if eigenstrain.bulk_points() > maximum_eigenstrain_points {
            return Err(invalid(
                "CutFEM dynamics exceeds the eigenstrain-point allowance",
            ));
        }
        let mut model = Self::new(source, mass, occupancy, maximum_native_dofs, gate)?;
        let native_force = source.eigenstrain_load(eigenstrain, gate)?;
        let force = model.restrict_force(&native_force)?;
        poll(gate)?;
        model.eigenstrain = Some(HeldEigenstrain3 {
            prepared: eigenstrain,
            force,
            maximum_kinematic_strain,
            maximum_points: maximum_eigenstrain_points,
        });
        Ok(model)
    }

    /// Sorted independent native DOFs. Adaptive entries are master coordinates.
    pub fn free_dofs(&self) -> &[usize] {
        &self.free
    }

    /// Expand an algebraic free field into the original operator coordinates.
    /// Fixed entries are exactly zero. Adaptive physical-node reconstruction
    /// remains `AdaptiveElasticity3::physical_displacements` on this result.
    pub fn expand_field(&self, field: &[f64]) -> Result<Vec<f64>, ElasticityError3> {
        self.admit_vector(field, self.free.len())?;
        let mut full = vec![0.0; self.mass.n()];
        self.scatter(field, &mut full)?;
        Ok(full)
    }

    /// Restrict a native force [N] with the same virtual-work transpose as M/K.
    /// Loads on stationary fixed DOFs do no work; callers retain those original
    /// entries if they need support reactions, which this adapter does not infer.
    /// Adaptive forces must already be condensed by their native load owner.
    pub fn restrict_force(&self, force: &[f64]) -> Result<Vec<f64>, ElasticityError3> {
        self.admit_vector(force, self.mass.n())?;
        Ok(self.free.iter().map(|&dof| force[dof]).collect())
    }

    /// Advance using the existing generalized-alpha method and this same gate.
    /// `native_stage_force_n` contains external forces only, evaluated at the
    /// method's force stage time `t_n+(1-alpha_f)h`. A configured eigenstrain
    /// already contributes through internal resistance and must not be added
    /// to this input. The original integrator owns convergence, physical-state
    /// admission and atomic publication; no second accepted clock is kept here.
    pub fn advance(
        &self,
        method: &OperatorGeneralizedAlpha,
        state: &mut SecondOrderState,
        native_stage_force_n: &[f64],
    ) -> Result<ImplicitStepTelemetry, TimeSolveError> {
        let force = self
            .restrict_force(native_stage_force_n)
            .map_err(|error| match error {
                ElasticityError3::Cancelled => TimeSolveError::Cancelled,
                other => TimeSolveError::ModelRefused {
                    phase: "native stage force",
                    reason: other.to_string(),
                },
            })?;
        method.step_controlled(state, self, &force, &mut || self.gate.is_requested())
    }

    fn admit_vector(&self, values: &[f64], n: usize) -> Result<(), ElasticityError3> {
        poll(self.gate)?;
        if values.len() != n {
            return Err(invalid("CutFEM dynamics field dimension mismatch"));
        }
        for (i, value) in values.iter().enumerate() {
            if i % 192 == 0 {
                poll(self.gate)?;
            }
            if !value.is_finite() {
                return Err(invalid("nonfinite CutFEM dynamics field"));
            }
        }
        Ok(())
    }

    fn scatter(&self, reduced: &[f64], full: &mut [f64]) -> Result<(), ElasticityError3> {
        for (i, (&dof, &value)) in self.free.iter().zip(reduced).enumerate() {
            if i % 192 == 0 {
                poll(self.gate)?;
            }
            full[dof] = value;
        }
        poll(self.gate)
    }

    fn action(
        &self,
        input: &[f64],
        output: &mut [f64],
        action: Action,
    ) -> Result<(), ElasticityError3> {
        self.admit_vector(input, self.free.len())?;
        if output.len() != self.free.len() {
            return Err(invalid("CutFEM dynamics output dimension mismatch"));
        }
        let mut scratch = self.scratch.borrow_mut();
        let Scratch {
            input: full,
            output: result,
            geometric,
        } = &mut *scratch;
        self.scatter(input, full)?;
        match action {
            Action::Mass => self
                .mass
                .apply_into(self.occupancy, full, result, || checkpoint(self.gate))?,
            Action::TangentTranspose => self.source.operator().apply_transpose(full, result),
            Action::Tangent | Action::Resistance => self.source.operator().apply(full, result),
        }
        let resistance = matches!(action, Action::Resistance);
        let prestress = if matches!(action, Action::Mass) {
            None
        } else {
            self.prestress.as_ref()
        };
        if let Some(held) = prestress {
            if resistance {
                // Center before multiplication: G(q-q_0), never Gq-Gq_0.
                // Strong rows stay zero because the reference was admitted.
                for (i, (value, reference)) in full.iter_mut().zip(&held.reference).enumerate() {
                    if i % 192 == 0 {
                        poll(self.gate)?;
                    }
                    *value -= reference;
                }
            }
            // The native initial-stress operator is symmetric, including its
            // selected Nitsche flux and exact adaptive constraint transpose.
            held.prepared
                .apply_into(full, geometric, || checkpoint(self.gate))?;
        }
        poll(self.gate)?;
        for (i, (&dof, value)) in self.free.iter().zip(output.iter_mut()).enumerate() {
            if i % 192 == 0 {
                poll(self.gate)?;
            }
            *value = result[dof];
            // Resolve the thermal offset before adding a weak geometric force.
            if resistance && let Some(held) = &self.eigenstrain {
                *value -= held.force[i];
            }
            if prestress.is_some() {
                *value += geometric[dof];
            }
            if !value.is_finite() {
                return Err(invalid("CutFEM dynamics operator overflow"));
            }
        }
        poll(self.gate)
    }
}

impl SecondOrderProblem for CutSecondOrderProblem3<'_> {
    fn dimension(&self) -> usize {
        self.free.len()
    }

    fn validate_state(&self, q: &[f64], v: &[f64]) -> Result<(), String> {
        let admitted = || -> Result<(), ElasticityError3> {
            self.admit_vector(q, self.free.len())?;
            self.admit_vector(v, self.free.len())?;
            let mut scratch = self.scratch.borrow_mut();
            self.scatter(q, &mut scratch.input)?;
            match &self.eigenstrain {
                Some(held) => self
                    .source
                    .admit_eigenstrain(held, &scratch.input, self.gate),
                None => self.source.admit(&scratch.input, self.gate),
            }
        };
        admitted().map_err(|error| error.to_string())
    }

    fn mass_apply(&self, input: &[f64], output: &mut [f64]) {
        if self.action(input, output, Action::Mass).is_err() {
            output.fill(f64::NAN);
        }
    }

    fn damping_apply(&self, input: &[f64], output: &mut [f64]) {
        if self.admit_vector(input, self.free.len()).is_err() || output.len() != self.free.len() {
            output.fill(f64::NAN);
        } else {
            output.fill(0.0);
        }
    }

    fn internal_force(&self, q: &[f64], output: &mut [f64]) {
        if self.action(q, output, Action::Resistance).is_err() {
            output.fill(f64::NAN);
        }
    }

    fn tangent_apply(&self, _: &[f64], direction: &[f64], output: &mut [f64]) {
        // Both held free strain and the geometric reference are affine offsets.
        if self.action(direction, output, Action::Tangent).is_err() {
            output.fill(f64::NAN);
        }
    }

    fn resistance_balance(
        &self,
        q: &[f64],
        dq: &[f64],
        fraction: f64,
        v: &[f64],
        force: &[f64],
        output: &mut [f64],
    ) {
        fs_time::galpha::linear_resistance_balance(self, q, dq, fraction, v, force, output);
    }

    fn preconditioner_apply(
        &self,
        _: &[f64],
        weights: SecondOrderOperatorWeights,
        _: usize,
        _: usize,
        residual: &[f64],
        output: &mut [f64],
    ) {
        // A diagonal preconditioner only; the physical mass action stays consistent.
        for ((value, rhs), diagonal) in output.iter_mut().zip(residual).zip(&self.mass_diagonal) {
            *value = rhs / (weights.mass * diagonal);
        }
        if self.gate.is_requested() {
            output.fill(f64::NAN);
        }
    }
}

/// State/force pullbacks at this fixed physical configuration. No design,
/// geometry, occupancy, free-strain, reference-prestress or material parameter
/// derivative is claimed. A selected geometric tangent remains frozen.
impl SecondOrderVjp for CutSecondOrderProblem3<'_> {
    fn parameter_count(&self) -> usize {
        0
    }

    fn mass_transpose_apply(&self, seed: &[f64], output: &mut [f64]) -> Result<(), String> {
        self.action(seed, output, Action::Mass)
            .map_err(|error| error.to_string())
    }

    fn damping_transpose_apply(&self, seed: &[f64], output: &mut [f64]) -> Result<(), String> {
        self.admit_vector(seed, self.free.len())
            .map_err(|error| error.to_string())?;
        if output.len() != self.free.len() {
            return Err("CutFEM dynamics output dimension mismatch".into());
        }
        output.fill(0.0);
        Ok(())
    }

    fn tangent_transpose_apply(
        &self,
        _: &[f64],
        seed: &[f64],
        output: &mut [f64],
    ) -> Result<(), String> {
        self.action(seed, output, Action::TangentTranspose)
            .map_err(|error| error.to_string())
    }

    fn residual_parameter_vjp(
        &self,
        _: &[f64],
        _: &[f64],
        _: &[f64],
        _: &[f64],
        parameters: &mut [f64],
    ) -> Result<(), String> {
        poll(self.gate).map_err(|error| error.to_string())?;
        if parameters.is_empty() {
            Ok(())
        } else {
            Err("CutFEM dynamics has no varying model parameters".into())
        }
    }
}
