//! Physics/IR boundary tests. Analytic providers isolate composition and refusal
//! behavior; mesh-based primal/adjoint integration is exercised by fs-ascent.
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_opt::reverse::physics::{PhysicsBinding, PhysicsError, PhysicsModel, PhysicsSample, PhysicsSignature};
use fs_opt::reverse::{HessianError, ReverseError, ReverseLimits, ReverseProgram};
use fs_opt::{ConstraintKind, Manifold, NodeId, OptError, Problem, ProblemBuilder, ReverseProblem, Sense};
use fs_qty::Dims;
use std::cell::Cell;

const LIMITS: ReverseLimits = ReverseLimits { max_nodes: 256, max_scalar_slots: 4096 };

#[derive(Debug)]
struct Model<'a> {
    name: &'a str,
    dimension: usize,
    units: Dims,
    calls: Cell<usize>,
    failure: Option<PhysicsError>,
    malformed: u8,
    cancel_after: Option<&'a CancelGate>,
}
impl Model<'_> {
    fn normal() -> Self {
        Self { name: "quadratic-study", dimension: 2, units: Dims::NONE,
            calls: Cell::new(0), failure: None, malformed: 0, cancel_after: None }
    }
}
impl PhysicsModel for Model<'_> {
    fn signature(&self) -> PhysicsSignature<'_> {
        PhysicsSignature { study: self.name, parameter_count: self.dimension,
            parameter_dims: self.units, residual_dims: Dims::NONE }
    }
    fn value_gradient(&self, x: &[f64], _cx: Option<&Cx<'_>>) -> Result<PhysicsSample, PhysicsError> {
        self.calls.set(self.calls.get() + 1);
        if let Some(error) = &self.failure { return Err(error.clone()); }
        let mut sample = PhysicsSample { value: x[0]*x[0] + 3.0*x[1], gradient: vec![2.0*x[0], 3.0] };
        match self.malformed {
            1 => { sample.gradient.pop(); }
            2 => sample.value = f64::INFINITY,
            3 => sample.gradient[1] = f64::NAN,
            _ => {}
        }
        if let Some(gate) = self.cancel_after { gate.request(); }
        Ok(sample)
    }
}

fn fixture(adjoint: bool) -> (Problem, NodeId, NodeId) {
    let mut b = ProblemBuilder::new();
    let x = b.var("x", Manifold::Rn { dim: 2 }, Dims::NONE).unwrap();
    let xr = b.var_ref(x).unwrap();
    let norm = b.norm_sq(xr).unwrap();
    let pde = b.pde_residual("quadratic-study", x, adjoint, Dims::NONE).unwrap();
    let same = b.pde_residual("quadratic-study", x, adjoint, Dims::NONE).unwrap();
    assert_eq!(same, pde, "one shared physics expression");
    let square = b.powi(pde, 2).unwrap();
    let half = b.konst(0.5, Dims::NONE).unwrap();
    let regularizer = b.mul(half, norm).unwrap();
    let total = b.add(square, regularizer).unwrap();
    let two = b.konst(2.0, Dims::NONE).unwrap();
    let residual = b.sub(pde, two).unwrap();
    b.objective(total, Sense::Minimize, 1.0).unwrap();
    b.constraint(residual, ConstraintKind::LeZero, "study cap").unwrap();
    (b.finish(), pde, norm)
}

#[test]
fn physical_gradient_composes_with_algebra_and_constraint_pullbacks_once() {
    let (problem, node, _) = fixture(true);
    let model = Model::normal();
    let oracle = ReverseProblem::new_with_physics(&problem, LIMITS,
        &[PhysicsBinding { node, model: &model }], None).unwrap();
    let mut input = vec![vec![2.0, -1.0]];
    let tape = oracle.evaluate(&input).unwrap();
    input[0][0] = 100.0;
    assert_eq!(tape.objective_value(), 3.5);
    assert_eq!(tape.constraint_values(), &[-1.0]);
    assert_eq!(tape.objective_gradient().unwrap(), vec![vec![10.0, 5.0]]);
    assert_eq!(tape.constraint_pullback(&[2.0]).unwrap(), vec![vec![8.0, 6.0]]);
    assert_eq!(tape.lagrangian_gradient(&[2.0]).unwrap(), vec![vec![18.0, 11.0]]);
    assert_eq!(tape.objective_gradient().unwrap(), vec![vec![10.0, 5.0]]);
    assert_eq!(model.calls.get(), 1, "pullbacks must not rerun physics");
}

#[test]
fn duplicate_roots_accumulate_and_separate_models_bind_exact_variables() {
    let (problem, node, _) = fixture(true);
    let model = Model::normal();
    let program = ReverseProgram::new_with_physics(&problem, &[node, node], LIMITS,
        &[PhysicsBinding { node, model: &model }]).unwrap();
    let tape = program.evaluate(&[vec![2.0, -1.0]], None).unwrap();
    assert_eq!(tape.pullback(&[1.0, 2.0], None).unwrap(), vec![vec![12.0, 9.0]]);
    assert_eq!(model.calls.get(), 1);

    let mut b = ProblemBuilder::new();
    let x = b.var("x", Manifold::Rn { dim: 2 }, Dims::NONE).unwrap();
    let y = b.var("y", Manifold::Rn { dim: 2 }, Dims::NONE).unwrap();
    let px = b.pde_residual("quadratic-study", x, true, Dims::NONE).unwrap();
    let py = b.pde_residual("quadratic-study", y, true, Dims::NONE).unwrap();
    let problem = b.finish();
    let program = ReverseProgram::new_with_physics(&problem, &[px, py], LIMITS,
        &[PhysicsBinding { node: py, model: &model }, PhysicsBinding { node: px, model: &model }]).unwrap();
    let tape = program.evaluate(&[vec![2.0, -1.0], vec![-1.0, 4.0]], None).unwrap();
    assert_eq!(tape.values(), &[1.0, 13.0]);
    assert_eq!(tape.pullback(&[2.0, -1.0], None).unwrap(), vec![vec![8.0, 6.0], vec![2.0, -3.0]]);
    assert_eq!(model.calls.get(), 3);
}

#[test]
fn default_unbound_and_unreachable_semantics_are_preserved() {
    let (problem, node, norm) = fixture(true);
    assert!(matches!(ReverseProgram::new(&problem, &[node], LIMITS),
        Err(ReverseError::Evaluation(OptError::Unevaluable { .. }))));
    let pure = ReverseProgram::new(&problem, &[norm], LIMITS).unwrap();
    let model = Model::normal();
    let unused = ReverseProgram::new_with_physics(&problem, &[norm], LIMITS,
        &[PhysicsBinding { node, model: &model }]).unwrap();
    assert_eq!(unused.scalar_slots(), pure.scalar_slots());
    let input = [vec![2.0, -1.0]];
    let a = pure.evaluate(&input, None).unwrap();
    let b = unused.evaluate(&input, None).unwrap();
    assert_eq!(a.values()[0].to_bits(), b.values()[0].to_bits());
    assert_eq!(a.pullback(&[1.0], None).unwrap(), b.pullback(&[1.0], None).unwrap());
    assert_eq!(model.calls.get(), 0);
}

#[test]
fn signatures_duplicates_and_unavailable_adjoint_refuse_before_execution() {
    let (problem, node, norm) = fixture(true);
    for case in 0..3 {
        let mut model = Model::normal();
        match case {
            0 => model.name = "different-scientific-study",
            1 => model.dimension = 3,
            _ => model.units = Dims([1, 0, 0, 0, 0, 0]),
        }
        assert!(matches!(ReverseProgram::new_with_physics(&problem, &[node], LIMITS,
            &[PhysicsBinding { node, model: &model }]), Err(ReverseError::PhysicsBinding { .. })));
        assert_eq!(model.calls.get(), 0);
    }
    let model = Model::normal();
    let binding = PhysicsBinding { node, model: &model };
    assert!(matches!(ReverseProgram::new_with_physics(&problem, &[node], LIMITS,
        &[binding, binding]), Err(ReverseError::PhysicsBinding { what: "duplicate physics binding", .. })));
    assert!(ReverseProgram::new_with_physics(&problem, &[norm], LIMITS,
        &[PhysicsBinding { node: norm, model: &model }]).is_err());
    assert!(ReverseProgram::new_with_physics(&problem, &[norm], LIMITS,
        &[PhysicsBinding { node: NodeId(u32::MAX), model: &model }]).is_err());
    let (problem, node, _) = fixture(false);
    assert!(matches!(ReverseProgram::new_with_physics(&problem, &[node], LIMITS,
        &[PhysicsBinding { node, model: &model }]), Err(ReverseError::PhysicsBinding { .. })));
}

#[test]
fn cached_physics_gradients_are_in_the_scalar_storage_cap() {
    let mut b = ProblemBuilder::new();
    let x = b.var("x", Manifold::Rn { dim: 2 }, Dims::NONE).unwrap();
    let node = b.pde_residual("quadratic-study", x, true, Dims::NONE).unwrap();
    let problem = b.finish();
    let model = Model::normal();
    for max_scalar_slots in [3, 4] {
        assert!(matches!(ReverseProgram::new_with_physics(&problem, &[node],
            ReverseLimits { max_scalar_slots, ..LIMITS }, &[PhysicsBinding { node, model: &model }]),
            Err(ReverseError::Evaluation(OptError::CapExceeded { .. }))));
    }
    let exact = ReverseProgram::new_with_physics(&problem, &[node],
        ReverseLimits { max_scalar_slots: 5, ..LIMITS }, &[PhysicsBinding { node, model: &model }]).unwrap();
    assert_eq!(exact.scalar_slots(), 5); // two inputs + scalar value + two derivatives
    assert_eq!(model.calls.get(), 0);
}

#[test]
fn malformed_or_nonfinite_provider_outputs_never_publish_a_tape() {
    let (problem, node, _) = fixture(true);
    for malformed in 1..=3 {
        let mut model = Model::normal(); model.malformed = malformed;
        let program = ReverseProgram::new_with_physics(&problem, &[node], LIMITS,
            &[PhysicsBinding { node, model: &model }]).unwrap();
        let error = program.evaluate(&[vec![2.0, -1.0]], None).unwrap_err();
        match malformed {
            1 => assert!(matches!(error, ReverseError::PhysicsGradientLength { expected: 2, actual: 1, .. })),
            2 => assert!(matches!(error, ReverseError::Evaluation(OptError::EvalNonFinite { .. }))),
            _ => assert!(matches!(error, ReverseError::NonFiniteAdjoint { component: 1, .. })),
        }
    }
}

#[test]
fn original_budget_solver_and_domain_errors_keep_node_attribution() {
    let (problem, node, _) = fixture(true);
    for failure in [
        PhysicsError::Budget { resource: "matvecs", used: 10, limit: 10 },
        PhysicsError::NotConverged { phase: "adjoint", iterations: 7,
            residual_bits: 0.1f64.to_bits(), tolerance_bits: 1e-8f64.to_bits() },
        PhysicsError::Domain("coefficient must be positive"),
    ] {
        let mut model = Model::normal(); model.failure = Some(failure.clone());
        let program = ReverseProgram::new_with_physics(&problem, &[node], LIMITS,
            &[PhysicsBinding { node, model: &model }]).unwrap();
        assert_eq!(program.evaluate(&[vec![2.0, -1.0]], None).unwrap_err(),
            ReverseError::Physics { node, source: failure });
    }
}

#[test]
fn first_order_physics_does_not_fabricate_a_second_order_rule() {
    let (problem, node, norm) = fixture(true);
    let model = Model::normal();
    let program = ReverseProgram::new_with_physics(&problem, &[norm, node], LIMITS,
        &[PhysicsBinding { node, model: &model }]).unwrap();
    let tape = program.evaluate(&[vec![2.0, -1.0]], None).unwrap();
    assert!(matches!(tape.hessian_vector_product(&[1.0, 1.0], &[vec![1.0, 2.0]], None),
        Err(HessianError::Reverse(ReverseError::Evaluation(OptError::Unevaluable { node: id, .. }))) if id == node.0));
    assert_eq!(tape.hessian_vector_product(&[1.0, 0.0], &[vec![1.0, 2.0]], None).unwrap(), vec![vec![2.0, 4.0]]);
    assert_eq!(model.calls.get(), 1);
}

#[test]
fn cancellation_brackets_provider_execution_and_preserves_original_channel() {
    let (problem, node, _) = fixture(true);
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 0, kernel_id: 1, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        let mut model = Model::normal(); model.cancel_after = Some(&gate);
        let program = ReverseProgram::new_with_physics(&problem, &[node], LIMITS,
            &[PhysicsBinding { node, model: &model }]).unwrap();
        assert_eq!(program.evaluate(&[vec![2.0, -1.0]], Some(&cx)).unwrap_err(),
            ReverseError::Evaluation(OptError::Cancelled));
        assert_eq!(model.calls.get(), 1);
        assert_eq!(program.evaluate(&[vec![2.0, -1.0]], Some(&cx)).unwrap_err(),
            ReverseError::Evaluation(OptError::Cancelled));
        assert_eq!(model.calls.get(), 1);
    });
    let mut model = Model::normal(); model.failure = Some(PhysicsError::Cancelled);
    let program = ReverseProgram::new_with_physics(&problem, &[node], LIMITS,
        &[PhysicsBinding { node, model: &model }]).unwrap();
    assert_eq!(program.evaluate(&[vec![2.0, -1.0]], None).unwrap_err(),
        ReverseError::Evaluation(OptError::Cancelled));
}
