//! G0/G1/G3: held thermal strain, consistent inertia and the original time adjoint.
use std::ops::ControlFlow;

use fs_blake3::ContentHash;
use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
use fs_cutfem::elastic3::{CutElasticity3, ElasticMaterialField3, ElasticityError3};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_exec::CancelGate;
use fs_ivl::Interval;
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_solid::cutfront::dynamics3::{CutElasticityRef3, CutSecondOrderProblem3};
use fs_solver::{FlexiblePreconditioner, LinearOp, NewtonKrylovConfig};
use fs_time::galpha::initialization::{
    InitialSolveConfig, second_order_acceleration, second_order_acceleration_vjp,
};
use fs_time::galpha::second_order_adjoint::{SecondOrderAdjointConfig, SecondOrderVjp};
use fs_time::galpha::{
    ImplicitSolveConfig, OperatorGeneralizedAlpha, SecondOrderProblem, SecondOrderState,
    TimeSolveError,
};

const EXPANSION: f64 = 0.004;
const POINTS: usize = 100_000;

#[path = "cutfem_eigenstrain_dynamics/prestress.rs"]
mod prestress;

struct Top;
impl CutSdf3 for Top {
    fn value(&self, p: [f64; 3]) -> f64 {
        p[2] - 0.73
    }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        Interval::new(lo[2] - 0.73, hi[2] - 0.73)
    }
    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], axis: HeightAxis) -> Interval {
        let value = if axis == HeightAxis::Z { 1.0 } else { 0.0 };
        Interval::new(value, value)
    }
}
fn poll() -> ControlFlow<()> {
    ControlFlow::Continue(())
}
fn material(limit: f64) -> ElasticAcousticState {
    let mut tensor = std::array::from_fn(|i| {
        std::array::from_fn(|j| if i == j { 15.0 + i as f64 } else { 0.0 })
    });
    tensor[0][3] = 1.2;
    tensor[3][0] = 1.2;
    ElasticAcousticState::try_new_mandel(tensor, 3.0, ContentHash([42; 32]), limit).unwrap()
}
fn control<T>(action: impl FnOnce(&mut QuadratureControl3<'_>) -> T) -> T {
    let mut callback = |_| poll();
    let mut control = QuadratureControl3::new(
        QuadratureOptions3 {
            depth: 0,
            ..Default::default()
        },
        &mut callback,
    )
    .unwrap();
    action(&mut control)
}
fn cartesian(limit: f64) -> CutElasticity3 {
    control(|control| {
        CutElasticity3::build_with_material_field(
            HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(),
            [1; 3],
            &Top,
            ElasticMaterialField3::Uniform(&material(limit)),
            &|p| p[0] == 0.0,
            Default::default(),
            control,
        )
        .unwrap()
    })
}
fn adaptive() -> AdaptiveElasticity3 {
    let tree = Octree3::uniform(1, 4, 4096).unwrap();
    let tree = tree
        .refined(&[*tree.leaves().iter().next().unwrap()], poll)
        .unwrap();
    control(|control| {
        AdaptiveElasticity3::build_with_material_field(
            HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(),
            &tree,
            &Top,
            ElasticMaterialField3::Uniform(&material(0.1)),
            &|p| p[0] == 0.0,
            Default::default(),
            control,
        )
        .unwrap()
    })
}
fn free_strain(_: usize, _: [f64; 3]) -> Result<[f64; 6], ElasticityError3> {
    Ok([EXPANSION, 0.0, 0.0, 0.0, 0.0, 0.0])
}
fn equilibrium(nodes: &[[f64; 3]]) -> Vec<f64> {
    nodes
        .iter()
        .flat_map(|p| [EXPANSION * p[0], 0.0, 0.0])
        .collect()
}
struct Identity;
impl FlexiblePreconditioner for Identity {
    fn apply(&self, _: usize, residual: &[f64], output: &mut [f64]) {
        output.copy_from_slice(residual);
    }
}
fn initial_config(n: usize) -> InitialSolveConfig {
    InitialSolveConfig {
        restart: n,
        max_cycles: 4,
        tolerance: 1e-12,
    }
}
fn acceleration(model: &CutSecondOrderProblem3<'_>, q: &[f64], v: &[f64], f: &[f64]) -> Vec<f64> {
    second_order_acceleration(
        model,
        q,
        v,
        f,
        &Identity,
        initial_config(model.dimension()),
        1_000_000,
        &mut || false,
    )
    .unwrap()
    .value
}
fn method(n: usize, dt: f64, rho_inf: f64) -> OperatorGeneralizedAlpha {
    OperatorGeneralizedAlpha::new(
        n,
        dt,
        rho_inf,
        ImplicitSolveConfig {
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1e-12,
                relative_tolerance: 1e-12,
                linear_restart: n,
                max_linear_cycles: 4,
                forcing_maximum: 1e-3,
                ..Default::default()
            },
            max_newton_iterations: 12,
        },
    )
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn close(a: f64, b: f64, relative: f64, absolute: f64) {
    assert!(
        (a - b).abs() <= absolute + relative * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}

#[test]
fn hot_equilibrium_uses_elastic_calibration_and_a_separate_total_strain_budget() {
    let op = cartesian(1e-6);
    let mass = op.prepare_material_mass(10_000, poll).unwrap();
    let occupancy = vec![1.0; op.cells()];
    let prepared = op.prepare_eigenstrain(&free_strain, POINTS, poll).unwrap();
    let gate = CancelGate::new_clock_free();
    let source = CutElasticityRef3::Cartesian(&op);
    // q=0 is outside the narrow elastic domain; construction must not assume it.
    let hot = CutSecondOrderProblem3::new_with_eigenstrain(
        source,
        &mass,
        &occupancy,
        &prepared,
        0.02,
        POINTS,
        op.n(),
        &gate,
    )
    .unwrap();
    let q = hot.restrict_force(&equilibrium(op.nodes())).unwrap();
    let zero = vec![0.0; hot.dimension()];
    assert!(hot.validate_state(&zero, &zero).is_err());
    hot.validate_state(&q, &zero).unwrap();
    let cold = CutSecondOrderProblem3::new(source, &mass, &occupancy, op.n(), &gate).unwrap();
    assert!(cold.validate_state(&q, &zero).is_err());
    let a = acceleration(&hot, &q, &zero, &zero);
    assert!(a.iter().all(|v| v.abs() < 1e-11));
    let mut state = SecondOrderState::new(0.0, &q, &zero, &a);
    hot.advance(
        &method(hot.dimension(), 0.02, 1.0),
        &mut state,
        &vec![0.0; op.n()],
    )
    .unwrap();
    for (actual, expected) in state.q.iter().zip(q) {
        close(*actual, expected, 1e-10, 1e-13);
    }
    let tight = CutSecondOrderProblem3::new_with_eigenstrain(
        source,
        &mass,
        &occupancy,
        &prepared,
        EXPANSION / 2.0,
        POINTS,
        op.n(),
        &gate,
    )
    .unwrap();
    assert!(tight.validate_state(&state.q, &zero).is_err());
}

#[test]
fn initial_thermal_acceleration_uses_scaled_free_force_and_consistent_physical_mass() {
    let mut op = cartesian(0.1);
    let prepared = op.prepare_eigenstrain(&free_strain, POINTS, poll).unwrap();
    op.set_scales(&[0.37]).unwrap();
    let mass = op.prepare_material_mass(10_000, poll).unwrap();
    let occupancy = vec![0.61; op.cells()];
    let gate = CancelGate::new_clock_free();
    let hot = CutSecondOrderProblem3::new_with_eigenstrain(
        CutElasticityRef3::Cartesian(&op),
        &mass,
        &occupancy,
        &prepared,
        0.02,
        POINTS,
        op.n(),
        &gate,
    )
    .unwrap();
    let n = hot.dimension();
    let zero = vec![0.0; n];
    let a = acceleration(&hot, &zero, &zero, &zero);
    let b = hot
        .restrict_force(&op.eigenstrain_load(&prepared, poll).unwrap())
        .unwrap();
    assert!(b.iter().any(|v| v.abs() > 1e-6));
    let mut dense_mass = vec![0.0; n * n];
    let mut full = vec![0.0; op.n()];
    let mut output = full.clone();
    for (column, &dof) in hot.free_dofs().iter().enumerate() {
        full.fill(0.0);
        full[dof] = 1.0;
        mass.apply_into(&occupancy, &full, &mut output, poll)
            .unwrap();
        for (row, &native) in hot.free_dofs().iter().enumerate() {
            dense_mass[row * n + column] = output[native];
        }
    }
    assert!((0..n).any(|i| (0..n).any(|j| i != j && dense_mass[i * n + j] != 0.0)));
    let mut reference = b.clone();
    fs_la::factor::lu(&dense_mass, n)
        .unwrap()
        .solve(&mut reference);
    for (actual, expected) in a.iter().zip(reference) {
        close(*actual, expected, 1e-10, 1e-12);
    }
    let mut resistance = zero.clone();
    hot.internal_force(&zero, &mut resistance);
    for (actual, force) in resistance.iter().zip(b) {
        close(*actual, -force, 1e-13, 1e-14);
    }
}

#[test]
fn thermal_oscillation_is_a_translated_cold_trajectory_on_cartesian_and_hanging_grids() {
    let cart = cartesian(0.1);
    let tree = adaptive();
    assert!(tree.physical_nodes().len() > tree.nodes().len());
    for source in [
        CutElasticityRef3::Cartesian(&cart),
        CutElasticityRef3::Adaptive(&tree),
    ] {
        let (mass, prepared, nodes, native_n) = match source {
            CutElasticityRef3::Cartesian(op) => (
                op.prepare_material_mass(POINTS, poll).unwrap(),
                op.prepare_eigenstrain(&free_strain, POINTS, poll).unwrap(),
                op.nodes(),
                op.n(),
            ),
            CutElasticityRef3::Adaptive(op) => (
                op.prepare_material_mass(POINTS, poll).unwrap(),
                op.prepare_eigenstrain(&free_strain, POINTS, poll).unwrap(),
                op.nodes(),
                op.n(),
            ),
        };
        let occupancy = vec![0.8; mass.cells()];
        let gate = CancelGate::new_clock_free();
        let hot = CutSecondOrderProblem3::new_with_eigenstrain(
            source, &mass, &occupancy, &prepared, 0.02, POINTS, native_n, &gate,
        )
        .unwrap();
        let cold = CutSecondOrderProblem3::new(source, &mass, &occupancy, native_n, &gate).unwrap();
        let n = hot.dimension();
        let offset = hot.restrict_force(&equilibrium(nodes)).unwrap();
        let x: Vec<_> = (0..n).map(|i| 1e-5 * ((i % 5) as f64 - 2.0)).collect();
        let v: Vec<_> = (0..n).map(|i| 2e-5 * ((i % 7) as f64 - 3.0)).collect();
        let q: Vec<_> = offset.iter().zip(&x).map(|(q, x)| q + x).collect();
        let zero = vec![0.0; n];
        let cold_a = acceleration(&cold, &x, &v, &zero);
        let hot_a = acceleration(&hot, &q, &v, &zero);
        for (a, b) in cold_a.iter().zip(&hot_a) {
            close(*a, *b, 1e-8, 1e-11);
        }
        for rho_inf in [1.0, 0.6] {
            let mut hot_state = SecondOrderState::new(0.0, &q, &v, &hot_a);
            let mut cold_state = SecondOrderState::new(0.0, &x, &v, &cold_a);
            let stepper = method(n, 0.02, rho_inf);
            for _ in 0..3 {
                hot.advance(&stepper, &mut hot_state, &vec![0.0; native_n])
                    .unwrap();
                cold.advance(&stepper, &mut cold_state, &vec![0.0; native_n])
                    .unwrap();
                for ((actual, shift), expected) in
                    hot_state.q.iter().zip(&offset).zip(&cold_state.q)
                {
                    close(actual - shift, *expected, 2e-7, 2e-11);
                }
                for (actual, expected) in hot_state
                    .v
                    .iter()
                    .zip(&cold_state.v)
                    .chain(hot_state.a.iter().zip(&cold_state.a))
                {
                    close(*actual, *expected, 2e-7, 2e-10);
                }
                assert_eq!(hot_state.t, cold_state.t);
                assert_eq!(hot_state.steps, cold_state.steps);
            }
        }
        if let CutElasticityRef3::Adaptive(op) = source {
            let physical = op
                .physical_displacements(&hot.expand_field(&offset).unwrap())
                .unwrap();
            for (p, u) in op.physical_nodes().iter().zip(physical.chunks_exact(3)) {
                close(u[0], EXPANSION * p[0], 1e-12, 1e-14);
            }
        }
    }
}

#[test]
fn affine_tangent_and_existing_initial_and_step_adjoints_have_no_free_force_offset() {
    let op = cartesian(1e-5);
    let mass = op.prepare_material_mass(10_000, poll).unwrap();
    let occupancy = vec![1.0; op.cells()];
    let prepared = op.prepare_eigenstrain(&free_strain, POINTS, poll).unwrap();
    let gate = CancelGate::new_clock_free();
    let hot = CutSecondOrderProblem3::new_with_eigenstrain(
        CutElasticityRef3::Cartesian(&op),
        &mass,
        &occupancy,
        &prepared,
        0.02,
        POINTS,
        op.n(),
        &gate,
    )
    .unwrap();
    let n = hot.dimension();
    let q = hot.restrict_force(&equilibrium(op.nodes())).unwrap();
    let zero = vec![0.0; n];
    let direction: Vec<_> = (0..n).map(|i| ((i % 5) as f64 - 2.0) / 3.0).collect();
    let seed: Vec<_> = (0..n).map(|i| 0.2 + 0.03 * i as f64).collect();
    let mut tangent = zero.clone();
    let mut transpose = zero.clone();
    hot.tangent_apply(&q, &direction, &mut tangent);
    hot.tangent_transpose_apply(&q, &seed, &mut transpose)
        .unwrap();
    close(
        dot(&seed, &tangent),
        dot(&direction, &transpose),
        1e-12,
        1e-12,
    );
    let epsilon = 1e-8;
    let plus: Vec<_> = q
        .iter()
        .zip(&direction)
        .map(|(q, d)| q + epsilon * d)
        .collect();
    let minus: Vec<_> = q
        .iter()
        .zip(&direction)
        .map(|(q, d)| q - epsilon * d)
        .collect();
    let mut fp = zero.clone();
    let mut fm = zero.clone();
    hot.internal_force(&plus, &mut fp);
    hot.internal_force(&minus, &mut fm);
    for ((p, m), k) in fp.iter().zip(&fm).zip(&tangent) {
        close((p - m) / (2.0 * epsilon), *k, 1e-8, 2e-9);
    }
    hot.tangent_apply(&q, &zero, &mut tangent);
    assert!(tangent.iter().all(|v| *v == 0.0));

    let initial_gradient = second_order_acceleration_vjp(
        &hot,
        &q,
        &zero,
        &zero,
        &seed,
        &Identity,
        &Identity,
        initial_config(n),
        1_000_000,
        &mut || false,
    )
    .unwrap();
    let numerical = (dot(&seed, &acceleration(&hot, &plus, &zero, &zero))
        - dot(&seed, &acceleration(&hot, &minus, &zero, &zero)))
        / (2.0 * epsilon);
    close(
        dot(&initial_gradient.initial_q, &direction),
        numerical,
        2e-8,
        2e-7,
    );
    assert!(initial_gradient.parameters.is_empty());

    let initial = SecondOrderState::new(0.0, &q, &zero, &initial_gradient.acceleration);
    let stepper = method(n, 0.03, 0.6);
    let terminal_v: Vec<_> = seed.iter().map(|s| 0.07 * s).collect();
    let terminal_a: Vec<_> = seed.iter().map(|s| -0.003 * s).collect();
    let gradient = stepper
        .step_vjp(
            &initial,
            &hot,
            &zero,
            (&seed, &terminal_v, &terminal_a),
            &Identity,
            SecondOrderAdjointConfig {
                restart: n,
                max_cycles: 4,
                tolerance: 1e-12,
            },
            1_000_000,
            &mut || false,
        )
        .unwrap();
    let dq: Vec<_> = direction.iter().map(|d| 1e-4 * d).collect();
    let dv: Vec<_> = direction.iter().map(|d| -2e-3 * d).collect();
    let da: Vec<_> = direction.iter().map(|d| 0.01 * d).collect();
    let df: Vec<_> = direction.iter().map(|d| 0.03 * d).collect();
    let h = 1e-4;
    let sample = |sign: f64| {
        let mut state = initial.clone();
        for (x, d) in state
            .q
            .iter_mut()
            .zip(&dq)
            .chain(state.v.iter_mut().zip(&dv))
            .chain(state.a.iter_mut().zip(&da))
        {
            *x += sign * h * d;
        }
        let force: Vec<_> = df.iter().map(|d| sign * h * d).collect();
        stepper.step(&mut state, &hot, &force).unwrap();
        dot(&seed, &state.q) + dot(&terminal_v, &state.v) + dot(&terminal_a, &state.a)
    };
    let analytic = dot(&gradient.initial_q, &dq)
        + dot(&gradient.initial_v, &dv)
        + dot(&gradient.initial_a, &da)
        + dot(&gradient.forcing, &df);
    close(
        analytic,
        (sample(1.0) - sample(-1.0)) / (2.0 * h),
        2e-7,
        2e-9,
    );
    assert!(gradient.parameters.is_empty());
}

#[test]
fn thermal_owner_budgets_endpoint_refusal_and_late_cancellation_preserve_state() {
    let op = cartesian(1e-6);
    let foreign = cartesian(1e-6);
    let mass = op.prepare_material_mass(10_000, poll).unwrap();
    let occupancy = vec![1.0; op.cells()];
    let prepared = op.prepare_eigenstrain(&free_strain, POINTS, poll).unwrap();
    let other = foreign
        .prepare_eigenstrain(&free_strain, POINTS, poll)
        .unwrap();
    let gate = CancelGate::new_clock_free();
    let source = CutElasticityRef3::Cartesian(&op);
    assert!(
        CutSecondOrderProblem3::new_with_eigenstrain(
            source,
            &mass,
            &occupancy,
            &other,
            0.02,
            POINTS,
            op.n(),
            &gate,
        )
        .is_err()
    );
    assert!(
        CutSecondOrderProblem3::new_with_eigenstrain(
            source,
            &mass,
            &occupancy,
            &prepared,
            0.02,
            prepared.bulk_points() - 1,
            op.n(),
            &gate,
        )
        .is_err()
    );
    for budget in [0.0, 1.0, f64::NAN] {
        assert!(
            CutSecondOrderProblem3::new_with_eigenstrain(
                source,
                &mass,
                &occupancy,
                &prepared,
                budget,
                POINTS,
                op.n(),
                &gate,
            )
            .is_err()
        );
    }
    let hot = CutSecondOrderProblem3::new_with_eigenstrain(
        source,
        &mass,
        &occupancy,
        &prepared,
        0.02,
        POINTS,
        op.n(),
        &gate,
    )
    .unwrap();
    let n = hot.dimension();
    let q = hot.restrict_force(&equilibrium(op.nodes())).unwrap();
    let zero = vec![0.0; n];
    let dt = 1e-4;
    // Midpoint strain stays below the 1e-6 elastic ceiling; the candidate
    // endpoint reaches approximately 1.5e-6 and must not be published.
    let velocity: Vec<_> = q.iter().map(|q| q / EXPANSION * (1.5e-6 / dt)).collect();
    let a = acceleration(&hot, &q, &velocity, &zero);
    let mut state = SecondOrderState::new(0.0, &q, &velocity, &a);
    let before = state.clone();
    assert!(matches!(
        hot.advance(&method(n, dt, 1.0), &mut state, &vec![0.0; op.n()]),
        Err(TimeSolveError::ModelRefused {
            phase: "accepted endpoint",
            ..
        })
    ));
    assert_eq!(state, before);

    let mut state = SecondOrderState::new(0.0, &q, &zero, &acceleration(&hot, &q, &zero, &zero));
    hot.advance(&method(n, dt, 1.0), &mut state, &vec![0.0; op.n()])
        .unwrap();
    let accepted = state.clone();
    let mut calls = 0;
    let result = method(n, dt, 1.0).step_controlled(&mut state, &hot, &zero, &mut || {
        calls += 1;
        if calls == 6 {
            gate.request();
        }
        gate.is_requested()
    });
    assert!(matches!(result, Err(TimeSolveError::Cancelled)));
    assert!(calls >= 6);
    assert_eq!(state, accepted);
    assert!(matches!(
        CutSecondOrderProblem3::new_with_eigenstrain(
            source,
            &mass,
            &occupancy,
            &prepared,
            0.02,
            POINTS,
            op.n(),
            &gate,
        ),
        Err(ElasticityError3::Cancelled)
    ));
}
