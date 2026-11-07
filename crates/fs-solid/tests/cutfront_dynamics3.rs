//! G1/G3: original CutFEM material inertia, native motion, work and pullbacks.
use std::ops::ControlFlow;

use fs_blake3::ContentHash;
use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
use fs_cutfem::elastic3::{CutElasticity3, ElasticMaterialField3};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_exec::CancelGate;
use fs_ivl::Interval;
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_solid::cutfront::dynamics3::{CutElasticityRef3, CutSecondOrderProblem3};
use fs_solver::{FlexiblePreconditioner, LinearOp, NewtonKrylovConfig};
use fs_time::galpha::initialization::{InitialSolveConfig, second_order_acceleration};
use fs_time::galpha::second_order_adjoint::SecondOrderAdjointConfig;
use fs_time::galpha::{
    GeneralizedAlpha, ImplicitSolveConfig, OperatorGeneralizedAlpha, SecondOrderProblem,
    SecondOrderState, TimeSolveError, galpha_step,
};

struct Top;
impl CutSdf3 for Top {
    fn value(&self, p: [f64; 3]) -> f64 {
        p[2] - 0.73
    }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        Interval::new(lo[2], hi[2]) - Interval::new(0.73, 0.73)
    }
    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], axis: HeightAxis) -> Interval {
        let value = if axis == HeightAxis::Z { 1.0 } else { 0.0 };
        Interval::new(value, value)
    }
}
fn poll() -> ControlFlow<()> {
    ControlFlow::Continue(())
}
fn material() -> ElasticAcousticState {
    let mut c = std::array::from_fn(|i| {
        std::array::from_fn(|j| if i == j { 15.0 + i as f64 } else { 0.0 })
    });
    c[0][3] = 1.2;
    c[3][0] = 1.2;
    ElasticAcousticState::try_new_mandel(c, 3.0, ContentHash([42; 32]), 0.1).unwrap()
}
fn cartesian(clamped: bool) -> CutElasticity3 {
    let mut callback = |_| poll();
    let mut control = QuadratureControl3::new(
        QuadratureOptions3 {
            depth: 0,
            ..Default::default()
        },
        &mut callback,
    )
    .unwrap();
    CutElasticity3::build_with_material_field(
        HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(),
        [1; 3],
        &Top,
        ElasticMaterialField3::Uniform(&material()),
        &|p| clamped && p[0] == 0.0,
        Default::default(),
        &mut control,
    )
    .unwrap()
}
fn adaptive() -> AdaptiveElasticity3 {
    let tree = Octree3::uniform(1, 4, 4096).unwrap();
    let tree = tree
        .refined(&[*tree.leaves().iter().next().unwrap()], poll)
        .unwrap();
    let mut callback = |_| poll();
    let mut control = QuadratureControl3::new(
        QuadratureOptions3 {
            depth: 0,
            ..Default::default()
        },
        &mut callback,
    )
    .unwrap();
    AdaptiveElasticity3::build_with_material_field(
        HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(),
        &tree,
        &Top,
        ElasticMaterialField3::Uniform(&material()),
        &|p| p[0] == 0.0,
        Default::default(),
        &mut control,
    )
    .unwrap()
}
struct Identity;
impl FlexiblePreconditioner for Identity {
    fn apply(&self, _: usize, residual: &[f64], output: &mut [f64]) {
        output.copy_from_slice(residual);
    }
}
fn method(n: usize, dt: f64) -> OperatorGeneralizedAlpha {
    OperatorGeneralizedAlpha::new(
        n,
        dt,
        1.0,
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
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1.0),
        "{a} != {b}"
    );
}
fn energy(model: &CutSecondOrderProblem3<'_>, state: &SecondOrderState) -> f64 {
    let mut work = vec![0.0; model.dimension()];
    model.mass_apply(&state.v, &mut work);
    let kinetic = 0.5 * dot(&state.v, &work);
    model.internal_force(&state.q, &mut work);
    kinetic + 0.5 * dot(&state.q, &work)
}

#[test]
fn native_material_body_force_drives_free_cut_body_with_exact_acceleration_and_work() {
    let op = cartesian(false);
    let mass = op.prepare_material_mass(10_000, poll).unwrap();
    let occupancy = vec![1.0; op.cells()];
    let gate = CancelGate::new_clock_free();
    let model = CutSecondOrderProblem3::new(
        CutElasticityRef3::Cartesian(&op),
        &mass,
        &occupancy,
        op.n(),
        &gate,
    )
    .unwrap();
    let acceleration = [0.7, -0.3, 0.2];
    let body = op
        .prepare_material_body_load(&|_| acceleration, 10_000, poll)
        .unwrap();
    let native_force = body.load(&occupancy, poll).unwrap();
    let force = model.restrict_force(&native_force).unwrap();
    let n = model.dimension();
    let q0 = vec![0.0; n];
    let v0: Vec<_> = model
        .free_dofs()
        .iter()
        .map(|dof| [0.02, -0.03, 0.01][dof % 3])
        .collect();
    let initial = second_order_acceleration(
        &model,
        &q0,
        &v0,
        &force,
        &Identity,
        InitialSolveConfig {
            restart: n,
            max_cycles: 4,
            tolerance: 1e-12,
        },
        100_000,
        &mut || gate.is_requested(),
    )
    .unwrap();
    for (&dof, &value) in model.free_dofs().iter().zip(&initial.value) {
        close(value, acceleration[dof % 3], 1e-10);
    }
    let mut state = SecondOrderState::new(0.0, &q0, &v0, &initial.value);
    let e0 = energy(&model, &state);
    let stepper = method(n, 0.005);
    for _ in 0..6 {
        model.advance(&stepper, &mut state, &native_force).unwrap();
        for (i, &dof) in model.free_dofs().iter().enumerate() {
            close(
                state.q[i],
                v0[i] * state.t + 0.5 * acceleration[dof % 3] * state.t.powi(2),
                1e-10,
            );
            close(state.v[i], v0[i] + acceleration[dof % 3] * state.t, 1e-10);
        }
    }
    close(energy(&model, &state) - e0, dot(&force, &state.q), 1e-11);
    let speed2: f64 = [0.02, -0.03, 0.01]
        .iter()
        .zip(acceleration)
        .map(|(v, a)| (v + a * state.t).powi(2))
        .sum();
    close(energy(&model, &state), 0.5 * 3.0 * 0.73 * speed2, 1e-11);
    assert_eq!(state.steps, 6);
}

#[test]
fn clamped_cartesian_and_hanging_node_motion_match_original_consistent_dense_pencils() {
    let cart = cartesian(true);
    let adaptive = adaptive();
    assert!(adaptive.physical_nodes().len() > adaptive.nodes().len());
    for source in [
        CutElasticityRef3::Cartesian(&cart),
        CutElasticityRef3::Adaptive(&adaptive),
    ] {
        let (operator, mass, native_force): (&dyn LinearOp, _, _) = match source {
            CutElasticityRef3::Cartesian(op) => (
                op,
                op.prepare_material_mass(100_000, poll).unwrap(),
                op.body_load(&|_| [0.001, -0.002, 0.003], poll).unwrap(),
            ),
            CutElasticityRef3::Adaptive(op) => (
                op,
                op.prepare_material_mass(100_000, poll).unwrap(),
                op.body_load(&|_| [0.001, -0.002, 0.003], poll).unwrap(),
            ),
        };
        let occupancy = vec![1.0; mass.cells()];
        let gate = CancelGate::new_clock_free();
        let model =
            CutSecondOrderProblem3::new(source, &mass, &occupancy, operator.n(), &gate).unwrap();
        let n = model.dimension();
        assert!(n < operator.n());
        let mut m = vec![0.0; n * n];
        let mut k = vec![0.0; n * n];
        let mut full = vec![0.0; operator.n()];
        let mut output = full.clone();
        for (column, &dof) in model.free_dofs().iter().enumerate() {
            full.fill(0.0);
            full[dof] = 1.0;
            mass.apply_into(&occupancy, &full, &mut output, poll)
                .unwrap();
            for (row, &i) in model.free_dofs().iter().enumerate() {
                m[row * n + column] = output[i];
            }
            operator.apply(&full, &mut output);
            for (row, &i) in model.free_dofs().iter().enumerate() {
                k[row * n + column] = output[i];
            }
        }
        assert!((0..n).any(|i| (0..n).any(|j| i != j && m[i * n + j] != 0.0)));
        let force = model.restrict_force(&native_force).unwrap();
        let mut q = vec![0.0; n];
        let mut v = q.clone();
        let mut a = force.clone();
        fs_la::factor::lu(&m, n).unwrap().solve(&mut a);
        let mut state = SecondOrderState::new(0.0, &q, &v, &a);
        let dense = GeneralizedAlpha::new(&m, &vec![0.0; n * n], &k, n, 0.03, 1.0);
        let stepper = method(n, 0.03);
        for _ in 0..3 {
            model.advance(&stepper, &mut state, &native_force).unwrap();
            galpha_step(&dense, &mut q, &mut v, &mut a, &force);
            for (actual, reference) in state
                .q
                .iter()
                .zip(&q)
                .chain(state.v.iter().zip(&v))
                .chain(state.a.iter().zip(&a))
            {
                close(*actual, *reference, 2e-9);
            }
        }
        close(energy(&model, &state), dot(&force, &state.q), 2e-10);
        let expanded = model.expand_field(&state.q).unwrap();
        let fixed = match source {
            CutElasticityRef3::Cartesian(op) => op.fixed(),
            CutElasticityRef3::Adaptive(op) => op.fixed(),
        };
        for (node, fixed) in fixed.iter().enumerate() {
            if *fixed {
                assert_eq!(&expanded[3 * node..3 * node + 3], &[0.0; 3]);
            }
        }
    }
}

#[test]
fn fixed_configuration_force_pullback_matches_a_real_cut_body_step() {
    let op = cartesian(true);
    let mass = op.prepare_material_mass(10_000, poll).unwrap();
    let occupancy = vec![1.0; op.cells()];
    let gate = CancelGate::new_clock_free();
    let model = CutSecondOrderProblem3::new(
        CutElasticityRef3::Cartesian(&op),
        &mass,
        &occupancy,
        op.n(),
        &gate,
    )
    .unwrap();
    let n = model.dimension();
    let zero = vec![0.0; n];
    let initial = SecondOrderState::new(0.0, &zero, &zero, &zero);
    let force = vec![0.001; n];
    let seed: Vec<_> = (0..n).map(|i| 0.2 + 0.03 * i as f64).collect();
    let direction: Vec<_> = (0..n).map(|i| ((i % 5) as f64 - 2.0) / 3.0).collect();
    let stepper = method(n, 0.03);
    let gradient = stepper
        .step_vjp(
            &initial,
            &model,
            &force,
            (&seed, &zero, &zero),
            &Identity,
            SecondOrderAdjointConfig {
                restart: n,
                max_cycles: 4,
                tolerance: 1e-12,
            },
            100_000,
            &mut || gate.is_requested(),
        )
        .unwrap();
    assert!(gradient.parameters.is_empty());
    let sample = |sign: f64| {
        let input: Vec<_> = force
            .iter()
            .zip(&direction)
            .map(|(f, d)| f + sign * 1e-4 * d)
            .collect();
        let mut state = initial.clone();
        stepper.step(&mut state, &model, &input).unwrap();
        dot(&seed, &state.q)
    };
    close(
        dot(&gradient.forcing, &direction),
        (sample(1.0) - sample(-1.0)) / 2e-4,
        1e-8,
    );
}

#[test]
fn foreign_mass_invalid_force_strain_and_cancellation_preserve_the_native_trajectory() {
    let op = cartesian(true);
    let foreign = cartesian(true);
    let mass = op.prepare_material_mass(10_000, poll).unwrap();
    let other_mass = foreign.prepare_material_mass(10_000, poll).unwrap();
    let occupancy = vec![1.0; op.cells()];
    let gate = CancelGate::new_clock_free();
    let source = CutElasticityRef3::Cartesian(&op);
    assert!(CutSecondOrderProblem3::new(source, &other_mass, &occupancy, op.n(), &gate).is_err());
    assert!(CutSecondOrderProblem3::new(source, &mass, &occupancy, op.n() - 1, &gate).is_err());
    assert!(
        CutSecondOrderProblem3::new(source, &mass, &vec![0.0; op.cells()], op.n(), &gate).is_err()
    );
    let model = CutSecondOrderProblem3::new(source, &mass, &occupancy, op.n(), &gate).unwrap();
    let n = model.dimension();
    let zero = vec![0.0; n];
    let mut state = SecondOrderState::new(0.0, &zero, &zero, &zero);
    let force = vec![0.001; op.n()];
    let stepper = method(n, 0.01);
    model.advance(&stepper, &mut state, &force).unwrap();
    let accepted = state.clone();
    let mut invalid = force.clone();
    invalid[0] = f64::NAN;
    assert!(model.advance(&stepper, &mut state, &invalid).is_err());
    assert_eq!(state, accepted);
    assert!(model.validate_state(&vec![1.0; n], &zero).is_err());
    gate.request();
    assert!(matches!(
        model.advance(&stepper, &mut state, &force),
        Err(TimeSolveError::Cancelled)
    ));
    assert_eq!(state, accepted);
}
