//! G1/G3/G4: full physical inertia survives initialization, time stepping and AD.
//! The small manufactured displacement/slope pencil has nonzero mass cross
//! terms; independent dense LU solves provide the forward and gradient oracle.
use fs_la::factor::lu;
use fs_plate::PlateModel;
use fs_plate_transient::{
    NoPlateLoad, PlateDynamics, PlateDynamicsBudget, PlateDynamicsError, PlateDynamicsParameters,
    PlateLoad,
};
use fs_solver::{FlexiblePreconditioner, NewtonKrylovConfig};
use fs_sparse::{Coo, Csr, LdltOptions};
use fs_time::galpha::{
    GeneralizedAlpha, ImplicitSolveConfig, OperatorGeneralizedAlpha, SecondOrderProblem,
    SecondOrderState, galpha_step,
    initialization::{
        InitialSolveConfig, initial_workspace_components, second_order_acceleration,
        second_order_acceleration_vjp,
    },
    second_order_adjoint::{
        SecondOrderAdjointConfig, SecondOrderVjp,
        trajectory::{
            RecordedStructural, StructuralRecordingConfig, StructuralReplayBudget,
            StructuralTrajectoryModel, samples::StructuralSampleObjective,
        },
    },
};

const N: usize = 3;
const STEP: f64 = 0.01;
const STEPS: usize = 24;
const RHO: f64 = 0.8;
const SAMPLES: [usize; 4] = [0, 8, 16, STEPS];
const Q0: [f64; N] = [0.001, -0.003, 0.002];
const V0: [f64; N] = [0.01, -0.02, 0.03];
const SENSOR: [f64; N] = [1.0, 0.2, -0.3];

fn matrix(values: [f64; 9]) -> Csr {
    let mut coo = Coo::new(N, N);
    for row in 0..N {
        for column in 0..N {
            if values[N * row + column] != 0.0 {
                coo.push(row, column, values[N * row + column]);
            }
        }
    }
    coo.assemble()
}
fn plate() -> PlateModel {
    PlateModel {
        k: matrix([40.0, -6.0, 3.0, -6.0, 14.0, -2.0, 3.0, -2.0, 26.0]),
        m: matrix([2.0, 0.12, -0.09, 0.12, 0.3, 0.025, -0.09, 0.025, 0.4]),
        dof_map: vec![Some(0), Some(1), Some(2)],
        free: N,
    }
}
fn parameters(p: [f64; 5]) -> PlateDynamicsParameters {
    PlateDynamicsParameters {
        stiffness_scale: p[0],
        mass_scale: p[1],
        mass_damping_per_s: p[2],
        stiffness_damping_s: p[3],
    }
}
fn budget(plate: &PlateModel) -> PlateDynamicsBudget {
    PlateDynamicsBudget {
        max_dofs: N,
        max_nonzeros: plate.k.nnz() + plate.m.nnz(),
        max_full_dofs: N,
        max_load_parameters: 1,
    }
}
fn initial_config() -> InitialSolveConfig {
    InitialSolveConfig {
        restart: N,
        max_cycles: 4,
        tolerance: 1e-13,
    }
}
fn method(rho: f64) -> OperatorGeneralizedAlpha {
    OperatorGeneralizedAlpha::new(
        N,
        STEP,
        rho,
        ImplicitSolveConfig {
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1e-13,
                relative_tolerance: 1e-12,
                linear_restart: N,
                max_linear_cycles: 4,
                forcing_minimum: 1e-14,
                forcing_maximum: 1e-12,
                ..NewtonKrylovConfig::default()
            },
            max_newton_iterations: 8,
        },
    )
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn close(actual: &[f64], expected: &[f64]) {
    for (&a, &b) in actual.iter().zip(expected) {
        assert!((a - b).abs() < 2e-10 * b.abs().max(1e-3), "{a:e} != {b:e}");
    }
}
struct Load(f64);
impl PlateLoad for Load {
    fn parameter_count(&self) -> usize {
        1
    }
    fn forcing(&self, t: f64, out: &mut [f64]) -> Result<(), String> {
        for (value, direction) in
            out.iter_mut()
                .zip([1.0 + t, -0.3 * (1.0 + 0.5 * t), 0.2 * (1.0 - 0.2 * t)])
        {
            *value = self.0 * direction;
        }
        Ok(())
    }
    fn forcing_vjp(&self, t: f64, seed: &[f64], out: &mut [f64]) -> Result<(), String> {
        out[0] = dot(
            seed,
            &[1.0 + t, -0.3 * (1.0 + 0.5 * t), 0.2 * (1.0 - 0.2 * t)],
        );
        Ok(())
    }
}
struct Observation;
impl StructuralSampleObjective for Observation {
    fn evaluate(
        &self,
        _: usize,
        state: &SecondOrderState,
        bars: (&mut [f64], &mut [f64], &mut [f64]),
        parameters: &mut [f64],
    ) -> Result<f64, String> {
        let value = dot(&state.q, &SENSOR);
        for (bar, weight) in bars.0.iter_mut().zip(SENSOR) {
            *bar = value * weight;
        }
        bars.1.fill(0.0);
        bars.2.fill(0.0);
        parameters.fill(0.0);
        Ok(0.5 * value * value)
    }
}
fn initial<L: PlateLoad>(model: &PlateDynamics<'_, L>) -> SecondOrderState {
    let mut force = vec![0.0; N];
    model.forcing(0.0, &mut force).unwrap();
    let cfg = initial_config();
    let a = second_order_acceleration(
        model,
        &Q0,
        &V0,
        &force,
        &model.mass_preconditioner(),
        cfg,
        initial_workspace_components(N, model.parameter_count(), cfg).unwrap(),
        &mut || false,
    )
    .unwrap();
    SecondOrderState::new(0.0, &Q0, &V0, &a.value)
}
fn dense_run(plate: &PlateModel, p: [f64; 5]) -> (SecondOrderState, f64) {
    let mass: Vec<_> = plate.m.to_dense().iter().map(|v| p[1] * v).collect();
    let stiffness: Vec<_> = plate.k.to_dense().iter().map(|v| p[0] * v).collect();
    let damping: Vec<_> = mass
        .iter()
        .zip(&stiffness)
        .map(|(m, k)| p[2] * m + p[3] * k)
        .collect();
    let load = Load(p[4]);
    let mut force = vec![0.0; N];
    load.forcing(0.0, &mut force).unwrap();
    let mut acceleration: Vec<_> = (0..N)
        .map(|i| {
            force[i]
                - dot(&stiffness[N * i..N * (i + 1)], &Q0)
                - dot(&damping[N * i..N * (i + 1)], &V0)
        })
        .collect();
    lu(&mass, N).unwrap().solve(&mut acceleration);
    let mut state = SecondOrderState::new(0.0, &Q0, &V0, &acceleration);
    let method = GeneralizedAlpha::new(&mass, &damping, &stiffness, N, STEP, RHO);
    let mut values = vec![0.5 * dot(&state.q, &SENSOR).powi(2)];
    for step in 1..=STEPS {
        load.forcing(state.t + STEP / (1.0 + RHO), &mut force)
            .unwrap();
        galpha_step(&method, &mut state.q, &mut state.v, &mut state.a, &force);
        state.t += STEP;
        state.steps += 1;
        if SAMPLES.contains(&step) {
            values.push(0.5 * dot(&state.q, &SENSOR).powi(2));
        }
    }
    (state, values.into_iter().rev().sum())
}

#[test]
fn consistent_inertia_initialization_trajectory_and_all_parameter_adjoints_match_dense() {
    let plate = plate();
    let p = [1.1, 0.9, 0.4, 0.03, 1.7];
    let load = Load(p[4]);
    let dynamics =
        PlateDynamics::new(&plate, parameters(p), &load, budget(&plate), &mut || false).unwrap();
    let start = initial(&dynamics);
    let cfg = SecondOrderAdjointConfig {
        restart: N,
        max_cycles: 4,
        tolerance: 1e-12,
    };
    let method = method(RHO);
    let mut tape = RecordedStructural::new(
        method,
        &dynamics,
        &start,
        StructuralRecordingConfig {
            steps: STEPS,
            adjoint: cfg,
            max_workspace_components: method.adjoint_workspace_components(5, cfg).unwrap(),
        },
    )
    .unwrap();
    tape.advance(STEPS, STEPS, &mut || false).unwrap();
    let result = tape
        .pullback_samples(
            &SAMPLES,
            SAMPLES.len(),
            &Observation,
            &dynamics
                .effective_preconditioner(STEP, RHO, &mut || false)
                .unwrap(),
            StructuralReplayBudget {
                checkpoints: tape.required_checkpoints(),
                forward_steps: 256,
            },
            &mut || false,
        )
        .unwrap();
    let mass = dynamics.mass_preconditioner();
    let cfg = initial_config();
    let mut force = vec![0.0; N];
    dynamics.forcing(0.0, &mut force).unwrap();
    let initial_gradient = second_order_acceleration_vjp(
        &dynamics,
        &start.q,
        &start.v,
        &force,
        &result.gradient.initial_a,
        &mass,
        &mass,
        cfg,
        initial_workspace_components(N, 5, cfg).unwrap(),
        &mut || false,
    )
    .unwrap();
    let mut force_gradient = vec![0.0; 5];
    dynamics
        .forcing_vjp(0.0, &initial_gradient.forcing, &mut force_gradient)
        .unwrap();
    let gradient: Vec<_> = (0..5)
        .map(|i| result.gradient.parameters[i] + initial_gradient.parameters[i] + force_gradient[i])
        .collect();
    let (expected, value) = dense_run(&plate, p);
    close(&tape.state().q, &expected.q);
    close(&tape.state().v, &expected.v);
    close(&tape.state().a, &expected.a);
    assert!((result.value - value).abs() < 2e-11 * value);
    for i in 0..5 {
        let eps = 2e-6;
        let (mut plus, mut minus) = (p, p);
        plus[i] += eps;
        minus[i] -= eps;
        let fd = (dense_run(&plate, plus).1 - dense_run(&plate, minus).1) / (2.0 * eps);
        assert!(fd.abs() > 1e-7, "uninformative parameter {i}");
        assert!(
            (gradient[i] - fd).abs() < 2e-6 * fd.abs(),
            "parameter {i}: {} vs {fd}",
            gradient[i]
        );
    }
    // The diagonal is a preconditioner, never a replacement inertia operator.
    let mut jacobi = vec![0.0; N];
    let mut reconstructed = jacobi.clone();
    mass.apply(0, &force, &mut jacobi);
    dynamics.mass_apply(&jacobi, &mut reconstructed);
    assert!(
        force
            .iter()
            .zip(&reconstructed)
            .any(|(f, r)| (f - r).abs() > 0.01)
    );
}

#[test]
fn full_inertia_midpoint_vibration_conserves_or_dissipates_physical_energy() {
    let plate = plate();
    for damping in [[0.0, 0.0], [0.4, 0.03]] {
        let p = [1.1, 0.9, damping[0], damping[1], 0.0];
        let dynamics = PlateDynamics::new(
            &plate,
            parameters(p),
            &NoPlateLoad,
            budget(&plate),
            &mut || false,
        )
        .unwrap();
        let mut state = initial(&dynamics);
        let method = method(1.0);
        let energy = |state: &SecondOrderState| {
            let mut inertia = vec![0.0; N];
            let mut stiffness = inertia.clone();
            plate.m.spmv(&state.v, &mut inertia);
            plate.k.spmv(&state.q, &mut stiffness);
            0.5 * (p[1] * dot(&state.v, &inertia) + p[0] * dot(&state.q, &stiffness))
        };
        let start = energy(&state);
        let mut previous = start;
        for _ in 0..240 {
            method.step(&mut state, &dynamics, &[0.0; N]).unwrap();
            let current = energy(&state);
            assert!(current >= 0.0 && current <= previous + 2e-10 * start);
            if damping[0] == 0.0 {
                assert!((current - start).abs() < 2e-9 * start);
            }
            previous = current;
        }
        if damping[0] > 0.0 {
            assert!(previous < 0.2 * start);
        }
    }
}

#[test]
fn nonphysical_consistent_mass_and_factor_budget_are_refused_before_publication() {
    let p = parameters([1.0, 1.0, 0.0, 0.0, 0.0]);
    for values in [
        [1.0, 2.0, 0.0, 2.0, 1.0, 0.0, 0.0, 0.0, 1.0], // positive diagonal, indefinite
        [1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0], // singular
        [1.0, 0.4, 0.0, 0.1, 1.0, 0.0, 0.0, 0.0, 1.0], // nonsymmetric
        [1.0, f64::NAN, 0.0, f64::NAN, 1.0, 0.0, 0.0, 0.0, 1.0],
    ] {
        let mut model = plate();
        model.m = matrix(values);
        assert!(
            PlateDynamics::new(&model, p, &NoPlateLoad, budget(&model), &mut || false).is_err()
        );
    }
    let model = plate();
    let limited = LdltOptions {
        max_factor_bytes: Some(0),
        ..LdltOptions::default()
    };
    assert!(matches!(
        PlateDynamics::new_with_mass_options(
            &model,
            p,
            &NoPlateLoad,
            budget(&model),
            &limited,
            &mut || false
        ),
        Err(PlateDynamicsError::MassFactor(_))
    ));
    let mut diagonal = model.clone();
    diagonal.m = matrix([2.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0, 0.4]);
    PlateDynamics::new_with_mass_options(
        &diagonal,
        p,
        &NoPlateLoad,
        budget(&diagonal),
        &limited,
        &mut || false,
    )
    .unwrap(); // no unnecessary factor for lumped mass
    let mut polls = 0;
    PlateDynamics::new(&model, p, &NoPlateLoad, budget(&model), &mut || {
        polls += 1;
        false
    })
    .unwrap();
    for stop in 1..=polls {
        let mut count = 0;
        assert!(matches!(
            PlateDynamics::new(&model, p, &NoPlateLoad, budget(&model), &mut || {
                count += 1;
                count == stop
            }),
            Err(PlateDynamicsError::Cancelled)
        ));
    }
}
