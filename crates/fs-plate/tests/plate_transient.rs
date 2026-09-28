//! G3/G4: actual DKT spatial dynamics, independent dense stepping and adjoints.
#[path = "../examples/moving_plate.rs"]
#[allow(dead_code)] // Share the concrete spatial experiment, not the dense oracle.
mod moving;

use fs_plate::PlateModel;
use fs_plate_transient::{NoPlateLoad, PlateDynamics, PlateDynamicsError, PlateLoad};
use fs_time::galpha::{
    GeneralizedAlpha, SecondOrderProblem, SecondOrderState, galpha_step,
    initialization::{initial_workspace_components, second_order_acceleration},
    second_order_adjoint::trajectory::{
        RecordedStructural, StructuralRecordingConfig, StructuralRecordingStatus,
        StructuralTrajectoryError,
    },
};
use std::{cell::Cell, collections::BTreeSet};

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn close_vector(actual: &[f64], expected: &[f64]) {
    let scale = expected
        .iter()
        .fold(1e-8_f64, |largest, v| largest.max(v.abs()));
    let error = actual
        .iter()
        .zip(expected)
        .fold(0.0_f64, |largest, (a, b)| largest.max((a - b).abs()));
    assert!(
        error <= 2e-9 * scale,
        "vector error {error:.9e}, scale {scale:.9e}"
    );
}

// Independent dense operator construction and LU stepping on the very same
// assembled spatial model. No PlateDynamics action or time-adjoint is used.
fn dense_run(
    mesh: &fs_plate::PlateMesh,
    plate: &PlateModel,
    p: [f64; 5],
) -> (SecondOrderState, f64) {
    let n = plate.free;
    let m = plate
        .m
        .to_dense()
        .iter()
        .map(|value| p[1] * value)
        .collect::<Vec<_>>();
    let k = plate
        .k
        .to_dense()
        .iter()
        .map(|value| p[0] * value)
        .collect::<Vec<_>>();
    let c = m
        .iter()
        .zip(&k)
        .map(|(m, k)| p[2] * m + p[3] * k)
        .collect::<Vec<_>>();
    let dense = GeneralizedAlpha::new(&m, &c, &k, n, moving::STEP, moving::RHO);
    let load = moving::MovingLoad {
        mesh,
        plate,
        amplitude: p[4],
    };
    let mut force = vec![0.0; n];
    load.forcing(0.0, &mut force).unwrap();
    let mut state = SecondOrderState::new(
        0.0,
        &vec![0.0; n],
        &vec![0.0; n],
        &(0..n).map(|i| force[i] / m[i * n + i]).collect::<Vec<_>>(),
    );
    let weights = moving::sensor_weights(mesh, plate).unwrap();
    let alpha_f = moving::RHO / (1.0 + moving::RHO);
    let mut values = vec![0.0];
    for step in 1..=moving::STEPS {
        load.forcing(moving::STEP.mul_add(1.0 - alpha_f, state.t), &mut force)
            .unwrap();
        galpha_step(&dense, &mut state.q, &mut state.v, &mut state.a, &force);
        state.t += moving::STEP;
        state.steps += 1;
        if moving::SAMPLES.contains(&step) {
            values.push(0.5 * (dot(&weights, &state.q) / moving::DEFLECTION_SCALE_M).powi(2));
        }
    }
    (state, values.into_iter().rev().sum())
}

#[test]
fn meshed_moving_force_matches_dense_forward_and_all_five_parameter_derivatives() {
    let (mesh, plate) = moving::build_plate().unwrap();
    assert_eq!(mesh.tris.len(), 32);
    assert_eq!(plate.free, 27);
    let p = [1.1, 0.9, 3.0, 2e-4, 1.7];
    let actual = moving::evaluate(&mesh, &plate, p).unwrap();
    let (expected, value) = dense_run(&mesh, &plate, p);
    close_vector(&actual.final_state.q, &expected.q);
    close_vector(&actual.final_state.v, &expected.v);
    close_vector(&actual.final_state.a, &expected.a);
    assert!((actual.loss - value).abs() < 1e-11 * value);
    assert_eq!(actual.checkpoints, 4);
    assert_eq!(actual.replays, 36);
    let load = moving::MovingLoad {
        mesh: &mesh,
        plate: &plate,
        amplitude: p[4],
    };
    let triangles = (0..=moving::STEPS)
        .map(|i| load.stencil(i as f64 * moving::STEP).unwrap().triangle())
        .collect::<BTreeSet<_>>();
    assert!(
        triangles.len() >= 4,
        "force must cross several actual spatial elements"
    );
    let perturbations = [2e-5, 2e-5, 5e-3, 2e-8, 2e-5];
    for i in 0..5 {
        let (mut plus, mut minus) = (p, p);
        plus[i] += perturbations[i];
        minus[i] -= perturbations[i];
        let fd = (dense_run(&mesh, &plate, plus).1 - dense_run(&mesh, &plate, minus).1)
            / (2.0 * perturbations[i]);
        assert!(
            fd.abs() > 1e-8,
            "parameter {i} must have an informative gradient"
        );
        assert!(
            (actual.gradient[i] - fd).abs() < 5e-5 * fd.abs(),
            "parameter {i}: adjoint {:.12e}, dense FD {fd:.12e}",
            actual.gradient[i]
        );
    }
    // Common scaling of K/M/F leaves this response unchanged. This also states
    // why an experiment must fix mass (or another scale) before fitting all three.
    let common = p[0] * actual.gradient[0] + p[1] * actual.gradient[1] + p[4] * actual.gradient[4];
    assert!(common.abs() < 1e-10);
}

fn energy(plate: &PlateModel, p: [f64; 5], state: &SecondOrderState) -> f64 {
    let mut stiffness = vec![0.0; plate.free];
    let mut mass = vec![0.0; plate.free];
    plate.k.spmv(&state.q, &mut stiffness);
    plate.m.spmv(&state.v, &mut mass);
    0.5 * (p[0] * dot(&state.q, &stiffness) + p[1] * dot(&state.v, &mass))
}

#[test]
fn rho_one_unforced_plate_conserves_or_dissipates_assembled_mechanical_energy() {
    let (mesh, plate) = moving::build_plate().unwrap();
    let mut q = vec![0.0; plate.free];
    // A smooth clamped shape with its actual nodal slopes, not independent
    // arbitrary rotations. The section/mesh still supplies the full DKT energy.
    for (full, reduced) in plate.dof_map.iter().enumerate() {
        if let Some(reduced) = reduced {
            let (x, y) = mesh.nodes[full / 3];
            let shape_x = 16.0 * x * x * (1.0 - x) * (1.0 - x);
            let shape_y = 16.0 * y * y * (0.6 - y) * (0.6 - y) / 0.6_f64.powi(4);
            q[*reduced] = 2e-5
                * match full % 3 {
                    0 => shape_x * shape_y,
                    1 => 32.0 * x * (1.0 - x) * (1.0 - 2.0 * x) * shape_y,
                    _ => shape_x * 32.0 * y * (0.6 - y) * (0.6 - 2.0 * y) / 0.6_f64.powi(4),
                };
        }
    }
    let zero = vec![0.0; plate.free];
    for damping in [[0.0, 0.0], [3.0, 2e-4]] {
        let p = [1.1, 0.9, damping[0], damping[1], 0.0];
        let dynamics = PlateDynamics::new(
            &plate,
            moving::dynamics_parameters(p),
            &NoPlateLoad,
            moving::dynamics_budget(&plate),
            &mut || false,
        )
        .unwrap();
        let mass = dynamics.mass_preconditioner();
        let initial = second_order_acceleration(
            &dynamics,
            &q,
            &zero,
            &zero,
            &mass,
            moving::initial_config(),
            initial_workspace_components(plate.free, 4, moving::initial_config()).unwrap(),
            &mut || false,
        )
        .unwrap();
        let mut state = SecondOrderState::new(0.0, &q, &zero, &initial.value);
        let method = moving::method(plate.free, 1.0);
        let initial_energy = energy(&plate, p, &state);
        assert!(initial_energy > 0.0);
        let mut previous = initial_energy;
        for _ in 0..48 {
            method.step(&mut state, &dynamics, &zero).unwrap();
            let current = energy(&plate, p, &state);
            assert!(current.is_finite() && current >= 0.0);
            assert!(
                current <= previous + 2e-7 * initial_energy,
                "unforced plate gained energy: {previous:.12e} -> {current:.12e}"
            );
            if damping[0] == 0.0 {
                assert!((current - initial_energy).abs() < 2e-7 * initial_energy);
            }
            previous = current;
        }
        if damping[0] > 0.0 {
            assert!(previous < 0.9 * initial_energy);
        }
    }
}

struct InterruptibleLoad(Cell<bool>);
impl PlateLoad for InterruptibleLoad {
    fn parameter_count(&self) -> usize {
        0
    }
    fn forcing(&self, _time: f64, out: &mut [f64]) -> Result<(), String> {
        if self.0.get() {
            return Err("load source unavailable".into());
        }
        out.fill(0.0);
        Ok(())
    }
    fn forcing_vjp(&self, _time: f64, _seed: &[f64], out: &mut [f64]) -> Result<(), String> {
        out.fill(0.0);
        Ok(())
    }
}

#[test]
fn plate_admission_cancellation_and_failed_loading_preserve_the_accepted_state() {
    let (_, plate) = moving::build_plate().unwrap();
    let p = moving::dynamics_parameters([1.1, 0.9, 3.0, 2e-4, 0.0]);
    let budget = moving::dynamics_budget(&plate);
    for refused in [
        fs_plate_transient::PlateDynamicsBudget {
            max_dofs: plate.free - 1,
            ..budget
        },
        fs_plate_transient::PlateDynamicsBudget {
            max_nonzeros: budget.max_nonzeros - 1,
            ..budget
        },
        fs_plate_transient::PlateDynamicsBudget {
            max_full_dofs: budget.max_full_dofs - 1,
            ..budget
        },
    ] {
        assert!(matches!(
            PlateDynamics::new(&plate, p, &NoPlateLoad, refused, &mut || false),
            Err(PlateDynamicsError::Budget(_))
        ));
    }
    let mut polls = 0;
    PlateDynamics::new(&plate, p, &NoPlateLoad, budget, &mut || {
        polls += 1;
        false
    })
    .unwrap();
    for stop in 1..=polls {
        let mut count = 0;
        assert!(matches!(
            PlateDynamics::new(&plate, p, &NoPlateLoad, budget, &mut || {
                count += 1;
                count == stop
            }),
            Err(PlateDynamicsError::Cancelled)
        ));
    }
    assert!(matches!(
        PlateDynamics::new(
            &plate,
            fs_plate_transient::PlateDynamicsParameters {
                mass_scale: 0.0,
                ..p
            },
            &NoPlateLoad,
            budget,
            &mut || false
        ),
        Err(PlateDynamicsError::InvalidInput(_))
    ));
    let load = InterruptibleLoad(Cell::new(false));
    let dynamics = PlateDynamics::new(&plate, p, &load, budget, &mut || false).unwrap();
    let n = plate.free;
    let method = moving::method(n, moving::RHO);
    let adjoint = moving::adjoint_config(n);
    assert!(matches!(
        dynamics.effective_preconditioner(moving::STEP, moving::RHO, &mut || true),
        Err(PlateDynamicsError::Cancelled)
    ));
    let initial = moving::resting_state(&dynamics).unwrap();
    let mut tape = RecordedStructural::new(
        method,
        &dynamics,
        &initial,
        StructuralRecordingConfig {
            steps: 3,
            adjoint,
            max_workspace_components: method.adjoint_workspace_components(4, adjoint).unwrap(),
        },
    )
    .unwrap();
    tape.advance(1, 3, &mut || false).unwrap();
    let before = tape.state().clone();
    load.0.set(true);
    assert!(matches!(
        tape.advance(1, 3, &mut || false),
        Err(StructuralTrajectoryError::Forcing(_))
    ));
    assert_eq!(tape.state(), &before);
    load.0.set(false);
    assert_eq!(
        tape.advance(1, 3, &mut || true).unwrap().status,
        StructuralRecordingStatus::Cancelled
    );
    assert_eq!(tape.state(), &before);
    assert_eq!(
        tape.advance(2, 3, &mut || false).unwrap().status,
        StructuralRecordingStatus::ReachedEnd
    );
    assert_eq!(tape.accepted_steps(), 3);
    let mut action = vec![0.0; n];
    dynamics.mass_apply(&initial.a, &mut action);
    assert!(action.iter().all(|value| *value == 0.0));
}
