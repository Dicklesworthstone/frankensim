//! Authored contact parameters through the production discrete adjoint.
#[path = "../examples/calibrate_plate_impact.rs"]
#[allow(dead_code)]
mod example;
#[path = "../examples/material_plate_impact.rs"]
#[allow(dead_code)]
mod sourced;

use example::{
    Parameters, add_initial_pullback, adjoint_config, authored_impact, consistent_state,
    driving_load, evaluate, method,
};
use fs_plate::{PlateMesh, PlateModel};
use fs_plate_transient::NoPlateLoad;
use fs_sparse::Coo;
use fs_time::galpha::SecondOrderProblem;
use fs_time::galpha::second_order_adjoint::{
    SecondOrderVjp, trajectory::StructuralTrajectoryModel,
};

const PARAMETERS: Parameters = [1.1, 0.9, 0.04, 0.02, 0.7, 0.2, 700.0, 3.0, 0.005];

fn fixture() -> (PlateMesh, PlateModel, Vec<f64>, Vec<f64>) {
    let mesh =
        PlateMesh::from_unstructured(vec![(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)], vec![[0, 1, 2]])
            .unwrap();
    let mut k = Coo::new(9, 9);
    let mut m = Coo::new(9, 9);
    for i in 0..9 {
        k.push(i, i, 4.0 * (1.0 + i as f64 / 10.0));
        m.push(i, i, 2.0 + i as f64 / 10.0);
    }
    let plate = PlateModel {
        k: k.assemble(),
        m: m.assemble(),
        dof_map: (0..9).map(Some).collect(),
        free: 9,
    };
    let mut q: Vec<_> = (0..10).map(|i| 0.001 * i as f64).collect();
    let mut v: Vec<_> = (0..10).map(|i| 0.002 * i as f64).collect();
    q[9] = 0.04;
    v[9] = 0.15;
    (mesh, plate, q, v)
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 1.0e-8 + 2.0e-5 * a.abs().max(b.abs()),
        "{a:.16e} != {b:.16e}"
    );
}

#[test]
fn step_pullback_includes_mass_stiffness_loss_gap_and_initial_consistency() {
    let (mesh, plate, q, v) = fixture();
    let load = driving_load(&mesh, &plate, PARAMETERS[4]).unwrap();
    let model = authored_impact(&mesh, &plate, PARAMETERS, &load)
        .unwrap()
        .with_contact_calibration()
        .unwrap();
    let initial = consistent_state(&model, &q, &v).unwrap();
    let step = method(10, 0.001);
    let time = step.forcing_time(0.0).unwrap();
    let mut force = vec![0.0; 10];
    model.forcing(time, &mut force).unwrap();
    let (q_seed, v_seed, a_seed) = ([0.2; 10], [-0.3; 10], [0.04; 10]);
    let config = adjoint_config(10);
    let pulled = step
        .step_vjp(
            &initial,
            &model,
            &force,
            (&q_seed, &v_seed, &a_seed),
            &model.mass_preconditioner(),
            config,
            step.adjoint_workspace_components(model.parameter_count(), config)
                .unwrap(),
            &mut || false,
        )
        .unwrap();
    let unchained = pulled.parameters.clone();
    let mut gradient = pulled.parameters;
    let mut applied = vec![f64::NAN; model.parameter_count()];
    model
        .forcing_vjp(time, &pulled.forcing, &mut applied)
        .unwrap();
    for (g, f) in gradient.iter_mut().zip(applied) {
        *g += f;
    }
    add_initial_pullback(&model, &initial, &pulled.initial_a, &mut gradient).unwrap();
    assert!(
        (gradient[5] - unchained[5]).abs() > 1.0e-4,
        "the mass chain through consistent acceleration must be observable"
    );
    let loss = |p: Parameters| {
        let load = driving_load(&mesh, &plate, p[4]).unwrap();
        let model = authored_impact(&mesh, &plate, p, &load)
            .unwrap()
            .with_contact_calibration()
            .unwrap();
        let mut end = consistent_state(&model, &q, &v).unwrap();
        let mut force = vec![0.0; 10];
        model.forcing(time, &mut force).unwrap();
        step.step(&mut end, &model, &force).unwrap();
        dot(&q_seed, &end.q) + dot(&v_seed, &end.v) + dot(&a_seed, &end.a)
    };
    for index in 5..9 {
        let epsilon = 1.0e-5 * PARAMETERS[index].abs().max(1.0e-3);
        let (mut plus, mut minus) = (PARAMETERS, PARAMETERS);
        plus[index] += epsilon;
        minus[index] -= epsilon;
        close(
            gradient[index],
            (loss(plus) - loss(minus)) / (2.0 * epsilon),
        );
    }
}

#[test]
fn trajectory_preserves_applied_load_prefix_and_zero_forcing_suffix() {
    let (mesh, plate, q, v) = fixture();
    let load = driving_load(&mesh, &plate, PARAMETERS[4]).unwrap();
    let model = authored_impact(&mesh, &plate, PARAMETERS, &load)
        .unwrap()
        .with_contact_calibration()
        .unwrap();
    assert_eq!(model.parameter_count(), 9);
    let mut force_seed = vec![1.0; 10];
    force_seed[9] = 19.0;
    let mut force_gradient = vec![f64::NAN; 9];
    model
        .forcing_vjp(0.013, &force_seed, &mut force_gradient)
        .unwrap();
    assert_eq!(&force_gradient[..4], &[0.0; 4]);
    close(force_gradient[4], 1.26);
    assert_eq!(&force_gradient[5..], &[0.0; 4]);
    let targets = [0.1, 0.05, 0.0];
    let run = evaluate(&mesh, &plate, PARAMETERS, &q, &v, 0.001, &targets, true).unwrap();
    let gradient = run.gradient.unwrap();
    for index in 4..9 {
        let epsilon = 1.0e-5 * PARAMETERS[index].abs().max(1.0e-3);
        let (mut plus, mut minus) = (PARAMETERS, PARAMETERS);
        plus[index] += epsilon;
        minus[index] -= epsilon;
        let a = evaluate(&mesh, &plate, plus, &q, &v, 0.001, &targets, false)
            .unwrap()
            .loss;
        let b = evaluate(&mesh, &plate, minus, &q, &v, 0.001, &targets, false)
            .unwrap()
            .loss;
        close(gradient[index], (a - b) / (2.0 * epsilon));
    }
    assert!(
        gradient[4].abs() > 1.0e-7,
        "the original applied load must affect the trajectory"
    );
}

#[test]
fn calibration_is_opt_in_and_material_source_coefficients_stay_owned() {
    let (mesh, plate, q, v) = fixture();
    let load = driving_load(&mesh, &plate, PARAMETERS[4]).unwrap();
    let model = authored_impact(&mesh, &plate, PARAMETERS, &load).unwrap();
    assert_eq!(model.parameter_count(), 5);
    let mut before = vec![0.0; 10];
    model.resistance(&q, &v, &mut before);
    let model = model.with_contact_calibration().unwrap();
    assert_eq!(model.parameter_count(), 9);
    let mut after = vec![0.0; 10];
    model.resistance(&q, &v, &mut after);
    assert_eq!(before, after);
    let source = sourced::synthetic_source("calibration refusal", 1.0e7, 900.0, 0.4).unwrap();
    let (mesh, plate) = sourced::build_plate(&source).unwrap();
    let no_load = NoPlateLoad;
    let model = sourced::build_impact(&mesh, &plate, source, &no_load).unwrap();
    assert_eq!(model.parameter_count(), 4);
    let Err(error) = model.with_contact_calibration() else {
        panic!("source-bound calibration must require material re-resolution");
    };
    assert!(error.to_string().contains("material re-resolution"));
}

#[test]
fn contact_stiffness_gradient_reduces_actual_synthetic_observation_loss() {
    let (mesh, plate, q, v) = fixture();
    let targets = evaluate(&mesh, &plate, PARAMETERS, &q, &v, 0.001, &[0.0; 3], false)
        .unwrap()
        .samples;
    let mut guess = PARAMETERS;
    guess[6] *= 0.7;
    let current = evaluate(&mesh, &plate, guess, &q, &v, 0.001, &targets, true).unwrap();
    let log_gradient = guess[6] * current.gradient.unwrap()[6];
    assert!(log_gradient.is_finite() && log_gradient != 0.0);
    guess[6] *= (-0.05 * log_gradient.signum()).exp();
    let trial = evaluate(&mesh, &plate, guess, &q, &v, 0.001, &targets, false).unwrap();
    assert!(
        trial.loss < current.loss,
        "{} >= {}",
        trial.loss,
        current.loss
    );
}
