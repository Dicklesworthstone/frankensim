//! G1/G3/G4: actual free plate vibration, exact midpoint oscillator and damage.
#[path = "../../fs-plate/tests/fatigue/fixture.rs"]
#[allow(dead_code)]
mod fixture;
use fixture::*;

use fs_plate::fatigue::{FatiguePlate, PlateFatigueError};
use fs_plate::{AssemblyOptions, EdgeSupport, PlateMesh};
use fs_plate_transient::fatigue::{PlateFatigueInitial, PlateFatigueRun, PlateFatigueRunError};
use fs_plate_transient::{
    NoPlateLoad, PlateDynamics, PlateDynamicsBudget, PlateDynamicsParameters,
};
use fs_solver::NewtonKrylovConfig;
use fs_time::galpha::{
    ImplicitSolveConfig, OperatorGeneralizedAlpha, SecondOrderState,
    initialization::InitialSolveConfig,
};

fn plate(samples: u64) -> FatiguePlate {
    let f = fixture(E, 0.2);
    let mesh = PlateMesh::rectangle(0.02, 0.02, 2, 2);
    let element = mesh.tris.iter().position(|t| t.contains(&4)).unwrap();
    FatiguePlate::assemble(
        &mesh,
        vec![section(&f); mesh.tris.len()],
        &PlateMesh::rectangle_boundary(2, 2),
        &AssemblyOptions {
            pretension: 0.0,
            support: EdgeSupport::Clamped,
        },
        vec![
            gauge(&f, element, H / 2.0, [1.0, 0.0]),
            gauge(&f, element, -H / 2.0, [1.0, 0.0]),
        ],
        config(samples),
        &mut || false,
    )
    .unwrap()
}
fn parameters(damping: f64) -> PlateDynamicsParameters {
    PlateDynamicsParameters {
        stiffness_scale: 1.0,
        mass_scale: 1.0,
        mass_damping_per_s: damping,
        stiffness_damping_s: 0.0,
    }
}
fn budget(plate: &FatiguePlate) -> PlateDynamicsBudget {
    let p = plate.model();
    PlateDynamicsBudget {
        max_dofs: p.free,
        max_nonzeros: p.k.nnz() + p.m.nnz(),
        max_full_dofs: p.dof_map.len(),
        max_load_parameters: 0,
    }
}
fn dynamics(plate: &FatiguePlate, damping: f64) -> PlateDynamics<'_, NoPlateLoad> {
    PlateDynamics::new(
        plate.model(),
        parameters(damping),
        &NoPlateLoad,
        budget(plate),
        &mut || false,
    )
    .unwrap()
}
fn method(n: usize, step: f64) -> OperatorGeneralizedAlpha {
    OperatorGeneralizedAlpha::new(
        n,
        step,
        1.0,
        ImplicitSolveConfig {
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1e-12,
                relative_tolerance: 1e-12,
                linear_restart: n,
                max_linear_cycles: 4,
                forcing_minimum: 1e-13,
                forcing_maximum: 1e-11,
                ..NewtonKrylovConfig::default()
            },
            max_newton_iterations: 8,
        },
    )
}
fn initial<'a>(q: &'a [f64], v: &'a [f64]) -> PlateFatigueInitial<'a> {
    PlateFatigueInitial {
        time_s: 0.0,
        displacement: q,
        velocity: v,
        solve: InitialSolveConfig {
            restart: 3,
            max_cycles: 3,
            tolerance: 1e-13,
        },
        maximum_workspace_components: 10000,
    }
}

#[test]
fn g1_free_vibration_recovers_exact_signed_cycles_and_miner_damage() {
    let plate = plate(1000);
    let p = plate.model();
    assert_eq!(p.free, 3);
    let w = p.dof_map[3 * 4].unwrap();
    let m = p.m.to_dense();
    let k = p.k.to_dense();
    // Symmetry decouples the sole center displacement from its two slopes.
    for j in 0..3 {
        if j != w {
            assert!(k[3 * w + j].abs() < 1e-12 * k[3 * w + w]);
        }
    }
    let omega = (k[3 * w + w] / m[3 * w + w]).sqrt();
    let steps_per_cycle = 64;
    // Midpoint advances a free oscillator by theta=2 atan(omega*h/2).
    let h = 2.0 * (std::f64::consts::PI / steps_per_cycle as f64).tan() / omega;
    let stepper = method(3, h);
    let mut q = vec![0.0; 3];
    q[w] = 1e-7;
    let zero = vec![0.0; 3];
    let mut run = PlateFatigueRun::new(
        &plate,
        dynamics(&plate, 0.0),
        stepper,
        initial(&q, &zero),
        &mut || false,
    )
    .unwrap();
    let amplitude = run.observation().stresses_pa()[0].abs();
    assert!(amplitude > 1.0);
    close(run.state().a[w], -omega * omega * q[w], 1e-12);
    let mut native = run.state().clone();
    let mut checkpoint = None;
    for i in 1..=3 * steps_per_cycle {
        run.advance(&mut || false).unwrap();
        stepper
            .step(&mut native, &dynamics(&plate, 0.0), &zero)
            .unwrap();
        native.history.clear();
        assert_eq!(run.state(), &native);
        let expected =
            q[w] * (2.0 * std::f64::consts::PI * i as f64 / steps_per_cycle as f64).cos();
        assert!((run.state().q[w] - expected).abs() < 2e-9 * q[w]);
        close(
            run.observation().stresses_pa()[0],
            -run.observation().stresses_pa()[1],
            1e-13,
        );
        if i == steps_per_cycle + 11 {
            checkpoint = Some(run.clone());
        }
    }
    let before = run.summaries();
    let preview = run.preview_finish(&mut || false).unwrap();
    assert_eq!(run.summaries(), before);
    for summary in &preview {
        close(summary.miner_damage, 3.0 * amplitude / LIFE_STRESS, 3e-9);
    }
    let mut replay = checkpoint.unwrap();
    while replay.state().steps < run.state().steps {
        replay.advance(&mut || false).unwrap();
    }
    assert_eq!(replay.state(), run.state());
    let finished = run.finish(&mut || false).unwrap();
    assert_eq!(preview, finished);
    assert_eq!(replay.finish(&mut || false).unwrap(), finished);
    assert_eq!(run.finish(&mut || false).unwrap(), finished);
    assert!(matches!(
        run.advance(&mut || false),
        Err(PlateFatigueRunError::Fatigue(PlateFatigueError::Finished))
    ));
}

#[test]
fn g3_physical_damping_reduces_the_same_native_vibration_damage() {
    let plate = plate(1000);
    let p = plate.model();
    let w = p.dof_map[12].unwrap();
    let omega = (p.k.to_dense()[3 * w + w] / p.m.to_dense()[3 * w + w]).sqrt();
    let h = 2.0 * (std::f64::consts::PI / 64.0).tan() / omega;
    let mut q = vec![0.0; 3];
    q[w] = 1e-7;
    let v = vec![0.0; 3];
    let mut bare = PlateFatigueRun::new(
        &plate,
        dynamics(&plate, 0.0),
        method(3, h),
        initial(&q, &v),
        &mut || false,
    )
    .unwrap();
    let mut damped = PlateFatigueRun::new(
        &plate,
        dynamics(&plate, 0.1 * omega),
        method(3, h),
        initial(&q, &v),
        &mut || false,
    )
    .unwrap();
    for _ in 0..192 {
        bare.advance(&mut || false).unwrap();
        damped.advance(&mut || false).unwrap();
    }
    assert!(damped.state().q[w].abs() < 0.5 * bare.state().q[w].abs());
    let a = bare.finish(&mut || false).unwrap()[0].miner_damage;
    let b = damped.finish(&mut || false).unwrap()[0].miner_damage;
    assert!(b > 0.3 * a && b < 0.85 * a, "damped={b}, bare={a}");
}

#[test]
fn g4_sample_budget_and_late_cancellation_rollback_mechanics_and_cycles() {
    let plate = plate(2);
    let mut q = vec![0.0; 3];
    q[plate.model().dof_map[12].unwrap()] = 1e-7;
    let v = vec![0.0; 3];
    let mut run = PlateFatigueRun::new(
        &plate,
        dynamics(&plate, 0.0),
        method(3, 1e-6),
        initial(&q, &v),
        &mut || false,
    )
    .unwrap();
    let mut trial = run.clone();
    let mut calls = 0;
    trial
        .advance(&mut || {
            calls += 1;
            false
        })
        .unwrap();
    let old = run.state().clone();
    let summaries = run.summaries();
    let mut seen = 0;
    assert!(
        run.advance(&mut || {
            seen += 1;
            seen == calls
        })
        .is_err()
    );
    assert_eq!(run.state(), &old);
    assert_eq!(run.summaries(), summaries);
    run.advance(&mut || false).unwrap();
    let old = run.state().clone();
    let summaries = run.summaries();
    assert!(matches!(
        run.advance(&mut || false),
        Err(PlateFatigueRunError::Fatigue(
            PlateFatigueError::Gauge { .. }
        ))
    ));
    assert_eq!(run.state(), &old);
    assert_eq!(run.summaries(), summaries);
    assert!(run.finish(&mut || true).is_err());
    assert_eq!(run.state(), &old);
    assert_eq!(run.summaries(), summaries);
}

#[test]
fn g4_foreign_or_rescaled_mechanics_cannot_change_the_observed_material() {
    let plate = plate(20);
    let q = vec![0.0; 3];
    let v = q.clone();
    let mut scaled = parameters(0.0);
    scaled.stiffness_scale = 2.0;
    let d = PlateDynamics::new(
        plate.model(),
        scaled,
        &NoPlateLoad,
        budget(&plate),
        &mut || false,
    )
    .unwrap();
    assert!(matches!(
        PlateFatigueRun::new(&plate, d, method(3, 1e-6), initial(&q, &v), &mut || false),
        Err(PlateFatigueRunError::Invalid(_))
    ));
    let foreign = plate.model().clone();
    let d = PlateDynamics::new(
        &foreign,
        parameters(0.0),
        &NoPlateLoad,
        budget(&plate),
        &mut || false,
    )
    .unwrap();
    assert!(matches!(
        PlateFatigueRun::new(&plate, d, method(3, 1e-6), initial(&q, &v), &mut || false),
        Err(PlateFatigueRunError::Invalid(_))
    ));
    // A native motion checkpoint remains a plain inspectable physical state.
    let run = PlateFatigueRun::new(
        &plate,
        dynamics(&plate, 0.0),
        method(3, 1e-6),
        initial(&q, &v),
        &mut || false,
    )
    .unwrap();
    assert_eq!(run.state(), &SecondOrderState::new(0.0, &q, &v, &v));
}

#[path = "fatigue/findley.rs"]
mod findley;
