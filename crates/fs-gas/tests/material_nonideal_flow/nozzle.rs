//! Real-gas momentum in a stationary nozzle, using the original NASA/PR owner.

use super::*;
use fs_gas::body_force::{BodyForce, ForceControl};

fn mesh(scale: f64) -> Mesh1d {
    let positions = vec![0.0, 0.08, 0.25, 0.5, 0.75, 0.92, 1.0];
    // A(x)=scale*(0.04+0.03*x^2); retain its exact integrated volumes.
    let areas = positions
        .iter()
        .map(|x| scale * (0.04 + 0.03 * x * x))
        .collect();
    let volumes = positions
        .windows(2)
        .map(|x| scale * (0.04 * (x[1] - x[0]) + 0.01 * (x[1] * x[1] * x[1] - x[0] * x[0] * x[0])))
        .collect();
    Mesh1d::new(positions, areas, volumes).unwrap()
}

fn ends(boundary: RealGasBoundary) -> RealGasBoundaries {
    RealGasBoundaries::Ends {
        left: boundary,
        right: boundary,
    }
}

fn muscl() -> StepControl {
    StepControl {
        reconstruction: Reconstruction::Muscl,
        ..controls()
    }
}

#[test]
fn resting_real_gas_preserves_state_and_balances_end_traction_against_side_wall() {
    let gas = gas();
    let q = state(&gas, 10.0, 0.0, 400.0);
    for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
        let mut flow = RealGasEuler1d::new_on_mesh(
            gas.clone(),
            mesh(1.0),
            vec![q; 6],
            ends(RealGasBoundary::Reflecting),
            || false,
        )
        .unwrap();
        let pressure = flow.primitives()[0].pressure;
        for _ in 0..6 {
            let report = flow
                .step(
                    Time::new(2e-6),
                    StepControl {
                        reconstruction,
                        ..controls()
                    },
                    || false,
                )
                .unwrap();
            check_balance(report);
            assert_eq!(flow.states(), &[q; 6]);
            assert_eq!(report.before, report.after);
            let side = report.elapsed.value() * pressure * (0.07 - 0.04);
            near(report.area_momentum_impulse, side, 2e-15);
            near(report.boundary_change.momentum, -side, 2e-15);
            assert_eq!(report.boundary_change.mass, 0.0);
            assert_eq!(report.boundary_change.energy, 0.0);
            assert_eq!(report.wall_work, 0.0);
        }
    }
}

#[test]
fn closed_compressing_nozzle_closes_gas_impulse_without_side_wall_energy_work() {
    let gas = gas();
    let states = (0..6)
        .map(|i| {
            state(
                &gas,
                10.0 + 0.2 * i as f64,
                if i < 3 { 35.0 } else { -20.0 },
                430.0,
            )
        })
        .collect();
    let mut flow = RealGasEuler1d::new_on_mesh(
        gas,
        mesh(1.0),
        states,
        ends(RealGasBoundary::Reflecting),
        || false,
    )
    .unwrap();
    let before = flow.totals(|| false).unwrap();
    let mut impulse = 0.0;
    for _ in 0..12 {
        let report = flow.step(Time::new(1e-6), muscl(), || false).unwrap();
        check_balance(report);
        impulse += report.boundary_change.momentum + report.area_momentum_impulse;
        assert_eq!(report.boundary_change.mass, 0.0);
        assert_eq!(report.boundary_change.energy, 0.0);
        assert!(report.area_momentum_impulse > 0.0);
    }
    let after = flow.totals(|| false).unwrap();
    near(after.mass, before.mass, 2e-14);
    near(after.energy, before.energy, 3e-9);
    near(after.momentum - before.momentum, impulse, 2e-13);
}

#[test]
fn common_area_scale_preserves_forced_nozzle_states_and_scales_extensive_momentum_work() {
    let gas = gas();
    let states = (0..6)
        .map(|i| {
            state(
                &gas,
                10.0 + 0.1 * i as f64,
                13.0 + 2.0 * i as f64,
                400.0 + 2.0 * i as f64,
            )
        })
        .collect::<Vec<_>>();
    let make = |scale| {
        RealGasEuler1d::new_on_mesh(
            gas.clone(),
            mesh(scale),
            states.clone(),
            ends(RealGasBoundary::Transmissive),
            || false,
        )
        .unwrap()
    };
    let mut a = make(1.0);
    let mut b = make(8.0);
    for _ in 0..4 {
        let step = |flow: &mut RealGasEuler1d| {
            flow.step_with_body_force(
                Time::new(2e-6),
                muscl(),
                BodyForce::UniformAcceleration(1e5),
                ForceControl { maximum_mach: 4.0 },
                || false,
            )
            .unwrap()
        };
        let x = step(&mut a);
        let y = step(&mut b);
        check_balance(x);
        check_balance(y);
        assert_eq!(a.states(), b.states());
        assert_eq!(a.primitives(), b.primitives());
        assert_eq!(x.elapsed, y.elapsed);
        assert_eq!(x.maximum_cfl, y.maximum_cfl);
        for (small, large) in [
            (x.area_momentum_impulse, y.area_momentum_impulse),
            (x.boundary_change.momentum, y.boundary_change.momentum),
            (x.source_change.impulse, y.source_change.impulse),
            (x.source_change.work, y.source_change.work),
            (x.after.energy, y.after.energy),
        ] {
            assert_eq!(8.0 * small, large);
        }
    }
}

#[test]
fn periodic_nozzle_requires_one_shared_seam_and_boundary_refusal_is_atomic() {
    let gas = gas();
    let q = state(&gas, 10.0, 15.0, 400.0);
    assert!(
        RealGasEuler1d::new_on_mesh(
            gas.clone(),
            mesh(1.0),
            vec![q; 6],
            RealGasBoundaries::Periodic,
            || false,
        )
        .is_err()
    );
    let original = ends(RealGasBoundary::Transmissive);
    let mut flow =
        RealGasEuler1d::new_on_mesh(gas.clone(), mesh(1.0), vec![q; 6], original, || false)
            .unwrap();
    let primitive = flow.primitives().to_vec();
    assert!(
        flow.set_boundaries(RealGasBoundaries::Periodic, || false)
            .is_err()
    );
    assert_eq!(flow.boundaries(), original);
    assert_eq!(flow.states(), &[q; 6]);
    assert_eq!(flow.primitives(), primitive);
    assert_eq!(flow.time().value(), 0.0);

    let periodic = Mesh1d::new(
        vec![0.0, 0.5, 1.0],
        vec![0.04, 0.06, 0.04],
        vec![0.025, 0.025],
    )
    .unwrap();
    let mut flow = RealGasEuler1d::new_on_mesh(
        gas,
        periodic,
        vec![q; 2],
        RealGasBoundaries::Periodic,
        || false,
    )
    .unwrap();
    let report = flow.step(Time::new(1e-6), muscl(), || false).unwrap();
    check_balance(report);
    assert_eq!(report.boundary_change, Default::default());
}

#[test]
fn nozzle_cancellation_preserves_the_original_field_and_retry_replays_the_same_step() {
    let gas = gas();
    let states = (0..6)
        .map(|i| state(&gas, 10.0, 10.0 + 3.0 * i as f64, 410.0))
        .collect::<Vec<_>>();
    let make = || {
        RealGasEuler1d::new_on_mesh(
            gas.clone(),
            mesh(1.0),
            states.clone(),
            ends(RealGasBoundary::Reflecting),
            || false,
        )
        .unwrap()
    };
    let mut baseline = make();
    let mut polls = 0;
    let expected = baseline
        .step(Time::new(2e-6), muscl(), || {
            polls += 1;
            false
        })
        .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut flow = make();
        let primitive = flow.primitives().to_vec();
        let mut calls = 0;
        assert_eq!(
            flow.step(Time::new(2e-6), muscl(), || {
                calls += 1;
                calls == stop
            }),
            Err(RealGasFlowError::Cancelled)
        );
        assert_eq!(flow.states(), states);
        assert_eq!(flow.primitives(), primitive);
        assert_eq!(flow.time().value(), 0.0);
        let actual = flow.step(Time::new(2e-6), muscl(), || false).unwrap();
        assert_eq!(actual, expected);
        assert_eq!(flow.states(), baseline.states());
    }
}
