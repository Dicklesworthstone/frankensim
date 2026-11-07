//! Original-source real-gas motion, pressure work and whole-state publication.

use super::*;
use fs_gas::body_force::{BodyForce, ForceControl};
use fs_gas::dynamics::nonideal::RealGasWallControl;
use fs_gas::material::nonideal_shock::RealGasShockControl;

fn wall_controls() -> RealGasWallControl {
    RealGasWallControl {
        compression: RealGasShockControl {
            energy: domain().energy,
            compression_ratio_bracket: [1.0 + 1e-8, 3.1],
            mass_flux_tolerance_kg_m2_s: 1e-8,
            momentum_tolerance_pa: 1e-4,
            jump_energy_tolerance_j_kg: 1e-5,
            entropy_tolerance_j_kg_k: 1e-8,
            maximum_density_iterations: 128,
        },
        minimum_rarefaction_density_kg_m3: 6.0,
        velocity_tolerance_m_s: 2e-6,
        integral_tolerance_m_s: 1e-7,
        maximum_quadrature_depth: 10,
    }
}

fn boundaries(physical: bool) -> RealGasBoundaries {
    let boundary = if physical {
        RealGasBoundary::PhysicalReflecting(wall_controls())
    } else {
        RealGasBoundary::Reflecting
    };
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

fn field(gas: RealGas, states: Vec<Conserved>, physical: bool) -> RealGasEuler1d {
    RealGasEuler1d::new(
        gas,
        Length::new(1.0 / states.len() as f64),
        0.07,
        states,
        boundaries(physical),
        || false,
    )
    .unwrap()
}

#[test]
fn original_real_gas_moving_fluxes_are_galilean_and_wall_work_is_pressure_times_speed() {
    let gas = gas();
    let make = |boost| {
        [
            state(&gas, 8.0, 30.0 + boost, 500.0),
            state(&gas, 12.0, -20.0 + boost, 480.0),
        ]
    };
    let q = make(0.0);
    let f = gas
        .moving_hllc_flux(q[0], q[1], Velocity::new(17.0), || false)
        .unwrap();
    let boost = 85.0;
    let shifted = make(boost);
    let g = gas
        .moving_hllc_flux(shifted[0], shifted[1], Velocity::new(17.0 + boost), || {
            false
        })
        .unwrap();
    near(g.mass, f.mass, 1e-7);
    near(g.momentum, f.momentum + boost * f.mass, 1e-3);
    near(
        g.energy,
        f.energy + boost * f.momentum + 0.5 * boost * boost * f.mass,
        0.1,
    );
    near(g.signal_speed, f.signal_speed, 1e-7);

    let p = gas
        .reflecting_wall_flux(state(&gas, 10.0, 30.0, 400.0), 1.0, wall_controls(), || {
            false
        })
        .unwrap();
    let w = gas
        .moving_reflecting_wall_flux(
            state(&gas, 10.0, 35.0, 400.0),
            1.0,
            Velocity::new(5.0),
            wall_controls(),
            || false,
        )
        .unwrap();
    assert_eq!(w.mass, 0.0);
    near(w.momentum, p.momentum, 1e-6);
    near(w.signal_speed, p.signal_speed, 1e-8);
    assert_eq!(w.energy, 5.0 * w.momentum);
}

#[test]
fn comoving_uniform_real_gas_keeps_state_and_translates_nonuniform_geometry() {
    for physical in [false, true] {
        let gas = gas();
        let positions = vec![0.0, 0.08, 0.25, 0.5, 0.75, 0.92, 1.0];
        let volumes = positions.windows(2).map(|x| 0.07 * (x[1] - x[0])).collect();
        let mesh = Mesh1d::new(positions, vec![0.07; 7], volumes).unwrap();
        let q = state(&gas, 10.0, 14.0, 400.0);
        let mut flow =
            RealGasEuler1d::new_on_mesh(gas, mesh, vec![q; 6], boundaries(physical), || false)
                .unwrap();
        let volumes = flow.mesh().volumes().to_vec();
        for _ in 0..6 {
            let old = flow.mesh().positions().to_vec();
            let report = flow
                .step_with_moving_walls(Time::new(2e-6), muscl(), [14.0; 2], || false)
                .unwrap();
            check_balance(report);
            assert_eq!(flow.states(), &[q; 6]);
            assert_eq!(flow.mesh().volumes(), volumes);
            assert_eq!(report.boundary_change, Default::default());
            for (a, b) in flow.mesh().positions().iter().zip(old) {
                assert_eq!(*a, b + 14.0 * report.elapsed.value());
            }
        }
    }
}

#[test]
fn physical_moving_compression_closes_mass_lab_momentum_and_pressure_work() {
    let gas = gas();
    let q = state(&gas, 10.0, 0.0, 450.0);
    let mut flow = field(gas, vec![q; 8], true);
    let before = flow.totals(|| false).unwrap();
    let mut impulse = 0.0;
    let mut work = 0.0;
    for _ in 0..8 {
        let report = flow
            .step_with_moving_walls(Time::new(2e-6), muscl(), [0.0, -12.0], || false)
            .unwrap();
        check_balance(report);
        assert_eq!(report.boundary_change.mass, 0.0);
        assert!(report.boundary_change.energy > 0.0);
        near(
            report.right_boundary_change.energy,
            -12.0 * report.right_boundary_change.momentum,
            2e-12,
        );
        assert_eq!(report.left_boundary_change.energy, 0.0);
        assert_eq!(report.wall_heat, 0.0);
        assert_eq!(report.wall_work, 0.0);
        impulse += report.boundary_change.momentum;
        work += report.boundary_change.energy;
    }
    let after = flow.totals(|| false).unwrap();
    near(after.mass, before.mass, 2e-14);
    near(after.momentum - before.momentum, impulse, 2e-12);
    near(after.energy - before.energy, work, 4e-9);
    near(
        flow.mesh().positions()[8],
        1.0 - 12.0 * flow.time().value(),
        5e-15,
    );
}

#[test]
fn moving_force_density_uses_both_original_and_displaced_volumes() {
    let gas = gas();
    let states = (0..8)
        .map(|i| state(&gas, 10.0 + 0.1 * i as f64, 3.0, 430.0))
        .collect();
    let mut flow = field(gas, states, false);
    let v0: f64 = flow.mesh().volumes().iter().sum();
    let force = 900_000.0;
    let report = flow
        .step_with_moving_walls_and_body_force(
            Time::new(1e-5),
            muscl(),
            [0.0, -3.0],
            BodyForce::UniformForceDensity(force),
            ForceControl { maximum_mach: 4.0 },
            || false,
        )
        .unwrap();
    let v1: f64 = flow.mesh().volumes().iter().sum();
    check_balance(report);
    near(
        report.source_change.impulse,
        0.5 * report.elapsed.value() * force * (v0 + v1),
        3e-12,
    );
    near(
        report.after.energy - report.before.energy,
        report.boundary_change.energy + report.source_change.work,
        5e-9,
    );
    assert!(v1 < v0);
    assert_eq!(report.wall_heat, 0.0);
}

#[test]
fn zero_motion_reuses_original_forced_and_unforced_paths() {
    for forced in [false, true] {
        let gas = gas();
        let states = (0..6)
            .map(|i| state(&gas, 10.0, if i < 3 { 15.0 } else { -5.0 }, 430.0))
            .collect::<Vec<_>>();
        let mut a = field(gas.clone(), states.clone(), false);
        let mut b = field(gas, states, false);
        let (x, y) = if forced {
            let force = BodyForce::UniformAcceleration(1e4);
            let f = ForceControl { maximum_mach: 4.0 };
            (
                a.step_with_body_force(Time::new(2e-6), muscl(), force, f, || false)
                    .unwrap(),
                b.step_with_moving_walls_and_body_force(
                    Time::new(2e-6),
                    muscl(),
                    [0.0; 2],
                    force,
                    f,
                    || false,
                )
                .unwrap(),
            )
        } else {
            (
                a.step(Time::new(2e-6), muscl(), || false).unwrap(),
                b.step_with_moving_walls(Time::new(2e-6), muscl(), [0.0; 2], || false)
                    .unwrap(),
            )
        };
        assert_eq!(x, y);
        assert_eq!(a.states(), b.states());
        assert_eq!(a.primitives(), b.primitives());
        assert_eq!(a.mesh(), b.mesh());
    }
}

#[test]
fn moving_original_energy_gauge_is_invariant_to_large_formation_datums() {
    let original = gas();
    let shifted = bind(
        common::model([1e22, -3e21], [2e18, -7e17], [200.0, 2000.0]),
        domain(),
    );
    let states = (0..6)
        .map(|i| {
            state(
                &original,
                10.0 + 0.1 * i as f64,
                5.0 - 2.0 * i as f64,
                430.0,
            )
        })
        .collect::<Vec<_>>();
    let mut a = field(original, states.clone(), false);
    let mut b = field(shifted, states, false);
    for _ in 0..4 {
        let x = a
            .step_with_moving_walls(Time::new(2e-6), muscl(), [2.0, -4.0], || false)
            .unwrap();
        let y = b
            .step_with_moving_walls(Time::new(2e-6), muscl(), [2.0, -4.0], || false)
            .unwrap();
        assert_eq!(x, y);
        assert_eq!(a.states(), b.states());
        assert_eq!(a.primitives(), b.primitives());
        assert_eq!(a.mesh(), b.mesh());
    }
}

#[test]
fn cancellation_after_prior_motion_preserves_geometry_cache_sources_and_retry() {
    let gas = gas();
    let states = (0..6)
        .map(|i| state(&gas, 10.0 + 0.1 * i as f64, 2.0, 430.0))
        .collect::<Vec<_>>();
    let make = || {
        let mut flow = field(gas.clone(), states.clone(), false);
        flow.step_with_moving_walls(Time::new(2e-6), muscl(), [1.0, -4.0], || false)
            .unwrap();
        flow
    };
    let run = |flow: &mut RealGasEuler1d, cancelled: &mut dyn FnMut() -> bool| {
        flow.step_with_moving_walls_and_body_force(
            Time::new(2e-6),
            muscl(),
            [1.0, -4.0],
            BodyForce::UniformAcceleration(1e4),
            ForceControl { maximum_mach: 4.0 },
            cancelled,
        )
    };
    let mut baseline = make();
    let mut polls = 0;
    let expected = run(&mut baseline, &mut || {
        polls += 1;
        false
    })
    .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut flow = make();
        let q = flow.states().to_vec();
        let p = flow.primitives().to_vec();
        let mesh = flow.mesh().clone();
        let time = flow.time();
        let mut calls = 0;
        assert_eq!(
            run(&mut flow, &mut || {
                calls += 1;
                calls == stop
            }),
            Err(RealGasFlowError::Cancelled)
        );
        assert_eq!(flow.states(), q);
        assert_eq!(flow.primitives(), p);
        assert_eq!(flow.mesh(), &mesh);
        assert_eq!(flow.time(), time);
        let actual = run(&mut flow, &mut || false).unwrap();
        assert_eq!(actual, expected);
        assert_eq!(flow.states(), baseline.states());
        assert_eq!(flow.primitives(), baseline.primitives());
        assert_eq!(flow.mesh(), baseline.mesh());
    }
    let mut unsupported = RealGasEuler1d::new(
        gas,
        Length::new(1.0 / 6.0),
        0.07,
        states.clone(),
        RealGasBoundaries::Periodic,
        || false,
    )
    .unwrap();
    assert!(
        unsupported
            .step_with_moving_walls(Time::new(1e-6), muscl(), [1.0, -1.0], || false)
            .is_err()
    );
    assert_eq!(unsupported.states(), states);
    assert_eq!(unsupported.time().value(), 0.0);
}
