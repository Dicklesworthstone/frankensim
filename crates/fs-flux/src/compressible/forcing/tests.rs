use super::*;
use crate::compressible::tests::{close, gas};

const PERIODIC: Boundaries2d = Boundaries2d {
    x: AxisBoundaries::Periodic,
    y: AxisBoundaries::Periodic,
};
const FORCE_CONTROL: ForceControl2d = ForceControl2d { maximum_mach: 8.0 };

fn uniform(density: f64, velocity: [f64; 2]) -> Euler2d {
    let gas = gas(0.4);
    let grid = CartesianGrid::new(5, 4, 0.3, 0.4, 0.7).unwrap();
    let state = State::from_primitive(&gas, density, velocity, 2.3).unwrap();
    Euler2d::new(gas, grid, vec![state; grid.len()], PERIODIC, || false).unwrap()
}

fn workspace(flow: &Euler2d) -> ForcingWorkspace2d {
    flow.forcing_workspace(flow.required_forcing_workspace_bytes().unwrap())
        .unwrap()
}

fn balance(report: &StepReport) {
    assert!(report.maximum_cfl <= 0.40000000000001);
    for residual in report.conservation_residual.array() {
        close(residual, 0.0, 2e-12);
    }
}

#[test]
fn g1_uniform_acceleration_and_force_density_have_distinct_inertia_and_signed_work() {
    // Unequal mass densities must accelerate identically under a, but not F/V.
    // Nonconstant Cp rules out a constant-gamma shortcut in the energy check.
    let velocity = [0.8, -0.3];
    for density in [0.7, 3.4] {
        for forcing in [
            BodyForce2d::UniformAcceleration([-3.0, 1.0]),
            BodyForce2d::UniformForceDensity([-3.0, 1.0]),
        ] {
            let mut flow = uniform(density, velocity);
            let initial = flow.states()[0].primitive(flow.gas(), || false).unwrap();
            let mass = flow.totals().mass;
            let mut scratch = workspace(&flow);
            let report = flow
                .step_with_forcing(
                    Time::new(0.01),
                    StepControl::default(),
                    forcing,
                    FORCE_CONTROL,
                    &mut scratch,
                    || false,
                )
                .unwrap();
            assert_eq!(report.elapsed.value(), 0.01);
            let scale = if matches!(forcing, BodyForce2d::UniformAcceleration(_)) {
                1.0
            } else {
                density.recip()
            };
            let dv = [-0.03 * scale, 0.01 * scale];
            let expected_work = mass
                * (0..2)
                    .map(|axis| dv[axis] * (velocity[axis] + 0.5 * dv[axis]))
                    .sum::<f64>();
            assert!(expected_work < 0.0);
            close(report.source_change.work, expected_work, 3e-13);
            close(
                report.after.energy - report.before.energy,
                expected_work,
                3e-13,
            );
            for axis in 0..2 {
                close(report.source_change.impulse[axis], mass * dv[axis], 3e-13);
            }
            for state in flow.states() {
                let final_state = state.primitive(flow.gas(), || false).unwrap();
                close(final_state.density, density, 3e-14);
                close(final_state.pressure, initial.pressure, 3e-12);
                close(final_state.temperature, initial.temperature, 3e-12);
                for axis in 0..2 {
                    close(final_state.velocity[axis], velocity[axis] + dv[axis], 3e-14);
                }
            }
            balance(&report);
        }
    }
}

#[test]
fn g1_spatial_acceleration_uses_each_current_cell_density() {
    let gas = gas(0.0);
    let grid = CartesianGrid::new(6, 4, 0.25, 0.25, 0.4).unwrap();
    let states: Vec<_> = (0..grid.len())
        .map(|cell| {
            let density = if cell / 6 < 2 { 0.5 } else { 2.0 };
            State::from_primitive(&gas, density, [0.2, 0.0], 1.0).unwrap()
        })
        .collect();
    let forces: Vec<_> = states.iter().map(|q| [q.density * 2.0, 0.0]).collect();
    let accelerations = vec![[2.0, 0.0]; states.len()];
    let mut a = Euler2d::new(gas.clone(), grid, states.clone(), PERIODIC, || false).unwrap();
    let mut f = Euler2d::new(gas, grid, states.clone(), PERIODIC, || false).unwrap();
    let mut wa = workspace(&a);
    let mut wf = workspace(&f);
    let ra = a
        .step_with_forcing(
            Time::new(0.02),
            StepControl::default(),
            BodyForce2d::Acceleration(&accelerations),
            FORCE_CONTROL,
            &mut wa,
            || false,
        )
        .unwrap();
    let rf = f
        .step_with_forcing(
            Time::new(0.02),
            StepControl::default(),
            BodyForce2d::ForceDensity(&forces),
            FORCE_CONTROL,
            &mut wf,
            || false,
        )
        .unwrap();
    for (cell, (actual_a, actual_f)) in a.states().iter().zip(f.states()).enumerate() {
        close(actual_a.density, states[cell].density, 1e-14);
        close(actual_a.momentum[0] / actual_a.density, 0.24, 1e-14);
        close(actual_a.momentum[1], 0.0, 1e-14);
        for (left, right) in actual_a.array().into_iter().zip(actual_f.array()) {
            close(left, right, 2e-14);
        }
    }
    balance(&ra);
    balance(&rf);

    // Now advect the density jump along the accelerated direction. The second
    // kick must use the NEW density, not a force precomputed from the first.
    let gas = crate::compressible::tests::gas(0.0);
    let initial: Vec<_> = (0..grid.len())
        .map(|cell| {
            let density = if cell % 6 < 3 { 0.5 } else { 2.0 };
            State::from_primitive(&gas, density, [0.2, 0.0], 1.0).unwrap()
        })
        .collect();
    let mut moving = Euler2d::new(gas, grid, initial.clone(), PERIODIC, || false).unwrap();
    let mut work = workspace(&moving);
    let report = moving
        .step_with_forcing(
            Time::new(0.02),
            StepControl::default(),
            BodyForce2d::Acceleration(&accelerations),
            FORCE_CONTROL,
            &mut work,
            || false,
        )
        .unwrap();
    assert!(
        moving
            .states()
            .iter()
            .zip(&initial)
            .any(|(new, old)| { (new.density - old.density).abs() > 1e-6 })
    );
    for q in moving.states() {
        close(q.momentum[0] / q.density, 0.24, 2e-13);
        close(q.momentum[1], 0.0, 2e-13);
    }
    close(
        report.source_change.impulse[0],
        report.before.mass * 0.04,
        2e-13,
    );
    balance(&report);
}

#[test]
fn g3_zero_source_preserves_native_bits_and_workspace_budget() {
    let make = || {
        let gas = gas(0.2);
        let grid = CartesianGrid::new(7, 5, 0.2, 0.2, 0.6).unwrap();
        let states = (0..grid.len())
            .map(|cell| {
                let density = if cell % 7 < 3 { 0.8 } else { 1.3 };
                State::from_primitive(&gas, density, [0.2, -0.1], 1.4).unwrap()
            })
            .collect();
        Euler2d::new(gas, grid, states, PERIODIC, || false).unwrap()
    };
    let mut ordinary = make();
    let mut forced = make();
    let bytes = forced.required_forcing_workspace_bytes().unwrap();
    assert!(forced.forcing_workspace(bytes - 1).is_err());
    let mut scratch = forced.forcing_workspace(bytes).unwrap();
    assert_eq!(scratch.workspace_bytes(), bytes);
    let expected = ordinary
        .step(Time::new(0.01), StepControl::default(), || false)
        .unwrap();
    let actual = forced
        .step_with_forcing(
            Time::new(0.01),
            StepControl::default(),
            BodyForce2d::UniformAcceleration([0.0; 2]),
            FORCE_CONTROL,
            &mut scratch,
            || false,
        )
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(forced.states(), ordinary.states());
    for side in [
        BoundarySide::XMin,
        BoundarySide::XMax,
        BoundarySide::YMin,
        BoundarySide::YMax,
    ] {
        assert_eq!(
            forced.boundary_impulses(side),
            ordinary.boundary_impulses(side)
        );
    }
    assert_eq!(actual.source_change, ForceChange2d::default());
}

#[test]
fn g1_acceleration_retry_admits_endpoint_cfl_and_mach() {
    for maximum_mach in [8.0, 0.1] {
        let mut flow = uniform(1.0, [0.0; 2]);
        let mut scratch = workspace(&flow);
        let report = flow
            .step_with_forcing(
                Time::new(0.1),
                StepControl::default(),
                BodyForce2d::UniformAcceleration([80.0, -40.0]),
                ForceControl2d { maximum_mach },
                &mut scratch,
                || false,
            )
            .unwrap();
        assert!(report.halvings > 0);
        let dt = report.elapsed.value();
        for state in flow.states() {
            let p = state.primitive(flow.gas(), || false).unwrap();
            close(p.velocity[0], 80.0 * dt, 2e-14);
            close(p.velocity[1], -40.0 * dt, 2e-14);
            assert!(p.velocity[0].hypot(p.velocity[1]) / p.sound_speed <= maximum_mach);
        }
        balance(&report);
    }
}

#[test]
fn g4_invalid_sources_and_late_cancellation_preserve_complete_accepted_state() {
    let mut flow = uniform(1.2, [0.2, -0.1]);
    let mut scratch = workspace(&flow);
    flow.step(Time::new(0.01), StepControl::default(), || false)
        .unwrap();
    let states = flow.states().to_vec();
    let time = flow.time();
    let report = flow.last_step().copied();
    let boundaries = flow.boundary_impulses.clone();
    for source in [
        BodyForce2d::UniformAcceleration([f64::NAN, 0.0]),
        BodyForce2d::Acceleration(&[]),
        BodyForce2d::UniformForceDensity([f64::INFINITY, 0.0]),
        BodyForce2d::Coriolis {
            angular_velocity: [1.0, 0.0, 0.0],
            acceleration: [0.0; 2],
        },
        BodyForce2d::Coriolis {
            angular_velocity: [0.0, 0.0, f64::NAN],
            acceleration: [0.0; 2],
        },
        BodyForce2d::CoriolisField {
            angular_velocity: [0.0, 0.0, 1.0],
            acceleration: &[],
        },
    ] {
        assert!(
            flow.step_with_forcing(
                Time::new(0.01),
                StepControl::default(),
                source,
                FORCE_CONTROL,
                &mut scratch,
                || false,
            )
            .is_err()
        );
        assert_eq!(flow.states(), states);
        assert_eq!(flow.time(), time);
        assert_eq!(flow.last_step().copied(), report);
        assert_eq!(flow.boundary_impulses, boundaries);
    }
    assert!(matches!(
        flow.step_with_forcing(
            Time::new(0.01),
            StepControl {
                maximum_halvings: 0,
                ..StepControl::default()
            },
            BodyForce2d::UniformAcceleration([200.0, 0.0]),
            ForceControl2d { maximum_mach: 0.5 },
            &mut scratch,
            || false,
        ),
        Err(Error::StepRetries { .. })
    ));
    assert_eq!(flow.states(), states);
    assert_eq!(flow.time(), time);
    assert_eq!(flow.last_step().copied(), report);
    assert_eq!(flow.boundary_impulses, boundaries);
    let source = BodyForce2d::UniformAcceleration([0.4, -0.3]);
    let mut probe = uniform(1.2, [0.2, -0.1]);
    probe
        .step(Time::new(0.01), StepControl::default(), || false)
        .unwrap();
    let mut probe_scratch = workspace(&probe);
    let mut checkpoints = 0;
    let expected = probe
        .step_with_forcing(
            Time::new(0.01),
            StepControl::default(),
            source,
            FORCE_CONTROL,
            &mut probe_scratch,
            || {
                checkpoints += 1;
                false
            },
        )
        .unwrap();
    let mut observed = 0;
    assert!(matches!(
        flow.step_with_forcing(
            Time::new(0.01),
            StepControl::default(),
            source,
            FORCE_CONTROL,
            &mut scratch,
            || {
                observed += 1;
                observed == checkpoints
            },
        ),
        Err(Error::Cancelled)
    ));
    assert_eq!(observed, checkpoints);
    assert_eq!(flow.states(), states);
    assert_eq!(flow.time(), time);
    assert_eq!(flow.last_step().copied(), report);
    assert_eq!(flow.boundary_impulses, boundaries);
    let retry = flow
        .step_with_forcing(
            Time::new(0.01),
            StepControl::default(),
            source,
            FORCE_CONTROL,
            &mut scratch,
            || false,
        )
        .unwrap();
    assert_eq!(retry, expected);
    assert_eq!(flow.states(), probe.states());
}

#[test]
fn g1_coriolis_and_prescribed_acceleration_preserve_uniform_geostrophic_balance() {
    let mut flow = uniform(2.0, [0.5, -0.25]);
    let original = flow.states().to_vec();
    let mut scratch = workspace(&flow);
    // a = 2*Omega cross u: the held pressure-drive acceleration exactly
    // balances the Coriolis acceleration. No centrifugal field is implied.
    let force = BodyForce2d::Coriolis {
        angular_velocity: [0.0, 0.0, 2.0],
        acceleration: [1.0, 2.0],
    };
    for _ in 0..8 {
        let report = flow
            .step_with_forcing(
                Time::new(0.015625),
                StepControl::default(),
                force,
                FORCE_CONTROL,
                &mut scratch,
                || false,
            )
            .unwrap();
        assert_eq!(flow.states(), original);
        assert_eq!(report.source_change, ForceChange2d::default());
        balance(&report);
    }
}

#[test]
fn g3_zero_spin_retains_the_prescribed_acceleration_arithmetic() {
    let acceleration = [0.4, -0.2];
    let mut original = uniform(1.3, [0.2, -0.1]);
    let mut rotating = uniform(1.3, [0.2, -0.1]);
    let mut original_work = workspace(&original);
    let mut rotating_work = workspace(&rotating);
    let expected = original
        .step_with_forcing(
            Time::new(0.01),
            StepControl::default(),
            BodyForce2d::UniformAcceleration(acceleration),
            FORCE_CONTROL,
            &mut original_work,
            || false,
        )
        .unwrap();
    let actual = rotating
        .step_with_forcing(
            Time::new(0.01),
            StepControl::default(),
            BodyForce2d::Coriolis {
                angular_velocity: [0.0; 3],
                acceleration,
            },
            FORCE_CONTROL,
            &mut rotating_work,
            || false,
        )
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(rotating.states(), original.states());
}

#[test]
fn g1_beta_plane_source_uses_each_local_rate_without_heating() {
    let flow = uniform(1.0, [0.6, -0.2]);
    let rotation: Vec<_> = (0..flow.states().len())
        .map(|cell| [0.0, 0.0, 0.75 * ((cell / 5) as f64 - 1.0)])
        .collect();
    let initial: Vec<_> = (0..flow.states().len())
        .map(|cell| {
            State::from_primitive(flow.gas(), 0.7 + 0.1 * cell as f64, [0.6, -0.2], 2.3).unwrap()
        })
        .collect();
    let force = BodyForce2d::CoriolisSpatial {
        angular_velocity: &rotation,
        acceleration: &[[0.0; 2]],
    };
    assert!(force.admit::<Error>(initial.len(), &mut || false).unwrap());
    let mut candidate = initial.clone();
    let dt = 0.07;
    let (change, _) = kick_field(
        flow.gas(),
        flow.grid(),
        &mut candidate,
        force,
        FORCE_CONTROL,
        dt,
        true,
        &mut || false,
    )
    .unwrap();
    assert_eq!(change.work, 0.0);
    for (cell, (before, after)) in initial.iter().zip(&candidate).enumerate() {
        // f(y)=2*Omega_z(y). An independent planar rotation checks both
        // hemispheric signs, the factor two and the exactly stationary equator.
        let angle = -2.0 * (dt * rotation[cell][2]).atan();
        let (sine, cosine) = angle.sin_cos();
        let expected = [0.6 * cosine + 0.2 * sine, 0.6 * sine - 0.2 * cosine];
        assert_eq!(after.density, before.density);
        assert_eq!(after.energy.to_bits(), before.energy.to_bits());
        for axis in 0..2 {
            close(after.momentum[axis] / after.density, expected[axis], 4e-15);
        }
        let speed2: f64 = after
            .momentum
            .into_iter()
            .map(|m| (m / after.density).powi(2))
            .sum();
        close(speed2, 0.4, 4e-15);
        if rotation[cell][2] == 0.0 {
            assert_eq!(after, before);
        }
    }
}

#[test]
fn g3_spatial_coriolis_broadcasts_without_changing_uniform_source_arithmetic() {
    let acceleration = [[0.4, -0.2]];
    for omega in [[0.0; 3], [0.0, 0.0, -1.25]] {
        let rotation = [omega];
        let rotations = vec![omega; 20];
        let accelerations = vec![acceleration[0]; 20];
        for (angular_velocity, acceleration_field) in [
            (&rotation[..], &acceleration[..]),
            (&rotation[..], &accelerations[..]),
            (&rotations[..], &acceleration[..]),
            (&rotations[..], &accelerations[..]),
        ] {
            let mut original = uniform(1.3, [0.2, -0.1]);
            let mut spatial = uniform(1.3, [0.2, -0.1]);
            let mut original_work = workspace(&original);
            let mut spatial_work = workspace(&spatial);
            let pointer = spatial_work.original.as_ptr();
            let expected = original
                .step_with_forcing(
                    Time::new(0.01),
                    StepControl::default(),
                    BodyForce2d::Coriolis {
                        angular_velocity: omega,
                        acceleration: acceleration[0],
                    },
                    FORCE_CONTROL,
                    &mut original_work,
                    || false,
                )
                .unwrap();
            let actual = spatial
                .step_with_forcing(
                    Time::new(0.01),
                    StepControl::default(),
                    BodyForce2d::CoriolisSpatial {
                        angular_velocity,
                        acceleration: acceleration_field,
                    },
                    FORCE_CONTROL,
                    &mut spatial_work,
                    || false,
                )
                .unwrap();
            assert_eq!(spatial.states(), original.states());
            assert_eq!(actual, expected);
            assert_eq!(spatial_work.original.as_ptr(), pointer);
        }
    }
}

#[test]
fn g4_spatial_coriolis_admits_complete_fields_and_cancels_atomically() {
    let mut flow = uniform(1.2, [0.2, -0.1]);
    let rotation: Vec<_> = (0..20)
        .map(|cell| [0.0, 0.0, 0.3 * ((cell / 5) as f64 - 1.0)])
        .collect();
    let mut tilted = rotation.clone();
    tilted[19][0] = 0.1;
    let mut nonfinite = rotation.clone();
    nonfinite[19][2] = f64::NAN;
    let mut bad_acceleration = vec![[0.0; 2]; 20];
    bad_acceleration[19][1] = f64::INFINITY;
    let states = flow.states().to_vec();
    let time = flow.time();
    let history = flow.last_step().copied();
    let boundaries = flow.boundary_impulses.clone();
    let mut scratch = workspace(&flow);
    for (angular_velocity, acceleration) in [
        (&rotation[..19], &[[0.0; 2]][..]),
        (&[][..], &[[0.0; 2]][..]),
        (&tilted[..], &[[0.0; 2]][..]),
        (&nonfinite[..], &[[0.0; 2]][..]),
        (&rotation[..], &[][..]),
        (&rotation[..], &bad_acceleration[..19]),
        (&rotation[..], &bad_acceleration[..]),
    ] {
        assert!(matches!(
            flow.step_with_forcing(
                Time::new(0.01),
                StepControl::default(),
                BodyForce2d::CoriolisSpatial {
                    angular_velocity,
                    acceleration
                },
                FORCE_CONTROL,
                &mut scratch,
                || false,
            ),
            Err(Error::Invalid(_))
        ));
        assert_eq!(flow.states(), states);
        assert_eq!(flow.time(), time);
        assert_eq!(flow.last_step().copied(), history);
        assert_eq!(flow.boundary_impulses, boundaries);
    }
    let source = BodyForce2d::CoriolisSpatial {
        angular_velocity: &rotation,
        acceleration: &[[0.0; 2]],
    };
    let mut probe = uniform(1.2, [0.2, -0.1]);
    let mut probe_work = workspace(&probe);
    let mut calls = 0;
    let expected = probe
        .step_with_forcing(
            Time::new(0.01),
            StepControl::default(),
            source,
            FORCE_CONTROL,
            &mut probe_work,
            || {
                calls += 1;
                false
            },
        )
        .unwrap();
    let mut seen = 0;
    assert!(matches!(
        flow.step_with_forcing(
            Time::new(0.01),
            StepControl::default(),
            source,
            FORCE_CONTROL,
            &mut scratch,
            || {
                seen += 1;
                seen == calls
            },
        ),
        Err(Error::Cancelled)
    ));
    assert_eq!(flow.states(), states);
    assert_eq!(flow.time(), time);
    assert_eq!(flow.last_step().copied(), history);
    assert_eq!(flow.boundary_impulses, boundaries);
    let actual = flow
        .step_with_forcing(
            Time::new(0.01),
            StepControl::default(),
            source,
            FORCE_CONTROL,
            &mut scratch,
            || false,
        )
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(flow.states(), probe.states());
}
