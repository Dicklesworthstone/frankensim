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
