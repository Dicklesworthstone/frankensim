use super::*;
use crate::compressible::tests::{close, gas};

const PERIODIC: Boundaries3d = Boundaries3d {
    x: AxisBoundaries::Periodic,
    y: AxisBoundaries::Periodic,
    z: AxisBoundaries::Periodic,
};
const OPEN: AxisBoundaries = AxisBoundaries::Ends {
    lower: BoundaryCondition::Transmissive,
    upper: BoundaryCondition::Transmissive,
};
const OPENS: Boundaries3d = Boundaries3d {
    x: OPEN,
    y: OPEN,
    z: OPEN,
};

fn field(
    grid: CartesianGrid,
    boundaries: Boundaries3d,
    sample: impl Fn(usize) -> (f64, [f64; 3], f64),
) -> Euler3d {
    let source = gas(0.0);
    let states = (0..grid.len())
        .map(|cell| {
            let (rho, u, p) = sample(cell);
            State::from_primitive(&source, rho, u, p).unwrap()
        })
        .collect();
    Euler3d::new(source, grid, states, boundaries, || false).unwrap()
}
fn cap() -> ForceControl3d {
    ForceControl3d {
        maximum_mach: 100.0,
    }
}
fn translation(u: [f64; 3]) -> MeshMotion {
    MeshMotion {
        lower_velocity: u,
        upper_velocity: u,
    }
}
fn check(report: StepReport) {
    assert!(report.maximum_cfl <= 0.40000000000001);
    for value in report.conservation_residual.array() {
        close(value, 0.0, 2e-11);
    }
}

#[test]
fn g1_moving_body_acceleration_has_exact_vector_impulse_kinetic_work_and_boost_law() {
    let grid = CartesianGrid::new(3, 2, 2, 0.3, 0.4, 0.5).unwrap();
    let u = [0.23, -0.17, 0.11];
    let a = [0.7, -0.4, 0.2];
    let boost = [4.0, -3.0, 2.0];
    let mut reports = Vec::new();
    for shift in [[0.0; 3], boost] {
        let velocity = std::array::from_fn(|d| u[d] + shift[d]);
        let mut flow = field(grid, PERIODIC, |_| (1.7, velocity, 3.4));
        let initial = flow.states[0].primitive(&flow.gas, || false).unwrap();
        let mut work = flow.forcing_workspace(usize::MAX).unwrap();
        let report = flow
            .step_moving_with_forcing(
                Time::new(0.005),
                StepControl::default(),
                translation(velocity),
                BodyForce3d::UniformAcceleration(a),
                cap(),
                &mut work,
                || false,
            )
            .unwrap();
        check(report);
        assert_eq!(report.elapsed, Time::new(0.005));
        let dt = report.elapsed.value();
        for d in 0..3 {
            close(
                report.source_change.impulse[d],
                report.before.mass * a[d] * dt,
                5e-13,
            );
        }
        let kinetic = report.before.mass
            * (0..3)
                .map(|d| (velocity[d] + 0.5 * a[d] * dt) * a[d] * dt)
                .sum::<f64>();
        close(report.source_change.work, kinetic, 5e-13);
        for q in &flow.states {
            let p = q.primitive(&flow.gas, || false).unwrap();
            close(p.temperature, initial.temperature, 2e-10);
            close(p.density, initial.density, 3e-13);
            for d in 0..3 {
                close(p.velocity[d], velocity[d] + a[d] * dt, 5e-13);
            }
        }
        for q in report.boundary_change.array() {
            close(q, 0.0, 2e-12);
        }
        reports.push(report);
    }
    for d in 0..3 {
        close(
            reports[0].source_change.impulse[d],
            reports[1].source_change.impulse[d],
            3e-13,
        );
    }
    close(
        reports[1].source_change.work - reports[0].source_change.work,
        (0..3)
            .map(|d| boost[d] * reports[0].source_change.impulse[d])
            .sum(),
        3e-13,
    );
    close(reports[0].maximum_cfl, reports[1].maximum_cfl, 3e-12);

    // A large uniform laboratory translation does not consume a lab-frame
    // transport CFL budget; the independent physical Mach cap still applies.
    let mut flow = field(grid, PERIODIC, |_| (1.0, boost, 2.0));
    let mut work = flow.forcing_workspace(usize::MAX).unwrap();
    let report = flow
        .step_moving_with_forcing(
            Time::new(0.025),
            StepControl::default(),
            translation(boost),
            BodyForce3d::UniformAcceleration(a),
            cap(),
            &mut work,
            || false,
        )
        .unwrap();
    assert_eq!(report.elapsed, Time::new(0.025));
    check(report);
}

#[test]
fn g1_moving_force_density_integrates_actual_old_and_new_volumes_on_every_axis() {
    let grid = CartesianGrid::new(3, 2, 2, 0.3, 0.4, 0.5).unwrap();
    let rho = 1.3;
    let u = [0.23, -0.17, 0.11];
    let f = [0.7, -0.4, 0.2];
    let per_cell = vec![f; grid.len()];
    for axis in 0..3 {
        for dilation in [-0.3, 0.3] {
            for force in [
                BodyForce3d::UniformForceDensity(f),
                BodyForce3d::ForceDensity(&per_cell),
            ] {
                let mut flow = field(grid, OPENS, |_| (rho, u, 2.0 * rho));
                let initial = flow.states[0].primitive(&flow.gas, || false).unwrap();
                let mut motion = translation([0.07, -0.02, 0.04]);
                motion.upper_velocity[axis] += dilation;
                let mut work = flow.forcing_workspace(usize::MAX).unwrap();
                let report = flow
                    .step_moving_with_forcing(
                        Time::new(0.005),
                        StepControl::default(),
                        motion,
                        force,
                        cap(),
                        &mut work,
                        || false,
                    )
                    .unwrap();
                check(report);
                let dt = report.elapsed.value();
                let v0 = grid.len() as f64 * grid.cell_volume();
                let v1 = flow.grid.len() as f64 * flow.grid.cell_volume();
                for d in 0..3 {
                    close(
                        report.source_change.impulse[d],
                        0.5 * dt * f[d] * (v0 + v1),
                        5e-13,
                    );
                }
                let expected_work = (0..3)
                    .map(|d| {
                        0.5 * dt
                            * f[d]
                            * (v0 * (u[d] + 0.25 * dt * f[d] / rho)
                                + v1 * (u[d] + 0.75 * dt * f[d] / rho))
                    })
                    .sum::<f64>();
                close(report.source_change.work, expected_work, 5e-13);
                for q in &flow.states {
                    let p = q.primitive(&flow.gas, || false).unwrap();
                    close(p.temperature, initial.temperature, 3e-11);
                    for d in 0..3 {
                        close(p.velocity[d], u[d] + dt * f[d] / rho, 5e-13);
                    }
                }
                // The open, swept boundary carries the first-kicked uniform
                // state; its transport is distinct from both force half-kicks.
                close(report.boundary_change.mass, rho * (v1 - v0), 3e-13);
                for d in 0..3 {
                    close(
                        report.boundary_change.momentum[d],
                        rho * (u[d] + 0.5 * dt * f[d] / rho) * (v1 - v0),
                        3e-13,
                    );
                }
            }
        }
    }
}

#[test]
fn g0_moving_body_zero_force_and_zero_motion_keep_original_transactions_exact() {
    let grid = CartesianGrid::new(4, 2, 2, 0.3, 0.4, 0.5).unwrap();
    let sample = |cell: usize| (1.0 + 0.01 * cell as f64, [0.2, -0.1, 0.3], 2.0);
    for stationary in [false, true] {
        let mut baseline = field(grid, PERIODIC, sample);
        let mut forced = field(grid, PERIODIC, sample);
        let motion = if stationary {
            MeshMotion::default()
        } else {
            translation([0.2, -0.1, 0.3])
        };
        let force = BodyForce3d::UniformAcceleration(if stationary {
            [0.4, -0.2, 0.1]
        } else {
            [0.0; 3]
        });
        let mut a = baseline.forcing_workspace(usize::MAX).unwrap();
        let mut b = forced.forcing_workspace(usize::MAX).unwrap();
        let expected = if stationary {
            baseline.step_with_forcing(
                Time::new(0.005),
                StepControl::default(),
                force,
                cap(),
                &mut a,
                || false,
            )
        } else {
            baseline.step_moving(Time::new(0.005), StepControl::default(), motion, || false)
        }
        .unwrap();
        let actual = forced
            .step_moving_with_forcing(
                Time::new(0.005),
                StepControl::default(),
                motion,
                force,
                cap(),
                &mut b,
                || false,
            )
            .unwrap();
        assert_eq!(actual, expected);
        assert_eq!(forced.states, baseline.states);
        assert_eq!(forced.grid, baseline.grid);
        assert_eq!(forced.boundary_impulses, baseline.boundary_impulses);
        assert_eq!(forced.last_step, baseline.last_step);
    }
}

#[test]
fn g1_moving_body_checks_reconstructed_and_final_kicked_physical_mach() {
    let grid = CartesianGrid::new(3, 1, 1, 1.0, 1.0, 1.0).unwrap();
    let mut flow = field(grid, OPENS, |cell| {
        let u = 0.5 + 0.5 * cell as f64;
        (1.0, [u, 0.0, 0.0], u * u)
    });
    let control = ForceControl3d { maximum_mach: 0.87 };
    for q in &flow.states {
        let p = q.primitive(&flow.gas, || false).unwrap();
        assert!(p.velocity[0] / p.sound_speed < control.maximum_mach);
    }
    // The centre cell's left MC trace is u=0.75, p=0.5:
    // Mach=0.75/sqrt(1.4*0.5)>0.87 despite every centre being admitted.
    assert!(0.75 / (1.4_f64 * 0.5).sqrt() > control.maximum_mach);
    let mut work = flow.forcing_workspace(usize::MAX).unwrap();
    let report = flow
        .step_moving_with_forcing(
            Time::new(0.001),
            StepControl::default(),
            translation([0.2, 0.0, 0.0]),
            BodyForce3d::UniformAcceleration([0.0, 1e-8, 0.0]),
            control,
            &mut work,
            || false,
        )
        .unwrap();
    assert!(report.first_order_retry);
    check(report);

    let make = || field(grid, PERIODIC, |_| (1.0, [0.0; 3], 2.0));
    let mut flow = make();
    let original = flow.states.clone();
    let mut work = flow.forcing_workspace(usize::MAX).unwrap();
    let cap = ForceControl3d { maximum_mach: 0.2 };
    let force = BodyForce3d::UniformAcceleration([20.0, 0.0, 0.0]);
    let motion = translation([0.1, 0.0, 0.0]);
    let refusal = flow.step_moving_with_forcing(
        Time::new(0.02),
        StepControl {
            maximum_halvings: 0,
            ..StepControl::default()
        },
        motion,
        force,
        cap,
        &mut work,
        || false,
    );
    assert!(refusal.is_err());
    assert_eq!(flow.states, original);
    assert_eq!(flow.time, 0.0);
    assert_eq!(flow.grid, grid);
    assert!(flow.last_step.is_none());
    let report = flow
        .step_moving_with_forcing(
            Time::new(0.02),
            StepControl::default(),
            motion,
            force,
            cap,
            &mut work,
            || false,
        )
        .unwrap();
    assert!(report.halvings > 0);
    for q in &flow.states {
        let p = q.primitive(&flow.gas, || false).unwrap();
        cap.admit_velocity::<Error>(p.velocity, p.sound_speed)
            .unwrap();
    }
    check(report);
}

#[test]
fn g1_moving_body_final_cancellation_preserves_geometry_loads_and_retry_receipt() {
    let grid = CartesianGrid::new(3, 2, 2, 0.3, 0.4, 0.5).unwrap();
    let make = || {
        field(grid, OPENS, |cell| {
            (1.0 + 0.01 * cell as f64, [0.1, -0.2, 0.3], 2.0)
        })
    };
    let mut motion = translation([0.03, -0.02, 0.01]);
    motion.upper_velocity[0] += 0.2;
    let force = BodyForce3d::UniformForceDensity([0.3, -0.2, 0.4]);
    let mut probe = make();
    let mut work = probe.forcing_workspace(usize::MAX).unwrap();
    let mut calls = 0;
    let expected = probe
        .step_moving_with_forcing(
            Time::new(0.005),
            StepControl::default(),
            motion,
            force,
            cap(),
            &mut work,
            || {
                calls += 1;
                false
            },
        )
        .unwrap();
    for cutoff in [calls / 2, calls] {
        let mut flow = make();
        let states = flow.states.clone();
        let loads = flow.boundary_impulses.clone();
        let mut work = flow.forcing_workspace(usize::MAX).unwrap();
        let mut count = 0;
        assert_eq!(
            flow.step_moving_with_forcing(
                Time::new(0.005),
                StepControl::default(),
                motion,
                force,
                cap(),
                &mut work,
                || {
                    count += 1;
                    count == cutoff
                },
            ),
            Err(Error::Cancelled)
        );
        assert_eq!(flow.states, states);
        assert_eq!(flow.grid, grid);
        assert_eq!(flow.boundary_impulses, loads);
        assert_eq!(flow.time, 0.0);
        assert!(flow.last_step.is_none());
        let replay = flow
            .step_moving_with_forcing(
                Time::new(0.005),
                StepControl::default(),
                motion,
                force,
                cap(),
                &mut work,
                || false,
            )
            .unwrap();
        assert_eq!(replay, expected);
        assert_eq!(flow.states, probe.states);
    }
}
