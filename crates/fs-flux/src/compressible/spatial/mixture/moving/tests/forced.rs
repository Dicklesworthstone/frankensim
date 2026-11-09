use super::*;

fn material(reference: f64) -> MixtureGas {
    gas(vec![
        source("A", R, 3.5, 0.0, reference),
        source("B", 2.0 * R, 3.5, 0.0, -2.0 * reference),
        source("absent", 0.5 * R, 3.5, 0.0, 3.0 * reference),
    ])
}

fn workspace(flow: &MixtureEuler3d) -> ForcingWorkspace3d {
    flow.forcing_workspace(flow.required_forcing_workspace_bytes().unwrap())
        .unwrap()
}

fn control() -> ForceControl3d {
    ForceControl3d { maximum_mach: 20.0 }
}

fn closed() -> MixtureEuler3d {
    field(
        material(0.0),
        CartesianGrid::new(1, 1, 1, 0.7, 0.6, 0.5).unwrap(),
        WALLS,
        |_| ([0.6, 0.9, 0.0], [0.08, -0.05, 0.03], 1.3),
    )
}

fn source_balance(report: StepReport) {
    assert!(report.maximum_cfl <= StepControl::default().cfl * (1.0 + 1e-13));
    for value in report.conservation_residual.species_masses {
        close(value, 0.0, 5e-12);
    }
    for value in report.conservation_residual.momentum {
        close(value, 0.0, 5e-12);
    }
    close(report.conservation_residual.energy, 0.0, 2e-10);
    close(report.absolute_energy_residual_j, 0.0, 3e-9);
}

#[test]
fn g0_zero_force_and_zero_motion_retain_the_original_transactions() {
    let acceleration = [0.12, -0.09, 0.07];
    let motion = MeshMotion {
        lower_velocity: [-0.02, 0.0, 0.0],
        upper_velocity: [0.04, 0.0, 0.0],
    };
    let mut expected = closed();
    let mut actual = closed();
    let mut work = workspace(&actual);
    let a = expected
        .step_moving(Time::new(0.003), StepControl::default(), motion, || false)
        .unwrap();
    let b = actual
        .step_moving_with_forcing(
            Time::new(0.003),
            StepControl::default(),
            motion,
            BodyForce3d::UniformAcceleration([0.0; 3]),
            control(),
            &mut work,
            || false,
        )
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(expected.states(), actual.states());
    assert_eq!(expected.grid(), actual.grid());
    assert_eq!(expected.boundary_impulses, actual.boundary_impulses);

    let mut expected = closed();
    let mut actual = closed();
    let mut expected_work = workspace(&expected);
    let mut actual_work = workspace(&actual);
    let mut a_polls = 0;
    let a = expected
        .step_with_forcing(
            Time::new(0.003),
            StepControl::default(),
            BodyForce3d::UniformAcceleration(acceleration),
            control(),
            &mut expected_work,
            || {
                a_polls += 1;
                false
            },
        )
        .unwrap();
    let mut b_polls = 0;
    let b = actual
        .step_moving_with_forcing(
            Time::new(0.003),
            StepControl::default(),
            MeshMotion::default(),
            BodyForce3d::UniformAcceleration(acceleration),
            control(),
            &mut actual_work,
            || {
                b_polls += 1;
                false
            },
        )
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(a_polls, b_polls);
    assert_eq!(expected.states(), actual.states());
    assert_eq!(expected.boundary_impulses, actual.boundary_impulses);
}

#[test]
fn g1_translating_forced_mixture_uses_relative_cfl_exact_work_and_original_gauges() {
    let acceleration = [0.12, -0.09, 0.07];
    let initial_velocity = [8.0, -3.0, 1.5];
    let motion = MeshMotion {
        lower_velocity: initial_velocity,
        upper_velocity: initial_velocity,
    };
    for reference in [0.0, 1e6] {
        for mode in 0..4 {
            let mut flow = field(
                material(reference),
                CartesianGrid::new(2, 2, 2, 0.8, 0.9, 1.0).unwrap(),
                PERIODIC,
                |_| ([0.6, 0.9, 0.0], initial_velocity, 1.3),
            );
            let mut work = workspace(&flow);
            let bytes = flow.workspace_bytes() + work.workspace_bytes();
            let initial = flow
                .cell(0)
                .unwrap()
                .primitive(flow.gas(), || false)
                .unwrap();
            let accelerations = vec![acceleration; flow.cell_count()];
            let densities = vec![acceleration.map(|a| 1.5 * a); flow.cell_count()];
            let force = match mode {
                0 => BodyForce3d::UniformAcceleration(acceleration),
                1 => BodyForce3d::Acceleration(&accelerations),
                2 => BodyForce3d::UniformForceDensity(acceleration.map(|a| 1.5 * a)),
                _ => BodyForce3d::ForceDensity(&densities),
            };
            let report = flow
                .step_moving_with_forcing(
                    Time::new(0.05),
                    StepControl::default(),
                    motion,
                    force,
                    control(),
                    &mut work,
                    || false,
                )
                .unwrap();
            source_balance(report);
            assert_eq!(report.elapsed, Time::new(0.05));
            assert_eq!(report.halvings, 0);
            assert_eq!(flow.workspace_bytes() + work.workspace_bytes(), bytes);
            let dt = report.elapsed.value();
            let mass = report.before.mass();
            let expected_work: f64 = (0..3)
                .map(|axis| {
                    mass * acceleration[axis]
                        * dt
                        * (initial_velocity[axis] + 0.5 * acceleration[axis] * dt)
                })
                .sum();
            close(report.source_change.work, expected_work, 3e-11);
            for axis in 0..3 {
                close(
                    report.source_change.impulse[axis],
                    mass * acceleration[axis] * dt,
                    3e-12,
                );
                close(
                    flow.grid().origin()[axis],
                    initial_velocity[axis] * dt,
                    3e-15,
                );
            }
            for i in 0..flow.cell_count() {
                let q = flow.cell(i).unwrap();
                let p = q.primitive(flow.gas(), || false).unwrap();
                for axis in 0..3 {
                    close(
                        p.velocity[axis],
                        initial_velocity[axis] + acceleration[axis] * dt,
                        3e-10,
                    );
                }
                close(p.temperature, initial.temperature, 3e-9);
                close(p.pressure, initial.pressure, 3e-9);
                for (actual, expected) in q.partial_densities.iter().zip([0.6, 0.9, 0.0]) {
                    close(*actual, expected, 3e-13);
                }
                assert_eq!(q.partial_densities[2], 0.0);
            }
        }
    }
}

#[test]
fn g2_all_axis_deformation_distinguishes_mass_acceleration_from_volume_force() {
    let acceleration = [0.12, -0.09, 0.07];
    let force_density = acceleration.map(|a| 1.5 * a);
    for axis in 0..3 {
        for sign in [-1.0, 1.0] {
            for density_force in [false, true] {
                let mut flow = closed();
                let mut work = workspace(&flow);
                let old_volume = flow.grid().cell_volume();
                let initial = flow
                    .cell(0)
                    .unwrap()
                    .primitive(flow.gas(), || false)
                    .unwrap();
                let mut motion = MeshMotion::default();
                motion.lower_velocity[axis] = -0.02 * sign;
                motion.upper_velocity[axis] = 0.04 * sign;
                let force = if density_force {
                    BodyForce3d::UniformForceDensity(force_density)
                } else {
                    BodyForce3d::UniformAcceleration(acceleration)
                };
                let report = flow
                    .step_moving_with_forcing(
                        Time::new(0.003),
                        StepControl::default(),
                        motion,
                        force,
                        control(),
                        &mut work,
                        || false,
                    )
                    .unwrap();
                source_balance(report);
                let dt = report.elapsed.value();
                let new_volume = flow.grid().cell_volume();
                let mass = report.before.mass();
                close(report.after.mass(), mass, 3e-13);
                let first: [f64; 3] = std::array::from_fn(|a| {
                    0.5 * dt
                        * if density_force {
                            force_density[a] * old_volume
                        } else {
                            acceleration[a] * mass
                        }
                });
                let second: [f64; 3] = std::array::from_fn(|a| {
                    0.5 * dt
                        * if density_force {
                            force_density[a] * new_volume
                        } else {
                            acceleration[a] * report.after.mass()
                        }
                });
                let exact_work: f64 = (0..3)
                    .map(|a| {
                        first[a] * (initial.velocity[a] + 0.5 * first[a] / mass)
                            + second[a] * (report.after.momentum[a] / mass - 0.5 * second[a] / mass)
                    })
                    .sum();
                close(report.source_change.work, exact_work, 3e-13);
                for a in 0..3 {
                    close(report.source_change.impulse[a], first[a] + second[a], 3e-14);
                    close(
                        report.after.momentum[a]
                            - report.before.momentum[a]
                            - report.boundary_change.momentum[a],
                        first[a] + second[a],
                        3e-13,
                    );
                }
                if density_force {
                    assert!(
                        (report.source_change.impulse[0] - force_density[0] * old_volume * dt)
                            .abs()
                            > 1e-9
                    );
                }
                for side in BoundarySide::ALL {
                    for face in flow.boundary_impulses(side) {
                        assert_eq!(face.species_masses, &[0.0; 3]);
                    }
                }
            }
        }
    }
}

#[test]
fn g4_forced_moving_refusal_and_late_cancellation_preserve_every_accepted_quantity() {
    let force = BodyForce3d::UniformAcceleration([0.12, -0.09, 0.07]);
    let motion = MeshMotion {
        lower_velocity: [-0.02, 0.0, 0.0],
        upper_velocity: [0.04, 0.0, 0.0],
    };
    let make = || {
        let mut flow = closed();
        let mut work = workspace(&flow);
        flow.step_moving_with_forcing(
            Time::new(0.001),
            StepControl::default(),
            motion,
            force,
            control(),
            &mut work,
            || false,
        )
        .unwrap();
        (flow, work)
    };
    let (mut expected, mut expected_work) = make();
    let mut polls = 0;
    let accepted = expected
        .step_moving_with_forcing(
            Time::new(0.003),
            StepControl::default(),
            motion,
            force,
            control(),
            &mut expected_work,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    for cutoff in [1, polls / 3, 2 * polls / 3, polls] {
        let (mut flow, mut work) = make();
        let states = flow.states().to_vec();
        let grid = *flow.grid();
        let time = flow.time();
        let report = *flow.last_step().unwrap();
        let loads = flow.boundary_impulses.clone();
        let bytes = flow.workspace_bytes() + work.workspace_bytes();
        let mut seen = 0;
        assert_eq!(
            flow.step_moving_with_forcing(
                Time::new(0.003),
                StepControl::default(),
                motion,
                force,
                control(),
                &mut work,
                || {
                    seen += 1;
                    seen == cutoff
                },
            ),
            Err(Error::Cancelled)
        );
        assert_eq!(flow.states(), states);
        assert_eq!(*flow.grid(), grid);
        assert_eq!(flow.time(), time);
        assert_eq!(flow.last_step(), Some(&report));
        assert_eq!(flow.boundary_impulses, loads);
        assert_eq!(flow.workspace_bytes() + work.workspace_bytes(), bytes);
        assert!(
            flow.step_moving_with_forcing(
                Time::new(0.003),
                StepControl::default(),
                motion,
                force,
                ForceControl3d { maximum_mach: 1e-5 },
                &mut work,
                || false,
            )
            .is_err()
        );
        assert_eq!(flow.states(), states);
        assert_eq!(*flow.grid(), grid);
        assert_eq!(flow.last_step(), Some(&report));
        let replay = flow
            .step_moving_with_forcing(
                Time::new(0.003),
                StepControl::default(),
                motion,
                force,
                control(),
                &mut work,
                || false,
            )
            .unwrap();
        assert_eq!(replay, accepted);
        assert_eq!(flow.states(), expected.states());
        assert_eq!(flow.boundary_impulses, expected.boundary_impulses);
    }
}
