use super::*;
use crate::compressible::tests::{close, gas};

mod timed;

const PERIODIC: Boundaries3d = Boundaries3d {
    x: AxisBoundaries::Periodic,
    y: AxisBoundaries::Periodic,
    z: AxisBoundaries::Periodic,
};
fn field(
    shape: [usize; 3],
    spacing: [f64; 3],
    cp_linear: f64,
    initial: impl Fn(usize) -> [f64; 5],
) -> Euler3d {
    let gas = gas(cp_linear);
    let grid = CartesianGrid::new(
        shape[0], shape[1], shape[2], spacing[0], spacing[1], spacing[2],
    )
    .unwrap();
    let states = (0..grid.len())
        .map(|cell| {
            let p = initial(cell);
            State::from_primitive(&gas, p[0], [p[1], p[2], p[3]], p[4]).unwrap()
        })
        .collect();
    Euler3d::new(gas, grid, states, PERIODIC, || false).unwrap()
}
fn uniform() -> Euler3d {
    field([3; 3], [1.0; 3], 0.0, |_| [1.0, 0.0, 0.0, 0.0, 2.0])
}
fn cap() -> ForceControl3d {
    ForceControl3d { maximum_mach: 10.0 }
}
fn balance(report: &StepReport) {
    for r in report.conservation_residual.array() {
        close(r, 0.0, 5e-12);
    }
}
fn impulses(flow: &Euler3d) -> [Vec<BoundaryImpulse>; 6] {
    BoundarySide::ALL.map(|side| flow.boundary_impulses(side).to_vec())
}

#[test]
fn g1_uniform_acceleration_and_force_density_supply_only_momentum_and_kinetic_work() {
    let rho = 1.7;
    let velocity = [0.23, -0.17, 0.11];
    let acceleration = [2.0, -1.0, 0.5];
    let density_force = acceleration.map(|a| rho * a);
    let accelerations = vec![acceleration; 36];
    let force_densities = vec![density_force; 36];
    for force in [
        BodyForce3d::UniformAcceleration(acceleration),
        BodyForce3d::Acceleration(&accelerations),
        BodyForce3d::UniformForceDensity(density_force),
        BodyForce3d::ForceDensity(&force_densities),
    ] {
        let mut flow = field([3, 4, 3], [0.2, 0.3, 0.4], 0.3, |_| {
            [rho, velocity[0], velocity[1], velocity[2], 2.0 * rho]
        });
        let original = flow.states()[0].primitive(flow.gas(), || false).unwrap();
        let before = flow.totals();
        let mut work = flow.forcing_workspace(usize::MAX).unwrap();
        let pointer = work.original.as_ptr();
        let report = flow
            .step_with_forcing(
                Time::new(0.01),
                StepControl::default(),
                force,
                cap(),
                &mut work,
                || false,
            )
            .unwrap();
        assert_eq!(report.elapsed, Time::new(0.01));
        let dt = report.elapsed.value();
        let change = acceleration.map(|a| a * dt);
        let specific_work: f64 = (0..3)
            .map(|a| velocity[a] * change[a] + 0.5 * change[a] * change[a])
            .sum();
        for q in flow.states() {
            let p = q.primitive(flow.gas(), || false).unwrap();
            close(p.density, rho, 2e-13);
            close(p.temperature, original.temperature, 5e-11);
            for axis in 0..3 {
                close(p.velocity[axis], velocity[axis] + change[axis], 2e-13);
            }
        }
        for axis in 0..3 {
            close(
                report.source_change.impulse[axis],
                before.mass * change[axis],
                2e-13,
            );
        }
        close(
            report.source_change.work,
            before.mass * specific_work,
            2e-13,
        );
        close(
            report.after.energy - report.before.energy,
            report.source_change.work,
            3e-12,
        );
        for value in report.boundary_change.array() {
            close(value, 0.0, 2e-13);
        }
        assert_eq!(pointer, work.original.as_ptr());
        balance(&report);
    }
}

#[test]
fn g1_heterogeneous_kick_uses_each_density_and_does_not_turn_rounded_impulse_into_heat() {
    let flow = field([2; 3], [0.2, 0.3, 0.4], 0.2, |i| {
        let rho = 1.0 + 0.1 * i as f64;
        [
            rho,
            0.03 * (i as f64 - 3.0),
            0.02,
            -0.04,
            rho * (2.0 + 0.05 * i as f64),
        ]
    });
    let loads: Vec<_> = (0..8)
        .map(|i| [0.2 + 0.01 * i as f64, -0.1, 0.3 - 0.02 * i as f64])
        .collect();
    for acceleration in [false, true] {
        let force = if acceleration {
            BodyForce3d::Acceleration(&loads)
        } else {
            BodyForce3d::ForceDensity(&loads)
        };
        let mut candidate = flow.states().to_vec();
        let (source, rate) = kick_field(
            flow.gas(),
            flow.grid(),
            &mut candidate,
            force,
            cap(),
            0.07,
            true,
            &mut || false,
        )
        .unwrap();
        assert!(rate > 0.0);
        let mut expected_impulse = [0.0; 3];
        let mut expected_work = 0.0;
        for (i, (before, after)) in flow.states().iter().zip(&candidate).enumerate() {
            let initial = before.primitive(flow.gas(), || false).unwrap();
            let final_state = after.primitive(flow.gas(), || false).unwrap();
            let next: [f64; 3] = core::array::from_fn(|axis| {
                initial.velocity[axis]
                    + 0.07 * loads[i][axis] / if acceleration { 1.0 } else { before.density }
            });
            let kinetic_change = 0.5
                * before.density
                * (0..3)
                    .map(|axis| {
                        next[axis] * next[axis] - initial.velocity[axis] * initial.velocity[axis]
                    })
                    .sum::<f64>();
            close(after.energy - before.energy, kinetic_change, 3e-13);
            close(final_state.temperature, initial.temperature, 5e-11);
            assert_eq!(after.density, before.density);
            for axis in 0..3 {
                close(final_state.velocity[axis], next[axis], 2e-13);
                expected_impulse[axis] += before.density
                    * (next[axis] - initial.velocity[axis])
                    * flow.grid().cell_volume();
            }
            expected_work += kinetic_change * flow.grid().cell_volume();
        }
        for axis in 0..3 {
            close(source.impulse[axis], expected_impulse[axis], 2e-13);
        }
        close(source.work, expected_work, 2e-13);
    }
    let rounded = kick::<Error, 3>(
        BodyForce3d::UniformForceDensity([0.25, 0.0, 0.0]),
        0,
        1.0,
        [1e16, 0.0, 0.0],
        1e32,
        1.0,
    )
    .unwrap();
    assert_eq!(rounded.momentum, [1e16, 0.0, 0.0]);
    assert_eq!(rounded.impulse_density, [0.0; 3]);
    assert_eq!(rounded.work_density, 0.0);
    assert_eq!(rounded.energy.to_bits(), 1e32_f64.to_bits());
}

#[test]
fn g0_zero_force_retains_original_bits_and_admission_respects_shape_and_byte_limits() {
    let make = || {
        field([3; 3], [0.2; 3], 0.2, |i| {
            [1.0 + 0.01 * i as f64, 0.1, -0.04, 0.02, 2.0]
        })
    };
    let zero = vec![[0.0, -0.0, 0.0]; 27];
    for force in [
        BodyForce3d::UniformAcceleration([0.0; 3]),
        BodyForce3d::Acceleration(&zero),
        BodyForce3d::UniformForceDensity([-0.0; 3]),
        BodyForce3d::ForceDensity(&zero),
    ] {
        let mut actual = make();
        let mut original = make();
        let mut work = actual.forcing_workspace(usize::MAX).unwrap();
        for _ in 0..2 {
            let expected = original
                .step(Time::new(0.001), StepControl::default(), || false)
                .unwrap();
            let report = actual
                .step_with_forcing(
                    Time::new(0.001),
                    StepControl::default(),
                    force,
                    cap(),
                    &mut work,
                    || false,
                )
                .unwrap();
            assert_eq!(report, expected);
            assert_eq!(actual.states(), original.states());
            assert_eq!(actual.time(), original.time());
            assert_eq!(impulses(&actual), impulses(&original));
        }
    }
    let mut flow = make();
    let required = flow.required_forcing_workspace_bytes().unwrap();
    assert_eq!(required, 40 * 27);
    assert!(flow.forcing_workspace(required - 1).is_err());
    let mut work = flow.forcing_workspace(required).unwrap();
    assert_eq!(work.workspace_bytes(), required);
    let short = vec![[0.1; 3]; 26];
    let mut nonfinite = vec![[0.0; 3]; 27];
    nonfinite[13][1] = f64::INFINITY;
    let before = flow.states().to_vec();
    for force in [
        BodyForce3d::Acceleration(&short),
        BodyForce3d::ForceDensity(&nonfinite),
        BodyForce3d::UniformAcceleration([f64::NAN, 0.0, 0.0]),
        BodyForce3d::Coriolis {
            angular_velocity: [0.0, f64::INFINITY, 0.0],
            acceleration: [0.0; 3],
        },
        BodyForce3d::CoriolisField {
            angular_velocity: [0.1, 0.2, 0.3],
            acceleration: &short,
        },
        BodyForce3d::CoriolisField {
            angular_velocity: [0.1, 0.2, 0.3],
            acceleration: &nonfinite,
        },
    ] {
        assert!(matches!(
            flow.step_with_forcing(
                Time::new(0.001),
                StepControl::default(),
                force,
                cap(),
                &mut work,
                || false,
            ),
            Err(Error::Invalid(_))
        ));
    }
    for maximum_mach in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            flow.step_with_forcing(
                Time::new(0.001),
                StepControl::default(),
                BodyForce3d::UniformAcceleration([0.1, 0.0, 0.0]),
                ForceControl3d { maximum_mach },
                &mut work,
                || false,
            )
            .is_err()
        );
    }
    let other = field([2; 3], [1.0; 3], 0.0, |_| [1.0, 0.0, 0.0, 0.0, 2.0]);
    let mut wrong_work = other.forcing_workspace(usize::MAX).unwrap();
    assert!(
        flow.step_with_forcing(
            Time::new(0.001),
            StepControl::default(),
            BodyForce3d::UniformAcceleration([0.1; 3]),
            cap(),
            &mut wrong_work,
            || false,
        )
        .is_err()
    );
    assert_eq!(flow.states(), before);
    assert_eq!(flow.time().value(), 0.0);
    assert!(flow.last_step().is_none());
}

#[test]
fn g1_forcing_retries_the_original_state_for_endpoint_cfl_and_authored_mach() {
    let force = BodyForce3d::UniformAcceleration([20.0, 0.0, 0.0]);
    let mut flow = uniform();
    let mut work = flow.forcing_workspace(usize::MAX).unwrap();
    let report = flow
        .step_with_forcing(
            Time::new(1.0),
            StepControl::default(),
            force,
            cap(),
            &mut work,
            || false,
        )
        .unwrap();
    assert!(report.halvings > 0);
    let p = flow.states()[0].primitive(flow.gas(), || false).unwrap();
    let endpoint_rate: f64 = (0..3)
        .map(|axis| (p.velocity[axis].abs() + p.sound_speed) / flow.grid().spacing()[axis])
        .sum();
    close(
        report.maximum_cfl,
        endpoint_rate * report.elapsed.value(),
        2e-12,
    );
    assert!(report.maximum_cfl <= StepControl::default().cfl * (1.0 + 1e-12));
    close(
        report.source_change.impulse[0],
        report.before.mass * 20.0 * report.elapsed.value(),
        3e-13,
    );
    close(p.temperature, 2.0, 5e-11);
    balance(&report);
    let mut fresh = uniform();
    let mut fresh_work = fresh.forcing_workspace(usize::MAX).unwrap();
    let replay = fresh
        .step_with_forcing(
            report.elapsed,
            StepControl::default(),
            force,
            cap(),
            &mut fresh_work,
            || false,
        )
        .unwrap();
    assert_eq!(fresh.states(), flow.states());
    assert_eq!(replay.source_change, report.source_change);
    assert_eq!(replay.boundary_change, report.boundary_change);
    assert_eq!(replay.halvings, 0);

    let mut limited = uniform();
    let mut limited_work = limited.forcing_workspace(usize::MAX).unwrap();
    let before = limited.states().to_vec();
    assert!(matches!(
        limited.step_with_forcing(
            Time::new(0.01),
            StepControl {
                maximum_halvings: 0,
                ..StepControl::default()
            },
            force,
            ForceControl3d { maximum_mach: 0.01 },
            &mut limited_work,
            || false,
        ),
        Err(Error::StepRetries { .. })
    ));
    assert_eq!(limited.states(), before);
    assert_eq!(limited.time().value(), 0.0);
    let result = limited
        .step_with_forcing(
            Time::new(1.0),
            StepControl::default(),
            force,
            ForceControl3d { maximum_mach: 0.3 },
            &mut limited_work,
            || false,
        )
        .unwrap();
    assert!(result.halvings > 0);
    for q in limited.states() {
        let p = q.primitive(limited.gas(), || false).unwrap();
        let speed = p.velocity.into_iter().map(|u| u * u).sum::<f64>().sqrt();
        assert!(speed / p.sound_speed <= 0.3 * (1.0 + 1e-12));
    }
    balance(&result);
}

#[test]
fn g4_g5_cancelled_force_trials_keep_accepted_history_and_retry_without_duplicate_work() {
    let force = BodyForce3d::UniformForceDensity([0.4, -0.2, 0.1]);
    let mut reference = uniform();
    let mut work = reference.forcing_workspace(usize::MAX).unwrap();
    reference
        .step_with_forcing(
            Time::new(0.005),
            StepControl::default(),
            force,
            cap(),
            &mut work,
            || false,
        )
        .unwrap();
    let mut calls = 0;
    let expected = reference
        .step_with_forcing(
            Time::new(0.01),
            StepControl::default(),
            force,
            cap(),
            &mut work,
            || {
                calls += 1;
                false
            },
        )
        .unwrap();
    for cutoff in [1, calls / 2, calls] {
        let mut candidate = uniform();
        let mut work = candidate.forcing_workspace(usize::MAX).unwrap();
        candidate
            .step_with_forcing(
                Time::new(0.005),
                StepControl::default(),
                force,
                cap(),
                &mut work,
                || false,
            )
            .unwrap();
        let before = candidate.states().to_vec();
        let time = candidate.time();
        let history = candidate.last_step().copied();
        let loads = impulses(&candidate);
        let mut seen = 0;
        assert!(matches!(
            candidate.step_with_forcing(
                Time::new(0.01),
                StepControl::default(),
                force,
                cap(),
                &mut work,
                || {
                    seen += 1;
                    seen == cutoff
                },
            ),
            Err(Error::Cancelled)
        ));
        assert_eq!(candidate.states(), before);
        assert_eq!(candidate.time(), time);
        assert_eq!(candidate.last_step().copied(), history);
        assert_eq!(impulses(&candidate), loads);
        let replay = candidate
            .step_with_forcing(
                Time::new(0.01),
                StepControl::default(),
                force,
                cap(),
                &mut work,
                || false,
            )
            .unwrap();
        assert_eq!(replay, expected);
        assert_eq!(candidate.states(), reference.states());
        assert_eq!(impulses(&candidate), impulses(&reference));
    }
}

#[test]
fn g1_arbitrary_axis_coriolis_rotates_without_heating_the_periodic_gas() {
    let velocity = [0.23, -0.17, 0.11];
    let mut flow = field([3; 3], [0.3; 3], 0.3, |_| {
        [1.7, velocity[0], velocity[1], velocity[2], 3.4]
    });
    let initial = flow.states()[0].primitive(flow.gas(), || false).unwrap();
    let energy = flow.states()[0].energy;
    let mut workspace = flow.forcing_workspace(usize::MAX).unwrap();
    let pointer = workspace.original.as_ptr();
    let mut angle = 0.0_f64;
    for _ in 0..256 {
        let report = flow
            .step_with_forcing(
                Time::new(0.002),
                StepControl::default(),
                BodyForce3d::Coriolis {
                    angular_velocity: [1.0, 2.0, 2.0],
                    acceleration: [0.0; 3],
                },
                cap(),
                &mut workspace,
                || false,
            )
            .unwrap();
        // Two midpoint source half steps each turn by -2 atan(|Omega|*dt/2).
        angle -= 4.0 * (1.5 * report.elapsed.value()).atan();
        assert_eq!(report.source_change.work, 0.0);
        balance(&report);
    }
    assert_eq!(workspace.original.as_ptr(), pointer);
    // Independent Rodrigues rotation about the unit vector (1,2,2)/3.
    let axis = [1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0];
    let projection: f64 = (0..3).map(|i| axis[i] * velocity[i]).sum();
    let cross = [
        axis[1] * velocity[2] - axis[2] * velocity[1],
        axis[2] * velocity[0] - axis[0] * velocity[2],
        axis[0] * velocity[1] - axis[1] * velocity[0],
    ];
    let (sine, cosine) = angle.sin_cos();
    for q in flow.states() {
        let p = q.primitive(flow.gas(), || false).unwrap();
        assert_eq!(q.energy.to_bits(), energy.to_bits());
        close(p.temperature, initial.temperature, 5e-11);
        close(p.pressure, initial.pressure, 5e-11);
        for i in 0..3 {
            let expected =
                cosine * velocity[i] + sine * cross[i] + (1.0 - cosine) * projection * axis[i];
            close(p.velocity[i], expected, 2e-12);
        }
    }
}

#[test]
fn g1_coriolis_field_uses_current_density_and_only_mechanical_work() {
    let flow = field([2; 3], [0.2, 0.3, 0.4], 0.2, |i| {
        let rho = 1.0 + 0.1 * i as f64;
        [rho, 0.03 * (i as f64 - 3.0), 0.02, -0.04, 2.0 * rho]
    });
    let acceleration: Vec<_> = (0..8)
        .map(|i| [0.2 + 0.01 * i as f64, -0.1, 0.3])
        .collect();
    let omega = [0.7, -0.3, 0.2];
    let force = BodyForce3d::CoriolisField {
        angular_velocity: omega,
        acceleration: &acceleration,
    };
    assert!(force.admit::<Error>(flow.states().len(), &mut || false).unwrap());
    let mut candidate = flow.states().to_vec();
    let (source, _) = kick_field(
        flow.gas(),
        flow.grid(),
        &mut candidate,
        force,
        cap(),
        0.07,
        true,
        &mut || false,
    )
    .unwrap();
    let mut expected_work = 0.0;
    for (i, (before, after)) in flow.states().iter().zip(&candidate).enumerate() {
        let p0 = before.primitive(flow.gas(), || false).unwrap();
        let p1 = after.primitive(flow.gas(), || false).unwrap();
        let midpoint: [f64; 3] =
            core::array::from_fn(|axis| 0.5 * (p0.velocity[axis] + p1.velocity[axis]));
        let spin = [
            omega[1] * midpoint[2] - omega[2] * midpoint[1],
            omega[2] * midpoint[0] - omega[0] * midpoint[2],
            omega[0] * midpoint[1] - omega[1] * midpoint[0],
        ];
        let work: f64 = (0..3)
            .map(|axis| 0.07 * before.density * acceleration[i][axis] * midpoint[axis])
            .sum();
        expected_work += work * flow.grid().cell_volume();
        assert_eq!(after.density, before.density);
        close(p1.temperature, p0.temperature, 5e-11);
        close(after.energy - before.energy, work, 3e-13);
        for axis in 0..3 {
            close(
                p1.velocity[axis] - p0.velocity[axis],
                0.07 * (acceleration[i][axis] - 2.0 * spin[axis]),
                3e-13,
            );
        }
    }
    close(source.work, expected_work, 3e-13);

    // Vary both the axis and the drive. Density is an independent cell field;
    // the midpoint momentum equation must use all three at the SAME index.
    let rotations: Vec<_> = (0..8)
        .map(|i| [0.1 * i as f64 - 0.3, 0.2, 0.4 - 0.05 * i as f64])
        .collect();
    let force = BodyForce3d::CoriolisSpatial {
        angular_velocity: &rotations,
        acceleration: &acceleration,
    };
    assert!(force.admit::<Error>(8, &mut || false).unwrap());
    candidate.copy_from_slice(flow.states());
    let (source, _) = kick_field(
        flow.gas(),
        flow.grid(),
        &mut candidate,
        force,
        cap(),
        0.07,
        true,
        &mut || false,
    )
    .unwrap();
    let mut expected_work = 0.0;
    for (i, (before, after)) in flow.states().iter().zip(&candidate).enumerate() {
        let p0 = before.primitive(flow.gas(), || false).unwrap();
        let p1 = after.primitive(flow.gas(), || false).unwrap();
        let midpoint: [f64; 3] =
            core::array::from_fn(|axis| 0.5 * (p0.velocity[axis] + p1.velocity[axis]));
        let omega = rotations[i];
        let spin = [
            omega[1] * midpoint[2] - omega[2] * midpoint[1],
            omega[2] * midpoint[0] - omega[0] * midpoint[2],
            omega[0] * midpoint[1] - omega[1] * midpoint[0],
        ];
        let work: f64 = (0..3)
            .map(|axis| 0.07 * before.density * acceleration[i][axis] * midpoint[axis])
            .sum();
        expected_work += work * flow.grid().cell_volume();
        for axis in 0..3 {
            close(
                p1.velocity[axis] - p0.velocity[axis],
                0.07 * (acceleration[i][axis] - 2.0 * spin[axis]),
                3e-13,
            );
        }
        assert_eq!(after.density, before.density);
        close(after.energy - before.energy, work, 3e-13);
        close(p1.temperature, p0.temperature, 5e-11);
    }
    close(source.work, expected_work, 3e-13);
}

#[test]
fn g1_spatial_coriolis_axes_drive_one_conservative_periodic_trajectory() {
    let mut flow = field([3; 3], [1.0; 3], 0.2, |_| [1.3, 0.23, -0.17, 0.11, 2.0]);
    let rotations: Vec<_> = (0..flow.states().len())
        .map(|cell| match cell % 3 {
            0 => [0.7, -0.3, 0.2],
            1 => [-0.1, 0.4, 0.6],
            _ => [0.0; 3],
        })
        .collect();
    let before = flow.totals();
    let original = flow.states().to_vec();
    let force = BodyForce3d::CoriolisSpatial {
        angular_velocity: &rotations,
        acceleration: &[[0.0; 3]],
    };
    let mut workspace = flow.forcing_workspace(usize::MAX).unwrap();
    let pointer = workspace.original.as_ptr();
    for _ in 0..8 {
        let report = flow
            .step_with_forcing(
                Time::new(0.002),
                StepControl::default(),
                force,
                cap(),
                &mut workspace,
                || false,
            )
            .unwrap();
        assert_eq!(report.source_change.work, 0.0);
        close(report.after.mass, before.mass, 2e-12);
        close(report.after.energy, before.energy, 2e-12);
        balance(&report);
    }
    assert!(flow.states().iter().zip(&original).any(|(a, b)| {
        a.momentum
            .iter()
            .zip(b.momentum)
            .any(|(x, y)| (x - y).abs() > 1e-5)
    }));
    assert_eq!(workspace.original.as_ptr(), pointer);
}
