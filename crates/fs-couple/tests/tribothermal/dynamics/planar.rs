//! G1/G3 physical paths through the original material and caloric owners.
use super::*;
use fs_couple::tribothermal::dynamics::planar::PlanarMotionSolve;

fn planar_controls(step: f64) -> PlanarMotionSolve {
    PlanarMotionSolve {
        sliding: motion_controls().sliding,
        maximum_step_s: step,
        maximum_segments: 512,
        maximum_motion_trials: 10_000,
        maximum_halvings: 512,
        maximum_bisections: 96,
        impulse_tolerance_n_s: 1e-8,
        event_time_tolerance_s: 1e-10,
        total_energy_tolerance_j: 1e-6,
    }
}

#[test]
fn arbitrary_tangent_orientation_stops_with_full_physical_heat_and_retained_masses() {
    for direction in [[0.6, 0.8], [-0.8, 0.6], [0.0, -1.0]] {
        let (a, ca, body_a) = sourced_body(2.0, 2.0);
        let (b, cb, body_b) = sourced_body(3.0, 3.0);
        let cards = [a, b];
        let curves = [ca, cb];
        let point = context();
        let thermal = inventory_with_geometry(&curves, [2.0, 3.0], [1.0, 1.0]);
        let card = interface_with_static(
            &cards,
            0.6,
            "speed",
            SPEED,
            [(0.0, 0.2), (2.0, 0.2)],
            2000.0,
        );
        let contact = MaterialSlidingContact::new(
            &thermal,
            binding(
                &cards,
                &curves,
                &card,
                &point,
                HeatPartition::new(0.4, 0.6, 0.0).unwrap(),
            ),
            10.0,
            1.0,
            InputAuthority::SyntheticFixture,
        )
        .unwrap();
        let pair = SlidingPair::from_material_bodies(
            &thermal,
            contact,
            [body_a, body_b],
            [direction, [0.0; 2]],
        )
        .unwrap();
        let next = pair
            .step_planar_candidate(
                &thermal,
                1.0,
                MotionDrive {
                    normal_load_n: 10.0,
                    body_forces_n: [[0.0; 2]; 2],
                },
                planar_controls(1.0),
                &[0.0; 2],
                other,
                || false,
            )
            .unwrap();
        let report = next.report();
        assert_eq!(report.segments.len(), 2);
        assert!(report.segments[0].reached_stop);
        close(report.segments[0].end_time_s, 0.6);
        assert_eq!(
            report.segments[1].sliding.friction.regime,
            fs_tribo::FrictionRegime::Sticking
        );
        for velocity in next.pair().velocities_mps() {
            for axis in 0..2 {
                close(velocity[axis], 0.4 * direction[axis]);
            }
        }
        close(report.kinetic_energy_change_j, -0.6);
        close(report.deposited_heat_j, 0.6);
        close(report.thermal_energy_change_j, 0.6);
        close(report.external_work_j, 0.0);
        close(report.total_energy_residual_j, 0.0);
        close(next.pair().contact().accumulated_distance_m(), 0.3);
        for axis in 0..2 {
            close(
                next.pair().contact().state().slip_m()[axis],
                0.3 * direction[axis],
            );
        }
        close(next.thermal().temperature_k(CELLS[0]).unwrap(), 300.06);
        close(next.thermal().temperature_k(CELLS[1]).unwrap(), 300.04);
        assert!(next.pair().material_bodies().is_some());
        assert_eq!(pair.velocities_mps(), [direction, [0.0; 2]]);
        assert_eq!(pair.contact().state().slip_m(), [0.0; 2]);
    }
}

#[test]
fn oblique_stop_then_reversal_preserves_signed_motion_and_unsigned_heat() {
    let (a, ca) = common::material(1.0, 2.0, 0.0);
    let (b, cb) = common::material(1.0, 3.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let point = context();
    let thermal = inventory(&curves, [1.0, 1.0]);
    let card = interface_with_static(
        &cards,
        0.2,
        "speed",
        SPEED,
        [(0.0, 0.2), (4.0, 0.2)],
        2000.0,
    );
    let contact = MaterialSlidingContact::new(
        &thermal,
        binding(
            &cards,
            &curves,
            &card,
            &point,
            HeatPartition::new(0.5, 0.5, 0.0).unwrap(),
        ),
        10.0,
        1.0,
        InputAuthority::SyntheticFixture,
    )
    .unwrap();
    let direction = [0.6, 0.8];
    let pair = SlidingPair::new(&thermal, contact, [direction, [0.0; 2]]).unwrap();
    let next = pair
        .step_planar_candidate(
            &thermal,
            0.5,
            MotionDrive {
                normal_load_n: 10.0,
                body_forces_n: [[-6.0 * direction[0], -6.0 * direction[1]], [0.0; 2]],
            },
            planar_controls(0.5),
            &[0.0; 2],
            other,
            || false,
        )
        .unwrap();
    let report = next.report();
    assert_eq!(report.segments.len(), 2);
    assert!(report.segments[0].reached_stop);
    close(report.segments[0].end_time_s, 0.1);
    for axis in 0..2 {
        close(
            next.pair().velocities_mps()[0][axis],
            -1.4 * direction[axis],
        );
        close(
            next.pair().velocities_mps()[1][axis],
            -0.6 * direction[axis],
        );
        close(
            report.body_displacements_m[0][axis],
            -0.18 * direction[axis],
        );
        close(
            report.body_displacements_m[1][axis],
            -0.07 * direction[axis],
        );
        close(
            next.pair().contact().state().slip_m()[axis],
            -0.11 * direction[axis],
        );
    }
    close(report.external_work_j, 1.08);
    close(report.kinetic_energy_change_j, 0.66);
    close(report.deposited_heat_j, 0.42);
    close(next.pair().contact().accumulated_distance_m(), 0.21);
    close(report.total_energy_residual_j, 0.0);
}

#[test]
fn planar_stiction_uses_the_vector_load_and_distinct_kinetic_capacity() {
    let (a, ca) = common::material(1.0, 2.0, 0.0);
    let (b, cb) = common::material(1.0, 3.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let point = context();
    let thermal = inventory(&curves, [1.0, 1.0]);
    let card = interface_with_static(
        &cards,
        0.6,
        "speed",
        SPEED,
        [(0.0, 0.2), (4.0, 0.2)],
        2000.0,
    );
    let contact = MaterialSlidingContact::new(
        &thermal,
        binding(
            &cards,
            &curves,
            &card,
            &point,
            HeatPartition::new(0.5, 0.5, 0.0).unwrap(),
        ),
        10.0,
        0.0,
        InputAuthority::SyntheticFixture,
    )
    .unwrap();
    let pair = SlidingPair::new(&thermal, contact, [[0.0; 2]; 2]).unwrap();
    let held = pair
        .step_planar_candidate(
            &thermal,
            0.1,
            MotionDrive {
                normal_load_n: 10.0,
                body_forces_n: [[3.0, 4.0], [-3.0, -4.0]],
            },
            planar_controls(0.1),
            &[0.0; 2],
            other,
            || false,
        )
        .unwrap();
    assert_eq!(held.pair().velocities_mps(), [[0.0; 2]; 2]);
    assert_eq!(held.report().deposited_heat_j, 0.0);
    assert_eq!(
        held.report().segments[0].sliding.friction.traction_n,
        [3.0, 4.0]
    );
    let next = pair
        .step_planar_candidate(
            &thermal,
            0.1,
            MotionDrive {
                normal_load_n: 10.0,
                body_forces_n: [[4.0, 5.0], [-4.0, -5.0]],
            },
            planar_controls(0.1),
            &[0.0; 2],
            other,
            || false,
        )
        .unwrap();
    let force = 41.0_f64.sqrt();
    let acceleration = 2.0 * (force - 2.0);
    for axis in 0..2 {
        let direction = [4.0, 5.0][axis] / force;
        close(
            next.pair().velocities_mps()[0][axis],
            0.05 * acceleration * direction,
        );
        close(
            next.pair().velocities_mps()[1][axis],
            -0.05 * acceleration * direction,
        );
        close(
            next.report().segments[0].sliding.friction.traction_n[axis],
            2.0 * direction,
        );
    }
    close(next.report().deposited_heat_j, acceleration * 0.01);
    close(next.report().total_energy_residual_j, 0.0);
}

#[test]
fn grazing_superkinetic_turn_cannot_be_promoted_to_a_static_branch() {
    let (a, ca) = common::material(1.0, 2.0, 0.0);
    let (b, cb) = common::material(1.0, 3.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let point = context();
    let thermal = inventory(&curves, [1.0, 1.0]);
    let card = interface_with_static(
        &cards,
        0.6,
        "speed",
        SPEED,
        [(0.0, 0.2), (4.0, 0.2)],
        2000.0,
    );
    let contact = MaterialSlidingContact::new(
        &thermal,
        binding(
            &cards,
            &curves,
            &card,
            &point,
            HeatPartition::new(0.5, 0.5, 0.0).unwrap(),
        ),
        10.0,
        1.0_f64.hypot(1e-4),
        InputAuthority::SyntheticFixture,
    )
    .unwrap();
    let pair = SlidingPair::new(&thermal, contact, [[0.5, 0.5e-4], [-0.5, -0.5e-4]]).unwrap();
    let mut controls = planar_controls(0.05);
    // The fictitious arrest's transverse impulse is smaller than this gate.
    // It still cannot stick: 2 N kinetic < 3 N drive < 6 N static, with a
    // nonzero transverse velocity. The continuous path turns, then accelerates.
    controls.impulse_tolerance_n_s = 1e-3;
    let next = pair
        .step_planar_candidate(
            &thermal,
            0.2,
            MotionDrive {
                normal_load_n: 10.0,
                body_forces_n: [[-3.0, 0.0], [3.0, 0.0]],
            },
            controls,
            &[0.0; 2],
            other,
            || false,
        )
        .unwrap();
    let report = next.report();
    assert!(report.segments.iter().all(|segment| !segment.reached_stop));
    assert!(report.halvings > 0);
    let velocity = next.pair().velocities_mps();
    assert!(velocity[0][0] - velocity[1][0] < -0.19);
    assert!(velocity[0][1] - velocity[1][1] >= -1e-14);
    assert!(report.deposited_heat_j > 0.1);
    close(report.total_energy_residual_j, 0.0);
}

#[test]
fn curved_arrest_refines_to_independent_planar_coulomb_path_and_heat() {
    let (a, ca) = common::material(1.0, 2.0, 0.0);
    let (b, cb) = common::material(1.0, 3.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let point = context();
    let thermal = inventory(&curves, [1.0, 1.0]);
    let card = interface_with_static(
        &cards,
        0.4,
        "speed",
        SPEED,
        [(0.0, 0.2), (4.0, 0.2)],
        2000.0,
    );
    let contact = MaterialSlidingContact::new(
        &thermal,
        binding(
            &cards,
            &curves,
            &card,
            &point,
            HeatPartition::new(0.5, 0.5, 0.0).unwrap(),
        ),
        10.0,
        1.0,
        InputAuthority::SyntheticFixture,
    )
    .unwrap();
    let pair = SlidingPair::new(&thermal, contact, [[0.0, 0.5], [0.0, -0.5]]).unwrap();
    // Independent integration of dr/dt=A*cos(theta)-c, dtheta/dt=-A*sin(theta)/r.
    // Here A=2, c=4, r0=1, theta0=pi/2; the finite stopping integrals are
    // t*=c/(c²-A²), x*=1/[8(c-A)]-1/[8(c+A)],
    // y*=1/[2(2c-A)]+1/[2(2c+A)], s*=1/[8(c-A)]+1/(4c)+1/[8(c+A)].
    let exact = [1.0 / 3.0, 1.0 / 24.0, 2.0 / 15.0, 7.0 / 48.0, 7.0 / 24.0];
    let mut errors = Vec::new();
    for h in [0.04, 0.02, 0.01] {
        let next = pair
            .step_planar_candidate(
                &thermal,
                0.5,
                MotionDrive {
                    normal_load_n: 10.0,
                    body_forces_n: [[1.0, 0.0], [-1.0, 0.0]],
                },
                planar_controls(h),
                &[0.0; 2],
                other,
                || false,
            )
            .unwrap();
        let report = next.report();
        assert_eq!(next.pair().velocities_mps(), [[0.0; 2]; 2]);
        let stop = report
            .segments
            .iter()
            .find(|segment| segment.reached_stop)
            .unwrap();
        let actual = [
            stop.end_time_s,
            report.body_displacements_m[0][0] - report.body_displacements_m[1][0],
            report.body_displacements_m[0][1] - report.body_displacements_m[1][1],
            next.pair().contact().accumulated_distance_m(),
            report.deposited_heat_j,
        ];
        let error = actual
            .into_iter()
            .zip(exact)
            .map(|(a, e)| (a - e).abs())
            .fold(0.0_f64, f64::max);
        errors.push(error);
        close(report.total_energy_residual_j, 0.0);
        close(report.deposited_heat_j, 0.25 + report.external_work_j);
        for residual in report.body_momentum_residual_n_s {
            assert!(residual[0].hypot(residual[1]) <= planar_controls(h).impulse_tolerance_n_s);
        }
        assert!(report.segments.iter().any(|s| {
            s.relative_velocity_before_mps[0] > 0.0 && s.relative_velocity_before_mps[1] > 0.0
        }));
    }
    assert!(errors[1] < 0.4 * errors[0], "{errors:?}");
    assert!(errors[2] < 0.4 * errors[1], "{errors:?}");
    assert!(errors[2] < 1e-4, "{errors:?}");
}

fn turning_reference() -> [f64; 5] {
    // Independent continuous ODE integration: mr=.5, F=(0,3), friction=2 N.
    // Positions and arc length are integrated, not reconstructed from energy.
    fn rhs(x: [f64; 5]) -> [f64; 5] {
        let speed = x[0].hypot(x[1]);
        [
            -4.0 * x[0] / speed,
            6.0 - 4.0 * x[1] / speed,
            x[0],
            x[1],
            speed,
        ]
    }
    let h = 0.2 / 32_768.0;
    let mut x = [1.0, 0.0, 0.0, 0.0, 0.0];
    for _ in 0..32_768 {
        let a = rhs(x);
        let b = rhs(std::array::from_fn(|i| x[i] + 0.5 * h * a[i]));
        let c = rhs(std::array::from_fn(|i| x[i] + 0.5 * h * b[i]));
        let d = rhs(std::array::from_fn(|i| x[i] + h * c[i]));
        for i in 0..5 {
            x[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    x
}

#[test]
fn driven_turning_matches_independent_time_refinement_and_keeps_tangent_frame_objectivity() {
    let (a, ca) = common::material(1.0, 2.0, 0.0);
    let (b, cb) = common::material(1.0, 3.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let point = context();
    let thermal = inventory(&curves, [1.0, 1.0]);
    let card = interface_with_static(
        &cards,
        0.6,
        "speed",
        SPEED,
        [(0.0, 0.2), (4.0, 0.2)],
        2000.0,
    );
    let contact = MaterialSlidingContact::new(
        &thermal,
        binding(
            &cards,
            &curves,
            &card,
            &point,
            HeatPartition::new(0.5, 0.5, 0.0).unwrap(),
        ),
        10.0,
        1.0,
        InputAuthority::SyntheticFixture,
    )
    .unwrap();
    let expected = turning_reference();
    let mut errors = Vec::new();
    let rotate = |v: [f64; 2]| [0.6 * v[0] - 0.8 * v[1], 0.8 * v[0] + 0.6 * v[1]];
    for h in [0.04, 0.02, 0.01] {
        let mut reference = None;
        for rotated in [false, true] {
            let transform = |v| if rotated { rotate(v) } else { v };
            let pair = SlidingPair::new(
                &thermal,
                contact.clone(),
                [transform([0.5, 0.0]), transform([-0.5, 0.0])],
            )
            .unwrap();
            let next = pair
                .step_planar_candidate(
                    &thermal,
                    0.2,
                    MotionDrive {
                        normal_load_n: 10.0,
                        body_forces_n: [transform([0.0, 3.0]), transform([0.0, -3.0])],
                    },
                    planar_controls(h),
                    &[0.0; 2],
                    other,
                    || false,
                )
                .unwrap();
            let report = next.report();
            // The target is reached by the real final segment, including when
            // repeated 0.02 s additions would otherwise leave one clock ULP.
            assert_eq!(report.segments.len(), (0.2 / h).round() as usize);
            assert_eq!(report.segments.last().unwrap().end_time_s, 0.2);
            let velocities = next.pair().velocities_mps();
            let actual = [
                velocities[0][0] - velocities[1][0],
                velocities[0][1] - velocities[1][1],
                report.body_displacements_m[0][0] - report.body_displacements_m[1][0],
                report.body_displacements_m[0][1] - report.body_displacements_m[1][1],
                next.pair().contact().accumulated_distance_m(),
            ];
            if rotated {
                let original: [f64; 5] = reference.unwrap();
                let velocity = rotate([original[0], original[1]]);
                let displacement = rotate([original[2], original[3]]);
                for (actual, expected) in actual.into_iter().zip([
                    velocity[0],
                    velocity[1],
                    displacement[0],
                    displacement[1],
                    original[4],
                ]) {
                    close(actual, expected);
                }
            } else {
                errors.push(
                    actual
                        .into_iter()
                        .zip(expected)
                        .map(|(a, e)| (a - e).abs())
                        .fold(0.0_f64, f64::max),
                );
                reference = Some(actual);
            }
            close(report.deposited_heat_j, 2.0 * actual[4]);
            close(report.total_energy_residual_j, 0.0);
            assert!(report.segments.iter().all(|s| !s.reached_stop));
        }
    }
    assert!(errors[1] < 0.3 * errors[0], "{errors:?}");
    assert!(errors[2] < 0.3 * errors[1], "{errors:?}");
    assert!(errors[2] < 3e-4, "{errors:?}");
}

#[test]
fn material_temperature_feedback_changes_the_same_planar_motion_and_heat() {
    let direction = [0.6, 0.8];
    let mut speeds = Vec::new();
    for capacity_a in [0.1, 1.0] {
        let (a, ca) = common::material(1.0, capacity_a, 0.0);
        let (b, cb) = common::material(1.0, 3.0, 0.0);
        let cards = [a, b];
        let curves = [ca, cb];
        let point = context();
        let thermal = inventory(&curves, [1.0, 1.0]);
        let card = interface_with_static(
            &cards,
            0.6,
            "Ta",
            Temperature::DIMS,
            [(299.0, 0.1), (310.0, 1.2)],
            2000.0,
        );
        let contact = MaterialSlidingContact::new(
            &thermal,
            binding(
                &cards,
                &curves,
                &card,
                &point,
                HeatPartition::new(0.5, 0.5, 0.0).unwrap(),
            ),
            10.0,
            1.0,
            InputAuthority::SyntheticFixture,
        )
        .unwrap();
        let pair = SlidingPair::new(&thermal, contact, [direction, [0.0; 2]]).unwrap();
        let next = pair
            .step_planar_candidate(
                &thermal,
                0.1,
                MotionDrive {
                    normal_load_n: 10.0,
                    body_forces_n: [[0.0; 2]; 2],
                },
                planar_controls(0.1),
                &[0.0; 2],
                other,
                || false,
            )
            .unwrap();
        // Independently eliminated source equation: mu=.2+.05*mu*(1-mu)/Ca.
        let k = 0.05 / capacity_a;
        let mu = 0.4 / ((1.0 - k) + ((1.0 - k) * (1.0 - k) + 0.8 * k).sqrt());
        let heat = mu * (1.0 - mu);
        let v = next.pair().velocities_mps();
        let speed = (v[0][0] - v[1][0]).hypot(v[0][1] - v[1][1]);
        close(speed, 1.0 - 2.0 * mu);
        close(next.report().deposited_heat_j, heat);
        close(next.report().kinetic_energy_change_j, -heat);
        close(
            next.thermal().temperature_k(CELLS[0]).unwrap(),
            300.0 + heat / (2.0 * capacity_a),
        );
        close(
            next.report().segments[0].sliding.material.speed_mps,
            1.0 - mu,
        );
        speeds.push(speed);
    }
    assert!(speeds[1] > speeds[0] + 0.1);
}

#[test]
fn discarded_arrest_source_errors_cancellation_and_budgets_leave_every_owner_unchanged() {
    let (a, ca) = common::material(1.0, 2.0, 0.0);
    let (b, cb) = common::material(1.0, 3.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let point = context();
    let thermal = inventory(&curves, [1.0, 1.0]);
    let card = interface_with_static(
        &cards,
        0.4,
        "speed",
        SPEED,
        [(0.0, 0.2), (4.0, 0.2)],
        2000.0,
    );
    let contact = MaterialSlidingContact::new(
        &thermal,
        binding(
            &cards,
            &curves,
            &card,
            &point,
            HeatPartition::new(0.5, 0.5, 0.0).unwrap(),
        ),
        10.0,
        1.0,
        InputAuthority::SyntheticFixture,
    )
    .unwrap();
    let pair = SlidingPair::new(&thermal, contact, [[0.6, 0.8], [0.0; 2]]).unwrap();
    let drive = MotionDrive {
        normal_load_n: 10.0,
        body_forces_n: [[0.0; 2]; 2],
    };
    let old_cells = thermal.cells().clone();
    let mut limited = planar_controls(1.0);
    limited.maximum_segments = 1;
    assert!(matches!(
        pair.step_planar_candidate(&thermal, 1.0, drive, limited, &[0.0; 2], other, || false),
        Err(SlidingError::MotionBudget)
    ));
    // A coarse oblique midpoint trial would reverse the velocity without
    // resolving its turning path. Exhaustion must refuse, never publish it.
    let turning =
        SlidingPair::new(&thermal, pair.contact().clone(), [[0.0, 0.5], [0.0, -0.5]]).unwrap();
    let mut limited = planar_controls(0.5);
    limited.maximum_halvings = 0;
    assert!(matches!(
        turning.step_planar_candidate(
            &thermal,
            0.5,
            MotionDrive {
                normal_load_n: 10.0,
                body_forces_n: [[1.0, 0.0], [-1.0, 0.0]],
            },
            limited,
            &[0.0; 2],
            other,
            || false,
        ),
        Err(SlidingError::MotionBudget)
    ));
    assert_eq!(turning.velocities_mps(), [[0.0, 0.5], [0.0, -0.5]]);
    let mut limited = planar_controls(1.0);
    limited.maximum_motion_trials = 1;
    assert!(matches!(
        pair.step_planar_candidate(&thermal, 1.0, drive, limited, &[0.0; 2], other, || false),
        Err(SlidingError::MotionBudget)
    ));
    assert!(matches!(
        pair.step_planar_candidate(
            &thermal,
            1.0,
            drive,
            planar_controls(1.0),
            &[0.0; 2],
            |_state| Err::<StepCoefficients<()>, _>("source unavailable"),
            || false
        ),
        Err(SlidingError::Source("source unavailable"))
    ));
    let polls = Cell::new(0usize);
    let good = pair
        .step_planar_candidate(
            &thermal,
            1.0,
            drive,
            planar_controls(1.0),
            &[0.0; 2],
            other,
            || {
                polls.set(polls.get() + 1);
                false
            },
        )
        .unwrap();
    let last = polls.get();
    polls.set(0);
    assert!(matches!(
        pair.step_planar_candidate(
            &thermal,
            1.0,
            drive,
            planar_controls(1.0),
            &[0.0; 2],
            other,
            || {
                polls.set(polls.get() + 1);
                polls.get() == last
            }
        ),
        Err(SlidingError::Cancelled)
    ));
    assert_eq!(thermal.cells(), &old_cells);
    assert_eq!(pair.contact().time_s(), 0.0);
    assert_eq!(pair.contact().state().slip_m(), [0.0; 2]);
    assert_eq!(pair.contact().work().dissipated_work_j(), 0.0);
    assert_eq!(pair.velocities_mps(), [[0.6, 0.8], [0.0; 2]]);
    let retry = pair
        .step_planar_candidate(
            &thermal,
            1.0,
            drive,
            planar_controls(1.0),
            &[0.0; 2],
            other,
            || false,
        )
        .unwrap();
    assert_eq!(retry.pair().velocities_mps(), good.pair().velocities_mps());
    assert_eq!(
        retry.pair().contact().state(),
        good.pair().contact().state()
    );
    assert_eq!(retry.thermal().cells(), good.thermal().cells());
}
