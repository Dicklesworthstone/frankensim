//! Pressure feedback and work transfer between two original real-gas fields.

use super::*;
use fs_gas::dynamics::Primitive;

const AREA: f64 = 0.07;

fn chamber(
    gas: RealGas,
    start: f64,
    length: f64,
    count: usize,
    temperature: f64,
    physical: bool,
) -> RealGasEuler1d {
    let positions: Vec<_> = (0..=count)
        .map(|i| start + length * (i as f64 / count as f64))
        .collect();
    let volumes = positions.windows(2).map(|x| AREA * (x[1] - x[0])).collect();
    let mesh = Mesh1d::new(positions, vec![AREA; count + 1], volumes).unwrap();
    let boundary = if physical {
        RealGasBoundary::PhysicalReflecting(RealGasWallControl {
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
        })
    } else {
        RealGasBoundary::Reflecting
    };
    let q = state(&gas, 10.0, 0.0, temperature);
    RealGasEuler1d::new_on_mesh(
        gas,
        mesh,
        vec![q; count],
        RealGasBoundaries::Ends {
            left: boundary,
            right: boundary,
        },
        || false,
    )
    .unwrap()
}

fn pair(temperatures: [f64; 2], physical: bool) -> (RealGasEuler1d, RealGasEuler1d) {
    (
        chamber(gas(), 0.0, 1.0, 3, temperatures[0], physical),
        chamber(gas(), 1.0, 0.7, 5, temperatures[1], physical),
    )
}

#[derive(Debug, PartialEq)]
struct Snapshot {
    states: Vec<Conserved>,
    primitives: Vec<Primitive>,
    positions: Vec<f64>,
    volumes: Vec<f64>,
    time: Time,
}
fn snapshot(flow: &RealGasEuler1d) -> Snapshot {
    Snapshot {
        states: flow.states().to_vec(),
        primitives: flow.primitives().to_vec(),
        positions: flow.mesh().positions().to_vec(),
        volumes: flow.mesh().volumes().to_vec(),
        time: flow.time(),
    }
}

#[test]
fn equal_pressure_retains_both_original_fields_on_unequal_grids() {
    for physical in [false, true] {
        let (mut left, mut right) = pair([450.0; 2], physical);
        let before = [snapshot(&left), snapshot(&right)];
        let pressure = left.primitives()[0].pressure;
        let mut piston = Piston::new(5.0, 0.0).unwrap();
        let report = left
            .step_with_shared_piston(
                &mut right,
                &mut piston,
                Time::new(2e-6),
                muscl(),
                0.0,
                mechanical_controls(),
                || false,
            )
            .unwrap();
        assert_eq!(piston.velocity_m_s(), 0.0);
        assert_eq!(report.displacement_m, 0.0);
        assert_eq!(report.exterior_pressure_work_j, 0.0);
        assert_eq!(report.total_energy_residual_j, 0.0);
        assert_eq!(report.gas_trials, 1);
        assert_eq!(left.time(), right.time());
        for (flow, original) in [(&left, &before[0]), (&right, &before[1])] {
            assert_eq!(flow.states(), original.states);
            assert_eq!(flow.primitives(), original.primitives);
            assert_eq!(flow.mesh().positions(), original.positions);
            assert_eq!(flow.mesh().volumes(), original.volumes);
        }
        near(
            report.gas.left.right_boundary_change.momentum,
            -pressure * AREA * report.gas.left.elapsed.value(),
            1e-12,
        );
        near(
            report.gas.right.left_boundary_change.momentum,
            pressure * AREA * report.gas.right.elapsed.value(),
            1e-12,
        );
    }
}

#[test]
fn pressure_difference_drives_both_directions_and_transfers_work_between_real_gases() {
    for temperatures in [[550.0, 400.0], [400.0, 550.0]] {
        let (mut left, mut right) = pair(temperatures, false);
        let before = [
            left.totals(|| false).unwrap(),
            right.totals(|| false).unwrap(),
        ];
        let mut piston = Piston::new(5.0, 0.0).unwrap();
        let acceleration = AREA * (left.primitives()[0].pressure - right.primitives()[0].pressure)
            / piston.mass_kg();
        let first = left
            .step_with_shared_piston(
                &mut right,
                &mut piston,
                Time::new(1e-7),
                muscl(),
                0.0,
                mechanical_controls(),
                || false,
            )
            .unwrap();
        near(
            piston.velocity_m_s() / first.gas.left.elapsed.value(),
            acceleration,
            3e-5 * acceleration.abs(),
        );
        let mut fixed_impulse = first.fixed_wall_impulse_n_s;
        for _ in 0..4 {
            let report = left
                .step_with_shared_piston(
                    &mut right,
                    &mut piston,
                    Time::new(2e-6),
                    muscl(),
                    0.0,
                    mechanical_controls(),
                    || false,
                )
                .unwrap();
            for side in [report.gas.left, report.gas.right] {
                check_balance(side);
                assert_eq!(side.boundary_change.mass, 0.0);
                assert_eq!(side.wall_heat, 0.0);
            }
            assert!(report.displacement_m * acceleration > 0.0);
            assert!(report.gas.left.right_boundary_change.energy * acceleration < 0.0);
            assert!(report.gas.right.left_boundary_change.energy * acceleration > 0.0);
            assert_eq!(report.gas.left.left_boundary_change.energy, 0.0);
            assert_eq!(report.gas.right.right_boundary_change.energy, 0.0);
            near(
                report.gas.left.right_boundary_change.energy
                    + report.gas.right.left_boundary_change.energy,
                -report.midpoint_velocity_m_s * report.gas_impulse_n_s,
                2e-12,
            );
            assert_eq!(left.mesh().positions()[3], right.mesh().positions()[0]);
            near(
                left.mesh().volumes().iter().sum::<f64>()
                    + right.mesh().volumes().iter().sum::<f64>(),
                1.7 * AREA,
                1e-14,
            );
            fixed_impulse += report.fixed_wall_impulse_n_s;
        }
        let after = [
            left.totals(|| false).unwrap(),
            right.totals(|| false).unwrap(),
        ];
        for i in 0..2 {
            near(after[i].mass, before[i].mass, 1e-13);
        }
        near(
            after[0].energy - before[0].energy + after[1].energy - before[1].energy
                + 0.5 * piston.mass_kg() * piston.velocity_m_s().powi(2),
            0.0,
            2e-6,
        );
        near(
            after[0].momentum - before[0].momentum + after[1].momentum - before[1].momentum
                + piston.mass_kg() * piston.velocity_m_s()
                + fixed_impulse,
            0.0,
            2e-9,
        );
        assert_eq!(left.mesh().positions()[0], 0.0);
        assert_eq!(right.mesh().positions()[5], 1.7);
    }
}

#[test]
fn the_shorter_chamber_sets_one_duration_and_both_candidates_are_native_ale_fields() {
    for physical in [false, true] {
        let make = || {
            (
                chamber(gas(), 0.0, 1.0, 3, 500.0, physical),
                chamber(gas(), 1.0, 0.01, 4, 400.0, physical),
            )
        };
        let (mut left, mut right) = make();
        let (mut prescribed_left, mut prescribed_right) = make();
        let mut piston = Piston::new(5.0, 0.0).unwrap();
        let report = left
            .step_with_shared_piston(
                &mut right,
                &mut piston,
                Time::new(1.0),
                muscl(),
                0.0,
                mechanical_controls(),
                || false,
            )
            .unwrap();
        assert!(report.gas.left.elapsed.value() < 1e-5);
        assert_eq!(report.gas.left.elapsed, report.gas.right.elapsed);
        assert_eq!(left.time(), right.time());
        let left_report = prescribed_left
            .step_with_moving_walls(
                report.gas.left.elapsed,
                muscl(),
                [0.0, report.midpoint_velocity_m_s],
                || false,
            )
            .unwrap();
        let right_report = prescribed_right
            .step_with_moving_walls(
                report.gas.right.elapsed,
                muscl(),
                [report.midpoint_velocity_m_s, 0.0],
                || false,
            )
            .unwrap();
        assert_eq!(report.gas.left, left_report);
        assert_eq!(report.gas.right, right_report);
        assert_eq!(snapshot(&left), snapshot(&prescribed_left));
        assert_eq!(snapshot(&right), snapshot(&prescribed_right));
    }
}

#[test]
fn independent_source_datums_spring_storage_and_damping_keep_the_two_gas_work_balance() {
    let (mut left, mut right) = pair([460.0, 440.0], false);
    let mut shifted_left = chamber(
        bind(
            common::model([1e22, -3e21], [2e18, -7e17], [200.0, 2000.0]),
            domain(),
        ),
        0.0,
        1.0,
        3,
        460.0,
        false,
    );
    let mut shifted_right = chamber(
        bind(
            common::model([-7e21, 5e21], [-3e18, 2e18], [200.0, 2000.0]),
            domain(),
        ),
        1.0,
        0.7,
        5,
        440.0,
        false,
    );
    let mut piston = Piston::new(3.0, 2.0).unwrap();
    let mut shifted_piston = piston;
    let spring = PistonSuspension::new(8000.0, 80.0, 0.95).unwrap();
    let initial = [
        left.totals(|| false).unwrap(),
        right.totals(|| false).unwrap(),
    ];
    let before_spring = spring.stored_energy_j(1.0).unwrap();
    let (mut applied, mut damping) = (0.0, 0.0);
    for _ in 0..4 {
        let report = left
            .step_with_shared_suspended_piston(
                &mut right,
                &mut piston,
                Time::new(2e-6),
                muscl(),
                700.0,
                spring,
                mechanical_controls(),
                || false,
            )
            .unwrap();
        let shifted = shifted_left
            .step_with_shared_suspended_piston(
                &mut shifted_right,
                &mut shifted_piston,
                Time::new(2e-6),
                muscl(),
                700.0,
                spring,
                mechanical_controls(),
                || false,
            )
            .unwrap();
        assert_eq!(report, shifted);
        assert_eq!(snapshot(&left), snapshot(&shifted_left));
        assert_eq!(snapshot(&right), snapshot(&shifted_right));
        assert_eq!(piston, shifted_piston);
        assert!(report.damping_energy_j > 0.0);
        near(
            report.applied_force_work_j,
            700.0 * report.piston.displacement_m,
            1e-15,
        );
        near(report.total_energy_residual_j, 0.0, 3e-7);
        applied += report.applied_force_work_j;
        damping += report.damping_energy_j;
    }
    let after = [
        left.totals(|| false).unwrap(),
        right.totals(|| false).unwrap(),
    ];
    near(
        after[0].energy - initial[0].energy + after[1].energy - initial[1].energy
            + 1.5 * (piston.velocity_m_s().powi(2) - 4.0)
            + spring.stored_energy_j(left.mesh().positions()[3]).unwrap()
            - before_spring
            + damping
            - applied,
        0.0,
        2e-6,
    );
}

#[test]
fn a_refused_root_and_final_cancellation_preserve_both_complete_owners() {
    let (mut reference_left, mut reference_right) = pair([550.0, 400.0], false);
    let mut expected_piston = Piston::new(5.0, 0.0).unwrap();
    let mut polls = 0;
    let expected = reference_left
        .step_with_shared_piston(
            &mut reference_right,
            &mut expected_piston,
            Time::new(2e-6),
            muscl(),
            0.0,
            mechanical_controls(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let (mut left, mut right) = pair([550.0, 400.0], false);
    let before = [snapshot(&left), snapshot(&right)];
    let mut piston = Piston::new(5.0, 0.0).unwrap();
    let original_piston = piston;
    let error = left
        .step_with_shared_piston(
            &mut right,
            &mut piston,
            Time::new(2e-6),
            StepControl {
                maximum_halvings: 0,
                ..muscl()
            },
            0.0,
            PistonControl {
                midpoint_velocity_bracket_m_s: [-1e-6, 1e-6],
                ..mechanical_controls()
            },
            || false,
        )
        .unwrap_err();
    assert!(matches!(error, PistonError::StepRetries { last, .. }
        if matches!(*last, PistonError::OutsideBracket { .. })));
    assert_eq!([snapshot(&left), snapshot(&right)], before);
    assert_eq!(piston, original_piston);
    let mut observed = 0;
    assert_eq!(
        left.step_with_shared_piston(
            &mut right,
            &mut piston,
            Time::new(2e-6),
            muscl(),
            0.0,
            mechanical_controls(),
            || {
                observed += 1;
                observed == polls
            }
        ),
        Err(PistonError::Cancelled)
    );
    assert_eq!([snapshot(&left), snapshot(&right)], before);
    assert_eq!(piston, original_piston);
    let actual = left
        .step_with_shared_piston(
            &mut right,
            &mut piston,
            Time::new(2e-6),
            muscl(),
            0.0,
            mechanical_controls(),
            || false,
        )
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(snapshot(&left), snapshot(&reference_left));
    assert_eq!(snapshot(&right), snapshot(&reference_right));
    assert_eq!(piston, expected_piston);
}
