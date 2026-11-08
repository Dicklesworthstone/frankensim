//! Work-conjugate finite inertia on the original nonideal pressure/energy owner.

use super::*;
use fs_gas::dynamics::nonideal::RealGasWallControl;
use fs_gas::dynamics::piston::{Piston, PistonControl, PistonError, PistonLoad, PistonSuspension};
use fs_gas::material::nonideal_shock::RealGasShockControl;

#[path = "piston/forcing.rs"]
mod forcing;

fn mechanical_controls() -> PistonControl {
    PistonControl {
        midpoint_velocity_bracket_m_s: [-20.0, 20.0],
        momentum_tolerance_n_s: 2e-10,
        energy_tolerance_j: 3e-7,
        maximum_iterations: 64,
    }
}

fn field(gas: RealGas, physical: bool) -> RealGasEuler1d {
    let q = state(&gas, 10.0, 0.0, 450.0);
    let end = if physical {
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
    RealGasEuler1d::new(
        gas,
        Length::new(0.125),
        0.07,
        vec![q; 8],
        RealGasBoundaries::Ends {
            left: end,
            right: end,
        },
        || false,
    )
    .unwrap()
}

fn muscl() -> StepControl {
    StepControl {
        reconstruction: Reconstruction::Muscl,
        ..controls()
    }
}

#[test]
fn equal_original_pressure_holds_a_finite_piston_at_rest_without_work() {
    for physical in [false, true] {
        let mut flow = field(gas(), physical);
        let states = flow.states().to_vec();
        let positions = flow.mesh().positions().to_vec();
        let volumes = flow.mesh().volumes().to_vec();
        let primitive = flow.primitives().to_vec();
        let mut piston = Piston::new(5.0, 0.0).unwrap();
        let report = flow
            .step_with_piston(
                &mut piston,
                Time::new(2e-6),
                muscl(),
                PistonLoad {
                    exterior_pressure_pa: primitive[0].pressure,
                    external_force_n: 0.0,
                },
                mechanical_controls(),
                || false,
            )
            .unwrap();
        assert_eq!(flow.states(), states);
        assert_eq!(flow.primitives(), primitive);
        assert_eq!(flow.mesh().positions(), positions);
        assert_eq!(flow.mesh().volumes(), volumes);
        assert_eq!(piston.velocity_m_s(), 0.0);
        assert_eq!(report.displacement_m, 0.0);
        assert_eq!(report.gas.boundary_change.energy, 0.0);
        assert_eq!(report.total_energy_residual_j, 0.0);
        assert_eq!(report.gas_trials, 1);
    }
}

#[test]
fn dense_gas_expansion_and_compression_exchange_actual_piston_kinetic_energy() {
    for physical in [false, true] {
        for initial_velocity in [0.0, -12.0] {
            let mut flow = field(gas(), physical);
            let before = flow.totals(|| false).unwrap();
            let mut piston = Piston::new(5.0, initial_velocity).unwrap();
            let initial_kinetic = 0.5 * piston.mass_kg() * initial_velocity * initial_velocity;
            let mut pressure_work = 0.0;
            let mut support_impulse = 0.0;
            for _ in 0..6 {
                let report = flow
                    .step_with_piston(
                        &mut piston,
                        Time::new(1e-6),
                        muscl(),
                        PistonLoad::default(),
                        mechanical_controls(),
                        || false,
                    )
                    .unwrap();
                check_balance(report.gas);
                assert_eq!(report.gas.left_boundary_change.energy, 0.0);
                assert_eq!(report.gas.boundary_change.mass, 0.0);
                assert_eq!(report.gas.wall_heat, 0.0);
                assert_eq!(report.gas.wall_work, 0.0);
                near(report.total_energy_residual_j, 0.0, 3e-7);
                near(report.combined_momentum_residual_n_s, 0.0, 2e-10);
                pressure_work += report.gas.boundary_change.energy;
                support_impulse += report.fixed_wall_impulse_n_s;
            }
            let after = flow.totals(|| false).unwrap();
            let final_kinetic = 0.5 * piston.mass_kg() * piston.velocity_m_s().powi(2);
            near(after.mass, before.mass, 2e-14);
            near(
                (after.energy - before.energy) + (final_kinetic - initial_kinetic),
                0.0,
                2e-6,
            );
            near(after.energy - before.energy, pressure_work, 1e-8);
            near(
                after.momentum - before.momentum
                    + piston.mass_kg() * (piston.velocity_m_s() - initial_velocity)
                    + support_impulse,
                0.0,
                2e-9,
            );
            if initial_velocity == 0.0 {
                assert!(piston.velocity_m_s() > 0.0);
                assert!(flow.mesh().positions()[8] > 1.0);
                assert!(after.energy < before.energy);
            } else {
                assert!(piston.velocity_m_s() > initial_velocity && piston.velocity_m_s() < 0.0);
                assert!(flow.mesh().positions()[8] < 1.0);
                assert!(after.energy > before.energy);
            }
        }
    }
}

#[test]
fn coupled_candidate_is_the_native_real_gas_ale_field_at_its_solved_speed() {
    for physical in [false, true] {
        let source = gas();
        let mut coupled = field(source.clone(), physical);
        let mut prescribed = field(source, physical);
        let mut piston = Piston::new(3.0, -7.0).unwrap();
        let report = coupled
            .step_with_piston(
                &mut piston,
                Time::new(1e-6),
                muscl(),
                PistonLoad::default(),
                mechanical_controls(),
                || false,
            )
            .unwrap();
        let original = prescribed
            .step_with_moving_walls(
                report.gas.elapsed,
                muscl(),
                [0.0, report.midpoint_velocity_m_s],
                || false,
            )
            .unwrap();
        assert_eq!(report.gas.elapsed, original.elapsed);
        assert_eq!(coupled.states(), prescribed.states());
        assert_eq!(coupled.primitives(), prescribed.primitives());
        assert_eq!(coupled.mesh().positions(), prescribed.mesh().positions());
        assert_eq!(coupled.mesh().volumes(), prescribed.mesh().volumes());
        assert_eq!(report.gas.boundary_change, original.boundary_change);
        assert!(piston.velocity_m_s() > -7.0);
    }
}

#[test]
fn external_drive_spring_and_damper_close_the_combined_real_gas_work_balance() {
    let mut flow = field(gas(), false);
    let mut piston = Piston::new(2.0, 3.0).unwrap();
    let spring = PistonSuspension::new(8e3, 80.0, 0.95).unwrap();
    let before = flow.totals(|| false).unwrap();
    let kinetic_before = 0.5 * piston.mass_kg() * piston.velocity_m_s().powi(2);
    let spring_before = spring.stored_energy_j(1.0).unwrap();
    let load = PistonLoad {
        exterior_pressure_pa: flow.primitives()[0].pressure,
        external_force_n: 700.0,
    };
    let (mut input_work, mut exterior_work, mut damping) = (0.0, 0.0, 0.0);
    for _ in 0..8 {
        let report = flow
            .step_with_suspended_piston(
                &mut piston,
                Time::new(2e-6),
                muscl(),
                load,
                spring,
                mechanical_controls(),
                || false,
            )
            .unwrap();
        assert!(report.damping_energy_j > 0.0);
        near(report.total_energy_residual_j, 0.0, 3e-7);
        near(
            report.applied_force_work_j,
            load.external_force_n * report.piston.displacement_m,
            1e-15,
        );
        input_work += report.applied_force_work_j;
        exterior_work += report.piston.exterior_pressure_work_j;
        damping += report.damping_energy_j;
    }
    let after = flow.totals(|| false).unwrap();
    let kinetic_after = 0.5 * piston.mass_kg() * piston.velocity_m_s().powi(2);
    let spring_after = spring.stored_energy_j(flow.mesh().positions()[8]).unwrap();
    near(
        (after.energy - before.energy)
            + (kinetic_after - kinetic_before)
            + (spring_after - spring_before)
            + damping
            + exterior_work
            - input_work,
        0.0,
        3e-6,
    );
}

#[test]
fn source_datum_changes_cannot_change_real_gas_piston_motion_or_work() {
    let shifted = bind(
        common::model([1e22, -3e21], [2e18, -7e17], [200.0, 2000.0]),
        domain(),
    );
    let mut original = field(gas(), false);
    let mut translated = field(shifted, false);
    let mut a = Piston::new(3.0, -8.0).unwrap();
    let mut b = a;
    for _ in 0..4 {
        let x = original
            .step_with_piston(
                &mut a,
                Time::new(2e-6),
                muscl(),
                PistonLoad::default(),
                mechanical_controls(),
                || false,
            )
            .unwrap();
        let y = translated
            .step_with_piston(
                &mut b,
                Time::new(2e-6),
                muscl(),
                PistonLoad::default(),
                mechanical_controls(),
                || false,
            )
            .unwrap();
        assert_eq!(x, y);
        assert_eq!(a, b);
        assert_eq!(original.states(), translated.states());
        assert_eq!(original.primitives(), translated.primitives());
    }
}

#[test]
fn final_poll_cancellation_and_failed_bracket_publish_no_gas_cache_mesh_or_piston() {
    let source = gas();
    let mut reference = field(source.clone(), false);
    let mut probe = Piston::new(4.0, 0.0).unwrap();
    let mut polls = 0;
    let expected = reference
        .step_with_piston(
            &mut probe,
            Time::new(2e-6),
            muscl(),
            PistonLoad::default(),
            mechanical_controls(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(polls > 100);
    let mut flow = field(source, false);
    let states = flow.states().to_vec();
    let primitive = flow.primitives().to_vec();
    let positions = flow.mesh().positions().to_vec();
    let volumes = flow.mesh().volumes().to_vec();
    let time = flow.time();
    let mut piston = Piston::new(4.0, 0.0).unwrap();
    let initial_piston = piston;
    let mut observed = 0;
    let error = flow.step_with_piston(
        &mut piston,
        Time::new(2e-6),
        muscl(),
        PistonLoad::default(),
        mechanical_controls(),
        || {
            observed += 1;
            observed == polls
        },
    );
    assert_eq!(error, Err(PistonError::Cancelled));
    assert_eq!(flow.states(), states);
    assert_eq!(flow.primitives(), primitive);
    assert_eq!(flow.mesh().positions(), positions);
    assert_eq!(flow.mesh().volumes(), volumes);
    assert_eq!(flow.time(), time);
    assert_eq!(piston, initial_piston);
    let error = flow
        .step_with_piston(
            &mut piston,
            Time::new(2e-6),
            StepControl {
                maximum_halvings: 0,
                ..muscl()
            },
            PistonLoad::default(),
            PistonControl {
                midpoint_velocity_bracket_m_s: [-1e-6, 1e-6],
                ..mechanical_controls()
            },
            || false,
        )
        .unwrap_err();
    assert!(matches!(error, PistonError::StepRetries { last, .. }
        if matches!(*last, PistonError::OutsideBracket { .. })));
    assert_eq!(flow.states(), states);
    assert_eq!(flow.primitives(), primitive);
    assert_eq!(flow.mesh().positions(), positions);
    assert_eq!(flow.mesh().volumes(), volumes);
    assert_eq!(flow.time(), time);
    assert_eq!(piston, initial_piston);
    let actual = flow
        .step_with_piston(
            &mut piston,
            Time::new(2e-6),
            muscl(),
            PistonLoad::default(),
            mechanical_controls(),
            || false,
        )
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(flow.states(), reference.states());
    assert_eq!(flow.primitives(), reference.primitives());
    assert_eq!(piston, probe);
}
