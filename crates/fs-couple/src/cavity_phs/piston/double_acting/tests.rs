//! Actual material-owner integration; gas fixture data are synthetic, not fits.
use super::*;
use crate::cavity_phs::piston::tests::{controls as single_controls, drive, gas};

fn controls() -> DoublePistonSolve {
    DoublePistonSolve {
        piston: single_controls(),
        opposing_gas: single_controls().gas,
        volume_tolerance_m3: 1e-18,
    }
}
fn drive_pair(heat_j: [f64; 2]) -> DoublePistonInput {
    DoublePistonInput {
        duration_s: drive(0.0).duration_s,
        heat_j,
        external_force_n: 0.0,
    }
}
fn pair(x: f64, slope: f64, opposing_area: f64, velocity: f64) -> DoubleActingGasPiston {
    DoubleActingGasPiston::new(
        GasPiston::new(gas(0.75, 0.001), 1.0, 0.01, velocity).unwrap(),
        gas(x, slope),
        opposing_area,
    )
    .unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.16e} != {b:.16e}");
}
fn check(report: &DoublePistonStep) {
    assert!(report.momentum_residual_n_s.abs() <= controls().piston.momentum_tolerance_n_s);
    assert!(report.total_energy_residual_j.abs() <= controls().piston.energy_tolerance_j);
    assert!(report.mechanical_energy_residual_j.abs() <= controls().piston.energy_tolerance_j);
    assert!(report.volume_exchange_residual_m3.abs() <= controls().volume_tolerance_m3);
    assert!(report.material_trials <= 2 * (controls().piston.maximum_iterations + 3));
    assert_eq!(report.material_trials % 2, 0);
}

#[test]
fn distinct_materials_at_equal_pressure_remain_bitwise_stationary() {
    let mut piston = pair(0.25, 0.0, 0.01, 0.0);
    let before = piston.clone();
    let report = piston
        .step(drive_pair([0.0, 0.0]), controls(), || false)
        .unwrap();
    assert_eq!(piston, before);
    assert_eq!(report.displacement_m, 0.0);
    assert_eq!(report.total_energy_residual_j, 0.0);
    assert_eq!(report.material_trials, 2);
    assert_ne!(
        report.midpoint[0].state().sound_speed,
        report.midpoint[1].state().sound_speed
    );
    check(&report);
}

#[test]
fn primary_heating_compresses_and_heats_the_finite_opposing_material() {
    let mut piston = pair(0.25, 0.001, 0.01, 0.0);
    let before = piston.clone();
    let report = piston
        .step(drive_pair([20.0, 0.0]), controls(), || false)
        .unwrap();
    assert!(piston.primary().velocity_m_s() > 0.0);
    assert!(piston.opposing().volume_m3() < before.opposing().volume_m3());
    assert!(
        piston.opposing().state().state().temperature
            > before.opposing().state().state().temperature
    );
    assert!(
        piston.opposing().state().state().pressure > before.opposing().state().state().pressure
    );
    assert!(
        piston.opposing().state().state().sound_speed
            > before.opposing().state().state().sound_speed
    );
    assert_ne!(
        piston.opposing().state().state().dynamic_viscosity,
        before.opposing().state().state().dynamic_viscosity
    );
    assert!(report.gas[0].work_by_gas_j > 0.0);
    assert!(report.gas[1].work_by_gas_j < 0.0);
    near(
        report.kinetic_energy_change_j,
        report.gas[0].work_by_gas_j + report.gas[1].work_by_gas_j,
        1e-7,
    );
    near(
        piston.primary().gas().volume_m3() + piston.opposing().volume_m3(),
        0.002,
        1e-18,
    );
    for (old, new) in [
        (before.primary().gas(), piston.primary().gas()),
        (before.opposing(), piston.opposing()),
    ] {
        assert_eq!(new.mass_kg(), old.mass_kg());
        assert_eq!(new.point().axes().get("specimen-lot"), Some(&42.0));
        assert_eq!(new.state().components().len(), 2);
    }
    check(&report);
}

#[test]
fn finite_opposing_gas_is_not_a_fixed_pressure_reservoir() {
    let mut pair = pair(0.75, 0.001, 0.01, 0.0);
    let mut reservoir = pair.primary().clone();
    let report = pair
        .step(drive_pair([20.0, 0.0]), controls(), || false)
        .unwrap();
    reservoir
        .step(drive(20.0), single_controls(), || false)
        .unwrap();
    assert!(pair.primary().velocity_m_s() < reservoir.velocity_m_s());
    assert!(report.midpoint[1].state().pressure > 100000.0);
    check(&report);
}

#[test]
fn swapping_equal_area_chambers_reverses_motion_without_changing_energy() {
    let mut a = pair(0.75, 0.001, 0.01, 0.0);
    let mut b = a.clone();
    let ra = a
        .step(drive_pair([20.0, 0.0]), controls(), || false)
        .unwrap();
    let rb = b
        .step(drive_pair([0.0, 20.0]), controls(), || false)
        .unwrap();
    assert!(rb.displacement_m < 0.0);
    near(ra.displacement_m, -rb.displacement_m, 1e-12);
    near(
        a.primary().velocity_m_s(),
        -b.primary().velocity_m_s(),
        1e-9,
    );
    near(ra.kinetic_energy_change_j, rb.kinetic_energy_change_j, 1e-9);
    near(
        a.primary().gas().state().state().temperature,
        b.opposing().state().state().temperature,
        1e-7,
    );
    check(&ra);
    check(&rb);
}

#[test]
fn unequal_faces_use_one_displacement_and_both_true_pressure_areas() {
    let mut piston = pair(0.25, 0.001, 0.005, 0.1);
    let before = piston.clone();
    let input = DoublePistonInput {
        external_force_n: -500.0,
        ..drive_pair([5.0, 3.0])
    };
    let report = piston.step(input, controls(), || false).unwrap();
    near(
        piston.primary().gas().volume_m3() - before.primary().gas().volume_m3(),
        0.01 * report.displacement_m,
        1e-18,
    );
    near(
        piston.opposing().volume_m3() - before.opposing().volume_m3(),
        -0.005 * report.displacement_m,
        1e-18,
    );
    let force = 0.01 * report.midpoint[0].state().pressure
        - 0.005 * report.midpoint[1].state().pressure
        - 500.0;
    near(
        piston.primary().velocity_m_s() - before.primary().velocity_m_s(),
        input.duration_s * force,
        1e-10,
    );
    near(
        report.displacement_m,
        0.5 * input.duration_s * (report.velocity_m_s[0] + report.velocity_m_s[1]),
        1e-15,
    );
    check(&report);
}

#[test]
fn opposing_composition_data_changes_both_gas_and_mechanical_trajectory() {
    let run = |x, slope| {
        let mut piston = pair(x, slope, 0.01, 0.0);
        piston
            .step(drive_pair([20.0, 0.0]), controls(), || false)
            .unwrap();
        piston
    };
    let a = run(0.75, 0.001);
    let b = run(0.25, 0.001);
    let c = run(0.75, 0.0);
    // One step with 20 J: the composition changes the step-mean acceleration by
    // more than 2e-5 m/s^2 (1e-8 m/s on the authored 0.0005 s step).
    let dv = (a.primary().velocity_m_s() - b.primary().velocity_m_s()).abs();
    assert!(dv > 2e-5 * drive(0.0).duration_s, "{dv}");
    assert_ne!(
        a.primary().gas().state().state().temperature,
        b.primary().gas().state().state().temperature
    );
    assert_ne!(
        a.opposing().state().state().sound_speed,
        b.opposing().state().state().sound_speed
    );
    assert_ne!(a.primary().velocity_m_s(), c.primary().velocity_m_s());
}

#[test]
fn joint_tangent_includes_opposing_discrete_thermal_response() {
    let piston = pair(0.25, 0.001, 0.006, -0.1);
    let input = drive_pair([20.0, 7.0]);
    let volume = 0.00102;
    let delta = 1e-8;
    let center = piston
        .trial(volume, input, controls(), &mut || false)
        .unwrap();
    let plus = piston
        .trial(volume + delta, input, controls(), &mut || false)
        .unwrap();
    let minus = piston
        .trial(volume - delta, input, controls(), &mut || false)
        .unwrap();
    let fd = (plus.residual - minus.residual) / (2.0 * delta);
    assert!((center.tangent / fd - 1.0).abs() < 2e-6);
    assert!(center.tangent > center.value.primary.tangent);
}

#[test]
fn long_trajectory_has_one_kinetic_storage_and_no_external_pressure_supply() {
    let mut piston = pair(0.25, 0.001, 0.01, 0.03);
    let initial = piston.clone();
    let mut supplied = 0.0;
    for i in 0..80 {
        let input = DoublePistonInput {
            external_force_n: if i < 20 { 0.5 } else { 0.0 },
            ..drive_pair([
                if i < 10 { 0.2 } else { 0.0 },
                if i < 5 { 0.1 } else { 0.0 },
            ])
        };
        let report = piston.step(input, controls(), || false).unwrap();
        supplied += input.heat_j.iter().sum::<f64>() + report.external_force_work_j;
        check(&report);
    }
    let gas_change: f64 = [
        (initial.primary().gas(), piston.primary().gas()),
        (initial.opposing(), piston.opposing()),
    ]
    .iter()
    .map(|(old, new)| {
        old.mass_kg()
            * (new.state().specific_internal_energy_j_kg()
                - old.state().specific_internal_energy_j_kg())
    })
    .sum();
    let kinetic = 0.5
        * piston.primary().mass_kg()
        * (piston.primary().velocity_m_s().powi(2) - initial.primary().velocity_m_s().powi(2));
    near(
        gas_change + kinetic,
        supplied,
        80.0 * controls().piston.energy_tolerance_j,
    );
    near(
        piston.primary().gas().volume_m3() + piston.opposing().volume_m3(),
        0.002,
        1e-16,
    );
}

#[test]
fn both_acoustic_bulk_moduli_predict_shared_small_signal_period() {
    let mut piston = pair(0.25, 0.0, 0.01, 0.001);
    let stiffness: f64 = [
        (piston.primary().gas(), 0.01_f64),
        (piston.opposing(), 0.01_f64),
    ]
    .iter()
    .map(|(gas, area)| {
        let medium = gas.state().state();
        medium.density * medium.sound_speed.powi(2) * area.powi(2) / gas.volume_m3()
    })
    .sum();
    let omega = (stiffness / piston.primary().mass_kg()).sqrt();
    let input = DoublePistonInput {
        duration_s: core::f64::consts::TAU / omega / 256.0,
        ..drive_pair([0.0, 0.0])
    };
    for _ in 0..256 {
        check(&piston.step(input, controls(), || false).unwrap());
    }
    assert!((piston.primary().velocity_m_s() / 0.001 - 1.0).abs() < 2e-3);
}

#[test]
fn source_geometry_and_numerical_refusals_leave_both_owners_unchanged() {
    let mut piston = pair(0.25, 0.001, 0.01, 0.0);
    let before = piston.clone();
    for input in [
        drive_pair([1e9, 0.0]),
        drive_pair([0.0, 1e9]),
        DoublePistonInput {
            duration_s: 0.0,
            ..drive_pair([1.0, 0.0])
        },
        DoublePistonInput {
            external_force_n: f64::NAN,
            ..drive_pair([0.0, 1.0])
        },
    ] {
        assert!(piston.step(input, controls(), || false).is_err());
        assert_eq!(piston, before);
    }
    for ctrl in [
        DoublePistonSolve {
            volume_tolerance_m3: 0.0,
            ..controls()
        },
        DoublePistonSolve {
            piston: PistonSolve {
                volume_bracket_m3: [0.0007, 0.003],
                ..single_controls()
            },
            ..controls()
        },
        DoublePistonSolve {
            opposing_gas: GasEnergySolve {
                maximum_temperature_k: 300.0,
                ..single_controls().gas
            },
            ..controls()
        },
    ] {
        assert!(
            piston
                .step(drive_pair([20.0, 0.0]), ctrl, || false)
                .is_err()
        );
        assert_eq!(piston, before);
    }
    // At the common 50 us tick the 20 J trial converges in one Newton
    // correction. Resolve the same piston at 200 us (omega*dt < 0.04) and
    // supply enough heat to require another correction, still inside both
    // material temperature domains. Prove that this is work exhaustion,
    // rather than a source, geometry or time-resolution refusal.
    let nonlinear = DoublePistonInput {
        duration_s: 0.0002,
        ..drive_pair([300.0, 0.0])
    };
    let omega2: f64 = [piston.primary().gas(), piston.opposing()]
        .into_iter()
        .map(|gas| {
            let medium = gas.state().state();
            medium.density * medium.sound_speed.powi(2) * 0.01_f64.powi(2)
                / (gas.volume_m3() * piston.primary().mass_kg())
        })
        .sum();
    assert!(omega2.sqrt() * nonlinear.duration_s < 0.04);
    let mut resolved = before.clone();
    let accepted = resolved.step(nonlinear, controls(), || false).unwrap();
    assert!(accepted.iterations > 1);
    check(&accepted);
    let exhausted = DoublePistonSolve {
        piston: PistonSolve {
            maximum_iterations: 1,
            ..single_controls()
        },
        ..controls()
    };
    assert!(matches!(
        piston.step(nonlinear, exhausted, || false),
        Err(PistonError::NotConverged { iterations: 1, .. })
    ));
    assert_eq!(piston, before);
    let tiny_bracket = DoublePistonSolve {
        piston: PistonSolve {
            volume_bracket_m3: [0.001 - 1e-10, 0.001 + 1e-10],
            ..single_controls()
        },
        ..controls()
    };
    assert!(matches!(
        piston.step(drive_pair([20.0, 0.0]), tiny_bracket, || false),
        Err(PistonError::OutsideBracket { .. })
    ));
    assert_eq!(piston, before);
}

#[test]
fn cancellation_including_last_publication_gate_is_atomic_for_the_pair() {
    let original = pair(0.25, 0.001, 0.01, 0.0);
    let mut probe = original.clone();
    let mut polls = 0;
    probe
        .step(drive_pair([20.0, 5.0]), controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    for cancel_at in [1, polls / 2, polls - 1, polls] {
        let mut piston = original.clone();
        let mut count = 0;
        assert!(matches!(
            piston.step(drive_pair([20.0, 5.0]), controls(), || {
                count += 1;
                count == cancel_at
            }),
            Err(PistonError::Cancelled)
        ));
        assert_eq!(piston, original);
    }
}

#[test]
fn invalid_opposing_face_is_refused_without_inventing_geometry() {
    for area in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            DoubleActingGasPiston::new(
                GasPiston::new(gas(0.75, 0.001), 1.0, 0.01, 0.0).unwrap(),
                gas(0.25, 0.0),
                area,
            )
            .is_err()
        );
    }
}
