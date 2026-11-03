//! Original contact-source changes preserve real native motion and work.
use super::*;
#[path = "../../../fs-contact/tests/adhesive_temperature/common.rs"]
mod thermal;
use fs_contact::interface_binding::adhesive_shear::SourceAdhesiveShear;
use fs_contact::normal_patch::adhesive::temperature::AdhesiveTemperatures;
use fs_contact::normal_patch::dmt::coupled::DmtCoupledError;
use fs_material::state_point::MaterialPropertySelection;

fn temperatures() -> AdhesiveTemperatures {
    AdhesiveTemperatures {
        surface_a_k: 310.0,
        surface_b_k: 290.0,
        interface_k: 315.0,
    }
}

// Direct dimensional contact/gap/separated potential, independently using
// the manufactured original material curves rather than the refresh receipt.
fn potential(delta: f64, t: AdhesiveTemperatures) -> f64 {
    let ea = 1e11 * (1.0 - 0.01 * (t.surface_a_k - 300.0));
    let eb = 2e11 * (1.0 - 0.01 * (t.surface_b_k - 300.0));
    let e = 1.0 / ((1.0 - 0.25_f64.powi(2)) / ea + (1.0 - 0.25_f64.powi(2)) / eb);
    let w = 0.01 * (1.0 - 0.03 * (t.interface_k - 300.0));
    let r: f64 = 1e-6;
    let pull = 2.0 * core::f64::consts::PI * r * w;
    if delta >= 0.0 {
        (8.0 / 15.0) * e * r.sqrt() * delta.powf(2.5) - pull * (delta + 0.5 * Z)
    } else if delta > -Z {
        -0.5 * pull * (delta + Z).powi(2) / Z
    } else {
        0.0
    }
}

fn mean_area(a: f64, b: f64) -> f64 {
    let approach = if a >= 0.0 && b >= 0.0 {
        0.5 * (a + b)
    } else if a <= 0.0 && b <= 0.0 {
        0.0
    } else {
        0.5 * a.max(b).powi(2) / (a - b).abs()
    };
    core::f64::consts::PI * 1e-6 * approach
}

#[test]
fn prescribed_temperature_preserves_native_history_and_drives_actual_area_shear() {
    let fixture = thermal::Fixture::with_options(1e5, 0.01, true, 320.0);
    let normal = fixture.dmt(Z);
    let shear = SourceAdhesiveShear::resolve(
        normal.normal_state(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let cold = AdhesiveTemperatures {
        surface_a_k: 300.0,
        surface_b_k: 300.0,
        interface_k: 300.0,
    };
    for initial_delta in [0.2 * Z, -0.4 * Z, -2.0 * Z] {
        let mut system = make(
            &normal,
            &shear,
            [0.0, 1e-3, -0.5e-3],
            [initial_delta, 0.0, 0.0],
            true,
            RATE,
        );
        system.step(&[0.0; 6], [0.0; 3]).unwrap();
        let checkpoint = (
            states(&system),
            system.state().position_m(),
            system.state().regime(),
            system.samples_rendered(),
            *system.last_frame().unwrap(),
            system.work_columns().clone(),
            system.compliance_m_per_n(),
            system.config(),
        );
        let delta = system.state().position_m()[0];
        let before = system.total_energy_j().unwrap();
        let refreshed = normal
            .refresh_at_held_indentation(delta, fixture.cards(), temperatures())
            .unwrap();
        let hot_shear = shear.rebind_dmt_temperature(&normal, &refreshed).unwrap();
        let report = system
            .apply_source_refresh(&refreshed, &hot_shear, &CancelGate::new_clock_free())
            .unwrap();
        assert_eq!(
            (
                states(&system),
                system.state().position_m(),
                system.state().regime(),
                system.samples_rendered(),
                *system.last_frame().unwrap(),
                system.work_columns().clone(),
                system.compliance_m_per_n(),
                system.config()
            ),
            checkpoint,
        );
        assert!(core::ptr::eq(system.normal(), refreshed.source()));
        assert!(core::ptr::eq(system.shear(), &hot_shear));
        assert_eq!(report.sample, 1);
        assert_eq!(report.energy_before_j, before);
        assert_eq!(report.energy_after_j, system.total_energy_j().unwrap());
        close(
            report.surface_energy_before_j,
            potential(delta, cold),
            1e-30,
            2e-13,
        );
        close(
            report.surface_energy_after_j,
            potential(delta, temperatures()),
            1e-30,
            2e-13,
        );
        close(
            report.parameter_work_j,
            potential(delta, temperatures()) - potential(delta, cold),
            1e-30,
            2e-13,
        );
        assert!(report.energy_residual_j.abs() <= report.energy_tolerance_j);

        let energy = system.total_energy_j().unwrap();
        let mut dissipation = 0.0;
        for _ in 0..8 {
            let old = system.state().position_m();
            let frame = *system.step(&[0.0; 6], [0.0; 3]).unwrap();
            let area = mean_area(old[0], frame.position_m[0]);
            close(frame.path.mean_real_area_m2, area, 1e-30, 2e-13);
            close(frame.path.static_capacity_n, 2550.0 * area, 1e-30, 2e-13);
            close(frame.path.kinetic_capacity_n, 1700.0 * area, 1e-30, 2e-13);
            close(
                frame.friction_dissipation_j,
                1700.0 * area * (frame.position_m[1] - old[1]).hypot(frame.position_m[2] - old[2]),
                1e-30,
                2e-8,
            );
            dissipation += frame.friction_dissipation_j + frame.network_dissipation_j;
            close(
                system.total_energy_j().unwrap() + dissipation,
                energy,
                1e-25,
                2e-8,
            );
        }
        // A second transfer follows the retained source, at the now actual gap;
        // no contact reset or reconstructed modal network is needed to cool.
        let delta = system.state().position_m()[0];
        let cooling = system
            .normal()
            .refresh_at_held_indentation(delta, fixture.cards(), cold)
            .unwrap();
        let cold_shear = system
            .shear()
            .rebind_dmt_temperature(system.normal(), &cooling)
            .unwrap();
        let old_native = states(&system);
        let old_position = system.state().position_m();
        let old_frame = *system.last_frame().unwrap();
        let cooled = system
            .apply_source_refresh(&cooling, &cold_shear, &CancelGate::new_clock_free())
            .unwrap();
        assert_eq!(states(&system), old_native);
        assert_eq!(system.state().position_m(), old_position);
        assert_eq!(*system.last_frame().unwrap(), old_frame);
        assert_eq!(cooled.sample, 9);
        close(
            cooled.parameter_work_j,
            potential(delta, cold) - potential(delta, temperatures()),
            1e-30,
            2e-13,
        );
        assert!(cooled.energy_residual_j.abs() <= cooled.energy_tolerance_j);
        system.step(&[0.0; 6], [0.0; 3]).unwrap();
    }
}

#[test]
fn rejected_temperature_transfers_leave_native_motion_sources_and_retry_unchanged() {
    let fixture = thermal::Fixture::with_options(1e5, 0.01, false, 320.0);
    let normal = fixture.dmt(Z);
    let shear = SourceAdhesiveShear::resolve(
        normal.normal_state(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let mut system = make(
        &normal,
        &shear,
        [0.0, 1e-3, 0.0],
        [0.2 * Z, 0.0, 0.0],
        false,
        RATE,
    );
    let mut replay = make(
        &normal,
        &shear,
        [0.0, 1e-3, 0.0],
        [0.2 * Z, 0.0, 0.0],
        false,
        RATE,
    );
    system.step(&[0.0; 6], [0.0; 3]).unwrap();
    replay.step(&[0.0; 6], [0.0; 3]).unwrap();
    let checkpoint = (
        states(&system),
        system.state(),
        system.samples_rendered(),
        *system.last_frame().unwrap(),
        system.total_energy_j().unwrap(),
    );
    let delta = system.state().position_m()[0];
    let refreshed = normal
        .refresh_at_held_indentation(delta, fixture.cards(), temperatures())
        .unwrap();
    let hot_shear = shear.rebind_dmt_temperature(&normal, &refreshed).unwrap();
    let stale = normal
        .refresh_at_held_indentation(delta + 0.01 * Z, fixture.cards(), temperatures())
        .unwrap();
    let held_load = normal
        .refresh_at_held_load(
            normal.evaluate(delta).unwrap().normal_force_n,
            fixture.cards(),
            temperatures(),
        )
        .unwrap();
    for (refresh, next_shear) in [
        (&stale, &hot_shear),
        (&refreshed, &shear),
        (&held_load, &hot_shear),
    ] {
        assert!(matches!(
            system.apply_source_refresh(refresh, next_shear, &CancelGate::new_clock_free()),
            Err(AdhesiveDmtModalError::Contact(
                DmtCoupledError::StateMismatch
            ))
        ));
        assert_eq!(
            (
                states(&system),
                system.state(),
                system.samples_rendered(),
                *system.last_frame().unwrap(),
                system.total_energy_j().unwrap()
            ),
            checkpoint
        );
        assert!(core::ptr::eq(system.normal(), &normal));
        assert!(core::ptr::eq(system.shear(), &shear));
    }
    let gate = CancelGate::new_clock_free();
    gate.request();
    assert!(matches!(
        system.apply_source_refresh(&refreshed, &hot_shear, &gate),
        Err(AdhesiveDmtModalError::Network(
            ModalCouplingError::Cancelled
        ))
    ));
    assert_eq!(
        (
            states(&system),
            system.state(),
            system.samples_rendered(),
            *system.last_frame().unwrap(),
            system.total_energy_j().unwrap()
        ),
        checkpoint
    );
    let report = system
        .apply_source_refresh(&refreshed, &hot_shear, &CancelGate::new_clock_free())
        .unwrap();
    let replay_report = replay
        .apply_source_refresh(&refreshed, &hot_shear, &CancelGate::new_clock_free())
        .unwrap();
    assert_eq!(report, replay_report);
    assert_eq!(
        system.step(&[0.0; 6], [0.0; 3]).unwrap(),
        replay.step(&[0.0; 6], [0.0; 3]).unwrap()
    );
    assert_eq!(states(&system), states(&replay));
}
