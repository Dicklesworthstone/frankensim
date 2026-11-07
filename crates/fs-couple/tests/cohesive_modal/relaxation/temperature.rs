//! G2 source refresh changes future native force without resetting bond memory.
use super::*;
use fs_couple::cohesive_contact::CohesiveLoadLine;
use fs_couple::render::schedule::force::coupled::contact::cohesive::relaxation::temperature::ViscoelasticBondModalTemperatureReport;

fn thermal_source(card: &fs_matdb::InterfaceSystemCard) -> ResolvedViscoelasticBond {
    resolve_viscoelastic_bond(
        card,
        &QueryPoint::new().with("T", 280.0).unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
        material::limits(),
    )
    .unwrap()
}

#[test]
fn source_refresh_preserves_native_motion_and_uses_retained_memory_on_next_step() {
    let card = material::temperature_card();
    let cold = thermal_source(&card);
    let hot = cold.at_temperature(&card, "T", 320.0).unwrap();
    let mut system = build(&cold, [0.001, -0.0005, 0.0003], false);
    for _ in 0..40 {
        system.step(&[0.0; 6], [0.0; 3]).unwrap();
    }
    let native = native_states(&system);
    let old_point = point_observation(&system);
    let last = system.last_frame().copied();
    let before = system.total_energy_j().unwrap();
    let expected = system
        .point()
        .transfer_temperature(hot.law(), 1e-8)
        .unwrap();
    let report = system
        .apply_temperature_successor(&hot, "T", &CancelGate::new())
        .unwrap();
    assert_eq!(native_states(&system), native);
    assert_eq!(system.samples_rendered(), 40);
    assert_eq!(system.last_frame().copied(), last);
    assert_eq!(system.point().jump_m(), old_point.0);
    assert_eq!(system.point().elapsed_s(), old_point.3);
    assert!(core::ptr::eq(system.source(), &hot));
    assert!(report.parameter_work_j < 0.0);
    close(
        report.parameter_work_j,
        AREA * expected.parameter_work_j_m2,
        1e-15,
    );
    close(
        system.total_energy_j().unwrap() - before,
        report.parameter_work_j,
        1e-13,
    );
    assert!(report.energy_residual_j.abs() <= report.energy_tolerance_j);
    assert_ne!(system.point().traction_pa(), old_point.1);

    // Independently apply the original load line to retained transferred branch
    // memory and the two native bodies' actual unforced free flight.
    let dt = system.sample_period_s();
    let relative_velocity: [f64; 3] = core::array::from_fn(|i| {
        native[3 + i].velocity_m_sqrt_kg_per_s / det::sqrt(MB)
            - native[i].velocity_m_sqrt_kg_per_s / det::sqrt(MA)
    });
    let free = core::array::from_fn(|i| old_point.0[i] + dt * relative_velocity[i]);
    let mut line =
        CohesiveLoadLine::new(AREA, system.compliance_m_per_n(), config().contact).unwrap();
    let ramp = hot.law().prepare(dt).unwrap();
    let anticipated = line
        .trial_relaxation(&expected.point, &ramp, free, || true)
        .unwrap();
    let frame = *system.step(&[0.0; 6], [0.0; 3]).unwrap();
    for i in 0..3 {
        close(
            frame.jump_m[i],
            anticipated.material.point.jump_m()[i],
            1e-13,
        );
        close(frame.mean_force_n[i], anticipated.mean_force_n[i], 1e-9);
        close(
            frame.endpoint_traction_pa[i],
            anticipated.material.traction_pa[i],
            1e-4,
        );
    }
    close(
        frame.relaxation_heat_j,
        AREA * anticipated.material.viscous_dissipated_j_m2,
        1e-13,
    );
}

#[test]
fn temperature_receiver_and_foreign_source_refusals_are_atomic_and_replayable() {
    let card = material::temperature_card();
    let cold = thermal_source(&card);
    let hot = cold.at_temperature(&card, "T", 320.0).unwrap();
    let foreign = source();
    let mut system = build(&cold, [0.001, -0.0005, 0.0003], false);
    let mut replay = build(&cold, [0.001, -0.0005, 0.0003], false);
    for _ in 0..20 {
        for model in [&mut system, &mut replay] {
            model.step(&[0.0; 6], [0.0; 3]).unwrap();
        }
    }
    let native = native_states(&system);
    let point = point_observation(&system);
    let last = system.last_frame().copied();
    let gate = CancelGate::new();
    let refused: Result<(), CohesiveModalError> = system
        .apply_temperature_successor_with_admission(&hot, "T", &gate, |report, candidate| {
            assert!(report.parameter_work_j < 0.0);
            assert_eq!(candidate.jump_m(), point.0);
            Err(CohesiveModalError::Network(
                ModalCouplingError::EnergyBalance {
                    residual_j: 1.0,
                    tolerance_j: 0.0,
                },
            ))
        });
    assert!(refused.is_err());
    let cancelled: Result<ViscoelasticBondModalTemperatureReport, CohesiveModalError> = system
        .apply_temperature_successor_with_admission(&hot, "T", &gate, |report, _| {
            gate.request();
            Ok(*report)
        });
    assert!(cancelled.is_err());
    assert!(
        system
            .apply_temperature_successor(&foreign, "T", &CancelGate::new())
            .is_err()
    );
    assert_eq!(native_states(&system), native);
    assert_eq!(point_observation(&system), point);
    assert_eq!(system.last_frame().copied(), last);
    assert!(core::ptr::eq(system.source(), &cold));
    for model in [&mut system, &mut replay] {
        model
            .apply_temperature_successor(&hot, "T", &CancelGate::new())
            .unwrap();
    }
    for _ in 0..50 {
        assert_eq!(
            *system.step(&[0.0; 6], [0.0; 3]).unwrap(),
            *replay.step(&[0.0; 6], [0.0; 3]).unwrap()
        );
        assert_eq!(native_states(&system), native_states(&replay));
    }
}
