use super::*;

#[test]
fn pressure_source_advances_the_original_zone_and_live_bore_on_both_orientations() {
    let gate = CancelGate::new_clock_free();
    let reference = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100000.0).unwrap();
    for high in [false, true] {
        for substeps in [1, 3] {
            let mut coupled = system(100000.0, high, substeps, DgBoundary::Matched);
            let mut original = zone(GasSpec::dry_air_ussa1976(), reference, high, substeps);
            let side = usize::from(high);
            for pressure in [2.0, 2.0, -1.0, 0.0] {
                let drive = NsfBoundaryDrive::PressureOffset(pressure);
                let frame = coupled.step_with_boundary(drive, &gate).unwrap();
                assert_eq!(frame.drive, NsfDgDrive::Boundary(drive));
                let mut controls = [drive; 2];
                controls[side] =
                    NsfBoundaryDrive::IncomingCharacteristic(frame.incoming_characteristic_m_s);
                let trial = original
                    .preview_boundary_interval(controls, substeps, CELLS * substeps, || false)
                    .unwrap();
                assert_eq!(trial.boundary_exchange(), &frame.faces);
                trial.accept();
                assert_eq!(original.state_hash(), coupled.zone().state_hash());
                let source = frame.faces[1 - side];
                assert!((source.mean_acoustic_pressure_pa - pressure).abs() < 1e-9);
                let work =
                    pressure * source.normal_volume_flow_m3_s * coupled.load().propagator().dt();
                assert!(
                    (source.acoustic_pressure_work_out_j - work).abs()
                        <= 2e-12 * work.abs().max(1e-24)
                );
                assert!(frame.pressure_residual_pa.abs() < 1e-6);
            }
        }
    }
}

#[test]
fn pressure_is_fixed_while_returning_bore_waves_and_material_change_source_flow() {
    let gate = CancelGate::new_clock_free();
    let mut matched = system(100000.0, true, 1, DgBoundary::Matched);
    let mut closed = system(100000.0, true, 1, DgBoundary::Closed);
    let mut dense = system(200000.0, true, 1, DgBoundary::Matched);
    let mut reflection_difference: f64 = 0.0;
    let mut material_difference: f64 = 0.0;
    for _ in 0..96 {
        let drive = NsfBoundaryDrive::PressureOffset(2.0);
        let a = matched.step_with_boundary(drive, &gate).unwrap();
        let b = closed.step_with_boundary(drive, &gate).unwrap();
        let c = dense.step_with_boundary(drive, &gate).unwrap();
        for frame in [&a, &b, &c] {
            assert!((frame.faces[0].mean_acoustic_pressure_pa - 2.0).abs() < 1e-9);
        }
        reflection_difference = reflection_difference
            .max((a.faces[0].normal_volume_flow_m3_s - b.faces[0].normal_volume_flow_m3_s).abs());
        material_difference = material_difference
            .max((a.faces[0].normal_volume_flow_m3_s - c.faces[0].normal_volume_flow_m3_s).abs());
    }
    assert!(
        reflection_difference > 1e-12,
        "source flow must respond to live bore feedback"
    );
    assert!(
        material_difference > 1e-12,
        "fixed pressure must use actual material impedance"
    );
    assert_ne!(matched.zone().state_hash(), closed.zone().state_hash());
    assert_ne!(matched.load().state(), dense.load().state());
}

#[test]
fn typed_playback_preserves_old_wave_bits_and_pressure_callback_identity() {
    let gate = Arc::new(CancelGate::new_clock_free());
    let mut old = NsfDgPerformance::new(
        system(100000.0, true, 1, DgBoundary::Matched),
        ApertureObservation::Inlet,
        RATE,
        vec![0.0; 8],
        8,
        Arc::clone(&gate),
    )
    .unwrap();
    let mut typed = NsfDgPerformance::new_with_boundary(
        system(100000.0, true, 1, DgBoundary::Matched),
        ApertureObservation::Inlet,
        RATE,
        vec![NsfBoundaryDrive::IncomingCharacteristic(0.0); 8],
        8,
        Arc::clone(&gate),
    )
    .unwrap();
    let mut a = [0.0; 8];
    let mut b = [0.0; 8];
    old.block(&mut a).unwrap();
    typed.block(&mut b).unwrap();
    assert_eq!(a, b);
    assert_eq!(
        old.system().zone().state_hash(),
        typed.system().zone().state_hash()
    );
    assert_eq!(old.system().load().state(), typed.system().load().state());

    let make = |block| {
        let sys = system(100000.0, true, 1, DgBoundary::Matched);
        let count = source_intervals_for(sys.load().propagator().dt(), 144000, 12).unwrap();
        NsfDgPerformance::new_resampled_with_boundary(
            sys,
            ApertureObservation::Inlet,
            144000,
            vec![NsfBoundaryDrive::PressureOffset(2.0); count],
            12,
            block,
            Arc::clone(&gate),
        )
        .unwrap()
    };
    let mut whole = make(12);
    let mut parts = make(5);
    let mut a = [0.0; 12];
    let mut b = [0.0; 12];
    whole.block(&mut a).unwrap();
    for chunk in b.chunks_mut(5) {
        parts.block(chunk).unwrap();
    }
    assert_eq!(a, b);
    assert_eq!(whole.drive_history(), parts.drive_history());
    assert_eq!(
        whole.system().zone().state_hash(),
        parts.system().zone().state_hash()
    );
    assert_eq!(whole.system().load().state(), parts.system().load().state());
}

#[test]
fn pressure_domain_and_cancellation_refusals_preserve_both_owners_and_drive() {
    let gate = CancelGate::new_clock_free();
    let mut coupled = system(100000.0, true, 1, DgBoundary::Matched);
    coupled
        .step_with_boundary(NsfBoundaryDrive::PressureOffset(2.0), &gate)
        .unwrap();
    let before = (
        coupled.zone().state_hash(),
        coupled.load().state().clone(),
        coupled.zone().boundary_drives(),
        coupled.accepted_steps(),
    );
    for pressure in [f64::NAN, 10000.0, -100000.0] {
        assert!(
            coupled
                .step_with_boundary(NsfBoundaryDrive::PressureOffset(pressure), &gate)
                .is_err()
        );
        assert_eq!(
            before,
            (
                coupled.zone().state_hash(),
                coupled.load().state().clone(),
                coupled.zone().boundary_drives(),
                coupled.accepted_steps()
            )
        );
    }
    gate.request();
    assert!(matches!(
        coupled.step_with_boundary(NsfBoundaryDrive::PressureOffset(3.0), &gate),
        Err(NsfDgError::Cancelled)
    ));
    assert_eq!(
        before,
        (
            coupled.zone().state_hash(),
            coupled.load().state().clone(),
            coupled.zone().boundary_drives(),
            coupled.accepted_steps()
        )
    );
}
