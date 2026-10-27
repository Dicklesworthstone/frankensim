use super::*;

#[test]
fn independent_exterior_drives_both_feedback_and_receiver_without_replacing_bore_gas() {
    let text = SOURCE.replace("reed 0.0004 0.013 6000 0 0 0.35",
        "reed 0.0004 0.013 6000 0.00001 100 0.35");
    // The receiver is compared with its lossless Rayleigh owner, so the
    // comparison declares the lossless observer law (ruling 6283).
    let explicit = text.replace(RECEIVER,
        &format!("{RECEIVER}\nexterior 303.15 100000 0.4\nobserver_law none"));
    let mut original = ReedPerformance::from_bytes(text.as_bytes(), 37).unwrap();
    let mut changed = ReedPerformance::from_bytes(explicit.as_bytes(), 37).unwrap();
    let external = GasState::try_new_moist_air(303.15, 100000.0, 0.4).unwrap();
    let info = changed.outlet_receiver().unwrap();
    assert_eq!(info.gas, external);
    assert_ne!(info.gas, changed.segment_gases()[1]);
    assert_eq!(original.segment_gases(), changed.segment_gases());
    assert_eq!(voice(&original).rho, voice(&changed).rho);
    assert_eq!(voice(&original).zc, voice(&changed).zc);
    assert_eq!(original.outlet_receiver().unwrap().axial_delay_samples, info.axial_delay_samples);
    let duct = Duct { segments: vec![
        Segment::Cylinder { radius: 0.0015, length: 0.2 },
        Segment::Cylinder { radius: 0.0015, length: 0.3 },
    ] };
    let bound = RegionalDuct::new(&duct, changed.segment_gases()).unwrap()
        .with_exterior_gas(external).unwrap();
    let mut reference = bound.baffled_outlet(48000, info.fft_samples, info.receiver, None).unwrap();
    reference.step(voice(&changed).p_plus_prev).unwrap();
    let mut motion_changed = false;
    for _ in 0..701 {
        original.pressure_block(&mut [0.0]).unwrap();
        let mut actual = [0.0];
        changed.pressure_block(&mut actual).unwrap();
        assert_eq!(actual[0].to_bits(), reference.step(voice(&changed).p_plus_prev).unwrap().to_bits());
        motion_changed |= voice(&original).reed_y.to_bits() != voice(&changed).reed_y.to_bits();
    }
    assert!(motion_changed, "exterior must change the feedback load, not only microphone gain");
    assert_eq!(raw(&explicit, 37), raw(&explicit, 113));
}

#[test]
fn identical_exterior_is_exact_and_unsupported_interfaces_refuse_with_source_lines() {
    let same = SOURCE.replace(RECEIVER, &format!(
        "{RECEIVER}\nexterior mixture 293.15 100000 mole 3 N2 0.4 O2 0.2 CO2 0.4"));
    assert_eq!(raw(&SOURCE, 37), raw(&same, 113));
    for record in [
        "exterior mixture 300 100000 mole 1 Xe 1",
        "exterior mixture 300 100000 mole 2 N2 70 O2 30",
        "exterior mixture 300 100000 mole 1 H2O 1",
    ] {
        let text = SOURCE.replace(RECEIVER, &format!("{RECEIVER}\n{record}"));
        let error = ReedPerformance::from_bytes(text.as_bytes(), 37).err().unwrap();
        assert!(matches!(&error, ReedPerformanceError::Gas { line: 6, .. }));
        assert!(std::error::Error::source(&error).is_some());
    }
    for text in [
        SOURCE.replace(RECEIVER, &format!("{RECEIVER}\nexterior 300 90000 0")),
        SOURCE.replace(RECEIVER, "listener_m 0.5\nexterior 300 100000 0"),
        SOURCE.replace(RECEIVER, &format!("{RECEIVER}\nexterior 300 100000 0\nexterior 300 100000 0")),
        SOURCE.replace(RECEIVER, &format!("{RECEIVER}\nexterior mixture NaN 100000 mole 1 N2 1")),
    ] { assert!(ReedPerformance::from_bytes(text.as_bytes(), 37).is_err()); }
}

#[test]
fn g3_two_medium_open_and_vented_chimneys_reach_feedback_and_receiver() {
    for sigma in [1.0, 0.3] {
        let geometry = Duct { segments: vec![
            Segment::Cylinder { radius: 0.0015, length: 0.1 },
            Segment::ToneHole { hole_radius: 0.0005, chimney_height: 0.002,
                bore_radius: 0.0015, state: HoleState::Vent(sigma) },
            Segment::Cylinder { radius: 0.0015, length: 0.1 },
            Segment::Cylinder { radius: 0.0015, length: 0.3 },
        ] };
        let text = SOURCE.replace("segments 2\ncylinder 0.0015 0.2",
            &format!("segments 4\ncylinder 0.0015 0.1\nhole 0.0005 0.002 0.0015 {sigma}\ncylinder 0.0015 0.1"))
            .replace("reed 0.0004 0.013 6000 0 0 0.35", "reed 0.0004 0.013 6000 0.00001 100 0.35")
            .replace(RECEIVER, &format!("{RECEIVER}\nobserver_law none"));
        let explicit = text.replace(RECEIVER, &format!("{RECEIVER}\nexterior 303.15 100000 0.4"));
        let mut baseline = ReedPerformance::from_bytes(text.as_bytes(), 37).unwrap();
        let mut actual = ReedPerformance::from_bytes(explicit.as_bytes(), 37).unwrap();
        let info = actual.outlet_receiver().unwrap();
        assert_eq!(baseline.segment_gases(), actual.segment_gases());
        assert_ne!(info.gas, actual.segment_gases()[1]);
        assert_eq!(voice(&baseline).rho, voice(&actual).rho);
        assert_eq!(voice(&baseline).zc, voice(&actual).zc);
        let bound = RegionalDuct::new(&geometry, actual.segment_gases()).unwrap()
            .with_exterior_gas(info.gas).unwrap();
        let (mut bell, mut holes) = bound.baffled_outlet_and_open_holes_delayed_at(
            48000, info.fft_samples, info.receiver, fs_duct::WallModel::Rigid, info.gas.sound_speed,
        ).unwrap();
        assert_eq!(holes.len(), 1);
        assert_eq!(holes[0].info().gas, info.gas);
        bell.step(voice(&actual).p_plus_prev).unwrap();
        holes[0].step(voice(&actual).p_plus_prev).unwrap();
        let (mut motion_changed, mut sound_changed, mut hole_audible) = (false, false, false);
        for _ in 0..1201 {
            let mut original = [0.0];
            baseline.pressure_block(&mut original).unwrap();
            let mut sample = [0.0];
            actual.pressure_block(&mut sample).unwrap();
            assert!(original[0].is_finite() && sample[0].is_finite());
            let outgoing = voice(&actual).p_plus_prev;
            let bell_pa = bell.step(outgoing).unwrap();
            let hole_pa = holes[0].step(outgoing).unwrap();
            assert_eq!(sample[0].to_bits(), (bell_pa + hole_pa).to_bits());
            motion_changed |= voice(&baseline).reed_y.to_bits() != voice(&actual).reed_y.to_bits();
            sound_changed |= original[0].to_bits() != sample[0].to_bits();
            hole_audible |= hole_pa.abs() > 1e-9;
        }
        assert!(motion_changed && sound_changed, "room material must alter feedback and sound");
        assert!(hole_audible, "the independently loaded chimney must reach the microphone");
        assert!(actual.last_outlet_error().is_none());
        assert_eq!(raw(&explicit, 37), raw(&explicit, 113));
    }
}
