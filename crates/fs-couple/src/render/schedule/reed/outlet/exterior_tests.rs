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
fn g3_two_medium_open_chimney_reaches_the_scheduled_reed_and_microphone() {
    let holes = SOURCE.replace("segments 2\ncylinder 0.0015 0.2",
        "segments 4\ncylinder 0.0015 0.1\nhole 0.0005 0.002 0.0015 1\ncylinder 0.0015 0.1")
        .replace("reed 0.0004 0.013 6000 0 0 0.35",
            "reed 0.0004 0.013 6000 0.00001 100 0.35");
    let changed = holes.replace(RECEIVER,
        &format!("{RECEIVER}\nexterior 303.15 100000 0.4\nobserver_law none"));
    let baseline = holes.replace(RECEIVER, &format!("{RECEIVER}\nobserver_law none"));
    let mut a = ReedPerformance::from_bytes(baseline.as_bytes(), 37).unwrap();
    let mut b = ReedPerformance::from_bytes(changed.as_bytes(), 37).unwrap();
    assert_eq!(a.segment_gases(), b.segment_gases());
    assert_eq!(voice(&a).rho, voice(&b).rho);
    assert_eq!(voice(&a).zc, voice(&b).zc);
    let mut motion_changed = false;
    let mut sound_changed = false;
    for _ in 0..701 {
        let (mut pa, mut pb) = ([0.0], [0.0]);
        a.pressure_block(&mut pa).unwrap(); b.pressure_block(&mut pb).unwrap();
        assert!(pa[0].is_finite() && pb[0].is_finite());
        motion_changed |= voice(&a).reed_y.to_bits() != voice(&b).reed_y.to_bits();
        sound_changed |= pa[0].to_bits() != pb[0].to_bits();
    }
    assert!(motion_changed && sound_changed);
    assert_eq!(raw(&changed, 37), raw(&changed, 113));
}
