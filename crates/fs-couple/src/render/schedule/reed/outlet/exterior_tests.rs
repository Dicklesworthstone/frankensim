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

fn regional_mouth_fixture(law: &str) -> (String, Duct) {
    let text = SOURCE
        .replace("segments 2\ncylinder 0.0015 0.2",
            "segments 5\ncylinder 0.0015 0.1\nhole 0.0005 0.002 0.0015 1\ncylinder 0.0015 0.1\nhole 0.0008 0.003 0.0015 0.25")
        .replace("reed 0.0004 0.013 6000 0 0 0.35",
            "reed 0.0004 0.013 6000 0.00001 100 0.35")
        .replace(RECEIVER,
            &format!("{RECEIVER}\nexterior 303.15 100000 0.4\nobserver_law {law}"));
    let duct = Duct { segments: vec![
        Segment::Cylinder { radius: 0.0015, length: 0.1 },
        Segment::ToneHole { hole_radius: 0.0005, chimney_height: 0.002,
            bore_radius: 0.0015, state: HoleState::Open },
        Segment::Cylinder { radius: 0.0015, length: 0.1 },
        Segment::ToneHole { hole_radius: 0.0008, chimney_height: 0.003,
            bore_radius: 0.0015, state: HoleState::Vent(0.25) },
        Segment::Cylinder { radius: 0.0015, length: 0.3 },
    ] };
    (text, duct)
}

#[test]
fn regional_mouth_receipts_bind_open_and_vented_geometry_to_their_own_exterior_ranges() {
    let (text, _) = regional_mouth_fixture("stokes-kirchhoff");
    let performance = ReedPerformance::from_bytes(text.as_bytes(), 37).unwrap();
    let exterior = GasState::try_new_moist_air(303.15, 100000.0, 0.4).unwrap();
    let bell = performance.outlet_receiver().unwrap();
    let mouths: Vec<_> = performance.open_hole_observers().collect();
    assert_eq!(mouths.len(), 2);
    assert_eq!(performance.segment_gases().len(), 5);
    assert_ne!(performance.segment_gases()[1], exterior);
    assert_ne!(performance.segment_gases()[4], exterior);
    let legacy: Vec<_> = performance.molecular_open_hole_observers().collect();
    assert_eq!(legacy.len(), mouths.len());
    for (((info, receipt), (old_info, old_receipt)), (axial, radius)) in mouths.iter()
        .zip(&legacy).zip([(0.1, 0.0005), (0.2, 0.0008 * 0.25_f64.sqrt())]) {
        let position: [f64; 3] = [0.5 - axial, 0.0, 0.5];
        let range = position[0].hypot(position[2]);
        assert_eq!(info.gas, exterior);
        assert_eq!(info.radius_m, radius);
        assert_eq!(info.receiver.position_m, position);
        assert_eq!(info.sample_rate_hz, bell.sample_rate_hz);
        assert_eq!(info.fft_samples, bell.fft_samples);
        assert_eq!(info.pre_arrival_latency_samples, bell.pre_arrival_latency_samples);
        assert!(info.axial_delay_samples < bell.axial_delay_samples);
        assert_eq!(receipt.filter.range_m, range);
        assert!(range > performance.observer().filter.range_m);
        assert_eq!(receipt.consumer_delay_s, None);
        assert!(!receipt.internal_part_unpropagated);
        let expected = observer_filter(Some(&DeclaredObserverLaw::StokesKirchhoff),
            &exterior, performance.gas_records(), range, info.sample_rate_hz).unwrap();
        assert_eq!(receipt.filter, *expected.receipt());
        assert_eq!(receipt.range_spread_m, Some(range_spread(position, radius)));
        assert_eq!(old_info.receiver.position_m, info.receiver.position_m);
        assert_eq!(old_info.radius_m, info.radius_m);
        assert_eq!(old_receipt, receipt);
    }
    assert_ne!(mouths[0].1.filter.range_m, mouths[1].1.filter.range_m);
}

fn check_mouth_delay_once(info: OutletReceiverInfo, receipt: ReedObserverReceipt,
    filter: &ObserverFilter) {
    use crate::pcm_wav::baffled::{BaffledPressure, RayleighMedium};
    let moment = |values: &[f64]| {
        let (first, total) = values.iter().enumerate().fold((0.0, 0.0),
            |(first, total), (n, value)| (first + n as f64 * value, total + value));
        first / total
    };
    let speed = receipt.filter.delay_speed_m_s.unwrap();
    let mut microphone = BaffledPressure::circular_outlet(info.radius_m, info.sample_rate_hz,
        info.receiver, RayleighMedium { density: info.gas.density, sound_speed: speed }).unwrap();
    assert_eq!(info.exterior_delay_samples, microphone.delay_samples);
    let kernel: Vec<_> = (0..microphone.delay_samples.1 + 4)
        .map(|_| microphone.step(&[1.0]).unwrap()).collect();
    let rate = f64::from(info.sample_rate_hz);
    let chain_delay = (moment(&kernel) + moment(filter.taps())) / rate;
    let physical_delay = receipt.filter.range_m / info.gas.sound_speed;
    let bound = receipt.range_spread_m.unwrap() / speed + 1e-6 / rate;
    assert!((chain_delay - physical_delay).abs() <= bound,
        "mouth delay {chain_delay:e}, physical {physical_delay:e}, bound {bound:e}");
    assert!(receipt.filter.low_frequency_delay_s > 0.0);
    assert!(matches!(receipt.filter.absorption.law, GasAbsorptionModel::Iso9613Air { .. }));
}

#[test]
fn regional_iso_pressure_is_the_sum_of_independent_mouth_paths_without_feedback_changes() {
    let (text, duct) = regional_mouth_fixture("iso9613");
    let mut actual = ReedPerformance::from_bytes(text.as_bytes(), 37).unwrap();
    let mut lossless = ReedPerformance::from_bytes(
        text.replace("observer_law iso9613", "observer_law none").as_bytes(), 37).unwrap();
    assert_eq!(actual.segment_gases(), lossless.segment_gases());
    let bell_info = actual.outlet_receiver().unwrap();
    let bell_receipt = actual.observer();
    let mouths: Vec<_> = actual.open_hole_observers().collect();
    let mut bell_filter = observer_filter(Some(&DeclaredObserverLaw::Iso9613), &bell_info.gas,
        actual.gas_records(), bell_receipt.filter.range_m, bell_info.sample_rate_hz).unwrap();
    let mut filters: Vec<_> = mouths.iter().map(|(info, receipt)| {
        observer_filter(Some(&DeclaredObserverLaw::Iso9613), &info.gas,
            actual.gas_records(), receipt.filter.range_m, info.sample_rate_hz).unwrap()
    }).collect();
    assert_eq!(mouths.len(), 2);
    assert_eq!(bell_receipt.filter, *bell_filter.receipt());
    check_mouth_delay_once(bell_info, bell_receipt, &bell_filter);
    for ((info, receipt), filter) in mouths.iter().zip(&filters) {
        assert_eq!(receipt.filter, *filter.receipt());
        check_mouth_delay_once(*info, *receipt, filter);
    }
    let hole_speeds: Vec<_> = filters.iter()
        .map(|filter| filter.receipt().delay_speed_m_s.unwrap()).collect();
    let bound = RegionalDuct::new(&duct, actual.segment_gases()).unwrap()
        .with_exterior_gas(bell_info.gas).unwrap();
    let (mut bell, mut holes, _) = bound.baffled_outlet_and_open_holes_and_shunt_with_delay_speeds(
        bell_info.sample_rate_hz, voice(&actual).line.history().len(), bell_info.receiver,
        fs_duct::WallModel::Rigid, bell_filter.receipt().delay_speed_m_s.unwrap(),
        &hole_speeds, 1.0).unwrap();
    assert_eq!(holes.len(), filters.len());
    let prime = voice(&actual).p_plus_prev;
    bell.step(prime).unwrap();
    for hole in &mut holes { hole.step(prime).unwrap(); }
    let mut emitted = Vec::with_capacity(actual.info().samples as usize);
    let (mut hole_audible, mut absorption_audible) = (false, false);
    for _ in 0..actual.info().samples {
        let (mut received, mut unfiltered) = ([0.0], [0.0]);
        actual.pressure_block(&mut received).unwrap();
        lossless.pressure_block(&mut unfiltered).unwrap();
        let a = voice(&actual);
        let b = voice(&lossless);
        assert_eq!(a.reed_y.to_bits(), b.reed_y.to_bits());
        assert_eq!(a.reed_v.to_bits(), b.reed_v.to_bits());
        assert_eq!(a.p_plus_prev.to_bits(), b.p_plus_prev.to_bits());
        assert_eq!(a.line.incoming().to_bits(), b.line.incoming().to_bits());
        let mut expected = bell_filter.process(bell.step(a.p_plus_prev).unwrap());
        for (hole, filter) in holes.iter_mut().zip(&mut filters) {
            let contribution = filter.process(hole.step(a.p_plus_prev).unwrap());
            expected += contribution;
            hole_audible |= contribution.abs() > 1e-9;
        }
        assert_eq!(received[0].to_bits(), expected.to_bits());
        absorption_audible |= (received[0] - unfiltered[0]).abs() > 1e-9;
        emitted.push(received[0]);
    }
    assert!(hole_audible, "side-mouth propagation must contribute pressure");
    assert!(absorption_audible, "the declared exterior law must reach the receiver");
    assert_eq!(voice(&actual).line.history(), voice(&lossless).line.history());
    assert!(actual.last_outlet_error().is_none());
    assert_eq!((actual.remaining_samples(), lossless.remaining_samples()), (0, 0));
    assert_eq!(emitted, raw(&text, 37));
    assert_eq!(emitted, raw(&text, 113));
}
