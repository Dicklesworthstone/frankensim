use super::*;

const SOURCE: &str = include_str!("../../../../../examples/reed-baffled-mixture.performance");
const RECEIVER: &str = "baffled_outlet 0 0 0.5 4 16 1000";

fn voice(p: &ReedPerformance) -> &ReedBoreVoice {
    let RenderVoice::ReedBore(voice) = &p.renderer.context().voices[0] else { panic!("not a reed") };
    voice
}
fn raw(text: &str, block: usize) -> Vec<f64> {
    let mut p = ReedPerformance::from_bytes(text.as_bytes(), block).unwrap();
    let mut out = vec![0.0; p.info().samples as usize];
    for part in out.chunks_mut(block) { p.pressure_block(part).unwrap(); }
    assert_eq!(p.remaining_samples(), 0);
    assert!(out.iter().all(|p| p.is_finite()));
    assert!(out.iter().any(|p| *p != 0.0));
    out
}

#[test]
fn terminal_observation_matches_its_owner_without_changing_either_reed_model() {
    let duct = Duct { segments: vec![
        Segment::Cylinder { radius: 0.0015, length: 0.2 },
        Segment::Cylinder { radius: 0.0015, length: 0.3 },
    ] };
    for replacement in ["reed 0.0004 0.013 6000 0 0 0.35", "reed 0.0004 0.013 6000 0.00001 100 0.35"] {
        let text = SOURCE.replace("reed 0.0004 0.013 6000 0 0 0.35", replacement);
        let mut physical = ReedPerformance::from_bytes(text.as_bytes(), 113).unwrap();
        let proxy_text = text.replace(RECEIVER, "listener_m 0.5");
        let mut proxy = ReedPerformance::from_bytes(proxy_text.as_bytes(), 113).unwrap();
        assert!(proxy.outlet_receiver().is_none());
        let info = physical.outlet_receiver().unwrap();
        let mut expected = RegionalDuct::new(&duct, physical.segment_gases()).unwrap()
            .baffled_outlet(48000, info.fft_samples, info.receiver, None).unwrap();
        expected.step(voice(&proxy).p_plus_prev).unwrap();
        for _ in 0..701 {
            proxy.pressure_block(&mut [0.0]).unwrap();
            let mut actual = [0.0];
            physical.pressure_block(&mut actual).unwrap();
            assert_eq!(actual[0].to_bits(), expected.step(voice(&proxy).p_plus_prev).unwrap().to_bits());
            assert_eq!(voice(&proxy).reed_y.to_bits(), voice(&physical).reed_y.to_bits());
            assert_eq!(voice(&proxy).reed_v.to_bits(), voice(&physical).reed_v.to_bits());
        }
        assert_eq!(voice(&proxy).line.history(), voice(&physical).line.history());
        assert_eq!(proxy.samples_rendered(), physical.samples_rendered());
        assert_eq!(info.gas, physical.segment_gases()[1]);
        assert!(physical.last_outlet_error().is_none());
    }
}

#[test]
fn observer_choice_leaves_reed_and_bore_state_bitwise_identical_over_a_full_note() {
    // G3 (X-Obs1): the receiver only observes. Over the whole performance, the
    // reed opening, both travelling waves and the bore volume flow carry the
    // same bits under the legacy listener and the baffled outlet.
    for replacement in ["reed 0.0004 0.013 6000 0 0 0.35", "reed 0.0004 0.013 6000 0.00001 100 0.35"] {
        let text = SOURCE.replace("reed 0.0004 0.013 6000 0 0 0.35", replacement);
        let mut baffled = ReedPerformance::from_bytes(text.as_bytes(), 1).unwrap();
        let mut legacy = ReedPerformance::from_bytes(text.replace(RECEIVER, "listener_m 0.5").as_bytes(), 1).unwrap();
        assert!(baffled.outlet_receiver().is_some() && legacy.outlet_receiver().is_none());
        let samples = baffled.info().samples;
        assert_eq!(samples, legacy.info().samples);
        let flow = |v: &ReedBoreVoice| (v.p_plus_prev - v.line.incoming()) / v.zc;
        let mut moved = false;
        for _ in 0..samples {
            baffled.pressure_block(&mut [0.0]).unwrap();
            legacy.pressure_block(&mut [0.0]).unwrap();
            let (a, b) = (voice(&baffled), voice(&legacy));
            assert_eq!(a.reed_y.to_bits(), b.reed_y.to_bits());
            assert_eq!(a.reed_v.to_bits(), b.reed_v.to_bits());
            assert_eq!(a.p_plus_prev.to_bits(), b.p_plus_prev.to_bits());
            assert_eq!(a.line.incoming().to_bits(), b.line.incoming().to_bits());
            assert_eq!(flow(a).to_bits(), flow(b).to_bits());
            moved |= a.p_plus_prev != 0.0 && flow(a) != 0.0;
        }
        assert!(moved, "the note must actually drive the bore");
        assert_eq!(voice(&baffled).line.history(), voice(&legacy).line.history());
        assert_eq!((baffled.remaining_samples(), legacy.remaining_samples()), (0, 0));
    }
}

#[test]
fn delays_partitioning_and_request_refusals_preserve_the_complete_source() {
    let mut source = ReedPerformance::from_bytes(SOURCE.as_bytes(), 113).unwrap();
    let info = source.outlet_receiver().unwrap();
    let first = info.axial_delay_samples + info.exterior_delay_samples.0 - 1;
    let mut oversized = [123.0; 114];
    assert!(source.pressure_block(&mut oversized).is_err());
    assert_eq!(oversized, [123.0; 114]);
    assert!(source.pressure_block(&mut []).is_err());
    assert_eq!(source.samples_rendered(), 0);
    let mut resumed = vec![0.0; source.info().samples as usize];
    for block in resumed.chunks_mut(113) { source.pressure_block(block).unwrap(); }
    assert_eq!(resumed, raw(SOURCE, 37));
    assert!(resumed[..first].iter().all(|p| *p == 0.0));
    let mut past_end = [321.0];
    assert!(source.pressure_block(&mut past_end).is_err());
    assert_eq!(past_end, [321.0]);
}

#[test]
fn gas_changes_reach_physical_pressure_and_invalid_observation_retains_its_error() {
    let base = raw(SOURCE, 37);
    for text in [
        SOURCE.replace("gas mixture 293.15 100000 mole 3 N2 0.4 O2 0.2 CO2 0.4",
            "gas mixture 293.15 100000 mole 2 N2 0.7 O2 0.3"),
        SOURCE.replace("mixture 293.15 100000", "mixture 293.15 50000"),
    ] { assert_ne!(base, raw(&text, 113)); }
    for terminal in ["closed", "ideal-open", "unflanged"] {
        assert!(ReedPerformance::from_bytes(SOURCE.replace("termination flanged",
            &format!("termination {terminal}")).as_bytes(), 37).is_err());
    }
    let error = ReedPerformance::from_bytes(SOURCE.replace(RECEIVER,
        "baffled_outlet 0 0 0.01 4 16 1000").as_bytes(), 37).err().unwrap();
    assert!(matches!(&error, ReedPerformanceError::Outlet(_)));
    assert!(std::error::Error::source(&error).is_some());
}

#[path = "exterior_tests.rs"]
mod exterior;
