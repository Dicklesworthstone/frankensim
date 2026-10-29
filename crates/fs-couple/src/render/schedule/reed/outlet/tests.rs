use super::*;

const SOURCE_FIXTURE: &str = include_str!("../../../../../examples/reed-baffled-mixture.performance");
/// The fixture's physics without its sidecar declaration: these library tests
/// add their own observation records, and the library never refuses an
/// undeclared law (ruling 6532 puts that in the sidecar writer).
static SOURCE: std::sync::LazyLock<String> =
    std::sync::LazyLock::new(|| SOURCE_FIXTURE.replace("observer_law stokes-kirchhoff\n", ""));
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
        // The owner is the lossless Rayleigh transfer, so the comparison
        // declares the lossless observer law (ruling 6283).
        let text = SOURCE.replace("reed 0.0004 0.013 6000 0 0 0.35", replacement)
            .replace(RECEIVER, &format!("{RECEIVER}\nobserver_law none"));
        let mut physical = ReedPerformance::from_bytes(text.as_bytes(), 113).unwrap();
        assert_eq!(physical.observer().filter.absorption.law, crate::air_path::GasAbsorptionModel::Lossless);
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

/// Ruling 6324, conditions 3 and 4, for a dispersive law (the ISO companion
/// on a legacy exterior) at 0.5 m and 30 m. Group delays at DC are first
/// moments: `sum n k[n] / sum k[n]`.
/// - Baffled: the Rayleigh kernel, delayed at the receipt's speed, then the
///   observer FIR, totals `r / c_0`. The bound is the element-range spread over
///   `c_ref` (the kernel's DC delay is its harmonic-mean element range) plus
///   1e-6 sample of rounding. The receipt's delay IS the FIR's first moment,
///   so `d` and `D` enter exactly, with no quantization.
/// - Proxy: the FIR delayed by the typed consumer delay totals `r / c_0` to
///   rounding. A consumer that applies it therefore places the onset where
///   the baffled path does, within the sum of both bounds.
#[test]
fn dispersive_observer_delay_is_counted_once_on_both_observers() {
    use crate::pcm_wav::baffled::{BaffledPressure, RayleighMedium};
    let moment = |values: &[f64]| {
        let (first, zeroth) = values.iter().enumerate()
            .fold((0.0, 0.0), |(a, b), (n, v)| (a + n as f64 * v, b + v));
        first / zeroth
    };
    let legacy = "303.15 100000 0.4";
    for range in [0.5, 30.0] {
        let receiver_line = format!("baffled_outlet 0 0 {range} 4 16 1000");
        let baffled_text = SOURCE.replace(RECEIVER, &format!("{receiver_line}\nexterior {legacy}"));
        let baffled = ReedPerformance::from_bytes(baffled_text.as_bytes(), 37).unwrap();
        let proxy_text = SOURCE
            .replace("ambient mixture 293.15 100000 mole 2 N2 0.7 O2 0.3", &format!("ambient {legacy}"))
            .replace(RECEIVER, &format!("listener_m {range}"));
        let proxy = ReedPerformance::from_bytes(proxy_text.as_bytes(), 37).unwrap();
        let (b, q) = (baffled.observer(), proxy.observer());
        assert!(matches!(b.filter.absorption.law, crate::air_path::GasAbsorptionModel::Iso9613Air { .. }));
        assert_eq!(b.filter.taps, q.filter.taps, "the same law, gas, range and clock");
        assert_eq!(b.consumer_delay_s, None);
        let consumer_delay_s = q.consumer_delay_s.expect("the proxy types its consumer delay");
        let info = baffled.outlet_receiver().unwrap();
        let rate = f64::from(info.sample_rate_hz);
        let c0 = info.gas.sound_speed;
        let mut mic = BaffledPressure::circular_outlet(info.radius_m, info.sample_rate_hz, info.receiver,
            RayleighMedium { density: info.gas.density, sound_speed: b.filter.delay_speed_m_s.unwrap() }).unwrap();
        let kernel: Vec<f64> = (0..mic.delay_samples.1 + 4).map(|_| mic.step(&[1.0]).unwrap()).collect();
        let fir = baffled.reference_observer(&info.gas);
        let fir_delay = moment(fir.taps());
        let baffled_chain = (moment(&kernel) + fir_delay) / rate;
        let proxy_chain = consumer_delay_s + fir_delay / rate;
        let spread = b.range_spread_m.unwrap();
        let bound_baffled = spread / b.filter.delay_speed_m_s.unwrap() + 1e-6 / rate;
        let bound_proxy = 1e-6 / rate;
        println!(
            "{{\"bead\":\"pexfp.18\",\"test\":\"observer_delay_counted_once\",\"range_m\":{range},\"alignment_samples\":{},\"padding_samples\":{},\"latency_s\":{:e},\"r_over_c0_s\":{:e},\"baffled_chain_s\":{baffled_chain:e},\"proxy_chain_s\":{proxy_chain:e},\"bound_baffled_s\":{bound_baffled:e},\"bound_proxy_s\":{bound_proxy:e}}}",
            b.filter.alignment_samples, b.filter.padding_samples, b.filter.latency_s, range / c0
        );
        assert!((baffled_chain - range / c0).abs() <= bound_baffled, "baffled {baffled_chain} vs r/c0 {}", range / c0);
        assert!((proxy_chain - range / c0).abs() <= bound_proxy, "proxy {proxy_chain} vs r/c0 {}", range / c0);
        assert!((proxy_chain - baffled_chain).abs() <= bound_baffled + bound_proxy);
    }
}

/// Ruling 6332: the receiver-delay admission on explicit numbers refuses
/// with them, and admits at equality. Planted negative at admitted geometry:
/// the shortest admitted receiver (5 cm on axis, 48 kHz, a legacy exterior
/// under the ISO companion) cannot absorb its filter's delay and refuses with
/// its numbers; at 50 cm it is admitted.
#[test]
fn receiver_delay_admission_refuses_with_its_numbers() {
    let refused = super::super::check_receiver_delay(0.05, 7.0, 7.25, 1e-10).unwrap_err();
    assert_eq!(refused, AcousticRealizeError::ObserverLatency {
        range_m: 0.05, available_samples: 7.0, needed_samples: 7.25, budget_energy: 1e-10 });
    assert!(refused.to_string().contains("observer_law none"), "{refused}");
    assert!(super::super::check_receiver_delay(0.05, 7.0, 7.0, 1e-10).is_ok());
    let at = |z: f64| SOURCE.replace(RECEIVER,
        &format!("baffled_outlet 0 0 {z} 4 16 1000\nexterior 303.15 100000 0.4"));
    match ReedPerformance::from_bytes(at(0.05).as_bytes(), 37) {
        Err(ReedPerformanceError::Model(AcousticRealizeError::ObserverLatency { range_m, available_samples, needed_samples, budget_energy })) => {
            println!(
                "{{\"bead\":\"pexfp.18\",\"test\":\"receiver_delay_admission\",\"range_m\":{range_m},\"outcome\":\"refused\",\"available_samples\":{available_samples},\"needed_samples\":{needed_samples},\"budget_energy\":{budget_energy:e}}}"
            );
            assert!(available_samples < needed_samples);
        }
        Ok(_) => panic!("a 5 cm receiver at 48 kHz must refuse its observer filter's delay"),
        Err(other) => panic!("unexpected refusal {other}"),
    }
    let admitted = ReedPerformance::from_bytes(at(0.5).as_bytes(), 37).unwrap();
    let receipt = admitted.observer().filter;
    println!(
        "{{\"bead\":\"pexfp.18\",\"test\":\"receiver_delay_admission\",\"range_m\":0.5,\"outcome\":\"admitted\",\"alignment_samples\":{},\"padding_samples\":{},\"budget_energy\":{:e},\"low_frequency_delay_samples\":{}}}",
        receipt.alignment_samples, receipt.padding_samples, receipt.budget_energy, receipt.low_frequency_delay_s * 48_000.0
    );
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
    assert_eq!(resumed, raw(&SOURCE, 37));
    assert!(resumed[..first].iter().all(|p| *p == 0.0));
    let mut past_end = [321.0];
    assert!(source.pressure_block(&mut past_end).is_err());
    assert_eq!(past_end, [321.0]);
}

#[test]
fn gas_changes_reach_physical_pressure_and_invalid_observation_retains_its_error() {
    let base = raw(&SOURCE, 37);
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

#[test]
fn open_tone_hole_radiates_to_listener_before_bell_without_changing_in_loop_state() {
    // The hole goes between the fixture's two cylinders; the last keeps its regional gas suffix.
    let with_hole = SOURCE
        .replace(
            "segments 2\ncylinder 0.0015 0.2\ncylinder 0.0015 0.3",
            "segments 3\ncylinder 0.0015 0.2\nhole 0.0010 0.004 0.0015 1\ncylinder 0.0015 0.3",
        )
        .replace(RECEIVER, &format!("{RECEIVER}\nobserver_law none"));
    assert!(with_hole.contains("\nhole 0.0010 0.004 0.0015 1\n"), "the fixture's bore changed");
    let mut baffled = ReedPerformance::from_bytes(with_hole.as_bytes(), 1).unwrap();
    let mut legacy = ReedPerformance::from_bytes(
        with_hole.replace(RECEIVER, "listener_m 0.5").as_bytes(),
        1,
    )
    .unwrap();
    let duct = Duct {
        segments: vec![
            Segment::Cylinder {
                radius: 0.0015,
                length: 0.2,
            },
            Segment::ToneHole {
                hole_radius: 0.0010,
                chimney_height: 0.004,
                bore_radius: 0.0015,
                state: HoleState::Open,
            },
            Segment::Cylinder {
                radius: 0.0015,
                length: 0.3,
            },
        ],
    };
    let info = baffled.outlet_receiver().unwrap();
    let regional = RegionalDuct::new(&duct, baffled.segment_gases()).unwrap();
    let mut bell_only = regional
        .baffled_outlet(48000, info.fft_samples, info.receiver, None)
        .unwrap();
    let (mut bell, mut holes) = regional
        .baffled_outlet_and_open_holes_delayed_at(
            48000,
            info.fft_samples,
            info.receiver,
            fs_duct::WallModel::from_pin(None),
            info.gas.sound_speed,
        )
        .unwrap();
    assert_eq!(holes.len(), 1);
    assert!(holes[0].info().axial_delay_samples < bell.info().axial_delay_samples);
    bell_only.step(voice(&legacy).p_plus_prev).unwrap();
    bell.step(voice(&legacy).p_plus_prev).unwrap();
    holes[0].step(voice(&legacy).p_plus_prev).unwrap();
    let mut hole_differs_from_bell_only = false;
    for _ in 0..701 {
        legacy.pressure_block(&mut [0.0]).unwrap();
        let mut actual = [0.0];
        baffled.pressure_block(&mut actual).unwrap();
        let p_plus = voice(&legacy).p_plus_prev;
        let bell_pa = bell_only.step(p_plus).unwrap();
        let expected_sum = bell.step(p_plus).unwrap() + holes[0].step(p_plus).unwrap();
        assert_eq!(actual[0].to_bits(), expected_sum.to_bits());
        assert_eq!(voice(&legacy).reed_y.to_bits(), voice(&baffled).reed_y.to_bits());
        assert_eq!(voice(&legacy).reed_v.to_bits(), voice(&baffled).reed_v.to_bits());
        hole_differs_from_bell_only |= (actual[0] - bell_pa).abs() > 1e-9;
    }
    assert!(hole_differs_from_bell_only, "open hole must radiate to the listener");
}

#[path = "exterior_tests.rs"]
mod exterior;
