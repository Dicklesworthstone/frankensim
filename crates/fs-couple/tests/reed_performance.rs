//! Authored pressure phrases use the existing reed/TMM solver, not a fixture envelope.
use std::io::Cursor;
use fs_couple::pcm_wav::{encode_pcm16_wav, decimate::Decimator};
use fs_couple::pcm_wav::observation::{DecimatedRenderer, PressureRenderer};
use fs_couple::pcm_wav::stream::{Pcm16WavStream, ScheduledWavProgress, render_pressure_pcm16};
use fs_couple::render::ReedBoreVoice;
use fs_couple::render::schedule::reed::{ReedPerformance, MAX_REED_PERFORMANCE_BYTES};
use fs_couple::thin_plate::PlateBank;
use fs_duct::{Duct, HoleState, Segment, Termination};
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use fs_scenario::{BeatingReed, gesture::{GestureEvent, GestureSchedule, GestureTarget, GestureTrack, GestureValue}};

const INPUT: &str = include_str!("../examples/reed-duct.performance");
fn schedule() -> GestureSchedule {
    let value = |pressure_pa| GestureValue::PressurePa(pressure_pa);
    GestureSchedule::try_new(700, vec![GestureTrack {
        id: "blow".into(), target: GestureTarget::BlowingPressure,
        initial: value(0.0), events: vec![
            GestureEvent { time_s: 0.0, transition_s: 0.01, value: value(2800.0) },
            GestureEvent { time_s: 0.005, transition_s: 0.01, value: value(3500.0) },
            GestureEvent { time_s: 0.02, transition_s: 0.0, value: value(0.0) },
            GestureEvent { time_s: 0.04, transition_s: 0.005, value: value(2800.0) },
            GestureEvent { time_s: 0.06, transition_s: 0.0, value: value(0.0) },
        ],
    }]).unwrap()
}
fn geometry() -> Duct { Duct { segments: vec![Segment::Cylinder { radius: 0.0022, length: 0.5 }] } }
fn reference(rate: u32, samples: usize, duct: &Duct, temp: f64, mass: f64, stiffness: f64) -> Vec<f64> {
    let air = GasState::try_new_moist_air(temp, 101_325.0, 0.0).unwrap();
    let reed = BeatingReed { rest_opening_m: 0.0004, width_m: 0.013,
        closing_pressure_pa: 6000.0, mass_kg: mass, stiffness_n_m: stiffness,
        damping_ratio: 0.35, blowing_pressure_pa: 0.0, attack_s: 0.0 };
    let mut voice = ReedBoreVoice::new(duct, &air, reed, Termination::UnflangedOpen,
        PlateBank::default(), 1.0, rate, samples, None).unwrap();
    let gesture = schedule();
    let mut out = Vec::with_capacity(samples);
    // Independent sample-by-sample application, not the adapter's compiler or
    // callback splitter. Non-divisor control clocks exercise ceil placement.
    for sample in 0..samples as u64 {
        voice.set_blowing_pressure(gesture.sample("blow", sample * 700 / u64::from(rate)).unwrap());
        let mut pressure = [0.0]; voice.step_block(&mut pressure).unwrap();
        out.push(0.0 + pressure[0]); // same additive identity as the common mixer
    }
    out
}
fn render(text: &str, block: usize) -> Vec<f64> {
    let mut source = ReedPerformance::from_bytes(text.as_bytes(), block).unwrap();
    let mut out = vec![0.0; source.info().samples as usize];
    for chunk in out.chunks_mut(block) { source.block(chunk).unwrap(); }
    out
}
fn bits(values: &[f64]) -> Vec<u64> { values.iter().map(|v| v.to_bits()).collect() }

#[test]
fn canonical_pressure_phrase_matches_solver_and_preserves_release_across_callbacks() {
    assert_eq!(INPUT.split_once("\nschedule\n").unwrap().1.as_bytes(), schedule().to_canonical_bytes());
    let expected = reference(48_000, 4801, &geometry(), 293.15, 0.0, 0.0);
    assert!(expected[3000..].iter().any(|v| v.abs() > 1e-12), "release must not erase bore history");
    for partition in [1, 37, 512] {
        assert_eq!(bits(&render(INPUT, partition)), bits(&expected));
    }
}

#[test]
fn admitted_massive_reed_uses_the_same_retained_structural_dynamics() {
    let input = INPUT.replace("6000 0 0 0.35", "6000 0.00001 500 0.35");
    let source = ReedPerformance::from_bytes(input.as_bytes(), 37).unwrap();
    assert!(source.info().massive_reed);
    let expected = reference(48_000, 4801, &geometry(), 293.15, 1e-5, 500.0);
    assert_eq!(bits(&render(&input, 37)), bits(&expected));
    assert_ne!(bits(&expected), bits(&render(INPUT, 37)), "massive mechanics must not become a static aperture");
}

#[test]
fn supplied_segmented_bore_and_moist_air_owner_control_the_actual_trajectory() {
    let duct = Duct { segments: vec![
        Segment::Cylinder { radius: 0.0022, length: 0.2 },
        Segment::ToneHole { hole_radius: 0.0005, chimney_height: 0.001,
            bore_radius: 0.0022, state: HoleState::Vent(0.5) },
        Segment::Cone { inlet_radius: 0.0022, outlet_radius: 0.0023, length: 0.3 },
    ] };
    let input = INPUT.replace("segments 1\ncylinder 0.0022 0.5",
        "segments 3\ncylinder 0.0022 0.2\nhole 0.0005 0.001 0.0022 0.5\ncone 0.0022 0.0023 0.3")
        .replace("ambient 293.15", "ambient 303.15");
    let source = ReedPerformance::from_bytes(input.as_bytes(), 37).unwrap();
    assert_eq!(source.info().segments, 3); assert_eq!(source.info().tone_holes, 1);
    let expected = reference(48_000, 4801, &duct, 303.15, 0.0, 0.0);
    assert_eq!(bits(&render(&input, 37)), bits(&expected));
    assert_ne!(bits(&expected), bits(&render(INPUT, 37)));
}

#[test]
fn complete_window_and_callback_admission_precede_pressure_changes() {
    let mut source = ReedPerformance::from_bytes(INPUT.as_bytes(), 37).unwrap();
    assert_eq!(source.samples_rendered(), 0); assert!(source.renderer().applied_controls().is_empty());
    assert!(source.validate_sample_count(4802).is_err());
    let mut sentinel = [12345.0; 38];
    assert!(source.block(&mut sentinel).is_err()); assert_eq!(sentinel, [12345.0; 38]);
    assert!(source.renderer().applied_controls().is_empty());
    let mut out = vec![0.0; 4800];
    for chunk in out.chunks_mut(37) { source.block(chunk).unwrap(); }
    assert_eq!(source.remaining_samples(), 1);
    assert!(source.block(&mut sentinel[..2]).is_err()); assert_eq!(source.remaining_samples(), 1);
    source.block(&mut sentinel[..1]).unwrap();
    assert!(source.validate_sample_count(1).is_err());
    assert_eq!(source.renderer().pending_controls().len(), 0);
}

#[test]
fn rate_conversion_cancellation_and_retry_retain_the_finite_physical_window() {
    let input = INPUT.replace("audio 48000 4801", "audio 96000 9602");
    let pressure = reference(96_000, 9602, &geometry(), 293.15, 0.0, 0.0);
    let mut filter = Decimator::new(2, 1).unwrap();
    let expected: Vec<f64> = pressure.chunks_exact(2).map(|block| {
        let value = filter.preview(block).unwrap()[0]; filter.commit(); value
    }).collect();
    let (wav, clips) = encode_pcm16_wav(&expected, 48_000, 20000.0).unwrap();
    let finite = ReedPerformance::from_bytes(input.as_bytes(), 37).unwrap();
    let mut source = DecimatedRenderer::new(finite, 96_000, 48_000, 37).unwrap();
    assert!(source.validate_sample_count(4802).is_err());
    let mut stream = Pcm16WavStream::new(Cursor::new(Vec::new()), 48_000, 20000.0, 37).unwrap();
    let mut scratch = [0.0; 37];
    render_pressure_pcm16(&mut source, &mut stream, &CancelGate::new(), &mut scratch, 37).unwrap();
    let pending = source.source().renderer().pending_controls().to_vec();
    let before = stream.sink().get_ref().clone();
    let gate = CancelGate::new(); gate.request(); scratch.fill(12345.0);
    assert_eq!(render_pressure_pcm16(&mut source, &mut stream, &gate, &mut scratch, 4764).unwrap(),
        ScheduledWavProgress::Cancelled { samples: 0 });
    assert_eq!(scratch, [12345.0; 37]); assert_eq!(stream.sink().get_ref(), &before);
    assert_eq!(source.source().samples_rendered(), 74);
    assert_eq!(source.source().renderer().pending_controls(), pending.as_slice());
    render_pressure_pcm16(&mut source, &mut stream, &CancelGate::new(), &mut scratch, 4764).unwrap();
    let (output, summary) = stream.finish().unwrap();
    assert_eq!(output.into_inner(), wav); assert_eq!(summary.clipped_samples, clips as u64);
    assert_eq!(source.source().remaining_samples(), 0);
}

#[test]
fn malformed_or_unobserved_physics_and_resource_requests_refuse_before_realization() {
    for text in [
        INPUT.replace("segments 1", "segments 18446744073709551615"),
        INPUT.replace("6000 0 0 0.35", "6000 -1 0 0.35"),
        INPUT.replace("ambient 293.15", "ambient NaN"),
        INPUT.replace("ambient 293.15 101325 0", "ambient 293.15 101325 2"),
        INPUT.replace("cylinder 0.0022 0.5", "cylinder 0.0022 0"),
        INPUT.replace("cylinder 0.0022 0.5", "hole 0.0005 0.001 0.0022 1.1"),
        INPUT.replace("termination unflanged", "termination imaginary"),
        INPUT.replace("compile_limits 100000 1024", "compile_limits 0 1024"),
        INPUT.replace("compile_limits 100000 1024", "compile_limits 100000 1"),
        INPUT.replace("events\t5", "events\t18446744073709551615"),
        INPUT.replace("tracks\t1", "tracks\t18446744073709551615"),
        INPUT.replace("control_rate_hz\t700", "control_rate_hz\t0"),
        INPUT.replace("control_rate_hz\t700", "control_rate_hz\t96000"),
        INPUT.replace("event\t6e-2", "event\t2e-1"),
        INPUT.replace("listener_m 1", "listener_m 1 ignored"),
        format!("{INPUT}ignored\n"),
    ] { assert!(ReedPerformance::from_bytes(text.as_bytes(), 37).is_err(), "accepted {text}"); }
    assert!(ReedPerformance::from_bytes(INPUT.as_bytes(), 0).is_err());
    assert!(ReedPerformance::from_bytes(&vec![b'x'; MAX_REED_PERFORMANCE_BYTES + 1], 37).is_err());
}
