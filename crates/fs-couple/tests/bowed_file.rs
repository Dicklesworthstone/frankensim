//! Authored bow mechanics and gestures through the existing solver and real CLI.
use fs_couple::bowed_string::{BowGesture, BowedRunConfig, BowedStringCard, FrictionIsland, Termination};
use fs_couple::bowed_string::runtime::BowedStringState;
use fs_couple::bowed_string::runtime::schedule::file::BowedPerformance;
use fs_couple::pcm_wav::{decimate::Decimator, encode_pcm16_wav};
use fs_couple::render::plate::file::PlatePerformance;
use fs_couple::stribeck_friction::StribeckFriction;
use fs_couple::thin_plate::CompactBody;
use fs_material::gas::GasState;
use fs_scenario::{RadiatingPlate, gesture::{GestureEvent, GestureSchedule, GestureTarget, GestureTrack, GestureValue}};

const INPUT: &str = include_str!("../examples/bowed-string.performance");
const PLATE: &str = include_str!("../examples/plate-mesh.performance");
fn schedule() -> GestureSchedule {
    let value = |v, f, s| GestureValue::Bow { velocity_m_per_s: v, normal_force_n: f, station: s };
    GestureSchedule::try_new(700, vec![GestureTrack { id: "bow".into(), target: GestureTarget::BowStroke { string: 0 },
        initial: value(0.45, 3.9, 0.11), events: vec![
            GestureEvent { time_s: 0.0, transition_s: 0.01, value: value(-0.3, 2.0, 0.2) },
            GestureEvent { time_s: 0.005, transition_s: 0.01, value: value(0.25, 1.0, 0.15) },
            GestureEvent { time_s: 0.02, transition_s: 0.0, value: value(0.0, 0.0, 0.15) },
            GestureEvent { time_s: 0.04, transition_s: 0.0, value: value(-0.3, 2.0, 0.12) },
            GestureEvent { time_s: 0.06, transition_s: 0.0, value: value(0.0, 0.0, 0.12) },
        ] }]).unwrap()
}
fn reference(rate: u32, samples: usize, tension: f64) -> (Vec<f64>, f64) {
    let config = BowedRunConfig {
        card: BowedStringCard { length_m: 0.65, tension_n: tension, linear_density_kg_m: 0.0006,
            bending_stiffness_n_m2: 0.0, viscous_bending_n_m2_s: 0.0, mode_count: 16,
            zetas: vec![0.001; 16], sample_rate_hz: rate },
        island: FrictionIsland::Stribeck(StribeckFriction::try_new(0.8, 0.4, 0.04).unwrap()),
        gesture: BowGesture::admit(0.45, 3.9, 0.11).unwrap(), steps: samples, subsamples: 16,
        termination: Termination::PlateOnePort {
            body: Box::new(CompactBody::from_radiator(RadiatingPlate {
                area_m2: 0.003, mass_kg: 0.15, frequency_hz: 280.0, damping_ratio: 0.02,
            }).unwrap()), ambient: GasState::try_new_moist_air(293.15, 101_325.0, 0.0).unwrap(),
        }, listener_m: 1.0,
    };
    let mut state = BowedStringState::new(&config, 1).unwrap();
    let schedule = schedule();
    let mut pressure = Vec::new();
    let mut previous = None;
    for sample in 0..samples as u64 {
        let GestureValue::Bow { velocity_m_per_s: v, normal_force_n: f, station: s } =
            schedule.sample_value("bow", sample * 700 / u64::from(rate)).unwrap() else { panic!("bow value") };
        let key = [v.to_bits(), f.to_bits(), s.to_bits()];
        if previous != Some(key) { state.set_bow(v, f, s).unwrap(); previous = Some(key); }
        pressure.push(0.0 + state.step().unwrap().radiated_pressure_pa.unwrap());
    }
    (pressure, state.total_modal_energy_j())
}
fn bits(values: &[f64]) -> Vec<u64> { values.iter().map(|v| v.to_bits()).collect() }

#[test]
fn canonical_file_matches_physical_solver_and_retains_release_ringdown_across_callbacks() {
    assert_eq!(INPUT.split_once("\nschedule\n").unwrap().1.as_bytes(), schedule().to_canonical_bytes());
    let (expected, energy) = reference(48_000, 4801, 60.0);
    assert!(expected[3000..].iter().any(|p| p.abs() > 1e-12), "released string/body must keep ringing");
    for partition in [1, 37, 512] {
        let performance = BowedPerformance::from_bytes(INPUT.as_bytes(), partition).unwrap();
        assert_eq!(performance.info().gesture_events, 5);
        assert!(performance.info().compiled_controls > 5);
        let mut renderer = performance.into_renderer();
        assert_eq!(renderer.samples_rendered(), 0);
        let mut actual = vec![0.0; 4801];
        for block in actual.chunks_mut(partition) { renderer.block(block).unwrap(); }
        assert_eq!(bits(&actual), bits(&expected));
        let bow = renderer.context().bowed_performance(0).unwrap();
        assert_eq!(bow.remaining_samples(), 0);
        assert_eq!(bow.state().total_modal_energy_j().to_bits(), energy.to_bits());
        assert!(bow.applied_controls().iter().any(|event| event.velocity_m_s < 0.0));
        assert_eq!(bow.applied_controls().last().unwrap().normal_force_n, 0.0);
        assert!(renderer.validate_sample_count(1).is_err());
    }
}

#[test]
fn malformed_physics_unknown_tracks_unobserved_commands_and_budgets_refuse() {
    for text in [
        INPUT.replace("0 0 16\n", "0 0 18446744073709551615\n"),
        INPUT.replace("stribeck 0.8 0.4", "stribeck 0.8 0.9"),
        INPUT.replace("subsamples 16", "subsamples 0"),
        INPUT.replace("body 0.003", "body 0"),
        INPUT.replace("ambient 293.15", "ambient NaN"),
        INPUT.replace("compile_limits 100000 1024", "compile_limits 0 1024"),
        INPUT.replace("compile_limits 100000 1024", "compile_limits 100000 1"),
        INPUT.replace("track\tbow\tbow\t0", "track\tbow\tbow\t1"),
        INPUT.replace("events\t5", "events\t18446744073709551615"),
        INPUT.replace("tracks\t1", "tracks\t18446744073709551615"),
        INPUT.replace("event\t6e-2", "event\t2e-1"),
        INPUT.replace("listener_m 1", "listener_m 1 ignored"),
        format!("{INPUT}ignored\n"),
    ] { assert!(BowedPerformance::from_bytes(text.as_bytes(), 37).is_err(), "accepted {text}"); }
    assert!(BowedPerformance::from_bytes(INPUT.as_bytes(), 0).is_err());
    assert!(BowedPerformance::from_bytes(&vec![b'x'; 1024 * 1024 + 1], 37).is_err());
}

#[test]
fn bow_and_mesh_files_render_together_with_native_clock_gestures_and_tension_changes() {
    let stamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let dir = std::env::temp_dir().join(format!("fs-bowed-file-{}-{stamp}", std::process::id()));
    std::fs::create_dir(&dir).unwrap();
    let plate_path = dir.join("plate.performance"); std::fs::write(&plate_path, PLATE).unwrap();
    let mut plate = PlatePerformance::from_bytes(PLATE.as_bytes(), 37).unwrap().into_renderer();
    let mut plate_pressure = vec![0.0; 4801];
    for block in plate_pressure.chunks_mut(37) { plate.block(block).unwrap(); }
    let mut outputs = Vec::new();
    for (index, (ratio, tension, block)) in [(1, 60.0, 37), (1, 66.0, 512), (2, 60.0, 1)].into_iter().enumerate() {
        let rate = 48_000 * ratio;
        let samples = 4801 * ratio as usize;
        let text = INPUT.replace("audio 48000 4801 1", &format!("audio {rate} {samples} 1"))
            .replace("string 0.65 60 ", &format!("string 0.65 {tension} "));
        let bow_path = dir.join(format!("bow-{index}.performance")); std::fs::write(&bow_path, text).unwrap();
        let out = dir.join(format!("mixed-{index}.wav"));
        let mut command = std::process::Command::new(env!("CARGO_BIN_EXE_music_render"));
        command.arg("ensemble").arg(&out).arg("--bow").arg(&bow_path).arg("--plate").arg(&plate_path)
            .args(["--full-scale-pa", "1", "--block", &block.to_string()]);
        if ratio != 1 { command.arg("--decimate"); }
        let result = command.output().unwrap();
        assert!(result.status.success(), "{}", String::from_utf8_lossy(&result.stdout));
        let mut filter = Decimator::new(ratio as usize, 1).unwrap();
        let delay = filter.delay_output_frames() as usize;
        let mut expected: Vec<_> = reference(rate, samples, tension).0.chunks_exact(ratio as usize).map(|group| {
            let p = filter.preview(group).unwrap()[0]; filter.commit(); p
        }).collect();
        for (p, q) in expected[delay..].iter_mut().zip(&plate_pressure) { *p += q; }
        let actual = std::fs::read(&out).unwrap();
        assert_eq!(actual, encode_pcm16_wav(&expected, 48_000, 1.0).unwrap().0);
        let metadata = std::fs::read_to_string(out.with_extension("provenance.json")).unwrap();
        assert!(metadata.contains("\"kind\":\"bow\""));
        assert!(metadata.contains("frankensim-bowed-performance-v1"));
        assert!(metadata.contains(&format!("\"common_delay_output_samples\":{delay}")));
        outputs.push(actual);
    }
    assert_ne!(&outputs[0][44..], &outputs[1][44..], "string tension must affect the actual mixed sound");
}
