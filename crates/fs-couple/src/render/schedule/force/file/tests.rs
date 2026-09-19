use super::*;

const INPUT: &str = "frankensim-modal-performance-v1\n\
sample_rate_hz 48000\nsamples 257\nfull_scale_pa 0.02\n\
limits 0.9 1 100000 1e9 1e7\ncompile_limits 32 100\nvoices 1\n\
voice retain-state 2 2\nmode 1000 0.03 1 0 0 0\nmode 2000 0.04 0.75 0 0 0\n\
port 0.5 1 -0.5\nport -0.25 0.25 2\nevents 3\n\
force 71 0 1 0.5\nforce 13 0 0 0\nforce 120 0 1 0\n";

fn render(input: &str, block: usize) -> Vec<f64> {
    let parsed = ModalPerformance::from_bytes(input.as_bytes(), block).unwrap();
    let mut pressure = vec![0.0; parsed.info().samples as usize];
    let mut renderer = parsed.into_renderer();
    for chunk in pressure.chunks_mut(block) { renderer.block(chunk).unwrap(); }
    pressure
}

#[test]
fn imported_physical_forces_match_direct_modes_at_every_callback_partition() {
    let modes = [(1000.0, 0.03, 1.0), (2000.0, 0.04, 0.75)].into_iter()
        .map(|(omega, damping, transfer)| ModalAcousticMode {
            angular_frequency_rad_s: omega, damping_ratio: damping,
            pressure_per_modal_velocity: C64::new(transfer, 0.0),
        }).collect();
    let mut model = ModalAcousticTimeModel::try_new(48000, modes,
        ModalAcousticTimeBudget::audible_reference()).unwrap();
    let expected: Vec<_> = (0..257).map(|sample| {
        let a = if sample < 13 { 0.5 } else { 0.0 };
        let b = if sample < 71 { -0.25 } else if sample < 120 { 0.5 } else { 0.0 };
        model.step(&[a + 0.25*b, -0.5*a + 2.0*b]).unwrap().observer_pressure_pa.to_bits()
    }).collect();
    assert!(expected.iter().any(|&p| p != 0.0_f64.to_bits()));
    for block in [1, 7, 37, 64, 257] {
        assert_eq!(render(INPUT, block).into_iter().map(f64::to_bits).collect::<Vec<_>>(), expected);
    }
}

#[test]
fn preload_release_and_supplied_history_both_reach_the_live_runtime() {
    let preload = INPUT.replace("retain-state", "static-preload");
    let samples = render(&preload, 37);
    assert!(samples[..13].iter().all(|p| p.abs() < 1e-14));
    assert!(samples[13..].iter().any(|p| p.abs() > 1e-6));
    let history = INPUT.replace("mode 1000 0.03 1 0 0 0", "mode 1000 0.03 1 0 0.0001 0.2");
    assert_ne!(render(&history, 37), render(INPUT, 37));
    assert!(ModalPerformance::from_bytes(history.replace("retain-state", "static-preload").as_bytes(), 37).is_err());
}

#[test]
fn counts_fields_clocks_and_whole_input_are_admitted_before_rendering() {
    for (from, to) in [
        ("sample_rate_hz 48000", "sample_rate_hz 0"),
        ("samples 257", "samples 0"),
        ("samples 257", "samples 28800001"),
        ("full_scale_pa 0.02", "full_scale_pa NaN"),
        ("voices 1", "voices 18446744073709551615"),
        ("voice retain-state 2 2", "voice retain-state 4097 2"),
        ("voice retain-state 2 2", "voice retain-state 2 65536"),
        ("events 3", "events 65537"),
        ("force 71 0 1 0.5", "force 257 0 1 0.5"),
        ("force 71 0 1 0.5", "force 71 9 1 0.5"),
        ("force 71 0 1 0.5", "force 71 0 9 0.5"),
        ("port 0.5 1 -0.5", "port 0.5 1 -0.5 ignored"),
        ("port 0.5 1 -0.5", "port 0.5 1"),
        ("mode 1000 0.03 1 0 0 0", "mode 1000 -1 1 0 0 0"),
        ("mode 1000 0.03 1 0 0 0", "mode 1000000 0.03 1 0 0 0"),
        ("mode 1000 0.03 1 0 0 0", "mode 1000 0.03 1 0 NaN 0"),
        ("limits 0.9 1", "limits 2 1"),
        ("compile_limits 32 100", "compile_limits 0 100"),
        ("compile_limits 32 100", "compile_limits 32 1"),
    ] {
        assert!(ModalPerformance::from_bytes(INPUT.replace(from, to).as_bytes(), 64).is_err(), "{to}");
    }
    for block in [0, 65537] { assert!(ModalPerformance::from_bytes(INPUT.as_bytes(), block).is_err()); }
    assert!(ModalPerformance::from_bytes(format!("{INPUT}unknown 1\n").as_bytes(), 64).is_err());
    assert!(ModalPerformance::from_bytes(b"\xff", 64).is_err());
    assert!(ModalPerformance::from_bytes(&vec![b'x'; MAX_MODAL_PERFORMANCE_BYTES+1], 64).is_err());
}

#[test]
fn every_truncated_record_refuses_and_errors_name_the_line() {
    for (index, _) in INPUT.match_indices('\n') {
        if index+1 < INPUT.len() {
            assert!(ModalPerformance::from_bytes(&INPUT.as_bytes()[..=index], 64).is_err());
        }
    }
    assert!(matches!(ModalPerformance::from_bytes(b"wrong\n", 64),
        Err(ModalPerformanceError::Input { line: 1, .. })));
    assert!(matches!(ModalPerformance::from_bytes(INPUT.replace("samples 257", "samples 257 extra").as_bytes(), 64),
        Err(ModalPerformanceError::Input { line: 3, .. })));
}

#[test]
fn complete_source_bytes_bind_parameters_but_callback_size_does_not() {
    let a = ModalPerformance::from_bytes(INPUT.as_bytes(), 37).unwrap().info();
    let b = ModalPerformance::from_bytes(INPUT.as_bytes(), 64).unwrap().info();
    assert_eq!(a, b);
    assert_eq!(a.input_hash, hash_domain(MODAL_PERFORMANCE_HASH_DOMAIN, INPUT.as_bytes()));
    let changed = INPUT.replace("mode 2000", "mode 2500");
    assert_ne!(ModalPerformance::from_bytes(changed.as_bytes(), 37).unwrap().info().input_hash, a.input_hash);
    assert_ne!(render(&changed, 37), render(INPUT, 37));
    assert_eq!((a.voices, a.modes, a.force_events), (1, 2, 3));
}

#[test]
fn multiple_imported_voices_sum_without_fallback_to_a_fixture() {
    let body = "voice retain-state 2 2\nmode 1000 0.03 1 0 0 0\nmode 2000 0.04 0.75 0 0 0\nport 0.5 1 -0.5\nport -0.25 0.25 2\n";
    let one = INPUT.replace("events 3\nforce 71 0 1 0.5\nforce 13 0 0 0\nforce 120 0 1 0\n", "events 0\n");
    let two = one.replace("voices 1", "voices 2").replace(body, &format!("{body}{body}"));
    let expected = render(&one, 37);
    let actual = render(&two, 64);
    for (a, b) in actual.into_iter().zip(expected) { assert_eq!(a.to_bits(), (b+b).to_bits()); }
    let over_modes = two.replacen("voice retain-state 2 2", "voice retain-state 4096 1", 1);
    // Truncated modal rows, and aggregate limits, cannot create a partial runtime.
    assert!(ModalPerformance::from_bytes(over_modes.as_bytes(), 64).is_err());
}
