//! G3 source phase, retained coupled state, and accepted-clock file integration.
use super::*;

const EXAMPLE: &[u8] =
    include_bytes!("../../../../../../../../examples/coupled_swept_tones.performance");

fn integrated_phase(knots: &[(u64, f64)], rate: f64) -> f64 {
    knots
        .windows(2)
        .map(|p| 0.5 * (p[0].1 + p[1].1) * (p[1].0 - p[0].0) as f64 / rate)
        .sum()
}

#[test]
fn coupled_sweeps_preserve_integrated_phase_through_silence_and_callback_partitions() {
    let mut whole = CoupledTonePerformance::from_bytes(EXAMPLE, 4800).unwrap();
    let mut split = CoupledTonePerformance::from_bytes(EXAMPLE, 37).unwrap();
    assert_eq!(whole.frequency_sweep_count(), 3);
    assert_eq!(whole.envelope_count(), 3);
    assert!(whole.maximum_sweep_force_error_bound_n() > 0.0);
    assert!(whole.maximum_sweep_force_error_bound_n() <= 1e-4);
    let mut a = vec![0.0; 4800];
    let mut b = vec![0.0; 4800];
    whole.block(&mut a).unwrap();
    for block in b.chunks_mut(37) {
        split.block(block).unwrap();
    }
    assert_eq!(a, b);
    assert_eq!(whole.cosine, split.cosine);
    assert_eq!(whole.sine, split.sine);
    assert_eq!(whole.network().last_frame(), split.network().last_frame());
    for (a, b) in whole
        .network()
        .components()
        .iter()
        .zip(split.network().components())
    {
        assert_eq!(a.states(), b.states());
    }
    assert!(a.iter().any(|p| p.abs() > 1e-7));
    let phases = [
        integrated_phase(
            &[
                (0, 3769.9111843077517),
                (960, 5026.548245743669),
                (1440, 4084.070449666731),
                (3000, 4398.22971502571),
                (4800, 4398.22971502571),
            ],
            48000.0,
        ),
        integrated_phase(
            &[
                (0, 4712.38898038469),
                (1800, 4084.070449666731),
                (3200, 4398.22971502571),
                (4800, 4398.22971502571),
            ],
            48000.0,
        ),
        integrated_phase(
            &[
                (0, 4411.424404170788),
                (1800, 4712.38898038469),
                (2800, 4084.070449666731),
                (4800, 4411.424404170788),
            ],
            48000.0,
        ),
    ];
    // Two rows share a prepared carrier but have distinct physical sweeps;
    // original row identity and signed physical port projection must survive.
    for (tone, mode, c, s) in [
        (0, 0, -0.15, -0.05),
        (1, 1, 0.14, -0.035),
        (2, 0, 0.1, -0.04),
    ] {
        let (sn, cs) = phases[tone].sin_cos();
        let i = 4 * mode + tone + 1;
        assert!((whole.cosine[i] - (c * cs + s * sn)).abs() < 2e-12);
        assert!((whole.sine[i] - (s * cs - c * sn)).abs() < 2e-12);
    }
    assert_eq!(
        whole.info().input_hash,
        fs_blake3::hash_domain(MODAL_COUPLED_PERFORMANCE_HASH_DOMAIN, EXAMPLE)
    );
}

#[test]
fn swept_release_keeps_coupled_and_receiver_ringdown_without_source_work() {
    let text = std::str::from_utf8(EXAMPLE)
        .unwrap()
        .replace("force 613 1 0 0.04", "force 613 1 0 0");
    let mut source = CoupledTonePerformance::from_bytes(text.as_bytes(), 3200).unwrap();
    source.block(&mut vec![0.0; 3200]).unwrap();
    let before = source.network().total_energy_j().unwrap();
    let mut tail = vec![0.0; 1600];
    source.block(&mut tail).unwrap();
    let after = source.network().total_energy_j().unwrap();
    assert!(before > after && after > 0.0);
    assert!(tail.iter().any(|p| p.abs() > 1e-9));
    let frame = source.network().last_frame().unwrap();
    assert_eq!(frame.external_work_j, 0.0);
    assert_eq!(frame.numerical_dissipation_j, 0.0);
    assert!(frame.total_dissipation_j() > 0.0);
}

#[test]
fn frequency_admission_and_failed_sample_preserve_all_accepted_source_state() {
    let text = std::str::from_utf8(EXAMPLE).unwrap();
    for invalid in [
        text.replace("sweep 0 5 0.0001", "sweep 0 5 1e-30"),
        text.replace("omega 960 5026.548245743669", "omega 960 140000"),
        text.replace("sweep 1 4 0.0001", "sweep 0 4 0.0001"),
        text.replace(
            "omega 4800 4411.424404170788",
            "omega 4801 4411.424404170788",
        ),
    ] {
        assert!(CoupledTonePerformance::from_bytes(invalid.as_bytes(), 240).is_err());
    }
    let text = text.replace("events 3\n", "events 4\nforce 240 0 0 1e20\n");
    let mut source = CoupledTonePerformance::from_bytes(text.as_bytes(), 240).unwrap();
    source.block(&mut vec![0.0; 240]).unwrap();
    let phase = (source.cosine.clone(), source.sine.clone());
    let held = source.held.clone();
    let event = source.next_event;
    let frame = source.network().last_frame().cloned();
    let states: Vec<_> = source
        .network()
        .components()
        .iter()
        .map(|model| model.states().to_vec())
        .collect();
    let next_force = |source: &CoupledTonePerformance| {
        let gain = source
            .envelopes
            .as_ref()
            .unwrap()
            .endpoints(0, 240)
            .unwrap();
        source
            .frequency_sweeps
            .as_ref()
            .unwrap()
            .step(0, 240, gain)
            .unwrap()
            .unwrap()
            .endpoints(1.0, 0.0)
            .unwrap()
    };
    let expected = next_force(&source);
    // Exercise the candidate operation directly so the second attempt also
    // proves retry stability beneath the file renderer's fatal-error policy.
    for _ in 0..2 {
        assert!(source.step().is_err());
        assert_eq!(source.completed, 240);
        assert_eq!(source.network().samples_rendered(), 240);
        assert_eq!(source.next_event, event);
        assert_eq!(source.held, held);
        assert_eq!((source.cosine.clone(), source.sine.clone()), phase);
        assert_eq!(next_force(&source), expected);
        assert_eq!(source.network().last_frame(), frame.as_ref());
        for (model, state) in source.network().components().iter().zip(&states) {
            assert_eq!(model.states(), state);
        }
    }
}

#[test]
fn empty_sweep_section_retains_existing_envelope_path_exactly() {
    let text = include_str!("../../../../../../../../examples/coupled_enveloped_tones.performance");
    let explicit = text.replace(
        "coupled_broadband ",
        "frequency_sweeps 0\ncoupled_broadband ",
    );
    let mut original = CoupledTonePerformance::from_bytes(text.as_bytes(), 4800).unwrap();
    let mut empty = CoupledTonePerformance::from_bytes(explicit.as_bytes(), 4800).unwrap();
    let mut a = vec![0.0; 4800];
    let mut b = vec![0.0; 4800];
    original.block(&mut a).unwrap();
    empty.block(&mut b).unwrap();
    assert_eq!(a, b);
    assert_eq!(empty.frequency_sweep_count(), 0);
    assert_eq!(empty.maximum_sweep_force_error_bound_n(), 0.0);
    assert_eq!(
        original.network().last_frame(),
        empty.network().last_frame()
    );
}
