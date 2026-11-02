use super::*;

const EXAMPLE: &[u8] =
    include_bytes!("../../../../../../../../examples/coupled_continuous_tones.performance");

fn analytic_file() -> String {
    "frankensim-modal-performance-v2\nsample_rate_hz 48000\nsamples 128\nfull_scale_pa 1\nlimits 0.9 1 100 100 1000000\ncompile_limits 32 4096\nvoices 2\nvoice retain-state 1 1\nmode 2000 0 2 0 0 0\nport 0 -2\nvoice retain-state 1 1\nmode 2000 0 -1 0 0 0\nport 0 1.2\ncoupling_limits 1 100000 0.9 100 1000000 100000000 1e-11 1e-12 1e-8\nconnections 1\nconnection 7000000 0 0\nleft 0 1\nright 1 1\nevents 1\nforce 0 0 0 0.1\ncoupled_tone_limits 8 128 134217728 10000000000\ntones 2\ntone 0 0 60000 1.5 0.25\ntone 1 0 60000 -0.4 0.1\n".into()
}

#[test]
fn v2_example_keeps_network_preload_phase_events_and_frame_across_callbacks() {
    let mut whole = CoupledTonePerformance::from_bytes(EXAMPLE, 4800).unwrap();
    let mut split = CoupledTonePerformance::from_bytes(EXAMPLE, 37).unwrap();
    assert_eq!(whole.continuous_force_tones(), 3);
    assert_eq!(whole.compiled_controls(), 3);
    let k = 200000.0;
    let a = 1256.6370614359173_f64.powi(2) + k;
    let d = 1884.9555921538758_f64.powi(2) + k;
    let expected = [
        (d * 0.075 + k * 0.014) / (a * d - k * k),
        (k * 0.075 + a * 0.014) / (a * d - k * k),
    ];
    for (model, q) in whole.network().components().iter().zip(expected) {
        assert!((model.states()[0].displacement_m_sqrt_kg - q).abs() < 1e-18);
        assert_eq!(model.states()[0].velocity_m_sqrt_kg_per_s, 0.0);
    }
    let mut one = vec![0.0; 4800];
    let mut many = vec![0.0; 4800];
    whole.block(&mut one).unwrap();
    for chunk in many.chunks_mut(37) {
        split.block(chunk).unwrap();
    }
    assert_eq!(one, many);
    assert!(one.iter().any(|p| p.abs() > 1e-6));
    assert_eq!(whole.cosine, split.cosine);
    assert_eq!(whole.sine, split.sine);
    assert_eq!(whole.network().last_frame(), split.network().last_frame());
    assert_eq!(whole.network().samples_rendered(), 4800);
    assert_eq!(whole.samples_rendered(), 4800);
    assert_eq!(whole.held, vec![0.0, -0.7 * 0.04]);
    assert_eq!(whole.next_event, whole.events.len());
    let frame = whole.network().last_frame().unwrap();
    assert!(frame.component_dissipation_j > 0.0);
    assert!(frame.connection_dissipation_j > 0.0);
    assert_eq!(frame.numerical_dissipation_j, 0.0);
    assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
    assert_eq!(
        whole.info().input_hash,
        fs_blake3::hash_domain(MODAL_COUPLED_PERFORMANCE_HASH_DOMAIN, EXAMPLE)
    );
    assert!(whole.block(&mut [0.0]).is_err());
}

// Independent exact forced-oscillator solution. The two identical connected
// masses separate into symmetric/antisymmetric coordinates, so this tests the
// actual file projection and coupling against physics outside the dense owner.
fn motion(w: f64, omega: f64, fc: f64, fs: f64, dc: f64, t: f64) -> (f64, f64) {
    let den = w * w - omega * omega;
    let q = fc / den * ((omega * t).cos() - (w * t).cos())
        + fs / den * ((omega * t).sin() - omega / w * (w * t).sin())
        + dc / (w * w) * (1.0 - (w * t).cos());
    let v = fc / den * (-omega * (omega * t).sin() + w * (w * t).sin())
        + fs / den * omega * ((omega * t).cos() - (w * t).cos())
        + dc / w * (w * t).sin();
    (q, v)
}

#[test]
fn signed_tones_follow_continuous_coupled_solution_and_physical_work() {
    let text = analytic_file();
    let mut performance = CoupledTonePerformance::from_bytes(text.as_bytes(), 1).unwrap();
    let mut work = 0.0;
    for sample in 1..=128 {
        let mut pressure = [0.0];
        performance.block(&mut pressure).unwrap();
        let t = f64::from(sample) / 48000.0;
        let plus = motion(
            2000.0,
            60000.0,
            (-3.0 - 0.48) / 2.0,
            (-0.5 + 0.12) / 2.0,
            -0.1,
            t,
        );
        let minus = motion(
            18_000_000.0_f64.sqrt(),
            60000.0,
            (-3.0 + 0.48) / 2.0,
            (-0.5 - 0.12) / 2.0,
            -0.1,
            t,
        );
        let expected = [
            (plus.0 + minus.0, plus.1 + minus.1),
            (plus.0 - minus.0, plus.1 - minus.1),
        ];
        for (model, (q, v)) in performance.network().components().iter().zip(expected) {
            assert!((model.states()[0].displacement_m_sqrt_kg - q).abs() < 2e-16);
            assert!((model.states()[0].velocity_m_sqrt_kg_per_s - v).abs() < 2e-12);
        }
        assert!((pressure[0] - (plus.1 + 3.0 * minus.1)).abs() < 1e-11);
        let frame = performance.network().last_frame().unwrap();
        work += frame.external_work_j;
        assert_eq!(frame.total_dissipation_j(), 0.0);
        assert_eq!(frame.numerical_dissipation_j, 0.0);
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
    }
    assert!((work - performance.network().total_energy_j().unwrap()).abs() < 1e-13);
    // First-interval ZOH would substitute Fc+DC for the entire 1.25-rad tone
    // arc, producing measurably different velocity even on this fine clock.
    let dt: f64 = 1.0 / 48000.0;
    let zp = (-3.48 / 2.0 - 0.1) / 2000.0 * (2000.0 * dt).sin();
    let wm = 18_000_000.0_f64.sqrt();
    let zm = (-2.52 / 2.0 - 0.1) / wm * (wm * dt).sin();
    let p = motion(2000.0, 60000.0, -1.74, -0.19, -0.1, dt).1;
    let m = motion(wm, 60000.0, -1.26, -0.31, -0.1, dt).1;
    assert!(((zp + 3.0 * zm) / (p + 3.0 * m) - 1.0).abs() > 0.01);
}

#[test]
fn source_refusal_keeps_phase_event_cursor_and_complete_network_unchanged() {
    let text = analytic_file().replace("0.9 100 1000000 100000000", "0.9 100 1e-30 100000000");
    let mut performance = CoupledTonePerformance::from_bytes(text.as_bytes(), 1).unwrap();
    let phase = (performance.cosine.clone(), performance.sine.clone());
    let states: Vec<_> = performance
        .network()
        .components()
        .iter()
        .map(|m| m.states().to_vec())
        .collect();
    assert!(performance.block(&mut [0.0]).is_err());
    assert_eq!(performance.samples_rendered(), 0);
    assert_eq!(performance.network().samples_rendered(), 0);
    assert_eq!(performance.network().last_frame(), None);
    assert_eq!(performance.next_event, 0);
    assert_eq!(performance.held, vec![0.0, 0.0]);
    assert_eq!(
        (performance.cosine.clone(), performance.sine.clone()),
        phase
    );
    for (model, old) in performance.network().components().iter().zip(states) {
        assert_eq!(model.states(), old);
    }
    assert!(performance.validate_sample_rate(48000).is_err());
}

#[test]
fn coupled_file_rejects_unresolved_sources_budgets_and_unsupported_models() {
    let valid = analytic_file();
    for bad in [
        valid.replace("8 128 134217728 10000000000", "2 128 134217728 10000000000"),
        valid.replace("8 128 134217728 10000000000", "8 128 1 10000000000"),
        valid.replace("8 128 134217728 10000000000", "8 128 134217728 1"),
        valid.replace("8 128 134217728 10000000000", "8 1 134217728 10000000000"),
        valid.replace("compile_limits 32 4096", "compile_limits 32 2"),
        valid.replace("tone 1 0 60000", "tone 1 9 60000"),
        valid.replace("tone 1 0 60000", "tone 1 0 140000"),
        valid.replace("tone 1 0 60000", "tone 1 0 NaN"),
        valid.replace("mode 2000 0 2 0", "mode 2000 0 2 0.1"),
        valid.replacen("voice retain-state", "voice static-preload", 1),
        valid.replace(
            "frankensim-modal-performance-v2",
            "frankensim-modal-performance-v3",
        ),
        format!("{valid}broadband 2\n"),
    ] {
        assert!(
            CoupledTonePerformance::from_bytes(bad.as_bytes(), 128).is_err(),
            "admitted {bad}"
        );
    }
    let performance = CoupledTonePerformance::from_bytes(valid.as_bytes(), 16).unwrap();
    assert!(performance.validate_sample_rate(44100).is_err());
    assert!(performance.validate_sample_count(129).is_err());
}
