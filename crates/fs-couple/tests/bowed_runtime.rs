//! Stateful bow performance: real modal physics, no synthetic signal generator.
use fs_couple::bowed_string::{
    BowGesture, BowedRunConfig, BowedRunError, BowedStringCard, FrictionIsland,
    Termination, run_bowed,
};
use fs_couple::bowed_string::runtime::{BowedRenderOutcome, BowedSample, BowedStringState};
use fs_couple::stribeck_friction::StribeckFriction;
use fs_exec::CancelGate;

fn config() -> BowedRunConfig {
    BowedRunConfig {
        card: BowedStringCard {
            length_m: 0.65,
            tension_n: 60.0,
            linear_density_kg_m: 6.0e-4,
            bending_stiffness_n_m2: 0.0,
            viscous_bending_n_m2_s: 0.0,
            mode_count: 16,
            zetas: (0..16).map(|k| 1.0e-3 * (1.0 + 0.55 * f64::from(k))).collect(),
            sample_rate_hz: 48_000,
        },
        island: FrictionIsland::Stribeck(StribeckFriction::try_new(0.8, 0.4, 0.04).unwrap()),
        gesture: BowGesture::admit(0.45, 3.9, 0.11).unwrap(),
        steps: 1024,
        subsamples: 16,
        termination: Termination::Rigid,
        listener_m: 1.0,
    }
}

fn mechanical_bits(sample: BowedSample) -> [u64; 4] {
    [sample.bow_point_velocity_m_s.to_bits(), sample.relative_velocity_m_s.to_bits(),
        sample.bridge_force_n.to_bits(), sample.total_modal_energy_j.to_bits()]
}

#[test]
fn batch_and_block_collectors_preserve_every_endpoint() {
    let cfg = config();
    let batch = run_bowed(&cfg).unwrap();
    assert!(batch.bridge_force_n.iter().any(|v| v.abs() > 1e-8));
    for partition in [1, 37, 256, 1024] {
        let mut state = BowedStringState::new(&cfg, 1024).unwrap();
        let mut samples = vec![BowedSample::default(); cfg.steps];
        for block in samples.chunks_mut(partition) { state.block(block).unwrap(); }
        for (i, sample) in samples.iter().enumerate() {
            assert_eq!(sample.bow_point_velocity_m_s.to_bits(), batch.bow_point_velocity_m_s[i].to_bits());
            assert_eq!(sample.relative_velocity_m_s.to_bits(), batch.relative_velocity_m_s[i].to_bits());
            assert_eq!(sample.bridge_force_n.to_bits(), batch.bridge_force_n[i].to_bits());
            assert_eq!(sample.radiated_pressure_pa, None);
        }
        assert_eq!(state.samples_rendered(), cfg.steps as u64);
        assert_eq!(state.total_modal_energy_j().to_bits(), batch.final_total_energy_j.to_bits());
        assert_eq!(state.peak_modal_energy_j().to_bits(), batch.peak_total_energy_j.to_bits());
    }
}

#[test]
fn first_sample_matches_independent_oscillator_solution() {
    let mut cfg = config();
    cfg.card.mode_count = 1;
    cfg.card.zetas = vec![0.01];
    cfg.subsamples = 1;
    cfg.gesture = BowGesture::admit(0.2, 1.0, 0.2).unwrap();
    cfg.island = FrictionIsland::ViscousOnly { viscous_n_s_per_m: 0.01 };
    let kappa = core::f64::consts::PI / cfg.card.length_m;
    let omega = kappa * (cfg.card.tension_n / cfg.card.linear_density_kg_m).sqrt();
    let decay = 0.01 * omega;
    let wd = (omega * omega - decay * decay).sqrt();
    let norm = (cfg.card.linear_density_kg_m * cfg.card.length_m * 0.5).sqrt();
    let phi = (core::f64::consts::PI * 0.2).sin() / norm;
    let force = 0.002 * phi;
    let dt = 1.0 / 48_000.0;
    let q = force / (omega * omega)
        * (1.0 - (-decay * dt).exp() * ((wd * dt).cos() + decay / wd * (wd * dt).sin()));
    let v = force * (-decay * dt).exp() * (wd * dt).sin() / wd;
    let sample = BowedStringState::new(&cfg, 1).unwrap().step().unwrap();
    assert!((sample.bow_point_velocity_m_s - phi * v).abs() < 1e-12);
    assert!((sample.bridge_force_n - cfg.card.tension_n * kappa / norm * q).abs() < 1e-12);
    assert!((sample.total_modal_energy_j - 0.5 * (v * v + omega * omega * q * q)).abs() < 1e-16);
}

#[test]
fn lifting_bow_preserves_ringdown_and_removes_even_viscous_contact() {
    let mut cfg = config();
    cfg.island = FrictionIsland::ViscousOnly { viscous_n_s_per_m: 0.01 };
    let mut a = BowedStringState::new(&cfg, 512).unwrap();
    let mut b = BowedStringState::new(&cfg, 512).unwrap();
    for _ in 0..512 { a.step().unwrap(); b.step().unwrap(); }
    let energy = a.total_modal_energy_j();
    assert!(energy > 1e-10);
    a.set_bow(9.0, 0.0, 0.11).unwrap();
    b.set_bow(-9.0, 0.0, 0.11).unwrap();
    assert_eq!(a.total_modal_energy_j().to_bits(), energy.to_bits());
    assert_eq!(a.samples_rendered(), 512);
    let mut previous = energy;
    let mut audible_motion = false;
    for _ in 0..512 {
        let x = a.step().unwrap();
        let y = b.step().unwrap();
        assert_eq!(x.bridge_force_n.to_bits(), y.bridge_force_n.to_bits());
        assert_eq!(x.bow_point_velocity_m_s.to_bits(), y.bow_point_velocity_m_s.to_bits());
        assert_eq!(x.total_modal_energy_j.to_bits(), y.total_modal_energy_j.to_bits());
        assert!(x.total_modal_energy_j <= previous + 1e-12 * energy);
        previous = x.total_modal_energy_j;
        audible_motion |= x.bow_point_velocity_m_s.abs() > 1e-8;
    }
    assert!(audible_motion, "release must not silence the string");
}

#[test]
fn bow_reversal_station_change_and_reentry_do_not_reset_physics() {
    let render = |partition: usize| {
        let mut state = BowedStringState::new(&config(), 512).unwrap();
        let mut samples = Vec::new();
        for (v, force, station) in [(0.45, 3.9, 0.11), (0.0, 0.0, 0.11), (-0.3, 2.0, 0.2)] {
            let before = state.total_modal_energy_j().to_bits();
            state.set_bow(v, force, station).unwrap();
            assert_eq!(before, state.total_modal_energy_j().to_bits());
            let mut phase = vec![BowedSample::default(); 512];
            for block in phase.chunks_mut(partition) { state.block(block).unwrap(); }
            samples.extend(phase.into_iter().map(mechanical_bits));
        }
        samples
    };
    assert_eq!(render(512), render(37));
}

#[test]
fn cancelled_state_resumes_without_restarting_contact() {
    let cfg = config();
    let mut state = BowedStringState::new(&cfg, 37).unwrap();
    let mut prefix = vec![BowedSample::default(); 83];
    assert_eq!(state.render_under_gate(&CancelGate::new(), &mut prefix, 37).unwrap(),
        BowedRenderOutcome::Completed { samples: 83 });
    let before = state.total_modal_energy_j().to_bits();
    let sentinel = BowedSample { bridge_force_n: 12345.0, ..BowedSample::default() };
    let mut suffix = vec![sentinel; 119];
    let cancelled = CancelGate::new();
    cancelled.request();
    assert_eq!(state.render_under_gate(&cancelled, &mut suffix, 37).unwrap(),
        BowedRenderOutcome::Cancelled { samples: 0 });
    assert_eq!(suffix, vec![sentinel; 119]);
    assert_eq!(state.samples_rendered(), 83);
    assert_eq!(state.total_modal_energy_j().to_bits(), before);
    state.render_under_gate(&CancelGate::new(), &mut suffix, 37).unwrap();
    let mut direct = BowedStringState::new(&cfg, 202).unwrap();
    let mut expected = vec![BowedSample::default(); 202];
    direct.block(&mut expected).unwrap();
    prefix.extend(suffix);
    assert_eq!(prefix.into_iter().map(mechanical_bits).collect::<Vec<_>>(),
        expected.into_iter().map(mechanical_bits).collect::<Vec<_>>());
}

#[test]
fn invalid_controls_and_blocks_are_transactional() {
    let mut state = BowedStringState::new(&config(), 8).unwrap();
    let mut direct = BowedStringState::new(&config(), 8).unwrap();
    state.step().unwrap(); direct.step().unwrap();
    for (v, force, station) in [(f64::NAN, 1.0, 0.1), (0.1, -1.0, 0.1),
        (0.1, f64::INFINITY, 0.1), (0.1, 1.0, 1.0), (0.1, 1.0, 0.0)] {
        assert!(state.set_bow(v, force, station).is_err());
    }
    assert!(state.block(&mut []).is_err());
    let sentinel = BowedSample { bridge_force_n: 123.0, ..BowedSample::default() };
    let mut out = [sentinel; 9];
    assert!(state.block(&mut out).is_err());
    assert_eq!(out, [sentinel; 9]);
    let mut pressure = [123.0; 8];
    assert!(state.pressure_block(&mut pressure).is_err());
    assert_eq!(pressure, [123.0; 8]);
    for _ in 0..8 { assert_eq!(mechanical_bits(state.step().unwrap()), mechanical_bits(direct.step().unwrap())); }
}

#[test]
fn failed_substep_poisons_voice_and_prevents_partial_resume() {
    let mut cfg = config();
    cfg.island = FrictionIsland::ViscousOnly { viscous_n_s_per_m: f64::MAX };
    let mut state = BowedStringState::new(&cfg, 8).unwrap();
    assert!(state.step().is_err());
    assert!(matches!(state.step(), Err(BowedRunError::Poisoned)));
    assert!(matches!(state.set_bow(0.0, 0.0, 0.11), Err(BowedRunError::Poisoned)));
    let mut untouched = [BowedSample::default(); 8];
    assert!(matches!(state.block(&mut untouched), Err(BowedRunError::Poisoned)));
    assert_eq!(untouched, [BowedSample::default(); 8]);
}
