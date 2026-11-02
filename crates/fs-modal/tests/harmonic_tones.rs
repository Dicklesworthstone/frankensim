//! Independent physical power and interference checks for simultaneous tones.

use fs_math::c64::C64;
use fs_modal::time::{
    MAX_HARMONIC_TONES, ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget,
    ModalAcousticTimeModel, ModalAcousticWorkspace, ModalToneBudget, PreparedModalHarmonic,
    PreparedModalTones,
};

fn state(q: f64, v: f64) -> ModalAcousticState {
    ModalAcousticState {
        displacement_m_sqrt_kg: q,
        velocity_m_sqrt_kg_per_s: v,
    }
}
fn model(modes: &[(f64, f64)]) -> ModalAcousticTimeModel {
    ModalAcousticTimeModel::try_new_with_free_coordinates(
        48_000,
        modes
            .iter()
            .map(|&(w, zeta)| ModalAcousticMode {
                angular_frequency_rad_s: w,
                damping_ratio: zeta,
                pressure_per_modal_velocity: C64::ZERO,
            })
            .collect(),
        ModalAcousticTimeBudget {
            maximum_abs_displacement_m_sqrt_kg: 100.0,
            maximum_abs_velocity_m_sqrt_kg_per_s: 1.0e6,
            maximum_total_energy_j: 1.0e20,
            ..ModalAcousticTimeBudget::audible_reference()
        },
    )
    .unwrap()
}
fn prepare(m: &ModalAcousticTimeModel, frequencies: &[f64], h: f64) -> PreparedModalTones {
    PreparedModalTones::try_new(m, frequencies, h, ModalToneBudget::realtime_reference()).unwrap()
}
fn near(actual: f64, expected: f64, tolerance: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * scale.max(1.0e-300),
        "actual={actual:.17e} expected={expected:.17e} scale={scale:.17e}"
    );
}

// Elementary homogeneous solution and a complex particular solution form an
// oracle independent of the production entire-series/composed power factors.
fn free(w: f64, zeta: f64, q: f64, v: f64, t: f64) -> (f64, f64) {
    if w == 0.0 {
        return (q + v * t, v);
    }
    if zeta < 1.0 {
        let decay = zeta * w;
        let wd = w * (1.0 - zeta * zeta).sqrt();
        let b = (v + decay * q) / wd;
        let (s, c) = (wd * t).sin_cos();
        let envelope = (-decay * t).exp();
        (
            envelope * (q * c + b * s),
            envelope * ((wd * b - decay * q) * c - (wd * q + decay * b) * s),
        )
    } else if zeta == 1.0 {
        let b = v + w * q;
        let envelope = (-w * t).exp();
        (envelope * (q + b * t), envelope * (v - w * b * t))
    } else {
        let root = (zeta * zeta - 1.0).sqrt();
        let slow = -w / (zeta + root);
        let fast = -w * (zeta + root);
        let a = (v - fast * q) / (slow - fast);
        let b = q - a;
        (
            a * (slow * t).exp() + b * (fast * t).exp(),
            slow * a * (slow * t).exp() + fast * b * (fast * t).exp(),
        )
    }
}
fn analytic(
    w: f64,
    zeta: f64,
    frequencies: &[f64],
    fc: &[f64],
    fs: &[f64],
    initial: ModalAcousticState,
    t: f64,
) -> ModalAcousticState {
    let (mut q, mut v) = free(
        w,
        zeta,
        initial.displacement_m_sqrt_kg,
        initial.velocity_m_sqrt_kg_per_s,
        t,
    );
    for ((&drive, &c), &s) in frequencies.iter().zip(fc).zip(fs) {
        if drive == 0.0 && w == 0.0 {
            q += 0.5 * c * t * t;
            v += c * t;
        } else if w == 0.0 {
            let (sn, cs) = (drive * t).sin_cos();
            q += (c * (1.0 - cs) + s * (drive * t - sn)) / (drive * drive);
            v += (c * sn + s * (1.0 - cs)) / drive;
        } else if drive == w && zeta == 0.0 {
            let (sn, cs) = (w * t).sin_cos();
            q += c * t * sn / (2.0 * w) + s * (sn - w * t * cs) / (2.0 * w * w);
            v += c * (sn + w * t * cs) / (2.0 * w) + s * t * sn / 2.0;
        } else {
            let amplitude = C64::new(c, if drive == 0.0 { 0.0 } else { s })
                / C64::new(w * w - drive * drive, -2.0 * zeta * w * drive);
            let (sn, cs) = (drive * t).sin_cos();
            let transient = free(w, zeta, -amplitude.re, -drive * amplitude.im, t);
            q += amplitude.re * cs + amplitude.im * sn + transient.0;
            v += drive * (-amplitude.re * sn + amplitude.im * cs) + transient.1;
        }
    }
    state(q, v)
}
fn oracle_power(
    w: f64,
    zeta: f64,
    frequencies: &[f64],
    fc: &[f64],
    fs: &[f64],
    initial: ModalAcousticState,
    h: f64,
) -> (f64, f64) {
    let (mut work, mut heat) = (0.0, 0.0);
    let n = 16_384;
    for i in 0..=n {
        let t = h * i as f64 / n as f64;
        let v = analytic(w, zeta, frequencies, fc, fs, initial, t).velocity_m_sqrt_kg_per_s;
        let force = frequencies
            .iter()
            .zip(fc)
            .zip(fs)
            .map(|((&w, &c), &s)| c * (w * t).cos() + s * (w * t).sin())
            .sum::<f64>();
        let weight = if i == 0 || i == n {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        work += weight * force * v;
        heat += weight * 2.0 * zeta * w * v * v;
    }
    (work * h / (3.0 * n as f64), heat * h / (3.0 * n as f64))
}

#[test]
fn g1_complete_power_matches_continuous_inharmonic_dc_and_transient_oracle() {
    let frequencies = [0.0, 13.0, 13.4, 38.7];
    let fc = [0.3, 1.3, -0.7, 0.5];
    let fs = [99.0, -0.4, 0.9, 0.2];
    let initial = state(0.03, -0.4);
    for (w, zeta) in [
        (20.0, 0.07),
        (20.0, 1.0),
        (20.0, 3.0),
        (50.0, 1.0e-12),
        (0.0, 0.0),
    ] {
        let mut m = model(&[(w, zeta)]);
        m.restore_states(&[initial]).unwrap();
        let p = prepare(&m, &frequencies, 0.3);
        let preview = p.advance_mode(0, initial, &fc, &fs).unwrap();
        let mut workspace = ModalAcousticWorkspace::new(&m);
        let frame = m.step_tones_into(&p, &fc, &fs, &mut workspace).unwrap();
        let (work, heat) = oracle_power(w, zeta, &frequencies, &fc, &fs, initial, 0.3);
        near(frame.input_work_j, work, 5.0e-10, work.abs());
        near(frame.viscous_dissipation_j, heat, 5.0e-10, heat.abs());
        assert!(frame.viscous_dissipation_j >= 0.0);
        let expected = analytic(w, zeta, &frequencies, &fc, &fs, initial, 0.3);
        near(
            m.states()[0].displacement_m_sqrt_kg,
            expected.displacement_m_sqrt_kg,
            3.0e-11,
            0.1,
        );
        near(
            m.states()[0].velocity_m_sqrt_kg_per_s,
            expected.velocity_m_sqrt_kg_per_s,
            3.0e-11,
            1.0,
        );
        assert_eq!(preview, m.states()[0]);
        assert_eq!(frame.observer_pressure_pa, 0.0);
    }
}

#[test]
fn g1_close_opposing_tones_keep_small_beat_work_and_positive_heat() {
    let frequencies = [13.0, 13.000001];
    let (fc, fs) = ([1.0, -1.0], [0.3, -0.3]);
    let mut m = model(&[(20.0, 0.07)]);
    let p = prepare(&m, &frequencies, 0.7);
    let mut workspace = ModalAcousticWorkspace::new(&m);
    let frame = m.step_tones_into(&p, &fc, &fs, &mut workspace).unwrap();
    let (work, heat) = oracle_power(20.0, 0.07, &frequencies, &fc, &fs, state(0.0, 0.0), 0.7);
    assert!(work > 0.0 && heat > 0.0);
    near(frame.input_work_j, work, 2.0e-7, work);
    near(frame.viscous_dissipation_j, heat, 2.0e-7, heat);
    let mut separate_heat = 0.0;
    for index in 0..2 {
        let mut one = model(&[(20.0, 0.07)]);
        let single = PreparedModalHarmonic::try_new(&one, frequencies[index], 0.7).unwrap();
        let mut ws = ModalAcousticWorkspace::new(&one);
        separate_heat += one
            .step_harmonic_into(&single, &[fc[index]], &[fs[index]], &mut ws)
            .unwrap()
            .viscous_dissipation_j;
    }
    assert!(
        separate_heat > 1.0e10 * heat,
        "single-tone heat addition would erase destructive interference"
    );
}

#[test]
fn g1_duplicate_tones_cancel_before_motion_power_and_static_equilibrium() {
    let mut m = model(&[(7.0, 0.1)]);
    m.restore_states(&[state(0.03, -0.4)]).unwrap();
    let mut held = m.clone();
    let frequencies = [0.0, 13.0, 13.0, -0.0, 13.0];
    let p = prepare(&m, &frequencies, 0.3);
    let mut workspace = ModalAcousticWorkspace::new(&m);
    let frame = m
        .step_tones_into(
            &p,
            &[0.3, 1.0e100, -1.0e100, 1.0, 0.0],
            &[1.0e100, 1.0e100, -1.0e100, -1.0e100, 0.0],
            &mut workspace,
        )
        .unwrap();
    let expected = held.step_duration(&[1.3], 0.3).unwrap();
    near(frame.input_work_j, expected.input_work_j, 3.0e-11, 0.1);
    near(
        frame.viscous_dissipation_j,
        expected.viscous_dissipation_j,
        3.0e-11,
        expected.viscous_dissipation_j,
    );
    near(
        m.states()[0].velocity_m_sqrt_kg_per_s,
        held.states()[0].velocity_m_sqrt_kg_per_s,
        3.0e-11,
        1.0,
    );
    // Compensation preserves a small physical source between two large
    // exactly equal opposing tones; no near-frequency coalescence occurs.
    let mut compensated = model(&[(7.0, 0.1)]);
    let p = prepare(&compensated, &[13.0, 13.0, 13.0], 0.3);
    let mut reference = compensated.clone();
    let single = PreparedModalHarmonic::try_new(&reference, 13.0, 0.3).unwrap();
    let mut ws = ModalAcousticWorkspace::new(&compensated);
    let mut sws = ModalAcousticWorkspace::new(&reference);
    let frame = compensated
        .step_tones_into(&p, &[1.0e16, 1.0, -1.0e16], &[0.0; 3], &mut ws)
        .unwrap();
    let expected = reference
        .step_harmonic_into(&single, &[1.0], &[0.0], &mut sws)
        .unwrap();
    near(
        frame.input_work_j,
        expected.input_work_j,
        3.0e-11,
        expected.input_work_j,
    );
    near(
        frame.viscous_dissipation_j,
        expected.viscous_dissipation_j,
        3.0e-11,
        expected.viscous_dissipation_j,
    );
}

#[test]
fn g1_resonant_and_detuned_tones_share_one_homogeneous_state() {
    let frequencies = [0.0, 11.0, 17.2];
    let (fc, fs) = ([0.2, 0.7, -0.3], [0.0, 0.4, 0.6]);
    let initial = state(0.03, -0.2);
    let h = 5.0 * core::f64::consts::TAU / 11.0;
    let mut m = model(&[(11.0, 0.0)]);
    m.restore_states(&[initial]).unwrap();
    let p = prepare(&m, &frequencies, h);
    let mut workspace = ModalAcousticWorkspace::new(&m);
    let frame = m.step_tones_into(&p, &fc, &fs, &mut workspace).unwrap();
    let expected = analytic(11.0, 0.0, &frequencies, &fc, &fs, initial, h);
    near(
        m.states()[0].displacement_m_sqrt_kg,
        expected.displacement_m_sqrt_kg,
        3.0e-11,
        0.1,
    );
    near(
        m.states()[0].velocity_m_sqrt_kg_per_s,
        expected.velocity_m_sqrt_kg_per_s,
        3.0e-11,
        1.0,
    );
    let (work, _) = oracle_power(11.0, 0.0, &frequencies, &fc, &fs, initial, h);
    near(frame.input_work_j, work, 2.0e-10, work.abs());
    assert_eq!(frame.viscous_dissipation_j, 0.0);
}

fn run(
    m: &mut ModalAcousticTimeModel,
    p: &PreparedModalTones,
    c: &mut [f64],
    s: &mut [f64],
    steps: usize,
) -> (f64, f64) {
    let mut ws = ModalAcousticWorkspace::new(m);
    let (mut work, mut heat) = (0.0, 0.0);
    for _ in 0..steps {
        let frame = m.step_tones_into(p, c, s, &mut ws).unwrap();
        work += frame.input_work_j;
        heat += frame.viscous_dissipation_j;
        for (index, (c, s)) in c.iter_mut().zip(s.iter_mut()).enumerate() {
            (*c, *s) = p
                .advance_force_phase(index % p.tone_count(), *c, *s)
                .unwrap();
        }
    }
    (work, heat)
}

#[test]
fn g3_multimode_multitone_work_heat_subdivision_and_replay() {
    let mut source = model(&[(7.0, 0.05), (17.0, 1.0), (0.0, 0.0)]);
    source
        .restore_states(&[state(0.03, -0.1), state(-0.01, 0.2), state(0.1, 0.3)])
        .unwrap();
    let frequencies = [0.0, 9.0, 9.3];
    let forces = (
        [0.7, -0.2, 0.3, 0.4, 0.1, -0.1, -0.2, 0.5, -0.7],
        [0.1, 0.4, -0.5, 0.8, -0.3, 0.4, -0.7, 0.2, 0.4],
    );
    let mut whole = source.clone();
    let p = prepare(&whole, &frequencies, 0.31);
    let (mut c, mut s) = forces;
    let expected = run(&mut whole, &p, &mut c, &mut s, 1);
    for count in [4, 16] {
        let mut m = source.clone();
        let p = prepare(&m, &frequencies, 0.31 / count as f64);
        let (mut c, mut s) = forces;
        let actual = run(&mut m, &p, &mut c, &mut s, count);
        near(actual.0, expected.0, 5.0e-11, expected.0.abs());
        near(actual.1, expected.1, 5.0e-11, expected.1);
        for (actual, expected) in m.states().iter().zip(whole.states()) {
            near(
                actual.displacement_m_sqrt_kg,
                expected.displacement_m_sqrt_kg,
                5.0e-11,
                0.2,
            );
            near(
                actual.velocity_m_sqrt_kg_per_s,
                expected.velocity_m_sqrt_kg_per_s,
                5.0e-11,
                1.0,
            );
        }
        let checkpoint = (m.states().to_vec(), c, s);
        let forward = run(&mut m, &p, &mut c, &mut s, 7);
        let endpoint = m.states().to_vec();
        m.restore_states(&checkpoint.0).unwrap();
        let (mut c, mut s) = (checkpoint.1, checkpoint.2);
        assert_eq!(run(&mut m, &p, &mut c, &mut s, 7), forward);
        assert_eq!(m.states(), endpoint);
    }
}

#[test]
fn g1_near_dc_positive_heat_survives_another_lossless_mode() {
    let mut m = model(&[(1.0, 0.1), (400.0, 0.0)]);
    m.restore_states(&[state(1.0, 0.0), state(0.1, 100.0)])
        .unwrap();
    let p = prepare(&m, &[0.0, 1.0e-8, 3.1], 1.0e-3);
    let mut ws = ModalAcousticWorkspace::new(&m);
    let frame = m
        .step_tones_into(&p, &[0.4, 0.6, 0.0, 0.0, 0.0, 0.0], &[0.0; 6], &mut ws)
        .unwrap();
    // Only the 0.6*cos(Omega*t) part leaves initial static equilibrium.
    let expected = 0.2 * 0.6 * 0.6 * 1.0e-32 * 1.0e-21 / 252.0;
    near(
        frame.modal_viscous_losses[0].energy_j,
        expected,
        3.0e-4,
        expected,
    );
    assert_eq!(frame.modal_viscous_losses[1].energy_j, 0.0);
}

#[test]
fn g0_bounded_preparation_and_transactional_runtime_refusals() {
    let mut m = model(&[(7.0, 0.1), (13.0, 0.2)]);
    let b = ModalToneBudget::realtime_reference();
    for budget in [
        ModalToneBudget {
            maximum_tones: 1,
            ..b
        },
        ModalToneBudget {
            maximum_storage_bytes: 1,
            ..b
        },
        ModalToneBudget {
            maximum_preparation_work: 1,
            ..b
        },
        ModalToneBudget {
            maximum_step_work: 1,
            ..b
        },
    ] {
        assert!(PreparedModalTones::try_new(&m, &[9.0, 9.1], 0.2, budget).is_err());
    }
    for frequencies in [
        vec![],
        vec![1.0; MAX_HARMONIC_TONES + 1],
        vec![f64::NAN],
        vec![-1.0],
        vec![200_000.0],
    ] {
        assert!(PreparedModalTones::try_new(&m, &frequencies, 0.2, b).is_err());
    }
    let p = prepare(&m, &[9.0, 9.1], 0.2);
    let initial = m.states().to_vec();
    let mut ws = ModalAcousticWorkspace::new(&m);
    for c in [
        vec![1.0],
        vec![1.0, 0.0, 0.0, f64::NAN],
        vec![1.0, 0.0, 0.0, 1.0e6],
    ] {
        assert!(m.step_tones_into(&p, &c, &[0.0; 4], &mut ws).is_err());
        assert_eq!(m.states(), initial);
    }
    let mut unrelated = model(&[(7.0, 0.1), (13.0, 0.2)]);
    assert!(
        unrelated
            .step_tones_into(&p, &[0.0; 4], &[0.0; 4], &mut ws)
            .is_err()
    );
    assert!(
        p.advance_mode(2, state(0.0, 0.0), &[0.0; 2], &[0.0; 2])
            .is_err()
    );
    assert!(p.advance_force_phase(2, 0.0, 0.0).is_err());
    let mut narrow = ModalAcousticTimeModel::try_new(
        48_000,
        vec![ModalAcousticMode {
            angular_frequency_rad_s: 7.0,
            damping_ratio: 0.1,
            pressure_per_modal_velocity: C64::new(1.0, 0.2),
        }],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    assert!(PreparedModalTones::try_new(&narrow, &[9.0], 0.1, b).is_err());
    let mechanical = narrow.without_pressure_transfers();
    let p = prepare(&mechanical, &[9.0], 0.1);
    let mut ws = ModalAcousticWorkspace::new(&narrow);
    assert!(narrow.step_tones_into(&p, &[0.0], &[0.0], &mut ws).is_err());
}

#[test]
fn g0_highest_admitted_tone_is_retained_without_truncation() {
    let frequencies = (0..MAX_HARMONIC_TONES)
        .map(|i| 100.0 + 17.3 * i as f64)
        .collect::<Vec<_>>();
    let mut m = model(&[(700.0, 0.02)]);
    let p = prepare(&m, &frequencies, 1.0 / 48_000.0);
    assert_eq!(p.tone_count(), MAX_HARMONIC_TONES);
    assert_eq!(p.mode_count(), 1);
    assert_eq!(p.angular_frequencies_rad_s(), frequencies);
    let mut c = [0.0; MAX_HARMONIC_TONES];
    c[MAX_HARMONIC_TONES - 1] = 2.0;
    let s = [0.0; MAX_HARMONIC_TONES];
    let mut reference = m.clone();
    let single = PreparedModalHarmonic::try_new(
        &reference,
        frequencies[MAX_HARMONIC_TONES - 1],
        p.duration_s(),
    )
    .unwrap();
    let mut ws = ModalAcousticWorkspace::new(&m);
    let mut sws = ModalAcousticWorkspace::new(&reference);
    let actual = m.step_tones_into(&p, &c, &s, &mut ws).unwrap();
    let expected = reference
        .step_harmonic_into(&single, &[2.0], &[0.0], &mut sws)
        .unwrap();
    near(
        actual.input_work_j,
        expected.input_work_j,
        3.0e-11,
        expected.input_work_j,
    );
    near(
        actual.viscous_dissipation_j,
        expected.viscous_dissipation_j,
        3.0e-11,
        expected.viscous_dissipation_j,
    );
    near(
        m.states()[0].velocity_m_sqrt_kg_per_s,
        reference.states()[0].velocity_m_sqrt_kg_per_s,
        3.0e-11,
        1.0e-4,
    );
}
