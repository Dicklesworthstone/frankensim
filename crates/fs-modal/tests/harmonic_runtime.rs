//! G1/G3 physical harmonic loading through the actual shared modal runtime.

use fs_math::c64::C64;
use fs_modal::time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
    ModalAcousticWorkspace, PreparedModalHarmonic,
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

fn near(actual: f64, expected: f64, relative: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= relative * scale.max(1.0e-300),
        "actual={actual:.17e}, expected={expected:.17e}, scale={scale:.17e}"
    );
}

// Independent particular solution plus elementary homogeneous solutions.
// No production transition or power integral enters this oracle.
fn analytic(
    w: f64,
    zeta: f64,
    drive: f64,
    fc: f64,
    fs: f64,
    initial: ModalAcousticState,
    t: f64,
) -> ModalAcousticState {
    let (s, c) = (drive * t).sin_cos();
    if w == 0.0 {
        return state(
            initial.displacement_m_sqrt_kg
                + initial.velocity_m_sqrt_kg_per_s * t
                + (fc * (1.0 - c) + fs * (drive * t - s)) / (drive * drive),
            initial.velocity_m_sqrt_kg_per_s + (fc * s + fs * (1.0 - c)) / drive,
        );
    }
    let q = C64::new(fc, fs) / C64::new(w * w - drive * drive, -2.0 * zeta * w * drive);
    let a = initial.displacement_m_sqrt_kg - q.re;
    let v = initial.velocity_m_sqrt_kg_per_s - drive * q.im;
    let (free_q, free_v) = if zeta < 1.0 {
        let decay = zeta * w;
        let wd = w * (1.0 - zeta * zeta).sqrt();
        let b = (v + decay * a) / wd;
        let (s, c) = (wd * t).sin_cos();
        let envelope = (-decay * t).exp();
        (
            envelope * (a * c + b * s),
            envelope * ((wd * b - decay * a) * c - (wd * a + decay * b) * s),
        )
    } else if zeta == 1.0 {
        let b = v + w * a;
        let envelope = (-w * t).exp();
        (envelope * (a + b * t), envelope * (v - w * b * t))
    } else {
        let root = (zeta * zeta - 1.0).sqrt();
        let slow = -w / (zeta + root);
        let fast = -w * (zeta + root);
        let p = (v - fast * a) / (slow - fast);
        let r = a - p;
        (
            p * (slow * t).exp() + r * (fast * t).exp(),
            slow * p * (slow * t).exp() + fast * r * (fast * t).exp(),
        )
    };
    state(
        q.re * c + q.im * s + free_q,
        drive * (-q.re * s + q.im * c) + free_v,
    )
}

#[test]
fn g1_work_and_heat_match_independent_continuous_power_quadrature() {
    for (w, zeta, drive, h) in [
        (20.0, 0.07, 13.0, 0.3),
        (20.0, 1.0, 13.0, 0.3),
        (20.0, 3.0, 13.0, 0.3),
        (100.0, 1.0e-12, 97.0, 0.1),
        (0.0, 0.0, 13.0, 0.3),
    ] {
        let (fc, fs) = (1.3, -0.7);
        let initial = state(0.03, -0.4);
        let mut m = model(&[(w, zeta)]);
        m.restore_states(&[initial]).unwrap();
        let p = PreparedModalHarmonic::try_new(&m, drive, h).unwrap();
        let mut workspace = ModalAcousticWorkspace::new(&m);
        let frame = m
            .step_harmonic_into(&p, &[fc], &[fs], &mut workspace)
            .unwrap();
        let (mut work, mut heat) = (0.0, 0.0);
        let n = 16_384;
        for i in 0..=n {
            let t = h * i as f64 / n as f64;
            let v = analytic(w, zeta, drive, fc, fs, initial, t).velocity_m_sqrt_kg_per_s;
            let force = fc * (drive * t).cos() + fs * (drive * t).sin();
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
        work *= h / (3.0 * n as f64);
        heat *= h / (3.0 * n as f64);
        near(frame.input_work_j, work, 2.0e-10, work.abs());
        near(frame.viscous_dissipation_j, heat, 2.0e-10, heat.abs());
        assert!(frame.viscous_dissipation_j >= 0.0);
        assert_eq!(frame.observer_pressure_pa, 0.0);
        let expected = analytic(w, zeta, drive, fc, fs, initial, h);
        near(
            m.states()[0].displacement_m_sqrt_kg,
            expected.displacement_m_sqrt_kg,
            2.0e-11,
            0.1,
        );
        near(
            m.states()[0].velocity_m_sqrt_kg_per_s,
            expected.velocity_m_sqrt_kg_per_s,
            2.0e-11,
            1.0,
        );
    }
}

#[test]
fn g1_resonant_cycle_work_is_not_endpoint_force_times_displacement() {
    let (w, fc) = (11.0, 0.7);
    let h = 3.0 * core::f64::consts::TAU / w;
    let mut m = model(&[(w, 0.0)]);
    let p = PreparedModalHarmonic::try_new(&m, w, h).unwrap();
    let mut workspace = ModalAcousticWorkspace::new(&m);
    let frame = m
        .step_harmonic_into(&p, &[fc], &[0.0], &mut workspace)
        .unwrap();
    let expected = fc * fc * h * h / 8.0;
    near(frame.input_work_j, expected, 1.0e-11, expected);
    near(frame.total_modal_energy_j, expected, 1.0e-11, expected);
    assert_eq!(frame.viscous_dissipation_j, 0.0);
    assert_eq!(frame.modal_viscous_losses[0].energy_j, 0.0);
    assert!(frame.input_work_j > 1.0e6 * (fc * m.states()[0].displacement_m_sqrt_kg).abs());
}

#[test]
fn g1_quadrature_interference_is_retained_in_positive_heat() {
    let run = |fc, fs| {
        let mut m = model(&[(13.0, 0.1)]);
        let p = PreparedModalHarmonic::try_new(&m, 7.0, 0.2).unwrap();
        let mut workspace = ModalAcousticWorkspace::new(&m);
        let frame = m
            .step_harmonic_into(&p, &[fc], &[fs], &mut workspace)
            .unwrap();
        (frame.viscous_dissipation_j, frame.input_work_j)
    };
    let combined = run(0.8, 0.9);
    let cosine = run(0.8, 0.0);
    let sine = run(0.0, 0.9);
    assert!((combined.0 - cosine.0 - sine.0).abs() > 0.01 * combined.0);
    assert!((combined.1 - cosine.1 - sine.1).abs() > 0.01 * combined.1.abs());
}

#[test]
fn g1_dc_work_heat_and_motion_match_the_existing_physical_zoh_owner() {
    for (w, zeta) in [
        (0.0, 0.0),
        (7.0, 0.0),
        (7.0, 0.05),
        (7.0, 1.0),
        (7.0, 4.0),
        (100.0, 1.0e12),
    ] {
        let mut harmonic = model(&[(w, zeta)]);
        harmonic.restore_states(&[state(0.03, -0.4)]).unwrap();
        let mut held = harmonic.clone();
        let p = PreparedModalHarmonic::try_new(&harmonic, 0.0, 0.3).unwrap();
        let mut workspace = ModalAcousticWorkspace::new(&harmonic);
        let actual = harmonic
            .step_harmonic_into(&p, &[1.3], &[-7.0], &mut workspace)
            .unwrap();
        let expected = held.step_duration(&[1.3], 0.3).unwrap();
        near(actual.input_work_j, expected.input_work_j, 2.0e-11, 0.1);
        near(
            actual.viscous_dissipation_j,
            expected.viscous_dissipation_j,
            2.0e-11,
            expected.viscous_dissipation_j.abs(),
        );
        near(
            actual.total_modal_energy_j,
            expected.total_modal_energy_j,
            2.0e-11,
            expected.total_modal_energy_j.abs(),
        );
        near(
            harmonic.states()[0].displacement_m_sqrt_kg,
            held.states()[0].displacement_m_sqrt_kg,
            2.0e-11,
            0.1,
        );
        near(
            harmonic.states()[0].velocity_m_sqrt_kg_per_s,
            held.states()[0].velocity_m_sqrt_kg_per_s,
            2.0e-11,
            1.0,
        );
    }
}

fn run_steps(
    m: &mut ModalAcousticTimeModel,
    p: &PreparedModalHarmonic,
    fc: &mut [f64],
    fs: &mut [f64],
    steps: usize,
) -> (f64, f64) {
    let mut workspace = ModalAcousticWorkspace::new(m);
    let (mut work, mut heat) = (0.0, 0.0);
    for _ in 0..steps {
        let frame = m.step_harmonic_into(p, fc, fs, &mut workspace).unwrap();
        work += frame.input_work_j;
        heat += frame.viscous_dissipation_j;
        for (c, s) in fc.iter_mut().zip(fs.iter_mut()) {
            (*c, *s) = p.advance_force_phase(*c, *s).unwrap();
        }
    }
    (work, heat)
}

#[test]
fn g3_multimode_work_and_heat_compose_and_checkpoint_replay() {
    let mut source = model(&[(7.0, 0.05), (17.0, 1.0), (0.0, 0.0)]);
    source
        .restore_states(&[state(0.03, -0.1), state(-0.01, 0.2), state(0.1, 0.3)])
        .unwrap();
    let forces = ([0.7, -0.2, 0.3], [0.1, 0.4, -0.5]);
    let mut whole = source.clone();
    let full = PreparedModalHarmonic::try_new(&whole, 9.0, 0.31).unwrap();
    let (mut full_fc, mut full_fs) = forces;
    let expected = run_steps(&mut whole, &full, &mut full_fc, &mut full_fs, 1);
    for count in [4, 16] {
        let mut m = source.clone();
        let p = PreparedModalHarmonic::try_new(&m, 9.0, 0.31 / count as f64).unwrap();
        let (mut fc, mut fs) = forces;
        let actual = run_steps(&mut m, &p, &mut fc, &mut fs, count);
        near(actual.0, expected.0, 3.0e-11, expected.0.abs());
        near(actual.1, expected.1, 3.0e-11, expected.1.abs());
        for (actual, expected) in m.states().iter().zip(whole.states()) {
            near(
                actual.displacement_m_sqrt_kg,
                expected.displacement_m_sqrt_kg,
                3.0e-11,
                0.2,
            );
            near(
                actual.velocity_m_sqrt_kg_per_s,
                expected.velocity_m_sqrt_kg_per_s,
                3.0e-11,
                1.0,
            );
        }
        let checkpoint = (m.states().to_vec(), fc, fs);
        let forward = run_steps(&mut m, &p, &mut fc, &mut fs, 7);
        let endpoint = m.states().to_vec();
        m.restore_states(&checkpoint.0).unwrap();
        let (mut fc, mut fs) = (checkpoint.1, checkpoint.2);
        let replay = run_steps(&mut m, &p, &mut fc, &mut fs, 7);
        assert_eq!(endpoint, m.states());
        assert_eq!(forward, replay);
    }
}

#[test]
fn g1_resolved_near_dc_heat_survives_large_other_mode_energy() {
    let mut m = model(&[(1.0, 0.1), (400.0, 0.0)]);
    m.restore_states(&[state(1.0, 0.0), state(0.1, 100.0)])
        .unwrap();
    let p = PreparedModalHarmonic::try_new(&m, 1.0e-8, 1.0e-3).unwrap();
    let mut workspace = ModalAcousticWorkspace::new(&m);
    let frame = m
        .step_harmonic_into(&p, &[1.0, 0.0], &[0.0, 0.0], &mut workspace)
        .unwrap();
    // At static equilibrium v=-F*Omega^2*t^3/6 + O(t^4), giving
    // heat=c*F^2*Omega^4*h^7/252. This is far below total-energy ulps.
    let expected = 0.2 * 1.0e-32 * 1.0e-21 / 252.0;
    assert!(frame.modal_viscous_losses[0].energy_j > 0.0);
    near(
        frame.modal_viscous_losses[0].energy_j,
        expected,
        3.0e-4,
        expected,
    );
    assert_eq!(frame.modal_viscous_losses[1].energy_j, 0.0);
    assert_eq!(
        frame.viscous_dissipation_j,
        frame.modal_viscous_losses[0].energy_j
    );
}

#[test]
fn g0_runtime_refusals_preserve_every_accepted_state_and_observation_scope() {
    let mut m = model(&[(7.0, 0.1), (13.0, 0.2)]);
    let p = PreparedModalHarmonic::try_new(&m, 9.0, 0.2).unwrap();
    let initial = m.states().to_vec();
    let mut workspace = ModalAcousticWorkspace::new(&m);
    for (fc, fs) in [
        (vec![1.0], vec![0.0, 0.0]),
        (vec![1.0, f64::NAN], vec![0.0, 0.0]),
        (vec![1.0, 1.0e6], vec![0.0, 0.0]),
        (vec![1.0, 1.0e308], vec![0.0, 0.0]),
    ] {
        assert!(m.step_harmonic_into(&p, &fc, &fs, &mut workspace).is_err());
        assert_eq!(m.states(), initial);
    }
    let mut unrelated = model(&[(7.0, 0.1), (13.0, 0.2)]);
    assert!(
        unrelated
            .step_harmonic_into(&p, &[1.0, 1.0], &[0.0, 0.0], &mut workspace)
            .is_err()
    );
    assert!(PreparedModalHarmonic::try_new(&m, 200_000.0, 0.1).is_err());
    assert!(PreparedModalHarmonic::try_new(&m, 9.0, 0.0).is_err());
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
    assert!(PreparedModalHarmonic::try_new(&narrow, 9.0, 0.1).is_err());
    let mechanical = narrow.without_pressure_transfers();
    let p = PreparedModalHarmonic::try_new(&mechanical, 9.0, 0.1).unwrap();
    let mut workspace = ModalAcousticWorkspace::new(&narrow);
    assert!(
        narrow
            .step_harmonic_into(&p, &[1.0], &[0.0], &mut workspace)
            .is_err()
    );
    assert_eq!(narrow.states(), &[ModalAcousticState::default()]);
}
