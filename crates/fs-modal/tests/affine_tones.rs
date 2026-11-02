//! Physical attack/release trajectories, independent work and positive heat.
use fs_math::c64::C64;
use fs_modal::time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
    ModalAcousticWorkspace, ModalToneBudget, PreparedModalTones,
};

fn model(w: f64, zeta: f64) -> ModalAcousticTimeModel {
    ModalAcousticTimeModel::try_new_with_free_coordinates(
        48_000,
        vec![ModalAcousticMode {
            angular_frequency_rad_s: w,
            damping_ratio: zeta,
            pressure_per_modal_velocity: C64::ZERO,
        }],
        ModalAcousticTimeBudget {
            maximum_abs_displacement_m_sqrt_kg: 100.0,
            maximum_abs_velocity_m_sqrt_kg_per_s: 1e5,
            maximum_total_energy_j: 1e12,
            ..ModalAcousticTimeBudget::audible_reference()
        },
    )
    .unwrap()
}
fn prepare(model: &ModalAcousticTimeModel, w: &[f64], h: f64) -> PreparedModalTones {
    PreparedModalTones::try_new_affine(model, w, h, ModalToneBudget::realtime_reference()).unwrap()
}
fn near(a: f64, b: f64, scale: f64) {
    assert!(
        (a - b).abs() < 3e-11 * scale.max(1e-300),
        "actual={a:.17e}, expected={b:.17e}"
    );
}

#[test]
fn affine_dc_and_exact_resonance_match_independent_closed_forms() {
    let h = 0.37;
    let mut free = model(0.0, 0.0);
    let p = prepare(&free, &[0.0], h);
    let mut workspace = ModalAcousticWorkspace::new(&free);
    let frame = free
        .step_tones_affine_into(&p, &[2.0], &[8.0], &[6.0], &[-9.0], &mut workspace)
        .unwrap();
    let q = h * h * (1.0 + 4.0 / 6.0);
    let v = 4.0 * h;
    near(free.states()[0].displacement_m_sqrt_kg, q, q);
    near(free.states()[0].velocity_m_sqrt_kg_per_s, v, v);
    near(frame.input_work_j, 0.5 * v * v, v * v);
    assert_eq!(frame.viscous_dissipation_j, 0.0);

    let w = 11.0;
    let mut resonant = model(w, 0.0);
    let p = prepare(&resonant, &[w], h);
    let mut workspace = ModalAcousticWorkspace::new(&resonant);
    let frame = resonant
        .step_tones_affine_into(&p, &[0.0], &[0.0], &[1.0], &[0.0], &mut workspace)
        .unwrap();
    let (s, c) = (w * h).sin_cos();
    let q = h * s / (4.0 * w) + c / (4.0 * w * w) - s / (4.0 * w * w * w * h);
    let v = h * c / 4.0 + s / (4.0 * w);
    near(resonant.states()[0].displacement_m_sqrt_kg, q, h * h);
    near(resonant.states()[0].velocity_m_sqrt_kg_per_s, v, h);
    near(frame.input_work_j, 0.5 * (v * v + w * w * q * q), h * h);
    assert_eq!(frame.viscous_dissipation_j, 0.0);
}

fn force(t: f64, h: f64, omega: &[f64], c0: &[f64], s0: &[f64], c1: &[f64], s1: &[f64]) -> f64 {
    (0..omega.len())
        .map(|i| {
            let (s, c) = (omega[i] * t).sin_cos();
            (c0[i] + (c1[i] - c0[i]) * t / h) * c + (s0[i] + (s1[i] - s0[i]) * t / h) * s
        })
        .sum()
}
// Independent physical ODE integration, including work and heat as states.
fn reference(
    w: f64,
    zeta: f64,
    h: f64,
    omega: &[f64],
    c0: &[f64],
    s0: &[f64],
    c1: &[f64],
    s1: &[f64],
) -> [f64; 4] {
    let mut x = [0.01, -0.03, 0.0, 0.0];
    let dt = h / 8192.0;
    let derivative = |t: f64, x: [f64; 4]| {
        let f = force(t, h, omega, c0, s0, c1, s1);
        [
            x[1],
            f - w * w * x[0] - 2.0 * zeta * w * x[1],
            f * x[1],
            2.0 * zeta * w * x[1] * x[1],
        ]
    };
    for i in 0..8192 {
        let t = i as f64 * dt;
        let a = derivative(t, x);
        let b = derivative(
            t + dt / 2.0,
            core::array::from_fn(|j| x[j] + dt * a[j] / 2.0),
        );
        let c = derivative(
            t + dt / 2.0,
            core::array::from_fn(|j| x[j] + dt * b[j] / 2.0),
        );
        let d = derivative(t + dt, core::array::from_fn(|j| x[j] + dt * c[j]));
        for j in 0..4 {
            x[j] += dt * (a[j] + 2.0 * b[j] + 2.0 * c[j] + d[j]) / 6.0;
        }
    }
    x
}

#[test]
fn nearby_envelopes_retain_cross_work_and_heat_across_damping_regimes() {
    let omega = [13.0, 13.0007, 0.0];
    let (c0, s0, c1, s1) = (
        [1.0, -0.8, 0.2],
        [0.3, -0.1, 0.0],
        [0.2, 0.7, -0.1],
        [-0.2, 0.4, 0.0],
    );
    for (w, zeta) in [(11.0, 0.08), (11.0, 1.0), (11.0, 1.7), (0.0, 0.0)] {
        let mut m = model(w, zeta);
        m.restore_states(&[ModalAcousticState {
            displacement_m_sqrt_kg: 0.01,
            velocity_m_sqrt_kg_per_s: -0.03,
        }])
        .unwrap();
        let p = prepare(&m, &omega, 0.31);
        let mut scratch = ModalAcousticWorkspace::new(&m);
        let frame = m
            .step_tones_affine_into(&p, &c0, &s0, &c1, &s1, &mut scratch)
            .unwrap();
        let expected = reference(w, zeta, 0.31, &omega, &c0, &s0, &c1, &s1);
        near(m.states()[0].displacement_m_sqrt_kg, expected[0], 0.1);
        near(m.states()[0].velocity_m_sqrt_kg_per_s, expected[1], 1.0);
        near(frame.input_work_j, expected[2], 0.1);
        near(frame.viscous_dissipation_j, expected[3], 0.1);
        assert!(frame.viscous_dissipation_j >= 0.0);
    }
}

#[test]
fn envelope_partition_preserves_motion_work_heat_and_same_basis_phase() {
    let omega = [9.0, 9.00003];
    let h = 0.43;
    let (c0, s0, c1, s1) = ([1.0, -0.6], [0.7, -0.3], [-0.2, 0.0], [0.0, 0.2]);
    let mut whole = model(8.0, 0.03);
    let mut split = whole.clone();
    let p = prepare(&whole, &omega, h);
    let mut scratch = ModalAcousticWorkspace::new(&whole);
    let expected = whole
        .step_tones_affine_into(&p, &c0, &s0, &c1, &s1, &mut scratch)
        .unwrap()
        .clone();
    let p = prepare(&split, &omega, h / 17.0);
    let mut work = 0.0;
    let mut heat = 0.0;
    for segment in 0..17 {
        let f0 = segment as f64 / 17.0;
        let f1 = (segment + 1) as f64 / 17.0;
        let mut rows = [[0.0; 2]; 4];
        for i in 0..2 {
            let (s, c) = (omega[i] * h * f0).sin_cos();
            for (end, f) in [(0, f0), (2, f1)] {
                let cc = c0[i] + (c1[i] - c0[i]) * f;
                let ss = s0[i] + (s1[i] - s0[i]) * f;
                rows[end][i] = cc * c + ss * s;
                rows[end + 1][i] = ss * c - cc * s;
            }
        }
        let frame = split
            .step_tones_affine_into(&p, &rows[0], &rows[1], &rows[2], &rows[3], &mut scratch)
            .unwrap();
        work += frame.input_work_j;
        heat += frame.viscous_dissipation_j;
    }
    near(
        split.states()[0].displacement_m_sqrt_kg,
        whole.states()[0].displacement_m_sqrt_kg,
        0.1,
    );
    near(
        split.states()[0].velocity_m_sqrt_kg_per_s,
        whole.states()[0].velocity_m_sqrt_kg_per_s,
        1.0,
    );
    near(work, expected.input_work_j, 0.1);
    near(heat, expected.viscous_dissipation_j, 0.1);
}

#[test]
fn preparation_family_duplicate_cancellation_and_rejection_preserve_state() {
    let mut m = model(8.0, 0.1);
    let p = prepare(&m, &[9.0, 9.0], 0.1);
    let mut scratch = ModalAcousticWorkspace::new(&m);
    let frame = m
        .step_tones_affine_into(
            &p,
            &[1e308, -1e308],
            &[0.0; 2],
            &[-1e308, 1e308],
            &[0.0; 2],
            &mut scratch,
        )
        .unwrap();
    assert_eq!(frame.input_work_j, 0.0);
    assert_eq!(frame.viscous_dissipation_j, 0.0);
    assert_eq!(m.states()[0], ModalAcousticState::default());
    assert!(
        m.step_tones_into(&p, &[0.0; 2], &[0.0; 2], &mut scratch)
            .is_err()
    );
    assert!(
        m.step_tones_affine_into(
            &p,
            &[1e200; 2],
            &[0.0; 2],
            &[1e200; 2],
            &[0.0; 2],
            &mut scratch
        )
        .is_err()
    );
    assert_eq!(m.states()[0], ModalAcousticState::default());
    let old = PreparedModalTones::try_new(&m, &[9.0], 0.1, ModalToneBudget::realtime_reference())
        .unwrap();
    assert!(!old.is_affine());
    assert!(
        m.step_tones_affine_into(&old, &[0.0], &[0.0], &[0.0], &[0.0], &mut scratch)
            .is_err()
    );
    assert!(
        PreparedModalTones::try_new_affine(
            &m,
            &[9.0],
            0.1,
            ModalToneBudget {
                maximum_storage_bytes: 1,
                ..ModalToneBudget::realtime_reference()
            }
        )
        .is_err()
    );
}
