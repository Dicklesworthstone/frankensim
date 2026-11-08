//! G1/G3/G4 physical-rate banks: interference, ramps, rigid drag and ownership.
use fs_math::c64::C64;
use fs_modal::time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeError,
    ModalAcousticTimeModel, ModalAcousticWorkspace, ModalToneBudget, PreparedModalTones,
};

fn state(q: f64, v: f64) -> ModalAcousticState {
    ModalAcousticState {
        displacement_m_sqrt_kg: q,
        velocity_m_sqrt_kg_per_s: v,
    }
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 3e-11, "actual={a:e}, expected={b:e}");
}
fn prepare(omega: &[f64], rates: &[f64], frequencies: &[f64], h: f64) -> PreparedModalTones {
    PreparedModalTones::try_new_with_damping_rates(
        omega,
        rates,
        frequencies,
        h,
        ModalToneBudget::realtime_reference(),
    )
    .unwrap()
}

// Independent physical ODE integration of displacement, velocity, work and heat.
// The damping coefficient is c, not the amplitude decay c/2.
fn reference(
    omega: f64,
    rate: f64,
    h: f64,
    frequencies: &[f64],
    c0: &[f64],
    s0: &[f64],
    c1: &[f64],
    s1: &[f64],
    initial: ModalAcousticState,
) -> [f64; 4] {
    let mut x = [
        initial.displacement_m_sqrt_kg,
        initial.velocity_m_sqrt_kg_per_s,
        0.0,
        0.0,
    ];
    let dt = h / 4096.0;
    let derivative = |t: f64, x: [f64; 4]| {
        let force: f64 = (0..frequencies.len())
            .map(|i| {
                let (sin, cos) = (frequencies[i] * t).sin_cos();
                (c0[i] + (c1[i] - c0[i]) * t / h) * cos + (s0[i] + (s1[i] - s0[i]) * t / h) * sin
            })
            .sum();
        [
            x[1],
            force - omega * omega * x[0] - rate * x[1],
            force * x[1],
            rate * x[1] * x[1],
        ]
    };
    for i in 0..4096 {
        let t = i as f64 * dt;
        let a = derivative(t, x);
        let b = derivative(
            t + 0.5 * dt,
            core::array::from_fn(|j| x[j] + 0.5 * dt * a[j]),
        );
        let c = derivative(
            t + 0.5 * dt,
            core::array::from_fn(|j| x[j] + 0.5 * dt * b[j]),
        );
        let d = derivative(t + dt, core::array::from_fn(|j| x[j] + dt * c[j]));
        for j in 0..4 {
            x[j] += dt * (a[j] + 2.0 * b[j] + 2.0 * c[j] + d[j]) / 6.0;
        }
    }
    x
}

#[test]
fn exact_duplicates_cancel_before_power_but_near_tones_retain_beats() {
    let zero = state(0.0, 0.0);
    let exact = prepare(&[2.0], &[0.5], &[3.0, 3.0], 0.7);
    let cancelled = exact
        .advance_mode_with_energy(0, zero, &[1.0, -1.0], &[0.2, -0.2])
        .unwrap();
    assert_eq!(cancelled.state, zero);
    assert_eq!(cancelled.input_work_j, 0.0);
    assert_eq!(cancelled.viscous_dissipation_j, 0.0);
    let frequencies = [3.0, 3.001];
    let c = [1.0, -1.0];
    let s = [0.2, -0.2];
    let beats = prepare(&[2.0], &[0.5], &frequencies, 0.7);
    assert_eq!(beats.angular_frequencies_rad_s(), frequencies);
    let frame = beats.advance_mode_with_energy(0, zero, &c, &s).unwrap();
    let expected = reference(2.0, 0.5, 0.7, &frequencies, &c, &s, &c, &s, zero);
    near(frame.state.displacement_m_sqrt_kg, expected[0]);
    near(frame.state.velocity_m_sqrt_kg_per_s, expected[1]);
    assert!(frame.viscous_dissipation_j > 1e-12);
    assert!((frame.viscous_dissipation_j - expected[3]).abs() < expected[3] * 1e-7);
    assert!((frame.input_work_j - expected[2]).abs() < expected[2].abs() * 1e-7);
}

#[test]
fn affine_envelopes_keep_all_cross_terms_across_physical_damping_regimes() {
    let omega = [4.0, 4.0, 4.0, 0.0];
    let rates = [0.4, 8.0, 12.0, 2.5];
    let frequencies = [9.0, 9.002, 0.0];
    let (c0, s0, c1, s1) = (
        [1.0, -0.8, 0.2],
        [0.3, -0.1, 0.0],
        [0.2, 0.7, -0.1],
        [-0.2, 0.4, 0.0],
    );
    let initial = state(0.01, -0.03);
    let bank = PreparedModalTones::try_new_affine_with_damping_rates(
        &omega,
        &rates,
        &frequencies,
        0.31,
        ModalToneBudget::realtime_reference(),
    )
    .unwrap();
    assert_eq!(bank.mode_count(), 4);
    assert!(bank.is_affine());
    for i in 0..omega.len() {
        let frame = bank
            .advance_mode_affine_with_energy(i, initial, &c0, &s0, &c1, &s1)
            .unwrap();
        let expected = reference(
            omega[i],
            rates[i],
            0.31,
            &frequencies,
            &c0,
            &s0,
            &c1,
            &s1,
            initial,
        );
        near(frame.state.displacement_m_sqrt_kg, expected[0]);
        near(frame.state.velocity_m_sqrt_kg_per_s, expected[1]);
        near(frame.input_work_j, expected[2]);
        near(frame.viscous_dissipation_j, expected[3]);
        assert!(frame.viscous_dissipation_j > 0.0);
    }
    assert!(bank.advance_mode_with_energy(0, initial, &c0, &s0).is_err());
}

#[test]
fn free_drag_uses_physical_rate_without_a_sampling_clock() {
    let (h, rate, force) = (0.7_f64, 2.5_f64, 1.2_f64);
    let initial = state(-0.2, 0.3);
    let p = prepare(&[0.0], &[rate], &[0.0], h);
    let frame = p
        .advance_mode_with_energy(0, initial, &[force], &[0.0])
        .unwrap();
    let terminal = force / rate;
    let tail = initial.velocity_m_sqrt_kg_per_s - terminal;
    let i0 = -(-rate * h).exp_m1() / rate;
    let i2 = -(-2.0 * rate * h).exp_m1() / (2.0 * rate);
    let distance = terminal * h + tail * i0;
    near(
        frame.state.displacement_m_sqrt_kg,
        initial.displacement_m_sqrt_kg + distance,
    );
    near(
        frame.state.velocity_m_sqrt_kg_per_s,
        terminal + tail * (-rate * h).exp(),
    );
    near(frame.input_work_j, force * distance);
    near(
        frame.viscous_dissipation_j,
        rate * (terminal * terminal * h + 2.0 * terminal * tail * i0 + tail * tail * i2),
    );
    let high = prepare(&[2e6], &[1e3], &[1e6], 2e-7);
    assert!(
        high.advance_mode_with_energy(0, state(0.0, 0.0), &[1.0], &[0.0])
            .unwrap()
            .state
            .velocity_m_sqrt_kg_per_s
            .is_finite()
    );
}

#[test]
fn model_binding_remains_required_and_legacy_frames_share_the_same_owner() {
    let mut model = ModalAcousticTimeModel::try_new_with_free_coordinates(
        48_000,
        vec![ModalAcousticMode {
            angular_frequency_rad_s: 4.0,
            damping_ratio: 0.1,
            pressure_per_modal_velocity: C64::ZERO,
        }],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    let initial = state(0.01, -0.03);
    model.restore_states(&[initial]).unwrap();
    let budget = ModalToneBudget::realtime_reference();
    let bound = PreparedModalTones::try_new(&model, &[3.0], 0.1, budget).unwrap();
    let pure = bound
        .advance_mode_with_energy(0, initial, &[0.8], &[0.2])
        .unwrap();
    let mut scratch = ModalAcousticWorkspace::new(&model);
    {
        let accepted = model
            .step_tones_into(&bound, &[0.8], &[0.2], &mut scratch)
            .unwrap();
        assert_eq!(accepted.input_work_j, pure.input_work_j);
        assert_eq!(accepted.viscous_dissipation_j, pure.viscous_dissipation_j);
    }
    assert_eq!(model.states()[0], pure.state);
    let unbound = prepare(&[4.0], &[0.8], &[3.0], 0.1);
    let before = model.states().to_vec();
    assert!(
        model
            .step_tones_into(&unbound, &[0.8], &[0.2], &mut scratch)
            .is_err()
    );
    assert_eq!(model.states(), before);
    assert!(PreparedModalTones::try_new(&model, &[1e6], 1e-6, budget).is_err());
    assert!(
        PreparedModalTones::try_new_with_damping_rates(&[4.0], &[0.8], &[1e6], 1e-6, budget)
            .is_ok()
    );
}

#[test]
fn aggregate_admission_and_preparation_cancellation_refuse_cleanly() {
    let budget = ModalToneBudget::realtime_reference();
    assert!(PreparedModalTones::try_new_with_damping_rates(&[], &[], &[1.0], 0.1, budget).is_err());
    assert!(
        PreparedModalTones::try_new_with_damping_rates(&[1.0], &[-1e-30], &[1.0], 0.1, budget)
            .is_err()
    );
    assert!(
        PreparedModalTones::try_new_with_damping_rates(&[1.0, 2.0], &[0.5], &[1.0], 0.1, budget)
            .is_err()
    );
    for limited in [
        ModalToneBudget {
            maximum_storage_bytes: 0,
            ..budget
        },
        ModalToneBudget {
            maximum_preparation_work: 0,
            ..budget
        },
        ModalToneBudget {
            maximum_step_work: 0,
            ..budget
        },
    ] {
        assert!(matches!(
            PreparedModalTones::try_new_with_damping_rates(&[1.0], &[0.5], &[1.0], 0.1, limited),
            Err(ModalAcousticTimeError::BudgetExceeded { .. })
        ));
    }
    let mut polls = 0;
    let result = PreparedModalTones::try_new_with_damping_rates_controlled(
        &[1.0, 2.0],
        &[0.5, 0.7],
        &[3.0, 4.0],
        0.1,
        budget,
        true,
        || {
            polls += 1;
            polls == 9
        },
    );
    assert!(matches!(
        result,
        Err(ModalAcousticTimeError::InvalidInput {
            what: "physical tone preparation cancelled"
        })
    ));
    assert_eq!(polls, 9);
}
