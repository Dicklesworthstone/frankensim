//! Authoritative viscous-rate banks, including exactly free damped coordinates.

use fs_math::c64::C64;
use fs_modal::time::{
    LinearForceTransition, ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget,
    ModalAcousticTimeModel, ModalAcousticWorkspace, ModalToneBudget, PreparedHarmonicForce,
    PreparedModalHarmonic, PreparedModalTones, advance_exact_zoh,
};

fn state(q: f64, v: f64) -> ModalAcousticState {
    ModalAcousticState {
        displacement_m_sqrt_kg: q,
        velocity_m_sqrt_kg_per_s: v,
    }
}

fn mode(omega: f64, ratio: f64) -> ModalAcousticMode {
    ModalAcousticMode {
        angular_frequency_rad_s: omega,
        damping_ratio: ratio,
        pressure_per_modal_velocity: C64::ZERO,
    }
}

fn budget() -> ModalAcousticTimeBudget {
    ModalAcousticTimeBudget {
        maximum_abs_displacement_m_sqrt_kg: 100.0,
        maximum_total_energy_j: 1.0e20,
        ..ModalAcousticTimeBudget::audible_reference()
    }
}

fn bank(sample_rate: u32, omega: &[f64], rates: &[f64]) -> ModalAcousticTimeModel {
    ModalAcousticTimeModel::try_new_with_damping_rates(
        sample_rate,
        omega.iter().map(|&w| mode(w, 0.0)).collect(),
        rates,
        budget(),
    )
    .unwrap()
}

fn near(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 3.0e-11 * scale.max(1.0e-300),
        "actual={actual:.17e}, expected={expected:.17e}, scale={scale:.17e}"
    );
}

#[test]
fn physical_metadata_is_complete_and_unrepresentable_ratios_fail_visibly() {
    let rates = [2.5, 4.0, f64::from_bits(1), 0.0, 0.25];
    let model = bank(48_000, &[0.0, 1.0e-310, 1.0e3, 0.0, 4.0], &rates);
    assert_eq!(model.damping_rates_per_s(), rates);
    for (i, &rate) in rates.iter().enumerate() {
        assert_eq!(model.damping_rate_per_s(i), rate);
    }
    // Undefined, overflowing and underflowing positive ratios all remain
    // distinguishable from physical zero loss through the legacy metadata.
    for m in &model.modes()[..3] {
        assert_eq!(m.damping_ratio, f64::INFINITY);
    }
    assert_eq!(model.modes()[3].damping_ratio, 0.0);
    assert_eq!(model.modes()[4].damping_ratio, 0.03125);
    let wrong_legacy = advance_exact_zoh(model.modes()[0], state(0.0, 1.0), 0.0, 0.1);
    assert!(!wrong_legacy.velocity_m_sqrt_kg_per_s.is_finite());
    let correct = model
        .advance_exact_zoh_duration_step(0, state(0.0, 1.0), 0.0, 0.1)
        .unwrap();
    near(correct.velocity_m_sqrt_kg_per_s, (-0.25_f64).exp(), 1.0);

    let copied = ModalAcousticTimeModel::try_new_with_damping_rates(
        48_000,
        model.modes().to_vec(),
        model.damping_rates_per_s(),
        budget(),
    )
    .unwrap();
    assert_eq!(copied.modes(), model.modes());
    assert_eq!(copied.damping_rates_per_s(), model.damping_rates_per_s());
    assert_eq!(
        model.without_pressure_transfers().damping_rates_per_s(),
        rates
    );
}

#[test]
fn exact_free_drag_zoh_retains_motion_work_heat_and_port_mobility_at_every_duration() {
    let rate = 2.5_f64;
    let initial = state(-0.2, 0.3);
    let force = 1.2;
    let mut source = bank(4, &[0.0], &[rate]);
    source.restore_states(&[initial]).unwrap();
    for h in [0.25, 0.125, 0.7] {
        let mut model = source.clone();
        let mut workspace = ModalAcousticWorkspace::new(&model);
        let terminal = force / rate;
        let tail = initial.velocity_m_sqrt_kg_per_s - terminal;
        let i0 = -(-rate * h).exp_m1() / rate;
        let i2 = -(-2.0 * rate * h).exp_m1() / (2.0 * rate);
        let distance = terminal * h + tail * i0;
        let expected = state(
            initial.displacement_m_sqrt_kg + distance,
            terminal + tail * (-rate * h).exp(),
        );
        let heat = rate * (terminal * terminal * h + 2.0 * terminal * tail * i0 + tail * tail * i2);
        let preview = model
            .advance_exact_zoh_duration_step(0, initial, force, h)
            .unwrap();
        let frame = model
            .step_duration_into(&[force], h, &mut workspace)
            .unwrap();
        near(
            frame.input_work_j,
            force * distance,
            (force * distance).abs(),
        );
        near(frame.viscous_dissipation_j, heat, heat);
        assert_eq!(
            frame.viscous_dissipation_j,
            frame.modal_viscous_losses[0].energy_j
        );
        assert_eq!(model.states()[0], preview);
        near(
            preview.displacement_m_sqrt_kg,
            expected.displacement_m_sqrt_kg,
            1.0,
        );
        near(
            preview.velocity_m_sqrt_kg_per_s,
            expected.velocity_m_sqrt_kg_per_s,
            1.0,
        );
        let target = -0.13;
        let port_force = source
            .held_force_for_port_velocity(&[0.4], target, h)
            .unwrap();
        let endpoint = source
            .advance_exact_zoh_duration_step(0, initial, 0.4 * port_force, h)
            .unwrap();
        near(0.4 * endpoint.velocity_m_sqrt_kg_per_s, target, 1.0);
    }
    assert_eq!(source.states(), &[initial]);
}

#[test]
fn mixed_free_soft_and_elastic_ramps_route_every_physical_loss_once() {
    let frequencies = [0.0, 1.0e-310, 4.0];
    let rates = [2.5, 4.0, 0.7];
    let initial = [state(-0.2, 0.3), state(0.1, -0.2), state(-0.02, 0.4)];
    let start = [1.2, -0.4, 0.3];
    let end = [-0.4, 0.7, 1.1];
    let mut model = bank(8, &frequencies, &rates);
    model.restore_states(&initial).unwrap();
    model.prepare_linear_forcing().unwrap();
    let expected = (0..3)
        .map(|i| {
            LinearForceTransition::try_new_with_damping_rate(frequencies[i], rates[i], 0.125)
                .unwrap()
                .advance_with_energy(initial[i], start[i], end[i])
                .unwrap()
        })
        .collect::<Vec<_>>();
    let mut workspace = ModalAcousticWorkspace::new(&model);
    let frame = model
        .step_linear_into(&start, &end, &mut workspace)
        .unwrap();
    assert_eq!(
        frame.input_work_j,
        expected.iter().map(|f| f.input_work_j).sum::<f64>()
    );
    assert_eq!(
        frame.viscous_dissipation_j,
        expected
            .iter()
            .map(|f| f.viscous_dissipation_j)
            .sum::<f64>()
    );
    for (loss, expected) in frame.modal_viscous_losses.iter().zip(&expected) {
        assert_eq!(loss.energy_j, expected.viscous_dissipation_j);
        assert!(loss.energy_j > 0.0);
    }
    for (actual, expected) in model.states().iter().zip(expected) {
        assert_eq!(*actual, expected.state);
    }
}

#[test]
fn free_decay_retains_tiny_heat_and_stiff_drag_without_a_loss_floor() {
    let h = 0.25;
    let initial = state(0.0, 0.7);
    for rate in [1.0e-18_f64, 1.0e-5, 2.5, 1.0e6] {
        let mut model = bank(4, &[0.0], &[rate]);
        model.restore_states(&[initial]).unwrap();
        let mut workspace = ModalAcousticWorkspace::new(&model);
        let frame = model.step_duration_into(&[0.0], h, &mut workspace).unwrap();
        let heat = -0.5 * initial.velocity_m_sqrt_kg_per_s.powi(2) * (-2.0 * rate * h).exp_m1();
        near(frame.viscous_dissipation_j, heat, heat);
        assert!(frame.viscous_dissipation_j > 0.0);
        assert_eq!(frame.input_work_j, 0.0);
        let q = -initial.velocity_m_sqrt_kg_per_s * (-rate * h).exp_m1() / rate;
        let v = initial.velocity_m_sqrt_kg_per_s * (-rate * h).exp();
        near(model.states()[0].displacement_m_sqrt_kg, q, q.abs());
        near(model.states()[0].velocity_m_sqrt_kg_per_s, v, v.abs());
        if rate == 1.0e-18 {
            assert_eq!(
                model.states()[0].velocity_m_sqrt_kg_per_s,
                initial.velocity_m_sqrt_kg_per_s
            );
        }
    }
}

#[test]
fn harmonic_bank_uses_the_same_rate_motion_work_and_heat_as_the_scalar_owner() {
    let (rate, drive, h) = (2.5, 3.0, 0.7);
    let initial = state(-0.2, 0.3);
    let mut model = bank(48_000, &[0.0], &[rate]);
    model.restore_states(&[initial]).unwrap();
    let prepared = PreparedModalHarmonic::try_new(&model, drive, h).unwrap();
    let expected = PreparedHarmonicForce::try_new_with_damping_rate(0.0, rate, drive, h)
        .unwrap()
        .advance_with_energy(initial, 1.2, -0.4)
        .unwrap();
    let mut workspace = ModalAcousticWorkspace::new(&model);
    let frame = model
        .step_harmonic_into(&prepared, &[1.2], &[-0.4], &mut workspace)
        .unwrap();
    assert_eq!(frame.input_work_j, expected.input_work_j);
    assert_eq!(frame.viscous_dissipation_j, expected.viscous_dissipation_j);
    assert_eq!(
        frame.modal_viscous_losses[0].roundoff_tolerance_j,
        expected.dissipation_roundoff_tolerance_j
    );
    assert_eq!(model.states(), &[expected.state]);
}

// Independently integrate q'=v, v'=F(t)-c*v and the two physical powers.
// No production transition or power quadrature enters this RK4 reference.
fn free_reference(
    rate: f64,
    initial: ModalAcousticState,
    h: f64,
    force: impl Fn(f64) -> f64,
) -> [f64; 4] {
    let mut x = [
        initial.displacement_m_sqrt_kg,
        initial.velocity_m_sqrt_kg_per_s,
        0.0,
        0.0,
    ];
    let dt = h / 8192.0;
    let derivative = |t: f64, y: [f64; 4]| {
        let f = force(t);
        [y[1], f - rate * y[1], f * y[1], rate * y[1] * y[1]]
    };
    for step in 0..8192 {
        let t = dt * step as f64;
        let a = derivative(t, x);
        let b = derivative(
            t + 0.5 * dt,
            core::array::from_fn(|i| x[i] + 0.5 * dt * a[i]),
        );
        let c = derivative(
            t + 0.5 * dt,
            core::array::from_fn(|i| x[i] + 0.5 * dt * b[i]),
        );
        let d = derivative(t + dt, core::array::from_fn(|i| x[i] + dt * c[i]));
        for i in 0..4 {
            x[i] += dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    x
}

#[test]
fn simultaneous_and_affine_tones_keep_free_drag_cross_terms_and_power() {
    let (rate, h) = (2.5, 0.7);
    let initial = state(-0.2, 0.3);
    let frequencies = [0.0, 3.0, 7.25, 3.0];
    let c0 = [0.4, 1.2, -0.3, 0.1];
    let s0 = [17.0, -0.4, 0.2, 0.05];
    let c1 = [-0.2, 0.3, 0.8, -0.05];
    let s1 = [31.0, 0.6, -0.1, 0.2];
    for affine in [false, true] {
        let mut model = bank(48_000, &[0.0], &[rate]);
        model.restore_states(&[initial]).unwrap();
        let prepared = if affine {
            PreparedModalTones::try_new_affine(
                &model,
                &frequencies,
                h,
                ModalToneBudget::realtime_reference(),
            )
        } else {
            PreparedModalTones::try_new(
                &model,
                &frequencies,
                h,
                ModalToneBudget::realtime_reference(),
            )
        }
        .unwrap();
        let expected = free_reference(rate, initial, h, |t| {
            (0..frequencies.len())
                .map(|i| {
                    let u = if affine { t / h } else { 0.0 };
                    let (s, c) = (frequencies[i] * t).sin_cos();
                    (c0[i] + u * (c1[i] - c0[i])) * c + (s0[i] + u * (s1[i] - s0[i])) * s
                })
                .sum()
        });
        let mut workspace = ModalAcousticWorkspace::new(&model);
        let frame = if affine {
            model.step_tones_affine_into(&prepared, &c0, &s0, &c1, &s1, &mut workspace)
        } else {
            model.step_tones_into(&prepared, &c0, &s0, &mut workspace)
        }
        .unwrap();
        near(frame.input_work_j, expected[2], expected[2].abs());
        near(frame.viscous_dissipation_j, expected[3], expected[3]);
        assert!(frame.viscous_dissipation_j > 0.0);
        near(
            model.states()[0].displacement_m_sqrt_kg,
            expected[0],
            expected[0].abs(),
        );
        near(
            model.states()[0].velocity_m_sqrt_kg_per_s,
            expected[1],
            expected[1].abs(),
        );
    }
}

#[test]
fn physical_rate_updates_preserve_state_and_invalidate_every_changed_cache() {
    let initial = state(-0.2, 0.3);
    let mut source = bank(4, &[0.0], &[2.5]);
    source.restore_states(&[initial]).unwrap();
    source.prepare_linear_forcing().unwrap();
    let harmonic = PreparedModalHarmonic::try_new(&source, 3.0, 0.25).unwrap();
    let tones =
        PreparedModalTones::try_new(&source, &[3.0], 0.25, ModalToneBudget::realtime_reference())
            .unwrap();
    let mut same = source.with_damping_rates(&[2.5]).unwrap();
    let mut workspace = ModalAcousticWorkspace::new(&same);
    same.step_harmonic_into(&harmonic, &[1.2], &[-0.4], &mut workspace)
        .unwrap();
    let mut changed = source.with_damping_rates(&[5.0]).unwrap();
    // Both legacy diagnostic ratios are infinite, so checking modes alone
    // cannot identify this physical coefficient change.
    assert_eq!(changed.modes(), source.modes());
    assert_eq!(changed.states(), source.states());
    assert_eq!(changed.damping_rates_per_s(), &[5.0]);
    assert!(
        changed
            .step_harmonic_into(&harmonic, &[1.2], &[-0.4], &mut workspace)
            .is_err()
    );
    assert!(
        changed
            .step_tones_into(&tones, &[1.2], &[-0.4], &mut workspace)
            .is_err()
    );
    assert!(
        changed
            .step_linear_into(&[1.2], &[-0.4], &mut workspace)
            .is_err()
    );
    assert_eq!(changed.states(), &[initial]);
    changed.prepare_linear_forcing().unwrap();
    changed
        .step_linear_into(&[1.2], &[-0.4], &mut workspace)
        .unwrap();
    let reset = source.with_damping_ratios(&[0.0]).unwrap();
    assert_eq!(reset.damping_rates_per_s(), &[0.0]);
    assert_eq!(reset.modes()[0].damping_ratio, 0.0);
    assert_eq!(reset.states(), source.states());
    assert_eq!(source.damping_rates_per_s(), &[2.5]);
}

#[test]
fn rate_admission_and_rejected_steps_do_not_mutate_accepted_state() {
    let mut model = bank(4, &[0.0], &[2.5]);
    let initial = state(-0.2, 0.3);
    model.restore_states(&[initial]).unwrap();
    for rates in [
        vec![],
        vec![-1.0],
        vec![f64::NAN],
        vec![f64::INFINITY],
        vec![1.0, 2.0],
    ] {
        assert!(
            ModalAcousticTimeModel::try_new_with_damping_rates(
                4,
                vec![mode(0.0, 0.0)],
                &rates,
                budget()
            )
            .is_err()
        );
        assert!(model.with_damping_rates(&rates).is_err());
        assert_eq!(model.states(), &[initial]);
    }
    let mut radiating_free = mode(0.0, 0.0);
    radiating_free.pressure_per_modal_velocity = C64::from_re(1.0);
    assert!(
        ModalAcousticTimeModel::try_new_with_damping_rates(
            4,
            vec![radiating_free],
            &[2.5],
            budget()
        )
        .is_err()
    );
    let mut workspace = ModalAcousticWorkspace::new(&model);
    assert!(
        model
            .step_duration_into(&[f64::MAX], 0.7, &mut workspace)
            .is_err()
    );
    assert!(
        model
            .advance_exact_zoh_duration_step(1, initial, 0.0, 0.1)
            .is_err()
    );
    assert!(
        model
            .advance_exact_zoh_duration_step(0, initial, 0.0, 0.0)
            .is_err()
    );
    assert_eq!(model.states(), &[initial]);
}

#[test]
fn legacy_ratio_zoh_keeps_its_exact_arithmetic_and_noop_rate_updates() {
    let initial = state(0.1, -0.2);
    for ratio in [0.0, 1.0e-12, 0.125, 1.0, 3.0, 1.0e12] {
        let original_mode = mode(4.0, ratio);
        let mut model =
            ModalAcousticTimeModel::try_new(48_000, vec![original_mode], budget()).unwrap();
        model.restore_states(&[initial]).unwrap();
        let same = model
            .with_damping_rates(model.damping_rates_per_s())
            .unwrap();
        for h in [1.0 / 48_000.0, 0.5 / 48_000.0, 0.017] {
            let expected = advance_exact_zoh(original_mode, initial, 1.2, h);
            assert_eq!(
                model
                    .advance_exact_zoh_duration_step(0, initial, 1.2, h)
                    .unwrap(),
                expected
            );
            assert_eq!(
                same.advance_exact_zoh_duration_step(0, initial, 1.2, h)
                    .unwrap(),
                expected
            );
        }
    }
}

#[test]
fn ratio_updates_keep_authoritative_rates_finite_and_refuse_transactionally() {
    let mut model =
        ModalAcousticTimeModel::try_new(48_000, vec![mode(4.0, 0.125)], budget()).unwrap();
    let initial = state(0.1, -0.2);
    model.restore_states(&[initial]).unwrap();
    assert!(model.with_damping_ratios(&[f64::MAX]).is_err());
    assert_eq!(model.damping_rates_per_s(), &[1.0]);
    assert_eq!(model.states(), &[initial]);
    assert!(ModalAcousticTimeModel::try_new(48_000, vec![mode(4.0, f64::MAX)], budget()).is_err());

    // The physical coefficient is finite even though evaluating 2*zeta first
    // overflows. Ratio admission and updates use the same guarded product.
    let soft =
        ModalAcousticTimeModel::try_new(48_000, vec![mode(1.0e-100, 0.0)], budget()).unwrap();
    let changed = soft.with_damping_ratios(&[1.0e308]).unwrap();
    assert!(changed.damping_rate_per_s(0).is_finite());
    assert_eq!(changed.damping_rate_per_s(0), 2.0 * (1.0e308 * 1.0e-100));
    assert_eq!(soft.damping_rates_per_s(), &[0.0]);
}
