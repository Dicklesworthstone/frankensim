//! Original affine envelopes, exact clock variations and coherent work/heat.
use super::*;

fn affine_bank(natural: f64, damping: f64, frequencies: &[f64], h: f64) -> PreparedModalTones {
    PreparedModalTones::try_new_affine_with_damping_rates(
        &[natural],
        &[damping],
        frequencies,
        h,
        ModalToneBudget::realtime_reference(),
    )
    .unwrap()
}

// Independently integrate the ORIGINAL physical ODE and each frequency
// variational equation. The envelope is interpolated in the start carrier basis.
fn ode(
    natural: f64,
    damping: f64,
    frequencies: [f64; 3],
    h: f64,
    initial: ModalAcousticState,
    endpoints: [[f64; 3]; 4],
) -> [f64; 16] {
    let rhs = |t: f64, x: [f64; 16]| {
        let u = t / h;
        let mut force = 0.0;
        let mut partial = [0.0; 3];
        for j in 0..3 {
            let (s, c) = (frequencies[j] * t).sin_cos();
            let fc = (1.0 - u) * endpoints[0][j] + u * endpoints[2][j];
            let fs = (1.0 - u) * endpoints[1][j] + u * endpoints[3][j];
            force += fc * c + fs * s;
            partial[j] = t * (fs * c - fc * s);
        }
        let mut out = [0.0; 16];
        out[..4].copy_from_slice(&[
            x[1],
            force - natural * natural * x[0] - damping * x[1],
            force * x[1],
            damping * x[1] * x[1],
        ]);
        for j in 0..3 {
            let i = 4 + 4 * j;
            out[i] = x[i + 1];
            out[i + 1] = partial[j] - natural * natural * x[i] - damping * x[i + 1];
            out[i + 2] = partial[j] * x[1] + force * x[i + 1];
            out[i + 3] = 2.0 * damping * x[1] * x[i + 1];
        }
        out
    };
    let add = |x: [f64; 16], y: [f64; 16], a: f64| core::array::from_fn(|i| x[i] + a * y[i]);
    let mut x = [0.0; 16];
    x[0] = initial.displacement_m_sqrt_kg;
    x[1] = initial.velocity_m_sqrt_kg_per_s;
    let dt = h / 8192.0;
    for j in 0..8192 {
        let t = dt * j as f64;
        let a = rhs(t, x);
        let b = rhs(t + dt * 0.5, add(x, a, dt * 0.5));
        let c = rhs(t + dt * 0.5, add(x, b, dt * 0.5));
        let d = rhs(t + dt, add(x, c, dt));
        for i in 0..16 {
            x[i] += dt / 6.0 * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]);
        }
    }
    x
}

#[test]
fn affine_clock_power_matches_original_variational_odes_and_constant_reduction() {
    let endpoints = [
        [0.7, -0.4, 0.2],
        [-0.2, 0.6, 0.5],
        [-0.3, 0.8, 0.4],
        [0.9, -0.5, 0.1],
    ];
    let initial = state(0.3, -0.2);
    for (natural, damping, frequencies, h) in [
        (0.0, 0.0, [0.0, 1.3, 2.1], 0.7),
        (0.0, 3.0, [0.0, 1.3, 2.1], 0.7),
        (2.0, 0.0, [2.0, 2.0001, 3.3], 2.3),
        (2.0, 4.0, [0.0, 2.0, 2.0], 1.7),
        (2.0, 20.0, [0.001, 1.3, 2.1], 1.7),
        (20.0, 0.03, [20.0, 20.0001, 32.0], 0.7),
        (2.0, 0.4, [0.0, 1e-10, 1e-10], 0.7),
    ] {
        let bank = affine_bank(natural, damping, &frequencies, h);
        let prepared = bank
            .prepare_power_frequency_gradients(
                ModalToneBudget::realtime_reference(),
                ModalTonePowerGradientBudget::realtime_reference(),
                || false,
            )
            .unwrap();
        let result = prepared
            .gradients_mode_affine(
                0,
                initial,
                &endpoints[0],
                &endpoints[1],
                &endpoints[2],
                &endpoints[3],
                || false,
            )
            .unwrap();
        let reference = ode(natural, damping, frequencies, h, initial, endpoints);
        for (j, gradient) in result.source_derivatives().iter().enumerate() {
            close(gradient.input_work, reference[6 + 4 * j], 3e-10);
            close(gradient.viscous_dissipation, reference[7 + 4 * j], 3e-10);
        }
        let constant = super::bank(&[natural], &[damping], &frequencies, h);
        let plain = constant
            .prepare_power_frequency_gradients(
                ModalToneBudget::realtime_reference(),
                ModalTonePowerGradientBudget::realtime_reference(),
                || false,
            )
            .unwrap()
            .gradients_mode(0, initial, &endpoints[0], &endpoints[1], || false)
            .unwrap();
        let reduced = prepared
            .gradients_mode_affine(
                0,
                initial,
                &endpoints[0],
                &endpoints[1],
                &endpoints[0],
                &endpoints[1],
                || false,
            )
            .unwrap();
        for (a, b) in reduced
            .source_derivatives()
            .iter()
            .zip(plain.source_derivatives())
        {
            close(a.input_work, b.input_work, 3e-12);
            close(a.viscous_dissipation, b.viscous_dissipation, 3e-12);
        }
    }
}

#[test]
fn independent_affine_duplicates_survive_overflowing_individual_deltas() {
    let h = 0.1;
    let initial = state(0.3, 0.2);
    let bank = affine_bank(0.0, 0.0, &[0.0, 0.0], h);
    let prepared = bank
        .prepare_power_frequency_gradients(
            ModalToneBudget::realtime_reference(),
            ModalTonePowerGradientBudget::realtime_reference(),
            || false,
        )
        .unwrap();
    let magnitude: f64 = 1e308;
    assert!((magnitude - (-magnitude)).is_infinite());
    let result = prepared
        .gradients_mode_affine(
            0,
            initial,
            &[0.0; 2],
            &[-magnitude, magnitude],
            &[0.0; 2],
            &[magnitude, -magnitude],
            || false,
        )
        .unwrap();
    // Nominal force is zero, v is constant. Integral v*t*S(t) is exact.
    let expected = (magnitude / 6.0) * h * h * initial.velocity_m_sqrt_kg_per_s;
    relative(result.source_derivatives()[0].input_work, expected, 3e-13);
    relative(result.source_derivatives()[1].input_work, -expected, 3e-13);
    assert_eq!(result.source_derivatives()[0].viscous_dissipation, 0.0);
    assert_eq!(result.source_derivatives()[1].viscous_dissipation, 0.0);
    assert_eq!(
        result.source_derivatives()[0].input_work + result.source_derivatives()[1].input_work,
        0.0
    );
}

#[test]
fn long_decay_keeps_the_complementary_start_clock_work() {
    let h = 1e6;
    let bank = affine_bank(2.0, 4.0, &[0.0], h);
    let prepared = bank
        .prepare_power_frequency_gradients(
            ModalToneBudget::realtime_reference(),
            ModalTonePowerGradientBudget::realtime_reference(),
            || false,
        )
        .unwrap();
    let result = prepared
        .gradients_mode_affine(0, state(1.0, 0.0), &[4.0], &[1.0], &[4.0], &[0.0], || false)
        .unwrap();
    // The nominal state is the exact stationary equilibrium. The clock force
    // is t*(1-t/h), so W_Omega=4*q_Omega(h). Its polynomial particular gives
    // q_Omega(h)=1/4-3/(8h); all initial transients underflow at this interval.
    relative(
        result.source_derivatives()[0].input_work,
        1.0 - 1.5 / h,
        3e-8,
    );
    assert_eq!(result.source_derivatives()[0].viscous_dissipation, 0.0);
}

#[test]
fn full_quadratic_dimension_units_budgets_and_reuse_are_admitted() {
    let frequencies: Vec<_> = (0..32).map(|j| 0.4 + 0.03 * j as f64).collect();
    let h = 0.2;
    let bank = affine_bank(2.0, 0.3, &frequencies, h);
    let cold = ModalToneBudget::realtime_reference();
    let query = ModalTonePowerGradientBudget::realtime_reference();
    assert!(
        bank.prepare_power_frequency_gradients(
            ModalToneBudget {
                maximum_storage_bytes: 1_000_000,
                ..cold
            },
            query,
            || false,
        )
        .is_err()
    );
    assert!(
        bank.prepare_power_frequency_gradients(
            cold,
            ModalTonePowerGradientBudget {
                maximum_query_work: 1_000_000,
                ..query
            },
            || false,
        )
        .is_err()
    );
    assert!(
        bank.prepare_power_frequency_gradients(
            cold,
            ModalTonePowerGradientBudget {
                maximum_scratch_bytes: 1,
                ..query
            },
            || false,
        )
        .is_err()
    );
    let prepared = bank
        .prepare_power_frequency_gradients(cold, query, || false)
        .unwrap();
    let initial = state(0.3, -0.2);
    let (mut c0, mut s0, mut c1, mut s1) = ([0.0; 32], [0.0; 32], [0.0; 32], [0.0; 32]);
    c0[31] = 0.7;
    s0[31] = -0.2;
    c1[31] = -0.3;
    s1[31] = 0.9;
    let result = prepared
        .gradients_mode_affine(0, initial, &c0, &s0, &c1, &s1, || false)
        .unwrap();
    let single = affine_bank(2.0, 0.3, &[frequencies[31]], h);
    let reference = single
        .prepare_power_frequency_gradients(cold, query, || false)
        .unwrap()
        .gradients_mode_affine(0, initial, &[0.7], &[-0.2], &[-0.3], &[0.9], || false)
        .unwrap();
    close(
        result.source_derivatives()[31].input_work,
        reference.source_derivatives()[0].input_work,
        2e-12,
    );
    close(
        result.source_derivatives()[31].viscous_dissipation,
        reference.source_derivatives()[0].viscous_dissipation,
        2e-12,
    );
    assert!(
        result.source_derivatives()[..31]
            .iter()
            .all(|g| g.input_work == 0.0 && g.viscous_dissipation == 0.0)
    );
    let mut polls = 0;
    assert!(
        prepared
            .gradients_mode_affine(0, initial, &c0, &s0, &c1, &s1, || {
                polls += 1;
                polls > 5
            })
            .is_err()
    );
    assert!(
        prepared
            .gradients_mode_affine(0, initial, &c0[..31], &s0, &c1, &s1, || false)
            .is_err()
    );
    assert_eq!(
        result,
        prepared
            .gradients_mode_affine(0, initial, &c0, &s0, &c1, &s1, || false)
            .unwrap()
    );
    // Physical coordinate/time conversion, not only an admission-scale change.
    let (a, b) = (1e3, 1e-2);
    let converted = affine_bank(2.0 / b, 0.3 / b, &[frequencies[31] / b], h * b);
    let scaled = converted
        .prepare_power_frequency_gradients(cold, query, || false)
        .unwrap()
        .gradients_mode_affine(
            0,
            state(0.3 * a, -0.2 * a / b),
            &[0.7 * a / b / b],
            &[-0.2 * a / b / b],
            &[-0.3 * a / b / b],
            &[0.9 * a / b / b],
            || false,
        )
        .unwrap();
    close(
        scaled.source_derivatives()[0].input_work * b / a / a,
        reference.source_derivatives()[0].input_work,
        3e-13,
    );
    close(
        scaled.source_derivatives()[0].viscous_dissipation * b / a / a,
        reference.source_derivatives()[0].viscous_dissipation,
        3e-13,
    );
}
