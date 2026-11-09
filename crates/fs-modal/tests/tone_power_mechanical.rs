//! G1/G3: coherent mechanical work/heat tangents and physical parameter chains.

use fs_modal::time::{
    ModalAcousticState, ModalToneBudget, ModalTonePowerGradientBudget,
    ModalTonePowerMechanicalGradient, ModalToneSensitivityBudget, PreparedHarmonicForce,
    PreparedModalTones,
};

fn state(q: f64, v: f64) -> ModalAcousticState {
    ModalAcousticState {
        displacement_m_sqrt_kg: q,
        velocity_m_sqrt_kg_per_s: v,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual:.17e} != {expected:.17e}"
    );
}

fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        actual.is_finite()
            && expected != 0.0
            && (actual - expected).abs() <= tolerance * expected.abs(),
        "{actual:.17e} != {expected:.17e}"
    );
}

fn bank(
    natural: f64,
    damping: f64,
    frequencies: &[f64],
    h: f64,
    affine: bool,
) -> PreparedModalTones {
    let create = if affine {
        PreparedModalTones::try_new_affine_with_damping_rates
    } else {
        PreparedModalTones::try_new_with_damping_rates
    };
    create(
        &[natural],
        &[damping],
        frequencies,
        h,
        ModalToneBudget::realtime_reference(),
    )
    .unwrap()
}

fn motion_budget() -> ModalToneSensitivityBudget {
    ModalToneSensitivityBudget {
        maximum_storage_bytes: 8 * 1024 * 1024,
        maximum_preparation_work: 100_000_000,
        maximum_query_work: 100_000,
    }
}

fn entries(g: ModalTonePowerMechanicalGradient) -> [f64; 4] {
    [
        g.stiffness.input_work,
        g.stiffness.viscous_dissipation,
        g.damping_rate.input_work,
        g.damping_rate.viscous_dissipation,
    ]
}

// Independent physical RK4: [q-q0,v,W,Q,q_k,v_k,W_k,Q_k,q_c,v_c,W_c,Q_c].
// The variational acceleration sources are -q and -v. The heat c derivative
// includes v^2 explicitly even when c=0. cos-1 and q-q0 preserve a tiny DC
// perturbation of an initially balanced equilibrium in the reference too.
#[allow(clippy::too_many_arguments)]
fn reference(
    natural: f64,
    damping: f64,
    frequencies: &[f64],
    h: f64,
    initial: ModalAcousticState,
    knots: [&[f64]; 4],
) -> [f64; 12] {
    let [c0, s0, c1, s1] = knots;
    let q0 = initial.displacement_m_sqrt_kg;
    let k = natural * natural;
    let held = c0.iter().sum::<f64>();
    let imbalance = held - k * q0;
    let rhs = |t: f64, x: [f64; 12]| {
        let u = t / h;
        let mut increment = 0.0;
        for j in 0..frequencies.len() {
            let dc = c1[j] - c0[j];
            let cosine = c0[j] + u * dc;
            let sine = s0[j] + u * (s1[j] - s0[j]);
            let angle = frequencies[j] * t;
            increment += u * dc - 2.0 * cosine * (0.5 * angle).sin().powi(2) + sine * angle.sin();
        }
        let force = held + increment;
        [
            x[1],
            imbalance + increment - k * x[0] - damping * x[1],
            force * x[1],
            damping * x[1] * x[1],
            x[5],
            -(q0 + x[0]) - k * x[4] - damping * x[5],
            force * x[5],
            2.0 * damping * x[1] * x[5],
            x[9],
            -x[1] - k * x[8] - damping * x[9],
            force * x[9],
            x[1] * x[1] + 2.0 * damping * x[1] * x[9],
        ]
    };
    let add = |x: [f64; 12], y: [f64; 12], a: f64| core::array::from_fn(|i| x[i] + a * y[i]);
    let mut x = [0.0; 12];
    x[1] = initial.velocity_m_sqrt_kg_per_s;
    let steps = 16384;
    let dt = h / steps as f64;
    for step in 0..steps {
        let t = step as f64 * dt;
        let a = rhs(t, x);
        let b = rhs(t + dt * 0.5, add(x, a, dt * 0.5));
        let d = rhs(t + dt * 0.5, add(x, b, dt * 0.5));
        let e = rhs(t + dt, add(x, d, dt));
        for i in 0..12 {
            x[i] += (dt / 6.0) * (a[i] + 2.0 * b[i] + 2.0 * d[i] + e[i]);
        }
    }
    x[0] += q0;
    x
}

#[test]
fn constant_and_affine_mechanics_match_independent_power_variational_odes() {
    let (c0, s0) = ([0.7, -0.4, 0.2], [-0.2, 0.6, 0.5]);
    let initial = state(0.3, -0.2);
    for affine in [false, true] {
        let (c1, s1) = if affine {
            ([-0.2, 0.4, 0.8], [0.2, -0.3, 0.1])
        } else {
            (c0, s0)
        };
        for (natural, damping, frequencies, h) in [
            (0.0, 0.0, [0.0, 1.3, 2.1], 0.7),
            (0.0, 3.0, [0.0, 1.3, 2.1], 0.7),
            (2.0, 0.0, [2.0, 2.0001, 3.3], 2.3),
            (2.0, 4.0, [0.0, 2.0, 2.0], 1.7),
            (2.0, 20.0, [0.001, 1.3, 2.1], 1.7),
            (20.0, 0.03, [20.0, 20.0001, 32.0], 0.7),
            (2.0, 0.4, [0.0, 1e-10, 1e-10], 0.7),
        ] {
            let bank = bank(natural, damping, &frequencies, h, affine);
            let power = bank
                .prepare_power_mechanical_gradients(
                    ModalToneBudget::realtime_reference(),
                    ModalTonePowerGradientBudget::realtime_reference(),
                    || false,
                )
                .unwrap();
            let query = || {
                if affine {
                    power.gradients_mode_affine(0, initial, &c0, &s0, &c1, &s1, || false)
                } else {
                    power.gradients_mode(0, initial, &c0, &s0, || false)
                }
                .unwrap()
            };
            let actual = query();
            let expected = reference(
                natural,
                damping,
                &frequencies,
                h,
                initial,
                [&c0, &s0, &c1, &s1],
            );
            for (a, b) in entries(actual).into_iter().zip([
                expected[6],
                expected[7],
                expected[10],
                expected[11],
            ]) {
                close(a, b, 5e-10);
            }
            assert_eq!(actual, query());
            if damping == 0.0 {
                assert_eq!(actual.stiffness.viscous_dissipation, 0.0);
            }
            // Independent energy bookkeeping: test the resulting derivative,
            // not a production path that defines heat by this subtraction.
            let motion = bank
                .prepare_sensitivities(motion_budget(), || false)
                .unwrap();
            let frame = if affine {
                motion.advance_mode_affine(0, initial, &c0, &s0, &c1, &s1, || false)
            } else {
                motion.advance_mode(0, initial, &c0, &s0, || false)
            }
            .unwrap();
            let q = frame.state.displacement_m_sqrt_kg;
            let v = frame.state.velocity_m_sqrt_kg_per_s;
            let p = frame.parameters;
            close(
                actual.stiffness.input_work - actual.stiffness.viscous_dissipation,
                natural * natural * q * p.displacement_by_stiffness
                    + v * p.velocity_by_stiffness
                    + 0.5 * (q * q - initial.displacement_m_sqrt_kg.powi(2)),
                2e-11,
            );
            close(
                actual.damping_rate.input_work - actual.damping_rate.viscous_dissipation,
                natural * natural * q * p.displacement_by_damping_rate
                    + v * p.velocity_by_damping_rate,
                2e-11,
            );
        }
    }
}

#[test]
fn static_dc_force_only_work_and_zero_drag_heat_derivative_are_retained() {
    // A nominally motionless forced equilibrium has nonzero dW/dk. A QR
    // tangent that discarded the original force-only residual would miss it.
    for affine in [false, true] {
        for damping in [0.0, 4.0] {
            let h = 0.7;
            let bank = bank(2.0, damping, &[0.0, 0.0], h, affine);
            let prepared = bank
                .prepare_power_mechanical_gradients(
                    ModalToneBudget::realtime_reference(),
                    ModalTonePowerGradientBudget::realtime_reference(),
                    || false,
                )
                .unwrap();
            let actual = if affine {
                prepared.gradients_mode_affine(
                    0,
                    state(0.25, 0.0),
                    &[0.7, 0.3],
                    &[0.2, -0.4],
                    &[0.7, 0.3],
                    &[0.2, -0.4],
                    || false,
                )
            } else {
                prepared.gradients_mode(0, state(0.25, 0.0), &[0.7, 0.3], &[0.2, -0.4], || false)
            }
            .unwrap();
            let retention = if damping == 0.0 {
                (2.0 * h).cos()
            } else {
                (1.0 + 2.0 * h) * (-2.0 * h).exp()
            };
            relative(
                actual.stiffness.input_work,
                -0.0625 * (1.0 - retention),
                2e-13,
            );
            assert_eq!(actual.stiffness.viscous_dissipation, 0.0);
            assert_eq!(actual.damping_rate.input_work, 0.0);
            assert_eq!(actual.damping_rate.viscous_dissipation, 0.0);
        }
    }
    let (h, q, v, force) = (1.7_f64, 0.2_f64, -0.3_f64, 0.6_f64);
    let free = bank(0.0, 0.0, &[0.0, 0.0], h, false);
    let actual = free
        .prepare_power_mechanical_gradients(
            ModalToneBudget::realtime_reference(),
            ModalTonePowerGradientBudget::realtime_reference(),
            || false,
        )
        .unwrap()
        .gradients_mode(0, state(q, v), &[0.7, -0.1], &[0.8, -0.8], || false)
        .unwrap();
    close(
        actual.stiffness.input_work,
        -force * (q * h.powi(2) / 2.0 + v * h.powi(3) / 6.0 + force * h.powi(4) / 24.0),
        2e-14,
    );
    close(
        actual.damping_rate.input_work,
        -force * (v * h.powi(2) / 2.0 + force * h.powi(3) / 6.0),
        2e-14,
    );
    assert_eq!(actual.stiffness.viscous_dissipation, 0.0);
    relative(
        actual.damping_rate.viscous_dissipation,
        v * v * h + v * force * h * h + force * force * h.powi(3) / 3.0,
        2e-13,
    );
}

#[test]
fn tiny_near_dc_heat_and_exact_duplicate_cancellation_preserve_mechanics() {
    let initial = state(0.25, 0.0);
    for affine in [false, true] {
        let bank = bank(2.0, 0.7, &[1e-6, 1e-6], 0.7, affine);
        let prepared = bank
            .prepare_power_mechanical_gradients(
                ModalToneBudget::realtime_reference(),
                ModalTonePowerGradientBudget::realtime_reference(),
                || false,
            )
            .unwrap();
        let (c, s) = ([0.7, 0.3], [0.2, -0.2]);
        let actual = if affine {
            prepared.gradients_mode_affine(0, initial, &c, &s, &c, &s, || false)
        } else {
            prepared.gradients_mode(0, initial, &c, &s, || false)
        }
        .unwrap();
        let r = reference(
            2.0,
            0.7,
            &[1e-6],
            0.7,
            initial,
            [&[1.0], &[0.0], &[1.0], &[0.0]],
        );
        relative(actual.stiffness.viscous_dissipation, r[7], 2e-10);
        relative(actual.damping_rate.viscous_dissipation, r[11], 2e-10);
        // Both nominal sources cancel at exactly equal frequency; mechanical
        // derivatives must be those of the single free initial response.
        let zero = if affine {
            prepared.gradients_mode_affine(
                0,
                state(0.2, -0.3),
                &[0.7, -0.7],
                &[0.4, -0.4],
                &[-0.2, 0.2],
                &[-0.6, 0.6],
                || false,
            )
        } else {
            prepared.gradients_mode(0, state(0.2, -0.3), &[0.7, -0.7], &[0.4, -0.4], || false)
        }
        .unwrap();
        let scalar = PreparedHarmonicForce::try_new_with_damping_rate(2.0, 0.7, 0.0, 0.7)
            .unwrap()
            .prepare_power_sensitivities()
            .unwrap()
            .gradients(state(0.2, -0.3), 0.0, 0.0)
            .unwrap();
        assert_eq!(zero.stiffness.input_work, 0.0);
        assert_eq!(zero.damping_rate.input_work, 0.0);
        close(
            zero.stiffness.viscous_dissipation,
            scalar.viscous_dissipation.stiffness,
            2e-13,
        );
        close(
            zero.damping_rate.viscous_dissipation,
            scalar.viscous_dissipation.damping_rate,
            2e-13,
        );
    }
}

#[test]
fn single_tone_scalar_agreement_and_all_thirty_two_source_columns() {
    for (natural, damping, drive, h) in [
        (0.0, 3.0, 1.3, 0.7),
        (2.0, 0.0, 2.0, 2.3),
        (2.0, 4.0, 2.0, 7.0),
        (2.0, 20.0, 0.001, 1.7),
    ] {
        let bank = bank(natural, damping, &[drive], h, false);
        let actual = bank
            .prepare_power_mechanical_gradients(
                ModalToneBudget::realtime_reference(),
                ModalTonePowerGradientBudget::realtime_reference(),
                || false,
            )
            .unwrap()
            .gradients_mode(0, state(0.3, -0.2), &[0.7], &[-0.4], || false)
            .unwrap();
        let scalar = PreparedHarmonicForce::try_new_with_damping_rate(natural, damping, drive, h)
            .unwrap()
            .prepare_power_sensitivities()
            .unwrap()
            .gradients(state(0.3, -0.2), 0.7, -0.4)
            .unwrap();
        for (a, b) in entries(actual).into_iter().zip([
            scalar.input_work.stiffness,
            scalar.viscous_dissipation.stiffness,
            scalar.input_work.damping_rate,
            scalar.viscous_dissipation.damping_rate,
        ]) {
            close(a, b, 2e-12);
        }
    }
    let frequencies: [f64; 32] = core::array::from_fn(|i| 1.8 + 0.01 * i as f64);
    let bank = bank(2.0, 0.3, &frequencies, 0.7, false);
    let mut cosine = [0.0; 32];
    let mut sine = [0.0; 32];
    cosine[31] = 0.7;
    sine[31] = -0.4;
    let actual = bank
        .prepare_power_mechanical_gradients(
            ModalToneBudget::realtime_reference(),
            ModalTonePowerGradientBudget::realtime_reference(),
            || false,
        )
        .unwrap()
        .gradients_mode(0, state(0.3, -0.2), &cosine, &sine, || false)
        .unwrap();
    let scalar = PreparedHarmonicForce::try_new_with_damping_rate(2.0, 0.3, frequencies[31], 0.7)
        .unwrap()
        .prepare_power_sensitivities()
        .unwrap()
        .gradients(state(0.3, -0.2), 0.7, -0.4)
        .unwrap();
    for (a, b) in entries(actual).into_iter().zip([
        scalar.input_work.stiffness,
        scalar.viscous_dissipation.stiffness,
        scalar.input_work.damping_rate,
        scalar.viscous_dissipation.damping_rate,
    ]) {
        close(a, b, 2e-12);
    }
}

#[test]
fn unequal_affine_intervals_reverse_the_complete_mechanical_power_objective() {
    let frequencies: [f64; 2] = [1.8, 2.1];
    let (c0, s0, c1, s1) = ([0.7, -0.4], [-0.2, 0.6], [-0.2, 0.4], [0.2, -0.3]);
    let initial = state(0.3, -0.2);
    let mut history = Vec::new();
    let mut x = initial;
    let mut start = 0.0;
    for h in [0.13, 0.27, 0.41, 0.19] {
        let end = start + h;
        let mut knots = [[0.0; 2]; 4];
        for j in 0..2 {
            let (s, c) = (frequencies[j] * start).sin_cos();
            for (offset, t) in [(0, start), (2, end)] {
                let fc = (1.0 - t) * c0[j] + t * c1[j];
                let fs = (1.0 - t) * s0[j] + t * s1[j];
                knots[offset][j] = fc * c + fs * s;
                knots[offset + 1][j] = -fc * s + fs * c;
            }
        }
        let bank = bank(2.0, 0.3, &frequencies, h, true);
        let next = bank
            .advance_mode_affine_with_energy(0, x, &knots[0], &knots[1], &knots[2], &knots[3])
            .unwrap()
            .state;
        history.push((bank, x, knots));
        x = next;
        start = end;
    }
    let (mut qbar, mut vbar) = (0.0, 0.0);
    let (mut kbar, mut cbar) = (0.0, 0.0);
    for (bank, x, knots) in history.iter().rev() {
        let [fc, fs, ec, es] = knots;
        let mechanical = bank
            .prepare_power_mechanical_gradients(
                ModalToneBudget::realtime_reference(),
                ModalTonePowerGradientBudget::realtime_reference(),
                || false,
            )
            .unwrap()
            .gradients_mode_affine(0, *x, fc, fs, ec, es, || false)
            .unwrap();
        let inputs = bank
            .prepare_power_gradients(ModalTonePowerGradientBudget::realtime_reference(), || false)
            .unwrap()
            .gradients_mode_affine(0, *x, fc, fs, ec, es, || false)
            .unwrap();
        let motion = bank
            .prepare_sensitivities(motion_budget(), || false)
            .unwrap()
            .pullback_mode_affine(0, *x, fc, fs, ec, es, qbar, vbar, || false)
            .unwrap();
        kbar += 0.7 * mechanical.stiffness.input_work
            - 0.2 * mechanical.stiffness.viscous_dissipation
            + motion.stiffness;
        cbar += 0.7 * mechanical.damping_rate.input_work
            - 0.2 * mechanical.damping_rate.viscous_dissipation
            + motion.damping_rate;
        qbar = 0.7 * inputs.initial_displacement.input_work
            - 0.2 * inputs.initial_displacement.viscous_dissipation
            + motion.initial_displacement;
        vbar = 0.7 * inputs.initial_velocity.input_work
            - 0.2 * inputs.initial_velocity.viscous_dissipation
            + motion.initial_velocity;
    }
    let full = bank(2.0, 0.3, &frequencies, 1.0, true);
    let expected = full
        .prepare_power_mechanical_gradients(
            ModalToneBudget::realtime_reference(),
            ModalTonePowerGradientBudget::realtime_reference(),
            || false,
        )
        .unwrap()
        .gradients_mode_affine(0, initial, &c0, &s0, &c1, &s1, || false)
        .unwrap();
    close(
        kbar,
        0.7 * expected.stiffness.input_work - 0.2 * expected.stiffness.viscous_dissipation,
        2e-12,
    );
    close(
        cbar,
        0.7 * expected.damping_rate.input_work - 0.2 * expected.damping_rate.viscous_dissipation,
        2e-12,
    );
}

#[test]
fn physical_units_and_bounded_refusals_preserve_original_bank() {
    let cold = ModalToneBudget::realtime_reference();
    let query = ModalTonePowerGradientBudget::realtime_reference();
    let base = bank(2.0, 0.3, &[1.8, 2.1], 0.7, false);
    for bad in [
        ModalToneBudget {
            maximum_tones: 1,
            ..cold
        },
        ModalToneBudget {
            maximum_storage_bytes: 0,
            ..cold
        },
        ModalToneBudget {
            maximum_preparation_work: 0,
            ..cold
        },
        ModalToneBudget {
            maximum_step_work: 0,
            ..cold
        },
    ] {
        assert!(
            base.prepare_power_mechanical_gradients(bad, query, || false)
                .is_err()
        );
    }
    for bad in [
        ModalTonePowerGradientBudget {
            maximum_scratch_bytes: 0,
            ..query
        },
        ModalTonePowerGradientBudget {
            maximum_query_work: 0,
            ..query
        },
    ] {
        assert!(
            base.prepare_power_mechanical_gradients(cold, bad, || false)
                .is_err()
        );
    }
    let mut polls = 0;
    assert!(
        base.prepare_power_mechanical_gradients(cold, query, || {
            polls += 1;
            polls == 20
        })
        .is_err()
    );
    let prepared = base
        .prepare_power_mechanical_gradients(cold, query, || false)
        .unwrap();
    let (c, s) = ([0.7, -0.4], [-0.2, 0.6]);
    let initial = state(0.3, -0.2);
    let original = prepared
        .gradients_mode(0, initial, &c, &s, || false)
        .unwrap();
    let forward = base.advance_mode_with_energy(0, initial, &c, &s).unwrap();
    assert!(
        prepared
            .gradients_mode(1, initial, &c, &s, || false)
            .is_err()
    );
    assert!(
        prepared
            .gradients_mode(0, initial, &c[..1], &s, || false)
            .is_err()
    );
    assert!(
        prepared
            .gradients_mode(0, initial, &[f64::INFINITY, 0.0], &s, || false)
            .is_err()
    );
    assert!(
        prepared
            .gradients_mode(0, state(f64::NAN, 0.0), &c, &s, || false)
            .is_err()
    );
    assert!(
        prepared
            .gradients_mode_affine(0, initial, &c, &s, &c, &s, || false)
            .is_err()
    );
    let mut polls = 0;
    assert!(
        prepared
            .gradients_mode(0, initial, &c, &s, || {
                polls += 1;
                polls == 5
            })
            .is_err()
    );
    assert_eq!(
        original,
        prepared
            .gradients_mode(0, initial, &c, &s, || false)
            .unwrap()
    );
    assert_eq!(
        forward,
        base.advance_mode_with_energy(0, initial, &c, &s).unwrap()
    );
    for affine in [false, true] {
        let source = bank(2.0, 0.3, &[1.8, 2.1], 0.7, affine);
        let prepared = source
            .prepare_power_mechanical_gradients(cold, query, || false)
            .unwrap();
        if affine {
            assert!(
                prepared
                    .gradients_mode(0, initial, &c, &s, || false)
                    .is_err()
            );
        }
        let (ec, es) = ([-0.2, 0.4], [0.2, -0.3]);
        let actual = if affine {
            prepared.gradients_mode_affine(0, initial, &c, &s, &ec, &es, || false)
        } else {
            prepared.gradients_mode(0, initial, &c, &s, || false)
        }
        .unwrap();
        // q'=a*q, t'=b*t gives dJ'/dk'=a²*dJ/dk, dJ'/dc'=a²/b*dJ/dc.
        let (a, b) = (1e3, 1e-2);
        let changed = bank(2.0 / b, 0.3 / b, &[1.8 / b, 2.1 / b], 0.7 * b, affine);
        let scaled = changed
            .prepare_power_mechanical_gradients(cold, query, || false)
            .unwrap();
        let sc = c.map(|f| a / b / b * f);
        let ss = s.map(|f| a / b / b * f);
        let sec = ec.map(|f| a / b / b * f);
        let ses = es.map(|f| a / b / b * f);
        let si = state(a * 0.3, a / b * (-0.2));
        let result = if affine {
            scaled.gradients_mode_affine(0, si, &sc, &ss, &sec, &ses, || false)
        } else {
            scaled.gradients_mode(0, si, &sc, &ss, || false)
        }
        .unwrap();
        let actual = entries(actual);
        let result = entries(result);
        for i in 0..4 {
            close(
                result[i] / a / a * if i < 2 { 1.0 } else { b },
                actual[i],
                2e-12,
            );
        }
    }
}
