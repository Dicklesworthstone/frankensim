//! G1/G3: physical coherent power covectors, including envelope and phase controls.

use fs_modal::time::{
    ModalAcousticState, ModalAffineTonePowerGradient, ModalToneBudget, ModalTonePowerDerivative,
    ModalTonePowerGradient, ModalTonePowerGradientBudget, ModalToneSensitivityBudget,
    PreparedHarmonicForce, PreparedModalTones,
};

fn state(q: f64, v: f64) -> ModalAcousticState {
    ModalAcousticState {
        displacement_m_sqrt_kg: q,
        velocity_m_sqrt_kg_per_s: v,
    }
}

fn motion_budget() -> ModalToneSensitivityBudget {
    ModalToneSensitivityBudget {
        maximum_storage_bytes: 8 * 1024 * 1024,
        maximum_preparation_work: 100_000_000,
        maximum_query_work: 100_000,
    }
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        a.is_finite() && (a - b).abs() <= tolerance * (1.0 + b.abs()),
        "actual={a:.17e}, expected={b:.17e}"
    );
}

fn pair_close(a: ModalTonePowerDerivative, b: ModalTonePowerDerivative, tolerance: f64) {
    close(a.input_work, b.input_work, tolerance);
    close(a.viscous_dissipation, b.viscous_dissipation, tolerance);
}

fn constant_entries(g: &ModalTonePowerGradient) -> Vec<ModalTonePowerDerivative> {
    let mut result = vec![g.initial_displacement, g.initial_velocity];
    for source in g.source_gradients() {
        result.extend([source.cosine_force, source.sine_force, source.initial_phase]);
    }
    result
}

fn affine_entries(g: &ModalAffineTonePowerGradient) -> Vec<ModalTonePowerDerivative> {
    let mut result = vec![g.initial_displacement, g.initial_velocity];
    for source in g.source_gradients() {
        result.extend([
            source.cosine_start,
            source.sine_start,
            source.cosine_end,
            source.sine_end,
            source.initial_phase,
        ]);
    }
    result
}

fn bank(natural: f64, c: f64, frequencies: &[f64], h: f64, affine: bool) -> PreparedModalTones {
    if affine {
        PreparedModalTones::try_new_affine_with_damping_rates(
            &[natural],
            &[c],
            frequencies,
            h,
            ModalToneBudget::realtime_reference(),
        )
        .unwrap()
    } else {
        PreparedModalTones::try_new_with_damping_rates(
            &[natural],
            &[c],
            frequencies,
            h,
            ModalToneBudget::realtime_reference(),
        )
        .unwrap()
    }
}

// Independent physical q/v/work/heat ODE with one variational equation for
// every input. No production factor, transition or energy difference is used.
#[allow(clippy::too_many_arguments)]
fn reference(
    natural: f64,
    c: f64,
    frequencies: &[f64],
    h: f64,
    initial: ModalAcousticState,
    c0: &[f64],
    s0: &[f64],
    c1: &[f64],
    s1: &[f64],
    affine: bool,
) -> Vec<f64> {
    let width = if affine { 4 } else { 2 };
    let controls = 2 + (width + 1) * frequencies.len();
    let rhs = |t: f64, x: &[f64]| {
        let u = if affine { t / h } else { 0.0 };
        let mut force = 0.0;
        let mut df = vec![0.0; controls];
        for (j, &omega) in frequencies.iter().enumerate() {
            let (sine, cosine) = (omega * t).sin_cos();
            let ca = (1.0 - u) * c0[j] + u * c1[j];
            let sa = (1.0 - u) * s0[j] + u * s1[j];
            force += ca * cosine + sa * sine;
            let p = 2 + (width + 1) * j;
            df[p] = (1.0 - u) * cosine;
            df[p + 1] = (1.0 - u) * sine;
            if affine {
                df[p + 2] = u * cosine;
                df[p + 3] = u * sine;
            }
            df[p + width] = sa * cosine - ca * sine;
        }
        let mut out = vec![0.0; x.len()];
        out[..4].copy_from_slice(&[
            x[1],
            force - natural * natural * x[0] - c * x[1],
            force * x[1],
            c * x[1] * x[1],
        ]);
        for (p, &df) in df.iter().enumerate() {
            let i = 4 + 4 * p;
            out[i] = x[i + 1];
            out[i + 1] = df - natural * natural * x[i] - c * x[i + 1];
            out[i + 2] = df * x[1] + force * x[i + 1];
            out[i + 3] = 2.0 * c * x[1] * x[i + 1];
        }
        out
    };
    let add = |x: &[f64], y: &[f64], scale: f64| {
        x.iter()
            .zip(y)
            .map(|(x, y)| x + scale * y)
            .collect::<Vec<_>>()
    };
    let mut x = vec![0.0; 4 + 4 * controls];
    x[0] = initial.displacement_m_sqrt_kg;
    x[1] = initial.velocity_m_sqrt_kg_per_s;
    x[4] = 1.0;
    x[9] = 1.0;
    let n = 8192;
    let dt = h / n as f64;
    for j in 0..n {
        let t = dt * j as f64;
        let a = rhs(t, &x);
        let b = rhs(t + 0.5 * dt, &add(&x, &a, 0.5 * dt));
        let d = rhs(t + 0.5 * dt, &add(&x, &b, 0.5 * dt));
        let e = rhs(t + dt, &add(&x, &d, dt));
        for i in 0..x.len() {
            x[i] += (dt / 6.0) * (a[i] + 2.0 * b[i] + 2.0 * d[i] + e[i]);
        }
    }
    x
}

#[test]
fn all_coherent_controls_match_independent_power_variational_dynamics() {
    let (c0, s0, c1, s1) = (
        [0.7, -0.4, 0.2],
        [-0.2, 0.6, 0.5],
        [0.1, 0.2, -0.5],
        [0.4, -0.3, 0.1],
    );
    let initial = state(0.3, -0.2);
    for (natural, c, frequencies, h) in [
        (0.0, 0.0, [0.0, 1.3, 2.1], 0.7),
        (0.0, 3.0, [0.0, 1.3, 2.1], 0.7),
        (2.0, 0.0, [2.0, 2.0001, 3.3], 2.3),
        (2.0, 4.0, [0.0, 2.0, 2.0], 1.7),
        (2.0, 20.0, [0.001, 1.3, 2.1], 1.7),
        (20.0, 0.03, [20.0, 20.0001, 32.0], 0.7),
    ] {
        for affine in [false, true] {
            let bank = bank(natural, c, &frequencies, h, affine);
            let gradients = bank
                .prepare_power_gradients(ModalTonePowerGradientBudget::realtime_reference(), || {
                    false
                })
                .unwrap();
            let (forward, values) = if affine {
                (
                    bank.advance_mode_affine_with_energy(0, initial, &c0, &s0, &c1, &s1)
                        .unwrap(),
                    affine_entries(
                        &gradients
                            .gradients_mode_affine(0, initial, &c0, &s0, &c1, &s1, || false)
                            .unwrap(),
                    ),
                )
            } else {
                (
                    bank.advance_mode_with_energy(0, initial, &c0, &s0).unwrap(),
                    constant_entries(
                        &gradients
                            .gradients_mode(0, initial, &c0, &s0, || false)
                            .unwrap(),
                    ),
                )
            };
            let r = reference(
                natural,
                c,
                &frequencies,
                h,
                initial,
                &c0,
                &s0,
                &c1,
                &s1,
                affine,
            );
            close(forward.input_work_j, r[2], 3e-10);
            close(forward.viscous_dissipation_j, r[3], 3e-10);
            for (p, g) in values.iter().enumerate() {
                close(g.input_work, r[4 + 4 * p + 2], 3e-10);
                close(g.viscous_dissipation, r[4 + 4 * p + 3], 3e-10);
                if c == 0.0 {
                    assert_eq!(g.viscous_dissipation, 0.0);
                }
            }
        }
    }
}

#[test]
fn duplicate_sources_keep_independent_power_controls_at_exact_cancellation() {
    let (c, h, v) = (3.0_f64, 0.7_f64, 2.0_f64);
    let bank = bank(0.0, c, &[0.0, 0.0], h, false);
    let prepared = bank
        .prepare_power_gradients(ModalTonePowerGradientBudget::realtime_reference(), || false)
        .unwrap();
    let forward = bank
        .advance_mode_with_energy(0, state(0.2, v), &[0.7, -0.7], &[0.4, -0.4])
        .unwrap();
    assert_eq!(forward.input_work_j, 0.0);
    let g = prepared
        .gradients_mode(0, state(0.2, v), &[0.7, -0.7], &[0.4, -0.4], || false)
        .unwrap();
    let expected_work = v * (1.0 - (-c * h).exp()) / c;
    let expected_heat =
        2.0 * v * ((1.0 - (-c * h).exp()) / c - (1.0 - (-2.0 * c * h).exp()) / (2.0 * c));
    for (i, source) in g.source_gradients().iter().enumerate() {
        close(source.cosine_force.input_work, expected_work, 3e-15);
        close(
            source.cosine_force.viscous_dissipation,
            expected_heat,
            3e-15,
        );
        assert_eq!(source.sine_force, ModalTonePowerDerivative::default());
        close(
            source.initial_phase.input_work,
            if i == 0 {
                0.4 * expected_work
            } else {
                -0.4 * expected_work
            },
            3e-15,
        );
    }
    // Each authored delta overflows, but coalesced endpoint pairs and their
    // independent physical covectors are representable and remain admitted.
    let affine = self::bank(0.0, c, &[0.0, 0.0], 0.01, true);
    let gradient = affine
        .prepare_power_gradients(ModalTonePowerGradientBudget::realtime_reference(), || false)
        .unwrap()
        .gradients_mode_affine(
            0,
            state(0.2, v),
            &[1.0e308, -1.0e308],
            &[0.0, 0.0],
            &[-1.0e308, 1.0e308],
            &[0.0, 0.0],
            || false,
        )
        .unwrap();
    assert_eq!(
        gradient.source_gradients()[0],
        gradient.source_gradients()[1]
    );
    assert!(gradient.source_gradients()[0].cosine_start.input_work > 0.0);
    assert!(gradient.source_gradients()[0].cosine_end.input_work > 0.0);
}

#[test]
fn affine_constant_reduction_and_complete_energy_transpose_agree() {
    let (natural, c, h) = (2.0, 0.4, 1.7);
    let frequencies = [1.8, 3.1];
    let (c0, s0) = ([0.7, -0.4], [-0.2, 0.6]);
    let initial = state(0.3, -0.2);
    let constant = bank(natural, c, &frequencies, h, false);
    let affine = bank(natural, c, &frequencies, h, true);
    let cg = constant
        .prepare_power_gradients(ModalTonePowerGradientBudget::realtime_reference(), || false)
        .unwrap()
        .gradients_mode(0, initial, &c0, &s0, || false)
        .unwrap();
    let ag = affine
        .prepare_power_gradients(ModalTonePowerGradientBudget::realtime_reference(), || false)
        .unwrap()
        .gradients_mode_affine(0, initial, &c0, &s0, &c0, &s0, || false)
        .unwrap();
    pair_close(cg.initial_displacement, ag.initial_displacement, 1e-13);
    pair_close(cg.initial_velocity, ag.initial_velocity, 1e-13);
    for (a, b) in cg.source_gradients().iter().zip(ag.source_gradients()) {
        pair_close(
            a.cosine_force,
            ModalTonePowerDerivative {
                input_work: b.cosine_start.input_work + b.cosine_end.input_work,
                viscous_dissipation: b.cosine_start.viscous_dissipation
                    + b.cosine_end.viscous_dissipation,
            },
            1e-13,
        );
        pair_close(
            a.sine_force,
            ModalTonePowerDerivative {
                input_work: b.sine_start.input_work + b.sine_end.input_work,
                viscous_dissipation: b.sine_start.viscous_dissipation
                    + b.sine_end.viscous_dissipation,
            },
            1e-13,
        );
        pair_close(a.initial_phase, b.initial_phase, 1e-13);
    }
    let (c1, s1) = ([0.1, 0.2], [0.4, -0.3]);
    let next = affine
        .advance_mode_affine_with_energy(0, initial, &c0, &s0, &c1, &s1)
        .unwrap();
    let g = affine
        .prepare_power_gradients(ModalTonePowerGradientBudget::realtime_reference(), || false)
        .unwrap()
        .gradients_mode_affine(0, initial, &c0, &s0, &c1, &s1, || false)
        .unwrap();
    let motion = affine
        .prepare_sensitivities(motion_budget(), || false)
        .unwrap()
        .pullback_mode_affine(
            0,
            initial,
            &c0,
            &s0,
            &c1,
            &s1,
            natural * natural * next.state.displacement_m_sqrt_kg,
            next.state.velocity_m_sqrt_kg_per_s,
            || false,
        )
        .unwrap();
    close(
        g.initial_displacement.input_work - g.initial_displacement.viscous_dissipation,
        motion.initial_displacement - natural * natural * initial.displacement_m_sqrt_kg,
        1e-13,
    );
    close(
        g.initial_velocity.input_work - g.initial_velocity.viscous_dissipation,
        motion.initial_velocity - initial.velocity_m_sqrt_kg_per_s,
        1e-13,
    );
    for (g, m) in g.source_gradients().iter().zip(motion.source_cotangents()) {
        for (g, m) in [
            g.cosine_start,
            g.sine_start,
            g.cosine_end,
            g.sine_end,
            g.initial_phase,
        ]
        .into_iter()
        .zip([
            m.cosine_start,
            m.sine_start,
            m.cosine_end,
            m.sine_end,
            m.initial_phase,
        ]) {
            close(g.input_work - g.viscous_dissipation, m, 1e-13);
        }
    }
}

#[test]
fn full_tone_capacity_preserves_the_last_source_and_original_scalar_gradients() {
    let frequencies: Vec<_> = (0..32).map(|i| 0.1 * (i + 1) as f64).collect();
    let bank = bank(2.0, 0.3, &frequencies, 0.7, false);
    let (mut c, mut s) = ([0.0; 32], [0.0; 32]);
    c[31] = 0.7;
    s[31] = -0.4;
    let g = bank
        .prepare_power_gradients(ModalTonePowerGradientBudget::realtime_reference(), || false)
        .unwrap()
        .gradients_mode(0, state(0.3, -0.2), &c, &s, || false)
        .unwrap();
    assert_eq!(g.source_gradients().len(), 32);
    let scalar =
        PreparedHarmonicForce::try_new_with_damping_rate(2.0, 0.3, frequencies[31], 0.7).unwrap();
    let reference = scalar
        .prepare_power_sensitivities()
        .unwrap()
        .gradients(state(0.3, -0.2), 0.7, -0.4)
        .unwrap();
    pair_close(
        g.initial_displacement,
        ModalTonePowerDerivative {
            input_work: reference.input_work.initial_displacement,
            viscous_dissipation: reference.viscous_dissipation.initial_displacement,
        },
        2e-13,
    );
    pair_close(
        g.initial_velocity,
        ModalTonePowerDerivative {
            input_work: reference.input_work.initial_velocity,
            viscous_dissipation: reference.viscous_dissipation.initial_velocity,
        },
        2e-13,
    );
    for (g, w, h) in [
        g.source_gradients()[31].cosine_force,
        g.source_gradients()[31].sine_force,
        g.source_gradients()[31].initial_phase,
    ]
    .into_iter()
    .zip([
        reference.input_work.cosine_force,
        reference.input_work.sine_force,
        reference.input_work.initial_phase,
    ])
    .zip([
        reference.viscous_dissipation.cosine_force,
        reference.viscous_dissipation.sine_force,
        reference.viscous_dissipation.initial_phase,
    ])
    .map(|((g, w), h)| (g, w, h))
    {
        pair_close(
            g,
            ModalTonePowerDerivative {
                input_work: w,
                viscous_dissipation: h,
            },
            2e-13,
        );
    }
}

#[test]
fn unequal_step_affine_reverse_sweep_matches_the_complete_multimode_interval() {
    let frequencies = [1.8_f64, 3.1];
    let initial = [state(0.3, -0.2), state(-0.4, 0.1)];
    let (c0, s0, c1, s1) = ([0.7, -0.4], [-0.2, 0.6], [0.1, 0.2], [0.4, -0.3]);
    let intervals = [0.13, 0.27, 0.41, 0.19];
    let total_h: f64 = intervals.iter().sum();
    let prepare = |h| {
        PreparedModalTones::try_new_affine_with_damping_rates(
            &[2.0, 5.0],
            &[0.3, 1.0],
            &frequencies,
            h,
            ModalToneBudget::realtime_reference(),
        )
        .unwrap()
    };
    let mut states = initial;
    let mut elapsed = 0.0;
    let mut history = Vec::new();
    for h in intervals {
        let bank = prepare(h);
        let mut forces = [[0.0; 2]; 4];
        for j in 0..2 {
            let (sine, cosine) = (frequencies[j] * elapsed).sin_cos();
            for (endpoint, u) in [elapsed / total_h, (elapsed + h) / total_h]
                .into_iter()
                .enumerate()
            {
                let c = (1.0 - u) * c0[j] + u * c1[j];
                let s = (1.0 - u) * s0[j] + u * s1[j];
                forces[2 * endpoint][j] = cosine * c + sine * s;
                forces[2 * endpoint + 1][j] = -sine * c + cosine * s;
            }
        }
        let before = states;
        for mode in 0..2 {
            states[mode] = bank
                .advance_mode_affine(
                    mode,
                    states[mode],
                    &forces[0],
                    &forces[1],
                    &forces[2],
                    &forces[3],
                )
                .unwrap();
        }
        history.push((bank, elapsed, h, forces, before));
        elapsed += h;
    }
    let mut state_bars = [[0.0; 2]; 2];
    let mut knot_bars = [[0.0; 2]; 4];
    let weights = [1.0, 0.6];
    for (bank, elapsed, h, forces, states) in history.iter().rev() {
        let power = bank
            .prepare_power_gradients(ModalTonePowerGradientBudget::realtime_reference(), || false)
            .unwrap();
        let motion = bank
            .prepare_sensitivities(motion_budget(), || false)
            .unwrap();
        for mode in 0..2 {
            let objective = |g: ModalTonePowerDerivative| {
                weights[mode] * (0.7 * g.input_work - 0.2 * g.viscous_dissipation)
            };
            let g = power
                .gradients_mode_affine(
                    mode,
                    states[mode],
                    &forces[0],
                    &forces[1],
                    &forces[2],
                    &forces[3],
                    || false,
                )
                .unwrap();
            let m = motion
                .pullback_mode_affine(
                    mode,
                    states[mode],
                    &forces[0],
                    &forces[1],
                    &forces[2],
                    &forces[3],
                    state_bars[mode][0],
                    state_bars[mode][1],
                    || false,
                )
                .unwrap();
            state_bars[mode] = [
                objective(g.initial_displacement) + m.initial_displacement,
                objective(g.initial_velocity) + m.initial_velocity,
            ];
            for (j, (g, m)) in g
                .source_gradients()
                .iter()
                .zip(m.source_cotangents())
                .enumerate()
            {
                let local = [
                    objective(g.cosine_start) + m.cosine_start,
                    objective(g.sine_start) + m.sine_start,
                    objective(g.cosine_end) + m.cosine_end,
                    objective(g.sine_end) + m.sine_end,
                ];
                let (sine, cosine) = (frequencies[j] * elapsed).sin_cos();
                for (endpoint, u) in [elapsed / total_h, (elapsed + h) / total_h]
                    .into_iter()
                    .enumerate()
                {
                    let cb = cosine * local[2 * endpoint] - sine * local[2 * endpoint + 1];
                    let sb = sine * local[2 * endpoint] + cosine * local[2 * endpoint + 1];
                    knot_bars[0][j] += (1.0 - u) * cb;
                    knot_bars[1][j] += (1.0 - u) * sb;
                    knot_bars[2][j] += u * cb;
                    knot_bars[3][j] += u * sb;
                }
            }
        }
    }
    let bank = prepare(total_h);
    let power = bank
        .prepare_power_gradients(ModalTonePowerGradientBudget::realtime_reference(), || false)
        .unwrap();
    let mut expected_knots = [[0.0; 2]; 4];
    let mut expected_phases = [0.0; 2];
    for mode in 0..2 {
        let objective = |g: ModalTonePowerDerivative| {
            weights[mode] * (0.7 * g.input_work - 0.2 * g.viscous_dissipation)
        };
        let g = power
            .gradients_mode_affine(mode, initial[mode], &c0, &s0, &c1, &s1, || false)
            .unwrap();
        close(
            state_bars[mode][0],
            objective(g.initial_displacement),
            4e-13,
        );
        close(state_bars[mode][1], objective(g.initial_velocity), 4e-13);
        for (j, g) in g.source_gradients().iter().enumerate() {
            for (endpoint, v) in [g.cosine_start, g.sine_start, g.cosine_end, g.sine_end]
                .into_iter()
                .enumerate()
            {
                expected_knots[endpoint][j] += objective(v);
            }
            expected_phases[j] += objective(g.initial_phase);
        }
    }
    for (a, b) in knot_bars
        .into_iter()
        .flatten()
        .zip(expected_knots.into_iter().flatten())
    {
        close(a, b, 4e-13);
    }
    for j in 0..2 {
        close(
            s0[j] * knot_bars[0][j] - c0[j] * knot_bars[1][j] + s1[j] * knot_bars[2][j]
                - c1[j] * knot_bars[3][j],
            expected_phases[j],
            4e-13,
        );
    }
}

#[test]
fn physical_units_budgets_and_cancelled_queries_preserve_coherent_preparation() {
    let bank = bank(2.0, 0.3, &[1.8, 2.1], 0.7, true);
    let budget = ModalTonePowerGradientBudget::realtime_reference();
    assert!(
        bank.prepare_power_gradients(
            ModalTonePowerGradientBudget {
                maximum_scratch_bytes: 0,
                ..budget
            },
            || false
        )
        .is_err()
    );
    assert!(
        bank.prepare_power_gradients(
            ModalTonePowerGradientBudget {
                maximum_query_work: 0,
                ..budget
            },
            || false
        )
        .is_err()
    );
    assert!(bank.prepare_power_gradients(budget, || true).is_err());
    let prepared = bank.prepare_power_gradients(budget, || false).unwrap();
    let (c0, s0, c1, s1) = ([0.7, -0.4], [-0.2, 0.6], [0.1, 0.2], [0.4, -0.3]);
    let original = prepared
        .gradients_mode_affine(0, state(0.3, -0.2), &c0, &s0, &c1, &s1, || false)
        .unwrap();
    assert!(
        prepared
            .gradients_mode(0, state(0.3, -0.2), &c0, &s0, || false)
            .is_err()
    );
    assert!(
        prepared
            .gradients_mode_affine(1, state(0.3, -0.2), &c0, &s0, &c1, &s1, || false)
            .is_err()
    );
    assert!(
        prepared
            .gradients_mode_affine(0, state(0.3, -0.2), &c0[..1], &s0, &c1, &s1, || false)
            .is_err()
    );
    assert!(
        prepared
            .gradients_mode_affine(0, state(f64::NAN, -0.2), &c0, &s0, &c1, &s1, || false)
            .is_err()
    );
    assert!(
        prepared
            .gradients_mode_affine(
                0,
                state(0.3, -0.2),
                &c0,
                &s0,
                &[f64::INFINITY, 0.0],
                &s1,
                || false
            )
            .is_err()
    );
    let mut polls = 0;
    assert!(
        prepared
            .gradients_mode_affine(0, state(0.3, -0.2), &c0, &s0, &c1, &s1, || {
                polls += 1;
                polls == 5
            })
            .is_err()
    );
    assert_eq!(
        original,
        prepared
            .gradients_mode_affine(0, state(0.3, -0.2), &c0, &s0, &c1, &s1, || false)
            .unwrap()
    );
    let (a, b) = (1.0e3, 1.0e-2);
    let changed = self::bank(2.0 / b, 0.3 / b, &[1.8 / b, 2.1 / b], 0.7 * b, true);
    let scale_force = |f: [f64; 2]| f.map(|f| a / b / b * f);
    let scaled = changed
        .prepare_power_gradients(budget, || false)
        .unwrap()
        .gradients_mode_affine(
            0,
            state(a * 0.3, a / b * (-0.2)),
            &scale_force(c0),
            &scale_force(s0),
            &scale_force(c1),
            &scale_force(s1),
            || false,
        )
        .unwrap();
    let mut scales = vec![a, a / b];
    for _ in 0..2 {
        scales.extend([a / b / b, a / b / b, a / b / b, a / b / b, 1.0]);
    }
    let energy = a * a / b / b;
    for ((old, new), scale) in affine_entries(&original)
        .into_iter()
        .zip(affine_entries(&scaled))
        .zip(scales)
    {
        pair_close(
            old,
            ModalTonePowerDerivative {
                input_work: new.input_work / energy * scale,
                viscous_dissipation: new.viscous_dissipation / energy * scale,
            },
            2e-13,
        );
    }
}
