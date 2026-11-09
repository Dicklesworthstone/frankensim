//! G1/G3: attack/release envelopes, their independent clocks and reverse knots.

use fs_modal::time::{
    ModalAcousticState, ModalToneBudget, ModalToneSensitivityBudget, PreparedModalTones,
};
use std::cell::Cell;

fn state(q: f64, v: f64) -> ModalAcousticState {
    ModalAcousticState {
        displacement_m_sqrt_kg: q,
        velocity_m_sqrt_kg_per_s: v,
    }
}
fn budget() -> ModalToneSensitivityBudget {
    ModalToneSensitivityBudget {
        maximum_storage_bytes: 8 * 1024 * 1024,
        maximum_preparation_work: 100_000_000,
        maximum_query_work: 100_000,
    }
}
fn bank(omega: &[f64], c: &[f64], frequencies: &[f64], h: f64) -> PreparedModalTones {
    PreparedModalTones::try_new_affine_with_damping_rates(
        omega,
        c,
        frequencies,
        h,
        ModalToneBudget::realtime_reference(),
    )
    .unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        a.is_finite() && (a - b).abs() <= tolerance * (1.0 + b.abs()),
        "{a:.17e} vs {b:.17e}"
    );
}

// Independent physical variational ODEs. Keep q-q0 and cos-1 explicitly so a
// tiny near-DC perturbation of static equilibrium does not subtract two O(1)
// accelerations. No prepared transition or sensitivity code enters this oracle.
#[allow(clippy::too_many_arguments)]
fn reference(
    omega: f64,
    c: f64,
    frequencies: &[f64],
    h: f64,
    initial: ModalAcousticState,
    endpoints: [&[f64]; 4],
    source: usize,
) -> [f64; 10] {
    let [c0, s0, c1, s1] = endpoints;
    let q0 = initial.displacement_m_sqrt_kg;
    let held = c0.iter().sum::<f64>() - omega * omega * q0;
    let rhs = |t: f64, x: [f64; 10]| {
        let u = t / h;
        let mut force = held;
        for j in 0..c0.len() {
            let dc = c1[j] - c0[j];
            let fc = c0[j] + u * dc;
            let fs = s0[j] + u * (s1[j] - s0[j]);
            let angle = frequencies[j] * t;
            force += u * dc - 2.0 * fc * (0.5 * angle).sin().powi(2) + fs * angle.sin();
        }
        let (sine, cosine) = (frequencies[source] * t).sin_cos();
        let phase = (s0[source] + u * (s1[source] - s0[source])) * cosine
            - (c0[source] + u * (c1[source] - c0[source])) * sine;
        [
            x[1],
            force - omega * omega * x[0] - c * x[1],
            x[3],
            -(q0 + x[0]) - omega * omega * x[2] - c * x[3],
            x[5],
            -x[1] - omega * omega * x[4] - c * x[5],
            x[7],
            t * phase - omega * omega * x[6] - c * x[7],
            x[9],
            phase - omega * omega * x[8] - c * x[9],
        ]
    };
    let add = |a: [f64; 10], b: [f64; 10], scale: f64| std::array::from_fn(|i| a[i] + scale * b[i]);
    let steps = 16384;
    let dt = h / steps as f64;
    let mut x = [0.0; 10];
    x[1] = initial.velocity_m_sqrt_kg_per_s;
    for j in 0..steps {
        let t = j as f64 * dt;
        let a = rhs(t, x);
        let b = rhs(t + dt * 0.5, add(x, a, dt * 0.5));
        let d = rhs(t + dt * 0.5, add(x, b, dt * 0.5));
        let e = rhs(t + dt, add(x, d, dt));
        x = std::array::from_fn(|i| x[i] + dt / 6.0 * (a[i] + 2.0 * b[i] + 2.0 * d[i] + e[i]));
    }
    x[0] += q0;
    x
}

#[test]
fn affine_mechanics_clocks_and_endpoint_transpose_match_independent_odes() {
    let frequencies = [0.0, 2.0, 2.03, 5.0];
    let (c0, s0) = ([0.3, 0.7, -0.2, 0.1], [-0.4, 0.2, 0.3, -0.15]);
    let (c1, s1) = ([-0.2, 0.4, 0.8, -0.4], [0.2, -0.3, 0.1, 0.6]);
    let initial = state(0.3, -0.2);
    for (omega, c, h) in [
        (0.0, 0.0, 0.7),
        (0.0, 1.3, 0.7),
        (2.0, 0.0, 1.1),
        (2.0, 4.0, 0.8),
        (2.0, 20.0, 0.8),
        (2.0, 1e-3, 4.0),
    ] {
        let bank = bank(&[omega], &[c], &frequencies, h);
        let prepared = bank.prepare_sensitivities(budget(), || false).unwrap();
        let frame = prepared
            .advance_mode_affine(0, initial, &c0, &s0, &c1, &s1, || false)
            .unwrap();
        assert_eq!(
            frame.state,
            bank.advance_mode_affine(0, initial, &c0, &s0, &c1, &s1)
                .unwrap()
        );
        let pb = prepared
            .pullback_mode_affine(0, initial, &c0, &s0, &c1, &s1, 0.2, -0.8, || false)
            .unwrap();
        let p = frame.parameters;
        for j in 0..frequencies.len() {
            let s = frame.source_partials()[j];
            let actual = [
                frame.state.displacement_m_sqrt_kg,
                frame.state.velocity_m_sqrt_kg_per_s,
                p.displacement_by_stiffness,
                p.velocity_by_stiffness,
                p.displacement_by_damping_rate,
                p.velocity_by_damping_rate,
                s.displacement_by_forcing_angular_frequency,
                s.velocity_by_forcing_angular_frequency,
                s.displacement_by_initial_phase,
                s.velocity_by_initial_phase,
            ];
            let expected = reference(omega, c, &frequencies, h, initial, [&c0, &s0, &c1, &s1], j);
            for (a, b) in actual.into_iter().zip(expected) {
                close(a, b, 3e-10);
            }
            close(pb.stiffness, 0.2 * expected[2] - 0.8 * expected[3], 3e-10);
            close(
                pb.damping_rate,
                0.2 * expected[4] - 0.8 * expected[5],
                3e-10,
            );
            close(
                pb.source_cotangents()[j].forcing_angular_frequency,
                0.2 * expected[6] - 0.8 * expected[7],
                3e-10,
            );
        }
        let dc0 = [0.2, -0.3, 0.4, -0.1];
        let ds0 = [0.7, 0.2, -0.4, 0.3];
        let dc1 = [-0.3, 0.1, 0.5, 0.2];
        let ds1 = [-0.1, 0.4, 0.3, -0.7];
        let direction = bank
            .advance_mode_affine(0, state(0.4, -0.7), &dc0, &ds0, &dc1, &ds1)
            .unwrap();
        let mut contracted = 0.4 * pb.initial_displacement - 0.7 * pb.initial_velocity;
        for j in 0..4 {
            let s = pb.source_cotangents()[j];
            contracted += dc0[j] * s.cosine_start
                + ds0[j] * s.sine_start
                + dc1[j] * s.cosine_end
                + ds1[j] * s.sine_end;
        }
        close(
            contracted,
            0.2 * direction.displacement_m_sqrt_kg - 0.8 * direction.velocity_m_sqrt_kg_per_s,
            3e-13,
        );
    }
}

#[test]
fn identical_endpoints_reduce_to_constant_tones_with_all_force_columns() {
    let (frequencies, c0, s0) = ([0.0, 3.0, 3.001], [0.5, 0.7, -0.4], [0.2, -0.3, 0.1]);
    let constant = PreparedModalTones::try_new_with_damping_rates(
        &[2.0],
        &[4.0],
        &frequencies,
        0.7,
        ModalToneBudget::realtime_reference(),
    )
    .unwrap();
    let affine = bank(&[2.0], &[4.0], &frequencies, 0.7);
    let a = affine.prepare_sensitivities(budget(), || false).unwrap();
    let b = constant.prepare_sensitivities(budget(), || false).unwrap();
    let initial = state(0.3, -0.4);
    let fa = a
        .advance_mode_affine(0, initial, &c0, &s0, &c0, &s0, || false)
        .unwrap();
    let fb = b.advance_mode(0, initial, &c0, &s0, || false).unwrap();
    assert_eq!(fa.state, fb.state);
    assert_eq!(fa.parameters, fb.parameters);
    let pa = a
        .pullback_mode_affine(0, initial, &c0, &s0, &c0, &s0, 0.2, -0.8, || false)
        .unwrap();
    let pb = b
        .pullback_mode(0, initial, &c0, &s0, 0.2, -0.8, || false)
        .unwrap();
    for (a, b) in pa.source_cotangents().iter().zip(pb.source_cotangents()) {
        close(a.cosine_start + a.cosine_end, b.cosine_force, 2e-14);
        close(a.sine_start + a.sine_end, b.sine_force, 2e-14);
        close(a.initial_phase, b.initial_phase, 2e-14);
        close(
            a.forcing_angular_frequency,
            b.forcing_angular_frequency,
            2e-14,
        );
    }
}

#[test]
fn cancelled_duplicate_envelopes_keep_independent_clock_columns_without_delta_overflow() {
    let bank = bank(&[2.0], &[0.2], &[2.0, 2.0, 2.001], 0.01);
    let prepared = bank.prepare_sensitivities(budget(), || false).unwrap();
    for amplitude in [1.0, 1e308] {
        let c0 = [amplitude, -amplitude, 0.0];
        let c1 = [-amplitude, amplitude, 0.0];
        let s0 = [0.4 * amplitude, -0.4 * amplitude, 0.0];
        let s1 = [-0.7 * amplitude, 0.7 * amplitude, 0.0];
        let frame = prepared
            .advance_mode_affine(0, state(0.0, 0.0), &c0, &s0, &c1, &s1, || false)
            .unwrap();
        assert_eq!(frame.state, state(0.0, 0.0));
        assert_eq!(frame.parameters.displacement_by_stiffness, 0.0);
        assert_eq!(frame.parameters.velocity_by_damping_rate, 0.0);
        let p = prepared
            .pullback_mode_affine(0, state(0.0, 0.0), &c0, &s0, &c1, &s1, 0.2, -0.8, || false)
            .unwrap();
        let s = p.source_cotangents();
        assert_eq!(s[0].cosine_start, s[1].cosine_start);
        assert_eq!(s[0].sine_end, s[1].sine_end);
        assert!(s[0].forcing_angular_frequency.abs() > 1e-7);
        assert_eq!(
            s[0].forcing_angular_frequency,
            -s[1].forcing_angular_frequency
        );
        assert_eq!(s[0].initial_phase, -s[1].initial_phase);
        assert_eq!(s[2].forcing_angular_frequency, 0.0);
        assert_ne!(s[2].cosine_end, 0.0);
    }
}

#[test]
fn near_dc_mechanics_and_long_decay_complementary_clock_moments_keep_small_terms() {
    let frequencies = [0.0, 1e-8, 2e-8];
    let fc = [0.5, 1.25, 2.25];
    let s0 = [0.0; 3];
    let s1 = [0.0, 1e-8, -2e-8];
    let near = bank(&[2.0], &[4.0], &frequencies, 1.0);
    let p = near.prepare_sensitivities(budget(), || false).unwrap();
    let frame = p
        .advance_mode_affine(0, state(1.0, 0.0), &fc, &s0, &fc, &s1, || false)
        .unwrap();
    let expected = reference(
        2.0,
        4.0,
        &frequencies,
        1.0,
        state(1.0, 0.0),
        [&fc, &s0, &fc, &s1],
        0,
    );
    for (actual, expected) in [
        (frame.parameters.displacement_by_damping_rate, expected[4]),
        (frame.parameters.velocity_by_damping_rate, expected[5]),
    ] {
        assert!(expected.abs() > 1e-20 && expected.abs() < 1e-16);
        assert!(
            (actual / expected - 1.0).abs() < 3e-10,
            "{actual:e} vs {expected:e}"
        );
    }
    // Critical DC kernel g(t)=t*exp(-t). For h >> 1, a start envelope
    // (1-t/h) gives q=2/h, v=-1/h. Its frequency derivative under a unit
    // sine quadrature is driven by t*(1-t/h), giving q=2, v=-1. Subtracting
    // full and end gains at h=1e18 would lose these complementary moments.
    let long = bank(&[1.0], &[2.0], &[0.0], 1e18);
    let p = long.prepare_sensitivities(budget(), || false).unwrap();
    let frame = p
        .advance_mode_affine(0, state(0.0, 0.0), &[0.0], &[1.0], &[0.0], &[0.0], || false)
        .unwrap();
    close(
        frame.source_partials()[0].displacement_by_forcing_angular_frequency,
        2.0,
        2e-12,
    );
    close(
        frame.source_partials()[0].velocity_by_forcing_angular_frequency,
        -1.0,
        2e-12,
    );
    let pull = p
        .pullback_mode_affine(
            0,
            state(0.0, 0.0),
            &[0.0],
            &[1.0],
            &[0.0],
            &[0.0],
            1.0,
            0.0,
            || false,
        )
        .unwrap();
    assert!((pull.source_cotangents()[0].cosine_start / 2e-18 - 1.0).abs() < 2e-12);
}

#[test]
fn reverse_affine_trajectory_preserves_shared_knots_and_absolute_source_clock() {
    let intervals = [0.03, 0.2, 0.07, 0.19];
    let h: f64 = intervals.iter().sum();
    let frequencies = [2.0, 2.0, 2.03];
    let omega = [0.0, 3.0];
    let c = [0.4, 0.6];
    let initial = [state(0.3, -0.2), state(-0.4, 0.1)];
    let endpoints = [
        [
            [0.7, -0.7, 0.3],
            [-0.4, 0.4, 0.2],
            [-0.1, 0.5, -0.3],
            [0.2, -0.1, 0.5],
        ],
        [
            [0.5, 0.2, -0.4],
            [0.3, -0.1, 0.5],
            [-0.3, -0.6, 0.1],
            [0.4, 0.3, -0.2],
        ],
    ];
    let banks: Vec<_> = intervals
        .iter()
        .map(|&dt| bank(&omega, &c, &frequencies, dt))
        .collect();
    let sensitivities: Vec<_> = banks
        .iter()
        .map(|b| b.prepare_sensitivities(budget(), || false).unwrap())
        .collect();
    let whole = bank(&omega, &c, &frequencies, h);
    let direct = whole.prepare_sensitivities(budget(), || false).unwrap();
    let endpoint = [(0.2, -0.8), (-0.4, 0.7)];
    for mode in 0..2 {
        let e = endpoints[mode];
        let mut x = initial[mode];
        let mut t = 0.0;
        let mut history = Vec::new();
        for (step, &dt) in intervals.iter().enumerate() {
            let mut local = [[0.0; 3]; 4];
            for (u, row) in [(t / h, 0), ((t + dt) / h, 2)] {
                for j in 0..3 {
                    let fc = (1.0 - u) * e[0][j] + u * e[2][j];
                    let fs = (1.0 - u) * e[1][j] + u * e[3][j];
                    let (sine, cosine) = (frequencies[j] * t).sin_cos();
                    local[row][j] = cosine * fc + sine * fs;
                    local[row + 1][j] = cosine * fs - sine * fc;
                }
            }
            history.push((x, t, local));
            x = banks[step]
                .advance_mode_affine(mode, x, &local[0], &local[1], &local[2], &local[3])
                .unwrap();
            t += dt;
        }
        let full = whole
            .advance_mode_affine(mode, initial[mode], &e[0], &e[1], &e[2], &e[3])
            .unwrap();
        close(x.displacement_m_sqrt_kg, full.displacement_m_sqrt_kg, 3e-14);
        close(
            x.velocity_m_sqrt_kg_per_s,
            full.velocity_m_sqrt_kg_per_s,
            3e-14,
        );
        let mut sb = endpoint[mode];
        let (mut kb, mut cb) = (0.0, 0.0);
        let mut forces = [[0.0; 3]; 4];
        let mut frequency = [0.0; 3];
        let mut phase = [0.0; 3];
        for step in (0..intervals.len()).rev() {
            let (state, t, local) = history[step];
            let p = sensitivities[step]
                .pullback_mode_affine(
                    mode,
                    state,
                    &local[0],
                    &local[1],
                    &local[2],
                    &local[3],
                    sb.0,
                    sb.1,
                    || false,
                )
                .unwrap();
            sb = (p.initial_displacement, p.initial_velocity);
            kb += p.stiffness;
            cb += p.damping_rate;
            for (j, s) in p.source_cotangents().iter().enumerate() {
                frequency[j] += s.forcing_angular_frequency + t * s.initial_phase;
                phase[j] += s.initial_phase;
                let (sine, cosine) = (frequencies[j] * t).sin_cos();
                for (u, bc, bs) in [
                    (t / h, s.cosine_start, s.sine_start),
                    ((t + intervals[step]) / h, s.cosine_end, s.sine_end),
                ] {
                    let c_bar = cosine * bc - sine * bs;
                    let s_bar = sine * bc + cosine * bs;
                    forces[0][j] += (1.0 - u) * c_bar;
                    forces[1][j] += (1.0 - u) * s_bar;
                    forces[2][j] += u * c_bar;
                    forces[3][j] += u * s_bar;
                }
            }
        }
        let p = direct
            .pullback_mode_affine(
                mode,
                initial[mode],
                &e[0],
                &e[1],
                &e[2],
                &e[3],
                endpoint[mode].0,
                endpoint[mode].1,
                || false,
            )
            .unwrap();
        for (a, b) in [
            (sb.0, p.initial_displacement),
            (sb.1, p.initial_velocity),
            (kb, p.stiffness),
            (cb, p.damping_rate),
        ] {
            close(a, b, 4e-13);
        }
        for (j, s) in p.source_cotangents().iter().enumerate() {
            for (a, b) in [
                (forces[0][j], s.cosine_start),
                (forces[1][j], s.sine_start),
                (forces[2][j], s.cosine_end),
                (forces[3][j], s.sine_end),
                (frequency[j], s.forcing_angular_frequency),
                (phase[j], s.initial_phase),
            ] {
                close(a, b, 4e-13);
            }
            let expected = reference(
                omega[mode],
                c[mode],
                &frequencies,
                h,
                initial[mode],
                [&e[0], &e[1], &e[2], &e[3]],
                j,
            );
            close(
                frequency[j],
                endpoint[mode].0 * expected[6] + endpoint[mode].1 * expected[7],
                3e-10,
            );
        }
    }
}

#[test]
fn affine_budget_law_and_cancellation_refusals_are_reusable() {
    let bank = bank(&[2.0], &[0.3], &[0.0, 2.0], 0.3);
    for b in [
        ModalToneSensitivityBudget {
            maximum_storage_bytes: 0,
            ..budget()
        },
        ModalToneSensitivityBudget {
            maximum_preparation_work: 0,
            ..budget()
        },
        ModalToneSensitivityBudget {
            maximum_query_work: 0,
            ..budget()
        },
    ] {
        assert!(bank.prepare_sensitivities(b, || false).is_err());
    }
    let polls = Cell::new(0);
    assert!(
        bank.prepare_sensitivities(budget(), || {
            polls.set(polls.get() + 1);
            polls.get() == 5
        })
        .is_err()
    );
    let p = bank.prepare_sensitivities(budget(), || false).unwrap();
    let (c0, s0, c1, s1) = ([0.7, -0.3], [0.0, 0.4], [0.2, 0.8], [0.3, -0.1]);
    let initial = state(0.3, -0.2);
    let expected = p
        .advance_mode_affine(0, initial, &c0, &s0, &c1, &s1, || false)
        .unwrap();
    polls.set(0);
    assert!(
        p.advance_mode_affine(0, initial, &c0, &s0, &c1, &s1, || {
            polls.set(polls.get() + 1);
            polls.get() == 4
        })
        .is_err()
    );
    assert!(p.advance_mode(0, initial, &c0, &s0, || false).is_err());
    assert!(
        p.advance_mode_affine(1, initial, &c0, &s0, &c1, &s1, || false)
            .is_err()
    );
    assert!(
        p.advance_mode_affine(0, initial, &c0, &s0, &c1[..1], &s1, || false)
            .is_err()
    );
    assert!(
        p.pullback_mode_affine(0, initial, &c0, &s0, &c1, &s1, f64::NAN, 0.8, || false)
            .is_err()
    );
    assert!(
        p.pullback_mode_affine(0, initial, &c0, &s0, &c1, &s1, 0.2, 0.8, || true)
            .is_err()
    );
    assert_eq!(
        expected,
        p.advance_mode_affine(0, initial, &c0, &s0, &c1, &s1, || false)
            .unwrap()
    );
}
