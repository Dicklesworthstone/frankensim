//! G1/G3: coherent bank gradients retain every original source and one state.

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
    PreparedModalTones::try_new_with_damping_rates(
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

// Independent physical variational RK4. Integrate displacement relative to q0
// and evaluate cos-1 by -2*sin(theta/2)^2 to retain a near-DC equilibrium's tiny
// residual. Coordinates are [q-q0,v,q_k,v_k,q_c,v_c,q_Omega,v_Omega,q_phi,v_phi].
#[allow(clippy::too_many_arguments)]
fn reference(
    omega: f64,
    c: f64,
    frequencies: &[f64],
    h: f64,
    initial: ModalAcousticState,
    fc: &[f64],
    fs: &[f64],
    source: usize,
) -> [f64; 10] {
    let q0 = initial.displacement_m_sqrt_kg;
    let held = fc.iter().sum::<f64>() - omega * omega * q0;
    let rhs = |t: f64, x: [f64; 10]| {
        let mut force = held;
        for j in 0..fc.len() {
            let angle = frequencies[j] * t;
            force += -2.0 * fc[j] * (0.5 * angle).sin().powi(2) + fs[j] * angle.sin();
        }
        let (sine, cosine) = (frequencies[source] * t).sin_cos();
        let phase = fs[source] * cosine - fc[source] * sine;
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
fn simultaneous_mechanics_and_each_source_match_independent_variational_odes() {
    let frequencies = [0.0, 2.0, 2.03, 5.0];
    let (fc, fs) = ([0.3, 0.7, -0.2, 0.1], [-0.4, 0.2, 0.3, -0.15]);
    let initial = state(0.3, -0.2);
    for (omega, c, h) in [
        (0.0, 1.3, 0.7),
        (2.0, 0.0, 1.1),
        (2.0, 4.0, 0.8),
        (2.0, 20.0, 0.8),
        (2.0, 1e-3, 4.0),
    ] {
        let bank = bank(&[omega], &[c], &frequencies, h);
        let prepared = bank.prepare_sensitivities(budget(), || false).unwrap();
        let frame = prepared
            .advance_mode(0, initial, &fc, &fs, || false)
            .unwrap();
        assert_eq!(
            frame.state,
            bank.advance_mode(0, initial, &fc, &fs).unwrap()
        );
        let pb = prepared
            .pullback_mode(0, initial, &fc, &fs, 0.2, -0.8, || false)
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
            let expected = reference(omega, c, &frequencies, h, initial, &fc, &fs, j);
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
        // State and all force columns form the exact transpose of the original
        // coherent map, even though its independent source clocks are nonlinear.
        let dc = [0.2, -0.3, 0.4, -0.1];
        let ds = [0.7, 0.2, -0.4, 0.3];
        let direction = bank.advance_mode(0, state(0.4, -0.7), &dc, &ds).unwrap();
        let mut contracted = 0.4 * pb.initial_displacement - 0.7 * pb.initial_velocity;
        for j in 0..4 {
            contracted += dc[j] * pb.source_cotangents()[j].cosine_force
                + ds[j] * pb.source_cotangents()[j].sine_force;
        }
        close(
            contracted,
            0.2 * direction.displacement_m_sqrt_kg - 0.8 * direction.velocity_m_sqrt_kg_per_s,
            3e-13,
        );
    }
}

#[test]
fn cancelled_duplicates_keep_individual_detuning_phase_and_force_gradients() {
    let bank = bank(&[2.0], &[0.2], &[2.0, 2.0, 2.001, 0.0], 0.7);
    let prepared = bank.prepare_sensitivities(budget(), || false).unwrap();
    let fc = [1.0, -1.0, 0.0, 0.0];
    let fs = [0.4, -0.4, 0.0, 0.0];
    let frame = prepared
        .advance_mode(0, state(0.0, 0.0), &fc, &fs, || false)
        .unwrap();
    assert_eq!(frame.state, state(0.0, 0.0));
    assert_eq!(frame.parameters.displacement_by_stiffness, 0.0);
    assert_eq!(frame.parameters.velocity_by_damping_rate, 0.0);
    let pb = prepared
        .pullback_mode(0, state(0.0, 0.0), &fc, &fs, 0.2, -0.8, || false)
        .unwrap();
    let rows = pb.source_cotangents();
    assert_eq!(rows.len(), 4);
    assert_eq!(rows[0].cosine_force, rows[1].cosine_force);
    assert_eq!(rows[0].sine_force, rows[1].sine_force);
    assert!(rows[0].cosine_force.abs() > 1e-3 && rows[0].forcing_angular_frequency.abs() > 1e-3);
    assert!(rows[0].initial_phase.abs() > 1e-3);
    close(
        rows[0].forcing_angular_frequency + rows[1].forcing_angular_frequency,
        0.0,
        1e-15,
    );
    close(rows[0].initial_phase + rows[1].initial_phase, 0.0, 1e-15);
    assert_ne!(rows[0].cosine_force, rows[2].cosine_force);
    assert_eq!(rows[2].forcing_angular_frequency, 0.0);
    assert_ne!(rows[2].cosine_force, 0.0);
}

#[test]
fn aggregate_near_dc_equilibrium_retains_tiny_damping_derivatives() {
    let frequencies = [0.0, 1e-8, 2e-8];
    let fc = [0.5, 1.25, 2.25];
    let fs = [0.0; 3];
    let bank = bank(&[2.0], &[4.0], &frequencies, 1.0);
    let prepared = bank.prepare_sensitivities(budget(), || false).unwrap();
    let frame = prepared
        .advance_mode(0, state(1.0, 0.0), &fc, &fs, || false)
        .unwrap();
    let expected = reference(2.0, 4.0, &frequencies, 1.0, state(1.0, 0.0), &fc, &fs, 0);
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
    assert_eq!(
        frame.state,
        bank.advance_mode(0, state(1.0, 0.0), &fc, &fs).unwrap()
    );
}

#[test]
fn reverse_bank_trajectory_preserves_original_sources_across_modes_and_steps() {
    let intervals = [0.03, 0.2, 0.07, 0.19];
    let duration: f64 = intervals.iter().sum();
    let frequencies = [2.0, 2.0, 2.03];
    let omegas = [0.0, 3.0];
    let rates = [0.4, 0.6];
    let banks: Vec<_> = intervals
        .iter()
        .map(|&h| bank(&omegas, &rates, &frequencies, h))
        .collect();
    let prepared: Vec<_> = banks
        .iter()
        .map(|b| b.prepare_sensitivities(budget(), || false).unwrap())
        .collect();
    let initial = [state(0.3, -0.2), state(-0.4, 0.1)];
    let initial_c = [[0.7, -0.7, 0.3], [0.5, 0.2, -0.4]];
    let initial_s = [[-0.4, 0.4, 0.2], [0.3, -0.1, 0.5]];
    let mut x = initial;
    let mut fc = initial_c;
    let mut fs = initial_s;
    let mut history = Vec::new();
    for bank in &banks {
        history.push((x, fc, fs));
        for mode in 0..2 {
            x[mode] = bank
                .advance_mode(mode, x[mode], &fc[mode], &fs[mode])
                .unwrap();
            for j in 0..3 {
                (fc[mode][j], fs[mode][j]) = bank
                    .advance_force_phase(j, fc[mode][j], fs[mode][j])
                    .unwrap();
            }
        }
    }
    let endpoint = [(0.2, -0.8), (-0.4, 0.7)];
    let mut state_bar = endpoint;
    let mut cb = [[0.0; 3]; 2];
    let mut sb = cb;
    let mut frequency_bar = [0.0; 3];
    let mut kb = [0.0; 2];
    let mut damping_bar = kb;
    for step in (0..intervals.len()).rev() {
        let (before, fc, fs) = history[step];
        for mode in 0..2 {
            let body = prepared[step]
                .pullback_mode(
                    mode,
                    before[mode],
                    &fc[mode],
                    &fs[mode],
                    state_bar[mode].0,
                    state_bar[mode].1,
                    || false,
                )
                .unwrap();
            kb[mode] += body.stiffness;
            damping_bar[mode] += body.damping_rate;
            state_bar[mode] = (body.initial_displacement, body.initial_velocity);
            for j in 0..3 {
                let source = body.source_cotangents()[j];
                let clock = banks[step]
                    .pullback_force_phase(j, fc[mode][j], fs[mode][j], cb[mode][j], sb[mode][j])
                    .unwrap();
                frequency_bar[j] +=
                    source.forcing_angular_frequency + clock.forcing_angular_frequency;
                cb[mode][j] = source.cosine_force + clock.cosine_force;
                sb[mode][j] = source.sine_force + clock.sine_force;
            }
        }
    }
    let whole = bank(&omegas, &rates, &frequencies, duration);
    let complete = whole.prepare_sensitivities(budget(), || false).unwrap();
    let mut direct_frequency = [0.0; 3];
    let mut phase_bar = [0.0; 3];
    let mut oracle_phase = [0.0; 3];
    for mode in 0..2 {
        let direct = complete
            .pullback_mode(
                mode,
                initial[mode],
                &initial_c[mode],
                &initial_s[mode],
                endpoint[mode].0,
                endpoint[mode].1,
                || false,
            )
            .unwrap();
        close(kb[mode], direct.stiffness, 3e-13);
        close(damping_bar[mode], direct.damping_rate, 3e-13);
        close(state_bar[mode].0, direct.initial_displacement, 3e-13);
        close(state_bar[mode].1, direct.initial_velocity, 3e-13);
        for j in 0..3 {
            let source = direct.source_cotangents()[j];
            direct_frequency[j] += source.forcing_angular_frequency;
            phase_bar[j] += initial_s[mode][j] * cb[mode][j] - initial_c[mode][j] * sb[mode][j];
            close(cb[mode][j], source.cosine_force, 3e-13);
            close(sb[mode][j], source.sine_force, 3e-13);
            let oracle = reference(
                omegas[mode],
                rates[mode],
                &frequencies,
                duration,
                initial[mode],
                &initial_c[mode],
                &initial_s[mode],
                j,
            );
            close(
                source.forcing_angular_frequency,
                endpoint[mode].0 * oracle[6] + endpoint[mode].1 * oracle[7],
                3e-10,
            );
            oracle_phase[j] += endpoint[mode].0 * oracle[8] + endpoint[mode].1 * oracle[9];
        }
    }
    for j in 0..3 {
        close(frequency_bar[j], direct_frequency[j], 3e-13);
        close(phase_bar[j], oracle_phase[j], 3e-10);
    }
}

#[test]
fn additional_budgets_affine_laws_and_cancelled_queries_refuse_without_mutation() {
    let bank = bank(&[2.0, 3.0], &[0.2, 0.3], &[0.0, 2.0, 2.01], 0.3);
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
            polls.get() == 4
        })
        .is_err()
    );
    let prepared = bank.prepare_sensitivities(budget(), || false).unwrap();
    let fc = [0.7, -0.3, 0.2];
    let fs = [0.0, 0.4, -0.1];
    let expected = prepared
        .advance_mode(0, state(0.3, -0.2), &fc, &fs, || false)
        .unwrap();
    polls.set(0);
    assert!(
        prepared
            .advance_mode(0, state(0.3, -0.2), &fc, &fs, || {
                polls.set(polls.get() + 1);
                polls.get() == 3
            })
            .is_err()
    );
    assert!(
        prepared
            .pullback_mode(0, state(0.3, -0.2), &fc, &fs, 0.2, 0.8, || true)
            .is_err()
    );
    assert!(
        prepared
            .advance_mode(7, state(0.0, 0.0), &fc, &fs, || false)
            .is_err()
    );
    assert!(
        prepared
            .advance_mode(0, state(0.0, 0.0), &fc[..2], &fs, || false)
            .is_err()
    );
    assert!(
        prepared
            .pullback_mode(0, state(0.0, 0.0), &fc, &fs, f64::NAN, 0.8, || false)
            .is_err()
    );
    assert_eq!(
        expected,
        prepared
            .advance_mode(0, state(0.3, -0.2), &fc, &fs, || false)
            .unwrap()
    );
    let affine = PreparedModalTones::try_new_affine_with_damping_rates(
        &[2.0],
        &[0.2],
        &[2.0],
        0.3,
        ModalToneBudget::realtime_reference(),
    )
    .unwrap();
    let affine = affine.prepare_sensitivities(budget(), || false).unwrap();
    assert!(
        affine
            .advance_mode(0, state(0.0, 0.0), &[1.0], &[0.0], || false)
            .is_err()
    );
}
