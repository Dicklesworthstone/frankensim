//! G1/G3: coherent work/heat source clocks, independent detuning and phase history.

use fs_modal::time::{
    ModalAcousticState, ModalToneBudget, ModalTonePowerDerivative, ModalTonePowerGradientBudget,
    ModalToneSensitivityBudget, PreparedHarmonicForce, PreparedModalTones,
};

fn state(q: f64, v: f64) -> ModalAcousticState {
    ModalAcousticState {
        displacement_m_sqrt_kg: q,
        velocity_m_sqrt_kg_per_s: v,
    }
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        a.is_finite() && (a - b).abs() <= tolerance * (1.0 + b.abs()),
        "{a:.17e} != {b:.17e}"
    );
}

fn relative(a: f64, b: f64, tolerance: f64) {
    assert!(
        a.is_finite() && b != 0.0 && (a - b).abs() <= tolerance * b.abs(),
        "{a:.17e} != {b:.17e}"
    );
}

fn bank(natural: &[f64], c: &[f64], frequencies: &[f64], h: f64) -> PreparedModalTones {
    PreparedModalTones::try_new_with_damping_rates(
        natural,
        c,
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

// Independent physical ODE and three source-frequency variational systems.
// State is [q,v,W,Q, q_Omega0,v_Omega0,W_Omega0,Q_Omega0, ...].
#[allow(clippy::too_many_arguments)]
fn reference(
    natural: f64,
    c: f64,
    frequencies: [f64; 3],
    h: f64,
    initial: ModalAcousticState,
    cosines: [f64; 3],
    sines: [f64; 3],
) -> [f64; 16] {
    let rhs = |t: f64, x: [f64; 16]| {
        let mut force = 0.0;
        let mut partial = [0.0; 3];
        for j in 0..3 {
            let (s, z) = (frequencies[j] * t).sin_cos();
            force += cosines[j] * z + sines[j] * s;
            partial[j] = t * (sines[j] * z - cosines[j] * s);
        }
        let mut out = [0.0; 16];
        out[..4].copy_from_slice(&[
            x[1],
            force - natural * natural * x[0] - c * x[1],
            force * x[1],
            c * x[1] * x[1],
        ]);
        for j in 0..3 {
            let i = 4 + 4 * j;
            out[i] = x[i + 1];
            out[i + 1] = partial[j] - natural * natural * x[i] - c * x[i + 1];
            out[i + 2] = partial[j] * x[1] + force * x[i + 1];
            out[i + 3] = 2.0 * c * x[1] * x[i + 1];
        }
        out
    };
    let add = |x: [f64; 16], y: [f64; 16], a: f64| core::array::from_fn(|i| x[i] + a * y[i]);
    let mut x = [0.0; 16];
    x[0] = initial.displacement_m_sqrt_kg;
    x[1] = initial.velocity_m_sqrt_kg_per_s;
    let n = 8192;
    let dt = h / n as f64;
    for j in 0..n {
        let t = dt * j as f64;
        let a = rhs(t, x);
        let b = rhs(t + 0.5 * dt, add(x, a, 0.5 * dt));
        let d = rhs(t + 0.5 * dt, add(x, b, 0.5 * dt));
        let e = rhs(t + dt, add(x, d, dt));
        for i in 0..16 {
            x[i] += (dt / 6.0) * (a[i] + 2.0 * b[i] + 2.0 * d[i] + e[i]);
        }
    }
    x
}

#[test]
fn coherent_source_frequency_power_derivatives_match_variational_odes() {
    let (fc, fs) = ([0.7, -0.4, 0.2], [-0.2, 0.6, 0.5]);
    let initial = state(0.3, -0.2);
    for (natural, c, frequencies, h) in [
        (0.0, 0.0, [0.0, 1.3, 2.1], 0.7),
        (0.0, 3.0, [0.0, 1.3, 2.1], 0.7),
        (2.0, 0.0, [2.0, 2.0001, 3.3], 2.3),
        (2.0, 4.0, [0.0, 2.0, 2.0], 1.7),
        (2.0, 20.0, [0.001, 1.3, 2.1], 1.7),
        (20.0, 0.03, [20.0, 20.0001, 32.0], 0.7),
        (2.0, 0.4, [0.0, 1e-10, 1e-10], 0.7),
    ] {
        let bank = bank(&[natural], &[c], &frequencies, h);
        let before = bank.advance_mode_with_energy(0, initial, &fc, &fs).unwrap();
        let prepared = bank
            .prepare_power_frequency_gradients(
                ModalToneBudget::realtime_reference(),
                ModalTonePowerGradientBudget::realtime_reference(),
                || false,
            )
            .unwrap();
        let g = prepared
            .gradients_mode(0, initial, &fc, &fs, || false)
            .unwrap();
        let r = reference(natural, c, frequencies, h, initial, fc, fs);
        for (j, derivative) in g.source_derivatives().iter().enumerate() {
            close(derivative.input_work, r[4 + 4 * j + 2], 4e-10);
            close(derivative.viscous_dissipation, r[4 + 4 * j + 3], 4e-10);
            if c == 0.0 {
                assert_eq!(derivative.viscous_dissipation, 0.0);
            }
        }
        assert_eq!(
            before,
            bank.advance_mode_with_energy(0, initial, &fc, &fs).unwrap()
        );
        assert_eq!(
            g,
            prepared
                .gradients_mode(0, initial, &fc, &fs, || false)
                .unwrap()
        );
    }
}

#[test]
fn cancelling_duplicates_retain_opposite_nonzero_independent_detuning() {
    let (c, h, v) = (3.0_f64, 0.7_f64, 2.0_f64);
    let bank = bank(&[0.0], &[c], &[0.0, 0.0], h);
    let prepared = bank
        .prepare_power_frequency_gradients(
            ModalToneBudget::realtime_reference(),
            ModalTonePowerGradientBudget::realtime_reference(),
            || false,
        )
        .unwrap();
    let (fc, fs) = ([0.7, -0.7], [0.4, -0.4]);
    let original = bank
        .advance_mode_with_energy(0, state(0.2, v), &fc, &fs)
        .unwrap();
    assert_eq!(original.input_work_j, 0.0);
    let g = prepared
        .gradients_mode(0, state(0.2, v), &fc, &fs, || false)
        .unwrap();
    // At DC, df/dOmega_j = S_j*t. Nominal F=0, v=v0*exp(-c*t).
    let expected = fs[0] * v * (1.0 - (1.0 + c * h) * (-c * h).exp()) / (c * c);
    relative(g.source_derivatives()[0].input_work, expected, 4e-14);
    relative(g.source_derivatives()[1].input_work, -expected, 4e-14);
    assert_ne!(g.source_derivatives()[0].viscous_dissipation, 0.0);
    assert_eq!(
        g.source_derivatives()[0].input_work + g.source_derivatives()[1].input_work,
        0.0
    );
    assert_eq!(
        g.source_derivatives()[0].viscous_dissipation
            + g.source_derivatives()[1].viscous_dissipation,
        0.0
    );
    // Even a virtually overflowing h*C direction is harmless when the original
    // cancelled cosine source has exactly zero DC clock partial.
    let long = self::bank(&[0.0], &[c], &[0.0, 0.0], 2.0);
    let long_clock = long
        .prepare_power_frequency_gradients(
            ModalToneBudget::realtime_reference(),
            ModalTonePowerGradientBudget::realtime_reference(),
            || false,
        )
        .unwrap();
    let zero = long_clock
        .gradients_mode(0, state(0.2, v), &[1e308, -1e308], &[0.0, 0.0], || false)
        .unwrap();
    assert!(
        zero.source_derivatives()
            .iter()
            .all(|p| *p == ModalTonePowerDerivative::default())
    );
}

#[test]
fn scalar_power_agrees_including_near_dc_equilibrium_and_last_source() {
    for (natural, c, drive, h, initial, fc, fs) in [
        (0.0, 0.0, 0.0, 0.7, state(0.3, -0.2), 0.7, -0.4),
        (2.0, 0.0, 2.0, 2.3, state(0.3, -0.2), 0.7, -0.4),
        (2.0, 4.0, 2.0, 7.0, state(0.3, -0.2), 0.7, -0.4),
        (2.0, 0.7, 1e-6, 0.7, state(0.25, 0.0), 1.0, 0.0),
    ] {
        let bank = bank(&[natural], &[c], &[drive], h);
        let actual = bank
            .prepare_power_frequency_gradients(
                ModalToneBudget::realtime_reference(),
                ModalTonePowerGradientBudget::realtime_reference(),
                || false,
            )
            .unwrap()
            .gradients_mode(0, initial, &[fc], &[fs], || false)
            .unwrap();
        let scalar =
            PreparedHarmonicForce::try_new_with_damping_rate(natural, c, drive, h).unwrap();
        let expected = scalar
            .prepare_power_sensitivities()
            .unwrap()
            .gradients(initial, fc, fs)
            .unwrap();
        let actual = actual.source_derivatives()[0];
        if drive == 1e-6 {
            relative(
                actual.input_work,
                expected.input_work.forcing_angular_frequency,
                3e-11,
            );
            relative(
                actual.viscous_dissipation,
                expected.viscous_dissipation.forcing_angular_frequency,
                3e-11,
            );
        } else {
            close(
                actual.input_work,
                expected.input_work.forcing_angular_frequency,
                3e-12,
            );
            close(
                actual.viscous_dissipation,
                expected.viscous_dissipation.forcing_angular_frequency,
                3e-12,
            );
        }
    }
    let bank = bank(&[2.0], &[0.3], &[1.8; 32], 0.7);
    let (mut fc, mut fs) = ([0.0; 32], [0.0; 32]);
    fc[31] = 0.7;
    fs[31] = -0.4;
    let actual = bank
        .prepare_power_frequency_gradients(
            ModalToneBudget::realtime_reference(),
            ModalTonePowerGradientBudget::realtime_reference(),
            || false,
        )
        .unwrap()
        .gradients_mode(0, state(0.3, -0.2), &fc, &fs, || false)
        .unwrap();
    assert_eq!(actual.source_derivatives().len(), 32);
    assert!(
        actual.source_derivatives()[..31]
            .iter()
            .all(|p| *p == ModalTonePowerDerivative::default())
    );
    let scalar = PreparedHarmonicForce::try_new_with_damping_rate(2.0, 0.3, 1.8, 0.7).unwrap();
    let expected = scalar
        .prepare_power_sensitivities()
        .unwrap()
        .gradients(state(0.3, -0.2), 0.7, -0.4)
        .unwrap();
    close(
        actual.source_derivatives()[31].input_work,
        expected.input_work.forcing_angular_frequency,
        3e-13,
    );
    close(
        actual.source_derivatives()[31].viscous_dissipation,
        expected.viscous_dissipation.forcing_angular_frequency,
        3e-13,
    );
}

#[test]
fn unequal_step_multimode_reverse_includes_the_complete_source_phase_clock() {
    let frequencies = [1.8, 3.6];
    let initial = [state(0.3, -0.2), state(-0.4, 0.1)];
    let (original_c, original_s) = ([0.7, -0.4], [-0.2, 0.6]);
    let intervals = [0.13, 0.27, 0.41, 0.19];
    let prepare = |h| bank(&[2.0, 5.0], &[0.3, 1.0], &frequencies, h);
    let mut states = initial;
    let (mut fc, mut fs) = (original_c, original_s);
    let mut history = Vec::new();
    for h in intervals {
        let bank = prepare(h);
        let before = states;
        for mode in 0..2 {
            states[mode] = bank.advance_mode(mode, states[mode], &fc, &fs).unwrap();
        }
        let (mut next_c, mut next_s) = ([0.0; 2], [0.0; 2]);
        for source in 0..2 {
            (next_c[source], next_s[source]) = bank
                .advance_force_phase(source, fc[source], fs[source])
                .unwrap();
        }
        history.push((bank, before, fc, fs));
        fc = next_c;
        fs = next_s;
    }
    let mut state_bars = [[0.0; 2]; 2];
    let mut force_bars = [[0.0; 2]; 2];
    let mut frequency_bars = [0.0; 2];
    let weights = [1.0, 0.6];
    for (bank, states, fc, fs) in history.iter().rev() {
        let input = bank
            .prepare_power_gradients(ModalTonePowerGradientBudget::realtime_reference(), || false)
            .unwrap();
        let clock = bank
            .prepare_power_frequency_gradients(
                ModalToneBudget::realtime_reference(),
                ModalTonePowerGradientBudget::realtime_reference(),
                || false,
            )
            .unwrap();
        let motion = bank
            .prepare_sensitivities(motion_budget(), || false)
            .unwrap();
        let mut beginning_forces = [[0.0; 2]; 2];
        // Phase rotates once per shared source, not once per driven mode.
        for source in 0..2 {
            let phase = bank
                .pullback_force_phase(
                    source,
                    fc[source],
                    fs[source],
                    force_bars[source][0],
                    force_bars[source][1],
                )
                .unwrap();
            beginning_forces[source] = [phase.cosine_force, phase.sine_force];
            frequency_bars[source] += phase.forcing_angular_frequency;
        }
        for mode in 0..2 {
            let objective = |p: ModalTonePowerDerivative| {
                weights[mode] * (0.7 * p.input_work - 0.2 * p.viscous_dissipation)
            };
            let power = input
                .gradients_mode(mode, states[mode], fc, fs, || false)
                .unwrap();
            let frequency = clock
                .gradients_mode(mode, states[mode], fc, fs, || false)
                .unwrap();
            let end = motion
                .pullback_mode(
                    mode,
                    states[mode],
                    fc,
                    fs,
                    state_bars[mode][0],
                    state_bars[mode][1],
                    || false,
                )
                .unwrap();
            state_bars[mode] = [
                objective(power.initial_displacement) + end.initial_displacement,
                objective(power.initial_velocity) + end.initial_velocity,
            ];
            for source in 0..2 {
                let p = power.source_gradients()[source];
                let e = end.source_cotangents()[source];
                beginning_forces[source][0] += objective(p.cosine_force) + e.cosine_force;
                beginning_forces[source][1] += objective(p.sine_force) + e.sine_force;
                frequency_bars[source] +=
                    objective(frequency.source_derivatives()[source]) + e.forcing_angular_frequency;
            }
        }
        force_bars = beginning_forces;
    }
    let full = prepare(intervals.iter().sum());
    let clock = full
        .prepare_power_frequency_gradients(
            ModalToneBudget::realtime_reference(),
            ModalTonePowerGradientBudget::realtime_reference(),
            || false,
        )
        .unwrap();
    let mut expected = [0.0; 2];
    for mode in 0..2 {
        let g = clock
            .gradients_mode(mode, initial[mode], &original_c, &original_s, || false)
            .unwrap();
        for source in 0..2 {
            expected[source] += weights[mode]
                * (0.7 * g.source_derivatives()[source].input_work
                    - 0.2 * g.source_derivatives()[source].viscous_dissipation);
        }
    }
    for (a, b) in frequency_bars.into_iter().zip(expected) {
        close(a, b, 6e-13);
    }
    // Frequencies are Omega0 and 2*Omega0; this is an authored parameter chain.
    close(
        frequency_bars[0] + 2.0 * frequency_bars[1],
        expected[0] + 2.0 * expected[1],
        8e-13,
    );
}

#[test]
fn physical_units_and_all_admissions_leave_original_bank_reusable() {
    let bank = bank(&[2.0], &[0.3], &[1.8, 2.1], 0.7);
    let cold = ModalToneBudget::realtime_reference();
    let query = ModalTonePowerGradientBudget::realtime_reference();
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
            bank.prepare_power_frequency_gradients(bad, query, || false)
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
            bank.prepare_power_frequency_gradients(cold, bad, || false)
                .is_err()
        );
    }
    assert!(
        bank.prepare_power_frequency_gradients(cold, query, || true)
            .is_err()
    );
    let mut polls = 0;
    assert!(
        bank.prepare_power_frequency_gradients(cold, query, || {
            polls += 1;
            polls == 4
        })
        .is_err()
    );
    let affine = PreparedModalTones::try_new_affine_with_damping_rates(
        &[2.0],
        &[0.3],
        &[1.8, 2.1],
        0.7,
        cold,
    )
    .unwrap();
    assert!(
        affine
            .prepare_power_frequency_gradients(cold, query, || false)
            .is_err()
    );
    let prepared = bank
        .prepare_power_frequency_gradients(cold, query, || false)
        .unwrap();
    let (fc, fs) = ([0.7, -0.4], [-0.2, 0.6]);
    let initial = state(0.3, -0.2);
    let original = prepared
        .gradients_mode(0, initial, &fc, &fs, || false)
        .unwrap();
    assert!(
        prepared
            .gradients_mode(1, initial, &fc, &fs, || false)
            .is_err()
    );
    assert!(
        prepared
            .gradients_mode(0, initial, &fc[..1], &fs, || false)
            .is_err()
    );
    assert!(
        prepared
            .gradients_mode(0, state(f64::NAN, 0.0), &fc, &fs, || false)
            .is_err()
    );
    assert!(
        prepared
            .gradients_mode(0, initial, &[f64::INFINITY, 0.0], &fs, || false)
            .is_err()
    );
    let mut polls = 0;
    assert!(
        prepared
            .gradients_mode(0, initial, &fc, &fs, || {
                polls += 1;
                polls == 5
            })
            .is_err()
    );
    assert_eq!(
        original,
        prepared
            .gradients_mode(0, initial, &fc, &fs, || false)
            .unwrap()
    );
    // q'=a*q,t'=b*t => J'=a²/b²*J and dJ'/dOmega'=a²/b*dJ/dOmega.
    let (a, b) = (1e3, 1e-2);
    let changed = self::bank(&[2.0 / b], &[0.3 / b], &[1.8 / b, 2.1 / b], 0.7 * b);
    let scaled = changed
        .prepare_power_frequency_gradients(cold, query, || false)
        .unwrap()
        .gradients_mode(
            0,
            state(a * 0.3, a / b * (-0.2)),
            &fc.map(|f| a / b / b * f),
            &fs.map(|f| a / b / b * f),
            || false,
        )
        .unwrap();
    for (old, new) in original
        .source_derivatives()
        .iter()
        .zip(scaled.source_derivatives())
    {
        close(new.input_work * b / a / a, old.input_work, 4e-13);
        close(
            new.viscous_dissipation * b / a / a,
            old.viscous_dissipation,
            4e-13,
        );
    }
}
