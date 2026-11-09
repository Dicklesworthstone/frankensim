//! G1/G3: independent power variational ODEs, energy balance and trajectory chains.

use fs_modal::time::{HarmonicForcePullback, ModalAcousticState, PreparedHarmonicForce};

fn state(q: f64, v: f64) -> ModalAcousticState {
    ModalAcousticState {
        displacement_m_sqrt_kg: q,
        velocity_m_sqrt_kg_per_s: v,
    }
}

fn entries(g: HarmonicForcePullback) -> [f64; 8] {
    [
        g.initial_displacement,
        g.initial_velocity,
        g.cosine_force,
        g.sine_force,
        g.stiffness,
        g.damping_rate,
        g.forcing_angular_frequency,
        g.initial_phase,
    ]
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        a.is_finite() && (a - b).abs() <= tolerance * (1.0 + b.abs()),
        "{a:.17e} != {b:.17e}"
    );
}

// Physical [q,v,W,Q] and their eight variational systems. This RK4 oracle does
// not use the production exponential, quadratic forms or power composition.
fn reference(k: f64, c: f64, omega: f64, h: f64, q: f64, v: f64, fc: f64, fs: f64) -> [f64; 36] {
    let rhs = |t: f64, x: [f64; 36]| {
        let (s, z) = (omega * t).sin_cos();
        let force = fc * z + fs * s;
        let phase = fs * z - fc * s;
        let mut result = [0.0; 36];
        result[..4].copy_from_slice(&[
            x[1],
            force - k * x[0] - c * x[1],
            force * x[1],
            c * x[1] * x[1],
        ]);
        for p in 0..8 {
            let df = match p {
                2 => z,
                3 => s,
                6 => t * phase,
                7 => phase,
                _ => 0.0,
            };
            let i = 4 + 4 * p;
            result[i] = x[i + 1];
            result[i + 1] = df
                - k * x[i]
                - c * x[i + 1]
                - if p == 4 { x[0] } else { 0.0 }
                - if p == 5 { x[1] } else { 0.0 };
            result[i + 2] = df * x[1] + force * x[i + 1];
            result[i + 3] = 2.0 * c * x[1] * x[i + 1] + if p == 5 { x[1] * x[1] } else { 0.0 };
        }
        result
    };
    let add = |x: [f64; 36], y: [f64; 36], f: f64| core::array::from_fn(|i| x[i] + f * y[i]);
    let mut x = [0.0; 36];
    x[0] = q;
    x[1] = v;
    x[4] = 1.0;
    x[9] = 1.0;
    let n = 16384;
    let dt = h / n as f64;
    for j in 0..n {
        let t = dt * j as f64;
        let a = rhs(t, x);
        let b = rhs(t + 0.5 * dt, add(x, a, 0.5 * dt));
        let d = rhs(t + 0.5 * dt, add(x, b, 0.5 * dt));
        let e = rhs(t + dt, add(x, d, dt));
        for i in 0..36 {
            x[i] += (dt / 6.0) * (a[i] + 2.0 * b[i] + 2.0 * d[i] + e[i]);
        }
    }
    x
}

#[test]
fn physical_work_and_heat_gradients_match_independent_variational_dynamics() {
    for (natural, c, drive, h) in [
        (0.0, 0.0, 0.0, 0.7),
        (0.0, 3.0, 2.0, 0.7),
        (2.0, 0.0, 2.0, 2.3),
        (2.0, 4.0, 2.0, 1.7),
        (2.0, 20.0, 0.001, 1.7),
        (10.0, 0.1, 10.000001, 2.3),
        (2.0, 0.4, 9.0, 7.0),
    ] {
        let base = PreparedHarmonicForce::try_new_with_damping_rate(natural, c, drive, h).unwrap();
        let before = base
            .advance_with_energy(state(0.3, -0.2), 0.7, -0.4)
            .unwrap();
        let prepared = base.prepare_power_sensitivities().unwrap();
        let g = prepared.gradients(state(0.3, -0.2), 0.7, -0.4).unwrap();
        let r = reference(natural * natural, c, drive, h, 0.3, -0.2, 0.7, -0.4);
        close(before.input_work_j, r[2], 3e-10);
        close(before.viscous_dissipation_j, r[3], 3e-10);
        for (p, (work, heat)) in entries(g.input_work)
            .into_iter()
            .zip(entries(g.viscous_dissipation))
            .enumerate()
        {
            close(work, r[4 + 4 * p + 2], 3e-9);
            close(heat, r[4 + 4 * p + 3], 3e-9);
        }
        assert_eq!(
            before,
            base.advance_with_energy(state(0.3, -0.2), 0.7, -0.4)
                .unwrap()
        );
        assert_eq!(g, prepared.gradients(state(0.3, -0.2), 0.7, -0.4).unwrap());
    }
}

#[test]
fn zero_loss_retains_heat_derivative_and_dc_equilibrium_remains_exact() {
    let (h, v, force) = (1.7_f64, -0.3_f64, 0.6_f64);
    let base = PreparedHarmonicForce::try_new_with_damping_rate(0.0, 0.0, 0.0, h).unwrap();
    let g = base
        .prepare_power_sensitivities()
        .unwrap()
        .gradients(state(0.2, v), force, 0.8)
        .unwrap();
    let heat = entries(g.viscous_dissipation);
    for (i, x) in heat.into_iter().enumerate() {
        if i != 5 {
            assert_eq!(x, 0.0);
        }
    }
    close(
        heat[5],
        v * v * h + v * force * h * h + force * force * h.powi(3) / 3.0,
        3e-15,
    );
    assert!(heat[5] > 0.0);
    // Latent sine amplitude contributes to frequency/phase even at exact DC.
    assert_ne!(g.input_work.forcing_angular_frequency, 0.0);
    assert_ne!(g.input_work.initial_phase, 0.0);
    let equilibrium =
        PreparedHarmonicForce::try_new_with_damping_rate(2.0, 0.7, 0.0, 23.0).unwrap();
    let g = equilibrium
        .prepare_power_sensitivities()
        .unwrap()
        .gradients(state(0.25, 0.0), 1.0, 0.0)
        .unwrap();
    assert!(entries(g.viscous_dissipation).into_iter().all(|v| v == 0.0));
}

#[test]
fn every_input_satisfies_differentiated_independent_energy_balance() {
    let (natural, c, drive, h) = (2.0, 0.35, 2.0, 8.7);
    let initial = state(0.3, -0.2);
    let base = PreparedHarmonicForce::try_new_with_damping_rate(natural, c, drive, h).unwrap();
    let forward = base.advance_with_energy(initial, 0.7, -0.4).unwrap();
    let g = base
        .prepare_power_sensitivities()
        .unwrap()
        .gradients(initial, 0.7, -0.4)
        .unwrap();
    let endpoint = base
        .motion()
        .prepare_sensitivities()
        .unwrap()
        .pullback(
            initial,
            0.7,
            -0.4,
            natural * natural * forward.state.displacement_m_sqrt_kg,
            forward.state.velocity_m_sqrt_kg_per_s,
        )
        .unwrap();
    let mut expected = entries(endpoint);
    expected[0] -= natural * natural * initial.displacement_m_sqrt_kg;
    expected[1] -= initial.velocity_m_sqrt_kg_per_s;
    expected[4] += 0.5
        * (forward.state.displacement_m_sqrt_kg.powi(2) - initial.displacement_m_sqrt_kg.powi(2));
    for ((work, heat), balance) in entries(g.input_work)
        .into_iter()
        .zip(entries(g.viscous_dissipation))
        .zip(expected)
    {
        close(work - heat, balance, 2e-11);
    }
}

#[test]
fn changing_step_trajectory_accumulates_power_state_and_source_clock_chains() {
    let intervals = [0.13, 0.27, 0.41, 0.19];
    let mut history = Vec::new();
    let (mut x, mut fc, mut fs) = (state(0.3, -0.2), 0.7, -0.4);
    for h in intervals {
        let base = PreparedHarmonicForce::try_new_with_damping_rate(2.0, 0.3, 1.8, h).unwrap();
        let next = base.advance_with_energy(x, fc, fs).unwrap();
        let source = base.motion().advance_force_phase(fc, fs).unwrap();
        history.push((base, x, fc, fs));
        x = next.state;
        fc = source.0;
        fs = source.1;
    }
    let mut total = [0.0; 8];
    for (base, x, fc, fs) in history.iter().rev() {
        let local = base
            .prepare_power_sensitivities()
            .unwrap()
            .gradients(*x, *fc, *fs)
            .unwrap();
        let work = entries(local.input_work);
        let heat = entries(local.viscous_dissipation);
        let motion = entries(
            base.motion()
                .prepare_sensitivities()
                .unwrap()
                .pullback(*x, *fc, *fs, total[0], total[1])
                .unwrap(),
        );
        let phase = base
            .motion()
            .pullback_force_phase(*fc, *fs, total[2], total[3])
            .unwrap();
        for p in 0..4 {
            total[p] = 0.7 * work[p] - 0.2 * heat[p] + motion[p];
        }
        total[2] += phase.cosine_force;
        total[3] += phase.sine_force;
        for p in 4..7 {
            total[p] += 0.7 * work[p] - 0.2 * heat[p] + motion[p];
        }
        total[6] += phase.forcing_angular_frequency;
    }
    total[7] = -0.4 * total[2] - 0.7 * total[3];
    let base =
        PreparedHarmonicForce::try_new_with_damping_rate(2.0, 0.3, 1.8, intervals.iter().sum())
            .unwrap();
    let full = base
        .prepare_power_sensitivities()
        .unwrap()
        .gradients(state(0.3, -0.2), 0.7, -0.4)
        .unwrap();
    for ((a, w), heat) in total
        .into_iter()
        .zip(entries(full.input_work))
        .zip(entries(full.viscous_dissipation))
    {
        close(a, 0.7 * w - 0.2 * heat, 3e-13);
    }
}

#[test]
fn physical_units_and_failed_queries_preserve_the_original_energy_cache() {
    let base = PreparedHarmonicForce::try_new_with_damping_rate(2.0, 0.3, 1.8, 0.7).unwrap();
    let prepared = base.prepare_power_sensitivities().unwrap();
    let original = prepared.gradients(state(0.3, -0.2), 0.7, -0.4).unwrap();
    assert!(prepared.gradients(state(f64::NAN, 0.0), 0.7, -0.4).is_err());
    assert!(
        prepared
            .gradients(state(0.3, -0.2), f64::INFINITY, 0.0)
            .is_err()
    );
    assert_eq!(
        original,
        prepared.gradients(state(0.3, -0.2), 0.7, -0.4).unwrap()
    );
    // q'=a*q, t'=b*t => energy'=a^2/b^2 * energy.
    let (a, b) = (1.0e3_f64, 1.0e-2_f64);
    let changed =
        PreparedHarmonicForce::try_new_with_damping_rate(2.0 / b, 0.3 / b, 1.8 / b, 0.7 * b)
            .unwrap();
    let scaled = changed
        .prepare_power_sensitivities()
        .unwrap()
        .gradients(
            state(a * 0.3, a / b * (-0.2)),
            a / b / b * 0.7,
            a / b / b * (-0.4),
        )
        .unwrap();
    let energy = a * a / b / b;
    let inputs = [
        a,
        a / b,
        a / b / b,
        a / b / b,
        1.0 / b / b,
        1.0 / b,
        1.0 / b,
        1.0,
    ];
    for (old, new) in [
        (original.input_work, scaled.input_work),
        (original.viscous_dissipation, scaled.viscous_dissipation),
    ] {
        for ((x, y), scale) in entries(old).into_iter().zip(entries(new)).zip(inputs) {
            close(y / energy * scale, x, 3e-13);
        }
    }
}
