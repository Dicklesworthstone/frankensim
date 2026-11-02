//! Continuous sources retain original coupled motion, work and positive heat.
//! The references below use elementary oscillator solutions and a separately
//! written physical ODE, not the prepared augmented exponential or its factors.
#![allow(unsafe_code)]

use fs_phs::linear_exact::{
    DenseExponentialBudget, LinearSplitBudget, PreparedExactTones, PreparedExactZoh,
};
use fs_phs::{
    PortHamiltonian, PreparedStepError, QuadraticStorage, RelaxationBranch, ResistanceWitness,
    duffing_oscillator, mass_spring_damper,
};
use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::Cell;

thread_local! { static TRACK: Cell<(bool, usize)> = const { Cell::new((false, 0)) }; }
struct Counter;
fn count() {
    let _ = TRACK.try_with(|cell| {
        let (enabled, n) = cell.get();
        if enabled {
            cell.set((true, n + 1));
        }
    });
}
// SAFETY: unchanged layouts/pointers are forwarded to System; the thread-local
// counter only observes allocation, including reallocations in a warm step.
unsafe impl GlobalAlloc for Counter {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        count();
        unsafe { System.alloc(layout) }
    }
    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        unsafe { System.dealloc(ptr, layout) }
    }
    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, size: usize) -> *mut u8 {
        count();
        unsafe { System.realloc(ptr, layout, size) }
    }
}
#[global_allocator]
static ALLOCATOR: Counter = Counter;

fn split() -> LinearSplitBudget {
    LinearSplitBudget {
        max_modes: 32,
        max_omega_dt: 0.1,
    }
}
fn dense() -> DenseExponentialBudget {
    DenseExponentialBudget {
        max_dimension: 64,
        max_dense_bytes: 128 * 1024 * 1024,
        max_multiply_terms: 20_000_000_000,
    }
}
fn near(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "{actual:.17e} versus {expected:.17e}; difference {:.4e}",
        actual - expected
    );
}
fn force(t: f64, frequencies: &[f64], cosine: &[f64], sine: &[f64]) -> f64 {
    frequencies
        .iter()
        .zip(cosine)
        .zip(sine)
        .map(|((&w, &c), &s)| c * (w * t).cos() + s * (w * t).sin())
        .sum()
}
fn simpson(f: impl Fn(f64) -> f64, dt: f64) -> f64 {
    const COUNT: usize = 8192;
    let h = dt / COUNT as f64;
    let mut value = f(0.0) + f(dt);
    for i in 1..COUNT {
        value += if i % 2 == 0 { 2.0 } else { 4.0 } * f(i as f64 * h);
    }
    value * h / 3.0
}

fn oscillator(
    t: f64,
    m: f64,
    k: f64,
    damping: f64,
    initial: [f64; 2],
    frequencies: &[f64],
    cosine: &[f64],
    sine: &[f64],
) -> [f64; 2] {
    let mut particular = [0.0; 2];
    let mut start = [0.0; 2];
    for ((&w, &c), &s) in frequencies.iter().zip(cosine).zip(sine) {
        let real = k - m * w * w;
        let imag = damping * w;
        let denom = real * real + imag * imag;
        let a = (real * c - imag * s) / denom;
        let b = (imag * c + real * s) / denom;
        particular[0] += a * (w * t).cos() + b * (w * t).sin();
        particular[1] += w * (-a * (w * t).sin() + b * (w * t).cos());
        start[0] += a;
        start[1] += w * b;
    }
    let a = initial[0] - start[0];
    let v = initial[1] / m - start[1];
    let gamma = damping / (2.0 * m);
    let delta = k / m - gamma * gamma;
    let (q, v) = if delta.abs() < 1e-13 * (k / m) {
        let b = v + gamma * a;
        let e = (-gamma * t).exp();
        ((a + b * t) * e, (b - gamma * (a + b * t)) * e)
    } else if delta > 0.0 {
        let beta = delta.sqrt();
        let b = (v + gamma * a) / beta;
        let (s, c) = (beta * t).sin_cos();
        let e = (-gamma * t).exp();
        (
            e * (a * c + b * s),
            e * ((b * beta - gamma * a) * c + (-a * beta - gamma * b) * s),
        )
    } else {
        let b = (-delta).sqrt();
        let r1 = -gamma + b;
        let r2 = -gamma - b;
        let a1 = (v - r2 * a) / (r1 - r2);
        let a2 = a - a1;
        (
            a1 * (r1 * t).exp() + a2 * (r2 * t).exp(),
            r1 * a1 * (r1 * t).exp() + r2 * a2 * (r2 * t).exp(),
        )
    };
    [q + particular[0], m * (v + particular[1])]
}

#[test]
fn inharmonic_oscillator_matches_closed_motion_and_independent_power_across_damping() {
    for (m, k, zeta, dt) in [
        (1.3_f64, 3.7_f64, 0.0, 0.8),
        (1.3, 3.7, 0.02, 0.8),
        (1.3, 3.7, 1.0, 0.8),
        (1.3, 3.7, 1.7, 0.8),
        (0.7, 0.7 * 320.0 * 320.0, 1e-7, 0.025),
    ] {
        let w = (k / m).sqrt();
        let damping = 2.0 * zeta * (m * k).sqrt();
        let frequencies = [0.0, 0.71 * w, 1.013 * w, 2.31 * w];
        let cosine = [0.7, -0.2, 0.4, 0.1];
        let sine = [12.0, 0.3, -0.17, 0.05]; // A DC sine coefficient must have no effect.
        let initial = [0.03, -0.07];
        let system = mass_spring_damper(m, k, damping).unwrap();
        let mut prepared =
            PreparedExactTones::new(&system, &frequencies, dt, split(), dense()).unwrap();
        let mut next = [0.0; 2];
        let mut out = [0.0];
        let record = prepared
            .step_into(&initial, &cosine, &sine, &mut next, &mut out)
            .unwrap();
        let endpoint = oscillator(dt, m, k, damping, initial, &frequencies, &cosine, &sine);
        for i in 0..2 {
            near(next[i], endpoint[i], 3e-11, 5e-12);
        }
        let velocity =
            |t| oscillator(t, m, k, damping, initial, &frequencies, &cosine, &sine)[1] / m;
        let work = simpson(|t| force(t, &frequencies, &cosine, &sine) * velocity(t), dt);
        let heat = simpson(|t| damping * velocity(t).powi(2), dt);
        near(record.supplied, work, 5e-10, 3e-13);
        near(record.dissipated, heat, 5e-10, 1e-17);
        near(out[0], (endpoint[0] - initial[0]) / dt, 3e-11, 3e-12);
        assert!(record.dissipated >= 0.0);
        near(
            record.delta_h + record.dissipated,
            record.supplied,
            1e-11,
            1e-12,
        );
    }
}

#[test]
fn exact_resonance_has_secular_motion_and_true_positive_work_despite_zero_mean_output() {
    let (m, w, amplitude) = (1.7_f64, 3.0_f64, 0.8_f64);
    let dt = 4.0 * std::f64::consts::PI / w;
    let system = mass_spring_damper(m, m * w * w, 0.0).unwrap();
    let mut prepared = PreparedExactTones::new(&system, &[w], dt, split(), dense()).unwrap();
    let mut next = [0.0; 2];
    let mut out = [0.0];
    let record = prepared
        .step_into(&[0.0; 2], &[amplitude], &[0.0], &mut next, &mut out)
        .unwrap();
    near(next[0], 0.0, 0.0, 2e-13);
    near(next[1], amplitude * dt / 2.0, 2e-13, 2e-13);
    near(
        record.supplied,
        amplitude * amplitude * dt * dt / (8.0 * m),
        5e-13,
        2e-13,
    );
    near(out[0], 0.0, 0.0, 2e-13);
    assert_eq!(record.dissipated, 0.0);
    assert!(record.supplied > 0.5);
    assert!((amplitude * out[0] * dt - record.supplied).abs() > 0.5);
}

fn coupled() -> PortHamiltonian {
    let mut q = vec![0.0; 16];
    q[0] = 6.0;
    q[1] = -2.0;
    q[4] = -2.0;
    q[5] = 9.0;
    q[10] = 1.0 / 1.2;
    q[15] = 1.0 / 0.8;
    let mut j = vec![0.0; 16];
    j[2] = 1.0;
    j[7] = 1.0;
    j[8] = -1.0;
    j[13] = -1.0;
    let mut r = vec![0.0; 16];
    r[10] = 0.08 + 0.05;
    r[15] = 0.12 + 0.05;
    r[11] = -0.05;
    r[14] = -0.05;
    PortHamiltonian::new(
        4,
        2,
        j,
        r,
        vec![0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0],
        Box::new(QuadraticStorage::new(q, 4).unwrap()),
    )
    .unwrap()
    .with_relaxation_branches(vec![RelaxationBranch {
        projection: vec![1.0, -1.0, 0.0, 0.0],
        stiffness: 3.0,
        relaxation_time_s: 0.07,
    }])
    .unwrap()
}
fn witness(n: usize, rows: Vec<f64>) -> ResistanceWitness {
    let mut matrix = vec![0.0; n * n];
    for row in rows.chunks_exact(n) {
        for i in 0..n {
            for j in 0..n {
                matrix[i * n + j] += row[i] * row[j];
            }
        }
    }
    ResistanceWitness::try_from_source_rows(n, matrix, rows).unwrap()
}
fn channels() -> Vec<ResistanceWitness> {
    vec![
        witness(
            5,
            vec![
                0.0,
                0.0,
                0.08_f64.sqrt(),
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.12_f64.sqrt(),
                0.0,
            ],
        ),
        witness(5, vec![0.0, 0.0, 0.05_f64.sqrt(), -0.05_f64.sqrt(), 0.0]),
        witness(5, vec![0.0, 0.0, 0.0, 0.0, 1.0 / 0.07_f64.sqrt()]),
    ]
}

fn physical_rhs(
    t: f64,
    x: &[f64; 11],
    frequencies: &[f64],
    cosine: &[f64],
    sine: &[f64],
) -> [f64; 11] {
    let k = frequencies.len();
    let u1 = force(t, frequencies, &cosine[..k], &sine[..k]);
    let u2 = force(t, frequencies, &cosine[k..], &sine[k..]);
    let v1 = x[2] / 1.2;
    let v2 = x[3] / 0.8;
    let memory = 3.0_f64.sqrt() * (x[0] - x[1]) - x[4];
    [
        v1,
        v2,
        -4.0 * x[0] - 2.0 * (x[0] - x[1]) - 3.0_f64.sqrt() * memory - 0.08 * v1 - 0.05 * (v1 - v2)
            + u1,
        -7.0 * x[1] + 2.0 * (x[0] - x[1]) + 3.0_f64.sqrt() * memory - 0.12 * v2
            + 0.05 * (v1 - v2)
            + u2,
        memory / 0.07,
        u1 * v1 + u2 * v2,
        0.08 * v1 * v1 + 0.12 * v2 * v2,
        0.05 * (v1 - v2).powi(2),
        memory * memory / 0.07,
        v1,
        v2,
    ]
}
fn rk4(
    initial: &[f64; 5],
    dt: f64,
    frequencies: &[f64],
    cosine: &[f64],
    sine: &[f64],
) -> [f64; 11] {
    let mut x = [0.0; 11];
    x[..5].copy_from_slice(initial);
    const COUNT: usize = 4096;
    let h = dt / COUNT as f64;
    for step in 0..COUNT {
        let t = step as f64 * h;
        let a = physical_rhs(t, &x, frequencies, cosine, sine);
        let b = physical_rhs(
            t + h / 2.0,
            &std::array::from_fn(|i| x[i] + h * a[i] / 2.0),
            frequencies,
            cosine,
            sine,
        );
        let c = physical_rhs(
            t + h / 2.0,
            &std::array::from_fn(|i| x[i] + h * b[i] / 2.0),
            frequencies,
            cosine,
            sine,
        );
        let d = physical_rhs(
            t + h,
            &std::array::from_fn(|i| x[i] + h * c[i]),
            frequencies,
            cosine,
            sine,
        );
        for i in 0..11 {
            x[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    x
}

#[test]
fn complete_coupled_memory_and_original_heat_channels_match_independent_physical_ode() {
    let frequencies = [0.0, 2.0, 2.031]; // Close but unequal tones keep their beat cross terms.
    let cosine = [0.2, 0.5, -0.31, -0.1, 0.22, 0.4];
    let sine = [0.0, -0.2, 0.17, 0.0, 0.1, -0.21];
    let initial = [0.1, -0.03, 0.07, -0.11, 0.02];
    let dt = 0.43;
    let mut prepared = PreparedExactTones::new_owned_with_resistance_channels(
        coupled(),
        &frequencies,
        dt,
        split(),
        dense(),
        &channels(),
    )
    .unwrap();
    let mut next = [0.0; 5];
    let mut output = [0.0; 2];
    let mut losses = [0.0; 3];
    let record = prepared
        .step_observed_into(
            &initial,
            &cosine,
            &sine,
            &mut next,
            &mut output,
            &mut losses,
        )
        .unwrap();
    let expected = rk4(&initial, dt, &frequencies, &cosine, &sine);
    for i in 0..5 {
        near(next[i], expected[i], 5e-11, 3e-12);
    }
    near(record.supplied, expected[5], 5e-11, 3e-13);
    for i in 0..3 {
        near(losses[i], expected[i + 6], 5e-11, 3e-14);
        assert!(losses[i] > 0.0);
    }
    near(losses.iter().sum(), record.dissipated, 5e-12, 5e-14);
    for i in 0..2 {
        near(output[i], expected[i + 9] / dt, 5e-11, 3e-12);
    }
    assert_eq!(prepared.state_dim(), 5);
    assert_eq!(prepared.port_dim(), 2);
    assert_eq!(prepared.resistance_channel_count(), 3);
}

#[test]
fn dc_and_exact_duplicate_sources_retain_free_motion_and_small_residual_force() {
    for stiffness in [0.0, 4.0] {
        let system = mass_spring_damper(1.3, stiffness, 0.17).unwrap();
        let mut prepared =
            PreparedExactTones::new(&system, &[0.0, -0.0, 0.0], 0.27, split(), dense()).unwrap();
        let held = PreparedExactZoh::new(&system, 0.27, split()).unwrap();
        let initial = [0.03, -0.02];
        let expected = held.step(&initial, &[1.0]).unwrap();
        let mut next = [0.0; 2];
        let mut out = [0.0];
        let record = prepared
            .step_into(
                &initial,
                &[1e16, 1.0, -1e16],
                &[7.0, -2.0, 0.0],
                &mut next,
                &mut out,
            )
            .unwrap();
        for i in 0..2 {
            near(next[i], expected.x[i], 2e-12, 2e-14);
        }
        near(record.supplied, expected.supplied, 3e-12, 2e-14);
        near(record.dissipated, expected.dissipated, 3e-12, 1e-16);
        near(out[0], expected.y[0], 2e-12, 2e-14);
    }
    let system = mass_spring_damper(1.0, 4.0, 0.1).unwrap();
    let mut prepared =
        PreparedExactTones::new(&system, &[2.0, 2.0], 0.2, split(), dense()).unwrap();
    let mut next = [9.0; 2];
    let mut out = [9.0];
    let record = prepared
        .step_into(
            &[0.0; 2],
            &[1e100, -1e100],
            &[-1e100, 1e100],
            &mut next,
            &mut out,
        )
        .unwrap();
    assert_eq!(next, [0.0; 2]);
    assert_eq!(out, [0.0]);
    assert_eq!(record.supplied, 0.0);
    assert_eq!(record.dissipated, 0.0);
}

#[test]
fn phase_checkpoint_and_subdivision_preserve_cross_tone_work_and_heat() {
    let system = mass_spring_damper(0.9, 5.0, 0.08).unwrap();
    let frequencies = [0.0, 2.31, 2.32, 4.7];
    let mut cosine = [0.2, 0.7, -0.4, 0.17];
    let mut sine = [0.0, -0.2, 0.3, 0.07];
    let initial = [0.03, -0.08];
    let mut long = PreparedExactTones::new(&system, &frequencies, 1.1, split(), dense()).unwrap();
    let mut next = [0.0; 2];
    let mut output = [0.0];
    let reference = long
        .step_into(&initial, &cosine, &sine, &mut next, &mut output)
        .unwrap();
    let endpoint = next;
    let mut short =
        PreparedExactTones::new(&system, &frequencies, 1.1 / 16.0, split(), dense()).unwrap();
    let mut x = initial;
    let mut work = 0.0;
    let mut heat = 0.0;
    for _ in 0..16 {
        let saved = (x, cosine, sine);
        let record = short
            .step_into(&x, &cosine, &sine, &mut next, &mut output)
            .unwrap();
        let first = (next, output, record.clone());
        let repeated = short
            .step_into(&saved.0, &saved.1, &saved.2, &mut next, &mut output)
            .unwrap();
        assert_eq!(first, (next, output, repeated));
        work += record.supplied;
        heat += record.dissipated;
        x = next;
        for tone in 0..frequencies.len() {
            (cosine[tone], sine[tone]) = short
                .advance_force_phase(tone, cosine[tone], sine[tone])
                .unwrap();
        }
    }
    for i in 0..2 {
        near(x[i], endpoint[i], 3e-11, 3e-13);
    }
    near(work, reference.supplied, 3e-11, 3e-13);
    near(heat, reference.dissipated, 3e-11, 3e-14);
}

#[test]
fn static_balance_and_tiny_motion_keep_original_heat_scale() {
    for damping in [1e-200, 1e-8, 0.2, 2.0, 20.0] {
        let system = mass_spring_damper(1.0, 1.0, damping).unwrap();
        let mut prepared =
            PreparedExactTones::new(&system, &[0.0, 2.3], 1.0, split(), dense()).unwrap();
        let mut next = [0.0; 2];
        let mut out = [0.0];
        let steady = prepared
            .step_into(&[1.0, 0.0], &[1.0, 0.0], &[0.0; 2], &mut next, &mut out)
            .unwrap();
        assert_eq!(steady.dissipated, 0.0);
        let q = 1.0 + 1e-10;
        let record = prepared
            .step_into(&[q, 0.0], &[1.0, 0.0], &[0.0; 2], &mut next, &mut out)
            .unwrap();
        // Solve only the true represented physical perturbation about the held
        // equilibrium, so reference heat contains no subtraction of O(1) work.
        let expected = simpson(
            |t| {
                let v = oscillator(t, 1.0, 1.0, damping, [q - 1.0, 0.0], &[], &[], &[])[1];
                damping * v * v
            },
            1.0,
        );
        assert!(record.dissipated > 0.0);
        near(record.dissipated, expected, 4e-10, 0.0);
    }
}

#[test]
fn warm_multitone_channels_and_phase_allocate_nothing() {
    let frequencies = [0.0, 2.0, 2.03];
    let mut cosine = [0.2, 0.5, -0.31, -0.1, 0.22, 0.4];
    let mut sine = [0.0, -0.2, 0.17, 0.0, 0.1, -0.21];
    let mut prepared = PreparedExactTones::new_owned_with_resistance_channels(
        coupled(),
        &frequencies,
        0.003,
        split(),
        dense(),
        &channels(),
    )
    .unwrap();
    let mut x = [0.1, -0.03, 0.07, -0.11, 0.02];
    let mut next = [0.0; 5];
    let mut output = [0.0; 2];
    let mut loss = [0.0; 3];
    prepared
        .step_observed_into(&x, &cosine, &sine, &mut next, &mut output, &mut loss)
        .unwrap();
    TRACK.with(|cell| cell.set((true, 0)));
    for _ in 0..1000 {
        prepared
            .step_observed_into(&x, &cosine, &sine, &mut next, &mut output, &mut loss)
            .unwrap();
        x = next;
        for p in 0..2 {
            for t in 0..3 {
                let i = p * 3 + t;
                (cosine[i], sine[i]) = prepared.advance_force_phase(t, cosine[i], sine[i]).unwrap();
            }
        }
    }
    let allocations = TRACK.with(|cell| {
        let (_, n) = cell.get();
        cell.set((false, n));
        n
    });
    assert_eq!(
        allocations, 0,
        "prepared complete-state/channel stepping and cached phase must allocate nothing"
    );
}

#[test]
fn preparation_and_runtime_refusals_publish_no_partial_outputs() {
    let system = mass_spring_damper(1.0, 4.0, 0.1).unwrap();
    for frequencies in [&[][..], &[-1.0][..], &[f64::NAN][..], &[f64::INFINITY][..]] {
        assert!(PreparedExactTones::new(&system, frequencies, 0.1, split(), dense()).is_err());
    }
    for dt in [0.0, -0.1, f64::NAN] {
        assert!(PreparedExactTones::new(&system, &[2.0], dt, split(), dense()).is_err());
    }
    for budget in [
        DenseExponentialBudget {
            max_dimension: 5,
            ..dense()
        },
        DenseExponentialBudget {
            max_dense_bytes: 1,
            ..dense()
        },
        DenseExponentialBudget {
            max_multiply_terms: 1,
            ..dense()
        },
    ] {
        assert!(PreparedExactTones::new(&system, &[2.0], 0.1, split(), budget).is_err());
    }
    assert!(matches!(
        PreparedExactTones::new_owned_with_resistance_channels_controlled(
            mass_spring_damper(1.0, 4.0, 0.1).unwrap(),
            &[2.0],
            0.1,
            split(),
            dense(),
            &[],
            || true
        ),
        Err(PreparedStepError::Cancelled)
    ));
    let nonlinear = duffing_oscillator(1.0, 4.0, 0.2, 0.1).unwrap();
    assert!(PreparedExactTones::new(&nonlinear, &[2.0], 0.1, split(), dense()).is_err());
    let wrong = witness(1, vec![1.0]);
    assert!(
        PreparedExactTones::new_owned_with_resistance_channels(
            mass_spring_damper(1.0, 4.0, 0.1).unwrap(),
            &[2.0],
            0.1,
            split(),
            dense(),
            &[wrong]
        )
        .is_err()
    );
    let mut prepared = PreparedExactTones::new_owned_with_resistance_channels(
        coupled(),
        &[2.0],
        0.1,
        split(),
        dense(),
        &channels(),
    )
    .unwrap();
    let mut next = [7.0; 5];
    let mut output = [8.0; 2];
    let mut loss = [9.0; 3];
    for cosine in [[f64::NAN, 0.2], [f64::MAX, f64::MAX]] {
        assert!(
            prepared
                .step_observed_into(
                    &[0.1; 5],
                    &cosine,
                    &[0.0; 2],
                    &mut next,
                    &mut output,
                    &mut loss
                )
                .is_err()
        );
        assert_eq!(next, [7.0; 5]);
        assert_eq!(output, [8.0; 2]);
        assert_eq!(loss, [9.0; 3]);
    }
    assert!(
        prepared
            .step_observed_into(
                &[0.1; 5],
                &[0.1; 2],
                &[0.0; 2],
                &mut next,
                &mut output,
                &mut loss[..2]
            )
            .is_err()
    );
    assert_eq!(next, [7.0; 5]);
    assert_eq!(output, [8.0; 2]);
    assert_eq!(loss, [9.0; 3]);
    assert!(prepared.advance_force_phase(1, 0.0, 0.0).is_err());
    assert!(prepared.advance_force_phase(0, f64::NAN, 0.0).is_err());
}
