//! Affine carrier envelopes retain actual motion, power, heat and observation.
#![allow(unsafe_code)]
use fs_phs::linear_exact::{
    DenseExponentialBudget, LinearSplitBudget, PreparedExactTones, PreparedToneObserver,
};
use fs_phs::{
    PortHamiltonian, PreparedStepError, RelaxationBranch, ResistanceWitness, mass_spring_damper,
};
use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::Cell;
thread_local! { static TRACK:Cell<(bool,usize)>=const { Cell::new((false,0)) }; }
struct Counter;
fn count() {
    let _ = TRACK.try_with(|c| {
        let (on, n) = c.get();
        if on {
            c.set((true, n + 1));
        }
    });
}
// SAFETY: all pointers/layouts are passed unchanged to System. The thread-local
// counter only observes allocations and never changes allocation behavior.
unsafe impl GlobalAlloc for Counter {
    unsafe fn alloc(&self, l: Layout) -> *mut u8 {
        count();
        unsafe { System.alloc(l) }
    }
    unsafe fn dealloc(&self, p: *mut u8, l: Layout) {
        unsafe { System.dealloc(p, l) }
    }
    unsafe fn realloc(&self, p: *mut u8, l: Layout, n: usize) -> *mut u8 {
        count();
        unsafe { System.realloc(p, l, n) }
    }
}
#[global_allocator]
static ALLOCATOR: Counter = Counter;
fn budget() -> DenseExponentialBudget {
    DenseExponentialBudget {
        max_dimension: 64,
        max_dense_bytes: 128 * 1024 * 1024,
        max_multiply_terms: 20_000_000_000,
    }
}
fn split() -> LinearSplitBudget {
    LinearSplitBudget {
        max_modes: 16,
        max_omega_dt: 0.1,
    }
}
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() < 4e-11 * (1.0 + b.abs()),
        "{a:.17e} versus {b:.17e}"
    );
}
fn rk4<const N: usize>(
    mut x: [f64; N],
    duration: f64,
    steps: usize,
    rhs: impl Fn(f64, &[f64; N]) -> [f64; N],
) -> [f64; N] {
    let h = duration / steps as f64;
    for i in 0..steps {
        let t = i as f64 * h;
        let a = rhs(t, &x);
        let b = rhs(t + h / 2.0, &std::array::from_fn(|j| x[j] + h * a[j] / 2.0));
        let c = rhs(t + h / 2.0, &std::array::from_fn(|j| x[j] + h * b[j] / 2.0));
        let d = rhs(t + h, &std::array::from_fn(|j| x[j] + h * c[j]));
        for j in 0..N {
            x[j] += h * (a[j] + 2.0 * b[j] + 2.0 * c[j] + d[j]) / 6.0;
        }
    }
    x
}
fn affine_force(t: f64, h: f64, w: &[f64], c0: &[f64], s0: &[f64], c1: &[f64], s1: &[f64]) -> f64 {
    (0..w.len())
        .map(|i| {
            (c0[i] + t / h * (c1[i] - c0[i])) * (w[i] * t).cos()
                + (s0[i] + t / h * (s1[i] - s0[i])) * (w[i] * t).sin()
        })
        .sum()
}

#[test]
fn resonant_affine_attack_has_independent_secular_motion_and_signed_work() {
    let (m, w, slope) = (1.7_f64, 3.0_f64, 0.7_f64);
    let h = std::f64::consts::TAU / w;
    let mut source = PreparedExactTones::new_owned_affine(
        mass_spring_damper(m, m * w * w, 0.0).unwrap(),
        &[w],
        h,
        split(),
        budget(),
    )
    .unwrap();
    let mut next = [0.0; 2];
    let mut mean = [0.0];
    let record = source
        .step_affine_into(
            &[0.0; 2],
            &[0.0],
            &[0.0],
            &[slope * h],
            &[0.0],
            &mut next,
            &mut mean,
        )
        .unwrap();
    // Convolution of t*cos(w*t) with sin(w*t)/(m*w), evaluated at one period.
    let q = slope * h / (4.0 * m * w * w);
    let p = slope * h * h / 4.0;
    close(next[0], q);
    close(next[1], p);
    close(mean[0], q / h);
    let independent = rk4([0.0; 3], h, 8192, |t, x| {
        let f = slope * t * (w * t).cos();
        [x[1] / m, -m * w * w * x[0] + f, f * x[1] / m]
    });
    close(record.supplied, independent[2]);
    close(record.supplied, (p * p / m + m * w * w * q * q) / 2.0);
    assert_eq!(record.dissipated, 0.0);
    assert!(source.is_affine());
}

fn material() -> PortHamiltonian {
    mass_spring_damper(0.8, 5.0, 0.15)
        .unwrap()
        .with_relaxation_branches(vec![RelaxationBranch {
            projection: vec![1.0, 0.0],
            stiffness: 2.5,
            relaxation_time_s: 0.08,
        }])
        .unwrap()
}
fn channel(n: usize, row: Vec<f64>) -> ResistanceWitness {
    let matrix = (0..n)
        .flat_map(|i| {
            (0..n).map({
                let row = &row;
                move |j| row[i] * row[j]
            })
        })
        .collect();
    ResistanceWitness::try_from_source_rows(n, matrix, row).unwrap()
}
fn channels() -> Vec<ResistanceWitness> {
    vec![
        channel(3, vec![0.0, 0.15_f64.sqrt(), 0.0]),
        channel(3, vec![0.0, 0.0, 1.0 / 0.08_f64.sqrt()]),
    ]
}

#[test]
fn affine_beats_material_channels_and_causal_pressure_follow_one_original_trajectory() {
    let frequencies = [0.0, 1.99, 2.01];
    let c0 = [0.1, 0.7, -0.3];
    let s0 = [0.0, 0.2, -0.1];
    let c1 = [0.4, 0.0, 0.2];
    let s1 = [0.0, -0.3, 0.4];
    let h = 0.37;
    let mut source = PreparedExactTones::new_owned_affine_with_resistance_channels(
        material(),
        &frequencies,
        h,
        split(),
        budget(),
        &channels(),
    )
    .unwrap();
    let mut observer = PreparedToneObserver::new(
        &source,
        1,
        1,
        &[-6.0],
        &[0.0, 1.0 / 0.8, 0.0],
        &[-1.3],
        &[0.0, 0.2 / 0.8, 0.0],
        budget(),
    )
    .unwrap();
    let initial = [0.03, -0.02, 0.01];
    let filter = [0.04];
    let mut next = [0.0; 3];
    let mut mean = [0.0];
    let mut losses = [0.0; 2];
    let record = source
        .step_affine_observed_into(
            &initial,
            &c0,
            &s0,
            &c1,
            &s1,
            &mut next,
            &mut mean,
            &mut losses,
        )
        .unwrap();
    let mut filter_next = [0.0];
    let mut pressure = [0.0];
    observer
        .preview_affine_into(
            &source,
            &initial,
            &filter,
            &c0,
            &s0,
            &c1,
            &s1,
            &mut filter_next,
            &mut pressure,
        )
        .unwrap();
    let expected = rk4(
        [
            initial[0], initial[1], initial[2], filter[0], 0.0, 0.0, 0.0, 0.0,
        ],
        h,
        8192,
        |t, x| {
            let v = x[1] / 0.8;
            let memory = 2.5_f64.sqrt() * x[0] - x[2];
            let force = affine_force(t, h, &frequencies, &c0, &s0, &c1, &s1);
            [
                v,
                -5.0 * x[0] - 0.15 * v - 2.5_f64.sqrt() * memory + force,
                memory / 0.08,
                -6.0 * x[3] + v,
                force * v,
                0.15 * v * v,
                memory * memory / 0.08,
                v,
            ]
        },
    );
    for i in 0..3 {
        close(next[i], expected[i]);
    }
    close(filter_next[0], expected[3]);
    close(pressure[0], -1.3 * expected[3] + 0.2 * expected[1] / 0.8);
    close(record.supplied, expected[4]);
    for i in 0..2 {
        close(losses[i], expected[5 + i]);
        assert!(losses[i] > 0.0);
    }
    close(record.dissipated, losses.iter().sum());
    close(mean[0], expected[7] / h);
}

#[test]
fn dc_free_mass_and_duplicate_endpoint_coalescing_retain_small_ramps() {
    let (m, h) = (2.0_f64, 0.4_f64);
    let mut source = PreparedExactTones::new_owned_affine(
        mass_spring_damper(m, 0.0, 0.0).unwrap(),
        &[0.0, -0.0, 0.0],
        h,
        split(),
        budget(),
    )
    .unwrap();
    let initial = [0.1, -0.2];
    let mut next = [0.0; 2];
    let mut mean = [0.0];
    // Endpoint sums are exactly 1 and 2. Subtracting the individual endpoints
    // first can erase a small affine source beside cancelling large terms.
    let record = source
        .step_affine_into(
            &initial,
            &[1e16, 1.0, -1e16],
            &[0.0; 3],
            &[1.0, 4.0, -3.0],
            &[0.0; 3],
            &mut next,
            &mut mean,
        )
        .unwrap();
    let p = initial[1] + 1.5 * h;
    let q = initial[0] + initial[1] / m * h + h * h / m * (0.5 + 1.0 / 6.0);
    close(next[0], q);
    close(next[1], p);
    close(mean[0], (q - initial[0]) / h);
    close(
        record.supplied,
        (p * p - initial[1] * initial[1]) / (2.0 * m),
    );
    assert_eq!(record.dissipated, 0.0);
}

#[test]
fn zero_envelope_delta_matches_existing_constant_motion_work_heat_and_observer() {
    let system = material();
    let frequencies = [0.0, 1.9, 2.1];
    let h = 0.19;
    let mut constant =
        PreparedExactTones::new(&system, &frequencies, h, split(), budget()).unwrap();
    let mut affine =
        PreparedExactTones::new_affine(&system, &frequencies, h, split(), budget()).unwrap();
    assert!(!constant.is_affine());
    let mut a = PreparedToneObserver::new(
        &constant,
        1,
        1,
        &[-4.0],
        &[0.0, 1.0, 0.0],
        &[1.0],
        &[0.0; 3],
        budget(),
    )
    .unwrap();
    let mut b = PreparedToneObserver::new(
        &affine,
        1,
        1,
        &[-4.0],
        &[0.0, 1.0, 0.0],
        &[1.0],
        &[0.0; 3],
        budget(),
    )
    .unwrap();
    let initial = [0.1, -0.03, 0.02];
    let c = [0.2, -0.1, 0.4];
    let s = [0.0, 0.2, -0.3];
    let mut next_a = [0.0; 3];
    let mut next_b = [0.0; 3];
    let mut out_a = [0.0];
    let mut out_b = [0.0];
    let ra = constant
        .step_into(&initial, &c, &s, &mut next_a, &mut out_a)
        .unwrap();
    let rb = affine
        .step_affine_into(&initial, &c, &s, &c, &s, &mut next_b, &mut out_b)
        .unwrap();
    for i in 0..3 {
        close(next_a[i], next_b[i]);
    }
    close(out_a[0], out_b[0]);
    close(ra.supplied, rb.supplied);
    close(ra.dissipated, rb.dissipated);
    close(ra.delta_h, rb.delta_h);
    let mut za = [0.0];
    let mut zb = [0.0];
    a.preview_into(&constant, &initial, &[0.04], &c, &s, &mut za, &mut out_a)
        .unwrap();
    b.preview_affine_into(
        &affine,
        &initial,
        &[0.04],
        &c,
        &s,
        &c,
        &s,
        &mut zb,
        &mut out_b,
    )
    .unwrap();
    close(za[0], zb[0]);
    close(out_a[0], out_b[0]);
}

#[test]
fn affine_phase_subdivision_preserves_work_heat_filter_history_with_zero_warm_allocations() {
    let frequencies = [0.0, 1.7, 1.701];
    let h = 0.8;
    const N: usize = 1000;
    let c0 = [0.2, 0.6, -0.1];
    let s0 = [0.0, 0.2, -0.3];
    let c1 = [0.0, 0.0, 0.4];
    let s1 = [0.0, 0.0, 0.1];
    let system = mass_spring_damper(0.9, 5.0, 0.1).unwrap();
    let mut long =
        PreparedExactTones::new_affine(&system, &frequencies, h, split(), budget()).unwrap();
    let mut long_observer = PreparedToneObserver::new(
        &long,
        1,
        1,
        &[-7.0],
        &[0.0, 1.0],
        &[1.0],
        &[0.0, 0.2],
        budget(),
    )
    .unwrap();
    let initial = [0.03, 0.07];
    let initial_filter = [0.2];
    let mut expected = [0.0; 2];
    let mut mean = [0.0];
    let record = long
        .step_affine_into(&initial, &c0, &s0, &c1, &s1, &mut expected, &mut mean)
        .unwrap();
    let mut expected_filter = [0.0];
    let mut expected_pressure = [0.0];
    long_observer
        .preview_affine_into(
            &long,
            &initial,
            &initial_filter,
            &c0,
            &s0,
            &c1,
            &s1,
            &mut expected_filter,
            &mut expected_pressure,
        )
        .unwrap();
    let mut short =
        PreparedExactTones::new_affine(&system, &frequencies, h / N as f64, split(), budget())
            .unwrap();
    let mut observer = PreparedToneObserver::new(
        &short,
        1,
        1,
        &[-7.0],
        &[0.0, 1.0],
        &[1.0],
        &[0.0, 0.2],
        budget(),
    )
    .unwrap();
    let mut c = c0;
    let mut s = s0;
    let mut dc = std::array::from_fn::<_, 3, _>(|i| (c1[i] - c0[i]) / N as f64);
    let mut ds = std::array::from_fn::<_, 3, _>(|i| (s1[i] - s0[i]) / N as f64);
    let mut x = initial;
    let mut z = initial_filter;
    let mut next = [0.0; 2];
    let mut zn = [0.0];
    let mut pressure = [0.0];
    let mut work = 0.0;
    let mut heat = 0.0;
    TRACK.with(|cell| cell.set((true, 0)));
    for _ in 0..N {
        let ce = std::array::from_fn::<_, 3, _>(|i| c[i] + dc[i]);
        let se = std::array::from_fn::<_, 3, _>(|i| s[i] + ds[i]);
        observer
            .preview_affine_into(&short, &x, &z, &c, &s, &ce, &se, &mut zn, &mut pressure)
            .unwrap();
        let step = short
            .step_affine_into(&x, &c, &s, &ce, &se, &mut next, &mut mean)
            .unwrap();
        x = next;
        z = zn;
        work += step.supplied;
        heat += step.dissipated;
        for i in 0..3 {
            (c[i], s[i]) = short.advance_force_phase(i, ce[i], se[i]).unwrap();
            (dc[i], ds[i]) = short.advance_force_phase(i, dc[i], ds[i]).unwrap();
        }
    }
    let allocations = TRACK.with(|cell| {
        let (_, n) = cell.get();
        cell.set((false, n));
        n
    });
    assert_eq!(allocations, 0);
    for i in 0..2 {
        close(x[i], expected[i]);
    }
    close(z[0], expected_filter[0]);
    close(pressure[0], expected_pressure[0]);
    close(work, record.supplied);
    close(heat, record.dissipated);
}

#[test]
fn smooth_attack_release_has_second_order_affine_reconstruction_ladder() {
    let duration = 1.3;
    let carrier = 2.03;
    let envelope = |t: f64| (std::f64::consts::PI * t / duration).sin().powi(2);
    let expected = rk4([0.0; 2], duration, 16384, |t, x| {
        [
            x[1],
            -4.0 * x[0] - 0.08 * x[1] + envelope(t) * (carrier * t).cos(),
        ]
    });
    let mut errors = Vec::new();
    for steps in [8, 16, 32] {
        let h = duration / steps as f64;
        let mut source = PreparedExactTones::new_owned_affine(
            mass_spring_damper(1.0, 4.0, 0.08).unwrap(),
            &[carrier],
            h,
            split(),
            budget(),
        )
        .unwrap();
        let mut x = [0.0; 2];
        let mut next = [0.0; 2];
        let mut mean = [0.0];
        for i in 0..steps {
            let t = i as f64 * h;
            let (s, c) = (carrier * t).sin_cos();
            source
                .step_affine_into(
                    &x,
                    &[envelope(t) * c],
                    &[-envelope(t) * s],
                    &[envelope(t + h) * c],
                    &[-envelope(t + h) * s],
                    &mut next,
                    &mut mean,
                )
                .unwrap();
            x = next;
        }
        errors.push((x[0] - expected[0]).hypot(x[1] - expected[1]));
    }
    for pair in errors.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!(
            ratio > 3.7 && ratio < 4.3,
            "affine refinement errors {errors:?}, ratio {ratio}"
        );
    }
}

#[test]
fn affine_budget_wrong_family_overflow_and_cancel_refusals_publish_nothing() {
    let system = mass_spring_damper(1.0, 4.0, 0.1).unwrap();
    let small = DenseExponentialBudget {
        max_dimension: 7,
        ..budget()
    };
    assert!(PreparedExactTones::new(&system, &[2.0], 0.1, split(), small).is_ok());
    assert!(PreparedExactTones::new_affine(&system, &[2.0], 0.1, split(), small).is_err());
    assert!(matches!(
        PreparedExactTones::new_owned_affine_with_resistance_channels_controlled(
            mass_spring_damper(1.0, 4.0, 0.1).unwrap(),
            &[2.0],
            0.1,
            split(),
            budget(),
            &[],
            || true
        ),
        Err(PreparedStepError::Cancelled)
    ));
    let mut source =
        PreparedExactTones::new_affine(&system, &[2.0], 0.1, split(), budget()).unwrap();
    let mut constant = PreparedExactTones::new(&system, &[2.0], 0.1, split(), budget()).unwrap();
    let mut next = [7.0; 2];
    let mut mean = [8.0];
    assert!(
        source
            .step_into(&[0.0; 2], &[0.0], &[0.0], &mut next, &mut mean)
            .is_err()
    );
    assert!(
        constant
            .step_affine_into(
                &[0.0; 2],
                &[0.0],
                &[0.0],
                &[0.0],
                &[0.0],
                &mut next,
                &mut mean
            )
            .is_err()
    );
    for (start, end) in [(0.0, f64::NAN), (-f64::MAX, f64::MAX), (0.0, f64::MAX)] {
        assert!(
            source
                .step_affine_into(
                    &[0.0; 2],
                    &[start],
                    &[0.0],
                    &[end],
                    &[0.0],
                    &mut next,
                    &mut mean
                )
                .is_err()
        );
        assert_eq!((next, mean), ([7.0; 2], [8.0]));
    }
    let mut observer = PreparedToneObserver::new(
        &source,
        1,
        1,
        &[-0.2],
        &[0.0, 1.0],
        &[2.0],
        &[0.0; 2],
        budget(),
    )
    .unwrap();
    let mut filter = [9.0];
    let mut pressure = [10.0];
    assert!(
        observer
            .preview_affine_into(
                &source,
                &[0.0; 2],
                &[1e308],
                &[0.0],
                &[0.0],
                &[0.0],
                &[0.0],
                &mut filter,
                &mut pressure
            )
            .is_err()
    );
    assert_eq!((filter, pressure), ([9.0], [10.0]));
    let mut polls = 0;
    assert!(matches!(
        observer.preview_affine_into_controlled(
            &source,
            &[0.0; 2],
            &[0.0],
            &[0.0],
            &[0.0],
            &[1.0],
            &[0.0],
            &mut filter,
            &mut pressure,
            || {
                polls += 1;
                polls == 3
            }
        ),
        Err(PreparedStepError::Cancelled)
    ));
    assert_eq!((filter, pressure), ([9.0], [10.0]));
    let record = source
        .step_affine_into(
            &[0.0; 2],
            &[0.0],
            &[0.0],
            &[1.0],
            &[0.0],
            &mut next,
            &mut mean,
        )
        .unwrap();
    let accepted = (next, mean, record.clone());
    let repeated = source
        .step_affine_into(
            &[0.0; 2],
            &[0.0],
            &[0.0],
            &[1.0],
            &[0.0],
            &mut next,
            &mut mean,
        )
        .unwrap();
    assert_eq!(accepted, (next, mean, repeated));
}
