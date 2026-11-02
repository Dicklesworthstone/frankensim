//! Observe actual within-step motion, including material memory and close tones.
#![allow(unsafe_code)]
use fs_phs::linear_exact::{
    DenseExponentialBudget, LinearSplitBudget, PreparedExactTones, PreparedToneObserver,
};
use fs_phs::{PreparedStepError, RelaxationBranch, mass_spring_damper};
use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::Cell;

thread_local! { static TRACK: Cell<(bool,usize)> = const { Cell::new((false,0)) }; }
struct Counter;
fn count() {
    let _ = TRACK.try_with(|c| {
        let (on, n) = c.get();
        if on {
            c.set((true, n + 1));
        }
    });
}
// SAFETY: pointers and layouts pass unchanged to System. The counter merely
// observes thread-local allocations and has no influence on returned memory.
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
        max_dense_bytes: 64 * 1024 * 1024,
        max_multiply_terms: 10_000_000_000,
    }
}
fn split() -> LinearSplitBudget {
    LinearSplitBudget {
        max_modes: 16,
        max_omega_dt: 0.1,
    }
}
fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() < 2e-11 * (1.0 + b.abs()),
        "{a:.17e} versus {b:.17e}"
    );
}
fn source(dt: f64, frequencies: &[f64]) -> PreparedExactTones<'static> {
    PreparedExactTones::new_owned(
        mass_spring_damper(1.0, 4.0, 0.0).unwrap(),
        frequencies,
        dt,
        split(),
        budget(),
    )
    .unwrap()
}
fn force(t: f64, w: &[f64], fc: &[f64], fs: &[f64]) -> f64 {
    w.iter()
        .zip(fc)
        .zip(fs)
        .map(|((&w, &fc), &fs)| fc * (w * t).cos() + fs * (w * t).sin())
        .sum()
}
fn analytic_velocity(t: f64, w: &[f64], fc: &[f64], fs: &[f64]) -> f64 {
    w.iter()
        .zip(fc)
        .zip(fs)
        .map(|((&w, &fc), &fs)| {
            (fc * (-w * (w * t).sin() + 2.0 * (2.0 * t).sin())
                + fs * w * ((w * t).cos() - (2.0 * t).cos()))
                / (4.0 - w * w)
        })
        .sum()
}
fn quadrature(dt: f64, f: impl Fn(f64) -> f64) -> f64 {
    const N: usize = 8192;
    let h = dt / N as f64;
    let mut value = f(0.0) + f(dt);
    for i in 1..N {
        value += if i % 2 == 0 { 2.0 } else { 4.0 } * f(i as f64 * h);
    }
    value * h / 3.0
}

#[test]
fn causal_filter_matches_independent_continuous_convolution_and_preserves_signed_transfer() {
    let dt = 0.73;
    let frequencies = [0.0, 1.3, 1.3001];
    let cosine = [0.2, 0.7, -0.4];
    let sine = [17.0, -0.2, 0.3];
    let source = source(dt, &frequencies);
    let mut observer = PreparedToneObserver::new(
        &source,
        1,
        2,
        &[-9.0],
        &[0.0, 1.0],
        &[-1.7, 1.7],
        &[0.0, 0.23, 0.0, -0.23],
        budget(),
    )
    .unwrap();
    let initial_filter = 0.04;
    let mut filter = [0.0];
    let mut output = [0.0; 2];
    observer
        .preview_into(
            &source,
            &[0.0; 2],
            &[initial_filter],
            &cosine,
            &sine,
            &mut filter,
            &mut output,
        )
        .unwrap();
    let expected = (-9.0 * dt).exp() * initial_filter
        + quadrature(dt, |t| {
            (-9.0 * (dt - t)).exp() * analytic_velocity(t, &frequencies, &cosine, &sine)
        });
    near(filter[0], expected);
    near(
        output[0],
        -1.7 * expected + 0.23 * analytic_velocity(dt, &frequencies, &cosine, &sine),
    );
    assert_eq!(output[1], -output[0]);
    // A held initial velocity is zero and loses the complete causal source term.
    assert!((filter[0] - (-9.0 * dt).exp() * initial_filter).abs() > 0.01);
    assert_eq!(observer.observer_state_dim(), 1);
    assert_eq!(observer.output_dim(), 2);
}

fn physical_rhs(t: f64, x: &[f64; 5], w: &[f64], fc: &[f64], fs: &[f64]) -> [f64; 5] {
    let velocity = x[1] / 1.4;
    let memory = 3.0_f64.sqrt() * x[0] - x[2];
    [
        velocity,
        -5.0 * x[0] - 0.12 * velocity - 3.0_f64.sqrt() * memory + force(t, w, fc, fs),
        memory / 0.09,
        -3.0 * x[3] + x[4] + 1.4 * x[1] - 0.3 * x[2],
        -0.7 * x[3] - 5.0 * x[4] + 0.2 * x[0] - 0.4 * x[1] + 0.8 * x[2],
    ]
}
fn reference(mut x: [f64; 5], dt: f64, w: &[f64], fc: &[f64], fs: &[f64]) -> [f64; 5] {
    let h = dt / 4096.0;
    for i in 0..4096 {
        let t = i as f64 * h;
        let a = physical_rhs(t, &x, w, fc, fs);
        let b = physical_rhs(
            t + h / 2.0,
            &std::array::from_fn(|j| x[j] + h * a[j] / 2.0),
            w,
            fc,
            fs,
        );
        let c = physical_rhs(
            t + h / 2.0,
            &std::array::from_fn(|j| x[j] + h * b[j] / 2.0),
            w,
            fc,
            fs,
        );
        let d = physical_rhs(t + h, &std::array::from_fn(|j| x[j] + h * c[j]), w, fc, fs);
        for j in 0..5 {
            x[j] += h * (a[j] + 2.0 * b[j] + 2.0 * c[j] + d[j]) / 6.0;
        }
    }
    x
}

#[test]
fn original_coupling_and_material_history_drive_observer_without_physical_reaction() {
    let system = mass_spring_damper(1.4, 5.0, 0.12)
        .unwrap()
        .with_relaxation_branches(vec![RelaxationBranch {
            projection: vec![1.0, 0.0],
            stiffness: 3.0,
            relaxation_time_s: 0.09,
        }])
        .unwrap();
    let frequencies = [0.0, 1.7, 2.21];
    let cosine = [0.2, 0.5, -0.3];
    let sine = [0.0, -0.2, 0.1];
    let dt = 0.37;
    let mut source =
        PreparedExactTones::new_owned(system, &frequencies, dt, split(), budget()).unwrap();
    let mut observer = PreparedToneObserver::new(
        &source,
        2,
        2,
        &[-3.0, 1.0, -0.7, -5.0],
        &[0.0, 1.4, -0.3, 0.2, -0.4, 0.8],
        &[1.0, -0.3, -0.2, 1.1],
        &[0.1, 0.7, 0.2, -0.2, -0.1, 0.3],
        budget(),
    )
    .unwrap();
    let initial = [0.03, -0.07, 0.2];
    let initial_filter = [0.04, -0.06];
    let mut endpoint = [0.0; 3];
    let mut mean = [0.0];
    let before = source
        .step_into(&initial, &cosine, &sine, &mut endpoint, &mut mean)
        .unwrap();
    let accepted = (endpoint, mean);
    let mut filter = [0.0; 2];
    let mut output = [0.0; 2];
    observer
        .preview_into(
            &source,
            &initial,
            &initial_filter,
            &cosine,
            &sine,
            &mut filter,
            &mut output,
        )
        .unwrap();
    let expected = reference(
        [
            initial[0],
            initial[1],
            initial[2],
            initial_filter[0],
            initial_filter[1],
        ],
        dt,
        &frequencies,
        &cosine,
        &sine,
    );
    for i in 0..3 {
        near(endpoint[i], expected[i]);
    }
    for i in 0..2 {
        near(filter[i], expected[3 + i]);
    }
    near(
        output[0],
        expected[3] - 0.3 * expected[4] + 0.1 * expected[0] + 0.7 * expected[1] + 0.2 * expected[2],
    );
    near(
        output[1],
        -0.2 * expected[3] + 1.1 * expected[4] - 0.2 * expected[0] - 0.1 * expected[1]
            + 0.3 * expected[2],
    );
    let after = source
        .step_into(&initial, &cosine, &sine, &mut endpoint, &mut mean)
        .unwrap();
    assert_eq!((endpoint, mean), accepted);
    assert_eq!(
        before, after,
        "an observer may not alter physical motion, storage, work or heat"
    );
}

#[test]
fn direct_feedthrough_and_duplicate_cancellation_use_same_source_and_survive_owner_move() {
    let mut source = source(0.23, &[1.7, 1.7, 1.7]);
    let mut observer =
        PreparedToneObserver::new(&source, 0, 1, &[], &[], &[], &[0.0, 2.0], budget()).unwrap();
    let cosine = [1e16, 1.0, -1e16];
    let sine = [-1e16, 0.0, 1e16];
    let initial = [0.07, -0.2];
    let mut endpoint = [0.0; 2];
    let mut mean = [0.0];
    source
        .step_into(&initial, &cosine, &sine, &mut endpoint, &mut mean)
        .unwrap();
    let moved = source;
    let mut output = [0.0];
    observer
        .preview_into(&moved, &initial, &[], &cosine, &sine, &mut [], &mut output)
        .unwrap();
    near(output[0], 2.0 * endpoint[1]);
    assert_eq!(observer.observer_state_dim(), 0);
}

#[test]
fn source_identity_capacity_and_late_refusals_preserve_both_caller_outputs() {
    let original = source(0.1, &[1.7]);
    let other = source(0.1, &[1.7]);
    let mut observer = PreparedToneObserver::new(
        &original,
        1,
        1,
        &[-0.2],
        &[0.0, 1.0],
        &[2.0],
        &[0.0; 2],
        budget(),
    )
    .unwrap();
    for bad in [
        DenseExponentialBudget {
            max_dimension: 4,
            ..budget()
        },
        DenseExponentialBudget {
            max_dense_bytes: 1,
            ..budget()
        },
        DenseExponentialBudget {
            max_multiply_terms: 1,
            ..budget()
        },
    ] {
        assert!(
            PreparedToneObserver::new(
                &original,
                1,
                1,
                &[-0.2],
                &[0.0, 1.0],
                &[2.0],
                &[0.0; 2],
                bad
            )
            .is_err()
        );
    }
    assert!(
        PreparedToneObserver::new(
            &original,
            1,
            1,
            &[f64::NAN],
            &[0.0, 1.0],
            &[2.0],
            &[0.0; 2],
            budget()
        )
        .is_err()
    );
    assert!(
        PreparedToneObserver::new(
            &original,
            1,
            1,
            &[-0.2],
            &[1.0],
            &[2.0],
            &[0.0; 2],
            budget()
        )
        .is_err()
    );
    assert!(matches!(
        PreparedToneObserver::new_controlled(
            &original,
            1,
            1,
            &[-0.2],
            &[0.0, 1.0],
            &[2.0],
            &[0.0; 2],
            budget(),
            || true
        ),
        Err(PreparedStepError::Cancelled)
    ));
    let mut filter = [7.0];
    let mut output = [8.0];
    assert!(
        observer
            .preview_into(
                &other,
                &[0.0; 2],
                &[0.0],
                &[0.2],
                &[0.0],
                &mut filter,
                &mut output
            )
            .is_err()
    );
    assert_eq!((filter, output), ([7.0], [8.0]));
    // State candidate is finite; the final pressure alone overflows.
    assert!(
        observer
            .preview_into(
                &original,
                &[0.0; 2],
                &[1e308],
                &[0.0],
                &[0.0],
                &mut filter,
                &mut output
            )
            .is_err()
    );
    assert_eq!((filter, output), ([7.0], [8.0]));
    let mut calls = 0;
    assert!(matches!(
        observer.preview_into_controlled(
            &original,
            &[0.0; 2],
            &[0.0],
            &[0.2],
            &[0.0],
            &mut filter,
            &mut output,
            || {
                calls += 1;
                calls == 3
            }
        ),
        Err(PreparedStepError::Cancelled)
    ));
    assert_eq!((filter, output), ([7.0], [8.0]));
    observer
        .preview_into(
            &original,
            &[0.0; 2],
            &[0.0],
            &[0.2],
            &[0.0],
            &mut filter,
            &mut output,
        )
        .unwrap();
    let accepted = (filter, output);
    observer
        .preview_into(
            &original,
            &[0.0; 2],
            &[0.0],
            &[0.2],
            &[0.0],
            &mut filter,
            &mut output,
        )
        .unwrap();
    assert_eq!((filter, output), accepted);
}

#[test]
fn continuous_observer_subdivision_and_warm_phase_steps_allocate_nothing() {
    let frequencies = [0.0, 1.7, 1.71];
    let mut cosine = [0.2, 0.5, -0.3];
    let mut sine = [0.0, -0.2, 0.1];
    let long = source(1.0, &frequencies);
    let mut reference = PreparedToneObserver::new(
        &long,
        1,
        1,
        &[-3.0],
        &[0.0, 1.0],
        &[1.0],
        &[0.0, 0.2],
        budget(),
    )
    .unwrap();
    let initial = [0.07, -0.2];
    let initial_filter = [0.03];
    let mut expected_filter = [0.0];
    let mut expected_output = [0.0];
    reference
        .preview_into(
            &long,
            &initial,
            &initial_filter,
            &cosine,
            &sine,
            &mut expected_filter,
            &mut expected_output,
        )
        .unwrap();
    let mut short = source(0.001, &frequencies);
    let mut observer = PreparedToneObserver::new(
        &short,
        1,
        1,
        &[-3.0],
        &[0.0, 1.0],
        &[1.0],
        &[0.0, 0.2],
        budget(),
    )
    .unwrap();
    let mut x = initial;
    let mut filter = initial_filter;
    let mut next = [0.0; 2];
    let mut mean = [0.0];
    let mut next_filter = [0.0];
    let mut output = [0.0];
    TRACK.with(|c| c.set((true, 0)));
    for _ in 0..1000 {
        observer
            .preview_into(
                &short,
                &x,
                &filter,
                &cosine,
                &sine,
                &mut next_filter,
                &mut output,
            )
            .unwrap();
        short
            .step_into(&x, &cosine, &sine, &mut next, &mut mean)
            .unwrap();
        x = next;
        filter = next_filter;
        for i in 0..frequencies.len() {
            (cosine[i], sine[i]) = short.advance_force_phase(i, cosine[i], sine[i]).unwrap();
        }
    }
    let allocations = TRACK.with(|c| {
        let (_, n) = c.get();
        c.set((false, n));
        n
    });
    assert_eq!(allocations, 0);
    near(filter[0], expected_filter[0]);
    near(output[0], expected_output[0]);
}
