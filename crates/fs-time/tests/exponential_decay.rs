//! G1/G3: analytical forced decay, spectral scaling, and transactional stepping.

use fs_time::stiff::{ExpEuler, ExpEulerError};

fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-300),
        "{actual:e} != {expected:e}"
    );
}

fn zero(_: &[f64], output: &mut [f64]) {
    output.fill(0.0);
}

#[test]
fn constant_source_has_the_exact_decay_and_equilibrium() {
    let method = ExpEuler::try_new(&[-0.7], 1, 0.2).unwrap();
    let mut workspace = method.workspace();
    let source = |_: &[f64], out: &mut [f64]| out[0] = 0.3;
    let mut state = [1.4];
    for step in 1..=40 {
        method
            .step_checked(&mut state, &source, &mut workspace)
            .unwrap();
        let equilibrium = 0.3 / 0.7;
        let exact = equilibrium + (1.4 - equilibrium) * (-0.14 * f64::from(step)).exp();
        relative(state[0], exact, 8e-15);
    }
    let mut convenience = [1.4];
    let mut checked = convenience;
    method.step(&mut convenience, &source);
    method
        .step_checked(&mut checked, &source, &mut workspace)
        .unwrap();
    assert_eq!(convenience, checked);
}

#[test]
fn overflowing_rate_duration_preserves_the_finite_forced_equilibrium() {
    // The old h*expm1(h*lambda)/(h*lambda) evaluated this coefficient as
    // zero because h*lambda=-infinity, although -source/lambda is finite.
    let method = ExpEuler::try_new(&[-1e200], 1, 1e200).unwrap();
    let mut state = [3.0];
    method
        .step_checked(
            &mut state,
            &|_, out| out[0] = 1e200,
            &mut method.workspace(),
        )
        .unwrap();
    relative(state[0], 1.0, 3e-15);

    // Very short times retain the finite source impulse as well.
    let method = ExpEuler::try_new(&[-1e300], 1, 1e-300).unwrap();
    state = [0.0];
    method.step(&mut state, &|_, out| out[0] = 1e300);
    relative(state[0], -(-1.0_f64).exp_m1(), 3e-15);
}

#[test]
fn homogeneous_tail_survives_an_underflowing_retention_coefficient() {
    // Independent 90-digit decimal values of 1e200*exp(-800) and exp(-720).
    for (duration, expected) in [
        (800.0, 3.667_874_584_177_687e-148),
        (720.0, 2.032_230_802_424_293_2e-113),
    ] {
        let method = ExpEuler::try_new(&[-1.0], 1, duration).unwrap();
        for sign in [-1.0, 1.0] {
            let mut state = [sign * 1e200];
            method.step(&mut state, &zero);
            relative(state[0], sign * expected, 4e-15);
        }
    }
}

#[test]
fn large_finite_coupled_operator_keeps_both_physical_decay_rates() {
    // Sum/difference modes have rates 0.5e308 and 1.5e308. An unscaled
    // Jacobi diagonal-sum overflow incorrectly skipped this rotation.
    let method = ExpEuler::try_new(&[-1e308, 5e307, 5e307, -1e308], 2, 1e-308).unwrap();
    let mut state = [1.0, 0.0];
    method.step(&mut state, &zero);
    let slow = (-0.5_f64).exp();
    let fast = (-1.5_f64).exp();
    relative(state[0], 0.5 * (slow + fast), 5e-15);
    relative(state[1], 0.5 * (slow - fast), 5e-15);
}

#[test]
fn zero_rates_and_tiny_rate_duration_retain_the_source_limit() {
    let method = ExpEuler::try_new(&[0.0, 0.0, 0.0, -2.0], 2, 0.25).unwrap();
    let mut state = [1.0, 3.0];
    method.step(&mut state, &|_, out| out.copy_from_slice(&[4.0, 2.0]));
    relative(state[0], 2.0, 3e-15);
    relative(state[1], 1.0 + 2.0 * (-0.5_f64).exp(), 3e-15);

    // h*lambda underflows, but h*source is an ordinary representable value.
    let method = ExpEuler::try_new(&[-1e-200], 1, 1e-200).unwrap();
    let mut state = [0.0];
    method.step(&mut state, &|_, out| out[0] = 1e200);
    relative(state[0], 1.0, 3e-15);
}

#[test]
fn constant_forcing_is_invariant_under_time_partition() {
    let a = [-2.0, 0.7, 0.7, -3.0];
    let source = |_: &[f64], out: &mut [f64]| out.copy_from_slice(&[0.3, -0.7]);
    let whole = ExpEuler::try_new(&a, 2, 1.2).unwrap();
    let part = ExpEuler::try_new(&a, 2, 0.03).unwrap();
    let mut full = [0.8, -0.2];
    let mut split = full;
    whole.step(&mut full, &source);
    let mut workspace = part.workspace();
    for _ in 0..40 {
        part.step_checked(&mut split, &source, &mut workspace)
            .unwrap();
    }
    for (actual, expected) in split.into_iter().zip(full) {
        relative(actual, expected, 3e-14);
    }
}

#[test]
fn preparation_refuses_invalid_or_unresolved_spectra() {
    for h in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            ExpEuler::try_new(&[-1.0], 1, h),
            Err(ExpEulerError::InvalidInput(_))
        ));
    }
    assert!(matches!(
        ExpEuler::try_new(&[-1.0], 2, 1.0),
        Err(ExpEulerError::Dimension {
            expected: 4,
            actual: 1
        })
    ));
    assert!(ExpEuler::try_new(&[], usize::MAX, 1.0).is_err());
    assert!(ExpEuler::try_new(&[f64::NAN], 1, 1.0).is_err());
    assert!(ExpEuler::try_new(&[f64::INFINITY], 1, 1.0).is_err());
    assert!(matches!(
        ExpEuler::try_new(&[-1.0, 0.3, 0.1, -2.0], 2, 1.0),
        Err(ExpEulerError::NonSymmetric { row: 0, column: 1 })
    ));
    assert!(matches!(
        ExpEuler::try_new(&[1000.0], 1, 1.0),
        Err(ExpEulerError::Spectrum(_))
    ));
    assert!(matches!(
        ExpEuler::try_new(&[-1e300, 0.0, 0.0, -1e-300], 2, 1.0),
        Err(ExpEulerError::Spectrum(_))
    ));
}

#[test]
fn checked_refusals_leave_the_entire_state_unchanged() {
    let method = ExpEuler::try_new(&[-1.0, 0.0, 0.0, -2.0], 2, 0.2).unwrap();
    let mut workspace = method.workspace();
    let initial = [0.7, -0.4];
    let mut state = initial;
    assert_eq!(
        method.step_checked(&mut state, &|_, out| out[0] = 0.0, &mut workspace),
        Err(ExpEulerError::NonFiniteSource { index: 1 })
    );
    assert_eq!(state, initial);
    assert_eq!(
        method.step_checked(
            &mut state,
            &|_, out| out.fill(f64::INFINITY),
            &mut workspace
        ),
        Err(ExpEulerError::NonFiniteSource { index: 0 })
    );
    assert_eq!(state, initial);
    let other = ExpEuler::new(&[-1.0], 1, 0.2);
    assert!(matches!(
        method.step_checked(&mut state, &zero, &mut other.workspace()),
        Err(ExpEulerError::Dimension {
            expected: 2,
            actual: 1
        })
    ));
    assert_eq!(state, initial);
    let mut short = [0.7];
    assert!(matches!(
        method.step_checked(&mut short, &zero, &mut workspace),
        Err(ExpEulerError::Dimension {
            expected: 2,
            actual: 1
        })
    ));
    assert_eq!(short, [0.7]);

    let method = ExpEuler::new(&[0.0; 4], 2, 2.0);
    assert_eq!(
        method.step_checked(
            &mut state,
            &|_, out| out.copy_from_slice(&[1.0, f64::MAX]),
            &mut method.workspace(),
        ),
        Err(ExpEulerError::Arithmetic)
    );
    assert_eq!(state, initial);
    state[1] = f64::NAN;
    assert_eq!(
        method.step_checked(&mut state, &zero, &mut method.workspace()),
        Err(ExpEulerError::NonFiniteState { index: 1 })
    );
    assert_eq!(state[0], initial[0]);
    assert!(state[1].is_nan());
}
