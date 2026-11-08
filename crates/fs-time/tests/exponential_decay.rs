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

fn path_rates(n: usize, rate: f64) -> Vec<f64> {
    let mut rates = vec![0.0; n * n];
    for i in 0..n - 1 {
        rates[i * n + i + 1] = rate;
        rates[(i + 1) * n + i] = rate;
    }
    rates
}

#[test]
fn insulated_diffusion_retains_equilibrium_in_long_stiff_windows() {
    for n in [3, 5, 9] {
        let rates = path_rates(n, 1.0);
        let mut generator = rates.clone();
        for i in 0..n {
            generator[i * n + i] = -rates[i * n..(i + 1) * n].iter().sum::<f64>();
        }
        // A numerical zero eigenvalue used to become spurious growth or loss
        // here. Both the existing matrix API and authored-rate API must retain
        // the analytic constant mode, independently of the stiff duration.
        let matrix = ExpEuler::try_new(&generator, n, 1e16).unwrap();
        let exchange = ExpEuler::try_new_diffusion(&rates, n, 1e16).unwrap();
        let mut from_matrix: Vec<_> = (0..n).map(|i| (i + 1) as f64).collect();
        let mut from_rates = from_matrix.clone();
        matrix.step(&mut from_matrix, &zero);
        exchange.step(&mut from_rates, &zero);
        assert_eq!(from_matrix, from_rates);
        let mean = 0.5 * (n + 1) as f64;
        for &value in &from_matrix {
            relative(value, mean, 3e-14);
        }
        relative(from_matrix.iter().sum(), n as f64 * mean, 3e-14);

        let mut uniform = vec![2.5; n];
        matrix.step(&mut uniform, &zero);
        for value in uniform {
            relative(value, 2.5, 3e-14);
        }
    }
}

#[test]
fn disconnected_insulated_regions_keep_separate_material_balances() {
    let n = 6;
    let mut rates = vec![0.0; n * n];
    for (i, j) in [(0, 2), (1, 3), (3, 5)] {
        rates[i * n + j] = 1.0;
        rates[j * n + i] = 1.0;
    }
    let method = ExpEuler::try_new_diffusion(&rates, n, 1e16).unwrap();
    let mut state = [2.0, 1.0, 6.0, 3.0, -7.0, 5.0];
    method.step(&mut state, &zero);
    for (value, expected) in state.into_iter().zip([4.0, 3.0, 4.0, 3.0, -7.0, 3.0]) {
        relative(value, expected, 3e-14);
    }
}

#[test]
fn balanced_sources_do_not_create_a_mean_in_stiff_diffusion() {
    let method = ExpEuler::try_new_diffusion(&path_rates(3, 1.0), 3, 1e16).unwrap();
    let mut state = [3.0; 3];
    method.step(&mut state, &|_, out| out.copy_from_slice(&[0.4, -0.4, 0.0]));
    // At equilibrium u0-u1=0.4, u1=u2, and total material remains nine.
    let lower = (9.0 - 0.4) / 3.0;
    for (value, expected) in state.into_iter().zip([lower + 0.4, lower, lower]) {
        relative(value, expected, 3e-14);
    }
    relative(state.iter().sum(), 9.0, 3e-14);

    // A real component source advances its mean by exactly the held source
    // integral; it must neither disappear nor be treated as an eigenvalue.
    state = [1.0; 3];
    method.step(&mut state, &|_, out| out.fill(1e-16));
    for value in state {
        relative(value, 2.0, 3e-14);
    }
}

#[test]
fn multiscale_component_balances_survive_stiff_source_integration() {
    let n = 6;
    let large = (1_u64 << 53) as f64;
    let tiny = 1.0 / large;
    let mut rates = vec![(1_u64 << 60) as f64; n * n];
    for i in 0..n {
        rates[i * n + i] = 0.0;
    }
    let method = ExpEuler::try_new_diffusion(&rates, n, large).unwrap();
    let mut workspace = method.workspace();
    // All six values are exactly binary64. Neumaier's compensated sum
    // leaves -2^-53 in this order despite their exactly zero real sum.
    // Integrating that residual over h=2^53 used to remove one unit.
    let balanced = [large, 1.0, tiny, -large, -1.0, -tiny];
    for offset in 0..n {
        for net in [0.0, tiny, -tiny] {
            let mut source = balanced;
            source[5] += net;
            source.rotate_left(offset);
            let mut state = [0.0; 6];
            method
                .step_checked(
                    &mut state,
                    &|_, out| out.copy_from_slice(&source),
                    &mut workspace,
                )
                .unwrap();
            let expected = large * net;
            assert!(
                (state.iter().sum::<f64>() - expected).abs() <= 2e-14,
                "offset={offset}, net={net:e}, state={state:?}",
            );
        }
        let mut state = balanced;
        state.rotate_left(offset);
        method
            .step_checked(&mut state, &zero, &mut workspace)
            .unwrap();
        assert_eq!(state, [0.0; 6]);
    }

    // An intermediate expansion that overflows is refused atomically. Scratch
    // from that refusal must not contaminate a later representable balance.
    let mut state = [1.0; 6];
    assert_eq!(
        method.step_checked(
            &mut state,
            &|_, out| out.copy_from_slice(&[f64::MAX, f64::MAX, -f64::MAX, -f64::MAX, 0.0, 0.0],),
            &mut workspace,
        ),
        Err(ExpEulerError::Arithmetic),
    );
    assert_eq!(state, [1.0; 6]);
    method
        .step_checked(
            &mut state,
            &|_, out| out.copy_from_slice(&balanced),
            &mut workspace,
        )
        .unwrap();
    relative(state.iter().sum(), 6.0, 3e-14);
}

#[test]
fn driven_exchange_matches_the_independent_two_cell_solution() {
    let rate: f64 = 3.0;
    let h: f64 = 0.2;
    let method = ExpEuler::try_new_diffusion(&path_rates(2, rate), 2, h).unwrap();
    let mut state = [1.0, 0.0];
    let source = |_: &[f64], out: &mut [f64]| out.copy_from_slice(&[0.7, -0.1]);
    let mut workspace = method.workspace();
    for step in 1..=20 {
        method
            .step_checked(&mut state, &source, &mut workspace)
            .unwrap();
        let time = f64::from(step) * h;
        let mean = 0.5 + 0.3 * time;
        let difference =
            0.8 / (2.0 * rate) + (1.0 - 0.8 / (2.0 * rate)) * (-2.0 * rate * time).exp();
        relative(state[0], mean + 0.5 * difference, 8e-14);
        relative(state[1], mean - 0.5 * difference, 8e-14);
    }
}

#[test]
fn genuine_small_diffusion_growth_and_decay_rates_are_retained() {
    let method = ExpEuler::try_new_diffusion(&path_rates(2, 1e-18), 2, 1e18).unwrap();
    let mut state = [1.0, 0.0];
    method.step(&mut state, &zero);
    let difference = (-2.0_f64).exp();
    relative(state[0], 0.5 * (1.0 + difference), 8e-15);
    relative(state[1], 0.5 * (1.0 - difference), 8e-15);

    for small in [-1e-18, 1e-18] {
        let generic = ExpEuler::try_new(&[-1.0, 0.0, 0.0, small], 2, 1e18).unwrap();
        state = [0.0, 1.0];
        generic.step(&mut state, &zero);
        assert_eq!(state[0], 0.0);
        relative(state[1], (small * 1e18).exp(), 4e-15);
    }
}

#[test]
fn diffusion_admission_rejects_nonreciprocal_or_invalid_models() {
    for rates in [
        [0.0, -1.0, -1.0, 0.0],
        [0.0, 1.0, 2.0, 0.0],
        [-1.0, 1.0, 1.0, -1.0],
        [0.0, f64::NAN, f64::NAN, 0.0],
        [0.0, f64::INFINITY, f64::INFINITY, 0.0],
    ] {
        assert!(ExpEuler::try_new_diffusion(&rates, 2, 1.0).is_err());
    }
    assert!(ExpEuler::try_new_diffusion(&[], 0, 1.0).is_err());
    assert!(ExpEuler::try_new_diffusion(&[0.0], 2, 1.0).is_err());
    assert!(ExpEuler::try_new_diffusion(&[], usize::MAX, 1.0).is_err());
    assert!(ExpEuler::try_new_diffusion(&[0.0], 1, 0.0).is_err());
    let isolated = ExpEuler::try_new_diffusion(&[0.0], 1, 0.5).unwrap();
    let mut state = [2.0];
    isolated.step(&mut state, &|_, out| out[0] = 3.0);
    assert_eq!(state, [3.5]);
}

#[test]
fn unequal_reservoirs_match_the_analytic_heat_and_power_balance() {
    // C=(2,8) J/K, G=3 W/K and P=(7,-2) W. The mean heats at
    // (P0+P1)/(C0+C1)=0.5 K/s; the difference obeys d'=-1.875*d+3.75.
    let capacities = [2.0, 8.0];
    let conductances = [0.0, 3.0, 3.0, 0.0];
    let method = ExpEuler::try_new_capacity_diffusion(&conductances, &capacities, 0.2).unwrap();
    let initial = [350.0, 290.0];
    let mut state = initial;
    // Workspaces remain interchangeable with ordinary same-dimension steppers.
    let ordinary = ExpEuler::try_new(&[-1.0, 0.0, 0.0, -2.0], 2, 0.2).unwrap();
    let mut workspace = ordinary.workspace();
    let power = |physical: &[f64], out: &mut [f64]| {
        assert!(physical.iter().all(|value| (280.0..360.0).contains(value)));
        out.copy_from_slice(&[7.0, -2.0]);
    };
    for step in 1..=20 {
        method
            .step_checked(&mut state, &power, &mut workspace)
            .unwrap();
        let time = 0.2 * f64::from(step);
        let mean = 302.0 + 0.5 * time;
        let difference = 2.0 + 58.0 * (-1.875 * time).exp();
        relative(state[0], mean + 0.8 * difference, 8e-14);
        relative(state[1], mean - 0.2 * difference, 8e-14);
        relative(2.0 * state[0] + 8.0 * state[1], 3020.0 + 5.0 * time, 8e-14);
    }
    let whole = ExpEuler::try_new_capacity_diffusion(&conductances, &capacities, 4.0).unwrap();
    let mut direct = initial;
    whole.step(&mut direct, &power);
    for (a, b) in direct.into_iter().zip(state) {
        relative(a, b, 8e-14);
    }
}

#[test]
fn capacity_diffusion_retains_separate_weighted_equilibria() {
    let capacities = [1.0, 3.0, 8.0, 2.0, 5.0, 11.0];
    let n = capacities.len();
    let mut conductances = vec![0.0; n * n];
    for (i, j, value) in [(0, 2, 1.0), (1, 3, 2.0), (3, 5, 4.0)] {
        conductances[i * n + j] = value;
        conductances[j * n + i] = value;
    }
    let method = ExpEuler::try_new_capacity_diffusion(&conductances, &capacities, 1e16).unwrap();
    let mut state = [10.0, 20.0, 30.0, 40.0, 7.0, 60.0];
    method.step(&mut state, &zero);
    for (value, expected) in
        state
            .into_iter()
            .zip([250.0 / 9.0, 50.0, 250.0 / 9.0, 50.0, 7.0, 50.0])
    {
        relative(value, expected, 8e-14);
    }
    relative(state[0] + 8.0 * state[2], 250.0, 8e-14);
    relative(
        3.0 * state[1] + 2.0 * state[3] + 11.0 * state[5],
        800.0,
        8e-14,
    );

    // The isolated cell's source is an extensive power, not a temperature rate.
    let isolated = ExpEuler::try_new_capacity_diffusion(&[0.0], &[5.0], 0.5).unwrap();
    let mut one = [7.0];
    isolated.step(&mut one, &|_, out| out[0] = 3.0);
    relative(one[0], 7.3, 8e-15);

    // h*P alone underflows, but the transformed impulse and physical state do not.
    let tiny = ExpEuler::try_new_capacity_diffusion(&[0.0], &[1e-200], 1e-200).unwrap();
    let mut one = [0.0];
    tiny.step(&mut one, &|_, out| out[0] = 1e-200);
    relative(one[0], 1e-200, 8e-15);
}

#[test]
fn weighted_totals_retain_product_tails_and_balanced_physical_power() {
    let epsilon = 2.0_f64.powi(-27);
    let capacities = [1.0 + epsilon, 1.0];
    let method =
        ExpEuler::try_new_capacity_diffusion(&[0.0, 1.0, 1.0, 0.0], &capacities, 1e16).unwrap();
    // The rounded product C0*u0 is 1.0, but its exact value is 1-2^-54.
    // Retaining the FMA residual is necessary to keep this real stored amount.
    let mut state = [1.0 - epsilon, -1.0];
    assert_eq!(capacities[0] * state[0], 1.0);
    method.step(&mut state, &zero);
    let expected = -(epsilon * epsilon) / (capacities[0] + capacities[1]);
    for value in state {
        relative(value, expected, 8e-14);
    }

    let large = (1_u64 << 53) as f64;
    let tiny = 1.0 / large;
    let capacities = [1.0, 2.0, 3.0, 5.0, 7.0, 11.0];
    let mut conductances = vec![(1_u64 << 60) as f64; 36];
    for i in 0..6 {
        conductances[6 * i + i] = 0.0;
    }
    let method = ExpEuler::try_new_capacity_diffusion(&conductances, &capacities, large).unwrap();
    let mut workspace = method.workspace();
    for net in [0.0, tiny, -tiny] {
        let power = [large, 1.0, tiny, -large, -1.0, -tiny + net];
        let mut state = [0.0; 6];
        method
            .step_checked(
                &mut state,
                &|_, out| out.copy_from_slice(&power),
                &mut workspace,
            )
            .unwrap();
        let amount: f64 = capacities.iter().zip(state).map(|(c, u)| c * u).sum();
        assert!(
            (amount - large * net).abs() < 2e-13,
            "amount={amount:e}, net={net:e}"
        );
    }
}

#[test]
fn capacity_units_and_unit_storage_leave_physical_trajectories_unchanged() {
    let conductances = path_rates(3, 2.0);
    let ordinary = ExpEuler::try_new_diffusion(&conductances, 3, 0.4).unwrap();
    let unit = ExpEuler::try_new_capacity_diffusion(&conductances, &[1.0; 3], 0.4).unwrap();
    let mut a = [2.0, 4.0, 7.0];
    let mut b = a;
    let forcing = |_: &[f64], out: &mut [f64]| out.copy_from_slice(&[0.7, -0.2, 0.1]);
    ordinary.step(&mut a, &forcing);
    unit.step(&mut b, &forcing);
    assert_eq!(a, b);

    for scale in [1e-120, 1.0, 1e120] {
        let method = ExpEuler::try_new_capacity_diffusion(
            &[0.0, 3.0 * scale, 3.0 * scale, 0.0],
            &[2.0 * scale, 8.0 * scale],
            0.2,
        )
        .unwrap();
        let mut state = [350.0, 290.0];
        method.step(&mut state, &|_, out| {
            out.copy_from_slice(&[7.0 * scale, -2.0 * scale])
        });
        let mean = 302.1;
        let difference = 2.0 + 58.0 * (-0.375_f64).exp();
        relative(state[0], mean + 0.8 * difference, 5e-14);
        relative(state[1], mean - 0.2 * difference, 5e-14);
    }
}

#[test]
fn capacity_diffusion_refuses_invalid_or_unrepresentable_physics_atomically() {
    let conductances = [0.0, 3.0, 3.0, 0.0];
    for capacities in [
        [0.0, 1.0],
        [-1.0, 1.0],
        [f64::NAN, 1.0],
        [f64::INFINITY, 1.0],
        [f64::MAX, f64::MAX],
    ] {
        assert!(ExpEuler::try_new_capacity_diffusion(&conductances, &capacities, 0.2).is_err());
    }
    assert!(ExpEuler::try_new_capacity_diffusion(&[], &[], 0.2).is_err());
    assert!(ExpEuler::try_new_capacity_diffusion(&[0.0], &[2.0, 8.0], 0.2).is_err());
    assert!(ExpEuler::try_new_capacity_diffusion(&conductances, &[2.0, 8.0], 0.0).is_err());
    for conductances in [
        [1.0, 3.0, 3.0, 0.0],
        [0.0, -3.0, -3.0, 0.0],
        [0.0, 3.0, 4.0, 0.0],
        [0.0, f64::INFINITY, f64::INFINITY, 0.0],
    ] {
        assert!(ExpEuler::try_new_capacity_diffusion(&conductances, &[2.0, 8.0], 0.2).is_err());
    }
    assert!(
        ExpEuler::try_new_capacity_diffusion(&[0.0, 1e-100, 1e-100, 0.0], &[1e300, 1e300], 0.2,)
            .is_err()
    );

    let method = ExpEuler::try_new_capacity_diffusion(&conductances, &[0.25, 1.0], 0.01).unwrap();
    let mut workspace = method.workspace();
    let mut state = [2.0, 3.0];
    let before = state;
    assert_eq!(
        method.step_checked(
            &mut state,
            &|_, out| out.copy_from_slice(&[f64::MAX, 0.0]),
            &mut workspace,
        ),
        Err(ExpEulerError::Arithmetic)
    );
    assert_eq!(state, before);
    method
        .step_checked(&mut state, &zero, &mut workspace)
        .unwrap();
    relative(0.25 * state[0] + state[1], 3.5, 8e-14);

    let method = ExpEuler::try_new_capacity_diffusion(&conductances, &[2.0, 8.0], 0.2).unwrap();
    let mut overflowing = [f64::MAX, 1.0];
    assert_eq!(
        method.step_checked(&mut overflowing, &zero, &mut method.workspace(),),
        Err(ExpEulerError::Arithmetic)
    );
    assert_eq!(overflowing, [f64::MAX, 1.0]);

    // The product itself is normal, but its exact low bit is below 2^-1074.
    // Losing that residual must not silently change a conserved amount.
    let tiny_capacity = (4.0 * f64::MIN_POSITIVE) * (1.0 + f64::EPSILON);
    let method = ExpEuler::try_new_capacity_diffusion(&[0.0], &[tiny_capacity], 1.0).unwrap();
    let mut unresolved = [1.0 + f64::EPSILON];
    assert_eq!(
        method.step_checked(&mut unresolved, &zero, &mut method.workspace(),),
        Err(ExpEulerError::Arithmetic)
    );
    assert_eq!(unresolved, [1.0 + f64::EPSILON]);
}
