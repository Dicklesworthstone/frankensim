//! Independent scalar and full forced-trajectory checks of source loss.
use fs_phs::{
    PreparedStepError,
    linear_exact::{DenseExponentialBudget, PreparedSourceLoss},
};

fn budget() -> DenseExponentialBudget {
    DenseExponentialBudget {
        max_dimension: 256,
        max_dense_bytes: 64 * 1024 * 1024,
        max_multiply_terms: 100_000_000_000,
    }
}
fn close(a: f64, b: f64, scale: f64) {
    assert!((a - b).abs() <= 3e-12 * scale.max(1e-300), "{a:e} != {b:e}");
}

#[test]
fn positive_source_integral_matches_forced_analytic_heat_and_exact_silent_modes() {
    let (a, r, dt, x, u) = (3.0_f64, 0.7_f64, 0.4_f64, 3.5_f64, 0.9_f64);
    let prepared =
        PreparedSourceLoss::new_controlled(1, 1, &[-a], &[1.0], &[r.sqrt()], dt, budget(), || {
            false
        })
        .unwrap();
    let equilibrium = u / a;
    let amplitude = x - equilibrium;
    let expected = r
        * (equilibrium * equilibrium * dt
            + 2.0 * equilibrium * amplitude * (-(-a * dt).exp_m1()) / a
            + amplitude * amplitude * (-(-2.0 * a * dt).exp_m1()) / (2.0 * a));
    close(
        prepared.evaluate_controlled(&[x], &[u], || false).unwrap(),
        expected,
        expected,
    );
    assert_eq!(prepared.duration(), dt);

    // One conserved translation and one relaxing relative coordinate. The
    // original positive source nullspace survives a large nonzero reference.
    let relative = PreparedSourceLoss::new_controlled(
        2,
        0,
        &[-a, a, 0.0, 0.0],
        &[],
        &[r.sqrt(), -r.sqrt()],
        dt,
        budget(),
        || false,
    )
    .unwrap();
    assert_eq!(
        relative
            .evaluate_controlled(&[4096.0, 4096.0], &[], || false)
            .unwrap(),
        0.0
    );
    let tiny = [1e-100, 0.0];
    let expected = r * tiny[0] * tiny[0] * (-(-2.0 * a * dt).exp_m1()) / (2.0 * a);
    close(
        relative.evaluate_controlled(&tiny, &[], || false).unwrap(),
        expected,
        expected,
    );
}

#[test]
fn independent_loss_channels_follow_the_same_nonnormal_forced_motion() {
    let (r1, r2, omega, force, dt) = (0.2_f64, 0.8_f64, 4.0_f64, 0.9_f64, 0.7_f64);
    let a = [-r1, omega, -omega, -r2];
    let first = PreparedSourceLoss::new_controlled(
        2,
        1,
        &a,
        &[1.0, 0.0],
        &[r1.sqrt(), 0.0],
        dt,
        budget(),
        || false,
    )
    .unwrap();
    let second = PreparedSourceLoss::new_controlled(
        2,
        1,
        &a,
        &[1.0, 0.0],
        &[0.0, r2.sqrt()],
        dt,
        budget(),
        || false,
    )
    .unwrap();
    let initial = [0.7, -0.4];
    let actual = [
        first
            .evaluate_controlled(&initial, &[force], || false)
            .unwrap(),
        second
            .evaluate_controlled(&initial, &[force], || false)
            .unwrap(),
    ];
    // Independent RK4 integrates BOTH state and physical powers, never a
    // production cached form or a share of an inferred total decrement.
    let rhs = |y: [f64; 5]| {
        [
            a[0] * y[0] + a[1] * y[1] + force,
            a[2] * y[0] + a[3] * y[1],
            r1 * y[0] * y[0],
            r2 * y[1] * y[1],
            force * y[0],
        ]
    };
    let mut y = [initial[0], initial[1], 0.0, 0.0, 0.0];
    let h = dt / 4096.0;
    for _ in 0..4096 {
        let k1 = rhs(y);
        let k2 = rhs(core::array::from_fn(|i| y[i] + 0.5 * h * k1[i]));
        let k3 = rhs(core::array::from_fn(|i| y[i] + 0.5 * h * k2[i]));
        let k4 = rhs(core::array::from_fn(|i| y[i] + h * k3[i]));
        for i in 0..5 {
            y[i] += h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
        }
    }
    close(actual[0], y[2], y[2]);
    close(actual[1], y[3], y[3]);
    close(
        0.5 * (y[0] * y[0] + y[1] * y[1]) + actual.iter().sum::<f64>(),
        0.5 * (initial[0] * initial[0] + initial[1] * initial[1]) + y[4],
        1.0,
    );
}

#[test]
fn original_source_layout_budget_and_cancellation_cannot_publish_partial_heat() {
    let make =
        |b| PreparedSourceLoss::new_controlled(1, 1, &[-1.0], &[1.0], &[1.0], 0.1, b, || false);
    assert!(
        make(DenseExponentialBudget {
            max_dimension: 1,
            ..budget()
        })
        .is_err()
    );
    assert!(
        make(DenseExponentialBudget {
            max_dense_bytes: 0,
            ..budget()
        })
        .is_err()
    );
    assert!(
        make(DenseExponentialBudget {
            max_multiply_terms: 0,
            ..budget()
        })
        .is_err()
    );
    assert!(
        PreparedSourceLoss::new_controlled(
            1,
            1,
            &[-1.0],
            &[1.0],
            &[1.0, 0.0],
            0.1,
            budget(),
            || false
        )
        .is_err()
    );
    assert!(matches!(
        PreparedSourceLoss::new_controlled(1, 1, &[-1.0], &[1.0], &[1.0], 0.1, budget(), || true),
        Err(PreparedStepError::Cancelled)
    ));
    let mut polls = 0;
    assert!(matches!(
        PreparedSourceLoss::new_controlled(1, 1, &[-1.0], &[1.0], &[1.0], 0.1, budget(), || {
            polls += 1;
            polls == 2
        }),
        Err(PreparedStepError::Cancelled)
    ));
    let prepared = make(budget()).unwrap();
    assert!(
        prepared
            .evaluate_controlled(&[1.0], &[f64::NAN], || false)
            .is_err()
    );
    assert!(prepared.evaluate_controlled(&[], &[1.0], || false).is_err());
    assert!(matches!(
        prepared.evaluate_controlled(&[1.0], &[1.0], || true),
        Err(PreparedStepError::Cancelled)
    ));
    let expected = prepared
        .evaluate_controlled(&[1.0], &[1.0], || false)
        .unwrap();
    close(expected, 0.1, 0.1); // Forced equilibrium still dissipates actual input work.
    let empty =
        PreparedSourceLoss::new_controlled(1, 0, &[-1.0], &[], &[], 0.1, budget(), || false)
            .unwrap();
    assert_eq!(
        empty.evaluate_controlled(&[1.0], &[], || false).unwrap(),
        0.0
    );
}

#[test]
fn source_resistor_heat_uses_the_voltage_difference_including_near_equilibrium() {
    let (rate, conductance, dt) = (0.1_f64, 0.09_f64, 0.7_f64);
    let source = PreparedSourceLoss::new_with_feedthrough_controlled(
        1,
        1,
        &[-rate],
        &[rate],
        &[conductance.sqrt()],
        &[-conductance.sqrt()],
        dt,
        budget(),
        || false,
    )
    .unwrap();
    for (initial, voltage) in [(0.0, 4.0), (-3.0, -1.0), (1e16 + 2.0, 1e16)] {
        // C*dV/dt=G*(Vs-V). This physical resistor heat is independent of
        // the capacitor's energy datum and is not the Norton shunt G*V².
        let difference = initial - voltage;
        let expected =
            conductance * difference * difference * (-(-2.0 * rate * dt).exp_m1()) / (2.0 * rate);
        let mut scratch = [f64::NAN; 2];
        let actual = source
            .evaluate_into_controlled(&[initial], &[voltage], &mut scratch, || false)
            .unwrap();
        close(actual, expected, expected);
        assert_eq!(
            actual,
            source
                .evaluate_controlled(&[initial], &[voltage], || false)
                .unwrap()
        );
    }
    assert_eq!(
        source
            .evaluate_controlled(&[1e16], &[1e16], || false)
            .unwrap(),
        0.0
    );
}

#[test]
fn multiple_physical_source_rows_include_all_state_input_cross_terms() {
    let dt = 0.4;
    let source = PreparedSourceLoss::new_with_feedthrough_controlled(
        2,
        2,
        &[0.0; 4],
        &[1.0, -2.0, 0.5, 1.0],
        &[1.0, 2.0, -3.0, 0.5],
        &[-2.0, 1.0, 0.25, -4.0],
        dt,
        budget(),
        || false,
    )
    .unwrap();
    let (initial, inputs) = ([2.0, -4.0], [3.0, -1.0]);
    // x(t)=[2+5t, -4+0.5t]. The original two physical loss rows are
    // [-13+6t, -3.25-14.75t]; integrate their squares independently.
    let expected: f64 = [(-13.0, 6.0), (-3.25, -14.75)]
        .into_iter()
        .map(|(a, b)| dt * (a * a + a * b * dt + b * b * dt * dt / 3.0))
        .sum();
    close(
        source
            .evaluate_controlled(&initial, &inputs, || false)
            .unwrap(),
        expected,
        expected,
    );

    for (row, (a, b)) in [(-13.0, 6.0), (-3.25, -14.75)].into_iter().enumerate() {
        let signed = dt * (a + 0.5 * b * dt);
        close(
            source
                .signed_integral_into_controlled(row, &initial, &inputs, &mut [0.0; 4], || false)
                .unwrap(),
            signed,
            signed.abs(),
        );
    }

    let pure_source = PreparedSourceLoss::new_with_feedthrough_controlled(
        1,
        2,
        &[-2.0],
        &[0.0, 0.0],
        &[0.0],
        &[2.0, -3.0],
        dt,
        budget(),
        || false,
    )
    .unwrap();
    close(
        pure_source
            .evaluate_controlled(&[100.0], &[4.0, 5.0], || false)
            .unwrap(),
        49.0 * dt,
        49.0 * dt,
    );
    close(
        pure_source
            .signed_integral_into_controlled(0, &[100.0], &[4.0, 5.0], &mut [0.0; 2], || false)
            .unwrap(),
        -7.0 * dt,
        7.0 * dt,
    );
    assert_eq!(
        pure_source
            .evaluate_controlled(&[100.0], &[0.0, 0.0], || false)
            .unwrap(),
        0.0
    );
}

#[test]
fn source_feedthrough_preserves_zero_mapping_and_refuses_incomplete_or_cancelled_work() {
    let original =
        PreparedSourceLoss::new_controlled(1, 1, &[-3.0], &[1.0], &[0.7], 0.4, budget(), || false)
            .unwrap();
    let zero = PreparedSourceLoss::new_with_feedthrough_controlled(
        1,
        1,
        &[-3.0],
        &[1.0],
        &[0.7],
        &[0.0],
        0.4,
        budget(),
        || false,
    )
    .unwrap();
    assert_eq!(
        original
            .evaluate_controlled(&[3.5], &[0.9], || false)
            .unwrap()
            .to_bits(),
        zero.evaluate_controlled(&[3.5], &[0.9], || false)
            .unwrap()
            .to_bits()
    );
    for invalid in [&[][..], &[1.0, 2.0][..], &[f64::NAN][..]] {
        assert!(
            PreparedSourceLoss::new_with_feedthrough_controlled(
                1,
                1,
                &[-3.0],
                &[1.0],
                &[0.7],
                invalid,
                0.4,
                budget(),
                || false,
            )
            .is_err()
        );
    }
    assert!(
        PreparedSourceLoss::new_with_feedthrough_controlled(
            1,
            1,
            &[-3.0],
            &[1.0],
            &[0.7],
            &[0.2],
            0.4,
            DenseExponentialBudget {
                max_dense_bytes: 0,
                ..budget()
            },
            || false,
        )
        .is_err()
    );
    assert!(matches!(
        PreparedSourceLoss::new_with_feedthrough_controlled(
            1,
            1,
            &[-3.0],
            &[1.0],
            &[0.7],
            &[0.2],
            0.4,
            budget(),
            || true,
        ),
        Err(PreparedStepError::Cancelled)
    ));
    let source = PreparedSourceLoss::new_with_feedthrough_controlled(
        1,
        1,
        &[-3.0],
        &[1.0],
        &[0.7],
        &[0.2],
        0.4,
        budget(),
        || false,
    )
    .unwrap();
    let expected = source
        .evaluate_controlled(&[3.5], &[0.9], || false)
        .unwrap();
    let mut scratch = [0.0; 2];
    let mut polls = 0;
    assert!(matches!(
        source.evaluate_into_controlled(&[3.5], &[0.9], &mut scratch, || {
            polls += 1;
            polls == 3
        }),
        Err(PreparedStepError::Cancelled)
    ));
    assert_eq!(
        source
            .evaluate_into_controlled(&[3.5], &[0.9], &mut scratch, || false)
            .unwrap(),
        expected
    );
    assert!(
        source
            .evaluate_into_controlled(&[3.5], &[0.9], &mut [], || false)
            .is_err()
    );
}

#[test]
fn signed_source_integral_retains_small_regenerative_current_about_large_bus_voltage() {
    let dt = 0.3_f64;
    // Unit RC: I(t)=Vs-V(t)=-(V0-Vs)*exp(-t). Both the signed
    // charge and positive heat have independent, cancellation-safe primitives.
    let source = PreparedSourceLoss::new_with_feedthrough_controlled(
        1,
        1,
        &[-1.0],
        &[1.0],
        &[-1.0],
        &[1.0],
        dt,
        budget(),
        || false,
    )
    .unwrap();
    for (voltage, offset) in [(24.0, 1e-7), (24.0, 1e-10), (1e16, 2.0)] {
        let initial = voltage + offset;
        let difference = initial - voltage;
        let charge = difference * (-dt).exp_m1();
        let heat = -0.5 * difference * difference * (-2.0 * dt).exp_m1();
        let actual = source
            .signed_integral_into_controlled(0, &[initial], &[voltage], &mut [0.0; 2], || false)
            .unwrap();
        assert!(actual < 0.0);
        close(actual, charge, charge.abs());
        close(
            source
                .evaluate_controlled(&[initial], &[voltage], || false)
                .unwrap(),
            heat,
            heat,
        );
        assert_eq!(source.signed_integral_into_controlled(
            0, &[voltage], &[voltage], &mut [0.0; 2], || false,
        ).unwrap(), 0.0);
    }
}

#[test]
fn signed_source_selection_excludes_unrelated_overflow_and_preserves_retry() {
    // Signed state-only differences also need compensated initialization;
    // rounding these two products independently would return 0.25 instead.
    let difference = PreparedSourceLoss::new_controlled(
        2,
        0,
        &[0.0; 4],
        &[],
        &[0.1, -0.1],
        1.0,
        budget(),
        || false,
    )
    .unwrap();
    close(
        difference
            .signed_integral_into_controlled(0, &[1e16 + 2.0, 1e16], &[], &mut [0.0; 4], || false)
            .unwrap(),
        0.2,
        0.2,
    );
    let source = PreparedSourceLoss::new_controlled(
        2,
        0,
        &[0.0; 4],
        &[],
        &[1.0, 0.0, 0.0, 1.0],
        100.0,
        budget(),
        || false,
    )
    .unwrap();
    let initial = [1.0, 1e308];
    let mut scratch = [0.0; 4];
    // The other source has an unrepresentable integral; an exact zero
    // coefficient must prevent its factor row from contaminating row zero.
    close(
        source
            .signed_integral_into_controlled(0, &initial, &[], &mut scratch, || false)
            .unwrap(),
        100.0,
        100.0,
    );
    assert!(
        source
            .signed_integral_into_controlled(1, &initial, &[], &mut scratch, || false,)
            .is_err()
    );
    assert!(
        source
            .signed_integral_into_controlled(2, &initial, &[], &mut scratch, || false,)
            .is_err()
    );
    assert!(matches!(
        source.signed_integral_into_controlled(0, &initial, &[], &mut scratch, || true,),
        Err(PreparedStepError::Cancelled)
    ));
    close(
        source
            .signed_integral_into_controlled(0, &initial, &[], &mut scratch, || false)
            .unwrap(),
        100.0,
        100.0,
    );
}

#[test]
fn g1_weak_damper_heat_survives_large_lossless_motion_and_rescaling() {
    // Unit masses with one physical dashpot along g=[1,1,-1]:
    // v'=-c*g*(g.v), power=c*(g.v)^2. The exact relative velocity
    // r=g.v obeys r'=-3*c*r while both orthogonal translations are lossless.
    // No dense eigensolver or production integral supplies this reference.
    let g = [1.0_f64, 1.0, -1.0];
    for c in [1.0_f64 / 1024.0, 1.0, 1024.0] {
        let flow: Vec<_> = g
            .iter()
            .flat_map(|a| g.iter().map(move |b| -c * a * b))
            .collect();
        let source = g.map(|value| c.sqrt() * value);
        let dt = 0.125 / c;
        let prepared =
            PreparedSourceLoss::new_controlled(3, 0, &flow, &[], &source, dt, budget(), || false)
                .unwrap();
        for (background, relative) in [
            (1.0, 1.0e-12),
            (1.0, 2.0_f64.powi(-55)),
            (1.0, -2.0_f64.powi(-55)),
            (1.0e16, 1.0),
        ] {
            let state = [background, relative, background];
            let expected = relative * relative * (-(-6.0 * c * dt).exp_m1()) / 6.0;
            let mut scratch = [f64::NAN; 6];
            let heat = prepared
                .evaluate_into_controlled(&state, &[], &mut scratch, || false)
                .unwrap();
            assert!(
                heat > 0.0,
                "a resolved physical velocity difference lost its heat"
            );
            close(heat, expected, expected);
            assert_eq!(
                heat,
                prepared.evaluate_controlled(&state, &[], || false).unwrap()
            );
        }
        assert_eq!(
            prepared
                .evaluate_controlled(&[1.0, 0.0, 1.0], &[], || false)
                .unwrap(),
            0.0
        );
    }
}

#[test]
fn g1_weak_conservative_coupling_must_reach_the_physical_dashpot() {
    // H=|x|^2/2; the first three coordinates couple conservatively through
    // g=[1,1,-1] to the fourth, whose damping rate is one:
    // r'= -3*y, y'=r-y, r=g.x. Only y^2 is dissipated.
    // Initially y=0, so all heat must come from the small acceleration that
    // an ordinary g.x sum can erase between large cancelling coordinates.
    let flow = [
        0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, -1.0, -1.0,
    ];
    let dt = 1.0_f64;
    let prepared = PreparedSourceLoss::new_controlled(
        4,
        0,
        &flow,
        &[],
        &[0.0, 0.0, 0.0, 1.0],
        dt,
        budget(),
        || false,
    )
    .unwrap();
    let frequency = (3.0_f64 - 0.25).sqrt();
    // y(t)=r0*exp(-t/2)*sin(frequency*t)/frequency. Integrate y^2
    // analytically; no subtraction of the large conserved background energy.
    let twice_frequency = 2.0 * frequency;
    let cosine_integral = (1.0
        + (-dt).exp()
            * (-(twice_frequency * dt).cos() + twice_frequency * (twice_frequency * dt).sin()))
        / (1.0 + twice_frequency * twice_frequency);
    let unit_heat = 0.5 * (-(-dt).exp_m1() - cosine_integral) / (frequency * frequency);
    assert!(unit_heat > 0.0);
    for amplitude in [1.0e-12, 2.0_f64.powi(-55), -2.0_f64.powi(-55)] {
        let expected = amplitude * amplitude * unit_heat;
        let heat = prepared
            .evaluate_into_controlled(&[1.0, amplitude, 1.0, 0.0], &[], &mut [0.0; 8], || false)
            .unwrap();
        assert!(
            heat > 0.0,
            "conservative cross coupling must drive the actual damper"
        );
        close(heat, expected, expected);
    }
    assert_eq!(
        prepared
            .evaluate_controlled(&[1.0, 0.0, 1.0, 0.0], &[], || false)
            .unwrap(),
        0.0
    );
}

#[test]
fn g3_state_only_heat_preserves_channel_sums_and_explicit_zero_feedthrough() {
    // Two physical dashpots share the same g.v, with c1=1/4 and c2=3/4.
    // Their sum produces A=-g*g^T; each heat is its own c_i times
    // the same analytic squared-velocity integral.
    let g = [1.0_f64, 1.0, -1.0];
    let flow: Vec<_> = g
        .iter()
        .flat_map(|a| g.iter().map(move |b| -a * b))
        .collect();
    let dt = 0.125_f64;
    let sources: Vec<_> = [0.25_f64, 0.75]
        .iter()
        .flat_map(|c| g.iter().map(move |value| c.sqrt() * value))
        .collect();
    let combined =
        PreparedSourceLoss::new_controlled(3, 1, &flow, &[0.0; 3], &sources, dt, budget(), || {
            false
        })
        .unwrap();
    let explicit_zero = PreparedSourceLoss::new_with_feedthrough_controlled(
        3,
        1,
        &flow,
        &[0.0; 3],
        &sources,
        &[0.0; 2],
        dt,
        budget(),
        || false,
    )
    .unwrap();
    let relative = 2.0_f64.powi(-55);
    let state = [1.0, relative, 1.0];
    let expected = relative * relative * (-(-6.0 * dt).exp_m1()) / 6.0;
    let heat = combined
        .evaluate_controlled(&state, &[7.0], || false)
        .unwrap();
    close(heat, expected, expected);
    assert_eq!(
        heat,
        explicit_zero
            .evaluate_controlled(&state, &[7.0], || false)
            .unwrap()
    );
    let mut separate = 0.0;
    for (index, c) in [0.25_f64, 0.75].into_iter().enumerate() {
        let channel = PreparedSourceLoss::new_controlled(
            3,
            0,
            &flow,
            &[],
            &sources[3 * index..3 * (index + 1)],
            dt,
            budget(),
            || false,
        )
        .unwrap();
        let value = channel.evaluate_controlled(&state, &[], || false).unwrap();
        close(value, c * expected, c * expected);
        separate += value;
    }
    close(heat, separate, expected);
}
