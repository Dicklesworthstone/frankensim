//! Exact reduced-time ramps and atomic section updates: Gauntlet G2
//! analytic convolution references and G3 partition/refusal metamorphics.

use fs_airfoil::indicial::{IndicialKernel, KUSSNER_2POLE, MAX_DS, WAGNER_JONES};
use fs_airfoil::unsteady::{
    DuhamelState, UNSTEADY_SECTION_V1, UnsteadySectionSpec, UnsteadySectionState,
};

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}, tolerance {tolerance:.3e}"
    );
}

// Integrate the analytic indicial step response against a constant input
// slope. This reference uses platform exp_m1 and the integrated response,
// independently of the production deficiency recurrence.
fn ramp_response(kernel: &IndicialKernel, start: f64, end: f64, duration: f64) -> f64 {
    let mut integral = duration;
    for i in 0..2 {
        let lost = kernel.a[i] / kernel.b[i] * -(-kernel.b[i] * duration).exp_m1();
        integral -= lost;
    }
    start + (end - start) / duration * integral
}

fn state_bits(state: &UnsteadySectionState) -> [u64; 7] {
    [
        state.motion.x[0].to_bits(),
        state.motion.x[1].to_bits(),
        state.motion.input_prev.to_bits(),
        state.gust.x[0].to_bits(),
        state.gust.x[1].to_bits(),
        state.gust.input_prev.to_bits(),
        state.separation.f.to_bits(),
    ]
}

#[test]
fn both_kernels_match_the_integrated_step_response_for_linear_ramps() {
    for kernel in [WAGNER_JONES, KUSSNER_2POLE] {
        for duration in [0.025, 0.5, 3.0, MAX_DS] {
            let mut state = DuhamelState::trim(-0.05);
            state.advance_linear(&kernel, duration, 0.15).unwrap();
            close(
                state.effective(),
                ramp_response(&kernel, -0.05, 0.15, duration),
                3e-15,
            );
            assert_eq!(state.input_prev.to_bits(), 0.15_f64.to_bits());
        }
        // This is a distinct physical interpolation from the legacy
        // midpoint input step; a large admitted interval distinguishes them.
        let mut ramp = DuhamelState::trim(0.0);
        let mut midpoint = ramp;
        ramp.advance_linear(&kernel, 3.0, 0.2).unwrap();
        midpoint.advance(&kernel, 3.0, 0.2).unwrap();
        assert!((ramp.effective() - midpoint.effective()).abs() > 1e-4);
    }
}

#[test]
fn unequal_ramp_partitions_compose_with_history_and_a_slope_reversal() {
    for kernel in [WAGNER_JONES, KUSSNER_2POLE] {
        let mut whole = DuhamelState {
            x: [0.03125, -0.015625],
            input_prev: -0.0625,
        };
        let mut split = whole;
        for (end, duration) in [(0.25, 8.0), (-0.125, 2.0), (0.0, 0.5)] {
            let start = whole.input_prev;
            whole.advance_linear(&kernel, duration, end).unwrap();
            let mut previous = 0.0;
            for fraction in [0.125, 0.375, 0.5, 1.0] {
                split
                    .advance_linear(
                        &kernel,
                        duration * (fraction - previous),
                        start + fraction * (end - start),
                    )
                    .unwrap();
                previous = fraction;
            }
            close(split.x[0], whole.x[0], 5e-16);
            close(split.x[1], whole.x[1], 5e-16);
            close(split.effective(), whole.effective(), 1e-15);
            assert_eq!(split.input_prev.to_bits(), whole.input_prev.to_bits());
        }
    }
}

#[test]
fn attached_section_sums_independent_motion_and_gust_ramps() {
    let spec = UNSTEADY_SECTION_V1;
    let mut section = UnsteadySectionState::trim(&spec, 0.01);
    let out = section.advance_linear(&spec, 1.25, 0.08, 0.03).unwrap();
    let expected = ramp_response(&spec.motion_kernel, 0.01, 0.08, 1.25)
        + ramp_response(&spec.gust_kernel, 0.0, 0.03, 1.25);
    close(out.alpha_eff_rad, expected, 2e-15);
    close(out.cl, spec.cl_alpha * expected, 2e-14);
    assert_eq!(out.cm_quarter, spec.cm0);
    assert_eq!(out.f, 1.0);
}

#[test]
fn zero_duration_freezes_all_bits_while_explicit_jumps_have_instantaneous_gain() {
    for kernel in [WAGNER_JONES, KUSSNER_2POLE] {
        let mut state = DuhamelState::trim(0.0);
        state.x[1] = -0.0;
        let original = state;
        for zero in [0.0, -0.0] {
            state.advance_linear(&kernel, zero, 0.25).unwrap();
            assert_eq!(state.x.map(f64::to_bits), original.x.map(f64::to_bits));
            assert_eq!(state.input_prev.to_bits(), original.input_prev.to_bits());
        }
        state.apply_jump(&kernel, 0.25).unwrap();
        close(state.effective(), 0.25 * kernel.phi(0.0), 2e-17);
        for duration in [0.125, 1.0, 4.0] {
            let mut held = state;
            held.advance_linear(&kernel, duration, 0.25).unwrap();
            close(held.effective(), 0.25 * kernel.phi(duration), 1e-16);
        }
    }
    let spec = UNSTEADY_SECTION_V1;
    for f in [0.2, 1e-300] {
        let mut state = UnsteadySectionState::trim(&spec, 0.0);
        state.motion.x[1] = -0.0;
        state.separation.f = f;
        let original = state_bits(&state);
        for zero in [0.0, -0.0] {
            state.advance_linear(&spec, zero, 0.4, -0.2).unwrap();
            assert_eq!(state_bits(&state), original);
        }
        state.apply_jump(&spec, 0.4, -0.2).unwrap();
        assert_eq!(state.separation.f.to_bits(), f.to_bits());
        assert_eq!(state.motion.input_prev, 0.4);
        assert_eq!(state.gust.input_prev, -0.2);
    }
}

#[test]
fn tiny_intervals_and_extreme_admitted_poles_keep_their_ramp_limits() {
    for kernel in [WAGNER_JONES, KUSSNER_2POLE] {
        for duration in [f64::from_bits(1), 1e-14] {
            let mut state = DuhamelState::trim(0.0);
            state.advance_linear(&kernel, duration, 0.25).unwrap();
            for i in 0..2 {
                let z = kernel.b[i] * duration;
                // Analytic small-z limit of the integrated pole.
                let expected = kernel.a[i] * 0.25 * (1.0 - z / 2.0);
                close(state.x[i], expected, 3e-17);
            }
            assert!(state.effective().is_finite());
        }
    }
    let extreme = IndicialKernel {
        kernel_id: "test-finite-pole-product-overflow",
        a: [0.5, 0.5],
        b: [f64::MAX, f64::MAX],
    };
    extreme.admit().unwrap();
    let mut state = DuhamelState::trim(0.0);
    state.advance_linear(&extreme, 2.0, 1.0).unwrap();
    let expected = (0.5 / f64::MAX) / 2.0;
    assert!(expected > 0.0);
    close(state.x[0], expected, f64::from_bits(1));
    close(state.x[1], expected, f64::from_bits(1));
}

fn advance_kind(
    kind: usize,
    state: &mut UnsteadySectionState,
    spec: &UnsteadySectionSpec,
    motion: f64,
    gust: f64,
) -> Result<(), fs_airfoil::Refusal> {
    match kind {
        0 => state.advance(spec, 1.0, motion, gust),
        1 => state.advance_linear(spec, 1.0, motion, gust),
        _ => state.apply_jump(spec, motion, gust),
    }
    .map(|_| ())
}

#[test]
fn every_section_update_refuses_atomically_including_second_channel_and_output_failure() {
    let spec = UNSTEADY_SECTION_V1;
    let mut baseline = UnsteadySectionState::trim(&spec, 0.02);
    baseline.advance_linear(&spec, 0.25, 0.4, 0.05).unwrap();
    for kind in 0..3 {
        let mut state = baseline;
        let before = state_bits(&state);
        let error = advance_kind(kind, &mut state, &spec, 0.1, f64::NAN).unwrap_err();
        assert_eq!(error.code, "non-finite-input");
        assert_eq!(state_bits(&state), before);

        let mut invalid = spec;
        invalid.gust_kernel.b[0] = -1.0;
        assert_eq!(
            advance_kind(kind, &mut state, &invalid, 0.1, 0.2)
                .unwrap_err()
                .code,
            "kernel-params-invalid"
        );
        assert_eq!(state_bits(&state), before);

        invalid = spec;
        invalid.separation.t_f = 0.0;
        assert_eq!(
            advance_kind(kind, &mut state, &invalid, 0.1, 0.2)
                .unwrap_err()
                .code,
            "separation-model-invalid"
        );
        assert_eq!(state_bits(&state), before);

        invalid = spec;
        invalid.cm0 = f64::NAN;
        assert_eq!(
            advance_kind(kind, &mut state, &invalid, 0.1, 0.2)
                .unwrap_err()
                .code,
            "unsteady-section-spec-invalid"
        );
        assert_eq!(state_bits(&state), before);

        invalid = spec;
        invalid.cl_alpha = f64::MAX;
        assert_eq!(
            advance_kind(kind, &mut state, &invalid, 8.0, 0.0)
                .unwrap_err()
                .code,
            "unsteady-state-invalid"
        );
        assert_eq!(state_bits(&state), before);

        state.gust.x[1] = f64::INFINITY;
        let corrupt_before = state_bits(&state);
        assert_eq!(
            advance_kind(kind, &mut state, &spec, 0.1, 0.2)
                .unwrap_err()
                .code,
            "unsteady-state-invalid"
        );
        assert_eq!(state_bits(&state), corrupt_before);
    }
}

#[test]
fn ramp_clock_boundaries_and_unrepresentable_updates_refuse_without_mutation() {
    let spec = UNSTEADY_SECTION_V1;
    let initial = UnsteadySectionState::trim(&spec, 0.0);
    for invalid_ds in [
        -f64::from_bits(1),
        MAX_DS.next_up(),
        f64::INFINITY,
        f64::NAN,
    ] {
        let mut state = initial;
        let error = state
            .advance_linear(&spec, invalid_ds, 0.2, 0.1)
            .unwrap_err();
        assert_eq!(error.code, "reduced-time-increment-invalid");
        assert_eq!(state_bits(&state), state_bits(&initial));
    }
    let mut cap = initial;
    assert!(cap.advance_linear(&spec, MAX_DS, 0.2, 0.1).is_ok());

    let mut huge = DuhamelState::trim(-f64::MAX);
    let before = huge;
    assert_eq!(
        huge.advance_linear(&WAGNER_JONES, 1.0, f64::MAX)
            .unwrap_err()
            .code,
        "unsteady-state-invalid"
    );
    assert_eq!(huge, before);
    assert_eq!(
        huge.apply_jump(&WAGNER_JONES, f64::MAX).unwrap_err().code,
        "unsteady-state-invalid"
    );
    assert_eq!(huge, before);

    let mut invalid_separation = initial;
    invalid_separation.separation.f = 1.0_f64.next_up();
    let before = state_bits(&invalid_separation);
    assert_eq!(
        invalid_separation
            .advance_linear(&spec, 0.0, 0.2, 0.1)
            .unwrap_err()
            .code,
        "unsteady-state-invalid"
    );
    assert_eq!(state_bits(&invalid_separation), before);
}
