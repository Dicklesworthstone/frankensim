//! Physical-midpoint history acceptance and cambered read-only response.
//! References integrate the continuous lag coordinates independently of the
//! production deficiency recurrence, with a changing chordwise clock.

use fs_airfoil::indicial::{IndicialKernel, MAX_DS};
use fs_airfoil::unsteady::{
    SectionMidpointInputs, UNSTEADY_SECTION_V1, UnsteadySectionSpec, UnsteadySectionState,
};

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}, tolerance {tolerance:.3e}"
    );
}

fn bits(state: &UnsteadySectionState) -> [u64; 7] {
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

fn inputs() -> SectionMidpointInputs {
    SectionMidpointInputs {
        motion_start_rad: 0.11,
        gust_start_rad: -0.02,
        motion_midpoint_rad: 0.15,
        gust_midpoint_rad: 0.01,
        effective_midpoint_rad: 0.42,
    }
}

#[test]
fn camber_shifts_loads_without_changing_geometric_history_or_separation() {
    let spec = UnsteadySectionSpec {
        cm0: 0.01 - core::f64::consts::PI * 0.02,
        ..UNSTEADY_SECTION_V1
    };
    let mut state = UnsteadySectionState::trim(&spec, 0.32);
    let accepted = state.advance_linear(&spec, 0.4, 0.36, 0.02).unwrap();
    let before = bits(&state);
    for zero in [0.0, -0.0] {
        let response = state.response(&spec, zero).unwrap();
        assert_eq!(response.cl.to_bits(), accepted.cl.to_bits());
        assert_eq!(response.cm_quarter.to_bits(), accepted.cm_quarter.to_bits());
        assert_eq!(
            response.alpha_eff_rad.to_bits(),
            accepted.alpha_eff_rad.to_bits()
        );
        assert_eq!(response.f.to_bits(), accepted.f.to_bits());
    }
    let cambered = state.response(&spec, -0.04).unwrap();
    let factor = ((1.0 + accepted.f.sqrt()) / 2.0).powi(2);
    let lift_shift = spec.cl_alpha * 0.04 * factor;
    let x_ac = 0.25 + spec.ac_shift * (1.0 - accepted.f);
    close(cambered.cl, accepted.cl + lift_shift, 1e-15);
    close(
        cambered.cm_quarter,
        accepted.cm_quarter + lift_shift * (0.25 - x_ac),
        1e-16,
    );
    assert_eq!(
        cambered.alpha_eff_rad.to_bits(),
        accepted.alpha_eff_rad.to_bits()
    );
    assert_eq!(cambered.f.to_bits(), accepted.f.to_bits());
    assert_eq!(bits(&state), before);
    assert!(state.response(&spec, f64::NAN).is_err());
    assert_eq!(bits(&state), before);
}

// z_i = a_i*u - x_i is continuous across input jumps. Solve
// dz_i/ds = b_i*(a_i*u-z_i) for a ramp from the ACTUAL start input;
// the incoming history may still be referenced to an extrapolated input.
fn continuous_lag_end(
    kernel: &IndicialKernel,
    pole: usize,
    z_start: f64,
    actual_start: f64,
    extrapolated_end: f64,
    ds: f64,
) -> f64 {
    let b = kernel.b[pole];
    let lost = -(-b * ds).exp_m1();
    let slope = (extrapolated_end - actual_start) / ds;
    z_start * (-b * ds).exp() + kernel.a[pole] * (actual_start * lost + slope * (ds - lost / b))
}

#[test]
fn midpoint_rebases_continuous_memory_and_uses_the_solved_midpoint_separation() {
    let spec = UNSTEADY_SECTION_V1;
    let mut state = UnsteadySectionState::trim(&spec, 0.08);
    state.apply_jump(&spec, 0.16, -0.03).unwrap();
    state.advance_linear(&spec, 0.3, 0.12, 0.05).unwrap();
    state.separation.f = 0.58;
    let before = state;
    let sample = inputs();
    let ds = 0.8;
    state.advance_midpoint(&spec, ds, sample).unwrap();

    for (old, new, kernel, start, midpoint) in [
        (
            before.motion,
            state.motion,
            spec.motion_kernel,
            sample.motion_start_rad,
            sample.motion_midpoint_rad,
        ),
        (
            before.gust,
            state.gust,
            spec.gust_kernel,
            sample.gust_start_rad,
            sample.gust_midpoint_rad,
        ),
    ] {
        let end = midpoint + (midpoint - start);
        assert_eq!(new.input_prev.to_bits(), end.to_bits());
        for pole in 0..2 {
            let z_start = kernel.a[pole] * old.input_prev - old.x[pole];
            let expected = continuous_lag_end(&kernel, pole, z_start, start, end, ds);
            let actual = kernel.a[pole] * new.input_prev - new.x[pole];
            close(actual, expected, 3e-16);
        }
    }
    let target = (-(sample.effective_midpoint_rad - spec.separation.alpha_break_rad)
        / spec.separation.decay_rad)
        .exp();
    let expected_f = target + (before.separation.f - target) * (-ds / spec.separation.t_f).exp();
    close(state.separation.f, expected_f, 2e-16);
    // The extrapolated response stays attached here; using it as the target
    // would relax in the wrong direction.
    assert!(state.response(&spec, 0.0).unwrap().alpha_eff_rad.abs() < 0.3);
    assert!(state.separation.f < before.separation.f);
}

#[test]
fn midpoint_zero_clock_freezes_bits_and_refusals_are_transactional() {
    let spec = UNSTEADY_SECTION_V1;
    let mut seed = UnsteadySectionState::trim(&spec, -0.0);
    seed.motion.x[0] = -0.0;
    seed.gust.input_prev = -0.0;
    seed.gust.x[1] = -0.0;
    seed.separation.f = 0.37;
    let before = bits(&seed);
    for ds in [0.0, -0.0] {
        let mut state = seed;
        state
            .advance_midpoint(
                &spec,
                ds,
                SectionMidpointInputs {
                    motion_start_rad: f64::MAX,
                    motion_midpoint_rad: -f64::MAX,
                    ..inputs()
                },
            )
            .unwrap();
        assert_eq!(bits(&state), before);
    }
    let bad_samples = [
        SectionMidpointInputs {
            motion_start_rad: f64::NAN,
            ..inputs()
        },
        SectionMidpointInputs {
            gust_start_rad: f64::INFINITY,
            ..inputs()
        },
        SectionMidpointInputs {
            motion_midpoint_rad: f64::NAN,
            ..inputs()
        },
        SectionMidpointInputs {
            gust_midpoint_rad: f64::MAX,
            ..inputs()
        },
        SectionMidpointInputs {
            effective_midpoint_rad: f64::NAN,
            ..inputs()
        },
    ];
    for sample in bad_samples {
        let mut state = seed;
        assert!(state.advance_midpoint(&spec, 0.4, sample).is_err());
        assert_eq!(bits(&state), before);
    }
    for ds in [-1.0, f64::NAN, f64::INFINITY, MAX_DS + 1.0] {
        let mut state = seed;
        assert!(state.advance_midpoint(&spec, ds, inputs()).is_err());
        assert_eq!(bits(&state), before);
    }
    let mut invalid = seed;
    invalid.separation.f = 1.01;
    let invalid_before = bits(&invalid);
    assert!(invalid.advance_midpoint(&spec, 0.0, inputs()).is_err());
    assert!(invalid.response(&spec, 0.0).is_err());
    assert_eq!(bits(&invalid), invalid_before);
    let mut bad_zero = seed;
    assert!(
        bad_zero
            .advance_midpoint(&spec, 0.0, bad_samples[0])
            .is_err()
    );
    assert_eq!(bits(&bad_zero), before);
}

#[test]
fn midpoint_accepts_history_without_claiming_a_physical_endpoint_load() {
    let spec = UnsteadySectionSpec {
        cl_alpha: f64::MAX,
        ..UNSTEADY_SECTION_V1
    };
    let mut state = UnsteadySectionState::trim(&spec, 0.0);
    state
        .advance_midpoint(
            &spec,
            0.3,
            SectionMidpointInputs {
                motion_start_rad: 0.0,
                gust_start_rad: 0.0,
                motion_midpoint_rad: 4.0,
                gust_midpoint_rad: 0.0,
                effective_midpoint_rad: 0.1,
            },
        )
        .unwrap();
    assert_eq!(state.motion.input_prev, 8.0);
    // Only an explicit load query can reject overflowing coefficients. The
    // extrapolated endpoint used to advance valid memory is not a force solve.
    assert!(state.response(&spec, 0.0).is_err());
}

fn reduced_time(t: f64) -> f64 {
    t + t * t / 4.0
}

fn variable_clock_circulation_error(steps: u32) -> f64 {
    let spec = UNSTEADY_SECTION_V1;
    let duration = 1.7;
    let dt = duration / f64::from(steps);
    let motion = |s: f64| 0.03 + 0.12 * s;
    let gust = |s: f64| -0.015 - 0.04 * s;
    let mut state = UnsteadySectionState::trim(&spec, motion(0.0));
    state.gust = fs_airfoil::unsteady::DuhamelState::trim(gust(0.0));
    for step in 0..steps {
        let t = f64::from(step) * dt;
        let t_mid = t + dt / 2.0;
        let start = reduced_time(t);
        let midpoint = reduced_time(t_mid);
        let mut predictor = state;
        predictor
            .apply_jump(&spec, motion(start), gust(start))
            .unwrap();
        predictor
            .advance_linear(
                &spec,
                0.5 * dt * (1.0 + 0.5 * t),
                motion(start),
                gust(start),
            )
            .unwrap();
        let midpoint_response = predictor
            .apply_jump(&spec, motion(midpoint), gust(midpoint))
            .unwrap();
        state
            .advance_midpoint(
                &spec,
                dt * (1.0 + 0.5 * t_mid),
                SectionMidpointInputs {
                    motion_start_rad: motion(start),
                    gust_start_rad: gust(start),
                    motion_midpoint_rad: motion(midpoint),
                    gust_midpoint_rad: gust(midpoint),
                    effective_midpoint_rad: midpoint_response.alpha_eff_rad,
                },
            )
            .unwrap();
    }
    let total_s = reduced_time(duration);
    state
        .apply_jump(&spec, motion(total_s), gust(total_s))
        .unwrap();
    let mut error: f64 = 0.0;
    // Both inputs are linear in exact reduced time, so each deficiency has a
    // closed form even though the physical chordwise speed changes in time.
    for (history, kernel, slope) in [
        (state.motion, spec.motion_kernel, 0.12),
        (state.gust, spec.gust_kernel, -0.04),
    ] {
        for pole in 0..2 {
            let expected =
                kernel.a[pole] * slope / kernel.b[pole] * -(-kernel.b[pole] * total_s).exp_m1();
            error = error.max((history.x[pole] - expected).abs());
        }
    }
    error
}

#[test]
fn variable_midpoint_clock_converges_quadratically_for_both_circulatory_channels() {
    let coarse = variable_clock_circulation_error(16);
    let medium = variable_clock_circulation_error(32);
    let fine = variable_clock_circulation_error(64);
    assert!(fine < 3e-6, "circulation error {fine:e}");
    for ratio in [coarse / medium, medium / fine] {
        assert!((3.6..4.4).contains(&ratio), "refinement ratio {ratio}");
    }
}

fn variable_clock_separation_error(steps: u32) -> f64 {
    let spec = UNSTEADY_SECTION_V1;
    let duration = 1.7;
    let dt = duration / f64::from(steps);
    let f_start = 0.35;
    let slope = 0.08;
    let mut state = UnsteadySectionState::trim(&spec, 0.0);
    state.separation.f = f_start;
    for step in 0..steps {
        let t_mid = (f64::from(step) + 0.5) * dt;
        let target = f_start + slope * reduced_time(t_mid);
        // Prescribe a smooth solved midpoint incidence whose static
        // separation target is linear in exact reduced time.
        let alpha = spec.separation.alpha_break_rad - spec.separation.decay_rad * target.ln();
        state
            .advance_midpoint(
                &spec,
                dt * (1.0 + 0.5 * t_mid),
                SectionMidpointInputs {
                    motion_start_rad: 0.0,
                    gust_start_rad: 0.0,
                    motion_midpoint_rad: 0.0,
                    gust_midpoint_rad: 0.0,
                    effective_midpoint_rad: alpha,
                },
            )
            .unwrap();
    }
    let total_s = reduced_time(duration);
    let expected = f_start
        + slope * (total_s - spec.separation.t_f * -(-total_s / spec.separation.t_f).exp_m1());
    (state.separation.f - expected).abs()
}

#[test]
fn actual_midpoint_separation_target_converges_quadratically_with_changing_speed() {
    let coarse = variable_clock_separation_error(16);
    let medium = variable_clock_separation_error(32);
    let fine = variable_clock_separation_error(64);
    assert!(fine < 5e-6, "separation error {fine:e}");
    for ratio in [coarse / medium, medium / fine] {
        assert!((3.6..4.4).contains(&ratio), "refinement ratio {ratio}");
    }
}
