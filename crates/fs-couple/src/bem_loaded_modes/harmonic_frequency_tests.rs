use super::super::super::tests::oscillator;
use super::*;
use fs_phs::RelaxationBranch;
use std::cell::Cell;

fn options() -> BemHarmonicFrequencyOptions {
    BemHarmonicFrequencyOptions {
        response: BemHarmonicOptions {
            max_load_cases: 8,
            max_response_bytes: 1_000_000,
            max_work: 1_000_000,
            power_relative_tolerance: 1e-10,
        },
        max_additional_response_bytes: 1_000_000,
        max_additional_work: 1_000_000,
    }
}

fn damped() -> BemLoadedModes {
    let mut model = oscillator(0.2);
    model.settings.max_dimension = 12;
    model
        .with_relaxation_branches(
            &[RelaxationBranch {
                projection: vec![0.5, -0.3],
                stiffness: 9_000.0,
                relaxation_time_s: 0.003,
            }],
            1,
            || false,
        )
        .unwrap()
}

fn solve(model: &BemLoadedModes, w: f64, force: &[C64]) -> BemHarmonicResponse {
    model
        .harmonic_response_batch(w, &[force], &options().response, || false)
        .unwrap()
        .remove(0)
}

fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(f64::MIN_POSITIVE),
        "{actual:?} != {expected:?}"
    );
}

fn scalar(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(f64::MIN_POSITIVE),
        "{actual:e} != {expected:e}"
    );
}

fn same_response(actual: &BemHarmonicResponse, expected: &BemHarmonicResponse) {
    for (a, b) in [
        (&actual.modal_displacement, &expected.modal_displacement),
        (&actual.modal_velocity, &expected.modal_velocity),
        (&actual.surface_pressure_pa, &expected.surface_pressure_pa),
        (
            &actual.surface_normal_velocity_m_s,
            &expected.surface_normal_velocity_m_s,
        ),
        (
            &actual.elastic_memory_sqrt_j,
            &expected.elastic_memory_sqrt_j,
        ),
    ] {
        assert_eq!(a.len(), b.len());
        for (a, b) in a.iter().zip(b) {
            assert_eq!(
                (a.re.to_bits(), a.im.to_bits()),
                (b.re.to_bits(), b.im.to_bits())
            );
        }
    }
    for (a, b) in [
        (
            actual.angular_frequency_rad_s,
            expected.angular_frequency_rad_s,
        ),
        (actual.mechanical_residual, expected.mechanical_residual),
        (actual.boundary_residual, expected.boundary_residual),
        (actual.memory_residual, expected.memory_residual),
        (actual.power.input_w, expected.power.input_w),
        (actual.power.kelvin_w, expected.power.kelvin_w),
        (actual.power.relaxation_w, expected.power.relaxation_w),
        (actual.power.radiated_w, expected.power.radiated_w),
        (actual.power.balance_w, expected.power.balance_w),
        (
            actual.power.apparent_power_scale_w,
            expected.power.apparent_power_scale_w,
        ),
        (
            actual.power.relative_balance,
            expected.power.relative_balance,
        ),
    ] {
        assert_eq!(a.to_bits(), b.to_bits());
    }
}

#[test]
fn g1_dark_maxwell_tangent_matches_exact_modulus_and_phased_source_work() {
    let mut base = oscillator(0.0);
    base.settings.max_dimension = 12;
    base.stiffness[3] = 60_000.0;
    base.viscosity[3] = 0.0;
    let tau = 1.0 / 300.0;
    let stiffness = 50_000.0_f64;
    let model = base
        .with_relaxation_branches(
            &[RelaxationBranch {
                projection: vec![0.0, 1.0],
                stiffness,
                relaxation_time_s: tau,
            }],
            1,
            || false,
        )
        .unwrap();
    let df = C64::new(0.0004, -0.0007);
    let mut saw_negative_heat_slope = false;
    for w in [10.0, 100.0, 295.8, 600.0, 1200.0] {
        let f = C64::new(0.7, -0.4) + df.scale(w - 300.0);
        let force = [C64::ZERO, f];
        let slope = [C64::ZERO, df];
        let r = model
            .harmonic_response_batch_with_frequency_derivative(
                w,
                &[&force],
                Some(&[&slope]),
                &options(),
                || false,
            )
            .unwrap()
            .remove(0);
        // Independent scalar Maxwell modulus; this coordinate has no skin
        // opening, Kelvin loss or mechanical coupling to the radiating mode.
        let denominator = C64::new(1.0, -w * tau);
        let g = C64::new(0.0, -w * tau) / denominator;
        let dg = C64::new(0.0, -tau) / (denominator * denominator);
        let dynamic = C64::from_re(60_000.0 - w * w) + g.scale(stiffness);
        let dynamic_slope = C64::from_re(-2.0 * w) + dg.scale(stiffness);
        let x = f / dynamic;
        let dx = (df - dynamic_slope * x) / dynamic;
        let v = C64::new(0.0, -w) * x;
        let dv = C64::new(0.0, -1.0) * (x + dx.scale(w));
        let memory = (g * x).scale(stiffness.sqrt());
        let dmemory = (dg * x + g * dx).scale(stiffness.sqrt());
        close(r.modal_displacement_derivative[1], dx, 2e-10);
        close(r.modal_velocity_derivative[1], dv, 2e-10);
        close(r.elastic_memory_derivative_sqrt_j[0], dmemory, 2e-10);
        let heat_slope = (memory.conj() * dmemory).re / tau;
        let input_slope = 0.5 * (df.conj() * v + f.conj() * dv).re;
        scalar(r.power_derivative.relaxation_w_per_rad_s, heat_slope, 2e-9);
        scalar(r.power_derivative.input_w_per_rad_s, input_slope, 2e-9);
        scalar(input_slope, heat_slope, 2e-10);
        saw_negative_heat_slope |= r.power_derivative.relaxation_w_per_rad_s < 0.0;
        assert!(r.power_derivative.kelvin_w_per_rad_s.abs() < 1e-25);
        assert!(r.power_derivative.radiated_w_per_rad_s.abs() < 1e-25);
        assert!(
            r.surface_pressure_derivative_pa
                .iter()
                .all(|p| p.abs() < 1e-20)
        );
        assert!(r.power_derivative.relative_balance < 1e-10);
        same_response(&r.response, &solve(&model, w, &force));
    }
    assert!(
        saw_negative_heat_slope,
        "frequency slopes must not be clipped as dissipation"
    );
}

fn difference_error(exact: &[C64], plus: &[C64], minus: &[C64], h: f64) -> f64 {
    assert_eq!(exact.len(), plus.len());
    assert_eq!(exact.len(), minus.len());
    let mut error = 0.0_f64;
    let mut scale = 0.0_f64;
    for ((&derivative, &a), &b) in exact.iter().zip(plus).zip(minus) {
        error = error.max(((a - b).scale(0.5 / h) - derivative).abs());
        scale = scale.max(derivative.abs());
    }
    error / scale.max(f64::MIN_POSITIVE)
}

#[test]
fn g3_full_coupled_state_and_power_tangents_match_refined_frequency_differences() {
    let model = damped();
    let source = [C64::new(1.0, 0.3), C64::new(-0.4, 0.7)];
    let slope = [C64::new(0.001, -0.002), C64::new(-0.0003, 0.0008)];
    for w in [280.0, 400.0, 540.0] {
        let r = model
            .harmonic_response_batch_with_frequency_derivative(
                w,
                &[&source],
                Some(&[&slope]),
                &options(),
                || false,
            )
            .unwrap()
            .remove(0);
        same_response(&r.response, &solve(&model, w, &source));
        assert!(r.mechanical_derivative_residual < 1e-10);
        assert!(r.boundary_derivative_residual < 1e-10);
        assert!(r.memory_derivative_residual < 1e-10);
        assert!(r.power_derivative.relative_balance < 1e-10);
        assert!(r.power_derivative.apparent_scale_w_per_rad_s > 0.0);
        assert!(
            r.power_derivative.balance_w_per_rad_s.abs()
                <= 1e-10 * r.power_derivative.apparent_scale_w_per_rad_s
        );
        let mut first_error = 0.0;
        for (level, h) in [0.1_f64, 0.025, 0.00625].into_iter().enumerate() {
            let plus_source = [source[0] + slope[0].scale(h), source[1] + slope[1].scale(h)];
            let minus_source = [source[0] - slope[0].scale(h), source[1] - slope[1].scale(h)];
            let plus = solve(&model, w + h, &plus_source);
            let minus = solve(&model, w - h, &minus_source);
            let mut error = 0.0_f64;
            for (exact, a, b) in [
                (
                    &r.modal_displacement_derivative,
                    &plus.modal_displacement,
                    &minus.modal_displacement,
                ),
                (
                    &r.modal_velocity_derivative,
                    &plus.modal_velocity,
                    &minus.modal_velocity,
                ),
                (
                    &r.surface_pressure_derivative_pa,
                    &plus.surface_pressure_pa,
                    &minus.surface_pressure_pa,
                ),
                (
                    &r.surface_normal_velocity_derivative_m_s,
                    &plus.surface_normal_velocity_m_s,
                    &minus.surface_normal_velocity_m_s,
                ),
                (
                    &r.elastic_memory_derivative_sqrt_j,
                    &plus.elastic_memory_sqrt_j,
                    &minus.elastic_memory_sqrt_j,
                ),
            ] {
                error = error.max(difference_error(exact, a, b, h));
            }
            let tolerance = 2e-4 * (h / 0.1).powi(2) + 3e-8;
            assert!(
                error < tolerance,
                "w={w}, h={h}, state tangent error={error:e}"
            );
            if level == 0 {
                first_error = error;
            }
            if level == 2 {
                assert!(error < first_error / 8.0 + 1e-7);
            }
            for (exact, a, b) in [
                (
                    r.power_derivative.input_w_per_rad_s,
                    plus.power.input_w,
                    minus.power.input_w,
                ),
                (
                    r.power_derivative.kelvin_w_per_rad_s,
                    plus.power.kelvin_w,
                    minus.power.kelvin_w,
                ),
                (
                    r.power_derivative.relaxation_w_per_rad_s,
                    plus.power.relaxation_w,
                    minus.power.relaxation_w,
                ),
                (
                    r.power_derivative.radiated_w_per_rad_s,
                    plus.power.radiated_w,
                    minus.power.radiated_w,
                ),
            ] {
                scalar((a - b) / (2.0 * h), exact, 4.0 * tolerance);
            }
        }
    }
}

#[test]
fn g5_one_assembly_batch_keeps_primal_bits_linearity_and_zero_drive_tangent() {
    let mut model = damped();
    model.settings.max_operator_evaluations = 1;
    let f = [C64::new(1.0, 0.3), C64::ZERO];
    let g = [C64::ZERO, C64::new(-0.4, 0.7)];
    let h = [f[0] + g[0], f[1] + g[1]];
    let df = [C64::new(0.001, -0.002), C64::ZERO];
    let dg = [C64::ZERO, C64::new(0.0003, 0.0008)];
    let dh = [df[0] + dg[0], df[1] + dg[1]];
    let zero = [C64::ZERO; 2];
    let forces: [&[C64]; 4] = [&f, &g, &h, &zero];
    let derivatives: [&[C64]; 4] = [&df, &dg, &dh, &f];
    let batch = model
        .harmonic_response_batch_with_frequency_derivative(
            350.0,
            &forces,
            Some(&derivatives),
            &options(),
            || false,
        )
        .unwrap();
    let primal = model
        .harmonic_response_batch(350.0, &forces, &options().response, || false)
        .unwrap();
    for (actual, expected) in batch.iter().zip(&primal) {
        same_response(&actual.response, expected);
    }
    for (a, b, sum) in [
        (
            &batch[0].modal_displacement_derivative,
            &batch[1].modal_displacement_derivative,
            &batch[2].modal_displacement_derivative,
        ),
        (
            &batch[0].surface_pressure_derivative_pa,
            &batch[1].surface_pressure_derivative_pa,
            &batch[2].surface_pressure_derivative_pa,
        ),
        (
            &batch[0].elastic_memory_derivative_sqrt_j,
            &batch[1].elastic_memory_derivative_sqrt_j,
            &batch[2].elastic_memory_derivative_sqrt_j,
        ),
    ] {
        for ((&a, &b), &sum) in a.iter().zip(b).zip(sum) {
            close(sum, a + b, 2e-10);
        }
    }
    // F=0 does not imply F'=0: its state tangent is the ordinary response to
    // F', while the quadratic source/heat/radiation power slopes remain zero.
    for (tangent, driven) in [
        (
            &batch[3].modal_displacement_derivative,
            &primal[0].modal_displacement,
        ),
        (
            &batch[3].modal_velocity_derivative,
            &primal[0].modal_velocity,
        ),
        (
            &batch[3].surface_pressure_derivative_pa,
            &primal[0].surface_pressure_pa,
        ),
        (
            &batch[3].elastic_memory_derivative_sqrt_j,
            &primal[0].elastic_memory_sqrt_j,
        ),
    ] {
        for (&a, &b) in tangent.iter().zip(driven) {
            close(a, b, 2e-10);
        }
    }
    let power = &batch[3].power_derivative;
    assert_eq!(
        [
            power.input_w_per_rad_s,
            power.kelvin_w_per_rad_s,
            power.relaxation_w_per_rad_s,
            power.radiated_w_per_rad_s
        ],
        [0.0; 4]
    );
    let fixed = model
        .harmonic_response_batch_with_frequency_derivative(350.0, &[&f], None, &options(), || false)
        .unwrap();
    let explicit = model
        .harmonic_response_batch_with_frequency_derivative(
            350.0,
            &[&f],
            Some(&[&zero]),
            &options(),
            || false,
        )
        .unwrap();
    assert_eq!(
        fixed[0].modal_displacement_derivative,
        explicit[0].modal_displacement_derivative
    );
    assert_eq!(
        fixed[0].surface_pressure_derivative_pa,
        explicit[0].surface_pressure_derivative_pa
    );
}

#[test]
fn g0_g4_bad_sources_limits_and_late_cancellation_leave_model_reusable() {
    let model = damped();
    let f = [C64::ONE, C64::ZERO];
    let df = [C64::new(0.002, -0.001), C64::ZERO];
    let malformed = [C64::new(f64::NAN, 0.0), C64::ZERO];
    for derivatives in [
        Some(&[][..]),
        Some(&[&df[..1]][..]),
        Some(&[&malformed[..]][..]),
    ] {
        assert!(
            model
                .harmonic_response_batch_with_frequency_derivative(
                    350.0,
                    &[&f],
                    derivatives,
                    &options(),
                    || false,
                )
                .is_err()
        );
    }
    assert!(
        model
            .harmonic_response_batch_with_frequency_derivative(
                350.0,
                &[&f[..1]],
                None,
                &options(),
                || false,
            )
            .is_err()
    );
    for defect in 0..6 {
        let mut limits = options();
        match defect {
            0 => limits.response.max_load_cases = 0,
            1 => limits.response.max_response_bytes = 0,
            2 => limits.response.max_work = 0,
            3 => limits.max_additional_response_bytes = 0,
            4 => limits.max_additional_work = 0,
            _ => limits.response.power_relative_tolerance = f64::NAN,
        }
        assert!(
            model
                .harmonic_response_batch_with_frequency_derivative(
                    350.0,
                    &[&f],
                    Some(&[&df]),
                    &limits,
                    || false,
                )
                .is_err()
        );
    }
    let polls = Cell::new(0);
    let before = model
        .harmonic_response_batch_with_frequency_derivative(
            350.0,
            &[&f, &f],
            Some(&[&df, &df]),
            &options(),
            || {
                polls.set(polls.get() + 1);
                false
            },
        )
        .unwrap();
    let total = polls.get();
    assert!(total > 4);
    for cut in [1, total / 2, total - 1, total] {
        let count = Cell::new(0);
        assert!(
            model
                .harmonic_response_batch_with_frequency_derivative(
                    350.0,
                    &[&f, &f],
                    Some(&[&df, &df]),
                    &options(),
                    || {
                        count.set(count.get() + 1);
                        count.get() >= cut
                    },
                )
                .is_err(),
            "cancellation at {cut}/{total}"
        );
    }
    let after = model
        .harmonic_response_batch_with_frequency_derivative(
            350.0,
            &[&f, &f],
            Some(&[&df, &df]),
            &options(),
            || false,
        )
        .unwrap();
    same_response(&before[0].response, &after[0].response);
    assert_eq!(
        before[0].modal_displacement_derivative,
        after[0].modal_displacement_derivative
    );
    assert_eq!(
        before[0].surface_pressure_derivative_pa,
        after[0].surface_pressure_derivative_pa
    );
    assert_eq!(
        before[0].elastic_memory_derivative_sqrt_j,
        after[0].elastic_memory_derivative_sqrt_j
    );
}
