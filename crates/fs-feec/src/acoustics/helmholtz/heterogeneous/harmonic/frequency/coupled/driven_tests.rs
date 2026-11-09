//! G1 actual driven peaks and bandwidths of the original complete operator.
use super::*;
use fs_phs::harmonic::resonance::{
    ResonanceControls, ResonanceSearchError, ResonanceSurveyControls, search_resonance,
    survey_resonances,
};

fn scalar(value: f64) -> Csr {
    let mut matrix = Coo::new(1, 1);
    matrix.push(0, 0, value);
    matrix.assemble()
}
fn controls() -> ResonanceControls {
    ResonanceControls {
        absolute_frequency_tolerance_rad_s: 2e-9,
        relative_frequency_tolerance: 0.0,
        maximum_evaluations: 180,
        measure_half_power_bandwidth: true,
    }
}
fn sample(
    omega: f64,
    tangent: bool,
    k: f64,
    r: f64,
    arms: &[(f64, f64)],
    fluid: C64,
    traction: C64,
    coordinate: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(C64, C64, f64), PressureError> {
    let blocks = [scalar(0.0), scalar(1.0), scalar(0.0), scalar(0.0)];
    let owner = CoupledPencil::prepare(
        omega,
        blocks.each_ref(),
        &[0],
        vec![Wall {
            dofs: vec![0],
            mass: scalar(1.0),
            area_m2: 1.0,
            surface_density: 1.0,
            stiffness: k,
            resistance: r,
            arms: arms.to_vec(),
        }],
        MechanicalHarmonicLimits::default(),
        cancelled,
    )?;
    let source = owner.prepare_traction(0, &[1.0], cancelled)?;
    let (response, _) = owner.solve_with_tractions(&[fluid], &[(&source, traction)], cancelled)?;
    if tangent {
        let slope = owner.frequency_derivative_with_tractions(
            &response.pressure,
            &[fluid],
            &[C64::ZERO],
            &[],
            cancelled,
        )?;
        Ok((
            response.pressure[coordinate],
            slope.pressure[coordinate],
            response.backward_error.max(slope.backward_error),
        ))
    } else {
        Ok((
            response.pressure[coordinate],
            C64::ZERO,
            response.backward_error,
        ))
    }
}

#[test]
fn actual_pressure_displacement_velocity_peaks_have_their_own_exact_widths() {
    let k = 4.0_f64; // wall stiffness 3 PLUS the original fluid loading 1.
    let r = 0.4_f64;
    let expected_displacement = (k - 0.5 * r * r).sqrt();
    let half_delta = r * (k - 0.25 * r * r).sqrt();
    let lower = (expected_displacement * expected_displacement - half_delta).sqrt();
    let upper = (expected_displacement * expected_displacement + half_delta).sqrt();
    for coordinate in [0, 1, 2] {
        let result = search_resonance(
            (1.4, 2.5),
            controls(),
            || false,
            |omega, tangent, cancelled| {
                sample(
                    omega,
                    tangent,
                    3.0,
                    r,
                    &[],
                    C64::ZERO,
                    C64::new(0.7, -0.3),
                    coordinate,
                    cancelled,
                )
            },
        )
        .unwrap();
        let expected_peak = if coordinate == 2 {
            k.sqrt()
        } else {
            expected_displacement
        };
        let expected_width = if coordinate == 2 { r } else { upper - lower };
        assert!((result.peak.angular_frequency_rad_s - expected_peak).abs() < 3e-9);
        let bandwidth = result.half_power.unwrap();
        assert!((bandwidth.bandwidth_rad_s - expected_width).abs() < 4e-9);
        // The observed peak differs from the free pole's damped frequency.
        assert!((expected_peak - (k - 0.25 * r * r).sqrt()).abs() > 0.005);
        assert!(result.maximum_relative_residual < 5e-13);
    }
}

#[test]
fn hereditary_full_state_peak_and_coherent_drive_notch_are_not_conflated() {
    let target = 2.0_f64;
    let rate = 3.0_f64;
    let ka = 2.0_f64;
    let r = 0.2;
    let denominator = rate * rate + target * target;
    let imaginary = r * target + ka * rate * target / denominator;
    let imaginary_slope = r + ka * rate * (rate * rate - target * target) / denominator.powi(2);
    let real_slope = -2.0 * target + 2.0 * ka * rate * rate * target / denominator.powi(2);
    let real = -imaginary * imaginary_slope / real_slope;
    let k = real + target * target - ka * target * target / denominator - 1.0;
    let arms = [(ka.sqrt(), rate)];
    let peak = search_resonance(
        (1.0, 3.0),
        controls(),
        || false,
        |omega, tangent, cancelled| {
            sample(
                omega,
                tangent,
                k,
                r,
                &arms,
                C64::ZERO,
                C64::ONE,
                1,
                cancelled,
            )
        },
    )
    .unwrap();
    assert!((peak.peak.angular_frequency_rad_s - target).abs() < 3e-9);
    let width = peak.half_power.unwrap();
    for crossing in [width.lower, width.upper] {
        let value = sample(
            crossing.angular_frequency_rad_s,
            false,
            k,
            r,
            &arms,
            C64::ZERO,
            C64::ONE,
            1,
            &mut || false,
        )
        .unwrap()
        .0;
        assert!((value.norm_sq() - 0.5 * peak.observation.norm_sq()).abs() < 3e-9);
    }
    // At omega=2, fixed coherent b=-2i and traction=1 cancel in the ORIGINAL
    // displacement drive b/(i*omega)+traction; pressure itself need not vanish.
    let fluid = C64::new(0.0, -target);
    let notch = sample(target, true, k, r, &arms, fluid, C64::ONE, 1, &mut || false).unwrap();
    assert!(notch.0.abs() < 2e-13);
    assert!(notch.1.abs() > 0.1);
    let refused = search_resonance(
        (1.95, 2.05),
        controls(),
        || false,
        |omega, tangent, cancelled| {
            sample(omega, tangent, k, r, &arms, fluid, C64::ONE, 1, cancelled)
        },
    );
    assert!(matches!(
        refused,
        Err(ResonanceSearchError::NotBracketed(_))
    ));
    let survey = survey_resonances(
        &[1.95, 2.0, 2.05],
        ResonanceSurveyControls {
            absolute_frequency_tolerance_rad_s: 2e-9,
            relative_frequency_tolerance: 0.0,
            maximum_evaluations: 100,
            maximum_peaks: 4,
        },
        || false,
        |omega, tangent, cancelled| {
            sample(omega, tangent, k, r, &arms, fluid, C64::ONE, 1, cancelled)
        },
    )
    .unwrap();
    assert!(survey.peaks.is_empty());
}
