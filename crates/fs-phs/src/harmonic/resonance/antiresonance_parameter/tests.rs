//! Numerical range gates for the zero-safe implicit notch derivative.
use super::*;
use crate::harmonic::resonance::ResonanceFrequency;

fn minimum(value: C64, frequency: C64) -> DrivenAntiresonance {
    DrivenAntiresonance {
        minimum: ResonanceFrequency {
            angular_frequency_rad_s: 1.0,
            bracket_rad_s: (0.9, 1.1),
        },
        observation: value,
        observation_derivative: frequency,
        recovery: None,
        evaluations: 3,
        maximum_relative_residual: 0.0,
    }
}

#[test]
fn exact_zero_ignores_zero_weight_huge_jets_and_retains_tiny_frequency_motion() {
    // At H=0 and H_omega=1, implicit differentiation gives omega_p=-H_p
    // for these real jets. H_omegaomega and H_omega_p have exact zero weights,
    // regardless of magnitude. Scaling by either one would erase H_p.
    let mut queries = 0;
    let mut cancelled = || false;
    let result = resolve_parameter_response(
        minimum(C64::ZERO, C64::ONE),
        |omega, second_order, _| {
            queries += 1;
            assert_eq!(omega, 1.0);
            assert!(second_order);
            Ok(ParameterObservationJet {
                value: C64::ZERO,
                frequency: C64::ONE,
                parameter: C64::from_re(1e-200),
                second_frequency: C64::from_re(1e200),
                mixed: C64::from_re(1e200),
                relative_residual: 0.0,
            })
        },
        &mut cancelled,
    )
    .unwrap();
    assert_ne!(result.minimum_frequency_derivative_rad_s, 0.0);
    assert!(
        (result.minimum_frequency_derivative_rad_s / 1e-200 + 1.0).abs()
            <= 8.0 * f64::EPSILON
    );
    assert_eq!(result.observation_derivative, C64::ZERO);
    assert_eq!(result.scaled_power_curvature_s2, 2.0);
    assert_eq!(result.scaled_stationarity_residual_s, 0.0);
    assert!(result.recovery_derivative.is_none());
    assert_eq!(queries, 1);
    assert_eq!(result.evaluations, 4);
}

#[test]
fn nonzero_minimum_refuses_underflowed_normalized_products_instead_of_false_zero() {
    // H=i, H_omega=1 is stationary. Its true power curvature is
    // 2*(1+1e200), and omega_p=-1/(1+1e200). A common amplitude scale of
    // 1e200 leaves nonzero normalized factors 1e-200 whose product underflows.
    // Losing that product must refuse, not publish omega_p=0.
    let mut cancelled = || false;
    let result = resolve_parameter_response(
        minimum(C64::new(0.0, 1.0), C64::ONE),
        |_, _, _| {
            Ok(ParameterObservationJet {
                value: C64::new(0.0, 1.0),
                frequency: C64::ONE,
                parameter: C64::ONE,
                second_frequency: C64::new(0.0, 1e200),
                mixed: C64::from_re(1e200),
                relative_residual: 0.0,
            })
        },
        &mut cancelled,
    );
    assert!(matches!(result, Err(DrivenResonanceError::Invalid(_))));
}

#[test]
fn higher_order_response_zero_has_no_admitted_simple_notch_derivative() {
    // H=(omega-1)^2+p has H=H_omega=0, H_omegaomega=2 at p=0.
    // Positive curvature of H itself does not give positive curvature of |H|².
    let mut cancelled = || false;
    let result = resolve_parameter_response(
        minimum(C64::ZERO, C64::ZERO),
        |_, _, _| {
            Ok(ParameterObservationJet {
                value: C64::ZERO,
                frequency: C64::ZERO,
                parameter: C64::ONE,
                second_frequency: C64::from_re(2.0),
                mixed: C64::ZERO,
                relative_residual: 0.0,
            })
        },
        &mut cancelled,
    );
    assert!(matches!(result, Err(DrivenResonanceError::Invalid(_))));
}

#[test]
fn nonzero_minimum_refuses_subnormal_normalization_that_erases_cancellation_residue() {
    // All nonzero physical jet entries below are normal. Only normalization
    // by P creates subnormals: H_p.re/P and H_omega_p.im/P round to +t/-t,
    // erasing the representable physical numerator 2^-448. Refuse this range
    // instead of reporting zero parameter motion.
    let t = f64::from_bits(1_u64 << 14); // exactly 2^-1060, without powi underflow
    let p = 2.0_f64.powi(664);
    let hp = C64::new((p * t) * (1.0 + f64::EPSILON), p);
    let hwp = C64::new(0.0, -(p * t));
    let exact_numerator = hp.re + hwp.im;
    assert_eq!(exact_numerator, 2.0_f64.powi(-448));
    assert!(exact_numerator > 0.0);
    let mut cancelled = || false;
    let result = resolve_parameter_response(
        minimum(C64::new(0.0, 1.0), C64::ONE),
        |_, _, _| {
            Ok(ParameterObservationJet {
                value: C64::new(0.0, 1.0),
                frequency: C64::ONE,
                parameter: hp,
                second_frequency: C64::ZERO,
                mixed: hwp,
                relative_residual: 0.0,
            })
        },
        &mut cancelled,
    );
    assert!(matches!(result, Err(DrivenResonanceError::Invalid(_))));
}

#[test]
fn minimum_refuses_lost_low_product_word_despite_normal_high_words() {
    // All normalized components and product high words are normal, yet the
    // exact product tail lies at 2^-1104. Losing it before opposite terms
    // cancel erases a physical -2^-440 numerator and a shift near +2^60.
    let u = 2.0_f64.powi(-500);
    let p = 2.0_f64.powi(664);
    let h = C64::new(0.0, u);
    let hw = C64::from_re(u * (1.0 - f64::EPSILON));
    let hp = C64::new((p * u) * (1.0 + f64::EPSILON), p);
    let hwp = C64::new(0.0, -(p * u));
    // An independent fused physical contraction retains the exact remainder.
    let numerator = hp.re.mul_add(hw.re, h.im * hwp.im);
    assert_eq!(numerator, -2.0_f64.powi(-440));
    assert!(numerator.is_normal());
    let mut cancelled = || false;
    let result = resolve_parameter_response(
        minimum(h, hw),
        |_, _, _| {
            Ok(ParameterObservationJet {
                value: h,
                frequency: hw,
                parameter: hp,
                second_frequency: C64::new(0.0, 1.0),
                mixed: hwp,
                relative_residual: 0.0,
            })
        },
        &mut cancelled,
    );
    assert!(matches!(result, Err(DrivenResonanceError::Invalid(_))));
}
