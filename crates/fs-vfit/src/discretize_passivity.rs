//! Global positive-real admission of the ACTUAL stored digital coefficients.
//!
//! With t = sin²(theta/2), each section's real part is a quadratic numerator
//! over its strictly positive squared denominator. Direct Bernstein products
//! retain DC/Nyquist endpoint structure without an expanded power polynomial.
//! Exact differentiators and integrators have separate boundary-pole storage
//! conditions at Nyquist and DC, respectively.

use super::{DigitalFilter, DiscretizeError, is_lossless_differentiator, is_lossless_integrator};
use crate::passivity::global::{
    MAX_GLOBAL_ORDER, bernstein_add, bernstein_mul, check_bernstein_nonnegative,
    interval_add as add, interval_mul as mul, value as bernstein_value,
};
use crate::passivity::{
    GlobalPassivityBudget, GlobalPassivityReport, GlobalPassivityUnresolved, GlobalPassivityVerdict,
};
use fs_ivl::Interval;
use fs_ivl::expansion::{expansion_sign, fast_expansion_sum_zeroelim};
use fs_math::det;

/// Certify impedance passivity over the entire digital frequency circle.
///
/// All stored binary64 coefficients are treated as exact real inputs. Ordinary
/// sections must satisfy strict Jury stability. Exact Tustin differentiators
/// and integrators are summed separately before checking their nonnegative
/// boundary-pole coefficients. Positive inertance cannot hide negative compliance.
/// Neither the direct slot nor an individual section must be nonnegative.
///
/// The report concerns the emitted coefficients, not the continuous fit or a
/// finite frequency grid. It does not prove roundoff-free runtime work. A
/// negative summed storage-element coefficient yields a violation without a
/// frequency witness: its failure is the boundary-pole storage condition.
/// Unestablished strict stability or unrepresentable arithmetic yields an
/// unresolved NumericRange report, never a pass. A frequency witness is omitted
/// if its rounded physical frequency cannot be rechecked with a negative
/// interval enclosure at the supplied sample period.
///
/// The budget's max_order bounds both the number of sections and the accumulated
/// ordinary polynomial degree (one for a first-order section, otherwise two),
/// with the shared hard cap of 64. Constant sections and exact storage elements need no
/// polynomial degree. The original filter and its state are never modified.
///
/// # Errors
/// Invalid sample interval or nonfinite stored coefficients.
pub fn check_digital_passivity(
    filter: &DigitalFilter,
    budget: GlobalPassivityBudget,
) -> Result<GlobalPassivityReport, DiscretizeError> {
    if !(filter.t_s > 0.0 && filter.t_s.is_finite()) {
        return Err(DiscretizeError::BadSampleInterval);
    }
    if !filter.direct.is_finite() {
        return Err(DiscretizeError::NonFiniteRuntimeValue {
            field: "digital_direct",
            index: None,
        });
    }
    let limit = budget.max_order.min(MAX_GLOBAL_ORDER);
    if filter.sections.len() > limit {
        return Ok(unresolved(GlobalPassivityUnresolved::OrderLimit));
    }
    let mut maximum = filter.direct.abs();
    for (index, section) in filter.sections.iter().enumerate() {
        if section.a.iter().chain(&section.b).any(|x| !x.is_finite()) {
            return Err(DiscretizeError::NonFiniteRuntimeValue {
                field: "digital_section",
                index: Some(index),
            });
        }
        for &value in &section.b {
            maximum = maximum.max(value.abs());
        }
    }

    // A common POWER-OF-TWO change of impedance units bounds the exact-input
    // sums. Round-trip admission below refuses any coefficient lost to
    // underflow; it never substitutes a rounded normalized source.
    let unit = if maximum == 0.0 {
        1.0
    } else {
        let exponent = maximum.to_bits() & 0x7ff0_0000_0000_0000;
        f64::from_bits(exponent.max(1_u64 << 52))
    };
    let Some(direct) = scaled_coefficient(filter.direct, unit) else {
        return Ok(unresolved(GlobalPassivityUnresolved::NumericRange));
    };
    let mut constants = vec![direct];
    let mut differentiators = Vec::new();
    let mut integrators = Vec::new();
    let mut factors = Vec::with_capacity(filter.sections.len());
    let mut degree = 0usize;

    for section in &filter.sections {
        let mut b = [0.0; 3];
        for (output, &input) in b.iter_mut().zip(&section.b) {
            let Some(scaled) = scaled_coefficient(input, unit) else {
                return Ok(unresolved(GlobalPassivityUnresolved::NumericRange));
            };
            *output = scaled;
        }
        if is_lossless_differentiator(section) {
            differentiators.push(b[0]);
            continue;
        }
        if is_lossless_integrator(section) {
            integrators.push(b[0]);
            continue;
        }
        let [a1, a2] = section.a;
        // These necessary bounds also keep the following exact sums finite.
        if a2.abs() >= 1.0 || a1.abs() >= 2.0 {
            return Ok(unresolved(GlobalPassivityUnresolved::NumericRange));
        }
        let plus = input_sum(&[1.0, a1, a2]);
        let minus = input_sum(&[1.0, -a1, a2]);
        if plus.lo() <= 0.0 || minus.lo() <= 0.0 {
            return Ok(unresolved(GlobalPassivityUnresolved::NumericRange));
        }
        if a1 == 0.0 && a2 == 0.0 && b[1] == 0.0 && b[2] == 0.0 {
            constants.push(b[0]);
            continue;
        }
        degree += if a2 == 0.0 && b[2] == 0.0 { 1 } else { 2 };
        if degree > limit {
            return Ok(unresolved(GlobalPassivityUnresolved::OrderLimit));
        }
        factors.push(section_factors(section.a, b, plus, minus));
    }

    // For a differentiator, u_n = w_(n+1)+w_n and p_n = g*(w_(n+1)-w_n).
    // Thus T*p_n*u_n = T*g*(w_(n+1)^2-w_n^2). Parallel copies have the
    // same zero-state recurrence; their SUMMED g must be nonnegative.
    // For an integrator, u_n=w_(n+1)-w_n and p_n=g*(w_(n+1)+w_n),
    // giving the same storage identity with an independent DC coordinate.
    // Expansion arithmetic retains even a weak negative remainder between
    // individually large positive and negative gains.
    if expansion_sign(&exact_sum(&differentiators)) < 0
        || expansion_sign(&exact_sum(&integrators)) < 0
    {
        return Ok(GlobalPassivityReport {
            verdict: GlobalPassivityVerdict::Violated,
            subdivisions: 0,
            violation_frequency_rad_s: None,
            unresolved: None,
        });
    }

    let mut numerator = vec![input_sum(&constants)];
    let mut denominator = vec![Interval::point(1.0)];
    for (part_numerator, part_denominator) in factors {
        numerator = bernstein_add(
            &bernstein_mul(&numerator, &part_denominator),
            &bernstein_mul(&denominator, &part_numerator),
        );
        denominator = bernstein_mul(&denominator, &part_denominator);
    }
    let certificate = check_bernstein_nonnegative(numerator.clone(), budget);
    let frequency = certificate
        .witness_parameter
        .and_then(|t| finite_witness(&numerator, filter.t_s, t));
    Ok(GlobalPassivityReport {
        verdict: certificate.verdict,
        subdivisions: certificate.subdivisions,
        violation_frequency_rad_s: frequency,
        unresolved: certificate.unresolved,
    })
}

// The report promises a witness at the RETURNED rounded physical frequency,
// not merely at the polynomial parameter that generated the candidate.
fn finite_witness(coefficients: &[Interval], sample_interval: f64, t: f64) -> Option<f64> {
    let angle = 2.0 * det::atan2(det::sqrt(t), det::sqrt(1.0 - t));
    let omega = angle / sample_interval;
    if !omega.is_finite() {
        return None;
    }
    let half_angle = mul(
        Interval::point(0.5),
        mul(Interval::point(omega), Interval::point(sample_interval)),
    );
    let sine = half_angle.sin();
    let parameter = mul(sine, sine).intersect(Interval::new(0.0, 1.0))?;
    (bernstein_value(coefficients, parameter).hi() < 0.0).then_some(omega)
}

fn unresolved(reason: GlobalPassivityUnresolved) -> GlobalPassivityReport {
    GlobalPassivityReport {
        verdict: GlobalPassivityVerdict::Unresolved,
        subdivisions: 0,
        violation_frequency_rad_s: None,
        unresolved: Some(reason),
    }
}

// This identity test is about the represented source, not a tolerance.
#[allow(clippy::float_cmp)]
fn scaled_coefficient(value: f64, unit: f64) -> Option<f64> {
    let scaled = value / unit;
    (scaled.is_finite() && scaled * unit == value).then_some(scaled)
}

// Inputs have already been bounded by normalization or strict Jury's necessary
// coefficient bounds. Their finite expansions cannot overflow. Empty is the
// exact zero, including cancellation between different represented inputs.
fn exact_sum(values: &[f64]) -> Vec<f64> {
    let mut result = Vec::new();
    for &value in values {
        if value != 0.0 {
            result = fast_expansion_sum_zeroelim(&result, &[value]);
        }
    }
    result
}

fn input_sum(values: &[f64]) -> Interval {
    let expansion = exact_sum(values);
    let mut sum = Interval::point(0.0);
    for value in expansion {
        sum = add(sum, Interval::point(value));
    }
    sum
}

fn section_factors(
    [a1, a2]: [f64; 2],
    [b0, b1, b2]: [f64; 3],
    plus: Interval,
    minus: Interval,
) -> (Vec<Interval>, Vec<Interval>) {
    let n_plus = input_sum(&[b0, b1, b2]);
    let n_minus = input_sum(&[b0, -b1, b2]);
    if a1 == 0.0 && a2 == 0.0 {
        // An FIR has denominator one. Its middle Bernstein coefficient is
        // b0-3*b2; form that exact input sum without a cancellation-prone
        // conversion from an already rounded power polynomial.
        let numerator = if b2 == 0.0 {
            vec![n_plus, n_minus]
        } else {
            vec![n_plus, input_sum(&[b0, -b2, -b2, -b2]), n_minus]
        };
        return (numerator, vec![Interval::point(1.0)]);
    }
    let numerator_ends = [mul(n_plus, plus), mul(n_minus, minus)];
    let denominator_ends = [mul(plus, plus), mul(minus, minus)];
    if a2 == 0.0 && b2 == 0.0 {
        return (numerator_ends.to_vec(), denominator_ends.to_vec());
    }

    let two = Interval::point(2.0);
    let one_minus_a2 = input_sum(&[1.0, -a2]);
    let b0_minus_b2 = input_sum(&[b0, -b2]);
    let denominator_middle = add(-mul(plus, minus), mul(two, mul(one_minus_a2, one_minus_a2)));
    let numerator_middle = add(
        -mul(
            Interval::point(0.5),
            add(mul(n_plus, minus), mul(n_minus, plus)),
        ),
        mul(two, mul(b0_minus_b2, one_minus_a2)),
    );
    (
        vec![numerator_ends[0], numerator_middle, numerator_ends[1]],
        vec![denominator_ends[0], denominator_middle, denominator_ends[1]],
    )
}
