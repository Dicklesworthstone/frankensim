//! All-finite-time proofs for declared, signed total mass flows.
//!
//! Signals clamp outside their own domains. The union of their knots therefore
//! covers every nonconstant piece, without adding a second scenario clock.
//! Exact opposite histories use IEEE sign symmetry. Other linear/held histories
//! use outward-rounded affine endpoint bounds and a conservative evaluation-
//! roundoff envelope. A failed bound is Unknown, never evidence of balance.

use crate::scenario::ValidationError;
use crate::{
    bc::{BcValue, BoundaryCondition},
    signal::{Interp, TimeSignal},
};
use fs_ivl::Interval as I;

#[derive(Clone, Copy)]
pub(crate) struct FluxTerm<'a> {
    pub bc: &'a BoundaryCondition,
    pub factor: f64,
}

fn finite(value: I) -> bool {
    value.lo().is_finite() && value.hi().is_finite()
}
fn magnitude(value: I) -> f64 {
    value.lo().abs().max(value.hi().abs())
}

impl FluxTerm<'_> {
    fn constant(self) -> Option<f64> {
        let value = match &self.bc.value {
            Some(BcValue::Uniform(q) | BcValue::Signal(TimeSignal::Constant(q))) => q.value,
            Some(BcValue::Signal(TimeSignal::Chebfun(p))) if p.cheb.coeffs().len() == 1 => {
                0.5 * p.cheb.coeffs()[0]
            }
            _ => return None,
        };
        let value = self.factor * value;
        value.is_finite().then_some(value)
    }

    fn bound(
        self,
        checkpoint: &mut impl FnMut(&'static str) -> Result<(), ValidationError>,
    ) -> Result<Option<I>, ValidationError> {
        if !self.factor.is_finite() {
            return Ok(None);
        }
        if let Some(value) = self.constant() {
            return Ok(Some(I::point(value.abs())));
        }
        let Some(BcValue::Signal(signal)) = &self.bc.value else {
            return Ok(None);
        };
        let mut bound = I::point(0.0);
        match signal {
            TimeSignal::Ramp {
                from,
                to,
                t_start,
                t_end,
            } => {
                if !from.value.is_finite()
                    || !to.value.is_finite()
                    || !t_start.is_finite()
                    || !t_end.is_finite()
                    || t_start >= t_end
                {
                    return Ok(None);
                }
                bound = I::point(from.value.abs()) + I::point(to.value.abs());
            }
            TimeSignal::Table { values, times, .. } => {
                if values.is_empty() || values.len() != times.len() {
                    return Ok(None);
                }
                for value in values {
                    checkpoint("net-flux all-time magnitude")?;
                    if !value.is_finite() {
                        return Ok(None);
                    }
                    bound = bound.hull(I::point(value.abs()));
                }
                // Each interpolated sample is the sum of two weighted values.
                bound = bound * I::point(2.0);
            }
            TimeSignal::Chebfun(profile) => {
                // For clamped finite physical coordinates, the stable affine
                // map lies in [-2,2], including endpoint-rounding slack. This
                // deliberately loose Clenshaw bound proves finite intermediate
                // arithmetic without assuming |T_k|<=1 for a rounded argument.
                let (mut b1, mut b2) = (I::point(0.0), I::point(0.0));
                for &coefficient in profile.cheb.coeffs().iter().skip(1).rev() {
                    checkpoint("net-flux all-time Clenshaw bound")?;
                    if !coefficient.is_finite() {
                        return Ok(None);
                    }
                    let b0 = I::point(4.0) * b1 + I::point(coefficient.abs()) + b2;
                    b2 = b1;
                    b1 = b0;
                    if !finite(b0) {
                        return Ok(None);
                    }
                }
                let Some(&c0) = profile.cheb.coeffs().first() else {
                    return Ok(None);
                };
                if !c0.is_finite() {
                    return Ok(None);
                }
                bound = I::point(2.0) * b1 + I::point(0.5) * I::point(c0.abs()) + b2;
            }
            TimeSignal::Constant(_) => return Ok(None),
        }
        let bound = bound * I::point(self.factor.abs());
        Ok(finite(bound).then_some(bound))
    }

    fn opposes(
        self,
        other: Self,
        checkpoint: &mut impl FnMut(&'static str) -> Result<(), ValidationError>,
    ) -> Result<bool, ValidationError> {
        if let (Some(a), Some(b)) = (self.constant(), other.constant()) {
            return Ok(a == -b);
        }
        if self.factor.abs() != other.factor.abs() {
            return Ok(false);
        }
        let negate = self.factor.is_sign_negative() == other.factor.is_sign_negative();
        let values_match = |a: f64, b: f64| if negate { a == -b } else { a == b };
        let (Some(BcValue::Signal(a)), Some(BcValue::Signal(b))) =
            (&self.bc.value, &other.bc.value)
        else {
            return Ok(false);
        };
        match (a, b) {
            (
                TimeSignal::Ramp {
                    t_start: a0,
                    t_end: a1,
                    from: af,
                    to: at,
                },
                TimeSignal::Ramp {
                    t_start: b0,
                    t_end: b1,
                    from: bf,
                    to: bt,
                },
            ) => Ok(a0 == b0
                && a1 == b1
                && values_match(af.value, bf.value)
                && values_match(at.value, bt.value)),
            (
                TimeSignal::Table {
                    times: at,
                    values: av,
                    interp: ai,
                    ..
                },
                TimeSignal::Table {
                    times: bt,
                    values: bv,
                    interp: bi,
                    ..
                },
            ) => {
                if ai != bi || at.len() != bt.len() || av.len() != bv.len() {
                    return Ok(false);
                }
                for (&a, &b) in at.iter().zip(bt) {
                    checkpoint("net-flux opposite knots")?;
                    if a != b {
                        return Ok(false);
                    }
                }
                for (&a, &b) in av.iter().zip(bv) {
                    checkpoint("net-flux opposite values")?;
                    if !values_match(a, b) {
                        return Ok(false);
                    }
                }
                Ok(true)
            }
            (TimeSignal::Chebfun(a), TimeSignal::Chebfun(b)) => {
                if a.cheb.domain() != b.cheb.domain()
                    || a.cheb.coeffs().len() != b.cheb.coeffs().len()
                {
                    return Ok(false);
                }
                for (&a, &b) in a.cheb.coeffs().iter().zip(b.cheb.coeffs()) {
                    checkpoint("net-flux opposite coefficients")?;
                    if !values_match(a, b) {
                        return Ok(false);
                    }
                }
                Ok(true)
            }
            _ => Ok(false),
        }
    }

    fn affine_piece(self, left: f64, right: f64) -> Option<(I, I, I)> {
        if let Some(value) = self.constant() {
            return Some((I::point(value), I::point(value), I::point(0.0)));
        }
        let Some(BcValue::Signal(signal)) = &self.bc.value else {
            return None;
        };
        let (start, end, from, to) = match signal {
            TimeSignal::Ramp {
                t_start,
                t_end,
                from,
                to,
            } => (*t_start, *t_end, from.value, to.value),
            TimeSignal::Table {
                times,
                values,
                interp,
                ..
            } => {
                if times.is_empty() || times.len() != values.len() {
                    return None;
                }
                let position = times.partition_point(|&t| t <= left);
                if position == 0 {
                    return self.constant_piece(values[0]);
                }
                if position == times.len() {
                    return self.constant_piece(*values.last()?);
                }
                if *interp == Interp::Hold {
                    return self.constant_piece(values[position - 1]);
                }
                (
                    times[position - 1],
                    times[position],
                    values[position - 1],
                    values[position],
                )
            }
            _ => return None,
        };
        if right <= start {
            return self.constant_piece(from);
        }
        if left >= end {
            return self.constant_piece(to);
        }
        if left < start || right > end || start >= end {
            return None;
        }
        let a = affine_endpoint(left, start, end, from, to)? * I::point(self.factor);
        let b = affine_endpoint(right, start, end, from, to)? * I::point(self.factor);
        // Fraction formation uses at most seven rounded operations, with the
        // overflow-safe time scaling only when opposite-sign width overflows.
        // For 0<=alpha<=1, absolute alpha error is <32 epsilon. Weighted-value
        // products/sum and final factor multiplication are covered by 64
        // epsilon and eight minimum subnormals BEFORE applying the factor.
        // The latter is essential when a large factor amplifies underflow.
        let error = ((I::point(from.abs()) + I::point(to.abs())) * I::point(64.0 * f64::EPSILON)
            + I::point(8.0 * f64::from_bits(1)))
            * I::point(self.factor.abs())
            + I::point(f64::from_bits(1));
        (finite(a) && finite(b) && finite(error)).then_some((a, b, error))
    }

    fn constant_piece(self, value: f64) -> Option<(I, I, I)> {
        let value = self.factor * value;
        value
            .is_finite()
            .then_some((I::point(value), I::point(value), I::point(0.0)))
    }
}

fn affine_endpoint(t: f64, start: f64, end: f64, from: f64, to: f64) -> Option<I> {
    if t == start {
        return Some(I::point(from));
    }
    if t == end {
        return Some(I::point(to));
    }
    let alpha = if (end - start).is_finite() {
        (I::point(t) - I::point(start)) / (I::point(end) - I::point(start))
    } else {
        let scale = I::point(start.abs().max(end.abs()));
        (I::point(t) / scale - I::point(start) / scale)
            / (I::point(end) / scale - I::point(start) / scale)
    };
    if !finite(alpha) {
        return None;
    }
    let value = (I::point(1.0) - alpha) * I::point(from) + alpha * I::point(to);
    finite(value).then_some(value)
}

/// An affirmative result proves the actual f64 flux sum meets the declared
/// relative tolerance for every finite f64 time, including clamped tails. The
/// declaration (knots, values, domains, factors and term order) is the witness;
/// replay reruns this same bounded arithmetic instead of trusting a verdict.
pub(crate) fn certifies_all_time(
    terms: &[FluxTerm<'_>],
    knots: &[f64],
    tolerance: f64,
    checkpoint: &mut impl FnMut(&'static str) -> Result<(), ValidationError>,
) -> Result<bool, ValidationError> {
    let mut gross_bound = I::point(0.0);
    for &term in terms {
        checkpoint("net-flux all-time bound")?;
        let Some(bound) = term.bound(checkpoint)? else {
            return Ok(false);
        };
        gross_bound = gross_bound + bound;
    }
    let nu = I::point(terms.len() as f64) * I::point(f64::EPSILON);
    if nu.hi() >= 0.25 {
        return Ok(false);
    }
    let gamma = nu / (I::point(1.0) - nu);
    if !finite(gross_bound * (I::point(1.0) + gamma)) {
        return Ok(false);
    }

    let mut paired = Vec::new();
    paired
        .try_reserve_exact(terms.len())
        .map_err(|_| ValidationError::AllocationRefused {
            resource: "net-flux opposite pairs",
            requested: terms.len(),
        })?;
    paired.resize(terms.len(), false);
    let mut exact_pairs = true;
    for (i, &term) in terms.iter().enumerate() {
        checkpoint("net-flux opposite pairing")?;
        if paired[i] {
            continue;
        }
        if term.opposes(term, checkpoint)? {
            paired[i] = true;
            continue;
        }
        let mut opposite = None;
        for j in i + 1..terms.len() {
            checkpoint("net-flux opposite candidate")?;
            if !paired[j] && term.opposes(terms[j], checkpoint)? {
                opposite = Some(j);
                break;
            }
        }
        if let Some(j) = opposite {
            paired[i] = true;
            paired[j] = true;
        } else {
            exact_pairs = false;
        }
    }
    // For finite f64 inputs whose exact sum is zero, sequential addition has
    // |fl(sum)|<=gamma_n*sum(abs). Subnormal sums are exact integer multiples
    // of the minimum subnormal. The computed gross sum is >=(1-gamma_n)*gross.
    // Leave a factor-two margin for the rounded tolerance*gross comparison.
    // If that product rounds to zero, gross is so small that every partial
    // sum is subnormal and addition is exact, so an opposite set sums to zero.
    let exact_sum_roundoff_fits =
        (I::point(2.0) * gamma / (I::point(1.0) - gamma)).hi() <= tolerance;
    if exact_pairs && exact_sum_roundoff_fits {
        return Ok(true);
    }

    for pair in knots.windows(2) {
        checkpoint("net-flux affine piece")?;
        let (left, right) = (pair[0], pair[1]);
        if left == right || left.next_up() >= right {
            continue;
        }
        let (mut net_left, mut net_right, mut gross_min, mut error) =
            (I::point(0.0), I::point(0.0), I::point(0.0), I::point(0.0));
        let mut all_zero = true;
        for (index, &term) in terms.iter().enumerate() {
            checkpoint("net-flux affine enclosure")?;
            // Exact pairs contribute zero to the true sum. Omitting their
            // absolute values only lowers the gross bound; the gamma term
            // still covers actual summation of every original provider.
            if paired[index] {
                continue;
            }
            let Some((a, b, roundoff)) = term.affine_piece(left, right) else {
                return Ok(false);
            };
            all_zero &= a == I::point(0.0) && b == a && roundoff == a;
            net_left = net_left + a;
            net_right = net_right + b;
            gross_min = gross_min + I::point(a.hull(b).abs().lo());
            error = error + roundoff;
        }
        if all_zero {
            if !exact_sum_roundoff_fits {
                return Ok(false);
            }
            continue;
        }
        let error = error + gamma * gross_bound;
        let upper_net = I::point(magnitude(net_left.hull(net_right))) + error;
        let lower_gross = gross_min - error;
        if !finite(upper_net)
            || lower_gross.lo() <= 0.0
            || upper_net.hi() > (I::point(tolerance) * lower_gross).lo()
        {
            return Ok(false);
        }
    }
    Ok(true)
}

#[cfg(test)]
mod tests {
    use super::{FluxTerm, certifies_all_time};
    use crate::scenario::ValidationError;
    use crate::{BcKind, BcValue, BoundaryCondition, ChebProfile, Compat, Physics, TimeSignal};
    use fs_cheb::Cheb1;
    use fs_qty::Dims;

    #[test]
    fn cancellation_inside_the_all_time_proof_returns_no_partial_certificate() {
        let conditions = [1.0, -1.0].map(|sign| BoundaryCondition {
            region: "port".into(),
            physics: Physics::IncompressibleFlow,
            kind: BcKind::MassFlowInlet,
            compatibility: Some(Compat::Incompressible),
            frame: 0,
            value: Some(BcValue::Signal(TimeSignal::Chebfun(ChebProfile {
                cheb: Cheb1::from_coeffs(0.0, 1.0, vec![sign, 2.0 * sign, -sign]),
                dims: Dims([0, 1, -1, 0, 0, 0]),
            }))),
        });
        let terms = conditions.each_ref().map(|bc| FluxTerm { bc, factor: 1.0 });
        let result = certifies_all_time(&terms, &[0.0, 1.0], 1e-9, &mut |phase| {
            if phase == "net-flux opposite coefficients" {
                Err(ValidationError::Cancelled {
                    phase,
                    completed: 0,
                    planned: 100,
                })
            } else {
                Ok(())
            }
        });
        assert!(matches!(
            result,
            Err(ValidationError::Cancelled {
                phase: "net-flux opposite coefficients",
                ..
            })
        ));
    }
}
