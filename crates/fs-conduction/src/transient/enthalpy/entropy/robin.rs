//! Analytic entropy of affine Robin reservoir fields on one P1 triangle.
//!
//! The native heat rows remain q_i = integral N_i*h*(T-B) dA. This module
//! returns C = integral h*(T-B)*(1/B - 1/Bmax) dA, so the original rows give
//! bath entropy -sum(q_i)/Bmax-C and production
//! sum(q_i*(T_i-Bmax)/(T_i*Bmax))+C. Bmax is only an algebraic reference.
//! No mean bath, quadrature, diagonal face approximation or extra heat is used.
//!
//! Splitting along the median bath isotherm gives at most two triangles with
//! two equal bath values. Transverse integration is polynomial and exact.
//! The remaining rational cubic Bernstein moments have closed logarithmic
//! forms; a positive, geometrically bounded series avoids cancellation near
//! equal temperatures. Work and storage are bounded independently of contrast.

use super::{EnthalpyError, Sum, finite, poll};
use fs_exec::Cx;

const MOMENT: [f64; 4] = [0.25, 1.0 / 12.0, 1.0 / 12.0, 0.25];
const MAX_SERIES_TERMS: usize = 128;

/// Integral x^k(1-x)^(3-k)/(r+(1-r)x) dx minus its r=1 value.
/// Every returned remainder is nonnegative; no numerical sign is clamped.
fn reciprocal_excess(
    cx: &Cx<'_>,
    low: f64,
    high: f64,
) -> Result<[f64; 4], EnthalpyError> {
    let gap = (high - low) / high;
    if gap == 0.0 {
        return Ok([0.0; 4]);
    }
    if gap <= 0.5 {
        // 1/(1-gap*(1-x)) = sum_m gap^m*(1-x)^m.
        // The m=0 moment is removed analytically, preserving a tiny contrast.
        let mut terms = MOMENT;
        let mut sums = [0.0; 4];
        let mut lost = [0.0; 4];
        for m in 0..MAX_SERIES_TERMS {
            if m % 16 == 0 {
                poll(cx, m)?;
            }
            let mut converged = true;
            for k in 0..4 {
                terms[k] *= gap * (4 - k + m) as f64 / (5 + m) as f64;
                let next = sums[k] + terms[k];
                lost[k] += (sums[k] - next) + terms[k];
                sums[k] = next;
                // Every subsequent ratio is <gap, so this bounds the
                // omitted positive analytic tail, not floating-point error.
                converged &= terms[k] * gap / (1.0 - gap)
                    <= 0.25 * f64::EPSILON * (sums[k] + lost[k]);
            }
            if converged {
                return Ok(std::array::from_fn(|k| sums[k] + lost[k]));
            }
        }
        return Err(EnthalpyError::InvalidInput(
            "Robin reciprocal moments exceeded their fixed analytic-series work",
        ));
    }
    let ratio = low / high;
    let d = 1.0 - ratio;
    let d2 = d * d;
    let d3 = d2 * d;
    let d4 = d2 * d2;
    // Avoid a rounded or underflowed subnormal ratio inside log. Ordinary
    // ratios use one log, avoiding subtraction of two large logarithms.
    let logarithm = if ratio >= f64::MIN_POSITIVE {
        -ratio.ln()
    } else {
        high.ln() - low.ln()
    };
    let r2 = ratio * ratio;
    let moments = [
        (logarithm - d - 0.5 * d2 - d3 / 3.0) / d4,
        (d3 / 3.0 + 1.5 * ratio * d2 + r2 * d - ratio * logarithm) / d4,
        (d3 / 6.0 - 0.5 * ratio * d2 - r2 * d + r2 * logarithm) / d4,
        (d3 / 3.0 - 0.5 * ratio * d2 + r2 * d - r2 * ratio * logarithm) / d4,
    ];
    let excess = std::array::from_fn(|k| moments[k] - MOMENT[k]);
    if excess.iter().any(|value| !value.is_finite() || *value < 0.0) {
        return Err(EnthalpyError::InvalidInput(
            "Robin reciprocal moment is not representable",
        ));
    }
    Ok(excess)
}

/// One corner has bath value corner_bath; both opposite nodes have edge_bath.
#[allow(clippy::too_many_arguments)]
fn corner(
    cx: &Cx<'_>,
    area: f64,
    fraction: f64,
    htc: [f64; 3],
    delta: [f64; 3],
    corner_bath: f64,
    edge_bath: f64,
    reference: f64,
) -> Result<f64, EnthalpyError> {
    let high = corner_bath.max(edge_bath);
    let low = corner_bath.min(edge_bath);
    let extra = reciprocal_excess(cx, low, high)?;
    let shift = (reference - high) / reference;
    // With t at the exceptional corner, the slice Jacobian is 2A*(1-t).
    // We need t^n*(1-t)^(3-n), n=0,1,2. Reversing t maps n to 3-n.
    let weights: [f64; 3] = std::array::from_fn(|n| {
        let k = if corner_bath >= edge_bath { n } else { 3 - n };
        extra[k] + MOMENT[k] * shift
    });
    let block = [
        [weights[2], 0.5 * weights[1], 0.5 * weights[1]],
        [0.5 * weights[1], weights[0] / 3.0, weights[0] / 6.0],
        [0.5 * weights[1], weights[0] / 6.0, weights[0] / 3.0],
    ];
    let mut result = Sum::default();
    for i in 0..3 {
        for j in 0..3 {
            result.add(finite((htc[i] * block[i][j]) * (delta[j] / high))?)?;
        }
    }
    finite((2.0 * area * fraction) * result.total()?)
}

/// Exact rational face correction in W/K, evaluated in bounded f64 arithmetic.
/// The caller has admitted positive bath temperatures and the native fields.
pub(super) fn correction(
    cx: &Cx<'_>,
    area: f64,
    htc: [f64; 3],
    temperature: [f64; 3],
    bath: [f64; 3],
) -> Result<f64, EnthalpyError> {
    poll(cx, 0)?;
    let mut order = [0, 1, 2];
    for (i, j) in [(0, 1), (1, 2), (0, 1)] {
        if bath[order[i]] > bath[order[j]] {
            order.swap(i, j);
        }
    }
    let [a, b, c] = order;
    let (low, middle, high) = (bath[a], bath[b], bath[c]);
    if low == high {
        return Ok(0.0);
    }
    let delta = std::array::from_fn::<_, 3, _>(|i| temperature[i] - bath[i]);
    if middle == low {
        return corner(cx, area, 1.0, [htc[c], htc[a], htc[b]],
            [delta[c], delta[a], delta[b]], high, low, high);
    }
    if middle == high {
        return corner(cx, area, 1.0, [htc[a], htc[b], htc[c]],
            [delta[a], delta[b], delta[c]], low, high, high);
    }
    let low_fraction = (middle - low) / (high - low);
    let high_fraction = (high - middle) / (high - low);
    if low_fraction == 0.0 || high_fraction == 0.0 {
        return Err(EnthalpyError::InvalidInput(
            "Robin isotherm subdivision is not representable",
        ));
    }
    // Separate positive fractions retain a thin triangle even if subtracting
    // its fraction from one would round it away.
    let htc_at_cut = finite(high_fraction * htc[a] + low_fraction * htc[c])?;
    let delta_at_cut = finite(high_fraction * delta[a] + low_fraction * delta[c])?;
    let cold = corner(cx, area, low_fraction, [htc[a], htc[b], htc_at_cut],
        [delta[a], delta[b], delta_at_cut], low, middle, high)?;
    let hot = corner(cx, area, high_fraction, [htc[c], htc_at_cut, htc[b]],
        [delta[c], delta_at_cut, delta[b]], high, middle, high)?;
    finite(cold + hot)
}
