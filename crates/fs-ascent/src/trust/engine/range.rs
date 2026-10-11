//! Range-safe fallbacks for the same Euclidean quadratic subproblem.
//! Ordinary well-conditioned reductions retain the established operation tree.
//! These are floating-point repairs, not interval certificates or a guarantee
//! for arbitrary conditioning, callback overflow, or changes of variable units.

use fs_math::dd::Dd;

// x = m * 2^e exactly; |m| in [1,2) for nonzero finite x. Splitting the
// subnormal mantissa directly avoids flushing it while finding its exponent.
fn split(x: f64) -> (f64, i32) {
    let bits = x.to_bits();
    let sign = bits & (1u64 << 63);
    let fraction = bits & ((1u64 << 52) - 1);
    let exponent = ((bits >> 52) & 0x7ff) as i32;
    if exponent != 0 {
        (f64::from_bits(sign | (1023u64 << 52) | fraction), exponent - 1023)
    } else if fraction == 0 {
        (x, 0)
    } else {
        let shift = fraction.leading_zeros() - 11;
        let fraction = (fraction << shift) & ((1u64 << 52) - 1);
        (f64::from_bits(sign | (1023u64 << 52) | fraction), -1022 - shift as i32)
    }
}

fn power(exponent: i32) -> f64 {
    debug_assert!((-1022..=1023).contains(&exponent));
    f64::from_bits(((exponent + 1023) as u64) << 52)
}

// m is bounded by 4. Split extreme powers so a representable subnormal or
// top-binade answer is not lost by materializing 2^e first.
fn scale(m: f64, exponent: i32) -> f64 {
    if exponent > 1023 {
        if exponent > 1025 { return m.signum() * f64::INFINITY; }
        (m * power(exponent - 1023)) * power(1023)
    } else if exponent < -1022 {
        if exponent < -1077 { return 0.0f64.copysign(m); }
        (m * power(exponent + 1022)) * power(-1022)
    } else {
        m * power(exponent)
    }
}

/// Sign of d'Hd and rr/(d'Hd). Finite callback vectors can have an overflowing
/// dot product while the Newton step is representable. Only exceptional dots
/// use product-wise binary exponents and compensated significand accumulation;
/// scaling whole vectors by their maxima would lose crossed large/small terms.
pub(super) fn curvature(rr: f64, d: &[f64], hd: &[f64]) -> (bool, f64) {
    let ordinary: f64 = d.iter().zip(hd).map(|(a, b)| a * b).sum();
    if ordinary.is_normal() {
        return (ordinary <= 0.0, rr / ordinary);
    }
    let largest = d.iter().zip(hd).filter(|(a, b)| **a != 0.0 && **b != 0.0)
        .map(|(&a, &b)| split(a).1 + split(b).1).max();
    let Some(largest) = largest else { return (true, 0.0); };
    let mut sum = Dd::ZERO;
    for (&a, &b) in d.iter().zip(hd) {
        if a == 0.0 || b == 0.0 { continue; }
        let (a, ea) = split(a);
        let (b, eb) = split(b);
        let product = Dd::from_f64(a) * Dd::from_f64(b);
        let shift = ea + eb - largest;
        sum = sum + Dd::from_pair(scale(product.hi, shift), scale(product.lo, shift));
    }
    if !Dd::ZERO.lt(sum) { return (true, 0.0); }
    let (numerator, ne) = split(rr);
    let (denominator, de) = split(sum.to_f64());
    (false, scale(numerator / denominator, ne - largest - de))
}

/// Follow the positive ray to the sphere, updating p only after all output
/// coordinates are finite. A dimensionless fallback returns the point directly:
/// its ray parameter may overflow even when every coordinate is representable.
pub(super) fn to_boundary(p: &mut [f64], d: &[f64], delta: f64) -> bool {
    if p.len() != d.len() || !delta.is_finite() || delta <= 0.0
        || p.iter().chain(d).any(|v| !v.is_finite()) { return false; }
    let pd: f64 = p.iter().zip(d).map(|(a, b)| a * b).sum();
    let dd: f64 = d.iter().map(|v| v * v).sum();
    let pp: f64 = p.iter().map(|v| v * v).sum();
    let disc = pd.mul_add(pd, dd * (delta * delta - pp));
    let root = fs_math::det::sqrt(disc.max(0.0));
    let tau = (-pd + root) / dd;
    if dd.is_normal() && disc.is_normal() && (delta * delta - pp) > 1e-8 * (delta * delta)
        && tau.is_finite() && tau > 0.0 && (pd <= 0.0 || root >= 2.0 * pd)
    {
        let next: Vec<_> = p.iter().zip(d).map(|(&a, &b)| tau.mul_add(b, a)).collect();
        if next.iter().all(|v| v.is_finite()) {
            p.copy_from_slice(&next);
            return true;
        }
    }
    let magnitude = d.iter().map(|v| v.abs()).fold(0.0f64, f64::max);
    if magnitude == 0.0 { return false; }
    let radius = Dd::from_f64(delta);
    let direction_scale = Dd::from_f64(magnitude);
    let mut aa = Dd::ZERO;
    let mut bb = Dd::ZERO;
    let mut norm_squared = Dd::ZERO;
    for (&point, &direction) in p.iter().zip(d) {
        let z = Dd::from_f64(point) / radius;
        let u = Dd::from_f64(direction) / direction_scale;
        aa = aa + u * u;
        bb = bb + z * u;
        norm_squared = norm_squared + z * z;
    }
    let gap = Dd::ONE - norm_squared;
    if !gap.hi.is_finite() || gap.lt(Dd::ZERO) { return false; }
    let root = (bb * bb + aa * gap).sqrt();
    // The rationalized positive root avoids subtracting nearly equal numbers
    // when p already points outwards near the boundary. The DD radial gap also
    // preserves a small transverse coordinate in a near-tangential crossing.
    let distance = if Dd::ZERO.lt(bb) { gap / (root + bb) }
        else { (root - bb) / aa };
    let next: Vec<_> = p.iter().zip(d).map(|(&point, &direction)| {
        let z = Dd::from_f64(point) / radius;
        let u = Dd::from_f64(direction) / direction_scale;
        ((z + distance * u) * radius).to_f64()
    }).collect();
    if !distance.hi.is_finite() || next.iter().any(|v| !v.is_finite()) { return false; }
    p.copy_from_slice(&next);
    true
}

#[cfg(test)]
#[path = "range_tests.rs"]
mod tests;
