//! Exact affine-pulse means for a passive linear/two-photon cell cascade.
//!
//! Normalize the incident ramp by its interval maximum, x in [r, 1]. Every
//! prefix has the form A f_b(x), f_b=x/[b+(1-b)x], with A its peak output and
//! 0<=b<=1. A cell maps p to E p/(1+D p), hence b_next=b/(1+D A).
//! Its positive removal is L f_b + Q f_b f_b_next. Here L is the peak
//! linear-only removal and Q is total peak removal minus L. Cancellation in
//! Q cannot dominate total heat: both moments lie in [1/3,1], and L remains
//! positive. The original cell kernels supply all dimensional powers.
//!
//! These are exact frozen-coefficient time integrals up to binary64 arithmetic;
//! temperature evolution still needs the native temporal refinement owner.

use fs_math::det;

#[derive(Clone, Copy, Debug)]
pub(super) struct CellMean {
    pub absorbed_w: f64,
    pub transmitted_w: f64,
    pub next_shape: f64,
}

pub(super) fn mean_cell(
    minimum_fraction: f64,
    shape: f64,
    absorbed_peak_w: f64,
    transmitted_peak_w: f64,
    linear_absorbed_peak_w: f64,
    linear_transmitted_peak_w: f64,
) -> Result<CellMean, &'static str> {
    if ![minimum_fraction, shape]
        .iter()
        .all(|x| x.is_finite() && (0.0..=1.0).contains(x))
        || ![
            absorbed_peak_w,
            transmitted_peak_w,
            linear_absorbed_peak_w,
            linear_transmitted_peak_w,
        ]
        .iter()
        .all(|x| x.is_finite() && *x >= 0.0)
    {
        return Err("invalid affine-pulse cell partition");
    }
    let allowance = 64.0 * f64::EPSILON;
    if linear_absorbed_peak_w - absorbed_peak_w > allowance * linear_absorbed_peak_w {
        return Err("nonlinear peak absorption is below its linear partition");
    }
    let inverse = if linear_transmitted_peak_w == 0.0 {
        if transmitted_peak_w != 0.0 {
            return Err("nonlinear transmission exceeds an opaque linear cell");
        }
        0.0
    } else {
        let value = transmitted_peak_w / linear_transmitted_peak_w;
        if !value.is_finite() || value > 1.0 + allowance {
            return Err("nonlinear transmission exceeds the passive linear cell");
        }
        value.min(1.0)
    };
    let next_shape = shape * inverse;
    let first = first_moment(shape, minimum_fraction);
    let lower_first = 0.5 + 0.5 * minimum_fraction;
    let lower_second = (1.0 + minimum_fraction + minimum_fraction * minimum_fraction) / 3.0;
    let extra = (absorbed_peak_w - linear_absorbed_peak_w).max(0.0);
    let second = if extra == 0.0 {
        lower_second
    } else {
        product_moment(shape, next_shape, minimum_fraction)
    };
    let remaining = if shape == next_shape {
        first
    } else {
        first_moment(next_shape, minimum_fraction)
    };
    for (value, lower) in [
        (first, lower_first),
        (second, lower_second),
        (remaining, lower_first),
    ] {
        if !value.is_finite() || value < lower - allowance || value > 1.0 + allowance {
            return Err("unrepresentable affine-pulse moment");
        }
    }
    let absorbed_w = linear_absorbed_peak_w * first.clamp(lower_first, 1.0)
        + extra * second.clamp(lower_second, 1.0);
    let transmitted_w = transmitted_peak_w * remaining.clamp(lower_first, 1.0);
    if !absorbed_w.is_finite() || !transmitted_w.is_finite() {
        return Err("unrepresentable affine-pulse power");
    }
    Ok(CellMean {
        absorbed_w,
        transmitted_w,
        next_shape,
    })
}

fn first_moment(b: f64, r: f64) -> f64 {
    if b == 1.0 {
        return 0.5 + 0.5 * r;
    }
    1.0 - deficit(b, r)
}

fn deficit(b: f64, r: f64) -> f64 {
    if b == 0.0 || r == 1.0 {
        return 0.0;
    }
    let width = 1.0 - r;
    let q = (1.0 - b) * width;
    // Evaluate 1-q without losing a tiny positive b at r=0.
    let denominator = b + (1.0 - b) * r;
    b * width * first_remainder(q, denominator)
}

// Integral_0^1 y/(1-q*y) dy. The positive series resolves q=0.
fn first_remainder(q: f64, denominator: f64) -> f64 {
    if q == 0.0 {
        return 0.5;
    }
    if q < 0.25 {
        let mut power = 1.0;
        let mut sum = 0.0;
        for n in 0..32 {
            sum += power / f64::from(n + 2);
            power *= q;
        }
        sum
    } else {
        (-det::ln(denominator) - q) / (q * q)
    }
}

fn product_moment(a: f64, b: f64, r: f64) -> f64 {
    if a == 1.0 && b == 1.0 {
        return (1.0 + r + r * r) / 3.0;
    }
    if a == 0.0 {
        return first_moment(b, r);
    }
    if b == 0.0 {
        return first_moment(a, r);
    }
    if r == 1.0 {
        return 1.0;
    }
    let da = deficit(a, r);
    let db = deficit(b, r);
    // The omitted positive cross term is <=min(da,db), below a fraction of
    // one ulp of this moment (which is >=1/3). This also avoids 0*infinity
    // when an opaque nonlinear prefix leaves a subnormal shape parameter.
    if da.min(db) < f64::EPSILON / 32.0 {
        return 1.0 - da - db;
    }
    let width = 1.0 - r;
    let smaller = a.min(b);
    let larger = a.max(b);
    let hi = (1.0 - smaller) * width;
    let lo = (1.0 - larger) * width;
    let denominator_hi = smaller + (1.0 - smaller) * r;
    let denominator_lo = larger + (1.0 - larger) * r;
    // Integral y²/[(1-hi*y)(1-lo*y)]. This form never subtracts nearby
    // first_remainder values. The logarithmic divided difference has an
    // explicit equal-argument limit, including nearly identical layers.
    let cross_integral = if hi < 0.25 {
        let mut homogeneous = 1.0;
        let mut power_hi = 1.0;
        let mut sum = 0.0;
        for n in 0..32 {
            sum += homogeneous / f64::from(n + 3);
            power_hi *= hi;
            homogeneous = lo * homogeneous + power_hi;
        }
        sum
    } else {
        let delta = denominator_lo - denominator_hi;
        let logarithmic_mean_inverse = if delta == 0.0 {
            1.0 / denominator_hi
        } else {
            log1p_ratio(delta / denominator_hi) / denominator_hi
        };
        (logarithmic_mean_inverse - 1.0 - (hi + lo) * first_remainder(lo, denominator_lo))
            / (hi * hi)
    };
    1.0 - da - db + (a * width) * (b * width) * cross_integral
}

// ln(1+x)/x with its removable x=0 limit. The atanh series avoids the
// rounded 1+x that otherwise destroys almost equal denominator queries.
fn log1p_ratio(x: f64) -> f64 {
    if x < 1.0 {
        let z = x / (2.0 + x);
        let square = z * z;
        let mut power = 1.0;
        let mut sum = 0.0;
        for n in 0..20 {
            sum += power / f64::from(2 * n + 1);
            power *= square;
        }
        2.0 * sum / (2.0 + x)
    } else {
        det::ln(1.0 + x) / x
    }
}

#[cfg(test)]
#[path = "moments/tests.rs"]
mod tests;
