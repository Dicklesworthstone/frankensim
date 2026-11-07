//! Integrate the actual intersection of a positive-friction disk and Hertz disk.
//!
//! Local-slip polar radius q has speed a*|omega|*q and coefficient
//! mu_zero*(1-q/k). Rays begin at the zero-slip point, distance d from the
//! patch centre. Their physical pressure roots give the radial bounds; circle
//! intersections split the angular bounds. No whole-patch node needs to find
//! a small crescent, and no large Coulomb response is subtracted.

use super::{
    Compensated, GAUSS, HertzFrictionQuadrature, HertzSlidingPatchParameters, HertzSlidingWrench,
    MAXIMUM_PANELS, PI, TriboError, det, finish,
};

#[allow(clippy::too_many_arguments)]
pub(super) fn intersecting_disk(
    parameters: HertzSlidingPatchParameters,
    slip: [f64; 2],
    spin: f64,
    speed: f64,
    edge_speed: f64,
    cutoff_speed: f64,
    quadrature: HertzFrictionQuadrature,
) -> Result<HertzSlidingWrench, TriboError> {
    let d = speed / edge_speed;
    let k = cutoff_speed / edge_speed;
    let inside = d < 1.0;
    let end = if inside { PI } else { det::asin(d.recip()) };
    let crossing = if k >= d + 1.0 {
        0.0
    } else {
        crossing_angle(d, k)?
    };
    // An external zero-slip point sees only rays subtended by the Hertz disk.
    // Before its tangent distance, the active circle terminates that view.
    let end = if !inside && k < det::sqrt((d - 1.0) * (d + 1.0)) {
        crossing
    } else {
        end
    };
    let half = (PI / 4.0) / quadrature.radial_panels as f64;
    let mut radial = [[0.0; 3]; 4 * MAXIMUM_PANELS];
    for panel in 0..quadrature.radial_panels {
        for (i, (node, weight)) in GAUSS.into_iter().enumerate() {
            let theta = (2.0 * panel as f64 + 1.0 + node) * half;
            let (s, c) = det::sin_cos(theta);
            radial[4 * panel + i] = [s * s, c * c, 2.0 * s * c * half * weight];
        }
    }
    let angular_half = (PI / 4.0) / quadrature.angular_panels as f64;
    let scale = (3.0 / PI) * parameters.normal_force_n * parameters.friction_coefficient;
    let mut integral = Compensated::default();
    for (start, stop) in [(0.0, crossing.min(end)), (crossing.min(end), end)] {
        if start == stop {
            continue;
        }
        for panel in 0..quadrature.angular_panels {
            for (node, weight) in GAUSS {
                let theta = (2.0 * panel as f64 + 1.0 + node) * angular_half;
                let (s, c) = det::sin_cos(theta);
                let phi = start + (stop - start) * s * s;
                let angular_weight = 2.0 * (stop - start) * s * c * angular_half * weight;
                let (sine, cosine) = det::sin_cos(phi);
                let transverse = d * sine;
                let root = if d == 1.0 {
                    // The near root is exactly zero on the rim. Preserve that
                    // geometry even when sqrt(1-sin^2(phi)) and cos(phi) round
                    // differently by much more than a microscopic cutoff k.
                    cosine.max(0.0)
                } else {
                    det::sqrt(((1.0 - transverse) * (1.0 + transverse)).max(0.0))
                };
                let middle = d * cosine;
                let (lower, width, pressure_span, cutoff_span) = if inside {
                    let height = (1.0 - d) * (1.0 + d);
                    let upper = if middle >= 0.0 {
                        root + middle
                    } else {
                        height / (root - middle)
                    };
                    let negative_root = if middle >= 0.0 {
                        height / (root + middle)
                    } else {
                        root - middle
                    };
                    (negative_root, upper.min(k), upper, k)
                } else {
                    // Recover the near root without subtracting two nearly
                    // equal values when the stationary point approaches the rim.
                    let near = (d - 1.0) * (d + 1.0) / (middle + root);
                    let available = if near <= 0.5 * k {
                        k - near
                    } else {
                        (k - middle) + root
                    };
                    (near, available.min(2.0 * root), 2.0 * root, available)
                };
                if width <= 0.0 {
                    continue;
                }
                for &[fraction, complement, radial_weight] in
                    &radial[..4 * quadrature.radial_panels]
                {
                    let offset = width * fraction;
                    let (q, moment_arm, pressure) = if inside {
                        (
                            offset,
                            offset - middle,
                            det::sqrt((pressure_span - offset).max(0.0))
                                * det::sqrt(offset + lower),
                        )
                    } else {
                        (
                            lower + offset,
                            offset - root,
                            det::sqrt(offset) * det::sqrt((pressure_span - offset).max(0.0)),
                        )
                    };
                    // Keep the small remaining coefficient and radial pressure
                    // factors directly, including on a narrow exterior crescent.
                    let coefficient = ((cutoff_span - width) + width * complement) / k;
                    let ring = if q >= 1.0 {
                        scale * (q * angular_weight)
                    } else {
                        (scale * q) * angular_weight
                    };
                    let load = (((ring * width) * radial_weight) * coefficient) * pressure;
                    integral.add([load * cosine, load * moment_arm, load * (edge_speed * q)]);
                }
            }
        }
    }
    if integral.total[2] == 0.0 {
        return Err(TriboError::InvalidInput {
            field: "hertz_friction_unresolved_active_region",
        });
    }
    finish(
        slip.map(|v| -integral.total[0] * (v / speed)),
        -(integral.total[1] * parameters.contact_radius_m) * spin.signum(),
        integral.total[2],
        slip,
        spin,
    )
}

fn crossing_angle(d: f64, k: f64) -> Result<f64, TriboError> {
    // Heron's sorted-side form resolves both a microscopic disk at d=1 and
    // a narrow angular view with d,k >> 1. acos alone loses the latter when
    // its cosine rounds to one; subtracting d +/- k loses the former.
    let mut sides = [d, k, 1.0];
    sides.sort_by(f64::total_cmp);
    let [c, b, a] = sides;
    let gap = c - (a - b);
    if gap <= 0.0 {
        return Err(TriboError::InvalidInput {
            field: "hertz_friction_unresolved_active_region",
        });
    }
    let sine = 0.5
        * ((det::sqrt(a + (b + c)) * det::sqrt(gap)) / d)
        * ((det::sqrt(c + (a - b)) * det::sqrt(a + (b - c))) / k);
    let cosine = if (d - 1.0).abs() <= (k - 1.0).abs() {
        ((d - 1.0) / k) * ((d + 1.0) / (2.0 * d)) + k / (2.0 * d)
    } else {
        ((k - 1.0) / d) * ((k + 1.0) / (2.0 * k)) + d / (2.0 * k)
    };
    let angle = det::atan2(sine, cosine);
    if angle.is_finite() && angle >= 0.0 {
        Ok(angle)
    } else {
        Err(TriboError::InvalidInput {
            field: "hertz_friction_active_geometry",
        })
    }
}
