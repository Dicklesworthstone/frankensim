//! Original orthotropic Coulomb traction over an actual elliptical Hertz patch.
//!
//! The footprint and material texture have independent axes. Local slip is
//! `u + omega J r`; the original kinetic ellipse supplies its associated support
//! point, not a scalar coefficient multiplying the slip direction. The same
//! pressure integral determines force, centre moment and positive physical heat.
//!
//! In normalized disk coordinates, unoriented lines through the zero-slip point
//! make the local traction constant on each signed half-line. The two required
//! radial Hertz moments integrate exactly. Only one bounded angular quadrature
//! remains. Outside the disk, `p=sin(theta)` removes the rim square root, while
//! symmetric directions retain the small force difference responsible for torque
//! at very small representable spin. A rank-one kinetic ellipse splits at its direction
//! discontinuity instead of relying on a quadrature node to find it.
//!
//! This is prescribed instantaneous gross sliding. Rest selects zero traction;
//! it is not a static reaction solve. No shear history, normal-contact evolution,
//! velocity evolution, adhesion, rolling-deformation loss, heat partition or
//! measured interface calibration is implied. Angular refinement is a numerical
//! convergence requirement, not an error certificate.

use core::f64::consts::FRAC_PI_2;

use fs_math::{dd::Dd, det};

use super::{HertzSlidingPatch, HertzSlidingPatchParameters, HertzSlidingWrench, finite};
use crate::tangential::elliptic::{EllipticCoulomb, from_principal, support, to_principal};
use crate::{TangentialLaw, TriboError, nonnegative_finite, positive_finite};

const GAUSS: [(f64, f64); 4] = [
    (-0.861_136_311_594_052_6, 0.347_854_845_137_453_85),
    (-0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
    (0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
    (0.861_136_311_594_052_6, 0.347_854_845_137_453_85),
];

/// Original normal-owner geometry in the current tangent coordinates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EllipticHertzPatchParameters {
    pub normal_force_n: f64,
    /// Positive semi-axis lengths [m], in their original supplied order.
    pub semi_axes_m: [f64; 2],
    /// Direction of the first footprint axis; normalized at construction.
    pub footprint_axis: [f64; 2],
}

/// Bounded resolution of the remaining one-dimensional angular integral.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EllipticHertzQuadrature {
    /// 1..=256 panels per smooth half-angle interval, four Gauss nodes per panel.
    /// Rank-one friction may split that interval once: at most 4096 original
    /// directional support evaluations, with no allocation or adaptive growth.
    pub angular_panels: usize,
}

/// Original constitutive refusal or cancellation, with no partial response.
#[derive(Clone, Debug, PartialEq)]
pub enum EllipticHertzError {
    Law(TriboError),
    Cancelled,
}
impl From<TriboError> for EllipticHertzError {
    fn from(error: TriboError) -> Self {
        Self::Law(error)
    }
}
impl core::fmt::Display for EllipticHertzError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Law(error) => error.fmt(f),
            Self::Cancelled => write!(f, "elliptical Hertz sliding cancelled"),
        }
    }
}
impl std::error::Error for EllipticHertzError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Law(error) => Some(error),
            Self::Cancelled => None,
        }
    }
}

/// Immutable elliptical pressure field and unchanged directional interface law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EllipticHertzSlidingPatch {
    parameters: EllipticHertzPatchParameters,
    law: EllipticCoulomb,
    texture: [f64; 2],
    capacities: [f64; 2],
    circular_isotropic: Option<HertzSlidingPatch>,
}

impl EllipticHertzSlidingPatch {
    pub fn new(
        mut parameters: EllipticHertzPatchParameters,
        law: EllipticCoulomb,
    ) -> Result<Self, TriboError> {
        nonnegative_finite(parameters.normal_force_n, "elliptic_hertz_normal_force_n")?;
        for a in parameters.semi_axes_m {
            positive_finite(a, "elliptic_hertz_semi_axis_m")?;
        }
        for x in parameters.footprint_axis {
            finite(x, "elliptic_hertz_footprint_axis")?;
        }
        let norm = det::hypot(parameters.footprint_axis[0], parameters.footprint_axis[1]);
        positive_finite(norm, "elliptic_hertz_footprint_axis_norm")?;
        parameters.footprint_axis = parameters.footprint_axis.map(|x| x / norm);
        let mu = law.kinetic_coefficients();
        let capacities = mu.map(|x| x * parameters.normal_force_n);
        for x in capacities {
            finite(x, "elliptic_hertz_force_capacity_n")?;
        }
        finite(
            capacities[0].max(capacities[1])
                * parameters.semi_axes_m[0].max(parameters.semi_axes_m[1]),
            "elliptic_hertz_moment_bound_nm",
        )?;
        let circular_isotropic =
            if parameters.semi_axes_m[0] == parameters.semi_axes_m[1] && mu[0] == mu[1] {
                Some(HertzSlidingPatch::new(HertzSlidingPatchParameters {
                    normal_force_n: parameters.normal_force_n,
                    contact_radius_m: parameters.semi_axes_m[0],
                    friction_coefficient: mu[0],
                })?)
            } else {
                None
            };
        Ok(Self {
            parameters,
            law,
            texture: law.initial_state().principal_axis(),
            capacities,
            circular_isotropic,
        })
    }

    /// Admitted geometry, with the footprint direction normalized once.
    #[must_use]
    pub const fn parameters(&self) -> EllipticHertzPatchParameters {
        self.parameters
    }

    /// The original ordered static/kinetic data and supplied material texture.
    #[must_use]
    pub const fn friction_law(&self) -> &EllipticCoulomb {
        &self.law
    }

    pub fn evaluate(
        &self,
        slip_velocity_m_s: [f64; 2],
        spin_rad_s: f64,
        quadrature: EllipticHertzQuadrature,
    ) -> Result<HertzSlidingWrench, EllipticHertzError> {
        self.evaluate_cancellable(slip_velocity_m_s, spin_rad_s, quadrature, || false)
    }

    /// Integrate the original local kinetic support under prescribed rigid
    /// translation and spin. Positive scaling of both velocities leaves the
    /// wrench unchanged and scales physical power. Cancellation is polled at
    /// every angular node, including before and after every exact fast path.
    pub fn evaluate_cancellable(
        &self,
        slip_velocity_m_s: [f64; 2],
        spin_rad_s: f64,
        quadrature: EllipticHertzQuadrature,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HertzSlidingWrench, EllipticHertzError> {
        if !(1..=256).contains(&quadrature.angular_panels) {
            return Err(invalid("elliptic_hertz_angular_panels").into());
        }
        poll(&mut cancelled)?;
        for x in slip_velocity_m_s {
            finite(x, "slip_velocity_m_s")?;
        }
        finite(spin_rad_s, "spin_rad_s")?;
        let at_rest = slip_velocity_m_s == [0.0; 2] && spin_rad_s == 0.0;
        let response = if let Some(circle) = self.circular_isotropic {
            circle.evaluate(slip_velocity_m_s, spin_rad_s)?
        } else if at_rest || self.capacities == [0.0; 2] {
            finish([Dd::ZERO; 4], slip_velocity_m_s, spin_rad_s, at_rest)?
        } else if spin_rad_s == 0.0 {
            let velocity = to_principal(self.texture, slip_velocity_m_s)?;
            let resistive = support(self.capacities, velocity)?;
            let force = from_principal(self.texture, resistive.map(|x| -x))?;
            let power = dot(resistive, velocity);
            finish(
                [
                    Dd::from_f64(force[0]),
                    Dd::from_f64(force[1]),
                    Dd::ZERO,
                    power,
                ],
                slip_velocity_m_s,
                spin_rad_s,
                false,
            )?
        } else {
            self.integrate(slip_velocity_m_s, spin_rad_s, quadrature, &mut cancelled)?
        };
        poll(&mut cancelled)?;
        Ok(response)
    }

    fn map(&self, z: [f64; 2]) -> Result<[f64; 2], TriboError> {
        from_principal(
            self.parameters.footprint_axis,
            [
                self.parameters.semi_axes_m[0] * z[0],
                self.parameters.semi_axes_m[1] * z[1],
            ],
        )
    }

    fn integrate(
        &self,
        slip: [f64; 2],
        spin: f64,
        controls: EllipticHertzQuadrature,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<HertzSlidingWrench, EllipticHertzError> {
        // q = A^-1 J u. Neither d=qnorm/|omega| outside the disk nor its square
        // is formed, so weak spin cannot create a huge displaced origin.
        let mut q = to_principal(self.parameters.footprint_axis, turn(slip))?;
        for (i, x) in q.iter_mut().enumerate() {
            *x = finite(
                *x / self.parameters.semi_axes_m[i],
                "elliptic_hertz_inverse_slip",
            )?;
        }
        let qnorm = finite(det::hypot(q[0], q[1]), "elliptic_hertz_inverse_speed")?;
        let inside = qnorm <= spin.abs();
        let ratio = if inside {
            qnorm / spin.abs()
        } else {
            spin.abs() / qnorm
        };
        if qnorm > 0.0 && ratio == 0.0 {
            return Err(invalid("elliptic_hertz_unrepresentable_motion_ratio").into());
        }
        let e1 = if qnorm == 0.0 {
            [1.0, 0.0]
        } else {
            q.map(|x| spin.signum() * (x / qnorm))
        };
        let e2 = turn(e1);
        let a = self.map(e1)?;
        let b = self.map(e2)?;
        let va = turn(a).map(|x| -spin.signum() * x);
        let vb = turn(b).map(|x| spin.signum() * x);
        let split = self.rank_one_split(va, vb, inside, ratio)?;
        let boundaries = [0.0, split.unwrap_or(FRAC_PI_2), FRAC_PI_2];
        let mut sum = [Dd::ZERO; 4];
        for interval in boundaries.windows(2) {
            if interval[0] == interval[1] {
                continue;
            }
            let width = (interval[1] - interval[0]) / controls.angular_panels as f64;
            for panel in 0..controls.angular_panels {
                let centre = interval[0] + (panel as f64 + 0.5) * width;
                for (node, gauss_weight) in GAUSS {
                    poll(cancelled)?;
                    let theta = centre + 0.5 * width * node;
                    let (s, c) = (det::sin(theta), det::cos(theta));
                    let (base, delta, weight, first, second, speed) = if inside {
                        let d = ratio;
                        let r2 = (1.0 - d) * (1.0 + d) + (d * c) * (d * c);
                        (
                            va.map(|x| c * x),
                            vb.map(|x| s * x),
                            1.5 * r2,
                            c * (5.0 * (d * s) * (d * s) - 1.0) / 4.0,
                            s * ((d * c) * (d * c) + r2 / 4.0),
                            spin.abs() * ((d * c) * (d * c) + r2 / 4.0),
                        )
                    } else {
                        let k = ratio;
                        let radial_cos = if k == 1.0 {
                            c
                        } else {
                            det::sqrt((1.0 - k) * (1.0 + k) + (k * c) * (k * c))
                        };
                        let extra = k * (c * c) / (4.0 * radial_cos);
                        (
                            va.map(|x| radial_cos * x),
                            vb.map(|x| (k * s) * x),
                            1.5 * c * c * c,
                            k * (5.0 * s * s - 1.0) / 4.0,
                            s * (radial_cos + k * extra),
                            qnorm * radial_cos + spin.abs() * extra,
                        )
                    };
                    let (mean, difference) =
                        paired_support(self.capacities, self.texture, base, delta)?;
                    let force_weight = if inside { weight * ratio * c } else { weight };
                    let values = [
                        force_weight * mean[0],
                        force_weight * mean[1],
                        weight * (first * cross(a, mean) + second * cross(b, difference)),
                        -weight * speed * (dot(mean, base) + dot(difference, delta)).to_f64(),
                    ];
                    let quadrature_weight = Dd::from_f64(0.5 * width * gauss_weight);
                    for (total, value) in sum.iter_mut().zip(values) {
                        finite(value, "elliptic_hertz_integrand")?;
                        *total = *total + quadrature_weight * Dd::from_f64(value);
                    }
                }
            }
        }
        finish(sum, slip, spin, false).map_err(Into::into)
    }

    fn rank_one_split(
        &self,
        base: [f64; 2],
        delta: [f64; 2],
        inside: bool,
        ratio: f64,
    ) -> Result<Option<f64>, TriboError> {
        let active = match self.capacities {
            [a, 0.0] if a > 0.0 => 0,
            [0.0, b] if b > 0.0 => 1,
            _ => return Ok(None),
        };
        let x = to_principal(self.texture, base)?[active].abs();
        let y = to_principal(self.texture, delta)?[active].abs();
        let angle = if inside {
            det::atan2(x, y)
        } else {
            let scale = x.max(y);
            if scale == 0.0 || ratio == 0.0 {
                return Ok(None);
            }
            let sine = (x / scale) / det::hypot(x / scale, y / scale) / ratio;
            if sine >= 1.0 {
                return Ok(None);
            }
            det::asin(sine)
        };
        Ok((angle > 0.0 && angle < FRAC_PI_2).then_some(angle))
    }
}

// Half-sum and half-difference of the ORIGINAL acting support at base +/- delta.
// Rationalize n_- - n_+ using -4 b.d/(n_+ + n_-), retaining O(delta) traction
// when the two direct supports round to the same number, and O(base) mean
// traction when opposite supports nearly cancel. Shared scaling admits the
// original zero coefficient directions without squaring dimensional forces.
fn paired_support(
    capacities: [f64; 2],
    axis: [f64; 2],
    base: [f64; 2],
    delta: [f64; 2],
) -> Result<([f64; 2], [f64; 2]), TriboError> {
    let base = to_principal(axis, base)?;
    let delta = to_principal(axis, delta)?;
    let plus = support(capacities, [base[0] + delta[0], base[1] + delta[1]])?;
    let minus = support(capacities, [base[0] - delta[0], base[1] - delta[1]])?;
    let cap_scale = capacities[0].max(capacities[1]);
    let motion_scale = base[0]
        .abs()
        .max(base[1].abs())
        .max(delta[0].abs())
        .max(delta[1].abs());
    if cap_scale == 0.0 || motion_scale == 0.0 {
        return Ok(([0.0; 2], [0.0; 2]));
    }
    let mut b =
        std::array::from_fn::<_, 2, _>(|i| (capacities[i] / cap_scale) * (base[i] / motion_scale));
    let mut d =
        std::array::from_fn::<_, 2, _>(|i| (capacities[i] / cap_scale) * (delta[i] / motion_scale));
    let scale = b[0].abs().max(b[1].abs()).max(d[0].abs()).max(d[1].abs());
    if scale == 0.0 {
        return Ok(([0.0; 2], [0.0; 2]));
    }
    for x in &mut b {
        *x /= scale;
    }
    for x in &mut d {
        *x /= scale;
    }
    let np = det::hypot(b[0] + d[0], b[1] + d[1]);
    let nm = det::hypot(b[0] - d[0], b[1] - d[1]);
    let (mean, difference) = if np.min(nm) < 0.25 * np.max(nm) || np == 0.0 || nm == 0.0 {
        // Near a one-sided zero, a rounded numerator divided by the tiny norm
        // would invent a finite difference. Direct native supports are bounded
        // and well resolved here; only comparable norms need rationalization.
        (
            std::array::from_fn(|i| -0.5 * plus[i] - 0.5 * minus[i]),
            std::array::from_fn(|i| -0.5 * plus[i] + 0.5 * minus[i]),
        )
    } else {
        let sum = np + nm;
        let gamma = (Dd::from_f64(4.0) * dot(b, d) / Dd::from_f64(sum * sum)).to_f64();
        let normalized = |first: [f64; 2], second: [f64; 2]| {
            std::array::from_fn(|i| {
                let numerator = (Dd::from_f64(first[i])
                    - Dd::from_f64(second[i]) * Dd::from_f64(gamma))
                .to_f64();
                -capacities[i] * (0.5 * (numerator / np) + 0.5 * (numerator / nm))
            })
        };
        (normalized(b, d), normalized(d, b))
    };
    Ok((
        from_principal(axis, mean)?,
        from_principal(axis, difference)?,
    ))
}

fn turn(v: [f64; 2]) -> [f64; 2] {
    [-v[1], v[0]]
}
fn dot(a: [f64; 2], b: [f64; 2]) -> Dd {
    Dd::from_f64(a[0]) * Dd::from_f64(b[0]) + Dd::from_f64(a[1]) * Dd::from_f64(b[1])
}
fn cross(a: [f64; 2], b: [f64; 2]) -> f64 {
    dot([-a[1], a[0]], b).to_f64()
}
fn invalid(field: &'static str) -> TriboError {
    TriboError::InvalidInput { field }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), EllipticHertzError> {
    if cancelled() {
        Err(EllipticHertzError::Cancelled)
    } else {
        Ok(())
    }
}

fn finish(
    sum: [Dd; 4],
    slip: [f64; 2],
    spin: f64,
    at_rest: bool,
) -> Result<HertzSlidingWrench, TriboError> {
    let mut values = [0.0; 4];
    for (value, sum) in values.iter_mut().zip(sum) {
        *value = finite(sum.to_f64(), "elliptic_hertz_integral")?;
    }
    if values[3] < 0.0 {
        return Err(invalid("elliptic_hertz_negative_heat"));
    }
    let force = [values[0], values[1]];
    let residual =
        dot(force, slip) + Dd::from_f64(values[2]) * Dd::from_f64(spin) + Dd::from_f64(values[3]);
    Ok(HertzSlidingWrench {
        force_n: force,
        spin_moment_nm: values[2],
        dissipated_power_w: values[3],
        power_residual_w: finite(residual.to_f64(), "elliptic_hertz_power_residual")?,
        at_rest,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn paired_support_near_a_one_sided_zero_retains_the_original_rank_one_law() {
        for epsilon in [1e-8, 1e-12, 1e-14, f64::EPSILON] {
            let same =
                paired_support([2.0, 0.0], [1.0, 0.0], [1.0, 0.0], [1.0 - epsilon, 0.0]).unwrap();
            assert_eq!(same, ([-2.0, 0.0], [0.0, 0.0]));
            let reversed =
                paired_support([2.0, 0.0], [1.0, 0.0], [1.0, 0.0], [1.0 + epsilon, 0.0]).unwrap();
            assert_eq!(reversed, ([0.0, 0.0], [-2.0, 0.0]));
        }
    }
}
