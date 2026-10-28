//! Physical Kelvin--Voigt viscosity of the surviving cohesive bonds.
//!
//! The explicitly selected constitutive choice is t_v,i = (1-d) eta_i gdot_i.
//! Eta has units Pa s/m; it is an interface viscosity, not a bulk viscosity or
//! numerical damage-relaxation parameter. The original intrinsic fracture law
//! still owns damage and recoverable energy. Compression has a normal dashpot
//! while bonds remain; after complete failure only the original independent
//! compression penalty survives. Viscous dissipation never heals damage.
//!
//! A step declares a straight jump path at constant rate. Endpoint traction is
//! used for implicit equilibrium; viscous work integrates the evolving intact
//! fraction along the entire path. Multiplying endpoint traction by the jump is
//! NOT that integral when damage changes. The degradation and frozen-onset mixed
//! extension are explicit modeling choices, not a universal adhesive calibration.

use super::{
    BilinearCohesive, CohesiveError, CohesivePoint, CohesiveStep, active, det, dot, finite, norm,
    positive,
};

/// Explicit physical viscosity and caller-declared local rate applicability.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CohesiveViscosity {
    normal_pa_s_m: f64,
    tangential_pa_s_m: f64,
    maximum_rate_m_s: f64,
}
impl CohesiveViscosity {
    /// Both normal and isotropic shear viscosities must be positive, finite SI
    /// values. The rate bound applies to the Euclidean jump-rate norm in m/s.
    pub fn new(
        normal_pa_s_m: f64,
        tangential_pa_s_m: f64,
        maximum_rate_m_s: f64,
    ) -> Result<Self, CohesiveError> {
        for x in [normal_pa_s_m, tangential_pa_s_m, maximum_rate_m_s] {
            positive(x)?;
        }
        positive(normal_pa_s_m * maximum_rate_m_s)?;
        positive(tangential_pa_s_m * maximum_rate_m_s)?;
        Ok(Self {
            normal_pa_s_m,
            tangential_pa_s_m,
            maximum_rate_m_s,
        })
    }
    pub const fn normal_pa_s_m(self) -> f64 {
        self.normal_pa_s_m
    }
    pub const fn tangential_pa_s_m(self) -> f64 {
        self.tangential_pa_s_m
    }
    pub const fn maximum_rate_m_s(self) -> f64 {
        self.maximum_rate_m_s
    }
    fn coefficients(self) -> [f64; 3] {
        [
            self.normal_pa_s_m,
            self.tangential_pa_s_m,
            self.tangential_pa_s_m,
        ]
    }
}
impl BilinearCohesive {
    /// Attach physical viscosity to the original intrinsic equilibrium envelope.
    /// Intrinsic Gc must exclude this viscous loss; fitting Gc to total finite-
    /// rate work and then adding eta would count the same loss twice.
    #[must_use]
    pub fn with_viscosity(mut self, viscosity: CohesiveViscosity) -> Self {
        self.viscosity = Some(viscosity);
        self
    }
    #[must_use]
    pub const fn viscosity(&self) -> Option<CohesiveViscosity> {
        self.viscosity
    }
}

/// Timed response preserving the original equilibrium/fracture account.
#[derive(Clone, Debug, PartialEq)]
pub struct CohesiveRateStep<'a> {
    /// Original irreversible point, equilibrium traction and intrinsic energy.
    pub equilibrium: CohesiveStep<'a>,
    /// Total endpoint traction, including the physical dashpots, Pa.
    pub traction_pa: [f64; 3],
    /// Exact smooth derivative of total endpoint traction, at fixed old history
    /// and fixed dt. It is generally nonsymmetric during mixed-mode initiation.
    pub tangent_pa_m: Option<[[f64; 3]; 3]>,
    /// Same explicitly selected semismooth branch as the equilibrium owner.
    pub branch_tangent_pa_m: [[f64; 3]; 3],
    /// Nonnegative path integral of sum eta_i (1-d) gdot_i^2, J/m2.
    pub viscous_dissipated_j_m2: f64,
}
impl<'a> From<CohesiveStep<'a>> for CohesiveRateStep<'a> {
    fn from(equilibrium: CohesiveStep<'a>) -> Self {
        Self {
            traction_pa: equilibrium.traction_pa,
            tangent_pa_m: equilibrium.tangent_pa_m,
            branch_tangent_pa_m: equilibrium.branch_tangent_pa_m,
            equilibrium,
            viscous_dissipated_j_m2: 0.0,
        }
    }
}
impl CohesiveRateStep<'_> {
    /// Total analytically integrated path work, distinct from endpoint quadrature.
    #[must_use]
    pub fn total_work_j_m2(&self) -> f64 {
        self.equilibrium.energy.work + self.viscous_dissipated_j_m2
    }
}

impl<'a> CohesivePoint<'a> {
    /// Advance a straight jump path at constant rate over an explicit positive dt.
    /// No dashpot memory is invented; the accepted point retains only the original
    /// separation and damage. Its `traction_pa()` remains the equilibrium traction.
    /// The returned timed traction includes the rate for this interval only.
    pub fn trial_over_time(
        &self,
        jump_m: [f64; 3],
        dt_s: f64,
        energy_tolerance_j_m2: f64,
    ) -> Result<CohesiveRateStep<'a>, CohesiveError> {
        positive(dt_s)?;
        let equilibrium = self.trial_equilibrium(jump_m, energy_tolerance_j_m2)?;
        let Some(viscosity) = self.law.viscosity else {
            return Ok(equilibrium.into());
        };
        let delta = core::array::from_fn::<_, 3, _>(|i| jump_m[i] - self.jump[i]);
        let rate = delta.map(|v| v / dt_s);
        if norm(rate)? > viscosity.maximum_rate_m_s {
            return Err(CohesiveError::Invalid(
                "cohesive jump rate exceeds declared applicability",
            ));
        }
        let intact = integrity(&equilibrium.point);
        let derivative = integrity_gradient(&equilibrium)?;
        let eta = viscosity.coefficients();
        let mut traction = equilibrium.traction_pa;
        let mut tangent = equilibrium.branch_tangent_pa_m;
        for i in 0..3 {
            let dashpot = finite(eta[i] * rate[i])?;
            traction[i] = finite(traction[i] + intact * dashpot)?;
            for j in 0..3 {
                tangent[i][j] = finite(
                    tangent[i][j]
                        + dashpot * derivative[j]
                        + if i == j {
                            (intact * eta[i]) / dt_s
                        } else {
                            0.0
                        },
                )?;
            }
        }
        let mean_intact = mean_integrity(self, &equilibrium.point)?;
        let mut loss = 0.0;
        for i in 0..3 {
            loss = finite(loss + eta[i] * rate[i].abs() * delta[i].abs() * mean_intact)?;
        }
        finite(equilibrium.energy.work + loss)?;
        Ok(CohesiveRateStep {
            tangent_pa_m: equilibrium.tangent_pa_m.map(|_| tangent),
            branch_tangent_pa_m: tangent,
            traction_pa: traction,
            viscous_dissipated_j_m2: loss,
            equilibrium,
        })
    }
}

fn integrity(point: &CohesivePoint<'_>) -> f64 {
    if point.law.mixed_mode.is_some() {
        point.mixed_integrity()
    } else {
        point.law.intact(point.maximum)
    }
}

// Recover the same integrity derivative from the original complete tangent,
// including sensitivity of the first mixed-mode initiation. This avoids a
// second implementation of the damage evolution or an approximate derivative.
fn integrity_gradient(step: &CohesiveStep<'_>) -> Result<[f64; 3], CohesiveError> {
    let p = &step.point;
    let intact = integrity(p);
    if intact == 0.0 || intact == 1.0 {
        return Ok([0.0; 3]);
    }
    let k = match p.law.mixed_mode {
        Some(m) => [
            m.normal_stiffness_pa_m,
            m.tangential_stiffness_pa_m,
            m.tangential_stiffness_pa_m,
        ],
        None => [p.law.stiffness; 3],
    };
    let g = active(p.jump);
    let mut row = 0;
    for i in 1..3 {
        if (k[i] * g[i]).abs() > (k[row] * g[row]).abs() {
            row = i;
        }
    }
    let divisor = k[row] * g[row];
    if divisor == 0.0 {
        return Ok([0.0; 3]);
    }
    let mut derivative = [0.0; 3];
    for j in 0..3 {
        derivative[j] = finite(
            (step.branch_tangent_pa_m[row][j] - if row == j { intact * k[row] } else { 0.0 })
                / divisor,
        )?;
    }
    Ok(derivative)
}

// The active energy norm is convex along a straight jump path. Damage grows only
// beyond its final outward crossing of the retained maximum. Split at that
// crossing, complete failure and the normal compression switch. The softening
// integral is analytic: a/(b-a) * integral (b/r(s)-1) ds.
fn mean_integrity(old: &CohesivePoint<'_>, end: &CohesivePoint<'_>) -> Result<f64, CohesiveError> {
    let old_intact = integrity(old);
    if old_intact == 0.0 {
        return Ok(0.0);
    }
    let (weights, onset, failure, maximum, lower) = match old.law.mixed_mode {
        Some(p) => {
            let Some(h) = end.mixed_history else {
                return Ok(1.0);
            };
            let lower = if old.mixed_history.is_none() {
                super::mixed_mode::initiate(p, old.jump, end.jump)?.1
            } else {
                0.0
            };
            (
                [
                    det::sqrt(p.normal_stiffness_pa_m),
                    det::sqrt(p.tangential_stiffness_pa_m),
                    det::sqrt(p.tangential_stiffness_pa_m),
                ],
                h.onset_sqrt_energy,
                h.failure_sqrt_energy,
                if old.mixed_history.is_some() {
                    old.maximum
                } else {
                    h.onset_sqrt_energy
                },
                lower,
            )
        }
        None => ([1.0; 3], old.law.onset, old.law.failure, old.maximum, 0.0),
    };
    let delta = core::array::from_fn::<_, 3, _>(|i| end.jump[i] - old.jump[i]);
    let at = |s: f64| core::array::from_fn::<_, 3, _>(|i| old.jump[i] + s * delta[i]);
    let radius = |s: f64| norm(core::array::from_fn(|i| weights[i] * active(at(s))[i]));
    let threshold = maximum.max(onset);
    if radius(1.0)? <= threshold {
        return Ok(old_intact);
    }
    let crossing = |target: f64, start: f64| -> Result<f64, CohesiveError> {
        let (mut a, mut b) = (start, 1.0);
        for _ in 0..64 {
            let m = 0.5 * a + 0.5 * b;
            if m == a || m == b {
                break;
            }
            if radius(m)? <= target {
                a = m;
            } else {
                b = m;
            }
        }
        Ok(0.5 * a + 0.5 * b)
    };
    let start = crossing(threshold, lower)?;
    let finish = if radius(1.0)? >= failure {
        crossing(failure, start)?
    } else {
        1.0
    };
    let mut result = old_intact * start;
    let switch = if delta[0] != 0.0 {
        -old.jump[0] / delta[0]
    } else {
        f64::INFINITY
    };
    let cut = if switch > start && switch < finish {
        switch
    } else {
        finish
    };
    for (a, b) in [(start, cut), (cut, finish)] {
        if b > a {
            let u = core::array::from_fn(|i| weights[i] * active(at(a))[i]);
            let v = core::array::from_fn(|i| weights[i] * active(at(b))[i]);
            result += (b - a) * softening_average(u, v, onset, failure)?;
        }
    }
    finite(result)?;
    let roundoff = 256.0 * f64::EPSILON;
    if result < -roundoff || result > old_intact + roundoff {
        return Err(CohesiveError::Invalid(
            "nonrepresentable viscous path integral",
        ));
    }
    Ok(result.clamp(0.0, old_intact))
}

fn softening_average(
    a: [f64; 3],
    b: [f64; 3],
    onset: f64,
    failure: f64,
) -> Result<f64, CohesiveError> {
    let r0 = positive(norm(a)?)?;
    let u = a.map(|v| v / r0);
    let d = core::array::from_fn::<_, 3, _>(|i| (b[i] - a[i]) / r0);
    let l = norm(d)?;
    let ud = dot(u, d);
    let dd = dot(d, d);
    let inverse_change = if 2.0 * ud.abs() + dd <= 0.001 {
        // Analytic binomial integral of (1+2 u.d t+d.d t^2)^(-1/2)-1.
        // Eight terms leave < 2e-28 absolute truncation here and retain the
        // small difference from unity near complete failure.
        let mut sum = 0.0;
        let mut binomial = 1.0;
        for k in 1..=8 {
            binomial *= (-0.5 - f64::from(k - 1)) / f64::from(k);
            let mut moment = 0.0;
            let mut choose = 1.0;
            for j in 0..=k {
                if j > 0 {
                    choose *= f64::from(k - j + 1) / f64::from(j);
                }
                moment += choose * powi(2.0 * ud, k - j) * powi(dd, j) / f64::from(k + j + 1);
            }
            sum += binomial * moment;
        }
        sum
    } else {
        let end_radius = norm(core::array::from_fn(|i| u[i] + d[i]))?;
        let radial_change = (2.0 * ud + dd) / (end_radius + 1.0);
        let logarithm = det::ln(1.0 + (l + radial_change) / (1.0 + ud / l));
        logarithm / l - 1.0
    };
    let value =
        (onset / (failure - onset)) * ((failure - r0) / r0 + (failure / r0) * inverse_change);
    finite(value)
}
fn powi(x: f64, n: i32) -> f64 {
    let mut r = 1.0;
    for _ in 0..n {
        r *= x;
    }
    r
}

#[cfg(test)]
mod tests;
