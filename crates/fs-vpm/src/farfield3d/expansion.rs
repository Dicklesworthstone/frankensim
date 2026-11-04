//! Continuous-filament Cartesian moments and their absolute Taylor remainder.
//!
//! For `F(x,c)=x/(|x|²+c)^(3/2)`, the retained polynomial is
//! `F - DxF δ + DxxF[δ,δ]/2 + DcF η`. Split its remainder into spatial
//! order three, the spatial change in `η DcF`, and core order two. The
//! directional/operator derivative bounds are respectively `24/q^(5/2)`,
//! `6/q^(5/2)`, and `15/(4 q³)`, with `q=|x|²+c`. Integrating those bounds
//! along the Taylor path gives the expressions used below. Multiplication
//! by the unsigned circulation-length measure bounds the complete line
//! integrals, without treating a finite segment as a point at its midpoint.
//!
//! The spatial constant 24 follows directly from `z=x/sqrt(q)`, `t=z·e`:
//! `q^(5/2) D_e³F = (-9+45t²)e + (45t-105t³)z`. Since `|z|≤1`, its squared
//! norm is at most `81+405t²-1485t⁴+1575t⁶ ≤ 576`; the difference from 576
//! factors as `(1-t²)(495+90t²+1575t⁴)`. Also,
//! `DxDcF = -3/(2q^(5/2)) (I-5zzᵀ)` has norm at most `6/q^(5/2)`, and
//! `DccF = 15x/(4q^(7/2))` has norm at most `15/(4q³)`.

use super::{cross, dist, eps_pairs};
use crate::filament3d::WakeSegment;

/// A fixed-size cell cache, shared by retained far cells and transient trees.
/// Source storage and exact-leaf evaluation belong to the caller.
#[derive(Clone, Debug)]
pub(crate) struct Expansion {
    pub centroid: [f64; 3],
    pub radius_m: f64,
    pub w_sum: [f64; 3],
    pub w_abs_sum: f64,
    pub dip: [[f64; 3]; 3],
    pub quad: [[[f64; 3]; 3]; 3],
    pub core2_m2: f64,
    pub core_moment: [f64; 3],
    pub max_abs_core_delta_m2: f64,
    pub max_source_core2_m2: f64,
    pub source_count: usize,
}

impl Expansion {
    /// Repeated source passes preserve the retained far-cell arithmetic/order.
    /// A clone of the iterator must produce the same finite-core sources.
    pub(crate) fn from_sources(sources: impl Clone + Iterator<Item = WakeSegment>) -> Self {
        let mut weight_scale = 0.0_f64;
        let mut core_scale = 0.0_f64;
        let mut source_count = 0;
        for source in sources.clone() {
            assert!(
                source.core2_m2.is_finite() && source.core2_m2 >= 0.0,
                "finite nonnegative source cores"
            );
            weight_scale =
                weight_scale.max(source.circulation.abs() * dist(source.start, source.end));
            core_scale = core_scale.max(source.core2_m2);
            source_count += 1;
        }
        let minimum_core = sources
            .clone()
            .map(|s| s.core2_m2)
            .fold(core_scale, f64::min);
        let normalized_weight = |s: WakeSegment| {
            if weight_scale == 0.0 {
                0.0
            } else {
                (s.circulation.abs() * dist(s.start, s.end)) / weight_scale
            }
        };
        let weight_sum: f64 = sources.clone().map(normalized_weight).sum();
        let mut centroid = [0.0; 3];
        let mut scaled_core_sum = 0.0;
        if weight_sum > 0.0 {
            for s in sources.clone() {
                let weight = normalized_weight(s) / weight_sum;
                for k in 0..3 {
                    centroid[k] += weight * (0.5 * s.start[k] + 0.5 * s.end[k]);
                }
                if core_scale > 0.0 {
                    scaled_core_sum += weight * (s.core2_m2 / core_scale);
                }
            }
        }
        let core2_m2 = if minimum_core == core_scale {
            core_scale
        } else if weight_sum == 0.0 {
            sources.clone().next().map_or(0.0, |s| s.core2_m2)
        } else {
            core_scale * scaled_core_sum.min(1.0)
        };
        let mut result = Self {
            centroid,
            radius_m: 0.0,
            w_sum: [0.0; 3],
            w_abs_sum: 0.0,
            dip: [[0.0; 3]; 3],
            quad: [[[0.0; 3]; 3]; 3],
            core2_m2,
            core_moment: [0.0; 3],
            max_abs_core_delta_m2: 0.0,
            max_source_core2_m2: core_scale,
            source_count,
        };
        for s in sources {
            let dl = [
                s.end[0] - s.start[0],
                s.end[1] - s.start[1],
                s.end[2] - s.start[2],
            ];
            let w = [
                s.circulation * dl[0],
                s.circulation * dl[1],
                s.circulation * dl[2],
            ];
            let d = [
                0.5 * s.start[0] + 0.5 * s.end[0] - centroid[0],
                0.5 * s.start[1] + 0.5 * s.end[1] - centroid[1],
                0.5 * s.start[2] + 0.5 * s.end[2] - centroid[2],
            ];
            let delta_core = s.core2_m2 - core2_m2;
            result.max_abs_core_delta_m2 = result.max_abs_core_delta_m2.max(delta_core.abs());
            result.w_abs_sum += s.circulation.abs() * dist(s.start, s.end);
            for j in 0..3 {
                result.w_sum[j] += w[j];
                result.core_moment[j] += w[j] * delta_core;
                for l in 0..3 {
                    result.dip[j][l] += w[j] * d[l];
                    for m in 0..3 {
                        // Uniform line parameter s∈[-1/2,1/2] has mean(s²)=1/12.
                        result.quad[j][l][m] += w[j] * (d[l] * d[m] + dl[l] * dl[m] / 12.0);
                    }
                }
            }
            result.radius_m = result
                .radius_m
                .max(dist(s.start, centroid))
                .max(dist(s.end, centroid));
        }
        result
    }

    pub(crate) fn is_finite(&self) -> bool {
        [
            self.radius_m,
            self.w_abs_sum,
            self.core2_m2,
            self.max_abs_core_delta_m2,
            self.max_source_core2_m2,
        ]
        .iter()
        .all(|v| v.is_finite() && *v >= 0.0)
            && self
                .centroid
                .iter()
                .chain(&self.w_sum)
                .chain(&self.core_moment)
                .chain(self.dip.iter().flatten())
                .chain(self.quad.iter().flatten().flatten())
                .all(|v| v.is_finite())
    }


    /// Relative arithmetic allowances require normal unsigned scales for every
    /// retained moment order. Otherwise even nonzero geometry can have a zero
    /// cached weight/quadrupole after underflow while its exact field is finite.
    /// Such caches remain valid source containers but cannot be approximated.
    fn resolved_moment_scales(&self) -> bool {
        let w = self.w_abs_sum;
        if !(w > 0.0 && w.is_normal()) {
            return false;
        }
        if self.radius_m > 0.0 {
            let dipole_scale = w * self.radius_m;
            let quadrupole_scale = dipole_scale * self.radius_m;
            if !dipole_scale.is_normal() || !quadrupole_scale.is_normal() {
                return false;
            }
        }
        self.max_abs_core_delta_m2 == 0.0
            || (w * self.max_abs_core_delta_m2).is_normal()
    }

    /// Pure quadrupole plus leading core-variation polynomial. The caller
    /// admits its domain/error before use; this never scans exact sources.
    pub(crate) fn velocity(&self, p: [f64; 3]) -> [f64; 3] {
        let r = [
            p[0] - self.centroid[0],
            p[1] - self.centroid[1],
            p[2] - self.centroid[2],
        ];
        let distance2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
        let rho2 = distance2 + self.core2_m2;
        let rho = rho2.sqrt();
        if rho == 0.0 {
            return [0.0; 3];
        }
        let rho3 = rho2 * rho;
        let rho5 = rho3 * rho2;
        let rho7 = rho5 * rho2;
        let f0 = [r[0] / rho3, r[1] / rho3, r[2] / rho3];
        let kd = |i: usize, j: usize| -> f64 { if i == j { 1.0 } else { 0.0 } };
        let jmat = |k: usize, l: usize| -> f64 { kd(k, l) / rho3 - 3.0 * r[k] * r[l] / rho5 };
        let hmat = |i: usize, j: usize, k: usize| -> f64 {
            -3.0 * (kd(i, j) * r[k] + kd(i, k) * r[j] + kd(j, k) * r[i]) / rho5
                + 15.0 * r[i] * r[j] * r[k] / rho7
        };
        let mut v = cross(self.w_sum, f0);
        for i in 0..3 {
            let (j1, k1, j2, k2) = eps_pairs(i);
            let mut acc = 0.0;
            for l in 0..3 {
                acc += self.dip[j1][l] * jmat(k1, l) - self.dip[j2][l] * jmat(k2, l);
            }
            v[i] -= acc;
        }
        for i in 0..3 {
            let (j1, k1, j2, k2) = eps_pairs(i);
            let mut acc = 0.0;
            for l in 0..3 {
                for m in 0..3 {
                    acc +=
                        self.quad[j1][l][m] * hmat(k1, l, m) - self.quad[j2][l][m] * hmat(k2, l, m);
                }
            }
            v[i] += 0.5 * acc;
        }
        let core_correction = cross(self.core_moment, r);
        for i in 0..3 {
            v[i] -= 1.5 * core_correction[i] / rho5;
        }
        let k = 1.0 / (4.0 * core::f64::consts::PI);
        [k * v[0], k * v[1], k * v[2]]
    }


    /// Analytic derivative of exactly the velocity polynomial above. Acceptance
    /// must satisfy both the velocity and gradient remainder bounds.
    pub(crate) fn kinematics(&self, p: [f64; 3]) -> crate::filament3d::WakeKinematics {
        let r: [f64; 3] = core::array::from_fn(|i| p[i] - self.centroid[i]);
        let rho2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2] + self.core2_m2;
        let rho3 = rho2 * rho2.sqrt();
        let rho5 = rho3 * rho2;
        let rho7 = rho5 * rho2;
        let rho9 = rho7 * rho2;
        let kd = |i: usize, j: usize| -> f64 { if i == j { 1.0 } else { 0.0 } };
        let jmat = |i: usize, j: usize| kd(i, j) / rho3 - 3.0 * r[i] * r[j] / rho5;
        let hmat = |i: usize, j: usize, k: usize| {
            -3.0 * (kd(i, j) * r[k] + kd(i, k) * r[j] + kd(j, k) * r[i]) / rho5
                + 15.0 * r[i] * r[j] * r[k] / rho7
        };
        let tmat = |i: usize, j: usize, k: usize, l: usize| {
            -3.0 * (kd(i, j) * kd(k, l) + kd(i, k) * kd(j, l) + kd(i, l) * kd(j, k))
                / rho5
                + 15.0
                    * (kd(i, j) * r[k] * r[l]
                        + kd(i, k) * r[j] * r[l]
                        + kd(i, l) * r[j] * r[k]
                        + kd(j, k) * r[i] * r[l]
                        + kd(j, l) * r[i] * r[k]
                        + kd(k, l) * r[i] * r[j])
                    / rho7
                - 105.0 * r[i] * r[j] * r[k] * r[l] / rho9
        };
        let scale = 1.0 / (4.0 * core::f64::consts::PI);
        let core_cross_r = cross(self.core_moment, r);
        let mut gradient = [[0.0; 3]; 3];
        for n in 0..3 {
            let mut column = cross(self.w_sum, core::array::from_fn(|i| jmat(i, n)));
            for i in 0..3 {
                let (j1, k1, j2, k2) = eps_pairs(i);
                for l in 0..3 {
                    column[i] -=
                        self.dip[j1][l] * hmat(k1, l, n) - self.dip[j2][l] * hmat(k2, l, n);
                    for m in 0..3 {
                        column[i] += 0.5
                            * (self.quad[j1][l][m] * tmat(k1, l, m, n)
                                - self.quad[j2][l][m] * tmat(k2, l, m, n));
                    }
                }
            }
            let mut direction = [0.0; 3];
            direction[n] = 1.0;
            let core_cross_direction = cross(self.core_moment, direction);
            for i in 0..3 {
                gradient[i][n] = scale
                    * (column[i] - 1.5 * core_cross_direction[i] / rho5
                        + 7.5 * core_cross_r[i] * r[n] / rho7);
            }
        }
        crate::filament3d::WakeKinematics {
            velocity: self.velocity(p),
            gradient,
        }
    }

    /// Absolute Euclidean velocity-error bound [m/s] for the retained Taylor
    /// polynomial. Infinity requires opening the node or evaluating sources.
    ///
    /// The complete real-arithmetic remainder is bounded analytically. A
    /// separate conservative allowance accounts for unsigned accumulation and
    /// world-coordinate conditioning; this is not a directed-rounding enclosure
    /// of every floating-point operation or a time-integration error bound.
    pub(crate) fn error_bound(&self, p: [f64; 3]) -> f64 {
        if !self.is_finite() {
            return f64::INFINITY;
        }
        self.admitted_error_bound(p)
    }

    /// Same bound after the caller has admitted this immutable cache once.
    pub(crate) fn admitted_error_bound(&self, p: [f64; 3]) -> f64 {
        if !self.resolved_moment_scales() || p.iter().any(|v| !v.is_finite()) {
            return f64::INFINITY;
        }
        let r = [
            p[0] - self.centroid[0],
            p[1] - self.centroid[1],
            p[2] - self.centroid[2],
        ];
        let distance2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
        let rho2 = distance2 + self.core2_m2;
        // Match the existing far-cell exact-fallback guards. A tree must not
        // charge one cell evaluation while silently evaluating all its sources.
        if distance2 <= self.radius_m * self.radius_m
            || self.max_abs_core_delta_m2 > rho2 / 16.0
            || !rho2.is_finite()
        {
            return f64::INFINITY;
        }
        let rho = rho2.sqrt();
        let rho3 = rho2 * rho;
        let rho5 = rho3 * rho2;
        let rho7 = rho5 * rho2;
        if !rho3.is_normal() || !rho5.is_normal() || !rho7.is_normal() {
            return f64::INFINITY;
        }
        let distance = r[0].hypot(r[1]).hypot(r[2]);
        let d = distance - self.radius_m;
        let theta = self.radius_m / rho;
        if d <= 0.0 || !(0.0 <= theta && theta < 1.0) {
            return f64::INFINITY;
        }
        let one_minus = 1.0 - theta;
        let core_delta = self.max_abs_core_delta_m2;
        let spatial_q = d * d + self.core2_m2;
        let core_q = d * d + (self.core2_m2 - core_delta).max(0.0);
        let a_over_s = self.radius_m / spatial_q.sqrt();
        let spatial_derivative = 4.0 * a_over_s.powi(3) * (rho2 / spatial_q);
        // Integral of 12*a³*(1-t)²/(rho-t*a)^5 over 0≤t≤1.
        let spatial_integral = theta.powi(3) * (4.0 - 3.0 * theta) / one_minus.powi(2);
        let spatial = spatial_derivative.min(spatial_integral);
        let (mixed, core) = if core_delta == 0.0 {
            (0.0, 0.0)
        } else {
            let mixed_derivative = 6.0 * a_over_s * (core_delta / spatial_q) * (rho2 / spatial_q);
            // Integral of 6*Delta*a/(rho-t*a)^5, without cancellation
            // in the equivalent difference (1-theta)^(-4)-1.
            let polynomial = 4.0 - 6.0 * theta + 4.0 * theta * theta - theta.powi(3);
            let mixed_integral = 1.5 * (core_delta / rho2) * theta * polynomial / one_minus.powi(4);
            let core = 1.875 * (core_delta / core_q).powi(2) * (rho2 / core_q);
            (mixed_derivative.min(mixed_integral), core)
        };
        let scale = (self.w_abs_sum / rho2) / (4.0 * core::f64::consts::PI);
        if !scale.is_normal() {
            return f64::INFINITY;
        }
        // Unsigned moments bound every polynomial term even under strong
        // signed cancellation. Include source-count summation and cancellation
        // of large world coordinates when forming midpoint/target offsets.
        let eps = f64::EPSILON;
        let roundoff = 64.0 * eps * (self.source_count as f64 + 8.0);
        if roundoff >= 0.125 {
            return f64::INFINITY;
        }
        let coordinate_scale = p.iter().map(|v| v.abs()).fold(0.0_f64, f64::max)
            + self
                .centroid
                .iter()
                .map(|v| v.abs())
                .fold(0.0_f64, f64::max)
            + self.radius_m;
        let conditioned = 64.0 * eps * (coordinate_scale / rho);
        let delta = core_delta / rho2;
        let arithmetic = roundoff / (1.0 - roundoff)
            * (1.0 + 2.0 * theta + 3.0 * theta * theta + 1.5 * delta)
            + conditioned * (2.0 + 6.0 * theta + 12.0 * theta * theta + 6.0 * delta);
        let bound = scale * (spatial + mixed + core + arithmetic) * (1.0 + 64.0 * eps);
        if bound.is_finite() && bound > 0.0 {
            bound
        } else {
            f64::INFINITY
        }
    }
    /// Absolute Frobenius Jacobian-error bound [1/s], after cache admission.
    ///
    /// Differentiate the retained Taylor polynomial before bounding its
    /// remainder. With J = DxF, unit direction e, z = x/sqrt(q), t = z.e,
    /// u = t^2 and s = |z|^2, 0 <= u <= s <= 1:
    ///
    /// q^3 D_e^3 J = (45t-105t^3)I + 90t ee^T
    ///              + (45-315t^2)(ze^T+ez^T) + (-315t+945t^3)zz^T.
    ///
    /// Its squared Frobenius norm is convex in s. At s=1 it is
    /// 4050+8100u-28350u^2+37800u^3, whose difference from 21600 is
    /// (1-u)(17550+9450u+37800u^2). At s=u the difference is
    /// 225(1-u)[(63u^2-77u/2)^2 + 2219u^2/4 -147u+96] >= 0.
    /// Thus the spatial remainder coefficient is 60*sqrt(6)/6.
    ///
    /// D_e DcJ = (15/2)q^-3 [tI+ez^T+ze^T-7tzz^T].
    /// The normalized squared norm is likewise convex in s; its endpoints
    /// are 2+16u and 11u-42u^2+49u^3, both at most 18.
    /// Finally DccJ=(15/4)q^-7/2(I-7zz^T) has Frobenius norm at most
    /// (15/4)*sqrt(38)/q^(7/2). Core Taylor remainder contributes half.
    ///
    /// The resulting coefficients are 10*sqrt(6), (45/2)*sqrt(2), and
    /// (15/8)*sqrt(38). Unsigned line measure bounds all source integrals.
    /// A separate arithmetic allowance is conservative, not directed rounding.
    pub(crate) fn admitted_gradient_error_bound(&self, p: [f64; 3]) -> f64 {
        if !self.resolved_moment_scales() || p.iter().any(|v| !v.is_finite()) {
            return f64::INFINITY;
        }
        let r: [f64; 3] = core::array::from_fn(|i| p[i] - self.centroid[i]);
        let distance2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
        let rho2 = distance2 + self.core2_m2;
        if distance2 <= self.radius_m * self.radius_m
            || self.max_abs_core_delta_m2 > rho2 / 16.0
            || !rho2.is_finite()
        {
            return f64::INFINITY;
        }
        let rho = rho2.sqrt();
        let rho3 = rho2 * rho;
        let rho5 = rho3 * rho2;
        let rho7 = rho5 * rho2;
        let rho9 = rho7 * rho2;
        if [rho3, rho5, rho7, rho9].iter().any(|v| !v.is_normal()) {
            return f64::INFINITY;
        }
        let d = r[0].hypot(r[1]).hypot(r[2]) - self.radius_m;
        let theta = self.radius_m / rho;
        if d <= 0.0 || !(0.0 <= theta && theta < 1.0) {
            return f64::INFINITY;
        }
        let core_delta = self.max_abs_core_delta_m2;
        let spatial_q = d * d + self.core2_m2;
        let core_q = d * d + (self.core2_m2 - core_delta).max(0.0);
        if !spatial_q.is_normal() || !core_q.is_normal() {
            return f64::INFINITY;
        }
        // Form dimensionless ratios: q^3 can overflow even when a finite
        // physical remainder is representable. Such overflow must not look
        // like a zero error bound.
        let a_over_s = self.radius_m / spatial_q.sqrt();
        let rho_over_s3 = (rho / spatial_q.sqrt()).powi(3);
        let spatial = 10.0 * 6.0_f64.sqrt() * a_over_s.powi(3) * rho_over_s3;
        let (mixed, core) = if core_delta == 0.0 {
            (0.0, 0.0)
        } else {
            (
                22.5 * 2.0_f64.sqrt() * a_over_s * (core_delta / spatial_q) * rho_over_s3,
                1.875
                    * 38.0_f64.sqrt()
                    * (core_delta / core_q).powi(2)
                    * (rho / core_q.sqrt()).powi(3),
            )
        };
        let scale = ((self.w_abs_sum / rho2) / rho) / (4.0 * core::f64::consts::PI);
        if !scale.is_normal() {
            return f64::INFINITY;
        }
        let eps = f64::EPSILON;
        let roundoff = 256.0 * eps * (self.source_count as f64 + 16.0);
        if roundoff >= 0.125 {
            return f64::INFINITY;
        }
        let coordinate_scale = p.iter().map(|v| v.abs()).fold(0.0_f64, f64::max)
            + self
                .centroid
                .iter()
                .map(|v| v.abs())
                .fold(0.0_f64, f64::max)
            + self.radius_m;
        let conditioned = 256.0 * eps * (coordinate_scale / rho);
        let delta = core_delta / rho2;
        // Operator bounds for the four retained terms and their next spatial
        // derivatives; sqrt(3) converts this arithmetic allowance to Frobenius.
        let arithmetic = 3.0_f64.sqrt()
            * (roundoff / (1.0 - roundoff)
                * (2.0 + 24.0 * theta + 102.0 * theta * theta + 6.0 * delta)
                + conditioned * (24.0 + 204.0 * theta + 1110.0 * theta * theta + 75.0 * delta));
        let bound = scale * (spatial + mixed + core + arithmetic) * (1.0 + 256.0 * eps);
        if bound.is_finite() && bound > 0.0 {
            bound
        } else {
            f64::INFINITY
        }
    }

}
