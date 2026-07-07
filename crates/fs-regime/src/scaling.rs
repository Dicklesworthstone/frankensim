//! Canonical nondimensionalization: recommended base scales derived from
//! the problem's roles, applied to any dimensioned quantity through its
//! SI exponent vector. Scale-aware transforms exist to improve
//! CONDITIONING — the crate ships a small exact condition-number probe so
//! the improvement is a measured, ledgerable fact rather than folklore.

use crate::RegimeError;
use crate::groups::{Role, RoleInput};
use fs_math::det;
use fs_qty::QtyAny;

/// Recommended base scales (SI): every quantity nondimensionalizes as
/// `value / (m_scale^a · kg_scale^b · s_scale^c · K_scale^d · A_scale^e)`
/// with `(a..e)` its dimension exponents.
#[derive(Debug, Clone, PartialEq)]
pub struct ScalingMap {
    /// Length scale L* (m).
    pub length: f64,
    /// Mass scale M* (kg) — ρL*³ when density is known.
    pub mass: f64,
    /// Time scale T* (s) — L*/U* when a velocity is known.
    pub time: f64,
    /// Temperature scale (K); 1 when thermal effects are untracked.
    pub temperature: f64,
    /// Current scale (A); 1 unless electromagnetics enters.
    pub current: f64,
}

impl ScalingMap {
    /// Derive the canonical scales from role-tagged inputs: `L*` from
    /// Length, `T* = L*/U*` from Velocity, `M* = ρL*³` from Density.
    ///
    /// # Errors
    /// [`RegimeError::MissingRole`] without a Length; [`RegimeError::
    /// BadValue`] on non-positive scales.
    pub fn recommend(inputs: &[RoleInput]) -> Result<ScalingMap, RegimeError> {
        let find = |role: Role| inputs.iter().find(|i| i.role == role).map(|i| i.qty.value);
        let length = find(Role::Length).ok_or(RegimeError::MissingRole {
            role: "length",
            context: "scaling recommendation".to_string(),
        })?;
        let time = find(Role::Velocity).map_or_else(
            || find(Role::Frequency).map_or(1.0, f64::recip),
            |u| length / u,
        );
        let mass = find(Role::Density).map_or(1.0, |rho| rho * length * length * length);
        for (name, v) in [("length", length), ("time", time), ("mass", mass)] {
            if !(v.is_finite() && v > 0.0) {
                return Err(RegimeError::BadValue {
                    what: format!("{name} scale must be positive and finite, got {v}"),
                });
            }
        }
        Ok(ScalingMap {
            length,
            mass,
            time,
            temperature: 1.0,
            current: 1.0,
        })
    }

    /// The scale factor for a dimension vector.
    #[must_use]
    pub fn factor(&self, dims: [i8; 5]) -> f64 {
        let base = [
            self.length,
            self.mass,
            self.time,
            self.temperature,
            self.current,
        ];
        base.iter()
            .zip(dims)
            .map(|(&b, d)| b.powi(i32::from(d)))
            .product()
    }

    /// Nondimensionalize one quantity.
    #[must_use]
    pub fn apply(&self, q: QtyAny) -> f64 {
        q.value / self.factor(q.dims.0)
    }

    /// Redimensionalize a unit-free number back to SI.
    #[must_use]
    pub fn unapply(&self, value: f64, dims: [i8; 5]) -> QtyAny {
        QtyAny::new(value * self.factor(dims), fs_qty::Dims(dims))
    }
}

/// Exact 2-norm condition number of a small dense matrix (row-major,
/// n×n) via cyclic-Jacobi eigenvalues of AᵀA. Meant for fixture-scale
/// systems (n ≲ 32) — this is a measurement probe, not an LA kernel.
///
/// # Errors
/// [`RegimeError::BadValue`] for empty/singular-to-working-precision
/// input.
pub fn condition_number(a: &[f64], n: usize) -> Result<f64, RegimeError> {
    if n == 0 || a.len() != n * n {
        return Err(RegimeError::BadValue {
            what: format!("matrix shape {} vs n²={}", a.len(), n * n),
        });
    }
    // G = AᵀA (symmetric positive semidefinite).
    let mut g = vec![0.0f64; n * n];
    for i in 0..n {
        for j in 0..n {
            let mut s = 0.0;
            for k in 0..n {
                s += a[k * n + i] * a[k * n + j];
            }
            g[i * n + j] = s;
        }
    }
    // Cyclic Jacobi sweeps.
    for _ in 0..64 {
        let mut off = 0.0f64;
        for p in 0..n {
            for q in (p + 1)..n {
                off += g[p * n + q] * g[p * n + q];
            }
        }
        if off < 1e-300 {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                let gpq = g[p * n + q];
                if gpq.abs() < 1e-300 {
                    continue;
                }
                let theta = (g[q * n + q] - g[p * n + p]) / (2.0 * gpq);
                let t = theta.signum() / (theta.abs() + det::sqrt(theta * theta + 1.0));
                let c = 1.0 / det::sqrt(t * t + 1.0);
                let s = t * c;
                for k in 0..n {
                    let (gpk, gqk) = (g[p * n + k], g[q * n + k]);
                    g[p * n + k] = c * gpk - s * gqk;
                    g[q * n + k] = s * gpk + c * gqk;
                }
                for k in 0..n {
                    let (gkp, gkq) = (g[k * n + p], g[k * n + q]);
                    g[k * n + p] = c * gkp - s * gkq;
                    g[k * n + q] = s * gkp + c * gkq;
                }
            }
        }
    }
    let mut min_ev = f64::INFINITY;
    let mut max_ev = 0.0f64;
    for i in 0..n {
        let ev = g[i * n + i].max(0.0);
        min_ev = min_ev.min(ev);
        max_ev = max_ev.max(ev);
    }
    if min_ev <= 0.0 || !max_ev.is_finite() {
        return Err(RegimeError::BadValue {
            what: "matrix is singular to working precision".to_string(),
        });
    }
    Ok(det::sqrt(max_ev / min_ev))
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_qty::Dims;

    #[test]
    fn condition_of_identity_is_one() {
        let eye = [1.0, 0.0, 0.0, 1.0];
        let c = condition_number(&eye, 2).expect("cond");
        assert!((c - 1.0).abs() < 1e-12);
    }

    #[test]
    fn condition_matches_known_diagonal() {
        let d = [100.0, 0.0, 0.0, 0.5];
        let c = condition_number(&d, 2).expect("cond");
        assert!((c - 200.0).abs() < 1e-9, "got {c}");
    }

    #[test]
    fn apply_unapply_round_trip() {
        let map = ScalingMap {
            length: 0.02,
            mass: 9.6e-6,
            time: 0.066_667,
            temperature: 1.0,
            current: 1.0,
        };
        let q = QtyAny::new(101_325.0, Dims([-1, 1, -2, 0, 0]));
        let nd = map.apply(q);
        let back = map.unapply(nd, [-1, 1, -2, 0, 0]);
        assert!((back.value - q.value).abs() / q.value < 1e-12);
        assert_eq!(back.dims, q.dims);
    }
}
