//! Conformal GA Cl(4,1) (plan §7.7): points, spheres, circles, and planes
//! as blades, used where CGA is genuinely superior — tangency- and
//! sphere-rich construction (sphere chains, fillet/lip scaffolding,
//! Apollonius-type incidence). Basis: (e1, e2, e3, e+, e−); null frame
//! `n∞ = e− + e+`, `n_o = (e− − e+)/2` with `n_o · n∞ = −1`.

use crate::GaError;
use crate::facade::Vec3;
use crate::mv::Cga;

const E1: usize = 1;
const E2: usize = 2;
const E3: usize = 4;
const EP: usize = 8;
const EM: usize = 16;

/// The point at infinity `n∞`.
#[must_use]
pub fn n_inf() -> Cga {
    let mut v = Cga::zero();
    v.0[EM] = 1.0;
    v.0[EP] = 1.0;
    v
}

/// The origin null vector `n_o`.
#[must_use]
pub fn n_o() -> Cga {
    let mut v = Cga::zero();
    v.0[EM] = 0.5;
    v.0[EP] = -0.5;
    v
}

/// Conformal embedding `up(p) = n_o + p + ½‖p‖² n∞` (a null vector).
#[must_use]
pub fn up(p: Vec3) -> Cga {
    let half_sq = 0.5 * p.dot(p);
    let mut v = Cga::zero();
    v.0[E1] = p.x;
    v.0[E2] = p.y;
    v.0[E3] = p.z;
    v.0[EP] = half_sq - 0.5;
    v.0[EM] = half_sq + 0.5;
    v
}

/// Inverse embedding: normalize so `X · n∞ = −1`, then read the Euclidean
/// components.
///
/// # Errors
/// [`GaError::ZeroWeight`] when `X · n∞ = 0` (an ideal/flat element with
/// no finite representative).
pub fn down(x: &Cga) -> Result<Vec3, GaError> {
    let dot = x.gp(&n_inf()).scalar_part();
    if dot == 0.0 {
        return Err(GaError::ZeroWeight {
            context: "down(): element has no finite representative",
        });
    }
    let s = -1.0 / dot;
    Ok(Vec3::new(x.0[E1] * s, x.0[E2] * s, x.0[E3] * s))
}

/// Inverse pseudoscalar (for duality; the CGA metric is non-degenerate).
#[must_use]
fn pseudo_inverse() -> Cga {
    let i = Cga::blade(Cga::PSEUDO, 1.0);
    let denom = i.gp(&i.reverse()).scalar_part();
    i.reverse().scale(1.0 / denom)
}

/// Metric dual `X I⁻¹` (maps direct ↔ dual representations).
#[must_use]
pub fn dual(x: &Cga) -> Cga {
    x.gp(&pseudo_inverse())
}

/// Dual sphere `s = up(c) − ½r² n∞` (grade 1; `X · s = 0` iff `X = up(p)`
/// lies on the sphere).
#[must_use]
pub fn dual_sphere(center: Vec3, radius: f64) -> Cga {
    up(center).sub(&n_inf().scale(0.5 * radius * radius))
}

/// Direct sphere through four points: `P₁ ∧ P₂ ∧ P₃ ∧ P₄` (grade 4;
/// degenerate — coplanar/coincident points — yields a flat or zero blade).
#[must_use]
pub fn sphere_through(p1: Vec3, p2: Vec3, p3: Vec3, p4: Vec3) -> Cga {
    up(p1).wedge(&up(p2)).wedge(&up(p3)).wedge(&up(p4))
}

/// Direct circle through three points (grade 3).
#[must_use]
pub fn circle_through(p1: Vec3, p2: Vec3, p3: Vec3) -> Cga {
    up(p1).wedge(&up(p2)).wedge(&up(p3))
}

/// Direct plane through three points (grade 4 flat: wedge with `n∞`).
#[must_use]
pub fn plane_through(p1: Vec3, p2: Vec3, p3: Vec3) -> Cga {
    circle_through(p1, p2, p3).wedge(&n_inf())
}

/// Incidence measure of a point with a direct round/flat `X`: the ∞-norm
/// of `up(p) ∧ X` (zero iff the point lies on it).
#[must_use]
pub fn incidence(p: Vec3, x: &Cga) -> f64 {
    up(p).wedge(x).max_abs()
}

/// Center and radius of a direct (grade-4) sphere.
///
/// # Errors
/// [`GaError::ZeroWeight`] for flat/degenerate input (coplanar points).
pub fn sphere_center_radius(sphere: &Cga) -> Result<(Vec3, f64), GaError> {
    let s = dual(sphere); // grade-1 dual sphere, arbitrary scale
    let weight = -s.gp(&n_inf()).scalar_part();
    if weight == 0.0 {
        return Err(GaError::ZeroWeight {
            context: "sphere_center_radius(): flat or degenerate blade",
        });
    }
    let s = s.scale(1.0 / weight); // now s = up(c) − ½r² n∞
    let r_sq = s.gp(&s).scalar_part();
    let center_mv = s.add(&n_inf().scale(0.5 * r_sq));
    let center = down(&center_mv)?;
    Ok((center, fs_math::det::sqrt(r_sq.max(0.0))))
}

/// Tangency residual for two dual spheres: `(s₁·s₂)² − s₁² s₂²` scaled to
/// the spheres' magnitude — zero (to tolerance) iff internally or
/// externally tangent.
#[must_use]
pub fn tangency_residual(s1: &Cga, s2: &Cga) -> f64 {
    let dot = s1.gp(s2).scalar_part();
    let n1 = s1.gp(s1).scalar_part();
    let n2 = s2.gp(s2).scalar_part();
    let scale = (n1 * n2).abs().max(dot * dot).max(f64::MIN_POSITIVE);
    (dot * dot - n1 * n2) / scale
}
