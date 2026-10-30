//! Analytic continuation of the existing constant-panel exterior operator.
//!
//! For `exp(-i omega t)`, `k=omega/c` and `A(k) p = B(k) q`, with
//! `q=i rho omega v`. Both signs of Im(k) are admitted: outgoing resonances
//! grow in space and decay in time. The matrices below retain the parent's
//! triangle-CBIE or centroid/disc Burton--Miller discretization, including its
//! static row closures. The derivative is analytic, never a frequency difference.
//!
//! A boundary impedance `A^-1 B` is meromorphic. A contour method requiring a
//! holomorphic operator must retain the boundary unknowns in an augmented pencil
//! (and exclude k=0 for Burton--Miller), rather than count zeros of that quotient.
//! Burton--Miller removes real irregular frequencies; it does NOT certify that
//! an arbitrary complex contour contains no fictitious representation roots.
//! Spatial/formulation convergence and physical mode identification remain with
//! the caller. No fitting, damping floor, symmetrization or passivity repair.
use super::super::*;
use super::DampedWaveNumber;

// The continuation evaluates both values and derivatives on one fixed geometric
// traversal. The real-frequency owner is left untouched; tests compare the
// resulting real-axis matrices with that independent original assembly bitwise.
#[path = "helmholtz_complex/quadrature.rs"]
mod quadrature;
use quadrature::{triangle_remainder_rule, triangle_self_rule};

/// Explicit logical storage and deterministic assembly-work ceilings.
#[derive(Debug, Clone, Copy)]
pub struct ComplexOptions {
    /// Additional panel cap, no greater than the parent's dense hard cap.
    pub max_panels: usize,
    /// Payload ceiling: four dense complex matrices (64*n*n bytes) with
    /// derivatives, or two (32*n*n bytes) for [`assemble_values_with_cancel`].
    /// Excludes caller geometry, allocator overhead and any later factorization.
    pub max_bytes: usize,
    /// Conservative kernel-evaluation units: 512*n*n for triangle CBIE,
    /// 8*n*n for centroid/disc rows. This is a deterministic admission measure,
    /// not a wall-clock promise or a floating-point-operation measurement.
    pub max_work: u64,
}

/// Typed admission, numerical and cancellation failures; no partial operator.
#[derive(Debug, Clone, PartialEq)]
pub enum ComplexError {
    /// Existing surface, resolution or physical input refusal.
    Boundary(HelmholtzError),
    /// Invalid complex frequency, options, or nonfinite arithmetic.
    Invalid { what: &'static str },
    /// Requested number of panels exceeds the caller's cap.
    PanelBudget { required: usize, maximum: usize },
    /// Four-matrix payload or declared kernel work exceeds its ceiling.
    Budget {
        required_bytes: usize,
        maximum_bytes: usize,
        required_work: u64,
        maximum_work: u64,
    },
    /// A fallible bounded allocation failed.
    Allocation,
    /// Cancellation was observed before allocation or between target rows.
    Cancelled,
}
impl core::fmt::Display for ComplexError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Boundary(e) => e.fmt(f),
            Self::Invalid { what } => write!(f, "FS-BEM-COMPLEX-INVALID: {what}"),
            Self::PanelBudget { required, maximum } => write!(
                f,
                "FS-BEM-COMPLEX-PANEL-BUDGET: {required} panels exceeds {maximum}"
            ),
            Self::Budget {
                required_bytes,
                maximum_bytes,
                required_work,
                maximum_work,
            } => write!(
                f,
                "FS-BEM-COMPLEX-BUDGET: bytes {required_bytes}/{maximum_bytes}, kernel work {required_work}/{maximum_work}"
            ),
            Self::Allocation => write!(
                f,
                "FS-BEM-COMPLEX-ALLOCATION: matrix payload could not be reserved"
            ),
            Self::Cancelled => write!(f, "FS-BEM-COMPLEX-CANCELLED: assembly cancelled"),
        }
    }
}
impl std::error::Error for ComplexError {}
impl From<HelmholtzError> for ComplexError {
    fn from(value: HelmholtzError) -> Self {
        Self::Boundary(value)
    }
}

/// Original row-major boundary equations and derivatives with respect to k
/// [1/m]. In an omega [rad/s] pencil multiply both derivatives by 1/c.
#[derive(Debug, Clone)]
pub struct ComplexOperator {
    /// Boundary panel count.
    pub dimension: usize,
    /// Left matrix in A p=B q.
    pub a: Vec<C64>,
    /// Right matrix in A p=B q; its area integration is already included.
    pub b: Vec<C64>,
    /// Analytic dA/dk at the same k and on the same geometric quadrature.
    pub da_dk: Vec<C64>,
    /// Analytic dB/dk at the same k and on the same geometric quadrature.
    pub db_dk: Vec<C64>,
}

/// The same boundary pair without derivative storage or derivative evaluation.
/// Frequency-domain material solves need only A p=B q.
#[derive(Debug, Clone)]
pub struct ComplexBoundaryPair {
    pub dimension: usize,
    pub a: Vec<C64>,
    pub b: Vec<C64>,
}

fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}
fn invalid(what: &'static str) -> ComplexError {
    ComplexError::Invalid { what }
}
fn zeros(n: usize) -> Result<Vec<C64>, ComplexError> {
    let mut out = Vec::new();
    out.try_reserve_exact(n)
        .map_err(|_| ComplexError::Allocation)?;
    out.resize(n, C64::ZERO);
    Ok(out)
}

/// Assemble the original dense boundary pair and its analytic derivative.
///
/// Triangle `PlainCbie` shares the parent's exact static integrals, projected
/// near-face fan, 8x8 near/Duffy and 4x4 far remainder rules. Burton--Miller
/// shares its established centroid/disc rule; the separate real-frequency
/// `PanelRule::TriangleNear` is not silently substituted. Positive real k values
/// reproduce the original default assembly bitwise. Complex frequencies use
/// the same outgoing Green function continued into either half-plane.
///
/// Cancellation is polled before admission/allocation, at every target row,
/// and before publication. No LU or boundary-mode elimination occurs here.
/// # Errors
/// Zero/nonfinite k, invalid surface/options, insufficient resolution at |k|,
/// budgets, allocation, nonfinite kernels or cancellation.
pub fn assemble_with_cancel(
    surface: &SpherePanels,
    k: C64,
    formulation: Formulation,
    options: &ComplexOptions,
    cancelled: impl FnMut() -> bool,
) -> Result<ComplexOperator, ComplexError> {
    assemble::<true>(surface, k, formulation, options, cancelled)
}

/// Assemble only the original A/B values, on identical quadrature and with
/// identical admission/cancellation to [`assemble_with_cancel`]. The requested
/// matrix payload is halved; the conservative work ceiling is unchanged.
pub fn assemble_values_with_cancel(
    surface: &SpherePanels,
    k: C64,
    formulation: Formulation,
    options: &ComplexOptions,
    cancelled: impl FnMut() -> bool,
) -> Result<ComplexBoundaryPair, ComplexError> {
    let out = assemble::<false>(surface, k, formulation, options, cancelled)?;
    Ok(ComplexBoundaryPair {
        dimension: out.dimension,
        a: out.a,
        b: out.b,
    })
}

fn assemble<const DERIVATIVES: bool>(
    surface: &SpherePanels,
    k: C64,
    formulation: Formulation,
    options: &ComplexOptions,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ComplexOperator, ComplexError> {
    if cancelled() {
        return Err(ComplexError::Cancelled);
    }
    if !finite(k) || !k.abs().is_finite() || k.abs() == 0.0 {
        return Err(invalid("wavenumber must be finite and nonzero"));
    }
    if options.max_panels == 0 || options.max_panels > MAX_DENSE_PANELS {
        return Err(invalid("panel budget must be in 1..=MAX_DENSE_PANELS"));
    }
    surface
        .validate()
        .map_err(|_| invalid("invalid closed source surface"))?;
    admit_resolution(surface, k.abs())?;
    let n = surface.centroids().len();
    if n > options.max_panels {
        return Err(ComplexError::PanelBudget {
            required: n,
            maximum: options.max_panels,
        });
    }
    let triangles = if formulation == Formulation::PlainCbie {
        surface.triangles()
    } else {
        None
    };
    let entries = n
        .checked_mul(n)
        .ok_or_else(|| invalid("matrix dimension overflow"))?;
    let required_bytes = entries
        .checked_mul((if DERIVATIVES { 4 } else { 2 }) * core::mem::size_of::<C64>())
        .ok_or_else(|| invalid("matrix payload overflow"))?;
    let required_work = (entries as u64)
        .checked_mul(if triangles.is_some() { 512 } else { 8 })
        .ok_or_else(|| invalid("matrix work overflow"))?;
    if required_bytes > options.max_bytes || required_work > options.max_work {
        return Err(ComplexError::Budget {
            required_bytes,
            maximum_bytes: options.max_bytes,
            required_work,
            maximum_work: options.max_work,
        });
    }
    let mut out = ComplexOperator {
        dimension: n,
        a: zeros(entries)?,
        b: zeros(entries)?,
        da_dk: if DERIVATIVES {
            zeros(entries)?
        } else {
            Vec::new()
        },
        db_dk: if DERIVATIVES {
            zeros(entries)?
        } else {
            Vec::new()
        },
    };
    let sign = match formulation {
        Formulation::PlainCbie => 0.0,
        Formulation::BurtonMiller => 1.0,
        Formulation::BurtonMillerWrongAlphaSign => -1.0,
    };
    // Preserve the original real-axis alpha operation, including signed zeros.
    let alpha = if k.im == 0.0 {
        alpha_for(formulation, k.re)
    } else {
        C64::new(0.0, sign) * k.recip()
    };
    let dalpha = if sign == 0.0 {
        C64::ZERO
    } else {
        -(alpha * k.recip())
    };
    let half = C64::from_re(0.5);
    let centroids = surface.centroids();
    let normals = surface.normals();
    let areas = surface.areas();
    for i in 0..n {
        if cancelled() {
            return Err(ComplexError::Cancelled);
        }
        let (xi, ni) = (centroids[i], normals[i]);
        let (mut n0_row, mut d0_row) = (C64::ZERO, C64::ZERO);
        if sign != 0.0 {
            for j in 0..n {
                if j != i {
                    let (_, d0, _, n0) = kernels(0.0, xi, ni, centroids[j], normals[j]);
                    n0_row = n0_row + n0.scale(areas[j]);
                    d0_row = d0_row + d0.scale(areas[j]);
                }
            }
        }
        for j in 0..n {
            let (value, derivative) = if i == j {
                let (s, nr, ds, dnr) = if let Some(ts) = triangles {
                    let (s, ds) = triangle_single_self::<DERIVATIVES>(k, xi, ts[i]);
                    (s, C64::ZERO, ds, C64::ZERO)
                } else {
                    disc_with_derivative::<DERIVATIVES>(k, areas[i])
                };
                let dself = if sign == 0.0 {
                    C64::ZERO
                } else {
                    -half - d0_row
                };
                (
                    [s, dself, dself, nr - n0_row],
                    [ds, C64::ZERO, C64::ZERO, dnr],
                )
            } else if let Some(ts) = triangles {
                let (s, d, ds, dd) = if triangle_is_near(xi, centroids[j], ts[j]) {
                    triangle_weak::<DERIVATIVES, 8>(k, xi, ts[j], normals[j], (&GL8_X, &GL8_W))
                } else {
                    triangle_weak::<DERIVATIVES, 4>(k, xi, ts[j], normals[j], (&GL4_X, &GL4_W))
                };
                ([s, d, C64::ZERO, C64::ZERO], [ds, dd, C64::ZERO, C64::ZERO])
            } else {
                let (v, d) =
                    point_with_derivative::<DERIVATIVES>(k, xi, ni, centroids[j], normals[j]);
                (v.map(|x| x.scale(areas[j])), d.map(|x| x.scale(areas[j])))
            };
            let [s, d, dp, hyp] = value;
            let [ds, dd, ddp, dhyp] = derivative;
            let mut av = d + alpha * hyp;
            let mut bv = s + alpha * dp;
            let da = dd + dalpha * hyp + alpha * dhyp;
            let mut db = ds + dalpha * dp + alpha * ddp;
            if i == j {
                av = av - half;
                bv = bv + alpha * half;
                db = db + dalpha * half;
            }
            if !finite(av) || !finite(bv) || (DERIVATIVES && (!finite(da) || !finite(db))) {
                return Err(invalid("nonfinite complex boundary kernel or derivative"));
            }
            let ij = i * n + j;
            out.a[ij] = av;
            out.b[ij] = bv;
            if DERIVATIVES {
                out.da_dk[ij] = da;
                out.db_dk[ij] = db;
            }
        }
    }
    if cancelled() {
        return Err(ComplexError::Cancelled);
    }
    Ok(out)
}

fn wave(k: C64) -> DampedWaveNumber {
    DampedWaveNumber {
        k,
        phase_origin_m: 0.0,
    }
}

fn point_with_derivative<const DERIVATIVES: bool>(
    k: C64,
    x: [f64; 3],
    nx: [f64; 3],
    y: [f64; 3],
    ny: [f64; 3],
) -> ([C64; 4], [C64; 4]) {
    let (g, dy, dx, hyp) = if k.im == 0.0 {
        kernels(k.re, x, nx, y, ny)
    } else {
        kernels(wave(k), x, nx, y, ny)
    };
    if !DERIVATIVES {
        return ([g, dy, dx, hyp], [C64::ZERO; 4]);
    }
    let d = sub(x, y);
    let r = norm(d);
    let imaginary = C64::new(0.0, 1.0);
    let radial = C64::new(-1.0 / r - k.im, k.re);
    let gp = g * radial;
    let dg = g * imaginary.scale(r);
    let dgp = dg * radial + g * imaginary;
    let dgpp = dgp * radial + gp * imaginary + dg.scale(1.0 / (r * r));
    let nyd = dot(ny, d);
    let nxd = dot(nx, d);
    let dradial = dgpp.scale(1.0 / (r * r)) - dgp.scale(1.0 / (r * r * r));
    (
        [g, dy, dx, hyp],
        [
            dg,
            dgp.scale(-nyd / r),
            dgp.scale(nxd / r),
            dradial.scale(-nxd * nyd) - dgp.scale(dot(nx, ny) / r),
        ],
    )
}

fn exp_minus_one(z: C64) -> C64 {
    if z.abs() < 0.25 {
        let mut term = z;
        let mut result = term;
        for m in 2..=24 {
            term = (term * z).scale(1.0 / f64::from(m));
            result = result + term;
        }
        result
    } else {
        expik(z.im).scale(det::exp(z.re)) - C64::ONE
    }
}

fn disc_with_derivative<const DERIVATIVES: bool>(k: C64, area: f64) -> (C64, C64, C64, C64) {
    let (s, nr) = if k.im == 0.0 {
        self_terms(k.re, area)
    } else {
        self_terms(wave(k), area)
    };
    if !DERIVATIVES {
        return (s, nr, C64::ZERO, C64::ZERO);
    }
    let a = det::sqrt(area / core::f64::consts::PI);
    let imaginary = C64::new(0.0, 1.0);
    let z = (imaginary * k).scale(a);
    let phi_prime = if z.abs() < 0.25 {
        let mut sum = C64::from_re(0.5);
        let mut term = sum;
        for m in 1..24 {
            // m z^(m-1)/(m+1)! -> (m+1) z^m/(m+2)!.
            term = (term * z).scale(f64::from(m + 1) / f64::from(m * (m + 2)));
            sum = sum + term;
        }
        sum
    } else {
        (z + (z - C64::ONE) * exp_minus_one(z)) * (z * z).recip()
    };
    let ds = (imaginary * phi_prime).scale(0.5 * a * a);
    let dnr = (imaginary * exp_minus_one(z)).scale(-0.5);
    (s, nr, ds, dnr)
}

// Exact static subtraction makes the integrands regular on close opposite skins.
// d/dk of the single remainder is i exp(ikr)/(4pi); d/dk of the double
// remainder is k times the single remainder. No differencing scale is needed.
fn single_remainder(k: C64, r: f64) -> C64 {
    if k.im == 0.0 {
        regularized_single_layer(k.re, r)
    } else {
        green_remainder(wave(k), r)
    }
}
fn single_derivative(k: C64, r: f64) -> C64 {
    (C64::new(0.0, 1.0) * wave(k).outgoing(r)).scale(1.0 / (4.0 * core::f64::consts::PI))
}
fn double_remainder(k: C64, r: f64) -> C64 {
    if k.im == 0.0 {
        return regularized_double_layer_remainder(k.re, r);
    }
    let ik = C64::new(0.0, 1.0) * k;
    let z = ik.scale(r);
    let four_pi = 4.0 * core::f64::consts::PI;
    if z.abs() < 0.25 {
        // [exp(z)(1-z)-1+z^2/2]/r^3 = SUM_(m>=3) (1-m) z^m/(m! r^3).
        let mut term = (ik * ik * ik).scale(-1.0 / 3.0);
        let mut sum = term;
        for m in 3..24 {
            term = (term * z).scale(f64::from(m) / f64::from((m + 1) * (m - 1)));
            sum = sum + term;
        }
        sum.scale(1.0 / four_pi)
    } else {
        (wave(k).outgoing(r) * (C64::ONE - z) - C64::ONE + (z * z).scale(0.5))
            .scale(1.0 / (four_pi * r * r * r))
    }
}

fn triangle_single_self<const DERIVATIVES: bool>(
    k: C64,
    x: [f64; 3],
    triangle: [[f64; 3]; 3],
) -> (C64, C64) {
    let mut s = C64::from_re(triangle_self_static(x, triangle));
    let mut ds = C64::ZERO;
    triangle_self_rule(x, triangle, |r, weight| {
        s = s + single_remainder(k, r).scale(weight);
        if DERIVATIVES {
            ds = ds + single_derivative(k, r).scale(weight);
        }
    });
    (s, ds)
}
fn triangle_weak<const DERIVATIVES: bool, const N: usize>(
    k: C64,
    x: [f64; 3],
    triangle: [[f64; 3]; 3],
    normal: [f64; 3],
    rule: (&[f64; N], &[f64; N]),
) -> (C64, C64, C64, C64) {
    let (ss, sd) = triangle_static_influence(x, triangle, normal);
    let normal_offset = -dot(sub(triangle[0], x), normal);
    let mut s = C64::from_re(ss);
    let mut ds = C64::ZERO;
    let mut d = if k.im == 0.0 {
        C64::from_re(sd + 0.5 * k.re * k.re * normal_offset * ss)
    } else {
        C64::from_re(sd) + (k * k).scale(0.5 * normal_offset * ss)
    };
    triangle_remainder_rule(x, triangle, normal, rule, |r, weight| {
        let g = single_remainder(k, r);
        let dgdny = double_remainder(k, r).scale(normal_offset);
        s = s + g.scale(weight);
        d = d + dgdny.scale(weight);
        if DERIVATIVES {
            ds = ds + single_derivative(k, r).scale(weight);
        }
    });
    (s, d, ds, (k * s).scale(normal_offset))
}

#[cfg(test)]
#[path = "helmholtz_complex/tests.rs"]
mod tests;
