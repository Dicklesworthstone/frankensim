//! Three-dimensional energy-based magnetic history with a true vector return map.
//!
//! Each original branch minimizes
//! `|z|²/(2 chi) - H.z + hc |z-z_old|` on the ball `|z| <= zs`.
//! The Euclidean friction norm couples all three components. The unsaturated
//! minimizer is exact shrinkage; an active saturation ball is solved through
//! its scalar Lagrange multiplier, not by projecting an unconstrained update.
//! This is the quadratic-potential specialization of the variational dry-friction
//! construction in Prigozhin et al., IEEE Trans. Magn. 52 (2016), §§II–III,
//! <https://arxiv.org/abs/1605.08748>. Its saturation ball is an explicit model
//! choice; the paper's nonlinear anhysteretic fit is not inferred.
//!
//! `B = mu_b H + sum(z)` and `W = mu_b |H|²/2 + sum(|z|²/(2 chi))`.
//! The incremental pinning heat is `Q = sum(hc |delta z|)`. Endpoint work
//! `H_new.delta B = delta W + Q + D` includes a separately reported nonnegative
//! discretization remainder D; D is NOT thermal heat. Unlike a scalar monotone
//! segment, a finite rotating-field increment does not determine an exact
//! continuous path integral. Resolve turns/reversals and refine the history.
//!
//! An explicitly fitted branch viscosity `eta >= 0` [m s/H] adds
//! `eta |z-z_old|²/(2 dt)` to that SAME incremental objective. Its backward
//! Euler force is `eta (z-z_old)/dt`; its physical switching heat is the FULL
//! `eta |z-z_old|²/dt`, while the original stored energy is unchanged. This
//! intrinsic polarization kinetics does not represent Ohmic eddy-current loss.
//!
//! Parameters are explicit caller fits at fixed temperature/stress. This law
//! supplies neither empirical material authority, eddy/excess loss, texture
//! anisotropy, magnetic damage, nor automatic spatial/circuit/thermal coupling.

use super::{FieldStorageError, HysteresisBranch, Parameters, ScalarHysteresis};
use fs_la::factor::cholesky;
use std::sync::Arc;

type Vector = [f64; 3];
type Matrix = [[f64; 3]; 3];

/// Physical residual allowances and bounded local work. All tolerances are SI.
#[derive(Clone, Copy, Debug)]
pub struct VectorHysteresisControls {
    pub relative_tolerance: f64,
    pub absolute_field_tolerance_a_m: f64,
    pub absolute_flux_tolerance_t: f64,
    pub work_tolerance_j_m3: f64,
    pub maximum_branch_iterations: usize,
    pub maximum_inverse_iterations: usize,
    pub maximum_backtracks: usize,
}

#[derive(Debug)]
pub enum VectorHysteresisError {
    Material(FieldStorageError),
    Invalid(&'static str),
    Cancelled,
    /// Requested |B| exceeds even saturated polarization plus the supported
    /// background field. A loading continuation may retry a smaller increment.
    FluxRange {
        observed_t: f64,
        maximum_t: f64,
    },
    BranchNonconvergence {
        branch: usize,
        iterations: usize,
    },
    FluxNonconvergence {
        iterations: usize,
        residual_t: f64,
    },
}
impl From<FieldStorageError> for VectorHysteresisError {
    fn from(value: FieldStorageError) -> Self {
        Self::Material(value)
    }
}
impl core::fmt::Display for VectorHysteresisError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "vector magnetic history: {self:?}")
    }
}
impl std::error::Error for VectorHysteresisError {}
type Result<T> = core::result::Result<T, VectorHysteresisError>;

/// Fixed isotropic branch law; the original scalar parameter admission is reused.
#[derive(Clone, Debug, PartialEq)]
pub struct VectorHysteresis {
    parameters: Arc<Parameters>,
    viscosities_m_s_h: Option<Arc<[f64]>>,
}

/// Private accepted three-dimensional polarization memory and its physical state.
#[derive(Clone, Debug, PartialEq)]
pub struct VectorHysteresisState {
    parameters: Arc<Parameters>,
    viscosities_m_s_h: Option<Arc<[f64]>>,
    field_a_m: Vector,
    flux_density_t: Vector,
    polarizations_t: Vec<Vector>,
    stored_energy_j_m3: f64,
}
impl VectorHysteresisState {
    pub const fn field_a_m(&self) -> Vector {
        self.field_a_m
    }
    pub const fn flux_density_t(&self) -> Vector {
        self.flux_density_t
    }
    pub fn polarizations_t(&self) -> &[Vector] {
        &self.polarizations_t
    }
    pub const fn stored_energy_j_m3(&self) -> f64 {
        self.stored_energy_j_m3
    }
}

/// Endpoint response and separately owned irreversible and numerical transfers.
#[derive(Clone, Debug, PartialEq)]
pub struct VectorHysteresisStep {
    pub state: VectorHysteresisState,
    /// Total physical switching heat, pinning plus viscosity; deposit once [J/m³].
    pub irreversible_heat_j_m3: f64,
    /// Nonnegative dry-friction heat `sum(hc |delta z|)` [J/m³].
    pub pinning_heat_j_m3: f64,
    /// Nonnegative viscous heat `sum(eta |delta z|²/dt)` [J/m³].
    /// This is intrinsic switching heat, excluding resolved Ohmic eddy loss.
    pub viscous_heat_j_m3: f64,
    /// Stable difference of the complete stored energy [J/m³].
    pub energy_change_j_m3: f64,
    /// H_new dot the actually retained change of B [J/m³].
    pub endpoint_work_j_m3: f64,
    /// Nonnegative endpoint-discretization remainder, never physical heat [J/m³].
    pub numerical_dissipation_j_m3: f64,
    /// Endpoint work minus stored change, physical heat and numerical remainder.
    pub energy_residual_j_m3: f64,
    /// Frozen-history dB/dH [H/m], symmetric positive definite. At pinning and
    /// saturation kinks the selected active set is not a smoothness claim.
    pub incremental_permeability_h_m: Matrix,
    pub maximum_constitutive_residual_a_m: f64,
    pub inverse_iterations: usize,
    pub branch_iterations: usize,
    pub inverse_residual_t: f64,
}

impl VectorHysteresisStep {
    /// Frozen-history potential `W(B,z_new) + Q_pinning + Q_viscous/2`.
    /// Its B-gradient is the endpoint H and its Hessian is the incremental
    /// reluctivity. A field Newton solve must minimize this potential rather
    /// than reversible storage alone. This is an incremental variational
    /// functional; it is not path-independent physical stored energy.
    pub fn incremental_potential_j_m3(&self) -> Result<f64> {
        finite(
            self.state.stored_energy_j_m3 + self.pinning_heat_j_m3 + 0.5 * self.viscous_heat_j_m3,
        )
    }

    /// Same frozen-history potential relative to the accepted stored energy.
    /// This stable form avoids subtracting two large state energies when
    /// globalizing a small increment about a magnetized operating point.
    pub fn incremental_potential_change_j_m3(&self) -> Result<f64> {
        finite(self.energy_change_j_m3 + self.pinning_heat_j_m3 + 0.5 * self.viscous_heat_j_m3)
    }

    /// Consistent frozen-history `dH/dB` [m/H] for a B-driven field solver.
    /// The positive background supplies invertibility even while branches are
    /// pinned or saturated. Kinks retain the original selected active set.
    pub fn incremental_reluctivity_m_h(&self) -> Result<Matrix> {
        inverse_tangent(self.incremental_permeability_h_m)
    }
}

/// One native SPD inverse for magnetic and electric constitutive facades.
pub(crate) fn inverse_tangent(tangent: Matrix) -> Result<Matrix> {
    let flat: [f64; 9] = std::array::from_fn(|k| tangent[k / 3][k % 3]);
    let factor = cholesky(&flat, 3)
        .map_err(|_| invalid("vector incremental permeability is not representably SPD"))?;
    let mut inverse = [[0.0; 3]; 3];
    for column in 0..3 {
        let mut unit = [0.0; 3];
        unit[column] = 1.0;
        factor.solve(&mut unit);
        for row in 0..3 {
            inverse[row][column] = finite(unit[row])?;
        }
    }
    // The exact Hessian is symmetric; erase only triangular-solve roundoff.
    for row in 0..3 {
        for column in 0..row {
            let value = f64::midpoint(inverse[row][column], inverse[column][row]);
            inverse[row][column] = value;
            inverse[column][row] = value;
        }
    }
    Ok(inverse)
}

fn invalid(s: &'static str) -> VectorHysteresisError {
    VectorHysteresisError::Invalid(s)
}
fn finite(x: f64) -> Result<f64> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(invalid("vector magnetic arithmetic is unrepresentable"))
    }
}
fn norm(x: Vector) -> f64 {
    x[0].hypot(x[1]).hypot(x[2])
}
fn dot(x: Vector, y: Vector) -> f64 {
    x[0].mul_add(y[0], x[1].mul_add(y[1], x[2] * y[2]))
}
fn sub(a: Vector, b: Vector) -> Vector {
    std::array::from_fn(|i| a[i] - b[i])
}
fn mv(a: Matrix, x: Vector) -> Vector {
    a.map(|row| dot(row, x))
}
// These analytic Hessians are symmetric. Evaluate each pair only once so
// floating-point association cannot violate the native assembly's exact
// symmetry contract; do not average or repair an arbitrary input tensor.
fn symmetric_matrix(value: impl Fn(usize, usize) -> f64) -> Matrix {
    let mut result = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..=i {
            let entry = value(i, j);
            result[i][j] = entry;
            result[j][i] = entry;
        }
    }
    result
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(VectorHysteresisError::Cancelled)
    } else {
        Ok(())
    }
}
fn field_target(c: VectorHysteresisControls, scale: f64) -> f64 {
    c.absolute_field_tolerance_a_m
        .max(c.relative_tolerance * scale)
}
fn admit_controls(c: VectorHysteresisControls) -> Result<()> {
    if !c.relative_tolerance.is_finite()
        || !(0.0..1.0).contains(&c.relative_tolerance)
        || [
            c.absolute_field_tolerance_a_m,
            c.absolute_flux_tolerance_t,
            c.work_tolerance_j_m3,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
    {
        return Err(invalid(
            "positive finite physical tolerances and relative tolerance below one required",
        ));
    }
    Ok(())
}

struct BranchResponse {
    z: Vector,
    tangent: Matrix,
    multiplier_m_h: f64,
    residual_a_m: f64,
    iterations: usize,
}

// Exact unconstrained minimizer at a = 1/chi + lambda. Its H derivative is
// D=(I-hc/|r|*(I-n*n^T))/a, r=H-a*z_old. Also dz/da=-D*z.
fn unconstrained(old: Vector, h: Vector, a: f64, hc: f64) -> Result<(Vector, Matrix)> {
    let r = std::array::from_fn(|i| (-a).mul_add(old[i], h[i]));
    shrinkage(old, r, a, hc)
}

fn shrinkage(old: Vector, r: Vector, a: f64, hc: f64) -> Result<(Vector, Matrix)> {
    if hc == 0.0 {
        // No pinning cusp exists at r=0. A viscous branch still has compliance
        // 1/a there, including during an initially unpolarized field hold.
        let z = std::array::from_fn(|i| (r[i] / a) + old[i]);
        let tangent = symmetric_matrix(|i, j| if i == j { a.recip() } else { 0.0 });
        for value in z.into_iter().chain(tangent.into_iter().flatten()) {
            finite(value)?;
        }
        return Ok((z, tangent));
    }
    let length = finite(norm(r))?;
    if length <= hc {
        return Ok((old, [[0.0; 3]; 3]));
    }
    let direction = r.map(|x| x / length);
    let fraction = hc / length;
    let z = std::array::from_fn(|i| ((1.0 - fraction) / a).mul_add(r[i], old[i]));
    let tangent = symmetric_matrix(|i, j| {
        ((if i == j { 1.0 - fraction } else { 0.0 }) + fraction * direction[i] * direction[j]) / a
    });
    for value in z.into_iter().chain(tangent.into_iter().flatten()) {
        finite(value)?;
    }
    Ok((z, tangent))
}

// The same shrinkage and constrained return map serve both laws. For eta=0,
// retain the original a=1/chi+lambda arithmetic exactly. For eta>0 the root
// coordinate is lambda itself: forming H+k*z_old would lose the physical drive
// and inflate residual tolerances when k=eta/dt is large.
fn branch_unconstrained(
    old: Vector,
    h: Vector,
    inverse_chi: f64,
    viscous_stiffness: f64,
    normal_coordinate: f64,
    hc: f64,
) -> Result<(Vector, Matrix)> {
    if viscous_stiffness == 0.0 {
        unconstrained(old, h, normal_coordinate, hc)
    } else {
        let drive = std::array::from_fn(|i| {
            (-normal_coordinate).mul_add(old[i], (-inverse_chi).mul_add(old[i], h[i]))
        });
        let stiffness = finite(inverse_chi + viscous_stiffness + normal_coordinate)?;
        shrinkage(old, drive, stiffness, hc)
    }
}

fn stationarity_residual(
    old: Vector,
    h: Vector,
    z: Vector,
    chi: f64,
    hc: f64,
    lambda: f64,
    viscous_stiffness: f64,
) -> Result<f64> {
    let increment = sub(z, old);
    let distance = finite(norm(increment))?;
    Ok(if distance == 0.0 {
        (norm(std::array::from_fn(|i| h[i] - z[i] / chi - lambda * z[i])) - hc).max(0.0)
    } else {
        norm(std::array::from_fn(|i| {
            let dry = h[i] - z[i] / chi - lambda * z[i] - hc * increment[i] / distance;
            if viscous_stiffness == 0.0 {
                dry
            } else {
                (-viscous_stiffness).mul_add(increment[i], dry)
            }
        }))
    })
}

fn branch_response(
    branch: &HysteresisBranch,
    old: Vector,
    h: Vector,
    viscous_stiffness: f64,
    c: VectorHysteresisControls,
    index: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<BranchResponse> {
    poll(cancelled)?;
    let chi = branch.susceptibility_h_m;
    let hc = branch.pinning_field_a_m;
    let zs = branch.saturation_polarization_t;
    let a0 = chi.recip();
    let old_radius = norm(old);
    let target = field_target(c, norm(h) + old_radius / chi + hc);
    if hc == 0.0 && viscous_stiffness == 0.0 {
        // A reversible branch has no finite pinning region, even at H=0 or on
        // the saturation sphere. Its transverse derivative must remain active.
        let magnitude = norm(h);
        let saturated = magnitude > 0.0 && chi * magnitude >= zs;
        let scale = if saturated { zs / magnitude } else { chi };
        let z = h.map(|v| scale * v);
        let tangent = symmetric_matrix(|i, j| {
            scale
                * ((if i == j { 1.0 } else { 0.0 })
                    - if saturated {
                        (h[i] / magnitude) * (h[j] / magnitude)
                    } else {
                        0.0
                    })
        });
        let multiplier = if saturated { magnitude / zs - a0 } else { 0.0 };
        let residual = norm(std::array::from_fn(|i| {
            h[i] - z[i] / chi - multiplier * z[i]
        }));
        if !residual.is_finite() || residual > target {
            return Err(invalid(
                "reversible vector branch missed its constitutive residual",
            ));
        }
        return Ok(BranchResponse {
            z,
            tangent,
            multiplier_m_h: multiplier,
            residual_a_m: residual,
            iterations: 0,
        });
    }
    let drive = std::array::from_fn(|i| h[i] - old[i] / chi);
    let mut pinned_remainder = drive;
    let mut normal = 0.0;
    // At the spherical constraint a normal reaction can balance outward H.
    // This is the actual pinned normal-cone condition, not component clipping.
    if old_radius >= zs * (1.0 - 8.0 * f64::EPSILON) && old_radius > 0.0 {
        let unit = old.map(|x| x / old_radius);
        normal = dot(drive, unit).max(0.0);
        pinned_remainder = std::array::from_fn(|i| drive[i] - normal * unit[i]);
    }
    let pinned_residual = (norm(pinned_remainder) - hc).max(0.0);
    if hc > 0.0 && pinned_residual <= target {
        // At a numerical pinning boundary, an unrepresentable displacement
        // cannot supply a reliable direction. Admit the retained state only
        // against the user's unchanged, directly measured field residual.
        return Ok(BranchResponse {
            z: old,
            tangent: [[0.0; 3]; 3],
            multiplier_m_h: if normal == 0.0 {
                0.0
            } else {
                normal / old_radius
            },
            residual_a_m: pinned_residual,
            iterations: 0,
        });
    }
    let normal_origin = if viscous_stiffness == 0.0 { a0 } else { 0.0 };
    let (free, free_tangent) =
        branch_unconstrained(old, h, a0, viscous_stiffness, normal_origin, hc)?;
    let (mut z, mut tangent, mut lambda, mut iterations) = (free, free_tangent, 0.0, 0);
    if norm(free) > zs {
        let mut lo = normal_origin;
        // At an active sphere z.(z-z_old)>=0 and the pinning normal work is
        // nonnegative. Thus lambda <= |H|/zs; this bracket needs no k*z_old.
        let mut hi = finite((norm(h) + hc) / zs)?.max(normal_origin);
        let mut a = hi;
        let mut accepted = false;
        for iteration in 0..c.maximum_branch_iterations {
            poll(cancelled)?;
            iterations = iteration + 1;
            let (point, d) = branch_unconstrained(old, h, a0, viscous_stiffness, a, hc)?;
            let radius = finite(norm(point))?;
            let defect = radius - zs;
            // Radius error converted through the local inverse scale to A/m.
            let stiffness = if viscous_stiffness == 0.0 {
                a
            } else {
                finite(a0 + viscous_stiffness + a)?
            };
            if defect.abs() * stiffness <= target && radius > 0.0 {
                let normalized = point.map(|x| x * (zs / radius));
                let candidate_lambda = if viscous_stiffness == 0.0 { a - a0 } else { a };
                // Radial normalization also rotates the pinning direction.
                // Its inverse scale hc/|delta z| can exceed a greatly, so a
                // small radius defect alone does not establish stationarity.
                // Keep refining this SAME multiplier/bracket until the actual
                // retained point passes the original field residual as well.
                if stationarity_residual(
                    old,
                    h,
                    normalized,
                    chi,
                    hc,
                    candidate_lambda,
                    viscous_stiffness,
                )? <= target
                {
                    z = normalized;
                    let dz = mv(d, z);
                    let denominator = finite(dot(z, dz))?;
                    if denominator <= 0.0 {
                        return Err(invalid(
                            "saturated vector tangent lost its positive radial compliance",
                        ));
                    }
                    tangent = symmetric_matrix(|i, j| d[i][j] - dz[i] * dz[j] / denominator);
                    lambda = candidate_lambda;
                    accepted = true;
                    break;
                }
            }
            if defect > 0.0 {
                lo = a;
            } else {
                hi = a;
            }
            let radial = if radius > 0.0 {
                dot(point, mv(d, point)) / radius
            } else {
                0.0
            };
            let proposed = if radial > 0.0 {
                a + defect / radial
            } else {
                f64::NAN
            };
            let next = if proposed.is_finite() && proposed > lo && proposed < hi {
                proposed
            } else {
                f64::midpoint(lo, hi)
            };
            if next == a {
                break;
            }
            a = next;
        }
        if !accepted {
            return Err(VectorHysteresisError::BranchNonconvergence {
                branch: index,
                iterations,
            });
        }
    }
    let mut residual = stationarity_residual(old, h, z, chi, hc, lambda, viscous_stiffness)?;
    if lambda == 0.0 && residual.is_finite() && residual > target && norm(z) <= zs {
        // Adding a tiny rotating slip to a large remanence quantizes its biased
        // component. The actual retained difference can then have a direction
        // error amplified by hc/|delta z|, despite exact analytic shrinkage.
        // Correct that ACTUAL KKT residual, never an ideal unretained increment.
        // Coordinate Newton steps let a well-resolved component absorb this
        // rounding without requiring an unrepresentable change in another.
        for _ in iterations..c.maximum_branch_iterations {
            poll(cancelled)?;
            iterations += 1;
            let increment = sub(z, old);
            let distance = finite(norm(increment))?;
            if distance == 0.0 {
                break;
            }
            let direction = increment.map(|value| value / distance);
            let base = finite(a0 + viscous_stiffness)?;
            let mut best = (z, residual);
            for axis in 0..3 {
                let dry = h[axis] - z[axis] / chi - hc * increment[axis] / distance;
                let force = if viscous_stiffness == 0.0 {
                    dry
                } else {
                    (-viscous_stiffness).mul_add(increment[axis], dry)
                };
                let curvature = base + hc * (1.0 - direction[axis] * direction[axis]) / distance;
                if !curvature.is_finite() || curvature <= 0.0 {
                    continue;
                }
                let mut trial = z;
                trial[axis] += force / curvature;
                if !trial[axis].is_finite() || trial[axis] == z[axis] || norm(trial) > zs {
                    continue;
                }
                let defect = stationarity_residual(old, h, trial, chi, hc, 0.0, viscous_stiffness)?;
                if defect < best.1 {
                    best = (trial, defect);
                }
            }
            if best.1 >= residual {
                break;
            }
            (z, residual) = best;
            if residual <= target {
                // Invert the original local Hessian using the final retained
                // slip. Its radial and transverse compliances are analytic;
                // this retains the same consistent field-Newton derivative.
                let increment = sub(z, old);
                let distance = finite(norm(increment))?;
                let radial = base.recip();
                tangent = if distance == 0.0 {
                    symmetric_matrix(|i, j| if hc == 0.0 && i == j { radial } else { 0.0 })
                } else {
                    let direction = increment.map(|value| value / distance);
                    let transverse = if hc == 0.0 {
                        radial
                    } else {
                        distance / finite(base * distance + hc)?
                    };
                    symmetric_matrix(|i, j| {
                        (if i == j { transverse } else { 0.0 })
                            + (radial - transverse) * direction[i] * direction[j]
                    })
                };
                break;
            }
        }
    }
    if !residual.is_finite() || residual > target || norm(z) > zs * (1.0 + 8.0 * f64::EPSILON) {
        return Err(invalid(
            "vector branch failed its original stationarity or saturation constraint",
        ));
    }
    Ok(BranchResponse {
        z,
        tangent,
        multiplier_m_h: lambda,
        residual_a_m: residual,
        iterations,
    })
}

struct Response {
    b: Vector,
    tangent: Matrix,
    maximum_residual: f64,
    branch_iterations: usize,
}

impl VectorHysteresis {
    /// Original positive reversible background [H/m], excluding all branch
    /// slopes. This is a coercivity lower bound, not the full permeability.
    pub fn background_permeability_h_m(&self) -> f64 {
        self.parameters.background_permeability_h_m
    }

    /// Original rotationally invariant field-magnitude support.
    pub(crate) fn field_limit_a_m(&self) -> f64 {
        self.parameters.field_limit_a_m
    }

    /// Reuse the scalar owner's parameter meanings and at-most-1024-branch
    /// admission. All three spatial components share each branch and norm.
    pub fn new(
        background_permeability_h_m: f64,
        field_limit_a_m: f64,
        branches: Vec<HysteresisBranch>,
    ) -> Result<Self> {
        let law = ScalarHysteresis::new(background_permeability_h_m, field_limit_a_m, branches)?;
        Self::from_isotropic_scalar(&law)
    }

    // The source adapter explicitly selects isotropy. Sharing the original
    // immutable parameters preserves exact coefficient and state identity.
    pub(crate) fn from_isotropic_scalar(law: &ScalarHysteresis) -> Result<Self> {
        if law
            .parameters
            .branches
            .iter()
            .any(|b| !b.susceptibility_h_m.recip().is_finite())
        {
            return Err(invalid(
                "vector branch inverse susceptibility is unrepresentable",
            ));
        }
        Ok(Self {
            parameters: Arc::clone(&law.parameters),
            viscosities_m_s_h: None,
        })
    }

    /// Bind one intrinsic polarization viscosity [m s/H] to each original
    /// branch. All values must be finite and nonnegative. A zero profile is
    /// exactly the original rate-independent law, without a new numerical path.
    ///
    /// A positive profile requires explicit durations for every candidate. Its
    /// fitted kinetics must exclude losses already owned by an Ohmic eddy model.
    /// Changing this profile does not convert an existing accepted history.
    pub fn with_viscosities_m_s_h(mut self, viscosities_m_s_h: Vec<f64>) -> Result<Self> {
        if viscosities_m_s_h.len() != self.parameters.branches.len()
            || viscosities_m_s_h
                .iter()
                .any(|eta| !eta.is_finite() || *eta < 0.0)
        {
            return Err(invalid(
                "one finite nonnegative intrinsic viscosity per original branch required",
            ));
        }
        self.viscosities_m_s_h = if viscosities_m_s_h.iter().any(|eta| *eta > 0.0) {
            Some(viscosities_m_s_h.into())
        } else {
            None
        };
        Ok(self)
    }

    /// Original immutable intrinsic kinetic fit, or None for exact zero rate
    /// dependence. Values have units [m s/H] and retain original branch order.
    pub fn viscosities_m_s_h(&self) -> Option<&[f64]> {
        self.viscosities_m_s_h.as_deref()
    }

    pub const fn is_rate_dependent(&self) -> bool {
        self.viscosities_m_s_h.is_some()
    }

    pub fn initial_state(&self) -> VectorHysteresisState {
        VectorHysteresisState {
            parameters: Arc::clone(&self.parameters),
            viscosities_m_s_h: self.viscosities_m_s_h.clone(),
            field_a_m: [0.0; 3],
            flux_density_t: [0.0; 3],
            polarizations_t: vec![[0.0; 3]; self.parameters.branches.len()],
            stored_energy_j_m3: 0.0,
        }
    }
    fn admit_history(&self, previous: &VectorHysteresisState) -> Result<()> {
        if !Arc::ptr_eq(&self.parameters, &previous.parameters)
            && self.parameters != previous.parameters
            || self.viscosities_m_s_h != previous.viscosities_m_s_h
        {
            return Err(invalid(
                "vector history must retain its original constitutive parameters",
            ));
        }
        Ok(())
    }
    fn admit_duration(&self, duration_s: Option<f64>) -> Result<()> {
        if let Some(duration) = duration_s {
            if !duration.is_finite() || duration <= 0.0 {
                return Err(invalid("positive finite switching duration required"));
            }
        } else if self.is_rate_dependent() {
            return Err(invalid(
                "rate-dependent vector history requires an explicit switching duration",
            ));
        }
        Ok(())
    }

    fn viscous_stiffness(&self, index: usize, duration_s: Option<f64>) -> Result<f64> {
        let Some(viscosities) = &self.viscosities_m_s_h else {
            return Ok(0.0);
        };
        let duration = duration_s
            .ok_or_else(|| invalid("rate-dependent vector history requires a duration"))?;
        let eta = viscosities[index];
        let stiffness = finite(eta / duration)?;
        if eta > 0.0 && stiffness == 0.0 {
            return Err(invalid("intrinsic viscous stiffness is unrepresentable"));
        }
        Ok(stiffness)
    }

    fn admit(
        &self,
        previous: &VectorHysteresisState,
        c: VectorHysteresisControls,
        duration_s: Option<f64>,
    ) -> Result<()> {
        admit_controls(c)?;
        self.admit_duration(duration_s)?;
        self.admit_history(previous)
    }

    /// Stable difference of the frozen-history H-controlled coenergy
    /// `mu_b |H|²/2 + sum(H.z - |z|²/(2chi) - hc|z-z_old|)`.
    /// `from` and `to` are candidate states evaluated from the SAME `previous`;
    /// this function never rebases their pinning norm to a Newton trial.
    /// The field solver owns candidate stationarity. This is a variational
    /// objective difference, not a physical work or heat increment.
    pub(crate) fn coenergy_difference_j_m3(
        &self,
        previous: &VectorHysteresisState,
        from: &VectorHysteresisState,
        to: &VectorHysteresisState,
    ) -> Result<f64> {
        self.admit_duration(None)?;
        self.admit_history(previous)?;
        self.admit_history(from)?;
        self.admit_history(to)?;
        let dh = sub(to.field_a_m, from.field_a_m);
        let hm = std::array::from_fn(|i| f64::midpoint(to.field_a_m[i], from.field_a_m[i]));
        let mut change = finite(self.parameters.background_permeability_h_m * dot(dh, hm))?;
        let mut correction = 0.0;
        for (index, branch) in self.parameters.branches.iter().enumerate() {
            let a = from.polarizations_t[index];
            let b = to.polarizations_t[index];
            let old = previous.polarizations_t[index];
            let dz = sub(b, a);
            let zm = std::array::from_fn(|i| f64::midpoint(a[i], b[i]));
            let da = sub(a, old);
            let db = sub(b, old);
            let na = finite(norm(da))?;
            let nb = finite(norm(db))?;
            let scale = na.max(nb);
            // Factored norm difference; scaling avoids overflow in na+nb and
            // preserves tiny changes beneath a large common pinning distance.
            let distance_change = if scale == 0.0 {
                0.0
            } else {
                let denominator = na / scale + nb / scale;
                dot(
                    dz,
                    std::array::from_fn(|i| (da[i] / scale + db[i] / scale) / denominator),
                )
            };
            let value = finite(
                dot(dh, zm)
                    + dot(
                        dz,
                        std::array::from_fn(|i| hm[i] - zm[i] / branch.susceptibility_h_m),
                    )
                    - branch.pinning_field_a_m * distance_change,
            )?;
            let next = finite(change + value)?;
            correction = finite(
                correction
                    + if change.abs() >= value.abs() {
                        (change - next) + value
                    } else {
                        (value - next) + change
                    },
            )?;
            change = next;
        }
        finite(change + correction)
    }
    fn response(
        &self,
        previous: &VectorHysteresisState,
        h: Vector,
        duration_s: Option<f64>,
        c: VectorHysteresisControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Response> {
        self.response_with(previous, h, duration_s, c, cancelled, &mut |_, _, _, _| {})
    }
    fn response_with(
        &self,
        previous: &VectorHysteresisState,
        h: Vector,
        duration_s: Option<f64>,
        c: VectorHysteresisControls,
        cancelled: &mut impl FnMut() -> bool,
        observe: &mut impl FnMut(&HysteresisBranch, Vector, f64, &BranchResponse),
    ) -> Result<Response> {
        poll(cancelled)?;
        let magnitude = norm(h);
        if !magnitude.is_finite() || magnitude > self.parameters.field_limit_a_m {
            return Err(invalid(
                "vector magnetic field exceeds its original magnitude domain",
            ));
        }
        let mu = self.parameters.background_permeability_h_m;
        let mut result = Response {
            b: h.map(|x| mu * x),
            tangent: std::array::from_fn(|i| {
                std::array::from_fn(|j| if i == j { mu } else { 0.0 })
            }),
            maximum_residual: 0.0,
            branch_iterations: 0,
        };
        for (index, (branch, &old)) in self
            .parameters
            .branches
            .iter()
            .zip(&previous.polarizations_t)
            .enumerate()
        {
            let stiffness = self.viscous_stiffness(index, duration_s)?;
            let r = branch_response(branch, old, h, stiffness, c, index, cancelled)?;
            for i in 0..3 {
                result.b[i] = finite(result.b[i] + r.z[i])?;
                for j in 0..3 {
                    result.tangent[i][j] = finite(result.tangent[i][j] + r.tangent[i][j])?;
                }
            }
            result.maximum_residual = result.maximum_residual.max(r.residual_a_m);
            result.branch_iterations += r.iterations;
            observe(branch, old, stiffness, &r);
        }
        Ok(result)
    }
    /// One backward variational increment from the SAME accepted history.
    /// Refine changes of field direction to resolve the continuous rotating path.
    pub fn candidate_field(
        &self,
        previous: &VectorHysteresisState,
        field_a_m: Vector,
        controls: VectorHysteresisControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<VectorHysteresisStep> {
        self.candidate_field_increment(previous, field_a_m, None, controls, cancelled)
    }

    /// Backward Euler switching over one positive physical duration [s]. Every
    /// candidate starts from the same accepted memory. Holds can continue to
    /// relax while the applied field is unchanged; time refinement is required
    /// to resolve the continuous kinetics. Original storage and physical dry,
    /// viscous, and numerical transfers remain separately identifiable.
    pub fn candidate_field_timed(
        &self,
        previous: &VectorHysteresisState,
        field_a_m: Vector,
        duration_s: f64,
        controls: VectorHysteresisControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<VectorHysteresisStep> {
        self.candidate_field_increment(previous, field_a_m, Some(duration_s), controls, cancelled)
    }

    fn candidate_field_increment(
        &self,
        previous: &VectorHysteresisState,
        field_a_m: Vector,
        duration_s: Option<f64>,
        controls: VectorHysteresisControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<VectorHysteresisStep> {
        self.admit(previous, controls, duration_s)?;
        let mu = self.parameters.background_permeability_h_m;
        let dh = sub(field_a_m, previous.field_a_m);
        let mut energy = 0.5 * mu * dot(field_a_m, field_a_m);
        let mut change = mu
            * dot(
                dh,
                std::array::from_fn(|i| f64::midpoint(field_a_m[i], previous.field_a_m[i])),
            );
        let mut dissipation = 0.5 * mu * dot(dh, dh);
        let mut pinning_heat = 0.0;
        let mut viscous_heat = 0.0;
        let mut polarizations = Vec::with_capacity(self.parameters.branches.len());
        // Each branch return map is evaluated once. Inverse trials use the
        // same kernel without allocating or retaining candidate branch memory.
        let response = self.response_with(
            previous,
            field_a_m,
            duration_s,
            controls,
            cancelled,
            &mut |branch, old, stiffness, r| {
                let dz = sub(r.z, old);
                let radius = norm(old);
                let zs = branch.saturation_polarization_t;
                // The radial term is nonnegative by the original ball constraint.
                // A <=8-ulp old-radius excess is only normalization roundoff; full
                // endpoint work is checked independently below.
                let radial = (zs - radius).max(0.0) * (zs + radius);
                let normal_work = 0.5 * r.multiplier_m_h * (dot(dz, dz) + radial);
                energy += 0.5 * dot(r.z, r.z) / branch.susceptibility_h_m;
                change += dot(dz, std::array::from_fn(|i| f64::midpoint(r.z[i], old[i])))
                    / branch.susceptibility_h_m;
                pinning_heat += branch.pinning_field_a_m * norm(dz);
                // Form force before work to avoid squaring a tiny increment
                // before multiplying by a large, still representable stiffness.
                viscous_heat += dot(dz, dz.map(|value| stiffness * value));
                dissipation += 0.5 * dot(dz, dz) / branch.susceptibility_h_m + normal_work;
                polarizations.push(r.z);
            },
        )?;
        let heat = pinning_heat + viscous_heat;
        let work = dot(field_a_m, sub(response.b, previous.flux_density_t));
        let residual = work - change - heat - dissipation;
        let allowance = controls
            .work_tolerance_j_m3
            .max(controls.relative_tolerance * (work.abs() + change.abs() + heat + dissipation));
        for x in [
            energy,
            change,
            pinning_heat,
            viscous_heat,
            heat,
            dissipation,
            work,
            residual,
        ] {
            finite(x)?;
        }
        if heat < 0.0 || dissipation < 0.0 || residual.abs() > allowance {
            return Err(invalid(
                "vector endpoint work, storage, physical heat and numerical dissipation do not close",
            ));
        }
        poll(cancelled)?;
        Ok(VectorHysteresisStep {
            state: VectorHysteresisState {
                parameters: Arc::clone(&self.parameters),
                viscosities_m_s_h: self.viscosities_m_s_h.clone(),
                field_a_m,
                flux_density_t: response.b,
                polarizations_t: polarizations,
                stored_energy_j_m3: energy,
            },
            irreversible_heat_j_m3: heat,
            pinning_heat_j_m3: pinning_heat,
            viscous_heat_j_m3: viscous_heat,
            energy_change_j_m3: change,
            endpoint_work_j_m3: work,
            numerical_dissipation_j_m3: dissipation,
            energy_residual_j_m3: residual,
            incremental_permeability_h_m: response.tangent,
            maximum_constitutive_residual_a_m: response.maximum_residual,
            inverse_iterations: 0,
            branch_iterations: response.branch_iterations,
            inverse_residual_t: 0.0,
        })
    }

    /// B-driven inverse using the consistent SPD frozen-history permeability.
    /// Every Newton/line-search trial starts from `previous`; no trial commits
    /// internal polarization. The returned B is recomputed, never overwritten
    /// with the requested value to hide a constitutive residual.
    pub fn candidate_flux_density(
        &self,
        previous: &VectorHysteresisState,
        flux_density_t: Vector,
        controls: VectorHysteresisControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<VectorHysteresisStep> {
        self.candidate_flux_increment(previous, flux_density_t, None, controls, cancelled)
    }

    /// Timed B-driven inverse of the same backward Euler return map. All Newton
    /// trials use this unchanged duration and original accepted polarization;
    /// the returned flux, tangent, work, and both physical heats are recomputed.
    pub fn candidate_flux_density_timed(
        &self,
        previous: &VectorHysteresisState,
        flux_density_t: Vector,
        duration_s: f64,
        controls: VectorHysteresisControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<VectorHysteresisStep> {
        self.candidate_flux_increment(
            previous,
            flux_density_t,
            Some(duration_s),
            controls,
            cancelled,
        )
    }

    fn candidate_flux_increment(
        &self,
        previous: &VectorHysteresisState,
        flux_density_t: Vector,
        duration_s: Option<f64>,
        controls: VectorHysteresisControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<VectorHysteresisStep> {
        self.admit(previous, controls, duration_s)?;
        if flux_density_t.iter().any(|x| !x.is_finite()) {
            return Err(invalid("finite vector flux density required"));
        }
        let magnitude = finite(norm(flux_density_t))?;
        let maximum = self.parameters.background_permeability_h_m * self.parameters.field_limit_a_m
            + self
                .parameters
                .branches
                .iter()
                .map(|b| b.saturation_polarization_t)
                .sum::<f64>();
        if magnitude > maximum {
            return Err(VectorHysteresisError::FluxRange {
                observed_t: magnitude,
                maximum_t: maximum,
            });
        }
        let target = finite(
            controls
                .absolute_flux_tolerance_t
                .max(controls.relative_tolerance * magnitude),
        )?;
        let mut h = previous.field_a_m;
        let mut response = self.response(previous, h, duration_s, controls, cancelled)?;
        let mut residual = finite(norm(sub(response.b, flux_density_t)))?;
        let mut branch_iterations = response.branch_iterations;
        for iteration in 0..=controls.maximum_inverse_iterations {
            poll(cancelled)?;
            if residual <= target {
                let mut result =
                    self.candidate_field_increment(previous, h, duration_s, controls, cancelled)?;
                result.inverse_iterations = iteration;
                result.branch_iterations += branch_iterations;
                result.inverse_residual_t =
                    finite(norm(sub(result.state.flux_density_t, flux_density_t)))?;
                if result.inverse_residual_t > target {
                    return Err(invalid(
                        "retained vector inverse missed its physical flux residual",
                    ));
                }
                return Ok(result);
            }
            if iteration == controls.maximum_inverse_iterations {
                break;
            }
            let flat: [f64; 9] = std::array::from_fn(|k| response.tangent[k / 3][k % 3]);
            let factor = cholesky(&flat, 3)
                .map_err(|_| invalid("vector incremental permeability is not representably SPD"))?;
            let mut direction = sub(response.b, flux_density_t);
            factor.solve(&mut direction);
            if direction.iter().any(|x| !x.is_finite()) {
                return Err(invalid(
                    "vector inverse Newton direction is unrepresentable",
                ));
            }
            let mut scale: f64 = 1.0;
            let mut accepted = None;
            for _ in 0..controls.maximum_backtracks {
                poll(cancelled)?;
                let trial = std::array::from_fn(|i| (-scale).mul_add(direction[i], h[i]));
                if norm(trial) <= self.parameters.field_limit_a_m {
                    let next = self.response(previous, trial, duration_s, controls, cancelled)?;
                    branch_iterations += next.branch_iterations;
                    let error = finite(norm(sub(next.b, flux_density_t)))?;
                    if error <= target || error <= (1.0 - 1e-4 * scale) * residual {
                        accepted = Some((trial, next, error));
                        break;
                    }
                }
                scale *= 0.5;
            }
            let Some((next_h, next_response, next_residual)) = accepted else {
                return Err(VectorHysteresisError::FluxNonconvergence {
                    iterations: iteration + 1,
                    residual_t: residual,
                });
            };
            h = next_h;
            response = next_response;
            residual = next_residual;
        }
        Err(VectorHysteresisError::FluxNonconvergence {
            iterations: controls.maximum_inverse_iterations,
            residual_t: residual,
        })
    }
}

#[cfg(test)]
mod tests;
