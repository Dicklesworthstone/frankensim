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
    BranchNonconvergence { branch: usize, iterations: usize },
    FluxNonconvergence { iterations: usize, residual_t: f64 },
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
}

/// Private accepted three-dimensional polarization memory and its physical state.
#[derive(Clone, Debug, PartialEq)]
pub struct VectorHysteresisState {
    parameters: Arc<Parameters>,
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
    /// Nonnegative incremental dry-friction heat; deposit at most once [J/m³].
    pub irreversible_heat_j_m3: f64,
    /// Stable difference of the complete stored energy [J/m³].
    pub energy_change_j_m3: f64,
    /// H_new dot the actually retained change of B [J/m³].
    pub endpoint_work_j_m3: f64,
    /// Nonnegative endpoint-discretization remainder, never physical heat [J/m³].
    pub numerical_dissipation_j_m3: f64,
    /// Endpoint work minus stored change, pinning heat and numerical remainder.
    pub energy_residual_j_m3: f64,
    /// Frozen-history dB/dH [H/m], symmetric positive definite. At pinning and
    /// saturation kinks the selected active set is not a smoothness claim.
    pub incremental_permeability_h_m: Matrix,
    pub maximum_constitutive_residual_a_m: f64,
    pub inverse_iterations: usize,
    pub branch_iterations: usize,
    pub inverse_residual_t: f64,
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
    let length = finite(norm(r))?;
    if length <= hc {
        return Ok((old, [[0.0; 3]; 3]));
    }
    let direction = r.map(|x| x / length);
    let fraction = hc / length;
    let z = std::array::from_fn(|i| ((1.0 - fraction) / a).mul_add(r[i], old[i]));
    let tangent = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            ((if i == j { 1.0 - fraction } else { 0.0 }) + fraction * direction[i] * direction[j])
                / a
        })
    });
    for value in z.into_iter().chain(tangent.into_iter().flatten()) {
        finite(value)?;
    }
    Ok((z, tangent))
}

fn branch_response(
    branch: &HysteresisBranch,
    old: Vector,
    h: Vector,
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
    if hc == 0.0 {
        // A reversible branch has no finite pinning region, even at H=0 or on
        // the saturation sphere. Its transverse derivative must remain active.
        let magnitude = norm(h);
        let saturated = magnitude > 0.0 && chi * magnitude >= zs;
        let scale = if saturated { zs / magnitude } else { chi };
        let z = h.map(|v| scale * v);
        let tangent = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                scale
                    * ((if i == j { 1.0 } else { 0.0 })
                        - if saturated {
                            (h[i] / magnitude) * (h[j] / magnitude)
                        } else {
                            0.0
                        })
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
    if pinned_residual <= target {
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
    let (free, free_tangent) = unconstrained(old, h, a0, hc)?;
    let (mut z, mut tangent, mut lambda, mut iterations) = (free, free_tangent, 0.0, 0);
    if norm(free) > zs {
        let mut lo = a0;
        let mut hi = finite((norm(h) + hc) / zs)?.max(a0);
        let mut a = hi;
        let mut accepted = false;
        for iteration in 0..c.maximum_branch_iterations {
            poll(cancelled)?;
            iterations = iteration + 1;
            let (point, d) = unconstrained(old, h, a, hc)?;
            let radius = finite(norm(point))?;
            let defect = radius - zs;
            // Radius error converted through the local inverse scale to A/m.
            if defect.abs() * a <= target && radius > 0.0 {
                // Only a converged multiplier is normalized for roundoff.
                z = point.map(|x| x * (zs / radius));
                let (_, d) = unconstrained(old, h, a, hc)?;
                let dz = mv(d, z);
                let denominator = finite(dot(z, dz))?;
                if denominator <= 0.0 {
                    return Err(invalid(
                        "saturated vector tangent lost its positive radial compliance",
                    ));
                }
                tangent = std::array::from_fn(|i| {
                    std::array::from_fn(|j| d[i][j] - dz[i] * dz[j] / denominator)
                });
                lambda = a - a0;
                accepted = true;
                break;
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
    let increment = sub(z, old);
    let distance = finite(norm(increment))?;
    let residual = if distance == 0.0 {
        (norm(std::array::from_fn(|i| h[i] - z[i] / chi - lambda * z[i])) - hc).max(0.0)
    } else {
        norm(std::array::from_fn(|i| {
            h[i] - z[i] / chi - lambda * z[i] - hc * increment[i] / distance
        }))
    };
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
    /// Reuse the scalar owner's parameter meanings and at-most-1024-branch
    /// admission. All three spatial components share each branch and norm.
    pub fn new(
        background_permeability_h_m: f64,
        field_limit_a_m: f64,
        branches: Vec<HysteresisBranch>,
    ) -> Result<Self> {
        let law = ScalarHysteresis::new(background_permeability_h_m, field_limit_a_m, branches)?;
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
            parameters: law.parameters,
        })
    }
    pub fn initial_state(&self) -> VectorHysteresisState {
        VectorHysteresisState {
            parameters: Arc::clone(&self.parameters),
            field_a_m: [0.0; 3],
            flux_density_t: [0.0; 3],
            polarizations_t: vec![[0.0; 3]; self.parameters.branches.len()],
            stored_energy_j_m3: 0.0,
        }
    }
    fn admit(&self, previous: &VectorHysteresisState, c: VectorHysteresisControls) -> Result<()> {
        admit_controls(c)?;
        if !Arc::ptr_eq(&self.parameters, &previous.parameters)
            && self.parameters != previous.parameters
        {
            return Err(invalid(
                "vector history must retain its original constitutive parameters",
            ));
        }
        Ok(())
    }
    fn response(
        &self,
        previous: &VectorHysteresisState,
        h: Vector,
        c: VectorHysteresisControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Response> {
        self.response_with(previous, h, c, cancelled, &mut |_, _, _| {})
    }
    fn response_with(
        &self,
        previous: &VectorHysteresisState,
        h: Vector,
        c: VectorHysteresisControls,
        cancelled: &mut impl FnMut() -> bool,
        observe: &mut impl FnMut(&HysteresisBranch, Vector, &BranchResponse),
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
            let r = branch_response(branch, old, h, c, index, cancelled)?;
            for i in 0..3 {
                result.b[i] = finite(result.b[i] + r.z[i])?;
                for j in 0..3 {
                    result.tangent[i][j] = finite(result.tangent[i][j] + r.tangent[i][j])?;
                }
            }
            result.maximum_residual = result.maximum_residual.max(r.residual_a_m);
            result.branch_iterations += r.iterations;
            observe(branch, old, &r);
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
        self.admit(previous, controls)?;
        let mu = self.parameters.background_permeability_h_m;
        let dh = sub(field_a_m, previous.field_a_m);
        let mut energy = 0.5 * mu * dot(field_a_m, field_a_m);
        let mut change = mu
            * dot(
                dh,
                std::array::from_fn(|i| f64::midpoint(field_a_m[i], previous.field_a_m[i])),
            );
        let mut dissipation = 0.5 * mu * dot(dh, dh);
        let mut heat = 0.0;
        let mut polarizations = Vec::with_capacity(self.parameters.branches.len());
        // Each branch return map is evaluated once. Inverse trials use the
        // same kernel without allocating or retaining candidate branch memory.
        let response = self.response_with(
            previous,
            field_a_m,
            controls,
            cancelled,
            &mut |branch, old, r| {
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
                heat += branch.pinning_field_a_m * norm(dz);
                dissipation += 0.5 * dot(dz, dz) / branch.susceptibility_h_m + normal_work;
                polarizations.push(r.z);
            },
        )?;
        let work = dot(field_a_m, sub(response.b, previous.flux_density_t));
        let residual = work - change - heat - dissipation;
        let allowance = controls
            .work_tolerance_j_m3
            .max(controls.relative_tolerance * (work.abs() + change.abs() + heat + dissipation));
        for x in [energy, change, heat, dissipation, work, residual] {
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
                field_a_m,
                flux_density_t: response.b,
                polarizations_t: polarizations,
                stored_energy_j_m3: energy,
            },
            irreversible_heat_j_m3: heat,
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
        self.admit(previous, controls)?;
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
            return Err(invalid(
                "vector flux requires a field outside the original model domain",
            ));
        }
        let target = finite(
            controls
                .absolute_flux_tolerance_t
                .max(controls.relative_tolerance * magnitude),
        )?;
        let mut h = previous.field_a_m;
        let mut response = self.response(previous, h, controls, cancelled)?;
        let mut residual = finite(norm(sub(response.b, flux_density_t)))?;
        let mut branch_iterations = response.branch_iterations;
        for iteration in 0..=controls.maximum_inverse_iterations {
            poll(cancelled)?;
            if residual <= target {
                let mut result = self.candidate_field(previous, h, controls, cancelled)?;
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
                    let next = self.response(previous, trial, controls, cancelled)?;
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
