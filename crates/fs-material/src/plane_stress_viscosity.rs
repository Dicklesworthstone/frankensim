//! Source-bound, reduced plane-stress Kelvin--Voigt viscosity.
//!
//! In the declared material plane, [sigma_11, sigma_22, sigma_12] =
//! V [epsilon_dot_11, epsilon_dot_22, gamma_dot_12]. V is symmetric,
//! in Pa s, and positive semidefinite, including directional loss and zero. These are
//! REDUCED in-plane coefficients, not three-dimensional bulk viscosities:
//! independently condensing an arbitrary 3-D elastic/viscous law is generally
//! frequency dependent and is not performed here. No elastic Poisson ratio,
//! isotropy, loss factor, material name or damping time is used to invent V.
//!
//! The engineering-shear convention makes stress dot strain-rate the physical
//! dissipation density. Rotation uses the original tensor owner. A homogeneous
//! plate's bending viscosity is V*h^3/12; this does not add membrane prestress
//! loss, a heat-transport law, hereditary memory or a calibrated material claim.

use crate::state_point::ResolvedMaterialStatePoint;
use crate::{MaterialError, tensor};
use fs_qty::{DynViscosity, QuantitySpec};

/// Six independent reduced coefficients in order 11,22,66,12,16,26.
/// Zero couplings must be explicitly supplied; missing entries never mean zero.
pub const PLANE_STRESS_VISCOSITY_PROPERTIES: [&str; 6] = [
    "plane_stress_viscosity_11", "plane_stress_viscosity_22",
    "plane_stress_viscosity_66", "plane_stress_viscosity_12",
    "plane_stress_viscosity_16", "plane_stress_viscosity_26",
];

/// Material axes spanning the characterized plane, in ordered local x,y roles.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MaterialPlane {
    /// Material axes 1,2; axis 3 is normal.
    Axes12,
    /// Material axes 2,1 in that order; coefficients must use that ordered basis.
    Axes21,
    /// Material axes 1,3; axis 2 is normal.
    Axes13,
    /// Material axes 2,3; axis 1 is normal.
    Axes23,
}

/// Frozen source state and its explicitly selected reduced dissipative tensor.
/// No independent material-query point, geometry or mutable history is hidden.
#[derive(Clone, Debug, PartialEq)]
pub struct PlaneStressViscosity {
    source: ResolvedMaterialStatePoint,
    properties: [String; 6],
    plane: MaterialPlane,
    matrix: [f64; 9],
    // Only semidefinite sources need a Gram evaluation. Keeping the established
    // SPD path preserves its stress, power and rotation arithmetic exactly.
    gram: Option<[[f64; 6]; 3]>,
    omega_band: (f64, f64),
}

fn invalid(what: &str) -> MaterialError {
    MaterialError::Parameters { what: what.into() }
}
/// Factor a passive symmetric three-component work-conjugate resultant as
/// `V = L L^T`. This is bilinear-form admission, not a three-dimensional
/// constitutive reduction. Zero eigenvalues and directional loss are allowed.
///
/// Uses the shared original-unit compensated work-matrix factor below. The
/// `64*EPSILON` admission allowance is relative to each pair of diagonal row
/// scales; no positive loss floor or rounded normalized source is introduced.
/// # Errors
/// Nonfinite, asymmetric or active data, or an unrepresentable factor.
pub fn factor_passive_resultant(v: &[f64; 9]) -> Result<[[f64; 6]; 3], MaterialError> {
    if v.iter().any(|x| !x.is_finite())
        || (0..3).any(|i| (0..3).any(|j| v[3*i+j] != v[3*j+i])) {
        return Err(invalid("passive resultant must be finite and symmetric"));
    }
    for i in 0..3 {
        if v[3*i+i] < 0.0
            || (v[3*i+i] == 0.0 && (0..3).any(|j| v[3*i+j] != 0.0)) {
            return Err(invalid("passive resultant has a negative diagonal or coupled zero diagonal"));
        }
    }
    let mut block = [[0.0; 6]; 6];
    for i in 0..3 { for j in 0..3 { block[i][j] = v[3*i+j]; } }
    let factor = factor_passive_work_matrix(&block)?;
    Ok(core::array::from_fn(|i| factor[i]))
}

/// Factor a finite symmetric passive six-component work matrix as `L L^T`.
/// The supplied coordinates are work-conjugate: Mandel tensors and shell
/// membrane/curvature resultants may use different dimensional row scales.
/// Zero and rank-deficient inputs preserve their resolved null directions.
///
/// Pivoted LDL Schur updates use original dimensional coefficients and the
/// shared double-double owner. Square roots occur only in output channels.
/// Positive pivots are retained down to local introduced-arithmetic uncertainty,
/// independently of the largest diagonal. The final factor reconstructs each
/// source entry within `64*EPSILON` of its pair of diagonal row scales.
/// This matrix-only admission cannot identify error already present in input
/// data, and supplies no empirical material or discretization certificate.
///
/// Fixed six-component storage needs no heap allocation or eigensolve.
pub fn factor_passive_work_matrix(block: &[[f64; 6]; 6])
    -> Result<[[f64; 6]; 6], MaterialError>
{
    use fs_math::dd::Dd;
    let bad = || invalid("material work matrix is not finite, symmetric, representable and passive");
    if block.iter().flatten().any(|x| !x.is_finite())
        || (0..6).any(|i| (0..6).any(|j| block[i][j] != block[j][i])) {
        return Err(bad());
    }
    let mut scale = [0.0; 6];
    for i in 0..6 {
        if block[i][i] < 0.0
            || (block[i][i] == 0.0 && block[i].iter().any(|v| *v != 0.0)) {
            return Err(bad());
        }
        scale[i] = if block[i][i] == 0.0 { 1.0 } else { fs_math::det::sqrt(block[i][i]) };
    }
    const ADMISSION: f64 = 64.0 * f64::EPSILON;
    // Error starts at zero for the exact input f64 coefficients. An untouched
    // small diagonal therefore has no tolerance tied to a stronger direction.
    const ROUND: f64 = 16.0 * f64::EPSILON * f64::EPSILON;
    let mut residual = block.map(|row| row.map(Dd::from_f64));
    let mut uncertainty = [[0.0; 6]; 6];
    let mut used = [false; 6];
    let mut out = [[0.0; 6]; 6];
    for column in 0..6 {
        let mut pivot = None;
        let mut diagonal = Dd::ZERO;
        for i in 0..6 {
            if used[i] { continue; }
            let value = residual[i][i];
            let normalized = value.to_f64() / scale[i] / scale[i];
            if !normalized.is_finite() || normalized < -ADMISSION { return Err(bad()); }
            if diagonal.lt(value) && Dd::from_f64(uncertainty[i][i]).lt(value) {
                pivot = Some(i);
                diagonal = value;
            }
        }
        let Some(p) = pivot else { break; };
        used[p] = true;
        let lower = diagonal.to_f64() - uncertainty[p][p];
        if !lower.is_finite() || lower <= 0.0 { return Err(bad()); }
        let root = diagonal.sqrt();
        out[p][column] = root.to_f64();
        let mut ratios = [Dd::ZERO; 6];
        let mut ratio_error = [0.0; 6];
        for i in 0..6 {
            if used[i] { continue; }
            out[i][column] = (residual[i][p] / root).to_f64();
            ratios[i] = residual[i][p] / diagonal;
            let ratio = ratios[i].abs().to_f64();
            ratio_error[i] = (uncertainty[i][p] + ratio * uncertainty[p][p]) / lower
                + ROUND * ratio;
        }
        if out[p][column] == 0.0 || out.iter().any(|row| !row[column].is_finite()) {
            return Err(bad());
        }
        for i in 0..6 {
            if used[i] { continue; }
            for j in i..6 {
                if used[j] { continue; }
                // Dividing the larger numerator first retains representable
                // products even when the smaller individual ratio underflows.
                let (a, b) = if residual[i][p].abs().lt(residual[j][p].abs()) {
                    (j, i)
                } else { (i, j) };
                let update = ratios[a] * residual[b][p];
                let error = uncertainty[i][j]
                    + ratios[a].abs().to_f64() * uncertainty[b][p]
                    + ratio_error[a] * residual[b][p].abs().to_f64()
                    + ratio_error[a] * uncertainty[b][p]
                    + ROUND * residual[i][j].abs().to_f64()
                    + ROUND * update.abs().to_f64();
                let next = residual[i][j] - update;
                if !error.is_finite() || !next.hi.is_finite() || !next.lo.is_finite() {
                    return Err(bad());
                }
                residual[i][j] = next;
                residual[j][i] = next;
                uncertainty[i][j] = error;
                uncertainty[j][i] = error;
            }
        }
    }
    // An active cross term cannot disappear behind an unrelated strong axis.
    for i in 0..6 { for j in i..6 {
        let mut actual = Dd::ZERO;
        for column in 0..6 {
            actual = actual + Dd::from_f64(out[i][column]) * Dd::from_f64(out[j][column]);
        }
        let difference = (actual - Dd::from_f64(block[i][j])).abs().to_f64();
        let normalized = difference / scale[i] / scale[j];
        if !normalized.is_finite() || normalized > ADMISSION { return Err(bad()); }
    } }
    Ok(out)
}

fn admit_passive(v: &[f64; 9]) -> Result<Option<[[f64; 6]; 3]>, MaterialError> {
    if admit_matrix(v).is_ok() { Ok(None) }
    else { factor_passive_resultant(v).map(Some) }
}

fn admit_matrix(v: &[f64; 9]) -> Result<(), MaterialError> {
    if v.iter().any(|x| !x.is_finite()) {
        return Err(invalid("plane-stress viscosity must be finite"));
    }
    if v.iter().all(|x| *x == 0.0) { return Ok(()); }
    if [v[0], v[4], v[8]].iter().any(|x| *x <= 0.0) {
        return Err(invalid("viscosity requires positive definiteness or an explicit all-zero tensor"));
    }
    // Diagonal congruence avoids products of dimensional coefficients. Refuse
    // unresolved/indefinite pivots; do not clip eigenvalues or repair source data.
    let s = [v[0].sqrt(), v[4].sqrt(), v[8].sqrt()];
    let a = v[1] / s[0] / s[1];
    let b = v[2] / s[0] / s[2];
    let c = v[5] / s[1] / s[2];
    let pivot = 1.0 - a*a;
    if !pivot.is_finite() || pivot <= 0.0 {
        return Err(invalid("plane-stress viscosity has a nonpositive or unresolved pivot"));
    }
    let l = (c - a*b) / pivot.sqrt();
    let last = 1.0 - b*b - l*l;
    if !last.is_finite() || last <= 0.0 {
        return Err(invalid("plane-stress viscosity has a nonpositive or unresolved pivot"));
    }
    Ok(())
}

impl PlaneStressViscosity {
    /// Bind six source-resolved scalar coefficients at the SAME frozen state.
    /// Every selected property must have an explicit frequency-validity band;
    /// their intersection is retained even for zero entries. Only constants or
    /// temperature-only curves at a fixed physical temperature are admitted.
    /// Names are in 11,22,66,12,16,26 order and must be distinct/nonempty.
    ///
    /// Rank-deficient tensors retain their undamped directions; no weak loss
    /// is inserted. The original full-rank evaluation remains unchanged.
    /// Source uncertainty and provenance remain those of the original claims.
    pub fn from_material(source: ResolvedMaterialStatePoint, plane: MaterialPlane,
        properties: [&str; 6]) -> Result<Self, MaterialError>
    {
        let mut values = [0.0; 6];
        let mut band = (0.0_f64, f64::INFINITY);
        for (i, name) in properties.iter().enumerate() {
            if name.is_empty() || name.len() > 256 || properties[..i].contains(name) {
                return Err(invalid("viscosity needs six distinct bounded source-property names"));
            }
            let p = source.property(name).ok_or_else(|| invalid("missing plane-stress viscosity coefficient"))?;
            if p.requirement().quantity() != QuantitySpec::dimensional(DynViscosity::DIMS)
                || !p.is_constant_at_fixed_environment() || !p.value_si().is_finite() {
                return Err(invalid("reduced viscosity needs fixed-state constants with quantity Pa s"));
            }
            let (lo, hi) = p.angular_frequency_band_rad_s()
                .map_err(|e| invalid(&format!("viscosity frequency context: {e}")))?
                .ok_or_else(|| invalid("each viscosity coefficient requires an explicit frequency band"))?;
            band = (band.0.max(lo), band.1.min(hi));
            values[i] = p.value_si();
        }
        if !band.0.is_finite() || !band.1.is_finite() || band.0 < 0.0 || band.1 < band.0 {
            return Err(invalid("viscosity coefficient frequency bands do not intersect"));
        }
        let [a, b, c, d, e, f] = values;
        let matrix = [a, d, e, d, b, f, e, f, c];
        let gram = admit_passive(&matrix)?;
        Ok(Self { source, plane, properties: properties.map(str::to_string), matrix, gram, omega_band: band })
    }

    /// Exact original material bundle, including original quantity/source records.
    #[must_use] pub const fn source(&self) -> &ResolvedMaterialStatePoint { &self.source }
    /// The explicitly characterized in-plane material axes.
    #[must_use] pub const fn plane(&self) -> MaterialPlane { self.plane }
    /// Selected source property names in 11,22,66,12,16,26 order.
    #[must_use] pub const fn properties(&self) -> &[String; 6] { &self.properties }
    /// Dense symmetric material-frame V in engineering-shear convention [Pa s].
    #[must_use] pub const fn matrix_pa_s(&self) -> &[f64; 9] { &self.matrix }
    /// Intersected angular-frequency applicability [rad/s], not a low-pass filter.
    #[must_use] pub const fn omega_band_rad_s(&self) -> (f64, f64) { self.omega_band }
    /// Check an actual modal frequency, not merely the source query frequency.
    pub fn validate_frequency(&self, omega: f64) -> Result<(), MaterialError> {
        if !omega.is_finite() || omega < self.omega_band.0 || omega > self.omega_band.1 {
            return Err(invalid("actual frequency lies outside the reduced viscosity source band"));
        }
        Ok(())
    }

    /// Stress [Pa] and irreversible power density [W/m3] at a finite local
    /// engineering strain rate [1/s]. No internal history or extra stored energy.
    /// Unrepresentable stress/power refuses; negative roundoff is never clamped.
    pub fn stress_and_power(&self, strain_rate: [f64; 3]) -> Result<([f64; 3], f64), MaterialError> {
        if strain_rate.iter().any(|x| !x.is_finite()) {
            return Err(invalid("viscous strain rate must be finite"));
        }
        let (stress, power) = if let Some(factor) = &self.gram {
            let projected: [f64; 6] = core::array::from_fn(|k|
                (0..3).map(|i| factor[i][k]*strain_rate[i]).sum());
            let stress = core::array::from_fn(|i|
                (0..6).map(|k| factor[i][k]*projected[k]).sum());
            (stress, projected.iter().map(|x| x*x).sum())
        } else {
            let stress: [f64; 3] = core::array::from_fn(|i|
                (0..3).map(|j| self.matrix[3*i+j]*strain_rate[j]).sum());
            let power = stress.iter().zip(strain_rate).map(|(s,e)| s*e).sum::<f64>();
            (stress, power)
        };
        if stress.iter().any(|x| !x.is_finite()) || !power.is_finite() || power < 0.0 {
            return Err(invalid("viscous stress or dissipated power is unrepresentable"));
        }
        Ok((stress, power))
    }

    /// Tensor in mesh coordinates when material in-plane x is rotated CCW by
    /// angle about the positive plane normal. Geometry is not rotated. Reuse the
    /// original tensor transforms for engineering strain and tensor stress.
    pub fn rotated_matrix_pa_s(&self, angle: f64) -> Result<[f64; 9], MaterialError> {
        if !angle.is_finite() { return Err(invalid("viscosity orientation must be finite radians")); }
        if angle == 0.0 || self.matrix.iter().all(|x| *x == 0.0) { return Ok(self.matrix); }
        let q = tensor::rotation([0.0, 0.0, 1.0], angle);
        let qt = std::array::from_fn(|i| std::array::from_fn(|j| q[j][i]));
        let mut out = [0.0; 9];
        for j in 0..3 {
            let mut e = [0.0; 6];
            e[if j == 2 { 3 } else { j }] = if j == 2 { 0.5 } else { 1.0 };
            let local = tensor::rotate(&e, &qt);
            let (s, _) = self.stress_and_power([local[0],local[1],2.0*local[3]])?;
            let s = tensor::rotate(&[s[0],s[1],0.0,s[2],0.0,0.0], &q);
            let column = [s[0],s[1],s[3]];
            // One rounded value per symmetric pair preserves work conjugacy.
            for i in 0..=j { out[3*i+j] = column[i]; out[3*j+i] = column[i]; }
        }
        admit_passive(&out)?;
        Ok(out)
    }

    /// Homogeneous-thickness bending viscosity [N m s] in mesh coordinates.
    /// Both stiffness and viscosity must use this same engineering-curvature
    /// convention; no elastic stiffness or geometric prestress is damped here.
    pub fn bending_matrix_n_m_s(&self, thickness_m: f64, angle: f64) -> Result<[f64; 9], MaterialError> {
        if !thickness_m.is_finite() || thickness_m <= 0.0 {
            return Err(invalid("viscous plate thickness must be positive and finite"));
        }
        let factor = thickness_m*thickness_m*(thickness_m/12.0);
        if !factor.is_finite() || factor <= 0.0 {
            return Err(invalid("viscous section thickness integral is unrepresentable"));
        }
        let original = self.rotated_matrix_pa_s(angle)?;
        let result = original.map(|v| v*factor);
        if result.iter().zip(original).any(|(r,v)| !r.is_finite() || (v != 0.0 && *r == 0.0)) {
            return Err(invalid("viscous section reduction lost a coefficient"));
        }
        admit_passive(&result)?;
        Ok(result)
    }
}
