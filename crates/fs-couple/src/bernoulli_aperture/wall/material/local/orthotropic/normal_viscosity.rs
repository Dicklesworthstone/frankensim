//! Exact zero-axial-stress elimination of a two-dimensional Kelvin--Voigt law.
//!
//! With normal strain e=[e_h,e_a], sigma=C*e+D*e', eliminating sigma_a=0
//! means eliminating C+sD AS A WHOLE. Condensing C and D independently drops
//! a physical relaxation mode whenever their normal couplings differ:
//!
//! M(s) = C_hh+s D_hh - (C_ha+s D_ha)^2/(C_aa+s D_aa)
//!      = E_h + eta*s + delta_E*s*tau/(1+s*tau),
//! eta=D_hh-D_ha^2/D_aa, tau=D_aa/C_aa,
//! delta_E=C_aa*(C_ha/C_aa-D_ha/D_aa)^2.
//!
//! The existing Maxwell state and Kelvin resistance realize this identity;
//! no frequency fit, effective material card or second time integrator exists.
use super::*;
use fs_material::visco::GeneralizedMaxwell;
use fs_math::det;

/// Schur complement of a positive-semidefinite normal viscosity matrix.
/// Exact rank deficiency is admissible; no negative value is clipped to zero.
fn residual_viscosity(d: [f64; 3]) -> Result<f64, AcousticRealizeError> {
    let [hh, ha, aa] = d;
    if d.iter().any(|x| !x.is_finite()) || hh < 0.0 || aa < 0.0 {
        return Err(invalid("normal wall viscosity needs finite nonnegative diagonal entries"));
    }
    if hh == 0.0 || aa == 0.0 {
        if ha != 0.0 {
            return Err(invalid("a zero diagonal viscosity requires zero cross viscosity"));
        }
        return Ok(hh);
    }
    let correlation = (ha / det::sqrt(hh)) / det::sqrt(aa);
    if !correlation.is_finite() || correlation.abs() > 1.0
        || (ha != 0.0 && correlation == 0.0)
    {
        return Err(invalid("normal wall viscosity must be positive semidefinite and representable"));
    }
    let determinant = (1.0 - correlation.abs()) * (1.0 + correlation.abs());
    let residual = hh * determinant;
    if !residual.is_finite() || (determinant > 0.0 && residual == 0.0) {
        return Err(invalid("normal wall viscous Schur complement is not representable"));
    }
    Ok(residual)
}

impl LocalMaterialWall {
    /// Requirements for an independent symmetric normal-plane viscosity D.
    /// Names are ordered [D_hh,D_ha,D_aa], all in Pa s in the declared local
    /// hoop/axial frame. These are stiffness-viscosity tensor entries, NOT
    /// directional Young viscosities or a three-dimensional tensor awaiting
    /// radial plane-stress elimination. The cross term may be negative.
    pub fn orthotropic_normal_viscosity_requirements(
        axes: [usize; 2],
        viscosity_properties: [&str; 3],
    ) -> Result<Vec<ScalarPropertyRequirement>, AcousticRealizeError> {
        let mut requirements = Self::orthotropic_requirements(axes, None)?;
        for name in viscosity_properties {
            if name.trim() != name {
                return Err(invalid("normal wall viscosity needs exact source property names"));
            }
            if let Some(prior) = requirements.iter().find(|r| r.name() == name) {
                if prior.dims() != DynViscosity::DIMS {
                    return Err(invalid("normal wall viscosity property conflicts with an elastic quantity"));
                }
            } else {
                requirements.push(ScalarPropertyRequirement::try_new(
                    name, DynViscosity::DIMS, ScalarAdmissibility::Finite,
                ).map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?);
            }
        }
        Ok(requirements)
    }

    /// Couple independently sourced normal Kelvin--Voigt viscosity to the
    /// original directional wall, including its exact axial relaxation state.
    ///
    /// C is the original reciprocal orthotropic normal-plane elasticity, D is
    /// the supplied symmetric positive-semidefinite two-dimensional viscosity.
    /// Fixed axial strain selects C_hh+s D_hh. Free axial stress eliminates
    /// C+sD together, retaining the resulting positive Maxwell arm AND the
    /// residual Kelvin resistance in the same gas/wall DG operator. With
    /// D_aa=0, passivity requires D_ha=0 and no internal axial state is needed.
    /// Proportional coupling has exactly zero additional Maxwell strength.
    ///
    /// All original material queries, axes, density and thickness are retained.
    /// The generated spectrum records the source state and named coefficients;
    /// it is an algebraic constitutive reduction, not another material card.
    /// Its memory retains the original tensor claims and axial-strain identity.
    /// The live wall wrapper preserves that physical strain as coefficients
    /// change, rather than falsely treating the generated arm as a source pair.
    ///
    /// The explicit temporal limits have the same domain as WallRelaxation::new;
    /// exact-DG propagation does not impose those midpoint-network limits.
    /// Existing thin-wall, fixed-state, principal-aligned local geometry applies.
    /// There is no axial wave inertia or full 3D viscoelastic condensation here.
    #[allow(clippy::too_many_arguments)]
    pub fn new_orthotropic_normal_kelvin_voigt(
        material: ResolvedMaterialStatePoint,
        thickness_m: f64,
        axial: AxialConstraint,
        axes: [usize; 2],
        viscosity_properties: [&str; 3],
        band_hz: [f64; 2],
        max_dt_over_tau: f64,
        max_angular_step: f64,
    ) -> Result<Self, AcousticRealizeError> {
        Self::orthotropic_normal_viscosity_requirements(axes, viscosity_properties)?;
        if !max_dt_over_tau.is_finite() || !(0.0 < max_dt_over_tau && max_dt_over_tau <= 2.0)
            || !max_angular_step.is_finite() || !(0.0 < max_angular_step && max_angular_step <= 1.0)
        {
            return Err(invalid("normal wall viscosity needs explicit bounded relaxation clock limits"));
        }
        // Reuse the original material/band/positive-compliance admission. The
        // fixed-strain factor is needed even when the ACTUAL constraint is free.
        let mut wall = Self::new_orthotropic(material, thickness_m,
            AxialConstraint::StrainFixed, axes, None, band_hz)?;
        let mut d = [0.0; 3];
        for (value, name) in d.iter_mut().zip(viscosity_properties) {
            *value = property(&wall.material, name, DynViscosity::DIMS, band_hz)?;
        }
        let residual = residual_viscosity(d)?;
        let constrained_hoop = wall.young * wall.axial_factor;
        let equilibrium = wall.young;
        wall.axial = axial;
        wall.viscosity = d[0];
        wall.loss = WallMaterialLoss::KelvinVoigt;
        // Coefficients below are already reduced in the physical hoop frame.
        // The radius/thickness/area map remains the ORIGINAL base_patch owner.
        let elastic_factor = wall.axial_factor;
        wall.axial_factor = 1.0;
        if axial == AxialConstraint::StrainFixed {
            wall.young = constrained_hoop;
            return Ok(wall);
        }
        if d[2] == 0.0 {
            return Ok(wall); // D_ha=0 has already been established, not defaulted.
        }
        wall.viscosity = residual;
        let (_, axial_name, poisson_name) = axes_names(axes)?;
        let axial_young = property(&wall.material, axial_name, Pressure::DIMS, band_hz)?;
        let nu = property(&wall.material, poisson_name, Dims::NONE, band_hz)?;
        let axial_modulus = axial_young * elastic_factor;
        let elastic_ratio = if axes[0] < axes[1] { nu } else { (nu / axial_young) * equilibrium };
        let viscous_ratio = d[1] / d[2];
        if !axial_modulus.is_finite() || axial_modulus <= 0.0
            || !elastic_ratio.is_finite() || (nu != 0.0 && elastic_ratio == 0.0)
            || !viscous_ratio.is_finite() || (d[1] != 0.0 && viscous_ratio == 0.0)
        {
            return Err(invalid("normal wall constitutive coupling is not representable"));
        }
        let mismatch = elastic_ratio - viscous_ratio;
        if mismatch == 0.0 {
            return Ok(wall); // The exact constitutive identity has no memory arm.
        }
        let root = det::sqrt(axial_modulus) * mismatch;
        let modulus = root * root;
        let time = d[2] / axial_modulus;
        if !modulus.is_finite() || modulus <= 0.0 || !time.is_finite() || time <= 0.0
            || !time.recip().is_finite()
        {
            return Err(invalid("normal wall condensation lost a positive relaxation coefficient"));
        }
        let law = GeneralizedMaxwell::new(equilibrium, vec![(modulus, time)])
            .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        let spectrum = super::super::super::WallRelaxation::new(law, format!(
            "exact normal-plane Kelvin-Voigt elimination; material-state:{}; axes:{axes:?}; viscosity-properties:{viscosity_properties:?}",
            wall.material.identity()), band_hz, max_dt_over_tau, max_angular_step)?;
        let spectrum = spectrum.with_normal_plane_history(&wall.material, axes,
            viscosity_properties, elastic_ratio, mismatch)?;
        wall.loss = WallMaterialLoss::Maxwell(spectrum);
        Ok(wall)
    }
}
