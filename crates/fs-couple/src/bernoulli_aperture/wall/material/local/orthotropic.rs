//! Principal-aligned plane-stress membrane, not an equivalent isotropic card.
//!
//! In the hoop/axial plane, S_ii=1/E_i and S_ij=-nu_ij/E_i for i<j.
//! Radial traction is neglected in the existing thin-wall limit. Free axial
//! stress gives E_hoop; fixed axial strain gives 1/(S_hh-S_ha^2/S_aa).
//! The symmetric normal compliance block must be positive definite. Shear and
//! the unused radial coefficients do not enter this axisymmetric local problem.
//! Off-axis grain, axial wave inertia, bending and thick walls remain excluded.

use super::*;
use fs_material::state_point::{
    ORTHOTROPIC_POISSON_RATIO_PROPERTIES, ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES,
    ScalarAdmissibility, ScalarPropertyRequirement,
};

const TIME_DIMS: Dims = Dims([0, 0, 1, 0, 0, 0]);

fn axes_names(axes: [usize; 2]) -> Result<(&'static str, &'static str, &'static str), AcousticRealizeError> {
    let [hoop, axial] = axes;
    if hoop >= 3 || axial >= 3 || hoop == axial {
        return Err(invalid("orthotropic wall needs distinct hoop and axial principal axes in 0..3"));
    }
    let pair = match (hoop.min(axial), hoop.max(axial)) {
        (0, 1) => 0,
        (0, 2) => 1,
        (1, 2) => 2,
        _ => unreachable!("distinct principal axes were checked"),
    };
    Ok((ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES[hoop],
        ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES[axial],
        ORTHOTROPIC_POISSON_RATIO_PROPERTIES[pair]))
}

impl LocalMaterialWall {
    /// Original scalar requirements for the selected principal membrane plane.
    /// Axes are zero-based material coordinates, ordered [hoop, axial].
    /// A supplied time property explicitly declares proportional Kelvin--Voigt
    /// response C_viscous=tau*C_elastic on this plane; it is not a fitted loss
    /// factor or permission to combine unrelated directional viscosities.
    pub fn orthotropic_requirements(
        axes: [usize; 2],
        proportional_viscosity_time_property: Option<&str>,
    ) -> Result<Vec<ScalarPropertyRequirement>, AcousticRealizeError> {
        let (hoop, axial, poisson) = axes_names(axes)?;
        let mut fields = vec![
            (DENSITY_PROPERTY, Density::DIMS, ScalarAdmissibility::StrictlyPositive),
            (hoop, Pressure::DIMS, ScalarAdmissibility::StrictlyPositive),
            (axial, Pressure::DIMS, ScalarAdmissibility::StrictlyPositive),
            (poisson, Dims::NONE, ScalarAdmissibility::Finite),
        ];
        if let Some(name) = proportional_viscosity_time_property {
            if name.trim() != name || fields.iter().any(|(key, _, _)| *key == name) {
                return Err(invalid("proportional wall viscosity needs a distinct exact source property"));
            }
            fields.push((name, TIME_DIMS, ScalarAdmissibility::NonNegative));
        }
        fields.into_iter().map(|(name, dims, domain)| {
            ScalarPropertyRequirement::try_new(name, dims, domain)
                .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))
        }).collect()
    }

    /// Source-backed orthotropic radial motion through the original wall owner.
    ///
    /// Material axes must follow the local cylindrical hoop/axial directions;
    /// this does not reinterpret a fixed Cartesian grain direction as a radial
    /// texture. The original density, two directional moduli and major Poisson
    /// ratio retain their receipts. No effective material or missing coefficient
    /// is invented. Positivity is checked for the represented normal compliance
    /// block, not asserted for an unspecified full three-dimensional material.
    ///
    /// With an explicit time property, C_viscous=tau*C_elastic shares the same
    /// reciprocal reduction for either axial restraint. General nonproportional
    /// anisotropic viscosity/relaxation needs additional constitutive states and
    /// is NOT approximated by this constructor. Gas losses remain independent.
    pub fn new_orthotropic(
        material: ResolvedMaterialStatePoint,
        thickness_m: f64,
        axial: AxialConstraint,
        axes: [usize; 2],
        proportional_viscosity_time_property: Option<&str>,
        band_hz: [f64; 2],
    ) -> Result<Self, AcousticRealizeError> {
        // The same source requirements and geometry admission as the native
        // input path; direct callers cannot bypass them with a partial bundle.
        Self::orthotropic_requirements(axes, proportional_viscosity_time_property)?;
        if !thickness_m.is_finite() || thickness_m <= 0.0
            || !band_hz.iter().all(|f| f.is_finite() && *f > 0.0)
            || band_hz[0] >= band_hz[1]
            || !(core::f64::consts::TAU * band_hz[1]).is_finite()
        {
            return Err(invalid("orthotropic wall requires positive thickness and an ordered positive frequency band"));
        }
        let (hoop_name, axial_name, poisson_name) = axes_names(axes)?;
        let density = property(&material, DENSITY_PROPERTY, Density::DIMS, band_hz)?;
        let young = property(&material, hoop_name, Pressure::DIMS, band_hz)?;
        let axial_young = property(&material, axial_name, Pressure::DIMS, band_hz)?;
        let nu = property(&material, poisson_name, Dims::NONE, band_hz)?;
        if density <= 0.0 || young <= 0.0 || axial_young <= 0.0 {
            return Err(invalid("orthotropic wall needs positive density and both principal moduli"));
        }
        // nu_ij is strain in j under stress in i, always i<j. Reciprocity
        // gives nu_ji=nu_ij*E_j/E_i, NOT nu_ji=nu_ij. Normalize before the
        // determinant check instead of forming tiny dimensional compliances.
        let (e_i, e_j) = if axes[0] < axes[1] { (young, axial_young) }
            else { (axial_young, young) };
        let coupling = (nu / e_i.sqrt()) * e_j.sqrt();
        if !coupling.is_finite() || coupling.abs() >= 1.0
            || (nu != 0.0 && coupling == 0.0)
        {
            return Err(invalid("orthotropic wall normal compliance must be positive definite and representable"));
        }
        let determinant = (1.0 - coupling.abs()) * (1.0 + coupling.abs());
        let axial_factor = match axial {
            AxialConstraint::StressFree => 1.0,
            AxialConstraint::StrainFixed => determinant.recip(),
        };
        if !axial_factor.is_finite() || !(young * axial_factor).is_finite() {
            return Err(invalid("orthotropic constrained hoop stiffness is not representable"));
        }
        let (viscosity, loss) = match proportional_viscosity_time_property {
            None => (0.0, WallMaterialLoss::Elastic),
            Some(name) => {
                let time = property(&material, name, TIME_DIMS, band_hz)?;
                let viscosity = young * time;
                if time < 0.0 || !viscosity.is_finite()
                    || !(viscosity * axial_factor).is_finite()
                    || (time > 0.0 && viscosity == 0.0)
                {
                    return Err(invalid("orthotropic proportional viscosity must be nonnegative and representable"));
                }
                (viscosity, WallMaterialLoss::KelvinVoigt)
            }
        };
        Ok(Self {
            material, thickness_m, axial, loss, band_hz, density, young,
            viscosity, axial_factor, orthotropic_axes: Some(axes),
            proportional_viscosity_time_property: proportional_viscosity_time_property.map(str::to_owned),
        })
    }

    /// Explicit zero-based [hoop, axial] source axes; None is the isotropic path.
    #[must_use]
    pub const fn orthotropic_axes(&self) -> Option<[usize; 2]> {
        self.orthotropic_axes
    }

    /// Original time-property name for C_viscous=tau*C_elastic, when requested.
    #[must_use]
    pub fn proportional_viscosity_time_property(&self) -> Option<&str> {
        self.proportional_viscosity_time_property.as_deref()
    }
}
