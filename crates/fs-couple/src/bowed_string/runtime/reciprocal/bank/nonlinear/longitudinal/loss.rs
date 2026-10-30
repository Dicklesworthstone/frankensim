//! Select source-bound axial constitutive loss without borrowing fluid damping.
use super::{ImpactError, MaterialStringSource, invalid};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::section::axial::ConcentricAxialViscosity;

/// Explicit constitutive selection, independent of the nonlinear strain energy.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum LongitudinalMaterialLoss {
    /// Original lossless longitudinal coordinates. No viscosity is inferred.
    Elastic,
    /// Linear axial strain-rate loss from each original solid state [Pa s].
    /// This does not include the finite-slope y_x*y_xt contribution to strain rate.
    KelvinVoigt {
        /// Explicit uniaxial property; every layer must supply it in its bundle.
        viscosity_property: String,
    },
    /// The same source viscosity applied to the complete retained strain rate
    /// `u_xt + y_x*y_xt + z_x*z_xt`. The original joint nonlinear/contact solve
    /// owns the resulting force and loss; linear modal axial damping is absent.
    /// This remains the declared moderate-slope strain reduction, not a fully
    /// geometrically exact material law or a thermal evolution model.
    FiniteStrainKelvinVoigt {
        /// Explicit uniaxial viscosity [Pa s] from every original solid layer.
        viscosity_property: String,
    },
}
impl LongitudinalMaterialLoss {
    pub(super) fn finite_strain(&self) -> bool {
        matches!(self, Self::FiniteStrainKelvinVoigt { .. })
    }
    pub(in crate::bowed_string::runtime::reciprocal::bank::nonlinear) fn bind(&self, source: &MaterialStringSource)
        -> Result<Option<ConcentricAxialViscosity>, ImpactError>
    {
        let viscosity_property = match self {
            Self::Elastic => return Ok(None),
            Self::KelvinVoigt { viscosity_property }
            | Self::FiniteStrainKelvinVoigt { viscosity_property } => viscosity_property,
        };
        let section = match source {
            MaterialStringSource::WithPolarizations(p) => return self.bind(p.source()),
            MaterialStringSource::WithTorsion(t) => return self.bind(t.bending()),
            MaterialStringSource::Homogeneous(s) => ConcentricElasticSection::try_new(0.0,
                vec![ConcentricMaterialLayer {
                    outer_radius_m: s.radius_m(), material: s.material().clone(),
                }]).map_err(|e| ImpactError::Owner(e.to_string()))?,
            MaterialStringSource::Concentric { section, .. } => section.elastic().clone(),
            MaterialStringSource::ConcentricProny { section, .. } => section.elastic().clone(),
            MaterialStringSource::Wound(_) => return Err(invalid("winding axial loss needs its own resolved constitutive model")),
            MaterialStringSource::Measured(_) => return Err(invalid("whole-string measurements contain no axial viscosity law")),
        };
        section.with_axial_viscosity(viscosity_property)
            .map(Some).map_err(|e| ImpactError::Owner(e.to_string()))
    }
}
