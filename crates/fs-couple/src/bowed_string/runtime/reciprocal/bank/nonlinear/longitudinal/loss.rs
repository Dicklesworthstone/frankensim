//! Select linear axial loss without borrowing transverse or fluid damping.
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
}
impl LongitudinalMaterialLoss {
    pub(in crate::bowed_string::runtime::reciprocal::bank::nonlinear) fn bind(&self, source: &MaterialStringSource)
        -> Result<Option<ConcentricAxialViscosity>, ImpactError>
    {
        let Self::KelvinVoigt { viscosity_property } = self else { return Ok(None); };
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
        };
        section.with_axial_viscosity(viscosity_property)
            .map(Some).map_err(|e| ImpactError::Owner(e.to_string()))
    }
}
