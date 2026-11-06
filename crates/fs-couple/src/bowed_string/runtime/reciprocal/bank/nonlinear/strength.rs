//! Original material strengths applied to the actual nonlinear string state.

use super::{
    BodyPotential, ImpactError, MAX_IMPACT_MODES, MaterialStringSource, NonlinearMaterialString,
    invalid,
};
use crate::render::plate::impact::string::StringTensileAssessment;
use fs_material::section::strength::AxialTensileCapacity;
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};

impl MaterialStringSource {
    fn tensile_capacity(&self, property: &str) -> Result<AxialTensileCapacity, ImpactError> {
        let capacity = match self {
            Self::WithPolarizations(source) => return source.source().tensile_capacity(property),
            Self::WithTorsion(source) => return source.bending().tensile_capacity(property),
            Self::Homogeneous(specimen) => {
                if specimen.string().relaxation_bending.is_some() {
                    return Err(invalid(
                        "hereditary string tensile admission needs actual layer memory stresses",
                    ));
                }
                ConcentricElasticSection::try_new(
                    0.0,
                    vec![ConcentricMaterialLayer {
                        outer_radius_m: specimen.radius_m(),
                        material: specimen.material().clone(),
                    }],
                )
                .and_then(|section| section.axial_tensile_capacity(property))
            }
            Self::Concentric { section, .. } => section.elastic().axial_tensile_capacity(property),
            Self::ConcentricProny { .. } => {
                return Err(invalid(
                    "hereditary string tensile admission needs actual layer memory stresses",
                ));
            }
            Self::ConcentricModalLoss { section, .. } => {
                section.elastic().axial_tensile_capacity(property)
            }
            Self::Measured(_) | Self::Wound(_) => {
                return Err(invalid(
                    "measured or helical strings need an explicit constituent stress model for tensile admission",
                ));
            }
        };
        capacity.map_err(|error| ImpactError::Owner(error.to_string()))
    }
}

impl NonlinearMaterialString {
    /// Select a named sourced tensile limit on every load-bearing layer.
    ///
    /// The actual Kirchhoff--Carrier tension and a continuous-span bending
    /// curvature envelope give `sigma_i <= E_i (T/EA + r_i curvature_bound)`.
    /// Both transverse polarizations contribute to one tension and one vector
    /// curvature. Initial and future candidate states exceeding this envelope
    /// refuse through the original native acceptance gate, preserving accepted
    /// state, time and histories. Forces and energy are never clipped or changed.
    ///
    /// This explicitly selects a small-strain elastic normal-stress criterion.
    /// Source values are not converted to allowables. Viscous bending stress,
    /// torsional/shear failure, rate effects, necking and fracture are not modeled.
    /// A conservative modal envelope can refuse a cancelling shape. Every source
    /// condition stays fixed; no temperature evolution or material refresh occurs.
    ///
    /// # Errors
    /// Missing/nonstatic/invalid limits, measured, wound or hereditary sections,
    /// duplicate selection, inconsistent geometry and overlimit initial state.
    pub fn with_tensile_limit(mut self, property: &str) -> Result<Self, ImpactError> {
        let capacity = self.source.tensile_capacity(property)?;
        let BodyPotential::String(potential) = &self.bodies[0].potential else {
            return Err(invalid(
                "material tensile admission requires its original nonlinear string",
            ));
        };
        let potential = potential
            .clone()
            .with_tensile_capacity(capacity, self.card.length_m)?;
        let mut q = [0.0; MAX_IMPACT_MODES];
        for (index, state) in self.bodies[0].initial.iter().enumerate() {
            q[index] = state.displacement_m_sqrt_kg;
        }
        potential.observe(&q[..self.transverse_mode_count()])?;
        self.bodies[0].potential = BodyPotential::String(potential);
        Ok(self)
    }

    /// Complete original section and per-layer strength receipts, when selected.
    #[must_use]
    pub fn tensile_capacity(&self) -> Option<&AxialTensileCapacity> {
        match &self.bodies[0].potential {
            BodyPotential::String(potential) => potential.tensile_capacity(),
            _ => None,
        }
    }

    /// Diagnose elastic tensile demand in the complete transverse coordinate
    /// prefix (one or two planes, excluding torsion). Does not enforce the limit.
    pub fn tensile_assessment(
        &self,
        transverse_coordinates: &[f64],
    ) -> Result<Option<StringTensileAssessment>, ImpactError> {
        match &self.bodies[0].potential {
            BodyPotential::String(potential) => {
                potential.tensile_assessment(transverse_coordinates)
            }
            _ => Err(invalid(
                "material tensile assessment requires its original nonlinear string",
            )),
        }
    }
}

impl super::runtime::MaterialStringRun {
    /// Elastic tensile envelope of the accepted native state, including both
    /// transverse planes. Connected receiver and torsion coordinates are excluded.
    /// Uses bounded stack scratch and the original immutable strength binding.
    pub fn tensile_assessment(&self) -> Result<Option<StringTensileAssessment>, ImpactError> {
        let count = self.source().transverse_mode_count();
        let mut q = [0.0; MAX_IMPACT_MODES];
        for (index, value) in q[..count].iter_mut().enumerate() {
            *value = self.system().state()[2 * index];
        }
        self.source().tensile_assessment(&q[..count])
    }
}
