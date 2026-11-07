//! Material-bound tensile assessments of actual small-strain solid stresses.
//!
//! The shared `fs-material` maximum-principal-tension criterion is evaluated
//! separately at every retained point. Neither compression nor averaging can
//! cancel a tensile demand. These are observations of the accepted mechanics,
//! not constitutive fracture, element deletion, or continuum safety bounds.

use crate::linear3::{QuadraticThermalResponse, TetElasticMaterial, TetThermalStressReport};
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_material::state_point::{
    IsotropicElasticStatePoint, IsotropicSolidStatePoint, MaterialStatePointError,
    ResolvedMaterialStatePoint,
};
use fs_material::strength::{IsotropicTensileLimit, TensileAssessment};

/// Direction-dependent tensile limits in the original material-axis frame.
pub mod orthotropic;

/// A tensile criterion and the source-derived elastic state it may assess.
///
/// Elasticity and strength can have different property selections, but must
/// describe the same immutable card, material condition and physical query.
/// Construct the mechanics with [`Self::material`] or the same resolved state.
/// A limit at another temperature, moisture content or process is refused.
#[derive(Debug, Clone, PartialEq)]
pub struct TensileMaterialBinding {
    material: TetElasticMaterial,
    limit: IsotropicTensileLimit,
}

impl TensileMaterialBinding {
    /// Bind a separately selected strength claim to the native elastic state.
    pub fn from_isotropic_elastic(
        state: &IsotropicElasticStatePoint,
        limit: IsotropicTensileLimit,
    ) -> Result<Self, TensileAssessmentError> {
        admit_same_state(state.resolved(), limit.resolved())?;
        Ok(Self {
            material: TetElasticMaterial::from_resolved_elastic_state(state),
            limit,
        })
    }

    /// Bind strength to an isotropic state that also retains a yield source.
    /// This does not replace its yield surface with the tensile criterion.
    pub fn from_isotropic_solid(
        state: &IsotropicSolidStatePoint,
        limit: IsotropicTensileLimit,
    ) -> Result<Self, TensileAssessmentError> {
        admit_same_state(state.resolved(), limit.resolved())?;
        Ok(Self {
            material: TetElasticMaterial::from_resolved_state(state),
            limit,
        })
    }

    /// The source-derived state to consume in the original solid operator.
    #[must_use]
    pub const fn material(&self) -> &TetElasticMaterial {
        &self.material
    }

    /// Complete strength-source receipt and explicit scalar criterion.
    #[must_use]
    pub const fn limit(&self) -> &IsotropicTensileLimit {
        &self.limit
    }
}

pub(crate) fn admit_same_state(
    mechanics: &ResolvedMaterialStatePoint,
    strength: &ResolvedMaterialStatePoint,
) -> Result<(), TensileAssessmentError> {
    if mechanics.card_identity() != strength.card_identity()
        || mechanics.material() != strength.material()
        || mechanics.query_point() != strength.query_point()
        || mechanics.axis_quantities() != strength.axis_quantities()
    {
        return Err(TensileAssessmentError::Invalid(
            "elasticity and tensile limit must name the same material card and physical state",
        ));
    }
    Ok(())
}

/// Source bindings in the original response's element/point order.
#[derive(Debug, Clone, Copy)]
pub enum TensileMaterialField<'a> {
    Uniform(&'a TensileMaterialBinding),
    PerElement(&'a [TensileMaterialBinding]),
    /// Flattened element-major, native-quadrature order. This permits local
    /// temperature-dependent limits alongside local thermoelastic stiffness.
    PerPoint(&'a [TensileMaterialBinding]),
}

impl<'a> TensileMaterialField<'a> {
    fn validate(self, elements: usize, points: usize) -> Result<(), TensileAssessmentError> {
        let valid = match self {
            Self::Uniform(_) => true,
            Self::PerElement(values) => values.len() == elements,
            Self::PerPoint(values) => values.len() == points,
        };
        if valid {
            Ok(())
        } else {
            Err(TensileAssessmentError::Invalid(
                "tensile bindings must exactly cover the original response field",
            ))
        }
    }

    fn at(self, element: usize, point: usize) -> &'a TensileMaterialBinding {
        match self {
            Self::Uniform(value) => value,
            Self::PerElement(values) => &values[element],
            Self::PerPoint(values) => &values[point],
        }
    }
}

/// Strength result at one actual native stress sample.
#[derive(Debug, Clone, PartialEq)]
pub struct TensileStressSample {
    element: usize,
    quadrature_point: Option<usize>,
    material_state_identity: ContentHash,
    limit_source_identity: ContentHash,
    assessment: TensileAssessment,
}

impl TensileStressSample {
    /// Index in the native element connectivity.
    #[must_use]
    pub const fn element(&self) -> usize {
        self.element
    }

    /// Native stiffness-point index; None means a constant P1 element stress.
    #[must_use]
    pub const fn quadrature_point(&self) -> Option<usize> {
        self.quadrature_point
    }

    #[must_use]
    pub const fn material_state_identity(&self) -> ContentHash {
        self.material_state_identity
    }

    /// Identity of the complete source bundle behind the selected limit.
    #[must_use]
    pub const fn limit_source_identity(&self) -> ContentHash {
        self.limit_source_identity
    }

    #[must_use]
    pub const fn assessment(&self) -> &TensileAssessment {
        &self.assessment
    }
}

/// Transactional tensile comparison of one accepted physical response.
#[derive(Debug, Clone, PartialEq)]
pub struct SolidTensileAssessment {
    response_identity: ContentHash,
    samples: Vec<TensileStressSample>,
    governing_sample: Option<usize>,
    identity: ContentHash,
}

impl SolidTensileAssessment {
    #[must_use]
    pub const fn response_identity(&self) -> ContentHash {
        self.response_identity
    }

    #[must_use]
    pub fn samples(&self) -> &[TensileStressSample] {
        &self.samples
    }

    /// First sample attaining the largest utilization, for deterministic ties.
    #[must_use]
    pub fn governing_sample(&self) -> Option<&TensileStressSample> {
        self.governing_sample.map(|index| &self.samples[index])
    }

    #[must_use]
    pub fn reaches_limit(&self) -> bool {
        self.governing_sample()
            .is_some_and(|sample| sample.assessment.reaches_limit())
    }

    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
}

/// Native source, size or cancellation refusal; no partial result is returned.
#[derive(Debug, Clone, PartialEq)]
pub enum TensileAssessmentError {
    Invalid(&'static str),
    Budget,
    SourceMismatch { element: usize, point: Option<usize> },
    Material(MaterialStatePointError),
    Elastic(crate::linear3::TetElasticError),
    Cancelled,
}
impl core::fmt::Display for TensileAssessmentError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "solid tensile assessment refused: {self:?}")
    }
}
impl std::error::Error for TensileAssessmentError {}

impl TetThermalStressReport {
    /// Compare every constant P1 stress with its material's tensile limit.
    ///
    /// This also consumes the original response inside transient runs. Mandel
    /// shear components are converted exactly once to tensor-shear components.
    /// `maximum_samples` caps allocation and work before any assessment.
    pub fn assess_tensile_strength(
        &self,
        materials: TensileMaterialField<'_>,
        maximum_samples: usize,
        cx: &Cx<'_>,
    ) -> Result<SolidTensileAssessment, TensileAssessmentError> {
        assess_samples(
            self.identity(),
            self.elements().len(),
            self.elements().len(),
            self.elements().iter().enumerate().map(|(element, stress)| {
                (element, None, *stress.stress_mandel_pa(), stress.material_state_identity())
            }),
            materials,
            maximum_samples,
            cx,
        )
    }
}

impl QuadraticThermalResponse {
    /// Compare every original stiffness-quadrature stress, without averaging.
    ///
    /// Q1/Q2 bending and thermal gradients can exceed a limit even when their
    /// element mean vanishes. These retained samples do not bound stresses
    /// between quadrature points or at an unsampled surface. B-bar uses its
    /// admitted assumed strain, consistently with the native force and energy.
    pub fn assess_tensile_strength(
        &self,
        materials: TensileMaterialField<'_>,
        maximum_samples: usize,
        cx: &Cx<'_>,
    ) -> Result<SolidTensileAssessment, TensileAssessmentError> {
        checkpoint(cx)?;
        let mut points = 0_usize;
        for element in self.elements() {
            checkpoint(cx)?;
            points = points.checked_add(element.stress_points().len())
                .ok_or(TensileAssessmentError::Budget)?;
        }
        assess_samples(
            self.identity(),
            self.elements().len(),
            points,
            self.elements().iter().enumerate().flat_map(|(element, response)| {
                response.stress_points().iter().enumerate().map(move |(point, stress)| {
                    (element, Some(point), *stress.stress_mandel_pa(), stress.material_state_identity())
                })
            }),
            materials,
            maximum_samples,
            cx,
        )
    }
}

#[allow(clippy::too_many_arguments)]
fn assess_samples(
    response_identity: ContentHash,
    elements: usize,
    points: usize,
    stresses: impl Iterator<Item = (usize, Option<usize>, [f64; 6], ContentHash)>,
    materials: TensileMaterialField<'_>,
    maximum_samples: usize,
    cx: &Cx<'_>,
) -> Result<SolidTensileAssessment, TensileAssessmentError> {
    checkpoint(cx)?;
    if points > maximum_samples {
        return Err(TensileAssessmentError::Budget);
    }
    materials.validate(elements, points)?;
    let mut samples: Vec<TensileStressSample> = Vec::with_capacity(points);
    let mut governing_sample = None;
    let mut largest = -1.0;
    let mut identity = DomainHasher::new("org.frankensim.fs-solid.tensile-assessment.v1");
    identity.update(response_identity.as_bytes());
    for (flat, (element, point, stress, material_state_identity)) in stresses.enumerate() {
        checkpoint(cx)?;
        let binding = materials.at(element, flat);
        if material_state_identity.0 != binding.material.material_state_identity {
            return Err(TensileAssessmentError::SourceMismatch { element, point });
        }
        let assessment = binding.limit.assess(&fs_material::tensor::from_mandel(&stress))
            .map_err(TensileAssessmentError::Material)?;
        if assessment.utilization() > largest {
            largest = assessment.utilization();
            governing_sample = Some(flat);
        }
        let limit_source_identity = binding.limit.resolved().identity();
        identity.update(material_state_identity.as_bytes());
        identity.update(limit_source_identity.as_bytes());
        identity.update(&binding.limit.limit_pa().to_bits().to_le_bytes());
        for value in assessment.principal_stresses_pa() {
            identity.update(&value.to_bits().to_le_bytes());
        }
        samples.push(TensileStressSample {
            element, quadrature_point: point, material_state_identity,
            limit_source_identity, assessment,
        });
    }
    checkpoint(cx)?;
    Ok(SolidTensileAssessment {
        response_identity, samples, governing_sample, identity: identity.finalize(),
    })
}

pub(crate) fn checkpoint(cx: &Cx<'_>) -> Result<(), TensileAssessmentError> {
    cx.checkpoint().map_err(|_| TensileAssessmentError::Cancelled)
}
