//! Direction-dependent strength observation on the native solid stress field.
//!
//! The material axes are the same proper rotation used to assemble stiffness.
//! Each actual stress sample is pulled back to those axes before applying the
//! explicitly selected shared strength criterion. Stress principal axes are
//! not material axes. Neither spatial averaging nor an isotropic eigenvalue
//! criterion can substitute for this directional comparison.
//!
//! Maximum material-normal tension and the full 3-D Tsai--Wu interaction use
//! one sampling path. These small-strain observations do not evolve damage,
//! change mechanics or bound unsampled continuum peaks. Named strengths must
//! be applicable to the chosen material condition and stress convention.

mod criterion;
pub use criterion::{OrthotropicStrengthCriterion, OrthotropicStrengthObservation};

use super::{TensileAssessmentError, admit_same_state, checkpoint};
use crate::linear3::{QuadraticThermalResponse, TetElasticMaterial, TetThermalStressReport};
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_material::state_point::OrthotropicElasticStatePoint;
use fs_material::strength::{
    OrthotropicTensileAssessment, OrthotropicTensileLimit, OrthotropicTsaiWu,
    OrthotropicTsaiWuAssessment,
};

/// Original tensile-only material-axis criterion.
pub type OrthotropicTensileMaterialBinding =
    OrthotropicStrengthMaterialBinding<OrthotropicTensileLimit>;
/// Full normal/shear interaction with original source-defined strengths.
pub type OrthotropicTsaiWuMaterialBinding = OrthotropicStrengthMaterialBinding<OrthotropicTsaiWu>;
pub type OrthotropicTensileMaterialField<'a> =
    OrthotropicStrengthMaterialField<'a, OrthotropicTensileLimit>;
pub type OrthotropicTsaiWuMaterialField<'a> =
    OrthotropicStrengthMaterialField<'a, OrthotropicTsaiWu>;
pub type OrthotropicTensileStressSample =
    OrthotropicStrengthStressSample<OrthotropicTensileAssessment>;
pub type OrthotropicTsaiWuStressSample =
    OrthotropicStrengthStressSample<OrthotropicTsaiWuAssessment>;
pub type OrthotropicSolidTensileAssessment =
    OrthotropicSolidStrengthAssessment<OrthotropicTensileAssessment>;
pub type OrthotropicSolidTsaiWuAssessment =
    OrthotropicSolidStrengthAssessment<OrthotropicTsaiWuAssessment>;

/// Source-matched directional strength and the exact oriented elastic state.
#[derive(Debug, Clone, PartialEq)]
pub struct OrthotropicStrengthMaterialBinding<C: OrthotropicStrengthCriterion> {
    material: TetElasticMaterial,
    limit: C,
    world_to_principal: [[f64; 3]; 3],
}

impl<C: OrthotropicStrengthCriterion> OrthotropicStrengthMaterialBinding<C> {
    /// Bind the named limits to their source-derived orthotropic stiffness.
    ///
    /// `principal_to_world` has the material-axis unit vectors as columns in
    /// the native mesh frame, exactly as in the solid constitutive constructor.
    /// Its validation and transpose are performed once at admission. The
    /// complete oriented stiffness identity is retained by native mechanics.
    pub fn try_from_resolved_orthotropic_state(
        state: &OrthotropicElasticStatePoint,
        limit: C,
        principal_to_world: [[f64; 3]; 3],
        orientation_identity: ContentHash,
    ) -> Result<Self, TensileAssessmentError> {
        admit_same_state(state.resolved(), limit.resolved())?;
        let material = TetElasticMaterial::try_from_resolved_orthotropic_state(
            state,
            principal_to_world,
            orientation_identity,
        )
        .map_err(TensileAssessmentError::Elastic)?;
        Ok(Self {
            material,
            limit,
            world_to_principal: transpose(principal_to_world),
        })
    }

    /// The original oriented source-derived state for the solid operator.
    #[must_use]
    pub const fn material(&self) -> &TetElasticMaterial {
        &self.material
    }

    /// Explicitly selected criterion and its complete original source receipts.
    #[must_use]
    pub const fn limit(&self) -> &C {
        &self.limit
    }

    /// Columns are the admitted material axes in the native mesh frame.
    #[must_use]
    pub fn principal_to_world(&self) -> [[f64; 3]; 3] {
        transpose(self.world_to_principal)
    }
}

fn transpose(matrix: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    std::array::from_fn(|i| std::array::from_fn(|j| matrix[j][i]))
}

/// Bindings in the original native response order, without homogenization.
#[derive(Debug)]
pub enum OrthotropicStrengthMaterialField<'a, C: OrthotropicStrengthCriterion> {
    Uniform(&'a OrthotropicStrengthMaterialBinding<C>),
    PerElement(&'a [OrthotropicStrengthMaterialBinding<C>]),
    /// Flattened element-major, native-quadrature-point order. Each local
    /// temperature, process and material orientation must match the response.
    PerPoint(&'a [OrthotropicStrengthMaterialBinding<C>]),
}

impl<C: OrthotropicStrengthCriterion> Copy for OrthotropicStrengthMaterialField<'_, C> {}
impl<C: OrthotropicStrengthCriterion> Clone for OrthotropicStrengthMaterialField<'_, C> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<'a, C: OrthotropicStrengthCriterion> OrthotropicStrengthMaterialField<'a, C> {
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
                "orthotropic strength bindings must exactly cover the original response field",
            ))
        }
    }

    fn at(self, element: usize, point: usize) -> &'a OrthotropicStrengthMaterialBinding<C> {
        match self {
            Self::Uniform(value) => value,
            Self::PerElement(values) => &values[element],
            Self::PerPoint(values) => &values[point],
        }
    }
}

/// Selected strength assessment of one unaveraged native stress sample.
#[derive(Debug, Clone, PartialEq)]
pub struct OrthotropicStrengthStressSample<A: OrthotropicStrengthObservation> {
    element: usize,
    quadrature_point: Option<usize>,
    material_state_identity: ContentHash,
    limit_source_identity: ContentHash,
    assessment: A,
}

impl<A: OrthotropicStrengthObservation> OrthotropicStrengthStressSample<A> {
    #[must_use]
    pub const fn element(&self) -> usize {
        self.element
    }

    /// None denotes constant P1 stress; Some is the original quadrature index.
    #[must_use]
    pub const fn quadrature_point(&self) -> Option<usize> {
        self.quadrature_point
    }

    /// Original oriented constitutive state, including its material-axis frame.
    #[must_use]
    pub const fn material_state_identity(&self) -> ContentHash {
        self.material_state_identity
    }

    #[must_use]
    pub const fn limit_source_identity(&self) -> ContentHash {
        self.limit_source_identity
    }

    #[must_use]
    pub const fn assessment(&self) -> &A {
        &self.assessment
    }
}

/// Atomic comparison of one original solid response with directional limits.
#[derive(Debug, Clone, PartialEq)]
pub struct OrthotropicSolidStrengthAssessment<A: OrthotropicStrengthObservation> {
    response_identity: ContentHash,
    samples: Vec<OrthotropicStrengthStressSample<A>>,
    governing_sample: Option<usize>,
    identity: ContentHash,
}

impl<A: OrthotropicStrengthObservation> OrthotropicSolidStrengthAssessment<A> {
    #[must_use]
    pub fn samples(&self) -> &[OrthotropicStrengthStressSample<A>] {
        &self.samples
    }

    /// First sample with largest utilization of its selected criterion.
    #[must_use]
    pub fn governing_sample(&self) -> Option<&OrthotropicStrengthStressSample<A>> {
        self.governing_sample.map(|i| &self.samples[i])
    }

    #[must_use]
    pub fn reaches_limit(&self) -> bool {
        self.governing_sample()
            .is_some_and(|s| s.assessment.reaches_limit())
    }

    #[must_use]
    pub const fn response_identity(&self) -> ContentHash {
        self.response_identity
    }

    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
}

impl TetThermalStressReport {
    /// Compare every constant P1 stress in its original material-axis frame.
    /// This includes the same stress reports returned by inertial trajectories.
    pub fn assess_orthotropic_tensile_strength(
        &self,
        materials: OrthotropicTensileMaterialField<'_>,
        maximum_samples: usize,
        cx: &Cx<'_>,
    ) -> Result<OrthotropicSolidTensileAssessment, TensileAssessmentError> {
        self.assess_orthotropic_strength(materials, maximum_samples, cx)
    }

    /// Assess the explicitly selected material-axis criterion on every native
    /// P1 stress, retaining its original material and orientation identity.
    /// The criterion type determines whether only normal tension or the full
    /// normal/shear interaction is assessed; no criterion is inferred.
    pub fn assess_orthotropic_strength<C: OrthotropicStrengthCriterion>(
        &self,
        materials: OrthotropicStrengthMaterialField<'_, C>,
        maximum_samples: usize,
        cx: &Cx<'_>,
    ) -> Result<OrthotropicSolidStrengthAssessment<C::Assessment>, TensileAssessmentError> {
        assess_samples(
            self.identity(),
            self.elements().len(),
            self.elements().len(),
            self.elements().iter().enumerate().map(|(element, stress)| {
                (
                    element,
                    None,
                    *stress.stress_mandel_pa(),
                    stress.material_state_identity(),
                )
            }),
            materials,
            maximum_samples,
            cx,
        )
    }
}

impl QuadraticThermalResponse {
    /// Compare every native stiffness-quadrature stress without averaging.
    ///
    /// Rotated material axes and local strength conditions follow each original
    /// point. Standard/B-bar stress semantics are preserved. No stress bound is
    /// inferred between retained points or on an unsampled surface.
    pub fn assess_orthotropic_tensile_strength(
        &self,
        materials: OrthotropicTensileMaterialField<'_>,
        maximum_samples: usize,
        cx: &Cx<'_>,
    ) -> Result<OrthotropicSolidTensileAssessment, TensileAssessmentError> {
        self.assess_orthotropic_strength(materials, maximum_samples, cx)
    }

    /// Apply the selected strength surface to every original quadrature
    /// stress. Each local temperature and material-axis orientation must match
    /// its stiffness state. No element-average stress replaces these samples.
    pub fn assess_orthotropic_strength<C: OrthotropicStrengthCriterion>(
        &self,
        materials: OrthotropicStrengthMaterialField<'_, C>,
        maximum_samples: usize,
        cx: &Cx<'_>,
    ) -> Result<OrthotropicSolidStrengthAssessment<C::Assessment>, TensileAssessmentError> {
        checkpoint(cx)?;
        let mut points = 0_usize;
        for element in self.elements() {
            checkpoint(cx)?;
            points = points
                .checked_add(element.stress_points().len())
                .ok_or(TensileAssessmentError::Budget)?;
        }
        assess_samples(
            self.identity(),
            self.elements().len(),
            points,
            self.elements()
                .iter()
                .enumerate()
                .flat_map(|(element, response)| {
                    response
                        .stress_points()
                        .iter()
                        .enumerate()
                        .map(move |(point, stress)| {
                            (
                                element,
                                Some(point),
                                *stress.stress_mandel_pa(),
                                stress.material_state_identity(),
                            )
                        })
                }),
            materials,
            maximum_samples,
            cx,
        )
    }
}

#[allow(clippy::too_many_arguments)]
fn assess_samples<C: OrthotropicStrengthCriterion>(
    response_identity: ContentHash,
    elements: usize,
    points: usize,
    stresses: impl Iterator<Item = (usize, Option<usize>, [f64; 6], ContentHash)>,
    materials: OrthotropicStrengthMaterialField<'_, C>,
    maximum_samples: usize,
    cx: &Cx<'_>,
) -> Result<OrthotropicSolidStrengthAssessment<C::Assessment>, TensileAssessmentError> {
    checkpoint(cx)?;
    if points > maximum_samples {
        return Err(TensileAssessmentError::Budget);
    }
    materials.validate(elements, points)?;
    let mut samples: Vec<OrthotropicStrengthStressSample<C::Assessment>> = Vec::new();
    samples
        .try_reserve_exact(points)
        .map_err(|_| TensileAssessmentError::Budget)?;
    let mut governing_sample = None;
    let mut largest = -1.0;
    let mut identity = DomainHasher::new(C::identity_domain());
    identity.update(response_identity.as_bytes());
    for (flat, (element, point, stress, material_state_identity)) in stresses.enumerate() {
        checkpoint(cx)?;
        let binding = materials.at(element, flat);
        if material_state_identity.0 != binding.material.material_state_identity {
            return Err(TensileAssessmentError::SourceMismatch { element, point });
        }
        let local_stress = fs_material::tensor::rotate(
            &fs_material::tensor::from_mandel(&stress),
            &binding.world_to_principal,
        );
        let assessment = binding
            .limit
            .assess(&local_stress)
            .map_err(TensileAssessmentError::Material)?;
        if assessment.utilization() > largest {
            largest = assessment.utilization();
            governing_sample = Some(flat);
        }
        let limit_source_identity = binding.limit.resolved().identity();
        identity.update(material_state_identity.as_bytes());
        identity.update(limit_source_identity.as_bytes());
        binding.limit.update_identity(&mut identity);
        for value in local_stress {
            identity.update(&value.to_bits().to_le_bytes());
        }
        samples.push(OrthotropicStrengthStressSample {
            element,
            quadrature_point: point,
            material_state_identity,
            limit_source_identity,
            assessment,
        });
    }
    checkpoint(cx)?;
    Ok(OrthotropicSolidStrengthAssessment {
        response_identity,
        samples,
        governing_sample,
        identity: identity.finalize(),
    })
}
