//! Direction-dependent tensile observation on the native solid stress field.
//!
//! The material axes are the same proper rotation used to assemble stiffness.
//! Each actual stress sample is pulled back to those axes before applying the
//! shared maximum material-normal-stress criterion. Stress principal axes are
//! not material axes. Neither spatial averaging nor an isotropic eigenvalue
//! criterion can substitute for this directional comparison.
//!
//! These small-strain observations do not implement a mixed-mode failure
//! surface, shear/compressive strength, fracture evolution or a continuum peak
//! bound. Named strengths must be applicable to the chosen material condition.

use super::{TensileAssessmentError, admit_same_state, checkpoint};
use crate::linear3::{QuadraticThermalResponse, TetElasticMaterial, TetThermalStressReport};
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_material::state_point::OrthotropicElasticStatePoint;
use fs_material::strength::{OrthotropicTensileAssessment, OrthotropicTensileLimit};

/// Source-matched directional strength and the exact oriented elastic state.
#[derive(Debug, Clone, PartialEq)]
pub struct OrthotropicTensileMaterialBinding {
    material: TetElasticMaterial,
    limit: OrthotropicTensileLimit,
    world_to_principal: [[f64; 3]; 3],
}

impl OrthotropicTensileMaterialBinding {
    /// Bind the named limits to their source-derived orthotropic stiffness.
    ///
    /// `principal_to_world` has the material-axis unit vectors as columns in
    /// the native mesh frame, exactly as in the solid constitutive constructor.
    /// Its validation and transpose are performed once at admission. The
    /// complete oriented stiffness identity is retained by native mechanics.
    pub fn try_from_resolved_orthotropic_state(
        state: &OrthotropicElasticStatePoint,
        limit: OrthotropicTensileLimit,
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

    /// Three explicitly selected strengths and their complete source receipts.
    #[must_use]
    pub const fn limit(&self) -> &OrthotropicTensileLimit {
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
#[derive(Debug, Clone, Copy)]
pub enum OrthotropicTensileMaterialField<'a> {
    Uniform(&'a OrthotropicTensileMaterialBinding),
    PerElement(&'a [OrthotropicTensileMaterialBinding]),
    /// Flattened element-major, native-quadrature-point order. Each local
    /// temperature, process and material orientation must match the response.
    PerPoint(&'a [OrthotropicTensileMaterialBinding]),
}

impl<'a> OrthotropicTensileMaterialField<'a> {
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
                "orthotropic tensile bindings must exactly cover the original response field",
            ))
        }
    }

    fn at(self, element: usize, point: usize) -> &'a OrthotropicTensileMaterialBinding {
        match self {
            Self::Uniform(value) => value,
            Self::PerElement(values) => &values[element],
            Self::PerPoint(values) => &values[point],
        }
    }
}

/// Directional tensile assessment of one unaveraged native stress sample.
#[derive(Debug, Clone, PartialEq)]
pub struct OrthotropicTensileStressSample {
    element: usize,
    quadrature_point: Option<usize>,
    material_state_identity: ContentHash,
    limit_source_identity: ContentHash,
    assessment: OrthotropicTensileAssessment,
}

impl OrthotropicTensileStressSample {
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
    pub const fn assessment(&self) -> &OrthotropicTensileAssessment {
        &self.assessment
    }
}

/// Atomic comparison of one original solid response with directional limits.
#[derive(Debug, Clone, PartialEq)]
pub struct OrthotropicSolidTensileAssessment {
    response_identity: ContentHash,
    samples: Vec<OrthotropicTensileStressSample>,
    governing_sample: Option<usize>,
    identity: ContentHash,
}

impl OrthotropicSolidTensileAssessment {
    #[must_use]
    pub fn samples(&self) -> &[OrthotropicTensileStressSample] {
        &self.samples
    }

    /// First sample with largest utilization; each sample also names its axis.
    #[must_use]
    pub fn governing_sample(&self) -> Option<&OrthotropicTensileStressSample> {
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
fn assess_samples(
    response_identity: ContentHash,
    elements: usize,
    points: usize,
    stresses: impl Iterator<Item = (usize, Option<usize>, [f64; 6], ContentHash)>,
    materials: OrthotropicTensileMaterialField<'_>,
    maximum_samples: usize,
    cx: &Cx<'_>,
) -> Result<OrthotropicSolidTensileAssessment, TensileAssessmentError> {
    checkpoint(cx)?;
    if points > maximum_samples {
        return Err(TensileAssessmentError::Budget);
    }
    materials.validate(elements, points)?;
    let mut samples: Vec<OrthotropicTensileStressSample> = Vec::new();
    samples
        .try_reserve_exact(points)
        .map_err(|_| TensileAssessmentError::Budget)?;
    let mut governing_sample = None;
    let mut largest = -1.0;
    let mut identity =
        DomainHasher::new("org.frankensim.fs-solid.orthotropic-tensile-assessment.v1");
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
            .assess_material_stress(&local_stress)
            .map_err(TensileAssessmentError::Material)?;
        if assessment.utilization() > largest {
            largest = assessment.utilization();
            governing_sample = Some(flat);
        }
        let limit_source_identity = binding.limit.resolved().identity();
        identity.update(material_state_identity.as_bytes());
        identity.update(limit_source_identity.as_bytes());
        for name in binding.limit.limit_properties() {
            identity.update(&(name.len() as u64).to_le_bytes());
            identity.update(name.as_bytes());
        }
        for value in binding.limit.limits_pa() {
            identity.update(&value.to_bits().to_le_bytes());
        }
        for value in local_stress {
            identity.update(&value.to_bits().to_le_bytes());
        }
        samples.push(OrthotropicTensileStressSample {
            element,
            quadrature_point: point,
            material_state_identity,
            limit_source_identity,
            assessment,
        });
    }
    checkpoint(cx)?;
    Ok(OrthotropicSolidTensileAssessment {
        response_identity,
        samples,
        governing_sample,
        identity: identity.finalize(),
    })
}
