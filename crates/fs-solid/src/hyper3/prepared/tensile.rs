//! Sourced tensile limits on the original finite-strain P1 response.
//!
//! The same source-bound hyperelastic owner supplies mechanics and its explicitly
//! selected Cauchy-strength limit. Native recovery pushes P forward with the
//! actual total F and J; nominal engineering strength is never converted by an
//! inferred area change. These observations do not advance damage or certify
//! equilibrium, fracture initiation, or an unsampled continuum stress maximum.

use super::super::{HyperTetError, HyperTetProblem, checkpoint, response};
use super::PreparedHyperTetProblem;
use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_material::state_point::MaterialStatePointError;
use fs_material::strength::{ResolvedTensileHyperelastic, TensileAssessment};

/// One original element response and its selected true-stress tensile limit.
#[derive(Debug, Clone, PartialEq)]
pub struct HyperTetTensileSample {
    response: response::HyperTetElementResponse,
    material_binding_identity: ContentHash,
    limit_source_identity: ContentHash,
    limit_pa: f64,
    assessment: TensileAssessment,
}

impl HyperTetTensileSample {
    /// Unmodified P, F, J and Cauchy stress from the native response owner.
    pub const fn response(&self) -> &response::HyperTetElementResponse {
        &self.response
    }

    /// Complete material-card, model selection and Cauchy-limit binding.
    pub const fn material_binding_identity(&self) -> ContentHash {
        self.material_binding_identity
    }

    /// Original source bundle for the selected tensile property and conditions.
    pub const fn limit_source_identity(&self) -> ContentHash {
        self.limit_source_identity
    }

    /// Explicitly selected current-area tensile limit, Pa.
    pub const fn limit_pa(&self) -> f64 {
        self.limit_pa
    }

    pub const fn assessment(&self) -> &TensileAssessment {
        &self.assessment
    }
}

/// Immutable observations in original element order.
#[derive(Debug, Clone, PartialEq)]
pub struct HyperTetTensileAssessment {
    samples: Vec<HyperTetTensileSample>,
    governing_sample: Option<usize>,
}

impl HyperTetTensileAssessment {
    pub fn samples(&self) -> &[HyperTetTensileSample] {
        &self.samples
    }

    /// First sample attaining the largest utilization; ties are deterministic.
    pub fn governing_sample(&self) -> Option<&HyperTetTensileSample> {
        self.governing_sample.map(|index| &self.samples[index])
    }

    pub fn reaches_limit(&self) -> bool {
        self.governing_sample()
            .is_some_and(|sample| sample.assessment.reaches_limit())
    }
}

/// Native mechanics, source, numerical or resource refusal; no partial field.
#[derive(Debug, Clone, PartialEq)]
pub enum HyperTetTensileError {
    Native(HyperTetError),
    Material(MaterialStatePointError),
    SampleBudget {
        required: usize,
        maximum: usize,
    },
    Allocation,
    /// An objective constitutive response must be symmetric within roundoff.
    NonsymmetricCauchy {
        element: usize,
        relative_skew: f64,
    },
}

impl core::fmt::Display for HyperTetTensileError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite-strain tensile assessment refused: {self:?}")
    }
}
impl std::error::Error for HyperTetTensileError {}
impl From<HyperTetError> for HyperTetTensileError {
    fn from(error: HyperTetError) -> Self {
        Self::Native(error)
    }
}
impl From<MaterialStatePointError> for HyperTetTensileError {
    fn from(error: MaterialStatePointError) -> Self {
        Self::Material(error)
    }
}

impl HyperTetProblem<'_, ResolvedTensileHyperelastic> {
    /// Compare actual current-area stresses to the original material limits.
    ///
    /// The displacement can come from the native nonlinear solve or another
    /// admitted observation. This method does not itself establish equilibrium.
    /// Prepare once and use the prepared method for repeated trajectory samples.
    pub fn assess_tensile_strength(
        &self,
        displacement_m: &[f64],
        maximum_samples: usize,
        cx: &Cx<'_>,
    ) -> Result<HyperTetTensileAssessment, HyperTetTensileError> {
        checkpoint(cx)?;
        admit_sample_count(self.tetrahedra.len(), maximum_samples)?;
        self.prepare(cx)?
            .assess_tensile_strength(displacement_m, maximum_samples, cx)
    }
}

impl PreparedHyperTetProblem<'_, ResolvedTensileHyperelastic> {
    /// Assess each original source-bound element with its cached geometry.
    ///
    /// Only stress and potential are recovered: no tangent, assembled operator,
    /// second equilibrium solve or intermediate full response field is built.
    /// Work and output allocation are bounded by maximum_samples, with a
    /// cancellation checkpoint at every element and before result publication.
    pub fn assess_tensile_strength(
        &self,
        displacement_m: &[f64],
        maximum_samples: usize,
        cx: &Cx<'_>,
    ) -> Result<HyperTetTensileAssessment, HyperTetTensileError> {
        checkpoint(cx)?;
        admit_sample_count(self.elements.len(), maximum_samples)?;
        self.problem.check_len(displacement_m)?;
        let mut samples: Vec<HyperTetTensileSample> = Vec::new();
        samples
            .try_reserve_exact(self.elements.len())
            .map_err(|_| HyperTetTensileError::Allocation)?;
        let mut governing_sample = None;
        let mut largest = -1.0;
        for (element, geometry) in self.elements.iter().enumerate() {
            checkpoint(cx)?;
            let deformation = HyperTetProblem::<ResolvedTensileHyperelastic>::deformation_gradient(
                geometry,
                displacement_m,
            );
            let response = response::recover(&self.problem, element, geometry, deformation)?;
            let material = self.problem.material.at(element);
            let stress = symmetric_cauchy(&response.cauchy_pa, element)?;
            let assessment = material.limit().assess(&stress)?;
            if assessment.utilization() > largest {
                largest = assessment.utilization();
                governing_sample = Some(element);
            }
            samples.push(HyperTetTensileSample {
                response,
                material_binding_identity: material.identity(),
                limit_source_identity: material.limit().resolved().identity(),
                limit_pa: material.limit().limit_pa(),
                assessment,
            });
        }
        checkpoint(cx)?;
        Ok(HyperTetTensileAssessment {
            samples,
            governing_sample,
        })
    }
}

fn admit_sample_count(required: usize, maximum: usize) -> Result<(), HyperTetTensileError> {
    if required > maximum {
        Err(HyperTetTensileError::SampleBudget { required, maximum })
    } else {
        Ok(())
    }
}

/// Admit only roundoff-scale skew, retaining the original tensor in the sample.
/// The symmetric part is used solely for the isotropic principal-stress check.
/// This numerical tolerance is not an interval bound or a material allowance.
fn symmetric_cauchy(stress: &[f64; 9], element: usize) -> Result<[f64; 6], HyperTetTensileError> {
    let scale = stress.iter().fold(0.0_f64, |value, x| value.max(x.abs()));
    if scale == 0.0 {
        return Ok([0.0; 6]);
    }
    let relative_skew = [(1, 3), (5, 7), (6, 2)]
        .into_iter()
        .fold(0.0_f64, |value, (i, j)| {
            value.max((stress[i] / scale - stress[j] / scale).abs())
        });
    if relative_skew > 64.0 * f64::EPSILON {
        return Err(HyperTetTensileError::NonsymmetricCauchy {
            element,
            relative_skew,
        });
    }
    // After skew admission this cannot overflow. Equal subnormal entries stay
    // intact, unlike averaging separately halved components.
    let mean = |a: f64, b: f64| a + 0.5 * (b - a);
    Ok([
        stress[0],
        stress[4],
        stress[8],
        mean(stress[1], stress[3]),
        mean(stress[5], stress[7]),
        mean(stress[6], stress[2]),
    ])
}
