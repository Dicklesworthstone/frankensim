//! Sealed criteria sharing the native directional stress sampling path.

use fs_blake3::DomainHasher;
use fs_material::Voigt;
use fs_material::state_point::{MaterialStatePointError, ResolvedMaterialStatePoint};
use fs_material::strength::{
    OrthotropicTensileAssessment, OrthotropicTensileLimit, OrthotropicTsaiWu,
    OrthotropicTsaiWuAssessment,
};

mod sealed {
    pub trait Criterion {}
    pub trait Assessment {}
}

/// Admitted fixed-material-axis criteria supported by the native solid owner.
/// Sealed so source matching and criterion identity cannot be omitted by an
/// arbitrary callback. Both variants reuse the original material arithmetic.
pub trait OrthotropicStrengthCriterion:
    sealed::Criterion + Clone + core::fmt::Debug + PartialEq
{
    type Assessment: OrthotropicStrengthObservation;
    fn resolved(&self) -> &ResolvedMaterialStatePoint;
    fn assess(&self, stress: &Voigt) -> Result<Self::Assessment, MaterialStatePointError>;
    #[doc(hidden)]
    fn identity_domain() -> &'static str;
    #[doc(hidden)]
    fn update_identity(&self, identity: &mut DomainHasher);
}

/// Common nonnegative utilization without erasing criterion-specific outputs.
pub trait OrthotropicStrengthObservation:
    sealed::Assessment + Clone + core::fmt::Debug + PartialEq
{
    fn utilization(&self) -> f64;
    fn reaches_limit(&self) -> bool {
        self.utilization() >= 1.0
    }
}

impl sealed::Criterion for OrthotropicTensileLimit {}
impl OrthotropicStrengthCriterion for OrthotropicTensileLimit {
    type Assessment = OrthotropicTensileAssessment;
    fn resolved(&self) -> &ResolvedMaterialStatePoint {
        self.resolved()
    }
    fn assess(&self, stress: &Voigt) -> Result<Self::Assessment, MaterialStatePointError> {
        self.assess_material_stress(stress)
    }
    fn identity_domain() -> &'static str {
        "org.frankensim.fs-solid.orthotropic-tensile-assessment.v1"
    }
    fn update_identity(&self, identity: &mut DomainHasher) {
        for name in self.limit_properties() {
            hash_name(identity, name);
        }
        for value in self.limits_pa() {
            identity.update(&value.to_bits().to_le_bytes());
        }
    }
}
impl sealed::Assessment for OrthotropicTensileAssessment {}
impl OrthotropicStrengthObservation for OrthotropicTensileAssessment {
    fn utilization(&self) -> f64 {
        self.utilization()
    }
}

impl sealed::Criterion for OrthotropicTsaiWu {}
impl OrthotropicStrengthCriterion for OrthotropicTsaiWu {
    type Assessment = OrthotropicTsaiWuAssessment;
    fn resolved(&self) -> &ResolvedMaterialStatePoint {
        self.resolved()
    }
    fn assess(&self, stress: &Voigt) -> Result<Self::Assessment, MaterialStatePointError> {
        self.assess_material_stress(stress)
    }
    fn identity_domain() -> &'static str {
        "org.frankensim.fs-solid.orthotropic-tsai-wu-assessment.v1"
    }
    fn update_identity(&self, identity: &mut DomainHasher) {
        for name in self
            .strength_properties()
            .into_iter()
            .chain(self.interaction_properties())
        {
            hash_name(identity, name);
        }
        for value in self
            .strengths_pa()
            .into_iter()
            .chain(self.normalized_interactions())
        {
            identity.update(&value.to_bits().to_le_bytes());
        }
    }
}
impl sealed::Assessment for OrthotropicTsaiWuAssessment {}
impl OrthotropicStrengthObservation for OrthotropicTsaiWuAssessment {
    fn utilization(&self) -> f64 {
        self.utilization()
    }
}

fn hash_name(identity: &mut DomainHasher, name: &str) {
    identity.update(&(name.len() as u64).to_le_bytes());
    identity.update(name.as_bytes());
}
