//! Adhesive rolling memory on an original, attached JKR contact.
//!
//! Dominik and Tielens (1995), DOI 10.1080/01418619508243800, equation 25,
//! supplies the elastic rolling moment. The original normal owner supplies
//! the actual attached radius, reduced radius/modulus and WORK OF ADHESION W;
//! the shared tribology owner uses W = 2 gamma in the paper's convention.
//! Critical rolling displacement is a separate, explicitly authored parameter.
//!
//! This is a scalar fixed-axis, frozen-normal-state approximation. The selected
//! elastic/perfectly plastic continuation is explicit, not an inferred atomistic
//! law. Rolling storage is distinct from irreversible loss and normal JKR
//! storage. `advance_normal` explicitly follows a prescribed attached JKR path
//! at held rolling angle, retains the original memory and returns the signed
//! work needed to change rolling stiffness. It does not add reciprocal normal
//! loading or solve a relaxed mixed-mode contact area. Temperature, detachment
//! or a changing rolling axis need other history/work owners; rebuilding this
//! owner cannot adopt its predecessor's history. No net-load cutoff is applied.

use super::{JkrError, JkrResponse, SourceJkrSphereContact};
use fs_blake3::{ContentHash, DomainHasher};
use fs_math::dd::Dd;
use fs_tribo::rolling_loss::adhesive::{
    AdhesiveRollingContact, AdhesiveRollingError, AdhesiveRollingParameters, AdhesiveRollingState,
    AdhesiveRollingStep,
};
use std::sync::Arc;

/// Independent rolling calibration and small-shift approximation ceiling.
/// Neither field is inferred from bulk elasticity or work of adhesion.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrRollingControls {
    /// Elastic contact shift at irreversible rolling onset, m.
    pub critical_rolling_displacement_m: f64,
    /// Maximum critical shift divided by the actual attached contact radius.
    /// This authored approximation limit must lie strictly between zero and one.
    pub maximum_critical_displacement_to_contact_radius: f64,
}

/// Original normal-source or scalar rolling-history refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum SourceJkrRollingError {
    Normal(JkrError),
    Rolling(AdhesiveRollingError),
    /// Accepted memory belongs to another source, normal state or calibration.
    StateMismatch,
}
impl From<JkrError> for SourceJkrRollingError {
    fn from(error: JkrError) -> Self {
        Self::Normal(error)
    }
}
impl From<AdhesiveRollingError> for SourceJkrRollingError {
    fn from(error: AdhesiveRollingError) -> Self {
        Self::Rolling(error)
    }
}
impl core::fmt::Display for SourceJkrRollingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source-bound JKR rolling: {self:?}")
    }
}
impl std::error::Error for SourceJkrRollingError {}

/// Sealed scalar memory and its complete original normal/source identity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SourceJkrRollingState {
    identity: ContentHash,
    constitutive: AdhesiveRollingState,
}
impl SourceJkrRollingState {
    pub fn elastic_shift_m(&self) -> f64 {
        self.constitutive.elastic_shift_m()
    }
    pub fn accumulated_angle_rad(&self) -> f64 {
        self.constitutive.accumulated_angle_rad()
    }
}

/// Immutable trial with distinct reversible storage, rolling loss and work.
/// Only `state()` supplies the source-bound successor; public diagnostics
/// cannot be rewritten into an accepted physical history.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SourceJkrRollingStep {
    state: SourceJkrRollingState,
    rolling: AdhesiveRollingStep,
}
impl SourceJkrRollingStep {
    pub const fn state(&self) -> SourceJkrRollingState {
        self.state
    }
    pub const fn rolling(&self) -> &AdhesiveRollingStep {
        &self.rolling
    }
}

/// Original-source attached normal continuation with held rolling memory.
/// This is an immutable proposal; publishing source and state together commits it.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceJkrRollingNormalStep {
    source: SourceJkrRollingContact,
    state: SourceJkrRollingState,
    normal_contact_work_j: f64,
    rolling_storage_before_j: f64,
    rolling_storage_after_j: f64,
    rolling_parameter_work_j: f64,
    applied_work_j: f64,
}
impl SourceJkrRollingNormalStep {
    pub const fn source(&self) -> &SourceJkrRollingContact {
        &self.source
    }
    pub const fn state(&self) -> SourceJkrRollingState {
        self.state
    }
    /// Original integral F_JKR d(delta), with no branch jump or rolling work.
    pub const fn normal_contact_work_j(&self) -> f64 {
        self.normal_contact_work_j
    }
    pub const fn rolling_storage_before_j(&self) -> f64 {
        self.rolling_storage_before_j
    }
    pub const fn rolling_storage_after_j(&self) -> f64 {
        self.rolling_storage_after_j
    }
    /// Signed change of rolling spring storage at fixed elastic shift. The
    /// prescribed normal controller supplies or receives this work, never heat.
    pub const fn rolling_parameter_work_j(&self) -> f64 {
        self.rolling_parameter_work_j
    }
    /// Normal contact work plus signed rolling parameter work, J.
    pub const fn applied_work_j(&self) -> f64 {
        self.applied_work_j
    }
    pub fn into_parts(self) -> (SourceJkrRollingContact, SourceJkrRollingState) {
        (self.source, self.state)
    }
}

/// Frozen original JKR geometry plus one explicitly selected rolling law.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceJkrRollingContact {
    // Source bundles own material strings/vectors. A normal-area proposal must
    // share these immutable receipts, not deep-copy them at every time step.
    normal: Arc<SourceJkrSphereContact>,
    response: JkrResponse,
    law: AdhesiveRollingContact,
    identity: ContentHash,
    controls: JkrRollingControls,
}
impl SourceJkrRollingContact {
    /// Evaluate the original attached branch at signed approach, then retain
    /// that normal state throughout scalar rolling. Negative and zero normal
    /// resultants keep their real adhesive area. No response fields supplied
    /// by a caller are trusted as constitutive inputs.
    pub fn new(
        normal: SourceJkrSphereContact,
        indentation_m: f64,
        controls: JkrRollingControls,
    ) -> Result<Self, SourceJkrRollingError> {
        let response = normal.evaluate_attached(indentation_m)?;
        let law = AdhesiveRollingContact::new(AdhesiveRollingParameters {
            effective_radius_m: normal.effective_radius_m(),
            work_of_adhesion_j_m2: normal.normal_state().adhesion_energy_j_per_m2(),
            reduced_modulus_pa: normal.normal_state().elastic().reduced_modulus_pa(),
            contact_radius_m: response.contact_radius_m,
            critical_rolling_displacement_m: controls.critical_rolling_displacement_m,
            maximum_critical_displacement_to_contact_radius: controls
                .maximum_critical_displacement_to_contact_radius,
        })?;
        let identity = identity(&normal, indentation_m, controls);
        Ok(Self {
            normal: Arc::new(normal),
            response,
            law,
            identity,
            controls,
        })
    }

    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    pub fn normal_source(&self) -> &SourceJkrSphereContact {
        &self.normal
    }
    /// Original normal response. Its storage is NOT included in rolling storage.
    pub const fn normal_response(&self) -> &JkrResponse {
        &self.response
    }
    pub const fn law(&self) -> &AdhesiveRollingContact {
        &self.law
    }
    pub const fn controls(&self) -> JkrRollingControls {
        self.controls
    }
    pub fn initial_state(&self) -> SourceJkrRollingState {
        SourceJkrRollingState {
            identity: self.identity,
            constitutive: self.law.initial_state(),
        }
    }
    pub fn stored_energy_j(
        &self,
        state: &SourceJkrRollingState,
    ) -> Result<f64, SourceJkrRollingError> {
        self.check(state)?;
        Ok(self.law.stored_energy_j(&state.constitutive)?)
    }
    /// Exact monotone scalar increment, including elastic unloading and yield
    /// reversal. Constant evaluation cost; no allocations or source re-query.
    /// Loss is unpartitioned: a thermal consumer must declare its destinations.
    pub fn advance(
        &self,
        state: &SourceJkrRollingState,
        angle_increment_rad: f64,
    ) -> Result<SourceJkrRollingStep, SourceJkrRollingError> {
        self.check(state)?;
        let rolling = self.law.advance(&state.constitutive, angle_increment_rad)?;
        Ok(SourceJkrRollingStep {
            state: SourceJkrRollingState {
                identity: self.identity,
                constitutive: rolling.state,
            },
            rolling,
        })
    }

    /// Prescribe one attached normal-approach segment at fixed rolling angle,
    /// source, calibration and frame. Split reversals into separate calls.
    /// Actual area changes stiffness without resetting elastic/plastic memory.
    /// Original JKR work and signed rolling-storage work stay separate; neither
    /// creates loss. A target below the attached fold refuses without release.
    /// No reciprocal normal loading, mixed-mode relaxation or thermal transfer
    /// is supplied. Bounded arithmetic shares the original immutable receipts.
    pub fn advance_normal(
        &self,
        state: &SourceJkrRollingState,
        indentation_m: f64,
    ) -> Result<SourceJkrRollingNormalStep, SourceJkrRollingError> {
        self.check(state)?;
        let response = self.normal.evaluate_attached(indentation_m)?;
        let refresh = self
            .law
            .refresh_contact_radius(&state.constitutive, response.contact_radius_m)?;
        let normal_contact_work_j = self
            .normal
            .reversible_work_between(self.response.indentation_m, indentation_m)?;
        let rolling_storage_before_j = refresh.storage_before_j();
        let rolling_storage_after_j = refresh.storage_after_j();
        let rolling_parameter_work_j = refresh.storage_change_j();
        let applied_work_j =
            (Dd::from_f64(normal_contact_work_j) + Dd::from_f64(rolling_parameter_work_j)).to_f64();
        if !applied_work_j.is_finite() {
            return Err(AdhesiveRollingError::Unrepresentable {
                field: "normal and rolling parameter work",
            }
            .into());
        }
        let (law, constitutive) = refresh.into_parts();
        let identity = identity(&self.normal, indentation_m, self.controls);
        Ok(SourceJkrRollingNormalStep {
            source: Self {
                normal: self.normal.clone(),
                response,
                law,
                identity,
                controls: self.controls,
            },
            state: SourceJkrRollingState {
                identity,
                constitutive,
            },
            normal_contact_work_j,
            rolling_storage_before_j,
            rolling_storage_after_j,
            rolling_parameter_work_j,
            applied_work_j,
        })
    }
    fn check(&self, state: &SourceJkrRollingState) -> Result<(), SourceJkrRollingError> {
        if state.identity != self.identity {
            return Err(SourceJkrRollingError::StateMismatch);
        }
        Ok(())
    }
}

fn identity(
    normal: &SourceJkrSphereContact,
    indentation_m: f64,
    controls: JkrRollingControls,
) -> ContentHash {
    let mut hash = DomainHasher::new("org.frankensim.fs-contact.jkr-rolling.v1");
    hash.update(normal.identity().as_bytes());
    for value in [
        indentation_m,
        controls.critical_rolling_displacement_m,
        controls.maximum_critical_displacement_to_contact_radius,
    ] {
        hash.update(&value.to_bits().to_le_bytes());
    }
    hash.finalize()
}
