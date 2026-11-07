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
//! storage. A changing normal area, temperature, detachment or rolling axis
//! needs a separate history/work transfer; rebuilding this owner cannot adopt
//! its predecessor's history. No net-normal-load friction cutoff is applied.

use super::{JkrError, JkrResponse, SourceJkrSphereContact};
use fs_blake3::{ContentHash, DomainHasher};
use fs_tribo::rolling_loss::adhesive::{
    AdhesiveRollingContact, AdhesiveRollingError, AdhesiveRollingParameters, AdhesiveRollingState,
    AdhesiveRollingStep,
};

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

/// Frozen original JKR geometry plus one explicitly selected rolling law.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceJkrRollingContact {
    normal: SourceJkrSphereContact,
    response: JkrResponse,
    law: AdhesiveRollingContact,
    identity: ContentHash,
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
        let mut hash = DomainHasher::new("org.frankensim.fs-contact.jkr-rolling.v1");
        hash.update(normal.identity().as_bytes());
        for value in [
            indentation_m,
            controls.critical_rolling_displacement_m,
            controls.maximum_critical_displacement_to_contact_radius,
        ] {
            hash.update(&value.to_bits().to_le_bytes());
        }
        Ok(Self {
            normal,
            response,
            law,
            identity: hash.finalize(),
        })
    }

    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    pub const fn normal_source(&self) -> &SourceJkrSphereContact {
        &self.normal
    }
    /// Original normal response. Its storage is NOT included in rolling storage.
    pub const fn normal_response(&self) -> &JkrResponse {
        &self.response
    }
    pub const fn law(&self) -> &AdhesiveRollingContact {
        &self.law
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
    fn check(&self, state: &SourceJkrRollingState) -> Result<(), SourceJkrRollingError> {
        if state.identity != self.identity {
            return Err(SourceJkrRollingError::StateMismatch);
        }
        Ok(())
    }
}
