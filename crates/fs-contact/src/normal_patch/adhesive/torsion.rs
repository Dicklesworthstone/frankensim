//! Bonded torsion on the original attached JKR area and two elastic bodies.
//!
//! The normal response supplies actual radius, including at tensile resultant
//! load. Both original shear moduli supply K_phi = 16 G_bar a^3/3, with
//! 1/G_bar = 1/G_a + 1/G_b. Critical twist is independent interface calibration;
//! neither work of adhesion nor a Coulomb coefficient determines it.
//!
//! This freezes the original mode-I JKR area. The caller explicitly limits the
//! omitted mode-III release at the largest stored twist:
//! G_III = d(K_phi phi^2/2)/d(pi a^2) = 4 G_bar a phi^2/pi.
//! G_III/W is an approximation-envelope screen, not a bound on radius error or
//! a solution of mixed-mode Griffith equilibrium, especially near a fold.
//! No plastic loss, normal storage or adhesion release is counted as heat here.

use super::{JkrError, JkrResponse, SourceJkrSphereContact};
use core::f64::consts::PI;
use fs_blake3::{ContentHash, DomainHasher};
use fs_tribo::adhesive::torsion::{
    AdhesiveTorsionContact, AdhesiveTorsionError, AdhesiveTorsionParameters, AdhesiveTorsionState,
    AdhesiveTorsionStep,
};

/// Independent irreversible threshold and explicit frozen-area approximation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrTorsionControls {
    pub critical_twist_angle_rad: f64,
    /// Critical rim displacement / contact radius; strictly between 0 and 1.
    pub maximum_critical_twist_angle_rad: f64,
    /// Largest admitted G_III/W at the critical elastic twist, between 0 and 1.
    /// This is an authored model ceiling, not an accuracy certificate.
    pub maximum_mode_iii_energy_release_to_adhesion: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum SourceJkrTorsionError {
    Normal(JkrError),
    Torsion(AdhesiveTorsionError),
    Invalid { field: &'static str },
    OutsideFrozenArea { ratio: f64, limit: f64 },
    StateMismatch,
}
impl From<JkrError> for SourceJkrTorsionError {
    fn from(error: JkrError) -> Self {
        Self::Normal(error)
    }
}
impl From<AdhesiveTorsionError> for SourceJkrTorsionError {
    fn from(error: AdhesiveTorsionError) -> Self {
        Self::Torsion(error)
    }
}
impl core::fmt::Display for SourceJkrTorsionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source-bound JKR torsion: {self:?}")
    }
}
impl std::error::Error for SourceJkrTorsionError {}

/// Sealed twist memory, including all normal source, calibration and limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SourceJkrTorsionState {
    identity: ContentHash,
    constitutive: AdhesiveTorsionState,
}
impl SourceJkrTorsionState {
    pub fn elastic_twist_rad(&self) -> f64 {
        self.constitutive.elastic_twist_rad()
    }
    pub fn accumulated_angle_rad(&self) -> f64 {
        self.constitutive.accumulated_angle_rad()
    }
}

/// Unmodifiable source-bound successor and its distinct storage/work/loss.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SourceJkrTorsionStep {
    state: SourceJkrTorsionState,
    torsion: AdhesiveTorsionStep,
}
impl SourceJkrTorsionStep {
    pub const fn state(&self) -> SourceJkrTorsionState {
        self.state
    }
    pub const fn torsion(&self) -> &AdhesiveTorsionStep {
        &self.torsion
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct SourceJkrTorsionContact {
    normal: SourceJkrSphereContact,
    response: JkrResponse,
    law: AdhesiveTorsionContact,
    identity: ContentHash,
    maximum_mode_iii_release: f64,
}
impl SourceJkrTorsionContact {
    /// Re-evaluate the original normal source; no caller-supplied radius or
    /// response fields enter stiffness. The selected plastic plateau is a
    /// lumped bonded-junction model, not the local torsional partial-slip law.
    pub fn new(
        normal: SourceJkrSphereContact,
        indentation_m: f64,
        controls: JkrTorsionControls,
    ) -> Result<Self, SourceJkrTorsionError> {
        let limit = controls.maximum_mode_iii_energy_release_to_adhesion;
        if !limit.is_finite() || limit <= 0.0 || limit >= 1.0 {
            return Err(SourceJkrTorsionError::Invalid {
                field: "mode-III release/adhesion ceiling",
            });
        }
        let response = normal.evaluate_attached(indentation_m)?;
        let elastic = normal.normal_state().elastic();
        let law = AdhesiveTorsionContact::new(AdhesiveTorsionParameters {
            contact_radius_m: response.contact_radius_m,
            shear_moduli_pa: [elastic.surface_a_state(), elastic.surface_b_state()]
                .map(|s| s.elastic_state().shear_modulus_pa()),
            critical_twist_angle_rad: controls.critical_twist_angle_rad,
            maximum_critical_twist_angle_rad: controls.maximum_critical_twist_angle_rad,
        })?;
        let angle = controls.critical_twist_angle_rad;
        let maximum_mode_iii_release = (4.0 / PI)
            * (law.reduced_shear_modulus_pa() * angle)
            * (response.contact_radius_m * angle);
        let ratio = maximum_mode_iii_release / normal.normal_state().adhesion_energy_j_per_m2();
        if !ratio.is_finite() || ratio < 0.0 {
            return Err(SourceJkrTorsionError::Invalid {
                field: "unrepresentable mode-III release",
            });
        }
        if ratio > limit {
            return Err(SourceJkrTorsionError::OutsideFrozenArea { ratio, limit });
        }
        let mut hash = DomainHasher::new("org.frankensim.fs-contact.jkr-torsion.v1");
        hash.update(normal.identity().as_bytes());
        for value in [
            indentation_m,
            controls.critical_twist_angle_rad,
            controls.maximum_critical_twist_angle_rad,
            limit,
        ] {
            hash.update(&value.to_bits().to_le_bytes());
        }
        Ok(Self {
            normal,
            response,
            law,
            identity: hash.finalize(),
            maximum_mode_iii_release,
        })
    }
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    pub const fn normal_source(&self) -> &SourceJkrSphereContact {
        &self.normal
    }
    /// Normal storage remains separate from the torsional storage below.
    pub const fn normal_response(&self) -> &JkrResponse {
        &self.response
    }
    pub const fn law(&self) -> &AdhesiveTorsionContact {
        &self.law
    }
    /// Omitted mode-III release evaluated at the critical elastic twist [J/m2].
    pub const fn maximum_mode_iii_release_j_m2(&self) -> f64 {
        self.maximum_mode_iii_release
    }
    pub fn initial_state(&self) -> SourceJkrTorsionState {
        SourceJkrTorsionState {
            identity: self.identity,
            constitutive: self.law.initial_state(),
        }
    }
    pub fn stored_energy_j(
        &self,
        state: &SourceJkrTorsionState,
    ) -> Result<f64, SourceJkrTorsionError> {
        self.check(state)?;
        Ok(self.law.stored_energy_j(&state.constitutive)?)
    }
    pub fn advance(
        &self,
        state: &SourceJkrTorsionState,
        angle_increment_rad: f64,
    ) -> Result<SourceJkrTorsionStep, SourceJkrTorsionError> {
        self.check(state)?;
        let torsion = self.law.advance(&state.constitutive, angle_increment_rad)?;
        Ok(SourceJkrTorsionStep {
            state: SourceJkrTorsionState {
                identity: self.identity,
                constitutive: torsion.state,
            },
            torsion,
        })
    }
    fn check(&self, state: &SourceJkrTorsionState) -> Result<(), SourceJkrTorsionError> {
        if state.identity != self.identity {
            return Err(SourceJkrTorsionError::StateMismatch);
        }
        Ok(())
    }
}
