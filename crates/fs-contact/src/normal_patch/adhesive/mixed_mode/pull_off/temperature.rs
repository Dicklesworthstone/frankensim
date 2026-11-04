//! Original-source thermal continuation of a finite-stiffness adhesive fixture.
//!
//! Both apparatus stages stay fixed while the original bulk moduli, work of
//! adhesion and junction strengths change. The actual mixed-mode equilibrium
//! therefore moves; its contact AND both spring energies enter parameter work.
//! A vanished attached equilibrium or newly exceeded slip boundary refuses.
//! Detached history and previously released energy never become fresh adhesion.
//!
//! This endpoint parameter update does not resolve an intervening thermal path,
//! thermal expansion, snap dynamics, entropy production or a caloric heat law.

use super::{
    JkrPullOffBranch, JkrPullOffError, JkrPullOffFixture, JkrPullOffResponse, JkrPullOffState,
    SourceJkrMixedMode, Work,
};
use crate::normal_patch::adhesive::temperature::{AdhesiveTemperatureError, AdhesiveTemperatures};
use fs_blake3::ContentHash;
use fs_exec::CancelGate;
use fs_matdb::MaterialCard;

/// Failed original-source query or unavailable retained apparatus equilibrium.
#[derive(Clone, Debug, PartialEq)]
pub enum JkrPullOffTemperatureError {
    Source(AdhesiveTemperatureError),
    Fixture(JkrPullOffError),
}
impl From<AdhesiveTemperatureError> for JkrPullOffTemperatureError {
    fn from(error: AdhesiveTemperatureError) -> Self {
        Self::Source(error)
    }
}
impl From<JkrPullOffError> for JkrPullOffTemperatureError {
    fn from(error: JkrPullOffError) -> Self {
        Self::Fixture(error)
    }
}
impl core::fmt::Display for JkrPullOffTemperatureError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "mixed-mode adhesive fixture temperature refused: {self:?}"
        )
    }
}
impl std::error::Error for JkrPullOffTemperatureError {}

/// Immutable original-source successor at the exact accepted apparatus stages.
/// Adopting `response().state()` on `source()` is the caller's atomic commit.
#[derive(Clone, Debug, PartialEq)]
pub struct JkrPullOffThermalRefresh {
    original_fixture_identity: ContentHash,
    source: JkrPullOffFixture,
    before: JkrPullOffResponse,
    response: JkrPullOffResponse,
    parameter_work_j: f64,
}
impl JkrPullOffThermalRefresh {
    pub const fn original_fixture_identity(&self) -> ContentHash {
        self.original_fixture_identity
    }
    pub const fn source(&self) -> &JkrPullOffFixture {
        &self.source
    }
    pub fn into_source(self) -> JkrPullOffFixture {
        self.source
    }
    pub const fn before(&self) -> &JkrPullOffResponse {
        &self.before
    }
    pub const fn response(&self) -> &JkrPullOffResponse {
        &self.response
    }
    /// Change in contact-plus-apparatus potential, J, with either sign.
    pub const fn potential_change_j(&self) -> f64 {
        self.parameter_work_j
    }
    /// Both external stages are held, even though the elastic tip can move.
    pub const fn mechanical_work_j(&self) -> f64 {
        0.0
    }
    /// Signed supply by the prescribed parameter change; never friction heat.
    pub const fn parameter_work_j(&self) -> f64 {
        self.parameter_work_j
    }
}

impl JkrPullOffFixture {
    /// Refresh the original normal and junction sources at an accepted state.
    ///
    /// The ordered bulk cards must be the originals. Only their independent
    /// absolute temperatures and the interface temperature change; original
    /// claims, geometry, both spring stiffnesses, stages and solve limits stay
    /// fixed. Actual mixed-mode area, tip position and reactions are recomputed.
    /// The original maximum-radius bound is retained, not used as a new preload.
    ///
    /// Attached history remains attached only when the held stage lies on the
    /// successor's admitted full-stick branch. Lost stability or static capacity
    /// refuses instead of producing an unowned snap, slip or energy release.
    /// Detached history remains detached, including its accumulated release.
    /// A no-op preserves the complete fixture and canonical response exactly.
    ///
    /// Construction and equilibrium queries retain their original bounded
    /// evaluation budgets. No candidate changes an accepted caller-owned state.
    pub fn refresh_at_held_stages(
        &self,
        previous: JkrPullOffState,
        material_cards: [&MaterialCard; 2],
        temperatures: AdhesiveTemperatures,
        gate: &CancelGate,
    ) -> Result<JkrPullOffThermalRefresh, JkrPullOffTemperatureError> {
        if previous.identity != self.identity {
            return Err(JkrPullOffError::StateMismatch.into());
        }
        let mut budget = Work::new(self.limits.maximum_evaluations, gate);
        budget.poll()?;
        let mut before = match previous.branch {
            JkrPullOffBranch::Attached => self.attached(previous.normal_stage_m, &mut budget)?,
            JkrPullOffBranch::Detached => {
                self.detached(previous.normal_stage_m, previous.released_energy_j)?
            }
        };
        before.state.released_energy_j = previous.released_energy_j;

        let normal = self
            .contact
            .normal_source()
            .at_temperatures(material_cards, temperatures)?;
        // `at_temperatures` establishes original ordered normal-source
        // continuity; the retained shear owner supplies its pinned claims.
        let shear = self
            .shear
            .rebind_thermal_normal(normal.normal_state())
            .map_err(JkrPullOffError::from)?;
        if normal == *self.contact.normal_source() && shear == self.shear {
            budget.poll()?;
            return Ok(JkrPullOffThermalRefresh {
                original_fixture_identity: self.identity,
                source: self.clone(),
                before,
                response: before,
                parameter_work_j: 0.0,
            });
        }

        let contact = SourceJkrMixedMode::new(normal, self.contact.envelope)
            .map_err(JkrPullOffError::from)?;
        let source = Self::new(
            contact,
            shear,
            self.normal_stiffness,
            self.tangential_stiffness,
            self.tangential_stage,
            self.limits,
            gate,
        )?;
        let mut response = match previous.branch {
            JkrPullOffBranch::Attached => source.attached(previous.normal_stage_m, &mut budget)?,
            JkrPullOffBranch::Detached => {
                source.detached(previous.normal_stage_m, previous.released_energy_j)?
            }
        };
        response.state.released_energy_j = previous.released_energy_j;
        let parameter_work_j = response.total_stored_energy_j - before.total_stored_energy_j;
        super::finite(parameter_work_j)?;
        budget.poll()?;
        Ok(JkrPullOffThermalRefresh {
            original_fixture_identity: self.identity,
            source,
            before,
            response,
            parameter_work_j,
        })
    }
}
