//! Explicit external work for held-coordinate angular temperature continuation.
//!
//! Actual finite-cell temperatures select the original bulk material claims.
//! The independently prescribed interface temperature selects its own claims.
//! Signed normal and angular storage changes belong to a named external
//! parameter-work source, never to caloric heat or inferred contact entropy.

use super::JkrSnapHeatError;
use fs_contact::normal_patch::adhesive::{
    rolling::SourceJkrRollingError,
    temperature::{AdhesiveTemperatureError, AdhesiveTemperatures},
    torsion::SourceJkrTorsionError,
};
use fs_exec::CancelGate;
use fs_mbd::tangential_friction::adhesive::{
    rolling::RigidAdhesiveRollingError, torsion::RigidAdhesiveTorsionError,
};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};

#[derive(Clone, Copy, Debug)]
pub struct JkrAngularTemperatureHeatControls<'a> {
    /// Physical driver responsible for signed coefficient work, including return.
    pub parameter_work_source: &'a str,
    /// Finite nonnegative bound on the absolute total normal-plus-angular work.
    pub maximum_absolute_parameter_work_j: f64,
}

#[derive(Debug)]
pub enum JkrAngularTemperatureHeatError {
    Invalid(&'static str),
    HistoryMismatch,
    ParameterWorkLimit { work_j: f64, maximum_j: f64 },
    Heat(JkrSnapHeatError),
    Temperature(AdhesiveTemperatureError),
    Rolling(SourceJkrRollingError),
    Torsion(SourceJkrTorsionError),
    NativeRolling(RigidAdhesiveRollingError),
    NativeTorsion(RigidAdhesiveTorsionError),
}
impl From<JkrSnapHeatError> for JkrAngularTemperatureHeatError {
    fn from(error: JkrSnapHeatError) -> Self {
        Self::Heat(error)
    }
}
impl From<DepositionError> for JkrAngularTemperatureHeatError {
    fn from(error: DepositionError) -> Self {
        Self::Heat(error.into())
    }
}
impl From<AdhesiveTemperatureError> for JkrAngularTemperatureHeatError {
    fn from(error: AdhesiveTemperatureError) -> Self {
        Self::Temperature(error)
    }
}
impl From<SourceJkrRollingError> for JkrAngularTemperatureHeatError {
    fn from(error: SourceJkrRollingError) -> Self {
        Self::Rolling(error)
    }
}
impl From<SourceJkrTorsionError> for JkrAngularTemperatureHeatError {
    fn from(error: SourceJkrTorsionError) -> Self {
        Self::Torsion(error)
    }
}
impl From<RigidAdhesiveRollingError> for JkrAngularTemperatureHeatError {
    fn from(error: RigidAdhesiveRollingError) -> Self {
        Self::NativeRolling(error)
    }
}
impl From<RigidAdhesiveTorsionError> for JkrAngularTemperatureHeatError {
    fn from(error: RigidAdhesiveTorsionError) -> Self {
        Self::NativeTorsion(error)
    }
}
impl core::fmt::Display for JkrAngularTemperatureHeatError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adhesive angular temperature continuation: {self:?}")
    }
}
impl std::error::Error for JkrAngularTemperatureHeatError {}

impl JkrAngularTemperatureHeatControls<'_> {
    pub(super) fn validate(self) -> Result<(), JkrAngularTemperatureHeatError> {
        if self.parameter_work_source.trim().is_empty()
            || !self.maximum_absolute_parameter_work_j.is_finite()
            || self.maximum_absolute_parameter_work_j < 0.0
        {
            return Err(JkrAngularTemperatureHeatError::Invalid(
                "named external parameter-work source and finite nonnegative work bound required",
            ));
        }
        Ok(())
    }
    pub(super) fn admit_work(self, work_j: f64) -> Result<(), JkrAngularTemperatureHeatError> {
        if !work_j.is_finite() || work_j.abs() > self.maximum_absolute_parameter_work_j {
            return Err(JkrAngularTemperatureHeatError::ParameterWorkLimit {
                work_j,
                maximum_j: self.maximum_absolute_parameter_work_j,
            });
        }
        Ok(())
    }
}

pub(super) fn temperatures(
    thermal: &DepositionSimulation<'_>,
    cells: [CellId; 2],
    interface_k: f64,
) -> Result<AdhesiveTemperatures, JkrAngularTemperatureHeatError> {
    if !interface_k.is_finite() || interface_k <= 0.0 {
        return Err(JkrAngularTemperatureHeatError::Invalid(
            "positive finite prescribed interface temperature required",
        ));
    }
    Ok(AdhesiveTemperatures {
        surface_a_k: thermal.temperature_k(cells[0])?,
        surface_b_k: thermal.temperature_k(cells[1])?,
        interface_k,
    })
}

pub(super) fn check(
    gate: &CancelGate,
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<(), JkrAngularTemperatureHeatError> {
    if gate.is_requested() {
        return Err(DepositionError::Cancelled.into());
    }
    checkpoint()?;
    if gate.is_requested() {
        return Err(DepositionError::Cancelled.into());
    }
    Ok(())
}
