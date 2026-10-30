//! Keep the complete original ordered interface attached during frame transport.

use super::{
    ObjectiveCohesiveConfig, ObjectiveCohesiveError, ObjectiveCohesivePoint, ObjectiveCohesiveStep,
};
use fs_matdb::MaterialStateId;
use fs_material::{
    cohesive::source::ResolvedBilinearCohesive, state_point::MaterialStatePointError,
};
use fs_mbd::{UnitQuaternion, Vec3};

/// Source identity or numerical transport refusal; no original history changes.
#[derive(Clone, Debug, PartialEq)]
pub enum SourceBoundObjectiveCohesiveError {
    MaterialAssignment(MaterialStatePointError),
    DifferentOriginalLaw,
    Numerical(ObjectiveCohesiveError),
}
impl From<ObjectiveCohesiveError> for SourceBoundObjectiveCohesiveError {
    fn from(value: ObjectiveCohesiveError) -> Self {
        Self::Numerical(value)
    }
}
impl core::fmt::Display for SourceBoundObjectiveCohesiveError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source-bound objective cohesive trial refused: {self:?}")
    }
}
impl std::error::Error for SourceBoundObjectiveCohesiveError {}

/// Borrowed original receipt, law and history; no material reselection or clone.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SourceBoundObjectiveCohesivePoint<'a> {
    point: ObjectiveCohesivePoint<'a>,
    source: &'a ResolvedBilinearCohesive,
}

/// Numerical work ports and their still-bound material successor.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceBoundObjectiveCohesiveStep<'a> {
    numerical: ObjectiveCohesiveStep<'a>,
    source: &'a ResolvedBilinearCohesive,
}
impl<'a> SourceBoundObjectiveCohesiveStep<'a> {
    pub const fn numerical(&self) -> &ObjectiveCohesiveStep<'a> {
        &self.numerical
    }
    pub const fn point(&self) -> SourceBoundObjectiveCohesivePoint<'a> {
        SourceBoundObjectiveCohesivePoint {
            point: self.numerical.point,
            source: self.source,
        }
    }
}

impl<'a> SourceBoundObjectiveCohesivePoint<'a> {
    /// Retain an existing history only if it borrows this exact source's law.
    /// Equal coefficients from another source are insufficient. Receiving
    /// surface materials must match the original ordered interface assignment.
    pub fn new(
        point: ObjectiveCohesivePoint<'a>,
        source: &'a ResolvedBilinearCohesive,
        receiving_materials: [&MaterialStateId; 2],
    ) -> Result<Self, SourceBoundObjectiveCohesiveError> {
        source
            .require_surface_materials(receiving_materials[0], receiving_materials[1])
            .map_err(SourceBoundObjectiveCohesiveError::MaterialAssignment)?;
        if !core::ptr::eq(point.material().law(), source.law()) {
            return Err(SourceBoundObjectiveCohesiveError::DifferentOriginalLaw);
        }
        Ok(Self { point, source })
    }

    pub const fn numerical(&self) -> &ObjectiveCohesivePoint<'a> {
        &self.point
    }
    pub const fn source(&self) -> &'a ResolvedBilinearCohesive {
        self.source
    }

    /// Transport the frame without selecting coefficients or changing history.
    pub fn transport_frame(
        self,
        frame: UnitQuaternion,
    ) -> Result<Self, SourceBoundObjectiveCohesiveError> {
        Ok(Self {
            point: self.point.transport_frame(frame)?,
            ..self
        })
    }

    /// Physical trial under the original fixed external source conditions.
    pub fn trial(
        &self,
        frame: UnitQuaternion,
        separation_world_m: Vec3,
        dt_s: f64,
        config: ObjectiveCohesiveConfig,
    ) -> Result<SourceBoundObjectiveCohesiveStep<'a>, SourceBoundObjectiveCohesiveError> {
        Ok(SourceBoundObjectiveCohesiveStep {
            numerical: self.point.trial(frame, separation_world_m, dt_s, config)?,
            source: self.source,
        })
    }

    /// Physical local jump with independent frame transport, without projection.
    pub fn trial_local_jump(
        &self,
        frame: UnitQuaternion,
        jump_m: [f64; 3],
        dt_s: f64,
        config: ObjectiveCohesiveConfig,
    ) -> Result<SourceBoundObjectiveCohesiveStep<'a>, SourceBoundObjectiveCohesiveError> {
        Ok(SourceBoundObjectiveCohesiveStep {
            numerical: self.point.trial_local_jump(frame, jump_m, dt_s, config)?,
            source: self.source,
        })
    }

    /// Preserve the original source while evaluating the same objective trial
    /// and its analytic damaged-area compression measure [Pa].
    pub fn trial_with_damage_weighted_compression(
        &self,
        frame: UnitQuaternion,
        separation_world_m: Vec3,
        dt_s: f64,
        config: ObjectiveCohesiveConfig,
    ) -> Result<(SourceBoundObjectiveCohesiveStep<'a>, f64), SourceBoundObjectiveCohesiveError>
    {
        let (numerical, pressure) = self.point.trial_with_damage_weighted_compression(
            frame,
            separation_world_m,
            dt_s,
            config,
        )?;
        Ok((
            SourceBoundObjectiveCohesiveStep {
                numerical,
                source: self.source,
            },
            pressure,
        ))
    }
}

#[cfg(test)]
mod tests;
