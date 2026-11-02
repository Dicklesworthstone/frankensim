//! Source-owned temperature transfer on the existing three-component bond.
//!
//! The declared ordering changes temperature at the accepted local separation
//! and frame, then optionally advances the original spatial mechanical trial.
//! Geometry is retained bit-for-bit during the temperature stage. Its signed
//! parameter work is external to the mechanical ports; additional damage and
//! subsequent intrinsic fracture and viscosity remain separate physical losses.
//! Common and mixed-mode histories use the original material transfer owner.
//! This prescribed-temperature, staggered composition does not solve heat
//! transport, assign fracture heat, infer cure, or reconstruct an unseen path.

use super::super::{
    ObjectiveCohesiveConfig, ObjectiveCohesiveError, ObjectiveCohesivePoint, ObjectiveCohesiveStep,
    finite_scalar, validate_config,
};
use super::{
    ResolvedBilinearCohesive, SourceBoundObjectiveCohesiveError, SourceBoundObjectiveCohesivePoint,
    SourceBoundObjectiveCohesiveStep,
};
use fs_material::cohesive::temperature::CohesiveTemperatureStep;
use fs_mbd::{UnitQuaternion, Vec3};

/// Fixed-geometry successor and independently scaled extensive energy accounts.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceBoundObjectiveCohesiveTemperatureStep<'a> {
    point: SourceBoundObjectiveCohesivePoint<'a>,
    /// Original per-area constitutive history transfer; no mechanical time step.
    pub material: CohesiveTemperatureStep<'a>,
    /// Signed coefficient work supplied by the external temperature source [J].
    pub parameter_work_j: f64,
    /// Additional irreversible damage release at fixed local separation [J].
    pub damage_dissipated_j: f64,
    /// Actual new-minus-old extensive bond storage at unchanged geometry [J].
    pub stored_change_j: f64,
    /// Parameter work minus storage change and new damage release [J].
    pub energy_residual_j: f64,
}
impl<'a> SourceBoundObjectiveCohesiveTemperatureStep<'a> {
    /// Original geometry/history, now borrowing the explicitly admitted source.
    pub const fn point(&self) -> SourceBoundObjectiveCohesivePoint<'a> {
        self.point
    }
}

/// Complete temperature-plus-mechanical proposal under one old-to-new balance.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceBoundObjectiveCohesiveThermalStep<'a> {
    temperature: SourceBoundObjectiveCohesiveTemperatureStep<'a>,
    mechanical: SourceBoundObjectiveCohesiveStep<'a>,
    /// Actual final-minus-original extensive storage across both stages [J].
    pub total_stored_change_j: f64,
    /// Fixed-jump thermal damage plus subsequent intrinsic fracture loss [J].
    pub total_fracture_dissipated_j: f64,
    /// Physical bond-dashpot loss during the subsequent mechanical path [J].
    pub viscous_dissipated_j: f64,
    /// Spatial input work plus parameter work minus storage and both losses [J].
    pub energy_residual_j: f64,
}
impl<'a> SourceBoundObjectiveCohesiveThermalStep<'a> {
    /// Accepted candidate only; the original source, geometry and history remain
    /// untouched unless the external physical owner publishes this successor.
    pub const fn point(&self) -> SourceBoundObjectiveCohesivePoint<'a> {
        self.mechanical.point()
    }
    /// Original spatial force pair, balancing frame couple and exact path work.
    pub const fn numerical(&self) -> &ObjectiveCohesiveStep<'a> {
        self.mechanical.numerical()
    }
    /// Separate coefficient work and fixed-separation damage release.
    pub const fn temperature(&self) -> &SourceBoundObjectiveCohesiveTemperatureStep<'a> {
        &self.temperature
    }
}

impl SourceBoundObjectiveCohesivePoint<'_> {
    /// Transfer to an explicitly resolved temperature of the SAME ordered source
    /// and model family. Frame, reference offset, actual world separation and
    /// bonded area remain exactly unchanged. No force/couple displacement work
    /// or viscous rate is attributed to this instantaneous parameter stage.
    pub fn transfer_temperature<'b>(
        &self,
        next_source: &'b ResolvedBilinearCohesive,
        temperature_axis: &str,
        config: ObjectiveCohesiveConfig,
    ) -> Result<SourceBoundObjectiveCohesiveTemperatureStep<'b>, SourceBoundObjectiveCohesiveError>
    {
        validate_config(config)?;
        self.source
            .require_temperature_successor(next_source, temperature_axis)
            .map_err(SourceBoundObjectiveCohesiveError::MaterialAssignment)?;
        let material = self
            .point
            .material()
            .transfer_temperature(next_source.law(), config.material_energy_tolerance_j_m2)
            .map_err(ObjectiveCohesiveError::Material)?;
        // Reconstructing through `new` would unnecessarily round the accepted
        // world separation again. Only the explicitly transferred history moves.
        let numerical = ObjectiveCohesivePoint {
            material: material.point,
            ..self.point
        };
        let area = self.point.area_m2();
        let before = finite_scalar(area * self.point.material().stored_energy_j_m2())?;
        let after = finite_scalar(area * numerical.material().stored_energy_j_m2())?;
        let parameter_work = finite_scalar(area * material.parameter_work_j_m2)?;
        let released = finite_scalar(area * material.damage_dissipated_j_m2)?;
        let stored_change = finite_scalar(after - before)?;
        let defect = finite_scalar(parameter_work - stored_change - released)?;
        check_energy(defect, config.spatial_work_tolerance_j)?;
        Ok(SourceBoundObjectiveCohesiveTemperatureStep {
            point: SourceBoundObjectiveCohesivePoint {
                point: numerical,
                source: next_source,
            },
            material,
            parameter_work_j: parameter_work,
            damage_dissipated_j: released,
            stored_change_j: stored_change,
            energy_residual_j: defect,
        })
    }

    /// Stage the temperature transfer, then the original three-component bond
    /// path between actual world separations and material frames. The original
    /// spatial force/couple identity remains responsible for mechanical work.
    /// An independent complete energy gate includes the temperature stage before
    /// a successor is returned; a later rotation/rate/work refusal publishes
    /// neither the temperature change nor new damage. Subdivide combined paths
    /// for accuracy: this is an explicit first-order thermal composition.
    #[allow(clippy::too_many_arguments)]
    pub fn trial_with_temperature<'b>(
        &self,
        next_source: &'b ResolvedBilinearCohesive,
        temperature_axis: &str,
        frame: UnitQuaternion,
        separation_world_m: Vec3,
        dt_s: f64,
        config: ObjectiveCohesiveConfig,
    ) -> Result<SourceBoundObjectiveCohesiveThermalStep<'b>, SourceBoundObjectiveCohesiveError>
    {
        let temperature = self.transfer_temperature(next_source, temperature_axis, config)?;
        let mechanical = temperature
            .point()
            .trial(frame, separation_world_m, dt_s, config)?;
        let area = self.point.area_m2();
        let before = finite_scalar(area * self.point.material().stored_energy_j_m2())?;
        let after = finite_scalar(
            area * mechanical
                .point()
                .numerical()
                .material()
                .stored_energy_j_m2(),
        )?;
        let stored_change = finite_scalar(after - before)?;
        let fracture = finite_scalar(
            temperature.damage_dissipated_j
                + area
                    * mechanical
                        .numerical()
                        .material
                        .equilibrium
                        .energy
                        .dissipated,
        )?;
        let viscous =
            finite_scalar(area * mechanical.numerical().material.viscous_dissipated_j_m2)?;
        let defect = finite_scalar(
            mechanical.numerical().spatial_work_j + temperature.parameter_work_j
                - stored_change
                - fracture
                - viscous,
        )?;
        check_energy(defect, config.spatial_work_tolerance_j)?;
        Ok(SourceBoundObjectiveCohesiveThermalStep {
            temperature,
            mechanical,
            total_stored_change_j: stored_change,
            total_fracture_dissipated_j: fracture,
            viscous_dissipated_j: viscous,
            energy_residual_j: defect,
        })
    }
}

fn check_energy(defect_j: f64, tolerance_j: f64) -> Result<(), ObjectiveCohesiveError> {
    if defect_j.abs() > tolerance_j {
        Err(ObjectiveCohesiveError::WorkDefect {
            defect_j,
            tolerance_j,
        })
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests;
