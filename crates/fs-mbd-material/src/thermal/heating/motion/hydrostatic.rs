//! Continuing buoyancy under source-resolved free thermal expansion.
//!
//! The existing clipped-volume hydrostatic owner supplies every force and
//! pressure potential. The existing thermal scheduler and momentum-preserving
//! configuration transfer advance the same material inventory. Pressure work
//! caused by changing the displacement envelope is explicit external parameter
//! work; it is not silently removed from or deposited into caloric heat.

use super::*;
use crate::hydrostatic::{
    DisplacementMesh, HydrostaticLoad, HydrostaticMedium, HydrostaticStepControls,
    MAX_DISPLACEMENT_CELLS,
};
use fs_material::fluid_transport::FluidTransport;
use std::collections::BTreeSet;

/// A non-overlapping part of the original fluid-excluding envelope, attached
/// to one named material region. Sealed voids may be included. The complete
/// partition and compatibility of different regions' free expansions remain
/// the caller's geometry obligations, as in the original hydrostatic owner.
#[derive(Clone, Debug)]
pub struct ThermalDisplacementRegion {
    pub region: String,
    pub envelope: DisplacementMesh,
}

/// One accepted thermal, mechanical and prescribed-fluid accounting window.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HeatedHydrostaticAudit {
    pub thermal: EnergyAudit,
    pub mechanical_before: DynamicsDiagnostics,
    pub mechanical_after: DynamicsDiagnostics,
    pub load_before: HydrostaticLoad,
    pub load_after: HydrostaticLoad,
    /// Signed inertia-remap plus fluid-displacement configuration work [J].
    pub parameter_work_j: f64,
    /// Absolute contributions are summed BEFORE cancellation between accounts.
    pub absolute_parameter_work_j: f64,
    /// The pressure-potential component of parameter_work_j [J].
    pub displacement_parameter_work_j: f64,
    /// Native fixed-configuration K + U_gravity + U_pressure defect [J].
    pub numerical_energy_change_j: f64,
    pub absolute_numerical_energy_change_j: f64,
    /// Thermal + all represented mechanical/fluid storage minus external
    /// heat, configuration work and measured numerical change [J].
    pub combined_residual_j: f64,
}

impl<'a> HeatedRigidBody<'a> {
    /// Continue the same heated body under gravity and hydrostatic buoyancy.
    ///
    /// Each original envelope vertex follows its own source region's accepted
    /// free strain and fixed point. Geometry is rebuilt FROM THE REFERENCE,
    /// never compounded, and the body's mass remains unchanged. The medium's
    /// source density, prescribed free surface and gravity are frozen during
    /// this advance; its gravity must match the original mechanical integrator.
    /// Continue with the same medium and reference envelope. Replacing those
    /// inputs is a model change whose reservoir/exchange work is not inferred.
    ///
    /// Native kick-drift motion precedes each thermal configuration update.
    /// This first-order split requires timestep refinement. No waves, fluid
    /// inertia/drag, moving reservoir, flooding, phase change, expansion-induced
    /// fluid kinetic energy or reciprocal thermoelastic heat is introduced.
    /// The existing absolute work/error budgets include pressure-potential
    /// contributions. Any late source, budget or cancellation failure discards
    /// the entire proposed thermal and mechanical history.
    pub fn advance_hydrostatic<F: FluidTransport + ?Sized>(
        &self,
        medium: &HydrostaticMedium<'_, F>,
        envelope: &[ThermalDisplacementRegion],
        target_s: f64,
        config: StepConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, HeatedHydrostaticAudit), HeatedMotionError> {
        checkpoint()?;
        if medium.gravity() != self.motion.integrator.gravity() {
            return Err(HeatedMotionError::Invalid(
                "hydrostatic and mechanical gravity must agree",
            ));
        }
        let mut parts = Vec::new();
        let mut keys = BTreeSet::new();
        let mut count = 0usize;
        for part in envelope {
            checkpoint()?;
            let index = self
                .reference()
                .regions()
                .binary_search_by(|r| r.key().cmp(part.region.as_str()))
                .map_err(|_| {
                    HeatedMotionError::Invalid("displacement region is not in the original body")
                })?;
            if !keys.insert(index) {
                return Err(HeatedMotionError::Invalid("duplicate displacement region"));
            }
            count = count
                .checked_add(part.envelope.cells().len())
                .ok_or(HeatedMotionError::Budget("displacement cell count"))?;
            if count > MAX_DISPLACEMENT_CELLS {
                return Err(HeatedMotionError::Budget("displacement cell count"));
            }
            parts.push((index, &part.envelope));
        }
        if parts.is_empty() {
            return Err(HeatedMotionError::Invalid("empty displacement envelope"));
        }
        parts.sort_by_key(|&(index, _)| index);
        let observer = HydrostaticExpansion {
            motion: self.motion.clone(),
            medium,
            parts,
            displacement_work_j: 0.0,
        };
        let load_before = observer.load(&mut checkpoint)?;
        let (thermal, observer, thermal_audit) =
            self.thermal
                .advance_observed(&observer, target_s, config, &mut checkpoint)?;
        let load_after = observer.load(&mut checkpoint)?;
        let work = observer.motion.work;
        let displacement_parameter_work_j = observer.displacement_work_j;
        let next = Self {
            thermal,
            motion: observer.motion,
        };
        let mechanical_before = self.diagnostics()?;
        let mechanical_after = next.diagnostics()?;
        let combined_residual_j = finite(
            thermal_audit.residual_j
                + mechanical_change(mechanical_before, mechanical_after)?
                + (load_after.pressure_potential_j() - load_before.pressure_potential_j())
                - work.parameter_j
                - work.numerical_j,
        )?;
        checkpoint()?;
        Ok((
            next,
            HeatedHydrostaticAudit {
                thermal: thermal_audit,
                mechanical_before,
                mechanical_after,
                load_before,
                load_after,
                parameter_work_j: work.parameter_j,
                absolute_parameter_work_j: work.absolute_parameter_j,
                displacement_parameter_work_j,
                numerical_energy_change_j: work.numerical_j,
                absolute_numerical_energy_change_j: work.absolute_numerical_j,
                combined_residual_j,
            },
        ))
    }
}

struct HydrostaticExpansion<'a, 'm, 'f, F: FluidTransport + ?Sized> {
    motion: RigidExpansion<'a>,
    medium: &'m HydrostaticMedium<'f, F>,
    parts: Vec<(usize, &'m DisplacementMesh)>,
    displacement_work_j: f64,
}
impl<F: FluidTransport + ?Sized> Clone for HydrostaticExpansion<'_, '_, '_, F> {
    fn clone(&self) -> Self {
        Self {
            motion: self.motion.clone(),
            medium: self.medium,
            parts: self.parts.clone(),
            displacement_work_j: self.displacement_work_j,
        }
    }
}
impl<F: FluidTransport + ?Sized> HydrostaticExpansion<'_, '_, '_, F> {
    fn mesh(
        &self,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<DisplacementMesh, HeatedMotionError> {
        let mut cells = Vec::new();
        for &(index, reference) in &self.parts {
            let region = &self.motion.expansion.configuration.regions()[index];
            let strain = region.expansion().free_linear_strain();
            for cell in reference.cells() {
                checkpoint()?;
                // Same reference-coordinate map as the original mass owner;
                // an exact zero strain retains the original geometry bits.
                cells.push(cell.map(|point| {
                    if strain == 0.0 {
                        point
                    } else {
                        point + (point - region.fixed_point_m()).scale(strain)
                    }
                }));
            }
        }
        checked(checkpoint, |cancelled| {
            DisplacementMesh::new(cells, MAX_DISPLACEMENT_CELLS, cancelled)
        })
    }
    fn load(
        &self,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<HydrostaticLoad, HeatedMotionError> {
        let mesh = self.mesh(checkpoint)?;
        checked(checkpoint, |cancelled| {
            self.medium
                .load_frame(&mesh, &self.motion.principal, self.motion.state, cancelled)
        })
    }
}
impl<F: FluidTransport + ?Sized> DepositionObserver for HydrostaticExpansion<'_, '_, '_, F> {
    type Error = HeatedMotionError;
    fn begin(&mut self) -> Result<(), Self::Error> {
        self.displacement_work_j = 0.0;
        self.motion.begin()
    }
    fn event(
        &mut self,
        event: &Event,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), Self::Error> {
        self.motion.event(event, thermal)
    }
    fn limit_timestep_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        proposed_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Self::Error> {
        self.motion
            .limit_timestep_s(thermal, proposed_s, checkpoint)
    }
    fn step(
        &mut self,
        start_s: f64,
        dt_s: f64,
        old: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        if start_s != self.motion.time_s {
            return Err(HeatedMotionError::Invalid(
                "hydrostatic and thermal clocks disagree",
            ));
        }
        let mesh = self.mesh(checkpoint)?;
        let step = checked(checkpoint, |cancelled| {
            self.medium.step_frame(
                &mesh,
                &self.motion.principal,
                self.motion.state,
                HydrostaticStepControls {
                    duration_s: dt_s,
                    maximum_energy_defect_j: self.motion.controls.maximum_numerical_energy_change_j,
                },
                cancelled,
            )
        })?;
        self.motion.finish_step(
            start_s,
            dt_s,
            old,
            thermal,
            step.dynamics.state_after,
            step.energy_defect_j,
            checkpoint,
        )?;
        let displacement_work_j = finite(
            self.load(checkpoint)?.pressure_potential_j() - step.load_after.pressure_potential_j(),
        )?;
        self.displacement_work_j = finite(self.displacement_work_j + displacement_work_j)?;
        self.motion.work.parameter_j = finite(self.motion.work.parameter_j + displacement_work_j)?;
        self.motion.work.absolute_parameter_j =
            finite(self.motion.work.absolute_parameter_j + displacement_work_j.abs())?;
        if self.motion.work.absolute_parameter_j
            > self.motion.controls.transfer.maximum_parameter_work_j
        {
            return Err(HeatedMotionError::Budget(
                "total absolute configuration work including displacement",
            ));
        }
        checkpoint()?;
        Ok(())
    }
}
