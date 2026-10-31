//! Native material bushings, finite caloric storage and continuing motion.
//!
//! Each accepted midpoint damping loss enters its declared original cells once.
//! Their actual endpoint temperatures refresh the same interface claims before
//! the next mechanical step. Geometry and carrier mass remain prescribed.
//! Signed elastic parameter work is external work, not additional damping heat
//! or a complete reversible thermoelastic law. The thermal/mechanical split is
//! first order; callers retain timestep and source-band qualification duties.

use super::{AttachmentError, MaterialBushing, step_bushings_implicit};
use crate::articulated::ResolvedArticulated;
use crate::thermal::heating::{
    HeatingError, LinkHeating, ThermalRegion, initialize_articulated,
};
use fs_matdb::QueryPoint;
use fs_material::thermal::{
    interface::ThermalInterfacePoint, source::ThermalMaterialError,
};
use fs_mbd::articulated::motion::FreeFloatingMotionState;
use fs_mbd::articulated::motion::spatial_loads::{
    MAX_SPATIAL_MOTION_LOADS, SpatialMotionLoads,
    implicit::{ImplicitSpatialConfig, ImplicitSpatialStep},
};
use fs_thermal::enthalpy::{
    CellId, DepositionBudget, DepositionError, DepositionObserver, DepositionSimulation,
    EnergyAudit, Event, StepConfig, TimedEvent,
};
use fs_thermal::observers::{CoupledObservers, ObserverCouplingError};
use std::collections::BTreeMap;

/// Original ordered interface and explicit partition of its physical loss.
#[derive(Clone, Debug)]
pub struct BushingHeatBinding<'a> {
    pub bushing: MaterialBushing<'a>,
    pub point: ThermalInterfacePoint<'a>,
    /// Cells of the actual A/B material regions named by the attachment.
    pub cells: [CellId; 2],
    /// Caller-declared heat partition in [0,1]; no effusivity model is inferred.
    pub fraction_to_a: f64,
}

/// Per-advance physical budgets, in addition to native solver/thermal budgets.
#[derive(Clone, Copy, Debug)]
pub struct BushingHeatControls {
    pub mechanics: ImplicitSpatialConfig,
    /// Bound on the sum of absolute elastic parameter work over this advance.
    pub maximum_parameter_work_j: f64,
    /// Gate on each native work defect and the complete coupled energy residual.
    pub energy_tolerance_j: f64,
}

/// Thermal source failure, native mechanical refusal or a coupling defect.
#[derive(Debug)]
pub enum BushingHeatError {
    Thermal(DepositionError),
    Heating(HeatingError),
    Attachment(AttachmentError),
    Source(ThermalMaterialError),
    Invalid(&'static str),
}
impl From<DepositionError> for BushingHeatError {
    fn from(e: DepositionError) -> Self { Self::Thermal(e) }
}
impl From<HeatingError> for BushingHeatError {
    fn from(e: HeatingError) -> Self { Self::Heating(e) }
}
impl From<AttachmentError> for BushingHeatError {
    fn from(e: AttachmentError) -> Self { Self::Attachment(e) }
}
impl From<ThermalMaterialError> for BushingHeatError {
    fn from(e: ThermalMaterialError) -> Self { Self::Source(e) }
}
impl core::fmt::Display for BushingHeatError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "heated material bushings: {self:?}")
    }
}
impl std::error::Error for BushingHeatError {}

/// All energy terms cover the same returned mechanical and thermal histories.
pub struct BushingHeatReport {
    /// Raw native thermal audit; its external heat already includes damping.
    pub thermal: EnergyAudit,
    pub dissipated_heat_j: f64,
    /// Boundary/source and composed-observer heat, excluding bushing damping.
    /// A composed domain's own storage/work balance remains with that observer.
    pub external_heat_j: f64,
    /// Signed U_new(q_after)-U_old(q_after), supplied by the parameter controller.
    pub elastic_parameter_work_j: f64,
    pub absolute_parameter_work_j: f64,
    /// Native actuator and external-wrench/gravity work.
    pub mechanical_supplied_work_j: f64,
    /// Actual kinetic plus bushing-storage change, including parameter work.
    pub mechanical_energy_change_j: f64,
    /// Native work defects remain numerical error and are never deposited.
    pub mechanical_work_residual_j: f64,
    pub total_energy_residual_j: f64,
}

/// One restartable accepted history. A failed advance changes neither owner.
#[derive(Clone)]
pub struct HeatedBushings<'a> {
    reference: &'a ResolvedArticulated,
    thermal: DepositionSimulation<'a>,
    regions: Vec<Vec<ThermalRegion<'a>>>,
    bindings: Vec<BushingHeatBinding<'a>>,
    state: FreeFloatingMotionState,
}
impl<'a> HeatedBushings<'a> {
    /// Compile fixed reference masses and original caloric cards, then bind all
    /// sources to their actual paired initial temperatures. Expansion data in
    /// LinkHeating is not consumed: this owner declares rigid carrier geometry.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        reference: &'a ResolvedArticulated,
        links: Vec<LinkHeating<'a>>,
        mut bindings: Vec<BushingHeatBinding<'a>>,
        state: FreeFloatingMotionState,
        events: Vec<TimedEvent>,
        budget: DepositionBudget,
        config: StepConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, BushingHeatError> {
        checkpoint()?;
        if bindings.is_empty() || bindings.len() > MAX_SPATIAL_MOTION_LOADS {
            return Err(BushingHeatError::Invalid("bounded nonempty bushing inventory required"));
        }
        let (thermal, regions) = initialize_articulated(
            reference, links, events, budget, config, &mut checkpoint,
        )?;
        for binding in &mut bindings {
            checkpoint()?;
            if !core::ptr::eq(reference, binding.bushing.mechanism())
                || !binding.fraction_to_a.is_finite()
                || !(0.0..=1.0).contains(&binding.fraction_to_a)
            {
                return Err(BushingHeatError::Invalid("original mechanism and heat partition required"));
            }
            let point = actual_point(binding, &regions, &thermal)?;
            let original = binding.bushing.axes()[0].resolved();
            if original.axis_quantities() != point.axis_quantities()
                || original.query_point().len() != point.axes().len()
                || !original.query_point().iter().zip(point.axes()).all(|((a, x), (b, y))| {
                    a == b && x.to_bits() == y.to_bits()
                })
            {
                return Err(BushingHeatError::Invalid("initial interface must use the actual caloric state"));
            }
            let rebound = checked(&mut checkpoint, |cancel| {
                binding.bushing.rebind_at_state(binding.point.interface, &point, &state, 0.0, cancel)
            })?;
            binding.bushing = rebound.binding;
        }
        checkpoint()?;
        Ok(Self { reference, thermal, regions, bindings, state })
    }

    pub const fn thermal(&self) -> &DepositionSimulation<'a> { &self.thermal }
    pub const fn state(&self) -> &FreeFloatingMotionState { &self.state }
    pub fn bindings(&self) -> &[BushingHeatBinding<'a>] { &self.bindings }

    /// Advance on the original event/CFL-clipped thermal clock. All sourced
    /// bushings act in ONE native implicit mechanical solve per substep.
    /// Additional axial connections refuse because they have no heat routes.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_to(
        &self,
        target_s: f64,
        loads: SpatialMotionLoads<'_>,
        thermal_config: StepConfig,
        controls: BushingHeatControls,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, BushingHeatReport), BushingHeatError> {
        checkpoint()?;
        let observer = self.observer(loads, controls)?;
        let (thermal, observer, audit) = self.thermal.advance_observed(
            &observer, target_s, thermal_config, &mut checkpoint,
        )?;
        self.finish(thermal, observer, audit, &mut checkpoint)
    }

    /// Compose another physical source/transport owner on these SAME substeps.
    ///
    /// Its heat is added to native damping once; conductance overrides, event
    /// boundaries and step limits remain with the original thermal scheduler.
    /// The supplied observer runs first at each callback, and bushing sources
    /// and coefficient refresh see the same start/end thermal states. Both
    /// histories publish only after the coupled energy and final cancellation
    /// checks succeed. Clone must retain independent mutable physical history.
    ///
    /// The report treats non-bushing heat as external to this mechanical and
    /// caloric subsystem. If the observer owns electrical, chemical or other
    /// storage, its returned history must account for that supply; this method
    /// does not infer a total-energy law for an arbitrary composed domain.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_with<O: DepositionObserver>(
        &self,
        observer: &O,
        target_s: f64,
        loads: SpatialMotionLoads<'_>,
        thermal_config: StepConfig,
        controls: BushingHeatControls,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, O, BushingHeatReport), ObserverCouplingError<O::Error, BushingHeatError>> {
        checkpoint()?;
        let mechanics = self.observer(loads, controls).map_err(ObserverCouplingError::Second)?;
        let pair = CoupledObservers::new(observer.clone(), mechanics);
        let (thermal, pair, audit) = self.thermal.advance_observed(
            &pair, target_s, thermal_config, &mut checkpoint,
        )?;
        let (observer, mechanics) = pair.into_parts();
        let (next, report) = self.finish(thermal, mechanics, audit, &mut checkpoint)
            .map_err(|error| match error {
                BushingHeatError::Thermal(error) => ObserverCouplingError::Thermal(error),
                error => ObserverCouplingError::Second(error),
            })?;
        Ok((next, observer, report))
    }

    fn observer<'b>(
        &'b self,
        loads: SpatialMotionLoads<'b>,
        controls: BushingHeatControls,
    ) -> Result<EvolvingBushings<'a, 'b>, BushingHeatError> {
        if !loads.connections.is_empty()
            || !nonnegative(controls.maximum_parameter_work_j)
            || !nonnegative(controls.energy_tolerance_j)
        {
            return Err(BushingHeatError::Invalid("unrouted connections or invalid work budgets"));
        }
        Ok(EvolvingBushings {
            regions: &self.regions,
            bindings: self.bindings.clone(),
            state: self.state.clone(),
            time_s: self.thermal.time_s(),
            loads,
            controls,
            pending: None,
            heat: 0.0,
            parameter: 0.0,
            parameter_used: 0.0,
            supplied: 0.0,
            mechanical_change: 0.0,
            mechanical_residual: 0.0,
        })
    }

    fn finish(
        &self,
        thermal: DepositionSimulation<'a>,
        observer: EvolvingBushings<'a, '_>,
        audit: EnergyAudit,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, BushingHeatReport), BushingHeatError> {
        let external = finite(audit.external_heat_j - observer.heat)?;
        let residual = finite(observer.mechanical_change + audit.stored_energy_change_j
            - observer.supplied - observer.parameter - external)?;
        if residual.abs() > observer.controls.energy_tolerance_j || observer.pending.is_some() {
            return Err(BushingHeatError::Invalid("coupled bushing energy balance failed"));
        }
        checkpoint()?;
        let report = BushingHeatReport {
            thermal: audit,
            dissipated_heat_j: observer.heat,
            external_heat_j: external,
            elastic_parameter_work_j: observer.parameter,
            absolute_parameter_work_j: observer.parameter_used,
            mechanical_supplied_work_j: observer.supplied,
            mechanical_energy_change_j: observer.mechanical_change,
            mechanical_work_residual_j: observer.mechanical_residual,
            total_energy_residual_j: residual,
        };
        Ok((Self {
            reference: self.reference,
            thermal,
            regions: self.regions.clone(),
            bindings: observer.bindings,
            state: observer.state,
        }, report))
    }
}

#[derive(Clone)]
struct EvolvingBushings<'a, 'b> {
    regions: &'b [Vec<ThermalRegion<'a>>],
    bindings: Vec<BushingHeatBinding<'a>>,
    state: FreeFloatingMotionState,
    time_s: f64,
    loads: SpatialMotionLoads<'b>,
    controls: BushingHeatControls,
    pending: Option<(ImplicitSpatialStep, f64)>,
    heat: f64,
    parameter: f64,
    parameter_used: f64,
    supplied: f64,
    mechanical_change: f64,
    mechanical_residual: f64,
}
impl DepositionObserver for EvolvingBushings<'_, '_> {
    type Error = BushingHeatError;
    fn requires_step_interval(&self) -> bool { true }

    fn event(&mut self, event: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        if matches!(event, Event::Birth(_)) {
            return Err(BushingHeatError::Invalid("fixed mechanical inventory cannot accept births"));
        }
        Ok(())
    }

    fn step_sources(
        &mut self,
        _: &DepositionSimulation<'_>,
        start_s: f64,
        dt_s: f64,
        end_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        checkpoint()?;
        if start_s.to_bits() != self.time_s.to_bits() || self.pending.is_some() {
            return Err(BushingHeatError::Invalid("mechanical and thermal intervals differ"));
        }
        let refs = self.bindings.iter().map(|b| &b.bushing).collect::<Vec<_>>();
        let step = checked(checkpoint, |cancel| {
            step_bushings_implicit(&refs, &self.state, self.loads, dt_s, self.controls.mechanics, cancel)
        })?;
        if step.spatial.mechanical_work_residual.abs() > self.controls.energy_tolerance_j {
            return Err(BushingHeatError::Invalid("native mechanical work defect exceeds budget"));
        }
        let mut sources = BTreeMap::new();
        let mut heat = 0.0;
        for binding in &self.bindings {
            let power = checked(checkpoint, |cancel| {
                binding.bushing.evaluate(&step.midpoint_state, cancel)
            })?.dissipation_power;
            let loss = finite(dt_s * power)?;
            if loss < 0.0 {
                return Err(BushingHeatError::Invalid("negative native damping heat"));
            }
            let a = finite(loss * binding.fraction_to_a)?;
            for (cell, energy) in binding.cells.into_iter().zip([a, loss - a]) {
                let watts = finite(energy / dt_s)?;
                let value = sources.entry(cell).or_insert(0.0);
                *value = finite(*value + watts)?;
            }
            heat = finite(heat + loss)?;
        }
        let routed = finite(dt_s * sources.values().sum::<f64>())?;
        if (heat - step.spatial.dissipated_energy).abs() > self.controls.energy_tolerance_j
            || (routed - heat).abs() > self.controls.energy_tolerance_j
        {
            return Err(BushingHeatError::Invalid("native loss and explicit heat routes differ"));
        }
        self.pending = Some((step, end_s));
        Ok(sources)
    }

    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        checkpoint()?;
        let (step, end_s) = self.pending.take()
            .ok_or(BushingHeatError::Invalid("missing native mechanical candidate"))?;
        if thermal.time_s().to_bits() != end_s.to_bits() {
            return Err(BushingHeatError::Invalid("accepted mechanical and thermal clocks differ"));
        }
        let mut parameter = 0.0;
        for binding in &mut self.bindings {
            checkpoint()?;
            let point = actual_point(binding, self.regions, thermal)?;
            let allowance = self.controls.maximum_parameter_work_j - self.parameter_used;
            let rebound = checked(checkpoint, |cancel| {
                binding.bushing.rebind_at_state(
                    binding.point.interface, &point, &step.spatial.motion.state_after,
                    allowance, cancel,
                )
            })?;
            parameter = finite(parameter + rebound.parameter_work_j)?;
            self.parameter_used = finite(self.parameter_used + rebound.parameter_work_j.abs())?;
            if self.parameter_used > self.controls.maximum_parameter_work_j {
                return Err(BushingHeatError::Invalid("elastic parameter work exceeds advance budget"));
            }
            binding.bushing = rebound.binding;
        }
        let native = &step.spatial;
        self.heat = finite(self.heat + native.dissipated_energy)?;
        self.parameter = finite(self.parameter + parameter)?;
        self.supplied = finite(self.supplied + native.motion.external_work + native.motion.actuator_work)?;
        self.mechanical_change = finite(self.mechanical_change
            + native.motion.momentum_after.kinetic_energy - native.motion.momentum_before.kinetic_energy
            + native.stored_energy_after - native.stored_energy_before + parameter)?;
        self.mechanical_residual = finite(self.mechanical_residual + native.mechanical_work_residual)?;
        self.state = native.motion.state_after.clone();
        self.time_s = end_s;
        Ok(())
    }
}

/// Resolve original ordered regions before asking their source-coordinate owner.
fn actual_point(
    binding: &BushingHeatBinding<'_>,
    regions: &[Vec<ThermalRegion<'_>>],
    thermal: &DepositionSimulation<'_>,
) -> Result<QueryPoint, BushingHeatError> {
    let side = |side: usize| {
        let link = binding.bushing.geometry.frames[side].link;
        let region = regions.get(link)
            .and_then(|regions| regions.iter().find(|r| r.key == binding.bushing.regions[side]))
            .ok_or(BushingHeatError::Invalid("attachment lost its original thermal region"))?;
        if region.cell_id != binding.cells[side] {
            return Err(BushingHeatError::Invalid("heat cell differs from attachment material region"));
        }
        let cell = thermal.cells().get(&region.cell_id)
            .ok_or(BushingHeatError::Invalid("attachment thermal cell is absent"))?;
        let state = region.caloric.curve().state_at_specific_enthalpy(cell.specific_h_j_kg)
            .map_err(HeatingError::from)?;
        if !state.is_fully_solid() {
            return Err(BushingHeatError::Invalid("rigid bushing carrier must remain fully solid"));
        }
        Ok(state)
    };
    Ok(binding.point.at([side(0)?, side(1)?])?)
}

fn checked<T>(
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    operation: impl FnOnce(&mut dyn FnMut() -> bool) -> Result<T, AttachmentError>,
) -> Result<T, BushingHeatError> {
    let mut failure = None;
    let result = operation(&mut || match checkpoint() {
        Ok(()) => false,
        Err(error) => { failure = Some(error); true }
    });
    if let Some(error) = failure { return Err(error.into()); }
    Ok(result?)
}
fn nonnegative(value: f64) -> bool { value.is_finite() && value >= 0.0 }
fn finite(value: f64) -> Result<f64, BushingHeatError> {
    if value.is_finite() { Ok(value) }
    else { Err(BushingHeatError::Invalid("unrepresentable bushing energy or heat")) }
}
