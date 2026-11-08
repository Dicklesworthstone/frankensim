//! Accepted thermal liquid states driving the original finite-bath owner.
//!
//! One well-mixed liquid cell retains its original caloric chart and exact tank
//! mass. Each native accepted substep refreshes the SAME pinned liquid source,
//! then the original tank solves its waterline at held displacement geometry.
//! Reference thermal exchange geometry remains the caller's lumped model.
//!
//! Expansion changes complete liquid gravitational storage, including the
//! empty-bath datum. This is separately supplied external parameter work, with
//! a whole-call sum-absolute allowance. It is never silently free, relabeled
//! caloric heat, or subtracted again as ambient pressure work. No closed-system
//! energy balance, moving-body impulse, mixing, sloshing or fluid inertia is claimed.

use super::{FiniteHydrostaticTank, FiniteTankDensityChange, HydrostaticError, finite};
use crate::hydrostatic::DisplacementMesh;
use fs_ga::Se3;
use fs_material::{
    fluid_transport::FluidTransportState,
    liquid::ResolvedLiquidState,
    liquid_temperature::{LiquidCaloricEndpoint, LiquidCaloricError, LiquidTemperatureLaw},
    phase::EquilibriumEnthalpyPhaseCurve,
};
use fs_thermal::enthalpy::{
    CellId, DepositedCell, DepositionError, DepositionObserver, DepositionSimulation, Event,
};
use std::collections::BTreeMap;

/// Explicit source agreement and external gravitational-work admission.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteTankThermalControls {
    /// Original chart/source density disagreement allowed at either endpoint.
    pub maximum_source_density_defect_kg_m3: f64,
    /// Sum of absolute COMPLETE gravitational changes over this advance [J].
    /// Cooling/heating cancellation cannot hide the amount of parameter work.
    pub maximum_absolute_parameter_work_j: f64,
}

/// Per-advance accounts for this liquid only; native heat retains its own audit.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FiniteTankThermalAudit {
    /// Actual m*(h_final-h_initial) of the bound original thermal cell [J].
    /// This may include internal exchange with other cells, not only external heat.
    pub caloric_storage_change_j: f64,
    /// Signed complete liquid gravitational-storage change [J].
    /// Positive means the external parameter owner supplied work.
    pub gravitational_parameter_work_j: f64,
    /// Sum of absolute per-update complete gravitational changes [J].
    pub absolute_gravitational_parameter_work_j: f64,
    pub source_updates: usize,
}

/// A refusal discards the thermal, liquid-source and hydrostatic candidate together.
#[derive(Debug)]
pub enum FiniteTankThermalError {
    Thermal(DepositionError),
    Caloric(LiquidCaloricError),
    Hydrostatic(HydrostaticError),
    Invalid(&'static str),
    ParameterWork {
        observed_absolute_j: f64,
        maximum_absolute_j: f64,
    },
}
impl From<DepositionError> for FiniteTankThermalError {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl From<LiquidCaloricError> for FiniteTankThermalError {
    fn from(error: LiquidCaloricError) -> Self {
        Self::Caloric(error)
    }
}
impl From<HydrostaticError> for FiniteTankThermalError {
    fn from(error: HydrostaticError) -> Self {
        Self::Hydrostatic(error)
    }
}
impl core::fmt::Display for FiniteTankThermalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "thermal finite tank: {self:?}")
    }
}
impl std::error::Error for FiniteTankThermalError {}

/// Native thermal observation with owned source endpoints and borrowed fixed geometry.
///
/// Use the returned observer with DepositionSimulation::advance_observed or an
/// existing observed implicit owner. Its clones contain independent mutable
/// history. The endpoint owns its source bundle; native tank views borrow that
/// endpoint without a self-reference or a hidden material query.
#[derive(Clone)]
pub struct FiniteTankThermalObserver<'t, 's, 'a> {
    template: &'t FiniteHydrostaticTank<'s, ResolvedLiquidState>,
    caloric: &'a EquilibriumEnthalpyPhaseCurve,
    reference_cell: DepositedCell,
    endpoint: LiquidCaloricEndpoint<'a>,
    transport: FluidTransportState,
    volume_m3: f64,
    mesh: &'t DisplacementMesh,
    world_from_body: Se3,
    controls: FiniteTankThermalControls,
    time_s: f64,
    initial_specific_enthalpy_j_kg: f64,
    audit: FiniteTankThermalAudit,
    last_density_change: Option<FiniteTankDensityChange>,
}

impl<'s> FiniteHydrostaticTank<'s, ResolvedLiquidState> {
    /// Observe one original liquid inventory through actual accepted heat steps.
    ///
    /// The cell must use the SAME caloric object, exact liquid mass and coherent
    /// reference liquid volume. The source law must resolve the tank's complete
    /// original liquid state at this initial enthalpy, retaining its selected
    /// claims and nonthermal context. Solid/coexisting phases refuse.
    ///
    /// Tank, displacement envelope and its pose stay fixed during the window.
    /// Native volume/source controls and overflow refusals remain authoritative.
    /// Heat-area/Biot and quasistatic/refinement assumptions stay explicit with
    /// their owners. This adapter neither adds heat nor moves the held body.
    pub fn thermal_density_observer<'t, 'a>(
        &'t self,
        thermal: &DepositionSimulation<'_>,
        liquid_cell: CellId,
        caloric: &'a EquilibriumEnthalpyPhaseCurve,
        law: &LiquidTemperatureLaw<'a>,
        mesh: &'t DisplacementMesh,
        world_from_body: Se3,
        controls: FiniteTankThermalControls,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<FiniteTankThermalObserver<'t, 's, 'a>, FiniteTankThermalError> {
        checkpoint()?;
        if ![
            controls.maximum_source_density_defect_kg_m3,
            controls.maximum_absolute_parameter_work_j,
        ]
        .iter()
        .all(|v| v.is_finite() && *v >= 0.0)
        {
            return Err(FiniteTankThermalError::Invalid(
                "finite nonnegative source/work allowances required",
            ));
        }
        let cell = thermal
            .cells()
            .get(&liquid_cell)
            .ok_or(FiniteTankThermalError::Invalid(
                "liquid thermal cell is not active",
            ))?;
        let volume = self.liquid_volume_m3();
        if cell.mass_kg != self.liquid_mass_kg()
            || !thermal.uses_material_law(cell.material_id, caloric)
            || (cell.geometry.volume_m3 - volume).abs()
                > 128.0 * f64::EPSILON * cell.geometry.volume_m3.max(volume)
        {
            return Err(FiniteTankThermalError::Invalid(
                "thermal and tank inventories need the same mass, original caloric and reference volume",
            ));
        }
        let state = caloric
            .state_at_specific_enthalpy(cell.specific_h_j_kg)
            .map_err(LiquidCaloricError::from)?;
        let endpoint = law.prepare_caloric_state_endpoint(
            caloric,
            state,
            state,
            controls.maximum_source_density_defect_kg_m3,
        )?;
        if endpoint.law().reference_state() != self.source() {
            return Err(FiniteTankThermalError::Invalid(
                "thermal law does not retain the tank's original complete liquid source",
            ));
        }
        hydrostatic_call(&mut checkpoint, |cancelled| {
            self.fill_at_pose(mesh, world_from_body, cancelled)
        })?;
        checkpoint()?;
        Ok(FiniteTankThermalObserver {
            template: self,
            caloric,
            reference_cell: cell.clone(),
            endpoint,
            transport: self.transport(),
            volume_m3: volume,
            mesh,
            world_from_body,
            controls,
            time_s: thermal.time_s(),
            initial_specific_enthalpy_j_kg: cell.specific_h_j_kg,
            audit: FiniteTankThermalAudit::default(),
            last_density_change: None,
        })
    }
}

impl FiniteTankThermalObserver<'_, '_, '_> {
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn audit(&self) -> FiniteTankThermalAudit {
        self.audit
    }
    pub const fn endpoint(&self) -> &LiquidCaloricEndpoint<'_> {
        &self.endpoint
    }
    pub const fn last_density_change(&self) -> Option<&FiniteTankDensityChange> {
        self.last_density_change.as_ref()
    }

    /// Continue original loads/steps using the accepted source and exact inventory.
    /// Reuses the original validated tank and the retained transport snapshot.
    /// No source is queried, volume solve performed, or new mass reconstructed.
    pub fn current_tank(
        &self,
    ) -> Result<FiniteHydrostaticTank<'_, ResolvedLiquidState>, HydrostaticError> {
        self.template.with_snapshot_inventory(
            self.endpoint.law().reference_state(),
            self.transport,
            self.volume_m3,
            self.template.liquid_mass_kg(),
        )
    }

    fn admit_cell(
        &self,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<f64, FiniteTankThermalError> {
        let old = &self.reference_cell;
        let cell = thermal
            .cells()
            .get(&old.id)
            .ok_or(FiniteTankThermalError::Invalid(
                "original liquid cell disappeared",
            ))?;
        if cell.mass_kg != old.mass_kg
            || cell.material_id != old.material_id
            || cell.body_id != old.body_id
            || cell.segment_id != old.segment_id
            || cell.geometry != old.geometry
            || !thermal.uses_material_law(cell.material_id, self.caloric)
        {
            return Err(FiniteTankThermalError::Invalid(
                "original liquid mass, caloric identity or reference geometry changed",
            ));
        }
        Ok(cell.specific_h_j_kg)
    }

    fn refresh(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), FiniteTankThermalError> {
        checkpoint()?;
        let h = self.admit_cell(thermal)?;
        if h.to_bits() == self.endpoint.state().specific_enthalpy_j_kg().to_bits() {
            checkpoint()?;
            self.time_s = thermal.time_s();
            return Ok(());
        }
        let state = self
            .caloric
            .state_at_specific_enthalpy(h)
            .map_err(LiquidCaloricError::from)?;
        let endpoint = self.endpoint.law().prepare_caloric_state_endpoint(
            self.caloric,
            self.endpoint.state(),
            state,
            self.controls.maximum_source_density_defect_kg_m3,
        )?;
        let (transport, volume, change) = {
            let previous = self.current_tank()?;
            let (next, change) = hydrostatic_call(checkpoint, |cancelled| {
                previous.rebind_source_at_pose(
                    endpoint.law().reference_state(),
                    self.mesh,
                    self.world_from_body,
                    cancelled,
                )
            })?;
            (next.transport(), next.liquid_volume_m3(), change)
        };
        let work = finite(
            self.audit.gravitational_parameter_work_j + change.liquid_gravitational_energy_change_j,
        )?;
        let absolute = finite(
            self.audit.absolute_gravitational_parameter_work_j
                + change.liquid_gravitational_energy_change_j.abs(),
        )?;
        if absolute > self.controls.maximum_absolute_parameter_work_j {
            return Err(FiniteTankThermalError::ParameterWork {
                observed_absolute_j: absolute,
                maximum_absolute_j: self.controls.maximum_absolute_parameter_work_j,
            });
        }
        let caloric_change =
            finite(self.reference_cell.mass_kg * finite(h - self.initial_specific_enthalpy_j_kg)?)?;
        let updates = self
            .audit
            .source_updates
            .checked_add(1)
            .ok_or(DepositionError::Arithmetic)?;
        checkpoint()?;
        self.endpoint = endpoint;
        self.transport = transport;
        self.volume_m3 = volume;
        self.time_s = thermal.time_s();
        self.audit = FiniteTankThermalAudit {
            caloric_storage_change_j: caloric_change,
            gravitational_parameter_work_j: work,
            absolute_gravitational_parameter_work_j: absolute,
            source_updates: updates,
        };
        self.last_density_change = Some(change);
        Ok(())
    }
}

impl DepositionObserver for FiniteTankThermalObserver<'_, '_, '_> {
    type Error = FiniteTankThermalError;

    fn begin(&mut self) -> Result<(), Self::Error> {
        self.initial_specific_enthalpy_j_kg = self.endpoint.state().specific_enthalpy_j_kg();
        self.audit = FiniteTankThermalAudit::default();
        self.last_density_change = None;
        Ok(())
    }

    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        checkpoint()?;
        if thermal.time_s() != self.time_s
            || self.admit_cell(thermal)?.to_bits()
                != self.endpoint.state().specific_enthalpy_j_kg().to_bits()
        {
            return Err(FiniteTankThermalError::Invalid(
                "thermal and finite-tank source epochs differ",
            ));
        }
        Ok(BTreeMap::new())
    }

    fn event(
        &mut self,
        _event: &Event,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), Self::Error> {
        if thermal.time_s() != self.time_s {
            return Err(FiniteTankThermalError::Invalid(
                "finite-tank event skipped an accepted thermal epoch",
            ));
        }
        self.refresh(thermal, &mut || Ok(()))
    }

    fn step(
        &mut self,
        start_s: f64,
        _duration_s: f64,
        _old_temperatures_k: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        if start_s != self.time_s || thermal.time_s() <= start_s {
            return Err(FiniteTankThermalError::Invalid(
                "finite-tank step must continue its accepted thermal epoch",
            ));
        }
        self.refresh(thermal, checkpoint)
    }
}

fn hydrostatic_call<T>(
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    call: impl FnOnce(&mut dyn FnMut() -> bool) -> Result<T, HydrostaticError>,
) -> Result<T, FiniteTankThermalError> {
    let mut stopped = None;
    let result = call(&mut || {
        if stopped.is_some() {
            return true;
        }
        match checkpoint() {
            Ok(()) => false,
            Err(error) => {
                stopped = Some(error);
                true
            }
        }
    });
    if let Some(error) = stopped {
        return Err(error.into());
    }
    Ok(result?)
}
