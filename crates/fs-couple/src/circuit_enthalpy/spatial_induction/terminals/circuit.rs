//! Original electrode circuits, finite caloric receivers and physical losses.
//!
//! Nonlinear and hysteretic field histories share the native thermal residual
//! and these graph-source/heat admissions. Only authored graph resistors have
//! receivers here. Numerical interface reference resistors are not physical
//! bodies and cannot supply heat to the thermal state.

use super::super::{
    CircuitEnthalpyError, Error, MagneticThermalError, MaterialNetworkSample,
    SpatialInductionThermalReport, SpatialMagneticRegion, SpatialMagneticThermalNetwork, check,
    finite, invalid,
};
use crate::circuit_enthalpy::{
    ResistorCell,
    contacts::heat::{HeatNetworkError, HeatNetworkSample},
    material::CellResistor,
};
use fs_em::{
    spatial::vector_potential::induction::terminals::circuit::SpatialTerminalCircuitError,
    winding::{nonlinear::resistors::GraphResistorMaterial, ohmic::WindingGraphResistorHeat},
};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::{
    cell::RefCell,
    collections::{BTreeMap, BTreeSet},
};

pub(in crate::circuit_enthalpy::spatial_induction) use super::super::nonlinear::winding::{
    admit_initial_graph_sources, pinned_graph_sources,
};

/// Original field/material, native electrode-circuit or complete heat refusal.
#[derive(Debug)]
pub enum SpatialTerminalCircuitThermalError<E = HeatNetworkError> {
    Native(MagneticThermalError),
    Circuit(SpatialTerminalCircuitError),
    Heat(E),
}

impl<E> From<MagneticThermalError> for SpatialTerminalCircuitThermalError<E> {
    fn from(error: MagneticThermalError) -> Self {
        Self::Native(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for SpatialTerminalCircuitThermalError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spatial electrode circuit/finite heat: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for SpatialTerminalCircuitThermalError<E> {}

pub(in crate::circuit_enthalpy::spatial_induction) type TerminalError<E = HeatNetworkError> =
    CircuitEnthalpyError<SpatialTerminalCircuitThermalError<E>>;

pub(in crate::circuit_enthalpy::spatial_induction) fn electrical<E>(
    error: SpatialTerminalCircuitError,
) -> TerminalError<E> {
    match error {
        SpatialTerminalCircuitError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(SpatialTerminalCircuitThermalError::Circuit(other)),
    }
}

pub(in crate::circuit_enthalpy::spatial_induction) fn heat_error<E: From<HeatNetworkError>>(
    error: HeatNetworkError,
) -> TerminalError<E> {
    match error {
        HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => {
            CircuitEnthalpyError::Source(SpatialTerminalCircuitThermalError::Heat(other.into()))
        }
    }
}

/// Original magnetic bodies and explicit authored graph-resistor receivers.
pub struct SpatialTerminalCircuitThermalNetwork<'a> {
    pub magnetic: SpatialMagneticThermalNetwork<'a>,
    /// Exactly one receiver for each authored graph resistor. With material
    /// bindings, each resistor occupies its own whole-volume caloric body.
    pub graph_resistors: &'a [ResistorCell],
}

/// One jointly accepted electrode circuit, original field and finite heat state.
pub struct SpatialTerminalCircuitThermalCandidate<'t, M, D, S = MaterialNetworkSample> {
    pub(in crate::circuit_enthalpy::spatial_induction) thermal: DepositionSimulation<'t>,
    pub(in crate::circuit_enthalpy::spatial_induction) circuit: M,
    pub(in crate::circuit_enthalpy::spatial_induction) report: SpatialInductionThermalReport<D, S>,
}

impl<'t, M, D, S> SpatialTerminalCircuitThermalCandidate<'t, M, D, S> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn circuit(&self) -> &M {
        &self.circuit
    }
    pub fn report(&self) -> &SpatialInductionThermalReport<D, S> {
        &self.report
    }
    /// Physical field and authored resistor heat already paid to finite cells.
    pub fn deposited_heat_j(&self) -> f64 {
        self.report.coupled.deposited_joule_j
    }
    /// Publish the accepted circuit/field and thermal histories together.
    /// Physical heat is already deposited. Numerical loss and the canceled
    /// interface-reference work must never be deposited by the caller.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        M,
        SpatialInductionThermalReport<D, S>,
    ) {
        (self.thermal, self.circuit, self.report)
    }
}

pub(in crate::circuit_enthalpy::spatial_induction) type GraphRoutes = BTreeMap<usize, CellId>;

impl SpatialTerminalCircuitThermalNetwork<'_> {
    pub(in crate::circuit_enthalpy::spatial_induction) fn routes(
        &self,
        thermal: &DepositionSimulation<'_>,
        accepted_time_s: f64,
        graph: &CircuitGraph,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<GraphRoutes, Error> {
        check(cancel)?;
        if accepted_time_s.to_bits() != thermal.time_s().to_bits() {
            return Err(invalid("electrode circuit and thermal clocks must agree"));
        }
        let mut routes = BTreeMap::new();
        for route in self.graph_resistors {
            check(cancel)?;
            if !thermal.cells().contains_key(&route.cell)
                || !matches!(
                    graph.branches.get(route.branch),
                    Some((_, _, Branch::Resistor { .. }))
                )
                || routes.insert(route.branch, route.cell).is_some()
            {
                return Err(invalid(
                    "authored graph resistor heat needs unique existing branch/cell correspondence",
                ));
            }
        }
        let mut expected = 0;
        for (_, _, branch) in &graph.branches {
            check(cancel)?;
            if matches!(branch, Branch::Resistor { .. }) {
                expected += 1;
            }
        }
        if routes.len() != expected {
            return Err(invalid(
                "every authored graph resistor needs a finite thermal receiver",
            ));
        }
        Ok(routes)
    }

    pub(in crate::circuit_enthalpy::spatial_induction) fn admit_graph_material_bindings(
        &self,
        retained: &[GraphResistorMaterial],
        bindings: Option<&[CellResistor]>,
        routes: &GraphRoutes,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<(), Error> {
        check(cancel)?;
        let Some(bindings) = bindings else {
            if !retained.is_empty() {
                return Err(invalid(
                    "retained graph resistor sources require their complete thermal bindings",
                ));
            }
            return Ok(());
        };
        if bindings.len() != routes.len() {
            return Err(invalid("complete graph material/thermal bindings required"));
        }
        if bindings.is_empty() {
            return Ok(());
        }
        let mut magnetic_cells = BTreeSet::new();
        for region in self.magnetic.magnetic_regions {
            check(cancel)?;
            if let SpatialMagneticRegion::Cell { cell, .. } = *region {
                magnetic_cells.insert(cell);
            }
        }
        let mut branches = BTreeSet::new();
        let mut cells = BTreeSet::new();
        for binding in bindings {
            check(cancel)?;
            if routes.get(&binding.destination.branch) != Some(&binding.destination.cell)
                || !branches.insert(binding.destination.branch)
                || !cells.insert(binding.destination.cell)
                || magnetic_cells.contains(&binding.destination.cell)
            {
                return Err(invalid(
                    "graph material bindings need exact routes and separate whole-cell bodies",
                ));
            }
        }
        Ok(())
    }
}

pub(in crate::circuit_enthalpy::spatial_induction) fn deposit_graph_heat(
    losses: &[WindingGraphResistorHeat],
    routes: &GraphRoutes,
    powers: &mut BTreeMap<CellId, f64>,
    dt: f64,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<(), Error> {
    if losses.len() != routes.len() || !dt.is_finite() || dt <= 0.0 {
        return Err(invalid(
            "complete authored graph heat receipts and positive duration required",
        ));
    }
    let mut seen = BTreeSet::new();
    for loss in losses {
        checkpoint()?;
        let receiver = *routes
            .get(&loss.branch)
            .ok_or_else(|| invalid("authored graph heat lacks its original finite receiver"))?;
        if !seen.insert(loss.branch) {
            return Err(invalid("repeated authored graph heat receipt"));
        }
        let heat = finite(loss.heat_j)?;
        let power = finite(heat / dt)?;
        if heat < 0.0 || heat != 0.0 && power == 0.0 {
            return Err(invalid("negative or unrepresentable authored graph heat"));
        }
        let previous = powers.get(&receiver).copied().unwrap_or(0.0);
        let combined = finite(previous + power)?;
        if power != 0.0 && combined == previous {
            return Err(invalid("authored graph heat disappeared in receiver sum"));
        }
        powers.insert(receiver, combined);
    }
    Ok(())
}

pub(in crate::circuit_enthalpy::spatial_induction) fn close_heat_balance<
    D,
    S: AsRef<HeatNetworkSample>,
    E: From<HeatNetworkError>,
>(
    report: &mut SpatialInductionThermalReport<D, S>,
    dt: f64,
    tolerance: f64,
) -> Result<(), TerminalError<E>> {
    super::super::heat::close_heat_balance::<D, S, E>(report, dt, tolerance).map_err(|error| {
        error.map_source(|source| match source {
            super::super::heat::SpatialInductionHeatError::Native(error) => {
                SpatialTerminalCircuitThermalError::Native(error)
            }
            super::super::heat::SpatialInductionHeatError::Heat(error) => {
                SpatialTerminalCircuitThermalError::Heat(error)
            }
        })
    })
}
