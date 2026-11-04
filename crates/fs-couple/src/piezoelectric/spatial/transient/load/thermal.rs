//! Spatial piezoelectric loads on the existing source-pinned thermal network.
//!
//! Each nonlinear trial restarts the original spatial midpoint solver from its
//! accepted history. The original material network resolves resistor and thermal
//! conductivities at the proposed enthalpies; the original thermal driver owns
//! deposition, conduction, Newton iteration and joint energy acceptance.

use super::{ElectrodeDrive, LoadedPiezoelectricStep, ResistiveElectrodeLoad};
use crate::circuit_enthalpy::{
    CircuitEnthalpyError, CircuitEnthalpyReport, CircuitEnthalpySolve, DrivenCoefficients,
    DriverAccount, ResistorCell, ThermalDriver, admit_driver_unknowns, cell_driver_candidate,
    material::{MaterialNetworkError, MaterialNetworkSample, PreparedMaterialNetwork},
};
use crate::material_circuit::MaterialResistor;
use crate::piezoelectric::spatial::transient::{
    TetPiezoelectricTransient, TetPiezoelectricTransientCandidate, TetPiezoelectricTransientStep,
};
use crate::piezoelectric::spatial::{TetPiezoelectricControls, TetPiezoelectricError};
use fs_exec::Cx;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};
use std::collections::BTreeMap;
use std::sync::Arc;

/// Connect an original material-network resistor to two spatial electrodes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ThermalElectrodeLoad {
    /// Stable resistor branch identifier retained by the prepared network.
    pub resistor_branch: usize,
    /// Positive load current flows from the first electrode to the second.
    pub electrodes: [usize; 2],
}

/// Preserve the original material and spatial-field refusals through Newton.
#[derive(Debug)]
pub enum LoadThermalError {
    /// Caloric/source identity, transport-domain or geometry refusal.
    Material(MaterialNetworkError),
    /// Original spatial electromechanical solve or physical work refusal.
    Piezoelectric(TetPiezoelectricError),
}
impl core::fmt::Display for LoadThermalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spatial piezoelectric load heating: {self:?}")
    }
}
impl std::error::Error for LoadThermalError {}
type Result<T> = core::result::Result<T, CircuitEnthalpyError<LoadThermalError>>;

/// Native motion, electrode work and the exact physical heat receivers.
#[derive(Clone, Debug)]
pub struct SpatialPiezoelectricThermalStep {
    /// Original spatial midpoint step, including momentum and charge residuals.
    pub native: TetPiezoelectricTransientStep,
    /// Actual external electrical work and source-bound resistor losses.
    pub loads: LoadedPiezoelectricStep,
    /// Original material-network heat routes, in the same order as `loads.loads`.
    pub destinations: Vec<ResistorCell>,
}

/// Original coupled thermal account with complete accepted material receipts.
pub type SpatialPiezoelectricThermalReport =
    CircuitEnthalpyReport<MaterialNetworkSample, SpatialPiezoelectricThermalStep>;

/// Both histories after one complete admitted interval; no input is mutated.
pub struct SpatialPiezoelectricThermalCandidate<'thermal, 'model, 'mesh> {
    thermal: DepositionSimulation<'thermal>,
    body: TetPiezoelectricTransient<'model, 'mesh>,
    report: SpatialPiezoelectricThermalReport,
}
impl<'thermal, 'model, 'mesh> SpatialPiezoelectricThermalCandidate<'thermal, 'model, 'mesh> {
    /// Candidate finite caloric state, including passive conductive neighbours.
    pub const fn thermal(&self) -> &DepositionSimulation<'thermal> {
        &self.thermal
    }
    /// Candidate displacement, momentum, electrode charge and dielectric field.
    pub const fn body(&self) -> &TetPiezoelectricTransient<'model, 'mesh> {
        &self.body
    }
    /// Physical work/heat closure and the original accepted material evaluation.
    pub const fn report(&self) -> &SpatialPiezoelectricThermalReport {
        &self.report
    }
    /// Publish the complete mechanical/electrical/thermal candidate together.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'thermal>,
        TetPiezoelectricTransient<'model, 'mesh>,
        SpatialPiezoelectricThermalReport,
    ) {
        (self.thermal, self.body, self.report)
    }
}

struct LoadDriver<'model, 'mesh> {
    candidate: TetPiezoelectricTransientCandidate<'model, 'mesh>,
    loads: LoadedPiezoelectricStep,
    destinations: Vec<ResistorCell>,
}
impl<'model, 'mesh> ThermalDriver for LoadDriver<'model, 'mesh> {
    type Model = TetPiezoelectricTransient<'model, 'mesh>;
    type Step = SpatialPiezoelectricThermalStep;
    fn account(&self) -> DriverAccount<'_> {
        let native = self.candidate.step();
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: self.loads.joule_heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: native.kinetic_energy_change_j + native.stored_energy_change_j,
            supplied_j: native.mechanical_work_j + self.loads.electrical_source_work_j,
        }
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        let (body, native) = self.candidate.into_parts();
        (
            body,
            SpatialPiezoelectricThermalStep {
                native,
                loads: self.loads,
                destinations: self.destinations,
            },
        )
    }
}

fn bind_loads(
    resistors: Vec<MaterialResistor>,
    bindings: &[ThermalElectrodeLoad],
) -> Result<Vec<ResistiveElectrodeLoad>> {
    let mut sources = BTreeMap::new();
    for resistor in resistors {
        if sources
            .insert(resistor.branch, resistor.conductor)
            .is_some()
        {
            return Err(CircuitEnthalpyError::Invalid(
                "duplicate thermal electrode resistor",
            ));
        }
    }
    let loads =
        bindings
            .iter()
            .map(|binding| {
                let conductor = sources.remove(&binding.resistor_branch).ok_or(
                    CircuitEnthalpyError::Invalid("missing thermal electrode resistor"),
                )?;
                Ok(ResistiveElectrodeLoad {
                    electrodes: binding.electrodes,
                    conductor: Arc::new(conductor),
                })
            })
            .collect::<Result<Vec<_>>>()?;
    if !sources.is_empty() {
        return Err(CircuitEnthalpyError::Invalid(
            "unconnected material-network resistor",
        ));
    }
    Ok(loads)
}

impl<'model, 'mesh> TetPiezoelectricTransient<'model, 'mesh> {
    /// Evolve the original spatial body and finite material resistor cells at
    /// one joint thermal endpoint. Each prepared resistor must be connected
    /// exactly once; its original network binding alone selects the heat cell.
    /// All native conductive links, caloric laws and authored boundaries remain
    /// in the same enthalpy residual. No ceramic or artificial damping heat is
    /// added. Source selections stay pinned while R(T) and conduction refresh.
    ///
    /// Spatial and thermal clocks must agree. Every trial restarts the accepted
    /// displacement, velocity and charge, and actual midpoint load heat is
    /// deposited exactly once. The returned candidate passes the common thermal
    /// and combined storage/work audits before either history can be published.
    ///
    /// The original fixed-geometry, constant-tensor, small-strain ceramic and
    /// scalar stationary Ohmic load scopes apply. Only the external resistors
    /// and thermal network change temperature here. Thermal strain, dielectric
    /// loss, depoling and changing ceramic coefficients are not inferred.
    /// Midpoint motion with endpoint thermal coefficients is first order in
    /// general; energy closure does not certify temporal or spatial accuracy.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_thermal_loads_candidate<'thermal>(
        &self,
        thermal: &DepositionSimulation<'thermal>,
        target_s: f64,
        material: &PreparedMaterialNetwork<'_, '_>,
        bindings: &[ThermalElectrodeLoad],
        drives: &[ElectrodeDrive],
        midpoint_nodal_forces_n: &[[f64; 3]],
        endpoint_nodal_free_charge_c: &[f64],
        field_controls: TetPiezoelectricControls,
        thermal_controls: CircuitEnthalpySolve,
        guess_j: &[f64],
        cx: &Cx<'_>,
    ) -> Result<SpatialPiezoelectricThermalCandidate<'thermal, 'model, 'mesh>> {
        let cancelled = || cx.checkpoint().is_err();
        if cancelled() {
            return Err(CircuitEnthalpyError::Cancelled);
        }
        admit_driver_unknowns(thermal, thermal_controls, guess_j)?;
        let dt = target_s - thermal.time_s();
        if self.state().time_s() != thermal.time_s() || !dt.is_finite() || dt <= 0.0 {
            return Err(CircuitEnthalpyError::Invalid(
                "spatial piezoelectric and thermal clocks must share one positive interval",
            ));
        }
        let mut routes = material
            .resistor_destinations()
            .into_iter()
            .map(|route| (route.branch, route))
            .collect::<BTreeMap<_, _>>();
        if bindings.len() != routes.len() {
            return Err(CircuitEnthalpyError::Invalid(
                "every material-network resistor requires one thermal electrode binding",
            ));
        }
        let mut destinations = Vec::with_capacity(bindings.len());
        for binding in bindings {
            if cancelled() {
                return Err(CircuitEnthalpyError::Cancelled);
            }
            destinations.push(routes.remove(&binding.resistor_branch).ok_or(
                CircuitEnthalpyError::Invalid(
                    "duplicate or foreign thermal electrode load binding",
                ),
            )?);
        }
        if !routes.is_empty() {
            return Err(CircuitEnthalpyError::Invalid(
                "unconnected material-network resistor",
            ));
        }
        let resolve = |state: &DepositionSimulation<'_>| {
            material
                .resolve(state, cancelled)
                .map_err(|error| match error {
                    MaterialNetworkError::Thermal(DepositionError::Cancelled) => {
                        CircuitEnthalpyError::Cancelled
                    }
                    error => CircuitEnthalpyError::Source(LoadThermalError::Material(error)),
                })
        };
        // Admit the existing source/geometry at the accepted state as well as
        // the proposed endpoint; a later in-domain trial cannot repair it.
        bind_loads(resolve(thermal)?.resistors, bindings)?;
        let (thermal, body, report) = cell_driver_candidate(
            thermal,
            target_s,
            thermal_controls,
            guess_j,
            |trial, _remaining| {
                let coefficients = resolve(trial)?;
                let loads = bind_loads(coefficients.resistors, bindings)?;
                let (mut candidate, loads) = self
                    .step_with_loads_candidate(
                        dt,
                        drives,
                        &loads,
                        midpoint_nodal_forces_n,
                        endpoint_nodal_free_charge_c,
                        field_controls,
                        cx,
                    )
                    .map_err(|error| match error {
                        TetPiezoelectricError::Cancelled => CircuitEnthalpyError::Cancelled,
                        error => {
                            CircuitEnthalpyError::Source(LoadThermalError::Piezoelectric(error))
                        }
                    })?;
                // The native duration solve retains its exact dt and physical
                // result. Stamp the shared absolute endpoint instead of the
                // separately rounded old_time + dt, which can differ by one ULP
                // and otherwise make the next coupled interval refuse.
                candidate.body.state.time_s = target_s;
                Ok(DrivenCoefficients {
                    circuit: LoadDriver {
                        candidate,
                        loads,
                        destinations: destinations.clone(),
                    },
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    source: coefficients.source,
                    circuit_trials: 1,
                })
            },
            |driver, dt, checkpoint| {
                let mut powers = BTreeMap::new();
                for (heat, route) in driver.loads.loads.iter().zip(&driver.destinations) {
                    checkpoint()?;
                    let power = heat.joule_heat_j / dt;
                    let sum = powers.entry(route.cell).or_insert(0.0);
                    *sum += power;
                    if !power.is_finite()
                        || power < 0.0
                        || !sum.is_finite()
                        || (heat.joule_heat_j != 0.0 && power == 0.0)
                    {
                        return Err(DepositionError::Arithmetic.into());
                    }
                }
                Ok(powers)
            },
            cancelled,
        )?;
        Ok(SpatialPiezoelectricThermalCandidate {
            thermal,
            body,
            report,
        })
    }
}
