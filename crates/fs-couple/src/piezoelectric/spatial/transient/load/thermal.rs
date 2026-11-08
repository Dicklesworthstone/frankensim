//! Spatial piezoelectric loads on the existing source-pinned thermal network.
//!
//! Each nonlinear trial restarts the original spatial midpoint solver from its
//! accepted history. The original material network resolves resistor and thermal
//! conductivities at the proposed enthalpies; the original thermal driver owns
//! deposition, conduction, Newton iteration and joint energy acceptance.
//! Explicit ceramic routes also deposit intrinsic material losses into their
//! original finite caloric inventory, with unchanged source-backed coefficients.
//! Incident and outgoing guide waves share that same solve. Their signed work
//! crosses the thermal assembly boundary; it is never deposited as local heat.

mod ceramic;
pub use ceramic::ThermalCeramicCell;
/// Full tensor temperature feedback and reversible heat on a finite ceramic body.
pub mod tensor;

use super::{ElectrodeDrive, LoadedPiezoelectricStep, ResistiveElectrodeLoad};
use crate::circuit_enthalpy::{
    CircuitEnthalpyError, CircuitEnthalpyReport, CircuitEnthalpySolve, DrivenCoefficients,
    DriverAccount, ResistorCell, ThermalDriver, admit_driver_unknowns, cell_driver_candidate,
    material::{MaterialNetworkError, MaterialNetworkSample, PreparedMaterialNetwork},
};
use crate::material_circuit::MaterialResistor;
use crate::piezoelectric::spatial::transient::{
    TetPiezoelectricTransient, TetPiezoelectricTransientCandidate, TetPiezoelectricTransientStep,
    acoustic::{TransientAcousticPort, TransientAcousticPortResponse},
};
use crate::piezoelectric::spatial::{TetPiezoelectricControls, TetPiezoelectricError};
use fs_exec::Cx;
use fs_math::dd::Dd;
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
    /// Midpoint support reactions, including inertia and intrinsic memory.
    /// Free mechanical coordinates have zero entries. Their work is retained
    /// separately in `native.support_work_j` and included in the joint supply.
    pub midpoint_support_force_n: Vec<[f64; 3]>,
    /// Incoming/outgoing guide waves and signed work, in input port order.
    /// This energy crosses the assembly boundary, separately from deposited heat.
    pub acoustic_ports: Vec<TransientAcousticPortResponse>,
    /// Original material-network heat routes, in the same order as `loads.loads`.
    pub destinations: Vec<ResistorCell>,
    /// Complete original-tetrahedron to ceramic caloric-cell correspondence.
    /// Heat amounts remain in the native body's per-cell relaxation reports.
    pub ceramic_destinations: Vec<ThermalCeramicCell>,
}

/// Original coupled thermal account with complete accepted material receipts.
/// The common report's `deposited_joule_j` includes both conductor and intrinsic
/// ceramic losses here. Its native/load subreports retain their separate owners.
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
    candidate: Arc<TetPiezoelectricTransientCandidate<'model, 'mesh>>,
    loads: LoadedPiezoelectricStep,
    destinations: Vec<ResistorCell>,
    ceramic_destinations: Arc<[ThermalCeramicCell]>,
}
impl<'model, 'mesh> ThermalDriver for LoadDriver<'model, 'mesh> {
    type Model = TetPiezoelectricTransient<'model, 'mesh>;
    type Step = SpatialPiezoelectricThermalStep;
    fn account(&self) -> DriverAccount<'_> {
        let native = self.candidate.step();
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: self.loads.joule_heat_j + native.material_heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: native.retained_energy_change_j().to_f64(),
            supplied_j: (native.retained_mechanical_work_j()
                + self.loads.retained_electrical_source_work_j()
                - Dd::from_f64(native.acoustic_outward_work_j)).to_f64(),
        }
    }
    fn retained_energy_balance_j(&self) -> Dd {
        let native = self.candidate.step();
        native.retained_energy_change_j() - native.retained_mechanical_work_j()
            - self.loads.retained_electrical_source_work_j()
            + Dd::from_f64(native.acoustic_outward_work_j)
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        let mut candidate = Arc::unwrap_or_clone(self.candidate);
        let midpoint_support_force_n = std::mem::take(&mut candidate.midpoint_support_force_n);
        let acoustic_ports = std::mem::take(&mut candidate.acoustic_ports);
        let (body, native) = candidate.into_parts();
        (
            body,
            SpatialPiezoelectricThermalStep {
                native,
                loads: self.loads,
                midpoint_support_force_n,
                acoustic_ports,
                destinations: self.destinations,
                ceramic_destinations: self.ceramic_destinations.to_vec(),
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
        self.step_with_material_thermal_loads_candidate(
            thermal,
            target_s,
            material,
            bindings,
            &[],
            "T",
            drives,
            midpoint_nodal_forces_n,
            endpoint_nodal_free_charge_c,
            field_controls,
            thermal_controls,
            guess_j,
            cx,
        )
    }

    /// Prescribed-support excitation with physical resistor heating in the
    /// original finite thermal network. Each original constrained DOF must have
    /// one endpoint velocity, including stationary supports. Its displacement
    /// follows the native midpoint kinematics; the full consistent mass and
    /// reciprocal electroelastic response determine the actual support work.
    ///
    /// The native support reaction and work survive every thermal trial.
    /// Resistor temperature feedback, heat conduction and joint acceptance use
    /// that mechanical supply without treating it as an extra heat source.
    /// Both accepted histories remain unchanged until the candidate is consumed.
    /// The same fixed-coefficient ceramic scope as
    /// [`Self::step_with_thermal_loads_candidate`] applies.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_thermal_support_motion_candidate<'thermal>(
        &self,
        thermal: &DepositionSimulation<'thermal>,
        target_s: f64,
        endpoint_support_velocities_m_s: &[(usize, f64)],
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
        self.advance_material_thermal_candidate(
            thermal,
            target_s,
            Some(endpoint_support_velocities_m_s),
            &[],
            material,
            bindings,
            &[],
            "T",
            drives,
            midpoint_nodal_forces_n,
            endpoint_nodal_free_charge_c,
            field_controls,
            thermal_controls,
            guess_j,
            cx,
        )
    }

    /// Evolve intrinsic ceramic heating, physical resistor losses and the
    /// original finite thermal network in one jointly accepted candidate.
    /// Every original tetrahedron needs one explicit ceramic destination when
    /// intrinsic memories are attached. Several tetrahedra may share one cell
    /// only when their complete original volume and mass match that cell.
    ///
    /// Ceramic and caloric owners must retain the same card, complete nonthermal
    /// context, density claim and actual borrowed caloric law. The full caloric
    /// chart remains solid with constant reference density. Every accepted and
    /// proposed temperature must lie inside all selected piezoelectric and
    /// relaxation source claims; changed temperature refuses any coefficient
    /// that depends on the owned temperature axis. This admits finite heating
    /// of a fixed-coefficient ceramic, without inventing parameter work or
    /// changing the native histories. Thermal strain, changing tensors, phase
    /// transition and depoling require separately accounted models.
    ///
    /// Ceramic heat enters the native deposition callback exactly once and is
    /// internal transfer, not an additional external heat source. All existing
    /// conductivity, resistor feedback, Newton budgets and physical audits stay
    /// authoritative. A refusal leaves both supplied histories unchanged.
    /// With no resistor feedback, one immutable spatial candidate is reused
    /// across thermal trials; every trial still validates its actual sources.
    /// Endpoint-dependent resistor and transport coefficients retain the
    /// original generally first-order thermal time discretization.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_material_thermal_loads_candidate<'thermal>(
        &self,
        thermal: &DepositionSimulation<'thermal>,
        target_s: f64,
        material: &PreparedMaterialNetwork<'_, '_>,
        bindings: &[ThermalElectrodeLoad],
        ceramic_bindings: &[ThermalCeramicCell],
        temperature_axis: &str,
        drives: &[ElectrodeDrive],
        midpoint_nodal_forces_n: &[[f64; 3]],
        endpoint_nodal_free_charge_c: &[f64],
        field_controls: TetPiezoelectricControls,
        thermal_controls: CircuitEnthalpySolve,
        guess_j: &[f64],
        cx: &Cx<'_>,
    ) -> Result<SpatialPiezoelectricThermalCandidate<'thermal, 'model, 'mesh>> {
        self.advance_material_thermal_candidate(
            thermal,
            target_s,
            None,
            &[],
            material,
            bindings,
            ceramic_bindings,
            temperature_axis,
            drives,
            midpoint_nodal_forces_n,
            endpoint_nodal_free_charge_c,
            field_controls,
            thermal_controls,
            guess_j,
            cx,
        )
    }

    /// Prescribed-support excitation with both intrinsic ceramic losses and
    /// physical shunt losses deposited into their original finite caloric cells.
    ///
    /// Support input has the same complete endpoint-velocity convention as
    /// [`Self::step_with_support_motion_candidate`]. Each nonlinear proposal
    /// restarts the accepted displacement, velocity, charge and material memory.
    /// The original native solve retains support inertia and memory forces;
    /// actual support work enters the combined electrical/mechanical/thermal
    /// audit exactly once. Ceramic loss is internal transfer to its named cell.
    ///
    /// All source, inventory and unchanged-coefficient restrictions of
    /// [`Self::step_with_material_thermal_loads_candidate`] still apply. A
    /// cancelled, out-of-domain or unconverged proposal publishes neither body
    /// nor heat state. When no resistors are present, thermal trials
    /// reuse one immutable native moving-support candidate.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_material_thermal_support_motion_candidate<'thermal>(
        &self,
        thermal: &DepositionSimulation<'thermal>,
        target_s: f64,
        endpoint_support_velocities_m_s: &[(usize, f64)],
        material: &PreparedMaterialNetwork<'_, '_>,
        bindings: &[ThermalElectrodeLoad],
        ceramic_bindings: &[ThermalCeramicCell],
        temperature_axis: &str,
        drives: &[ElectrodeDrive],
        midpoint_nodal_forces_n: &[[f64; 3]],
        endpoint_nodal_free_charge_c: &[f64],
        field_controls: TetPiezoelectricControls,
        thermal_controls: CircuitEnthalpySolve,
        guess_j: &[f64],
        cx: &Cx<'_>,
    ) -> Result<SpatialPiezoelectricThermalCandidate<'thermal, 'model, 'mesh>> {
        self.advance_material_thermal_candidate(
            thermal,
            target_s,
            Some(endpoint_support_velocities_m_s),
            &[],
            material,
            bindings,
            ceramic_bindings,
            temperature_axis,
            drives,
            midpoint_nodal_forces_n,
            endpoint_nodal_free_charge_c,
            field_controls,
            thermal_controls,
            guess_j,
            cx,
        )
    }

    /// Receive and emit sound while evolving source-bound load and ceramic heat.
    /// The original guide traction, electrode currents, consistent inertia and
    /// intrinsic material memories participate in every joint thermal trial.
    /// Signed outward acoustic work is an external energy transfer, counted once;
    /// only resistor and intrinsic material losses enter their finite heat cells.
    ///
    /// Use `None` for stationary supports or supply every original support DOF's
    /// endpoint velocity. Ceramic routes, source domains and thermal feedback
    /// retain the restrictions of [`Self::step_with_material_thermal_loads_candidate`].
    /// The guide's fluid owner is fixed during this interval; fluid heating and
    /// temperature-dependent ceramic tensors require separate coupled models.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_acoustic_thermal_loads_candidate<'thermal>(
        &self,
        thermal: &DepositionSimulation<'thermal>,
        target_s: f64,
        acoustic_ports: &[TransientAcousticPort<'_>],
        endpoint_support_velocities_m_s: Option<&[(usize, f64)]>,
        material: &PreparedMaterialNetwork<'_, '_>,
        bindings: &[ThermalElectrodeLoad],
        ceramic_bindings: &[ThermalCeramicCell],
        temperature_axis: &str,
        drives: &[ElectrodeDrive],
        midpoint_nodal_forces_n: &[[f64; 3]],
        endpoint_nodal_free_charge_c: &[f64],
        field_controls: TetPiezoelectricControls,
        thermal_controls: CircuitEnthalpySolve,
        guess_j: &[f64],
        cx: &Cx<'_>,
    ) -> Result<SpatialPiezoelectricThermalCandidate<'thermal, 'model, 'mesh>> {
        self.advance_material_thermal_candidate(
            thermal, target_s, endpoint_support_velocities_m_s, acoustic_ports,
            material, bindings, ceramic_bindings, temperature_axis, drives,
            midpoint_nodal_forces_n, endpoint_nodal_free_charge_c, field_controls,
            thermal_controls, guess_j, cx,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_material_thermal_candidate<'thermal>(
        &self,
        thermal: &DepositionSimulation<'thermal>,
        target_s: f64,
        endpoint_support_velocities_m_s: Option<&[(usize, f64)]>,
        acoustic_ports: &[TransientAcousticPort<'_>],
        material: &PreparedMaterialNetwork<'_, '_>,
        bindings: &[ThermalElectrodeLoad],
        ceramic_bindings: &[ThermalCeramicCell],
        temperature_axis: &str,
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
        let ceramic = ceramic::CeramicPlan::prepare(
            self,
            thermal,
            material,
            ceramic_bindings,
            temperature_axis,
            cx,
        )?;
        let ceramic_destinations: Arc<[ThermalCeramicCell]> = ceramic.destinations().into();
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
        let original = resolve(thermal)?;
        ceramic.admit(self, thermal, &original.source, cx)?;
        let original_loads = bind_loads(original.resistors, bindings)?;
        let spatial_candidate = |loads: &[ResistiveElectrodeLoad]| {
            let (mut candidate, loads) = self
                .advance_loaded_with_acoustic_candidate(
                    dt,
                    endpoint_support_velocities_m_s,
                    acoustic_ports,
                    drives,
                    loads,
                    midpoint_nodal_forces_n,
                    endpoint_nodal_free_charge_c,
                    field_controls,
                    cx,
                )
                .map_err(|error| match error {
                    TetPiezoelectricError::Cancelled => CircuitEnthalpyError::Cancelled,
                    error => CircuitEnthalpyError::Source(LoadThermalError::Piezoelectric(error)),
                })?;
            // Preserve the caller's absolute endpoint rather than a separately
            // rounded old_time + dt. Field integration still uses the exact dt.
            candidate.body.state.time_s = target_s;
            Ok::<_, CircuitEnthalpyError<LoadThermalError>>((Arc::new(candidate), loads))
        };
        // Held ceramic coefficients have no thermal feedback. When no resistor
        // is present, source/heat-network Newton trials share one native field
        // candidate; only final publication may clone its accepted state once.
        let fixed_spatial = if original_loads.is_empty() {
            Some(spatial_candidate(&[])?)
        } else {
            None
        };
        let (thermal, body, report) = cell_driver_candidate(
            thermal,
            target_s,
            thermal_controls,
            guess_j,
            |trial, _remaining| {
                let coefficients = resolve(trial)?;
                ceramic.admit(self, trial, &coefficients.source, cx)?;
                let loads = bind_loads(coefficients.resistors, bindings)?;
                let (candidate, loads) = if let Some((candidate, loads)) = &fixed_spatial {
                    (Arc::clone(candidate), loads.clone())
                } else {
                    spatial_candidate(&loads)?
                };
                Ok(DrivenCoefficients {
                    circuit: LoadDriver {
                        candidate,
                        loads,
                        destinations: destinations.clone(),
                        ceramic_destinations: Arc::clone(&ceramic_destinations),
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
                ceramic.deposit(driver.candidate.body(), dt, &mut powers, checkpoint)?;
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
