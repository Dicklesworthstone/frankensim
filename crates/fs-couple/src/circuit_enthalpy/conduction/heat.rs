//! Spatial current and complete sourced heat transport on one physical endpoint.
//! The prepared heat assembly supplies the SAME cells/cards used by the field
//! material refresh. Its original interface and radiosity owners supply heat;
//! the original spatial electrical solve alone supplies Joule dissipation.

mod assembly;
pub use assembly::{ConductionHeatError, step_response_field_circuit_candidate};

use super::circuit::response::{
    ConductionResponseCircuit, ConductionResponseCircuitEnthalpyCandidate,
    thermal::ConductionResponseFieldThermalCandidate,
};
use super::response::ConductionResponseEnthalpyCandidate;
use super::*;
use crate::circuit_enthalpy::contacts::heat::{
    HeatNetworkError, HeatNetworkSample, PreparedHeatNetwork, closed_balance,
};
use crate::circuit_enthalpy::{ResistorCell, field_storage::FieldThermalBinding};
use fs_em::conduction::response::{ConductionResponseControls, PlanarConductionResponse};
use fs_exec::Cx;

fn heat_error(error: HeatNetworkError) -> Error {
    match error {
        HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(ConductionThermalError::Heat(other)),
    }
}

impl PreparedHeatNetwork<'_, '_> {
    /// Advance nonlinear spatial current, native electrical contacts and sourced
    /// R(T)/C(T)/L(T) with this complete heat assembly at ONE finite endpoint.
    ///
    /// Bulk and surface Joule heat, reversible field-entropy transfer, original
    /// thermal contacts, source-integrated layers and closed-enclosure radiosity
    /// all enter the existing circuit/enthalpy residual once. Each constitutive
    /// trial starts from the same accepted charge/flux. The returned source
    /// sample retains every original thermal/interface/radiation calculation.
    ///
    /// Graph resistors come from this assembly's original material bindings;
    /// the explicit spatial contact bindings retain their native two-sided heat
    /// owner. Series graph electrical contacts refuse because this spatial
    /// circuit has no owner for their work. Every field branch needs an original
    /// whole-volume caloric binding, with the existing affine field law,
    /// constant-density solid chart and positive capacity/coordinate envelope.
    ///
    /// Physical field U = psi + T s and caloric storage close the original
    /// external supply. Closed radiation's signed numerical net remains visible
    /// in the heat sample and is excluded from external_heat_j; it is never
    /// redistributed into cells. E/J/B/H-dependent thermal boundaries refuse
    /// without a resolved interface field. Geometry and enclosure view factors
    /// remain fixed. Coupling is first order and requires timestep refinement;
    /// no phase motion, hysteresis, dielectric relaxation or fitted loss is added.
    #[allow(clippy::too_many_arguments)]
    pub fn step_conduction_response_field_circuit_candidate<'t>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &ConductionResponseCircuit,
        ports: &[f64],
        contacts: &[contacts::ConductionContactThermal<'_>],
        fields: &[FieldThermalBinding],
        regions: &[ConductionRegion],
        volume_tolerance_m3: f64,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        response: ConductionResponseControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
    ) -> Result<ConductionResponseFieldThermalCandidate<'t, HeatNetworkSample>, Error> {
        step_response_field_circuit_candidate(
            self,
            cx,
            thermal,
            target_s,
            circuit,
            ports,
            contacts,
            fields,
            regions,
            volume_tolerance_m3,
            assembly,
            field_controls,
            response,
            controls,
            initial_energy_increments_j,
        )
        .map_err(|error| {
            error.map_source(|source| match source {
                ConductionHeatError::Conduction(error) => error,
                ConductionHeatError::Heat(error) => ConductionThermalError::Heat(error),
            })
        })
    }

    /// Advance the original nonlinear spatial circuit with this complete
    /// material heat assembly. Electrical charge/flux, bulk E/J response,
    /// original electrical-contact loss, thermal-contact conduction and
    /// enclosure radiation see the SAME endpoint enthalpies and source claims.
    /// Every field trial starts from the accepted circuit history, and every
    /// heat trial starts from the accepted finite caloric history.
    ///
    /// Original thermal surfaces and layer laws are evaluated by this prepared
    /// network; their samples remain in the returned report. Bulk and contact
    /// Joule heat enter once through the native circuit field. Other fixed graph
    /// resistors require explicit receiving cells in `graph_resistors`; this
    /// route does not silently omit material graph-resistor/contact bindings.
    /// All boundaries must be independent of E/J in response and applicability.
    /// Closed-enclosure radiation is internal exchange and its numerical net is
    /// checked against the original combined-energy allowance.
    #[allow(clippy::too_many_arguments)]
    pub fn step_conduction_response_contact_circuit_candidate<'t>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &ConductionResponseCircuit,
        ports: &[f64],
        graph_resistors: &[ResistorCell],
        contacts: &[contacts::ConductionContactThermal<'_>],
        regions: &[ConductionRegion],
        volume_tolerance_m3: f64,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        response: ConductionResponseControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
    ) -> Result<ConductionResponseCircuitEnthalpyCandidate<'t, HeatNetworkSample>, Error> {
        let cancel = RefCell::new(|| cx.checkpoint().is_err());
        check(&cancel)?;
        let bulk = self.thermal_material_network().map_err(heat_error)?;
        let layout = ConductionThermalNetwork {
            cells: bulk.cells,
            links: bulk.links,
            regions,
            volume_tolerance_m3,
        };
        if regions.len() != circuit.sources().len() || regions.len() > assembly.max_cells {
            return Err(invalid(
                "one bounded heat route per original circuit triangle required",
            ));
        }
        let mut field_axes = Vec::new();
        for region in regions {
            check(&cancel)?;
            let source = circuit
                .sources()
                .get(region.triangle)
                .ok_or_else(|| invalid("thermal route names an absent circuit conductor"))?;
            for (axis, _) in source.axes() {
                field_axes.push((region.cell, axis.as_str()));
            }
        }
        self.admit_field_independent_boundaries(cx, &field_axes)
            .map_err(heat_error)?;
        let next = layout.step_response_circuit_resolved_candidate(
            thermal,
            target_s,
            circuit,
            ports,
            graph_resistors,
            contacts,
            assembly,
            field_controls,
            response,
            controls,
            initial_energy_increments_j,
            |trial| {
                let coefficients = self.resolve(cx, trial).map_err(heat_error)?;
                if !coefficients.contacts.is_empty() || !coefficients.bulk.resistors.is_empty() {
                    return Err(invalid(
                        "spatial response circuit cannot omit other source-bound electrical work",
                    ));
                }
                Ok(coefficients.bulk)
            },
            |sample| &sample.material,
            &cancel,
        )?;
        let report = next.report();
        closed_balance(
            target_s - thermal.time_s(),
            report.total_energy_residual_j,
            &report.source,
            controls.total_energy_tolerance_j,
        )
        .map_err(heat_error)?;
        check(&cancel)?;
        Ok(next)
    }

    /// Close nonlinear material current, thermal contacts and radiation together.
    ///
    /// Each trial resolves the original source sigma/rho(T, phase, |E|, |J|)
    /// through its native spatial field owner and re-evaluates this complete heat
    /// assembly at the SAME endpoint enthalpies. Native local E.J is deposited
    /// once. Interfaces, emissivity and optional source-integrated layer fluxes
    /// retain their original owners, selected claims and receiving cells.
    ///
    /// Voltage, distributed current density and equipotential total-current
    /// boundaries keep their original meanings. Every native field iteration,
    /// including discarded thermal trials and final fresh admission, consumes
    /// the common `maximum_circuit_trials` budget. Constitutive defects and
    /// iteration error are never heat. Returned field and caloric histories
    /// must be published together; refusal/cancellation publishes neither.
    ///
    /// Caloric storage, thermal transport, interface resistance and emissivity
    /// must be independent of E/J in both response and source applicability.
    /// No interface electric field or field-dependent radiating temperature is
    /// inferred. Geometry, trace areas/view factors, fixed external conditions
    /// and confined phase-chart densities retain the original heat/field limits.
    /// This is first-order endpoint coupling, not a uniqueness, unconditional
    /// stability, phase-motion or measured-material validation claim.
    /// Closed-enclosure radiation is internal exchange; external heat is
    /// `non_electrical_heat_j - dt * source.radiation_net_power_w()`.
    #[allow(clippy::too_many_arguments)]
    pub fn step_conduction_response_candidate<'t>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &PlanarConductionResponse,
        drive: ConductionDrive<'_>,
        regions: &[ConductionRegion],
        volume_tolerance_m3: f64,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        response: ConductionResponseControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
    ) -> Result<ConductionResponseEnthalpyCandidate<'t, HeatNetworkSample>, Error> {
        self.step_conduction_response_contacts_candidate(
            cx,
            thermal,
            target_s,
            problem,
            drive,
            &[],
            regions,
            volume_tolerance_m3,
            assembly,
            field_controls,
            response,
            controls,
            initial_energy_increments_j,
        )
    }

    /// Nonlinear source-driven bulk current and source-refreshed electrical
    /// contacts share this complete thermal-contact/radiation network. Native
    /// volume and surface Joule losses enter their original receiving cells
    /// exactly once; numerical constitutive defects never become heat.
    #[allow(clippy::too_many_arguments)]
    pub fn step_conduction_response_contacts_candidate<'t>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &PlanarConductionResponse,
        drive: ConductionDrive<'_>,
        contacts: &[contacts::ConductionContactThermal<'_>],
        regions: &[ConductionRegion],
        volume_tolerance_m3: f64,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        response: ConductionResponseControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
    ) -> Result<ConductionResponseEnthalpyCandidate<'t, HeatNetworkSample>, Error> {
        let cancel = RefCell::new(|| cx.checkpoint().is_err());
        check(&cancel)?;
        let bulk = self.thermal_material_network().map_err(heat_error)?;
        let layout = ConductionThermalNetwork {
            cells: bulk.cells,
            links: bulk.links,
            regions,
            volume_tolerance_m3,
        };
        if regions.len() != problem.mesh().faces().len() || regions.len() > assembly.max_cells {
            return Err(invalid(
                "one bounded heat route per original electrical triangle required",
            ));
        }
        let mut field_axes = Vec::new();
        for region in regions {
            check(&cancel)?;
            let source = problem
                .sources()
                .get(region.triangle)
                .ok_or_else(|| invalid("thermal route names an absent spatial conductor"))?;
            for (axis, _) in source.axes() {
                field_axes.push((region.cell, axis.as_str()));
            }
        }
        self.admit_field_independent_boundaries(cx, &field_axes)
            .map_err(heat_error)?;
        let next = layout.step_response_resolved_candidate(
            thermal,
            target_s,
            problem,
            drive,
            contacts,
            assembly,
            field_controls,
            response,
            controls,
            initial_energy_increments_j,
            |trial| {
                let coefficients = self.resolve(cx, trial).map_err(heat_error)?;
                if !coefficients.contacts.is_empty() || !coefficients.bulk.resistors.is_empty() {
                    return Err(invalid(
                        "spatial current cannot silently omit graph electrical work",
                    ));
                }
                Ok(coefficients.bulk)
            },
            |sample| &sample.material,
            &cancel,
        )?;
        let report = next.report();
        closed_balance(
            target_s - thermal.time_s(),
            report.total_energy_residual_j,
            &report.source,
            controls.total_energy_tolerance_j,
        )
        .map_err(heat_error)?;
        check(&cancel)?;
        Ok(next)
    }

    /// Joint spatial current, material-temperature/phase feedback, thermal
    /// contact resistance and closed-enclosure radiation.
    ///
    /// The field's actual E.J integral heats the explicitly mapped original
    /// cells once. At every nonlinear trial the existing prepared heat network
    /// re-evaluates pinned bulk/interface/emissivity sources, and the existing
    /// field owner re-evaluates its pinned scalar/oriented conductivity at that
    /// SAME thermal state. No graph resistor, second timeline or heat law is
    /// introduced. Current-controlled electrodes still solve their voltage
    /// from the full spatial operator on every trial.
    ///
    /// `regions` and `volume_tolerance_m3` have the original
    /// [`ConductionThermalNetwork`] whole-volume semantics. Material/caloric
    /// identities, fixed coordinates, constant phase-chart density, and source
    /// support retain the original strict checks. The heat assembly must have
    /// no electrical graph bindings: those require their own circuit owner.
    /// A thermal-interface resistance is not an electrical contact resistance.
    ///
    /// Geometry, enclosure view factors and contact/radiating areas remain
    /// fixed. This does not model material motion, volume-changing melting or
    /// evolving interface topology. In particular a non-solid phase chart must
    /// meet the same explicitly confined constant-density assumptions as the
    /// original spatial conduction lane; no changing surface is inferred.
    ///
    /// The returned report retains the complete material/contact/radiation
    /// samples. Closed-enclosure radiation is internal exchange: its signed
    /// numerical net is separately checked against the ORIGINAL total-energy
    /// allowance. The physically external non-electrical heat is
    /// `report.non_electrical_heat_j - dt * report.source.radiation_net_power_w()`.
    /// First-order backward-Euler coupling requires time/spatial refinement;
    /// this composition does not establish experimental material validation.
    #[allow(clippy::too_many_arguments)]
    pub fn step_conduction_candidate<'t, C: PlanarConductorMaterial>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &PlanarConduction<C>,
        drive: ConductionDrive<'_>,
        regions: &[ConductionRegion],
        volume_tolerance_m3: f64,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
    ) -> Result<ConductionEnthalpyCandidate<'t, C, HeatNetworkSample>, Error> {
        self.step_conduction_contacts_candidate(
            cx,
            thermal,
            target_s,
            problem,
            drive,
            &[],
            regions,
            volume_tolerance_m3,
            assembly,
            field_controls,
            controls,
            initial_energy_increments_j,
        )
    }

    /// Advance source-backed electrical contacts through this SAME complete
    /// thermal-contact/radiation network. Surface voltage-jump losses enter
    /// their actual adjacent material cells under each explicit heat fraction.
    /// Interface and bulk electrical laws refresh at the same endpoint as all
    /// heat-transport sources. No separate thermal step or loss account is added.
    /// All restrictions of [`Self::step_conduction_candidate`] and the native
    /// stationary contact-field owner remain; graph electrical work requires
    /// its own circuit owner and cannot be omitted from this field-only lane.
    #[allow(clippy::too_many_arguments)]
    pub fn step_conduction_contacts_candidate<'t, C: PlanarConductorMaterial>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        problem: &PlanarConduction<C>,
        drive: ConductionDrive<'_>,
        contacts: &[contacts::ConductionContactThermal<'_>],
        regions: &[ConductionRegion],
        volume_tolerance_m3: f64,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
    ) -> Result<ConductionEnthalpyCandidate<'t, C, HeatNetworkSample>, Error> {
        let cancel = RefCell::new(|| cx.checkpoint().is_err());
        check(&cancel)?;
        let bulk = self.thermal_material_network().map_err(heat_error)?;
        let layout = ConductionThermalNetwork {
            cells: bulk.cells,
            links: bulk.links,
            regions,
            volume_tolerance_m3,
        };
        let next = layout.step_resolved_candidate(
            thermal,
            target_s,
            problem,
            drive,
            contacts,
            assembly,
            field_controls,
            controls,
            initial_energy_increments_j,
            |trial| {
                let coefficients = self.resolve(cx, trial).map_err(heat_error)?;
                // The original assembly has already excluded graph bindings;
                // preserve a fail-closed boundary if another owner is extended.
                if !coefficients.contacts.is_empty() || !coefficients.bulk.resistors.is_empty() {
                    return Err(invalid(
                        "spatial current cannot silently omit graph electrical work",
                    ));
                }
                Ok(coefficients.bulk)
            },
            |sample| &sample.material,
            &cancel,
        )?;
        let report = next.report();
        closed_balance(
            target_s - thermal.time_s(),
            report.total_energy_residual_j,
            &report.source,
            controls.total_energy_tolerance_j,
        )
        .map_err(heat_error)?;
        check(&cancel)?;
        Ok(next)
    }
}
