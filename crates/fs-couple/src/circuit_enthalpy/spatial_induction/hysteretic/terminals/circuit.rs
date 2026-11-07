//! Native circuit-driven solid conductors, vector magnetic memory and finite heat.
//!
//! Graph storage, electrode reactions, A–phi and vector play all restart from
//! accepted histories on every thermal trial. The original field thermal driver
//! owns core heat; only authored graph resistor losses are added here. Numerical
//! reference loss and numerical magnetic dissipation never heat a finite body.

pub use super::super::super::terminals::circuit::{
    SpatialTerminalCircuitThermalCandidate, SpatialTerminalCircuitThermalError,
    SpatialTerminalCircuitThermalNetwork,
};
use super::super::super::{
    SpatialMagneticRegion,
    terminals::circuit::{
        TerminalError, admit_initial_graph_sources, close_heat_balance, deposit_graph_heat,
        electrical, heat_error, pinned_graph_sources,
    },
};
use super::super::{
    CircuitEnthalpyError, CircuitEnthalpySolve, CircuitThermalCoefficients, DriverAccount,
    HystereticInductionDriver, MaterialNetwork, MaterialNetworkSample, RefCell,
    SpatialHystereticInductionStep, SpatialMagneticThermalNetwork, ThermalDriver, check, invalid,
    material, native_error,
};
use crate::circuit_enthalpy::{ResistorCell, contacts::heat::HeatAssembly, material::CellResistor};
use fs_em::spatial::vector_potential::induction::{
    hysteretic::terminals::circuit::{
        SpatialHystereticTerminalCircuit, SpatialHystereticTerminalCircuitCandidate,
        SpatialHystereticTerminalCircuitControls, SpatialHystereticTerminalCircuitStep,
    },
    terminals::circuit::SpatialTerminalCircuitError,
};
use fs_exec::Cx;
use fs_material::conductor::{ElectricalPropertyBasis, spatial::SpatialConductorMaterial};
use fs_thermal::enthalpy::DepositionSimulation;

impl<'m, C: SpatialConductorMaterial> ThermalDriver
    for SpatialHystereticTerminalCircuitCandidate<'m, C>
{
    type Model = SpatialHystereticTerminalCircuit<'m, C>;
    type Step = SpatialHystereticTerminalCircuitStep;

    fn account(&self) -> DriverAccount<'_> {
        Self::receipt_account(self.step())
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        SpatialHystereticTerminalCircuitCandidate::into_parts(self)
    }
}

impl<C: SpatialConductorMaterial> HystereticInductionDriver
    for SpatialHystereticTerminalCircuitCandidate<'_, C>
{
    fn induction_step(&self) -> &SpatialHystereticInductionStep {
        &self.step().field
    }
    fn receipt_account(step: &Self::Step) -> DriverAccount<'_> {
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: step.heat_j,
            numerical_dissipation_j: step.numerical_dissipation_j,
            stored_change_j: step.stored_change_j,
            supplied_j: step.supplied_work_j,
        }
    }
    fn field_trials(&self) -> usize {
        self.step().field_trials
    }
}

impl SpatialTerminalCircuitThermalNetwork<'_> {
    /// Advance native RLC storage, physical electrode current, vector memory and
    /// finite caloric inventory together. The field determines bulk current and
    /// conductor loss; authored graph resistance retains its fixed value.
    ///
    /// Each original graph resistor requires a heat destination. Magnetic fit
    /// coefficients remain fixed over the complete admitted caloric interval;
    /// core conductivity follows the actual thermal endpoint. Every nested
    /// circuit/field trial consumes the existing shared thermal solve budget.
    #[allow(clippy::too_many_arguments)]
    pub fn step_hysteretic_candidate<'t, 'm, C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &SpatialHystereticTerminalCircuit<'m, C>,
        external_inputs: &[f64],
        endpoint_impressed_current_a_m2: &[[f64; 3]],
        controls: CircuitEnthalpySolve,
        circuit_controls: SpatialHystereticTerminalCircuitControls,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        SpatialTerminalCircuitThermalCandidate<
            't,
            SpatialHystereticTerminalCircuit<'m, C>,
            SpatialHystereticTerminalCircuitStep,
        >,
        TerminalError,
    > {
        self.step_hysteretic_bound(
            thermal,
            target_s,
            circuit,
            external_inputs,
            endpoint_impressed_current_a_m2,
            None,
            controls,
            circuit_controls,
            initial_energy_increments_j,
            cancelled,
        )
    }

    /// Additionally resolve every original graph resistor from its own finite
    /// material body. Resistance refresh preserves circuit fluxes, charges,
    /// topology, source choices and the full accepted vector magnetic history.
    #[allow(clippy::too_many_arguments)]
    pub fn step_hysteretic_material_candidate<'t, 'm, C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &SpatialHystereticTerminalCircuit<'m, C>,
        external_inputs: &[f64],
        endpoint_impressed_current_a_m2: &[[f64; 3]],
        resistors: &[CellResistor],
        electrical_bases: &[(usize, ElectricalPropertyBasis)],
        controls: CircuitEnthalpySolve,
        circuit_controls: SpatialHystereticTerminalCircuitControls,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        SpatialTerminalCircuitThermalCandidate<
            't,
            SpatialHystereticTerminalCircuit<'m, C>,
            SpatialHystereticTerminalCircuitStep,
        >,
        TerminalError,
    > {
        self.step_hysteretic_bound(
            thermal,
            target_s,
            circuit,
            external_inputs,
            endpoint_impressed_current_a_m2,
            Some((resistors, electrical_bases)),
            controls,
            circuit_controls,
            initial_energy_increments_j,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_hysteretic_bound<'t, 'm, C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &SpatialHystereticTerminalCircuit<'m, C>,
        external_inputs: &[f64],
        impressed: &[[f64; 3]],
        graph_material: Option<(&[CellResistor], &[(usize, ElectricalPropertyBasis)])>,
        controls: CircuitEnthalpySolve,
        circuit_controls: SpatialHystereticTerminalCircuitControls,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        SpatialTerminalCircuitThermalCandidate<
            't,
            SpatialHystereticTerminalCircuit<'m, C>,
            SpatialHystereticTerminalCircuitStep,
        >,
        TerminalError,
    > {
        let cancel = RefCell::new(cancelled);
        check(&cancel).map_err(native_error)?;
        crate::circuit_enthalpy::admit_driver_unknowns(
            thermal,
            controls,
            initial_energy_increments_j,
        )?;
        let network = MaterialNetwork {
            cells: self.magnetic.cells,
            resistors: graph_material.map_or(&[], |(resistors, _)| resistors),
            links: self.magnetic.links,
        };
        let prepared = network
            .prepare_with_bases(
                thermal,
                graph_material.map_or(&[], |(_, bases)| bases),
                || (*cancel.borrow_mut())(),
            )
            .map_err(material)
            .map_err(native_error)?;
        self.step_hysteretic_resolved_candidate(
            thermal,
            target_s,
            circuit,
            external_inputs,
            impressed,
            graph_material.map(|(resistors, _)| resistors),
            controls,
            circuit_controls,
            initial_energy_increments_j,
            |state| {
                prepared
                    .resolve(state, || (*cancel.borrow_mut())())
                    .map_err(material)
                    .map_err(native_error)
            },
            |sample| sample,
            &cancel,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_hysteretic_resolved_candidate<'t, 'm, S, E, C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &SpatialHystereticTerminalCircuit<'m, C>,
        external_inputs: &[f64],
        impressed: &[[f64; 3]],
        graph_material: Option<&[CellResistor]>,
        controls: CircuitEnthalpySolve,
        circuit_controls: SpatialHystereticTerminalCircuitControls,
        initial_energy_increments_j: &[f64],
        resolve: impl Fn(
            &DepositionSimulation<'_>,
        ) -> Result<CircuitThermalCoefficients<S>, TerminalError<E>>,
        sample: impl Fn(&S) -> &MaterialNetworkSample,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<
        SpatialTerminalCircuitThermalCandidate<
            't,
            SpatialHystereticTerminalCircuit<'m, C>,
            SpatialHystereticTerminalCircuitStep,
            S,
        >,
        TerminalError<E>,
    > {
        check(cancel).map_err(native_error)?;
        let routes = self
            .routes(thermal, circuit.time_s(), circuit.graph(), cancel)
            .map_err(native_error)?;
        self.admit_graph_material_bindings(
            circuit.graph_resistor_materials(),
            graph_material,
            &routes,
            cancel,
        )
        .map_err(native_error)?;
        let initial = resolve(thermal)?;
        admit_initial_graph_sources(
            circuit.graph(),
            circuit.graph_resistor_materials(),
            graph_material.map_or(0, |bindings| bindings.len()),
            &initial.resistors,
            cancel,
        )
        .map_err(native_error)?;
        let dt = target_s - thermal.time_s();
        let (thermal, circuit, report) = self.magnetic.step_hysteretic_driver_candidate(
            thermal,
            target_s,
            circuit.field(),
            controls,
            initial_energy_increments_j,
            resolve,
            &sample,
            |_, coefficients, conductors, remaining| {
                let mut inner = circuit_controls;
                inner.maximum_field_trials = inner.maximum_field_trials.min(remaining);
                if inner.maximum_field_trials == 0 {
                    return Err(CircuitEnthalpyError::CircuitBudget);
                }
                let refreshed_graph;
                let driven = if graph_material.is_some() {
                    let sources =
                        pinned_graph_sources(&initial.resistors, &coefficients.resistors, cancel)
                            .map_err(native_error)?;
                    refreshed_graph = circuit
                        .with_graph_resistors(sources, &mut || (*cancel.borrow_mut())())
                        .map_err(electrical)?;
                    &refreshed_graph
                } else {
                    if !coefficients.resistors.is_empty() {
                        return Err(native_error(invalid(
                            "unowned graph resistor source appeared in an electrode thermal trial",
                        )));
                    }
                    circuit
                };
                driven
                    .step_with_conductors_candidate(
                        external_inputs,
                        impressed,
                        dt,
                        conductors,
                        inner,
                        &mut || (*cancel.borrow_mut())(),
                    )
                    .map_err(|error| match error {
                        SpatialTerminalCircuitError::FieldTrialBudget
                            if remaining < circuit_controls.maximum_field_trials =>
                        {
                            CircuitEnthalpyError::CircuitBudget
                        }
                        other => electrical(other),
                    })
            },
            |candidate, powers, dt, checkpoint| {
                deposit_graph_heat(
                    &candidate.step().graph_resistors,
                    &routes,
                    powers,
                    dt,
                    checkpoint,
                )
                .map_err(native_error)
            },
            cancel,
        )?;
        Ok(SpatialTerminalCircuitThermalCandidate {
            thermal,
            circuit,
            report,
        })
    }
}

/// Continue the same native circuit, solid-electrode field and vector memory
/// with the complete original heat provider. Thermal contact, layer, radiation
/// and coolant payloads retain their own source and finite-history owners.
#[allow(clippy::too_many_arguments)]
pub fn heat_assembly_candidate<'t, 'm, A: HeatAssembly + ?Sized, C: SpatialConductorMaterial>(
    provider: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    circuit: &SpatialHystereticTerminalCircuit<'m, C>,
    external_inputs: &[f64],
    endpoint_impressed_current_a_m2: &[[f64; 3]],
    regions: &[SpatialMagneticRegion],
    graph_resistors: &[ResistorCell],
    volume_tolerance_m3: f64,
    controls: CircuitEnthalpySolve,
    circuit_controls: SpatialHystereticTerminalCircuitControls,
    initial_energy_increments_j: &[f64],
) -> Result<
    SpatialTerminalCircuitThermalCandidate<
        't,
        SpatialHystereticTerminalCircuit<'m, C>,
        SpatialHystereticTerminalCircuitStep,
        A::Sample,
    >,
    TerminalError<A::Error>,
> {
    let cancel = RefCell::new(|| cx.checkpoint().is_err());
    check(&cancel).map_err(native_error)?;
    crate::circuit_enthalpy::admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
    let prepared = provider.prepared_heat_network();
    let bulk = prepared.field_material_network().map_err(heat_error)?;
    prepared
        .admit_field_independent_boundaries(cx, &[])
        .map_err(heat_error)?;
    let layout = SpatialTerminalCircuitThermalNetwork {
        magnetic: SpatialMagneticThermalNetwork {
            cells: bulk.cells,
            links: bulk.links,
            magnetic_regions: regions,
            volume_tolerance_m3,
        },
        graph_resistors,
    };
    let mut next = layout.step_hysteretic_resolved_candidate(
        thermal,
        target_s,
        circuit,
        external_inputs,
        endpoint_impressed_current_a_m2,
        (!bulk.resistors.is_empty()).then_some(bulk.resistors),
        controls,
        circuit_controls,
        initial_energy_increments_j,
        |state| {
            let coefficients = provider
                .resolve_heat(cx, state)
                .map_err(|error| error.map_source(SpatialTerminalCircuitThermalError::Heat))?;
            if !coefficients.contacts.is_empty() {
                return Err(native_error(invalid(
                    "hysteretic electrode heat cannot discard electrical contact work",
                )));
            }
            Ok(coefficients.bulk)
        },
        |sample: &A::Sample| &sample.as_ref().material,
        &cancel,
    )?;
    close_heat_balance(
        &mut next.report,
        target_s - thermal.time_s(),
        controls.total_energy_tolerance_j,
    )?;
    check(&cancel).map_err(native_error)?;
    Ok(next)
}
