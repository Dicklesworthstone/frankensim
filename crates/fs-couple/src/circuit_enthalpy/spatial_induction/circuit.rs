//! Native RLC circuits, solid-conductor induction and original finite enthalpy.
//!
//! Every thermal proposal refreshes the original conductor and magnetic sources,
//! then closes the actual electrode field against the original circuit history.
//! Conductor-frame Joule heat and authored resistor heat enter explicit finite
//! receivers once. Magnetic entropy/parameter work stays with the existing
//! caloric owner; numerical reference resistors never become physical heat.
//!
//! The native circuit graph, its R/L/C coefficients, and geometry remain fixed.
//! Endpoint field coefficients and backward-Euler enthalpy give first-order
//! coupled accuracy in general. Prescribed motion transfers external mechanical
//! work; it does not move the material interface or advance a mechanical body.

use super::terminals::circuit::{
    SpatialTerminalCircuitThermalCandidate, SpatialTerminalCircuitThermalNetwork,
};
use super::*;
use crate::circuit_enthalpy::ResistorCell;
use fs_em::spatial::vector_potential::induction::terminals::circuit::{
    SpatialTerminalCircuit, SpatialTerminalCircuitCandidate, SpatialTerminalCircuitControls,
    SpatialTerminalCircuitError, SpatialTerminalCircuitStep,
};

/// The shared accepted electrode-circuit/caloric candidate, specialized to the
/// linear magnetic field. No trial updates any of its three histories.
pub type SpatialCircuitThermalCandidate<
    't,
    'm,
    C = ResolvedConductorState,
    S = MaterialNetworkSample,
> = SpatialTerminalCircuitThermalCandidate<
    't,
    SpatialTerminalCircuit<'m, C>,
    SpatialTerminalCircuitStep,
    S,
>;

impl<C: SpatialConductorMaterial> LinearInductionHistory<C> for SpatialTerminalCircuit<'_, C> {
    fn time_s(&self) -> f64 {
        SpatialTerminalCircuit::time_s(self)
    }
    fn problem(&self) -> &SpatialVectorMagnetostatics {
        self.field().problem()
    }
    fn conductors(&self) -> &[Option<C>] {
        self.field().conductors()
    }
    fn flux_density_t(
        &self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<[f64; 3]>, SpatialMagneticError> {
        Ok(self.field().field(cancelled)?.flux_density_t().to_vec())
    }
}

impl<'m, C: SpatialConductorMaterial> ThermalDriver for SpatialTerminalCircuitCandidate<'m, C> {
    type Model = SpatialTerminalCircuit<'m, C>;
    type Step = SpatialTerminalCircuitStep;

    fn account(&self) -> DriverAccount<'_> {
        let step = self.step();
        let field = &step.field;
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: field.heat_j + step.resistor_heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: field.stored_after_j - field.stored_before_j - field.parameter_work_j
                + step.circuit_stored_change_j,
            supplied_j: field.impressed_work_j - field.mechanical_work_j + step.external_work_j,
        }
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        SpatialTerminalCircuitCandidate::into_parts(self)
    }
}

impl<C: SpatialConductorMaterial> LinearInductionThermalDriver
    for SpatialTerminalCircuitCandidate<'_, C>
{
    fn induction_step(&self) -> &SpatialInductionStep {
        &self.step().field
    }
    fn induction_receipt(step: &Self::Step) -> &SpatialInductionStep {
        &step.field
    }
    fn field_trials(&self) -> usize {
        self.step().field_trials
    }
    fn additional_stored_change_j(step: &Self::Step) -> f64 {
        step.circuit_stored_change_j
    }
    fn additional_supplied_j(step: &Self::Step) -> f64 {
        step.external_work_j
    }
}

fn circuit_error(error: SpatialTerminalCircuitError) -> Error {
    match error {
        SpatialTerminalCircuitError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(MagneticThermalError::SpatialTerminalCircuit(other)),
    }
}

impl SpatialMagneticThermalNetwork<'_> {
    /// Solve original field, circuit and finite-temperature feedback together.
    ///
    /// Every tetrahedral conductor retains its original full-volume caloric
    /// receiver. Every original graph resistor needs one explicit destination;
    /// its fixed authored resistance is retained. Passive thermal cells outside
    /// the field may receive circuit loss and participate in heat transport.
    /// The original magnetic thermal laws own field entropy, effective-capacity
    /// admission and reversible heat separately from both Ohmic contributions.
    ///
    /// All field-response probes count against maximum_circuit_trials, including
    /// final replay, and each native field solve retains its independent bounds.
    /// Every proposal starts from the same accepted A, circuit fluxes/charges
    /// and caloric history. Publication returns all three endpoints together.
    #[allow(clippy::too_many_arguments)]
    pub fn step_terminal_circuit_candidate<'t, 'm, C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &SpatialTerminalCircuit<'m, C>,
        external_inputs: &[f64],
        endpoint_impressed_current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        resistor_destinations: &[ResistorCell],
        bindings: &[SpatialMagneticThermalLaw],
        controls: CircuitEnthalpySolve,
        circuit_controls: SpatialTerminalCircuitControls,
        assembly: WeightedAssemblyLimits,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialCircuitThermalCandidate<'t, 'm, C>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        super::super::admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
        let receiving = SpatialTerminalCircuitThermalNetwork {
            magnetic: SpatialMagneticThermalNetwork {
                cells: self.cells,
                links: self.links,
                magnetic_regions: self.magnetic_regions,
                volume_tolerance_m3: self.volume_tolerance_m3,
            },
            graph_resistors: resistor_destinations,
        };
        let destinations = receiving.routes(thermal, circuit.time_s(), circuit.graph(), &cancel)?;
        let material_network = MaterialNetwork {
            cells: self.cells,
            resistors: &[],
            links: self.links,
        };
        let prepared = material_network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        let (thermal, circuit, report) = self.step_linear_driver_candidate(
            thermal,
            target_s,
            circuit,
            bindings,
            controls,
            initial_energy_increments_j,
            |trial| {
                prepared
                    .resolve(trial, || (*cancel.borrow_mut())())
                    .map_err(material)
            },
            |sample| sample,
            |_, _, materials, conductors, remaining| {
                let mut limited = circuit_controls;
                limited.maximum_field_trials = limited.maximum_field_trials.min(remaining);
                circuit
                    .step_with_materials_candidate(
                        external_inputs,
                        endpoint_impressed_current_a_m2,
                        velocity_m_s,
                        target_s - thermal.time_s(),
                        materials,
                        conductors,
                        assembly,
                        limited,
                        &mut || (*cancel.borrow_mut())(),
                    )
                    .map_err(|error| {
                        if matches!(error, SpatialTerminalCircuitError::FieldTrialBudget)
                            && remaining < circuit_controls.maximum_field_trials
                        {
                            CircuitEnthalpyError::CircuitBudget
                        } else {
                            circuit_error(error)
                        }
                    })
            },
            |proposal, powers, dt, checkpoint| {
                for (&branch, &cell) in &destinations {
                    checkpoint()?;
                    let loss = proposal
                        .step()
                        .branches
                        .get(branch)
                        .filter(|loss| loss.branch == branch)
                        .ok_or_else(|| invalid("native circuit omitted an original resistor"))?;
                    let heat = loss.absorbed_work_j;
                    let power = finite(heat / dt)?;
                    if heat < 0.0 || heat != 0.0 && power == 0.0 {
                        return Err(invalid(
                            "circuit heat is negative or lost in power conversion",
                        ));
                    }
                    let previous = *powers.get(&cell).unwrap_or(&0.0);
                    let next = finite(previous + power)?;
                    if power != 0.0 && next == previous {
                        return Err(invalid("circuit heat disappeared in its thermal receiver"));
                    }
                    powers.insert(cell, next);
                }
                Ok(())
            },
            &cancel,
        )?;
        Ok(SpatialCircuitThermalCandidate {
            thermal,
            circuit,
            report,
        })
    }
}
