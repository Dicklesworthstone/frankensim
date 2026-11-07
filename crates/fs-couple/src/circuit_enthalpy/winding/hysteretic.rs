//! Native winding circuits, vector magnetic memory and finite heat storage.
//!
//! Every thermal trial restarts the SAME accepted circuit, field potential and
//! branch polarizations. Original source cards supply winding/core conductivity
//! and optional graph resistance at the proposed thermal endpoint. Spatial
//! pinning, core eddy, stranded copper and graph losses reach their original
//! cells once; numerical field/end loss remains separate from physical heat.
//! The original fitted magnetic coefficients and stationary geometry stay fixed.

use super::{
    CircuitEnthalpyError, CircuitEnthalpyReport, CircuitEnthalpySolve, DrivenCoefficients,
    DriverAccount, Error, ThermalDriver, WindingThermalNetwork, WindingThermalSample,
    cell_driver_candidate, check, invalid, material, nonlinear::electrical, same_material,
};
use crate::circuit_enthalpy::{
    CircuitThermalCoefficients,
    eddy::hysteretic::deposit_constitutive_heat,
    material::{CellResistor, MaterialNetwork, MaterialNetworkSample},
};
use fs_em::{
    planar::eddy::EddyConductor,
    winding::{
        hysteretic::{
            HystereticWindingCandidate, HystereticWindingCircuit, HystereticWindingControls,
            HystereticWindingError, HystereticWindingStep,
        },
        nonlinear::resistors::GraphResistorMaterial,
        ohmic::WindingConductor,
    },
};
use fs_material::conductor::ElectricalPropertyBasis;
use fs_phs::circuit::Branch;
use fs_thermal::enthalpy::DepositionSimulation;
use std::{cell::RefCell, collections::BTreeMap};

/// One circuit, spatial vector-memory field and finite caloric history.
pub struct HystereticWindingEnthalpyCandidate<'t, 'm, S = MaterialNetworkSample> {
    thermal: DepositionSimulation<'t>,
    circuit: HystereticWindingCircuit<'m>,
    report: CircuitEnthalpyReport<WindingThermalSample<S>, HystereticWindingStep>,
}

impl<'t, 'm, S> HystereticWindingEnthalpyCandidate<'t, 'm, S> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn circuit(&self) -> &HystereticWindingCircuit<'m> {
        &self.circuit
    }
    pub fn report(&self) -> &CircuitEnthalpyReport<WindingThermalSample<S>, HystereticWindingStep> {
        &self.report
    }
    /// Commit all histories together. Physical pinning and Joule heat is
    /// already deposited; numerical loss must not be deposited again.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        HystereticWindingCircuit<'m>,
        CircuitEnthalpyReport<WindingThermalSample<S>, HystereticWindingStep>,
    ) {
        (self.thermal, self.circuit, self.report)
    }
}

impl<'m> ThermalDriver for HystereticWindingCandidate<'m> {
    type Model = HystereticWindingCircuit<'m>;
    type Step = HystereticWindingStep;

    fn account(&self) -> DriverAccount<'_> {
        let step = self.step();
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: step.heat_j,
            numerical_dissipation_j: step.numerical_dissipation_j,
            // Preserve small resolved work beside a large magnetic bias.
            // This includes stable native graph, field and end-mode changes.
            stored_change_j: step.stored_change_j,
            // Private winding-port work cancels inside the native field/circuit
            // owner. Fixed fitted coefficients add no magnetic parameter work.
            supplied_j: step.supplied_work_j,
        }
    }

    fn into_parts(self) -> (Self::Model, Self::Step) {
        HystereticWindingCandidate::into_parts(self)
    }
}

impl WindingThermalNetwork<'_> {
    /// Evolve native winding current, remanence, core diffusion and finite
    /// material heat on the same original enthalpy residual.
    ///
    /// All thermal/Newton trials start from the same accepted circuit and
    /// vector-memory field. Conductivity follows each proposed cell temperature;
    /// fitted magnetic coefficients must remain constant on their source-admitted
    /// temperature intervals. A dissipating magnetic triangle needs an explicit
    /// finite receiver, including electrically insulating material. Core eddy,
    /// stranded copper, native graph and magnetic pinning/viscous heat are deposited once.
    /// Backward-Euler numerical loss participates only in the work balance.
    ///
    /// Geometry and mapped wholly solid, constant-density caloric charts remain
    /// fixed. This is first-order coupled stepping. Missing source support,
    /// cancellation or exhausted nested field budgets publish no partial history.
    #[allow(clippy::too_many_arguments)]
    pub fn step_hysteretic_candidate<'t, 'm>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &HystereticWindingCircuit<'m>,
        external_inputs: &[f64],
        controls: CircuitEnthalpySolve,
        winding_controls: HystereticWindingControls,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<HystereticWindingEnthalpyCandidate<'t, 'm>, Error> {
        if !circuit.graph_resistor_materials().is_empty() {
            return Err(invalid(
                "source-bound graph resistors require their complete thermal material bindings",
            ));
        }
        self.step_hysteretic_bound(
            thermal,
            target_s,
            circuit,
            external_inputs,
            None,
            controls,
            winding_controls,
            initial_energy_increments_j,
            cancelled,
        )
    }

    /// The same circuit/field/enthalpy transaction with every ordinary graph
    /// resistor resolved from its original whole-cell material and endpoint.
    ///
    /// Bindings must match all existing graph heat destinations exactly. Their
    /// whole-cell volumes remain separate from the magnetic geometry. Accepted
    /// native resistance, retained source choices and conductor geometry must
    /// agree with the original thermal state before any trial is attempted.
    /// Resistivity and conductivity sources use the existing material owner;
    /// changing either electrical response does not reset magnetic memory.
    #[allow(clippy::too_many_arguments)]
    pub fn step_hysteretic_material_candidate<'t, 'm>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &HystereticWindingCircuit<'m>,
        external_inputs: &[f64],
        resistors: &[CellResistor],
        electrical_bases: &[(usize, ElectricalPropertyBasis)],
        controls: CircuitEnthalpySolve,
        winding_controls: HystereticWindingControls,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<HystereticWindingEnthalpyCandidate<'t, 'm>, Error> {
        self.step_hysteretic_bound(
            thermal,
            target_s,
            circuit,
            external_inputs,
            Some((resistors, electrical_bases)),
            controls,
            winding_controls,
            initial_energy_increments_j,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_hysteretic_bound<'t, 'm>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &HystereticWindingCircuit<'m>,
        external_inputs: &[f64],
        graph_material: Option<(&[CellResistor], &[(usize, ElectricalPropertyBasis)])>,
        controls: CircuitEnthalpySolve,
        winding_controls: HystereticWindingControls,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<HystereticWindingEnthalpyCandidate<'t, 'm>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        crate::circuit_enthalpy::admit_driver_unknowns(
            thermal,
            controls,
            initial_energy_increments_j,
        )?;
        if circuit.time_s().to_bits() != thermal.time_s().to_bits()
            || circuit.field().time_s().to_bits() != thermal.time_s().to_bits()
            || self.cells.len() != thermal.cells().len()
        {
            return Err(invalid(
                "matching hysteretic field/circuit/thermal clocks and complete cell bindings required",
            ));
        }
        let network = MaterialNetwork {
            cells: self.cells,
            resistors: graph_material.map_or(&[], |(resistors, _)| resistors),
            links: self.links,
        };
        let prepared = network
            .prepare_with_bases(
                thermal,
                graph_material.map_or(&[], |(_, bases)| bases),
                || (*cancel.borrow_mut())(),
            )
            .map_err(material)?;
        self.step_hysteretic_resolved(
            thermal,
            target_s,
            circuit,
            external_inputs,
            graph_material.map(|(resistors, _)| resistors),
            controls,
            winding_controls,
            initial_energy_increments_j,
            |state| {
                prepared
                    .resolve(state, || (*cancel.borrow_mut())())
                    .map_err(material)
            },
            |sample| sample,
            &cancel,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_hysteretic_resolved<'t, 'm, S>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &HystereticWindingCircuit<'m>,
        external_inputs: &[f64],
        graph_resistors: Option<&[CellResistor]>,
        controls: CircuitEnthalpySolve,
        winding_controls: HystereticWindingControls,
        initial_energy_increments_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'_>) -> Result<CircuitThermalCoefficients<S>, Error>,
        sample: impl Fn(&S) -> &MaterialNetworkSample,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<HystereticWindingEnthalpyCandidate<'t, 'm, S>, Error> {
        check(cancel)?;
        crate::circuit_enthalpy::admit_driver_unknowns(
            thermal,
            controls,
            initial_energy_increments_j,
        )?;
        if circuit.time_s().to_bits() != thermal.time_s().to_bits()
            || circuit.field().time_s().to_bits() != thermal.time_s().to_bits()
            || self.cells.len() != thermal.cells().len()
        {
            return Err(invalid(
                "matching hysteretic field/circuit/thermal clocks and complete cell bindings required",
            ));
        }
        if graph_resistors.is_none() && !circuit.graph_resistor_materials().is_empty() {
            return Err(invalid(
                "source-bound graph resistors require their complete thermal material bindings",
            ));
        }
        let mut conductors = BTreeMap::new();
        let mut core_regions = BTreeMap::new();
        for source in circuit.field().conductors() {
            check(cancel)?;
            if conductors
                .insert(source.triangle, &source.conductor)
                .is_some()
                || core_regions
                    .insert(source.triangle, source.region)
                    .is_some()
            {
                return Err(invalid("duplicate original core conductor triangle"));
            }
        }
        for source in circuit.conductors() {
            check(cancel)?;
            if conductors
                .insert(source.triangle, &source.conductor)
                .is_some()
            {
                return Err(invalid(
                    "stranded winding and core conductor triangles overlap",
                ));
            }
        }
        let triangles: Vec<_> = conductors.keys().copied().collect();
        let routes = self.admit_routes_for_mesh(
            thermal,
            circuit.field().problem().mesh(),
            &triangles,
            circuit.graph(),
            cancel,
        )?;
        let layout = self.magnetic_layout();
        layout.admit_hysteretic_temperature_support(
            circuit.field().problem().materials(),
            &routes.magnetic,
            cancel,
        )?;
        if let Some(resistors) = graph_resistors {
            if resistors.len() != routes.graph.len()
                || resistors.iter().any(|source| {
                    routes.graph.get(&source.destination.branch) != Some(&source.destination.cell)
                        || routes
                            .magnetic
                            .magnetic
                            .values()
                            .any(|cell| *cell == Some(source.destination.cell))
                })
            {
                return Err(invalid(
                    "sourced graph resistors need exact existing heat routes and distinct whole-cell volumes",
                ));
            }
        }
        let accepted = resolve(thermal)?;
        if graph_resistors.is_some() {
            for source in &accepted.resistors {
                check(cancel)?;
                let Some((_, _, Branch::Resistor { ohms })) =
                    circuit.graph().branches.get(source.branch)
                else {
                    return Err(invalid("material graph source names a non-resistor branch"));
                };
                if ohms.to_bits() != source.conductor.resistance_ohm().to_bits() {
                    return Err(invalid(
                        "accepted graph resistance does not match its actual thermal source",
                    ));
                }
                if !circuit.graph_resistor_materials().is_empty()
                    && !circuit.graph_resistor_materials().iter().any(|old| {
                        old.branch == source.branch
                            && old.conductor.source_basis() == source.conductor.source_basis()
                            && old.conductor.source_key() == source.conductor.source_key()
                            && old.conductor.length_m().to_bits()
                                == source.conductor.length_m().to_bits()
                            && old.conductor.area_m2().to_bits()
                                == source.conductor.area_m2().to_bits()
                            && same_material(old.conductor.material(), source.conductor.material())
                    })
                {
                    return Err(invalid(
                        "graph resistor changed its accepted material source or geometry",
                    ));
                }
            }
        }
        let resolve_materials = |state: &DepositionSimulation<'_>,
                                 source: &MaterialNetworkSample| {
            layout.hysteretic_materials_at(
                state,
                circuit.field().problem().materials(),
                &conductors,
                source,
                &routes.magnetic,
                cancel,
            )
        };
        let (initial, _) = resolve_materials(thermal, sample(&accepted.source))?;
        if initial.len() != conductors.len()
            || initial.iter().any(|(triangle, source)| {
                conductors.get(triangle).is_none_or(|original| {
                    source.source_basis() != original.source_basis()
                        || source.source_key() != original.source_key()
                        || source.resistivity_ohm_m().to_bits()
                            != original.resistivity_ohm_m().to_bits()
                        || !same_material(source.material(), original.material())
                })
            })
        {
            return Err(invalid(
                "accepted hysteretic winding/core sources must match the actual thermal endpoint",
            ));
        }
        let (thermal, circuit, report) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            initial_energy_increments_j,
            |state, remaining| {
                check(cancel)?;
                // Bound all nested native field trials before starting any
                // circuit solve, including discarded stencils and final replay.
                let mut inner = winding_controls;
                inner.maximum_field_trials = inner.maximum_field_trials.min(remaining);
                if inner.maximum_field_trials == 0 {
                    return Err(CircuitEnthalpyError::CircuitBudget);
                }
                let outer_limited = remaining < winding_controls.maximum_field_trials;
                let coefficients = resolve(state)?;
                let (sources, magnetic) = resolve_materials(state, sample(&coefficients.source))?;
                let mut core = Vec::new();
                let mut winding = Vec::new();
                for (triangle, conductor) in sources {
                    check(cancel)?;
                    if let Some(&region) = core_regions.get(&triangle) {
                        core.push(EddyConductor {
                            triangle,
                            region,
                            conductor,
                        });
                    } else {
                        winding.push(WindingConductor {
                            triangle,
                            conductor,
                        });
                    }
                }
                // Rebind only endpoint source coefficients on a clone of the
                // original accepted circuit. Never advance from a thermal trial.
                let rebound;
                let driven = if graph_resistors.is_some() {
                    rebound = circuit
                        .with_graph_resistors(
                            coefficients
                                .resistors
                                .iter()
                                .map(|source| GraphResistorMaterial {
                                    branch: source.branch,
                                    conductor: source.conductor.clone(),
                                })
                                .collect(),
                            &mut || (*cancel.borrow_mut())(),
                        )
                        .map_err(electrical)?;
                    &rebound
                } else {
                    circuit
                };
                let proposal = driven
                    .step_to_candidate(target_s, core, winding, external_inputs, inner, &mut || {
                        (*cancel.borrow_mut())()
                    })
                    .map_err(|error| match error {
                        HystereticWindingError::TrialBudget if outer_limited => {
                            CircuitEnthalpyError::CircuitBudget
                        }
                        other => electrical(other),
                    })?;
                let circuit_trials = proposal.step().field_trials;
                Ok(DrivenCoefficients {
                    circuit: proposal,
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    source: WindingThermalSample {
                        network: coefficients.source,
                        magnetic,
                        magnetic_caloric: None,
                        magnetic_balance: None,
                    },
                    circuit_trials,
                })
            },
            |proposal: &HystereticWindingCandidate<'_>, dt, checkpoint| {
                let step = proposal.step();
                let mut powers = routes.deposit_losses(
                    &step.cells,
                    &step.field.cells,
                    &step.graph_resistors,
                    dt,
                    checkpoint,
                )?;
                deposit_constitutive_heat(
                    &step.field.magnetic.cells,
                    &routes.magnetic,
                    dt,
                    &mut powers,
                    checkpoint,
                )?;
                Ok(powers)
            },
            || (*cancel.borrow_mut())(),
        )?;
        check(cancel)?;
        Ok(HystereticWindingEnthalpyCandidate {
            thermal,
            circuit,
            report,
        })
    }
}
