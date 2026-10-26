//! Actual voltage-driven nonlinear fields and their original finite thermal cells.
//!
//! The native circuit and nonlinear A_z histories are solved afresh in every
//! thermal trial. The original source cards supply resistivity at each proposed
//! endpoint; the same magnetic laws, geometry and caloric charts remain pinned.
//! Physical winding/core/graph Joule heat is deposited once. Numerical field
//! dissipation and private circuit-port work are never mistaken for heat.

use super::{
    CircuitEnthalpyError, CircuitEnthalpyReport, CircuitEnthalpySolve, DrivenCoefficients,
    DriverAccount, Error, ThermalDriver, WindingThermalError, WindingThermalNetwork,
    WindingThermalSample, cell_driver_candidate, check, invalid, material, same_material,
};
use crate::circuit_enthalpy::material::MaterialNetwork;
use fs_em::{
    planar::eddy::EddyConductor,
    winding::{
        nonlinear::{
            NonlinearWindingCandidate, NonlinearWindingCircuit, NonlinearWindingControls,
            NonlinearWindingError, NonlinearWindingStep,
        },
        ohmic::WindingConductor,
    },
};
use fs_thermal::enthalpy::DepositionSimulation;
use std::{cell::RefCell, collections::BTreeMap};

/// The circuit, field and finite thermal histories accepted over one interval.
pub struct NonlinearWindingEnthalpyCandidate<'t, 'm> {
    thermal: DepositionSimulation<'t>,
    circuit: NonlinearWindingCircuit<'m>,
    report: CircuitEnthalpyReport<WindingThermalSample, NonlinearWindingStep>,
}
impl<'t, 'm> NonlinearWindingEnthalpyCandidate<'t, 'm> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn circuit(&self) -> &NonlinearWindingCircuit<'m> {
        &self.circuit
    }
    pub fn report(&self) -> &CircuitEnthalpyReport<WindingThermalSample, NonlinearWindingStep> {
        &self.report
    }
    /// Commit all histories together; reported physical heat is already deposited.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        NonlinearWindingCircuit<'m>,
        CircuitEnthalpyReport<WindingThermalSample, NonlinearWindingStep>,
    ) {
        (self.thermal, self.circuit, self.report)
    }
}
impl<'m> ThermalDriver for NonlinearWindingCandidate<'m> {
    type Model = NonlinearWindingCircuit<'m>;
    type Step = NonlinearWindingStep;

    fn account(&self) -> DriverAccount<'_> {
        let step = self.step();
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: step.heat_j,
            numerical_dissipation_j: step.numerical_dissipation_j,
            stored_change_j: step.stored_after_j - step.stored_before_j,
            // Native work also contains private winding-port exchanges. Those
            // cancel against the field/winding work in this combined owner.
            supplied_j: step.supplied_work_j,
        }
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        NonlinearWindingCandidate::into_parts(self)
    }
}
fn electrical(error: NonlinearWindingError) -> Error {
    match error {
        NonlinearWindingError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(WindingThermalError::NonlinearWinding(other)),
    }
}

impl WindingThermalNetwork<'_> {
    /// Solve source-defined saturation, native circuit currents and finite
    /// temperature-dependent Joule heating on the original enthalpy residual.
    ///
    /// Conductivity is re-resolved at the actual proposed thermal endpoint.
    /// Every nested field/circuit solve starts from the same accepted histories.
    /// All conducting triangles need explicit original heat destinations; core
    /// eddy conductors and stranded windings remain disjoint. Every native graph
    /// resistor also deposits its authored physical heat.
    ///
    /// The complete H(B) law and geometry remain fixed; temperature-dependent
    /// magnetic storage refuses until its entropy/parameter work is represented.
    /// Mapped cells use their original wholly solid, constant-density caloric
    /// charts. This is first-order coupled stepping. Backward-Euler numerical
    /// dissipation participates in the global work balance but never heats cells.
    /// No accepted history changes on refusal, budget exhaustion or cancellation.
    #[allow(clippy::too_many_arguments)]
    pub fn step_nonlinear_candidate<'t, 'm>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &NonlinearWindingCircuit<'m>,
        external_inputs: &[f64],
        controls: CircuitEnthalpySolve,
        winding_controls: NonlinearWindingControls,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<NonlinearWindingEnthalpyCandidate<'t, 'm>, Error> {
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
                "matching nonlinear field/circuit/thermal clocks and complete cell bindings required",
            ));
        }
        let network = MaterialNetwork {
            cells: self.cells,
            resistors: &[],
            links: self.links,
        };
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        let mut conductors = BTreeMap::new();
        let mut core_regions = BTreeMap::new();
        for source in circuit.field().conductors() {
            check(&cancel)?;
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
            check(&cancel)?;
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
            &cancel,
        )?;
        let layout = self.magnetic_layout();
        let accepted = prepared
            .resolve(thermal, || (*cancel.borrow_mut())())
            .map_err(material)?;
        let (initial, _) = layout.nonlinear_materials_at(
            thermal,
            circuit.field().problem(),
            &conductors,
            &accepted.source,
            &routes.magnetic,
            &cancel,
        )?;
        if initial.len() != conductors.len()
            || initial.iter().any(|(triangle, source)| {
                conductors.get(triangle).is_none_or(|original| {
                    source.source_basis() != original.source_basis()
                        || source.resistivity_ohm_m().to_bits()
                            != original.resistivity_ohm_m().to_bits()
                        || !same_material(source.material(), original.material())
                })
            })
        {
            return Err(invalid(
                "accepted nonlinear winding/core sources must match the actual thermal endpoint",
            ));
        }
        let (thermal, circuit, report) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            initial_energy_increments_j,
            |state, _remaining| {
                check(&cancel)?;
                let coefficients = prepared
                    .resolve(state, || (*cancel.borrow_mut())())
                    .map_err(material)?;
                let (sources, magnetic) = layout.nonlinear_materials_at(
                    state,
                    circuit.field().problem(),
                    &conductors,
                    &coefficients.source,
                    &routes.magnetic,
                    &cancel,
                )?;
                let mut core = Vec::new();
                let mut winding = Vec::new();
                for (triangle, conductor) in sources {
                    check(&cancel)?;
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
                let proposal = circuit
                    .step_to_candidate(
                        target_s,
                        core,
                        winding,
                        external_inputs,
                        winding_controls,
                        &mut || (*cancel.borrow_mut())(),
                    )
                    .map_err(electrical)?;
                Ok(DrivenCoefficients {
                    circuit: proposal,
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    source: WindingThermalSample {
                        network: coefficients.source,
                        magnetic,
                    },
                    circuit_trials: 1,
                })
            },
            |proposal: &NonlinearWindingCandidate<'_>, dt, checkpoint| {
                let step = proposal.step();
                routes.deposit_losses(
                    &step.cells,
                    &step.field.cells,
                    &step.graph_resistors,
                    dt,
                    checkpoint,
                )
            },
            || (*cancel.borrow_mut())(),
        )?;
        check(&cancel)?;
        Ok(NonlinearWindingEnthalpyCandidate {
            thermal,
            circuit,
            report,
        })
    }
}
