//! Actual voltage-driven nonlinear fields and their original finite thermal cells.
//!
//! The native circuit and nonlinear A_z histories are solved afresh in every
//! thermal trial. The original source cards supply resistivity at each proposed
//! endpoint; the same magnetic laws, geometry and caloric charts remain pinned.
//! Physical winding/core/graph Joule heat is deposited once. Numerical field
//! dissipation and private circuit-port work are never mistaken for heat.

/// Reversible magnetic temperature feedback on this same circuit/field solve.
pub mod thermal;

use super::{
    CircuitEnthalpyError, CircuitEnthalpyReport, CircuitEnthalpySolve, DrivenCoefficients,
    DriverAccount, Error, ThermalDriver, WindingThermalError, WindingThermalNetwork,
    WindingThermalSample, cell_driver_candidate, check, invalid, material, same_material,
};
use crate::circuit_enthalpy::magnetic::caloric::{
    MagneticEnergyCell, MagneticThermalLaw, PreparedMagneticCaloric,
};
use crate::circuit_enthalpy::material::{CellResistor, MaterialNetwork};
use fs_em::planar::nonlinear::IsotropicMagneticMaterial;
use fs_em::{
    planar::eddy::EddyConductor,
    winding::{
        nonlinear::{
            NonlinearWindingCandidate, NonlinearWindingCircuit, NonlinearWindingControls,
            NonlinearWindingError, NonlinearWindingStep, resistors::GraphResistorMaterial,
        },
        ohmic::WindingConductor,
    },
};
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_material::conductor::ElectricalPropertyBasis;
use fs_phs::circuit::Branch;
use fs_thermal::enthalpy::DepositionSimulation;
use std::{cell::RefCell, collections::BTreeMap};
use thermal::MagneticCaloricBalance;

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
    /// Independent physical internal-energy balance for a magnetic thermal
    /// step. The ordinary report also retains the native free-energy/parameter
    /// work account; only this balance treats reversible field heat internally.
    pub fn caloric_balance(&self) -> Option<&MagneticCaloricBalance> {
        self.report.source.magnetic_balance.as_ref()
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
            // Magnetic parameter work is explicit at this intermediate
            // free-energy boundary; the caloric owner closes its destination.
            supplied_j: step.supplied_work_j + step.parameter_work_j,
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
        if !circuit.graph_resistor_materials().is_empty() {
            return Err(invalid(
                "source-bound graph resistors require their complete thermal material bindings",
            ));
        }
        self.step_nonlinear_bound(
            thermal,
            target_s,
            circuit,
            external_inputs,
            None,
            None,
            controls,
            winding_controls,
            initial_energy_increments_j,
            cancelled,
        )
    }

    /// The same finite circuit/field/enthalpy residual with every ordinary graph
    /// resistor resolved from its OWN whole-cell material and endpoint state.
    /// The accepted native R and any retained source snapshot must match the
    /// initial thermal state. Original claim choices are pinned throughout the
    /// interval; conductivity and resistivity inputs share the existing owner.
    ///
    /// Bindings must match the existing graph heat destinations exactly. These
    /// separate bulk-resistor bodies cannot also be volumes covered by magnetic
    /// triangles. Their source, geometry, phase and material-conduction links
    /// are admitted by PreparedMaterialNetwork; no second heat or R(T) law is
    /// introduced. Source-domain and late thermal refusals preserve all histories.
    #[allow(clippy::too_many_arguments)]
    pub fn step_nonlinear_material_candidate<'t, 'm>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &NonlinearWindingCircuit<'m>,
        external_inputs: &[f64],
        resistors: &[CellResistor],
        electrical_bases: &[(usize, ElectricalPropertyBasis)],
        controls: CircuitEnthalpySolve,
        winding_controls: NonlinearWindingControls,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<NonlinearWindingEnthalpyCandidate<'t, 'm>, Error> {
        self.step_nonlinear_bound(
            thermal,
            target_s,
            circuit,
            external_inputs,
            Some((resistors, electrical_bases)),
            None,
            controls,
            winding_controls,
            initial_energy_increments_j,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_nonlinear_bound<'t, 'm>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &NonlinearWindingCircuit<'m>,
        external_inputs: &[f64],
        graph_material: Option<(&[CellResistor], &[(usize, ElectricalPropertyBasis)])>,
        magnetic_thermal: Option<(&[MagneticThermalLaw], WeightedAssemblyLimits)>,
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
        if graph_material.is_some()
            && (network.resistors.len() != routes.graph.len()
                || network.resistors.iter().any(|source| {
                    routes.graph.get(&source.destination.branch) != Some(&source.destination.cell)
                        || routes
                            .magnetic
                            .magnetic
                            .values()
                            .any(|cell| *cell == Some(source.destination.cell))
                }))
        {
            return Err(invalid(
                "sourced graph resistors need exact existing heat routes and distinct whole-cell volumes",
            ));
        }
        let layout = self.magnetic_layout();
        let caloric = if let Some((bindings, _)) = magnetic_thermal {
            if circuit
                .field()
                .problem()
                .materials()
                .iter()
                .any(|m| !matches!(m, IsotropicMagneticMaterial::Linear(_)))
            {
                return Err(invalid(
                    "magnetic caloric winding feedback currently requires sourced linear permeability in every region",
                ));
            }
            let field = circuit
                .field()
                .field(&mut || (*cancel.borrow_mut())())
                .map_err(|e| match e {
                    fs_em::planar::eddy::nonlinear::NonlinearEddyError::Cancelled => {
                        CircuitEnthalpyError::Cancelled
                    }
                    other => {
                        CircuitEnthalpyError::Source(WindingThermalError::NonlinearInduction(other))
                    }
                })?;
            let originals: Vec<_> = circuit
                .field()
                .problem()
                .materials()
                .iter()
                .map(IsotropicMagneticMaterial::material)
                .collect();
            Some(PreparedMagneticCaloric::new(
                &layout,
                thermal,
                circuit.field().problem().mesh(),
                &originals,
                field.magnetic_flux_density_t(),
                bindings,
                &routes.magnetic,
                &cancel,
            )?)
        } else {
            None
        };
        let accepted = prepared
            .resolve(thermal, || (*cancel.borrow_mut())())
            .map_err(material)?;
        if graph_material.is_some() {
            for source in &accepted.resistors {
                check(&cancel)?;
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
        let resolve_materials =
            |state: &DepositionSimulation<'_>,
             sample: &crate::circuit_enthalpy::material::MaterialNetworkSample| {
                if caloric.is_some() {
                    thermal::materials_at(
                        &layout,
                        state,
                        circuit.field().problem(),
                        &conductors,
                        sample,
                        &routes.magnetic,
                        &cancel,
                    )
                } else {
                    layout.nonlinear_materials_at(
                        state,
                        circuit.field().problem(),
                        &conductors,
                        sample,
                        &routes.magnetic,
                        &cancel,
                    )
                }
            };
        let (initial, initial_magnetic) = resolve_materials(thermal, &accepted.source)?;
        if let Some(caloric) = &caloric {
            caloric.validate_sources(thermal, &initial_magnetic, &cancel)?;
        }
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
        let (thermal, circuit, mut report) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            initial_energy_increments_j,
            |state, _remaining| {
                check(&cancel)?;
                let mut coefficients = prepared
                    .resolve(state, || (*cancel.borrow_mut())())
                    .map_err(material)?;
                let (sources, magnetic) = resolve_materials(state, &coefficients.source)?;
                if let Some(caloric) = &caloric {
                    caloric.validate_sources(state, &magnetic, &cancel)?;
                }
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
                let driven = if graph_material.is_some() {
                    circuit
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
                        .map_err(electrical)?
                } else {
                    circuit.clone()
                };
                let proposal = if let Some((_, assembly)) = magnetic_thermal {
                    driven
                        .step_with_materials_candidate(
                            target_s,
                            magnetic
                                .iter()
                                .map(|s| IsotropicMagneticMaterial::Linear(s.material.clone()))
                                .collect(),
                            core,
                            winding,
                            external_inputs,
                            assembly,
                            winding_controls,
                            &mut || (*cancel.borrow_mut())(),
                        )
                        .map_err(electrical)?
                } else {
                    driven
                        .step_to_candidate(
                            target_s,
                            core,
                            winding,
                            external_inputs,
                            winding_controls,
                            &mut || (*cancel.borrow_mut())(),
                        )
                        .map_err(electrical)?
                };
                let magnetic_caloric = if let Some(caloric) = &caloric {
                    Some(caloric.deposit(
                        state,
                        proposal.step().field.magnetic_cells.iter().map(|cell| {
                            MagneticEnergyCell {
                                triangle: cell.triangle,
                                flux_density_t: cell.flux_density_t,
                                free_energy_j: cell.free_energy_j,
                                parameter_work_j: cell.parameter_work_j,
                            }
                        }),
                        &mut coefficients.sources_w,
                        target_s - thermal.time_s(),
                        controls.total_energy_tolerance_j,
                        &cancel,
                    )?)
                } else {
                    None
                };
                Ok(DrivenCoefficients {
                    circuit: proposal,
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    source: WindingThermalSample {
                        network: coefficients.source,
                        magnetic,
                        magnetic_caloric,
                        magnetic_balance: None,
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
        report.source.magnetic_balance =
            thermal::balance(&report, controls.total_energy_tolerance_j)?;
        check(&cancel)?;
        Ok(NonlinearWindingEnthalpyCandidate {
            thermal,
            circuit,
            report,
        })
    }
}
