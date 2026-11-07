//! Native nonlinear electrode circuits on the original finite heat residual.
//!
//! Each thermal proposal restarts the accepted circuit and physical A–phi
//! electrode owner. Core conductivity and optional magnetic entropy follow the
//! actual endpoint temperatures; authored graph resistors may retain numerical
//! resistance or refresh their original material R(T) on separate finite bodies.
//! Only physical field and authored graph heat reach the thermal state.

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
use super::super::*;
use crate::circuit_enthalpy::{ResistorCell, contacts::heat::HeatAssembly, material::CellResistor};
use fs_em::spatial::vector_potential::induction::nonlinear::terminals::circuit::{
    SpatialNonlinearTerminalCircuit, SpatialNonlinearTerminalCircuitCandidate,
    SpatialNonlinearTerminalCircuitControls, SpatialNonlinearTerminalCircuitError,
    SpatialNonlinearTerminalCircuitStep,
};
use fs_exec::Cx;
use fs_material::conductor::ElectricalPropertyBasis;

impl<'m, C: SpatialConductorMaterial> ThermalDriver
    for SpatialNonlinearTerminalCircuitCandidate<'m, C>
{
    type Model = SpatialNonlinearTerminalCircuit<'m, C>;
    type Step = SpatialNonlinearTerminalCircuitStep;

    fn account(&self) -> DriverAccount<'_> {
        Self::receipt_account(self.step())
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        SpatialNonlinearTerminalCircuitCandidate::into_parts(self)
    }
}

impl<C: SpatialConductorMaterial> InductionThermalDriver
    for SpatialNonlinearTerminalCircuitCandidate<'_, C>
{
    fn induction_step(&self) -> &SpatialNonlinearInductionStep {
        &self.step().field
    }
    fn induction_receipt(step: &Self::Step) -> &SpatialNonlinearInductionStep {
        &step.field
    }
    fn receipt_account(step: &Self::Step) -> DriverAccount<'_> {
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: step.heat_j,
            numerical_dissipation_j: step.numerical_dissipation_j,
            stored_change_j: step.stored_change_j - step.parameter_work_j,
            supplied_j: step.supplied_work_j - step.mechanical_work_j,
        }
    }
    fn field_trials(&self) -> usize {
        self.step().field_trials
    }
}

impl SpatialTerminalCircuitThermalNetwork<'_> {
    /// Advance a physical electrode circuit and its finite caloric state as one
    /// candidate. Native circuit storage, electrode work and the spatial field
    /// are solved together at every thermal trial. Numerical reference-resistor
    /// work cancels inside the native circuit and is never deposited as heat.
    ///
    /// Core conductivity responds to actual temperature. Explicit magnetic laws
    /// enable the original linear mu(T) or separable H(B,T) entropy exchange;
    /// otherwise magnetic sources must remain fixed over the caloric interval.
    /// This stationary adapter retains original solid, whole-volume receivers.
    #[allow(clippy::too_many_arguments)]
    pub fn step_nonlinear_candidate<'t, 'm, C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &SpatialNonlinearTerminalCircuit<'m, C>,
        external_inputs: &[f64],
        endpoint_impressed_current_a_m2: &[[f64; 3]],
        magnetic_laws: Option<&[SpatialMagneticThermalLaw<MagneticConstitutiveThermalLaw>]>,
        controls: CircuitEnthalpySolve,
        circuit_controls: SpatialNonlinearTerminalCircuitControls,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        SpatialTerminalCircuitThermalCandidate<
            't,
            SpatialNonlinearTerminalCircuit<'m, C>,
            SpatialNonlinearTerminalCircuitStep,
        >,
        TerminalError,
    > {
        self.step_nonlinear_bound(
            thermal,
            target_s,
            circuit,
            external_inputs,
            endpoint_impressed_current_a_m2,
            None,
            magnetic_laws,
            controls,
            circuit_controls,
            initial_energy_increments_j,
            cancelled,
        )
    }

    /// Resolve all authored graph resistors from their own finite material
    /// bodies. Accepted source identity, geometry and R(T) must agree before
    /// stepping; every trial refreshes the same original sources without
    /// restarting circuit storage or electrode/field history.
    #[allow(clippy::too_many_arguments)]
    pub fn step_nonlinear_material_candidate<'t, 'm, C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &SpatialNonlinearTerminalCircuit<'m, C>,
        external_inputs: &[f64],
        endpoint_impressed_current_a_m2: &[[f64; 3]],
        resistors: &[CellResistor],
        electrical_bases: &[(usize, ElectricalPropertyBasis)],
        magnetic_laws: Option<&[SpatialMagneticThermalLaw<MagneticConstitutiveThermalLaw>]>,
        controls: CircuitEnthalpySolve,
        circuit_controls: SpatialNonlinearTerminalCircuitControls,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        SpatialTerminalCircuitThermalCandidate<
            't,
            SpatialNonlinearTerminalCircuit<'m, C>,
            SpatialNonlinearTerminalCircuitStep,
        >,
        TerminalError,
    > {
        self.step_nonlinear_bound(
            thermal,
            target_s,
            circuit,
            external_inputs,
            endpoint_impressed_current_a_m2,
            Some((resistors, electrical_bases)),
            magnetic_laws,
            controls,
            circuit_controls,
            initial_energy_increments_j,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_nonlinear_bound<'t, 'm, C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &SpatialNonlinearTerminalCircuit<'m, C>,
        external_inputs: &[f64],
        endpoint_impressed_current_a_m2: &[[f64; 3]],
        graph_material: Option<(&[CellResistor], &[(usize, ElectricalPropertyBasis)])>,
        magnetic_laws: Option<&[SpatialMagneticThermalLaw<MagneticConstitutiveThermalLaw>]>,
        controls: CircuitEnthalpySolve,
        circuit_controls: SpatialNonlinearTerminalCircuitControls,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        SpatialTerminalCircuitThermalCandidate<
            't,
            SpatialNonlinearTerminalCircuit<'m, C>,
            SpatialNonlinearTerminalCircuitStep,
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
        self.step_nonlinear_resolved_candidate(
            thermal,
            target_s,
            circuit,
            external_inputs,
            endpoint_impressed_current_a_m2,
            graph_material.map(|(resistors, _)| resistors),
            magnetic_laws,
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
    fn step_nonlinear_resolved_candidate<'t, 'm, S, E, C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &SpatialNonlinearTerminalCircuit<'m, C>,
        external_inputs: &[f64],
        endpoint_impressed_current_a_m2: &[[f64; 3]],
        graph_material: Option<&[CellResistor]>,
        magnetic_laws: Option<&[SpatialMagneticThermalLaw<MagneticConstitutiveThermalLaw>]>,
        controls: CircuitEnthalpySolve,
        circuit_controls: SpatialNonlinearTerminalCircuitControls,
        initial_energy_increments_j: &[f64],
        resolve: impl Fn(
            &DepositionSimulation<'_>,
        ) -> Result<CircuitThermalCoefficients<S>, TerminalError<E>>,
        sample: impl Fn(&S) -> &MaterialNetworkSample,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<
        SpatialTerminalCircuitThermalCandidate<
            't,
            SpatialNonlinearTerminalCircuit<'m, C>,
            SpatialNonlinearTerminalCircuitStep,
            S,
        >,
        TerminalError<E>,
    > {
        check(cancel).map_err(native_error)?;
        crate::circuit_enthalpy::admit_driver_unknowns(
            thermal,
            controls,
            initial_energy_increments_j,
        )?;
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
        let velocity = vec![[0.0; 3]; circuit.field().problem().positions_m().len()];
        let dt = target_s - thermal.time_s();
        let (thermal, circuit, report) = self.magnetic.step_constitutive_driver_candidate(
            thermal,
            target_s,
            circuit.field(),
            magnetic_laws,
            controls,
            initial_energy_increments_j,
            resolve,
            sample,
            |_, coefficients, materials, conductors, remaining| {
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
                let proposal = if let Some(materials) = materials {
                    driven.step_with_materials_candidate(
                        external_inputs,
                        endpoint_impressed_current_a_m2,
                        &velocity,
                        dt,
                        materials,
                        conductors,
                        inner,
                        &mut || (*cancel.borrow_mut())(),
                    )
                } else {
                    driven.step_with_conductors_candidate(
                        external_inputs,
                        endpoint_impressed_current_a_m2,
                        &velocity,
                        dt,
                        conductors,
                        inner,
                        &mut || (*cancel.borrow_mut())(),
                    )
                };
                proposal.map_err(|error| match error {
                    SpatialNonlinearTerminalCircuitError::FieldTrialBudget
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

/// Couple the same physical electrode circuit to the complete original heat
/// assembly, including thermal interfaces, radiation and finite material coolant.
/// Graph material bodies resolve their original R(T); every native graph loss
/// retains its explicit destination. Electrical contacts need a separate native
/// work owner and cannot be silently discarded. Boundaries are field independent.
#[allow(clippy::too_many_arguments)]
pub fn heat_assembly_candidate<'t, 'm, A: HeatAssembly + ?Sized, C: SpatialConductorMaterial>(
    provider: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    circuit: &SpatialNonlinearTerminalCircuit<'m, C>,
    external_inputs: &[f64],
    endpoint_impressed_current_a_m2: &[[f64; 3]],
    regions: &[SpatialMagneticRegion],
    graph_resistors: &[ResistorCell],
    volume_tolerance_m3: f64,
    magnetic_laws: Option<&[SpatialMagneticThermalLaw<MagneticConstitutiveThermalLaw>]>,
    controls: CircuitEnthalpySolve,
    circuit_controls: SpatialNonlinearTerminalCircuitControls,
    initial_energy_increments_j: &[f64],
) -> Result<
    SpatialTerminalCircuitThermalCandidate<
        't,
        SpatialNonlinearTerminalCircuit<'m, C>,
        SpatialNonlinearTerminalCircuitStep,
        A::Sample,
    >,
    TerminalError<A::Error>,
> {
    let cancel = RefCell::new(|| cx.checkpoint().is_err());
    check(&cancel).map_err(native_error)?;
    crate::circuit_enthalpy::admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
    let prepared = provider.prepared_heat_network();
    let bulk = prepared.field_material_network().map_err(heat_error)?;
    let mut magnetic_axes = Vec::new();
    for region in regions {
        check(&cancel).map_err(native_error)?;
        if let SpatialMagneticRegion::Cell { tetrahedron, cell } = *region {
            if let Some(SpatialNonlinearMaterial::Reversible(law)) =
                circuit.field().problem().materials().get(tetrahedron)
            {
                magnetic_axes.push((cell, law.source_path().axis.as_str()));
            }
        }
    }
    prepared
        .admit_field_independent_boundaries(cx, &magnetic_axes)
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
    let mut next = layout.step_nonlinear_resolved_candidate(
        thermal,
        target_s,
        circuit,
        external_inputs,
        endpoint_impressed_current_a_m2,
        (!bulk.resistors.is_empty()).then_some(bulk.resistors),
        magnetic_laws,
        controls,
        circuit_controls,
        initial_energy_increments_j,
        |state| {
            let coefficients = provider
                .resolve_heat(cx, state)
                .map_err(|error| error.map_source(SpatialTerminalCircuitThermalError::Heat))?;
            if !coefficients.contacts.is_empty() {
                return Err(native_error(invalid(
                    "electrode circuit heat cannot discard electrical contact work",
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
