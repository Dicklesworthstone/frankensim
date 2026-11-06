//! Native three-dimensional winding circuits, magnetic memory and finite heat.
//!
//! Every electrical/thermal trial restarts from the accepted circuit, A–phi,
//! vector play and caloric states. The original spatial winding owner supplies
//! terminal work, end storage and local strand losses; the field owner supplies
//! eddy and pinning heat. They enter the original enthalpy residual together.

pub use super::super::nonlinear::winding::{
    SpatialWindingThermalError, SpatialWindingThermalNetwork,
};
use super::super::{
    SpatialMagneticRegion,
    nonlinear::winding::{
        WindingError, admit_initial_graph_sources, deposit_winding_heat, electrical, heat_error,
        pinned_graph_sources,
    },
};
use super::{
    CircuitEnthalpyError, CircuitEnthalpySolve, CircuitThermalCoefficients, DriverAccount,
    HystereticInductionDriver, MaterialNetwork, MaterialNetworkSample, RefCell,
    SpatialHystereticInductionStep, SpatialInductionThermalReport, SpatialMagneticThermalNetwork,
    ThermalDriver, check, finite, invalid, material, native_error, same_material,
};
use crate::circuit_enthalpy::{
    ResistorCell,
    contacts::heat::{HeatAssembly, HeatNetworkSample, closed_balance},
    material::CellResistor,
};
use fs_em::spatial::vector_potential::induction::hysteretic::winding::{
    SpatialHystereticWindingCandidate, SpatialHystereticWindingCircuit,
    SpatialHystereticWindingControls, SpatialHystereticWindingError, SpatialHystereticWindingStep,
};
use fs_exec::Cx;
use fs_material::conductor::{
    ElectricalPropertyBasis, ResolvedConductorState, spatial::SpatialConductorMaterial,
};
use fs_thermal::enthalpy::DepositionSimulation;

/// One jointly accepted native circuit, vector memory and finite thermal state.
pub struct SpatialHystereticWindingEnthalpyCandidate<
    't,
    'm,
    S = MaterialNetworkSample,
    C: SpatialConductorMaterial = ResolvedConductorState,
> {
    thermal: DepositionSimulation<'t>,
    circuit: SpatialHystereticWindingCircuit<'m, C>,
    report: SpatialInductionThermalReport<SpatialHystereticWindingStep, S>,
}

impl<'t, 'm, S, C: SpatialConductorMaterial>
    SpatialHystereticWindingEnthalpyCandidate<'t, 'm, S, C>
{
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn circuit(&self) -> &SpatialHystereticWindingCircuit<'m, C> {
        &self.circuit
    }
    pub fn report(&self) -> &SpatialInductionThermalReport<SpatialHystereticWindingStep, S> {
        &self.report
    }
    /// Physical core, strand and graph heat already paid to finite receivers.
    pub fn deposited_heat_j(&self) -> f64 {
        self.report.coupled.deposited_joule_j
    }
    /// Publish all histories together. Physical heat is already deposited;
    /// numerical endpoint dissipation remains a separate energy export.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        SpatialHystereticWindingCircuit<'m, C>,
        SpatialInductionThermalReport<SpatialHystereticWindingStep, S>,
    ) {
        (self.thermal, self.circuit, self.report)
    }
}

impl<'m, C: SpatialConductorMaterial> ThermalDriver for SpatialHystereticWindingCandidate<'m, C> {
    type Model = SpatialHystereticWindingCircuit<'m, C>;
    type Step = SpatialHystereticWindingStep;

    fn account(&self) -> DriverAccount<'_> {
        Self::receipt_account(self.step())
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        SpatialHystereticWindingCandidate::into_parts(self)
    }
}

impl<C: SpatialConductorMaterial> HystereticInductionDriver
    for SpatialHystereticWindingCandidate<'_, C>
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

impl SpatialWindingThermalNetwork<'_> {
    /// Advance the native 3D winding/circuit, vector hysteresis and enthalpy as
    /// one candidate. Original strand resistivity and core conductivity follow
    /// the actual endpoint temperatures; magnetic fit coefficients remain fixed
    /// throughout each receiver's admitted caloric interval.
    ///
    /// Core eddy heat, irreversible pinning heat, strand copper heat and ordinary
    /// graph-resistor heat have explicit receivers. Numerical field/end-storage
    /// dissipation is not heat. Fixed geometry and fully solid whole-volume
    /// receivers are required, including electrically insulating lossy cores.
    /// Source changes and physical reversals need resolved time endpoints.
    #[allow(clippy::too_many_arguments)]
    pub fn step_hysteretic_candidate<'t, 'm, C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &SpatialHystereticWindingCircuit<'m, C>,
        external_inputs: &[f64],
        controls: CircuitEnthalpySolve,
        winding_controls: SpatialHystereticWindingControls,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        SpatialHystereticWindingEnthalpyCandidate<'t, 'm, MaterialNetworkSample, C>,
        WindingError,
    > {
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

    /// Resolve every ordinary graph resistor from its own finite material body.
    /// Original geometry/source identity and the accepted R(T) must agree;
    /// source refresh retains all accepted circuit and magnetic histories.
    #[allow(clippy::too_many_arguments)]
    pub fn step_hysteretic_material_candidate<'t, 'm, C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &SpatialHystereticWindingCircuit<'m, C>,
        external_inputs: &[f64],
        resistors: &[CellResistor],
        electrical_bases: &[(usize, ElectricalPropertyBasis)],
        controls: CircuitEnthalpySolve,
        winding_controls: SpatialHystereticWindingControls,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        SpatialHystereticWindingEnthalpyCandidate<'t, 'm, MaterialNetworkSample, C>,
        WindingError,
    > {
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
    fn step_hysteretic_bound<'t, 'm, C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &SpatialHystereticWindingCircuit<'m, C>,
        external_inputs: &[f64],
        graph_material: Option<(&[CellResistor], &[(usize, ElectricalPropertyBasis)])>,
        controls: CircuitEnthalpySolve,
        winding_controls: SpatialHystereticWindingControls,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        SpatialHystereticWindingEnthalpyCandidate<'t, 'm, MaterialNetworkSample, C>,
        WindingError,
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
            graph_material.map(|(resistors, _)| resistors),
            controls,
            winding_controls,
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
        circuit: &SpatialHystereticWindingCircuit<'m, C>,
        external_inputs: &[f64],
        graph_material: Option<&[CellResistor]>,
        controls: CircuitEnthalpySolve,
        winding_controls: SpatialHystereticWindingControls,
        initial_energy_increments_j: &[f64],
        resolve: impl Fn(
            &DepositionSimulation<'_>,
        ) -> Result<CircuitThermalCoefficients<S>, WindingError<E>>,
        sample: impl Fn(&S) -> &MaterialNetworkSample,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<SpatialHystereticWindingEnthalpyCandidate<'t, 'm, S, C>, WindingError<E>> {
        check(cancel).map_err(native_error)?;
        let problem = circuit.field().problem();
        let volumes = (0..problem.mesh().tets.len())
            .map(|cell| {
                check(cancel)?;
                problem
                    .cell_volume_m3(cell)
                    .ok_or_else(|| invalid("hysteretic winding tetrahedron lacks physical volume"))
            })
            .collect::<Result<Vec<_>, super::Error>>()
            .map_err(native_error)?;
        let originals = problem
            .materials()
            .iter()
            .map(|law| law.material())
            .collect::<Vec<_>>();
        let routes = self
            .routes_for_sources(
                thermal,
                circuit.time_s(),
                &volumes,
                &originals,
                circuit.conductors(),
                circuit.graph(),
                cancel,
            )
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
        let accepted = self
            .refresh_winding_sources(
                circuit.conductors(),
                &routes,
                sample(&initial.source),
                cancel,
            )
            .map_err(native_error)?;
        for (original, refreshed) in circuit.conductors().iter().zip(&accepted) {
            check(cancel).map_err(native_error)?;
            if !same_material(
                refreshed.conductor.material(),
                original.conductor.material(),
            ) || refreshed.conductor.resistivity_ohm_m().to_bits()
                != original.conductor.resistivity_ohm_m().to_bits()
            {
                return Err(native_error(invalid(
                    "accepted hysteretic winding source differs from its actual caloric state",
                )));
            }
        }
        let (thermal, circuit, report) = self.magnetic.step_hysteretic_driver_candidate(
            thermal,
            target_s,
            circuit.field(),
            controls,
            initial_energy_increments_j,
            resolve,
            &sample,
            |_, coefficients, conductors, remaining| {
                let sources = self
                    .refresh_winding_sources(
                        circuit.conductors(),
                        &routes,
                        sample(&coefficients.source),
                        cancel,
                    )
                    .map_err(native_error)?;
                let mut inner = winding_controls;
                inner.maximum_field_trials = inner.maximum_field_trials.min(remaining);
                if inner.maximum_field_trials == 0 {
                    return Err(CircuitEnthalpyError::CircuitBudget);
                }
                let refreshed_graph;
                let driven = if graph_material.is_some() {
                    let graph_sources =
                        pinned_graph_sources(&initial.resistors, &coefficients.resistors, cancel)
                            .map_err(native_error)?;
                    refreshed_graph = circuit
                        .with_graph_resistors(graph_sources, &mut || (*cancel.borrow_mut())())
                        .map_err(electrical)?;
                    &refreshed_graph
                } else {
                    if !coefficients.resistors.is_empty() {
                        return Err(native_error(invalid(
                            "unowned graph resistor source appeared in a hysteretic thermal trial",
                        )));
                    }
                    circuit
                };
                driven
                    .step_with_conductors_candidate(
                        target_s,
                        conductors,
                        sources,
                        external_inputs,
                        inner,
                        &mut || (*cancel.borrow_mut())(),
                    )
                    .map_err(|error| match error {
                        SpatialHystereticWindingError::TrialBudget
                            if remaining < winding_controls.maximum_field_trials =>
                        {
                            CircuitEnthalpyError::CircuitBudget
                        }
                        other => electrical(other),
                    })
            },
            |candidate, powers, dt, checkpoint| {
                let step = candidate.step();
                deposit_winding_heat(
                    &step.cells,
                    &step.graph_resistors,
                    circuit.conductors().len(),
                    &routes,
                    powers,
                    dt,
                    checkpoint,
                )
                .map_err(native_error)
            },
            cancel,
        )?;
        Ok(SpatialHystereticWindingEnthalpyCandidate {
            thermal,
            circuit,
            report,
        })
    }
}

/// Couple retained 3D winding hysteresis to the original complete heat assembly.
/// Contacts, radiation and coolant use their existing finite histories and source
/// payloads. Boundaries are field independent; electrical contact work requires
/// its own electrical owner and cannot be silently discarded by this adapter.
#[allow(clippy::too_many_arguments)]
pub fn heat_assembly_candidate<'t, 'm, A: HeatAssembly + ?Sized, C: SpatialConductorMaterial>(
    provider: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    circuit: &SpatialHystereticWindingCircuit<'m, C>,
    external_inputs: &[f64],
    regions: &[SpatialMagneticRegion],
    graph_resistors: &[ResistorCell],
    volume_tolerance_m3: f64,
    fill_tolerance: f64,
    controls: CircuitEnthalpySolve,
    winding_controls: SpatialHystereticWindingControls,
    initial_energy_increments_j: &[f64],
) -> Result<SpatialHystereticWindingEnthalpyCandidate<'t, 'm, A::Sample, C>, WindingError<A::Error>>
{
    let cancel = RefCell::new(|| cx.checkpoint().is_err());
    let prepared = provider.prepared_heat_network();
    let bulk = prepared.field_material_network().map_err(heat_error)?;
    prepared
        .admit_field_independent_boundaries(cx, &[])
        .map_err(heat_error)?;
    let layout = SpatialWindingThermalNetwork {
        magnetic: SpatialMagneticThermalNetwork {
            cells: bulk.cells,
            links: bulk.links,
            magnetic_regions: regions,
            volume_tolerance_m3,
        },
        graph_resistors,
        fill_tolerance,
    };
    let mut next = layout.step_hysteretic_resolved_candidate(
        thermal,
        target_s,
        circuit,
        external_inputs,
        (!bulk.resistors.is_empty()).then_some(bulk.resistors),
        controls,
        winding_controls,
        initial_energy_increments_j,
        |state| {
            let coefficients = provider
                .resolve_heat(cx, state)
                .map_err(|error| error.map_source(SpatialWindingThermalError::Heat))?;
            if !coefficients.contacts.is_empty() {
                return Err(native_error(invalid(
                    "hysteretic winding heat cannot discard electrical contact work",
                )));
            }
            Ok(coefficients.bulk)
        },
        |sample: &A::Sample| &sample.as_ref().material,
        &cancel,
    )?;
    let dt = target_s - thermal.time_s();
    let sample: &HeatNetworkSample = next.report.coupled.source.network.as_ref();
    closed_balance(
        dt,
        next.report.coupled.total_energy_residual_j,
        sample,
        controls.total_energy_tolerance_j,
    )
    .map_err(heat_error)?;
    closed_balance(
        dt,
        next.report.total_energy_residual_j,
        sample,
        controls.total_energy_tolerance_j,
    )
    .map_err(heat_error)?;
    let radiation_net = finite(dt * sample.radiation_net_power_w())?;
    next.report.external_heat_j = finite(next.report.external_heat_j - radiation_net)?;
    next.report.total_energy_residual_j =
        finite(next.report.total_energy_residual_j + radiation_net)?;
    check(&cancel).map_err(native_error)?;
    Ok(next)
}
