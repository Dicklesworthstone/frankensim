//! Native three-dimensional winding terminals, saturation and finite heat.
//!
//! One original A–phi/circuit candidate is evaluated at each proposed thermal
//! endpoint. Source resistivity, eddy conductivity and optional magnetic entropy
//! use the same material cards and caloric histories. The existing enthalpy
//! owner deposits each physical loss once; numerical magnetic loss stays out.

use super::*;
use crate::circuit_enthalpy::{
    ResistorCell,
    contacts::heat::{HeatAssembly, HeatNetworkError, HeatNetworkSample, closed_balance},
    magnetic::MagneticRoutes,
};
use fs_em::{
    spatial::vector_potential::induction::nonlinear::winding::{
        SpatialNonlinearWindingCandidate, SpatialNonlinearWindingCircuit,
        SpatialNonlinearWindingError, SpatialNonlinearWindingStep, SpatialWindingConductor,
    },
    winding::nonlinear::NonlinearWindingControls,
};
use fs_exec::Cx;
use fs_phs::circuit::Branch;
use std::collections::BTreeSet;

#[derive(Debug)]
pub enum SpatialWindingThermalError<E = HeatNetworkError> {
    Native(MagneticThermalError),
    Winding(SpatialNonlinearWindingError),
    Heat(E),
}
impl<E> From<MagneticThermalError> for SpatialWindingThermalError<E> {
    fn from(error: MagneticThermalError) -> Self {
        Self::Native(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for SpatialWindingThermalError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spatial winding/finite heat: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for SpatialWindingThermalError<E> {}
type WindingError<E = HeatNetworkError> = CircuitEnthalpyError<SpatialWindingThermalError<E>>;

fn electrical<E>(error: SpatialNonlinearWindingError) -> WindingError<E> {
    match error {
        SpatialNonlinearWindingError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(SpatialWindingThermalError::Winding(other)),
    }
}
fn heat_error<E: From<HeatNetworkError>>(error: HeatNetworkError) -> WindingError<E> {
    match error {
        HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(SpatialWindingThermalError::Heat(other.into())),
    }
}

/// Whole-volume caloric ownership for stationary distributed winding cells.
/// Multiple insulated strands may share one tetrahedron only when they use
/// the same original material and their declared fill fractions sum to one.
/// Partial-fill composites require a mixture caloric owner and refuse here.
pub struct SpatialWindingThermalNetwork<'a> {
    pub magnetic: SpatialMagneticThermalNetwork<'a>,
    /// Every authored numeric graph resistor needs exactly one heat receiver.
    /// Its resistance stays authored; source-resolved R(T) acts in the spatial
    /// winding/conductor cells, not by silently assigning graph material.
    pub graph_resistors: &'a [ResistorCell],
    /// Absolute, dimensionless allowance for summation roundoff, at most 1e-10.
    /// No fill fraction or source volume is rescaled to meet this allowance.
    pub fill_tolerance: f64,
}

/// Joint accepted circuit, full A–phi history and original finite caloric state.
pub struct SpatialWindingEnthalpyCandidate<
    't,
    'm,
    S = MaterialNetworkSample,
    C: SpatialConductorMaterial = ResolvedConductorState,
> {
    thermal: DepositionSimulation<'t>,
    circuit: SpatialNonlinearWindingCircuit<'m, C>,
    report: SpatialInductionThermalReport<SpatialNonlinearWindingStep, S>,
}
impl<'t, 'm, S, C: SpatialConductorMaterial> SpatialWindingEnthalpyCandidate<'t, 'm, S, C> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn circuit(&self) -> &SpatialNonlinearWindingCircuit<'m, C> {
        &self.circuit
    }
    pub fn report(&self) -> &SpatialInductionThermalReport<SpatialNonlinearWindingStep, S> {
        &self.report
    }
    /// Publish all histories together. Physical heat and reversible field heat
    /// have already entered the caloric stores and must not be deposited again.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        SpatialNonlinearWindingCircuit<'m, C>,
        SpatialInductionThermalReport<SpatialNonlinearWindingStep, S>,
    ) {
        (self.thermal, self.circuit, self.report)
    }
}

impl<'m, C: SpatialConductorMaterial> ThermalDriver for SpatialNonlinearWindingCandidate<'m, C> {
    type Model = SpatialNonlinearWindingCircuit<'m, C>;
    type Step = SpatialNonlinearWindingStep;
    fn account(&self) -> DriverAccount<'_> {
        Self::receipt_account(self.step())
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        SpatialNonlinearWindingCandidate::into_parts(self)
    }
}
impl<C: SpatialConductorMaterial> InductionThermalDriver
    for SpatialNonlinearWindingCandidate<'_, C>
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

struct Routes {
    winding: BTreeMap<usize, CellId>,
    graph: BTreeMap<usize, CellId>,
}

impl SpatialWindingThermalNetwork<'_> {
    /// Joint stationary 3D field/circuit/material/finite-enthalpy step.
    ///
    /// With `magnetic_laws=None`, the original H(B) law stays temperature
    /// independent. Explicit bindings admit the existing linear mu(T) or
    /// source-separable H(B,T) law and close its reversible entropy exchange.
    /// Original A, circuit flux/charge and thermal histories restart unchanged
    /// in every Newton trial. The work budget counts actual nested field trials.
    ///
    /// Eddy heat, stranded heat and graph heat each have their original receiver.
    /// This adapter is stationary: it does not omit heat advection from a moving
    /// solid. Fixed geometry, wholly solid constant-density caloric charts and
    /// the existing first-order coupled scheme remain the admitted scope.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate<'t, 'm, C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &SpatialNonlinearWindingCircuit<'m, C>,
        external_inputs: &[f64],
        magnetic_laws: Option<&[SpatialMagneticThermalLaw<MagneticConstitutiveThermalLaw>]>,
        controls: CircuitEnthalpySolve,
        winding_controls: NonlinearWindingControls<SpatialNonlinearInductionControls>,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialWindingEnthalpyCandidate<'t, 'm, MaterialNetworkSample, C>, WindingError>
    {
        let cancel = RefCell::new(cancelled);
        let network = MaterialNetwork {
            cells: self.magnetic.cells,
            resistors: &[],
            links: self.magnetic.links,
        };
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)
            .map_err(native_error)?;
        self.step_resolved_candidate(
            thermal,
            target_s,
            circuit,
            external_inputs,
            magnetic_laws,
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

    fn routes<C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'_>,
        circuit: &SpatialNonlinearWindingCircuit<'_, C>,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<Routes, Error> {
        check(cancel)?;
        if circuit.time_s().to_bits() != thermal.time_s().to_bits()
            || !self.fill_tolerance.is_finite()
            || self.fill_tolerance < 0.0
            || self.fill_tolerance > 1e-10
        {
            return Err(invalid(
                "matching circuit/thermal clocks and bounded fill roundoff required",
            ));
        }
        let problem = circuit.field().problem();
        let regions = self.magnetic.regions();
        let layout = MagneticThermalNetwork {
            cells: self.magnetic.cells,
            links: self.magnetic.links,
            magnetic_regions: &regions,
            volume_tolerance_m3: self.magnetic.volume_tolerance_m3,
        };
        let volumes = (0..problem.mesh().tets.len())
            .map(|cell| {
                check(cancel)?;
                problem
                    .cell_volume_m3(cell)
                    .ok_or_else(|| invalid("winding tetrahedron lacks original physical volume"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let mut fills = BTreeMap::<usize, f64>::new();
        let mut material_by_cell = BTreeMap::new();
        for source in circuit.conductors() {
            check(cancel)?;
            let fill = fills.entry(source.cell).or_default();
            *fill = finite(*fill + source.fill_fraction)?;
            let original = problem
                .materials()
                .get(source.cell)
                .ok_or_else(|| invalid("winding thermal source names a foreign tetrahedron"))?;
            if source.conductor.material().card_identity() != original.material().card_identity()
                || material_by_cell
                    .get(&source.cell)
                    .is_some_and(|previous| !same_material(source.conductor.material(), previous))
            {
                return Err(invalid(
                    "winding caloric volume requires one original material source",
                ));
            }
            material_by_cell.insert(source.cell, source.conductor.material().clone());
        }
        if fills
            .values()
            .any(|fill| (*fill - 1.0).abs() > self.fill_tolerance)
        {
            return Err(invalid(
                "partially filled winding cells require explicit mixture caloric volumes",
            ));
        }
        let MagneticRoutes {
            conductors: winding,
            ..
        } = layout.routes_for_volumes(
            thermal,
            &volumes,
            &fills.keys().copied().collect::<Vec<_>>(),
            cancel,
        )?;
        let mut graph = BTreeMap::new();
        for route in self.graph_resistors {
            check(cancel)?;
            if !thermal.cells().contains_key(&route.cell)
                || !matches!(
                    circuit.graph().branches.get(route.branch),
                    Some((_, _, Branch::Resistor { .. }))
                )
                || graph.insert(route.branch, route.cell).is_some()
            {
                return Err(invalid(
                    "graph resistor heat needs unique existing branch/cell correspondence",
                ));
            }
        }
        if circuit
            .graph()
            .branches
            .iter()
            .filter(|(_, _, branch)| matches!(branch, Branch::Resistor { .. }))
            .count()
            != graph.len()
        {
            return Err(invalid(
                "every native graph resistor requires a finite thermal receiver",
            ));
        }
        Ok(Routes { winding, graph })
    }

    #[allow(clippy::too_many_arguments)]
    fn step_resolved_candidate<'t, 'm, S, E, C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &SpatialNonlinearWindingCircuit<'m, C>,
        external_inputs: &[f64],
        magnetic_laws: Option<&[SpatialMagneticThermalLaw<MagneticConstitutiveThermalLaw>]>,
        controls: CircuitEnthalpySolve,
        winding_controls: NonlinearWindingControls<SpatialNonlinearInductionControls>,
        initial_energy_increments_j: &[f64],
        resolve: impl Fn(
            &DepositionSimulation<'_>,
        ) -> Result<CircuitThermalCoefficients<S>, WindingError<E>>,
        sample: impl Fn(&S) -> &MaterialNetworkSample,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<SpatialWindingEnthalpyCandidate<'t, 'm, S, C>, WindingError<E>> {
        let routes = self
            .routes(thermal, circuit, cancel)
            .map_err(native_error)?;
        let refresh = |sources: &MaterialNetworkSample| {
            circuit
                .conductors()
                .iter()
                .map(|original| {
                    check(cancel)?;
                    let receiver = routes.winding[&original.cell];
                    let binding = self
                        .magnetic
                        .cells
                        .iter()
                        .find(|b| b.cell == receiver)
                        .ok_or_else(|| {
                            invalid("winding receiver lacks its original material binding")
                        })?;
                    let phase = sources
                        .cells
                        .iter()
                        .find(|c| c.cell == receiver)
                        .ok_or_else(|| {
                            invalid("winding receiver lacks its resolved endpoint phase")
                        })?
                        .phase;
                    let query = binding.point.at(phase).map_err(|e| {
                        material(
                            crate::circuit_enthalpy::material::MaterialNetworkError::Material(e),
                        )
                    })?;
                    let refreshed = original
                        .conductor
                        .material()
                        .at_query_point(binding.point.material, &query)
                        .map_err(|e| {
                            material(
                                crate::circuit_enthalpy::material::MaterialNetworkError::Material(
                                    fs_material::thermal::source::ThermalMaterialError::Material(e),
                                ),
                            )
                        })?;
                    let conductor = original.conductor.with_material(refreshed).map_err(|e| {
                        CircuitEnthalpyError::Source(MagneticThermalError::Conductor(e))
                    })?;
                    Ok(SpatialWindingConductor {
                        conductor,
                        ..original.clone()
                    })
                })
                .collect::<Result<Vec<_>, Error>>()
        };
        let initial = resolve(thermal)?;
        let accepted = refresh(sample(&initial.source)).map_err(native_error)?;
        for (original, refreshed) in circuit.conductors().iter().zip(&accepted) {
            check(cancel).map_err(native_error)?;
            if !same_material(
                refreshed.conductor.material(),
                original.conductor.material(),
            ) || refreshed.conductor.resistivity_ohm_m().to_bits()
                != original.conductor.resistivity_ohm_m().to_bits()
            {
                return Err(native_error(invalid(
                    "accepted winding source differs from its actual caloric state",
                )));
            }
        }
        let velocity = vec![[0.0; 3]; circuit.field().problem().mesh().vertex_count];
        let (thermal, circuit, report) = self.magnetic.step_constitutive_driver_candidate(
            thermal,
            target_s,
            circuit.field(),
            magnetic_laws,
            controls,
            initial_energy_increments_j,
            resolve,
            sample,
            |_, samples, materials, conductors, remaining| {
                let sources = refresh(samples).map_err(native_error)?;
                let mut inner = winding_controls;
                inner.maximum_field_trials = inner.maximum_field_trials.min(remaining);
                if inner.maximum_field_trials == 0 {
                    return Err(CircuitEnthalpyError::CircuitBudget);
                }
                let candidate = if let Some(materials) = materials {
                    circuit.step_with_materials_candidate(
                        target_s,
                        &velocity,
                        materials,
                        conductors,
                        sources,
                        external_inputs,
                        inner,
                        &mut || (*cancel.borrow_mut())(),
                    )
                } else {
                    circuit.step_with_conductors_candidate(
                        target_s,
                        &velocity,
                        conductors,
                        sources,
                        external_inputs,
                        inner,
                        &mut || (*cancel.borrow_mut())(),
                    )
                };
                candidate.map_err(|error| match error {
                    SpatialNonlinearWindingError::TrialBudget
                        if remaining < winding_controls.maximum_field_trials =>
                    {
                        CircuitEnthalpyError::CircuitBudget
                    }
                    other => electrical(other),
                })
            },
            |candidate, powers, dt, checkpoint| {
                let step = candidate.step();
                let mut seen = BTreeSet::new();
                for loss in &step.cells {
                    checkpoint()?;
                    let receiver = *routes.winding.get(&loss.cell).ok_or_else(|| {
                        native_error(invalid("spatial strand heat lacks its original receiver"))
                    })?;
                    if !seen.insert((loss.winding, loss.cell)) {
                        return Err(native_error(invalid(
                            "repeated spatial strand heat receipt",
                        )));
                    }
                    add_heat(powers, receiver, loss.heat_j, dt).map_err(native_error)?;
                }
                if seen.len() != circuit.conductors().len()
                    || step.graph_resistors.len() != routes.graph.len()
                {
                    return Err(native_error(invalid(
                        "native winding/graph heat receipt is incomplete",
                    )));
                }
                let mut seen_graph = BTreeSet::new();
                for loss in &step.graph_resistors {
                    checkpoint()?;
                    let receiver = *routes.graph.get(&loss.branch).ok_or_else(|| {
                        native_error(invalid("native graph heat lacks its original receiver"))
                    })?;
                    if !seen_graph.insert(loss.branch) {
                        return Err(native_error(invalid("repeated native graph heat receipt")));
                    }
                    add_heat(powers, receiver, loss.heat_j, dt).map_err(native_error)?;
                }
                Ok(())
            },
            cancel,
        )?;
        Ok(SpatialWindingEnthalpyCandidate {
            thermal,
            circuit,
            report,
        })
    }
}

fn add_heat(
    powers: &mut BTreeMap<CellId, f64>,
    cell: CellId,
    heat_j: f64,
    dt: f64,
) -> Result<(), Error> {
    let power = finite(heat_j / dt)?;
    if heat_j < 0.0 || heat_j != 0.0 && power == 0.0 {
        return Err(invalid(
            "spatial winding heat is negative or lost in power conversion",
        ));
    }
    let old = powers.get(&cell).copied().unwrap_or(0.0);
    powers.insert(cell, finite(old + power)?);
    Ok(())
}

/// Continue the same native field/circuit/caloric history through a complete
/// original heat assembly, including thermal contacts, radiation and coolant.
/// Thermal boundaries remain field independent; numeric graph resistor heat
/// routes are explicit, while spatial conductor sources follow temperature.
#[allow(clippy::too_many_arguments)]
pub fn heat_assembly_candidate<'t, 'm, A: HeatAssembly + ?Sized, C: SpatialConductorMaterial>(
    provider: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    circuit: &SpatialNonlinearWindingCircuit<'m, C>,
    external_inputs: &[f64],
    regions: &[super::super::SpatialMagneticRegion],
    graph_resistors: &[ResistorCell],
    volume_tolerance_m3: f64,
    fill_tolerance: f64,
    magnetic_laws: Option<&[SpatialMagneticThermalLaw<MagneticConstitutiveThermalLaw>]>,
    controls: CircuitEnthalpySolve,
    winding_controls: NonlinearWindingControls<SpatialNonlinearInductionControls>,
    initial_energy_increments_j: &[f64],
) -> Result<SpatialWindingEnthalpyCandidate<'t, 'm, A::Sample, C>, WindingError<A::Error>> {
    let cancel = RefCell::new(|| cx.checkpoint().is_err());
    let prepared = provider.prepared_heat_network();
    let bulk = prepared.thermal_material_network().map_err(heat_error)?;
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
    let mut next = layout.step_resolved_candidate(
        thermal,
        target_s,
        circuit,
        external_inputs,
        magnetic_laws,
        controls,
        winding_controls,
        initial_energy_increments_j,
        |state| {
            let coefficients = provider
                .resolve_heat(cx, state)
                .map_err(|e| e.map_source(SpatialWindingThermalError::Heat))?;
            if !coefficients.contacts.is_empty() || !coefficients.bulk.resistors.is_empty() {
                return Err(native_error(invalid(
                    "spatial winding heat cannot discard sourced graph/contact electrical work",
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
