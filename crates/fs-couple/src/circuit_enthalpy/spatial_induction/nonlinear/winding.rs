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
    material::CellResistor,
};
use fs_em::{
    spatial::vector_potential::induction::nonlinear::winding::{
        SpatialNonlinearWindingCandidate, SpatialNonlinearWindingCircuit,
        SpatialNonlinearWindingError, SpatialNonlinearWindingStep, SpatialWindingCellHeat,
        SpatialWindingConductor,
    },
    winding::{
        nonlinear::{NonlinearWindingControls, resistors::GraphResistorMaterial},
        ohmic::WindingGraphResistorHeat,
    },
};
use fs_exec::Cx;
use fs_material::{
    conductor::ElectricalPropertyBasis, field_storage::admit_source_temperature_refresh,
    state_point::ResolvedMaterialStatePoint,
};
use fs_phs::circuit::{Branch, CircuitGraph};
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
pub(in crate::circuit_enthalpy::spatial_induction) type WindingError<E = HeatNetworkError> =
    CircuitEnthalpyError<SpatialWindingThermalError<E>>;

pub(in crate::circuit_enthalpy::spatial_induction) fn electrical<E>(
    error: SpatialNonlinearWindingError,
) -> WindingError<E> {
    match error {
        SpatialNonlinearWindingError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(SpatialWindingThermalError::Winding(other)),
    }
}
pub(in crate::circuit_enthalpy::spatial_induction) fn heat_error<E: From<HeatNetworkError>>(
    error: HeatNetworkError,
) -> WindingError<E> {
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
    /// Every original graph resistor needs exactly one heat receiver. Numeric
    /// resistance stays authored unless its complete material binding is given
    /// to `step_material_candidate` or the complete material heat assembly.
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

pub(in crate::circuit_enthalpy::spatial_induction) struct Routes {
    pub(in crate::circuit_enthalpy::spatial_induction) winding: BTreeMap<usize, CellId>,
    pub(in crate::circuit_enthalpy::spatial_induction) graph: BTreeMap<usize, CellId>,
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
        self.step_bound(
            thermal,
            target_s,
            circuit,
            external_inputs,
            None,
            magnetic_laws,
            controls,
            winding_controls,
            initial_energy_increments_j,
            cancelled,
        )
    }

    /// Resolve every ordinary graph resistor from its own finite material body.
    ///
    /// Original length, area, source choices and accepted resistance must agree
    /// before stepping. Graph bodies cannot also be whole magnetic volumes.
    /// Each endpoint uses its own R(T), retained alongside native branch heat;
    /// no circuit, magnetic or thermal history is restarted by source refresh.
    #[allow(clippy::too_many_arguments)]
    pub fn step_material_candidate<'t, 'm, C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &SpatialNonlinearWindingCircuit<'m, C>,
        external_inputs: &[f64],
        resistors: &[CellResistor],
        electrical_bases: &[(usize, ElectricalPropertyBasis)],
        magnetic_laws: Option<&[SpatialMagneticThermalLaw<MagneticConstitutiveThermalLaw>]>,
        controls: CircuitEnthalpySolve,
        winding_controls: NonlinearWindingControls<SpatialNonlinearInductionControls>,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialWindingEnthalpyCandidate<'t, 'm, MaterialNetworkSample, C>, WindingError>
    {
        self.step_bound(
            thermal,
            target_s,
            circuit,
            external_inputs,
            Some((resistors, electrical_bases)),
            magnetic_laws,
            controls,
            winding_controls,
            initial_energy_increments_j,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_bound<'t, 'm, C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &SpatialNonlinearWindingCircuit<'m, C>,
        external_inputs: &[f64],
        graph_material: Option<(&[CellResistor], &[(usize, ElectricalPropertyBasis)])>,
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
        self.step_resolved_candidate(
            thermal,
            target_s,
            circuit,
            external_inputs,
            graph_material.map(|(resistors, _)| resistors),
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
        let problem = circuit.field().problem();
        let volumes = (0..problem.mesh().tets.len())
            .map(|cell| {
                check(cancel)?;
                problem
                    .cell_volume_m3(cell)
                    .ok_or_else(|| invalid("winding tetrahedron lacks original physical volume"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let originals = problem
            .materials()
            .iter()
            .map(SpatialNonlinearMaterial::material)
            .collect::<Vec<_>>();
        self.routes_for_sources(
            thermal,
            circuit.time_s(),
            &volumes,
            &originals,
            circuit.conductors(),
            circuit.graph(),
            cancel,
        )
    }

    // Reversible and hysteretic circuit owners share this original physical
    // volume/source correspondence; neither supplies another thermal topology.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::circuit_enthalpy::spatial_induction) fn routes_for_sources(
        &self,
        thermal: &DepositionSimulation<'_>,
        accepted_time_s: f64,
        volumes: &[f64],
        originals: &[&ResolvedMaterialStatePoint],
        strands: &[SpatialWindingConductor],
        circuit_graph: &CircuitGraph,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<Routes, Error> {
        check(cancel)?;
        if accepted_time_s.to_bits() != thermal.time_s().to_bits()
            || !self.fill_tolerance.is_finite()
            || self.fill_tolerance < 0.0
            || self.fill_tolerance > 1e-10
        {
            return Err(invalid(
                "matching circuit/thermal clocks and bounded fill roundoff required",
            ));
        }
        if originals.len() != volumes.len() {
            return Err(invalid(
                "complete original winding geometry/material maps required",
            ));
        }
        let regions = self.magnetic.regions();
        let layout = MagneticThermalNetwork {
            cells: self.magnetic.cells,
            links: self.magnetic.links,
            magnetic_regions: &regions,
            volume_tolerance_m3: self.magnetic.volume_tolerance_m3,
        };
        let mut fills = BTreeMap::<usize, f64>::new();
        let mut material_by_cell = BTreeMap::<usize, &ResolvedMaterialStatePoint>::new();
        for source in strands {
            check(cancel)?;
            let fill = fills.entry(source.cell).or_default();
            *fill = finite(*fill + source.fill_fraction)?;
            let original = originals
                .get(source.cell)
                .ok_or_else(|| invalid("winding thermal source names a foreign tetrahedron"))?;
            if source.conductor.material().card_identity() != original.card_identity()
                || material_by_cell
                    .get(&source.cell)
                    .is_some_and(|previous| !same_material(source.conductor.material(), previous))
            {
                return Err(invalid(
                    "winding caloric volume requires one original material source",
                ));
            }
            material_by_cell.insert(source.cell, source.conductor.material());
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
            volumes,
            &fills.keys().copied().collect::<Vec<_>>(),
            cancel,
        )?;
        let mut graph = BTreeMap::new();
        for route in self.graph_resistors {
            check(cancel)?;
            if !thermal.cells().contains_key(&route.cell)
                || !matches!(
                    circuit_graph.branches.get(route.branch),
                    Some((_, _, Branch::Resistor { .. }))
                )
                || graph.insert(route.branch, route.cell).is_some()
            {
                return Err(invalid(
                    "graph resistor heat needs unique existing branch/cell correspondence",
                ));
            }
        }
        if circuit_graph
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

    pub(in crate::circuit_enthalpy::spatial_induction) fn admit_graph_material_bindings(
        &self,
        retained: &[GraphResistorMaterial],
        bindings: Option<&[CellResistor]>,
        routes: &Routes,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<(), Error> {
        check(cancel)?;
        if bindings.is_none() && !retained.is_empty() {
            return Err(invalid(
                "retained graph resistor sources require their complete thermal bindings",
            ));
        }
        if let Some(bindings) = bindings {
            if bindings.len() != routes.graph.len() {
                return Err(invalid(
                    "graph material bindings need exact heat routes and separate whole-cell volumes",
                ));
            }
            if bindings.is_empty() {
                return Ok(());
            }
            let mut magnetic_cells = BTreeSet::new();
            for region in self.magnetic.magnetic_regions {
                check(cancel)?;
                if let super::super::SpatialMagneticRegion::Cell { cell, .. } = *region {
                    magnetic_cells.insert(cell);
                }
            }
            for binding in bindings {
                check(cancel)?;
                if routes.graph.get(&binding.destination.branch) != Some(&binding.destination.cell)
                    || magnetic_cells.contains(&binding.destination.cell)
                {
                    return Err(invalid(
                        "graph material bindings need exact heat routes and separate whole-cell volumes",
                    ));
                }
            }
        }
        Ok(())
    }

    // Resolve each original stranded source at the same endpoint phase as the
    // field. Build the thermal lookups once, not once per strand/cell pair.
    pub(in crate::circuit_enthalpy::spatial_induction) fn refresh_winding_sources(
        &self,
        strands: &[SpatialWindingConductor],
        routes: &Routes,
        sample: &MaterialNetworkSample,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<Vec<SpatialWindingConductor>, Error> {
        check(cancel)?;
        if strands.is_empty() {
            return Ok(Vec::new());
        }
        let mut bindings = BTreeMap::new();
        for binding in self.magnetic.cells {
            check(cancel)?;
            bindings.entry(binding.cell).or_insert(binding);
        }
        let mut phases = BTreeMap::new();
        for source in &sample.cells {
            check(cancel)?;
            phases.entry(source.cell).or_insert(source.phase);
        }
        strands
            .iter()
            .map(|original| {
                check(cancel)?;
                let receiver = *routes
                    .winding
                    .get(&original.cell)
                    .ok_or_else(|| invalid("winding source lacks its original thermal receiver"))?;
                let binding = bindings.get(&receiver).ok_or_else(|| {
                    invalid("winding receiver lacks its original material binding")
                })?;
                let phase = *phases
                    .get(&receiver)
                    .ok_or_else(|| invalid("winding receiver lacks its resolved endpoint phase"))?;
                let query = binding.point.at(phase).map_err(|e| {
                    material(crate::circuit_enthalpy::material::MaterialNetworkError::Material(e))
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
            .collect()
    }

    #[allow(clippy::too_many_arguments)]
    fn step_resolved_candidate<'t, 'm, S, E, C: SpatialConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        circuit: &SpatialNonlinearWindingCircuit<'m, C>,
        external_inputs: &[f64],
        graph_material: Option<&[CellResistor]>,
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
        self.admit_graph_material_bindings(
            circuit.graph_resistor_materials(),
            graph_material,
            &routes,
            cancel,
        )
        .map_err(native_error)?;
        let refresh = |sources: &MaterialNetworkSample| {
            self.refresh_winding_sources(circuit.conductors(), &routes, sources, cancel)
        };
        let initial = resolve(thermal)?;
        admit_initial_graph_sources(
            circuit.graph(),
            circuit.graph_resistor_materials(),
            graph_material.map_or(0, |bindings| bindings.len()),
            &initial.resistors,
            cancel,
        )
        .map_err(native_error)?;
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
            &sample,
            |_, coefficients, materials, conductors, remaining| {
                let sources = refresh(sample(&coefficients.source)).map_err(native_error)?;
                let mut inner = winding_controls;
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
                            "unowned graph resistor source appeared in a thermal trial",
                        )));
                    }
                    circuit
                };
                let candidate = if let Some(materials) = materials {
                    driven.step_with_materials_candidate(
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
                    driven.step_with_conductors_candidate(
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
        Ok(SpatialWindingEnthalpyCandidate {
            thermal,
            circuit,
            report,
        })
    }
}

pub(in crate::circuit_enthalpy::spatial_induction) fn admit_initial_graph_sources(
    graph: &CircuitGraph,
    retained: &[GraphResistorMaterial],
    expected_count: usize,
    sources: &[crate::material_circuit::MaterialResistor],
    cancel: &RefCell<impl FnMut() -> bool>,
) -> Result<(), Error> {
    if sources.len() != expected_count {
        return Err(invalid("graph thermal source coverage changed"));
    }
    for source in sources {
        check(cancel)?;
        let Some((_, _, Branch::Resistor { ohms })) = graph.branches.get(source.branch) else {
            return Err(invalid("material graph source names a non-resistor branch"));
        };
        if ohms.to_bits() != source.conductor.resistance_ohm().to_bits()
            || !retained.is_empty()
                && !retained.iter().any(|old| {
                    old.branch == source.branch
                        && same_graph_geometry(&old.conductor, &source.conductor)
                        && same_material(old.conductor.material(), source.conductor.material())
                })
        {
            return Err(invalid(
                "accepted graph resistance/source does not match its original thermal state",
            ));
        }
    }
    Ok(())
}

// Both magnetic histories consume the exact same original native strand and
// graph loss receipts. Core Joule/pinning heat belongs to their field adapter.
#[allow(clippy::too_many_arguments)]
pub(in crate::circuit_enthalpy::spatial_induction) fn deposit_winding_heat(
    cells: &[SpatialWindingCellHeat],
    graph_losses: &[WindingGraphResistorHeat],
    expected_strands: usize,
    routes: &Routes,
    powers: &mut BTreeMap<CellId, f64>,
    dt: f64,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<(), Error> {
    let mut seen = BTreeSet::new();
    for loss in cells {
        checkpoint()?;
        let receiver = *routes
            .winding
            .get(&loss.cell)
            .ok_or_else(|| invalid("spatial strand heat lacks its original receiver"))?;
        if !seen.insert((loss.winding, loss.cell)) {
            return Err(invalid("repeated spatial strand heat receipt"));
        }
        add_heat(powers, receiver, loss.heat_j, dt)?;
    }
    if seen.len() != expected_strands || graph_losses.len() != routes.graph.len() {
        return Err(invalid("native winding/graph heat receipt is incomplete"));
    }
    let mut seen_graph = BTreeSet::new();
    for loss in graph_losses {
        checkpoint()?;
        let receiver = *routes
            .graph
            .get(&loss.branch)
            .ok_or_else(|| invalid("native graph heat lacks its original receiver"))?;
        if !seen_graph.insert(loss.branch) {
            return Err(invalid("repeated native graph heat receipt"));
        }
        add_heat(powers, receiver, loss.heat_j, dt)?;
    }
    Ok(())
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

pub(in crate::circuit_enthalpy::spatial_induction) fn same_graph_geometry(
    original: &fs_material::conductor::ResolvedConductor,
    proposed: &fs_material::conductor::ResolvedConductor,
) -> bool {
    original.source_basis() == proposed.source_basis()
        && original.source_key() == proposed.source_key()
        && original.length_m().to_bits() == proposed.length_m().to_bits()
        && original.area_m2().to_bits() == proposed.area_m2().to_bits()
}

// A heat provider may choose transport coefficients, but it cannot replace an
// accepted resistor's original material source while solving its temperature.
pub(in crate::circuit_enthalpy::spatial_induction) fn pinned_graph_sources(
    original: &[crate::material_circuit::MaterialResistor],
    proposed: &[crate::material_circuit::MaterialResistor],
    cancel: &RefCell<impl FnMut() -> bool>,
) -> Result<Vec<GraphResistorMaterial>, Error> {
    if original.len() != proposed.len() {
        return Err(invalid("graph resistor thermal source coverage changed"));
    }
    let mut originals = BTreeMap::new();
    for source in original {
        check(cancel)?;
        originals.entry(source.branch).or_insert(source);
    }
    let mut seen = BTreeSet::new();
    proposed
        .iter()
        .map(|source| {
            check(cancel)?;
            let old = originals
                .get(&source.branch)
                .ok_or_else(|| invalid("thermal trial changed the graph resistor branch"))?;
            if !seen.insert(source.branch)
                || !same_graph_geometry(&old.conductor, &source.conductor)
            {
                return Err(invalid(
                    "thermal trial changed graph resistor geometry or source choice",
                ));
            }
            admit_source_temperature_refresh(old.conductor.material(), source.conductor.material())
                .map_err(|e| CircuitEnthalpyError::Source(MagneticThermalError::Magnetic(e)))?;
            if old
                .conductor
                .with_material(source.conductor.material().clone())
                .map_err(|e| CircuitEnthalpyError::Source(MagneticThermalError::Conductor(e)))?
                != source.conductor
            {
                return Err(invalid(
                    "graph resistor left its original electrical source law",
                ));
            }
            Ok(GraphResistorMaterial {
                branch: source.branch,
                conductor: source.conductor.clone(),
            })
        })
        .collect()
}

/// Continue the same native field/circuit/caloric history through a complete
/// original heat assembly, including thermal contacts, radiation and coolant.
/// Thermal boundaries remain field independent. Complete graph material
/// bindings refresh their original R(T) on separate finite bodies; otherwise
/// graph resistance stays authored. Every heat route remains explicit.
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
    let mut next = layout.step_resolved_candidate(
        thermal,
        target_s,
        circuit,
        external_inputs,
        (!bulk.resistors.is_empty()).then_some(bulk.resistors),
        magnetic_laws,
        controls,
        winding_controls,
        initial_energy_increments_j,
        |state| {
            let coefficients = provider
                .resolve_heat(cx, state)
                .map_err(|e| e.map_source(SpatialWindingThermalError::Heat))?;
            if !coefficients.contacts.is_empty() {
                return Err(native_error(invalid(
                    "spatial winding heat cannot discard electrical contact work",
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
