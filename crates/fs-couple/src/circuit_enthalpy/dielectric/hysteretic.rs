//! Spatial dielectric switching and finite caloric storage on one timeline.
//!
//! The existing Gauss/electrode owner advances one endpoint or adaptively
//! refined affine path from accepted memory. Its fixed fitted coefficients and
//! prescribed electrical loading do not depend on temperature, so thermal
//! Newton trials reuse that SAME immutable candidate and physical heat. They
//! refresh source admission and native thermal transport, not polarization or
//! the electrical path. The
//! accepted temperature queries and both histories publish together.
//!
//! Original field storage, electrode/volume-charge work and physical switching
//! loss enter the existing enthalpy driver's energy gate. Endpoint numerical
//! dissipation is never heat. Fixed regions may have no physical switching loss;
//! every dissipating triangle needs an explicitly mapped finite receiver.
//!
//! Fits and mesh geometry stay fixed, and mapped caloric charts remain wholly
//! solid with constant density. Thermal transport is backward Euler and the
//! prescribed field increment uses endpoint work. Resolve electrical turns,
//! thermal events and path increments; this supplies no temperature-dependent
//! switching, electrocaloric heat, leakage, geometry evolution or measured fit.
//! Electrical refinement accumulates per-cell Q, deposited as Q/dt over one
//! backward-Euler thermal interval. It does not refine heat-release chronology
//! or increase the coupled thermal time order. [`trajectory`] instead refines
//! both histories together, resolving when that physical heat reaches each cell.

pub mod trajectory;

use crate::circuit_enthalpy::{
    CircuitEnthalpyError, CircuitEnthalpyReport, CircuitEnthalpySolve, CircuitThermalCoefficients,
    DrivenCoefficients, DriverAccount, ThermalDriver, admit_driver_unknowns, cell_driver_candidate,
    contacts::heat::{HeatAssembly, HeatNetworkError, closed_balance},
    material::{
        MaterialCell, MaterialNetwork, MaterialNetworkError, MaterialNetworkSample,
        PreparedMaterialNetwork,
    },
};
use fs_em::electrostatic::{
    ElectrostaticError,
    nonlinear::{
        DielectricMaterial, ElectrodeControl, NonlinearElectrostaticControls,
        NonlinearElectrostaticError,
        hysteretic::{
            HystereticElectrostaticStep, PlanarHystereticElectrostatics,
            path::{HystereticPathControls, HystereticPathError, HystereticPathSummary},
        },
    },
};
use fs_exec::Cx;
use fs_matdb::QueryPoint;
use fs_material::field_storage::{
    FieldStorageError, dielectric::hysteresis::vector::ResolvedVectorDielectricHysteresis,
};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::{cell::RefCell, collections::BTreeMap, sync::Arc};

/// Exactly one explicit physical receiver or fixed context per field triangle.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DielectricRegion {
    Cell {
        triangle: usize,
        cell: CellId,
    },
    /// Region outside this thermal network. Any nonzero switching heat refuses.
    Fixed {
        triangle: usize,
    },
}
impl DielectricRegion {
    fn triangle(self) -> usize {
        match self {
            Self::Cell { triangle, .. } | Self::Fixed { triangle } => triangle,
        }
    }
    fn receiver(self) -> Option<CellId> {
        match self {
            Self::Cell { cell, .. } => Some(cell),
            Self::Fixed { .. } => None,
        }
    }
}

#[derive(Debug)]
pub enum HystereticDielectricThermalError<E = HeatNetworkError> {
    Material(MaterialNetworkError),
    Dielectric(FieldStorageError),
    Field(NonlinearElectrostaticError),
    Path(HystereticPathError),
    Heat(E),
}
impl<E: core::fmt::Debug> core::fmt::Display for HystereticDielectricThermalError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spatial dielectric hysteresis/enthalpy: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for HystereticDielectricThermalError<E> {}
type Error<E = HeatNetworkError> = CircuitEnthalpyError<HystereticDielectricThermalError<E>>;

fn invalid<E>(reason: &'static str) -> Error<E> {
    CircuitEnthalpyError::Invalid(reason)
}
fn finite<E>(value: f64) -> Result<f64, Error<E>> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(
            "spatial dielectric/thermal arithmetic is unrepresentable",
        ))
    }
}
fn material<E>(error: MaterialNetworkError) -> Error<E> {
    match error {
        MaterialNetworkError::Thermal(DepositionError::Cancelled) => {
            CircuitEnthalpyError::Cancelled
        }
        error => CircuitEnthalpyError::Source(HystereticDielectricThermalError::Material(error)),
    }
}
fn dielectric<E>(error: FieldStorageError) -> Error<E> {
    CircuitEnthalpyError::Source(HystereticDielectricThermalError::Dielectric(error))
}
fn field<E>(error: NonlinearElectrostaticError) -> Error<E> {
    match error {
        NonlinearElectrostaticError::Electric(ElectrostaticError::Cancelled) => {
            CircuitEnthalpyError::Cancelled
        }
        error => CircuitEnthalpyError::Source(HystereticDielectricThermalError::Field(error)),
    }
}
fn path<E>(error: HystereticPathError) -> Error<E> {
    match error {
        HystereticPathError::Cancelled => CircuitEnthalpyError::Cancelled,
        error => CircuitEnthalpyError::Source(HystereticDielectricThermalError::Path(error)),
    }
}
fn heat<E: From<HeatNetworkError>>(error: HeatNetworkError) -> Error<E> {
    match error {
        HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
        error => CircuitEnthalpyError::Source(HystereticDielectricThermalError::Heat(error.into())),
    }
}
fn check<E>(cancelled: &RefCell<impl FnMut() -> bool>) -> Result<(), Error<E>> {
    if (*cancelled.borrow_mut())() {
        Err(CircuitEnthalpyError::Cancelled)
    } else {
        Ok(())
    }
}
fn source<E>(
    field: &PlanarHystereticElectrostatics,
    triangle: usize,
) -> Result<&ResolvedVectorDielectricHysteresis, Error<E>> {
    match field.problem().materials().get(triangle) {
        Some(DielectricMaterial::Hysteretic(material)) => Ok(material.source()),
        _ => Err(invalid(
            "retained spatial dielectric lost its original material correspondence",
        )),
    }
}
fn add<E>((sum, correction): &mut (f64, f64), value: f64) -> Result<(), Error<E>> {
    let value = finite(value)?;
    let next = finite(*sum + value)?;
    *correction = finite(
        *correction
            + if sum.abs() >= value.abs() {
                (*sum - next) + value
            } else {
                (value - next) + *sum
            },
    )?;
    *sum = next;
    Ok(())
}

struct Receiver<'n, 'a> {
    binding: &'n MaterialCell<'a>,
    volume_m3: f64,
    density_kg_m3: f64,
    interval_k: [f64; 2],
}

/// Original electrical work/accuracy policy, evaluated once per thermal step.
/// Affine paths retain the native global solve/operator budgets and empirical
/// polarization, switching-heat and work error estimates. They are not thermal
/// time-error estimates or automatic waveform-turn detection.
#[derive(Clone, Copy, Debug)]
pub enum DielectricFieldSolve {
    Endpoint(NonlinearElectrostaticControls),
    AffinePath(HystereticPathControls),
}
impl From<NonlinearElectrostaticControls> for DielectricFieldSolve {
    fn from(value: NonlinearElectrostaticControls) -> Self {
        Self::Endpoint(value)
    }
}
impl From<HystereticPathControls> for DielectricFieldSolve {
    fn from(value: HystereticPathControls) -> Self {
        Self::AffinePath(value)
    }
}

/// Accepted whole-mesh polarization with original finite-cell heat ownership.
/// Prepared geometry, heat source selections and routing persist across steps.
pub struct ThermalHystereticElectrostatics<'n, 'a> {
    field: PlanarHystereticElectrostatics,
    time_s: f64,
    network: &'n MaterialNetwork<'a>,
    prepared: Arc<PreparedMaterialNetwork<'n, 'a>>,
    regions: Arc<[DielectricRegion]>,
    receivers: Arc<BTreeMap<CellId, Receiver<'n, 'a>>>,
    volume_tolerance_m3: f64,
}
impl<'n, 'a> ThermalHystereticElectrostatics<'n, 'a> {
    /// Bind a compatible prepared field to the actual thermal clock and sources.
    /// Mapped triangles partition whole finite receiver volumes; many triangles
    /// may share one receiver. Each original constant fit must explicitly cover
    /// its receiver's complete wholly solid caloric temperature interval.
    /// `volume_tolerance_m3` is an absolute geometry/mass agreement allowance,
    /// never an inferred fill fraction or a rescaling of physical heat.
    pub fn bind(
        field: PlanarHystereticElectrostatics,
        thermal: &DepositionSimulation<'_>,
        network: &'n MaterialNetwork<'a>,
        regions: &[DielectricRegion],
        volume_tolerance_m3: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        Self::bind_with_link_conductances(
            field,
            thermal,
            network,
            regions,
            volume_tolerance_m3,
            &BTreeMap::new(),
            cancelled,
        )
    }

    /// Additional native heat links belong to a complete HeatAssembly provider.
    /// Continue with `step_candidate_with_heat_assembly`; the ordinary path
    /// refuses such unowned links rather than retaining their initial values.
    #[allow(clippy::too_many_arguments)]
    pub fn bind_with_link_conductances(
        field: PlanarHystereticElectrostatics,
        thermal: &DepositionSimulation<'_>,
        network: &'n MaterialNetwork<'a>,
        regions: &[DielectricRegion],
        volume_tolerance_m3: f64,
        initial_link_conductances_w_k: &BTreeMap<u64, f64>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        let cells = field.problem().mesh().faces().len();
        if regions.len() != cells
            || !volume_tolerance_m3.is_finite()
            || volume_tolerance_m3 < 0.0
            || !network.resistors.is_empty()
        {
            return Err(invalid(
                "complete dielectric regions, a finite volume tolerance and a thermal-only material network are required",
            ));
        }
        let prepared = network
            .prepare_with_link_conductances(thermal, &[], initial_link_conductances_w_k, || {
                (*cancel.borrow_mut())()
            })
            .map_err(material)?;
        let bindings: BTreeMap<_, _> = network
            .cells
            .iter()
            .map(|binding| (binding.cell, binding))
            .collect();
        let mut ordered = vec![None; cells];
        let mut volumes = BTreeMap::new();
        for &region in regions {
            check(&cancel)?;
            let triangle = region.triangle();
            if triangle >= cells || ordered[triangle].replace(region).is_some() {
                return Err(invalid("duplicate or foreign dielectric triangle"));
            }
            if let Some(cell) = region.receiver() {
                if !bindings.contains_key(&cell) || !thermal.cells().contains_key(&cell) {
                    return Err(invalid(
                        "dielectric receiver lacks an active original caloric binding",
                    ));
                }
                let volume = field
                    .problem()
                    .mesh()
                    .face_measure(triangle)
                    .ok_or_else(|| invalid("missing original dielectric triangle volume"))?;
                if !volume.is_finite() || volume <= 0.0 {
                    return Err(invalid(
                        "positive finite dielectric triangle volume required",
                    ));
                }
                add(volumes.entry(cell).or_insert((0.0, 0.0)), volume)?;
            }
        }
        let mut receivers = BTreeMap::new();
        for (cell, (volume, correction)) in volumes {
            check(&cancel)?;
            let volume = finite(volume + correction)?;
            let binding = bindings[&cell];
            let knots = binding.curve.knots();
            let (Some(first), Some(last)) = (knots.first(), knots.last()) else {
                return Err(invalid("dielectric receiver has no caloric interval"));
            };
            let density = first.bulk_density_kg_m3;
            if density <= 0.0
                || !density.is_finite()
                || first.temperature_k >= last.temperature_k
                || knots.iter().any(|knot| {
                    knot.liquid_mass_fraction != 0.0
                        || (knot.bulk_density_kg_m3 / density - 1.0).abs() > 1e-10
                })
            {
                return Err(invalid(
                    "fixed dielectric geometry requires a wholly solid constant-density caloric chart",
                ));
            }
            receivers.insert(
                cell,
                Receiver {
                    binding,
                    volume_m3: volume,
                    density_kg_m3: density,
                    interval_k: [first.temperature_k, last.temperature_k],
                },
            );
        }
        let regions: Vec<_> = ordered.into_iter().map(|region| region.unwrap()).collect();
        for region in &regions {
            check(&cancel)?;
            if let Some(cell) = region.receiver() {
                let receiver = &receivers[&cell];
                let original = source(&field, region.triangle())?;
                if original.material().card_identity()
                    != receiver.binding.point.material.content_hash()
                {
                    return Err(invalid(
                        "dielectric and caloric receiver must use the same original material card",
                    ));
                }
                for property in original.material().properties() {
                    check(&cancel)?;
                    receiver
                        .binding
                        .point
                        .admit_static(property)
                        .map_err(MaterialNetworkError::Material)
                        .map_err(material)?;
                }
                original
                    .scalar_source()
                    .admit_temperature_interval(
                        receiver.binding.point.temperature_axis,
                        receiver.interval_k,
                    )
                    .map_err(dielectric)?;
            }
        }
        let result = Self {
            field,
            time_s: thermal.time_s(),
            network,
            prepared: Arc::new(prepared),
            regions: regions.into(),
            receivers: Arc::new(receivers),
            volume_tolerance_m3,
        };
        let queries = result.source_queries(thermal, &cancel)?;
        result.admit_accepted_queries(&queries, &cancel)?;
        check(&cancel)?;
        Ok(result)
    }

    pub fn field(&self) -> &PlanarHystereticElectrostatics {
        &self.field
    }
    pub fn time_s(&self) -> f64 {
        self.time_s
    }
    pub fn regions(&self) -> &[DielectricRegion] {
        &self.regions
    }

    // Bind admitted every immutable constant fit over this receiver's complete
    // interval at the exact original context. Only the owned temperature varies.
    // Trial admission is therefore per receiver, not a repeated per-triangle
    // material selection; final publication still refreshes every source receipt.
    fn receiver_queries<E>(
        &self,
        thermal: &DepositionSimulation<'_>,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<BTreeMap<CellId, QueryPoint>, Error<E>> {
        let mut points = BTreeMap::new();
        for (&id, receiver) in self.receivers.iter() {
            check(cancel)?;
            let cell = thermal
                .cells()
                .get(&id)
                .ok_or_else(|| invalid("dielectric heat receiver is inactive"))?;
            let mass_volume = finite(cell.mass_kg / receiver.density_kg_m3)?;
            if (cell.geometry.volume_m3 - receiver.volume_m3).abs() > self.volume_tolerance_m3
                || (mass_volume - receiver.volume_m3).abs() > self.volume_tolerance_m3
            {
                return Err(invalid(
                    "dielectric triangles must retain their whole fixed receiving-cell volume and mass",
                ));
            }
            let phase = receiver
                .binding
                .curve
                .state_at_specific_enthalpy(cell.specific_h_j_kg)
                .map_err(MaterialNetworkError::Phase)
                .map_err(material)?;
            let temperature = phase.temperature_k();
            if !temperature.is_finite()
                || temperature < receiver.interval_k[0]
                || temperature > receiver.interval_k[1]
            {
                return Err(invalid(
                    "dielectric receiver left its original caloric temperature interval",
                ));
            }
            let point = receiver
                .binding
                .point
                .at(phase)
                .map_err(MaterialNetworkError::Material)
                .map_err(material)?;
            points.insert(id, point);
        }
        Ok(points)
    }

    fn source_queries<E>(
        &self,
        thermal: &DepositionSimulation<'_>,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<Vec<QueryPoint>, Error<E>> {
        let points = self.receiver_queries(thermal, cancel)?;
        let mut queries = Vec::with_capacity(self.regions.len());
        for &region in self.regions.iter() {
            check(cancel)?;
            let original = source(&self.field, region.triangle())?;
            let query = match region.receiver() {
                Some(cell) => points[&cell].clone(),
                None => original.material().query(),
            };
            original
                .scalar_source()
                .at_query_point(&query)
                .map_err(dielectric)?;
            queries.push(query);
        }
        Ok(queries)
    }

    fn admit_accepted_queries<E>(
        &self,
        queries: &[QueryPoint],
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<(), Error<E>> {
        for (triangle, query) in queries.iter().enumerate() {
            check(cancel)?;
            if source(&self.field, triangle)?.material().query() != *query {
                return Err(invalid(
                    "accepted dielectric queries must match the complete actual thermal state",
                ));
            }
        }
        Ok(())
    }

    /// One native endpoint or bounded affine-path calculation, then finite heat.
    /// Physical switching heat is already deposited in the returned endpoint;
    /// the numerical loss remains only in the original field account.
    /// A refined electrical path deposits its total per-cell Q/dt over this one
    /// thermal interval; thermal time accuracy still requires interval refinement.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        free_charge_density_c_m3: &[f64],
        boundary: &[ElectrodeControl],
        controls: CircuitEnthalpySolve,
        field_controls: impl Into<DielectricFieldSolve>,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<HystereticDielectricEnthalpyCandidate<'t, 'n, 'a>, Error> {
        let cancel = RefCell::new(cancelled);
        self.step_resolved_candidate(
            thermal,
            target_s,
            free_charge_density_c_m3,
            boundary,
            controls,
            field_controls.into(),
            initial_energy_increments_j,
            |state| {
                self.prepared
                    .resolve(state, || (*cancel.borrow_mut())())
                    .map_err(material)
            },
            &cancel,
        )
    }

    /// Original finite conduction, contacts, layers, radiation and coolant share
    /// this thermal endpoint. Electrical graph branches and field-dependent
    /// boundary material fits refuse: this field owner cannot supply their work.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate_with_heat_assembly<'t, A: HeatAssembly + ?Sized>(
        &self,
        cx: &Cx<'_>,
        assembly: &A,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        free_charge_density_c_m3: &[f64],
        boundary: &[ElectrodeControl],
        controls: CircuitEnthalpySolve,
        field_controls: impl Into<DielectricFieldSolve>,
        initial_energy_increments_j: &[f64],
    ) -> Result<HystereticDielectricEnthalpyCandidate<'t, 'n, 'a, A::Sample>, Error<A::Error>> {
        cx.checkpoint()
            .map_err(|_| CircuitEnthalpyError::Cancelled)?;
        let prepared = assembly.prepared_heat_network();
        let network = prepared
            .thermal_material_network()
            .map_err(heat::<A::Error>)?;
        if !std::ptr::eq(network, self.network) {
            return Err(invalid(
                "dielectric heat assembly must retain the same original material network",
            ));
        }
        prepared
            .admit_field_independent_boundaries(cx, &[])
            .map_err(heat::<A::Error>)?;
        let cancel = RefCell::new(|| cx.checkpoint().is_err());
        let candidate = self.step_resolved_candidate(
            thermal,
            target_s,
            free_charge_density_c_m3,
            boundary,
            controls,
            field_controls.into(),
            initial_energy_increments_j,
            |state| {
                let coefficients = assembly
                    .resolve_heat(cx, state)
                    .map_err(|error| error.map_source(HystereticDielectricThermalError::Heat))?;
                if !coefficients.contacts.is_empty() {
                    return Err(invalid(
                        "spatial dielectric heating cannot discard electrical contact work",
                    ));
                }
                Ok(coefficients.bulk)
            },
            &cancel,
        )?;
        closed_balance(
            target_s - thermal.time_s(),
            candidate.report.total_energy_residual_j,
            candidate.report.source.as_ref(),
            controls.total_energy_tolerance_j,
        )
        .map_err(heat::<A::Error>)?;
        cx.checkpoint()
            .map_err(|_| CircuitEnthalpyError::Cancelled)?;
        Ok(candidate)
    }

    #[allow(clippy::too_many_arguments)]
    fn step_resolved_candidate<'t, S, E>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        free_charge_density_c_m3: &[f64],
        boundary: &[ElectrodeControl],
        controls: CircuitEnthalpySolve,
        field_controls: DielectricFieldSolve,
        initial_energy_increments_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<CircuitThermalCoefficients<S>, Error<E>>,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<HystereticDielectricEnthalpyCandidate<'t, 'n, 'a, S>, Error<E>> {
        check(cancel)?;
        admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
        if thermal.time_s().to_bits() != self.time_s.to_bits()
            || !target_s.is_finite()
            || target_s <= self.time_s
        {
            return Err(invalid(
                "matching accepted dielectric/thermal clocks and a later finite target are required",
            ));
        }
        let queries = self.source_queries(thermal, cancel)?;
        self.admit_accepted_queries(&queries, cancel)?;
        if !resolve(thermal)?.resistors.is_empty() {
            return Err(invalid(
                "spatial dielectric heating cannot discard native resistor work",
            ));
        }
        // Constant fits and prescribed loading make this electrical
        // proposal temperature-independent. Never redo Gauss or advance its
        // memory during a discarded thermal Newton stencil.
        let (proposal, path_summary) = match field_controls {
            DielectricFieldSolve::Endpoint(controls) => (
                self.field
                    .solve_candidate(free_charge_density_c_m3, boundary, controls, &mut || {
                        (*cancel.borrow_mut())()
                    })
                    .map_err(field)?,
                None,
            ),
            DielectricFieldSolve::AffinePath(controls) => {
                let (proposal, summary) = self
                    .field
                    .solve_affine_path_candidate(
                        free_charge_density_c_m3,
                        boundary,
                        controls,
                        &mut || (*cancel.borrow_mut())(),
                    )
                    .map_err(path)?
                    .into_electrostatic_candidate();
                (proposal, Some(summary))
            }
        };
        let powers = self.switching_powers(proposal.step(), target_s - self.time_s, cancel)?;
        let (thermal, (), report) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            initial_energy_increments_j,
            |state, _remaining| {
                check(cancel)?;
                let _ = self.receiver_queries::<E>(state, cancel)?;
                let coefficients = resolve(state)?;
                if !coefficients.resistors.is_empty() {
                    return Err(invalid(
                        "spatial dielectric heating cannot discard native resistor work",
                    ));
                }
                Ok(DrivenCoefficients {
                    circuit: FrozenFieldDriver(proposal.step()),
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    source: coefficients.source,
                    // The existing budget counts complete source/driver
                    // evaluations. The electrical path is not recomputed.
                    circuit_trials: 1,
                })
            },
            |_, _dt, checkpoint| {
                let mut result = BTreeMap::new();
                for (&cell, &power) in &powers {
                    checkpoint()?;
                    result.insert(cell, power);
                }
                Ok(result)
            },
            || (*cancel.borrow_mut())(),
        )?;
        let queries = self.source_queries(&thermal, cancel)?;
        let (field, step) = proposal
            .with_source_queries(&queries, &mut || (*cancel.borrow_mut())())
            .map_err(field)?
            .into_parts();
        // Move the original account, not a second constitutive receipt or a
        // per-trial copy of the spatial history. Only its source queries changed.
        let report = CircuitEnthalpyReport {
            electrical: step,
            thermal: report.thermal,
            source: report.source,
            nonlinear: report.nonlinear,
            maximum_cell_residual_j: report.maximum_cell_residual_j,
            circuit_trials: report.circuit_trials,
            deposited_joule_j: report.deposited_joule_j,
            non_electrical_heat_j: report.non_electrical_heat_j,
            total_energy_residual_j: report.total_energy_residual_j,
        };
        check(cancel)?;
        Ok(HystereticDielectricEnthalpyCandidate {
            thermal,
            field: Self {
                field,
                time_s: target_s,
                network: self.network,
                prepared: Arc::clone(&self.prepared),
                regions: Arc::clone(&self.regions),
                receivers: Arc::clone(&self.receivers),
                volume_tolerance_m3: self.volume_tolerance_m3,
            },
            report,
            path_summary,
        })
    }

    fn switching_powers<E>(
        &self,
        step: &HystereticElectrostaticStep,
        dt_s: f64,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<BTreeMap<CellId, f64>, Error<E>> {
        if step.cells.len() != self.regions.len() {
            return Err(invalid(
                "switching heat lost original triangle correspondence",
            ));
        }
        let mut heat = BTreeMap::new();
        for (triangle, loss) in step.cells.iter().enumerate() {
            check(cancel)?;
            if loss.triangle != triangle
                || !loss.switching_heat_j.is_finite()
                || loss.switching_heat_j < 0.0
            {
                return Err(invalid(
                    "invalid original dielectric switching heat receipt",
                ));
            }
            match self.regions[triangle].receiver() {
                Some(cell) => add(
                    heat.entry(cell).or_insert((0.0, 0.0)),
                    loss.switching_heat_j,
                )?,
                None if loss.switching_heat_j != 0.0 => {
                    return Err(invalid(
                        "every switching dielectric triangle needs a finite heat receiver",
                    ));
                }
                None => {}
            }
        }
        let mut powers = BTreeMap::new();
        for (cell, (sum, correction)) in heat {
            check(cancel)?;
            let energy = finite(sum + correction)?;
            let power = finite(energy / dt_s)?;
            if energy > 0.0 && power == 0.0 {
                return Err(invalid(
                    "physical dielectric heat underflowed during power conversion",
                ));
            }
            powers.insert(cell, power);
        }
        Ok(powers)
    }
}

// The original accepted spatial proposal remains borrowed and immutable until
// the heat owner accepts. No field state or native receipt is cloned per trial.
struct FrozenFieldDriver<'a>(&'a HystereticElectrostaticStep);
impl ThermalDriver for FrozenFieldDriver<'_> {
    type Model = ();
    type Step = ();
    fn account(&self) -> DriverAccount<'_> {
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: self.0.switching_heat_j,
            numerical_dissipation_j: self.0.numerical_dissipation_j,
            stored_change_j: self.0.stored_change_j,
            supplied_j: self.0.endpoint_source_work_j,
        }
    }
    fn into_parts(self) -> ((), ()) {
        ((), ())
    }
}

/// Complete field and finite thermal histories, publishable together.
pub struct HystereticDielectricEnthalpyCandidate<'t, 'n, 'a, S = MaterialNetworkSample> {
    thermal: DepositionSimulation<'t>,
    field: ThermalHystereticElectrostatics<'n, 'a>,
    report: CircuitEnthalpyReport<S, HystereticElectrostaticStep>,
    path_summary: Option<HystereticPathSummary>,
}
impl<'t, 'n, 'a, S> HystereticDielectricEnthalpyCandidate<'t, 'n, 'a, S> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn field(&self) -> &ThermalHystereticElectrostatics<'n, 'a> {
        &self.field
    }
    /// `circuit_trials` counts thermal source/driver evaluations. One native
    /// field/path proposal supplies `electrical`; numerical dissipation is not heat.
    pub fn report(&self) -> &CircuitEnthalpyReport<S, HystereticElectrostaticStep> {
        &self.report
    }
    pub fn deposited_heat_j(&self) -> f64 {
        self.report.deposited_joule_j
    }
    /// Global attempted field work and accepted-only refinement estimates, when
    /// an affine path was requested. No thermal trial resets these counters.
    pub fn path_summary(&self) -> Option<&HystereticPathSummary> {
        self.path_summary.as_ref()
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        ThermalHystereticElectrostatics<'n, 'a>,
        CircuitEnthalpyReport<S, HystereticElectrostaticStep>,
        Option<HystereticPathSummary>,
    ) {
        (self.thermal, self.field, self.report, self.path_summary)
    }
}
