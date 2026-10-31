//! Backward-Euler conduction on the existing thermal network.
//!
//! The unknown is each cell's energy increment in joules. Residuals, reservoir
//! heat, source heat, entropy and final acceptance remain owned by the original
//! `BackwardEulerStep`. This module differentiates the SAME geometric Fourier
//! conductances and solves their tridiagonal Newton system with the L1 pivoted
//! linear-algebra owner. Latent plateaus retain zero dT/dh; no apparent heat
//! capacity or independently integrated temperature field is introduced.
//!
//! A bounded line search enforces constitutive domains and residual decrease.
//! Failed complete intervals may be halved, without publishing their states or
//! heat. This is nonlinear globalization, not automatic time-error control.
//! Backward Euler remains first order and requires independent time refinement.
//! The direct fast path admits fixed-mass nearest-neighbor chains in increasing
//! CellId order, with all links and reservoirs supplied by the geometric
//! conduction observer. Scheduled boundary changes occur at their exact times.
//! `ImplicitConductionGraph` admits branched, cyclic and disconnected networks
//! through the same nonlinear/time owner with a bounded sparse Krylov solve.

use super::*;
use crate::enthalpy::backward_euler::{BackwardEulerStep, BackwardEulerTrial, StepCoefficients};
use crate::enthalpy::entropy::EntropyAudit;
use fs_la::tridiagonal::{TridiagonalError, TridiagonalWorkspace};

/// Accuracy-controlled comparisons of the same original implicit intervals.
#[path = "implicit/adaptive.rs"]
pub mod adaptive;

#[path = "implicit/graph.rs"]
mod graph;
use graph::LinearTopology;
pub use graph::{ImplicitConductionGraph, ImplicitConductionLinearControls};

/// One fused constitutive query, with derivatives on the same enthalpy chart.
/// At a genuine corner the material owner explicitly selects its one-sided
/// generalized derivative. The original residual still decides acceptance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImplicitConductionState {
    pub temperature_k: f64,
    pub temperature_slope_kg_k_j: f64,
    pub conductivity_w_m_k: f64,
    pub conductivity_slope_kg_w_m_k_j: f64,
}

/// Source-owned error coordinates for adaptive time integration. A phase
/// fraction is optional because a caloric law need not represent phase change.
/// When supplied, its phase identity must be fixed by the material owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImplicitConductionObservables {
    pub temperature_k: f64,
    pub phase_mass_fraction: Option<f64>,
}

/// Immutable original caloric/conductivity objects and their analytic tangents.
/// Bounds must cover a continuous admitted caloric domain. Transport data,
/// pressure and phase closure remain explicit material-owner responsibilities.
pub trait ImplicitConductionLaw {
    fn caloric_law(&self) -> &dyn ThermalLaw;
    fn conductivity_law(&self) -> &dyn ConductivityLaw;
    fn enthalpy_bounds_j_kg(&self) -> [f64; 2];
    fn state_at_specific_enthalpy(
        &self,
        specific_h_j_kg: f64,
    ) -> Result<ImplicitConductionState, DepositionError>;

    /// Resolve temperature and any declared phase observable on this SAME
    /// caloric chart. Owners can override this to share one original inverse.
    /// Requesting phase error control without an observable refuses explicitly.
    fn time_observables(
        &self,
        specific_h_j_kg: f64,
    ) -> Result<ImplicitConductionObservables, DepositionError> {
        Ok(ImplicitConductionObservables {
            temperature_k: self
                .state_at_specific_enthalpy(specific_h_j_kg)?
                .temperature_k,
            phase_mass_fraction: None,
        })
    }
}

/// Physical interval and bounded nonlinear work, distinct from PDE accuracy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImplicitConductionControls {
    pub maximum_dt_s: f64,
    /// Smallest retried interval. A final/event-clipped interval may be shorter.
    pub minimum_dt_s: f64,
    pub maximum_steps: usize,
    pub maximum_newton_iterations: usize,
    pub maximum_backtracks: usize,
    pub maximum_halvings: usize,
    /// Whole-advance ceiling, including rejected trials and final entropy queries.
    pub maximum_residual_evaluations: usize,
    pub cell_tolerance_j: f64,
    /// Gates each interval and the independently reduced whole-call balance.
    pub total_tolerance_j: f64,
}
impl Default for ImplicitConductionControls {
    fn default() -> Self {
        Self {
            maximum_dt_s: 1.0,
            minimum_dt_s: 1e-9,
            maximum_steps: 65_536,
            maximum_newton_iterations: 32,
            maximum_backtracks: 24,
            maximum_halvings: 24,
            maximum_residual_evaluations: 1_048_576,
            cell_tolerance_j: 1e-8,
            total_tolerance_j: 1e-7,
        }
    }
}

/// Original physical accounts plus bounded nonlinear work diagnostics.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImplicitConductionReport {
    pub energy: EnergyAudit,
    pub entropy: EntropyAudit,
    pub nonlinear_iterations: usize,
    /// Krylov iterations, including work on rejected nonlinear/time candidates.
    /// Zero on the direct tridiagonal chain path.
    pub linear_iterations: usize,
    pub residual_evaluations: usize,
    pub rejected_intervals: usize,
    pub maximum_cell_residual_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum ImplicitConductionError {
    Thermal(DepositionError),
    Linear(TridiagonalError),
    Invalid(&'static str),
    NotConverged {
        maximum_residual_j: f64,
    },
    /// A graph Newton direction exhausted its restart cycles or broke down.
    KrylovNotConverged {
        /// Last completed cycle's true relative Euclidean linear residual.
        relative_residual: f64,
    },
}
impl From<DepositionError> for ImplicitConductionError {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl From<TridiagonalError> for ImplicitConductionError {
    fn from(error: TridiagonalError) -> Self {
        Self::Linear(error)
    }
}
impl core::fmt::Display for ImplicitConductionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "implicit conduction: {self:?}")
    }
}
impl std::error::Error for ImplicitConductionError {}

/// Reusable Newton storage: O(cells) for the direct chain and O(cells + links)
/// for graph coefficients. The graph's Krylov owner also uses bounded temporary
/// Arnoldi storage. Failed calls may change scratch but not the accepted timeline.
#[derive(Debug)]
pub struct ImplicitConductionWorkspace {
    linear: TridiagonalWorkspace,
    energy: Vec<f64>,
    proposal: Vec<f64>,
    direction: Vec<f64>,
    lower: Vec<f64>,
    diagonal: Vec<f64>,
    upper: Vec<f64>,
    rhs: Vec<f64>,
    // Graph off-diagonals only. Reused across Newton evaluations; chains retain
    // their three-diagonal storage and direct solve without Krylov allocation.
    graph_entries: Option<Vec<(usize, usize, f64)>>,
}
impl ImplicitConductionWorkspace {
    pub fn new(cells: usize) -> Result<Self, ImplicitConductionError> {
        fn allocate(n: usize) -> Result<Vec<f64>, ImplicitConductionError> {
            let mut v = Vec::new();
            v.try_reserve_exact(n)
                .map_err(|_| TridiagonalError::Capacity)?;
            v.resize(n, 0.0);
            Ok(v)
        }
        Ok(Self {
            linear: TridiagonalWorkspace::new(cells)?,
            energy: allocate(cells)?,
            proposal: allocate(cells)?,
            direction: allocate(cells)?,
            diagonal: allocate(cells)?,
            lower: allocate(cells.saturating_sub(1))?,
            upper: allocate(cells.saturating_sub(1))?,
            rhs: allocate(cells)?,
            graph_entries: None,
        })
    }
    pub fn cells(&self) -> usize {
        self.energy.len()
    }
}

/// Prepared structured solve on the original graph and exact physical bindings.
/// The output timeline's lifetime follows its own original caloric references,
/// not the lifetime of this optional tangent adapter.
pub struct ImplicitConductionChain<'a> {
    conduction: ConductionObserver<'a>,
    laws: BTreeMap<MaterialId, &'a dyn ImplicitConductionLaw>,
    cells: BTreeMap<CellId, BoundCell>,
    ids: Vec<CellId>,
    index: BTreeMap<CellId, usize>,
    topology: LinearTopology,
}

struct Endpoint {
    states: Vec<ImplicitConductionState>,
    link_values: BTreeMap<u64, f64>,
    reservoir_values: BTreeMap<CellId, f64>,
}

impl<'a> ImplicitConductionChain<'a> {
    pub fn new(
        thermal: &DepositionSimulation<'_>,
        conduction: &ConductionObserver<'a>,
        laws: &[(MaterialId, &'a dyn ImplicitConductionLaw)],
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ImplicitConductionError> {
        Self::prepare(thermal, conduction, laws, LinearTopology::Chain, checkpoint)
    }

    fn prepare(
        thermal: &DepositionSimulation<'_>,
        conduction: &ConductionObserver<'a>,
        laws: &[(MaterialId, &'a dyn ImplicitConductionLaw)],
        topology: LinearTopology,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ImplicitConductionError> {
        checkpoint()?;
        let ids: Vec<_> = thermal.cells.keys().copied().collect();
        if ids.is_empty() || laws.len() != thermal.materials.len() {
            return Err(ImplicitConductionError::Invalid(
                "active thermal cells and one tangent law per caloric material required",
            ));
        }
        let mut bindings = BTreeMap::new();
        for &(id, law) in laws {
            checkpoint()?;
            if bindings.insert(id, law).is_some() {
                return Err(ImplicitConductionError::Invalid(
                    "duplicate implicit material",
                ));
            }
        }
        let result = Self {
            conduction: conduction.clone(),
            laws: bindings,
            cells: thermal
                .cells
                .iter()
                .map(|(&id, cell)| {
                    (
                        id,
                        BoundCell {
                            material_id: cell.material_id,
                            body_id: cell.body_id,
                            segment_id: cell.segment_id,
                            mass_kg: cell.mass_kg,
                            geometry: cell.geometry,
                        },
                    )
                })
                .collect(),
            index: ids.iter().enumerate().map(|(i, &id)| (id, i)).collect(),
            ids,
            topology,
        };
        result.admit(thermal, &mut checkpoint)?;
        result.endpoint(thermal, &mut checkpoint)?;
        checkpoint()?;
        Ok(result)
    }

    pub fn cell_ids(&self) -> &[CellId] {
        &self.ids
    }

    fn admit(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), ImplicitConductionError> {
        checkpoint()?;
        if thermal.cells.len() != self.ids.len()
            || self.conduction.links.len() != thermal.links.len()
            || (matches!(self.topology, LinearTopology::Chain)
                && (thermal.links.len() != self.ids.len() - 1 || !thermal.covered_areas.is_empty()))
        {
            return Err(ImplicitConductionError::Invalid(
                "thermal graph bindings changed or direct-chain topology is invalid",
            ));
        }
        let mut neighbors = vec![false; self.ids.len().saturating_sub(1)];
        for (&id, binding) in &self.conduction.links {
            checkpoint()?;
            let link = thermal
                .links
                .get(&id)
                .ok_or(ImplicitConductionError::Invalid("chain link disappeared"))?;
            if link.a != binding.a || link.b != binding.b || link.area_m2 != binding.area_m2 {
                return Err(ImplicitConductionError::Invalid(
                    "chain link binding changed",
                ));
            }
            let a = *self
                .index
                .get(&link.a)
                .ok_or(ImplicitConductionError::Invalid("unknown chain endpoint"))?;
            let b = *self
                .index
                .get(&link.b)
                .ok_or(ImplicitConductionError::Invalid("unknown chain endpoint"))?;
            if matches!(self.topology, LinearTopology::Graph(_)) {
                continue;
            }
            if a.abs_diff(b) != 1 {
                return Err(ImplicitConductionError::Invalid(
                    "chain links must join consecutive CellIds",
                ));
            }
            if neighbors[a.min(b)] {
                return Err(ImplicitConductionError::Invalid(
                    "repeated nearest-neighbor link",
                ));
            }
            neighbors[a.min(b)] = true;
        }
        if matches!(self.topology, LinearTopology::Chain)
            && neighbors.iter().any(|&present| !present)
        {
            return Err(ImplicitConductionError::Invalid(
                "each neighboring cell pair needs exactly one link",
            ));
        }
        for &id in &self.ids {
            checkpoint()?;
            let cell = thermal
                .cells
                .get(&id)
                .ok_or(ImplicitConductionError::Invalid("chain cell disappeared"))?;
            let law = self
                .laws
                .get(&cell.material_id)
                .ok_or(ImplicitConductionError::Invalid(
                    "missing implicit material binding",
                ))?;
            if !thermal.uses_material_law(cell.material_id, law.caloric_law()) {
                return Err(ImplicitConductionError::Invalid(
                    "implicit caloric object differs from original timeline",
                ));
            }
            if let Some(k) = self.conduction.laws.get(&cell.material_id) {
                if !std::ptr::addr_eq(
                    *k as *const dyn ConductivityLaw,
                    law.conductivity_law() as *const dyn ConductivityLaw,
                ) {
                    return Err(ImplicitConductionError::Invalid(
                        "implicit conductivity object differs from original transport",
                    ));
                }
            } else if thermal.links.len() != 0 || cell.boundary.reservoir.is_some() {
                return Err(ImplicitConductionError::Invalid(
                    "missing original conductivity binding",
                ));
            }
            if let Some(bound) = self.cells.get(&id) {
                if cell.material_id != bound.material_id
                    || cell.mass_kg != bound.mass_kg
                    || cell.body_id != bound.body_id
                    || cell.segment_id != bound.segment_id
                    || cell.geometry != bound.geometry
                {
                    return Err(ImplicitConductionError::Invalid(
                        "implicit reference cell binding changed",
                    ));
                }
            }
            let [lo, hi] = law.enthalpy_bounds_j_kg();
            if !lo.is_finite()
                || !hi.is_finite()
                || lo >= hi
                || !(lo..=hi).contains(&cell.specific_h_j_kg)
            {
                return Err(DepositionError::MaterialDomain.into());
            }
            match (cell.boundary.reservoir, self.conduction.reservoirs.get(&id)) {
                (Some(r), Some(b)) if r.area_m2 == b.area_m2 => {}
                (None, None) => {}
                _ => {
                    return Err(ImplicitConductionError::Invalid(
                        "implicit reservoir binding changed",
                    ));
                }
            }
        }
        checkpoint()?;
        Ok(())
    }

    fn endpoint(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Endpoint, ImplicitConductionError> {
        let mut states = Vec::with_capacity(self.ids.len());
        for &id in &self.ids {
            checkpoint()?;
            let cell = &thermal.cells[&id];
            let law = self.laws[&cell.material_id];
            let [lo, hi] = law.enthalpy_bounds_j_kg();
            if !(lo..=hi).contains(&cell.specific_h_j_kg) {
                return Err(DepositionError::MaterialDomain.into());
            }
            let state = law.state_at_specific_enthalpy(cell.specific_h_j_kg)?;
            let max_slope = thermal.material_temperature_slopes[&cell.material_id];
            if !positive(state.temperature_k)
                || !nonnegative(state.temperature_slope_kg_k_j)
                || state.temperature_slope_kg_k_j > max_slope
                || !nonnegative(state.conductivity_w_m_k)
                || !state.conductivity_slope_kg_w_m_k_j.is_finite()
                || (state.conductivity_w_m_k == 0.0 && state.conductivity_slope_kg_w_m_k_j != 0.0)
            {
                return Err(DepositionError::MaterialDomain.into());
            }
            states.push(state);
        }
        let mut link_values = BTreeMap::new();
        for (&id, binding) in &self.conduction.links {
            checkpoint()?;
            let a = states[self.index[&binding.a]];
            let b = states[self.index[&binding.b]];
            link_values.insert(
                id,
                series_contact_conductance_w_k(
                    binding.area_m2,
                    binding.geometry.half_distances_m,
                    [a.conductivity_w_m_k, b.conductivity_w_m_k],
                    binding.geometry.contact_resistance_m2_k_w,
                )?,
            );
        }
        let mut reservoir_values = BTreeMap::new();
        for (&id, binding) in &self.conduction.reservoirs {
            checkpoint()?;
            let k = states[self.index[&id]].conductivity_w_m_k;
            let g = if k == 0.0 || binding.area_m2 == 0.0 {
                0.0
            } else {
                let bulk = finite(binding.geometry.half_distance_m / k)?;
                let resistance = finite(bulk + binding.geometry.surface_resistance_m2_k_w)?;
                let value = finite(binding.area_m2 / resistance)?;
                if !positive(bulk) || !positive(resistance) || !positive(value) {
                    return Err(DepositionError::Arithmetic.into());
                }
                value
            };
            reservoir_values.insert(id, g);
        }
        checkpoint()?;
        Ok(Endpoint {
            states,
            link_values,
            reservoir_values,
        })
    }

    fn evaluate<'m>(
        &self,
        step: &BackwardEulerStep<'_, 'm>,
        energy: &[f64],
        entropy: bool,
        controls: ImplicitConductionControls,
        report: &mut ImplicitConductionReport,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<BackwardEulerTrial<'m, Endpoint>, ImplicitConductionError> {
        if report.residual_evaluations >= controls.maximum_residual_evaluations {
            return Err(DepositionError::Budget("implicit residual evaluations").into());
        }
        report.residual_evaluations += 1;
        // Share the bounded checkpoint with the original owner and its nested
        // constitutive callback. Calls are sequential, with no retained borrow.
        let checkpoint = std::cell::RefCell::new(checkpoint);
        let resolve = |thermal: &DepositionSimulation<'m>| {
            let endpoint = self.endpoint(thermal, &mut || (checkpoint.borrow_mut())())?;
            let reservoirs = endpoint.reservoir_values.clone();
            Ok::<_, ImplicitConductionError>((
                StepCoefficients {
                    sources_w: BTreeMap::new(),
                    conductances_w_k: endpoint.link_values.clone(),
                    payload: endpoint,
                },
                reservoirs,
            ))
        };
        if entropy {
            step.evaluate_with_entropy_and_reservoirs(energy, resolve, || {
                (checkpoint.borrow_mut())()
            })
        } else {
            step.evaluate_with_reservoirs(energy, resolve, || (checkpoint.borrow_mut())())
        }
    }

    fn jacobian(
        &self,
        thermal: &DepositionSimulation<'_>,
        endpoint: &Endpoint,
        dt: f64,
        workspace: &mut ImplicitConductionWorkspace,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), ImplicitConductionError> {
        workspace.diagonal.fill(1.0);
        workspace.lower.fill(0.0);
        workspace.upper.fill(0.0);
        match self.topology {
            LinearTopology::Chain => workspace.graph_entries = None,
            LinearTopology::Graph(_) => {
                let entries = workspace.graph_entries.get_or_insert_with(Vec::new);
                entries.clear();
                entries
                    .try_reserve(
                        self.conduction
                            .links
                            .len()
                            .checked_mul(2)
                            .ok_or(TridiagonalError::Capacity)?,
                    )
                    .map_err(|_| TridiagonalError::Capacity)?;
            }
        }
        let mut covered: Vec<f64> = self
            .ids
            .iter()
            .map(|id| {
                thermal.cells[id]
                    .boundary
                    .reservoir
                    .map_or(0.0, |r| r.area_m2)
            })
            .collect();
        for (&id, link) in &self.conduction.links {
            checkpoint()?;
            let ia = self.index[&link.a];
            let ib = self.index[&link.b];
            let a = endpoint.states[ia];
            let b = endpoint.states[ib];
            let g = endpoint.link_values[&id];
            let ga = conductance_derivative(g, link.area_m2, link.geometry.half_distances_m[0], a)?;
            let gb = conductance_derivative(g, link.area_m2, link.geometry.half_distances_m[1], b)?;
            let difference = b.temperature_k - a.temperature_k;
            let qa = finite(
                (ga * difference - g * a.temperature_slope_kg_k_j) / thermal.cells[&link.a].mass_kg,
            )?;
            let qb = finite(
                (gb * difference + g * b.temperature_slope_kg_k_j) / thermal.cells[&link.b].mass_kg,
            )?;
            add_entry(workspace, ia, ia, -dt * qa)?;
            add_entry(workspace, ia, ib, -dt * qb)?;
            add_entry(workspace, ib, ia, dt * qa)?;
            add_entry(workspace, ib, ib, dt * qb)?;
            let coverage = thermal.contact_areas(&thermal.links[&id]);
            covered[ia] = finite(covered[ia] + coverage[0])?;
            covered[ib] = finite(covered[ib] + coverage[1])?;
        }
        for (i, &id) in self.ids.iter().enumerate() {
            checkpoint()?;
            let cell = &thermal.cells[&id];
            let state = endpoint.states[i];
            let area = finite(cell.geometry.surface_area_m2 - covered[i])?;
            if area < 0.0 {
                return Err(DepositionError::Arithmetic.into());
            }
            let mut derivative =
                -(area * cell.boundary.htc_w_m2_k) * state.temperature_slope_kg_k_j;
            if let Some(reservoir) = cell.boundary.reservoir {
                let binding = self.conduction.reservoirs[&id];
                let g = endpoint.reservoir_values[&id];
                let dg = conductance_derivative(
                    g,
                    binding.area_m2,
                    binding.geometry.half_distance_m,
                    state,
                )?;
                derivative += dg * (reservoir.temperature_k - state.temperature_k)
                    - g * state.temperature_slope_kg_k_j;
            }
            add_entry(workspace, i, i, -dt * derivative / cell.mass_kg)?;
        }
        checkpoint()?;
        Ok(())
    }

    fn solve_step<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        controls: ImplicitConductionControls,
        workspace: &mut ImplicitConductionWorkspace,
        report: &mut ImplicitConductionReport,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'m>, EnergyAudit, EntropyAudit, f64), ImplicitConductionError>
    {
        let step = BackwardEulerStep::new(thermal, target_s, &mut *checkpoint)?;
        workspace.energy.fill(0.0);
        let mut trial = self.evaluate(
            &step,
            &workspace.energy,
            false,
            controls,
            report,
            checkpoint,
        )?;
        for iteration in 0..=controls.maximum_newton_iterations {
            checkpoint()?;
            let maximum = trial.maximum_cell_residual_j();
            if maximum <= controls.cell_tolerance_j
                && trial.audit().residual_j.abs() <= controls.total_tolerance_j
            {
                let final_trial =
                    self.evaluate(&step, &workspace.energy, true, controls, report, checkpoint)?;
                let maximum = final_trial.maximum_cell_residual_j();
                let entropy =
                    *final_trial
                        .entropy_audit()
                        .ok_or(ImplicitConductionError::Invalid(
                            "missing final entropy account",
                        ))?;
                let (next, _, audit) = final_trial.accept(
                    controls.cell_tolerance_j,
                    controls.total_tolerance_j,
                    &mut *checkpoint,
                )?;
                return Ok((next, audit, entropy, maximum));
            }
            if iteration == controls.maximum_newton_iterations {
                return Err(ImplicitConductionError::NotConverged {
                    maximum_residual_j: maximum,
                });
            }
            self.jacobian(
                thermal,
                trial.payload(),
                step.duration_s(),
                workspace,
                checkpoint,
            )?;
            for (rhs, &residual) in workspace.rhs.iter_mut().zip(trial.residual_j()) {
                *rhs = -residual;
            }
            match self.topology {
                LinearTopology::Chain => {
                    let mut cancellation_error = None;
                    let solved = workspace.linear.solve_into(
                        &workspace.lower,
                        &workspace.diagonal,
                        &workspace.upper,
                        &workspace.rhs,
                        &mut workspace.direction,
                        || match checkpoint() {
                            Ok(()) => false,
                            Err(e) => {
                                cancellation_error = Some(e);
                                true
                            }
                        },
                    );
                    if let Some(error) = cancellation_error {
                        return Err(error.into());
                    }
                    solved?;
                }
                LinearTopology::Graph(linear) => {
                    graph::solve(workspace, linear, report, checkpoint)?;
                }
            }
            report.nonlinear_iterations += 1;
            let mut accepted = None;
            let mut fraction = 1.0;
            for _ in 0..=controls.maximum_backtracks {
                checkpoint()?;
                let mut representable = false;
                let mut admitted = true;
                for (i, &id) in self.ids.iter().enumerate() {
                    let mut value = workspace.energy[i] + fraction * workspace.direction[i];
                    let cell = &thermal.cells[&id];
                    let h = cell.specific_h_j_kg + value / cell.mass_kg;
                    let [lo, hi] = self.laws[&cell.material_id].enthalpy_bounds_j_kg();
                    admitted &= value.is_finite() && h.is_finite() && (lo..=hi).contains(&h);
                    // A diffusion resolvent has arbitrarily small far-field
                    // increments. If this proposal cannot change the stored
                    // coordinate, propose its exactly representable zero.
                    // This does NOT discard a heat/source term: the original
                    // residual still charges all endpoint heat to this cell,
                    // and both local and whole-call joule gates remain intact.
                    if h == cell.specific_h_j_kg {
                        value = 0.0;
                    }
                    representable |= value != workspace.energy[i];
                    workspace.proposal[i] = value;
                }
                if admitted && representable {
                    match self.evaluate(
                        &step,
                        &workspace.proposal,
                        false,
                        controls,
                        report,
                        checkpoint,
                    ) {
                        Ok(candidate) => {
                            let value = candidate.maximum_cell_residual_j();
                            if value <= (1.0 - 1e-4 * fraction) * maximum
                                || (value <= controls.cell_tolerance_j
                                    && candidate.audit().residual_j.abs()
                                        <= controls.total_tolerance_j)
                            {
                                accepted = Some(candidate);
                                break;
                            }
                        }
                        Err(ImplicitConductionError::Thermal(
                            DepositionError::MaterialDomain | DepositionError::Arithmetic,
                        )) => {}
                        Err(error) => return Err(error),
                    }
                }
                fraction *= 0.5;
            }
            let Some(candidate) = accepted else {
                return Err(ImplicitConductionError::NotConverged {
                    maximum_residual_j: maximum,
                });
            };
            workspace.energy.copy_from_slice(&workspace.proposal);
            trial = candidate;
        }
        unreachable!("bounded nonlinear loop returns at its limit")
    }

    /// Advance the original timeline with endpoint-evaluated nonlinear heat.
    /// Whole-interval retries and final joule gates are transactional. No
    /// explicit diffusion CFL restriction is imposed on an accepted interval.
    pub fn advance_to<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        controls: ImplicitConductionControls,
        workspace: &mut ImplicitConductionWorkspace,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'m>, ImplicitConductionReport), ImplicitConductionError> {
        checkpoint()?;
        validate_advance(
            thermal.time_s,
            target_s,
            self.ids.len(),
            controls,
            workspace.cells(),
        )?;
        self.admit(thermal, &mut checkpoint)?;
        let mut next = thermal.clone();
        let mut observer = self.conduction.clone();
        let mut report = empty_report();
        let mut supplied = Sum::default();
        loop {
            checkpoint()?;
            // The original owner applies due events in exact authored order.
            let (with_events, observed, events) = next.advance_observed(
                &observer,
                next.time_s,
                StepConfig {
                    max_dt_s: controls.maximum_dt_s,
                    energy_tolerance_j: controls.total_tolerance_j,
                },
                &mut checkpoint,
            )?;
            next = with_events;
            observer = observed;
            report.energy.events += events.events;
            self.admit(&next, &mut checkpoint)?;
            if next.time_s == target_s {
                break;
            }
            if report.energy.substeps >= controls.maximum_steps.min(next.budget.max_substeps) {
                return Err(DepositionError::Budget("implicit thermal intervals").into());
            }
            let stop = next
                .schedule
                .get(next.cursor)
                .map_or(target_s, |e| e.time_s.min(target_s));
            let mut dt = controls.maximum_dt_s.min(stop - next.time_s);
            let mut halvings = 0;
            loop {
                checkpoint()?;
                let end = if dt == stop - next.time_s {
                    stop
                } else {
                    next.time_s + dt
                };
                if !end.is_finite() || end <= next.time_s {
                    return Err(DepositionError::Arithmetic.into());
                }
                match self.solve_step(
                    &next,
                    end,
                    controls,
                    workspace,
                    &mut report,
                    &mut checkpoint,
                ) {
                    Ok((accepted, energy, entropy, maximum)) => {
                        supplied.add(energy.external_heat_j)?;
                        accumulate_entropy(&mut report.entropy, entropy)?;
                        report.maximum_cell_residual_j =
                            report.maximum_cell_residual_j.max(maximum);
                        report.energy.substeps += 1;
                        next = accepted;
                        break;
                    }
                    Err(
                        error @ (ImplicitConductionError::NotConverged { .. }
                        | ImplicitConductionError::KrylovNotConverged { .. }
                        | ImplicitConductionError::Linear(TridiagonalError::Singular {
                            ..
                        })),
                    ) => {
                        report.rejected_intervals += 1;
                        if halvings >= controls.maximum_halvings || dt * 0.5 < controls.minimum_dt_s
                        {
                            return Err(error);
                        }
                        halvings += 1;
                        dt *= 0.5;
                    }
                    Err(error) => return Err(error),
                }
            }
        }
        let mut storage = Sum::default();
        for &id in &self.ids {
            checkpoint()?;
            storage.add(
                thermal.cells[&id].mass_kg
                    * (next.cells[&id].specific_h_j_kg - thermal.cells[&id].specific_h_j_kg),
            )?;
        }
        report.energy.external_heat_j = supplied.total()?;
        report.energy.stored_energy_change_j = storage.total()?;
        report.energy.residual_j =
            finite(report.energy.stored_energy_change_j - report.energy.external_heat_j)?;
        if report.energy.residual_j.abs() > controls.total_tolerance_j {
            return Err(DepositionError::EnergyBalance(report.energy.residual_j).into());
        }
        checkpoint()?;
        Ok((next, report))
    }
}

fn validate_advance(
    start_s: f64,
    target_s: f64,
    cells: usize,
    controls: ImplicitConductionControls,
    workspace_cells: usize,
) -> Result<(), ImplicitConductionError> {
    if !target_s.is_finite()
        || target_s < start_s
        || !positive(controls.maximum_dt_s)
        || !positive(controls.minimum_dt_s)
        || controls.minimum_dt_s > controls.maximum_dt_s
        || controls.maximum_steps == 0
        || controls.maximum_newton_iterations == 0
        || controls.maximum_newton_iterations > 128
        || controls.maximum_backtracks > 32
        || controls.maximum_halvings > 32
        || controls.maximum_residual_evaluations == 0
        || !nonnegative(controls.cell_tolerance_j)
        || !nonnegative(controls.total_tolerance_j)
        || workspace_cells != cells
    {
        return Err(ImplicitConductionError::Invalid(
            "invalid implicit interval, work, accuracy or workspace controls",
        ));
    }
    Ok(())
}

fn conductance_derivative(
    g: f64,
    area: f64,
    distance: f64,
    state: ImplicitConductionState,
) -> Result<f64, DepositionError> {
    if g == 0.0 || state.conductivity_slope_kg_w_m_k_j == 0.0 {
        return Ok(0.0);
    }
    // G*(bulk resistance / total resistance)*(k'/k), keeping the geometric
    // resistance ratio bounded rather than directly forming G²/k².
    finite(
        g * finite((g / area) * (distance / state.conductivity_w_m_k))?
            * (state.conductivity_slope_kg_w_m_k_j / state.conductivity_w_m_k),
    )
}
fn add_entry(
    w: &mut ImplicitConductionWorkspace,
    row: usize,
    col: usize,
    value: f64,
) -> Result<(), DepositionError> {
    let target = if row == col {
        &mut w.diagonal[row]
    } else if let Some(entries) = &mut w.graph_entries {
        entries.push((row, col, finite(value)?));
        return Ok(());
    } else if row == col + 1 {
        &mut w.lower[col]
    } else if col == row + 1 {
        &mut w.upper[row]
    } else {
        return Err(DepositionError::InvalidInput(
            "nontridiagonal conduction Jacobian",
        ));
    };
    *target = finite(*target + finite(value)?)?;
    Ok(())
}
fn empty_report() -> ImplicitConductionReport {
    ImplicitConductionReport {
        energy: EnergyAudit {
            incoming_mass_kg: 0.0,
            incoming_energy_j: 0.0,
            external_heat_j: 0.0,
            stored_energy_change_j: 0.0,
            residual_j: 0.0,
            substeps: 0,
            events: 0,
        },
        entropy: EntropyAudit::default(),
        nonlinear_iterations: 0,
        linear_iterations: 0,
        residual_evaluations: 0,
        rejected_intervals: 0,
        maximum_cell_residual_j: 0.0,
    }
}
fn accumulate_entropy(total: &mut EntropyAudit, next: EntropyAudit) -> Result<(), DepositionError> {
    total.thermal_storage_change_j_k =
        finite(total.thermal_storage_change_j_k + next.thermal_storage_change_j_k)?;
    total.reservoir_entropy_input_j_k =
        finite(total.reservoir_entropy_input_j_k + next.reservoir_entropy_input_j_k)?;
    total.conduction_generation_j_k =
        finite(total.conduction_generation_j_k + next.conduction_generation_j_k)?;
    total.unclassified_source_at_cell_j_k =
        finite(total.unclassified_source_at_cell_j_k + next.unclassified_source_at_cell_j_k)?;
    total.has_unclassified_sources |= next.has_unclassified_sources;
    total.temporal_defect_j_k = finite(
        total.thermal_storage_change_j_k
            - total.reservoir_entropy_input_j_k
            - total.conduction_generation_j_k
            - total.unclassified_source_at_cell_j_k,
    )?;
    Ok(())
}
#[derive(Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) -> Result<(), DepositionError> {
        let value = finite(value)?;
        let next = finite(self.value + value)?;
        let lost = if self.value.abs() >= value.abs() {
            (self.value - next) + value
        } else {
            (value - next) + self.value
        };
        self.correction = finite(self.correction + lost)?;
        self.value = next;
        Ok(())
    }
    fn total(self) -> Result<f64, DepositionError> {
        finite(self.value + self.correction)
    }
}

#[cfg(test)]
#[path = "implicit/tests.rs"]
mod tests;
