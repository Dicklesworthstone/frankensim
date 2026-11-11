//! Implicit enthalpy conduction on the existing thermal network.
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
//! conduction observer. Selected ambient faces share its material/film law;
//! unselected ambient channels retain authored htc. Bath events occur at their exact times.
//! `ImplicitConductionGraph` admits branched, cyclic and disconnected networks
//! through the same nonlinear/time owner with a bounded sparse Krylov solve.
//! The opt-in `advance_sdirk2_to` path shares this residual and Newton owner,
//! with second-order L-stable quadrature for fixed networks, piecewise constant
//! source/bath events and continuous prescribed source and bath histories,
//! including native Robin coefficients.
//! Selected planar layers may expand or contract at fixed mass and footprint.
//! Their original density and specific-volume tangent update every native
//! Fourier/contact/surface resistance and its Jacobian, including on latent
//! plateaus. Custom flux owners may explicitly supply the same geometry chain
//! rule; owners that supply only fixed-geometry tangents still refuse.
//! Independently supplied current mechanical face measures are also supported,
//! held fixed through the interval's Newton trials and SDIRK stages. Reference
//! mass and coverage stay authoritative; mechanical and density geometry do
//! not compose implicitly, and coupling still requires time refinement.

use super::*;
use crate::enthalpy::backward_euler::{
    BackwardEulerStep, BackwardEulerTrial, BoundaryConductances, BoundaryExchangeStepCoefficients,
    BoundaryExchanges, BoundaryStepCoefficients, StageQuadrature, StepCoefficients,
    TransportStepCoefficients, resolvable_increment_j,
};
use crate::enthalpy::entropy::EntropyAudit;
use crate::enthalpy::SourcePowers;
use fs_la::tridiagonal::{TridiagonalError, TridiagonalWorkspace};

/// Continuous small-signal thermal impedance at an admitted equilibrium.
#[path = "implicit/harmonic.rs"]
pub mod harmonic;

/// Direct stationary heat balance on the same material and transport owners.
#[path = "implicit/equilibrium.rs"]
pub mod equilibrium;

/// Accuracy-controlled comparisons of the same original implicit intervals.
#[path = "implicit/adaptive.rs"]
pub mod adaptive;

/// Second-order native stages for fixed networks and scheduled boundary data.
#[path = "implicit/sdirk.rs"]
mod sdirk;

#[path = "implicit/graph.rs"]
mod graph;
use graph::LinearTopology;
pub use graph::{ImplicitConductionGraph, ImplicitConductionLinearControls};

/// IMEX sources and conservative transport with transactional observer histories.
#[path = "implicit/observed.rs"]
mod observed;
pub use observed::ImplicitConductionObservedError;

/// Local heat-generation feedback evaluated with the conductive endpoint.
#[path = "implicit/source.rs"]
pub mod source;
use source::{ImplicitHeatSourceLaw, ImplicitHeatSourceState};

/// Sparse cross-cell heat-generation feedback in the same native residual.
#[path = "implicit/coupled_source.rs"]
pub mod coupled_source;
use coupled_source::{CoupledHeatSource, ImplicitCoupledHeatSourceState};

/// Source-owned state-dependent contact resistance in the same implicit solve.
#[path = "implicit/contact.rs"]
pub mod contact;
use contact::{ImplicitContactLaw, ImplicitContactState};

/// Integrated half-cell fluxes and their complete two-sided tangents.
#[path = "implicit/link.rs"]
pub mod link;
use link::{ImplicitLinkLaw, ImplicitLinkState};

/// Integrated material-to-bath fluxes on original reservoir and ambient faces.
#[path = "implicit/boundary.rs"]
pub mod boundary;
use boundary::{ImplicitBoundaryGeometry, ImplicitBoundaryLaw, ImplicitBoundaryState,
    ImplicitBoundaryExchangeLaw, ImplicitBoundaryExchangeState};

/// Endpoint upwind enthalpy in the same conservative nonlinear solve.
#[path = "implicit/advection.rs"]
mod advection;
use advection::{ClassifiedAdvectionEntropy, ImplicitAdvection, ImplicitAdvectionStage};

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

    /// Original specific-volume derivative `d(1/rho)/dh` [m³/J]. This is
    /// required only for material-density planar layers; `None` explicitly
    /// leaves that geometry coupling unsupported. The density value still
    /// comes from the SAME original [`ThermalLaw::bulk_density_kg_m3`] owner.
    ///
    /// Use the admitted caloric chart's one-sided derivative at corners.
    /// A latent plateau may have nonzero volume slope despite zero `dT/dh`.
    /// This changes resistance at conserved mass, not storage or pressure work.
    fn specific_volume_derivative_at_specific_enthalpy(
        &self,
        _specific_h_j_kg: f64,
    ) -> Result<Option<f64>, DepositionError> {
        Ok(None)
    }

    /// Admit a two-sided local harmonic operating point on the original law.
    /// A finite-step one-sided tangent at a phase boundary need not define a
    /// harmonic derivative. Owners with known corners must refuse them here.
    /// The default preserves the existing caller obligation to choose a smooth
    /// branch; it is not a differentiability certificate. Time stepping and
    /// stationary root finding do not use this additional response gate.
    fn admit_harmonic_operating_point(
        &self,
        _specific_h_j_kg: f64,
    ) -> Result<(), DepositionError> {
        Ok(())
    }

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
    /// The whole-call joule budget, gating the independently reduced whole-call
    /// balance. Each interval's Newton stop and gate take its allowance: the
    /// budget accrued through the interval's end in proportion to time,
    /// `total_tolerance_j * (end - call start) / call span`, less the residual
    /// magnitudes already accepted in the call. N accepted interval residuals
    /// thus sum in magnitude to at most this budget times `1 + 2 * 2^-53` to
    /// first order in 2^-53 (the rounding of the allowances and of their
    /// compensated running sum), and a short interval inherits what the
    /// intervals before it left unspent.
    pub total_tolerance_j: f64,
}
impl ImplicitConductionControls {
    // An interval's controls: the same limits, with its allowance (see
    // `total_tolerance_j`) as the joule tolerance. `spent_j` is the compensated
    // sum of the residual magnitudes accepted earlier in the call; at the call's
    // target the accrued budget is `total_tolerance_j` exactly. A negative
    // difference, possible only by rounding, is an allowance of zero.
    fn interval_allowance(
        self,
        call_start_s: f64,
        call_target_s: f64,
        end_s: f64,
        spent_j: f64,
    ) -> Result<Self, DepositionError> {
        let fraction =
            finite(finite(end_s - call_start_s)? / finite(call_target_s - call_start_s)?)?;
        let accrued = finite(self.total_tolerance_j * fraction)?;
        Ok(Self {
            total_tolerance_j: finite(accrued - spent_j)?.max(0.0),
            ..self
        })
    }
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
    /// A passive SDIRK2 interval loses cell-plus-bath entropy beyond its
    /// local joule-residual and floating-point allowance. Retrying may reduce
    /// the physical interval; neither state nor entropy is repaired.
    EntropyDecrease {
        entropy_change_j_k: f64,
        allowance_j_k: f64,
    },
    /// A graph Newton direction exhausted its restart cycles or broke down.
    KrylovNotConverged {
        /// Last completed cycle's true relative Euclidean linear residual.
        relative_residual: f64,
    },
    /// An interval's allowance of `total_tolerance_j` is below the rounding floor
    /// of its own energy audit, with every cell already within `cell_tolerance_j`:
    /// f64 cannot resolve what is left of the whole-call budget at this interval.
    /// Raised at that interval, which publishes nothing; it is not retried.
    AllowanceBelowAuditFloor {
        /// The interval's allowance of `total_tolerance_j`.
        allowance_j: f64,
        /// The rounding floor of the refused trial's energy audit.
        floor_j: f64,
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

/// Reusable Newton storage: O(cells) for the direct chain and
/// O(cells + links + source dependencies) for graph coefficients. The graph's Krylov owner also uses bounded temporary
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
    contacts: BTreeMap<u64, &'a dyn ImplicitContactLaw>,
    link_laws: BTreeMap<u64, &'a dyn ImplicitLinkLaw>,
    reservoir_laws: BTreeMap<CellId, &'a dyn ImplicitBoundaryLaw>,
    ambient_laws: BTreeMap<CellId, &'a dyn ImplicitBoundaryLaw>,
    reservoir_exchanges: BTreeMap<CellId, &'a dyn ImplicitBoundaryExchangeLaw>,
    ambient_exchanges: BTreeMap<CellId, &'a dyn ImplicitBoundaryExchangeLaw>,
    reservoir_histories: BTreeMap<CellId, boundary::BoundBathHistory<'a>>,
    ambient_histories: BTreeMap<CellId, boundary::BoundBathHistory<'a>>,
    heat_sources: BTreeMap<CellId, &'a dyn ImplicitHeatSourceLaw>,
    coupled_heat_source: Option<CoupledHeatSource<'a>>,
    advection: Option<ImplicitAdvection>,
}

struct Endpoint {
    states: Vec<ImplicitConductionState>,
    // Empty for fixed geometry. Selected layers share one density/tangent
    // query per cell across all normal links and external faces.
    layer_distances: Vec<Option<ImplicitHalfDistance>>,
    link_values: BTreeMap<u64, f64>,
    boundary_values: BoundaryConductances,
    exchange_values: BoundaryExchanges,
    contacts: BTreeMap<u64, ImplicitContactState>,
    link_states: BTreeMap<u64, ImplicitLinkState>,
    reservoir_states: BTreeMap<CellId, ImplicitBoundaryState>,
    ambient_states: BTreeMap<CellId, ImplicitBoundaryState>,
    reservoir_exchange_states: BTreeMap<CellId, ImplicitBoundaryExchangeState>,
    ambient_exchange_states: BTreeMap<CellId, ImplicitBoundaryExchangeState>,
    heat_source_states: BTreeMap<CellId, ImplicitHeatSourceState>,
    coupled_heat_source_state: Option<ImplicitCoupledHeatSourceState>,
    transport_tangent: Vec<(usize, usize, f64)>,
    advection_multiplier: f64,
    advection_entropy: Option<ClassifiedAdvectionEntropy>,
    // Set only by frozen-clock stationary/response preparation. Evolution
    // continues to evaluate and account its own native endpoint fluxes.
    stationary_advection: Option<advection::StationaryAdvection>,
}

#[derive(Clone, Copy)]
struct ImplicitHalfDistance {
    value_m: f64,
    /// `(dd/dh)/d = rho * d(1/rho)/dh` at fixed mass and footprint.
    relative_enthalpy_slope_kg_j: f64,
}

fn endpoint_distance(
    layers: &[Option<ImplicitHalfDistance>],
    index: usize,
    reference_m: f64,
) -> ImplicitHalfDistance {
    layers
        .get(index)
        .copied()
        .flatten()
        .unwrap_or(ImplicitHalfDistance {
            value_m: reference_m,
            relative_enthalpy_slope_kg_j: 0.0,
        })
}

fn require_geometry_tangent<T>(state: Option<T>) -> Result<T, ImplicitConductionError> {
    state.ok_or(ImplicitConductionError::Invalid(
        "integrated flux laws do not supply material-density geometry derivatives",
    ))
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
        let mut result = Self {
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
                            phase_layer: None,
                        },
                    )
                })
                .collect(),
            index: ids.iter().enumerate().map(|(i, &id)| (id, i)).collect(),
            ids,
            topology,
            contacts: BTreeMap::new(),
            link_laws: BTreeMap::new(),
            reservoir_laws: BTreeMap::new(),
            ambient_laws: BTreeMap::new(),
            reservoir_exchanges: BTreeMap::new(),
            ambient_exchanges: BTreeMap::new(),
            reservoir_histories: BTreeMap::new(),
            ambient_histories: BTreeMap::new(),
            heat_sources: BTreeMap::new(),
            coupled_heat_source: None,
            advection: None,
        };
        if conduction.has_phase_layers || conduction.current_geometry.is_some() {
            // Retain original geometric admission, including reference mesh,
            // material, inventory and coverage, before fixing this solve's mass.
            result.conduction.prepare(thermal, &mut checkpoint)?;
        }
        result.admit(thermal, &mut checkpoint)?;
        result.endpoint(thermal, &mut checkpoint)?;
        checkpoint()?;
        Ok(result)
    }

    pub fn cell_ids(&self) -> &[CellId] {
        &self.ids
    }

    fn has_material_density_geometry(&self, id: CellId) -> bool {
        self.conduction.has_material_density_layers
            && self
                .conduction
                .cells
                .get(&id)
                .and_then(|cell| cell.phase_layer)
                .is_some_and(|layer| layer.uses_material_density)
    }

    fn admit(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), ImplicitConductionError> {
        checkpoint()?;
        // Every contact is bound to a conductance or named as carried outside.
        if thermal.cells.len() != self.ids.len()
            || self.conduction.links.len() + self.conduction.outside.len() != thermal.links.len()
            || (matches!(self.topology, LinearTopology::Chain)
                && (thermal.links.len() != self.ids.len() - 1 || !thermal.covered_areas.is_empty()))
        {
            return Err(ImplicitConductionError::Invalid(
                "thermal graph bindings changed or direct-chain topology is invalid",
            ));
        }
        self.conduction.admit_ambient_areas(thermal, checkpoint)?;
        for (&id, named) in &self.conduction.outside {
            checkpoint()?;
            named.admit(thermal, id)?;
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
            if self.conduction.current_geometry.is_some()
                && thermal.contact_areas(link) != binding.covered_area_m2
            {
                return Err(ImplicitConductionError::Invalid(
                    "current conduction reference contact coverage changed",
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
            } else if thermal.links.len() != 0
                || cell.boundary.reservoir.is_some()
                || self.conduction.ambient_faces.contains_key(&id)
            {
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
        self.admit_bath_histories(thermal, checkpoint)?;
        checkpoint()?;
        Ok(())
    }

    fn endpoint(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Endpoint, ImplicitConductionError> {
        let mut states = Vec::with_capacity(self.ids.len());
        let mut layer_distances = if self.conduction.has_phase_layers {
            Vec::with_capacity(self.ids.len())
        } else {
            Vec::new()
        };
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
            {
                return Err(DepositionError::MaterialDomain.into());
            }
            states.push(state);
            if self.conduction.has_phase_layers {
                let distance = if let Some(phase) = self.conduction.cells.get(&id)
                    .and_then(|binding| binding.phase_layer)
                {
                    let mut layer = phase.layer;
                    let relative_enthalpy_slope_kg_j = if phase.uses_material_density {
                        layer.density_kg_m3 = material_density(thermal, cell)?;
                        let volume_slope = law
                            .specific_volume_derivative_at_specific_enthalpy(cell.specific_h_j_kg)?
                            .ok_or(ImplicitConductionError::Invalid(
                                "implicit material-density layers require an original specific-volume tangent",
                            ))?;
                        finite(layer.density_kg_m3 * volume_slope)?
                    } else {
                        0.0
                    };
                    Some(ImplicitHalfDistance {
                        value_m: layer.half_distance(cell.mass_kg)?,
                        relative_enthalpy_slope_kg_j,
                    })
                } else {
                    None
                };
                layer_distances.push(distance);
            }
        }
        let mut link_values = BTreeMap::new();
        let mut contacts = BTreeMap::new();
        let mut link_states = BTreeMap::new();
        for (&id, binding) in &self.conduction.links {
            checkpoint()?;
            let ia = self.index[&binding.a];
            let ib = self.index[&binding.b];
            let metric = self.conduction.link_metrics(binding);
            let distances = [
                endpoint_distance(&layer_distances, ia, metric.half_distances_m[0]),
                endpoint_distance(&layer_distances, ib, metric.half_distances_m[1]),
            ];
            let geometry = ConductiveLink {
                half_distances_m: distances.map(|distance| distance.value_m),
                ..binding.geometry
            };
            if let Some(law) = self.link_laws.get(&id) {
                let enthalpies = [
                    thermal.cells[&binding.a].specific_h_j_kg,
                    thermal.cells[&binding.b].specific_h_j_kg,
                ];
                let mut patch = thermal
                    .contact_patch(id)
                    .ok_or(ImplicitConductionError::Invalid(
                        "implicit flux link disappeared",
                    ))?;
                // Constitutive flux uses the current common area; per-side
                // material-reference coverage and original identities stay put.
                patch.link.area_m2 = metric.area_m2;
                let state = if self.has_material_density_geometry(binding.a)
                    || self.has_material_density_geometry(binding.b)
                {
                    require_geometry_tangent(
                        law.state_at_specific_enthalpies_with_geometry_tangent(
                            enthalpies,
                            geometry,
                            patch,
                            distances.map(|distance| distance.relative_enthalpy_slope_kg_j),
                            checkpoint,
                        )?,
                    )?
                } else {
                    law.state_at_specific_enthalpies(enthalpies, geometry, patch, checkpoint)?
                };
                state.admit()?;
                link_values.insert(id, state.conductance_w_k);
                link_states.insert(id, state);
                continue;
            }
            let a = states[ia];
            let b = states[ib];
            let resistance = if let Some(law) = self.contacts.get(&id) {
                let enthalpies = [
                    thermal.cells[&binding.a].specific_h_j_kg,
                    thermal.cells[&binding.b].specific_h_j_kg,
                ];
                let bounds = law.enthalpy_bounds_j_kg();
                if enthalpies
                    .iter()
                    .zip(bounds)
                    .any(|(&h, [lo, hi])| !(lo..=hi).contains(&h))
                {
                    return Err(DepositionError::MaterialDomain.into());
                }
                let state = law.state_at_specific_enthalpies(enthalpies)?;
                state.admit()?;
                contacts.insert(id, state);
                state.resistance_m2_k_w
            } else {
                binding.geometry.contact_resistance_m2_k_w
            };
            link_values.insert(
                id,
                series_contact_conductance_w_k(
                    metric.area_m2,
                    geometry.half_distances_m,
                    [a.conductivity_w_m_k, b.conductivity_w_m_k],
                    resistance,
                )?,
            );
        }
        // Contacts carried outside: the residual owner takes exactly zero in
        // place of their authored conductances and keeps the areas they
        // cover. Written by the binding, as for the explicit timeline.
        self.conduction.carry_outside(thermal, &mut link_values, checkpoint)?;
        let mut boundary_values = BoundaryConductances::default();
        let mut exchange_values = BoundaryExchanges::default();
        let mut reservoir_states = BTreeMap::new();
        let mut ambient_states = BTreeMap::new();
        let mut reservoir_exchange_states = BTreeMap::new();
        let mut ambient_exchange_states = BTreeMap::new();
        for (&id, binding) in &self.conduction.reservoirs {
            checkpoint()?;
            let cell = &thermal.cells[&id];
            let reservoir = cell.boundary.reservoir.ok_or(ImplicitConductionError::Invalid(
                "implicit reservoir disappeared",
            ))?;
            let bath_temperature_k = boundary::temperature_at(
                &self.reservoir_histories, id, reservoir.temperature_k, thermal.time_s,
            )?;
            let index = self.index[&id];
            let metric = self.conduction.reservoir_metrics(binding);
            let distance =
                endpoint_distance(&layer_distances, index, metric.half_distance_m);
            let geometry = ImplicitBoundaryGeometry {
                area_m2: metric.area_m2,
                half_distance_m: distance.value_m,
                surface_resistance_m2_k_w: binding.geometry.surface_resistance_m2_k_w,
            };
            if let Some(law) = self.reservoir_exchanges.get(&id) {
                let state = if self.has_material_density_geometry(id) {
                    require_geometry_tangent(
                        law.state_at_specific_enthalpy_with_geometry_tangent_and_time(
                            cell.specific_h_j_kg,
                            bath_temperature_k,
                            geometry,
                            distance.relative_enthalpy_slope_kg_j,
                            thermal.time_s,
                            checkpoint,
                        )?,
                    )?
                } else {
                    law.state_at_specific_enthalpy_and_time(
                        cell.specific_h_j_kg,
                        bath_temperature_k,
                        geometry,
                        thermal.time_s,
                        checkpoint,
                    )?
                }
                .validate()?;
                exchange_values.reservoir.insert(id, state.baths);
                reservoir_exchange_states.insert(id, state);
                continue;
            }
            if let Some(law) = self.reservoir_laws.get(&id) {
                let state = if self.has_material_density_geometry(id) {
                    require_geometry_tangent(
                        law.state_at_specific_enthalpy_with_geometry_tangent(
                            cell.specific_h_j_kg,
                            bath_temperature_k,
                            geometry,
                            distance.relative_enthalpy_slope_kg_j,
                            checkpoint,
                        )?,
                    )?
                } else {
                    law.state_at_specific_enthalpy(
                        cell.specific_h_j_kg,
                        bath_temperature_k,
                        geometry,
                        checkpoint,
                    )?
                };
                state.admit()?;
                boundary_values
                    .reservoir_w_k
                    .insert(id, state.conductance_w_k);
                reservoir_states.insert(id, state);
                continue;
            }
            let k = states[index].conductivity_w_m_k;
            let g = surface_conductance_w_k(
                metric.area_m2,
                distance.value_m,
                k,
                binding.geometry.surface_resistance_m2_k_w,
            )?;
            boundary_values.reservoir_w_k.insert(id, g);
        }
        for (&id, binding) in &self.conduction.ambient_faces {
            checkpoint()?;
            let cell = &thermal.cells[&id];
            let bath_temperature_k = boundary::temperature_at(
                &self.ambient_histories, id, cell.boundary.ambient_k, thermal.time_s,
            )?;
            let index = self.index[&id];
            let metric = self.conduction.ambient_metrics(binding);
            let distance =
                endpoint_distance(&layer_distances, index, metric.half_distance_m);
            let geometry = ImplicitBoundaryGeometry {
                area_m2: metric.area_m2,
                half_distance_m: distance.value_m,
                surface_resistance_m2_k_w: binding.geometry.surface_resistance_m2_k_w,
            };
            if let Some(law) = self.ambient_exchanges.get(&id) {
                let state = if self.has_material_density_geometry(id) {
                    require_geometry_tangent(
                        law.state_at_specific_enthalpy_with_geometry_tangent_and_time(
                            cell.specific_h_j_kg,
                            bath_temperature_k,
                            geometry,
                            distance.relative_enthalpy_slope_kg_j,
                            thermal.time_s,
                            checkpoint,
                        )?,
                    )?
                } else {
                    law.state_at_specific_enthalpy_and_time(
                        cell.specific_h_j_kg,
                        bath_temperature_k,
                        geometry,
                        thermal.time_s,
                        checkpoint,
                    )?
                }
                .validate()?;
                exchange_values.ambient.insert(id, state.baths);
                ambient_exchange_states.insert(id, state);
                continue;
            }
            if let Some(law) = self.ambient_laws.get(&id) {
                let state = if self.has_material_density_geometry(id) {
                    require_geometry_tangent(
                        law.state_at_specific_enthalpy_with_geometry_tangent(
                            cell.specific_h_j_kg,
                            bath_temperature_k,
                            geometry,
                            distance.relative_enthalpy_slope_kg_j,
                            checkpoint,
                        )?,
                    )?
                } else {
                    law.state_at_specific_enthalpy(
                        cell.specific_h_j_kg,
                        bath_temperature_k,
                        geometry,
                        checkpoint,
                    )?
                };
                state.admit()?;
                boundary_values
                    .ambient_w_k
                    .insert(id, state.conductance_w_k);
                ambient_states.insert(id, state);
                continue;
            }
            let g = surface_conductance_w_k(
                metric.area_m2,
                distance.value_m,
                states[index].conductivity_w_m_k,
                binding.geometry.surface_resistance_m2_k_w,
            )?;
            boundary_values.ambient_w_k.insert(id, g);
        }
        let mut heat_source_states = BTreeMap::new();
        for (&id, law) in &self.heat_sources {
            checkpoint()?;
            let state = law.state_at_specific_enthalpy_and_time(
                thermal.cells[&id].specific_h_j_kg,
                thermal.time_s,
                checkpoint,
            )?;
            state.admit()?;
            heat_source_states.insert(id, state);
        }
        let coupled_heat_source_state = self.coupled_heat_source.as_ref()
            .map(|source| source.endpoint(thermal, checkpoint))
            .transpose()?;
        checkpoint()?;
        let mut endpoint = Endpoint {
            states,
            layer_distances,
            link_values,
            boundary_values,
            exchange_values,
            contacts,
            link_states,
            reservoir_states,
            ambient_states,
            reservoir_exchange_states,
            ambient_exchange_states,
            heat_source_states,
            coupled_heat_source_state,
            transport_tangent: Vec::new(),
            advection_multiplier: 1.0,
            advection_entropy: None,
            stationary_advection: None,
        };
        self.apply_bath_histories(thermal, &mut endpoint, checkpoint)?;
        Ok(endpoint)
    }

    #[cfg(test)]
    fn evaluate<'m>(
        &self,
        step: &BackwardEulerStep<'_, 'm>,
        energy: &[f64],
        sources_w: &BTreeMap<CellId, f64>,
        transfers: &[HeatTransfer],
        entropy: bool,
        controls: ImplicitConductionControls,
        report: &mut ImplicitConductionReport,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<BackwardEulerTrial<'m, Endpoint>, ImplicitConductionError> {
        let mut sources = SourcePowers::default();
        for (&id, &power) in sources_w {
            checkpoint()?;
            sources.powers.entry(id).or_default().add(power)?;
            sources.has_nonzero |= power != 0.0;
        }
        let advection_stage = self.advection.as_ref()
            .map(|advection| advection.stage(step.target_s(), checkpoint))
            .transpose()?;
        self.evaluate_stage(
            step, energy, &sources, transfers, entropy, None,
            advection_stage.as_ref(), controls, report, checkpoint,
        )
    }

    fn evaluate_stage<'m>(
        &self,
        step: &BackwardEulerStep<'_, 'm>,
        energy: &[f64],
        sources: &SourcePowers,
        transfers: &[HeatTransfer],
        entropy: bool,
        quadrature: Option<StageQuadrature<'_>>,
        advection_stage: Option<&ImplicitAdvectionStage<'_>>,
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
            let mut endpoint = self.endpoint(thermal, &mut || (checkpoint.borrow_mut())())?;
            let boundaries = endpoint.boundary_values.clone();
            // The residual owner consumes the bath powers; the analytic
            // Jacobian retains only the already prepared exchange states.
            let exchanges = std::mem::take(&mut endpoint.exchange_values);
            let mut sources = sources.clone();
            for (&id, state) in &endpoint.heat_source_states {
                (checkpoint.borrow_mut())()?;
                sources.add(thermal, id, state.power_w)?;
            }
            if let (Some(source), Some(state)) =
                (&self.coupled_heat_source, &endpoint.coupled_heat_source_state)
            {
                for (&id, &power) in source.cells.iter().zip(&state.powers_w) {
                    (checkpoint.borrow_mut())()?;
                    sources.add(thermal, id, power)?;
                }
            }
            let mut transfers = transfers.to_vec();
            if let Some(advection) = &self.advection {
                let stage = advection_stage.ok_or(DepositionError::Arithmetic)?;
                endpoint.advection_multiplier = stage.multiplier;
                let (transport, flow_sources, tangent, flow_entropy) = advection.endpoint(
                    thermal, stage, &self.index, &endpoint.states, entropy,
                    &mut || (checkpoint.borrow_mut())(),
                )?;
                endpoint.transport_tangent = tangent;
                if let Some(flow_entropy) = flow_entropy {
                    let mut classified = ClassifiedAdvectionEntropy::sample(
                        thermal, &sources, &transfers, flow_entropy,
                        quadrature.map_or(step.duration_s(), |stage| stage.implicit_dt_s),
                        &mut || (checkpoint.borrow_mut())(),
                    )?;
                    // Preserve actual local source ownership even when these
                    // endpoint watts cancel an independent observed source.
                    classified.has_unclassified_sources |= endpoint.heat_source_states
                        .values().any(|state| state.power_w != 0.0)
                        || endpoint.coupled_heat_source_state.as_ref()
                            .is_some_and(|state| state.powers_w.iter().any(|&power| power != 0.0));
                    endpoint.advection_entropy = Some(classified);
                }
                for (id, watts) in transport.heat.bulk.sources_w {
                    (checkpoint.borrow_mut())()?;
                    sources.add(thermal, id, watts)?;
                }
                flow_sources.visit(&mut |id, watts| {
                    (checkpoint.borrow_mut())()?;
                    sources.add(thermal, id, watts)
                })?;
                transfers.extend(transport.transfers);
            }
            Ok::<_, ImplicitConductionError>((BoundaryExchangeStepCoefficients {
                transport: TransportStepCoefficients {
                    heat: BoundaryStepCoefficients {
                        bulk: StepCoefficients {
                            sources_w: BTreeMap::new(),
                            conductances_w_k: endpoint.link_values.clone(),
                            payload: endpoint,
                        },
                        boundaries,
                    },
                    // Observed transfers stay at the accepted start; an explicitly
                    // bound flow graph is evaluated at this actual endpoint.
                    // Both use the original paired debit/credit owner exactly once.
                    transfers,
                },
                exchanges,
            }, sources))
        };
        let mut trial = step.evaluate_exchange_source_option(
            energy, resolve, entropy, quadrature, || (checkpoint.borrow_mut())(),
        )?;
        if let Some(classified) = trial.payload().advection_entropy {
            trial.classify_endpoint_advection(
                classified.advection, classified.source_entropy_j_k,
                classified.has_unclassified_sources, quadrature,
            )?;
        }
        if trial.entropy_audit().is_some() {
            let owned_sources = trial.payload().heat_source_states
                .values().any(|state| state.power_w != 0.0)
                || trial.payload().coupled_heat_source_state.as_ref()
                    .is_some_and(|state| state.powers_w.iter().any(|&power| power != 0.0));
            // This also updates any retained first-stage heat. Applying it
            // after advection classification preserves every original owner.
            trial.preserve_unclassified_sources(owned_sources);
        }
        Ok(trial)
    }

    fn jacobian(
        &self,
        thermal: &DepositionSimulation<'_>,
        endpoint: &Endpoint,
        dt: f64,
        workspace: &mut ImplicitConductionWorkspace,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), ImplicitConductionError> {
        self.assemble_tangent(thermal, endpoint, dt, 1.0, workspace, checkpoint)
    }

    // Shared exact heat-rate tangent. Harmonic analysis requests zero identity
    // shift directly; subtracting I from a time-step Jacobian would erase slow
    // physical relaxation rates through cancellation.
    fn assemble_tangent(
        &self,
        thermal: &DepositionSimulation<'_>,
        endpoint: &Endpoint,
        dt: f64,
        identity: f64,
        workspace: &mut ImplicitConductionWorkspace,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), ImplicitConductionError> {
        workspace.diagonal.fill(identity);
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
                            .and_then(|n| n.checked_add(self.advection.as_ref()
                                .map_or(0, |flow| flow.tangent.len())))
                            .and_then(|n| n.checked_add(endpoint.transport_tangent.len()))
                            .and_then(|n| n.checked_add(self.coupled_heat_source.as_ref()
                                .map_or(0, |source| source.dependencies.len())))
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
            let (qa, qb) = if let Some(state) = endpoint.link_states.get(&id) {
                // The supplied Q leaves A. The original residual uses heat
                // INTO A and energy, not specific enthalpy, coordinates.
                let [a, b] = state.power_enthalpy_derivatives_w_kg_j;
                (
                    finite(-a / thermal.cells[&link.a].mass_kg)?,
                    finite(-b / thermal.cells[&link.b].mass_kg)?,
                )
            } else {
                let a = endpoint.states[ia];
                let b = endpoint.states[ib];
                let g = endpoint.link_values[&id];
                let metric = self.conduction.link_metrics(link);
                let difference = b.temperature_k - a.temperature_k;
                if metric.area_m2 != 0.0
                    && a.conductivity_w_m_k == 0.0
                    && b.conductivity_w_m_k == 0.0
                    && a.conductivity_slope_kg_w_m_k_j != 0.0
                    && b.conductivity_slope_kg_w_m_k_j != 0.0
                    && difference != 0.0
                {
                    // Two independently opening insulators give a
                    // direction-dependent harmonic-mean derivative. At equal
                    // temperatures the complete heat flux is second order;
                    // otherwise a source-owned complete link tangent is needed.
                    return Err(ImplicitConductionError::Invalid(
                        "simultaneous insulating onsets at unequal temperatures require a complete link law",
                    ));
                }
                let mut ga = if b.conductivity_w_m_k == 0.0 {
                    0.0
                } else {
                    conductance_derivative(g, metric.area_m2, endpoint_distance(
                        &endpoint.layer_distances, ia, metric.half_distances_m[0],
                    ), a)?
                };
                let mut gb = if a.conductivity_w_m_k == 0.0 {
                    0.0
                } else {
                    conductance_derivative(g, metric.area_m2, endpoint_distance(
                        &endpoint.layer_distances, ib, metric.half_distances_m[1],
                    ), b)?
                };
                if let Some(contact) = endpoint.contacts.get(&id) {
                    ga = finite(ga + contact.conductance_derivative(g, metric.area_m2, 0)?)?;
                    gb = finite(gb + contact.conductance_derivative(g, metric.area_m2, 1)?)?;
                }
                (
                    finite(
                        (ga * difference - g * a.temperature_slope_kg_k_j)
                            / thermal.cells[&link.a].mass_kg,
                    )?,
                    finite(
                        (gb * difference + g * b.temperature_slope_kg_k_j)
                            / thermal.cells[&link.b].mass_kg,
                    )?,
                )
            };
            add_entry(workspace, ia, ia, -dt * qa)?;
            add_entry(workspace, ia, ib, -dt * qb)?;
            add_entry(workspace, ib, ia, dt * qa)?;
            add_entry(workspace, ib, ib, dt * qb)?;
            let coverage = thermal.contact_areas(&thermal.links[&id]);
            covered[ia] = finite(covered[ia] + coverage[0])?;
            covered[ib] = finite(covered[ib] + coverage[1])?;
        }
        for (&id, named) in &self.conduction.outside {
            checkpoint()?;
            // No conduction entry; the contact still covers both faces, as
            // in the residual owner's exposed areas.
            let coverage = thermal.contact_areas(&thermal.links[&id]);
            let (ia, ib) = (self.index[&named.a], self.index[&named.b]);
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
            let mut derivative = if let Some(boundary) = endpoint.ambient_exchange_states.get(&id) {
                -boundary.power_enthalpy_derivative_w_kg_j
            } else if let Some(boundary) = endpoint.ambient_states.get(&id) {
                -boundary.power_enthalpy_derivative_w_kg_j
            } else if let Some(binding) = self.conduction.ambient_faces.get(&id) {
                let g = endpoint.boundary_values.ambient_w_k[&id];
                let metric = self.conduction.ambient_metrics(binding);
                let dg = conductance_derivative(
                    g,
                    metric.area_m2,
                    endpoint_distance(&endpoint.layer_distances, i, metric.half_distance_m),
                    state,
                )?;
                dg * (cell.boundary.ambient_k - state.temperature_k)
                    - g * state.temperature_slope_kg_k_j
            } else {
                -(area * cell.boundary.htc_w_m2_k) * state.temperature_slope_kg_k_j
            };
            if let Some(reservoir) = cell.boundary.reservoir {
                if let Some(boundary) = endpoint.reservoir_exchange_states.get(&id) {
                    derivative -= boundary.power_enthalpy_derivative_w_kg_j;
                } else if let Some(boundary) = endpoint.reservoir_states.get(&id) {
                    derivative -= boundary.power_enthalpy_derivative_w_kg_j;
                } else {
                    let binding = self.conduction.reservoirs[&id];
                    let g = endpoint.boundary_values.reservoir_w_k[&id];
                    let metric = self.conduction.reservoir_metrics(&binding);
                    let dg = conductance_derivative(
                        g,
                        metric.area_m2,
                        endpoint_distance(&endpoint.layer_distances, i, metric.half_distance_m),
                        state,
                    )?;
                    derivative += dg * (reservoir.temperature_k - state.temperature_k)
                        - g * state.temperature_slope_kg_k_j;
                }
            }
            if let Some(source) = endpoint.heat_source_states.get(&id) {
                derivative = finite(derivative + source.power_enthalpy_derivative_w_kg_j)?;
            }
            add_entry(workspace, i, i, -dt * derivative / cell.mass_kg)?;
        }
        if let (Some(source), Some(state)) =
            (&self.coupled_heat_source, &endpoint.coupled_heat_source_state)
        {
            for (&dependency, &derivative) in source.dependencies.iter()
                .zip(&state.power_enthalpy_derivatives_w_kg_j)
            {
                checkpoint()?;
                add_entry(
                    workspace, dependency.row, dependency.column,
                    dependency.energy_jacobian_entry(dt, derivative)?,
                )?;
            }
        }
        if let Some(advection) = &self.advection {
            for &(row, col, derivative) in &advection.tangent {
                checkpoint()?;
                let derivative = advection::scaled_tangent(derivative, endpoint.advection_multiplier)?;
                add_entry(workspace, row, col, -dt * derivative)?;
            }
            // Liquid/reconstruction entries were evaluated from this same
            // stage's scaled rates and must not be multiplied a second time.
            for &(row, col, derivative) in &endpoint.transport_tangent {
                checkpoint()?;
                add_entry(workspace, row, col, -dt * derivative)?;
            }
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
        self.solve_step_with_sources(
            thermal,
            target_s,
            &SourcePowers::default(),
            &[],
            controls,
            workspace,
            report,
            checkpoint,
        )
    }

    // The same bounded linear solve serves transient Newton and stationary
    // heat balance; only their physical residual and identity shift differ.
    fn solve_direction(
        &self,
        workspace: &mut ImplicitConductionWorkspace,
        report: &mut ImplicitConductionReport,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), ImplicitConductionError> {
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
                        Err(error) => { cancellation_error = Some(error); true }
                    },
                );
                if let Some(error) = cancellation_error { return Err(error.into()); }
                solved?;
            }
            LinearTopology::Graph(linear) => graph::solve(workspace, linear, report, checkpoint)?,
        }
        Ok(())
    }

    fn solve_step_with_sources<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        sources: &SourcePowers,
        transfers: &[HeatTransfer],
        controls: ImplicitConductionControls,
        workspace: &mut ImplicitConductionWorkspace,
        report: &mut ImplicitConductionReport,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'m>, EnergyAudit, EntropyAudit, f64), ImplicitConductionError>
    {
        let trial = self.solve_stage(
            thermal, target_s, sources, transfers, None, controls, workspace, report, checkpoint,
        )?;
        accept_trial(trial, controls, checkpoint)
    }

    // Shared nonlinear owner. The optional stage weight changes only the
    // native residual's rate quadrature and its analytic Jacobian; all trial
    // enthalpies remain relative to this actual physical initial state.
    fn solve_stage<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        sources: &SourcePowers,
        transfers: &[HeatTransfer],
        quadrature: Option<StageQuadrature<'_>>,
        controls: ImplicitConductionControls,
        workspace: &mut ImplicitConductionWorkspace,
        report: &mut ImplicitConductionReport,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<BackwardEulerTrial<'m, Endpoint>, ImplicitConductionError> {
        let step = BackwardEulerStep::new(thermal, target_s, &mut *checkpoint)?;
        let advection_stage = self.advection.as_ref()
            .map(|advection| advection.stage(target_s, &mut *checkpoint))
            .transpose()?;
        let implicit_dt_s = quadrature.map_or(step.duration_s(), |stage| stage.implicit_dt_s);
        workspace.energy.fill(0.0);
        let mut trial = self.evaluate_stage(
            &step,
            &workspace.energy,
            sources,
            transfers,
            false,
            quadrature,
            advection_stage.as_ref(),
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
                return self.evaluate_stage(
                    &step,
                    &workspace.energy,
                    sources,
                    transfers,
                    true,
                    quadrature,
                    advection_stage.as_ref(),
                    controls,
                    report,
                    checkpoint,
                );
            }
            if maximum <= controls.cell_tolerance_j {
                let floor_j = self.audit_floor_j(thermal, &workspace.energy, trial.residual_j())?;
                if controls.total_tolerance_j < floor_j {
                    return Err(ImplicitConductionError::AllowanceBelowAuditFloor {
                        allowance_j: controls.total_tolerance_j,
                        floor_j,
                    });
                }
            }
            if iteration == controls.maximum_newton_iterations {
                return Err(ImplicitConductionError::NotConverged {
                    maximum_residual_j: maximum,
                });
            }
            self.jacobian(
                thermal,
                trial.payload(),
                implicit_dt_s,
                workspace,
                checkpoint,
            )?;
            for (rhs, &residual) in workspace.rhs.iter_mut().zip(trial.residual_j()) {
                *rhs = -residual;
            }
            self.solve_direction(workspace, report, checkpoint)?;
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
                    value = resolvable_increment_j(cell.specific_h_j_kg, value, cell.mass_kg);
                    representable |= value != workspace.energy[i];
                    workspace.proposal[i] = value;
                }
                if admitted && representable {
                    match self.evaluate_stage(
                        &step,
                        &workspace.proposal,
                        sources,
                        transfers,
                        false,
                        quadrature,
                        advection_stage.as_ref(),
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

    // The rounding floor of a trial's energy audit. A cell's stored joules move
    // on its enthalpy grid, in steps of m * ulp(h), so without cancellation
    // between cells the audit resolves no finer than half a step per cell. The
    // audit is a compensated sum of at most K = 2n + 3 rounded terms (storage and
    // heat per cell, and the external account's three parts), whose rounding is
    // at most gamma_K = K u / (1 - K u) times their magnitudes, u = 2^-53.
    // Both parts come from this trial; neither is a chosen scale.
    fn audit_floor_j(
        &self,
        thermal: &DepositionSimulation<'_>,
        energy_j: &[f64],
        residual_j: &[f64],
    ) -> Result<f64, DepositionError> {
        let unit = f64::EPSILON / 2.0;
        let terms = (2 * self.ids.len() + 3) as f64;
        let gamma = finite(terms * unit / (1.0 - terms * unit))?;
        let (mut grid, mut magnitude) = (Sum::default(), Sum::default());
        for ((&id, &stored), &residual) in self.ids.iter().zip(energy_j).zip(residual_j) {
            let cell = &thermal.cells[&id];
            let h = finite(cell.specific_h_j_kg + stored / cell.mass_kg)?.abs();
            grid.add(0.5 * cell.mass_kg * (h.next_up() - h))?;
            magnitude.add(stored.abs() + (stored - residual).abs())?;
        }
        finite(grid.total()? + gamma * magnitude.total()?)
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
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'m>, ImplicitConductionReport), ImplicitConductionError> {
        let (next, (), report) = self
            .advance_observed_to(thermal, &(), target_s, controls, workspace, checkpoint)
            .map_err(|error| match error {
                ImplicitConductionObservedError::Implicit(error) => error,
                ImplicitConductionObservedError::Observer(error) => error.into(),
            })?;
        Ok((next, report))
    }
}

fn accept_trial<'m>(
    trial: BackwardEulerTrial<'m, Endpoint>,
    controls: ImplicitConductionControls,
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<(DepositionSimulation<'m>, EnergyAudit, EntropyAudit, f64), ImplicitConductionError> {
    let maximum = trial.maximum_cell_residual_j();
    let entropy = *trial
        .entropy_audit()
        .ok_or(ImplicitConductionError::Invalid("missing final entropy account"))?;
    let (next, _, audit) = trial.accept(
        controls.cell_tolerance_j,
        controls.total_tolerance_j,
        checkpoint,
    )?;
    Ok((next, audit, entropy, maximum))
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
    distance: ImplicitHalfDistance,
    state: ImplicitConductionState,
) -> Result<f64, DepositionError> {
    if area == 0.0 {
        return Ok(0.0);
    }
    if state.conductivity_w_m_k == 0.0 {
        if state.conductivity_slope_kg_w_m_k_j == 0.0 {
            return Ok(0.0);
        }
        // G=A*k/(d+k*R_other): at k=0 its admitted one-sided derivative
        // is A*k'/d, evaluated at the current density-driven thickness.
        // The d' term vanishes at k=0. For a link the caller first excludes
        // an insulating opposite side.
        return finite((area / distance.value_m) * state.conductivity_slope_kg_w_m_k_j);
    }
    if g == 0.0
        || (state.conductivity_slope_kg_w_m_k_j == 0.0
            && distance.relative_enthalpy_slope_kg_j == 0.0)
    {
        return Ok(0.0);
    }
    // G*(bulk resistance / total resistance)*(k'/k - d'/d), keeping the geometric
    // resistance ratio bounded rather than directly forming G²/k².
    finite(
        g * finite((g / area) * (distance.value_m / state.conductivity_w_m_k))?
            * finite(state.conductivity_slope_kg_w_m_k_j / state.conductivity_w_m_k
                - distance.relative_enthalpy_slope_kg_j)?,
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
    total.accumulate(next)
}
#[derive(Clone, Copy, Default)]
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

#[cfg(test)]
#[path = "implicit/density_tests.rs"]
mod density_tests;
