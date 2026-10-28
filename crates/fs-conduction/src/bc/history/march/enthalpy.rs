//! Scheduled melting/solidification through the existing reference-mass solver.
//!
//! Specific enthalpy, not temperature, is the persistent state. Latent heat is
//! stored once in each original phase chart; neither apparent heat capacity nor
//! equilibrium-density changes are substituted into the reference masses.
//! Boundary/source histories use the same explicit endpoint quadrature as the
//! temperature runner. Prescribed temperatures use each node's original chart;
//! ambiguous plateau temperatures and out-of-domain values remain refused.
//!
//! This is stationary P1 transport, not melt flow, nucleation, expansion, moving
//! interfaces or remeshing. Constitutive charts, conductivity assignments,
//! phase multipliers and matching-face contact remain with their existing owners.

/// Radiative reservoir histories driving the same persistent enthalpy state.
pub mod radiation;

/// Local time-error control in temperature, enthalpy and liquid fraction.
pub mod adaptive;

/// Material-backed thermal contact feedback retained through scheduled evolution.
pub mod contact;

/// Moving absorbed-power footprints with persistent phase and optional radiation.
pub mod scanning;

use super::{ThermalHistoryProblem, poll};
use super::events;
use crate::{ConductionError, ConductionProblem, ScalarField, ThermalBc, ThermalBoundary};
use crate::material::LiquidMassFractionConductivity;
use crate::transient::enthalpy::{EnthalpyBudget, EnthalpyError, EnthalpyStepConfig, EnthalpyStepSolution};
use crate::transient::enthalpy::heterogeneous::{
    HeterogeneousEnthalpyBackwardEuler, HeterogeneousEnthalpyError, ReferenceEnthalpyMaterial,
};
use fs_exec::Cx;

/// Explicit constitutive identities and frozen reference mass for this mesh.
/// Conductivity still uses `ThermalHistoryProblem::element_materials`; the two
/// assignments are never guessed from one another or averaged at shared nodes.
pub struct EnthalpyHistoryMaterials<'s, 'c> {
    /// Original phase charts and separately supplied reference densities.
    pub table: &'s [ReferenceEnthalpyMaterial<'c>],
    /// One index into `table` per tetrahedron. Shared vertices must have one ID.
    pub of_element: &'s [usize],
    /// Optional multiplier on each original conductivity tensor, evaluated at
    /// the mean nodal liquid mass fraction on EVERY nonlinear candidate.
    pub phase_conductivity: Option<Vec<LiquidMassFractionConductivity>>,
}

/// Explicit spatial admission and total lifetime physical-step allowance.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EnthalpyHistoryBudget {
    /// Original enthalpy storage/assembly caps, not a substitute solver budget.
    pub spatial: EnthalpyBudget,
    /// Total accepted endpoints; retries never consume this allowance.
    pub max_steps: usize,
}

/// Resumable scenario trajectory retaining the full caloric state of each node.
///
/// Construction supports one or several materials through the same checked
/// heterogeneous storage owner. The constructor copies assignments and retains
/// the original chart borrows; the input table need not outlive construction.
/// Only read-only storage inspection is exposed. A step cannot replace material
/// identities, reset a liquid fraction, or omit the selected phase conductivity.
pub struct EnthalpyHistoryRun<'h, 'm, 'c> {
    problem: ThermalHistoryProblem<'h, 'm>,
    storage: HeterogeneousEnthalpyBackwardEuler<'m, 'c>,
    specific_h: Vec<f64>,
    temperature: Vec<f64>,
    liquid_mass_fraction: Vec<f64>,
    time_s: f64,
    accepted_steps: usize,
    max_steps: usize,
}

impl<'h, 'm, 'c> EnthalpyHistoryRun<'h, 'm, 'c> {
    /// Bind the original mesh, reference materials, schedules and initial J/kg.
    ///
    /// Initial temperature and liquid fraction are resolved from the supplied
    /// enthalpy at each node, never inverted from an ambiguous temperature.
    /// No physical step or boundary-storage reset occurs during admission.
    ///
    /// # Errors
    /// Original material/mesh/budget/chart refusals, malformed source or boundary
    /// data, zero step allowance, cancellation, or a prescribed temperature with
    /// no unique admissible enthalpy in its original material chart.
    pub fn new(
        cx: &Cx<'_>,
        problem: ThermalHistoryProblem<'h, 'm>,
        materials: EnthalpyHistoryMaterials<'_, 'c>,
        initial_specific_h: &[f64],
        budget: EnthalpyHistoryBudget,
    ) -> Result<Self, HeterogeneousEnthalpyError> {
        poll(cx, 0).map_err(conduction)?;
        if budget.max_steps == 0 {
            return Err(invalid("enthalpy history needs a positive accepted-step budget"));
        }
        let mesh = problem.boundary.mesh();
        if initial_specific_h.len() != mesh.vertex_count() {
            return Err(invalid("initial specific enthalpy must contain one value per mesh vertex"));
        }
        let mut storage = HeterogeneousEnthalpyBackwardEuler::new(
            cx, mesh, materials.table, materials.of_element, budget.spatial,
        )?;
        if let Some(laws) = materials.phase_conductivity {
            storage = storage.with_phase_conductivity(cx, laws)?;
        }
        if let Some(assigned) = problem.element_materials {
            assigned.validate_for(mesh).map_err(conduction)?;
        }
        let time_s = problem.boundary.initial_time_s();
        let boundary = problem.boundary.at(cx, time_s).map_err(conduction)?;
        storage.admit_prescribed_boundary(cx, &boundary)?;
        problem.source.at(time_s, mesh.vertex_count()).map_err(conduction)?;
        let mut temperature = Vec::with_capacity(initial_specific_h.len());
        let mut liquid_mass_fraction = Vec::with_capacity(initial_specific_h.len());
        for (vertex, &h) in initial_specific_h.iter().enumerate() {
            if vertex % 512 == 0 { poll(cx, vertex).map_err(conduction)?; }
            let curve = storage.phase_curve_at(vertex)
                .ok_or_else(|| invalid("reference chart is missing for an admitted vertex"))?;
            let state = curve.state_at_specific_enthalpy(h)
                .map_err(|source| EnthalpyError::Phase { vertex, source })?;
            temperature.push(state.temperature_k());
            liquid_mass_fraction.push(state.liquid_mass_fraction());
        }
        let specific_h = initial_specific_h.to_vec();
        poll(cx, initial_specific_h.len()).map_err(conduction)?;
        Ok(Self { problem, storage, specific_h, temperature, liquid_mass_fraction,
            time_s, accepted_steps: 0, max_steps: budget.max_steps })
    }

    /// Accepted nodal specific enthalpies [J/kg], including latent storage.
    #[must_use]
    pub fn specific_enthalpy_j_kg(&self) -> &[f64] { &self.specific_h }

    /// Chart-resolved accepted nodal temperatures [K].
    #[must_use]
    pub fn temperature(&self) -> &[f64] { &self.temperature }

    /// Chart-resolved accepted liquid mass fractions; no separate phase update.
    #[must_use]
    pub fn liquid_mass_fraction(&self) -> &[f64] { &self.liquid_mass_fraction }

    /// Original masses, material assignment and per-node chart inspection.
    #[must_use]
    pub fn storage(&self) -> &HeterogeneousEnthalpyBackwardEuler<'m, 'c> { &self.storage }

    /// Physical time of the accepted state [s].
    #[must_use]
    pub const fn time_s(&self) -> f64 { self.time_s }

    /// Successfully published endpoints, excluding failed attempts.
    #[must_use]
    pub const fn accepted_steps(&self) -> usize { self.accepted_steps }

    /// Lifetime accepted-step allowance, retained on resume.
    #[must_use]
    pub const fn max_steps(&self) -> usize { self.max_steps }

    /// Increase the allowance without changing any caloric state or schedule.
    ///
    /// # Errors
    /// A new limit not strictly larger than the current limit.
    pub fn extend_step_budget(&mut self, max_steps: usize) -> Result<(), HeterogeneousEnthalpyError> {
        if max_steps <= self.max_steps { return Err(invalid("enthalpy step budget must increase")); }
        self.max_steps = max_steps;
        Ok(())
    }

    /// Solve one explicit endpoint with the original nonlinear enthalpy owner.
    ///
    /// Natural boundaries, generation, k(T), phase-dependent tensor multipliers
    /// and contact use the SAME old enthalpy on every trial. The production
    /// response retains Newton work, phase fields and the actual joule balance.
    /// A failed endpoint changes neither state nor time nor the accepted budget.
    ///
    /// # Errors
    /// Budget/time/history/domain/cancellation or original nonlinear/energy
    /// refusal. Caller-owned time grids are not automatically event-split.
    pub fn advance_to(&mut self, cx: &Cx<'_>, time_s: f64, config: EnthalpyStepConfig)
        -> Result<EnthalpyStepSolution, HeterogeneousEnthalpyError>
    {
        let (dt_s, boundary, source) = self.prepare_endpoint(cx, time_s)?;
        let solved = self.storage.advance(cx, ConductionProblem {
            mesh: self.problem.boundary.mesh(), boundary: &boundary,
            material: self.problem.material, element_materials: self.problem.element_materials,
            source: &source,
        }, self.problem.interfaces, &self.specific_h, dt_s, config)?;
        self.publish(cx, time_s, &solved)?;
        Ok(solved)
    }

    /// Advance ONE interval ending at the next load/boundary knot or step cap.
    /// Held data use their exact interval-end left limits, so short heat pulses
    /// change latent enthalpy even when their switch-on and switch-off both lie
    /// before the requested target. Repeat until `time_s() == target_s`.
    /// The original Newton, phase, energy and natural-boundary checks remain.
    ///
    /// # Errors
    /// Invalid time/cap, exhausted budget, cancellation or any original solve
    /// refusal. All caloric fields and time remain unchanged on a failed call.
    pub fn advance_event_step(&mut self, cx: &Cx<'_>, target_s: f64,
        max_step_s: f64, config: EnthalpyStepConfig)
        -> Result<events::HistoryEventStep<EnthalpyStepSolution>, HeterogeneousEnthalpyError>
    {
        poll(cx, self.accepted_steps).map_err(conduction)?;
        if self.accepted_steps >= self.max_steps {
            return Err(invalid("enthalpy history step budget exhausted"));
        }
        let start_time_s = self.time_s;
        let end_time_s = events::next_endpoint(cx, self.problem.boundary, &self.problem.source,
            std::iter::empty(), start_time_s, target_s, max_step_s).map_err(conduction)?;
        let boundary = events::boundary_before(cx, self.problem.boundary, end_time_s).map_err(conduction)?;
        self.storage.admit_prescribed_boundary(cx, &boundary)?;
        let source = events::source_before(&self.problem.source, end_time_s, self.specific_h.len())
            .map_err(conduction)?;
        let report = self.storage.advance(cx, ConductionProblem {
            mesh: self.problem.boundary.mesh(), boundary: &boundary,
            material: self.problem.material, element_materials: self.problem.element_materials,
            source: &source,
        }, self.problem.interfaces, &self.specific_h, end_time_s-start_time_s, config)?;
        self.publish(cx, end_time_s, &report)?;
        Ok(events::HistoryEventStep { start_time_s, end_time_s, report })
    }

    fn prepare_endpoint(&self, cx: &Cx<'_>, time_s: f64)
        -> Result<(f64, ThermalBoundary, ScalarField), HeterogeneousEnthalpyError>
    {
        poll(cx, self.accepted_steps).map_err(conduction)?;
        if self.accepted_steps >= self.max_steps {
            return Err(invalid("enthalpy history step budget exhausted"));
        }
        let dt_s = time_s - self.time_s;
        if !time_s.is_finite() || !dt_s.is_finite() || dt_s <= 0.0 {
            return Err(invalid("enthalpy endpoint must be strictly later with a finite positive interval"));
        }
        let boundary = self.problem.boundary.at(cx, time_s).map_err(conduction)?;
        self.storage.admit_prescribed_boundary(cx, &boundary)?;
        let source = self.problem.source.at(time_s, self.problem.boundary.mesh().vertex_count())
            .map_err(conduction)?;
        poll(cx, self.accepted_steps).map_err(conduction)?;
        Ok((dt_s, boundary, source))
    }

    fn publish(&mut self, cx: &Cx<'_>, time_s: f64, solved: &EnthalpyStepSolution)
        -> Result<(), HeterogeneousEnthalpyError>
    {
        // Prepare every owned array before changing any accepted field.
        let h = solved.specific_enthalpy_j_kg.clone();
        let temperature = solved.temperature.clone();
        let fraction = solved.liquid_mass_fraction.clone();
        poll(cx, self.accepted_steps).map_err(conduction)?;
        self.specific_h = h;
        self.temperature = temperature;
        self.liquid_mass_fraction = fraction;
        self.time_s = time_s;
        self.accepted_steps += 1;
        Ok(())
    }
}

fn natural_only(boundary: &ThermalBoundary) -> Result<(), HeterogeneousEnthalpyError> {
    if boundary.conditions().iter().any(|row| matches!(row, ThermalBc::Dirichlet { .. })) {
        return Err(EnthalpyError::UnsupportedDirichlet.into());
    }
    Ok(())
}

fn conduction(error: ConductionError) -> HeterogeneousEnthalpyError { EnthalpyError::from(error).into() }
fn invalid(what: &'static str) -> HeterogeneousEnthalpyError { EnthalpyError::InvalidInput(what).into() }
