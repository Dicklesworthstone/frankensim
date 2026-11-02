//! Scheduled black reservoirs with original material emissivity on P1 traces.
//!
//! Each interval uses the existing material radiation and primary-enthalpy
//! owners. Their full h feedback, pointwise epsilon*T^4 integration, physical
//! energy gates and source-domain refusals remain authoritative. Histories own
//! only time sampling and atomic publication; no phase or heat solver is added.

use super::super::{
    EnthalpyHistoryRun, EnthalpyStepConfig, EnthalpyStepSolution,
    HeterogeneousEnthalpyError, conduction, invalid, poll,
    adaptive::{PhaseErrorTolerance, phase_error},
};
use super::super::super::events::{
    self, HistoryEventStep,
    adaptive::{AdaptiveHistoryConfig, AdaptiveHistoryStep, admit_proposal, select_step},
};
use crate::{
    ConductionError, ConductionProblem, MaterialAmbientRadiationPatch,
    MaterialPointwiseAmbientRadiation, SURFACE_EMISSIVITY_PROPERTY, ScalarField,
    SurfaceEmissivity, TEMPERATURE_DIMS, ThermalBoundary,
};
use fs_exec::Cx;
use fs_material::{
    phase::EquilibriumEnthalpyPhaseCurve,
    state_point::MaterialPropertySelection,
    thermal::source::ThermalMaterialPoint,
};
use fs_scenario::signal::TimeSignal;

/// One original material trace and its independent black-reservoir schedule.
/// Different caloric objects require distinct traces, exactly as in the direct
/// material radiation owner. Temperature and liquid fraction evolve from h;
/// pressure, finish and other material conditions remain explicitly fixed.
#[derive(Debug, Clone)]
pub struct MaterialRadiationRegionHistory<'a> {
    pub region: String,
    /// Exact original caloric object used by every participating storage node.
    pub curve: &'a EquilibriumEnthalpyPhaseCurve,
    pub material: ThermalMaterialPoint<'a>,
    /// Resolved at construction; the selected claim remains pinned on resume.
    pub selection: MaterialPropertySelection,
    /// Absolute reservoir temperature in kelvin, independent of convection.
    pub ambient_temperature: TimeSignal,
}

/// Persistent h/T/f trajectory with material pointwise ambient radiation.
///
/// The original material declarations and selected source claims remain attached
/// across every endpoint, retry and resume. Radiation adds to the original
/// natural boundary; its signed heat remains separate from convection/source
/// heat. The reservoir has view factor one, with no enclosure or spectral model.
pub struct MaterialRadiativeEnthalpyHistoryRun<'h, 'm, 'c> {
    run: EnthalpyHistoryRun<'h, 'm, 'c>,
    patches: Vec<MaterialRadiationRegionHistory<'c>>,
    reference_emissivities: Vec<SurfaceEmissivity>,
}

impl<'h, 'm, 'c> MaterialRadiativeEnthalpyHistoryRun<'h, 'm, 'c> {
    /// Attach original material sources before the first accepted endpoint.
    ///
    /// The original radiation owner admits geometry, source interpolation and
    /// the initial nodal h states. Exact chart-object identity is checked against
    /// storage, so another chart with the same card or temperature cannot replace
    /// its latent state. Construction takes no physical step.
    ///
    /// # Errors
    /// Advanced trajectory, malformed reservoir/region, mismatched chart, source
    /// support/interpolation refusal or cancellation. No input state is advanced.
    pub fn new(
        cx: &Cx<'_>, run: EnthalpyHistoryRun<'h, 'm, 'c>,
        patches: Vec<MaterialRadiationRegionHistory<'c>>,
    ) -> Result<Self, HeterogeneousEnthalpyError> {
        poll(cx, 0).map_err(conduction)?;
        if run.accepted_steps() != 0 {
            return Err(invalid("material radiation must be attached before the first enthalpy step"));
        }
        let mut reference_emissivities = Vec::with_capacity(patches.len());
        for (index, patch) in patches.iter().enumerate() {
            poll(cx, index).map_err(conduction)?;
            let (emissivity, _) = SurfaceEmissivity::from_condition(
                &patch.region, patch.material, patch.selection.clone(),
            ).map_err(conduction)?;
            reference_emissivities.push(emissivity);
        }
        let prepared = Self { run, patches, reference_emissivities };
        let boundary = prepared.run.problem.boundary.at(cx, prepared.run.time_s())
            .map_err(conduction)?;
        let patches = prepared.sample(cx, prepared.run.time_s(), false).map_err(conduction)?;
        let radiation = MaterialPointwiseAmbientRadiation::new(
            cx, prepared.run.problem.boundary.mesh(), &boundary, &patches,
        ).map_err(conduction)?;
        if !radiation.matches_curves(cx, |vertex| prepared.run.storage.phase_curve_at(vertex)
            .expect("admitted history vertex has an original caloric chart"))
            .map_err(conduction)?
        {
            return Err(invalid("material radiation history requires the exact storage caloric objects"));
        }
        radiation.load_and_report(cx, &prepared.run.specific_h, &prepared.run.temperature)
            .map_err(conduction)?;
        poll(cx, prepared.patches.len()).map_err(conduction)?;
        Ok(prepared)
    }

    /// Accepted enthalpy, temperature, liquid fractions, material state and time.
    #[must_use]
    pub fn conduction(&self) -> &EnthalpyHistoryRun<'h, 'm, 'c> { &self.run }

    /// Original immutable material bindings, selections and reservoir schedules.
    #[must_use]
    pub fn patches(&self) -> &[MaterialRadiationRegionHistory<'c>] { &self.patches }

    /// Preparation-time answers and pinned claims, in patch order. Actual nodal
    /// emissivities remain in each accepted step's pointwise radiation report.
    #[must_use]
    pub fn reference_emissivities(&self) -> &[SurfaceEmissivity] {
        &self.reference_emissivities
    }

    /// Extend the lifetime endpoint allowance without dropping material sources.
    pub fn extend_step_budget(&mut self, maximum: usize) -> Result<(), HeterogeneousEnthalpyError> {
        self.run.extend_step_budget(maximum)
    }

    /// Advance one explicit endpoint using right-continuous schedule values.
    ///
    /// Original source emissivity is refreshed at every Newton h candidate,
    /// including on a constant-temperature latent plateau. The returned native
    /// response keeps pointwise radiation, convection and generation separate.
    /// Caller-owned time grids are not event-split by this method.
    ///
    /// # Errors
    /// Original time/budget/history/source/solver/energy/cancellation refusal.
    /// Every refusal preserves all accepted caloric fields and the clock.
    pub fn advance_to(&mut self, cx: &Cx<'_>, time_s: f64, config: EnthalpyStepConfig)
        -> Result<EnthalpyStepSolution, HeterogeneousEnthalpyError>
    {
        let (dt_s, boundary, source) = self.run.prepare_endpoint(cx, time_s)?;
        let solved = self.solve(cx, &self.run.specific_h, &boundary, &source,
            dt_s, time_s, false, config)?;
        self.run.publish(cx, time_s, &solved)?;
        Ok(solved)
    }

    /// Advance one interval clipped to source, boundary and reservoir events.
    /// Held schedules use their exact interval-end left limits: a bath switching
    /// off at the endpoint cannot erase heat transferred during that interval.
    /// Repeat until the accepted time reaches the requested target.
    ///
    /// # Errors
    /// Invalid target/cap or an original physical/budget/cancellation refusal.
    /// A failed interval changes none of the accepted state or source bindings.
    pub fn advance_event_step(&mut self, cx: &Cx<'_>, target_s: f64,
        max_step_s: f64, config: EnthalpyStepConfig)
        -> Result<HistoryEventStep<EnthalpyStepSolution>, HeterogeneousEnthalpyError>
    {
        poll(cx, self.run.accepted_steps()).map_err(conduction)?;
        if self.run.accepted_steps() >= self.run.max_steps() {
            return Err(invalid("enthalpy history step budget exhausted"));
        }
        let start_time_s = self.run.time_s();
        let end_time_s = events::next_endpoint(cx, self.run.problem.boundary,
            &self.run.problem.source, self.patches.iter().map(|p| &p.ambient_temperature),
            start_time_s, target_s, max_step_s).map_err(conduction)?;
        let report = self.solve_before(cx, &self.run.specific_h, start_time_s, end_time_s, config)?;
        self.run.publish(cx, end_time_s, &report)?;
        Ok(HistoryEventStep { start_time_s, end_time_s, report })
    }

    /// Adapt the same scheduled material response in h, T and liquid fraction.
    ///
    /// Coarse and half-step proposals share immutable accepted history. The
    /// second half uses the first half's complete h state; its source emissivity
    /// is never frozen to the coarse or initial phase. All three error gates and
    /// original physical gates must pass. Both accepted reports are returned,
    /// and the two endpoint slots publish atomically through the original run.
    /// This estimates local time error, not a global or continuum error bound.
    ///
    /// # Errors
    /// Original solver/source refusal, invalid policy, insufficient two-step
    /// budget, cancellation or an unmet error target at the floor/attempt cap.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_adaptive_step(&mut self, cx: &Cx<'_>, target_s: f64, proposed_step_s: f64,
        step: EnthalpyStepConfig, adaptive: AdaptiveHistoryConfig, phase: PhaseErrorTolerance)
        -> Result<AdaptiveHistoryStep<EnthalpyStepSolution>, HeterogeneousEnthalpyError>
    {
        poll(cx, self.run.accepted_steps()).map_err(conduction)?;
        phase.validate()?;
        let start = self.run.time_s();
        let cap = admit_proposal(cx, self.run.accepted_steps(), self.run.max_steps(),
            proposed_step_s, adaptive).map_err(conduction)?;
        let end = events::next_endpoint(cx, self.run.problem.boundary, &self.run.problem.source,
            self.patches.iter().map(|p| &p.ambient_temperature), start, target_s, cap)
            .map_err(conduction)?;
        let selected = select_step(cx, start, end, adaptive, conduction, |middle, end| {
            let coarse = self.solve_before(cx, &self.run.specific_h, start, end, step)?;
            let first = self.solve_before(cx, &self.run.specific_h, start, middle, step)?;
            let second = self.solve_before(cx, &first.specific_enthalpy_j_kg, middle, end, step)?;
            let error = phase_error(cx, &self.run, &coarse, &second, adaptive, phase)?;
            Ok((first, second, error))
        })?;
        self.run.publish(cx, selected.end_time_s, &selected.second_half)?;
        self.run.accepted_steps += 1; // The other reserved half step, without partial publication.
        Ok(selected)
    }

    fn sample(&self, cx: &Cx<'_>, time_s: f64, before: bool)
        -> Result<Vec<MaterialAmbientRadiationPatch<'_>>, ConductionError>
    {
        let mut sampled = Vec::with_capacity(self.patches.len());
        for (index, (patch, reference)) in self.patches.iter()
            .zip(&self.reference_emissivities).enumerate()
        {
            poll(cx, index)?;
            let value = if before {
                events::value_before(&patch.ambient_temperature, time_s)?
            } else {
                patch.ambient_temperature.eval(time_s).map_err(|error|
                    ConductionError::Radiation { surface: patch.region.clone(),
                        what: format!("radiative reservoir history: {error}"),
                        fix: "supply a valid absolute-temperature history".into() })?
            };
            if value.dims != TEMPERATURE_DIMS {
                return Err(ConductionError::Dimensions {
                    context: format!("material radiative reservoir on {:?}", patch.region),
                    expected: TEMPERATURE_DIMS.0, found: value.dims.0,
                });
            }
            sampled.push(MaterialAmbientRadiationPatch {
                region: &patch.region, ambient_temperature_k: value.value,
                curve: patch.curve, material: patch.material,
                selection: MaterialPropertySelection::PinnedByProperty(vec![(
                    SURFACE_EMISSIVITY_PROPERTY.into(), reference.receipt().selected,
                )]),
            });
        }
        Ok(sampled)
    }

    fn solve_before(&self, cx: &Cx<'_>, old: &[f64], from: f64, to: f64,
        config: EnthalpyStepConfig) -> Result<EnthalpyStepSolution, HeterogeneousEnthalpyError>
    {
        let boundary = events::boundary_before(cx, self.run.problem.boundary, to).map_err(conduction)?;
        self.run.storage.admit_prescribed_boundary(cx, &boundary)?;
        let source = events::source_before(&self.run.problem.source, to, old.len()).map_err(conduction)?;
        self.solve(cx, old, &boundary, &source, to - from, to, true, config)
    }

    #[allow(clippy::too_many_arguments)]
    fn solve(&self, cx: &Cx<'_>, old: &[f64], boundary: &ThermalBoundary,
        source: &ScalarField, dt_s: f64, time_s: f64, before: bool, config: EnthalpyStepConfig)
        -> Result<EnthalpyStepSolution, HeterogeneousEnthalpyError>
    {
        let patches = self.sample(cx, time_s, before).map_err(conduction)?;
        let mesh = self.run.problem.boundary.mesh();
        let radiation = MaterialPointwiseAmbientRadiation::new(cx, mesh, boundary, &patches)
            .map_err(conduction)?;
        Ok(self.run.storage.advance_with_material_ambient_radiation(cx,
            ConductionProblem { mesh, boundary, source, material: self.run.problem.material,
                element_materials: self.run.problem.element_materials },
            self.run.problem.interfaces, old, dt_s, config, &radiation)?)
    }
}
