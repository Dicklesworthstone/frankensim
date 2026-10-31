//! Fatigue observation on the ORIGINAL quasistatic tetrahedral load solve.
//!
//! The path borrows one immutable reference mesh, constitutive material field,
//! support set and thermal state. Each accepted load is solved by `loads3` and
//! the actual recovered tensor feeds fixed-plane S-N and/or Findley gauges.
//! No displacement/stress solver, mesh interpolation or constitutive law is
//! duplicated. Independent S-N gauges retain rainflow residues across calls.
//!
//! Findley searches explicitly declared fixed planes on a closed stress window,
//! retaining both shear-vector directions and independent normal-stress peaks.
//! It reports utilization of a calibrated fatigue limit, not cycles to failure
//! or Miner damage. Closing a window does not establish periodic stabilization.
//! Restart anchors a new window at the last actual accepted native stress; it
//! does not erase fixed-plane history or refund the global sample budget.
//!
//! Both models require fixed-state, small-strain, elastic, rate-independent
//! applicability. Sampling must resolve stress extrema and rotating shear paths;
//! time labels do not turn the quasistatic solve into dynamics. No stiffness
//! degradation, heat or remeshing history is inferred from empirical fatigue.

use std::sync::Arc;

use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_material::fatigue::critical_plane::findley::{
    FindleyCalibration, FindleyCycle, FindleyCycleResult, FindleyError, FindleyLimits,
    StressClosureTolerance,
};
use fs_material::fatigue::{
    FatigueError, FatigueSummary, MeanStressCorrection, NormalStressDirection, RainflowLimits,
    RainflowMiner, SnCurve,
};
use fs_material::state_point::ResolvedMaterialStatePoint;

use crate::linear3::{
    ElasticTensorOrder, StressTensorBasis, StressTensorNotation, TetLinearElasticProblem,
    TetMaterialField, TetStaticSolveConfig, TetThermalStrainField, TetThermalStressReport,
};
use crate::loads3::transient::TransientError3;
use crate::loads3::{LoadError3, TetLoadCase3, TetLoadSolution3};

/// An explicitly calibrated fixed-plane observation at one original element.
#[derive(Debug, Clone)]
pub struct TetFatigueGauge {
    pub element: usize,
    /// Components in the path's declared reference-mesh frame.
    pub normal_in_mesh_frame: [f64; 3],
    /// Exact elastic state to which the caller applies this calibration.
    pub material_state_identity: [u8; 32],
    pub curve: Arc<SnCurve>,
    pub mean_stress: MeanStressCorrection,
}

/// A shear-sensitive fatigue-limit search on an actual closed stress window.
/// The finite plane set uses the same fixed mesh frame as the native response.
#[derive(Debug, Clone)]
pub struct TetFindleyGauge {
    pub element: usize,
    pub material_state_identity: [u8; 32],
    /// Original immutable material bundle whose identity the elastic field
    /// actually carries. Mandatory for a source-resolved calibration; its card,
    /// named state and complete typed query must match that calibration's source.
    /// Derived/oriented material identities need their own lineage receipt and
    /// are not accepted merely because a material name or temperature agrees.
    pub material_source: Option<Arc<ResolvedMaterialStatePoint>>,
    pub calibration: Arc<FindleyCalibration>,
    pub normals_in_mesh_frame: Vec<[f64; 3]>,
    /// Separate per-gauge limits bound retained tensors, plane projection and
    /// minimum-enclosing-circle work. Transactional trials need another copy.
    pub limits: FindleyLimits,
    pub closure: StressClosureTolerance,
}

/// Explicit mechanics and bounded history controls; no invented life defaults.
#[derive(Debug, Clone, Copy)]
pub struct TetFatigueConfig {
    pub solver: TetStaticSolveConfig,
    pub maximum_absolute_mandel_strain: f64,
    /// The sample limit covers the entire path, including Findley-only paths.
    pub history: RainflowLimits,
    /// Combined fixed-plane S-N and Findley gauge count.
    pub maximum_gauges: usize,
    /// Aggregate fixed-plane reversal slots per retained copy. Findley gauges
    /// declare separate tensor-storage and work limits in `TetFindleyGauge`.
    pub maximum_total_residue: usize,
}

/// Observation limits shared by quasistatic and native inertial fatigue paths.
/// Mechanics keeps its original solver controls; these limits govern only
/// stress recovery and retained empirical histories.
#[derive(Debug, Clone, Copy)]
pub struct TetFatigueHistoryConfig {
    pub maximum_absolute_mandel_strain: f64,
    /// The sample limit is global and is not reset by a new Findley window.
    pub history: RainflowLimits,
    /// Combined fixed-plane S-N and Findley gauge count.
    pub maximum_gauges: usize,
    /// Aggregate fixed-plane reversal slots per retained copy. A transactional
    /// candidate needs another bounded copy; Findley limits are per gauge.
    pub maximum_total_residue: usize,
}

impl TetFatigueConfig {
    fn history_config(self) -> TetFatigueHistoryConfig {
        TetFatigueHistoryConfig {
            maximum_absolute_mandel_strain: self.maximum_absolute_mandel_strain,
            history: self.history,
            maximum_gauges: self.maximum_gauges,
            maximum_total_residue: self.maximum_total_residue,
        }
    }
}

#[derive(Debug)]
pub enum TetFatigueError {
    Invalid(&'static str),
    Budget,
    SourceMismatch {
        gauge: usize,
        what: &'static str,
    },
    Gauge {
        gauge: usize,
        error: FatigueError,
    },
    /// Index within the declared Findley bank, independent of S-N gauge count.
    Findley {
        gauge: usize,
        error: FindleyError,
    },
    Load(LoadError3),
    Transient(TransientError3),
    Cancelled,
    Finished,
}
impl core::fmt::Display for TetFatigueError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "tetrahedral fatigue path refused: {self:?}")
    }
}
impl std::error::Error for TetFatigueError {}
impl From<LoadError3> for TetFatigueError {
    fn from(value: LoadError3) -> Self {
        Self::Load(value)
    }
}
impl From<TransientError3> for TetFatigueError {
    fn from(value: TransientError3) -> Self {
        Self::Transient(value)
    }
}

#[derive(Debug, Clone, PartialEq)]
struct GaugeState {
    element: usize,
    direction: NormalStressDirection,
    counter: RainflowMiner,
}

#[derive(Debug, Clone)]
struct FindleyGaugeState {
    element: usize,
    material_source: Option<Arc<ResolvedMaterialStatePoint>>,
    counter: FindleyCycle,
    /// Only the private native-report acceptance path can write this anchor.
    last_stress_pa: Option<[f64; 6]>,
}

/// Native mechanical response and the same accepted signed S-N samples.
/// Findley window results are available through the owning path's accessor.
#[derive(Debug)]
pub struct TetFatigueStep {
    pub time_s: f64,
    pub mechanics: TetLoadSolution3,
    /// In the original S-N gauge order; positive means tension on that plane.
    pub signed_stress_pa: Vec<f64>,
    /// Open-history damage: unresolved final half cycles are not included yet.
    pub fatigue: Vec<FatigueSummary>,
}

/// Mesh/material/support/temperature binding cannot be replaced during a path.
/// A failed solve, source refusal or cancellation publishes no partial gauges.
pub struct TetFatiguePath<'p, 'm> {
    problem: &'p TetLinearElasticProblem<'m>,
    thermal: TetThermalStrainField<'p>,
    solver: TetStaticSolveConfig,
    observation: TetFatigueObserver,
}

/// Shared private observer. Its callers own the original mechanical binding;
/// accepting arbitrary external stress arrays is deliberately not public API.
#[derive(Debug, Clone)]
pub(crate) struct TetFatigueObserver {
    mesh_frame: ContentHash,
    config: TetFatigueHistoryConfig,
    element_count: usize,
    gauges: Vec<GaugeState>,
    findley: Vec<FindleyGaugeState>,
    samples: u64,
    last_time_s: Option<f64>,
    finished: bool,
}

impl<'p, 'm> TetFatiguePath<'p, 'm> {
    pub fn new(
        problem: &'p TetLinearElasticProblem<'m>,
        thermal: TetThermalStrainField<'p>,
        mesh_frame: ContentHash,
        gauges: Vec<TetFatigueGauge>,
        config: TetFatigueConfig,
        cx: &Cx<'_>,
    ) -> Result<Self, TetFatigueError> {
        Self::new_with_findley(problem, thermal, mesh_frame, gauges, Vec::new(), config, cx)
    }

    /// Observe S-N gauges, Findley gauges or both on the same actual response.
    /// At least one gauge must be declared; no fatigue law is chosen implicitly.
    pub fn new_with_findley(
        problem: &'p TetLinearElasticProblem<'m>,
        thermal: TetThermalStrainField<'p>,
        mesh_frame: ContentHash,
        gauges: Vec<TetFatigueGauge>,
        findley: Vec<TetFindleyGauge>,
        config: TetFatigueConfig,
        cx: &Cx<'_>,
    ) -> Result<Self, TetFatigueError> {
        poll(cx)?;
        if !config.solver.relative_residual_tolerance.is_finite()
            || config.solver.relative_residual_tolerance <= 0.0
            || config.solver.relative_residual_tolerance >= 1.0
            || config.solver.maximum_iterations == 0
        {
            return Err(TetFatigueError::Invalid(
                "declare valid quasistatic solve limits",
            ));
        }
        let observation = TetFatigueObserver::new(
            problem,
            thermal,
            mesh_frame,
            gauges,
            findley,
            config.history_config(),
            cx,
        )?;
        Ok(Self {
            problem,
            thermal,
            solver: config.solver,
            observation,
        })
    }

    pub fn last_time_s(&self) -> Option<f64> {
        self.observation.last_time_s()
    }

    /// Number of actual native samples; restarting a window does not change it.
    pub fn observed_samples(&self) -> u64 {
        self.observation.observed_samples()
    }

    /// Complete original S-N calibration, including its declared scope.
    pub fn curve(&self, gauge: usize) -> Option<&SnCurve> {
        self.observation.curve(gauge)
    }

    pub fn summaries(&self) -> Vec<FatigueSummary> {
        self.observation.summaries()
    }

    pub fn findley_results(&self) -> Vec<Option<FindleyCycleResult>> {
        self.observation.findley_results()
    }

    pub fn findley_calibration(&self, gauge: usize) -> Option<&FindleyCalibration> {
        self.observation.findley_calibration(gauge)
    }

    pub fn findley_material_source(&self, gauge: usize) -> Option<&ResolvedMaterialStatePoint> {
        self.observation.findley_material_source(gauge)
    }

    /// Close just the accepted Findley window. S-N histories remain open.
    /// Every Findley gauge must accept before any result is published.
    pub fn close_findley_cycle(
        &mut self,
        cx: &Cx<'_>,
    ) -> Result<Vec<FindleyCycleResult>, TetFatigueError> {
        self.observation.close_findley_cycle(cx)
    }

    /// Start a fresh Findley window from the last actual accepted native stress.
    /// An open warm-up window may be discarded explicitly. Fixed-plane histories,
    /// accepted time and the global sample budget remain continuous.
    pub fn restart_findley_cycle(&mut self, cx: &Cx<'_>) -> Result<(), TetFatigueError> {
        self.observation.restart_findley_cycle(cx)
    }

    /// Solve one load and accept every gauge together. The FIRST call is the
    /// first sample, not an inferred zero-load state. Times must strictly grow.
    pub fn advance(
        &mut self,
        time_s: f64,
        load: &TetLoadCase3<'_, '_>,
        cx: &Cx<'_>,
    ) -> Result<TetFatigueStep, TetFatigueError> {
        self.advance_inner(time_s, load, false, cx)
    }

    /// Accept the endpoint and close all Findley windows in one transaction.
    /// Closure or circle-work refusal leaves time, S-N histories and every
    /// Findley window at the previous actual accepted sample.
    pub fn advance_and_close_findley_cycle(
        &mut self,
        time_s: f64,
        load: &TetLoadCase3<'_, '_>,
        cx: &Cx<'_>,
    ) -> Result<TetFatigueStep, TetFatigueError> {
        self.advance_inner(time_s, load, true, cx)
    }

    fn advance_inner(
        &mut self,
        time_s: f64,
        load: &TetLoadCase3<'_, '_>,
        close_findley: bool,
        cx: &Cx<'_>,
    ) -> Result<TetFatigueStep, TetFatigueError> {
        self.observation.check_next(time_s, cx)?;
        if close_findley {
            self.observation.require_findley(cx)?;
        }
        // Existing owner checks exact mesh equality, source identities, solver
        // residual and strain domain before returning this accepted response.
        let mechanics = load.solve(
            self.problem,
            self.thermal,
            self.solver,
            self.observation.mesh_frame(),
            self.observation.maximum_strain(),
            cx,
        )?;
        let (signed_stress_pa, fatigue) =
            self.observation
                .accept(time_s, &mechanics.stress, close_findley, cx)?;
        Ok(TetFatigueStep {
            time_s,
            mechanics,
            signed_stress_pa,
            fatigue,
        })
    }

    /// Close the entire observation. Both S-N half cycles and Findley closure
    /// must accept before any history is finalized. Results remain inspectable.
    pub fn finish(&mut self, cx: &Cx<'_>) -> Result<Vec<FatigueSummary>, TetFatigueError> {
        self.observation.finish(cx)
    }
}

impl TetFatigueObserver {
    pub(crate) fn new(
        problem: &TetLinearElasticProblem<'_>,
        thermal: TetThermalStrainField<'_>,
        mesh_frame: ContentHash,
        gauges: Vec<TetFatigueGauge>,
        findley: Vec<TetFindleyGauge>,
        config: TetFatigueHistoryConfig,
        cx: &Cx<'_>,
    ) -> Result<Self, TetFatigueError> {
        poll(cx)?;
        if gauges
            .len()
            .checked_add(findley.len())
            .is_none_or(|n| n == 0 || n > config.maximum_gauges)
            || config.history.maximum_samples == 0
            || gauges
                .len()
                .checked_mul(config.history.maximum_residue)
                .is_none_or(|n| n > config.maximum_total_residue)
        {
            return Err(TetFatigueError::Budget);
        }
        if mesh_frame == ContentHash([0; 32])
            || !config.maximum_absolute_mandel_strain.is_finite()
            || config.maximum_absolute_mandel_strain <= 0.0
        {
            return Err(TetFatigueError::Invalid(
                "declare a nonzero frame and valid strain limit",
            ));
        }
        let count = problem.tetrahedra.len();
        if matches!(thermal, TetThermalStrainField::PerElement(s) if s.len() != count)
            || matches!(problem.materials, TetMaterialField::PerElement(s) if s.len() != count)
        {
            return Err(TetFatigueError::Invalid(
                "material/temperature field must follow the original connectivity",
            ));
        }
        let mut states = Vec::with_capacity(gauges.len());
        for (index, gauge) in gauges.into_iter().enumerate() {
            poll(cx)?;
            bind_material(
                problem,
                thermal,
                index,
                gauge.element,
                gauge.material_state_identity,
                gauge.curve.scope().temperature_k,
            )?;
            let direction =
                NormalStressDirection::new(gauge.normal_in_mesh_frame).map_err(|error| {
                    TetFatigueError::Gauge {
                        gauge: index,
                        error,
                    }
                })?;
            let counter = RainflowMiner::new(gauge.curve, gauge.mean_stress, config.history)
                .map_err(|error| TetFatigueError::Gauge {
                    gauge: index,
                    error,
                })?;
            states.push(GaugeState {
                element: gauge.element,
                direction,
                counter,
            });
        }
        let fixed_count = states.len();
        let mut searches = Vec::with_capacity(findley.len());
        for (index, gauge) in findley.into_iter().enumerate() {
            poll(cx)?;
            let source_index = fixed_count + index;
            bind_material(
                problem,
                thermal,
                source_index,
                gauge.element,
                gauge.material_state_identity,
                gauge.calibration.scope().temperature_k,
            )?;
            if let Some(material) = gauge.material_source.as_deref() {
                if material.identity().0 != gauge.material_state_identity {
                    return Err(TetFatigueError::SourceMismatch {
                        gauge: source_index,
                        what: "Findley material receipt is not the actual elastic state identity",
                    });
                }
            }
            if let Some(source) = gauge.calibration.source() {
                let material = gauge.material_source.as_deref().ok_or(
                    TetFatigueError::SourceMismatch {
                        gauge: source_index,
                        what: "source-backed Findley requires the original immutable elastic material bundle",
                    },
                )?;
                if !source.matches_material_state(
                    material.card_identity(),
                    material.material(),
                    &material.query(),
                ) {
                    return Err(TetFatigueError::SourceMismatch {
                        gauge: source_index,
                        what: "Findley source card, material state or typed query differs from the actual elastic bundle",
                    });
                }
            }
            let counter = FindleyCycle::new(
                gauge.calibration,
                gauge.normals_in_mesh_frame,
                gauge.limits,
                gauge.closure,
            )
            .map_err(|error| TetFatigueError::Findley {
                gauge: index,
                error,
            })?;
            searches.push(FindleyGaugeState {
                element: gauge.element,
                material_source: gauge.material_source,
                counter,
                last_stress_pa: None,
            });
        }
        poll(cx)?;
        Ok(Self {
            mesh_frame,
            config,
            element_count: count,
            gauges: states,
            findley: searches,
            samples: 0,
            last_time_s: None,
            finished: false,
        })
    }

    pub(crate) fn mesh_frame(&self) -> ContentHash {
        self.mesh_frame
    }

    pub(crate) fn maximum_strain(&self) -> f64 {
        self.config.maximum_absolute_mandel_strain
    }

    pub(crate) fn last_time_s(&self) -> Option<f64> {
        self.last_time_s
    }

    pub(crate) fn observed_samples(&self) -> u64 {
        self.samples
    }

    pub(crate) fn curve(&self, gauge: usize) -> Option<&SnCurve> {
        self.gauges.get(gauge).map(|g| g.counter.curve())
    }

    pub(crate) fn summaries(&self) -> Vec<FatigueSummary> {
        self.gauges.iter().map(|g| g.counter.summary()).collect()
    }

    pub(crate) fn findley_results(&self) -> Vec<Option<FindleyCycleResult>> {
        self.findley
            .iter()
            .map(|g| g.counter.result().copied())
            .collect()
    }

    pub(crate) fn findley_calibration(&self, gauge: usize) -> Option<&FindleyCalibration> {
        self.findley.get(gauge).map(|g| g.counter.calibration())
    }

    pub(crate) fn findley_material_source(
        &self,
        gauge: usize,
    ) -> Option<&ResolvedMaterialStatePoint> {
        self.findley
            .get(gauge)
            .and_then(|g| g.material_source.as_deref())
    }

    pub(crate) fn require_findley(&self, cx: &Cx<'_>) -> Result<(), TetFatigueError> {
        poll(cx)?;
        if self.finished {
            return Err(TetFatigueError::Finished);
        }
        if self.findley.is_empty() {
            return Err(TetFatigueError::Invalid(
                "no Findley cycle gauges were declared",
            ));
        }
        Ok(())
    }

    pub(crate) fn close_findley_cycle(
        &mut self,
        cx: &Cx<'_>,
    ) -> Result<Vec<FindleyCycleResult>, TetFatigueError> {
        self.require_findley(cx)?;
        finish_findley_atomically(&mut self.findley, cx)
    }

    pub(crate) fn restart_findley_cycle(&mut self, cx: &Cx<'_>) -> Result<(), TetFatigueError> {
        self.require_findley(cx)?;
        // Build empty windows directly: copying an old tensor history merely
        // to discard it would make a restart unnecessarily proportional to H.
        let mut candidate = Vec::with_capacity(self.findley.len());
        for (index, gauge) in self.findley.iter().enumerate() {
            poll(cx)?;
            let sample = gauge.last_stress_pa.ok_or(TetFatigueError::Invalid(
                "a Findley restart requires an actual accepted native sample",
            ))?;
            let mut counter = gauge
                .counter
                .fresh()
                .map_err(|error| TetFatigueError::Findley {
                    gauge: index,
                    error,
                })?;
            counter
                .push_with_cancel(sample, &mut || cx.checkpoint().is_err())
                .map_err(|error| TetFatigueError::Findley {
                    gauge: index,
                    error,
                })?;
            candidate.push(FindleyGaugeState {
                element: gauge.element,
                material_source: gauge.material_source.clone(),
                counter,
                last_stress_pa: Some(sample),
            });
        }
        poll(cx)?;
        self.findley = candidate;
        Ok(())
    }

    pub(crate) fn check_next(&self, time_s: f64, cx: &Cx<'_>) -> Result<(), TetFatigueError> {
        poll(cx)?;
        if self.finished {
            return Err(TetFatigueError::Finished);
        }
        if !time_s.is_finite() || time_s < 0.0 || self.last_time_s.is_some_and(|t| time_s <= t) {
            return Err(TetFatigueError::Invalid(
                "sample times must be finite, nonnegative and strictly increasing",
            ));
        }
        if self.samples >= self.config.history.maximum_samples {
            return Err(TetFatigueError::Budget);
        }
        if self.findley.iter().any(|g| g.counter.result().is_some()) {
            return Err(TetFatigueError::Invalid(
                "restart a closed Findley cycle before advancing",
            ));
        }
        Ok(())
    }

    /// The owning path must recover this response from its original bound
    /// mechanics. All gauge changes and optional closure commit together only
    /// after the final cancellation poll. External tensors cannot enter here.
    pub(crate) fn accept(
        &mut self,
        time_s: f64,
        stress: &TetThermalStressReport,
        close_findley: bool,
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, Vec<FatigueSummary>), TetFatigueError> {
        self.check_next(time_s, cx)?;
        if close_findley {
            self.require_findley(cx)?;
        }
        if stress.mesh_frame() != self.mesh_frame || stress.elements().len() != self.element_count {
            return Err(TetFatigueError::Invalid(
                "stress report differs from the bound frame or element field",
            ));
        }
        let mut candidate = self.gauges.clone();
        let mut signed_stress_pa = Vec::with_capacity(candidate.len());
        for (index, gauge) in candidate.iter_mut().enumerate() {
            poll(cx)?;
            let tensor = stress_tensor(stress, gauge.element, self.mesh_frame)?;
            let signed =
                gauge
                    .direction
                    .project(&tensor)
                    .map_err(|error| TetFatigueError::Gauge {
                        gauge: index,
                        error,
                    })?;
            gauge
                .counter
                .push_with_cancel(signed, &mut || cx.checkpoint().is_err())
                .map_err(|error| TetFatigueError::Gauge {
                    gauge: index,
                    error,
                })?;
            signed_stress_pa.push(signed);
        }
        let fatigue = candidate.iter().map(|g| g.counter.summary()).collect();
        let mut tensors = Vec::with_capacity(self.findley.len());
        let mut checkpoints = Vec::with_capacity(self.findley.len());
        for gauge in &self.findley {
            poll(cx)?;
            tensors.push(stress_tensor(stress, gauge.element, self.mesh_frame)?);
            checkpoints.push(gauge.counter.sample_count());
        }
        // Append tentatively without cloning H previous tensors every time.
        // Every pre-step Findley result is None by check_next. Reverting each
        // original prefix also clears a result finalized before a later refusal.
        let proposal = (|| {
            for (index, (gauge, &tensor)) in self.findley.iter_mut().zip(&tensors).enumerate() {
                gauge
                    .counter
                    .push_with_cancel(tensor, &mut || cx.checkpoint().is_err())
                    .map_err(|error| TetFatigueError::Findley {
                        gauge: index,
                        error,
                    })?;
            }
            if close_findley {
                finish_findley(&mut self.findley, cx)?;
            }
            poll(cx)
        })();
        if let Err(error) = proposal {
            for (gauge, checkpoint) in self.findley.iter_mut().zip(checkpoints) {
                gauge
                    .counter
                    .discard_after(checkpoint)
                    .expect("an accepted prefix remains valid after a tentative append");
            }
            return Err(error);
        }
        for (gauge, tensor) in self.findley.iter_mut().zip(tensors) {
            gauge.last_stress_pa = Some(tensor);
        }
        self.gauges = candidate;
        self.last_time_s = Some(time_s);
        self.samples += 1;
        Ok((signed_stress_pa, fatigue))
    }

    pub(crate) fn finish(&mut self, cx: &Cx<'_>) -> Result<Vec<FatigueSummary>, TetFatigueError> {
        poll(cx)?;
        if self.finished {
            return Ok(self.summaries());
        }
        let mut candidate = self.gauges.clone();
        for (index, gauge) in candidate.iter_mut().enumerate() {
            gauge
                .counter
                .finish_with_cancel(&mut || cx.checkpoint().is_err())
                .map_err(|error| TetFatigueError::Gauge {
                    gauge: index,
                    error,
                })?;
        }
        let summaries = candidate.iter().map(|g| g.counter.summary()).collect();
        // This is the last fallible operation, including its final poll. Its
        // cache-only transaction never copies the retained tensor trajectories.
        finish_findley_atomically(&mut self.findley, cx)?;
        self.gauges = candidate;
        self.finished = true;
        Ok(summaries)
    }
}

fn bind_material(
    problem: &TetLinearElasticProblem<'_>,
    thermal: TetThermalStrainField<'_>,
    gauge: usize,
    element: usize,
    material_state_identity: [u8; 32],
    calibration_temperature_k: f64,
) -> Result<(), TetFatigueError> {
    if element >= problem.tetrahedra.len() {
        return Err(TetFatigueError::Invalid(
            "fatigue gauge is outside the original element field",
        ));
    }
    let material = match problem.materials {
        TetMaterialField::Uniform(m) => m,
        TetMaterialField::PerElement(m) => &m[element],
    };
    let temperature = match thermal {
        TetThermalStrainField::Uniform(t) => t,
        TetThermalStrainField::PerElement(t) => &t[element],
    };
    if material_state_identity == [0; 32]
        || material_state_identity != material.material_state_identity
        || temperature.elastic_material_state_identity() != material.material_state_identity
    {
        return Err(TetFatigueError::SourceMismatch {
            gauge,
            what: "fatigue/thermal binding differs from the actual elastic material state",
        });
    }
    if calibration_temperature_k != temperature.temperature_k() {
        return Err(TetFatigueError::SourceMismatch {
            gauge,
            what: "fatigue calibration temperature differs from the retained solid state",
        });
    }
    Ok(())
}

fn stress_tensor(
    stress: &TetThermalStressReport,
    element: usize,
    frame: ContentHash,
) -> Result<[f64; 6], TetFatigueError> {
    // Original stress-basis owner: Mandel shear is never treated as physical
    // shear, and unsigned equivalent stress never replaces a tensor trajectory.
    let transformed = stress
        .stress_in_frame(
            element,
            StressTensorBasis {
                notation: StressTensorNotation::Tensor,
                order: ElasticTensorOrder::XxYyZzXyYzZx,
                frame,
            },
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        )
        .map_err(|error| TetFatigueError::Load(LoadError3::Solver(error)))?;
    Ok(*transformed.stress_pa())
}

fn finish_findley(
    gauges: &mut [FindleyGaugeState],
    cx: &Cx<'_>,
) -> Result<Vec<FindleyCycleResult>, TetFatigueError> {
    let mut results = Vec::with_capacity(gauges.len());
    for (index, gauge) in gauges.iter_mut().enumerate() {
        results.push(
            gauge
                .counter
                .finish_with_cancel(&mut || cx.checkpoint().is_err())
                .map_err(|error| TetFatigueError::Findley {
                    gauge: index,
                    error,
                })?,
        );
    }
    Ok(results)
}

/// Finalization changes only cached results. Preserve any previously closed
/// windows, and roll newly finalized ones back to their original open prefix
/// if another gauge or the final cancellation poll refuses.
fn finish_findley_atomically(
    gauges: &mut [FindleyGaugeState],
    cx: &Cx<'_>,
) -> Result<Vec<FindleyCycleResult>, TetFatigueError> {
    let checkpoints: Vec<_> = gauges
        .iter()
        .map(|g| (g.counter.sample_count(), g.counter.result().is_none()))
        .collect();
    let proposal = finish_findley(gauges, cx).and_then(|results| {
        poll(cx)?;
        Ok(results)
    });
    if proposal.is_err() {
        for (gauge, (checkpoint, was_open)) in gauges.iter_mut().zip(checkpoints) {
            if was_open {
                gauge
                    .counter
                    .discard_after(checkpoint)
                    .expect("finalization cannot remove an accepted stress sample");
            }
        }
    }
    proposal
}

fn poll(cx: &Cx<'_>) -> Result<(), TetFatigueError> {
    cx.checkpoint().map_err(|_| TetFatigueError::Cancelled)
}
