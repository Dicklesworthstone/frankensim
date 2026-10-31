//! Stress-life observation on the ORIGINAL quasistatic tetrahedral load solve.
//!
//! The path borrows one immutable reference mesh, constitutive material field,
//! support set and thermal state. Each accepted load is solved by `loads3`, and
//! the actual recovered stress is projected onto explicitly chosen fixed planes.
//! No displacement/stress solver, mesh interpolation, or constitutive law is
//! duplicated. Independent gauges retain separate rainflow residues across calls.
//!
//! This is a fixed-state, small-strain, elastic, rate-independent stress-life
//! baseline. Sampling must resolve load reversals; time labels do not turn the
//! quasistatic solve into dynamics. Gauge selection is NOT a critical-plane or
//! multiaxial failure criterion. Calibration applicability remains the caller's
//! responsibility; matching a material identity and temperature is necessary,
//! not sufficient evidence of component life. No stiffness degradation or heat
//! is inferred from Miner damage, and no history is transferred through remeshing.

use std::sync::Arc;

use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_material::fatigue::{
    FatigueError, FatigueSummary, MeanStressCorrection, NormalStressDirection, RainflowLimits,
    RainflowMiner, SnCurve,
};

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

/// Explicit mechanics and bounded history controls; no invented life defaults.
#[derive(Debug, Clone, Copy)]
pub struct TetFatigueConfig {
    pub solver: TetStaticSolveConfig,
    pub maximum_absolute_mandel_strain: f64,
    pub history: RainflowLimits,
    pub maximum_gauges: usize,
    /// Aggregate logical reversal slots per retained copy; candidates need additional copies.
    pub maximum_total_residue: usize,
}

/// Observation limits shared by quasistatic and native inertial fatigue paths.
/// Mechanics keeps its original solver controls; these limits govern only
/// stress recovery and the retained empirical histories.
#[derive(Debug, Clone, Copy)]
pub struct TetFatigueHistoryConfig {
    pub maximum_absolute_mandel_strain: f64,
    pub history: RainflowLimits,
    pub maximum_gauges: usize,
    /// Aggregate logical reversal slots per retained copy; a transactional
    /// candidate needs an additional bounded copy.
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
    SourceMismatch { gauge: usize, what: &'static str },
    Gauge { gauge: usize, error: FatigueError },
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

/// Native mechanical response and the same accepted signed fatigue samples.
#[derive(Debug)]
pub struct TetFatigueStep {
    pub time_s: f64,
    pub mechanics: TetLoadSolution3,
    /// In the original gauge order; positive means tension on that plane.
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

    /// Complete original calibration for a gauge, including its declared scope.
    pub fn curve(&self, gauge: usize) -> Option<&SnCurve> {
        self.observation.curve(gauge)
    }

    pub fn summaries(&self) -> Vec<FatigueSummary> {
        self.observation.summaries()
    }

    /// Solve one load and accept every gauge together. The FIRST call is the
    /// first sample, not an inferred zero-load state. Times must strictly grow.
    pub fn advance(
        &mut self,
        time_s: f64,
        load: &TetLoadCase3<'_, '_>,
        cx: &Cx<'_>,
    ) -> Result<TetFatigueStep, TetFatigueError> {
        self.observation.check_next(time_s, cx)?;
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
        let (signed_stress_pa, fatigue) = self.observation.accept(time_s, &mechanics.stress, cx)?;
        Ok(TetFatigueStep {
            time_s,
            mechanics,
            signed_stress_pa,
            fatigue,
        })
    }

    /// Close the entire observation, NOT an individual load block. Every gauge
    /// must admit its endpoint half cycles before any finalization is published.
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
        config: TetFatigueHistoryConfig,
        cx: &Cx<'_>,
    ) -> Result<Self, TetFatigueError> {
        poll(cx)?;
        if gauges.is_empty()
            || gauges.len() > config.maximum_gauges
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
            if gauge.element >= count {
                return Err(TetFatigueError::Invalid(
                    "fatigue gauge is outside the original element field",
                ));
            }
            let material = match problem.materials {
                TetMaterialField::Uniform(m) => m,
                TetMaterialField::PerElement(m) => &m[gauge.element],
            };
            let temperature = match thermal {
                TetThermalStrainField::Uniform(t) => t,
                TetThermalStrainField::PerElement(t) => &t[gauge.element],
            };
            if gauge.material_state_identity == [0; 32]
                || gauge.material_state_identity != material.material_state_identity
                || temperature.elastic_material_state_identity() != material.material_state_identity
            {
                return Err(TetFatigueError::SourceMismatch {
                    gauge: index,
                    what: "fatigue/thermal binding differs from the actual elastic material state",
                });
            }
            if gauge.curve.scope().temperature_k != temperature.temperature_k() {
                return Err(TetFatigueError::SourceMismatch {
                    gauge: index,
                    what: "fatigue calibration temperature differs from the retained solid state",
                });
            }
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
        poll(cx)?;
        Ok(Self {
            mesh_frame,
            config,
            element_count: count,
            gauges: states,
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

    /// Complete original calibration for a gauge, including its declared scope.
    pub(crate) fn curve(&self, gauge: usize) -> Option<&SnCurve> {
        self.gauges.get(gauge).map(|g| g.counter.curve())
    }

    pub(crate) fn summaries(&self) -> Vec<FatigueSummary> {
        self.gauges.iter().map(|g| g.counter.summary()).collect()
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
        if self
            .gauges
            .iter()
            .any(|g| g.counter.summary().samples >= self.config.history.maximum_samples)
        {
            return Err(TetFatigueError::Budget);
        }
        Ok(())
    }

    /// The owning path must recover this response from its original bound
    /// mechanics. Time, field size and frame checks additionally guard internal
    /// misuse. All gauge changes are committed together after the final poll.
    pub(crate) fn accept(
        &mut self,
        time_s: f64,
        stress: &TetThermalStressReport,
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, Vec<FatigueSummary>), TetFatigueError> {
        self.check_next(time_s, cx)?;
        if stress.mesh_frame() != self.mesh_frame || stress.elements().len() != self.element_count {
            return Err(TetFatigueError::Invalid(
                "stress report differs from the bound frame or element field",
            ));
        }
        let mut candidate = self.gauges.clone();
        let mut signed_stress_pa = Vec::with_capacity(candidate.len());
        for (index, gauge) in candidate.iter_mut().enumerate() {
            poll(cx)?;
            // Use the original stress-basis owner; do not mistake Mandel shear
            // for physical shear or replace reversal by unsigned von Mises.
            let stress = stress
                .stress_in_frame(
                    gauge.element,
                    StressTensorBasis {
                        notation: StressTensorNotation::Tensor,
                        order: ElasticTensorOrder::XxYyZzXyYzZx,
                        frame: self.mesh_frame,
                    },
                    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                )
                .map_err(|error| TetFatigueError::Load(LoadError3::Solver(error)))?;
            let signed = gauge
                .direction
                .project(stress.stress_pa())
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
        poll(cx)?;
        self.gauges = candidate;
        self.last_time_s = Some(time_s);
        Ok((signed_stress_pa, fatigue))
    }

    /// Close the entire observation, NOT an individual load block. Every gauge
    /// must admit its endpoint half cycles before any finalization is published.
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
        poll(cx)?;
        self.gauges = candidate;
        self.finished = true;
        Ok(summaries)
    }
}

fn poll(cx: &Cx<'_>) -> Result<(), TetFatigueError> {
    cx.checkpoint().map_err(|_| TetFatigueError::Cancelled)
}
