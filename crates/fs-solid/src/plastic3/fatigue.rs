//! Bounded strain-life observation on the native source-bound J2 load path.
//!
//! Each accepted load factor uses the original `solve_from` owner and retains
//! its plastic strain and isotropic/kinematic history. Gauges observe signed
//! TOTAL MECHANICAL strain in a fixed reference plane, not accumulated plastic
//! strain or unsigned equivalent stress. `ResolvedJ2Plasticity` has no thermal
//! eigenstrain: its source temperature and coefficients stay fixed throughout.
//! Thermal-wrapped laws and changing-temperature paths need their own admitted
//! mechanical-strain transfer and cannot enter this fixed-state adapter.
//!
//! A Coffin-Manson-Basquin calibration is either an explicit caller declaration
//! or an immutable source-card fit. Matching the original material state,
//! process, modulus and temperature prevents accidental rebinding; source fits
//! additionally pin the original card, external query and modulus claim. These
//! checks do not validate the fit or predict component life. Morrow/SWT use
//! actual paired stress/strain and refuse unsupported misaligned histories.
//! Miner damage changes neither stiffness, heat nor failure topology. Native
//! small-strain/P1 and linear-hardening limits remain in force. Resolve every
//! physical turning point in `load_path`; no extrema between targets are inferred.

use std::sync::Arc;

use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_material::fatigue::strain_life::{
    PairedStrainLifeMiner, StrainLifeCurve, StrainLifeMeanStress, StrainLifeMiner,
    StrainLifeSummary, StrainStressSample,
};
use fs_material::fatigue::{FatigueError, NormalStressDirection, RainflowLimits};
use fs_material::plastic::J2State;
use fs_material::plastic::source::ResolvedJ2Plasticity;
use fs_material::state_point::YOUNG_MODULUS_PROPERTY;

use super::{
    SmallStrainLawField, SmallStrainTetError, SmallStrainTetProblem, SmallStrainTetSolution,
};

/// Fixed-plane observation at one original tetrahedron, with explicit source
/// identity and total mechanical strain-life calibration.
#[derive(Debug, Clone)]
pub struct J2FatigueGauge {
    pub element: usize,
    pub normal_in_mesh_frame: [f64; 3],
    pub material_state_identity: ContentHash,
    pub curve: Arc<StrainLifeCurve>,
    pub mean_stress: StrainLifeMeanStress,
}

/// Logical storage/work limits; transactional proposals need an additional
/// copy of the latest native solution and these bounded gauge histories.
#[derive(Debug, Clone, Copy)]
pub struct J2FatigueConfig {
    pub history: RainflowLimits,
    pub maximum_gauges: usize,
    pub maximum_total_residue: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub enum J2FatigueError {
    Invalid(&'static str),
    Budget,
    SourceMismatch { gauge: usize, what: &'static str },
    Gauge { gauge: usize, error: FatigueError },
    Mechanics(SmallStrainTetError),
    Cancelled,
    Finished,
    PathIncomplete,
}
impl core::fmt::Display for J2FatigueError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "cyclic J2 strain-life path refused: {self:?}")
    }
}
impl std::error::Error for J2FatigueError {}
impl From<SmallStrainTetError> for J2FatigueError {
    fn from(error: SmallStrainTetError) -> Self {
        match error {
            SmallStrainTetError::Cancelled => Self::Cancelled,
            error => Self::Mechanics(error),
        }
    }
}

/// Actual accepted observations in one fixed plane. Named Morrow/SWT policies
/// consume this same signed stress together with total mechanical strain.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct J2FatigueSample {
    pub mechanical_strain: f64,
    pub normal_stress_pa: f64,
}

#[derive(Debug, Clone)]
struct GaugeState {
    element: usize,
    normal: NormalStressDirection,
    counter: GaugeCounter,
    sample: J2FatigueSample,
}

#[derive(Debug, Clone)]
enum GaugeCounter {
    Scalar(StrainLifeMiner),
    Paired(PairedStrainLifeMiner),
}
impl GaugeCounter {
    fn new(
        curve: Arc<StrainLifeCurve>,
        policy: StrainLifeMeanStress,
        limits: RainflowLimits,
    ) -> Result<Self, FatigueError> {
        if policy == StrainLifeMeanStress::Uncorrected {
            Ok(Self::Scalar(StrainLifeMiner::new(curve, policy, limits)?))
        } else {
            Ok(Self::Paired(PairedStrainLifeMiner::new(
                curve, policy, limits,
            )?))
        }
    }
    fn curve(&self) -> &StrainLifeCurve {
        match self {
            Self::Scalar(c) => c.curve(),
            Self::Paired(c) => c.curve(),
        }
    }
    fn summary(&self) -> StrainLifeSummary {
        match self {
            Self::Scalar(c) => c.summary(),
            Self::Paired(c) => c.summary(),
        }
    }
    fn push(
        &mut self,
        sample: J2FatigueSample,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), FatigueError> {
        match self {
            Self::Scalar(c) => c.push_with_cancel(sample.mechanical_strain, cancelled),
            Self::Paired(c) => c.push_with_cancel(
                StrainStressSample {
                    mechanical_strain: sample.mechanical_strain,
                    normal_stress_pa: sample.normal_stress_pa,
                },
                cancelled,
            ),
        }
    }
    fn finish(
        &mut self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<StrainLifeSummary, FatigueError> {
        match self {
            Self::Scalar(c) => c.finish_with_cancel(cancelled),
            Self::Paired(c) => c.finish_with_cancel(cancelled),
        }
    }
}

/// Common private source binding and accepted-solution observer for the
/// quasistatic and physical-time J2 owners. No external stress arrays enter.
#[derive(Debug, Clone)]
pub(super) struct J2FatigueObserver {
    mesh_frame: ContentHash,
    gauges: Vec<GaugeState>,
    maximum_samples: u64,
    finished: bool,
}

/// Resumable observation of one immutable native load path. `advance` commits
/// a requested block only after every mechanical increment and gauge accepts.
/// Cloning retains the complete current plastic and unresolved cycle histories.
/// Only the latest native increment record is retained, not the whole trace.
#[derive(Debug, Clone)]
pub struct J2FatiguePath<'p, 'm> {
    problem: &'p SmallStrainTetProblem<'m, ResolvedJ2Plasticity>,
    solution: SmallStrainTetSolution<J2State>,
    observer: J2FatigueObserver,
    accepted_increments: usize,
}

impl<'p, 'm> J2FatiguePath<'p, 'm> {
    /// Solve the real virgin, zero-load state and accept it as the first
    /// observation. The declared load path then supplies absolute targets.
    /// Source identities, labels, fixed positive T and modulus must agree with
    /// the actual resolved constitutive field. No thermal strain is subtracted
    /// or guessed, and no pre-existing plastic/fatigue history is fabricated.
    pub fn new(
        problem: &'p SmallStrainTetProblem<'m, ResolvedJ2Plasticity>,
        mesh_frame: ContentHash,
        gauges: Vec<J2FatigueGauge>,
        config: J2FatigueConfig,
        cx: &Cx<'_>,
    ) -> Result<Self, J2FatigueError> {
        poll(cx)?;
        if mesh_frame == ContentHash([0; 32]) || problem.load_path.is_empty() {
            return Err(J2FatigueError::Invalid(
                "declare a nonzero mesh frame and finite nonempty absolute load path",
            ));
        }
        let samples = u64::try_from(problem.load_path.len())
            .ok()
            .and_then(|n| n.checked_add(1));
        if samples.is_none_or(|n| n > config.history.maximum_samples) {
            return Err(J2FatigueError::Budget);
        }
        for chunk in problem.load_path.chunks(256) {
            poll(cx)?;
            if chunk.iter().any(|x| !x.is_finite()) {
                return Err(J2FatigueError::Invalid("load factors must be finite"));
            }
        }
        let mut observer = J2FatigueObserver::new(problem, mesh_frame, gauges, config, cx)?;
        let solution = with_path(problem, &[0.0]).solve(cx)?;
        observer.observe(problem, &solution, cx)?;
        poll(cx)?;
        Ok(Self {
            problem,
            solution,
            observer,
            accepted_increments: 0,
        })
    }

    pub fn mesh_frame(&self) -> ContentHash {
        self.observer.mesh_frame()
    }
    pub fn accepted_increments(&self) -> usize {
        self.accepted_increments
    }
    pub fn remaining_increments(&self) -> usize {
        self.problem.load_path.len() - self.accepted_increments
    }
    pub fn solution(&self) -> &SmallStrainTetSolution<J2State> {
        &self.solution
    }
    pub fn summaries(&self) -> Vec<StrainLifeSummary> {
        self.observer.summaries()
    }
    pub fn samples(&self) -> Vec<J2FatigueSample> {
        self.observer.samples()
    }
    pub fn curve(&self, gauge: usize) -> Option<&StrainLifeCurve> {
        self.observer.curve(gauge)
    }
    /// Original complete property receipts and native constitutive metadata.
    pub fn material(&self, gauge: usize) -> Option<&ResolvedJ2Plasticity> {
        self.observer.material(self.problem, gauge)
    }

    /// Accept exactly `count` planned increments as one transaction. Every
    /// original load factor is observed; callback boundaries close no cycles.
    /// Nonconvergence, cancellation, strain/calibration refusal or any gauge
    /// failure preserves the complete prior mechanics and fatigue histories.
    pub fn advance(&mut self, count: usize, cx: &Cx<'_>) -> Result<(), J2FatigueError> {
        poll(cx)?;
        if self.observer.finished {
            return Err(J2FatigueError::Finished);
        }
        if count == 0 || count > self.remaining_increments() {
            return Err(J2FatigueError::Invalid(
                "advance needs a nonzero count within the remaining load path",
            ));
        }
        let end = self.accepted_increments + count;
        let mut solution = self.solution.clone();
        let mut observer = self.observer.clone();
        for &factor in &self.problem.load_path[self.accepted_increments..end] {
            poll(cx)?;
            solution = with_path(self.problem, &[factor]).solve_from(
                solution.displacement_m.as_flattened(),
                &solution.element_states,
                cx,
            )?;
            observer.observe(self.problem, &solution, cx)?;
        }
        poll(cx)?;
        self.solution = solution;
        self.observer = observer;
        self.accepted_increments = end;
        Ok(())
    }

    /// Finish the entire declared observation once, after all planned targets.
    /// Open endpoint half cycles are finalized together for every gauge.
    pub fn finish(&mut self, cx: &Cx<'_>) -> Result<Vec<StrainLifeSummary>, J2FatigueError> {
        poll(cx)?;
        if self.remaining_increments() != 0 {
            return Err(J2FatigueError::PathIncomplete);
        }
        self.observer.finish(cx)
    }
}

impl J2FatigueObserver {
    pub(super) fn new(
        problem: &SmallStrainTetProblem<'_, ResolvedJ2Plasticity>,
        mesh_frame: ContentHash,
        gauges: Vec<J2FatigueGauge>,
        config: J2FatigueConfig,
        cx: &Cx<'_>,
    ) -> Result<Self, J2FatigueError> {
        poll(cx)?;
        if mesh_frame == ContentHash([0; 32]) {
            return Err(J2FatigueError::Invalid("declare a nonzero mesh frame"));
        }
        if gauges.is_empty()
            || gauges.len() > config.maximum_gauges
            || gauges
                .len()
                .checked_mul(config.history.maximum_residue)
                .is_none_or(|n| n > config.maximum_total_residue)
        {
            return Err(J2FatigueError::Budget);
        }
        problem.laws.validate_count(problem.tetrahedra.len())?;
        let mut prepared = Vec::with_capacity(gauges.len());
        for (index, gauge) in gauges.into_iter().enumerate() {
            poll(cx)?;
            if gauge.element >= problem.tetrahedra.len() {
                return Err(J2FatigueError::Invalid(
                    "gauge is outside the original tetrahedral field",
                ));
            }
            let law = problem.laws.at(gauge.element);
            let material = law.material();
            let scope = gauge.curve.scope();
            if gauge.material_state_identity != material.identity()
                || scope.material != material.material().chemistry
                || scope.process != material.material().process
            {
                return Err(J2FatigueError::SourceMismatch {
                    gauge: index,
                    what: "strain-life identity/material/process differs from the original J2 state",
                });
            }
            let temperature_axis = gauge
                .curve
                .source()
                .map_or("T", |source| source.request().temperature_axis.as_str());
            if material.state_coordinate(temperature_axis) != Some(scope.temperature_k)
                || gauge.curve.parameters().youngs_modulus_pa != law.law().elastic.youngs
            {
                return Err(J2FatigueError::SourceMismatch {
                    gauge: index,
                    what: "strain-life temperature/modulus differs from the frozen J2 law",
                });
            }
            if let Some(source) = gauge.curve.source() {
                let youngs = material.property(YOUNG_MODULUS_PROPERTY).ok_or(
                    J2FatigueError::SourceMismatch {
                        gauge: index,
                        what: "original J2 modulus receipt is missing or ambiguous",
                    },
                )?;
                if !source.matches_material_state(
                    material.card_identity(),
                    material.material(),
                    &material.query(),
                    youngs.answer(),
                ) {
                    return Err(J2FatigueError::SourceMismatch {
                        gauge: index,
                        what: "CMB card, state, external conditions or modulus claim differs from the original J2 source",
                    });
                }
            }
            let normal =
                NormalStressDirection::new(gauge.normal_in_mesh_frame).map_err(|error| {
                    J2FatigueError::Gauge {
                        gauge: index,
                        error,
                    }
                })?;
            let counter = GaugeCounter::new(gauge.curve, gauge.mean_stress, config.history)
                .map_err(|error| J2FatigueError::Gauge {
                    gauge: index,
                    error,
                })?;
            prepared.push(GaugeState {
                element: gauge.element,
                normal,
                counter,
                sample: J2FatigueSample {
                    mechanical_strain: 0.0,
                    normal_stress_pa: 0.0,
                },
            });
        }
        poll(cx)?;
        Ok(Self {
            mesh_frame,
            gauges: prepared,
            maximum_samples: config.history.maximum_samples,
            finished: false,
        })
    }
    pub(super) fn mesh_frame(&self) -> ContentHash {
        self.mesh_frame
    }
    pub(super) fn summaries(&self) -> Vec<StrainLifeSummary> {
        self.gauges.iter().map(|g| g.counter.summary()).collect()
    }
    pub(super) fn samples(&self) -> Vec<J2FatigueSample> {
        self.gauges.iter().map(|g| g.sample).collect()
    }
    pub(super) fn curve(&self, gauge: usize) -> Option<&StrainLifeCurve> {
        self.gauges.get(gauge).map(|g| g.counter.curve())
    }
    pub(super) fn material<'a>(
        &self,
        problem: &'a SmallStrainTetProblem<'_, ResolvedJ2Plasticity>,
        gauge: usize,
    ) -> Option<&'a ResolvedJ2Plasticity> {
        self.gauges.get(gauge).map(|g| problem.laws.at(g.element))
    }
    pub(super) fn check_sample(&self, cx: &Cx<'_>) -> Result<(), J2FatigueError> {
        poll(cx)?;
        if self.finished {
            return Err(J2FatigueError::Finished);
        }
        if self
            .gauges
            .iter()
            .any(|g| g.counter.summary().samples >= self.maximum_samples)
        {
            return Err(J2FatigueError::Budget);
        }
        Ok(())
    }
    /// Called only on an uncommitted candidate owned by the mechanical path.
    pub(super) fn observe(
        &mut self,
        problem: &SmallStrainTetProblem<'_, ResolvedJ2Plasticity>,
        solution: &SmallStrainTetSolution<J2State>,
        cx: &Cx<'_>,
    ) -> Result<(), J2FatigueError> {
        self.check_sample(cx)?;
        observe(problem, solution, &mut self.gauges, cx)
    }
    pub(super) fn finish(&mut self, cx: &Cx<'_>) -> Result<Vec<StrainLifeSummary>, J2FatigueError> {
        poll(cx)?;
        if self.finished {
            return Ok(self.summaries());
        }
        let mut gauges = self.gauges.clone();
        for (index, gauge) in gauges.iter_mut().enumerate() {
            gauge
                .counter
                .finish(&mut || cx.checkpoint().is_err())
                .map_err(|error| J2FatigueError::Gauge {
                    gauge: index,
                    error,
                })?;
        }
        poll(cx)?;
        self.gauges = gauges;
        self.finished = true;
        Ok(self.summaries())
    }
}

pub(super) fn poll(cx: &Cx<'_>) -> Result<(), J2FatigueError> {
    cx.checkpoint().map_err(|_| J2FatigueError::Cancelled)
}

pub(super) fn with_path<'a>(
    problem: &'a SmallStrainTetProblem<'_, ResolvedJ2Plasticity>,
    path: &'a [f64],
) -> SmallStrainTetProblem<'a, ResolvedJ2Plasticity> {
    SmallStrainTetProblem {
        nodes_m: problem.nodes_m,
        tetrahedra: problem.tetrahedra,
        laws: match &problem.laws {
            SmallStrainLawField::Uniform(law) => SmallStrainLawField::Uniform(law),
            SmallStrainLawField::PerElement(laws) => SmallStrainLawField::PerElement(laws),
        },
        prescribed_m: problem.prescribed_m,
        nodal_forces_n: problem.nodal_forces_n,
        load_path: path,
        budget: problem.budget,
        settings: problem.settings,
    }
}

fn observe(
    problem: &SmallStrainTetProblem<'_, ResolvedJ2Plasticity>,
    solution: &SmallStrainTetSolution<J2State>,
    gauges: &mut [GaugeState],
    cx: &Cx<'_>,
) -> Result<(), J2FatigueError> {
    let count = problem.tetrahedra.len();
    if solution.element_strain.len() != count
        || solution.element_stress.len() != count
        || solution.element_states.len() != count
        || solution.displacement_m.len() != problem.nodes_m.len()
    {
        return Err(J2FatigueError::Invalid(
            "native response does not cover the original mesh",
        ));
    }
    // Enforce the original law's actual componentwise tensor-Voigt envelope at
    // EVERY element, including unobserved ones. The generic native solver
    // deliberately leaves this model-domain obligation to its owning caller.
    for (index, strain) in solution.element_strain.iter().enumerate() {
        if index % 256 == 0 {
            poll(cx)?;
        }
        let limit = problem.laws.at(index).law().elastic.strain_limit.min(0.05);
        let state = &solution.element_states[index];
        if strain.iter().any(|x| !x.is_finite() || x.abs() > limit)
            || solution.element_stress[index]
                .iter()
                .any(|x| !x.is_finite())
            || !state.alpha.is_finite()
            || state.alpha < 0.0
            || state.plastic_strain.iter().any(|x| !x.is_finite())
        {
            return Err(J2FatigueError::Invalid(
                "native response exceeds the admitted J2 small-strain/finite-state domain",
            ));
        }
    }
    for (index, gauge) in gauges.iter_mut().enumerate() {
        let strain = gauge
            .normal
            .project_strain(&solution.element_strain[gauge.element])
            .map_err(|error| J2FatigueError::Gauge {
                gauge: index,
                error,
            })?;
        let stress = gauge
            .normal
            .project(&solution.element_stress[gauge.element])
            .map_err(|error| J2FatigueError::Gauge {
                gauge: index,
                error,
            })?;
        let sample = J2FatigueSample {
            mechanical_strain: strain,
            normal_stress_pa: stress,
        };
        gauge
            .counter
            .push(sample, &mut || cx.checkpoint().is_err())
            .map_err(|error| J2FatigueError::Gauge {
                gauge: index,
                error,
            })?;
        gauge.sample = sample;
    }
    poll(cx)?;
    Ok(())
}
