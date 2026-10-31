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
//! checks do not validate the fit or predict component life. Fixed-plane
//! Morrow/SWT rainflow uses actual paired stress/strain and refuses unsupported
//! misaligned histories. Optional closed-cycle SWT and Fatemi-Socie searches
//! use the full native tensors and independently sampled normal-stress peaks.
//! Fatemi-Socie uses a separately calibrated torsional strain-life relation,
//! retaining engineering-shear amplitude and its named normal-stress policy. Local tensor closure
//! alone does not establish plastic stabilization or repeatable component life.
//! Fatigue damage changes neither stiffness, heat nor failure topology. Native
//! small-strain/P1 and linear-hardening limits remain in force. Resolve every
//! physical turning point in `load_path`; no extrema between targets are inferred.

use std::sync::Arc;

use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_material::fatigue::critical_plane::fatemi_socie::{
    FatemiSocieCalibration, FatemiSocieCycle, FatemiSocieCycleResult, FatemiSocieError,
    FatemiSocieLimits,
};
use fs_material::fatigue::critical_plane::miner::{
    CriticalPlaneMiner, CriticalPlaneMinerError, CriticalPlaneMinerLimits, CriticalPlaneMinerResult,
};
use fs_material::fatigue::critical_plane::{
    CriticalPlaneError, CriticalPlaneLimits, CycleClosureTolerance, SwtCycle, SwtCycleResult,
    TensorFatigueSample,
};
use fs_material::fatigue::strain_life::{
    PairedStrainLifeMiner, StrainLifeCurve, StrainLifeMeanStress, StrainLifeMiner,
    StrainLifeSummary, StrainStressSample,
};
use fs_material::fatigue::{FatigueError, NormalStressDirection, RainflowLimits};
use fs_material::plastic::J2State;
use fs_material::plastic::source::ResolvedJ2Plasticity;
use fs_material::state_point::{
    POISSON_RATIO_PROPERTY, YIELD_STRESS_PROPERTY, YOUNG_MODULUS_PROPERTY,
};

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

/// A declared finite plane search on one CLOSED representative tensor cycle.
/// Plane directions use the same fixed mesh frame as the native solution.
/// This retains the source-bound CMB fit; no repeated plastic state is inferred.
#[derive(Debug, Clone)]
pub struct J2CriticalPlaneGauge {
    pub element: usize,
    pub material_state_identity: ContentHash,
    pub curve: Arc<StrainLifeCurve>,
    pub normals_in_mesh_frame: Vec<[f64; 3]>,
    pub limits: CriticalPlaneLimits,
    pub closure: CycleClosureTolerance,
}

/// A shear-strain life search on one closed representative native tensor cycle.
/// Plane directions share the original mesh frame. The torsional calibration
/// is distinct from an axial CMB fit and includes its initial-yield reference.
#[derive(Debug, Clone)]
pub struct J2FatemiSocieGauge {
    pub element: usize,
    pub material_state_identity: ContentHash,
    pub calibration: Arc<FatemiSocieCalibration>,
    pub normals_in_mesh_frame: Vec<[f64; 3]>,
    pub limits: FatemiSocieLimits,
    pub closure: CycleClosureTolerance,
}

/// Variable-amplitude finite-plane Miner accumulation on one native element.
/// Every accepted J2 endpoint enters every declared plane; no representative
/// cycle closure or stabilization is fabricated.
#[derive(Debug, Clone)]
pub struct J2VariableAmplitudeGauge {
    pub element: usize,
    pub material_state_identity: ContentHash,
    pub curve: Arc<StrainLifeCurve>,
    pub mean_stress: StrainLifeMeanStress,
    pub normals_in_mesh_frame: Vec<[f64; 3]>,
    pub limits: CriticalPlaneMinerLimits,
}

/// Logical storage/work limits. Native proposals copy the mechanical state
/// and bounded fixed/SWT state; shear histories use accepted-prefix rollback.
#[derive(Debug, Clone, Copy)]
pub struct J2FatigueConfig {
    pub history: RainflowLimits,
    /// Combined fixed-plane, SWT, Fatemi-Socie and variable-amplitude gauge count.
    pub maximum_gauges: usize,
    /// Total bounded rainflow residue slots across fixed-plane gauges and every
    /// plane of every variable-amplitude gauge. Representative-cycle banks have
    /// their own explicit plane/sample/work limits.
    pub maximum_total_residue: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub enum J2FatigueError {
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
    CriticalPlane {
        gauge: usize,
        error: CriticalPlaneError,
    },
    FatemiSocie {
        gauge: usize,
        error: FatemiSocieError,
    },
    VariableAmplitude {
        gauge: usize,
        error: CriticalPlaneMinerError,
    },
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
struct CriticalGaugeState {
    element: usize,
    counter: SwtCycle,
}

#[derive(Debug, Clone)]
struct FatemiSocieGaugeState {
    element: usize,
    counter: FatemiSocieCycle,
}

#[derive(Debug, Clone)]
struct VariableAmplitudeGaugeState {
    element: usize,
    counter: CriticalPlaneMiner,
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
    critical: Vec<CriticalGaugeState>,
    fatemi_socie: Vec<FatemiSocieGaugeState>,
    variable_amplitude: Vec<VariableAmplitudeGaugeState>,
    samples: u64,
    maximum_samples: u64,
    finished: bool,
}

/// Only bounded fixed/SWT state is copied. Each shear-history checkpoint is
/// an accepted prefix count; rollback never duplicates the retained tensors.
pub(super) struct J2FatigueCheckpoint {
    gauges: Vec<GaugeState>,
    critical: Vec<CriticalGaugeState>,
    shear_prefixes: Vec<(u64, bool)>,
    variable_amplitude: Vec<VariableAmplitudeGaugeState>,
    samples: u64,
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
        Self::new_with_critical_planes(problem, mesh_frame, gauges, Vec::new(), config, cx)
    }

    /// Observe fixed-plane rainflow, closed-cycle SWT searches, or both on the
    /// same native trajectory. Every bank starts at the real zero-load state.
    /// After conditioning, explicitly restart its cycle at an accepted native
    /// checkpoint; no supplied tensor or imported constitutive state is accepted.
    pub fn new_with_critical_planes(
        problem: &'p SmallStrainTetProblem<'m, ResolvedJ2Plasticity>,
        mesh_frame: ContentHash,
        gauges: Vec<J2FatigueGauge>,
        critical: Vec<J2CriticalPlaneGauge>,
        config: J2FatigueConfig,
        cx: &Cx<'_>,
    ) -> Result<Self, J2FatigueError> {
        Self::new_with_fatemi_socie(
            problem,
            mesh_frame,
            gauges,
            critical,
            Vec::new(),
            config,
            cx,
        )
    }

    /// Observe fixed CMB, tensile SWT and/or shear Fatemi-Socie gauges on the
    /// same actual J2 trajectory. Every declared bank begins at the native
    /// zero-load state and participates in the existing atomic cycle lifecycle.
    pub fn new_with_fatemi_socie(
        problem: &'p SmallStrainTetProblem<'m, ResolvedJ2Plasticity>,
        mesh_frame: ContentHash,
        gauges: Vec<J2FatigueGauge>,
        critical: Vec<J2CriticalPlaneGauge>,
        fatemi_socie: Vec<J2FatemiSocieGauge>,
        config: J2FatigueConfig,
        cx: &Cx<'_>,
    ) -> Result<Self, J2FatigueError> {
        Self::new_with_variable_amplitude(
            problem,
            mesh_frame,
            gauges,
            critical,
            fatemi_socie,
            Vec::new(),
            config,
            cx,
        )
    }

    /// Observe variable-amplitude finite-plane Miner damage alongside all
    /// existing fixed-plane and representative-cycle fatigue banks.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_variable_amplitude(
        problem: &'p SmallStrainTetProblem<'m, ResolvedJ2Plasticity>,
        mesh_frame: ContentHash,
        gauges: Vec<J2FatigueGauge>,
        critical: Vec<J2CriticalPlaneGauge>,
        fatemi_socie: Vec<J2FatemiSocieGauge>,
        variable_amplitude: Vec<J2VariableAmplitudeGauge>,
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
        let mut observer = J2FatigueObserver::new(
            problem,
            mesh_frame,
            gauges,
            critical,
            fatemi_socie,
            variable_amplitude,
            config,
            cx,
        )?;
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

    /// One entry per critical gauge; an unfinished window has no life result.
    pub fn critical_plane_results(&self) -> Vec<Option<SwtCycleResult>> {
        self.observer.critical_plane_results()
    }
    pub fn critical_plane_curve(&self, gauge: usize) -> Option<&StrainLifeCurve> {
        self.observer.critical_plane_curve(gauge)
    }

    /// Separate shear-mode results; no axial-life or damage substitution.
    pub fn fatemi_socie_results(&self) -> Vec<Option<FatemiSocieCycleResult>> {
        self.observer.fatemi_socie_results()
    }
    pub fn fatemi_socie_calibration(&self, gauge: usize) -> Option<&FatemiSocieCalibration> {
        self.observer.fatemi_socie_calibration(gauge)
    }
    pub fn variable_amplitude_results(&self) -> Vec<Option<CriticalPlaneMinerResult>> {
        self.observer.variable_amplitude_results()
    }
    pub fn variable_amplitude_curve(&self, gauge: usize) -> Option<&StrainLifeCurve> {
        self.observer.variable_amplitude_curve(gauge)
    }

    /// Discard only the previous critical cycle window/result and begin at the
    /// CURRENT accepted native tensors, for example after physical warm-up.
    /// Mechanics, accepted increments and every fixed-gauge history are retained.
    /// Discarded windows are not silently credited as counted fatigue cycles.
    pub fn restart_critical_plane_cycle(&mut self, cx: &Cx<'_>) -> Result<(), J2FatigueError> {
        self.observer
            .restart_critical_plane_cycle(&self.solution, cx)
    }

    /// Accept exactly `count` planned increments as one transaction. Every
    /// original load factor is observed; callback boundaries close no cycles.
    /// Nonconvergence, cancellation, strain/calibration refusal or any gauge
    /// failure preserves the complete prior mechanics and fatigue histories.
    pub fn advance(&mut self, count: usize, cx: &Cx<'_>) -> Result<(), J2FatigueError> {
        self.advance_block(count, false, cx)
    }

    /// Advance native increments and require all critical windows to be closed
    /// representative cycles BEFORE committing mechanics or any gauge. Neither
    /// ordinary advance boundaries nor individual Newton steps close a cycle.
    /// Results describe the observed cycle, without repeated-history scaling.
    pub fn advance_cycle(&mut self, count: usize, cx: &Cx<'_>) -> Result<(), J2FatigueError> {
        self.observer.require_critical(cx)?;
        self.advance_block(count, true, cx)
    }

    fn advance_block(
        &mut self,
        count: usize,
        close_cycle: bool,
        cx: &Cx<'_>,
    ) -> Result<(), J2FatigueError> {
        poll(cx)?;
        if self.observer.finished {
            return Err(J2FatigueError::Finished);
        }
        if count == 0 || count > self.remaining_increments() {
            return Err(J2FatigueError::Invalid(
                "advance needs a nonzero count within the remaining load path",
            ));
        }
        self.observer.check_sample(cx)?;
        let end = self.accepted_increments + count;
        let mut solution = self.solution.clone();
        let checkpoint = self.observer.checkpoint();
        let proposal = (|| {
            for &factor in &self.problem.load_path[self.accepted_increments..end] {
                poll(cx)?;
                solution = with_path(self.problem, &[factor]).solve_from(
                    solution.displacement_m.as_flattened(),
                    &solution.element_states,
                    cx,
                )?;
                self.observer.observe(self.problem, &solution, cx)?;
            }
            if close_cycle {
                self.observer.close_critical_plane_cycle(cx)?;
            }
            poll(cx)
        })();
        if let Err(error) = proposal {
            self.observer.restore(checkpoint);
            return Err(error);
        }
        self.solution = solution;
        self.accepted_increments = end;
        Ok(())
    }

    /// Finish the entire declared observation once, after all planned targets.
    /// Fixed-plane endpoint half cycles and every critical cycle are finalized
    /// atomically. A critical bank must contain a resolved, closed single cycle.
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
        critical: Vec<J2CriticalPlaneGauge>,
        fatemi_socie: Vec<J2FatemiSocieGauge>,
        variable_amplitude: Vec<J2VariableAmplitudeGauge>,
        config: J2FatigueConfig,
        cx: &Cx<'_>,
    ) -> Result<Self, J2FatigueError> {
        poll(cx)?;
        if mesh_frame == ContentHash([0; 32]) {
            return Err(J2FatigueError::Invalid("declare a nonzero mesh frame"));
        }
        let fixed_residue = gauges
            .len()
            .checked_mul(config.history.maximum_residue)
            .ok_or(J2FatigueError::Budget)?;
        let variable_residue = variable_amplitude.iter().try_fold(0usize, |total, gauge| {
            gauge
                .normals_in_mesh_frame
                .len()
                .checked_mul(gauge.limits.rainflow.maximum_residue)
                .and_then(|n| total.checked_add(n))
                .ok_or(J2FatigueError::Budget)
        })?;
        if gauges
            .len()
            .checked_add(critical.len())
            .and_then(|n| n.checked_add(fatemi_socie.len()))
            .and_then(|n| n.checked_add(variable_amplitude.len()))
            .is_none_or(|count| count == 0 || count > config.maximum_gauges)
            || fixed_residue
                .checked_add(variable_residue)
                .is_none_or(|n| n > config.maximum_total_residue)
        {
            return Err(J2FatigueError::Budget);
        }
        problem.laws.validate_count(problem.tetrahedra.len())?;
        let mut prepared = Vec::with_capacity(gauges.len());
        for (index, gauge) in gauges.into_iter().enumerate() {
            poll(cx)?;
            bind_curve(
                problem,
                index,
                gauge.element,
                gauge.material_state_identity,
                &gauge.curve,
            )?;
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
        let fixed_count = prepared.len();
        let mut searches = Vec::with_capacity(critical.len());
        for (index, gauge) in critical.into_iter().enumerate() {
            poll(cx)?;
            bind_curve(
                problem,
                fixed_count + index,
                gauge.element,
                gauge.material_state_identity,
                &gauge.curve,
            )?;
            let counter = SwtCycle::new(
                gauge.curve,
                gauge.normals_in_mesh_frame,
                gauge.limits,
                gauge.closure,
            )
            .map_err(|error| J2FatigueError::CriticalPlane {
                gauge: index,
                error,
            })?;
            searches.push(CriticalGaugeState {
                element: gauge.element,
                counter,
            });
        }
        let shear_offset = fixed_count + searches.len();
        let mut shear_searches = Vec::with_capacity(fatemi_socie.len());
        for (index, gauge) in fatemi_socie.into_iter().enumerate() {
            poll(cx)?;
            bind_fatemi_socie(
                problem,
                shear_offset + index,
                gauge.element,
                gauge.material_state_identity,
                &gauge.calibration,
            )?;
            let counter = FatemiSocieCycle::new(
                gauge.calibration,
                gauge.normals_in_mesh_frame,
                gauge.limits,
                gauge.closure,
            )
            .map_err(|error| J2FatigueError::FatemiSocie {
                gauge: index,
                error,
            })?;
            shear_searches.push(FatemiSocieGaugeState {
                element: gauge.element,
                counter,
            });
        }
        let variable_offset = shear_offset + shear_searches.len();
        let mut variable_searches = Vec::with_capacity(variable_amplitude.len());
        for (index, gauge) in variable_amplitude.into_iter().enumerate() {
            poll(cx)?;
            bind_curve(
                problem,
                variable_offset + index,
                gauge.element,
                gauge.material_state_identity,
                &gauge.curve,
            )?;
            let counter = CriticalPlaneMiner::new(
                gauge.curve,
                gauge.mean_stress,
                gauge.normals_in_mesh_frame,
                gauge.limits,
            )
            .map_err(|error| J2FatigueError::VariableAmplitude {
                gauge: index,
                error,
            })?;
            variable_searches.push(VariableAmplitudeGaugeState {
                element: gauge.element,
                counter,
            });
        }
        poll(cx)?;
        Ok(Self {
            mesh_frame,
            gauges: prepared,
            critical: searches,
            fatemi_socie: shear_searches,
            variable_amplitude: variable_searches,
            samples: 0,
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
    pub(super) fn critical_plane_results(&self) -> Vec<Option<SwtCycleResult>> {
        self.critical
            .iter()
            .map(|g| g.counter.result().cloned())
            .collect()
    }
    pub(super) fn critical_plane_curve(&self, gauge: usize) -> Option<&StrainLifeCurve> {
        self.critical.get(gauge).map(|g| g.counter.curve())
    }
    pub(super) fn fatemi_socie_results(&self) -> Vec<Option<FatemiSocieCycleResult>> {
        self.fatemi_socie
            .iter()
            .map(|g| g.counter.result().copied())
            .collect()
    }
    pub(super) fn fatemi_socie_calibration(&self, gauge: usize) -> Option<&FatemiSocieCalibration> {
        self.fatemi_socie
            .get(gauge)
            .map(|g| g.counter.calibration())
    }
    pub(super) fn variable_amplitude_results(&self) -> Vec<Option<CriticalPlaneMinerResult>> {
        self.variable_amplitude
            .iter()
            .map(|g| g.counter.result().cloned())
            .collect()
    }
    pub(super) fn variable_amplitude_curve(&self, gauge: usize) -> Option<&StrainLifeCurve> {
        self.variable_amplitude.get(gauge).map(|g| g.counter.curve())
    }
    pub(super) fn checkpoint(&self) -> J2FatigueCheckpoint {
        J2FatigueCheckpoint {
            gauges: self.gauges.clone(),
            critical: self.critical.clone(),
            shear_prefixes: self
                .fatemi_socie
                .iter()
                .map(|g| (g.counter.sample_count(), g.counter.result().is_none()))
                .collect(),
            variable_amplitude: self.variable_amplitude.clone(),
            samples: self.samples,
            finished: self.finished,
        }
    }
    /// Internal accepted prefixes are always valid. Rollback cannot allocate,
    /// poll cancellation or discard a previously closed shear-cycle result.
    pub(super) fn restore(&mut self, checkpoint: J2FatigueCheckpoint) {
        self.gauges = checkpoint.gauges;
        self.critical = checkpoint.critical;
        self.variable_amplitude = checkpoint.variable_amplitude;
        self.samples = checkpoint.samples;
        self.finished = checkpoint.finished;
        for (gauge, (prefix, was_open)) in
            self.fatemi_socie.iter_mut().zip(checkpoint.shear_prefixes)
        {
            if was_open {
                gauge
                    .counter
                    .discard_after(prefix)
                    .expect("native observer rollback retains a previously accepted prefix");
            }
        }
    }

    pub(super) fn require_critical(&self, cx: &Cx<'_>) -> Result<(), J2FatigueError> {
        poll(cx)?;
        if self.finished {
            return Err(J2FatigueError::Finished);
        }
        if self.critical.is_empty() && self.fatemi_socie.is_empty() {
            return Err(J2FatigueError::Invalid(
                "no critical-plane cycle gauges were declared",
            ));
        }
        Ok(())
    }

    pub(super) fn restart_critical_plane_cycle(
        &mut self,
        solution: &SmallStrainTetSolution<J2State>,
        cx: &Cx<'_>,
    ) -> Result<(), J2FatigueError> {
        self.require_critical(cx)?;
        let mut critical = Vec::with_capacity(self.critical.len());
        for (index, gauge) in self.critical.iter().enumerate() {
            poll(cx)?;
            let counter = gauge
                .counter
                .fresh()
                .map_err(|error| J2FatigueError::CriticalPlane {
                    gauge: index,
                    error,
                })?;
            let mut next = CriticalGaugeState {
                element: gauge.element,
                counter,
            };
            push_critical(&mut next, index, solution, cx)?;
            critical.push(next);
        }
        let mut shear = Vec::with_capacity(self.fatemi_socie.len());
        for (index, gauge) in self.fatemi_socie.iter().enumerate() {
            poll(cx)?;
            let counter = gauge
                .counter
                .fresh()
                .map_err(|error| J2FatigueError::FatemiSocie {
                    gauge: index,
                    error,
                })?;
            let mut next = FatemiSocieGaugeState {
                element: gauge.element,
                counter,
            };
            push_fatemi_socie(&mut next, index, solution, cx)?;
            shear.push(next);
        }
        poll(cx)?;
        self.critical = critical;
        self.fatemi_socie = shear;
        Ok(())
    }

    pub(super) fn close_critical_plane_cycle(&mut self, cx: &Cx<'_>) -> Result<(), J2FatigueError> {
        self.require_critical(cx)?;
        let checkpoint = self.checkpoint();
        let proposal = (|| {
            for (index, gauge) in self.critical.iter_mut().enumerate() {
                gauge
                    .counter
                    .finish_with_cancel(&mut || cx.checkpoint().is_err())
                    .map_err(|error| J2FatigueError::CriticalPlane {
                        gauge: index,
                        error,
                    })?;
            }
            for (index, gauge) in self.fatemi_socie.iter_mut().enumerate() {
                gauge
                    .counter
                    .finish_with_cancel(&mut || cx.checkpoint().is_err())
                    .map_err(|error| J2FatigueError::FatemiSocie {
                        gauge: index,
                        error,
                    })?;
            }
            poll(cx)
        })();
        if proposal.is_err() {
            self.restore(checkpoint);
        }
        proposal
    }

    pub(super) fn check_sample(&self, cx: &Cx<'_>) -> Result<(), J2FatigueError> {
        poll(cx)?;
        if self.finished {
            return Err(J2FatigueError::Finished);
        }
        if self.samples >= self.maximum_samples {
            return Err(J2FatigueError::Budget);
        }
        if self.critical.iter().any(|g| g.counter.result().is_some())
            || self
                .fatemi_socie
                .iter()
                .any(|g| g.counter.result().is_some())
        {
            return Err(J2FatigueError::Invalid(
                "restart a closed critical-plane cycle before advancing",
            ));
        }
        Ok(())
    }
    /// Called inside the owning native block transaction, or before the initial
    /// state is published. The native owner restores its checkpoint on refusal.
    pub(super) fn observe(
        &mut self,
        problem: &SmallStrainTetProblem<'_, ResolvedJ2Plasticity>,
        solution: &SmallStrainTetSolution<J2State>,
        cx: &Cx<'_>,
    ) -> Result<(), J2FatigueError> {
        self.check_sample(cx)?;
        observe(problem, solution, &mut self.gauges, cx)?;
        for (index, gauge) in self.critical.iter_mut().enumerate() {
            push_critical(gauge, index, solution, cx)?;
        }
        for (index, gauge) in self.fatemi_socie.iter_mut().enumerate() {
            push_fatemi_socie(gauge, index, solution, cx)?;
        }
        for (index, gauge) in self.variable_amplitude.iter_mut().enumerate() {
            push_variable_amplitude(gauge, index, solution, cx)?;
        }
        poll(cx)?;
        self.samples += 1;
        Ok(())
    }
    pub(super) fn finish(&mut self, cx: &Cx<'_>) -> Result<Vec<StrainLifeSummary>, J2FatigueError> {
        poll(cx)?;
        if self.finished {
            return Ok(self.summaries());
        }
        let checkpoint = self.checkpoint();
        let proposal = (|| {
            for (index, gauge) in self.gauges.iter_mut().enumerate() {
                gauge
                    .counter
                    .finish(&mut || cx.checkpoint().is_err())
                    .map_err(|error| J2FatigueError::Gauge {
                        gauge: index,
                        error,
                    })?;
            }
            if !self.critical.is_empty() || !self.fatemi_socie.is_empty() {
                self.close_critical_plane_cycle(cx)?;
            }
            for (index, gauge) in self.variable_amplitude.iter_mut().enumerate() {
                gauge
                    .counter
                    .finish_with_cancel(&mut || cx.checkpoint().is_err())
                    .map_err(|error| J2FatigueError::VariableAmplitude {
                        gauge: index,
                        error,
                    })?;
            }
            poll(cx)
        })();
        if let Err(error) = proposal {
            self.restore(checkpoint);
            return Err(error);
        }
        self.finished = true;
        Ok(self.summaries())
    }
}

fn push_critical(
    gauge: &mut CriticalGaugeState,
    index: usize,
    solution: &SmallStrainTetSolution<J2State>,
    cx: &Cx<'_>,
) -> Result<(), J2FatigueError> {
    gauge
        .counter
        .push_with_cancel(
            TensorFatigueSample {
                stress_pa: solution.element_stress[gauge.element],
                total_mechanical_strain: solution.element_strain[gauge.element],
            },
            &mut || cx.checkpoint().is_err(),
        )
        .map_err(|error| J2FatigueError::CriticalPlane {
            gauge: index,
            error,
        })
}

fn push_variable_amplitude(
    gauge: &mut VariableAmplitudeGaugeState,
    index: usize,
    solution: &SmallStrainTetSolution<J2State>,
    cx: &Cx<'_>,
) -> Result<(), J2FatigueError> {
    gauge
        .counter
        .push_with_cancel(
            TensorFatigueSample {
                stress_pa: solution.element_stress[gauge.element],
                total_mechanical_strain: solution.element_strain[gauge.element],
            },
            &mut || cx.checkpoint().is_err(),
        )
        .map_err(|error| J2FatigueError::VariableAmplitude {
            gauge: index,
            error,
        })
}

fn push_fatemi_socie(
    gauge: &mut FatemiSocieGaugeState,
    index: usize,
    solution: &SmallStrainTetSolution<J2State>,
    cx: &Cx<'_>,
) -> Result<(), J2FatigueError> {
    gauge
        .counter
        .push_with_cancel(
            TensorFatigueSample {
                stress_pa: solution.element_stress[gauge.element],
                total_mechanical_strain: solution.element_strain[gauge.element],
            },
            &mut || cx.checkpoint().is_err(),
        )
        .map_err(|error| J2FatigueError::FatemiSocie {
            gauge: index,
            error,
        })
}

// SourceMismatch uses fixed, SWT, Fatemi-Socie, then variable-amplitude indexing.

fn bind_curve(
    problem: &SmallStrainTetProblem<'_, ResolvedJ2Plasticity>,
    index: usize,
    element: usize,
    identity: ContentHash,
    curve: &StrainLifeCurve,
) -> Result<(), J2FatigueError> {
    if element >= problem.tetrahedra.len() {
        return Err(J2FatigueError::Invalid(
            "gauge is outside the original tetrahedral field",
        ));
    }
    let law = problem.laws.at(element);
    let material = law.material();
    let scope = curve.scope();
    if identity != material.identity()
        || scope.material != material.material().chemistry
        || scope.process != material.material().process
    {
        return Err(J2FatigueError::SourceMismatch {
            gauge: index,
            what: "strain-life identity/material/process differs from the original J2 state",
        });
    }
    let temperature_axis = curve
        .source()
        .map_or("T", |source| source.request().temperature_axis.as_str());
    if material.state_coordinate(temperature_axis) != Some(scope.temperature_k)
        || curve.parameters().youngs_modulus_pa != law.law().elastic.youngs
    {
        return Err(J2FatigueError::SourceMismatch {
            gauge: index,
            what: "strain-life temperature/modulus differs from the frozen J2 law",
        });
    }
    if let Some(source) = curve.source() {
        let youngs =
            material
                .property(YOUNG_MODULUS_PROPERTY)
                .ok_or(J2FatigueError::SourceMismatch {
                    gauge: index,
                    what: "original J2 modulus receipt is missing or ambiguous",
                })?;
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
    Ok(())
}

fn bind_fatemi_socie(
    problem: &SmallStrainTetProblem<'_, ResolvedJ2Plasticity>,
    index: usize,
    element: usize,
    identity: ContentHash,
    calibration: &FatemiSocieCalibration,
) -> Result<(), J2FatigueError> {
    if element >= problem.tetrahedra.len() {
        return Err(J2FatigueError::Invalid(
            "gauge is outside the original tetrahedral field",
        ));
    }
    let law = problem.laws.at(element);
    let material = law.material();
    let scope = calibration.scope();
    let parameters = calibration.parameters();
    if identity != material.identity()
        || scope.material != material.material().chemistry
        || scope.process != material.material().process
    {
        return Err(J2FatigueError::SourceMismatch {
            gauge: index,
            what: "shear-life identity/material/process differs from the original J2 state",
        });
    }
    let temperature_axis = calibration
        .source()
        .map_or("T", |source| source.request().temperature_axis.as_str());
    if material.state_coordinate(temperature_axis) != Some(scope.temperature_k)
        || parameters.shear_modulus_pa != law.law().elastic.lame().1
        || parameters.initial_yield_stress_pa != law.law().yield_stress
    {
        return Err(J2FatigueError::SourceMismatch {
            gauge: index,
            what: "shear-life temperature, shear modulus or initial-yield reference differs from the frozen J2 law",
        });
    }
    if let Some(source) = calibration.source() {
        let answer = |name| {
            material
                .property(name)
                .map(|p| p.answer())
                .ok_or(J2FatigueError::SourceMismatch {
                    gauge: index,
                    what: "original J2 elastic/initial-yield receipt is missing or ambiguous",
                })
        };
        if !source.matches_material_state(
            material.card_identity(),
            material.material(),
            &material.query(),
            answer(YOUNG_MODULUS_PROPERTY)?,
            answer(POISSON_RATIO_PROPERTY)?,
            answer(YIELD_STRESS_PROPERTY)?,
        ) {
            return Err(J2FatigueError::SourceMismatch {
                gauge: index,
                what: "shear-life card, state, typed conditions or original elastic/initial-yield source differs from J2",
            });
        }
    }
    Ok(())
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
