//! Variable-amplitude strain-life on original CST/DKT shell material points.
//!
//! Gauges retain the actual facet, DKT surface point and thickness quadrature
//! point. Each accepted endpoint supplies the complete native plane-stress
//! strain, including transverse contraction and plastic history. Stress is the
//! original elastic map of the accepted elastic strain; observation performs
//! no return map or additional material update. Normals use the fixed local
//! facet frame, not global coordinates or a thickness-averaged resultant.
//!
//! Source-bound J2, Voce and finite-branch Chaboche laws retain their original
//! accepted histories, including every backstress branch. Backstress translates
//! the yield surface; it is not subtracted from the physical fatigue stress.
//! The material law and its external conditions remain fixed.
//! Thermal eigenstrain, changing-temperature calibration and time-dependent
//! constitutive laws are outside this adapter. A finite plane bank uses the
//! existing bounded rainflow/Miner owner and explicit uncorrected, Morrow or
//! SWT policy. Paired policies retain its coincident-extrema admission rules.
//! Separate Fatemi-Socie windows observe representative shear cycles with their
//! original torsional calibration and full-tensor closure. Optional shear Miner
//! accumulation counts each explicitly closed window once on every fixed plane;
//! no arbitrary shear reversal count or repeated exposure is inferred.
//! These empirical estimates do not change stiffness, deposit heat, create
//! cracks, establish cyclic stabilization or certify component service life.

mod loading;
mod shear;

pub use loading::ShellFatigueLoading;
pub use shear::{ShellFatemiSocieConfig, ShellFatemiSocieGauge};

use std::sync::Arc;

use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_material::fatigue::critical_plane::TensorFatigueSample;
use fs_material::fatigue::critical_plane::fatemi_socie::{
    FatemiSocieCalibration, FatemiSocieCycleResult, FatemiSocieError, FatemiSocieMinerSummary,
    FatemiSociePlaneDamage,
};
use fs_material::fatigue::critical_plane::miner::{
    CriticalPlaneMiner, CriticalPlaneMinerError, CriticalPlaneMinerLimits, CriticalPlaneMinerResult,
};
use fs_material::fatigue::cyclic::SourceBoundFatigueLaw;
use fs_material::fatigue::strain_life::{StrainLifeCurve, StrainLifeMeanStress, StrainLifeSummary};
use fs_material::plastic::source::ResolvedJ2Plasticity;
use fs_material::state_point::{ResolvedMaterialStatePoint, YOUNG_MODULUS_PROPERTY};
use fs_qty::semantic::{QuantityKind, ValueForm};
use fs_qty::{QuantitySpec, Temperature};

use super::dynamics::{PlasticShellReferenceInertia, ShellInertialEnergy};
use super::{
    Element, PlasticShellSections, SmallStrainShellProblem, SmallStrainShellSolution,
    SmallStrainTetError,
};
use shear::ShearBank;

/// One original material integration point and a declared finite plane bank.
#[derive(Debug, Clone)]
pub struct ShellFatigueGauge {
    pub element: usize,
    /// Original DKT midside quadrature index, in `0..3`.
    pub surface_point: usize,
    /// Index into the assigned section's original `points()` array.
    pub thickness_point: usize,
    pub material_state_identity: ContentHash,
    pub curve: Arc<StrainLifeCurve>,
    pub mean_stress: StrainLifeMeanStress,
    /// Directions in the native fixed local facet frame, including its normal.
    pub normals_in_facet_frame: Vec<[f64; 3]>,
    pub limits: CriticalPlaneMinerLimits,
}

/// Normal-strain bank bounds and the combined normal/shear gauge count.
/// Representative shear windows have separate aggregate sample/work limits.
#[derive(Debug, Clone, Copy)]
pub struct ShellFatigueConfig {
    /// Combined normal-strain and representative shear-cycle gauge count.
    pub maximum_gauges: usize,
    /// Sum of normal-strain plane counts times each plane's rainflow residue bound.
    pub maximum_total_residue: usize,
    /// Sum of normal-strain planes times samples for the complete declared load path,
    /// including its initial zero-load observation.
    pub maximum_total_projection_work: u64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum ShellFatigueError {
    Invalid(&'static str),
    Budget,
    SourceMismatch {
        /// Combined index: normal-strain gauges, then representative shear gauges.
        gauge: usize,
        what: &'static str,
    },
    Gauge {
        gauge: usize,
        error: CriticalPlaneMinerError,
    },
    FatemiSocie {
        gauge: usize,
        error: FatemiSocieError,
    },
    Mechanics(SmallStrainTetError),
    Cancelled,
    Finished,
    PathIncomplete,
}
impl core::fmt::Display for ShellFatigueError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native layered-shell fatigue refused: {self:?}")
    }
}
impl std::error::Error for ShellFatigueError {}
impl From<SmallStrainTetError> for ShellFatigueError {
    fn from(error: SmallStrainTetError) -> Self {
        match error {
            SmallStrainTetError::Cancelled => Self::Cancelled,
            error => Self::Mechanics(error),
        }
    }
}

#[derive(Debug, Clone)]
struct GaugeState {
    element: usize,
    surface_point: usize,
    thickness_point: usize,
    counter: CriticalPlaneMiner,
}

/// Resumable native shell mechanics with transactional fatigue observations.
///
/// Reference geometry is prepared once and shared by clones. Only the latest
/// native increment record, constitutive histories and bounded rainflow state
/// are retained. Every absolute factor in the original `load_path` is visited;
/// callers must resolve physical reversals with those targets. A failed native
/// solve or any later gauge refusal leaves the entire accepted step unchanged.
#[derive(Clone)]
pub struct ShellFatiguePath<'p, 'm, L: SourceBoundFatigueLaw = ResolvedJ2Plasticity> {
    problem: &'p SmallStrainShellProblem<'m, L>,
    elements: Arc<Vec<Element>>,
    solution: SmallStrainShellSolution<L::State>,
    velocity: Vec<[f64; 6]>,
    gauges: Vec<GaugeState>,
    samples: Vec<TensorFatigueSample>,
    shear: ShearBank,
    accepted_increments: usize,
    inertia: Option<&'p PlasticShellReferenceInertia<'m>>,
    elapsed_time_s: f64,
    finished: bool,
}
impl<L: SourceBoundFatigueLaw> core::fmt::Debug for ShellFatiguePath<'_, '_, L> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("ShellFatiguePath")
            .field("accepted_increments", &self.accepted_increments)
            .field("gauges", &self.gauges)
            .field("fatemi_socie_gauges", &self.shear.len())
            .field("elapsed_time_s", &self.elapsed_time_s)
            .field("finished", &self.finished)
            .finish_non_exhaustive()
    }
}

impl<'p, 'm, L: SourceBoundFatigueLaw> ShellFatiguePath<'p, 'm, L> {
    /// Observe the native virgin zero-load solution once, before the declared
    /// path. The shell retains its original material sources throughout; no
    /// external tensor history or unverified preload is imported.
    pub fn new(
        problem: &'p SmallStrainShellProblem<'m, L>,
        gauges: Vec<ShellFatigueGauge>,
        config: ShellFatigueConfig,
        cx: &Cx<'_>,
    ) -> Result<Self, ShellFatigueError> {
        Self::new_with_fatemi_socie(
            problem,
            gauges,
            Vec::new(),
            config,
            ShellFatemiSocieConfig::default(),
            cx,
        )
    }

    /// Observe normal-strain Miner damage, separately calibrated representative
    /// Fatemi-Socie cycles, or both on the same native shell trajectory. Shear
    /// windows start at the actual zero-load state; explicit restart selects a
    /// later accepted checkpoint after conditioning. Their limits apply to each
    /// declared window, independently of the whole-path normal-strain history.
    pub fn new_with_fatemi_socie(
        problem: &'p SmallStrainShellProblem<'m, L>,
        gauges: Vec<ShellFatigueGauge>,
        shear_gauges: Vec<ShellFatemiSocieGauge>,
        config: ShellFatigueConfig,
        shear_config: ShellFatemiSocieConfig,
        cx: &Cx<'_>,
    ) -> Result<Self, ShellFatigueError> {
        poll(cx)?;
        if gauges
            .len()
            .checked_add(shear_gauges.len())
            .is_none_or(|count| count == 0 || count > config.maximum_gauges)
            || !gauges.is_empty()
                && (config.maximum_total_residue == 0 || config.maximum_total_projection_work == 0)
        {
            return Err(ShellFatigueError::Budget);
        }
        let samples = u64::try_from(problem.load_path.len())
            .ok()
            .and_then(|count| count.checked_add(1))
            .ok_or(ShellFatigueError::Budget)?;
        let mut residue = 0usize;
        let mut work = 0u64;
        for gauge in &gauges {
            poll(cx)?;
            let planes = gauge.normals_in_facet_frame.len();
            residue = planes
                .checked_mul(gauge.limits.rainflow.maximum_residue)
                .and_then(|count| residue.checked_add(count))
                .ok_or(ShellFatigueError::Budget)?;
            let gauge_work = u64::try_from(planes)
                .ok()
                .and_then(|count| count.checked_mul(samples))
                .ok_or(ShellFatigueError::Budget)?;
            work = work
                .checked_add(gauge_work)
                .ok_or(ShellFatigueError::Budget)?;
            if residue > config.maximum_total_residue
                || work > config.maximum_total_projection_work
                || samples > gauge.limits.rainflow.maximum_samples
                || gauge_work > gauge.limits.maximum_projection_work
            {
                return Err(ShellFatigueError::Budget);
            }
        }
        let elements = Arc::new(problem.prepare(cx)?);
        let mut shear = ShearBank::new(problem, shear_gauges, gauges.len(), shear_config, cx)?;
        let mut bank = Vec::with_capacity(gauges.len());
        for (index, gauge) in gauges.into_iter().enumerate() {
            poll(cx)?;
            bind_curve(problem, index, &gauge)?;
            let counter = CriticalPlaneMiner::new(
                gauge.curve,
                gauge.mean_stress,
                gauge.normals_in_facet_frame,
                gauge.limits,
            )
            .map_err(|error| gauge_error(index, error))?;
            bank.push(GaugeState {
                element: gauge.element,
                surface_point: gauge.surface_point,
                thickness_point: gauge.thickness_point,
                counter,
            });
        }
        let zero = endpoint(problem, problem.prescribed, problem.nodal_loads, &[0.0]);
        let solution = zero.solve_prepared(&elements, None, None, cx)?;
        let initial_samples = observe(problem, &solution, &mut bank, cx)?;
        shear.observe(problem, &solution, false, cx)?;
        Ok(Self {
            problem,
            elements,
            solution,
            velocity: vec![[0.0; 6]; problem.mesh.nodes.len()],
            gauges: bank,
            samples: initial_samples,
            shear,
            accepted_increments: 0,
            inertia: None,
            elapsed_time_s: 0.0,
            finished: false,
        })
    }

    pub fn solution(&self) -> &SmallStrainShellSolution<L::State> {
        &self.solution
    }
    pub fn velocity(&self) -> &[[f64; 6]] {
        &self.velocity
    }
    /// Latest accepted full tensors, in gauge order and each native facet frame.
    pub fn samples(&self) -> &[TensorFatigueSample] {
        &self.samples
    }
    pub const fn next_step(&self) -> usize {
        self.accepted_increments
    }
    /// Includes the zero-load observation, even if the path starts at zero too.
    pub fn sample_count(&self) -> u64 {
        self.accepted_increments as u64 + 1
    }
    pub const fn elapsed_time_s(&self) -> f64 {
        self.elapsed_time_s
    }
    pub const fn is_finished(&self) -> bool {
        self.finished
    }
    pub fn is_complete(&self) -> bool {
        self.accepted_increments == self.problem.load_path.len()
    }
    pub fn curve(&self, gauge: usize) -> Option<&StrainLifeCurve> {
        self.gauges.get(gauge).map(|g| g.counter.curve())
    }
    pub fn material_source(&self, gauge: usize) -> Option<&ResolvedMaterialStatePoint> {
        let gauge = self.gauges.get(gauge)?;
        Some(point_law(self.problem, gauge.element, gauge.thickness_point).fatigue_source())
    }
    /// Closed-cycle damage accumulated so far, before unresolved endpoint halves.
    pub fn plane_summaries(&self, gauge: usize) -> Option<Vec<StrainLifeSummary>> {
        self.gauges.get(gauge).map(|g| g.counter.plane_summaries())
    }
    pub fn results(&self) -> Vec<Option<CriticalPlaneMinerResult>> {
        self.gauges
            .iter()
            .map(|g| g.counter.result().cloned())
            .collect()
    }

    /// Completed representative shear windows. Open windows have no life result;
    /// optional fixed-plane Miner totals are exposed separately.
    pub fn fatemi_socie_results(&self) -> Vec<Option<FatemiSocieCycleResult>> {
        self.shear.results()
    }

    /// Enable bounded Miner accumulation before observing a shear cycle,
    /// including after explicit conditioning. All gauges are admitted together;
    /// each accepted closed window counts once, without inferred repetitions.
    pub fn enable_fatemi_socie_miner(
        &mut self,
        maximum_cycles: u64,
        cx: &Cx<'_>,
    ) -> Result<(), ShellFatigueError> {
        self.require_shear(cx)?;
        self.shear.enable_miner(maximum_cycles, cx)
    }

    /// Greatest accumulated damage on a fixed plane in each opted-in gauge.
    pub fn fatemi_socie_miner_summaries(&self) -> Vec<Option<FatemiSocieMinerSummary>> {
        self.shear.miner_summaries()
    }

    /// Per-plane totals retain the original facet-frame normal ordering.
    pub fn fatemi_socie_plane_damage(&self, gauge: usize) -> Option<&[FatemiSociePlaneDamage]> {
        self.shear.miner_plane_damage(gauge)
    }

    /// Latest accepted physical tensors in shear-gauge order and facet frames.
    pub fn fatemi_socie_samples(&self) -> &[TensorFatigueSample] {
        self.shear.samples()
    }

    pub fn fatemi_socie_calibration(&self, gauge: usize) -> Option<&FatemiSocieCalibration> {
        self.shear.calibration(gauge)
    }

    pub fn fatemi_socie_material_source(
        &self,
        gauge: usize,
    ) -> Option<&ResolvedMaterialStatePoint> {
        self.shear.material_source(self.problem, gauge)
    }

    /// Begin every shear window at the CURRENT accepted native tensor state.
    /// This discards only the previous representative windows/results. Completed
    /// shear Miner totals, normal rainflow history, mechanics, velocity and
    /// physical time are preserved; discarded conditioning windows receive no
    /// implicit fatigue credit.
    pub fn restart_fatemi_socie_cycle(&mut self, cx: &Cx<'_>) -> Result<(), ShellFatigueError> {
        self.require_shear(cx)?;
        self.shear.restart(self.problem, &self.solution, cx)
    }

    /// Close the current representative shear windows together. Full stress
    /// and strain closure is required; one late refusal preserves every prior
    /// window and shear damage total. Each opted-in plane receives one cycle's
    /// damage, once. This changes no accepted mechanics or normal-strain history.
    pub fn close_fatemi_socie_cycle(
        &mut self,
        cx: &Cx<'_>,
    ) -> Result<Vec<FatemiSocieCycleResult>, ShellFatigueError> {
        self.require_shear(cx)?;
        self.shear.close(cx)
    }

    /// Accept one native quasistatic target. `false` means the path is complete.
    /// Static preload may precede dynamics at zero velocity; static continuation
    /// after a dynamic step refuses rather than silently discarding momentum.
    pub fn advance_next(&mut self, cx: &Cx<'_>) -> Result<bool, ShellFatigueError> {
        self.advance_next_with_fields(self.problem.prescribed, self.problem.nodal_loads, false, cx)
    }

    /// Accept the next static endpoint and close all representative shear
    /// cycles as one transaction. Failed closure publishes neither mechanics
    /// nor any normal/shear observation from the proposed endpoint.
    pub fn advance_next_cycle(&mut self, cx: &Cx<'_>) -> Result<bool, ShellFatigueError> {
        self.require_shear(cx)?;
        self.advance_next_with_fields(self.problem.prescribed, self.problem.nodal_loads, true, cx)
    }

    fn advance_next_with_fields(
        &mut self,
        prescribed: &[(usize, f64)],
        nodal_loads: &[(usize, f64)],
        close_shear: bool,
        cx: &Cx<'_>,
    ) -> Result<bool, ShellFatigueError> {
        self.admit_advance(cx)?;
        if self.inertia.is_some() {
            return Err(ShellFatigueError::Invalid(
                "quasistatic continuation cannot discard accepted shell motion",
            ));
        }
        let Some(&factor) = self.problem.load_path.get(self.accepted_increments) else {
            return Ok(false);
        };
        let factors = [factor];
        let problem = endpoint(self.problem, prescribed, nodal_loads, &factors);
        let solution = problem.solve_prepared(
            &self.elements,
            Some((&self.solution.coordinates, &self.solution.section_states)),
            None,
            cx,
        )?;
        let mut gauges = self.gauges.clone();
        let samples = observe(self.problem, &solution, &mut gauges, cx)?;
        self.shear
            .observe(self.problem, &solution, close_shear, cx)?;
        self.solution = solution;
        self.gauges = gauges;
        self.samples = samples;
        self.accepted_increments += 1;
        Ok(true)
    }

    /// Advance the next absolute target with original backward-Euler inertia.
    /// Both prescribed motion and nodal loads are scaled by its factor. Retain
    /// the same physical reference-inertia owner across every dynamic step.
    /// Returned kinetic and numerical-loss accounts are the native solver's;
    /// fatigue adds no force, material heat or numerical damping.
    pub fn advance_next_backward_euler(
        &mut self,
        reference: &'p PlasticShellReferenceInertia<'m>,
        dt_s: f64,
        cx: &Cx<'_>,
    ) -> Result<Option<ShellInertialEnergy>, ShellFatigueError> {
        self.advance_next_backward_euler_with_fields(
            reference,
            dt_s,
            self.problem.prescribed,
            self.problem.nodal_loads,
            false,
            cx,
        )
    }

    /// Accept the next physical endpoint and close all representative shear
    /// cycles atomically, including velocity, time and original inertia receipt.
    pub fn advance_next_backward_euler_cycle(
        &mut self,
        reference: &'p PlasticShellReferenceInertia<'m>,
        dt_s: f64,
        cx: &Cx<'_>,
    ) -> Result<Option<ShellInertialEnergy>, ShellFatigueError> {
        self.require_shear(cx)?;
        self.advance_next_backward_euler_with_fields(
            reference,
            dt_s,
            self.problem.prescribed,
            self.problem.nodal_loads,
            true,
            cx,
        )
    }

    fn advance_next_backward_euler_with_fields(
        &mut self,
        reference: &'p PlasticShellReferenceInertia<'m>,
        dt_s: f64,
        prescribed: &[(usize, f64)],
        nodal_loads: &[(usize, f64)],
        close_shear: bool,
        cx: &Cx<'_>,
    ) -> Result<Option<ShellInertialEnergy>, ShellFatigueError> {
        self.admit_advance(cx)?;
        if self
            .inertia
            .is_some_and(|old| !std::ptr::eq(old, reference))
        {
            return Err(ShellFatigueError::Invalid(
                "dynamic fatigue must retain the original shell reference inertia",
            ));
        }
        let Some(&factor) = self.problem.load_path.get(self.accepted_increments) else {
            return Ok(None);
        };
        let time = self.elapsed_time_s + dt_s;
        if !dt_s.is_finite() || dt_s <= 0.0 || !time.is_finite() || time <= self.elapsed_time_s {
            return Err(ShellFatigueError::Invalid(
                "dynamic fatigue needs a positive representable physical interval",
            ));
        }
        let prescribed: Vec<_> = prescribed
            .iter()
            .map(|&(dof, value)| (dof, factor * value))
            .collect();
        let loads: Vec<_> = nodal_loads
            .iter()
            .map(|&(dof, value)| (dof, factor * value))
            .collect();
        let problem = endpoint(self.problem, &prescribed, &loads, &[1.0]);
        let mut step = problem.advance_backward_euler_prepared(
            &self.elements,
            Some(&self.solution),
            &self.velocity,
            reference,
            dt_s,
            cx,
        )?;
        let mut gauges = self.gauges.clone();
        let samples = observe(self.problem, &step.solution, &mut gauges, cx)?;
        // The native bridge uses a factor-one endpoint with scaled fields;
        // retain the caller's actual absolute path factor in the public record.
        step.solution.increments[0].load_factor = factor;
        self.shear
            .observe(self.problem, &step.solution, close_shear, cx)?;
        self.solution = step.solution;
        self.velocity = step.velocity;
        self.gauges = gauges;
        self.samples = samples;
        self.inertia = Some(reference);
        self.elapsed_time_s = time;
        self.accepted_increments += 1;
        Ok(Some(step.inertia))
    }

    /// Read-only hypothetical end of the normal-strain history, including its
    /// residual half cycles. Does not consume reversals, complete the path, or
    /// close/evaluate any representative Fatemi-Socie window.
    pub fn preview_finish(
        &self,
        cx: &Cx<'_>,
    ) -> Result<Vec<CriticalPlaneMinerResult>, ShellFatigueError> {
        poll(cx)?;
        let mut results = Vec::with_capacity(self.gauges.len());
        for (index, gauge) in self.gauges.iter().enumerate() {
            results.push(
                gauge
                    .counter
                    .preview_finish_with_cancel(&mut || cx.checkpoint().is_err())
                    .map_err(|error| gauge_error(index, error))?,
            );
        }
        poll(cx)?;
        Ok(results)
    }

    /// Finish every normal-strain bank and representative shear cycle atomically
    /// after all declared targets are accepted. A shear closure/calibration
    /// refusal preserves the prior normal endpoint state and every shear window.
    pub fn finish(
        &mut self,
        cx: &Cx<'_>,
    ) -> Result<Vec<CriticalPlaneMinerResult>, ShellFatigueError> {
        poll(cx)?;
        if !self.is_complete() {
            return Err(ShellFatigueError::PathIncomplete);
        }
        let mut gauges = self.gauges.clone();
        let mut results = Vec::with_capacity(gauges.len());
        for (index, gauge) in gauges.iter_mut().enumerate() {
            results.push(
                gauge
                    .counter
                    .finish_with_cancel(&mut || cx.checkpoint().is_err())
                    .map_err(|error| gauge_error(index, error))?,
            );
        }
        self.shear.close(cx)?;
        self.gauges = gauges;
        self.finished = true;
        Ok(results)
    }

    fn admit_advance(&self, cx: &Cx<'_>) -> Result<(), ShellFatigueError> {
        poll(cx)?;
        if self.finished {
            return Err(ShellFatigueError::Finished);
        }
        self.shear.admit_advance()?;
        Ok(())
    }

    fn require_shear(&self, cx: &Cx<'_>) -> Result<(), ShellFatigueError> {
        poll(cx)?;
        if self.finished {
            return Err(ShellFatigueError::Finished);
        }
        if self.shear.len() == 0 {
            return Err(ShellFatigueError::Invalid(
                "no representative Fatemi-Socie gauges were declared",
            ));
        }
        Ok(())
    }
}

fn poll(cx: &Cx<'_>) -> Result<(), ShellFatigueError> {
    cx.checkpoint().map_err(|_| ShellFatigueError::Cancelled)
}
fn gauge_error(gauge: usize, error: CriticalPlaneMinerError) -> ShellFatigueError {
    match error {
        CriticalPlaneMinerError::Cancelled => ShellFatigueError::Cancelled,
        error => ShellFatigueError::Gauge { gauge, error },
    }
}

fn endpoint<'a, L: SourceBoundFatigueLaw>(
    problem: &'a SmallStrainShellProblem<'_, L>,
    prescribed: &'a [(usize, f64)],
    nodal_loads: &'a [(usize, f64)],
    load_path: &'a [f64],
) -> SmallStrainShellProblem<'a, L> {
    SmallStrainShellProblem {
        mesh: problem.mesh,
        sections: match &problem.sections {
            PlasticShellSections::Uniform(section) => PlasticShellSections::Uniform(section),
            PlasticShellSections::PerElement(sections) => {
                PlasticShellSections::PerElement(sections)
            }
            PlasticShellSections::PerElementBorrowed(sections) => {
                PlasticShellSections::PerElementBorrowed(sections)
            }
        },
        coordinate_length_m: problem.coordinate_length_m,
        rotation_limit_rad: problem.rotation_limit_rad,
        drilling_alpha: problem.drilling_alpha,
        prescribed,
        nodal_loads,
        load_path,
        budget: problem.budget,
        settings: problem.settings,
    }
}

fn point_law<'a, L: SourceBoundFatigueLaw>(
    problem: &'a SmallStrainShellProblem<'_, L>,
    element: usize,
    thickness_point: usize,
) -> &'a L {
    let section = problem.sections.at(element);
    &section.layers()[section.points()[thickness_point].layer].law
}

fn accepted_sample<L: SourceBoundFatigueLaw>(
    problem: &SmallStrainShellProblem<'_, L>,
    solution: &SmallStrainShellSolution<L::State>,
    element: usize,
    surface_point: usize,
    thickness_point: usize,
) -> Result<TensorFatigueSample, ShellFatigueError> {
    let point = &solution.section_states[3 * element + surface_point].points[thickness_point];
    point_law(problem, element, thickness_point)
        .accepted_fatigue_sample(&point.strain, &point.history)
        .map_err(|error| {
            SmallStrainTetError::InvalidInput {
                what: format!(
                    "shell fatigue element {element}, surface point {surface_point}, thickness point {thickness_point}: {error}"
                ),
            }
            .into()
        })
}

fn observe<L: SourceBoundFatigueLaw>(
    problem: &SmallStrainShellProblem<'_, L>,
    solution: &SmallStrainShellSolution<L::State>,
    gauges: &mut [GaugeState],
    cx: &Cx<'_>,
) -> Result<Vec<TensorFatigueSample>, ShellFatigueError> {
    let mut samples = Vec::with_capacity(gauges.len());
    for (index, gauge) in gauges.iter_mut().enumerate() {
        poll(cx)?;
        let sample = accepted_sample(
            problem,
            solution,
            gauge.element,
            gauge.surface_point,
            gauge.thickness_point,
        )?;
        gauge
            .counter
            .push_with_cancel(sample, &mut || cx.checkpoint().is_err())
            .map_err(|error| gauge_error(index, error))?;
        samples.push(sample);
    }
    Ok(samples)
}

fn bind_curve<L: SourceBoundFatigueLaw>(
    problem: &SmallStrainShellProblem<'_, L>,
    index: usize,
    gauge: &ShellFatigueGauge,
) -> Result<(), ShellFatigueError> {
    if gauge.element >= problem.mesh.tris.len() || gauge.surface_point >= 3 {
        return Err(ShellFatigueError::Invalid(
            "gauge is outside the original shell surface quadrature",
        ));
    }
    let section = problem.sections.at(gauge.element);
    let point = section
        .points()
        .get(gauge.thickness_point)
        .ok_or(ShellFatigueError::Invalid(
            "gauge is outside the original shell thickness quadrature",
        ))?;
    let law = &section.layers()[point.layer].law;
    let material = law.fatigue_source();
    let scope = gauge.curve.scope();
    let mismatch = |what| ShellFatigueError::SourceMismatch { gauge: index, what };
    if gauge.material_state_identity != material.identity()
        || scope.material != material.material().chemistry
        || scope.process != material.material().process
    {
        return Err(mismatch(
            "CMB identity/material/process differs from the native shell layer",
        ));
    }
    let axis = gauge
        .curve
        .source()
        .map_or("T", |source| source.request().temperature_axis.as_str());
    let absolute_temperature = material
        .axis_quantities()
        .get(axis)
        .is_some_and(|quantity| {
            *quantity == QuantitySpec::dimensional(Temperature::DIMS)
                || quantity.semantic_type().is_some_and(|semantic| {
                    semantic.kind() == QuantityKind::AbsoluteTemperature
                        && semantic.form() == ValueForm::Static
                })
        });
    if !absolute_temperature
        || material.state_coordinate(axis) != Some(scope.temperature_k)
        || gauge.curve.parameters().youngs_modulus_pa != law.fatigue_backbone().elastic.youngs
    {
        return Err(mismatch(
            "CMB absolute temperature/modulus differs from the frozen shell law",
        ));
    }
    if let Some(source) = gauge.curve.source() {
        let youngs = material
            .property(YOUNG_MODULUS_PROPERTY)
            .ok_or_else(|| mismatch("original shell modulus receipt is missing or ambiguous"))?;
        if !source.matches_material_state(
            material.card_identity(),
            material.material(),
            &material.query(),
            youngs.answer(),
        ) {
            return Err(mismatch(
                "CMB card, state, external conditions or modulus claim differs from the shell source",
            ));
        }
    }
    Ok(())
}
