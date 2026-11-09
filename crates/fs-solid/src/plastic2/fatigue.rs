//! Variable-amplitude fatigue on accepted plane-stress/plane-strain histories.
//!
//! Each original P1/Q1 integration point retains its complete three-dimensional
//! constitutive strain and physical stress. Plane stress therefore retains the
//! solved thickness contraction; plane strain retains constrained transverse
//! stress, and B-bar uses its actual assumed material strain. Neither principal
//! nor unsigned equivalent stress replaces this signed tensor history.
//!
//! This observer calls the original nonlinear solve at every load entry. It
//! never performs a second return map, changes stiffness, or infers material
//! heat from Miner damage. J2, Voce and Chaboche histories remain native. Source
//! conditions stay fixed; thermal/rate laws, remeshing and service-life
//! validation require separate models and evidence. Optional Fatemi-Socie
//! gauges share the same solve and accepted tensors, while retaining explicit
//! closed representative-cycle windows alongside the live normal-strain Miner
//! histories. Optional shear Miner accumulation counts each explicitly closed
//! window once on every fixed plane. No repetition count or stabilized response
//! is inferred.

pub mod shear;
pub use shear::{PlaneFatemiSocieGauge, PlaneShearFatigueBudget};

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

use super::{Element, SmallStrainPlaneError, SmallStrainPlaneProblem, SmallStrainPlaneSolution};
use crate::plastic3::SmallStrainLawField;

/// One native quadrature point and a caller-declared finite plane bank.
#[derive(Debug, Clone)]
pub struct PlaneFatigueGauge {
    pub element: usize,
    /// Local index in the original `mesh2::quad_points` order.
    pub quadrature_point: usize,
    pub material_state_identity: ContentHash,
    pub curve: Arc<StrainLifeCurve>,
    pub mean_stress: StrainLifeMeanStress,
    /// Fixed mesh x/y axes, with z the out-of-plane direction.
    pub normals_in_mesh_frame: Vec<[f64; 3]>,
    pub limits: CriticalPlaneMinerLimits,
}

/// Normal-strain history bounds and combined gauge count, including the real
/// initial zero-load sample. Shear windows have separate aggregate budgets.
#[derive(Debug, Clone, Copy)]
pub struct PlaneFatigueConfig {
    /// Combined normal-strain and shear representative-cycle gauges.
    pub maximum_gauges: usize,
    /// Sum of each gauge's planes times its per-plane rainflow residue limit.
    pub maximum_total_residue: usize,
    /// Sum of normal-strain plane counts times the declared path's sample count.
    pub maximum_total_projection_work: u64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum PlaneFatigueError {
    Invalid(&'static str),
    Budget,
    SourceMismatch {
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
    Mechanics(SmallStrainPlaneError),
    Cancelled,
    Finished,
    PathIncomplete,
}
impl core::fmt::Display for PlaneFatigueError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native plane-continuum fatigue refused: {self:?}")
    }
}
impl std::error::Error for PlaneFatigueError {}
impl From<SmallStrainPlaneError> for PlaneFatigueError {
    fn from(error: SmallStrainPlaneError) -> Self {
        match error {
            SmallStrainPlaneError::Cancelled => Self::Cancelled,
            error => Self::Mechanics(error),
        }
    }
}

#[derive(Debug, Clone)]
struct GaugeState {
    element: usize,
    point: usize,
    counter: CriticalPlaneMiner,
}

/// Resumable native mechanics and bounded fatigue, published transactionally.
/// Keep the latest solution, bounded normal residues and explicitly bounded
/// representative shear windows. Clones share immutable original quadrature
/// geometry and calibration, and retain complete open observation histories.
#[derive(Debug, Clone)]
pub struct PlaneFatiguePath<'p, 'm, L: SourceBoundFatigueLaw = ResolvedJ2Plasticity> {
    problem: &'p SmallStrainPlaneProblem<'m, L>,
    elements: Arc<Vec<Element>>,
    solution: SmallStrainPlaneSolution<L::State>,
    gauges: Vec<GaugeState>,
    samples: Vec<TensorFatigueSample>,
    fatemi_socie: Vec<shear::GaugeState>,
    shear_samples: Vec<TensorFatigueSample>,
    accepted_increments: usize,
    finished: bool,
}

impl<'p, 'm, L: SourceBoundFatigueLaw> PlaneFatiguePath<'p, 'm, L> {
    /// Begin at the original virgin zero-load equilibrium. Every subsequent
    /// absolute load-path entry must resolve the physical reversals of interest.
    pub fn new(
        problem: &'p SmallStrainPlaneProblem<'m, L>,
        gauges: Vec<PlaneFatigueGauge>,
        config: PlaneFatigueConfig,
        cx: &Cx<'_>,
    ) -> Result<Self, PlaneFatigueError> {
        Self::new_with_shear(
            problem,
            gauges,
            Vec::new(),
            config,
            PlaneShearFatigueBudget::default(),
            cx,
        )
    }

    /// Observe normal-strain Miner and/or a separately calibrated shear-cycle
    /// bank on the same native solve. Shear windows may be restarted after
    /// conditioning without resetting the normal-strain cycle history.
    pub fn new_with_shear(
        problem: &'p SmallStrainPlaneProblem<'m, L>,
        gauges: Vec<PlaneFatigueGauge>,
        fatemi_socie: Vec<PlaneFatemiSocieGauge>,
        config: PlaneFatigueConfig,
        shear_budget: PlaneShearFatigueBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, PlaneFatigueError> {
        poll(cx)?;
        let gauge_count = gauges
            .len()
            .checked_add(fatemi_socie.len())
            .ok_or(PlaneFatigueError::Budget)?;
        if gauge_count == 0 || gauge_count > config.maximum_gauges {
            return Err(PlaneFatigueError::Budget);
        }
        shear::admit_budget(&fatemi_socie, shear_budget, cx)?;
        let count = u64::try_from(problem.load_path.len())
            .ok()
            .and_then(|n| n.checked_add(1))
            .ok_or(PlaneFatigueError::Budget)?;
        let (mut residue, mut work) = (0usize, 0u64);
        for gauge in &gauges {
            poll(cx)?;
            let planes = gauge.normals_in_mesh_frame.len();
            residue = planes
                .checked_mul(gauge.limits.rainflow.maximum_residue)
                .and_then(|n| residue.checked_add(n))
                .ok_or(PlaneFatigueError::Budget)?;
            let gauge_work = u64::try_from(planes)
                .ok()
                .and_then(|n| n.checked_mul(count))
                .ok_or(PlaneFatigueError::Budget)?;
            work = work
                .checked_add(gauge_work)
                .ok_or(PlaneFatigueError::Budget)?;
            if residue > config.maximum_total_residue
                || work > config.maximum_total_projection_work
                || count > gauge.limits.rainflow.maximum_samples
                || gauge_work > gauge.limits.maximum_projection_work
            {
                return Err(PlaneFatigueError::Budget);
            }
        }
        let elements = Arc::new(problem.prepare(cx)?);
        let mut bank = Vec::with_capacity(gauges.len());
        for (index, gauge) in gauges.into_iter().enumerate() {
            poll(cx)?;
            if elements
                .get(gauge.element)
                .is_none_or(|element| gauge.quadrature_point >= element.points.len())
            {
                return Err(PlaneFatigueError::Invalid(
                    "gauge is outside the original plane quadrature",
                ));
            }
            bind_curve(
                problem.law(gauge.element, gauge.quadrature_point),
                index,
                &gauge,
            )?;
            bank.push(GaugeState {
                element: gauge.element,
                point: gauge.quadrature_point,
                counter: CriticalPlaneMiner::new(
                    gauge.curve,
                    gauge.mean_stress,
                    gauge.normals_in_mesh_frame,
                    gauge.limits,
                )
                .map_err(|error| gauge_error(index, error))?,
            });
        }
        let mut shear = shear::prepare(problem, &elements, fatemi_socie, bank.len(), cx)?;
        let solution = endpoint(problem, &[0.0]).solve_prepared(&elements, None, cx)?;
        let samples = observe(problem, &solution, &mut bank, cx)?;
        let shear_samples = shear::observe(problem, &solution, &mut shear, cx)?;
        poll(cx)?;
        Ok(Self {
            problem,
            elements,
            solution,
            gauges: bank,
            samples,
            fatemi_socie: shear,
            shear_samples,
            accepted_increments: 0,
            finished: false,
        })
    }

    pub fn solution(&self) -> &SmallStrainPlaneSolution<L::State> {
        &self.solution
    }
    /// Current accepted tensors in the normal-strain gauge order.
    pub fn samples(&self) -> &[TensorFatigueSample] {
        &self.samples
    }
    pub const fn accepted_increments(&self) -> usize {
        self.accepted_increments
    }
    pub fn remaining_increments(&self) -> usize {
        self.problem.load_path.len() - self.accepted_increments
    }
    pub fn sample_count(&self) -> u64 {
        self.accepted_increments as u64 + 1
    }
    pub const fn is_finished(&self) -> bool {
        self.finished
    }
    pub fn curve(&self, gauge: usize) -> Option<&StrainLifeCurve> {
        self.gauges.get(gauge).map(|g| g.counter.curve())
    }
    pub fn material_source(&self, gauge: usize) -> Option<&ResolvedMaterialStatePoint> {
        self.gauges
            .get(gauge)
            .map(|g| self.problem.law(g.element, g.point).fatigue_source())
    }
    /// Damage on each declared plane before unresolved endpoint half cycles.
    pub fn plane_summaries(&self, gauge: usize) -> Option<Vec<StrainLifeSummary>> {
        self.gauges.get(gauge).map(|g| g.counter.plane_summaries())
    }
    pub fn results(&self) -> Vec<Option<CriticalPlaneMinerResult>> {
        self.gauges
            .iter()
            .map(|g| g.counter.result().cloned())
            .collect()
    }

    /// Current accepted tensors in the declared shear-gauge order.
    pub fn fatemi_socie_samples(&self) -> &[TensorFatigueSample] {
        &self.shear_samples
    }
    pub fn fatemi_socie_sample_count(&self, gauge: usize) -> Option<u64> {
        self.fatemi_socie
            .get(gauge)
            .map(|g| g.counter.sample_count())
    }
    pub fn fatemi_socie_results(&self) -> Vec<Option<FatemiSocieCycleResult>> {
        self.fatemi_socie
            .iter()
            .map(|g| g.counter.result().copied())
            .collect()
    }

    /// Select bounded Miner accumulation for explicitly closed shear cycles.
    /// Enable before observing a cycle, including after an explicit conditioning
    /// restart. Every gauge is admitted together; failed admission or cancellation
    /// changes neither mechanics nor any existing fatigue state.
    pub fn enable_fatemi_socie_miner(
        &mut self,
        maximum_cycles: u64,
        cx: &Cx<'_>,
    ) -> Result<(), PlaneFatigueError> {
        poll(cx)?;
        if self.finished {
            return Err(PlaneFatigueError::Finished);
        }
        if self.fatemi_socie.is_empty() {
            return Err(PlaneFatigueError::Invalid(
                "no shear-cycle gauges were declared",
            ));
        }
        let mut gauges = Vec::new();
        gauges
            .try_reserve_exact(self.fatemi_socie.len())
            .map_err(|_| PlaneFatigueError::Budget)?;
        for (index, gauge) in self.fatemi_socie.iter().enumerate() {
            poll(cx)?;
            if gauge.counter.sample_count() > 1 {
                return Err(PlaneFatigueError::Invalid(
                    "enable FS/Miner at the first observation of a fresh cycle window",
                ));
            }
            let mut next = gauge.clone();
            next.counter.enable_miner(maximum_cycles).map_err(|error| {
                PlaneFatigueError::FatemiSocie {
                    gauge: index,
                    error,
                }
            })?;
            gauges.push(next);
        }
        poll(cx)?;
        self.fatemi_socie = gauges;
        Ok(())
    }

    /// Greatest accumulated damage on a fixed plane in each opted-in gauge.
    /// Unfinished or discarded windows contribute no cycle exposure.
    pub fn fatemi_socie_miner_summaries(&self) -> Vec<Option<FatemiSocieMinerSummary>> {
        self.fatemi_socie
            .iter()
            .map(|g| g.counter.miner_summary())
            .collect()
    }

    /// Per-plane totals retain the original mesh-frame normal ordering.
    pub fn fatemi_socie_plane_damage(&self, gauge: usize) -> Option<&[FatemiSociePlaneDamage]> {
        self.fatemi_socie.get(gauge)?.counter.miner_plane_damage()
    }
    pub fn fatemi_socie_calibration(&self, gauge: usize) -> Option<&FatemiSocieCalibration> {
        self.fatemi_socie
            .get(gauge)
            .map(|g| g.counter.calibration())
    }
    pub fn fatemi_socie_material_source(
        &self,
        gauge: usize,
    ) -> Option<&ResolvedMaterialStatePoint> {
        self.fatemi_socie
            .get(gauge)
            .map(|g| self.problem.law(g.element, g.point).fatigue_source())
    }

    /// Begin a new declared representative cycle at the actual current native
    /// state. Conditioning is explicit; this does not reset mechanical history,
    /// normal-strain or completed shear Miner totals, or the accepted load-path
    /// position. Discarded open windows contribute no shear damage.
    pub fn restart_fatemi_socie_cycle(&mut self, cx: &Cx<'_>) -> Result<(), PlaneFatigueError> {
        poll(cx)?;
        if self.finished {
            return Err(PlaneFatigueError::Finished);
        }
        if self.fatemi_socie.is_empty() {
            return Err(PlaneFatigueError::Invalid(
                "no shear-cycle gauges were declared",
            ));
        }
        let (gauges, samples) =
            shear::restart(self.problem, &self.solution, &self.fatemi_socie, cx)?;
        self.fatemi_socie = gauges;
        self.shear_samples = samples;
        Ok(())
    }

    /// Finalize the current closed shear window at any accepted endpoint.
    /// Closure of all six stress and strain components is required. The result
    /// predicts the declared representative cycle's life. Opted-in Miner totals
    /// separately count this window once on each fixed plane; finalizing the
    /// same window again adds no exposure.
    /// Further mechanical advancement requires an explicit cycle restart.
    pub fn finish_fatemi_socie_cycle(
        &mut self,
        cx: &Cx<'_>,
    ) -> Result<Vec<FatemiSocieCycleResult>, PlaneFatigueError> {
        if self.fatemi_socie.is_empty() {
            return Err(PlaneFatigueError::Invalid(
                "no shear-cycle gauges were declared",
            ));
        }
        shear::finish_cycles(&mut self.fatemi_socie, cx)
    }

    /// Accept a chronological block. All intermediate endpoints are observed;
    /// a late mechanics/gauge refusal rolls back this entire requested block.
    /// Only current mechanics and bounded normal residues are cloned. Shear
    /// windows restore saved prefixes without copying their retained tensors.
    pub fn advance(&mut self, count: usize, cx: &Cx<'_>) -> Result<(), PlaneFatigueError> {
        poll(cx)?;
        if self.finished {
            return Err(PlaneFatigueError::Finished);
        }
        if count > self.remaining_increments() {
            return Err(PlaneFatigueError::Invalid(
                "requested block exceeds the declared load path",
            ));
        }
        if count == 0 {
            return Ok(());
        }
        if self
            .fatemi_socie
            .iter()
            .any(|g| g.counter.result().is_some())
        {
            return Err(PlaneFatigueError::Invalid(
                "restart a closed shear cycle before advancing",
            ));
        }
        let mut solution = self.solution.clone();
        let mut gauges = self.gauges.clone();
        let mut samples = self.samples.clone();
        let mut shear_samples = self.shear_samples.clone();
        let checkpoint = shear::checkpoint(&self.fatemi_socie);
        let attempt = (|| -> Result<(), PlaneFatigueError> {
            for offset in 0..count {
                poll(cx)?;
                let factor = self.problem.load_path[self.accepted_increments + offset];
                let next = endpoint(self.problem, &[factor]).solve_prepared(
                    &self.elements,
                    Some((
                        solution.displacement_m.as_flattened(),
                        &solution.quadrature_states,
                    )),
                    cx,
                )?;
                samples = observe(self.problem, &next, &mut gauges, cx)?;
                shear_samples = shear::observe(self.problem, &next, &mut self.fatemi_socie, cx)?;
                solution = next;
            }
            poll(cx)
        })();
        if let Err(error) = attempt {
            shear::restore(&mut self.fatemi_socie, &checkpoint);
            return Err(error);
        }
        self.solution = solution;
        self.gauges = gauges;
        self.samples = samples;
        self.shear_samples = shear_samples;
        self.accepted_increments += count;
        Ok(())
    }

    /// Read-only damage if the accepted prefix ended now. The native solution,
    /// continuation, unresolved residues and sample allowances are unchanged.
    pub fn preview_finish(
        &self,
        cx: &Cx<'_>,
    ) -> Result<Vec<CriticalPlaneMinerResult>, PlaneFatigueError> {
        poll(cx)?;
        self.gauges
            .iter()
            .enumerate()
            .map(|(index, gauge)| {
                gauge
                    .counter
                    .preview_finish_with_cancel(&mut || cx.checkpoint().is_err())
                    .map_err(|error| gauge_error(index, error))
            })
            .collect()
    }

    /// Close normal endpoint half cycles and any declared shear cycle only at
    /// the end of the load path. Repeated finalization is idempotent; a failed
    /// gauge seals neither bank. Shear results remain available separately.
    pub fn finish(
        &mut self,
        cx: &Cx<'_>,
    ) -> Result<Vec<CriticalPlaneMinerResult>, PlaneFatigueError> {
        poll(cx)?;
        if self.remaining_increments() != 0 {
            return Err(PlaneFatigueError::PathIncomplete);
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
        poll(cx)?;
        // This includes the final cancellation poll and atomic shear rollback.
        // Publishing the staged normal bank and finished flag cannot then fail.
        shear::finish_cycles(&mut self.fatemi_socie, cx)?;
        self.gauges = gauges;
        self.finished = true;
        Ok(results)
    }
}

fn poll(cx: &Cx<'_>) -> Result<(), PlaneFatigueError> {
    cx.checkpoint().map_err(|_| PlaneFatigueError::Cancelled)
}
fn gauge_error(gauge: usize, error: CriticalPlaneMinerError) -> PlaneFatigueError {
    match error {
        CriticalPlaneMinerError::Cancelled => PlaneFatigueError::Cancelled,
        error => PlaneFatigueError::Gauge { gauge, error },
    }
}
fn endpoint<'a, L: SourceBoundFatigueLaw>(
    problem: &'a SmallStrainPlaneProblem<'_, L>,
    load_path: &'a [f64],
) -> SmallStrainPlaneProblem<'a, L> {
    SmallStrainPlaneProblem {
        mesh: problem.mesh,
        laws: match &problem.laws {
            SmallStrainLawField::Uniform(law) => SmallStrainLawField::Uniform(law),
            SmallStrainLawField::PerElement(laws) => SmallStrainLawField::PerElement(laws),
            SmallStrainLawField::PerQuadrature(laws) => SmallStrainLawField::PerQuadrature(laws),
        },
        plane: problem.plane,
        formulation: problem.formulation,
        thickness_m: problem.thickness_m,
        reduction: problem.reduction,
        rotation_limit_rad: problem.rotation_limit_rad,
        prescribed_m: problem.prescribed_m,
        nodal_forces_n: problem.nodal_forces_n,
        load_path,
        budget: problem.budget,
        settings: problem.settings,
    }
}
fn observe<L: SourceBoundFatigueLaw>(
    problem: &SmallStrainPlaneProblem<'_, L>,
    solution: &SmallStrainPlaneSolution<L::State>,
    gauges: &mut [GaugeState],
    cx: &Cx<'_>,
) -> Result<Vec<TensorFatigueSample>, PlaneFatigueError> {
    let mut samples = Vec::with_capacity(gauges.len());
    for (index, gauge) in gauges.iter_mut().enumerate() {
        poll(cx)?;
        let sample = sample_at(problem, solution, gauge.element, gauge.point)?;
        gauge
            .counter
            .push_with_cancel(sample, &mut || cx.checkpoint().is_err())
            .map_err(|error| gauge_error(index, error))?;
        samples.push(sample);
    }
    Ok(samples)
}
fn sample_at<L: SourceBoundFatigueLaw>(
    problem: &SmallStrainPlaneProblem<'_, L>,
    solution: &SmallStrainPlaneSolution<L::State>,
    element: usize,
    quadrature_point: usize,
) -> Result<TensorFatigueSample, PlaneFatigueError> {
    let point = solution.element_offsets[element] + quadrature_point;
    let state = &solution.quadrature_states[point];
    problem
        .law(element, quadrature_point)
        .validate_fatigue_history(&state.strain, &state.history)
        .map_err(|_| {
            PlaneFatigueError::Invalid("native fatigue material history is inadmissible")
        })?;
    Ok(TensorFatigueSample {
        stress_pa: solution.quadrature_stress_pa[point],
        total_mechanical_strain: state.strain,
    })
}

fn absolute_temperature(material: &ResolvedMaterialStatePoint, axis: &str) -> bool {
    material
        .axis_quantities()
        .get(axis)
        .is_some_and(|quantity| {
            *quantity == QuantitySpec::dimensional(Temperature::DIMS)
                || quantity.semantic_type().is_some_and(|semantic| {
                    semantic.kind() == QuantityKind::AbsoluteTemperature
                        && semantic.form() == ValueForm::Static
                })
        })
}

fn bind_curve<L: SourceBoundFatigueLaw>(
    law: &L,
    index: usize,
    gauge: &PlaneFatigueGauge,
) -> Result<(), PlaneFatigueError> {
    let material = law.fatigue_source();
    let scope = gauge.curve.scope();
    let mismatch = |what| PlaneFatigueError::SourceMismatch { gauge: index, what };
    if gauge.material_state_identity != material.identity()
        || scope.material != material.material().chemistry
        || scope.process != material.material().process
    {
        return Err(mismatch(
            "CMB identity/material/process differs from the native plane law",
        ));
    }
    let axis = gauge
        .curve
        .source()
        .map_or("T", |s| s.request().temperature_axis.as_str());
    if !absolute_temperature(material, axis)
        || material.state_coordinate(axis) != Some(scope.temperature_k)
        || gauge.curve.parameters().youngs_modulus_pa != law.fatigue_backbone().elastic.youngs
    {
        return Err(mismatch(
            "CMB absolute temperature/modulus differs from the frozen plane law",
        ));
    }
    if let Some(source) = gauge.curve.source() {
        let youngs = material
            .property(YOUNG_MODULUS_PROPERTY)
            .ok_or_else(|| mismatch("original modulus receipt is missing or ambiguous"))?;
        if !source.matches_material_state(
            material.card_identity(),
            material.material(),
            &material.query(),
            youngs.answer(),
        ) {
            return Err(mismatch(
                "CMB card/state/conditions/modulus claim differs from the native source",
            ));
        }
    }
    Ok(())
}
