//! Shear-controlled life from the same accepted native plane-continuum tensors.
//!
//! The material owner evaluates one explicitly declared closed representative
//! cycle with its independent torsional strain-life calibration. Engineering
//! shear amplitude comes from the full two-component shear path on each plane;
//! the original maximum-shear-plane selection and full-cycle normal-stress
//! maximum remain in fs-material. No axial coefficient is relabeled, no cycle
//! repetitions are inferred. Explicitly enabled FS/Miner accumulates one cycle
//! per closed window on each fixed plane through that same material owner.

use std::sync::Arc;

use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_material::fatigue::critical_plane::fatemi_socie::{
    FatemiSocieCalibration, FatemiSocieCycle, FatemiSocieCycleResult, FatemiSocieError,
    FatemiSocieLimits,
};
use fs_material::fatigue::critical_plane::{CycleClosureTolerance, TensorFatigueSample};
use fs_material::fatigue::cyclic::SourceBoundFatigueLaw;
use fs_material::state_point::{
    POISSON_RATIO_PROPERTY, YIELD_STRESS_PROPERTY, YOUNG_MODULUS_PROPERTY,
};

use super::{
    Element, PlaneFatigueError, SmallStrainPlaneProblem, SmallStrainPlaneSolution,
    absolute_temperature, poll, sample_at,
};

/// One original P1/Q1 integration point and its separate shear-life fit.
#[derive(Debug, Clone)]
pub struct PlaneFatemiSocieGauge {
    pub element: usize,
    pub quadrature_point: usize,
    pub material_state_identity: ContentHash,
    pub calibration: Arc<FatemiSocieCalibration>,
    /// Fixed mesh x/y axes, with z the out-of-plane direction.
    pub normals_in_mesh_frame: Vec<[f64; 3]>,
    pub limits: FatemiSocieLimits,
    pub closure: CycleClosureTolerance,
}

/// Aggregate bounds across all retained representative-cycle windows.
/// Sample bounds count complete paired stress/strain tensors, not components.
/// Projection and circle bounds sum each gauge's declared per-window allowance.
/// Default admits no shear gauges; normal-strain Miner budgets stay separate.
#[derive(Debug, Clone, Copy, Default)]
pub struct PlaneShearFatigueBudget {
    pub maximum_retained_samples: u64,
    pub maximum_planes: usize,
    pub maximum_projection_work: u64,
    pub maximum_circle_work: u64,
}

#[derive(Debug, Clone)]
pub(super) struct GaugeState {
    pub(super) element: usize,
    pub(super) point: usize,
    pub(super) counter: FatemiSocieCycle,
}

pub(super) fn admit_budget(
    gauges: &[PlaneFatemiSocieGauge],
    budget: PlaneShearFatigueBudget,
    cx: &Cx<'_>,
) -> Result<(), PlaneFatigueError> {
    let (mut samples, mut planes, mut projection, mut circle) = (0u64, 0usize, 0u64, 0u64);
    for gauge in gauges {
        poll(cx)?;
        samples = samples
            .checked_add(gauge.limits.maximum_samples)
            .ok_or(PlaneFatigueError::Budget)?;
        planes = planes
            .checked_add(gauge.normals_in_mesh_frame.len())
            .ok_or(PlaneFatigueError::Budget)?;
        projection = projection
            .checked_add(gauge.limits.maximum_projection_work)
            .ok_or(PlaneFatigueError::Budget)?;
        circle = circle
            .checked_add(gauge.limits.maximum_circle_work)
            .ok_or(PlaneFatigueError::Budget)?;
        if samples > budget.maximum_retained_samples
            || planes > budget.maximum_planes
            || projection > budget.maximum_projection_work
            || circle > budget.maximum_circle_work
        {
            return Err(PlaneFatigueError::Budget);
        }
    }
    Ok(())
}

pub(super) fn prepare<L: SourceBoundFatigueLaw>(
    problem: &SmallStrainPlaneProblem<'_, L>,
    elements: &[Element],
    gauges: Vec<PlaneFatemiSocieGauge>,
    source_index_offset: usize,
    cx: &Cx<'_>,
) -> Result<Vec<GaugeState>, PlaneFatigueError> {
    let mut bank = Vec::with_capacity(gauges.len());
    for (index, gauge) in gauges.into_iter().enumerate() {
        poll(cx)?;
        if elements
            .get(gauge.element)
            .is_none_or(|element| gauge.quadrature_point >= element.points.len())
        {
            return Err(PlaneFatigueError::Invalid(
                "shear gauge is outside the original plane quadrature",
            ));
        }
        bind(
            problem.law(gauge.element, gauge.quadrature_point),
            source_index_offset + index,
            &gauge,
        )?;
        bank.push(GaugeState {
            element: gauge.element,
            point: gauge.quadrature_point,
            counter: FatemiSocieCycle::new(
                gauge.calibration,
                gauge.normals_in_mesh_frame,
                gauge.limits,
                gauge.closure,
            )
            .map_err(|error| gauge_error(index, error))?,
        });
    }
    Ok(bank)
}

pub(super) fn observe<L: SourceBoundFatigueLaw>(
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

pub(super) fn checkpoint(gauges: &[GaugeState]) -> Vec<(u64, bool)> {
    gauges
        .iter()
        .map(|g| (g.counter.sample_count(), g.counter.result().is_none()))
        .collect()
}

/// Previously accepted prefixes are valid by construction. Restoration neither
/// allocates nor polls cancellation, and preserves already finalized results.
pub(super) fn restore(gauges: &mut [GaugeState], checkpoint: &[(u64, bool)]) {
    for (gauge, &(prefix, was_open)) in gauges.iter_mut().zip(checkpoint) {
        if was_open {
            gauge
                .counter
                .discard_after(prefix)
                .expect("native plane rollback retains a previously accepted shear prefix");
        }
    }
}

pub(super) fn restart<L: SourceBoundFatigueLaw>(
    problem: &SmallStrainPlaneProblem<'_, L>,
    solution: &SmallStrainPlaneSolution<L::State>,
    gauges: &[GaugeState],
    cx: &Cx<'_>,
) -> Result<(Vec<GaugeState>, Vec<TensorFatigueSample>), PlaneFatigueError> {
    let mut next = Vec::with_capacity(gauges.len());
    for (index, gauge) in gauges.iter().enumerate() {
        poll(cx)?;
        next.push(GaugeState {
            element: gauge.element,
            point: gauge.point,
            counter: gauge.counter.fresh().map_err(|e| gauge_error(index, e))?,
        });
    }
    let samples = observe(problem, solution, &mut next, cx)?;
    poll(cx)?;
    Ok((next, samples))
}

pub(super) fn finish_cycles(
    gauges: &mut [GaugeState],
    cx: &Cx<'_>,
) -> Result<Vec<FatemiSocieCycleResult>, PlaneFatigueError> {
    poll(cx)?;
    let previous = checkpoint(gauges);
    let result = (|| {
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
        Ok(results)
    })();
    if result.is_err() {
        restore(gauges, &previous);
    }
    result
}

fn gauge_error(gauge: usize, error: FatemiSocieError) -> PlaneFatigueError {
    match error {
        FatemiSocieError::Cancelled => PlaneFatigueError::Cancelled,
        error => PlaneFatigueError::FatemiSocie { gauge, error },
    }
}

fn bind<L: SourceBoundFatigueLaw>(
    law: &L,
    index: usize,
    gauge: &PlaneFatemiSocieGauge,
) -> Result<(), PlaneFatigueError> {
    let material = law.fatigue_source();
    let calibration = &gauge.calibration;
    let scope = calibration.scope();
    let parameters = calibration.parameters();
    let mismatch = |what| PlaneFatigueError::SourceMismatch { gauge: index, what };
    if gauge.material_state_identity != material.identity()
        || scope.material != material.material().chemistry
        || scope.process != material.material().process
    {
        return Err(mismatch(
            "shear-life identity/material/process differs from the native plane law",
        ));
    }
    let axis = calibration
        .source()
        .map_or("T", |source| source.request().temperature_axis.as_str());
    if !absolute_temperature(material, axis)
        || material.state_coordinate(axis) != Some(scope.temperature_k)
        || parameters.shear_modulus_pa != law.fatigue_backbone().elastic.lame().1
        || parameters.initial_yield_stress_pa != law.fatigue_backbone().yield_stress
    {
        return Err(mismatch(
            "shear-life absolute temperature, G or initial-yield reference differs from the native plane law",
        ));
    }
    if let Some(source) = calibration.source() {
        let answer = |name| {
            material.property(name).map(|p| p.answer()).ok_or_else(|| {
                mismatch("original elastic/initial-yield receipt is missing or ambiguous")
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
            return Err(mismatch(
                "shear-life card/state/conditions or original elastic/initial-yield source differs from the plane law",
            ));
        }
    }
    Ok(())
}
