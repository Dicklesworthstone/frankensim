//! Representative shear windows on the same accepted native shell points.
//! Prefix rollback keeps accumulated tensor storage out of ordinary step copies.

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
    POISSON_RATIO_PROPERTY, ResolvedMaterialStatePoint, YIELD_STRESS_PROPERTY,
    YOUNG_MODULUS_PROPERTY,
};
use fs_qty::semantic::{QuantityKind, ValueForm};
use fs_qty::{QuantitySpec, Temperature};

use super::{
    ShellFatigueError, SmallStrainShellProblem, SmallStrainShellSolution, accepted_sample,
    point_law, poll,
};

/// One original material point and a finite-plane representative shear cycle.
/// The torsional calibration is independent of an axial strain-life fit. Plane
/// directions use the same fixed local facet frame as the native point tensors.
#[derive(Debug, Clone)]
pub struct ShellFatemiSocieGauge {
    pub element: usize,
    /// Original DKT midside quadrature index, in `0..3`.
    pub surface_point: usize,
    /// Index in the assigned section's original thickness `points()` array.
    pub thickness_point: usize,
    pub material_state_identity: ContentHash,
    pub calibration: Arc<FatemiSocieCalibration>,
    pub normals_in_facet_frame: Vec<[f64; 3]>,
    pub limits: FatemiSocieLimits,
    pub closure: CycleClosureTolerance,
}

/// Aggregate bounds on the sum of declared representative-window limits.
/// Samples are retained once per gauge, not once per plane. Each explicit
/// restart opens a new bounded window; none of these budgets implies repeated
/// cycle exposure or a lifetime damage estimate. Zero bounds admit no shear bank.
#[derive(Debug, Clone, Copy, Default)]
pub struct ShellFatemiSocieConfig {
    pub maximum_total_samples: usize,
    pub maximum_total_projection_work: u64,
    pub maximum_total_circle_work: u64,
}

#[derive(Debug, Clone)]
struct ShearGaugeState {
    element: usize,
    surface_point: usize,
    thickness_point: usize,
    counter: FatemiSocieCycle,
}

#[derive(Debug, Clone)]
pub(super) struct ShearBank {
    gauges: Vec<ShearGaugeState>,
    samples: Vec<TensorFatigueSample>,
}

impl ShearBank {
    pub(super) fn new<L: SourceBoundFatigueLaw>(
        problem: &SmallStrainShellProblem<'_, L>,
        gauges: Vec<ShellFatemiSocieGauge>,
        source_offset: usize,
        config: ShellFatemiSocieConfig,
        cx: &Cx<'_>,
    ) -> Result<Self, ShellFatigueError> {
        let (mut samples, mut projection, mut circle) = (0usize, 0u64, 0u64);
        for gauge in &gauges {
            poll(cx)?;
            samples = usize::try_from(gauge.limits.maximum_samples)
                .ok()
                .and_then(|count| samples.checked_add(count))
                .ok_or(ShellFatigueError::Budget)?;
            projection = projection
                .checked_add(gauge.limits.maximum_projection_work)
                .ok_or(ShellFatigueError::Budget)?;
            circle = circle
                .checked_add(gauge.limits.maximum_circle_work)
                .ok_or(ShellFatigueError::Budget)?;
            if samples > config.maximum_total_samples
                || projection > config.maximum_total_projection_work
                || circle > config.maximum_total_circle_work
            {
                return Err(ShellFatigueError::Budget);
            }
        }
        let mut prepared = Vec::with_capacity(gauges.len());
        for (index, gauge) in gauges.into_iter().enumerate() {
            poll(cx)?;
            bind_calibration(problem, source_offset + index, &gauge)?;
            let counter = FatemiSocieCycle::new(
                gauge.calibration,
                gauge.normals_in_facet_frame,
                gauge.limits,
                gauge.closure,
            )
            .map_err(|error| shear_error(index, error))?;
            prepared.push(ShearGaugeState {
                element: gauge.element,
                surface_point: gauge.surface_point,
                thickness_point: gauge.thickness_point,
                counter,
            });
        }
        poll(cx)?;
        Ok(Self {
            gauges: prepared,
            samples: Vec::new(),
        })
    }

    pub(super) fn len(&self) -> usize {
        self.gauges.len()
    }

    pub(super) fn results(&self) -> Vec<Option<FatemiSocieCycleResult>> {
        self.gauges
            .iter()
            .map(|g| g.counter.result().copied())
            .collect()
    }

    pub(super) fn samples(&self) -> &[TensorFatigueSample] {
        &self.samples
    }

    pub(super) fn calibration(&self, gauge: usize) -> Option<&FatemiSocieCalibration> {
        self.gauges.get(gauge).map(|g| g.counter.calibration())
    }

    pub(super) fn material_source<'a, L: SourceBoundFatigueLaw>(
        &self,
        problem: &'a SmallStrainShellProblem<'_, L>,
        gauge: usize,
    ) -> Option<&'a ResolvedMaterialStatePoint> {
        let gauge = self.gauges.get(gauge)?;
        Some(point_law(problem, gauge.element, gauge.thickness_point).fatigue_source())
    }

    pub(super) fn admit_advance(&self) -> Result<(), ShellFatigueError> {
        if self.gauges.iter().any(|g| g.counter.result().is_some()) {
            return Err(ShellFatigueError::Invalid(
                "restart closed Fatemi-Socie windows before advancing native motion",
            ));
        }
        Ok(())
    }

    fn checkpoint(&self) -> Vec<(u64, bool)> {
        self.gauges
            .iter()
            .map(|g| (g.counter.sample_count(), g.counter.result().is_some()))
            .collect()
    }

    fn restore(&mut self, checkpoint: &[(u64, bool)]) {
        for (gauge, &(prefix, was_closed)) in self.gauges.iter_mut().zip(checkpoint) {
            if !was_closed {
                gauge
                    .counter
                    .discard_after(prefix)
                    .expect("shell shear rollback retains a previously accepted prefix");
            }
        }
    }

    /// Last fallible stage of a native endpoint proposal. On refusal, every
    /// appended sample/candidate result is removed; latest accepted samples
    /// remain untouched. No previously accumulated tensor history is copied.
    pub(super) fn observe<L: SourceBoundFatigueLaw>(
        &mut self,
        problem: &SmallStrainShellProblem<'_, L>,
        solution: &SmallStrainShellSolution<L::State>,
        close_cycle: bool,
        cx: &Cx<'_>,
    ) -> Result<(), ShellFatigueError> {
        self.admit_advance()?;
        let checkpoint = self.checkpoint();
        let proposal = (|| {
            let mut samples = Vec::with_capacity(self.gauges.len());
            for (index, gauge) in self.gauges.iter_mut().enumerate() {
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
                    .map_err(|error| shear_error(index, error))?;
                samples.push(sample);
            }
            if close_cycle {
                self.close_inner(cx)?;
            }
            poll(cx)?;
            Ok(samples)
        })();
        match proposal {
            Ok(samples) => {
                self.samples = samples;
                Ok(())
            }
            Err(error) => {
                self.restore(&checkpoint);
                Err(error)
            }
        }
    }

    pub(super) fn restart<L: SourceBoundFatigueLaw>(
        &mut self,
        problem: &SmallStrainShellProblem<'_, L>,
        solution: &SmallStrainShellSolution<L::State>,
        cx: &Cx<'_>,
    ) -> Result<(), ShellFatigueError> {
        let mut gauges = Vec::with_capacity(self.gauges.len());
        let mut samples = Vec::with_capacity(self.gauges.len());
        for (index, gauge) in self.gauges.iter().enumerate() {
            poll(cx)?;
            let sample = accepted_sample(
                problem,
                solution,
                gauge.element,
                gauge.surface_point,
                gauge.thickness_point,
            )?;
            let mut counter = gauge
                .counter
                .fresh()
                .map_err(|error| shear_error(index, error))?;
            counter
                .push_with_cancel(sample, &mut || cx.checkpoint().is_err())
                .map_err(|error| shear_error(index, error))?;
            gauges.push(ShearGaugeState {
                element: gauge.element,
                surface_point: gauge.surface_point,
                thickness_point: gauge.thickness_point,
                counter,
            });
            samples.push(sample);
        }
        poll(cx)?;
        self.gauges = gauges;
        self.samples = samples;
        Ok(())
    }

    pub(super) fn close(
        &mut self,
        cx: &Cx<'_>,
    ) -> Result<Vec<FatemiSocieCycleResult>, ShellFatigueError> {
        let checkpoint = self.checkpoint();
        let result = self.close_inner(cx);
        if result.is_err() {
            self.restore(&checkpoint);
        }
        result
    }

    fn close_inner(
        &mut self,
        cx: &Cx<'_>,
    ) -> Result<Vec<FatemiSocieCycleResult>, ShellFatigueError> {
        let mut results = Vec::with_capacity(self.gauges.len());
        for (index, gauge) in self.gauges.iter_mut().enumerate() {
            results.push(
                gauge
                    .counter
                    .finish_with_cancel(&mut || cx.checkpoint().is_err())
                    .map_err(|error| shear_error(index, error))?,
            );
        }
        poll(cx)?;
        Ok(results)
    }
}

fn shear_error(gauge: usize, error: FatemiSocieError) -> ShellFatigueError {
    match error {
        FatemiSocieError::Cancelled => ShellFatigueError::Cancelled,
        error => ShellFatigueError::FatemiSocie { gauge, error },
    }
}

fn bind_calibration<L: SourceBoundFatigueLaw>(
    problem: &SmallStrainShellProblem<'_, L>,
    index: usize,
    gauge: &ShellFatemiSocieGauge,
) -> Result<(), ShellFatigueError> {
    if gauge.element >= problem.mesh.tris.len() || gauge.surface_point >= 3 {
        return Err(ShellFatigueError::Invalid(
            "shear gauge is outside the original shell surface quadrature",
        ));
    }
    let section = problem.sections.at(gauge.element);
    let point = section
        .points()
        .get(gauge.thickness_point)
        .ok_or(ShellFatigueError::Invalid(
            "shear gauge is outside the original shell thickness quadrature",
        ))?;
    let law = &section.layers()[point.layer].law;
    let material = law.fatigue_source();
    let scope = gauge.calibration.scope();
    let parameters = gauge.calibration.parameters();
    let mismatch = |what| ShellFatigueError::SourceMismatch { gauge: index, what };
    if gauge.material_state_identity != material.identity()
        || scope.material != material.material().chemistry
        || scope.process != material.material().process
    {
        return Err(mismatch(
            "shear-life identity/material/process differs from the native shell layer",
        ));
    }
    let axis = gauge
        .calibration
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
        || parameters.shear_modulus_pa != law.fatigue_backbone().elastic.lame().1
        || parameters.initial_yield_stress_pa != law.fatigue_backbone().yield_stress
    {
        return Err(mismatch(
            "shear-life absolute temperature, shear modulus or initial-yield reference differs from the frozen shell law",
        ));
    }
    if let Some(source) = gauge.calibration.source() {
        let answer = |name| {
            material
                .property(name)
                .map(|property| property.answer())
                .ok_or_else(|| {
                    mismatch("original shell elastic/initial-yield receipt is missing or ambiguous")
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
                "shear-life card, state, typed conditions or elastic/initial-yield source differs from the shell layer",
            ));
        }
    }
    Ok(())
}
