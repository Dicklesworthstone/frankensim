//! Multiaxial fatigue-limit usage on closed native elastic plate stress blocks.
//!
//! The shared Findley owner retains the full tensor path and evaluates its
//! two-dimensional shear-traction envelope on each declared fixed plane.
//! Its usage factor is separate from normal-stress Miner damage. No cycle
//! count, repetition count, finite life or cyclic stabilization is inferred.

use std::sync::Arc;

use fs_material::fatigue::critical_plane::findley::{
    FindleyCalibration, FindleyCycle, FindleyCycleResult, FindleyError, FindleyLimits,
    StressClosureTolerance,
};
use fs_qty::semantic::{QuantityKind, ValueForm};
use fs_qty::{QuantitySpec, Temperature};

use super::{PlateFatigueError, PlateFatigueSection, dot, local, poll};

/// One physical material point with a caller-declared finite plane search.
#[derive(Debug, Clone)]
pub struct PlateFindleyGauge {
    pub element: usize,
    /// Original DKT midside quadrature point, in `0..3`.
    pub surface_point: usize,
    /// Signed distance from the centered midsurface, in `[-h/2, h/2]` [m].
    pub depth_m: f64,
    pub temperature_axis: String,
    pub calibration: Arc<FindleyCalibration>,
    /// Three-dimensional normals in the original fixed mesh x/y/z frame.
    pub normals_in_mesh_frame: Vec<[f64; 3]>,
    pub limits: FindleyLimits,
    pub closure: StressClosureTolerance,
}

/// Aggregate upper bounds for all Findley gauges, in addition to native mesh
/// and gauge admission. One retained tensor sample contains six `f64` values.
/// Projection/circle limits sum the corresponding per-gauge declared ceilings.
/// The default admits only an empty Findley bank.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct PlateFindleyBudget {
    pub maximum_total_samples: u64,
    pub maximum_total_planes: usize,
    pub maximum_total_projection_work: u64,
    pub maximum_total_circle_work: u64,
}

pub(super) fn admit_budget(
    gauges: &[PlateFindleyGauge],
    budget: PlateFindleyBudget,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), PlateFatigueError> {
    let (mut samples, mut planes, mut projection, mut circle) = (0_u64, 0_usize, 0_u64, 0_u64);
    for gauge in gauges {
        poll(cancelled)?;
        samples = samples
            .checked_add(gauge.limits.maximum_samples)
            .ok_or(PlateFatigueError::Budget)?;
        planes = planes
            .checked_add(gauge.normals_in_mesh_frame.len())
            .ok_or(PlateFatigueError::Budget)?;
        projection = projection
            .checked_add(gauge.limits.maximum_projection_work)
            .ok_or(PlateFatigueError::Budget)?;
        circle = circle
            .checked_add(gauge.limits.maximum_circle_work)
            .ok_or(PlateFatigueError::Budget)?;
        if samples > budget.maximum_total_samples
            || planes > budget.maximum_total_planes
            || projection > budget.maximum_total_projection_work
            || circle > budget.maximum_total_circle_work
        {
            return Err(PlateFatigueError::Budget);
        }
    }
    Ok(())
}

pub(super) fn error(gauge: usize, error: FindleyError) -> PlateFatigueError {
    if error == FindleyError::Cancelled {
        PlateFatigueError::Cancelled
    } else {
        PlateFatigueError::FindleyGauge { gauge, error }
    }
}

/// The original compiled plane search is shared by fresh observation windows.
#[derive(Debug)]
pub(super) struct PreparedGauge {
    pub(super) binding: PlateFindleyGauge,
    dofs: [Option<usize>; 9],
    /// [sigma_xx, sigma_yy, tau_xy] rows, with engineering twist converted once.
    stress: [[f64; 9]; 3],
    membrane_stress_pa: f64,
    seed: FindleyCycle,
}

pub(super) fn admit_gauge(
    index: usize,
    gauge: &PlateFindleyGauge,
    sections: &[PlateFatigueSection],
) -> Result<FindleyCycle, PlateFatigueError> {
    let invalid = PlateFatigueError::Invalid;
    let section = sections
        .get(gauge.element)
        .ok_or_else(|| invalid("Findley element is outside the original plate"))?;
    if gauge.surface_point >= 3
        || !gauge.depth_m.is_finite()
        || gauge.depth_m.abs() > 0.5 * section.thickness_m
    {
        return Err(invalid(
            "Findley point must lie in its original DKT section",
        ));
    }
    let material = section.material.resolved();
    let scope = gauge.calibration.scope();
    let mismatch = |what| PlateFatigueError::FindleySourceMismatch { gauge: index, what };
    let absolute = material
        .axis_quantities()
        .get(&gauge.temperature_axis)
        .is_some_and(|quantity| {
            *quantity == QuantitySpec::dimensional(Temperature::DIMS)
                || quantity.semantic_type().is_some_and(|semantic| {
                    semantic.kind() == QuantityKind::AbsoluteTemperature
                        && semantic.form() == ValueForm::Static
                })
        });
    if !absolute
        || material.state_coordinate(&gauge.temperature_axis) != Some(scope.temperature_k)
        || scope.material != material.material().chemistry
        || scope.process != material.material().process
    {
        return Err(mismatch(
            "Findley material/process/absolute temperature differs from the section",
        ));
    }
    if gauge.calibration.source().is_some_and(|source| {
        source.request().temperature_axis != gauge.temperature_axis
            || !source.matches_material_state(
                material.card_identity(),
                material.material(),
                &material.query(),
            )
    }) {
        return Err(mismatch(
            "Findley original card/state/fixed source context differs from the section",
        ));
    }
    FindleyCycle::new(
        gauge.calibration.clone(),
        gauge.normals_in_mesh_frame.clone(),
        gauge.limits,
        gauge.closure,
    )
    .map_err(|e| error(index, e))
}

impl PreparedGauge {
    pub(super) fn new(
        binding: PlateFindleyGauge,
        seed: FindleyCycle,
        section: &PlateFatigueSection,
        dofs: [Option<usize>; 9],
        curvature: &[[f64; 9]; 3],
        pretension: f64,
    ) -> Result<Self, PlateFatigueError> {
        let e = section.material.young_modulus_pa();
        let nu = section.material.poisson_ratio();
        let scale = e / (1.0 - nu * nu);
        let z = binding.depth_m;
        let stress = [
            std::array::from_fn(|j| z * scale * (curvature[0][j] + nu * curvature[1][j])),
            std::array::from_fn(|j| z * scale * (nu * curvature[0][j] + curvature[1][j])),
            std::array::from_fn(|j| z * (e / (2.0 * (1.0 + nu))) * curvature[2][j]),
        ];
        let membrane_stress_pa = pretension / section.thickness_m;
        if stress.iter().flatten().any(|s| !s.is_finite()) || !membrane_stress_pa.is_finite() {
            return Err(PlateFatigueError::Invalid(
                "unrepresentable DKT Findley stress recovery",
            ));
        }
        Ok(Self {
            binding,
            dofs,
            stress,
            membrane_stress_pa,
            seed,
        })
    }

    /// Called only after the parent's complete native displacement/domain check.
    pub(super) fn recover(&self, q: &[f64]) -> Result<[f64; 6], PlateFatigueError> {
        let u = local(&self.dofs, q);
        let value = [
            self.membrane_stress_pa + dot(&self.stress[0], &u),
            self.membrane_stress_pa + dot(&self.stress[1], &u),
            0.0,
            dot(&self.stress[2], &u),
            0.0,
            0.0,
        ];
        if value.iter().any(|s| !s.is_finite()) {
            return Err(PlateFatigueError::Invalid("plate Findley stress overflow"));
        }
        Ok(value)
    }
}

/// All gauges share the same observation window and accepted sample clock.
#[derive(Debug, Clone)]
pub(super) struct Bank {
    counters: Vec<FindleyCycle>,
}

impl Bank {
    pub(super) fn new(gauges: &[PreparedGauge]) -> Result<Self, PlateFatigueError> {
        Ok(Self {
            counters: gauges
                .iter()
                .enumerate()
                .map(|(i, g)| g.seed.fresh().map_err(|e| error(i, e)))
                .collect::<Result<_, _>>()?,
        })
    }
    pub(super) fn sample_count(&self, gauge: usize) -> Option<u64> {
        self.counters.get(gauge).map(FindleyCycle::sample_count)
    }
    pub(super) fn results(&self) -> Vec<Option<FindleyCycleResult>> {
        self.counters.iter().map(|c| c.result().copied()).collect()
    }
    pub(super) fn push(
        &mut self,
        samples: &[[f64; 6]],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), PlateFatigueError> {
        poll(cancelled)?;
        if samples.len() != self.counters.len() {
            return Err(PlateFatigueError::Invalid(
                "incomplete Findley gauge sample",
            ));
        }
        // Never invalidate a previously finalized window as a rollback side effect.
        if let Some((index, _)) = self
            .counters
            .iter()
            .enumerate()
            .find(|(_, c)| c.result().is_some())
        {
            return Err(error(index, FindleyError::Finished));
        }
        let previous = self.counters.first().map_or(0, FindleyCycle::sample_count);
        let attempt = (|| {
            for (i, (counter, sample)) in self.counters.iter_mut().zip(samples).enumerate() {
                counter
                    .push_with_cancel(*sample, cancelled)
                    .map_err(|e| error(i, e))?;
            }
            poll(cancelled)
        })();
        if let Err(failure) = attempt {
            for (i, counter) in self.counters.iter_mut().enumerate() {
                counter.discard_after(previous).map_err(|e| error(i, e))?;
            }
            return Err(failure);
        }
        Ok(())
    }
    pub(super) fn preview(
        &self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<FindleyCycleResult>, PlateFatigueError> {
        poll(cancelled)?;
        let mut results = Vec::with_capacity(self.counters.len());
        for (i, counter) in self.counters.iter().enumerate() {
            poll(cancelled)?;
            results.push(match counter.result() {
                Some(result) => *result,
                None => counter
                    .clone()
                    .finish_with_cancel(cancelled)
                    .map_err(|e| error(i, e))?,
            });
        }
        poll(cancelled)?;
        Ok(results)
    }
    pub(super) fn finish(
        &mut self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<FindleyCycleResult>, PlateFatigueError> {
        poll(cancelled)?;
        let already_finished: Vec<_> = self.counters.iter().map(|c| c.result().is_some()).collect();
        let attempt = (|| {
            let mut results = Vec::with_capacity(self.counters.len());
            for (i, counter) in self.counters.iter_mut().enumerate() {
                results.push(
                    counter
                        .finish_with_cancel(cancelled)
                        .map_err(|e| error(i, e))?,
                );
            }
            poll(cancelled)?;
            Ok(results)
        })();
        if attempt.is_err() {
            for (i, (counter, was_finished)) in
                self.counters.iter_mut().zip(already_finished).enumerate()
            {
                if !was_finished {
                    counter
                        .discard_after(counter.sample_count())
                        .map_err(|e| error(i, e))?;
                }
            }
        }
        attempt
    }
    pub(super) fn restart(
        &mut self,
        samples: &[[f64; 6]],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), PlateFatigueError> {
        poll(cancelled)?;
        let mut next = Self {
            counters: self
                .counters
                .iter()
                .enumerate()
                .map(|(i, c)| c.fresh().map_err(|e| error(i, e)))
                .collect::<Result<_, _>>()?,
        };
        next.push(samples, cancelled)?;
        *self = next;
        Ok(())
    }
}
