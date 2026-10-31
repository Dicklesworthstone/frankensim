//! Variable-amplitude critical-plane strain-life accumulation.
//!
//! Each accepted tensor sample is projected onto every caller-declared fixed
//! plane, then counted by the existing bounded rainflow/Palmgren-Miner owner.
//! The selected plane is the one with greatest accumulated Miner damage over
//! the COMPLETE observation. This is deliberately different from the single-
//! representative-cycle SWT and Fatemi-Socie searches.
//!
//! Mean-stress correction is explicit. Uncorrected counting uses only the
//! signed normal mechanical strain history. Morrow/SWT use actual coincident
//! normal stress/strain observations and retain the paired-history refusal when
//! phase lag or relaxation would make scalar endpoint pairing ambiguous.
//!
//! This supplies a finite-plane approximation to multiaxial variable-amplitude
//! fatigue. It does not optimize over the continuum sphere, infer missed
//! extrema, model load-sequence interaction beyond linear Miner summation,
//! establish cyclic stabilization, or convert damage into stiffness/cracks.

use std::sync::Arc;

use super::{SearchPlane, TensorFatigueSample};
use crate::fatigue::strain_life::{
    PairedStrainLifeMiner, StrainLifeCurve, StrainLifeMeanStress, StrainLifeMiner,
    StrainLifeSummary, StrainStressSample,
};
use crate::fatigue::{FatigueError, RainflowLimits};

/// Bounded finite-plane search plus the existing per-plane rainflow storage.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CriticalPlaneMinerLimits {
    pub rainflow: RainflowLimits,
    pub maximum_planes: usize,
    /// Total accepted samples times declared planes.
    pub maximum_projection_work: u64,
}

/// Final accumulated damage on the most damaging declared plane.
#[derive(Debug, Clone, PartialEq)]
pub struct CriticalPlaneMinerResult {
    pub samples: u64,
    pub evaluated_planes: usize,
    pub selected_plane: usize,
    pub normal_in_frame: [f64; 3],
    pub summary: StrainLifeSummary,
    pub projection_work: u64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum CriticalPlaneMinerError {
    Invalid(&'static str),
    Budget,
    Fatigue { plane: usize, error: FatigueError },
    Cancelled,
    Finished,
}
impl core::fmt::Display for CriticalPlaneMinerError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "critical-plane variable-amplitude fatigue refused: {self:?}")
    }
}
impl std::error::Error for CriticalPlaneMinerError {}

#[derive(Debug, Clone, PartialEq)]
enum PlaneCounter {
    Scalar(StrainLifeMiner),
    Paired(PairedStrainLifeMiner),
}
impl PlaneCounter {
    fn new(
        curve: Arc<StrainLifeCurve>,
        correction: StrainLifeMeanStress,
        limits: RainflowLimits,
    ) -> Result<Self, FatigueError> {
        if correction == StrainLifeMeanStress::Uncorrected {
            Ok(Self::Scalar(StrainLifeMiner::new(curve, correction, limits)?))
        } else {
            Ok(Self::Paired(PairedStrainLifeMiner::new(
                curve, correction, limits,
            )?))
        }
    }

    fn push(
        &mut self,
        strain: f64,
        stress_pa: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), FatigueError> {
        match self {
            Self::Scalar(counter) => counter.push_with_cancel(strain, cancelled),
            Self::Paired(counter) => counter.push_with_cancel(
                StrainStressSample {
                    mechanical_strain: strain,
                    normal_stress_pa: stress_pa,
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
            Self::Scalar(counter) => counter.finish_with_cancel(cancelled),
            Self::Paired(counter) => counter.finish_with_cancel(cancelled),
        }
    }

    fn summary(&self) -> StrainLifeSummary {
        match self {
            Self::Scalar(counter) => counter.summary(),
            Self::Paired(counter) => counter.summary(),
        }
    }
}

/// Streaming finite-plane variable-amplitude Miner owner.
///
/// Pushes are transactional across the ENTIRE plane bank: if any later plane
/// refuses a sample, all earlier-plane candidates are discarded. Finish is
/// similarly atomic, so a life-domain or mean-stress refusal on one plane does
/// not publish half-cycle damage on another.
#[derive(Debug, Clone)]
pub struct CriticalPlaneMiner {
    curve: Arc<StrainLifeCurve>,
    correction: StrainLifeMeanStress,
    planes: Arc<Vec<SearchPlane>>,
    counters: Vec<PlaneCounter>,
    limits: CriticalPlaneMinerLimits,
    samples: u64,
    projection_work: u64,
    result: Option<CriticalPlaneMinerResult>,
}

impl CriticalPlaneMiner {
    pub fn new(
        curve: Arc<StrainLifeCurve>,
        correction: StrainLifeMeanStress,
        normals: Vec<[f64; 3]>,
        limits: CriticalPlaneMinerLimits,
    ) -> Result<Self, CriticalPlaneMinerError> {
        if normals.is_empty() {
            return Err(CriticalPlaneMinerError::Invalid(
                "declare at least one search plane",
            ));
        }
        if limits.maximum_planes == 0
            || limits.rainflow.maximum_samples == 0
            || limits.rainflow.maximum_residue == 0
            || limits.maximum_projection_work == 0
            || normals.len() > limits.maximum_planes
        {
            return Err(CriticalPlaneMinerError::Budget);
        }
        let plane_count =
            u64::try_from(normals.len()).map_err(|_| CriticalPlaneMinerError::Budget)?;
        if plane_count > limits.maximum_projection_work {
            return Err(CriticalPlaneMinerError::Budget);
        }

        let mut planes = Vec::new();
        planes
            .try_reserve_exact(normals.len())
            .map_err(|_| CriticalPlaneMinerError::Budget)?;
        let mut counters = Vec::new();
        counters
            .try_reserve_exact(normals.len())
            .map_err(|_| CriticalPlaneMinerError::Budget)?;

        for (index, normal) in normals.into_iter().enumerate() {
            let plane = SearchPlane::new(normal).map_err(|error| {
                CriticalPlaneMinerError::Fatigue {
                    plane: index,
                    error,
                }
            })?;
            let counter = PlaneCounter::new(curve.clone(), correction, limits.rainflow).map_err(
                |error| CriticalPlaneMinerError::Fatigue {
                    plane: index,
                    error,
                },
            )?;
            planes.push(plane);
            counters.push(counter);
        }

        Ok(Self {
            curve,
            correction,
            planes: Arc::new(planes),
            counters,
            limits,
            samples: 0,
            projection_work: 0,
            result: None,
        })
    }

    pub fn curve(&self) -> &StrainLifeCurve {
        &self.curve
    }

    pub fn mean_stress(&self) -> StrainLifeMeanStress {
        self.correction
    }

    pub fn sample_count(&self) -> u64 {
        self.samples
    }

    pub fn plane_count(&self) -> usize {
        self.planes.len()
    }

    pub fn projection_work(&self) -> u64 {
        self.projection_work
    }

    pub fn limits(&self) -> CriticalPlaneMinerLimits {
        self.limits
    }

    pub fn result(&self) -> Option<&CriticalPlaneMinerResult> {
        self.result.as_ref()
    }

    pub fn plane_summaries(&self) -> Vec<StrainLifeSummary> {
        self.counters.iter().map(PlaneCounter::summary).collect()
    }

    pub fn push(&mut self, sample: TensorFatigueSample) -> Result<(), CriticalPlaneMinerError> {
        self.push_with_cancel(sample, &mut || false)
    }

    pub fn push_with_cancel(
        &mut self,
        sample: TensorFatigueSample,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), CriticalPlaneMinerError> {
        if cancelled() {
            return Err(CriticalPlaneMinerError::Cancelled);
        }
        if self.result.is_some() {
            return Err(CriticalPlaneMinerError::Finished);
        }
        let next_samples = self
            .samples
            .checked_add(1)
            .ok_or(CriticalPlaneMinerError::Budget)?;
        let plane_count =
            u64::try_from(self.planes.len()).map_err(|_| CriticalPlaneMinerError::Budget)?;
        let next_work = self
            .projection_work
            .checked_add(plane_count)
            .ok_or(CriticalPlaneMinerError::Budget)?;
        if next_samples > self.limits.rainflow.maximum_samples
            || next_work > self.limits.maximum_projection_work
        {
            return Err(CriticalPlaneMinerError::Budget);
        }

        // Clone only bounded rainflow residues/counters; immutable planes/curve
        // remain shared. No full tensor history is retained.
        let mut candidate = self.counters.clone();
        for (index, (plane, counter)) in self
            .planes
            .iter()
            .copied()
            .zip(candidate.iter_mut())
            .enumerate()
        {
            if cancelled() {
                return Err(CriticalPlaneMinerError::Cancelled);
            }
            let (strain, stress_pa) = plane.project(sample).map_err(|error| {
                CriticalPlaneMinerError::Fatigue {
                    plane: index,
                    error,
                }
            })?;
            counter
                .push(strain, stress_pa, cancelled)
                .map_err(|error| CriticalPlaneMinerError::Fatigue {
                    plane: index,
                    error,
                })?;
        }
        if cancelled() {
            return Err(CriticalPlaneMinerError::Cancelled);
        }
        self.counters = candidate;
        self.samples = next_samples;
        self.projection_work = next_work;
        Ok(())
    }

    pub fn push_block(
        &mut self,
        samples: &[TensorFatigueSample],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), CriticalPlaneMinerError> {
        if cancelled() {
            return Err(CriticalPlaneMinerError::Cancelled);
        }
        let count = u64::try_from(samples.len()).map_err(|_| CriticalPlaneMinerError::Budget)?;
        let plane_count =
            u64::try_from(self.planes.len()).map_err(|_| CriticalPlaneMinerError::Budget)?;
        let next_samples = self
            .samples
            .checked_add(count)
            .ok_or(CriticalPlaneMinerError::Budget)?;
        let extra_work = count
            .checked_mul(plane_count)
            .ok_or(CriticalPlaneMinerError::Budget)?;
        let next_work = self
            .projection_work
            .checked_add(extra_work)
            .ok_or(CriticalPlaneMinerError::Budget)?;
        if next_samples > self.limits.rainflow.maximum_samples
            || next_work > self.limits.maximum_projection_work
        {
            return Err(CriticalPlaneMinerError::Budget);
        }
        let mut candidate = self.clone();
        for &sample in samples {
            candidate.push_with_cancel(sample, cancelled)?;
        }
        if cancelled() {
            return Err(CriticalPlaneMinerError::Cancelled);
        }
        *self = candidate;
        Ok(())
    }

    pub fn finish(&mut self) -> Result<CriticalPlaneMinerResult, CriticalPlaneMinerError> {
        self.finish_with_cancel(&mut || false)
    }

    pub fn finish_with_cancel(
        &mut self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<CriticalPlaneMinerResult, CriticalPlaneMinerError> {
        if cancelled() {
            return Err(CriticalPlaneMinerError::Cancelled);
        }
        if let Some(result) = &self.result {
            return Ok(result.clone());
        }
        if self.samples == 0 {
            return Err(CriticalPlaneMinerError::Invalid(
                "at least one tensor observation is required",
            ));
        }

        let mut candidate = self.counters.clone();
        let mut selected: Option<(usize, StrainLifeSummary)> = None;
        for (index, counter) in candidate.iter_mut().enumerate() {
            if cancelled() {
                return Err(CriticalPlaneMinerError::Cancelled);
            }
            let summary =
                counter
                    .finish(cancelled)
                    .map_err(|error| CriticalPlaneMinerError::Fatigue {
                        plane: index,
                        error,
                    })?;
            if selected.as_ref().is_none_or(|(_, best)| {
                summary.miner_damage > best.miner_damage
                    || (summary.miner_damage == best.miner_damage
                        && summary.maximum_strain_amplitude > best.maximum_strain_amplitude)
            }) {
                selected = Some((index, summary));
            }
        }
        if cancelled() {
            return Err(CriticalPlaneMinerError::Cancelled);
        }
        let (selected_plane, summary) =
            selected.ok_or(CriticalPlaneMinerError::Invalid("missing search planes"))?;
        let result = CriticalPlaneMinerResult {
            samples: self.samples,
            evaluated_planes: self.planes.len(),
            selected_plane,
            normal_in_frame: self.planes[selected_plane].normal,
            summary,
            projection_work: self.projection_work,
        };
        self.counters = candidate;
        self.result = Some(result.clone());
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fatigue::FatigueScope;
    use crate::fatigue::strain_life::CmbParameters;

    fn curve() -> Arc<StrainLifeCurve> {
        Arc::new(
            StrainLifeCurve::new(
                FatigueScope {
                    source_reference: "synthetic critical-plane Miner fixture".into(),
                    material: "synthetic coupon".into(),
                    process: "fixture".into(),
                    surface_and_geometry: "smooth coupon".into(),
                    environment: "fixed".into(),
                    load_basis: "variable-amplitude projected normal strain".into(),
                    temperature_k: 293.15,
                },
                CmbParameters {
                    youngs_modulus_pa: 200e9,
                    fatigue_strength_coefficient_pa: 1e9,
                    fatigue_strength_exponent: -0.5,
                    fatigue_ductility_coefficient: 0.002,
                    fatigue_ductility_exponent: -0.5,
                },
                0.5,
                1e8,
            )
            .unwrap(),
        )
    }

    fn limits(samples: u64, planes: usize) -> CriticalPlaneMinerLimits {
        CriticalPlaneMinerLimits {
            rainflow: RainflowLimits {
                maximum_samples: samples,
                maximum_residue: 64,
            },
            maximum_planes: planes,
            maximum_projection_work: samples * planes as u64,
        }
    }

    fn axial(strain: f64, stress: f64) -> TensorFatigueSample {
        TensorFatigueSample {
            total_mechanical_strain: [strain, 0.0, 0.0, 0.0, 0.0, 0.0],
            stress_pa: [stress, 0.0, 0.0, 0.0, 0.0, 0.0],
        }
    }

    fn close(a: f64, b: f64) {
        assert!(
            (a - b).abs() <= 8e-12 * a.abs().max(b.abs()).max(1e-30),
            "{a:e} != {b:e}"
        );
    }

    #[test]
    fn uniaxial_variable_amplitude_matches_the_original_fixed_plane_counter() {
        let history = [
            axial(-0.001, -200e6),
            axial(0.001, 200e6),
            axial(-0.0005, -100e6),
            axial(0.0015, 300e6),
            axial(-0.001, -200e6),
        ];
        let mut critical = CriticalPlaneMiner::new(
            curve(),
            StrainLifeMeanStress::SmithWatsonTopper,
            vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            limits(32, 8),
        )
        .unwrap();
        let mut fixed = PairedStrainLifeMiner::new(
            curve(),
            StrainLifeMeanStress::SmithWatsonTopper,
            RainflowLimits {
                maximum_samples: 32,
                maximum_residue: 64,
            },
        )
        .unwrap();
        for sample in history {
            critical.push(sample).unwrap();
            fixed
                .push(StrainStressSample {
                    mechanical_strain: sample.total_mechanical_strain[0],
                    normal_stress_pa: sample.stress_pa[0],
                })
                .unwrap();
        }
        let expected = fixed.finish().unwrap();
        let actual = critical.finish().unwrap();
        assert_eq!(actual.selected_plane, 0);
        assert_eq!(actual.summary, expected);
        assert_eq!(actual.samples, history.len() as u64);
        assert_eq!(actual.projection_work, 2 * history.len() as u64);
    }

    #[test]
    fn oblique_axial_history_selects_the_declared_physical_plane() {
        let n = [1.0 / 3.0_f64.sqrt(), 1.0 / 3.0_f64.sqrt(), 1.0 / 3.0_f64.sqrt()];
        let tensor = |strain: f64, stress: f64| {
            let weights = [
                n[0] * n[0],
                n[1] * n[1],
                n[2] * n[2],
                n[0] * n[1],
                n[1] * n[2],
                n[2] * n[0],
            ];
            TensorFatigueSample {
                total_mechanical_strain: weights.map(|w| w * strain),
                stress_pa: weights.map(|w| w * stress),
            }
        };
        let mut miner = CriticalPlaneMiner::new(
            curve(),
            StrainLifeMeanStress::Uncorrected,
            vec![[1.0, 0.0, 0.0], n, [0.0, 0.0, 1.0]],
            limits(16, 8),
        )
        .unwrap();
        for sign in [-1.0, 1.0, -0.5, 1.5, -1.0] {
            miner.push(tensor(0.001 * sign, 200e6 * sign)).unwrap();
        }
        let result = miner.finish().unwrap();
        assert_eq!(result.selected_plane, 1);
        assert!(result.summary.miner_damage > miner.plane_summaries()[0].miner_damage);
    }

    #[test]
    fn arbitrary_multiaxial_phase_lag_is_explicitly_refused_for_paired_mean_stress() {
        let samples = [
            TensorFatigueSample {
                total_mechanical_strain: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                stress_pa: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            },
            TensorFatigueSample {
                total_mechanical_strain: [0.001, 0.0, 0.0, 0.0, 0.0, 0.0],
                stress_pa: [-100e6, 0.0, 0.0, 0.0, 0.0, 0.0],
            },
        ];
        let mut paired = CriticalPlaneMiner::new(
            curve(),
            StrainLifeMeanStress::Morrow,
            vec![[1.0, 0.0, 0.0]],
            limits(8, 2),
        )
        .unwrap();
        paired.push(samples[0]).unwrap();
        let before = paired.plane_summaries();
        assert!(matches!(
            paired.push(samples[1]),
            Err(CriticalPlaneMinerError::Fatigue {
                error: FatigueError::UnsupportedPairedHistory { .. },
                ..
            })
        ));
        assert_eq!(paired.sample_count(), 1);
        assert_eq!(paired.plane_summaries(), before);

        let mut uncorrected = CriticalPlaneMiner::new(
            curve(),
            StrainLifeMeanStress::Uncorrected,
            vec![[1.0, 0.0, 0.0]],
            limits(8, 2),
        )
        .unwrap();
        uncorrected.push_block(&samples, &mut || false).unwrap();
        assert_eq!(uncorrected.sample_count(), 2);
    }

    #[test]
    fn chunk_boundaries_do_not_change_multiplane_damage() {
        let history = [
            axial(-0.001, -200e6),
            axial(0.001, 200e6),
            axial(-0.00025, -50e6),
            axial(0.00075, 150e6),
            axial(-0.00125, -250e6),
            axial(0.001, 200e6),
        ];
        let make = || {
            CriticalPlaneMiner::new(
                curve(),
                StrainLifeMeanStress::Uncorrected,
                vec![[1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
                limits(32, 8),
            )
            .unwrap()
        };
        let mut whole = make();
        whole.push_block(&history, &mut || false).unwrap();
        let expected = whole.finish().unwrap();
        for split in 0..=history.len() {
            let mut chunked = make();
            chunked
                .push_block(&history[..split], &mut || false)
                .unwrap();
            chunked
                .push_block(&history[split..], &mut || false)
                .unwrap();
            assert_eq!(chunked.finish().unwrap(), expected);
            assert_eq!(chunked.plane_summaries(), whole.plane_summaries());
        }
    }

    #[test]
    fn projection_and_sample_budgets_are_global_and_transactional() {
        let mut miner = CriticalPlaneMiner::new(
            curve(),
            StrainLifeMeanStress::Uncorrected,
            vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            CriticalPlaneMinerLimits {
                rainflow: RainflowLimits {
                    maximum_samples: 3,
                    maximum_residue: 8,
                },
                maximum_planes: 2,
                maximum_projection_work: 6,
            },
        )
        .unwrap();
        for x in [-0.001, 0.001, -0.001] {
            miner.push(axial(x, 200e9 * x)).unwrap();
        }
        let before = miner.plane_summaries();
        assert_eq!(
            miner.push(axial(0.0, 0.0)),
            Err(CriticalPlaneMinerError::Budget)
        );
        assert_eq!(miner.sample_count(), 3);
        assert_eq!(miner.plane_summaries(), before);
        assert_eq!(miner.finish().unwrap().projection_work, 6);
    }

    #[test]
    fn finish_is_atomic_and_idempotent_across_plane_specific_life_refusal() {
        let mut miner = CriticalPlaneMiner::new(
            curve(),
            StrainLifeMeanStress::SmithWatsonTopper,
            vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            limits(16, 4),
        )
        .unwrap();
        for sign in [-1.0, 1.0, -1.0] {
            miner
                .push(TensorFatigueSample {
                    total_mechanical_strain: [0.001 * sign, 0.0005 * sign, 0.0, 0.0, 0.0, 0.0],
                    stress_pa: [200e6 * sign, -150e6 + 50e6 * sign, 0.0, 0.0, 0.0, 0.0],
                })
                .unwrap();
        }
        let before = miner.plane_summaries();
        assert!(matches!(
            miner.finish(),
            Err(CriticalPlaneMinerError::Fatigue {
                error: FatigueError::NoTensileStress { .. },
                ..
            })
        ));
        assert_eq!(miner.result(), None);
        assert_eq!(miner.plane_summaries(), before);

        let mut clean = CriticalPlaneMiner::new(
            curve(),
            StrainLifeMeanStress::Uncorrected,
            vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            limits(16, 4),
        )
        .unwrap();
        for sign in [-1.0, 1.0, -1.0] {
            clean.push(axial(0.001 * sign, 200e6 * sign)).unwrap();
        }
        let result = clean.finish().unwrap();
        assert_eq!(clean.finish().unwrap(), result);
        assert_eq!(clean.push(axial(0.0, 0.0)), Err(CriticalPlaneMinerError::Finished));
        close(result.summary.maximum_strain_amplitude, 0.001);
    }
}
