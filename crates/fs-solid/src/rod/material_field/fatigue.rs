//! Signed stress histories at fixed material points of the original Cosserat rod.
//!
//! Native intrinsic-relative strains supply axial/bending stress; the original
//! circular-section law supplies torsional shear. Each layer retains its own E,
//! G, source state and integrated thermal expansion. Directions and points are
//! fixed in the co-rotating material frame, never moved to an instantaneous
//! stress maximum. Finite rigid rotation therefore does not create a new cycle.
//!
//! Scalar S-N rainflow and representative-cycle Findley evaluation reuse the
//! material owners. Transverse-force shear is either explicitly unassessed or
//! recovered through an opt-in circular Saint-Venant owner, homogeneous or
//! bonded layered with a common original Poisson ratio; Q/A is
//! never used as a local stress field. No plasticity, variable-temperature calibration transfer,
//! notch correction, life qualification or fatigue-to-stiffness/heat law is
//! inferred. Caller-selected load endpoints must resolve physical extrema.

use std::sync::Arc;

use fs_exec::CancelGate;
use fs_material::fatigue::critical_plane::findley::{
    FindleyCalibration, FindleyCycle, FindleyCycleResult, FindleyError, FindleyLimits,
    StressClosureTolerance,
};
use fs_material::fatigue::{
    FatigueError, FatigueScope, FatigueSummary, MeanStressCorrection, NormalStressDirection,
    RainflowLimits, RainflowMiner, SnCurve,
};
use fs_material::section::ConcentricShearSection;
use fs_material::section::transverse_stress::CircularTransverseShearRecovery;
use fs_material::section::transverse_stress::TransverseShearStressPoint;
use fs_material::section::transverse_stress::layered::LayeredCircularTransverseShearRecovery;
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_qty::semantic::{QuantityKind, ValueForm};
use fs_qty::{QuantitySpec, Temperature};

use super::{SpatialMaterialRod, SpatialRodError, checkpoint};
use crate::rod::{TipLoad, WorldTipLoad};

/// A permanent point in one original material annulus, in local metres.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RodFatiguePoint {
    pub segment: usize,
    pub layer: usize,
    /// Local [y,z]; +x is the member axis. Both sides of an interface may be
    /// observed separately by explicitly selecting their different layers.
    pub coordinates_yz_m: [f64; 2],
}

/// Signed normal-stress life on one fixed material plane at one fixed point.
#[derive(Debug, Clone)]
pub struct RodFatigueGauge {
    pub point: RodFatiguePoint,
    pub normal_in_material_frame: NormalStressDirection,
    /// Exact native absolute-temperature axis; no alias or ambient is guessed.
    pub temperature_axis: String,
    pub curve: Arc<SnCurve>,
    pub mean_stress: MeanStressCorrection,
}

/// Shear-sensitive fatigue-limit usage of an explicitly closed stress window.
#[derive(Debug, Clone)]
pub struct RodFindleyGauge {
    pub point: RodFatiguePoint,
    pub temperature_axis: String,
    pub calibration: Arc<FindleyCalibration>,
    pub normals_in_material_frame: Vec<[f64; 3]>,
    pub limits: FindleyLimits,
    pub closure: StressClosureTolerance,
}

/// Observation storage and physical applicability, independent of native Newton.
#[derive(Debug, Clone, Copy)]
pub struct RodFatigueConfig {
    /// Global accepted sample ceiling, including the initial pose exactly once.
    pub history: RainflowLimits,
    /// Combined scalar and Findley gauge count.
    pub maximum_gauges: usize,
    /// Aggregate declared scalar reversal slots per transaction copy.
    pub maximum_total_residue: usize,
    /// Sum of the Findley windows' declared retained tensor-sample ceilings.
    pub maximum_total_findley_samples: usize,
    /// Whole selected annulus axial/bending strain bound, in (0, 0.05].
    /// Both total and mechanical strain are checked; thermal bounds also apply.
    pub maximum_absolute_axial_strain: f64,
    /// Whole selected annulus engineering-shear bound, in (0, 0.05]. With
    /// transverse recovery, the sum of the separate transverse/torsion peaks
    /// conservatively bounds their combined physical strain everywhere.
    pub maximum_engineering_shear: f64,
    /// Maximum omitted native transverse-force magnitude [N], finite and >=0.
    /// Zero is strict. A nonzero allowance is the caller's explicit decision
    /// that transverse shear may be neglected, not a local stress-error bound.
    pub maximum_unassessed_transverse_force_n: f64,
}

/// Latest actual local observation; no von Mises or spatial maximum substitution.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RodFatigueSample {
    pub point: RodFatiguePoint,
    pub axial_mechanical_strain: f64,
    /// Tensor-Voigt [xx,yy,zz,xy,yz,zx] [Pa], in the local material frame.
    /// Axial/bending and circular torsion, plus transverse Saint-Venant shear
    /// when the observation was explicitly constructed with that model.
    pub stress_pa: [f64; 6],
    /// Native transverse resultants [Qy,Qz] [N] whose local stress is absent.
    /// Exactly zero when transverse recovery accounts for this contribution;
    /// otherwise even a small force is retained without a significance cutoff.
    pub unassessed_transverse_force_n: [f64; 2],
}

#[derive(Debug, Clone, PartialEq)]
pub enum RodFatigueError {
    Invalid(&'static str),
    Budget,
    /// Combined declaration order: scalar gauges, then Findley gauges.
    SourceMismatch {
        gauge: usize,
        what: &'static str,
    },
    Gauge {
        gauge: usize,
        error: FatigueError,
    },
    /// Index within the Findley bank.
    Findley {
        gauge: usize,
        error: FindleyError,
    },
    Mechanics(SpatialRodError),
    Cancelled,
    Finished,
}
impl core::fmt::Display for RodFatigueError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native rod fatigue refused: {self:?}")
    }
}
impl std::error::Error for RodFatigueError {}
impl From<SpatialRodError> for RodFatigueError {
    fn from(error: SpatialRodError) -> Self {
        match error {
            SpatialRodError::Cancelled => Self::Cancelled,
            error => Self::Mechanics(error),
        }
    }
}

#[derive(Debug, Clone)]
struct PointBinding {
    point: RodFatiguePoint,
    segment_slot: usize,
    youngs_modulus_pa: f64,
    thermal_strain: f64,
    radius_m: f64,
    strain_limit: f64,
}
#[derive(Debug, Clone)]
struct SegmentBinding {
    segment: usize,
    section: Arc<ConcentricShearSection>,
    transverse: Option<TransverseRecovery>,
}

#[derive(Clone, Copy)]
pub(super) enum TransverseModel {
    Omitted,
    Homogeneous,
    LayeredCommonPoisson,
}

#[derive(Debug, Clone)]
enum TransverseRecovery {
    Homogeneous(CircularTransverseShearRecovery),
    Layered(LayeredCircularTransverseShearRecovery),
}
impl TransverseRecovery {
    fn evaluate(
        &self,
        layer: usize,
        point: [f64; 2],
        force: [f64; 2],
        bound: f64,
    ) -> Result<TransverseShearStressPoint, fs_material::state_point::MaterialStatePointError> {
        match self {
            Self::Homogeneous(source) => source.evaluate(layer, point, force, bound),
            Self::Layered(source) => source.evaluate(layer, point, force, bound),
        }
    }
}
#[derive(Debug, Clone, Copy, Default)]
struct SegmentSample {
    gamma: [f64; 3],
    curvature: [f64; 3],
    transverse_force_n: [f64; 2],
}
#[derive(Debug, Clone)]
struct NormalGauge {
    direction: NormalStressDirection,
    counter: RainflowMiner,
}

/// Reusable observation bank, separate from native mechanical advancement.
/// Immutable admissions and selected-segment recovery are prepared once. The
/// native path is its acceptance owner; queries never advance or re-solve it.
#[derive(Debug, Clone)]
pub struct RodFatigueObservation {
    // A motion bridge may pin the immutable prepared source and admit it in O(1).
    pub(super) motion_source: Option<Arc<SpatialMaterialRod>>,
    segments: Vec<SegmentBinding>,
    points: Vec<PointBinding>,
    normal: Vec<NormalGauge>,
    findley: Vec<FindleyCycle>,
    samples: Vec<RodFatigueSample>,
    recovery: Vec<SegmentSample>,
    config: RodFatigueConfig,
    sample_count: u64,
    finished: bool,
}

struct ObservationCheckpoint {
    normal: Vec<NormalGauge>,
    samples: Vec<RodFatigueSample>,
    sample_count: u64,
    findley: Vec<(u64, bool)>,
}

impl RodFatigueObservation {
    /// Bind original source sections and observe the supplied accepted pose
    /// once. A preloaded specimen starts at its actual stress, not fictitious zero.
    pub fn new(
        rod: &SpatialMaterialRod,
        gauges: Vec<RodFatigueGauge>,
        findley: Vec<RodFindleyGauge>,
        config: RodFatigueConfig,
        gate: &CancelGate,
    ) -> Result<Self, RodFatigueError> {
        let mut bank = Self::prepare(rod, gauges, findley, config, TransverseModel::Omitted, gate)?;
        bank.observe_native(rod, gate)?;
        Ok(bank)
    }

    /// Observe the same native history with transverse-force shear included at
    /// every fixed point. Each selected segment must retain an explicitly
    /// isotropic, elastically homogeneous solid circle or concentric tube.
    /// The source-owned Saint-Venant field satisfies the actual section force
    /// and traction-free contours; it does not change the native authored
    /// shear correction, mechanics or energy. End effects and cross-sectional
    /// inertia remain outside this quasistatic section recovery.
    pub fn new_with_transverse_shear(
        rod: &SpatialMaterialRod,
        gauges: Vec<RodFatigueGauge>,
        findley: Vec<RodFindleyGauge>,
        config: RodFatigueConfig,
        gate: &CancelGate,
    ) -> Result<Self, RodFatigueError> {
        let mut bank = Self::prepare(
            rod,
            gauges,
            findley,
            config,
            TransverseModel::Homogeneous,
            gate,
        )?;
        bank.observe_native(rod, gate)?;
        Ok(bank)
    }

    /// Resolve bonded circular layers with different original E/G but one
    /// common explicitly isotropic Poisson ratio. Exact warping/traction
    /// transmission supplies each layer's own stress; no homogenized G is
    /// inserted into its fatigue tensor. The original mechanics and atomic
    /// history retain their existing owners and small-strain bounds.
    pub fn new_with_layered_transverse_shear(
        rod: &SpatialMaterialRod,
        gauges: Vec<RodFatigueGauge>,
        findley: Vec<RodFindleyGauge>,
        config: RodFatigueConfig,
        gate: &CancelGate,
    ) -> Result<Self, RodFatigueError> {
        let mut bank = Self::prepare(
            rod,
            gauges,
            findley,
            config,
            TransverseModel::LayeredCommonPoisson,
            gate,
        )?;
        bank.observe_native(rod, gate)?;
        Ok(bank)
    }

    pub(super) fn prepare(
        rod: &SpatialMaterialRod,
        gauges: Vec<RodFatigueGauge>,
        findley: Vec<RodFindleyGauge>,
        config: RodFatigueConfig,
        transverse_shear: TransverseModel,
        gate: &CancelGate,
    ) -> Result<Self, RodFatigueError> {
        poll(gate)?;
        let count = gauges
            .len()
            .checked_add(findley.len())
            .ok_or(RodFatigueError::Budget)?;
        if count == 0
            || count > config.maximum_gauges
            || config.history.maximum_samples == 0
            || config.history.maximum_residue < 3
            || !valid_strain_limit(config.maximum_absolute_axial_strain)
            || !valid_strain_limit(config.maximum_engineering_shear)
            || !config.maximum_unassessed_transverse_force_n.is_finite()
            || config.maximum_unassessed_transverse_force_n < 0.0
        {
            return Err(RodFatigueError::Invalid(
                "declare nonempty bounded gauges and small-strain domains",
            ));
        }
        if gauges
            .len()
            .checked_mul(config.history.maximum_residue)
            .is_none_or(|n| n > config.maximum_total_residue)
        {
            return Err(RodFatigueError::Budget);
        }
        let retained = findley
            .iter()
            .try_fold(0usize, |sum, gauge| {
                usize::try_from(gauge.limits.maximum_samples)
                    .ok()
                    .and_then(|n| sum.checked_add(n))
            })
            .ok_or(RodFatigueError::Budget)?;
        if retained > config.maximum_total_findley_samples {
            return Err(RodFatigueError::Budget);
        }
        let mut bank = Self {
            motion_source: None,
            segments: Vec::new(),
            points: Vec::with_capacity(count),
            normal: Vec::with_capacity(gauges.len()),
            findley: Vec::with_capacity(findley.len()),
            samples: Vec::with_capacity(count),
            recovery: Vec::new(),
            config,
            sample_count: 0,
            finished: false,
        };
        for (index, gauge) in gauges.into_iter().enumerate() {
            let material = layer_material(rod, gauge.point)?;
            let named = material.material();
            admit_scope(
                material,
                &gauge.temperature_axis,
                gauge.curve.scope(),
                &format!("{}; {}", named.chemistry, named.phase),
                index,
            )?;
            if let Some(source) = gauge.curve.source() {
                if gauge.temperature_axis != source.request().temperature_axis
                    || !source.matches_material_state(
                        material.card_identity(),
                        material.material(),
                        &material.query(),
                    )
                {
                    return Err(source_error(
                        index,
                        "S-N source differs from the actual layer card, state or typed context",
                    ));
                }
            }
            bank.add_point(rod, gauge.point, transverse_shear, gate)?;
            bank.normal.push(NormalGauge {
                direction: gauge.normal_in_material_frame,
                counter: RainflowMiner::new(gauge.curve, gauge.mean_stress, config.history)
                    .map_err(|error| gauge_error(index, error))?,
            });
        }
        for (index, gauge) in findley.into_iter().enumerate() {
            let source_index = bank.normal.len() + index;
            let material = layer_material(rod, gauge.point)?;
            admit_scope(
                material,
                &gauge.temperature_axis,
                gauge.calibration.scope(),
                &material.material().chemistry,
                source_index,
            )?;
            if let Some(source) = gauge.calibration.source() {
                if gauge.temperature_axis != source.request().temperature_axis
                    || !source.matches_material_state(
                        material.card_identity(),
                        material.material(),
                        &material.query(),
                    )
                {
                    return Err(source_error(
                        source_index,
                        "Findley source differs from the actual layer card, state or typed context",
                    ));
                }
            }
            bank.add_point(rod, gauge.point, transverse_shear, gate)?;
            bank.findley.push(
                FindleyCycle::new(
                    gauge.calibration,
                    gauge.normals_in_material_frame,
                    gauge.limits,
                    gauge.closure,
                )
                .map_err(|error| findley_error(index, error))?,
            );
        }
        bank.recovery
            .resize(bank.segments.len(), SegmentSample::default());
        Ok(bank)
    }

    fn add_point(
        &mut self,
        rod: &SpatialMaterialRod,
        point: RodFatiguePoint,
        transverse_shear: TransverseModel,
        gate: &CancelGate,
    ) -> Result<(), RodFatigueError> {
        poll(gate)?;
        let section = &rod.sections[point.segment]; // Scope admission checked both indices.
        section
            .torsional_stress_at_twist_bounded(
                point.layer,
                point.coordinates_yz_m,
                0.0,
                self.config.maximum_engineering_shear,
            )
            .map_err(SpatialRodError::from)?;
        let segment_slot = match self
            .segments
            .iter()
            .position(|s| s.segment == point.segment)
        {
            Some(slot) => slot,
            None => {
                let section = Arc::new(section.clone());
                let transverse = match transverse_shear {
                    TransverseModel::Omitted => None,
                    TransverseModel::Homogeneous => Some(TransverseRecovery::Homogeneous(
                        CircularTransverseShearRecovery::try_new(Arc::clone(&section))
                            .map_err(SpatialRodError::from)?,
                    )),
                    TransverseModel::LayeredCommonPoisson => Some(TransverseRecovery::Layered(
                        LayeredCircularTransverseShearRecovery::try_new(Arc::clone(&section))
                            .map_err(SpatialRodError::from)?,
                    )),
                };
                self.segments.push(SegmentBinding {
                    segment: point.segment,
                    section,
                    transverse,
                });
                self.segments.len() - 1
            }
        };
        let thermal_strain = rod.thermal_sections().map_or(0.0, |sections| {
            sections[point.segment].expansion_paths()[point.layer].free_linear_strain()
        });
        let strain_limit =
            rod.thermal
                .as_ref()
                .map_or(self.config.maximum_absolute_axial_strain, |thermal| {
                    thermal
                        .linear_strain_limit()
                        .min(self.config.maximum_absolute_axial_strain)
                });
        self.points.push(PointBinding {
            point,
            segment_slot,
            thermal_strain,
            strain_limit,
            youngs_modulus_pa: section.elastic().layer_integrals()[point.layer].young_modulus_pa(),
            radius_m: section.elastic().layers()[point.layer].outer_radius_m,
        });
        self.samples.push(RodFatigueSample {
            point,
            axial_mechanical_strain: 0.0,
            stress_pa: [0.0; 6],
            unassessed_transverse_force_n: [0.0; 2],
        });
        Ok(())
    }

    pub fn sample_count(&self) -> u64 {
        self.sample_count
    }
    pub fn is_finished(&self) -> bool {
        self.finished
    }
    pub fn normal_samples(&self) -> &[RodFatigueSample] {
        &self.samples[..self.normal.len()]
    }
    pub fn findley_samples(&self) -> &[RodFatigueSample] {
        &self.samples[self.normal.len()..]
    }
    pub fn curve(&self, gauge: usize) -> Option<&SnCurve> {
        self.normal.get(gauge).map(|g| g.counter.curve())
    }
    pub fn findley_calibration(&self, gauge: usize) -> Option<&FindleyCalibration> {
        self.findley.get(gauge).map(FindleyCycle::calibration)
    }
    /// Original source in combined scalar-then-Findley declaration order.
    pub fn material_source(&self, gauge: usize) -> Option<&ResolvedMaterialStatePoint> {
        let p = self.points.get(gauge)?;
        Some(&self.segments[p.segment_slot].section.elastic().layers()[p.point.layer].material)
    }
    pub fn summaries(&self) -> Vec<FatigueSummary> {
        self.normal.iter().map(|g| g.counter.summary()).collect()
    }
    pub fn findley_results(&self) -> Vec<Option<FindleyCycleResult>> {
        self.findley.iter().map(|g| g.result().copied()).collect()
    }

    /// Include pending endpoint half cycles without finalizing a live counter,
    /// adding a native observation, or closing a representative Findley window.
    pub fn preview_fatigue(
        &self,
        gate: &CancelGate,
    ) -> Result<Vec<FatigueSummary>, RodFatigueError> {
        poll(gate)?;
        self.normal
            .iter()
            .enumerate()
            .map(|(index, gauge)| {
                gauge
                    .counter
                    .preview_finish_with_cancel(&mut || cancelled(gate))
                    .map_err(|error| gauge_error(index, error))
            })
            .collect()
    }

    /// Close all Findley windows together. Scalar rainflow stays open.
    pub fn close_findley_cycle(
        &mut self,
        gate: &CancelGate,
    ) -> Result<Vec<FindleyCycleResult>, RodFatigueError> {
        self.require_findley(gate)?;
        self.finish_findley(gate)
    }

    /// Begin a fresh window at the latest actual sample. This explicit operation
    /// may discard an open warm-up window; scalar histories and global budget
    /// remain continuous. Prior stress traces are not cloned or replayed.
    pub fn restart_findley_cycle(&mut self, gate: &CancelGate) -> Result<(), RodFatigueError> {
        self.require_findley(gate)?;
        let mut windows = Vec::with_capacity(self.findley.len());
        for (index, old) in self.findley.iter().enumerate() {
            poll(gate)?;
            let mut fresh = old.fresh().map_err(|error| findley_error(index, error))?;
            fresh
                .push_with_cancel(
                    self.samples[self.normal.len() + index].stress_pa,
                    &mut || cancelled(gate),
                )
                .map_err(|error| findley_error(index, error))?;
            windows.push(fresh);
        }
        poll(gate)?;
        self.findley = windows;
        Ok(())
    }

    /// Finalize scalar half cycles and all declared Findley windows atomically.
    /// A nonclosed representative window leaves the complete observation live.
    pub fn finish(&mut self, gate: &CancelGate) -> Result<Vec<FatigueSummary>, RodFatigueError> {
        poll(gate)?;
        if self.finished {
            return Ok(self.summaries());
        }
        let mut normal = self.normal.clone();
        let mut results = Vec::with_capacity(normal.len());
        for (index, gauge) in normal.iter_mut().enumerate() {
            results.push(
                gauge
                    .counter
                    .finish_with_cancel(&mut || cancelled(gate))
                    .map_err(|error| gauge_error(index, error))?,
            );
        }
        self.finish_findley(gate)?;
        self.normal = normal;
        self.finished = true;
        Ok(results)
    }

    fn require_findley(&self, gate: &CancelGate) -> Result<(), RodFatigueError> {
        poll(gate)?;
        if self.finished {
            return Err(RodFatigueError::Finished);
        }
        if self.findley.is_empty() {
            return Err(RodFatigueError::Invalid("no Findley gauges were declared"));
        }
        Ok(())
    }

    fn admit_steps(&self, additional: usize, gate: &CancelGate) -> Result<(), RodFatigueError> {
        poll(gate)?;
        if self.finished {
            return Err(RodFatigueError::Finished);
        }
        let additional = u64::try_from(additional).map_err(|_| RodFatigueError::Budget)?;
        if additional == 0
            || self
                .sample_count
                .checked_add(additional)
                .is_none_or(|n| n > self.config.history.maximum_samples)
        {
            return Err(RodFatigueError::Budget);
        }
        for (index, gauge) in self.findley.iter().enumerate() {
            if gauge.result().is_some() {
                return Err(findley_error(index, FindleyError::Finished));
            }
            let count = gauge
                .sample_count()
                .checked_add(additional)
                .ok_or(RodFatigueError::Budget)?;
            if count > gauge.limits().maximum_samples
                || count
                    .checked_mul(gauge.plane_count() as u64)
                    .is_none_or(|n| n > gauge.limits().maximum_projection_work)
            {
                return Err(findley_error(index, FindleyError::Budget));
            }
        }
        Ok(())
    }

    fn checkpoint(&self) -> ObservationCheckpoint {
        ObservationCheckpoint {
            normal: self.normal.clone(),
            samples: self.samples.clone(),
            sample_count: self.sample_count,
            findley: self
                .findley
                .iter()
                .map(|g| (g.sample_count(), g.result().is_none()))
                .collect(),
        }
    }
    fn restore(&mut self, old: ObservationCheckpoint) {
        self.normal = old.normal;
        self.samples = old.samples;
        self.sample_count = old.sample_count;
        for (gauge, (count, was_open)) in self.findley.iter_mut().zip(old.findley) {
            if was_open {
                gauge
                    .discard_after(count)
                    .expect("native observation only appends Findley samples");
            }
        }
    }

    fn finish_findley(
        &mut self,
        gate: &CancelGate,
    ) -> Result<Vec<FindleyCycleResult>, RodFatigueError> {
        let open: Vec<_> = self.findley.iter().map(|g| g.result().is_none()).collect();
        let outcome = (|| {
            let mut results = Vec::with_capacity(self.findley.len());
            for (index, gauge) in self.findley.iter_mut().enumerate() {
                results.push(
                    gauge
                        .finish_with_cancel(&mut || cancelled(gate))
                        .map_err(|error| findley_error(index, error))?,
                );
            }
            poll(gate)?;
            Ok(results)
        })();
        if outcome.is_err() {
            for (gauge, was_open) in self.findley.iter_mut().zip(open) {
                if was_open {
                    gauge
                        .discard_after(gauge.sample_count())
                        .expect("finalization retains all Findley samples");
                }
            }
        }
        outcome
    }

    // Caller owns the complete native/observation transaction. These same
    // selected segment strains can be supplied by an accepted motion pose.
    fn observe_native(
        &mut self,
        rod: &SpatialMaterialRod,
        gate: &CancelGate,
    ) -> Result<(), RodFatigueError> {
        self.admit_steps(1, gate)?;
        for (binding, sample) in self.segments.iter().zip(&mut self.recovery) {
            poll(gate)?;
            let (gamma, curvature) = rod.segment_strains(binding.segment)?;
            let (force, _) = rod.resultants_from_strains(binding.segment, gamma, curvature)?;
            *sample = SegmentSample {
                gamma,
                curvature,
                transverse_force_n: [force[1], force[2]],
            };
        }
        self.observe_recovered(gate)
    }

    fn observe_recovered(&mut self, gate: &CancelGate) -> Result<(), RodFatigueError> {
        for (point, sample) in self.points.iter().zip(&mut self.samples) {
            poll(gate)?;
            let native = self.recovery[point.segment_slot];
            let binding = &self.segments[point.segment_slot];
            let transverse =
                fs_math::det::hypot(native.transverse_force_n[0], native.transverse_force_n[1]);
            if !transverse.is_finite()
                || (binding.transverse.is_none()
                    && transverse > self.config.maximum_unassessed_transverse_force_n)
            {
                return Err(RodFatigueError::Invalid(
                    "native transverse force exceeds the declared omitted-shear allowance",
                ));
            }
            let bending =
                point.radius_m * fs_math::det::hypot(native.curvature[1], native.curvature[2]);
            let mechanical_center = native.gamma[0] - point.thermal_strain;
            if !bending.is_finite()
                || !mechanical_center.is_finite()
                || native.gamma[0].abs() + bending > point.strain_limit
                || mechanical_center.abs() + bending > point.strain_limit
                || point.thermal_strain.abs() > point.strain_limit
            {
                return Err(RodFatigueError::Invalid(
                    "selected rod layer exceeds its total/mechanical axial-bending strain domain",
                ));
            }
            let [y, z] = point.point.coordinates_yz_m;
            let strain = mechanical_center + native.curvature[1] * z - native.curvature[2] * y;
            let stress = point.youngs_modulus_pa * strain;
            if !strain.is_finite() || !stress.is_finite() || strain != 0.0 && stress == 0.0 {
                return Err(RodFatigueError::Invalid(
                    "rod point normal stress is unrepresentable",
                ));
            }
            let torsion = binding
                .section
                .torsional_stress_at_twist_bounded(
                    point.point.layer,
                    point.point.coordinates_yz_m,
                    native.curvature[0],
                    self.config.maximum_engineering_shear,
                )
                .map_err(SpatialRodError::from)?;
            let mut stress_pa = torsion.stress_tensor_pa();
            stress_pa[0] = stress;
            let unassessed_transverse_force_n = if let Some(recovery) = &binding.transverse {
                let recovered = recovery
                    .evaluate(
                        point.point.layer,
                        point.point.coordinates_yz_m,
                        native.transverse_force_n,
                        self.config.maximum_engineering_shear,
                    )
                    .map_err(SpatialRodError::from)?;
                let combined_bound = recovered.maximum_layer_engineering_shear
                    + native.curvature[0].abs() * point.radius_m;
                if !combined_bound.is_finite()
                    || combined_bound > self.config.maximum_engineering_shear
                {
                    return Err(RodFatigueError::Invalid(
                        "combined transverse and torsional strain exceeds the selected layer domain",
                    ));
                }
                stress_pa[3] += recovered.shear_stress_pa[0];
                stress_pa[5] += recovered.shear_stress_pa[1];
                if !stress_pa[3].is_finite() || !stress_pa[5].is_finite() {
                    return Err(RodFatigueError::Invalid(
                        "combined rod shear stress overflows",
                    ));
                }
                [0.0; 2]
            } else {
                native.transverse_force_n
            };
            *sample = RodFatigueSample {
                point: point.point,
                axial_mechanical_strain: strain,
                stress_pa,
                unassessed_transverse_force_n,
            };
        }
        for (index, gauge) in self.normal.iter_mut().enumerate() {
            let normal = gauge
                .direction
                .project(&self.samples[index].stress_pa)
                .map_err(|error| gauge_error(index, error))?;
            gauge
                .counter
                .push_with_cancel(normal, &mut || cancelled(gate))
                .map_err(|error| gauge_error(index, error))?;
        }
        for (index, gauge) in self.findley.iter_mut().enumerate() {
            gauge
                .push_with_cancel(
                    self.samples[self.normal.len() + index].stress_pa,
                    &mut || cancelled(gate),
                )
                .map_err(|error| findley_error(index, error))?;
        }
        poll(gate)?;
        self.sample_count += 1;
        Ok(())
    }

    // Only the source-pinned motion bridge may call this; the public caller
    // cannot substitute arbitrary strain tensors or a refreshed material model.
    pub(super) fn accept_motion_strains(
        &mut self,
        rod: &SpatialMaterialRod,
        mut recover: impl FnMut(usize) -> Result<([f64; 3], [f64; 3]), SpatialRodError>,
        gate: &CancelGate,
    ) -> Result<(), RodFatigueError> {
        self.admit_steps(1, gate)?;
        let old = self.checkpoint();
        let outcome = (|| {
            for (binding, sample) in self.segments.iter().zip(&mut self.recovery) {
                poll(gate)?;
                let (gamma, curvature) = recover(binding.segment)?;
                let (force, _) = rod.resultants_from_strains(binding.segment, gamma, curvature)?;
                *sample = SegmentSample {
                    gamma,
                    curvature,
                    transverse_force_n: [force[1], force[2]],
                };
            }
            self.observe_recovered(gate)
        })();
        if outcome.is_err() {
            self.restore(old);
        }
        outcome
    }
}

/// Actual native equilibrium and all fatigue banks share one acceptance boundary.
/// Only immutable mechanics/observation views are exposed, so no material refresh,
/// pose edit or unobserved solve can bypass a history. Clone preserves continuation.
#[derive(Debug, Clone)]
pub struct SpatialRodFatiguePath {
    rod: SpatialMaterialRod,
    observation: RodFatigueObservation,
}
impl SpatialRodFatiguePath {
    pub fn new(
        rod: SpatialMaterialRod,
        gauges: Vec<RodFatigueGauge>,
        findley: Vec<RodFindleyGauge>,
        config: RodFatigueConfig,
        gate: &CancelGate,
    ) -> Result<Self, RodFatigueError> {
        let observation = RodFatigueObservation::new(&rod, gauges, findley, config, gate)?;
        Ok(Self { rod, observation })
    }

    /// Native equilibrium cycles including the original section's resolved
    /// transverse shear; see [`RodFatigueObservation::new_with_transverse_shear`].
    pub fn new_with_transverse_shear(
        rod: SpatialMaterialRod,
        gauges: Vec<RodFatigueGauge>,
        findley: Vec<RodFindleyGauge>,
        config: RodFatigueConfig,
        gate: &CancelGate,
    ) -> Result<Self, RodFatigueError> {
        let observation =
            RodFatigueObservation::new_with_transverse_shear(&rod, gauges, findley, config, gate)?;
        Ok(Self { rod, observation })
    }

    /// Original equilibrium cycles with bonded common-nu layer transmission.
    pub fn new_with_layered_transverse_shear(
        rod: SpatialMaterialRod,
        gauges: Vec<RodFatigueGauge>,
        findley: Vec<RodFindleyGauge>,
        config: RodFatigueConfig,
        gate: &CancelGate,
    ) -> Result<Self, RodFatigueError> {
        let observation = RodFatigueObservation::new_with_layered_transverse_shear(
            &rod, gauges, findley, config, gate,
        )?;
        Ok(Self { rod, observation })
    }
    pub fn rod(&self) -> &SpatialMaterialRod {
        &self.rod
    }
    pub fn observation(&self) -> &RodFatigueObservation {
        &self.observation
    }
    pub fn sample_count(&self) -> u64 {
        self.observation.sample_count()
    }
    pub fn summaries(&self) -> Vec<FatigueSummary> {
        self.observation.summaries()
    }
    pub fn preview_fatigue(
        &self,
        gate: &CancelGate,
    ) -> Result<Vec<FatigueSummary>, RodFatigueError> {
        self.observation.preview_fatigue(gate)
    }
    pub fn close_findley_cycle(
        &mut self,
        gate: &CancelGate,
    ) -> Result<Vec<FindleyCycleResult>, RodFatigueError> {
        self.observation.close_findley_cycle(gate)
    }
    pub fn restart_findley_cycle(&mut self, gate: &CancelGate) -> Result<(), RodFatigueError> {
        self.observation.restart_findley_cycle(gate)
    }
    pub fn finish(&mut self, gate: &CancelGate) -> Result<Vec<FatigueSummary>, RodFatigueError> {
        self.observation.finish(gate)
    }

    /// Explicit chronological load targets with world force and body-held tip
    /// couple. Each target uses ONE original native equilibrium step; Newton
    /// trials never enter the history. Supply intermediate targets to resolve
    /// extrema/rotation branches. Return native residual histories per target.
    /// A refusal anywhere rolls back the complete requested target block.
    pub fn advance(
        &mut self,
        loads: &[TipLoad],
        tolerance: f64,
        gate: &CancelGate,
    ) -> Result<Vec<Vec<f64>>, RodFatigueError> {
        self.advance_inner(loads.len(), false, gate, |rod, i| {
            rod.solve_static(&loads[i], 1, tolerance, gate)
        })
    }

    /// Same native trajectory with force AND couple held in world axes.
    pub fn advance_world(
        &mut self,
        loads: &[WorldTipLoad],
        tolerance: f64,
        gate: &CancelGate,
    ) -> Result<Vec<Vec<f64>>, RodFatigueError> {
        self.advance_inner(loads.len(), false, gate, |rod, i| {
            rod.solve_static_world(&loads[i], 1, tolerance, gate)
        })
    }

    /// Body-couple targets plus Findley closure in the SAME transaction.
    pub fn advance_and_close_findley_cycle(
        &mut self,
        loads: &[TipLoad],
        tolerance: f64,
        gate: &CancelGate,
    ) -> Result<Vec<Vec<f64>>, RodFatigueError> {
        self.advance_inner(loads.len(), true, gate, |rod, i| {
            rod.solve_static(&loads[i], 1, tolerance, gate)
        })
    }

    /// World-couple targets plus Findley closure in the SAME transaction.
    pub fn advance_world_and_close_findley_cycle(
        &mut self,
        loads: &[WorldTipLoad],
        tolerance: f64,
        gate: &CancelGate,
    ) -> Result<Vec<Vec<f64>>, RodFatigueError> {
        self.advance_inner(loads.len(), true, gate, |rod, i| {
            rod.solve_static_world(&loads[i], 1, tolerance, gate)
        })
    }

    fn advance_inner(
        &mut self,
        count: usize,
        close: bool,
        gate: &CancelGate,
        mut solve: impl FnMut(&mut SpatialMaterialRod, usize) -> Result<Vec<f64>, SpatialRodError>,
    ) -> Result<Vec<Vec<f64>>, RodFatigueError> {
        self.observation.admit_steps(count, gate)?;
        if close {
            self.observation.require_findley(gate)?;
        }
        let old = self.observation.checkpoint();
        // Native statics mutates only the accepted inner pose. Immutable source
        // sections, thermal paths and reference inventory need no per-block copy.
        let old_pose = self.rod.rod.clone();
        let outcome = (|| {
            let mut histories = Vec::with_capacity(count);
            for i in 0..count {
                poll(gate)?;
                histories.push(solve(&mut self.rod, i)?);
                self.observation.observe_native(&self.rod, gate)?;
            }
            if close {
                self.observation.close_findley_cycle(gate)?;
            }
            poll(gate)?;
            Ok(histories)
        })();
        if outcome.is_err() {
            self.rod.rod = old_pose;
            self.observation.restore(old);
        }
        outcome
    }
}

fn layer_material(
    rod: &SpatialMaterialRod,
    point: RodFatiguePoint,
) -> Result<&ResolvedMaterialStatePoint, RodFatigueError> {
    Ok(&rod
        .sections
        .get(point.segment)
        .and_then(|s| s.elastic().layers().get(point.layer))
        .ok_or(RodFatigueError::Invalid(
            "fatigue segment/layer is outside the original rod",
        ))?
        .material)
}

/// `scope_material` is the name the calibration's source gives this layer's
/// material: an S-N source writes "chemistry; phase" (fs-material
/// `fatigue::source`, as the plate and frame S-N gauges read it); a Findley
/// source writes the chemistry.
fn admit_scope(
    material: &ResolvedMaterialStatePoint,
    axis: &str,
    scope: &FatigueScope,
    scope_material: &str,
    index: usize,
) -> Result<(), RodFatigueError> {
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
    if axis.trim().is_empty()
        || !absolute_temperature
        || material.state_coordinate(axis) != Some(scope.temperature_k)
        || scope.material != scope_material
        || material.material().process != scope.process
    {
        return Err(source_error(
            index,
            "fatigue needs the original layer material/process and typed absolute temperature",
        ));
    }
    Ok(())
}
fn valid_strain_limit(value: f64) -> bool {
    value.is_finite() && value > 0.0 && value <= 0.05
}
fn poll(gate: &CancelGate) -> Result<(), RodFatigueError> {
    checkpoint(gate).map_err(Into::into)
}
fn cancelled(gate: &CancelGate) -> bool {
    checkpoint(gate).is_err()
}
fn source_error(gauge: usize, what: &'static str) -> RodFatigueError {
    RodFatigueError::SourceMismatch { gauge, what }
}
fn gauge_error(gauge: usize, error: FatigueError) -> RodFatigueError {
    match error {
        FatigueError::Cancelled => RodFatigueError::Cancelled,
        error => RodFatigueError::Gauge { gauge, error },
    }
}
fn findley_error(gauge: usize, error: FindleyError) -> RodFatigueError {
    match error {
        FindleyError::Cancelled => RodFatigueError::Cancelled,
        error => RodFatigueError::Findley { gauge, error },
    }
}
