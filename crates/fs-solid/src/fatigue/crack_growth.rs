//! Chronological Mode-I Paris growth from the original elastic solid response.
//!
//! Each gauge declares an UNCRACKED nominal-stress location, a fixed crack-plane
//! normal and a constant geometry factor over a bounded crack-length interval.
//! The recovered native stress supplies the tensile normal traction; appreciable
//! in-plane shear traction refuses instead of silently entering a Mode-I law.
//! The mesh does not contain the crack, and growth changes neither stiffness,
//! thermal energy nor topology. Small-scale yielding and the geometry factor's
//! applicability are caller declarations, not conclusions from this solver.
//!
//! A window begins at a physical stress extremum, contains one monotone outward
//! and one monotone return segment, and is closed EXPLICITLY. Plateaus are
//! accepted; additional reversals refuse. Every physical extremum must be
//! resolved by the supplied load samples. Closing a window advances exactly one
//! full cycle at the currently accepted crack length and starts the next window
//! at the same actual endpoint. No rainflow completion order, repeated spectrum,
//! crack-closure correction or fatigue-threshold extrapolation is inferred.
//!
//! Mesh, material, support and thermal fields stay fixed. Calibration data may
//! be caller-declared, or accompanied by the original immutable material bundle;
//! the latter additionally binds chemistry, process and its `T` coordinate.
//! Both are empirical estimates and establish no validated component lifetime.

use std::sync::Arc;

use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_material::fatigue::NormalStressDirection;
use fs_material::fatigue::crack_growth::{
    ConstantGeometryCrack, CrackGrowthError, CrackGrowthLimits, CrackGrowthState, CrackGrowthStep,
    ParisCrackGrowth, ParisLaw, StressCycle,
};
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_math::det;

use super::{TetFatigueError, bind_material, stress_tensor};
use crate::linear3::{
    TetLinearElasticProblem, TetMaterialField, TetStaticSolveConfig, TetThermalStrainField,
    TetThermalStressReport,
};
use crate::loads3::transient::TransientError3;
use crate::loads3::{LoadError3, TetLoadCase3, TetLoadSolution3};

/// Numerical allowance, in Pa plus a fraction of the named stress scale.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CrackStressTolerance {
    pub absolute_stress_pa: f64,
    pub relative: f64,
}

impl CrackStressTolerance {
    fn validate(self) -> Result<(), TetCrackGrowthError> {
        if !self.absolute_stress_pa.is_finite()
            || self.absolute_stress_pa < 0.0
            || !self.relative.is_finite()
            || !(0.0..1.0).contains(&self.relative)
        {
            return Err(TetCrackGrowthError::Invalid(
                "stress allowances need finite absolute >= 0 and relative in [0, 1)",
            ));
        }
        Ok(())
    }

    fn allowance(self, scale: f64) -> Result<f64, TetCrackGrowthError> {
        let value = self.relative.mul_add(scale.abs(), self.absolute_stress_pa);
        if value.is_finite() {
            Ok(value)
        } else {
            Err(TetCrackGrowthError::Invalid("stress allowance overflow"))
        }
    }
}

/// One explicitly calibrated nominal-stress/crack-geometry observation.
#[derive(Debug, Clone)]
pub struct TetCrackGrowthGauge {
    pub element: usize,
    pub normal_in_mesh_frame: [f64; 3],
    pub material_state_identity: [u8; 32],
    /// Optional original bundle. A supplied bundle must be the actual native
    /// material identity and agree with the Paris scope's chemistry/process/T.
    pub material_source: Option<Arc<ResolvedMaterialStatePoint>>,
    pub law: Arc<ParisLaw>,
    pub geometry: ConstantGeometryCrack,
    pub initial_length_m: f64,
    pub limits: CrackGrowthLimits,
    /// Closure compares first/last normal stresses using the window peak scale.
    pub closure: CrackStressTolerance,
    /// Shear traction must not exceed this allowance at any accepted sample.
    /// The relative scale is the absolute normal traction, not von Mises stress.
    pub mode_i: CrackStressTolerance,
}

/// Shared bounds for quasistatic and inertial observations. Retained storage
/// is constant per gauge; no sampled trajectory or rainflow residue is kept.
#[derive(Debug, Clone, Copy)]
pub struct TetCrackGrowthConfig {
    pub maximum_absolute_mandel_strain: f64,
    pub maximum_gauges: usize,
    pub maximum_samples: u64,
}

#[derive(Debug)]
pub enum TetCrackGrowthError {
    Invalid(&'static str),
    Budget,
    Binding(TetFatigueError),
    Gauge {
        gauge: usize,
        error: CrackGrowthError,
    },
    Cycle {
        gauge: usize,
        what: &'static str,
    },
    MixedMode {
        gauge: usize,
        shear_pa: f64,
        allowance_pa: f64,
    },
    Load(LoadError3),
    Transient(TransientError3),
    Cancelled,
    Finished,
}

impl core::fmt::Display for TetCrackGrowthError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native crack-growth observation refused: {self:?}")
    }
}
impl std::error::Error for TetCrackGrowthError {}
impl From<TetFatigueError> for TetCrackGrowthError {
    fn from(error: TetFatigueError) -> Self {
        match error {
            TetFatigueError::Cancelled => Self::Cancelled,
            error => Self::Binding(error),
        }
    }
}
impl From<LoadError3> for TetCrackGrowthError {
    fn from(error: LoadError3) -> Self {
        match error {
            LoadError3::Cancelled => Self::Cancelled,
            error => Self::Load(error),
        }
    }
}
impl From<TransientError3> for TetCrackGrowthError {
    fn from(error: TransientError3) -> Self {
        match error {
            TransientError3::Cancelled => Self::Cancelled,
            error => Self::Transient(error),
        }
    }
}

/// One accepted native solve. Open-window samples do not grow the crack until
/// their actual full cycle closes. `completed_cycles` preserves gauge order.
#[derive(Debug)]
pub struct TetCrackGrowthSample {
    pub time_s: f64,
    pub mechanics: TetLoadSolution3,
    pub normal_stress_pa: Vec<f64>,
    pub cracks: Vec<CrackGrowthState>,
    pub completed_cycles: Option<Vec<CrackGrowthStep>>,
}

/// Original quasistatic solves plus their complete chronological crack state.
/// A failed solve, cycle closure or later gauge publishes no time/history/growth.
pub struct TetCrackGrowthPath<'p, 'm> {
    problem: &'p TetLinearElasticProblem<'m>,
    thermal: TetThermalStrainField<'p>,
    solver: TetStaticSolveConfig,
    observation: TetCrackGrowthObserver,
}

impl<'p, 'm> TetCrackGrowthPath<'p, 'm> {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        problem: &'p TetLinearElasticProblem<'m>,
        thermal: TetThermalStrainField<'p>,
        mesh_frame: ContentHash,
        gauges: Vec<TetCrackGrowthGauge>,
        solver: TetStaticSolveConfig,
        config: TetCrackGrowthConfig,
        cx: &Cx<'_>,
    ) -> Result<Self, TetCrackGrowthError> {
        if !solver.relative_residual_tolerance.is_finite()
            || !(0.0..1.0).contains(&solver.relative_residual_tolerance)
            || solver.relative_residual_tolerance == 0.0
            || solver.maximum_iterations == 0
        {
            return Err(TetCrackGrowthError::Invalid("invalid native solve bounds"));
        }
        let observation =
            TetCrackGrowthObserver::new(problem, thermal, mesh_frame, gauges, config, cx)?;
        Ok(Self {
            problem,
            thermal,
            solver,
            observation,
        })
    }

    pub fn states(&self) -> Vec<CrackGrowthState> {
        self.observation.states()
    }
    pub fn law(&self, gauge: usize) -> Option<&ParisLaw> {
        self.observation.law(gauge)
    }
    pub fn material_source(&self, gauge: usize) -> Option<&ResolvedMaterialStatePoint> {
        self.observation.material_source(gauge)
    }
    pub fn observed_samples(&self) -> u64 {
        self.observation.samples
    }
    pub fn last_time_s(&self) -> Option<f64> {
        self.observation.last_time_s
    }

    /// The first call begins at the first actual stress, not an invented zero.
    pub fn advance(
        &mut self,
        time_s: f64,
        load: &TetLoadCase3<'_, '_>,
        cx: &Cx<'_>,
    ) -> Result<TetCrackGrowthSample, TetCrackGrowthError> {
        self.advance_inner(time_s, load, false, cx)
    }

    /// Accept the mechanical endpoint and exactly one closed cycle together.
    pub fn advance_and_close_cycle(
        &mut self,
        time_s: f64,
        load: &TetLoadCase3<'_, '_>,
        cx: &Cx<'_>,
    ) -> Result<TetCrackGrowthSample, TetCrackGrowthError> {
        self.advance_inner(time_s, load, true, cx)
    }

    fn advance_inner(
        &mut self,
        time_s: f64,
        load: &TetLoadCase3<'_, '_>,
        close_cycle: bool,
        cx: &Cx<'_>,
    ) -> Result<TetCrackGrowthSample, TetCrackGrowthError> {
        self.observation.check_next(time_s, cx)?;
        let mechanics = load.solve(
            self.problem,
            self.thermal,
            self.solver,
            self.observation.mesh_frame(),
            self.observation.maximum_strain(),
            cx,
        )?;
        let (normal_stress_pa, completed_cycles) =
            self.observation
                .accept(time_s, &mechanics.stress, close_cycle, cx)?;
        Ok(TetCrackGrowthSample {
            time_s,
            mechanics,
            normal_stress_pa,
            cracks: self.states(),
            completed_cycles,
        })
    }

    /// Close the already accepted window without adding a duplicate sample.
    pub fn close_cycle(
        &mut self,
        cx: &Cx<'_>,
    ) -> Result<Vec<CrackGrowthStep>, TetCrackGrowthError> {
        self.observation.close_cycle(cx)
    }

    /// Require all sampled excursions to have been explicitly closed. A lone
    /// initial/last-closure endpoint contributes no invented partial cycle.
    pub fn finish(&mut self, cx: &Cx<'_>) -> Result<Vec<CrackGrowthState>, TetCrackGrowthError> {
        self.observation.finish(cx)
    }
}

#[derive(Debug, Clone)]
struct GaugeState {
    element: usize,
    normal: NormalStressDirection,
    material_source: Option<Arc<ResolvedMaterialStatePoint>>,
    growth: ParisCrackGrowth,
    closure: CrackStressTolerance,
    mode_i: CrackStressTolerance,
    window: Option<CycleWindow>,
}

/// Shared only by owners that recover their original native stress response.
#[derive(Debug, Clone)]
pub(crate) struct TetCrackGrowthObserver {
    mesh_frame: ContentHash,
    element_count: usize,
    config: TetCrackGrowthConfig,
    gauges: Vec<GaugeState>,
    samples: u64,
    last_time_s: Option<f64>,
    finished: bool,
}

impl TetCrackGrowthObserver {
    pub(crate) fn new(
        problem: &TetLinearElasticProblem<'_>,
        thermal: TetThermalStrainField<'_>,
        mesh_frame: ContentHash,
        gauges: Vec<TetCrackGrowthGauge>,
        config: TetCrackGrowthConfig,
        cx: &Cx<'_>,
    ) -> Result<Self, TetCrackGrowthError> {
        poll(cx)?;
        if gauges.is_empty() || gauges.len() > config.maximum_gauges || config.maximum_samples == 0
        {
            return Err(TetCrackGrowthError::Budget);
        }
        if mesh_frame == ContentHash([0; 32])
            || !config.maximum_absolute_mandel_strain.is_finite()
            || config.maximum_absolute_mandel_strain <= 0.0
            || config.maximum_absolute_mandel_strain > 0.05
        {
            return Err(TetCrackGrowthError::Invalid(
                "declare frame and small-strain bound in (0, 0.05]",
            ));
        }
        let element_count = problem.tetrahedra.len();
        if matches!(thermal, TetThermalStrainField::PerElement(t) if t.len() != element_count)
            || matches!(problem.materials, TetMaterialField::PerElement(m) if m.len() != element_count)
        {
            return Err(TetCrackGrowthError::Invalid(
                "material/thermal field differs from original connectivity",
            ));
        }
        let mut states = Vec::with_capacity(gauges.len());
        for (index, gauge) in gauges.into_iter().enumerate() {
            poll(cx)?;
            let scope = gauge.law.scope();
            bind_material(
                problem,
                thermal,
                index,
                gauge.element,
                gauge.material_state_identity,
                scope.temperature_k,
            )?;
            if let Some(source) = gauge.material_source.as_deref() {
                if source.identity().0 != gauge.material_state_identity
                    || source.material().chemistry != scope.material
                    || source.material().process != scope.process
                    || source.state_coordinate("T") != Some(scope.temperature_k)
                {
                    return Err(TetCrackGrowthError::Binding(
                        TetFatigueError::SourceMismatch {
                            gauge: index,
                            what: "Paris scope/bundle differs from actual elastic identity, material, process or T",
                        },
                    ));
                }
            }
            gauge.closure.validate()?;
            gauge.mode_i.validate()?;
            let normal = NormalStressDirection::new(gauge.normal_in_mesh_frame)
                .map_err(|_| TetCrackGrowthError::Invalid("invalid fixed crack-plane normal"))?;
            let growth = ParisCrackGrowth::new(
                gauge.law,
                gauge.geometry,
                gauge.initial_length_m,
                gauge.limits,
            )
            .map_err(|error| TetCrackGrowthError::Gauge {
                gauge: index,
                error,
            })?;
            states.push(GaugeState {
                element: gauge.element,
                normal,
                material_source: gauge.material_source,
                growth,
                closure: gauge.closure,
                mode_i: gauge.mode_i,
                window: None,
            });
        }
        poll(cx)?;
        Ok(Self {
            mesh_frame,
            element_count,
            config,
            gauges: states,
            samples: 0,
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
    pub(crate) fn observed_samples(&self) -> u64 {
        self.samples
    }
    pub(crate) fn states(&self) -> Vec<CrackGrowthState> {
        self.gauges.iter().map(|g| g.growth.state()).collect()
    }
    pub(crate) fn law(&self, gauge: usize) -> Option<&ParisLaw> {
        self.gauges.get(gauge).map(|g| g.growth.law())
    }
    pub(crate) fn material_source(&self, gauge: usize) -> Option<&ResolvedMaterialStatePoint> {
        self.gauges
            .get(gauge)
            .and_then(|g| g.material_source.as_deref())
    }

    pub(crate) fn check_next(&self, time_s: f64, cx: &Cx<'_>) -> Result<(), TetCrackGrowthError> {
        poll(cx)?;
        if self.finished {
            return Err(TetCrackGrowthError::Finished);
        }
        if !time_s.is_finite()
            || time_s < 0.0
            || self.last_time_s.is_some_and(|last| time_s <= last)
        {
            return Err(TetCrackGrowthError::Invalid(
                "sample times must be finite, nonnegative and increasing",
            ));
        }
        if self.samples >= self.config.maximum_samples {
            return Err(TetCrackGrowthError::Budget);
        }
        Ok(())
    }

    pub(crate) fn accept(
        &mut self,
        time_s: f64,
        stress: &TetThermalStressReport,
        close_cycle: bool,
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, Option<Vec<CrackGrowthStep>>), TetCrackGrowthError> {
        self.check_next(time_s, cx)?;
        if stress.mesh_frame() != self.mesh_frame || stress.elements().len() != self.element_count {
            return Err(TetCrackGrowthError::Invalid(
                "stress differs from original frame/connectivity",
            ));
        }
        let mut candidate = self.gauges.clone();
        let mut values = Vec::with_capacity(candidate.len());
        for (index, gauge) in candidate.iter_mut().enumerate() {
            poll(cx)?;
            let tensor = stress_tensor(stress, gauge.element, self.mesh_frame)?;
            let normal_stress = gauge
                .normal
                .project(&tensor)
                .map_err(|_| TetCrackGrowthError::Invalid("nonfinite native normal stress"))?;
            let n = gauge.normal.unit_normal();
            let traction = [
                tensor[0] * n[0] + tensor[3] * n[1] + tensor[5] * n[2],
                tensor[3] * n[0] + tensor[1] * n[1] + tensor[4] * n[2],
                tensor[5] * n[0] + tensor[4] * n[1] + tensor[2] * n[2],
            ];
            let shear = det::hypot(
                det::hypot(
                    traction[0] - normal_stress * n[0],
                    traction[1] - normal_stress * n[1],
                ),
                traction[2] - normal_stress * n[2],
            );
            let allowance = gauge.mode_i.allowance(normal_stress)?;
            if !shear.is_finite() || shear > allowance {
                return Err(TetCrackGrowthError::MixedMode {
                    gauge: index,
                    shear_pa: shear,
                    allowance_pa: allowance,
                });
            }
            gauge
                .growth
                .stress_intensity(normal_stress)
                .map_err(|error| TetCrackGrowthError::Gauge {
                    gauge: index,
                    error,
                })?;
            if let Some(window) = gauge.window.as_mut() {
                window.push(normal_stress, gauge.closure, index)?;
            } else {
                gauge.window = Some(CycleWindow::new(normal_stress));
            }
            values.push(normal_stress);
        }
        let completed = if close_cycle {
            Some(close_gauges(&mut candidate, cx)?)
        } else {
            None
        };
        poll(cx)?;
        self.gauges = candidate;
        self.samples += 1;
        self.last_time_s = Some(time_s);
        Ok((values, completed))
    }

    pub(crate) fn close_cycle(
        &mut self,
        cx: &Cx<'_>,
    ) -> Result<Vec<CrackGrowthStep>, TetCrackGrowthError> {
        poll(cx)?;
        if self.finished {
            return Err(TetCrackGrowthError::Finished);
        }
        let mut candidate = self.gauges.clone();
        let steps = close_gauges(&mut candidate, cx)?;
        poll(cx)?;
        self.gauges = candidate;
        Ok(steps)
    }

    pub(crate) fn finish(
        &mut self,
        cx: &Cx<'_>,
    ) -> Result<Vec<CrackGrowthState>, TetCrackGrowthError> {
        poll(cx)?;
        for (index, gauge) in self.gauges.iter().enumerate() {
            if gauge.window.is_some_and(|w| w.samples > 1) {
                return Err(TetCrackGrowthError::Cycle {
                    gauge: index,
                    what: "close the complete final cycle; partial cycles cannot be inferred",
                });
            }
        }
        poll(cx)?;
        self.finished = true;
        Ok(self.states())
    }
}

#[derive(Debug, Clone, Copy)]
struct CycleWindow {
    first: f64,
    last: f64,
    minimum: f64,
    maximum: f64,
    direction: i8,
    reversed: bool,
    samples: u64,
}

impl CycleWindow {
    fn new(value: f64) -> Self {
        Self {
            first: value,
            last: value,
            minimum: value,
            maximum: value,
            direction: 0,
            reversed: false,
            samples: 1,
        }
    }

    fn push(
        &mut self,
        value: f64,
        tolerance: CrackStressTolerance,
        gauge: usize,
    ) -> Result<(), TetCrackGrowthError> {
        let direction = if value > self.last {
            1
        } else if value < self.last {
            -1
        } else {
            0
        };
        if direction != 0 {
            if self.direction == 0 {
                self.direction = direction;
            } else if direction != self.direction {
                self.reversed = true;
            } else if self.reversed {
                return Err(TetCrackGrowthError::Cycle {
                    gauge,
                    what: "additional reversal: close each chronological simple cycle explicitly",
                });
            }
        }
        let minimum = self.minimum.min(value);
        let maximum = self.maximum.max(value);
        let allowance = tolerance.allowance(maximum.abs().max(minimum.abs()))?;
        if self.reversed && f64::from(self.direction) * (value - self.first) < -allowance {
            return Err(TetCrackGrowthError::Cycle {
                gauge,
                what: "return segment crossed the starting extremum without closing the cycle",
            });
        }
        self.minimum = minimum;
        self.maximum = maximum;
        self.last = value;
        self.samples += 1;
        Ok(())
    }

    fn cycle(
        self,
        tolerance: CrackStressTolerance,
        gauge: usize,
    ) -> Result<StressCycle, TetCrackGrowthError> {
        if self.samples < 2 || (self.direction != 0 && !self.reversed) {
            return Err(TetCrackGrowthError::Cycle {
                gauge,
                what: "a complete outward/return excursion is required",
            });
        }
        let allowance = tolerance.allowance(self.maximum.abs().max(self.minimum.abs()))?;
        if (self.last - self.first).abs() > allowance {
            return Err(TetCrackGrowthError::Cycle {
                gauge,
                what: "normal stress has not returned to the original extremum",
            });
        }
        Ok(StressCycle {
            minimum_pa: self.minimum,
            maximum_pa: self.maximum,
            cycles: 1.0,
        })
    }
}

fn close_gauges(
    gauges: &mut [GaugeState],
    cx: &Cx<'_>,
) -> Result<Vec<CrackGrowthStep>, TetCrackGrowthError> {
    let mut steps = Vec::with_capacity(gauges.len());
    for (index, gauge) in gauges.iter_mut().enumerate() {
        poll(cx)?;
        let window = gauge.window.ok_or(TetCrackGrowthError::Cycle {
            gauge: index,
            what: "no native stress has been observed",
        })?;
        let cycle = window.cycle(gauge.closure, index)?;
        let step = gauge
            .growth
            .advance_with_cancel(cycle, &mut || cx.checkpoint().is_err())
            .map_err(|error| TetCrackGrowthError::Gauge {
                gauge: index,
                error,
            })?;
        gauge.window = Some(CycleWindow::new(window.last));
        steps.push(step);
    }
    Ok(steps)
}

fn poll(cx: &Cx<'_>) -> Result<(), TetCrackGrowthError> {
    cx.checkpoint().map_err(|_| TetCrackGrowthError::Cancelled)
}
