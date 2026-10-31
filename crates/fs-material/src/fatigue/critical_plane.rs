//! Tensor Smith-Watson-Topper (SWT) on one resolved, closed periodic cycle.
//!
//! Search a caller-declared finite set of FIXED planes for maximum normal
//! mechanical-strain range, then use the maximum normal stress on that plane
//! over the entire observed period. Exact strain-range ties select maximum
//! normal stress, then the earliest input plane. Maximizing the SWT product
//! instead is a different plane-selection convention under phase-lagged loading.
//! See COMSOL "Strain Based Fatigue Models" and NASA CR-182126, section 3.2.3:
//! <https://doc.comsol.com/6.4/doc/com.comsol.help.fatigue/fatigue_ug_sme.4.23.html>
//! <https://ntrs.nasa.gov/citations/19880015632>
//!
//! Stress extrema need not coincide with strain extrema. Every nonconstant
//! projected strain history must nevertheless contain exactly one maximum and
//! minimum around the observed period. Multi-cycle blocks and arbitrary
//! variable-amplitude multiaxial counting are outside this owner. This does not
//! predict shear-mode failure, nonproportional cyclic hardening or stabilization.
//! A stabilized, repeatable response and transferable tensile-mode calibration
//! remain applicability requirements; endpoint agreement alone proves neither.
//!
//! All tensors and planes share one fixed frame and tensor-Voigt order
//! [xx, yy, zz, xy, yz, zx]. Strain is total MECHANICAL strain, excluding thermal
//! eigenstrain, not engineering shear or equivalent accumulated plastic strain.
//! Observations must resolve physical extrema. No endpoint or between-sample
//! peak is fabricated. Closure tolerance explicitly admits endpoint mismatch
//! without changing any sample. The result is a finite-plane model estimate,
//! not a certified optimum over the sphere or a component-life qualification.

pub mod findley;

/// Shear-mode Fatemi-Socie life on a declared closed mechanical cycle.
pub mod fatemi_socie;

/// Variable-amplitude finite-plane strain-life accumulation.
pub mod miner;

mod shear_path;

use super::strain_life::{StrainLifeCurve, StrainLifeMeanStress};
use super::{FatigueError, NormalStressDirection, half_range};
use fs_math::det;
use std::sync::Arc;

/// Simultaneous signed tensors in the same fixed frame.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TensorFatigueSample {
    pub stress_pa: [f64; 6],
    pub total_mechanical_strain: [f64; 6],
}

/// A work unit projects one paired stress/strain observation onto one plane.
/// Storage is O(planes), including a preallocated transactional scratch buffer.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CriticalPlaneLimits {
    pub maximum_samples: u64,
    pub maximum_planes: usize,
    pub maximum_projection_work: u64,
}

/// Componentwise endpoint agreement: |last-first| <= absolute + relative *
/// max(|first|, |last|). These are caller limits, not error certificates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CycleClosureTolerance {
    pub absolute_strain: f64,
    pub absolute_stress_pa: f64,
    pub relative: f64,
}

/// A single-period, tensile-mode estimate on the declared finite plane set.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SwtCycleResult {
    pub samples: u64,
    pub evaluated_planes: usize,
    pub selected_plane: usize,
    pub normal_in_frame: [f64; 3],
    pub strain_amplitude: f64,
    pub minimum_normal_stress_pa: f64,
    pub maximum_normal_stress_pa: f64,
    pub cycles_to_failure: Option<f64>,
    /// One divided by predicted full cycles; zero only for zero strain range.
    pub damage_per_cycle: f64,
    pub closure: CycleClosureTolerance,
    pub maximum_strain_closure_defect: f64,
    pub maximum_stress_closure_defect_pa: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum CriticalPlaneError {
    Invalid(&'static str),
    Budget,
    InvalidSample {
        sample: u64,
        component: usize,
        quantity: &'static str,
    },
    NotClosed {
        component: usize,
        quantity: &'static str,
    },
    NotSingleCycle {
        plane: usize,
    },
    Fatigue(FatigueError),
    Allocation,
    Cancelled,
    Finished,
}
impl core::fmt::Display for CriticalPlaneError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "critical-plane SWT cycle refused: {self:?}")
    }
}
impl std::error::Error for CriticalPlaneError {}
impl From<FatigueError> for CriticalPlaneError {
    fn from(error: FatigueError) -> Self {
        match error {
            FatigueError::Cancelled => Self::Cancelled,
            FatigueError::Allocation => Self::Allocation,
            error => Self::Fatigue(error),
        }
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq)]
struct PlaneState {
    minimum_strain: f64,
    maximum_strain: f64,
    minimum_stress: f64,
    maximum_stress: f64,
    last_strain: f64,
    first_direction: i8,
    last_direction: i8,
    // Three means already invalid. Keep warm-up bounded and permit restart.
    direction_changes: u8,
}
impl PlaneState {
    fn observe(mut self, strain: f64, stress: f64, first: bool) -> Self {
        if first {
            return Self {
                minimum_strain: strain,
                maximum_strain: strain,
                minimum_stress: stress,
                maximum_stress: stress,
                last_strain: strain,
                ..Self::default()
            };
        }
        let direction = if strain > self.last_strain {
            1
        } else if strain < self.last_strain {
            -1
        } else {
            0
        };
        if direction != 0 {
            if self.first_direction == 0 {
                self.first_direction = direction;
            }
            if self.last_direction != 0 && direction != self.last_direction {
                self.direction_changes = (self.direction_changes + 1).min(3);
            }
            self.last_direction = direction;
        }
        self.last_strain = strain;
        self.minimum_strain = self.minimum_strain.min(strain);
        self.maximum_strain = self.maximum_strain.max(strain);
        self.minimum_stress = self.minimum_stress.min(stress);
        self.maximum_stress = self.maximum_stress.max(stress);
        self
    }

    fn is_single_cycle(self) -> bool {
        if self.minimum_strain == self.maximum_strain {
            return true;
        }
        let seam = u8::from(self.first_direction != self.last_direction);
        self.direction_changes + seam == 2
    }
}

#[derive(Debug, Clone, Copy)]
struct SearchPlane {
    normal: [f64; 3],
    weights: [f64; 6],
}
impl SearchPlane {
    fn new(direction: [f64; 3]) -> Result<Self, FatigueError> {
        let normal = NormalStressDirection::new(direction)?.unit_normal();
        let [x, y, z] = normal;
        Ok(Self {
            normal,
            weights: [x * x, y * y, z * z, 2.0 * x * y, 2.0 * y * z, 2.0 * z * x],
        })
    }

    fn project(self, sample: TensorFatigueSample) -> Result<(f64, f64), FatigueError> {
        let mut strain = 0.0;
        let mut stress = 0.0;
        for (component, weight) in self.weights.into_iter().enumerate() {
            strain += weight * sample.total_mechanical_strain[component];
            stress += weight * sample.stress_pa[component];
        }
        if !strain.is_finite() || !stress.is_finite() {
            return Err(FatigueError::Arithmetic);
        }
        Ok((strain, stress))
    }
}

/// Streaming O(planes) state. Pushes perform no allocation; a block proposal
/// uses one additional bounded snapshot. Refusal/cancellation never publishes
/// partial observations or results. The scratch buffer is not accepted history.
#[derive(Debug, Clone)]
pub struct SwtCycle {
    curve: Arc<StrainLifeCurve>,
    planes: Arc<Vec<SearchPlane>>,
    limits: CriticalPlaneLimits,
    closure: CycleClosureTolerance,
    states: Vec<PlaneState>,
    scratch: Vec<PlaneState>,
    first: Option<TensorFatigueSample>,
    last: Option<TensorFatigueSample>,
    samples: u64,
    projection_work: u64,
    result: Option<SwtCycleResult>,
}
impl SwtCycle {
    pub fn new(
        curve: Arc<StrainLifeCurve>,
        normals: Vec<[f64; 3]>,
        limits: CriticalPlaneLimits,
        closure: CycleClosureTolerance,
    ) -> Result<Self, CriticalPlaneError> {
        if limits.maximum_samples < 3
            || limits.maximum_planes == 0
            || limits.maximum_projection_work == 0
        {
            return Err(CriticalPlaneError::Invalid(
                "positive limits and at least three sample slots are required",
            ));
        }
        if normals.is_empty() {
            return Err(CriticalPlaneError::Invalid(
                "declare at least one search plane",
            ));
        }
        if normals.len() > limits.maximum_planes {
            return Err(CriticalPlaneError::Budget);
        }
        let count = u64::try_from(normals.len()).map_err(|_| CriticalPlaneError::Budget)?;
        if count > limits.maximum_projection_work {
            return Err(CriticalPlaneError::Budget);
        }
        if [
            closure.absolute_strain,
            closure.absolute_stress_pa,
            closure.relative,
        ]
        .into_iter()
        .any(|x| !x.is_finite() || x < 0.0)
        {
            return Err(CriticalPlaneError::Invalid(
                "closure tolerances must be finite and nonnegative",
            ));
        }
        // Reuse policy admission; a plastic-only CMB fit cannot supply SWT.
        curve.cycles_to_failure_with_stress(
            0.0,
            0.0,
            0.0,
            StrainLifeMeanStress::SmithWatsonTopper,
        )?;
        let mut planes = Vec::new();
        planes
            .try_reserve_exact(normals.len())
            .map_err(|_| CriticalPlaneError::Allocation)?;
        for normal in normals {
            planes.push(SearchPlane::new(normal)?);
        }
        let states = empty_states(planes.len())?;
        let scratch = empty_states(planes.len())?;
        Ok(Self {
            curve,
            planes: Arc::new(planes),
            limits,
            closure,
            states,
            scratch,
            first: None,
            last: None,
            samples: 0,
            projection_work: 0,
            result: None,
        })
    }

    pub fn curve(&self) -> &StrainLifeCurve {
        &self.curve
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
    pub fn limits(&self) -> CriticalPlaneLimits {
        self.limits
    }
    pub fn closure_tolerance(&self) -> CycleClosureTolerance {
        self.closure
    }
    pub fn result(&self) -> Option<&SwtCycleResult> {
        self.result.as_ref()
    }

    /// Start an empty observation window with the SAME immutable calibration,
    /// search planes, closure policy and limits. Does not mutate this window
    /// or silently copy its endpoint; the consumer supplies the new first sample.
    pub fn fresh(&self) -> Result<Self, CriticalPlaneError> {
        Ok(Self {
            curve: self.curve.clone(),
            planes: self.planes.clone(),
            limits: self.limits,
            closure: self.closure,
            states: empty_states(self.planes.len())?,
            scratch: empty_states(self.planes.len())?,
            first: None,
            last: None,
            samples: 0,
            projection_work: 0,
            result: None,
        })
    }

    pub fn push(&mut self, sample: TensorFatigueSample) -> Result<(), CriticalPlaneError> {
        self.push_with_cancel(sample, &mut || false)
    }

    pub fn push_with_cancel(
        &mut self,
        sample: TensorFatigueSample,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), CriticalPlaneError> {
        if cancelled() {
            return Err(CriticalPlaneError::Cancelled);
        }
        if self.result.is_some() {
            return Err(CriticalPlaneError::Finished);
        }
        let (next_samples, next_work) = self.admit_count(1)?;
        for (quantity, values) in [
            ("stress", &sample.stress_pa),
            ("mechanical strain", &sample.total_mechanical_strain),
        ] {
            for (component, value) in values.iter().enumerate() {
                if !value.is_finite() {
                    return Err(CriticalPlaneError::InvalidSample {
                        sample: self.samples,
                        component,
                        quantity,
                    });
                }
            }
        }
        for (index, normal) in self.planes.iter().enumerate() {
            if cancelled() {
                return Err(CriticalPlaneError::Cancelled);
            }
            let (strain, stress) = normal.project(sample)?;
            self.scratch[index] = self.states[index].observe(strain, stress, self.samples == 0);
        }
        if cancelled() {
            return Err(CriticalPlaneError::Cancelled);
        }
        std::mem::swap(&mut self.states, &mut self.scratch);
        if self.first.is_none() {
            self.first = Some(sample);
        }
        self.last = Some(sample);
        self.samples = next_samples;
        self.projection_work = next_work;
        Ok(())
    }

    pub fn push_block(
        &mut self,
        samples: &[TensorFatigueSample],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), CriticalPlaneError> {
        if cancelled() {
            return Err(CriticalPlaneError::Cancelled);
        }
        if self.result.is_some() {
            return Err(CriticalPlaneError::Finished);
        }
        self.admit_count(u64::try_from(samples.len()).map_err(|_| CriticalPlaneError::Budget)?)?;
        let mut candidate = self.fresh()?;
        candidate.states.copy_from_slice(&self.states);
        candidate.first = self.first;
        candidate.last = self.last;
        candidate.samples = self.samples;
        candidate.projection_work = self.projection_work;
        for &sample in samples {
            candidate.push_with_cancel(sample, cancelled)?;
        }
        if cancelled() {
            return Err(CriticalPlaneError::Cancelled);
        }
        *self = candidate;
        Ok(())
    }

    pub fn finish(&mut self) -> Result<SwtCycleResult, CriticalPlaneError> {
        self.finish_with_cancel(&mut || false)
    }

    pub fn finish_with_cancel(
        &mut self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SwtCycleResult, CriticalPlaneError> {
        if cancelled() {
            return Err(CriticalPlaneError::Cancelled);
        }
        if let Some(result) = self.result {
            return Ok(result);
        }
        if self.samples < 3 {
            return Err(CriticalPlaneError::Invalid(
                "a closed cycle requires at least three observations including both endpoints",
            ));
        }
        let first = self
            .first
            .ok_or(CriticalPlaneError::Invalid("missing initial observation"))?;
        let last = self
            .last
            .ok_or(CriticalPlaneError::Invalid("missing final observation"))?;
        let maximum_strain_closure_defect = closure_defect(
            &first.total_mechanical_strain,
            &last.total_mechanical_strain,
            self.closure.absolute_strain,
            self.closure.relative,
            "mechanical strain",
        )?;
        let maximum_stress_closure_defect_pa = closure_defect(
            &first.stress_pa,
            &last.stress_pa,
            self.closure.absolute_stress_pa,
            self.closure.relative,
            "stress",
        )?;
        let mut selected = 0;
        let mut amplitude = 0.0;
        for (index, state) in self.states.iter().copied().enumerate() {
            if cancelled() {
                return Err(CriticalPlaneError::Cancelled);
            }
            if !state.is_single_cycle() {
                return Err(CriticalPlaneError::NotSingleCycle { plane: index });
            }
            let candidate = half_range(state.minimum_strain, state.maximum_strain)?;
            if index == 0
                || candidate > amplitude
                || (candidate == amplitude
                    && state.maximum_stress > self.states[selected].maximum_stress)
            {
                selected = index;
                amplitude = candidate;
            }
        }
        let state = self.states[selected];
        let life = self.curve.cycles_to_failure_with_stress(
            amplitude,
            state.minimum_stress,
            state.maximum_stress,
            StrainLifeMeanStress::SmithWatsonTopper,
        )?;
        let damage = life.map_or(0.0, |n| 1.0 / n);
        if !damage.is_finite() {
            return Err(FatigueError::Arithmetic.into());
        }
        let result = SwtCycleResult {
            samples: self.samples,
            evaluated_planes: self.planes.len(),
            selected_plane: selected,
            normal_in_frame: self.planes[selected].normal,
            strain_amplitude: amplitude,
            minimum_normal_stress_pa: state.minimum_stress,
            maximum_normal_stress_pa: state.maximum_stress,
            cycles_to_failure: life,
            damage_per_cycle: damage,
            closure: self.closure,
            maximum_strain_closure_defect,
            maximum_stress_closure_defect_pa,
        };
        if cancelled() {
            return Err(CriticalPlaneError::Cancelled);
        }
        self.result = Some(result);
        Ok(result)
    }

    fn admit_count(&self, count: u64) -> Result<(u64, u64), CriticalPlaneError> {
        let next_samples = self
            .samples
            .checked_add(count)
            .ok_or(CriticalPlaneError::Budget)?;
        let planes = u64::try_from(self.planes.len()).map_err(|_| CriticalPlaneError::Budget)?;
        let next_work = count
            .checked_mul(planes)
            .and_then(|n| self.projection_work.checked_add(n))
            .ok_or(CriticalPlaneError::Budget)?;
        if next_samples > self.limits.maximum_samples
            || next_work > self.limits.maximum_projection_work
        {
            Err(CriticalPlaneError::Budget)
        } else {
            Ok((next_samples, next_work))
        }
    }
}

fn empty_states(count: usize) -> Result<Vec<PlaneState>, CriticalPlaneError> {
    let mut states = Vec::new();
    states
        .try_reserve_exact(count)
        .map_err(|_| CriticalPlaneError::Allocation)?;
    states.resize(count, PlaneState::default());
    Ok(states)
}

fn closure_defect(
    first: &[f64; 6],
    last: &[f64; 6],
    absolute: f64,
    relative: f64,
    quantity: &'static str,
) -> Result<f64, CriticalPlaneError> {
    let mut maximum = 0.0_f64;
    for component in 0..6 {
        let a = first[component];
        let b = last[component];
        if a == b {
            continue;
        }
        let difference = (a - b).abs();
        if !difference.is_finite() {
            return Err(FatigueError::Arithmetic.into());
        }
        let scale = a.abs().max(b.abs()).max(absolute);
        // Scaling keeps an otherwise overflowing tolerance sum harmless.
        let threshold = absolute / scale + relative * (a.abs().max(b.abs()) / scale);
        if (absolute == 0.0 && relative == 0.0) || difference / scale > threshold {
            return Err(CriticalPlaneError::NotClosed {
                component,
                quantity,
            });
        }
        maximum = maximum.max(difference);
    }
    Ok(maximum)
}

/// Nested polar/azimuth search on the unoriented hemisphere (n and -n are
/// equivalent planes). Resolution one is exactly the x, y, z axes, ordered
/// z, x, y. Doubling resolution retains the geometric coarse directions.
/// There are 4*q*q - 2*q + 1 normals. This is a finite search grid, not an
/// angular error certificate or proof of a continuum critical plane.
pub fn hemisphere_planes(
    resolution: usize,
    maximum_planes: usize,
) -> Result<Vec<[f64; 3]>, CriticalPlaneError> {
    if resolution == 0 {
        return Err(CriticalPlaneError::Invalid(
            "positive hemisphere resolution required",
        ));
    }
    let azimuths = resolution
        .checked_mul(4)
        .ok_or(CriticalPlaneError::Budget)?;
    let equator = resolution
        .checked_mul(2)
        .ok_or(CriticalPlaneError::Budget)?;
    let count = resolution
        .checked_sub(1)
        .and_then(|n| n.checked_mul(azimuths))
        .and_then(|n| n.checked_add(equator))
        .and_then(|n| n.checked_add(1))
        .ok_or(CriticalPlaneError::Budget)?;
    if count > maximum_planes {
        return Err(CriticalPlaneError::Budget);
    }
    let mut planes = Vec::new();
    planes
        .try_reserve_exact(count)
        .map_err(|_| CriticalPlaneError::Allocation)?;
    planes.push([0.0, 0.0, 1.0]);
    for band in 1..resolution {
        let theta = core::f64::consts::FRAC_PI_2 * (band as f64 / resolution as f64);
        let (sin_theta, cos_theta) = det::sin_cos(theta);
        for azimuth in 0..azimuths {
            let (sin_phi, cos_phi) = azimuth_direction(azimuth, resolution);
            planes.push([sin_theta * cos_phi, sin_theta * sin_phi, cos_theta]);
        }
    }
    for azimuth in 0..equator {
        let (sin_phi, cos_phi) = azimuth_direction(azimuth, resolution);
        planes.push([cos_phi, sin_phi, 0.0]);
    }
    Ok(planes)
}

fn azimuth_direction(index: usize, resolution: usize) -> (f64, f64) {
    if index.is_multiple_of(resolution) {
        return match index / resolution {
            0 => (0.0, 1.0),
            1 => (1.0, 0.0),
            2 => (0.0, -1.0),
            _ => (-1.0, 0.0),
        };
    }
    det::sin_cos(core::f64::consts::FRAC_PI_2 * (index as f64 / resolution as f64))
}

#[cfg(test)]
mod tests {
    use super::super::FatigueScope;
    use super::super::strain_life::CmbParameters;
    use super::*;

    fn curve() -> Arc<StrainLifeCurve> {
        Arc::new(
            StrainLifeCurve::new(
                FatigueScope {
                    source_reference: "synthetic analytical SWT references; not measured life"
                        .into(),
                    material: "fixture".into(),
                    process: "fixture".into(),
                    surface_and_geometry: "smooth tensile-mode coupon".into(),
                    environment: "fixed temperature".into(),
                    load_basis: "closed tensor mechanical-strain and stress cycle".into(),
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
    fn exact_closure() -> CycleClosureTolerance {
        CycleClosureTolerance {
            absolute_strain: 0.0,
            absolute_stress_pa: 0.0,
            relative: 0.0,
        }
    }
    fn limits() -> CriticalPlaneLimits {
        CriticalPlaneLimits {
            maximum_samples: 1000,
            maximum_planes: 2048,
            maximum_projection_work: 100000,
        }
    }
    fn cycle(planes: Vec<[f64; 3]>) -> SwtCycle {
        SwtCycle::new(curve(), planes, limits(), exact_closure()).unwrap()
    }
    fn axial(strain: f64, stress: f64) -> TensorFatigueSample {
        TensorFatigueSample {
            total_mechanical_strain: [strain, 0.0, 0.0, 0.0, 0.0, 0.0],
            stress_pa: [stress, 0.0, 0.0, 0.0, 0.0, 0.0],
        }
    }
    fn axial_cycle() -> [TensorFatigueSample; 3] {
        [
            axial(-0.001, -200e6),
            axial(0.001, 200e6),
            axial(-0.001, -200e6),
        ]
    }
    fn close(a: f64, b: f64) {
        assert!((a - b).abs() <= 8e-12 * b.abs().max(1e-20), "{a} != {b}");
    }
    fn same_accepted(a: &SwtCycle, b: &SwtCycle) {
        assert_eq!(a.states, b.states);
        assert_eq!(a.first, b.first);
        assert_eq!(a.last, b.last);
        assert_eq!(a.samples, b.samples);
        assert_eq!(a.projection_work, b.projection_work);
        assert_eq!(a.result, b.result);
    }
    // For b=c=-1/2, sigma_max*epsilon_a = 7e6/(2*N).
    fn analytical_life(strain_amplitude: f64, maximum_stress: f64) -> f64 {
        7e6 / (2.0 * strain_amplitude * maximum_stress)
    }

    #[test]
    fn g1_axial_reduction_uses_full_cycles_and_only_inverts_the_selected_plane() {
        let mut flow = cycle(vec![[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]);
        let mut samples = axial_cycle();
        for sample in &mut samples {
            // This irrelevant plane is below the finite calibration's support.
            sample.total_mechanical_strain[1] = sample.total_mechanical_strain[0] * 1e-20;
        }
        flow.push_block(&samples, &mut || false).unwrap();
        let result = flow.finish().unwrap();
        assert_eq!(result.selected_plane, 1);
        assert_eq!(result.normal_in_frame, [1.0, 0.0, 0.0]);
        assert_eq!(result.strain_amplitude, 0.001);
        close(
            result.cycles_to_failure.unwrap(),
            analytical_life(0.001, 200e6),
        );
        close(result.damage_per_cycle, 1.0 / analytical_life(0.001, 200e6));
        assert_eq!((result.samples, result.evaluated_planes), (3, 3));
        assert_eq!(flow.projection_work(), 9);
        assert_eq!(flow.result(), Some(&result));
        assert_eq!(flow.finish().unwrap(), result);
        assert_eq!(flow.push(samples[0]), Err(CriticalPlaneError::Finished));
    }

    #[test]
    fn g1_lagged_stress_peak_is_observed_away_from_both_strain_extrema() {
        let samples = [
            axial(-0.001, 0.0),
            axial(0.0, 300e6),
            axial(0.001, 0.0),
            axial(0.0, -300e6),
            axial(-0.001, 0.0),
        ];
        let mut flow = cycle(vec![[1.0, 0.0, 0.0]]);
        flow.push_block(&samples, &mut || false).unwrap();
        let result = flow.finish().unwrap();
        assert_eq!(result.maximum_normal_stress_pa, 300e6);
        close(
            result.cycles_to_failure.unwrap(),
            analytical_life(0.001, 300e6),
        );
        for split in 0..=samples.len() {
            let mut chunked = flow.fresh().unwrap();
            chunked
                .push_block(&samples[..split], &mut || false)
                .unwrap();
            chunked
                .push_block(&samples[split..], &mut || false)
                .unwrap();
            assert_eq!(chunked.finish().unwrap(), result);
        }
    }

    #[test]
    fn g1_plane_selection_is_maximum_strain_range_with_explicit_ties_not_maximum_product() {
        let mut flow = cycle(vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
        for sign in [-1.0, 1.0, -1.0] {
            flow.push(TensorFatigueSample {
                total_mechanical_strain: [0.002 * sign, 0.001 * sign, 0.0, 0.0, 0.0, 0.0],
                stress_pa: [100e6 * sign, 500e6 * sign, 0.0, 0.0, 0.0, 0.0],
            })
            .unwrap();
        }
        let result = flow.finish().unwrap();
        assert_eq!(result.selected_plane, 0);
        close(
            result.cycles_to_failure.unwrap(),
            analytical_life(0.002, 100e6),
        );
        assert!(0.001 * 500e6 > result.strain_amplitude * result.maximum_normal_stress_pa);

        let mut tied = cycle(vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
        for sign in [-1.0, 1.0, -1.0] {
            tied.push(TensorFatigueSample {
                total_mechanical_strain: [0.001 * sign, 0.001 * sign, 0.001 * sign, 0.0, 0.0, 0.0],
                stress_pa: [100e6 * sign, 300e6 * sign, 300e6 * sign, 0.0, 0.0, 0.0],
            })
            .unwrap();
        }
        assert_eq!(tied.finish().unwrap().selected_plane, 1);
    }

    #[test]
    fn g1_tensor_shear_projects_to_the_tensile_plane_without_engineering_factor_error() {
        let mut flow = cycle(vec![[1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]);
        for sign in [-1.0, 1.0, -1.0] {
            flow.push(TensorFatigueSample {
                total_mechanical_strain: [0.0, 0.0, 0.0, 0.001 * sign, 0.0, 0.0],
                stress_pa: [0.0, 0.0, 0.0, 100e6 * sign, 0.0, 0.0],
            })
            .unwrap();
        }
        let result = flow.finish().unwrap();
        assert_eq!(result.selected_plane, 1);
        close(result.strain_amplitude, 0.001);
        close(result.maximum_normal_stress_pa, 100e6);
        close(
            result.cycles_to_failure.unwrap(),
            analytical_life(0.001, 100e6),
        );
    }

    fn rotate_vector([x, y, z]: [f64; 3]) -> [f64; 3] {
        [0.6 * x - 0.8 * y, 0.8 * x + 0.6 * y, z]
    }
    fn rotate_tensor(a: [f64; 6]) -> [f64; 6] {
        let c = 0.6;
        let s = 0.8;
        [
            c * c * a[0] + s * s * a[1] - 2.0 * c * s * a[3],
            s * s * a[0] + c * c * a[1] + 2.0 * c * s * a[3],
            a[2],
            c * s * (a[0] - a[1]) + (c * c - s * s) * a[3],
            s * a[5] + c * a[4],
            c * a[5] - s * a[4],
        ]
    }

    #[test]
    fn g3_rotating_tensors_and_planes_together_preserves_life_and_plane_orientation() {
        let planes = vec![
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
        ];
        let mut original = cycle(planes.clone());
        let mut rotated = cycle(planes.into_iter().map(rotate_vector).collect());
        for sample in axial_cycle() {
            original.push(sample).unwrap();
            rotated
                .push(TensorFatigueSample {
                    stress_pa: rotate_tensor(sample.stress_pa),
                    total_mechanical_strain: rotate_tensor(sample.total_mechanical_strain),
                })
                .unwrap();
        }
        let a = original.finish().unwrap();
        let b = rotated.finish().unwrap();
        close(a.strain_amplitude, b.strain_amplitude);
        close(a.maximum_normal_stress_pa, b.maximum_normal_stress_pa);
        close(a.cycles_to_failure.unwrap(), b.cycles_to_failure.unwrap());
        for (actual, expected) in b
            .normal_in_frame
            .into_iter()
            .zip(rotate_vector(a.normal_in_frame))
        {
            close(actual, expected);
        }
        let mut opposite = cycle(vec![[-1.0, 0.0, 0.0]]);
        opposite.push_block(&axial_cycle(), &mut || false).unwrap();
        close(
            opposite.finish().unwrap().cycles_to_failure.unwrap(),
            a.cycles_to_failure.unwrap(),
        );
    }

    #[test]
    fn g3_cyclic_seam_origin_and_plateaus_work_but_nested_periods_require_restart() {
        let values = [-1.0, 0.0, 1.0, 0.0];
        for shift in 0..4 {
            let mut flow = cycle(vec![[1.0, 0.0, 0.0]]);
            for index in 0..=4 {
                let v = values[(index + shift) % 4];
                flow.push(axial(0.001 * v, 200e6 * v)).unwrap();
                flow.push(axial(0.001 * v, 200e6 * v)).unwrap();
            }
            close(
                flow.finish().unwrap().cycles_to_failure.unwrap(),
                analytical_life(0.001, 200e6),
            );
        }
        let mut warm = cycle(vec![[1.0, 0.0, 0.0]]);
        for index in 0..101 {
            let sign = if index % 2 == 0 { -1.0 } else { 1.0 };
            warm.push(axial(0.001 * sign, 200e6 * sign)).unwrap();
        }
        assert_eq!(warm.states[0].direction_changes, 3);
        let before = warm.clone();
        assert_eq!(
            warm.finish(),
            Err(CriticalPlaneError::NotSingleCycle { plane: 0 })
        );
        same_accepted(&warm, &before);
        let mut fresh = warm.fresh().unwrap();
        assert!(Arc::ptr_eq(&warm.curve, &fresh.curve));
        assert!(Arc::ptr_eq(&warm.planes, &fresh.planes));
        assert_eq!(fresh.sample_count(), 0);
        fresh.push_block(&axial_cycle(), &mut || false).unwrap();
        fresh.finish().unwrap();
    }

    #[test]
    fn g0_closure_is_full_tensor_explicit_and_reported_without_inserting_an_endpoint() {
        let mut exact = cycle(vec![[1.0, 0.0, 0.0]]);
        let mut samples = axial_cycle();
        samples[2].total_mechanical_strain[0] += 1e-8;
        samples[2].stress_pa[0] += 10.0;
        exact.push_block(&samples, &mut || false).unwrap();
        let before = exact.clone();
        assert!(matches!(
            exact.finish(),
            Err(CriticalPlaneError::NotClosed {
                quantity: "mechanical strain",
                ..
            })
        ));
        same_accepted(&exact, &before);
        let tolerance = CycleClosureTolerance {
            absolute_strain: 2e-8,
            absolute_stress_pa: 20.0,
            relative: 0.0,
        };
        let mut tolerant =
            SwtCycle::new(curve(), vec![[1.0, 0.0, 0.0]], limits(), tolerance).unwrap();
        tolerant.push_block(&samples, &mut || false).unwrap();
        let result = tolerant.finish().unwrap();
        assert_eq!(result.samples, 3);
        assert_eq!(result.closure, tolerance);
        close(
            result.maximum_strain_closure_defect,
            (samples[2].total_mechanical_strain[0] + 0.001).abs(),
        );
        assert_eq!(result.maximum_stress_closure_defect_pa, 10.0);
        assert_eq!(tolerant.last, Some(samples[2]));

        let mut hidden = cycle(vec![[1.0, 0.0, 0.0]]);
        let mut samples = axial_cycle();
        samples[2].stress_pa[4] = 1.0;
        hidden.push_block(&samples, &mut || false).unwrap();
        assert_eq!(
            hidden.finish(),
            Err(CriticalPlaneError::NotClosed {
                component: 4,
                quantity: "stress"
            })
        );
        let mut relative = SwtCycle::new(
            curve(),
            vec![[1.0, 0.0, 0.0]],
            limits(),
            CycleClosureTolerance {
                relative: 1e-4,
                ..exact_closure()
            },
        )
        .unwrap();
        let mut samples = axial_cycle();
        samples[2].stress_pa[0] += 10.0;
        relative.push_block(&samples, &mut || false).unwrap();
        assert_eq!(
            relative.finish().unwrap().maximum_stress_closure_defect_pa,
            10.0
        );
    }

    #[test]
    fn g4_push_block_and_finish_are_transactional_at_every_cancellation_boundary() {
        let mut base = cycle(vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
        base.push(axial_cycle()[0]).unwrap();
        for boundary in 0..5 {
            let mut attempt = base.clone();
            let mut calls = 0;
            let result = attempt.push_with_cancel(axial_cycle()[1], &mut || {
                calls += 1;
                calls > boundary
            });
            assert_eq!(result, Err(CriticalPlaneError::Cancelled));
            same_accepted(&attempt, &base);
            attempt.push(axial_cycle()[1]).unwrap();
            attempt.push(axial_cycle()[2]).unwrap();
            attempt.finish().unwrap();
        }
        let mut failed = base.clone();
        let mut invalid = axial_cycle()[2];
        invalid.stress_pa[2] = f64::NAN;
        assert!(matches!(
            failed.push_block(&[axial_cycle()[1], invalid], &mut || false),
            Err(CriticalPlaneError::InvalidSample { .. })
        ));
        same_accepted(&failed, &base);
        let mut calls = 0;
        assert_eq!(
            failed.push_block(&axial_cycle()[1..], &mut || {
                calls += 1;
                calls > 7
            }),
            Err(CriticalPlaneError::Cancelled)
        );
        same_accepted(&failed, &base);
        base.push_block(&axial_cycle()[1..], &mut || false).unwrap();
        for boundary in 0..5 {
            let mut attempt = base.clone();
            let mut calls = 0;
            assert_eq!(
                attempt.finish_with_cancel(&mut || {
                    calls += 1;
                    calls > boundary
                }),
                Err(CriticalPlaneError::Cancelled)
            );
            same_accepted(&attempt, &base);
        }
    }

    #[test]
    fn g0_model_support_and_limits_refuse_without_changing_the_observed_cycle() {
        for (amplitude, offset, expected_compression) in [(0.001, -300e6, true), (0.1, 0.0, false)]
        {
            let mut flow = cycle(vec![[1.0, 0.0, 0.0]]);
            for sign in [-1.0, 1.0, -1.0] {
                flow.push(axial(amplitude * sign, 200e6 * sign + offset))
                    .unwrap();
            }
            let before = flow.clone();
            let result = flow.finish();
            if expected_compression {
                assert!(matches!(
                    result,
                    Err(CriticalPlaneError::Fatigue(
                        FatigueError::NoTensileStress { .. }
                    ))
                ));
            } else {
                assert!(matches!(
                    result,
                    Err(CriticalPlaneError::Fatigue(
                        FatigueError::StrainAmplitudeOutsideCurve { .. }
                    ))
                ));
            }
            same_accepted(&flow, &before);
        }
        let mut zero = cycle(vec![[1.0, 0.0, 0.0]]);
        zero.push_block(&[axial(0.0, 0.0); 3], &mut || false)
            .unwrap();
        let result = zero.finish().unwrap();
        assert_eq!(result.cycles_to_failure, None);
        assert_eq!(result.damage_per_cycle, 0.0);
        for limited in [
            CriticalPlaneLimits {
                maximum_samples: 3,
                ..limits()
            },
            CriticalPlaneLimits {
                maximum_projection_work: 3,
                ..limits()
            },
        ] {
            let mut flow =
                SwtCycle::new(curve(), vec![[1.0, 0.0, 0.0]], limited, exact_closure()).unwrap();
            flow.push_block(&axial_cycle(), &mut || false).unwrap();
            let before = flow.clone();
            assert_eq!(flow.push(axial(0.0, 0.0)), Err(CriticalPlaneError::Budget));
            same_accepted(&flow, &before);
            flow.finish().unwrap();
        }
        assert!(SwtCycle::new(curve(), vec![], limits(), exact_closure()).is_err());
        assert!(SwtCycle::new(curve(), vec![[0.0; 3]], limits(), exact_closure()).is_err());
        for bad in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(
                SwtCycle::new(
                    curve(),
                    vec![[1.0, 0.0, 0.0]],
                    limits(),
                    CycleClosureTolerance {
                        relative: bad,
                        ..exact_closure()
                    }
                )
                .is_err()
            );
        }
        let mut short = cycle(vec![[1.0, 0.0, 0.0]]);
        assert!(short.finish().is_err());
        short.push(axial(0.0, 0.0)).unwrap();
        assert!(short.finish().is_err());
    }

    fn axial_tensor(direction: [f64; 3], strain: f64, stress: f64) -> TensorFatigueSample {
        let [x, y, z] = direction;
        let shape = [x * x, y * y, z * z, x * y, y * z, z * x];
        TensorFatigueSample {
            total_mechanical_strain: shape.map(|s| s * strain),
            stress_pa: shape.map(|s| s * stress),
        }
    }

    #[test]
    fn g1_hemisphere_grid_has_checked_counts_and_refines_an_oblique_axial_reference() {
        assert_eq!(
            hemisphere_planes(1, 3).unwrap(),
            vec![[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
        );
        assert_eq!(hemisphere_planes(2, 13).unwrap().len(), 13);
        assert_eq!(hemisphere_planes(2, 12), Err(CriticalPlaneError::Budget));
        assert!(hemisphere_planes(0, 10).is_err());
        assert_eq!(
            hemisphere_planes(usize::MAX, usize::MAX),
            Err(CriticalPlaneError::Budget)
        );
        let direction = NormalStressDirection::new([1.0, 2.0, 3.0])
            .unwrap()
            .unit_normal();
        let mut previous = 0.0;
        for resolution in [1, 2, 4, 8, 16] {
            let planes = hemisphere_planes(resolution, 2048).unwrap();
            let mut flow = cycle(planes);
            for sign in [-1.0, 1.0, -1.0] {
                flow.push(axial_tensor(direction, 0.001 * sign, 100e6 * sign))
                    .unwrap();
            }
            let result = flow.finish().unwrap();
            assert!(result.strain_amplitude + 1e-16 >= previous);
            assert!(result.strain_amplitude <= 0.001 * (1.0 + 1e-12));
            previous = result.strain_amplitude;
        }
        assert!(previous > 0.99 * 0.001);
    }
}
