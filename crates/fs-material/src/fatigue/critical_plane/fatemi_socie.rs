//! Fatemi-Socie shear-mode life for one declared closed mechanical cycle.
//!
//! Select maximum ENGINEERING shear-strain amplitude on the caller's fixed
//! plane set. On that plane, evaluate the maximum normal stress over the whole
//! cycle, not merely at a strain extremum. Exact amplitude ties select the
//! larger normal-stress maximum, then the earliest input plane. This is the
//! original maximum-shear-plane convention, not the extended maximum-product
//! criterion. The bounded shear calibration is
//!
//! `gamma_a * (1 + k*sigma_n_max/sigma_y0)
//!    = tau_f'/G * (2*N)^b_gamma + gamma_f' * (2*N)^c_gamma`.
//!
//! `N` is FULL CYCLES; the power-law argument `2*N` counts reversals. The shear
//! parameters are independently calibrated, never axial CMB values relabeled.
//! Signed and tensile-only normal-stress treatments are explicit alternatives.
//! See the original Fatemi/Socie plane convention and shear formulation:
//! <https://doi.org/10.1111/j.1460-2695.1988.tb01169.x>
//! <https://doc.comsol.com/6.4/doc/com.comsol.help.fatigue/fatigue_ug_sme.4.23.html>
//! The tensile-only treatment is also stated in:
//! <https://doi.org/10.1111/ffe.14153>.
//!
//! Input uses tensor-Voigt `[xx, yy, zz, xy, yz, zx]` in one fixed frame.
//! Strain is TOTAL MECHANICAL strain (elastic plus plastic), excluding thermal
//! eigenstrain and accumulated equivalent plastic strain. The engineering shear
//! vector on a plane is `2*(epsilon*n - n*(n^T*epsilon*n))`; its minimum enclosing
//! circle supplies the amplitude. A scalar shear magnitude or longest chord
//! would lose nonproportional loading information.
//!
//! This estimates life for a caller-declared resolved representative cycle.
//! It neither counts arbitrary multiaxial reversals nor infers repetitions or
//! accumulated damage. Full stress and strain endpoint closure is required but
//! does not establish cyclic stabilization, constitutive-model adequacy, or
//! calibration transfer. The native constitutive owner must supply the actual
//! stress response; no extra nonproportional hardening is invented here. The
//! finite plane search and floating-point circle remain numerical estimates.

pub mod source;

use super::shear_path::{
    CircleWork, SearchPlane, ShearPathError, dot6, enclosing_radius, normal_stress,
    tensor_normalization,
};
use super::{CriticalPlaneError, CycleClosureTolerance, TensorFatigueSample, closure_defect};
use crate::fatigue::strain_life::inverse;
use crate::fatigue::{FatigueError, FatigueScope};
use fs_math::det;
use std::sync::Arc;

/// Required convention for the crack-opening normal-stress term. These are
/// distinct calibrations, not an automatic fallback under compressive loading.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FatemiSocieNormalStress {
    /// Preserve the signed maximum. A nonpositive multiplier refuses; it is
    /// never reinterpreted as zero fatigue or infinite life.
    Signed,
    /// Only a tensile maximum augments the base shear parameter. Compression
    /// receives no additional fatigue-life credit.
    TensionOnly,
}

/// Shear strain-life coefficients and the separately calibrated FS correction.
/// G and the INITIAL yield reference are fixed material properties, not the
/// evolving tangent or hardened yield surface of a cyclic plasticity update.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FatemiSocieParameters {
    pub shear_modulus_pa: f64,
    pub shear_fatigue_strength_coefficient_pa: f64,
    pub shear_fatigue_strength_exponent: f64,
    pub shear_fatigue_ductility_coefficient: f64,
    pub shear_fatigue_ductility_exponent: f64,
    pub normal_stress_sensitivity: f64,
    pub initial_yield_stress_pa: f64,
    pub normal_stress: FatemiSocieNormalStress,
}

/// Immutable bounded shear calibration, preserving original source evidence.
#[derive(Debug, Clone, PartialEq)]
pub struct FatemiSocieCalibration {
    scope: FatigueScope,
    parameters: FatemiSocieParameters,
    minimum_cycles: f64,
    maximum_cycles: f64,
    elastic_coefficient: f64,
    log_minimum_reversals: f64,
    log_maximum_reversals: f64,
    minimum_parameter: f64,
    maximum_parameter: f64,
    source: Option<Arc<source::FatemiSocieCalibrationSource>>,
}
impl FatemiSocieCalibration {
    pub fn new(
        scope: FatigueScope,
        parameters: FatemiSocieParameters,
        minimum_cycles: f64,
        maximum_cycles: f64,
    ) -> Result<Self, FatemiSocieError> {
        scope.validate()?;
        let p = parameters;
        if !p.shear_modulus_pa.is_finite()
            || p.shear_modulus_pa <= 0.0
            || !p.initial_yield_stress_pa.is_finite()
            || p.initial_yield_stress_pa <= 0.0
            || !p.normal_stress_sensitivity.is_finite()
            || p.normal_stress_sensitivity < 0.0
            || !p.shear_fatigue_strength_coefficient_pa.is_finite()
            || p.shear_fatigue_strength_coefficient_pa < 0.0
            || !p.shear_fatigue_ductility_coefficient.is_finite()
            || p.shear_fatigue_ductility_coefficient < 0.0
            || (p.shear_fatigue_strength_coefficient_pa == 0.0
                && p.shear_fatigue_ductility_coefficient == 0.0)
            || !p.shear_fatigue_strength_exponent.is_finite()
            || p.shear_fatigue_strength_exponent >= 0.0
            || !p.shear_fatigue_ductility_exponent.is_finite()
            || p.shear_fatigue_ductility_exponent >= 0.0
        {
            return Err(FatemiSocieError::Invalid(
                "positive G/initial yield, nonnegative k/active shear coefficients, and negative shear exponents are required",
            ));
        }
        if !minimum_cycles.is_finite()
            || minimum_cycles < 0.5
            || !maximum_cycles.is_finite()
            || maximum_cycles <= minimum_cycles
        {
            return Err(FatemiSocieError::Invalid(
                "declare a finite shear-life interval beginning at or above half a cycle",
            ));
        }
        let elastic_coefficient = p.shear_fatigue_strength_coefficient_pa / p.shear_modulus_pa;
        if !elastic_coefficient.is_finite()
            || (p.shear_fatigue_strength_coefficient_pa > 0.0 && elastic_coefficient == 0.0)
        {
            return Err(FatemiSocieError::Invalid(
                "shear elastic fatigue coefficient is not representable",
            ));
        }
        let log_minimum_reversals = core::f64::consts::LN_2 + det::ln(minimum_cycles);
        let log_maximum_reversals = core::f64::consts::LN_2 + det::ln(maximum_cycles);
        let mut calibration = Self {
            scope,
            parameters,
            minimum_cycles,
            maximum_cycles,
            elastic_coefficient,
            log_minimum_reversals,
            log_maximum_reversals,
            minimum_parameter: 0.0,
            maximum_parameter: 0.0,
            source: None,
        };
        calibration.minimum_parameter = calibration.evaluate_log_reversals(log_maximum_reversals).0;
        calibration.maximum_parameter = calibration.evaluate_log_reversals(log_minimum_reversals).0;
        if log_minimum_reversals >= log_maximum_reversals
            || !calibration.minimum_parameter.is_finite()
            || calibration.minimum_parameter <= 0.0
            || !calibration.maximum_parameter.is_finite()
            || calibration.maximum_parameter <= calibration.minimum_parameter
        {
            return Err(FatemiSocieError::Invalid(
                "positive finite distinguishable shear-life endpoints are required",
            ));
        }
        Ok(calibration)
    }
    pub fn scope(&self) -> &FatigueScope {
        &self.scope
    }
    pub fn parameters(&self) -> FatemiSocieParameters {
        self.parameters
    }
    pub fn cycle_domain(&self) -> (f64, f64) {
        (self.minimum_cycles, self.maximum_cycles)
    }
    pub fn parameter_domain(&self) -> (f64, f64) {
        (self.minimum_parameter, self.maximum_parameter)
    }
    pub fn source(&self) -> Option<&source::FatemiSocieCalibrationSource> {
        self.source.as_deref()
    }

    /// The right-hand-side shear-life parameter, not a corrected observed
    /// strain amplitude. Under zero opening stress it equals engineering gamma_a.
    pub fn parameter_at_cycles(&self, cycles: f64) -> Result<f64, FatemiSocieError> {
        if !cycles.is_finite() || cycles < self.minimum_cycles || cycles > self.maximum_cycles {
            return Err(FatigueError::LifeOutsideCurve {
                cycles,
                minimum: self.minimum_cycles,
                maximum: self.maximum_cycles,
            }
            .into());
        }
        if cycles == self.minimum_cycles {
            return Ok(self.maximum_parameter);
        }
        if cycles == self.maximum_cycles {
            return Ok(self.minimum_parameter);
        }
        let value = self
            .evaluate_log_reversals(core::f64::consts::LN_2 + det::ln(cycles))
            .0;
        if !value.is_finite() || value <= 0.0 {
            return Err(FatigueError::Arithmetic.into());
        }
        Ok(value)
    }
    /// Invert the FS parameter only inside its original finite calibration.
    /// None means exactly zero modeled shear parameter, not an endurance limit
    /// or qualification of infinite component life.
    pub fn cycles_to_failure(&self, parameter: f64) -> Result<Option<f64>, FatemiSocieError> {
        inverse::invert_life(
            parameter,
            self.parameter_domain(),
            self.cycle_domain(),
            (self.log_minimum_reversals, self.log_maximum_reversals),
            |x| self.evaluate_log_reversals(x),
        )
        .map_err(Into::into)
    }
    fn evaluate_log_reversals(&self, x: f64) -> (f64, f64) {
        let p = self.parameters;
        let elastic = self.elastic_coefficient * det::exp(p.shear_fatigue_strength_exponent * x);
        let plastic = p.shear_fatigue_ductility_coefficient
            * det::exp(p.shear_fatigue_ductility_exponent * x);
        (
            elastic + plastic,
            p.shear_fatigue_strength_exponent * elastic
                + p.shear_fatigue_ductility_exponent * plastic,
        )
    }
    fn correction(&self, maximum_normal_stress_pa: f64) -> Result<f64, FatemiSocieError> {
        let p = self.parameters;
        let effective = match p.normal_stress {
            FatemiSocieNormalStress::Signed => maximum_normal_stress_pa,
            FatemiSocieNormalStress::TensionOnly => maximum_normal_stress_pa.max(0.0),
        };
        let increment = if p.normal_stress_sensitivity == 0.0 {
            0.0
        } else {
            let ratio = effective / p.initial_yield_stress_pa;
            if ratio.is_finite() && (ratio != 0.0 || effective == 0.0) {
                p.normal_stress_sensitivity * ratio
            } else {
                (p.normal_stress_sensitivity * effective) / p.initial_yield_stress_pa
            }
        };
        let factor = 1.0 + increment;
        if !factor.is_finite() {
            return Err(FatigueError::Arithmetic.into());
        }
        if factor <= 0.0 {
            return Err(FatemiSocieError::NonpositiveNormalCorrection {
                maximum_normal_stress_pa,
                factor,
            });
        }
        Ok(factor)
    }
}

/// Samples times planes is the projection work admitted before observation.
/// Circle work counts the same normalized point/predicate/construction units as
/// Findley. Storage is O(samples + planes); the bounded circle solver has a
/// cubic worst case and a deterministic fixed permutation for ordinary paths.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FatemiSocieLimits {
    pub maximum_samples: u64,
    pub maximum_planes: usize,
    pub maximum_projection_work: u64,
    pub maximum_circle_work: u64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FatemiSocieCycleResult {
    pub samples: u64,
    pub evaluated_planes: usize,
    pub selected_plane: usize,
    pub normal_in_frame: [f64; 3],
    pub engineering_shear_strain_amplitude: f64,
    pub minimum_normal_stress_pa: f64,
    /// Original signed maximum, retained even with TensionOnly treatment.
    pub maximum_normal_stress_pa: f64,
    pub normal_stress_correction: f64,
    pub fatemi_socie_parameter: f64,
    /// Full cycles on the bounded calibration; None only for zero parameter.
    pub cycles_to_failure: Option<f64>,
    pub projection_work: u64,
    pub circle_work: u64,
    pub closure: CycleClosureTolerance,
    pub maximum_strain_closure_defect: f64,
    pub maximum_stress_closure_defect_pa: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum FatemiSocieError {
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
    NotCyclic,
    NonpositiveNormalCorrection {
        maximum_normal_stress_pa: f64,
        factor: f64,
    },
    NumericalCircle,
    Fatigue(FatigueError),
    Closure(CriticalPlaneError),
    Allocation,
    Cancelled,
    Finished,
}
impl core::fmt::Display for FatemiSocieError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Fatemi-Socie cycle refused: {self:?}")
    }
}
impl std::error::Error for FatemiSocieError {}
impl From<FatigueError> for FatemiSocieError {
    fn from(error: FatigueError) -> Self {
        match error {
            FatigueError::Cancelled => Self::Cancelled,
            FatigueError::Allocation => Self::Allocation,
            error => Self::Fatigue(error),
        }
    }
}
impl From<ShearPathError> for FatemiSocieError {
    fn from(error: ShearPathError) -> Self {
        match error {
            ShearPathError::Budget => Self::Budget,
            ShearPathError::Cancelled => Self::Cancelled,
            ShearPathError::NotCyclic => Self::NotCyclic,
            ShearPathError::NumericalCircle => Self::NumericalCircle,
            ShearPathError::Fatigue(error) => error.into(),
        }
    }
}
impl From<CriticalPlaneError> for FatemiSocieError {
    fn from(error: CriticalPlaneError) -> Self {
        match error {
            CriticalPlaneError::NotClosed {
                component,
                quantity,
            } => Self::NotClosed {
                component,
                quantity,
            },
            CriticalPlaneError::Fatigue(error) => error.into(),
            error => Self::Closure(error),
        }
    }
}

/// Actual tensor observations stored once. Finalization projects one plane at
/// a time into reused scratch. Refusal never publishes a partial result; owning
/// native solvers use discard_after for cheap atomic gauge-bank rollback.
#[derive(Debug, Clone)]
pub struct FatemiSocieCycle {
    calibration: Arc<FatemiSocieCalibration>,
    planes: Arc<Vec<SearchPlane>>,
    limits: FatemiSocieLimits,
    closure: CycleClosureTolerance,
    samples: Vec<TensorFatigueSample>,
    result: Option<FatemiSocieCycleResult>,
}
impl FatemiSocieCycle {
    pub fn new(
        calibration: Arc<FatemiSocieCalibration>,
        normals: Vec<[f64; 3]>,
        limits: FatemiSocieLimits,
        closure: CycleClosureTolerance,
    ) -> Result<Self, FatemiSocieError> {
        if limits.maximum_samples < 3
            || limits.maximum_planes == 0
            || limits.maximum_projection_work == 0
            || limits.maximum_circle_work == 0
        {
            return Err(FatemiSocieError::Invalid(
                "positive work/plane limits and at least three observations are required",
            ));
        }
        if !closure.absolute_strain.is_finite()
            || closure.absolute_strain < 0.0
            || !closure.absolute_stress_pa.is_finite()
            || closure.absolute_stress_pa < 0.0
            || !closure.relative.is_finite()
            || closure.relative < 0.0
            || closure.relative > 1.0
        {
            return Err(FatemiSocieError::Invalid(
                "finite nonnegative closure tolerances and relative tolerance at most one are required",
            ));
        }
        if normals.is_empty() {
            return Err(FatemiSocieError::Invalid(
                "declare at least one search plane",
            ));
        }
        if normals.len() > limits.maximum_planes {
            return Err(FatemiSocieError::Budget);
        }
        let count = u64::try_from(normals.len()).map_err(|_| FatemiSocieError::Budget)?;
        if count > limits.maximum_projection_work {
            return Err(FatemiSocieError::Budget);
        }
        let mut planes = Vec::new();
        planes
            .try_reserve_exact(normals.len())
            .map_err(|_| FatemiSocieError::Allocation)?;
        for normal in normals {
            planes.push(SearchPlane::new(normal)?);
        }
        Ok(Self {
            calibration,
            planes: Arc::new(planes),
            limits,
            closure,
            samples: Vec::new(),
            result: None,
        })
    }
    pub fn calibration(&self) -> &FatemiSocieCalibration {
        &self.calibration
    }
    pub fn sample_count(&self) -> u64 {
        self.samples.len() as u64
    }
    pub fn plane_count(&self) -> usize {
        self.planes.len()
    }
    pub fn projection_work(&self) -> u64 {
        self.sample_count() * self.planes.len() as u64
    }
    pub fn circle_work(&self) -> u64 {
        self.result.map_or(0, |r| r.circle_work)
    }
    pub fn limits(&self) -> FatemiSocieLimits {
        self.limits
    }
    pub fn closure_tolerance(&self) -> CycleClosureTolerance {
        self.closure
    }
    pub fn result(&self) -> Option<&FatemiSocieCycleResult> {
        self.result.as_ref()
    }
    /// Empty window with the same immutable calibration, normals and policies.
    /// The consumer explicitly supplies the new first observation.
    pub fn fresh(&self) -> Result<Self, FatemiSocieError> {
        Ok(Self {
            calibration: self.calibration.clone(),
            planes: self.planes.clone(),
            limits: self.limits,
            closure: self.closure,
            samples: Vec::new(),
            result: None,
        })
    }
    /// Roll back an explicitly saved observed prefix and clear its candidate
    /// result. Ordinary native steps need not clone the accumulated history.
    pub fn discard_after(&mut self, sample_count: u64) -> Result<(), FatemiSocieError> {
        let count = usize::try_from(sample_count).map_err(|_| FatemiSocieError::Budget)?;
        if count > self.samples.len() {
            return Err(FatemiSocieError::Invalid(
                "cannot restore an unobserved prefix",
            ));
        }
        self.samples.truncate(count);
        self.result = None;
        Ok(())
    }
    pub fn push(&mut self, sample: TensorFatigueSample) -> Result<(), FatemiSocieError> {
        self.push_with_cancel(sample, &mut || false)
    }
    pub fn push_with_cancel(
        &mut self,
        sample: TensorFatigueSample,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), FatemiSocieError> {
        self.push_block(core::slice::from_ref(&sample), cancelled)
    }
    pub fn push_block(
        &mut self,
        samples: &[TensorFatigueSample],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), FatemiSocieError> {
        if cancelled() {
            return Err(FatemiSocieError::Cancelled);
        }
        if self.result.is_some() {
            return Err(FatemiSocieError::Finished);
        }
        self.admit_count(samples.len())?;
        let previous = self.samples.len();
        for (offset, sample) in samples.iter().enumerate() {
            if cancelled() {
                return Err(FatemiSocieError::Cancelled);
            }
            for (quantity, values) in [
                ("stress", sample.stress_pa),
                ("mechanical strain", sample.total_mechanical_strain),
            ] {
                for (component, value) in values.into_iter().enumerate() {
                    if !value.is_finite() {
                        return Err(FatemiSocieError::InvalidSample {
                            sample: (previous + offset) as u64,
                            component,
                            quantity,
                        });
                    }
                }
            }
        }
        self.samples
            .try_reserve(samples.len())
            .map_err(|_| FatemiSocieError::Allocation)?;
        for chunk in samples.chunks(64) {
            if cancelled() {
                self.samples.truncate(previous);
                return Err(FatemiSocieError::Cancelled);
            }
            self.samples.extend_from_slice(chunk);
        }
        if cancelled() {
            self.samples.truncate(previous);
            return Err(FatemiSocieError::Cancelled);
        }
        Ok(())
    }
    fn admit_count(&self, additional: usize) -> Result<(), FatemiSocieError> {
        let count = self
            .samples
            .len()
            .checked_add(additional)
            .ok_or(FatemiSocieError::Budget)?;
        let count = u64::try_from(count).map_err(|_| FatemiSocieError::Budget)?;
        let planes = u64::try_from(self.planes.len()).map_err(|_| FatemiSocieError::Budget)?;
        let work = count.checked_mul(planes).ok_or(FatemiSocieError::Budget)?;
        if count > self.limits.maximum_samples || work > self.limits.maximum_projection_work {
            return Err(FatemiSocieError::Budget);
        }
        Ok(())
    }
    pub fn finish(&mut self) -> Result<FatemiSocieCycleResult, FatemiSocieError> {
        self.finish_with_cancel(&mut || false)
    }
    pub fn finish_with_cancel(
        &mut self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FatemiSocieCycleResult, FatemiSocieError> {
        if cancelled() {
            return Err(FatemiSocieError::Cancelled);
        }
        if let Some(result) = self.result {
            return Ok(result);
        }
        if self.samples.len() < 3 {
            return Err(FatemiSocieError::Invalid(
                "a closed cycle requires at least three observations including both endpoints",
            ));
        }
        let first = self.samples[0];
        let last = self.samples[self.samples.len() - 1];
        let strain_defect = closure_defect(
            &first.total_mechanical_strain,
            &last.total_mechanical_strain,
            self.closure.absolute_strain,
            self.closure.relative,
            "mechanical strain",
        )?;
        let stress_defect = closure_defect(
            &first.stress_pa,
            &last.stress_pa,
            self.closure.absolute_stress_pa,
            self.closure.relative,
            "stress",
        )?;
        let (origin, scale) = tensor_normalization(
            self.samples.iter().map(|s| s.total_mechanical_strain),
            cancelled,
        )?;
        let mut path = Vec::new();
        path.try_reserve_exact(self.samples.len())
            .map_err(|_| FatemiSocieError::Allocation)?;
        let mut work = CircleWork {
            used: 0,
            maximum: self.limits.maximum_circle_work,
            cancelled,
        };
        let mut selected: Option<FatemiSocieCycleResult> = None;
        for (index, plane) in self.planes.iter().enumerate() {
            path.clear();
            let mut minimum = f64::INFINITY;
            let mut maximum = f64::NEG_INFINITY;
            for sample in &self.samples {
                work.check()?;
                let stress = normal_stress(plane.normal_weights, sample.stress_pa)?;
                minimum = minimum.min(stress);
                maximum = maximum.max(stress);
                let shifted = core::array::from_fn(|i| {
                    (sample.total_mechanical_strain[i] - origin[i]) / scale
                });
                let point = [
                    dot6(plane.shear_weights[0], shifted),
                    dot6(plane.shear_weights[1], shifted),
                ];
                if point.iter().any(|v| !v.is_finite()) {
                    return Err(FatigueError::Arithmetic.into());
                }
                path.push(point);
            }
            let tensor_amplitude = enclosing_radius(&mut path, &mut work)?;
            // Tensor-to-engineering shear is exactly a factor of two. Apply it
            // to the normalized radius before restoring physical strain scale.
            let amplitude = (2.0 * tensor_amplitude) * scale;
            if !amplitude.is_finite() || (tensor_amplitude > 0.0 && amplitude == 0.0) {
                return Err(FatigueError::Arithmetic.into());
            }
            if selected.as_ref().is_none_or(|r| {
                amplitude > r.engineering_shear_strain_amplitude
                    || (amplitude == r.engineering_shear_strain_amplitude
                        && maximum > r.maximum_normal_stress_pa)
            }) {
                selected = Some(FatemiSocieCycleResult {
                    samples: self.sample_count(),
                    evaluated_planes: self.planes.len(),
                    selected_plane: index,
                    normal_in_frame: plane.normal,
                    engineering_shear_strain_amplitude: amplitude,
                    minimum_normal_stress_pa: minimum,
                    maximum_normal_stress_pa: maximum,
                    normal_stress_correction: 0.0,
                    fatemi_socie_parameter: 0.0,
                    cycles_to_failure: None,
                    projection_work: self.projection_work(),
                    circle_work: 0,
                    closure: self.closure,
                    maximum_strain_closure_defect: strain_defect,
                    maximum_stress_closure_defect_pa: stress_defect,
                });
            }
        }
        work.check()?;
        let mut result = selected.ok_or(FatemiSocieError::Invalid("missing search planes"))?;
        result.normal_stress_correction = self
            .calibration
            .correction(result.maximum_normal_stress_pa)?;
        result.fatemi_socie_parameter =
            result.engineering_shear_strain_amplitude * result.normal_stress_correction;
        if !result.fatemi_socie_parameter.is_finite()
            || (result.engineering_shear_strain_amplitude > 0.0
                && result.fatemi_socie_parameter == 0.0)
        {
            return Err(FatigueError::Arithmetic.into());
        }
        result.cycles_to_failure = self
            .calibration
            .cycles_to_failure(result.fatemi_socie_parameter)?;
        result.circle_work = work.used;
        work.check()?;
        self.result = Some(result);
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scope() -> FatigueScope {
        FatigueScope {
            source_reference: "synthetic shear-life reference, not measured material".into(),
            material: "synthetic coupon".into(), process: "fixture".into(),
            surface_and_geometry: "smooth torsion specimen".into(), environment: "fixture".into(),
            load_basis: "total engineering shear strain amplitude; full-cycle life; explicit normal treatment".into(),
            temperature_k: 293.15,
        }
    }
    fn parameters(policy: FatemiSocieNormalStress) -> FatemiSocieParameters {
        FatemiSocieParameters {
            shear_modulus_pa: 80e9,
            shear_fatigue_strength_coefficient_pa: 200e6,
            shear_fatigue_strength_exponent: -0.5,
            shear_fatigue_ductility_coefficient: 0.3,
            shear_fatigue_ductility_exponent: -0.5,
            normal_stress_sensitivity: 0.5,
            initial_yield_stress_pa: 250e6,
            normal_stress: policy,
        }
    }
    fn calibration(policy: FatemiSocieNormalStress) -> Arc<FatemiSocieCalibration> {
        Arc::new(FatemiSocieCalibration::new(scope(), parameters(policy), 1.0, 1e6).unwrap())
    }
    fn limits() -> FatemiSocieLimits {
        FatemiSocieLimits {
            maximum_samples: 1000,
            maximum_planes: 1000,
            maximum_projection_work: 100_000,
            maximum_circle_work: 1_000_000,
        }
    }
    fn cycle(policy: FatemiSocieNormalStress, normals: Vec<[f64; 3]>) -> FatemiSocieCycle {
        FatemiSocieCycle::new(
            calibration(policy),
            normals,
            limits(),
            CycleClosureTolerance {
                absolute_strain: 0.0,
                absolute_stress_pa: 0.0,
                relative: 0.0,
            },
        )
        .unwrap()
    }
    fn sample(gamma_xy: f64, gamma_zx: f64, normal_stress: f64) -> TensorFatigueSample {
        TensorFatigueSample {
            stress_pa: [normal_stress, 0.0, 0.0, gamma_xy * 2e9, 0.0, gamma_zx * 2e9],
            total_mechanical_strain: [0.0, 0.0, 0.0, gamma_xy / 2.0, 0.0, gamma_zx / 2.0],
        }
    }
    fn line(amplitude: f64, normal_stress: f64) -> [TensorFatigueSample; 3] {
        [
            sample(-amplitude, 0.0, normal_stress),
            sample(amplitude, 0.0, normal_stress),
            sample(-amplitude, 0.0, normal_stress),
        ]
    }
    fn evaluate(
        model: &mut FatemiSocieCycle,
        samples: &[TensorFatigueSample],
    ) -> FatemiSocieCycleResult {
        model.push_block(samples, &mut || false).unwrap();
        model.finish().unwrap()
    }
    fn close(actual: f64, expected: f64, relative: f64) {
        assert!(
            (actual - expected).abs()
                <= relative * actual.abs().max(expected.abs()).max(f64::MIN_POSITIVE),
            "actual {actual:?}, expected {expected:?}"
        );
    }
    fn known_life(parameter: f64) -> f64 {
        0.5 * (0.3025 / parameter).powi(2)
    }

    #[test]
    fn shear_curve_inverse_matches_closed_form_in_full_cycles_and_strict_bounds() {
        let curve = calibration(FatemiSocieNormalStress::Signed);
        for cycles in [1.0, 10.0, 250.0, 10000.0, 1e6] {
            let expected = 0.3025 / det::sqrt(2.0 * cycles);
            let parameter = curve.parameter_at_cycles(cycles).unwrap();
            close(parameter, expected, 8e-15);
            close(
                curve.cycles_to_failure(parameter).unwrap().unwrap(),
                cycles,
                3e-13,
            );
        }
        let (minimum, maximum) = curve.parameter_domain();
        assert!(curve.cycles_to_failure(minimum / 2.0).is_err());
        assert!(curve.cycles_to_failure(maximum * 2.0).is_err());
        assert!(curve.parameter_at_cycles(0.5).is_err());
        assert_eq!(curve.cycles_to_failure(0.0).unwrap(), None);
        assert_eq!(curve.cycle_domain(), (1.0, 1e6));
    }

    #[test]
    fn inactive_shear_terms_and_invalid_calibrations_are_explicit() {
        for plastic_only in [false, true] {
            let mut p = parameters(FatemiSocieNormalStress::Signed);
            if plastic_only {
                p.shear_fatigue_strength_coefficient_pa = 0.0;
            } else {
                p.shear_fatigue_ductility_coefficient = 0.0;
            }
            let coefficient = if plastic_only { 0.3 } else { 0.0025 };
            let curve = FatemiSocieCalibration::new(scope(), p, 1.0, 1e6).unwrap();
            let parameter = coefficient / det::sqrt(200.0);
            close(
                curve.cycles_to_failure(parameter).unwrap().unwrap(),
                100.0,
                2e-13,
            );
        }
        for invalid in 0..5 {
            let mut p = parameters(FatemiSocieNormalStress::Signed);
            match invalid {
                0 => p.shear_modulus_pa = 0.0,
                1 => p.initial_yield_stress_pa = -1.0,
                2 => p.shear_fatigue_strength_exponent = 0.0,
                3 => p.normal_stress_sensitivity = -1.0,
                _ => {
                    p.shear_fatigue_strength_coefficient_pa = 0.0;
                    p.shear_fatigue_ductility_coefficient = 0.0;
                }
            }
            assert!(FatemiSocieCalibration::new(scope(), p, 1.0, 1e6).is_err());
        }
    }

    #[test]
    fn pure_torsion_preserves_tensor_to_engineering_factor_and_reversal_basis() {
        let mut model = cycle(FatemiSocieNormalStress::Signed, vec![[1.0, 0.0, 0.0]]);
        let result = evaluate(&mut model, &line(0.02, 0.0));
        close(result.engineering_shear_strain_amplitude, 0.02, 3e-14);
        close(result.fatemi_socie_parameter, 0.02, 3e-14);
        close(result.cycles_to_failure.unwrap(), known_life(0.02), 3e-13);
        assert_eq!(result.normal_stress_correction, 1.0);
        assert_eq!(result.samples, 3);
        assert_eq!(result.projection_work, 3);
        // Feeding tensor shear as engineering shear would quadruple this life.
        assert!(known_life(0.01) > 3.99 * result.cycles_to_failure.unwrap());
    }

    #[test]
    fn normal_peak_between_shear_extrema_shortens_life() {
        let samples = [
            sample(-0.02, 0.0, 0.0),
            sample(0.0, 0.0, 75e6),
            sample(0.02, 0.0, 0.0),
            sample(0.0, 0.0, 0.0),
            sample(-0.02, 0.0, 0.0),
        ];
        let mut model = cycle(FatemiSocieNormalStress::Signed, vec![[1.0, 0.0, 0.0]]);
        let result = evaluate(&mut model, &samples);
        assert_eq!(result.maximum_normal_stress_pa, 75e6);
        close(result.normal_stress_correction, 1.15, 3e-14);
        close(result.fatemi_socie_parameter, 0.023, 3e-14);
        close(result.cycles_to_failure.unwrap(), known_life(0.023), 3e-13);
        assert!(result.cycles_to_failure.unwrap() < known_life(0.02));
    }

    #[test]
    fn original_maximum_shear_plane_differs_from_maximum_product() {
        let mut p = parameters(FatemiSocieNormalStress::Signed);
        p.normal_stress_sensitivity = 1.0;
        let calibration = Arc::new(FatemiSocieCalibration::new(scope(), p, 1.0, 1e6).unwrap());
        let mut model = FatemiSocieCycle::new(
            calibration,
            vec![[det::sqrt(3.0), 1.0, 0.0], [1.0, 1.0, 0.0]],
            limits(),
            CycleClosureTolerance {
                absolute_strain: 0.0,
                absolute_stress_pa: 0.0,
                relative: 0.0,
            },
        )
        .unwrap();
        let state = |strain| TensorFatigueSample {
            stress_pa: [1e9, 0.0, 0.0, 0.0, 0.0, 0.0],
            total_mechanical_strain: [strain, 0.0, 0.0, 0.0, 0.0, 0.0],
        };
        let result = evaluate(&mut model, &[state(-0.02), state(0.02), state(-0.02)]);
        assert_eq!(result.selected_plane, 1);
        close(result.engineering_shear_strain_amplitude, 0.02, 3e-14);
        close(result.fatemi_socie_parameter, 0.06, 5e-14);
        // The 30-degree plane has smaller shear amplitude but greater product;
        // the original convention must nevertheless retain the 45-degree plane.
        let other_parameter = 0.02 * det::sqrt(3.0) / 2.0 * 4.0;
        assert!(other_parameter > result.fatemi_socie_parameter);
    }

    #[test]
    fn exact_shear_ties_select_greater_normal_stress_then_input_order() {
        let mut model = cycle(
            FatemiSocieNormalStress::Signed,
            vec![[1.0, -1.0, 0.0], [1.0, 1.0, 0.0], [1.0, 1.0, 0.0]],
        );
        let state = |strain| TensorFatigueSample {
            stress_pa: [0.0, 0.0, 0.0, 100e6, 0.0, 0.0],
            total_mechanical_strain: [strain, 0.0, 0.0, 0.0, 0.0, 0.0],
        };
        let result = evaluate(&mut model, &[state(-0.02), state(0.02), state(-0.02)]);
        assert_eq!(result.selected_plane, 1);
        close(result.maximum_normal_stress_pa, 100e6, 3e-14);
    }

    #[test]
    fn signed_and_tension_only_compression_are_distinct_and_never_zero_damage_clipping() {
        let mut signed = cycle(FatemiSocieNormalStress::Signed, vec![[1.0, 0.0, 0.0]]);
        let signed_result = evaluate(&mut signed, &line(0.02, -100e6));
        close(signed_result.normal_stress_correction, 0.8, 3e-14);
        close(
            signed_result.cycles_to_failure.unwrap(),
            known_life(0.016),
            3e-13,
        );
        let mut tensile = cycle(FatemiSocieNormalStress::TensionOnly, vec![[1.0, 0.0, 0.0]]);
        let tensile_result = evaluate(&mut tensile, &line(0.02, -1e9));
        assert_eq!(tensile_result.maximum_normal_stress_pa, -1e9);
        assert_eq!(tensile_result.normal_stress_correction, 1.0);
        close(
            tensile_result.cycles_to_failure.unwrap(),
            known_life(0.02),
            3e-13,
        );
        let mut refused = signed.fresh().unwrap();
        refused
            .push_block(&line(0.02, -500e6), &mut || false)
            .unwrap();
        assert!(matches!(
            refused.finish(),
            Err(FatemiSocieError::NonpositiveNormalCorrection { factor: 0.0, .. })
        ));
        assert!(refused.result().is_none());
    }

    #[test]
    fn rotating_and_triangular_engineering_shear_paths_use_shared_enclosing_circle() {
        let mut rotating = Vec::new();
        for i in 0..16 {
            let (s, c) = det::sin_cos(i as f64 * core::f64::consts::TAU / 16.0);
            rotating.push(sample(0.02 * c, 0.02 * s, 0.0));
        }
        rotating.push(rotating[0]);
        let triangle = [
            sample(0.02, 0.0, 0.0),
            sample(-0.01, 0.01 * det::sqrt(3.0), 0.0),
            sample(-0.01, -0.01 * det::sqrt(3.0), 0.0),
            sample(0.02, 0.0, 0.0),
        ];
        for samples in [&rotating[..], &triangle[..]] {
            let mut model = cycle(FatemiSocieNormalStress::Signed, vec![[1.0, 0.0, 0.0]]);
            let result = evaluate(&mut model, samples);
            close(result.engineering_shear_strain_amplitude, 0.02, 3e-14);
            close(result.cycles_to_failure.unwrap(), known_life(0.02), 3e-13);
        }
        let translated: Vec<_> = rotating
            .iter()
            .copied()
            .map(|mut s| {
                s.total_mechanical_strain[3] += 0.001;
                s.total_mechanical_strain[5] -= 0.004;
                s
            })
            .collect();
        let mut shifted = cycle(FatemiSocieNormalStress::Signed, vec![[1.0, 0.0, 0.0]]);
        close(
            evaluate(&mut shifted, &translated).engineering_shear_strain_amplitude,
            0.02,
            3e-14,
        );
    }

    fn rotate_tensor(t: [f64; 6], r: [[f64; 3]; 3]) -> [f64; 6] {
        let matrix = [[t[0], t[3], t[5]], [t[3], t[1], t[4]], [t[5], t[4], t[2]]];
        let mut result = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                for a in 0..3 {
                    for b in 0..3 {
                        result[i][j] += r[i][a] * matrix[a][b] * r[j][b];
                    }
                }
            }
        }
        [
            result[0][0],
            result[1][1],
            result[2][2],
            result[0][1],
            result[1][2],
            result[2][0],
        ]
    }
    #[test]
    fn simultaneous_frame_rotation_preserves_life() {
        let (s, c) = det::sin_cos(0.41);
        let rotation = [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]];
        let normals = vec![[1.0, 0.0, 0.0], [1.0, 1.0, 1.0], [0.0, 0.0, 1.0]];
        let rotated_normals = normals
            .iter()
            .map(|n| core::array::from_fn(|i| (0..3).map(|j| rotation[i][j] * n[j]).sum()))
            .collect();
        let history = [
            sample(-0.02, 0.0, 25e6),
            sample(0.0, 0.02, 75e6),
            sample(0.02, 0.0, 50e6),
            sample(0.0, -0.02, 0.0),
            sample(-0.02, 0.0, 25e6),
        ];
        let rotated: Vec<_> = history
            .iter()
            .map(|t| TensorFatigueSample {
                stress_pa: rotate_tensor(t.stress_pa, rotation),
                total_mechanical_strain: rotate_tensor(t.total_mechanical_strain, rotation),
            })
            .collect();
        let mut before = cycle(FatemiSocieNormalStress::Signed, normals);
        let mut after = cycle(FatemiSocieNormalStress::Signed, rotated_normals);
        let a = evaluate(&mut before, &history);
        let b = evaluate(&mut after, &rotated);
        assert_eq!(a.selected_plane, b.selected_plane);
        close(
            a.engineering_shear_strain_amplitude,
            b.engineering_shear_strain_amplitude,
            2e-13,
        );
        close(
            a.cycles_to_failure.unwrap(),
            b.cycles_to_failure.unwrap(),
            2e-12,
        );
    }

    #[test]
    fn full_tensor_closure_static_refusal_and_zero_shear_scope_are_explicit() {
        for strain_defect in [false, true] {
            let mut history = line(0.02, 0.0);
            if strain_defect {
                history[2].total_mechanical_strain[2] = 1e-5;
            } else {
                history[2].stress_pa[2] = 1.0;
            }
            let mut model = cycle(FatemiSocieNormalStress::Signed, vec![[1.0, 0.0, 0.0]]);
            model.push_block(&history, &mut || false).unwrap();
            assert!(matches!(
                model.finish(),
                Err(FatemiSocieError::NotClosed { component: 2, .. })
            ));
            assert!(model.result().is_none());
            model.closure.absolute_strain = 1e-5;
            model.closure.absolute_stress_pa = 1.0;
            let result = model.finish().unwrap();
            if strain_defect {
                assert_eq!(result.maximum_strain_closure_defect, 1e-5);
            } else {
                assert_eq!(result.maximum_stress_closure_defect_pa, 1.0);
            }
        }
        let mut constant = cycle(FatemiSocieNormalStress::Signed, vec![[1.0, 0.0, 0.0]]);
        constant
            .push_block(&[sample(0.01, 0.0, 0.0); 3], &mut || false)
            .unwrap();
        assert_eq!(constant.finish(), Err(FatemiSocieError::NotCyclic));
        let hydrostatic = |v| TensorFatigueSample {
            stress_pa: [v * 1e9, v * 1e9, v * 1e9, 0.0, 0.0, 0.0],
            total_mechanical_strain: [v, v, v, 0.0, 0.0, 0.0],
        };
        let mut hydro = constant.fresh().unwrap();
        let result = evaluate(
            &mut hydro,
            &[hydrostatic(-0.001), hydrostatic(0.001), hydrostatic(-0.001)],
        );
        assert_eq!(result.engineering_shear_strain_amplitude, 0.0);
        assert_eq!(result.cycles_to_failure, None);
    }

    #[test]
    fn exact_work_and_sample_limits_are_transactional() {
        let history = line(0.02, 0.0);
        let mut model = cycle(FatemiSocieNormalStress::Signed, vec![[1.0, 0.0, 0.0]]);
        let baseline = evaluate(&mut model, &history);
        let mut exact = model.fresh().unwrap();
        exact.limits.maximum_samples = 3;
        exact.limits.maximum_projection_work = 3;
        exact.limits.maximum_circle_work = baseline.circle_work;
        exact.push_block(&history, &mut || false).unwrap();
        assert_eq!(exact.push(history[0]), Err(FatemiSocieError::Budget));
        assert_eq!(exact.finish().unwrap(), baseline);
        let mut short = exact.fresh().unwrap();
        short.limits.maximum_circle_work -= 1;
        short.push_block(&history, &mut || false).unwrap();
        assert_eq!(short.finish(), Err(FatemiSocieError::Budget));
        assert!(short.result().is_none());
        let mut projection = exact.fresh().unwrap();
        projection.limits.maximum_projection_work = 2;
        assert_eq!(
            projection.push_block(&history, &mut || false),
            Err(FatemiSocieError::Budget)
        );
        assert_eq!(projection.sample_count(), 0);
    }

    #[test]
    fn cancellation_bad_samples_and_native_rollback_preserve_history() {
        let mut model = cycle(FatemiSocieNormalStress::Signed, vec![[1.0, 0.0, 0.0]]);
        let history = line(0.02, 0.0);
        model.push(history[0]).unwrap();
        let mut calls = 0;
        assert_eq!(
            model.push_block(&[history[1]; 10], &mut || {
                calls += 1;
                calls == 13
            }),
            Err(FatemiSocieError::Cancelled)
        );
        assert_eq!(model.sample_count(), 1);
        assert_eq!(model.samples[0], history[0]);
        let mut invalid = history[1];
        invalid.total_mechanical_strain[4] = f64::NAN;
        assert!(matches!(
            model.push(invalid),
            Err(FatemiSocieError::InvalidSample {
                sample: 1,
                component: 4,
                quantity: "mechanical strain"
            })
        ));
        assert_eq!(model.sample_count(), 1);
        model.push_block(&history[1..], &mut || false).unwrap();
        let mut calls = 0;
        assert_eq!(
            model.finish_with_cancel(&mut || {
                calls += 1;
                calls == 12
            }),
            Err(FatemiSocieError::Cancelled)
        );
        assert!(model.result().is_none());
        let original = model.finish().unwrap();
        assert_eq!(model.push(history[0]), Err(FatemiSocieError::Finished));
        assert!(model.discard_after(4).is_err());
        assert_eq!(model.result(), Some(&original));
        model.discard_after(2).unwrap();
        model.push(history[2]).unwrap();
        assert_eq!(model.finish().unwrap(), original);
        assert!(Arc::ptr_eq(
            &model.calibration,
            &model.fresh().unwrap().calibration
        ));
    }
}
