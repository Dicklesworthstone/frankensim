//! Bounded Mode-I fatigue crack growth under ordered elastic load cycles.
//!
//! The Paris relation is written in dimensionally explicit normalized form:
//! `da/dN = reference_rate * (delta_K / reference_delta_K)^exponent`.
//! `delta_K` is the FULL stress-intensity RANGE in Pa sqrt(m), not amplitude.
//! Normalization avoids an exponent-dependent and easily misconverted unit
//! for the conventional Paris coefficient C. Parameters are caller-supplied
//! empirical calibration, with the existing material/process/environment scope.
//!
//! For a declared constant geometry factor, `K = Y sigma sqrt(pi a)`. Each
//! constant-amplitude block integrates the evolving crack length analytically,
//! including the exponential exponent-two case. Work is constant per block,
//! independent of its cycle count. No frozen-K/Euler approximation is used.
//! Blocks stay in physical order; completed rainflow cycles must NOT be used
//! as a reordered crack-growth spectrum.
//!
//! This is the Paris region of long-crack, small-scale-yielding Mode-I LEFM.
//! Only nonnegative tensile stress ratios within declared support are admitted.
//! The lower fit bound is NOT a threshold: below-fit cycles refuse. No crack
//! closure, overload retardation, short-crack initiation, mixed-mode growth,
//! temperature transfer, evolving Y, mesh crack, or stiffness/heat law is
//! inferred. The upper K bound is a declared validity cutoff, not a fracture
//! toughness measurement or a failure certificate.
//!
//! References: NASA-HDBK-5010 Volume 2 Revision A, §7.3, equations 7.3-1–3;
//! <https://standards.nasa.gov/standard/NASA/NASA-HDBK-5010-VOLUME_2>;
//! AFGROW Damage Tolerant Design Handbook, §5.1.2 (Paris-region applicability).

use std::sync::Arc;

use fs_math::det;

use super::FatigueScope;

/// A normalized Paris fit; lengths are metres and N counts full cycles.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ParisParameters {
    /// Positive reference stress-intensity RANGE in Pa sqrt(m).
    pub reference_delta_k_pa_sqrt_m: f64,
    /// Positive crack extension at the reference range, in m/full cycle.
    pub reference_rate_m_per_cycle: f64,
    /// Positive dimensionless Paris exponent.
    pub exponent: f64,
}

/// Supported fit/LEFM domain at the fixed temperature in `FatigueScope`.
/// R support declares where this same fit is applicable; it does not introduce
/// a Walker/Forman correction or infer support between disjoint datasets.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ParisDomain {
    /// Inclusive lower Paris-region bound, NOT an arrest threshold.
    pub minimum_delta_k_pa_sqrt_m: f64,
    /// Inclusive upper Paris-region bound.
    pub maximum_delta_k_pa_sqrt_m: f64,
    /// Inclusive supported `K_min/K_max = sigma_min/sigma_max`, at least zero.
    pub minimum_load_ratio: f64,
    /// Inclusive supported R, strictly less than one.
    pub maximum_load_ratio: f64,
    /// Exclusive upper K_max validity limit in Pa sqrt(m).
    pub maximum_k_pa_sqrt_m: f64,
}

/// Immutable empirical calibration, never silently extrapolated.
#[derive(Debug, Clone, PartialEq)]
pub struct ParisLaw {
    scope: FatigueScope,
    parameters: ParisParameters,
    domain: ParisDomain,
}

/// Declared constant-Y geometry and the interval where it may be used.
/// The caller owns its relation to the nominal stress and crack geometry.
/// For a central through crack, `a` is the HALF length and Y=1 is the infinite-
/// plate approximation. A finite-width factor cannot be frozen as a crack
/// grows unless that approximation is justified over this complete interval.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConstantGeometryCrack {
    pub geometry_factor: f64,
    pub minimum_length_m: f64,
    pub maximum_length_m: f64,
}

/// Explicit cumulative work and cycle limits for one runtime crack.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CrackGrowthLimits {
    pub maximum_blocks: u64,
    pub maximum_cycles: f64,
}

/// One physical constant-amplitude block, in chronological spectrum order.
/// Fractional cycle equivalents are permitted for analytical continuation;
/// this interface does not assign half cycles to an unresolved stress trace.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StressCycle {
    pub minimum_pa: f64,
    pub maximum_pa: f64,
    pub cycles: f64,
}

/// Accepted internal crack state. No failure probability or stiffness is implied.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CrackGrowthState {
    pub crack_length_m: f64,
    pub cycles: f64,
    pub blocks: u64,
}

/// Original and advanced crack states plus actual endpoint stress intensities.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CrackGrowthStep {
    pub before: CrackGrowthState,
    pub after: CrackGrowthState,
    pub initial_delta_k_pa_sqrt_m: f64,
    pub final_delta_k_pa_sqrt_m: f64,
    pub initial_maximum_k_pa_sqrt_m: f64,
    pub final_maximum_k_pa_sqrt_m: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum CrackGrowthError {
    Invalid(&'static str),
    StressOutsideDomain {
        minimum_pa: f64,
        maximum_pa: f64,
    },
    LoadRatioOutsideDomain {
        ratio: f64,
        minimum: f64,
        maximum: f64,
    },
    DeltaKOutsideDomain {
        delta_k_pa_sqrt_m: f64,
        minimum: f64,
        maximum: f64,
    },
    MaximumKExceeded {
        maximum_k_pa_sqrt_m: f64,
        limit: f64,
    },
    CrackLengthOutsideDomain {
        length_m: f64,
        minimum: f64,
        maximum: f64,
    },
    /// The requested block reaches the mathematical finite-cycle singularity.
    /// This is outside the finite geometry/fit domain, not a rupture prediction.
    GrowthSingularity,
    Budget,
    Arithmetic,
    Cancelled,
}

impl core::fmt::Display for CrackGrowthError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Paris crack growth refused: {self:?}")
    }
}
impl std::error::Error for CrackGrowthError {}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

impl ParisLaw {
    pub fn new(
        scope: FatigueScope,
        parameters: ParisParameters,
        domain: ParisDomain,
    ) -> Result<Self, CrackGrowthError> {
        scope.validate().map_err(|_| {
            CrackGrowthError::Invalid(
                "declare the source, material, process, geometry, environment and temperature",
            )
        })?;
        if !positive(parameters.reference_delta_k_pa_sqrt_m)
            || !positive(parameters.reference_rate_m_per_cycle)
            || !positive(parameters.exponent)
            || !positive(domain.minimum_delta_k_pa_sqrt_m)
            || !positive(domain.maximum_delta_k_pa_sqrt_m)
            || !positive(domain.maximum_k_pa_sqrt_m)
            || domain.minimum_delta_k_pa_sqrt_m >= domain.maximum_delta_k_pa_sqrt_m
            || parameters.reference_delta_k_pa_sqrt_m < domain.minimum_delta_k_pa_sqrt_m
            || parameters.reference_delta_k_pa_sqrt_m > domain.maximum_delta_k_pa_sqrt_m
            || !domain.minimum_load_ratio.is_finite()
            || !domain.maximum_load_ratio.is_finite()
            || domain.minimum_load_ratio < 0.0
            || domain.minimum_load_ratio > domain.maximum_load_ratio
            || domain.maximum_load_ratio >= 1.0
            || parameters.reference_delta_k_pa_sqrt_m / (1.0 - domain.minimum_load_ratio)
                >= domain.maximum_k_pa_sqrt_m
        {
            return Err(CrackGrowthError::Invalid(
                "positive Paris parameters and nonempty tensile fit/LEFM support are required",
            ));
        }
        Ok(Self {
            scope,
            parameters,
            domain,
        })
    }

    pub fn scope(&self) -> &FatigueScope {
        &self.scope
    }
    pub fn parameters(&self) -> ParisParameters {
        self.parameters
    }
    pub fn domain(&self) -> ParisDomain {
        self.domain
    }

    /// da/dN in metres per full cycle, evaluated only in the declared support.
    /// A zero range is not evaluated through the positive-range Paris fit.
    pub fn growth_rate(
        &self,
        delta_k_pa_sqrt_m: f64,
        load_ratio: f64,
    ) -> Result<f64, CrackGrowthError> {
        self.check_ratio(load_ratio)?;
        self.check_delta_k(delta_k_pa_sqrt_m)?;
        self.check_maximum_k(delta_k_pa_sqrt_m / (1.0 - load_ratio))?;
        if delta_k_pa_sqrt_m == self.parameters.reference_delta_k_pa_sqrt_m {
            return Ok(self.parameters.reference_rate_m_per_cycle);
        }
        let log_ratio =
            det::ln(delta_k_pa_sqrt_m) - det::ln(self.parameters.reference_delta_k_pa_sqrt_m);
        let log_rate = det::ln(self.parameters.reference_rate_m_per_cycle)
            + self.parameters.exponent * log_ratio;
        let rate = det::exp(log_rate);
        if !positive(rate) {
            return Err(CrackGrowthError::Arithmetic);
        }
        Ok(rate)
    }

    fn check_ratio(&self, ratio: f64) -> Result<(), CrackGrowthError> {
        if !ratio.is_finite()
            || ratio < self.domain.minimum_load_ratio
            || ratio > self.domain.maximum_load_ratio
        {
            return Err(CrackGrowthError::LoadRatioOutsideDomain {
                ratio,
                minimum: self.domain.minimum_load_ratio,
                maximum: self.domain.maximum_load_ratio,
            });
        }
        Ok(())
    }

    fn check_delta_k(&self, delta_k_pa_sqrt_m: f64) -> Result<(), CrackGrowthError> {
        if !delta_k_pa_sqrt_m.is_finite()
            || delta_k_pa_sqrt_m < self.domain.minimum_delta_k_pa_sqrt_m
            || delta_k_pa_sqrt_m > self.domain.maximum_delta_k_pa_sqrt_m
        {
            return Err(CrackGrowthError::DeltaKOutsideDomain {
                delta_k_pa_sqrt_m,
                minimum: self.domain.minimum_delta_k_pa_sqrt_m,
                maximum: self.domain.maximum_delta_k_pa_sqrt_m,
            });
        }
        Ok(())
    }

    fn check_maximum_k(&self, maximum_k_pa_sqrt_m: f64) -> Result<(), CrackGrowthError> {
        if !maximum_k_pa_sqrt_m.is_finite()
            || maximum_k_pa_sqrt_m >= self.domain.maximum_k_pa_sqrt_m
        {
            return Err(CrackGrowthError::MaximumKExceeded {
                maximum_k_pa_sqrt_m,
                limit: self.domain.maximum_k_pa_sqrt_m,
            });
        }
        Ok(())
    }
}

impl ConstantGeometryCrack {
    fn validate(self) -> Result<(), CrackGrowthError> {
        if !positive(self.geometry_factor)
            || !positive(self.minimum_length_m)
            || !positive(self.maximum_length_m)
            || self.minimum_length_m >= self.maximum_length_m
        {
            return Err(CrackGrowthError::Invalid(
                "declare a positive constant Y and finite positive crack-length interval",
            ));
        }
        Ok(())
    }

    fn check_length(self, length_m: f64) -> Result<(), CrackGrowthError> {
        if !length_m.is_finite()
            || length_m < self.minimum_length_m
            || length_m > self.maximum_length_m
        {
            return Err(CrackGrowthError::CrackLengthOutsideDomain {
                length_m,
                minimum: self.minimum_length_m,
                maximum: self.maximum_length_m,
            });
        }
        Ok(())
    }

    fn intensity(self, length_m: f64, stress_pa: f64) -> Result<f64, CrackGrowthError> {
        // sqrt(pi)*sqrt(a) avoids overflowing pi*a for otherwise finite input.
        let value = (stress_pa * det::sqrt(length_m))
            * (self.geometry_factor * det::sqrt(std::f64::consts::PI));
        if !value.is_finite() || (stress_pa > 0.0 && value <= 0.0) {
            return Err(CrackGrowthError::Arithmetic);
        }
        Ok(value)
    }
}

/// Restartable analytical crack-growth owner; all updates are transactional.
/// The immutable calibration is shared and no load history is retained.
#[derive(Debug, Clone, PartialEq)]
pub struct ParisCrackGrowth {
    law: Arc<ParisLaw>,
    geometry: ConstantGeometryCrack,
    limits: CrackGrowthLimits,
    state: CrackGrowthState,
}

impl ParisCrackGrowth {
    pub fn new(
        law: Arc<ParisLaw>,
        geometry: ConstantGeometryCrack,
        initial_length_m: f64,
        limits: CrackGrowthLimits,
    ) -> Result<Self, CrackGrowthError> {
        geometry.validate()?;
        geometry.check_length(initial_length_m)?;
        if limits.maximum_blocks == 0 || !positive(limits.maximum_cycles) {
            return Err(CrackGrowthError::Invalid(
                "positive block and full-cycle budgets are required",
            ));
        }
        Ok(Self {
            law,
            geometry,
            limits,
            state: CrackGrowthState {
                crack_length_m: initial_length_m,
                cycles: 0.0,
                blocks: 0,
            },
        })
    }

    pub fn law(&self) -> &ParisLaw {
        &self.law
    }
    pub fn geometry(&self) -> ConstantGeometryCrack {
        self.geometry
    }
    pub fn limits(&self) -> CrackGrowthLimits {
        self.limits
    }
    pub fn state(&self) -> CrackGrowthState {
        self.state
    }

    /// Check one actual tensile nominal stress against the current crack and
    /// upper K limit, even if no cycle has closed yet. Static overload cannot
    /// be hidden in a constant or monotone observation window.
    pub fn stress_intensity(&self, stress_pa: f64) -> Result<f64, CrackGrowthError> {
        if !stress_pa.is_finite() || stress_pa < 0.0 {
            return Err(CrackGrowthError::StressOutsideDomain {
                minimum_pa: stress_pa,
                maximum_pa: stress_pa,
            });
        }
        let value = self
            .geometry
            .intensity(self.state.crack_length_m, stress_pa)?;
        self.law.check_maximum_k(value)?;
        Ok(value)
    }

    pub fn advance(&mut self, cycle: StressCycle) -> Result<CrackGrowthStep, CrackGrowthError> {
        self.advance_with_cancel(cycle, &mut || false)
    }

    /// Advance one complete constant-amplitude block, publishing nothing if any
    /// part of its evolving crack interval leaves the model's support.
    pub fn advance_with_cancel(
        &mut self,
        cycle: StressCycle,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<CrackGrowthStep, CrackGrowthError> {
        if cancelled() {
            return Err(CrackGrowthError::Cancelled);
        }
        let step = self.preview(cycle)?;
        if cancelled() {
            return Err(CrackGrowthError::Cancelled);
        }
        self.state = step.after;
        Ok(step)
    }

    /// Ordered spectrum transaction, with one bounded cancellation point before
    /// and after each analytical block. A late refusal rolls back every block.
    pub fn advance_block(
        &mut self,
        cycles: &[StressCycle],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<CrackGrowthState, CrackGrowthError> {
        if cancelled() {
            return Err(CrackGrowthError::Cancelled);
        }
        let count = u64::try_from(cycles.len()).map_err(|_| CrackGrowthError::Budget)?;
        if count > self.limits.maximum_blocks - self.state.blocks {
            return Err(CrackGrowthError::Budget);
        }
        let mut candidate = self.clone();
        for cycle in cycles {
            candidate.advance_with_cancel(*cycle, cancelled)?;
        }
        if cancelled() {
            return Err(CrackGrowthError::Cancelled);
        }
        self.state = candidate.state;
        Ok(self.state)
    }

    fn preview(&self, cycle: StressCycle) -> Result<CrackGrowthStep, CrackGrowthError> {
        if !cycle.minimum_pa.is_finite()
            || !cycle.maximum_pa.is_finite()
            || cycle.minimum_pa < 0.0
            || cycle.maximum_pa < cycle.minimum_pa
        {
            return Err(CrackGrowthError::StressOutsideDomain {
                minimum_pa: cycle.minimum_pa,
                maximum_pa: cycle.maximum_pa,
            });
        }
        if !positive(cycle.cycles) {
            return Err(CrackGrowthError::Invalid(
                "declare a positive finite full-cycle count",
            ));
        }
        let blocks = self
            .state
            .blocks
            .checked_add(1)
            .ok_or(CrackGrowthError::Budget)?;
        let cycles = self.state.cycles + cycle.cycles;
        if blocks > self.limits.maximum_blocks
            || !cycles.is_finite()
            || cycles > self.limits.maximum_cycles
            || cycles <= self.state.cycles
        {
            return Err(CrackGrowthError::Budget);
        }
        let a = self.state.crack_length_m;
        let initial_maximum_k_pa_sqrt_m = self.stress_intensity(cycle.maximum_pa)?;
        let range_pa = cycle.maximum_pa - cycle.minimum_pa;
        let initial_delta_k_pa_sqrt_m = self.geometry.intensity(a, range_pa)?;
        let next_length = if range_pa == 0.0 {
            // A static admissible stress is not a fatigue cycle. Still retain
            // the caller's exposure count and enforce its instantaneous K.
            a
        } else {
            let ratio = cycle.minimum_pa / cycle.maximum_pa;
            let rate = self.law.growth_rate(initial_delta_k_pa_sqrt_m, ratio)?;
            let t = (rate / a) * cycle.cycles;
            if !positive(t) {
                return Err(CrackGrowthError::Arithmetic);
            }
            let q = 1.0 - 0.5 * self.law.parameters.exponent;
            let x = q * t;
            if x <= -1.0 {
                return Err(CrackGrowthError::GrowthSingularity);
            }
            if !x.is_finite() {
                return Err(CrackGrowthError::Arithmetic);
            }
            // Integrating da/dN = rate*(a/a_initial)^(m/2) gives
            // ln(a_next/a_initial) = ln(1+q*t)/q. The continuous q=0 limit
            // is t; this form also stays accurate for exponents near two.
            let log_growth = t * log1p_over_x(x);
            let next = if log_growth < 0.5 {
                a + a * det::expm1(log_growth)
            } else {
                det::exp(det::ln(a) + log_growth)
            };
            self.geometry.check_length(next)?;
            if !next.is_finite() || next <= a {
                return Err(CrackGrowthError::Arithmetic);
            }
            next
        };
        let final_delta_k_pa_sqrt_m = self.geometry.intensity(next_length, range_pa)?;
        let final_maximum_k_pa_sqrt_m = self.geometry.intensity(next_length, cycle.maximum_pa)?;
        self.law.check_maximum_k(final_maximum_k_pa_sqrt_m)?;
        if range_pa > 0.0 {
            // a and both intensities increase monotonically throughout this
            // block, so endpoint checks admit the complete growth interval.
            self.law.check_delta_k(final_delta_k_pa_sqrt_m)?;
        }
        Ok(CrackGrowthStep {
            before: self.state,
            after: CrackGrowthState {
                crack_length_m: next_length,
                cycles,
                blocks,
            },
            initial_delta_k_pa_sqrt_m,
            final_delta_k_pa_sqrt_m,
            initial_maximum_k_pa_sqrt_m,
            final_maximum_k_pa_sqrt_m,
        })
    }
}

/// ln(1+x)/x with its continuous limit at zero. Correct the rounded 1+x
/// before division, so small cycle counts and m near two do not lose growth.
fn log1p_over_x(x: f64) -> f64 {
    let one_plus_x = 1.0 + x;
    if one_plus_x == 1.0 {
        1.0
    } else {
        det::ln(one_plus_x) / (one_plus_x - 1.0)
    }
}
