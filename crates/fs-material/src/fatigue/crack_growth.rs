//! Bounded Mode-I fatigue crack growth under ordered elastic load cycles.
//!
//! The Paris relation is written in dimensionally explicit normalized form:
//! `da/dN = reference_rate * (delta_K / reference_delta_K)^exponent`.
//! `delta_K` is the FULL stress-intensity RANGE in Pa sqrt(m), not amplitude.
//! Normalization avoids an exponent-dependent and easily misconverted unit
//! for the conventional Paris coefficient C. Parameters are caller-supplied
//! empirical calibration, with the existing material/process/environment scope.
//! An explicitly calibrated Walker option adds the normalized multiplier
//! `[(1-R_ref)/(1-R)]^[exponent*(1-gamma)]`. The reference rate then applies
//! at both the reference full range and reference load ratio. Gamma is the
//! exponent on `(1-R)` inside Walker's equivalent intensity, not the outer
//! crack-growth exponent. No Walker exponent is inferred from a Paris fit.
//!
//! Stress intensity is `K = Y(a) sigma sqrt(pi a)`. Constant geometry retains
//! an analytical O(1) fast path. Declared tabulated geometry uses log-log Y(a)
//! interpolation and analytically integrates each crossed power-law interval.
//! Work depends on crossed geometry intervals, not the cycle count. No
//! frozen-K/Euler approximation or cycle-by-cycle loop is used.
//! Blocks stay in physical order; completed rainflow cycles must NOT be used
//! as a reordered crack-growth spectrum.
//! A read-only inverse forecasts additional exposure to a declared larger crack
//! length under fixed tensile extrema, within the same domain and work limits.
//! The target is not an inferred failure length or a calibrated service life.
//!
//! This is the Paris region of long-crack, small-scale-yielding Mode-I LEFM.
//! Only nonnegative tensile stress ratios within declared support are admitted.
//! The lower fit bound is NOT a threshold: below-fit cycles refuse. No crack
//! closure, overload retardation, short-crack initiation, mixed-mode growth,
//! temperature transfer, resolved mesh crack, or stiffness/heat law is
//! inferred. The upper K bound is a declared validity cutoff, not a fracture
//! toughness measurement or a failure certificate.
//!
//! References: NASA-HDBK-5010 Volume 2 Revision A, §7.3, equations 7.3-1–3;
//! <https://standards.nasa.gov/standard/NASA/NASA-HDBK-5010-VOLUME_2>;
//! AFGROW Damage Tolerant Design Handbook, §5.1.2, equation 5.1.3, p. 5.1.9
//! (Walker convention and Paris-region applicability):
//! <https://www.afgrow.net/applications/dtdhandbook/pdfs/sec5_1_0.pdf>.

use std::sync::Arc;

use fs_math::det;

use super::FatigueScope;

/// Original material-card Paris calibration and retained query evidence.
pub mod source;

mod geometry;
pub use geometry::{
    CrackGeometry, CrackGeometryPoint, CrackGeometryTableContext, FiniteWidthCentreCrack,
    TabulatedCrackGeometry,
};

/// A normalized Paris fit; lengths are metres and N counts full cycles.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ParisParameters {
    /// Positive reference stress-intensity RANGE in Pa sqrt(m).
    pub reference_delta_k_pa_sqrt_m: f64,
    /// Positive crack extension at the reference range, in m/full cycle.
    /// With Walker correction, this is at its declared reference load ratio.
    pub reference_rate_m_per_cycle: f64,
    /// Positive dimensionless Paris exponent.
    pub exponent: f64,
}

/// Calibrated tensile Walker load-ratio dependence. AFGROW equation 5.1.3 is
/// `da/dN=C[(1-R)^gamma*K_max]^exponent`; equivalently, its intensity is
/// `Delta-K/(1-R)^(1-gamma)`. Only the declared `0 <= gamma <= 1` branch is
/// supported. Gamma=1 recovers Paris; gamma=0 depends on K_max alone within
/// the still-required positive-range fit domain. Neither is a static-growth law.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WalkerParameters {
    /// Load ratio at which `ParisParameters` declares its reference rate/range.
    pub reference_load_ratio: f64,
    /// Walker's inner `(1-R)` exponent; distinct from the Paris growth exponent.
    pub gamma: f64,
}

/// Supported fit/LEFM domain at the fixed temperature in `FatigueScope`.
/// The bounds apply to actual Delta-K and R, including when Walker is selected.
/// Declaring R support alone introduces no correction and does not infer support
/// between disjoint datasets. A Walker law needs its own calibrated parameters.
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
    walker: Option<WalkerParameters>,
    source: Option<Arc<source::ParisCalibrationSource>>,
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

/// Read-only analytical exposure to a caller-declared larger crack length.
/// Fractional full cycles are retained; neither failure nor service life is
/// inferred. Valid only for the fixed tensile extrema and existing calibration.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CrackGrowthForecast {
    pub before: CrackGrowthState,
    pub target_length_m: f64,
    pub additional_cycles: f64,
    pub initial_delta_k_pa_sqrt_m: f64,
    pub target_delta_k_pa_sqrt_m: f64,
    pub initial_maximum_k_pa_sqrt_m: f64,
    pub target_maximum_k_pa_sqrt_m: f64,
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
    /// The complete requested block would continue beyond the last declared
    /// geometry point. No crack length is extrapolated or partially published.
    GeometrySupportExhausted {
        length_m: f64,
        remaining_cycles: f64,
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

fn stress_range(minimum_pa: f64, maximum_pa: f64) -> Result<f64, CrackGrowthError> {
    if !minimum_pa.is_finite()
        || !maximum_pa.is_finite()
        || minimum_pa < 0.0
        || maximum_pa < minimum_pa
    {
        return Err(CrackGrowthError::StressOutsideDomain {
            minimum_pa,
            maximum_pa,
        });
    }
    Ok(maximum_pa - minimum_pa)
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
            walker: None,
            source: None,
        })
    }

    /// Declare a normalized Walker fit over the supplied actual-range/R domain.
    /// The reference rate/range must be admissible at the reference R itself.
    /// This creates caller-declared calibration; original material-card fits
    /// must be resolved through their dedicated source adapter instead.
    pub fn new_walker(
        scope: FatigueScope,
        parameters: ParisParameters,
        domain: ParisDomain,
        walker: WalkerParameters,
    ) -> Result<Self, CrackGrowthError> {
        let mut law = Self::new(scope, parameters, domain)?;
        if !walker.gamma.is_finite() || !(0.0..=1.0).contains(&walker.gamma) {
            return Err(CrackGrowthError::Invalid(
                "declare a calibrated Walker gamma in [0, 1]",
            ));
        }
        law.check_ratio(walker.reference_load_ratio)?;
        law.check_maximum_k(
            parameters.reference_delta_k_pa_sqrt_m / (1.0 - walker.reference_load_ratio),
        )?;
        law.walker = Some(walker);
        Ok(law)
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

    /// Explicit Walker calibration, or None for the unchanged plain Paris law.
    pub fn walker(&self) -> Option<WalkerParameters> {
        self.walker
    }

    /// Original immutable calibration, when resolved from a material card.
    pub fn source(&self) -> Option<&source::ParisCalibrationSource> {
        self.source.as_deref()
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
        let active_walker = self
            .walker
            .filter(|walker| walker.gamma != 1.0 && load_ratio != walker.reference_load_ratio);
        if delta_k_pa_sqrt_m == self.parameters.reference_delta_k_pa_sqrt_m
            && active_walker.is_none()
        {
            return Ok(self.parameters.reference_rate_m_per_cycle);
        }
        let mut log_ratio =
            det::ln(delta_k_pa_sqrt_m) - det::ln(self.parameters.reference_delta_k_pa_sqrt_m);
        if let Some(walker) = active_walker {
            // Form the small difference directly from R, before subtraction
            // from one could erase it. For separated ratios, subtract logs to
            // avoid loss near a relative change of -1.
            let relative = (load_ratio - walker.reference_load_ratio) / (1.0 - load_ratio);
            let log_ratio_shift = if relative.abs() < 0.5 {
                relative * log1p_over_x(relative)
            } else {
                det::ln(1.0 - walker.reference_load_ratio) - det::ln(1.0 - load_ratio)
            };
            log_ratio += (1.0 - walker.gamma) * log_ratio_shift;
        }
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
    geometry: CrackGeometry,
    limits: CrackGrowthLimits,
    state: CrackGrowthState,
}

impl ParisCrackGrowth {
    pub fn new(
        law: Arc<ParisLaw>,
        geometry: impl Into<CrackGeometry>,
        initial_length_m: f64,
        limits: CrackGrowthLimits,
    ) -> Result<Self, CrackGrowthError> {
        let geometry = geometry.into();
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
    pub fn geometry(&self) -> &CrackGeometry {
        &self.geometry
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

    /// Additional full-cycle equivalents to a strictly larger declared length,
    /// assuming these fixed tensile extrema and the declared geometry law.
    ///
    /// This read-only inverse checks every crossed geometry interval against
    /// the same R, Delta-K, maximum-K and geometry limits as advancement.
    /// It requires one remaining block and representable cumulative exposure
    /// within the cycle budget, without consuming either. A zero range refuses
    /// rather than implying infinite life. Nonrepresentable inverse arithmetic,
    /// including a rounded finite-cycle singularity, also refuses.
    pub fn cycles_to_length(
        &self,
        minimum_pa: f64,
        maximum_pa: f64,
        target_length_m: f64,
    ) -> Result<CrackGrowthForecast, CrackGrowthError> {
        self.cycles_to_length_with_cancel(minimum_pa, maximum_pa, target_length_m, &mut || false)
    }

    /// Cancellation-aware forecast, polling before/after the query and at each
    /// crossed table interval. Constant geometry retains constant work.
    pub fn cycles_to_length_with_cancel(
        &self,
        minimum_pa: f64,
        maximum_pa: f64,
        target_length_m: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<CrackGrowthForecast, CrackGrowthError> {
        if cancelled() {
            return Err(CrackGrowthError::Cancelled);
        }
        self.geometry.check_length(target_length_m)?;
        let a = self.state.crack_length_m;
        if target_length_m <= a {
            return Err(CrackGrowthError::Invalid(
                "forecast target must exceed the current crack length",
            ));
        }
        let range_pa = stress_range(minimum_pa, maximum_pa)?;
        if range_pa == 0.0 {
            return Err(CrackGrowthError::Invalid(
                "a growth forecast requires a positive tensile stress range",
            ));
        }
        let initial_maximum_k_pa_sqrt_m = self.stress_intensity(maximum_pa)?;
        let initial_delta_k_pa_sqrt_m = self.geometry.intensity(a, range_pa)?;
        let rate = self
            .law
            .growth_rate(initial_delta_k_pa_sqrt_m, minimum_pa / maximum_pa)?;
        let target_delta_k_pa_sqrt_m = self.geometry.intensity(target_length_m, range_pa)?;
        let target_maximum_k_pa_sqrt_m = self.geometry.intensity(target_length_m, maximum_pa)?;
        self.law.check_delta_k(target_delta_k_pa_sqrt_m)?;
        self.law.check_maximum_k(target_maximum_k_pa_sqrt_m)?;

        let additional_cycles = match &self.geometry {
            CrackGeometry::Constant(_) => {
                power_cycles(a, target_length_m, rate, 0.5 * self.law.parameters.exponent)?
            }
            CrackGeometry::Tabulated(table) => self.tabulated_cycles(
                table,
                minimum_pa,
                maximum_pa,
                target_length_m,
                rate,
                cancelled,
            )?,
        };
        self.checked_exposure(additional_cycles)?;
        if cancelled() {
            return Err(CrackGrowthError::Cancelled);
        }
        Ok(CrackGrowthForecast {
            before: self.state,
            target_length_m,
            additional_cycles,
            initial_delta_k_pa_sqrt_m,
            target_delta_k_pa_sqrt_m,
            initial_maximum_k_pa_sqrt_m,
            target_maximum_k_pa_sqrt_m,
        })
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
        let step = self.preview(cycle, cancelled)?;
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

    fn checked_exposure(&self, additional_cycles: f64) -> Result<(u64, f64), CrackGrowthError> {
        if !positive(additional_cycles) {
            return Err(CrackGrowthError::Invalid(
                "declare a positive finite full-cycle count",
            ));
        }
        let blocks = self
            .state
            .blocks
            .checked_add(1)
            .ok_or(CrackGrowthError::Budget)?;
        let cycles = self.state.cycles + additional_cycles;
        if blocks > self.limits.maximum_blocks
            || !cycles.is_finite()
            || cycles > self.limits.maximum_cycles
            || additional_cycles > self.limits.maximum_cycles - self.state.cycles
            || cycles <= self.state.cycles
        {
            return Err(CrackGrowthError::Budget);
        }
        Ok((blocks, cycles))
    }

    fn preview(
        &self,
        cycle: StressCycle,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<CrackGrowthStep, CrackGrowthError> {
        let range_pa = stress_range(cycle.minimum_pa, cycle.maximum_pa)?;
        let (blocks, cycles) = self.checked_exposure(cycle.cycles)?;
        let a = self.state.crack_length_m;
        let initial_maximum_k_pa_sqrt_m = self.stress_intensity(cycle.maximum_pa)?;
        let initial_delta_k_pa_sqrt_m = self.geometry.intensity(a, range_pa)?;
        let next_length = if range_pa == 0.0 {
            // A static admissible stress is not a fatigue cycle. Still retain
            // the caller's exposure count and enforce its instantaneous K.
            a
        } else {
            let ratio = cycle.minimum_pa / cycle.maximum_pa;
            let rate = self.law.growth_rate(initial_delta_k_pa_sqrt_m, ratio)?;
            let next = match &self.geometry {
                CrackGeometry::Constant(_) => {
                    power_length(a, rate, 0.5 * self.law.parameters.exponent, cycle.cycles)?
                }
                CrackGeometry::Tabulated(table) => {
                    self.tabulated_length(table, cycle, rate, cancelled)?
                }
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
            // Constant-Y intensities increase with length; table evolution has
            // additionally checked all crossed interval extrema, including any
            // interior intensity peak or valley.
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

/// Exact evolution for da/dN = rate_at_a * (length/a)^power. Table geometry
/// changes `power` to m*(1/2 + d ln(Y)/d ln(a)); the constant-Y path is m/2.
fn power_length(a: f64, rate: f64, power: f64, cycles: f64) -> Result<f64, CrackGrowthError> {
    let t = (rate / a) * cycles;
    let q = 1.0 - power;
    if !positive(t) || !q.is_finite() {
        return Err(CrackGrowthError::Arithmetic);
    }
    let x = q * t;
    if x <= -1.0 {
        return Err(CrackGrowthError::GrowthSingularity);
    }
    if !x.is_finite() {
        return Err(CrackGrowthError::Arithmetic);
    }
    let log_growth = t * log1p_over_x(x);
    Ok(if log_growth < 0.5 {
        a + a * det::expm1(log_growth)
    } else {
        det::exp(det::ln(a) + log_growth)
    })
}

/// Integral to an interval boundary. Positive infinity means a finite block
/// cannot reach that boundary; the caller can still evolve within the interval.
fn power_cycles_to_boundary(
    a: f64,
    target: f64,
    rate: f64,
    power: f64,
) -> Result<f64, CrackGrowthError> {
    let q = 1.0 - power;
    let log_growth = positive_log_ratio(a, target);
    let relative_rate = rate / a;
    if !q.is_finite() || !positive(log_growth) || !positive(relative_rate) {
        return Err(CrackGrowthError::Arithmetic);
    }
    let normalized = if q == 0.0 {
        log_growth
    } else {
        det::expm1(q * log_growth) / q
    };
    let cycles = normalized / relative_rate;
    if cycles.is_nan() || cycles <= 0.0 {
        return Err(CrackGrowthError::Arithmetic);
    }
    Ok(cycles)
}

fn power_cycles(a: f64, target: f64, rate: f64, power: f64) -> Result<f64, CrackGrowthError> {
    let cycles = power_cycles_to_boundary(a, target, rate, power)?;
    let forward_time = (rate / a) * cycles;
    let forward_argument = (1.0 - power) * forward_time;
    if !positive(cycles)
        || !positive(forward_time)
        || !forward_argument.is_finite()
        || forward_argument <= -1.0
    {
        return Err(CrackGrowthError::Arithmetic);
    }
    Ok(cycles)
}

/// ln(to/from), keeping adjacent representable values and wide ratios.
fn positive_log_ratio(from: f64, to: f64) -> f64 {
    let relative = (to - from) / from;
    if relative.is_finite() && relative > -1.0 {
        relative * log1p_over_x(relative)
    } else {
        det::ln(to) - det::ln(from)
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
