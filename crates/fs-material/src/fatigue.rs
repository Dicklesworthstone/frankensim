//! Streaming rainflow and stress-life fatigue (life-ladder-gahl, first rung).
//!
//! This is an empirical Palmgren-Miner baseline, not a constitutive softening
//! law, a multiaxial failure criterion, or a validated component-life claim.
//! Supply a SIGNED stress history in Pa, a fully reversed S-N calibration and
//! an explicitly named mean-stress policy. Never feed unsigned von Mises
//! stress as a substitute for a reversing tensile/compressive history.
//!
//! Monotone runs and plateaus are reduced incrementally. Unclosed reversals
//! survive arbitrary input blocks; call `finish` only at the end of the whole
//! observation, not at callback/window boundaries. Failed samples, blocks and
//! finalization do not alter accepted history or damage. Work per sample is
//! bounded by the explicit residue limit; no entire stress trace is retained.

use fs_math::det;
use std::sync::Arc;

/// Conditions identifying the supplied empirical calibration. These strings
/// are caller declarations, not independent verification of the source.
#[derive(Debug, Clone, PartialEq)]
pub struct FatigueScope {
    pub source_reference: String,
    pub material: String,
    pub process: String,
    pub surface_and_geometry: String,
    pub environment: String,
    pub load_basis: String,
    pub temperature_k: f64,
}

/// One fully reversed S-N observation: stress AMPLITUDE, not stress range.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SnPoint {
    pub amplitude_pa: f64,
    pub cycles_to_failure: f64,
}

/// Original calibration points and scope, retained without extrapolation.
#[derive(Debug, Clone, PartialEq)]
pub struct SnCurve {
    scope: FatigueScope,
    points: Vec<SnPoint>,
    endurance_limit_pa: Option<f64>,
}

/// Named treatment of the cycle's signed mean stress.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MeanStressCorrection {
    /// Explicit mean-independent approximation. Does not imply its validation.
    Uncorrected,
    /// Modified Goodman for tensile means; compressive means receive no credit.
    /// The declared ultimate strength also bounds the tensile cycle peak.
    GoodmanTensionOnly { ultimate_strength_pa: f64 },
}

/// Logical storage/work limits, checked before accepting a sample or reversal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RainflowLimits {
    pub maximum_samples: u64,
    pub maximum_residue: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub enum FatigueError {
    InvalidCalibration(&'static str),
    InvalidLimits,
    InvalidStress { sample: u64 },
    InvalidDirection,
    SampleLimit { limit: u64 },
    ResidueLimit { limit: usize },
    AmplitudeOutsideCurve { amplitude_pa: f64, minimum_pa: f64, maximum_pa: f64 },
    StaticStrengthExceeded { peak_pa: f64, ultimate_strength_pa: f64 },
    Arithmetic,
    Allocation,
    Cancelled,
    AlreadyFinished,
}

impl core::fmt::Display for FatigueError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "fatigue evaluation refused: {self:?}")
    }
}
impl std::error::Error for FatigueError {}

impl SnCurve {
    /// Points must have strictly increasing positive amplitudes and strictly
    /// decreasing finite lives >= 1 cycle. Interpolation is linear in log-log
    /// coordinates. No below/above-table extrapolation is ever performed.
    ///
    /// An optional endurance limit must be strictly below the lowest finite-
    /// life point. The gap between that limit and the first point remains
    /// UNKNOWN and refuses; it is not bridged by an invented infinite-life fit.
    pub fn new(
        scope: FatigueScope,
        points: Vec<SnPoint>,
        endurance_limit_pa: Option<f64>,
    ) -> Result<Self, FatigueError> {
        if [&scope.source_reference, &scope.material, &scope.process,
            &scope.surface_and_geometry, &scope.environment, &scope.load_basis]
            .iter().any(|s| s.trim().is_empty())
            || !scope.temperature_k.is_finite() || scope.temperature_k <= 0.0
        {
            return Err(FatigueError::InvalidCalibration("declare the source and all calibration conditions"));
        }
        if points.len() < 2 || points.iter().any(|p|
            !p.amplitude_pa.is_finite() || p.amplitude_pa <= 0.0
                || !p.cycles_to_failure.is_finite() || p.cycles_to_failure < 1.0)
        {
            return Err(FatigueError::InvalidCalibration("at least two positive finite S-N points with life >= 1 are required"));
        }
        if points.windows(2).any(|p|
            p[0].amplitude_pa >= p[1].amplitude_pa
                || p[0].cycles_to_failure <= p[1].cycles_to_failure
                || det::ln(p[0].amplitude_pa) >= det::ln(p[1].amplitude_pa))
        {
            return Err(FatigueError::InvalidCalibration("amplitudes must increase in resolvable log coordinates and lives must decrease"));
        }
        if endurance_limit_pa.is_some_and(|s|
            !s.is_finite() || s <= 0.0 || s >= points[0].amplitude_pa)
        {
            return Err(FatigueError::InvalidCalibration("endurance limit must be positive and strictly below finite-life data"));
        }
        Ok(Self { scope, points, endurance_limit_pa })
    }

    pub fn scope(&self) -> &FatigueScope { &self.scope }
    pub fn points(&self) -> &[SnPoint] { &self.points }
    pub fn endurance_limit_pa(&self) -> Option<f64> { self.endurance_limit_pa }

    /// `None` means zero amplitude or the explicitly declared endurance regime,
    /// not a demonstrated infinite service life. Everything else needs support
    /// from the original calibration table.
    pub fn cycles_to_failure(&self, amplitude_pa: f64) -> Result<Option<f64>, FatigueError> {
        if !amplitude_pa.is_finite() || amplitude_pa < 0.0 { return Err(FatigueError::Arithmetic); }
        if amplitude_pa == 0.0 || self.endurance_limit_pa.is_some_and(|s| amplitude_pa <= s) {
            return Ok(None);
        }
        let first = self.points[0];
        let last = self.points[self.points.len() - 1];
        if amplitude_pa < first.amplitude_pa || amplitude_pa > last.amplitude_pa {
            return Err(FatigueError::AmplitudeOutsideCurve {
                amplitude_pa, minimum_pa: first.amplitude_pa, maximum_pa: last.amplitude_pa,
            });
        }
        let upper = self.points.partition_point(|p| p.amplitude_pa < amplitude_pa);
        if self.points[upper].amplitude_pa == amplitude_pa {
            return Ok(Some(self.points[upper].cycles_to_failure));
        }
        let lo = self.points[upper - 1];
        let hi = self.points[upper];
        let t = (det::ln(amplitude_pa) - det::ln(lo.amplitude_pa))
            / (det::ln(hi.amplitude_pa) - det::ln(lo.amplitude_pa));
        let life = det::exp((1.0 - t) * det::ln(lo.cycles_to_failure)
            + t * det::ln(hi.cycles_to_failure));
        if !life.is_finite() || life < 1.0 { return Err(FatigueError::Arithmetic); }
        Ok(Some(life))
    }
}

impl MeanStressCorrection {
    fn validate(self) -> Result<(), FatigueError> {
        if let Self::GoodmanTensionOnly { ultimate_strength_pa: s } = self {
            if !s.is_finite() || s <= 0.0 {
                return Err(FatigueError::InvalidCalibration("Goodman requires a positive finite sourced ultimate strength"));
            }
        }
        Ok(())
    }

    /// Equivalent fully reversed amplitude. Compressive mean never increases
    /// predicted life in the tension-only Goodman policy.
    pub fn corrected_amplitude(self, amplitude_pa: f64, mean_pa: f64) -> Result<f64, FatigueError> {
        self.validate()?;
        if !amplitude_pa.is_finite() || amplitude_pa < 0.0 || !mean_pa.is_finite() {
            return Err(FatigueError::Arithmetic);
        }
        let value = match self {
            Self::Uncorrected => amplitude_pa,
            Self::GoodmanTensionOnly { ultimate_strength_pa: s } => {
                let peak = mean_pa + amplitude_pa;
                if !peak.is_finite() || peak >= s {
                    return Err(FatigueError::StaticStrengthExceeded { peak_pa: peak, ultimate_strength_pa: s });
                }
                amplitude_pa / (1.0 - mean_pa.max(0.0) / s)
            }
        };
        if !value.is_finite() { return Err(FatigueError::Arithmetic); }
        Ok(value)
    }
}

/// Damage already accounted for. Before `finish`, unresolved endpoint half
/// cycles are NOT included. `miner_damage >= 1` is only the empirical rule's
/// nominal exhaustion threshold; it is not a failure certificate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FatigueSummary {
    pub samples: u64,
    pub full_cycles: u64,
    pub half_cycles: u64,
    pub miner_damage: f64,
    pub maximum_corrected_amplitude_pa: f64,
    pub finished: bool,
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct TurningPoint { stress: f64 }

#[derive(Debug, Clone, Copy, PartialEq, Default)]
struct Tally {
    full: u64,
    half: u64,
    damage: f64,
    compensation: f64,
    maximum_amplitude: f64,
}

/// Incremental three-point rainflow reduction and Palmgren-Miner accumulation.
/// The shared curve is immutable; cloning a counter copies only bounded residue.
#[derive(Debug, Clone, PartialEq)]
pub struct RainflowMiner {
    curve: Arc<SnCurve>,
    correction: MeanStressCorrection,
    limits: RainflowLimits,
    residue: Vec<TurningPoint>,
    pending: Option<TurningPoint>,
    direction: i8,
    samples: u64,
    tally: Tally,
    finished: bool,
}

impl RainflowMiner {
    pub fn new(curve: Arc<SnCurve>, correction: MeanStressCorrection, limits: RainflowLimits) -> Result<Self, FatigueError> {
        correction.validate()?;
        if limits.maximum_samples == 0 || limits.maximum_residue < 2 { return Err(FatigueError::InvalidLimits); }
        Ok(Self { curve, correction, limits, residue: Vec::new(), pending: None,
            direction: 0, samples: 0, tally: Tally::default(), finished: false })
    }

    pub fn curve(&self) -> &SnCurve { &self.curve }
    pub fn pending_reversals(&self) -> usize { self.residue.len() }
    pub fn summary(&self) -> FatigueSummary {
        FatigueSummary { samples: self.samples, full_cycles: self.tally.full,
            half_cycles: self.tally.half, miner_damage: self.tally.damage,
            maximum_corrected_amplitude_pa: self.tally.maximum_amplitude, finished: self.finished }
    }

    pub fn push(&mut self, stress_pa: f64) -> Result<(), FatigueError> {
        self.push_with_cancel(stress_pa, &mut || false)
    }

    /// Transactional single sample. Cancellation is polled before acceptance
    /// and at every candidate cycle; existing accepted samples are untouched.
    pub fn push_with_cancel(&mut self, stress_pa: f64, cancelled: &mut impl FnMut() -> bool) -> Result<(), FatigueError> {
        if cancelled() { return Err(FatigueError::Cancelled); }
        if self.finished { return Err(FatigueError::AlreadyFinished); }
        if !stress_pa.is_finite() { return Err(FatigueError::InvalidStress { sample: self.samples }); }
        if self.samples >= self.limits.maximum_samples {
            return Err(FatigueError::SampleLimit { limit: self.limits.maximum_samples });
        }
        let point = TurningPoint { stress: stress_pa };
        let mut direction = self.direction;
        if let Some(previous) = self.pending {
            if point.stress != previous.stress {
                let next_direction = if point.stress > previous.stress { 1 } else { -1 };
                if direction == 0 || direction != next_direction {
                    self.append_reversal(previous, cancelled)?;
                }
                direction = next_direction;
            }
        }
        self.pending = Some(point);
        self.direction = direction;
        self.samples += 1;
        Ok(())
    }

    /// Atomic block ingestion; callback size does not change the reduction or
    /// damage summation order. Do not finish the counter between blocks.
    pub fn push_block(&mut self, stresses_pa: &[f64], cancelled: &mut impl FnMut() -> bool) -> Result<(), FatigueError> {
        if cancelled() { return Err(FatigueError::Cancelled); }
        if self.finished { return Err(FatigueError::AlreadyFinished); }
        let count = u64::try_from(stresses_pa.len()).map_err(|_| FatigueError::Arithmetic)?;
        if count > self.limits.maximum_samples - self.samples {
            return Err(FatigueError::SampleLimit { limit: self.limits.maximum_samples });
        }
        let mut candidate = self.clone();
        for &stress in stresses_pa { candidate.push_with_cancel(stress, cancelled)?; }
        if cancelled() { return Err(FatigueError::Cancelled); }
        *self = candidate;
        Ok(())
    }

    /// Account for the final endpoint and remaining half cycles ONCE. Repeated
    /// calls are idempotent. Refusal leaves the still-open counter unchanged.
    pub fn finish(&mut self) -> Result<FatigueSummary, FatigueError> {
        self.finish_with_cancel(&mut || false)
    }

    pub fn finish_with_cancel(&mut self, cancelled: &mut impl FnMut() -> bool) -> Result<FatigueSummary, FatigueError> {
        if cancelled() { return Err(FatigueError::Cancelled); }
        if self.finished { return Ok(self.summary()); }
        let mut candidate = self.clone();
        if let Some(last) = candidate.pending { candidate.append_reversal(last, cancelled)?; }
        let mut tally = candidate.tally;
        for pair in candidate.residue.windows(2) {
            if cancelled() { return Err(FatigueError::Cancelled); }
            tally = candidate.count_cycle(tally, pair[0], pair[1], false)?;
        }
        if cancelled() { return Err(FatigueError::Cancelled); }
        candidate.tally = tally;
        candidate.residue.clear();
        candidate.pending = None;
        candidate.finished = true;
        *self = candidate;
        Ok(self.summary())
    }

    fn count_cycle(&self, mut tally: Tally, a: TurningPoint, b: TurningPoint, full: bool) -> Result<Tally, FatigueError> {
        if a.stress == b.stress { return Ok(tally); }
        let amplitude = half_range(a.stress, b.stress)?;
        let mean = 0.5 * a.stress + 0.5 * b.stress;
        let corrected = self.correction.corrected_amplitude(amplitude, mean)?;
        let weight = if full { 1.0 } else { 0.5 };
        let damage = match self.curve.cycles_to_failure(corrected)? {
            Some(life) => weight / life,
            None => 0.0,
        };
        let corrected_term = damage - tally.compensation;
        let sum = tally.damage + corrected_term;
        let compensation = (sum - tally.damage) - corrected_term;
        if !sum.is_finite() || !compensation.is_finite() { return Err(FatigueError::Arithmetic); }
        if full { tally.full = tally.full.checked_add(1).ok_or(FatigueError::Arithmetic)?; }
        else { tally.half = tally.half.checked_add(1).ok_or(FatigueError::Arithmetic)?; }
        tally.damage = sum;
        tally.compensation = compensation;
        tally.maximum_amplitude = tally.maximum_amplitude.max(corrected);
        Ok(tally)
    }

    // Preview reductions and damage without mutating the accepted stack. A new
    // extreme can remove arbitrarily many nested pairs, but at most the explicit
    // residue limit. No clone/allocation proportional to trace length is needed.
    fn append_reversal(&mut self, point: TurningPoint, cancelled: &mut impl FnMut() -> bool) -> Result<(), FatigueError> {
        let mut keep = self.residue.len();
        let mut start_half = None;
        let mut tally = self.tally;
        while keep >= 2 {
            if cancelled() { return Err(FatigueError::Cancelled); }
            let a = self.residue[keep - 2];
            let b = self.residue[keep - 1];
            if half_range(b.stress, point.stress)? < half_range(a.stress, b.stress)? { break; }
            tally = self.count_cycle(tally, a, b, keep > 2)?;
            if keep == 2 { start_half = Some(b); keep = 0; break; }
            keep -= 2;
        }
        let required = if start_half.is_some() { 2 } else { keep.checked_add(1).ok_or(FatigueError::Arithmetic)? };
        if required > self.limits.maximum_residue {
            return Err(FatigueError::ResidueLimit { limit: self.limits.maximum_residue });
        }
        if cancelled() { return Err(FatigueError::Cancelled); }
        if required > self.residue.capacity() {
            self.residue.try_reserve_exact(required - self.residue.len()).map_err(|_| FatigueError::Allocation)?;
        }
        if let Some(start) = start_half {
            self.residue.clear();
            self.residue.push(start);
        } else { self.residue.truncate(keep); }
        self.residue.push(point);
        self.tally = tally;
        Ok(())
    }
}

fn half_range(a: f64, b: f64) -> Result<f64, FatigueError> {
    let difference = (a - b).abs();
    let result = if difference.is_finite() { 0.5 * difference } else { (0.5 * a - 0.5 * b).abs() };
    if !result.is_finite() || (a != b && result == 0.0) { return Err(FatigueError::Arithmetic); }
    Ok(result)
}

/// Fixed-plane signed normal stress, not a multiaxial fatigue criterion.
/// The direction and input tensor must be expressed in the SAME fixed frame.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NormalStressDirection([f64; 3]);

impl NormalStressDirection {
    pub fn new(direction: [f64; 3]) -> Result<Self, FatigueError> {
        if direction.iter().any(|x| !x.is_finite()) { return Err(FatigueError::InvalidDirection); }
        let scale = direction.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
        if scale == 0.0 { return Err(FatigueError::InvalidDirection); }
        let scaled = direction.map(|x| x / scale);
        let norm = det::sqrt(scaled.iter().map(|x| x * x).sum());
        Ok(Self(scaled.map(|x| x / norm)))
    }

    pub fn unit_normal(self) -> [f64; 3] { self.0 }

    /// Tensor-Voigt order `[xx, yy, zz, xy, yz, zx]`, NOT Mandel or engineering shear.
    pub fn project(self, stress_pa: &[f64; 6]) -> Result<f64, FatigueError> {
        if stress_pa.iter().any(|x| !x.is_finite()) { return Err(FatigueError::Arithmetic); }
        let [x, y, z] = self.0;
        let weights = [x*x, y*y, z*z, 2.0*x*y, 2.0*y*z, 2.0*z*x];
        let stress: f64 = weights.iter().zip(stress_pa).map(|(w, s)| w * s).sum();
        if !stress.is_finite() { return Err(FatigueError::Arithmetic); }
        Ok(stress)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scope() -> FatigueScope {
        FatigueScope { source_reference: "synthetic hand-check; not measured".into(),
            material: "synthetic coupon".into(), process: "fixture".into(),
            surface_and_geometry: "smooth axial coupon".into(), environment: "fixture".into(),
            load_basis: "fully reversed signed axial stress".into(), temperature_k: 293.15 }
    }
    fn curve() -> Arc<SnCurve> {
        Arc::new(SnCurve::new(scope(), vec![
            SnPoint { amplitude_pa: 1.0, cycles_to_failure: 1000.0 },
            SnPoint { amplitude_pa: 10.0, cycles_to_failure: 100.0 },
        ], None).unwrap())
    }
    fn counter() -> RainflowMiner {
        RainflowMiner::new(curve(), MeanStressCorrection::Uncorrected,
            RainflowLimits { maximum_samples: 10000, maximum_residue: 64 }).unwrap()
    }
    fn close(actual: f64, expected: f64) {
        assert!((actual - expected).abs() <= 1e-11 * expected.abs().max(1.0), "{actual} != {expected}");
    }

    #[test]
    fn sn_is_log_log_with_exact_knots_and_explicit_domain() {
        let curve = curve();
        assert_eq!(curve.cycles_to_failure(1.0).unwrap(), Some(1000.0));
        assert_eq!(curve.cycles_to_failure(10.0).unwrap(), Some(100.0));
        close(curve.cycles_to_failure(5.0).unwrap().unwrap(), 200.0);
        assert!(matches!(curve.cycles_to_failure(0.5), Err(FatigueError::AmplitudeOutsideCurve { .. })));
        assert!(curve.cycles_to_failure(11.0).is_err());
        assert!(curve.cycles_to_failure(f64::NAN).is_err());
        let endurance = SnCurve::new(scope(), curve.points().to_vec(), Some(0.25)).unwrap();
        assert_eq!(endurance.cycles_to_failure(0.25).unwrap(), None);
        assert!(endurance.cycles_to_failure(0.5).is_err());
        assert!(SnCurve::new(scope(), curve.points().to_vec(), Some(1.0)).is_err());
        let mut reversed = curve.points().to_vec(); reversed.reverse();
        assert!(SnCurve::new(scope(), reversed, None).is_err());
    }

    #[test]
    fn constant_and_monotone_records_have_no_fabricated_cycles() {
        let mut constant = counter();
        constant.push_block(&[3.0; 20], &mut || false).unwrap();
        assert_eq!(constant.finish().unwrap().miner_damage, 0.0);
        let mut ramp = counter();
        ramp.push_block(&[-2.0, -1.0, 0.0, 1.0, 2.0], &mut || false).unwrap();
        let summary = ramp.finish().unwrap();
        assert_eq!((summary.full_cycles, summary.half_cycles), (0, 1));
        close(summary.miner_damage, 0.001);
    }

    #[test]
    fn hand_counted_nested_history_and_chunking_are_identical() {
        // [-2,1] and [1,-3] are initial half cycles; [-1,3] is one
        // nested full cycle. Remaining half-cycle ranges are 8,9,8,6.
        let samples = [-2.0, 1.0, -3.0, 5.0, -1.0, 3.0, -4.0, 4.0, -2.0];
        let mut whole = counter();
        whole.push_block(&samples, &mut || false).unwrap();
        let summary = whole.finish().unwrap();
        assert_eq!((summary.full_cycles, summary.half_cycles), (1, 6));
        // N(a)=1000/a: D = (2 + .5*(1.5+2+4+4.5+4+3))/1000.
        close(summary.miner_damage, 0.0115);
        for split in 0..=samples.len() {
            let mut chunked = counter();
            chunked.push_block(&samples[..split], &mut || false).unwrap();
            chunked.push_block(&samples[split..], &mut || false).unwrap();
            assert_eq!(chunked.finish().unwrap(), summary);
        }
        let mut scalar = counter();
        for sample in samples { scalar.push(sample).unwrap(); }
        assert_eq!(scalar.finish().unwrap(), summary);
    }

    #[test]
    fn plateaus_do_not_double_count_and_finalization_is_once_only() {
        let mut c = counter();
        c.push_block(&[-2.0, -2.0, 2.0, 2.0, -2.0, -2.0], &mut || false).unwrap();
        let first = c.finish().unwrap();
        assert_eq!((first.full_cycles, first.half_cycles), (0, 2));
        close(first.miner_damage, 0.002);
        assert_eq!(c.finish().unwrap(), first);
        assert_eq!(c.push(1.0), Err(FatigueError::AlreadyFinished));
    }

    #[test]
    fn goodman_is_named_tensile_only_and_refuses_static_failure() {
        let rule = MeanStressCorrection::GoodmanTensionOnly { ultimate_strength_pa: 20.0 };
        close(rule.corrected_amplitude(3.0, 5.0).unwrap(), 4.0);
        close(rule.corrected_amplitude(3.0, -5.0).unwrap(), 3.0);
        assert!(matches!(rule.corrected_amplitude(3.0, 17.0), Err(FatigueError::StaticStrengthExceeded { .. })));
        let mut c = RainflowMiner::new(curve(), rule,
            RainflowLimits { maximum_samples: 100, maximum_residue: 16 }).unwrap();
        c.push_block(&[2.0, 8.0, 2.0], &mut || false).unwrap();
        close(c.finish().unwrap().miner_damage, 0.004);
    }

    #[test]
    fn unsupported_cycle_and_nonfinite_blocks_preserve_accepted_state() {
        let mut c = counter();
        c.push_block(&[-2.0, 2.0], &mut || false).unwrap();
        let before = c.clone();
        assert!(c.push_block(&[-2.0, f64::NAN], &mut || false).is_err());
        assert_eq!(c, before);
        let mut outside = counter();
        outside.push_block(&[-20.0, 20.0], &mut || false).unwrap();
        let before = outside.clone();
        assert!(matches!(outside.finish(), Err(FatigueError::AmplitudeOutsideCurve { .. })));
        assert_eq!(outside, before);
    }

    #[test]
    fn cancellation_during_reduction_and_finish_is_atomic() {
        let mut c = counter();
        c.push_block(&[-2.0, 1.0, -3.0, 5.0, -1.0, 3.0], &mut || false).unwrap();
        let before = c.clone();
        for boundary in 0..8 {
            let mut attempt = before.clone();
            let mut calls = 0;
            let result = attempt.finish_with_cancel(&mut || { calls += 1; calls > boundary });
            if result.is_err() { assert_eq!(attempt, before); }
        }
        let mut calls = 0;
        assert_eq!(c.push_block(&[-4.0, 4.0, -2.0], &mut || { calls += 1; calls > 3 }), Err(FatigueError::Cancelled));
        assert_eq!(c, before);
    }

    #[test]
    fn limits_refuse_without_forgetting_unclosed_history() {
        let mut c = RainflowMiner::new(curve(), MeanStressCorrection::Uncorrected,
            RainflowLimits { maximum_samples: 10, maximum_residue: 2 }).unwrap();
        c.push_block(&[-5.0, 5.0, -4.0], &mut || false).unwrap();
        let before = c.clone();
        assert_eq!(c.push(4.0), Err(FatigueError::ResidueLimit { limit: 2 }));
        assert_eq!(c, before);
        let mut capped = RainflowMiner::new(curve(), MeanStressCorrection::Uncorrected,
            RainflowLimits { maximum_samples: 2, maximum_residue: 8 }).unwrap();
        capped.push_block(&[-2.0, 2.0], &mut || false).unwrap();
        let before = capped.clone();
        assert_eq!(capped.push(0.0), Err(FatigueError::SampleLimit { limit: 2 }));
        assert_eq!(capped, before);
        assert!(capped.finish().unwrap().finished);
    }

    #[test]
    fn fixed_plane_projection_retains_reversal_and_tensor_shear_factor() {
        let axial = NormalStressDirection::new([1.0, 0.0, 0.0]).unwrap();
        assert_eq!(axial.project(&[-7.0, 2.0, 3.0, 0.0, 0.0, 0.0]).unwrap(), -7.0);
        let diagonal = NormalStressDirection::new([1.0, 1.0, 0.0]).unwrap();
        close(diagonal.project(&[0.0, 0.0, 0.0, 4.0, 0.0, 0.0]).unwrap(), 4.0);
        assert!(NormalStressDirection::new([0.0; 3]).is_err());
        assert!(NormalStressDirection::new([f64::NAN, 0.0, 1.0]).is_err());
    }

    #[test]
    fn long_periodic_history_keeps_bounded_residue_and_exact_cycle_counts() {
        let mut c = counter();
        c.push(-2.0).unwrap();
        for _ in 0..2000 { c.push(2.0).unwrap(); c.push(-2.0).unwrap(); assert!(c.pending_reversals() <= 2); }
        let summary = c.finish().unwrap();
        assert_eq!(summary.full_cycles * 2 + summary.half_cycles, 4000);
        close(summary.miner_damage, 4.0);
    }
}
