//! Bounded Coffin-Manson-Basquin total mechanical strain-life (MR62).
//!
//! `epsilon_a = sigma_f'/E * (2 N_f)^b + epsilon_f' * (2 N_f)^c`.
//! `epsilon_a` is half the signed total mechanical strain range; `N_f` is
//! FULL CYCLES, while `2 N_f` counts reversals. The elastic and plastic terms
//! describe one calibrated fatigue curve, not an additional constitutive law.
//! See COMSOL's "Strain Life Models", combined Basquin/Coffin-Manson relation:
//! <https://doc.comsol.com/6.4/doc/com.comsol.help.fatigue/fatigue_ug_sme.4.21.html>.
//!
//! Parameters and scope are caller-supplied calibration declarations. No fit is
//! inferred from yield strength, hardness or a plastic-hardening modulus. The
//! declared finite life interval is strict: neither a low-amplitude endurance
//! limit nor a short-life extrapolation is invented. `Uncorrected` explicitly
//! omits mean-stress effects; Morrow/SWT require the paired stress/strain owner.
//! This uniaxial/fixed-plane empirical Miner estimate is not a multiaxial
//! criterion, stiffness loss, dissipated heat, rupture or service-life proof.

use super::{CycleLife, FatigueError, FatigueScope, RainflowHistory, RainflowLimits};
use fs_math::det;
use std::sync::Arc;

/// Source-backed CMB calibration with explicit full-cycle/reversal fit conventions.
pub mod source;

pub(in crate::fatigue) mod inverse;
/// Paired stress/strain histories and named mean-stress corrections.
mod mean_stress;
pub use mean_stress::{PairedStrainLifeMiner, StrainStressSample};

/// Coefficients of the total mechanical strain-amplitude relation. Both
/// exponents are strictly negative. A zero coefficient admits the corresponding
/// single-term Basquin or Coffin-Manson limit; at least one term must be present.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CmbParameters {
    pub youngs_modulus_pa: f64,
    pub fatigue_strength_coefficient_pa: f64,
    pub fatigue_strength_exponent: f64,
    pub fatigue_ductility_coefficient: f64,
    pub fatigue_ductility_exponent: f64,
}

/// One immutable, explicitly bounded total mechanical strain-life calibration.
#[derive(Debug, Clone, PartialEq)]
pub struct StrainLifeCurve {
    scope: FatigueScope,
    parameters: CmbParameters,
    minimum_cycles: f64,
    maximum_cycles: f64,
    elastic_coefficient: f64,
    log_minimum_reversals: f64,
    log_maximum_reversals: f64,
    minimum_amplitude: f64,
    maximum_amplitude: f64,
    source: Option<Arc<source::CmbCalibrationSource>>,
}

impl StrainLifeCurve {
    /// The supplied domain must span at least two resolvable finite lives and
    /// start at no less than half a cycle (one reversal). Active coefficients
    /// and the strain amplitudes at both endpoints must be representable.
    pub fn new(
        scope: FatigueScope,
        parameters: CmbParameters,
        minimum_cycles: f64,
        maximum_cycles: f64,
    ) -> Result<Self, FatigueError> {
        scope.validate()?;
        let p = parameters;
        if !p.youngs_modulus_pa.is_finite()
            || p.youngs_modulus_pa <= 0.0
            || !p.fatigue_strength_coefficient_pa.is_finite()
            || p.fatigue_strength_coefficient_pa < 0.0
            || !p.fatigue_ductility_coefficient.is_finite()
            || p.fatigue_ductility_coefficient < 0.0
            || (p.fatigue_strength_coefficient_pa == 0.0 && p.fatigue_ductility_coefficient == 0.0)
            || !p.fatigue_strength_exponent.is_finite()
            || p.fatigue_strength_exponent >= 0.0
            || !p.fatigue_ductility_exponent.is_finite()
            || p.fatigue_ductility_exponent >= 0.0
        {
            return Err(FatigueError::InvalidCalibration(
                "CMB needs positive modulus, nonnegative active coefficients and negative exponents",
            ));
        }
        if !minimum_cycles.is_finite()
            || minimum_cycles < 0.5
            || !maximum_cycles.is_finite()
            || maximum_cycles <= minimum_cycles
        {
            return Err(FatigueError::InvalidCalibration(
                "declare a finite strain-life cycle interval starting at or above 0.5 cycles",
            ));
        }
        let elastic_coefficient = p.fatigue_strength_coefficient_pa / p.youngs_modulus_pa;
        if !elastic_coefficient.is_finite()
            || (p.fatigue_strength_coefficient_pa > 0.0 && elastic_coefficient == 0.0)
        {
            return Err(FatigueError::InvalidCalibration(
                "CMB elastic strain coefficient is not representable",
            ));
        }
        // Sum logs instead of forming 2*N: supported N may exceed MAX/2.
        let log_minimum_reversals = core::f64::consts::LN_2 + det::ln(minimum_cycles);
        let log_maximum_reversals = core::f64::consts::LN_2 + det::ln(maximum_cycles);
        let mut curve = Self {
            scope,
            parameters,
            minimum_cycles,
            maximum_cycles,
            elastic_coefficient,
            log_minimum_reversals,
            log_maximum_reversals,
            minimum_amplitude: 0.0,
            maximum_amplitude: 0.0,
            source: None,
        };
        curve.maximum_amplitude = curve.evaluate_log_reversals(log_minimum_reversals).0;
        curve.minimum_amplitude = curve.evaluate_log_reversals(log_maximum_reversals).0;
        if log_minimum_reversals >= log_maximum_reversals
            || !curve.maximum_amplitude.is_finite()
            || !curve.minimum_amplitude.is_finite()
            || curve.minimum_amplitude <= 0.0
            || curve.maximum_amplitude <= curve.minimum_amplitude
        {
            return Err(FatigueError::InvalidCalibration(
                "strain-life endpoints must be positive, finite and numerically distinguishable",
            ));
        }
        Ok(curve)
    }

    pub fn source(&self) -> Option<&source::CmbCalibrationSource> {
        self.source.as_deref()
    }

    pub fn scope(&self) -> &FatigueScope {
        &self.scope
    }
    pub fn parameters(&self) -> CmbParameters {
        self.parameters
    }
    pub fn cycle_domain(&self) -> (f64, f64) {
        (self.minimum_cycles, self.maximum_cycles)
    }
    pub fn amplitude_domain(&self) -> (f64, f64) {
        (self.minimum_amplitude, self.maximum_amplitude)
    }

    /// Total mechanical strain AMPLITUDE at a life in full cycles. Original
    /// endpoint values are retained exactly, avoiding inverse endpoint drift.
    pub fn amplitude_at_cycles(&self, cycles: f64) -> Result<f64, FatigueError> {
        if !cycles.is_finite() || cycles < self.minimum_cycles || cycles > self.maximum_cycles {
            return Err(FatigueError::LifeOutsideCurve {
                cycles,
                minimum: self.minimum_cycles,
                maximum: self.maximum_cycles,
            });
        }
        if cycles == self.minimum_cycles {
            return Ok(self.maximum_amplitude);
        }
        if cycles == self.maximum_cycles {
            return Ok(self.minimum_amplitude);
        }
        let amplitude = self
            .evaluate_log_reversals(core::f64::consts::LN_2 + det::ln(cycles))
            .0;
        if !amplitude.is_finite() || amplitude <= 0.0 {
            return Err(FatigueError::Arithmetic);
        }
        Ok(amplitude)
    }

    /// Invert only within the declared calibration. `None` is reserved for
    /// exactly zero amplitude; positive sub-domain amplitudes are unsupported.
    /// A bracketed Newton solve in log-reversals uses at most 64 iterations and
    /// deterministic elementary functions. Bisection safeguards every step.
    pub fn cycles_to_failure(&self, amplitude: f64) -> Result<Option<f64>, FatigueError> {
        self.invert_life(
            amplitude,
            self.minimum_amplitude,
            self.maximum_amplitude,
            |x| self.evaluate_log_reversals(x),
        )
    }

    // Every correction shares this bounded inverse and the original cycle
    // domain; only its amplitude relation and endpoint support change.
    fn invert_life(
        &self,
        amplitude: f64,
        minimum_amplitude: f64,
        maximum_amplitude: f64,
        evaluate: impl Fn(f64) -> (f64, f64),
    ) -> Result<Option<f64>, FatigueError> {
        inverse::invert_life(
            amplitude,
            (minimum_amplitude, maximum_amplitude),
            (self.minimum_cycles, self.maximum_cycles),
            (self.log_minimum_reversals, self.log_maximum_reversals),
            evaluate,
        )
    }

    fn evaluate_log_reversals(&self, log_reversals: f64) -> (f64, f64) {
        let p = self.parameters;
        let elastic =
            self.elastic_coefficient * det::exp(p.fatigue_strength_exponent * log_reversals);
        let plastic = p.fatigue_ductility_coefficient
            * det::exp(p.fatigue_ductility_exponent * log_reversals);
        (
            elastic + plastic,
            p.fatigue_strength_exponent * elastic + p.fatigue_ductility_exponent * plastic,
        )
    }
}

/// Required mean-stress convention. A scalar strain trace does not supply the
/// stress extrema needed by Morrow or Smith-Watson-Topper corrections.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StrainLifeMeanStress {
    /// Explicitly neglect mean stress; applicability must be calibrated.
    Uncorrected,
    /// Morrow modifies the elastic fatigue-strength term with the cycle's
    /// actual signed mean stress. Compressive mean receives model credit.
    Morrow,
    /// Fixed-plane Smith-Watson-Topper uses actual maximum tensile stress
    /// times total mechanical strain amplitude. This is not a plane search.
    SmithWatsonTopper,
}

/// Accounted cycles and dimensionless Miner damage. Open endpoint half cycles
/// are excluded until `finish`. Damage >= 1 is an empirical exhaustion estimate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StrainLifeSummary {
    pub samples: u64,
    pub full_cycles: u64,
    pub half_cycles: u64,
    pub miner_damage: f64,
    pub maximum_strain_amplitude: f64,
    pub finished: bool,
}

/// Signed total mechanical strain history using the SAME bounded rainflow
/// reduction and compensated damage summation as the stress-life owner.
#[derive(Debug, Clone, PartialEq)]
pub struct StrainLifeMiner {
    curve: Arc<StrainLifeCurve>,
    mean_stress: StrainLifeMeanStress,
    history: RainflowHistory,
}

impl StrainLifeMiner {
    pub fn new(
        curve: Arc<StrainLifeCurve>,
        mean_stress: StrainLifeMeanStress,
        limits: RainflowLimits,
    ) -> Result<Self, FatigueError> {
        match mean_stress {
            StrainLifeMeanStress::Uncorrected => {}
            StrainLifeMeanStress::Morrow | StrainLifeMeanStress::SmithWatsonTopper => {
                return Err(FatigueError::InvalidCalibration(
                    "mean-stress correction requires PairedStrainLifeMiner and actual stress/strain samples",
                ));
            }
        }
        Ok(Self {
            curve,
            mean_stress,
            history: RainflowHistory::new(limits)?,
        })
    }
    pub fn curve(&self) -> &StrainLifeCurve {
        &self.curve
    }
    pub fn mean_stress(&self) -> StrainLifeMeanStress {
        self.mean_stress
    }
    pub fn pending_reversals(&self) -> usize {
        self.history.residue.len()
    }
    pub fn summary(&self) -> StrainLifeSummary {
        StrainLifeSummary {
            samples: self.history.samples,
            full_cycles: self.history.tally.full,
            half_cycles: self.history.tally.half,
            miner_damage: self.history.tally.damage,
            maximum_strain_amplitude: self.history.tally.maximum_amplitude,
            finished: self.history.finished,
        }
    }
    pub fn push(&mut self, mechanical_strain: f64) -> Result<(), FatigueError> {
        self.push_with_cancel(mechanical_strain, &mut || false)
    }
    /// Failed sample or cycle admission leaves all accepted history untouched.
    pub fn push_with_cancel(
        &mut self,
        mechanical_strain: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), FatigueError> {
        self.history
            .push(mechanical_strain, self.curve.as_ref(), cancelled)
    }
    /// Blocks are atomic; never finalize at callback or solver-window boundaries.
    pub fn push_block(
        &mut self,
        mechanical_strains: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), FatigueError> {
        self.history
            .push_block(mechanical_strains, self.curve.as_ref(), cancelled)
    }
    pub fn finish(&mut self) -> Result<StrainLifeSummary, FatigueError> {
        self.finish_with_cancel(&mut || false)
    }
    /// Flush unresolved endpoint half cycles once, atomically and idempotently.
    pub fn finish_with_cancel(
        &mut self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<StrainLifeSummary, FatigueError> {
        self.history.finish(self.curve.as_ref(), cancelled)?;
        Ok(self.summary())
    }
}

impl CycleLife for StrainLifeCurve {
    fn admit_sample(&self, value: f64, _context: (), sample: u64) -> Result<(), FatigueError> {
        if value.is_finite() {
            Ok(())
        } else {
            Err(FatigueError::InvalidStrain { sample })
        }
    }
    fn evaluate_cycle(
        &self,
        amplitude: f64,
        _mean: f64,
        _contexts: [(); 2],
    ) -> Result<(Option<f64>, f64), FatigueError> {
        Ok((self.cycles_to_failure(amplitude)?, amplitude))
    }
}

#[cfg(test)]
mod tests {
    use super::super::{
        MeanStressCorrection, NormalStressDirection, RainflowMiner, SnCurve, SnPoint,
    };
    use super::*;

    fn scope() -> FatigueScope {
        FatigueScope {
            source_reference: "synthetic analytical CMB; no measured-life claim".into(),
            material: "fixture".into(),
            process: "fixture".into(),
            surface_and_geometry: "uniaxial smooth coupon".into(),
            environment: "isothermal".into(),
            load_basis: "signed axial total mechanical strain amplitude; full cycles".into(),
            temperature_k: 293.15,
        }
    }
    fn parameters() -> CmbParameters {
        CmbParameters {
            youngs_modulus_pa: 200e9,
            fatigue_strength_coefficient_pa: 1e9,
            fatigue_strength_exponent: -0.1,
            fatigue_ductility_coefficient: 0.5,
            fatigue_ductility_exponent: -0.6,
        }
    }
    fn curve() -> Arc<StrainLifeCurve> {
        Arc::new(StrainLifeCurve::new(scope(), parameters(), 0.5, 1e9).unwrap())
    }
    fn counter(curve: Arc<StrainLifeCurve>) -> StrainLifeMiner {
        StrainLifeMiner::new(
            curve,
            StrainLifeMeanStress::Uncorrected,
            RainflowLimits {
                maximum_samples: 10000,
                maximum_residue: 64,
            },
        )
        .unwrap()
    }
    fn close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 3e-12 * expected.abs(),
            "{actual} != {expected}"
        );
    }

    #[test]
    fn g1_cmb_roundtrips_full_cycles_and_each_physical_term() {
        let original = parameters();
        for p in [
            original,
            CmbParameters {
                fatigue_ductility_coefficient: 0.0,
                ..original
            },
            CmbParameters {
                fatigue_strength_coefficient_pa: 0.0,
                ..original
            },
        ] {
            let c = StrainLifeCurve::new(scope(), p, 0.5, 1e9).unwrap();
            for n in [0.5, 1.0, 2.0, 10.0, 100.0, 1000.0, 1e5, 1e8, 1e9] {
                // Independent platform-power arithmetic is a test oracle only.
                let expected = p.fatigue_strength_coefficient_pa / p.youngs_modulus_pa
                    * (2.0_f64 * n).powf(p.fatigue_strength_exponent)
                    + p.fatigue_ductility_coefficient
                        * (2.0_f64 * n).powf(p.fatigue_ductility_exponent);
                let amplitude = c.amplitude_at_cycles(n).unwrap();
                close(amplitude, expected);
                close(c.cycles_to_failure(amplitude).unwrap().unwrap(), n);
            }
        }
    }

    #[test]
    fn g1_equal_exponents_have_an_independent_closed_form_life() {
        let p = CmbParameters {
            fatigue_strength_exponent: -0.5,
            fatigue_ductility_exponent: -0.5,
            ..parameters()
        };
        let c = StrainLifeCurve::new(scope(), p, 0.5, 1e8).unwrap();
        // epsilon_a = (.005 + .5) / sqrt(2N), hence N = .5*(.505/a)^2.
        for amplitude in [0.01, 0.02, 0.1, 0.3] {
            close(
                c.cycles_to_failure(amplitude).unwrap().unwrap(),
                0.5 * (0.505 / amplitude).powi(2),
            );
        }
        let mut flow = counter(Arc::new(c));
        flow.push_block(&[-0.01, 0.01, -0.01], &mut || false)
            .unwrap();
        let result = flow.finish().unwrap();
        assert_eq!((result.full_cycles, result.half_cycles), (0, 2));
        assert_eq!(result.maximum_strain_amplitude, 0.01);
        close(
            result.miner_damage,
            1.0 / (0.5 * (0.505_f64 / 0.01).powi(2)),
        );
        assert_eq!(flow.finish().unwrap(), result);
        assert_eq!(flow.push(0.0), Err(FatigueError::AlreadyFinished));
    }

    #[test]
    fn g0_support_is_strict_and_zero_does_not_imply_a_fatigue_limit() {
        let c = curve();
        let (small, large) = c.amplitude_domain();
        assert_eq!(c.cycles_to_failure(0.0).unwrap(), None);
        assert_eq!(c.cycles_to_failure(small).unwrap(), Some(1e9));
        assert_eq!(c.cycles_to_failure(large).unwrap(), Some(0.5));
        for value in [small.next_down(), large.next_up()] {
            assert!(matches!(
                c.cycles_to_failure(value),
                Err(FatigueError::StrainAmplitudeOutsideCurve { .. })
            ));
        }
        for n in [
            0.0,
            0.5_f64.next_down(),
            1e9_f64.next_up(),
            f64::INFINITY,
            f64::NAN,
        ] {
            assert!(matches!(
                c.amplitude_at_cycles(n),
                Err(FatigueError::LifeOutsideCurve { .. })
            ));
        }
        let mut flow = counter(c);
        flow.push_block(&[-small / 4.0, small / 4.0], &mut || false)
            .unwrap();
        let before = flow.clone();
        assert!(matches!(
            flow.finish(),
            Err(FatigueError::StrainAmplitudeOutsideCurve { .. })
        ));
        assert_eq!(flow, before);
    }

    #[test]
    fn g3_strain_and_stress_models_use_identical_nested_cycle_reduction() {
        let p = CmbParameters {
            fatigue_strength_coefficient_pa: 0.0,
            fatigue_ductility_coefficient: 2.0,
            fatigue_ductility_exponent: -1.0,
            ..parameters()
        };
        // epsilon_a = 1/N. The S-N table is given the same numeric relation.
        let c = Arc::new(StrainLifeCurve::new(scope(), p, 10.0, 100000.0).unwrap());
        let sn = Arc::new(
            SnCurve::new(
                scope(),
                vec![
                    SnPoint {
                        amplitude_pa: 0.00001,
                        cycles_to_failure: 100000.0,
                    },
                    SnPoint {
                        amplitude_pa: 0.1,
                        cycles_to_failure: 10.0,
                    },
                ],
                None,
            )
            .unwrap(),
        );
        let samples = [
            -0.002, 0.001, -0.003, 0.005, -0.001, 0.003, -0.004, 0.004, -0.002,
        ];
        let mut whole = counter(c.clone());
        let mut stress = RainflowMiner::new(
            sn,
            MeanStressCorrection::Uncorrected,
            RainflowLimits {
                maximum_samples: 10000,
                maximum_residue: 64,
            },
        )
        .unwrap();
        whole.push_block(&samples, &mut || false).unwrap();
        stress.push_block(&samples, &mut || false).unwrap();
        let expected = whole.finish().unwrap();
        let stress_result = stress.finish().unwrap();
        assert_eq!((expected.full_cycles, expected.half_cycles), (1, 6));
        assert_eq!(
            (expected.full_cycles, expected.half_cycles),
            (stress_result.full_cycles, stress_result.half_cycles)
        );
        close(expected.miner_damage, 0.0115);
        close(expected.miner_damage, stress_result.miner_damage);
        for split in 0..=samples.len() {
            let mut chunked = counter(c.clone());
            chunked
                .push_block(&samples[..split], &mut || false)
                .unwrap();
            chunked
                .push_block(&samples[split..], &mut || false)
                .unwrap();
            assert_eq!(chunked.finish().unwrap(), expected);
        }
    }

    #[test]
    fn g0_plateaus_offsets_and_directional_tensor_strain_are_explicit() {
        let mut zero = counter(curve());
        zero.push_block(&[0.01; 12], &mut || false).unwrap();
        assert_eq!(zero.finish().unwrap().miner_damage, 0.0);
        let mut reversed = counter(curve());
        reversed
            .push_block(&[-0.01, -0.01, 0.01, 0.01, -0.01], &mut || false)
            .unwrap();
        let mut offset = counter(curve());
        offset
            .push_block(&[0.0, 0.0, 0.02, 0.02, 0.0], &mut || false)
            .unwrap();
        assert_eq!(reversed.finish().unwrap(), offset.finish().unwrap());
        let normal = NormalStressDirection::new([1.0, 1.0, 0.0]).unwrap();
        close(
            normal
                .project_strain(&[0.0, 0.0, 0.0, -0.004, 0.0, 0.0])
                .unwrap(),
            -0.004,
        );
    }

    #[test]
    fn g4_strain_refusals_and_cancellation_do_not_publish_partial_history() {
        let mut flow = counter(curve());
        flow.push_block(&[-0.01, 0.01], &mut || false).unwrap();
        let before = flow.clone();
        assert!(matches!(
            flow.push_block(&[-0.01, f64::NAN], &mut || false),
            Err(FatigueError::InvalidStrain { .. })
        ));
        assert_eq!(flow, before);
        for boundary in 0..6 {
            let mut attempt = before.clone();
            let mut calls = 0;
            let result = attempt.finish_with_cancel(&mut || {
                calls += 1;
                calls > boundary
            });
            if result.is_err() {
                assert_eq!(attempt, before);
            }
        }
        let mut calls = 0;
        assert_eq!(
            flow.push_block(&[-0.01, 0.01, -0.01], &mut || {
                calls += 1;
                calls > 3
            }),
            Err(FatigueError::Cancelled)
        );
        assert_eq!(flow, before);
    }

    #[test]
    fn g0_invalid_cmb_parameters_and_degenerate_support_refuse() {
        let p = parameters();
        for bad in [
            CmbParameters {
                youngs_modulus_pa: 0.0,
                ..p
            },
            CmbParameters {
                fatigue_strength_coefficient_pa: -1.0,
                ..p
            },
            CmbParameters {
                fatigue_ductility_coefficient: f64::NAN,
                ..p
            },
            CmbParameters {
                fatigue_strength_exponent: 0.0,
                ..p
            },
            CmbParameters {
                fatigue_ductility_exponent: 0.1,
                ..p
            },
            CmbParameters {
                fatigue_strength_coefficient_pa: 0.0,
                fatigue_ductility_coefficient: 0.0,
                ..p
            },
        ] {
            assert!(StrainLifeCurve::new(scope(), bad, 1.0, 1e6).is_err());
        }
        for (lo, hi) in [
            (0.0, 1e6),
            (0.49, 1e6),
            (1.0, 1.0),
            (2.0, 1.0),
            (1.0, f64::INFINITY),
        ] {
            assert!(StrainLifeCurve::new(scope(), p, lo, hi).is_err());
        }
    }

    #[test]
    fn g0_log_reversals_do_not_overflow_twice_the_supported_life() {
        let c = StrainLifeCurve::new(scope(), parameters(), 0.5, f64::MAX).unwrap();
        for cycles in [1e100, 1e200, f64::MAX] {
            let amplitude = c.amplitude_at_cycles(cycles).unwrap();
            close(c.cycles_to_failure(amplitude).unwrap().unwrap(), cycles);
        }
        let p = CmbParameters {
            fatigue_strength_coefficient_pa: f64::MIN_POSITIVE,
            youngs_modulus_pa: f64::MAX,
            ..parameters()
        };
        assert!(StrainLifeCurve::new(scope(), p, 0.5, 1e6).is_err());
    }
}
