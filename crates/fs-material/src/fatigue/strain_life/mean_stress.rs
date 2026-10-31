//! Named mean-stress corrections on coincident signed stress/strain reversals.
//!
//! Morrow: `epsilon_a = (sigma_f' - sigma_mean)/E * (2N)^b
//!                       + epsilon_f' * (2N)^c`.
//! SWT: `sigma_max * epsilon_a = sigma_f'^2/E * (2N)^(2b)
//!                             + sigma_f' * epsilon_f' * (2N)^(b+c)`.
//! These are empirical alternatives, not equivalent constitutive laws. See
//! COMSOL "Strain Life Models" / "Strain Based Fatigue Models", and NASA
//! CR-189670, section 3.2 (mean-stress correction of counted strain cycles).
//! <https://doc.comsol.com/6.4/doc/com.comsol.help.fatigue/fatigue_ug_sme.4.23.html>
//! <https://ntrs.nasa.gov/citations/19930003223>
//!
//! The paired counter admits a scalar fixed-plane history whose observed
//! stress and strain extrema coincide. Every observed increment must have
//! matching stress/strain direction; stress evolution at unchanged strain is
//! refused. A nested strain-cycle reduction may not discard stress extrema
//! outside the retained outer endpoint envelope. Accepted samples must resolve
//! reversals: these checks do not certify unsampled extrema or a multiaxial
//! loading path. Cycle stress extrema are the stress values paired with its
//! counted strain endpoints, not stresses reconstructed from an elastic modulus.
//! The constant-cycle correction applied to evolving/unstabilized cycles is an
//! explicit empirical approximation, not a hardening/relaxation or sequence law.

use super::super::{CycleLife, FatigueError, RainflowHistory, RainflowLimits, TurningPoint};
use super::{StrainLifeCurve, StrainLifeMeanStress, StrainLifeSummary};
use fs_math::det;
use std::sync::Arc;

/// Simultaneous signed normal components on one fixed physical plane.
/// Mechanical strain includes elastic and plastic parts, excludes thermal
/// eigenstrain, and is not accumulated equivalent plastic strain.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StrainStressSample {
    pub mechanical_strain: f64,
    pub normal_stress_pa: f64,
}
impl StrainStressSample {
    fn point(self) -> TurningPoint<f64> {
        TurningPoint {
            value: self.mechanical_strain,
            context: self.normal_stress_pa,
        }
    }
}

impl StrainLifeCurve {
    /// Evaluate one cycle using its actual stress extrema and total mechanical
    /// strain AMPLITUDE. Life remains inside the original calibrated cycle
    /// interval; the amplitude support is reevaluated for this cycle's mean/max.
    ///
    /// Morrow's signed compressive mean receives the model's life credit. A
    /// nonpositive elastic fatigue coefficient is refused, not called static
    /// failure: `sigma_f'` is not ultimate tensile strength. SWT refuses nonzero
    /// strain cycles without tensile maximum stress rather than reporting an
    /// invented infinite life for compression-only fatigue.
    pub fn cycles_to_failure_with_stress(
        &self,
        strain_amplitude: f64,
        stress_min_pa: f64,
        stress_max_pa: f64,
        correction: StrainLifeMeanStress,
    ) -> Result<Option<f64>, FatigueError> {
        validate_policy(self, correction)?;
        if !stress_min_pa.is_finite() || !stress_max_pa.is_finite() || stress_min_pa > stress_max_pa
        {
            return Err(FatigueError::InvalidStressExtrema {
                minimum_pa: stress_min_pa,
                maximum_pa: stress_max_pa,
            });
        }
        if !strain_amplitude.is_finite() || strain_amplitude < 0.0 {
            return Err(FatigueError::Arithmetic);
        }
        if strain_amplitude == 0.0 {
            return Ok(None);
        }
        let p = self.parameters;
        match correction {
            StrainLifeMeanStress::Uncorrected => self.cycles_to_failure(strain_amplitude),
            StrainLifeMeanStress::Morrow => {
                let mean = 0.5 * stress_min_pa + 0.5 * stress_max_pa;
                if mean >= p.fatigue_strength_coefficient_pa {
                    return Err(FatigueError::MeanStressOutsideStrainLife {
                        mean_pa: mean,
                        fatigue_strength_coefficient_pa: p.fatigue_strength_coefficient_pa,
                    });
                }
                let elastic_coefficient = self.elastic_coefficient - mean / p.youngs_modulus_pa;
                if !elastic_coefficient.is_finite() || elastic_coefficient <= 0.0 {
                    return Err(FatigueError::Arithmetic);
                }
                let evaluate = |x| {
                    let elastic = elastic_coefficient * det::exp(p.fatigue_strength_exponent * x);
                    let plastic = p.fatigue_ductility_coefficient
                        * det::exp(p.fatigue_ductility_exponent * x);
                    (
                        elastic + plastic,
                        p.fatigue_strength_exponent * elastic
                            + p.fatigue_ductility_exponent * plastic,
                    )
                };
                self.invert_life(
                    strain_amplitude,
                    evaluate(self.log_maximum_reversals).0,
                    evaluate(self.log_minimum_reversals).0,
                    evaluate,
                )
            }
            StrainLifeMeanStress::SmithWatsonTopper => {
                if stress_max_pa <= 0.0 {
                    return Err(FatigueError::NoTensileStress {
                        maximum_pa: stress_max_pa,
                    });
                }
                // Divide by sigma_f' first: avoid ever squaring a stress, and
                // return support in actual strain-amplitude units, not Pa.
                let stress_ratio = stress_max_pa / p.fatigue_strength_coefficient_pa;
                if !stress_ratio.is_finite() || stress_ratio <= 0.0 {
                    return Err(FatigueError::Arithmetic);
                }
                let b = 2.0 * p.fatigue_strength_exponent;
                let c = p.fatigue_strength_exponent + p.fatigue_ductility_exponent;
                let evaluate = |x| {
                    let elastic = self.elastic_coefficient * det::exp(b * x);
                    let plastic = p.fatigue_ductility_coefficient * det::exp(c * x);
                    (
                        (elastic + plastic) / stress_ratio,
                        (b * elastic + c * plastic) / stress_ratio,
                    )
                };
                self.invert_life(
                    strain_amplitude,
                    evaluate(self.log_maximum_reversals).0,
                    evaluate(self.log_minimum_reversals).0,
                    evaluate,
                )
            }
        }
    }
}

fn validate_policy(
    curve: &StrainLifeCurve,
    correction: StrainLifeMeanStress,
) -> Result<(), FatigueError> {
    match correction {
        StrainLifeMeanStress::Uncorrected => Ok(()),
        StrainLifeMeanStress::Morrow | StrainLifeMeanStress::SmithWatsonTopper => {
            let p = curve.parameters;
            if p.fatigue_strength_coefficient_pa <= 0.0
                || !(2.0 * p.fatigue_strength_exponent).is_finite()
                || !(p.fatigue_strength_exponent + p.fatigue_ductility_exponent).is_finite()
            {
                Err(FatigueError::InvalidCalibration(
                    "mean-stress correction needs a positive fatigue-strength coefficient and representable exponents",
                ))
            } else {
                Ok(())
            }
        }
    }
}

struct PairedLife<'a> {
    curve: &'a StrainLifeCurve,
    correction: StrainLifeMeanStress,
}
impl CycleLife<f64> for PairedLife<'_> {
    fn admit_sample(&self, value: f64, stress: f64, sample: u64) -> Result<(), FatigueError> {
        if !value.is_finite() {
            return Err(FatigueError::InvalidStrain { sample });
        }
        if !stress.is_finite() {
            return Err(FatigueError::InvalidStress { sample });
        }
        Ok(())
    }
    fn admit_transition(
        &self,
        previous: TurningPoint<f64>,
        next: TurningPoint<f64>,
        sample: u64,
    ) -> Result<(), FatigueError> {
        let invalid = if next.value > previous.value {
            next.context < previous.context
        } else if next.value < previous.value {
            next.context > previous.context
        } else {
            next.context != previous.context
        };
        if invalid {
            Err(FatigueError::UnsupportedPairedHistory {
                sample,
                what: "stress/strain extrema must coincide; phase lag, relaxation and opposite increments are unsupported",
            })
        } else {
            Ok(())
        }
    }
    fn admit_reduction(
        &self,
        outer_start: Option<TurningPoint<f64>>,
        a: TurningPoint<f64>,
        b: TurningPoint<f64>,
        end: TurningPoint<f64>,
        sample: u64,
    ) -> Result<(), FatigueError> {
        if let Some(start) = outer_start {
            self.admit_transition(start, end, sample)?;
            let low = start.context.min(end.context);
            let high = start.context.max(end.context);
            if [a.context, b.context]
                .into_iter()
                .any(|stress| stress < low || stress > high)
            {
                return Err(FatigueError::UnsupportedPairedHistory {
                    sample,
                    what: "nested strain reduction would hide an unmatched stress extreme",
                });
            }
        } else {
            self.admit_transition(b, end, sample)?;
        }
        Ok(())
    }
    fn evaluate_cycle(
        &self,
        amplitude: f64,
        _strain_mean: f64,
        stress: [f64; 2],
    ) -> Result<(Option<f64>, f64), FatigueError> {
        Ok((
            self.curve.cycles_to_failure_with_stress(
                amplitude,
                stress[0].min(stress[1]),
                stress[0].max(stress[1]),
                self.correction,
            )?,
            amplitude,
        ))
    }
}

/// Actual paired stress/strain observations with one bounded shared rainflow
/// owner. Failure and cancellation preserve the accepted pairs and all damage.
#[derive(Debug, Clone, PartialEq)]
pub struct PairedStrainLifeMiner {
    curve: Arc<StrainLifeCurve>,
    correction: StrainLifeMeanStress,
    history: RainflowHistory<f64>,
}
impl PairedStrainLifeMiner {
    pub fn new(
        curve: Arc<StrainLifeCurve>,
        correction: StrainLifeMeanStress,
        limits: RainflowLimits,
    ) -> Result<Self, FatigueError> {
        validate_policy(&curve, correction)?;
        Ok(Self {
            curve,
            correction,
            history: RainflowHistory::new(limits)?,
        })
    }
    pub fn curve(&self) -> &StrainLifeCurve {
        &self.curve
    }
    pub fn mean_stress(&self) -> StrainLifeMeanStress {
        self.correction
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
    pub fn push(&mut self, sample: StrainStressSample) -> Result<(), FatigueError> {
        self.push_with_cancel(sample, &mut || false)
    }
    pub fn push_with_cancel(
        &mut self,
        sample: StrainStressSample,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), FatigueError> {
        let law = PairedLife {
            curve: &self.curve,
            correction: self.correction,
        };
        self.history.push_point(sample.point(), &law, cancelled)
    }
    pub fn push_block(
        &mut self,
        samples: &[StrainStressSample],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), FatigueError> {
        let law = PairedLife {
            curve: &self.curve,
            correction: self.correction,
        };
        self.history
            .push_points(samples.iter().map(|sample| sample.point()), &law, cancelled)
    }
    pub fn finish(&mut self) -> Result<StrainLifeSummary, FatigueError> {
        self.finish_with_cancel(&mut || false)
    }
    pub fn finish_with_cancel(
        &mut self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<StrainLifeSummary, FatigueError> {
        let law = PairedLife {
            curve: &self.curve,
            correction: self.correction,
        };
        self.history.finish(&law, cancelled)?;
        Ok(self.summary())
    }
}

#[cfg(test)]
mod tests {
    use super::super::super::FatigueScope;
    use super::super::{CmbParameters, StrainLifeMiner};
    use super::*;

    fn curve() -> Arc<StrainLifeCurve> {
        Arc::new(
            StrainLifeCurve::new(
                FatigueScope {
                    source_reference: "synthetic closed-form references; not measured life".into(),
                    material: "fixture".into(),
                    process: "fixture".into(),
                    surface_and_geometry: "axial coupon".into(),
                    environment: "isothermal".into(),
                    load_basis: "total mechanical normal strain and paired normal stress".into(),
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
    fn counter(policy: StrainLifeMeanStress) -> PairedStrainLifeMiner {
        PairedStrainLifeMiner::new(
            curve(),
            policy,
            RainflowLimits {
                maximum_samples: 10000,
                maximum_residue: 64,
            },
        )
        .unwrap()
    }
    fn sample(strain: f64, stress: f64) -> StrainStressSample {
        StrainStressSample {
            mechanical_strain: strain,
            normal_stress_pa: stress,
        }
    }
    fn close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 4e-12 * expected.abs().max(1e-20),
            "{actual} != {expected}"
        );
    }

    #[test]
    fn g1_morrow_uses_the_actual_signed_stress_mean() {
        let c = curve();
        let amplitude = 0.001;
        let mut lives = Vec::new();
        for mean in [-200e6, 0.0, 200e6] {
            let life = c
                .cycles_to_failure_with_stress(
                    amplitude,
                    mean - 100e6,
                    mean + 100e6,
                    StrainLifeMeanStress::Morrow,
                )
                .unwrap()
                .unwrap();
            let expected = 0.5 * (((1e9 - mean) / 200e9 + 0.002) / amplitude).powi(2);
            close(life, expected);
            lives.push(life);
        }
        assert!(lives[0] > lives[1] && lives[1] > lives[2]);
        close(lives[1], c.cycles_to_failure(amplitude).unwrap().unwrap());
    }

    #[test]
    fn g1_swt_product_has_an_independent_closed_form_and_reversed_limit() {
        let c = curve();
        for maximum in [100e6, 300e6, 600e6] {
            let actual = c
                .cycles_to_failure_with_stress(
                    0.001,
                    maximum - 200e6,
                    maximum,
                    StrainLifeMeanStress::SmithWatsonTopper,
                )
                .unwrap()
                .unwrap();
            // b=c=-1/2: sigma_max*epsilon_a = (sigma_f'^2/E+sigma_f'*epsilon_f')/(2N).
            close(
                actual,
                (1e18 / 200e9 + 1e9 * 0.002) / (2.0 * maximum * 0.001),
            );
        }
        for n in [1.0, 10.0, 1000.0] {
            let amplitude = c.amplitude_at_cycles(n).unwrap();
            let stress = 1e9 / (2.0 * n).sqrt();
            close(
                c.cycles_to_failure_with_stress(
                    amplitude,
                    -stress,
                    stress,
                    StrainLifeMeanStress::SmithWatsonTopper,
                )
                .unwrap()
                .unwrap(),
                n,
            );
        }
    }

    #[test]
    fn g0_corrections_require_stress_and_admit_zero_but_not_compressive_swt_cycles() {
        let c = curve();
        for policy in [
            StrainLifeMeanStress::Morrow,
            StrainLifeMeanStress::SmithWatsonTopper,
        ] {
            assert!(
                StrainLifeMiner::new(
                    c.clone(),
                    policy,
                    RainflowLimits {
                        maximum_samples: 10,
                        maximum_residue: 8
                    }
                )
                .is_err()
            );
            let mut flow = counter(policy);
            flow.push(sample(0.0, 0.0)).unwrap();
            assert_eq!(flow.finish().unwrap().miner_damage, 0.0);
        }
        assert!(matches!(
            c.cycles_to_failure_with_stress(
                0.001,
                -200e6,
                0.0,
                StrainLifeMeanStress::SmithWatsonTopper
            ),
            Err(FatigueError::NoTensileStress { .. })
        ));
        assert!(matches!(
            c.cycles_to_failure_with_stress(0.001, 0.9e9, 1.1e9, StrainLifeMeanStress::Morrow),
            Err(FatigueError::MeanStressOutsideStrainLife { .. })
        ));
        assert!(matches!(
            c.cycles_to_failure_with_stress(0.001, 10.0, -10.0, StrainLifeMeanStress::Morrow),
            Err(FatigueError::InvalidStressExtrema { .. })
        ));
        // The mean-shifted support is recomputed at the original life bounds.
        assert!(matches!(
            c.cycles_to_failure_with_stress(
                c.amplitude_domain().1,
                0.8e9,
                1e9,
                StrainLifeMeanStress::Morrow
            ),
            Err(FatigueError::StrainAmplitudeOutsideCurve { .. })
        ));
    }

    #[test]
    fn g1_counted_mean_stress_changes_damage_without_changing_strain_cycles() {
        for policy in [
            StrainLifeMeanStress::Morrow,
            StrainLifeMeanStress::SmithWatsonTopper,
        ] {
            let mut low = counter(policy);
            let mut high = counter(policy);
            for strain in [-0.001, 0.001, -0.001] {
                low.push(sample(strain, 100e9 * strain)).unwrap();
                high.push(sample(strain, 100e9 * strain + 200e6)).unwrap();
            }
            let a = low.finish().unwrap();
            let b = high.finish().unwrap();
            assert_eq!(
                (a.full_cycles, a.half_cycles),
                (b.full_cycles, b.half_cycles)
            );
            assert_eq!(a.maximum_strain_amplitude, b.maximum_strain_amplitude);
            assert!(b.miner_damage > a.miner_damage);
            close(
                b.miner_damage,
                1.0 / curve()
                    .cycles_to_failure_with_stress(0.001, 100e6, 300e6, policy)
                    .unwrap()
                    .unwrap(),
            );
        }
    }

    #[test]
    fn g3_nested_paired_cycles_match_hand_counted_endpoints_across_every_block_split() {
        let strains = [-2.0, 1.0, -3.0, 5.0, -1.0, 3.0, -4.0, 4.0, -2.0];
        let samples = strains.map(|e| sample(e * 0.001, 100e6 * e + 300e6));
        let half_pairs = [
            (-2.0, 1.0),
            (1.0, -3.0),
            (-3.0, 5.0),
            (5.0, -4.0),
            (-4.0, 4.0),
            (4.0, -2.0),
        ];
        for policy in [
            StrainLifeMeanStress::Morrow,
            StrainLifeMeanStress::SmithWatsonTopper,
        ] {
            let c = curve();
            let damage = |a: f64, b: f64| {
                1.0 / c
                    .cycles_to_failure_with_stress(
                        0.0005 * (a - b).abs(),
                        100e6 * a.min(b) + 300e6,
                        100e6 * a.max(b) + 300e6,
                        policy,
                    )
                    .unwrap()
                    .unwrap()
            };
            let expected = damage(-1.0, 3.0)
                + half_pairs
                    .into_iter()
                    .map(|(a, b)| 0.5 * damage(a, b))
                    .sum::<f64>();
            let mut whole = counter(policy);
            whole.push_block(&samples, &mut || false).unwrap();
            let summary = whole.finish().unwrap();
            assert_eq!((summary.full_cycles, summary.half_cycles), (1, 6));
            close(summary.miner_damage, expected);
            for split in 0..=samples.len() {
                let mut chunked = counter(policy);
                chunked
                    .push_block(&samples[..split], &mut || false)
                    .unwrap();
                chunked
                    .push_block(&samples[split..], &mut || false)
                    .unwrap();
                assert_eq!(chunked.finish().unwrap(), summary);
            }
        }
    }

    #[test]
    fn g0_lag_relaxation_and_hidden_nested_stress_extrema_refuse_transactionally() {
        let mut flow = counter(StrainLifeMeanStress::Morrow);
        flow.push_block(&[sample(0.0, 0.0), sample(0.001, 100e6)], &mut || false)
            .unwrap();
        let before = flow.clone();
        for bad in [
            sample(0.001, 80e6),
            sample(0.002, 80e6),
            sample(0.0, 120e6),
            sample(f64::NAN, 0.0),
            sample(0.002, f64::NAN),
        ] {
            assert!(flow.push(bad).is_err());
            assert_eq!(flow, before);
        }
        let mut nested = counter(StrainLifeMeanStress::Morrow);
        nested
            .push_block(
                &[
                    sample(-0.005, -50e6),
                    sample(0.004, 100e6),
                    sample(0.001, 10e6),
                    sample(0.006, 60e6),
                ],
                &mut || false,
            )
            .unwrap();
        let before = nested.clone();
        assert!(matches!(
            nested.push(sample(-0.006, -60e6)),
            Err(FatigueError::UnsupportedPairedHistory { .. })
        ));
        assert_eq!(nested, before);
    }

    #[test]
    fn g4_paired_finish_and_blocks_keep_all_state_on_cancellation_or_model_refusal() {
        let mut flow = counter(StrainLifeMeanStress::Morrow);
        flow.push_block(&[sample(-0.001, -100e6), sample(0.001, 100e6)], &mut || {
            false
        })
        .unwrap();
        let before = flow.clone();
        for boundary in 0..8 {
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
            flow.push_block(&[sample(-0.001, -100e6), sample(0.001, 100e6)], &mut || {
                calls += 1;
                calls > 3
            }),
            Err(FatigueError::Cancelled)
        );
        assert_eq!(flow, before);
        let mut compressive = counter(StrainLifeMeanStress::SmithWatsonTopper);
        compressive
            .push_block(
                &[sample(-0.001, -200e6), sample(0.001, -100e6)],
                &mut || false,
            )
            .unwrap();
        let before = compressive.clone();
        assert!(matches!(
            compressive.finish(),
            Err(FatigueError::NoTensileStress { .. })
        ));
        assert_eq!(compressive, before);
    }
}
