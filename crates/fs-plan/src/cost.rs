//! Per-operator cost models (plan §11.4, Bet 12): quantile predictions of
//! wall cost from observed `tune`-table history — machine-specific ON
//! PURPOSE (kernel-level winners flip between the reference machines).
//!
//! Model: log-log power-law fit `cost ≈ exp(a)·size^b` by least squares,
//! with EMPIRICAL residual quantiles supplying the prediction bands — an
//! estimate is itself evidenced (P10/P50/P90 plus observation count and an
//! extrapolation flag; consumers see uncertainty, not a bare number).
//! Insufficient data is a structured refusal, never a guess.
//!
//! Determinism: fits are pure functions of the observation multiset
//! (sorted internally; nearest-rank quantiles with deterministic
//! tie-breaking) — identical ledger snapshots give identical models (P2).

/// One observed execution.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CostObservation {
    /// Problem-size feature (elements, DOF, rays…; > 0).
    pub size: f64,
    /// Measured wall cost in seconds (> 0).
    pub cost_s: f64,
}

/// A quantile cost prediction (an evidenced estimate).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CostPrediction {
    /// 10th-percentile cost, seconds.
    pub p10: f64,
    /// Median cost, seconds.
    pub p50: f64,
    /// 90th-percentile cost, seconds.
    pub p90: f64,
    /// Observations behind the fit.
    pub n_obs: usize,
    /// True when `size` lies outside the observed size range — the bands
    /// are then model extrapolation, not interpolation.
    pub extrapolated: bool,
}

/// Why a prediction was refused (Decalogue P10).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CostRefusal {
    /// Not enough observations for a quantile fit.
    InsufficientData {
        /// Observations available.
        have: usize,
        /// Observations required.
        need: usize,
    },
    /// A nonpositive size/cost was supplied (log-log domain).
    BadInput,
    /// The bounded observation budget is exhausted.
    ObservationLimit {
        /// Maximum observations retained by one model.
        limit: usize,
    },
    /// A freshness policy declared a non-positive horizon
    /// (bead jle3m).
    InvalidFreshnessHorizon {
        /// The refused horizon in ledger nanoseconds.
        horizon_ns: i64,
    },
    /// A fit or prediction produced a nonfinite intermediate/result.
    ArithmeticFailure {
        /// Failing numerical stage.
        stage: &'static str,
    },
    /// A calibration/error audit with no probes has no evidentiary meaning.
    EmptyEvaluation {
        /// Refused evaluation metric.
        metric: &'static str,
    },
    /// A held-out evaluation exceeded its explicit work budget.
    EvaluationLimit {
        /// Number of supplied probes.
        provided: usize,
        /// Maximum admitted probes.
        limit: usize,
    },
}

impl core::fmt::Display for CostRefusal {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            CostRefusal::InsufficientData { have, need } => write!(
                f,
                "cost model refuses to predict from {have} observation(s) (need ≥ {need}); \
                 run the operator or seed the tune table first"
            ),
            CostRefusal::BadInput => {
                write!(
                    f,
                    "sizes and costs must be positive finite (log-log model domain)"
                )
            }
            CostRefusal::ObservationLimit { limit } => {
                write!(f, "cost model observation limit {limit} is exhausted")
            }
            CostRefusal::InvalidFreshnessHorizon { horizon_ns } => write!(
                f,
                "freshness horizon must be positive; got {horizon_ns} ns                  (a zero horizon stales every receipt at its own recording instant)"
            ),
            CostRefusal::ArithmeticFailure { stage } => {
                write!(f, "cost model refused nonfinite arithmetic during {stage}")
            }
            CostRefusal::EmptyEvaluation { metric } => {
                write!(f, "cost model refuses empty {metric} evaluation")
            }
            CostRefusal::EvaluationLimit { provided, limit } => write!(
                f,
                "cost model evaluation supplied {provided} probes, exceeding limit {limit}"
            ),
        }
    }
}

impl core::error::Error for CostRefusal {}

/// Minimum observations before predictions are offered.
pub const MIN_OBS: usize = 3;

/// Maximum observations retained by one online model.
pub const MAX_COST_OBSERVATIONS: usize = 4_096;

/// Maximum held-out observations processed by one evaluation call.
pub const MAX_COST_EVALUATION_OBSERVATIONS: usize = 4_096;

struct DerivedFit {
    loglog: Option<(f64, f64)>,
    residuals: Vec<f64>,
}

/// A fitted per-(operator × shape-class × machine) cost model.
#[derive(Debug, Clone, Default)]
pub struct CostModel {
    obs: Vec<CostObservation>,
    /// (intercept a, slope b) of `ln cost = a + b·ln size`, when fitted.
    loglog: Option<(f64, f64)>,
    /// Sorted residuals `ln cost − (a + b·ln size)`.
    residuals: Vec<f64>,
}

impl CostModel {
    /// An empty model (predictions refuse until observations arrive).
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Fit from a batch of observations.
    ///
    /// # Errors
    /// [`CostRefusal::BadInput`] on nonpositive/nonfinite entries.
    pub fn fit(observations: &[CostObservation]) -> Result<Self, CostRefusal> {
        if observations.len() > MAX_COST_OBSERVATIONS {
            return Err(CostRefusal::ObservationLimit {
                limit: MAX_COST_OBSERVATIONS,
            });
        }
        for observation in observations {
            Self::validate_observation(*observation)?;
        }
        let mut obs = observations.to_vec();
        obs.sort_by(Self::observation_order);
        let DerivedFit { loglog, residuals } = Self::derive_fit(&obs)?;
        Ok(Self {
            obs,
            loglog,
            residuals,
        })
    }

    /// Record one execution and refit (online updating: predicted-vs-actual
    /// residuals feed the refit).
    ///
    /// # Errors
    /// [`CostRefusal::BadInput`] on nonpositive/nonfinite entries.
    pub fn observe(&mut self, o: CostObservation) -> Result<(), CostRefusal> {
        Self::validate_observation(o)?;
        if self.obs.len() >= MAX_COST_OBSERVATIONS {
            return Err(CostRefusal::ObservationLimit {
                limit: MAX_COST_OBSERVATIONS,
            });
        }
        let mut candidate = self.obs.clone();
        let insertion = candidate
            .binary_search_by(|existing| Self::observation_order(existing, &o))
            .unwrap_or_else(|index| index);
        candidate.insert(insertion, o);
        let DerivedFit { loglog, residuals } = Self::derive_fit(&candidate)?;
        self.obs = candidate;
        self.loglog = loglog;
        self.residuals = residuals;
        Ok(())
    }

    /// Observation count.
    #[must_use]
    pub fn n_obs(&self) -> usize {
        self.obs.len()
    }

    fn validate_observation(observation: CostObservation) -> Result<(), CostRefusal> {
        let valid = |value: f64| value.is_finite() && value > 0.0;
        if valid(observation.size) && valid(observation.cost_s) {
            Ok(())
        } else {
            Err(CostRefusal::BadInput)
        }
    }

    fn observation_order(left: &CostObservation, right: &CostObservation) -> core::cmp::Ordering {
        left.size
            .total_cmp(&right.size)
            .then(left.cost_s.total_cmp(&right.cost_s))
    }

    fn finite_sum(
        values: impl IntoIterator<Item = f64>,
        stage: &'static str,
    ) -> Result<f64, CostRefusal> {
        // Kahan summation keeps the centered regression stable while the
        // fixed iteration order preserves snapshot determinism.
        let mut sum = 0.0_f64;
        let mut correction = 0.0_f64;
        for value in values {
            let adjusted = value - correction;
            let next = sum + adjusted;
            correction = (next - sum) - adjusted;
            sum = next;
            if !sum.is_finite() || !correction.is_finite() {
                return Err(CostRefusal::ArithmeticFailure { stage });
            }
        }
        Ok(sum)
    }

    fn derive_fit(observations: &[CostObservation]) -> Result<DerivedFit, CostRefusal> {
        if observations.len() < MIN_OBS {
            return Ok(DerivedFit {
                loglog: None,
                residuals: Vec::new(),
            });
        }
        let xs = observations
            .iter()
            .map(|observation| observation.size.ln())
            .collect::<Vec<_>>();
        let ys = observations
            .iter()
            .map(|observation| observation.cost_s.ln())
            .collect::<Vec<_>>();
        if xs.iter().chain(&ys).any(|value| !value.is_finite()) {
            return Err(CostRefusal::ArithmeticFailure {
                stage: "log transform",
            });
        }
        let n = observations.len() as f64;
        let mean_x = Self::finite_sum(xs.iter().copied(), "mean size")? / n;
        let mean_y = Self::finite_sum(ys.iter().copied(), "mean cost")? / n;
        if !mean_x.is_finite() || !mean_y.is_finite() {
            return Err(CostRefusal::ArithmeticFailure { stage: "means" });
        }
        let centered_xx = Self::finite_sum(
            xs.iter().map(|x| {
                let centered = *x - mean_x;
                centered * centered
            }),
            "centered size variance",
        )?;
        let centered_xy = Self::finite_sum(
            xs.iter()
                .zip(&ys)
                .map(|(x, y)| (*x - mean_x) * (*y - mean_y)),
            "centered covariance",
        )?;
        let x_scale = xs.iter().map(|value| value.abs()).fold(1.0_f64, f64::max);
        let degeneracy_threshold = f64::EPSILON * n * x_scale * x_scale;
        if !centered_xx.is_finite() || !centered_xy.is_finite() || !degeneracy_threshold.is_finite()
        {
            return Err(CostRefusal::ArithmeticFailure {
                stage: "regression moments",
            });
        }
        let slope = if centered_xx <= degeneracy_threshold {
            0.0
        } else {
            centered_xy / centered_xx
        };
        let intercept = mean_y - slope * mean_x;
        if !slope.is_finite() || !intercept.is_finite() {
            return Err(CostRefusal::ArithmeticFailure {
                stage: "regression coefficients",
            });
        }
        let mut residuals = observations
            .iter()
            .zip(xs.iter().zip(&ys))
            .map(|(_, (x, y))| *y - (intercept + slope * *x))
            .collect::<Vec<_>>();
        if residuals.iter().any(|residual| !residual.is_finite()) {
            return Err(CostRefusal::ArithmeticFailure {
                stage: "fit residuals",
            });
        }
        residuals.sort_by(f64::total_cmp);
        Ok(DerivedFit {
            loglog: Some((intercept, slope)),
            residuals,
        })
    }

    fn residual_quantile(&self, p: f64) -> Result<f64, CostRefusal> {
        let last = self
            .residuals
            .len()
            .checked_sub(1)
            .ok_or(CostRefusal::ArithmeticFailure {
                stage: "residual quantile",
            })?;
        let idx = ((self.residuals.len() as f64 - 1.0) * p).round() as usize;
        self.residuals
            .get(idx.min(last))
            .copied()
            .ok_or(CostRefusal::ArithmeticFailure {
                stage: "residual quantile",
            })
    }

    /// Predict the cost distribution at `size`.
    ///
    /// # Errors
    /// [`CostRefusal::InsufficientData`] below [`MIN_OBS`];
    /// [`CostRefusal::BadInput`] for a nonpositive size.
    pub fn predict(&self, size: f64) -> Result<CostPrediction, CostRefusal> {
        let ok = |v: f64| v.is_finite() && v > 0.0;
        if !ok(size) {
            return Err(CostRefusal::BadInput);
        }
        let Some((a, b)) = self.loglog else {
            return Err(CostRefusal::InsufficientData {
                have: self.obs.len(),
                need: MIN_OBS,
            });
        };
        let mu = a + b * size.ln();
        if !mu.is_finite() {
            return Err(CostRefusal::ArithmeticFailure {
                stage: "prediction center",
            });
        }
        let q = |p: f64| -> Result<f64, CostRefusal> {
            let prediction = (mu + self.residual_quantile(p)?).exp();
            if prediction.is_finite() && prediction > 0.0 {
                Ok(prediction)
            } else {
                Err(CostRefusal::ArithmeticFailure {
                    stage: "prediction quantile",
                })
            }
        };
        let (p10, p50, p90) = (q(0.10)?, q(0.50)?, q(0.90)?);
        if p10 > p50 || p50 > p90 {
            return Err(CostRefusal::ArithmeticFailure {
                stage: "prediction quantile order",
            });
        }
        let observed_min = self
            .obs
            .first()
            .ok_or(CostRefusal::ArithmeticFailure {
                stage: "observed size range",
            })?
            .size;
        let observed_max = self
            .obs
            .last()
            .ok_or(CostRefusal::ArithmeticFailure {
                stage: "observed size range",
            })?
            .size;
        let extrapolated = size < observed_min || size > observed_max;
        Ok(CostPrediction {
            p10,
            p50,
            p90,
            n_obs: self.obs.len(),
            extrapolated,
        })
    }

    /// Calibration audit: the fraction of held-out observations whose cost
    /// falls inside this model's [p10, p90] band (the acceptance criterion
    /// "cost predictions within stated quantile bands empirically").
    ///
    /// # Errors
    /// Propagates prediction refusals.
    pub fn calibration(&self, held_out: &[CostObservation]) -> Result<f64, CostRefusal> {
        if held_out.is_empty() {
            return Err(CostRefusal::EmptyEvaluation {
                metric: "calibration",
            });
        }
        if held_out.len() > MAX_COST_EVALUATION_OBSERVATIONS {
            return Err(CostRefusal::EvaluationLimit {
                provided: held_out.len(),
                limit: MAX_COST_EVALUATION_OBSERVATIONS,
            });
        }
        let mut inside = 0usize;
        for o in held_out {
            Self::validate_observation(*o)?;
            let p = self.predict(o.size)?;
            if o.cost_s >= p.p10 && o.cost_s <= p.p90 {
                inside += 1;
            }
        }
        Ok(inside as f64 / held_out.len() as f64)
    }

    /// Mean absolute relative error of the median prediction on a probe set
    /// (the online-updating improvement metric).
    ///
    /// # Errors
    /// Propagates prediction refusals.
    pub fn median_rel_error(&self, probes: &[CostObservation]) -> Result<f64, CostRefusal> {
        if probes.is_empty() {
            return Err(CostRefusal::EmptyEvaluation {
                metric: "median relative error",
            });
        }
        if probes.len() > MAX_COST_EVALUATION_OBSERVATIONS {
            return Err(CostRefusal::EvaluationLimit {
                provided: probes.len(),
                limit: MAX_COST_EVALUATION_OBSERVATIONS,
            });
        }
        let mut relative_errors = Vec::with_capacity(probes.len());
        for o in probes {
            Self::validate_observation(*o)?;
            let p = self.predict(o.size)?;
            let relative_error = (p.p50 - o.cost_s).abs() / o.cost_s;
            if !relative_error.is_finite() {
                return Err(CostRefusal::ArithmeticFailure {
                    stage: "relative error",
                });
            }
            relative_errors.push(relative_error);
        }
        let result =
            Self::finite_sum(relative_errors, "relative error mean")? / probes.len() as f64;
        if result.is_finite() {
            Ok(result)
        } else {
            Err(CostRefusal::ArithmeticFailure {
                stage: "relative error mean",
            })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lcg(seed: &mut u64) -> f64 {
        *seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((*seed >> 11) as f64) / (1u64 << 53) as f64
    }

    /// Synthetic truth: cost = 3e-9 · size^1.5 · noise.
    fn synth(seed: &mut u64, size: f64) -> CostObservation {
        let noise = 0.8 + 0.4 * lcg(seed);
        CostObservation {
            size,
            cost_s: 3e-9 * size.powf(1.5) * noise,
        }
    }

    #[test]
    fn refuses_until_enough_data_then_predicts_with_bands() {
        let mut m = CostModel::new();
        assert!(matches!(
            m.predict(1e6),
            Err(CostRefusal::InsufficientData {
                have: 0,
                need: MIN_OBS
            })
        ));
        let mut seed = 0x5EED_C057_0000_0001u64;
        for i in 0..12 {
            m.observe(synth(&mut seed, 1e4 * f64::from(1 << (i % 6))))
                .unwrap();
        }
        let p = m.predict(3e4).unwrap();
        assert!(p.p10 <= p.p50 && p.p50 <= p.p90);
        assert!(!p.extrapolated);
        assert!(
            m.predict(1e9).unwrap().extrapolated,
            "outside observed sizes"
        );
        assert!(matches!(m.predict(-1.0), Err(CostRefusal::BadInput)));
    }

    #[test]
    fn fits_are_deterministic_regardless_of_arrival_order() {
        let mut seed = 0x5EED_C057_0000_0002u64;
        let obs: Vec<CostObservation> = (0..20)
            .map(|i| synth(&mut seed, 1e3 * f64::from(i + 1)))
            .collect();
        let m1 = CostModel::fit(&obs).unwrap();
        let mut rev = obs.clone();
        rev.reverse();
        let m2 = CostModel::fit(&rev).unwrap();
        let (p1, p2) = (m1.predict(5e3).unwrap(), m2.predict(5e3).unwrap());
        assert_eq!(p1, p2, "arrival order must not change the model");
    }

    #[test]
    fn invalid_online_observation_is_transactional() {
        let observations = [
            CostObservation {
                size: 1.0,
                cost_s: 2.0,
            },
            CostObservation {
                size: 2.0,
                cost_s: 4.0,
            },
            CostObservation {
                size: 4.0,
                cost_s: 8.0,
            },
        ];
        let mut model = CostModel::fit(&observations).unwrap();
        let before = model.predict(3.0).unwrap();
        assert_eq!(
            model.observe(CostObservation {
                size: 3.0,
                cost_s: f64::NAN,
            }),
            Err(CostRefusal::BadInput)
        );
        assert_eq!(model.n_obs(), observations.len());
        assert_eq!(model.predict(3.0).unwrap(), before);
    }

    #[test]
    fn centered_fit_handles_identical_sizes_without_unstable_slope() {
        let model = CostModel::fit(&[
            CostObservation {
                size: 1.0e200,
                cost_s: 1.0e-100,
            },
            CostObservation {
                size: 1.0e200,
                cost_s: 1.0,
            },
            CostObservation {
                size: 1.0e200,
                cost_s: 1.0e100,
            },
        ])
        .unwrap();
        let prediction = model.predict(1.0e200).unwrap();
        assert!(prediction.p10.is_finite());
        assert!(prediction.p50.is_finite());
        assert!(prediction.p90.is_finite());
        assert!(prediction.p10 <= prediction.p50 && prediction.p50 <= prediction.p90);
    }

    #[test]
    fn observations_accept_the_cap_and_refuse_limit_plus_one_transactionally() {
        let observations = (0..MAX_COST_OBSERVATIONS)
            .map(|index| CostObservation {
                size: (index + 1) as f64,
                cost_s: (index + 2) as f64,
            })
            .collect::<Vec<_>>();
        let mut model = CostModel::fit(&observations).unwrap();
        let before = model.predict(17.0).unwrap();
        assert_eq!(model.n_obs(), MAX_COST_OBSERVATIONS);
        assert_eq!(
            model.observe(CostObservation {
                size: 9_999.0,
                cost_s: 10_000.0,
            }),
            Err(CostRefusal::ObservationLimit {
                limit: MAX_COST_OBSERVATIONS,
            })
        );
        assert_eq!(model.n_obs(), MAX_COST_OBSERVATIONS);
        assert_eq!(model.predict(17.0).unwrap(), before);

        let mut over_limit = observations;
        over_limit.push(CostObservation {
            size: 10_001.0,
            cost_s: 10_002.0,
        });
        assert!(matches!(
            CostModel::fit(&over_limit),
            Err(CostRefusal::ObservationLimit {
                limit: MAX_COST_OBSERVATIONS
            })
        ));
    }

    #[test]
    fn empty_and_oversized_evaluations_refuse_without_vacuous_scores() {
        let model = CostModel::fit(&[
            CostObservation {
                size: 1.0,
                cost_s: 1.0,
            },
            CostObservation {
                size: 2.0,
                cost_s: 2.0,
            },
            CostObservation {
                size: 3.0,
                cost_s: 3.0,
            },
        ])
        .unwrap();
        assert_eq!(
            model.calibration(&[]),
            Err(CostRefusal::EmptyEvaluation {
                metric: "calibration"
            })
        );
        assert_eq!(
            model.median_rel_error(&[]),
            Err(CostRefusal::EmptyEvaluation {
                metric: "median relative error"
            })
        );
        let too_many = vec![
            CostObservation {
                size: 1.0,
                cost_s: 1.0,
            };
            MAX_COST_EVALUATION_OBSERVATIONS + 1
        ];
        assert!(matches!(
            model.calibration(&too_many),
            Err(CostRefusal::EvaluationLimit {
                provided,
                limit: MAX_COST_EVALUATION_OBSERVATIONS,
            }) if provided == MAX_COST_EVALUATION_OBSERVATIONS + 1
        ));
    }
}
