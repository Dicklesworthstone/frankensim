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
                write!(f, "sizes and costs must be positive finite (log-log model domain)")
            }
        }
    }
}

/// Minimum observations before predictions are offered.
pub const MIN_OBS: usize = 3;

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
        let mut model = Self::new();
        for &o in observations {
            model.observe(o)?;
        }
        Ok(model)
    }

    /// Record one execution and refit (online updating: predicted-vs-actual
    /// residuals feed the refit).
    ///
    /// # Errors
    /// [`CostRefusal::BadInput`] on nonpositive/nonfinite entries.
    pub fn observe(&mut self, o: CostObservation) -> Result<(), CostRefusal> {
        if !(o.size > 0.0) || !(o.cost_s > 0.0) || !o.size.is_finite() || !o.cost_s.is_finite() {
            return Err(CostRefusal::BadInput);
        }
        self.obs.push(o);
        // Deterministic order regardless of arrival: sort by (size, cost).
        self.obs
            .sort_by(|x, y| x.size.total_cmp(&y.size).then(x.cost_s.total_cmp(&y.cost_s)));
        self.refit();
        Ok(())
    }

    /// Observation count.
    #[must_use]
    pub fn n_obs(&self) -> usize {
        self.obs.len()
    }

    fn refit(&mut self) {
        if self.obs.len() < MIN_OBS {
            self.loglog = None;
            self.residuals.clear();
            return;
        }
        let n = self.obs.len() as f64;
        let (mut sx, mut sy, mut sxx, mut sxy) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        for o in &self.obs {
            let (x, y) = (o.size.ln(), o.cost_s.ln());
            sx += x;
            sy += y;
            sxx += x * x;
            sxy += x * y;
        }
        let denom = n * sxx - sx * sx;
        let (a, b) = if denom.abs() < 1e-30 {
            // All sizes identical: flat model at the mean log-cost.
            (sy / n, 0.0)
        } else {
            let b = (n * sxy - sx * sy) / denom;
            ((sy - b * sx) / n, b)
        };
        self.loglog = Some((a, b));
        self.residuals = self
            .obs
            .iter()
            .map(|o| o.cost_s.ln() - (a + b * o.size.ln()))
            .collect();
        self.residuals.sort_by(f64::total_cmp);
    }

    fn residual_quantile(&self, p: f64) -> f64 {
        let idx = ((self.residuals.len() as f64 - 1.0) * p).round() as usize;
        self.residuals[idx.min(self.residuals.len() - 1)]
    }

    /// Predict the cost distribution at `size`.
    ///
    /// # Errors
    /// [`CostRefusal::InsufficientData`] below [`MIN_OBS`];
    /// [`CostRefusal::BadInput`] for a nonpositive size.
    pub fn predict(&self, size: f64) -> Result<CostPrediction, CostRefusal> {
        if !(size > 0.0) || !size.is_finite() {
            return Err(CostRefusal::BadInput);
        }
        let Some((a, b)) = self.loglog else {
            return Err(CostRefusal::InsufficientData { have: self.obs.len(), need: MIN_OBS });
        };
        let mu = a + b * size.ln();
        let q = |p: f64| (mu + self.residual_quantile(p)).exp();
        let extrapolated = size < self.obs[0].size || size > self.obs[self.obs.len() - 1].size;
        Ok(CostPrediction {
            p10: q(0.10),
            p50: q(0.50),
            p90: q(0.90),
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
            return Ok(1.0);
        }
        let mut inside = 0usize;
        for o in held_out {
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
            return Ok(0.0);
        }
        let mut total = 0.0f64;
        for o in probes {
            let p = self.predict(o.size)?;
            total += (p.p50 - o.cost_s).abs() / o.cost_s;
        }
        Ok(total / probes.len() as f64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lcg(seed: &mut u64) -> f64 {
        *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((*seed >> 11) as f64) / (1u64 << 53) as f64
    }

    /// Synthetic truth: cost = 3e-9 · size^1.5 · noise.
    fn synth(seed: &mut u64, size: f64) -> CostObservation {
        let noise = 0.8 + 0.4 * lcg(seed);
        CostObservation { size, cost_s: 3e-9 * size.powf(1.5) * noise }
    }

    #[test]
    fn refuses_until_enough_data_then_predicts_with_bands() {
        let mut m = CostModel::new();
        assert!(matches!(
            m.predict(1e6),
            Err(CostRefusal::InsufficientData { have: 0, need: MIN_OBS })
        ));
        let mut seed = 0x5EED_C057_0000_0001u64;
        for i in 0..12 {
            m.observe(synth(&mut seed, 1e4 * f64::from(1 << (i % 6)))).unwrap();
        }
        let p = m.predict(3e4).unwrap();
        assert!(p.p10 <= p.p50 && p.p50 <= p.p90);
        assert!(!p.extrapolated);
        assert!(m.predict(1e9).unwrap().extrapolated, "outside observed sizes");
        assert!(matches!(m.predict(-1.0), Err(CostRefusal::BadInput)));
    }

    #[test]
    fn fits_are_deterministic_regardless_of_arrival_order() {
        let mut seed = 0x5EED_C057_0000_0002u64;
        let obs: Vec<CostObservation> =
            (0..20).map(|i| synth(&mut seed, 1e3 * f64::from(i + 1))).collect();
        let m1 = CostModel::fit(&obs).unwrap();
        let mut rev = obs.clone();
        rev.reverse();
        let m2 = CostModel::fit(&rev).unwrap();
        let (p1, p2) = (m1.predict(5e3).unwrap(), m2.predict(5e3).unwrap());
        assert_eq!(p1, p2, "arrival order must not change the model");
    }
}
