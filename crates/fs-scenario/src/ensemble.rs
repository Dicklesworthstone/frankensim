//! Stochastic scenario ensembles: seeded, replayable-from-seed generators
//! for environmental variability — wind gusts (Dryden spectra), ground
//! motions (Kanai–Tajimi spectral model), fluid-property bands (Carreau
//! parameter families). Realizations are BITWISE reproducible functions of
//! the complete canonical ensemble specification, member identity, and the
//! versioned Philox stream/synthesis semantics. The explicit recipe receipt
//! and seed-tree policy remain tracked by `frankensim-sj31i.39`.

use crate::ScenarioError;
use crate::scenario::Violation;
use fs_math::det;
use fs_qty::{Dims, QtyAny};
use fs_rand::StreamKey;
use std::fmt;

const TIME_DIMS: Dims = Dims([0, 0, 1, 0, 0, 0]);
/// fs-rand kernel ids for ensemble draws (stable across runs — part of
/// the logical identity, never a thread id).
const KERNEL_GUST: u32 = 0x5C01;
const KERNEL_GROUND: u32 = 0x5C02;
const KERNEL_CARREAU: u32 = 0x5C03;

/// The spectral / band model behind an ensemble.
#[derive(Debug, Clone, PartialEq)]
pub enum SpectrumModel {
    /// Dryden longitudinal gust spectrum:
    /// `S(ω) = σ²·(2L/(πV)) / (1 + (Lω/V)²)`.
    Dryden {
        /// Turbulence intensity σ (m/s).
        sigma: QtyAny,
        /// Length scale L (m).
        length_scale: QtyAny,
        /// Mean wind speed V (m/s).
        mean_speed: QtyAny,
    },
    /// Kanai–Tajimi ground-acceleration spectrum:
    /// `S(ω) = S₀·(1 + 4ζ²r²)/((1−r²)² + 4ζ²r²)`, `r = ω/ω_g`.
    KanaiTajimi {
        /// Bedrock intensity S₀ ((m/s²)²·s/rad, carried as a raw factor).
        s0: f64,
        /// Ground natural frequency ω_g (rad/s).
        omega_g: QtyAny,
        /// Ground damping ζ_g.
        zeta_g: f64,
    },
    /// A Carreau-fluid parameter band (the vessel robustness sweep):
    /// members sample each parameter uniformly inside its band.
    CarreauBand {
        /// Zero-shear viscosity band (Pa·s).
        eta_zero: [QtyAny; 2],
        /// Infinite-shear viscosity band (Pa·s).
        eta_inf: [QtyAny; 2],
        /// Relaxation-time band (s).
        lambda: [QtyAny; 2],
        /// Power-index band (dimensionless).
        n: [f64; 2],
    },
}

impl SpectrumModel {
    /// One-sided target PSD at angular frequency ω (rad/s). Parameter-band
    /// models are not spectra and are refused.
    ///
    /// # Errors
    /// Returns [`ScenarioError`] when the frequency or any model parameter is
    /// non-finite, dimensionally invalid, outside its physical domain, or the
    /// finite inputs produce an unrepresentable PSD value.
    pub fn try_psd(&self, omega: f64) -> Result<f64, ScenarioError> {
        if !omega.is_finite() || omega < 0.0 {
            return Err(ScenarioError::Evaluate {
                what: format!(
                    "one-sided PSD angular frequency {omega} must be finite and non-negative"
                ),
            });
        }
        let value = match self {
            SpectrumModel::Dryden {
                sigma,
                length_scale,
                mean_speed,
            } => {
                let velocity_dims = Dims([1, 0, -1, 0, 0, 0]);
                let length_dims = Dims([1, 0, 0, 0, 0, 0]);
                if sigma.dims != velocity_dims
                    || mean_speed.dims != velocity_dims
                    || length_scale.dims != length_dims
                {
                    return Err(ScenarioError::Evaluate {
                        what: "Dryden PSD parameters have invalid physical dimensions".to_string(),
                    });
                }
                let (s, l, v) = (sigma.value, length_scale.value, mean_speed.value);
                if !(s.is_finite()
                    && l.is_finite()
                    && v.is_finite()
                    && s > 0.0
                    && l > 0.0
                    && v > 0.0)
                {
                    return Err(ScenarioError::Evaluate {
                        what: "Dryden PSD parameters must be finite and positive".to_string(),
                    });
                }
                let x = l * omega / v;
                s * s * (2.0 * l / (core::f64::consts::PI * v)) / (1.0 + x * x)
            }
            SpectrumModel::KanaiTajimi {
                s0,
                omega_g,
                zeta_g,
            } => {
                if omega_g.dims != Dims([0, 0, -1, 0, 0, 0]) {
                    return Err(ScenarioError::Evaluate {
                        what: "Kanai–Tajimi ground frequency has invalid physical dimensions"
                            .to_string(),
                    });
                }
                if !(s0.is_finite()
                    && omega_g.value.is_finite()
                    && zeta_g.is_finite()
                    && *s0 > 0.0
                    && omega_g.value > 0.0
                    && *zeta_g > 0.0)
                {
                    return Err(ScenarioError::Evaluate {
                        what: "Kanai–Tajimi PSD parameters must be finite and positive".to_string(),
                    });
                }
                let r = omega / omega_g.value;
                let r2 = r * r;
                let four_z2_r2 = 4.0 * zeta_g * zeta_g * r2;
                s0 * (1.0 + four_z2_r2) / ((1.0 - r2) * (1.0 - r2) + four_z2_r2)
            }
            SpectrumModel::CarreauBand { .. } => {
                return Err(ScenarioError::Evaluate {
                    what: "Carreau parameter bands do not define a spectral density".to_string(),
                });
            }
        };
        if !value.is_finite() || value < 0.0 {
            return Err(ScenarioError::Evaluate {
                what: format!(
                    "finite model inputs produced an invalid PSD at angular frequency {omega}"
                ),
            });
        }
        Ok(value)
    }

    fn kernel(&self) -> u32 {
        match self {
            SpectrumModel::Dryden { .. } => KERNEL_GUST,
            SpectrumModel::KanaiTajimi { .. } => KERNEL_GROUND,
            SpectrumModel::CarreauBand { .. } => KERNEL_CARREAU,
        }
    }
}

/// A seeded ensemble specification.
#[derive(Debug, Clone, PartialEq)]
pub struct StochasticEnsemble {
    /// Ensemble name (IR identity).
    pub name: String,
    /// The study seed feeding the Philox streams.
    pub seed: u64,
    /// Member count.
    pub members: u32,
    /// Realization duration (s); ignored by band models.
    pub duration: QtyAny,
    /// Sample step (s); ignored by band models.
    pub dt: QtyAny,
    /// The model.
    pub model: SpectrumModel,
}

#[derive(Clone, Copy)]
struct EnsembleDiagnosticContext<'a>(&'a str);

impl fmt::Display for EnsembleDiagnosticContext<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "ensemble {:?}", self.0)
    }
}

/// One realized member: a sampled time series (spectral models) or a
/// parameter draw (band models).
#[derive(Debug, Clone, PartialEq)]
pub struct Realization {
    /// Sample times (s); empty for band models.
    pub times: Vec<f64>,
    /// Sampled values (spectral models) or parameters (band models).
    pub values: Vec<f64>,
}

/// Explicit admission budget for one realization.
///
/// `max_work` counts generated sample timestamps, coefficient pairs, and
/// sample-by-harmonic accumulation steps. It is deliberately independent of
/// wall-clock speed so the same request is admitted on every machine.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RealizationBudget {
    /// Maximum number of output samples.
    /// Band-model realizations have zero samples.
    pub max_samples: usize,
    /// Maximum deterministic work units.
    pub max_work: usize,
}

/// Conservative default used by [`StochasticEnsemble::realize`].
pub const DEFAULT_REALIZATION_BUDGET: RealizationBudget = RealizationBudget {
    max_samples: 65_536,
    max_work: 16_777_216,
};

impl Default for RealizationBudget {
    fn default() -> Self {
        DEFAULT_REALIZATION_BUDGET
    }
}

fn reserve_exact<T>(values: &mut Vec<T>, count: usize, context: &str) -> Result<(), ScenarioError> {
    values
        .try_reserve_exact(count)
        .map_err(|allocation_error| ScenarioError::Evaluate {
            what: format!(
                "{context}: allocation for {count} elements was refused: {allocation_error}"
            ),
        })
}

impl StochasticEnsemble {
    /// Realize member `member` — a bitwise-reproducible function of the full
    /// canonical ensemble specification, member identity, and implementation
    /// stream/synthesis semantics. The spectral representation is
    /// `x(t) = Σₖ √(S(ωₖ)Δω)·(aₖ cos ωₖt + bₖ sin ωₖt)`.
    ///
    /// # Errors
    /// [`ScenarioError`] for dimension/shape defects in the spec.
    pub fn realize(&self, member: u32) -> Result<Realization, ScenarioError> {
        self.realize_with_budget(member, RealizationBudget::default())
    }

    /// Realize a member under an explicit deterministic sample/work budget.
    ///
    /// Public fields make it possible to construct malformed ensembles
    /// without calling [`StochasticEnsemble::check`]. This method therefore
    /// validates the complete model independently before drawing or
    /// allocating anything.
    ///
    /// # Errors
    /// [`ScenarioError`] for invalid dimensions/domains, out-of-range members,
    /// arithmetic overflow, budget refusal, or allocation refusal.
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_precision_loss,
        clippy::cast_sign_loss
    )]
    pub fn realize_with_budget(
        &self,
        member: u32,
        budget: RealizationBudget,
    ) -> Result<Realization, ScenarioError> {
        let mut violations = Vec::new();
        self.check(&mut violations);
        if let Some(first) = violations.first() {
            return Err(ScenarioError::Evaluate {
                what: format!("ensemble {:?} is invalid: {}", self.name, first.what),
            });
        }
        if member >= self.members {
            return Err(ScenarioError::Evaluate {
                what: format!(
                    "member {member} out of range (ensemble {:?} has {})",
                    self.name, self.members
                ),
            });
        }
        let key = StreamKey {
            seed: self.seed,
            kernel: self.model.kernel(),
            tile: member,
        };
        let mut stream = key.stream();
        if let SpectrumModel::CarreauBand {
            eta_zero,
            eta_inf,
            lambda,
            n,
        } = &self.model
        {
            if budget.max_work < 4 {
                return Err(ScenarioError::Evaluate {
                    what: format!(
                        "ensemble {:?}: requested work 4 exceeds budget {}",
                        self.name, budget.max_work
                    ),
                });
            }
            let draw = |lo: f64, hi: f64, s: &mut fs_rand::Stream| lo + (hi - lo) * s.next_f64();
            let mut values = Vec::new();
            reserve_exact(&mut values, 4, "Carreau realization")?;
            values.push(draw(eta_zero[0].value, eta_zero[1].value, &mut stream));
            values.push(draw(eta_inf[0].value, eta_inf[1].value, &mut stream));
            values.push(draw(lambda[0].value, lambda[1].value, &mut stream));
            values.push(draw(n[0], n[1], &mut stream));
            if values.iter().any(|value| !value.is_finite()) {
                return Err(ScenarioError::Evaluate {
                    what: format!("ensemble {:?}: Carreau draw became non-finite", self.name),
                });
            }
            return Ok(Realization {
                times: Vec::new(),
                values,
            });
        }
        let ratio = self.duration.value / self.dt.value;
        if !ratio.is_finite() || ratio < 1.0 {
            return Err(ScenarioError::Evaluate {
                what: format!(
                    "ensemble {:?}: duration/dt ratio {ratio} must be finite and >= 1",
                    self.name
                ),
            });
        }
        let rounded_samples = ratio.round();
        if rounded_samples < 2.0 {
            return Err(ScenarioError::Evaluate {
                what: format!(
                    "ensemble {:?}: spectral grid rounds to {rounded_samples:.0} sample; at least 2 are required",
                    self.name
                ),
            });
        }
        if rounded_samples >= usize::MAX as f64 {
            return Err(ScenarioError::Evaluate {
                what: format!(
                    "ensemble {:?}: sample count cannot be represented as usize",
                    self.name
                ),
            });
        }
        if rounded_samples > budget.max_samples as f64 {
            return Err(ScenarioError::Evaluate {
                what: format!(
                    "ensemble {:?}: requested {:.0} samples exceeds budget {}",
                    self.name, rounded_samples, budget.max_samples
                ),
            });
        }
        let n_samples = rounded_samples as usize;
        let n_harmonics = n_samples / 2;
        let synthesis_work =
            n_samples
                .checked_mul(n_harmonics)
                .ok_or_else(|| ScenarioError::Evaluate {
                    what: format!("ensemble {:?}: realization work overflowed", self.name),
                })?;
        let work = synthesis_work
            .checked_add(n_harmonics)
            .and_then(|work| work.checked_add(n_samples))
            .ok_or_else(|| ScenarioError::Evaluate {
                what: format!("ensemble {:?}: realization work overflowed", self.name),
            })?;
        if work > budget.max_work {
            return Err(ScenarioError::Evaluate {
                what: format!(
                    "ensemble {:?}: requested work {work} exceeds budget {}",
                    self.name, budget.max_work
                ),
            });
        }
        let d_omega = 2.0 * core::f64::consts::PI / (n_samples as f64 * self.dt.value);
        if !d_omega.is_finite() || d_omega <= 0.0 {
            return Err(ScenarioError::Evaluate {
                what: format!("ensemble {:?}: frequency spacing is invalid", self.name),
            });
        }
        // Draw the Gaussian coefficient pairs in a fixed order (bitwise
        // determinism comes from the fixed draw and summation order).
        let mut coeffs = Vec::new();
        reserve_exact(&mut coeffs, n_harmonics, "spectral coefficients")?;
        for k in 1..=n_harmonics {
            let omega = k as f64 * d_omega;
            let spectral_density =
                self.model
                    .try_psd(omega)
                    .map_err(|error| ScenarioError::Evaluate {
                        what: format!(
                            "ensemble {:?}: PSD is invalid at angular frequency {omega}: {error}",
                            self.name
                        ),
                    })?;
            let amp = det::sqrt(spectral_density * d_omega);
            if !amp.is_finite() {
                return Err(ScenarioError::Evaluate {
                    what: format!(
                        "ensemble {:?}: spectral amplitude overflowed at angular frequency {omega}",
                        self.name
                    ),
                });
            }
            let cosine_coefficient = amp * stream.next_normal();
            let sine_coefficient = amp * stream.next_normal();
            if !cosine_coefficient.is_finite() || !sine_coefficient.is_finite() {
                return Err(ScenarioError::Evaluate {
                    what: format!(
                        "ensemble {:?}: spectral coefficient overflowed at angular frequency {omega}",
                        self.name
                    ),
                });
            }
            coeffs.push((omega, cosine_coefficient, sine_coefficient));
        }
        let mut times = Vec::new();
        reserve_exact(&mut times, n_samples, "realization times")?;
        for sample in 0..n_samples {
            let time = sample as f64 * self.dt.value;
            if !time.is_finite() {
                return Err(ScenarioError::Evaluate {
                    what: format!("ensemble {:?}: sample time overflowed", self.name),
                });
            }
            times.push(time);
        }
        let mut values = Vec::new();
        reserve_exact(&mut values, n_samples, "realization values")?;
        for &time in &times {
            let mut accumulated = 0.0f64;
            for &(omega, cosine_coefficient, sine_coefficient) in &coeffs {
                accumulated += cosine_coefficient * det::cos(omega * time)
                    + sine_coefficient * det::sin(omega * time);
            }
            if !accumulated.is_finite() {
                return Err(ScenarioError::Evaluate {
                    what: format!("ensemble {:?}: synthesized sample is non-finite", self.name),
                });
            }
            values.push(accumulated);
        }
        Ok(Realization { times, values })
    }

    /// Structural validation.
    pub fn check(&self, out: &mut Vec<Violation>) {
        let ctx = EnsembleDiagnosticContext(self.name.as_str());
        if self.name.is_empty() {
            out.push(Violation {
                code: "ensemble-name-empty",
                what: "ensemble identity is empty".to_string(),
                fix: "give every ensemble a nonempty exact UTF-8 name".to_string(),
            });
        }
        if self.members == 0 {
            out.push(Violation {
                code: "ensemble-empty",
                what: format!("{ctx}: zero members"),
                fix: "request at least one member".to_string(),
            });
        }
        // Band models do not sample this grid, but duration/dt still travel in
        // canonical IR and therefore remain part of the public artifact. They
        // must be finite, dimensionally coherent placeholders so an admitted
        // ensemble can always be serialized and reparsed losslessly.
        if self.dt.dims != TIME_DIMS || self.duration.dims != TIME_DIMS {
            out.push(Violation {
                code: "ensemble-time-dims",
                what: format!("{ctx}: duration/dt must be times (seconds)"),
                fix: "give duration and dt the SI dimensions of time".to_string(),
            });
        }
        let dt_ok = self.dt.value.is_finite() && self.dt.value > 0.0;
        let duration_ok = self.duration.value.is_finite()
            && self.duration.value > 0.0
            && self.duration.value >= self.dt.value;
        if !dt_ok || !duration_ok {
            out.push(Violation {
                code: "ensemble-time-range",
                what: format!(
                    "{ctx}: dt {} vs duration {}",
                    self.dt.value, self.duration.value
                ),
                fix: "choose finite values satisfying 0 < dt <= duration".to_string(),
            });
        }
        if dt_ok && duration_ok && !matches!(self.model, SpectrumModel::CarreauBand { .. }) {
            let ratio = self.duration.value / self.dt.value;
            let rounded = ratio.round();
            if !ratio.is_finite() || rounded < 2.0 || rounded >= usize::MAX as f64 {
                out.push(Violation {
                    code: "ensemble-spectral-grid",
                    what: format!(
                        "{ctx}: spectral duration/dt ratio {ratio} does not define a representable grid with at least two samples"
                    ),
                    fix: "choose finite duration/dt whose rounded sample count is in [2, usize::MAX)"
                        .to_string(),
                });
            }
        }
        match &self.model {
            SpectrumModel::Dryden {
                sigma,
                length_scale,
                mean_speed,
            } => {
                expect_dims(&ctx, "sigma", sigma, Dims([1, 0, -1, 0, 0, 0]), out);
                expect_dims(
                    &ctx,
                    "length scale",
                    length_scale,
                    Dims([1, 0, 0, 0, 0, 0]),
                    out,
                );
                expect_dims(
                    &ctx,
                    "mean speed",
                    mean_speed,
                    Dims([1, 0, -1, 0, 0, 0]),
                    out,
                );
                // `psd` divides by `mean_speed` (`x = l·ω/v`, and `2l/(πv)`), so
                // a zero/negative/non-finite speed makes every realization
                // inf/NaN — validate must reject it, not admit a NaN ensemble.
                // sigma (intensity) and length scale must likewise be positive.
                let positive = |v: f64| v.is_finite() && v > 0.0;
                if !positive(sigma.value)
                    || !positive(length_scale.value)
                    || !positive(mean_speed.value)
                {
                    out.push(Violation {
                        code: "ensemble-dryden-params",
                        what: format!(
                            "{ctx}: sigma, length scale, and mean speed must be positive"
                        ),
                        fix: "supply positive Dryden intensity, length scale, and mean speed"
                            .to_string(),
                    });
                }
            }
            SpectrumModel::KanaiTajimi {
                s0,
                omega_g,
                zeta_g,
            } => {
                expect_dims(&ctx, "omega_g", omega_g, Dims([0, 0, -1, 0, 0, 0]), out);
                // `psd` divides by `omega_g` (`r = ω/ω_g`), so a zero/negative
                // ground frequency makes every realization NaN — reject it
                // alongside S0 and zeta_g rather than admit a NaN ensemble.
                let positive = |v: f64| v.is_finite() && v > 0.0;
                if !positive(*s0) || !positive(*zeta_g) || !positive(omega_g.value) {
                    out.push(Violation {
                        code: "ensemble-kt-params",
                        what: format!("{ctx}: S0, zeta_g, and omega_g must be positive"),
                        fix:
                            "supply positive Kanai–Tajimi intensity, damping, and ground frequency"
                                .to_string(),
                    });
                }
            }
            SpectrumModel::CarreauBand {
                eta_zero,
                eta_inf,
                lambda,
                n,
            } => {
                let visc = Dims([-1, 1, -1, 0, 0, 0]);
                expect_dims(&ctx, "eta_zero lo", &eta_zero[0], visc, out);
                expect_dims(&ctx, "eta_zero hi", &eta_zero[1], visc, out);
                expect_dims(&ctx, "eta_inf lo", &eta_inf[0], visc, out);
                expect_dims(&ctx, "eta_inf hi", &eta_inf[1], visc, out);
                expect_dims(&ctx, "lambda lo", &lambda[0], TIME_DIMS, out);
                expect_dims(&ctx, "lambda hi", &lambda[1], TIME_DIMS, out);
                for (lo, hi, name) in [
                    (eta_zero[0].value, eta_zero[1].value, "eta_zero"),
                    (eta_inf[0].value, eta_inf[1].value, "eta_inf"),
                    (lambda[0].value, lambda[1].value, "lambda"),
                ] {
                    if !(lo.is_finite() && hi.is_finite() && lo > 0.0 && hi > 0.0) {
                        out.push(Violation {
                            code: "ensemble-carreau-positive-finite",
                            what: format!(
                                "{ctx}: {name} band [{lo}, {hi}] must contain finite positive values"
                            ),
                            fix: "use finite, strictly positive Carreau viscosity/time bounds"
                                .to_string(),
                        });
                    } else if lo > hi {
                        out.push(Violation {
                            code: "ensemble-band-order",
                            what: format!("{ctx}: {name} band [{lo}, {hi}] is inverted"),
                            fix: "order every band as [low, high]".to_string(),
                        });
                    }
                }
                if !(n[0].is_finite()
                    && n[1].is_finite()
                    && n[0] > 0.0
                    && n[1] > 0.0
                    && n[0] <= 1.0
                    && n[1] <= 1.0)
                {
                    out.push(Violation {
                        code: "ensemble-carreau-power-index",
                        what: format!(
                            "{ctx}: shear-thinning power-index band [{}, {}] must lie in (0, 1]",
                            n[0], n[1]
                        ),
                        fix:
                            "use finite Carreau power-index bounds satisfying 0 < low <= high <= 1"
                                .to_string(),
                    });
                } else if n[0] > n[1] {
                    out.push(Violation {
                        code: "ensemble-band-order",
                        what: format!("{ctx}: n band [{}, {}] is inverted", n[0], n[1]),
                        fix: "order every band as [low, high]".to_string(),
                    });
                }
                if eta_zero[0].value.is_finite()
                    && eta_inf[1].value.is_finite()
                    && eta_zero[0].value < eta_inf[1].value
                {
                    out.push(Violation {
                        code: "ensemble-carreau-viscosity-order",
                        what: format!(
                            "{ctx}: eta_zero low {} is below eta_inf high {}; independent draws could violate eta_zero >= eta_inf",
                            eta_zero[0].value, eta_inf[1].value
                        ),
                        fix: "separate the bands so every zero-shear viscosity is at least every infinite-shear viscosity"
                            .to_string(),
                    });
                }
            }
        }
    }
}

fn expect_dims<C: fmt::Display + ?Sized>(
    ctx: &C,
    name: &str,
    q: &QtyAny,
    expected: Dims,
    out: &mut Vec<Violation>,
) {
    if q.dims != expected {
        out.push(Violation {
            code: "ensemble-dims",
            what: format!(
                "{ctx}: {name} has dimensions {:?}, expected {:?}",
                q.dims.0, expected.0
            ),
            fix: format!("express {name} in coherent SI units"),
        });
    }
}

#[cfg(test)]
mod validation_internal_tests {
    use super::{EnsembleDiagnosticContext, SpectrumModel, StochasticEnsemble, TIME_DIMS};
    use fs_qty::{Dims, QtyAny};

    #[test]
    fn diagnostic_context_is_borrowed_and_output_stable() {
        let name = String::from("gust");
        let context = EnsembleDiagnosticContext(name.as_str());
        assert_eq!(format!("{context}"), "ensemble \"gust\"");

        let ensemble = StochasticEnsemble {
            name,
            seed: 1,
            members: 1,
            duration: QtyAny::new(2.0, TIME_DIMS),
            dt: QtyAny::new(1.0, TIME_DIMS),
            model: SpectrumModel::Dryden {
                sigma: QtyAny::new(1.0, Dims([1, 0, -1, 0, 0, 0])),
                length_scale: QtyAny::new(1.0, Dims([1, 0, 0, 0, 0, 0])),
                mean_speed: QtyAny::new(1.0, Dims([1, 0, -1, 0, 0, 0])),
            },
        };
        let mut findings = Vec::new();
        ensemble.check(&mut findings);
        assert!(findings.is_empty());
    }
}
