//! Stochastic scenario ensembles: seeded, replayable-from-seed generators
//! for environmental variability — wind gusts (Dryden spectra), ground
//! motions (Kanai–Tajimi spectral model), fluid-property bands (Carreau
//! parameter families). Realizations are BITWISE reproducible: member k
//! of ensemble (seed s) is a pure function of (s, model, k) via Philox
//! streams keyed by logical identity (fs-rand).

use crate::ScenarioError;
use crate::scenario::Violation;
use fs_math::det;
use fs_qty::{Dims, QtyAny};
use fs_rand::StreamKey;

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
    /// One-sided target PSD at angular frequency ω (rad/s); zero for
    /// band models.
    #[must_use]
    pub fn psd(&self, omega: f64) -> f64 {
        match self {
            SpectrumModel::Dryden {
                sigma,
                length_scale,
                mean_speed,
            } => {
                let (s, l, v) = (sigma.value, length_scale.value, mean_speed.value);
                let x = l * omega / v;
                s * s * (2.0 * l / (core::f64::consts::PI * v)) / (1.0 + x * x)
            }
            SpectrumModel::KanaiTajimi {
                s0,
                omega_g,
                zeta_g,
            } => {
                let r = omega / omega_g.value;
                let r2 = r * r;
                let four_z2_r2 = 4.0 * zeta_g * zeta_g * r2;
                s0 * (1.0 + four_z2_r2) / ((1.0 - r2) * (1.0 - r2) + four_z2_r2)
            }
            SpectrumModel::CarreauBand { .. } => 0.0,
        }
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

/// One realized member: a sampled time series (spectral models) or a
/// parameter draw (band models).
#[derive(Debug, Clone, PartialEq)]
pub struct Realization {
    /// Sample times (s); empty for band models.
    pub times: Vec<f64>,
    /// Sampled values (spectral models) or parameters (band models).
    pub values: Vec<f64>,
}

impl StochasticEnsemble {
    /// Realize member `member` — a pure, bitwise-reproducible function of
    /// `(seed, model, member)` via the spectral representation method with
    /// Gaussian coefficients: `x(t) = Σₖ √(S(ωₖ)Δω)·(aₖ cos ωₖt + bₖ sin ωₖt)`.
    ///
    /// # Errors
    /// [`ScenarioError`] for dimension/shape defects in the spec.
    pub fn realize(&self, member: u32) -> Result<Realization, ScenarioError> {
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
            let draw = |lo: f64, hi: f64, s: &mut fs_rand::Stream| lo + (hi - lo) * s.next_f64();
            let values = vec![
                draw(eta_zero[0].value, eta_zero[1].value, &mut stream),
                draw(eta_inf[0].value, eta_inf[1].value, &mut stream),
                draw(lambda[0].value, lambda[1].value, &mut stream),
                draw(n[0], n[1], &mut stream),
            ];
            return Ok(Realization {
                times: Vec::new(),
                values,
            });
        }
        if self.dt.dims != TIME_DIMS || self.duration.dims != TIME_DIMS {
            return Err(ScenarioError::Dimensions {
                context: format!("ensemble {:?} duration/dt", self.name),
                expected: TIME_DIMS.0,
                got: self.dt.dims.0,
            });
        }
        let dt_ok = self.dt.value.is_finite() && self.dt.value > 0.0;
        if !dt_ok || self.duration.value < self.dt.value {
            return Err(ScenarioError::Evaluate {
                what: format!("ensemble {:?}: need 0 < dt <= duration", self.name),
            });
        }
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let n_samples = (self.duration.value / self.dt.value).round() as usize;
        let n_harmonics = n_samples / 2;
        let d_omega = 2.0 * core::f64::consts::PI / (n_samples as f64 * self.dt.value);
        // Draw the Gaussian coefficient pairs in a fixed order (bitwise
        // determinism comes from the fixed draw and summation order).
        let mut coeffs = Vec::with_capacity(n_harmonics);
        for k in 1..=n_harmonics {
            let omega = k as f64 * d_omega;
            let amp = det::sqrt(self.model.psd(omega) * d_omega);
            coeffs.push((
                omega,
                amp * stream.next_normal(),
                amp * stream.next_normal(),
            ));
        }
        let times: Vec<f64> = (0..n_samples).map(|j| j as f64 * self.dt.value).collect();
        let values: Vec<f64> = times
            .iter()
            .map(|&t| {
                let mut acc = 0.0f64;
                for &(omega, a, b) in &coeffs {
                    acc += a * det::cos(omega * t) + b * det::sin(omega * t);
                }
                acc
            })
            .collect();
        Ok(Realization { times, values })
    }

    /// Structural validation.
    pub fn check(&self, out: &mut Vec<Violation>) {
        let ctx = format!("ensemble {:?}", self.name);
        if self.members == 0 {
            out.push(Violation {
                code: "ensemble-empty",
                what: format!("{ctx}: zero members"),
                fix: "request at least one member".to_string(),
            });
        }
        let spectral = !matches!(self.model, SpectrumModel::CarreauBand { .. });
        if spectral {
            if self.dt.dims != TIME_DIMS || self.duration.dims != TIME_DIMS {
                out.push(Violation {
                    code: "ensemble-time-dims",
                    what: format!("{ctx}: duration/dt must be times (seconds)"),
                    fix: "give duration and dt the SI dimensions of time".to_string(),
                });
            }
            let dt_ok = self.dt.value.is_finite() && self.dt.value > 0.0;
            if !dt_ok || self.duration.value < self.dt.value {
                out.push(Violation {
                    code: "ensemble-time-range",
                    what: format!(
                        "{ctx}: dt {} vs duration {}",
                        self.dt.value, self.duration.value
                    ),
                    fix: "choose 0 < dt <= duration".to_string(),
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
                if !positive(sigma.value) || !positive(length_scale.value) || !positive(mean_speed.value)
                {
                    out.push(Violation {
                        code: "ensemble-dryden-params",
                        what: format!("{ctx}: sigma, length scale, and mean speed must be positive"),
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
                        fix: "supply positive Kanai–Tajimi intensity, damping, and ground frequency"
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
                    (n[0], n[1], "n"),
                ] {
                    if lo > hi {
                        out.push(Violation {
                            code: "ensemble-band-order",
                            what: format!("{ctx}: {name} band [{lo}, {hi}] is inverted"),
                            fix: "order every band as [low, high]".to_string(),
                        });
                    }
                }
            }
        }
    }
}

fn expect_dims(ctx: &str, name: &str, q: &QtyAny, expected: Dims, out: &mut Vec<Violation>) {
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
