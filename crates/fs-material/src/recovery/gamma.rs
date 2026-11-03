//! Gamma-rate expectations without integrating a singular rate density.
//!
//! For shape a, rate scale theta and N equal CSTRs, Fubini's theorem applied
//! to two nonnegative integrands gives surviving fraction
//! S = E[(1 + b Y)^(-a)], Y ~ Gamma(N, 1), b = theta * mean_residence / N.
//! This combines the gamma Laplace transform (Vinnett and Waters, 2020,
//! https://doi.org/10.3390/min10100913, Table 1) with the Erlang RTD (Fogler,
//! https://websites.umich.edu/~elements/5e/17chap/Ch17_Web-Additional%20Material.pdf,
//! equations 17-4/17-5). Rate classes retain their rates through the train.
//!
//! In x=ln(Y), both integrands are smooth and log-concave even for a<1.
//! A tangent to log(f) bounds each infinite tail by f(x)/|d ln(f)/dx|.
//! Bounded adaptive Simpson integration estimates the remaining numerical
//! error. These are binary64 estimates, not certified interval enclosures.

use super::RecoveryError;
use fs_math::det;
use fs_qty::Dimensionless;

const MAX_EVALUATIONS: usize = 65_536;
const MAX_DEPTH: usize = 40;
const RELATIVE_QUADRATURE_TARGET: f64 = 2.0e-14;
const RELATIVE_TAIL_TARGET: f64 = 2.0e-16;

/// Numerical diagnostics for the gamma-rate/Erlang residence integral.
/// They describe numerical evaluation only, not parameter or model uncertainty.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GammaIntegrationReport {
    /// Combined quadrature, truncation and rounding estimate on each fraction.
    pub absolute_error_estimate: Dimensionless,
    /// Nominal analytic bound on the omitted infinite tails, evaluated in f64.
    pub omitted_tail_bound: Dimensionless,
    /// Independently integrated recovered plus surviving fraction minus one.
    pub independent_closure_residual: Dimensionless,
    /// Adaptive quadrature evaluations; at most 65,536 per response. Separate
    /// mode searches (96 iterations each) and tail searches (32 per side) are bounded.
    pub evaluations: usize,
}

pub(super) struct GammaFractions {
    pub converted: f64,
    pub surviving: f64,
    pub report: GammaIntegrationReport,
}

// Correct the rounded addition; retains x smaller than one ULP at one.
fn log_one_plus(x: f64) -> f64 {
    let sum = 1.0 + x;
    if sum == 1.0 {
        x
    } else {
        det::ln(sum) * (x / (sum - 1.0))
    }
}

fn softplus(x: f64) -> f64 {
    if x > 36.0 {
        x + log_one_plus(det::exp(-x))
    } else {
        log_one_plus(det::exp(x))
    }
}

fn log_softplus(x: f64) -> f64 {
    // Relative correction is below binary64 precision at this endpoint.
    if x < -36.0 { x } else { det::ln(softplus(x)) }
}

/// Exact gamma Laplace-transform response, evaluated without forming theta*t.
pub(super) fn batch(shape: f64, rate: f64, time: f64) -> Result<f64, RecoveryError> {
    let log_exposure = det::ln(shape) + log_softplus(det::ln(rate) + det::ln(time));
    let response = -det::expm1(-det::exp(log_exposure));
    if response.is_finite() && response > 0.0 && response <= 1.0 {
        Ok(response)
    } else {
        Err(RecoveryError::Arithmetic(
            "gamma batch response is unrepresentable",
        ))
    }
}

#[derive(Clone, Copy)]
struct Integrand {
    tanks: f64,
    log_shape: f64,
    log_scale_time: f64,
    log_factorial: f64,
    converted: bool,
}

impl Integrand {
    fn log_value_and_slope(self, x: f64) -> (f64, f64) {
        let y = det::exp(x);
        let z = self.log_scale_time + x;
        let log_exposure = self.log_shape + log_softplus(z);
        let exposure = det::exp(log_exposure);
        let base = self.tanks * x - y - self.log_factorial;
        if !self.converted {
            let slope = self.tanks - y - det::exp(self.log_shape - softplus(-z));
            return (base - exposure, slope);
        }
        let log_response = if log_exposure < -36.0 {
            log_exposure
        } else if exposure > 36.0 {
            log_one_plus(-det::exp(-exposure))
        } else {
            det::ln(-det::expm1(-exposure))
        };
        // d ln(1-exp(-a*ln(1+exp(z)))) / dz, written as two ratios
        // in [0,1]. Each decreases, proving log-concavity of this channel.
        let first = if z < -36.0 {
            1.0
        } else {
            det::exp(-softplus(-z) - log_softplus(z))
        };
        let second = if log_exposure < -36.0 {
            1.0
        } else if exposure > 700.0 {
            0.0
        } else {
            exposure / det::expm1(exposure)
        };
        (base + log_response, self.tanks - y + first * second)
    }

    fn mode(self) -> f64 {
        let log_n = det::ln(self.tanks);
        let (mut lo, mut hi) = if self.converted {
            (log_n, det::ln(self.tanks + 1.0))
        } else {
            (
                log_n - softplus(self.log_shape + self.log_scale_time),
                log_n,
            )
        };
        for _ in 0..96 {
            let mid = lo + (hi - lo) * 0.5;
            if mid == lo || mid == hi {
                break;
            }
            if self.log_value_and_slope(mid).1 > 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        lo + (hi - lo) * 0.5
    }
}

struct Integral {
    value: f64,
    error: f64,
    tail: f64,
}

struct Panel {
    lo: f64,
    hi: f64,
    values: [f64; 3],
    integral: f64,
    tolerance: f64,
    depth: usize,
}

fn simpson(lo: f64, hi: f64, values: [f64; 3]) -> f64 {
    (hi - lo) * (values[0] + 4.0 * values[1] + values[2]) / 6.0
}

fn integrate(model: Integrand, evaluations: &mut usize) -> Result<Integral, RecoveryError> {
    let mode = model.mode();
    let peak = model.log_value_and_slope(mode).0;
    if !peak.is_finite() {
        return Err(RecoveryError::Arithmetic(
            "gamma integral mode is unrepresentable",
        ));
    }
    let scaled = |x| det::exp(model.log_value_and_slope(x).0 - peak);
    let width = 1.0 / det::sqrt(model.tanks);
    // On each half of this interval the minimum is its non-mode endpoint.
    let lower_bound = width * (scaled(mode - width) + scaled(mode + width));
    if !lower_bound.is_finite() || lower_bound <= 0.0 {
        return Err(RecoveryError::Arithmetic("gamma integral lower bound"));
    }
    let tail_target = RELATIVE_TAIL_TARGET * lower_bound;
    let mut endpoints = [0.0; 2];
    let mut tail = 0.0;
    for (index, direction) in [-1.0, 1.0].into_iter().enumerate() {
        let mut distance = width;
        let mut found = false;
        for _ in 0..32 {
            let endpoint = mode + direction * distance;
            let (value, slope) = model.log_value_and_slope(endpoint);
            if direction * slope < 0.0 {
                let bound = det::exp(value - peak) / slope.abs();
                if bound.is_finite() && bound <= tail_target {
                    endpoints[index] = endpoint;
                    tail += bound;
                    found = true;
                    break;
                }
            }
            distance *= 2.0;
        }
        if !found {
            return Err(RecoveryError::Arithmetic(
                "gamma integral tail could not be bounded",
            ));
        }
    }
    let mut evaluate = |x| -> Result<f64, RecoveryError> {
        if *evaluations >= MAX_EVALUATIONS {
            return Err(RecoveryError::Arithmetic("gamma integration work limit"));
        }
        *evaluations += 1;
        let value = scaled(x);
        if value.is_finite() && value >= 0.0 {
            Ok(value)
        } else {
            Err(RecoveryError::Arithmetic("gamma integrand"))
        }
    };
    let tolerance = RELATIVE_QUADRATURE_TARGET * lower_bound;
    let mut panels = Vec::new();
    // The mode is an explicit endpoint: a wide interval cannot miss its peak.
    for (lo, hi) in [(endpoints[0], mode), (mode, endpoints[1])] {
        let values = [
            evaluate(lo)?,
            evaluate(lo + (hi - lo) * 0.5)?,
            evaluate(hi)?,
        ];
        panels.push(Panel {
            lo,
            hi,
            values,
            integral: simpson(lo, hi, values),
            tolerance: tolerance * 0.5,
            depth: 0,
        });
    }
    let (mut integral, mut error) = (0.0, 0.0);
    while let Some(panel) = panels.pop() {
        let mid = panel.lo + (panel.hi - panel.lo) * 0.5;
        let left_mid = panel.lo + (mid - panel.lo) * 0.5;
        let right_mid = mid + (panel.hi - mid) * 0.5;
        if left_mid == panel.lo || right_mid == panel.hi {
            return Err(RecoveryError::Arithmetic(
                "gamma quadrature interval is unresolved",
            ));
        }
        let left_values = [panel.values[0], evaluate(left_mid)?, panel.values[1]];
        let right_values = [panel.values[1], evaluate(right_mid)?, panel.values[2]];
        let left = simpson(panel.lo, mid, left_values);
        let right = simpson(mid, panel.hi, right_values);
        let difference = left + right - panel.integral;
        let estimate = difference.abs() / 15.0;
        // Subdivision cannot resolve fluctuations below the floating-point
        // conditioning of log(f)-log(f_mode), especially for tiny probabilities.
        // Admit that resolution explicitly and retain it in the reported error.
        let resolution = 16.0 * f64::EPSILON * (1.0 + peak.abs()) * (left + right);
        if estimate <= panel.tolerance + resolution {
            let candidate = left + right + difference / 15.0;
            if candidate < 0.0 || !candidate.is_finite() {
                return Err(RecoveryError::Arithmetic("gamma quadrature contribution"));
            }
            integral += candidate;
            error += estimate;
        } else {
            if panel.depth >= MAX_DEPTH {
                return Err(RecoveryError::Arithmetic(
                    "gamma quadrature did not converge",
                ));
            }
            for (lo, hi, values, part) in [
                (mid, panel.hi, right_values, right),
                (panel.lo, mid, left_values, left),
            ] {
                panels.push(Panel {
                    lo,
                    hi,
                    values,
                    integral: part,
                    tolerance: panel.tolerance * 0.5,
                    depth: panel.depth + 1,
                });
            }
        }
    }
    let value = det::exp(peak + det::ln(integral));
    if !value.is_finite() || value <= 0.0 {
        return Err(RecoveryError::Arithmetic(
            "positive gamma fraction is unrepresentable",
        ));
    }
    let rounding = 64.0 * f64::EPSILON * (1.0 + peak.abs());
    Ok(Integral {
        value,
        error: value * ((error + tail) / integral + rounding),
        tail: value * (tail / integral),
    })
}

pub(super) fn continuous(
    shape: f64,
    rate: f64,
    mean: f64,
    tanks: usize,
) -> Result<GammaFractions, RecoveryError> {
    let log_scale_time = det::ln(rate) + det::ln(mean) - det::ln(tanks as f64);
    let log_factorial = (1..tanks).map(|j| det::ln(j as f64)).sum();
    let mut model = Integrand {
        tanks: tanks as f64,
        log_shape: det::ln(shape),
        log_scale_time,
        log_factorial,
        converted: true,
    };
    let mut evaluations = 0;
    let recovered = integrate(model, &mut evaluations)?;
    model.converted = false;
    let surviving = integrate(model, &mut evaluations)?;
    let residual = recovered.value + surviving.value - 1.0;
    let error = recovered.error + surviving.error;
    if !residual.is_finite() || residual.abs() > error + 2.0e-13 {
        return Err(RecoveryError::Arithmetic(
            "independent gamma integral closure",
        ));
    }
    // R+S=1 is the defining probability identity, not a fitted stream balance.
    // Evaluate the larger probability by complementing the independently
    // resolved smaller one; this preserves tiny recoveries AND surviving tails.
    let (converted, remaining) = if recovered.value <= surviving.value {
        (recovered.value, 1.0 - recovered.value)
    } else {
        (1.0 - surviving.value, surviving.value)
    };
    if converted <= 0.0 || converted > 1.0 || remaining <= 0.0 || remaining > 1.0 {
        return Err(RecoveryError::Arithmetic(
            "gamma probabilities outside (0,1]",
        ));
    }
    Ok(GammaFractions {
        converted,
        surviving: remaining,
        report: GammaIntegrationReport {
            absolute_error_estimate: Dimensionless::new(error + residual.abs()),
            omitted_tail_bound: Dimensionless::new(recovered.tail + surviving.tail),
            independent_closure_residual: Dimensionless::new(residual),
            evaluations,
        },
    })
}
