//! Scalar stress-offset owner for a homogeneous regularized power-law passage.
//! All through-gap integrals are analytic. A short convergent power series
//! evaluates near-coincident endpoint divided differences without subtraction;
//! it is not a spatial quadrature or an artificial shear-rate floor.

use super::super::{Channel, CompiledLaw, PowerLawFilm, RegularizedPowerLaw, positive_or_log};
use super::{FilmError, Motion, WallResponse, fail};
use fs_math::det;

struct Law {
    material: RegularizedPowerLaw,
    compiled: CompiledLaw,
    core_viscosity: f64,
}

impl Law {
    fn rate(&self, stress: f64) -> Result<f64, FilmError> {
        if stress == 0.0 {
            return Ok(0.0);
        }
        if self.compiled.newtonian || stress <= self.compiled.transition_stress_pa {
            return positive_or_log(stress / self.core_viscosity, || {
                det::ln(stress) - det::ln(self.core_viscosity)
            });
        }
        let log_rate = self.compiled.log_minimum_shear
            + self.compiled.inverse_n * (det::ln(stress) - self.compiled.log_transition_stress);
        positive_or_log(
            self.material.minimum_shear_rate_per_s
                * det::pow(
                    stress / self.compiled.transition_stress_pa,
                    self.compiled.inverse_n,
                ),
            || log_rate,
        )
    }

    fn differential(&self, stress: f64) -> Result<f64, FilmError> {
        if self.compiled.newtonian || stress <= self.compiled.transition_stress_pa {
            return positive_or_log(1.0 / self.core_viscosity, || -det::ln(self.core_viscosity));
        }
        let rate = self.rate(stress)?;
        positive_or_log(self.compiled.inverse_n * (rate / stress), || {
            det::ln(self.compiled.inverse_n) + det::ln(rate) - det::ln(stress)
        })
    }

    fn stress(&self, rate: f64) -> Result<f64, FilmError> {
        if rate == 0.0 {
            return Ok(0.0);
        }
        if self.compiled.newtonian || rate <= self.material.minimum_shear_rate_per_s {
            return positive_or_log(self.core_viscosity * rate, || {
                det::ln(self.core_viscosity) + det::ln(rate)
            });
        }
        positive_or_log(
            self.material.consistency_pa_s_n * det::pow(rate, self.material.flow_index),
            || det::ln(self.material.consistency_pa_s_n) + self.material.flow_index * det::ln(rate),
        )
    }

    // Average g or g' on a POSITIVE stress interval. `span` retains a tiny
    // interval even if subtracting the rounded endpoint stresses loses it.
    fn average(&self, low: f64, high: f64, span: f64, derivative: bool) -> Result<f64, FilmError> {
        let core = self.compiled.transition_stress_pa;
        if self.compiled.newtonian || high <= core {
            return if derivative {
                self.differential(0.0)
            } else {
                self.rate(high - 0.5 * span)
            };
        }
        if low < core {
            let fraction = (core - low) / span;
            let inner = if derivative {
                self.differential(0.0)?
            } else {
                self.rate(0.5 * (low + core))?
            };
            let outer = self.average(core, high, high - core, derivative)?;
            let result = fraction * inner + (1.0 - fraction) * outer;
            return positive_or_log(result, || f64::NAN);
        }
        let exponent = self.compiled.inverse_n + if derivative { 0.0 } else { 1.0 };
        let distribution = power_moments(exponent, low, high, span)?;
        let endpoint = if derivative {
            self.differential(high)?
        } else {
            self.rate(high)?
        };
        positive_or_log(endpoint * distribution.average, || {
            det::ln(endpoint) + det::ln(distribution.average)
        })
    }

    // c and s are nonnegative here; odd symmetry is restored by the caller.
    // Cancelling the symmetric part of the stress interval BEFORE integrating
    // keeps arbitrarily small relative sliding in pressure-dominated flow.
    fn mean_rate(&self, c: f64, s: f64) -> Result<(f64, f64), FilmError> {
        if c == 0.0 {
            return Ok((0.0, self.rate(s)? / s));
        }
        let high = c + s;
        let low = (c - s).abs();
        let span = 2.0 * c.min(s);
        if !high.is_finite() || !span.is_finite() || span <= 0.0 {
            return fail("power-law wall stress interval is unrepresentable");
        }
        let mean = self.average(low, high, span, false)? * (c / s).min(1.0);
        let derivative = if c <= s {
            (0.5 * self.rate(high)? + 0.5 * self.rate(low)?) / s
        } else {
            self.average(low, high, span, true)?
        };
        if !mean.is_finite() || !derivative.is_finite() || derivative <= 0.0 {
            return fail("power-law wall constraint differential is unrepresentable");
        }
        Ok((mean, derivative))
    }

    fn offset(
        &self,
        s: f64,
        target: f64,
        maximum_iterations: usize,
        check: &mut impl FnMut() -> Result<(), FilmError>,
    ) -> Result<f64, FilmError> {
        let pure = self.stress(target)?;
        let mut lower = (pure - s).max(0.0);
        let mut upper = pure + s;
        if !upper.is_finite() {
            return fail("power-law wall stress bracket is unrepresentable");
        }
        let mut c = pure;
        for _ in 0..maximum_iterations {
            check()?;
            let (mean, derivative) = self.mean_rate(c, s)?;
            let residual = mean - target;
            if residual.abs() <= 16.0 * f64::EPSILON * mean.max(target) {
                return Ok(c);
            }
            if residual > 0.0 {
                upper = c;
            } else {
                lower = c;
            }
            let newton = c - residual / derivative;
            let next = if newton > lower && newton < upper {
                newton
            } else {
                0.5 * lower + 0.5 * upper
            };
            if next == c {
                return fail("power-law wall stress root stalled before velocity closure");
            }
            c = next;
        }
        fail("power-law wall stress iteration budget exhausted")
    }
}

// Moments of density (1-d*u)^(a-1) on u in [0,1]. The small-span branch
// retains the centered first moment directly, including tiny asymmetry.
struct Distribution {
    average: f64,
    mean: f64,
    second: f64,
    centered: f64,
    variance: f64,
}

fn power_moments(a: f64, low: f64, high: f64, span: f64) -> Result<Distribution, FilmError> {
    let d = span / high;
    let p = a - 1.0;
    let (average, mean, second, centered) = if d < 0.25 && p.abs() * d <= 0.5 {
        let (mut f0, mut f1, mut f2, mut fc) = (1.0, 0.5, 1.0 / 3.0, 0.0);
        let mut coefficient = 1.0;
        for k in 1..=48 {
            let kf = f64::from(k);
            coefficient *= (kf - 1.0 - p) * d / kf;
            f0 += coefficient / (kf + 1.0);
            f1 += coefficient / (kf + 2.0);
            f2 += coefficient / (kf + 3.0);
            fc += coefficient * kf / (2.0 * (kf + 1.0) * (kf + 2.0));
            if coefficient == 0.0 || (k >= 4 && coefficient.abs() < f64::EPSILON * f0.abs()) {
                break;
            }
        }
        (f0, f1 / f0, f2 / f0, fc / f0)
    } else {
        let r = 1.0 - d;
        let logarithm = if r > 0.0 && r < 1.0 {
            det::ln(r) * (d / (1.0 - r))
        } else if r == 1.0 {
            -d
        } else {
            det::ln(low) - det::ln(high)
        };
        let t = a * logarithm;
        let denominator = -det::expm1(t);
        let tail = det::exp(t) / denominator;
        let inverse1 = 1.0 / ((a + 1.0) * d);
        let inverse2 = 1.0 / ((a + 2.0) * d);
        let correction = (a * d) * tail;
        let m1 = (1.0 - correction) * inverse1;
        let m2 = (2.0 - correction * (2.0 + (a + 1.0) * d)) * inverse1 * inverse2;
        (denominator / (a * d), m1, m2, m1 - 0.5)
    };
    let variance = second - mean * mean;
    if [average, mean, second, centered, variance]
        .iter()
        .any(|v| !v.is_finite())
        || average <= 0.0
        || !(0.0..1.0).contains(&mean)
        || variance <= 0.0
        || second <= 0.0
        || second >= mean
    {
        return fail("power-law through-gap stress moments are unrepresentable");
    }
    Ok(Distribution {
        average,
        mean,
        second,
        centered,
        variance,
    })
}

// Unnormalized positive weighted variance is accumulated about local means,
// so near-wall concentration does not require J2 - J1*J1/J0 subtraction.
#[derive(Default)]
struct Moments {
    weight: f64,
    mean: f64,
    variance: f64,
    wall_weight: f64,
}

impl Moments {
    fn add(&mut self, weight: f64, mean: f64, variance: f64, wall_weight: f64) {
        let total = self.weight + weight;
        let delta = mean - self.mean;
        self.variance += weight * variance + (self.weight / total) * weight * delta * delta;
        self.mean += (weight / total) * delta;
        self.weight = total;
        self.wall_weight += weight * wall_weight;
    }
}

fn moments(law: &Law, c: f64, s: f64) -> Result<Moments, FilmError> {
    let core = law.compiled.transition_stress_pa;
    let (left, right) = (c + s, c - s);
    let mut cuts = [-1.0, 1.0, 1.0, 1.0];
    let mut count = 2;
    for transition in [core, -core] {
        if right < transition && transition < left {
            cuts[count] = (c - transition) / s;
            count += 1;
        }
    }
    cuts[..count].sort_by(f64::total_cmp);
    let mut result = Moments::default();
    for interval in cuts[..count].windows(2) {
        let (x0, x1) = (interval[0], interval[1]);
        let width = x1 - x0;
        if width == 0.0 {
            continue;
        }
        let middle = 0.5 * x0 + 0.5 * x1;
        let (t0, t1) = ((c - s * x0).abs(), (c - s * x1).abs());
        let (start, end, density, distribution) = if (c - s * middle).abs() <= core {
            (
                x0,
                x1,
                1.0 / law.core_viscosity,
                Distribution {
                    average: 1.0,
                    mean: 0.5,
                    second: 1.0 / 3.0,
                    centered: 0.0,
                    variance: 1.0 / 12.0,
                },
            )
        } else {
            let (start, end, high, low) = if t0 >= t1 {
                (x0, x1, t0, t1)
            } else {
                (x1, x0, t1, t0)
            };
            (
                start,
                end,
                law.differential(high)?,
                power_moments(law.compiled.inverse_n, low, high, s * width)?,
            )
        };
        let weight = (0.5 * width) * density * distribution.average;
        let mean = middle + (end - start) * distribution.centered;
        let variance = width * width * distribution.variance;
        // Bernstein form of (1-x)*(1+x) has only nonnegative terms.
        let u = distribution.mean;
        let uu = distribution.second;
        let wall_weight = (1.0 - start) * (1.0 + start) * (1.0 - 2.0 * u + uu)
            + ((1.0 - start) * (1.0 + end) + (1.0 - end) * (1.0 + start)) * (u - uu)
            + (1.0 - end) * (1.0 + end) * uu;
        if !weight.is_finite() || weight <= 0.0 || !wall_weight.is_finite() || wall_weight <= 0.0 {
            return fail("power-law differential stress weight is unrepresentable");
        }
        result.add(weight, mean, variance, wall_weight);
    }
    // Avoid cancellation of two symmetric weighted means for very small
    // relative sliding under a much larger pressure-induced stress.
    if c <= 0.5 * s {
        let (mean_rate, _) = law.mean_rate(c, s)?;
        let endpoint_mean = if c == 0.0 {
            0.0
        } else {
            c * law.average(s - c, s + c, 2.0 * c, true)?
        };
        result.mean = ((mean_rate - endpoint_mean) / s) / result.weight;
    }
    if [
        result.weight,
        result.mean,
        result.variance,
        result.wall_weight,
    ]
    .iter()
    .any(|v| !v.is_finite())
        || result.weight <= 0.0
        || result.variance <= 0.0
    {
        return fail("power-law differential stress variance is unrepresentable");
    }
    Ok(result)
}

#[allow(clippy::too_many_arguments)]
pub(in super::super) fn evaluate(
    film: &PowerLawFilm,
    index: usize,
    gap: f64,
    core_conductance: f64,
    dp: f64,
    motion: Motion,
    check: &mut impl FnMut() -> Result<(), FilmError>,
) -> Result<Channel, FilmError> {
    let edge = &film.film.channels[index];
    let law = Law {
        material: film.laws[index],
        compiled: film.compiled[index],
        core_viscosity: film.film.channel_viscosities_pa_s()[index],
    };
    let rate = motion.relative / gap;
    let s = dp * (0.5 * gap / edge.length_m);
    if !rate.is_finite() || !s.is_finite() || (dp != 0.0 && s == 0.0) {
        return fail("power-law moving-wall shear rate or stress is unrepresentable");
    }
    let pure_stress = law.stress(rate.abs())?.copysign(rate);
    let (c, flow, differential, aperture, relative_flow, cr) = if law.compiled.newtonian
        || pure_stress.abs() + s.abs() <= law.compiled.transition_stress_pa
    {
        let fixed = film.channel(index, gap, core_conductance, dp)?;
        let drag = law.core_viscosity / gap;
        (
            pure_stress,
            fixed.flow,
            fixed.differential,
            fixed.aperture_derivative,
            0.0,
            drag,
        )
    } else if motion.relative == 0.0 {
        let fixed = film.channel(index, gap, core_conductance, dp)?;
        let inverse_drag = law.rate(s.abs())? / s.abs();
        (
            0.0,
            fixed.flow,
            fixed.differential,
            fixed.aperture_derivative,
            0.0,
            1.0 / (gap * inverse_drag),
        )
    } else if dp == 0.0 {
        let viscosity = law.material.flow_index * (pure_stress / rate);
        let differential = (edge.width_m / edge.length_m) * gap * gap * gap / (12.0 * viscosity);
        (pure_stress, 0.0, differential, 0.0, 0.0, viscosity / gap)
    } else {
        let positive = law.offset(s.abs(), rate.abs(), film.controls.maximum_iterations, check)?;
        let c = positive.copysign(rate);
        let moments = moments(&law, positive, s.abs())?;
        let mean = moments.mean * rate.signum() * s.signum();
        let flow = (0.25 * edge.width_m * gap * gap) * s * moments.wall_weight;
        let differential =
            (0.25 * edge.width_m / edge.length_m) * gap * gap * gap * moments.variance;
        let aperture =
            2.0 * flow / gap + (0.5 * edge.width_m * gap) * (s * moments.variance + mean * rate);
        (
            c,
            flow,
            differential,
            aperture,
            -(0.5 * edge.width_m * gap) * mean,
            1.0 / (gap * moments.weight),
        )
    };
    let mean_effort = -(edge.width_m * gap) * dp;
    let relative_effort = (edge.width_m * edge.length_m) * c;
    let wall = WallResponse {
        mean_effort,
        relative_effort,
        relative_flow_derivative: relative_flow,
        relative_effort_derivative: (edge.width_m * edge.length_m) * cr,
        pressure_flow: flow,
        mid_stress: c,
        tangent_ambiguous: !law.compiled.newtonian
            && dp == 0.0
            && rate.abs() == law.material.minimum_shear_rate_per_s,
    };
    let result = Channel {
        flow: flow + motion.entrainment,
        differential,
        // Keep wall-transport derivatives separate until their incidence sum
        // cancels; otherwise large common translation can erase pressure terms.
        aperture_derivative: aperture,
        wall,
    };
    if [
        result.flow,
        differential,
        result.aperture_derivative,
        relative_flow,
        cr,
        mean_effort,
        relative_effort,
        -(c + s),
        c - s,
        wall.relative_effort_derivative,
        flow * dp,
        relative_effort * motion.relative,
        mean_effort * motion.mean + relative_effort * motion.relative,
    ]
    .iter()
    .any(|v| !v.is_finite())
        || differential <= 0.0
        || cr <= 0.0
        || flow * dp < 0.0
        || relative_effort * motion.relative < 0.0
    {
        return fail("power-law moving-wall response or work is unrepresentable");
    }
    Ok(result)
}
