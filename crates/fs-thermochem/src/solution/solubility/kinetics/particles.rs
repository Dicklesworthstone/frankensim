//! The separable closed-bath radius equation for a retained particle population.
//!
//! Writing s = r³ removes the vanishing area at exhaustion:
//! r' = -G_ref (r³ + B)/(3 V s_ref^(2/3)), B = capacity - total.
//! For B != 0 the radius is scaled by |B|^(1/3). The two rational
//! primitives are evaluated with their logarithmic singularity or far-field
//! tail isolated. A bounded safeguarded scalar solve finds elapsed time.
//! Small intervals use fixed eight-point Gauss quadrature of a smooth scaled
//! integrand to avoid subtracting almost equal primitive values. No area is
//! frozen during an accepted interval, and no explicit time-step clamp is used.

use super::{
    SolubilityErrorV1, SolubilityKineticStateV1, SolubilityKineticStepV1, SolubilityKineticsV1,
    product_ratio,
};
use fs_math::decay::{attenuate, decayed_fraction};
use fs_math::det::{atan2, exp, ln, pow};
use fs_qty::{Amount, Time, Volume};

const SQRT_3: f64 = 1.732_050_807_568_877_2;
const MAX_ITERATIONS: usize = 128;

pub(super) fn advance(
    model: &SolubilityKineticsV1,
    _reference_solid: Amount,
    solvent_volume: Volume,
    initial: &SolubilityKineticStateV1,
    duration: Time,
    mut result: SolubilityKineticStepV1,
) -> Result<SolubilityKineticStepV1, SolubilityErrorV1> {
    let solid = initial.solid.value();
    let dissolved = initial.dissolved.value();
    let capacity = result.equilibrium.capacity.value();
    let excess = capacity - initial.total_solute.value();

    if excess == 0.0 {
        // At capacity == total the finite-bath tail is algebraic, not
        // exponential. Its positive solid remains meaningful even when it is
        // below one ulp of the already dissolved stock.
        let g = model.conductance_at(initial.solid())?.value();
        let exposure = product_ratio(g, duration.value(), solvent_volume.value())?;
        let log_denominator = if exposure.is_finite() {
            log1p((2.0 / 3.0) * exposure)
        } else {
            ln(g) + ln(duration.value()) - ln(solvent_volume.value()) + ln(2.0 / 3.0)
        };
        let hazard = 1.5 * log_denominator;
        let remaining = attenuate(solid, hazard);
        let transfer = solid * decayed_fraction(hazard);
        return finish(initial, dissolved + transfer, remaining, transfer, result);
    }

    let scale = excess.abs();
    let ratio = solid / scale;
    let radius = if ratio.is_finite() && ratio > 0.0 {
        pow(ratio, 1.0 / 3.0)
    } else {
        exp((ln(solid) - ln(scale)) / 3.0)
    };
    if !radius.is_finite() || radius <= 0.0 {
        return Err(range("positive particle radius scale"));
    }
    let g_scale = model.conductance_at(Amount::new(scale))?.value();
    let elapsed = product_ratio(g_scale, duration.value(), solvent_volume.value())? / 3.0;
    if elapsed == 0.0 || elapsed.is_nan() {
        return Err(range("positive scaled particle interval"));
    }

    if excess > 0.0 {
        let exhaustion = positive_primitive(radius);
        if elapsed >= exhaustion {
            let event = product_ratio(solvent_volume.value(), 3.0 * exhaustion, g_scale)?;
            if !event.is_finite() || event <= 0.0 {
                return Err(range("finite particle exhaustion time"));
            }
            result.solid_exhaustion_time = Some(Time::new(event));
            return finish(initial, initial.total_solute.value(), 0.0, solid, result);
        }
        if elapsed > 0.5 * exhaustion {
            // Close to the finite event, solve for the small REMAINING radius
            // from the remaining time. Subtracting it from a radius fraction
            // near one would otherwise erase a representable residual solid.
            let remaining_time = exhaustion - elapsed;
            let endpoint = solve(
                remaining_time,
                radius,
                remaining_time,
                positive_primitive,
                |radius| {
                    if radius > 1.0 {
                        let inverse = 1.0 / radius;
                        let cube = inverse * inverse * inverse;
                        cube / (1.0 + cube)
                    } else {
                        1.0 / (1.0 + radius * radius * radius)
                    }
                },
            )?;
            let relative = endpoint / radius;
            let remaining = ((solid * relative) * relative) * relative;
            return finish(
                initial,
                dissolved + (solid - remaining),
                remaining,
                solid - remaining,
                result,
            );
        }
        // alpha is the lost fraction of the ORIGINAL radius. Keeping this
        // small change separately resolves transfer below one ulp of s itself.
        let derivative = positive_integrand(radius, 0.0);
        let alpha = solve(
            elapsed,
            1.0,
            elapsed / derivative,
            |alpha| positive_elapsed(radius, alpha),
            |alpha| positive_integrand(radius, alpha),
        )?;
        let transfer = solid * (alpha * (3.0 - alpha * (3.0 - alpha)));
        let remaining = if transfer <= 0.5 * solid {
            solid - transfer
        } else {
            let relative = 1.0 - alpha;
            ((solid * relative) * relative) * relative
        };
        return finish(initial, dissolved + transfer, remaining, transfer, result);
    }

    // B < 0: a positive solid equilibrium remains. q is the logarithmic
    // decrease in |r-r_eq|, so the asymptotic tail never crosses equilibrium.
    let driving = if dissolved > 0.5 * initial.total_solute.value() {
        excess + solid
    } else {
        capacity - dissolved
    };
    let delta = if (0.5..=2.0).contains(&radius) {
        (driving / scale) * inverse_q(radius)
    } else {
        radius - 1.0
    };
    if delta == 0.0 || !delta.is_finite() {
        return Err(range("nonzero particle radius departure"));
    }
    // This limit puts even a largest-finite amount's surviving deviation
    // below the smallest subnormal. It does not invent finite exhaustion.
    let q_limit = 800.0 + ln(driving.abs()).max(0.0);
    let limit_elapsed = negative_elapsed(radius, delta, q_limit);
    if elapsed >= limit_elapsed {
        return finish(initial, capacity, scale, driving, result);
    }
    let guess = if radius <= 1.0 {
        elapsed / inverse_q(radius)
    } else {
        let g = model.conductance_at(initial.solid())?.value();
        let exposure = product_ratio(g, duration.value(), solvent_volume.value())?;
        let inverse = 1.0 / radius;
        (exposure / 3.0) * (1.0 + inverse + inverse * inverse)
    };
    let q = solve(
        elapsed,
        q_limit,
        guess,
        |q| negative_elapsed(radius, delta, q),
        |q| inverse_q(radius_at(radius, delta, q)),
    )?;
    let endpoint_radius = radius_at(radius, delta, q);
    let fraction = decayed_fraction(q);
    let transfer = if driving > 0.0 {
        let alpha = (delta / radius) * fraction;
        solid * (alpha * (3.0 - alpha * (3.0 - alpha)))
    } else {
        let change = -delta * fraction;
        -scale
            * (change
                * (endpoint_radius * endpoint_radius + endpoint_radius * radius + radius * radius))
    };
    // Evaluate a small remaining dissolved surplus independently of transfer.
    // Subtracting a nearly exhausted dissolved donor would erase that trace.
    let q_ratio = exp(log_q(endpoint_radius) - log_q(radius));
    let remaining_driving = attenuate(driving, q) * q_ratio;
    let next_dissolved = if transfer.abs() <= 0.5 * dissolved {
        dissolved + transfer
    } else {
        capacity - remaining_driving
    };
    let next_solid = if transfer <= 0.5 * solid {
        solid - transfer
    } else {
        scale + remaining_driving
    };
    finish(initial, next_dissolved, next_solid, transfer, result)
}

fn finish(
    initial: &SolubilityKineticStateV1,
    dissolved: f64,
    solid: f64,
    transfer: f64,
    mut result: SolubilityKineticStepV1,
) -> Result<SolubilityKineticStepV1, SolubilityErrorV1> {
    if transfer == 0.0 || !transfer.is_finite() {
        return Err(range("nonzero finite particle phase transfer"));
    }
    result.state = SolubilityKineticStateV1::with_phases(initial.total_solute, dissolved, solid)?;
    result.transferred_to_solution = Amount::new(transfer);
    Ok(result)
}

fn range(message: &'static str) -> SolubilityErrorV1 {
    SolubilityErrorV1::Unrepresentable(message)
}

// Correct the rounding in 1+x before the deterministic logarithm. x >= 0.
fn log1p(x: f64) -> f64 {
    let shifted = 1.0 + x;
    ln(shifted) - ((shifted - 1.0) - x) / shifted
}

fn inverse_q(radius: f64) -> f64 {
    if radius > 1.0 {
        let inverse = 1.0 / radius;
        (inverse * inverse) / (1.0 + inverse + inverse * inverse)
    } else {
        1.0 / (1.0 + radius + radius * radius)
    }
}

fn log_q(radius: f64) -> f64 {
    if radius > 1.0 {
        let inverse = 1.0 / radius;
        2.0 * ln(radius) + log1p(inverse + inverse * inverse)
    } else {
        log1p(radius + radius * radius)
    }
}

fn radius_at(initial: f64, delta: f64, q: f64) -> f64 {
    if delta < 0.0 {
        initial + (-delta) * decayed_fraction(q)
    } else {
        1.0 + attenuate(delta, q)
    }
}

// Integral from radius to infinity of 1/(x³ +/- 1), radius >= 2.
// After x=1/u the series ratio is at most 1/8; 24 terms are ample.
fn tail(radius: f64, positive: bool) -> f64 {
    let inverse = 1.0 / radius;
    let cube = inverse * inverse * inverse;
    let ratio = if positive { -cube } else { cube };
    let mut power = 1.0;
    let mut sum = 0.5;
    for index in 1..24 {
        power *= ratio;
        sum += power / (3 * index + 2) as f64;
    }
    (inverse * inverse) * sum
}

// Integral from 0 to radius of 1/(1+x³). The angle is a difference
// evaluated as one atan2; the logarithm uses a positive log1p argument.
fn positive_primitive(radius: f64) -> f64 {
    if radius > 2.0 {
        return 2.0 * core::f64::consts::PI / (3.0 * SQRT_3) - tail(radius, true);
    }
    let denominator = radius * radius - radius + 1.0;
    log1p(3.0 * radius / denominator) / 6.0 + atan2(SQRT_3 * radius, 2.0 - radius) / SQRT_3
}

fn positive_integrand(initial: f64, alpha: f64) -> f64 {
    let radius = initial * (1.0 - alpha);
    if radius > 1.0 {
        let inverse = 1.0 / radius;
        ((initial / radius) * inverse * inverse) / (1.0 + inverse * inverse * inverse)
    } else {
        initial / (1.0 + radius * radius * radius)
    }
}

fn positive_elapsed(initial: f64, alpha: f64) -> f64 {
    if alpha <= 0.125 {
        return gauss8(alpha, |x| positive_integrand(initial, x));
    }
    let endpoint = initial * (1.0 - alpha);
    if endpoint > 2.0 {
        tail(endpoint, true) - tail(initial, true)
    } else {
        positive_primitive(initial) - positive_primitive(endpoint)
    }
}

fn negative_elapsed(initial: f64, delta: f64, q: f64) -> f64 {
    if q <= 0.125 {
        return gauss8(q, |x| inverse_q(radius_at(initial, delta, x)));
    }
    let endpoint = radius_at(initial, delta, q);
    if endpoint > 2.0 {
        return tail(endpoint, false) - tail(initial, false);
    }
    let angle = atan2(
        2.0 * SQRT_3 * (endpoint - initial),
        3.0 + (2.0 * endpoint + 1.0) * (2.0 * initial + 1.0),
    );
    q / 3.0 + (log_q(endpoint) - log_q(initial)) / 6.0 + angle / SQRT_3
}

fn gauss8(width: f64, integrand: impl Fn(f64) -> f64) -> f64 {
    let half = 0.5 * width;
    let mut sum = 0.0;
    for (node, weight) in [
        (0.183_434_642_495_649_8, 0.362_683_783_378_362),
        (0.525_532_409_916_329, 0.313_706_645_877_887_3),
        (0.796_666_477_413_626_7, 0.222_381_034_453_374_5),
        (0.960_289_856_497_536_3, 0.101_228_536_290_376_3),
    ] {
        sum += weight * (integrand(half * (1.0 - node)) + integrand(half * (1.0 + node)));
    }
    half * sum
}

fn solve(
    target: f64,
    upper: f64,
    guess: f64,
    elapsed: impl Fn(f64) -> f64,
    derivative: impl Fn(f64) -> f64,
) -> Result<f64, SolubilityErrorV1> {
    let mut low = 0.0;
    let mut high = upper;
    let mut candidate = if guess > low && guess < high {
        guess
    } else {
        0.5 * upper
    };
    for _ in 0..MAX_ITERATIONS {
        let value = elapsed(candidate);
        if !value.is_finite() || value < 0.0 {
            return Err(range("particle elapsed-time primitive"));
        }
        let residual = value - target;
        if residual.abs() <= 32.0 * f64::EPSILON * target {
            return Ok(candidate);
        }
        if residual > 0.0 {
            high = candidate;
        } else {
            low = candidate;
        }
        if high - low <= 16.0 * f64::EPSILON * candidate {
            return Ok(low + 0.5 * (high - low));
        }
        let newton = candidate - residual / derivative(candidate);
        candidate = if newton.is_finite() && newton > low && newton < high {
            newton
        } else {
            low + 0.5 * (high - low)
        };
    }
    Err(range(
        "bounded particle elapsed-time solve did not converge",
    ))
}
