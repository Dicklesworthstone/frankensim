//! Exact two-state exchange, including sinks and integrated event histories.
//!
//! For rates a:0->1 and b:1->0, P=S I+(1-S) Pi, where S=exp(-(a+b)t)
//! and each column of Pi is [b,a]/(a+b). The first and second normalized
//! integrals replace S by phi1 and phi2 and the stationary coefficient by
//! 1-phi1 and 1/2-phi2. Positive terms retain tiny transfers and source ages.

use super::{
    DecayError, DecayLimits, ReactionPlan, Sum, charge, det, poll, positive_or_zero, source,
};

pub(super) fn prepare(
    mut plan: ReactionPlan,
    limits: DecayLimits,
    with_source: bool,
    with_residence: bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(ReactionPlan, Option<source::SourceMaps>), DecayError> {
    // Fixed scalar work, matrix writes and both independent column audits.
    // This charge is independent of the number of exchange lifetimes.
    charge(
        &mut plan.work,
        384 + 128 * usize::from(with_source) + 32 * usize::from(with_residence),
        limits.maximum_work,
    )?;
    poll(cancelled)?;
    let [a, b] = [plan.rates[0], plan.rates[1]];
    let time = plan.duration_s;
    let scale = a.max(b);
    let scaled_sum = a / scale + b / scale;
    let wa = (a / scale) / scaled_sum;
    let wb = (b / scale) / scaled_sum;
    if wa == 0.0 || wb == 0.0 {
        return Err(DecayError::Unrepresentable("two-state exchange rate ratio"));
    }
    // Do not form a+b: that sum may overflow even for a finite exposure.
    let exposure = (scale * time) * scaled_sum;
    if exposure == 0.0 {
        return Err(DecayError::Unrepresentable("two-state exchange exposure"));
    }
    let survival = det::exp(-exposure);
    let loss = -det::expm1(-exposure);
    let (phi1, mean_loss, phi2, second_loss) = moments(exposure, scale, scaled_sum, time);
    let forward = wa * loss;
    let backward = wb * loss;
    if forward == 0.0 || backward == 0.0 {
        return Err(DecayError::Unrepresentable(
            "positive two-state exchange transfer",
        ));
    }
    plan.transition = vec![survival + backward, backward, forward, survival + forward];
    plan.reachable.fill(true);

    // k=a*b/(a+b), integrated without an overflowing product or underflowing
    // a subnormal k before multiplying it by a long duration. The same positive
    // cycle term belongs to BOTH event accounts. The additional direct term
    // retains the initial relaxation event even after millions of revisits.
    let (cycle_rate, cycle_weight) = if a <= b { (a, wb) } else { (b, wa) };
    let repeated = integrated_cycle(cycle_rate, time, cycle_weight, mean_loss)?;
    plan.integrated = vec![repeated + forward, repeated, repeated, repeated + backward];
    poll(cancelled)?;
    let maps = if with_source {
        let source_forward = wa * mean_loss;
        let source_backward = wb * mean_loss;
        let source_repeated = integrated_cycle(cycle_rate, time, cycle_weight, second_loss)?;
        let mut maps = source::SourceMaps::identity(2, with_residence);
        maps.population = vec![
            phi1 + source_backward,
            source_backward,
            source_forward,
            phi1 + source_forward,
        ];
        maps.removals = vec![
            source_repeated + source_forward,
            source_repeated,
            source_repeated,
            source_repeated + source_backward,
        ];
        if let Some(residence) = &mut maps.residence {
            let forward = wa * second_loss;
            let backward = wb * second_loss;
            residence.population = vec![phi2 + backward, backward, forward, phi2 + forward];
        }
        Some(maps)
    } else {
        None
    };
    for column in 0..2 {
        poll(cancelled)?;
        let initial = if column == 0 { [1.0, 0.0] } else { [0.0, 1.0] };
        let endpoint = [plan.transition[column], plan.transition[2 + column]];
        let removals = [plan.integrated[column], plan.integrated[2 + column]];
        plan.check_balance(&initial, None, &endpoint, &removals, column)?;
        if let Some(maps) = &maps {
            let endpoint = [maps.population[column], maps.population[2 + column]];
            let removals = [maps.removals[column], maps.removals[2 + column]];
            plan.check_balance(&initial, None, &endpoint, &removals, column)?;
        }
    }
    poll(cancelled)?;
    Ok((plan, maps))
}

/// Positive, independently integrated maps for a passive pair with separated
/// removal rates. Every coefficient is representable before replacing the
/// general owner; exhausted tails and exceptional scales retain its admission.
pub(super) struct DampedPair {
    transition: [f64; 4],
    integrated: [f64; 4],
    population: [f64; 4],
    removals: [f64; 4],
    residence: [f64; 4],
}

pub(super) fn damped_pair(
    plan: &mut ReactionPlan,
    limits: DecayLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Option<DampedPair>, DecayError> {
    if plan.nodes != 2 || plan.yields[0] != 0.0 || plan.yields[3] != 0.0 {
        return Ok(None);
    }
    let [a, b] = [plan.rates[0], plan.rates[1]];
    let [p, q] = [plan.yields[2], plan.yields[1]];
    if a == 0.0 || b == 0.0 || p <= 0.0 || q <= 0.0 || p > 1.0 || q > 1.0 {
        return Ok(None);
    }
    let scale = a.max(b);
    let small = a.min(b);
    let ratio = small / scale;
    let exposure = scale * plan.duration_s;
    // With separated rates and at least one fast lifetime, the two scalar
    // moment differences below are well conditioned. Near-degenerate modes
    // and tiny exposures continue through the original positive series.
    if ratio > 0.5 || exposure < 1.0 || !exposure.is_finite() {
        return Ok(None);
    }
    charge(&mut plan.work, 768, limits.maximum_work)?;
    poll(cancelled)?;
    let d = 1.0 - ratio;
    let discriminant = det::sqrt(d * d + 4.0 * ratio * p * q);
    // Rationalize the small root correction instead of subtracting two
    // nearly equal quantities. Neither a+b nor a*b needs to be representable.
    let correction = (2.0 * ratio * p * q) / (discriminant + d);
    let mut sink = Sum::default();
    sink.add(1.0 - p);
    sink.add(p * (1.0 - q));
    let slow = ((small * plan.duration_s) * sink.total()) / (1.0 + correction);
    let fast = exposure * (1.0 + correction);
    let gap = exposure * discriminant;
    if !fast.is_finite() || !gap.is_finite() {
        return Ok(None);
    }
    let slow_weight = correction / discriminant;
    let fast_weight = (d + correction) / discriminant;
    let (wa, wb) = if a >= b {
        (slow_weight, fast_weight)
    } else {
        (fast_weight, slow_weight)
    };
    let backward = ((b / scale) * q) / discriminant;
    let forward = ((a / scale) * p) / discriminant;
    let es = det::exp(-slow);
    let ef = det::exp(-fast);
    let difference = es * -det::expm1(-gap);
    let (first_s, _, second_s, _) = moments(slow, 1.0, 1.0, 1.0);
    let (first_f, _, second_f, _) = moments(fast, 1.0, 1.0, 1.0);
    let transition = [
        wa * es + wb * ef,
        backward * difference,
        forward * difference,
        wb * es + wa * ef,
    ];
    let population = [
        wa * first_s + wb * first_f,
        backward * (first_s - first_f),
        forward * (first_s - first_f),
        wb * first_s + wa * first_f,
    ];
    let residence = [
        wa * second_s + wb * second_f,
        backward * (second_s - second_f),
        forward * (second_s - second_f),
        wb * second_s + wa * second_f,
    ];
    let mut integrated = [0.0; 4];
    let mut removals = [0.0; 4];
    for i in 0..2 {
        let exposure = plan.rates[i] * plan.duration_s;
        for j in 0..2 {
            integrated[2 * i + j] = exposure * population[2 * i + j];
            removals[2 * i + j] = exposure * residence[2 * i + j];
        }
    }
    // This fast path never grants an underflow exemption or manufactures a
    // zero response. The general owner retains its existing tail accounting.
    if !transition
        .iter()
        .chain(&integrated)
        .chain(&population)
        .chain(&removals)
        .chain(&residence)
        .all(|&value| value.is_finite() && value > 0.0)
    {
        return Ok(None);
    }
    poll(cancelled)?;
    Ok(Some(DampedPair {
        transition,
        integrated,
        population,
        removals,
        residence,
    }))
}

impl DampedPair {
    pub(super) fn finish(
        self,
        mut plan: ReactionPlan,
        with_source: bool,
        with_residence: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(ReactionPlan, Option<source::SourceMaps>), DecayError> {
        plan.transition = self.transition.to_vec();
        plan.integrated = self.integrated.to_vec();
        plan.reachable.fill(true);
        let maps = with_source.then(|| {
            let mut maps = source::SourceMaps::identity(2, with_residence);
            maps.population = self.population.to_vec();
            maps.removals = self.removals.to_vec();
            if let Some(residence) = &mut maps.residence {
                residence.population = self.residence.to_vec();
            }
            maps
        });
        for column in 0..2 {
            poll(cancelled)?;
            let initial = if column == 0 { [1.0, 0.0] } else { [0.0, 1.0] };
            let endpoint = [plan.transition[column], plan.transition[2 + column]];
            let removals = [plan.integrated[column], plan.integrated[2 + column]];
            plan.check_balance(&initial, None, &endpoint, &removals, column)?;
            if let Some(maps) = &maps {
                let endpoint = [maps.population[column], maps.population[2 + column]];
                let removals = [maps.removals[column], maps.removals[2 + column]];
                plan.check_balance(&initial, None, &endpoint, &removals, column)?;
            }
        }
        poll(cancelled)?;
        Ok((plan, maps))
    }
}

// phi1=integral_0^1 exp(-x*s)ds; phi2=integral_0^1 (1-s)exp(-x*s)ds.
// Evaluate their small complements independently. In particular 1/2-phi2
// cannot be formed by subtraction when x is tiny: it owns repeated source
// events and cross-state residence, whose first nonzero term is x/6.
fn moments(x: f64, scale: f64, scaled_sum: f64, time: f64) -> (f64, f64, f64, f64) {
    if x <= 0.5 {
        let mut phi2 = Sum::default();
        let mut phi3 = Sum::default();
        let mut term2 = 0.5;
        let mut term3 = 1.0 / 6.0;
        phi2.add(term2);
        phi3.add(term3);
        // At x<=1/2 the first omitted terms are below 1e-22. This bounded
        // scalar series is not the network's uniformization series.
        for k in 1..=16 {
            term2 *= -x / (k + 2) as f64;
            term3 *= -x / (k + 3) as f64;
            phi2.add(term2);
            phi3.add(term3);
        }
        let phi2 = phi2.total();
        let mean_loss = x * phi2;
        (1.0 - mean_loss, mean_loss, phi2, x * phi3.total())
    } else {
        let inverse = if x.is_infinite() {
            ((1.0 / scale) / scaled_sum) / time
        } else {
            1.0 / x
        };
        let phi1 = -det::expm1(-x) * inverse;
        let mean_loss = 1.0 - phi1;
        let phi2 = mean_loss * inverse;
        (phi1, mean_loss, phi2, 0.5 - phi2)
    }
}

fn integrated_cycle(rate: f64, time: f64, weight: f64, factor: f64) -> Result<f64, DecayError> {
    let exposure = rate * time;
    let value = if exposure.is_finite() {
        (exposure * weight) * factor
    } else {
        // A bounded factor can make the complete count finite even when
        // rate*time alone overflows. No infinite intermediate is published.
        (rate * (weight * factor)) * time
    };
    if !positive_or_zero(value) {
        return Err(DecayError::Unrepresentable(
            "two-state exchange event count",
        ));
    }
    Ok(value)
}

#[cfg(test)]
#[path = "exchange/tests.rs"]
mod tests;
