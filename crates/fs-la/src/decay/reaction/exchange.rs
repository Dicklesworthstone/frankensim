//! Exact conservative two-state exchange, including integrated event histories.
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
