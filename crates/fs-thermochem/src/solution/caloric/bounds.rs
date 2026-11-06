//! Admission and reusable storage bounds from the original caloric sources.

use super::{
    BinaryCaloricErrorV1, BinaryLiquidCaloricModelV1, BinaryRegularSolutionV1,
    PureLiquidHeatCapacityV1, inventory, positive,
};
use fs_math::det::sqrt;
use fs_qty::Temperature;

// For a negative excess coefficient C, fixed-composition Cp is positive for
// every 0<=x<=1 iff -C < (sqrt(Cp_A)+sqrt(Cp_B))². At fixed x it is affine in
// temperature between source knots; checking the union of their knots and
// the common endpoints therefore covers the entire admitted domain.
pub(super) fn admit_positive_capacity(
    pair: &BinaryRegularSolutionV1,
    sources: &[PureLiquidHeatCapacityV1; 2],
    domain: [Temperature; 2],
) -> Result<(), BinaryCaloricErrorV1> {
    let c = pair.excess_heat_capacity_coefficient_j_per_mol_k();
    if c >= 0.0 {
        return Ok(());
    }
    for temperature in source_temperatures(sources, domain) {
        let a = lower_capacity_at(&sources[0], temperature);
        let b = lower_capacity_at(&sources[1], temperature);
        // Normalize before taking roots/squaring; neither huge nor subnormal
        // source magnitudes need overflow an otherwise valid comparison.
        // Directed rounding refuses an unresolved zero-capacity boundary.
        let scale = a.max(b).max(-c);
        let lower_root = |cp: f64| sqrt((cp / scale).next_down().max(0.0)).next_down().max(0.0);
        let roots = (lower_root(a) + lower_root(b)).next_down().max(0.0);
        let threshold = (roots * roots).next_down();
        if threshold <= ((-c) / scale).next_up() {
            return Err(BinaryCaloricErrorV1::Invalid(
                "excess heat capacity permits nonpositive total Cp inside the source domain",
            ));
        }
    }
    Ok(())
}

impl BinaryLiquidCaloricModelV1 {
    /// Conservative positive lower bound on total dH/dT [J/K] for this
    /// conserved inventory and inclusive temperature interval.
    ///
    /// Original piecewise-linear pure Cp curves are evaluated at every crossed
    /// source knot and both endpoints. A negative heat capacity of mixing is
    /// included at the feed composition: equilibrium splitting decreases
    /// x_A*x_B, so it cannot make that negative contribution more negative.
    /// Equilibrium redistribution adds a nonnegative response. A positive
    /// excess coefficient can safely be omitted from the bound. This is a
    /// storage/CFL bound, not a replacement for the physical local Cp.
    pub fn isobaric_heat_capacity_lower_bound_j_per_k(
        &self,
        component_moles: [f64; 2],
        domain: [Temperature; 2],
    ) -> Result<f64, BinaryCaloricErrorV1> {
        let (total, _) = inventory(component_moles)?;
        if domain.iter().any(|t| !t.value().is_finite())
            || domain[0].value() > domain[1].value()
            || domain[0].value() < self.owner.temperature_domain[0].value()
            || domain[1].value() > self.owner.temperature_domain[1].value()
        {
            return Err(BinaryCaloricErrorV1::Invalid("heat-capacity bound domain"));
        }
        let c = self
            .thermodynamics()
            .excess_heat_capacity_coefficient_j_per_mol_k();
        let excess_lower = if c < 0.0 && component_moles.iter().all(|&n| n > 0.0) {
            // n*x_A*x_B = min(n_A,n_B)*max(n_A,n_B)/n. This order avoids
            // overflowing the product of extensive amounts or losing a trace
            // fraction before its extensive contribution is represented.
            let smaller = component_moles[0].min(component_moles[1]);
            let larger = component_moles[0].max(component_moles[1]);
            let amount_upper = (smaller * (larger / total.next_down()).next_up()).next_up();
            (c * amount_upper).next_down()
        } else {
            0.0
        };
        let mut minimum = f64::INFINITY;
        for temperature in source_temperatures(&self.owner.heat_capacity, domain) {
            let mut capacity = 0.0;
            for (n, source) in component_moles.iter().zip(&self.owner.heat_capacity) {
                if *n > 0.0 {
                    let cp = lower_capacity_at(source, temperature);
                    let contribution = (n * cp).next_down().max(0.0);
                    capacity = (capacity + contribution).next_down();
                }
            }
            minimum = minimum.min((capacity + excess_lower).next_down());
        }
        if !positive(minimum) {
            return Err(BinaryCaloricErrorV1::Arithmetic(
                "positive heat-capacity lower bound is unrepresentable",
            ));
        }
        Ok(minimum)
    }
}

// At cropped domain endpoints or the other species' knots, bound the original
// affine interpolation itself. Rounding only n*Cp downward would not cover
// earlier upward roundoff in the temperature fraction/interpolation. Build
// upward from the smaller positive Cp using a downward temperature weight.
fn lower_capacity_at(source: &PureLiquidHeatCapacityV1, temperature: Temperature) -> f64 {
    let t = temperature.value();
    let knots = source.curve.knots_si();
    let index = knots
        .partition_point(|&(knot, _)| knot <= t)
        .saturating_sub(1)
        .min(knots.len() - 2);
    let [(a, ca), (b, cb)] = [knots[index], knots[index + 1]];
    if t == a {
        return ca;
    }
    if t == b {
        return cb;
    }
    let (base, increase, distance) = if ca <= cb {
        (
            ca,
            (cb - ca).next_down().max(0.0),
            (t - a).next_down().max(0.0),
        )
    } else {
        (
            cb,
            (ca - cb).next_down().max(0.0),
            (b - t).next_down().max(0.0),
        )
    };
    let weight = (distance / (b - a).next_up()).next_down().max(0.0);
    let increment = (increase * weight).next_down().max(0.0);
    (base + increment).next_down().max(0.0)
}

// No sampled caloric chart or extra allocation: duplicate shared knots are
// harmless, and admission/bounds are prepared once, outside thermal trials.
fn source_temperatures(
    sources: &[PureLiquidHeatCapacityV1; 2],
    domain: [Temperature; 2],
) -> impl Iterator<Item = Temperature> + '_ {
    domain.into_iter().chain(
        sources
            .iter()
            .flat_map(|source| source.curve.knots_si().iter())
            .filter(move |&&(t, _)| t > domain[0].value() && t < domain[1].value())
            .map(|&(t, _)| Temperature::new(t)),
    )
}
