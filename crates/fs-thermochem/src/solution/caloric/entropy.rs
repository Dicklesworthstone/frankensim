//! Relative entropy on the same pure-Cp and liquid phase states as enthalpy.

use super::{
    BinaryCaloricErrorV1, BinaryCaloricStateV1, BinaryLiquidCaloricModelV1,
    MolarHeatCapacityCurveV1, PureLiquidHeatCapacityV1, Sum, checked, interpolate, positive,
};
use fs_math::det::ln;
use fs_qty::Temperature;
use std::sync::Arc;

impl MolarHeatCapacityCurveV1 {
    /// Analytical integral of the original piecewise-linear Cp(T)/T [J/(mol K)].
    /// Only crossed source segments are visited. Local positive weights retain
    /// small temperature changes without subtracting reference entropies or
    /// nearly equal logarithms. Cooling reverses the identical ascending sum.
    pub fn entropy_change_j_per_mol_k(
        &self,
        from: Temperature,
        to: Temperature,
    ) -> Result<f64, BinaryCaloricErrorV1> {
        self.admit(from.value())?;
        self.admit(to.value())?;
        let (lo, hi, sign) = if from.value() <= to.value() {
            (from.value(), to.value(), 1.0)
        } else {
            (to.value(), from.value(), -1.0)
        };
        if lo == hi {
            return Ok(0.0);
        }
        let start = self
            .knots
            .partition_point(|&(t, _)| t <= lo)
            .saturating_sub(1);
        let mut sum = Sum::default();
        for pair in self.knots[start..].windows(2) {
            let a = lo.max(pair[0].0);
            let b = hi.min(pair[1].0);
            if b > a {
                let cp_a = interpolate(pair, a);
                let cp_b = interpolate(pair, b);
                let [cold, hot] = entropy_weights(a, b);
                // Scale positive capacities before the weighted sum: splitting
                // a subnormal Cp into two endpoint products can erase both.
                let scale = cp_a.max(cp_b);
                let increment = checked(
                    scale * ((cp_a / scale) * cold + (cp_b / scale) * hot),
                    "sensible entropy increment",
                )?;
                if !positive(increment) {
                    return Err(BinaryCaloricErrorV1::Arithmetic(
                        "positive sensible entropy increment underflows",
                    ));
                }
                sum.add(increment)?;
            }
            if pair[1].0 >= hi {
                break;
            }
        }
        checked(sign * sum.value()?, "sensible entropy integral")
    }
}

impl PureLiquidHeatCapacityV1 {
    /// Analytical integral of the original pure-liquid Cp(T)/T [J/(mol K)].
    pub fn entropy_change_j_per_mol_k(
        &self,
        from: Temperature,
        to: Temperature,
    ) -> Result<f64, BinaryCaloricErrorV1> {
        self.curve.entropy_change_j_per_mol_k(from, to)
    }
}

impl BinaryLiquidCaloricModelV1 {
    /// Entropy change [J/K] of the SAME conserved component inventory under
    /// this original caloric/source binding. Pure-component entropy datums
    /// cancel. Sensible Cp/T increments and the actual phase mixing-entropy
    /// change are included exactly once; a heat input alone cannot determine
    /// the entropy exchanged with an unspecified external reservoir.
    ///
    /// For equilibrated endpoints at fixed pressure/composition this is the
    /// reversible integral dH/T, including changing coexistence compositions.
    /// An explicitly homogeneous initial state retains its declared mixing
    /// entropy when evaluating the total change caused by adiabatic demixing.
    /// Phase/root precision remains that of the supplied accepted states.
    pub fn entropy_change_j_per_k(
        &self,
        initial: &BinaryCaloricStateV1,
        final_state: &BinaryCaloricStateV1,
    ) -> Result<f64, BinaryCaloricErrorV1> {
        let owns = |state: &BinaryCaloricStateV1| {
            Arc::ptr_eq(&self.owner, &state.owner) || self.owner == state.owner
        };
        if !owns(initial)
            || !owns(final_state)
            || initial.component_moles != final_state.component_moles
        {
            return Err(BinaryCaloricErrorV1::Invalid(
                "entropy increment requires the same caloric owner and component inventory",
            ));
        }
        self.entropy_increment(&[initial], final_state)
    }

    // The heat owner has already admitted every original parcel and closed the
    // component inventory. Evaluate only the accepted endpoint, never thermal
    // root trials. This also admits mixing parcels with different compositions
    // and temperatures without introducing absolute species entropy datums.
    pub(super) fn entropy_increment(
        &self,
        inputs: &[&BinaryCaloricStateV1],
        final_state: &BinaryCaloricStateV1,
    ) -> Result<f64, BinaryCaloricErrorV1> {
        let mut entropy = Sum::default();
        for input in inputs {
            for (n, source) in input.component_moles.iter().zip(&self.owner.heat_capacity) {
                if *n != 0.0 {
                    entropy.add(scaled_entropy(
                        *n,
                        source.entropy_change_j_per_mol_k(
                            input.temperature,
                            final_state.temperature,
                        )?,
                        "parcel sensible entropy increment",
                    )?)?;
                }
            }
            // Symmetric coexistence phases have identical mixing entropies.
            // Use the same actual accepted phase as the enthalpy owner, not
            // the unequilibrated feed composition or a transport Cp. This also
            // avoids lever-rule roundoff between identical phase entropies.
            // Weight the final-minus-original molar entropy by each incoming
            // parcel's conserved amount, rather than subtracting large totals.
            entropy.add(scaled_entropy(
                input.total_moles,
                final_state.mixing.mixing_entropy_j_per_mol_k()
                    - input.mixing.mixing_entropy_j_per_mol_k(),
                "parcel mixing entropy increment",
            )?)?;
        }
        entropy.value()
    }
}

fn scaled_entropy(
    moles: f64,
    molar_change: f64,
    quantity: &'static str,
) -> Result<f64, BinaryCaloricErrorV1> {
    let change = checked(moles * molar_change, quantity)?;
    if change == 0.0 && molar_change != 0.0 {
        return Err(BinaryCaloricErrorV1::Arithmetic(
            "nonzero entropy contribution underflows",
        ));
    }
    Ok(change)
}

// Integral of an affine capacity is Cp(a)*cold + Cp(b)*hot, where
// hot = 1-ln(1+r)/r and cold = ln(1+r)-hot, r=(b-a)/a. Both weights
// are positive. Evaluate their series directly near zero; fixed-length
// arithmetic plus fs-math's logarithm retains the original deterministic lane.
fn entropy_weights(low: f64, high: f64) -> [f64; 2] {
    let width = high - low;
    let r = width / low;
    if r < 0.125 {
        let mut cold = 0.0;
        let mut hot = 0.0;
        // At r<=1/8 the omitted terms are below double-precision roundoff.
        for n in (1..=18).rev() {
            let sign = if n % 2 == 0 { -1.0 } else { 1.0 };
            let n = f64::from(n);
            cold = sign / (n * (n + 1.0)) + r * cold;
            hot = sign / (n + 1.0) + r * hot;
        }
        [r * cold, r * hot]
    } else {
        let logarithm = if r.is_finite() {
            // log1p correction retains the rounding lost in forming 1+r.
            let shifted = 1.0 + r;
            ln(shifted) - ((shifted - 1.0) - r) / shifted
        } else {
            ln(high) - ln(low)
        };
        let hot = 1.0 - (low / width) * logarithm;
        [logarithm - hot, hot]
    }
}
