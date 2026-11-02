//! Inversion of the same uncorrected, source-specific Bond card.
//! No efficiency correction, recovery prediction or calibration is introduced.

use super::{
    BondGrindingCard, ComminutionError, Dimensionless, Length, MassFlowRate, Power,
    SpecificGrindingEnergy, det, positive, require,
};

impl BondGrindingCard {
    /// Solve for P80 from positive specific energy [J/kg of dry solids].
    ///
    /// The original work-index source and size support remain authoritative.
    /// Endpoint energies use this card's forward evaluator; their exact binary64
    /// values round-trip to the declared endpoint. Other requests are inverted
    /// algebraically and must pass normal size admission; nothing is clamped.
    pub fn product_size_for_specific_energy(
        &self,
        feed: Length,
        energy: SpecificGrindingEnergy,
    ) -> Result<Length, ComminutionError> {
        let energy = positive(
            energy.value(),
            "inverse Bond energy must be finite and positive",
        )?;
        let [minimum, maximum] = self.domain.product_bounds();
        self.domain.admit(feed, minimum)?;
        let maximum_energy = self.specific_energy(feed, minimum)?.value();
        require(
            energy <= maximum_energy,
            "energy implies P80 below the card's minimum",
        )?;
        if energy == maximum_energy {
            return Ok(minimum);
        }
        if maximum < feed {
            let minimum_energy = self.specific_energy(feed, maximum)?.value();
            require(
                energy >= minimum_energy,
                "energy implies P80 above the card's maximum",
            )?;
            if energy == minimum_energy {
                return Ok(maximum);
            }
        }
        let normalized = positive(
            energy / self.index.energy.value(),
            "inverse Bond normalized energy is unrepresentable",
        )?;
        let reference = 100.0e-6;
        let feed_term = det::sqrt(reference / feed.value());
        let product = if normalized <= feed_term {
            // P = F/(1+r)^2. Rationalizing F-P preserves tiny reductions
            // when rounding 1+r alone would erase the supplied energy.
            let r = normalized / feed_term;
            let relative_reduction = r * (2.0 + r) / ((1.0 + r) * (1.0 + r));
            feed.value().mul_add(-relative_reduction, feed.value())
        } else {
            // Split the two divisions so a very small squared denominator
            // cannot underflow when the final product size is representable.
            let denominator = feed_term + normalized;
            (reference / denominator) / denominator
        };
        let product = Length::new(positive(
            product,
            "inverse Bond product size is unrepresentable",
        )?);
        self.domain.admit(feed, product)?;
        Ok(product)
    }

    /// Solve for P80 at a positive dry-solids flow [kg/s] and pinion power [W].
    /// Utilization must be in (0,1]; downtime cannot supply a positive flow.
    /// Energy is `power * utilization / dry_mass_flow`, without an implied
    /// motor efficiency, circulating load or slurry-water correction.
    pub fn product_size_for_throughput(
        &self,
        feed: Length,
        net_pinion_power: Power,
        utilization: Dimensionless,
        dry_mass_flow: MassFlowRate,
    ) -> Result<Length, ComminutionError> {
        let power = positive(
            net_pinion_power.value(),
            "net pinion power must be finite and positive",
        )?;
        let flow = positive(
            dry_mass_flow.value(),
            "inverse Bond dry mass flow must be finite and positive",
        )?;
        let u = utilization.value();
        require(
            u.is_finite() && u > 0.0 && u <= 1.0,
            "inverse Bond utilization must be finite and in (0,1]",
        )?;
        let energy = positive_product_ratio(power, u, flow)?;
        self.product_size_for_specific_energy(feed, SpecificGrindingEnergy::new(energy))
    }
}

// Validated positive operands. Bounded binary scaling avoids destroying a
// representable energy through intermediate power underflow or flow overflow.
fn positive_product_ratio(a: f64, b: f64, divisor: f64) -> Result<f64, ComminutionError> {
    let (am, ae) = split(a);
    let (bm, be) = split(b);
    let (dm, de) = split(divisor);
    let mut mantissa = (am * bm) / dm;
    let mut exponent = ae + be - de;
    if mantissa >= 2.0 {
        mantissa *= 0.5;
        exponent += 1;
    } else if mantissa < 1.0 {
        mantissa *= 2.0;
        exponent -= 1;
    }
    let value = if (-1022..=1023).contains(&exponent) {
        mantissa * f64::from_bits(((exponent + 1023) as u64) << 52)
    } else if (-1075..-1022).contains(&exponent) {
        (mantissa * f64::MIN_POSITIVE) * f64::from_bits(((exponent + 2045) as u64) << 52)
    } else {
        return Err(ComminutionError::Invalid(
            "inverse Bond energy is unrepresentable",
        ));
    };
    positive(value, "inverse Bond energy is unrepresentable")
}

fn split(mut value: f64) -> (f64, i32) {
    let mut correction = 0;
    if value < f64::MIN_POSITIVE {
        value *= 4_503_599_627_370_496.0;
        correction = -52;
    }
    let bits = value.to_bits();
    let exponent = ((bits >> 52) & 0x7ff) as i32 - 1023 + correction;
    let mantissa = f64::from_bits((bits & 0x000f_ffff_ffff_ffff) | 0x3ff0_0000_0000_0000);
    (mantissa, exponent)
}

#[cfg(test)]
#[path = "inverse/tests.rs"]
mod tests;
