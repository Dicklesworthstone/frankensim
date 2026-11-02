//! Explicit dry-mass weighting of specific grinding energies at one product size.

use super::{
    BondGrindingCard, ComminutionError, SpecificGrindingEnergy, positive, require,
    throughput_from_net_power,
};
use fs_qty::{Dimensionless, Length, Mass, MassFlowRate, Power};

/// One ore/composite feed with its own measured or authored work-index context.
#[derive(Clone, Debug, PartialEq)]
pub struct GrindingFeed {
    /// Original source-specific grinding card, including its applicability window.
    pub card: BondGrindingCard,
    /// Dry-solids mass of this feed in the blend, in kilograms.
    pub dry_mass: Mass,
    /// This feed's 80%-passing size before the common grinding stage, in metres.
    pub feed_p80: Length,
}

/// A caller-declared additive-specific-energy approximation for an ore blend.
///
/// For one circuit and common target P80, evaluate every source card at its
/// own F80, then compute E_blend = sum(m_i E_i) / sum(m_i). Apply power only
/// AFTER this sum. Averaging work indices, feed sizes or throughputs does not
/// in general produce this result. The original feeds remain inspectable.
///
/// Constructing this type explicitly chooses energy additivity. It does not
/// establish that mixed ores grind independently in a real mill. Interaction,
/// classification, circulating load, recovery and calibration uncertainty need
/// separately supported models. Mineral names never supply a work index.
#[derive(Clone, Debug, PartialEq)]
pub struct AdditiveGrindingBlend {
    feeds: Vec<GrindingFeed>,
    scaled_masses: Vec<f64>,
    total_scaled_mass: f64,
    product_bounds: [Length; 2],
}

impl AdditiveGrindingBlend {
    /// Admit nonempty positive dry masses, valid feed sizes and a common product
    /// domain. The caller asserts the same circuit and additive-energy regime.
    /// No constituent is dropped or given an inferred source/test context.
    pub fn new(feeds: Vec<GrindingFeed>) -> Result<Self, ComminutionError> {
        require(
            !feeds.is_empty(),
            "grinding blend requires at least one feed",
        )?;
        let mut largest_mass = 0.0_f64;
        let mut minimum_product = 0.0_f64;
        let mut maximum_product = f64::INFINITY;
        for feed in &feeds {
            largest_mass = largest_mass.max(positive(
                feed.dry_mass.value(),
                "each blend feed requires finite positive dry mass",
            )?);
            let bounds = feed.card.domain().product_bounds();
            // This checks both the source F80 window and whether this feed can
            // reach any supported product size with actual reduction.
            feed.card.domain().admit(feed.feed_p80, bounds[0])?;
            minimum_product = minimum_product.max(bounds[0].value());
            maximum_product = maximum_product
                .min(bounds[1].value())
                .min(feed.feed_p80.value().next_down());
        }
        require(
            minimum_product <= maximum_product,
            "blend source cards have no common product-size support",
        )?;
        let scaled_masses = feeds
            .iter()
            .map(|feed| {
                positive(
                    feed.dry_mass.value() / largest_mass,
                    "blend mass ratio is unrepresentable",
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
        let total_scaled_mass = positive(
            compensated_sum(scaled_masses.iter().copied()),
            "blend total dry-mass ratio is unrepresentable",
        )?;
        Ok(Self {
            feeds,
            scaled_masses,
            total_scaled_mass,
            product_bounds: [Length::new(minimum_product), Length::new(maximum_product)],
        })
    }

    /// Original component masses, F80s, test methods and source provenance.
    #[must_use]
    pub fn feeds(&self) -> &[GrindingFeed] {
        &self.feeds
    }

    /// Inclusive common P80 window, also strictly smaller than every feed F80.
    #[must_use]
    pub const fn product_bounds(&self) -> [Length; 2] {
        self.product_bounds
    }

    /// Dry-mass-weighted energy in J/kg at a shared target P80.
    /// Every original card admits the requested point independently.
    pub fn specific_energy(
        &self,
        product: Length,
    ) -> Result<SpecificGrindingEnergy, ComminutionError> {
        let energies = self
            .feeds
            .iter()
            .map(|feed| feed.card.specific_energy(feed.feed_p80, product))
            .collect::<Result<Vec<_>, _>>()?;
        let largest_energy = energies.iter().map(|e| e.value()).fold(0.0_f64, f64::max);
        // Scale BOTH extensive factors before summing. Large masses or energies
        // must not overflow a weighted mean which is itself finite.
        let contributions = energies
            .iter()
            .zip(&self.scaled_masses)
            .map(|(energy, mass)| {
                positive(
                    (energy.value() / largest_energy) * mass,
                    "blend energy contribution is unrepresentable",
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
        let mean =
            (compensated_sum(contributions.into_iter()) / self.total_scaled_mass) * largest_energy;
        Ok(SpecificGrindingEnergy::new(positive(
            mean,
            "blend specific energy is unrepresentable",
        )?))
    }

    /// Blend dry-solids throughput in kg/s from energy and available pinion power.
    pub fn dry_throughput(
        &self,
        product: Length,
        net_pinion_power: Power,
        utilization: Dimensionless,
    ) -> Result<MassFlowRate, ComminutionError> {
        throughput_from_net_power(
            self.specific_energy(product)?,
            net_pinion_power,
            utilization,
        )
    }
}

fn compensated_sum(values: impl Iterator<Item = f64>) -> f64 {
    let mut sum = 0.0_f64;
    let mut correction = 0.0_f64;
    for value in values {
        let next = sum + value;
        correction += if sum >= value {
            (sum - next) + value
        } else {
            (value - next) + sum
        };
        sum = next;
    }
    sum + correction
}

#[cfg(test)]
#[path = "blend/tests.rs"]
mod tests;
