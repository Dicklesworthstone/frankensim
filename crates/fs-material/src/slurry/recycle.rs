//! Steady rougher/cleaner separation with an explicit cleaner-tail return.
//!
//! Fresh feed and returned cleaner tail enter the rougher. Rougher concentrate
//! feeds the cleaner; cleaner concentrate is final product. One NONSELECTIVE
//! fraction of cleaner tail returns, and the remainder is an explicit purge.
//! Rougher tails and purge form the external tailings stream. Partitions are
//! prescribed component MASS fractions, with separate carrier-liquid fractions.
//!
//! For a component with rougher fraction r, cleaner fraction c and return s,
//! its positive exit fraction is e=(1-r)+r*(c+(1-c)*(1-s)). Mixed feed is F/e.
//! This direct steady solution needs no convergence iteration. Existing stream
//! operations independently replay every separation, split and mixing balance.
//! No unexplained balancing stream is introduced.
//!
//! The original density evidence and both separators remain inspectable.
//! Constant prescribed partitions must apply to the resulting recycled blend;
//! this arithmetic does not establish that calibration or predict capacity,
//! breakage, liberation, reagent response or residence-dependent recovery.
//! This is a recovery circuit, not a grinding/population-balance model.

use super::separation::{SelectiveSeparator, SeparationError, SeparationResult};
use super::{SlurryComponent, SlurryError, SlurryStream};
use fs_matdb::Provenance;
use fs_qty::{Dimensionless, MassFlowRate};

const CLOSURE_TOLERANCE: f64 = 2e-12;

/// A refused prescription, nonunique closed loop or numerical material balance.
#[derive(Clone, Debug, PartialEq)]
pub enum RecycleError {
    /// Original source, component or stream admission failed.
    Slurry(SlurryError),
    /// An original separator refused its feed or prescribed partition coverage.
    Separation(SeparationError),
    /// Invalid circuit declaration or unrepresentable arithmetic.
    Invalid(&'static str),
    /// This identity has no external exit; no unique finite steady state exists.
    NoEscape {
        /// Exact solid identity or liquid-carrier identity with no exit.
        component: String,
    },
    /// An independently replayed balance exceeds the declared numeric tolerance.
    Closure {
        /// Identity whose independently replayed flow does not close.
        component: String,
        /// Either the external circuit boundary or the internal recycle closure.
        boundary: &'static str,
        /// Signed residual in kg/s; it is never used as an adjustment stream.
        residual: MassFlowRate,
    },
}
impl core::fmt::Display for RecycleError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rougher-cleaner recycle: {self:?}")
    }
}
impl std::error::Error for RecycleError {}
impl From<SlurryError> for RecycleError {
    fn from(value: SlurryError) -> Self {
        Self::Slurry(value)
    }
}
impl From<SeparationError> for RecycleError {
    fn from(value: SeparationError) -> Self {
        Self::Separation(value)
    }
}
type Result<T> = core::result::Result<T, RecycleError>;

/// Original prescribed stage partitions and return/purge assumption.
#[derive(Clone, Debug, PartialEq)]
pub struct RougherCleanerCircuit {
    rougher: SelectiveSeparator,
    cleaner: SelectiveSeparator,
    cleaner_tail_recycle: Dimensionless,
    provenance: Provenance,
}

/// Flow accounting for one exact solid identity or the separate liquid carrier.
#[derive(Clone, Debug, PartialEq)]
pub struct CircuitComponentBalance {
    /// Exact dry constituent or liquid-carrier identity.
    pub component: String,
    /// Fresh entering mass rate, in kg/s.
    pub fresh_mass_flow: MassFlowRate,
    /// Total mass rate entering the rougher, in kg/s.
    pub mixed_mass_flow: MassFlowRate,
    /// Cleaner-tail mass returned to the rougher, in kg/s.
    pub recycle_mass_flow: MassFlowRate,
    /// Final cleaner-concentrate mass rate, in kg/s.
    pub product_mass_flow: MassFlowRate,
    /// Mass leaving directly through rougher tails, in kg/s.
    pub rougher_tailings_mass_flow: MassFlowRate,
    /// Cleaner-tail mass leaving through the explicit purge, in kg/s.
    pub purge_mass_flow: MassFlowRate,
    /// Product + rougher tails + purge - fresh, in kg/s.
    pub external_conservation_residual: MassFlowRate,
    /// Mixed feed - independently mixed (fresh + recycle), in kg/s.
    pub recycle_closure_residual: MassFlowRate,
    /// Recycle/fresh mass rate; may exceed one, and is None for zero fresh feed.
    pub circulating_load: Option<Dimensionless>,
    /// Final product/fresh mass rate; None for zero fresh feed of this identity.
    pub overall_recovery: Option<Dimensionless>,
}

/// Physical streams and original source evidence for the complete steady circuit.
#[derive(Clone, Debug, PartialEq)]
pub struct RecycleResult {
    /// Original stage prescriptions, return fraction and circuit provenance.
    pub circuit: RougherCleanerCircuit,
    /// Original fresh feed, including its density sources and query conditions.
    pub fresh_feed: SlurryStream,
    /// Directly solved rougher feed, retaining each original constituent source.
    pub mixed_feed: SlurryStream,
    /// First stage; its tailings leave the recycle loop.
    pub rougher: SeparationResult,
    /// Second stage; its concentrate is the final circuit product.
    pub cleaner: SeparationResult,
    /// Replayed cleaner-tail stream returned to the rougher.
    pub recycle: SlurryStream,
    /// Complementary cleaner-tail stream leaving the circuit.
    pub purge: SlurryStream,
    /// Rougher tailings plus cleaner-tail purge, using the existing source mixer.
    pub combined_tailings: SlurryStream,
    /// Dry identities retain their existing canonical order and source density.
    pub components: Vec<CircuitComponentBalance>,
    /// Liquid-carrier balance; excluded from dry component grades and yields.
    pub carrier: CircuitComponentBalance,
}

impl RougherCleanerCircuit {
    /// Declare two existing separators and the fraction of every cleaner-tail
    /// component returned. The complement is purged; no hidden water recovery or
    /// fresh makeup stream is inferred. Supply makeup in the explicit fresh feed.
    pub fn new(
        rougher: SelectiveSeparator,
        cleaner: SelectiveSeparator,
        cleaner_tail_recycle: Dimensionless,
        provenance: Provenance,
    ) -> Result<Self> {
        let fraction = cleaner_tail_recycle.value();
        if !fraction.is_finite()
            || !(0.0..=1.0).contains(&fraction)
            || provenance.source.trim().is_empty()
            || provenance.license.trim().is_empty()
        {
            return Err(RecycleError::Invalid(
                "return fraction must be in [0,1] and circuit source/license must be nonblank",
            ));
        }
        Ok(Self {
            rougher,
            cleaner,
            cleaner_tail_recycle: Dimensionless::new(fraction.abs()),
            provenance,
        })
    }

    /// Original prescribed rougher partitions and their provenance.
    pub fn rougher(&self) -> &SelectiveSeparator {
        &self.rougher
    }
    /// Original prescribed cleaner partitions and their provenance.
    pub fn cleaner(&self) -> &SelectiveSeparator {
        &self.cleaner
    }
    /// Nonselective fraction of cleaner tail returned, on each component's mass.
    pub const fn cleaner_tail_recycle(&self) -> Dimensionless {
        self.cleaner_tail_recycle
    }
    /// Original source of the circuit topology and return/purge assumption.
    pub fn provenance(&self) -> &Provenance {
        &self.provenance
    }

    /// Solve directly and replay the original separator/split/mix operations.
    ///
    /// Both stages must cover the fresh solid identity set exactly, including
    /// zero-flow identities. A fully trapped identity refuses even at zero fresh
    /// flow: its circulating inventory would otherwise remain undetermined.
    /// No input is mutated. Work is bounded by the existing 128-solid limit.
    ///
    /// Both relative closure residuals must be at most 2e-12. External closure
    /// is scaled by FRESH flow, never by the potentially much larger recycle.
    /// An ill-conditioned near-trapped circuit can therefore refuse even when
    /// each separate stage is representable; no rounding deficit is patched.
    pub fn solve(&self, fresh: &SlurryStream) -> Result<RecycleResult> {
        coverage(&self.rougher, fresh)?;
        coverage(&self.cleaner, fresh)?;
        let s = self.cleaner_tail_recycle.value();
        let solids = fresh
            .solids()
            .iter()
            .zip(self.rougher.partitions())
            .zip(self.cleaner.partitions())
            .map(|((component, rougher), cleaner)| {
                mixed_component(
                    component,
                    rougher.to_concentrate().value(),
                    cleaner.to_concentrate().value(),
                    s,
                )
            })
            .collect::<Result<Vec<_>>>()?;
        let liquid = mixed_component(
            fresh.liquid(),
            self.rougher.carrier_to_concentrate().value(),
            self.cleaner.carrier_to_concentrate().value(),
            s,
        )?;
        let mixed_feed = SlurryStream::new_additive_volume(solids, liquid)?;
        let rougher = self.rougher.separate(&mixed_feed)?;
        let cleaner = self.cleaner.separate(&rougher.concentrate)?;
        let (recycle, purge) = cleaner.tailings.split(self.cleaner_tail_recycle)?;
        let combined_tailings = SlurryStream::mix(&[&rougher.tailings, &purge])?;
        let reconstructed_feed = SlurryStream::mix(&[fresh, &recycle])?;
        let streams = [
            fresh,
            &mixed_feed,
            &recycle,
            &cleaner.concentrate,
            &rougher.tailings,
            &purge,
            &combined_tailings,
            &reconstructed_feed,
        ];
        let components = (0..fresh.solids().len())
            .map(|i| balance(streams.map(|stream| &stream.solids()[i])))
            .collect::<Result<Vec<_>>>()?;
        let carrier = balance(streams.map(SlurryStream::liquid))?;
        Ok(RecycleResult {
            circuit: self.clone(),
            fresh_feed: fresh.clone(),
            mixed_feed,
            rougher,
            cleaner,
            recycle,
            purge,
            combined_tailings,
            components,
            carrier,
        })
    }
}

fn coverage(separator: &SelectiveSeparator, fresh: &SlurryStream) -> Result<()> {
    let missing: Vec<_> = fresh
        .solids()
        .iter()
        .filter(|solid| {
            separator
                .partitions()
                .binary_search_by(|p| p.component().cmp(solid.id()))
                .is_err()
        })
        .map(|solid| solid.id().to_owned())
        .collect();
    let unknown: Vec<_> = separator
        .partitions()
        .iter()
        .filter(|partition| {
            fresh
                .solids()
                .binary_search_by(|s| s.id().cmp(partition.component()))
                .is_err()
        })
        .map(|partition| partition.component().to_owned())
        .collect();
    if !missing.is_empty() || !unknown.is_empty() {
        Err(SeparationError::PartitionCoverage { missing, unknown }.into())
    } else {
        Ok(())
    }
}

fn mixed_component(component: &SlurryComponent, r: f64, c: f64, s: f64) -> Result<SlurryComponent> {
    // Sum positive exit paths instead of subtracting a rounded loop gain from
    // one. In particular r=s=1 and any positive c must retain the product exit.
    let escape = (1.0 - r) + r * (c + (1.0 - c) * (1.0 - s));
    if escape == 0.0 {
        return Err(RecycleError::NoEscape {
            component: component.id().to_owned(),
        });
    }
    if !escape.is_finite() || escape <= 0.0 || escape > 1.0 {
        return Err(RecycleError::Invalid(
            "component escape fraction is unrepresentable",
        ));
    }
    let mixed = component.mass_flow().value() / escape;
    if !mixed.is_finite() || (component.mass_flow().value() > 0.0 && mixed == 0.0) {
        return Err(RecycleError::Invalid(
            "circulating component flow is unrepresentable",
        ));
    }
    Ok(SlurryComponent::new(
        component.id(),
        MassFlowRate::new(mixed),
        component.density_source().clone(),
    )?)
}

fn balance(components: [&SlurryComponent; 8]) -> Result<CircuitComponentBalance> {
    let id = components[0].id();
    let [
        fresh,
        mixed,
        recycle,
        product,
        rougher_tail,
        purge,
        external_tail,
        rebuilt,
    ] = components.map(|c| c.mass_flow().value());
    let external_residual = (product + external_tail) - fresh;
    let recycle_residual = mixed - rebuilt;
    closure(id, "external", external_residual, fresh)?;
    closure(id, "recycle", recycle_residual, mixed)?;
    let circulating_load = ratio(recycle, fresh)?;
    let overall_recovery = ratio(product, fresh)?;
    if overall_recovery
        .is_some_and(|r| r.value() > 1.0 || (external_tail > 0.0 && r.value() == 1.0))
    {
        return Err(RecycleError::Invalid(
            "positive overall recovery or tail fraction is unrepresentable",
        ));
    }
    Ok(CircuitComponentBalance {
        component: id.to_owned(),
        fresh_mass_flow: MassFlowRate::new(fresh),
        mixed_mass_flow: MassFlowRate::new(mixed),
        recycle_mass_flow: MassFlowRate::new(recycle),
        product_mass_flow: MassFlowRate::new(product),
        rougher_tailings_mass_flow: MassFlowRate::new(rougher_tail),
        purge_mass_flow: MassFlowRate::new(purge),
        external_conservation_residual: MassFlowRate::new(external_residual),
        recycle_closure_residual: MassFlowRate::new(recycle_residual),
        circulating_load,
        overall_recovery,
    })
}

fn closure(id: &str, boundary: &'static str, residual: f64, scale: f64) -> Result<()> {
    if !residual.is_finite()
        || (scale == 0.0 && residual != 0.0)
        || (scale > 0.0 && (residual / scale).abs() > CLOSURE_TOLERANCE)
    {
        Err(RecycleError::Closure {
            component: id.to_owned(),
            boundary,
            residual: MassFlowRate::new(residual),
        })
    } else {
        Ok(())
    }
}

fn ratio(numerator: f64, denominator: f64) -> Result<Option<Dimensionless>> {
    if denominator == 0.0 {
        return Ok(None);
    }
    let value = numerator / denominator;
    if !value.is_finite() || (numerator > 0.0 && value == 0.0) {
        Err(RecycleError::Invalid(
            "circuit flow ratio is unrepresentable",
        ))
    } else {
        Ok(Some(Dimensionless::new(value)))
    }
}

#[cfg(test)]
#[path = "recycle_tests.rs"]
mod tests;
