//! Persistent shrinking cores competing for one finite, well-mixed reagent.
//!
//! All original local laws are first order in the same bulk concentration C.
//! The shared exposure dH=C dt therefore separates particle geometry from the
//! bath balance: F_i(u_i)-F_i(u_i,start)=H for every unfinished population.
//! F_i is the integral of the original film/surface/product-layer resistance.
//! A scalar time integral then closes the common reagent stock. No population
//! receives priority, and independently computed recovery curves are not clipped.
//!
//! Original radii, mineral amounts, source cards and initial solution are fixed.
//! Repeated advances retain each accepted core and the absolute source clock.
//! Original mineral-total loading, summed over its size bins, is the source
//! coordinate; merely subdividing a bin does not change that coordinate.
//! This closure assumes additive local rates in one homogeneous bath, constant
//! volume/temperature/coefficients, and no product inhibition or reaction feed-
//! back. Source applicability is required; mixture calibration is not inferred.

use super::{
    ShrinkingCoreError, ShrinkingCoreReaction, conversion, nonnegative, positive, quotient,
};
use fs_matdb::QueryPoint;
use fs_qty::{Amount, AmountConcentration, Dimensionless, Length, Qty, Temperature, Time, Volume};
use fs_thermochem::{
    ReactionId, SpeciesId,
    stoichiometry::{
        ConservedReaction,
        network::{NetworkAdvance, ReactionNetwork, ReactionNetworkError},
    },
};
use std::{
    collections::{BTreeMap, BTreeSet},
    sync::Arc,
};

#[path = "population/exposure.rs"]
mod exposure;
#[path = "population/integration.rs"]
mod integration;
use exposure::{CoreExposure, CoreProgress};

/// Integral of reagent concentration over time, in mol s/m³.
pub type ReagentExposure = Qty<-3, 0, 1, 0, 0, 1>;
/// Bound on distinct mineral/size populations in one batch.
pub const MAX_POPULATIONS: usize = 128;

/// A refusal leaves the accepted bath, particle geometry and clock unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum PopulationBatchError {
    /// Original local-law admission or bounded numerical refusal.
    Core(ShrinkingCoreError),
    /// Original exact reaction-network or finite-stock refusal.
    Network(ReactionNetworkError),
    /// Invalid population identity, layout or physical association.
    Invalid(&'static str),
    /// Nonfinite or unrepresentable positive stock, progress or history.
    Arithmetic(&'static str),
}
impl core::fmt::Display for PopulationBatchError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "shrinking-core population batch: {self:?}")
    }
}
impl std::error::Error for PopulationBatchError {}
impl From<ShrinkingCoreError> for PopulationBatchError {
    fn from(e: ShrinkingCoreError) -> Self {
        Self::Core(e)
    }
}
impl From<ReactionNetworkError> for PopulationBatchError {
    fn from(e: ReactionNetworkError) -> Self {
        Self::Network(e)
    }
}
type Result<T> = core::result::Result<T, PopulationBatchError>;

/// One originally virgin mineral/size population with an explicit source mapping.
#[derive(Clone, Debug)]
pub struct LeachPopulation {
    id: String,
    model: ShrinkingCoreReaction,
    radius: Length,
    original_amount: Amount,
    conditions: QueryPoint,
}
impl LeachPopulation {
    /// Declare a population. Radius is reactive-grain radius, not an inferred
    /// half-size of an unqualified bulk ore particle. Amount is original solid.
    pub fn new(
        id: impl Into<String>,
        model: ShrinkingCoreReaction,
        radius: Length,
        original_amount: Amount,
        conditions: QueryPoint,
    ) -> Result<Self> {
        let id = id.into();
        if id.trim().is_empty() || id.len() > 256 {
            return Err(PopulationBatchError::Invalid(
                "population IDs need 1..=256 bytes",
            ));
        }
        positive(radius.value(), "population radius must be positive")?;
        nonnegative(
            original_amount.value(),
            "original solid amount must be nonnegative",
        )?;
        Ok(Self {
            id,
            model,
            radius,
            original_amount,
            conditions,
        })
    }
    /// Original caller identity, independent of species and reaction IDs.
    pub fn id(&self) -> &str {
        &self.id
    }
    /// Original local law, exact equation, coefficients and source evidence.
    pub fn model(&self) -> &ShrinkingCoreReaction {
        &self.model
    }
    /// Fixed outer reactive-particle radius.
    pub const fn radius(&self) -> Length {
        self.radius
    }
    /// Original solid amount; it is never reset after partial conversion.
    pub const fn original_amount(&self) -> Amount {
        self.original_amount
    }
    /// Original additional source coordinates.
    pub fn conditions(&self) -> &QueryPoint {
        &self.conditions
    }
}

/// Independent work limits and estimated time-inversion tolerance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PopulationNumerics {
    /// Positive absolute time residual allowance.
    pub absolute_time_tolerance: Time,
    /// Relative time residual allowance in (0, 0.01].
    pub relative_time_tolerance: f64,
    /// Total scalar and inner-core evaluation budget, at most ten million.
    pub max_work: usize,
    /// Scalar time-inversion iterations, in 1..=128.
    pub max_iterations: usize,
    /// Iterations for any one intrinsic core inversion, in 1..=128.
    pub max_core_iterations: usize,
}
impl Default for PopulationNumerics {
    fn default() -> Self {
        Self {
            absolute_time_tolerance: Time::new(1e-8),
            relative_time_tolerance: 1e-9,
            max_work: 1_000_000,
            max_iterations: 96,
            max_core_iterations: 96,
        }
    }
}
impl PopulationNumerics {
    fn tolerance(self, t: f64) -> Result<f64> {
        if !self.absolute_time_tolerance.value().is_finite()
            || self.absolute_time_tolerance.value() <= 0.0
            || !self.relative_time_tolerance.is_finite()
            || !(0.0 < self.relative_time_tolerance && self.relative_time_tolerance <= 0.01)
            || !(1..=10_000_000).contains(&self.max_work)
            || !(1..=128).contains(&self.max_iterations)
            || !(1..=128).contains(&self.max_core_iterations)
        {
            return Err(PopulationBatchError::Invalid(
                "invalid population numerical budget",
            ));
        }
        Ok(positive(
            self.absolute_time_tolerance.value() + self.relative_time_tolerance * t,
            "population time tolerance range loss",
        )?)
    }
}

/// Geometry and molar account of one accepted population.
#[derive(Clone, Debug, PartialEq)]
pub struct PopulationState {
    /// Caller identity retained through all advances.
    pub id: String,
    /// Original solid species identity.
    pub solid: SpeciesId,
    /// Original reaction identity; multiple size bins may share it.
    pub reaction: ReactionId,
    /// Original amount, in mol of this solid formula unit.
    pub original_amount: Amount,
    /// Total written-equation extent for this population since initialization.
    pub cumulative_extent: Amount,
    /// Molar consumption from that written-equation extent.
    pub consumed: Amount,
    /// Original amount minus stoichiometric consumption.
    pub remaining: Amount,
    /// Spherical core conversion; not metal recovery or payability.
    pub conversion: Dimensionless,
    /// Remaining core radius; outer radius is unchanged.
    pub core_radius: Length,
    /// Geometry-derived remainder minus molar-account remainder; no correction
    /// is applied. Binary64 residuals are not conservation certificates.
    pub geometry_residual: Amount,
}

#[derive(Clone, Debug)]
struct Layout {
    populations: Vec<LeachPopulation>,
    curves: Vec<CoreExposure>,
    network: ReactionNetwork,
    slots: Vec<usize>,
    loadings: Vec<f64>,
    initial_inventory: Vec<(SpeciesId, Amount)>,
    reagent: SpeciesId,
    initial_reagent: f64,
    volume: Volume,
    temperature: Temperature,
}

/// Accepted private particle history and its shared finite bath.
#[derive(Clone, Debug)]
pub struct ShrinkingCoreBatch {
    layout: Arc<Layout>,
    penetrations: Vec<f64>,
    states: Vec<PopulationState>,
    inventory: Vec<(SpeciesId, Amount)>,
    elapsed: Time,
    exposure: ReagentExposure,
    accumulated_time_error: Time,
}

/// One atomic accepted advance, retaining the immutable original model layout.
#[derive(Clone, Debug)]
pub struct PopulationBatchAdvance {
    layout: Arc<Layout>,
    /// Time at the start of this call, measured from original batch initialization.
    pub previous_elapsed: Time,
    /// Absolute accepted batch clock, including a post-completion hold.
    pub elapsed: Time,
    /// Reagent exposure added by this call.
    pub exposure_increment: ReagentExposure,
    /// Total exposure since initialization.
    pub cumulative_exposure: ReagentExposure,
    /// Original-network account from ORIGINAL stocks and total reaction extents.
    /// These consumed/formed amounts are cumulative, not incremental step flows.
    pub cumulative: NetworkAdvance,
    /// Canonical population identities and their accepted core/molar accounts.
    pub populations: Vec<PopulationState>,
    /// Actual source conditions at the previous absolute time.
    pub initial_points: Vec<QueryPoint>,
    /// Actual source conditions at the accepted absolute time.
    pub final_points: Vec<QueryPoint>,
    /// Time residual plus adaptive quadrature estimate; not an interval bound.
    pub estimated_time_error: Time,
    /// Sum of reported step error estimates, not a certified global error bound.
    pub accumulated_time_error: Time,
    /// All charged work, including nested core inversions.
    pub work: usize,
    /// Number of time-integrand evaluations within that work.
    pub quadrature_evaluations: usize,
}
impl PopulationBatchAdvance {
    /// Original complete source-bound population declarations.
    pub fn models(&self) -> &[LeachPopulation] {
        &self.layout.populations
    }
    /// Original exact network, usable with the existing source-enthalpy consumer.
    pub fn reaction_network(&self) -> &ReactionNetwork {
        &self.layout.network
    }
}

impl ShrinkingCoreBatch {
    /// Initialize 1..=128 explicitly assigned populations in one closed bath.
    /// Supply the complete network inventory, including zero products. Each
    /// solid stock must equal the canonical sum of its assigned original bins.
    /// Duplicate reaction IDs require the same original conservation certificate.
    /// Products cannot feed any modeled reactant in this irreversible closure.
    pub fn new(
        mut populations: Vec<LeachPopulation>,
        reagent: SpeciesId,
        inventory: &[(SpeciesId, Amount)],
        volume: Volume,
        temperature: Temperature,
    ) -> Result<Self> {
        if populations.is_empty() || populations.len() > MAX_POPULATIONS {
            return Err(PopulationBatchError::Invalid(
                "requires 1..=128 populations",
            ));
        }
        positive(volume.value(), "batch liquid volume must be positive")?;
        positive(
            temperature.value(),
            "batch absolute temperature must be positive",
        )?;
        populations.sort_by(|a, b| a.id.cmp(&b.id));
        if populations.windows(2).any(|p| p[0].id == p[1].id) {
            return Err(PopulationBatchError::Invalid(
                "duplicate population identity",
            ));
        }
        let mut reactions: BTreeMap<ReactionId, ConservedReaction> = BTreeMap::new();
        let mut solids: BTreeMap<SpeciesId, Sum> = BTreeMap::new();
        let mut consumed = BTreeSet::from([reagent.clone()]);
        for population in &populations {
            let model = &population.model;
            if model.reagent != reagent {
                return Err(PopulationBatchError::Invalid(
                    "every population must consume the same reagent",
                ));
            }
            consumed.insert(model.solid.clone());
            solids
                .entry(model.solid.clone())
                .or_default()
                .add(population.original_amount.value())?;
            if let Some(previous) = reactions.get(model.reaction.id()) {
                if previous.terms() != model.reaction.terms()
                    || previous.conservation() != model.reaction.conservation()
                {
                    return Err(PopulationBatchError::Invalid(
                        "a repeated reaction ID changed its original equation",
                    ));
                }
            } else {
                reactions.insert(model.reaction.id().clone(), model.reaction.clone());
            }
        }
        if reactions.values().any(|r| {
            r.terms()
                .iter()
                .any(|t| t.coefficient > 0 && consumed.contains(&t.species))
        }) {
            return Err(PopulationBatchError::Invalid(
                "reaction products cannot feed modeled reactants",
            ));
        }
        let network = ReactionNetwork::new(reactions.into_values().collect())?;
        let zero: Vec<_> = network
            .reactions()
            .iter()
            .map(|r| (r.id().clone(), Amount::new(0.0)))
            .collect();
        let initial = network.advance_inventory(inventory, &zero)?.final_inventory;
        let mut totals = BTreeMap::new();
        for (species, sum) in solids {
            let total = sum.value()?;
            if total != stock(&initial, &species) {
                return Err(PopulationBatchError::Invalid(
                    "solid stock must equal its original assigned population sum",
                ));
            }
            totals.insert(species, total);
        }
        let initial_reagent = stock(&initial, &reagent);
        let mut curves = Vec::with_capacity(populations.len());
        let mut slots = Vec::with_capacity(populations.len());
        let mut loadings = Vec::with_capacity(populations.len());
        for population in &populations {
            curves.push(CoreExposure::new(&population.model, population.radius)?);
            slots.push(
                network
                    .reactions()
                    .binary_search_by(|r| r.id().cmp(population.model.reaction.id()))
                    .expect("original reaction retained"),
            );
            loadings.push(quotient(totals[&population.model.solid], volume.value())?);
        }
        let states = populations
            .iter()
            .map(|p| PopulationState {
                id: p.id.clone(),
                solid: p.model.solid.clone(),
                reaction: p.model.reaction.id().clone(),
                original_amount: p.original_amount,
                cumulative_extent: Amount::new(0.0),
                consumed: Amount::new(0.0),
                remaining: p.original_amount,
                conversion: Dimensionless::new(0.0),
                core_radius: p.radius,
                geometry_residual: Amount::new(0.0),
            })
            .collect();
        let n = populations.len();
        let layout = Arc::new(Layout {
            populations,
            curves,
            network,
            slots,
            loadings,
            initial_inventory: initial.clone(),
            reagent,
            initial_reagent,
            volume,
            temperature,
        });
        layout.points(Time::new(0.0), initial_reagent)?;
        Ok(Self {
            layout,
            penetrations: vec![0.0; n],
            states,
            inventory: initial,
            elapsed: Time::new(0.0),
            exposure: ReagentExposure::new(0.0),
            accumulated_time_error: Time::new(0.0),
        })
    }
    /// Original sources and fixed population assignments, in identity order.
    pub fn models(&self) -> &[LeachPopulation] {
        &self.layout.populations
    }
    /// Current accepted complete species inventory.
    pub fn inventory(&self) -> &[(SpeciesId, Amount)] {
        &self.inventory
    }
    /// Current accepted per-population core and molar accounts.
    pub fn populations(&self) -> &[PopulationState] {
        &self.states
    }
    /// Original exact network; repeated size bins are aggregated by reaction ID.
    pub fn reaction_network(&self) -> &ReactionNetwork {
        &self.layout.network
    }
    /// Absolute batch age. It is never reset on continuation.
    pub const fn elapsed(&self) -> Time {
        self.elapsed
    }
    /// Total accepted reagent exposure.
    pub const fn exposure(&self) -> ReagentExposure {
        self.exposure
    }

    /// Continue the existing particle histories for a nonnegative duration.
    /// No state changes until all source, stock, work and cancellation checks pass.
    pub fn advance(
        &mut self,
        duration: Time,
        numerics: PopulationNumerics,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PopulationBatchAdvance> {
        nonnegative(duration.value(), "advance duration must be nonnegative")?;
        let tolerance = numerics.tolerance(duration.value())?;
        let end = self.elapsed.value() + duration.value();
        if !end.is_finite() || (duration.value() > 0.0 && end == self.elapsed.value()) {
            return Err(PopulationBatchError::Arithmetic(
                "absolute batch clock loses range",
            ));
        }
        let reagent = stock(&self.inventory, &self.layout.reagent);
        let initial_points = self.layout.points(self.elapsed, reagent)?;
        self.layout.points(Time::new(end), reagent)?;
        let mut work = Work {
            used: 0,
            maximum: numerics.max_work,
            integrands: 0,
            cancelled: &mut cancelled,
        };
        work.tick()?;
        let mut evaluator = Evaluator::new(self, numerics.max_core_iterations);
        let (dose, error) = if duration.value() == 0.0 || reagent == 0.0 {
            (0.0, 0.0)
        } else if self.states.iter().all(|s| s.remaining.value() == 0.0) {
            (
                positive(
                    quotient(reagent, self.layout.volume.value())? * duration.value(),
                    "post-completion exposure loses range",
                )?,
                0.0,
            )
        } else {
            integration::solve(
                &mut evaluator,
                duration.value(),
                tolerance,
                numerics,
                &mut work,
            )?
        };
        evaluator.evaluate(dose, &mut work)?;
        let requested: Vec<_> = self
            .layout
            .network
            .reactions()
            .iter()
            .zip(&evaluator.extents)
            .map(|(r, &x)| (r.id().clone(), Amount::new(x)))
            .collect();
        let cumulative = self
            .layout
            .network
            .advance_inventory(&self.layout.initial_inventory, &requested)?;
        let remaining_reagent = stock(&cumulative.final_inventory, &self.layout.reagent);
        if reagent > 0.0 && remaining_reagent <= 0.0 {
            return Err(PopulationBatchError::Arithmetic(
                "finite time erased positive reagent stock",
            ));
        }
        let final_points = self.layout.points(Time::new(end), remaining_reagent)?;
        let states = self.layout.states(&evaluator.progress)?;
        let mut mineral_accounts: BTreeMap<SpeciesId, Sum> = BTreeMap::new();
        for state in &states {
            mineral_accounts
                .entry(state.solid.clone())
                .or_default()
                .add(state.remaining.value())?;
        }
        for (species, sum) in mineral_accounts {
            let reported = stock(&cumulative.final_inventory, &species);
            let original = stock(&self.layout.initial_inventory, &species);
            if (sum.value()? - reported).abs() > 64.0 * f64::EPSILON * original {
                return Err(PopulationBatchError::Arithmetic(
                    "population and network solid accounts disagree",
                ));
            }
        }
        let total_dose = self.exposure.value() + dose;
        if !total_dose.is_finite() || (dose > 0.0 && total_dose == self.exposure.value()) {
            return Err(PopulationBatchError::Arithmetic(
                "cumulative exposure loses range",
            ));
        }
        let accumulated = self.accumulated_time_error.value() + error;
        nonnegative(accumulated, "accumulated time-error estimate loses range")?;
        work.tick()?;
        let penetrations = evaluator.progress.iter().map(|p| p.penetration).collect();
        let report = PopulationBatchAdvance {
            layout: self.layout.clone(),
            previous_elapsed: self.elapsed,
            elapsed: Time::new(end),
            exposure_increment: ReagentExposure::new(dose),
            cumulative_exposure: ReagentExposure::new(total_dose),
            cumulative,
            populations: states.clone(),
            initial_points,
            final_points,
            estimated_time_error: Time::new(error),
            accumulated_time_error: Time::new(accumulated),
            work: work.used,
            quadrature_evaluations: work.integrands,
        };
        self.penetrations = penetrations;
        self.states = states;
        self.inventory = report.cumulative.final_inventory.clone();
        self.elapsed = report.elapsed;
        self.exposure = report.cumulative_exposure;
        self.accumulated_time_error = report.accumulated_time_error;
        Ok(report)
    }
}

impl Layout {
    fn points(&self, time: Time, reagent: f64) -> Result<Vec<QueryPoint>> {
        let concentration = AmountConcentration::new(quotient(reagent, self.volume.value())?);
        self.populations
            .iter()
            .zip(&self.loadings)
            .map(|(p, &loading)| {
                Ok(p.model.point(
                    self.volume,
                    p.radius,
                    concentration,
                    AmountConcentration::new(loading),
                    time,
                    self.temperature,
                    &p.conditions,
                )?)
            })
            .collect()
    }
    fn states(&self, progress: &[CoreProgress]) -> Result<Vec<PopulationState>> {
        self.populations
            .iter()
            .zip(progress)
            .map(|(p, core)| {
                let x = conversion(core.penetration);
                let extent = extent(p, x)?;
                let consumed = (p.model.coefficients[0] as f64) * extent;
                let remaining =
                    (-(p.model.coefficients[0] as f64)).mul_add(extent, p.original_amount.value());
                if !remaining.is_finite()
                    || remaining < 0.0
                    || (p.original_amount.value() > 0.0 && !core.complete && remaining == 0.0)
                    || (core.complete && remaining != 0.0)
                {
                    return Err(PopulationBatchError::Arithmetic(
                        "population extent lost its positive core stock",
                    ));
                }
                let radius = p.radius.value() * (1.0 - core.penetration);
                if !core.complete && radius == 0.0 {
                    return Err(PopulationBatchError::Arithmetic("core radius underflows"));
                }
                let geometric = p.original_amount.value() * core.remaining_fraction;
                let residual = geometric - remaining;
                if residual.abs() > 64.0 * f64::EPSILON * p.original_amount.value() {
                    return Err(PopulationBatchError::Arithmetic(
                        "geometric and molar core accounts disagree",
                    ));
                }
                Ok(PopulationState {
                    id: p.id.clone(),
                    solid: p.model.solid.clone(),
                    reaction: p.model.reaction.id().clone(),
                    original_amount: p.original_amount,
                    cumulative_extent: Amount::new(extent),
                    consumed: Amount::new(consumed),
                    remaining: Amount::new(remaining),
                    conversion: Dimensionless::new(x),
                    core_radius: Length::new(radius),
                    geometry_residual: Amount::new(residual),
                })
            })
            .collect()
    }
}

struct Evaluator<'a> {
    batch: &'a ShrinkingCoreBatch,
    progress: Vec<CoreProgress>,
    extents: Vec<f64>,
    sums: Vec<Sum>,
    core_iterations: usize,
}
impl<'a> Evaluator<'a> {
    fn new(batch: &'a ShrinkingCoreBatch, core_iterations: usize) -> Self {
        Self {
            batch,
            progress: batch
                .penetrations
                .iter()
                .map(|&u| CoreProgress {
                    penetration: u,
                    conversion_increment: 0.0,
                    remaining_fraction: (1.0 - u).powi(3),
                    complete: u == 1.0,
                })
                .collect(),
            extents: vec![0.0; batch.layout.network.reactions().len()],
            sums: vec![Sum::default(); batch.layout.network.reactions().len()],
            core_iterations,
        }
    }
    fn evaluate(&mut self, dose: f64, work: &mut Work<'_, impl FnMut() -> bool>) -> Result<f64> {
        work.tick()?;
        self.sums.fill(Sum::default());
        for (i, p) in self.batch.layout.populations.iter().enumerate() {
            let core = self.batch.layout.curves[i].advance(
                self.batch.penetrations[i],
                if p.original_amount.value() == 0.0 {
                    0.0
                } else {
                    dose
                },
                self.core_iterations,
                &mut || work.tick(),
            )?;
            if p.original_amount.value() > 0.0
                && dose > 0.0
                && self.batch.penetrations[i] < 1.0
                && core.conversion_increment <= 0.0
            {
                return Err(PopulationBatchError::Arithmetic(
                    "positive population progress disappeared",
                ));
            }
            let xi = extent(p, conversion(core.penetration))?;
            self.sums[self.batch.layout.slots[i]].add(xi)?;
            self.progress[i] = core;
        }
        let mut reagent = Sum::default();
        reagent.add(self.batch.layout.initial_reagent)?;
        for (i, reaction) in self.batch.layout.network.reactions().iter().enumerate() {
            let xi = self.sums[i].value()?;
            self.extents[i] = xi;
            let a = -(reaction
                .terms()
                .iter()
                .find(|t| t.species == self.batch.layout.reagent)
                .expect("shared reactant")
                .coefficient as f64);
            let used = a * xi;
            if !used.is_finite() || (xi > 0.0 && used == 0.0) {
                return Err(PopulationBatchError::Arithmetic(
                    "shared reagent demand loses range",
                ));
            }
            reagent.add(-used)?;
            reagent.add(-a.mul_add(xi, -used))?;
        }
        let amount = reagent.value()?;
        let concentration = amount / self.batch.layout.volume.value();
        if !concentration.is_finite() || (amount != 0.0 && concentration == 0.0) {
            return Err(PopulationBatchError::Arithmetic(
                "shared reagent concentration loses range",
            ));
        }
        Ok(concentration)
    }
}

fn extent(p: &LeachPopulation, x: f64) -> Result<f64> {
    let amount = p.original_amount.value() * x;
    let result = amount / p.model.coefficients[0] as f64;
    if !result.is_finite() || (p.original_amount.value() > 0.0 && x > 0.0 && result == 0.0) {
        return Err(PopulationBatchError::Arithmetic(
            "population extent loses range",
        ));
    }
    Ok(result)
}
fn stock(inventory: &[(SpeciesId, Amount)], species: &SpeciesId) -> f64 {
    inventory
        .iter()
        .find(|(id, _)| id == species)
        .expect("admitted network inventory")
        .1
        .value()
}

struct Work<'a, F: FnMut() -> bool> {
    used: usize,
    maximum: usize,
    integrands: usize,
    cancelled: &'a mut F,
}
impl<F: FnMut() -> bool> Work<'_, F> {
    fn tick(&mut self) -> core::result::Result<(), ShrinkingCoreError> {
        if (self.cancelled)() {
            return Err(ShrinkingCoreError::Cancelled);
        }
        if self.used >= self.maximum {
            return Err(ShrinkingCoreError::Budget);
        }
        self.used += 1;
        Ok(())
    }
}
#[derive(Clone, Copy, Default)]
struct Sum {
    sum: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) -> Result<()> {
        let next = self.sum + value;
        self.correction += if self.sum.abs() >= value.abs() {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        self.sum = next;
        self.value().map(|_| ())
    }
    fn value(self) -> Result<f64> {
        let value = self.sum + self.correction;
        if value.is_finite() {
            Ok(value)
        } else {
            Err(PopulationBatchError::Arithmetic(
                "population accumulation loses range",
            ))
        }
    }
}

#[cfg(test)]
#[path = "population/tests.rs"]
mod tests;
