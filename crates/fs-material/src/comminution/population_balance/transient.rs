//! Finite inventory in the original continuous perfectly mixed grinding volume.
//!
//! For class mass H_i, constant feed f_i and common residence tau, conservation
//! gives `dH_i/dt = f_i - (S_i + 1/tau) H_i + sum_j b_ij S_j H_j`.
//! This adds accumulation and flow to the original selection/breakage law; its
//! zero derivative is exactly the parent module's steady equation. The law and
//! assumptions are those cited there, not a measured transient calibration.
//! Every assayed element follows its parent dry mass, and water obeys
//! `dH_water/dt = f_water - H_water/tau` without breakage.
//!
//! Backward Euler uses a positive triangular solve at each declared substep.
//! Integrated discharge is the same substep's `dt * H_new/tau`, so inventory,
//! input and output close without replacing an account by a balance remainder.
//! This is a first-order time approximation, not an exact trajectory or error
//! enclosure. The caller supplies its step and work budget. Rates, daughters,
//! feed, residence and source conditions remain constant during an interval;
//! feed changes are separate calls carrying the returned inventory forward.
//!
//! This local model has no dynamic classifier/recycle, variable discharge law,
//! mineral liberation, rate-power calibration, Bond transient power estimate or
//! plant transfer claim. All physical inputs retain their original meaning.

use super::{
    ClassHoldup, ContinuousPopulationBalanceMill, MAX_SIZE_CLASSES, PopulationBalanceError, Sum,
    product, quotient, require,
};
use crate::separation::{ParticleClass, ParticleStream};
use fs_matdb::QueryPoint;
use fs_qty::{Dimensionless, Mass, MassFlowRate, Time, chemistry::ElementId};
use std::collections::BTreeMap;

/// Admission, numerical work, cancellation or a failed extensive balance.
#[derive(Clone, Debug, PartialEq)]
pub enum TransientMillError {
    /// Original continuous-mill or stream owner refused an input/account.
    Population(PopulationBalanceError),
    /// Requested interval needs more substeps than the explicit budget permits.
    Budget,
    /// Cancellation was observed; no partial inventory is published.
    Cancelled,
    /// A final unmodified extensive account fails the caller's tolerance.
    Balance {
        /// Dry solids, water, or an explicit assayed element identity.
        account: String,
        /// Initial inventory plus feed minus final inventory and discharge, kg.
        residual: Mass,
        /// Absolute plus relative tolerance at the actual account scale, kg.
        allowed: Mass,
    },
}
impl core::fmt::Display for TransientMillError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "transient grinding: {self:?}")
    }
}
impl std::error::Error for TransientMillError {}
impl From<PopulationBalanceError> for TransientMillError {
    fn from(value: PopulationBalanceError) -> Self {
        Self::Population(value)
    }
}
type Result<T> = core::result::Result<T, TransientMillError>;

/// Size-resolved extensive masses, also used for integrated input/discharge.
/// Element axes are explicit and identical in every class, including zeroes.
#[derive(Clone, Debug, PartialEq)]
pub struct MillInventory {
    classes: Vec<ClassHoldup>,
    water: Mass,
}
impl MillInventory {
    /// Admit original class masses in kg without normalization or inferred assays.
    pub fn new(mut classes: Vec<ClassHoldup>, water: Mass) -> Result<Self> {
        require(
            !classes.is_empty() && classes.len() <= MAX_SIZE_CLASSES,
            "inventory requires between one and 128 size classes",
        )?;
        classes.sort_by(|a, b| a.id.cmp(&b.id));
        let inventory = Self { classes, water };
        // A one-second drain delegates the identical nonnegative mass/assay
        // constraints to the established stream owner with coherent kg/s units.
        inventory.discharge_at(Time::new(1.0))?;
        inventory.dry_solids()?;
        Ok(inventory)
    }

    /// Empty inventory retaining the feed's explicit class layout and assay axes.
    pub fn empty_like(feed: &ParticleStream) -> Result<Self> {
        Self::integrate(feed, Time::new(0.0))
    }
    /// Original class identities, sizes and masses in canonical identity order.
    pub fn classes(&self) -> &[ClassHoldup] {
        &self.classes
    }
    /// Total carrier-water mass, kg; no liquid density or dissolved assay inferred.
    pub const fn water(&self) -> Mass {
        self.water
    }
    /// Sum of all dry class masses, including unassayed constituents, kg.
    pub fn dry_solids(&self) -> Result<Mass> {
        Ok(Mass::new(sum(self
            .classes
            .iter()
            .map(|c| c.dry_solids.value()))?))
    }
    /// Total mass of an explicitly assayed element; missing is not measured zero.
    pub fn element_mass(&self, element: &ElementId) -> Result<Mass> {
        require(
            self.classes[0].elements.contains_key(element),
            "element is not assayed in this inventory",
        )?;
        Ok(Mass::new(sum(self
            .classes
            .iter()
            .map(|c| c.elements[element].value()))?))
    }

    fn integrate(stream: &ParticleStream, elapsed: Time) -> Result<Self> {
        let classes = stream
            .classes()
            .iter()
            .map(|c| {
                Ok(ClassHoldup {
                    id: c.id().into(),
                    diameter: c.diameter(),
                    dry_solids: Mass::new(product(c.dry_solids().value(), elapsed.value())?),
                    elements: c
                        .elements()
                        .iter()
                        .map(|(id, flow)| {
                            Ok((
                                id.clone(),
                                Mass::new(product(flow.value(), elapsed.value())?),
                            ))
                        })
                        .collect::<Result<_>>()?,
                })
            })
            .collect::<Result<_>>()?;
        Self::new(
            classes,
            Mass::new(product(stream.water().value(), elapsed.value())?),
        )
    }

    fn discharge_at(&self, residence: Time) -> Result<ParticleStream> {
        let classes = self
            .classes
            .iter()
            .map(|c| {
                ParticleClass::new(
                    c.id.clone(),
                    c.diameter,
                    MassFlowRate::new(quotient(c.dry_solids.value(), residence.value())?),
                    c.elements
                        .iter()
                        .map(|(id, mass)| {
                            Ok((
                                id.clone(),
                                MassFlowRate::new(quotient(mass.value(), residence.value())?),
                            ))
                        })
                        .collect::<core::result::Result<_, PopulationBalanceError>>()?,
                )
                .map_err(PopulationBalanceError::from)
            })
            .collect::<core::result::Result<_, PopulationBalanceError>>()?;
        Ok(ParticleStream::new(
            classes,
            MassFlowRate::new(quotient(self.water.value(), residence.value())?),
        )
        .map_err(PopulationBalanceError::from)?)
    }
}

/// Explicit numerical resolution, work and conservation admission for one interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransientMillSettings {
    /// Largest permitted implicit substep, seconds; no adaptive accuracy claim.
    pub maximum_step: Time,
    /// Total substep budget in 1..=1_000_000, checked before advancing inventory.
    pub maximum_steps: usize,
    /// Absolute closure tolerance for every extensive account, kg.
    pub absolute_balance_tolerance: Mass,
    /// Relative closure tolerance, applied to each account's own magnitude.
    pub relative_balance_tolerance: Dimensionless,
}

/// Signed unmodified residuals: initial + integrated feed - final - discharge.
#[derive(Clone, Debug, PartialEq)]
pub struct TransientMillBalance {
    /// Dry mass closure, kg; individual classes change through breakage.
    pub dry_solids: Mass,
    /// Independent carrier-water closure, kg.
    pub water: Mass,
    /// Closure of every explicit solid-phase elemental mass, kg.
    pub elements: BTreeMap<ElementId, Mass>,
}

/// Accepted finite-time response retaining the original source and input state.
#[derive(Clone, Debug, PartialEq)]
pub struct TransientMillStep {
    model: ContinuousPopulationBalanceMill,
    initial: MillInventory,
    feed: ParticleStream,
    /// Actual source query, including feed flow and common residence.
    pub conditions: QueryPoint,
    /// Residence shared by dry classes and water for the complete interval.
    pub residence: Time,
    /// Advanced time, seconds; not the residence or an individual particle age.
    pub elapsed: Time,
    /// Final retained inventory; pass this into the next feed interval.
    pub inventory: MillInventory,
    /// Final instantaneous discharge rates, kg/s, obtained from final inventory.
    pub outlet: ParticleStream,
    /// Total feed mass entering over the interval, kg by class and element.
    pub integrated_feed: MillInventory,
    /// Total mass discharged over all substeps, kg by class and element.
    pub integrated_discharge: MillInventory,
    /// Actual extensive residuals; no account is replaced to make them vanish.
    pub balance: TransientMillBalance,
    /// Declared work, resolution and conservation limits.
    pub settings: TransientMillSettings,
    /// Completed equal nominal substeps (last step closes the interval). At most
    /// one extra step is used when rounding would otherwise exceed maximum_step.
    pub substeps: usize,
}
impl TransientMillStep {
    /// Original local source, including unmodified rates, daughters and validity.
    pub fn model(&self) -> &ContinuousPopulationBalanceMill {
        &self.model
    }
    /// Original initial inventory; failures never mutate the caller's state.
    pub fn initial(&self) -> &MillInventory {
        &self.initial
    }
    /// Original constant feed for this interval.
    pub fn feed(&self) -> &ParticleStream {
        &self.feed
    }
}

impl ContinuousPopulationBalanceMill {
    /// Advance one constant-input interval with positive backward-Euler substeps.
    ///
    /// Actual feed coordinates must satisfy the original source. Empty-feed
    /// washout requires zero-flow source support; a source requiring feed solids
    /// fraction refuses a completely empty feed because that ratio is undefined.
    /// Existing and incoming inventories must have identical class/assay axes.
    /// Cancellation is checked before work, at each class tile and before return.
    /// No refusal returns a partial step or changes either input.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_inventory(
        &self,
        initial: &MillInventory,
        feed: &ParticleStream,
        residence: Time,
        elapsed: Time,
        conditions: &QueryPoint,
        settings: TransientMillSettings,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<TransientMillStep> {
        poll(&mut cancelled)?;
        let substeps = admit_settings(elapsed, settings)?;
        let actual = self.admitted_feed(feed, residence, conditions, true)?;
        require(
            initial.classes.len() == feed.classes().len()
                && initial.classes.iter().zip(feed.classes()).all(|(a, b)| {
                    a.id == b.id()
                        && a.diameter == b.diameter()
                        && a.elements.keys().eq(b.elements().keys())
                }),
            "inventory and feed require identical size identities and explicit assays",
        )?;
        let elements: Vec<_> = initial.classes[0].elements.keys().cloned().collect();
        let channels = elements.len() + 1;
        let mut masses: Vec<Vec<f64>> = initial
            .classes
            .iter()
            .map(|c| {
                core::iter::once(c.dry_solids.value())
                    .chain(elements.iter().map(|e| c.elements[e].value()))
                    .collect()
            })
            .collect();
        let rates: Vec<Vec<f64>> = feed
            .classes()
            .iter()
            .map(|c| {
                core::iter::once(c.dry_solids().value())
                    .chain(elements.iter().map(|e| c.elements()[e].value()))
                    .collect()
            })
            .collect();
        let mut order: Vec<_> = (0..self.classes.len()).collect();
        order.sort_by(|&a, &b| {
            self.classes[b]
                .diameter
                .value()
                .total_cmp(&self.classes[a].diameter.value())
                .then_with(|| self.classes[a].id.cmp(&self.classes[b].id))
        });
        let mut discharged = vec![vec![Sum::default(); channels]; self.classes.len()];
        let mut water = initial.water.value();
        let mut discharged_water = Sum::default();
        let nominal_step = quotient(elapsed.value(), substeps as f64)?;
        let integrated_feed = MillInventory::integrate(feed, elapsed)?;
        for step in 0..substeps {
            poll(&mut cancelled)?;
            let dt = if step + 1 == substeps {
                elapsed.value() - nominal_step * step as f64
            } else {
                nominal_step
            };
            require(
                dt.is_finite() && dt > 0.0 && dt <= settings.maximum_step.value(),
                "requested time partition is not representable within the step limit",
            )?;
            let flushed = quotient(dt, residence.value())?;
            let mut incoming = vec![vec![Sum::default(); channels]; self.classes.len()];
            for (i, row) in incoming.iter_mut().enumerate() {
                for (channel, value) in row.iter_mut().enumerate() {
                    value.add(masses[i][channel])?;
                    value.add(product(dt, rates[i][channel])?)?;
                }
            }
            for &i in &order {
                poll(&mut cancelled)?;
                let selected = product(dt, self.kinetics[i].selection_rate.value())?;
                let denominator = sum([1.0, flushed, selected])?;
                for channel in 0..channels {
                    let mass = quotient(incoming[i][channel].value()?, denominator)?;
                    masses[i][channel] = mass;
                    discharged[i][channel].add(product(mass, flushed)?)?;
                    let broken = product(mass, selected)?;
                    for daughter in &self.kinetics[i].daughters {
                        let output = self
                            .classes
                            .binary_search_by(|c| c.id.cmp(&daughter.output))
                            .expect("original source admits every daughter identity");
                        incoming[output][channel]
                            .add(product(broken, daughter.mass_fraction.value())?)?;
                    }
                }
            }
            water = quotient(
                sum([water, product(dt, feed.water().value())?])?,
                sum([1.0, flushed])?,
            )?;
            discharged_water.add(product(water, flushed)?)?;
        }
        let inventory = inventory_from_channels(self, &elements, masses, water)?;
        let integrated_discharge = inventory_from_channels(
            self,
            &elements,
            discharged
                .into_iter()
                .map(|row| row.into_iter().map(Sum::value).collect())
                .collect::<core::result::Result<_, PopulationBalanceError>>()?,
            discharged_water.value()?,
        )?;
        let balance = close_accounts(
            initial,
            &integrated_feed,
            &inventory,
            &integrated_discharge,
            settings,
        )?;
        let outlet = inventory.discharge_at(residence)?;
        poll(&mut cancelled)?;
        Ok(TransientMillStep {
            model: self.clone(),
            initial: initial.clone(),
            feed: feed.clone(),
            conditions: actual,
            residence,
            elapsed,
            inventory,
            outlet,
            integrated_feed,
            integrated_discharge,
            balance,
            settings,
            substeps,
        })
    }
}

fn admit_settings(elapsed: Time, settings: TransientMillSettings) -> Result<usize> {
    require(
        elapsed.value().is_finite()
            && elapsed.value() > 0.0
            && settings.maximum_step.value().is_finite()
            && settings.maximum_step.value() > 0.0
            && (1..=1_000_000).contains(&settings.maximum_steps)
            && settings.absolute_balance_tolerance.value().is_finite()
            && settings.absolute_balance_tolerance.value() >= 0.0
            && settings.relative_balance_tolerance.value().is_finite()
            && settings.relative_balance_tolerance.value() >= 0.0,
        "require positive finite times, a work budget and nonnegative finite balance tolerances",
    )?;
    let needed = (elapsed.value() / settings.maximum_step.value())
        .ceil()
        .max(1.0);
    if !needed.is_finite() || needed > settings.maximum_steps as f64 || needed >= usize::MAX as f64
    {
        return Err(TransientMillError::Budget);
    }
    let mut count = needed as usize;
    let nominal = quotient(elapsed.value(), count as f64)?;
    let last = elapsed.value() - nominal * (count - 1) as f64;
    if last > settings.maximum_step.value() {
        count += 1;
        if count > settings.maximum_steps {
            return Err(TransientMillError::Budget);
        }
    }
    Ok(count)
}

fn inventory_from_channels(
    model: &ContinuousPopulationBalanceMill,
    elements: &[ElementId],
    values: Vec<Vec<f64>>,
    water: f64,
) -> Result<MillInventory> {
    MillInventory::new(
        model
            .classes
            .iter()
            .zip(values)
            .map(|(c, row)| ClassHoldup {
                id: c.id.clone(),
                diameter: c.diameter,
                dry_solids: Mass::new(row[0]),
                elements: elements
                    .iter()
                    .cloned()
                    .zip(row[1..].iter().copied().map(Mass::new))
                    .collect(),
            })
            .collect(),
        Mass::new(water),
    )
}

fn close_accounts(
    initial: &MillInventory,
    feed: &MillInventory,
    final_state: &MillInventory,
    discharge: &MillInventory,
    settings: TransientMillSettings,
) -> Result<TransientMillBalance> {
    let inventories = [initial, feed, final_state, discharge];
    let close = |account: String, values: [Mass; 4]| -> Result<Mass> {
        let input = sum([values[0].value(), values[1].value()])?;
        let output = sum([values[2].value(), values[3].value()])?;
        let residual = input - output;
        let allowed = settings.absolute_balance_tolerance.value()
            + settings.relative_balance_tolerance.value() * input.max(output);
        require(allowed.is_finite(), "balance tolerance overflows")?;
        if residual.abs() > allowed {
            return Err(TransientMillError::Balance {
                account,
                residual: Mass::new(residual),
                allowed: Mass::new(allowed),
            });
        }
        Ok(Mass::new(residual))
    };
    let dry_solids = close(
        "dry-solids".into(),
        [
            initial.dry_solids()?,
            feed.dry_solids()?,
            final_state.dry_solids()?,
            discharge.dry_solids()?,
        ],
    )?;
    let water = close("water".into(), inventories.map(|i| i.water))?;
    let mut elements = BTreeMap::new();
    for element in initial.classes[0].elements.keys() {
        elements.insert(
            element.clone(),
            close(
                format!("element:{element:?}"),
                [
                    initial.element_mass(element)?,
                    feed.element_mass(element)?,
                    final_state.element_mass(element)?,
                    discharge.element_mass(element)?,
                ],
            )?,
        );
    }
    Ok(TransientMillBalance {
        dry_solids,
        water,
        elements,
    })
}

fn sum(values: impl IntoIterator<Item = f64>) -> core::result::Result<f64, PopulationBalanceError> {
    let mut total = Sum::default();
    for value in values {
        total.add(value)?;
    }
    total.value()
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(TransientMillError::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
#[path = "transient/tests.rs"]
mod tests;
