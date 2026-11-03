//! Steady size-mass balance in one continuous, perfectly mixed grinding volume.
//!
//! The source declares a selection rate S_i [1/s] and the daughter dry-mass
//! fractions b_ij for each class. With a common residence time tau, discharge
//! rates satisfy p_i = (f_i + tau sum_j b_ij S_j p_j)/(1 + tau S_i).
//! Only strictly coarser j contribute, so the solve is triangular. This is the
//! continuous perfect-mixer balance in Minerals 12 (2022), 679, equation 2:
//! <https://doi.org/10.3390/min12060679>. The reference supplies the model
//! family, not rates or a calibration for any ore or mill.
//!
//! Every declared solid-phase element follows its parent mass congruently.
//! All size classes and carrier water share tau; class and water hold-ups are
//! their discharge rates times tau under this explicit perfect-mixing claim.
//! Daughter fractions are mass fractions, not particle-number fractions.
//! The finest classes have zero selection and collect all declared fragments.
//! A per-pass redistribution matrix is not silently converted into rates.
//! Daughter sums admit eight binary64 epsilons of source rounding, as in the
//! existing redistribution owner; account arithmetic is not an exact certificate.
//!
//! Rates and daughters retain their common original test context and validity.
//! They are held constant at the admitted operating point. No liberation,
//! agglomeration, size-dependent discharge, transient response, media trajectory,
//! power draw, wear, or transfer of laboratory rates to a plant is inferred.

use super::breakage::{BreakageTestContext, BreakageTransfer, MAX_SIZE_CLASSES, SizeClass};
pub use super::circuit::MILL_DRY_FLOW_AXIS;
use crate::separation::{ParticleClass, ParticleStream, SOLIDS_FRACTION_AXIS, SeparationError};
use fs_evidence::ValidityDomain;
use fs_matdb::{MatDbError, QueryPoint};
use fs_qty::{
    Dimensionless, Frequency, Length, Mass, MassFlowRate, QuantitySpec, Time, chemistry::ElementId,
};
use std::collections::BTreeMap;

/// Actual common mean residence, in seconds, rather than an individual age.
pub const MEAN_RESIDENCE_TIME_AXIS: &str = "mean-residence-time";

/// One original selection rate and its conditional daughter mass distribution.
#[derive(Clone, Debug, PartialEq)]
pub struct SelectionBreakageClass {
    /// Exactly one declaration is required per source size class.
    pub input: String,
    /// Nonnegative first-order selection frequency, in inverse seconds.
    pub selection_rate: Frequency,
    /// Positive-rate classes distribute all selected mass to smaller classes.
    /// Zero-rate classes have no daughters. Original fractions are not rescaled.
    pub daughters: Vec<BreakageTransfer>,
}

/// Invalid source, unsupported operating point, or unrepresentable account.
#[derive(Clone, Debug, PartialEq)]
pub enum PopulationBalanceError {
    /// Missing, inconsistent or nonphysical source/input declaration.
    Invalid(&'static str),
    /// The actual operating point leaves the original source support.
    OutsideValidity,
    /// Typed source-coordinate construction failed.
    Source(MatDbError),
    /// The existing particle-stream owner refused a material account.
    Stream(SeparationError),
    /// A required positive flow or inventory cannot be represented.
    Arithmetic(&'static str),
}
impl core::fmt::Display for PopulationBalanceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "continuous grinding: {self:?}")
    }
}
impl std::error::Error for PopulationBalanceError {}
impl From<MatDbError> for PopulationBalanceError {
    fn from(value: MatDbError) -> Self {
        Self::Source(value)
    }
}
impl From<SeparationError> for PopulationBalanceError {
    fn from(value: SeparationError) -> Self {
        Self::Stream(value)
    }
}
type Result<T> = core::result::Result<T, PopulationBalanceError>;

/// Immutable source declaration for a local continuous grinding law.
#[derive(Clone, Debug, PartialEq)]
pub struct ContinuousPopulationBalanceMill {
    classes: Vec<SizeClass>,
    kinetics: Vec<SelectionBreakageClass>,
    context: BreakageTestContext,
    validity: ValidityDomain,
}

/// Extensive in-mill accounts, in kg, with the original class identity.
#[derive(Clone, Debug, PartialEq)]
pub struct ClassHoldup {
    /// Original source class identity.
    pub id: String,
    /// Original representative particle diameter.
    pub diameter: Length,
    /// Total dry mass in this class, including unassayed constituents.
    pub dry_solids: Mass,
    /// Each explicitly assayed element's in-mill mass, without inferred assays.
    pub elements: BTreeMap<ElementId, Mass>,
}

/// Admitted steady outlet and its hold-up under the common-residence model.
#[derive(Clone, Debug, PartialEq)]
pub struct ContinuousMillSolution {
    model: ContinuousPopulationBalanceMill,
    /// Source query including the actual residence, flow and solids fraction.
    pub conditions: QueryPoint,
    /// Common residence of all solids classes and carrier water.
    pub mean_residence_time: Time,
    /// Size-resolved outgoing solids, original assayed elements and carrier water.
    pub outlet: ParticleStream,
    /// Canonical identity order, matching the outlet classes.
    pub class_holdup: Vec<ClassHoldup>,
    /// Sum of the retained class dry masses, in kg.
    pub dry_solids_holdup: Mass,
    /// Carrier water is also perfectly mixed with the same mean residence.
    pub water_holdup: Mass,
}
impl ContinuousMillSolution {
    /// Retained original rates, daughters, context and support.
    pub fn model(&self) -> &ContinuousPopulationBalanceMill {
        &self.model
    }
}

impl ContinuousPopulationBalanceMill {
    /// Admit source rates and strictly finer conservative daughter columns.
    /// The source must explicitly bound positive residence time and nonnegative
    /// dry feed rate. Every additional applicability axis must also be typed
    /// and finite. The same context applies to both rates and daughters.
    pub fn new(
        mut classes: Vec<SizeClass>,
        mut kinetics: Vec<SelectionBreakageClass>,
        context: BreakageTestContext,
        validity: ValidityDomain,
    ) -> Result<Self> {
        require(
            !classes.is_empty() && classes.len() <= MAX_SIZE_CLASSES,
            "declare between one and 128 size classes",
        )?;
        require(
            !context.specimen.trim().is_empty()
                && !context.method.trim().is_empty()
                && !context.provenance.source.trim().is_empty()
                && !context.provenance.license.trim().is_empty(),
            "rates and daughters require specimen, method, source and license",
        )?;
        for (axis, &(lo, hi)) in validity.bounds() {
            require(
                !axis.trim().is_empty()
                    && lo.is_finite()
                    && hi.is_finite()
                    && lo <= hi
                    && validity.axis_quantities().contains_key(axis),
                "source conditions require finite bounds and explicit quantities",
            )?;
        }
        for (axis, quantity, positive_lower) in [
            (
                MEAN_RESIDENCE_TIME_AXIS,
                QuantitySpec::dimensional(Time::DIMS),
                true,
            ),
            (
                MILL_DRY_FLOW_AXIS,
                QuantitySpec::dimensional(MassFlowRate::DIMS),
                false,
            ),
        ] {
            let (lo, _) = validity.bound(axis).ok_or(PopulationBalanceError::Invalid(
                "declare residence-time and dry-flow source support",
            ))?;
            require(
                validity.axis_quantities().get(axis) == Some(&quantity)
                    && if positive_lower { lo > 0.0 } else { lo >= 0.0 },
                "residence and dry-flow support require their physical units and bounds",
            )?;
        }
        if validity.bound(SOLIDS_FRACTION_AXIS).is_some() {
            require(
                validity.axis_quantities().get(SOLIDS_FRACTION_AXIS)
                    == Some(&QuantitySpec::dimensional(Dimensionless::DIMS)),
                "solids mass fraction must be dimensionless",
            )?;
        }
        classes.sort_by(|a, b| a.id.cmp(&b.id));
        for class in &classes {
            require(
                !class.id.trim().is_empty()
                    && !class.id.chars().any(char::is_control)
                    && class.diameter.value().is_finite()
                    && class.diameter.value() > 0.0,
                "source classes need unique usable identities and positive diameters",
            )?;
        }
        require(
            classes.windows(2).all(|p| p[0].id != p[1].id) && classes.len() == kinetics.len(),
            "declare one kinetic column per unique class",
        )?;
        kinetics.sort_by(|a, b| a.input.cmp(&b.input));
        for (class, column) in classes.iter().zip(&mut kinetics) {
            let rate = column.selection_rate.value();
            require(
                column.input == class.id && rate.is_finite() && rate >= 0.0,
                "each class requires exactly one finite nonnegative selection rate",
            )?;
            require(
                (rate == 0.0 && column.daughters.is_empty())
                    || (rate > 0.0
                        && !column.daughters.is_empty()
                        && column.daughters.len() < classes.len()),
                "positive rates need daughters; zero rates have none",
            )?;
            column.daughters.sort_by(|a, b| a.output.cmp(&b.output));
            require(
                column
                    .daughters
                    .windows(2)
                    .all(|p| p[0].output != p[1].output),
                "duplicate daughter class",
            )?;
            let mut total = Sum::default();
            for daughter in &column.daughters {
                let output = classes
                    .binary_search_by(|c| c.id.cmp(&daughter.output))
                    .map_err(|_| PopulationBalanceError::Invalid("unknown daughter class"))?;
                let fraction = daughter.mass_fraction.value();
                require(
                    classes[output].diameter < class.diameter
                        && fraction.is_finite()
                        && (0.0..=1.0).contains(&fraction),
                    "daughter mass fractions must enter strictly smaller classes",
                )?;
                total.add(fraction)?;
            }
            require(
                rate == 0.0 || (total.value()? - 1.0).abs() <= 8.0 * f64::EPSILON,
                "daughter mass fractions must sum to one",
            )?;
        }
        Ok(Self {
            classes,
            kinetics,
            context,
            validity,
        })
    }

    /// Source class layout in canonical identity order.
    pub fn classes(&self) -> &[SizeClass] {
        &self.classes
    }
    /// Original selection frequencies and daughter fractions, ordered by ID.
    pub fn kinetics(&self) -> &[SelectionBreakageClass] {
        &self.kinetics
    }
    /// Common original specimen, method and provenance for rates and daughters.
    pub fn test_context(&self) -> &BreakageTestContext {
        &self.context
    }
    /// Original applicability domain, never widened by evaluation.
    pub fn validity(&self) -> &ValidityDomain {
        &self.validity
    }

    /// Solve with the source-admitted mean residence in seconds. Dry solids,
    /// all measured elements and water are kept as separate extensive accounts.
    /// Caller conditions cannot replace the three internally computed axes.
    pub fn at_residence_time(
        &self,
        feed: &ParticleStream,
        residence: Time,
        conditions: &QueryPoint,
    ) -> Result<ContinuousMillSolution> {
        require(
            residence.value().is_finite() && residence.value() > 0.0,
            "mean residence must be finite and positive",
        )?;
        require(
            feed.classes().len() == self.classes.len()
                && feed
                    .classes()
                    .iter()
                    .zip(&self.classes)
                    .all(|(a, b)| a.id() == b.id && a.diameter() == b.diameter),
            "feed identities and diameters must exactly match the source layout",
        )?;
        let dry_flow = feed.dry_solids()?.value();
        require(
            dry_flow > 0.0,
            "continuous grinding requires positive dry feed",
        )?;
        for axis in [
            MEAN_RESIDENCE_TIME_AXIS,
            MILL_DRY_FLOW_AXIS,
            SOLIDS_FRACTION_AXIS,
        ] {
            require(
                !conditions.axes().contains_key(axis),
                "caller conditions cannot replace actual residence, flow or solids fraction",
            )?;
        }
        let actual = conditions
            .clone()
            .with_quantity(
                MEAN_RESIDENCE_TIME_AXIS,
                QuantitySpec::dimensional(Time::DIMS),
                residence.value(),
            )?
            .with_quantity(
                MILL_DRY_FLOW_AXIS,
                QuantitySpec::dimensional(MassFlowRate::DIMS),
                dry_flow,
            )?
            .with_quantity(
                SOLIDS_FRACTION_AXIS,
                QuantitySpec::dimensional(Dimensionless::DIMS),
                feed.solids_mass_fraction()?.value(),
            )?;
        if !self
            .validity
            .contains_typed(actual.axes(), actual.axis_quantities())
        {
            return Err(PopulationBalanceError::OutsideValidity);
        }
        let elements: Vec<_> = feed.classes()[0].elements().keys().cloned().collect();
        let channels = elements.len() + 1;
        let mut incoming = vec![vec![Sum::default(); channels]; self.classes.len()];
        let mut discharge = vec![vec![0.0; channels]; self.classes.len()];
        for (i, class) in feed.classes().iter().enumerate() {
            incoming[i][0].add(class.dry_solids().value())?;
            for (j, element) in elements.iter().enumerate() {
                incoming[i][j + 1].add(class.elements()[element].value())?;
            }
        }
        let mut order: Vec<_> = (0..self.classes.len()).collect();
        order.sort_by(|&a, &b| {
            self.classes[b]
                .diameter
                .value()
                .total_cmp(&self.classes[a].diameter.value())
                .then_with(|| self.classes[a].id.cmp(&self.classes[b].id))
        });
        for i in order {
            let selected = product(residence.value(), self.kinetics[i].selection_rate.value())?;
            let denominator = 1.0 + selected;
            for channel in 0..channels {
                discharge[i][channel] = quotient(incoming[i][channel].value()?, denominator)?;
                let broken = product(discharge[i][channel], selected)?;
                for daughter in &self.kinetics[i].daughters {
                    let output = self
                        .classes
                        .binary_search_by(|c| c.id.cmp(&daughter.output))
                        .expect("constructor admits every daughter identity");
                    incoming[output][channel]
                        .add(product(broken, daughter.mass_fraction.value())?)?;
                }
            }
        }
        let mut outlets = Vec::with_capacity(self.classes.len());
        let mut class_holdup = Vec::with_capacity(self.classes.len());
        let mut dry_holdup = Sum::default();
        for (class, amounts) in self.classes.iter().zip(discharge) {
            let assays = elements
                .iter()
                .cloned()
                .zip(amounts[1..].iter().copied().map(MassFlowRate::new))
                .collect();
            outlets.push(ParticleClass::new(
                class.id.clone(),
                class.diameter,
                MassFlowRate::new(amounts[0]),
                assays,
            )?);
            let dry_solids = product(amounts[0], residence.value())?;
            dry_holdup.add(dry_solids)?;
            let mut held_elements = BTreeMap::new();
            for (element, amount) in elements.iter().zip(&amounts[1..]) {
                held_elements.insert(
                    element.clone(),
                    Mass::new(product(*amount, residence.value())?),
                );
            }
            class_holdup.push(ClassHoldup {
                id: class.id.clone(),
                diameter: class.diameter,
                dry_solids: Mass::new(dry_solids),
                elements: held_elements,
            });
        }
        Ok(ContinuousMillSolution {
            model: self.clone(),
            conditions: actual,
            mean_residence_time: residence,
            outlet: ParticleStream::new(outlets, feed.water())?,
            class_holdup,
            dry_solids_holdup: Mass::new(dry_holdup.value()?),
            water_holdup: Mass::new(product(feed.water().value(), residence.value())?),
        })
    }

    /// Infer common residence as positive dry hold-up / actual dry feed rate.
    /// At fixed hold-up, increasing throughput shortens residence and changes
    /// the product PSD; no relation between hold-up, power or mill geometry is
    /// invented. Source support is checked at the resulting operating point.
    pub fn at_dry_holdup(
        &self,
        feed: &ParticleStream,
        dry_holdup: Mass,
        conditions: &QueryPoint,
    ) -> Result<ContinuousMillSolution> {
        require(
            dry_holdup.value().is_finite() && dry_holdup.value() > 0.0,
            "dry hold-up must be finite and positive",
        )?;
        let dry_flow = feed.dry_solids()?.value();
        require(dry_flow > 0.0, "hold-up requires positive dry feed")?;
        self.at_residence_time(
            feed,
            Time::new(quotient(dry_holdup.value(), dry_flow)?),
            conditions,
        )
    }
}

fn require(condition: bool, reason: &'static str) -> Result<()> {
    if condition {
        Ok(())
    } else {
        Err(PopulationBalanceError::Invalid(reason))
    }
}
fn product(a: f64, b: f64) -> Result<f64> {
    let value = a * b;
    if !value.is_finite() || (a > 0.0 && b > 0.0 && value == 0.0) {
        Err(PopulationBalanceError::Arithmetic(
            "positive product is unrepresentable",
        ))
    } else {
        Ok(value)
    }
}
fn quotient(a: f64, b: f64) -> Result<f64> {
    let value = a / b;
    if !value.is_finite() || (a > 0.0 && value == 0.0) {
        Err(PopulationBalanceError::Arithmetic(
            "positive quotient is unrepresentable",
        ))
    } else {
        Ok(value)
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
        if !next.is_finite() {
            return Err(PopulationBalanceError::Arithmetic("account sum overflows"));
        }
        self.correction += if self.sum >= value {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        self.sum = next;
        Ok(())
    }
    fn value(self) -> Result<f64> {
        let value = self.sum + self.correction;
        if value.is_finite() {
            Ok(value)
        } else {
            Err(PopulationBalanceError::Arithmetic("account sum overflows"))
        }
    }
}

#[cfg(test)]
#[path = "population_balance/tests.rs"]
mod tests;
