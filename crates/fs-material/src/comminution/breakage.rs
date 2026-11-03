//! Source-declared, congruent size redistribution for one pass through a mill.
//!
//! Each column supplies fractions of one input class's dry mass entering named
//! output classes. Positive off-diagonal transfers must strictly decrease size.
//! The same coefficients carry every explicitly assayed element; this is size
//! redistribution, not mineral liberation, selective chemistry, or a predicted
//! breakage law. Water passes through unchanged. No coefficients are inferred
//! from Bond energy, mineral names, residence time, or a population-balance fit.
//!
//! Original fractions are retained without renormalization. Column closure is
//! admitted within eight binary64 epsilons; resulting conservation has that
//! source-rounding allowance plus arithmetic rounding. Compensated sums retain
//! small contributions, and positive transfer underflow refuses atomically.

use crate::separation::{ParticleClass, ParticleStream, SeparationError};
use fs_evidence::ValidityDomain;
use fs_matdb::{Provenance, QueryPoint};
use fs_qty::{Dimensionless, Length, MassFlowRate};
use std::collections::BTreeMap;

/// Maximum declared size classes, bounding the per-pass matrix work.
pub const MAX_SIZE_CLASSES: usize = 128;

/// One source-declared particle class; identity is distinct from its diameter.
#[derive(Clone, Debug, PartialEq)]
pub struct SizeClass {
    /// Unique, nonempty class identity, preserved in the outgoing stream.
    pub id: String,
    /// Positive representative diameter in metres.
    pub diameter: Length,
}

/// One original entry in a dry-mass redistribution column.
#[derive(Clone, Debug, PartialEq)]
pub struct BreakageTransfer {
    /// Identity of a declared output class.
    pub output: String,
    /// Fraction of this column's input dry mass, in [0,1].
    pub mass_fraction: Dimensionless,
}

/// Sparse source column: omitted output classes receive exactly zero transfer.
#[derive(Clone, Debug, PartialEq)]
pub struct BreakageColumn {
    /// Identity of the input class whose complete dry mass this column allocates.
    pub input: String,
    /// Unique known outputs; fractions sum to one within the admission allowance.
    pub transfers: Vec<BreakageTransfer>,
}

/// Source of the empirical per-pass coefficients and their applicability.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BreakageTestContext {
    /// Original ore/composite and class-layout identity.
    pub specimen: String,
    /// Test procedure, preparation, pass definition, and operating regime.
    pub method: String,
    /// Original data source and license; no calibration is manufactured.
    pub provenance: Provenance,
}

/// Invalid source matrix, incompatible feed, unsupported conditions, or range loss.
#[derive(Clone, Debug, PartialEq)]
pub enum BreakageError {
    /// An input or source declaration violates the stated matrix contract.
    Invalid(
        /// The failed admission rule.
        &'static str,
    ),
    /// At least one declared source condition is missing, mistyped, or unsupported.
    OutsideValidity,
    /// A positive transfer or accumulated flow cannot be represented.
    Arithmetic(
        /// The failed numerical operation.
        &'static str,
    ),
    /// The original particle-stream owner refused a proposed output.
    Stream(
        /// Original stream admission failure, without replacing its meaning.
        SeparationError,
    ),
}

impl core::fmt::Display for BreakageError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "size redistribution: {self:?}")
    }
}
impl std::error::Error for BreakageError {}
impl From<SeparationError> for BreakageError {
    fn from(error: SeparationError) -> Self {
        Self::Stream(error)
    }
}

/// Immutable, source-bound redistribution matrix with canonical class identities.
#[derive(Clone, Debug, PartialEq)]
pub struct SizeRedistributionCard {
    classes: Vec<SizeClass>,
    columns: Vec<BreakageColumn>,
    context: BreakageTestContext,
    validity: ValidityDomain,
}

impl SizeRedistributionCard {
    /// Admit exactly one conservative column per class. Positive transfers may
    /// retain the same ID or enter a strictly smaller diameter; equally sized
    /// classes cannot exchange mass. Explicit zero entries carry no material.
    /// Every declared source-condition axis requires finite bounds and an exact
    /// quantity convention. An unconstrained domain is an explicit source claim.
    pub fn new(
        mut classes: Vec<SizeClass>,
        mut columns: Vec<BreakageColumn>,
        context: BreakageTestContext,
        validity: ValidityDomain,
    ) -> Result<Self, BreakageError> {
        require(
            !classes.is_empty() && classes.len() <= MAX_SIZE_CLASSES,
            "declare between one and 128 size classes",
        )?;
        require(
            !context.specimen.trim().is_empty()
                && !context.method.trim().is_empty()
                && !context.provenance.source.trim().is_empty()
                && !context.provenance.license.trim().is_empty(),
            "specimen, method, source and license are required",
        )?;
        for (axis, &(lo, hi)) in validity.bounds() {
            require(
                !axis.trim().is_empty()
                    && lo.is_finite()
                    && hi.is_finite()
                    && lo <= hi
                    && validity.axis_quantities().contains_key(axis),
                "source conditions require named finite bounds and explicit quantities",
            )?;
        }
        classes.sort_by(|a, b| a.id.cmp(&b.id));
        for class in &classes {
            require(
                !class.id.trim().is_empty() && !class.id.chars().any(char::is_control),
                "invalid size-class identity",
            )?;
            require(
                class.diameter.value().is_finite() && class.diameter.value() > 0.0,
                "class diameter must be finite and positive",
            )?;
        }
        require(
            classes.windows(2).all(|pair| pair[0].id != pair[1].id),
            "duplicate size-class identity",
        )?;
        require(
            columns.len() == classes.len(),
            "supply exactly one column for every input class",
        )?;
        columns.sort_by(|a, b| a.input.cmp(&b.input));
        for (input, column) in classes.iter().zip(&mut columns) {
            require(
                input.id == column.input,
                "missing, duplicate or unknown input column",
            )?;
            require(
                !column.transfers.is_empty() && column.transfers.len() <= classes.len(),
                "invalid transfer count",
            )?;
            column.transfers.sort_by(|a, b| a.output.cmp(&b.output));
            require(
                column
                    .transfers
                    .windows(2)
                    .all(|pair| pair[0].output != pair[1].output),
                "duplicate output transfer",
            )?;
            let mut total = Accumulator::default();
            for transfer in &column.transfers {
                let index = classes
                    .binary_search_by(|class| class.id.cmp(&transfer.output))
                    .map_err(|_| BreakageError::Invalid("unknown output class"))?;
                let fraction = transfer.mass_fraction.value();
                require(
                    fraction.is_finite() && (0.0..=1.0).contains(&fraction),
                    "mass fractions must be finite and in [0,1]",
                )?;
                require(
                    fraction == 0.0
                        || transfer.output == input.id
                        || classes[index].diameter < input.diameter,
                    "positive transfer must retain identity or strictly reduce diameter",
                )?;
                total.add(fraction)?;
            }
            require(
                (total.value()? - 1.0).abs() <= 8.0 * f64::EPSILON,
                "each input column must conserve dry mass",
            )?;
        }
        Ok(Self {
            classes,
            columns,
            context,
            validity,
        })
    }

    /// Original class layout in canonical identity order, not diameter order.
    pub fn classes(&self) -> &[SizeClass] {
        &self.classes
    }
    /// Original fractions, with input columns and output entries sorted by ID.
    pub fn columns(&self) -> &[BreakageColumn] {
        &self.columns
    }
    /// Original test context and source provenance.
    pub fn test_context(&self) -> &BreakageTestContext {
        &self.context
    }
    /// Original source support; queries never widen it.
    pub fn validity(&self) -> &ValidityDomain {
        &self.validity
    }

    /// Read a source coefficient. Omitted entries between known classes mean
    /// zero; unknown input or output identities refuse. This getter does not
    /// admit operating conditions: consumers must still call `transform`.
    pub fn fraction(&self, output: &str, input: &str) -> Result<Dimensionless, BreakageError> {
        self.class_index(output)?;
        let column = &self.columns[self.class_index(input)?];
        Ok(column
            .transfers
            .binary_search_by(|entry| entry.output.as_str().cmp(output))
            .map(|index| column.transfers[index].mass_fraction)
            .unwrap_or(Dimensionless::new(0.0)))
    }

    /// Original diagonal coefficient; no complement is inferred or rounded away.
    pub fn retained_fraction(&self, input: &str) -> Result<Dimensionless, BreakageError> {
        self.fraction(input, input)
    }

    /// Independently sum the off-diagonal coefficients, preserving tiny breakage
    /// even when the retained coefficient rounds to one. The column's original
    /// closure allowance applies; this is never replaced by `1 - retained`.
    pub fn fragmented_fraction(&self, input: &str) -> Result<Dimensionless, BreakageError> {
        let column = &self.columns[self.class_index(input)?];
        let mut total = Accumulator::default();
        for transfer in &column.transfers {
            if transfer.output != input {
                total.add(transfer.mass_fraction.value())?;
            }
        }
        Ok(Dimensionless::new(total.value()?))
    }

    /// Apply this pass to an exactly matching class layout and source conditions.
    /// Every transferred element follows its input dry mass congruently; no
    /// unknown assay becomes zero. Water and the original input remain unchanged.
    /// Zero-flow output classes retain their explicit size/element layout.
    pub fn transform(
        &self,
        feed: &ParticleStream,
        conditions: &QueryPoint,
    ) -> Result<ParticleStream, BreakageError> {
        if !self
            .validity
            .contains_typed(conditions.axes(), conditions.axis_quantities())
        {
            return Err(BreakageError::OutsideValidity);
        }
        require(
            feed.classes().len() == self.classes.len()
                && feed
                    .classes()
                    .iter()
                    .zip(&self.classes)
                    .all(|(actual, expected)| {
                        actual.id() == expected.id && actual.diameter() == expected.diameter
                    }),
            "feed identities and diameters must match the source class layout exactly",
        )?;
        let elements: Vec<_> = feed.classes()[0].elements().keys().cloned().collect();
        let mut totals = vec![vec![Accumulator::default(); elements.len() + 1]; self.classes.len()];
        for (input, column) in feed.classes().iter().zip(&self.columns) {
            for transfer in &column.transfers {
                let fraction = transfer.mass_fraction.value();
                let output = self.class_index(&transfer.output)?;
                totals[output][0].add(product(input.dry_solids().value(), fraction)?)?;
                for (index, element) in elements.iter().enumerate() {
                    totals[output][index + 1]
                        .add(product(input.elements()[element].value(), fraction)?)?;
                }
            }
        }
        let mut outputs = Vec::with_capacity(self.classes.len());
        for (class, amounts) in self.classes.iter().zip(totals) {
            let mut assays = BTreeMap::new();
            for (element, amount) in elements.iter().zip(&amounts[1..]) {
                assays.insert(element.clone(), MassFlowRate::new(amount.value()?));
            }
            outputs.push(ParticleClass::new(
                class.id.clone(),
                class.diameter,
                MassFlowRate::new(amounts[0].value()?),
                assays,
            )?);
        }
        Ok(ParticleStream::new(outputs, feed.water())?)
    }

    fn class_index(&self, id: &str) -> Result<usize, BreakageError> {
        self.classes
            .binary_search_by(|class| class.id.as_str().cmp(id))
            .map_err(|_| BreakageError::Invalid("unknown size-class identity"))
    }
}

fn require(condition: bool, reason: &'static str) -> Result<(), BreakageError> {
    if condition {
        Ok(())
    } else {
        Err(BreakageError::Invalid(reason))
    }
}

fn product(flow: f64, fraction: f64) -> Result<f64, BreakageError> {
    let value = flow * fraction;
    if !value.is_finite() || (flow > 0.0 && fraction > 0.0 && value == 0.0) {
        Err(BreakageError::Arithmetic(
            "positive transfer is unrepresentable",
        ))
    } else {
        Ok(value)
    }
}

#[derive(Clone, Copy, Default)]
struct Accumulator {
    sum: f64,
    correction: f64,
}

impl Accumulator {
    fn add(&mut self, value: f64) -> Result<(), BreakageError> {
        let next = self.sum + value;
        if !next.is_finite() {
            return Err(BreakageError::Arithmetic("accumulated flow overflows"));
        }
        self.correction += if self.sum >= value {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        self.sum = next;
        Ok(())
    }
    fn value(self) -> Result<f64, BreakageError> {
        let value = self.sum + self.correction;
        if value.is_finite() {
            Ok(value)
        } else {
            Err(BreakageError::Arithmetic("compensated total overflows"))
        }
    }
}

#[cfg(test)]
#[path = "breakage/tests.rs"]
mod tests;
