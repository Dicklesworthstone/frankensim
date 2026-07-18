//! Deterministic weighted-empirical consumption of structured backlash output.
//!
//! This module consumes an already-admitted [`StructuredPropagationReceipt`]
//! and publishes its complete weighted support, caller-requested
//! finite-population quantiles, and mode shares. It does not refit a
//! distribution, infer units, or establish gear, clearance, reliability,
//! coverage, or physical-validity authority.

use core::{cmp::Ordering, fmt};

use std::{collections::BTreeMap, num::NonZeroU64};

use crate::{
    STRUCTURED_PROPAGATION_SCHEMA_V1, StructuredLawId, StructuredModelIdentity,
    StructuredPopulationModel, StructuredPropagationReceipt,
};

/// Receipt schema for the gear-backlash consumer.
pub const GEAR_BACKLASH_CONSUMER_SCHEMA_V1: u32 = 1;
/// Maximum distinct empirical quantiles retained by one receipt.
pub const MAX_GEAR_BACKLASH_QUANTILES_V1: usize = 128;

/// Explicit caller-declared unit of structured outputs interpreted as backlash.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum GearBacklashLengthUnitV1 {
    /// Metres.
    Metre = 1,
    /// Millimetres.
    Millimetre = 2,
    /// Micrometres.
    Micrometre = 3,
    /// Nanometres.
    Nanometre = 4,
}

impl GearBacklashLengthUnitV1 {
    /// Stable receipt tag.
    #[must_use]
    pub const fn tag(self) -> u8 {
        self as u8
    }

    /// Binary64 multiplier from the submitted unit to metres.
    #[must_use]
    pub const fn metres_per_unit(self) -> f64 {
        match self {
            Self::Metre => 1.0,
            Self::Millimetre => 1.0e-3,
            Self::Micrometre => 1.0e-6,
            Self::Nanometre => 1.0e-9,
        }
    }

    /// Stable unit spelling.
    #[must_use]
    pub const fn symbol(self) -> &'static str {
        match self {
            Self::Metre => "m",
            Self::Millimetre => "mm",
            Self::Micrometre => "um",
            Self::Nanometre => "nm",
        }
    }
}

/// Refusal from constructing one exact rational probability.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GearBacklashProbabilityErrorV1 {
    /// Denominator was zero.
    ZeroDenominator,
    /// Numerator exceeded the denominator.
    AboveOne,
}

impl GearBacklashProbabilityErrorV1 {
    /// Stable diagnostic code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::ZeroDenominator => "GearBacklashProbabilityZeroDenominator",
            Self::AboveOne => "GearBacklashProbabilityAboveOne",
        }
    }
}

impl fmt::Display for GearBacklashProbabilityErrorV1 {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::ZeroDenominator => "gear-backlash probability denominator must be nonzero",
            Self::AboveOne => "gear-backlash probability must lie in the closed unit interval",
        })
    }
}

impl std::error::Error for GearBacklashProbabilityErrorV1 {}

/// Canonically reduced rational probability in the closed unit interval.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct GearBacklashProbabilityV1 {
    numerator: u64,
    denominator: NonZeroU64,
}

impl GearBacklashProbabilityV1 {
    /// Construct and reduce one exact rational probability.
    ///
    /// # Errors
    /// Refuses a zero denominator or a value above one.
    pub fn try_new(
        numerator: u64,
        denominator: u64,
    ) -> Result<Self, GearBacklashProbabilityErrorV1> {
        let Some(denominator) = NonZeroU64::new(denominator) else {
            return Err(GearBacklashProbabilityErrorV1::ZeroDenominator);
        };
        if numerator > denominator.get() {
            return Err(GearBacklashProbabilityErrorV1::AboveOne);
        }
        let divisor = greatest_common_divisor(numerator, denominator.get());
        let Some(reduced_denominator) = NonZeroU64::new(denominator.get() / divisor) else {
            return Err(GearBacklashProbabilityErrorV1::ZeroDenominator);
        };
        Ok(Self {
            numerator: numerator / divisor,
            denominator: reduced_denominator,
        })
    }

    /// Reduced numerator.
    #[must_use]
    pub const fn numerator(self) -> u64 {
        self.numerator
    }

    /// Reduced nonzero denominator.
    #[must_use]
    pub const fn denominator(self) -> NonZeroU64 {
        self.denominator
    }

    fn is_reached_by(self, cumulative_weight: u64, total_weight: u64) -> bool {
        self.numerator == 0
            || u128::from(cumulative_weight) * u128::from(self.denominator.get())
                >= u128::from(self.numerator) * u128::from(total_weight)
    }
}

/// One distinct weighted support value retained in ascending binary64 order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct GearBacklashSupportPointV1 {
    value: GearBacklashValueV1,
    relative_weight: u64,
    cumulative_before: u64,
    cumulative_at: u64,
}

impl GearBacklashSupportPointV1 {
    /// Distinct signed response value in source and coherent-SI units.
    #[must_use]
    pub const fn value(self) -> GearBacklashValueV1 {
        self.value
    }

    /// Exact multiplicity at this support value.
    #[must_use]
    pub const fn relative_weight(self) -> u64 {
        self.relative_weight
    }

    /// Exact cumulative multiplicity strictly below this support value.
    #[must_use]
    pub const fn cumulative_before(self) -> u64 {
        self.cumulative_before
    }

    /// Exact cumulative multiplicity at or below this support value.
    #[must_use]
    pub const fn cumulative_at(self) -> u64 {
        self.cumulative_at
    }
}

impl PartialOrd for GearBacklashProbabilityV1 {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for GearBacklashProbabilityV1 {
    fn cmp(&self, other: &Self) -> Ordering {
        (u128::from(self.numerator) * u128::from(other.denominator.get()))
            .cmp(&(u128::from(other.numerator) * u128::from(self.denominator.get())))
    }
}

/// Signed backlash value retaining submitted-unit and coherent-SI bits.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct GearBacklashValueV1 {
    source_bits: u64,
    unit: GearBacklashLengthUnitV1,
    metres_bits: u64,
}

impl GearBacklashValueV1 {
    /// Value in the caller-declared structured-output unit.
    #[must_use]
    pub fn source_value(self) -> f64 {
        f64::from_bits(self.source_bits)
    }

    /// Exact caller-declared unit.
    #[must_use]
    pub const fn unit(self) -> GearBacklashLengthUnitV1 {
        self.unit
    }

    /// Deterministic binary64 conversion to metres.
    #[must_use]
    pub fn metres(self) -> f64 {
        f64::from_bits(self.metres_bits)
    }
}

/// Backlash population variance retaining source-unit and square-metre bits.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct GearBacklashVarianceV1 {
    source_bits: u64,
    unit: GearBacklashLengthUnitV1,
    square_metres_bits: u64,
}

impl GearBacklashVarianceV1 {
    /// Population variance in squared caller-declared units.
    #[must_use]
    pub fn source_variance(self) -> f64 {
        f64::from_bits(self.source_bits)
    }

    /// Unit whose square applies to the source variance.
    #[must_use]
    pub const fn unit(self) -> GearBacklashLengthUnitV1 {
        self.unit
    }

    /// Deterministic binary64 conversion to square metres.
    #[must_use]
    pub fn square_metres(self) -> f64 {
        f64::from_bits(self.square_metres_bits)
    }
}

/// One left-continuous inverse-CDF result over the weighted finite population.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GearBacklashQuantileV1 {
    probability: GearBacklashProbabilityV1,
    support_index: usize,
    value: GearBacklashValueV1,
    cumulative_before: u64,
    cumulative_at: u64,
}

impl GearBacklashQuantileV1 {
    /// Exact reduced requested probability.
    #[must_use]
    pub const fn probability(&self) -> GearBacklashProbabilityV1 {
        self.probability
    }

    /// Ordinal of the selected row in the report's ascending support table.
    #[must_use]
    pub const fn support_index(&self) -> usize {
        self.support_index
    }

    /// Selected weighted empirical support value.
    #[must_use]
    pub const fn value(&self) -> GearBacklashValueV1 {
        self.value
    }

    /// Exact cumulative multiplicity strictly below the selected value.
    #[must_use]
    pub const fn cumulative_before(&self) -> u64 {
        self.cumulative_before
    }

    /// Exact cumulative multiplicity at or below the selected value.
    #[must_use]
    pub const fn cumulative_at(&self) -> u64 {
        self.cumulative_at
    }
}

/// Retained population share of one structured response mode.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GearBacklashModeShareV1 {
    law: StructuredLawId,
    law_key: Box<str>,
    piece: usize,
    mode_key: Box<str>,
    relative_weight: u64,
    leaf_count: usize,
}

impl GearBacklashModeShareV1 {
    /// Structured law identity.
    #[must_use]
    pub const fn law(&self) -> StructuredLawId {
        self.law
    }

    /// Stable structured law key.
    #[must_use]
    pub fn law_key(&self) -> &str {
        &self.law_key
    }

    /// Law-local response-piece ordinal.
    #[must_use]
    pub const fn piece(&self) -> usize {
        self.piece
    }

    /// Stable response-mode key.
    #[must_use]
    pub fn mode_key(&self) -> &str {
        &self.mode_key
    }

    /// Exact multiplicity selecting this mode.
    #[must_use]
    pub const fn relative_weight(&self) -> u64 {
        self.relative_weight
    }

    /// Distinct retained leaves selecting this mode.
    #[must_use]
    pub const fn leaf_count(&self) -> usize {
        self.leaf_count
    }
}

/// Caller declaration for one deterministic backlash-consumption report.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GearBacklashConsumerDraftV1 {
    /// Explicit unit assigned by the caller to all structured output values.
    pub output_unit: GearBacklashLengthUnitV1,
    /// Requested probabilities in non-semantic caller order.
    pub quantiles: Vec<GearBacklashProbabilityV1>,
}

impl GearBacklashConsumerDraftV1 {
    /// Consume one admitted structured receipt without reevaluating its model.
    ///
    /// # Errors
    /// Refuses resource/alias, retained-invariant, or SI-conversion gaps.
    pub fn consume(
        self,
        structured: &StructuredPropagationReceipt,
    ) -> Result<GearBacklashConsumptionReceiptV1, GearBacklashConsumptionErrorV1> {
        consume_structured_backlash(self, structured)
    }
}

/// Deterministic weighted-empirical backlash report.
#[derive(Debug, Clone, PartialEq)]
pub struct GearBacklashConsumptionReceiptV1 {
    schema_version: u32,
    structured_schema_version: u32,
    model: StructuredPopulationModel,
    output_unit: GearBacklashLengthUnitV1,
    total_weight: u64,
    mean: GearBacklashValueV1,
    variance: GearBacklashVarianceV1,
    standard_deviation: GearBacklashValueV1,
    support: Box<[GearBacklashSupportPointV1]>,
    quantiles: Box<[GearBacklashQuantileV1]>,
    modes: Box<[GearBacklashModeShareV1]>,
}

impl GearBacklashConsumptionReceiptV1 {
    /// Fixed consumer schema version.
    #[must_use]
    pub const fn schema_version(&self) -> u32 {
        self.schema_version
    }

    /// Exact upstream structured-propagation schema consumed by this report.
    #[must_use]
    pub const fn structured_schema_version(&self) -> u32 {
        self.structured_schema_version
    }

    /// Complete admitted model retained for replay and semantic comparison.
    #[must_use]
    pub const fn model(&self) -> &StructuredPopulationModel {
        &self.model
    }

    /// Exact external identity retained by the structured model owner.
    #[must_use]
    pub const fn model_identity(&self) -> &StructuredModelIdentity {
        &self.model.identity
    }

    /// Caller-declared interpretation unit.
    #[must_use]
    pub const fn output_unit(&self) -> GearBacklashLengthUnitV1 {
        self.output_unit
    }

    /// Exact finite-population multiplicity.
    #[must_use]
    pub const fn total_weight(&self) -> u64 {
        self.total_weight
    }

    /// Structured population mean with explicit unit conversion.
    #[must_use]
    pub const fn mean(&self) -> GearBacklashValueV1 {
        self.mean
    }

    /// Structured population variance with squared-unit conversion.
    #[must_use]
    pub const fn variance(&self) -> GearBacklashVarianceV1 {
        self.variance
    }

    /// Structured population standard deviation with explicit conversion.
    #[must_use]
    pub const fn standard_deviation(&self) -> GearBacklashValueV1 {
        self.standard_deviation
    }

    /// Complete ascending weighted support used by every reported quantile.
    #[must_use]
    pub fn support(&self) -> &[GearBacklashSupportPointV1] {
        &self.support
    }

    /// Quantiles in ascending exact-probability order.
    #[must_use]
    pub fn quantiles(&self) -> &[GearBacklashQuantileV1] {
        &self.quantiles
    }

    /// Law-major, piece-minor mode shares retained from the structured receipt.
    #[must_use]
    pub fn modes(&self) -> &[GearBacklashModeShareV1] {
        &self.modes
    }
}

/// Numeric quantity whose unit conversion failed closed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GearBacklashNumericQuantityV1 {
    /// Population mean.
    Mean,
    /// Population variance.
    Variance,
    /// Population standard deviation.
    StandardDeviation,
    /// One row of the complete ascending support table.
    Support(usize),
}

/// Stable unit-conversion failure class.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GearBacklashNumericIssueV1 {
    /// Source or converted value was NaN or infinite.
    NonFinite,
    /// A variance or standard deviation was negative.
    Negative,
    /// A nonzero source value vanished during SI conversion.
    SiUnderflow,
}

/// Structured refusal from backlash consumption.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GearBacklashConsumptionErrorV1 {
    /// At least one quantile must be requested.
    NoQuantiles,
    /// Raw quantile requests exceeded the fixed cap.
    QuantileLimit {
        /// Submitted request count.
        actual: usize,
        /// Maximum admitted count.
        max: usize,
    },
    /// One normalized probability appeared more than once.
    DuplicateQuantile {
        /// Repeated exact probability.
        probability: GearBacklashProbabilityV1,
    },
    /// The upstream evaluator schema was not the one this consumer implements.
    StructuredSchema {
        /// Supported upstream schema.
        expected: u32,
        /// Supplied upstream schema.
        actual: u32,
    },
    /// Privately constructed structured support was unexpectedly empty.
    ReceiptSupportEmpty,
    /// Retained multiplicity accumulation overflowed.
    ReceiptWeightOverflow,
    /// Leaf weights disagreed with the structured receipt total.
    ReceiptWeightMismatch {
        /// Total declared by the structured receipt.
        declared: u64,
        /// Total reconstructed from retained leaves.
        reconstructed: u64,
    },
    /// One retained mode could not resolve its law/piece key.
    ReceiptModeInvariantGap {
        /// Missing structured law.
        law: StructuredLawId,
        /// Missing law-local piece ordinal.
        piece: usize,
    },
    /// One retained mode disagreed with aggregation over the retained leaves.
    ReceiptModeLeafMismatch {
        /// Structured law identity.
        law: StructuredLawId,
        /// Law-local piece ordinal.
        piece: usize,
        /// Mode receipt's exact multiplicity.
        declared_weight: u64,
        /// Multiplicity independently reconstructed from leaves.
        reconstructed_weight: u64,
        /// Mode receipt's distinct-leaf count.
        declared_leaf_count: usize,
        /// Distinct-leaf count independently reconstructed from leaves.
        reconstructed_leaf_count: usize,
    },
    /// Mode weights disagreed with the structured receipt total.
    ReceiptModeWeightMismatch {
        /// Total declared by the structured receipt.
        declared: u64,
        /// Total reconstructed from retained modes.
        reconstructed: u64,
    },
    /// Two distinct source support values collapsed to one coherent-SI value.
    SiSupportAliasing {
        /// Lower source support bits.
        lower_source_bits: u64,
        /// Upper source support bits.
        upper_source_bits: u64,
        /// Colliding coherent-SI bits.
        metres_bits: u64,
    },
    /// A retained source or coherent-SI value could not be represented.
    Numeric {
        /// Quantity being converted.
        quantity: GearBacklashNumericQuantityV1,
        /// Exact failure class.
        issue: GearBacklashNumericIssueV1,
    },
}

impl GearBacklashConsumptionErrorV1 {
    /// Stable machine-actionable diagnostic code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::NoQuantiles => "GearBacklashNoQuantiles",
            Self::QuantileLimit { .. } => "GearBacklashQuantileLimit",
            Self::DuplicateQuantile { .. } => "GearBacklashDuplicateQuantile",
            Self::StructuredSchema { .. } => "GearBacklashStructuredSchema",
            Self::ReceiptSupportEmpty => "GearBacklashReceiptSupportEmpty",
            Self::ReceiptWeightOverflow => "GearBacklashReceiptWeightOverflow",
            Self::ReceiptWeightMismatch { .. } => "GearBacklashReceiptWeightMismatch",
            Self::ReceiptModeInvariantGap { .. } => "GearBacklashReceiptModeInvariantGap",
            Self::ReceiptModeLeafMismatch { .. } => "GearBacklashReceiptModeLeafMismatch",
            Self::ReceiptModeWeightMismatch { .. } => "GearBacklashReceiptModeWeightMismatch",
            Self::SiSupportAliasing { .. } => "GearBacklashSiSupportAliasing",
            Self::Numeric { .. } => "GearBacklashNumeric",
        }
    }
}

impl fmt::Display for GearBacklashConsumptionErrorV1 {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NoQuantiles => formatter.write_str("gear-backlash report requires a quantile"),
            Self::QuantileLimit { actual, max } => {
                write!(
                    formatter,
                    "gear-backlash report has {actual} quantiles; maximum is {max}"
                )
            }
            Self::DuplicateQuantile { probability } => write!(
                formatter,
                "gear-backlash probability {}/{} is repeated",
                probability.numerator(),
                probability.denominator()
            ),
            Self::StructuredSchema { expected, actual } => write!(
                formatter,
                "gear-backlash consumer expects structured schema {expected}, got {actual}"
            ),
            Self::ReceiptSupportEmpty => {
                formatter.write_str("structured receipt retained no backlash support")
            }
            Self::ReceiptWeightOverflow => {
                formatter.write_str("structured backlash weight accumulation overflowed")
            }
            Self::ReceiptWeightMismatch {
                declared,
                reconstructed,
            } => write!(
                formatter,
                "structured backlash leaf weight {reconstructed} differs from declared {declared}"
            ),
            Self::ReceiptModeInvariantGap { law, piece } => write!(
                formatter,
                "structured backlash mode law {} piece {piece} is unresolved",
                law.0
            ),
            Self::ReceiptModeLeafMismatch {
                law,
                piece,
                declared_weight,
                reconstructed_weight,
                declared_leaf_count,
                reconstructed_leaf_count,
            } => write!(
                formatter,
                "structured backlash mode law {} piece {piece} declares weight/count \
                 {declared_weight}/{declared_leaf_count}, reconstructed \
                 {reconstructed_weight}/{reconstructed_leaf_count}",
                law.0
            ),
            Self::ReceiptModeWeightMismatch {
                declared,
                reconstructed,
            } => write!(
                formatter,
                "structured backlash mode weight {reconstructed} differs from declared {declared}"
            ),
            Self::SiSupportAliasing {
                lower_source_bits,
                upper_source_bits,
                metres_bits,
            } => write!(
                formatter,
                "distinct backlash support bits {lower_source_bits:#018x} and \
                 {upper_source_bits:#018x} collapse to SI bits {metres_bits:#018x}"
            ),
            Self::Numeric { quantity, issue } => {
                write!(
                    formatter,
                    "gear-backlash {quantity:?} conversion failed: {issue:?}"
                )
            }
        }
    }
}

impl std::error::Error for GearBacklashConsumptionErrorV1 {}

#[derive(Debug, Clone, Copy)]
struct SupportPoint {
    value: f64,
    relative_weight: u64,
}

/// Consume an admitted structured finite population as weighted backlash data.
///
/// The quantile convention is `inf { x | F_weighted(x) >= p }`, with `p = 0`
/// explicitly selecting the minimum support value. Equal binary64 support
/// values are grouped before cumulative brackets are reported.
///
/// # Errors
/// Refuses resource/alias, retained-invariant, or SI-conversion gaps.
#[allow(clippy::too_many_lines)] // One atomic audit: support, exact ranks, modes, then publication.
pub fn consume_structured_backlash(
    draft: GearBacklashConsumerDraftV1,
    structured: &StructuredPropagationReceipt,
) -> Result<GearBacklashConsumptionReceiptV1, GearBacklashConsumptionErrorV1> {
    if draft.quantiles.is_empty() {
        return Err(GearBacklashConsumptionErrorV1::NoQuantiles);
    }
    if draft.quantiles.len() > MAX_GEAR_BACKLASH_QUANTILES_V1 {
        return Err(GearBacklashConsumptionErrorV1::QuantileLimit {
            actual: draft.quantiles.len(),
            max: MAX_GEAR_BACKLASH_QUANTILES_V1,
        });
    }
    let mut probabilities = draft.quantiles;
    probabilities.sort_unstable();
    if let Some(pair) = probabilities.windows(2).find(|pair| pair[0] == pair[1]) {
        return Err(GearBacklashConsumptionErrorV1::DuplicateQuantile {
            probability: pair[0],
        });
    }
    if structured.schema_version() != STRUCTURED_PROPAGATION_SCHEMA_V1 {
        return Err(GearBacklashConsumptionErrorV1::StructuredSchema {
            expected: STRUCTURED_PROPAGATION_SCHEMA_V1,
            actual: structured.schema_version(),
        });
    }

    let mut leaves = structured
        .leaves()
        .iter()
        .map(|leaf| (leaf.output(), leaf.relative_weight().get()))
        .collect::<Vec<_>>();
    leaves.sort_by(|left, right| left.0.total_cmp(&right.0));
    let mut support = Vec::<SupportPoint>::with_capacity(leaves.len());
    for (value, relative_weight) in leaves {
        if let Some(last) = support.last_mut()
            && last.value.total_cmp(&value) == Ordering::Equal
        {
            last.relative_weight = last
                .relative_weight
                .checked_add(relative_weight)
                .ok_or(GearBacklashConsumptionErrorV1::ReceiptWeightOverflow)?;
        } else {
            support.push(SupportPoint {
                value,
                relative_weight,
            });
        }
    }
    if support.is_empty() {
        return Err(GearBacklashConsumptionErrorV1::ReceiptSupportEmpty);
    }
    let reconstructed_weight = support.iter().try_fold(0_u64, |sum, point| {
        sum.checked_add(point.relative_weight)
            .ok_or(GearBacklashConsumptionErrorV1::ReceiptWeightOverflow)
    })?;
    let total_weight = structured.total_weight();
    if reconstructed_weight != total_weight {
        return Err(GearBacklashConsumptionErrorV1::ReceiptWeightMismatch {
            declared: total_weight,
            reconstructed: reconstructed_weight,
        });
    }

    let mut published_support = Vec::<GearBacklashSupportPointV1>::with_capacity(support.len());
    let mut cumulative = 0_u64;
    for (index, point) in support.into_iter().enumerate() {
        let cumulative_before = cumulative;
        cumulative = cumulative
            .checked_add(point.relative_weight)
            .ok_or(GearBacklashConsumptionErrorV1::ReceiptWeightOverflow)?;
        let value = convert_value(
            point.value,
            draft.output_unit,
            GearBacklashNumericQuantityV1::Support(index),
            false,
        )?;
        if let Some(previous) = published_support.last()
            && previous.value().metres().total_cmp(&value.metres()) == Ordering::Equal
        {
            return Err(GearBacklashConsumptionErrorV1::SiSupportAliasing {
                lower_source_bits: previous.value().source_value().to_bits(),
                upper_source_bits: value.source_value().to_bits(),
                metres_bits: value.metres().to_bits(),
            });
        }
        published_support.push(GearBacklashSupportPointV1 {
            value,
            relative_weight: point.relative_weight,
            cumulative_before,
            cumulative_at: cumulative,
        });
    }

    let mut quantiles = Vec::with_capacity(probabilities.len());
    for probability in probabilities {
        let selected = published_support
            .iter()
            .enumerate()
            .find(|(_, point)| probability.is_reached_by(point.cumulative_at(), total_weight));
        let Some((support_index, point)) = selected else {
            return Err(GearBacklashConsumptionErrorV1::ReceiptWeightMismatch {
                declared: total_weight,
                reconstructed: published_support
                    .last()
                    .map_or(0, |point| point.cumulative_at()),
            });
        };
        quantiles.push(GearBacklashQuantileV1 {
            probability,
            support_index,
            value: point.value(),
            cumulative_before: point.cumulative_before(),
            cumulative_at: point.cumulative_at(),
        });
    }

    let mut leaf_modes = BTreeMap::<(StructuredLawId, usize), (u64, usize)>::new();
    for leaf in structured.leaves() {
        let aggregate = leaf_modes
            .entry((leaf.law(), leaf.selected_piece()))
            .or_default();
        aggregate.0 = aggregate
            .0
            .checked_add(leaf.relative_weight().get())
            .ok_or(GearBacklashConsumptionErrorV1::ReceiptWeightOverflow)?;
        aggregate.1 = aggregate
            .1
            .checked_add(1)
            .ok_or(GearBacklashConsumptionErrorV1::ReceiptWeightOverflow)?;
    }
    let mut modes = Vec::with_capacity(structured.modes().len());
    let mut mode_weight = 0_u64;
    for mode in structured.modes() {
        let Some(law) = structured.model().laws.get(mode.law().index()) else {
            return Err(GearBacklashConsumptionErrorV1::ReceiptModeInvariantGap {
                law: mode.law(),
                piece: mode.piece(),
            });
        };
        let Some(piece) = law.pieces.get(mode.piece()) else {
            return Err(GearBacklashConsumptionErrorV1::ReceiptModeInvariantGap {
                law: mode.law(),
                piece: mode.piece(),
            });
        };
        let (reconstructed_weight, reconstructed_leaf_count) = leaf_modes
            .remove(&(mode.law(), mode.piece()))
            .unwrap_or_default();
        if mode.relative_weight() != reconstructed_weight
            || mode.leaf_count() != reconstructed_leaf_count
        {
            return Err(GearBacklashConsumptionErrorV1::ReceiptModeLeafMismatch {
                law: mode.law(),
                piece: mode.piece(),
                declared_weight: mode.relative_weight(),
                reconstructed_weight,
                declared_leaf_count: mode.leaf_count(),
                reconstructed_leaf_count,
            });
        }
        mode_weight = mode_weight
            .checked_add(mode.relative_weight())
            .ok_or(GearBacklashConsumptionErrorV1::ReceiptWeightOverflow)?;
        modes.push(GearBacklashModeShareV1 {
            law: mode.law(),
            law_key: law.key.clone().into_boxed_str(),
            piece: mode.piece(),
            mode_key: piece.mode_key.clone().into_boxed_str(),
            relative_weight: mode.relative_weight(),
            leaf_count: mode.leaf_count(),
        });
    }
    if mode_weight != total_weight {
        return Err(GearBacklashConsumptionErrorV1::ReceiptModeWeightMismatch {
            declared: total_weight,
            reconstructed: mode_weight,
        });
    }
    if let Some(((law, piece), (reconstructed_weight, reconstructed_leaf_count))) =
        leaf_modes.into_iter().next()
    {
        return Err(GearBacklashConsumptionErrorV1::ReceiptModeLeafMismatch {
            law,
            piece,
            declared_weight: 0,
            reconstructed_weight,
            declared_leaf_count: 0,
            reconstructed_leaf_count,
        });
    }

    Ok(GearBacklashConsumptionReceiptV1 {
        schema_version: GEAR_BACKLASH_CONSUMER_SCHEMA_V1,
        structured_schema_version: structured.schema_version(),
        model: structured.model().clone(),
        output_unit: draft.output_unit,
        total_weight,
        mean: convert_value(
            structured.mean(),
            draft.output_unit,
            GearBacklashNumericQuantityV1::Mean,
            false,
        )?,
        variance: convert_variance(structured.variance(), draft.output_unit)?,
        standard_deviation: convert_value(
            structured.standard_deviation(),
            draft.output_unit,
            GearBacklashNumericQuantityV1::StandardDeviation,
            true,
        )?,
        support: published_support.into_boxed_slice(),
        quantiles: quantiles.into_boxed_slice(),
        modes: modes.into_boxed_slice(),
    })
}

fn convert_value(
    source: f64,
    unit: GearBacklashLengthUnitV1,
    quantity: GearBacklashNumericQuantityV1,
    require_nonnegative: bool,
) -> Result<GearBacklashValueV1, GearBacklashConsumptionErrorV1> {
    if !source.is_finite() {
        return Err(GearBacklashConsumptionErrorV1::Numeric {
            quantity,
            issue: GearBacklashNumericIssueV1::NonFinite,
        });
    }
    if require_nonnegative && source < 0.0 {
        return Err(GearBacklashConsumptionErrorV1::Numeric {
            quantity,
            issue: GearBacklashNumericIssueV1::Negative,
        });
    }
    let metres = source * unit.metres_per_unit();
    if !metres.is_finite() {
        return Err(GearBacklashConsumptionErrorV1::Numeric {
            quantity,
            issue: GearBacklashNumericIssueV1::NonFinite,
        });
    }
    if source != 0.0 && metres == 0.0 {
        return Err(GearBacklashConsumptionErrorV1::Numeric {
            quantity,
            issue: GearBacklashNumericIssueV1::SiUnderflow,
        });
    }
    Ok(GearBacklashValueV1 {
        source_bits: source.to_bits(),
        unit,
        metres_bits: metres.to_bits(),
    })
}

fn convert_variance(
    source: f64,
    unit: GearBacklashLengthUnitV1,
) -> Result<GearBacklashVarianceV1, GearBacklashConsumptionErrorV1> {
    let quantity = GearBacklashNumericQuantityV1::Variance;
    if !source.is_finite() {
        return Err(GearBacklashConsumptionErrorV1::Numeric {
            quantity,
            issue: GearBacklashNumericIssueV1::NonFinite,
        });
    }
    if source < 0.0 {
        return Err(GearBacklashConsumptionErrorV1::Numeric {
            quantity,
            issue: GearBacklashNumericIssueV1::Negative,
        });
    }
    let scale = unit.metres_per_unit();
    let square_metres = source * scale * scale;
    if !square_metres.is_finite() {
        return Err(GearBacklashConsumptionErrorV1::Numeric {
            quantity,
            issue: GearBacklashNumericIssueV1::NonFinite,
        });
    }
    if source != 0.0 && square_metres == 0.0 {
        return Err(GearBacklashConsumptionErrorV1::Numeric {
            quantity,
            issue: GearBacklashNumericIssueV1::SiUnderflow,
        });
    }
    Ok(GearBacklashVarianceV1 {
        source_bits: source.to_bits(),
        unit,
        square_metres_bits: square_metres.to_bits(),
    })
}

const fn greatest_common_divisor(mut left: u64, mut right: u64) -> u64 {
    while right != 0 {
        let remainder = left % right;
        left = right;
        right = remainder;
    }
    if left == 0 { 1 } else { left }
}
