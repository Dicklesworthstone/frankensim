//! Source-bound particle classification and conservative two-product streams.
//!
//! The Rosin–Rammler corrected partition is
//! `c(d) = 1 - exp(-ln(2) * (d/d50)^m)`. Actual solids recovery to underflow
//! is `b + (1-b)c(d)` under the explicit water-bypass approximation.
//! Source: Kawatra et al. (2005), DOE report DE-FC26-01NT41062, equation 4,
//! <https://doi.org/10.2172/887498>; Kujawa (2011), *Cycloning of Tailings for
//! the Production of Sand as TSF Construction Material*, equation 1,
//! <https://www.patersoncooke.com/wp-content/uploads/TMW2011-Christian-Kujawa-881.pdf>.
//! `ln(2)` retains the exact mathematical 50% definition of corrected d50.
//!
//! Cut size, sharpness and bypass are supplied test parameters, not predictions
//! from geometry. There is no inferred fishhook, density correction, liberation,
//! breakage, dissolution or chemical reaction. A class is compositionally
//! homogeneous for the declared split; different mineral responses need separate
//! classes and explicitly supplied partitions. All flows are coherent SI.

use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, MatDbError, QueryPoint};
use fs_math::det;
use fs_qty::{Dimensionless, Dims, Length, QuantitySpec};

#[path = "separation/stream.rs"]
mod stream;
pub use stream::{ParticleClass, ParticleStream, SeparatedStreams};

/// Version-1 stateless empirical classifier identity in an ordinary model card.
pub const CLASSIFIER_LAW_ID: &str = "rosin-rammler-water-bypass-classifier";
/// Owned query axis, in metres; must have finite positive declared bounds.
pub const PARTICLE_DIAMETER_AXIS: &str = "particle-diameter";
/// Owned query axis, dry solids / (dry solids + water) by mass.
pub const SOLIDS_FRACTION_AXIS: &str = "solids-mass-fraction";
const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);

/// Invalid source, unsupported operating point or unrepresentable material flow.
#[derive(Clone, Debug, PartialEq)]
pub enum SeparationError {
    /// Model-card or typed coordinate admission failed in fs-matdb.
    Source(MatDbError),
    /// Invalid physical input or unsupported model schema.
    Invalid(&'static str),
    /// Missing, mistyped or out-of-range source operating conditions.
    OutsideValidity,
    /// A nonzero flow or finite sum cannot be represented.
    Arithmetic(&'static str),
}
impl core::fmt::Display for SeparationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "particle separation: {self:?}")
    }
}
impl std::error::Error for SeparationError {}
impl From<MatDbError> for SeparationError {
    fn from(error: MatDbError) -> Self {
        Self::Source(error)
    }
}

/// Two complementary fractions with both small tails retained explicitly.
/// No normalization or inferred material recovery is performed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Partition {
    product: Dimensionless,
    reject: Dimensionless,
}

impl Partition {
    /// Supply the product fraction and calculate its complement.
    pub fn from_product(product: Dimensionless) -> Result<Self, SeparationError> {
        fraction(product.value())?;
        Self::from_fractions(product, Dimensionless::new(1.0 - product.value()))
    }

    /// Retain independently evaluated complementary fractions, e.g. a tiny
    /// exponential tail alongside a product fraction rounded to one. Their sum
    /// must equal one within eight binary64 epsilons; neither is rescaled.
    pub fn from_fractions(
        product: Dimensionless,
        reject: Dimensionless,
    ) -> Result<Self, SeparationError> {
        fraction(product.value())?;
        fraction(reject.value())?;
        require(
            (product.value() + reject.value() - 1.0).abs() <= 8.0 * f64::EPSILON,
            "partition fractions must sum to one",
        )?;
        Ok(Self { product, reject })
    }

    /// Fraction entering the named product or classifier underflow.
    pub const fn product(self) -> Dimensionless {
        self.product
    }
    /// Fraction entering the reject or classifier overflow, retaining small tails.
    pub const fn reject(self) -> Dimensionless {
        self.reject
    }
}

/// Immutable executable classifier retaining its complete source model card.
#[derive(Clone, Debug, PartialEq)]
pub struct ClassifierCard {
    source: ConstitutiveModelCard,
    cut_size: Length,
    sharpness: f64,
    water: Partition,
}

impl ClassifierCard {
    /// Admit the exact law/version with three dimensioned parameters:
    /// `corrected-cut-size` [m], `sharpness` [1], and
    /// `water-recovery-to-underflow` [1]. Both owned query axes require declared
    /// finite bounds; all other source conditions are supplied at evaluation.
    pub fn from_model_card(source: &ConstitutiveModelCard) -> Result<Self, SeparationError> {
        source.validate()?;
        require(
            source.law.0 == CLASSIFIER_LAW_ID && source.law_version == 1,
            "unsupported classifier law or version",
        )?;
        require(
            source.state_schema_version == 1
                && source.initial_state == InitialStatePolicy::ZeroInternalState,
            "classifier requires stateless schema version one",
        )?;
        require(
            !source.sources.is_empty(),
            "classifier requires source artifact identities",
        )?;
        require(
            source.parameters.len() == 3,
            "classifier requires exactly three parameters",
        )?;
        let parameter = |name: &str, dims: Dims| -> Result<f64, SeparationError> {
            let value = source
                .parameters
                .get(name)
                .ok_or(SeparationError::Invalid("missing classifier parameter"))?;
            require(
                value.dims == dims,
                "classifier parameter has wrong dimensions",
            )?;
            Ok(value.value)
        };
        let cut_size = positive(
            parameter("corrected-cut-size", LENGTH)?,
            "cut size must be positive",
        )?;
        let sharpness = positive(
            parameter("sharpness", Dims::NONE)?,
            "sharpness must be positive",
        )?;
        let water = Partition::from_product(Dimensionless::new(parameter(
            "water-recovery-to-underflow",
            Dims::NONE,
        )?))?;
        for (axis, dims) in [
            (PARTICLE_DIAMETER_AXIS, LENGTH),
            (SOLIDS_FRACTION_AXIS, Dims::NONE),
        ] {
            let (lo, hi) = source
                .validity
                .bound(axis)
                .ok_or(SeparationError::Invalid("missing classifier validity axis"))?;
            require(
                lo.is_finite() && hi.is_finite() && lo <= hi,
                "unusable classifier validity bounds",
            )?;
            if axis == PARTICLE_DIAMETER_AXIS {
                positive(lo, "size support must be positive")?;
            } else {
                fraction(lo)?;
                fraction(hi)?;
            }
            if let Some(quantity) = source.validity.axis_quantities().get(axis) {
                require(
                    *quantity == QuantitySpec::dimensional(dims),
                    "classifier axis has wrong quantity convention",
                )?;
            }
        }
        Ok(Self {
            source: source.clone(),
            cut_size: Length::new(cut_size),
            sharpness,
            water,
        })
    }

    /// Original immutable parameters, source identities and applicability domain.
    pub const fn model_card(&self) -> &ConstitutiveModelCard {
        &self.source
    }
    /// Declared water split, also used as the fine-particle bypass approximation.
    pub const fn water_partition(&self) -> Partition {
        self.water
    }

    /// Evaluate actual solids partition to underflow/overflow at one class size.
    /// Extra source conditions (e.g. pressure or test identity) must be supplied
    /// explicitly. The two owned axes cannot be overridden through `conditions`.
    pub fn partition(
        &self,
        diameter: Length,
        solids_fraction: Dimensionless,
        conditions: &QueryPoint,
    ) -> Result<Partition, SeparationError> {
        positive(
            diameter.value(),
            "particle size must be finite and positive",
        )?;
        fraction(solids_fraction.value())?;
        require(
            !conditions.axes().contains_key(PARTICLE_DIAMETER_AXIS)
                && !conditions.axes().contains_key(SOLIDS_FRACTION_AXIS),
            "caller conditions override owned stream axes",
        )?;
        let mut point = conditions.clone();
        for (axis, value) in [
            (PARTICLE_DIAMETER_AXIS, diameter.value()),
            (SOLIDS_FRACTION_AXIS, solids_fraction.value()),
        ] {
            point = if let Some(quantity) = self.source.validity.axis_quantities().get(axis) {
                point.with_quantity(axis, *quantity, value)?
            } else {
                point.with(axis, value)?
            };
        }
        if !self
            .source
            .validity
            .contains_typed(point.axes(), point.axis_quantities())
        {
            return Err(SeparationError::OutsideValidity);
        }
        // Logarithms avoid an overflowing diameter ratio. The exponential tails
        // are evaluated separately so fine-product and coarse-reject flows are
        // not erased by subtraction from a rounded unit fraction.
        let relative_size = (diameter.value() - self.cut_size.value()) / self.cut_size.value();
        let log_ratio = if relative_size.abs() <= 0.5 {
            // Recover rounding in 1+x before the logarithm. Subtracting two
            // separately rounded logarithms loses adjacent sizes near d50.
            let rounded = 1.0 + relative_size;
            if rounded == 1.0 {
                relative_size
            } else {
                det::ln(rounded) * (relative_size / (rounded - 1.0))
            }
        } else {
            det::ln(diameter.value()) - det::ln(self.cut_size.value())
        };
        let log_power = self.sharpness * log_ratio;
        let exponent = core::f64::consts::LN_2 * det::exp(log_power);
        let corrected = -det::expm1(-exponent);
        let underflow = self.water.product.value() + self.water.reject.value() * corrected;
        let overflow = self.water.reject.value() * det::exp(-exponent);
        Partition::from_fractions(Dimensionless::new(underflow), Dimensionless::new(overflow))
    }

    /// Classify all solid classes and water using the stream's ACTUAL solids
    /// mass fraction. A single supplied card applies to every class here;
    /// component-specific laws use `ParticleStream::split` with explicit splits.
    pub fn classify(
        &self,
        feed: &ParticleStream,
        conditions: &QueryPoint,
    ) -> Result<SeparatedStreams, SeparationError> {
        let solids_fraction = feed.solids_mass_fraction()?;
        let partitions = feed
            .classes()
            .iter()
            .map(|class| {
                Ok((
                    class.id().to_owned(),
                    self.partition(class.diameter(), solids_fraction, conditions)?,
                ))
            })
            .collect::<Result<Vec<_>, SeparationError>>()?;
        feed.split(&partitions, self.water)
    }
}

fn require(condition: bool, reason: &'static str) -> Result<(), SeparationError> {
    if condition {
        Ok(())
    } else {
        Err(SeparationError::Invalid(reason))
    }
}
fn positive(value: f64, reason: &'static str) -> Result<f64, SeparationError> {
    require(value.is_finite() && value > 0.0, reason)?;
    Ok(value)
}
fn nonnegative(value: f64, reason: &'static str) -> Result<f64, SeparationError> {
    require(value.is_finite() && value >= 0.0, reason)?;
    Ok(value)
}
fn fraction(value: f64) -> Result<f64, SeparationError> {
    require(
        value.is_finite() && (0.0..=1.0).contains(&value),
        "fraction must be finite and in [0,1]",
    )?;
    Ok(value)
}

#[cfg(test)]
#[path = "separation/tests.rs"]
mod tests;
