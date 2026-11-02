//! Selective, prescribed partitions of existing solid constituents and carrier.
//!
//! A caller supplies the fraction of each exact component identity reporting to
//! concentrate, plus a separate carrier-liquid fraction. The complementary flow
//! reports to tailings. These are declared separation inputs, not inferred
//! flotation recovery, cyclone efficiency, liberation, kinetics or calibration.
//! Identities may describe minerals, size classes, or mineral-by-size classes;
//! reported dry fractions refer only to those identities and are not elemental
//! assays. Source density evidence and the additive-volume stream assumption
//! remain unchanged. No reaction, dissolution or new constituent is introduced.

use super::{
    MAX_COMPONENT_ID_BYTES, MAX_SOLID_COMPONENTS, SlurryComponent, SlurryError, SlurryStream,
};
use fs_matdb::Provenance;
use fs_qty::{Dimensionless, MassFlowRate};

/// Invalid partitions, incomplete component coverage, or a refused outlet stream.
#[derive(Clone, Debug, PartialEq)]
pub enum SeparationError {
    /// Original component/stream admission or finite-arithmetic refusal.
    Slurry(SlurryError),
    /// Invalid physical input or unrepresentable separation metric.
    Invalid(&'static str),
    /// The prescribed solid partitions do not exactly cover the feed identities.
    PartitionCoverage {
        /// Feed identities without a prescribed partition.
        missing: Vec<String>,
        /// Prescribed identities absent from the feed's solid set.
        unknown: Vec<String>,
    },
}

impl core::fmt::Display for SeparationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Slurry(error) => write!(f, "selective separation: {error}"),
            Self::Invalid(what) => write!(f, "selective separation: {what}"),
            Self::PartitionCoverage { missing, unknown } => write!(
                f,
                "selective separation: missing partitions {missing:?}; unknown partitions {unknown:?}"
            ),
        }
    }
}
impl std::error::Error for SeparationError {}
impl From<SlurryError> for SeparationError {
    fn from(value: SlurryError) -> Self {
        Self::Slurry(value)
    }
}

type Result<T> = core::result::Result<T, SeparationError>;

fn require(condition: bool, message: &'static str) -> Result<()> {
    if condition {
        Ok(())
    } else {
        Err(SeparationError::Invalid(message))
    }
}

fn fraction(value: Dimensionless) -> Result<Dimensionless> {
    require(
        value.value().is_finite() && (0.0..=1.0).contains(&value.value()),
        "partition fraction must be finite and in [0, 1]",
    )?;
    Ok(Dimensionless::new(value.value().abs()))
}

/// One exact solid identity's prescribed fraction to the concentrate outlet.
#[derive(Clone, Debug, PartialEq)]
pub struct ComponentPartition {
    component: String,
    to_concentrate: Dimensionless,
}

impl ComponentPartition {
    /// Declare a fraction in [0, 1], without clipping or a default for other solids.
    pub fn new(component: impl Into<String>, to_concentrate: Dimensionless) -> Result<Self> {
        let component = component.into();
        require(
            !component.trim().is_empty() && component.len() <= MAX_COMPONENT_ID_BYTES,
            "component identity must contain 1 to 256 bytes and cannot be blank",
        )?;
        Ok(Self {
            component,
            to_concentrate: fraction(to_concentrate)?,
        })
    }

    /// Exact feed-solid identity; no mineral-name or size-class inference.
    #[must_use]
    pub fn component(&self) -> &str {
        &self.component
    }

    /// Prescribed fraction of this component's feed mass reporting to concentrate.
    #[must_use]
    pub const fn to_concentrate(&self) -> Dimensionless {
        self.to_concentrate
    }
}

/// Bounded caller-prescribed separation with retained source/assumption provenance.
#[derive(Clone, Debug, PartialEq)]
pub struct SelectiveSeparator {
    partitions: Vec<ComponentPartition>,
    carrier_to_concentrate: Dimensionless,
    provenance: Provenance,
}

impl SelectiveSeparator {
    /// Admit up to 128 distinct solid partitions and one carrier-liquid partition.
    ///
    /// Provenance describes these prescribed fractions. Supplying a citation or
    /// artifact alone does not establish a calibrated cyclone/flotation model.
    pub fn new(
        mut partitions: Vec<ComponentPartition>,
        carrier_to_concentrate: Dimensionless,
        provenance: Provenance,
    ) -> Result<Self> {
        require(
            partitions.len() <= MAX_SOLID_COMPONENTS,
            "too many component partitions",
        )?;
        partitions.sort_by(|left, right| left.component.cmp(&right.component));
        require(
            partitions
                .windows(2)
                .all(|pair| pair[0].component != pair[1].component),
            "component partitions must have distinct identities",
        )?;
        require(
            !provenance.source.trim().is_empty() && !provenance.license.trim().is_empty(),
            "prescribed partitions require nonblank source and license provenance",
        )?;
        Ok(Self {
            partitions,
            carrier_to_concentrate: fraction(carrier_to_concentrate)?,
            provenance,
        })
    }

    /// Canonically ordered prescribed partitions for the exact solid identity set.
    #[must_use]
    pub fn partitions(&self) -> &[ComponentPartition] {
        &self.partitions
    }

    /// Prescribed carrier-liquid fraction, independent of all solid partitions.
    #[must_use]
    pub const fn carrier_to_concentrate(&self) -> Dimensionless {
        self.carrier_to_concentrate
    }

    /// Original provenance of the caller's separation assumptions or source data.
    #[must_use]
    pub fn provenance(&self) -> &Provenance {
        &self.provenance
    }

    /// Produce concentrate and tailings without mutating the feed.
    ///
    /// Every feed-solid identity, including zero-flow identities, must have
    /// exactly one partition; extra identities refuse. A zero-flow component has
    /// no computed recovery. Dry-component grades are undefined in an outlet
    /// with no dry solids. Positive flows erased by underflow refuse atomically.
    pub fn separate(&self, feed: &SlurryStream) -> Result<SeparationResult> {
        let missing: Vec<_> = feed
            .solids()
            .iter()
            .filter(|solid| {
                self.partitions
                    .binary_search_by(|partition| partition.component.as_str().cmp(solid.id()))
                    .is_err()
            })
            .map(|solid| solid.id().to_owned())
            .collect();
        let unknown: Vec<_> = self
            .partitions
            .iter()
            .filter(|partition| {
                feed.solids()
                    .binary_search_by(|solid| solid.id().cmp(partition.component.as_str()))
                    .is_err()
            })
            .map(|partition| partition.component.clone())
            .collect();
        if !missing.is_empty() || !unknown.is_empty() {
            return Err(SeparationError::PartitionCoverage { missing, unknown });
        }
        let mut concentrate_solids = Vec::with_capacity(self.partitions.len());
        let mut tailing_solids = Vec::with_capacity(self.partitions.len());
        let mut balances = Vec::with_capacity(self.partitions.len());
        for (solid, partition) in feed.solids().iter().zip(&self.partitions) {
            let (concentrate, tailings, balance) = split(solid, partition.to_concentrate)?;
            concentrate_solids.push(concentrate);
            tailing_solids.push(tailings);
            balances.push(balance);
        }
        let (concentrate_carrier, tailing_carrier, carrier) =
            split(feed.liquid(), self.carrier_to_concentrate)?;
        let concentrate =
            SlurryStream::new_additive_volume(concentrate_solids, concentrate_carrier)?;
        let tailings = SlurryStream::new_additive_volume(tailing_solids, tailing_carrier)?;
        let streams = [feed, &concentrate, &tailings];
        let dry_rates = streams.map(|stream| stream.state().dry_solids_mass_flow.value());
        let positive_components = streams.map(|stream| {
            stream
                .solids()
                .iter()
                .filter(|solid| solid.mass_flow().value() > 0.0)
                .count()
        });
        let mut components = Vec::with_capacity(balances.len());
        for balance in balances {
            let rates = [
                balance.feed_mass_flow.value(),
                balance.concentrate_mass_flow.value(),
                balance.tailings_mass_flow.value(),
            ];
            let mut grades = [None; 3];
            for index in 0..3 {
                grades[index] =
                    dry_fraction(rates[index], dry_rates[index], positive_components[index])?;
            }
            components.push(DryComponentSeparation {
                balance,
                feed_dry_fraction: grades[0],
                concentrate_dry_fraction: grades[1],
                tailings_dry_fraction: grades[2],
            });
        }
        Ok(SeparationResult {
            dry_mass_yield: recovery(dry_rates[0], dry_rates[1], dry_rates[2])?,
            concentrate,
            tailings,
            separator: self.clone(),
            components,
            carrier,
        })
    }
}

/// Independently inspectable component flow accounting for one separation.
#[derive(Clone, Debug, PartialEq)]
pub struct ComponentSplit {
    /// Exact feed constituent identity.
    pub component: String,
    /// Feed mass rate, in kg/s.
    pub feed_mass_flow: MassFlowRate,
    /// Mass rate reporting to concentrate, in kg/s.
    pub concentrate_mass_flow: MassFlowRate,
    /// Complementary mass rate reporting to tailings, in kg/s.
    pub tailings_mass_flow: MassFlowRate,
    /// `(concentrate + tailings) - feed`, in kg/s, including arithmetic rounding.
    pub conservation_residual: MassFlowRate,
    /// Computed concentrate/feed ratio; `None` for zero feed of this constituent.
    /// This is a result of the prescribed partition, not measured plant recovery.
    pub concentrate_recovery: Option<Dimensionless>,
}

/// Dry mass fractions for one exact supplied solid identity, not elemental grade.
#[derive(Clone, Debug, PartialEq)]
pub struct DryComponentSeparation {
    /// Component flow and recovery accounting, independent of carrier mass.
    pub balance: ComponentSplit,
    /// Component mass / total dry-solid mass in the feed, or `None` if dry-empty.
    pub feed_dry_fraction: Option<Dimensionless>,
    /// Component mass / total dry-solid mass in concentrate, or `None` if dry-empty.
    pub concentrate_dry_fraction: Option<Dimensionless>,
    /// Component mass / total dry-solid mass in tailings, or `None` if dry-empty.
    pub tailings_dry_fraction: Option<Dimensionless>,
}

/// Conserved outlet streams, retained partition assumptions and component metrics.
#[derive(Clone, Debug, PartialEq)]
pub struct SeparationResult {
    /// Concentrate stream with original component density sources and contexts.
    pub concentrate: SlurryStream,
    /// Tailings stream with original component density sources and contexts.
    pub tailings: SlurryStream,
    /// Original caller-prescribed partitions and their provenance.
    pub separator: SelectiveSeparator,
    /// Per-solid flow accounting and dry fractions, in canonical identity order.
    pub components: Vec<DryComponentSeparation>,
    /// Independent carrier-liquid flow and recovery accounting.
    pub carrier: ComponentSplit,
    /// Concentrate dry mass / feed dry mass; `None` for zero dry feed.
    pub dry_mass_yield: Option<Dimensionless>,
}

fn split(
    component: &SlurryComponent,
    to_concentrate: Dimensionless,
) -> Result<(SlurryComponent, SlurryComponent, ComponentSplit)> {
    let feed = component.mass_flow().value();
    let fraction = to_concentrate.value();
    let concentrate = feed * fraction;
    let tailings = feed - concentrate;
    require(
        feed == 0.0 || fraction == 0.0 || fraction == 1.0 || (concentrate > 0.0 && tailings > 0.0),
        "positive component outlet is unrepresentable",
    )?;
    let residual = (concentrate + tailings) - feed;
    require(
        residual.is_finite(),
        "component balance residual is unrepresentable",
    )?;
    let balance = ComponentSplit {
        component: component.id().to_owned(),
        feed_mass_flow: MassFlowRate::new(feed),
        concentrate_mass_flow: MassFlowRate::new(concentrate),
        tailings_mass_flow: MassFlowRate::new(tailings),
        conservation_residual: MassFlowRate::new(residual),
        concentrate_recovery: recovery(feed, concentrate, tailings)?,
    };
    Ok((
        component.with_flow(concentrate)?,
        component.with_flow(tailings)?,
        balance,
    ))
}

fn recovery(feed: f64, concentrate: f64, tailings: f64) -> Result<Option<Dimensionless>> {
    if feed == 0.0 {
        return Ok(None);
    }
    let value = concentrate / feed;
    require(
        value.is_finite()
            && (0.0..=1.0).contains(&value)
            && (concentrate == 0.0 || value > 0.0)
            && (tailings == 0.0 || value < 1.0),
        "positive recovery or dry-mass yield is unrepresentable",
    )?;
    Ok(Some(Dimensionless::new(value)))
}

fn dry_fraction(
    component: f64,
    dry_total: f64,
    positive_components: usize,
) -> Result<Option<Dimensionless>> {
    if dry_total == 0.0 {
        return Ok(None);
    }
    let value = component / dry_total;
    require(
        value.is_finite()
            && (0.0..=1.0).contains(&value)
            && (component == 0.0 || value > 0.0)
            && (positive_components <= 1 || value < 1.0),
        "positive dry-component mass fraction is unrepresentable",
    )?;
    Ok(Some(Dimensionless::new(value)))
}

#[cfg(test)]
#[path = "separation_tests.rs"]
mod tests;
