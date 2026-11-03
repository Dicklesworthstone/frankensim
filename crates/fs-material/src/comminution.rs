//! Bounded Bond grinding energy, dry throughput, and ideal mill critical speed.
//!
//! The Bond relation uses a caller's test-specific work index, not a universal
//! mineral constant. The card retains specimen, preparation/method, closing
//! screen and provenance; its declared size window is not proof of calibration.
//! Outputs are empirical estimates, not measured plant performance.
//!
//! Formula and scope: Global Mining Guidelines Group (2021), *Determining the
//! Bond Efficiency of Industrial Grinding Circuits*, GMG01-MP-2021, §§2, 4–5:
//! <https://gmggroup.org/wp-content/uploads/2024/07/GUIDELINE__Determining-the-Bond-Efficiency-of-Industrial-Grinding-Circuits_2021.pdf>.
//! This uncorrected baseline refuses P80 below 70 µm. No circuit efficiency,
//! SAG/AG model, liberation, recovery, wear, slurry, or particle trajectory
//! model is inferred. Net power means power at the grinding drive pinion,
//! excluding auxiliaries; motor nameplate power needs explicit conversion.
//!
//! Single-source operations use bounded scalar arithmetic with correctly
//! rounded square roots. Explicit additive-energy blends traverse their declared
//! feed list; no operation performs IO or maintains mutable material history.

use fs_matdb::Provenance;
use fs_math::det;
use fs_qty::{Acceleration, Dimensionless, Frequency, Length, MassFlowRate, Power, Qty};

#[path = "comminution/blend.rs"]
mod blend;
pub use blend::{AdditiveGrindingBlend, GrindingFeed};

#[path = "comminution/inverse.rs"]
mod inverse;

/// Source-declared congruent per-pass size redistribution with conserved elements.
#[path = "comminution/breakage.rs"]
pub mod breakage;

/// Closed grinding/classification balances with explicit source applicability.
#[path = "comminution/circuit.rs"]
pub mod circuit;

/// Continuous perfectly mixed grinding with source-declared selection rates.
#[path = "comminution/population_balance.rs"]
pub mod population_balance;

/// Grinding energy per dry mass, in coherent SI joules per kilogram.
pub type SpecificGrindingEnergy = Qty<2, 0, -2, 0, 0, 0>;

const JOULES_PER_KWH: f64 = 3_600_000.0;
const KG_PER_SHORT_TON: f64 = 907.18474;
const MINIMUM_PRODUCT_M: f64 = 70.0e-6;

/// Structured admission failure; no extrapolated or silently clamped result.
#[derive(Debug, Clone, PartialEq)]
pub enum ComminutionError {
    /// Invalid physical input, missing source context, or unrepresentable result.
    Invalid(&'static str),
    /// The requested size lies outside the caller's declared interval.
    OutsideValidity {
        /// Either F80 (feed) or P80 (product).
        axis: &'static str,
        /// Requested 80%-passing size in metres.
        value_m: f64,
        /// Inclusive lower size bound in metres.
        minimum_m: f64,
        /// Inclusive upper size bound in metres.
        maximum_m: f64,
    },
}

impl core::fmt::Display for ComminutionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(what) => write!(f, "comminution: {what}"),
            Self::OutsideValidity {
                axis,
                value_m,
                minimum_m,
                maximum_m,
            } => write!(
                f,
                "comminution: {axis}={value_m} m outside [{minimum_m}, {maximum_m}] m"
            ),
        }
    }
}
impl std::error::Error for ComminutionError {}

fn require(condition: bool, what: &'static str) -> Result<(), ComminutionError> {
    if condition {
        Ok(())
    } else {
        Err(ComminutionError::Invalid(what))
    }
}
fn positive(value: f64, what: &'static str) -> Result<f64, ComminutionError> {
    require(value.is_finite() && value > 0.0, what)?;
    Ok(value)
}

/// The original dry-mass convention remains inspectable after normalization.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WorkIndexMassBasis {
    /// 1,000 kg of dry solids.
    MetricTonne,
    /// 907.18474 kg of dry solids.
    ShortTon,
}

/// A reported Bond work index with an explicit original unit convention.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BondWorkIndex {
    source_value: f64,
    source_basis: WorkIndexMassBasis,
    energy: SpecificGrindingEnergy,
}

impl BondWorkIndex {
    /// Normalize kWh per dry metric tonne (1,000 kg) to joules per kilogram.
    pub fn from_kwh_per_metric_tonne(value: f64) -> Result<Self, ComminutionError> {
        Self::from_source(value, WorkIndexMassBasis::MetricTonne, 1_000.0)
    }

    /// Normalize kWh per dry US short ton (exactly 907.18474 kg).
    pub fn from_kwh_per_short_ton(value: f64) -> Result<Self, ComminutionError> {
        Self::from_source(value, WorkIndexMassBasis::ShortTon, KG_PER_SHORT_TON)
    }

    fn from_source(
        value: f64,
        source_basis: WorkIndexMassBasis,
        kg: f64,
    ) -> Result<Self, ComminutionError> {
        positive(value, "work index must be finite and positive")?;
        let energy = positive(
            value * (JOULES_PER_KWH / kg),
            "work index conversion is unrepresentable",
        )?;
        Ok(Self {
            source_value: value,
            source_basis,
            energy: SpecificGrindingEnergy::new(energy),
        })
    }

    /// Original numeric work index, in kWh per declared source mass basis.
    #[must_use]
    pub const fn source_value(self) -> f64 {
        self.source_value
    }

    /// Original denominator, retained to prevent tonne/short-ton ambiguity.
    #[must_use]
    pub const fn source_basis(self) -> WorkIndexMassBasis {
        self.source_basis
    }

    /// Normalized reference work index in J/kg of dry solids.
    #[must_use]
    pub const fn specific_energy(self) -> SpecificGrindingEnergy {
        self.energy
    }
}

/// Explicit applicability window on the two 80%-passing sizes, in metres.
/// This baseline has no fine-grinding correction below 70 µm.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BondSizeDomain {
    feed: [Length; 2],
    product: [Length; 2],
}

impl BondSizeDomain {
    /// Admit inclusive finite feed/product bounds with some possible size reduction.
    pub fn new(feed: [Length; 2], product: [Length; 2]) -> Result<Self, ComminutionError> {
        for range in [feed, product] {
            positive(range[0].value(), "minimum size must be finite and positive")?;
            positive(range[1].value(), "maximum size must be finite and positive")?;
            require(range[0] <= range[1], "size domain has reversed bounds")?;
        }
        require(
            product[0].value() >= MINIMUM_PRODUCT_M,
            "uncorrected Bond model requires P80 >= 70 micrometres",
        )?;
        require(
            feed[1] > product[0],
            "size domain contains no reduction with F80 > P80",
        )?;
        Ok(Self { feed, product })
    }

    /// Inclusive feed-size bounds, in metres.
    #[must_use]
    pub const fn feed_bounds(self) -> [Length; 2] {
        self.feed
    }

    /// Inclusive product-size bounds, in metres.
    #[must_use]
    pub const fn product_bounds(self) -> [Length; 2] {
        self.product
    }

    fn admit(self, feed: Length, product: Length) -> Result<(), ComminutionError> {
        for (axis, value, bounds) in [("F80", feed, self.feed), ("P80", product, self.product)] {
            positive(value.value(), "F80 and P80 must be finite and positive")?;
            if value < bounds[0] || value > bounds[1] {
                return Err(ComminutionError::OutsideValidity {
                    axis,
                    value_m: value.value(),
                    minimum_m: bounds[0].value(),
                    maximum_m: bounds[1].value(),
                });
            }
        }
        require(
            feed > product,
            "grinding requires F80 > P80; zero and inverse reduction refuse",
        )
    }
}

/// Original work-index test context. No mineral name supplies a default Wi.
#[derive(Clone, Debug, PartialEq)]
pub struct BondTestContext {
    /// Located composite/sample identity, or an explicitly labelled authored fixture.
    pub specimen: String,
    /// Test procedure and feed preparation; identifies the kind of work-index test.
    pub method: String,
    /// Test closing-screen aperture, distinct from the product's 80%-passing size.
    pub closing_screen: Length,
    /// Original lab report or authored scenario, usage terms, and optional artifact hash.
    pub provenance: Provenance,
}

/// Immutable empirical grinding card with a retained original test declaration.
#[derive(Clone, Debug, PartialEq)]
pub struct BondGrindingCard {
    index: BondWorkIndex,
    test: BondTestContext,
    domain: BondSizeDomain,
}

impl BondGrindingCard {
    /// Admit a reported work index and retain its declared test and size support.
    pub fn new(
        index: BondWorkIndex,
        test: BondTestContext,
        domain: BondSizeDomain,
    ) -> Result<Self, ComminutionError> {
        require(
            !test.specimen.trim().is_empty(),
            "work index requires a specimen identity",
        )?;
        require(
            !test.method.trim().is_empty(),
            "work index requires its test method and preparation",
        )?;
        require(
            !test.provenance.source.trim().is_empty(),
            "work index requires source provenance",
        )?;
        require(
            !test.provenance.license.trim().is_empty(),
            "work index requires source usage terms",
        )?;
        positive(
            test.closing_screen.value(),
            "test closing screen must be finite and positive",
        )?;
        Ok(Self {
            index,
            test,
            domain,
        })
    }

    /// Original work index and source mass basis.
    #[must_use]
    pub const fn work_index(&self) -> BondWorkIndex {
        self.index
    }

    /// Immutable original specimen, preparation, screen and source declaration.
    #[must_use]
    pub fn test_context(&self) -> &BondTestContext {
        &self.test
    }

    /// Caller-declared size support; not automatically inferred from the screen.
    #[must_use]
    pub const fn domain(&self) -> BondSizeDomain {
        self.domain
    }

    /// Evaluate the uncorrected Bond baseline in J/kg of dry solids.
    ///
    /// With sizes in micrometres, W = 10 Wi (1/sqrt(P80) - 1/sqrt(F80)).
    /// Here sizes are SI metres, so the dimensionless factor is
    /// sqrt(100e-6/P80) * (1-sqrt(P80/F80)). The latter difference is
    /// rationalized to retain small positive reductions when F80 approaches P80.
    /// The caller must choose test/size applicability; this is no assertion that
    /// a measured Wi transfers to another ore, blend, circuit or preparation.
    pub fn specific_energy(
        &self,
        feed: Length,
        product: Length,
    ) -> Result<SpecificGrindingEnergy, ComminutionError> {
        self.domain.admit(feed, product)?;
        let f = feed.value();
        let p = product.value();
        let factor = det::sqrt(100.0e-6 / p) * (((f - p) / f) / (1.0 + det::sqrt(p / f)));
        let energy = positive(
            self.index.energy.value() * factor,
            "grinding energy is unrepresentable",
        )?;
        Ok(SpecificGrindingEnergy::new(energy))
    }

    /// Estimate dry-solids mass flow from pinion power, excluding auxiliaries.
    /// Utilization is an explicit fraction of the reporting interval; 0 means
    /// downtime. No efficiency, circulating-load or slurry-water factor is inferred.
    pub fn dry_throughput(
        &self,
        feed: Length,
        product: Length,
        net_pinion_power: Power,
        utilization: Dimensionless,
    ) -> Result<MassFlowRate, ComminutionError> {
        let energy = self.specific_energy(feed, product)?;
        throughput_from_net_power(energy, net_pinion_power, utilization)
    }
}

/// Dry mass flow (kg/s) = available grinding power (W) / specific energy (J/kg).
/// The caller explicitly supplies pinion power and reporting-period utilization.
pub fn throughput_from_net_power(
    energy: SpecificGrindingEnergy,
    net_pinion_power: Power,
    utilization: Dimensionless,
) -> Result<MassFlowRate, ComminutionError> {
    positive(
        energy.value(),
        "specific grinding energy must be finite and positive",
    )?;
    positive(
        net_pinion_power.value(),
        "net pinion power must be finite and positive",
    )?;
    let u = utilization.value();
    require(
        u.is_finite() && (0.0..=1.0).contains(&u),
        "utilization must be finite and in [0,1]",
    )?;
    if u == 0.0 {
        return Ok(MassFlowRate::new(0.0));
    }
    let rate = positive(
        net_pinion_power.value() * u / energy.value(),
        "dry throughput is unrepresentable",
    )?;
    Ok(MassFlowRate::new(rate))
}

/// Ideal centrifuging threshold, in revolutions per second (multiply by 60 for rpm).
///
/// A single ball touching the inside liner has orbit radius (D-d)/2. Balancing
/// gravity and centripetal acceleration gives n_c = sqrt(2g/(D-d))/(2*pi).
/// D is the free diameter INSIDE the liner, d is the ball diameter, and g is
/// supplied explicitly. This ideal threshold is not an operating-speed
/// recommendation and omits lifters, charge interactions, slip and slurry.
pub fn ball_mill_critical_speed(
    inside_diameter: Length,
    ball_diameter: Length,
    gravity: Acceleration,
) -> Result<Frequency, ComminutionError> {
    let d_mill = positive(
        inside_diameter.value(),
        "inside mill diameter must be finite and positive",
    )?;
    let d_ball = positive(
        ball_diameter.value(),
        "ball diameter must be finite and positive",
    )?;
    let g = positive(gravity.value(), "gravity must be finite and positive")?;
    require(
        d_ball < d_mill,
        "ball diameter must be smaller than the free mill diameter",
    )?;
    // Split the square roots before division: 2*g or g/(D-d) can overflow
    // even when the final frequency is finite.
    let frequency = (det::sqrt(g) / core::f64::consts::TAU * core::f64::consts::SQRT_2)
        / det::sqrt(d_mill - d_ball);
    Ok(Frequency::new(positive(
        frequency,
        "critical speed is unrepresentable",
    )?))
}

#[cfg(test)]
#[path = "comminution/tests.rs"]
mod tests;
