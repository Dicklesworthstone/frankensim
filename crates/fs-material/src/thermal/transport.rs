//! Prepared thermal transport along one fixed-context temperature path.
//!
//! Admission uses the canonical material source-line evaluator, retaining every
//! source breakpoint and usage receipt. Runtime evaluation performs one binary
//! search and interpolates rho, c and k separately before deriving rho*c,
//! diffusivity and effusivity. It performs no allocation, hashing or database
//! query. Products and ratios are never independently linearly interpolated.
//!
//! The caller selects the specific-heat constraint and heat-flow direction.
//! Scalar k does not establish isotropy; c is not converted between Cp, Cv or
//! fixed-strain heat capacities. No phase change, enthalpy integration, moisture
//! evolution or uncertainty propagation is inferred. Existing caloric and
//! conservative thermal owners remain responsible for those operations.

/// Local thermal response from directly supplied volumetric heat capacity.
#[path = "transport/volumetric.rs"]
pub mod volumetric;

/// Constant-mass density from explicitly selected free length-change sources.
#[path = "transport/expansion.rs"]
pub mod expansion;

/// Sourced heat transport using constant-mass free-expansion density.
#[path = "transport/expanding.rs"]
pub mod expanding;

/// Additive sensible heat storage for explicitly co-temperature constituents.
#[path = "transport/storage.rs"]
pub mod storage;

/// Source-qualified ordered multilayer conductivity and contact response.
#[path = "transport/stack.rs"]
pub mod stack;

// Shared affine-capacity integration and inversion, independent of its units.
#[path = "transport/capacity.rs"]
mod capacity;

use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{
    ClaimSelection, MaterialCard, PropertyKey, PropertyLineSlice, PropertySupportError, QueryPoint,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec};

const DENSITY: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const SPECIFIC_HEAT: Dims = Dims([2, 0, -2, -1, 0, 0]);
const CONDUCTIVITY: Dims = Dims([1, 1, -3, -1, 0, 0]);
const TEMPERATURE: Dims = Dims([0, 0, 0, 1, 0, 0]);
const INVERSE_TIME: Dims = Dims([0, 0, -1, 0, 0, 0]);

/// Bound inherited from three canonical source lines, without adaptive resampling.
pub const MAX_THERMAL_TRANSPORT_KNOTS: usize = 3 * fs_matdb::MAX_PROPERTY_LINE_KNOTS;

/// Exact source roles, fixed conditions and independent claim selections.
#[derive(Clone, Copy, Debug)]
pub struct ThermalTransportSourceQuery<'a> {
    /// Bulk density [kg/m3].
    pub density: &'a PropertyKey,
    /// Specific heat appropriate to the consumer's caloric constraint [J/(kg K)].
    pub specific_heat: &'a PropertyKey,
    /// Scalar conductivity along the selected heat-flow direction [W/(m K)].
    pub conductivity: &'a PropertyKey,
    /// The only changing source coordinate, explicitly typed positive absolute kelvin.
    pub temperature_axis: &'a str,
    /// Lower temperature and every required source condition.
    pub lower: &'a QueryPoint,
    /// Higher temperature; all other coordinates and descriptors stay fixed.
    pub upper: &'a QueryPoint,
    /// Density, specific heat, conductivity selections, in that order.
    pub selections: [ClaimSelection; 3],
}

/// Admission or evaluation failure; no partially prepared path is published.
#[derive(Clone, Debug, PartialEq)]
pub enum ThermalTransportError {
    /// Invalid quantity role, source support, temperature or numerical result.
    Invalid(&'static str),
    /// Original canonical source refusal, including a witness query point.
    Source {
        /// Required property whose source could not support the entire path.
        property: String,
        /// Unmodified source-support diagnosis.
        error: Box<PropertySupportError>,
    },
}

impl core::fmt::Display for ThermalTransportError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(reason) => write!(f, "thermal transport: {reason}"),
            Self::Source { property, error } => write!(f, "thermal transport {property}: {error}"),
        }
    }
}
impl std::error::Error for ThermalTransportError {}

/// One coherent numerical evaluation of the selected source models.
/// Original uncertainty and receipts remain on the prepared path; this is not
/// a newly measured property claim or a new canonical point-query receipt.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalTransportSample {
    path_identity: ContentHash,
    temperature_k: f64,
    values: [f64; 3],
    volumetric_heat_capacity: f64,
    diffusivity: f64,
    effusivity: f64,
}

impl ThermalTransportSample {
    /// Original prepared model and fixed source-condition identity.
    #[must_use]
    pub const fn path_identity(&self) -> ContentHash {
        self.path_identity
    }
    /// Evaluation temperature [K].
    #[must_use]
    pub const fn temperature_k(&self) -> f64 {
        self.temperature_k
    }
    /// Bulk density [kg/m3].
    #[must_use]
    pub const fn density_kg_m3(&self) -> f64 {
        self.values[0]
    }
    /// Selected specific heat [J/(kg K)].
    #[must_use]
    pub const fn specific_heat_j_kg_k(&self) -> f64 {
        self.values[1]
    }
    /// Conductivity along the selected source direction [W/(m K)].
    #[must_use]
    pub const fn conductivity_w_m_k(&self) -> f64 {
        self.values[2]
    }
    /// Local rho*c [J/(m3 K)], without mechanical or latent-heat terms.
    #[must_use]
    pub const fn volumetric_heat_capacity_j_m3_k(&self) -> f64 {
        self.volumetric_heat_capacity
    }
    /// Local k/(rho*c) [m2/s].
    #[must_use]
    pub const fn diffusivity_m2_s(&self) -> f64 {
        self.diffusivity
    }
    /// Local sqrt(k*rho*c) [J/(m2 K sqrt(s))].
    #[must_use]
    pub const fn effusivity_j_m2_k_sqrt_s(&self) -> f64 {
        self.effusivity
    }
}

/// Source values on the union grid; only the three independent inputs interpolate.
#[derive(Clone, Copy, Debug, PartialEq)]
struct TransportKnot {
    temperature: f64,
    values: [f64; 3],
}

/// Immutable prepared model usable across repeated thermal evaluations.
#[derive(Clone, Debug, PartialEq)]
pub struct PreparedThermalTransport {
    material_card: ContentHash,
    lower: QueryPoint,
    upper: QueryPoint,
    sources: [PropertyLineSlice; 3],
    knots: Vec<TransportKnot>,
    identity: ContentHash,
}

impl PreparedThermalTransport {
    /// Admit the complete fixed-context path before preparing its numerical view.
    ///
    /// Reuses canonical source selection, explicit pins and continuous support,
    /// including connected fixed-context cuts through sampled surfaces. Unknown
    /// temperature support, source gaps, discrete-only spans, changing selection
    /// and conflicting observations refuse. Each selected claim must explicitly
    /// bound the temperature axis: an unconditioned constant is not a thermal law.
    ///
    /// Construction is bounded by the existing source-line budgets. All original
    /// receipts, marginal uncertainties and card identity survive; no covariance
    /// or propagated uncertainty is invented. Re-prepare when material identity,
    /// any non-temperature condition or the selected source model changes.
    pub fn try_new(
        card: &MaterialCard,
        query: ThermalTransportSourceQuery<'_>,
    ) -> Result<Self, ThermalTransportError> {
        let invalid = ThermalTransportError::Invalid;
        let keys = [query.density, query.specific_heat, query.conductivity];
        let conductivity_quantity = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::ThermalConductivity,
            ValueForm::Static,
        ));
        for (index, (key, dims)) in keys
            .iter()
            .zip([DENSITY, SPECIFIC_HEAT, CONDUCTIVITY])
            .enumerate()
        {
            if (key.quantity() != QuantitySpec::dimensional(dims)
                && !(index == 2 && key.quantity() == conductivity_quantity))
                || key.hardness_test().is_some()
                || key.elastic_component().is_some()
                || key.strain_component().is_some()
                || key.stress_component().is_some()
            {
                return Err(invalid(
                    "source roles require scalar density, specific heat and static thermal conductivity",
                ));
            }
        }
        for point in [query.lower, query.upper] {
            admit_temperature_point(point, query.temperature_axis)?;
        }
        if query.lower.axes()[query.temperature_axis] >= query.upper.axes()[query.temperature_axis]
        {
            return Err(invalid(
                "the prepared temperature interval must have positive width",
            ));
        }
        let resolve = |index: usize| {
            let source = source_line(
                card,
                keys[index],
                query.temperature_axis,
                query.lower,
                query.upper,
                query.selections[index],
            )?;
            if source
                .knots
                .iter()
                .any(|(_, answer)| !positive(answer.evidence.value.value))
            {
                return Err(invalid(
                    "density, specific heat and conductivity must remain finite and positive",
                ));
            }
            Ok(source)
        };
        let sources = [resolve(0)?, resolve(1)?, resolve(2)?];
        let mut temperatures: Vec<_> = sources
            .iter()
            .flat_map(|s| s.knots.iter().map(|(t, _)| *t))
            .collect();
        temperatures.sort_by(f64::total_cmp);
        temperatures.dedup();
        if temperatures.len() > MAX_THERMAL_TRANSPORT_KNOTS {
            return Err(invalid(
                "thermal source union exceeds its breakpoint budget",
            ));
        }
        let mut identity =
            DomainHasher::new("org.frankensim.fs-material.prepared-thermal-transport.v1");
        let material_card = card.content_hash();
        identity.update(material_card.as_bytes());
        for source in &sources {
            identity.update(&(source.knots.len() as u64).to_le_bytes());
            for (_, answer) in &source.knots {
                identity.update(answer.receipt.content_hash().as_bytes());
            }
        }
        let identity = identity.finalize();
        let mut knots = Vec::with_capacity(temperatures.len());
        for temperature in temperatures {
            let values = std::array::from_fn(|i| source_value(&sources[i], temperature));
            checked_sample(identity, temperature, values)?;
            knots.push(TransportKnot {
                temperature,
                values,
            });
        }
        Ok(Self {
            material_card,
            lower: query.lower.clone(),
            upper: query.upper.clone(),
            sources,
            knots,
            identity,
        })
    }

    /// Evaluate the retained piecewise-affine inputs with one binary search.
    /// No query, allocation, mutation or extrapolation occurs in this hot path.
    /// Roundoff can differ from a direct source query after union-grid splitting;
    /// this preserves the same mathematical source model, not bitwise equality.
    /// Preparation proves source support and checks numerical values at knots.
    /// An interior product or ratio can still overflow or underflow for extreme
    /// inputs; evaluation refuses that sample without changing the prepared path.
    pub fn at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<ThermalTransportSample, ThermalTransportError> {
        let span = self.temperature_span_k();
        if !temperature_k.is_finite() || temperature_k < span[0] || temperature_k > span[1] {
            return Err(ThermalTransportError::Invalid(
                "temperature is outside the prepared source interval",
            ));
        }
        let index = self
            .knots
            .partition_point(|knot| knot.temperature < temperature_k);
        let upper = self.knots[index];
        let values = if upper.temperature == temperature_k {
            upper.values
        } else {
            let lower = self.knots[index - 1];
            let weight =
                (temperature_k - lower.temperature) / (upper.temperature - lower.temperature);
            std::array::from_fn(|i| lower.values[i] + weight * (upper.values[i] - lower.values[i]))
        };
        checked_sample(self.identity, temperature_k, values)
    }

    /// Prepared path identity, including original source roles and receipts.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    /// Immutable material-card identity used at admission.
    #[must_use]
    pub const fn material_card_identity(&self) -> ContentHash {
        self.material_card
    }
    /// Exact fixed source context at the lower and upper temperatures.
    #[must_use]
    pub fn query_endpoints(&self) -> [&QueryPoint; 2] {
        [&self.lower, &self.upper]
    }
    /// Source lines in density, specific-heat, conductivity order, with original receipts.
    #[must_use]
    pub const fn source_lines(&self) -> &[PropertyLineSlice; 3] {
        &self.sources
    }
    /// Inclusive prepared temperature range [K].
    #[must_use]
    pub fn temperature_span_k(&self) -> [f64; 2] {
        [
            self.knots[0].temperature,
            self.knots[self.knots.len() - 1].temperature,
        ]
    }
    /// Every retained source breakpoint in increasing temperature order.
    pub fn temperatures_k(&self) -> impl ExactSizeIterator<Item = f64> + '_ {
        self.knots.iter().map(|knot| knot.temperature)
    }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

fn checked_sample(
    identity: ContentHash,
    temperature: f64,
    values: [f64; 3],
) -> Result<ThermalTransportSample, ThermalTransportError> {
    let volumetric_heat_capacity = values[0] * values[1];
    let [diffusivity, effusivity] = transport_values(volumetric_heat_capacity, values[2])?;
    if !values.into_iter().all(positive) {
        return Err(ThermalTransportError::Invalid(
            "thermal inputs or derived transport are not positive finite representable values",
        ));
    }
    Ok(ThermalTransportSample {
        path_identity: identity,
        temperature_k: temperature,
        values,
        volumetric_heat_capacity,
        diffusivity,
        effusivity,
    })
}

fn transport_values(capacity: f64, conductivity: f64) -> Result<[f64; 2], ThermalTransportError> {
    let diffusivity = conductivity / capacity;
    // Avoid overflow in k*C when its square root is representable.
    let effusivity = conductivity.sqrt() * capacity.sqrt();
    if ![capacity, conductivity, diffusivity, effusivity]
        .into_iter()
        .all(positive)
    {
        return Err(ThermalTransportError::Invalid(
            "thermal inputs or derived transport are not positive finite representable values",
        ));
    }
    Ok([diffusivity, effusivity])
}

fn admit_temperature_point(point: &QueryPoint, axis: &str) -> Result<(), ThermalTransportError> {
    let absolute = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    if !point.axes().get(axis).is_some_and(|&t| positive(t))
        || !point
            .axis_quantities()
            .get(axis)
            .is_some_and(|q| *q == absolute || *q == QuantitySpec::dimensional(TEMPERATURE))
    {
        return Err(ThermalTransportError::Invalid(
            "the source path requires positive absolute temperatures in kelvin",
        ));
    }
    if point.axes().iter().any(|(name, value)| {
        *value != 0.0
            && (name == "omega"
                || name == "frequency"
                || point
                    .axis_quantities()
                    .get(name)
                    .is_some_and(|q| q.dims() == INVERSE_TIME))
    }) {
        return Err(ThermalTransportError::Invalid(
            "nonzero frequency or rate coordinates require an explicitly coupled thermal owner",
        ));
    }
    Ok(())
}

fn source_line(
    card: &MaterialCard,
    key: &PropertyKey,
    axis: &str,
    lower: &QueryPoint,
    upper: &QueryPoint,
    selection: ClaimSelection,
) -> Result<PropertyLineSlice, ThermalTransportError> {
    let invalid = ThermalTransportError::Invalid;
    let source = card
        .claims()
        .query_line_typed(key, axis, lower, upper, selection)
        .map_err(|error| ThermalTransportError::Source {
            property: key.name().to_owned(),
            error: Box::new(error),
        })?;
    let claim = card
        .claims()
        .claim(source.knots[0].1.receipt.selected)
        .expect("canonical source line selected an existing claim");
    if !claim.validity.bounds().contains_key(axis) {
        return Err(invalid(
            "each thermal source must explicitly bound the temperature coordinate",
        ));
    }
    let temperature_difference = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::TemperatureDifference,
        ValueForm::Static,
    ));
    if claim.validity.bounds().keys().any(|name| {
        name != axis
            && (name == "T"
                || lower
                    .axis_quantities()
                    .get(name)
                    .is_some_and(|q| q.dims() == TEMPERATURE && *q != temperature_difference))
    }) {
        return Err(invalid(
            "a second source temperature requires an explicitly coupled thermal owner",
        ));
    }
    if !claim.value.axes().all(|(name, dims)| {
        if name == axis {
            dims == TEMPERATURE
        } else {
            lower.axis_quantities().get(name).is_some_and(|q| {
                q.dims() == dims && (dims != TEMPERATURE || *q == temperature_difference)
            }) && (dims != INVERSE_TIME || lower.axes().get(name) == Some(&0.0))
        }
    }) {
        return Err(invalid(
            "source value axes must match the typed temperature and fixed static conditions",
        ));
    }
    Ok(source)
}

fn source_value(source: &PropertyLineSlice, temperature: f64) -> f64 {
    let index = source.knots.partition_point(|(t, _)| *t < temperature);
    let (hi, answer) = &source.knots[index];
    let high = answer.evidence.value.value;
    if *hi == temperature {
        return high;
    }
    let (lo, answer) = &source.knots[index - 1];
    let low = answer.evidence.value.value;
    low + (temperature - lo) / (hi - lo) * (high - low)
}
