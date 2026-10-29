//! Additive sensible heat storage at one explicitly shared temperature.
//!
//! Components retain their own material, conditioning, caloric constraint and
//! source evidence. Fixed masses multiply specific capacities; fixed reference
//! volumes multiply native volumetric capacities. Neither density nor thermal
//! conductivity is needed. The caller declares compatible work constraints and
//! a co-temperature approximation; this is not a homogenized conductivity,
//! chemical mixture law, spatial equilibration model or phase-change model.

#[path = "storage/timeline.rs"]
mod timeline;
pub use timeline::HeatStorageThermalLaw;

use super::{
    SPECIFIC_HEAT, ThermalTransportError, admit_temperature_point, capacity::CapacityCurve,
    positive, source_line, source_value, volumetric::VOLUMETRIC_HEAT_CAPACITY_DIMS,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{ClaimSelection, MaterialCard, PropertyKey, PropertyLineSlice, QueryPoint};
use fs_qty::QuantitySpec;

/// Bound preparation and retention independently of the number of queries.
pub const MAX_HEAT_STORAGE_COMPONENTS: usize = 64;
/// Every retained breakpoint comes from an original admitted source line.
pub const MAX_HEAT_STORAGE_KNOTS: usize =
    MAX_HEAT_STORAGE_COMPONENTS * fs_matdb::MAX_PROPERTY_LINE_KNOTS;

/// Explicit physical amount and matching source-capacity basis.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum HeatStorageAmount {
    /// Fixed constituent mass multiplying capacity in J/(kg K).
    MassKg(f64),
    /// Fixed reference volume multiplying native capacity in J/(m3 K).
    ReferenceVolumeM3(f64),
}

impl HeatStorageAmount {
    fn value(self) -> f64 {
        match self {
            Self::MassKg(value) | Self::ReferenceVolumeM3(value) => value,
        }
    }
    fn tag(self) -> u8 {
        match self {
            Self::MassKg(_) => 0,
            Self::ReferenceVolumeM3(_) => 1,
        }
    }
}

/// One constituent and its complete fixed-context source selection.
#[derive(Clone, Copy, Debug)]
pub struct HeatStorageComponentQuery<'a> {
    /// Stable constituent identity; duplicate IDs refuse.
    pub id: u64,
    pub card: &'a MaterialCard,
    pub capacity: &'a PropertyKey,
    pub amount: HeatStorageAmount,
    pub temperature_axis: &'a str,
    pub lower: &'a QueryPoint,
    pub upper: &'a QueryPoint,
    pub selection: ClaimSelection,
}

/// Original source evidence and amount, never an independent measured mixture.
#[derive(Clone, Debug, PartialEq)]
pub struct HeatStorageComponent {
    id: u64,
    material_card: ContentHash,
    amount: HeatStorageAmount,
    temperature_axis: String,
    endpoints: [QueryPoint; 2],
    source: PropertyLineSlice,
}

impl HeatStorageComponent {
    #[must_use]
    pub const fn id(&self) -> u64 {
        self.id
    }
    #[must_use]
    pub const fn material_card_identity(&self) -> ContentHash {
        self.material_card
    }
    #[must_use]
    pub const fn amount(&self) -> HeatStorageAmount {
        self.amount
    }
    #[must_use]
    pub fn temperature_axis(&self) -> &str {
        &self.temperature_axis
    }
    #[must_use]
    pub const fn query_endpoints(&self) -> &[QueryPoint; 2] {
        &self.endpoints
    }
    #[must_use]
    pub const fn source_line(&self) -> &PropertyLineSlice {
        &self.source
    }
}

/// Prepared total capacity C(T) = sum(m_i c_i(T)) + sum(V_j Cvol_j(T)).
///
/// Components are sorted by stable ID before summation and hashing. The path
/// retains every breakpoint inside the common support interval; no interval
/// extrapolation or cross-population interpolation occurs. Preparation costs
/// O(N K log K_source + K log K), with N bounded by 64 and K by the source-line
/// budget. Runtime capacity uses one binary search. Energy, inverse and entropy
/// use bounded breakpoint traversal, independent of the number of components,
/// with no allocation, database queries, hashing or numerical root iteration.
#[derive(Clone, Debug, PartialEq)]
pub struct PreparedHeatStorage {
    components: Vec<HeatStorageComponent>,
    // Absolute kelvin, total capacity J/K; positive and ordered.
    knots: Vec<[f64; 2]>,
    identity: ContentHash,
}

impl PreparedHeatStorage {
    pub fn try_new(
        queries: &[HeatStorageComponentQuery<'_>],
    ) -> Result<Self, ThermalTransportError> {
        if queries.is_empty() || queries.len() > MAX_HEAT_STORAGE_COMPONENTS {
            return Err(invalid(
                "heat storage needs one to 64 explicit constituents",
            ));
        }
        let mut ordered: Vec<_> = queries.iter().collect();
        ordered.sort_by_key(|query| query.id);
        if ordered.windows(2).any(|w| w[0].id == w[1].id) {
            return Err(invalid(
                "heat storage constituent identities must be distinct",
            ));
        }
        let mut components = Vec::with_capacity(ordered.len());
        let mut low: f64 = 0.0;
        let mut high = f64::INFINITY;
        for query in ordered {
            admit(query)?;
            let source = source_line(
                query.card,
                query.capacity,
                query.temperature_axis,
                query.lower,
                query.upper,
                query.selection,
            )?;
            if source
                .knots
                .iter()
                .any(|(_, answer)| !positive(answer.evidence.value.value))
            {
                return Err(invalid(
                    "all selected heat capacities must be finite and positive",
                ));
            }
            low = low.max(query.lower.axes()[query.temperature_axis]);
            high = high.min(query.upper.axes()[query.temperature_axis]);
            components.push(HeatStorageComponent {
                id: query.id,
                material_card: query.card.content_hash(),
                amount: query.amount,
                temperature_axis: query.temperature_axis.to_owned(),
                endpoints: [query.lower.clone(), query.upper.clone()],
                source,
            });
        }
        if low >= high {
            return Err(invalid(
                "constituents have no positive-width common temperature support",
            ));
        }
        let mut temperatures: Vec<_> = components
            .iter()
            .flat_map(|c| c.source.knots.iter().map(|(t, _)| *t))
            .filter(|t| *t >= low && *t <= high)
            .collect();
        // Each intersection endpoint is an original selected source endpoint.
        temperatures.sort_by(f64::total_cmp);
        temperatures.dedup();
        if temperatures.len() < 2 || temperatures.len() > MAX_HEAT_STORAGE_KNOTS {
            return Err(invalid(
                "heat storage exceeds its admitted breakpoint budget",
            ));
        }
        let mut knots = Vec::with_capacity(temperatures.len());
        for temperature in temperatures {
            let mut capacity = 0.0;
            let mut correction = 0.0;
            for component in &components {
                let contribution =
                    component.amount.value() * source_value(&component.source, temperature);
                if !positive(contribution) {
                    return Err(invalid(
                        "constituent capacity is not positive and representable",
                    ));
                }
                let corrected = contribution - correction;
                let next = capacity + corrected;
                if !positive(next) {
                    return Err(invalid(
                        "assembled capacity is not finite and representable",
                    ));
                }
                correction = (next - capacity) - corrected;
                capacity = next;
            }
            knots.push([temperature, capacity]);
        }
        let mut hash = DomainHasher::new("org.frankensim.fs-material.prepared-heat-storage.v1");
        hash.update(&(components.len() as u64).to_le_bytes());
        for component in &components {
            hash.update(&component.id.to_le_bytes());
            hash.update(component.material_card.as_bytes());
            hash.update(&[component.amount.tag()]);
            hash.update(&component.amount.value().to_bits().to_le_bytes());
            hash.update(&(component.temperature_axis.len() as u64).to_le_bytes());
            hash.update(component.temperature_axis.as_bytes());
            hash.update(&(component.source.knots.len() as u64).to_le_bytes());
            for (_, answer) in &component.source.knots {
                hash.update(answer.receipt.content_hash().as_bytes());
            }
        }
        Ok(Self {
            components,
            knots,
            identity: hash.finalize(),
        })
    }

    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    #[must_use]
    pub fn components(&self) -> &[HeatStorageComponent] {
        &self.components
    }
    #[must_use]
    pub fn temperature_span_k(&self) -> [f64; 2] {
        [self.knots[0][0], self.knots[self.knots.len() - 1][0]]
    }
    pub fn temperatures_k(&self) -> impl ExactSizeIterator<Item = f64> + '_ {
        self.knots.iter().map(|knot| knot[0])
    }
    pub fn heat_capacity_j_k(&self, temperature_k: f64) -> Result<f64, ThermalTransportError> {
        CapacityCurve::new(&self.knots).capacity(temperature_k)
    }
    /// Signed sum of constituent sensible storage changes in joules.
    pub fn energy_change_j(
        &self,
        initial_k: f64,
        final_k: f64,
    ) -> Result<f64, ThermalTransportError> {
        CapacityCurve::new(&self.knots).energy_change(initial_k, final_k)
    }
    /// Unique temperature under the same co-temperature, fixed-amount model.
    pub fn temperature_after_energy_change(
        &self,
        initial_k: f64,
        change_j: f64,
    ) -> Result<f64, ThermalTransportError> {
        CapacityCurve::new(&self.knots).temperature_after_energy_change(initial_k, change_j)
    }
    /// Reversible sensible thermal entropy change, integral C(T)/T dT [J/K].
    /// This explicit interpretation supplies no mixing entropy, absolute datum,
    /// chemical reaction, latent heat or measurement uncertainty propagation.
    pub fn sensible_entropy_change_j_k(
        &self,
        initial_k: f64,
        final_k: f64,
    ) -> Result<f64, ThermalTransportError> {
        CapacityCurve::new(&self.knots).entropy_change(initial_k, final_k)
    }
}

fn admit(query: &HeatStorageComponentQuery<'_>) -> Result<(), ThermalTransportError> {
    let dims = match query.amount {
        HeatStorageAmount::MassKg(_) => SPECIFIC_HEAT,
        HeatStorageAmount::ReferenceVolumeM3(_) => VOLUMETRIC_HEAT_CAPACITY_DIMS,
    };
    if !positive(query.amount.value())
        || query.capacity.quantity() != QuantitySpec::dimensional(dims)
        || query.capacity.hardness_test().is_some()
        || query.capacity.elastic_component().is_some()
        || query.capacity.strain_component().is_some()
        || query.capacity.stress_component().is_some()
    {
        return Err(invalid(
            "heat capacity must match its positive physical mass or volume basis",
        ));
    }
    for point in [query.lower, query.upper] {
        admit_temperature_point(point, query.temperature_axis)?;
    }
    if query.lower.axis_quantities() != query.upper.axis_quantities()
        || !query.lower.axes().keys().eq(query.upper.axes().keys())
        || query.lower.axes().iter().any(|(axis, value)| {
            axis != query.temperature_axis && query.upper.axes()[axis] != *value
        })
        || query.lower.axes()[query.temperature_axis] >= query.upper.axes()[query.temperature_axis]
    {
        return Err(invalid(
            "each capacity path must retain its complete fixed source context",
        ));
    }
    Ok(())
}

fn invalid(reason: &'static str) -> ThermalTransportError {
    ThermalTransportError::Invalid(reason)
}
