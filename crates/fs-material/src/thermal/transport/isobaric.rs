//! Source Cp(T) on the existing fixed-mass, constant-pressure heat timeline.
//!
//! This caloric adapter needs neither conductivity nor a density-temperature law.
//! The original thermal owner evolves m*h and exchanges heat. Geometry remains
//! caller-authored heat-transfer geometry; this law supplies no current volume,
//! pressure evolution, phase change, mechanical storage or separate p*dV account.
//! The caller selects true single-phase isobaric Cp at fixed composition. A
//! same-dimension Cv, apparent capacity or strained-body capacity cannot substitute.
//!
//! The shared capacity kernel owns integration, entropy and inversion. Original
//! source knots, receipts and uncertainty survive; no resampled caloric chart,
//! second integrator or material-name branch is introduced.

use super::{
    SPECIFIC_HEAT, ThermalTransportError, admit_temperature_point, capacity::CapacityCurve,
    positive, source_line,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{ClaimSelection, MaterialCard, PropertyKey, PropertyLineSlice, QueryPoint};
use fs_qty::{Pressure, QuantitySpec};
use fs_thermal::enthalpy::{DepositionError, ThermalLaw, ThermalStorageBasis};

/// Explicit Cp role and one source-supported isobar at fixed composition.
#[derive(Clone, Copy, Debug)]
pub struct IsobaricHeatCapacityQuery<'a> {
    pub specific_heat_cp: &'a PropertyKey,
    pub temperature_axis: &'a str,
    /// Positive absolute pressure, typed in pascals and explicitly bounded by
    /// the selected source. A free extra query coordinate proves no support.
    pub pressure_axis: &'a str,
    pub lower: &'a QueryPoint,
    pub upper: &'a QueryPoint,
    pub selection: ClaimSelection,
}

/// Immutable single-phase caloric law. Zero h is at the caller's datum T.
/// Source Cp alone does not identify an absolute formation-energy datum.
/// Materials exchanging mass must use compatible externally established datums.
#[derive(Clone, Debug, PartialEq)]
pub struct IsobaricHeatCapacityLaw {
    material_card: ContentHash,
    endpoints: [QueryPoint; 2],
    source: PropertyLineSlice,
    knots: Vec<[f64; 2]>,
    pressure_pa: f64,
    datum_temperature_k: f64,
    enthalpy_span: [f64; 2],
    maximum_slope: f64,
    identity: ContentHash,
}

impl IsobaricHeatCapacityLaw {
    /// Admit the complete typed Cp path and source pressure, preserving all
    /// breakpoints. Cp must be positive throughout every affine segment.
    /// Pressure, composition and all other source coordinates stay fixed.
    pub fn try_new(
        card: &MaterialCard,
        query: IsobaricHeatCapacityQuery<'_>,
        datum_temperature_k: f64,
    ) -> Result<Self, ThermalTransportError> {
        let invalid = ThermalTransportError::Invalid;
        let key = query.specific_heat_cp;
        if key.quantity() != QuantitySpec::dimensional(SPECIFIC_HEAT)
            || key.hardness_test().is_some()
            || key.elastic_component().is_some()
            || key.strain_component().is_some()
            || key.stress_component().is_some()
            || query.temperature_axis == query.pressure_axis
        {
            return Err(invalid(
                "isobaric calorimetry requires a scalar static specific Cp key",
            ));
        }
        for point in [query.lower, query.upper] {
            admit_temperature_point(point, query.temperature_axis)?;
            if point.axis_quantities().get(query.pressure_axis)
                != Some(&QuantitySpec::dimensional(Pressure::DIMS))
                || !point
                    .axes()
                    .get(query.pressure_axis)
                    .is_some_and(|&p| positive(p))
            {
                return Err(invalid(
                    "isobaric Cp requires explicit positive absolute pressure in pascals",
                ));
            }
        }
        if query.lower.axis_quantities() != query.upper.axis_quantities()
            || !query.lower.axes().keys().eq(query.upper.axes().keys())
            || query.lower.axes().iter().any(|(axis, value)| {
                axis != query.temperature_axis && query.upper.axes()[axis] != *value
            })
            || query.lower.axes()[query.temperature_axis]
                >= query.upper.axes()[query.temperature_axis]
        {
            return Err(invalid(
                "the Cp interval must retain its complete fixed source context",
            ));
        }
        let source = source_line(
            card,
            key,
            query.temperature_axis,
            query.lower,
            query.upper,
            query.selection,
        )?;
        let claim = card
            .claims()
            .claim(source.knots[0].1.receipt.selected)
            .ok_or(invalid("selected Cp source is missing"))?;
        if !claim.validity.bounds().contains_key(query.pressure_axis) {
            return Err(invalid(
                "selected Cp source must explicitly support the declared pressure",
            ));
        }
        let knots: Vec<_> = source
            .knots
            .iter()
            .map(|(t, answer)| [*t, answer.evidence.value.value])
            .collect();
        if knots.len() < 2
            || knots.len() > fs_matdb::MAX_PROPERTY_LINE_KNOTS
            || knots.iter().any(|[t, cp]| !positive(*t) || !positive(*cp))
            || knots.windows(2).any(|pair| pair[0][0] >= pair[1][0])
        {
            return Err(invalid(
                "source Cp requires bounded ordered temperatures and positive capacities",
            ));
        }
        let curve = CapacityCurve::new(&knots);
        curve.capacity(datum_temperature_k)?;
        let enthalpy_span = [
            curve.energy_change(datum_temperature_k, knots[0][0])?,
            curve.energy_change(datum_temperature_k, knots[knots.len() - 1][0])?,
        ];
        let minimum = knots.iter().map(|k| k[1]).fold(f64::INFINITY, f64::min);
        let maximum_slope = (1.0 / minimum).next_up();
        if !positive(maximum_slope)
            || !enthalpy_span.iter().all(|h| h.is_finite())
            || enthalpy_span[0] >= enthalpy_span[1]
        {
            return Err(invalid(
                "isobaric caloric range or inverse slope is unrepresentable",
            ));
        }
        let pressure_pa = query.lower.axes()[query.pressure_axis];
        let mut hash = DomainHasher::new("org.frankensim.fs-material.isobaric-cp-law.v1");
        hash.update(card.content_hash().as_bytes());
        for axis in [query.temperature_axis, query.pressure_axis] {
            hash.update(&(axis.len() as u64).to_le_bytes());
            hash.update(axis.as_bytes());
        }
        hash.update(&(knots.len() as u64).to_le_bytes());
        for (_, answer) in &source.knots {
            hash.update(answer.receipt.content_hash().as_bytes());
        }
        hash.update(&pressure_pa.to_bits().to_le_bytes());
        hash.update(&datum_temperature_k.to_bits().to_le_bytes());
        Ok(Self {
            material_card: card.content_hash(),
            endpoints: [query.lower.clone(), query.upper.clone()],
            source,
            knots,
            pressure_pa,
            datum_temperature_k,
            enthalpy_span,
            maximum_slope,
            identity: hash.finalize(),
        })
    }

    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    #[must_use]
    pub const fn material_card_identity(&self) -> ContentHash {
        self.material_card
    }
    /// Unmodified original answers, receipts and uncertainty at every source knot.
    #[must_use]
    pub const fn source_line(&self) -> &PropertyLineSlice {
        &self.source
    }
    #[must_use]
    pub const fn query_endpoints(&self) -> &[QueryPoint; 2] {
        &self.endpoints
    }
    #[must_use]
    pub const fn pressure_pa(&self) -> f64 {
        self.pressure_pa
    }
    #[must_use]
    pub const fn datum_temperature_k(&self) -> f64 {
        self.datum_temperature_k
    }
    #[must_use]
    pub fn temperature_span_k(&self) -> [f64; 2] {
        [self.knots[0][0], self.knots[self.knots.len() - 1][0]]
    }
    #[must_use]
    pub const fn specific_enthalpy_span_j_kg(&self) -> [f64; 2] {
        self.enthalpy_span
    }

    pub fn specific_heat_cp_j_kg_k(
        &self,
        temperature_k: f64,
    ) -> Result<f64, ThermalTransportError> {
        CapacityCurve::new(&self.knots).capacity(temperature_k)
    }

    /// Integral Cp dT from the selected datum; no density or volume enters.
    pub fn specific_enthalpy_j_kg(&self, temperature_k: f64) -> Result<f64, ThermalTransportError> {
        self.specific_enthalpy_change_j_kg(self.datum_temperature_k, temperature_k)
    }

    /// Direct crossed-segment integration retains small changes near a large datum.
    pub fn specific_enthalpy_change_j_kg(
        &self,
        initial_k: f64,
        final_k: f64,
    ) -> Result<f64, ThermalTransportError> {
        CapacityCurve::new(&self.knots).energy_change(initial_k, final_k)
    }

    /// Reversible integral Cp/T dT on this same isobar, not absolute source entropy.
    pub fn specific_entropy_change_at_temperatures_j_kg_k(
        &self,
        initial_k: f64,
        final_k: f64,
    ) -> Result<f64, ThermalTransportError> {
        CapacityCurve::new(&self.knots).entropy_change(initial_k, final_k)
    }

    /// Same analytic capacity inverse for heating or cooling, without extrapolation.
    pub fn temperature_after_specific_heat(
        &self,
        initial_k: f64,
        heat_j_kg: f64,
    ) -> Result<f64, ThermalTransportError> {
        CapacityCurve::new(&self.knots).temperature_after_energy_change(initial_k, heat_j_kg)
    }
}

impl ThermalLaw for IsobaricHeatCapacityLaw {
    fn storage_basis(&self) -> ThermalStorageBasis {
        ThermalStorageBasis::Enthalpy
    }

    fn temperature_k(&self, specific_h_j_kg: f64) -> Result<f64, DepositionError> {
        let [low, high] = self.enthalpy_span;
        if !specific_h_j_kg.is_finite() || specific_h_j_kg < low || specific_h_j_kg > high {
            return Err(DepositionError::MaterialDomain);
        }
        let [cold, hot] = self.temperature_span_k();
        if specific_h_j_kg == low {
            return Ok(cold);
        }
        if specific_h_j_kg == high {
            return Ok(hot);
        }
        self.temperature_after_specific_heat(self.datum_temperature_k, specific_h_j_kg)
            .map_err(|_| DepositionError::MaterialDomain)
    }

    fn maximum_temperature_slope(&self) -> f64 {
        self.maximum_slope
    }

    fn specific_entropy_change_j_kg_k(
        &self,
        initial_h_j_kg: f64,
        final_h_j_kg: f64,
    ) -> Result<Option<f64>, DepositionError> {
        let initial_k = self.temperature_k(initial_h_j_kg)?;
        let final_k = self.temperature_k(final_h_j_kg)?;
        if initial_h_j_kg == final_h_j_kg {
            return Ok(Some(0.0));
        }
        let dh = final_h_j_kg - initial_h_j_kg;
        let change = if initial_k == final_k {
            dh / initial_k
        } else {
            let ds = self
                .specific_entropy_change_at_temperatures_j_kg_k(initial_k, final_k)
                .map_err(|_| DepositionError::MaterialDomain)?;
            let represented_dh = self
                .specific_enthalpy_change_j_kg(initial_k, final_k)
                .map_err(|_| DepositionError::MaterialDomain)?;
            // Preserve actual stored dh when inversion rounds the endpoints.
            dh * (ds / represented_dh)
        };
        if !change.is_finite()
            || change == 0.0
            || change.is_sign_positive() != (final_h_j_kg > initial_h_j_kg)
        {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(Some(change))
    }
}
