//! Source Cv(T) on an explicitly fixed-density, fixed-composition path.
//!
//! At fixed mass and volume, du = Cv dT and ds = Cv/T dT. These integrals use
//! the same affine-capacity kernel as isobaric and volumetric storage, but the
//! timeline stores internal energy, not stream enthalpy. The original thermal
//! owner supplies heat exchange, time integration, entropy and energy balances.
//!
//! The caller selects true single-phase Cv data. The source must support the
//! actual density throughout the temperature interval; Cv sampled only along
//! an isobar is insufficient. This supplies neither an EOS/pressure trajectory
//! nor phase, reaction, strain or compositional evolution. In particular it
//! does not convert Cp or volumetric capacity to Cv by changing a key's name.
//! Thermodynamic convention: NASA Glenn, "Specific Heat – Cp & Cv",
//! https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-heat-cp-cv/.

use super::{
    DENSITY, SPECIFIC_HEAT, ThermalTransportError, admit_temperature_point,
    capacity::{CapacityCurve, CapacityInverse},
    positive, source_line,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{ClaimSelection, MaterialCard, PropertyKey, PropertyLineSlice, QueryPoint};
use fs_qty::{Pressure, QuantitySpec};
use fs_thermal::enthalpy::{DepositedCell, DepositionError, ThermalLaw, ThermalStorageBasis};

/// Original Cv source and complete endpoints of one fixed-density path.
#[derive(Clone, Copy, Debug)]
pub struct IsochoricHeatCapacityQuery<'a> {
    /// Mass-specific constant-volume heat capacity, J/(kg K).
    pub specific_heat_cv: &'a PropertyKey,
    /// Positive absolute temperature, in kelvin.
    pub temperature_axis: &'a str,
    /// Fixed positive mass density, kg/m³, explicitly supported by the source.
    pub density_axis: &'a str,
    /// Lower temperature with the complete fixed typed context.
    pub lower: &'a QueryPoint,
    /// Upper temperature with identical remaining typed context.
    pub upper: &'a QueryPoint,
    /// Original source selection, pinned over the complete path.
    pub selection: ClaimSelection,
}

/// Immutable source-backed internal energy with zero u at the declared datum T.
/// The original source receipts remain available; a numerical query creates
/// no new source authority. Native thermal cells must match the source density.
#[derive(Clone, Debug, PartialEq)]
pub struct IsochoricHeatCapacityLaw {
    material_card: ContentHash,
    endpoints: [QueryPoint; 2],
    temperature_axis: String,
    density_axis: String,
    source: PropertyLineSlice,
    knots: Vec<[f64; 2]>,
    density_kg_m3: f64,
    datum_temperature_k: f64,
    energy_span: [f64; 2],
    inverse: Option<CapacityInverse>,
    maximum_slope: f64,
    identity: ContentHash,
}

impl IsochoricHeatCapacityLaw {
    /// Prepare exact affine Cv integrals and a logarithmic-time inverse index.
    /// Every original source breakpoint is retained. A failed optional search
    /// index retains bounded direct inversion; it never changes the Cv law.
    ///
    /// Explicit pressure coordinates refuse: pressure generally changes along
    /// an isochore, and a fixed-pressure source slice cannot establish this path.
    pub fn try_new(
        card: &MaterialCard,
        query: IsochoricHeatCapacityQuery<'_>,
        datum_temperature_k: f64,
    ) -> Result<Self, ThermalTransportError> {
        let invalid = ThermalTransportError::Invalid;
        let key = query.specific_heat_cv;
        if key.quantity() != QuantitySpec::dimensional(SPECIFIC_HEAT)
            || key.hardness_test().is_some()
            || key.elastic_component().is_some()
            || key.strain_component().is_some()
            || key.stress_component().is_some()
            || query.temperature_axis == query.density_axis
        {
            return Err(invalid(
                "isochoric calorimetry requires a scalar static specific Cv key",
            ));
        }
        for point in [query.lower, query.upper] {
            admit_temperature_point(point, query.temperature_axis)?;
            if point.axis_quantities().get(query.density_axis)
                != Some(&QuantitySpec::dimensional(DENSITY))
                || !point.axes().get(query.density_axis).is_some_and(|&rho| positive(rho))
                || point.axis_quantities().values().any(|q| q.dims() == Pressure::DIMS)
            {
                return Err(invalid(
                    "isochoric Cv requires positive typed mass density without a fixed pressure/stress coordinate",
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
                "the Cv interval must retain its complete fixed source context",
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
        let claim = card.claims().claim(source.knots[0].1.receipt.selected)
            .ok_or(invalid("selected Cv source is missing"))?;
        if !claim.validity.bounds().contains_key(query.density_axis) {
            return Err(invalid(
                "selected Cv source must explicitly support the declared density",
            ));
        }
        let knots: Vec<_> = source.knots.iter()
            .map(|(t, answer)| [*t, answer.evidence.value.value])
            .collect();
        if knots.len() < 2
            || knots.len() > fs_matdb::MAX_PROPERTY_LINE_KNOTS
            || knots.iter().any(|[t, cv]| !positive(*t) || !positive(*cv))
            || knots.windows(2).any(|pair| pair[0][0] >= pair[1][0])
        {
            return Err(invalid(
                "source Cv requires ordered temperatures and positive capacities",
            ));
        }
        let curve = CapacityCurve::new(&knots);
        curve.capacity(datum_temperature_k)?;
        let energy_span = [
            curve.energy_change(datum_temperature_k, knots[0][0])?,
            curve.energy_change(datum_temperature_k, knots[knots.len() - 1][0])?,
        ];
        let minimum = knots.iter().map(|k| k[1]).fold(f64::INFINITY, f64::min);
        let maximum_slope = minimum.recip().next_up();
        if !positive(maximum_slope)
            || !energy_span.iter().all(|u| u.is_finite())
            || energy_span[0] >= energy_span[1]
        {
            return Err(invalid(
                "isochoric energy range or inverse slope is unrepresentable",
            ));
        }
        let inverse = CapacityInverse::new(&curve, datum_temperature_k).ok();
        let density_kg_m3 = query.lower.axes()[query.density_axis];
        let mut hash = DomainHasher::new("org.frankensim.fs-material.isochoric-cv-law.v2");
        hash.update(card.content_hash().as_bytes());
        for axis in [query.temperature_axis, query.density_axis] {
            hash.update(&(axis.len() as u64).to_le_bytes());
            hash.update(axis.as_bytes());
        }
        hash.update(&(knots.len() as u64).to_le_bytes());
        for (_, answer) in &source.knots {
            hash.update(answer.receipt.content_hash().as_bytes());
        }
        hash.update(&density_kg_m3.to_bits().to_le_bytes());
        hash.update(&datum_temperature_k.to_bits().to_le_bytes());
        Ok(Self {
            material_card: card.content_hash(),
            endpoints: [query.lower.clone(), query.upper.clone()],
            temperature_axis: query.temperature_axis.to_owned(),
            density_axis: query.density_axis.to_owned(),
            source,
            knots,
            density_kg_m3,
            datum_temperature_k,
            energy_span,
            inverse,
            maximum_slope,
            identity: hash.finalize(),
        })
    }

    /// Identity of this source path, density and energy datum.
    pub const fn identity(&self) -> ContentHash { self.identity }
    /// Original material card, with no replacement source selection.
    pub const fn material_card_identity(&self) -> ContentHash { self.material_card }
    /// Original answers, uncertainty and source receipts.
    pub const fn source_line(&self) -> &PropertyLineSlice { &self.source }
    /// Original complete typed source conditions.
    pub const fn query_endpoints(&self) -> &[QueryPoint; 2] { &self.endpoints }
    /// Original absolute-temperature coordinate.
    pub fn temperature_axis(&self) -> &str { &self.temperature_axis }
    /// Original fixed-density coordinate.
    pub fn density_axis(&self) -> &str { &self.density_axis }
    /// Source-supported fixed density [kg/m³].
    pub const fn density_kg_m3(&self) -> f64 { self.density_kg_m3 }
    /// Temperature at which the declared sensible internal energy is zero.
    pub const fn datum_temperature_k(&self) -> f64 { self.datum_temperature_k }
    /// Inclusive original temperature support [K].
    pub fn temperature_span_k(&self) -> [f64; 2] {
        [self.knots[0][0], self.knots[self.knots.len() - 1][0]]
    }
    /// Inclusive internal-energy support in the declared datum [J/kg].
    pub const fn specific_internal_energy_span_j_kg(&self) -> [f64; 2] { self.energy_span }

    /// Source-interpolated Cv at the actual temperature [J/(kg K)].
    pub fn specific_heat_cv_j_kg_k(&self, temperature_k: f64) -> Result<f64, ThermalTransportError> {
        CapacityCurve::new(&self.knots).capacity(temperature_k)
    }

    /// Integral Cv dT from the declared datum. No p/rho flow work is included.
    /// Reuse the caloric index and one original affine segment, retaining the
    /// canonical support endpoints and direct integration at numerical edges.
    pub fn specific_internal_energy_j_kg(
        &self, temperature_k: f64,
    ) -> Result<f64, ThermalTransportError> {
        let [cold, hot] = self.temperature_span_k();
        let [low, high] = self.energy_span;
        if temperature_k == cold { return Ok(low); }
        if temperature_k == hot { return Ok(high); }
        let curve = CapacityCurve::new(&self.knots);
        let direct = || curve.energy_change(self.datum_temperature_k, temperature_k);
        let energy = self.inverse.as_ref().map_or_else(
            direct,
            |inverse| inverse.energy_change_from_datum(&curve, temperature_k),
        )?;
        if energy < low || energy > high { return direct(); }
        Ok(energy)
    }

    /// Signed original Cv integral. Datum-aligned changes reuse the caloric
    /// index; other changes never subtract large reference energies.
    pub fn specific_internal_energy_change_j_kg(
        &self, initial_k: f64, final_k: f64,
    ) -> Result<f64, ThermalTransportError> {
        if initial_k == self.datum_temperature_k {
            return self.specific_internal_energy_j_kg(final_k);
        }
        if final_k == self.datum_temperature_k {
            return self.specific_internal_energy_j_kg(initial_k).map(|u| -u);
        }
        CapacityCurve::new(&self.knots).energy_change(initial_k, final_k)
    }

    /// Reversible entropy on this same fixed-density, fixed-composition path.
    pub fn specific_entropy_change_at_temperatures_j_kg_k(
        &self, initial_k: f64, final_k: f64,
    ) -> Result<f64, ThermalTransportError> {
        CapacityCurve::new(&self.knots).entropy_change(initial_k, final_k)
    }

    /// Analytic segment inversion of signed heat, without a frozen Cv update.
    pub fn temperature_after_specific_heat(
        &self, initial_k: f64, heat_j_kg: f64,
    ) -> Result<f64, ThermalTransportError> {
        CapacityCurve::new(&self.knots).temperature_after_energy_change(initial_k, heat_j_kg)
    }
}

impl ThermalLaw for IsochoricHeatCapacityLaw {
    fn storage_basis(&self) -> ThermalStorageBasis { ThermalStorageBasis::InternalEnergy }
    fn caloric_identity(&self) -> Option<[u8; 32]> { Some(*self.identity.as_bytes()) }

    fn validate_cell(&self, cell: &DepositedCell) -> Result<(), DepositionError> {
        let expected = self.density_kg_m3 * cell.geometry.volume_m3;
        if !positive(cell.mass_kg)
            || !positive(cell.geometry.volume_m3)
            || !positive(expected)
            || (cell.mass_kg - expected).abs()
                > 64.0 * f64::EPSILON * cell.mass_kg.max(expected)
        {
            return Err(DepositionError::InvalidInput(
                "cell mass/volume differs from its Cv source density",
            ));
        }
        self.temperature_k(cell.specific_h_j_kg)?;
        Ok(())
    }

    fn temperature_k(&self, u: f64) -> Result<f64, DepositionError> {
        let [low, high] = self.energy_span;
        if !u.is_finite() || u < low || u > high {
            return Err(DepositionError::MaterialDomain);
        }
        let [cold, hot] = self.temperature_span_k();
        if u == low { return Ok(cold); }
        if u == high { return Ok(hot); }
        self.inverse.as_ref().map_or_else(
            || self.temperature_after_specific_heat(self.datum_temperature_k, u),
            |inverse| inverse.temperature(&CapacityCurve::new(&self.knots), u),
        ).map_err(|_| DepositionError::MaterialDomain)
    }

    fn maximum_temperature_slope(&self) -> f64 { self.maximum_slope }

    fn specific_entropy_change_j_kg_k(
        &self, initial_u: f64, final_u: f64,
    ) -> Result<Option<f64>, DepositionError> {
        let initial_k = self.temperature_k(initial_u)?;
        let final_k = self.temperature_k(final_u)?;
        if initial_u == final_u { return Ok(Some(0.0)); }
        let du = final_u - initial_u;
        let entropy = if initial_k == final_k {
            du / initial_k
        } else {
            let ds = self.specific_entropy_change_at_temperatures_j_kg_k(initial_k, final_k)
                .map_err(|_| DepositionError::MaterialDomain)?;
            let represented_du = self.specific_internal_energy_change_j_kg(initial_k, final_k)
                .map_err(|_| DepositionError::MaterialDomain)?;
            du * (ds / represented_du)
        };
        if !entropy.is_finite()
            || entropy == 0.0
            || entropy.is_sign_positive() != (final_u > initial_u)
        {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(Some(entropy))
    }
}

#[cfg(test)]
#[path = "isochoric/tests.rs"]
mod tests;
