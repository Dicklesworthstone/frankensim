//! Exact source Cp branches on the existing finite-rate crystallization owner.
//!
//! At one fixed pressure and composition, dh/dT = Cp. A selected source value
//! L_ref = h_parent(T_ref) - h_crystal(T_ref) therefore gives
//! L(T) = L_ref + integral(Cp_parent - Cp_crystal) dT. This is the Kirchhoff
//! relation, following the IUPAC heat-capacity definition:
//! https://goldbook.iupac.org/terms/view/H02753.
//!
//! Both original Cp laws are integrated without a resampled enthalpy chart.
//! The parent has zero h at T_ref; the crystal has -L_ref there. The native
//! kinetic owner evolves their mass-weighted enthalpy and crystal inventory.
//! No additional latent-heat source is introduced. The supplied branches must
//! represent the same composition's parent and crystalline phases, including
//! any metastable support actually used. These roles are caller declarations.
//!
//! Cp and one enthalpy difference do not establish a relative phase entropy,
//! Gibbs stability, melting temperature, density or pressure-work law. Kinetic
//! applicability and the melting limit remain explicit on the original owner.

use super::super::{
    ThermalTransportError, admit_temperature_point,
    capacity::{CapacityCurve, CapacityInverse},
    positive,
};
use super::IsobaricHeatCapacityLaw;
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{
    ClaimSelection, MaterialAnswer, MaterialCard, PropertyKey, PropertySupportError, QueryPoint,
};
use fs_qty::{Dims, QuantitySpec};
use fs_thermal::crystallization::{
    CrystallizationCaloric, CrystallizationConfig, CrystallizationError, CrystallizationLaw,
    CrystallizationRate,
};
use fs_thermal::enthalpy::{DepositionError, ThermalLaw, ThermalStorageBasis};

const SPECIFIC_ENERGY: Dims = Dims([2, 0, -2, 0, 0, 0]);

/// Explicit original phase sources and one common-pressure enthalpy datum.
#[derive(Clone, Copy, Debug)]
pub struct IsobaricCrystallizationQuery<'a> {
    pub parent: &'a IsobaricHeatCapacityLaw,
    pub crystal: &'a IsobaricHeatCapacityLaw,
    /// Positive parent-minus-crystal enthalpy, J/kg, at `reference`.
    pub enthalpy_gap: &'a PropertyKey,
    /// Complete typed reference state supported by both Cp sources. It may
    /// lie outside the kinetic interval, e.g. a melting-point calorimetry datum.
    pub reference: &'a QueryPoint,
    pub selection: ClaimSelection,
    /// Runtime interval, continuously supported by BOTH original Cp laws.
    pub temperature_range_k: [f64; 2],
}

/// A borrowed original Cp law in the pair's common relative enthalpy chart.
///
/// Register this wrapper itself on a crystallizing thermal timeline. Its
/// identity includes the phase role, reference and admitted interval; the
/// unre-based Cp source is not the same caloric object.
#[derive(Clone, Debug)]
pub struct IsobaricCrystallizationBranch<'a> {
    source: &'a IsobaricHeatCapacityLaw,
    reference_k: f64,
    offset_j_kg: f64,
    range_k: [f64; 2],
    span_j_kg: [f64; 2],
    minimum_cp: f64,
    maximum_slope: f64,
    inverse: Option<CapacityInverse>,
    identity: ContentHash,
}

impl<'a> IsobaricCrystallizationBranch<'a> {
    fn new(
        source: &'a IsobaricHeatCapacityLaw,
        reference_k: f64,
        offset_j_kg: f64,
        range_k: [f64; 2],
        phase_role: u8,
    ) -> Result<Self, ThermalTransportError> {
        source.specific_heat_cp_j_kg_k(reference_k)?;
        let sensible = [
            source.specific_enthalpy_change_j_kg(reference_k, range_k[0])?,
            source.specific_enthalpy_change_j_kg(reference_k, range_k[1])?,
        ];
        let span_j_kg = sensible.map(|change| offset_j_kg + change);
        if !span_j_kg.iter().all(|h| h.is_finite())
            || span_j_kg[0] >= span_j_kg[1]
            || sensible
                .iter()
                .zip(span_j_kg)
                .any(|(&dh, h)| dh != 0.0 && h == offset_j_kg)
        {
            return Err(invalid("rebased phase enthalpy is not representable"));
        }
        let mut minimum_cp = source
            .specific_heat_cp_j_kg_k(range_k[0])?
            .min(source.specific_heat_cp_j_kg_k(range_k[1])?);
        for &[t, cp] in &source.knots {
            if range_k[0] < t && t < range_k[1] {
                minimum_cp = minimum_cp.min(cp);
            }
        }
        let maximum_slope = minimum_cp.recip().next_up();
        if !positive(maximum_slope) {
            return Err(invalid("phase inverse-capacity bound is not representable"));
        }
        let inverse = CapacityInverse::new(&CapacityCurve::new(&source.knots), reference_k).ok();
        let mut hash = DomainHasher::new("org.frankensim.fs-material.exact-cp-phase-branch.v1");
        hash.update(source.identity().as_bytes());
        hash.update(&[phase_role]);
        for value in [reference_k, offset_j_kg, range_k[0], range_k[1]] {
            hash.update(&value.to_bits().to_le_bytes());
        }
        Ok(Self {
            source,
            reference_k,
            offset_j_kg,
            range_k,
            span_j_kg,
            minimum_cp,
            maximum_slope,
            inverse,
            identity: hash.finalize(),
        })
    }

    /// Unmodified original source, including its receipts and uncertainty.
    pub const fn source(&self) -> &'a IsobaricHeatCapacityLaw {
        self.source
    }
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    pub const fn temperature_range_k(&self) -> [f64; 2] {
        self.range_k
    }
    pub const fn specific_enthalpy_span_j_kg(&self) -> [f64; 2] {
        self.span_j_kg
    }

    pub fn enthalpy_at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<f64, ThermalTransportError> {
        admit_temperature(temperature_k, self.range_k)?;
        if temperature_k == self.range_k[0] {
            return Ok(self.span_j_kg[0]);
        }
        if temperature_k == self.range_k[1] {
            return Ok(self.span_j_kg[1]);
        }
        let mut change = match &self.inverse {
            Some(inverse) => inverse
                .energy_change_from_datum(&CapacityCurve::new(&self.source.knots), temperature_k)?,
            None => self
                .source
                .specific_enthalpy_change_j_kg(self.reference_k, temperature_k)?,
        };
        let mut h = self.offset_j_kg + change;
        // The index is only an optimization. A final-bit accumulation
        // difference must not exclude a directly representable branch state.
        if !h.is_finite() || h < self.span_j_kg[0] || h > self.span_j_kg[1] {
            change = self
                .source
                .specific_enthalpy_change_j_kg(self.reference_k, temperature_k)?;
            h = self.offset_j_kg + change;
        }
        if !h.is_finite()
            || h < self.span_j_kg[0]
            || h > self.span_j_kg[1]
            || (change != 0.0 && h == self.offset_j_kg)
        {
            return Err(invalid("phase enthalpy increment is not representable"));
        }
        Ok(h)
    }
}

impl CrystallizationCaloric for IsobaricCrystallizationBranch<'_> {
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some(*self.identity.as_bytes())
    }
    fn minimum_heat_capacity_j_kg_k(&self) -> Option<f64> {
        Some(self.minimum_cp)
    }
    fn specific_enthalpy_j_kg(&self, temperature_k: f64) -> Result<f64, CrystallizationError> {
        self.enthalpy_at_temperature(temperature_k).map_err(|_| {
            CrystallizationError::Invalid("temperature outside the original exact-Cp phase branch")
        })
    }
    // The trait's default None deliberately leaves density/pressure work unavailable.
}

impl ThermalLaw for IsobaricCrystallizationBranch<'_> {
    fn storage_basis(&self) -> ThermalStorageBasis {
        ThermalStorageBasis::Enthalpy
    }
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some(*self.identity.as_bytes())
    }
    fn maximum_temperature_slope(&self) -> f64 {
        self.maximum_slope
    }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !h.is_finite() || h < self.span_j_kg[0] || h > self.span_j_kg[1] {
            return Err(DepositionError::MaterialDomain);
        }
        if h == self.span_j_kg[0] {
            return Ok(self.range_k[0]);
        }
        if h == self.span_j_kg[1] {
            return Ok(self.range_k[1]);
        }
        let sensible = h - self.offset_j_kg;
        if !sensible.is_finite() {
            return Err(DepositionError::MaterialDomain);
        }
        let t = self
            .inverse
            .as_ref()
            .map_or_else(
                || {
                    self.source
                        .temperature_after_specific_heat(self.reference_k, sensible)
                },
                |inverse| inverse.temperature(&CapacityCurve::new(&self.source.knots), sensible),
            )
            .map_err(|_| DepositionError::MaterialDomain)?;
        admit_temperature(t, self.range_k).map_err(|_| DepositionError::MaterialDomain)?;
        Ok(t)
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        initial_h: f64,
        final_h: f64,
    ) -> Result<Option<f64>, DepositionError> {
        let initial = self.temperature_k(initial_h)?;
        let final_t = self.temperature_k(final_h)?;
        if initial_h == final_h {
            return Ok(Some(0.0));
        }
        let dh = final_h - initial_h;
        let entropy = if initial == final_t {
            dh / initial
        } else {
            let ds = self
                .source
                .specific_entropy_change_at_temperatures_j_kg_k(initial, final_t)
                .map_err(|_| DepositionError::MaterialDomain)?;
            let represented_dh = self
                .source
                .specific_enthalpy_change_j_kg(initial, final_t)
                .map_err(|_| DepositionError::MaterialDomain)?;
            dh * (ds / represented_dh)
        };
        if !entropy.is_finite()
            || entropy == 0.0
            || entropy.is_sign_positive() != (final_h > initial_h)
        {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(Some(entropy))
    }
}

/// Source-bound parent/crystal calorics for the original kinetic and thermal owners.
#[derive(Clone, Debug)]
pub struct IsobaricCrystallizationPair<'a> {
    parent: IsobaricCrystallizationBranch<'a>,
    crystal: IsobaricCrystallizationBranch<'a>,
    gap_source: MaterialAnswer,
    identity: ContentHash,
}

impl<'a> IsobaricCrystallizationPair<'a> {
    /// Preserve both original Cp laws and source the relative phase datum once.
    /// Preparation checks all affine-Cp intervals, including interior minima
    /// of their quadratic enthalpy difference. Numerical checks use f64; this
    /// is not an interval enclosure or experimental stability certificate.
    pub fn try_new(
        card: &MaterialCard,
        query: IsobaricCrystallizationQuery<'a>,
    ) -> Result<Self, ThermalTransportError> {
        let parent = query.parent;
        let crystal = query.crystal;
        let axis = parent.temperature_axis();
        let pressure = parent.pressure_axis();
        let context = &parent.query_endpoints()[0];
        if card.content_hash() != parent.material_card_identity()
            || card.content_hash() != crystal.material_card_identity()
            || crystal.temperature_axis() != axis
            || crystal.pressure_axis() != pressure
            || !same_context(context, &crystal.query_endpoints()[0], axis)
            || !same_context(context, query.reference, axis)
        {
            return Err(invalid(
                "phase Cp sources and gap datum must retain one card and fixed typed context",
            ));
        }
        admit_temperature_point(query.reference, axis)?;
        let reference_k = query.reference.axes()[axis];
        let [low, high] = query.temperature_range_k;
        if !positive(low) || !positive(high) || low >= high {
            return Err(invalid(
                "crystallization requires an ordered positive temperature interval",
            ));
        }
        for source in [parent, crystal] {
            let span = source.temperature_span_k();
            if low < span[0] || high > span[1] || reference_k < span[0] || reference_k > span[1] {
                return Err(invalid(
                    "both Cp sources must cover the phase interval and reference temperature",
                ));
            }
        }
        let key = query.enthalpy_gap;
        if key.quantity() != QuantitySpec::dimensional(SPECIFIC_ENERGY)
            || key.hardness_test().is_some()
            || key.elastic_component().is_some()
            || key.strain_component().is_some()
            || key.stress_component().is_some()
        {
            return Err(invalid(
                "phase datum requires scalar static mass-specific enthalpy",
            ));
        }
        let gap_source = match query.selection {
            ClaimSelection::Policy(policy) => {
                card.claims().query_typed(key, query.reference, policy)
            }
            ClaimSelection::Pinned(pin) => {
                card.claims().query_pinned_typed(key, query.reference, pin)
            }
        }
        .map_err(|error| ThermalTransportError::Source {
            property: key.name().to_owned(),
            error: Box::new(PropertySupportError::Evaluation {
                point: query.reference.clone(),
                error: Box::new(error),
            }),
        })?;
        let selected = card
            .claims()
            .claim(gap_source.receipt.selected)
            .ok_or(invalid("selected phase enthalpy-gap source is missing"))?;
        if !selected.validity.bounds().contains_key(axis)
            || !selected.validity.bounds().contains_key(pressure)
            || !positive(gap_source.evidence.value.value)
        {
            return Err(invalid(
                "positive phase enthalpy gap requires explicit source temperature and pressure support",
            ));
        }
        let parent = IsobaricCrystallizationBranch::new(parent, reference_k, 0.0, [low, high], 0)?;
        let crystal = IsobaricCrystallizationBranch::new(
            crystal,
            reference_k,
            -gap_source.evidence.value.value,
            [low, high],
            1,
        )?;
        let mut hash =
            DomainHasher::new("org.frankensim.fs-material.exact-cp-crystallization-pair.v1");
        hash.update(parent.identity.as_bytes());
        hash.update(crystal.identity.as_bytes());
        hash.update(gap_source.receipt.content_hash().as_bytes());
        let pair = Self {
            parent,
            crystal,
            gap_source,
            identity: hash.finalize(),
        };
        pair.check_gap()?;
        Ok(pair)
    }

    pub const fn parent(&self) -> &IsobaricCrystallizationBranch<'a> {
        &self.parent
    }
    pub const fn crystal(&self) -> &IsobaricCrystallizationBranch<'a> {
        &self.crystal
    }
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    pub const fn enthalpy_gap_source(&self) -> &MaterialAnswer {
        &self.gap_source
    }
    pub const fn temperature_range_k(&self) -> [f64; 2] {
        self.parent.range_k
    }
    pub fn maximum_temperature_slope(&self) -> f64 {
        self.parent.maximum_slope.max(self.crystal.maximum_slope)
    }

    /// Parent-minus-crystal enthalpy at the actual temperature; no constant-L assumption.
    pub fn latent_enthalpy_j_kg(&self, temperature_k: f64) -> Result<f64, ThermalTransportError> {
        let gap = self.parent.enthalpy_at_temperature(temperature_k)?
            - self.crystal.enthalpy_at_temperature(temperature_k)?;
        if !positive(gap) {
            return Err(invalid(
                "phase enthalpy gap must remain positive and representable",
            ));
        }
        Ok(gap)
    }

    /// Bind these exact branches to the existing finite-rate kinetic integrator.
    /// Kinetic identity must describe the actual supplied rate and applicability;
    /// the caller supplies the melting limit, which Cp and L do not determine.
    pub fn law<'b>(
        &'b self,
        rate: &'b dyn CrystallizationRate,
        kinetic_source_identity: ContentHash,
        avrami_exponent: f64,
        maximum_crystalline_mass_fraction: f64,
        melting_temperature_k: f64,
    ) -> Result<CrystallizationLaw<'b>, CrystallizationError> {
        let mut hash =
            DomainHasher::new("org.frankensim.fs-material.exact-cp-crystallization-law.v1");
        hash.update(self.identity.as_bytes());
        hash.update(kinetic_source_identity.as_bytes());
        for value in [
            avrami_exponent,
            maximum_crystalline_mass_fraction,
            melting_temperature_k,
        ] {
            hash.update(&value.to_bits().to_le_bytes());
        }
        CrystallizationLaw::new(
            &self.parent,
            &self.crystal,
            rate,
            CrystallizationConfig {
                identity: *hash.finalize().as_bytes(),
                avrami_exponent,
                maximum_crystalline_mass_fraction,
                temperature_range_k: self.temperature_range_k(),
                melting_temperature_k,
            },
        )
    }

    fn check_gap(&self) -> Result<(), ThermalTransportError> {
        let [low, high] = self.temperature_range_k();
        let mut temperatures = vec![low, high];
        for source in [self.parent.source, self.crystal.source] {
            temperatures.extend(
                source
                    .knots
                    .iter()
                    .map(|k| k[0])
                    .filter(|t| low < *t && *t < high),
            );
        }
        temperatures.sort_by(f64::total_cmp);
        temperatures.dedup();
        for &t in &temperatures {
            self.latent_enthalpy_j_kg(t)?;
        }
        for interval in temperatures.windows(2) {
            let difference = |t| -> Result<f64, ThermalTransportError> {
                Ok(self.parent.source.specific_heat_cp_j_kg_k(t)?
                    - self.crystal.source.specific_heat_cp_j_kg_k(t)?)
            };
            let left = difference(interval[0])?;
            let right = difference(interval[1])?;
            if left < 0.0 && right > 0.0 {
                // L'(T) is affine. This zero is the only interior minimum.
                // Scaling avoids overflow from subtracting opposite large slopes.
                let scale = (-left).max(right);
                let weight = (-left / scale) / ((-left / scale) + right / scale);
                let t = interval[0] + weight * (interval[1] - interval[0]);
                if !(interval[0] < t && t < interval[1]) {
                    return Err(invalid(
                        "interior phase-enthalpy minimum is not representable",
                    ));
                }
                self.latent_enthalpy_j_kg(t)?;
            }
        }
        Ok(())
    }
}

fn same_context(a: &QueryPoint, b: &QueryPoint, axis: &str) -> bool {
    a.axis_quantities() == b.axis_quantities()
        && a.axes().keys().eq(b.axes().keys())
        && a.axes()
            .iter()
            .all(|(name, value)| name == axis || b.axes()[name] == *value)
}

fn admit_temperature(t: f64, range: [f64; 2]) -> Result<(), ThermalTransportError> {
    if !t.is_finite() || t < range[0] || t > range[1] {
        return Err(invalid("temperature outside original phase Cp support"));
    }
    Ok(())
}

fn invalid(reason: &'static str) -> ThermalTransportError {
    ThermalTransportError::Invalid(reason)
}

#[cfg(test)]
mod tests;
