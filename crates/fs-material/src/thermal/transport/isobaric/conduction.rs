//! Source-exact isobaric storage and conduction on the existing thermal owner.
//!
//! The borrowed caloric law integrates its original affine Cp(T) exactly. This
//! adapter compiles conductivity on the same material card, temperature span,
//! pressure and fixed composition. It supplies dT/dh=1/Cp and dk/dh=(dk/dT)/Cp
//! to the existing implicit chain/graph and their adaptive time integrators.
//! No resampled enthalpy chart, effective capacity or second solver is created.
//!
//! Conductivity is a positive scalar for the caller's declared heat-flow
//! direction. Fixed reference geometry, single phase and the original isobaric
//! interpretation still apply. Spatial resolution and temporal accuracy are
//! properties of the thermal solve, not of this constitutive adapter.

use super::super::{CONDUCTIVITY, ThermalTransportError, positive, source_line};
use super::IsobaricHeatCapacityLaw;
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{ClaimSelection, MaterialCard, PropertyKey, PropertyLineSlice};
use fs_qty::{
    QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermal::enthalpy::conduction::ConductivityLaw;
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionLaw, ImplicitConductionObservables, ImplicitConductionState,
};
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};

/// One original caloric owner and its continuously supported conductivity path.
/// Construction resolves source evidence once; runtime queries allocate nothing.
#[derive(Clone, Debug)]
pub struct IsobaricConductionLaw<'a> {
    caloric: &'a IsobaricHeatCapacityLaw,
    source: PropertyLineSlice,
    knots: Vec<[f64; 2]>,
    slopes: Vec<f64>,
    identity: ContentHash,
}

impl<'a> IsobaricConductionLaw<'a> {
    /// Select conductivity from the exact original card and fixed source context.
    /// Every conductivity breakpoint is retained independently of the Cp knots.
    /// Scalar, explicitly directional static conductivity is admitted; tensor,
    /// dynamic, pressure-unbounded and out-of-support sources refuse.
    pub fn try_new(
        card: &MaterialCard,
        caloric: &'a IsobaricHeatCapacityLaw,
        conductivity: &PropertyKey,
        selection: ClaimSelection,
    ) -> Result<Self, ThermalTransportError> {
        let invalid = ThermalTransportError::Invalid;
        let quantity = conductivity.quantity();
        let static_conductivity = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::ThermalConductivity,
            ValueForm::Static,
        ));
        if card.content_hash() != caloric.material_card_identity() {
            return Err(invalid(
                "conduction must retain the original isobaric material card",
            ));
        }
        if (quantity != QuantitySpec::dimensional(CONDUCTIVITY)
            && quantity != static_conductivity)
            || conductivity.hardness_test().is_some()
            || conductivity.elastic_component().is_some()
            || conductivity.strain_component().is_some()
            || conductivity.stress_component().is_some()
        {
            return Err(invalid(
                "isobaric conduction requires scalar static conductivity",
            ));
        }
        let [lower, upper] = caloric.query_endpoints();
        let source = source_line(
            card,
            conductivity,
            caloric.temperature_axis(),
            lower,
            upper,
            selection,
        )?;
        if source.knots.len() < 2 {
            return Err(invalid("conductivity source requires a nonempty temperature interval"));
        }
        let claim = card
            .claims()
            .claim(source.knots[0].1.receipt.selected)
            .ok_or(invalid("selected conductivity source is missing"))?;
        if !claim.validity.bounds().contains_key(caloric.pressure_axis()) {
            return Err(invalid(
                "conductivity source must support the original isobaric pressure",
            ));
        }
        let knots: Vec<_> = source
            .knots
            .iter()
            .map(|(temperature, answer)| [*temperature, answer.evidence.value.value])
            .collect();
        if knots.iter().any(|knot| !positive(knot[1])) {
            return Err(invalid(
                "source conductivity must remain finite and positive",
            ));
        }
        let mut slopes = Vec::with_capacity(knots.len() - 1);
        for pair in knots.windows(2) {
            let width = pair[1][0] - pair[0][0];
            let difference = pair[1][1] - pair[0][1];
            if !positive(width) {
                return Err(invalid("source conductivity temperatures must increase"));
            }
            let slope = difference / width;
            if !slope.is_finite() || (difference != 0.0 && slope == 0.0) {
                return Err(invalid(
                    "source conductivity derivative is unrepresentable",
                ));
            }
            slopes.push(slope);
        }
        let mut hash = DomainHasher::new("org.frankensim.fs-material.isobaric-conduction.v1");
        hash.update(caloric.identity().as_bytes());
        hash.update(&(source.knots.len() as u64).to_le_bytes());
        for (_, answer) in &source.knots {
            hash.update(answer.receipt.content_hash().as_bytes());
        }
        Ok(Self {
            caloric,
            source,
            knots,
            slopes,
            identity: hash.finalize(),
        })
    }

    /// The exact caloric law registered on the thermal timeline.
    #[must_use]
    pub const fn caloric(&self) -> &'a IsobaricHeatCapacityLaw {
        self.caloric
    }

    /// Original material-source receipts and uncertainty, without resampling.
    #[must_use]
    pub const fn source_line(&self) -> &PropertyLineSlice {
        &self.source
    }

    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    // Share original k(T) interpolation between explicit and implicit consumers.
    // The explicit route needs no Cp query or enthalpy tangent evaluation.
    fn conductivity_at_temperature(
        &self,
        temperature: f64,
    ) -> Result<(f64, f64), ThermalTransportError> {
        let invalid = ThermalTransportError::Invalid;
        if !temperature.is_finite()
            || temperature < self.knots[0][0]
            || temperature > self.knots[self.knots.len() - 1][0]
        {
            return Err(invalid("temperature outside original conductivity source"));
        }
        let index = self
            .knots
            .partition_point(|knot| knot[0] <= temperature)
            .min(self.knots.len() - 1)
            - 1;
        let [a, ka] = self.knots[index];
        let [b, kb] = self.knots[index + 1];
        let conductivity = if temperature == a {
            ka
        } else if temperature == b {
            kb
        } else {
            let fraction = (temperature - a) / (b - a);
            (1.0 - fraction) * ka + fraction * kb
        };
        if !positive(conductivity) {
            return Err(invalid("source conductivity is unrepresentable"));
        }
        Ok((conductivity, self.slopes[index]))
    }

    /// One original inverse and analytic derivatives on the same temperature.
    /// Interior conductivity knots select their right derivative; the upper
    /// endpoint uses the left derivative. Cp remains continuous at its knots.
    pub fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitConductionState, ThermalTransportError> {
        let invalid = ThermalTransportError::Invalid;
        let temperature = self
            .caloric
            .temperature_k(h)
            .map_err(|_| invalid("enthalpy outside original isobaric law"))?;
        let cp = self.caloric.specific_heat_cp_j_kg_k(temperature)?;
        let (conductivity, dk_dt) = self.conductivity_at_temperature(temperature)?;
        let temperature_slope = 1.0 / cp;
        let conductivity_slope = dk_dt / cp;
        if !positive(temperature_slope)
            || !conductivity_slope.is_finite()
            || (dk_dt != 0.0 && conductivity_slope == 0.0)
        {
            return Err(invalid(
                "isobaric conduction state or tangent is unrepresentable",
            ));
        }
        Ok(ImplicitConductionState {
            temperature_k: temperature,
            temperature_slope_kg_k_j: temperature_slope,
            conductivity_w_m_k: conductivity,
            conductivity_slope_kg_w_m_k_j: conductivity_slope,
        })
    }
}

impl ConductivityLaw for IsobaricConductionLaw<'_> {
    fn conductivity_w_m_k(&self, h: f64) -> Result<f64, DepositionError> {
        let temperature = self.caloric.temperature_k(h)?;
        self.conductivity_at_temperature(temperature)
            .map(|(conductivity, _)| conductivity)
            .map_err(|_| DepositionError::MaterialDomain)
    }
}

impl ImplicitConductionLaw for IsobaricConductionLaw<'_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.caloric
    }

    fn conductivity_law(&self) -> &dyn ConductivityLaw {
        self
    }

    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        self.caloric.specific_enthalpy_span_j_kg()
    }

    fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        Self::state_at_specific_enthalpy(self, h).map_err(|_| DepositionError::MaterialDomain)
    }

    fn time_observables(&self, h: f64) -> Result<ImplicitConductionObservables, DepositionError> {
        Ok(ImplicitConductionObservables {
            temperature_k: self.caloric.temperature_k(h)?,
            phase_mass_fraction: None,
        })
    }
}

#[cfg(test)]
mod tests;
