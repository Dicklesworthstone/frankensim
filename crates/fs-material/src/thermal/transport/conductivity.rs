//! Shared scalar source preparation; caloric constraints remain with the owner.

use super::{CONDUCTIVITY, ThermalTransportError, positive, source_line};
use fs_blake3::DomainHasher;
use fs_matdb::{ClaimSelection, MaterialCard, PropertyKey, PropertyLineSlice, QueryPoint};
use fs_qty::{
    QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermal::enthalpy::conduction::implicit::ImplicitConductionState;

#[derive(Clone, Debug)]
pub(super) struct PreparedScalarConductivity {
    source: PropertyLineSlice,
    knots: Vec<[f64; 2]>,
    slopes: Vec<f64>,
}

impl PreparedScalarConductivity {
    // The wrapper supplies its original caloric endpoints and separately
    // requires support for the thermodynamic constraint (pressure or density).
    pub(super) fn try_new(
        card: &MaterialCard,
        conductivity: &PropertyKey,
        temperature_axis: &str,
        lower: &QueryPoint,
        upper: &QueryPoint,
        selection: ClaimSelection,
    ) -> Result<Self, ThermalTransportError> {
        let invalid = ThermalTransportError::Invalid;
        let quantity = conductivity.quantity();
        let static_conductivity = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::ThermalConductivity,
            ValueForm::Static,
        ));
        if (quantity != QuantitySpec::dimensional(CONDUCTIVITY) && quantity != static_conductivity)
            || conductivity.hardness_test().is_some()
            || conductivity.elastic_component().is_some()
            || conductivity.strain_component().is_some()
            || conductivity.stress_component().is_some()
        {
            return Err(invalid(
                "caloric conduction requires scalar static conductivity",
            ));
        }
        let source = source_line(
            card,
            conductivity,
            temperature_axis,
            lower,
            upper,
            selection,
        )?;
        if source.knots.len() < 2 {
            return Err(invalid(
                "conductivity source requires a nonempty temperature interval",
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
                return Err(invalid("source conductivity derivative is unrepresentable"));
            }
            slopes.push(slope);
        }
        Ok(Self {
            source,
            knots,
            slopes,
        })
    }

    pub(super) const fn source_line(&self) -> &PropertyLineSlice {
        &self.source
    }

    // Optional integrated-layer preparation retains every original k knot.
    // Its separate bounded kernel never narrows admission of the point law.
    pub(super) fn integrated_curve(
        &self,
    ) -> Result<crate::thermal::conduction::ConductivityCurve, fs_thermal::enthalpy::DepositionError>
    {
        if self.knots.len() > crate::thermal::conduction::MAX_CONDUCTIVITY_KNOTS {
            return Err(fs_thermal::enthalpy::DepositionError::Budget(
                "integrated conductivity source knots",
            ));
        }
        crate::thermal::conduction::ConductivityCurve::new(self.knots.clone())
    }

    // Preserve the existing isobaric identity preimage: original knot count
    // followed by each original receipt, without a new wrapper hash or tag.
    pub(super) fn append_identity(&self, hash: &mut DomainHasher) {
        hash.update(&(self.source.knots.len() as u64).to_le_bytes());
        for (_, answer) in &self.source.knots {
            hash.update(answer.receipt.content_hash().as_bytes());
        }
    }

    // Right derivative at interior source knots; left derivative at the upper
    // endpoint. No source query, allocation or capacity lookup during evaluation.
    pub(super) fn at_temperature(
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

    // The caller owns de/dT, choosing Cp for h or Cv for u on its exact chart.
    // Historical native h-labelled tangent fields carry the same selected e.
    pub(super) fn state(
        &self,
        temperature: f64,
        capacity: f64,
    ) -> Result<ImplicitConductionState, ThermalTransportError> {
        let invalid = ThermalTransportError::Invalid;
        if !positive(capacity) {
            return Err(invalid(
                "caloric conduction requires positive finite capacity",
            ));
        }
        let (conductivity, dk_dt) = self.at_temperature(temperature)?;
        let temperature_slope = 1.0 / capacity;
        let conductivity_slope = dk_dt / capacity;
        if !positive(temperature_slope)
            || !conductivity_slope.is_finite()
            || (dk_dt != 0.0 && conductivity_slope == 0.0)
        {
            return Err(invalid(
                "caloric conduction state or tangent is unrepresentable",
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
