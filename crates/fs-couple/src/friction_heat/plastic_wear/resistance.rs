//! Original contact resistance at solved interface temperature and loaded pressure.
//!
//! The two uniform thermal receivers use one massless heat-generation plane.
//! A declared split of total area-specific resistance locates that plane; heat
//! partition still follows the solved temperatures and each body's finite
//! inventory. Pressure means force divided by ACTUAL loaded area on the native
//! equal-area grid, not nominal whole-face pressure or an inferred real area.
//! A source calibrated to another pressure/temperature convention is unsuitable.
//! This does not resolve lateral constriction or a finite third-body heat store.

use fs_material::state_point::{MaterialStatePointError, ResolvedInterfaceStatePoint};
use fs_material::thermal::interface::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
};
use fs_math::dd::Dd;
use fs_qty::QuantitySpec;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

use super::*;

/// Original selected resistance claim and its explicitly owned pressure axis.
#[derive(Clone, Debug)]
pub struct PlasticThermalResistance {
    source: ResolvedInterfaceStatePoint,
    pressure_axis: String,
    temperature_axis: String,
    fraction_to_first: f64,
}

/// Exact resistance source selected at the accepted coupled endpoint.
#[derive(Clone, Debug, PartialEq)]
pub struct PlasticThermalResistanceSample {
    pub source: ResolvedInterfaceStatePoint,
    /// Native resultant divided by actual positive-pressure grid area [Pa].
    pub mean_loaded_pressure_pa: f64,
    /// Resistance from the massless heat-generation plane to ordered A/B [m² K/W].
    pub area_resistance_m2_k_w: [f64; 2],
}

impl<E> From<MaterialStatePointError> for PlasticSlidingHeatError<E> {
    fn from(error: MaterialStatePointError) -> Self {
        Self::ResistanceSource(error)
    }
}

impl PlasticThermalResistance {
    pub const fn original_source(&self) -> &ResolvedInterfaceStatePoint {
        &self.source
    }

    /// Fraction of total contact resistance between the generation plane and A.
    pub const fn fraction_to_first(&self) -> f64 {
        self.fraction_to_first
    }

    pub(super) fn sample<E>(
        &self,
        contact: &SourcePlasticSlidingWear,
        state: &SourcePlasticSlidingWearState,
        interface_temperature_k: f64,
        cx: &Cx<'_>,
    ) -> Result<PlasticThermalResistanceSample, PlasticSlidingHeatError<E>> {
        checkpoint(cx)?;
        let pressure = state
            .wear_state()
            .contact_state()
            .contact_state()
            .pressure_pa();
        let count = pressure.iter().filter(|p| **p > 0.0).count();
        if count == 0 || !interface_temperature_k.is_finite() || interface_temperature_k <= 0.0 {
            return Err(PlasticSlidingHeatError::Invalid(
                "contact resistance needs a loaded area and positive interface temperature",
            ));
        }
        // Equal physical cell areas cancel in force/loaded-area. Scale before
        // compensated summation so a representable mean cannot overflow a sum.
        let mut mean = Dd::from_f64(0.0);
        for (i, &p) in pressure.iter().enumerate() {
            if i % 64 == 0 {
                checkpoint(cx)?;
            }
            if !p.is_finite() || p < 0.0 {
                return Err(PlasticSlidingHeatError::Invalid(
                    "finite native compressive pressure required",
                ));
            }
            if p > 0.0 {
                mean = mean + Dd::from_f64(p / count as f64);
            }
        }
        let mean_loaded_pressure_pa = finite(mean.to_f64())?;
        if mean_loaded_pressure_pa <= 0.0 {
            return Err(DepositionError::Arithmetic.into());
        }
        let query = self
            .source
            .query()
            .with(&self.pressure_axis, mean_loaded_pressure_pa)
            .and_then(|query| query.with(&self.temperature_axis, interface_temperature_k))
            .map_err(|_| {
                PlasticSlidingHeatError::Invalid("finite original resistance coordinates required")
            })?;
        let source = self.source.at_query_point(
            contact
                .wear_source()
                .contact()
                .normal_source()
                .elastic()
                .interface()
                .card(),
            &query,
        )?;
        let resistance = source
            .property_by_key(&resistance_key())
            .ok_or(PlasticSlidingHeatError::Invalid(
                "original thermal resistance property disappeared",
            ))?
            .value_si();
        if !resistance.is_finite() || resistance < 0.0 {
            return Err(PlasticSlidingHeatError::Invalid(
                "finite nonnegative sourced contact resistance required",
            ));
        }
        let first = resistance * self.fraction_to_first;
        let area_resistance_m2_k_w = [first, resistance - first];
        if (self.fraction_to_first > 0.0 && resistance > 0.0 && first == 0.0)
            || (self.fraction_to_first < 1.0
                && resistance > 0.0
                && area_resistance_m2_k_w[1] == 0.0)
        {
            return Err(DepositionError::Arithmetic.into());
        }
        checkpoint(cx)?;
        Ok(PlasticThermalResistanceSample {
            source,
            mean_loaded_pressure_pa,
            area_resistance_m2_k_w,
        })
    }
}

fn resistance_key() -> PropertyKey {
    PropertyKey::new(
        AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
        AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
    )
}

impl PlasticSlidingHeatContact<'_> {
    /// Bind original R(T_interface, mean_loaded_pressure, fixed conditions).
    ///
    /// `source` contains the canonical area-specific thermal contact resistance
    /// from the SAME ordered interface card. Pressure must be a distinct typed
    /// static-pressure axis outside the normal query. All fixed conditions and
    /// temperature roles must agree with the existing contact. Every endpoint
    /// pins the original claim; a support failure returns no physical candidate.
    ///
    /// `fraction_to_first` in [0,1] explicitly locates the heat-generation plane
    /// within total resistance. It is not an inferred heat fraction. The current
    /// route must already carry the same source-resolved split. This setting
    /// survives accepted thermal steps, adaptive refinement and material export.
    pub fn with_sourced_contact_resistance(
        mut self,
        source: ResolvedInterfaceStatePoint,
        pressure_axis: String,
        fraction_to_first: f64,
        cx: &Cx<'_>,
    ) -> Result<Self, PlasticSlidingHeatError> {
        checkpoint(cx)?;
        let normal = self.source.wear_source().contact().normal_source();
        let fixed = normal.interface_state();
        let temperature_axis = self.source.friction_source().axes().temperature.clone();
        let pressure_quantity =
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static));
        if !fraction_to_first.is_finite()
            || !(0.0..=1.0).contains(&fraction_to_first)
            || pressure_axis.is_empty()
            || fixed.state_coordinate(&pressure_axis).is_some()
            || source.axis_quantities().get(&pressure_axis) != Some(&pressure_quantity)
            || source.card_identity() != fixed.card_identity()
            || source.surface_materials() != fixed.surface_materials()
            || !source
                .query_point()
                .iter()
                .filter(|(axis, _)| axis != &pressure_axis)
                .eq(fixed.query_point().iter())
            || !source
                .axis_quantities()
                .iter()
                .filter(|(axis, _)| *axis != &pressure_axis)
                .eq(fixed.axis_quantities().iter())
            || source.state_coordinate(&temperature_axis)
                != Some(self.source.interface_temperature_k())
        {
            return Err(PlasticSlidingHeatError::Invalid(
                "resistance requires the original ordered interface/fixed conditions, owned temperature/pressure, and an explicit plane split",
            ));
        }
        let mut owned = vec![(pressure_axis.as_str(), pressure_quantity)];
        if let Some(&quantity) = source.axis_quantities().get(&temperature_axis) {
            owned.push((temperature_axis.as_str(), quantity));
        }
        if source.properties().len() != 1
            || source.property_by_key(&resistance_key()).is_none()
            || !source.properties()[0]
                .has_only_owned_coordinates_and_fixed_external_conditions(&owned)
        {
            return Err(PlasticSlidingHeatError::Invalid(
                "one original canonical thermal resistance claim with only owned or fixed external coordinates required",
            ));
        }
        let resistance = PlasticThermalResistance {
            source,
            pressure_axis,
            temperature_axis,
            fraction_to_first,
        };
        let sample = resistance.sample(
            &self.source,
            &self.state,
            self.source.interface_temperature_k(),
            cx,
        )?;
        let accepted = self.route.area_resistance_m2_k_w();
        if (0..2).any(|i| {
            (accepted[i] - sample.area_resistance_m2_k_w[i]).abs()
                > 64.0 * f64::EPSILON * accepted[i].max(sample.area_resistance_m2_k_w[i])
        }) {
            return Err(PlasticSlidingHeatError::Invalid(
                "accepted heat route must retain the original source-resolved resistance split",
            ));
        }
        self.resistance = Some(resistance);
        Ok(self)
    }

    pub const fn contact_resistance(&self) -> Option<&PlasticThermalResistance> {
        self.resistance.as_ref()
    }
}
