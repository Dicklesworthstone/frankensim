//! Reciprocal scalar thermoelectric transport on an original material span.
//!
//! The convention is J = sigma (-grad(V) - S grad(T)) and
//! q = -k grad(T) + S T J, at zero magnetic field. S is ABSOLUTE thermopower;
//! a coefficient measured relative to another conductor is not this property.
//! Electrical resistance, Fourier conduction and caloric storage retain their
//! existing owners. This module supplies only the reversible cross terms.

use fs_matdb::{MaterialCard, PropertyValue, QueryPoint};
use fs_qty::{Dims, QuantitySpec, Temperature};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

use crate::state_point::{MaterialPropertySelection, MaterialStatePointError,
    ResolvedMaterialStatePoint, ScalarAdmissibility, ScalarPropertyRequirement,
    integrate_piecewise_linear, pinned_temperature_line, resolve_material_state_point};

/// Source-bound Kelvin transport through the public constitutive graph.
pub mod graph;

/// Absolute scalar Seebeck coefficient under the convention in this module.
pub const ABSOLUTE_SEEBECK_PROPERTY: &str = "absolute_seebeck_coefficient";
/// Volts per kelvin, with signed values permitted.
pub const SEEBECK_DIMS: Dims = Dims([2, 1, -3, -1, -1, 0]);

/// Original material or unsupported constitutive-domain refusal.
#[derive(Debug)]
pub enum ThermoelectricError {
    /// Canonical source, units, selection or validity refusal.
    Material(MaterialStatePointError),
    /// Unsupported state, source form or unrepresentable arithmetic.
    Invalid(&'static str),
}
impl From<MaterialStatePointError> for ThermoelectricError {
    fn from(error: MaterialStatePointError) -> Self { Self::Material(error) }
}
impl core::fmt::Display for ThermoelectricError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "thermoelectric material: {self:?}")
    }
}
impl std::error::Error for ThermoelectricError {}

/// One source claim, admitted throughout the temperature span at fixed context.
/// Complete endpoint and bulk receipts remain available. Constant scalar or
/// supported piecewise-linear temperature paths supply exact integrals on the two
/// thermal edges a-to-bulk and bulk-to-b. Their integral means define a discrete
/// Kelvin pair: energy and entropy transfer close for all three thermal stores.
/// This is a spatial reduction; endpoint Peltier flux and distributed Thomson
/// heat approach the continuum law as the thermal cells are refined.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedSeebeckSpan {
    material: [ResolvedMaterialStatePoint; 2],
    bulk_material: ResolvedMaterialStatePoint,
    temperature_k: [f64; 2],
    bulk_temperature_k: f64,
    edge_seebeck_v_k: [f64; 2],
    voltage_drop_v: f64,
}

/// Reversible terminal transfer, excluding Fourier flux and Joule dissipation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermoelectricTransfer {
    /// Voltage DROP from terminal a to b, integral(S(T) dT) [V].
    pub voltage_drop_v: f64,
    /// Signed current from a to b [A].
    pub current_a: f64,
    /// Heat INTO terminal a and b: [-S_left T_a I, S_right T_b I] [W].
    pub terminal_heat_w: [f64; 2],
    /// Heat INTO the original bulk store, I T_bulk (S_left - S_right) [W].
    /// This signed variable-S exchange is separate from positive Joule loss.
    pub bulk_heat_w: f64,
    /// Work supplied TO the circuit, -I integral(S(T) dT) [W].
    pub electrical_supply_w: f64,
}

impl ResolvedSeebeckSpan {
    /// Resolve the SAME claim at both absolute-temperature endpoints
    /// AND at the bulk temperature in `context` (which can be hotter through
    /// Joule heating). All three actual material states must remain supported.
    /// Only the named temperature coordinate changes; all external coordinates
    /// and their quantity descriptors are retained. The canonical line query
    /// checks continuous source support through the entire path, including a
    /// fixed-context slice through a two-coordinate law. Every crossed cell
    /// breakpoint is integrated; unsupported gaps and discrete-only spans
    /// refuse even when the three endpoint queries individually succeed.
    pub fn resolve(
        card: &MaterialCard, context: &QueryPoint, temperature_axis: &str,
        temperature_k: [f64; 2], selection: MaterialPropertySelection,
    ) -> Result<Self, ThermoelectricError> {
        let invalid = ThermoelectricError::Invalid;
        let quantity = context.axis_quantities().get(temperature_axis);
        if !context.axes().contains_key(temperature_axis)
            || context.axes().get(temperature_axis).is_some_and(|t| *t <= 0.0)
            || quantity.is_some_and(|q| *q != QuantitySpec::dimensional(Temperature::DIMS)
                && *q != QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::AbsoluteTemperature, ValueForm::Static)))
            || temperature_k.iter().any(|t| !t.is_finite() || *t <= 0.0) {
            return Err(invalid("existing absolute-temperature coordinate and positive endpoints required"));
        }
        // This reduction has no magnetic-parity or frequency-dependent tensor.
        for (axis, value) in context.axes() {
            let dims = context.axis_quantities().get(axis).map(|q| q.dims());
            let field_or_rate = ["B", "H", "omega", "f"].contains(&axis.as_str())
                || dims.is_some_and(|d| [Dims([0, 1, -2, 0, -1, 0]),
                    Dims([-1, 0, 0, 0, 1, 0]), fs_qty::Frequency::DIMS].contains(&d));
            if field_or_rate && *value != 0.0 {
                return Err(invalid("scalar Kelvin coupling requires zero magnetic field and static source coordinates"));
            }
        }
        let requirement = ScalarPropertyRequirement::try_new(
            ABSOLUTE_SEEBECK_PROPERTY, SEEBECK_DIMS, ScalarAdmissibility::Finite)?;
        let point = |temperature| context.clone().with(temperature_axis, temperature)
            .map_err(|_| invalid("temperature update refused by the original query owner"));
        let a = resolve_material_state_point(card, &point(temperature_k[0])?,
            std::slice::from_ref(&requirement), selection)?;
        let property = a.property(ABSOLUTE_SEEBECK_PROPERTY)
            .ok_or(invalid("absolute Seebeck source is absent"))?;
        let selected = property.answer().receipt.selected;
        let claim = card.claims().claim(selected)
            .ok_or(invalid("selected absolute Seebeck claim is absent"))?;
        let pinned = MaterialPropertySelection::PinnedByProperty(vec![
            (ABSOLUTE_SEEBECK_PROPERTY.to_owned(), selected)]);
        let bulk_material = resolve_material_state_point(card, context,
            std::slice::from_ref(&requirement), pinned.clone())?;
        let b = resolve_material_state_point(card, &point(temperature_k[1])?, &[requirement], pinned)?;
        let bulk_temperature_k = context.axes()[temperature_axis];
        let lower = temperature_k[0].min(temperature_k[1]).min(bulk_temperature_k);
        let upper = temperature_k[0].max(temperature_k[1]).max(bulk_temperature_k);
        let path = pinned_temperature_line(card, &claim.key, temperature_axis,
            &point(lower)?, &point(upper)?, selected)?;
        let average = |start: f64, end: f64, at_start: f64| -> Result<f64, ThermoelectricError> {
            let value = if start == end { at_start } else {
                match &claim.value {
                    PropertyValue::Scalar { value, .. } => *value,
                    PropertyValue::Curve { abscissa, knots, .. } if abscissa == temperature_axis => {
                        let integral = integrate_piecewise_linear(knots, start, end)?;
                        let mean = integral / (end - start);
                        if integral != 0.0 && mean == 0.0 {
                            return Err(invalid("Seebeck integral mean is not representable"));
                        }
                        mean
                    }
                    PropertyValue::Curve { .. } | PropertyValue::Surface2D { .. } => {
                        let integral = integrate_piecewise_linear(&path, start, end)?;
                        let mean = integral / (end - start);
                        if integral != 0.0 && mean == 0.0 {
                            return Err(invalid("Seebeck integral mean is not representable"));
                        }
                        mean
                    }
                }
            };
            if !value.is_finite() || value != 0.0 && start != end && value * (end - start) == 0.0 {
                return Err(invalid("Seebeck integral mean is not representable"));
            }
            Ok(value)
        };
        let bulk_s = bulk_material.property(ABSOLUTE_SEEBECK_PROPERTY)
            .ok_or(invalid("bulk Seebeck source is absent"))?.value_si();
        let edge_seebeck_v_k = [
            average(temperature_k[0], bulk_temperature_k, property.value_si())?,
            average(bulk_temperature_k, temperature_k[1], bulk_s)?,
        ];
        let voltage_drop_v = edge_seebeck_v_k[0] * (bulk_temperature_k - temperature_k[0])
            + edge_seebeck_v_k[1] * (temperature_k[1] - bulk_temperature_k);
        if !voltage_drop_v.is_finite() {
            return Err(invalid("integrated Seebeck voltage is not representable"));
        }
        Ok(Self { material: [a, b], bulk_material, temperature_k, bulk_temperature_k,
            edge_seebeck_v_k, voltage_drop_v })
    }

    /// Exact original endpoint source bundles and their usage receipts.
    pub const fn material(&self) -> &[ResolvedMaterialStatePoint; 2] { &self.material }
    /// Same selected claim at the actual bulk/context temperature.
    pub const fn bulk_material(&self) -> &ResolvedMaterialStatePoint { &self.bulk_material }
    /// Temperatures in declared a-to-b terminal order.
    pub const fn temperature_k(&self) -> [f64; 2] { self.temperature_k }
    /// Actual temperature of the original bulk caloric store.
    pub const fn bulk_temperature_k(&self) -> f64 { self.bulk_temperature_k }
    /// Signed integral-mean thermopower [V/K] on a-to-bulk and bulk-to-b.
    /// An isothermal edge uses the original source at that temperature.
    pub const fn edge_seebeck_v_k(&self) -> [f64; 2] { self.edge_seebeck_v_k }

    /// Reciprocal Kelvin transfer at a current supplied by the electrical owner.
    /// Sum(heat) + electrical_supply = 0 and Sum(heat/T) = 0 analytically,
    /// including the bulk store. Freezing a variable S at one point, or omitting
    /// bulk_heat_w, breaks these relations. Constant S retains zero bulk heat.
    /// The coupled solver must independently admit their numerical balances.
    pub fn transfer(&self, current_a: f64) -> Result<ThermoelectricTransfer, ThermoelectricError> {
        let [a, b] = self.temperature_k;
        let [left, right] = self.edge_seebeck_v_k;
        let voltage_drop_v = self.voltage_drop_v;
        let terminal_heat_w = [-left * a * current_a, right * b * current_a];
        let bulk_heat_w = (left - right) * self.bulk_temperature_k * current_a;
        let electrical_supply_w = -voltage_drop_v * current_a;
        if [current_a, voltage_drop_v, terminal_heat_w[0], terminal_heat_w[1], bulk_heat_w, electrical_supply_w]
            .iter().any(|x| !x.is_finite())
            || current_a != 0.0 && (left != 0.0 && terminal_heat_w[0] == 0.0
                || right != 0.0 && terminal_heat_w[1] == 0.0
                || left != right && bulk_heat_w == 0.0)
            || voltage_drop_v != 0.0 && current_a != 0.0 && electrical_supply_w == 0.0 {
            return Err(ThermoelectricError::Invalid("thermoelectric transfer is not representable"));
        }
        Ok(ThermoelectricTransfer { voltage_drop_v, current_a, terminal_heat_w, bulk_heat_w, electrical_supply_w })
    }
}
