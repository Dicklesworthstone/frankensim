//! Reciprocal scalar thermoelectric transport for one constant-Seebeck span.
//!
//! The convention is J = sigma (-grad(V) - S grad(T)) and
//! q = -k grad(T) + S T J, at zero magnetic field. S is ABSOLUTE thermopower;
//! a coefficient measured relative to another conductor is not this property.
//! Electrical resistance, Fourier conduction and caloric storage retain their
//! existing owners. This module supplies only the reversible cross terms.

use fs_matdb::{EvaluationDecision, MaterialCard, QueryPoint};
use fs_qty::{Dims, QuantitySpec, Temperature};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

use crate::state_point::{MaterialPropertySelection, MaterialStatePointError,
    ResolvedMaterialStatePoint, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point};

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
/// Complete endpoint receipts remain available; a source curve is NEVER frozen
/// to create the constant-S model. Its Thomson coefficient T dS/dT is zero.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedSeebeckSpan {
    material: [ResolvedMaterialStatePoint; 2],
    bulk_material: ResolvedMaterialStatePoint,
    temperature_k: [f64; 2],
    seebeck_v_k: f64,
}

/// Reversible terminal transfer, excluding Fourier flux and Joule dissipation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermoelectricTransfer {
    /// Voltage DROP from terminal a to b, S (T_b - T_a) [V].
    pub voltage_drop_v: f64,
    /// Signed current from a to b [A].
    pub current_a: f64,
    /// Heat INTO terminal a and b, respectively: [-S T_a I, S T_b I] [W].
    pub terminal_heat_w: [f64; 2],
    /// Work supplied TO the circuit, -I S (T_b - T_a) [W].
    pub electrical_supply_w: f64,
}

impl ResolvedSeebeckSpan {
    /// Resolve the SAME constant claim at both absolute-temperature endpoints
    /// AND at the bulk temperature in `context` (which can be hotter through
    /// Joule heating). All three actual material states must remain supported.
    /// Only the named temperature coordinate changes; all external coordinates
    /// and their quantity descriptors are retained. Rectangular source validity
    /// then covers the full fixed-context interval between these endpoints.
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
        if property.answer().receipt.decision != EvaluationDecision::ConstantWithinValidity {
            return Err(invalid("constant-Seebeck span requires a validity-wide scalar claim, not a sampled curve"));
        }
        let seebeck_v_k = property.value_si();
        let pinned = MaterialPropertySelection::PinnedByProperty(vec![
            (ABSOLUTE_SEEBECK_PROPERTY.to_owned(), property.answer().receipt.selected)]);
        let bulk_material = resolve_material_state_point(card, context,
            std::slice::from_ref(&requirement), pinned.clone())?;
        let b = resolve_material_state_point(card, &point(temperature_k[1])?, &[requirement], pinned)?;
        Ok(Self { material: [a, b], bulk_material, temperature_k, seebeck_v_k })
    }

    /// Exact original endpoint source bundles and their usage receipts.
    pub const fn material(&self) -> &[ResolvedMaterialStatePoint; 2] { &self.material }
    /// Same selected claim at the actual bulk/context temperature.
    pub const fn bulk_material(&self) -> &ResolvedMaterialStatePoint { &self.bulk_material }
    /// Temperatures in declared a-to-b terminal order.
    pub const fn temperature_k(&self) -> [f64; 2] { self.temperature_k }
    /// Signed absolute thermopower [V/K], never a material-name preset.
    pub const fn seebeck_v_k(&self) -> f64 { self.seebeck_v_k }

    /// Reciprocal Kelvin transfer at a current supplied by the electrical owner.
    /// Sum(heat) + electrical_supply = 0 and Sum(heat/T) = 0 analytically.
    /// The coupled solver must independently admit their numerical balances.
    pub fn transfer(&self, current_a: f64) -> Result<ThermoelectricTransfer, ThermoelectricError> {
        let [a, b] = self.temperature_k;
        let voltage_drop_v = self.seebeck_v_k * (b - a);
        let terminal_heat_w = [-self.seebeck_v_k * a * current_a,
            self.seebeck_v_k * b * current_a];
        let electrical_supply_w = -voltage_drop_v * current_a;
        if [current_a, voltage_drop_v, terminal_heat_w[0], terminal_heat_w[1], electrical_supply_w]
            .iter().any(|x| !x.is_finite())
            || self.seebeck_v_k != 0.0 && a != b && voltage_drop_v == 0.0
            || self.seebeck_v_k != 0.0 && current_a != 0.0 && terminal_heat_w.contains(&0.0)
            || voltage_drop_v != 0.0 && current_a != 0.0 && electrical_supply_w == 0.0 {
            return Err(ThermoelectricError::Invalid("thermoelectric transfer is not representable"));
        }
        Ok(ThermoelectricTransfer { voltage_drop_v, current_a, terminal_heat_w, electrical_supply_w })
    }
}
