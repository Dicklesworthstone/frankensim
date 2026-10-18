//! Ordered interface conditions from two original material/enthalpy states.
//! No interfacial temperature, pressure, phase fraction or mixture is inferred.
use super::source::{ThermalMaterialError, ThermalMaterialPoint};
use crate::phase::EquilibriumPhaseState;
use fs_matdb::{InterfaceSystemCard, MaterialCard, QueryPoint};
use fs_qty::{Dims, QuantitySpec, Temperature};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

/// Original ordered surfaces, interface source and explicit evolving coordinates.
#[derive(Clone, Copy, Debug)]
pub struct ThermalInterfacePoint<'a> {
    /// Original ordered surface/medium/history system, not an unordered pair.
    pub interface: &'a InterfaceSystemCard,
    /// Actual bulk cards in interface A/B order. Their content identities must
    /// agree with the two phase states; named surface IDs must also agree.
    pub materials: [&'a MaterialCard; 2],
    /// Full interface conditions. Undeclared axes (pressure, gap, lot) stay fixed.
    pub context: &'a QueryPoint,
    /// Distinct existing absolute-temperature axes in A/B order.
    pub temperature_axes: [&'a str; 2],
    /// Optional distinct existing liquid MASS fraction axes, in A/B order.
    pub liquid_fraction_axes: [Option<&'a str>; 2],
    /// Optional existing Kelvin axis explicitly defined as (T_A + T_B)/2.
    /// This is a declared source coordinate, not thermal equilibrium or a
    /// solved interface temperature. None supplies no mean-temperature law.
    pub mean_temperature_axis: Option<&'a str>,
}

impl ThermalInterfacePoint<'_> {
    /// Rebind only the declared axes, using the SAME bulk-state coordinate owner.
    /// Both phase identities and ordered surface identities are checked before
    /// any property is queried. Aliased axes refuse instead of losing a side.
    pub fn at(&self, states: [EquilibriumPhaseState; 2]) -> Result<QueryPoint, ThermalMaterialError> {
        let invalid = ThermalMaterialError::Invalid;
        let axes: Vec<_> = self.temperature_axes.iter().copied()
            .chain(self.liquid_fraction_axes.iter().copied().flatten())
            .chain(self.mean_temperature_axis).collect();
        for (i, axis) in axes.iter().enumerate() {
            if axis.trim().is_empty() || axes[..i].contains(axis) {
                return Err(invalid("thermal interface coordinates must be nonblank and distinct"));
            }
        }
        let surfaces = [self.interface.surface_a(), self.interface.surface_b()];
        let mut point = self.context.clone();
        for side in 0..2 {
            if &surfaces[side].material != self.materials[side].id() {
                return Err(invalid("thermal interface must retain ordered bulk material identities"));
            }
            point = ThermalMaterialPoint {
                material: self.materials[side], context: &point,
                temperature_axis: self.temperature_axes[side],
                liquid_fraction_axis: self.liquid_fraction_axes[side],
            }.at(states[side])?;
        }
        if let Some(axis) = self.mean_temperature_axis {
            let absolute = QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature, ValueForm::Static));
            if !point.axes().contains_key(axis) || point.axis_quantities().get(axis)
                .is_some_and(|q| *q != absolute && *q != QuantitySpec::dimensional(Temperature::DIMS)) {
                return Err(invalid("declared interface mean needs an existing absolute-temperature axis"));
            }
            let mean = 0.5 * states[0].temperature_k() + 0.5 * states[1].temperature_k();
            point = point.with(axis, mean)
                .map_err(|e| ThermalMaterialError::Coordinate(e.to_string()))?;
        }
        Ok(point)
    }

    /// Whether a one-dimensional static source follows a declared evolving
    /// thermal coordinate. Other conditions remain explicit validity inputs;
    /// this does not admit frequency-domain or history-dependent response laws.
    pub fn is_thermal_axis(&self, name: &str, dims: Dims) -> bool {
        (dims == Temperature::DIMS && (self.temperature_axes.contains(&name)
            || self.mean_temperature_axis == Some(name)))
            || (dims == Dims::NONE && self.liquid_fraction_axes.contains(&Some(name)))
    }
}
