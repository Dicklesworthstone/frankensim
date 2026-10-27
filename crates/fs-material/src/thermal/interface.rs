//! Ordered interface conditions from two original material/enthalpy states.
//! No interfacial temperature, pressure, phase fraction or mixture is inferred.
use super::source::{ThermalMaterialError, ThermalMaterialPoint};
use crate::phase::EquilibriumPhaseState;
use fs_matdb::{InterfaceSystemCard, MaterialCard, QueryPoint};
use fs_qty::{Dims, QuantitySpec, Temperature};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

/// Canonical ordered-interface resistance property, shared by spatial and
/// finite-cell thermal consumers [m² K/W].
pub const AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY: &str =
    "area-specific-thermal-contact-resistance";
/// SI dimensions of area-specific thermal resistance [m² K/W].
pub const AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS: Dims = Dims([0, -1, 3, 1, 0, 0]);

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
    /// Resolve a static resistance from the original ordered interface at the
    /// actual two caloric states. Temperature and phase remain distinct on the
    /// two sides; fixed pressure/medium/history conditions are never inferred.
    /// The returned bundle can refresh the same selected claim with
    /// `at_query_point`, without selecting a replacement at a support boundary.
    pub fn resolve_resistance(
        &self,
        states: [EquilibriumPhaseState; 2],
        selection: crate::state_point::MaterialPropertySelection,
    ) -> Result<crate::state_point::ResolvedInterfaceStatePoint, ThermalMaterialError> {
        use crate::state_point::{
            ScalarAdmissibility, ScalarPropertyRequirement, resolve_interface_state_point,
        };
        let point = self.at(states)?;
        let requirement = ScalarPropertyRequirement::try_new(
            AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
            AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .map_err(ThermalMaterialError::Material)?;
        let resolved =
            resolve_interface_state_point(self.interface, &point, &[requirement], selection)
                .map_err(ThermalMaterialError::Material)?;
        let property = &resolved.properties()[0];
        if !property
            .source_axes()
            .iter()
            .all(|(axis, dims)| self.is_response_axis(axis, *dims))
            || property
                .angular_frequency_band_rad_s()
                .map_err(ThermalMaterialError::Material)?
                .is_some_and(|(low, _)| low > 0.0)
            || point.axes().iter().any(|(axis, value)| {
                (axis == "omega"
                    || point
                        .axis_quantities()
                        .get(axis)
                        .and_then(|q| q.semantic_type())
                        .is_some_and(|s| matches!(s.kind(), QuantityKind::Frequency(_))))
                    && *value != 0.0
            })
        {
            return Err(ThermalMaterialError::Invalid(
                "contact resistance needs static response in owned thermal/phase or typed fixed external coordinates",
            ));
        }
        Ok(resolved)
    }

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

    /// Whether a source axis follows a declared evolving thermal coordinate.
    /// Fixed external response axes are admitted separately by `is_response_axis`.
    pub fn is_thermal_axis(&self, name: &str, dims: Dims) -> bool {
        (dims == Temperature::DIMS && (self.temperature_axes.contains(&name)
            || self.mean_temperature_axis == Some(name)))
            || (dims == Dims::NONE && self.liquid_fraction_axes.contains(&Some(name)))
    }

    /// Whether a static source response may vary along this coordinate.
    /// Declared thermal/phase coordinates follow the actual endpoint states.
    /// Typed external temperature, pressure, composition, moisture and
    /// concentration coordinates stay at their explicitly supplied values.
    /// Every axis of a curve or surface must pass this check; the original
    /// property evaluator still owns interpolation and supported-domain checks.
    /// Field, frequency, rate, geometry and history response axes require their
    /// own evolving owners and are not made static by appearing in the context.
    pub fn is_response_axis(&self, name: &str, dims: Dims) -> bool {
        if self.is_thermal_axis(name, dims) {
            return true;
        }
        let quantity = self.context.axis_quantities().get(name).copied();
        quantity.is_some_and(|quantity| quantity.dims() == dims)
            && crate::state_point::is_fixed_external_coordinate(
                name,
                quantity,
                self.context.axes().get(name).copied(),
            )
    }
}
