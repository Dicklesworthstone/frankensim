//! Source-preserving local response of an ordered resistive interface.
//!
//! The field owner supplies the actual voltage jump and normal current density.
//! This adapter does not replace either with bulk fields or total face current.

use super::*;
use std::sync::Arc;

/// Local electrical coordinate, oriented from source surface A to surface B.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ElectricalInterfaceCoordinate {
    VoltageJump,
    VoltageJumpMagnitude,
    CurrentDensity,
    CurrentDensityMagnitude,
}
impl ElectricalInterfaceCoordinate {
    pub const fn quantity(self) -> QuantitySpec {
        QuantitySpec::dimensional(match self {
            Self::VoltageJump | Self::VoltageJumpMagnitude => Dims([2, 1, -3, 0, -1, 0]),
            Self::CurrentDensity | Self::CurrentDensityMagnitude => Dims([-2, 0, 0, 0, 1, 0]),
        })
    }
    fn value(self, voltage: f64, density: f64) -> f64 {
        match self {
            Self::VoltageJump => voltage,
            Self::VoltageJumpMagnitude => voltage.abs(),
            Self::CurrentDensity => density,
            Self::CurrentDensityMagnitude => density.abs(),
        }
    }
    fn admits(self, value: f64) -> bool {
        value.is_finite()
            && match self {
                Self::VoltageJumpMagnitude | Self::CurrentDensityMagnitude => value >= 0.0,
                _ => true,
            }
    }
}

/// One complete pinned interface source and its explicitly owned coordinates.
/// Temperature and declared side mass fractions are refreshed by other owners;
/// the electrical owner closes only local voltage jump and current density.
#[derive(Clone, Debug, PartialEq)]
pub struct ScalarElectricalInterfaceResponse {
    card: Arc<InterfaceSystemCard>,
    original: ResolvedElectricalInterface,
    axes: Vec<(String, ElectricalInterfaceCoordinate)>,
    phase_axes: [Option<String>; 2],
}

fn invalid(quantity: &'static str) -> MaterialStatePointError {
    MaterialStatePointError::InvalidDerived { quantity }
}

impl ScalarElectricalInterfaceResponse {
    pub fn new(
        card: &InterfaceSystemCard,
        original: ResolvedElectricalInterface,
        axes: &[(&str, ElectricalInterfaceCoordinate)],
    ) -> Result<Self, MaterialStatePointError> {
        Self::new_with_liquid_mass_fraction_axes(card, original, axes, [None, None])
    }

    /// Bind optional existing dimensional mass-fraction coordinates on A/B.
    /// This declares ownership, not a phase mixture or heat-partition law.
    pub fn new_with_liquid_mass_fraction_axes(
        card: &InterfaceSystemCard,
        original: ResolvedElectricalInterface,
        axes: &[(&str, ElectricalInterfaceCoordinate)],
        phase_axes: [Option<String>; 2],
    ) -> Result<Self, MaterialStatePointError> {
        if axes.len() > 2 || card.content_hash() != original.resolved().card_identity() {
            return Err(invalid(
                "original interface card and at most two local electrical coordinates required",
            ));
        }
        let law = Self {
            card: Arc::new(card.clone()),
            original,
            axes: axes
                .iter()
                .map(|(name, coordinate)| ((*name).into(), *coordinate))
                .collect(),
            phase_axes,
        };
        law.admit(&law.original)?;
        let replay = law.at_query_point(&law.original.resolved().query())?;
        if replay
            .resolved()
            .properties()
            .iter()
            .zip(law.original.resolved().properties())
            .any(|(a, b)| a.value_si().to_bits() != b.value_si().to_bits())
        {
            return Err(invalid(
                "interface response must replay its original selected sources",
            ));
        }
        Ok(law)
    }

    pub fn original(&self) -> &ResolvedElectricalInterface {
        &self.original
    }
    pub fn axes(&self) -> &[(String, ElectricalInterfaceCoordinate)] {
        &self.axes
    }
    pub fn liquid_mass_fraction_axes(&self) -> &[Option<String>; 2] {
        &self.phase_axes
    }

    /// Admit a refresh without changing source roles, pins, names or quantities.
    /// Non-electrical context may change only through this explicit source query.
    pub fn admit(
        &self,
        source: &ResolvedElectricalInterface,
    ) -> Result<(), MaterialStatePointError> {
        let old = self.original.resolved();
        let state = source.resolved();
        if source.surfaces() != self.original.surfaces()
            || state.card_identity() != self.card.content_hash()
            || state.axis_quantities() != old.axis_quantities()
            || state
                .query_point()
                .iter()
                .map(|(n, _)| n)
                .ne(old.query_point().iter().map(|(n, _)| n))
            || state.properties().len() != old.properties().len()
            || state
                .properties()
                .iter()
                .zip(old.properties())
                .any(|(a, b)| {
                    a.requirement() != b.requirement()
                        || a.source_axes() != b.source_axes()
                        || a.answer().receipt.selected != b.answer().receipt.selected
                })
            || !(source.resistance_ohm_m2().is_finite() && source.resistance_ohm_m2() > 0.0)
        {
            return Err(invalid(
                "local interface response requires complete original sources and positive resistance",
            ));
        }
        let mut owned = Vec::with_capacity(4);
        for (index, (name, coordinate)) in self.axes.iter().enumerate() {
            if name.trim().is_empty()
                || self.axes[..index]
                    .iter()
                    .any(|(n, c)| n == name || c.quantity() == coordinate.quantity())
                || state.axis_quantities().get(name) != Some(&coordinate.quantity())
                || !state
                    .state_coordinate(name)
                    .is_some_and(|v| coordinate.admits(v))
            {
                return Err(invalid(
                    "unique existing typed local interface electrical coordinates required",
                ));
            }
            owned.push((name.as_str(), coordinate.quantity()));
        }
        for name in self.phase_axes.iter().flatten() {
            let quantity = QuantitySpec::dimensional(Dims::NONE);
            if name.trim().is_empty()
                || owned.iter().any(|(n, _)| *n == name)
                || state.axis_quantities().get(name) != Some(&quantity)
                || !state
                    .state_coordinate(name)
                    .is_some_and(|v| v.is_finite() && (0.0..=1.0).contains(&v))
            {
                return Err(invalid(
                    "interface phase ownership requires unique dimensional mass fractions",
                ));
            }
            owned.push((name.as_str(), quantity));
        }
        if !state
            .property(ELECTRICAL_CONTACT_RESISTANCE_PROPERTY)
            .is_some_and(|p| p.has_only_owned_coordinates_and_fixed_external_conditions(&owned))
        {
            return Err(invalid(
                "every interface response and applicability coordinate needs an owner",
            ));
        }
        let electrical = |dims| {
            [
                Dims([2, 1, -3, 0, -1, 0]),
                Dims([-2, 0, 0, 0, 1, 0]),
                Dims([0, 0, 0, 0, 1, 0]),
                Dims([1, 1, -3, 0, -1, 0]),
            ]
            .contains(&dims)
        };
        let covered = |name: &str, dims| {
            !electrical(dims)
                || self
                    .axes
                    .iter()
                    .any(|(n, c)| n == name && c.quantity().dims() == dims)
        };
        if state
            .axis_quantities()
            .iter()
            .any(|(n, q)| !covered(n, q.dims()))
            || state
                .properties()
                .iter()
                .any(|p| p.source_axes().iter().any(|(n, d)| !covered(n, *d)))
        {
            return Err(invalid(
                "interface electrical coordinates cannot remain frozen or imply a bulk field",
            ));
        }
        Ok(())
    }

    pub fn at_query_point(
        &self,
        point: &QueryPoint,
    ) -> Result<ResolvedElectricalInterface, MaterialStatePointError> {
        let source = self.original.at_query_point(&self.card, point)?;
        self.admit(&source)?;
        Ok(source)
    }

    /// Query actual local A-to-B jump [V] and normal current density [A/m2].
    /// No area, bulk electric field or total contact current is substituted.
    pub fn at_operating_point(
        &self,
        seed: &ResolvedElectricalInterface,
        voltage_jump_v: f64,
        current_density_a_m2: f64,
    ) -> Result<ResolvedElectricalInterface, MaterialStatePointError> {
        self.admit(seed)?;
        if !voltage_jump_v.is_finite() || !current_density_a_m2.is_finite() {
            return Err(invalid(
                "finite actual interface jump and current density required",
            ));
        }
        let mut query = seed.resolved().query();
        for (name, coordinate) in &self.axes {
            query = query
                .with(name, coordinate.value(voltage_jump_v, current_density_a_m2))
                .map_err(|_| invalid("unrepresentable interface operating coordinate"))?;
        }
        self.at_query_point(&query)
    }

    /// Requery at damped electrical coordinates without blending resistance,
    /// source authority, temperature or phase. Original support remains binding.
    pub fn relax_operating_point(
        &self,
        applied: &ResolvedElectricalInterface,
        evaluated: &ResolvedElectricalInterface,
        weight: f64,
    ) -> Result<ResolvedElectricalInterface, MaterialStatePointError> {
        if !(weight.is_finite() && weight > 0.0 && weight <= 1.0) {
            return Err(invalid("interface relaxation must be in (0,1]"));
        }
        self.admit(applied)?;
        self.admit(evaluated)?;
        let mut query = applied.resolved().query();
        for ((name, a), (_, b)) in applied
            .resolved()
            .query_point()
            .iter()
            .zip(evaluated.resolved().query_point())
        {
            if self.axes.iter().any(|(n, _)| n == name) {
                let value = if weight == 1.0 {
                    *b
                } else if a.is_sign_negative() == b.is_sign_negative() {
                    weight.mul_add(b - a, *a)
                } else {
                    (1.0 - weight) * a + weight * b
                };
                query = query
                    .with(name, value)
                    .map_err(|_| invalid("unrepresentable relaxed interface coordinate"))?;
            } else if a.to_bits() != b.to_bits() {
                return Err(invalid(
                    "interface relaxation cannot change non-electrical context",
                ));
            }
        }
        if weight == 1.0 || self.axes.is_empty() {
            Ok(evaluated.clone())
        } else {
            self.at_query_point(&query)
        }
    }
}
