//! Original source-region contributions to native diagonal modal resistance.
//! These are physical operator rows, never fractions of a measured total loss.
use super::*;
use std::collections::BTreeSet;

/// An immutable physical region's original mass-normalized viscous operator.
/// Geometry producers select one region within one original modal component.
/// A row `(mode, rate)` contributes `rate * velocity[mode]^2` watts, where the
/// native velocity is in m sqrt(kg)/s. All rows of a region share its ONE caloric
/// inventory. Material and ambient resistance remain independently accounted.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialViscousHeatSource {
    component: usize,
    region: usize,
    rates: Vec<(usize, f64)>,
    // Original property name, geometry/mass multiplier and accepted viscosity.
    // Kept private so a temperature binder cannot invent a source selector.
    properties: Vec<(String, f64, f64)>,
}

impl MaterialViscousHeatSource {
    pub(crate) fn from_material_region(
        component: usize,
        region: usize,
        rates: Vec<(usize, f64)>,
    ) -> Self {
        Self {
            component,
            region,
            rates,
            properties: Vec::new(),
        }
    }

    /// Original mechanical component, not a newly created material body.
    pub const fn component(&self) -> usize {
        self.component
    }
    /// Original physical region within that component (e.g. an annular layer).
    pub const fn region(&self) -> usize {
        self.region
    }
    /// Increasing original global mode indices and nonnegative source rates [1/s].
    /// A zero row retains a source law that may become positive as its body warms.
    pub fn rates_per_s(&self) -> &[(usize, f64)] {
        &self.rates
    }

    pub(crate) fn with_property_rows(mut self, properties: Vec<(String, f64, f64)>)
        -> Result<Self, ModalCouplingError> {
        if properties.len() != self.rates.len() || properties.iter().zip(&self.rates)
            .any(|((name, multiplier, value), (_, rate))| name.trim().is_empty()
                || !multiplier.is_finite() || *multiplier <= 0.0
                || !value.is_finite() || *value < 0.0
                || !(*multiplier * *value).is_finite()
                || (*rate - *multiplier * *value).abs()
                    > 128.0 * f64::EPSILON * rate.abs().max((*multiplier * *value).abs())) {
            return Err(invalid("viscosity selectors must retain each original nonnegative source coefficient and geometry multiplier"));
        }
        self.properties = properties;
        Ok(self)
    }

    pub(crate) fn property_names(&self) -> impl Iterator<Item = &str> {
        self.properties.iter().map(|(name, _, _)| name.as_str())
    }

    pub(crate) fn at_material_source(&self, source: &fs_material::state_point::ResolvedMaterialStatePoint)
        -> Result<Self, ModalCouplingError> {
        let mut next = self.clone();
        for ((name, multiplier, accepted), (_, rate)) in next.properties.iter_mut().zip(&mut next.rates) {
            let property = source.property(name).ok_or_else(|| invalid("viscosity mobility lost its original property selector"))?;
            let value = property.value_si();
            if property.requirement().quantity() != fs_qty::QuantitySpec::dimensional(fs_qty::DynViscosity::DIMS)
                || !property.is_constant_at_fixed_environment() || !value.is_finite() || value < 0.0 {
                return Err(invalid("viscosity mobility requires its nonnegative source law at held external conditions"));
            }
            if value.to_bits() != accepted.to_bits() {
                *rate = finite(value * *multiplier)?;
                if value > 0.0 && *rate == 0.0 { return Err(invalid("source viscosity lost its positive original modal resistance")); }
                *accepted = value;
            }
        }
        Ok(next)
    }
}

impl CoupledModalSystem {
    /// Check the exact union of routed regions against the ORIGINAL diagonal
    /// resistance. Region identity prevents routing the same layer twice even
    /// when other materials/gas contribute additional damping to its modes.
    /// Different physical regions may legitimately act on the same velocity.
    pub(crate) fn admit_material_viscosity<'a>(
        &self,
        sources: impl IntoIterator<Item = &'a MaterialViscousHeatSource>,
        gate: &CancelGate,
    ) -> Result<(), ModalCouplingError> {
        let mut owners = BTreeSet::new();
        let mut total = vec![0.0; self.mode_count()];
        for source in sources {
            poll(Some(gate))?;
            let component = self
                .models
                .get(source.component)
                .ok_or_else(|| invalid("viscous heat needs its original modal component"))?;
            if source.rates.is_empty() || !owners.insert((source.component, source.region)) {
                return Err(invalid(
                    "each original material viscosity region may heat its body only once",
                ));
            }
            let start = self.offsets[source.component];
            let end = start + component.modes().len();
            let mut previous = None;
            for &(mode, rate) in &source.rates {
                poll(Some(gate))?;
                if mode < start
                    || mode >= end
                    || previous.is_some_and(|p| mode <= p)
                    || !rate.is_finite()
                    || rate < 0.0
                {
                    return Err(invalid(
                        "viscous heat needs ordered nonnegative original component resistance rows",
                    ));
                }
                previous = Some(mode);
                total[mode] = finite(total[mode] + rate)?;
            }
        }
        for (rate, mode) in total.iter().zip(self.models.iter().flat_map(|m| m.modes())) {
            poll(Some(gate))?;
            let native = finite(2.0 * mode.damping_ratio * mode.angular_frequency_rad_s)?;
            let roundoff = 64.0 * f64::EPSILON * native.max(*rate);
            if *rate - native > roundoff {
                return Err(invalid(
                    "routed source viscosities exceed the original native modal resistance",
                ));
            }
        }
        Ok(())
    }

    /// Replace only explicitly retained source contributions. The remainder of
    /// each ORIGINAL modal operator (gas and other material regions) survives;
    /// no measured decay, total-loss fraction or reconstructed damping is used.
    pub(crate) fn material_viscosity_updates<'a>(
        &self, pairs: impl IntoIterator<Item = (&'a MaterialViscousHeatSource, &'a MaterialViscousHeatSource)>,
        gate: &CancelGate,
    ) -> Result<Vec<super::super::rates::mobility::ModalDampingRatio>, ModalCouplingError> {
        let mut old = vec![0.0; self.mode_count()];
        let mut new = vec![0.0; self.mode_count()];
        let mut changed = vec![false; self.mode_count()];
        let mut originals = Vec::new();
        for (before, after) in pairs {
            poll(Some(gate))?;
            if (before.component, before.region) != (after.component, after.region)
                || before.rates.len() != after.rates.len()
                || before.properties.len() != after.properties.len()
                || before.properties.iter().zip(&after.properties).any(|(a, b)| a.0 != b.0 || a.1.to_bits() != b.1.to_bits()) {
                return Err(invalid("viscosity continuation must retain the original source regions, selectors and geometry"));
            }
            originals.push(before);
            for (&(mode, a), &(index, b)) in before.rates.iter().zip(&after.rates) {
                if mode != index || mode >= old.len() || !b.is_finite() || b < 0.0 {
                    return Err(invalid("viscosity continuation lost an original nonnegative resistance row"));
                }
                old[mode] = finite(old[mode] + a)?;
                new[mode] = finite(new[mode] + b)?;
                changed[mode] |= a != b;
            }
        }
        self.admit_material_viscosity(originals, gate)?;
        let mut updates = Vec::new();
        for (index, mode) in self.models.iter().flat_map(|m| m.modes()).enumerate() {
            if !changed[index] { continue; }
            poll(Some(gate))?;
            let native = finite(2.0 * mode.damping_ratio * mode.angular_frequency_rad_s)?;
            let remainder = native - old[index];
            let tolerance = 64.0 * f64::EPSILON * native.max(old[index]);
            if mode.angular_frequency_rad_s <= 0.0 || remainder < -tolerance {
                return Err(invalid("source viscosity cannot remove unowned native modal resistance"));
            }
            // Only admission-sized floating cancellation can reach this zero.
            let next_rate = finite(remainder.max(0.0) + new[index])?;
            let next = finite(next_rate / (2.0 * mode.angular_frequency_rad_s))?;
            if next_rate > 0.0 && next == 0.0 { return Err(invalid("updated native viscosity is unrepresentable")); }
            updates.push(super::super::rates::mobility::ModalDampingRatio {
                mode: index, previous: mode.damping_ratio, next,
            });
        }
        Ok(updates)
    }
}
