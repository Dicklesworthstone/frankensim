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
    /// Increasing original global mode indices and positive source rates [1/s].
    pub fn rates_per_s(&self) -> &[(usize, f64)] {
        &self.rates
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
                    || rate <= 0.0
                {
                    return Err(invalid(
                        "viscous heat needs ordered positive original component resistance rows",
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
}
