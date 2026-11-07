//! Fixed source roles shared by receiving dielectric/thermal owners.
//!
//! The complete original snapshot and every coefficient name remain visible.
//! Receiving heat owners must admit constancy over the actual temperature path;
//! this interface does not refresh a pole or manufacture parameter work.

use super::sealed;
use fs_material::{
    field_storage::dielectric::{bulk::ResolvedDebyeMaterial, lorentz::ResolvedLorentzMaterial},
    state_point::ResolvedMaterialStatePoint,
};

pub trait DielectricMaterialState: sealed::Sealed {
    fn material(&self) -> &ResolvedMaterialStatePoint;
    fn coefficient_properties(&self) -> impl Iterator<Item = &str>;
}

impl sealed::Sealed for ResolvedDebyeMaterial {}
impl DielectricMaterialState for ResolvedDebyeMaterial {
    fn material(&self) -> &ResolvedMaterialStatePoint {
        ResolvedDebyeMaterial::material(self)
    }
    fn coefficient_properties(&self) -> impl Iterator<Item = &str> {
        std::iter::once(self.high_frequency_property()).chain(
            self.pole_properties().iter().flat_map(|pole| {
                [
                    pole.permittivity_increment.as_str(),
                    pole.relaxation_time.as_str(),
                ]
            }),
        )
    }
}

impl sealed::Sealed for ResolvedLorentzMaterial {}
impl DielectricMaterialState for ResolvedLorentzMaterial {
    fn material(&self) -> &ResolvedMaterialStatePoint {
        ResolvedLorentzMaterial::material(self)
    }
    fn coefficient_properties(&self) -> impl Iterator<Item = &str> {
        std::iter::once(self.high_frequency_property()).chain(
            self.pole_properties().iter().flat_map(|pole| {
                [
                    pole.permittivity_increment.as_str(),
                    pole.angular_frequency.as_str(),
                    pole.damping_rate.as_str(),
                ]
            }),
        )
    }
}
