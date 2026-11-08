//! Original passive material/joint channels on the same resolved moving wall.

use super::{ThermoviscousInterior, ThermoviscousSourceError, invalid};
use crate::plate_viscosity::ModalPlateViscosity;
use fs_exec::CancelGate;
use fs_feec::acoustics::thermoviscous::ThermoviscousError;

#[cfg(test)]
mod tests;

impl<'c, S> ThermoviscousInterior<'c, S> {
    /// Attach the original full mass-normalized material/joint viscosity to
    /// this resolved fluid's existing structural wall basis. The passive
    /// factor is reused, including signed cross terms and exact loss null
    /// spaces; it is never replaced by one damping ratio per mode.
    ///
    /// Source shapes and mode order must match the viscosity's original basis.
    /// Existing wall-mode damping must be zero to avoid counting it twice.
    /// Resolved fluid viscosity, caloric relaxation and wall conduction retain
    /// their separate owners. The attachment applies to harmonic response,
    /// causal motion, original heat accounts and existing duct connections.
    ///
    /// The native cache consumes sqrt(R) times each original passive channel,
    /// so arbitrary factorization reference resistance cannot change C or the
    /// physical response. No source property is refitted or queried per step.
    pub fn with_structural_viscosity(
        mut self,
        viscosity: &ModalPlateViscosity,
        gate: &CancelGate,
    ) -> Result<Self, ThermoviscousSourceError> {
        if gate.is_requested() {
            return Err(ThermoviscousError::Cancelled.into());
        }
        let count = self
            .fe
            .moving_wall()
            .map(|wall| wall.structure.modes().len());
        if count != Some(viscosity.mode_count()) {
            return Err(invalid(
                "structural viscosity requires the complete original moving-wall basis",
            ));
        }
        let root_r = viscosity.reference_resistance_n_s_m().sqrt();
        let mut channels = Vec::with_capacity(viscosity.channel_count());
        for original in viscosity.channel_shapes() {
            if gate.is_requested() {
                return Err(ThermoviscousError::Cancelled.into());
            }
            channels.push(original.iter().map(|&shape| root_r * shape).collect());
        }
        self.fe = self
            .fe
            .with_moving_wall_damping_channels(channels, &mut || gate.is_requested())?;
        Ok(self)
    }
}
