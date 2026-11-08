//! Laboratory energy transport bound to the original moving array's fluid/frame.

use super::{
    CompactSourceStrength, MovingCompactObserver, MovingMultipoleHistory, MovingMultipoleSource,
    MovingSourceArrayField, PreparedMovingSourceArray, poll,
};
use crate::AeroacError;
use crate::curle3d::energy::{self, AcousticEnergyTransport};

/// The unchanged coherent field plus its quadratic perturbation transport.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MovingSourceArrayEnergy {
    pub field: MovingSourceArrayField,
    pub transport: AcousticEnergyTransport,
}

impl PreparedMovingSourceArray<'_> {
    /// Full Q/F fields and laboratory energy at one original microphone event.
    /// Only the original receiver's accepted mean-flow frame supplies U; its
    /// translational velocity never becomes mean advection. Work, cancellation
    /// and causal/source-history admission remain the original field owner's.
    pub fn force_energy_at(
        &self,
        receiver_time_s: f64,
        source_at: &mut impl FnMut(usize, f64) -> Result<CompactSourceStrength, AeroacError>,
        force_history_at: &mut impl FnMut(
            usize,
            usize,
            &MovingCompactObserver,
            f64,
        ) -> Result<[f64; 3], AeroacError>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MovingSourceArrayEnergy, AeroacError> {
        let field = self.force_field_at(receiver_time_s, source_at, force_history_at, cancelled)?;
        self.with_energy(field, cancelled)
    }

    /// Full coherent volume/force/stress fields and transported energy.
    /// Every source and image contributes to p/u BEFORE quadratic reduction;
    /// source energy or intensity values are never summed independently.
    /// The complete particle velocity, including retained vortical history,
    /// contributes. Laboratory transport through a moving control surface
    /// separately needs its velocity; this result does not substitute the
    /// microphone motion for that flux or for the original fluid mean flow.
    pub fn multipole_energy_at(
        &self,
        receiver_time_s: f64,
        source_at: &mut impl FnMut(usize, f64) -> Result<MovingMultipoleSource, AeroacError>,
        history_at: &mut impl FnMut(
            usize,
            usize,
            &MovingCompactObserver,
            f64,
        ) -> Result<MovingMultipoleHistory, AeroacError>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MovingSourceArrayEnergy, AeroacError> {
        let field = self.multipole_field_at(receiver_time_s, source_at, history_at, cancelled)?;
        self.with_energy(field, cancelled)
    }

    fn with_energy(
        &self,
        field: MovingSourceArrayField,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MovingSourceArrayEnergy, AeroacError> {
        poll(cancelled)?;
        let receiver = self.receiver();
        let velocity = receiver
            .flow_frame()
            .map_or([0.0; 3], |frame| frame.velocity_m_s());
        let transport = energy::instantaneous(
            receiver.medium(),
            velocity,
            field.pressure_pa,
            field.particle_velocity_m_s,
        )?;
        poll(cancelled)?;
        Ok(MovingSourceArrayEnergy { field, transport })
    }
}
