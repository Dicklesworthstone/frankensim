//! Read the same mixture face stresses without advancing species or energy.
use super::*;

impl<L: TransportLaw> MolecularTransport3d<L> {
    pub(super) fn prepare_stress_field(
        &mut self,
        flow: &MixtureEuler3d,
        controls: DissipationControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), DissipationError> {
        poll(cancelled)?;
        controls.validate()?;
        if self.gas() != flow.gas() || self.grid != flow.grid || self.boundaries != flow.boundaries {
            return Err(DissipationError::Invalid("stress observer and accepted mixture differ"));
        }
        reservoir::admit(self.boundaries, &flow.reservoirs, self.pressure_outlets)?;
        spatial(
            &mut self.source, &self.grid, self.boundaries, self.pressure_outlets,
            flow.reservoirs.accepted(), &self.walls, &flow.states,
            &mut self.work.primitive, &mut self.work.gradient,
            &mut self.work.viscosity, &mut self.work.flux,
            #[cfg(feature = "experimental-sgs")]
            self.sgs.as_mut(),
            controls, cancelled,
        )?;
        poll(cancelled)
    }
}
