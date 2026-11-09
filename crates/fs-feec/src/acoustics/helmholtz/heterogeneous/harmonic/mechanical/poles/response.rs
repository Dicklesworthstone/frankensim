//! Complete source-to-mode participation on the original acoustic FE owner.

use super::{
    MechanicalPoleError, MechanicalPoleLimits, MechanicalResonancePole,
    PreparedMechanicalHarmonic,
};
use super::super::{
    BoundaryVelocitySource, HeterogeneousError, PreparedHeterogeneous,
    PreparedMechanicalHarmonicTraction, VelocityTrace, VolumeVelocitySource,
};
use super::super::frequency::coupled::poles::PreparedCoupledResonance;
pub use super::super::frequency::coupled::poles::{
    MechanicalPoleResponseOptions, MechanicalResonanceResidue, MechanicalWallResidue,
};
use fs_math::c64::C64;

/// A source projection can fail independently of the accepted physical pole.
#[derive(Debug)]
pub enum MechanicalResonanceResponseError {
    Source(HeterogeneousError),
    Pole(MechanicalPoleError),
}
impl core::fmt::Display for MechanicalResonanceResponseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Source(error) => write!(f, "acoustic resonance source refused: {error}"),
            Self::Pole(error) => write!(f, "{error}"),
        }
    }
}
impl std::error::Error for MechanicalResonanceResponseError {}

/// One reusable bordered factor of the complete cavity/wall/memory equations.
/// Its lifetime binds every source projection and traction to the same FE owner.
pub struct PreparedMechanicalResonance<'p, 'm, 'f, 'c> {
    owner: &'p PreparedMechanicalHarmonic<'m, 'f, 'c>,
    response: PreparedCoupledResonance<'p>,
}

impl<'m, 'f, 'c> PreparedMechanicalHarmonic<'m, 'f, 'c> {
    /// Prepare true source participation for an already resolved simple pole.
    ///
    /// Every original equation and the pole correction are checked again.
    /// Repeated/defective roots can refuse; no inverse wall impedance, modal
    /// damping ratio or symmetry approximation is introduced. The factor is
    /// retained for subsequent coherent source queries at quadratic cost.
    pub fn prepare_resonance_response(
        &self,
        pole: &MechanicalResonancePole,
        options: MechanicalPoleResponseOptions,
        limits: MechanicalPoleLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedMechanicalResonance<'_, 'm, 'f, 'c>, MechanicalPoleError> {
        let response = self.pencil.prepare_resonance_response(
            pole, options, limits, cancelled,
        )?;
        Ok(PreparedMechanicalResonance { owner: self, response })
    }
}

impl<'m, 'f, 'c> PreparedMechanicalResonance<'_, 'm, 'f, 'c> {
    /// Original physical mode, independently checked during preparation.
    /// Its common amplitude is arbitrary; source residues below are absolute.
    #[must_use]
    pub fn mode(&self) -> &MechanicalResonancePole {
        self.response.mode()
    }

    /// Residue for coherent inward motion, boundary/volume injection and wall
    /// traction, all evaluated AT the complex Laplace pole. For a prescribed
    /// frequency-dependent source, the caller supplies its analytic value there;
    /// a real-frequency-only material/source law is not silently continued.
    ///
    /// Native projection, units, original source binding and pressure-release
    /// refusal are identical to the complete harmonic field. Fluid loading is
    /// s*b and mechanical loading is W*traction. Signed phase interference and
    /// exact dark modes are retained; source powers are never added here.
    ///
    /// For exp(+i*omega*t), the angular-frequency residue is -i times this
    /// Laplace residue. It describes one pole's local contribution; obtaining
    /// an actual driven peak, linewidth or power requires the full response.
    #[allow(clippy::too_many_arguments)]
    pub fn forcing_residue(
        &self,
        inward_velocity: &[VelocityTrace],
        boundary_sources: &[(&BoundaryVelocitySource<'_, '_, '_>, C64)],
        volume_sources: &[(&VolumeVelocitySource<'_, '_, '_>, C64)],
        tractions: &[(&PreparedMechanicalHarmonicTraction<'_, 'm, 'f, 'c>, C64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MechanicalResonanceResidue, MechanicalResonanceResponseError> {
        let columns = self.owner.traction_columns(tractions, cancelled)
            .map_err(MechanicalResonanceResponseError::Source)?;
        let load = PreparedHeterogeneous::project_sources_for(
            self.owner.model,
            &self.owner.pressure_release,
            inward_velocity,
            boundary_sources,
            volume_sources,
            cancelled,
        ).map_err(MechanicalResonanceResponseError::Source)?;
        self.response.forcing_residue(&load, &columns, cancelled)
            .map_err(MechanicalResonanceResponseError::Pole)
    }
}

#[cfg(test)]
mod tests;
