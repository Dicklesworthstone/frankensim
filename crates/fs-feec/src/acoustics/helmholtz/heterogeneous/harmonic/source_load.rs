//! Original source loads for complete pressure-potential multiphysics pencils.

use super::{
    BoundaryVelocitySource, HeterogeneousError, HeterogeneousHelmholtzFe, PreparedHeterogeneous,
    VelocityTrace, VolumeVelocitySource, finite, invalid, poll,
};

use crate::acoustics::helmholtz::pressure_release::EmptyFreePressure;
use fs_math::c64::C64;

impl HeterogeneousHelmholtzFe<'_, '_> {
    /// Combine physical injection in the complete original pressure basis.
    ///
    /// For `beta * p_dot + div(v) = s`, the pressure potential `p = psi_dot`
    /// receives `b_i = integral(phi_i * s dV) + integral(phi_i * v_in dA)`.
    /// Each returned entry is a generalized volume flow [m³/s]. There is NO
    /// density or frequency factor: the pressure-potential equation receives
    /// `b`, whereas a pressure equation receives its time derivative.
    /// The conjugate peak-phasor work is `Re(sum(conj(p_i) * b_i))/2`.
    ///
    /// Prepared shapes are real spatial columns with no phasor convention.
    /// This projection preserves the supplied complex amplitudes unchanged;
    /// a receiving exp(-i omega t) solver requires exp(-i omega t) amplitudes.
    /// Conjugate an amplitude when converting an existing exp(+i omega t)
    /// drive, rather than adding a sign, omega, density or area here.
    ///
    /// Source count, original material/geometry identity, finite amplitudes,
    /// complete boundary support and pressure-release exclusions use the same
    /// admission as the native pressure solve, even for zero drive amplitudes.
    /// The receiving pencil must retain the stated pressure-release condition;
    /// output remains in the complete original coordinate order. Empty faces
    /// admit a full pressure-potential pencil with no eliminated coordinates.
    /// The native signed P2 and volume quadrature is reused without resampling,
    /// normalization or inversion of a bare acoustic operator.
    pub fn project_source_load(
        &self,
        pressure_release_faces: &[[u32; 3]],
        inward_velocity: &[VelocityTrace],
        boundary_sources: &[(&BoundaryVelocitySource<'_, '_, '_>, C64)],
        volume_sources: &[(&VolumeVelocitySource<'_, '_, '_>, C64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<C64>, HeterogeneousError> {
        let constraints = self.geometry.pressure_release_constraints(
            pressure_release_faces,
            EmptyFreePressure::Refuse,
            cancelled,
        )?;
        let load = PreparedHeterogeneous::project_sources_for(
            self,
            &constraints,
            inward_velocity,
            boundary_sources,
            volume_sources,
            cancelled,
        )?;
        for &value in &load {
            poll(cancelled)?;
            if !finite(value) {
                return Err(invalid("combined physical source load overflow").into());
            }
        }
        Ok(load)
    }
}
