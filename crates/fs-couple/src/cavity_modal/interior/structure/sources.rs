//! Original distributed injection driving complete structural/acoustic fields.

use super::{
    PreparedStructuralInterior, StructuralInterior, StructuralInteriorError,
    StructuralInteriorResponse, poll,
};
use fs_exec::Cx;
use fs_feec::acoustics::helmholtz::heterogeneous::{
    BoundaryVelocitySource, VelocityTrace, VolumeVelocitySource,
};
use fs_math::c64::C64;

impl StructuralInterior<'_, '_, '_> {
    /// Project linear boundary motion, original P2 boundary sources and volume
    /// injection in the complete pressure basis. Sources must belong to this
    /// exact original heterogeneous material/geometry owner, including sources
    /// with zero amplitude. Every physical shape and signed coefficient survives.
    ///
    /// All drives use peak exp(-i omega t) amplitudes. The real spatial source
    /// columns have no phasor convention; conjugate amplitudes imported from an
    /// exp(+i omega t) drive. The pressure-potential RHS receives this [m³/s]
    /// load directly, with no additional time derivative or density scaling.
    /// The original pencil retains every pressure coordinate, so no
    /// pressure-release face is declared. Prescribed inward motion is additive
    /// relative to the coupled structural motion and any passive wall flow.
    pub fn project_sources(
        &self,
        inward_velocity: &[VelocityTrace],
        boundary_sources: &[(&BoundaryVelocitySource<'_, '_, '_>, C64)],
        volume_sources: &[(&VolumeVelocitySource<'_, '_, '_>, C64)],
        cx: &Cx<'_>,
    ) -> Result<Vec<C64>, StructuralInteriorError> {
        poll(cx)?;
        Ok(self.interior.project_source_load(
            &[],
            inward_velocity,
            boundary_sources,
            volume_sources,
            &mut || cx.checkpoint().is_err(),
        )?)
    }
}

impl PreparedStructuralInterior<'_, '_, '_, '_> {
    /// Drive the existing complete factor with original physical injection and
    /// structural forces together. Source work is already included ONCE in
    /// `average_acoustic_input_power_w`, paired with the solved pressure rather
    /// than the pressure potential. It is separate from material/wall heat.
    pub fn solve_with_sources(
        &self,
        structural_force: &[C64],
        inward_velocity: &[VelocityTrace],
        boundary_sources: &[(&BoundaryVelocitySource<'_, '_, '_>, C64)],
        volume_sources: &[(&VolumeVelocitySource<'_, '_, '_>, C64)],
        cx: &Cx<'_>,
    ) -> Result<StructuralInteriorResponse, StructuralInteriorError> {
        let load =
            self.model
                .project_sources(inward_velocity, boundary_sources, volume_sources, cx)?;
        self.solve(structural_force, &load, cx)
    }
}
