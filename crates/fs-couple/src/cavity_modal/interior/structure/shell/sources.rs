//! Distributed original pressure sources driving the original material shell.

use super::{PreparedTensorShellInterior, TensorShellInteriorError, TensorShellInteriorResponse};
use fs_exec::Cx;
use fs_feec::acoustics::helmholtz::heterogeneous::{
    BoundaryVelocitySource, VelocityTrace, VolumeVelocitySource,
};
use fs_math::c64::C64;

impl PreparedTensorShellInterior<'_, '_, '_, '_, '_> {
    /// Drive all original material-shell and pressure coordinates with physical
    /// source columns belonging to the retained heterogeneous interior.
    /// Amplitudes use exp(-i omega t), including volume flow [m³/s] and P2 inward
    /// velocity [m/s]. Linear traces are inward velocity [m/s] as well.
    /// Original full-cycle material strain and reciprocal-work admission runs
    /// unchanged. Total signed source work is already part of the returned
    /// `field.average_acoustic_input_power_w`; it is not extra material heat.
    pub fn solve_with_sources(
        &self,
        structural_force: &[C64],
        inward_velocity: &[VelocityTrace],
        boundary_sources: &[(&BoundaryVelocitySource<'_, '_, '_>, C64)],
        volume_sources: &[(&VolumeVelocitySource<'_, '_, '_>, C64)],
        cx: &Cx<'_>,
    ) -> Result<TensorShellInteriorResponse, TensorShellInteriorError> {
        let load = self.owner.model.project_sources(
            inward_velocity,
            boundary_sources,
            volume_sources,
            cx,
        )?;
        self.solve(structural_force, &load, cx)
    }
}
