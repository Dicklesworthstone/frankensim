//! Prepared molecular rheology in the original moving-wall VOF transaction.

use super::*;
use crate::d3q19::freesurface3::{
    PlanarWallFieldControls3, PlanarWallFieldReport3, PlanarWallVelocityField3,
};

impl FreeSurfaceRheology3 {
    /// Advance the actual non-Newtonian liquid under a tangential wall field.
    ///
    /// The accepted local constitutive solution supplies the even relaxation
    /// in the same native collision that reflects these moving-wall links.
    /// BGK remains BGK; TRT retains its authored even/odd product at every wet
    /// cell. Original forcing, gas reconstruction, capillarity, mass exchange
    /// and conversion advance once on one private native candidate.
    ///
    /// The actual endpoint populations, including wall-generated stress and
    /// the Guo half-force correction, then resolve the same material law before
    /// publication. No stationary-wall or gas velocity is substituted into a
    /// gradient. Newly wet cells use their original global field coordinates.
    /// Failed wall, native, constitutive or cancellation checks preserve the
    /// accepted flow, relaxation map and constitutive statistics together.
    ///
    /// Returns the constitutive statistics APPLIED in this collision and the
    /// original face/total impulse, angular impulse, couple and actuator work.
    /// `stats()` describes the accepted endpoint, ready for the next call.
    /// The two relaxation maps reuse existing entries; native wall storage and
    /// its explicit domain/face budgets retain their original behavior.
    ///
    /// Velocities, wall geometry and exchange retain native lattice units.
    /// The plane must remain fixed and completely flooded, with tangent-only
    /// material motion. This is an explicit isothermal constitutive step; it
    /// adds no caloric heating, implicit wall/material iteration, moving contact
    /// line, normal displacement, or independently actuated gas pressure.
    pub fn step_with_planar_wall_velocity_field(
        &mut self,
        motion: PlanarWallVelocityField3<'_>,
        controls: PlanarWallFieldControls3,
        gate: &CancelGate,
    ) -> Result<(RheologyStats, PlanarWallFieldReport3), FreeSurfaceRheologyError3> {
        check_cancel(gate)?;
        let (candidate, exchange) = self.lattice.planar_wall_field_candidate_with_relaxation(
            motion,
            controls,
            Some(&self.tau),
            gate,
        )?;
        let applied = self.accept_candidate(candidate, gate)?;
        Ok((applied, exchange))
    }
}
