//! Explicit full-gradient SGS in the original moving-wall collision.

use super::{
    CancelGate, FreeSurface3, FreeSurfaceSgsError3, FreeSurfaceSgsReport3,
    FreeSurfaceSgsWorkspace3, SgsConfig, check_cancel,
};
use crate::d3q19::freesurface3::{
    FreeSurfaceKinematicsWorkspace3, PlanarWallFieldControls3, PlanarWallFieldReport3,
    PlanarWallVelocityField3,
};

impl FreeSurface3 {
    /// Propose one wall-driven VOF tick with explicit WALE/Vreman viscosity.
    ///
    /// The accepted interior liquid velocity supplies the full gradient, held
    /// for this tick; no wall-law gradient or subgrid wall stress is inferred.
    /// Empty interface carriers keep molecular relaxation. Occupied cells retain
    /// the gradient closure's second-order support and viscosity-ratio gates.
    ///
    /// The SAME collision supplies the populations reflected by the original
    /// face-field wall operator. Its force, torque, couple and actuator work
    /// therefore belong to this modeled-viscosity trajectory. The returned SGS
    /// report describes the input state, not thermal energy deposited by the wall.
    /// Native wall geometry, speed, face-mass and whole-liquid checks still apply.
    ///
    /// The wall/domain budgets and gradient workspace budget precede the single
    /// native candidate clone. Zero coefficient takes the original None-tau wall
    /// path, preserving native populations and complete wall receipts exactly.
    pub fn planar_wall_velocity_field_gradient_sgs_candidate(
        &self,
        motion: PlanarWallVelocityField3<'_>,
        controls: PlanarWallFieldControls3,
        gate: &CancelGate,
        config: SgsConfig,
        kinematics: &mut FreeSurfaceKinematicsWorkspace3,
        workspace: &mut FreeSurfaceSgsWorkspace3,
    ) -> Result<
        (Self, PlanarWallFieldReport3, FreeSurfaceSgsReport3),
        FreeSurfaceSgsError3,
    > {
        check_cancel(gate)?;
        controls.admit_shape(
            motion.wall,
            self.domain_cells(),
            motion.tangential_velocity.len(),
        )?;
        let sgs = self.prepare_gradient_sgs(gate, config, kinematics, workspace)?;
        let tau = (config.model.coefficient() != 0.0).then_some(&workspace.tau);
        let (candidate, wall) =
            self.planar_wall_field_candidate_with_relaxation(motion, controls, tau, gate)?;
        check_cancel(gate)?;
        Ok((candidate, wall, sgs))
    }

    /// Publish the combined native wall/SGS candidate only after all checks.
    /// Returns actual wall exchange followed by pre-step liquid SGS diagnostics.
    pub fn step_with_planar_wall_velocity_field_and_gradient_sgs(
        &mut self,
        motion: PlanarWallVelocityField3<'_>,
        controls: PlanarWallFieldControls3,
        gate: &CancelGate,
        config: SgsConfig,
        kinematics: &mut FreeSurfaceKinematicsWorkspace3,
        workspace: &mut FreeSurfaceSgsWorkspace3,
    ) -> Result<(PlanarWallFieldReport3, FreeSurfaceSgsReport3), FreeSurfaceSgsError3> {
        let (candidate, wall, sgs) = self.planar_wall_velocity_field_gradient_sgs_candidate(
            motion,
            controls,
            gate,
            config,
            kinematics,
            workspace,
        )?;
        check_cancel(gate)?;
        *self = candidate;
        Ok((wall, sgs))
    }
}
