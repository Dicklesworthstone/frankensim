//! Body-source impulse and work in the original reacting piston composition.

use super::*;

impl ReactingEuler3d {
    /// Couple reacting ALE gas, held body forcing and both finite pistons.
    ///
    /// The original source kicks surround the central moving advection.
    /// Gas body impulse and kinetic work remain in the reacting report and
    /// are removed from the closed gas/piston mechanical residual. Boundary
    /// impulse/work still drive the same simultaneous finite-mass solve.
    ///
    /// `force` acts on the gas. Forces on the pistons remain explicit in
    /// `loads`; gravity or rotation of a complete assembly is not inferred.
    /// Coefficients are held in original cell order through each complete
    /// trial. Failed source, wave or mechanical proposals publish nothing.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_inertial_walls_and_forcing(
        &mut self,
        axis: usize,
        walls: &mut [Piston; 2],
        maximum_duration: Time,
        transport_controls: ReactiveStepControl,
        loads: [InertialWallLoad; 2],
        controls: InertialWallsControl,
        force: BodyForce3d<'_>,
        force_controls: ForceControl3d,
        workspace: &mut ForcingWorkspace3d,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveInertialWallsStep, ReactiveInertialWallsError> {
        self.piston_step::<ConstantTransport>(
            None,
            axis,
            walls,
            maximum_duration,
            transport_controls,
            DissipationControl::default(),
            super::super::super::MAX_MIXTURE_WORKSPACE_BYTES,
            Some((force, force_controls, workspace)),
            loads,
            controls,
            cancelled,
        )
    }

    /// Add native molecular stress, diffusion and heat to the forced reacting
    /// piston step. All source halves, geometry and mechanical proposals share
    /// one accepted physical interval and reuse their prepared workspaces.
    ///
    /// Chemical reference changes, gas body work, wall heat, boundary work,
    /// fixed-side support impulse and guide reactions retain their independent
    /// accounts. The explicit byte budget includes forcing scratch as well as
    /// both reacting frames, the reactor and the molecular operator. All gas,
    /// source-binding and piston publication follows the final checkpoint.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_dissipation_forcing_and_inertial_walls<L: TransportLaw>(
        &mut self,
        operator: &mut MolecularTransport3d<L>,
        axis: usize,
        walls: &mut [Piston; 2],
        maximum_duration: Time,
        transport_controls: ReactiveStepControl,
        dissipation_controls: DissipationControl,
        maximum_workspace_bytes: usize,
        loads: [InertialWallLoad; 2],
        controls: InertialWallsControl,
        force: BodyForce3d<'_>,
        force_controls: ForceControl3d,
        workspace: &mut ForcingWorkspace3d,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveInertialWallsStep, ReactiveInertialWallsError> {
        self.piston_step(
            Some(operator),
            axis,
            walls,
            maximum_duration,
            transport_controls,
            dissipation_controls,
            maximum_workspace_bytes,
            Some((force, force_controls, workspace)),
            loads,
            controls,
            cancelled,
        )
    }
}
