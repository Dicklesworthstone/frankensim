//! Original body impulse/work inside the moving reacting transaction.

use super::*;

impl ReactingEuler3d {
    /// Advance C/2, K/2, ALE, K/2, C/2 with the original held body source.
    ///
    /// Acceleration uses actual stage mass; force density uses actual old/new
    /// physical volume. The original source map supplies represented impulse
    /// and kinetic work without changing species or inventing thermal energy.
    /// The same entry consumes the original midpoint or exact Coriolis source.
    /// Coefficients are held in the original cell order throughout one trial;
    /// they are not resampled as an Eulerian field on the moving grid.
    ///
    /// Physical Mach and mesh-relative acoustic admission remain separate.
    /// Chemistry, both kicks and geometry restart together on shortening, and
    /// accepted cells, clock and source/boundary reports publish atomically.
    /// A zero source follows the existing moving composition after admission.
    #[allow(clippy::too_many_arguments)]
    pub fn step_moving_with_forcing(
        &mut self,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        motion: MeshMotion,
        force: BodyForce3d<'_>,
        force_controls: ForceControl3d,
        workspace: &mut ForcingWorkspace3d,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MovingReactiveStepReport3d, ReactiveError> {
        self.step_moving_impl::<ConstantTransport>(
            maximum_duration,
            controls,
            None,
            DissipationControl::default(),
            super::super::super::MAX_MIXTURE_WORKSPACE_BYTES,
            motion,
            Some((force, force_controls, workspace)),
            cancelled,
        )
    }

    /// Add original molecular half intervals around the complete forced
    /// reacting ALE split. Molecular stress work, wall heat, body work and the
    /// chemical reference-energy source retain separate physical accounts.
    ///
    /// The explicit budget covers both reacting fields, the prepared reactor,
    /// the molecular operator and forcing scratch together. Refusal preserves
    /// both accepted owners; numerical scratch may be reused. Closed/periodic
    /// moving-domain admission is unchanged. This source splitting is not a
    /// hydrostatic well-balanced method or a moving reference-frame transform.
    #[allow(clippy::too_many_arguments)]
    pub fn step_moving_with_dissipation_and_forcing<L: TransportLaw>(
        &mut self,
        operator: &mut MolecularTransport3d<L>,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        dissipation_controls: DissipationControl,
        maximum_workspace_bytes: usize,
        motion: MeshMotion,
        force: BodyForce3d<'_>,
        force_controls: ForceControl3d,
        workspace: &mut ForcingWorkspace3d,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MovingReactiveStepReport3d, ReactiveError> {
        self.step_moving_impl(
            maximum_duration,
            controls,
            Some(operator),
            dissipation_controls,
            maximum_workspace_bytes,
            motion,
            Some((force, force_controls, workspace)),
            cancelled,
        )
    }

    /// Admit and normalize the held source before any physical proposal.
    /// Retain its allocated bytes in the budget even for a zero source.
    pub(in crate::compressible::spatial::mixture::reactive) fn prepare_moving_forcing(
        &self,
        forcing: &mut Option<forcing::Context<'_, '_>>,
        maximum_workspace_bytes: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<usize, ReactiveError> {
        let bytes = forcing
            .as_ref()
            .map_or(0, |(_, _, work)| work.workspace_bytes());
        let remaining = maximum_workspace_bytes
            .min(super::super::super::MAX_MIXTURE_WORKSPACE_BYTES)
            .checked_sub(bytes)
            .ok_or(ReactiveError::Invalid(
                "moving reactive forcing workspace budget",
            ))?;
        let active = if let Some((force, control, work)) = forcing.as_mut() {
            control.admit::<super::super::super::Error>()?;
            work.admit(&self.accepted)?;
            force.admit::<super::super::super::Error>(self.accepted.cell_count(), cancelled)?
        } else {
            false
        };
        if !active {
            *forcing = None;
        }
        Ok(remaining)
    }
}
