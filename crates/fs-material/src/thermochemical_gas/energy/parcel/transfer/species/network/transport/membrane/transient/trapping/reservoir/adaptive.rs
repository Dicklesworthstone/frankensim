//! Every temporal probe includes original gas sources and EOS reconstruction.
use super::super::super::slab::{
    adaptive::{MembraneDiffusionAdvanceError, MembraneDiffusionAdvanceWorkspace, Reservoirs},
    trapping::adaptive::{
        self as kinetics, KineticNative, TrappingMembraneAdvanceControls,
        TrappingMembraneAdvanceReport,
    },
};
use super::*;

#[derive(Clone, Debug, PartialEq)]
pub struct IsothermalTrappingGasAdvanceReport {
    pub diffusion: TrappingMembraneAdvanceReport,
    /// Net accepted fixed-temperature streams, including reversals. The
    /// individual original gas caloric balances were checked on every probe.
    /// This account is not a summed or whole-film energy residual certificate.
    pub accounting: IsothermalTrappingGasAccount,
}
impl<'a> KineticNative for IsothermalTrappingGasMembrane<'a> {
    type Error = IsothermalTrappingGasError;
    fn kinetics(&self) -> &TrappingMembraneDiffusionState {
        self.membrane()
    }
    fn reservoirs(&self) -> Option<Reservoirs> {
        Some(Reservoirs {
            capacities: self.finite.capacities().map(|c| c.value_si()),
            pressures: self.finite.reservoir_pressures_pa(),
        })
    }
    fn is_cancelled(error: &Self::Error) -> bool {
        *error == IsothermalTrappingGasError::Cancelled
    }
    fn kinetic_step(
        &self,
        _: [f64; 2],
        duration: f64,
        controls: TrappingMembraneSolve,
        work: &mut TrappingMembraneWorkspace,
        cancelled: &mut impl FnMut() -> bool,
    ) -> core::result::Result<(Self, TrappingMembraneStep), Self::Error> {
        let (state, report) = self
            .step_with_controls(duration, controls, work, cancelled)?
            .into_parts();
        Ok((state, report.diffusion))
    }
}
impl<'a> IsothermalTrappingGasMembrane<'a> {
    pub fn adaptive_workspace(
        &self,
    ) -> Result<MembraneDiffusionAdvanceWorkspace<TrappingMembraneWorkspace>> {
        Ok(self.finite.adaptive_workspace()?)
    }
    /// Advance actual gas parcels and both membrane stores atomically. Every
    /// coarse/fine probe calls the SAME source-relative gas reconstruction;
    /// gas source failures cannot be concealed by time-step adaptation. The
    /// shared algebraic audit includes actual EOS-reseeded reservoir pressures.
    pub fn advance_candidate(
        &self,
        duration_s: f64,
        controls: TrappingMembraneAdvanceControls,
        work: &mut MembraneDiffusionAdvanceWorkspace<TrappingMembraneWorkspace>,
        mut cancelled: impl FnMut() -> bool,
    ) -> core::result::Result<
        IsothermalTrappingGasCandidate<'a, IsothermalTrappingGasAdvanceReport>,
        MembraneDiffusionAdvanceError<IsothermalTrappingGasError>,
    > {
        let (state, diffusion) =
            kinetics::advance(self, [0.0; 2], duration_s, controls, work, &mut cancelled)?
                .into_parts();
        let amounts = [
            -diffusion.transport.face_amounts_mol[0],
            *diffusion
                .transport
                .face_amounts_mol
                .last()
                .expect("nonempty membrane faces"),
        ];
        let accounting = self
            .account(
                &state,
                amounts,
                diffusion.mobile_change_mol,
                diffusion.trapped_change_mol,
            )
            .map_err(MembraneDiffusionAdvanceError::Native)?;
        if cancelled() {
            return Err(MembraneDiffusionAdvanceError::Cancelled);
        }
        Ok(IsothermalTrappingGasCandidate {
            state,
            report: IsothermalTrappingGasAdvanceReport {
                diffusion,
                accounting,
            },
        })
    }
}
