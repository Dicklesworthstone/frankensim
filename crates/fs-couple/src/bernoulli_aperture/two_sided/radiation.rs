//! The shared source-bound MIMO exterior on the original moving-valve owner.
//!
//! Each acoustic side can own one ordered radiation group. Its pressure returns
//! through that side's original network to the single, paired or modal valve;
//! no alternate reed engine, extra feedback delay or output coloration is added.
//! Cross-side exterior coupling is not represented by two independent groups.

use super::{TwoSidedAperture, ValveSide, invalid, poll, wave};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::radiation::environment::multiport::{
    MimoEnvironment, MimoRadiationLimits, MimoRadiationModel,
};
use crate::pcm_wav::baffled::openings::duct::DuctRadiationSpec;
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use fs_math::c64::C64;
use fs_vfit::mimo::MimoPassivityReport;
use fs_vfit::waveguide::network::multiport::MimoPortFrame;

impl TwoSidedAperture {
    /// Attach the original source-bound exterior to one unadvanced acoustic side.
    /// `radiation` addresses this side's local graph: every opening's `network`
    /// must be zero, in the source matrix's exact row/column order. A later
    /// `TwoSidedPressure` uses global indices 0 (upstream) and 1 (downstream).
    ///
    /// Original valve materials, wall/supply histories, regional gases, loss
    /// subdivisions and fractional delays remain with their existing owners.
    /// Both sides may carry a group; the groups do not mutually couple across
    /// the valve. Each later receiver must match the retained source gas,
    /// aperture geometry/order and mechanics clock, including Diagnostic mode.
    ///
    /// This consumes an externally produced proper undelayed model; it does not
    /// generate or refit radiation data. The existing continuous and GridOnly
    /// sampled checks do not establish rational storage or complete M2.5-W3c
    /// qualification. Signed exchanged work is separately audited in each frame.
    ///
    /// # Errors
    /// An advanced valve, a second group, mismatched clock/gas/geometry, a
    /// preexisting reactive terminal, failed original model/resource admission
    /// or cancellation publishes no partially attached system.
    pub fn with_mimo_environment(
        mut self,
        side: ValveSide,
        source: MimoRadiationModel,
        radiation: DuctRadiationSpec,
        limits: MimoRadiationLimits,
        gate: &CancelGate,
    ) -> Result<Self, AcousticRealizeError> {
        poll(gate)?;
        let mechanics = radiation.reconstruction.mechanics_rate_hz();
        if self.valve.accepted_steps() != 0 || mechanics == 0
            || self.valve.time_step_s().to_bits() != (1.0 / f64::from(mechanics)).to_bits()
            || self.sides[side.index()].mimo_environment.is_some()
        {
            return Err(invalid("MIMO exterior requires an unadvanced valve, its exact clock and an unbound side"));
        }
        let side = &mut self.sides[side.index()];
        let (environment, terminal) = MimoEnvironment::prepare(
            &side.graph, &side.gases, source, radiation, limits, gate,
        )?;
        side.network.install_mimo_terminals(terminal).map_err(wave)?;
        poll(gate)?;
        side.mimo_environment = Some(environment);
        Ok(self)
    }

    /// Original exterior, aperture order and source identities on this side.
    #[must_use]
    pub fn mimo_environment(&self, side: ValveSide) -> Option<&MimoEnvironment> {
        self.sides[side.index()].mimo_environment.as_ref()
    }

    /// Accepted pressure, volume flow and signed work in the source matrix order.
    #[must_use]
    pub fn mimo_port_frames(&self, side: ValveSide) -> Option<&[MimoPortFrame]> {
        self.sides[side.index()].network.mimo_port_frames()
    }

    /// Continuous report and separate actual sampled-response diagnostic.
    #[must_use]
    pub fn mimo_passivity(&self, side: ValveSide)
        -> Option<(&MimoPassivityReport, &MimoPassivityReport)>
    {
        self.sides[side.index()].network.mimo_passivity()
    }

    /// Accepted original rational pole coordinates; no copied or reset history.
    #[must_use]
    pub fn mimo_pole_states(&self, side: ValveSide) -> Option<&[C64]> {
        self.sides[side.index()].network.mimo_pole_states()
    }

    // Remap only addresses. The common observer's medium, geometry, reconstruction
    // and per-side ordering are then checked by the same owner used by flute.
    pub(super) fn validate_mimo_receiver(
        &self,
        radiation: &DuctRadiationSpec,
        exterior: GasState,
    ) -> Result<(), AcousticRealizeError> {
        for (network, side) in self.sides.iter().enumerate() {
            let Some(environment) = &side.mimo_environment else { continue; };
            let mut local = radiation.clone();
            local.openings.retain(|opening| opening.network == network);
            for opening in &mut local.openings { opening.network = 0; }
            environment.validate_receiver(&side.graph, &local, exterior)?;
        }
        Ok(())
    }
}
