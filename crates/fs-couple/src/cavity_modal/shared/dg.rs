//! An additional resolved duct downstream of the existing compact cavity neck.
//!
//! The original neck keeps its bore length, inner end inertia and selected loss;
//! the DG geometry is additional length. There is no abrupt area/gas jump, outer
//! neck correction, or second radiation load at this connection. Outlet/branch
//! radiation belongs to the supplied original `DgPortLoad`.

use super::*;
pub mod harmonic;
pub use harmonic::SharedCavityDgHarmonicResponse;
use crate::dg_port::DgPortLoad;
use crate::render::schedule::force::coupled::dg::PreparedModalDg;
pub use crate::render::schedule::force::coupled::dg::{
    ModalDgError as SharedCavityDgError, ModalDgFrame as SharedCavityDgFrame,
};

/// One clock and one transaction for the retained cavity and distributed duct.
/// Neither accepted owner can be mutated outside the coupled step.
pub struct SharedCavityDg {
    cavity: SharedCavity,
    load: DgPortLoad,
    bridge: PreparedModalDg,
    frame: Option<SharedCavityDgFrame>,
}
impl SharedCavity {
    /// Attach an ADDITIONAL downstream duct to the existing neck's outward flow.
    /// Requires `ExteriorMouth`, an ordinary DG inlet, the same inlet gas/area,
    /// an exact clock match, no accepted steps, and no direct exterior radiation
    /// already attached to this neck. An area difference requires a separately
    /// modeled junction; only floating-point roundoff in area is tolerated here.
    /// Existing cavity/neck thermal and viscous history remains unchanged.
    ///
    /// Each step solves the two original mean-port response laws and checks
    /// equal/opposite work before either owner publishes. Pressure is held on
    /// the cavity; its mean swept flow is held on the DG inlet. This is NOT the
    /// continuous joint exponential: coupling-clock refinement remains necessary
    /// independently of DG spatial order and wall-band accuracy. No formal time
    /// order or experimental validity is claimed by the energy check.
    /// Impulsive excitation can make coarse clock refinement nonmonotone;
    /// resolve the retained duct spectrum as well as the compact cavity modes.
    ///
    /// The original modal energy ceiling applies to combined stored energy,
    /// and its pressure ceiling also bounds mean DG inlet pressure. Rational
    /// terminal exchange is reported separately, without a storage certificate.
    /// Fixed source mobility is prepared once; stepping adds no matrix setup
    /// beyond the existing owners' bounded propagation and candidate allocations.
    pub fn with_downstream_duct(
        self,
        load: DgPortLoad,
        gate: &CancelGate,
    ) -> Result<SharedCavityDg, SharedCavityDgError> {
        if gate.is_requested() {
            return Err(ModalCouplingError::Cancelled.into());
        }
        let area = self.source.neck_area_m2();
        if self.source.loss() != CavityNeckLoss::ExteriorMouth
            || self.radiation_neck_port.is_some()
            || self.source.gas() != load.duct().inlet_gas()
            || (area - load.duct().inlet_area_m2()).abs() > 8.0 * f64::EPSILON * area
        {
            return Err(network_invalid("downstream DG requires an internal-only neck, identical inlet gas/area, and no direct exterior neck load").into());
        }
        let mut column = vec![0.0; self.system.mode_count()];
        column[self.surface_modes] = self.source.neck_volume_shape_m2_per_sqrt_kg();
        let bridge = PreparedModalDg::new(&self.system, &load, column, gate)?;
        Ok(SharedCavityDg {
            cavity: self,
            load,
            bridge,
            frame: None,
        })
    }
}
impl SharedCavityDg {
    /// Original cavity geometry, components, material memory and accepted ports.
    /// Its harmonic response describes this cavity subsystem only; it does not
    /// include the downstream DG duct attached by this owning adapter.
    pub const fn cavity(&self) -> &SharedCavity {
        &self.cavity
    }
    /// Original duct, branches, wall states and accepted terminal exchanges.
    pub const fn load(&self) -> &DgPortLoad {
        &self.load
    }
    /// Last fully accepted joint frame, absent before the first step.
    pub const fn last_frame(&self) -> Option<&SharedCavityDgFrame> {
        self.frame.as_ref()
    }
    /// Joint storage represented by the original owners [J].
    pub fn stored_energy_j(&self) -> Result<f64, SharedCavityDgError> {
        Ok(self.cavity.network().total_energy_j()? + self.load.stored_energy_j())
    }
    /// Advance the original surface forces and an additional prescribed neck
    /// pressure. Positive pressure acts inward and its work remains external;
    /// never feed the returned DG reaction back as an additional source.
    /// All accepted mechanical/acoustic/material histories, clocks and frames
    /// survive a rejected step unchanged, including cancellation before commit.
    pub fn step(
        &mut self,
        surface_force: &[f64],
        incident_neck_pressure_pa: f64,
        gate: &CancelGate,
    ) -> Result<&SharedCavityDgFrame, SharedCavityDgError> {
        let before =
            self.cavity
                .prepare_external_step(surface_force, incident_neck_pressure_pa, gate)?;
        let frame = self.bridge.step(
            &mut self.cavity.system,
            &mut self.load,
            &self.cavity.external,
            gate,
        )?;
        self.cavity
            .external
            .copy_from_slice(self.bridge.applied_forces());
        self.cavity.record_accepted_step(before);
        self.frame = Some(frame);
        Ok(self.frame.as_ref().expect("joint step accepted"))
    }
}

#[cfg(test)]
mod tests;
