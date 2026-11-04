//! Instantaneous momentum transfer through the original gas-piston wall.
//!
//! Only wall velocity changes. Volume/position, gas mass, composition, calorics,
//! source state and the preceding finite-time port account stay unchanged. No
//! gas pressure work, heat or constitutive proposal is assigned to zero elapsed
//! time. Subsequent motion uses the existing source-backed gas/piston step.

use super::{GasPiston, PistonGas, PistonSupport};

mod kernel;
pub use kernel::{PistonImpulseControls, PistonImpulseError, PistonImpulseReceipt};

/// Reciprocal frictionless impact with a native rigid-body material point.
pub mod rigid;
pub use rigid::{
    RigidPistonImpactContact, RigidPistonImpactControls, RigidPistonImpactError,
    RigidPistonImpactReceipt,
};

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), PistonImpulseError> {
    if cancelled() {
        Err(PistonImpulseError::Cancelled)
    } else {
        Ok(())
    }
}

impl<G: PistonGas> GasPiston<G> {
    /// Apply an already resolved external impulse [N s] to this original wall.
    /// Positive impulse points toward primary-gas expansion. A free wall uses
    /// its actual conserved mass; a fixed wall retains zero velocity and returns
    /// the opposite support reaction with zero work. No contact law is inferred.
    ///
    /// The event allocates nothing and makes no gas/material query. All gas and
    /// prior port state remains intact, including on cancellation or refusal.
    pub fn apply_impulse(
        &mut self,
        impulse_n_s: f64,
        controls: PistonImpulseControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PistonImpulseReceipt, PistonImpulseError> {
        poll(&mut cancelled)?;
        let receipt = self.stage_impulse(impulse_n_s, controls)?;
        poll(&mut cancelled)?;
        self.publish_impulse(&receipt);
        Ok(receipt)
    }

    pub(super) fn stage_impulse(
        &self,
        impulse_n_s: f64,
        controls: PistonImpulseControls,
    ) -> Result<PistonImpulseReceipt, PistonImpulseError> {
        kernel::candidate(
            self.mass_kg,
            self.velocity_m_s,
            self.support == PistonSupport::Fixed,
            impulse_n_s,
            controls,
        )
    }

    pub(super) fn publish_impulse(&mut self, receipt: &PistonImpulseReceipt) {
        self.velocity_m_s = receipt.velocity_m_s[1];
    }
}

#[cfg(test)]
mod tests;
