//! Fixed support for the original gas, material and mechanical state owner.
//!
//! This is an ideal stationary constraint. It supplies a reaction impulse,
//! never removes kinetic energy, invents an inertial mass or credits heat.

use super::{GasPiston, PistonError, PistonGas, PistonGasState, PistonInput, PistonSolve, Trial, invalid, poll};

/// Mechanical boundary data retained through gas reaction, flow and heat steps.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PistonSupport {
    /// The original finite-mass translational momentum equation.
    Free,
    /// An ideal support holds the accepted wall position exactly.
    Fixed,
}

impl<G: PistonGas> GasPiston<G> {
    /// Attach an ideal fixed support to a wall already at rest.
    ///
    /// The physical wall mass, exposed area, source gas and volume are retained.
    /// Moving walls refuse: stopping one requires an independently accounted
    /// collision or braking impulse and must not silently discard kinetic energy.
    /// This constraint is retained by every existing reaction/flow/heat candidate.
    pub fn with_fixed_support(mut self) -> Result<Self, PistonError> {
        if self.velocity_m_s != 0.0 {
            return Err(invalid("fixed support requires a wall at rest; account for its stopping impulse and energy first"));
        }
        self.support = PistonSupport::Fixed;
        Ok(self)
    }

    /// Accepted mechanical constraint; no separate gas or geometry owner.
    #[must_use]
    pub const fn support(&self) -> PistonSupport { self.support }

    pub(super) fn fixed_trial(
        &self, input: PistonInput, controls: PistonSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Trial<G>, PistonError> {
        poll(cancelled)?;
        if self.velocity_m_s != 0.0 {
            return Err(invalid("fixed support cannot carry nonzero wall velocity"));
        }
        let gas = self.gas.midpoint_heat_volume(
            self.gas.volume_m3(), input.heat_j, controls.gas, &mut *cancelled,
        )?;
        let force = self.area_m2*(gas.midpoint().pressure_pa()-input.exterior_pressure_pa)
            + input.external_force_n;
        let support_impulse = -input.duration_s*force;
        if !force.is_finite() || !support_impulse.is_finite() {
            return Err(invalid("fixed-support reaction impulse is unrepresentable"));
        }
        // No root or volume perturbation exists for this constraint. The unused
        // zero tangent never enters the freely moving wall's impulse solver.
        Ok(Trial { gas, displacement: 0.0, velocity: 0.0, residual: 0.0,
            tangent: 0.0, support_impulse })
    }
}

#[cfg(test)]
mod tests;
