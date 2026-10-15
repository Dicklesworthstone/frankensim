//! Material-bound piston inertia from the existing rigid-body compiler.
//!
//! The body owns regional geometry, density queries and spatial inertia. The
//! piston owns the translational constraint and thermal/momentum time step.
//! Binding them does not create a second mass calculator or rigid-body solver.

use super::{ClosedGasParcel, GasPiston, PistonError, PistonInput, PistonSolve, PistonStep};
use fs_mbd_material::ResolvedBody;

/// A rigid material-bound moving wall and its coupled closed gas.
///
/// All regions of `body` translate together at the piston velocity. Rotations
/// are constrained; centre-of-mass offsets and rotational inertia therefore do
/// not replace its total translational mass with an ad hoc effective mass.
/// Rigidity, geometry coverage and guide reactions remain model assumptions.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialGasPiston {
    body: ResolvedBody,
    piston: GasPiston,
}
impl MaterialGasPiston {
    /// Bind a source-resolved rigid specimen to the existing coupled dynamics.
    ///
    /// `body` is produced by `fs_mbd_material::resolve_body`, retaining density
    /// claims, region geometry and all original material state conditions. The
    /// total mass comes ONLY from its canonical spatial inertia. There is no
    /// competing numeric mass parameter or material-name selection here.
    ///
    /// `area_m2` is the pressure-exposed projected area of the translating wall,
    /// not its volume or the gas volume. Geometry must supply that independent
    /// boundary quantity. The gas and solid may have different temperatures.
    /// This binds a new physical model; it does not replace a live body's mass,
    /// conserve momentum across a material swap or infer thermal expansion.
    pub fn new(
        gas: ClosedGasParcel, body: ResolvedBody, area_m2: f64, velocity_m_s: f64,
    ) -> Result<Self, PistonError> {
        let piston = GasPiston::new(gas, body.inertia().mass(), area_m2, velocity_m_s)?;
        Ok(Self { body, piston })
    }

    /// Original source-bound specimen, including every density usage receipt.
    /// Immutable access prevents changing mass behind the dynamics owner's back.
    #[must_use]
    pub const fn body(&self) -> &ResolvedBody { &self.body }
    /// Existing accepted dynamics and gas state; no duplicate coordinates.
    #[must_use]
    pub const fn piston(&self) -> &GasPiston { &self.piston }
    /// Accepted common thermal/flow/acoustic gas state and source model.
    #[must_use]
    pub const fn gas(&self) -> &ClosedGasParcel { self.piston.gas() }

    /// Advance the original coupled gas and momentum equations unchanged.
    /// Failed input, source support, numerical balance or cancellation leaves
    /// body, gas, velocity and volume intact. The solid specimen is rigid and
    /// its source state is frozen: gas heating is not silently solid heating.
    pub fn step(
        &mut self, input: PistonInput, controls: PistonSolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PistonStep, PistonError> {
        self.piston.step(input, controls, cancelled)
    }
}

#[cfg(test)]
mod tests;
