//! Free decay of the same complete pressure/wall realization as forced response.

use super::PreparedMechanicalHarmonic;
use fs_math::c64::C64;

mod response;
pub use response::{
    MechanicalPoleResponseOptions, MechanicalResonanceResidue,
    MechanicalResonanceResponseError, MechanicalWallResidue, PreparedMechanicalResonance,
};
pub use super::frequency::coupled::poles::{
    MechanicalPoleError, MechanicalPoleLimits, MechanicalPoleOptions, MechanicalPoleWall,
    MechanicalResonancePole, MechanicalResonanceRegion, MechanicalResonanceSpectrum,
    MechanicalSpectrumOptions,
};

impl PreparedMechanicalHarmonic<'_, '_, '_> {
    /// Full finite-pencil coordinates: free pressure, then each wall's original
    /// displacement, velocity and Maxwell-extension blocks in wall-DOF order.
    /// The Maxwell extension is `sqrt(k_j)*w-z_j`, without eliminating any arm.
    #[must_use]
    pub fn resonance_coordinate_count(&self) -> usize {
        self.pencil.resonance_coordinate_count()
    }

    /// Original pressure coordinates occupying the beginning of a pole seed.
    #[must_use]
    pub fn resonance_free_pressure_dofs(&self) -> &[usize] {
        &self.pressure_release.free_dofs
    }

    /// Refine one nonzero simple finite mode of the original complete cavity
    /// and causal wall equations under `exp(s*t)`. Pure real relaxation and
    /// unobservable memory modes remain available; their rates are not forced
    /// peak frequencies or half-power Q values.
    ///
    /// The source was prepared at a nonsingular real reference frequency, but
    /// the query evaluates its ORIGINAL polynomial at complex `s`, with every
    /// Maxwell arm, consistent mass and reciprocal face coupling intact. Only
    /// the numerical coordinate scaling is reused from that reference. No wall
    /// impedance is inverted and no material loss is frozen at a trial rate.
    ///
    /// `seed` follows `resonance_coordinate_count` in physical coordinates.
    /// Fixed pressure rows lift to exact zeros. Each original free equation is
    /// checked independently after the shared Newton solve. Exact zero-rate
    /// pressure/gauge classification is outside this API; repeated/defective
    /// roots may refuse. All returned modes describe the retained spatial model,
    /// without a continuum or experimental validation claim.
    pub fn resonance_pole(
        &self,
        rate_per_s: C64,
        seed: &[C64],
        options: &MechanicalPoleOptions,
        limits: MechanicalPoleLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MechanicalResonancePole, MechanicalPoleError> {
        self.pencil.resonance_pole(rate_per_s, seed, options, limits, cancelled)
    }

    /// Discover numerical nonzero finite modes inside a declared Laplace-plane
    /// ellipse, using the shared Beyn/Newton owners. The ellipse must exclude
    /// zero. Counts/rank agreement and nested quadrature are numerical checks,
    /// not an interval completeness certificate. No partial spectrum is
    /// returned if any original mode, budget or cancellation check fails.
    pub fn resonance_spectrum(
        &self,
        region: MechanicalResonanceRegion,
        options: &MechanicalSpectrumOptions,
        limits: MechanicalPoleLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MechanicalResonanceSpectrum, MechanicalPoleError> {
        self.pencil.resonance_spectrum(region, options, limits, cancelled)
    }
}

#[cfg(test)]
mod tests;
