//! Frozen homogeneous two-phase compressibility (Wood's relation).
//!
//! The two phases move together and share a pressure perturbation. Their masses
//! are frozen over that perturbation; each phase supplies its own adiabatic
//! modulus. Surface tension, bubble inertia/resonance, slip, phase transfer and
//! thermal exchange between phases are absent. This is a long-wavelength,
//! low-frequency mechanical-equilibrium limit within those explicit assumptions,
//! not the equilibrium boiling sound speed or a broadband cavitation model.
//!
//! Brennen, *Cavitation and Bubble Dynamics*, section 6.2, equation (6.13):
//! <https://media.library.caltech.edu/CaltechBOOK:1995.001/chap6.htm>.

use super::{FluidAcousticState, FluidAcoustics};
use crate::MaterialError;

/// Two explicitly supplied phase responses and their actual vapor/void volume
/// fraction. Source identity, common pressure and phase state remain with the
/// original owners; this numerical view does not infer equilibrium from names.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrozenHomogeneousMixture {
    liquid: FluidAcousticState,
    dispersed: FluidAcousticState,
    void_fraction: f64,
    state: FluidAcousticState,
}

impl FrozenHomogeneousMixture {
    /// Use volume fraction alpha, NOT vapor mass quality. The dispersed phase
    /// can be vapor or a separately modeled noncondensable gas. Its own response
    /// is never replaced by an assumed gamma*p law.
    pub fn try_new(
        liquid: FluidAcousticState,
        dispersed: FluidAcousticState,
        void_fraction: f64,
    ) -> Result<Self, MaterialError> {
        if !void_fraction.is_finite() || !(0.0..=1.0).contains(&void_fraction) {
            return Err(MaterialError::Parameters {
                what: "homogeneous mixture requires a finite void volume fraction in [0,1]".into(),
            });
        }
        let state = if void_fraction == 0.0 {
            liquid
        } else if void_fraction == 1.0 {
            dispersed
        } else {
            let alpha = void_fraction;
            let density =
                (1.0 - alpha) * liquid.density_kg_m3() + alpha * dispersed.density_kg_m3();
            let compliance = (1.0 - alpha) * liquid.adiabatic_compressibility_per_pa()
                + alpha * dispersed.adiabatic_compressibility_per_pa();
            FluidAcousticState::from_adiabatic_bulk_modulus(density, 1.0 / compliance)?
        };
        Ok(Self {
            liquid,
            dispersed,
            void_fraction,
            state,
        })
    }

    /// Read the existing material/EOS acoustic trait once for each phase.
    pub fn from_phases(
        liquid: &impl FluidAcoustics,
        dispersed: &impl FluidAcoustics,
        void_fraction: f64,
    ) -> Result<Self, MaterialError> {
        Self::try_new(
            liquid.acoustic_state()?,
            dispersed.acoustic_state()?,
            void_fraction,
        )
    }

    /// Convert the phase-transfer owner's vapor MASS fraction using additive
    /// specific volumes. This is not alpha=x or a linear density interpolation.
    pub fn from_mass_fraction(
        liquid: FluidAcousticState,
        dispersed: FluidAcousticState,
        dispersed_mass_fraction: f64,
    ) -> Result<Self, MaterialError> {
        let x = dispersed_mass_fraction;
        if !x.is_finite() || !(0.0..=1.0).contains(&x) {
            return Err(MaterialError::Parameters {
                what: "homogeneous mixture requires a finite dispersed mass fraction in [0,1]"
                    .into(),
            });
        }
        if x == 0.0 || x == 1.0 {
            return Self::try_new(liquid, dispersed, x);
        }
        let liquid_volume = (1.0 - x) / liquid.density_kg_m3();
        let dispersed_volume = x / dispersed.density_kg_m3();
        let volume = liquid_volume + dispersed_volume;
        let alpha = dispersed_volume / volume;
        if !volume.is_finite() || volume <= 0.0 || alpha <= 0.0 || alpha >= 1.0 {
            return Err(MaterialError::Parameters {
                what: "two-phase mass-to-volume fraction conversion is not representable".into(),
            });
        }
        Self::try_new(liquid, dispersed, alpha)
    }

    /// Original liquid acoustic response, without mixture changes.
    pub const fn liquid(self) -> FluidAcousticState {
        self.liquid
    }
    /// Original dispersed-phase acoustic response.
    pub const fn dispersed(self) -> FluidAcousticState {
        self.dispersed
    }
    /// Actual dispersed-phase volume fraction alpha.
    pub const fn void_fraction(self) -> f64 {
        self.void_fraction
    }
    /// Homogeneous mixture response from phase volume and compliance sums.
    pub const fn state(self) -> FluidAcousticState {
        self.state
    }
    /// Actual dispersed mass fraction corresponding to the retained alpha.
    pub fn dispersed_mass_fraction(self) -> f64 {
        self.void_fraction * self.dispersed.density_kg_m3() / self.state.density_kg_m3()
    }
}

impl FluidAcoustics for FrozenHomogeneousMixture {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        Ok(self.state)
    }
}
