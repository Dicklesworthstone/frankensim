//! Finite-rate compression of one retained regional liquid/vapor water state.
//!
//! The original frozen and equilibrium acoustic owners supply K_f and K_eq.
//! A caller-declared positive compression-relaxation time supplies ONE Maxwell
//! arm: K(omega) = K_eq + (K_f-K_eq) s*tau/(1+s*tau), s=-i*omega.
//! The original real density then gives compressibility 1/K, density response
//! d(rho)/dp = rho/K, wavenumber omega*sqrt(rho/K), and impedance sqrt(rho*K).
//! These are one coherent response, not independently interpolated speeds.
//!
//! For compressive strain e and an internal relaxed strain q, the law is
//! p = K_eq*e + (K_f-K_eq)*(e-q), tau*dq/dt = e-q. Stored energy per reference
//! volume is [K_eq*e^2 + (K_f-K_eq)*(e-q)^2]/2; dissipation is
//! (K_f-K_eq)*(e-q)^2/tau >= 0. This is the existing positive compression arm,
//! not a bulk viscosity inferred from a propagation attenuation coefficient.
//!
//! IF97 supplies the endpoints, NOT a phase-transfer rate. The declared tau is
//! the relaxation time at fixed compression. The equivalent compliance time
//! at fixed pressure is tau*K_f/K_eq. This single-time constitutive assumption
//! needs caller-established rate data and a validity band; it does not derive
//! heat/mass-transfer kinetics, phase-specific temperatures or a quality history.
//! Perturbations must be small and remain strictly inside coexistence. Common
//! pressure/velocity and homogenization must be justified, with negligible slip,
//! bubble inertia, interface energy and noncondensable gas. No molecular,
//! viscothermal or wall losses are added. Source fit/root defects remain with
//! the original endpoint owners and are never converted into dissipation.

/// Exact compression histories with fixed-mass density and physical arm heat.
pub mod transient;

use super::super::regional::RegionalWaterState;
use super::thermal::{CoexistenceAcousticClosure, RegionalWaterAcoustics};
use crate::MaterialError;
use crate::relaxation::bulk::BulkRelaxationArm;
use fs_math::c64::C64;

/// One explicitly declared compression timescale on original water endpoints.
///
/// Construction and frequency evaluation allocate no successful-path storage,
/// repeat no EOS/root solve, and do not change the thermal or phase inventory.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RelaxingRegionalWater {
    equilibrium: RegionalWaterAcoustics,
    frozen: RegionalWaterAcoustics,
    arm: BulkRelaxationArm,
}

impl RelaxingRegionalWater {
    /// Retain the same interior coexistence state for both acoustic limits.
    ///
    /// Refuses nonpositive/nonfinite time, homogeneous/exact-endpoint states,
    /// unavailable original tangents, and endpoints without K_f > K_eq.
    /// The supplied time is a constitutive input, not an IF97-derived rate.
    pub fn try_new(
        source: RegionalWaterState,
        compression_relaxation_time_s: f64,
    ) -> Result<Self, MaterialError> {
        if !compression_relaxation_time_s.is_finite() || compression_relaxation_time_s <= 0.0 {
            return Err(MaterialError::Parameters {
                what: "water compression relaxation requires a finite positive time".into(),
            });
        }
        if !matches!(source.vapor_quality(), Some(x) if x > 0.0 && x < 1.0) {
            return Err(invalid(
                "water compression relaxation requires interior liquid/vapor coexistence",
            ));
        }
        let equilibrium = RegionalWaterAcoustics::try_new(
            source,
            CoexistenceAcousticClosure::HomogeneousEquilibrium,
        )?;
        let frozen =
            RegionalWaterAcoustics::try_new(source, CoexistenceAcousticClosure::FrozenPhaseMasses)?;
        let delta_bulk_modulus_pa = frozen.state().adiabatic_bulk_modulus_pa()
            - equilibrium.state().adiabatic_bulk_modulus_pa();
        if !delta_bulk_modulus_pa.is_finite() || delta_bulk_modulus_pa <= 0.0 {
            return Err(invalid(
                "water compression relaxation requires a positive frozen-to-equilibrium modulus gap",
            ));
        }
        Ok(Self {
            equilibrium,
            frozen,
            arm: BulkRelaxationArm {
                delta_bulk_modulus_pa,
                relaxation_time_s: compression_relaxation_time_s,
            },
        })
    }

    /// Complete unchanged thermal state, including its original regional roots.
    pub const fn source(&self) -> &RegionalWaterState {
        self.equilibrium.source()
    }

    /// Original equilibrium response, source tangent and source-fit defects.
    pub const fn equilibrium(&self) -> &RegionalWaterAcoustics {
        &self.equilibrium
    }

    /// Original frozen phase-mass response with the same thermal density.
    pub const fn frozen(&self) -> &RegionalWaterAcoustics {
        &self.frozen
    }

    /// Positive modulus increment and explicit fixed-compression time.
    pub const fn arm(&self) -> BulkRelaxationArm {
        self.arm
    }

    /// Coherent passive local response for finite nonnegative omega [rad/s].
    ///
    /// The convention is exp(i*k*x - i*omega*t): Im(k) >= 0 attenuates
    /// pressure amplitude; Im(K) <= 0 dissipates compression work. DC returns
    /// equilibrium compliance, zero wavenumber and finite limiting impedance.
    /// An unrepresentable response refuses without clamping or changing rho.
    pub fn response(&self, omega_rad_s: f64) -> Result<RelaxingWaterResponse, MaterialError> {
        let bulk_modulus_pa = C64::from_re(self.equilibrium.state().adiabatic_bulk_modulus_pa())
            + self.arm.response(omega_rad_s)?;
        let density_kg_m3 = self.equilibrium.state().density_kg_m3();
        let compressibility_per_pa = bulk_modulus_pa.recip();
        let density_susceptibility_kg_m3_pa = compressibility_per_pa.scale(density_kg_m3);
        // Form K/rho first: rho*K need not be representable even when Z is.
        let speed = bulk_modulus_pa.scale(density_kg_m3.recip()).sqrt();
        let wavenumber_rad_m = speed.recip().scale(omega_rad_s);
        let specific_impedance_pa_s_m = speed.scale(density_kg_m3);
        let positive_real = |z: C64| z.re.is_finite() && z.re > 0.0 && z.im.is_finite();
        if ![
            bulk_modulus_pa,
            compressibility_per_pa,
            density_susceptibility_kg_m3_pa,
            specific_impedance_pa_s_m,
        ]
        .into_iter()
        .all(positive_real)
            || bulk_modulus_pa.im > 0.0
            || compressibility_per_pa.im < 0.0
            || density_susceptibility_kg_m3_pa.im < 0.0
            || specific_impedance_pa_s_m.im > 0.0
            || !wavenumber_rad_m.re.is_finite()
            || !wavenumber_rad_m.im.is_finite()
            || wavenumber_rad_m.im < 0.0
            || (omega_rad_s > 0.0 && wavenumber_rad_m.re <= 0.0)
        {
            return Err(invalid("water compression response is unrepresentable"));
        }
        Ok(RelaxingWaterResponse {
            bulk_modulus_pa,
            compressibility_per_pa,
            density_susceptibility_kg_m3_pa,
            wavenumber_rad_m,
            specific_impedance_pa_s_m,
            density_kg_m3,
        })
    }
}

/// Frequency response of the original real-density state; no complex EOS.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RelaxingWaterResponse {
    bulk_modulus_pa: C64,
    compressibility_per_pa: C64,
    density_susceptibility_kg_m3_pa: C64,
    wavenumber_rad_m: C64,
    specific_impedance_pa_s_m: C64,
    density_kg_m3: f64,
}

impl RelaxingWaterResponse {
    /// Pressure/compressive-strain response [Pa].
    pub const fn bulk_modulus_pa(&self) -> C64 {
        self.bulk_modulus_pa
    }
    /// Compressive-strain/pressure response [1/Pa].
    pub const fn compressibility_per_pa(&self) -> C64 {
        self.compressibility_per_pa
    }
    /// Density-perturbation/pressure response for fixed total mass [kg/(m^3 Pa)].
    pub const fn density_susceptibility_kg_m3_pa(&self) -> C64 {
        self.density_susceptibility_kg_m3_pa
    }
    /// Forward wavenumber [rad/m]; imaginary part is amplitude loss [Np/m].
    pub const fn wavenumber_rad_m(&self) -> C64 {
        self.wavenumber_rad_m
    }
    /// Pressure/particle-velocity impedance [Pa s/m], also finite at DC.
    pub const fn specific_impedance_pa_s_m(&self) -> C64 {
        self.specific_impedance_pa_s_m
    }
    /// Unmodified real thermal density [kg/m^3], independent of frequency.
    pub const fn density_kg_m3(&self) -> f64 {
        self.density_kg_m3
    }
}

fn invalid(what: &str) -> MaterialError {
    MaterialError::State { what: what.into() }
}

#[cfg(test)]
mod tests;
