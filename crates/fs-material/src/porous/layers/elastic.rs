//! Refraction and P/S mode conversion into an isotropic elastic half-space.
//!
//! The original admitted elastic tensor supplies density and longitudinal/shear
//! speeds. The incident fluid and finite porous/fluid layers are the existing
//! prepared owners. The receiver is a frozen, unstressed, lossless small-strain
//! solid; no attenuation, poroelastic frame motion or anisotropic projection is
//! inferred. Static source moduli retain their original response regime.
//!
//! For exp(+i omega t), outgoing waves use exp(-i k_n z), with real positive
//! or imaginary negative k_n. Tangential slowness is shared by the fluid, P and
//! SV waves. At the interface normal velocity and normal traction are continuous
//! and shear traction vanishes; tangential fluid/solid velocity may slip.
//!
//! With P/S velocity amplitudes A/B and s=sin(theta_S), shear freedom is
//! `2 s cos(theta_P) A + (1-2s^2) B = 0`. Its normalized null vector supplies
//! the pressure/velocity boundary pair without dividing by a critical cosine,
//! the zero of cos(2 theta_S) or the surface impedance. Both outgoing modal
//! powers are evaluated independently, including the P-evanescent/S-radiating
//! interval; an evanescent mode retains its field but carries no normal power.
//!
//! Reference: Martire et al., SPECFEM2D-DG, Appendix A, Eqs. (A1)-(A5):
//! <https://doi.org/10.1093/gji/ggab308>. Here amplitudes are particle velocity
//! relative to incident particle velocity, and pressure is positive compression.

pub mod plate;
pub use plate::{ElasticPlateLayers, LayeredElasticPlateTransmission};
pub mod stack;
pub use stack::{ElasticPlateStack, ElasticPlateStackLayer};

use super::transmission::transmit_state;
use super::{
    Layer, RigidBackedLayers, bad, divided, finite, magnitude_scale, normal_wavenumber_ratio,
};
use crate::MaterialError;
use crate::elastic::acoustic::ElasticAcousticState;
use crate::fluid_acoustics::thermoviscous::ThermoviscousFluid;
use fs_math::{c64::C64, det};

/// Original equivalent-fluid layers terminated by an explicitly isotropic
/// elastic source. Preparation evaluates the two Christoffel speeds once;
/// angle queries allocate nothing and retain the original material identity.
pub struct ElasticTransmittingLayers<'a> {
    prepared: RigidBackedLayers<'a>,
    receiver: &'a ElasticAcousticState,
    longitudinal_speed: f64,
    shear_speed: f64,
    longitudinal_ratio: f64,
    shear_ratio: f64,
    impedance_ratio: f64,
    lossless: bool,
}

/// Coherent fields at the entrance and solid surface, and normal energy flux.
/// The P/S velocity coefficients multiply incident `pressure/(rho*c)`.
/// In the incidence plane their outgoing polarizations are respectively
/// `(sin(theta_P), cos(theta_P))` and `(cos(theta_S), -sin(theta_S))`,
/// continued analytically for evanescent waves. The solid surface pressure
/// coefficient is minus normal stress divided by incident pressure.
#[derive(Clone, Copy, Debug)]
pub struct LayeredElasticTransmission {
    pub pressure_reflection: C64,
    pub surface_pressure_transmission: C64,
    pub longitudinal_velocity_transmission: C64,
    pub shear_velocity_transmission: C64,
    pub longitudinal_normal_wavenumber_rad_m: C64,
    pub shear_normal_wavenumber_rad_m: C64,
    pub reflected_fraction: f64,
    pub longitudinal_transmitted_fraction: f64,
    pub shear_transmitted_fraction: f64,
    /// Loss in the finite porous layers; the elastic half-space is lossless.
    pub absorbed_fraction: f64,
}

impl<'a> ElasticTransmittingLayers<'a> {
    /// The receiver must originate from `ElasticAcousticState::from_isotropic`
    /// or `from_bulk_and_shear`, including the existing source-selected K/G
    /// adapter. Arbitrary Mandel/cubic/orthotropic states are refused rather
    /// than fitted to an isotropic surrogate. The source strain ceiling and
    /// original identity remain available on the borrowed state; coefficients
    /// describe infinitesimal waves, not admission of a finite drive amplitude.
    pub fn new(
        incident: &'a dyn ThermoviscousFluid,
        receiver: &'a ElasticAcousticState,
        omega: f64,
        layers: &[Layer<'a>],
    ) -> Result<Self, MaterialError> {
        if !receiver.is_explicitly_isotropic() {
            return Err(bad(
                "elastic layer receiver requires an explicitly isotropic source",
            ));
        }
        let prepared = RigidBackedLayers::new(incident, omega, layers)?;
        let modes = receiver.directional_waves([0.0, 0.0, 1.0])?.modes;
        let longitudinal_speed = modes[2].phase_speed_m_s;
        let shear_speed = modes[0].phase_speed_m_s;
        let incident_state = prepared.fluid().acoustic();
        let longitudinal_ratio = incident_state.sound_speed_m_s() / longitudinal_speed;
        let shear_ratio = incident_state.sound_speed_m_s() / shear_speed;
        let impedance_ratio =
            incident_state.characteristic_impedance_pa_s_m() / modes[2].impedance_pa_s_m;
        if [longitudinal_ratio, shear_ratio, impedance_ratio]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(bad("elastic receiving contrast is unrepresentable"));
        }
        let lossless = prepared
            .layers
            .iter()
            .all(|layer| layer.porous.is_none() || layer.source.thickness_m() == 0.0);
        Ok(Self {
            prepared,
            receiver,
            longitudinal_speed,
            shear_speed,
            longitudinal_ratio,
            shear_ratio,
            impedance_ratio,
            lossless,
        })
    }

    pub const fn prepared_layers(&self) -> &RigidBackedLayers<'a> {
        &self.prepared
    }

    pub const fn receiving_state(&self) -> &'a ElasticAcousticState {
        self.receiver
    }

    /// Longitudinal and shear phase speeds, in that order, from the same
    /// original elastic tensor and density used by the mechanical owner.
    pub const fn receiving_phase_speeds_m_s(&self) -> [f64; 2] {
        [self.longitudinal_speed, self.shear_speed]
    }

    /// Coherent reflection, refraction, mode conversion and finite-layer loss.
    /// Critical/grazing receiver modes use their finite traction limit; exact
    /// grazing incidence in the *incident* fluid is outside the normal-power
    /// normalization, as it is for the existing fluid transmission owner.
    pub fn scattering(&self, cosine: f64) -> Result<LayeredElasticTransmission, MaterialError> {
        if !cosine.is_finite() || cosine <= 0.0 || cosine > 1.0 {
            return Err(bad(
                "elastic layer transmission needs incidence cosine in (0,1]",
            ));
        }
        let incident_speed = self.prepared.fluid().acoustic().sound_speed_m_s();
        let qp = normal_wavenumber_ratio(
            self.longitudinal_ratio,
            cosine,
            incident_speed,
            self.longitudinal_speed,
        );
        let qs =
            normal_wavenumber_ratio(self.shear_ratio, cosine, incident_speed, self.shear_speed);
        let p_kn = qp.scale(self.prepared.k0);
        let s_kn = qs.scale(self.prepared.k0);
        if !finite(p_kn)
            || !finite(s_kn)
            || (qp != C64::ZERO && p_kn == C64::ZERO)
            || (qs != C64::ZERO && s_kn == C64::ZERO)
        {
            return Err(bad("elastic receiving wavenumber is unrepresentable"));
        }
        let cp = divided(qp, self.longitudinal_ratio);
        let cs = divided(qs, self.shear_ratio);
        let sine = det::sqrt((1.0 - cosine) * (1.0 + cosine));
        let ss = sine / self.shear_ratio;
        let a = (-2.0 * ss).mul_add(ss, 1.0);
        let b = cp.scale(2.0 * ss);
        if !a.is_finite() || !finite(b) || !finite(cp) || !finite(cs) {
            return Err(bad("elastic refraction direction is unrepresentable"));
        }
        let scale = a.abs().max(magnitude_scale(b));
        let (p_mode, s_mode) = if scale == 0.0 {
            // Simultaneous P grazing and zero shear traction (lambda=0).
            // The grazing P amplitude is unconstrained and carries no normal
            // flux; the outgoing SV mode alone supplies the finite load.
            (C64::ZERO, -C64::ONE)
        } else {
            (C64::from_re(a / scale), -divided(b, scale))
        };
        let speed_ratio = self.shear_speed / self.longitudinal_speed;
        // Normalize solid velocity by rho_s*c_P. Keep pressure AND velocity;
        // neither the Rayleigh pressure zero nor a critical velocity zero is
        // divided out. Every finite layer sees the same boundary work pair.
        let pressure = p_mode.scale(a) - (cs * s_mode).scale(2.0 * speed_ratio * ss);
        let normal_velocity = cp * p_mode - s_mode.scale(ss);
        let normal_impedance_ratio = self.impedance_ratio / cosine;
        if !normal_impedance_ratio.is_finite() || normal_impedance_ratio <= 0.0 {
            return Err(bad("elastic receiving normal impedance is unrepresentable"));
        }
        let velocity = normal_velocity.scale(normal_impedance_ratio);
        let state = transmit_state(&self.prepared, cosine, pressure, velocity)?;
        let log_ratio = det::ln(self.impedance_ratio);
        let longitudinal = if cp.re == 0.0 {
            0.0
        } else {
            state.power(p_mode, log_ratio - det::ln(cosine) + det::ln(cp.re))
        };
        let shear = if cs.re == 0.0 {
            0.0
        } else {
            state.power(
                s_mode,
                log_ratio - det::ln(cosine) + det::ln(speed_ratio) + det::ln(cs.re),
            )
        };
        let reflected = state.reflection.norm_sq();
        let absorbed = if self.lossless {
            0.0
        } else {
            state.entering - longitudinal - shear
        };
        let tolerance = 1024.0 * f64::EPSILON;
        if !longitudinal.is_finite()
            || longitudinal < 0.0
            || !shear.is_finite()
            || shear < 0.0
            || !absorbed.is_finite()
            || absorbed < -tolerance
            || (reflected + longitudinal + shear + absorbed - 1.0).abs() > tolerance
        {
            return Err(bad(
                "elastic layer transmission failed passive modal-power balance",
            ));
        }
        Ok(LayeredElasticTransmission {
            pressure_reflection: state.reflection,
            surface_pressure_transmission: state.coefficient(pressure)?,
            longitudinal_velocity_transmission: state.scaled_coefficient(p_mode, log_ratio)?,
            shear_velocity_transmission: state.scaled_coefficient(s_mode, log_ratio)?,
            longitudinal_normal_wavenumber_rad_m: p_kn,
            shear_normal_wavenumber_rad_m: s_kn,
            reflected_fraction: reflected,
            longitudinal_transmitted_fraction: longitudinal,
            shear_transmitted_fraction: shear,
            absorbed_fraction: absorbed,
        })
    }
}
