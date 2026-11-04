//! Directional pressure from the same full flanged-mouth radial basis.
//!
//! For the original modal volume flows U_n, the half-space Rayleigh integral
//! has far-field amplitude
//! F(mu) = -i rho omega/pi sum_n U_n sign(J0(alpha_n)) D_n(ka sqrt(1-mu^2)),
//! where p(r,mu) = F(mu) exp(i k r)/r + O(1/r^2), exp(-i omega t).
//! The original radiation matrix's analytic Hankel transforms, removable
//! singularity limits and signed normalization are reused without a surface
//! discretization. Every retained mode participates coherently; only U_0 is
//! net volume flow, and zero-net-flow modes still radiate off axis.
//!
//! This is a homogeneous nondispersive receiving half-space and an infinite
//! rigid baffle. It is a far-field observation, not a near-field Rayleigh
//! quadrature, an atmospheric propagation law or an additional radiation load.

use super::{FLANGED_MAX_KA, FLANGED_MAX_MODES, RULES, fill_d};
use crate::DuctError;
use crate::fluid::FluidDuctError;
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_math::{c64::C64, dd::Dd, det};

/// One coherent outgoing far-field sample, in peak-amplitude phasors.
#[derive(Debug, Clone, Copy)]
pub struct FlangedFarField {
    /// F(theta) [Pa m], referenced to the circular mouth centre.
    /// The leading outgoing pressure is F exp(i k r)/r.
    pub amplitude_pa_m: C64,
    /// Mean power per unit solid angle [W/sr], |F|^2/(2 rho c).
    /// Integrating over the outward hemisphere gives the full modal radiation
    /// power, including interference between different mouth modes.
    pub power_per_steradian_w: f64,
}

/// Original-material, frequency-bound flanged-mouth observation.
///
/// Preparation samples the receiving material once and shares the original
/// certified radial-root table. Directional queries allocate nothing and take
/// O(modes) work. The supplied flows must use the same physical mouth radius,
/// frequency and signed basis as the duct solve; no condensed plane-mode load
/// or independent piston is substituted for the retained higher modes.
#[derive(Debug, Clone)]
pub struct FlangedModalObserver {
    exterior: FluidAcousticState,
    wavenumber_per_m: f64,
    ka: f64,
    pressure_scale: f64,
    roots: &'static [(f64, f64)],
}

impl FlangedModalObserver {
    /// Prepare the receiving medium for the original circular mouth.
    ///
    /// The modes are psi_n(r/a) = J0(alpha_n r/a)/|J0(alpha_n)|,
    /// U_n = integral(v_out psi_n dA), including the plane mode n=0.
    /// Geometry, frequency and retained-mode admission match the original
    /// full flanged radiation matrix. No interior-medium property enters.
    ///
    /// # Errors
    /// Original material refusals, invalid geometry/frequency/mode count, the
    /// shared exterior ka limit, or an unrepresentable physical pressure scale.
    pub fn from_material<M: FluidAcoustics + ?Sized>(
        exterior: &M,
        radius_m: f64,
        omega: f64,
        modes: usize,
    ) -> Result<Self, FluidDuctError> {
        let exterior = exterior.acoustic_state()?;
        if !(radius_m > 0.0 && radius_m.is_finite() && omega > 0.0 && omega.is_finite()) {
            return Err(DuctError::BadParameter {
                what: "flanged observation needs positive finite radius and frequency",
            }
            .into());
        }
        if !(1..=FLANGED_MAX_MODES).contains(&modes) {
            return Err(DuctError::BadParameter {
                what: "flanged observation mode count outside 1..=FLANGED_MAX_MODES",
            }
            .into());
        }
        let wavenumber_per_m = omega / exterior.sound_speed_m_s();
        let ka = wavenumber_per_m * radius_m;
        if !(ka <= FLANGED_MAX_KA) {
            return Err(DuctError::RadiationKaTooLarge { ka }.into());
        }
        let volume_impedance = exterior.characteristic_impedance_pa_s_m()
            / (core::f64::consts::PI * radius_m * radius_m);
        let pressure_scale = (omega / core::f64::consts::PI) * exterior.density_kg_m3();
        if !(volume_impedance > 0.0
            && volume_impedance.is_finite()
            && pressure_scale > 0.0
            && pressure_scale.is_finite())
        {
            return Err(DuctError::BadParameter {
                what: "flanged observation physical scales are unrepresentable",
            }
            .into());
        }
        let rules = RULES.as_ref().map_err(Clone::clone)?;
        Ok(Self {
            exterior,
            wavenumber_per_m,
            ka,
            pressure_scale,
            roots: &rules.roots[..modes - 1],
        })
    }

    /// Coherent amplitude and angular power in the outward hemisphere.
    ///
    /// cos_theta is the cosine from the outward mouth axis: one is on axis
    /// and zero is grazing the baffle. Aperture modes need not propagate inside
    /// the duct to contribute to this exterior field. The amplitude is the
    /// exact far-field limit of the admitted finite radial expansion.
    ///
    /// # Errors
    /// A mismatched or nonfinite flow vector, a direction outside the outward
    /// hemisphere, or unrepresentable pressure/power.
    pub fn far_field(
        &self,
        outward_flow_m3_s: &[C64],
        cos_theta: f64,
    ) -> Result<FlangedFarField, DuctError> {
        let modes = self.roots.len() + 1;
        if outward_flow_m3_s.len() != modes
            || outward_flow_m3_s
                .iter()
                .any(|z| !z.re.is_finite() || !z.im.is_finite())
            || !cos_theta.is_finite()
            || !(0.0..=1.0).contains(&cos_theta)
        {
            return Err(DuctError::BadParameter {
                what: "flanged observation requires original finite modal flows and outward direction",
            });
        }
        let x = self.ka * det::sqrt((1.0 - cos_theta) * (1.0 + cos_theta));
        let mut d = [0.0; FLANGED_MAX_MODES];
        fill_d(
            x,
            fs_math::special::bessel_j1(x),
            self.roots,
            &mut d[..modes],
        );
        let mut re = Dd::ZERO;
        let mut im = Dd::ZERO;
        for (n, (&flow, &hankel)) in outward_flow_m3_s.iter().zip(&d).enumerate() {
            let sign = if n == 0 || self.roots[n - 1].1 > 0.0 {
                1.0
            } else {
                -1.0
            };
            let weight = Dd::from_f64(sign * hankel);
            re = re + Dd::from_f64(flow.re) * weight;
            im = im + Dd::from_f64(flow.im) * weight;
        }
        let amplitude = C64::new(im.to_f64(), -re.to_f64()).scale(self.pressure_scale);
        let normalized =
            amplitude.abs() / det::sqrt(self.exterior.characteristic_impedance_pa_s_m());
        let power = 0.5 * normalized * normalized;
        if !amplitude.re.is_finite() || !amplitude.im.is_finite() || !power.is_finite() {
            return Err(DuctError::BadParameter {
                what: "flanged far-field pressure or power is unrepresentable",
            });
        }
        Ok(FlangedFarField {
            amplitude_pa_m: amplitude,
            power_per_steradian_w: power,
        })
    }

    /// Leading far-zone pressure [Pa] at distance r from the mouth centre.
    ///
    /// This applies the original outward exp(i k r)/r propagation once.
    /// The caller must place the receiver in the far zone: both a/r and
    /// k a^2/(2r) must be small at the required accuracy. Passing a finite
    /// positive distance alone does not establish that approximation.
    /// Near-field curvature and evanescent exterior pressure are not included.
    ///
    /// # Errors
    /// As [Self::far_field], plus an invalid distance or propagation phase.
    pub fn far_field_pressure(
        &self,
        outward_flow_m3_s: &[C64],
        cos_theta: f64,
        distance_m: f64,
    ) -> Result<C64, DuctError> {
        let phase = self.wavenumber_per_m * distance_m;
        if !(distance_m > 0.0 && distance_m.is_finite() && phase.is_finite()) {
            return Err(DuctError::BadParameter {
                what: "flanged far-field pressure needs finite positive distance and phase",
            });
        }
        let amplitude = self.far_field(outward_flow_m3_s, cos_theta)?.amplitude_pa_m;
        let (sin, cos) = det::sin_cos(phase);
        let pressure = amplitude.scale(1.0 / distance_m) * C64::new(cos, sin);
        if !pressure.re.is_finite() || !pressure.im.is_finite() {
            return Err(DuctError::BadParameter {
                what: "flanged far-field pressure is unrepresentable",
            });
        }
        Ok(pressure)
    }
}

#[cfg(test)]
mod tests;
