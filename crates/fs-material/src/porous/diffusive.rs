//! Positive causal memory of the original rigid-frame JCA/JCAL fluid law.
//!
//! The original density law (Johnson, Koplik & Dashen, JFM 176, 1987,
//! doi:10.1017/S0022112087000727) gives
//! `Z(s)/rho_inf = s + d sqrt(1 + tau*s)`. The original compressibility
//! (Lafarge et al., JASA 102, 1997, doi:10.1121/1.419690, Eqs. 22--30) gives
//! `Y(s)/C_inf = s + (gamma-1) s F(s)`, with
//! `F(s) = A sqrt(1+B*s)/(s+A sqrt(1+B*s))`.
//!
//! The following positive measures follow from their negative-real-axis cuts:
//! `sqrt(1+tau*s)-1 = integral s/(s+r) sqrt(tau*r-1)/(pi*r) dr`,
//! over `r > 1/tau`. For `F`, there is one pole `-r0` below `-1/B`,
//! `r0 = 2*A/(sqrt((A*B)^2+4)+A*B)`, with positive residue
//! `h0 = 2*r0^2/(2*r0+A^2*B)`. Its remaining Stieltjes density is
//! `A*r*sqrt(B*r-1)/(pi*(r^2+A^2*(B*r-1)))`, for `r > 1/B`.
//! In particular the isolated thermal pole must not be dropped.
//!
//! Substitute `r = sec(theta)^2/tau` or `sec(theta)^2/B`, then apply
//! positive midpoint quadrature on `0 < theta < pi/2`. This yields positive
//! Foster arms `g*s/(s+p)`; no fitted residues, loss clipping, or artificial
//! tail load enter the model. With energy coordinates x,y each arm has
//! storage `(x^2+y^2)/2` and loss `(sqrt(g)*x-sqrt(p)*y)^2`.
//! Positive diffusive variables for rigid porous equivalent fluids and their
//! energy interpretation are also developed by Moufid et al., JCP 451 (2022),
//! 110888, doi:10.1016/j.jcp.2021.110888.
//!
//! This is a coefficient producer, not another time integrator. The finite
//! quadrature is checked against the ORIGINAL complex material law on 129
//! logarithmic frequencies, separately for complex response and real loss.
//! Those sampled maxima are not continuous-band error bounds. Signals outside
//! the declared source-supported band, rigid-frame assumptions, and frozen
//! fluid state have no accuracy or source qualification from this admission.

use super::PorousFluidLaw;
use crate::MaterialError;
use crate::fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState};
use fs_math::{c64::C64, det};

/// Bounded number of quadrature poles per continuum. Thermal relaxation also
/// retains its one isolated real pole. Consumers admit their own total state.
pub const MAX_POROUS_DIFFUSIVE_POLES: usize = 64;

/// One positive normalized resistance/conductance `g*s/(s+p)`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PorousDiffusiveTerm {
    /// Resistance divided by high-frequency inertia, or conductance divided
    /// by high-frequency compliance [1/s].
    pub gain_per_s: f64,
    /// Original constitutive relaxation rate [1/s].
    pub rate_per_s: f64,
}

/// Frozen source-bound finite-band realization of the original porous law.
/// Porosity already occurs in rho_inf and C_inf; use superficial volume flow
/// through the full specimen area, with no second porosity multiplier.
#[derive(Clone, Debug)]
pub struct PorousDiffusive {
    material: PorousFluidLaw,
    fluid: ThermoviscousFluidState,
    band_rad_s: [f64; 2],
    density_infinite: f64,
    bulk_infinite: f64,
    series_drag_per_s: f64,
    series: Vec<PorousDiffusiveTerm>,
    shunt: Vec<PorousDiffusiveTerm>,
    sampled_max_complex_error: f64,
    sampled_max_real_loss_error: f64,
}

fn invalid(what: &'static str) -> MaterialError {
    MaterialError::Parameters { what: what.into() }
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn finite(value: C64) -> bool {
    value.re.is_finite() && value.im.is_finite()
}
fn term(gain: f64, rate: f64) -> Result<PorousDiffusiveTerm, MaterialError> {
    if !positive(gain) || !positive(rate) {
        return Err(invalid(
            "porous diffusive pole or positive weight is unrepresentable",
        ));
    }
    Ok(PorousDiffusiveTerm {
        gain_per_s: gain,
        rate_per_s: rate,
    })
}

impl PorousDiffusive {
    /// Resolve the fluid ONCE, retain the original skeleton, and admit a
    /// positive causal approximation against that same material's complex law.
    /// `poles` controls each continuous measure; Cp=Cv omits all thermal arms.
    /// `maximum_relative_error` applies separately to both complex responses
    /// and their nonzero real losses. Weak thermal loss is not hidden behind
    /// the much larger reversible acoustic compliance.
    ///
    /// No state propagation, material heating or source-domain extension is
    /// performed. The band is angular frequency [rad/s], inside the original
    /// card's measured/identified frequency support.
    pub fn new<M: ThermoviscousFluid + ?Sized>(
        material: &PorousFluidLaw,
        fluid: &M,
        poles: usize,
        band_rad_s: [f64; 2],
        maximum_relative_error: f64,
    ) -> Result<Self, MaterialError> {
        if !(1..=MAX_POROUS_DIFFUSIVE_POLES).contains(&poles)
            || !positive(band_rad_s[0])
            || !positive(band_rad_s[1])
            || band_rad_s[1] < band_rad_s[0]
            || !positive(maximum_relative_error)
            || maximum_relative_error >= 1.0
        {
            return Err(invalid(
                "porous diffusion needs bounded poles, positive ordered band and error below one",
            ));
        }
        let fluid = fluid.thermoviscous_state()?;
        let first = material.effective_fluid(&fluid, band_rad_s[0])?;
        material.effective_fluid(&fluid, band_rad_s[1])?;
        let card = material.reference_card();
        let acoustic = fluid.acoustic();
        let transport = fluid.transport();
        let rho = acoustic.density_kg_m3();
        let mu = transport.dynamic_viscosity_pa_s();
        let density_infinite = card.tortuosity * rho / card.porosity;
        let bulk_infinite = acoustic.adiabatic_bulk_modulus_pa() / card.porosity;
        let drag = first.flow_resistivity_pa_s_m2() / density_infinite;
        let length = density_infinite / first.flow_resistivity_pa_s_m2() / card.viscous_length_m;
        let tau = 4.0 * (mu / rho) * length * length;
        if [
            density_infinite,
            bulk_infinite,
            1.0 / bulk_infinite,
            drag,
            tau,
        ]
        .into_iter()
        .any(|value| !positive(value))
        {
            return Err(invalid(
                "porous high-frequency storage or viscous timescale is unrepresentable",
            ));
        }
        let mut series = Vec::with_capacity(poles);
        for i in 0..poles {
            let theta = core::f64::consts::FRAC_PI_2 * ((i as f64 + 0.5) / poles as f64);
            let (sin, cos) = det::sin_cos(theta);
            let tan = sin / cos;
            series.push(term(
                drag * (tan * tan / poles as f64),
                (1.0 / cos / cos) / tau,
            )?);
        }

        let gamma_excess = fluid.heat_capacity_ratio() - 1.0;
        let mut shunt = Vec::new();
        if gamma_excess > 0.0 {
            let k0 = card
                .thermal_permeability_m2
                .unwrap_or(card.porosity * card.thermal_length_m * card.thermal_length_m / 8.0);
            let pr = transport.prandtl();
            let a = mu * card.porosity / (rho * pr * k0);
            let length = k0 / (card.thermal_length_m * card.porosity);
            let b = 4.0 * (rho * pr / mu) * length * length;
            let q = a * b;
            if [a, b, q].into_iter().any(|value| !positive(value)) {
                return Err(invalid(
                    "porous thermal constitutive timescales are unrepresentable",
                ));
            }
            // Cancellation-free positive root of t*t + q*t = 1.
            let root = det::hypot(q, 2.0);
            let t = 1.0 / (0.5 * root + 0.5 * q);
            let rate = a * t;
            let residue = rate * ((2.0 * t) / root);
            shunt.reserve(poles + 1);
            shunt.push(term(gamma_excess * residue, rate)?);
            for i in 0..poles {
                let theta = core::f64::consts::FRAC_PI_2 * ((i as f64 + 0.5) / poles as f64);
                let (sin, cos) = det::sin_cos(theta);
                let tan = sin / cos;
                // Divide the sec^4 numerator and denominator by sec^4.
                let scaled = q * sin * cos;
                let divisor = 1.0 + scaled * scaled;
                let gain = gamma_excess * a * (tan * tan / poles as f64) / divisor;
                shunt.push(term(gain, (1.0 / cos / cos) / b)?);
            }
        }
        let mut result = Self {
            material: material.clone(),
            fluid,
            band_rad_s,
            density_infinite,
            bulk_infinite,
            series_drag_per_s: drag,
            series,
            shunt,
            sampled_max_complex_error: 0.0,
            sampled_max_real_loss_error: 0.0,
        };
        let log_low = det::ln(band_rad_s[0]);
        let log_high = det::ln(band_rad_s[1]);
        for i in 0..=128 {
            let omega = match i {
                0 => band_rad_s[0],
                128 => band_rad_s[1],
                _ => det::exp(log_low + (log_high - log_low) * (i as f64 / 128.0)),
            };
            let original = material.effective_fluid(&result.fluid, omega)?;
            let s = C64::new(0.0, omega);
            let z = s * original
                .dynamic_density_kg_m3()
                .scale(1.0 / density_infinite);
            let y = s / original
                .dynamic_bulk_modulus_pa()
                .scale(1.0 / bulk_infinite);
            let (realized_z, realized_y) = result.normalized_at(omega);
            for (approximate, exact) in [(realized_z, z), (realized_y, y)] {
                if !finite(approximate) || !finite(exact) || exact.re < 0.0 {
                    return Err(invalid(
                        "porous diffusion reference response is unrepresentable",
                    ));
                }
                let error = (approximate - exact).abs() / exact.abs();
                let loss_error = if exact.re > 0.0 {
                    (approximate.re - exact.re).abs() / exact.re
                } else if approximate.re == 0.0 {
                    0.0
                } else {
                    return Err(invalid(
                        "porous diffusion creates loss in an inactive source",
                    ));
                };
                if !error.is_finite()
                    || !loss_error.is_finite()
                    || error > maximum_relative_error
                    || loss_error > maximum_relative_error
                {
                    return Err(invalid(
                        "porous diffusive quadrature is unresolved on the declared band",
                    ));
                }
                result.sampled_max_complex_error = result.sampled_max_complex_error.max(error);
                result.sampled_max_real_loss_error =
                    result.sampled_max_real_loss_error.max(loss_error);
            }
        }
        Ok(result)
    }

    /// Original skeleton and characterization-fluid viscosity.
    pub const fn material(&self) -> &PorousFluidLaw {
        &self.material
    }
    /// Original coherent fluid snapshot resolved once at admission.
    pub const fn fluid(&self) -> ThermoviscousFluidState {
        self.fluid
    }
    /// Admitted original source interval [rad/s].
    pub const fn band_rad_s(&self) -> [f64; 2] {
        self.band_rad_s
    }
    /// High-frequency effective inertia density alpha_inf*rho/phi [kg/m3].
    pub const fn density_infinite_kg_m3(&self) -> f64 {
        self.density_infinite
    }
    /// High-frequency effective modulus K_s/phi [Pa].
    pub const fn bulk_infinite_pa(&self) -> f64 {
        self.bulk_infinite
    }
    /// Exact static Darcy drag divided by high-frequency inertia [1/s].
    pub const fn series_drag_per_s(&self) -> f64 {
        self.series_drag_per_s
    }
    /// Positive viscous-memory coefficients, excluding the static drag.
    pub fn series_terms(&self) -> &[PorousDiffusiveTerm] {
        &self.series
    }
    /// Positive thermal-memory coefficients, isolated real pole first.
    pub fn shunt_terms(&self) -> &[PorousDiffusiveTerm] {
        &self.shunt
    }
    /// Maximum relative complex response discrepancy on the admission samples.
    /// This is not an interval-enclosed continuous-band error bound.
    pub const fn sampled_max_complex_error(&self) -> f64 {
        self.sampled_max_complex_error
    }
    /// Maximum relative real-loss discrepancy on the admission samples.
    pub const fn sampled_max_real_loss_error(&self) -> f64 {
        self.sampled_max_real_loss_error
    }

    fn normalized_at(&self, omega: f64) -> (C64, C64) {
        let s = C64::new(0.0, omega);
        let mut series = s + C64::from_re(self.series_drag_per_s);
        let mut shunt = s;
        for (terms, response) in [(&self.series, &mut series), (&self.shunt, &mut shunt)] {
            for term in terms {
                *response =
                    *response + (s / (s + C64::from_re(term.rate_per_s))).scale(term.gain_per_s);
            }
        }
        (series, shunt)
    }

    /// Realized normalized series impedance and shunt admittance, including
    /// their high-frequency acoustic `s` terms, in exp(+i omega t) convention.
    /// The band check rejects queries not covered by this approximation.
    pub fn normalized_response(&self, omega: f64) -> Result<(C64, C64), MaterialError> {
        if !omega.is_finite() || omega < self.band_rad_s[0] || omega > self.band_rad_s[1] {
            return Err(invalid(
                "porous diffusive query is outside its admitted band",
            ));
        }
        let response = self.normalized_at(omega);
        if !finite(response.0) || !finite(response.1) {
            return Err(invalid("porous diffusive response is unrepresentable"));
        }
        Ok(response)
    }

    /// Realized effective density [kg/m3] and bulk modulus [Pa], using the
    /// same porosity and exp(+i omega t) convention as the original material.
    pub fn effective_density_bulk(&self, omega: f64) -> Result<(C64, C64), MaterialError> {
        let (z, y) = self.normalized_response(omega)?;
        let s = C64::new(0.0, omega);
        let density = (z / s).scale(self.density_infinite);
        let bulk = (s / y).scale(self.bulk_infinite);
        if !finite(density) || !finite(bulk) {
            return Err(invalid(
                "porous diffusive effective density or modulus is unrepresentable",
            ));
        }
        Ok((density, bulk))
    }
}

#[cfg(test)]
mod tests;
