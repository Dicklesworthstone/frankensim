//! Linear gas and surrounding-liquid heat diffusion in the original bubble.
//!
//! For exp(-i omega t), a uniform gas pressure and an isothermal interface give
//! `kappa(omega) = gamma / [1 + 3(gamma-1)(z coth(z)-1)/z^2]`,
//! `z^2 = -i omega R^2/D`, and `D = k_g/(rho_g cp_g)`. This is the spherical
//! heat-equation solution, not a fitted polytropic exponent or damping rate.
//! See Commander and Prosperetti (1989), equations 19--28 (opposite Fourier
//! convention):
//! <https://pages.jh.edu/aprospe1/publications/PapersPublished/BubblyLiquids/Commander.pdf>.
//!
//! Static continuation keeps the original ISOTHERMAL gas inventory; the
//! material's gamma enters only the harmonic response. Gas coefficients refer
//! to the current equilibrium pressure and liquid/interface temperature.
//! [`IsothermalWallBubble`] imposes an interface temperature reservoir;
//! [`ConductingLiquidBubble`] instead solves the surrounding infinite liquid's
//! spherical heat equation with continuous interface temperature and heat flux.
//! Changing vapor inventory, latent heat, molecular relaxation and nonlinear
//! excursions are excluded. Continuum gas conduction, a spatially
//! uniform incident pressure and small fractional-radius amplitudes remain
//! caller applicability conditions. Gas compactness, negligible gas inertia
//! and a small vapor fraction are checked separately; liquid compactness alone
//! cannot establish uniform gas pressure.

use super::{BubbleLinearization, CavitationError, KellerMiksisBubble, positive};
use core::f64::consts::PI;
use fs_material::{
    MaterialError,
    fluid_transport::{FluidTransport, ThermalFluidTransportState},
    gas::GasState,
};
use fs_math::{c64::C64, det};

/// Conservative admission ceilings, not calibrated error bounds.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalBubbleLimits {
    /// Largest omega R/c_g, in (0, 0.1], for the uniform-pressure gas.
    pub maximum_gas_wavenumber_radius: f64,
    /// Largest rho_g/rho_l, in (0, 0.1], neglecting gas inertia.
    pub maximum_gas_liquid_density_ratio: f64,
    /// Largest p_v/(p_g+p_v), in (0, 0.1], neglecting vapor heat/mass transfer.
    pub maximum_vapor_pressure_fraction: f64,
}

/// Original material/bubble errors or an explicit thermal-model boundary.
#[derive(Clone, Debug, PartialEq)]
pub enum ThermalBubbleError {
    Material(MaterialError),
    Bubble(CavitationError),
    InvalidInput(&'static str),
    StateMismatch(&'static str),
    OutsideModel {
        quantity: &'static str,
        value: f64,
        maximum: f64,
    },
    NonFiniteDerived(&'static str),
}

impl core::fmt::Display for ThermalBubbleError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "thermal bubble calculation refused: {self:?}")
    }
}
impl std::error::Error for ThermalBubbleError {}
impl From<MaterialError> for ThermalBubbleError {
    fn from(error: MaterialError) -> Self {
        Self::Material(error)
    }
}
impl From<CavitationError> for ThermalBubbleError {
    fn from(error: CavitationError) -> Self {
        Self::Bubble(error)
    }
}

fn bound(quantity: &'static str, value: f64, maximum: f64) -> Result<(), ThermalBubbleError> {
    if !value.is_finite() || value < 0.0 || value > maximum {
        return Err(ThermalBubbleError::OutsideModel {
            quantity,
            value,
            maximum,
        });
    }
    Ok(())
}

fn finite_complex(value: C64) -> bool {
    value.re.is_finite() && value.im.is_finite()
}

impl ThermalBubbleLimits {
    /// Validate controls even when an enclosing population contains no gas.
    pub fn validate(self) -> Result<(), ThermalBubbleError> {
        for (name, value) in [
            ("gas compactness limit", self.maximum_gas_wavenumber_radius),
            (
                "gas/liquid density-ratio limit",
                self.maximum_gas_liquid_density_ratio,
            ),
            (
                "vapor pressure-fraction limit",
                self.maximum_vapor_pressure_fraction,
            ),
        ] {
            if !positive(value) || value > 0.1 {
                return Err(ThermalBubbleError::InvalidInput(name));
            }
        }
        Ok(())
    }
}

/// Thermal response about one stable equilibrium of the ORIGINAL gas law.
/// Construction prepares all state-dependent coefficients; queries allocate
/// nothing and neither evolve nor re-equilibrate the retained gas inventory.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsothermalWallBubble {
    bubble: KellerMiksisBubble,
    linearization: BubbleLinearization,
    gas: GasState,
    limits: ThermalBubbleLimits,
    gas_pressure_pa: f64,
    thermal_diffusivity_m2_s: f64,
    thermal_diffusion_time_s: f64,
    inertia_pa_s2: f64,
}

/// Harmonic coefficients for a fractional-radius amplitude eta = delta R/R.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalBubbleResponse {
    /// Complex gas exponent: 1 at DC, gamma at high thermal frequency.
    pub effective_polytropic_exponent: C64,
    /// omega R^2/D, using cp-based gas thermal diffusivity.
    pub thermal_diffusion_parameter: f64,
    /// omega R/c_g, independently checked against uniform gas pressure.
    pub gas_wavenumber_radius: f64,
    /// Native KM numerator/denominator [s^2]. Fractional-radius susceptibility
    /// is -pressure_transfer_s2/(rho_l R^2), with exp(-i omega t).
    pub pressure_transfer_s2: C64,
    /// 3 p_g kappa - 2 sigma/R - i 4 mu_l omega [Pa].
    pub wall_stiffness_pa: C64,
    /// Mean heat-conduction dissipation divided by |eta|^2 [W]. This includes
    /// both gas and liquid for [`ConductingLiquidBubble`], and gas alone for
    /// [`IsothermalWallBubble`]. It excludes viscosity and outgoing acoustic
    /// power, and is zero at DC.
    pub thermal_power_per_fractional_radius_squared_w: f64,
}

impl IsothermalWallBubble {
    /// Retain the original bubble and a coherent dry-gas material snapshot.
    /// Temperature must equal the liquid's; gas pressure must match p_g(R)
    /// within 64 machine epsilons. No reference radius, mass or EOS is reset.
    pub fn new(
        bubble: KellerMiksisBubble,
        equilibrium_radius_m: f64,
        gas: GasState,
        limits: ThermalBubbleLimits,
    ) -> Result<Self, ThermalBubbleError> {
        limits.validate()?;
        // This is the material owner's coherent gas admission and its original
        // overflow-safe cp-based diffusivity, not a second transport closure.
        let gas_transport = gas.transport_state()?;
        if bubble.gas().polytropic_exponent() != 1.0 {
            return Err(ThermalBubbleError::InvalidInput(
                "isothermal static gas inventory",
            ));
        }
        if gas.water_mole_fraction != 0.0 {
            return Err(ThermalBubbleError::InvalidInput(
                "noncondensable gas contains water",
            ));
        }
        let liquid = bubble.liquid();
        if gas.temperature != liquid.temperature_k() {
            return Err(ThermalBubbleError::StateMismatch(
                "gas and interface temperature",
            ));
        }
        let linearization = bubble.linearization_at_radius(equilibrium_radius_m)?;
        // Static wall balance makes the local jump zero even after a large
        // ambient shift. An explicit-radius caller still must keep the frozen
        // liquid near its ORIGINAL operating point, as population continuation
        // does; repeated material preparation cannot reset that reference.
        if equilibrium_radius_m != bubble.gas().reference_radius_m() {
            let reference =
                super::inception::static_response(&bubble, bubble.gas().reference_radius_m())?;
            let pressure_ratio = (linearization.equilibrium_absolute_pressure_pa
                - reference.equilibrium_absolute_pressure_pa)
                .abs()
                / liquid.liquid_density_kg_m3()
                / liquid.liquid_sound_speed_m_s()
                / liquid.liquid_sound_speed_m_s();
            let maximum_pressure_ratio = bubble.compressibility_limits().maximum_pressure_ratio;
            if !pressure_ratio.is_finite() || pressure_ratio > maximum_pressure_ratio {
                return Err(CavitationError::PressureRatioExceeded {
                    pressure_ratio,
                    maximum_pressure_ratio,
                }
                .into());
            }
        }
        let gas_pressure_pa = bubble.gas().pressure_at_radius(equilibrium_radius_m)?;
        if (gas.pressure - gas_pressure_pa).abs()
            > 64.0 * f64::EPSILON * gas.pressure.max(gas_pressure_pa)
        {
            return Err(ThermalBubbleError::StateMismatch(
                "gas equilibrium partial pressure",
            ));
        }
        bound(
            "gas/liquid density ratio",
            gas.density / liquid.liquid_density_kg_m3(),
            limits.maximum_gas_liquid_density_ratio,
        )?;
        bound(
            "vapor pressure fraction",
            liquid.saturation_pressure_pa() / (gas_pressure_pa + liquid.saturation_pressure_pa()),
            limits.maximum_vapor_pressure_fraction,
        )?;
        let thermal_diffusivity_m2_s = gas_transport.thermal_diffusivity_m2_s();
        let diffusion_length = equilibrium_radius_m / det::sqrt(thermal_diffusivity_m2_s);
        let thermal_diffusion_time_s = diffusion_length * diffusion_length;
        let inertia_pa_s2 =
            liquid.liquid_density_kg_m3() * equilibrium_radius_m * equilibrium_radius_m;
        if !positive(thermal_diffusivity_m2_s)
            || !positive(thermal_diffusion_time_s)
            || !positive(inertia_pa_s2)
        {
            return Err(ThermalBubbleError::NonFiniteDerived(
                "gas heat-diffusion coefficients",
            ));
        }
        Ok(Self {
            bubble,
            linearization,
            gas,
            limits,
            gas_pressure_pa,
            thermal_diffusivity_m2_s,
            thermal_diffusion_time_s,
            inertia_pa_s2,
        })
    }

    pub const fn bubble(&self) -> KellerMiksisBubble {
        self.bubble
    }
    pub const fn linearization(&self) -> BubbleLinearization {
        self.linearization
    }
    pub const fn gas(&self) -> GasState {
        self.gas
    }
    pub const fn limits(&self) -> ThermalBubbleLimits {
        self.limits
    }
    pub const fn thermal_diffusivity_m2_s(&self) -> f64 {
        self.thermal_diffusivity_m2_s
    }

    /// Steady harmonic response at a nonnegative angular frequency [rad/s].
    /// An enclosing acoustic owner must also admit the local amplitude and
    /// incident-field/liquid compactness before using these coefficients.
    pub fn response(
        &self,
        angular_frequency_rad_s: f64,
    ) -> Result<ThermalBubbleResponse, ThermalBubbleError> {
        let omega = angular_frequency_rad_s;
        let (x, gas_wavenumber_radius) = self.frequency_parameters(omega)?;
        let kappa = effective_polytropic_exponent(self.gas.gamma, x);
        self.response_with_exponent(omega, x, gas_wavenumber_radius, kappa)
    }

    fn frequency_parameters(&self, omega: f64) -> Result<(f64, f64), ThermalBubbleError> {
        if !omega.is_finite() || omega < 0.0 {
            return Err(ThermalBubbleError::InvalidInput(
                "nonnegative angular frequency",
            ));
        }
        let radius = self.linearization.equilibrium_radius_m;
        let gas_wavenumber_radius = omega * (radius / self.gas.sound_speed);
        bound(
            "gas compactness",
            gas_wavenumber_radius,
            self.limits.maximum_gas_wavenumber_radius,
        )?;
        let x = omega * self.thermal_diffusion_time_s;
        if !x.is_finite() {
            return Err(ThermalBubbleError::NonFiniteDerived(
                "thermal diffusion parameter",
            ));
        }
        Ok((x, gas_wavenumber_radius))
    }

    // Both thermal boundary conditions use this SAME native KM tangent.
    // There is no second radial equation, resonance or radiation coefficient.
    fn response_with_exponent(
        &self,
        omega: f64,
        x: f64,
        gas_wavenumber_radius: f64,
        kappa: C64,
    ) -> Result<ThermalBubbleResponse, ThermalBubbleError> {
        let radius = self.linearization.equilibrium_radius_m;
        let liquid = self.bubble.liquid();
        let wall_stiffness_pa = C64::new(
            3.0 * self.gas_pressure_pa * kappa.re - 2.0 * liquid.surface_tension_n_m() / radius,
            3.0 * self.gas_pressure_pa * kappa.im - 4.0 * liquid.dynamic_viscosity_pa_s() * omega,
        );
        // Divide the entire native KM operator by its forcing derivative.
        // W includes viscous wall stress, so this retains the native viscous
        // acceleration multiplier 1+b_visc*tau exactly once. Unlike directly
        // dividing the expanded numerator/denominator, all loss terms have
        // the same sign; tiny radiation is not a difference of O(omega) terms.
        let pressure_transfer_s2 = if omega == 0.0 {
            C64::from_re(1.0 / self.linearization.stiffness_coefficient_per_s2)
        } else {
            let derivative = C64::new(1.0, -omega * self.linearization.pressure_derivative_time_s);
            let inertia = C64::from_re(omega * omega) / derivative;
            let impedance = wall_stiffness_pa.scale(1.0 / self.inertia_pa_s2) - inertia;
            impedance.recip()
        };
        let thermal_power_per_fractional_radius_squared_w = if omega == 0.0 {
            0.0
        } else {
            6.0 * PI * radius * radius * radius * self.gas_pressure_pa * (omega * -kappa.im)
        };
        if !finite_complex(kappa)
            || !finite_complex(wall_stiffness_pa)
            || !finite_complex(pressure_transfer_s2)
            || (pressure_transfer_s2.re == 0.0 && pressure_transfer_s2.im == 0.0)
            || !thermal_power_per_fractional_radius_squared_w.is_finite()
            || thermal_power_per_fractional_radius_squared_w < 0.0
        {
            return Err(ThermalBubbleError::NonFiniteDerived(
                "thermal bubble response",
            ));
        }
        Ok(ThermalBubbleResponse {
            effective_polytropic_exponent: kappa,
            thermal_diffusion_parameter: x,
            gas_wavenumber_radius,
            pressure_transfer_s2,
            wall_stiffness_pa,
            thermal_power_per_fractional_radius_squared_w,
        })
    }
}

/// Exact harmonic conduction in the gas and an infinite surrounding liquid.
///
/// The far-field liquid temperature is fixed. Temperature and conductive heat
/// flux are continuous at a thermally perfect interface; no fitted heat-transfer
/// coefficient or new material law is introduced. Liquid transport belongs to
/// the original equilibrium point and remains frozen. Thermal expansion of the
/// liquid, pressure-driven liquid temperature changes, interfacial thermal
/// resistance, phase change and interactions between neighboring temperature
/// fields are outside this single-bubble closure. A small vapor pressure
/// fraction alone does not certify negligible latent heat: the caller must
/// establish the noncondensable-gas approximation for the actual liquid.
/// These are the gas/liquid heat equations and interface conditions of
/// Prosperetti (1977), Sec. IV, restricted to spatially uniform gas pressure:
/// <https://pages.jh.edu/aprospe1/publications/PapersPublished/GasBubbles/JASA_Thermal_77.pdf>.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConductingLiquidBubble {
    base: IsothermalWallBubble,
    liquid: ThermalFluidTransportState,
    liquid_thermal_diffusion_time_s: f64,
}

/// Harmonic interface heat exchange and the independent conduction losses.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConductingLiquidBubbleResponse {
    /// Original KM response with the coupled gas/liquid complex exponent.
    pub bubble: ThermalBubbleResponse,
    /// omega R^2/D_l, using the same liquid's cp-based diffusivity.
    pub liquid_thermal_diffusion_parameter: f64,
    /// Interface temperature perturbation divided by eta = delta R/R [K].
    pub interface_temperature_per_fractional_radius_k: C64,
    /// Heat rate from gas into liquid divided by eta [W], with exp(-i omega t).
    /// This oscillatory heat rate is distinct from mean irreversible loss.
    pub outward_heat_per_fractional_radius_w: C64,
    /// Mean irreversible gas conduction loss divided by |eta|^2 [W].
    pub gas_thermal_power_per_fractional_radius_squared_w: f64,
    /// Mean irreversible liquid conduction loss divided by |eta|^2 [W].
    pub liquid_thermal_power_per_fractional_radius_squared_w: f64,
}

impl ConductingLiquidBubble {
    /// Retain the admitted original bubble and a coherent liquid thermal view.
    /// Temperature, density and dynamic viscosity must exactly match the SAME
    /// equilibrium liquid snapshot. Sources remain with the caller's material
    /// owner; raw numerical views do not create a material-data receipt.
    pub fn new(
        base: IsothermalWallBubble,
        liquid: ThermalFluidTransportState,
    ) -> Result<Self, ThermalBubbleError> {
        let original = base.bubble.liquid();
        let transport = liquid.transport();
        if liquid.temperature_k() != original.temperature_k()
            || transport.density_kg_m3() != original.liquid_density_kg_m3()
            || transport.dynamic_viscosity_pa_s() != original.dynamic_viscosity_pa_s()
        {
            return Err(ThermalBubbleError::StateMismatch(
                "surrounding liquid thermal state",
            ));
        }
        let diffusion_length = base.linearization.equilibrium_radius_m
            / det::sqrt(transport.thermal_diffusivity_m2_s());
        let liquid_thermal_diffusion_time_s = diffusion_length * diffusion_length;
        if !positive(liquid_thermal_diffusion_time_s) {
            return Err(ThermalBubbleError::NonFiniteDerived(
                "liquid heat-diffusion coefficients",
            ));
        }
        Ok(Self {
            base,
            liquid,
            liquid_thermal_diffusion_time_s,
        })
    }

    pub const fn base(&self) -> &IsothermalWallBubble {
        &self.base
    }

    pub const fn liquid(&self) -> ThermalFluidTransportState {
        self.liquid
    }

    /// Solve the two spherical heat equations at one admitted frequency.
    ///
    /// For exp(-i omega t), set z_j^2=-i omega R^2/D_j, B=z_g coth(z_g)-1,
    /// F=3B/z_g^2, Y_g=4 pi R k_g B, Y_l=4 pi R k_l(1+z_l). Continuity gives
    /// T_i/T_ad=Y_g/(Y_g+Y_l), and the original gas mass/energy balance gives
    /// kappa=gamma/[1+(gamma-1)F Y_l/(Y_g+Y_l)]. The liquid solution is
    /// T_l(r)=T_i (R/r) exp[-z_l(r/R-1)]. Queries are allocation-free.
    pub fn response(
        &self,
        angular_frequency_rad_s: f64,
    ) -> Result<ConductingLiquidBubbleResponse, ThermalBubbleError> {
        let omega = angular_frequency_rad_s;
        let (x, gas_wavenumber_radius) = self.base.frequency_parameters(omega)?;
        let liquid_x = omega * self.liquid_thermal_diffusion_time_s;
        if !liquid_x.is_finite() {
            return Err(ThermalBubbleError::NonFiniteDerived(
                "liquid thermal diffusion parameter",
            ));
        }
        let gas_f = spherical_heat_average(x);
        // Recover B from F, retaining Re(B)=O(x^2) without subtracting one
        // from z*coth(z). This is essential for gas heat loss near DC.
        let gas_b = C64::new((x / 3.0) * gas_f.im, -(x / 3.0) * gas_f.re);
        let liquid_a = det::sqrt(0.5 * liquid_x);
        // Remove the common 4*pi*R factor before dividing admittances.
        let gas_admittance = gas_b.scale(self.base.gas.thermal_conductivity);
        let liquid_admittance = C64::new(1.0 + liquid_a, -liquid_a)
            .scale(self.liquid.transport().thermal_conductivity_w_m_k());
        if !finite_complex(gas_admittance) || !finite_complex(liquid_admittance) {
            return Err(ThermalBubbleError::NonFiniteDerived(
                "spherical heat admittance",
            ));
        }
        // Scaling avoids overflowing the sum at high conductivity. Evaluate
        // BOTH shares directly: 1 minus a nearly unit share erases real loss.
        let scale = gas_admittance
            .re
            .max(-gas_admittance.im)
            .max(liquid_admittance.re)
            .max(-liquid_admittance.im);
        let gas_scaled = C64::new(gas_admittance.re / scale, gas_admittance.im / scale);
        let liquid_scaled = C64::new(liquid_admittance.re / scale, liquid_admittance.im / scale);
        let sum = gas_scaled + liquid_scaled;
        // A large share's ordinary division subtracts two O(sqrt(omega))
        // imaginary terms to recover an O(omega) phase. The smaller term
        // may already have disappeared when forming the denominator. Form
        // the cross determinant from the ORIGINAL two admittances instead.
        // Its opposite signs are shared exactly by the two fractions.
        let share_imaginary = (gas_scaled * liquid_scaled.conj()).im / sum.norm_sq();
        let interface_share = C64::new((gas_scaled / sum).re, share_imaginary);
        let gas_drop_share = C64::new((liquid_scaled / sum).re, -share_imaginary);
        let effective_f = gas_f * gas_drop_share;
        let gamma = self.base.gas.gamma;
        let kappa = C64::from_re(gamma) / (C64::ONE + effective_f.scale(gamma - 1.0));
        let bubble = self
            .base
            .response_with_exponent(omega, x, gas_wavenumber_radius, kappa)?;
        let temperature = self.liquid.temperature_k();
        let adiabatic_temperature = kappa.scale(-3.0 * temperature * ((gamma - 1.0) / gamma));
        let interface_temperature_per_fractional_radius_k = adiabatic_temperature * interface_share;
        let gas_drop = adiabatic_temperature * gas_drop_share;
        let area_length = 4.0 * PI * self.base.linearization.equilibrium_radius_m;
        // Either side gives the same flux. Use the gas side when its
        // admittance is smaller to avoid multiplying huge and tiny values.
        let outward_heat_per_fractional_radius_w = if gas_admittance.abs()
            <= liquid_admittance.abs()
        {
            (gas_admittance * gas_drop).scale(area_length)
        } else {
            (liquid_admittance * interface_temperature_per_fractional_radius_k).scale(area_length)
        };
        // Green's identity gives P_j=Re(Y_j)|temperature drop|^2/(2*T0).
        // Compute both nonnegative entropy integrals independently, without
        // subtracting liquid loss from total gas p-V work. The square-root
        // form also postpones squaring tiny temperature amplitudes.
        let gas_loss_amplitude =
            det::sqrt((0.5 * area_length / temperature) * gas_admittance.re) * gas_drop.abs();
        let liquid_loss_amplitude =
            det::sqrt((0.5 * area_length / temperature) * liquid_admittance.re)
                * interface_temperature_per_fractional_radius_k.abs();
        let gas_thermal_power_per_fractional_radius_squared_w =
            gas_loss_amplitude * gas_loss_amplitude;
        let liquid_thermal_power_per_fractional_radius_squared_w =
            liquid_loss_amplitude * liquid_loss_amplitude;
        if !finite_complex(interface_temperature_per_fractional_radius_k)
            || !finite_complex(outward_heat_per_fractional_radius_w)
            || !gas_thermal_power_per_fractional_radius_squared_w.is_finite()
            || !liquid_thermal_power_per_fractional_radius_squared_w.is_finite()
        {
            return Err(ThermalBubbleError::NonFiniteDerived(
                "coupled thermal bubble response",
            ));
        }
        Ok(ConductingLiquidBubbleResponse {
            bubble,
            liquid_thermal_diffusion_parameter: liquid_x,
            interface_temperature_per_fractional_radius_k,
            outward_heat_per_fractional_radius_w,
            gas_thermal_power_per_fractional_radius_squared_w,
            liquid_thermal_power_per_fractional_radius_squared_w,
        })
    }
}

// Positive real x only. The Bernoulli expansion resolves thermal losses at
// arbitrarily low frequencies without subtracting z*coth(z) from one. Its
// first omitted term at x=1/4 is below 5e-20. The middle branch evaluates coth
// with exp(-2 Re z); the high-frequency branch drops corrections below 1e-17.
fn effective_polytropic_exponent(gamma: f64, x: f64) -> C64 {
    if x == 0.0 {
        return C64::ONE;
    }
    let f = spherical_heat_average(x);
    C64::from_re(gamma) / (C64::ONE + f.scale(gamma - 1.0))
}

// F=3(z*coth(z)-1)/z^2 for z^2=-i*x. Shared by both boundary conditions.
fn spherical_heat_average(x: f64) -> C64 {
    if x <= 0.25 {
        let x2 = x * x;
        let re = 1.0
            + x2 * (-2.0 / 315.0
                + x2 * (2.0 / 31185.0
                    + x2 * (-4.0 / 6081075.0
                        + x2 * (6.7523539550426975e-9 - x2 * 6.931929779700787e-11))));
        let im = x
            * (1.0 / 15.0
                + x2 * (-1.0 / 1575.0
                    + x2 * (1382.0 / 212837625.0
                        + x2 * (-6.664382636993903e-8
                            + x2 * (6.841545361377655e-10 - x2 * 7.023512045947465e-12)))));
        C64::new(re, im)
    } else {
        let a = det::sqrt(0.5 * x);
        let scale = 1.5 / a;
        if a >= 20.0 {
            C64::new(scale, scale * (1.0 - 1.0 / a))
        } else {
            let e = det::exp(-2.0 * a);
            let (sin, cos) = det::sin_cos(a);
            let one_minus_e = -det::expm1(-2.0 * a);
            let denominator = one_minus_e * one_minus_e + 4.0 * e * sin * sin;
            let coth_re = one_minus_e * (1.0 + e) / denominator;
            let coth_im = 4.0 * e * sin * cos / denominator;
            C64::new(
                scale * (coth_re - coth_im),
                3.0 * a.mul_add(coth_re + coth_im, -1.0) / x,
            )
        }
    }
}

#[cfg(test)]
#[path = "thermal/tests.rs"]
mod tests;
