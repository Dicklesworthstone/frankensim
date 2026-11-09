//! Physical absorption of the original homogeneous acoustic mode.
//!
//! Its disturbance-energy flux is Re(t u* + T' q'*/T0)/2, with total
//! compressive traction t and Fourier heat flux q'=-kappa grad(T'). The
//! positive losses are eta |grad(u)|²/2 and kappa |grad(T')|²/(2 T0).
//! The latter is T0 times entropy production, not another caloric heat source.
//! See Myers (1991), doi:10.1017/S0022112091002434, and the thermoviscous
//! intensity/dissipation definitions in
//! <https://doc.comsol.com/6.4/doc/com.comsol.help.aco/aco_ug_thermo.09.74.html>.

use super::{BulkThermoviscousWave, invalid};
use crate::MaterialError;
use crate::fluid_acoustics::thermoviscous::absorption::BulkThermoviscousEquilibrium;
use fs_math::{c64::C64, dd::Dd, det};

mod coherent;
pub use coherent::{BulkCoherentCellAbsorption, BulkCoherentWaveSegment};

/// Prepared one-way bulk-wave fields and absorption at ONE original equilibrium.
/// Peak phasors use exp(i k x - i omega t), with x >= 0. No material is queried
/// during observation; constant-work field and cell queries allocate nothing.
/// This is the original acoustic eigenmode, without a reflected wave, independent
/// thermal mode, boundary condition, molecular law or mean-temperature evolution.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BulkWaveAbsorption {
    equilibrium: BulkThermoviscousEquilibrium,
    wave: BulkThermoviscousWave,
    omega_rad_s: f64,
    temperature_over_velocity: C64,
    longitudinal_viscosity: f64,
}

/// Local physical phasors and cycle-mean work per physical cross-sectional area.
/// Thermal availability flux is part of acoustic energy flux; the complex
/// conductive heat flux is a first-order oscillation, not mean deposited heat.
/// Both power fluxes are signed along +x; physical dissipation is nonnegative.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BulkWaveField {
    pub velocity_m_s: C64,
    /// Original longitudinal strain-rate phasor du/dx [1/s]. Retained before
    /// squaring so coherent frequencies can supply physical finite-time heat.
    pub velocity_gradient_per_s: C64,
    pub pressure_pa: C64,
    pub normal_traction_pa: C64,
    pub temperature_perturbation_k: C64,
    pub conductive_heat_flux_w_m2: C64,
    pub mechanical_power_flux_w_m2: f64,
    pub thermal_availability_flux_w_m2: f64,
    /// Direct irreversible shear-plus-bulk viscous heating [W/m³].
    pub viscous_heating_w_m3: f64,
    /// Conductive entropy generation [W/(m³ K)], distinct from caloric heat.
    pub conductive_entropy_w_m3_k: f64,
    /// T0 times conductive entropy generation [W/m³].
    pub thermal_availability_loss_w_m3: f64,
}
impl BulkWaveField {
    /// Full disturbance-energy flux, including the original thermal port.
    pub fn acoustic_power_flux_w_m2(self) -> f64 {
        self.mechanical_power_flux_w_m2 + self.thermal_availability_flux_w_m2
    }
}

/// Exact spatial integrals within the admitted linear constant-medium model.
/// Only viscous_heat_w is a direct caloric heat source. Thermal entropy and
/// its availability loss stay separate; a mean thermal solver owns heat transport.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BulkWaveCellAbsorption {
    pub incoming_acoustic_power_w: f64,
    pub outgoing_acoustic_power_w: f64,
    pub viscous_heat_w: f64,
    pub conductive_entropy_w_k: f64,
    pub thermal_availability_loss_w: f64,
    pub acoustic_power_loss_w: f64,
    /// Incoming minus outgoing minus independently integrated physical loss.
    /// Retained even when the true attenuation is below endpoint-power roundoff.
    pub balance_residual_w: f64,
}

impl BulkThermoviscousEquilibrium {
    /// Prepare physical absorption of the native traveling acoustic mode.
    /// The SAME signed expansion and temperature determine the thermal field;
    /// no ideal-gas inference from Cp/Cv is made. Positive finite frequency is
    /// required for cycle means. Original branch, source-band and continuum
    /// limits remain; this is not a finite-speed or experimentally qualified law.
    pub fn prepare_wave_absorption(
        &self,
        omega_rad_s: f64,
    ) -> Result<BulkWaveAbsorption, MaterialError> {
        if !omega_rad_s.is_finite() || omega_rad_s <= 0.0 {
            return Err(invalid(
                "bulk wave absorption needs positive finite frequency",
            ));
        }
        let absorption = self.absorption();
        let state = absorption.state();
        let wave = absorption.thermoviscous_wave(omega_rad_s)?;
        // Let q=(omega/(k*c))² and y=omega*kappa/(rho*Cp*c²).
        // T'/u = alpha_p*T0/(rho*Cp) * Z_t/(q+i*gamma*y).
        // The native root retains the denominator BEFORE cancellation; using
        // rounded k in the heat equation loses it at coincident diffusion scales.
        let temperature_over_velocity = if self.isobaric_expansion_per_k() == 0.0 {
            C64::ZERO
        } else {
            (wave.traction_impedance_pa_s_m() / wave.thermal_denominator).scale(
                (self.isobaric_expansion_per_k() * self.temperature_k()
                    / state.transport().specific_heat_j_kg_k())
                    / state.acoustic().density_kg_m3(),
            )
        };
        let longitudinal_viscosity = (4.0 / 3.0) * state.transport().dynamic_viscosity_pa_s()
            + absorption.bulk_viscosity_pa_s();
        if !finite(temperature_over_velocity) || !longitudinal_viscosity.is_finite() {
            return Err(invalid(
                "bulk thermal wave coefficients are unrepresentable",
            ));
        }
        Ok(BulkWaveAbsorption {
            equilibrium: *self,
            wave,
            omega_rad_s,
            temperature_over_velocity,
            longitudinal_viscosity,
        })
    }
}

impl BulkWaveAbsorption {
    /// Original source-coherent acoustic, caloric and thermal snapshot.
    pub const fn equilibrium(&self) -> &BulkThermoviscousEquilibrium {
        &self.equilibrium
    }
    /// Original k and distinct pressure/traction impedances.
    pub const fn wave(&self) -> BulkThermoviscousWave {
        self.wave
    }
    pub const fn angular_frequency_rad_s(&self) -> f64 {
        self.omega_rad_s
    }

    /// Observe one forward wave driven by the stated peak velocity at x=0.
    /// Distance is the physical length along propagation.
    pub fn at_distance(
        &self,
        origin_velocity_m_s: C64,
        distance_m: f64,
    ) -> Result<BulkWaveField, MaterialError> {
        admit(origin_velocity_m_s, distance_m)?;
        let velocity = propagate(
            origin_velocity_m_s,
            self.wave.wavenumber_rad_m(),
            distance_m,
        )?;
        self.field(velocity)
    }

    fn field(&self, velocity: C64) -> Result<BulkWaveField, MaterialError> {
        let k = self.wave.wavenumber_rad_m();
        let ik = C64::new(-k.im, k.re);
        let t0 = self.equilibrium.temperature_k();
        let kappa = self
            .equilibrium
            .absorption()
            .state()
            .transport()
            .thermal_conductivity_w_m_k();
        let temperature = self.temperature_over_velocity * velocity;
        let gradient_t = ik * temperature;
        let pressure = self.wave.pressure_impedance_pa_s_m() * velocity;
        let traction = self.wave.traction_impedance_pa_s_m() * velocity;
        let heat_flux = -gradient_t.scale(kappa);
        let gradient_u = ik * velocity;
        // Weight before squaring; small physical phasors need not have a
        // representable unweighted square. Every loss is an original positive
        // constitutive form, not a pressure-work or endpoint-power subtraction.
        let viscous = gradient_u
            .scale(det::sqrt(0.5 * self.longitudinal_viscosity))
            .norm_sq();
        let conductive = gradient_t.scale(det::sqrt(0.5 * kappa / t0)).norm_sq();
        let entropy = conductive / t0;
        let mechanical = velocity
            .scale(det::sqrt(0.5 * self.wave.traction_impedance_pa_s_m().re))
            .norm_sq();
        let thermal = temperature
            .scale(det::sqrt((0.5 * kappa / t0) * k.im))
            .norm_sq();
        if ![
            velocity,
            gradient_u,
            pressure,
            traction,
            temperature,
            heat_flux,
        ]
        .into_iter()
        .all(finite)
            || ![
                viscous,
                conductive,
                entropy,
                mechanical,
                thermal,
                mechanical + thermal,
            ]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(invalid(
                "bulk wave field or physical work is unrepresentable",
            ));
        }
        Ok(BulkWaveField {
            velocity_m_s: velocity,
            velocity_gradient_per_s: gradient_u,
            pressure_pa: pressure,
            normal_traction_pa: traction,
            temperature_perturbation_k: temperature,
            conductive_heat_flux_w_m2: heat_flux,
            mechanical_power_flux_w_m2: mechanical,
            thermal_availability_flux_w_m2: thermal,
            viscous_heating_w_m3: viscous,
            conductive_entropy_w_m3_k: entropy,
            thermal_availability_loss_w_m3: conductive,
        })
    }

    /// Integrate the same physical losses over a constant-area cell [start,end].
    /// An expm1 integral retains weak absorption below endpoint-power roundoff,
    /// and only decaying exponentials are formed in opaque cells. Common phase
    /// cancels in every work term, so this query performs no phase evaluation,
    /// spatial quadrature, material query or allocation. Equal endpoints give
    /// zero loss. Positive finite physical area is explicit and applied once.
    pub fn absorption_between(
        &self,
        origin_velocity_m_s: C64,
        start_m: f64,
        end_m: f64,
        area_m2: f64,
    ) -> Result<BulkWaveCellAbsorption, MaterialError> {
        admit(origin_velocity_m_s, start_m)?;
        if !end_m.is_finite() || end_m < start_m || !area_m2.is_finite() || area_m2 <= 0.0 {
            return Err(invalid(
                "bulk absorption needs ordered finite depths and positive finite area",
            ));
        }
        let alpha = self.wave.wavenumber_rad_m().im;
        let width = end_m - start_m;
        let depth = alpha * width;
        let left_velocity = origin_velocity_m_s.scale(det::exp(-alpha * start_m));
        let left = self.field(left_velocity)?;
        let right = self.field(left_velocity.scale(det::exp(-depth)))?;
        let length = if depth == 0.0 {
            width
        } else if depth < 0.5 {
            width * (-det::expm1(-2.0 * depth) / (2.0 * depth))
        } else {
            (-0.5 * det::expm1(-2.0 * depth)) / alpha
        };
        let viscous = (left.viscous_heating_w_m3 * length) * area_m2;
        let conductive = (left.thermal_availability_loss_w_m3 * length) * area_m2;
        let entropy = conductive / self.equilibrium.temperature_k();
        let absorbed = viscous + conductive;
        let incoming = left.acoustic_power_flux_w_m2() * area_m2;
        let outgoing = right.acoustic_power_flux_w_m2() * area_m2;
        let residual =
            (Dd::from_f64(incoming) - Dd::from_f64(outgoing) - Dd::from_f64(absorbed)).to_f64();
        if ![
            viscous, conductive, entropy, absorbed, incoming, outgoing, residual,
        ]
        .into_iter()
        .all(f64::is_finite)
        {
            return Err(invalid("bulk wave cell absorption is unrepresentable"));
        }
        Ok(BulkWaveCellAbsorption {
            incoming_acoustic_power_w: incoming,
            outgoing_acoustic_power_w: outgoing,
            viscous_heat_w: viscous,
            conductive_entropy_w_k: entropy,
            thermal_availability_loss_w: conductive,
            acoustic_power_loss_w: absorbed,
            balance_residual_w: residual,
        })
    }
}

fn propagate(origin: C64, k: C64, distance: f64) -> Result<C64, MaterialError> {
    let amplitude = det::exp(-k.im * distance);
    if amplitude == 0.0 || origin == C64::ZERO {
        return Ok(C64::ZERO);
    }
    let angle = k.re * distance;
    if !angle.is_finite() {
        return Err(invalid("surviving bulk wave phase is unrepresentable"));
    }
    let (sin, cos) = det::sin_cos(angle);
    Ok(origin * C64::new(cos, sin).scale(amplitude))
}

fn finite(value: C64) -> bool {
    value.re.is_finite() && value.im.is_finite()
}
fn admit(velocity: C64, distance: f64) -> Result<(), MaterialError> {
    if finite(velocity) && distance.is_finite() && distance >= 0.0 {
        Ok(())
    } else {
        Err(invalid(
            "bulk wave needs finite velocity and finite nonnegative distance",
        ))
    }
}

#[cfg(test)]
mod tests;
