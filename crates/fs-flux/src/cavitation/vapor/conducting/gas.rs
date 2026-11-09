//! Bounded homogeneous gas-vapor contents for the resolved-liquid bubble.
//!
//! The immutable original dry inventory and saturated water share ONE volume
//! and ONE temperature. Partial pressures add; the original dry Cv and IF97
//! phase calorics add in the first law. This additive real-vapor/ideal-dry
//! approximation neglects cross-virial and interspecies fugacity corrections.
//! It is a dilute, common-temperature, rapidly mixed closure, not a general
//! real-mixture EOS or an internal gas PDE.
//! Noncondensable dissolution, chemistry, gas-side concentration layers,
//! internal thermal boundary layers and interspecies fugacity corrections
//! require separate physical owners when the declared approximation fails.
//!
//! The caller supplies lower bounds for the WET mixture's conductivity and
//! water/noncondensable diffusivity over an explicit T/p/composition domain.
//! No transport number is inferred from dry conductivity, species names or an
//! unrelated kinetic fit. Bounds can be source-derived or explicitly authored
//! problem data; their physical authority remains the caller's.
//!
//! Zero gas flux gives u_g(R)=Rdot and u_v(R)=Rdot-j/rho_v. The species momentum
//! flux therefore retains j^2/rho_v in the existing recoil law. Fast diffusion
//! is screened using dY_v/dr|R=(1-Y_v)j/(rho_mix D). The reported gradient and
//! relaxation scales are applicability estimates, not rigorous error bounds.
//!
//! The internal Fourier power is Q_inside=Q_liquid-L*mdot_v. The estimate
//! |Q_inside|/(4 pi k_mix R) screens the temperature variation needed to carry
//! it; it is distinct from the dry component's heat m_g Cv Tdot+p_g Vdot.
//!
//! Physical context: Hao, Zhang & Prosperetti, Phys. Rev. Fluids 2, 034303
//! (2017), equations 2--6 and 17--19, explicitly resolve internal diffusion:
//! <https://link.aps.org/accepted/10.1103/PhysRevFluids.2.034303>.

use super::super::{VaporBubbleObservation, VaporBubbleThermalState, bound};
use super::{ConductingVaporBubbleError, GasVaporObservation};
use crate::cavitation::positive;
use core::f64::consts::PI;
use fs_material::cavitation::vapor::gas::{BubbleGasSample, ConservedBubbleGas};

/// Explicit wet-mixture transport lower bounds and their complete local domain.
/// Source identity/uncertainty is retained by the caller's material selection;
/// this numerical view does not claim that authored values are measured data.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasVaporTransportBounds {
    pub minimum_temperature_k: f64,
    pub maximum_temperature_k: f64,
    pub minimum_total_pressure_pa: f64,
    pub maximum_total_pressure_pa: f64,
    pub minimum_vapor_mass_fraction: f64,
    pub maximum_vapor_mass_fraction: f64,
    pub minimum_thermal_conductivity_w_m_k: f64,
    /// Lower bound for relevant composition relaxation. For a dry blend it
    /// must also bound the slowest active mixing mode, not one selected pair.
    pub minimum_binary_diffusivity_m2_s: f64,
}

/// Actual-mixture inertia, equilibration and spatial-uniformity limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasVaporLimits {
    pub transport: GasVaporTransportBounds,
    /// Largest (rho_v+rho_g)/rho_l, in (0,0.1].
    pub maximum_total_density_ratio: f64,
    /// Largest mixture-relative/wall speed divided by frozen mixture sound.
    pub maximum_mixture_mach_number: f64,
    /// Largest (R/c_mix)*|p_b'/p_b|, in (0,0.1].
    pub maximum_pressure_relaxation_ratio: f64,
    /// Largest estimated center/interface temperature variation [K].
    pub maximum_internal_temperature_difference_k: f64,
    /// Relative composition-gradient scale, in (0,0.1].
    pub maximum_relative_composition_gradient: f64,
    /// Thermal relaxation time times admitted change/forcing rate, in (0,0.1].
    pub maximum_thermal_response_ratio: f64,
    /// Diffusion relaxation time times admitted change/forcing rate, in (0,0.1].
    pub maximum_diffusion_response_ratio: f64,
    /// Authored upper rate [1/s] covering imposed bandwidth and relative
    /// pressure/temperature change rates, including extrema where a local
    /// first derivative vanishes. Zero asserts genuinely constant forcing.
    /// The temporal step ceiling must independently resolve that forcing.
    pub maximum_imposed_rate_per_s: f64,
}

fn invalid(message: &'static str) -> ConductingVaporBubbleError {
    ConductingVaporBubbleError::InvalidInput(message)
}

fn in_domain(
    quantity: &'static str,
    value: f64,
    minimum: f64,
    maximum: f64,
) -> Result<(), ConductingVaporBubbleError> {
    if !value.is_finite() || value < minimum || value > maximum {
        return Err(ConductingVaporBubbleError::OutsideDomain {
            quantity,
            value,
            minimum,
            maximum,
        });
    }
    Ok(())
}

impl GasVaporTransportBounds {
    fn validate(self) -> Result<(), ConductingVaporBubbleError> {
        if !positive(self.minimum_temperature_k)
            || !positive(self.maximum_temperature_k)
            || self.minimum_temperature_k >= self.maximum_temperature_k
            || !positive(self.minimum_total_pressure_pa)
            || !positive(self.maximum_total_pressure_pa)
            || self.minimum_total_pressure_pa >= self.maximum_total_pressure_pa
            || !self.minimum_vapor_mass_fraction.is_finite()
            || !self.maximum_vapor_mass_fraction.is_finite()
            || self.minimum_vapor_mass_fraction < 0.0
            || self.maximum_vapor_mass_fraction > 1.0
            || self.minimum_vapor_mass_fraction >= self.maximum_vapor_mass_fraction
            || !positive(self.minimum_thermal_conductivity_w_m_k)
            || !positive(self.minimum_binary_diffusivity_m2_s)
        {
            return Err(invalid(
                "wet-mixture transport bounds and explicit source domain",
            ));
        }
        Ok(())
    }

    fn admit(
        self,
        temperature: f64,
        pressure: f64,
        vapor_fraction: f64,
    ) -> Result<(), ConductingVaporBubbleError> {
        in_domain(
            "wet-mixture transport temperature",
            temperature,
            self.minimum_temperature_k,
            self.maximum_temperature_k,
        )?;
        in_domain(
            "wet-mixture transport total pressure",
            pressure,
            self.minimum_total_pressure_pa,
            self.maximum_total_pressure_pa,
        )?;
        in_domain(
            "wet-mixture transport vapor mass fraction",
            vapor_fraction,
            self.minimum_vapor_mass_fraction,
            self.maximum_vapor_mass_fraction,
        )
    }
}

impl GasVaporLimits {
    pub fn validate(self) -> Result<(), ConductingVaporBubbleError> {
        self.transport.validate()?;
        for value in [
            self.maximum_total_density_ratio,
            self.maximum_mixture_mach_number,
            self.maximum_pressure_relaxation_ratio,
            self.maximum_relative_composition_gradient,
            self.maximum_thermal_response_ratio,
            self.maximum_diffusion_response_ratio,
        ] {
            if !positive(value) || value > 0.1 {
                return Err(invalid(
                    "gas-vapor dimensionless limits must lie in (0,0.1]",
                ));
            }
        }
        if !positive(self.maximum_internal_temperature_difference_k)
            || self.maximum_internal_temperature_difference_k
                > 0.05 * self.transport.minimum_temperature_k
            || !self.maximum_imposed_rate_per_s.is_finite()
            || self.maximum_imposed_rate_per_s < 0.0
        {
            return Err(invalid(
                "bounded internal temperature variation and imposed rate",
            ));
        }
        Ok(())
    }
}

#[derive(Clone, Debug, PartialEq)]
pub(super) struct GasBinding {
    pub(super) inventory: ConservedBubbleGas,
    limits: GasVaporLimits,
}

impl GasBinding {
    pub(super) fn new(
        inventory: ConservedBubbleGas,
        limits: GasVaporLimits,
        thermal: &VaporBubbleThermalState,
    ) -> Result<Self, ConductingVaporBubbleError> {
        limits.validate()?;
        let reference = inventory.reference_parcel();
        if reference.state().state().temperature != thermal.state.temperature_k()
            || (reference.volume_m3() - thermal.volume_m3).abs()
                > 64.0 * f64::EPSILON * reference.volume_m3().max(thermal.volume_m3)
        {
            return Err(invalid(
                "gas inventory must match accepted bubble temperature and volume",
            ));
        }
        let sample = inventory.sample(thermal.state.temperature_k(), thermal.volume_m3)?;
        let density = thermal.vapor.vapor_density_kg_m3() + sample.density_kg_m3;
        if !positive(density) {
            return Err(invalid("represented initial mixture density"));
        }
        limits.transport.admit(
            thermal.state.temperature_k(),
            thermal.vapor.pressure_pa() + sample.pressure_pa,
            thermal.vapor.vapor_density_kg_m3() / density,
        )?;
        // Direct original EOS/caloric/transport qualification. This creates no
        // second accepted parcel and does not change the reference inventory.
        inventory.qualify(thermal.state.temperature_k(), thermal.volume_m3)?;
        Ok(Self { inventory, limits })
    }

    pub(super) fn sample(
        &self,
        temperature_k: f64,
        volume_m3: f64,
    ) -> Result<BubbleGasSample, ConductingVaporBubbleError> {
        let controls = self.inventory.controls();
        // A narrow gas bracket may be crossed only by an RK trial. Retain a
        // typed refinable domain error before the original source samples it.
        in_domain(
            "noncondensable caloric temperature",
            temperature_k,
            controls.minimum_temperature_k,
            controls.maximum_temperature_k,
        )?;
        Ok(self.inventory.sample(temperature_k, volume_m3)?)
    }

    pub(super) fn observe(
        &self,
        sample: BubbleGasSample,
        bubble: &VaporBubbleObservation,
        liquid_density_kg_m3: f64,
    ) -> Result<GasVaporObservation, ConductingVaporBubbleError> {
        let radius = bubble.state.radius_m();
        let temperature = bubble.state.temperature_k();
        let vapor_density = bubble.vapor.vapor_density_kg_m3();
        let density = vapor_density + sample.density_kg_m3;
        let vapor_fraction = vapor_density / density;
        if !positive(density) || !positive(vapor_fraction) || vapor_fraction >= 1.0 {
            return Err(invalid("represented positive gas and vapor mass fractions"));
        }
        self.limits
            .transport
            .admit(temperature, bubble.bubble_pressure_pa, vapor_fraction)?;
        bound(
            "gas-vapor/liquid density ratio",
            density / liquid_density_kg_m3,
            self.limits.maximum_total_density_ratio,
        )?;

        // Frozen-composition acoustic derivative of the declared additive
        // free-energy approximation, using ORIGINAL component derivatives:
        // c^2 = p_rho|T + T p_T|rho^2/(rho C_v,volume).
        // The IF97 endpoint remains the vapor partial-density state. No
        // equilibrium-phase gamma or dry-component sound speed is substituted.
        let water = bubble.vapor.saturation().vapor();
        let compressibility = water.isothermal_compressibility().value();
        let expansion = water.thermal_expansion().value();
        let vapor_cv = water.cv().value();
        let vapor_cp = water.cp().value();
        let cv_volume = vapor_density * vapor_cv + sample.density_kg_m3 * sample.cv_j_kg_k;
        // The sum of component rho*Cp is a conservative screening capacity
        // for this additive model (Cauchy-Schwarz bounds the mixture Cp).
        // It is not presented as an independently measured wet-mixture Cp.
        let cp_volume = vapor_density * vapor_cp + sample.density_kg_m3 * sample.cp_j_kg_k;
        if !positive(compressibility)
            || !expansion.is_finite()
            || !positive(cv_volume)
            || !positive(cp_volume)
        {
            return Err(invalid(
                "positive original mixture acoustic/thermal derivatives",
            ));
        }
        let pressure_temperature = expansion / compressibility + sample.pressure_pa / temperature;
        let sound_squared = (1.0 / compressibility
            + sample.pressure_pa
            + temperature * pressure_temperature * pressure_temperature / cv_volume)
            / density;
        if !positive(sound_squared) {
            return Err(invalid("represented frozen gas-vapor acoustic derivative"));
        }
        let sound = fs_math::det::sqrt(sound_squared);
        let mixture_relative_velocity = bubble.vapor_mass_flux_kg_m2_s / density;
        let mixture_velocity = bubble.wall_velocity_m_s - mixture_relative_velocity;
        let mach = bubble
            .wall_velocity_m_s
            .abs()
            .max(mixture_relative_velocity.abs())
            .max(mixture_velocity.abs())
            / sound;
        bound(
            "gas-vapor mixture Mach",
            mach,
            self.limits.maximum_mixture_mach_number,
        )?;

        let volume_rate = bubble.volume_rate_m3_s / bubble.volume_m3;
        let temperature_rate = bubble.temperature_rate_k_s / temperature;
        let pressure_rate = bubble.vapor.pressure_temperature_derivative_pa_k()
            * bubble.temperature_rate_k_s
            + sample.pressure_pa * (temperature_rate - volume_rate);
        let relative_pressure_rate = (pressure_rate / bubble.bubble_pressure_pa).abs();
        let pressure_relaxation = (radius / sound) * relative_pressure_rate;
        bound(
            "gas-vapor pressure relaxation",
            pressure_relaxation,
            self.limits.maximum_pressure_relaxation_ratio,
        )?;

        let transport = self.limits.transport;
        let internal_heat =
            bubble.heat_power_w - bubble.vapor.latent_enthalpy_j_kg() * bubble.vapor_mass_rate_kg_s;
        let temperature_difference = internal_heat.abs()
            / (4.0 * PI * transport.minimum_thermal_conductivity_w_m_k * radius);
        bound(
            "estimated internal gas-vapor temperature difference",
            temperature_difference,
            self.limits.maximum_internal_temperature_difference_k,
        )?;

        // max(1,Y_g/Y_v) screens relative variations in BOTH species,
        // including vapor depletion next to a condensing gas-rich interface.
        let composition_factor = (sample.density_kg_m3 / vapor_density).max(1.0);
        let flux = bubble.vapor_mass_flux_kg_m2_s.abs();
        let concentration_gradient = (flux * radius
            / (density * transport.minimum_binary_diffusivity_m2_s))
            * composition_factor;
        bound(
            "gas-vapor relative composition gradient",
            concentration_gradient,
            self.limits.maximum_relative_composition_gradient,
        )?;
        let composition_rate = (3.0 * flux / (density * radius)) * composition_factor;
        let rate = self
            .limits
            .maximum_imposed_rate_per_s
            .max(temperature_rate.abs())
            .max(volume_rate.abs())
            .max(relative_pressure_rate)
            .max(composition_rate);
        let thermal_time =
            radius * radius * cp_volume / transport.minimum_thermal_conductivity_w_m_k;
        let diffusion_time = radius * radius / transport.minimum_binary_diffusivity_m2_s;
        if !positive(thermal_time) || !positive(diffusion_time) || !rate.is_finite() {
            return Err(invalid("represented wet-mixture equilibration times"));
        }
        let thermal_response = thermal_time * rate;
        let diffusion_response = diffusion_time * rate;
        bound(
            "gas-vapor thermal response",
            thermal_response,
            self.limits.maximum_thermal_response_ratio,
        )?;
        bound(
            "gas-vapor species diffusion response",
            diffusion_response,
            self.limits.maximum_diffusion_response_ratio,
        )?;

        let energy_rate = sample.heat_capacity_j_k * bubble.temperature_rate_k_s;
        let work = sample.pressure_pa * bubble.volume_rate_m3_s;
        let heat = energy_rate + work;
        if ![internal_heat, energy_rate, work, heat]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(invalid("represented original gas energy/heat/work rates"));
        }
        Ok(GasVaporObservation {
            mass_kg: sample.mass_kg,
            density_kg_m3: sample.density_kg_m3,
            partial_pressure_pa: sample.pressure_pa,
            cv_j_kg_k: sample.cv_j_kg_k,
            energy_rate_w: energy_rate,
            pressure_work_power_w: work,
            heat_power_w: heat,
            mixture_density_kg_m3: density,
            mixture_sound_speed_m_s: sound,
            vapor_mass_fraction: vapor_fraction,
            internal_heat_power_w: internal_heat,
            estimated_internal_temperature_difference_k: temperature_difference,
            relative_composition_gradient_scale: concentration_gradient,
            thermal_response_ratio: thermal_response,
            diffusion_response_ratio: diffusion_response,
            pressure_relaxation_ratio: pressure_relaxation,
        })
    }
}

#[cfg(test)]
#[path = "gas/tests.rs"]
mod tests;
