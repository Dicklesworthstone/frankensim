//! Nonlinear pure-water vapor bubbles with phase mass and latent heat.
//!
//! The original IF97 saturation source supplies p(T), rho_v(T), u_v(T), h_l(T)
//! and their analytic tangents. Vapor occupies the geometric bubble volume:
//! m_v=rho_v V, U_v=m_v u_v. Condensate returns to the surrounding liquid;
//! liquid-droplet volume is not inserted into a nominal vapor cavity.
//!
//! Positive j evaporates into vapor. Mass conservation requires
//! Rdot=U_l+j/rho_l, where U_l is the liquid velocity at the interface. Its
//! exact incompressible momentum equation retains the Stefan slip and recoil:
//! R U_l'+2 Rdot U_l-U_l^2/2 =
//! (p_v-p_inf-2 sigma/R-4 mu U_l/R+j^2(1/rho_v-1/rho_l))/rho_l.
//! This is Eq. 2.1 of Sullivan et al., JFM 948 (2022), with their J=-j:
//! <https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/inertiothermal-vapour-bubble-growth/20819A1EC619FA2D8DA80FBDDD91103C>.
//!
//! The leading thermal interface law is U_v'=Q_l+h_l m_v'-p_v V'. Together
//! with saturation mass it gives a two-by-two scalar solve for Tdot and j;
//! no fitted evaporation time, accommodation coefficient or heat coefficient
//! is introduced. Q_l is an explicit heat-power port on the LIQUID side of
//! the interface. Neglected kinetic, viscous and liquid-pressure enthalpy
//! contributions, and the flat-interface equilibrium approximation, receive
//! independent stage admission. Those scales are applicability screens, not
//! rigorous error bounds on a nonequilibrium interface.
//!
//! The bulk liquid, explicit surface tension and external heat reservoir remain
//! frozen. Actual original-water density, viscosity and sound-speed variation
//! is checked at every interface temperature. Vapor pressure is homogeneous;
//! vapor flow/phase-flux Mach and pressure-relaxation time scales are bounded.
//! No noncondensable gas, kinetic phase resistance, finite liquid inventory,
//! resolved liquid heat equation, Kelvin equilibrium, liquid compressibility,
//! shock, jet or erosion law is claimed. In particular this is not a
//! phase-changing Keller--Miksis equation.

use super::{
    BubbleLimits, BubbleNumerics, CavitationError, finite, positive, radial_liquid_acceleration,
    stepper,
};
use core::f64::consts::PI;
use fs_material::MaterialError;
use fs_material::cavitation::{SaturatedWaterCavitation, vapor::SaturatedWaterVapor};

/// Geometric radius and LIQUID velocity, which differs from Rdot under phase change.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VaporBubbleState {
    time_s: f64,
    radius_m: f64,
    liquid_velocity_m_s: f64,
    temperature_k: f64,
}

impl VaporBubbleState {
    pub fn try_new(
        time_s: f64,
        radius_m: f64,
        liquid_velocity_m_s: f64,
        temperature_k: f64,
    ) -> Result<Self, VaporBubbleError> {
        if !time_s.is_finite()
            || !positive(radius_m)
            || !liquid_velocity_m_s.is_finite()
            || !positive(temperature_k)
        {
            return Err(VaporBubbleError::InvalidInput("vapor bubble state"));
        }
        Ok(Self {
            time_s,
            radius_m,
            liquid_velocity_m_s,
            temperature_k,
        })
    }
    pub const fn time_s(self) -> f64 {
        self.time_s
    }
    pub const fn radius_m(self) -> f64 {
        self.radius_m
    }
    pub const fn liquid_velocity_m_s(self) -> f64 {
        self.liquid_velocity_m_s
    }
    pub const fn temperature_k(self) -> f64 {
        self.temperature_k
    }
}

/// Declared source, spherical-flow and thermal-interface applicability domain.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VaporBubbleLimits {
    /// Radius, liquid/interface Mach, and actual vapor/liquid density limits.
    pub bubble: BubbleLimits,
    pub minimum_temperature_k: f64,
    pub maximum_temperature_k: f64,
    /// Largest absolute pressure departure from the ORIGINAL liquid pressure,
    /// divided by rho_l c_l^2; in (0,0.1].
    pub maximum_pressure_ratio: f64,
    /// Bounds vapor velocity, phase-relative/interface Mach, and R/c_v times
    /// |p_v'/p_v|, independently of the liquid Mach check; in (0,0.1].
    pub maximum_vapor_mach_number: f64,
    /// Largest source rho_l, mu_l, c_l or T fractional change; in (0,0.05].
    pub maximum_frozen_liquid_relative_change: f64,
    /// Largest omitted interface specific-energy scale / latent enthalpy,
    /// and liquid-pressure chemical-potential scale / p_v v_v; in (0,0.01].
    pub maximum_neglected_enthalpy_fraction: f64,
}

impl VaporBubbleLimits {
    pub fn validate(self) -> Result<(), VaporBubbleError> {
        self.bubble.validate()?;
        if !self.minimum_temperature_k.is_finite()
            || !self.maximum_temperature_k.is_finite()
            || self.minimum_temperature_k < 273.15
            || self.maximum_temperature_k > 623.15
            || self.minimum_temperature_k >= self.maximum_temperature_k
        {
            return Err(VaporBubbleError::InvalidInput(
                "IF97 vapor temperature interval",
            ));
        }
        for (name, value, ceiling) in [
            (
                "liquid pressure-ratio limit",
                self.maximum_pressure_ratio,
                0.1,
            ),
            ("vapor Mach limit", self.maximum_vapor_mach_number, 0.1),
            (
                "frozen-liquid variation limit",
                self.maximum_frozen_liquid_relative_change,
                0.05,
            ),
            (
                "neglected interface enthalpy limit",
                self.maximum_neglected_enthalpy_fraction,
                0.01,
            ),
        ] {
            if !positive(value) || value > ceiling {
                return Err(VaporBubbleError::InvalidInput(name));
            }
        }
        Ok(())
    }
}

/// Original adaptive radial controls plus physical thermal and balance scales.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VaporBubbleNumerics {
    /// Velocity tolerance applies to U_l; the geometric wall velocity is derived.
    pub radial: BubbleNumerics,
    pub temperature_absolute_tolerance_k: f64,
    pub energy_absolute_tolerance_j: f64,
    pub mass_absolute_tolerance_kg: f64,
    pub energy_balance_tolerance_j: f64,
    pub mass_balance_tolerance_kg: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum VaporBubbleError {
    Material(MaterialError),
    Dynamics(CavitationError),
    InvalidInput(&'static str),
    OutsideModel {
        quantity: &'static str,
        value: f64,
        maximum: f64,
    },
    TemperatureOutsideDomain {
        temperature_k: f64,
        minimum_k: f64,
        maximum_k: f64,
    },
    EnergyBalance {
        residual_j: f64,
        tolerance_j: f64,
    },
    MassBalance {
        residual_kg: f64,
        tolerance_kg: f64,
    },
    Cancelled,
}
impl core::fmt::Display for VaporBubbleError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "saturated-vapor bubble refused: {self:?}")
    }
}
impl std::error::Error for VaporBubbleError {}
impl From<MaterialError> for VaporBubbleError {
    fn from(error: MaterialError) -> Self {
        Self::Material(error)
    }
}
impl From<CavitationError> for VaporBubbleError {
    fn from(error: CavitationError) -> Self {
        if error == CavitationError::Cancelled {
            Self::Cancelled
        } else {
            Self::Dynamics(error)
        }
    }
}

fn bound(quantity: &'static str, value: f64, maximum: f64) -> Result<(), VaporBubbleError> {
    if !value.is_finite() || value < 0.0 || value > maximum {
        return Err(VaporBubbleError::OutsideModel {
            quantity,
            value,
            maximum,
        });
    }
    Ok(())
}

fn volume(radius: f64) -> Result<f64, VaporBubbleError> {
    let result = (4.0 * PI / 3.0) * radius * radius * radius;
    if !positive(result) {
        return Err(VaporBubbleError::InvalidInput("representable vapor volume"));
    }
    Ok(result)
}

/// Source snapshot supplied to the caller's deterministic liquid-side heat port.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VaporBubbleThermalState {
    pub state: VaporBubbleState,
    pub vapor: SaturatedWaterVapor,
    pub volume_m3: f64,
    pub vapor_mass_kg: f64,
}

/// Instantaneous phase flux, mechanical flow and independently reported limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VaporBubbleObservation {
    pub state: VaporBubbleState,
    /// Original phase endpoints and raw IF97 Gibbs discrepancy remain available.
    pub vapor: SaturatedWaterVapor,
    pub wall_velocity_m_s: f64,
    pub liquid_acceleration_m_s2: f64,
    pub temperature_rate_k_s: f64,
    pub vapor_mass_kg: f64,
    /// Positive for evaporation, negative for condensation.
    pub vapor_mass_rate_kg_s: f64,
    pub vapor_mass_flux_kg_m2_s: f64,
    pub volume_m3: f64,
    /// Geometric dV_b/dt, distinct from the displaced-liquid acoustic source.
    pub volume_rate_m3_s: f64,
    /// Actual incompressible liquid flux, 4 pi R^2 U_l [m3/s].
    pub liquid_volume_flow_m3_s: f64,
    /// Analytic derivative of actual liquid flux; no derivative of heat needed.
    pub liquid_volume_flow_rate_m3_s2: f64,
    pub liquid_interface_pressure_pa: f64,
    pub heat_power_w: f64,
    pub vapor_recoil_pressure_pa: f64,
    pub neglected_enthalpy_fraction: f64,
    pub neglected_specific_enthalpy_j_kg: f64,
    pub flat_interface_pressure_fraction: f64,
    pub frozen_liquid_relative_change: f64,
}

/// Whole-window source mass and energy accounts, committed atomically.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VaporBubbleStep {
    pub observation: VaporBubbleObservation,
    pub accepted_steps: usize,
    pub rejected_steps: usize,
    pub rhs_evaluations: usize,
    pub maximum_normalized_local_error: f64,
    pub vapor_mass_change_kg: f64,
    pub transferred_vapor_mass_kg: f64,
    pub mass_residual_kg: f64,
    /// Actual Delta(m_v u_v) in the released IF97 energy datum.
    pub vapor_energy_change_j: f64,
    pub supplied_heat_j: f64,
    /// Integral h_l dm_v in the SAME original source datum.
    pub transferred_enthalpy_j: f64,
    /// Integral p_v dV_b, positive for work done by vapor.
    pub pressure_work_j: f64,
    /// Original initial liquid enthalpy used only as a constant energy anchor.
    pub reference_liquid_enthalpy_j_kg: f64,
    pub reference_invariant_energy_change_j: f64,
    pub reference_invariant_transferred_enthalpy_j: f64,
    /// Delta(U-h_ref m) - Q - integral(h_l-h_ref)dm + integral p_v dV.
    /// This removes arbitrary common mass-energy datums from balance admission.
    pub energy_residual_j: f64,
    /// Same balance in the released absolute datum, including h_ref*mass residual.
    pub original_datum_energy_residual_j: f64,
    /// Integral |dm/dt| times the admitted omitted specific-enthalpy scale.
    /// This is a model-approximation indicator, not an error certificate.
    pub neglected_interface_energy_scale_j: f64,
}

/// Native phase-changing bubble in its original frozen IF97 surrounding liquid.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SaturatedVaporBubble {
    source: SaturatedWaterCavitation,
    state: VaporBubbleState,
    limits: VaporBubbleLimits,
}

/// Borrowed uncommitted candidate, allowing acoustic admission before publication.
pub struct VaporBubbleProposal<'a> {
    target: &'a mut SaturatedVaporBubble,
    step: VaporBubbleStep,
}
impl VaporBubbleProposal<'_> {
    pub const fn step(&self) -> &VaporBubbleStep {
        &self.step
    }
    pub fn commit(self) -> VaporBubbleStep {
        self.target.state = self.step.observation.state;
        self.step
    }
}

impl SaturatedVaporBubble {
    /// Bind vapor calorics and surrounding mechanics to the SAME original water.
    /// Surface tension stays explicitly authored on the retained water source.
    pub fn new(
        liquid: SaturatedWaterCavitation,
        state: VaporBubbleState,
        limits: VaporBubbleLimits,
    ) -> Result<Self, VaporBubbleError> {
        limits.validate()?;
        let result = Self {
            source: liquid,
            state,
            limits,
        };
        result.thermal_state(state)?;
        Ok(result)
    }
    pub const fn source(&self) -> &SaturatedWaterCavitation {
        &self.source
    }
    pub const fn state(&self) -> VaporBubbleState {
        self.state
    }
    pub const fn limits(&self) -> VaporBubbleLimits {
        self.limits
    }

    fn thermal_state(
        &self,
        state: VaporBubbleState,
    ) -> Result<(VaporBubbleThermalState, f64), VaporBubbleError> {
        let limits = self.limits;
        if state.radius_m < limits.bubble.minimum_radius_m
            || state.radius_m > limits.bubble.maximum_radius_m
        {
            return Err(CavitationError::RadiusOutsideDomain {
                radius_m: state.radius_m,
                minimum_radius_m: limits.bubble.minimum_radius_m,
                maximum_radius_m: limits.bubble.maximum_radius_m,
            }
            .into());
        }
        if state.temperature_k < limits.minimum_temperature_k
            || state.temperature_k > limits.maximum_temperature_k
        {
            return Err(VaporBubbleError::TemperatureOutsideDomain {
                temperature_k: state.temperature_k,
                minimum_k: limits.minimum_temperature_k,
                maximum_k: limits.maximum_temperature_k,
            });
        }
        let original = self.source.state();
        let (vapor, sampled_liquid) = if state.temperature_k == original.temperature_k() {
            (
                SaturatedWaterVapor::try_from_saturation(*self.source.saturation())?,
                original,
            )
        } else {
            let vapor = SaturatedWaterVapor::at_temperature_k(state.temperature_k)?;
            let sampled = SaturatedWaterCavitation::try_new(
                *vapor.saturation(),
                original.surface_tension_n_m(),
            )?
            .state();
            (vapor, sampled)
        };
        let variation = (sampled_liquid.liquid_density_kg_m3() / original.liquid_density_kg_m3()
            - 1.0)
            .abs()
            .max(
                (sampled_liquid.dynamic_viscosity_pa_s() / original.dynamic_viscosity_pa_s() - 1.0)
                    .abs(),
            )
            .max(
                (sampled_liquid.liquid_sound_speed_m_s() / original.liquid_sound_speed_m_s() - 1.0)
                    .abs(),
            )
            .max((state.temperature_k / original.temperature_k() - 1.0).abs());
        bound(
            "frozen liquid coefficient variation",
            variation,
            limits.maximum_frozen_liquid_relative_change,
        )?;
        bound(
            "vapor/liquid density ratio",
            vapor.vapor_density_kg_m3() / original.liquid_density_kg_m3(),
            limits.bubble.maximum_vapor_density_ratio,
        )?;
        bound(
            "liquid interface Mach",
            state.liquid_velocity_m_s.abs() / original.liquid_sound_speed_m_s(),
            limits.bubble.maximum_wall_mach_number,
        )?;
        let volume_m3 = volume(state.radius_m)?;
        let vapor_mass_kg = vapor.vapor_density_kg_m3() * volume_m3;
        if !positive(vapor_mass_kg) {
            return Err(VaporBubbleError::InvalidInput(
                "representable vapor inventory",
            ));
        }
        Ok((
            VaporBubbleThermalState {
                state,
                vapor,
                volume_m3,
                vapor_mass_kg,
            },
            variation,
        ))
    }

    pub fn observe(
        &self,
        absolute_pressure_pa: f64,
        heat_power_w: f64,
    ) -> Result<VaporBubbleObservation, VaporBubbleError> {
        let (thermal, variation) = self.thermal_state(self.state)?;
        self.observe_thermal(thermal, variation, absolute_pressure_pa, heat_power_w)
    }

    fn observe_thermal(
        &self,
        thermal: VaporBubbleThermalState,
        frozen_liquid_relative_change: f64,
        pressure: f64,
        heat: f64,
    ) -> Result<VaporBubbleObservation, VaporBubbleError> {
        if !pressure.is_finite() || !heat.is_finite() {
            return Err(VaporBubbleError::InvalidInput(
                "finite absolute pressure and interface heat power",
            ));
        }
        let state = thermal.state;
        let vapor = thermal.vapor;
        let liquid = self.source.state();
        let radius = state.radius_m;
        let velocity = state.liquid_velocity_m_s;
        let area = 4.0 * PI * radius * radius;
        let rho_l = liquid.liquid_density_kg_m3();
        let rho_v = vapor.vapor_density_kg_m3();
        let ratio = rho_v / rho_l;
        let latent = vapor.latent_enthalpy_j_kg();
        let density_tangent = (radius / 3.0) * vapor.vapor_density_temperature_derivative_kg_m3_k();
        let capacity = (radius / 3.0) * vapor.phase_temperature_capacity_j_m3_k();
        let volume_latent = rho_v * latent;
        let heat_flux = heat / area;
        let determinant = (1.0 - ratio) * capacity + (volume_latent / rho_l) * density_tangent;
        if !positive(area) || !positive(determinant) || !heat_flux.is_finite() {
            return Err(VaporBubbleError::InvalidInput(
                "representable saturated phase balance",
            ));
        }
        // Solve directly for Tdot and j. Computing j=rho_l*(Rdot-U_l)
        // would erase the small mass flux at ordinary water density ratios.
        let temperature_rate = ((1.0 - ratio) * heat_flux - volume_latent * velocity) / determinant;
        let mass_flux = (rho_v * capacity * velocity
            + density_tangent * (heat_flux - volume_latent * velocity))
            / determinant;
        let stefan_slip = mass_flux / rho_l;
        let wall_velocity = velocity + stefan_slip;
        let vapor_relative_velocity = mass_flux / rho_v;
        let vapor_fluid_velocity = -density_tangent * temperature_rate / rho_v;
        let vapor_mach = vapor_relative_velocity
            .abs()
            .max(vapor_fluid_velocity.abs())
            .max(wall_velocity.abs())
            / vapor.vapor_sound_speed_m_s();
        bound(
            "vapor or phase-relative Mach",
            vapor_mach,
            self.limits.maximum_vapor_mach_number,
        )?;
        bound(
            "vapor pressure relaxation ratio",
            (radius / vapor.vapor_sound_speed_m_s())
                * (vapor.pressure_temperature_derivative_pa_k() * temperature_rate
                    / vapor.pressure_pa())
                .abs(),
            self.limits.maximum_vapor_mach_number,
        )?;
        bound(
            "geometric interface Mach",
            wall_velocity.abs() / liquid.liquid_sound_speed_m_s(),
            self.limits.bubble.maximum_wall_mach_number,
        )?;
        let recoil = mass_flux * (vapor_relative_velocity - stefan_slip);
        let viscous_pressure = 4.0 * liquid.dynamic_viscosity_pa_s() * velocity / radius;
        let interface_departure =
            -2.0 * liquid.surface_tension_n_m() / radius - viscous_pressure + recoil;
        let interface_pressure = vapor.pressure_pa() + interface_departure;
        let bulk = rho_l * liquid.liquid_sound_speed_m_s() * liquid.liquid_sound_speed_m_s();
        bound(
            "original liquid pressure excursion",
            (pressure - liquid.saturation_pressure_pa())
                .abs()
                .max((interface_pressure - liquid.saturation_pressure_pa()).abs())
                / bulk,
            self.limits.maximum_pressure_ratio,
        )?;
        // The pressure term also includes the source liquid's thermal expansion:
        // |h_p| = |v-T v_T| <= v+T|v_T|. These remain first-order applicability
        // scales because the liquid interface may be tensile, outside IF97's EOS.
        let source_liquid = vapor.saturation().liquid();
        let pressure_enthalpy = interface_departure.abs()
            * (source_liquid.specific_volume().value()
                + state.temperature_k * source_liquid.dv_dt_at_pressure().value().abs());
        let omitted = pressure_enthalpy
            + viscous_pressure.abs() / rho_l
            + 0.5 * (vapor_relative_velocity * vapor_relative_velocity + stefan_slip * stefan_slip);
        let neglected_enthalpy_fraction = omitted / latent;
        bound(
            "neglected interface enthalpy / latent heat",
            neglected_enthalpy_fraction,
            self.limits.maximum_neglected_enthalpy_fraction,
        )?;
        let flat_interface_pressure_fraction = interface_departure.abs()
            * source_liquid.specific_volume().value()
            / (vapor.pressure_pa() / rho_v);
        bound(
            "flat-interface saturation pressure correction",
            flat_interface_pressure_fraction,
            self.limits.maximum_neglected_enthalpy_fraction,
        )?;
        // Preserve subtraction around the physical vapor/ambient pressure
        // before capillary stress, as in the original zero-mass-flux owner.
        let driving = -(pressure - vapor.pressure_pa()
            + 2.0 * liquid.surface_tension_n_m() / radius)
            - viscous_pressure
            + recoil;
        let acceleration =
            radial_liquid_acceleration(liquid, radius, velocity, stefan_slip, driving)?;
        let volume_rate = area * wall_velocity;
        let liquid_flow = area * velocity;
        let liquid_flow_rate = area * (acceleration + 2.0 * wall_velocity * velocity / radius);
        let mass_rate = area * mass_flux;
        for value in [
            temperature_rate,
            mass_flux,
            wall_velocity,
            recoil,
            interface_pressure,
            acceleration,
            volume_rate,
            liquid_flow,
            liquid_flow_rate,
            mass_rate,
            omitted,
        ] {
            finite(value, "phase-changing bubble response")?;
        }
        Ok(VaporBubbleObservation {
            state,
            vapor,
            wall_velocity_m_s: wall_velocity,
            liquid_acceleration_m_s2: acceleration,
            temperature_rate_k_s: temperature_rate,
            vapor_mass_kg: thermal.vapor_mass_kg,
            vapor_mass_rate_kg_s: mass_rate,
            vapor_mass_flux_kg_m2_s: mass_flux,
            volume_m3: thermal.volume_m3,
            volume_rate_m3_s: volume_rate,
            liquid_volume_flow_m3_s: liquid_flow,
            liquid_volume_flow_rate_m3_s2: liquid_flow_rate,
            liquid_interface_pressure_pa: interface_pressure,
            heat_power_w: heat,
            vapor_recoil_pressure_pa: recoil,
            neglected_enthalpy_fraction,
            neglected_specific_enthalpy_j_kg: omitted,
            flat_interface_pressure_fraction,
            frozen_liquid_relative_change,
        })
    }

    /// Integrate the native phase/flow state and independent heat, enthalpy,
    /// pressure-work and phase-mass quadratures on the SAME adaptive clock.
    /// Callbacks must be deterministic, bounded and resolved by maximum_step_s.
    pub fn advance_controlled<P, H, C>(
        &mut self,
        end_time_s: f64,
        numerics: &VaporBubbleNumerics,
        pressure: &P,
        heat: &H,
        cancelled: &mut C,
    ) -> Result<VaporBubbleStep, VaporBubbleError>
    where
        P: Fn(f64) -> f64,
        H: Fn(f64, &VaporBubbleThermalState) -> Result<f64, VaporBubbleError>,
        C: FnMut() -> bool,
    {
        Ok(self
            .prepare_advance_controlled(end_time_s, numerics, pressure, heat, cancelled)?
            .commit())
    }

    pub fn prepare_advance_controlled<'a, P, H, C>(
        &'a mut self,
        end_time_s: f64,
        numerics: &VaporBubbleNumerics,
        pressure: &P,
        heat: &H,
        cancelled: &mut C,
    ) -> Result<VaporBubbleProposal<'a>, VaporBubbleError>
    where
        P: Fn(f64) -> f64,
        H: Fn(f64, &VaporBubbleThermalState) -> Result<f64, VaporBubbleError>,
        C: FnMut() -> bool,
    {
        if cancelled() {
            return Err(VaporBubbleError::Cancelled);
        }
        if ![
            numerics.temperature_absolute_tolerance_k,
            numerics.energy_absolute_tolerance_j,
            numerics.mass_absolute_tolerance_kg,
            numerics.energy_balance_tolerance_j,
            numerics.mass_balance_tolerance_kg,
        ]
        .into_iter()
        .all(positive)
        {
            return Err(VaporBubbleError::InvalidInput(
                "positive thermal, mass and energy error budgets",
            ));
        }
        let (initial, _) = self.thermal_state(self.state)?;
        let anchor = initial.vapor.liquid_enthalpy_j_kg();
        let start = self.state;
        let result = stepper::advance_array(
            start.time_s,
            [
                start.radius_m,
                start.liquid_velocity_m_s,
                start.temperature_k,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ],
            end_time_s,
            &numerics.radial,
            [
                numerics.radial.radius_absolute_tolerance_m,
                numerics.radial.velocity_absolute_tolerance_m_s,
                numerics.temperature_absolute_tolerance_k,
                numerics.energy_absolute_tolerance_j,
                numerics.energy_absolute_tolerance_j,
                numerics.energy_absolute_tolerance_j,
                numerics.mass_absolute_tolerance_kg,
                numerics.energy_absolute_tolerance_j,
            ],
            &mut |time, values: &[f64; 8]| {
                let state = VaporBubbleState::try_new(time, values[0], values[1], values[2])?;
                let (thermal, variation) = self.thermal_state(state)?;
                let power = heat(time, &thermal)?;
                let observation =
                    self.observe_thermal(thermal, variation, pressure(time), power)?;
                let work = observation.vapor.pressure_pa() * observation.volume_rate_m3_s;
                let enthalpy = (observation.vapor.liquid_enthalpy_j_kg() - anchor)
                    * observation.vapor_mass_rate_kg_s;
                let neglected = observation.neglected_specific_enthalpy_j_kg
                    * observation.vapor_mass_rate_kg_s.abs();
                Ok([
                    observation.wall_velocity_m_s,
                    observation.liquid_acceleration_m_s2,
                    observation.temperature_rate_k_s,
                    power,
                    enthalpy,
                    work,
                    observation.vapor_mass_rate_kg_s,
                    neglected,
                ])
            },
            cancelled,
            &|error: &VaporBubbleError| {
                matches!(
                    error,
                    VaporBubbleError::OutsideModel { .. }
                        | VaporBubbleError::TemperatureOutsideDomain { .. }
                        | VaporBubbleError::Dynamics(
                            CavitationError::RadiusOutsideDomain { .. }
                                | CavitationError::NonFiniteDerived(_)
                        )
                )
            },
        )?;
        if cancelled() {
            return Err(VaporBubbleError::Cancelled);
        }
        let next = VaporBubbleState::try_new(
            result.time_s,
            result.values[0],
            result.values[1],
            result.values[2],
        )?;
        let (thermal, variation) = self.thermal_state(next)?;
        let power = heat(next.time_s, &thermal)?;
        let observation = self.observe_thermal(thermal, variation, pressure(next.time_s), power)?;
        // Factored volume and product differences preserve small physical
        // changes without subtracting two large phase-energy inventories.
        let delta_volume = (4.0 * PI / 3.0)
            * (next.radius_m - start.radius_m)
            * (next.radius_m * next.radius_m
                + next.radius_m * start.radius_m
                + start.radius_m * start.radius_m);
        let delta_mass = initial.volume_m3
            * (thermal.vapor.vapor_density_kg_m3() - initial.vapor.vapor_density_kg_m3())
            + thermal.vapor.vapor_density_kg_m3() * delta_volume;
        let mass_residual = finite(delta_mass - result.values[6], "vapor source mass balance")?;
        if mass_residual.abs() > numerics.mass_balance_tolerance_kg {
            return Err(VaporBubbleError::MassBalance {
                residual_kg: mass_residual,
                tolerance_kg: numerics.mass_balance_tolerance_kg,
            });
        }
        let anchored_change = initial.vapor_mass_kg
            * (thermal.vapor.vapor_internal_energy_j_kg()
                - initial.vapor.vapor_internal_energy_j_kg())
            + delta_mass * (thermal.vapor.vapor_internal_energy_j_kg() - anchor);
        let energy_residual = finite(
            anchored_change - result.values[3] - result.values[4] + result.values[5],
            "vapor source energy balance",
        )?;
        if energy_residual.abs() > numerics.energy_balance_tolerance_j {
            return Err(VaporBubbleError::EnergyBalance {
                residual_j: energy_residual,
                tolerance_j: numerics.energy_balance_tolerance_j,
            });
        }
        let full_change = finite(
            anchored_change + anchor * delta_mass,
            "vapor internal-energy change",
        )?;
        let full_enthalpy = finite(
            result.values[4] + anchor * result.values[6],
            "vapor transferred enthalpy",
        )?;
        let absolute_residual = finite(
            energy_residual + anchor * mass_residual,
            "original-datum vapor energy balance",
        )?;
        if !result.values[7].is_finite() || result.values[7] < 0.0 {
            return Err(VaporBubbleError::InvalidInput(
                "integrated omitted interface-energy scale",
            ));
        }
        if cancelled() {
            return Err(VaporBubbleError::Cancelled);
        }
        let step = VaporBubbleStep {
            observation,
            accepted_steps: result.accepted_steps,
            rejected_steps: result.rejected_steps,
            rhs_evaluations: result.rhs_evaluations,
            maximum_normalized_local_error: result.maximum_normalized_local_error,
            vapor_mass_change_kg: delta_mass,
            transferred_vapor_mass_kg: result.values[6],
            mass_residual_kg: mass_residual,
            vapor_energy_change_j: full_change,
            supplied_heat_j: result.values[3],
            transferred_enthalpy_j: full_enthalpy,
            pressure_work_j: result.values[5],
            reference_liquid_enthalpy_j_kg: anchor,
            reference_invariant_energy_change_j: anchored_change,
            reference_invariant_transferred_enthalpy_j: result.values[4],
            energy_residual_j: energy_residual,
            original_datum_energy_residual_j: absolute_residual,
            neglected_interface_energy_scale_j: result.values[7],
        };
        Ok(VaporBubbleProposal { target: self, step })
    }
}

#[cfg(test)]
#[path = "vapor/tests.rs"]
mod tests;
