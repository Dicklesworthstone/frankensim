//! Thermal gas diffusion in the original dilute acoustic populations.
//!
//! Each prepared native bubble retains its actual equilibrium radius and the
//! gas material snapshot evaluated at that gas pressure and interface temperature.
//! The gas diffusion equation supplies complex stiffness and thermal loss; the
//! original acoustic owner supplies all volume weights, scattering, physical
//! admission and duct propagation. Nothing is fitted to a target attenuation.
//!
//! Preparation borrows an immutable accepted population. A changed operating
//! pressure requires new gas snapshots from the original material law at each
//! continued class's gas pressure, followed by new thermal preparation. Stale
//! pressure/temperature coefficients refuse. No source is silently rescaled.
//!
//! The default interface is isothermal. `new_with_conducting_liquid` instead
//! resolves spherical heat diffusion in the surrounding liquid from the SAME
//! material source. Its interface temperature and gas/liquid dissipation feed
//! these original acoustic/population kernels without fitted thermal damping.
//! The independent liquid thermal layers must be small compared with the mean
//! bubble spacing; this model does not solve interbubble thermal interactions.
//!
//! This is a small-signal, ideal noncondensable gas with a homobaric interior,
//! not nonlinear thermal collapse or vapor kinetics.
//! The native gas-acoustic compactness, density and vapor-fraction limits apply
//! in addition to every original bubbly-acoustic limit. Gas continuum validity,
//! frozen liquid mechanics and negligible vapor heat/mass transfer remain
//! caller conditions; the interface temperature is solved in the conducting lane.
//!
//! Thermal and liquid-viscous absorption are reported separately from coherent
//! scattering attenuation. A local pressure amplitude gives local cycle-mean
//! heating; a pressure envelope alone is not a spatially integrated heat source.

use super::population::poll;
use super::*;
use fs_flux::cavitation::thermal::{
    ConductingLiquidBubble, IsothermalWallBubble, ThermalBubbleLimits, ThermalBubbleResponse,
};
use fs_material::fluid_transport::{ThermalFluidTransport, ThermalFluidTransportState};
use fs_material::gas::GasState;

/// Admission of separate spherical liquid thermal layers in the dilute medium.
/// These conservative ceilings are model boundaries, not calibrated errors.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConductingLiquidPopulationLimits {
    /// Largest sqrt(2 D_l/omega) divided by the mean bubble spacing, in (0,0.1].
    /// The exact DC limit has no oscillating interface temperature or heat flux.
    pub maximum_penetration_to_spacing: f64,
    /// Largest peak interface-temperature perturbation divided by T, in (0,0.05].
    pub maximum_fractional_interface_temperature: f64,
}

impl ConductingLiquidPopulationLimits {
    fn validate(self) -> Result<(), BubblyLiquidError> {
        for (name, value, maximum) in [
            (
                "liquid thermal penetration limit",
                self.maximum_penetration_to_spacing,
                0.1,
            ),
            (
                "interface temperature limit",
                self.maximum_fractional_interface_temperature,
                0.05,
            ),
        ] {
            if !value.is_finite() || value <= 0.0 || value > maximum {
                return Err(BubblyLiquidError::InvalidInput(name));
            }
        }
        Ok(())
    }
}

#[derive(Clone, Copy)]
enum PreparedBubble {
    Isothermal(IsothermalWallBubble),
    Conducting(ConductingLiquidBubble),
}

impl PreparedBubble {
    fn new(
        base: IsothermalWallBubble,
        liquid: Option<ThermalFluidTransportState>,
    ) -> Result<Self, fs_flux::cavitation::thermal::ThermalBubbleError> {
        match liquid {
            Some(liquid) => Ok(Self::Conducting(ConductingLiquidBubble::new(base, liquid)?)),
            None => Ok(Self::Isothermal(base)),
        }
    }

    fn base(&self) -> &IsothermalWallBubble {
        match self {
            Self::Isothermal(base) => base,
            Self::Conducting(model) => model.base(),
        }
    }

    fn response(&self, omega: f64) -> Result<PreparedResponse, BubblyLiquidError> {
        match self {
            Self::Isothermal(base) => {
                let bubble = base.response(omega)?;
                Ok(PreparedResponse {
                    gas_power_w: bubble.thermal_power_per_fractional_radius_squared_w,
                    bubble,
                    liquid_power_w: 0.0,
                    interface_temperature_per_fractional_radius_k: C64::ZERO,
                    liquid_diffusion_parameter: 0.0,
                })
            }
            Self::Conducting(model) => {
                let response = model.response(omega)?;
                Ok(PreparedResponse {
                    bubble: response.bubble,
                    gas_power_w: response.gas_thermal_power_per_fractional_radius_squared_w,
                    liquid_power_w: response.liquid_thermal_power_per_fractional_radius_squared_w,
                    interface_temperature_per_fractional_radius_k: response
                        .interface_temperature_per_fractional_radius_k,
                    liquid_diffusion_parameter: response.liquid_thermal_diffusion_parameter,
                })
            }
        }
    }
}

#[derive(Clone, Copy)]
struct PreparedResponse {
    bubble: ThermalBubbleResponse,
    gas_power_w: f64,
    liquid_power_w: f64,
    interface_temperature_per_fractional_radius_k: C64,
    liquid_diffusion_parameter: f64,
}

#[derive(Clone, Copy)]
struct LiquidBoundary {
    material: ThermalFluidTransportState,
    limits: ConductingLiquidPopulationLimits,
}

impl LiquidBoundary {
    fn admit_material(self, original: CavitationState) -> Result<(), BubblyLiquidError> {
        let transport = self.material.transport();
        if self.material.temperature_k() != original.temperature_k()
            || transport.density_kg_m3() != original.liquid_density_kg_m3()
            || transport.dynamic_viscosity_pa_s() != original.dynamic_viscosity_pa_s()
        {
            return Err(BubblyLiquidError::MaterialMismatch);
        }
        Ok(())
    }

    fn admit_spacing(self, omega: f64, number: f64) -> Result<(), BubblyLiquidError> {
        if omega == 0.0 || number == 0.0 {
            return Ok(());
        }
        let spacing = det::pow(number, -1.0 / 3.0);
        let penetration =
            det::sqrt(2.0 * self.material.transport().thermal_diffusivity_m2_s() / omega);
        bound(
            "liquid thermal penetration / mean spacing",
            penetration / spacing,
            self.limits.maximum_penetration_to_spacing,
        )
    }

    fn admit_temperature(
        self,
        response: PreparedResponse,
        single: SingleBubbleResponse,
    ) -> Result<f64, BubblyLiquidError> {
        let fraction = response.interface_temperature_per_fractional_radius_k.abs()
            * single.fractional_radius_amplitude
            / self.material.temperature_k();
        bound(
            "fractional interface temperature",
            fraction,
            self.limits.maximum_fractional_interface_temperature,
        )?;
        Ok(fraction)
    }
}

/// Native gas diffusion prepared once for a monodisperse acoustic owner.
pub struct ThermalDiluteBubblyLiquid<'p, 'm, F: CavitationProperties + ?Sized> {
    original: &'p DiluteBubblyLiquid<'m, F>,
    thermal: Option<PreparedBubble>,
    liquid: Option<LiquidBoundary>,
}

impl<'p, 'm, F: CavitationProperties + ?Sized> ThermalDiluteBubblyLiquid<'p, 'm, F> {
    /// An occupied population needs its current gas snapshot; exact zero void
    /// fraction needs `None`, without manufacturing gas material or an oscillator.
    pub fn new(
        original: &'p DiluteBubblyLiquid<'m, F>,
        gas: Option<GasState>,
        limits: ThermalBubbleLimits,
    ) -> Result<Self, BubblyLiquidError> {
        Self::prepare(original, gas, limits, None)
    }

    /// Resolve liquid heat capacity, conductivity and temperature from the
    /// original cavitation source. Density and viscosity must agree exactly
    /// with its retained liquid dynamics. No material is selected by name.
    pub fn new_with_conducting_liquid(
        original: &'p DiluteBubblyLiquid<'m, F>,
        gas: Option<GasState>,
        limits: ThermalBubbleLimits,
        liquid_limits: ConductingLiquidPopulationLimits,
    ) -> Result<Self, BubblyLiquidError>
    where
        F: ThermalFluidTransport,
    {
        liquid_limits.validate()?;
        let liquid = LiquidBoundary {
            material: original.source().thermal_transport_state()?,
            limits: liquid_limits,
        };
        Self::prepare(original, gas, limits, Some(liquid))
    }

    fn prepare(
        original: &'p DiluteBubblyLiquid<'m, F>,
        gas: Option<GasState>,
        limits: ThermalBubbleLimits,
        liquid: Option<LiquidBoundary>,
    ) -> Result<Self, BubblyLiquidError> {
        limits.validate()?;
        if let Some(liquid) = liquid {
            liquid.admit_material(original.material())?;
        }
        let thermal = match (original.linearization(), gas) {
            (Some(linear), Some(gas)) => Some(PreparedBubble::new(
                IsothermalWallBubble::new(
                    *original.bubble(),
                    linear.equilibrium_radius_m,
                    gas,
                    limits,
                )?,
                liquid.map(|boundary| boundary.material),
            )?),
            (None, None) => None,
            _ => {
                return Err(BubblyLiquidError::InvalidInput(
                    "thermal gas snapshot must match occupied bubble population",
                ));
            }
        };
        Ok(Self {
            original,
            thermal,
            liquid,
        })
    }

    pub const fn original(&self) -> &'p DiluteBubblyLiquid<'m, F> {
        self.original
    }
    pub fn thermal_bubble(&self) -> Option<&IsothermalWallBubble> {
        self.thermal.as_ref().map(PreparedBubble::base)
    }
    /// Same-source liquid thermal snapshot, or the imposed isothermal wall.
    pub fn liquid_thermal_state(&self) -> Option<ThermalFluidTransportState> {
        self.liquid.map(|boundary| boundary.material)
    }

    /// No source queries, root solves or allocations occur in a frequency query.
    pub fn response(
        &self,
        omega_rad_s: f64,
        maximum_local_pressure_amplitude_pa: f64,
    ) -> Result<ThermalBubblyResponse<BubblyLiquidResponse>, BubblyLiquidError> {
        let mut losses = Losses::default();
        let acoustic = self.original.response_with_kernel(
            omega_rad_s,
            maximum_local_pressure_amplitude_pa,
            &mut |linear, medium| {
                let model = self
                    .thermal
                    .as_ref()
                    .ok_or(BubblyLiquidError::InvalidInput(
                        "missing occupied thermal bubble",
                    ))?;
                let response = model.response(omega_rad_s)?;
                let thermal = response.bubble;
                let single = admit_thermal_response(
                    model.base(),
                    thermal,
                    linear,
                    medium,
                    omega_rad_s,
                    maximum_local_pressure_amplitude_pa,
                    self.original.limits,
                )?;
                let radius = linear.equilibrium_radius_m;
                let number = self.original.void_fraction()
                    / ((4.0 * core::f64::consts::PI / 3.0) * radius * radius * radius);
                let temperature_fraction = if let Some(liquid) = self.liquid {
                    liquid.admit_spacing(omega_rad_s, number)?;
                    liquid.admit_temperature(response, single)?
                } else {
                    0.0
                };
                losses.add(
                    number,
                    model.base(),
                    response,
                    single,
                    omega_rad_s,
                    temperature_fraction,
                )?;
                Ok(single)
            },
        )?;
        Ok(ThermalBubblyResponse {
            acoustic,
            diagnostics: losses.finish()?,
            maximum_local_pressure_amplitude_pa,
        })
    }
}

/// Per-class gas diffusion borrowed from the same immutable native distribution.
/// The prepared vector is allocated once; frequency queries allocate nothing.
pub struct ThermalPolydisperseBubblyLiquid<'p, 'm, F: CavitationProperties + ?Sized> {
    original: &'p PolydisperseBubblyLiquid<'m, F>,
    thermal: Vec<Option<PreparedBubble>>,
    liquid: Option<LiquidBoundary>,
}

impl<'p, 'm, F: CavitationProperties + ?Sized> ThermalPolydisperseBubblyLiquid<'p, 'm, F> {
    /// Snapshots use original class order, with `None` for each empty bin.
    /// Validate every class and cancellation before publishing the preparation.
    /// Failure preserves the borrowed acoustic owner and its original sources.
    pub fn new(
        original: &'p PolydisperseBubblyLiquid<'m, F>,
        gas_states: &[Option<GasState>],
        limits: ThermalBubbleLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, BubblyLiquidError> {
        Self::prepare(original, gas_states, limits, None, cancelled)
    }

    /// Prepare each ORIGINAL class with one coherent thermal view of the
    /// retained liquid source. Class radii and gas pressures remain distinct.
    pub fn new_with_conducting_liquid(
        original: &'p PolydisperseBubblyLiquid<'m, F>,
        gas_states: &[Option<GasState>],
        limits: ThermalBubbleLimits,
        liquid_limits: ConductingLiquidPopulationLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, BubblyLiquidError>
    where
        F: ThermalFluidTransport,
    {
        poll(cancelled)?;
        liquid_limits.validate()?;
        let liquid = LiquidBoundary {
            material: original.source().thermal_transport_state()?,
            limits: liquid_limits,
        };
        Self::prepare(original, gas_states, limits, Some(liquid), cancelled)
    }

    fn prepare(
        original: &'p PolydisperseBubblyLiquid<'m, F>,
        gas_states: &[Option<GasState>],
        limits: ThermalBubbleLimits,
        liquid: Option<LiquidBoundary>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, BubblyLiquidError> {
        poll(cancelled)?;
        limits.validate()?;
        if let Some(liquid) = liquid {
            liquid.admit_material(original.material())?;
        }
        if gas_states.len() != original.population().classes().len() {
            return Err(BubblyLiquidError::InvalidInput(
                "one thermal gas slot per original population class required",
            ));
        }
        let mut thermal = Vec::with_capacity(gas_states.len());
        for (class_index, (class, gas)) in original
            .population()
            .classes()
            .iter()
            .zip(gas_states)
            .enumerate()
        {
            poll(cancelled)?;
            let model = match (class.linearization(), gas) {
                (Some(linear), Some(gas)) => Some(
                    IsothermalWallBubble::new(
                        *class.bubble(),
                        linear.equilibrium_radius_m,
                        *gas,
                        limits,
                    )
                    .and_then(|base| {
                        PreparedBubble::new(base, liquid.map(|boundary| boundary.material))
                    })
                    .map_err(|cause| BubblyLiquidError::ThermalClass { class_index, cause })?,
                ),
                (None, None) => None,
                _ => {
                    return Err(BubblyLiquidError::InvalidInput(
                        "thermal gas slots must match occupied population classes",
                    ));
                }
            };
            thermal.push(model);
        }
        poll(cancelled)?;
        Ok(Self {
            original,
            thermal,
            liquid,
        })
    }

    pub const fn original(&self) -> &'p PolydisperseBubblyLiquid<'m, F> {
        self.original
    }
    /// Original native gas preparations in class order; no allocation or copy.
    pub fn thermal_bubbles(&self) -> impl ExactSizeIterator<Item = Option<&IsothermalWallBubble>> {
        self.thermal
            .iter()
            .map(|model| model.as_ref().map(PreparedBubble::base))
    }
    pub fn liquid_thermal_state(&self) -> Option<ThermalFluidTransportState> {
        self.liquid.map(|boundary| boundary.material)
    }

    pub fn response(
        &self,
        omega_rad_s: f64,
        maximum_local_pressure_amplitude_pa: f64,
    ) -> Result<ThermalBubblyResponse<PolydisperseBubblyLiquidResponse>, BubblyLiquidError> {
        self.response_controlled(
            omega_rad_s,
            maximum_local_pressure_amplitude_pa,
            &mut || false,
        )
    }

    /// Native class order and compensated sums are shared with the original
    /// polytropic path. Cancellation or any class refusal publishes no response.
    pub fn response_controlled(
        &self,
        omega_rad_s: f64,
        maximum_local_pressure_amplitude_pa: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ThermalBubblyResponse<PolydisperseBubblyLiquidResponse>, BubblyLiquidError> {
        poll(cancelled)?;
        if let Some(liquid) = self.liquid {
            liquid.admit_spacing(
                omega_rad_s,
                self.original.population().number_density_per_m3(),
            )?;
        }
        let mut losses = Losses::default();
        let acoustic =
            self.original.response_with_kernel(
                omega_rad_s,
                maximum_local_pressure_amplitude_pa,
                cancelled,
                &mut |class_index, class, linear, medium, limits| {
                    let model = self.thermal[class_index].as_ref().ok_or(
                        BubblyLiquidError::InvalidInput("missing occupied thermal class"),
                    )?;
                    let response = model.response(omega_rad_s).map_err(|error| match error {
                        BubblyLiquidError::Thermal(cause) => {
                            BubblyLiquidError::ThermalClass { class_index, cause }
                        }
                        other => other,
                    })?;
                    let thermal = response.bubble;
                    let single = admit_thermal_response(
                        model.base(),
                        thermal,
                        linear,
                        medium,
                        omega_rad_s,
                        maximum_local_pressure_amplitude_pa,
                        limits,
                    )?;
                    let temperature_fraction = if let Some(liquid) = self.liquid {
                        liquid.admit_temperature(response, single)?
                    } else {
                        0.0
                    };
                    losses.add(
                        class.number_density_per_m3(),
                        model.base(),
                        response,
                        single,
                        omega_rad_s,
                        temperature_fraction,
                    )?;
                    Ok(single)
                },
            )?;
        Ok(ThermalBubblyResponse {
            acoustic,
            diagnostics: losses.finish()?,
            maximum_local_pressure_amplitude_pa,
        })
    }
}

fn admit_thermal_response(
    model: &IsothermalWallBubble,
    thermal: ThermalBubbleResponse,
    linear: BubbleLinearization,
    medium: FluidAcousticState,
    omega: f64,
    amplitude: f64,
    limits: BubblyLiquidLimits,
) -> Result<SingleBubbleResponse, BubblyLiquidError> {
    admit_single_bubble_response(
        &model.bubble(),
        linear,
        medium,
        omega,
        amplitude,
        limits,
        SingleBubbleKernel {
            pressure_transfer_s2: thermal.pressure_transfer_s2,
            wall_stiffness_pa: thermal.wall_stiffness_pa,
        },
    )
}

/// Constitutive heat coefficients for peak (not RMS) local harmonic pressure.
/// They exclude coherent scattering, which redistributes acoustic energy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalBubblyDiagnostics {
    /// Combined gas and liquid conductive dissipation.
    pub thermal_loss_per_pressure_squared_w_m3_pa2: f64,
    pub gas_thermal_loss_per_pressure_squared_w_m3_pa2: f64,
    pub liquid_thermal_loss_per_pressure_squared_w_m3_pa2: f64,
    pub viscous_loss_per_pressure_squared_w_m3_pa2: f64,
    pub maximum_gas_wavenumber_radius: f64,
    pub maximum_thermal_diffusion_parameter: f64,
    pub maximum_liquid_thermal_diffusion_parameter: f64,
    pub maximum_fractional_interface_temperature: f64,
}

/// Local cycle-mean dissipated power per original mixture volume.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BubblyAbsorption {
    pub gas_thermal_w_m3: f64,
    pub liquid_thermal_w_m3: f64,
    pub liquid_viscous_w_m3: f64,
}

/// Original admitted propagation plus separate physical absorption.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalBubblyResponse<R> {
    acoustic: R,
    diagnostics: ThermalBubblyDiagnostics,
    maximum_local_pressure_amplitude_pa: f64,
}

impl<R> ThermalBubblyResponse<R> {
    /// Use this response's native segment-wave and duct-chain operations.
    pub const fn acoustic(&self) -> &R {
        &self.acoustic
    }
    pub const fn diagnostics(&self) -> ThermalBubblyDiagnostics {
        self.diagnostics
    }

    /// Evaluate heat only at a known local peak pressure within the admitted
    /// envelope. A duct field's actual amplitudes belong to its solved field;
    /// the input envelope alone does not prescribe their spatial distribution.
    pub fn mean_absorption_w_m3(
        &self,
        local_pressure_amplitude_pa: f64,
    ) -> Result<BubblyAbsorption, BubblyLiquidError> {
        bound(
            "local pressure amplitude",
            local_pressure_amplitude_pa,
            self.maximum_local_pressure_amplitude_pa,
        )?;
        let gas = self
            .diagnostics
            .gas_thermal_loss_per_pressure_squared_w_m3_pa2
            * local_pressure_amplitude_pa
            * local_pressure_amplitude_pa;
        let liquid = self
            .diagnostics
            .liquid_thermal_loss_per_pressure_squared_w_m3_pa2
            * local_pressure_amplitude_pa
            * local_pressure_amplitude_pa;
        let viscous = self.diagnostics.viscous_loss_per_pressure_squared_w_m3_pa2
            * local_pressure_amplitude_pa
            * local_pressure_amplitude_pa;
        nonnegative(gas)?;
        nonnegative(liquid)?;
        nonnegative(viscous)?;
        Ok(BubblyAbsorption {
            gas_thermal_w_m3: gas,
            liquid_thermal_w_m3: liquid,
            liquid_viscous_w_m3: viscous,
        })
    }
}

#[derive(Default)]
struct Losses {
    thermal: population::Sum,
    gas: population::Sum,
    liquid: population::Sum,
    viscous: population::Sum,
    gas_size: f64,
    diffusion: f64,
    liquid_diffusion: f64,
    temperature_fraction: f64,
}

fn nonnegative(value: f64) -> Result<(), BubblyLiquidError> {
    if !value.is_finite() || value < 0.0 {
        return Err(BubblyLiquidError::InvalidInput(
            "unrepresentable passive bubble absorption",
        ));
    }
    Ok(())
}

impl Losses {
    fn add(
        &mut self,
        number: f64,
        model: &IsothermalWallBubble,
        response: PreparedResponse,
        single: SingleBubbleResponse,
        omega: f64,
        temperature_fraction: f64,
    ) -> Result<(), BubblyLiquidError> {
        let thermal = response.bubble;
        let radius = model.linearization().equilibrium_radius_m;
        let magnitude = single.susceptibility.abs();
        let thermal_loss =
            number * thermal.thermal_power_per_fractional_radius_squared_w * magnitude * magnitude;
        let gas_loss = number * response.gas_power_w * magnitude * magnitude;
        let liquid_loss = number * response.liquid_power_w * magnitude * magnitude;
        // <16*pi*mu*R*Rdot^2> = 8*pi*mu*R^3*omega^2*|delta R/R|^2.
        let viscous_loss = number
            * (8.0 * core::f64::consts::PI)
            * model.bubble().liquid().dynamic_viscosity_pa_s()
            * radius
            * radius
            * radius
            * omega
            * omega
            * magnitude
            * magnitude;
        nonnegative(thermal_loss)?;
        nonnegative(gas_loss)?;
        nonnegative(liquid_loss)?;
        nonnegative(viscous_loss)?;
        self.thermal.add(thermal_loss);
        self.gas.add(gas_loss);
        self.liquid.add(liquid_loss);
        self.viscous.add(viscous_loss);
        self.gas_size = self.gas_size.max(thermal.gas_wavenumber_radius);
        self.diffusion = self.diffusion.max(thermal.thermal_diffusion_parameter);
        self.liquid_diffusion = self
            .liquid_diffusion
            .max(response.liquid_diffusion_parameter);
        self.temperature_fraction = self.temperature_fraction.max(temperature_fraction);
        Ok(())
    }

    fn finish(self) -> Result<ThermalBubblyDiagnostics, BubblyLiquidError> {
        let thermal = self.thermal.value();
        let gas = self.gas.value();
        let liquid = self.liquid.value();
        let viscous = self.viscous.value();
        nonnegative(thermal)?;
        nonnegative(gas)?;
        nonnegative(liquid)?;
        nonnegative(viscous)?;
        Ok(ThermalBubblyDiagnostics {
            thermal_loss_per_pressure_squared_w_m3_pa2: thermal,
            gas_thermal_loss_per_pressure_squared_w_m3_pa2: gas,
            liquid_thermal_loss_per_pressure_squared_w_m3_pa2: liquid,
            viscous_loss_per_pressure_squared_w_m3_pa2: viscous,
            maximum_gas_wavenumber_radius: self.gas_size,
            maximum_thermal_diffusion_parameter: self.diffusion,
            maximum_liquid_thermal_diffusion_parameter: self.liquid_diffusion,
            maximum_fractional_interface_temperature: self.temperature_fraction,
        })
    }
}

#[cfg(test)]
#[path = "thermal/tests.rs"]
mod tests;
