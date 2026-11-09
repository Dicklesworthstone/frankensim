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
//! This is a small-signal, ideal noncondensable gas with an isothermal interface
//! and homobaric interior, not nonlinear thermal collapse or vapor kinetics.
//! The native gas-acoustic compactness, density and vapor-fraction limits apply
//! in addition to every original bubbly-acoustic limit. Gas continuum validity
//! and an effectively isothermal liquid interface remain caller assumptions.
//!
//! Thermal and liquid-viscous absorption are reported separately from coherent
//! scattering attenuation. A local pressure amplitude gives local cycle-mean
//! heating; a pressure envelope alone is not a spatially integrated heat source.

use super::population::poll;
use super::*;
use fs_flux::cavitation::thermal::{
    IsothermalWallBubble, ThermalBubbleLimits, ThermalBubbleResponse,
};
use fs_material::gas::GasState;

/// Native gas diffusion prepared once for a monodisperse acoustic owner.
pub struct ThermalDiluteBubblyLiquid<'p, 'm, F: CavitationProperties + ?Sized> {
    original: &'p DiluteBubblyLiquid<'m, F>,
    thermal: Option<IsothermalWallBubble>,
}

impl<'p, 'm, F: CavitationProperties + ?Sized> ThermalDiluteBubblyLiquid<'p, 'm, F> {
    /// An occupied population needs its current gas snapshot; exact zero void
    /// fraction needs `None`, without manufacturing gas material or an oscillator.
    pub fn new(
        original: &'p DiluteBubblyLiquid<'m, F>,
        gas: Option<GasState>,
        limits: ThermalBubbleLimits,
    ) -> Result<Self, BubblyLiquidError> {
        limits.validate()?;
        let thermal = match (original.linearization(), gas) {
            (Some(linear), Some(gas)) => Some(IsothermalWallBubble::new(
                *original.bubble(),
                linear.equilibrium_radius_m,
                gas,
                limits,
            )?),
            (None, None) => None,
            _ => {
                return Err(BubblyLiquidError::InvalidInput(
                    "thermal gas snapshot must match occupied bubble population",
                ));
            }
        };
        Ok(Self { original, thermal })
    }

    pub const fn original(&self) -> &'p DiluteBubblyLiquid<'m, F> {
        self.original
    }
    pub const fn thermal_bubble(&self) -> Option<&IsothermalWallBubble> {
        self.thermal.as_ref()
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
                let thermal = model.response(omega_rad_s)?;
                let single = admit_thermal_response(
                    model,
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
                losses.add(number, model, thermal, single, omega_rad_s)?;
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
    thermal: Vec<Option<IsothermalWallBubble>>,
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
        poll(cancelled)?;
        limits.validate()?;
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
        Ok(Self { original, thermal })
    }

    pub const fn original(&self) -> &'p PolydisperseBubblyLiquid<'m, F> {
        self.original
    }
    pub fn thermal_bubbles(&self) -> &[Option<IsothermalWallBubble>] {
        &self.thermal
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
                    let thermal = model
                        .response(omega_rad_s)
                        .map_err(|cause| BubblyLiquidError::ThermalClass { class_index, cause })?;
                    let single = admit_thermal_response(
                        model,
                        thermal,
                        linear,
                        medium,
                        omega_rad_s,
                        maximum_local_pressure_amplitude_pa,
                        limits,
                    )?;
                    losses.add(
                        class.number_density_per_m3(),
                        model,
                        thermal,
                        single,
                        omega_rad_s,
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
    pub thermal_loss_per_pressure_squared_w_m3_pa2: f64,
    pub viscous_loss_per_pressure_squared_w_m3_pa2: f64,
    pub maximum_gas_wavenumber_radius: f64,
    pub maximum_thermal_diffusion_parameter: f64,
}

/// Local cycle-mean dissipated power per original mixture volume.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BubblyAbsorption {
    pub gas_thermal_w_m3: f64,
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
        let thermal = self.diagnostics.thermal_loss_per_pressure_squared_w_m3_pa2
            * local_pressure_amplitude_pa
            * local_pressure_amplitude_pa;
        let viscous = self.diagnostics.viscous_loss_per_pressure_squared_w_m3_pa2
            * local_pressure_amplitude_pa
            * local_pressure_amplitude_pa;
        nonnegative(thermal)?;
        nonnegative(viscous)?;
        Ok(BubblyAbsorption {
            gas_thermal_w_m3: thermal,
            liquid_viscous_w_m3: viscous,
        })
    }
}

#[derive(Default)]
struct Losses {
    thermal: population::Sum,
    viscous: population::Sum,
    gas_size: f64,
    diffusion: f64,
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
        thermal: ThermalBubbleResponse,
        single: SingleBubbleResponse,
        omega: f64,
    ) -> Result<(), BubblyLiquidError> {
        let radius = model.linearization().equilibrium_radius_m;
        let magnitude = single.susceptibility.abs();
        let thermal_loss =
            number * thermal.thermal_power_per_fractional_radius_squared_w * magnitude * magnitude;
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
        nonnegative(viscous_loss)?;
        self.thermal.add(thermal_loss);
        self.viscous.add(viscous_loss);
        self.gas_size = self.gas_size.max(thermal.gas_wavenumber_radius);
        self.diffusion = self.diffusion.max(thermal.thermal_diffusion_parameter);
        Ok(())
    }

    fn finish(self) -> Result<ThermalBubblyDiagnostics, BubblyLiquidError> {
        let thermal = self.thermal.value();
        let viscous = self.viscous.value();
        nonnegative(thermal)?;
        nonnegative(viscous)?;
        Ok(ThermalBubblyDiagnostics {
            thermal_loss_per_pressure_squared_w_m3_pa2: thermal,
            viscous_loss_per_pressure_squared_w_m3_pa2: viscous,
            maximum_gas_wavenumber_radius: self.gas_size,
            maximum_thermal_diffusion_parameter: self.diffusion,
        })
    }
}

#[cfg(test)]
#[path = "thermal/tests.rs"]
mod tests;
