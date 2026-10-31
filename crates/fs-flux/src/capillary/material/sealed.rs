//! Source-bound overdamped capillary filling against a closed gas pocket.
//!
//! Liquid viscosity/density and equilibrium wetting remain at their original
//! declared source state. The actual frozen-composition gas evolves through its
//! NASA-9 midpoint heat/volume owner, including pressure, temperature, caloric
//! storage and transport. The caller assigns the declared liquid/wall/exterior
//! pair; no gas identity is inferred from free-form material names.
//!
//! This is an implicit-midpoint time discretization, not an exact finite-time
//! Washburn solution. Reduce dt on boundary, source-domain or convergence
//! refusal. No end-event time or clipped motion is fabricated. Unheated steps
//! also refuse reversal of the overdamped endpoint driving pressure.
//! The gas is spatially uniform: gas through-flow resistance, inertia,
//! dynamic contact angles, liquid heating and phase change are outside this law.

use super::MaterialCapillaryTube;
use crate::capillary::sealed_mechanics as mechanics;
use core::f64::consts::PI;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::midpoint::GasHeatVolumeStep;
use fs_material::thermochemical_gas::energy::{GasEnergyError, GasEnergySolve};

/// A rigid circular tube, including any sealed downstream gas dead volume.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SealedCapillaryGeometry {
    pub radius_m: f64,
    pub tube_length_m: f64,
    pub wet_length_m: f64,
    pub dead_volume_m3: f64,
}
impl SealedCapillaryGeometry {
    fn require_volume(self, volume: f64) -> Result<(), SealedCapillaryError> {
        if ![self.radius_m, self.tube_length_m, volume]
            .into_iter()
            .all(|x| x.is_finite() && x > 0.0)
            || !self.wet_length_m.is_finite()
            || self.wet_length_m < 0.0
            || self.wet_length_m > self.tube_length_m
            || !self.dead_volume_m3.is_finite()
            || self.dead_volume_m3 < 0.0
        {
            return Err(SealedCapillaryError::Invalid(
                "positive tube/gas geometry and an in-tube meniscus required",
            ));
        }
        let area = PI * self.radius_m * self.radius_m;
        let expected = self.dead_volume_m3 + area * (self.tube_length_m - self.wet_length_m);
        if !area.is_finite()
            || area <= 0.0
            || !expected.is_finite()
            || expected <= 0.0
            || (expected - volume).abs() > 128.0 * f64::EPSILON * expected.max(volume)
        {
            return Err(SealedCapillaryError::Invalid(
                "original gas volume must equal dry tube plus declared dead volume",
            ));
        }
        Ok(())
    }
}

/// Full-step boundary conditions. Inlet pressure is ABSOLUTE, never gauge.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SealedCapillaryInput {
    pub inlet_pressure_pa: f64,
    pub duration_s: f64,
    /// Signed externally supplied thermal power into the gas [W].
    pub gas_heat_power_w: f64,
    /// Gravity component opposing filling [m/s²]; negative assists filling.
    pub opposing_gravity_m_s2: f64,
}
/// Explicit bounded work, pressure/length convergence and whole-system energy tolerances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SealedCapillarySolve {
    pub gas: GasEnergySolve,
    pub maximum_iterations: usize,
    pub pressure_absolute_tolerance_pa: f64,
    /// Maximum accepted local Newton correction [m]; pressure tolerance must also pass.
    pub length_absolute_tolerance_m: f64,
    pub energy_absolute_tolerance_j: f64,
    /// Applied to the sum of magnitudes of incremental work/storage, not absolute formation energy.
    pub energy_relative_tolerance: f64,
}
#[derive(Clone, Debug, PartialEq)]
pub enum SealedCapillaryError {
    Invalid(&'static str),
    Gas(GasEnergyError),
    Cancelled,
    /// Reduce dt or supply the physical boundary/event model; no partial step escaped.
    StepTooLarge(&'static str),
    NotConverged {
        iterations: usize,
        residual_pa: f64,
    },
    EnergyBalance {
        residual_j: f64,
        tolerance_j: f64,
    },
}
impl core::fmt::Display for SealedCapillaryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "sealed capillary: {self:?}")
    }
}
impl std::error::Error for SealedCapillaryError {}

/// Source-retaining liquid/wetting assignment plus the actual closed gas storage.
#[derive(Clone, Debug)]
pub struct SealedGasCapillary<'a> {
    tube: MaterialCapillaryTube<'a>,
    geometry: SealedCapillaryGeometry,
    gas: ClosedGasParcel,
}
/// Independent incremental accounts; gas heating appears exactly once.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SealedCapillaryStep {
    pub initial_wet_length_m: f64,
    pub wet_length_m: f64,
    pub elapsed_time_s: f64,
    pub transferred_liquid_volume_m3: f64,
    pub transferred_liquid_mass_kg: f64,
    pub inlet_pressure_work_j: f64,
    pub wall_surface_energy_change_j: f64,
    pub gravitational_energy_change_j: f64,
    pub liquid_viscous_dissipation_j: f64,
    /// Original gas owner's actual work, stored-energy increment, heat and residual.
    pub gas: GasHeatVolumeStep,
    pub midpoint_gas_pressure_pa: f64,
    pub pressure_residual_pa: f64,
    /// Inlet work + heat - gas storage - wall/gravity storage - liquid dissipation.
    pub energy_balance_residual_j: f64,
    pub iterations: usize,
}
#[derive(Clone, Debug)]
pub struct SealedCapillaryCandidate<'a> {
    endpoint: SealedGasCapillary<'a>,
    step: SealedCapillaryStep,
}
impl<'a> SealedCapillaryCandidate<'a> {
    #[must_use]
    pub const fn endpoint(&self) -> &SealedGasCapillary<'a> {
        &self.endpoint
    }
    #[must_use]
    pub const fn step(&self) -> &SealedCapillaryStep {
        &self.step
    }
    #[must_use]
    pub fn into_parts(self) -> (SealedGasCapillary<'a>, SealedCapillaryStep) {
        (self.endpoint, self.step)
    }
}
impl<'a> MaterialCapillaryTube<'a> {
    /// Bind existing closed gas storage without resampling species or copying scalar surrogates.
    pub fn seal(
        self,
        gas: ClosedGasParcel,
        geometry: SealedCapillaryGeometry,
    ) -> Result<SealedGasCapillary<'a>, SealedCapillaryError> {
        geometry.require_volume(gas.volume_m3())?;
        Ok(SealedGasCapillary {
            tube: self,
            geometry,
            gas,
        })
    }
}
impl<'a> SealedGasCapillary<'a> {
    #[must_use]
    pub const fn tube(&self) -> &MaterialCapillaryTube<'a> {
        &self.tube
    }
    #[must_use]
    pub const fn gas(&self) -> &ClosedGasParcel {
        &self.gas
    }
    #[must_use]
    pub const fn geometry(&self) -> SealedCapillaryGeometry {
        self.geometry
    }

    /// Propose one complete, immutable step. Cancellation, invalid source
    /// states, tube ends and energy refusal leave this owner exactly unchanged.
    /// Viscous loss is reported as dissipation, not implicitly injected into gas heat.
    pub fn step_candidate(
        &self,
        input: SealedCapillaryInput,
        controls: SealedCapillarySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SealedCapillaryCandidate<'a>, SealedCapillaryError> {
        controls.gas.validate().map_err(SealedCapillaryError::Gas)?;
        if !input.gas_heat_power_w.is_finite()
            || !controls.energy_absolute_tolerance_j.is_finite()
            || controls.energy_absolute_tolerance_j <= 0.0
            || !controls.energy_relative_tolerance.is_finite()
            || controls.energy_relative_tolerance < 0.0
        {
            return Err(SealedCapillaryError::Invalid(
                "finite heat power and explicit nonnegative energy tolerances required",
            ));
        }
        let heat = input.gas_heat_power_w * input.duration_s;
        if !heat.is_finite()
            || (input.gas_heat_power_w != 0.0 && input.duration_s > 0.0 && heat == 0.0)
        {
            return Err(SealedCapillaryError::Invalid(
                "full-step heat is unrepresentable",
            ));
        }
        let g = self.geometry;
        let density = self.tube.liquid.density_kg_m3();
        let capillary_pressure =
            2.0 * self.tube.surface_tension_n_m * self.tube.contact_angle_cosine / g.radius_m;
        let motion = mechanics::solve(
            mechanics::Input {
                radius: g.radius_m,
                length: g.tube_length_m,
                initial: g.wet_length_m,
                gas_volume: self.gas.volume_m3(),
                viscosity: self.tube.liquid.dynamic_viscosity_pa_s(),
                density,
                capillary_pressure,
                inlet_pressure: input.inlet_pressure_pa,
                gravity: input.opposing_gravity_m_s2,
                duration: input.duration_s,
            },
            mechanics::Controls {
                iterations: controls.maximum_iterations,
                pressure_tolerance: controls.pressure_absolute_tolerance_pa,
                length_tolerance: controls.length_absolute_tolerance_m,
            },
            |volume, cancel| {
                let gas = self
                    .gas
                    .midpoint_heat_volume(volume, heat, controls.gas, cancel)?;
                Ok(mechanics::Pressure {
                    value: gas.midpoint().state().pressure,
                    derivative: gas.step().pressure_final_volume_derivative_pa_m3,
                    payload: gas,
                })
            },
            &mut cancelled,
        )
        .map_err(|error| match error {
            mechanics::Error::Invalid(s) => SealedCapillaryError::Invalid(s),
            mechanics::Error::Gas(GasEnergyError::Cancelled) | mechanics::Error::Cancelled => {
                SealedCapillaryError::Cancelled
            }
            mechanics::Error::Gas(e) => SealedCapillaryError::Gas(e),
            mechanics::Error::Boundary => {
                SealedCapillaryError::StepTooLarge("tube/gas boundary reached; reduce dt")
            }
            mechanics::Error::NonMonotone => SealedCapillaryError::StepTooLarge(
                "cannot certify a stable increasing residual; reduce dt",
            ),
            mechanics::Error::NotConverged {
                iterations,
                residual,
            } => SealedCapillaryError::NotConverged {
                iterations,
                residual_pa: residual,
            },
        })?;
        let midpoint_pressure = motion.pressure.value;
        let (gas, gas_step) = motion.pressure.payload.into_parts();
        let geometry = SealedCapillaryGeometry {
            wet_length_m: motion.length,
            ..g
        };
        geometry.require_volume(gas.volume_m3())?;
        let endpoint_drive = (input.inlet_pressure_pa - gas.state().state().pressure)
            + capillary_pressure
            - density * input.opposing_gravity_m_s2 * motion.length;
        if heat == 0.0
            && ((motion.volume > 0.0 && endpoint_drive < -controls.pressure_absolute_tolerance_pa)
                || (motion.volume < 0.0
                    && endpoint_drive > controls.pressure_absolute_tolerance_pa))
        {
            return Err(SealedCapillaryError::StepTooLarge(
                "unheated overdamped equilibrium overshoot; reduce dt",
            ));
        }
        let mass = density * motion.volume;
        let residual = (motion.inlet_work + heat - gas_step.stored_energy_change_j)
            - motion.wall_energy
            - motion.gravity_energy
            - motion.dissipation;
        let scale = motion.inlet_work.abs()
            + heat.abs()
            + gas_step.stored_energy_change_j.abs()
            + motion.wall_energy.abs()
            + motion.gravity_energy.abs()
            + motion.dissipation;
        let tolerance =
            controls.energy_absolute_tolerance_j + controls.energy_relative_tolerance * scale;
        if !mass.is_finite()
            || (motion.volume != 0.0 && mass == 0.0)
            || !residual.is_finite()
            || !tolerance.is_finite()
        {
            return Err(SealedCapillaryError::Invalid(
                "mass or energy accounting is unrepresentable",
            ));
        }
        if residual.abs() > tolerance {
            return Err(SealedCapillaryError::EnergyBalance {
                residual_j: residual,
                tolerance_j: tolerance,
            });
        }
        if cancelled() {
            return Err(SealedCapillaryError::Cancelled);
        }
        Ok(SealedCapillaryCandidate {
            endpoint: SealedGasCapillary {
                tube: self.tube,
                geometry,
                gas,
            },
            step: SealedCapillaryStep {
                initial_wet_length_m: g.wet_length_m,
                wet_length_m: motion.length,
                elapsed_time_s: input.duration_s,
                transferred_liquid_volume_m3: motion.volume,
                transferred_liquid_mass_kg: mass,
                inlet_pressure_work_j: motion.inlet_work,
                wall_surface_energy_change_j: motion.wall_energy,
                gravitational_energy_change_j: motion.gravity_energy,
                liquid_viscous_dissipation_j: motion.dissipation,
                gas: gas_step,
                midpoint_gas_pressure_pa: midpoint_pressure,
                pressure_residual_pa: motion.pressure_residual,
                energy_balance_residual_j: residual,
                iterations: motion.iterations,
            },
        })
    }
}
