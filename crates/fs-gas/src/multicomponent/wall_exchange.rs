//! Stationary-wall momentum loss and conservative heat exchange.
//!
//! This source uses original mixture energy and Cv, preserves species, and
//! retains total gas energy during drag. Thermal conductance is h*A [W/K].
//! Linear/Darcy coefficients are explicit frozen engineering inputs; no
//! unsteady shock boundary-layer or universal heat-transfer closure is implied.

use super::{Cell, Error as FlowError, MixtureGas, positive};
use fs_math::det;
use fs_qty::{DynViscosity, Length, Time};

/// Original thermal laws in finite gas-wall storage, including latent plateaus.
mod caloric;
pub use caloric::{MaterialWallCaloric, WallCaloricConstraint, exchange_cell_with_heat_storage};

/// Stationary radial wall drag; its lost kinetic energy stays in the gas.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum WallDrag {
    /// No radial momentum sink.
    None,
    /// Frozen exponential relaxation, dm/dt = -rate*m.
    Linear {
        /// Nonnegative momentum relaxation rate [1/s].
        rate_per_s: f64,
    },
    /// Frozen Darcy-Weisbach factor, du/dt = -f*u*abs(u)/(2*Dh).
    Darcy {
        /// Nonnegative Darcy factor, four times the Fanning factor.
        friction_factor: f64,
        /// Positive hydraulic diameter [m].
        hydraulic_diameter_m: f64,
    },
}
impl WallDrag {
    /// Frozen fully developed circular Poiseuille relaxation, 32*mu/(rho*D²).
    ///
    /// The caller supplies source-consistent viscosity and density and owns
    /// the laminar/developed-flow applicability decision. This is not a model
    /// of a transient shock boundary layer.
    pub fn circular_laminar(
        viscosity: DynViscosity,
        density_kg_m3: f64,
        diameter: Length,
    ) -> Result<Self, WallExchangeError> {
        if !nonnegative(viscosity.value())
            || !positive(density_kg_m3)
            || !positive(diameter.value())
        {
            return Err(WallExchangeError::Invalid("laminar wall coefficient inputs"));
        }
        let rate = 32.0 * (viscosity.value() / density_kg_m3)
            / diameter.value() / diameter.value();
        if !nonnegative(rate) {
            return Err(WallExchangeError::Invalid("unrepresentable laminar wall rate"));
        }
        Ok(Self::Linear { rate_per_s: rate })
    }
    fn validate(self) -> Result<(), WallExchangeError> {
        let valid = match self {
            Self::None => true,
            Self::Linear { rate_per_s } => nonnegative(rate_per_s),
            Self::Darcy { friction_factor, hydraulic_diameter_m } => {
                nonnegative(friction_factor) && positive(hydraulic_diameter_m)
            }
        };
        if valid { Ok(()) } else { Err(WallExchangeError::Invalid("wall drag inputs")) }
    }
}

/// Thermal storage coupled to one gas cell; all temperatures are in kelvin.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum WallThermal {
    /// No heat exchange.
    Adiabatic,
    /// Infinite prescribed-temperature reservoir; its heat is external.
    Reservoir {
        /// Positive prescribed temperature [K].
        temperature_k: f64,
        /// Complete nonnegative exposed h*A [W/K].
        conductance_w_per_k: f64,
    },
    /// Finite wall; its heat is internal to gas plus wall. The authored path
    /// uses constant capacity. An explicit material heat-storage binding uses
    /// its original energy integral and keeps the capacity below as a tangent.
    Lumped {
        /// Current positive wall temperature [K].
        temperature_k: f64,
        /// Positive complete wall heat capacity [J/K], or the refreshed local
        /// tangent when this contact has a retained material storage binding.
        heat_capacity_j_per_k: f64,
        /// Complete nonnegative exposed h*A [W/K].
        conductance_w_per_k: f64,
    },
    /// Finite storage on an original immutable thermal law. Its source binding
    /// must be supplied to the caloric exchange entry or the network operator.
    /// The coordinate is h or u exactly as declared by that law.
    Caloric {
        /// Fixed wall mass [kg].
        mass_kg: f64,
        /// Current original specific thermal-storage coordinate [J/kg].
        specific_energy_j_kg: f64,
        /// Complete nonnegative exposed h*A [W/K].
        conductance_w_per_k: f64,
    },
}

/// One cell's stationary radial momentum and thermal contact.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WallContact {
    /// Radial momentum-loss law.
    pub drag: WallDrag,
    /// Thermal boundary or finite storage.
    pub thermal: WallThermal,
}
impl WallContact {
    pub(super) fn validate(self) -> Result<(), WallExchangeError> {
        self.drag.validate()?;
        let valid = match self.thermal {
            WallThermal::Adiabatic => true,
            WallThermal::Reservoir { temperature_k, conductance_w_per_k } => {
                positive(temperature_k) && nonnegative(conductance_w_per_k)
            }
            WallThermal::Lumped {
                temperature_k, heat_capacity_j_per_k, conductance_w_per_k,
            } => positive(temperature_k)
                && positive(heat_capacity_j_per_k)
                && nonnegative(conductance_w_per_k),
            WallThermal::Caloric {
                mass_kg, specific_energy_j_kg, conductance_w_per_k,
            } => positive(mass_kg) && specific_energy_j_kg.is_finite()
                && nonnegative(conductance_w_per_k),
        };
        if valid { Ok(()) } else { Err(WallExchangeError::Invalid("wall thermal inputs")) }
    }
}

/// Local original-equation residual and scalar nonlinear work budgets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WallControl {
    /// Absolute energy and implicit heat-equation allowance [J].
    pub absolute_energy_tolerance_j: f64,
    /// Relative allowance against exchanged energy, never a large energy gauge.
    pub relative_energy_tolerance: f64,
    /// Maximum safeguarded temperature iterations, in 1..=128.
    pub maximum_iterations: usize,
}
impl Default for WallControl {
    fn default() -> Self {
        Self {
            absolute_energy_tolerance_j: 1e-10,
            relative_energy_tolerance: 1e-10,
            maximum_iterations: 96,
        }
    }
}
impl WallControl {
    pub(super) fn validate(self) -> Result<(), WallExchangeError> {
        if !positive(self.absolute_energy_tolerance_j)
            || !positive(self.relative_energy_tolerance)
            || self.relative_energy_tolerance > 1e-3
            || !(1..=128).contains(&self.maximum_iterations)
        {
            return Err(WallExchangeError::Invalid("wall solve controls"));
        }
        Ok(())
    }
    fn tolerance(self, scale: f64) -> f64 {
        self.absolute_energy_tolerance_j + self.relative_energy_tolerance * scale
    }
}

/// Refusal returns no mutated cell or thermal storage.
#[derive(Clone, Debug, PartialEq)]
pub enum WallExchangeError {
    /// Incompatible or unrepresentable physical input.
    Invalid(&'static str),
    /// Original caloric/source admission refused.
    Flow(FlowError),
    /// Original source-bound molecular or material transport law refused.
    Transport(super::dissipation::DissipationError),
    /// Original fixed-amount wall heat-storage source refused.
    #[cfg(feature = "material")]
    Caloric(fs_material::thermal::transport::ThermalTransportError),
    /// Original finite-wall thermal source refused its state or geometry.
    Thermal(fs_thermal::enthalpy::DepositionError),
    /// Actual cell state is outside the declared developed-wall approximation.
    TransportDomain {
        /// Physical quantity whose applicability interval was exceeded.
        quantity: &'static str,
        /// Actual dimensionless Reynolds or Mach number.
        value: f64,
        /// Inclusive lower applicability bound.
        minimum: f64,
        /// Inclusive upper applicability bound.
        maximum: f64,
    },
    /// Original heat equation did not meet the requested iteration budget.
    TemperatureIterations,
    /// Actual published gas/wall storage or heat flux exceeded its budget.
    Balance {
        /// Signed actual gas-plus-wall energy residual [J].
        energy_residual_j: f64,
        /// Signed backward-Euler heat-equation residual [J].
        thermal_residual_j: f64,
    },
    /// Cancellation before publication.
    Cancelled,
}
impl core::fmt::Display for WallExchangeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "gas wall exchange: {self:?}")
    }
}
impl std::error::Error for WallExchangeError {}
impl From<fs_thermal::enthalpy::DepositionError> for WallExchangeError {
    fn from(value: fs_thermal::enthalpy::DepositionError) -> Self {
        match value {
            fs_thermal::enthalpy::DepositionError::Cancelled => Self::Cancelled,
            other => Self::Thermal(other),
        }
    }
}
impl From<FlowError> for WallExchangeError {
    fn from(value: FlowError) -> Self {
        match value {
            FlowError::Cancelled => Self::Cancelled,
            other => Self::Flow(other),
        }
    }
}
impl From<super::dissipation::DissipationError> for WallExchangeError {
    fn from(value: super::dissipation::DissipationError) -> Self {
        match value {
            super::dissipation::DissipationError::Cancelled => Self::Cancelled,
            other => Self::Transport(other),
        }
    }
}

#[cfg(feature = "material")]
impl From<fs_material::thermal::transport::ThermalTransportError> for WallExchangeError {
    fn from(value: fs_material::thermal::transport::ThermalTransportError) -> Self {
        Self::Caloric(value)
    }
}

/// Actual source exchanges of one complete cell interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WallCellReport {
    /// Physical source duration, with no flow-clock mutation.
    pub elapsed: Time,
    /// Opposite impulse on the fixed support, in the local gas axis [N s].
    pub wall_impulse_ns: f64,
    /// Lost gas kinetic energy retained as sensible energy [J].
    pub drag_heating_j: f64,
    /// Actual heat entering gas [J].
    pub heat_to_gas_j: f64,
    /// Actual represented finite-wall thermal-storage change [J], in the
    /// original law's declared h/u coordinate for a caloric wall.
    pub wall_energy_change_j: f64,
    /// Heat entering gas from a prescribed reservoir only [J].
    pub external_heat_to_gas_j: f64,
    /// Actual gas-plus-wall thermal-storage change minus reservoir heat [J].
    /// An enthalpy wall retains its original pressure-work convention.
    pub energy_residual_j: f64,
    /// Actual implicit heat-equation residual [J].
    pub thermal_residual_j: f64,
    /// Safeguarded thermal iterations.
    pub iterations: usize,
}

/// Returned momentum/energy and contact; caller retains all original species.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WallUpdate {
    /// Updated gas momentum density [kg/(m² s)].
    pub momentum: f64,
    /// Updated species-referenced total gas energy density [J/m³].
    pub energy: f64,
    /// Updated finite wall, or unchanged prescribed/adiabatic contact.
    pub contact: WallContact,
    /// Actual physical balances.
    pub report: WallCellReport,
}

/// Apply exact frozen drag followed by an implicit original-caloric heat step.
///
/// Species and the caller's fields are borrowed and never changed. Gas density
/// and cell volume remain fixed. A finite wall is eliminated in the scalar
/// temperature solve, then actual represented gas and wall storage are audited.
/// Heat integration is first order; symmetry around transport does not promote
/// it to a global second-order or entropy-stability guarantee.
pub fn exchange_cell(
    gas: &MixtureGas,
    state: Cell<'_>,
    volume_m3: f64,
    duration: Time,
    contact: WallContact,
    control: WallControl,
    mut cancelled: impl FnMut() -> bool,
) -> Result<WallUpdate, WallExchangeError> {
    poll(&mut cancelled)?;
    control.validate()?;
    contact.validate()?;
    if !positive(volume_m3) || !nonnegative(duration.value()) {
        return Err(WallExchangeError::Invalid("wall volume or interval"));
    }
    let initial = gas.primitive(state, &mut cancelled)?;
    let dt = duration.value();
    let factor = match contact.drag {
        WallDrag::None => 1.0,
        WallDrag::Linear { rate_per_s } => {
            let exponent = rate_per_s * dt;
            if !nonnegative(exponent) {
                return Err(WallExchangeError::Invalid("wall drag exponent"));
            }
            det::exp(-exponent)
        }
        WallDrag::Darcy { friction_factor, hydraulic_diameter_m } => {
            let damping = 0.5 * friction_factor * dt
                * (initial.velocity.abs() / hydraulic_diameter_m);
            if !nonnegative(damping) {
                return Err(WallExchangeError::Invalid("wall Darcy increment"));
            }
            1.0 / (1.0 + damping)
        }
    };
    let momentum = state.momentum * factor;
    let old_kinetic = 0.5 * state.momentum * (state.momentum / initial.density);
    let kinetic = 0.5 * momentum * (momentum / initial.density);
    let mut energy = state.energy;
    let heated = Cell { partial_densities: state.partial_densities, momentum, energy };
    let primitive = gas.primitive(heated, &mut cancelled)?;
    let mut updated = contact;
    let mut report = WallCellReport {
        elapsed: duration,
        wall_impulse_ns: (state.momentum - momentum) * volume_m3,
        drag_heating_j: (old_kinetic - kinetic) * volume_m3,
        heat_to_gas_j: 0.0,
        wall_energy_change_j: 0.0,
        external_heat_to_gas_j: 0.0,
        energy_residual_j: 0.0,
        thermal_residual_j: 0.0,
        iterations: 0,
    };
    let thermal = match contact.thermal {
        WallThermal::Adiabatic => None,
        WallThermal::Reservoir { temperature_k, conductance_w_per_k } => {
            Some((temperature_k, conductance_w_per_k, None))
        }
        WallThermal::Lumped {
            temperature_k, heat_capacity_j_per_k, conductance_w_per_k,
        } => Some((temperature_k, conductance_w_per_k, Some(heat_capacity_j_per_k))),
        WallThermal::Caloric { .. } => {
            return Err(WallExchangeError::Invalid("caloric wall requires its original thermal law"));
        }
    };
    if let Some((wall_t, conductance, capacity)) = thermal {
        let integrated = conductance * dt;
        if !nonnegative(integrated) {
            return Err(WallExchangeError::Invalid("integrated wall conductance"));
        }
        if integrated > 0.0 && wall_t != primitive.temperature {
            let effective = capacity.map_or(integrated, |c| {
                if integrated <= c {
                    integrated / (1.0 + integrated / c)
                } else {
                    c / (1.0 + c / integrated)
                }
            });
            let target = state.energy - kinetic;
            let domain = gas.domain();
            let lo = primitive.temperature.min(wall_t).max(domain.minimum_temperature);
            let hi = primitive.temperature.max(wall_t).min(domain.maximum_temperature);
            let equation = |t: f64, _: &mut dyn FnMut() -> bool|
                -> Result<(f64, Option<f64>, f64), WallExchangeError> {
                let (thermal, cv) = gas.caloric(state.partial_densities, t)?;
                let change = (thermal - target) * volume_m3;
                let heat = effective * (wall_t - t);
                let slope = cv * volume_m3 + effective;
                let residual = change - heat;
                if !residual.is_finite() || !positive(slope) || !heat.is_finite() {
                    return Err(WallExchangeError::Invalid("wall thermal equation"));
                }
                Ok((residual, Some(slope), change.abs().max(heat.abs())))
            };
            let (temperature, iterations) = solve_scalar(
                lo, hi, primitive.temperature, control, equation, &mut cancelled,
            )?;
            report.iterations = iterations;
            let (thermal, _) = gas.caloric(state.partial_densities, temperature)?;
            energy = state.energy + (thermal - target);
            let actual = gas.primitive(
                Cell { partial_densities: state.partial_densities, momentum, energy },
                &mut cancelled,
            )?;
            report.heat_to_gas_j = (energy - state.energy) * volume_m3;
            let actual_wall_t = if let Some(capacity) = capacity {
                let fraction = effective / capacity;
                let next = if fraction <= 0.5 {
                    wall_t + fraction * (actual.temperature - wall_t)
                } else {
                    actual.temperature + (1.0 - fraction) * (wall_t - actual.temperature)
                };
                report.wall_energy_change_j = capacity * (next - wall_t);
                updated.thermal = WallThermal::Lumped {
                    temperature_k: next,
                    heat_capacity_j_per_k: capacity,
                    conductance_w_per_k: conductance,
                };
                next
            } else {
                report.external_heat_to_gas_j = report.heat_to_gas_j;
                wall_t
            };
            report.thermal_residual_j =
                report.heat_to_gas_j - integrated * (actual_wall_t - actual.temperature);
            report.energy_residual_j = report.heat_to_gas_j + report.wall_energy_change_j
                - report.external_heat_to_gas_j;
            let scale = report.heat_to_gas_j.abs().max(report.wall_energy_change_j.abs());
            let tolerance = control.tolerance(scale);
            if report.energy_residual_j.abs() > tolerance
                || report.thermal_residual_j.abs() > tolerance
            {
                return Err(WallExchangeError::Balance {
                    energy_residual_j: report.energy_residual_j,
                    thermal_residual_j: report.thermal_residual_j,
                });
            }
        }
    }
    updated.validate()?;
    if ![
        momentum, energy, report.wall_impulse_ns, report.drag_heating_j,
        report.heat_to_gas_j, report.wall_energy_change_j,
        report.external_heat_to_gas_j, report.energy_residual_j,
        report.thermal_residual_j,
    ].iter().all(|x| x.is_finite()) {
        return Err(WallExchangeError::Invalid("wall extensive source balance"));
    }
    poll(&mut cancelled)?;
    Ok(WallUpdate { momentum, energy, contact: updated, report })
}

/// One safeguarded monotone equation owner for authored and original-source
/// storage. A source without an exact tangent uses the same bounded bisection.
fn solve_scalar(
    mut lo: f64,
    mut hi: f64,
    initial: f64,
    control: WallControl,
    mut equation: impl FnMut(f64, &mut dyn FnMut() -> bool)
        -> Result<(f64, Option<f64>, f64), WallExchangeError>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(f64, usize), WallExchangeError> {
    let lower = equation(lo, cancelled)?;
    let upper = equation(hi, cancelled)?;
    if lower.0 > control.tolerance(lower.2) || upper.0 < -control.tolerance(upper.2) {
        return Err(FlowError::EnergyDomain.into());
    }
    let mut value = initial.clamp(lo, hi);
    for iteration in 1..=control.maximum_iterations {
        poll(cancelled)?;
        let (residual, slope, scale) = equation(value, cancelled)?;
        if residual.abs() <= control.tolerance(scale) {
            return Ok((value, iteration));
        }
        if residual < 0.0 { lo = value; } else { hi = value; }
        let newton = slope.map(|slope| value - residual / slope);
        let next = match newton {
            Some(next) if next.is_finite() && next > lo && next < hi => next,
            _ => 0.5 * lo + 0.5 * hi,
        };
        if next == value { break; }
        value = next;
    }
    Err(WallExchangeError::TemperatureIterations)
}

fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), WallExchangeError> {
    if cancelled() { Err(WallExchangeError::Cancelled) } else { Ok(()) }
}
