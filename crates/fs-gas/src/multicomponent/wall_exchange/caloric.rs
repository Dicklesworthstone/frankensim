//! Original material storage in the existing gas/wall scalar exchange.

use super::*;
use fs_material::thermal::transport::storage::PreparedHeatStorage;
use std::sync::Arc;

/// Required physical interpretation of the source capacity. No Cp-to-Cv
/// conversion or unrepresented pressure/deformation work is inferred.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WallCaloricConstraint {
    /// The supplied capacity is dU/dT under the fixed wall-volume constraint.
    IsochoricInternalEnergy,
    /// Stationary solid approximation: thermal expansion, pressure-volume work
    /// and the Cp-Cv difference are explicitly neglected over this source span.
    StationarySolidNegligibleExpansion,
}

/// Immutable source and fixed constituent amounts, shared without copying tables.
/// Temperature remains exclusively in the original `WallContact` transaction.
/// A bound contact's capacity field is a refreshed tangent cache; its storage
/// change is ALWAYS the source integral, never that cache times delta T.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialWallCaloric {
    storage: Arc<PreparedHeatStorage>,
    constraint: WallCaloricConstraint,
}

impl MaterialWallCaloric {
    /// Bind an admitted source with an explicit compatible work constraint.
    /// This retains fixed amounts and source context; it adds no phase change,
    /// wall deformation, reaction, spatial wall conduction or external work.
    pub fn new(storage: Arc<PreparedHeatStorage>, constraint: WallCaloricConstraint) -> Self {
        Self {
            storage,
            constraint,
        }
    }
    pub fn storage(&self) -> &PreparedHeatStorage {
        &self.storage
    }
    pub const fn constraint(&self) -> WallCaloricConstraint {
        self.constraint
    }

    /// Construct the original finite-wall state from material data alone.
    pub fn contact(
        &self,
        temperature_k: f64,
        conductance_w_per_k: f64,
        drag: WallDrag,
    ) -> Result<WallContact, WallExchangeError> {
        let contact = WallContact {
            drag,
            thermal: WallThermal::Lumped {
                temperature_k,
                conductance_w_per_k,
                heat_capacity_j_per_k: self.storage.heat_capacity_j_k(temperature_k)?,
            },
        };
        contact.validate()?;
        Ok(contact)
    }

    pub(crate) fn bind_contact(
        &self,
        contact: WallContact,
    ) -> Result<WallContact, WallExchangeError> {
        let WallThermal::Lumped {
            temperature_k,
            conductance_w_per_k,
            ..
        } = contact.thermal
        else {
            return Err(WallExchangeError::Invalid(
                "material storage requires a finite wall",
            ));
        };
        self.contact(temperature_k, conductance_w_per_k, contact.drag)
    }
}

/// Source-backed nonlinear wall storage with the SAME gas calorics, radial
/// drag and safeguarded scalar owner as `exchange_cell`.
///
/// Solve Q = dt G (Tw(Q)-Tg(Q)), where each original energy owner supplies its
/// own temperature and tangent. The bracket intersects both source domains and
/// the passive temperature interval. Actual represented energy changes and
/// endpoint heat flux must pass the original budget before anything returns.
/// Heat integration is backward Euler; coefficients and geometry are held for
/// this source interval. The gas retains drag heating and all original species.
#[allow(clippy::too_many_arguments)]
pub fn exchange_cell_with_heat_storage(
    gas: &MixtureGas,
    state: Cell<'_>,
    volume_m3: f64,
    duration: Time,
    contact: WallContact,
    wall: &MaterialWallCaloric,
    control: WallControl,
    mut cancelled: impl FnMut() -> bool,
) -> Result<WallUpdate, WallExchangeError> {
    poll(&mut cancelled)?;
    let contact = wall.bind_contact(contact)?;
    let WallThermal::Lumped {
        temperature_k: initial_wall,
        conductance_w_per_k,
        ..
    } = contact.thermal
    else {
        unreachable!("binding admitted finite storage");
    };
    // Reuse the original radial loss and gas sensible-heating owner exactly.
    let mut update = exchange_cell(
        gas,
        state,
        volume_m3,
        duration,
        WallContact {
            drag: contact.drag,
            thermal: WallThermal::Adiabatic,
        },
        control,
        &mut cancelled,
    )?;
    update.contact = contact;
    let gas_at = |energy, cancelled: &mut dyn FnMut() -> bool| {
        gas.primitive(
            Cell {
                partial_densities: state.partial_densities,
                momentum: update.momentum,
                energy,
            },
            cancelled,
        )
    };
    let initial_gas = gas_at(update.energy, &mut cancelled)?;
    let integrated = conductance_w_per_k * duration.value();
    if !nonnegative(integrated) {
        return Err(WallExchangeError::Invalid(
            "integrated material wall conductance",
        ));
    }
    if integrated == 0.0 || initial_wall == initial_gas.temperature {
        poll(&mut cancelled)?;
        return Ok(update);
    }
    let initial_flux = integrated * (initial_wall - initial_gas.temperature);
    if !initial_flux.is_finite() {
        return Err(WallExchangeError::Invalid("material wall heat increment"));
    }
    if initial_flux.abs() <= control.tolerance(initial_flux.abs()) {
        update.report.thermal_residual_j = -initial_flux;
        poll(&mut cancelled)?;
        return Ok(update);
    }
    let storage = wall.storage();
    let [wall_low, wall_high] = storage.temperature_span_k();
    let domain = gas.domain();
    let kinetic = 0.5 * update.momentum * (update.momentum / initial_gas.density);
    let target = update.energy - kinetic;
    let gas_change_at = |temperature| -> Result<f64, WallExchangeError> {
        let energy = gas.caloric(state.partial_densities, temperature)?.0;
        let change = (energy - target) * volume_m3;
        if !change.is_finite() {
            return Err(WallExchangeError::Invalid("material wall gas energy range"));
        }
        Ok(change)
    };
    // Neither endpoint asks a source for an unsupported trial temperature.
    let (low, high) = if initial_wall > initial_gas.temperature {
        (
            0.0,
            gas_change_at(initial_wall.min(domain.maximum_temperature))?.min(
                -storage.energy_change_j(initial_wall, initial_gas.temperature.max(wall_low))?,
            ),
        )
    } else {
        (
            gas_change_at(initial_wall.max(domain.minimum_temperature))?.max(
                -storage.energy_change_j(initial_wall, initial_gas.temperature.min(wall_high))?,
            ),
            0.0,
        )
    };
    if !low.is_finite() || !high.is_finite() || low >= high {
        return Err(WallExchangeError::Invalid(
            "material wall has no supported heat interval",
        ));
    }
    let equation = |heat: f64, cancelled: &mut dyn FnMut() -> bool| {
        let gas_state = gas_at(update.energy + heat / volume_m3, cancelled)?;
        let wall_t = storage.temperature_after_energy_change(initial_wall, -heat)?;
        let wall_capacity = storage.heat_capacity_j_k(wall_t)?;
        let gas_capacity = gas
            .caloric(state.partial_densities, gas_state.temperature)?
            .1
            * volume_m3;
        let flux = integrated * (wall_t - gas_state.temperature);
        let residual = heat - flux;
        let slope = 1.0 + integrated * (1.0 / wall_capacity + 1.0 / gas_capacity);
        if !positive(gas_capacity) || !positive(slope) || !residual.is_finite() || !flux.is_finite()
        {
            return Err(WallExchangeError::Invalid("material wall heat equation"));
        }
        Ok((residual, slope, heat.abs().max(flux.abs())))
    };
    let (heat, iterations) = solve_scalar(low, high, 0.0, control, equation, &mut cancelled)?;
    let energy = update.energy + heat / volume_m3;
    let actual_heat = (energy - update.energy) * volume_m3;
    if !actual_heat.is_finite() || (heat != 0.0 && actual_heat == 0.0) {
        return Err(WallExchangeError::Invalid(
            "material wall heat is unrepresentable",
        ));
    }
    let actual_gas = gas_at(energy, &mut cancelled)?;
    let actual_wall = storage.temperature_after_energy_change(initial_wall, -actual_heat)?;
    let wall_change = storage.energy_change_j(initial_wall, actual_wall)?;
    let energy_residual = actual_heat + wall_change;
    let thermal_residual = actual_heat - integrated * (actual_wall - actual_gas.temperature);
    let tolerance = control.tolerance(actual_heat.abs().max(wall_change.abs()));
    if !energy_residual.is_finite()
        || !thermal_residual.is_finite()
        || energy_residual.abs() > tolerance
        || thermal_residual.abs() > tolerance
    {
        return Err(WallExchangeError::Balance {
            energy_residual_j: energy_residual,
            thermal_residual_j: thermal_residual,
        });
    }
    update.contact = wall.contact(actual_wall, conductance_w_per_k, contact.drag)?;
    update.energy = energy;
    update.report.heat_to_gas_j = actual_heat;
    update.report.wall_energy_change_j = wall_change;
    update.report.energy_residual_j = energy_residual;
    update.report.thermal_residual_j = thermal_residual;
    update.report.iterations = iterations;
    poll(&mut cancelled)?;
    Ok(update)
}
