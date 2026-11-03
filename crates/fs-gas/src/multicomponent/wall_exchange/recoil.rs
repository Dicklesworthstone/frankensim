//! One fixed gas amount exchanging frictional momentum with a finite support.
//!
//! The support slides along the gas axis without changing geometry. No external
//! drive maintains its velocity. Opposite impulse and reduced-mass slip decay
//! couple both momenta over one interval; an authored heat fraction partitions
//! the lost relative kinetic energy before the original caloric heat solve.

use super::{
    Cell, MixtureCaloric, MixtureGas, Time, WallContact, WallControl, WallDrag, WallExchangeError,
    WallGasCaloric, WallThermal, WallUpdate, nonnegative, poll, positive,
};

mod mechanics;
mod shared;
pub use shared::{
    SharedRecoilCell, SharedRecoilCellUpdate, SharedRecoilControl, SharedRecoilUpdate,
    SharedRecoilWorkspace,
};

/// Finite translating support, with no prescribed drive or moving gas volume.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SlidingSupport {
    /// Positive finite inertial mass [kg].
    pub mass_kg: f64,
    /// Axial laboratory velocity [m/s].
    pub velocity_m_s: f64,
}

/// Frozen resistance defined by gas acceleration against a fixed support.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RecoilLaw {
    /// Nonnegative linear coefficient a in du/dt = -a*s - b*abs(s)*s [1/s].
    pub linear_rate_per_s: f64,
    /// Nonnegative quadratic coefficient b [1/m].
    pub quadratic_per_m: f64,
    /// Authored fraction of slip dissipation retained by gas, in [0,1]. The
    /// remainder heats a finite `WallThermal::Lumped` support. No fraction is inferred.
    pub gas_heat_fraction: f64,
}

/// Independent dimensional acceptance budgets for the coupled interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RecoilControl {
    /// Original caloric solve and complete gas/support energy budget [J].
    pub wall: WallControl,
    /// Positive absolute gas-plus-support momentum residual budget [N s].
    pub momentum_tolerance_ns: f64,
}

/// Both accepted endpoints; the caller retains the original species densities.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RecoilUpdate {
    /// Gas state, returned thermal contact and complete heat/work exchanges.
    pub gas: WallUpdate,
    /// Support after the same source interval.
    pub support: SlidingSupport,
    /// Actual returned support kinetic-energy change [J].
    pub support_kinetic_energy_change_j: f64,
    /// Nonnegative relative-slip energy lost, before its authored partition [J].
    pub slip_dissipation_j: f64,
    /// Slip heat assigned to gas before subsequent thermal contact [J].
    pub gas_slip_heating_j: f64,
    /// Slip heat assigned to finite wall storage before thermal contact [J].
    pub wall_slip_heating_j: f64,
    /// Actual gas impulse plus actual support momentum change [N s].
    pub momentum_residual_ns: f64,
    /// Actual gas energy + support kinetic + finite-wall heat change, minus
    /// prescribed-reservoir heat [J]. No residual is repaired with an energy source.
    pub energy_residual_j: f64,
}

/// Exact frozen friction recoil followed by the original implicit heat owner.
///
/// Gas volume, species and mass remain fixed. A heat fraction below one requires
/// finite constant-capacity wall storage; original material caloric bindings
/// are not inferred. With all slip heat in gas, the usual adiabatic/reservoir/
/// lumped thermal choices remain available. Inputs are borrowed or copied and
/// no flow clock changes: a refusal publishes neither endpoint. One support
/// belongs to this one gas amount; this is not a sequential shared-wall network solve.
#[allow(clippy::too_many_arguments)]
pub fn exchange_cell_with_recoil(
    gas: &MixtureGas,
    state: Cell<'_>,
    volume_m3: f64,
    duration: Time,
    support: SlidingSupport,
    law: RecoilLaw,
    thermal: WallThermal,
    control: RecoilControl,
    cancelled: impl FnMut() -> bool,
) -> Result<RecoilUpdate, WallExchangeError> {
    exchange_caloric_with_recoil(
        &MixtureCaloric {
            gas,
            partial_densities: state.partial_densities,
        },
        state.momentum,
        state.energy,
        volume_m3,
        duration,
        support,
        law,
        thermal,
        control,
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
fn exchange_caloric_with_recoil(
    gas: &impl WallGasCaloric,
    initial_momentum: f64,
    initial_energy: f64,
    volume_m3: f64,
    duration: Time,
    support: SlidingSupport,
    law: RecoilLaw,
    thermal: WallThermal,
    control: RecoilControl,
    mut cancelled: impl FnMut() -> bool,
) -> Result<RecoilUpdate, WallExchangeError> {
    poll(&mut cancelled)?;
    control.wall.validate()?;
    let mut contact = WallContact {
        drag: WallDrag::None,
        thermal,
    };
    contact.validate()?;
    if !positive(control.momentum_tolerance_ns)
        || !law.gas_heat_fraction.is_finite()
        || !(0.0..=1.0).contains(&law.gas_heat_fraction)
        || (law.gas_heat_fraction < 1.0 && !matches!(thermal, WallThermal::Lumped { .. }))
    {
        return Err(WallExchangeError::Invalid(
            "recoil momentum budget or finite-wall heat partition",
        ));
    }
    let initial = gas.primitive(initial_momentum, initial_energy, &mut cancelled)?;
    let motion = mechanics::advance(
        initial_momentum,
        initial.density,
        initial.velocity,
        volume_m3,
        duration.value(),
        support,
        law,
    )?;
    if motion.momentum_residual_ns.abs() > control.momentum_tolerance_ns {
        return Err(WallExchangeError::Invalid(
            "finite-support momentum residual exceeds budget",
        ));
    }
    let gas_heating = law.gas_heat_fraction * motion.dissipation_j;
    let wall_heating = motion.dissipation_j - gas_heating;
    let mut deposited_wall_heat = 0.0;
    if let WallThermal::Lumped {
        temperature_k,
        heat_capacity_j_per_k,
        conductance_w_per_k,
    } = thermal
    {
        let next = temperature_k + wall_heating / heat_capacity_j_per_k;
        deposited_wall_heat = heat_capacity_j_per_k * (next - temperature_k);
        contact.thermal = WallThermal::Lumped {
            temperature_k: next,
            heat_capacity_j_per_k,
            conductance_w_per_k,
        };
    }
    let mechanical_energy = initial_energy + (motion.work_to_gas_j - wall_heating) / volume_m3;
    let mut update = super::exchange_caloric(
        gas,
        motion.momentum,
        mechanical_energy,
        volume_m3,
        duration,
        contact,
        control.wall,
        &mut cancelled,
    )?;
    update.report.wall_impulse_ns = -motion.gas_impulse_ns;
    update.report.mechanical_work_to_gas_j = motion.work_to_gas_j;
    // Preserve the original sensible-energy account: total slip loss heats
    // gas, then the authored wall share leaves it as heat. Explicit partition
    // fields below expose the net shares without counting that transfer twice.
    update.report.drag_heating_j = motion.dissipation_j;
    update.report.heat_to_gas_j -= wall_heating;
    update.report.wall_energy_change_j += deposited_wall_heat;
    let storage_change = (update.energy - initial_energy) * volume_m3
        + update.report.wall_energy_change_j
        - update.report.external_heat_to_gas_j;
    update.report.energy_residual_j = storage_change - motion.work_to_gas_j;
    let energy_residual_j = storage_change + motion.support_kinetic_energy_change_j;
    if ![
        update.report.heat_to_gas_j,
        update.report.wall_energy_change_j,
        update.report.energy_residual_j,
        energy_residual_j,
    ]
    .into_iter()
    .all(f64::is_finite)
    {
        return Err(WallExchangeError::Invalid(
            "finite-support extensive energy balance",
        ));
    }
    let scale = motion
        .work_to_gas_j
        .abs()
        .max(motion.dissipation_j)
        .max(update.report.heat_to_gas_j.abs())
        .max(update.report.wall_energy_change_j.abs());
    if energy_residual_j.abs() > control.wall.tolerance(scale)
        || update.report.energy_residual_j.abs() > control.wall.tolerance(scale)
    {
        return Err(WallExchangeError::Balance {
            energy_residual_j,
            thermal_residual_j: update.report.thermal_residual_j,
        });
    }
    poll(&mut cancelled)?;
    Ok(RecoilUpdate {
        gas: update,
        support: SlidingSupport {
            velocity_m_s: motion.support_velocity_m_s,
            ..support
        },
        support_kinetic_energy_change_j: motion.support_kinetic_energy_change_j,
        slip_dissipation_j: motion.dissipation_j,
        gas_slip_heating_j: gas_heating,
        wall_slip_heating_j: wall_heating,
        momentum_residual_ns: motion.momentum_residual_ns,
        energy_residual_j,
    })
}

#[cfg(test)]
mod mixture_tests;
#[cfg(test)]
mod tests;
