//! One fixed gas amount exchanging frictional momentum with a finite support.
//!
//! The support slides along the gas axis without changing geometry. No external
//! drive maintains its velocity. Opposite impulse and reduced-mass slip decay
//! couple both momenta over one interval; an authored heat fraction partitions
//! the lost relative kinetic energy before the original caloric heat solve.

use super::{
    Cell, MaterialWallCaloric, MixtureCaloric, MixtureGas, Time, WallCaloricConstraint,
    WallContact, WallControl, WallDrag, WallExchangeError, WallGasCaloric, WallThermal, WallUpdate,
    nonnegative, poll, positive,
};

pub(super) mod mechanics;
mod shared;
pub use shared::{
    SharedRecoilCell, SharedRecoilCellUpdate, SharedRecoilControl, SharedRecoilHeatCell,
    SharedRecoilHeatUpdate, SharedRecoilHeatWorkspace, SharedRecoilUpdate, SharedRecoilWorkspace,
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
    /// remainder heats a finite support through its authored constant capacity
    /// or explicitly supplied original material storage. No fraction is inferred.
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
        None,
        control,
        cancelled,
    )
}

/// The same exact recoil with the support's original temperature-dependent or
/// latent heat storage, followed by the existing implicit material heat owner.
///
/// The wall share of slip dissipation enters the original storage integral;
/// no frozen tangent is substituted for a finite energy change. The returned
/// wall state retains its temperature or native energy coordinate, including
/// latent inventory. Its complete endpoint storage change, the gas energy and
/// the support kinetic energy must close the original conservation budget.
///
/// The source must declare internal energy or the existing negligible-solid-
/// expansion approximation. Unqualified enthalpy or sensible heat alone does
/// not close this mechanical energy balance. The fixed thermal amount belongs
/// to the authored support inertia; it is not added to that inertia again.
/// Conductance, friction and geometry remain frozen over this split interval;
/// the heat solve remains backward Euler. A refusal publishes neither endpoint.
#[allow(clippy::too_many_arguments)]
pub fn exchange_cell_with_recoil_and_heat_storage(
    gas: &MixtureGas,
    state: Cell<'_>,
    volume_m3: f64,
    duration: Time,
    support: SlidingSupport,
    law: RecoilLaw,
    thermal: WallThermal,
    wall: &MaterialWallCaloric<'_>,
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
        Some(wall),
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
    wall: Option<&MaterialWallCaloric<'_>>,
    control: RecoilControl,
    mut cancelled: impl FnMut() -> bool,
) -> Result<RecoilUpdate, WallExchangeError> {
    let initial_contact = prepare_contact(
        thermal,
        wall,
        law.gas_heat_fraction,
        control,
        &mut cancelled,
    )?;
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
    complete_motion(
        gas,
        initial_energy,
        volume_m3,
        duration,
        support,
        law.gas_heat_fraction,
        initial_contact,
        wall,
        control,
        motion,
        &mut cancelled,
    )
}

/// One original wall-storage admission for frozen and live-source recoil.
pub(super) fn prepare_contact(
    thermal: WallThermal,
    wall: Option<&MaterialWallCaloric<'_>>,
    gas_heat_fraction: f64,
    control: RecoilControl,
    mut cancelled: impl FnMut() -> bool,
) -> Result<WallContact, WallExchangeError> {
    poll(&mut cancelled)?;
    control.wall.validate()?;
    let mut initial_contact = WallContact {
        drag: WallDrag::None,
        thermal,
    };
    initial_contact.validate()?;
    if let Some(wall) = wall {
        if !matches!(
            wall.constraint(),
            WallCaloricConstraint::IsochoricInternalEnergy
                | WallCaloricConstraint::StationarySolidNegligibleExpansion
        ) {
            return Err(WallExchangeError::Invalid(
                "finite-support recoil requires internal energy or explicit negligible expansion",
            ));
        }
        initial_contact = wall.bind_contact(initial_contact)?;
    }
    if !positive(control.momentum_tolerance_ns)
        || !gas_heat_fraction.is_finite()
        || !(0.0..=1.0).contains(&gas_heat_fraction)
        || (gas_heat_fraction < 1.0
            && wall.is_none()
            && !matches!(thermal, WallThermal::Lumped { .. }))
    {
        return Err(WallExchangeError::Invalid(
            "recoil momentum budget or finite-wall heat partition",
        ));
    }
    Ok(initial_contact)
}

/// Deposit the accepted physical slip loss through the existing original
/// caloric owner and audit actual complete gas, support and wall endpoints.
#[allow(clippy::too_many_arguments)]
pub(super) fn complete_motion(
    gas: &impl WallGasCaloric,
    initial_energy: f64,
    volume_m3: f64,
    duration: Time,
    support: SlidingSupport,
    gas_heat_fraction: f64,
    initial_contact: WallContact,
    wall: Option<&MaterialWallCaloric<'_>>,
    control: RecoilControl,
    motion: mechanics::Motion,
    mut cancelled: impl FnMut() -> bool,
) -> Result<RecoilUpdate, WallExchangeError> {
    if motion.momentum_residual_ns.abs() > control.momentum_tolerance_ns {
        return Err(WallExchangeError::Invalid(
            "finite-support momentum residual exceeds budget",
        ));
    }
    let gas_heating = gas_heat_fraction * motion.dissipation_j;
    let wall_heating = motion.dissipation_j - gas_heating;
    let mut contact = initial_contact;
    if let Some(wall) = wall {
        contact = wall.after(contact, wall_heating)?.contact;
    } else if let WallThermal::Lumped {
        temperature_k,
        heat_capacity_j_per_k,
        conductance_w_per_k,
    } = initial_contact.thermal
    {
        let next = temperature_k + wall_heating / heat_capacity_j_per_k;
        contact.thermal = WallThermal::Lumped {
            temperature_k: next,
            heat_capacity_j_per_k,
            conductance_w_per_k,
        };
    }
    let mechanical_energy = initial_energy + (motion.work_to_gas_j - wall_heating) / volume_m3;
    let mut update = match wall {
        Some(wall) => super::exchange_caloric_with_heat_storage(
            gas,
            motion.momentum,
            mechanical_energy,
            volume_m3,
            duration,
            contact,
            wall,
            None,
            control.wall,
            &mut cancelled,
        )?,
        None => super::exchange_caloric(
            gas,
            motion.momentum,
            mechanical_energy,
            volume_m3,
            duration,
            contact,
            control.wall,
            &mut cancelled,
        )?,
    };
    update.report.wall_impulse_ns = -motion.gas_impulse_ns;
    update.report.mechanical_work_to_gas_j = motion.work_to_gas_j;
    // Preserve the original sensible-energy account: total slip loss heats
    // gas, then the authored wall share leaves it as heat. Explicit partition
    // fields below expose the net shares without counting that transfer twice.
    update.report.drag_heating_j = motion.dissipation_j;
    update.report.heat_to_gas_j -= wall_heating;
    // Audit the complete returned storage endpoints, not a sum of separately
    // rounded frictional-heating and thermal-exchange differences.
    if let Some(wall) = wall {
        update.report.wall_energy_change_j =
            wall.energy_change_j(initial_contact, update.contact)?;
    } else if let (
        WallThermal::Lumped {
            temperature_k: initial,
            heat_capacity_j_per_k,
            ..
        },
        WallThermal::Lumped {
            temperature_k: final_temperature,
            ..
        },
    ) = (initial_contact.thermal, update.contact.thermal)
    {
        update.report.wall_energy_change_j = heat_capacity_j_per_k * (final_temperature - initial);
    }
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
