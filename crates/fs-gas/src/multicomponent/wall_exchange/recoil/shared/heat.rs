//! One shared mechanical support and one original finite thermal inventory.

use super::{
    RecoilControl, RecoilLaw, SharedRecoilCell, SharedRecoilControl, SharedRecoilUpdate,
    SharedRecoilWorkspace, SlidingSupport, Sum, Time, WallExchangeError, allocated, nonnegative,
    poll,
};
use crate::multicomponent::wall_exchange::{
    MaterialWallCaloric, SharedWallControl, SharedWallReport, WallCaloricConstraint, WallContact,
    WallDrag, WallThermal,
    shared::{Group, SharedWall, SharedWallGas, SharedWallPort, SharedWallStorage, Work, mixture},
};
use crate::multicomponent::{Cell, MixtureGas};
use std::cell::Cell as Counter;

/// One original gas amount coupled to the common support and common wall.
/// Species remain borrowed and fixed during this source interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SharedRecoilHeatCell<'a> {
    /// Original conservative gas state, including its own composition.
    pub state: Cell<'a>,
    /// Positive fixed gas volume [m³].
    pub volume_m3: f64,
    /// Frozen resistance and authored gas share of this cell's slip loss.
    pub law: RecoilLaw,
    /// Complete nonnegative conductance to the common thermal inventory [W/K].
    pub conductance_w_per_k: f64,
}

/// Complete candidate endpoints in the original input order. No caller state
/// is changed; retaining this result borrows the workspace until publication.
#[derive(Debug)]
pub struct SharedRecoilHeatUpdate<'a> {
    /// Original common-support motion, impulses, work and unpartitioned slip
    /// loss. Its energy residual covers mechanics alone.
    pub mechanics: SharedRecoilUpdate<'a>,
    /// Returned total gas energy densities [J/m³], paired with the momenta in
    /// `mechanics.cells`. Every original species density is retained unchanged.
    pub energies: &'a [f64],
    /// Heat entering each gas during the subsequent simultaneous thermal solve
    /// [J]. This excludes the frictional wall share deposited before that solve.
    pub heat_to_gas_j: &'a [f64],
    /// Each cell's slip loss assigned to the ONE wall inventory [J]. Its gas
    /// share is the corresponding mechanical slip loss minus this amount.
    pub wall_slip_heating_j: &'a [f64],
    /// One returned material state, including native latent energy if present.
    pub thermal: WallThermal,
    /// Original source's admitted mechanical storage interpretation.
    pub storage_constraint: WallCaloricConstraint,
    /// Actual complete returned gas energy change, summed over all cells [J].
    pub gas_energy_change_j: f64,
    /// Actual original-to-returned wall storage difference, counted once [J].
    pub wall_energy_change_j: f64,
    /// Actual gas + support kinetic + wall storage change [J]. No external
    /// source or numerical-residual correction is introduced.
    pub energy_residual_j: f64,
    /// Simultaneous thermal-stage balances and original source work. Source
    /// evaluations also include admission, frictional storage and final audits.
    pub thermal_report: SharedWallReport,
}

/// Prepared O(cells) scratch for shared recoil with original nonlinear or
/// latent wall storage. Mechanics and thermal contact reuse their existing
/// owners; there is no second caloric curve or per-cell copy of wall inventory.
pub struct SharedRecoilHeatWorkspace<'law> {
    mechanics: SharedRecoilWorkspace,
    inputs: Vec<SharedRecoilCell>,
    wall_heat: Vec<f64>,
    thermal: Group<'law>,
}

impl<'law> SharedRecoilHeatWorkspace<'law> {
    /// Bind one original finite storage and allocate all scratch once. Its own
    /// conductance must be zero: the individual interfaces supply conductance.
    /// `thermal` admits the initial storage representation, including its mass;
    /// each advance still receives the caller's actual current wall state.
    ///
    /// Only internal energy or explicit negligible solid expansion closes this
    /// mechanical balance. The thermal amount belongs to the supplied support
    /// inertia; this workspace never adds another copy of its mass.
    pub fn new(
        cell_count: usize,
        wall: MaterialWallCaloric<'law>,
        thermal: WallThermal,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, WallExchangeError> {
        poll(&mut cancelled)?;
        if cell_count == 0
            || !matches!(
                wall.constraint(),
                WallCaloricConstraint::IsochoricInternalEnergy
                    | WallCaloricConstraint::StationarySolidNegligibleExpansion
            )
        {
            return Err(WallExchangeError::Invalid(
                "shared recoil needs cells and internal energy or explicit negligible expansion",
            ));
        }
        let contact = bind_wall(&wall, thermal)?;
        let mut ports = allocated(
            cell_count,
            SharedWallPort {
                gas: SharedWallGas::Plenum(0),
                conductance_w_per_k: 0.0,
            },
            &mut cancelled,
        )?;
        for (i, port) in ports.iter_mut().enumerate() {
            // These are private local indices for the original shared thermal
            // owner, not additional plenums or network state.
            port.gas = SharedWallGas::Plenum(i);
        }
        let workspace = Self {
            mechanics: SharedRecoilWorkspace::new(cell_count, &mut cancelled)?,
            inputs: allocated(
                cell_count,
                SharedRecoilCell {
                    density_kg_m3: 0.0,
                    momentum_kg_m2_s: 0.0,
                    volume_m3: 0.0,
                    linear_rate_per_s: 0.0,
                    quadratic_per_m: 0.0,
                },
                &mut cancelled,
            )?,
            wall_heat: allocated(cell_count, 0.0, &mut cancelled)?,
            thermal: Group::new(
                SharedWall {
                    storage: SharedWallStorage {
                        thermal: contact.thermal,
                        heat_storage: Some(wall),
                    },
                    ports,
                },
                allocated(cell_count, 0.0, &mut cancelled)?,
            )?,
        };
        poll(&mut cancelled)?;
        Ok(workspace)
    }

    /// Advance common-support mechanics, deposit the sum of authored wall
    /// fractions ONCE, then solve all gas/wall heat interfaces simultaneously.
    /// The original gas Cv, material capacity integral and latent inventory
    /// remain active. Whole returned storage endpoints close the energy audit.
    ///
    /// Resistance, conductance and each gas volume are frozen over this source
    /// interval. The mechanical solver retains its own bounded time estimates;
    /// thermal splitting and the shared backward-Euler solve remain first order.
    /// No geometry, flow clock, or caller-owned state changes. Refusal or
    /// cancellation changes only scratch; retry uses the supplied inputs anew.
    #[allow(clippy::too_many_arguments)]
    pub fn advance<'a>(
        &'a mut self,
        gas: &MixtureGas,
        cells: &[SharedRecoilHeatCell<'_>],
        duration: Time,
        support: SlidingSupport,
        thermal: WallThermal,
        control: SharedRecoilControl,
        balance: RecoilControl,
        source_control: SharedWallControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SharedRecoilHeatUpdate<'a>, WallExchangeError> {
        poll(&mut cancelled)?;
        balance.wall.validate()?;
        if cells.len() != self.inputs.len() || source_control.maximum_source_evaluations == 0 {
            return Err(WallExchangeError::Invalid(
                "shared recoil thermal shape or source budget",
            ));
        }
        let work = Work {
            used: Counter::new(0),
            contacts: Counter::new(0),
            limit: source_control.maximum_source_evaluations,
        };
        work.charge()?;
        let source = self.thermal.accepted.storage.heat_storage.as_ref().unwrap();
        let initial_wall = bind_wall(source, thermal)?;
        let storage_constraint = source.constraint();
        for (i, cell) in cells.iter().enumerate() {
            poll(&mut cancelled)?;
            if !nonnegative(cell.conductance_w_per_k)
                || !cell.law.gas_heat_fraction.is_finite()
                || !(0.0..=1.0).contains(&cell.law.gas_heat_fraction)
            {
                return Err(WallExchangeError::Invalid(
                    "shared recoil conductance or heat partition",
                ));
            }
            work.charge()?;
            let primitive = gas.primitive(cell.state, &mut cancelled)?;
            self.inputs[i] = SharedRecoilCell {
                density_kg_m3: primitive.density,
                momentum_kg_m2_s: cell.state.momentum,
                volume_m3: cell.volume_m3,
                linear_rate_per_s: cell.law.linear_rate_per_s,
                quadratic_per_m: cell.law.quadratic_per_m,
            };
            self.thermal.volumes[i] = cell.volume_m3;
            self.thermal.accepted.ports[i].conductance_w_per_k = cell.conductance_w_per_k;
        }
        let motion = self.mechanics.advance(
            &self.inputs,
            duration,
            support,
            control,
            balance,
            &mut cancelled,
        )?;
        let mut wall_heat = Sum::default();
        for (i, (cell, mechanical)) in cells.iter().zip(motion.cells).enumerate() {
            poll(&mut cancelled)?;
            let gas_heat = cell.law.gas_heat_fraction * mechanical.slip_dissipation_j;
            let heat = mechanical.slip_dissipation_j - gas_heat;
            self.wall_heat[i] = heat;
            wall_heat.add(heat);
            let change = mechanical.work_to_gas_j - heat;
            let energy = cell.state.energy + change / cell.volume_m3;
            let residual = (energy - cell.state.energy) * cell.volume_m3 - change;
            if !energy.is_finite()
                || !residual.is_finite()
                || residual.abs()
                    > balance
                        .wall
                        .tolerance(change.abs().max(mechanical.slip_dissipation_j))
            {
                return Err(WallExchangeError::Balance {
                    energy_residual_j: residual,
                    thermal_residual_j: 0.0,
                });
            }
            self.thermal.initial[i] = [mechanical.momentum, energy];
        }
        let deposited = self
            .thermal
            .after_state(initial_wall.thermal, wall_heat.value(), &work)?;
        let deposit_residual = deposited.change - wall_heat.value();
        if !deposit_residual.is_finite()
            || deposit_residual.abs() > balance.wall.tolerance(wall_heat.value())
        {
            return Err(WallExchangeError::Balance {
                energy_residual_j: deposit_residual,
                thermal_residual_j: 0.0,
            });
        }
        let heat = self.thermal.solve(
            &|index, _| {
                let SharedWallGas::Plenum(i) = index else {
                    unreachable!("private local thermal indices");
                };
                mixture(gas, cells[i].state.partial_densities)
            },
            duration,
            balance.wall,
            &work,
            deposited.thermal,
            &mut cancelled,
        )?;
        work.charge()?;
        let source = self.thermal.accepted.storage.heat_storage.as_ref().unwrap();
        let wall_change = source.energy_change_j(
            initial_wall,
            WallContact {
                drag: WallDrag::None,
                thermal: self.thermal.candidate,
            },
        )?;
        let mut gas_change = Sum::default();
        let mut scale = Sum::default();
        for (i, (cell, mechanical)) in cells.iter().zip(motion.cells).enumerate() {
            poll(&mut cancelled)?;
            let change = (self.thermal.energies[i] - cell.state.energy) * cell.volume_m3;
            let mut residual = Sum::default();
            for value in [
                change,
                -mechanical.work_to_gas_j,
                self.wall_heat[i],
                -self.thermal.heat[i],
            ] {
                residual.add(value);
            }
            let cell_scale = mechanical.work_to_gas_j.abs()
                + mechanical.slip_dissipation_j
                + self.thermal.heat[i].abs();
            if !change.is_finite()
                || !cell_scale.is_finite()
                || !residual.value().is_finite()
                || residual.value().abs() > balance.wall.tolerance(cell_scale)
            {
                return Err(WallExchangeError::Balance {
                    energy_residual_j: residual.value(),
                    thermal_residual_j: heat.thermal_residual_j,
                });
            }
            gas_change.add(change);
            scale.add(cell_scale);
        }
        let mut residual = gas_change;
        residual.add(motion.support_kinetic_energy_change_j);
        residual.add(wall_change);
        let scale = scale
            .value()
            .max(wall_change.abs())
            .max(motion.support_kinetic_energy_change_j.abs());
        if !wall_change.is_finite()
            || !scale.is_finite()
            || !residual.value().is_finite()
            || residual.value().abs() > balance.wall.tolerance(scale)
        {
            return Err(WallExchangeError::Balance {
                energy_residual_j: residual.value(),
                thermal_residual_j: heat.thermal_residual_j,
            });
        }
        poll(&mut cancelled)?;
        Ok(SharedRecoilHeatUpdate {
            mechanics: motion,
            energies: &self.thermal.energies,
            heat_to_gas_j: &self.thermal.heat,
            wall_slip_heating_j: &self.wall_heat,
            thermal: self.thermal.candidate,
            storage_constraint,
            gas_energy_change_j: gas_change.value(),
            wall_energy_change_j: wall_change,
            energy_residual_j: residual.value(),
            thermal_report: SharedWallReport {
                elapsed: duration,
                energy_residual_j: heat.energy_residual_j,
                maximum_wall_residual_j: heat.energy_residual_j.abs(),
                maximum_contact_residual_j: heat.thermal_residual_j.abs(),
                source_evaluations: work.used.get(),
                wall_iterations: heat.iterations,
                contact_solves: work.contacts.get(),
            },
        })
    }
}

fn bind_wall(
    source: &MaterialWallCaloric<'_>,
    thermal: WallThermal,
) -> Result<WallContact, WallExchangeError> {
    let contact = WallContact {
        drag: WallDrag::None,
        thermal,
    };
    contact.validate()?;
    match thermal {
        WallThermal::Lumped {
            conductance_w_per_k: 0.0,
            ..
        }
        | WallThermal::Caloric {
            conductance_w_per_k: 0.0,
            ..
        } => source.bind_contact(contact),
        _ => Err(WallExchangeError::Invalid(
            "shared recoil requires one finite wall with conductance on its ports",
        )),
    }
}

#[cfg(test)]
mod tests;
