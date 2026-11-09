//! Finite per-cell support inertia on the original network source transaction.

use super::{
    DevelopedWallTransport, NetworkWallOperator, WallContact, WallControl, WallDrag,
    WallExchangeError, WallThermal, allocated, poll, positive, sum,
};
use crate::multicomponent::wall_exchange::{
    MaterialWallCaloric, RecoilControl, RecoilLaw, RecoilUpdate, SlidingSupport,
    WallCaloricConstraint, WallUpdate, exchange_cell_with_recoil,
    exchange_cell_with_recoil_and_heat_storage,
};
use crate::multicomponent::{Cell, MixtureGas};
use fs_qty::Time;

/// Mechanical exchanges with one finite support over all accepted source
/// halves. The support may belong to one cell or a complete branch. Units use
/// that branch's original axial coordinate. Heat transfer remains in the
/// original `NetworkWallReport`, separate from work.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SupportExchange {
    /// Sum of actual gas momentum changes from support friction [N s].
    pub gas_impulse_ns: f64,
    /// Actual support endpoint momentum change [N s].
    pub support_impulse_ns: f64,
    /// Mechanical work received by gas from this support [J].
    pub work_to_gas_j: f64,
    /// Actual gas energy changes during the friction stages only [J].
    pub gas_energy_change_j: f64,
    /// Actual support endpoint kinetic-energy change [J].
    pub kinetic_energy_change_j: f64,
    /// Relative-slip energy thermalized once across gas and finite wall [J].
    pub drag_heating_j: f64,
    /// Wall heating during friction stages, before subsequent thermal exchange
    /// [J]. Local walls report actual storage gain; shared supports report the
    /// assigned transfer, whose original shared inventory separately audits
    /// deposition and complete endpoints. This is an internal gas debit.
    pub wall_slip_heating_j: f64,
    /// Actual gas impulse plus support endpoint momentum change [N s].
    pub momentum_residual_ns: f64,
    /// Friction-stage gas, wall heating and actual support kinetic changes [J].
    /// Shared inventory representation and subsequent thermal residuals remain
    /// in the original network/shared-wall reports.
    pub energy_residual_j: f64,
}

pub(super) struct FiniteSupports {
    accepted: Vec<Option<SlidingSupport>>,
    trial: Vec<Option<SlidingSupport>>,
    exchanges: Vec<Option<SupportExchange>>,
    candidate_exchanges: Vec<Option<SupportExchange>>,
    gas_heat_fractions: Vec<f64>,
    momentum_tolerance_ns: f64,
}

impl NetworkWallOperator<'_> {
    /// Bind one independent finite sliding support to each selected branch cell.
    ///
    /// The nested shape is exactly the original branch/cell order. `None`
    /// retains the original prescribed/stationary wall. Plenums cannot acquire
    /// axial inertia because their original state has no bulk momentum.
    ///
    /// Each selected support replaces the contact's prescribed velocity while
    /// retaining its linear/Darcy/combined coefficients. Original circular
    /// material transport refreshes coefficients and checks Reynolds/Mach
    /// against the actual support-relative slip at BOTH endpoints.
    ///
    /// Mechanics reuses the exact reduced-mass recoil law. By default all slip
    /// loss stays in gas sensible energy; [`Self::with_finite_support_heat_partition`]
    /// can deposit a declared share into each support's finite wall. The
    /// original local/caloric/shared-wall/conduction owner then performs heat
    /// exchange. Shared THERMAL storage remains possible, but supports are
    /// mechanically independent: repeated descriptors do not create one shared
    /// mass. No wall geometry, mesh volume, external drive or clock is added.
    ///
    /// Bind before evolution. Support states and complete reaction receipts
    /// carry across both original source halves and publish only with the
    /// network; every complete retry restores them from accepted state.
    pub fn with_finite_supports(
        mut self,
        branches: Vec<Vec<Option<SlidingSupport>>>,
        momentum_tolerance_ns: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, WallExchangeError> {
        poll(&mut cancelled)?;
        if self.supports.is_some()
            || self.last_report.is_some()
            || self.trial_active
            || branches.len() != self.branches.len()
            || !positive(momentum_tolerance_ns)
        {
            return Err(WallExchangeError::Invalid(
                "finite supports must bind once to original branches before evolution",
            ));
        }
        let mut accepted = allocated(self.plenum_start, None)?;
        for (branch, (shape, input)) in self.branches.iter().zip(&branches).enumerate() {
            if input.len() != shape.mesh.cell_count() {
                return Err(WallExchangeError::Invalid("finite support cell axis"));
            }
            for (cell, support) in input.iter().copied().enumerate() {
                poll(&mut cancelled)?;
                if let Some(support) = support {
                    if self
                        .shared_supports
                        .as_ref()
                        .is_some_and(|supports| supports.has_branch(branch))
                    {
                        return Err(WallExchangeError::Invalid(
                            "independent and shared supports overlap",
                        ));
                    }
                    if !positive(support.mass_kg) || !support.velocity_m_s.is_finite() {
                        return Err(WallExchangeError::Invalid(
                            "finite support mass or velocity",
                        ));
                    }
                    self.contacts[shape.start + cell].drag =
                        at_velocity(self.contacts[shape.start + cell].drag, support.velocity_m_s);
                }
                accepted[shape.start + cell] = support;
            }
        }
        let count = accepted.len();
        let mut trial = allocated(count, None)?;
        trial.copy_from_slice(&accepted);
        let exchanges = allocated(count, None)?;
        let candidate_exchanges = allocated(count, None)?;
        let gas_heat_fractions = allocated(count, 1.0)?;
        poll(&mut cancelled)?;
        self.supports = Some(FiniteSupports {
            accepted,
            trial,
            exchanges,
            candidate_exchanges,
            gas_heat_fractions,
            momentum_tolerance_ns,
        });
        Ok(self)
    }

    /// Partition frictional slip heat between gas and each independent support's
    /// original finite wall. Each entry is the fraction retained in gas, in
    /// `[0, 1]`; a cell without an independent support must retain `1`.
    ///
    /// Bind after the supports and their local material storage, before any
    /// evolution. A fraction below one requires a finite local wall, with
    /// explicit internal energy or negligible expansion for material sources.
    /// One shared wall inventory cannot be substituted for these distinct local
    /// stores. The wall's thermal amount already belongs to the support inertia.
    ///
    /// Friction heats the original capacity integral or native energy coordinate
    /// before the existing thermal/conduction solve. No constant tangent or new
    /// heat-transfer model is introduced. All allocations occur at binding;
    /// retry and cancellation retain the original network transaction.
    pub fn with_finite_support_heat_partition(
        mut self,
        branches: Vec<Vec<f64>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, WallExchangeError> {
        poll(&mut cancelled)?;
        if self.last_report.is_some() || self.trial_active || branches.len() != self.branches.len()
        {
            return Err(WallExchangeError::Invalid(
                "finite support heat partition must bind before evolution",
            ));
        }
        let supports = self.supports.as_mut().ok_or(WallExchangeError::Invalid(
            "heat partition requires independent finite supports",
        ))?;
        for (shape, input) in self.branches.iter().zip(branches) {
            if input.len() != shape.mesh.cell_count() {
                return Err(WallExchangeError::Invalid(
                    "support heat partition cell axis",
                ));
            }
            for (cell, fraction) in input.into_iter().enumerate() {
                poll(&mut cancelled)?;
                let index = shape.start + cell;
                if !fraction.is_finite() || !(0.0..=1.0).contains(&fraction) {
                    return Err(WallExchangeError::Invalid("support gas heat fraction"));
                }
                if fraction < 1.0 {
                    if supports.accepted[index].is_none() {
                        return Err(WallExchangeError::Invalid(
                            "wall slip heat requires an independent finite support",
                        ));
                    }
                    friction_wall_conductance(self.contacts[index], self.caloric[index].as_ref())?;
                }
                supports.gas_heat_fractions[index] = fraction;
            }
        }
        poll(&mut cancelled)?;
        Ok(self)
    }

    /// Accepted finite supports in one original branch's cell order. This is
    /// read-only: changing inertia or velocity requires an explicit new owner.
    pub fn branch_supports(&self, branch: usize) -> Option<&[Option<SlidingSupport>]> {
        let shape = self.branches.get(branch)?;
        let supports = self.supports.as_ref()?;
        Some(&supports.accepted[shape.start..shape.start + shape.mesh.cell_count()])
    }

    /// Complete accepted mechanical exchanges in one branch's cell order.
    /// `None` cells retain their original prescribed/stationary wall reaction
    /// in `accepted_branch_impulses`. Different branch axes are never summed.
    pub fn accepted_branch_support_exchanges(
        &self,
        branch: usize,
    ) -> Option<&[Option<SupportExchange>]> {
        self.last_report?;
        let shape = self.branches.get(branch)?;
        let supports = self.supports.as_ref()?;
        Some(&supports.exchanges[shape.start..shape.start + shape.mesh.cell_count()])
    }
}

impl FiniteSupports {
    pub(super) fn any_selected(&self, start: usize, count: usize) -> bool {
        self.accepted[start..start + count]
            .iter()
            .any(Option::is_some)
    }

    pub(super) fn begin_trial(
        &mut self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), WallExchangeError> {
        for i in 0..self.accepted.len() {
            poll(cancelled)?;
            self.trial[i] = self.accepted[i];
            self.candidate_exchanges[i] = self.accepted[i].map(|_| SupportExchange::default());
        }
        Ok(())
    }

    pub(super) fn contact(&self, index: usize, mut contact: WallContact) -> WallContact {
        if let Some(support) = self.trial[index] {
            contact.drag = at_velocity(contact.drag, support.velocity_m_s);
        }
        contact
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn prepare(
        &mut self,
        index: usize,
        branch: usize,
        transport: Option<&mut DevelopedWallTransport>,
        gas: &MixtureGas,
        state: Cell<'_>,
        volume: f64,
        duration: Time,
        contact: WallContact,
        wall: Option<&MaterialWallCaloric<'_>>,
        controls: WallControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Option<RecoilUpdate>, WallExchangeError> {
        let Some(support) = self.trial[index] else {
            return Ok(None);
        };
        let mut law = coefficients(contact.drag);
        law.gas_heat_fraction = self.gas_heat_fractions[index];
        let balance = RecoilControl {
            wall: controls,
            momentum_tolerance_ns: self.momentum_tolerance_ns,
        };
        let conductance = if law.gas_heat_fraction == 1.0 {
            None
        } else {
            Some(friction_wall_conductance(contact, wall)?)
        };
        let mut thermal = contact.thermal;
        let wall = if conductance.is_some() {
            // Deposit friction only. The existing network owner applies its
            // material-conduction and heat exchange after the mechanical stage.
            set_finite_conductance(&mut thermal, 0.0)?;
            wall
        } else {
            thermal = WallThermal::Adiabatic;
            None
        };
        let mut motion = match transport.filter(|transport| transport.couples_drag(branch)) {
            Some(transport) => transport.advance_recoil(
                gas,
                branch,
                state,
                volume,
                duration,
                support,
                law.gas_heat_fraction,
                WallContact { thermal, ..contact },
                wall,
                balance,
                &mut *cancelled,
            )?,
            None => match wall {
                Some(wall) => exchange_cell_with_recoil_and_heat_storage(
                    gas,
                    state,
                    volume,
                    duration,
                    support,
                    law,
                    thermal,
                    wall,
                    balance,
                    &mut *cancelled,
                )?,
                None => exchange_cell_with_recoil(
                    gas,
                    state,
                    volume,
                    duration,
                    support,
                    law,
                    thermal,
                    balance,
                    &mut *cancelled,
                )?,
            },
        };
        if let Some(conductance) = conductance {
            set_finite_conductance(&mut motion.gas.contact.thermal, conductance)?;
        } else {
            motion.gas.contact.thermal = contact.thermal;
        }
        motion.gas.contact.drag = at_velocity(contact.drag, motion.support.velocity_m_s);
        let initial = self.accepted[index].ok_or(WallExchangeError::Invalid(
            "finite support lost its accepted inventory",
        ))?;
        let previous = self.candidate_exchanges[index].unwrap_or_default();
        let gas_impulse = sum([previous.gas_impulse_ns, -motion.gas.report.wall_impulse_ns]);
        let work = sum([
            previous.work_to_gas_j,
            motion.gas.report.mechanical_work_to_gas_j,
        ]);
        let gas_energy = sum([
            previous.gas_energy_change_j,
            (motion.gas.energy - state.energy) * volume,
        ]);
        // Audit actual complete support endpoints rather than summing kinetic
        // deltas from independently rounded half-step endpoints.
        let impulse = initial.mass_kg * (motion.support.velocity_m_s - initial.velocity_m_s);
        let kinetic = impulse * (0.5 * motion.support.velocity_m_s + 0.5 * initial.velocity_m_s);
        let heating = sum([previous.drag_heating_j, motion.slip_dissipation_j]);
        let wall_heating = sum([
            previous.wall_slip_heating_j,
            motion.gas.report.wall_energy_change_j,
        ]);
        let report = SupportExchange {
            gas_impulse_ns: gas_impulse,
            support_impulse_ns: impulse,
            work_to_gas_j: work,
            gas_energy_change_j: gas_energy,
            kinetic_energy_change_j: kinetic,
            drag_heating_j: heating,
            wall_slip_heating_j: wall_heating,
            momentum_residual_ns: gas_impulse + impulse,
            energy_residual_j: gas_energy + kinetic + wall_heating,
        };
        if ![
            gas_impulse,
            impulse,
            work,
            gas_energy,
            kinetic,
            heating,
            wall_heating,
            report.momentum_residual_ns,
            report.energy_residual_j,
        ]
        .into_iter()
        .all(f64::is_finite)
            || heating < 0.0
            || report.momentum_residual_ns.abs() > self.momentum_tolerance_ns
        {
            return Err(WallExchangeError::Invalid(
                "complete finite support momentum balance",
            ));
        }
        if report.energy_residual_j.abs() > controls.tolerance(work.abs().max(heating)) {
            return Err(WallExchangeError::Balance {
                energy_residual_j: report.energy_residual_j,
                thermal_residual_j: 0.0,
            });
        }
        self.trial[index] = Some(motion.support);
        self.candidate_exchanges[index] = Some(report);
        Ok(Some(motion))
    }

    pub(super) fn commit(&mut self) {
        core::mem::swap(&mut self.accepted, &mut self.trial);
        core::mem::swap(&mut self.exchanges, &mut self.candidate_exchanges);
    }
}

/// Retain exact one-cell mechanics in the common mechanical-update format.
pub(super) fn mechanical_update(motion: RecoilUpdate, drag: WallDrag) -> WallUpdate {
    let mut update = motion.gas;
    update.contact.drag = at_velocity(drag, motion.support.velocity_m_s);
    update.report.drag_heating_j = motion.slip_dissipation_j;
    update
}

/// Restore the coefficient source and accepted support velocity after the
/// original thermal owner ran with drag disabled. Its storage, tangent,
/// source identity and conduction admission remain authoritative.
pub(super) fn finish(
    mut thermal: WallUpdate,
    motion: WallUpdate,
    initial_contact: WallContact,
    wall: Option<&MaterialWallCaloric<'_>>,
    initial_energy: f64,
    volume: f64,
    controls: WallControl,
) -> Result<WallUpdate, WallExchangeError> {
    thermal.contact.drag = motion.contact.drag;
    thermal.report.wall_impulse_ns = motion.report.wall_impulse_ns;
    thermal.report.mechanical_work_to_gas_j = motion.report.mechanical_work_to_gas_j;
    thermal.report.drag_heating_j = motion.report.drag_heating_j;
    thermal.report.heat_to_gas_j += motion.report.heat_to_gas_j;
    // Audit complete endpoints, including friction and subsequent conduction,
    // rather than summing independently rounded changes of the wall store.
    thermal.report.wall_energy_change_j = match wall {
        Some(wall) => wall.energy_change_j(initial_contact, thermal.contact)?,
        None => match (initial_contact.thermal, thermal.contact.thermal) {
            (
                WallThermal::Lumped {
                    temperature_k: initial,
                    heat_capacity_j_per_k,
                    ..
                },
                WallThermal::Lumped {
                    temperature_k: final_temperature,
                    ..
                },
            ) => heat_capacity_j_per_k * (final_temperature - initial),
            // An adiabatic local contact may deliver its frictional wall
            // fraction to one separately owned shared inventory. Its exact
            // deposition and complete endpoint are audited by that owner.
            _ => thermal.report.wall_energy_change_j + motion.report.wall_energy_change_j,
        },
    };
    thermal.report.energy_residual_j = (thermal.energy - initial_energy) * volume
        + thermal.report.wall_energy_change_j
        - thermal.report.external_heat_to_gas_j
        - thermal.report.mechanical_work_to_gas_j;
    let scale = thermal
        .report
        .mechanical_work_to_gas_j
        .abs()
        .max(thermal.report.heat_to_gas_j.abs())
        .max(thermal.report.wall_energy_change_j.abs());
    if !thermal.report.heat_to_gas_j.is_finite()
        || !thermal.report.wall_energy_change_j.is_finite()
        || !thermal.report.energy_residual_j.is_finite()
        || thermal.report.energy_residual_j.abs() > controls.tolerance(scale)
    {
        return Err(WallExchangeError::Balance {
            energy_residual_j: thermal.report.energy_residual_j,
            thermal_residual_j: thermal.report.thermal_residual_j,
        });
    }
    Ok(thermal)
}

pub(super) fn friction_wall_conductance(
    contact: WallContact,
    wall: Option<&MaterialWallCaloric<'_>>,
) -> Result<f64, WallExchangeError> {
    if let Some(wall) = wall {
        if !matches!(
            wall.constraint(),
            WallCaloricConstraint::IsochoricInternalEnergy
                | WallCaloricConstraint::StationarySolidNegligibleExpansion
        ) {
            return Err(WallExchangeError::Invalid(
                "support wall heating requires internal energy or explicit negligible expansion",
            ));
        }
        wall.bind_contact(contact)?;
    }
    match contact.thermal {
        WallThermal::Lumped {
            conductance_w_per_k,
            ..
        } => Ok(conductance_w_per_k),
        WallThermal::Caloric {
            conductance_w_per_k,
            ..
        } if wall.is_some() => Ok(conductance_w_per_k),
        _ => Err(WallExchangeError::Invalid(
            "support slip heat requires original finite wall storage",
        )),
    }
}

fn set_finite_conductance(thermal: &mut WallThermal, value: f64) -> Result<(), WallExchangeError> {
    match thermal {
        WallThermal::Lumped {
            conductance_w_per_k,
            ..
        }
        | WallThermal::Caloric {
            conductance_w_per_k,
            ..
        } => *conductance_w_per_k = value,
        _ => {
            return Err(WallExchangeError::Invalid(
                "finite support wall storage was lost",
            ));
        }
    }
    Ok(())
}

pub(super) fn coefficients(drag: WallDrag) -> RecoilLaw {
    let (linear_rate_per_s, quadratic_per_m) = drag.resistance_coefficients();
    RecoilLaw {
        linear_rate_per_s,
        quadratic_per_m,
        gas_heat_fraction: 1.0,
    }
}

pub(super) fn at_velocity(drag: WallDrag, wall_velocity_m_s: f64) -> WallDrag {
    match drag {
        WallDrag::None => WallDrag::MovingLinear {
            rate_per_s: 0.0,
            wall_velocity_m_s,
        },
        WallDrag::Linear { rate_per_s } | WallDrag::MovingLinear { rate_per_s, .. } => {
            WallDrag::MovingLinear {
                rate_per_s,
                wall_velocity_m_s,
            }
        }
        WallDrag::Darcy {
            friction_factor,
            hydraulic_diameter_m,
        }
        | WallDrag::MovingDarcy {
            friction_factor,
            hydraulic_diameter_m,
            ..
        } => WallDrag::MovingDarcy {
            friction_factor,
            hydraulic_diameter_m,
            wall_velocity_m_s,
        },
        WallDrag::LinearQuadratic {
            linear_rate_per_s,
            quadratic_per_m,
            ..
        } => WallDrag::LinearQuadratic {
            linear_rate_per_s,
            quadratic_per_m,
            wall_velocity_m_s,
        },
    }
}
