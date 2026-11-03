//! Finite per-cell support inertia on the original network source transaction.

use super::{
    NetworkWallOperator, WallContact, WallControl, WallDrag, WallExchangeError, WallThermal,
    allocated, poll, positive, sum,
};
use crate::multicomponent::wall_exchange::{
    RecoilControl, RecoilLaw, RecoilUpdate, SlidingSupport, WallUpdate, exchange_cell_with_recoil,
};
use crate::multicomponent::{Cell, MixtureGas};
use fs_qty::Time;

/// Mechanical exchanges with one independent finite support over all accepted
/// source halves. Units use that branch's original axial coordinate. Heat
/// transfer remains in the original `NetworkWallReport`, separate from work.
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
    /// Relative-slip energy thermalized once in the gas [J].
    pub drag_heating_j: f64,
    /// Actual gas impulse plus support endpoint momentum change [N s].
    pub momentum_residual_ns: f64,
    /// Friction-stage gas energy changes plus actual support kinetic change [J].
    /// Original subsequent thermal storage/flux residuals remain separate.
    pub energy_residual_j: f64,
}

pub(super) struct FiniteSupports {
    accepted: Vec<Option<SlidingSupport>>,
    trial: Vec<Option<SlidingSupport>>,
    exchanges: Vec<Option<SupportExchange>>,
    candidate_exchanges: Vec<Option<SupportExchange>>,
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
    /// Mechanics reuses the exact reduced-mass recoil law. All slip loss stays
    /// in gas sensible energy; it is not deposited a second time as heat. The
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
        for (shape, input) in self.branches.iter().zip(&branches) {
            if input.len() != shape.mesh.cell_count() {
                return Err(WallExchangeError::Invalid("finite support cell axis"));
            }
            for (cell, support) in input.iter().copied().enumerate() {
                poll(&mut cancelled)?;
                if let Some(support) = support {
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
        poll(&mut cancelled)?;
        self.supports = Some(FiniteSupports {
            accepted,
            trial,
            exchanges,
            candidate_exchanges,
            momentum_tolerance_ns,
        });
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
        gas: &MixtureGas,
        state: Cell<'_>,
        volume: f64,
        duration: Time,
        drag: WallDrag,
        controls: WallControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Option<RecoilUpdate>, WallExchangeError> {
        let Some(support) = self.trial[index] else {
            return Ok(None);
        };
        let motion = exchange_cell_with_recoil(
            gas,
            state,
            volume,
            duration,
            support,
            coefficients(drag),
            WallThermal::Adiabatic,
            RecoilControl {
                wall: controls,
                momentum_tolerance_ns: self.momentum_tolerance_ns,
            },
            cancelled,
        )?;
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
        let report = SupportExchange {
            gas_impulse_ns: gas_impulse,
            support_impulse_ns: impulse,
            work_to_gas_j: work,
            gas_energy_change_j: gas_energy,
            kinetic_energy_change_j: kinetic,
            drag_heating_j: heating,
            momentum_residual_ns: gas_impulse + impulse,
            energy_residual_j: gas_energy + kinetic,
        };
        if ![
            gas_impulse,
            impulse,
            work,
            gas_energy,
            kinetic,
            heating,
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

/// Restore the original coefficient source and accepted support velocity after
/// the thermal owner ran with drag disabled. Its original heat storage, tangent,
/// source identity and conduction admission remain authoritative.
pub(super) fn finish(
    mut thermal: WallUpdate,
    motion: RecoilUpdate,
    drag: WallDrag,
    controls: WallControl,
) -> Result<WallUpdate, WallExchangeError> {
    thermal.contact.drag = at_velocity(drag, motion.support.velocity_m_s);
    thermal.report.wall_impulse_ns = motion.gas.report.wall_impulse_ns;
    thermal.report.mechanical_work_to_gas_j = motion.gas.report.mechanical_work_to_gas_j;
    thermal.report.drag_heating_j = motion.slip_dissipation_j;
    thermal.report.energy_residual_j += motion.gas.report.energy_residual_j;
    let scale = thermal
        .report
        .mechanical_work_to_gas_j
        .abs()
        .max(thermal.report.heat_to_gas_j.abs())
        .max(thermal.report.wall_energy_change_j.abs());
    if !thermal.report.energy_residual_j.is_finite()
        || thermal.report.energy_residual_j.abs() > controls.tolerance(scale)
    {
        return Err(WallExchangeError::Balance {
            energy_residual_j: thermal.report.energy_residual_j,
            thermal_residual_j: thermal.report.thermal_residual_j,
        });
    }
    Ok(thermal)
}

fn coefficients(drag: WallDrag) -> RecoilLaw {
    let (linear_rate_per_s, quadratic_per_m) = match drag {
        WallDrag::None => (0.0, 0.0),
        WallDrag::Linear { rate_per_s } | WallDrag::MovingLinear { rate_per_s, .. } => {
            (rate_per_s, 0.0)
        }
        WallDrag::Darcy {
            friction_factor,
            hydraulic_diameter_m,
        }
        | WallDrag::MovingDarcy {
            friction_factor,
            hydraulic_diameter_m,
            ..
        } => (0.0, 0.5 * friction_factor / hydraulic_diameter_m),
        WallDrag::LinearQuadratic {
            linear_rate_per_s,
            quadratic_per_m,
            ..
        } => (linear_rate_per_s, quadratic_per_m),
    };
    RecoilLaw {
        linear_rate_per_s,
        quadratic_per_m,
        gas_heat_fraction: 1.0,
    }
}

fn at_velocity(drag: WallDrag, wall_velocity_m_s: f64) -> WallDrag {
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
