//! Original finite thermal stores at the two impermeable moving-gas ends.

use super::super::{ReactingEuler1d, ReactiveError, ReactiveStepControl, ReactiveStepReport};
use crate::dynamics::piston::{Piston, PistonControl, PistonLoad, PistonStep};
use crate::multicomponent::dissipation::{DissipationControl, DissipationOperator, TransportLaw};
use crate::multicomponent::wall_exchange::{
    MaterialWallCaloric, WallCaloricConstraint, WallCellReport, WallContact, WallControl, WallDrag,
    WallExchangeError, WallThermal, exchange_cell, exchange_cell_with_heat_storage,
};
use crate::multicomponent::{
    Boundaries, Boundary, Mesh1d, MixtureEuler1d, MixtureGas, Totals, sum,
};
use fs_qty::Time;

/// One end-heat half interval, using the original local caloric solves.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PistonThermalReport {
    /// Actual physical duration, without a second flow-clock advance.
    pub elapsed: Time,
    /// Left and right end exchanges, in that order. A one-cell interval uses
    /// the same gas cell twice; the right solve sees the completed left solve.
    pub ends: [WallCellReport; 2],
}

/// Complete end-heat transfers surrounding molecular, chemical and Euler work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SplitPistonThermalReport {
    /// Heat at the original physical volumes, before molecular/chemical work.
    pub first: PistonThermalReport,
    /// Heat at the moved physical volumes, after molecular/chemical work.
    pub second: PistonThermalReport,
    /// Actual total conductive heat entering the gas [J].
    pub heat_to_gas_j: f64,
    /// Finite-wall internal-energy change between accepted endpoints [J].
    /// Original nonlinear source differences are used, never a frozen tangent.
    pub wall_energy_change_j: f64,
    /// Heat from prescribed infinite reservoirs only [J].
    pub external_heat_to_gas_j: f64,
    /// Gas heat plus finite-wall storage change minus reservoir heat [J].
    pub energy_residual_j: f64,
}

/// Complete thermal/reacting piston field report, with separate heat and work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalReactiveStepReport {
    /// The original D/C/A/C/D report between the heat halves. Its before/after
    /// states and conservation residual exclude the two outer heat intervals.
    pub reactive: ReactiveStepReport,
    /// Complete accepted field before either heat half.
    pub before: Totals,
    /// Complete endpoint field after both heat halves.
    pub after: Totals,
    /// Both original heat solves and the independent finite-wall storage change.
    pub heat: SplitPistonThermalReport,
    /// Actual physical gas chemical-energy change between full endpoints [J].
    /// This includes species-reference changes exactly once and includes heat.
    pub gas_energy_change_j: f64,
}

/// Two independent finite thermal stores or prescribed reservoirs at gas ends.
///
/// End order is `[left, right]`; either end may move. Conductance is the complete hA
/// in W/K, supplied by the caller; it is not a material or boundary-layer
/// correlation. Finite thermal masses are distinct, fixed amounts. An end's
/// thermal mass may represent part of its piston, whose total translating mass
/// remains with its mechanical owner; no second mass or work is added here.
/// Wall thermal expansion and deformation are outside this rigid-body model.
///
/// The original gas/source and accepted mesh are bound once. Candidate wall
/// states and one candidate mesh are reused through the existing piston root;
/// source laws are retained, not cloned into each velocity proposal. Only the
/// complete accepted gas/mechanical transaction publishes contacts and reports.
#[derive(Debug)]
pub struct PistonThermalEnds<'law> {
    gas: MixtureGas,
    mesh: Mesh1d,
    candidate_mesh: Mesh1d,
    contacts: [WallContact; 2],
    candidate_contacts: [WallContact; 2],
    caloric: [Option<MaterialWallCaloric<'law>>; 2],
    last_report: Option<SplitPistonThermalReport>,
}

impl<'law> PistonThermalEnds<'law> {
    /// Bind two original thermal contacts to a closed constant-area interval.
    ///
    /// Pass `[None, None]` for authored adiabatic/reservoir/constant-capacity
    /// walls. A source-backed contact must be constructed from its same caloric
    /// law. Internal-energy storage is admitted directly; the source's explicit
    /// stationary-solid/negligible-expansion convention is also admitted.
    /// Unqualified source enthalpy refuses because its missing pressure work
    /// cannot be inserted into a gas-plus-piston internal-energy balance.
    /// Radial drag refuses: longitudinal wall stress belongs to the original
    /// moving molecular operator and must not be applied a second time here.
    pub fn new(
        flow: &MixtureEuler1d,
        mut contacts: [WallContact; 2],
        caloric: [Option<MaterialWallCaloric<'law>>; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, WallExchangeError> {
        poll(&mut cancelled)?;
        closed(flow)?;
        for (contact, source) in contacts.iter_mut().zip(&caloric) {
            poll(&mut cancelled)?;
            contact.validate()?;
            if contact.drag != WallDrag::None {
                return Err(WallExchangeError::Invalid(
                    "thermal piston end cannot supply radial drag or duplicate wall stress",
                ));
            }
            match source {
                Some(source) => {
                    if source.constraint() == WallCaloricConstraint::SourceEnthalpy {
                        return Err(WallExchangeError::Invalid(
                            "thermal piston requires a rigid internal-energy storage convention",
                        ));
                    }
                    *contact = source.bind_contact(*contact)?;
                }
                None if matches!(contact.thermal, WallThermal::Caloric { .. }) => {
                    return Err(WallExchangeError::Invalid(
                        "thermal piston energy contact requires its original caloric source",
                    ));
                }
                None => {}
            }
        }
        let mesh = flow
            .mesh()
            .try_copy()
            .map_err(crate::multicomponent::moving::geometry_error)?;
        let candidate_mesh = mesh
            .try_copy()
            .map_err(crate::multicomponent::moving::geometry_error)?;
        poll(&mut cancelled)?;
        Ok(Self {
            gas: flow.gas().clone(),
            mesh,
            candidate_mesh,
            contacts,
            candidate_contacts: contacts,
            caloric,
            last_report: None,
        })
    }

    /// Accepted thermal contacts, including actual finite-wall temperatures/u.
    pub const fn contacts(&self) -> &[WallContact; 2] {
        &self.contacts
    }

    /// Original immutable source owners, in left/right end order.
    pub const fn heat_storage(&self) -> &[Option<MaterialWallCaloric<'law>>; 2] {
        &self.caloric
    }

    /// Exact accepted gas geometry; failed moved candidates never rebind it.
    pub const fn mesh(&self) -> &Mesh1d {
        &self.mesh
    }

    /// Heat report from the last jointly accepted complete piston step.
    pub const fn last_report(&self) -> Option<SplitPistonThermalReport> {
        self.last_report
    }

    pub(in crate::multicomponent::reactive) fn admit(
        &self,
        flow: &MixtureEuler1d,
    ) -> Result<(), WallExchangeError> {
        closed(flow)?;
        if self.gas != *flow.gas() || self.mesh != *flow.mesh() {
            return Err(WallExchangeError::Invalid(
                "thermal piston and accepted flow source or geometry differ",
            ));
        }
        Ok(())
    }

    pub(in crate::multicomponent::reactive) fn begin_trial(&mut self) {
        self.candidate_contacts = self.contacts;
    }

    pub(in crate::multicomponent::reactive) fn advance_trial(
        &mut self,
        flow: &mut MixtureEuler1d,
        duration: Time,
        controls: WallControl,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PistonThermalReport, WallExchangeError> {
        poll(cancelled)?;
        self.candidate_mesh
            .reset_motion_from(flow.mesh(), cancelled)
            .map_err(crate::multicomponent::moving::geometry_error)?;
        let n = self.gas.species_count();
        let right = self.mesh.cell_count() - 1;
        let mut exchange = |end: usize, cell: usize| {
            poll(cancelled)?;
            let contact = self.candidate_contacts[end];
            let state = flow.cell(cell).expect("admitted thermal end cell");
            let update = match &self.caloric[end] {
                Some(source) => exchange_cell_with_heat_storage(
                    &self.gas,
                    state,
                    flow.mesh.volumes()[cell],
                    duration,
                    contact,
                    source,
                    controls,
                    &mut *cancelled,
                )?,
                None => exchange_cell(
                    &self.gas,
                    state,
                    flow.mesh.volumes()[cell],
                    duration,
                    contact,
                    controls,
                    &mut *cancelled,
                )?,
            };
            // Drag is forbidden. Heat cannot invent impulse, mechanical work,
            // a species source or a kinetic-energy correction at either end.
            if update.momentum != state.momentum || update.report.wall_impulse_ns != 0.0 {
                return Err(WallExchangeError::Invalid("end heat changed gas momentum"));
            }
            flow.states[cell * (n + 2) + n + 1] = update.energy;
            self.candidate_contacts[end] = update.contact;
            Ok(update.report)
        };
        let ends = [exchange(0, 0)?, exchange(1, right)?];
        Ok(PistonThermalReport {
            elapsed: duration,
            ends,
        })
    }

    pub(in crate::multicomponent::reactive) fn complete(
        &self,
        first: PistonThermalReport,
        second: PistonThermalReport,
    ) -> Result<SplitPistonThermalReport, WallExchangeError> {
        let mut changes = [0.0; 2];
        for (end, change) in changes.iter_mut().enumerate() {
            let before = self.contacts[end];
            let after = self.candidate_contacts[end];
            *change = match &self.caloric[end] {
                Some(source) => source.energy_change_j(before, after)?,
                None => match (before.thermal, after.thermal) {
                    (
                        WallThermal::Lumped {
                            temperature_k: t0,
                            heat_capacity_j_per_k,
                            ..
                        },
                        WallThermal::Lumped {
                            temperature_k: t1, ..
                        },
                    ) => heat_capacity_j_per_k * (t1 - t0),
                    _ => 0.0,
                },
            };
        }
        let stages = first.ends.iter().chain(&second.ends);
        let heat_to_gas_j = sum(stages.clone().map(|s| s.heat_to_gas_j));
        let wall_energy_change_j = sum(changes);
        let external_heat_to_gas_j = sum(stages.map(|s| s.external_heat_to_gas_j));
        let energy_residual_j = sum([heat_to_gas_j, wall_energy_change_j, -external_heat_to_gas_j]);
        if ![
            heat_to_gas_j,
            wall_energy_change_j,
            external_heat_to_gas_j,
            energy_residual_j,
        ]
        .iter()
        .all(|v| v.is_finite())
        {
            return Err(WallExchangeError::Invalid(
                "non-finite complete piston heat balance",
            ));
        }
        Ok(SplitPistonThermalReport {
            first,
            second,
            heat_to_gas_j,
            wall_energy_change_j,
            external_heat_to_gas_j,
            energy_residual_j,
        })
    }

    pub(in crate::multicomponent::reactive) fn publish(
        &mut self,
        report: SplitPistonThermalReport,
    ) {
        self.contacts = self.candidate_contacts;
        core::mem::swap(&mut self.mesh, &mut self.candidate_mesh);
        self.last_report = Some(report);
    }
}

impl ReactingEuler1d {
    /// Advance finite end-wall heat, molecular transport, chemistry and piston motion.
    ///
    /// Every proposal in the original midpoint-speed root runs H(dt/2), D(dt/2),
    /// C(dt/2), A(dt), C(dt/2), D(dt/2), H(dt/2) from the same accepted gas and
    /// wall stores. H uses the original scalar caloric exchange at the old/new
    /// end-cell volumes. Both ends use left-then-right order, including for one
    /// gas cell. These backward-Euler heat stages do not establish second-order
    /// accuracy or a fully implicit conjugate heat-transfer solution.
    ///
    /// Pressure and longitudinal viscous impulse/work alone drive the mechanical
    /// root. Conductive heat changes the actual gas caloric state and subsequent
    /// pressure flux. The full accepted energy residual is gas chemical energy
    /// plus finite-wall storage plus piston kinetic change plus exterior work,
    /// minus additional-force work and prescribed-reservoir heat. Finite-wall
    /// heat is internal to that balance, never external work or duplicate heat.
    ///
    /// Shortened durations restart the entire root and all seven stages. Gas,
    /// mesh, clock, molecular geometry, piston, thermal contacts and report publish
    /// together after the final poll. Original source domains, chemistry and
    /// parabolic budgets remain authoritative; end hA is caller-prescribed.
    ///
    /// # Errors
    /// Any source/geometry, nonlinear, balance or cancellation refusal leaves all
    /// accepted owners unchanged and reusable, including finite thermal stores.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_piston_dissipation_and_heat<L: TransportLaw>(
        &mut self,
        operator: &mut DissipationOperator<L>,
        thermal: &mut PistonThermalEnds<'_>,
        piston: &mut Piston,
        maximum_duration: Time,
        controls: ReactiveStepControl,
        dissipation_controls: DissipationControl,
        thermal_controls: WallControl,
        load: PistonLoad,
        piston_controls: PistonControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PistonStep<ThermalReactiveStepReport>, ReactiveError> {
        self.step_with_piston_molecular(
            operator,
            Some((thermal, thermal_controls)),
            piston,
            maximum_duration,
            controls,
            dissipation_controls,
            load,
            piston_controls,
            cancelled,
        )
        .map(|report| {
            report.map_gas(|gas| ThermalReactiveStepReport {
                reactive: gas.flow,
                before: gas.before,
                after: gas.after,
                heat: gas
                    .heat
                    .expect("thermal piston entry always supplies heat stages"),
                gas_energy_change_j: gas.gas_energy_change_j,
            })
        })
    }
}

fn closed(flow: &MixtureEuler1d) -> Result<(), WallExchangeError> {
    if flow.mesh().uniform_cross_section().is_none()
        || flow.boundaries()
            != (Boundaries::Ends {
                left: Boundary::Reflecting,
                right: Boundary::Reflecting,
            })
    {
        return Err(WallExchangeError::Invalid(
            "thermal piston requires a constant-area closed interval",
        ));
    }
    Ok(())
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), WallExchangeError> {
    if cancelled() {
        Err(WallExchangeError::Cancelled)
    } else {
        Ok(())
    }
}
