//! One retained mechanical island and all of its finite thermal interfaces.

use super::{
    RigidFrictionHeatControls, RigidFrictionHeatError, compatible_law, numerical, validate_geometry,
};
use crate::friction_heat::FrictionHeatRoute;
use crate::friction_heat::network::{
    FrictionHeatNetworkInput, FrictionHeatNetworkReport, FrictionHeatNetworkStep,
    InterfaceHeatInput,
};
use fs_contact::tangential::nonsmooth::{NonsmoothContact, NonsmoothTangentialState};
use fs_exec::CancelGate;
use fs_mbd::tangential_friction::{RigidContactBinding, RigidFrictionReport, apply_rigid_friction};
use fs_mbd::{MassProperties, RigidBodyState};
use fs_solver::finite_difference::FallibleFiniteDifference;
use fs_solver::{NewtonKrylovState, NewtonReport, SolverRunProgress};
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};
use fs_tribo::partial_slip::PartialSlipLaw;
use std::cell::Cell;

/// An original loaded patch and its ordered physical body/thermal binding.
pub struct HeatedRigidContact<'a> {
    pub contact: &'a NonsmoothContact<'a>,
    pub binding: RigidContactBinding,
    pub route: &'a FrictionHeatRoute,
    /// Initial nonlinear guess only, never a temperature boundary condition.
    pub temperature_guess_k: f64,
}

/// One temperature-resolved original law, explicitly identified by thermal link.
pub struct RigidInterfaceHeatSource {
    pub link_id: u64,
    pub law: PartialSlipLaw,
    pub conductivity_w_m_k: [f64; 2],
}

pub struct RigidFrictionHeatNetworkSource<P> {
    /// Exactly the callback's canonical link order; each ID is checked.
    pub interfaces: Vec<RigidInterfaceHeatSource>,
    /// Outside heat and other original links only; friction is supplied here.
    pub transport: StepCoefficients<P>,
}

pub struct HeatedRigidFrictionNetworkReport<P> {
    pub mechanical: RigidFrictionReport,
    pub heat: FrictionHeatNetworkReport,
    /// Accepted original laws in canonical thermal link order.
    pub interfaces: Vec<RigidInterfaceHeatSource>,
    pub source: P,
    pub nonlinear: NewtonReport,
    pub trials: usize,
    pub combined_energy_residual_j: f64,
}

struct Staged<P> {
    bodies: Vec<RigidBodyState>,
    history: NonsmoothTangentialState,
    mechanical: RigidFrictionReport,
    interfaces: Vec<RigidInterfaceHeatSource>,
    source: P,
}

/// Accept all shared-body friction and shared-cell heat in one transaction.
///
/// Every Newton trial resolves all original laws at their distinct common
/// interface temperatures and then calls the native rigid multicontact owner
/// exactly once. Only each regional return map's physical heat enters that
/// contact's original thermal link. Elastic storage, numerical losses and
/// outside supplies retain separate accounts. Multiple interfaces may share
/// bodies or finite thermal cells; no sequential intermediate state is used.
///
/// `resolve` receives ascending original thermal link IDs and their matching
/// temperatures, independent of caller contact order. Its receipt must remain
/// staged until this function succeeds. Source identity, regional stiffness and
/// geometry retain the single-patch consumer's strict continuation rules.
/// Frozen poses, positive normal patches and coincident orthogonal thermal
/// faces are required; normal solving and geometry motion belong to their owners.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn apply_rigid_friction_heat_network<'m, P, E>(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    history: &mut NonsmoothTangentialState,
    contacts: &[HeatedRigidContact<'_>],
    thermal: &mut DepositionSimulation<'m>,
    dt_s: f64,
    controls: RigidFrictionHeatControls,
    energy_guess_j: &[f64],
    resolve: impl Fn(
        &DepositionSimulation<'m>,
        &[u64],
        &[f64],
    ) -> Result<RigidFrictionHeatNetworkSource<P>, E>,
    gate: &CancelGate,
) -> Result<HeatedRigidFrictionNetworkReport<P>, RigidFrictionHeatError<E>> {
    let checkpoint = || {
        if gate.is_requested() {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    };
    checkpoint()?;
    controls.newton.validate()?;
    let target = history.time_s() + dt_s;
    if history.time_s() != thermal.time_s()
        || !dt_s.is_finite()
        || dt_s <= 0.0
        || !target.is_finite()
        || target <= history.time_s()
        || bodies.len() != properties.len()
        || bodies.len() > controls.mechanical.maximum_bodies
        || contacts.is_empty()
        || contacts.len() > controls.mechanical.maximum_contacts
        || controls.maximum_cells == 0
        || thermal.cells().len() > controls.maximum_cells
        || controls.maximum_trials == 0
        || controls.maximum_newton_iterations == 0
        || !controls.difference_step_j.is_finite()
        || controls.difference_step_j <= 0.0
        || !controls.temperature_scale_j_per_k.is_finite()
        || controls.temperature_scale_j_per_k <= 0.0
        || !controls.combined_energy_tolerance_j.is_finite()
        || controls.combined_energy_tolerance_j < 0.0
        || energy_guess_j.len() != thermal.cells().len()
        || energy_guess_j.iter().any(|x| !x.is_finite())
    {
        return Err(RigidFrictionHeatError::Invalid(
            "synchronized bounded nonempty rigid/thermal network required",
        ));
    }
    let mut ordered: Vec<_> = contacts.iter().collect();
    ordered.sort_by_key(|contact| contact.route.patch.link.id);
    for input in &ordered {
        checkpoint()?;
        if input.binding.body_a >= bodies.len()
            || input.binding.body_b >= bodies.len()
            || input.binding.body_a == input.binding.body_b
            || !input.temperature_guess_k.is_finite()
            || input.temperature_guess_k <= 0.0
        {
            return Err(RigidFrictionHeatError::Invalid(
                "distinct valid body indices and positive Kelvin guesses required",
            ));
        }
        validate_geometry(
            bodies,
            properties,
            input.contact,
            input.binding,
            input.route,
            controls.mechanical.geometry_tolerance_m,
        )?;
        if let Some(old) = history.contact_checkpoint(input.contact.patch.patch_id()) {
            compatible_law(old.law(), input.contact.law)?;
        }
    }
    let routes: Vec<_> = ordered.iter().map(|input| input.route).collect();
    let step = FrictionHeatNetworkStep::new(
        thermal,
        target,
        &routes,
        controls.mechanical.maximum_contacts,
        checkpoint,
    )?;
    let cells = step.cell_ids().len();
    let dimension = cells
        .checked_add(contacts.len())
        .ok_or(RigidFrictionHeatError::Invalid(
            "thermal dimension overflow",
        ))?;
    let trials = Cell::new(0usize);
    let evaluate = |x: &[f64]| {
        checkpoint()?;
        if trials.get() >= controls.maximum_trials {
            return Err(RigidFrictionHeatError::TrialBudget);
        }
        trials.set(trials.get() + 1);
        if x.len() != dimension {
            return Err(RigidFrictionHeatError::Invalid(
                "incomplete rigid/thermal network iterate",
            ));
        }
        let temperatures: Vec<_> = ordered
            .iter()
            .enumerate()
            .map(|(i, input)| {
                input.temperature_guess_k + x[cells + i] / controls.temperature_scale_j_per_k
            })
            .collect();
        step.evaluate(
            &x[..cells],
            &temperatures,
            |endpoint, temperatures| {
                let RigidFrictionHeatNetworkSource {
                    interfaces,
                    transport,
                } = resolve(endpoint, step.link_ids(), temperatures)
                    .map_err(RigidFrictionHeatError::Source)?;
                if interfaces.len() != ordered.len() {
                    return Err(RigidFrictionHeatError::Invalid(
                        "incomplete resolved rigid interface laws",
                    ));
                }
                for (i, source) in interfaces.iter().enumerate() {
                    if source.link_id != step.link_ids()[i] {
                        return Err(RigidFrictionHeatError::Invalid(
                            "resolved rigid interface law order differs",
                        ));
                    }
                    compatible_law(ordered[i].contact.law, &source.law)?;
                }
                let resolved: Vec<_> = ordered
                    .iter()
                    .zip(&interfaces)
                    .map(|(input, source)| {
                        let contact = input.contact;
                        NonsmoothContact {
                            law: &source.law,
                            patch: contact.patch,
                            interface: contact.interface,
                            frame: contact.frame,
                            common_spin_rad_per_s: contact.common_spin_rad_per_s,
                            torsional_length_m: contact.torsional_length_m,
                            creepage_reference_speed_mps: contact.creepage_reference_speed_mps,
                            work: contact.work.clone(),
                        }
                    })
                    .collect();
                let bindings: Vec<_> = ordered.iter().map(|input| input.binding).collect();
                let mut next_bodies = bodies.to_vec();
                let mut next_history = history.clone();
                let mechanical = apply_rigid_friction(
                    &mut next_bodies,
                    properties,
                    &mut next_history,
                    &resolved,
                    &bindings,
                    dt_s,
                    controls.mechanical,
                    gate,
                )?;
                if mechanical.friction.coefficient_work_j != 0.0 {
                    return Err(RigidFrictionHeatError::UnbalancedCoefficientWork);
                }
                let mut heat_inputs = Vec::with_capacity(interfaces.len());
                for (input, source) in ordered.iter().zip(&interfaces) {
                    let index = mechanical
                        .contacts
                        .binary_search_by(|result| {
                            result.patch_id.as_str().cmp(input.contact.patch.patch_id())
                        })
                        .map_err(|_| {
                            RigidFrictionHeatError::Invalid("missing original mechanical contact")
                        })?;
                    heat_inputs.push(InterfaceHeatInput {
                        dissipated_j: mechanical.contacts[index].law.dissipation.heat_j,
                        conductivity_w_m_k: source.conductivity_w_m_k,
                    });
                }
                let StepCoefficients {
                    sources_w,
                    conductances_w_k,
                    payload,
                } = transport;
                Ok(FrictionHeatNetworkInput {
                    interfaces: heat_inputs,
                    transport: StepCoefficients {
                        sources_w,
                        conductances_w_k,
                        payload: Staged {
                            bodies: next_bodies,
                            history: next_history,
                            mechanical,
                            interfaces,
                            source: payload,
                        },
                    },
                })
            },
            checkpoint,
        )
    };
    let problem =
        FallibleFiniteDifference::new(dimension, controls.difference_step_j, |x: &[f64]| {
            evaluate(x).map(|trial| trial.residual_j().to_vec())
        })
        .map_err(numerical)?;
    let mut initial = energy_guess_j.to_vec();
    initial.resize(dimension, 0.0);
    let construction = NewtonKrylovState::new(&problem, initial, controls.newton);
    if let Some(error) = problem.take_failure() {
        return Err(numerical(error));
    }
    let mut solver = construction?;
    let result = solver.run_controlled(&problem, controls.maximum_newton_iterations, &mut || {
        gate.is_requested() || problem.has_failed()
    });
    if let Some(error) = problem.take_failure() {
        return Err(numerical(error));
    }
    checkpoint()?;
    if result.progress != SolverRunProgress::Complete {
        return Err(RigidFrictionHeatError::Interrupted(result.progress));
    }
    if !result.report.converged {
        return Err(RigidFrictionHeatError::NotConverged(result.report));
    }
    let trial = evaluate(&solver.x)?;
    let mechanical = &trial.payload().mechanical;
    let heat = trial.report();
    let combined = mechanical.kinetic_after_j - mechanical.kinetic_before_j
        + mechanical.friction.storage_change_j
        + heat.thermal.stored_energy_change_j
        - heat.external_heat_j
        + mechanical.friction.elastic_numerical_loss_j
        + mechanical.friction.mechanical_numerical_loss_j;
    if !combined.is_finite() || combined.abs() > controls.combined_energy_tolerance_j {
        return Err(RigidFrictionHeatError::Invalid(
            "joint rigid/thermal network work does not balance",
        ));
    }
    let (next_thermal, staged, heat) = trial.accept(controls.heat, checkpoint)?;
    checkpoint()?;
    drop(problem);
    bodies.copy_from_slice(&staged.bodies);
    *history = staged.history;
    *thermal = next_thermal;
    Ok(HeatedRigidFrictionNetworkReport {
        mechanical: staged.mechanical,
        heat,
        interfaces: staged.interfaces,
        source: staged.source,
        nonlinear: result.report,
        trials: trials.get(),
        combined_energy_residual_j: combined,
    })
}
