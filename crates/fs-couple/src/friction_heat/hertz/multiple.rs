//! Simultaneous source-backed Hertz patches and one finite caloric network.
//!
//! One immutable native island retains every shared-body cross response. Each
//! trial resolves every original coefficient source at that patch's own solved
//! interface temperature, then advances all mechanical patches together. The
//! original heat network sums their powers into shared cells exactly once.

use super::{RigidHertzHeatContact, validate_geometry};
use crate::friction_heat::network::{
    FrictionHeatNetworkInput, FrictionHeatNetworkReport, FrictionHeatNetworkStep,
    InterfaceHeatInput,
};
use crate::friction_heat::{FrictionHeatRoute, FrictionHeatTolerances};
use fs_contact::interface_binding::thermal_friction::hertz::{
    SourceHertzWrench, SourceHertzWrenchError,
};
use fs_exec::CancelGate;
use fs_mbd::tangential_friction::hertz::RigidHertzContact;
use fs_mbd::tangential_friction::hertz::multiple::{
    PreparedRigidHertzIsland, RigidHertzIslandCandidate, RigidHertzIslandContact,
    RigidHertzIslandControls, RigidHertzIslandError, RigidHertzIslandLaw,
};
use fs_mbd::{MassProperties, RigidBodyState};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};
use std::cell::Cell;

/// One original mechanical patch, pinned source and distinct thermal face.
/// The route's A/B order belongs to the contact binding's A/B order. Original
/// body/material identity and normal evolution remain with their owners.
#[derive(Clone, Copy)]
pub struct RigidHertzIslandHeatContact<'a> {
    pub patch_id: &'a str,
    pub contact: RigidHertzHeatContact<'a>,
    pub route: &'a FrictionHeatRoute,
    /// Positive initial guess only; never a prescribed physical temperature.
    pub temperature_guess_k: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct RigidHertzIslandHeatConductivity {
    pub link_id: u64,
    pub conductivity_w_m_k: [f64; 2],
}

pub struct RigidHertzIslandHeatSource<P> {
    /// Exactly the callback's ascending original link order, including IDs.
    pub interfaces: Vec<RigidHertzIslandHeatConductivity>,
    /// Other heat only. The bound friction and conductance already have owners.
    pub transport: StepCoefficients<P>,
}

#[derive(Clone, Copy, Debug)]
pub struct RigidHertzIslandHeatControls {
    pub mechanical: RigidHertzIslandControls,
    pub maximum_cells: usize,
    /// Bound cell increments plus distinct interface temperatures before allocation.
    pub maximum_unknowns: usize,
    pub maximum_trials: usize,
    pub maximum_newton_iterations: usize,
    pub difference_step_j: f64,
    /// Numerical scaling [J/K], not an added physical heat capacity.
    pub temperature_scale_j_per_k: f64,
    pub heat: FrictionHeatTolerances,
    pub combined_energy_tolerance_j: f64,
    pub newton: NewtonKrylovConfig,
    /// Final-root exact caloric entropy only; every active material must support it.
    pub audit_entropy: bool,
}

#[derive(Debug)]
pub enum RigidHertzIslandHeatError<E> {
    Invalid(&'static str),
    Mechanical(RigidHertzIslandError),
    FrictionSource(SourceHertzWrenchError),
    Thermal(DepositionError),
    Source(E),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    TrialBudget,
}
impl<E> From<RigidHertzIslandError> for RigidHertzIslandHeatError<E> {
    fn from(error: RigidHertzIslandError) -> Self {
        Self::Mechanical(error)
    }
}
impl<E> From<SourceHertzWrenchError> for RigidHertzIslandHeatError<E> {
    fn from(error: SourceHertzWrenchError) -> Self {
        Self::FrictionSource(error)
    }
}
impl<E> From<DepositionError> for RigidHertzIslandHeatError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl<E> From<NewtonError> for RigidHertzIslandHeatError<E> {
    fn from(error: NewtonError) -> Self {
        Self::Newton(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for RigidHertzIslandHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Hertz island friction/heat: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for RigidHertzIslandHeatError<E> {}

/// Original coefficient-use receipt paired with its stable mechanical identity.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidHertzIslandHeatReceipt {
    pub patch_id: String,
    pub friction_source: SourceHertzWrench,
}

/// Every native body, source receipt and finite thermal endpoint publishes together.
pub struct HeatedRigidHertzIslandCandidate<'m, P> {
    pub mechanical: RigidHertzIslandCandidate,
    pub thermal: DepositionSimulation<'m>,
    /// Canonical ascending original link order, with one heat owner per patch.
    pub heat: FrictionHeatNetworkReport,
    /// Same canonical thermal order; explicit patch IDs retain mechanical mapping.
    pub friction_sources: Vec<RigidHertzIslandHeatReceipt>,
    pub source: P,
    pub nonlinear: NewtonReport,
    pub trials: usize,
    pub target_time_s: f64,
    /// Direct actual kinetic plus stored caloric change, minus outside heat.
    pub combined_energy_residual_j: f64,
    /// Independent sum of native mechanics and original network audit residuals.
    pub summed_energy_residual_j: f64,
}

struct Staged<P> {
    mechanical: RigidHertzIslandCandidate,
    friction_sources: Vec<RigidHertzIslandHeatReceipt>,
    source: P,
}

/// Solve all source-resolved Hertz patches and shared finite cells simultaneously.
///
/// Input contacts may be permuted; their stable IDs and distinct original link
/// IDs preserve the same canonical mechanics and thermal solve. Each patch has
/// its own original source and solved common temperature. Nominal pressure is
/// explicit; T/p-only source admission refuses other evolving source axes.
/// The callback resolves original endpoint conductivities and other heat only.
/// It must stage any external receipts until this complete candidate is accepted.
///
/// Poses, normal loads, contact radii and cell geometries are frozen. The native
/// island's response-rank, effort and unresolved-endpoint-event refusals remain
/// in force. This stateless Hertz rung cannot consume elastic contact history
/// or rolling loss. The original orthogonal finite-cell/common-face model and
/// backward Euler heat need spatial/time refinement. Nonlinear convergence
/// does not establish a unique root or physical validity outside these models.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn propose_rigid_hertz_island_heat<'m, P, E>(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contacts: &[RigidHertzIslandHeatContact<'_>],
    thermal: &DepositionSimulation<'m>,
    accepted_time_s: f64,
    dt_s: f64,
    controls: RigidHertzIslandHeatControls,
    energy_guess_j: &[f64],
    resolve: impl Fn(
        &DepositionSimulation<'m>,
        &[u64],
        &[f64],
    ) -> Result<RigidHertzIslandHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<HeatedRigidHertzIslandCandidate<'m, P>, RigidHertzIslandHeatError<E>> {
    let checkpoint = || {
        if gate.is_requested() {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    };
    checkpoint()?;
    controls.newton.validate()?;
    controls.mechanical.validate()?;
    let count = contacts.len();
    let cells = thermal.cells().len();
    let dimension = cells
        .checked_add(count)
        .ok_or(RigidHertzIslandHeatError::Invalid(
            "thermal dimension overflow",
        ))?;
    let target = accepted_time_s + dt_s;
    if !accepted_time_s.is_finite()
        || accepted_time_s != thermal.time_s()
        || !dt_s.is_finite()
        || dt_s <= 0.0
        || !target.is_finite()
        || target <= accepted_time_s
        || count == 0
        || count > controls.mechanical.maximum_contacts
        || controls.maximum_cells == 0
        || cells > controls.maximum_cells
        || dimension > controls.maximum_unknowns
        || controls.maximum_trials == 0
        || controls.maximum_newton_iterations == 0
        || !controls.difference_step_j.is_finite()
        || controls.difference_step_j <= 0.0
        || !controls.temperature_scale_j_per_k.is_finite()
        || controls.temperature_scale_j_per_k <= 0.0
        || energy_guess_j.len() != cells
        || energy_guess_j.iter().any(|v| !v.is_finite())
        || contacts
            .iter()
            .any(|c| !c.temperature_guess_k.is_finite() || c.temperature_guess_k <= 0.0)
        || [
            controls.heat.cell_j,
            controls.heat.thermal_j,
            controls.heat.interface_j,
            controls.heat.total_j,
            controls.combined_energy_tolerance_j,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v < 0.0)
    {
        return Err(RigidHertzIslandHeatError::Invalid(
            "synchronized bounded bodies, patches and solved thermal coordinates required",
        ));
    }
    // Both owners use the same representable interval, including large clocks.
    let duration_s = target - accepted_time_s;
    let mut ordered: Vec<_> = contacts.iter().collect();
    ordered.sort_by_key(|c| c.route.link_id());
    if ordered
        .windows(2)
        .any(|p| p[0].route.link_id() == p[1].route.link_id())
    {
        return Err(RigidHertzIslandHeatError::Invalid(
            "each original thermal link requires one mechanical heat owner",
        ));
    }
    let mut initial = Vec::with_capacity(count);
    for c in &ordered {
        checkpoint()?;
        initial.push(c.contact.source.hertz_at(
            c.temperature_guess_k,
            c.contact.nominal_pressure_pa,
            c.contact.normal_force_n,
            c.contact.contact_radius_m,
        )?);
    }
    let native_contacts: Vec<_> = ordered
        .iter()
        .zip(&initial)
        .map(|(c, source)| RigidHertzIslandContact {
            patch_id: c.patch_id,
            contact: RigidHertzContact {
                binding: c.contact.binding,
                frame: c.contact.frame,
                law: source.law(),
            },
        })
        .collect();
    let prepared = PreparedRigidHertzIsland::new(
        bodies,
        properties,
        &native_contacts,
        duration_s,
        controls.mechanical,
        gate,
    )?;
    let mut mechanical_indices = Vec::with_capacity(count);
    for c in &ordered {
        checkpoint()?;
        let index = prepared
            .contacts()
            .binary_search_by(|p| p.patch_id.as_str().cmp(c.patch_id))
            .map_err(|_| RigidHertzIslandHeatError::Invalid("missing native patch identity"))?;
        mechanical_indices.push(index);
        let patch = &prepared.contacts()[index];
        let point = patch.contact_points_world[0];
        validate_geometry(
            [point.x, point.y, point.z],
            c.contact,
            c.route,
            controls.mechanical.mechanical.geometry_tolerance_m,
        )
        .map_err(RigidHertzIslandHeatError::Invalid)?;
    }
    let routes: Vec<_> = ordered.iter().map(|c| c.route).collect();
    let step = FrictionHeatNetworkStep::new(
        thermal,
        target,
        &routes,
        controls.mechanical.maximum_contacts,
        checkpoint,
    )?;
    let trials = Cell::new(0usize);
    let evaluate = |x: &[f64], entropy: bool| {
        checkpoint()?;
        if trials.get() >= controls.maximum_trials {
            return Err(RigidHertzIslandHeatError::TrialBudget);
        }
        trials.set(trials.get() + 1);
        if x.len() != dimension || x.iter().any(|v| !v.is_finite()) {
            return Err(RigidHertzIslandHeatError::Invalid(
                "complete finite thermal iterate required",
            ));
        }
        let temperatures: Vec<_> = ordered
            .iter()
            .enumerate()
            .map(|(i, c)| c.temperature_guess_k + x[cells + i] / controls.temperature_scale_j_per_k)
            .collect();
        let endpoint = |endpoint: &DepositionSimulation<'m>, temperatures: &[f64]| {
            let mut receipts = Vec::with_capacity(count);
            for (c, &temperature) in ordered.iter().zip(temperatures) {
                checkpoint()?;
                receipts.push(RigidHertzIslandHeatReceipt {
                    patch_id: c.patch_id.to_owned(),
                    friction_source: c.contact.source.hertz_at(
                        temperature,
                        c.contact.nominal_pressure_pa,
                        c.contact.normal_force_n,
                        c.contact.contact_radius_m,
                    )?,
                });
            }
            let laws: Vec<_> = receipts
                .iter()
                .map(|r| RigidHertzIslandLaw {
                    patch_id: &r.patch_id,
                    law: r.friction_source.law(),
                })
                .collect();
            let mechanical = prepared.advance(&laws, gate)?;
            let RigidHertzIslandHeatSource {
                interfaces,
                transport,
            } = resolve(endpoint, step.link_ids(), temperatures)
                .map_err(RigidHertzIslandHeatError::Source)?;
            if interfaces.len() != count
                || interfaces
                    .iter()
                    .zip(step.link_ids())
                    .any(|(source, id)| source.link_id != *id)
            {
                return Err(RigidHertzIslandHeatError::Invalid(
                    "conductivities must retain every canonical original link identity",
                ));
            }
            let mut inputs = Vec::with_capacity(count);
            for ((c, &index), source) in ordered.iter().zip(&mechanical_indices).zip(interfaces) {
                checkpoint()?;
                let native = &mechanical.report.contacts[index];
                if native.patch_id != c.patch_id {
                    return Err(RigidHertzIslandHeatError::Invalid(
                        "changed original Hertz heat owner",
                    ));
                }
                inputs.push(InterfaceHeatInput {
                    dissipated_j: native.dissipated_energy_j,
                    conductivity_w_m_k: source.conductivity_w_m_k,
                });
            }
            let StepCoefficients {
                sources_w,
                conductances_w_k,
                payload,
            } = transport;
            Ok(FrictionHeatNetworkInput {
                interfaces: inputs,
                transport: StepCoefficients {
                    sources_w,
                    conductances_w_k,
                    payload: Staged {
                        mechanical,
                        friction_sources: receipts,
                        source: payload,
                    },
                },
            })
        };
        if entropy {
            step.evaluate_with_entropy(&x[..cells], &temperatures, endpoint, checkpoint)
        } else {
            step.evaluate(&x[..cells], &temperatures, endpoint, checkpoint)
        }
    };
    let problem =
        FallibleFiniteDifference::new(dimension, controls.difference_step_j, |x: &[f64]| {
            evaluate(x, false).map(|trial| trial.residual_j().to_vec())
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
        return Err(RigidHertzIslandHeatError::Interrupted(result.progress));
    }
    if !result.report.converged {
        return Err(RigidHertzIslandHeatError::NotConverged(result.report));
    }
    let trial = evaluate(&solver.x, controls.audit_entropy)?;
    let native = &trial.payload().mechanical.report;
    let heat = trial.report();
    let combined = native.kinetic_after_j - native.kinetic_before_j
        + heat.thermal.stored_energy_change_j
        - heat.external_heat_j;
    let summed = native.total_energy_residual_j + heat.total_energy_residual_j;
    if [combined, summed]
        .iter()
        .any(|v| !v.is_finite() || v.abs() > controls.combined_energy_tolerance_j)
    {
        return Err(RigidHertzIslandHeatError::Invalid(
            "original island kinetic and finite-cell energy does not balance",
        ));
    }
    let (thermal, staged, heat) = trial.accept(controls.heat, checkpoint)?;
    checkpoint()?;
    Ok(HeatedRigidHertzIslandCandidate {
        mechanical: staged.mechanical,
        thermal,
        heat,
        friction_sources: staged.friction_sources,
        source: staged.source,
        nonlinear: result.report,
        trials: trials.get(),
        target_time_s: target,
        combined_energy_residual_j: combined,
        summed_energy_residual_j: summed,
    })
}

fn numerical<E>(
    error: FiniteDifferenceError<RigidHertzIslandHeatError<E>>,
) -> RigidHertzIslandHeatError<E> {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(field) => RigidHertzIslandHeatError::Invalid(field),
    }
}
