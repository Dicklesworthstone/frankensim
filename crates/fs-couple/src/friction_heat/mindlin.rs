//! Retained source Mindlin history coupled to the original finite heat network.
//!
//! The source continuation proves that the SAME selected elastic/friction
//! claims are constant throughout each temperature path. Native mechanics then
//! needs one solve, independent of thermal Newton iterations. Every heat trial
//! still resolves actual bulk and interface temperatures, admits the original
//! sources and transfers the complete retained history. A changed storage law
//! or friction bound refuses; no thermal query silently resets a contact.
//!
//! Each physical microslip loss enters its existing massless thermal face once.
//! Shared cells store energy once, while conductivities and caloric laws may
//! evolve at their actual endpoint temperatures. This fixed-normal/fixed-area
//! phase does not infer thermal expansion, wear geometry or moving coverage.

use super::network::{
    FrictionHeatNetworkInput, FrictionHeatNetworkReport, FrictionHeatNetworkStep,
    InterfaceHeatInput,
};
use super::{FrictionHeatRoute, FrictionHeatTolerances};
use fs_contact::interface_binding::mindlin::temperature::SourceMindlinTemperatureError;
use fs_contact::interface_binding::mindlin::{SourceMindlinContact, SourceMindlinState};
use fs_exec::CancelGate;
use fs_math::det;
use fs_mbd::tangential_friction::mindlin::multiple::{
    PreparedRigidMindlinIsland, RigidMindlinIslandCandidate, RigidMindlinIslandContact,
    RigidMindlinIslandControls, RigidMindlinIslandError,
};
use fs_mbd::{MassProperties, RigidBodyState, Wrench};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};
use std::cell::Cell;
use std::collections::BTreeMap;

/// One original source/history and its same physical circular thermal face.
#[derive(Clone, Copy)]
pub struct RigidMindlinIslandHeatContact<'a> {
    pub contact: RigidMindlinIslandContact<'a>,
    /// Ordered cells belong to the original mechanical contact's A/B bodies.
    pub route: &'a FrictionHeatRoute,
    /// Initial numerical guess, not a temperature boundary condition.
    pub temperature_guess_k: f64,
}

/// Refreshed original claims and transport, in canonical thermal-link order.
pub struct RigidMindlinIslandHeatInterface {
    pub link_id: u64,
    /// Same original selected claims at actual endpoint bulk/interface T.
    /// The source continuation independently proves invariant mechanics.
    pub law: SourceMindlinContact,
    pub conductivity_w_m_k: [f64; 2],
}

pub struct RigidMindlinIslandHeatSource<P> {
    pub interfaces: Vec<RigidMindlinIslandHeatInterface>,
    /// Other heat only; each bound interface already owns its friction/conduction.
    pub transport: StepCoefficients<P>,
}

#[derive(Clone, Copy, Debug)]
pub struct RigidMindlinIslandHeatControls {
    pub mechanical: RigidMindlinIslandControls,
    pub maximum_cells: usize,
    pub maximum_unknowns: usize,
    pub maximum_trials: usize,
    pub maximum_newton_iterations: usize,
    pub difference_step_j: f64,
    /// Numerical coordinate scaling [J/K], not a physical heat capacity.
    pub temperature_scale_j_per_k: f64,
    pub heat: FrictionHeatTolerances,
    pub combined_energy_tolerance_j: f64,
    pub newton: NewtonKrylovConfig,
    pub audit_entropy: bool,
}

#[derive(Debug)]
pub enum RigidMindlinIslandHeatError<E> {
    Invalid(&'static str),
    Mechanical(RigidMindlinIslandError),
    MindlinSource(SourceMindlinTemperatureError),
    Thermal(DepositionError),
    Source(E),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    TrialBudget,
}
impl<E> From<RigidMindlinIslandError> for RigidMindlinIslandHeatError<E> {
    fn from(value: RigidMindlinIslandError) -> Self {
        Self::Mechanical(value)
    }
}
impl<E> From<SourceMindlinTemperatureError> for RigidMindlinIslandHeatError<E> {
    fn from(value: SourceMindlinTemperatureError) -> Self {
        Self::MindlinSource(value)
    }
}
impl<E> From<DepositionError> for RigidMindlinIslandHeatError<E> {
    fn from(value: DepositionError) -> Self {
        Self::Thermal(value)
    }
}
impl<E> From<NewtonError> for RigidMindlinIslandHeatError<E> {
    fn from(value: NewtonError) -> Self {
        Self::Newton(value)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for RigidMindlinIslandHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "retained Mindlin island/heat: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for RigidMindlinIslandHeatError<E> {}

/// Endpoint source in canonical link order, explicitly mapped to native history.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidMindlinIslandHeatReceipt {
    pub link_id: u64,
    pub contact_id: String,
    pub law: SourceMindlinContact,
}

/// Publish original body momenta, retained histories, sources and heat together.
pub struct HeatedRigidMindlinIslandCandidate<'m, P> {
    pub mechanical: RigidMindlinIslandCandidate,
    pub thermal: DepositionSimulation<'m>,
    pub heat: FrictionHeatNetworkReport,
    pub sources: Vec<RigidMindlinIslandHeatReceipt>,
    pub source: P,
    pub nonlinear: NewtonReport,
    pub trials: usize,
    pub target_time_s: f64,
    /// Actual delta(K + retained shear storage + caloric energy), minus all
    /// external mechanical work and outside heat. Friction is internal once.
    pub combined_energy_residual_j: f64,
    pub summed_energy_residual_j: f64,
}

struct Staged<P> {
    sources: Vec<RigidMindlinIslandHeatReceipt>,
    histories: Vec<SourceMindlinState>,
    source: P,
}

/// Original shared-body mechanics and a simultaneous finite heat network.
#[allow(clippy::too_many_arguments)]
pub fn propose_rigid_mindlin_island_heat<'m, P, E>(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contacts: &[RigidMindlinIslandHeatContact<'_>],
    thermal: &DepositionSimulation<'m>,
    accepted_time_s: f64,
    dt_s: f64,
    controls: RigidMindlinIslandHeatControls,
    energy_guess_j: &[f64],
    resolve: impl Fn(
        &DepositionSimulation<'m>,
        &[u64],
        &[f64],
    ) -> Result<RigidMindlinIslandHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<HeatedRigidMindlinIslandCandidate<'m, P>, RigidMindlinIslandHeatError<E>> {
    propose(
        bodies,
        properties,
        contacts,
        None,
        thermal,
        accepted_time_s,
        dt_s,
        controls,
        energy_guess_j,
        resolve,
        gate,
    )
}

/// Held world COM forces and body-frame torques follow original global indices.
/// Once-per-body native work is an external supply, never frictional heat.
#[allow(clippy::too_many_arguments)]
pub fn propose_rigid_mindlin_island_driven_heat<'m, P, E>(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contacts: &[RigidMindlinIslandHeatContact<'_>],
    external: &[Wrench],
    thermal: &DepositionSimulation<'m>,
    accepted_time_s: f64,
    dt_s: f64,
    controls: RigidMindlinIslandHeatControls,
    energy_guess_j: &[f64],
    resolve: impl Fn(
        &DepositionSimulation<'m>,
        &[u64],
        &[f64],
    ) -> Result<RigidMindlinIslandHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<HeatedRigidMindlinIslandCandidate<'m, P>, RigidMindlinIslandHeatError<E>> {
    propose(
        bodies,
        properties,
        contacts,
        Some(external),
        thermal,
        accepted_time_s,
        dt_s,
        controls,
        energy_guess_j,
        resolve,
        gate,
    )
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn propose<'m, P, E>(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contacts: &[RigidMindlinIslandHeatContact<'_>],
    external: Option<&[Wrench]>,
    thermal: &DepositionSimulation<'m>,
    accepted_time_s: f64,
    dt_s: f64,
    controls: RigidMindlinIslandHeatControls,
    energy_guess_j: &[f64],
    resolve: impl Fn(
        &DepositionSimulation<'m>,
        &[u64],
        &[f64],
    ) -> Result<RigidMindlinIslandHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<HeatedRigidMindlinIslandCandidate<'m, P>, RigidMindlinIslandHeatError<E>> {
    let checkpoint = || {
        if gate.is_requested() {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    };
    checkpoint()?;
    controls.newton.validate()?;
    controls
        .mechanical
        .multiple
        .validate()
        .map_err(RigidMindlinIslandError::from)?;
    let count = contacts.len();
    let cells = thermal.cells().len();
    let dimension = cells
        .checked_add(count)
        .ok_or(RigidMindlinIslandHeatError::Invalid(
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
        || count > controls.mechanical.multiple.maximum_contacts
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
        || energy_guess_j.iter().any(|x| !x.is_finite())
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
        .any(|x| !x.is_finite() || *x < 0.0)
    {
        return Err(RigidMindlinIslandHeatError::Invalid(
            "bounded synchronized native bodies and finite heat required",
        ));
    }
    let duration_s = target - accepted_time_s;
    let mut ordered: Vec<_> = contacts.iter().collect();
    ordered.sort_by_key(|c| c.route.link_id());
    if ordered
        .windows(2)
        .any(|c| c[0].route.link_id() == c[1].route.link_id())
    {
        return Err(RigidMindlinIslandHeatError::Invalid(
            "one friction owner per original heat face",
        ));
    }
    let native: Vec<_> = ordered.iter().map(|c| c.contact).collect();
    let prepared = if let Some(external) = external {
        PreparedRigidMindlinIsland::new_driven(
            bodies,
            properties,
            &native,
            external,
            duration_s,
            controls.mechanical,
            &mut || gate.is_requested(),
        )?
    } else {
        PreparedRigidMindlinIsland::new(
            bodies,
            properties,
            &native,
            duration_s,
            controls.mechanical,
            &mut || gate.is_requested(),
        )?
    };
    let mut cell_owners = BTreeMap::new();
    for c in &ordered {
        checkpoint()?;
        validate_geometry(
            bodies,
            properties,
            c,
            controls.mechanical.geometry_tolerance_m,
        )
        .map_err(RigidMindlinIslandHeatError::Invalid)?;
        validate_temperatures(c.contact.law, thermal, c.route, None)?;
        let elastic = c
            .contact
            .law
            .normal_source()
            .model()
            .normal_state()
            .elastic();
        let bulk = [elastic.surface_a_state(), elastic.surface_b_state()];
        for (side, cell) in c.route.cells().into_iter().enumerate() {
            let owner = (
                c.contact.body_indices[side],
                bulk[side].resolved().card_identity(),
            );
            if cell_owners
                .insert(cell, owner)
                .is_some_and(|old| old != owner)
            {
                return Err(RigidMindlinIslandHeatError::Invalid(
                    "shared thermal cell changed original native body/material",
                ));
            }
        }
    }
    // Invariant selected source claims make this exact native solve independent
    // of all thermal coordinates. Each later source trial must prove that fact.
    let mut mechanical = prepared.advance(&mut || gate.is_requested())?;
    let mechanical_indices: Vec<_> = ordered
        .iter()
        .map(|c| {
            mechanical
                .contacts
                .binary_search_by(|p| p.contact_id.as_str().cmp(c.contact.contact_id))
                .map_err(|_| {
                    RigidMindlinIslandHeatError::Invalid("missing original native contact")
                })
        })
        .collect::<Result<_, _>>()?;
    let routes: Vec<_> = ordered.iter().map(|c| c.route).collect();
    let step = FrictionHeatNetworkStep::new(
        thermal,
        target,
        &routes,
        controls.mechanical.multiple.maximum_contacts,
        checkpoint,
    )?;
    let trials = Cell::new(0usize);
    let evaluate = |x: &[f64], entropy: bool| {
        checkpoint()?;
        if trials.get() >= controls.maximum_trials {
            return Err(RigidMindlinIslandHeatError::TrialBudget);
        }
        trials.set(trials.get() + 1);
        if x.len() != dimension || x.iter().any(|v| !v.is_finite()) {
            return Err(RigidMindlinIslandHeatError::Invalid(
                "complete finite thermal iterate required",
            ));
        }
        let temperatures: Vec<_> = ordered
            .iter()
            .enumerate()
            .map(|(i, c)| c.temperature_guess_k + x[cells + i] / controls.temperature_scale_j_per_k)
            .collect();
        let endpoint = |endpoint: &DepositionSimulation<'m>, temperatures: &[f64]| {
            let RigidMindlinIslandHeatSource {
                interfaces,
                transport,
            } = resolve(endpoint, step.link_ids(), temperatures)
                .map_err(RigidMindlinIslandHeatError::Source)?;
            if interfaces.len() != count
                || interfaces
                    .iter()
                    .zip(step.link_ids())
                    .any(|(s, id)| s.link_id != *id)
            {
                return Err(RigidMindlinIslandHeatError::Invalid(
                    "complete original link-ordered source claims required",
                ));
            }
            let mut inputs = Vec::with_capacity(count);
            let mut histories = Vec::with_capacity(count);
            let mut sources = Vec::with_capacity(count);
            for (((c, source), &index), &temperature) in ordered
                .iter()
                .zip(interfaces)
                .zip(&mechanical_indices)
                .zip(temperatures)
            {
                checkpoint()?;
                validate_temperatures(&source.law, endpoint, c.route, Some(temperature))?;
                let history = c.contact.law.continue_at_temperature(
                    &mechanical.contacts[index].contact.state,
                    &source.law,
                    &mut || gate.is_requested(),
                )?;
                inputs.push(InterfaceHeatInput {
                    dissipated_j: mechanical.contacts[index].contact.dissipated_j,
                    conductivity_w_m_k: source.conductivity_w_m_k,
                });
                histories.push(history);
                sources.push(RigidMindlinIslandHeatReceipt {
                    link_id: source.link_id,
                    contact_id: c.contact.contact_id.into(),
                    law: source.law,
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
                        sources,
                        histories,
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
        return Err(RigidMindlinIslandHeatError::Interrupted(result.progress));
    }
    if !result.report.converged {
        return Err(RigidMindlinIslandHeatError::NotConverged(result.report));
    }
    let trial = evaluate(&solver.x, controls.audit_entropy)?;
    let native = &mechanical.report;
    let heat = trial.report();
    let combined =
        native.kinetic_change_j + native.stored_change_j + heat.thermal.stored_energy_change_j
            - native.external_work_j
            - heat.external_heat_j;
    let summed = native.total_energy_residual_j + heat.total_energy_residual_j;
    if [combined, summed]
        .iter()
        .any(|v| !v.is_finite() || v.abs() > controls.combined_energy_tolerance_j)
    {
        return Err(RigidMindlinIslandHeatError::Invalid(
            "native storage/kinetic/caloric/external work does not close",
        ));
    }
    let (thermal, staged, heat) = trial.accept(controls.heat, checkpoint)?;
    for (&index, history) in mechanical_indices.iter().zip(staged.histories) {
        mechanical.contacts[index].contact.state = history;
    }
    checkpoint()?;
    Ok(HeatedRigidMindlinIslandCandidate {
        mechanical,
        thermal,
        heat,
        sources: staged.sources,
        source: staged.source,
        nonlinear: result.report,
        trials: trials.get(),
        target_time_s: target,
        combined_energy_residual_j: combined,
        summed_energy_residual_j: summed,
    })
}

fn validate_temperatures<E>(
    law: &SourceMindlinContact,
    thermal: &DepositionSimulation<'_>,
    route: &FrictionHeatRoute,
    interface: Option<f64>,
) -> Result<(), RigidMindlinIslandHeatError<E>> {
    if interface.is_some_and(|t| t.to_bits() != law.temperature_kelvin().to_bits()) {
        return Err(RigidMindlinIslandHeatError::Invalid(
            "source interface T differs from solved common face",
        ));
    }
    let elastic = law.normal_source().model().normal_state().elastic();
    for (bulk, cell) in [elastic.surface_a_state(), elastic.surface_b_state()]
        .into_iter()
        .zip(route.cells())
    {
        let temperature = thermal.temperature_k(cell)?;
        if !bulk
            .resolved()
            .query_point()
            .iter()
            .any(|(name, t)| name == "T" && t.to_bits() == temperature.to_bits())
        {
            return Err(RigidMindlinIslandHeatError::Invalid(
                "source bulk T differs from its actual original finite cell",
            ));
        }
    }
    Ok(())
}

fn validate_geometry(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    c: &RigidMindlinIslandHeatContact<'_>,
    tolerance: f64,
) -> Result<(), &'static str> {
    let geometry = c.route.geometry;
    let center = [geometry[0].start_m, geometry[0].end_m]
        .into_iter()
        .find(|p| *p == geometry[1].start_m || *p == geometry[1].end_m)
        .ok_or("original thermal face has no shared center")?;
    let radius = det::sqrt(c.route.contact_area_m2() / core::f64::consts::PI);
    if !radius.is_finite()
        || (radius - c.contact.law.normal_response().patch_radius_m).abs() > tolerance
    {
        return Err("original Mindlin footprint differs from circular finite heat face");
    }
    let normal = c.contact.frame.normal();
    for side in 0..2 {
        let index = c.contact.body_indices[side];
        let point = bodies[index]
            .point_kinematics(properties[index], c.contact.arms_body[side])
            .map_err(|_| "invalid original native heat attachment")?
            .point_world;
        if norm([
            point.x - center[0],
            point.y - center[1],
            point.z - center[2],
        ]) > tolerance
        {
            return Err("original mechanical contact and thermal face centers differ");
        }
        let axis =
            std::array::from_fn::<_, 3, _>(|i| geometry[side].end_m[i] - geometry[side].start_m[i]);
        if norm([
            axis[1] * normal[2] - axis[2] * normal[1],
            axis[2] * normal[0] - axis[0] * normal[2],
            axis[0] * normal[1] - axis[1] * normal[0],
        ]) > tolerance
        {
            return Err("original normal and finite-cell centrelines differ");
        }
    }
    Ok(())
}
fn norm(v: [f64; 3]) -> f64 {
    det::hypot(det::hypot(v[0], v[1]), v[2])
}
fn numerical<E>(
    error: FiniteDifferenceError<RigidMindlinIslandHeatError<E>>,
) -> RigidMindlinIslandHeatError<E> {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(field) => RigidMindlinIslandHeatError::Invalid(field),
    }
}
