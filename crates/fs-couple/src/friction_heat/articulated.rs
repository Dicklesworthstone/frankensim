//! Original articulated shear and finite-cell heat at one shared endpoint.
//!
//! Native contact geometry, full cross-contact mobility and objective history
//! transport are prepared once. Each thermal trial uses those same accepted
//! histories and original laws, with distinct solved interface temperatures.
//! Only their physical irreversible work enters the original thermal network.

use super::network::{
    FrictionHeatNetworkInput, FrictionHeatNetworkReport, FrictionHeatNetworkStep,
    InterfaceHeatInput,
};
use super::{FrictionHeatRoute, FrictionHeatTolerances};
use fs_exec::CancelGate;
use fs_math::det;
use fs_mbd::articulated_contacts::tangential::planar::multiple::ArticulatedMultiTangentialConfig;
use fs_mbd::articulated_contacts::tangential::planar::transport::{
    ArticulatedFramedMultiTangentialResult, ArticulatedFramedTangentialSite,
    PreparedArticulatedFramedInput, PreparedArticulatedFramedTangentialContacts,
};
use fs_mbd::articulated_contacts::{ArticulatedContactError, ArticulatedContactSystem};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};
use fs_tribo::tangential::transport::{ObjectiveTangentialLaw, TangentialTransportConfig};
use std::cell::Cell;

/// Explicit physical correspondence; thermal A belongs to the contact's point,
/// and thermal B to its dynamic counterpart or prescribed support.
#[derive(Clone, Copy)]
pub struct ArticulatedFrictionHeatRoute<'a> {
    pub contact_index: usize,
    pub route: &'a FrictionHeatRoute,
    /// Initial nonlinear guess only, never a prescribed contact temperature.
    pub temperature_guess_k: f64,
}

/// Original endpoint material conductivities in the route's A/B order.
#[derive(Clone, Copy, Debug)]
pub struct ArticulatedInterfaceHeatSource {
    pub link_id: u64,
    pub conductivity_w_m_k: [f64; 2],
}

pub struct ArticulatedFrictionHeatSource<P> {
    /// Exactly the callback's canonical ascending original link order.
    pub interfaces: Vec<ArticulatedInterfaceHeatSource>,
    /// Other heat only; friction heat and bound conductances have this owner.
    pub transport: StepCoefficients<P>,
}

/// Independent native mechanics, objective transport and finite-cell budgets.
#[derive(Clone, Copy, Debug)]
pub struct ArticulatedFrictionHeatControls {
    pub mechanical: ArticulatedMultiTangentialConfig,
    pub transport: TangentialTransportConfig,
    pub maximum_cells: usize,
    /// Bound cells plus distinct common interface temperatures before allocation.
    pub maximum_unknowns: usize,
    pub maximum_trials: usize,
    pub maximum_newton_iterations: usize,
    pub difference_step_j: f64,
    /// Numerical coordinate scale [J/K], not an added physical heat capacity.
    pub temperature_scale_j_per_k: f64,
    pub heat: FrictionHeatTolerances,
    pub combined_energy_tolerance_j: f64,
    pub newton: NewtonKrylovConfig,
    /// Final-candidate caloric entropy audit; all active materials must support it.
    pub audit_entropy: bool,
}

#[derive(Debug)]
pub enum ArticulatedFrictionHeatError<E> {
    Invalid(&'static str),
    Mechanical(ArticulatedContactError),
    Thermal(DepositionError),
    Source(E),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    TrialBudget,
}
impl<E> From<ArticulatedContactError> for ArticulatedFrictionHeatError<E> {
    fn from(error: ArticulatedContactError) -> Self {
        Self::Mechanical(error)
    }
}
impl<E> From<DepositionError> for ArticulatedFrictionHeatError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl<E> From<NewtonError> for ArticulatedFrictionHeatError<E> {
    fn from(error: NewtonError) -> Self {
        Self::Newton(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for ArticulatedFrictionHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "articulated friction/heat: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for ArticulatedFrictionHeatError<E> {}

/// Complete private candidate. The surrounding integrator must publish every
/// native body response, framed contact history and thermal endpoint together.
pub struct HeatedArticulatedFrictionCandidate<'m, S, P> {
    /// Original mechanical order; the matrices are owned only in this final report.
    pub mechanical: ArticulatedFramedMultiTangentialResult<S>,
    pub thermal: DepositionSimulation<'m>,
    /// Interfaces follow ascending original thermal link IDs.
    pub heat: FrictionHeatNetworkReport,
    /// Mechanical index corresponding to each canonical heat interface.
    pub contact_indices: Vec<usize>,
    pub source: P,
    pub nonlinear: NewtonReport,
    pub trials: usize,
    pub target_time_s: f64,
    /// Actual kinetic + accepted shear storage + caloric changes, less support
    /// work and outside heat. Physical friction cancels as one internal transfer.
    pub combined_energy_residual_j: f64,
}

struct Staged<'p, S, P> {
    mechanical: ArticulatedFramedMultiTangentialResult<S, &'p [f64]>,
    source: P,
}

/// Propose simultaneous articulated friction and finite-cell heat atomically.
///
/// Both accepted owners are borrowed. Mechanical sites retain their order;
/// bindings may be arbitrarily ordered, but must cover each site exactly once
/// with a distinct original link. The callback and heat report use ascending
/// link IDs. Cell-energy guesses follow sorted original CellId order.
///
/// All drive temperatures must be absent: they are solved coordinates. Every
/// trial calls the same original law with the same once-transported native
/// history. The callback resolves conductivities and other powers, never law
/// substitutions. It must stage source receipts without external publication
/// and must not add friction heat a second time. A source law owns applicability
/// and coefficient continuation; a frozen elastic source must refuse unsupported
/// temperatures. Variable elastic free energy needs a separate reversible
/// heat/entropy owner, which this adapter does not invent.
///
/// Contact points must coincide with the thermal shared face centre, and the
/// physical tangent normal must be parallel to both original cell centrelines.
/// Route A/B ownership is explicit, including prescribed supports with finite
/// thermal cells. The caller owns body/material identity and any supplied
/// nominal pressure. Poses, normal loads, faces and support velocities remain
/// frozen. Transport-then-shear mechanics and backward-Euler heat require time
/// and spatial refinement; nonlinear convergence does not establish uniqueness.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn propose_articulated_friction_heat<'m, L: ObjectiveTangentialLaw, P, E>(
    systems: &[ArticulatedContactSystem<'_>],
    sites: &[ArticulatedFramedTangentialSite<'_, L>],
    thermal: &DepositionSimulation<'m>,
    accepted_time_s: f64,
    bindings: &[ArticulatedFrictionHeatRoute<'_>],
    controls: ArticulatedFrictionHeatControls,
    energy_guess_j: &[f64],
    resolve: impl Fn(
        &DepositionSimulation<'m>,
        &[u64],
        &[f64],
    ) -> Result<ArticulatedFrictionHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<HeatedArticulatedFrictionCandidate<'m, L::State, P>, ArticulatedFrictionHeatError<E>> {
    let checkpoint = || {
        if gate.is_requested() {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    };
    checkpoint()?;
    controls.newton.validate()?;
    let contacts = sites.len();
    let cells = thermal.cells().len();
    let dimension = cells
        .checked_add(contacts)
        .ok_or(ArticulatedFrictionHeatError::Invalid(
            "thermal coordinate count overflow",
        ))?;
    let dt_s = sites.first().map_or(f64::NAN, |site| site.drive.dt_s);
    let target_time_s = accepted_time_s + dt_s;
    if !accepted_time_s.is_finite()
        || accepted_time_s != thermal.time_s()
        || !dt_s.is_finite()
        || dt_s <= 0.0
        || !target_time_s.is_finite()
        || target_time_s <= accepted_time_s
        || contacts == 0
        || contacts > controls.mechanical.coupling.maximum_contacts
        || bindings.len() != contacts
        || controls.maximum_cells == 0
        || cells > controls.maximum_cells
        || dimension > controls.maximum_unknowns
        || controls.maximum_trials == 0
        || controls.maximum_newton_iterations == 0
        || !controls.difference_step_j.is_finite()
        || controls.difference_step_j <= 0.0
        || !controls.temperature_scale_j_per_k.is_finite()
        || controls.temperature_scale_j_per_k <= 0.0
        || [
            controls.heat.cell_j,
            controls.heat.thermal_j,
            controls.heat.interface_j,
            controls.heat.total_j,
            controls.combined_energy_tolerance_j,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
        || energy_guess_j.len() != cells
        || energy_guess_j.iter().any(|x| !x.is_finite())
        || sites
            .iter()
            .any(|site| site.drive.contact_temperature_k.is_some() || site.drive.dt_s != dt_s)
    {
        return Err(ArticulatedFrictionHeatError::Invalid(
            "synchronized bounded cells/sites and solved contact temperatures required",
        ));
    }
    let mut seen = vec![false; contacts];
    for binding in bindings {
        checkpoint()?;
        if binding.contact_index >= contacts
            || seen[binding.contact_index]
            || !binding.temperature_guess_k.is_finite()
            || binding.temperature_guess_k <= 0.0
        {
            return Err(ArticulatedFrictionHeatError::Invalid(
                "each mechanical contact requires one heat route and positive Kelvin guess",
            ));
        }
        seen[binding.contact_index] = true;
    }
    let mut ordered: Vec<_> = bindings.iter().collect();
    ordered.sort_by_key(|binding| binding.route.link_id());
    if ordered
        .windows(2)
        .any(|pair| pair[0].route.link_id() == pair[1].route.link_id())
    {
        return Err(ArticulatedFrictionHeatError::Invalid(
            "each original thermal link requires exactly one mechanical heat owner",
        ));
    }
    let routes: Vec<_> = ordered.iter().map(|binding| binding.route).collect();
    let contact_indices: Vec<_> = ordered
        .iter()
        .map(|binding| binding.contact_index)
        .collect();
    let step = FrictionHeatNetworkStep::new(
        thermal,
        target_time_s,
        &routes,
        controls.mechanical.coupling.maximum_contacts,
        checkpoint,
    )?;
    let prepared = PreparedArticulatedFramedTangentialContacts::new(
        systems,
        sites,
        controls.mechanical,
        controls.transport,
        || gate.is_requested(),
    )?;
    for binding in &ordered {
        checkpoint()?;
        let point = prepared
            .mechanics()
            .contact_world_point(binding.contact_index)
            .ok_or(ArticulatedFrictionHeatError::Invalid(
                "missing prepared contact point",
            ))?;
        validate_geometry(
            [point.x, point.y, point.z],
            sites[binding.contact_index].contact.frame_world.normal(),
            binding.route,
            controls.mechanical.mechanical.geometry_tolerance_m,
        )?;
    }
    let trials = Cell::new(0usize);
    let evaluate = |x: &[f64], with_entropy: bool| {
        checkpoint()?;
        if trials.get() >= controls.maximum_trials {
            return Err(ArticulatedFrictionHeatError::TrialBudget);
        }
        trials.set(trials.get() + 1);
        if x.len() != dimension || x.iter().any(|x| !x.is_finite()) {
            return Err(ArticulatedFrictionHeatError::Invalid(
                "complete finite articulated/thermal iterate required",
            ));
        }
        let temperatures: Vec<_> = ordered
            .iter()
            .enumerate()
            .map(|(i, binding)| {
                binding.temperature_guess_k + x[cells + i] / controls.temperature_scale_j_per_k
            })
            .collect();
        let endpoint = |thermal: &DepositionSimulation<'m>, temperatures: &[f64]| {
            let ArticulatedFrictionHeatSource {
                interfaces,
                transport,
            } = resolve(thermal, step.link_ids(), temperatures)
                .map_err(ArticulatedFrictionHeatError::Source)?;
            if interfaces.len() != contacts
                || interfaces
                    .iter()
                    .zip(step.link_ids())
                    .any(|(source, id)| source.link_id != *id)
            {
                return Err(ArticulatedFrictionHeatError::Invalid(
                    "conductivities must retain every canonical original link identity",
                ));
            }
            let mut inputs: Vec<_> = sites
                .iter()
                .map(|site| PreparedArticulatedFramedInput {
                    law: site.law,
                    drive: site.drive,
                })
                .collect();
            for (i, &contact) in contact_indices.iter().enumerate() {
                inputs[contact].drive.contact_temperature_k = Some(temperatures[i]);
            }
            let mechanical = prepared.advance(&inputs, || gate.is_requested())?;
            let interfaces = interfaces
                .iter()
                .enumerate()
                .map(|(i, source)| InterfaceHeatInput {
                    dissipated_j: mechanical.mechanics.contacts[contact_indices[i]].dissipated_j,
                    conductivity_w_m_k: source.conductivity_w_m_k,
                })
                .collect();
            let StepCoefficients {
                sources_w,
                conductances_w_k,
                payload,
            } = transport;
            Ok(FrictionHeatNetworkInput {
                interfaces,
                transport: StepCoefficients {
                    sources_w,
                    conductances_w_k,
                    payload: Staged {
                        mechanical,
                        source: payload,
                    },
                },
            })
        };
        if with_entropy {
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
        return Err(ArticulatedFrictionHeatError::Interrupted(result.progress));
    }
    if !result.report.converged {
        return Err(ArticulatedFrictionHeatError::NotConverged(result.report));
    }
    let trial = evaluate(&solver.x, controls.audit_entropy)?;
    let combined = trial
        .payload()
        .mechanical
        .accepted_energy_balance_residual_j
        + trial.report().total_energy_residual_j;
    if !combined.is_finite() || combined.abs() > controls.combined_energy_tolerance_j {
        return Err(ArticulatedFrictionHeatError::Invalid(
            "joint articulated/contact/caloric energy does not balance",
        ));
    }
    let (thermal, staged, heat) = trial.accept(controls.heat, checkpoint)?;
    checkpoint()?;
    drop(problem);
    let mechanical = staged.mechanical.into_owned();
    checkpoint()?;
    Ok(HeatedArticulatedFrictionCandidate {
        mechanical,
        thermal,
        heat,
        contact_indices,
        source: staged.source,
        nonlinear: result.report,
        trials: trials.get(),
        target_time_s,
        combined_energy_residual_j: combined,
    })
}

fn validate_geometry<E>(
    point: [f64; 3],
    normal: [f64; 3],
    route: &FrictionHeatRoute,
    tolerance_m: f64,
) -> Result<(), ArticulatedFrictionHeatError<E>> {
    let geometry = route.geometry;
    let centre = [geometry[0].start_m, geometry[0].end_m]
        .into_iter()
        .find(|p| *p == geometry[1].start_m || *p == geometry[1].end_m)
        .ok_or(ArticulatedFrictionHeatError::Invalid(
            "thermal route has no common face centre",
        ))?;
    if norm(std::array::from_fn(|i| point[i] - centre[i])) > tolerance_m {
        return Err(ArticulatedFrictionHeatError::Invalid(
            "articulated contact must coincide with its thermal face centre",
        ));
    }
    for g in geometry {
        let a = std::array::from_fn::<_, 3, _>(|i| g.end_m[i] - g.start_m[i]);
        let cross = [
            a[1] * normal[2] - a[2] * normal[1],
            a[2] * normal[0] - a[0] * normal[2],
            a[0] * normal[1] - a[1] * normal[0],
        ];
        if norm(cross) > tolerance_m {
            return Err(ArticulatedFrictionHeatError::Invalid(
                "articulated tangent normal and thermal cell centrelines differ",
            ));
        }
    }
    Ok(())
}
fn norm(v: [f64; 3]) -> f64 {
    det::hypot(det::hypot(v[0], v[1]), v[2])
}
fn numerical<E>(
    error: FiniteDifferenceError<ArticulatedFrictionHeatError<E>>,
) -> ArticulatedFrictionHeatError<E> {
    match error {
        FiniteDifferenceError::Source(source) => source,
        FiniteDifferenceError::Invalid(message) => ArticulatedFrictionHeatError::Invalid(message),
    }
}
