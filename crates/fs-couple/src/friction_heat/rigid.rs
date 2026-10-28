//! Retained finite-patch rigid friction and the original thermal cells accept together.
//!
//! Every temperature trial starts from accepted body momenta and regional slip
//! memory. The native rigid owner derives all shear/spin mobility from masses,
//! inertias, poses and contact arms. Its physical plastic heat alone enters the
//! finite enthalpy network. Reversible storage and both numerical losses retain
//! their own accounts; no effective mass, heat fraction or thermal store is added.

/// Simultaneous native rigid contact islands with shared thermal cells.
pub mod network;

use super::{
    FrictionHeatInput, FrictionHeatReport, FrictionHeatRoute, FrictionHeatStep,
    FrictionHeatTolerances,
};
use fs_contact::tangential::nonsmooth::{NonsmoothContact, NonsmoothTangentialState};
use fs_exec::CancelGate;
use fs_mbd::tangential_friction::{
    RigidContactBinding, RigidFrictionControls, RigidFrictionError, RigidFrictionReport,
    apply_rigid_friction,
};
use fs_mbd::{MassProperties, RigidBodyState, Vec3};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};
use fs_tribo::partial_slip::PartialSlipLaw;
use std::cell::Cell;

/// Original source evaluated at the same candidate thermal endpoint.
pub struct RigidFrictionHeatSource<P> {
    /// Original retained law with temperature-resolved Coulomb coefficients.
    /// Stiffness, regional decomposition, source and model identity stay fixed.
    pub law: PartialSlipLaw,
    /// Conductivity from each original contacting material at its own cell state.
    pub conductivity_w_m_k: [f64; 2],
    /// Other heat/conduction only; physical friction work is supplied by mechanics.
    pub transport: StepCoefficients<P>,
}

/// Independent mechanical, nonlinear and thermal acceptance controls.
#[derive(Debug, Clone, Copy)]
pub struct RigidFrictionHeatControls {
    pub mechanical: RigidFrictionControls,
    pub maximum_cells: usize,
    pub maximum_trials: usize,
    pub maximum_newton_iterations: usize,
    pub difference_step_j: f64,
    /// Newton coordinate scale [J/K], never an additional heat capacity.
    pub temperature_scale_j_per_k: f64,
    pub heat: FrictionHeatTolerances,
    pub combined_energy_tolerance_j: f64,
    /// The inner forcing floor must respect finite-difference accuracy and the
    /// original mechanical solve tolerances. All physical gates remain separate.
    pub newton: NewtonKrylovConfig,
}

#[derive(Debug)]
pub enum RigidFrictionHeatError<E> {
    Invalid(&'static str),
    Mechanical(RigidFrictionError),
    Thermal(DepositionError),
    Source(E),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    TrialBudget,
    /// This consumer has no physical coefficient-work supply port. A changed
    /// storage law or any nonzero parameter work therefore refuses, even if a
    /// loose numerical energy tolerance would otherwise conceal it.
    UnbalancedCoefficientWork,
}
impl<E> From<DepositionError> for RigidFrictionHeatError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl<E> From<RigidFrictionError> for RigidFrictionHeatError<E> {
    fn from(error: RigidFrictionError) -> Self {
        Self::Mechanical(error)
    }
}
impl<E> From<NewtonError> for RigidFrictionHeatError<E> {
    fn from(error: NewtonError) -> Self {
        Self::Newton(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for RigidFrictionHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rigid-body friction/heat: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for RigidFrictionHeatError<E> {}

/// Joint receipt on the exact original rigid and finite-thermal endpoints.
pub struct HeatedRigidFrictionReport<P> {
    pub mechanical: RigidFrictionReport,
    pub heat: FrictionHeatReport,
    pub law: PartialSlipLaw,
    pub source: P,
    pub nonlinear: NewtonReport,
    pub trials: usize,
    /// Kinetic + retained elastic + caloric changes, minus outside heat, plus
    /// separately reported elastic/mechanical integration losses. Physical
    /// friction cancels internally once; no coefficient-work port is accepted.
    pub combined_energy_residual_j: f64,
}

struct Staged<P> {
    bodies: Vec<RigidBodyState>,
    history: NonsmoothTangentialState,
    mechanical: RigidFrictionReport,
    law: PartialSlipLaw,
    source: P,
}

/// Solve one loaded shear/spin patch and its finite thermal network atomically.
///
/// `route` A/B belongs to `binding` body A/B. Both original point arms must
/// meet the thermal route's shared face centre, and its area must equal the
/// supplied normal ellipse. The thermal cell coordinates are world coordinates
/// for this frozen-geometry phase. Poses, normal admission and any subsequent
/// material transport remain the surrounding mechanical owner's responsibility.
///
/// Only the original law's static/kinetic Coulomb coefficients may vary with
/// temperature. Retained history uses that law's coefficient continuation;
/// changing stiffness, decomposition or source is refused. There is no implicit
/// source for parameter work and no conversion of numerical damping into heat.
///
/// `resolve` must stage its receipt without publishing state, and must not add
/// the friction heat to its other powers. Energy guesses follow sorted CellId
/// order; temperature is a Newton guess, never an imposed boundary condition.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn apply_rigid_friction_heat<'m, P, E>(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    history: &mut NonsmoothTangentialState,
    contact: &NonsmoothContact<'_>,
    binding: RigidContactBinding,
    thermal: &mut DepositionSimulation<'m>,
    route: &FrictionHeatRoute,
    dt_s: f64,
    controls: RigidFrictionHeatControls,
    energy_guess_j: &[f64],
    temperature_guess_k: f64,
    resolve: impl Fn(&DepositionSimulation<'m>, f64) -> Result<RigidFrictionHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<HeatedRigidFrictionReport<P>, RigidFrictionHeatError<E>> {
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
        || binding.body_a >= bodies.len()
        || binding.body_b >= bodies.len()
        || binding.body_a == binding.body_b
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
        || !temperature_guess_k.is_finite()
        || temperature_guess_k <= 0.0
        || energy_guess_j.len() != thermal.cells().len()
        || energy_guess_j.iter().any(|x| !x.is_finite())
    {
        return Err(RigidFrictionHeatError::Invalid(
            "synchronized clocks and finite bounded rigid/thermal coordinates required",
        ));
    }
    validate_geometry(
        bodies,
        properties,
        contact,
        binding,
        route,
        controls.mechanical.geometry_tolerance_m,
    )?;
    if let Some(old) = history.contact_checkpoint(contact.patch.patch_id()) {
        compatible_law(old.law(), contact.law)?;
    }
    let step = FrictionHeatStep::new(thermal, target, route, checkpoint)?;
    let cells = step.cell_ids().len();
    let trials = Cell::new(0usize);
    let evaluate = |x: &[f64]| {
        checkpoint()?;
        if trials.get() >= controls.maximum_trials {
            return Err(RigidFrictionHeatError::TrialBudget);
        }
        trials.set(trials.get() + 1);
        if x.len() != cells + 1 {
            return Err(RigidFrictionHeatError::Invalid(
                "incomplete rigid/thermal iterate",
            ));
        }
        let temperature = temperature_guess_k + x[cells] / controls.temperature_scale_j_per_k;
        step.evaluate(
            &x[..cells],
            temperature,
            |endpoint, interface_temperature| {
                let RigidFrictionHeatSource {
                    law,
                    conductivity_w_m_k,
                    transport,
                } = resolve(endpoint, interface_temperature)
                    .map_err(RigidFrictionHeatError::Source)?;
                compatible_law(contact.law, &law)?;
                let StepCoefficients {
                    sources_w,
                    conductances_w_k,
                    payload,
                } = transport;
                let mut next_bodies = bodies.to_vec();
                let mut next_history = history.clone();
                let resolved = NonsmoothContact {
                    law: &law,
                    patch: contact.patch,
                    interface: contact.interface,
                    frame: contact.frame,
                    common_spin_rad_per_s: contact.common_spin_rad_per_s,
                    torsional_length_m: contact.torsional_length_m,
                    creepage_reference_speed_mps: contact.creepage_reference_speed_mps,
                    work: contact.work.clone(),
                };
                let mechanical = apply_rigid_friction(
                    &mut next_bodies,
                    properties,
                    &mut next_history,
                    &[resolved],
                    &[binding],
                    dt_s,
                    controls.mechanical,
                    gate,
                )?;
                if mechanical.friction.coefficient_work_j != 0.0 {
                    return Err(RigidFrictionHeatError::UnbalancedCoefficientWork);
                }
                Ok(FrictionHeatInput {
                    dissipated_j: mechanical.friction.physical_heat_j,
                    conductivity_w_m_k,
                    transport: StepCoefficients {
                        sources_w,
                        conductances_w_k,
                        payload: Staged {
                            bodies: next_bodies,
                            history: next_history,
                            mechanical,
                            law,
                            source: payload,
                        },
                    },
                })
            },
            checkpoint,
        )
    };
    let problem =
        FallibleFiniteDifference::new(cells + 1, controls.difference_step_j, |x: &[f64]| {
            evaluate(x).map(|trial| trial.residual_j().to_vec())
        })
        .map_err(numerical)?;
    let mut initial = energy_guess_j.to_vec();
    initial.push(0.0);
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
            "joint rigid/thermal work does not balance",
        ));
    }
    let (next_thermal, staged, heat) = trial.accept(controls.heat, checkpoint)?;
    checkpoint()?;
    drop(problem);
    // No fallible operation or user callback follows this publication boundary.
    bodies.copy_from_slice(&staged.bodies);
    *history = staged.history;
    *thermal = next_thermal;
    Ok(HeatedRigidFrictionReport {
        mechanical: staged.mechanical,
        heat,
        law: staged.law,
        source: staged.source,
        nonlinear: result.report,
        trials: trials.get(),
        combined_energy_residual_j: combined,
    })
}

fn compatible_law<E>(
    accepted: &PartialSlipLaw,
    candidate: &PartialSlipLaw,
) -> Result<(), RigidFrictionHeatError<E>> {
    let old = accepted.parameters();
    let mut next = candidate.parameters();
    if old.tangential_stiffness_n_per_m != next.tangential_stiffness_n_per_m
        || old.torsional_stiffness_nm_per_rad != next.torsional_stiffness_nm_per_rad
    {
        return Err(RigidFrictionHeatError::UnbalancedCoefficientWork);
    }
    next.static_mu = old.static_mu;
    next.kinetic_mu = old.kinetic_mu;
    if old != next
        || accepted.model_id() != candidate.model_id()
        || accepted.source_id() != candidate.source_id()
    {
        return Err(RigidFrictionHeatError::Invalid(
            "thermal continuation must retain the original source and regional storage law",
        ));
    }
    Ok(())
}

fn validate_geometry<E>(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contact: &NonsmoothContact<'_>,
    binding: RigidContactBinding,
    route: &FrictionHeatRoute,
    tolerance_m: f64,
) -> Result<(), RigidFrictionHeatError<E>> {
    if !tolerance_m.is_finite() || tolerance_m <= 0.0 {
        return Err(RigidFrictionHeatError::Invalid(
            "positive geometry tolerance required",
        ));
    }
    let area = core::f64::consts::PI
        * contact.patch.semi_axis_longitudinal_m()
        * contact.patch.semi_axis_lateral_m();
    if !area.is_finite() || (area - route.patch.link.area_m2).abs() > 128.0 * f64::EPSILON * area {
        return Err(RigidFrictionHeatError::Invalid(
            "normal ellipse and thermal contact must have the same physical area",
        ));
    }
    let geometry = route.geometry;
    let centre = [geometry[0].start_m, geometry[0].end_m]
        .into_iter()
        .find(|point| *point == geometry[1].start_m || *point == geometry[1].end_m)
        .ok_or(RigidFrictionHeatError::Invalid(
            "thermal route has no common face centre",
        ))?;
    let normal = Vec3::new(
        contact.frame.normal()[0],
        contact.frame.normal()[1],
        contact.frame.normal()[2],
    );
    let axis = Vec3::new(
        geometry[0].end_m[0] - geometry[0].start_m[0],
        geometry[0].end_m[1] - geometry[0].start_m[1],
        geometry[0].end_m[2] - geometry[0].start_m[2],
    );
    if norm(axis.cross(normal)) > tolerance_m {
        return Err(RigidFrictionHeatError::Invalid(
            "rigid normal and thermal cell centreline differ",
        ));
    }
    let centre = Vec3::new(centre[0], centre[1], centre[2]);
    for (body, arm) in [
        (binding.body_a, binding.arm_a_body),
        (binding.body_b, binding.arm_b_body),
    ] {
        let point = bodies[body]
            .point_kinematics(properties[body], arm)
            .map_err(RigidFrictionError::from)?
            .point_world;
        if norm(point.sub(centre)) > tolerance_m {
            return Err(RigidFrictionHeatError::Invalid(
                "rigid contact must coincide with the thermal face centre",
            ));
        }
    }
    Ok(())
}
fn norm(v: Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}
fn numerical<E>(
    error: FiniteDifferenceError<RigidFrictionHeatError<E>>,
) -> RigidFrictionHeatError<E> {
    match error {
        FiniteDifferenceError::Source(source) => source,
        FiniteDifferenceError::Invalid(message) => RigidFrictionHeatError::Invalid(message),
    }
}
