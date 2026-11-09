//! Original Hertz force/spin friction, finite recoil, and contact heat together.
//!
//! One prepared native pair supplies the full three-coordinate body response.
//! Every thermal trial refreshes the same pinned interface coefficient claims
//! and solves that original Hertz law. Its actual local dissipation enters one
//! existing finite thermal network through a massless common-temperature face.
//! Only complete mechanics/source/thermal candidates leave this owner.

use super::network::{
    FrictionHeatNetworkInput, FrictionHeatNetworkReport, FrictionHeatNetworkStep,
    InterfaceHeatInput,
};
use super::{FrictionHeatRoute, FrictionHeatTolerances};
use fs_contact::interface_binding::thermal_friction::SourceThermalCoulomb;
use fs_contact::interface_binding::thermal_friction::hertz::{
    SourceHertzWrench, SourceHertzWrenchError,
};
use fs_exec::CancelGate;
use fs_math::det;
use fs_mbd::tangential_friction::RigidContactBinding;
use fs_mbd::tangential_friction::hertz::{
    PreparedRigidHertzPair, RigidHertzContact, RigidHertzControls, RigidHertzError,
    RigidHertzReport,
};
use fs_mbd::{MassProperties, RigidBodyState};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};
use fs_tribo::partial_slip::TangentFrame;
use std::cell::Cell;

/// Explicit continuous force/spin drive, accepted mode and finite-cell heat.
pub mod driven;
/// Time-resolved thermal feedback through accepted original Hertz phases.
pub mod march;
/// Source-resolved simultaneous Hertz patches on one finite body/heat island.
pub mod multiple;
/// Original local-speed friction curves with finite-cell thermal feedback.
pub mod rate_dependent;
/// Original finite-body stopping trajectories and their finite-cell heat.
pub mod trajectory;
/// Curved original Hertz trajectories with finite-body thermal feedback.
pub mod transient;

/// The original ordered interface source and caller-owned closed normal patch.
pub struct RigidHertzHeatContact<'a, S = SourceThermalCoulomb> {
    pub binding: RigidContactBinding,
    pub frame: TangentFrame,
    pub source: &'a S,
    pub normal_force_n: f64,
    pub contact_radius_m: f64,
    /// Explicit original source coordinate, never inferred from N or patch area.
    pub nominal_pressure_pa: Option<f64>,
}
impl<S> Copy for RigidHertzHeatContact<'_, S> {}
impl<S> Clone for RigidHertzHeatContact<'_, S> {
    fn clone(&self) -> Self {
        *self
    }
}

/// Endpoint bulk conductivities and other heat from the original material owner.
pub struct RigidHertzHeatSource<P> {
    pub conductivity_w_m_k: [f64; 2],
    /// Bound-contact heat/conductance must not be added here a second time.
    pub transport: StepCoefficients<P>,
}

#[derive(Clone, Copy, Debug)]
pub struct RigidHertzHeatControls {
    pub mechanical: RigidHertzControls,
    pub maximum_cells: usize,
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
pub enum RigidHertzHeatError<E, M = RigidHertzError, F = SourceHertzWrenchError> {
    Invalid(&'static str),
    Mechanical(M),
    FrictionSource(F),
    Thermal(DepositionError),
    Source(E),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    TrialBudget,
}
impl<E> From<RigidHertzError> for RigidHertzHeatError<E> {
    fn from(e: RigidHertzError) -> Self {
        Self::Mechanical(e)
    }
}
impl<E, M> From<SourceHertzWrenchError> for RigidHertzHeatError<E, M> {
    fn from(e: SourceHertzWrenchError) -> Self {
        Self::FrictionSource(e)
    }
}
impl<E, M, F> From<DepositionError> for RigidHertzHeatError<E, M, F> {
    fn from(e: DepositionError) -> Self {
        Self::Thermal(e)
    }
}
impl<E, M, F> From<NewtonError> for RigidHertzHeatError<E, M, F> {
    fn from(e: NewtonError) -> Self {
        Self::Newton(e)
    }
}
impl<E: core::fmt::Debug, M: core::fmt::Debug, F: core::fmt::Debug> core::fmt::Display
    for RigidHertzHeatError<E, M, F>
{
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Hertz friction/heat: {self:?}")
    }
}
impl<E: std::error::Error + 'static, M: core::fmt::Debug, F: core::fmt::Debug> std::error::Error
    for RigidHertzHeatError<E, M, F>
{
}

/// Publish the two original body indices and thermal endpoint together.
/// Accepted bodies and the thermal owner remain borrowed throughout this solve.
pub struct HeatedRigidHertzCandidate<'m, P, M = RigidHertzReport, F = SourceHertzWrench> {
    pub binding: RigidContactBinding,
    /// Original binding's A/B order; every other body is unaffected.
    pub bodies_after: [RigidBodyState; 2],
    pub mechanical: M,
    pub thermal: DepositionSimulation<'m>,
    pub heat: FrictionHeatNetworkReport,
    /// The unchanged source owner and actual endpoint coefficient-use receipts.
    pub friction_source: F,
    pub source: P,
    pub nonlinear: NewtonReport,
    pub trials: usize,
    pub target_time_s: f64,
    /// Actual kinetic plus caloric change, minus outside heat and mechanical work.
    pub combined_energy_residual_j: f64,
}

/// Native work already checked by its original mechanical owner. This bridge
/// carries neither an alternative friction law nor another physical state.
pub(super) struct HertzHeatMechanics<M> {
    pub bodies_after: [RigidBodyState; 2],
    pub report: M,
    pub dissipated_j: f64,
    pub total_energy_residual_j: f64,
}

struct Staged<P, M, F> {
    mechanical: HertzHeatMechanics<M>,
    friction_source: F,
    source: P,
}

/// Solve one closed Hertz patch and its original finite heat network.
///
/// The source is immutable. Only its original selected coefficient claims are
/// queried at the solved common face temperature and supplied pressure. Normal
/// load/radius, poses, arms, masses and inertia remain fixed. The same original
/// geometry must define the circular thermal face. Thermal A/B belongs to the
/// binding's A/B; the caller owns body/material identity and normal evolution.
///
/// This stateless static/gross-slide/spin rung cannot consume stored elastic
/// contact memory. Rolling loss is a different work owner. No heat fraction,
/// reservoir, velocity cutoff, coefficient fit or numerical-loss heating is
/// supplied. The held-effort mechanics retains its coarse-reversal refusal.
/// The finite-cell/backward-Euler thermal approximation requires spatial/time
/// refinement; convergence of the nonlinear root does not prove uniqueness.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn propose_rigid_hertz_heat<'m, P, E>(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contact: RigidHertzHeatContact<'_>,
    thermal: &DepositionSimulation<'m>,
    route: &FrictionHeatRoute,
    accepted_time_s: f64,
    dt_s: f64,
    controls: RigidHertzHeatControls,
    energy_guess_j: &[f64],
    temperature_guess_k: f64,
    resolve: impl Fn(&DepositionSimulation<'m>, f64) -> Result<RigidHertzHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<HeatedRigidHertzCandidate<'m, P>, RigidHertzHeatError<E>> {
    solve_rigid_hertz_heat(
        contact,
        thermal,
        route,
        accepted_time_s,
        dt_s,
        controls,
        energy_guess_j,
        temperature_guess_k,
        |temperature| refresh_source(contact, temperature),
        |source, duration| {
            let prepared = PreparedRigidHertzPair::new(
                bodies,
                properties,
                RigidHertzContact {
                    binding: contact.binding,
                    frame: contact.frame,
                    law: source.law(),
                },
                duration,
                controls.mechanical,
                gate,
            )?;
            let point = prepared.contact_points_world()[0];
            Ok((prepared, [point.x, point.y, point.z]))
        },
        |prepared, source| {
            let candidate = prepared.advance(source.law(), gate)?;
            Ok(HertzHeatMechanics {
                bodies_after: candidate.bodies_after,
                dissipated_j: candidate.report.dissipated_energy_j,
                total_energy_residual_j: candidate.report.total_energy_residual_j,
                report: candidate.report,
            })
        },
        resolve,
        gate,
    )
}

/// One source/geometry/nonlinear/thermal transaction for original native
/// mechanical owners. A trial always restarts the same accepted native pair.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(super) fn solve_rigid_hertz_heat<'m, P, E, C, S, M, ME, F, FE>(
    contact: RigidHertzHeatContact<'_, C>,
    thermal: &DepositionSimulation<'m>,
    route: &FrictionHeatRoute,
    accepted_time_s: f64,
    dt_s: f64,
    controls: RigidHertzHeatControls,
    energy_guess_j: &[f64],
    temperature_guess_k: f64,
    refresh: impl Fn(f64) -> Result<F, FE>,
    prepare: impl FnOnce(&F, f64) -> Result<(S, [f64; 3]), ME>,
    mechanics: impl Fn(&S, &F) -> Result<HertzHeatMechanics<M>, ME>,
    resolve: impl Fn(&DepositionSimulation<'m>, f64) -> Result<RigidHertzHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<HeatedRigidHertzCandidate<'m, P, M, F>, RigidHertzHeatError<E, ME, FE>> {
    let checkpoint = || {
        if gate.is_requested() {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    };
    checkpoint()?;
    controls.newton.validate()?;
    let target = accepted_time_s + dt_s;
    let cells = thermal.cells().len();
    let dimension = cells
        .checked_add(1)
        .ok_or(RigidHertzHeatError::Invalid("thermal dimension overflow"))?;
    if !accepted_time_s.is_finite()
        || accepted_time_s != thermal.time_s()
        || !dt_s.is_finite()
        || dt_s <= 0.0
        || !target.is_finite()
        || target <= accepted_time_s
        || controls.maximum_cells == 0
        || cells > controls.maximum_cells
        || controls.maximum_trials == 0
        || controls.maximum_newton_iterations == 0
        || !controls.difference_step_j.is_finite()
        || controls.difference_step_j <= 0.0
        || !controls.temperature_scale_j_per_k.is_finite()
        || controls.temperature_scale_j_per_k <= 0.0
        || !temperature_guess_k.is_finite()
        || temperature_guess_k <= 0.0
        || [
            controls.heat.cell_j,
            controls.heat.thermal_j,
            controls.heat.interface_j,
            controls.heat.total_j,
            controls.combined_energy_tolerance_j,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v < 0.0)
        || energy_guess_j.len() != cells
        || energy_guess_j.iter().any(|v| !v.is_finite())
    {
        return Err(RigidHertzHeatError::Invalid(
            "synchronized finite bounded mechanical/thermal inputs required",
        ));
    }
    // Both owners advance over the same representable clock interval. At a
    // nonzero accepted time this subtraction need not reproduce requested dt.
    let duration_s = target - accepted_time_s;
    let initial_source =
        refresh(temperature_guess_k).map_err(RigidHertzHeatError::FrictionSource)?;
    let (prepared, point) =
        prepare(&initial_source, duration_s).map_err(RigidHertzHeatError::Mechanical)?;
    validate_geometry(
        point,
        contact,
        route,
        controls.mechanical.geometry_tolerance_m,
    )
    .map_err(RigidHertzHeatError::Invalid)?;
    let step = FrictionHeatNetworkStep::new(thermal, target, &[route], 1, checkpoint)?;
    let trials = Cell::new(0usize);
    let evaluate = |x: &[f64], entropy: bool| {
        checkpoint()?;
        if trials.get() >= controls.maximum_trials {
            return Err(RigidHertzHeatError::TrialBudget);
        }
        trials.set(trials.get() + 1);
        if x.len() != dimension || x.iter().any(|v| !v.is_finite()) {
            return Err(RigidHertzHeatError::Invalid(
                "complete finite thermal iterate required",
            ));
        }
        let temperature = temperature_guess_k + x[cells] / controls.temperature_scale_j_per_k;
        let endpoint = |endpoint: &DepositionSimulation<'m>, temperatures: &[f64]| {
            let temperature = temperatures[0];
            let friction_source =
                refresh(temperature).map_err(RigidHertzHeatError::FrictionSource)?;
            let mechanical =
                mechanics(&prepared, &friction_source).map_err(RigidHertzHeatError::Mechanical)?;
            let RigidHertzHeatSource {
                conductivity_w_m_k,
                transport,
            } = resolve(endpoint, temperature).map_err(RigidHertzHeatError::Source)?;
            let StepCoefficients {
                sources_w,
                conductances_w_k,
                payload,
            } = transport;
            Ok(FrictionHeatNetworkInput {
                interfaces: vec![InterfaceHeatInput {
                    dissipated_j: mechanical.dissipated_j,
                    conductivity_w_m_k,
                }],
                transport: StepCoefficients {
                    sources_w,
                    conductances_w_k,
                    payload: Staged {
                        mechanical,
                        friction_source,
                        source: payload,
                    },
                },
            })
        };
        if entropy {
            step.evaluate_with_entropy(&x[..cells], &[temperature], endpoint, checkpoint)
        } else {
            step.evaluate(&x[..cells], &[temperature], endpoint, checkpoint)
        }
    };
    let problem =
        FallibleFiniteDifference::new(dimension, controls.difference_step_j, |x: &[f64]| {
            evaluate(x, false).map(|trial| trial.residual_j().to_vec())
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
        return Err(RigidHertzHeatError::Interrupted(result.progress));
    }
    if !result.report.converged {
        return Err(RigidHertzHeatError::NotConverged(result.report));
    }
    let trial = evaluate(&solver.x, controls.audit_entropy)?;
    let combined =
        trial.payload().mechanical.total_energy_residual_j + trial.report().total_energy_residual_j;
    if !combined.is_finite() || combined.abs() > controls.combined_energy_tolerance_j {
        return Err(RigidHertzHeatError::Invalid(
            "original body and finite-cell energy does not balance",
        ));
    }
    let (thermal, staged, heat) = trial.accept(controls.heat, checkpoint)?;
    checkpoint()?;
    Ok(HeatedRigidHertzCandidate {
        binding: contact.binding,
        bodies_after: staged.mechanical.bodies_after,
        mechanical: staged.mechanical.report,
        thermal,
        heat,
        friction_source: staged.friction_source,
        source: staged.source,
        nonlinear: result.report,
        trials: trials.get(),
        target_time_s: target,
        combined_energy_residual_j: combined,
    })
}

pub(super) fn refresh_source(
    contact: RigidHertzHeatContact<'_>,
    temperature_k: f64,
) -> Result<SourceHertzWrench, SourceHertzWrenchError> {
    contact.source.hertz_at(
        temperature_k,
        contact.nominal_pressure_pa,
        contact.normal_force_n,
        contact.contact_radius_m,
    )
}

pub(super) fn validate_geometry<S>(
    point: [f64; 3],
    contact: RigidHertzHeatContact<'_, S>,
    route: &FrictionHeatRoute,
    tolerance: f64,
) -> Result<(), &'static str> {
    let geometry = route.geometry;
    let centre = [geometry[0].start_m, geometry[0].end_m]
        .into_iter()
        .find(|p| *p == geometry[1].start_m || *p == geometry[1].end_m)
        .ok_or("thermal route has no shared face")?;
    let radius = det::sqrt(route.contact_area_m2() / core::f64::consts::PI);
    if norm(std::array::from_fn(|i| point[i] - centre[i])) > tolerance
        || !radius.is_finite()
        || (radius - contact.contact_radius_m).abs() > tolerance
    {
        return Err("original Hertz centre/radius and thermal circular face must agree");
    }
    let n = contact.frame.normal();
    for g in geometry {
        let a = std::array::from_fn::<_, 3, _>(|i| g.end_m[i] - g.start_m[i]);
        if norm([
            a[1] * n[2] - a[2] * n[1],
            a[2] * n[0] - a[0] * n[2],
            a[0] * n[1] - a[1] * n[0],
        ]) > tolerance
        {
            return Err("Hertz normal and finite-cell centrelines differ");
        }
    }
    Ok(())
}
fn norm(v: [f64; 3]) -> f64 {
    det::hypot(det::hypot(v[0], v[1]), v[2])
}
fn numerical<E, M, F>(
    error: FiniteDifferenceError<RigidHertzHeatError<E, M, F>>,
) -> RigidHertzHeatError<E, M, F> {
    match error {
        FiniteDifferenceError::Source(e) => e,
        FiniteDifferenceError::Invalid(s) => RigidHertzHeatError::Invalid(s),
    }
}
