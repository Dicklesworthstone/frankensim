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
    PreparedRigidHertzPair, RigidHertzCandidate, RigidHertzContact, RigidHertzControls,
    RigidHertzError, RigidHertzReport,
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

/// The original ordered interface source and caller-owned closed normal patch.
#[derive(Clone, Copy)]
pub struct RigidHertzHeatContact<'a> {
    pub binding: RigidContactBinding,
    pub frame: TangentFrame,
    pub source: &'a SourceThermalCoulomb,
    pub normal_force_n: f64,
    pub contact_radius_m: f64,
    /// Explicit original source coordinate, never inferred from N or patch area.
    pub nominal_pressure_pa: Option<f64>,
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
pub enum RigidHertzHeatError<E> {
    Invalid(&'static str),
    Mechanical(RigidHertzError),
    FrictionSource(SourceHertzWrenchError),
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
impl<E> From<SourceHertzWrenchError> for RigidHertzHeatError<E> {
    fn from(e: SourceHertzWrenchError) -> Self {
        Self::FrictionSource(e)
    }
}
impl<E> From<DepositionError> for RigidHertzHeatError<E> {
    fn from(e: DepositionError) -> Self {
        Self::Thermal(e)
    }
}
impl<E> From<NewtonError> for RigidHertzHeatError<E> {
    fn from(e: NewtonError) -> Self {
        Self::Newton(e)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for RigidHertzHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Hertz friction/heat: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for RigidHertzHeatError<E> {}

/// Publish the two original body indices and thermal endpoint together.
/// Accepted bodies and the thermal owner remain borrowed throughout this solve.
pub struct HeatedRigidHertzCandidate<'m, P> {
    pub binding: RigidContactBinding,
    /// Original binding's A/B order; every other body is unaffected.
    pub bodies_after: [RigidBodyState; 2],
    pub mechanical: RigidHertzReport,
    pub thermal: DepositionSimulation<'m>,
    pub heat: FrictionHeatNetworkReport,
    /// The unchanged source owner and actual endpoint coefficient-use receipts.
    pub friction_source: SourceHertzWrench,
    pub source: P,
    pub nonlinear: NewtonReport,
    pub trials: usize,
    pub target_time_s: f64,
    /// Actual kinetic plus caloric change, minus outside heat.
    pub combined_energy_residual_j: f64,
}

struct Staged<P> {
    mechanical: RigidHertzCandidate,
    friction_source: SourceHertzWrench,
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
    let initial_source = contact.source.hertz_at(
        temperature_guess_k,
        contact.nominal_pressure_pa,
        contact.normal_force_n,
        contact.contact_radius_m,
    )?;
    let prepared = PreparedRigidHertzPair::new(
        bodies,
        properties,
        RigidHertzContact {
            binding: contact.binding,
            frame: contact.frame,
            law: initial_source.law(),
        },
        duration_s,
        controls.mechanical,
        gate,
    )?;
    let point = prepared.contact_points_world()[0];
    validate_geometry(
        [point.x, point.y, point.z],
        contact,
        route,
        controls.mechanical.geometry_tolerance_m,
    )?;
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
            let friction_source = contact.source.hertz_at(
                temperature,
                contact.nominal_pressure_pa,
                contact.normal_force_n,
                contact.contact_radius_m,
            )?;
            let mechanical = prepared.advance(friction_source.law(), gate)?;
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
                    dissipated_j: mechanical.report.dissipated_energy_j,
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
    let combined = trial.payload().mechanical.report.total_energy_residual_j
        + trial.report().total_energy_residual_j;
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

fn validate_geometry<E>(
    point: [f64; 3],
    contact: RigidHertzHeatContact<'_>,
    route: &FrictionHeatRoute,
    tolerance: f64,
) -> Result<(), RigidHertzHeatError<E>> {
    let geometry = route.geometry;
    let centre = [geometry[0].start_m, geometry[0].end_m]
        .into_iter()
        .find(|p| *p == geometry[1].start_m || *p == geometry[1].end_m)
        .ok_or(RigidHertzHeatError::Invalid(
            "thermal route has no shared face",
        ))?;
    let radius = det::sqrt(route.contact_area_m2() / core::f64::consts::PI);
    if norm(std::array::from_fn(|i| point[i] - centre[i])) > tolerance
        || !radius.is_finite()
        || (radius - contact.contact_radius_m).abs() > tolerance
    {
        return Err(RigidHertzHeatError::Invalid(
            "original Hertz centre/radius and thermal circular face must agree",
        ));
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
            return Err(RigidHertzHeatError::Invalid(
                "Hertz normal and finite-cell centrelines differ",
            ));
        }
    }
    Ok(())
}
fn norm(v: [f64; 3]) -> f64 {
    det::hypot(det::hypot(v[0], v[1]), v[2])
}
fn numerical<E>(error: FiniteDifferenceError<RigidHertzHeatError<E>>) -> RigidHertzHeatError<E> {
    match error {
        FiniteDifferenceError::Source(e) => e,
        FiniteDifferenceError::Invalid(s) => RigidHertzHeatError::Invalid(s),
    }
}
