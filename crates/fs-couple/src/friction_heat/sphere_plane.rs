//! A sourced spherical spin brake on the original finite thermal network.
//!
//! The native sphere owner resolves spin arrest, static retention, breakaway,
//! reversal and actual attitude. Its integrated physical Hertz loss enters the
//! original massless thermal face exactly once. Both finite body temperatures
//! and the common face temperature are solved together; that face temperature
//! refreshes the original pinned friction claims in every coupled trial.
//!
//! This stationary, circular, area-averaged contact-column model admits pure
//! normal-axis spin only. It does not transport heat through a moving footprint,
//! allocate collision loss, infer thermal expansion or resolve radial flash
//! temperature. Coefficients are held at each trial's endpoint temperature for
//! the whole native trajectory: thermal-window refinement remains necessary.

use super::network::{
    FrictionHeatNetworkInput, FrictionHeatNetworkReport, FrictionHeatNetworkStep,
    InterfaceHeatInput,
};
use super::{FrictionHeatRoute, FrictionHeatTolerances, finite};
use fs_contact::interface_binding::thermal_friction::SourceThermalCoulomb;
use fs_contact::interface_binding::thermal_friction::hertz::{
    SourceHertzWrench, SourceHertzWrenchError,
};
use fs_exec::CancelGate;
use fs_math::det;
use fs_mbd::Vec3;
use fs_mbd::sphere_plane::hertz_support::{
    SpherePlaneHertzState, SpherePlaneHertzSupportAdvance, SpherePlaneHertzSupportControls,
    SpherePlaneHertzSupportError, SpherePlaneHertzSupportPhase,
    advance_sphere_plane_hertz_support_with,
};
use fs_mbd::sphere_plane::{SpherePlane, SpherePlaneError, SpherePlaneEventKind};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};
use std::cell::Cell;

/// Original ordered sphere-A/plane-B interface and its separately admitted
/// circular footprint. The ideal normal support supplies the actual force.
#[derive(Clone, Copy)]
pub struct SpherePlaneHertzHeatContact<'a> {
    pub source: &'a SourceThermalCoulomb,
    pub contact_radius_m: f64,
    /// Original thermal body identities in interface A/B order. The surrounding
    /// material owner binds their local caloric regions to the interface card.
    pub thermal_body_ids: [u64; 2],
    /// Explicit original source coordinate, not N divided by patch area.
    pub nominal_pressure_pa: Option<f64>,
}

/// Each bulk's conductivity at its own actual candidate temperature, plus
/// other transport/heat only. Bound-contact conduction and friction are owned
/// by this composition and may not be added again by the callback.
pub struct SpherePlaneHertzHeatSource<P> {
    pub conductivity_w_m_k: [f64; 2],
    pub transport: StepCoefficients<P>,
}

#[derive(Clone, Copy, Debug)]
pub struct SpherePlaneHertzHeatControls {
    pub mechanical: SpherePlaneHertzSupportControls,
    pub maximum_cells: usize,
    pub maximum_trials: usize,
    pub maximum_newton_iterations: usize,
    pub difference_step_j: f64,
    /// Coordinate scaling for the nonlinear solve, never a physical capacity.
    pub temperature_scale_j_per_k: f64,
    pub geometry_tolerance_m: f64,
    pub heat: FrictionHeatTolerances,
    pub combined_energy_tolerance_j: f64,
    pub newton: NewtonKrylovConfig,
    pub audit_entropy: bool,
}

#[derive(Debug)]
pub enum SpherePlaneHertzHeatPhaseError {
    Invalid(&'static str),
    Friction(SourceHertzWrenchError),
}

#[derive(Debug)]
pub enum SpherePlaneHertzHeatError<E> {
    Invalid(&'static str),
    Mechanical(SpherePlaneHertzSupportError<SpherePlaneHertzHeatPhaseError>),
    Thermal(DepositionError),
    Source(E),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    TrialBudget,
}
impl<E> From<DepositionError> for SpherePlaneHertzHeatError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl<E> From<NewtonError> for SpherePlaneHertzHeatError<E> {
    fn from(error: NewtonError) -> Self {
        Self::Newton(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for SpherePlaneHertzHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spherical Hertz spin/finite heat: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for SpherePlaneHertzHeatError<E> {}

/// Original native motion, source receipt and finite thermal endpoint, accepted
/// together. No caller-owned state changes while any trial may still refuse.
pub struct HeatedSpherePlaneHertzCandidate<'m, P> {
    pub mechanical: SpherePlaneHertzSupportAdvance,
    pub thermal: DepositionSimulation<'m>,
    pub heat: FrictionHeatNetworkReport,
    pub friction_source: SourceHertzWrench,
    pub source: P,
    pub nonlinear: NewtonReport,
    pub trials: usize,
    pub target_time_s: f64,
    pub body_temperatures_before_k: [f64; 2],
    pub body_temperatures_after_k: [f64; 2],
    /// Actual kinetic plus caloric change, less outside heat, force/motor/plane
    /// work and the separately bounded native static constraint-work residual.
    /// The latter is a numerical account and is never deposited as heat.
    pub combined_energy_residual_j: f64,
}

struct Staged<P> {
    mechanical: SpherePlaneHertzSupportAdvance,
    friction_source: SourceHertzWrench,
    source: P,
}

/// Propose a complete heated spin-brake interval through the original owners.
///
/// The sphere COM and plane must be stationary, applied force normal and motor
/// couple parallel to the actual normal. Native support phases must retain
/// exactly zero centre slip and tangential drive; no velocity or momentum is
/// projected to enter this invariant subspace. The original circular thermal
/// face must match the actual centre, normal and radius. Route A lies on the
/// sphere side, route B on the plane side; both original body IDs are checked.
/// Caloric/conductivity source identity and contact-column applicability remain
/// with their material owner. No whole-sphere mass or thermal geometry is inferred.
///
/// Each nonlinear trial restarts the same opaque accepted mechanical state.
/// A single source snapshot at the solved face temperature supplies all native
/// stopping/holding/reversal events in that window. The native radius and actual
/// normal load remain unchanged across source phases. Only physical sustained
/// friction heats the two existing cells; native static-work defects remain
/// explicit, and any impact, flight or displaced footprint refuses.
///
/// Backward-Euler thermal transport and endpoint coefficient sampling are first
/// order. Exact held-coefficient mechanical stopping is not an exact stop under
/// continuously changing temperature. Refine thermal windows independently;
/// nonlinear convergence does not prove uniqueness or physical calibration.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn propose_sphere_plane_hertz_heat<'m, P, E>(
    scene: &SpherePlane,
    state: SpherePlaneHertzState,
    torque_world_nm: Vec3,
    contact: SpherePlaneHertzHeatContact<'_>,
    thermal: &DepositionSimulation<'m>,
    route: &FrictionHeatRoute,
    dt_s: f64,
    controls: SpherePlaneHertzHeatControls,
    energy_guess_j: &[f64],
    temperature_guess_k: f64,
    resolve: impl Fn(&DepositionSimulation<'m>, f64) -> Result<SpherePlaneHertzHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<HeatedSpherePlaneHertzCandidate<'m, P>, SpherePlaneHertzHeatError<E>> {
    let checkpoint = || {
        if gate.is_requested() {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    };
    checkpoint()?;
    controls.newton.validate()?;
    let origin = state.motion();
    let target = origin.time_s() + dt_s;
    let cells = thermal.cells().len();
    let dimension = cells
        .checked_add(1)
        .ok_or(SpherePlaneHertzHeatError::Invalid(
            "thermal dimension overflow",
        ))?;
    if origin.time_s() != thermal.time_s()
        || !dt_s.is_finite()
        || dt_s <= 0.0
        || !target.is_finite()
        || target <= origin.time_s()
        || scene.velocity_world != Vec3::ZERO
        || origin.body().linear_momentum_world() != Vec3::ZERO
        || !contact.contact_radius_m.is_finite()
        || contact.contact_radius_m <= 0.0
        || contact.contact_radius_m > scene.radius
        || contact.thermal_body_ids[0] == contact.thermal_body_ids[1]
        || cells == 0
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
            controls.geometry_tolerance_m,
            controls.combined_energy_tolerance_j,
            controls.heat.cell_j,
            controls.heat.thermal_j,
            controls.heat.interface_j,
            controls.heat.total_j,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
        || energy_guess_j.len() != cells
        || energy_guess_j.iter().any(|x| !x.is_finite())
    {
        return Err(SpherePlaneHertzHeatError::Invalid(
            "stationary spin geometry, synchronized clocks and finite bounded controls required",
        ));
    }
    let step = FrictionHeatNetworkStep::new(thermal, target, &[route], 1, checkpoint)?;
    let mut before_k = [0.0; 2];
    for (i, id) in route.cells().into_iter().enumerate() {
        if thermal.cells()[&id].body_id != contact.thermal_body_ids[i] {
            return Err(SpherePlaneHertzHeatError::Invalid(
                "original thermal bodies must retain ordered sphere-A/plane-B identity",
            ));
        }
        before_k[i] = thermal.temperature_k(id)?;
    }
    // Use the same representable interval for both native clocks.
    let duration = step.duration_s();
    let trials = Cell::new(0usize);
    let evaluate = |x: &[f64], entropy: bool| {
        checkpoint()?;
        if trials.get() >= controls.maximum_trials {
            return Err(SpherePlaneHertzHeatError::TrialBudget);
        }
        trials.set(trials.get() + 1);
        if x.len() != dimension || x.iter().any(|x| !x.is_finite()) {
            return Err(SpherePlaneHertzHeatError::Invalid(
                "complete finite thermal iterate required",
            ));
        }
        let temperature = temperature_guess_k + x[cells] / controls.temperature_scale_j_per_k;
        let endpoint = |endpoint: &DepositionSimulation<'m>, temperatures: &[f64]| {
            let temperature = temperatures[0];
            let mut snapshot: Option<SourceHertzWrench> = None;
            let mechanical = advance_sphere_plane_hertz_support_with(
                scene,
                state,
                torque_world_nm,
                &mut (),
                duration,
                controls.mechanical,
                |_, phase, _| {
                    admit_phase(
                        scene,
                        origin.time_s(),
                        torque_world_nm,
                        contact.contact_radius_m,
                        route,
                        phase,
                        controls.geometry_tolerance_m,
                    )
                    .map_err(SpherePlaneHertzHeatPhaseError::Invalid)?;
                    if let Some(source) = &snapshot {
                        if source.law().parameters().normal_force_n != phase.normal_force_n {
                            return Err(SpherePlaneHertzHeatPhaseError::Invalid(
                                "held normal support changed",
                            ));
                        }
                    } else {
                        if phase.time_s != origin.time_s() {
                            return Err(SpherePlaneHertzHeatPhaseError::Invalid(
                                "initial sustained contact required",
                            ));
                        }
                        snapshot = Some(
                            contact
                                .source
                                .hertz_at(
                                    temperature,
                                    contact.nominal_pressure_pa,
                                    phase.normal_force_n,
                                    contact.contact_radius_m,
                                )
                                .map_err(SpherePlaneHertzHeatPhaseError::Friction)?,
                        );
                    }
                    Ok((*snapshot.as_ref().expect("source initialized").law(), ()))
                },
                || gate.is_requested(),
            )
            .map_err(SpherePlaneHertzHeatError::Mechanical)?;
            if mechanical.impact_loss_j != 0.0
                || mechanical
                    .motion
                    .events
                    .iter()
                    .any(|e| e.kind == SpherePlaneEventKind::Impact)
                || mechanical.motion.state_after.body().pose().position_world()
                    != origin.body().pose().position_world()
                || mechanical.motion.state_after.body().linear_momentum_world() != Vec3::ZERO
                || mechanical.support_slip_m != [0.0; 2]
            {
                return Err(SpherePlaneHertzHeatError::Invalid(
                    "moving footprints and impacts require separate thermal owners",
                ));
            }
            let friction_source = snapshot.ok_or(SpherePlaneHertzHeatError::Invalid(
                "positive-load sustained contact required",
            ))?;
            let SpherePlaneHertzHeatSource {
                conductivity_w_m_k,
                transport,
            } = resolve(endpoint, temperature).map_err(SpherePlaneHertzHeatError::Source)?;
            let StepCoefficients {
                sources_w,
                conductances_w_k,
                payload,
            } = transport;
            Ok(FrictionHeatNetworkInput {
                interfaces: vec![InterfaceHeatInput {
                    dissipated_j: mechanical.sustained_friction_loss_j,
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
        return Err(SpherePlaneHertzHeatError::Interrupted(result.progress));
    }
    if !result.report.converged {
        return Err(SpherePlaneHertzHeatError::NotConverged(result.report));
    }
    let trial = evaluate(&solver.x, controls.audit_entropy)?;
    let m = &trial.payload().mechanical.motion;
    let delta_k = origin
        .body()
        .kinetic_energy_change_to(scene.properties, m.state_after.body())
        .map_err(|e| {
            SpherePlaneHertzHeatError::Mechanical(SpherePlaneHertzSupportError::Mechanics(
                SpherePlaneError::Dynamics(e),
            ))
        })?;
    let combined = finite(
        delta_k + trial.report().thermal.stored_energy_change_j
            - trial.report().external_heat_j
            - m.external_force_work
            - m.external_torque_work
            - m.plane_work
            - m.static_contact_work_residual,
    )?;
    if combined.abs() > controls.combined_energy_tolerance_j {
        return Err(SpherePlaneHertzHeatError::Invalid(
            "native kinetic and finite caloric work do not balance",
        ));
    }
    let (thermal, staged, heat) = trial.accept(controls.heat, checkpoint)?;
    let mut after_k = [0.0; 2];
    for (i, id) in route.cells().into_iter().enumerate() {
        after_k[i] = thermal.temperature_k(id)?;
    }
    checkpoint()?;
    Ok(HeatedSpherePlaneHertzCandidate {
        mechanical: staged.mechanical,
        thermal,
        heat,
        friction_source: staged.friction_source,
        source: staged.source,
        nonlinear: result.report,
        trials: trials.get(),
        target_time_s: target,
        body_temperatures_before_k: before_k,
        body_temperatures_after_k: after_k,
        combined_energy_residual_j: combined,
    })
}

fn admit_phase(
    scene: &SpherePlane,
    origin_s: f64,
    torque: Vec3,
    radius: f64,
    route: &FrictionHeatRoute,
    phase: &SpherePlaneHertzSupportPhase,
    tolerance: f64,
) -> Result<(), &'static str> {
    if phase.time_s < origin_s
        || phase.normal_force_n <= 0.0
        || phase.body.linear_momentum_world() != Vec3::ZERO
        || phase.slip_velocity_m_s != [0.0; 2]
        || phase.driving_tangent_acceleration_m_s2 != [0.0; 2]
        || scene.force_world.cross(phase.normal_world) != Vec3::ZERO
        || torque.cross(phase.normal_world) != Vec3::ZERO
    {
        return Err("stationary positive normal support and pure normal-axis spin required");
    }
    let geometry = route.geometry;
    let centre = [geometry[0].start_m, geometry[0].end_m]
        .into_iter()
        .find(|p| *p == geometry[1].start_m || *p == geometry[1].end_m)
        .ok_or("original thermal route has no shared face")?;
    let centre = Vec3::new(centre[0], centre[1], centre[2]);
    let thermal_radius = det::sqrt(route.contact_area_m2() / core::f64::consts::PI);
    if norm(phase.contact_point_world_m.sub(centre)) > tolerance
        || !thermal_radius.is_finite()
        || (thermal_radius - radius).abs() > tolerance
    {
        return Err("actual Hertz centre/radius and original thermal face must agree");
    }
    for (side, g) in geometry.into_iter().enumerate() {
        let a = Vec3::new(g.start_m[0], g.start_m[1], g.start_m[2]);
        let b = Vec3::new(g.end_m[0], g.end_m[1], g.end_m[2]);
        let offset = a.sub(centre).scale(0.5).add(b.sub(centre).scale(0.5));
        let signed = offset.dot(phase.normal_world) * if side == 0 { 1.0 } else { -1.0 };
        if norm(b.sub(a).cross(phase.normal_world)) > tolerance
            || !signed.is_finite()
            || signed <= 0.0
        {
            return Err("original thermal sphere-A/plane-B columns must follow the actual normal");
        }
    }
    Ok(())
}
fn norm(value: Vec3) -> f64 {
    det::hypot(det::hypot(value.x, value.y), value.z)
}
fn numerical<E>(
    error: FiniteDifferenceError<SpherePlaneHertzHeatError<E>>,
) -> SpherePlaneHertzHeatError<E> {
    match error {
        FiniteDifferenceError::Source(e) => e,
        FiniteDifferenceError::Invalid(s) => SpherePlaneHertzHeatError::Invalid(s),
    }
}
