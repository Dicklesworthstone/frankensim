//! Rigid bilateral distance loops on the original articulated tree.
//!
//! A central initial constraint impulse and the original implicit midpoint are
//! solved together for endpoint position closure. A second central impulse at
//! that endpoint enforces tangency with the original free-floating mass action.
//! Both impulses retain their actual quadratic work; no position is snapped and
//! no finite-step constraint work is labelled heat or silently discarded.
//!
//! This uses the SHAKE/RATTLE position/hidden-constraint construction described
//! by McLachlan et al., <https://arxiv.org/abs/1207.3367>. The underlying native
//! Lie midpoint is not asserted symplectic, so this composition makes no such
//! claim. Only regular independent distance constraints are admitted.

use super::super::MaterialPoint;
use super::super::bushing::{self, SpatialBushing};
use super::{
    ArticulatedModel, Cell, FallibleFiniteDifference, FiniteDifferenceError,
    FreeFloatingMotionError, FreeFloatingMotionState, ImplicitSpatialConfig, ImplicitSpatialError,
    MAX_FREE_FLOATING_STEP_LINKS, MidpointEquations, MidpointInitial, MomentumMap,
    NewtonKrylovConfig, NewtonKrylovState, NewtonReport, RefCell, SpatialMotionLoads,
    SpatialMotionStep, Vec3, advance, finish_midpoint, poll, validate_config,
    validate_external_wrenches, validate_forward_dynamics_header, validate_slice,
};
use crate::articulated::{
    ArticulatedError, BaseState, FreeFloatingImpulseResponse, JointKind, Kinematics,
    apply_free_floating_impulses, forward_kinematics, free_floating_impulse_delta,
    solve_floating_base_system,
};
use fs_ga::{GaError, Se3, Wrench};
use fs_la::factor::{Cholesky, cholesky};

/// Bound on simultaneous independently reacting distance loops.
pub const MAX_DISTANCE_LOOPS: usize = 32;

/// A positive rigid separation between two distinct original link material points.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DistanceLoop {
    /// First endpoint, unchanged in its link frame.
    pub a: MaterialPoint,
    /// Second endpoint on a different link.
    pub b: MaterialPoint,
    /// Required finite positive distance in metres.
    pub length_m: f64,
}

/// Nonlinear work, physical constraint tolerances and regularity admission.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DistanceLoopConfig {
    /// Existing scaled-rate solve policy. Its evaluation cap includes the
    /// coupled solve's stencils and its actual initial-impulse verification.
    pub nonlinear: ImplicitSpatialConfig,
    /// Maximum actual initial/final distance error, in metres.
    pub position_tolerance_m: f64,
    /// Maximum actual initial/final distance rate, in m/s.
    pub velocity_tolerance_mps: f64,
    /// Minimum mobility/free-body ratio and squared Cholesky pivot of the
    /// unit-diagonal mobility, in [4096*EPSILON,1).
    /// This is a floating regularity admission, not a certified condition bound.
    pub minimum_scaled_pivot: f64,
}

impl Default for DistanceLoopConfig {
    fn default() -> Self {
        Self {
            nonlinear: ImplicitSpatialConfig::default(),
            position_tolerance_m: 1.0e-9,
            velocity_tolerance_mps: 1.0e-9,
            minimum_scaled_pivot: 1.0e-10,
        }
    }
}

/// One position/velocity-admitted constrained step, including both impulse events.
#[derive(Debug, Clone, PartialEq)]
pub struct DistanceLoopStep {
    /// Canonical final constrained state; use this for the next step.
    pub state_after: FreeFloatingMotionState,
    /// Native smooth part between the initial and final constraint impulses.
    /// Its state is the intermediate endpoint BEFORE final tangent correction.
    pub midpoint: SpatialMotionStep,
    /// Actual native initial constraint impulse response, with quadratic work.
    pub initial_reaction: FreeFloatingImpulseResponse,
    /// Actual native final tangent impulse response, with quadratic work.
    pub final_reaction: FreeFloatingImpulseResponse,
    /// Initial impulse per loop, N s; positive pushes endpoints apart.
    pub initial_impulse_ns: Vec<f64>,
    /// Endpoint impulse per loop, N s, in the same convention.
    pub final_impulse_ns: Vec<f64>,
    /// Actual initial signed length residuals, metres.
    pub position_residual_before_m: Vec<f64>,
    /// Actual final signed length residuals, metres; positions are never projected.
    pub position_residual_after_m: Vec<f64>,
    /// Actual initial material-point separation rates, m/s.
    pub velocity_residual_before_mps: Vec<f64>,
    /// Actual final material-point separation rates, m/s.
    pub velocity_residual_after_mps: Vec<f64>,
    /// Minimum admitted unit-diagonal mobility pivot at the initial configuration.
    pub initial_scaled_pivot: f64,
    /// The corresponding endpoint admission; no constraint may silently disappear.
    pub final_scaled_pivot: f64,
    /// Coupled position/midpoint solve report, with the actual impulse rechecked.
    pub nonlinear: NewtonReport,
    /// All nonlinear stage evaluations including directional-difference probes.
    pub residual_evaluations: usize,
    /// Sum of the actual initial/final impulse work, J. Ideal continuous
    /// constraints do zero work, but this finite-step defect is not assumed zero.
    pub constraint_impulse_work: f64,
    /// Actual change of kinetic+finite-spring storage plus damping loss minus
    /// actuator/external work, J. Includes constraint-work and integration defects.
    pub mechanical_work_residual: f64,
    /// Actual net linear momentum change minus the external smooth impulse, N s.
    pub linear_momentum_residual_world: Vec3,
    /// Actual orbital-plus-spin change minus the external smooth impulse, N m s.
    pub angular_momentum_residual_world: Vec3,
}

/// Atomic refusal: no state, arbitrary dependent reaction or clipped pose is returned.
#[derive(Debug, Clone, PartialEq)]
pub enum DistanceLoopError {
    /// Malformed loop, controls, count or nonfinite arithmetic.
    Invalid(&'static str),
    /// Original force/dynamics/nonlinear/budget/cancellation refusal.
    Motion(ImplicitSpatialError),
    /// Input or endpoint is outside the admitted constraint position tolerance.
    Position {
        /// `initial` or `endpoint`.
        stage: &'static str,
        /// Original caller loop index.
        constraint: usize,
        /// Actual length minus declared length, m.
        residual_m: f64,
    },
    /// Initial or final velocity is not tangent within the requested tolerance.
    Velocity {
        /// `initial` or `endpoint`.
        stage: &'static str,
        /// Original caller loop index.
        constraint: usize,
        /// Actual material-point range rate, m/s.
        residual_mps: f64,
    },
    /// Dependent, singular, insufficiently regular or nonreciprocal mobility.
    /// This lane refuses instead of assigning an arbitrary redundant reaction.
    Mobility {
        /// `initial` or `endpoint`.
        stage: &'static str,
        /// First refused loop/pivot in the supplied order.
        constraint: usize,
        /// Concrete floating admission that failed.
        reason: &'static str,
    },
}

impl core::fmt::Display for DistanceLoopError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "articulated distance loop: {self:?}")
    }
}
impl std::error::Error for DistanceLoopError {}
impl From<ImplicitSpatialError> for DistanceLoopError {
    fn from(value: ImplicitSpatialError) -> Self {
        Self::Motion(value)
    }
}
impl From<FreeFloatingMotionError> for DistanceLoopError {
    fn from(value: FreeFloatingMotionError) -> Self {
        Self::Motion(value.into())
    }
}
impl From<ArticulatedError> for DistanceLoopError {
    fn from(value: ArticulatedError) -> Self {
        Self::Motion(value.into())
    }
}
impl From<GaError> for DistanceLoopError {
    fn from(value: GaError) -> Self {
        Self::Motion(value.into())
    }
}

/// Advance a closed tree with regular internal rigid distance loops.
///
/// Initial positions AND velocities must satisfy the declared constraints. The
/// first central impulses and native midpoint rates are simultaneous unknowns;
/// the endpoint tangent solve uses the full free-floating mobility, including
/// base recoil, all joints, off-centre arms and complete inertia tensors. Both
/// impulse applications use the original owner. No mass is clamped or replaced
/// by a penalty spring. Forces, finite springs and dampers retain their separate
/// storage/work account from the unconstrained lane.
///
/// At most 32 independent distance loops and 4096 links are admitted. Cached
/// initial inverse-mass columns cost O(loops*links), with O(loops*joints) storage;
/// each coupled residual adds one original FK pass to the native midpoint work.
/// The endpoint solves only a loops-by-loops SPD mobility. No dense generalized
/// mass/Jacobian matrix or nested nonlinear integration is introduced.
///
/// Returned endpoint residuals are floating observations, not a state enclosure.
/// Intermediate drift states need not lie on the manifold. There is no collision
/// detection, general joint/gear loop, nonholonomic rolling, automatic assembly,
/// redundancy resolution, branch certificate or adjoint. Numerical constraint
/// work remains visible; no general symplecticity or energy-exactness is claimed.
pub fn step_free_floating_distance_loops(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: SpatialMotionLoads<'_>,
    loops: &[DistanceLoop],
    duration: f64,
    config: DistanceLoopConfig,
    cancelled: impl FnMut() -> bool,
) -> Result<DistanceLoopStep, DistanceLoopError> {
    step_free_floating_distance_loops_with_bushings(
        model,
        state,
        loads,
        loops,
        &[],
        duration,
        config,
        cancelled,
    )
}

/// Advance the same rigid loops together with finite spatial bushings.
///
/// Both finite bushing strains/rates and rigid loop reactions use the original
/// tree at the simultaneous stage. All bushing storage, loss and internal work
/// is included in `midpoint` and the complete mechanical-work residual. A rigid
/// distance does not turn the bushing into a penalty constraint or erase its
/// tangential/rotational mechanics. An empty bushing slice is exactly the plain
/// distance-loop operation. All of that operation's limits and refusals apply.
pub fn step_free_floating_distance_loops_with_bushings(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    loads: SpatialMotionLoads<'_>,
    loops: &[DistanceLoop],
    bushings: &[SpatialBushing],
    duration: f64,
    config: DistanceLoopConfig,
    mut cancelled: impl FnMut() -> bool,
) -> Result<DistanceLoopStep, DistanceLoopError> {
    poll(&mut cancelled)?;
    validate_config(config.nonlinear)?;
    if !duration.is_finite() || duration <= 0.0 {
        return Err(FreeFloatingMotionError::InvalidDuration.into());
    }
    if loops.is_empty() || loops.len() > MAX_DISTANCE_LOOPS {
        return Err(DistanceLoopError::Invalid(
            "distance loop count must be in 1..=32",
        ));
    }
    if [
        config.position_tolerance_m,
        config.velocity_tolerance_mps,
        config.minimum_scaled_pivot,
    ]
    .iter()
    .any(|x| !x.is_finite() || *x <= 0.0)
        || config.minimum_scaled_pivot >= 1.0
        || config.minimum_scaled_pivot < 4096. * f64::EPSILON
    {
        return Err(DistanceLoopError::Invalid(
            "finite positive constraint tolerances and resolved mobility pivot required",
        ));
    }
    if model.link_count() > MAX_FREE_FLOATING_STEP_LINKS {
        return Err(FreeFloatingMotionError::TooManyLinks {
            got: model.link_count(),
            maximum: MAX_FREE_FLOATING_STEP_LINKS,
        }
        .into());
    }
    model.free_floating_complexity()?;
    validate_forward_dynamics_header(
        model,
        state.base.twist_body,
        &state.joint_position,
        &state.joint_velocity,
        loads.held.generalized_force,
        loads.held.gravity_world,
        loads.held.external_body_wrench,
    )?;
    validate_external_wrenches(loads.held.external_body_wrench)?;
    bushing::admit(model, loads, bushings, &mut cancelled)?;
    admit_loops(model, loops)?;
    let initial_map = MomentumMap::new(model, state, &mut cancelled)?;
    let initial_coordinates = initial_map.coordinates(state.base)?;
    let initial_load =
        bushing::evaluate(&initial_map, state.base, loads, bushings, &mut cancelled)?.load;
    let before = geometry(&initial_map.shape, loops)?;
    check_position(&before.gaps, config, "initial")?;
    check_velocity(&before.rates, config, "initial")?;
    let initial_mobility = Mobility::new(
        model,
        &state.joint_position,
        &before,
        config,
        "initial",
        true,
        &mut cancelled,
    )?;
    let n = model.dof_count();
    let count = loops.len();
    let mut scales = vec![config.nonlinear.angular_velocity_scale; n + 3];
    for (link, part) in model.links().iter().enumerate() {
        if part.joint.kind() == JointKind::Prismatic
            && let Some(dof) = model.dof_index(link)
        {
            scales[dof] = config.nonlinear.linear_velocity_scale;
        }
    }
    let impulse_scales: Vec<_> = initial_mobility
        .diagonal
        .iter()
        .map(|d| config.nonlinear.linear_velocity_scale / (d * d))
        .collect();
    let gap_scale = duration * config.nonlinear.linear_velocity_scale;
    let target = config
        .nonlinear
        .residual_tolerance
        .min(config.position_tolerance_m / gap_scale);
    if !gap_scale.is_finite()
        || gap_scale <= 0.0
        || !target.is_finite()
        || target <= 0.0
        || impulse_scales.iter().any(|s| !s.is_finite() || *s <= 0.0)
    {
        return Err(DistanceLoopError::Invalid(
            "unrepresentable constraint scaling",
        ));
    }
    let gravity_only_external = loads.world_points.is_empty()
        && loads
            .held
            .external_body_wrench
            .iter()
            .all(|w| w.force == Vec3::default() && w.torque == Vec3::default());
    let mut guess = state.joint_velocity.clone();
    let w = state.base.twist_body.angular;
    guess.extend([w.x, w.y, w.z]);
    for (x, s) in guess.iter_mut().zip(&scales) {
        *x /= s;
    }
    guess.resize(n + 3 + count, 0.0);
    let evaluations = Cell::new(0usize);
    let cancelled = RefCell::new(cancelled);
    let charge = || -> Result<(), DistanceLoopError> {
        let used = evaluations.get();
        if used >= config.nonlinear.maximum_evaluations {
            return Err(ImplicitSpatialError::EvaluationBudget {
                used,
                maximum: config.nonlinear.maximum_evaluations,
            }
            .into());
        }
        evaluations.set(used + 1);
        Ok(())
    };
    let residual = |x: &[f64]| -> Result<Vec<f64>, DistanceLoopError> {
        let mut cancel = cancelled.borrow_mut();
        poll(&mut *cancel)?;
        charge()?;
        let mut kicked = state.clone();
        for (j, column) in initial_mobility.joint_columns.iter().enumerate() {
            let magnitude = x[n + 3 + j] * impulse_scales[j];
            for (v, dv) in kicked.joint_velocity.iter_mut().zip(column) {
                *v += magnitude * dv;
            }
        }
        let midpoint = MidpointEquations {
            model,
            state: &kicked,
            loads,
            duration,
            initial: initial_coordinates,
            mass: initial_map.mass(),
            scales: &scales,
            gravity_only_external,
            bushings,
        }
        .evaluate(&x[..n + 3], &mut *cancel)?;
        let q = advance(
            &state.joint_position,
            &midpoint.state.joint_velocity,
            duration,
        );
        poll(&mut *cancel)?;
        let shape = forward_kinematics(
            model,
            BaseState::stationary(Se3::identity()),
            &q,
            &vec![0.; n],
        )?;
        let end = geometry(&shape, loops)?;
        let mut residual = midpoint.residual;
        residual.extend(end.gaps.iter().map(|g| g / gap_scale));
        validate_slice(&residual, "distance_loop.equations")?;
        poll(&mut *cancel)?;
        Ok(residual)
    };
    let problem =
        FallibleFiniteDifference::new(n + 3 + count, config.nonlinear.difference_step, residual)
            .map_err(difference_error)?;
    let solver_config = NewtonKrylovConfig {
        absolute_tolerance: target,
        relative_tolerance: 0.,
        linear_restart: config.nonlinear.krylov_restart.min(n + 3 + count),
        max_linear_cycles: 8,
        forcing_minimum: 1e-3,
        forcing_maximum: 0.1,
        ..Default::default()
    };
    let created = NewtonKrylovState::new(&problem, guess, solver_config);
    if let Some(error) = problem.take_failure() {
        return Err(difference_error(error));
    }
    let mut solver = created.map_err(ImplicitSpatialError::Solver)?;
    let run = solver.run_controlled(&problem, config.nonlinear.maximum_iterations, &mut || {
        problem.has_failed() || (*cancelled.borrow_mut())()
    });
    if let Some(error) = problem.take_failure() {
        return Err(difference_error(error));
    }
    match run.progress {
        fs_solver::nonlinear::SolverRunProgress::Complete => {}
        fs_solver::nonlinear::SolverRunProgress::Paused => {
            return Err(FreeFloatingMotionError::Cancelled.into());
        }
        other => return Err(ImplicitSpatialError::SolverProgress(other).into()),
    }
    let mut report = run.report;
    if !report.converged {
        return Err(ImplicitSpatialError::Unresolved(report).into());
    }
    let initial_impulse: Vec<_> = solver.x[n + 3..]
        .iter()
        .zip(&impulse_scales)
        .map(|(x, s)| x * s)
        .collect();
    let mut cancel = cancelled.borrow_mut();
    poll(&mut *cancel)?;
    let initial_reaction = apply(model, state, &before, &initial_impulse)?;
    let kicked = FreeFloatingMotionState {
        base: initial_reaction.base_after,
        joint_position: state.joint_position.clone(),
        joint_velocity: initial_reaction.joint_velocity_after.clone(),
    };
    // Re-evaluate with the ACTUAL native summed impulse, not cached unit-column
    // superposition or assumed zero floating momentum defect.
    charge()?;
    let kicked_map = MomentumMap::new(model, &kicked, &mut *cancel)?;
    let kicked_coordinates = kicked_map.coordinates(kicked.base)?;
    let actual_midpoint = MidpointEquations {
        model,
        state: &kicked,
        loads,
        duration,
        initial: kicked_coordinates,
        mass: kicked_map.mass(),
        scales: &scales,
        gravity_only_external,
        bushings,
    }
    .evaluate(&solver.x[..n + 3], &mut *cancel)?;
    let final_q = advance(
        &kicked.joint_position,
        &actual_midpoint.state.joint_velocity,
        duration,
    );
    let shape = forward_kinematics(
        model,
        BaseState::stationary(Se3::identity()),
        &final_q,
        &vec![0.; n],
    )?;
    let final_geometry = geometry(&shape, loops)?;
    let mut actual_residual = actual_midpoint.residual.clone();
    actual_residual.extend(final_geometry.gaps.iter().map(|g| g / gap_scale));
    report.residual_norm = fs_solver::norm2(&actual_residual);
    if !report.residual_norm.is_finite() || report.residual_norm > target {
        report.converged = false;
        return Err(ImplicitSpatialError::Unresolved(report).into());
    }
    let midpoint = finish_midpoint(
        model,
        &kicked,
        loads,
        duration,
        MidpointInitial {
            map: kicked_map,
            coordinates: kicked_coordinates,
            load: initial_load,
            momentum: initial_reaction.momentum_after,
        },
        actual_midpoint,
        bushings,
        &mut *cancel,
    )?;
    let endpoint = &midpoint.motion.state_after;
    let final_map = MomentumMap::new(model, endpoint, &mut *cancel)?;
    let end = geometry(&final_map.shape, loops)?;
    check_position(&end.gaps, config, "endpoint")?;
    let final_mobility = Mobility::new(
        model,
        &endpoint.joint_position,
        &end,
        config,
        "endpoint",
        false,
        &mut *cancel,
    )?;
    let final_impulse = final_mobility.cancel_rates(&end.rates)?;
    poll(&mut *cancel)?;
    let final_reaction = apply(model, endpoint, &end, &final_impulse)?;
    let state_after = FreeFloatingMotionState {
        base: final_reaction.base_after,
        joint_position: endpoint.joint_position.clone(),
        joint_velocity: final_reaction.joint_velocity_after.clone(),
    };
    let final_shape = forward_kinematics(
        model,
        BaseState::stationary(Se3::identity()),
        &state_after.joint_position,
        &state_after.joint_velocity,
    )?;
    let actual = geometry(&final_shape, loops)?;
    check_position(&actual.gaps, config, "endpoint")?;
    check_velocity(&actual.rates, config, "endpoint")?;
    let initial_momentum = initial_reaction.momentum_before;
    let final_momentum = final_reaction.momentum_after;
    let constraint_impulse_work = initial_reaction.impulse_work + final_reaction.impulse_work;
    let mechanical_work_residual = final_momentum.kinetic_energy - initial_momentum.kinetic_energy
        + midpoint.stored_energy_after
        - midpoint.stored_energy_before
        + midpoint.dissipated_energy
        - midpoint.motion.external_work
        - midpoint.motion.actuator_work;
    let linear_momentum_residual_world = final_momentum.linear_momentum_world
        - initial_momentum.linear_momentum_world
        - midpoint.motion.external_linear_impulse_world;
    let angular_momentum_residual_world = final_momentum.angular_momentum_world
        - initial_momentum.angular_momentum_world
        - midpoint.motion.external_angular_impulse_world;
    validate_slice(
        &[constraint_impulse_work, mechanical_work_residual],
        "distance_loop.work",
    )?;
    super::validate_vec3(
        linear_momentum_residual_world,
        "distance_loop.linear_momentum",
    )?;
    super::validate_vec3(
        angular_momentum_residual_world,
        "distance_loop.angular_momentum",
    )?;
    poll(&mut *cancel)?;
    Ok(DistanceLoopStep {
        state_after,
        midpoint,
        initial_reaction,
        final_reaction,
        initial_impulse_ns: initial_impulse,
        final_impulse_ns: final_impulse,
        position_residual_before_m: before.gaps,
        position_residual_after_m: actual.gaps,
        velocity_residual_before_mps: before.rates,
        velocity_residual_after_mps: actual.rates,
        initial_scaled_pivot: initial_mobility.minimum_pivot,
        final_scaled_pivot: final_mobility.minimum_pivot,
        nonlinear: report,
        residual_evaluations: evaluations.get(),
        constraint_impulse_work,
        mechanical_work_residual,
        linear_momentum_residual_world,
        angular_momentum_residual_world,
    })
}

struct Geometry {
    gaps: Vec<f64>,
    rates: Vec<f64>,
    terms: Vec<[(usize, Wrench); 2]>,
}

fn admit_loops(model: &ArticulatedModel, loops: &[DistanceLoop]) -> Result<(), DistanceLoopError> {
    for c in loops {
        if c.a.link >= model.link_count()
            || c.b.link >= model.link_count()
            || c.a.link == c.b.link
            || !c.length_m.is_finite()
            || c.length_m <= 0.
        {
            return Err(DistanceLoopError::Invalid("distance endpoints/length"));
        }
        for p in [c.a, c.b] {
            super::validate_vec3(p.point_body, "distance_loop.material_point")?;
        }
    }
    Ok(())
}

fn geometry(shape: &Kinematics, loops: &[DistanceLoop]) -> Result<Geometry, DistanceLoopError> {
    let mut gaps = Vec::with_capacity(loops.len());
    let mut rates = Vec::with_capacity(loops.len());
    let mut terms = Vec::with_capacity(loops.len());
    for c in loops {
        let a = shape.world_from_link[c.a.link];
        let b = shape.world_from_link[c.b.link];
        let delta = b.transform_point(c.b.point_body)? - a.transform_point(c.a.point_body)?;
        let length = delta.norm();
        if !length.is_finite() || length <= 0. {
            return Err(DistanceLoopError::Invalid(
                "coincident/nonfinite distance endpoints",
            ));
        }
        let normal = delta.scale(1. / length);
        let va = shape.body_twist[c.a.link];
        let vb = shape.body_twist[c.b.link];
        // shape is relative to a stationary base: bulk translation and the
        // identically zero n dot (omega cross delta) never pollute the rate.
        let relative = b
            .rotation()
            .rotate(vb.linear + vb.angular.cross(c.b.point_body))?
            - a.rotation()
                .rotate(va.linear + va.angular.cross(c.a.point_body))?;
        let fa = a.rotation().inverse().rotate(normal.scale(-1.))?;
        let fb = b.rotation().inverse().rotate(normal)?;
        gaps.push(length - c.length_m);
        rates.push(normal.dot(relative));
        terms.push([
            (c.a.link, Wrench::new(c.a.point_body.cross(fa), fa)),
            (c.b.link, Wrench::new(c.b.point_body.cross(fb), fb)),
        ]);
    }
    validate_slice(&gaps, "distance_loop.gaps")?;
    validate_slice(&rates, "distance_loop.rates")?;
    Ok(Geometry { gaps, rates, terms })
}

struct Mobility {
    diagonal: Vec<f64>,
    factor: Cholesky,
    minimum_pivot: f64,
    joint_columns: Vec<Vec<f64>>,
}
impl Mobility {
    fn new(
        model: &ArticulatedModel,
        q: &[f64],
        geometry: &Geometry,
        config: DistanceLoopConfig,
        stage: &'static str,
        retain_columns: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, DistanceLoopError> {
        let k = geometry.terms.len();
        let mut matrix = vec![0.; k * k];
        let mut joint_columns = Vec::with_capacity(if retain_columns { k } else { 0 });
        let joints = vec![0.; model.dof_count()];
        let mut body = vec![Wrench::default(); model.link_count()];
        for (column, terms) in geometry.terms.iter().enumerate() {
            poll(cancelled)?;
            for (link, wrench) in terms {
                body[*link] = *wrench;
            }
            let response = free_floating_impulse_delta(model, q, &joints, &body)?;
            for (row, row_terms) in geometry.terms.iter().enumerate() {
                matrix[row * k + column] = row_terms
                    .iter()
                    .map(|(link, w)| w.pairing(response.body[*link]))
                    .sum();
            }
            if retain_columns {
                joint_columns.push(response.joint);
            }
            for (link, _) in terms {
                body[*link] = Wrench::default();
            }
        }
        validate_slice(&matrix, "distance_loop.mobility")?;
        let mut diagonal = Vec::with_capacity(k);
        for i in 0..k {
            let d = matrix[i * k + i];
            let mut free_mobility = 0.;
            for (link, wrench) in geometry.terms[i] {
                let response =
                    solve_floating_base_system(&model.links()[link].inertia.matrix()?, wrench)?;
                free_mobility += wrench.pairing(response);
            }
            if d <= 0.
                || !d.is_finite()
                || !free_mobility.is_finite()
                || free_mobility <= 0.
                || d / free_mobility < config.minimum_scaled_pivot
            {
                return Err(DistanceLoopError::Mobility {
                    stage,
                    constraint: i,
                    reason: "constraint mobility vanishes relative to free endpoint bodies",
                });
            }
            diagonal.push(d.sqrt());
        }
        for i in 0..k {
            matrix[i * k + i] = 1.;
            for j in 0..i {
                let a = (matrix[i * k + j] / diagonal[i]) / diagonal[j];
                let b = (matrix[j * k + i] / diagonal[j]) / diagonal[i];
                if !a.is_finite() || !b.is_finite() || (a - b).abs() > 1e-10 {
                    return Err(DistanceLoopError::Mobility {
                        stage,
                        constraint: i,
                        reason: "constraint mobility reciprocity unresolved",
                    });
                }
                matrix[i * k + j] = 0.5 * (a + b);
                matrix[j * k + i] = matrix[i * k + j];
            }
        }
        poll(cancelled)?;
        let factor = cholesky(&matrix, k).map_err(|error| DistanceLoopError::Mobility {
            stage,
            constraint: match error {
                fs_la::factor::FactorError::NotSpd { index }
                | fs_la::factor::FactorError::Singular { index } => index,
            },
            reason: "dependent or indefinite constraint mobility",
        })?;
        let mut minimum_pivot = 1.0_f64;
        for i in 0..k {
            let pivot = factor.l(i, i) * factor.l(i, i);
            if !pivot.is_finite() || pivot < config.minimum_scaled_pivot {
                return Err(DistanceLoopError::Mobility {
                    stage,
                    constraint: i,
                    reason: "constraint mobility is insufficiently regular",
                });
            }
            minimum_pivot = minimum_pivot.min(pivot);
        }
        poll(cancelled)?;
        Ok(Self {
            diagonal,
            factor,
            minimum_pivot,
            joint_columns,
        })
    }
    fn cancel_rates(&self, rates: &[f64]) -> Result<Vec<f64>, DistanceLoopError> {
        let mut impulse: Vec<_> = rates
            .iter()
            .zip(&self.diagonal)
            .map(|(v, d)| -v / d)
            .collect();
        self.factor.solve(&mut impulse);
        for (j, d) in impulse.iter_mut().zip(&self.diagonal) {
            *j /= d;
        }
        validate_slice(&impulse, "distance_loop.endpoint_impulse")?;
        Ok(impulse)
    }
}

fn apply(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    geometry: &Geometry,
    magnitude: &[f64],
) -> Result<FreeFloatingImpulseResponse, DistanceLoopError> {
    let mut body = vec![Wrench::default(); model.link_count()];
    for (terms, scale) in geometry.terms.iter().zip(magnitude) {
        for (link, w) in terms {
            body[*link] = Wrench::new(
                body[*link].torque + w.torque.scale(*scale),
                body[*link].force + w.force.scale(*scale),
            );
        }
    }
    Ok(apply_free_floating_impulses(
        model,
        state.base,
        &state.joint_position,
        &state.joint_velocity,
        &vec![0.; model.dof_count()],
        &body,
    )?)
}
fn check_position(
    values: &[f64],
    config: DistanceLoopConfig,
    stage: &'static str,
) -> Result<(), DistanceLoopError> {
    for (constraint, &residual_m) in values.iter().enumerate() {
        if residual_m.abs() > config.position_tolerance_m {
            return Err(DistanceLoopError::Position {
                stage,
                constraint,
                residual_m,
            });
        }
    }
    Ok(())
}
fn check_velocity(
    values: &[f64],
    config: DistanceLoopConfig,
    stage: &'static str,
) -> Result<(), DistanceLoopError> {
    for (constraint, &residual_mps) in values.iter().enumerate() {
        if residual_mps.abs() > config.velocity_tolerance_mps {
            return Err(DistanceLoopError::Velocity {
                stage,
                constraint,
                residual_mps,
            });
        }
    }
    Ok(())
}
fn difference_error(error: FiniteDifferenceError<DistanceLoopError>) -> DistanceLoopError {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(reason) => ImplicitSpatialError::Difference(reason).into(),
    }
}
