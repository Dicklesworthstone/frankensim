//! Original nonlinear conductor response coupled to native rigid momentum.
//!
//! This is the same fixed-configuration force substep as the linear owner:
//! pose, geometry, imposed B and external source context are held. The original
//! source law closes at E'=E+u cross B and J at each proposed rigid midpoint.
//! A nonlinear field is never replaced by an affine motor/drag fit.

use super::*;
use fs_em::spatial::conduction::response::{
    SpatialConductionResponse, SpatialConductionResponseControls, SpatialConductionResponseError,
    motion::{
        SpatialMovingConductionResponseCandidate, SpatialMovingConductionResponseDrive,
        SpatialMovingConductionResponseReport,
    },
};
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_material::{conductor::response::ConductorMagnitude, state_point::MaterialStatePointError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
    finite_difference::{FallibleFiniteDifference, FiniteDifferenceError},
};
use std::cell::{Cell, RefCell};

/// Native source/field limits and independent physical momentum/work gates.
/// Newton coordinates are midpoint-velocity increments divided by the explicit
/// angular/linear scales; difference_step is dimensionless in those coordinates.
#[derive(Clone, Copy, Debug)]
pub struct SpatialConductionResponseKickControls {
    /// maximum_field_solves counts ALL native solves across numerical stencils,
    /// nonlinear constitutive iterations and final fresh admission in this call.
    pub physical: SpatialConductionKickControls,
    /// Per-field constitutive limits, also narrowed by the global remaining cap.
    pub response: SpatialConductionResponseControls,
    pub maximum_newton_iterations: usize,
    pub angular_velocity_scale_rad_s: f64,
    pub linear_velocity_scale_m_s: f64,
    pub difference_step: f64,
    pub newton: NewtonKrylovConfig,
}

/// Original source, material, native momentum and numerical refusals.
#[derive(Debug)]
pub enum SpatialConductionResponseKickError {
    Kick(SpatialConductionKickError),
    Response(SpatialConductionResponseError),
    Material(MaterialStatePointError),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    FieldSolveBudget {
        completed: usize,
        required: usize,
        limit: usize,
    },
    Cancelled,
}
impl core::fmt::Display for SpatialConductionResponseKickError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "nonlinear spatial conduction momentum kick: {self:?}")
    }
}
impl std::error::Error for SpatialConductionResponseKickError {}
impl From<SpatialConductionKickError> for SpatialConductionResponseKickError {
    fn from(error: SpatialConductionKickError) -> Self {
        match error {
            SpatialConductionKickError::Cancelled => Self::Cancelled,
            error => Self::Kick(error),
        }
    }
}
impl From<SpatialConductionResponseError> for SpatialConductionResponseKickError {
    fn from(error: SpatialConductionResponseError) -> Self {
        match error {
            SpatialConductionResponseError::Cancelled => Self::Cancelled,
            error => Self::Response(error),
        }
    }
}
impl From<MaterialStatePointError> for SpatialConductionResponseKickError {
    fn from(error: MaterialStatePointError) -> Self {
        Self::Material(error)
    }
}
type Result<T> = core::result::Result<T, SpatialConductionResponseKickError>;

/// Accepted nonlinear current/source state and actual native impulse response.
/// Joule heat is reported for its original thermal receivers, not deposited here.
#[derive(Debug)]
pub struct SpatialConductionResponseKick {
    pub mechanical: FreeFloatingImpulseResponse,
    pub field: SpatialMovingConductionResponseReport,
    pub midpoint_twist_body: Twist,
    pub lorentz_wrench_body: Wrench,
    pub terminal_work_j: f64,
    pub joule_heat_j: f64,
    pub mechanical_work_j: f64,
    pub external_work_j: f64,
    pub angular_velocity_residual_rad_s: f64,
    pub linear_velocity_residual_m_s: f64,
    pub mechanical_work_residual_j: f64,
    /// delta_kinetic + Joule - terminal - external [J].
    pub total_energy_residual_j: f64,
    /// All actual field solves, including discarded numerical/constitutive trials.
    pub field_solves: usize,
    pub source_evaluations: usize,
    pub nonlinear: NewtonReport,
}

/// Consume the accepted original source owner and momentum response together.
/// Failed work, source-domain queries and cancellation mutate neither input.
pub struct SpatialConductionResponseKickCandidate {
    problem: SpatialConductionResponse,
    report: SpatialConductionResponseKick,
}
impl SpatialConductionResponseKickCandidate {
    pub fn problem(&self) -> &SpatialConductionResponse {
        &self.problem
    }
    pub fn report(&self) -> &SpatialConductionResponseKick {
        &self.report
    }
    pub fn into_parts(self) -> (SpatialConductionResponse, SpatialConductionResponseKick) {
        (self.problem, self.report)
    }
}

/// Close actual original-source current and native rigid momentum at one
/// implicit midpoint, with the caller's held boundary/direct-electrode drive.
/// The same physical velocity produces motional EMF and receives Lorentz work.
/// Source/geometry snapshots are fixed during the force kick; no induction,
/// pose transport, gyroscopic evolution or thermal timeline is inferred.
#[allow(clippy::too_many_arguments)]
pub fn kick_response_candidate<'d>(
    problem: &SpatialConductionResponse,
    mass: &TetrahedralMass<'_>,
    drive: impl Into<SpatialMovingConductionResponseDrive<'d>>,
    flux_density_t: &[[f64; 3]],
    duration_s: f64,
    initial: FreeFloatingBaseState,
    external_wrench_body: Wrench,
    assembly: WeightedAssemblyLimits,
    controls: SpatialConductionResponseKickControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<SpatialConductionResponseKickCandidate> {
    kick_response_with_materials_candidate(
        problem,
        problem.conductors().to_vec(),
        mass,
        drive,
        flux_density_t,
        duration_s,
        initial,
        external_wrench_body,
        assembly,
        controls,
        cancelled,
    )
}

/// The same momentum/current closure at an externally owned thermal/phase
/// endpoint. The original density bundle is re-evaluated with its exact pins at
/// EVERY supplied source context and must retain the original density values.
/// Its response and applicability must be independent of the owned E'/J axes.
/// Thus changed electrical coordinates do not freeze an unowned inertia law.
///
/// A thermal caller additionally owns constant-density caloric admission and
/// support over its thermal step. This force substep changes neither inventory
/// nor inertia. It admits complete original contact-free tetrahedral bodies;
/// moving electrical interfaces need their own surface-work law.
///
/// Every residual/Jacobian evaluation starts from the SAME initial mechanical
/// history and supplied external material context. The final accepted current
/// is freshly solved, its original load is passed to native impulses once, and
/// independent midpoint, Lorentz/kinetic and total energy gates must all pass.
#[allow(clippy::too_many_arguments)]
pub fn kick_response_with_materials_candidate<'d>(
    problem: &SpatialConductionResponse,
    seeds: Vec<ResolvedConductorState>,
    mass: &TetrahedralMass<'_>,
    drive: impl Into<SpatialMovingConductionResponseDrive<'d>>,
    flux_density_t: &[[f64; 3]],
    duration_s: f64,
    initial: FreeFloatingBaseState,
    external_wrench_body: Wrench,
    assembly: WeightedAssemblyLimits,
    controls: SpatialConductionResponseKickControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<SpatialConductionResponseKickCandidate> {
    poll(cancelled)?;
    let physical = controls.physical;
    if ![
        duration_s,
        physical.angular_velocity_tolerance_rad_s,
        physical.linear_velocity_tolerance_m_s,
        physical.energy_tolerance_j,
        controls.angular_velocity_scale_rad_s,
        controls.linear_velocity_scale_m_s,
        controls.difference_step,
    ]
    .iter()
    .all(|v| v.is_finite() && *v > 0.0)
        || physical.maximum_field_solves == 0
        || controls.response.maximum_field_solves == 0
        || controls.maximum_newton_iterations == 0
        || !(1..=fs_em::conduction::response::MAX_RESPONSE_ITERATIONS)
            .contains(&controls.response.maximum_iterations)
    {
        return Err(invalid(
            "positive bounded nonlinear momentum/source controls required",
        ));
    }
    let cells = problem.mesh().tets.len();
    let vertices = problem.positions_m().len();
    if cells == 0
        || cells > physical.loads.max_cells
        || vertices > physical.loads.max_vertices
        || flux_density_t.len() != cells
        || !problem.contacts().is_empty()
        || seeds.len() != cells
        || mass.nodes_m().len() != vertices
        || mass.tetrahedra().len() != cells
    {
        return Err(invalid(
            "complete bounded contact-free nonlinear field and original mass required",
        ));
    }
    admit_mass(problem, &seeds, mass, cancelled)?;
    let model = ArticulatedModel::new(vec![Link::new(
        "original nonlinear tetrahedral conductor",
        None,
        Se3::identity(),
        JointModel::FIXED,
        mass.reference().inertia(),
    )])
    .map_err(SpatialConductionKickError::from)?;
    let drive = drive.into();
    let selected: Vec<_> = (0..cells).collect();
    let old = initial.twist_body.to_array();
    let external = external_wrench_body.to_array();
    for value in old.into_iter().chain(external) {
        finite(value)?;
    }
    let scales = [
        controls.angular_velocity_scale_rad_s,
        controls.angular_velocity_scale_rad_s,
        controls.angular_velocity_scale_rad_s,
        controls.linear_velocity_scale_m_s,
        controls.linear_velocity_scale_m_s,
        controls.linear_velocity_scale_m_s,
    ];
    let cancel = RefCell::new(cancelled);
    let field_solves = Cell::new(0usize);
    let source_evaluations = Cell::new(0usize);
    let evaluate = |x: &[f64]| -> Result<Trial> {
        poll(&mut || (*cancel.borrow_mut())())?;
        if x.len() != 6 {
            return Err(invalid("complete scaled midpoint coordinates required"));
        }
        let mut midpoint = [0.0; 6];
        for i in 0..6 {
            midpoint[i] = finite(old[i] + finite(scales[i] * x[i])?)?;
        }
        let midpoint_twist_body = twist(midpoint);
        let (nodal, cell_velocity) = rigid_velocity_trace(
            problem.positions_m(),
            &problem.mesh().tets,
            midpoint_twist_body,
            &mut || (*cancel.borrow_mut())(),
        )?;
        let remaining = physical.maximum_field_solves - field_solves.get();
        if remaining == 0 {
            return Err(SpatialConductionResponseKickError::FieldSolveBudget {
                completed: field_solves.get(),
                required: 1,
                limit: physical.maximum_field_solves,
            });
        }
        let bounded = SpatialConductionResponseControls {
            maximum_field_solves: controls.response.maximum_field_solves.min(remaining),
            ..controls.response
        };
        let candidate = problem
            .solve_with_motion_and_materials_and_contacts_candidate(
                seeds.clone(),
                Vec::new(),
                drive,
                SpatialConductionMotion {
                    velocity_m_s: &cell_velocity,
                    flux_density_t,
                },
                assembly,
                physical.field,
                bounded,
                &mut || (*cancel.borrow_mut())(),
            )
            .map_err(|error| match error {
                SpatialConductionResponseError::FieldSolveBudget {
                    completed,
                    required,
                    ..
                } if bounded.maximum_field_solves < controls.response.maximum_field_solves => {
                    SpatialConductionResponseKickError::FieldSolveBudget {
                        completed: field_solves.get() + completed,
                        required,
                        limit: physical.maximum_field_solves,
                    }
                }
                other => other.into(),
            })?;
        let solves = candidate.report().field_solves;
        if solves == 0 || solves > remaining {
            return Err(SpatialConductionResponseKickError::FieldSolveBudget {
                completed: field_solves.get(),
                required: solves,
                limit: physical.maximum_field_solves,
            });
        }
        field_solves.set(field_solves.get() + solves);
        source_evaluations.set(
            source_evaluations
                .get()
                .checked_add(candidate.report().source_evaluations)
                .ok_or_else(|| invalid("source evaluation count overflows"))?,
        );
        let loads = candidate
            .lorentz_loads_for_velocity(&selected, &nodal, 0.0, physical.loads, &mut || {
                (*cancel.borrow_mut())()
            })
            .map_err(SpatialConductionKickError::from)?;
        let lorentz_wrench_body = loads
            .wrench(&Se3::identity(), [0.0; 3])
            .map_err(SpatialConductionKickError::from)?;
        if finite(
            lorentz_wrench_body.pairing(midpoint_twist_body)
                - candidate.report().field.mechanical_power_w(),
        )?
        .abs()
            > physical.field.power_tolerance_w
        {
            return Err(invalid(
                "native nonlinear wrench and original-cell work disagree",
            ));
        }
        let force = lorentz_wrench_body.to_array();
        let mut impulse = [0.0; 6];
        for i in 0..6 {
            impulse[i] = finite(duration_s * finite(force[i] + external[i])?)?;
        }
        poll(&mut || (*cancel.borrow_mut())())?;
        let mechanical =
            apply_free_floating_impulses(&model, initial, &[], &[], &[], &[wrench(impulse)])
                .map_err(SpatialConductionKickError::from)?;
        let after = mechanical.base_after.twist_body.to_array();
        let actual_midpoint: [f64; 6] = core::array::from_fn(|i| 0.5 * old[i] + 0.5 * after[i]);
        let residual = (0..6)
            .map(|i| -> Result<f64> {
                Ok(finite(
                    finite(midpoint[i] - actual_midpoint[i])? / scales[i],
                )?)
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(Trial {
            candidate,
            mechanical,
            midpoint_twist_body,
            actual_midpoint,
            lorentz_wrench_body,
            residual,
        })
    };
    let numerical = FallibleFiniteDifference::new(6, controls.difference_step, |x: &[f64]| {
        evaluate(x).map(|trial| trial.residual)
    })
    .map_err(numerical_error)?;
    let construction = NewtonKrylovState::new(&numerical, vec![0.0; 6], controls.newton);
    if let Some(error) = numerical.take_failure() {
        return Err(numerical_error(error));
    }
    let mut solver = construction.map_err(SpatialConductionResponseKickError::Newton)?;
    let solved = solver.run_controlled(&numerical, controls.maximum_newton_iterations, &mut || {
        (*cancel.borrow_mut())() || numerical.has_failed()
    });
    if let Some(error) = numerical.take_failure() {
        return Err(numerical_error(error));
    }
    poll(&mut || (*cancel.borrow_mut())())?;
    if solved.progress != SolverRunProgress::Complete {
        return Err(SpatialConductionResponseKickError::Interrupted(
            solved.progress,
        ));
    }
    if !solved.report.converged {
        return Err(SpatialConductionResponseKickError::NotConverged(
            solved.report,
        ));
    }
    let trial = evaluate(&solver.x)?;
    let midpoint = trial.midpoint_twist_body.to_array();
    let mut angular_residual: f64 = 0.0;
    let mut linear_residual: f64 = 0.0;
    for i in 0..6 {
        let residual = finite(trial.actual_midpoint[i] - midpoint[i])?.abs();
        if i < 3 {
            angular_residual = angular_residual.max(residual);
        } else {
            linear_residual = linear_residual.max(residual);
        }
    }
    let field = &trial.candidate.report().field;
    let supplied_power_w = trial
        .candidate
        .report()
        .terminals
        .as_ref()
        .map(|terminals| terminals.supplied_power_w)
        .unwrap_or_else(|| field.terminal_power_w());
    let terminal_work_j = finite(duration_s * supplied_power_w)?;
    let joule_heat_j = finite(duration_s * field.joule_power_w())?;
    let mechanical_work_j = finite(duration_s * field.mechanical_power_w())?;
    let external_work_j =
        finite(duration_s * external_wrench_body.pairing(twist(trial.actual_midpoint)))?;
    let mechanical_work_residual_j = finite(
        finite(trial.mechanical.kinetic_energy_change - external_work_j)? - mechanical_work_j,
    )?;
    let total_energy_residual_j = finite(
        finite(finite(trial.mechanical.kinetic_energy_change + joule_heat_j)? - terminal_work_j)?
            - external_work_j,
    )?;
    if (field.joule_power_w() > 0.0 && joule_heat_j == 0.0)
        || angular_residual > physical.angular_velocity_tolerance_rad_s
        || linear_residual > physical.linear_velocity_tolerance_m_s
        || trial.mechanical.kinetic_work_residual.abs() > physical.energy_tolerance_j
        || mechanical_work_residual_j.abs() > physical.energy_tolerance_j
        || total_energy_residual_j.abs() > physical.energy_tolerance_j
    {
        return Err(SpatialConductionKickError::Balance {
            angular_velocity_residual_rad_s: angular_residual,
            linear_velocity_residual_m_s: linear_residual,
            mechanical_work_residual_j,
            total_energy_residual_j,
        }
        .into());
    }
    let (problem, field) = trial.candidate.into_parts();
    poll(&mut || (*cancel.borrow_mut())())?;
    Ok(SpatialConductionResponseKickCandidate {
        problem,
        report: SpatialConductionResponseKick {
            mechanical: trial.mechanical,
            field,
            midpoint_twist_body: trial.midpoint_twist_body,
            lorentz_wrench_body: trial.lorentz_wrench_body,
            terminal_work_j,
            joule_heat_j,
            mechanical_work_j,
            external_work_j,
            angular_velocity_residual_rad_s: angular_residual,
            linear_velocity_residual_m_s: linear_residual,
            mechanical_work_residual_j,
            total_energy_residual_j,
            field_solves: field_solves.get(),
            source_evaluations: source_evaluations.get(),
            nonlinear: solved.report,
        },
    })
}

struct Trial {
    candidate: SpatialMovingConductionResponseCandidate,
    mechanical: FreeFloatingImpulseResponse,
    midpoint_twist_body: Twist,
    actual_midpoint: [f64; 6],
    lorentz_wrench_body: Wrench,
    residual: Vec<f64>,
}

fn invalid(reason: &'static str) -> SpatialConductionResponseKickError {
    super::invalid(reason).into()
}
fn numerical_error(
    error: FiniteDifferenceError<SpatialConductionResponseKickError>,
) -> SpatialConductionResponseKickError {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(reason) => invalid(reason),
    }
}

// Original geometry and density pins are retained. Only explicitly owned
// electrical coordinates may vary during a mechanical residual evaluation.
// A caller may refresh external T/phase first, but cannot change the fixed
// body's inertia by silently sampling a different density at that endpoint.
fn admit_mass(
    problem: &SpatialConductionResponse,
    seeds: &[ResolvedConductorState],
    mass: &TetrahedralMass<'_>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<()> {
    for (original, current) in mass.nodes_m().iter().zip(problem.positions_m()) {
        poll(cancelled)?;
        if original != current {
            return Err(invalid(
                "mass and nonlinear field use different original vertices",
            ));
        }
    }
    for (original, current) in mass.tetrahedra().iter().zip(&problem.mesh().tets) {
        poll(cancelled)?;
        if *original != current.map(|v| v as usize) {
            return Err(invalid(
                "mass and nonlinear field use different original tetrahedra",
            ));
        }
    }
    let electrical = |dims| {
        [
            ConductorMagnitude::ElectricField.quantity().dims(),
            ConductorMagnitude::CurrentDensity.quantity().dims(),
        ]
        .contains(&dims)
    };
    for (region, &cell) in mass
        .reference()
        .regions()
        .iter()
        .zip(mass.region_elements())
    {
        poll(cancelled)?;
        let source = &problem.sources()[cell];
        let seed = &seeds[cell];
        source
            .admit(seed)
            .map_err(SpatialConductionResponseError::from)?;
        let original = region.material();
        if original.material() != seed.material().material()
            || original.card_identity() != seed.material().card_identity()
        {
            return Err(invalid(
                "electrical and inertial cells must retain the same material card",
            ));
        }
        for property in original.properties() {
            let receipt = &property.answer().receipt;
            let claim = source
                .card()
                .claims()
                .claim(receipt.selected)
                .ok_or_else(|| invalid("original inertial density claim is absent"))?;
            let field_axis = |axis: &str| {
                source.axes().iter().any(|(name, _)| name == axis)
                    || [
                        claim.validity.axis_quantities().get(axis),
                        receipt.axis_quantities.get(axis),
                        seed.material().axis_quantities().get(axis),
                    ]
                    .into_iter()
                    .flatten()
                    .any(|quantity| electrical(quantity.dims()))
            };
            if claim
                .value
                .axes()
                .any(|(axis, dims)| electrical(dims) || field_axis(axis))
                || claim.validity.bounds().keys().any(|axis| field_axis(axis))
            {
                return Err(invalid(
                    "inertial source response and applicability must be independent of electrical fields",
                ));
            }
        }
        let refreshed = original.at_query_point(source.card(), &seed.material().query())?;
        if original
            .properties()
            .iter()
            .zip(refreshed.properties())
            .any(|(a, b)| a.value_si().to_bits() != b.value_si().to_bits())
        {
            return Err(invalid(
                "refreshed conductor context changes the original mass density",
            ));
        }
    }
    poll(cancelled)?;
    Ok(())
}
