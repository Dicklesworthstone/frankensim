//! Two-axis retained shear coupled to the original articulated mass response.
//!
//! The two world-frame tangent rows share one native ABA workspace. Both
//! diagonal and cross terms of W = J M^-1 J^T enter C = dt^2 W / 2; the
//! shared constitutive embedding resolves du = -dt v_before - C F. One
//! accumulated impulse dt F then updates the original articulated systems.
//! Actual endpoint velocities independently recover the work-conjugate slip.
//!
//! This is one frozen-configuration shear substep. The original tangent frame
//! and the constitutive checkpoint must describe the same physical axes on
//! every call; no history rotation, contact discovery or normal-force solve is
//! inferred. In particular, a scalar Hertz--Mindlin law keeps its own refusal
//! of transverse motion. A two-axis layer or dry law supplies its own physics.

use super::{ArticulatedTangentialConfig, ArticulatedTangentialDrive, admit};
use crate::articulated::FreeFloatingImpulseResponse;
use crate::articulated_contacts::{
    ArticulatedContactError, ArticulatedContactSystem, ContactCounterpart, EventResponse,
    EventWorkspace, LinkContactPoint, Prepared, finite, magnitude, poll, prepare_point, vector,
};
use fs_ga::Vec3;
use fs_math::det;
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::tangential::matrix::{MatrixTangentialConfig, advance_with_matrix_compliance};
use fs_tribo::{TangentialLaw, TangentialStep, TangentialStepInput};

/// One common interface point with the original validated physical shear frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ArticulatedPlanarTangentialContact {
    /// Positive-effort articulated link material point.
    pub point: LinkContactPoint,
    /// Opposite dynamic point or translating support at the same world point.
    pub other: ContactCounterpart,
    /// Existing orthonormal world-frame tangent basis, in constitutive order.
    /// Accepted history must retain this physical frame association.
    pub frame_world: TangentFrame,
}

/// Mechanical admission and the independent bounded two-axis coupling solve.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ArticulatedPlanarTangentialConfig {
    /// Original force, time, speed, common-point and actual-work allowances.
    pub mechanical: ArticulatedTangentialConfig,
    /// Full matrix admission, constitutive-evaluation and compatibility budgets.
    pub coupling: MatrixTangentialConfig,
}

/// Private coupled candidate; publish all body responses and history together.
#[derive(Clone, Debug)]
pub struct ArticulatedPlanarTangentialResult<S> {
    /// Original ABA updates, in the supplied system order.
    pub responses: Vec<FreeFloatingImpulseResponse>,
    /// Original constitutive state, work-conjugate force and physical heat.
    pub contact: TangentialStep<S>,
    /// Full native contact mobility [1/kg], including both cross terms.
    pub mobility_kg_inverse: [[f64; 2]; 2],
    /// Held-force displacement response dt^2 W / 2 [m/N].
    pub compliance_m_per_n: [[f64; 2]; 2],
    /// Original point-minus-counterpart relative tangent speeds [m/s].
    pub relative_speed_before_m_s: [f64; 2],
    /// Speeds recomputed from the actual native endpoint states [m/s].
    pub relative_speed_after_m_s: [f64; 2],
    /// Actual counterpart-minus-point midpoint relative displacement [m].
    pub slip_increment_m: [f64; 2],
    /// Actual displacement minus the retained constitutive increment [m].
    pub compatibility_residual_m: [f64; 2],
    /// Outer matrix-coupling iterations and original-law evaluation count.
    pub coupling_iterations: usize,
    pub law_evaluations: usize,
    /// Original reversible contact energy before and after [J].
    pub contact_energy_before_j: f64,
    pub contact_energy_after_j: f64,
    /// Actual total articulated kinetic-energy change [J].
    pub kinetic_energy_change_j: f64,
    /// Signed work supplied by translating supports [J].
    pub support_work_j: f64,
    /// Contact-storage change + heat - F dot actual slip [J].
    pub contact_work_residual_j: f64,
    /// Actual kinetic + contact storage change + heat - support work [J].
    pub energy_balance_residual_j: f64,
    /// Actual world momentum changes after external support impulse [N s, N m s].
    pub linear_momentum_residual: Vec3,
    pub angular_momentum_residual: Vec3,
    /// Actual versus projected native point velocity discrepancy [m/s].
    pub endpoint_velocity_residual_m_s: f64,
}

/// Advance retained two-axis contact on at most 128 original articulated trees.
///
/// Both tangent impulses are solved together and applied together. The complete
/// compliance is passed to the shared `TangentialLaw` embedding; neither its
/// cross terms nor a transverse stored force are dropped. A prescribed support
/// supplies its signed physical work. Two dynamic counterparts exchange equal
/// and opposite impulses at the admitted common point.
///
/// The call only borrows accepted mechanics and constitutive history. Every
/// fallible solve, endpoint check and cancellation poll precedes return of the
/// complete candidate. The enclosing integrator owns pose/normal evolution,
/// physical frame transport and refinement across intrastep reversals. Dry
/// sticking is a condition on midpoint slip, not on endpoint velocity. Multiple
/// shear contacts and finite-patch spin/rolling are not inferred by this call.
pub fn advance_articulated_planar_tangential<L: TangentialLaw>(
    systems: &[ArticulatedContactSystem<'_>],
    contact: ArticulatedPlanarTangentialContact,
    law: &L,
    state: &L::State,
    drive: ArticulatedTangentialDrive,
    config: ArticulatedPlanarTangentialConfig,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ArticulatedPlanarTangentialResult<L::State>, ArticulatedContactError> {
    poll(&mut cancelled)?;
    let mechanical = config.mechanical;
    admit(systems.len(), drive, mechanical)?;
    config.coupling.validate()?;
    let before_energy = finite(law.stored_energy_j(state))?;
    if before_energy < 0.0 {
        return Err(ArticulatedContactError::Invalid("negative contact storage"));
    }
    let axes = admit_frame(contact.frame_world)?;
    let mut workspace = EventWorkspace::new(systems, &mut cancelled)?;
    let rows = [
        prepare_axis(
            contact,
            axes[0],
            &workspace,
            mechanical.geometry_tolerance_m,
        )?,
        prepare_axis(
            contact,
            axes[1],
            &workspace,
            mechanical.geometry_tolerance_m,
        )?,
    ];
    let before_speed = [
        finite(rows[0].velocity(&workspace.before))?,
        finite(rows[1].velocity(&workspace.before))?,
    ];
    let native = workspace.delassus(systems, &rows, &mut cancelled)?;
    let mobility = [[native[0], native[1]], [native[2], native[3]]];
    let mut compliance = [[0.0; 2]; 2];
    for i in 0..2 {
        for j in 0..2 {
            compliance[i][j] = finite(0.5 * drive.dt_s * drive.dt_s * mobility[i][j])?;
        }
    }
    let coupled = advance_with_matrix_compliance(
        law,
        state,
        &TangentialStepInput {
            normal_n: drive.normal_force_n,
            slip_increment_m: [
                finite(-drive.dt_s * before_speed[0])?,
                finite(-drive.dt_s * before_speed[1])?,
            ],
            dt_s: drive.dt_s,
            contact_temperature_k: drive.contact_temperature_k,
            nominal_pressure_pa: drive.nominal_pressure_pa,
            compliance_m_per_n: 0.0,
            solve: mechanical.solve,
        },
        compliance,
        config.coupling,
        &mut cancelled,
    )?;
    poll(&mut cancelled)?;
    let step = coupled.contact;
    let force = [finite(step.traction_n[0])?, finite(step.traction_n[1])?];
    if finite(det::hypot(force[0], force[1]))? > mechanical.solve.maximum_force_n
        || step.iterations > mechanical.solve.max_iterations
    {
        return Err(ArticulatedContactError::Unresolved(
            "constitutive effort/work cap exceeded",
        ));
    }
    for residual in step.compatibility_residual_m {
        if finite(residual)?.abs() > mechanical.displacement_tolerance_m {
            return Err(ArticulatedContactError::Unresolved(
                "constitutive two-axis shear compatibility",
            ));
        }
    }
    let after_energy = finite(law.stored_energy_j(&step.state))?;
    let stored_change = finite(after_energy - before_energy)?;
    let heat = finite(step.dissipated_j)?;
    if after_energy < 0.0 || heat < 0.0 {
        return Err(ArticulatedContactError::Invalid(
            "negative contact storage or heat",
        ));
    }
    for residual in [
        finite(stored_change - step.stored_change_j)?,
        finite(step.closure_residual_j)?,
        finite(coupled.body_work_residual_j)?,
    ] {
        if residual.abs() > mechanical.energy_tolerance_j {
            return Err(ArticulatedContactError::Unresolved(
                "constitutive two-axis contact energy",
            ));
        }
    }
    let impulse = [
        finite(drive.dt_s * force[0])?,
        finite(drive.dt_s * force[1])?,
    ];
    let predicted_speed = [
        finite(before_speed[0] + mobility[0][0] * impulse[0] + mobility[0][1] * impulse[1])?,
        finite(before_speed[1] + mobility[1][0] * impulse[0] + mobility[1][1] * impulse[1])?,
    ];
    let mechanics = workspace.apply(
        systems,
        &rows,
        EventResponse {
            impulses: &impulse,
            predicted_speeds: &predicted_speed,
            stored_energy_change: stored_change,
            dissipated_energy: heat,
            velocity_tolerance: mechanical.velocity_tolerance_m_s,
            energy_tolerance: mechanical.energy_tolerance_j,
        },
        &mut cancelled,
    )?;
    let after_speed = [mechanics.speeds_after[0], mechanics.speeds_after[1]];
    let mut slip = [0.0; 2];
    let mut compatibility = [0.0; 2];
    for i in 0..2 {
        slip[i] = finite(-drive.dt_s * (0.5 * before_speed[i] + 0.5 * after_speed[i]))?;
        compatibility[i] = finite(slip[i] - step.slip_increment_m[i])?;
        if compatibility[i].abs() > mechanical.displacement_tolerance_m {
            return Err(ArticulatedContactError::Unresolved(
                "actual two-axis shear displacement",
            ));
        }
    }
    let work_residual = finite(stored_change + heat - force[0] * slip[0] - force[1] * slip[1])?;
    if finite(det::hypot(slip[0], slip[1]) / drive.dt_s)? > mechanical.maximum_slip_speed_m_s
        || work_residual.abs() > mechanical.energy_tolerance_j
    {
        return Err(ArticulatedContactError::Unresolved(
            "actual two-axis shear speed/work",
        ));
    }
    poll(&mut cancelled)?;
    Ok(ArticulatedPlanarTangentialResult {
        responses: mechanics.responses,
        contact: step,
        mobility_kg_inverse: mobility,
        compliance_m_per_n: compliance,
        relative_speed_before_m_s: before_speed,
        relative_speed_after_m_s: after_speed,
        slip_increment_m: slip,
        compatibility_residual_m: compatibility,
        coupling_iterations: coupled.iterations,
        law_evaluations: coupled.law_evaluations,
        contact_energy_before_j: before_energy,
        contact_energy_after_j: after_energy,
        kinetic_energy_change_j: mechanics.kinetic_energy_change,
        support_work_j: mechanics.plane_work,
        contact_work_residual_j: work_residual,
        energy_balance_residual_j: mechanics.energy_balance_residual,
        linear_momentum_residual: mechanics.linear_momentum_residual,
        angular_momentum_residual: mechanics.angular_momentum_residual,
        endpoint_velocity_residual_m_s: mechanics.endpoint_velocity_residual,
    })
}

fn admit_frame(frame: TangentFrame) -> Result<[Vec3; 2], ArticulatedContactError> {
    let basis = [frame.longitudinal(), frame.lateral(), frame.normal()]
        .map(|v| Vec3::new(v[0], v[1], v[2]));
    // A nearly parallel constructor hint can lose orthogonality in its floating
    // projection. Admit arithmetic roundoff only; never reinterpret a skew or
    // partly normal coordinate as tangential work.
    for i in 0..3 {
        vector(basis[i])?;
        for j in 0..=i {
            let expected = if i == j { 1.0 } else { 0.0 };
            if (finite(basis[i].dot(basis[j]))? - expected).abs() > 256.0 * f64::EPSILON {
                return Err(ArticulatedContactError::Invalid(
                    "two-axis shear requires an orthonormal physical tangent frame",
                ));
            }
        }
    }
    Ok([basis[0], basis[1]])
}

fn prepare_axis(
    contact: ArticulatedPlanarTangentialContact,
    axis: Vec3,
    workspace: &EventWorkspace,
    geometry_tolerance_m: f64,
) -> Result<Prepared, ArticulatedContactError> {
    let (first, world_point) = prepare_point(contact.point, axis, &workspace.before)?;
    let (second, other_point, plane_speed) = match contact.other {
        ContactCounterpart::Body(point) => {
            if point.system_index == contact.point.system_index
                && point.link_index == contact.point.link_index
            {
                return Err(ArticulatedContactError::Invalid(
                    "two-axis shear requires distinct links",
                ));
            }
            let (term, position) = prepare_point(point, axis.scale(-1.0), &workspace.before)?;
            (Some(term), position, 0.0)
        }
        ContactCounterpart::PrescribedPlane {
            point_world,
            velocity_world,
        } => {
            vector(point_world)?;
            vector(velocity_world)?;
            (None, point_world, finite(axis.dot(velocity_world))?)
        }
    };
    if magnitude(world_point - other_point) > geometry_tolerance_m {
        return Err(ArticulatedContactError::Invalid(
            "two-axis shear material points do not coincide",
        ));
    }
    Ok(Prepared {
        first,
        second,
        axis_world: axis,
        world_point,
        plane_speed,
        minimum_speed: None,
    })
}
