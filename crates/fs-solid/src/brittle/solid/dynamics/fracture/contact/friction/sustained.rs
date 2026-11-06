//! Static directional-friction load transfer through the native midpoint.
//!
//! Held normal and tangent forces enter one original nonlinear material/mass
//! response. The shared unregularized force law resolves its complete coupled
//! response in work-conjugate elliptic coordinates. Only matched, resting pairs
//! may carry load for this interval; finite sliding needs spatial contact
//! geometry and is refused instead of creating a spurious internal couple.

use fs_contact::sustained::elliptic::EllipticForceCoordinates;
use fs_contact::sustained::solve::{MAX_CONTACTS, SolveReport};
use fs_contact::sustained::{ContactLawError, ContactPoint, NaturalMapScales, linearize};
use fs_exec::Cx;

use super::{
    BrittleFaceFriction, BrittleFaceFrictionError, BrittleFaceFrictionFrame, EllipticFriction,
};
use crate::brittle::solid::dynamics::fracture::contact::BrittleFaceImpactRow;
use crate::brittle::solid::dynamics::fracture::contact::approach::ReleasedFaceMotionView;
use crate::brittle::solid::dynamics::fracture::contact::sustained::{
    BrittleFaceSustainedControls, BrittleFaceSustainedError, NativeFaceForceResponse,
    NativeFaceSustainedTrial, Pair, SustainedContactWork, SustainedWorkAdmission,
    admit_released_face_sustained_work, balance, pair_view, solve_released_face_force_response,
    sustained_external_loads, validate_controls,
};
use crate::brittle::solid::dynamics::midpoint::BrittleTetMidpointEnergy;
use crate::brittle::solid::dynamics::{BrittleTetDynamicError, BrittleTetMotion};
use crate::brittle::solid::fracture::require_released_faces;
use crate::brittle::solid::path::Sum;
use crate::brittle::solid::{BrittleTetError, checkpoint, energy_gate};
use crate::plastic3::IncrementSettings;

/// Original normal/midpoint budgets plus static tangential admission.
#[derive(Clone, Copy, Debug)]
pub struct BrittleFaceSustainedFrictionControls {
    pub contact: BrittleFaceSustainedControls,
    /// Positive numerical scale in the native solver's tangent coordinates.
    /// This changes conditioning only, not physical stiffness or stored energy.
    pub tangent_scale_n_per_m: f64,
    /// Maximum physical two-tangent displacement during a loaded interval.
    /// Both actual endpoint relative velocities must also be at rest.
    pub slip_tolerance_m: f64,
}

/// Accepted native material/motion interval with physical and solver witnesses.
#[derive(Debug)]
pub struct BrittleFaceSustainedFrictionResult<Motion, Energy> {
    pub motion: Motion,
    /// Original external/support, storage and fracture accounts, with internal
    /// contact loads removed from dead-load work. Static contact supplies no heat.
    pub energy: Energy,
    /// Unmodified midpoint account including the applied contact nodal loads.
    pub native_energy: Energy,
    /// Original native solve in explicitly TRANSFORMED force/path coordinates.
    /// A converged local force representative does not prove global uniqueness.
    pub solve: SolveReport,
    pub rows: Vec<BrittleFaceImpactRow>,
    pub frames: Vec<BrittleFaceFrictionFrame>,
    /// Physical normal/tangent forces in the returned frames and row order.
    pub forces_n: Vec<[f64; 3]>,
    /// Actual endpoint gap and physical two-tangent displacement increment [m].
    pub paths_m: Vec<[f64; 3]>,
    pub coordinates: Vec<EllipticForceCoordinates>,
    pub normal_contact_work_j: f64,
    pub tangential_contact_work_j: f64,
    /// Sum of absolute normal and individual tangent work terms for every pair.
    pub absolute_contact_work_j: f64,
    /// Sum of individual absolute physical/solver coordinate-work differences.
    pub coordinate_work_error_j: f64,
    pub contact_work_balance_error_j: f64,
    /// Sum of each pair's larger initial/final absolute force couple, bounding
    /// its affine trace. Opposing unphysical couples cannot cancel this gate.
    pub absolute_contact_couple_n_m: f64,
    pub numerical_tangent_evaluations: usize,
}

pub type BrittleFaceSustainedFrictionStep<'a> =
    BrittleFaceSustainedFrictionResult<BrittleTetMotion<'a>, BrittleTetMidpointEnergy>;

pub(crate) type NativeFaceSustainedFrictionStep<Motion, Energy> =
    BrittleFaceSustainedFrictionResult<Motion, Energy>;

#[derive(Debug)]
pub enum BrittleFaceSustainedFrictionError {
    Sustained(BrittleFaceSustainedError),
    Frame(BrittleFaceFrictionError),
    Law {
        row: BrittleFaceImpactRow,
        error: ContactLawError,
    },
}
impl core::fmt::Display for BrittleFaceSustainedFrictionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "released brittle face static-friction interval refused: {self:?}"
        )
    }
}
impl std::error::Error for BrittleFaceSustainedFrictionError {}
impl From<BrittleFaceSustainedError> for BrittleFaceSustainedFrictionError {
    fn from(error: BrittleFaceSustainedError) -> Self {
        Self::Sustained(error)
    }
}
impl From<BrittleFaceFrictionError> for BrittleFaceSustainedFrictionError {
    fn from(error: BrittleFaceFrictionError) -> Self {
        Self::Frame(error)
    }
}
impl From<BrittleTetError> for BrittleFaceSustainedFrictionError {
    fn from(error: BrittleTetError) -> Self {
        Self::Sustained(error.into())
    }
}

impl<'a> BrittleTetMotion<'a> {
    /// Carry simultaneous normal and static frictional loads for one native
    /// physical interval, or continue unloaded separation of released pairs.
    ///
    /// Friction inputs follow retained face/vertex order. Two positive principal
    /// coefficients, or equal coefficients including zero, use the existing
    /// native ellipse. A single zero principal coefficient is explicitly refused
    /// by that sustained law; no artificial resistance fills its free direction.
    ///
    /// All non-tied pairs participate. Only initially active, matched pairs at
    /// rest can carry load. Loaded pairs must remain matched and at rest in all
    /// three relative endpoint velocity components, with negligible physical
    /// slip throughout the native affine displacement trace. An impact,
    /// open-to-loaded arrival or finite sliding requires another resolved mode.
    /// Unloaded pairs can separate without evaluating a friction law.
    ///
    /// The full nonlinear midpoint starts from the same accepted native history
    /// for every force trial, including all normal/shear and same-body coupling.
    /// After fresh force replay, physical cone, path, work and individual-couple
    /// checks gate one complete successor. No gap, position or velocity is
    /// projected, no penalty/compliance is added, and no contact work is assigned
    /// to heat. A refusal leaves the accepted body and step allowance unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_released_face_sustained_frictional(
        &self,
        dt_s: f64,
        friction: &[BrittleFaceFriction],
        prescribed_m: &[(usize, f64)],
        mean_nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        controls: BrittleFaceSustainedFrictionControls,
        cx: &Cx<'_>,
    ) -> Result<BrittleFaceSustainedFrictionStep<'a>, BrittleFaceSustainedFrictionError> {
        checkpoint(cx)?;
        if friction.len() != self.solid.released_faces.len() {
            return Err(BrittleFaceSustainedError::Invalid(
                "complete retained-face friction inputs required",
            )
            .into());
        }
        let directions: Vec<_> = friction.iter().map(|f| f.tangent_direction).collect();
        solve_released_face_sustained_frictional::<_, BrittleFaceSustainedFrictionError>(
            self.released_face_motion_view(),
            &directions,
            prescribed_m,
            mean_nodal_forces_n,
            controls,
            cx,
            |face, row, _, _, _| Ok(friction[face].vertices[row.vertex]),
            |total| {
                self.midpoint_candidate(
                    dt_s,
                    prescribed_m,
                    total,
                    settings,
                    controls.contact.work,
                    cx,
                )
            },
        )
    }
}

/// Shared native-body kernel. Resolve only initially active, matched, static
/// rows; open/separating/sliding rows keep zero friction and cannot acquire load
/// on acceptance. Source candidates stay private until this whole result passes.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(crate) fn solve_released_face_sustained_frictional<
    Step: NativeFaceSustainedTrial,
    E: From<BrittleFaceSustainedFrictionError>,
>(
    old_motion: ReleasedFaceMotionView<'_>,
    tangent_directions: &[[f64; 3]],
    prescribed_m: &[(usize, f64)],
    mean_nodal_forces_n: &[(usize, f64)],
    controls: BrittleFaceSustainedFrictionControls,
    cx: &Cx<'_>,
    mut resolve: impl FnMut(
        usize,
        &BrittleFaceImpactRow,
        &BrittleFaceFrictionFrame,
        [f64; 3],
        f64,
    ) -> Result<EllipticFriction, E>,
    candidate: impl FnMut(&[(usize, f64)]) -> Result<Step, BrittleTetDynamicError>,
) -> Result<NativeFaceSustainedFrictionStep<Step::Motion, Step::Energy>, E> {
    use BrittleFaceSustainedError as S;
    use BrittleFaceSustainedFrictionError as F;
    checkpoint(cx).map_err(F::from)?;
    validate_controls(controls.contact).map_err(F::from)?;
    if tangent_directions.len() != old_motion.released_faces.len()
        || !controls.tangent_scale_n_per_m.is_finite()
        || controls.tangent_scale_n_per_m <= 0.0
        || !controls.slip_tolerance_m.is_finite()
        || controls.slip_tolerance_m < 0.0
    {
        return Err(F::from(S::Invalid(
            "complete frame inputs, positive tangent scale and nonnegative slip gate required",
        ))
        .into());
    }
    let external = sustained_external_loads(
        old_motion.displacement_m.as_flattened().len(),
        prescribed_m,
        mean_nodal_forces_n,
        cx,
    )
    .map_err(F::from)?;
    let mut rows = Vec::new();
    let mut frames = Vec::new();
    let mut initial = Vec::new();
    let mut eligible = Vec::new();
    let mut coordinates = Vec::new();
    for (face_index, face) in old_motion.released_faces.iter().enumerate() {
        checkpoint(cx).map_err(F::from)?;
        let normal = super::unit(face.normal).map_err(F::from)?;
        let mut frame = None;
        for vertex in 0..3 {
            let nodes = [face.paired_nodes[0][vertex], face.paired_nodes[1][vertex]];
            if nodes[0] == nodes[1] {
                continue;
            }
            if rows.len() >= MAX_CONTACTS || rows.len() >= controls.contact.solve.maximum_contacts {
                return Err(F::from(S::Invalid(
                    "all released-face frames must fit the sustained solve budget",
                ))
                .into());
            }
            let row = BrittleFaceImpactRow {
                face_id: face.id,
                vertex,
                nodes,
            };
            let geometry = pair_view(old_motion, row, normal).map_err(F::from)?;
            if geometry.gap < 0.0 {
                return Err(F::from(S::Geometry {
                    face: face.id,
                    vertex,
                    reason: "negative initial gap cannot be repaired",
                })
                .into());
            }
            if geometry.gap <= controls.contact.activation_tolerance_m
                && geometry.speed < -controls.contact.velocity_tolerance_mps
            {
                return Err(F::from(S::ImpactRequired {
                    face: face.id,
                    vertex,
                    normal_speed_mps: geometry.speed,
                })
                .into());
            }
            let frame = match frame {
                Some(frame) => frame,
                None => {
                    let direction = super::unit(tangent_directions[face_index]).map_err(F::from)?;
                    let tangent2 = super::unit(super::cross(normal, direction)).map_err(F::from)?;
                    let tangent1 = super::unit(super::cross(tangent2, normal)).map_err(F::from)?;
                    let value = BrittleFaceFrictionFrame {
                        normal,
                        tangent1,
                        tangent2,
                    };
                    frame = Some(value);
                    value
                }
            };
            let velocity = super::velocity(old_motion.velocity_mps, row, frame).map_err(F::from)?;
            let active = geometry.gap <= controls.contact.activation_tolerance_m
                && geometry.tangent <= controls.contact.tangential_tolerance_m
                && velocity
                    .iter()
                    .all(|v| v.abs() <= controls.contact.velocity_tolerance_mps);
            let friction = if active {
                resolve(face_index, &row, &frame, velocity, geometry.gap)?
            } else {
                // This row can only be accepted unloaded. Its numerical normal
                // response still detects unexpected load or open-to-loaded arrival.
                EllipticFriction {
                    coefficients_mu: [0.0; 2],
                    principal_axis: [1.0, 0.0],
                }
            };
            coordinates.push(
                EllipticForceCoordinates::prepare(friction)
                    .map_err(|error| F::Law { row, error })?,
            );
            rows.push(row);
            frames.push(frame);
            initial.push(geometry);
            eligible.push(active);
        }
    }
    if rows.is_empty() {
        return Err(F::from(S::Invalid("no non-tied released-face contact rows")).into());
    }
    let friction: Vec<_> = coordinates.iter().map(|c| c.coefficient()).collect();
    let scales = vec![
        NaturalMapScales {
            normal_n_per_m: controls.contact.normal_scale_n_per_m,
            tangent_n_per_m: controls.tangent_scale_n_per_m,
        };
        rows.len()
    ];
    let seeds = vec![[0.0; 3]; rows.len()];
    let NativeFaceForceResponse {
        step,
        solve,
        contact_load_n,
        numerical_tangent_evaluations,
    } = solve_released_face_force_response(
        &external,
        &rows,
        &friction,
        &scales,
        &seeds,
        controls.contact.solve,
        cx,
        candidate,
        |i, q| {
            let force = coordinates[i]
                .force_to_physical(q)
                .map_err(response_law_error)?;
            spatial_force(frames[i], force)
        },
        |i, motion| {
            coordinates[i]
                .path_to_solver(physical_path(motion, rows[i], frames[i], initial[i])?)
                .map_err(response_law_error)
        },
    )
    .map_err(F::from)?;
    let final_motion = step.face_motion();
    let mut forces_n = Vec::with_capacity(rows.len());
    let mut paths_m = Vec::with_capacity(rows.len());
    let (mut normal_work, mut tangent_work, mut absolute_work, mut coordinate_work, mut couple) = (
        Sum::default(),
        Sum::default(),
        Sum::default(),
        Sum::default(),
        Sum::default(),
    );
    for (i, &row) in rows.iter().enumerate() {
        checkpoint(cx).map_err(F::from)?;
        let frame = frames[i];
        let before = initial[i];
        let after = pair_view(final_motion, row, frame.normal).map_err(F::from)?;
        let path = physical_path(final_motion, row, frame, before)
            .map_err(|error| F::from(S::Motion(error)))?;
        let q = solve.forces_n[i];
        let force = coordinates[i]
            .force_to_physical(q)
            .map_err(|error| F::Law { row, error })?;
        let mapped_path = coordinates[i]
            .path_to_solver(path)
            .map_err(|error| F::Law { row, error })?;
        if force[0] < 0.0
            || solve.paths_m[i] != mapped_path
            || (friction[i] == 0.0 && (force[1] != 0.0 || force[2] != 0.0))
        {
            return Err(F::from(S::Invalid(
                "final force/path identity or zero-friction traction refused",
            ))
            .into());
        }
        if after.gap < 0.0 {
            return Err(F::from(S::Geometry {
                face: row.face_id,
                vertex: row.vertex,
                reason: "negative final gap cannot be repaired",
            })
            .into());
        }
        let after_velocity =
            super::velocity(final_motion.velocity_mps, row, frame).map_err(F::from)?;
        let loaded = super::finite(super::norm(force)).map_err(F::from)?
            > controls.contact.force_tolerance_n;
        if loaded
            && (!eligible[i]
                || after.gap > controls.contact.activation_tolerance_m
                || after.tangent > controls.contact.tangential_tolerance_m
                || after_velocity
                    .iter()
                    .any(|v| v.abs() > controls.contact.velocity_tolerance_mps)
                || path[1].hypot(path[2]) > controls.slip_tolerance_m)
        {
            return Err(F::from(S::Mode { face: row.face_id, vertex: row.vertex,
                reason: "loaded friction requires matched static endpoints; sliding or arrival needs spatial contact geometry" }).into());
        }
        // Re-admit the actual physical force after its coordinate round trip.
        // Physical slip/speeds above remain independently bounded, so a weak
        // principal axis cannot conceal a moving point in solver coordinates.
        let measured_force = coordinates[i]
            .force_to_solver(force)
            .map_err(|error| F::Law { row, error })?;
        let law = linearize(
            ContactPoint {
                force_n: measured_force,
                gap_m: after.gap,
                slip_increment_m: [mapped_path[1], mapped_path[2]],
                friction_coefficient: friction[i],
            },
            scales[i],
        )
        .map_err(|error| F::Law { row, error })?;
        for defect in law.residual_n {
            balance(
                "actual sustained ellipse force residual",
                defect,
                controls.contact.solve.force_residual_tolerance_n,
            )
            .map_err(F::from)?;
        }
        balance(
            "actual sustained ellipse force capacity",
            (measured_force[1].hypot(measured_force[2]) - friction[i] * force[0]).max(0.0),
            controls.contact.force_tolerance_n,
        )
        .map_err(F::from)?;
        for (account, gap) in [
            ("sustained friction initial complementarity", before.gap),
            ("sustained friction final complementarity", after.gap),
        ] {
            energy_gate(
                account,
                force[0] * gap,
                controls.contact.complementarity_tolerance_j,
            )
            .map_err(F::from)?;
        }
        let normal = force[0] * (after.gap - before.gap);
        let tangent = [force[1] * path[1], force[2] * path[2]];
        normal_work.add(normal);
        tangent_work.add(tangent[0]);
        tangent_work.add(tangent[1]);
        absolute_work.add(normal.abs());
        absolute_work.add(tangent[0].abs());
        absolute_work.add(tangent[1].abs());
        let mut actual_work = Sum::default();
        let mut mapped_work = Sum::default();
        for value in [normal, tangent[0], tangent[1]] {
            actual_work.add(value);
        }
        for value in [
            q[0] * (after.gap - before.gap),
            q[1] * mapped_path[1],
            q[2] * mapped_path[2],
        ] {
            mapped_work.add(value);
        }
        coordinate_work.add((actual_work.value() - mapped_work.value()).abs());
        let world = spatial_force(frame, force).map_err(|error| F::from(S::Motion(error)))?;
        let start_couple =
            super::finite(super::norm(super::cross(before.separation, world))).map_err(F::from)?;
        let end_couple =
            super::finite(super::norm(super::cross(after.separation, world))).map_err(F::from)?;
        couple.add(start_couple.max(end_couple));
        forces_n.push(force);
        paths_m.push(path);
    }
    require_released_faces(
        old_motion.displacement_m.len(),
        old_motion.released_faces,
        final_motion.displacement_m,
        cx,
    )
    .map_err(F::from)?;
    let normal_contact_work_j = normal_work.value();
    let tangential_contact_work_j = tangent_work.value();
    let absolute_contact_work_j = absolute_work.value();
    let coordinate_work_error_j = coordinate_work.value();
    let absolute_contact_couple_n_m = couple.value();
    energy_gate(
        "sustained friction coordinate work",
        coordinate_work_error_j,
        controls.contact.contact_work_tolerance_j,
    )
    .map_err(F::from)?;
    let SustainedWorkAdmission {
        energy,
        native_energy,
        contact_work_balance_error_j,
    } = admit_released_face_sustained_work(
        old_motion,
        &step,
        mean_nodal_forces_n,
        &contact_load_n,
        SustainedContactWork {
            contact_work_j: normal_contact_work_j + tangential_contact_work_j,
            absolute_contact_work_j,
            absolute_contact_couple_n_m,
        },
        controls.contact,
        cx,
    )
    .map_err(F::from)?;
    checkpoint(cx).map_err(F::from)?;
    Ok(BrittleFaceSustainedFrictionResult {
        motion: step.into_motion(),
        energy,
        native_energy,
        solve,
        rows,
        frames,
        forces_n,
        paths_m,
        coordinates,
        normal_contact_work_j,
        tangential_contact_work_j,
        absolute_contact_work_j,
        coordinate_work_error_j,
        contact_work_balance_error_j,
        absolute_contact_couple_n_m,
        numerical_tangent_evaluations,
    })
}

fn response_law_error(error: ContactLawError) -> BrittleTetDynamicError {
    match error {
        ContactLawError::Invalid(reason) | ContactLawError::NonFinite(reason) => {
            BrittleTetDynamicError::Invalid(reason)
        }
    }
}

fn spatial_force(
    frame: BrittleFaceFrictionFrame,
    force: [f64; 3],
) -> Result<[f64; 3], BrittleTetDynamicError> {
    let axes = frame.axes();
    let mut result = [0.0; 3];
    for axis in 0..3 {
        result[axis] = super::product([axes[0][axis], axes[1][axis], axes[2][axis]], force)
            .map_err(|_| {
                BrittleTetDynamicError::Invalid("unrepresentable physical contact force")
            })?;
    }
    Ok(result)
}

fn physical_path(
    motion: ReleasedFaceMotionView<'_>,
    row: BrittleFaceImpactRow,
    frame: BrittleFaceFrictionFrame,
    before: Pair,
) -> Result<[f64; 3], BrittleTetDynamicError> {
    let separation = core::array::from_fn(|axis| {
        motion.displacement_m[row.nodes[1]][axis] - motion.displacement_m[row.nodes[0]][axis]
    });
    let increment = core::array::from_fn(|axis| separation[axis] - before.separation[axis]);
    let gap = super::super::dot(frame.normal, separation);
    if !gap.is_finite() {
        return Err(BrittleTetDynamicError::Invalid(
            "unrepresentable paired response gap",
        ));
    }
    Ok([
        gap,
        super::product(frame.tangent1, increment)
            .map_err(|_| BrittleTetDynamicError::Invalid("unrepresentable paired tangent path"))?,
        super::product(frame.tangent2, increment)
            .map_err(|_| BrittleTetDynamicError::Invalid("unrepresentable paired tangent path"))?,
    ])
}
