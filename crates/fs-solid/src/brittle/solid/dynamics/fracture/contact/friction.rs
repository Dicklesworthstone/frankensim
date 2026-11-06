//! Simultaneous directional friction on existing completed-fracture witnesses.
//!
//! Sparse paired point covectors enter the original consistent mass. Every
//! normal/shear and intercontact block enters the native elliptic Coulomb law;
//! one combined impulse then enters the original body momentum/work audit.
//! The existing solid elliptic owner independently checks the physical law at
//! the actual endpoint. No configuration, material history or clock changes.

/// Original interface sources sampled at actual incoming relative slip.
pub mod source;

/// Static directional-friction loads through the original nonlinear midpoint.
pub mod sustained;

/// One original physical-time driver for directional impact, static dwell and release.
pub mod time;

use std::sync::Arc;

pub use fs_contact::frictional_multicontact::elliptic::EllipticFriction;
use fs_contact::frictional_multicontact::elliptic::{
    EllipticContactResiduals, EllipticFrictionalContact, EllipticFrictionalContactError,
    EllipticFrictionalContactResult, solve_elliptic_frictional_contacts,
};
use fs_contact::frictional_multicontact::{
    FrictionalContactConfig, MAX_FRICTIONAL_CONTACTS, NormalVelocityLaw,
};
use fs_exec::Cx;
use fs_math::dd::Dd;

use super::{
    BrittleFaceImpactRow, ReleasedFaceImpactBody, SolidImpulseControls, SolidImpulseReport, cross,
    loads, norm,
};
use crate::brittle::solid::dynamics::{BrittleTetDynamicError, BrittleTetMotion};
use crate::brittle::solid::{BrittleTetError, BrittleTetSolid, checkpoint, energy_gate};
use crate::dynamics::frictional_impact::SolidFrictionalImpactError;
use crate::dynamics::frictional_impact::elliptic::actual_endpoint_residuals;
use crate::plastic3::SmallStrainTetError;

/// Explicit constitutive input for one retained face, in retained face order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BrittleFaceFriction {
    /// A finite nonzero world direction with a nonzero projection into the
    /// retained face plane. That projection defines tangent 1; n cross t1
    /// defines tangent 2. This constructs coordinates, not a new position.
    pub tangent_direction: [f64; 3],
    /// Principal coefficients and principal axis in that SAME tangent frame,
    /// one for each retained paired vertex. These are explicit model inputs;
    /// this entry point does not assign them material-source authority.
    pub vertices: [EllipticFriction; 3],
}

/// Right-handed contact coordinates retained with each physical row.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BrittleFaceFrictionFrame {
    pub normal: [f64; 3],
    pub tangent1: [f64; 3],
    pub tangent2: [f64; 3],
}
impl BrittleFaceFrictionFrame {
    fn axes(self) -> [[f64; 3]; 3] {
        [self.normal, self.tangent1, self.tangent2]
    }

    /// Coordinate construction shared by fixed-configuration events and dwell.
    /// The native retained normal has already been normalized by the caller.
    pub(crate) fn from_unit_normal(
        normal: [f64; 3],
        direction: [f64; 3],
    ) -> Result<Self, BrittleFaceFrictionError> {
        let direction = unit(direction)?;
        let tangent2 = unit(cross(normal, direction))?;
        let tangent1 = unit(cross(tangent2, normal))?;
        Ok(Self {
            normal,
            tangent1,
            tangent2,
        })
    }
}

/// Native contact/mass budgets and independent physical assembly allowances.
#[derive(Clone, Copy, Debug)]
pub struct BrittleFaceFrictionControls {
    pub contact: FrictionalContactConfig,
    pub impulse: SolidImpulseControls,
    pub activation_tolerance_m: f64,
    /// Active paired witnesses must lie on the same retained normal line.
    pub tangential_tolerance_m: f64,
    /// Actual versus predicted contact-velocity allowance. The original law's
    /// own velocity, impulse, complementarity and work limits also apply.
    pub velocity_tolerance_mps: f64,
    pub energy_tolerance_j: f64,
    /// All active non-tied pairs must fit this cap AND the original solver's
    /// contact and scalar-row caps. Contact scratch is quadratic in 3*rows;
    /// one native body mass response is retained at a time.
    pub maximum_contact_rows: usize,
}

/// Complete accepted velocity successor and original physical law witnesses.
#[derive(Debug)]
pub struct BrittleFaceFrictionStep<'a> {
    pub motion: BrittleTetMotion<'a>,
    /// Impulses and velocities use the returned [normal, tangent1, tangent2]
    /// frames. The inner transformed proof retains its native mapped basis.
    /// Velocity-only determinacy does not identify unique point loads or heat.
    pub contact: EllipticFrictionalContactResult,
    pub impulse: SolidImpulseReport,
    pub rows: Vec<BrittleFaceImpactRow>,
    pub frames: Vec<BrittleFaceFrictionFrame>,
    pub measured_residuals: EllipticContactResiduals,
    pub maximum_velocity_residual_mps: f64,
    /// Total physical midpoint contact work measured at the actual endpoint.
    /// This is a combined normal/tangential loss, not a heat allocation.
    pub relative_work_j: f64,
    /// Actual midpoint contact work minus the original law's relative work.
    pub contact_work_error_j: f64,
    /// Native kinetic change minus support work minus law relative work.
    pub relative_work_error_j: f64,
    /// Native kinetic change minus support work plus law dissipation.
    pub energy_balance_error_j: f64,
    /// Sum of individual pair |separation cross physical impulse|. Opposing
    /// artificial couples cannot cancel this independent angular gate.
    pub absolute_contact_couple_kg_m2_s: f64,
}

#[derive(Debug)]
pub enum BrittleFaceFrictionError {
    Invalid(&'static str),
    Motion(BrittleTetDynamicError),
    Contact(EllipticFrictionalContactError),
    NativeAudit(SolidFrictionalImpactError),
    Geometry {
        face: u64,
        vertex: usize,
        reason: &'static str,
    },
    Response {
        face: u64,
        vertex: usize,
        component: usize,
        residual_mps: f64,
        tolerance_mps: f64,
    },
    InternalImpulse {
        linear_norm_kg_m_s: f64,
        angular_norm_kg_m2_s: f64,
        absolute_couple_kg_m2_s: f64,
    },
}
impl core::fmt::Display for BrittleFaceFrictionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "released brittle face frictional impact refused: {self:?}"
        )
    }
}
impl std::error::Error for BrittleFaceFrictionError {}
impl From<BrittleTetError> for BrittleFaceFrictionError {
    fn from(error: BrittleTetError) -> Self {
        Self::Motion(error.into())
    }
}
impl From<SmallStrainTetError> for BrittleFaceFrictionError {
    fn from(error: SmallStrainTetError) -> Self {
        Self::Motion(error.into())
    }
}

impl<'a> BrittleTetMotion<'a> {
    /// Apply one coupled directional-friction event to every active released
    /// pair. Both input slices follow `solid().released_faces()` order.
    ///
    /// Shared crack-front nodes are omitted; every other geometrically active
    /// pair enters, including initially separating pairs. Negative opening or
    /// tangentially separated active witnesses refuse. This fixed-configuration
    /// law does not discover contacts or supply finite-sliding contact geometry.
    /// A caller tangent direction defines the contact coordinates; the native
    /// ellipse retains two principal coefficients, including an exact zero
    /// free-slip coefficient. Equal coefficients recover circular Coulomb
    /// friction, and two zeros recover the frictionless impulse law.
    ///
    /// Held velocity components retain their actual prescribed velocity and
    /// their signed support work. Normal and tangential impulses are solved
    /// together with all consistent-mass cross terms, then applied once. The
    /// actual endpoint, native momentum, each pair's couple and whole-body work
    /// gate publication. Positions, reference matter, material/fracture history,
    /// temperature, physical time and step allowance stay with the accepted
    /// motion. Any refusal leaves it available for retry.
    pub fn impact_released_faces_frictional(
        &self,
        restitution: &[f64],
        friction: &[BrittleFaceFriction],
        held_velocity_dofs: &[usize],
        origin_m: [f64; 3],
        controls: BrittleFaceFrictionControls,
        cx: &Cx<'_>,
    ) -> Result<BrittleFaceFrictionStep<'a>, BrittleFaceFrictionError> {
        checkpoint(cx)?;
        if friction.len() != self.solid.released_faces.len() {
            return Err(BrittleFaceFrictionError::Invalid(
                "complete retained-face friction inputs required",
            ));
        }
        let directions: Vec<_> = friction.iter().map(|f| f.tangent_direction).collect();
        let response = released_face_frictional_impact_with::<BrittleFaceFrictionError>(
            ReleasedFaceImpactBody {
                reference_nodes_m: self.solid.geometry.nodes_m(),
                displacement_m: &self.solid.solution.displacement_m,
                released_faces: &self.solid.released_faces,
                mass: &self.mass,
                velocity_mps: &self.velocity_mps,
                kinetic_energy_j: self.kinetic_energy_j,
            },
            restitution,
            &directions,
            held_velocity_dofs,
            origin_m,
            controls,
            cx,
            |face, row, _, _| Ok(friction[face].vertices[row.vertex]),
        )?;
        BrittleFaceFrictionStep::from_response(self, response, cx)
    }

    /// Publish an already audited native velocity without duplicating its fields.
    pub(crate) fn friction_velocity_successor(
        &self,
        velocity_mps: Vec<[f64; 3]>,
        kinetic_after_j: f64,
        cx: &Cx<'_>,
    ) -> Result<Self, BrittleFaceFrictionError> {
        let solid = BrittleTetSolid {
            geometry: self.solid.geometry,
            materials: self.solid.materials,
            bands: Arc::clone(&self.solid.bands),
            released_faces: Arc::clone(&self.solid.released_faces),
            budget: self.solid.budget,
            solution: self.solid.solution.clone(),
            stored_energy_j: self.solid.stored_energy_j,
            fracture_dissipation_j: self.solid.fracture_dissipation_j,
        };
        checkpoint(cx)?;
        Ok(Self {
            solid,
            velocity_mps,
            reference_density_kg_m3: Arc::clone(&self.reference_density_kg_m3),
            mass: Arc::clone(&self.mass),
            kinetic_energy_j: kinetic_after_j,
            time_s: self.time_s,
            maximum_steps: self.maximum_steps,
            accepted_steps: self.accepted_steps,
        })
    }
}

/// Native event before a mechanical or thermal caller publishes its unchanged
/// material owner. Source descendants retain staged samples until this succeeds.
pub(crate) struct ReleasedFaceFrictionResponse {
    pub(crate) velocity_mps: Vec<[f64; 3]>,
    pub(crate) contact: EllipticFrictionalContactResult,
    pub(crate) impulse: SolidImpulseReport,
    pub(crate) rows: Vec<BrittleFaceImpactRow>,
    pub(crate) frames: Vec<BrittleFaceFrictionFrame>,
    pub(crate) measured_residuals: EllipticContactResiduals,
    pub(crate) maximum_velocity_residual_mps: f64,
    pub(crate) relative_work_j: f64,
    pub(crate) contact_work_error_j: f64,
    pub(crate) relative_work_error_j: f64,
    pub(crate) energy_balance_error_j: f64,
    pub(crate) absolute_contact_couple_kg_m2_s: f64,
}

impl<'a> BrittleFaceFrictionStep<'a> {
    pub(super) fn from_response(
        motion: &BrittleTetMotion<'a>,
        response: ReleasedFaceFrictionResponse,
        cx: &Cx<'_>,
    ) -> Result<Self, BrittleFaceFrictionError> {
        Ok(Self {
            motion: motion.friction_velocity_successor(
                response.velocity_mps,
                response.impulse.kinetic_after_j,
                cx,
            )?,
            contact: response.contact,
            impulse: response.impulse,
            rows: response.rows,
            frames: response.frames,
            measured_residuals: response.measured_residuals,
            maximum_velocity_residual_mps: response.maximum_velocity_residual_mps,
            relative_work_j: response.relative_work_j,
            contact_work_error_j: response.contact_work_error_j,
            relative_work_error_j: response.relative_work_error_j,
            energy_balance_error_j: response.energy_balance_error_j,
            absolute_contact_couple_kg_m2_s: response.absolute_contact_couple_kg_m2_s,
        })
    }
}

/// One geometry scan and one original mass assembly for explicit or sourced
/// friction. The callback sees ONLY active rows at their actual incoming
/// velocity, in the exact frame subsequently used by the contact owner.
pub(crate) fn released_face_frictional_impact_with<E: From<BrittleFaceFrictionError>>(
    body: ReleasedFaceImpactBody<'_>,
    restitution: &[f64],
    tangent_directions: &[[f64; 3]],
    held_velocity_dofs: &[usize],
    origin_m: [f64; 3],
    controls: BrittleFaceFrictionControls,
    cx: &Cx<'_>,
    mut resolve: impl FnMut(
        usize,
        &BrittleFaceImpactRow,
        &BrittleFaceFrictionFrame,
        [f64; 3],
    ) -> Result<EllipticFriction, E>,
) -> Result<ReleasedFaceFrictionResponse, E> {
    use BrittleFaceFrictionError as F;
    checkpoint(cx).map_err(F::from)?;
    if restitution.len() != body.released_faces.len()
        || tangent_directions.len() != body.released_faces.len()
        || restitution
            .iter()
            .any(|e| !e.is_finite() || !(0.0..=1.0).contains(e))
        || !(1..=MAX_FRICTIONAL_CONTACTS).contains(&controls.maximum_contact_rows)
        || !(1..=MAX_FRICTIONAL_CONTACTS).contains(&controls.contact.maximum_contacts)
        || [
            controls.activation_tolerance_m,
            controls.tangential_tolerance_m,
            controls.velocity_tolerance_mps,
            controls.energy_tolerance_j,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
        || origin_m.iter().any(|x| !x.is_finite())
    {
        return Err(F::Invalid(
            "complete face inputs, bounded contacts and finite nonnegative gates required",
        )
        .into());
    }
    let capacity = controls
        .maximum_contact_rows
        .min(controls.contact.maximum_contacts)
        .min(controls.contact.matrix.max_variables / 3);
    let mut rows = Vec::with_capacity(capacity);
    let mut frames = Vec::with_capacity(capacity);
    let mut separations = Vec::with_capacity(capacity);
    let mut contacts = Vec::with_capacity(capacity);
    for (face_index, face) in body.released_faces.iter().enumerate() {
        checkpoint(cx).map_err(F::from)?;
        let normal = unit(face.normal)?;
        let mut frame = None;
        for vertex in 0..3 {
            let [a, b] = [face.paired_nodes[0][vertex], face.paired_nodes[1][vertex]];
            if a == b {
                continue;
            }
            // Corresponding reference coordinates coincide by topology. Direct
            // displacement differences retain gaps at large global origins.
            let separation = core::array::from_fn(|axis| {
                body.displacement_m[b][axis] - body.displacement_m[a][axis]
            });
            let gap = product(normal, separation)?;
            if gap < 0.0 {
                return Err(F::Geometry {
                    face: face.id,
                    vertex,
                    reason: "negative paired opening; penetration cannot be repaired",
                }
                .into());
            }
            if gap > controls.activation_tolerance_m {
                continue;
            }
            let tangent = core::array::from_fn(|axis| separation[axis] - gap * normal[axis]);
            if finite(norm(tangent))? > controls.tangential_tolerance_m {
                return Err(F::Geometry {
                    face: face.id,
                    vertex,
                    reason: "active paired witnesses are tangentially separated",
                }
                .into());
            }
            if rows.len() == capacity {
                return Err(F::Invalid(
                    "all active released-face frames must fit contact and scalar-row budgets",
                )
                .into());
            }
            let frame = match frame {
                Some(frame) => frame,
                None => {
                    let value = BrittleFaceFrictionFrame::from_unit_normal(
                        normal,
                        tangent_directions[face_index],
                    )?;
                    frame = Some(value);
                    value
                }
            };
            let row = BrittleFaceImpactRow {
                face_id: face.id,
                vertex,
                nodes: [a, b],
            };
            let velocity = velocity(body.velocity_mps, row, frame)?;
            let friction = resolve(face_index, &row, &frame, velocity)?;
            let id = u64::try_from(face_index)
                .ok()
                .and_then(|face| face.checked_mul(3))
                .and_then(|face| face.checked_add(vertex as u64))
                .ok_or(F::Invalid("unrepresentable stable contact identity"))?;
            contacts.push(EllipticFrictionalContact {
                id,
                gap,
                activation_tolerance: controls.activation_tolerance_m,
                relative_velocity: velocity,
                normal_law: NormalVelocityLaw::NewtonRestitution {
                    coefficient: restitution[face_index],
                },
                friction,
            });
            rows.push(row);
            frames.push(frame);
            separations.push(separation);
        }
    }
    if rows.is_empty() {
        return Err(F::Invalid("no active released-face contact rows").into());
    }
    let dimension = 3 * rows.len();
    let mut mobility = vec![0.0; dimension * dimension];
    for (column, (&row, &frame)) in rows.iter().zip(&frames).enumerate() {
        for (component, axis) in frame.axes().into_iter().enumerate() {
            checkpoint(cx).map_err(F::from)?;
            let (response, _) = body
                .mass
                .free_velocity_response(
                    held_velocity_dofs,
                    &loads(row, axis, 1.0),
                    controls.impulse,
                    cx,
                )
                .map_err(SmallStrainTetError::from)
                .map_err(F::from)?;
            for i in 0..rows.len() {
                let response = velocity(&response, rows[i], frames[i])?;
                for a in 0..3 {
                    mobility[(3 * i + a) * dimension + 3 * column + component] = response[a];
                }
            }
        }
    }
    let contact =
        solve_elliptic_frictional_contacts(&contacts, &mobility, controls.contact, || {
            cx.checkpoint().is_err()
        })
        .map_err(F::Contact)?;
    let mut impulses = Vec::with_capacity(6 * rows.len());
    let mut couples = Dd::ZERO;
    for (i, &amount) in contact.impulses.iter().enumerate() {
        checkpoint(cx).map_err(F::from)?;
        let axes = frames[i].axes();
        let mut physical = [0.0; 3];
        for axis in 0..3 {
            physical[axis] = product([axes[0][axis], axes[1][axis], axes[2][axis]], amount)?;
        }
        impulses.extend_from_slice(&loads(rows[i], physical, 1.0));
        couples = couples + Dd::from_f64(finite(norm(cross(separations[i], physical)))?);
    }
    let absolute_contact_couple_kg_m2_s = rounded(couples)?;
    let mut positions = Vec::with_capacity(body.velocity_mps.len());
    for (node, (reference, displacement)) in body
        .reference_nodes_m
        .iter()
        .zip(body.displacement_m)
        .enumerate()
    {
        if node % 256 == 0 {
            checkpoint(cx).map_err(F::from)?;
        }
        let position = core::array::from_fn(|axis| reference[axis] + displacement[axis]);
        if position.iter().any(|x| !x.is_finite()) {
            return Err(F::Invalid("unrepresentable physical position").into());
        }
        positions.push(position);
    }
    let (velocity_mps, impulse) = body
        .mass
        .impulse_candidate(
            &positions,
            body.velocity_mps,
            held_velocity_dofs,
            &impulses,
            origin_m,
            controls.impulse,
            cx,
        )
        .map_err(SmallStrainTetError::from)
        .map_err(F::from)?;
    let linear = finite(norm(impulse.applied.linear_kg_m_s))?;
    let angular = finite(norm(impulse.applied.angular_kg_m2_s))?;
    if linear > controls.impulse.linear_tolerance_kg_m_s
        || angular > controls.impulse.angular_tolerance_kg_m2_s
        || absolute_contact_couple_kg_m2_s > controls.impulse.angular_tolerance_kg_m2_s
    {
        return Err(F::InternalImpulse {
            linear_norm_kg_m_s: linear,
            angular_norm_kg_m2_s: angular,
            absolute_couple_kg_m2_s: absolute_contact_couple_kg_m2_s,
        }
        .into());
    }
    let mut actual = Vec::with_capacity(rows.len());
    let mut before = Vec::with_capacity(rows.len());
    let mut maximum_velocity_residual_mps = 0.0_f64;
    let (mut work, mut work_scale) = (Dd::ZERO, Dd::ZERO);
    for (i, &row) in rows.iter().enumerate() {
        checkpoint(cx).map_err(F::from)?;
        let after = velocity(&velocity_mps, row, frames[i])?;
        for a in 0..3 {
            let residual = finite(after[a] - contact.relative_velocities_after[i][a])?.abs();
            maximum_velocity_residual_mps = maximum_velocity_residual_mps.max(residual);
            if residual > controls.velocity_tolerance_mps {
                return Err(F::Response {
                    face: row.face_id,
                    vertex: row.vertex,
                    component: a,
                    residual_mps: residual,
                    tolerance_mps: controls.velocity_tolerance_mps,
                }
                .into());
            }
            let term = Dd::from_f64(contact.impulses[i][a])
                * (Dd::from_f64(0.5 * contacts[i].relative_velocity[a])
                    + Dd::from_f64(0.5 * after[a]));
            work = work + term;
            work_scale = work_scale + Dd::from_f64(rounded(term)?.abs());
        }
        before.push(contacts[i].relative_velocity);
        actual.push(after);
    }
    let order: Vec<_> = (0..rows.len()).collect();
    let measured_residuals = actual_endpoint_residuals(
        &contact,
        &before,
        &actual,
        &order,
        |i| contacts[i].friction,
        controls.contact,
        cx,
    )
    .map_err(F::NativeAudit)?;
    let relative_work_j = rounded(work)?;
    let allowed_work =
        finite(controls.contact.energy_tolerance + 64.0 * f64::EPSILON * rounded(work_scale)?)?;
    if relative_work_j > allowed_work {
        return Err(F::NativeAudit(SolidFrictionalImpactError::Work {
            defect_j: relative_work_j,
            allowance_j: allowed_work,
        })
        .into());
    }
    let contact_work_error_j = rounded(work - Dd::from_f64(contact.relative_work))?;
    let native_work = Dd::from_f64(impulse.kinetic_change_j) - Dd::from_f64(impulse.support_work_j);
    let relative_work_error_j = rounded(native_work - Dd::from_f64(contact.relative_work))?;
    let energy_balance_error_j = rounded(native_work + Dd::from_f64(contact.dissipated_energy))?;
    for (account, defect) in [
        (
            "frictional impact retained kinetic storage",
            impulse.kinetic_before_j - body.kinetic_energy_j,
        ),
        (
            "frictional impact actual/contact work",
            contact_work_error_j,
        ),
        (
            "frictional impact native/contact work",
            relative_work_error_j,
        ),
        (
            "frictional impact native/contact dissipation",
            energy_balance_error_j,
        ),
    ] {
        energy_gate(account, defect, controls.energy_tolerance_j).map_err(F::from)?;
    }
    checkpoint(cx).map_err(F::from)?;
    Ok(ReleasedFaceFrictionResponse {
        velocity_mps,
        contact,
        impulse,
        rows,
        frames,
        measured_residuals,
        maximum_velocity_residual_mps,
        relative_work_j,
        contact_work_error_j,
        relative_work_error_j,
        energy_balance_error_j,
        absolute_contact_couple_kg_m2_s,
    })
}

fn velocity(
    velocity: &[[f64; 3]],
    row: BrittleFaceImpactRow,
    frame: BrittleFaceFrictionFrame,
) -> Result<[f64; 3], BrittleFaceFrictionError> {
    let relative =
        core::array::from_fn(|axis| velocity[row.nodes[1]][axis] - velocity[row.nodes[0]][axis]);
    let mut result = [0.0; 3];
    for (i, axis) in frame.axes().into_iter().enumerate() {
        result[i] = product(axis, relative)?;
    }
    Ok(result)
}
fn product(a: [f64; 3], b: [f64; 3]) -> Result<f64, BrittleFaceFrictionError> {
    rounded(
        Dd::from_f64(a[0]) * Dd::from_f64(b[0])
            + Dd::from_f64(a[1]) * Dd::from_f64(b[1])
            + Dd::from_f64(a[2]) * Dd::from_f64(b[2]),
    )
}
fn unit(vector: [f64; 3]) -> Result<[f64; 3], BrittleFaceFrictionError> {
    if vector.iter().any(|x| !x.is_finite()) {
        return Err(BrittleFaceFrictionError::Invalid(
            "finite contact frame required",
        ));
    }
    let scale = vector.iter().fold(0.0_f64, |scale, x| scale.max(x.abs()));
    if scale == 0.0 {
        return Err(BrittleFaceFrictionError::Invalid(
            "nonzero normal and nonparallel tangent direction required",
        ));
    }
    let scaled = vector.map(|x| x / scale);
    let length = finite(norm(scaled))?;
    Ok(scaled.map(|x| x / length))
}
fn rounded(value: Dd) -> Result<f64, BrittleFaceFrictionError> {
    finite(value.hi)?;
    finite(value.lo)?;
    finite(value.to_f64())
}
fn finite(value: f64) -> Result<f64, BrittleFaceFrictionError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(BrittleFaceFrictionError::Invalid(
            "nonfinite contact geometry, velocity, impulse or work",
        ))
    }
}
