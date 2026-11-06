//! Simultaneous frictionless impacts on retained completed-fracture faces.
//!
//! Every active paired vertex contributes its signed physical normal row.
//! Native consistent-mass responses assemble the complete contact mobility,
//! including all same-body cross terms. The original contact owner solves one
//! coupled normal event; one summed native impulse updates the actual body.
//! Position, fracture history, source state, reference matter and time do not
//! change. This fixed-configuration event does not locate collisions, repair
//! penetration, supply persistent contact forces or infer impact heat.

/// Bounded native midpoint prefixes locating approaching released-face events.
pub mod approach;

/// Native nonlinear midpoint motion with sustained unilateral normal forces.
pub mod sustained;

/// Complete native physical intervals through simultaneous released-face impacts.
pub mod time;

use std::sync::Arc;

use fs_contact::multicontact::{
    MultiImpactConfig, MultiImpactError, MultiImpactResult, NormalContactImpact,
    solve_normal_impacts,
};
use fs_exec::Cx;

use crate::brittle::solid::dynamics::{BrittleTetDynamicError, BrittleTetMotion};
use crate::brittle::solid::fracture::ReleasedBrittleFace;
use crate::brittle::solid::{BrittleTetError, BrittleTetSolid, checkpoint, energy_gate};
use crate::inertia::P1ReferenceMass;
pub use crate::inertia::{SolidImpulseControls, SolidImpulseReport};
use crate::plastic3::SmallStrainTetError;

/// Original contact and mass-solve budgets, plus physical assembly gates.
#[derive(Clone, Copy, Debug)]
pub struct BrittleFaceImpactControls {
    pub impact: MultiImpactConfig,
    pub impulse: SolidImpulseControls,
    pub activation_tolerance_m: f64,
    /// Active paired witnesses must lie on the same reference normal line.
    pub tangential_tolerance_m: f64,
    /// Actual endpoint-speed, restitution and complementarity allowance.
    pub velocity_tolerance_mps: f64,
    pub energy_tolerance_j: f64,
    /// Hard cap in 1..=64 on ALL active non-tied paired vertices.
    /// New contact arrays are quadratic in this cap; native mass scratch uses
    /// `impulse.maximum_workspace_components` and existing body dimensions.
    pub maximum_contact_rows: usize,
}

/// Stable material-face and vertex identity for one native normal-contact row.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BrittleFaceImpactRow {
    pub face_id: u64,
    pub vertex: usize,
    /// A-side and B-side current nodal indices, respectively.
    pub nodes: [usize; 2],
}

/// Complete velocity successor with unchanged position, history and clock.
#[derive(Debug)]
pub struct BrittleFaceImpactStep<'a> {
    pub motion: BrittleTetMotion<'a>,
    pub contact: MultiImpactResult,
    pub impulse: SolidImpulseReport,
    pub rows: Vec<BrittleFaceImpactRow>,
    pub maximum_velocity_residual_mps: f64,
    /// Native kinetic change minus support work minus contact relative work.
    pub relative_work_error_j: f64,
    /// Native kinetic change minus support work plus contact dissipation.
    /// Admitted numerical positive work remains visible in this account.
    pub energy_balance_error_j: f64,
    /// Sum of individual |(x_B-x_A) cross impulse|. Opposing artificial
    /// contact couples cannot cancel the native angular-momentum gate.
    pub absolute_contact_couple_kg_m2_s: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum BrittleFaceImpactError {
    Invalid(&'static str),
    Motion(BrittleTetDynamicError),
    Contact(MultiImpactError),
    Geometry {
        face: u64,
        vertex: usize,
        reason: &'static str,
    },
    Response {
        face: u64,
        vertex: usize,
        residual_mps: f64,
        tolerance_mps: f64,
    },
    InternalImpulse {
        linear_norm_kg_m_s: f64,
        angular_norm_kg_m2_s: f64,
        absolute_couple_kg_m2_s: f64,
    },
}
impl core::fmt::Display for BrittleFaceImpactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "released brittle face impact refused: {self:?}")
    }
}
impl std::error::Error for BrittleFaceImpactError {}
impl From<BrittleTetError> for BrittleFaceImpactError {
    fn from(error: BrittleTetError) -> Self {
        Self::Motion(error.into())
    }
}
impl From<SmallStrainTetError> for BrittleFaceImpactError {
    fn from(error: SmallStrainTetError) -> Self {
        Self::Motion(error.into())
    }
}

impl<'a> BrittleTetMotion<'a> {
    /// Apply one simultaneous frictionless event to every active released pair.
    ///
    /// Restitution follows `solid().released_faces()` order. Shared crack-front
    /// vertices have identically zero relative motion and are omitted. Every
    /// other pair at or below the activation gap is included, even if initially
    /// separating. Penetration always refuses, including inside the solid's
    /// opening roundoff allowance. Active tangentially separated witnesses
    /// refuse: fixed paired vertices do not define finite-sliding contact.
    /// An empty active set returns a refusal without publishing a candidate.
    ///
    /// Held velocity components retain their current values through the native
    /// mass owner; their signed external reaction work is explicit. The full
    /// coupled contact law retains its rank, passivity and determinacy limits.
    /// Redundant impulses remain a representative, not a unique local pressure
    /// or heat allocation. Actual endpoint speeds, internal impulse moments and
    /// whole-body work must agree before the complete successor is returned.
    /// The physical clock, step allowance and cumulative fracture loss do not
    /// change. A rejected event leaves this accepted motion available for retry.
    pub fn impact_released_faces(
        &self,
        restitution: &[f64],
        held_velocity_dofs: &[usize],
        origin_m: [f64; 3],
        controls: BrittleFaceImpactControls,
        cx: &Cx<'_>,
    ) -> Result<BrittleFaceImpactStep<'a>, BrittleFaceImpactError> {
        let ReleasedFaceImpactResponse {
            velocity_mps,
            contact,
            impulse,
            rows,
            maximum_velocity_residual_mps,
            relative_work_error_j,
            energy_balance_error_j,
            absolute_contact_couple_kg_m2_s,
        } = released_face_impact(
            ReleasedFaceImpactBody {
                reference_nodes_m: self.solid.geometry.nodes_m(),
                displacement_m: &self.solid.solution.displacement_m,
                released_faces: &self.solid.released_faces,
                mass: &self.mass,
                velocity_mps: &self.velocity_mps,
                kinetic_energy_j: self.kinetic_energy_j,
            },
            restitution,
            held_velocity_dofs,
            origin_m,
            controls,
            cx,
        )?;
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
        Ok(BrittleFaceImpactStep {
            motion: BrittleTetMotion {
                solid,
                velocity_mps,
                reference_density_kg_m3: Arc::clone(&self.reference_density_kg_m3),
                mass: Arc::clone(&self.mass),
                kinetic_energy_j: impulse.kinetic_after_j,
                time_s: self.time_s,
                maximum_steps: self.maximum_steps,
                accepted_steps: self.accepted_steps,
            },
            contact,
            impulse,
            rows,
            maximum_velocity_residual_mps,
            relative_work_error_j,
            energy_balance_error_j,
            absolute_contact_couple_kg_m2_s,
        })
    }
}

/// Read-only native geometry and inertia; constitutive owners stay with callers.
/// Both thermal and isothermal bodies enter the identical physical impulse law.
pub(crate) struct ReleasedFaceImpactBody<'a> {
    pub(crate) reference_nodes_m: &'a [[f64; 3]],
    pub(crate) displacement_m: &'a [[f64; 3]],
    pub(crate) released_faces: &'a [ReleasedBrittleFace],
    pub(crate) mass: &'a P1ReferenceMass,
    pub(crate) velocity_mps: &'a [[f64; 3]],
    pub(crate) kinetic_energy_j: f64,
}

/// Native event result before its caller publishes the unchanged material owner.
pub(crate) struct ReleasedFaceImpactResponse {
    pub(crate) velocity_mps: Vec<[f64; 3]>,
    pub(crate) contact: MultiImpactResult,
    pub(crate) impulse: SolidImpulseReport,
    pub(crate) rows: Vec<BrittleFaceImpactRow>,
    pub(crate) maximum_velocity_residual_mps: f64,
    pub(crate) relative_work_error_j: f64,
    pub(crate) energy_balance_error_j: f64,
    pub(crate) absolute_contact_couple_kg_m2_s: f64,
}

pub(crate) fn released_face_impact(
    body: ReleasedFaceImpactBody<'_>,
    restitution: &[f64],
    held_velocity_dofs: &[usize],
    origin_m: [f64; 3],
    controls: BrittleFaceImpactControls,
    cx: &Cx<'_>,
) -> Result<ReleasedFaceImpactResponse, BrittleFaceImpactError> {
    use BrittleFaceImpactError as E;
    checkpoint(cx)?;
    if restitution.len() != body.released_faces.len()
        || restitution
            .iter()
            .any(|e| !e.is_finite() || !(0.0..=1.0).contains(e))
        || !(1..=64).contains(&controls.maximum_contact_rows)
        || [
            controls.activation_tolerance_m,
            controls.tangential_tolerance_m,
            controls.velocity_tolerance_mps,
            controls.energy_tolerance_j,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v < 0.0)
        || origin_m.iter().any(|v| !v.is_finite())
    {
        return Err(E::Invalid(
            "complete restitution, bounded rows and finite nonnegative gates required",
        ));
    }
    let mut rows = Vec::with_capacity(controls.maximum_contact_rows);
    let mut normals = Vec::with_capacity(controls.maximum_contact_rows);
    let mut separations = Vec::with_capacity(controls.maximum_contact_rows);
    let mut contacts = Vec::with_capacity(controls.maximum_contact_rows);
    for (face, &restitution) in body.released_faces.iter().zip(restitution) {
        checkpoint(cx)?;
        for vertex in 0..3 {
            let [a, b] = [face.paired_nodes[0][vertex], face.paired_nodes[1][vertex]];
            if a == b {
                continue;
            }
            // Corresponding reference coordinates are identical by the
            // original topology admission. Subtract displacements directly;
            // adding a large common spatial origin first can erase opening.
            let separation = core::array::from_fn(|axis| {
                body.displacement_m[b][axis] - body.displacement_m[a][axis]
            });
            let gap = finite(dot(face.normal, separation))?;
            if gap < 0.0 {
                return Err(E::Geometry {
                    face: face.id,
                    vertex,
                    reason: "negative paired opening; penetration cannot be repaired",
                });
            }
            if gap > controls.activation_tolerance_m {
                continue;
            }
            let tangent = core::array::from_fn(|axis| separation[axis] - gap * face.normal[axis]);
            if norm(tangent) > controls.tangential_tolerance_m {
                return Err(E::Geometry {
                    face: face.id,
                    vertex,
                    reason: "active paired witnesses are tangentially separated",
                });
            }
            if rows.len() == controls.maximum_contact_rows
                || rows.len() >= controls.impact.solver.max_variables
            {
                return Err(E::Invalid(
                    "all active released-face rows must fit the contact budget",
                ));
            }
            let row = BrittleFaceImpactRow {
                face_id: face.id,
                vertex,
                nodes: [a, b],
            };
            contacts.push(NormalContactImpact {
                gap,
                activation_tolerance: controls.activation_tolerance_m,
                relative_normal_speed: speed(body.velocity_mps, row, face.normal)?,
                restitution,
            });
            rows.push(row);
            normals.push(face.normal);
            separations.push(separation);
        }
    }
    if rows.is_empty() {
        return Err(E::Invalid("no active released-face contact rows"));
    }
    let count = rows.len();
    let mut mobility = vec![0.0; count * count];
    for (column, (&row, &normal)) in rows.iter().zip(&normals).enumerate() {
        checkpoint(cx)?;
        let (response, _) = body
            .mass
            .free_velocity_response(
                held_velocity_dofs,
                &loads(row, normal, 1.0),
                controls.impulse,
                cx,
            )
            .map_err(SmallStrainTetError::from)?;
        for i in 0..count {
            mobility[i * count + column] = speed(&response, rows[i], normals[i])?;
        }
    }
    let contact = solve_normal_impacts(&contacts, &mobility, controls.impact, || {
        cx.checkpoint().is_err()
    })
    .map_err(E::Contact)?;
    let mut impulses = Vec::with_capacity(6 * count);
    let mut absolute_contact_couple_kg_m2_s = 0.0;
    for (i, &amount) in contact.normal_impulses.iter().enumerate() {
        checkpoint(cx)?;
        impulses.extend_from_slice(&loads(rows[i], normals[i], amount));
        let impulse = normals[i].map(|n| n * amount);
        absolute_contact_couple_kg_m2_s =
            finite(absolute_contact_couple_kg_m2_s + norm(cross(separations[i], impulse)))?;
    }
    let mut positions = Vec::with_capacity(body.velocity_mps.len());
    for (node, (reference, displacement)) in body
        .reference_nodes_m
        .iter()
        .zip(body.displacement_m)
        .enumerate()
    {
        if node % 256 == 0 {
            checkpoint(cx)?;
        }
        let position = core::array::from_fn(|axis| reference[axis] + displacement[axis]);
        if position.iter().any(|x| !x.is_finite()) {
            return Err(E::Invalid("unrepresentable physical position"));
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
        .map_err(SmallStrainTetError::from)?;
    let linear = norm(impulse.applied.linear_kg_m_s);
    let angular = norm(impulse.applied.angular_kg_m2_s);
    if linear > controls.impulse.linear_tolerance_kg_m_s
        || angular > controls.impulse.angular_tolerance_kg_m2_s
        || absolute_contact_couple_kg_m2_s > controls.impulse.angular_tolerance_kg_m2_s
    {
        return Err(E::InternalImpulse {
            linear_norm_kg_m_s: linear,
            angular_norm_kg_m2_s: angular,
            absolute_couple_kg_m2_s: absolute_contact_couple_kg_m2_s,
        });
    }
    let mut maximum_velocity_residual_mps = 0.0_f64;
    for (i, &row) in rows.iter().enumerate() {
        checkpoint(cx)?;
        let actual = speed(&velocity_mps, row, normals[i])?;
        let stationarity = finite(actual - contact.restitution_targets[i])?;
        let kkt = if contact.normal_impulses[i] > 0.0 {
            stationarity.abs()
        } else {
            (-stationarity).max(0.0)
        };
        let residual = finite(actual - contact.relative_normal_speeds_after[i])?
            .abs()
            .max(kkt);
        maximum_velocity_residual_mps = maximum_velocity_residual_mps.max(residual);
        if residual > controls.velocity_tolerance_mps {
            return Err(E::Response {
                face: row.face_id,
                vertex: row.vertex,
                residual_mps: residual,
                tolerance_mps: controls.velocity_tolerance_mps,
            });
        }
    }
    let relative_work_error_j =
        (impulse.kinetic_change_j - impulse.support_work_j) - contact.relative_work;
    let energy_balance_error_j =
        (impulse.kinetic_change_j - impulse.support_work_j) + contact.dissipated_energy;
    for (account, defect) in [
        (
            "impact retained kinetic storage",
            impulse.kinetic_before_j - body.kinetic_energy_j,
        ),
        ("impact native/contact relative work", relative_work_error_j),
        ("impact native/contact dissipation", energy_balance_error_j),
    ] {
        energy_gate(account, defect, controls.energy_tolerance_j)?;
    }
    checkpoint(cx)?;
    Ok(ReleasedFaceImpactResponse {
        velocity_mps,
        contact,
        impulse,
        rows,
        maximum_velocity_residual_mps,
        relative_work_error_j,
        energy_balance_error_j,
        absolute_contact_couple_kg_m2_s,
    })
}

fn speed(
    velocity: &[[f64; 3]],
    row: BrittleFaceImpactRow,
    normal: [f64; 3],
) -> Result<f64, BrittleFaceImpactError> {
    finite(dot(
        normal,
        core::array::from_fn(|axis| velocity[row.nodes[1]][axis] - velocity[row.nodes[0]][axis]),
    ))
}
fn loads(row: BrittleFaceImpactRow, normal: [f64; 3], amount: f64) -> [(usize, f64); 6] {
    core::array::from_fn(|i| {
        (
            3 * row.nodes[i / 3] + i % 3,
            (if i < 3 { -amount } else { amount }) * normal[i % 3],
        )
    })
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0].mul_add(b[0], a[1].mul_add(b[1], a[2] * b[2]))
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn norm(a: [f64; 3]) -> f64 {
    a[0].hypot(a[1]).hypot(a[2])
}
fn finite(value: f64) -> Result<f64, BrittleFaceImpactError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(BrittleFaceImpactError::Invalid(
            "nonfinite contact geometry or arithmetic",
        ))
    }
}
