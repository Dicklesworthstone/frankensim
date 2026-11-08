//! Two-sided normal forces between retained, completely failed ductile faces.
//!
//! Every contact-force probe advances the SAME accepted native finite material,
//! original consistent mass, velocity and physical duration. Only the final
//! independently admitted endpoint is returned. The contact law adds no spring,
//! viscosity, density rescaling or replacement constitutive state.
//!
//! This bounded path resolves frictionless dwell and unloaded separation of
//! matching released triangles with fixed CURRENT release normals. Loaded
//! pairs must remain aligned, matched and normally at rest. Closing impacts,
//! finite sliding and rotation outside the release-frame tolerance refuse.
//! All paired gaps and the native tetrahedral orientation chord are checked;
//! unrelated fragment surfaces require a separate collision-detection owner.

use fs_contact::sustained::NaturalMapScales;
use fs_contact::sustained::solve::{self, SolveReport};
use fs_exec::Cx;
use fs_material::ductile::evolution::DifferentiableDuctileInitiationLaw;
use fs_math::dd::Dd;

use super::super::super::plane_contact::{
    PlaneContactControls, PlaneContactError, TetPlaneEndpoint, TetPlaneGeometry,
};
use super::super::dynamics::FiniteDuctileMotionControls;
use super::super::{DuctileStepConditions, FiniteDuctileTetError, IncrementSettings};
use super::{
    FiniteDuctileFractureAdvance, FiniteDuctileFractureError, FiniteDuctileFractureMotion,
};

/// Native complementarity, endpoint work and geometric admission budgets.
#[derive(Clone, Copy, Debug)]
pub struct FiniteDuctileFaceContactControls {
    /// Original force solver, signed/absolute normal-work gates and tet sweep.
    pub contact: PlaneContactControls,
    /// Natural-map conditioning only; this supplies no physical compliance.
    pub normal_scale_n_per_m: f64,
    pub activation_tolerance_m: f64,
    pub tangential_tolerance_m: f64,
    pub velocity_tolerance_mps: f64,
    pub force_tolerance_n: f64,
    pub complementarity_tolerance_j: f64,
    /// Sum of each pair's largest endpoint force couple; opposing couples
    /// cannot cancel one another in this angular-momentum gate.
    pub force_moment_tolerance_n_m: f64,
}

impl Default for FiniteDuctileFaceContactControls {
    fn default() -> Self {
        Self {
            contact: PlaneContactControls::default(),
            normal_scale_n_per_m: 1e6,
            activation_tolerance_m: 1e-10,
            tangential_tolerance_m: 1e-10,
            velocity_tolerance_mps: 1e-9,
            force_tolerance_n: 1e-8,
            complementarity_tolerance_j: 1e-8,
            force_moment_tolerance_n_m: 1e-8,
        }
    }
}

/// Stable released-face identity and its fixed current-frame normal.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteDuctileFaceContactRow {
    pub face_id: u64,
    pub vertex: usize,
    /// A and B nodes; a positive normal force acts -n on A and +n on B.
    pub nodes: [usize; 2],
    pub normal: [f64; 3],
}

/// Contact work is kept separate from physical plastic/damage loss and heat.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteDuctileFaceContactReport {
    pub rows: Vec<FiniteDuctileFaceContactRow>,
    pub solve: SolveReport,
    pub gap_before_m: Vec<f64>,
    /// Work of equal/opposite internal pair forces on the actual native step.
    pub normal_contact_work_j: f64,
    /// Negative contact work, retained as a signed endpoint quadrature defect.
    /// It is neither fracture energy nor an impact restitution law.
    pub normal_work_defect_j: f64,
    pub absolute_normal_work_defect_j: f64,
    pub absolute_contact_couple_n_m: f64,
    /// Original caller loads only; the native step's dead-load work includes
    /// internal contact loads and must not be added a second time.
    pub dead_load_work_j: f64,
    pub load_decomposition_error_j: f64,
    /// External/support work minus native physical/numerical accounts and the
    /// separate signed contact defect.
    pub work_balance_error_j: f64,
    pub minimum_gap_m: f64,
    pub path_windows: usize,
}

/// Native accepted material/momentum endpoint with the released faces retained.
#[derive(Debug)]
pub struct FiniteDuctileFaceContactStep<'g, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + core::fmt::Debug + PartialEq,
{
    pub step: FiniteDuctileFractureAdvance<'g, 'm, L>,
    pub contact: FiniteDuctileFaceContactReport,
}

#[derive(Debug)]
pub enum FiniteDuctileFaceContactError {
    Invalid(&'static str),
    Native(FiniteDuctileTetError),
    Fracture(FiniteDuctileFractureError),
    Contact(solve::SolveError<FiniteDuctileTetError>),
    Path(PlaneContactError<FiniteDuctileTetError>),
    Face {
        face: u64,
        vertex: usize,
        reason: &'static str,
    },
    ImpactRequired {
        face: u64,
        vertex: usize,
        normal_speed_mps: f64,
    },
    Balance {
        quantity: &'static str,
        value: f64,
        maximum_abs: f64,
    },
}
impl core::fmt::Display for FiniteDuctileFaceContactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite ductile released-face contact refused: {self:?}")
    }
}
impl std::error::Error for FiniteDuctileFaceContactError {}
impl From<FiniteDuctileTetError> for FiniteDuctileFaceContactError {
    fn from(error: FiniteDuctileTetError) -> Self {
        Self::Native(error)
    }
}
impl From<FiniteDuctileFractureError> for FiniteDuctileFaceContactError {
    fn from(error: FiniteDuctileFractureError) -> Self {
        Self::Fracture(error)
    }
}

impl<'g, 'm, L> FiniteDuctileFractureMotion<'g, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + core::fmt::Debug + PartialEq,
{
    /// Continue normal dwell or unloaded separation using the original finite
    /// backward-Euler solve. All non-tied released vertex pairs participate.
    ///
    /// Optional force seeds follow released-face/vertex order with tied pairs
    /// omitted; seeds are numerical guesses, not retained mechanical history.
    /// Every refused force trial leaves this accepted wrapper unchanged.
    /// Incoming closing impacts need a separate impulse/event owner. Neither
    /// penetration nor sliding is repaired by moving or resetting a vertex.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn advance_with_face_contact(
        &self,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        conditions: DuctileStepConditions<'_>,
        initial_normal_forces_n: Option<&[f64]>,
        settings: IncrementSettings,
        material_controls: FiniteDuctileMotionControls,
        controls: FiniteDuctileFaceContactControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteDuctileFaceContactStep<'g, 'm, L>, FiniteDuctileFaceContactError> {
        use FiniteDuctileFaceContactError as E;
        poll(cx)?;
        validate(controls)?;
        let solid = self.motion.solid();
        let nodes = solid.reference_nodes_m();
        let old_u = solid.displacement_m();
        self.admit_endpoint(old_u, cx)?;
        let geometry = TetPlaneGeometry {
            nodes_m: nodes,
            tetrahedra: solid.reference_tetrahedra(),
            elements: &solid.state.elements,
            displacement_m: old_u,
        };
        let mut rows = Vec::new();
        let mut before = Vec::new();
        for face in self.released_faces.iter() {
            poll(cx)?;
            for vertex in 0..3 {
                let nodes = [face.paired_nodes[0][vertex], face.paired_nodes[1][vertex]];
                if nodes[0] == nodes[1] {
                    continue;
                }
                if rows.len() >= solve::MAX_CONTACTS
                    || rows.len() >= controls.contact.solve.maximum_contacts
                {
                    return Err(E::Invalid(
                        "all released pairs must fit the native contact budget",
                    ));
                }
                let row = FiniteDuctileFaceContactRow {
                    face_id: face.id,
                    vertex,
                    nodes,
                    normal: face.normal,
                };
                let pair = pair_view(
                    solid.reference_nodes_m(),
                    old_u,
                    self.motion.velocity_mps(),
                    row,
                )?;
                if pair.gap < 0.0 {
                    return Err(face_error(row, "negative initial gap cannot be repaired"));
                }
                if pair.gap <= controls.activation_tolerance_m
                    && pair.speed < -controls.velocity_tolerance_mps
                {
                    return Err(E::ImpactRequired {
                        face: row.face_id,
                        vertex,
                        normal_speed_mps: pair.speed,
                    });
                }
                if !(0..3).any(|axis| {
                    row.normal[axis] != 0.0
                        && row.nodes.iter().any(|node| {
                            !prescribed_m.iter().any(|&(dof, _)| dof == 3 * node + axis)
                        })
                }) {
                    return Err(face_error(
                        row,
                        "a contact normal needs a free nodal coordinate",
                    ));
                }
                rows.push(row);
                before.push(pair);
            }
        }
        if rows.is_empty()
            || initial_normal_forces_n.is_some_and(|forces| {
                forces.len() != rows.len() || forces.iter().any(|f| !f.is_finite() || *f < 0.0)
            })
        {
            return Err(E::Invalid(
                "nonempty released rows and complete nonnegative seeds required",
            ));
        }
        let seeds: Vec<_> = (0..rows.len())
            .map(|i| [initial_normal_forces_n.map_or(0.0, |f| f[i]), 0.0, 0.0])
            .collect();
        let scales = vec![
            NaturalMapScales {
                normal_n_per_m: controls.normal_scale_n_per_m,
                tangent_n_per_m: controls.normal_scale_n_per_m,
            };
            rows.len()
        ];
        let friction = vec![0.0; rows.len()];
        let mut retained = None;
        let report = {
            let mut response = |forces: &[[f64; 3]],
                                stop: &mut dyn FnMut() -> bool|
             -> Result<Vec<[f64; 3]>, FiniteDuctileTetError> {
                drop(retained.take());
                if stop() {
                    return Err(super::super::SmallStrainTetError::Cancelled.into());
                }
                let mut applied = Vec::with_capacity(nodal_forces_n.len() + 6 * rows.len());
                applied.extend_from_slice(nodal_forces_n);
                for (row, force) in rows.iter().zip(forces) {
                    for axis in 0..3 {
                        let value = row.normal[axis] * force[0];
                        if !value.is_finite() {
                            return Err(super::super::invalid(
                                "finite paired normal force required",
                            )
                            .into());
                        }
                        if value != 0.0 {
                            applied.push((3 * row.nodes[0] + axis, -value));
                            applied.push((3 * row.nodes[1] + axis, value));
                        }
                    }
                }
                let step = self.motion.advance(
                    prescribed_m,
                    &applied,
                    body_force_n_m3,
                    conditions,
                    settings,
                    material_controls,
                    cx,
                )?;
                let mut paths = Vec::with_capacity(rows.len());
                for row in &rows {
                    let gap = opening(nodes, step.motion.solid().displacement_m(), *row);
                    if !gap.is_finite() {
                        return Err(
                            super::super::invalid("finite paired response gap required").into()
                        );
                    }
                    paths.push([gap, 0.0, 0.0]);
                }
                retained = Some((forces.to_vec(), paths.clone(), step));
                Ok(paths)
            };
            solve::solve(
                &friction,
                &scales,
                &seeds,
                controls.contact.solve,
                &mut response,
                &mut || cx.checkpoint().is_err(),
            )
            .map_err(E::Contact)?
        };
        let (forces, paths, step) =
            retained.ok_or(E::Invalid("missing final native force replay"))?;
        if forces != report.forces_n || paths != report.paths_m {
            return Err(E::Invalid(
                "final native response must match the contact receipt",
            ));
        }
        let new_u = step.motion.solid().displacement_m();
        self.admit_endpoint(new_u, cx)?;
        let path_windows = geometry
            .path::<FiniteDuctileTetError>(new_u, controls.contact, cx)
            .map_err(E::Path)?;
        let mut contact_work = Dd::ZERO;
        let mut absolute_work = Dd::ZERO;
        let mut absolute_couple = Dd::ZERO;
        let mut minimum_gap = f64::INFINITY;
        for (i, (row, old)) in rows.iter().zip(&before).enumerate() {
            poll(cx)?;
            let after = pair_view(nodes, new_u, step.motion.velocity_mps(), *row)?;
            let [force, tangent_1, tangent_2] = report.forces_n[i];
            if force < 0.0
                || tangent_1 != 0.0
                || tangent_2 != 0.0
                || report.paths_m[i] != [after.gap, 0.0, 0.0]
            {
                return Err(face_error(
                    *row,
                    "frictionless cone or final path identity refused",
                ));
            }
            if after.gap < 0.0 {
                return Err(face_error(*row, "negative final gap cannot be repaired"));
            }
            if force > controls.force_tolerance_n
                && (old.gap > controls.activation_tolerance_m
                    || after.gap > controls.activation_tolerance_m
                    || old.speed.abs() > controls.velocity_tolerance_mps
                    || after.speed.abs() > controls.velocity_tolerance_mps
                    || old.tangent > controls.tangential_tolerance_m
                    || after.tangent > controls.tangential_tolerance_m)
            {
                return Err(face_error(
                    *row,
                    "loaded pairs must remain matched, active and normally at rest",
                ));
            }
            gate(
                "initial normal complementarity",
                force * old.gap,
                controls.complementarity_tolerance_j,
            )?;
            gate(
                "final normal complementarity",
                force * after.gap,
                controls.complementarity_tolerance_j,
            )?;
            // Physical nodal-force work, evaluated before summing across pairs.
            let mut work = Dd::ZERO;
            for axis in 0..3 {
                let [a, b] = row.nodes;
                let relative_increment = (Dd::from_f64(new_u[b][axis])
                    - Dd::from_f64(old_u[b][axis]))
                    - (Dd::from_f64(new_u[a][axis]) - Dd::from_f64(old_u[a][axis]));
                work = work + Dd::from_f64(force * row.normal[axis]) * relative_increment;
            }
            contact_work = contact_work + work;
            absolute_work = absolute_work + Dd::from_f64(work.to_f64().abs());
            let physical_force = row.normal.map(|n| force * n);
            let couple = norm(cross(old.separation, physical_force))
                .max(norm(cross(after.separation, physical_force)));
            absolute_couple = absolute_couple + Dd::from_f64(couple);
            minimum_gap = minimum_gap.min(old.gap).min(after.gap);
        }
        let absolute_normal_work_defect_j = absolute_work.to_f64();
        gate(
            "absolute normal endpoint work defect",
            absolute_normal_work_defect_j,
            controls.contact.maximum_normal_work_defect_j,
        )?;
        gate(
            "absolute paired force couple",
            absolute_couple.to_f64(),
            controls.force_moment_tolerance_n_m,
        )?;
        let dead = geometry
            .dead_work::<FiniteDuctileTetError>(nodal_forces_n, body_force_n_m3, new_u, cx)
            .map_err(E::Path)?;
        let decomposition = (Dd::from_f64(step.dead_load_work_j()) - dead - contact_work).to_f64();
        gate(
            "native dead/contact load decomposition",
            decomposition,
            controls.contact.work_balance_tolerance_j,
        )?;
        let normal_defect = Dd::ZERO - contact_work;
        let balance = step
            .balance(dead + Dd::from_f64(step.support_work_j()), normal_defect)
            .to_f64();
        gate(
            "native external/material/contact work identity",
            balance,
            controls.contact.work_balance_tolerance_j,
        )?;
        poll(cx)?;
        Ok(FiniteDuctileFaceContactStep {
            step: self.retain_step(step),
            contact: FiniteDuctileFaceContactReport {
                rows,
                solve: report,
                gap_before_m: before.iter().map(|p| p.gap).collect(),
                normal_contact_work_j: contact_work.to_f64(),
                normal_work_defect_j: normal_defect.to_f64(),
                absolute_normal_work_defect_j,
                absolute_contact_couple_n_m: absolute_couple.to_f64(),
                dead_load_work_j: dead.to_f64(),
                load_decomposition_error_j: decomposition,
                work_balance_error_j: balance,
                minimum_gap_m: minimum_gap,
                path_windows,
            },
        })
    }
}

#[derive(Clone, Copy)]
struct Pair {
    gap: f64,
    speed: f64,
    tangent: f64,
    separation: [f64; 3],
}
fn separation(nodes: &[[f64; 3]], u: &[[f64; 3]], row: FiniteDuctileFaceContactRow) -> [f64; 3] {
    let [a, b] = row.nodes;
    core::array::from_fn(|axis| {
        (Dd::from_f64(nodes[b][axis]) - Dd::from_f64(nodes[a][axis]) + Dd::from_f64(u[b][axis])
            - Dd::from_f64(u[a][axis]))
        .to_f64()
    })
}
fn opening(nodes: &[[f64; 3]], u: &[[f64; 3]], row: FiniteDuctileFaceContactRow) -> f64 {
    dot(row.normal, separation(nodes, u, row))
}
fn pair_view(
    nodes: &[[f64; 3]],
    u: &[[f64; 3]],
    velocity: &[[f64; 3]],
    row: FiniteDuctileFaceContactRow,
) -> Result<Pair, FiniteDuctileFaceContactError> {
    let separation = separation(nodes, u, row);
    let gap = dot(row.normal, separation);
    let speed = dot(
        row.normal,
        core::array::from_fn(|axis| velocity[row.nodes[1]][axis] - velocity[row.nodes[0]][axis]),
    );
    let tangent = norm(core::array::from_fn(|axis| {
        separation[axis] - gap * row.normal[axis]
    }));
    if !gap.is_finite() || !speed.is_finite() || !tangent.is_finite() {
        return Err(face_error(
            row,
            "finite paired geometry and normal speed required",
        ));
    }
    Ok(Pair {
        gap,
        speed,
        tangent,
        separation,
    })
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    let mut value = Dd::ZERO;
    for axis in 0..3 {
        value = value + Dd::from_f64(a[axis]) * Dd::from_f64(b[axis]);
    }
    value.to_f64()
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn norm(v: [f64; 3]) -> f64 {
    v.into_iter().fold(0.0_f64, |n, x| n.hypot(x))
}
fn face_error(
    row: FiniteDuctileFaceContactRow,
    reason: &'static str,
) -> FiniteDuctileFaceContactError {
    FiniteDuctileFaceContactError::Face {
        face: row.face_id,
        vertex: row.vertex,
        reason,
    }
}
fn poll(cx: &Cx<'_>) -> Result<(), FiniteDuctileFaceContactError> {
    cx.checkpoint().map_err(|_| {
        FiniteDuctileFaceContactError::Native(super::super::SmallStrainTetError::Cancelled.into())
    })
}
fn gate(
    quantity: &'static str,
    value: f64,
    maximum_abs: f64,
) -> Result<(), FiniteDuctileFaceContactError> {
    if !value.is_finite() || value.abs() > maximum_abs {
        Err(FiniteDuctileFaceContactError::Balance {
            quantity,
            value,
            maximum_abs,
        })
    } else {
        Ok(())
    }
}
fn validate(c: FiniteDuctileFaceContactControls) -> Result<(), FiniteDuctileFaceContactError> {
    if !c.normal_scale_n_per_m.is_finite()
        || c.normal_scale_n_per_m <= 0.0
        || [
            c.activation_tolerance_m,
            c.tangential_tolerance_m,
            c.velocity_tolerance_mps,
            c.force_tolerance_n,
            c.complementarity_tolerance_j,
            c.force_moment_tolerance_n_m,
            c.contact.maximum_normal_work_defect_j,
            c.contact.work_balance_tolerance_j,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v < 0.0)
        || !(1..=1 << 24).contains(&c.contact.maximum_path_windows)
        || !(1..=1 << 20).contains(&c.contact.sweep.maximum_windows)
        || !c.contact.sweep.fraction_tolerance.is_finite()
        || c.contact.sweep.fraction_tolerance <= 0.0
        || c.contact.sweep.fraction_tolerance >= 1.0
    {
        return Err(FiniteDuctileFaceContactError::Invalid(
            "positive normal scale, finite nonnegative gates and bounded path controls required",
        ));
    }
    Ok(())
}
