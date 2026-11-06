//! Sustained frictionless forces on the original nonlinear brittle midpoint.
//!
//! The native contact owner solves normal force/configuration complementarity
//! through complete original material and inertia responses. Every force trial
//! starts from the same accepted history. Only the final, independently checked
//! nonpenetrating candidate advances time. Normal forces are interval means,
//! never impulses divided by time or artificial penalty springs.

/// Loaded-contact prefixes localizing arrival of other released faces.
pub mod approach;

use fs_contact::sustained::NaturalMapScales;
use fs_contact::sustained::solve::{MAX_CONTACTS, SolveControls, SolveError, SolveReport};
use fs_exec::Cx;

use super::BrittleFaceImpactRow;
use crate::brittle::solid::dynamics::midpoint::{
    BrittleTetMidpointControls, BrittleTetMidpointEnergy,
};
use crate::brittle::solid::dynamics::{BrittleTetDynamicError, BrittleTetMotion};
use crate::brittle::solid::path::Sum;
use crate::brittle::solid::{BrittleTetError, checkpoint, energy_gate};
use crate::plastic3::IncrementSettings;

/// Native nonlinear solve budgets and independent physical admission gates.
#[derive(Clone, Copy, Debug)]
pub struct BrittleFaceSustainedControls {
    pub solve: SolveControls,
    pub work: BrittleTetMidpointControls,
    /// Numerical conditioning only; no compliance or stored energy is added.
    pub normal_scale_n_per_m: f64,
    pub activation_tolerance_m: f64,
    pub tangential_tolerance_m: f64,
    pub velocity_tolerance_mps: f64,
    pub force_tolerance_n: f64,
    pub complementarity_tolerance_j: f64,
    pub contact_work_tolerance_j: f64,
    pub force_moment_tolerance_n_m: f64,
}

/// One accepted material/motion interval and its original coupled-force solve.
#[derive(Debug)]
pub struct BrittleFaceSustainedStep<'a> {
    pub motion: BrittleTetMotion<'a>,
    /// External dead-load and native support work only. Internal normal-force
    /// work must independently fit its numerical allowance; it is not a loss.
    pub energy: BrittleTetMidpointEnergy,
    /// Unmodified midpoint receipt including the applied contact nodal loads.
    pub native_energy: BrittleTetMidpointEnergy,
    pub solve: SolveReport,
    /// Participating pairs in retained face/vertex order. The standalone solve
    /// includes all non-tied pairs; constrained arrival prefixes include only
    /// the pairs that were initially resting, active and tangentially matched.
    pub rows: Vec<BrittleFaceImpactRow>,
    pub normal_contact_work_j: f64,
    pub absolute_contact_work_j: f64,
    pub contact_work_balance_error_j: f64,
    /// Sum of each pair's larger initial/final absolute force couple. This
    /// also bounds its affine displacement trace; opposing couples cannot hide.
    pub absolute_contact_couple_n_m: f64,
    /// Sum over all successfully completed original mechanical responses.
    pub numerical_tangent_evaluations: usize,
}

#[derive(Debug)]
pub enum BrittleFaceSustainedError {
    Invalid(&'static str),
    Motion(BrittleTetDynamicError),
    Solve(Box<SolveError<BrittleTetDynamicError>>),
    ImpactRequired {
        face: u64,
        vertex: usize,
        normal_speed_mps: f64,
    },
    Geometry {
        face: u64,
        vertex: usize,
        reason: &'static str,
    },
    Mode {
        face: u64,
        vertex: usize,
        reason: &'static str,
    },
    Balance {
        account: &'static str,
        defect: f64,
        tolerance: f64,
    },
}
impl core::fmt::Display for BrittleFaceSustainedError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "released brittle face sustained contact refused: {self:?}"
        )
    }
}
impl std::error::Error for BrittleFaceSustainedError {}
impl From<BrittleTetError> for BrittleFaceSustainedError {
    fn from(error: BrittleTetError) -> Self {
        Self::Motion(error.into())
    }
}

impl<'a> BrittleTetMotion<'a> {
    /// Resolve finite-duration frictionless dwell or unloaded separation.
    ///
    /// All non-tied released pairs participate, including currently open ones.
    /// Loaded pairs must remain active, normally at rest and tangentially
    /// matched at both endpoints. The actual native displacement trace is
    /// affine, so these endpoint conditions bound its opening and contact work.
    /// A closing impact, open-to-loaded arrival, finite sliding or unresolved
    /// normal work refuses; this is not a general continuous collision solver.
    /// No position, gap, material history or velocity is projected or reset.
    ///
    /// Original applied loads are interval means and prescribed values are
    /// endpoint displacements. Native reactions include contact loads on held
    /// coordinates. Failed solves publish no partial motion or consumed step.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_released_face_sustained(
        &self,
        dt_s: f64,
        prescribed_m: &[(usize, f64)],
        mean_nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        controls: BrittleFaceSustainedControls,
        cx: &Cx<'_>,
    ) -> Result<BrittleFaceSustainedStep<'a>, BrittleFaceSustainedError> {
        self.advance_released_face_sustained_from(
            dt_s,
            prescribed_m,
            mean_nodal_forces_n,
            None,
            settings,
            controls,
            cx,
        )
    }

    /// Previous interval forces are numerical seeds only; history is unchanged.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn advance_released_face_sustained_from(
        &self,
        dt_s: f64,
        prescribed_m: &[(usize, f64)],
        mean_nodal_forces_n: &[(usize, f64)],
        initial_normal_forces_n: Option<&[f64]>,
        settings: IncrementSettings,
        controls: BrittleFaceSustainedControls,
        cx: &Cx<'_>,
    ) -> Result<BrittleFaceSustainedStep<'a>, BrittleFaceSustainedError> {
        self.sustained_candidate(
            dt_s,
            prescribed_m,
            mean_nodal_forces_n,
            initial_normal_forces_n,
            SustainedScope::Complete,
            settings,
            controls,
            cx,
        )
    }

    /// Private constrained response. InitialResting defers opening admission
    /// only for unforced rows, which the arrival locator must check before
    /// exposing a motion. All participating-row and physical work gates remain.
    #[allow(clippy::too_many_arguments)]
    fn sustained_candidate(
        &self,
        dt_s: f64,
        prescribed_m: &[(usize, f64)],
        mean_nodal_forces_n: &[(usize, f64)],
        initial_normal_forces_n: Option<&[f64]>,
        scope: SustainedScope,
        settings: IncrementSettings,
        controls: BrittleFaceSustainedControls,
        cx: &Cx<'_>,
    ) -> Result<BrittleFaceSustainedStep<'a>, BrittleFaceSustainedError> {
        use BrittleFaceSustainedError as E;
        checkpoint(cx)?;
        validate_controls(controls)?;
        let old_u = self.solid.solution.displacement_m.as_flattened();
        if mean_nodal_forces_n.len() > old_u.len() || prescribed_m.len() > old_u.len() {
            return Err(E::Invalid("load entries exceed original DOF count"));
        }
        let mut external = vec![0.0; old_u.len()];
        for &(dof, value) in mean_nodal_forces_n {
            checkpoint(cx)?;
            if dof >= external.len() || !value.is_finite() {
                return Err(E::Invalid(
                    "finite external loads on existing DOFs required",
                ));
            }
            external[dof] += value;
            if !external[dof].is_finite() {
                return Err(E::Invalid("unrepresentable accumulated external load"));
            }
        }
        let mut rows = Vec::new();
        let mut normals = Vec::new();
        let mut initial = Vec::new();
        for face in self.solid.released_faces.iter() {
            checkpoint(cx)?;
            for vertex in 0..3 {
                let nodes = [face.paired_nodes[0][vertex], face.paired_nodes[1][vertex]];
                if nodes[0] == nodes[1] {
                    continue;
                }
                let row = BrittleFaceImpactRow {
                    face_id: face.id,
                    vertex,
                    nodes,
                };
                let geometry = pair(self, row, face.normal)?;
                if geometry.gap < 0.0 {
                    return Err(geometry_error(
                        row,
                        "negative initial gap cannot be repaired",
                    ));
                }
                if geometry.gap <= controls.activation_tolerance_m
                    && geometry.speed < -controls.velocity_tolerance_mps
                {
                    return Err(E::ImpactRequired {
                        face: face.id,
                        vertex,
                        normal_speed_mps: geometry.speed,
                    });
                }
                if scope == SustainedScope::InitialResting && !resting(geometry, controls) {
                    continue;
                }
                if rows.len() >= MAX_CONTACTS || rows.len() >= controls.solve.maximum_contacts {
                    return Err(E::Invalid(
                        "all participating released-face rows must fit the sustained solve budget",
                    ));
                }
                rows.push(row);
                normals.push(face.normal);
                initial.push(geometry);
            }
        }
        if rows.is_empty() {
            return Err(E::Invalid("no non-tied released-face contact rows"));
        }
        if initial_normal_forces_n.is_some_and(|forces| {
            forces.len() != rows.len() || forces.iter().any(|f| !f.is_finite() || *f < 0.0)
        }) {
            return Err(E::Invalid(
                "finite nonnegative force seeds in complete row order required",
            ));
        }
        let seeds: Vec<_> = (0..rows.len())
            .map(|i| [initial_normal_forces_n.map_or(0.0, |f| f[i]), 0.0, 0.0])
            .collect();
        let friction = vec![0.0; rows.len()];
        let scales = vec![
            NaturalMapScales {
                normal_n_per_m: controls.normal_scale_n_per_m,
                tangent_n_per_m: controls.normal_scale_n_per_m,
            };
            rows.len()
        ];
        let mut total: Vec<_> = external.iter().copied().enumerate().collect();
        let mut contact_load = vec![0.0; old_u.len()];
        let mut retained = None;
        let mut final_forces = Vec::new();
        let mut evaluations = 0usize;
        let mut tangent_evaluations = 0usize;
        // For mu=0 the two tangent equations are simply Ft=0. Their path
        // coordinates are dummy algebraic zeros: only normal forces are sent
        // to the body, and the returned tangent forces must be exactly zero.
        let solve = fs_contact::sustained::solve::solve(
            &friction,
            &scales,
            &seeds,
            controls.solve,
            &mut |forces, stop| {
                if stop() {
                    return Err(BrittleTetError::Cancelled.into());
                }
                retained = None;
                contact_load.fill(0.0);
                for ((row, normal), force) in rows.iter().zip(&normals).zip(forces) {
                    checkpoint(cx)?;
                    for axis in 0..3 {
                        let value = normal[axis] * force[0];
                        contact_load[3 * row.nodes[0] + axis] -= value;
                        contact_load[3 * row.nodes[1] + axis] += value;
                    }
                }
                for (i, ((_, value), &base)) in total.iter_mut().zip(&external).enumerate() {
                    if i % 256 == 0 {
                        checkpoint(cx)?;
                    }
                    *value = base + contact_load[i];
                    if !value.is_finite() {
                        return Err(BrittleTetDynamicError::Invalid(
                            "unrepresentable sustained total load",
                        ));
                    }
                }
                evaluations += 1;
                let step = self.midpoint_candidate(
                    dt_s,
                    prescribed_m,
                    &total,
                    settings,
                    controls.work,
                    cx,
                )?;
                tangent_evaluations = tangent_evaluations
                    .checked_add(step.numerical_tangent_evaluations)
                    .ok_or(BrittleTetDynamicError::Invalid(
                        "completed tangent count overflow",
                    ))?;
                let mut paths = Vec::with_capacity(rows.len());
                for (row, normal) in rows.iter().zip(&normals) {
                    let gap = opening(&step.motion, *row, *normal);
                    if !gap.is_finite() {
                        return Err(BrittleTetDynamicError::Invalid(
                            "unrepresentable paired response gap",
                        ));
                    }
                    paths.push([gap, 0.0, 0.0]);
                }
                final_forces.clear();
                final_forces.extend_from_slice(forces);
                retained = Some(step);
                Ok(paths)
            },
            &mut || cx.checkpoint().is_err(),
        )
        .map_err(|error| E::Solve(Box::new(error)))?;
        if evaluations != solve.response_evaluations || final_forces != solve.forces_n {
            return Err(E::Invalid("missing original final force replay"));
        }
        let step = retained.ok_or(E::Invalid("missing final native midpoint candidate"))?;
        let new_u = step.motion.solid.solution.displacement_m.as_flattened();
        let mut contact_work = Sum::default();
        let mut absolute_work = Sum::default();
        let mut absolute_couple = Sum::default();
        for (i, ((row, normal), before)) in rows.iter().zip(&normals).zip(&initial).enumerate() {
            checkpoint(cx)?;
            let after = pair(&step.motion, *row, *normal)?;
            let [force, tangent_1, tangent_2] = solve.forces_n[i];
            if force < 0.0
                || tangent_1 != 0.0
                || tangent_2 != 0.0
                || solve.paths_m[i] != [after.gap, 0.0, 0.0]
            {
                return Err(E::Invalid(
                    "final force/path identity or frictionless normal cone refused",
                ));
            }
            if after.gap < 0.0 {
                return Err(geometry_error(
                    *row,
                    "negative final gap cannot be repaired",
                ));
            }
            if force > controls.force_tolerance_n
                && (before.gap > controls.activation_tolerance_m
                    || after.gap > controls.activation_tolerance_m
                    || before.speed.abs() > controls.velocity_tolerance_mps
                    || after.speed.abs() > controls.velocity_tolerance_mps
                    || before.tangent > controls.tangential_tolerance_m
                    || after.tangent > controls.tangential_tolerance_m)
            {
                return Err(E::Mode {
                    face: row.face_id,
                    vertex: row.vertex,
                    reason: "loaded contact must remain matched, active and normally at rest",
                });
            }
            energy_gate(
                "sustained initial normal complementarity",
                force * before.gap,
                controls.complementarity_tolerance_j,
            )?;
            energy_gate(
                "sustained final normal complementarity",
                force * after.gap,
                controls.complementarity_tolerance_j,
            )?;
            let work = force * (after.gap - before.gap);
            contact_work.add(work);
            absolute_work.add(work.abs());
            let physical_force = normal.map(|value| force * value);
            absolute_couple.add(
                super::norm(super::cross(before.separation, physical_force))
                    .max(super::norm(super::cross(after.separation, physical_force))),
            );
        }
        if scope == SustainedScope::Complete {
            self.solid
                .require_released_faces(&step.motion.solid.solution.displacement_m, cx)?;
        }
        let mut external_work = Sum::default();
        for &(dof, force) in mean_nodal_forces_n {
            checkpoint(cx)?;
            external_work.add(force * (new_u[dof] - old_u[dof]));
        }
        let mut nodal_contact_work = Sum::default();
        let mut resultant = [Sum::default(); 3];
        for (dof, &force) in contact_load.iter().enumerate() {
            if dof % 256 == 0 {
                checkpoint(cx)?;
            }
            nodal_contact_work.add(force * (new_u[dof] - old_u[dof]));
            resultant[dof % 3].add(force);
        }
        let normal_contact_work_j = contact_work.value();
        let absolute_contact_work_j = absolute_work.value();
        let absolute_contact_couple_n_m = absolute_couple.value();
        let native_energy = step.energy;
        let contact_work_balance_error_j =
            (native_energy.dead_load_mean_work_j - external_work.value()) - normal_contact_work_j;
        let external = external_work.value() + native_energy.support_mean_work_j;
        let energy = BrittleTetMidpointEnergy {
            dead_load_mean_work_j: external_work.value(),
            equilibrium_work_residual_j: (external - native_energy.mean_internal_work_j)
                - native_energy.inertia.inertial_work_j,
            total_balance_error_j: ((external
                - (native_energy.stored_after_j - native_energy.stored_before_j))
                - native_energy.fracture_dissipation_j)
                - native_energy.inertia.kinetic_change_j,
            ..native_energy
        };
        for (account, defect, tolerance) in [
            (
                "sustained normal contact work",
                normal_contact_work_j,
                controls.contact_work_tolerance_j,
            ),
            (
                "sustained absolute normal contact work",
                absolute_contact_work_j,
                controls.contact_work_tolerance_j,
            ),
            (
                "sustained nodal/paired normal work",
                nodal_contact_work.value() - normal_contact_work_j,
                controls.contact_work_tolerance_j,
            ),
            (
                "sustained external/contact applied work",
                contact_work_balance_error_j,
                controls.contact_work_tolerance_j,
            ),
            (
                "sustained external-only equilibrium work",
                energy.equilibrium_work_residual_j,
                controls.work.equilibrium_work_tolerance_j,
            ),
            (
                "sustained external-only physical work",
                energy.total_balance_error_j,
                controls.work.total_energy_tolerance_j,
            ),
        ] {
            energy_gate(account, defect, tolerance)?;
        }
        balance(
            "internal contact resultant force",
            super::norm(resultant.map(Sum::value)),
            controls.force_tolerance_n,
        )?;
        balance(
            "absolute internal contact force couple",
            absolute_contact_couple_n_m,
            controls.force_moment_tolerance_n_m,
        )?;
        checkpoint(cx)?;
        Ok(BrittleFaceSustainedStep {
            motion: step.motion,
            energy,
            native_energy,
            solve,
            rows,
            normal_contact_work_j,
            absolute_contact_work_j,
            contact_work_balance_error_j,
            absolute_contact_couple_n_m,
            numerical_tangent_evaluations: tangent_evaluations,
        })
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum SustainedScope {
    Complete,
    InitialResting,
}

fn validate_controls(
    controls: BrittleFaceSustainedControls,
) -> Result<(), BrittleFaceSustainedError> {
    if !controls.normal_scale_n_per_m.is_finite()
        || controls.normal_scale_n_per_m <= 0.0
        || [
            controls.activation_tolerance_m,
            controls.tangential_tolerance_m,
            controls.velocity_tolerance_mps,
            controls.force_tolerance_n,
            controls.complementarity_tolerance_j,
            controls.contact_work_tolerance_j,
            controls.force_moment_tolerance_n_m,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v < 0.0)
    {
        return Err(BrittleFaceSustainedError::Invalid(
            "finite nonnegative physical gates and positive normal scale required",
        ));
    }
    Ok(())
}

fn resting(geometry: Pair, controls: BrittleFaceSustainedControls) -> bool {
    geometry.gap <= controls.activation_tolerance_m
        && geometry.speed.abs() <= controls.velocity_tolerance_mps
        && geometry.tangent <= controls.tangential_tolerance_m
}

#[derive(Clone, Copy)]
struct Pair {
    gap: f64,
    speed: f64,
    tangent: f64,
    separation: [f64; 3],
}
fn opening(motion: &BrittleTetMotion<'_>, row: BrittleFaceImpactRow, normal: [f64; 3]) -> f64 {
    super::dot(normal, separation(motion, row))
}
fn separation(motion: &BrittleTetMotion<'_>, row: BrittleFaceImpactRow) -> [f64; 3] {
    core::array::from_fn(|axis| {
        motion.solid.solution.displacement_m[row.nodes[1]][axis]
            - motion.solid.solution.displacement_m[row.nodes[0]][axis]
    })
}
fn pair(
    motion: &BrittleTetMotion<'_>,
    row: BrittleFaceImpactRow,
    normal: [f64; 3],
) -> Result<Pair, BrittleFaceSustainedError> {
    let separation = separation(motion, row);
    let gap = super::dot(normal, separation);
    let speed = super::dot(
        normal,
        core::array::from_fn(|axis| {
            motion.velocity_mps[row.nodes[1]][axis] - motion.velocity_mps[row.nodes[0]][axis]
        }),
    );
    let tangent = super::norm(core::array::from_fn(|axis| {
        separation[axis] - gap * normal[axis]
    }));
    if !gap.is_finite() || !speed.is_finite() || !tangent.is_finite() {
        return Err(geometry_error(
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
fn geometry_error(row: BrittleFaceImpactRow, reason: &'static str) -> BrittleFaceSustainedError {
    BrittleFaceSustainedError::Geometry {
        face: row.face_id,
        vertex: row.vertex,
        reason,
    }
}
fn balance(
    account: &'static str,
    defect: f64,
    tolerance: f64,
) -> Result<(), BrittleFaceSustainedError> {
    if !defect.is_finite() || defect.abs() > tolerance {
        Err(BrittleFaceSustainedError::Balance {
            account,
            defect,
            tolerance,
        })
    } else {
        Ok(())
    }
}
