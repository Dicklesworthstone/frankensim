//! Sustained Coulomb contact on fixed physical shell point/plane witnesses.
//!
//! The original gauge-free `fs-plate` point rows observe displacement and their
//! exact transpose applies forces and physical moments. These are the original
//! P1 interface traces, not enriched DKT interior reconstruction. The existing sustained
//! contact solver queries native plastic-shell equilibria or backward-Euler
//! momentum endpoints, always from the same accepted material histories.
//!
//! Geometry and contact frames remain fixed. The admitted path linearly
//! interpolates the two point/plane endpoints; nonnegative endpoint gaps exclude
//! penetration for that affine point path only. This is not a contact search,
//! finite-rotation, footprint, self-intersection or global CCD certificate.
//! Explicit rigid offsets add neither inertia nor compliance. Their physical
//! ownership and the supplied finite patch areas belong to the caller.

use super::drilling::{PlasticShellReferenceDrilling, ShellDrillingEnergy};
use super::dynamics::{PlasticShellReferenceInertia, ShellInertialEnergy};
use super::{
    PlasticShellSections, SmallStrainShellProblem, SmallStrainShellSolution, checkpoint, invalid,
};
use crate::plastic3::SmallStrainTetError;
use fs_contact::sustained::NaturalMapScales;
use fs_contact::sustained::solve::{self, SolveControls, SolveError, SolveReport};
use fs_exec::Cx;
use fs_material::SmallStrainLaw;
use fs_material::plastic::section::LayeredPlaneStressState;
use fs_plate::shell::ShellMesh;
use fs_plate::shell::point::triangle_offset_point_rows;
use std::fmt;

/// Original physical point and fixed oriented obstacle plane. The normal points
/// from the obstacle toward the shell; a positive normal force opens the gap.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShellContactWitness {
    pub id: u64,
    pub triangle: usize,
    pub barycentric: [f64; 3],
    /// Global reference arm from the barycentric midsurface point [m].
    pub offset_m: [f64; 3],
    pub reference_obstacle_point_m: [f64; 3],
    /// Right-handed orthonormal rows: normal, first tangent, second tangent.
    pub axes_world: [[f64; 3]; 3],
    /// Explicit physical patch area [m²], used only for contact pressure.
    pub area_m2: f64,
}

/// Inputs for one interval, in the prepared witnesses' original order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShellContactStep {
    /// Absolute obstacle translation relative to the reference witness [m].
    pub obstacle_start_displacement_m: [f64; 3],
    pub obstacle_end_displacement_m: [f64; 3],
    /// Finite nonnegative coefficient, held over this physical interval.
    pub friction_coefficient: f64,
    /// Numerical initial guess [normal, tangent 1, tangent 2], in N.
    pub initial_force_n: [f64; 3],
    pub scales: NaturalMapScales,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShellContactControls {
    pub solve: SolveControls,
    /// Independent absolute contact work/complementarity gate [J].
    pub work_tolerance_j: f64,
    /// Maximum aggregate fn*g0 [J], separate from the work-closure tolerance.
    /// Zero by default: contact may start open only when no closing force is
    /// required. Callers admitting endpoint approach must set a finite budget.
    pub maximum_normal_numerical_loss_j: f64,
}

impl Default for ShellContactControls {
    fn default() -> Self {
        Self {
            solve: SolveControls::default(),
            work_tolerance_j: 1e-9,
            maximum_normal_numerical_loss_j: 0.0,
        }
    }
}

/// Original source errors remain inspectable; no state is published on refusal.
#[derive(Debug)]
pub enum ShellContactError {
    Invalid(&'static str),
    Solid(SmallStrainTetError),
    Contact(SolveError<SmallStrainTetError>),
    Work {
        id: u64,
        residual_j: f64,
        tolerance_j: f64,
    },
    /// Endpoint contact would spend too much unresolved approach work.
    NormalApproachBudget {
        loss_j: f64,
        maximum_j: f64,
    },
}

impl From<SmallStrainTetError> for ShellContactError {
    fn from(error: SmallStrainTetError) -> Self {
        Self::Solid(error)
    }
}

impl fmt::Display for ShellContactError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid(reason) => write!(f, "invalid plastic-shell contact: {reason}"),
            Self::Solid(error) => write!(f, "original plastic shell: {error}"),
            Self::Contact(error) => write!(f, "original sustained contact: {error}"),
            Self::Work {
                id,
                residual_j,
                tolerance_j,
            } => write!(
                f,
                "shell contact {id} work residual {residual_j} J exceeds {tolerance_j} J"
            ),
            Self::NormalApproachBudget { loss_j, maximum_j } => write!(
                f,
                "shell contact approach loss {loss_j} J exceeds {maximum_j} J"
            ),
        }
    }
}

impl std::error::Error for ShellContactError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Solid(error) => Some(error),
            Self::Contact(error) => Some(error),
            _ => None,
        }
    }
}

#[derive(Clone, Debug)]
struct Point {
    nodes: [usize; 3],
    rows: [[[f64; 6]; 3]; 3],
    reference_gap_m: f64,
}

/// Immutable reference geometry, retaining caller order and IDs. Section laws
/// may change between thermal trials; the mesh and physical layer bounds may not.
#[derive(Clone, Debug)]
pub struct PlasticShellContacts<'m> {
    mesh: &'m ShellMesh,
    witnesses: Vec<ShellContactWitness>,
    points: Vec<Point>,
    layer_bounds: Vec<Vec<[f64; 2]>>,
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.iter().zip(b).fold(0.0, |sum, (x, y)| x.mul_add(y, sum))
}

fn section_count<L: SmallStrainLaw>(problem: &SmallStrainShellProblem<'_, L>) -> usize {
    match &problem.sections {
        PlasticShellSections::Uniform(_) => problem.mesh.tris.len(),
        PlasticShellSections::PerElement(sections) => sections.len(),
        PlasticShellSections::PerElementBorrowed(sections) => sections.len(),
    }
}

impl<'m> PlasticShellContacts<'m> {
    /// Prepare without invoking any material law. The returned lifetime borrows
    /// only the original mesh, independently of temporary thermal section laws.
    pub fn prepare<L: SmallStrainLaw>(
        mesh: &'m ShellMesh,
        problem: &SmallStrainShellProblem<'_, L>,
        witnesses: &[ShellContactWitness],
        maximum_contacts: usize,
        cx: &Cx<'_>,
    ) -> Result<Self, ShellContactError> {
        checkpoint(cx)?;
        if maximum_contacts == 0
            || maximum_contacts > solve::MAX_CONTACTS
            || witnesses.is_empty()
            || witnesses.len() > maximum_contacts
            || (!std::ptr::eq(mesh, problem.mesh) && mesh != problem.mesh)
            || mesh.nodes.len() < 3
            || mesh.nodes.len() > problem.budget.maximum_nodes
            || mesh.tris.is_empty()
            || mesh.tris.len() > problem.budget.maximum_elements
            || section_count(problem) != mesh.tris.len()
            || mesh.nodes.iter().flatten().any(|x| !x.is_finite())
        {
            return Err(ShellContactError::Invalid(
                "reference geometry, section count or contact budget",
            ));
        }
        let mut count = 0usize;
        for e in 0..mesh.tris.len() {
            checkpoint(cx)?;
            count = problem
                .sections
                .at(e)
                .point_count()
                .checked_mul(3)
                .and_then(|points| count.checked_add(points))
                .ok_or(ShellContactError::Invalid("material point count overflow"))?;
            if count > problem.budget.maximum_material_points {
                return Err(ShellContactError::Invalid(
                    "reference material point budget",
                ));
            }
        }
        let mut layer_bounds = Vec::with_capacity(mesh.tris.len());
        for e in 0..mesh.tris.len() {
            checkpoint(cx)?;
            layer_bounds.push(
                problem
                    .sections
                    .at(e)
                    .layers()
                    .iter()
                    .map(|layer| [layer.bottom_z_m, layer.top_z_m])
                    .collect(),
            );
        }
        let mut points = Vec::with_capacity(witnesses.len());
        for (i, witness) in witnesses.iter().enumerate() {
            checkpoint(cx)?;
            if witnesses[..i].iter().any(|old| old.id == witness.id)
                || !witness.area_m2.is_finite()
                || witness.area_m2 <= 0.0
                || witness
                    .reference_obstacle_point_m
                    .iter()
                    .any(|x| !x.is_finite())
                || witness.axes_world.iter().flatten().any(|x| !x.is_finite())
            {
                return Err(ShellContactError::Invalid(
                    "unique IDs, finite plane and positive patch area required",
                ));
            }
            let [n, t1, t2] = witness.axes_world;
            let cross = [
                n[1] * t1[2] - n[2] * t1[1],
                n[2] * t1[0] - n[0] * t1[2],
                n[0] * t1[1] - n[1] * t1[0],
            ];
            if (0..3).any(|a| {
                (0..3).any(|b| {
                    (dot3(witness.axes_world[a], witness.axes_world[b])
                        - if a == b { 1.0 } else { 0.0 })
                    .abs()
                        > 1e-12
                })
            }) || (dot3(cross, t2) - 1.0).abs() > 1e-12
            {
                return Err(ShellContactError::Invalid(
                    "contact frame must be right-handed and orthonormal",
                ));
            }
            let geometry = mesh
                .facet(witness.triangle)
                .map_err(|error| ShellContactError::Solid(invalid(error.to_string())))?;
            let nodes = mesh.tris[witness.triangle];
            let mut rows = [[[0.0; 6]; 3]; 3];
            for (axis, row) in rows.iter_mut().enumerate() {
                *row = triangle_offset_point_rows(
                    witness.barycentric,
                    witness.axes_world[axis],
                    witness.offset_m,
                    geometry.frame,
                    geometry.gradient,
                    [geometry.frame[2]; 3],
                )
                .map_err(|error| ShellContactError::Invalid(error.reason()))?;
            }
            let separation = core::array::from_fn(|axis| {
                (0..3)
                    .map(|a| witness.barycentric[a] * mesh.nodes[nodes[a]][axis])
                    .sum::<f64>()
                    + witness.offset_m[axis]
                    - witness.reference_obstacle_point_m[axis]
            });
            let reference_gap_m = dot3(n, separation);
            if !reference_gap_m.is_finite() {
                return Err(ShellContactError::Invalid(
                    "reference point/plane separation overflows",
                ));
            }
            points.push(Point {
                nodes,
                rows,
                reference_gap_m,
            });
        }
        Ok(Self {
            mesh,
            witnesses: witnesses.to_vec(),
            points,
            layer_bounds,
        })
    }

    /// Caller-supplied original order, including exact physical patch areas.
    pub fn witnesses(&self) -> &[ShellContactWitness] {
        &self.witnesses
    }

    fn admit<L: SmallStrainLaw>(
        &self,
        problem: &SmallStrainShellProblem<'_, L>,
        old: Option<&[[f64; 6]]>,
        steps: &[ShellContactStep],
        controls: ShellContactControls,
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, ShellContactError> {
        checkpoint(cx)?;
        if (!std::ptr::eq(self.mesh, problem.mesh) && self.mesh != problem.mesh)
            || section_count(problem) != self.layer_bounds.len()
            || steps.len() != self.points.len()
            || !controls.work_tolerance_j.is_finite()
            || controls.work_tolerance_j <= 0.0
            || !controls.maximum_normal_numerical_loss_j.is_finite()
            || controls.maximum_normal_numerical_loss_j < 0.0
            || old.is_some_and(|q| {
                q.len() != self.mesh.nodes.len() || q.iter().flatten().any(|x| !x.is_finite())
            })
        {
            return Err(ShellContactError::Invalid(
                "reference geometry, old coordinates, interval count or work gate",
            ));
        }
        for (e, bounds) in self.layer_bounds.iter().enumerate() {
            checkpoint(cx)?;
            let layers = problem.sections.at(e).layers();
            if layers.len() != bounds.len()
                || layers
                    .iter()
                    .zip(bounds)
                    .any(|(layer, bound)| [layer.bottom_z_m, layer.top_z_m] != *bound)
            {
                return Err(ShellContactError::Invalid(
                    "reference physical layer bounds changed",
                ));
            }
        }
        let mut gaps = Vec::with_capacity(steps.len());
        for ((point, witness), step) in self.points.iter().zip(&self.witnesses).zip(steps) {
            checkpoint(cx)?;
            if step
                .obstacle_start_displacement_m
                .iter()
                .chain(&step.obstacle_end_displacement_m)
                .any(|x| !x.is_finite())
            {
                return Err(ShellContactError::Invalid(
                    "finite obstacle endpoint translations required",
                ));
            }
            let gap = point.normal_gap(
                old,
                witness.axes_world[0],
                step.obstacle_start_displacement_m,
            );
            if !gap.is_finite() || gap < 0.0 {
                return Err(ShellContactError::Invalid(
                    "initial physical point/plane gap must be nonnegative",
                ));
            }
            gaps.push(gap);
        }
        Ok(gaps)
    }

    fn loads(
        &self,
        base: &[(usize, f64)],
        forces: &[[f64; 3]],
    ) -> Result<Vec<(usize, f64)>, SmallStrainTetError> {
        let capacity = base
            .len()
            .checked_add(self.points.len() * 18)
            .ok_or_else(|| invalid("shell contact load count overflow"))?;
        let mut loads = Vec::with_capacity(capacity);
        loads.extend_from_slice(base);
        for (point, force) in self.points.iter().zip(forces) {
            for a in 0..3 {
                for d in 0..6 {
                    let value = (0..3).fold(0.0, |sum, axis| {
                        point.rows[axis][a][d].mul_add(force[axis], sum)
                    });
                    if !value.is_finite() {
                        return Err(invalid("shell contact transpose force overflows"));
                    }
                    if value != 0.0 {
                        loads.push((6 * point.nodes[a] + d, value));
                    }
                }
            }
        }
        Ok(loads)
    }

    fn paths(
        &self,
        old: Option<&[[f64; 6]]>,
        new: &[[f64; 6]],
        steps: &[ShellContactStep],
    ) -> Vec<[f64; 3]> {
        self.points
            .iter()
            .zip(&self.witnesses)
            .zip(steps)
            .map(|((point, witness), step)| {
                let delta = point.observe(Some(new), old);
                let obstacle = core::array::from_fn(|axis| {
                    step.obstacle_end_displacement_m[axis]
                        - step.obstacle_start_displacement_m[axis]
                });
                let relative: [f64; 3] = core::array::from_fn(|axis| {
                    delta[axis] - dot3(witness.axes_world[axis], obstacle)
                });
                // Use the same absolute endpoint arithmetic as restart
                // admission. Incrementally accumulating g0 + B*dq - dplane
                // can accept a gap that recomputes as negative on continuation.
                let gap = point.normal_gap(
                    Some(new),
                    witness.axes_world[0],
                    step.obstacle_end_displacement_m,
                );
                [gap, relative[1], relative[2]]
            })
            .collect()
    }
}

impl Point {
    fn normal_gap(
        &self,
        coordinates: Option<&[[f64; 6]]>,
        normal: [f64; 3],
        obstacle_displacement_m: [f64; 3],
    ) -> f64 {
        self.reference_gap_m + self.observe_axis(0, coordinates, None)
            - dot3(normal, obstacle_displacement_m)
    }

    /// Physical metres/radians; rotations are never multiplied by solver length.
    fn observe(&self, new: Option<&[[f64; 6]]>, old: Option<&[[f64; 6]]>) -> [f64; 3] {
        core::array::from_fn(|axis| self.observe_axis(axis, new, old))
    }

    fn observe_axis(&self, axis: usize, new: Option<&[[f64; 6]]>, old: Option<&[[f64; 6]]>) -> f64 {
        (0..3).fold(0.0, |sum, a| {
            (0..6).fold(sum, |sum, d| {
                let delta = new.map_or(0.0, |q| q[self.nodes[a]][d])
                    - old.map_or(0.0, |q| q[self.nodes[a]][d]);
                self.rows[axis][a][d].mul_add(delta, sum)
            })
        })
    }
}

/// One solved contact's local work account. Numerical gap closing is never heat.
#[derive(Clone, Debug, PartialEq)]
pub struct ShellContactReport {
    pub id: u64,
    pub force_n: [f64; 3],
    pub initial_gap_m: f64,
    pub final_gap_m: f64,
    pub slip_increment_m: [f64; 2],
    /// Solved contact force divided by the explicit physical patch area [Pa].
    pub normal_pressure_pa: f64,
    /// Work of the held contact force on the physical shell point [J].
    pub shell_work_j: f64,
    /// Work supplied to the shell/contact pair by obstacle translation [J].
    pub obstacle_work_j: f64,
    /// Native Coulomb loss mu*fn*|slip| [J]; its equality to -ft·slip is gated.
    pub friction_loss_j: f64,
    /// fn*g0 [J]: endpoint approach loss, not restitution or material heat.
    pub normal_numerical_loss_j: f64,
    /// Signed normal contact work fn*(g1-g0) [J].
    pub normal_work_j: f64,
    /// fn*g1 [J], independently required to lie within the work gate.
    pub normal_constraint_residual_j: f64,
    /// Wshell-Wobstacle+friction_loss+normal_numerical_loss [J].
    pub work_residual_j: f64,
}

#[derive(Debug)]
pub struct ShellContactMotion {
    pub velocity: Vec<[f64; 6]>,
    pub inertia: ShellInertialEnergy,
}

#[derive(Debug)]
pub struct PlasticShellContactStep<S> {
    pub solution: SmallStrainShellSolution<S>,
    pub drilling: ShellDrillingEnergy,
    pub motion: Option<ShellContactMotion>,
    pub contacts: Vec<ShellContactReport>,
    pub contact_solve: SolveReport,
}

struct Candidate<S> {
    solution: SmallStrainShellSolution<S>,
    drilling: ShellDrillingEnergy,
    motion: Option<ShellContactMotion>,
    forces: Vec<[f64; 3]>,
    paths: Vec<[f64; 3]>,
}

impl<L: SmallStrainLaw> SmallStrainShellProblem<'_, L> {
    fn with_contact_loads<'a>(
        &'a self,
        loads: &'a [(usize, f64)],
    ) -> SmallStrainShellProblem<'a, L> {
        SmallStrainShellProblem {
            mesh: self.mesh,
            sections: match &self.sections {
                PlasticShellSections::Uniform(section) => PlasticShellSections::Uniform(section),
                PlasticShellSections::PerElement(sections) => {
                    PlasticShellSections::PerElement(sections)
                }
                PlasticShellSections::PerElementBorrowed(sections) => {
                    PlasticShellSections::PerElementBorrowed(sections)
                }
            },
            coordinate_length_m: self.coordinate_length_m,
            rotation_limit_rad: self.rotation_limit_rad,
            drilling_alpha: self.drilling_alpha,
            prescribed: self.prescribed,
            nodal_loads: loads,
            load_path: self.load_path,
            budget: self.budget,
            settings: self.settings,
        }
    }

    /// Solve one factor-one static endpoint with held contact forces. Every
    /// Newton/difference trial restarts from the same accepted history; timed
    /// material updates are not accumulated across contact iterations.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_from_with_contacts(
        &self,
        coordinates: &[[f64; 6]],
        states: &[LayeredPlaneStressState<L::State>],
        drilling: &PlasticShellReferenceDrilling<'_>,
        contacts: &PlasticShellContacts<'_>,
        steps: &[ShellContactStep],
        controls: ShellContactControls,
        cx: &Cx<'_>,
    ) -> Result<PlasticShellContactStep<L::State>, ShellContactError> {
        let gaps = contacts.admit(self, Some(coordinates), steps, controls, cx)?;
        if self.load_path != [1.0] {
            return Err(ShellContactError::Invalid(
                "contact requires one factor-one shell endpoint",
            ));
        }
        let elements = self.prepare_with_reference_drilling(Some(drilling), cx)?;
        contacts.solve(Some(coordinates), steps, &gaps, controls, cx, |forces| {
            let loads = contacts.loads(self.nodal_loads, forces)?;
            let problem = self.with_contact_loads(&loads);
            let solution =
                problem.solve_prepared(&elements, Some((coordinates, states)), None, cx)?;
            let drilling = problem.audit_reference_drilling(
                &elements,
                Some(coordinates),
                &solution.coordinates,
                cx,
            )?;
            Ok((solution, drilling, None))
        })
    }

    /// Advance on the original retained mass and drilling gauge. Material laws
    /// must describe the same dt. Contact does not turn backward-Euler loss,
    /// gauge work or the normal approach contribution into frictional heat.
    /// Incoming normal velocity is handled only by this endpoint discretization;
    /// its loss remains in the native BE numerical account. This does not resolve
    /// an impact time, restitution law or intermediate material/contact history.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_backward_euler_with_contacts(
        &self,
        initial: Option<&SmallStrainShellSolution<L::State>>,
        velocity: &[[f64; 6]],
        inertia: &PlasticShellReferenceInertia<'_>,
        drilling: &PlasticShellReferenceDrilling<'_>,
        contacts: &PlasticShellContacts<'_>,
        steps: &[ShellContactStep],
        dt_s: f64,
        controls: ShellContactControls,
        cx: &Cx<'_>,
    ) -> Result<PlasticShellContactStep<L::State>, ShellContactError> {
        let old = initial.map(|state| state.coordinates.as_slice());
        let gaps = contacts.admit(self, old, steps, controls, cx)?;
        let elements = self.prepare_with_reference_drilling(Some(drilling), cx)?;
        contacts.solve(old, steps, &gaps, controls, cx, |forces| {
            let loads = contacts.loads(self.nodal_loads, forces)?;
            let problem = self.with_contact_loads(&loads);
            let solved = problem
                .advance_backward_euler_prepared(&elements, initial, velocity, inertia, dt_s, cx)?;
            let drilling = problem.audit_reference_drilling(
                &elements,
                old,
                &solved.solution.coordinates,
                cx,
            )?;
            Ok((
                solved.solution,
                drilling,
                Some(ShellContactMotion {
                    velocity: solved.velocity,
                    inertia: solved.inertia,
                }),
            ))
        })
    }
}

impl PlasticShellContacts<'_> {
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    fn solve<S>(
        &self,
        old: Option<&[[f64; 6]]>,
        steps: &[ShellContactStep],
        gaps: &[f64],
        controls: ShellContactControls,
        cx: &Cx<'_>,
        mut mechanical: impl FnMut(
            &[[f64; 3]],
        ) -> Result<
            (
                SmallStrainShellSolution<S>,
                ShellDrillingEnergy,
                Option<ShellContactMotion>,
            ),
            SmallStrainTetError,
        >,
    ) -> Result<PlasticShellContactStep<S>, ShellContactError> {
        let coefficients: Vec<_> = steps.iter().map(|step| step.friction_coefficient).collect();
        let scales: Vec<_> = steps.iter().map(|step| step.scales).collect();
        let forces: Vec<_> = steps.iter().map(|step| step.initial_force_n).collect();
        let mut accepted = None;
        let contact_solve = solve::solve(
            &coefficients,
            &scales,
            &forces,
            controls.solve,
            &mut |force, cancelled| {
                if cancelled() {
                    return Err(SmallStrainTetError::Cancelled);
                }
                checkpoint(cx)?;
                let (solution, drilling, motion) = mechanical(force)?;
                let paths = self.paths(old, &solution.coordinates, steps);
                accepted = Some(Candidate {
                    solution,
                    drilling,
                    motion,
                    forces: force.to_vec(),
                    paths: paths.clone(),
                });
                Ok(paths)
            },
            &mut || cx.checkpoint().is_err(),
        )
        .map_err(ShellContactError::Contact)?;
        checkpoint(cx)?;
        let accepted = accepted.ok_or(ShellContactError::Invalid(
            "contact solver omitted its final mechanical replay",
        ))?;
        if accepted.forces != contact_solve.forces_n || accepted.paths != contact_solve.paths_m {
            return Err(ShellContactError::Invalid(
                "contact replay differs from the returned force/path",
            ));
        }
        let mut reports = Vec::with_capacity(steps.len());
        let mut normal_loss_j = 0.0;
        for i in 0..steps.len() {
            checkpoint(cx)?;
            let witness = &self.witnesses[i];
            let step = &steps[i];
            let force = contact_solve.forces_n[i];
            let path = contact_solve.paths_m[i];
            let delta = self.points[i].observe(Some(&accepted.solution.coordinates), old);
            let obstacle = core::array::from_fn(|axis| {
                step.obstacle_end_displacement_m[axis] - step.obstacle_start_displacement_m[axis]
            });
            let obstacle_axes = witness.axes_world.map(|axis| dot3(axis, obstacle));
            let shell_work_j = dot3(force, delta);
            let obstacle_work_j = dot3(force, obstacle_axes);
            let friction_loss_j = step.friction_coefficient * force[0] * path[1].hypot(path[2]);
            let tangent_work = force[1].mul_add(path[1], force[2] * path[2]);
            let normal_numerical_loss_j = force[0] * gaps[i];
            normal_loss_j += normal_numerical_loss_j;
            if !normal_loss_j.is_finite()
                || normal_loss_j > controls.maximum_normal_numerical_loss_j
            {
                return Err(ShellContactError::NormalApproachBudget {
                    loss_j: normal_loss_j,
                    maximum_j: controls.maximum_normal_numerical_loss_j,
                });
            }
            let normal_constraint_residual_j = force[0] * path[0];
            let normal_work_j = force[0] * (path[0] - gaps[i]);
            let normal_pressure_pa = force[0] / witness.area_m2;
            let work_residual_j =
                shell_work_j - obstacle_work_j + friction_loss_j + normal_numerical_loss_j;
            if [
                shell_work_j,
                obstacle_work_j,
                friction_loss_j,
                tangent_work,
                normal_numerical_loss_j,
                normal_constraint_residual_j,
                normal_work_j,
                normal_pressure_pa,
                work_residual_j,
            ]
            .iter()
            .any(|x| !x.is_finite())
                || force[0] < 0.0
                || path[0] < 0.0
            {
                return Err(ShellContactError::Invalid(
                    "unrepresentable or inadmissible contact work account",
                ));
            }
            // Neither a friction defect nor normal complementarity work may
            // disappear by cancelling the other in the aggregate identity.
            for residual_j in [
                tangent_work + friction_loss_j,
                normal_constraint_residual_j,
                work_residual_j,
            ] {
                if residual_j.abs() > controls.work_tolerance_j {
                    return Err(ShellContactError::Work {
                        id: witness.id,
                        residual_j,
                        tolerance_j: controls.work_tolerance_j,
                    });
                }
            }
            reports.push(ShellContactReport {
                id: witness.id,
                force_n: force,
                initial_gap_m: gaps[i],
                final_gap_m: path[0],
                slip_increment_m: [path[1], path[2]],
                normal_pressure_pa,
                shell_work_j,
                obstacle_work_j,
                friction_loss_j,
                normal_numerical_loss_j,
                normal_work_j,
                normal_constraint_residual_j,
                work_residual_j,
            });
        }
        Ok(PlasticShellContactStep {
            solution: accepted.solution,
            drilling: accepted.drilling,
            motion: accepted.motion,
            contacts: reports,
            contact_solve,
        })
    }
}
