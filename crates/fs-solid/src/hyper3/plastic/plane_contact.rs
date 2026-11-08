//! Shared rigid-plane contact for original finite tetrahedral endpoint owners.
//!
//! This module owns geometry, native normal complementarity and contact work.
//! The supplied endpoint retains all constitutive, momentum and physical-time
//! behavior. Every force probe starts from the same accepted specimen; only the
//! original final response matching the contact receipt can be published.
//!
//! All original vertices must stay on the admitted side of every tool. For P1
//! tetrahedra, the declared straight nodal/tool path then has affine gaps; the
//! native orientation sweep separately excludes collapse along that chord.
//! This is bounded frictionless tooling, not self-contact, an impact law,
//! curved tooling, generated crack surfaces or an exact-flow CCD claim.

use super::{Element, HyperTetError};
use fs_contact::sustained::NaturalMapScales;
use fs_contact::sustained::solve::{self, SolveControls, SolveError, SolveReport};
use fs_exec::Cx;
use fs_math::dd::Dd;
use fs_query::TriangleQueryError;
use fs_query::triangle::sweep::{
    LinearSweepControls, LinearSweepVerdict, linear_tetrahedron_orientation,
};

/// Translating rigid half-space. Its unit normal points into admissible space.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidPlaneTool {
    pub normal: [f64; 3],
    pub point_before_m: [f64; 3],
    /// Translation over the same physical interval as the native material.
    pub translation_m: [f64; 3],
}

/// One explicit original vertex/plane contact candidate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TetPlaneContact {
    pub node: usize,
    pub plane: usize,
    /// Nonnegative numerical initial guess; it does not advance history.
    pub initial_force_n: f64,
    /// Natural-map conditioning scale [N/m], not material stiffness.
    pub scale_n_per_m: f64,
}

/// Native force-solve budgets and independent geometric/work admission.
#[derive(Clone, Copy, Debug)]
pub struct PlaneContactControls {
    pub solve: SolveControls,
    /// Sum of absolute normal endpoint-work defects, never physical heat [J].
    pub maximum_normal_work_defect_j: f64,
    pub work_balance_tolerance_j: f64,
    pub sweep: LinearSweepControls,
    pub maximum_path_windows: usize,
}
impl Default for PlaneContactControls {
    fn default() -> Self {
        Self {
            solve: SolveControls::default(),
            maximum_normal_work_defect_j: 1e-6,
            work_balance_tolerance_j: 1e-7,
            sweep: LinearSweepControls::default(),
            maximum_path_windows: 1 << 20,
        }
    }
}

/// Original complementarity receipt, equal/opposite tool forces and work ports.
#[derive(Clone, Debug, PartialEq)]
pub struct PlaneContactReport {
    pub tools: Vec<RigidPlaneTool>,
    pub contacts: Vec<TetPlaneContact>,
    /// Original row order; both auxiliary tangent forces are zero.
    pub solve: SolveReport,
    pub gap_before_m: Vec<f64>,
    /// Held contact force dotted with the original nodal displacement change.
    pub body_work_j: f64,
    /// Held contact force dotted with the prescribed tool translation.
    pub tool_work_j: f64,
    /// Signed tool minus body work; never deposited as material heat.
    pub normal_work_defect_j: f64,
    pub absolute_normal_work_defect_j: f64,
    /// Caller dead nodal/body work, excluding the added contact forces.
    pub dead_load_work_j: f64,
    pub load_decomposition_error_j: f64,
    pub work_balance_error_j: f64,
    pub force_on_tool_n: Vec<[f64; 3]>,
    /// Compensated numerical minimum over every original vertex and plane.
    pub minimum_gap_m: f64,
    pub path_windows: usize,
}

/// Geometry/contact/work refusal preserving the original native error type.
#[derive(Debug)]
pub enum PlaneContactError<E> {
    Invalid(&'static str),
    Cancelled,
    Solid(E),
    Contact(SolveError<E>),
    Geometry(TriangleQueryError),
    Penetration {
        plane: usize,
        node: usize,
        at_end: bool,
        gap_m: f64,
    },
    PathBudget,
    PossibleElementCollapse {
        element: usize,
        fraction_start: f64,
        fraction_end: f64,
    },
    Energy {
        quantity: &'static str,
        value_j: f64,
        maximum_abs_j: f64,
    },
}
impl<E: core::fmt::Debug> core::fmt::Display for PlaneContactError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite plane contact refused: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for PlaneContactError<E> {}

/// Borrowed original geometry; no mesh, mass or constitutive state is rebuilt.
pub(super) struct TetPlaneGeometry<'a> {
    pub nodes_m: &'a [[f64; 3]],
    pub tetrahedra: &'a [[usize; 4]],
    pub elements: &'a [Element],
    pub displacement_m: &'a [[f64; 3]],
}

/// Each native endpoint subtracts its own complete physical/numerical accounts.
/// A nonpotential material need not invent any incremental scalar potential.
pub(super) trait TetPlaneEndpoint {
    fn displacement_m(&self) -> &[[f64; 3]];
    fn dead_load_work_j(&self) -> f64;
    fn support_work_j(&self) -> f64;
    fn balance(&self, external_work: Dd, normal_defect: Dd) -> Dd;
}

fn poll<E>(cx: &Cx<'_>) -> Result<(), PlaneContactError<E>> {
    cx.checkpoint().map_err(|_| PlaneContactError::Cancelled)
}
fn gate<E>(
    quantity: &'static str,
    value_j: f64,
    maximum_abs_j: f64,
) -> Result<(), PlaneContactError<E>> {
    if !value_j.is_finite() || value_j.abs() > maximum_abs_j {
        Err(PlaneContactError::Energy {
            quantity,
            value_j,
            maximum_abs_j,
        })
    } else {
        Ok(())
    }
}
fn native_invalid<E: From<HyperTetError>>(what: &'static str) -> E {
    HyperTetError::InvalidInput { what: what.into() }.into()
}

impl RigidPlaneTool {
    fn validate<E>(self) -> Result<(), PlaneContactError<E>> {
        let norm = self.normal.into_iter().fold(0.0_f64, |v, x| v.hypot(x));
        if self
            .normal
            .into_iter()
            .chain(self.point_before_m)
            .chain(self.translation_m)
            .any(|x| !x.is_finite())
            || (norm - 1.0).abs() > 128.0 * f64::EPSILON
        {
            return Err(PlaneContactError::Invalid(
                "finite plane point/translation and unit normal required",
            ));
        }
        Ok(())
    }
    fn gap(self, reference: [f64; 3], displacement: [f64; 3], at_end: bool) -> f64 {
        let mut result = Dd::ZERO;
        for axis in 0..3 {
            let relative = Dd::from_f64(reference[axis]) - Dd::from_f64(self.point_before_m[axis])
                + Dd::from_f64(displacement[axis])
                - Dd::from_f64(if at_end {
                    self.translation_m[axis]
                } else {
                    0.0
                });
            result = result + Dd::from_f64(self.normal[axis]) * relative;
        }
        result.to_f64()
    }
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(super) fn solve_tet_planes<S: TetPlaneEndpoint, E: From<HyperTetError>>(
    geometry: TetPlaneGeometry<'_>,
    prescribed: &[(usize, f64)],
    nodal: &[(usize, f64)],
    body: [f64; 3],
    tools: &[RigidPlaneTool],
    contacts: &[TetPlaneContact],
    controls: PlaneContactControls,
    mut endpoint: impl FnMut(&[(usize, f64)]) -> Result<S, E>,
    cx: &Cx<'_>,
) -> Result<(S, PlaneContactReport), PlaneContactError<E>> {
    poll(cx)?;
    geometry.admit(tools, contacts, prescribed, controls)?;
    let minimum_before = geometry.plane_gaps(tools, geometry.displacement_m, false, cx)?;
    let gaps_before = contacts
        .iter()
        .map(|row| {
            tools[row.plane].gap(
                geometry.nodes_m[row.node],
                geometry.displacement_m[row.node],
                false,
            )
        })
        .collect();
    let scales: Vec<_> = contacts
        .iter()
        .map(|row| NaturalMapScales {
            normal_n_per_m: row.scale_n_per_m,
            tangent_n_per_m: row.scale_n_per_m,
        })
        .collect();
    let initial: Vec<_> = contacts
        .iter()
        .map(|row| [row.initial_force_n, 0.0, 0.0])
        .collect();
    let friction = vec![0.0; contacts.len()];
    let mut retained = None;
    let report = {
        let mut response = |forces: &[[f64; 3]],
                            cancelled: &mut dyn FnMut() -> bool|
         -> Result<Vec<[f64; 3]>, E> {
            drop(retained.take());
            if cancelled() {
                return Err(HyperTetError::Cancelled.into());
            }
            let mut applied = Vec::with_capacity(nodal.len() + 3 * contacts.len());
            applied.extend_from_slice(nodal);
            for (row, force) in contacts.iter().zip(forces) {
                for axis in 0..3 {
                    let value = tools[row.plane].normal[axis] * force[0];
                    if !value.is_finite() {
                        return Err(native_invalid("finite contact nodal force required"));
                    }
                    if value != 0.0 {
                        applied.push((3 * row.node + axis, value));
                    }
                }
            }
            let step = endpoint(&applied)?;
            let mut paths = Vec::with_capacity(contacts.len());
            for row in contacts {
                let gap = tools[row.plane].gap(
                    geometry.nodes_m[row.node],
                    step.displacement_m()[row.node],
                    true,
                );
                if !gap.is_finite() {
                    return Err(native_invalid("finite plane contact gap required"));
                }
                paths.push([gap, 0.0, 0.0]);
            }
            retained = Some((forces.to_vec(), paths.clone(), step));
            Ok(paths)
        };
        solve::solve(
            &friction,
            &scales,
            &initial,
            controls.solve,
            &mut response,
            &mut || cx.checkpoint().is_err(),
        )
        .map_err(PlaneContactError::Contact)?
    };
    let (forces, paths, step) = retained.ok_or(PlaneContactError::Invalid(
        "missing final original response",
    ))?;
    if forces != report.forces_n
        || paths != report.paths_m
        || report
            .forces_n
            .iter()
            .any(|force| force[1] != 0.0 || force[2] != 0.0)
    {
        return Err(PlaneContactError::Invalid(
            "contact receipt must match the original frictionless replay",
        ));
    }
    let minimum_after = geometry.plane_gaps(tools, step.displacement_m(), true, cx)?;
    let path_windows = geometry.path(step.displacement_m(), controls, cx)?;
    let old_u = geometry.displacement_m;
    let new_u = step.displacement_m();
    let mut body_work = Dd::ZERO;
    let mut tool_work = Dd::ZERO;
    let mut normal_defect = Dd::ZERO;
    let mut absolute_defect = Dd::ZERO;
    let mut reactions = vec![[Dd::ZERO; 3]; tools.len()];
    for (row, force) in contacts.iter().zip(&report.forces_n) {
        poll(cx)?;
        let tool = tools[row.plane];
        let mut body = Dd::ZERO;
        let mut drive = Dd::ZERO;
        for axis in 0..3 {
            let value = Dd::from_f64(tool.normal[axis] * force[0]);
            body = body
                + value
                    * (Dd::from_f64(new_u[row.node][axis]) - Dd::from_f64(old_u[row.node][axis]));
            drive = drive + value * Dd::from_f64(tool.translation_m[axis]);
            reactions[row.plane][axis] = reactions[row.plane][axis] - value;
        }
        let defect = drive - body;
        body_work = body_work + body;
        tool_work = tool_work + drive;
        normal_defect = normal_defect + defect;
        absolute_defect = absolute_defect + Dd::from_f64(defect.to_f64().abs());
    }
    let absolute_normal_work_defect_j = absolute_defect.to_f64();
    gate(
        "absolute normal endpoint work defect",
        absolute_normal_work_defect_j,
        controls.maximum_normal_work_defect_j,
    )?;
    let dead = geometry.dead_work(nodal, body, new_u, cx)?;
    let decomposition = (Dd::from_f64(step.dead_load_work_j()) - dead - body_work).to_f64();
    gate(
        "native dead/contact load decomposition",
        decomposition,
        controls.work_balance_tolerance_j,
    )?;
    let balance = step.balance(
        dead + Dd::from_f64(step.support_work_j()) + tool_work,
        normal_defect,
    );
    let work_balance_error_j = balance.to_f64();
    gate(
        "tool/material/thermal work identity",
        work_balance_error_j,
        controls.work_balance_tolerance_j,
    )?;
    let force_on_tool_n: Vec<_> = reactions.into_iter().map(|v| v.map(Dd::to_f64)).collect();
    if force_on_tool_n
        .iter()
        .flatten()
        .any(|value| !value.is_finite())
    {
        return Err(PlaneContactError::Invalid(
            "finite tool resultants required",
        ));
    }
    poll(cx)?;
    Ok((
        step,
        PlaneContactReport {
            tools: tools.to_vec(),
            contacts: contacts.to_vec(),
            solve: report,
            gap_before_m: gaps_before,
            body_work_j: body_work.to_f64(),
            tool_work_j: tool_work.to_f64(),
            normal_work_defect_j: normal_defect.to_f64(),
            absolute_normal_work_defect_j,
            dead_load_work_j: dead.to_f64(),
            load_decomposition_error_j: decomposition,
            work_balance_error_j,
            force_on_tool_n,
            minimum_gap_m: minimum_before.min(minimum_after),
            path_windows,
        },
    ))
}

impl TetPlaneGeometry<'_> {
    fn admit<E>(
        &self,
        tools: &[RigidPlaneTool],
        contacts: &[TetPlaneContact],
        prescribed: &[(usize, f64)],
        controls: PlaneContactControls,
    ) -> Result<(), PlaneContactError<E>> {
        if contacts.is_empty()
            || contacts.len() > solve::MAX_CONTACTS
            || contacts.len() > controls.solve.maximum_contacts
            || tools.is_empty()
            || tools.len() > contacts.len()
            || !controls.maximum_normal_work_defect_j.is_finite()
            || controls.maximum_normal_work_defect_j < 0.0
            || !controls.work_balance_tolerance_j.is_finite()
            || controls.work_balance_tolerance_j < 0.0
            || !(1..=1 << 24).contains(&controls.maximum_path_windows)
            || !(1..=1 << 20).contains(&controls.sweep.maximum_windows)
            || !controls.sweep.fraction_tolerance.is_finite()
            || controls.sweep.fraction_tolerance <= 0.0
            || controls.sweep.fraction_tolerance >= 1.0
        {
            return Err(PlaneContactError::Invalid(
                "bounded nonempty tools/contacts, finite nonnegative work gates and valid path controls required",
            ));
        }
        for tool in tools {
            tool.validate()?;
        }
        let mut used = vec![false; tools.len()];
        for (index, row) in contacts.iter().enumerate() {
            if row.node >= self.nodes_m.len()
                || row.plane >= tools.len()
                || !row.initial_force_n.is_finite()
                || row.initial_force_n < 0.0
                || !row.scale_n_per_m.is_finite()
                || row.scale_n_per_m <= 0.0
                || contacts[..index]
                    .iter()
                    .any(|old| old.node == row.node && old.plane == row.plane)
            {
                return Err(PlaneContactError::Invalid(
                    "unique original node/plane rows, finite nonnegative guesses and positive numerical scales required",
                ));
            }
            if !(0..3).any(|axis| {
                tools[row.plane].normal[axis] != 0.0
                    && !prescribed
                        .iter()
                        .any(|&(dof, _)| dof == 3 * row.node + axis)
            }) {
                return Err(PlaneContactError::Invalid(
                    "a contact normal needs a free original nodal coordinate",
                ));
            }
            used[row.plane] = true;
        }
        if used.iter().any(|value| !value) {
            return Err(PlaneContactError::Invalid(
                "every tool needs at least one explicit contact candidate",
            ));
        }
        Ok(())
    }

    fn plane_gaps<E>(
        &self,
        tools: &[RigidPlaneTool],
        displacement: &[[f64; 3]],
        at_end: bool,
        cx: &Cx<'_>,
    ) -> Result<f64, PlaneContactError<E>> {
        let mut minimum = f64::INFINITY;
        for (plane, tool) in tools.iter().enumerate() {
            for (node, &reference) in self.nodes_m.iter().enumerate() {
                if node % 256 == 0 {
                    poll(cx)?;
                }
                let gap_m = tool.gap(reference, displacement[node], at_end);
                if !gap_m.is_finite() || gap_m < 0.0 {
                    return Err(PlaneContactError::Penetration {
                        plane,
                        node,
                        at_end,
                        gap_m,
                    });
                }
                minimum = minimum.min(gap_m);
            }
        }
        Ok(minimum)
    }

    fn path<E>(
        &self,
        after: &[[f64; 3]],
        controls: PlaneContactControls,
        cx: &Cx<'_>,
    ) -> Result<usize, PlaneContactError<E>> {
        let mut windows = 0;
        for (element, tet) in self.tetrahedra.iter().enumerate() {
            poll(cx)?;
            let remaining = controls.maximum_path_windows - windows;
            if remaining == 0 {
                return Err(PlaneContactError::PathBudget);
            }
            let relative = |u: &[[f64; 3]]| {
                tet.map(|node| {
                    core::array::from_fn(|axis| {
                        (Dd::from_f64(self.nodes_m[node][axis])
                            - Dd::from_f64(self.nodes_m[tet[0]][axis])
                            + Dd::from_f64(u[node][axis])
                            - Dd::from_f64(u[tet[0]][axis]))
                        .to_f64()
                    })
                })
            };
            let mut sweep = controls.sweep;
            sweep.maximum_windows = sweep.maximum_windows.min(remaining);
            let result = linear_tetrahedron_orientation(
                relative(self.displacement_m),
                relative(after),
                1.0,
                sweep,
                cx,
            )
            .map_err(PlaneContactError::Geometry)?;
            windows += result.examined_windows;
            if let LinearSweepVerdict::PossibleContact {
                fraction_start,
                fraction_end,
            } = result.verdict
            {
                return Err(PlaneContactError::PossibleElementCollapse {
                    element,
                    fraction_start,
                    fraction_end,
                });
            }
        }
        Ok(windows)
    }

    fn dead_work<E>(
        &self,
        nodal: &[(usize, f64)],
        body: [f64; 3],
        after: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<Dd, PlaneContactError<E>> {
        let delta = |node: usize, axis: usize| {
            Dd::from_f64(after[node][axis]) - Dd::from_f64(self.displacement_m[node][axis])
        };
        let mut work = Dd::ZERO;
        for (i, &(dof, force)) in nodal.iter().enumerate() {
            if i % 256 == 0 {
                poll(cx)?;
            }
            work = work + Dd::from_f64(force) * delta(dof / 3, dof % 3);
        }
        for (i, element) in self.elements.iter().enumerate() {
            if i % 256 == 0 {
                poll(cx)?;
            }
            for &node in &element.nodes {
                for axis in 0..3 {
                    work =
                        work + Dd::from_f64(body[axis] * element.volume / 4.0) * delta(node, axis);
                }
            }
        }
        Ok(work)
    }
}
