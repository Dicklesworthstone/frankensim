//! Frictionless rigid tooling on the original finite thermal plastic response.
//!
//! Fixed-normal planes translate over one physical interval. The original
//! `fs-contact::sustained` natural map solves nonnegative normal forces and
//! endpoint gaps together. Every force probe starts from the SAME accepted
//! material, temperature and momentum state and calls the native solid solve;
//! no penalty spring, surrogate compliance or second return map is introduced.
//!
//! The supplied node/plane rows are explicit candidates. All mesh vertices are
//! checked against every plane before and after the solve, so an omitted row
//! cannot permit an accepted penetration. For P1 cells and translating planes,
//! gap is affine in both material coordinates and the declared straight nodal
//! path. Endpoint gap checks admit that path to their compensated arithmetic
//! precision. The original interval geometry owner separately excludes element
//! collapse throughout this interpolation. This is not an exact-flow CCD claim,
//! a self-contact solver, a curved tool, or an automatic candidate generator.
//!
//! Contact forces are held endpoint forces, including for backward Euler.
//! Their work on the solid and on the moving tool are separate. The difference
//! `force * (gap_before - gap_after)` is a time-discretization defect, bounded
//! by the sum of its absolute values, never plastic heat. Refine the physical
//! interval at contact arrival; no impact or contact-duration law is inferred.

use super::{
    AffineThermalFiniteJ2, FiniteThermalJ2Material, FiniteThermoplasticControls,
    FiniteThermoplasticTetError, FiniteThermoplasticTetSolid, FiniteThermoplasticTetStep,
    HyperTetError, HyperTetSettings, invalid,
};
use fs_contact::sustained::NaturalMapScales;
use fs_contact::sustained::solve::{self, SolveControls, SolveError, SolveReport};
use fs_exec::Cx;
use fs_math::dd::Dd;
use fs_query::TriangleQueryError;
use fs_query::triangle::sweep::{
    LinearSweepControls, LinearSweepVerdict, linear_tetrahedron_orientation,
};

/// An infinite rigid half-space with explicitly prescribed translation.
/// Free space has `normal . (x - point_before - translation) >= 0` at the end.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RigidPlaneTool {
    /// Unit normal pointing from the rigid tool into the admissible solid side.
    pub normal: [f64; 3],
    /// A point on the plane at the accepted beginning of this interval, m.
    pub point_before_m: [f64; 3],
    /// Tool translation over the SAME physical interval as the solid, m.
    /// Tangential translation does no work in this frictionless model.
    pub translation_m: [f64; 3],
}

/// One original mesh vertex allowed to develop a reaction against one plane.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TetPlaneContact {
    pub node: usize,
    pub plane: usize,
    /// Nonnegative initial force guess, N; it does not change physical history.
    pub initial_force_n: f64,
    /// Positive natural-map conditioning scale, N/m; not a material stiffness.
    pub scale_n_per_m: f64,
}

/// Native force-solve budgets and independent geometric/work admission.
#[derive(Clone, Copy, Debug)]
pub struct PlaneContactControls {
    pub solve: SolveControls,
    /// Sum of absolute normal endpoint-work defects, J. Reduce dt/motion on refusal.
    pub maximum_normal_work_defect_j: f64,
    /// Absolute work decomposition and full mechanical/thermal identity gate, J.
    pub work_balance_tolerance_j: f64,
    /// Original interval tetrahedron-orientation controls, applied per element.
    pub sweep: LinearSweepControls,
    /// Total interval windows over the original mesh; no partial path is accepted.
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

/// Original contact solution, work ports, and complete declared path admission.
#[derive(Clone, Debug, PartialEq)]
pub struct PlaneContactReport {
    pub tools: Vec<RigidPlaneTool>,
    pub contacts: Vec<TetPlaneContact>,
    /// Original contact owner receipt, in caller row order. Tangent entries are zero.
    pub solve: SolveReport,
    pub gap_before_m: Vec<f64>,
    /// Actual held normal force dotted with each original nodal displacement change.
    pub body_work_j: f64,
    /// Actual held normal force dotted with prescribed tool translation.
    pub tool_work_j: f64,
    /// Signed tool work minus body work; never added to caloric heat.
    pub normal_work_defect_j: f64,
    /// Noncancelling gate quantity, including inactive-row complementarity error.
    pub absolute_normal_work_defect_j: f64,
    /// Work of the caller's original dead nodal/body loads, excluding contact.
    pub dead_load_work_j: f64,
    /// Native assembled dead-load work minus original-load and contact work.
    pub load_decomposition_error_j: f64,
    /// Complete external work minus native storage, heat, inertia and defects.
    pub work_balance_error_j: f64,
    /// Equal/opposite resultant on each tool, in original plane order, N.
    pub force_on_tool_n: Vec<[f64; 3]>,
    /// Smallest beginning/end gap over ALL original vertices and planes, m.
    /// This is a compensated numerical gap, not a rigorous distance enclosure.
    pub minimum_gap_m: f64,
    /// Interval windows spent proving positive element orientation on the chord.
    pub path_windows: usize,
}

/// One immutable complete candidate; the native material and contact accounts coexist.
#[derive(Debug)]
pub struct FiniteThermoplasticPlaneContactStep<
    'a,
    L: FiniteThermalJ2Material = AffineThermalFiniteJ2,
> {
    /// Original native result. Its dead-load work includes the applied contact
    /// forces; use `contact.dead_load_work_j` when counting tool work separately.
    pub step: FiniteThermoplasticTetStep<'a, L>,
    pub contact: PlaneContactReport,
}

/// Refusal leaves all accepted deformation, plastic, thermal and momentum state intact.
#[derive(Debug)]
pub enum FiniteThermoplasticPlaneContactError {
    Invalid(&'static str),
    Cancelled,
    Solid(FiniteThermoplasticTetError),
    Contact(SolveError<FiniteThermoplasticTetError>),
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

impl std::fmt::Display for FiniteThermoplasticPlaneContactError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "finite thermoplastic plane contact refused: {self:?}")
    }
}
impl std::error::Error for FiniteThermoplasticPlaneContactError {}

type Error = FiniteThermoplasticPlaneContactError;

fn poll(cx: &Cx<'_>) -> Result<(), Error> {
    cx.checkpoint().map_err(|_| Error::Cancelled)
}

fn gate(quantity: &'static str, value_j: f64, maximum_abs_j: f64) -> Result<(), Error> {
    if !value_j.is_finite() || value_j.abs() > maximum_abs_j {
        Err(Error::Energy {
            quantity,
            value_j,
            maximum_abs_j,
        })
    } else {
        Ok(())
    }
}

impl RigidPlaneTool {
    fn validate(self) -> Result<(), Error> {
        let norm = self.normal.into_iter().fold(0.0_f64, |v, x| v.hypot(x));
        if self
            .normal
            .into_iter()
            .chain(self.point_before_m)
            .chain(self.translation_m)
            .any(|x| !x.is_finite())
            || (norm - 1.0).abs() > 128.0 * f64::EPSILON
        {
            return Err(Error::Invalid(
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

impl<'a, L: FiniteThermalJ2Material> FiniteThermoplasticTetSolid<'a, L> {
    /// Propose finite-duration deformation against frictionless translating tools.
    ///
    /// All force responses retain the original accepted history, source families,
    /// physical dt and cached reference mass. Native force and material energy
    /// gates apply to every trial. A final fresh response must match the contact
    /// receipt exactly before full-mesh geometry and independent work admission.
    /// No accepted state, source, node or gap is projected or repaired.
    ///
    /// The caller supplies at least one candidate per tool, at most the original
    /// contact owner's 21-row cap. Omitted rows are permitted only while every
    /// original vertex remains on the free side. This is an endpoint method:
    /// source/history validity between endpoints and physical arrival timing
    /// require refinement, beyond the admitted affine geometry.
    ///
    /// # Errors
    /// Original solid/contact/geometry refusal, incomplete geometry candidates,
    /// malformed controls, exceeded work/path budgets, or cancellation.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn advance_with_plane_contact(
        &self,
        new_temperatures_k: &[f64],
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        dt_s: f64,
        tools: &[RigidPlaneTool],
        contacts: &[TetPlaneContact],
        settings: HyperTetSettings,
        material_controls: FiniteThermoplasticControls,
        controls: PlaneContactControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteThermoplasticPlaneContactStep<'a, L>, Error> {
        poll(cx)?;
        self.admit_plane_contact(tools, contacts, prescribed_m, controls)?;
        let minimum_before = self.plane_gaps(tools, self.displacement_m(), false, cx)?;
        let gaps_before: Vec<_> = contacts
            .iter()
            .map(|row| {
                tools[row.plane].gap(
                    self.state.nodes_m[row.node],
                    self.displacement_m()[row.node],
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
             -> Result<Vec<[f64; 3]>, FiniteThermoplasticTetError> {
                // Retain only the latest successful response, never a history
                // from a previous Newton force or a discarded material trial.
                drop(retained.take());
                if cancelled() {
                    return Err(HyperTetError::Cancelled.into());
                }
                let mut applied = Vec::with_capacity(nodal_forces_n.len() + 3 * contacts.len());
                applied.extend_from_slice(nodal_forces_n);
                for (row, force) in contacts.iter().zip(forces) {
                    for axis in 0..3 {
                        let value = tools[row.plane].normal[axis] * force[0];
                        if !value.is_finite() {
                            return Err(invalid("finite contact nodal force required").into());
                        }
                        if value != 0.0 {
                            applied.push((3 * row.node + axis, value));
                        }
                    }
                }
                let step = self.advance(
                    new_temperatures_k,
                    prescribed_m,
                    &applied,
                    body_force_n_m3,
                    dt_s,
                    settings,
                    material_controls,
                    cx,
                )?;
                let mut paths = Vec::with_capacity(contacts.len());
                for row in contacts {
                    let gap = tools[row.plane].gap(
                        self.state.nodes_m[row.node],
                        step.solid.displacement_m()[row.node],
                        true,
                    );
                    if !gap.is_finite() {
                        return Err(invalid("finite plane contact gap required").into());
                    }
                    // The two auxiliary tangent equations are exactly ft=0.
                    // They introduce no physical tangent load or slip law.
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
            .map_err(Error::Contact)?
        };
        let (forces, paths, step) =
            retained.ok_or(Error::Invalid("missing final original response"))?;
        if forces != report.forces_n
            || paths != report.paths_m
            || report
                .forces_n
                .iter()
                .any(|force| force[1] != 0.0 || force[2] != 0.0)
        {
            return Err(Error::Invalid(
                "contact receipt must match the original frictionless replay",
            ));
        }
        let minimum_after = self.plane_gaps(tools, step.solid.displacement_m(), true, cx)?;
        let path_windows = self.plane_contact_path(&step.solid, controls, cx)?;
        let old_u = self.displacement_m();
        let new_u = step.solid.displacement_m();
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
                // Use the same force product supplied to the original solid.
                let value = Dd::from_f64(tool.normal[axis] * force[0]);
                body = body
                    + value
                        * (Dd::from_f64(new_u[row.node][axis])
                            - Dd::from_f64(old_u[row.node][axis]));
                drive = drive + value * Dd::from_f64(tool.translation_m[axis]);
                reactions[row.plane][axis] = reactions[row.plane][axis] - value;
            }
            let defect = drive - body;
            body_work = body_work + body;
            tool_work = tool_work + drive;
            normal_defect = normal_defect + defect;
            absolute_defect = absolute_defect + Dd::from_f64(defect.to_f64().abs());
        }
        let body_work_j = body_work.to_f64();
        let tool_work_j = tool_work.to_f64();
        let normal_work_defect_j = normal_defect.to_f64();
        let absolute_normal_work_defect_j = absolute_defect.to_f64();
        gate(
            "absolute normal endpoint work defect",
            absolute_normal_work_defect_j,
            controls.maximum_normal_work_defect_j,
        )?;
        let dead = self.original_dead_work(nodal_forces_n, body_force_n_m3, new_u, cx)?;
        let dead_load_work_j = dead.to_f64();
        let decomposition =
            (Dd::from_f64(step.energy.dead_load_endpoint_work_j) - dead - body_work).to_f64();
        gate(
            "native dead/contact load decomposition",
            decomposition,
            controls.work_balance_tolerance_j,
        )?;
        let mut balance = dead + Dd::from_f64(step.energy.support_endpoint_work_j) + tool_work
            - Dd::from_f64(step.energy.internal_energy_change_j)
            - Dd::from_f64(step.energy.caloric_heat_j)
            - Dd::from_f64(step.energy.endpoint_work_defect_j)
            - normal_defect;
        if let Some(inertia) = &step.inertia {
            balance = balance
                - Dd::from_f64(inertia.kinetic_change_j)
                - Dd::from_f64(inertia.numerical_loss_j);
        }
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
            return Err(Error::Invalid("finite tool resultants required"));
        }
        poll(cx)?;
        Ok(FiniteThermoplasticPlaneContactStep {
            step,
            contact: PlaneContactReport {
                tools: tools.to_vec(),
                contacts: contacts.to_vec(),
                solve: report,
                gap_before_m: gaps_before,
                body_work_j,
                tool_work_j,
                normal_work_defect_j,
                absolute_normal_work_defect_j,
                dead_load_work_j,
                load_decomposition_error_j: decomposition,
                work_balance_error_j,
                force_on_tool_n,
                minimum_gap_m: minimum_before.min(minimum_after),
                path_windows,
            },
        })
    }

    fn admit_plane_contact(
        &self,
        tools: &[RigidPlaneTool],
        contacts: &[TetPlaneContact],
        prescribed: &[(usize, f64)],
        controls: PlaneContactControls,
    ) -> Result<(), Error> {
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
            return Err(Error::Invalid(
                "bounded nonempty tools/contacts, finite nonnegative work gates and valid path controls required",
            ));
        }
        for tool in tools {
            tool.validate()?;
        }
        let mut used = vec![false; tools.len()];
        for (index, row) in contacts.iter().enumerate() {
            if row.node >= self.state.nodes_m.len()
                || row.plane >= tools.len()
                || !row.initial_force_n.is_finite()
                || row.initial_force_n < 0.0
                || !row.scale_n_per_m.is_finite()
                || row.scale_n_per_m <= 0.0
                || contacts[..index]
                    .iter()
                    .any(|old| old.node == row.node && old.plane == row.plane)
            {
                return Err(Error::Invalid(
                    "unique original node/plane rows, finite nonnegative guesses and positive numerical scales required",
                ));
            }
            if !(0..3).any(|axis| {
                tools[row.plane].normal[axis] != 0.0
                    && !prescribed
                        .iter()
                        .any(|&(dof, _)| dof == 3 * row.node + axis)
            }) {
                return Err(Error::Invalid(
                    "a contact normal needs a free original nodal coordinate",
                ));
            }
            used[row.plane] = true;
        }
        if used.iter().any(|value| !value) {
            return Err(Error::Invalid(
                "every tool needs at least one explicit contact candidate",
            ));
        }
        Ok(())
    }

    fn plane_gaps(
        &self,
        tools: &[RigidPlaneTool],
        displacement: &[[f64; 3]],
        at_end: bool,
        cx: &Cx<'_>,
    ) -> Result<f64, Error> {
        let mut minimum = f64::INFINITY;
        for (plane, tool) in tools.iter().enumerate() {
            for (node, &reference) in self.state.nodes_m.iter().enumerate() {
                if node % 256 == 0 {
                    poll(cx)?;
                }
                let gap_m = tool.gap(reference, displacement[node], at_end);
                if !gap_m.is_finite() || gap_m < 0.0 {
                    return Err(Error::Penetration {
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

    fn plane_contact_path(
        &self,
        after: &Self,
        controls: PlaneContactControls,
        cx: &Cx<'_>,
    ) -> Result<usize, Error> {
        let mut windows = 0;
        for (element, tet) in self.state.tetrahedra.iter().enumerate() {
            poll(cx)?;
            let remaining = controls.maximum_path_windows - windows;
            if remaining == 0 {
                return Err(Error::PathBudget);
            }
            // Removing the same affine anchor from every vertex preserves
            // orientation and avoids losing edges to a large rigid translation.
            let relative = |u: &[[f64; 3]]| {
                tet.map(|node| {
                    core::array::from_fn(|axis| {
                        (Dd::from_f64(self.state.nodes_m[node][axis])
                            - Dd::from_f64(self.state.nodes_m[tet[0]][axis])
                            + Dd::from_f64(u[node][axis])
                            - Dd::from_f64(u[tet[0]][axis]))
                        .to_f64()
                    })
                })
            };
            let mut sweep = controls.sweep;
            sweep.maximum_windows = sweep.maximum_windows.min(remaining);
            // The original native tetrahedron owner admits positive reference
            // orientation only. The query retains that immutable sign.
            let result = linear_tetrahedron_orientation(
                relative(self.displacement_m()),
                relative(after.displacement_m()),
                1.0,
                sweep,
                cx,
            )
            .map_err(Error::Geometry)?;
            windows += result.examined_windows;
            if let LinearSweepVerdict::PossibleContact {
                fraction_start,
                fraction_end,
            } = result.verdict
            {
                return Err(Error::PossibleElementCollapse {
                    element,
                    fraction_start,
                    fraction_end,
                });
            }
        }
        Ok(windows)
    }

    fn original_dead_work(
        &self,
        nodal: &[(usize, f64)],
        body: [f64; 3],
        after: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<Dd, Error> {
        let before = self.displacement_m();
        let delta = |node: usize, axis: usize| {
            Dd::from_f64(after[node][axis]) - Dd::from_f64(before[node][axis])
        };
        let mut work = Dd::ZERO;
        for (i, &(dof, force)) in nodal.iter().enumerate() {
            if i % 256 == 0 {
                poll(cx)?;
            }
            // The original response has already admitted these node-major DOFs.
            work = work + Dd::from_f64(force) * delta(dof / 3, dof % 3);
        }
        for (i, element) in self.state.elements.iter().enumerate() {
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
