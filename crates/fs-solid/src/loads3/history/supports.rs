//! Typed support motion sampled at the same physical time as distributed loads.
//! Constraint components are WORLD xyz, fixed in the reference discretization.
//! Rotating payload frames rotate the prescribed vector, not the restraint axes.
//! This is a sequence of absolute quasi-static endpoints, not base-excited dynamics.
use super::*;
use std::collections::{BTreeMap, BTreeSet};

/// SI metres, with no implicit velocity integration or frame-translation term.
pub const DISPLACEMENT_DIMS: Dims = Dims([1, 0, 0, 0, 0, 0]);

/// One spatially uniform displacement history on explicitly selected nodes.
/// Only the selected WORLD components are prescribed. Others remain governed
/// by the original problem; an omitted fixed component stays zero. No component
/// may be assigned twice, even with equal values: ownership is unambiguous.
#[derive(Debug, Clone, Copy)]
pub struct SupportMotion3<'a> {
    pub nodes: &'a [usize],
    pub components: [bool; 3],
    pub displacement: &'a VectorPayload,
}

/// One retained endpoint with its exact problem, forces and support values.
/// The force and support histories cannot be resampled independently at solve
/// time, and material density cannot change after acceleration integration.
pub struct TetSupportSample3<'g, 'm, 'p, 's> {
    loads: TetLoadSample3<'g, 'm, 'p, 's>,
    prescribed_m: Vec<(usize, f64)>,
    maximum_displacement_m: f64,
}
impl TetSupportSample3<'_, '_, '_, '_> {
    pub fn time_s(&self) -> f64 { self.loads.time_s() }
    /// Sorted full xyz DOFs and absolute displacement values [m].
    pub fn prescribed_dofs_m(&self) -> &[(usize, f64)] { &self.prescribed_m }
    pub fn nodal_forces_n(&self) -> &[(usize, f64)] { self.loads.nodal_forces_n() }
    pub fn resultant(&self, origin_m: [f64; 3], cx: &Cx<'_>)
        -> Result<([f64; 3], [f64; 3]), LoadError3> { self.loads.resultant(origin_m, cx) }

    /// Preserve the original thermal/material state, support reactions and true
    /// residual gate; no endpoint becomes a correction-only displacement field.
    pub fn solve(&self, thermal: TetThermalStrainField<'_>, config: TetStaticSolveConfig,
        mesh_frame: ContentHash, maximum_strain: f64, cx: &Cx<'_>)
        -> Result<TetLoadSolution3, LoadError3> {
        self.loads.case.solve_with_prescribed_dofs(self.loads.problem, thermal,
            &self.prescribed_m, self.maximum_displacement_m, config, mesh_frame, maximum_strain, cx)
    }
}

impl TetLoadHistory3<'_> {
    /// Evaluate loads and support motion at ONE physical `time_s`, using the
    /// same declared frame tree and existing payload interpolation/endpoint rules.
    /// Even-parity metre vectors are required; spatial profiles, reference resets,
    /// distributions and stronger unhandled semantic contracts are not guessed.
    /// Components select world axes AFTER vector transformation, never moving
    /// constraint axes or an implicit projection. The bound applies separately
    /// to each selected displacement component, not to unselected components.
    ///
    /// At most 1024 support rows and at most the problem's fixed-DOF count of
    /// expanded assignments are admitted before binding allocation. Every node,
    /// DOF and source is checked; duplicate or unconstrained assignments refuse.
    /// Partial loads/supports never escape. Checkpoints precede publication and
    /// occur per row, node and fixed DOF; existing frame/signal kernels keep their
    /// own indivisible work boundaries. No time integration or clock is advanced.
    #[allow(clippy::too_many_arguments)]
    pub fn at_with_supports<'g, 'm, 'p, 's>(
        &self, geometry: &'g TetLoadGeometry3<'m>, problem: &'p TetLinearElasticProblem<'s>,
        supports: &[SupportMotion3<'_>], time_s: f64, maximum_displacement_m: f64, cx: &Cx<'_>,
    ) -> Result<TetSupportSample3<'g, 'm, 'p, 's>, LoadHistoryError3> {
        use LoadHistoryError3 as E;
        poll(cx)?;
        if !time_s.is_finite() || !maximum_displacement_m.is_finite() || maximum_displacement_m <= 0.0 {
            return Err(E::Invalid("finite time and positive finite displacement-component bound required"));
        }
        if supports.len() > 1024 || problem.nodes_m.len() > problem.budget.maximum_nodes {
            return Err(E::Invalid("support rows or problem nodes exceed the admitted envelope"));
        }
        let full_dofs = problem.nodes_m.len().checked_mul(3)
            .ok_or(E::Invalid("full support-DOF count overflows"))?;
        if problem.fixed_dofs.len() > full_dofs { return Err(E::Invalid("too many fixed components")); }
        let mut count = 0usize;
        for row in supports {
            poll(cx)?;
            let components = row.components.iter().filter(|v| **v).count();
            if row.nodes.is_empty() || components == 0 {
                return Err(E::Invalid("each support history needs nodes and at least one world component"));
            }
            count = row.nodes.len().checked_mul(components).and_then(|n| count.checked_add(n))
                .ok_or(E::Invalid("expanded support-DOF count overflows"))?;
            if count > problem.fixed_dofs.len() {
                return Err(E::Invalid("support assignments exceed the declared fixed-component count"));
            }
        }
        let mut fixed = BTreeSet::new();
        for &dof in problem.fixed_dofs {
            poll(cx)?;
            if dof >= full_dofs || !fixed.insert(dof) {
                return Err(E::Invalid("fixed components must be unique and inside the reference mesh"));
            }
        }
        let mut prescribed = BTreeMap::new();
        for row in supports {
            poll(cx)?;
            let value = sample(row.displacement, DISPLACEMENT_DIMS, self.frames, time_s, cx)?;
            for c in 0..3 {
                if row.components[c] && value[c].abs() > maximum_displacement_m {
                    return Err(E::Invalid("sampled support motion exceeds its displacement-component bound"));
                }
            }
            for &node in row.nodes {
                poll(cx)?;
                if node >= problem.nodes_m.len() { return Err(E::Invalid("support node is outside the reference mesh")); }
                for c in 0..3 {
                    if !row.components[c] { continue; }
                    let dof = 3 * node + c;
                    if !fixed.contains(&dof) { return Err(E::Invalid("support history names a component that is not fixed in the problem")); }
                    if prescribed.insert(dof, value[c]).is_some() {
                        return Err(E::Invalid("duplicate support-component assignment"));
                    }
                }
            }
        }
        let mut prescribed_m = Vec::with_capacity(count);
        for value in prescribed { poll(cx)?; prescribed_m.push(value); }
        let loads = self.at(geometry, problem, time_s, cx)?;
        poll(cx)?;
        Ok(TetSupportSample3 { loads, prescribed_m, maximum_displacement_m })
    }
}
