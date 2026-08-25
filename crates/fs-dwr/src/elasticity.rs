//! Goal-oriented compliance estimation for vector CutFEM elasticity.
//!
//! The compliance goal is the assembled algebraic work `J_h = b_h^T u_h`.
//! Because the vector operator is symmetric, the discrete compliance adjoint
//! is the primal solution itself.  The enriched weight is therefore
//! `w = u_{h/2} - u_h`, evaluated pointwise without missing-node fallbacks.
//!
//! Bulk, embedded-Nitsche, and outer-traction terms use the exact signs of the
//! coarse weak residual while integrating over the enriched cut partition. The
//! ghost term uses the coarse consistent-limit correction
//! `+g_h(u_h, u_h)`.  It is the coarse ghost residual against the limiting
//! smooth adjoint, whose normal-derivative jump vanishes.  This avoids
//! inventing enriched traces on inactive halves of a coarse ghost face.  It is
//! measured estimator evidence, not a certified error bound.

use fs_cutfem::quad::tensor_gauss;
use fs_cutfem::{
    BoundaryTraction, CellKey, CutElasticity, CutElasticitySolution, CutFemError, DesignBoxEdge,
    Quadtree, SharedFacePatch, supported_traction_segment,
};
use std::collections::{BTreeMap, BTreeSet};

/// The stabilization convention used by [`estimate_elasticity_compliance`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ElasticityGhostMethod {
    /// Evaluate `+g_h(u_h, u_h)` on the actual coarse ghost-face set and
    /// allocate each non-negative face contribution equally to its cells.
    CoarseConsistentEnergy,
}

impl ElasticityGhostMethod {
    /// Stable evidence label for JSONL rows and ledgers.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::CoarseConsistentEnergy => "coarse-consistent-energy",
        }
    }
}

/// Signed residual decomposition.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ElasticityResidualTerms {
    /// `sum_K integral(f dot w - sigma(u_h):epsilon(w))`.
    pub bulk: f64,
    /// Symmetric embedded-Nitsche residual, including inhomogeneous `g`.
    pub nitsche: f64,
    /// Dead-load work `integral(t_bar dot w)` on the design-box boundary.
    pub outer_traction: f64,
    /// Coarse consistent-limit stabilization correction described by the
    /// estimate's [`ElasticityGhostMethod`].
    pub ghost: f64,
}

impl ElasticityResidualTerms {
    /// Deterministic signed sum of all reported terms.
    #[must_use]
    pub fn total(self) -> f64 {
        ((self.bulk + self.nitsche) + self.outer_traction) + self.ghost
    }
}

/// Compliance DWR result on one uniform or 2:1-graded vector CutFEM grid.
#[derive(Debug, Clone)]
pub struct ElasticityDwrEstimate {
    /// Signed cell reconstruction of the estimated goal error.
    pub eta_signed: f64,
    /// `sum_K |eta_K|`, used only as deterministic marking mass.
    pub eta_abs: f64,
    /// Signed indicator per coarse active cell.
    pub indicators: BTreeMap<CellKey, f64>,
    /// Per-cell decomposition whose deterministic sums reconstruct [`Self::terms`].
    /// The aggregate [`Self::indicators`] remain the marking signal; this map
    /// makes term-specific scientific claims auditable without discarding any
    /// contribution from the actual estimator.
    pub cell_terms: BTreeMap<CellKey, ElasticityResidualTerms>,
    /// Coarse consistent-limit ghost correction per canonical coarse face.
    pub face_indicators: BTreeMap<(CellKey, CellKey), f64>,
    /// Signed decomposition before cell allocation.
    pub terms: ElasticityResidualTerms,
    /// Stabilization convention that produced `terms.ghost`.
    pub ghost_method: ElasticityGhostMethod,
    /// Coarse assembled-load compliance `b_h^T u_h`.
    pub j_primal: f64,
    /// Enriched assembled-load compliance `b_{h/2}^T u_{h/2}`.
    pub j_enriched: f64,
    /// Coarse vector displacement DOF count.
    pub dofs: usize,
    /// Enriched vector displacement DOF count.
    pub enriched_dofs: usize,
}

/// Estimate compliance error using self-adjoint coarse and enriched primal
/// solves with uncertified callback traction or no traction.
///
/// `J=b^T u` is physical dead-load compliance when essential data are
/// homogeneous. With nonzero embedded data `g`, it is explicitly the
/// algebraic assembled-load compliance, not a physical external-work claim.
///
/// # Errors
///
/// Propagates vector CutFEM build/solve refusals. Field evaluation fails closed
/// when a requested active cell or corner value is unavailable, and every
/// callback-derived or accumulated quantity must remain finite.
pub fn estimate_elasticity_compliance(
    problem: &CutElasticity<'_>,
    body: &dyn Fn(f64, f64) -> [f64; 2],
    embedded_data: &dyn Fn(f64, f64) -> [f64; 2],
) -> Result<ElasticityDwrEstimate, CutFemError> {
    estimate_elasticity_compliance_impl(problem, body, embedded_data, None)
}

/// Estimate compliance error using self-adjoint coarse and enriched primal
/// solves with typed boundary traction data.
///
/// Requires [`CutElasticity::boundary_traction`] to be `None`.
///
/// # Errors
///
/// Propagates vector CutFEM build/solve refusals or invalid traction inputs.
pub fn estimate_elasticity_compliance_with_boundary_traction(
    problem: &CutElasticity<'_>,
    body: &dyn Fn(f64, f64) -> [f64; 2],
    embedded_data: &dyn Fn(f64, f64) -> [f64; 2],
    boundary_traction: BoundaryTraction<'_>,
) -> Result<ElasticityDwrEstimate, CutFemError> {
    if problem.boundary_traction.is_some() {
        return Err(CutFemError::InvalidElasticityInput {
            what: "estimate_elasticity_compliance_with_boundary_traction requires CutElasticity::boundary_traction to be None"
                .to_string(),
        });
    }
    estimate_elasticity_compliance_impl(problem, body, embedded_data, Some(boundary_traction))
}

fn estimate_elasticity_compliance_impl(
    problem: &CutElasticity<'_>,
    body: &dyn Fn(f64, f64) -> [f64; 2],
    embedded_data: &dyn Fn(f64, f64) -> [f64; 2],
    boundary_traction: Option<BoundaryTraction<'_>>,
) -> Result<ElasticityDwrEstimate, CutFemError> {
    if problem.grid.max_level() >= 16 {
        return Err(invalid(
            "elasticity DWR enrichment requires one level of lattice headroom below level 16"
                .to_string(),
        ));
    }
    let coarse = match boundary_traction {
        Some(bt) => problem.solve_with_boundary_traction(body, embedded_data, bt)?,
        None => problem.solve(body, embedded_data)?,
    };
    let fine_grid = problem.grid.refined_once();
    let fine_problem = CutElasticity {
        grid: &fine_grid,
        sdf: problem.sdf,
        material: problem.material,
        nitsche_beta: problem.nitsche_beta,
        ghost_gamma: problem.ghost_gamma,
        stabilization_scaling: problem.stabilization_scaling,
        quad_depth: problem.quad_depth,
        clamp: problem.clamp,
        boundary_traction: if boundary_traction.is_some() {
            None
        } else {
            problem.boundary_traction
        },
        traction_free_interface: problem.traction_free_interface,
        solver_tol: problem.solver_tol,
        solver_max_iters: problem.solver_max_iters,
    };
    let fine = match boundary_traction {
        Some(bt) => fine_problem.solve_with_boundary_traction(body, embedded_data, bt)?,
        None => fine_problem.solve(body, embedded_data)?,
    };
    let effective_traction =
        boundary_traction.or_else(|| problem.boundary_traction.map(BoundaryTraction::Uncertified));
    let (lambda, mu) = problem.material.lame();
    let modulus = problem.effective_stabilization_modulus(lambda, mu)?;

    let coarse_active: BTreeSet<CellKey> = coarse.active_cells().iter().copied().collect();
    let fine_active: BTreeSet<CellKey> = fine.active_cells().iter().copied().collect();
    let mut indicators: BTreeMap<CellKey, f64> = coarse_active
        .iter()
        .copied()
        .map(|cell| (cell, 0.0))
        .collect();
    let mut cell_terms: BTreeMap<CellKey, ElasticityResidualTerms> = coarse_active
        .iter()
        .copied()
        .map(|cell| (cell, ElasticityResidualTerms::default()))
        .collect();
    let mut enriched_parents = BTreeSet::new();
    let mut terms = ElasticityResidualTerms::default();

    for fine_cell in fine_active.iter().copied() {
        let parent = coarse_parent(problem.grid, &fine_grid, fine_cell, &coarse_active)?;
        enriched_parents.insert(parent);
        let (lo, hi) = fine_grid.rect(fine_cell);
        let bulk_rule;
        let bulk: &[([f64; 2], f64)] = if let Some(rule) = fine.cut_rules().get(&fine_cell) {
            &rule.bulk
        } else {
            bulk_rule = {
                let mut points = Vec::with_capacity(9);
                tensor_gauss(lo, hi, &mut points);
                points
            };
            &bulk_rule
        };
        let mut cell_bulk = 0.0;
        for &(point, weight) in bulk {
            let (fine_value, fine_gradient) = fine.value_gradient(&fine_grid, fine_cell, point)?;
            let (coarse_value, coarse_gradient) =
                coarse.value_gradient(problem.grid, parent, point)?;
            let weight_value = sub2(fine_value, coarse_value);
            let weight_gradient = sub_gradient(fine_gradient, coarse_gradient);
            let force = body(point[0], point[1]);
            require_vector_finite(force, "body force in elasticity DWR bulk term")?;
            let stress = stress(lambda, mu, coarse_gradient);
            let value = weight
                * (dot2(force, weight_value) - stress_strain_contract(stress, weight_gradient));
            require_finite(value, "elasticity DWR bulk contribution")?;
            cell_bulk += value;
        }
        require_finite(cell_bulk, "elasticity DWR cell bulk sum")?;
        terms.bulk += cell_bulk;
        add_cell_term(&mut cell_terms, parent, |cell| cell.bulk += cell_bulk)?;
        add_indicator(&mut indicators, parent, cell_bulk)?;

        if !problem.traction_free_interface
            && let Some(rule) = fine.cut_rules().get(&fine_cell)
        {
            let coarse_h = problem.grid.cell_h(parent);
            let penalty = problem.nitsche_beta * modulus / coarse_h;
            require_finite(penalty, "elasticity DWR Nitsche penalty")?;
            let mut cell_nitsche = 0.0;
            for &(point, weight, normal) in &rule.iface {
                let (fine_value, fine_gradient) =
                    fine.value_gradient(&fine_grid, fine_cell, point)?;
                let (coarse_value, coarse_gradient) =
                    coarse.value_gradient(problem.grid, parent, point)?;
                let weight_value = sub2(fine_value, coarse_value);
                let weight_gradient = sub_gradient(fine_gradient, coarse_gradient);
                let data = embedded_data(point[0], point[1]);
                require_vector_finite(data, "embedded data in elasticity DWR Nitsche term")?;
                let displacement_gap = sub2(coarse_value, data);
                let coarse_traction = traction(stress(lambda, mu, coarse_gradient), normal);
                let weight_traction = traction(stress(lambda, mu, weight_gradient), normal);
                let value = weight
                    * (dot2(coarse_traction, weight_value)
                        + dot2(weight_traction, displacement_gap)
                        - penalty * dot2(displacement_gap, weight_value));
                require_finite(value, "elasticity DWR Nitsche contribution")?;
                cell_nitsche += value;
            }
            require_finite(cell_nitsche, "elasticity DWR cell Nitsche sum")?;
            terms.nitsche += cell_nitsche;
            add_cell_term(&mut cell_terms, parent, |cell| {
                cell.nitsche += cell_nitsche;
            })?;
            add_indicator(&mut indicators, parent, cell_nitsche)?;
        }

        let cell_traction = outer_traction_residual(
            problem,
            &coarse,
            &fine,
            &fine_grid,
            fine_cell,
            parent,
            effective_traction,
        )?;
        terms.outer_traction += cell_traction;
        add_cell_term(&mut cell_terms, parent, |cell| {
            cell.outer_traction += cell_traction;
        })?;
        add_indicator(&mut indicators, parent, cell_traction)?;
    }

    if let Some(missing) = coarse_active.difference(&enriched_parents).next() {
        return Err(invalid(format!(
            "coarse active cell {missing:?} has no enriched active child; non-nested active-domain loss is unsupported"
        )));
    }

    require_finite(terms.bulk, "elasticity DWR global bulk sum")?;
    require_finite(terms.nitsche, "elasticity DWR global Nitsche sum")?;
    require_finite(
        terms.outer_traction,
        "elasticity DWR global outer-traction sum",
    )?;

    let mut face_indicators = BTreeMap::new();
    if problem.ghost_gamma > 0.0 {
        for &(cell_a, cell_b) in coarse.ghost_faces() {
            let patch = problem
                .grid
                .shared_face_patch(cell_a, cell_b)
                .map_err(|error| {
                    invalid(format!(
                        "cannot reconstruct coarse elasticity ghost patch ({cell_a:?}, {cell_b:?}): {error}"
                    ))
                })?;
            let face_value = coarse_ghost_consistent_energy(problem, &coarse, patch, modulus)?;
            if face_indicators
                .insert((cell_a, cell_b), face_value)
                .is_some()
            {
                return Err(invalid(format!(
                    "duplicate canonical coarse ghost face ({cell_a:?}, {cell_b:?})"
                )));
            }
            terms.ghost += face_value;
            let half = 0.5 * face_value;
            add_cell_term(&mut cell_terms, cell_a, |cell| cell.ghost += half)?;
            add_cell_term(&mut cell_terms, cell_b, |cell| cell.ghost += half)?;
            add_indicator(&mut indicators, cell_a, half)?;
            add_indicator(&mut indicators, cell_b, half)?;
        }
    }
    require_finite(terms.ghost, "elasticity DWR coarse ghost-energy sum")?;

    for (&cell, &indicator) in &indicators {
        require_finite(
            indicator,
            &format!("elasticity DWR indicator for cell {cell:?}"),
        )?;
    }
    let eta_signed: f64 = indicators.values().sum();
    let eta_abs: f64 = indicators.values().map(|value| value.abs()).sum();
    require_finite(eta_signed, "elasticity DWR signed estimate")?;
    require_finite(eta_abs, "elasticity DWR marking mass")?;
    let reconstruction_magnitude = eta_abs.max(
        terms.bulk.abs() + terms.nitsche.abs() + terms.outer_traction.abs() + terms.ghost.abs(),
    );
    require_reconstruction(
        "elasticity DWR total cell allocation",
        eta_signed,
        terms.total(),
        reconstruction_magnitude,
    )?;
    validate_cell_term_reconstruction(&indicators, &cell_terms, terms)?;

    Ok(ElasticityDwrEstimate {
        eta_signed,
        eta_abs,
        indicators,
        cell_terms,
        face_indicators,
        terms,
        ghost_method: ElasticityGhostMethod::CoarseConsistentEnergy,
        j_primal: coarse.compliance(),
        j_enriched: fine.compliance(),
        dofs: coarse.dof_count(),
        enriched_dofs: fine.dof_count(),
    })
}

fn coarse_parent(
    coarse_grid: &Quadtree,
    fine_grid: &Quadtree,
    fine_cell: CellKey,
    coarse_active: &BTreeSet<CellKey>,
) -> Result<CellKey, CutFemError> {
    let (lo, hi) = fine_grid.rect(fine_cell);
    let center = [f64::midpoint(lo[0], hi[0]), f64::midpoint(lo[1], hi[1])];
    let parent = coarse_grid
        .find_leaf_at(center[0], center[1])
        .ok_or_else(|| {
            invalid(format!(
                "enriched active cell {fine_cell:?} has no containing coarse leaf"
            ))
        })?;
    if !coarse_active.contains(&parent) {
        return Err(invalid(format!(
            "enriched active cell {fine_cell:?} maps to inactive coarse leaf {parent:?}; non-nested active-domain variational crime is unsupported"
        )));
    }
    Ok(parent)
}

fn outer_traction_residual(
    problem: &CutElasticity<'_>,
    coarse: &CutElasticitySolution,
    fine: &CutElasticitySolution,
    fine_grid: &Quadtree,
    fine_cell: CellKey,
    parent: CellKey,
    boundary_traction: Option<BoundaryTraction<'_>>,
) -> Result<f64, CutFemError> {
    let Some(boundary_traction) = boundary_traction else {
        return Ok(0.0);
    };
    let (level, i, j) = fine_cell;
    let nmax = 1u32 << level;
    let corners = fine_grid.corner_nodes(fine_cell);
    let edges = [
        (j == 0, DesignBoxEdge::Bottom, [0usize, 1usize]),
        (i + 1 == nmax, DesignBoxEdge::Right, [1, 2]),
        (j + 1 == nmax, DesignBoxEdge::Top, [2, 3]),
        (i == 0, DesignBoxEdge::Left, [3, 0]),
    ];
    let mut eta = 0.0;
    for (on_boundary, edge, corner_indices) in edges {
        if !on_boundary {
            continue;
        }
        let point_a = fine_grid.node_pos(corners[corner_indices[0]]);
        let point_b = fine_grid.node_pos(corners[corner_indices[1]]);
        let dx = point_b[0] - point_a[0];
        let dy = point_b[1] - point_a[1];
        let Some(segment) = supported_traction_segment(boundary_traction, edge, point_a, point_b)
        else {
            continue;
        };
        let segment_lo = [
            segment.start[0].min(segment.end[0]),
            segment.start[1].min(segment.end[1]),
        ];
        let segment_hi = [
            segment.start[0].max(segment.end[0]),
            segment.start[1].max(segment.end[1]),
        ];
        let enclosure = problem.sdf.enclose(segment_lo, segment_hi);
        if !(enclosure.lo().is_finite() && enclosure.hi().is_finite()) {
            return Err(invalid(format!(
                "non-finite SDF enclosure on DWR supported {edge:?} segment {:?}--{:?}",
                segment.start, segment.end
            )));
        }
        if enclosure.lo() <= 0.0 && enclosure.hi() >= 0.0 {
            return Err(invalid(format!(
                "DWR supported {edge:?} boundary traction segment {:?}--{:?} is cut by the SDF; certified edge clipping is unavailable",
                segment.start, segment.end
            )));
        }
        if enclosure.lo() > 0.0 {
            continue;
        }
        let length = dx.hypot(dy);
        let midpoint = f64::midpoint(segment.local_start, segment.local_end);
        let half_span = 0.5 * (segment.local_end - segment.local_start);
        let gauss = half_span / 3.0f64.sqrt();
        for coordinate in [midpoint - gauss, midpoint + gauss] {
            let point = [point_a[0] + coordinate * dx, point_a[1] + coordinate * dy];
            let (fine_value, _) = fine.value_gradient(fine_grid, fine_cell, point)?;
            let (coarse_value, _) = coarse.value_gradient(problem.grid, parent, point)?;
            let weight_value = sub2(fine_value, coarse_value);
            let traction_value = (segment.value)(point[0], point[1]);
            require_vector_finite(
                traction_value,
                "outer traction in elasticity DWR boundary term",
            )?;
            let value = half_span * length * dot2(traction_value, weight_value);
            require_finite(value, "elasticity DWR outer-traction contribution")?;
            eta += value;
        }
    }
    require_finite(eta, "elasticity DWR cell outer-traction sum")?;
    Ok(eta)
}

fn coarse_ghost_consistent_energy(
    problem: &CutElasticity<'_>,
    coarse: &CutElasticitySolution,
    patch: SharedFacePatch,
    modulus: f64,
) -> Result<f64, CutFemError> {
    let (cell_a, cell_b) = patch.oriented_cells();
    let axis = patch.axis().index();
    let coordinate = patch.coordinate();
    let (tangent_lo, tangent_hi) = patch.tangent_interval();
    if !tangent_lo.is_finite() || !tangent_hi.is_finite() || tangent_hi <= tangent_lo {
        return Err(invalid(format!(
            "coarse ghost patch {patch:?} has non-finite or non-positive length"
        )));
    }
    let normal = patch.axis().normal();
    let gauss = 0.5 / 3.0f64.sqrt();
    let half_length = 0.5 * (tangent_hi - tangent_lo);
    let h = patch.h_f();
    let scale = problem.ghost_gamma * modulus * h * half_length;
    require_finite(scale, "elasticity DWR coarse ghost scale")?;
    let mut quadrature_sum = 0.0;
    for local in [0.5 - gauss, 0.5 + gauss] {
        let tangent = tangent_lo + local * (tangent_hi - tangent_lo);
        let point = if axis == 0 {
            [coordinate, tangent]
        } else {
            [tangent, coordinate]
        };
        let (_, coarse_gradient_a) = coarse.value_gradient(problem.grid, cell_a, point)?;
        let (_, coarse_gradient_b) = coarse.value_gradient(problem.grid, cell_b, point)?;
        let jump_u = normal_derivative_jump(coarse_gradient_a, coarse_gradient_b, normal);
        quadrature_sum += dot2(jump_u, jump_u);
    }
    let value = scale * quadrature_sum;
    require_finite(value, "elasticity DWR coarse ghost-energy contribution")?;
    Ok(value)
}

fn stress(lambda: f64, mu: f64, gradient: [[f64; 2]; 2]) -> [f64; 3] {
    let strain_xx = gradient[0][0];
    let strain_yy = gradient[1][1];
    let engineering_shear = gradient[0][1] + gradient[1][0];
    [
        (lambda + 2.0 * mu) * strain_xx + lambda * strain_yy,
        lambda * strain_xx + (lambda + 2.0 * mu) * strain_yy,
        mu * engineering_shear,
    ]
}

fn stress_strain_contract(stress: [f64; 3], gradient: [[f64; 2]; 2]) -> f64 {
    stress[0] * gradient[0][0]
        + stress[1] * gradient[1][1]
        + stress[2] * (gradient[0][1] + gradient[1][0])
}

fn traction(stress: [f64; 3], normal: [f64; 2]) -> [f64; 2] {
    [
        stress[0] * normal[0] + stress[2] * normal[1],
        stress[2] * normal[0] + stress[1] * normal[1],
    ]
}

fn normal_derivative_jump(
    gradient_a: [[f64; 2]; 2],
    gradient_b: [[f64; 2]; 2],
    normal: [f64; 2],
) -> [f64; 2] {
    [
        dot2(gradient_a[0], normal) - dot2(gradient_b[0], normal),
        dot2(gradient_a[1], normal) - dot2(gradient_b[1], normal),
    ]
}

fn sub2(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[0], a[1] - b[1]]
}

fn sub_gradient(a: [[f64; 2]; 2], b: [[f64; 2]; 2]) -> [[f64; 2]; 2] {
    [sub2(a[0], b[0]), sub2(a[1], b[1])]
}

fn dot2(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}

fn require_vector_finite(value: [f64; 2], what: &str) -> Result<(), CutFemError> {
    if value.iter().all(|component| component.is_finite()) {
        Ok(())
    } else {
        Err(invalid(format!("{what} is non-finite: {value:?}")))
    }
}

fn add_indicator(
    indicators: &mut BTreeMap<CellKey, f64>,
    cell: CellKey,
    contribution: f64,
) -> Result<(), CutFemError> {
    let indicator = indicators.get_mut(&cell).ok_or_else(|| {
        invalid(format!(
            "elasticity DWR allocation targets unknown coarse cell {cell:?}"
        ))
    })?;
    *indicator += contribution;
    require_finite(
        *indicator,
        &format!("elasticity DWR indicator for cell {cell:?}"),
    )
}

fn add_cell_term(
    cell_terms: &mut BTreeMap<CellKey, ElasticityResidualTerms>,
    cell: CellKey,
    add: impl FnOnce(&mut ElasticityResidualTerms),
) -> Result<(), CutFemError> {
    let terms = cell_terms.get_mut(&cell).ok_or_else(|| {
        invalid(format!(
            "elasticity DWR term allocation targets unknown coarse cell {cell:?}"
        ))
    })?;
    add(terms);
    for (name, value) in [
        ("bulk", terms.bulk),
        ("Nitsche", terms.nitsche),
        ("outer traction", terms.outer_traction),
        ("ghost", terms.ghost),
    ] {
        require_finite(
            value,
            &format!("elasticity DWR {name} allocation for cell {cell:?}"),
        )?;
    }
    Ok(())
}

fn require_finite(value: f64, what: &str) -> Result<(), CutFemError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(invalid(format!("{what} is non-finite: {value}")))
    }
}

fn require_reconstruction(
    what: &str,
    cell_sum: f64,
    term_sum: f64,
    accumulation_magnitude: f64,
) -> Result<(), CutFemError> {
    if !(cell_sum.is_finite() && term_sum.is_finite() && accumulation_magnitude.is_finite()) {
        return Err(invalid(format!(
            "{what} reconstruction inputs must be finite: cells={cell_sum}, terms={term_sum}, magnitude={accumulation_magnitude}"
        )));
    }
    let scale = cell_sum
        .abs()
        .max(term_sum.abs())
        .max(accumulation_magnitude.abs())
        .max(1.0);
    if (cell_sum - term_sum).abs() <= 64.0 * f64::EPSILON * scale {
        Ok(())
    } else {
        Err(invalid(format!(
            "{what} does not reconstruct the signed residual: cells={cell_sum:.17e}, terms={term_sum:.17e}"
        )))
    }
}

fn validate_cell_term_reconstruction(
    indicators: &BTreeMap<CellKey, f64>,
    cell_terms: &BTreeMap<CellKey, ElasticityResidualTerms>,
    global: ElasticityResidualTerms,
) -> Result<(), CutFemError> {
    if indicators.len() != cell_terms.len() {
        return Err(invalid(format!(
            "elasticity DWR cell evidence key counts disagree: indicators={}, term decompositions={}",
            indicators.len(),
            cell_terms.len()
        )));
    }

    let mut sums = ElasticityResidualTerms::default();
    let mut magnitudes = ElasticityResidualTerms::default();
    for (&cell, &cell_term) in cell_terms {
        let &indicator = indicators.get(&cell).ok_or_else(|| {
            invalid(format!(
                "elasticity DWR term decomposition targets unknown indicator cell {cell:?}"
            ))
        })?;
        require_reconstruction(
            &format!("elasticity DWR term decomposition for cell {cell:?}"),
            indicator,
            cell_term.total(),
            cell_term.bulk.abs()
                + cell_term.nitsche.abs()
                + cell_term.outer_traction.abs()
                + cell_term.ghost.abs(),
        )?;
        sums.bulk += cell_term.bulk;
        sums.nitsche += cell_term.nitsche;
        sums.outer_traction += cell_term.outer_traction;
        sums.ghost += cell_term.ghost;
        magnitudes.bulk += cell_term.bulk.abs();
        magnitudes.nitsche += cell_term.nitsche.abs();
        magnitudes.outer_traction += cell_term.outer_traction.abs();
        magnitudes.ghost += cell_term.ghost.abs();
    }

    for (name, cell_sum, term_sum, magnitude) in [
        ("bulk", sums.bulk, global.bulk, magnitudes.bulk),
        ("Nitsche", sums.nitsche, global.nitsche, magnitudes.nitsche),
        (
            "outer traction",
            sums.outer_traction,
            global.outer_traction,
            magnitudes.outer_traction,
        ),
        ("ghost", sums.ghost, global.ghost, magnitudes.ghost),
    ] {
        require_reconstruction(
            &format!("elasticity DWR per-cell {name} decomposition"),
            cell_sum,
            term_sum,
            magnitude,
        )?;
    }
    Ok(())
}

fn invalid(what: String) -> CutFemError {
    CutFemError::InvalidElasticityInput { what }
}

#[cfg(test)]
mod tests {
    use super::{
        ElasticityResidualTerms, require_reconstruction, validate_cell_term_reconstruction,
    };
    use std::collections::BTreeMap;

    #[test]
    fn reconstruction_refuses_nonfinite_sums() {
        assert!(require_reconstruction("test", f64::INFINITY, f64::INFINITY, 1.0).is_err());
        assert!(require_reconstruction("test", 1.0, f64::INFINITY, 1.0).is_err());
        assert!(require_reconstruction("test", f64::NAN, 1.0, 1.0).is_err());
        assert!(require_reconstruction("test", 1.0, 1.0, f64::INFINITY).is_err());
        assert!(require_reconstruction("test", 1.0, 1.0, 1.0).is_ok());
    }

    #[test]
    fn cell_term_reconstruction_refuses_inconsistent_evidence() {
        let cell = (0, 0, 0);
        let indicators = BTreeMap::from([(cell, 1.0)]);
        assert!(
            validate_cell_term_reconstruction(
                &indicators,
                &BTreeMap::new(),
                ElasticityResidualTerms::default(),
            )
            .is_err()
        );

        let cell_terms = BTreeMap::from([(
            cell,
            ElasticityResidualTerms {
                bulk: 0.5,
                ..ElasticityResidualTerms::default()
            },
        )]);
        assert!(
            validate_cell_term_reconstruction(
                &indicators,
                &cell_terms,
                ElasticityResidualTerms {
                    bulk: 0.5,
                    ..ElasticityResidualTerms::default()
                },
            )
            .is_err()
        );

        let consistent_terms = BTreeMap::from([(
            cell,
            ElasticityResidualTerms {
                bulk: 1.0,
                ..ElasticityResidualTerms::default()
            },
        )]);
        assert!(
            validate_cell_term_reconstruction(
                &indicators,
                &consistent_terms,
                ElasticityResidualTerms {
                    bulk: 2.0,
                    ..ElasticityResidualTerms::default()
                },
            )
            .is_err()
        );
        assert!(
            validate_cell_term_reconstruction(
                &indicators,
                &consistent_terms,
                ElasticityResidualTerms {
                    bulk: 1.0,
                    ..ElasticityResidualTerms::default()
                },
            )
            .is_ok()
        );
    }
}
