//! Three-dimensional Q1 elasticity on a fixed Cartesian background cut by an SDF.
//!
//! Integrates the retained reference elastic tensor over the `quad3` bulk rules.
//! Native material fields retain full Mandel stiffness, density and identity. The embedded surface is naturally traction free. Zero
//! displacement clamps act on selected background-box nodes; `dirichlet` adds
//! weak displacement conditions on the actual interface. Shared faces carry
//! a positive first-normal-derivative ghost penalty, independently of cut volume.
//!
//! Cell stiffness scales can change without repeating geometric integration.
//! Ghost coefficients average the two incident material moduli times their own
//! scales; the exact discrete scale pullback includes both bulk and ghost work.
//! Geometry derivatives, incompressible/mixed formulations,
//! continuum error bounds and mesh-convergence certification are not claimed.

use std::collections::{BTreeMap, BTreeSet};
use std::ops::ControlFlow;

use crate::fem::recomputed_euclidean_residual_claim;
use crate::quad3::{CutRules3, QuadratureControl3, QuadratureError3, cut_cell_rules3};
use crate::{CutSdf3, HexCell};
use fs_ivl::Interval;
use fs_material::{IsotropicElastic, elastic::acoustic::ElasticAcousticState};

mod material;
use fs_solver::krylov::{CgState, ResidualClaim};
use fs_solver::op::LinearOp;
use fs_sparse::precond::Precond;
pub use material::ElasticMaterialField3;
use material::{MaterialProvider3, ReferenceMaterial3, engineering_strain};

/// Assembly bounds and the dimensionless ghost coefficient.
#[derive(Debug, Clone, Copy)]
pub struct ElasticityOptions3 {
    /// Maximum background hexahedra, before any quadrature allocation.
    pub max_cells: usize,
    /// Maximum active displacement degrees of freedom, including clamps.
    pub max_dofs: usize,
    /// Coefficient multiplying the material modulus and normal-derivative jump.
    /// Legacy isotropy uses `mu`; tensor fields use a conservative Mandel row norm.
    pub ghost_gamma: f64,
}
impl Default for ElasticityOptions3 {
    fn default() -> Self {
        Self {
            max_cells: 32_768,
            max_dofs: 250_000,
            ghost_gamma: 0.1,
        }
    }
}

/// Typed refusal; interrupted assembly and solves publish no partial operator/field.
#[derive(Debug, Clone, PartialEq)]
pub enum ElasticityError3 {
    /// Invalid model, shape, or unrepresentable arithmetic.
    Invalid(&'static str),
    /// Geometric integration refused or consumed its shared budget.
    Quadrature(QuadratureError3),
    /// No positive-volume active cell was integrated.
    EmptyDomain,
    /// A possible cut domain had no quadrature support; refine instead of dropping it.
    UnresolvedSupport(HexCell),
    /// Caller stopped a load or solve operation.
    Cancelled,
    /// The explicit Euclidean residual did not pass within the shared solve cap.
    NotConverged {
        iterations: usize,
        relative_residual: f64,
    },
}
impl From<QuadratureError3> for ElasticityError3 {
    fn from(value: QuadratureError3) -> Self {
        Self::Quadrature(value)
    }
}
impl std::fmt::Display for ElasticityError3 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "3-D CutFEM elasticity refused: {self:?}")
    }
}
impl std::error::Error for ElasticityError3 {}

type Key3 = [usize; 3];
struct Cell3 {
    key: Key3,
    bounds: HexCell,
    nodes: [usize; 8],
    stiffness: [[f64; 24]; 24],
    material: ReferenceMaterial3,
    rules: CutRules3,
    cut: bool,
}
struct GhostPoint3 {
    cells: [usize; 2],
    nodes: Vec<usize>,
    jump: Vec<f64>,
    // Full reference coefficients on each side; interpolation halves each one.
    weights: [f64; 2],
}
impl GhostPoint3 {
    fn scaled_weight(&self, scales: &[f64]) -> f64 {
        let [left, right] = self.cells;
        if self.weights[0] == self.weights[1] {
            // Preserve the established uniform-material arithmetic.
            self.weights[0] * 0.5 * (scales[left] + scales[right])
        } else {
            0.5 * self.weights[0] * scales[left] + 0.5 * self.weights[1] * scales[right]
        }
    }
    fn scale_weights(&self) -> [f64; 2] {
        self.weights.map(|w| 0.5 * w)
    }
}

/// Geometry-integrated, matrix-free elasticity with independent per-cell scales.
pub struct CutElasticity3 {
    nodes: Vec<[f64; 3]>,
    fixed: Vec<bool>,
    // Full-material inventory, independent of stiffness interpolation or supports.
    reference_cell_masses_kg: Option<Vec<f64>>,
    cells: Vec<Cell3>,
    ghosts: Vec<GhostPoint3>,
    scales: Vec<f64>,
    volume_bounds: Interval,
    embedded: Option<dirichlet::DirichletData3>,
}

/// A solved discrete field admitted only by a recomputed `b-Au` residual.
#[derive(Debug, Clone)]
pub struct ElasticitySolution3 {
    coefficients: Vec<f64>,
    compliance: f64,
    iterations: usize,
    residual: ResidualClaim,
}
impl ElasticitySolution3 {
    /// Three displacement entries per active node, matching `nodes()`.
    #[must_use]
    pub fn coefficients(&self) -> &[f64] {
        &self.coefficients
    }
    /// Masked independent-load work, `b^T u`; not continuum-certified compliance.
    #[must_use]
    pub const fn compliance(&self) -> f64 {
        self.compliance
    }
    /// Total CG iterations including residual corrections.
    #[must_use]
    pub const fn iterations(&self) -> usize {
        self.iterations
    }
    /// Explicit true-Euclidean residual provenance.
    #[must_use]
    pub const fn residual_claim(&self) -> ResidualClaim {
        self.residual
    }
}

fn node_key(cell: Key3, corner: usize) -> Key3 {
    std::array::from_fn(|a| cell[a] + ((corner >> a) & 1))
}
fn position(key: Key3, counts: Key3, domain: HexCell) -> [f64; 3] {
    let (lo, hi) = (domain.lo(), domain.hi());
    std::array::from_fn(|a| {
        if key[a] == 0 {
            lo[a]
        } else if key[a] == counts[a] {
            hi[a]
        } else {
            lo[a] + (hi[a] - lo[a]) * (key[a] as f64 / counts[a] as f64)
        }
    })
}

/// Tensor-product values/physical gradients; corner bits select x/y/z endpoints.
fn q1(cell: HexCell, p: [f64; 3]) -> ([f64; 8], [[f64; 3]; 8]) {
    let (lo, hi) = (cell.lo(), cell.hi());
    let h: [f64; 3] = std::array::from_fn(|a| hi[a] - lo[a]);
    let t: [f64; 3] = std::array::from_fn(|a| (p[a] - lo[a]) / h[a]);
    let mut values = [0.0; 8];
    let mut gradients = [[0.0; 3]; 8];
    for corner in 0..8 {
        let factor: [f64; 3] = std::array::from_fn(|a| {
            if corner & (1 << a) == 0 {
                1.0 - t[a]
            } else {
                t[a]
            }
        });
        values[corner] = factor[0] * factor[1] * factor[2];
        for a in 0..3 {
            gradients[corner][a] = (if corner & (1 << a) == 0 { -1.0 } else { 1.0 }) / h[a]
                * factor[(a + 1) % 3]
                * factor[(a + 2) % 3];
        }
    }
    (values, gradients)
}

impl CutElasticity3 {
    /// Build one geometry-integrated operator with all cell scales initially one.
    /// The supplied field must be pure and have sound box/derivative enclosures.
    /// A clamp selects zero displacement only on background-box boundary nodes.
    /// The caller must support each connected component; a small algebraic
    /// residual alone does not prove uniqueness, coercivity or material validity.
    #[allow(clippy::too_many_arguments)]
    pub fn build(
        domain: HexCell,
        counts: [usize; 3],
        sdf: &dyn CutSdf3,
        material: &IsotropicElastic,
        clamp: &dyn Fn([f64; 3]) -> bool,
        options: ElasticityOptions3,
        control: &mut QuadratureControl3<'_>,
    ) -> Result<Self, ElasticityError3> {
        Self::build_core(domain, counts, sdf, material, clamp, options, true, control)
    }

    /// Build from admitted native tensors and densities in the global frame.
    /// Each active cut cell retains its original material identity and strain
    /// domain. Material assignment is fixed while stiffness scales change.
    /// This is Q1 displacement elasticity; admission of an SPD tensor does not
    /// certify freedom from locking or a sufficiently large boundary penalty.
    #[allow(clippy::too_many_arguments)]
    pub fn build_with_material_field(
        domain: HexCell,
        counts: [usize; 3],
        sdf: &dyn CutSdf3,
        materials: ElasticMaterialField3<'_>,
        clamp: &dyn Fn([f64; 3]) -> bool,
        options: ElasticityOptions3,
        control: &mut QuadratureControl3<'_>,
    ) -> Result<Self, ElasticityError3> {
        Self::build_material_core(
            domain,
            counts,
            sdf,
            MaterialProvider3::from_field(materials)?,
            clamp,
            options,
            true,
            control,
        )
    }

    // Only the embedded-boundary builder may defer support admission. It must
    // attach a nonempty boundary rule before publishing the resulting operator.
    #[allow(clippy::too_many_arguments)]
    fn build_core(
        domain: HexCell,
        counts: [usize; 3],
        sdf: &dyn CutSdf3,
        material: &IsotropicElastic,
        clamp: &dyn Fn([f64; 3]) -> bool,
        options: ElasticityOptions3,
        require_box_support: bool,
        control: &mut QuadratureControl3<'_>,
    ) -> Result<Self, ElasticityError3> {
        Self::build_material_core(
            domain,
            counts,
            sdf,
            MaterialProvider3::legacy(material)?,
            clamp,
            options,
            require_box_support,
            control,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn build_material_core(
        domain: HexCell,
        counts: [usize; 3],
        sdf: &dyn CutSdf3,
        materials: MaterialProvider3<'_>,
        clamp: &dyn Fn([f64; 3]) -> bool,
        options: ElasticityOptions3,
        require_box_support: bool,
        control: &mut QuadratureControl3<'_>,
    ) -> Result<Self, ElasticityError3> {
        control.poll()?;
        let count = counts
            .iter()
            .try_fold(1usize, |n, &v| n.checked_mul(v))
            .ok_or(ElasticityError3::Invalid("background size overflow"))?;
        if counts.contains(&0)
            || counts.iter().any(|&n| n > u32::MAX as usize)
            || count > options.max_cells
            || !options.ghost_gamma.is_finite()
            || options.ghost_gamma < 0.0
        {
            return Err(ElasticityError3::Invalid(
                "invalid background or assembly allowance",
            ));
        }
        let mut cells = Vec::new();
        let mut keys = BTreeSet::new();
        let mut volume_bounds = Interval::new(0.0, 0.0);
        for z in 0..counts[2] {
            for y in 0..counts[1] {
                for x in 0..counts[0] {
                    control.poll()?;
                    let key = [x, y, z];
                    let bounds = HexCell::try_new(
                        position(key, counts, domain),
                        position(node_key(key, 7), counts, domain),
                    )
                    .map_err(|_| ElasticityError3::Invalid("collapsed background cell"))?;
                    let rules = cut_cell_rules3(sdf, bounds, control)?;
                    volume_bounds = volume_bounds + rules.volume_bounds();
                    if rules.bulk().is_empty() {
                        if rules.cut_boxes() > 0 {
                            return Err(ElasticityError3::UnresolvedSupport(bounds));
                        }
                        continue;
                    }
                    // Separate assembly classification; the quadrature work counter
                    // counts quad3's own producer calls, not this extra box query.
                    let enclosure = sdf.enclose(bounds.lo(), bounds.hi());
                    control.poll()?;
                    if !enclosure.lo().is_finite()
                        || !enclosure.hi().is_finite()
                        || enclosure.lo() > enclosure.hi()
                    {
                        return Err(ElasticityError3::Invalid("invalid assembly enclosure"));
                    }
                    control.poll()?;
                    let material = materials.at(bounds)?;
                    control.poll()?;
                    let stiffness = assemble_cell_stiffness(bounds, &rules, &material, control)?;
                    for a in 0..8 {
                        keys.insert(node_key(key, a));
                        if keys.len() > options.max_dofs / 3 {
                            return Err(ElasticityError3::Invalid(
                                "active dof allowance exhausted",
                            ));
                        }
                    }
                    cells.push(Cell3 {
                        key,
                        bounds,
                        nodes: [0; 8],
                        stiffness,
                        material,
                        rules,
                        cut: enclosure.hi() >= 0.0,
                    });
                }
            }
        }
        if cells.is_empty() {
            return Err(ElasticityError3::EmptyDomain);
        }
        if !volume_bounds.hi().is_finite() {
            return Err(ElasticityError3::Invalid("domain volume overflow"));
        }
        volume_bounds = Interval::new(volume_bounds.lo().max(0.0), volume_bounds.hi());
        let ids: BTreeMap<_, _> = keys
            .iter()
            .enumerate()
            .map(|(id, &key)| (key, id))
            .collect();
        let nodes: Vec<_> = keys
            .iter()
            .map(|&key| position(key, counts, domain))
            .collect();
        let mut fixed = Vec::with_capacity(nodes.len());
        for (&key, &p) in keys.iter().zip(&nodes) {
            control.poll()?;
            let selected = clamp(p);
            control.poll()?;
            if selected && !(0..3).any(|a| key[a] == 0 || key[a] == counts[a]) {
                return Err(ElasticityError3::Invalid(
                    "clamp must select background-box boundary nodes",
                ));
            }
            fixed.push(selected);
        }
        if require_box_support && !fixed.iter().any(|b| *b) {
            return Err(ElasticityError3::Invalid(
                "no displacement support selected",
            ));
        }
        for cell in &mut cells {
            cell.nodes = std::array::from_fn(|a| ids[&node_key(cell.key, a)]);
        }
        let by_key: BTreeMap<_, _> = cells
            .iter()
            .enumerate()
            .map(|(id, c)| (c.key, id))
            .collect();
        let mut ghosts = Vec::new();
        if options.ghost_gamma > 0.0 {
            for (left, cell) in cells.iter().enumerate() {
                for axis in 0..3 {
                    let mut next = cell.key;
                    next[axis] += 1;
                    if let Some(&right) = by_key.get(&next) {
                        if cell.cut || cells[right].cut {
                            add_ghosts(
                                &cells,
                                left,
                                right,
                                axis,
                                options.ghost_gamma,
                                control,
                                &mut ghosts,
                            )?;
                        }
                    }
                }
            }
        }
        control.poll()?;
        let scales = vec![1.0; cells.len()];
        let reference_cell_masses_kg = collect_reference_masses(&cells)?;
        Ok(Self {
            nodes,
            fixed,
            reference_cell_masses_kg,
            cells,
            ghosts,
            scales,
            volume_bounds,
            embedded: None,
        })
    }

    /// Active node positions, in deterministic lattice-key order.
    #[must_use]
    pub fn nodes(&self) -> &[[f64; 3]] {
        &self.nodes
    }
    /// Number of retained active background cells / independent scale variables.
    #[must_use]
    pub fn cells(&self) -> usize {
        self.cells.len()
    }
    /// Fixed background-node mask; all three components are clamped together.
    #[must_use]
    pub fn fixed(&self) -> &[bool] {
        &self.fixed
    }
    /// Active logical grid keys, in the same order as `scales()` and `volumes()`.
    #[must_use]
    pub fn cell_keys(&self) -> Vec<[usize; 3]> {
        self.cells.iter().map(|c| c.key).collect()
    }
    /// Active cell-to-node connectivity, matching `nodes()`.
    #[must_use]
    pub fn cell_nodes(&self) -> Vec<[usize; 8]> {
        self.cells.iter().map(|c| c.nodes).collect()
    }
    /// Numerical cut volumes; these are not certified individual measures.
    #[must_use]
    pub fn volumes(&self) -> Vec<f64> {
        self.cells.iter().map(|c| c.rules.volume()).collect()
    }
    /// Original admitted material of an active cell, in scale/volume order.
    /// Legacy E/nu construction has no sourced density or material identity.
    #[must_use]
    pub fn cell_material_state(&self, cell: usize) -> Option<&ElasticAcousticState> {
        self.cells.get(cell)?.material.source_state()
    }
    /// Full-material cell masses [kg] from density times retained cut quadrature.
    /// Includes supported material; no SIMP stiffness scaling is applied here.
    /// These are numerical measures, not certified mass bounds. Legacy material
    /// construction returns `None` because it has no physical density.
    #[must_use]
    pub fn reference_cell_masses_kg(&self) -> Option<&[f64]> {
        self.reference_cell_masses_kg.as_deref()
    }
    /// Conservative domain-volume enclosure, independent of quadrature weights.
    #[must_use]
    pub const fn volume_bounds(&self) -> Interval {
        self.volume_bounds
    }
    /// Current positive multipliers of the reference material stiffness.
    #[must_use]
    pub fn scales(&self) -> &[f64] {
        &self.scales
    }
    /// Transactionally replace stiffness scales; no re-integration or remeshing.
    /// Scales in (0,1] represent removal relative to the reference material.
    pub fn set_scales(&mut self, scales: &[f64]) -> Result<(), ElasticityError3> {
        if scales.len() != self.cells.len()
            || !scales
                .iter()
                .all(|s| s.is_finite() && *s > 0.0 && *s <= 1.0)
        {
            return Err(ElasticityError3::Invalid(
                "one finite stiffness scale in (0,1] per active cell required",
            ));
        }
        self.scales.copy_from_slice(scales);
        Ok(())
    }

    /// Integrate one independent body-force density over the actual retained cuts.
    /// Force entries on clamped rows are zero (they do not do displacement work).
    pub fn body_load(
        &self,
        force: &dyn Fn([f64; 3]) -> [f64; 3],
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<f64>, ElasticityError3> {
        let mut rhs = vec![0.0; self.n()];
        for cell in &self.cells {
            for &(p, w) in cell.rules.bulk() {
                if checkpoint().is_break() {
                    return Err(ElasticityError3::Cancelled);
                }
                let f = force(p);
                if checkpoint().is_break() {
                    return Err(ElasticityError3::Cancelled);
                }
                if !f.iter().all(|v| v.is_finite()) {
                    return Err(ElasticityError3::Invalid("nonfinite body force"));
                }
                let (values, _) = q1(cell.bounds, p);
                for a in 0..8 {
                    if !self.fixed[cell.nodes[a]] {
                        for c in 0..3 {
                            rhs[3 * cell.nodes[a] + c] += w * values[a] * f[c];
                        }
                    }
                }
            }
        }
        if !rhs.iter().all(|v| v.is_finite()) {
            return Err(ElasticityError3::Invalid("load overflow"));
        }
        if checkpoint().is_break() {
            return Err(ElasticityError3::Cancelled);
        }
        Ok(rhs)
    }

    /// `u^T (dK/dscale_c) u` for each cell, including its own material coefficient
    /// in adjacent ghost energies and any embedded Nitsche terms. For a matching FIXED load
    /// the compliance derivative is the negative of this vector. Nonzero
    /// prescribed motion adds a density-dependent RHS: include its load
    /// derivative as well. This is a discrete scale derivative, NOT an
    /// SDF/shape derivative; arbitrary supplied vectors carry no solution claim.
    pub fn scale_quadratic_forms(&self, u: &[f64]) -> Result<Vec<f64>, ElasticityError3> {
        if u.len() != self.n() || !u.iter().all(|v| v.is_finite()) {
            return Err(ElasticityError3::Invalid(
                "invalid energy-contraction vector",
            ));
        }
        let mut energies = vec![0.0; self.cells.len()];
        for (id, cell) in self.cells.iter().enumerate() {
            let x: [f64; 24] = std::array::from_fn(|i| {
                if self.fixed[cell.nodes[i / 3]] {
                    0.0
                } else {
                    u[3 * cell.nodes[i / 3] + i % 3]
                }
            });
            for i in 0..24 {
                for j in 0..24 {
                    energies[id] += x[i] * cell.stiffness[i][j] * x[j];
                }
            }
        }
        for face in &self.ghosts {
            let mut energy = 0.0;
            for c in 0..3 {
                let jump: f64 = face
                    .nodes
                    .iter()
                    .zip(&face.jump)
                    .map(|(&node, &j)| {
                        if self.fixed[node] {
                            0.0
                        } else {
                            j * u[3 * node + c]
                        }
                    })
                    .sum();
                energy += jump * jump;
            }
            for (&cell, weight) in face.cells.iter().zip(face.scale_weights()) {
                energies[cell] += weight * energy;
            }
        }
        if !energies.iter().all(|e| e.is_finite()) {
            return Err(ElasticityError3::Invalid("energy contraction overflow"));
        }
        Ok(energies)
    }

    fn diagonal(&self) -> Result<Diagonal3, ElasticityError3> {
        let mut d = vec![0.0; self.n()];
        for (id, cell) in self.cells.iter().enumerate() {
            for i in 0..24 {
                d[3 * cell.nodes[i / 3] + i % 3] += self.scales[id] * cell.stiffness[i][i];
            }
        }
        for face in &self.ghosts {
            let weight = face.scaled_weight(&self.scales);
            for (&node, &j) in face.nodes.iter().zip(&face.jump) {
                for c in 0..3 {
                    d[3 * node + c] += weight * j * j;
                }
            }
        }
        for (i, value) in d.iter_mut().enumerate() {
            if self.fixed[i / 3] {
                *value = 1.0;
            }
            if !value.is_finite() || *value <= 0.0 {
                return Err(ElasticityError3::Invalid("nonpositive/nonfinite diagonal"));
            }
            *value = 1.0 / *value;
            if !value.is_finite() {
                return Err(ElasticityError3::Invalid("diagonal inverse overflow"));
            }
        }
        Ok(Diagonal3(d))
    }

    /// Jacobi-CG with bounded polling and the existing true-Euclidean residual
    /// correction gate. Iteration allowance is cumulative across corrections.
    /// Cancellation is distinct from nonconvergence; neither returns an iterate.
    pub fn solve_controlled(
        &self,
        force: &[f64],
        tolerance: f64,
        max_iterations: usize,
        poll_iterations: usize,
        mut checkpoint: impl FnMut(usize) -> ControlFlow<()>,
    ) -> Result<ElasticitySolution3, ElasticityError3> {
        if force.len() != self.n()
            || !force.iter().all(|f| f.is_finite())
            || !tolerance.is_finite()
            || tolerance <= 0.0
            || tolerance >= 1.0
            || poll_iterations == 0
        {
            return Err(ElasticityError3::Invalid(
                "invalid load, tolerance, or polling interval",
            ));
        }
        if checkpoint(0).is_break() {
            return Err(ElasticityError3::Cancelled);
        }
        let rhs: Vec<f64> = force
            .iter()
            .enumerate()
            .map(|(i, &f)| if self.fixed[i / 3] { 0.0 } else { f })
            .collect();
        let precond = self.diagonal()?;
        let mut x = vec![0.0; self.n()];
        let mut total = 0;
        loop {
            if checkpoint(total).is_break() {
                return Err(ElasticityError3::Cancelled);
            }
            let claim = recomputed_euclidean_residual_claim(self, &x, &rhs);
            let residual = claim.euclidean().expect("explicit Euclidean residual");
            if checkpoint(total).is_break() {
                return Err(ElasticityError3::Cancelled);
            }
            if !residual.is_finite() || !x.iter().all(|v| v.is_finite()) {
                return Err(ElasticityError3::NotConverged {
                    iterations: total,
                    relative_residual: residual,
                });
            }
            if residual < tolerance {
                let compliance: f64 = rhs.iter().zip(&x).map(|(b, u)| b * u).sum();
                if !compliance.is_finite() {
                    return Err(ElasticityError3::Invalid("compliance overflow"));
                }
                if checkpoint(total).is_break() {
                    return Err(ElasticityError3::Cancelled);
                }
                return Ok(ElasticitySolution3 {
                    coefficients: x,
                    compliance,
                    iterations: total,
                    residual: claim,
                });
            }
            if total >= max_iterations {
                return Err(ElasticityError3::NotConverged {
                    iterations: total,
                    relative_residual: residual,
                });
            }
            let mut ax = vec![0.0; self.n()];
            self.apply(&x, &mut ax);
            let r: Vec<f64> = rhs.iter().zip(ax).map(|(b, a)| b - a).collect();
            let mut correction = CgState::new(self, &precond, &r);
            let remaining = max_iterations - total;
            while correction.iters < remaining && !(correction.rel_residual() < tolerance) {
                if checkpoint(total + correction.iters).is_break() {
                    return Err(ElasticityError3::Cancelled);
                }
                let before = correction.iters;
                let _ = correction.run(
                    self,
                    &precond,
                    tolerance,
                    poll_iterations.min(remaining - before),
                );
                correction.history.clear();
                if checkpoint(total + correction.iters).is_break() {
                    return Err(ElasticityError3::Cancelled);
                }
                if correction.iters == before || !correction.rel_residual().is_finite() {
                    break;
                }
            }
            let completed = correction.iters;
            total += completed;
            for (value, delta) in x.iter_mut().zip(correction.x) {
                *value += delta;
            }
            if completed == 0 {
                return Err(ElasticityError3::NotConverged {
                    iterations: total,
                    relative_residual: residual,
                });
            }
        }
    }
}

impl LinearOp for CutElasticity3 {
    fn n(&self) -> usize {
        3 * self.nodes.len()
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        assert_eq!(x.len(), self.n());
        assert_eq!(y.len(), self.n());
        y.fill(0.0);
        for (id, cell) in self.cells.iter().enumerate() {
            let local: [f64; 24] = std::array::from_fn(|j| {
                if self.fixed[cell.nodes[j / 3]] {
                    0.0
                } else {
                    x[3 * cell.nodes[j / 3] + j % 3]
                }
            });
            for i in 0..24 {
                if !self.fixed[cell.nodes[i / 3]] {
                    let value: f64 = cell.stiffness[i]
                        .iter()
                        .zip(local)
                        .map(|(k, u)| k * u)
                        .sum();
                    y[3 * cell.nodes[i / 3] + i % 3] += self.scales[id] * value;
                }
            }
        }
        for face in &self.ghosts {
            let weight = face.scaled_weight(&self.scales);
            for c in 0..3 {
                let jump: f64 = face
                    .nodes
                    .iter()
                    .zip(&face.jump)
                    .map(|(&node, &j)| {
                        if self.fixed[node] {
                            0.0
                        } else {
                            j * x[3 * node + c]
                        }
                    })
                    .sum();
                for (&node, &j) in face.nodes.iter().zip(&face.jump) {
                    if !self.fixed[node] {
                        y[3 * node + c] += weight * j * jump;
                    }
                }
            }
        }
        for (node, &fixed) in self.fixed.iter().enumerate() {
            if fixed {
                y[3 * node..3 * node + 3].copy_from_slice(&x[3 * node..3 * node + 3]);
            }
        }
    }
}
struct Diagonal3(Vec<f64>);
impl Precond for Diagonal3 {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        for ((z, r), d) in z.iter_mut().zip(r).zip(&self.0) {
            *z = r * d;
        }
    }
}

// Both Cartesian and octree assembly use the same constitutive contraction.
// Tensor basis strains and stresses are prepared once per quadrature point,
// avoiding a 6x6 material multiplication for every one of the 300 entries.
fn assemble_cell_stiffness(
    bounds: HexCell,
    rules: &CutRules3,
    material: &ReferenceMaterial3,
    control: &mut QuadratureControl3<'_>,
) -> Result<[[f64; 24]; 24], ElasticityError3> {
    let mut stiffness = [[0.0; 24]; 24];
    for &(p, w) in rules.bulk() {
        control.poll()?;
        let (_, g) = q1(bounds, p);
        if let Some([lambda, mu]) = material.legacy_lame() {
            for i in 0..24 {
                for j in i..24 {
                    let (a, b, ci, cj) = (i / 3, j / 3, i % 3, j % 3);
                    let dot = g[a].iter().zip(g[b]).map(|(u, v)| u * v).sum::<f64>();
                    let value = lambda * g[a][ci] * g[b][cj]
                        + mu * g[a][cj] * g[b][ci]
                        + if ci == cj { mu * dot } else { 0.0 };
                    stiffness[i][j] += w * value;
                }
            }
        } else {
            let strains: [[f64; 6]; 24] =
                std::array::from_fn(|i| engineering_strain(g[i / 3], i % 3));
            let stresses: [[f64; 6]; 24] =
                std::array::from_fn(|i| material.basis_stress(g[i / 3], i % 3));
            for i in 0..24 {
                for j in i..24 {
                    let value: f64 = strains[i].iter().zip(stresses[j]).map(|(e, s)| e * s).sum();
                    stiffness[i][j] += w * value;
                }
            }
        }
    }
    for i in 0..24 {
        if stiffness[i][i] <= 0.0 {
            return Err(ElasticityError3::Invalid(
                "unrepresentable positive cell stiffness",
            ));
        }
        for j in i..24 {
            if !stiffness[i][j].is_finite() {
                return Err(ElasticityError3::Invalid("stiffness overflow"));
            }
            stiffness[j][i] = stiffness[i][j];
        }
    }
    Ok(stiffness)
}

fn collect_reference_masses(cells: &[Cell3]) -> Result<Option<Vec<f64>>, ElasticityError3> {
    if cells
        .iter()
        .all(|cell| cell.material.source_state().is_none())
    {
        return Ok(None);
    }
    cells
        .iter()
        .map(|cell| {
            let material = cell
                .material
                .source_state()
                .ok_or(ElasticityError3::Invalid(
                    "mixed legacy and native material inventory",
                ))?;
            let mass = material.density_kg_m3() * cell.rules.volume();
            if !mass.is_finite() || mass <= 0.0 {
                return Err(ElasticityError3::Invalid(
                    "unrepresentable material cell mass",
                ));
            }
            Ok(mass)
        })
        .collect::<Result<Vec<_>, _>>()
        .map(Some)
}

#[allow(clippy::too_many_arguments)]
fn add_ghosts(
    cells: &[Cell3],
    left: usize,
    right: usize,
    axis: usize,
    gamma: f64,
    control: &mut QuadratureControl3<'_>,
    out: &mut Vec<GhostPoint3>,
) -> Result<(), ElasticityError3> {
    let l = &cells[left];
    let r = &cells[right];
    let (lo, hi) = (l.bounds.lo(), l.bounds.hi());
    let (a, b) = ((axis + 1) % 3, (axis + 2) % 3);
    let h = (hi[axis] - lo[axis]).min(r.bounds.hi()[axis] - r.bounds.lo()[axis]);
    let gauss = [
        (-0.774_596_669_241_483_4, 5.0 / 9.0),
        (0.0, 8.0 / 9.0),
        (0.774_596_669_241_483_4, 5.0 / 9.0),
    ];
    for (u, wu) in gauss {
        for (v, wv) in gauss {
            control.poll()?;
            let mut p = lo;
            p[axis] = hi[axis];
            p[a] = f64::midpoint(lo[a], hi[a]) + 0.5 * (hi[a] - lo[a]) * u;
            p[b] = f64::midpoint(lo[b], hi[b]) + 0.5 * (hi[b] - lo[b]) * v;
            let (_, gl) = q1(l.bounds, p);
            let (_, gr) = q1(r.bounds, p);
            let mut jump = BTreeMap::new();
            for i in 0..8 {
                *jump.entry(l.nodes[i]).or_insert(0.0) += gl[i][axis];
                *jump.entry(r.nodes[i]).or_insert(0.0) -= gr[i][axis];
            }
            let weights = [l, r].map(|cell| {
                gamma
                    * cell.material.ghost_modulus()
                    * h
                    * 0.25
                    * (hi[a] - lo[a])
                    * (hi[b] - lo[b])
                    * wu
                    * wv
            });
            if weights.iter().any(|w| !w.is_finite() || *w <= 0.0)
                || !jump.values().all(|j: &f64| j.is_finite())
            {
                return Err(ElasticityError3::Invalid(
                    "unrepresentable ghost stabilization",
                ));
            }
            let (nodes, jump) = jump.into_iter().unzip();
            out.push(GhostPoint3 {
                cells: [left, right],
                nodes,
                jump,
                weights,
            });
        }
    }
    Ok(())
}

/// Q1-conforming locally refined octree execution of this 3-D kernel.
pub mod adaptive;

/// Oriented reference-surface traction and pressure loads.
pub mod surface;

/// Bulk quadrature stress observations and exact state/material pullbacks.
pub mod stress;

/// Weak prescribed displacement on selected zero-level surface patches.
pub mod dirichlet;
