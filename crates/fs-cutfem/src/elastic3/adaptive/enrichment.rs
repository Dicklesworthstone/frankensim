//! Conforming inter-grid transfer and localized weak residuals for 3-D goals.
//! These are numerical operators, not continuum-error certificates.
use super::*;

fn poll(checkpoint: &mut impl FnMut() -> ControlFlow<()>) -> Result<(), ElasticityError3> {
    if checkpoint().is_break() { Err(ElasticityError3::Cancelled) } else { Ok(()) }
}
fn field(op: &AdaptiveElasticity3, x: &[f64]) -> Result<(), ElasticityError3> {
    if x.len() != op.n() || !x.iter().all(|v| v.is_finite())
        || x.iter().enumerate().any(|(i, v)| op.fixed[i / 3] && *v != 0.0) {
        return Err(ElasticityError3::Invalid("expected finite homogeneous master field"));
    }
    Ok(())
}

/// Sparse Q1 transfer bound to the two actual operators it connects.
/// No field extrapolation, nearest-node substitution or missing-node zero fill.
/// Target active leaves must descend from source active leaves; every source
/// active leaf must retain at least one descendant. Geometry outside active
/// support is not inferred from this map. Both operators must describe the same
/// box, reference material and ghost coefficient. The caller owns the unchanged
/// implicit-domain assumption; no geometry equivalence certificate is minted.
pub struct AdaptiveTransfer3<'a> {
    coarse: &'a AdaptiveElasticity3,
    fine: &'a AdaptiveElasticity3,
    rows: Vec<Vec<(usize, f64)>>,
    parents: Vec<usize>,
}
impl<'a> AdaptiveTransfer3<'a> {
    /// Build bounded sparse interpolation without a dense/all-pairs search.
    /// Nodes are located by at most 21 levels of dyadic ancestor lookups, with
    /// at most eight incident boxes tested at a boundary. `max_terms` bounds
    /// the TOTAL retained scalar coefficients, not the number per row.
    pub fn new(coarse: &'a AdaptiveElasticity3, fine: &'a AdaptiveElasticity3,
        max_terms: usize, mut checkpoint: impl FnMut() -> ControlFlow<()>) -> Result<Self, ElasticityError3> {
        poll(&mut checkpoint)?;
        if coarse.domain != fine.domain || coarse.reference != fine.reference {
            return Err(ElasticityError3::Invalid("enrichment box/material/stabilization mismatch"));
        }
        let lookup: BTreeMap<_, _> = coarse.leaves.iter().enumerate()
            .map(|(i, c)| ((c.level(), c.index()), i)).collect();
        let mut parents = Vec::with_capacity(fine.cells());
        let mut covered = vec![false; coarse.cells()];
        for leaf in &fine.leaves {
            poll(&mut checkpoint)?;
            let parent = (0..=leaf.level()).rev().find_map(|level| {
                lookup.get(&(level, leaf.index().map(|v| v >> (leaf.level() - level)))).copied()
            }).ok_or(ElasticityError3::Invalid("fine active cell has no coarse active ancestor"))?;
            covered[parent] = true;
            parents.push(parent);
        }
        if covered.iter().any(|v| !v) {
            return Err(ElasticityError3::Invalid("enrichment lost a coarse active cell"));
        }
        let mut rows = Vec::with_capacity(fine.nodes.len());
        let mut terms = 0usize;
        for (target, key) in fine.master_lattice.iter().enumerate() {
            poll(&mut checkpoint)?;
            let mut found = None;
            'levels: for level in 0..=20u8 {
                let span = 1u32 << (20 - level);
                let side = 1u32 << level;
                let choices: [Vec<u32>; 3] = std::array::from_fn(|a| {
                    let i = key[a] / span;
                    let mut indices = Vec::with_capacity(2);
                    if key[a] % span == 0 && i > 0 { indices.push(i - 1); }
                    if i < side { indices.push(i); }
                    indices
                });
                for &x in &choices[0] { for &y in &choices[1] { for &z in &choices[2] {
                    if let Some(&cell) = lookup.get(&(level, [x, y, z])) {
                        found = Some(cell); break 'levels;
                    }
                } } }
            }
            let cell = &coarse.raw.cells[found.ok_or(ElasticityError3::Invalid("target master outside coarse active support"))?];
            let p = fine.nodes[target];
            if (0..3).any(|a| p[a] < cell.bounds.lo()[a] || p[a] > cell.bounds.hi()[a]) {
                return Err(ElasticityError3::Invalid("physical and dyadic transfer coordinates disagree"));
            }
            let (values, _) = q1(cell.bounds, p);
            let mut row = BTreeMap::new();
            for (corner, &value) in values.iter().enumerate() {
                if value == 0.0 { continue; }
                for &(master, weight) in &coarse.rows[cell.nodes[corner]] {
                    if coarse.fixed[master] { continue; }
                    *row.entry(master).or_insert(0.0) += value * weight;
                    if row.len() > max_terms.saturating_sub(terms) {
                        return Err(ElasticityError3::Invalid("transfer coefficient budget exhausted"));
                    }
                }
            }
            if row.values().any(|v: &f64| !v.is_finite() || *v < 0.0)
                || (fine.fixed[target] && !row.is_empty()) {
                return Err(ElasticityError3::Invalid("invalid weights or incompatible fine clamp"));
            }
            terms += row.len();
            rows.push(row.into_iter().collect());
        }
        poll(&mut checkpoint)?;
        Ok(Self { coarse, fine, rows, parents })
    }
    /// Source geometry and current constitutive scales (immutably borrowed).
    #[must_use] pub const fn coarse(&self) -> &AdaptiveElasticity3 { self.coarse }
    /// Enriched geometry and current constitutive scales (immutably borrowed).
    #[must_use] pub const fn fine(&self) -> &AdaptiveElasticity3 { self.fine }
    /// Coarse active-cell index for each fine active cell.
    #[must_use] pub fn parents(&self) -> &[usize] { &self.parents }
    /// Physical stiffness field inherited by children. Transfer these scales,
    /// not raw densities through a DIFFERENT filter, when estimating one design.
    #[must_use] pub fn inherited_scales(&self) -> Vec<f64> {
        self.parents.iter().map(|&parent| self.coarse.scales()[parent]).collect()
    }
    /// Interpolate a homogeneous coarse field into fine MASTER coordinates.
    /// Fine hanging values are subsequently reconstructed by the fine operator.
    pub fn prolongate(&self, x: &[f64], mut checkpoint: impl FnMut() -> ControlFlow<()>) -> Result<Vec<f64>, ElasticityError3> {
        poll(&mut checkpoint)?; field(self.coarse, x)?;
        let mut result = vec![0.0; self.fine.n()];
        for (i, row) in self.rows.iter().enumerate() {
            poll(&mut checkpoint)?;
            for c in 0..3 { result[3 * i + c] = row.iter().map(|&(m, w)| w * x[3 * m + c]).sum(); }
        }
        if !result.iter().all(|v| v.is_finite()) { return Err(ElasticityError3::Invalid("transfer overflow")); }
        poll(&mut checkpoint)?; Ok(result)
    }
}

/// Signed contributions to `f(w) - a(u,w)` on one retained numerical cut cell.
/// Ghost faces are apportioned half to each incident cell, exactly once in total.
#[derive(Debug, Clone, Copy, Default)]
pub struct CellResidual3 {
    /// Integrated body-force work against w.
    pub load: f64,
    /// Density-scaled elasticity bilinear form a_bulk(u,w).
    pub bulk: f64,
    /// Density-scaled ghost bilinear form a_ghost(u,w).
    pub ghost: f64,
}
impl CellResidual3 {
    /// Local weak residual, before taking any marking absolute value.
    #[must_use] pub fn residual(self) -> f64 { self.load - self.bulk - self.ghost }
}
impl AdaptiveElasticity3 {
    /// Recompute the true Euclidean residual for an externally retained field.
    /// Used to reject stale/partial fields before assigning them error evidence.
    pub fn field_residual(&self, x: &[f64], rhs: &[f64], mut checkpoint: impl FnMut() -> ControlFlow<()>) -> Result<f64, ElasticityError3> {
        poll(&mut checkpoint)?; field(self, x)?;
        if rhs.len() != self.n() || !rhs.iter().all(|v| v.is_finite()) {
            return Err(ElasticityError3::Invalid("invalid residual load"));
        }
        let b: Vec<_> = rhs.iter().enumerate().map(|(i, &v)| if self.fixed[i / 3] { 0.0 } else { v }).collect();
        let residual = recomputed_euclidean_residual_claim(self, x, &b).euclidean().expect("explicit Euclidean residual");
        poll(&mut checkpoint)?;
        if !residual.is_finite() { return Err(ElasticityError3::Invalid("nonfinite field residual")); }
        Ok(residual)
    }

    /// Localize the ACTUAL cut-operator weak residual using retained cell blocks
    /// and ghost jumps; there is no recovered stress or substituted stiffness.
    /// The body producer must be pure and match the solve's force density.
    /// This works on unconverged u (e.g. a prolonged coarse field) intentionally.
    pub fn cell_residuals(&self, u: &[f64], w: &[f64], body: &dyn Fn([f64; 3]) -> [f64; 3],
        mut checkpoint: impl FnMut() -> ControlFlow<()>) -> Result<Vec<CellResidual3>, ElasticityError3> {
        poll(&mut checkpoint)?; field(self, u)?; field(self, w)?;
        let u = self.physical_displacements(u)?;
        let w = self.physical_displacements(w)?;
        let mut result = vec![CellResidual3::default(); self.cells()];
        for (id, cell) in self.raw.cells.iter().enumerate() {
            poll(&mut checkpoint)?;
            let uc: [f64; 24] = std::array::from_fn(|i| u[3 * cell.nodes[i / 3] + i % 3]);
            let wc: [f64; 24] = std::array::from_fn(|i| w[3 * cell.nodes[i / 3] + i % 3]);
            for (i, &wi) in wc.iter().enumerate() {
                let applied: f64 = cell.stiffness[i].iter().zip(&uc).map(|(a, x)| a * x).sum();
                result[id].bulk += self.scales()[id] * wi * applied;
            }
            for &(p, weight) in cell.rules.bulk() {
                poll(&mut checkpoint)?;
                let f = body(p);
                poll(&mut checkpoint)?;
                if !f.iter().all(|v| v.is_finite()) { return Err(ElasticityError3::Invalid("nonfinite residual body force")); }
                let (n, _) = q1(cell.bounds, p);
                for a in 0..8 { for c in 0..3 { result[id].load += weight * n[a] * f[c] * wc[3 * a + c]; } }
            }
        }
        for face in &self.raw.ghosts {
            poll(&mut checkpoint)?;
            let scale = 0.5 * (self.scales()[face.cells[0]] + self.scales()[face.cells[1]]);
            let mut work = 0.0;
            for c in 0..3 {
                let ju: f64 = face.nodes.iter().zip(&face.jump).map(|(&n, &j)| j * u[3 * n + c]).sum();
                let jw: f64 = face.nodes.iter().zip(&face.jump).map(|(&n, &j)| j * w[3 * n + c]).sum();
                work += scale * face.weight * ju * jw;
            }
            for &cell in &face.cells { result[cell].ghost += 0.5 * work; }
        }
        if result.iter().any(|r| !r.load.is_finite() || !r.bulk.is_finite() || !r.ghost.is_finite() || !r.residual().is_finite()) {
            return Err(ElasticityError3::Invalid("local residual arithmetic overflow"));
        }
        poll(&mut checkpoint)?; Ok(result)
    }
}
