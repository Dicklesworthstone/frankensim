//! Geometry-owned von Karman stretching of prestressed planar films.
//!
//! The supplied plate pencil already contains bending and installed tension.
//! This reduction adds ONLY the incremental membrane energy
//! `min_u 1/2 integral (sym grad u + grad w tensor grad w / 2):C:(...)`.
//! Rim constraints act on in-plane displacement, not on every interior node.
//! Static condensation is exact for the supplied P1 in-plane space: one cold
//! elastic factorization and one RHS per symmetric pair of transverse modes.
//! Runtime retains positive facet strain energies instead of a dense quartic
//! tensor. It allocates nothing and introduces no new time integrator.
//!
//! Scope: moderate slopes, linear plane-stress material, quasi-static in-plane
//! relaxation. In-plane inertia, wrinkling, slackening and film viscoelasticity
//! are NOT inferred. Mesh/mode convergence and material calibration still matter.
use crate::{ModePair, PlateError, PlateMesh, PlateModel, PlateSection};

/// Explicit cold-work/storage ceilings, independent of instrument names.
#[derive(Debug, Clone, Copy)]
pub struct MembraneReductionBudget {
    /// All supplied modes must fit; none are silently dropped.
    pub max_modes: usize,
    /// Maximum geometry nodes, including constrained nodes.
    pub max_nodes: usize,
    /// Maximum facets times symmetric modal pairs.
    pub max_facet_pairs: usize,
    /// Maximum dense in-plane matrix plus all RHS scalar entries.
    pub max_solve_entries: usize,
    /// Relative mass, eigenpair and static-solve residual tolerance.
    pub relative_tolerance: f64,
}

#[derive(Debug, Clone)]
struct Facet {
    a: [f64; 9],
    gradients: Vec<[f64; 2]>,
    strains: Vec<[f64; 3]>,
}

/// Prestressed linear pencil plus positive, statically relaxed stretching.
/// Coordinates have unit modal mass. Dynamic state belongs to the time owner.
#[derive(Debug, Clone)]
pub struct MembraneReduction {
    omegas: Vec<f64>,
    linear: Vec<f64>,
    pairs: Vec<(usize, usize)>,
    facets: Vec<Facet>,
    solve_residual: f64,
}

fn bad(what: &'static str) -> PlateError { PlateError::BadSection { what } }
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn mul(a: &[f64; 9], x: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|i| a[3*i]*x[0] + a[3*i+1]*x[1] + a[3*i+2]*x[2])
}

impl MembraneReduction {
    /// Prepare a film of arbitrary planar triangulated shape and uniform section.
    /// `fixed_nodes` fix BOTH incremental in-plane translations at those nodes;
    /// out-of-plane supports and prestress remain those of `model`.
    ///
    /// The cold solver callback receives `(dimension, row_major_K, RHS_columns)`
    /// and returns solution columns for `K u = rhs`. Use the existing fs-la or
    /// fs-solver owner, factor ONCE, and solve all columns with that factor.
    /// Every returned solution is independently residual-checked here. The
    /// callback is not called when every in-plane DOF is fixed.
    ///
    /// # Errors
    /// Invalid geometry/material/basis, exceeded budgets or invalid static solve.
    pub fn new<F>(
        mesh: &PlateMesh, section: &PlateSection, model: &PlateModel,
        modes: &[ModePair], fixed_nodes: &[usize], budget: MembraneReductionBudget,
        solve: F,
    ) -> Result<Self, PlateError>
    where F: FnOnce(usize, &[f64], &[Vec<f64>]) -> Result<Vec<Vec<f64>>, PlateError> {
        let n = modes.len();
        let nodes = mesh.nodes.len();
        let count = n.checked_add(1).and_then(|v| n.checked_mul(v)).map(|v| v/2)
            .ok_or_else(|| bad("membrane pair count overflow"))?;
        let nn = n.checked_mul(n).ok_or_else(|| bad("membrane projection overflow"))?;
        let dofs = nodes.checked_mul(3).ok_or_else(|| bad("membrane DOF count overflow"))?;
        if n == 0 || n > budget.max_modes || nodes == 0 || nodes > budget.max_nodes
            || mesh.tris.is_empty()
            || count.checked_mul(mesh.tris.len()).is_none_or(|v| v > budget.max_facet_pairs)
            || !budget.relative_tolerance.is_finite() || budget.relative_tolerance <= 0.0
            || budget.relative_tolerance >= 1.0 {
            return Err(bad("membrane reduction exceeds declared budgets or has invalid tolerance"));
        }
        // Reuse the geometry owner's topology/orientation admission.
        let checked = PlateMesh::from_unstructured(mesh.nodes.clone(), mesh.tris.clone())?;
        section.validate()?;
        if model.dof_map.len() != dofs || model.k.nrows() != model.free
            || model.k.ncols() != model.free || model.m.nrows() != model.free
            || model.m.ncols() != model.free
            || model.dof_map.iter().flatten().any(|&i| i >= model.free)
            || modes.iter().any(|m| m.phi.len() != model.free || !m.lambda.is_finite()
                || m.lambda <= 0.0 || m.phi.iter().any(|v| !v.is_finite())) {
            return Err(bad("membrane modes do not match the finite plate pencil"));
        }
        let tol = budget.relative_tolerance;
        let mut linear = vec![0.0; nn];
        let mut mass = vec![0.0; model.free];
        let mut stiffness = mass.clone();
        for j in 0..n {
            model.m.spmv(&modes[j].phi, &mut mass);
            model.k.spmv(&modes[j].phi, &mut stiffness);
            let mut residual = 0.0_f64;
            let mut scale = f64::MIN_POSITIVE;
            for (&k, &m) in stiffness.iter().zip(&mass) {
                let expected = modes[j].lambda * m;
                if !k.is_finite() || !expected.is_finite() {
                    return Err(bad("nonfinite membrane eigenpair residual"));
                }
                residual = residual.max((k - expected).abs());
                scale = scale.max(k.abs()).max(expected.abs());
            }
            if residual > tol * scale { return Err(bad("membrane shape is not a pencil eigenvector")); }
            for i in 0..n {
                let m = modes[i].phi.iter().zip(&mass).map(|(a,b)| a*b).sum::<f64>();
                if !m.is_finite() || (m - if i == j {1.0} else {0.0}).abs() > tol {
                    return Err(bad("membrane basis is not mass orthonormal"));
                }
                let k = modes[i].phi.iter().zip(&stiffness).map(|(a,b)| a*b).sum::<f64>();
                let expected = if i == j { modes[i].lambda } else { 0.0 };
                let scale = modes[i].lambda.sqrt() * modes[j].lambda.sqrt();
                if !k.is_finite() || !scale.is_finite() || (k-expected).abs() > tol*scale {
                    return Err(bad("membrane frequencies disagree with the supplied pencil"));
                }
                linear[i*n+j] = k;
            }
        }
        // Preserve the actual projected quadratic energy with exact symmetry.
        for i in 0..n { for j in 0..i {
            let k = 0.5*linear[i*n+j] + 0.5*linear[j*n+i];
            linear[i*n+j] = k; linear[j*n+i] = k;
        }}
        let mut fixed = vec![false; nodes];
        for &node in fixed_nodes {
            if node >= nodes || fixed[node] { return Err(bad("invalid or duplicate in-plane support node")); }
            fixed[node] = true;
        }
        if fixed_nodes.len() < 2 { return Err(bad("in-plane membrane rigid motions require explicit supports")); }
        let mut free = 0;
        let mut map = vec![[None; 2]; nodes];
        for node in 0..nodes { if !fixed[node] {
            map[node] = [Some(free), Some(free+1)]; free += 2;
        }}
        let entries = free.checked_mul(free).and_then(|v|
            count.checked_mul(free).and_then(|r| v.checked_add(r)))
            .ok_or_else(|| bad("in-plane condensation storage overflow"))?;
        if entries > budget.max_solve_entries { return Err(bad("in-plane condensation exceeds solve budget")); }
        let pairs: Vec<_> = (0..n).flat_map(|i| (i..n).map(move |j| (i,j))).collect();
        let mut elastic = vec![0.0; free*free];
        let mut rhs = vec![vec![0.0; free]; count];
        let mut facets = Vec::with_capacity(checked.tris.len());
        let mut columns = Vec::with_capacity(checked.tris.len());
        for tri in &checked.tris {
            let (x0,y0) = checked.nodes[tri[0]];
            let (x1,y1) = checked.nodes[tri[1]];
            let (x2,y2) = checked.nodes[tri[2]];
            let twice_area = (x1-x0)*(y2-y0) - (x2-x0)*(y1-y0);
            if !twice_area.is_finite() || twice_area <= 0.0 { return Err(bad("membrane requires finite positive facet area")); }
            let grad = [[(y1-y2)/twice_area,(x2-x1)/twice_area],
                [(y2-y0)/twice_area,(x0-x2)/twice_area],
                [(y0-y1)/twice_area,(x1-x0)/twice_area]];
            // D=plane_stress*h^3/12; A here includes element area.
            let a = section.d.map(|d| 0.5*twice_area*12.0/(section.thickness*section.thickness)*d);
            if a.iter().any(|v| !v.is_finite()) { return Err(bad("membrane constitutive scaling overflow")); }
            let b: [[f64;3];6] = core::array::from_fn(|i| {
                let [dx,dy] = grad[i/2];
                if i%2 == 0 { [dx,0.0,dy] } else { [0.0,dy,dx] }
            });
            let ids: [Option<usize>;6] = core::array::from_fn(|i| map[tri[i/2]][i%2]);
            for i in 0..6 { if let Some(row) = ids[i] {
                for j in 0..6 { if let Some(col) = ids[j] {
                    elastic[row*free+col] += dot(b[i],mul(&a,b[j]));
                }}
            }}
            let mut gradients = Vec::with_capacity(n);
            for mode in modes {
                let mut g = [0.0;2];
                for i in 0..3 {
                    let w = model.dof_map[3*tri[i]].map_or(0.0, |d| mode.phi[d]);
                    g[0] += grad[i][0]*w; g[1] += grad[i][1]*w;
                }
                gradients.push(g);
            }
            let mut strains = Vec::with_capacity(count);
            for (pair, &(i,j)) in pairs.iter().enumerate() {
                let [ix,iy] = gradients[i]; let [jx,jy] = gradients[j];
                let eta = if i == j { [0.5*ix*ix,0.5*iy*iy,ix*iy] }
                    else { [ix*jx,iy*jy,ix*jy+jx*iy] };
                let stress = mul(&a,eta);
                for d in 0..6 { if let Some(row) = ids[d] { rhs[pair][row] -= dot(b[d],stress); }}
                strains.push(eta);
            }
            facets.push(Facet { a, gradients, strains });
            columns.push((ids,b));
        }
        if elastic.iter().chain(rhs.iter().flatten()).any(|v| !v.is_finite()) {
            return Err(bad("nonfinite in-plane stiffness or modal pair load"));
        }
        let solutions = if free == 0 { vec![Vec::new();count] } else { solve(free,&elastic,&rhs)? };
        if solutions.len() != count || solutions.iter().any(|u| u.len() != free || u.iter().any(|v| !v.is_finite())) {
            return Err(bad("static solver returned invalid displacement columns"));
        }
        let mut solve_residual = 0.0_f64;
        for (u,load) in solutions.iter().zip(&rhs) {
            let mut residual = 0.0_f64; let mut scale = f64::MIN_POSITIVE;
            for row in 0..free {
                let actual = elastic[row*free..(row+1)*free].iter().zip(u).map(|(k,u)| k*u).sum::<f64>();
                if !actual.is_finite() { return Err(bad("static membrane solve overflow")); }
                residual = residual.max((actual-load[row]).abs());
                scale = scale.max(actual.abs()).max(load[row].abs());
            }
            solve_residual = solve_residual.max(residual/scale);
            if residual > tol*scale { return Err(bad("static membrane relaxation residual exceeds tolerance")); }
        }
        for (facet,(ids,b)) in facets.iter_mut().zip(&columns) {
            for (strain,u) in facet.strains.iter_mut().zip(&solutions) {
                for d in 0..6 { if let Some(row) = ids[d] {
                    for c in 0..3 { strain[c] += b[d][c]*u[row]; }
                }}
                if strain.iter().any(|v| !v.is_finite()) { return Err(bad("condensed membrane strain overflow")); }
            }
        }
        Ok(Self { omegas: modes.iter().map(|m|m.lambda.sqrt()).collect(),
            linear, pairs, facets, solve_residual })
    }

    /// Supplied small-signal angular frequencies [rad/s].
    #[must_use]
    pub fn omegas(&self) -> &[f64] { &self.omegas }
    /// Retained mass-normalized coordinates.
    #[must_use]
    pub fn mode_count(&self) -> usize { self.omegas.len() }
    /// Worst independently checked cold static-solve relative residual.
    #[must_use]
    pub fn solve_residual(&self) -> f64 { self.solve_residual }
    fn strain(&self, facet: &Facet, q: &[f64]) -> [f64;3] {
        let mut strain = [0.0;3];
        for (&(i,j),s) in self.pairs.iter().zip(&facet.strains) {
            for c in 0..3 { strain[c] += q[i]*q[j]*s[c]; }
        }
        strain
    }
    /// Incremental stretching energy [J], nonnegative by construction.
    /// Invalid coordinates return NaN for the time owner's finite-set gate.
    #[must_use]
    pub fn stretching_energy(&self, q: &[f64]) -> f64 {
        if q.len() != self.mode_count() || q.iter().any(|v| !v.is_finite()) { return f64::NAN; }
        self.facets.iter().map(|f| { let s=self.strain(f,q); 0.5*dot(s,mul(&f.a,s)) }).sum()
    }
    /// Total potential [J], with installed tension and bending counted once.
    #[must_use]
    pub fn potential(&self, q: &[f64]) -> f64 {
        let n = self.mode_count();
        if q.len() != n || q.iter().any(|v| !v.is_finite()) { return f64::NAN; }
        let mut energy = self.stretching_energy(q);
        for i in 0..n { for j in 0..n { energy += 0.5*q[i]*self.linear[i*n+j]*q[j]; }}
        energy
    }
    /// Exact derivative of the same positive facet energy. No heap scratch.
    pub fn gradient(&self, q: &[f64], out: &mut [f64]) {
        let n = self.mode_count();
        if q.len() != n || out.len() != n || q.iter().any(|v| !v.is_finite()) { out.fill(f64::NAN); return; }
        for i in 0..n { out[i] = (0..n).map(|j| self.linear[i*n+j]*q[j]).sum(); }
        for f in &self.facets {
            let stress = mul(&f.a,self.strain(f,q));
            for (&(i,j),s) in self.pairs.iter().zip(&f.strains) {
                let force = dot(*s,stress);
                out[i] += q[j]*force; out[j] += q[i]*force;
            }
        }
    }
    /// Maximum P1 transverse slope, for caller-owned validity diagnostics.
    #[must_use]
    pub fn maximum_slope(&self, q: &[f64]) -> f64 {
        if q.len() != self.mode_count() || q.iter().any(|v| !v.is_finite()) { return f64::NAN; }
        let mut maximum = 0.0_f64;
        for f in &self.facets {
            let mut g = [0.0;2];
            for (d,&q) in f.gradients.iter().zip(q) { g[0] += d[0]*q; g[1] += d[1]*q; }
            let slope = g[0].hypot(g[1]);
            if !slope.is_finite() { return f64::NAN; }
            maximum = maximum.max(slope);
        }
        maximum
    }
}
