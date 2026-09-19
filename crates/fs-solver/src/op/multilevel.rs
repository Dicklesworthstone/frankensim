//! Sparse recursive Galerkin correction with a matrix-free finest operator.
//!
//! B_f = D_f^-1 + P_f V_1 P_f^T. Intermediate V-cycles use a forward
//! Gauss-Seidel sweep, a recursive residual correction, and a backward sweep.
//! Only the bounded LAST level is factored densely. For symmetric level
//! operators with positive diagonals and an SPD bottom inverse, each cycle
//! is fixed, linear and SPD: its smoother contribution is
//! (D+U)^-1 D (D+L)^-1, plus a congruence of the next SPD correction.
//! The additive finest diagonal covers coordinates absent from P_f.
//!
//! This is not a rediscretization hierarchy. The caller supplies the first
//! sparse P_f^T A_f P_f (e.g. from element contractions); subsequent operators
//! are formed here by bounded sparse Galerkin products. No fine matrix or fine
//! coordinate probing is needed. Fine SPD/first-Galerkin consistency remain
//! caller obligations; preparation is not a coercivity or convergence proof.
use std::collections::BTreeMap;
use std::ops::ControlFlow;
use fs_la::factor::{Cholesky, cholesky};
use fs_sparse::{Coo, Csr, ops::transpose, precond::Precond};
use super::LinearOp;

/// Structural limits, not a peak allocator/RSS bound. At most 20 levels and
/// 512 bottom coordinates are admitted, independent of caller-supplied caps.
#[derive(Debug, Clone, Copy)]
pub struct MultilevelBudget {
    /// Finest vector size; every subsequent level must strictly shrink.
    pub max_fine_dofs: usize,
    /// Number of levels INCLUDING the matrix-free finest level.
    pub max_levels: usize,
    /// Total retained prolongation entries (transposes have the same count).
    pub max_transfer_entries: usize,
    /// Total entries in all retained sparse coarse matrices.
    pub max_matrix_entries: usize,
    /// Cumulative upper-triangle Galerkin summands, including failed setup.
    pub max_galerkin_products: usize,
    /// Dense bottom factor dimension, NOT the first coarse dimension.
    pub max_coarsest_dofs: usize,
}
impl Default for MultilevelBudget {
    fn default() -> Self {
        Self { max_fine_dofs: 250_000, max_levels: 12,
            max_transfer_entries: 4_000_000, max_matrix_entries: 8_000_000,
            max_galerkin_products: 200_000_000, max_coarsest_dofs: 192 }
    }
}
/// Setup work. A product is one accumulated p_i*a_ij*p_j summand, not a
/// fine apply, scalar flop count, or wall-clock measurement.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct MultilevelWork {
    /// Attempted, admitted summands, including discarded candidates.
    pub galerkin_products: usize,
    /// Entries in completed admitted coarse matrices; partial candidates are
    /// not counted here (their products above are still charged).
    pub matrix_entries: usize,
}
/// No partial hierarchy is returned on any refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MultilevelError {
    Invalid(&'static str),
    Budget(&'static str),
    Cancelled,
    Nonsymmetric,
    NotPositiveDefinite,
}
impl std::fmt::Display for MultilevelError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "multilevel preparation refused: {self:?}")
    }
}
impl std::error::Error for MultilevelError {}

/// Shared by the physics owner's first Galerkin assembly and all deeper
/// products. Polls bracket row/level/factor work and each 1,024 summands.
/// Sparse transpose and the bounded bottom factorization are not preemptible.
pub struct MultilevelControl<'a> {
    budget: MultilevelBudget,
    work: MultilevelWork,
    checkpoint: &'a mut dyn FnMut(MultilevelWork) -> ControlFlow<()>,
}
impl<'a> MultilevelControl<'a> {
    pub fn new(budget: MultilevelBudget,
        checkpoint: &'a mut dyn FnMut(MultilevelWork) -> ControlFlow<()>) -> Result<Self, MultilevelError> {
        if !(2..=20).contains(&budget.max_levels) || !(1..=512).contains(&budget.max_coarsest_dofs) {
            return Err(MultilevelError::Invalid("level/bottom limits outside supported range"));
        }
        let mut result = Self { budget, work: MultilevelWork::default(), checkpoint };
        result.poll()?; Ok(result)
    }
    #[must_use] pub const fn work(&self) -> MultilevelWork { self.work }
    #[must_use] pub const fn budget(&self) -> MultilevelBudget { self.budget }
    pub fn poll(&mut self) -> Result<(), MultilevelError> {
        if (self.checkpoint)(self.work).is_break() { Err(MultilevelError::Cancelled) } else { Ok(()) }
    }
    fn product(&mut self) -> Result<(), MultilevelError> {
        if self.work.galerkin_products >= self.budget.max_galerkin_products {
            return Err(MultilevelError::Budget("Galerkin summands"));
        }
        self.work.galerkin_products += 1;
        if self.work.galerkin_products % 1024 == 0 { self.poll()?; }
        Ok(())
    }
    fn remaining_entries(&self) -> usize { self.budget.max_matrix_entries - self.work.matrix_entries }
    fn retain(&mut self, n: usize) -> Result<(), MultilevelError> {
        if n > self.remaining_entries() { return Err(MultilevelError::Budget("sparse hierarchy entries")); }
        self.work.matrix_entries += n; self.poll()
    }
}

/// Upper-triangle accumulator for a symmetric Galerkin form. Feed the FULL
/// bilinear expansion: entries with row>column are ignored, NOT folded/doubled.
/// Each admitted upper contribution consumes work. The finished matrix mirrors
/// the computed upper triangle exactly; callers must not pass nonsymmetric laws.
/// One in-progress accumulator is intended per control; candidate storage is
/// bounded by the remaining hierarchy-entry allowance.
pub struct SymmetricGalerkin {
    n: usize,
    upper: BTreeMap<(usize, usize), f64>,
    entries: usize,
    entry_cap: usize,
}
impl SymmetricGalerkin {
    pub fn new(n: usize, control: &mut MultilevelControl<'_>) -> Result<Self, MultilevelError> {
        control.poll()?;
        if n == 0 || n > control.budget.max_fine_dofs {
            return Err(MultilevelError::Budget("Galerkin dimension"));
        }
        Ok(Self { n, upper: BTreeMap::new(), entries: 0, entry_cap: control.remaining_entries() })
    }
    pub fn add(&mut self, i: usize, j: usize, value: f64,
        control: &mut MultilevelControl<'_>) -> Result<(), MultilevelError> {
        if i >= self.n || j >= self.n { return Err(MultilevelError::Invalid("Galerkin index")); }
        if i > j { return Ok(()); }
        control.product()?;
        if !value.is_finite() { return Err(MultilevelError::Invalid("nonfinite Galerkin summand")); }
        if value == 0.0 { return Ok(()); }
        if let Some(entry) = self.upper.get_mut(&(i, j)) {
            *entry += value;
            if !entry.is_finite() { return Err(MultilevelError::Invalid("Galerkin accumulation overflow")); }
        } else {
            let additional = if i == j { 1 } else { 2 };
            if additional > self.entry_cap.saturating_sub(self.entries) {
                return Err(MultilevelError::Budget("Galerkin candidate entries"));
            }
            self.entries += additional;
            self.upper.insert((i, j), value);
        }
        Ok(())
    }
    /// Finish one matrix. The hierarchy constructor accounts retained entries;
    /// this method does not double-charge a caller's first coarse matrix.
    pub fn finish(self, control: &mut MultilevelControl<'_>) -> Result<Csr, MultilevelError> {
        control.poll()?;
        if self.entries > control.remaining_entries() { return Err(MultilevelError::Budget("Galerkin candidate entries")); }
        let mut coo = Coo::new(self.n, self.n);
        for (entry, ((i, j), value)) in self.upper.into_iter().enumerate() {
            if entry % 256 == 0 { control.poll()?; }
            coo.push(i, j, value);
            if i != j { coo.push(j, i, value); }
        }
        let result = coo.assemble(); control.poll()?; Ok(result)
    }
}
fn transfer(p: &Csr, control: &mut MultilevelControl<'_>) -> Result<(), MultilevelError> {
    if p.nrows() == 0 || p.ncols() == 0 || p.ncols() >= p.nrows() {
        return Err(MultilevelError::Invalid("transfer must strictly reduce dimension"));
    }
    if p.nnz() > control.budget.max_transfer_entries { return Err(MultilevelError::Budget("transfer entries")); }
    for i in 0..p.nrows() {
        if i % 256 == 0 { control.poll()?; }
        let (c, v) = p.row(i);
        if c.len() != v.len() || c.iter().any(|&j| j >= p.ncols())
            || c.windows(2).any(|w| w[0] >= w[1]) || v.iter().any(|v| !v.is_finite()) {
            return Err(MultilevelError::Invalid("invalid sparse transfer row"));
        }
    }
    Ok(())
}
fn matrix(a: &Csr, control: &mut MultilevelControl<'_>) -> Result<Vec<f64>, MultilevelError> {
    if a.nrows() == 0 || a.nrows() != a.ncols() { return Err(MultilevelError::Invalid("coarse operator must be square")); }
    if a.nnz() > control.remaining_entries() { return Err(MultilevelError::Budget("sparse hierarchy entries")); }
    let mut diagonal = Vec::with_capacity(a.nrows());
    for i in 0..a.nrows() {
        if i % 64 == 0 { control.poll()?; }
        let (c, v) = a.row(i);
        if c.len() != v.len() || c.iter().any(|&j| j >= a.nrows())
            || c.windows(2).any(|w| w[0] >= w[1]) || v.iter().any(|v| !v.is_finite()) {
            return Err(MultilevelError::Invalid("invalid sparse operator row"));
        }
        for (&j, &value) in c.iter().zip(v) {
            // No silent symmetrization of a caller's nonsymmetric operator.
            if a.get(j, i) != value { return Err(MultilevelError::Nonsymmetric); }
        }
        let d = a.get(i, i);
        if !d.is_finite() || d <= 0.0 { return Err(MultilevelError::NotPositiveDefinite); }
        diagonal.push(d);
    }
    Ok(diagonal)
}

/// Form P^T A P by deterministic row traversal. No dense intermediate or
/// unbounded sparse-matrix multiplication is hidden in this operation.
/// A must already be a finite symmetric matrix with valid CSR storage.
pub fn sparse_galerkin(a: &Csr, p: &Csr, control: &mut MultilevelControl<'_>) -> Result<Csr, MultilevelError> {
    control.poll()?;
    if a.nrows() != p.nrows() || a.nrows() != a.ncols() {
        return Err(MultilevelError::Invalid("Galerkin operator/transfer shape"));
    }
    transfer(p, control)?;
    // The caller may use this primitive without constructing a hierarchy.
    // Validate symmetry/finite storage without retaining/counting A twice.
    if a.nrows() > control.budget.max_fine_dofs { return Err(MultilevelError::Budget("Galerkin source dimension")); }
    let mut out = SymmetricGalerkin::new(p.ncols(), control)?;
    for i in 0..a.nrows() {
        control.poll()?;
        let (ac, av) = a.row(i); let (pc, pv) = p.row(i);
        if ac.iter().any(|&j| j >= a.nrows()) || ac.windows(2).any(|w| w[0] >= w[1])
            || av.iter().any(|v| !v.is_finite()) { return Err(MultilevelError::Invalid("invalid Galerkin source row")); }
        for (&j, &value) in ac.iter().zip(av) {
            if a.get(j, i) != value { return Err(MultilevelError::Nonsymmetric); }
            if value == 0.0 { continue; }
            let (qc, qv) = p.row(j);
            for (&k, &u) in pc.iter().zip(pv) { for (&l, &v) in qc.iter().zip(qv) {
                if k <= l { out.add(k, l, (u * value) * v, control)?; }
            } }
        }
    }
    out.finish(control)
}
struct Transfer { p: Csr, pt: Csr }
struct Level { a: Csr, diagonal: Vec<f64> }

/// Matrix-free finest correction plus recursive symmetric sparse V-cycles.
/// Only the bottom level uses dense storage. The immutable borrow prevents
/// ordinary fine-operator mutation while this numerical preparation is alive.
/// Applies allocate temporary vectors and one full cycle is not preemptible;
/// this is a CPU reference path, not allocation-free or hard-real-time code.
pub struct SparseMultilevel<'a, A: LinearOp> {
    operator: &'a A,
    inverse_diagonal: Vec<f64>,
    transfers: Vec<Transfer>,
    levels: Vec<Level>,
    bottom: Cholesky,
    equilibration: Vec<f64>,
    sizes: Vec<usize>,
    work: MultilevelWork,
}
impl<'a, A: LinearOp> SparseMultilevel<'a, A> {
    /// Transfers are ordered finest-to-coarsest. `first_coarse` must be the
    /// ACTUAL P[0]^T A_f P[0], formed by the physics owner from the current
    /// stiffness. The source operator is never probed during this constructor.
    /// All remaining Galerkin matrices are formed internally. Only the last
    /// transfer's column dimension is subject to the small dense-factor cap.
    pub fn new(operator: &'a A, inverse_diagonal: &[f64], prolongations: Vec<Csr>, first_coarse: Csr,
        control: &mut MultilevelControl<'_>) -> Result<Self, MultilevelError> {
        control.poll()?;
        let n = operator.n();
        if n == 0 || inverse_diagonal.len() != n || inverse_diagonal.iter().any(|d| !d.is_finite() || *d <= 0.0) {
            return Err(MultilevelError::Invalid("positive fine diagonal/dimension required"));
        }
        if n > control.budget.max_fine_dofs || prolongations.is_empty()
            || prolongations.len() >= control.budget.max_levels {
            return Err(MultilevelError::Budget("finest size or hierarchy depth"));
        }
        let mut expected = n; let mut entries = 0usize; let mut sizes = vec![n];
        for p in &prolongations {
            if p.nrows() != expected { return Err(MultilevelError::Invalid("disconnected transfer dimensions")); }
            transfer(p, control)?;
            entries = entries.checked_add(p.nnz()).filter(|&v| v <= control.budget.max_transfer_entries)
                .ok_or(MultilevelError::Budget("total transfer entries"))?;
            expected = p.ncols(); sizes.push(expected);
        }
        if expected > control.budget.max_coarsest_dofs { return Err(MultilevelError::Budget("bottom factor size; provide another level")); }
        if first_coarse.nrows() != sizes[1] { return Err(MultilevelError::Invalid("first Galerkin dimension")); }
        let mut transfers = Vec::with_capacity(prolongations.len());
        for p in prolongations { control.poll()?; let pt = transpose(&p); transfers.push(Transfer { p, pt }); }
        let mut current = first_coarse; let mut levels = Vec::with_capacity(transfers.len());
        for level in 0..transfers.len() {
            let diagonal = matrix(&current, control)?;
            control.retain(current.nnz())?;
            let next = if level + 1 < transfers.len() {
                Some(sparse_galerkin(&current, &transfers[level+1].p, control)?)
            } else { None };
            levels.push(Level { a: current, diagonal });
            match next { Some(a) => current = a, None => break }
        }
        let last = levels.last().ok_or(MultilevelError::Invalid("empty coarse hierarchy"))?;
        let nb = last.a.nrows();
        let roots: Vec<f64> = last.diagonal.iter().map(|v| v.sqrt()).collect();
        let mut dense = vec![0.0; nb * nb];
        for i in 0..nb {
            control.poll()?;
            for j in 0..=i {
                let value = (last.a.get(i, j) / roots[i]) / roots[j];
                if !value.is_finite() { return Err(MultilevelError::Invalid("bottom equilibration overflow")); }
                dense[i*nb+j] = value;
            }
        }
        control.poll()?;
        let bottom = cholesky(&dense, nb).map_err(|_| MultilevelError::NotPositiveDefinite)?;
        let equilibration: Vec<f64> = roots.iter().map(|v| 1.0 / v).collect();
        if equilibration.iter().any(|v| !v.is_finite() || *v <= 0.0) { return Err(MultilevelError::Invalid("bottom inverse scale overflow")); }
        control.poll()?;
        Ok(Self { operator, inverse_diagonal: inverse_diagonal.to_vec(), transfers,
            levels, bottom, equilibration, sizes, work: control.work() })
    }
    #[must_use] pub fn operator(&self) -> &A { self.operator }
    /// Includes the matrix-free finest level first.
    #[must_use] pub fn level_sizes(&self) -> &[usize] { &self.sizes }
    #[must_use] pub const fn work(&self) -> MultilevelWork { self.work }
    fn cycle(&self, level: usize, rhs: &[f64]) -> Vec<f64> {
        let a = &self.levels[level].a; let n = a.nrows();
        if level + 1 == self.levels.len() {
            let mut x: Vec<f64> = rhs.iter().zip(&self.equilibration).map(|(r, d)| r*d).collect();
            self.bottom.solve(&mut x);
            for (v, d) in x.iter_mut().zip(&self.equilibration) { *v *= d; }
            return x;
        }
        let diagonal = &self.levels[level].diagonal;
        let mut x = vec![0.0; n];
        // Solve (D+L)x=r: one forward GS sweep from zero.
        for i in 0..n {
            let mut value = rhs[i]; let (cols, values) = a.row(i);
            for (&j, &v) in cols.iter().zip(values) { if j < i { value -= v*x[j]; } }
            x[i] = value / diagonal[i];
        }
        let mut residual = vec![0.0; n]; a.spmv(&x, &mut residual);
        for (v, r) in residual.iter_mut().zip(rhs) { *v = r - *v; }
        let t = &self.transfers[level+1];
        let mut coarse_rhs = vec![0.0; t.p.ncols()]; t.pt.spmv(&residual, &mut coarse_rhs);
        let coarse = self.cycle(level+1, &coarse_rhs);
        let mut delta = vec![0.0; n]; t.p.spmv(&coarse, &mut delta);
        for (v, d) in x.iter_mut().zip(&delta) { *v += d; }
        a.spmv(&x, &mut residual);
        for (v, r) in residual.iter_mut().zip(rhs) { *v = r - *v; }
        // Adjoint post-sweep: (D+U)delta = r-Ax, never another forward sweep.
        delta.fill(0.0);
        for i in (0..n).rev() {
            let mut value = residual[i]; let (cols, values) = a.row(i);
            for (&j, &v) in cols.iter().zip(values) { if j > i { value -= v*delta[j]; } }
            delta[i] = value / diagonal[i];
        }
        for (v, d) in x.iter_mut().zip(delta) { *v += d; }
        x
    }
}
impl<A: LinearOp> Precond for SparseMultilevel<'_, A> {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        assert_eq!(r.len(), self.operator.n()); assert_eq!(z.len(), r.len());
        let t = &self.transfers[0];
        let mut coarse_rhs = vec![0.0; t.p.ncols()]; t.pt.spmv(r, &mut coarse_rhs);
        let coarse = self.cycle(0, &coarse_rhs); t.p.spmv(&coarse, z);
        for ((z, r), d) in z.iter_mut().zip(r).zip(&self.inverse_diagonal) { *z += d*r; }
    }
}
