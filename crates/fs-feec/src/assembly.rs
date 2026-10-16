//! Assembly: materialize the exact integer incidence operators as
//! f64 CSR (±1 entries, per-row column sort for the CSR invariant)
//! and compose stiffness operators the FEEC way — K_k = d_kᵀ·M_{k+1}·d_k
//! through fs-sparse's deterministic transpose/spgemm. The exactness
//! of the sequence lives in the INTEGER d (fs-rep-mesh); everything
//! metric enters only through the mass matrix.

use fs_rep_mesh::Incidence;
use fs_sparse::{Csr, ops};

/// Materialize an integer incidence operator as an f64 CSR (entries
/// are exactly ±1.0; columns sorted per row for the canonical CSR
/// invariant — a pure reordering, no arithmetic).
#[must_use]
pub fn incidence_to_csr(inc: &Incidence) -> Csr {
    let mut row_ptr = Vec::with_capacity(inc.rows.len() + 1);
    let mut col_idx = Vec::new();
    let mut vals = Vec::new();
    row_ptr.push(0usize);
    for row in &inc.rows {
        let mut entries: Vec<(usize, i8)> = row.clone();
        entries.sort_unstable_by_key(|&(c, _)| c);
        for (c, s) in entries {
            col_idx.push(c);
            vals.push(f64::from(s));
        }
        row_ptr.push(col_idx.len());
    }
    Csr::from_parts(inc.rows.len(), inc.cols, row_ptr, col_idx, vals)
}

/// FEEC stiffness composition: K = dᵀ·M·d (e.g. the P₁ Poisson
/// stiffness is d₀ᵀ·M₁·d₀ — identical to classical FEM assembly, but
/// derived from the complex). Deterministic: Gustavson spgemm over
/// canonical CSR.
#[must_use]
pub fn stiffness(d: &Csr, mass_next: &Csr) -> Csr {
    let dt = ops::transpose(d);
    let md = ops::spgemm(mass_next, d);
    ops::spgemm(&dt, &md)
}
