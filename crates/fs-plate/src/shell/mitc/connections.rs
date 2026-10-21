//! Bilateral physical attachment resultants on the original shell coordinates.
//!
//! One signed row b measures displacement, and the SAME row returns its force.
//! Spring energy is k*(b.u)^2/2; viscous power is c*(b.v)^2. Stiffness is added
//! before the shared eigensolve, including all cross terms and support frames.
//! These finite resultants are neither exact ties nor a continuum mortar law.

use super::MitcShellModel;
use crate::PlateError;
use fs_sparse::Coo;

/// Largest number of physical modes one connection projection admits. It bounds the dense per-connection
/// rows; the caller's `max_work` still governs cost. Callers that retain more modes than this must refuse
/// before they reach it (fs-couple ties its own shell ceiling to this constant at compile time).
pub const MAX_CONNECTION_MODES: usize = 256;

/// Whether `n` supplied physical modes fit one projection: at least one, at most [`MAX_CONNECTION_MODES`].
#[must_use]
pub const fn connection_modes_admitted(n: usize) -> bool {
    n >= 1 && n <= MAX_CONNECTION_MODES
}

/// One contribution to a work-conjugate displacement coordinate.
#[derive(Clone, Debug, PartialEq)]
pub struct ShellConnectionTerm {
    /// Original physical midsurface node.
    pub node: usize,
    /// Global force fractions multiplying nodal translation.
    pub translation: [f64; 3],
    /// Global moment/force lever [m] multiplying physical axial rotation.
    pub rotation_m: [f64; 3],
}

/// One linear bilateral spring/dashpot at its stress-free reference state.
#[derive(Clone, Debug, PartialEq)]
pub struct ShellConnection {
    /// Complete signed footprint; equal nodes accumulate, never overwrite.
    pub terms: Vec<ShellConnectionTerm>,
    /// Nonnegative stiffness [N/m].
    pub stiffness_n_m: f64,
    /// Nonnegative viscous resistance [N s/m].
    pub damping_n_s_m: f64,
}

/// Cold connection work is independent of mesh and retained-mode allowances.
#[derive(Clone, Copy, Debug)]
pub struct ShellConnectionBudget {
    /// Maximum physical resultants, at most 256.
    pub max_connections: usize,
    /// Maximum total nodal footprint terms, at most 4096.
    pub max_terms: usize,
    /// Bound on scalar terms in each assembly or modal projection operation.
    pub max_work: usize,
}

fn invalid(what: &'static str) -> PlateError {
    PlateError::BadConstraint { what }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), PlateError> {
    if cancelled() { Err(PlateError::Cancelled) } else { Ok(()) }
}

fn admit(connections: &[ShellConnection], nodes: usize, budget: ShellConnectionBudget)
    -> Result<usize, PlateError>
{
    if budget.max_connections > 256 || budget.max_terms > 4096
        || connections.len() > budget.max_connections {
        return Err(invalid("shell connection count exceeds explicit allowance"));
    }
    let mut terms = 0usize;
    for link in connections {
        terms = terms.checked_add(link.terms.len())
            .ok_or_else(|| invalid("shell connection footprint extent overflow"))?;
        if link.terms.is_empty() || terms > budget.max_terms
            || !link.stiffness_n_m.is_finite() || link.stiffness_n_m < 0.0
            || !link.damping_n_s_m.is_finite() || link.damping_n_s_m < 0.0
            || (link.stiffness_n_m == 0.0 && link.damping_n_s_m == 0.0)
            || link.terms.iter().any(|t| t.node >= nodes
                || t.translation.iter().chain(&t.rotation_m).any(|x| !x.is_finite())) {
            return Err(invalid("shell connection needs finite passive coefficients and bounded physical terms"));
        }
    }
    Ok(terms)
}

fn physical_row(link: &ShellConnection, nodes: usize) -> Result<Vec<f64>, PlateError> {
    let mut row = vec![0.0; 6 * nodes];
    for term in &link.terms {
        for c in 0..3 {
            row[6 * term.node + c] += term.translation[c];
            row[6 * term.node + 3 + c] += term.rotation_m[c];
        }
    }
    if row.iter().any(|x| !x.is_finite()) || row.iter().all(|x| *x == 0.0) {
        return Err(invalid("shell connection physical footprint cancels or overflows"));
    }
    Ok(row)
}

impl MitcShellModel {
    /// Add k*b*b^T in the exact existing free/support basis, preserving M.
    /// b is the same physical work-conjugate row used by viscous projection.
    /// This is a cold, transactional operation. A connection completely fixed
    /// by Dirichlet constraints refuses instead of silently disappearing.
    /// No prestress, connection mass, geometric stiffness or rigid tie is added.
    pub fn with_connections(&self, connections: &[ShellConnection],
        budget: ShellConnectionBudget, mut cancelled: impl FnMut() -> bool)
        -> Result<Self, PlateError>
    {
        poll(&mut cancelled)?;
        let nodes = self.frames.len();
        admit(connections, nodes, budget)?;
        if connections.is_empty() { return Ok(self.clone()); }
        let work = connections.iter().try_fold(self.pencil.k.nnz(), |work, link| {
            let width = link.terms.len().checked_mul(6)?;
            work.checked_add(width.checked_mul(width)?)?
                .checked_add(nodes.checked_mul(72)?)
        });
        if work.is_none_or(|w| w > budget.max_work) {
            return Err(invalid("complete shell connection stiffness exceeds work allowance"));
        }
        let mut k = Coo::new(self.pencil.free, self.pencil.free);
        for i in 0..self.pencil.free {
            poll(&mut cancelled)?;
            let (columns, values) = self.pencil.k.row(i);
            for (&j, &value) in columns.iter().zip(values) { k.push(i, j, value); }
        }
        for link in connections {
            poll(&mut cancelled)?;
            let full = physical_row(link, nodes)?;
            let row = self.project_load(&full)?;
            let entries: Vec<_> = row.iter().enumerate().filter(|(_, v)| **v != 0.0).collect();
            if entries.is_empty() {
                return Err(invalid("shell connection is detached by its homogeneous supports"));
            }
            for &(i, &a) in &entries {
                poll(&mut cancelled)?;
                for &(j, &b) in &entries {
                    let value = link.stiffness_n_m * a * b;
                    if !value.is_finite() {
                        return Err(invalid("shell connection stiffness overflows"));
                    }
                    if value != 0.0 { k.push(i, j, value); }
                }
            }
        }
        poll(&mut cancelled)?;
        let mut out = self.clone();
        out.pencil.k = k.assemble();
        for i in 0..out.pencil.free {
            poll(&mut cancelled)?;
            if out.pencil.k.row(i).1.iter().any(|x| !x.is_finite()) {
                return Err(invalid("assembled shell connection stiffness overflows"));
            }
        }
        Ok(out)
    }

    /// Project every original signed footprint into supplied physical GLOBAL
    /// modes, preserving source order, sign and all coordinates. The resulting
    /// row realizes c*b*b^T with one passive channel, including semidefinite
    /// losses and exactly zero modal participation. No diagonal fit is applied.
    pub fn project_connections(&self, connections: &[ShellConnection],
        physical_modes: &[Vec<f64>], budget: ShellConnectionBudget,
        mut cancelled: impl FnMut() -> bool) -> Result<Vec<Vec<f64>>, PlateError>
    {
        poll(&mut cancelled)?;
        let nodes = self.frames.len();
        admit(connections, nodes, budget)?;
        if !connection_modes_admitted(physical_modes.len())
            || physical_modes.iter().any(|p| p.len() != 6 * nodes
                || p.iter().any(|v| !v.is_finite())) {
            return Err(invalid("shell connection projection needs bounded full physical modes"));
        }
        let work = nodes.checked_mul(6)
            .and_then(|w| w.checked_mul(physical_modes.len() + 1))
            .and_then(|w| w.checked_mul(connections.len()));
        if work.is_none_or(|w| w > budget.max_work) {
            return Err(invalid("complete shell connection projection exceeds work allowance"));
        }
        let mut out = Vec::with_capacity(connections.len());
        for link in connections {
            poll(&mut cancelled)?;
            let full = physical_row(link, nodes)?;
            let mut row = Vec::with_capacity(physical_modes.len());
            for mode in physical_modes {
                poll(&mut cancelled)?;
                let value = full.iter().zip(mode).map(|(a, b)| a * b).sum::<f64>();
                if !value.is_finite() {
                    return Err(invalid("shell connection modal projection overflows"));
                }
                row.push(value);
            }
            out.push(row);
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests;
