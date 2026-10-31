//! Shared algebra and electrode geometry for scalar-potential field owners.
//!
//! This layer has no constitutive units: the caller supplies its assembled
//! Hodge operator and load, then recovers its own charge/energy or current/heat.
//! Conductivity must never be passed off as dielectric permittivity.

use fs_rep_mesh::TriComplex2;
use fs_solver::{CheckedCgConfig, CheckedCgError, CsrOp, SolveReport, checked_cg};
use fs_sparse::{Csr, precond::IdentityPrecond};

use crate::planar::{MagneticError, PlanarMagnetostatics};

pub(crate) enum PotentialError {
    Invalid(&'static str),
    Cancelled,
    Nonconvergence(SolveReport),
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), PotentialError> {
    if cancelled() {
        Err(PotentialError::Cancelled)
    } else {
        Ok(())
    }
}

fn finite(value: f64) -> Result<f64, PotentialError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(PotentialError::Invalid(
            "scalar-potential arithmetic is not representable",
        ))
    }
}

/// Admit complete exterior traces and one anchor per mesh component through
/// the original geometric owner. This does not invoke a magnetic material law.
pub(crate) fn prescribed_boundary(
    mesh: &TriComplex2,
    boundary: &[(usize, f64)],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<Option<f64>>, PotentialError> {
    PlanarMagnetostatics::admit_mesh_boundary(mesh, boundary, cancelled).map_err(
        |error| match error {
            MagneticError::Cancelled => PotentialError::Cancelled,
            MagneticError::Invalid(reason) => PotentialError::Invalid(reason),
            _ => PotentialError::Invalid("invalid scalar-potential mesh"),
        },
    )
}

/// Eliminate the original prescribed coordinates and use the existing checked
/// CG recurrence, including its independently recomputed Euclidean residual.
pub(crate) fn solve(
    stiffness: &Csr,
    rhs: &[f64],
    prescribed: &[Option<f64>],
    controls: CheckedCgConfig,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<f64>, Option<SolveReport>), PotentialError> {
    poll(cancelled)?;
    let n = stiffness.nrows();
    if stiffness.ncols() != n
        || rhs.len() != n
        || prescribed.len() != n
        || rhs.iter().any(|v| !v.is_finite())
        || prescribed.iter().flatten().any(|v| !v.is_finite())
    {
        return Err(PotentialError::Invalid(
            "complete finite scalar-potential loads and prescribed values required",
        ));
    }
    let free: Vec<_> = (0..n).filter(|&v| prescribed[v].is_none()).collect();
    let mut index = vec![usize::MAX; n];
    for (i, &v) in free.iter().enumerate() {
        index[v] = i;
    }
    let mut row_ptr = vec![0];
    let mut columns = Vec::new();
    let mut values = Vec::new();
    let mut reduced_rhs = Vec::with_capacity(free.len());
    // Constant potentials lie in the full Hodge operator's nullspace. Solve
    // for differences before multiplying by that operator: a large common
    // voltage must not dominate either elimination or the CG residual scale.
    // The original absolute prescribed outputs remain bit-for-bit unchanged.
    let reference = prescribed.iter().flatten().copied().next().unwrap_or(0.0);
    for &v in &free {
        poll(cancelled)?;
        let mut load = rhs[v];
        let (cols, vals) = stiffness.row(v);
        for (&w, &value) in cols.iter().zip(vals) {
            if let Some(voltage) = prescribed[w] {
                let difference = finite(voltage - reference)?;
                load = finite((-value).mul_add(difference, load))?;
            } else {
                columns.push(index[w]);
                values.push(value);
            }
        }
        row_ptr.push(columns.len());
        reduced_rhs.push(load);
    }
    let matrix = Csr::from_parts(free.len(), free.len(), row_ptr, columns, values);
    let mut voltage: Vec<_> = prescribed.iter().map(|v| v.unwrap_or(0.0)).collect();
    let report = if free.is_empty() {
        None
    } else {
        let result = checked_cg(
            &CsrOp::symmetric(matrix),
            &IdentityPrecond,
            &reduced_rhs,
            controls,
            |_| poll(cancelled),
        )
        .map_err(|error| match error {
            CheckedCgError::InvalidInput(reason) => PotentialError::Invalid(reason),
            CheckedCgError::Interrupted(error) => error,
        })?;
        if !result.report.converged_euclidean() {
            return Err(PotentialError::Nonconvergence(result.report));
        }
        for (&v, &value) in free.iter().zip(&result.x) {
            voltage[v] = finite(reference + value)?;
        }
        Some(result.report)
    };
    poll(cancelled)?;
    Ok((voltage, report))
}

/// Electrode labels must refer to disjoint unions of complete exterior edges,
/// not point contacts. The field owner subsequently admits component anchoring.
pub(crate) fn terminal_boundary(
    mesh: &TriComplex2,
    terminals: &[&[usize]],
    voltages_v: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<(usize, f64)>, PotentialError> {
    let n = mesh.vertices().len();
    if terminals.len() < 2
        || terminals.len() > n
        || voltages_v.len() != terminals.len()
        || voltages_v.iter().any(|v| !v.is_finite())
    {
        return Err(PotentialError::Invalid(
            "complete finite electrode voltages and at least two terminals required",
        ));
    }
    let trace = mesh
        .boundary_trace()
        .map_err(|_| PotentialError::Invalid("invalid electrode boundary trace"))?;
    let mut boundary = Vec::new();
    let mut assigned = vec![false; n];
    for (&terminal, &value) in terminals.iter().zip(voltages_v) {
        if terminal.is_empty() || terminal.len() > n {
            return Err(PotentialError::Invalid(
                "nonempty bounded terminal traces required",
            ));
        }
        let mut selected = vec![false; n];
        for &vertex in terminal {
            poll(cancelled)?;
            if vertex >= n || assigned[vertex] {
                return Err(PotentialError::Invalid(
                    "invalid, overlapping or repeated terminal vertex",
                ));
            }
            assigned[vertex] = true;
            selected[vertex] = true;
            boundary.push((vertex, value));
        }
        let mut supported = vec![false; n];
        for edge in &trace.edges {
            poll(cancelled)?;
            let [a, b] = edge.oriented_vertices.map(|v| v as usize);
            if selected[a] && selected[b] {
                supported[a] = true;
                supported[b] = true;
            }
        }
        if terminal.iter().any(|&v| !supported[v]) {
            return Err(PotentialError::Invalid(
                "each terminal must cover exterior edges, not isolated points",
            ));
        }
    }
    poll(cancelled)?;
    Ok(boundary)
}

/// Physical components after ideal equipotential terminals join mesh pieces.
/// This owns connectivity only; the consumer must anchor every returned
/// component with a voltage control, never a guessed floating voltage.
pub(crate) fn terminal_components(
    mesh: &TriComplex2,
    terminals: &[&[usize]],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<usize>, PotentialError> {
    fn root(parent: &mut [usize], mut v: usize) -> usize {
        while parent[v] != v {
            parent[v] = parent[parent[v]];
            v = parent[v];
        }
        v
    }
    fn join(parent: &mut [usize], a: usize, b: usize) {
        let a = root(parent, a);
        let b = root(parent, b);
        parent[a.max(b)] = a.min(b);
    }
    let n = mesh.vertices().len();
    let mut parent: Vec<_> = (0..n).collect();
    for &[a, b] in mesh.edges() {
        poll(cancelled)?;
        join(&mut parent, a as usize, b as usize);
    }
    for terminal in terminals {
        if terminal.is_empty() || terminal.iter().any(|&v| v >= n) {
            return Err(PotentialError::Invalid("complete valid terminal traces required"));
        }
        for &v in *terminal {
            poll(cancelled)?;
            join(&mut parent, terminal[0], v);
        }
    }
    for v in 0..n {
        poll(cancelled)?;
        parent[v] = root(&mut parent, v);
    }
    Ok(parent)
}
