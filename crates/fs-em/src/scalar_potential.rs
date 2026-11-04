//! Shared algebra and electrode geometry for scalar-potential field owners.
//!
//! This layer has no constitutive units: the caller supplies its assembled
//! Hodge operator and load, then recovers its own charge/energy or current/heat.
//! Conductivity must never be passed off as dielectric permittivity.

use fs_rep_mesh::TriComplex2;
use fs_solver::{CheckedCgConfig, CheckedCgError, CsrOp, SolveReport, checked_cg};
use fs_sparse::{Csr, precond::Precond};

use crate::planar::{MagneticError, PlanarMagnetostatics};

pub(crate) enum PotentialError {
    Invalid(&'static str),
    Cancelled,
    Nonconvergence(SolveReport),
}

/// Native Jacobi application on the already admitted free-coordinate matrix.
/// The original matrix/RHS stay unscaled; checked_cg still accepts only its
/// independently recomputed Euclidean residual in the physical equation.
struct Diagonal(Vec<f64>);
impl Precond for Diagonal {
    fn apply(&self, residual: &[f64], output: &mut [f64]) {
        for ((out, &value), &inverse) in output.iter_mut().zip(residual).zip(&self.0) {
            *out = value * inverse;
        }
    }
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

/// One deterministic prescribed reference for each connected scalar operator.
/// Disconnected bodies have independent gauges: carrying another component's
/// common voltage into a reduced right-hand side can hide its physical load
/// below the Krylov residual scale. Nonzero entries include material contact
/// couplings, so electrically joined mesh pieces share the same reference.
pub(crate) fn component_references(
    stiffness: &Csr,
    prescribed: &[Option<f64>],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>, PotentialError> {
    fn root(parent: &mut [usize], mut vertex: usize) -> usize {
        while parent[vertex] != vertex {
            parent[vertex] = parent[parent[vertex]];
            vertex = parent[vertex];
        }
        vertex
    }
    let n = stiffness.nrows();
    if stiffness.ncols() != n
        || prescribed.len() != n
        || prescribed.iter().flatten().any(|value| !value.is_finite())
    {
        return Err(PotentialError::Invalid(
            "complete finite scalar-potential references required",
        ));
    }
    let mut components: Vec<_> = (0..n).collect();
    for vertex in 0..n {
        poll(cancelled)?;
        let (columns, values) = stiffness.row(vertex);
        for (&column, &value) in columns.iter().zip(values) {
            if column != vertex && value != 0.0 {
                let a = root(&mut components, vertex);
                let b = root(&mut components, column);
                components[a.max(b)] = a.min(b);
            }
        }
    }
    let mut anchors = vec![None; n];
    for vertex in 0..n {
        poll(cancelled)?;
        let component = root(&mut components, vertex);
        components[vertex] = component;
        if anchors[component].is_none() {
            anchors[component] = prescribed[vertex];
        }
    }
    let mut references = Vec::with_capacity(n);
    for component in components {
        poll(cancelled)?;
        references.push(anchors[component].ok_or(PotentialError::Invalid(
            "every scalar-potential component needs a prescribed reference",
        ))?);
    }
    Ok(references)
}

/// Eliminate the original prescribed coordinates and use the existing checked
/// Jacobi-CG recurrence, including its independently recomputed Euclidean
/// residual on the unscaled physical equation. Diagonal preparation shares
/// the original row traversal and retains one scalar per free coordinate.
pub(crate) fn solve(
    stiffness: &Csr,
    rhs: &[f64],
    prescribed: &[Option<f64>],
    controls: CheckedCgConfig,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<f64>, Option<SolveReport>), PotentialError> {
    let solved = solve_relative(stiffness, rhs, prescribed, controls, cancelled)?;
    let mut absolute = Vec::with_capacity(solved.relative.len());
    for (vertex, &relative) in solved.relative.iter().enumerate() {
        poll(cancelled)?;
        absolute.push(match prescribed[vertex] {
            Some(value) => value,
            None => finite(solved.references[vertex] + relative)?,
        });
    }
    Ok((absolute, solved.report))
}

/// Preserve physical potential differences separately from their gauges. An
/// induced free potential can be far below one ulp of a prescribed common
/// voltage; restoring that voltage before field recovery would erase it.
pub(crate) struct RelativePotential {
    pub relative: Vec<f64>,
    pub references: Vec<f64>,
    pub report: Option<SolveReport>,
}

pub(crate) fn solve_relative(
    stiffness: &Csr,
    rhs: &[f64],
    prescribed: &[Option<f64>],
    controls: CheckedCgConfig,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<RelativePotential, PotentialError> {
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
    let mut inverse_diagonal = Vec::with_capacity(free.len());
    // Constant potentials lie in the full Hodge operator's nullspace. Solve
    // for differences before multiplying by that operator: a large common
    // voltage must not dominate either elimination or the CG residual scale.
    // The original absolute prescribed outputs remain bit-for-bit unchanged.
    let references = component_references(stiffness, prescribed, cancelled)?;
    for &v in &free {
        poll(cancelled)?;
        let mut load = rhs[v];
        let mut diagonal = 0.0;
        let (cols, vals) = stiffness.row(v);
        for (&w, &value) in cols.iter().zip(vals) {
            if w == v {
                diagonal = value;
            }
            if let Some(voltage) = prescribed[w] {
                let difference = finite(voltage - references[v])?;
                load = finite((-value).mul_add(difference, load))?;
            } else {
                columns.push(index[w]);
                values.push(value);
            }
        }
        if !diagonal.is_finite() || diagonal <= 0.0 || !diagonal.recip().is_finite() {
            return Err(PotentialError::Invalid(
                "scalar-potential Jacobi diagonal must be representably positive",
            ));
        }
        inverse_diagonal.push(diagonal.recip());
        row_ptr.push(columns.len());
        reduced_rhs.push(load);
    }
    let matrix = Csr::from_parts(free.len(), free.len(), row_ptr, columns, values);
    let mut relative = prescribed
        .iter()
        .zip(&references)
        .map(|(value, &reference)| match value {
            Some(value) => finite(value - reference),
            None => Ok(0.0),
        })
        .collect::<Result<Vec<_>, _>>()?;
    let report = if free.is_empty() {
        None
    } else {
        let result = checked_cg(
            &CsrOp::symmetric(matrix),
            &Diagonal(inverse_diagonal),
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
            relative[v] = value;
        }
        Some(result.report)
    };
    poll(cancelled)?;
    Ok(RelativePotential {
        relative,
        references,
        report,
    })
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
            return Err(PotentialError::Invalid(
                "complete valid terminal traces required",
            ));
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

#[cfg(test)]
mod tests {
    use super::*;
    use fs_sparse::precond::IdentityPrecond;
    use std::convert::Infallible;

    fn controls() -> CheckedCgConfig {
        CheckedCgConfig {
            tolerance: 1e-12,
            max_iterations: 256,
            max_corrections: 2,
        }
    }

    /// Exact capacitor star: node zero is the common electrode and leaf i has
    /// independent C_i. The full operator has its physical constant nullspace;
    /// after grounding zero, Q_i=C_i V_i is an independent analytical solution.
    fn capacitor_star(capacitances: &[f64]) -> Csr {
        let n = capacitances.len() + 1;
        let mut row_ptr = vec![0];
        let mut columns = vec![0];
        let mut values = vec![capacitances.iter().sum()];
        for (leaf, &capacitance) in capacitances.iter().enumerate() {
            columns.push(leaf + 1);
            values.push(-capacitance);
        }
        row_ptr.push(columns.len());
        for (leaf, &capacitance) in capacitances.iter().enumerate() {
            columns.extend([0, leaf + 1]);
            values.extend([-capacitance, capacitance]);
            row_ptr.push(columns.len());
        }
        Csr::from_parts(n, n, row_ptr, columns, values)
    }

    #[test]
    fn jacobi_resolves_heterogeneous_capacitances_with_fewer_native_cg_iterations() {
        let capacitances: Vec<_> = (0..16).map(|i| 2_f64.powi(2 * i - 40)).collect();
        let matrix = capacitor_star(&capacitances);
        let mut rhs = vec![0.0];
        rhs.extend_from_slice(&capacitances);
        let mut prescribed = vec![None; rhs.len()];
        prescribed[0] = Some(0.0);
        let actual = solve_relative(&matrix, &rhs, &prescribed, controls(), &mut || false)
            .unwrap_or_else(|_| panic!("physical capacitor-star solve refused"));
        let report = actual.report.as_ref().unwrap();
        assert!(report.converged_euclidean());
        for &voltage in &actual.relative[1..] {
            assert!((voltage - 1.0).abs() < 1e-12);
        }
        // Check the unscaled physical equation independently of the
        // preconditioned recurrence or its recursive residual estimate.
        let residual_squared: f64 = capacitances
            .iter()
            .zip(&actual.relative[1..])
            .map(|(&c, &voltage)| (c - c * voltage).powi(2))
            .sum();
        let rhs_squared: f64 = capacitances.iter().map(|c| c * c).sum();
        assert!((residual_squared / rhs_squared).sqrt() < controls().tolerance);

        let count = capacitances.len();
        let reduced = Csr::from_parts(
            count,
            count,
            (0..=count).collect(),
            (0..count).collect(),
            capacitances,
        );
        let baseline = checked_cg(
            &CsrOp::symmetric(reduced),
            &IdentityPrecond,
            &rhs[1..],
            controls(),
            |_| Ok::<_, Infallible>(()),
        )
        .unwrap();
        assert!(baseline.report.converged_euclidean());
        assert!(report.iters <= 2, "Jacobi iterations: {}", report.iters);
        assert!(baseline.report.iters >= 8);
        assert!(report.iters < baseline.report.iters);

        // A common absolute bias must not alter either physical differences
        // or the deterministic iteration count and measured residual.
        prescribed[0] = Some(2_f64.powi(45));
        let shifted = solve_relative(&matrix, &rhs, &prescribed, controls(), &mut || false)
            .unwrap_or_else(|_| panic!("gauge-shifted capacitor-star solve refused"));
        assert_eq!(actual.relative, shifted.relative);
        assert_eq!(report.iters, shifted.report.as_ref().unwrap().iters);
        assert_eq!(
            report.rel_residual.to_bits(),
            shifted.report.as_ref().unwrap().rel_residual.to_bits()
        );
    }

    #[test]
    fn nonpositive_or_unrepresentable_jacobi_diagonals_refuse_before_solving() {
        for diagonal in [0.0, -1.0, f64::NAN, f64::INFINITY, 1e-320] {
            // Keep a connection to the reference even for malformed input so
            // this reaches diagonal admission rather than the gauge refusal.
            let matrix = Csr::from_parts(
                2,
                2,
                vec![0, 2, 4],
                vec![0, 1, 0, 1],
                vec![1.0, -1.0, -1.0, diagonal],
            );
            assert!(matches!(
                solve_relative(
                    &matrix,
                    &[0.0, 1.0],
                    &[Some(0.0), None],
                    controls(),
                    &mut || false,
                ),
                Err(PotentialError::Invalid(
                    "scalar-potential Jacobi diagonal must be representably positive"
                ))
            ));
        }
    }

    #[test]
    fn jacobi_setup_and_checked_iterations_remain_cancellable() {
        let matrix = capacitor_star(&[1e-12, 1e-9, 1e-6]);
        let rhs = [0.0, 1e-12, 2e-9, -1e-6];
        let prescribed = [Some(0.0), None, None, None];
        let mut calls = 0;
        let original = solve_relative(&matrix, &rhs, &prescribed, controls(), &mut || {
            calls += 1;
            false
        })
        .unwrap_or_else(|_| panic!("capacitor-star solve refused"));
        for stop in [1, 3, calls - 1, calls] {
            let mut count = 0;
            assert!(matches!(
                solve_relative(&matrix, &rhs, &prescribed, controls(), &mut || {
                    count += 1;
                    count == stop
                }),
                Err(PotentialError::Cancelled)
            ));
        }
        let replay = solve_relative(&matrix, &rhs, &prescribed, controls(), &mut || false)
            .unwrap_or_else(|_| panic!("capacitor-star replay refused"));
        assert_eq!(original.relative, replay.relative);
    }
}
