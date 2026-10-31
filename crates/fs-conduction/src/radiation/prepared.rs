//! Reuse the constant radiosity operator without retaining a mutable global cache.
//!
//! An immutable enclosure fixes both emissivities and view factors. Its matrix
//! is independent of the changing surface temperatures. Store its elimination
//! operations once; applying them to each new RHS preserves the one-shot
//! solver's pivot, subtraction, and FMA order, including interleaved row swaps.

use super::{
    ConductionError, Cx, GrayDiffuseEnclosure, MAX_ENCLOSURE_SURFACES, RadiosityReport,
    radiation_error,
};

/// A prepared radiosity operator borrowed from one immutable enclosure.
///
/// Preparation performs one bounded, cancellable dense factorization. Each
/// subsequent [`Self::solve`] performs only quadratic work for assembly,
/// triangular solution, and the ORIGINAL radiosity/heat report. Surface
/// temperatures and their emissivity validity are checked on every call.
///
/// The borrowed enclosure prevents detaching factors from their surfaces,
/// emissivities, or view-factor order. This object retains no temperature or
/// mutable iteration history and may be reused across independent solves.
/// Its storage is O(n^2), with the existing 256-surface enclosure limit; this
/// is not a bound on total conduction workspace. No physical certificate,
/// continuum error estimate, or temperature-dependent emissivity is added.
#[derive(Debug)]
pub struct PreparedRadiosity<'a> {
    enclosure: &'a GrayDiffuseEnclosure,
    factors: DenseRadiosityFactorization,
}

impl<'a> PreparedRadiosity<'a> {
    pub(super) fn new(
        cx: &Cx<'_>,
        enclosure: &'a GrayDiffuseEnclosure,
    ) -> Result<Self, ConductionError> {
        let n = enclosure.surfaces.len();
        let mut matrix = vec![vec![0.0; n]; n];
        for (i, row) in matrix.iter_mut().enumerate() {
            poll(cx, "radiation-radiosity-assemble", i)?;
            let emissivity = enclosure.surfaces[i].emissivity.value;
            for (j, entry) in row.iter_mut().enumerate() {
                *entry = if i == j { 1.0 } else { 0.0 }
                    - (1.0 - emissivity) * enclosure.view_factors.factors[i][j];
            }
        }
        let factors = DenseRadiosityFactorization::new(cx, &matrix)?;
        Ok(Self { enclosure, factors })
    }

    /// Evaluate a fresh temperature vector with the retained factorization.
    ///
    /// The same original matrix is reassembled for the residual calculation;
    /// a factor-space residual is never substituted for physical heat records.
    /// This preserves the ordinary one-shot report and its stated limitations.
    pub fn solve(
        &self,
        cx: &Cx<'_>,
        surface_temperatures_k: &[f64],
    ) -> Result<RadiosityReport, ConductionError> {
        self.enclosure.solve_with_factors(cx, surface_temperatures_k, Some(&self.factors))
    }

    pub(crate) fn enclosure(&self) -> &GrayDiffuseEnclosure {
        self.enclosure
    }

    // Load and tangent owners validate their physical inputs and independently
    // check the original equations. Reuse the same factorization for each RHS.
    pub(crate) fn solve_rhs(&self, cx: &Cx<'_>, rhs: &[f64]) -> Result<Vec<f64>, ConductionError> {
        self.factors.solve(cx, rhs)
    }
}

/// Keep the multipliers in elimination-step order, not in a conventional
/// row-permuted lower triangle: replaying the original interleaved swaps is
/// what preserves the original RHS arithmetic exactly.
#[derive(Debug)]
pub(super) struct DenseRadiosityFactorization {
    upper: Vec<Vec<f64>>,
    pivots: Vec<usize>,
    multipliers: Vec<Vec<f64>>,
}

impl DenseRadiosityFactorization {
    pub(super) fn new(cx: &Cx<'_>, matrix: &[Vec<f64>]) -> Result<Self, ConductionError> {
        let n = matrix.len();
        if !(2..=MAX_ENCLOSURE_SURFACES).contains(&n)
            || matrix.iter().any(|row| row.len() != n || row.iter().any(|v| !v.is_finite()))
        {
            return Err(radiation_error(
                "radiosity-system", "invalid prepared matrix dimensions or entries",
                "prepare an admitted finite enclosure within the existing surface limit",
            ));
        }
        let mut upper = matrix.to_vec();
        let mut pivots = Vec::with_capacity(n);
        let mut multipliers = Vec::with_capacity(n);
        let scale = upper.iter().flatten().copied().map(f64::abs)
            .fold(0.0f64, f64::max).max(f64::MIN_POSITIVE);
        for column in 0..n {
            poll(cx, "radiation-radiosity-factor", column)?;
            let mut pivot = column;
            let mut pivot_abs = upper[column][column].abs();
            for (row, values) in upper.iter().enumerate().skip(column + 1) {
                let candidate = values[column].abs();
                if candidate > pivot_abs {
                    pivot = row;
                    pivot_abs = candidate;
                }
            }
            if !pivot_abs.is_finite() || pivot_abs <= 64.0 * f64::EPSILON * scale {
                return Err(radiation_error(
                    "radiosity-system",
                    format!("pivot {column} has magnitude {pivot_abs}"),
                    "repair the emissivity/view-factor model; the radiosity operator is singular or ill-conditioned",
                ));
            }
            pivots.push(pivot);
            if pivot != column { upper.swap(column, pivot); }
            let pivot_row = upper[column].clone();
            let mut column_multipliers = Vec::with_capacity(n - column - 1);
            for values in upper.iter_mut().skip(column + 1) {
                let factor = values[column] / pivot_row[column];
                column_multipliers.push(factor);
                values[column] = 0.0;
                for entry in (column + 1)..n {
                    values[entry] = (-factor).mul_add(pivot_row[entry], values[entry]);
                }
            }
            multipliers.push(column_multipliers);
        }
        poll(cx, "radiation-radiosity-factor", n)?;
        Ok(Self { upper, pivots, multipliers })
    }

    pub(super) fn solve(&self, cx: &Cx<'_>, rhs: &[f64]) -> Result<Vec<f64>, ConductionError> {
        let n = self.upper.len();
        if rhs.len() != n || rhs.iter().any(|v| !v.is_finite()) {
            return Err(radiation_error(
                "radiosity-system", "invalid prepared right-hand side",
                "supply one finite thermal load per prepared surface",
            ));
        }
        let mut b = rhs.to_vec();
        for column in 0..n {
            poll(cx, "radiation-radiosity-apply", column)?;
            let pivot = self.pivots[column];
            if pivot != column { b.swap(column, pivot); }
            let pivot_rhs = b[column];
            for (offset, &factor) in self.multipliers[column].iter().enumerate() {
                let row = column + offset + 1;
                b[row] = (-factor).mul_add(pivot_rhs, b[row]);
            }
        }
        let mut x = vec![0.0; n];
        for row in (0..n).rev() {
            poll(cx, "radiation-radiosity-backsolve", row)?;
            let mut value = b[row];
            for (column, &known) in x.iter().enumerate().skip(row + 1) {
                value = (-self.upper[row][column]).mul_add(known, value);
            }
            x[row] = value / self.upper[row][row];
            if !x[row].is_finite() {
                return Err(radiation_error(
                    "radiosity-system", format!("solution component {row} is not finite"),
                    "repair the enclosure inputs",
                ));
            }
        }
        poll(cx, "radiation-radiosity-backsolve", n)?;
        Ok(x)
    }
}

fn poll(cx: &Cx<'_>, stage: &'static str, at: usize) -> Result<(), ConductionError> {
    cx.checkpoint().map_err(|_| ConductionError::Cancelled { stage, at })
}
