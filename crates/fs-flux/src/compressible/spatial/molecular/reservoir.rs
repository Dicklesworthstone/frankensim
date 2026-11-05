//! Tangential derivatives of the selected source-admitted bath traces.

use super::*;

/// Midpoint tangential derivatives in normal-frame component order [n, t1, t2].
/// The cell cache is indexed [velocity component][derivative axis]. A uniform
/// bath contributes zero tangential derivatives. Face profiles use centered
/// differences, periodic wrapping or quadratic one-sided end differences;
/// two faces use a linear difference and one face resolves no variation.
#[allow(clippy::too_many_arguments)]
pub(super) fn tangential_gradient(
    grid: &CartesianGrid,
    boundaries: Boundaries3d,
    reservoirs: super::super::reservoir::View<'_>,
    axis: usize,
    sign: f64,
    cell: usize,
    interior_gradient: &[[f64; 3]; 3],
) -> Result<[[f64; 3]; 2], MolecularError> {
    let [a, b] = tangents(axis);
    let c = grid.coordinates(cell);
    let mut result = [[0.0; 3]; 2];
    for (row, derivative) in [a, b].into_iter().enumerate() {
        let n = grid.shape[derivative];
        let i = c[derivative];
        let value = |index: usize| -> Result<[f64; 3], MolecularError> {
            let mut p = c;
            p[derivative] = index;
            reservoirs
                .get_at(grid, axis, sign, grid.index(p))
                .map(|q| q.primitive().velocity)
                .ok_or(MolecularError::Invalid(
                    "unbound tangential molecular reservoir trace",
                ))
        };
        let boundary_gradient = if n == 1 {
            [0.0; 3]
        } else if boundaries.axis(derivative) == AxisBoundaries::Periodic {
            let low = value(if i == 0 { n - 1 } else { i - 1 })?;
            let high = value(if i + 1 == n { 0 } else { i + 1 })?;
            std::array::from_fn(|d| (high[d] - low[d]) / (2.0 * grid.spacing[derivative]))
        } else if n == 2 {
            let low = value(0)?;
            let high = value(1)?;
            std::array::from_fn(|d| (high[d] - low[d]) / grid.spacing[derivative])
        } else if i == 0 {
            let at = value(0)?;
            let next = value(1)?;
            let far = value(2)?;
            std::array::from_fn(|d| {
                (4.0 * (next[d] - at[d]) - (far[d] - at[d])) / (2.0 * grid.spacing[derivative])
            })
        } else if i + 1 == n {
            let at = value(i)?;
            let next = value(i - 1)?;
            let far = value(i - 2)?;
            std::array::from_fn(|d| {
                (4.0 * (at[d] - next[d]) - (at[d] - far[d])) / (2.0 * grid.spacing[derivative])
            })
        } else {
            let low = value(i - 1)?;
            let high = value(i + 1)?;
            std::array::from_fn(|d| (high[d] - low[d]) / (2.0 * grid.spacing[derivative]))
        };
        for (column, component) in [axis, a, b].into_iter().enumerate() {
            result[row][column] =
                0.5 * interior_gradient[component][derivative] + 0.5 * boundary_gradient[component];
        }
    }
    if !result.iter().flatten().all(|x| x.is_finite()) {
        return Err(MolecularError::Invalid(
            "tangential molecular reservoir gradient overflow",
        ));
    }
    Ok(result)
}
