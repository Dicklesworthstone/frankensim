//! Original MAC stress work split between molecular heat and modeled SGS transfer.
//!
//! Each material subvolume sees the same harmonic edge stress. Its velocity
//! gradient is therefore inversely proportional to its TOTAL dynamic viscosity;
//! molecular and SGS losses share that local gradient, not two independently
//! harmonically averaged stress laws. No modeled transfer is deposited as heat.

use super::{ChtError, Solver, cell_of, finite_value, invalid_power};

/// Split one nonnegative original-volume stress power without subtracting
/// nearly equal extensive powers. A small positive molecular or SGS channel
/// therefore retains its own arithmetic.
pub(super) fn split_power(
    solver: &Solver<'_>,
    cell: usize,
    eddy: f64,
    power: f64,
) -> Result<[f64; 2], ChtError> {
    let molecular = solver.molecular_viscosity[cell];
    let modeled = solver.sgs_density(cell) * eddy;
    let total = solver.sgs_total_viscosity(cell, eddy);
    if !power.is_finite()
        || power < 0.0
        || !molecular.is_finite()
        || molecular <= 0.0
        || !modeled.is_finite()
        || modeled < 0.0
        || !total.is_finite()
        || total <= 0.0
    {
        return Err(invalid_power("invalid material SGS stress power partition"));
    }
    Ok([
        finite_value(
            "viscous_heat.molecular_partition",
            power * (molecular / total),
        )?,
        finite_value("viscous_heat.sgs_partition", power * (modeled / total))?,
    ])
}

/// Partition one original reciprocal edge power over its four material
/// subvolumes. The caller supplies the SAME native edge closure and power
/// used by both momentum components; areas retain the graded mesh geometry.
pub(super) fn edge_power(
    solver: &Solver<'_>,
    cells: [[usize; 3]; 4],
    axes: [usize; 2],
    eddy: f64,
    power: f64,
) -> Result<[[f64; 2]; 4], ChtError> {
    let [a, d] = axes;
    let areas = cells.map(|p| solver.domain.width(a, p[a]) * solver.domain.width(d, p[d]));
    let area: f64 = areas.iter().sum();
    let totals = cells.map(|p| solver.sgs_total_viscosity(cell_of(solver.n, p), eddy));
    let minimum = totals.iter().copied().fold(f64::INFINITY, f64::min);
    let weights: [f64; 4] = core::array::from_fn(|i| (areas[i] / area) * (minimum / totals[i]));
    let sum: f64 = weights.iter().sum();
    if !sum.is_finite() || sum <= 0.0 {
        return Err(invalid_power(
            "unrepresentable material SGS edge resistance partition",
        ));
    }
    let mut partition = [[0.0; 2]; 4];
    for i in 0..4 {
        partition[i] = split_power(
            solver,
            cell_of(solver.n, cells[i]),
            eddy,
            power * (weights[i] / sum),
        )?;
    }
    Ok(partition)
}

/// One reservoir-mouth half-edge has two material subvolumes. The common
/// normal half-width and third-axis width cancel from their resistance ratio;
/// only the original tangent widths and TOTAL dynamic viscosities remain.
pub(super) fn mouth_power(
    solver: &Solver<'_>,
    cells: [[usize; 3]; 2],
    tangent: usize,
    eddy: f64,
    power: f64,
) -> Result<[[f64; 2]; 2], ChtError> {
    let widths = cells.map(|p| solver.domain.width(tangent, p[tangent]));
    let width = widths[0] + widths[1];
    let totals = cells.map(|p| solver.sgs_total_viscosity(cell_of(solver.n, p), eddy));
    let minimum = totals[0].min(totals[1]);
    let weights: [f64; 2] = core::array::from_fn(|i| (widths[i] / width) * (minimum / totals[i]));
    let sum = weights[0] + weights[1];
    if !sum.is_finite() || sum <= 0.0 {
        return Err(invalid_power(
            "unrepresentable material SGS mouth resistance partition",
        ));
    }
    let mut partition = [[0.0; 2]; 2];
    for i in 0..2 {
        partition[i] = split_power(
            solver,
            cell_of(solver.n, cells[i]),
            eddy,
            power * (weights[i] / sum),
        )?;
    }
    Ok(partition)
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod pressure_tests;
