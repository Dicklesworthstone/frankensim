//! Reciprocal multi-terminal capacitance on original capacitor coordinates.
//!
//! The declared matrix relates branch charge to branch voltage, `q = C v`.
//! It is not a list of positive pair capacitors: either off-diagonal sign is
//! admissible when the complete reciprocal matrix is positive definite.
//! A reference/gauge choice belongs to the caller. Singular Maxwell matrices
//! must retain an explicit reference before entering independent charge states.

use super::{Branch, CapacitorCurve, CircuitDae, CircuitError, CircuitGraph, invalid};
use std::sync::Arc;

/// Bounds preparation, inverse verification and every dense block query.
pub const MAX_CAPACITANCE_BLOCK_COORDINATES: usize = 128;

/// Prepared `H(q) = q^T C^-1 q / 2`, bound to original graph branches.
///
/// Charge and voltage ordering is exactly [`Self::branches`], including each
/// branch's first-to-second-terminal orientation. The full supplied C owns
/// its diagonal too; the graph's positive scalar coefficients are topology
/// placeholders on these branches. Immutable coefficients and the admitted
/// inverse are shared across clones and branch remapping. The block's C is
/// never refactored when assembling another native circuit; ordinary native
/// graph and combined-storage admission still apply.
#[derive(Clone, Debug, PartialEq)]
pub struct CapacitanceBlock {
    branches: Arc<[usize]>,
    matrix_f: Arc<[f64]>,
    inverse_per_f: Arc<[f64]>,
}

impl CapacitanceBlock {
    /// Prepare a bounded, exactly reciprocal, positive-definite C [F].
    ///
    /// The matrix is row-major in the supplied branch order. Original values
    /// are retained without symmetrization, clipping or fitting. Diagonal
    /// equilibration makes the positivity and inverse checks independent of
    /// capacitance units and electrode scale. A tiny unresolved pivot refuses
    /// rather than turning a remaining gauge into an extremely stiff store.
    pub fn try_new(branches: Vec<usize>, matrix_f: Vec<f64>) -> Result<Self, CircuitError> {
        admit_branches(&branches)?;
        let n = branches.len();
        if matrix_f.len() != n * n || matrix_f.iter().any(|v| !v.is_finite()) {
            return Err(invalid("complete finite capacitance block matrix required"));
        }
        let mut scales = Vec::with_capacity(n);
        for i in 0..n {
            if matrix_f[i * n + i] <= 0.0 {
                return Err(invalid("capacitance block diagonals must be positive"));
            }
            scales.push(matrix_f[i * n + i].sqrt());
            for j in 0..i {
                if matrix_f[i * n + j] != matrix_f[j * n + i] {
                    return Err(invalid("capacitance block must be exactly reciprocal"));
                }
            }
        }
        let mut scaled = vec![0.0; n * n];
        let mut products = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..=i {
                let product = finite(scales[i] * scales[j])?;
                if product <= 0.0 {
                    return Err(invalid("capacitance block scaling is unrepresentable"));
                }
                let value = finite(matrix_f[i * n + j] / product)?;
                if matrix_f[i * n + j] != 0.0 && value == 0.0 {
                    return Err(invalid(
                        "capacitance block coupling underflowed during scaling",
                    ));
                }
                scaled[i * n + j] = value;
                scaled[j * n + i] = value;
                products[i * n + j] = product;
                products[j * n + i] = product;
            }
        }
        let factor = fs_la::factor::cholesky(&scaled, n)
            .map_err(|_| invalid("full capacitance block must be positive definite"))?;
        let pivot_floor = 16.0 * f64::EPSILON * n as f64;
        if (0..n).any(|i| factor.l(i, i).powi(2) <= pivot_floor) {
            return Err(invalid(
                "capacitance block gauge or storage is numerically unresolved",
            ));
        }
        let mut inverse = vec![0.0; n * n];
        let mut column = vec![0.0; n];
        for j in 0..n {
            column.fill(0.0);
            column[j] = 1.0;
            factor.solve(&mut column);
            for i in 0..n {
                inverse[i * n + j] = finite(column[i])?;
            }
        }
        // Share a single rounded reciprocal inverse coefficient, then check
        // the actual represented physical inverse against the ORIGINAL C.
        for i in 0..n {
            for j in 0..=i {
                let normalized = f64::midpoint(inverse[i * n + j], inverse[j * n + i]);
                let value = finite(normalized / products[i * n + j])?;
                if (i == j && value <= 0.0) || (normalized != 0.0 && value == 0.0) {
                    return Err(invalid("capacitance block inverse is unrepresentable"));
                }
                inverse[i * n + j] = value;
                inverse[j * n + i] = value;
            }
        }
        for i in 0..n {
            for j in 0..n {
                let mut residual = if i == j { -1.0 } else { 0.0 };
                for k in 0..n {
                    let represented = finite(inverse[k * n + j] * products[k * n + j])?;
                    residual = scaled[i * n + k].mul_add(represented, residual);
                }
                if !residual.is_finite() || residual.abs() > 1e-10 {
                    return Err(invalid(
                        "capacitance block inverse is numerically unresolved",
                    ));
                }
            }
        }
        Ok(Self {
            branches: branches.into(),
            matrix_f: matrix_f.into(),
            inverse_per_f: inverse.into(),
        })
    }

    /// Original graph branch indices, in matrix/charge/voltage order.
    pub fn branches(&self) -> &[usize] {
        &self.branches
    }

    /// Original reciprocal C, row-major, in farads.
    pub fn matrix_f(&self) -> &[f64] {
        &self.matrix_f
    }

    /// Admitted represented C inverse, row-major, in inverse farads.
    pub fn inverse_matrix_per_f(&self) -> &[f64] {
        &self.inverse_per_f
    }

    /// Rebind the same prepared constitutive matrix without refactorization.
    /// Graph element types and exclusive storage ownership are checked during
    /// assembly; remapping alone does not change charge/voltage orientation.
    pub fn with_branches(&self, branches: Vec<usize>) -> Result<Self, CircuitError> {
        admit_branches(&branches)?;
        if branches.len() != self.branches.len() {
            return Err(invalid(
                "capacitance remapping must preserve matrix dimension",
            ));
        }
        Ok(Self {
            branches: branches.into(),
            matrix_f: self.matrix_f.clone(),
            inverse_per_f: self.inverse_per_f.clone(),
        })
    }

    /// Electrode/branch charges C v [C], in this block's original order.
    pub fn charges_c(&self, voltage_v: &[f64]) -> Result<Vec<f64>, CircuitError> {
        multiply(&self.matrix_f, voltage_v, self.branches.len())
    }

    /// Actual endpoint constitutive voltages C^-1 q [V].
    pub fn voltages_v(&self, charge_c: &[f64]) -> Result<Vec<f64>, CircuitError> {
        multiply(&self.inverse_per_f, charge_c, self.branches.len())
    }

    /// Complete coupled energy [J], including every cross term.
    pub fn stored_energy_j(&self, charge_c: &[f64]) -> Result<f64, CircuitError> {
        let voltage = self.voltages_v(charge_c)?;
        let (mut sum, mut correction) = (0.0, 0.0);
        for (&q, &v) in charge_c.iter().zip(&voltage) {
            add(&mut sum, &mut correction, finite(q * v)?)?;
        }
        let energy = finite(0.5 * finite(sum + correction)?)?;
        if energy < 0.0 || (energy == 0.0 && charge_c.iter().any(|&q| q != 0.0)) {
            return Err(invalid(
                "nonzero coupled capacitive energy is not representable",
            ));
        }
        Ok(energy)
    }

    /// Stable change of the complete quadratic energy [J].
    /// Uses the work of the mean constitutive voltage against charge change,
    /// including midpoint rounding, instead of subtracting endpoint energies.
    pub fn energy_change_j(&self, before_c: &[f64], after_c: &[f64]) -> Result<f64, CircuitError> {
        let n = self.branches.len();
        admit_vector(before_c, n)?;
        admit_vector(after_c, n)?;
        let (mut sum, mut correction) = (0.0, 0.0);
        for i in 0..n {
            let change = finite(after_c[i] - before_c[i])?;
            if change == 0.0 {
                continue;
            }
            let (mut voltage, mut voltage_correction) = (0.0, 0.0);
            for j in 0..n {
                let midpoint = f64::midpoint(before_c[j], after_c[j]);
                let offset = (before_c[j] - midpoint) + (after_c[j] - midpoint);
                let inverse = self.inverse_per_f[i * n + j];
                add(
                    &mut voltage,
                    &mut voltage_correction,
                    finite(inverse * midpoint)?,
                )?;
                add(
                    &mut voltage,
                    &mut voltage_correction,
                    finite((0.5 * inverse) * offset)?,
                )?;
            }
            add(
                &mut sum,
                &mut correction,
                finite(change * finite(voltage + voltage_correction)?)?,
            )?;
        }
        finite(sum + correction)
    }
}

fn finite(value: f64) -> Result<f64, CircuitError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(
            "coupled capacitive arithmetic is not representable",
        ))
    }
}

fn admit_branches(branches: &[usize]) -> Result<(), CircuitError> {
    if !(1..=MAX_CAPACITANCE_BLOCK_COORDINATES).contains(&branches.len())
        || branches
            .iter()
            .enumerate()
            .any(|(i, &branch)| branch == usize::MAX || branches[..i].contains(&branch))
    {
        return Err(invalid(
            "bounded distinct capacitance block branches required",
        ));
    }
    Ok(())
}

fn admit_vector(values: &[f64], n: usize) -> Result<(), CircuitError> {
    if values.len() != n || values.iter().any(|v| !v.is_finite()) {
        return Err(invalid(
            "complete finite capacitance block coordinates required",
        ));
    }
    Ok(())
}

fn add(sum: &mut f64, correction: &mut f64, value: f64) -> Result<(), CircuitError> {
    let next = finite(*sum + value)?;
    *correction = finite(
        *correction
            + if sum.abs() >= value.abs() {
                (*sum - next) + value
            } else {
                (value - next) + *sum
            },
    )?;
    *sum = next;
    Ok(())
}

fn multiply(matrix: &[f64], values: &[f64], n: usize) -> Result<Vec<f64>, CircuitError> {
    admit_vector(values, n)?;
    let mut result = Vec::with_capacity(n);
    for i in 0..n {
        let (mut sum, mut correction) = (0.0, 0.0);
        for j in 0..n {
            add(
                &mut sum,
                &mut correction,
                finite(matrix[i * n + j] * values[j])?,
            )?;
        }
        result.push(finite(sum + correction)?);
    }
    if result.iter().all(|&v| v == 0.0) && values.iter().any(|&v| v != 0.0) {
        return Err(invalid(
            "nonzero coupled capacitive coordinates disappeared",
        ));
    }
    Ok(result)
}

/// Retain original charge indices and every other quadratic energy block.
pub(super) fn stamp(
    graph: &CircuitGraph,
    dae: &CircuitDae,
    blocks: &[CapacitanceBlock],
    curves: &[&CapacitorCurve],
    matrix: &mut [f64],
) -> Result<(), CircuitError> {
    if blocks.is_empty() {
        return Ok(());
    }
    let mut coordinates = vec![None; graph.branches.len()];
    let mut charges = dae.charge_index.iter();
    for (branch, &(_, _, element)) in graph.branches.iter().enumerate() {
        if matches!(element, Branch::Capacitor { .. }) {
            coordinates[branch] = charges.next().copied();
        }
    }
    let mut owned = vec![false; graph.branches.len()];
    for curve in curves {
        owned[curve.branch] = true;
    }
    let n = dae.system.state_dim();
    for block in blocks {
        let mut indices = Vec::with_capacity(block.branches.len());
        for &branch in block.branches.iter() {
            let Some(Some(coordinate)) = coordinates.get(branch) else {
                return Err(invalid(
                    "capacitance block must bind existing capacitor branches",
                ));
            };
            if owned[branch] {
                return Err(invalid(
                    "capacitor block and curve storage owners must not overlap",
                ));
            }
            owned[branch] = true;
            indices.push(*coordinate);
        }
        let count = indices.len();
        for (i, &row) in indices.iter().enumerate() {
            for (j, &column) in indices.iter().enumerate() {
                matrix[row * n + column] = block.inverse_per_f[i * count + j];
            }
        }
    }
    Ok(())
}
