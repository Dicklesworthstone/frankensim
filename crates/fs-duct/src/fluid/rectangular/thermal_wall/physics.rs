//! Original Q1 diffusion form on the physical rectangle.
//!
//! K + z M + sum(beta B) is complex symmetric. Its exact real lowering
//! [Re A, -Im A; -Im A, -Re A] uses the existing bounded sparse LDLT.
//! No dense matrix, iterative tolerance surrogate, or replacement thermal law.
use super::{
    RectangularThermalLimits, RectangularThermalWall as Wall, RectangularThermalWallError as Error,
    bad, finite, poll,
};
use fs_math::{c64::C64, dd::Dd, det};
use fs_sparse::{Coo, DirectOrdering, LdltOptions, SymbolicLdlt};
use std::cell::RefCell;

pub(super) const MAX_NODES: usize = 16_384;

#[derive(Debug)]
pub(super) struct ThermalField {
    pub dimensions: [f64; 2],
    pub cells: [usize; 2],
    stored: Vec<C64>,
    complement: bool,
    pub mean: C64,
    pub gradient_form: f64,
    pub field_form: f64,
    pub wall_forms: [f64; 4],
    pub relative_residual: f64,
    pub peak_payload_bytes: usize,
    pub factor_nonzeros: usize,
}

pub(super) fn prepare(
    dimensions: [f64; 2],
    z: C64,
    beta: [C64; 4],
    walls: [Wall; 4],
    limits: RectangularThermalLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<ThermalField, Error> {
    poll(cancelled, 0)?;
    let [nx, ny] = limits.cells;
    let nodes = nx
        .checked_add(1)
        .and_then(|x| ny.checked_add(1).and_then(|y| x.checked_mul(y)))
        .ok_or_else(|| bad("rectangular thermal mesh extent overflow"))?;
    if nx < 2
        || ny < 2
        || nodes > MAX_NODES
        || !dimensions.into_iter().all(|v| v.is_finite() && v > 0.0)
        || !(dimensions[0] * dimensions[1]).is_finite()
        || dimensions[0] * dimensions[1] <= 0.0
        || !finite(z)
        || z.re < 0.0
        || z.im <= 0.0
        || !beta
            .into_iter()
            .all(|v| finite(v) && v.re >= 0.0 && v.im >= 0.0)
        || !limits.maximum_relative_residual.is_finite()
        || !(0.0..1.0).contains(&limits.maximum_relative_residual)
        || limits.maximum_relative_residual == 0.0
        || !limits.maximum_cell_diffusion_phase.is_finite()
        || limits.maximum_cell_diffusion_phase <= 0.0
        || limits.maximum_cell_diffusion_phase > 1.0
    {
        return Err(bad(
            "rectangular thermal geometry, coefficients or numerical limits are invalid",
        )
        .into());
    }
    let h = [dimensions[0] / nx as f64, dimensions[1] / ny as f64];
    if !h.into_iter().all(|v| v.is_finite() && v > 0.0)
        || !h[0].mul_add(h[1], 0.0).is_finite()
        || h[0].max(h[1]) * det::sqrt(z.abs()) > limits.maximum_cell_diffusion_phase
    {
        return Err(
            bad("rectangular thermal mesh does not resolve its physical diffusion phase").into(),
        );
    }
    // COO staging plus canonical assembly, retained fields/rows and candidate
    // vectors. Native symbolic/factor storage is separately admitted below.
    let staged = nx
        .checked_mul(ny)
        .and_then(|v| v.checked_mul(64))
        .and_then(|v| (nx + ny).checked_mul(32).and_then(|b| v.checked_add(b)))
        .ok_or_else(|| bad("rectangular thermal staging extent overflow"))?;
    let fixed = staged
        .checked_mul(48)
        .and_then(|v| nodes.checked_mul(512).and_then(|n| v.checked_add(n)))
        .and_then(|v| v.checked_add(8192))
        .ok_or_else(|| bad("rectangular thermal payload extent overflow"))?;
    if fixed >= limits.maximum_payload_bytes {
        return Err(bad("rectangular thermal staging exceeds the declared payload").into());
    }
    let insulated = walls.into_iter().all(|w| w != Wall::Isothermal) && beta == [C64::ZERO; 4];
    let complement = insulated || z.abs() * dimensions[0].max(dimensions[1]).powi(2) >= 1.0;
    let mut stored = filled(nodes, if complement { C64::ONE } else { C64::ZERO })?;
    let mut map = filled(nodes, usize::MAX)?;
    let mut free = 0;
    for y in 0..=ny {
        poll(cancelled, y)?;
        for x in 0..=nx {
            let fixed_node = (x == 0 && walls[0] == Wall::Isothermal)
                || (x == nx && walls[1] == Wall::Isothermal)
                || (y == 0 && walls[2] == Wall::Isothermal)
                || (y == ny && walls[3] == Wall::Isothermal);
            if !fixed_node {
                map[y * (nx + 1) + x] = free;
                free += 1;
            }
        }
    }
    let (mut residual, mut peak, mut nnz) = (0.0_f64, fixed, 0);
    if insulated {
        stored.fill(C64::ZERO);
    } else {
        let mut coo = Coo::new(2 * free, 2 * free);
        let mut rhs = filled(2 * free, 0.0)?;
        for y in 0..ny {
            poll(cancelled, y)?;
            for x in 0..nx {
                let nodes = cell_nodes(nx, x, y);
                for a in 0..4 {
                    let row = map[nodes[a]];
                    if row != usize::MAX && !complement {
                        add_rhs(&mut rhs, row, z.scale(h[0] * h[1] * 0.25));
                    }
                    for b in 0..4 {
                        let (m, k) = element(h, a, b);
                        add_entry(
                            &mut coo,
                            &mut rhs,
                            &map,
                            nodes[a],
                            nodes[b],
                            C64::from_re(k) + z.scale(m),
                            complement,
                        );
                    }
                }
            }
        }
        for side in 0..4 {
            let (count, length) = if side < 2 { (ny, h[1]) } else { (nx, h[0]) };
            for edge in 0..count {
                if edge % 32 == 0 {
                    poll(cancelled, edge)?;
                }
                let pair = edge_nodes([nx, ny], side, edge);
                for a in 0..2 {
                    let row = map[pair[a]];
                    if row != usize::MAX && complement {
                        add_rhs(&mut rhs, row, beta[side].scale(length * 0.5));
                    }
                    for b in 0..2 {
                        add_entry(
                            &mut coo,
                            &mut rhs,
                            &map,
                            pair[a],
                            pair[b],
                            beta[side].scale(length * if a == b { 1.0 / 3.0 } else { 1.0 / 6.0 }),
                            complement,
                        );
                    }
                }
            }
        }
        poll(cancelled, staged)?;
        let original = coo.assemble();
        drop(coo);
        let solved = solve_assembled(&original, &rhs, free, limits, fixed, cancelled)?;
        residual = solved.relative_residual;
        peak = solved.peak_payload_bytes;
        nnz = solved.factor_nonzeros;
        let answer = solved.answer;
        for (node, &index) in map.iter().enumerate() {
            if index != usize::MAX {
                stored[node] = C64::new(answer[2 * index], answer[2 * index + 1]);
            }
        }
    }
    let mut result = ThermalField {
        dimensions,
        cells: [nx, ny],
        stored,
        complement,
        mean: C64::ZERO,
        gradient_form: 0.0,
        field_form: 0.0,
        wall_forms: [0.0; 4],
        relative_residual: residual,
        peak_payload_bytes: peak,
        factor_nonzeros: nnz,
    };
    result.integrate(cancelled)?;
    poll(cancelled, nodes)?;
    Ok(result)
}

/// Same scaled native solve for gas-only and conjugate gas/solid Q1 owners.
pub(super) struct SparseThermalSolution {
    pub answer: Vec<f64>,
    pub relative_residual: f64,
    pub peak_payload_bytes: usize,
    pub factor_nonzeros: usize,
}
pub(super) fn solve_assembled(
    original: &fs_sparse::Csr,
    rhs: &[f64],
    free: usize,
    limits: RectangularThermalLimits,
    fixed: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<SparseThermalSolution, Error> {
    let (mut residual, mut peak) = (0.0_f64, fixed);
    if (0..original.nrows()).any(|row| original.row(row).1.iter().any(|v| !v.is_finite()))
        || rhs.iter().any(|v| !v.is_finite())
    {
        return Err(bad("rectangular thermal matrix or forcing is unrepresentable").into());
    }
    let mut scaled = original.clone();
    let mut scale = filled(2 * free, 1.0)?;
    for i in 0..free {
        let diagonal = det::hypot(original.get(2 * i, 2 * i), original.get(2 * i, 2 * i + 1));
        if !diagonal.is_finite() || diagonal <= 0.0 {
            return Err(bad("rectangular thermal diagonal scale is unrepresentable").into());
        }
        let s = 1.0 / det::sqrt(diagonal);
        scale[2 * i] = s;
        scale[2 * i + 1] = s;
    }
    for i in 0..2 * free {
        if i % 64 == 0 {
            poll(cancelled, i)?;
        }
        let columns = scaled.row(i).0.to_vec();
        for (&j, value) in columns.iter().zip(scaled.row_values_mut(i)) {
            // Identical factor ordering preserves exact numeric symmetry.
            let (a, b) = if i <= j { (i, j) } else { (j, i) };
            *value = (*value * scale[a]) * scale[b];
        }
    }
    let mut direct = limits.sparse;
    direct.max_symbolic_entries = direct
        .max_symbolic_entries
        .min((limits.maximum_payload_bytes - fixed) / 64);
    let gate = RefCell::new(&mut *cancelled);
    let checkpoint = || (*gate.borrow_mut())();
    let symbolic =
        SymbolicLdlt::analyze_bounded(&scaled, DirectOrdering::Amd, &direct, &checkpoint)?;
    let symbolic_bytes = symbolic.symbolic_bytes();
    let remaining = limits
        .maximum_payload_bytes
        .checked_sub(fixed)
        .and_then(|v| v.checked_sub(symbolic_bytes))
        .ok_or_else(|| bad("rectangular thermal symbolic payload exceeds admission"))?;
    let factor = symbolic.factor_bounded(
        &scaled,
        &LdltOptions {
            max_factor_bytes: Some(remaining),
            ..LdltOptions::default()
        },
        &direct,
        &checkpoint,
    )?;
    if factor.inertia().positive != free || factor.inertia().negative != free {
        return Err(bad("rectangular complex thermal factor has inconsistent inertia").into());
    }
    let nnz = factor.stats().nnz_l;
    peak += symbolic_bytes + factor.stats().peak_factor_bytes;
    let scaled_rhs: Vec<_> = rhs.iter().zip(&scale).map(|(r, s)| r * s).collect();
    if checkpoint() {
        return Err(fs_sparse::LdltError::Cancelled.into());
    }
    let mut answer = factor.solve(&scaled_rhs);
    if checkpoint() {
        return Err(fs_sparse::LdltError::Cancelled.into());
    }
    for (v, s) in answer.iter_mut().zip(&scale) {
        *v *= *s;
    }
    for i in 0..2 * free {
        if i % 64 == 0 && checkpoint() {
            return Err(fs_sparse::LdltError::Cancelled.into());
        }
        let (columns, entries) = original.row(i);
        let mut difference = -Dd::from_f64(rhs[i]);
        let mut row_scale = Dd::from_f64(rhs[i].abs());
        for (&j, &a) in columns.iter().zip(entries) {
            difference = difference + Dd::from_f64(a) * Dd::from_f64(answer[j]);
            row_scale = row_scale + Dd::from_f64(a.abs()) * Dd::from_f64(answer[j].abs());
        }
        residual =
            residual.max(difference.to_f64().abs() / row_scale.to_f64().max(f64::MIN_POSITIVE));
    }
    if !residual.is_finite()
        || residual > limits.maximum_relative_residual
        || answer.iter().any(|v| !v.is_finite())
    {
        return Err(bad("rectangular thermal original-equation residual exceeds admission").into());
    }

    Ok(SparseThermalSolution {
        answer,
        relative_residual: residual,
        peak_payload_bytes: peak,
        factor_nonzeros: nnz,
    })
}
fn add_rhs(rhs: &mut [f64], row: usize, value: C64) {
    rhs[2 * row] += value.re;
    rhs[2 * row + 1] -= value.im;
}
fn add_entry(
    coo: &mut Coo,
    rhs: &mut [f64],
    map: &[usize],
    a: usize,
    b: usize,
    v: C64,
    complement: bool,
) {
    let row = map[a];
    let column = map[b];
    if row == usize::MAX {
        return;
    }
    if column == usize::MAX {
        if complement {
            add_rhs(rhs, row, v.scale(-1.0));
        }
    } else {
        coo.push(2 * row, 2 * column, v.re);
        coo.push(2 * row, 2 * column + 1, -v.im);
        coo.push(2 * row + 1, 2 * column, -v.im);
        coo.push(2 * row + 1, 2 * column + 1, -v.re);
    }
}
fn cell_nodes(nx: usize, x: usize, y: usize) -> [usize; 4] {
    let a = y * (nx + 1) + x;
    [a, a + 1, a + nx + 1, a + nx + 2]
}
fn edge_nodes([nx, ny]: [usize; 2], side: usize, edge: usize) -> [usize; 2] {
    match side {
        0 => [edge * (nx + 1), (edge + 1) * (nx + 1)],
        1 => [edge * (nx + 1) + nx, (edge + 1) * (nx + 1) + nx],
        2 => [edge, edge + 1],
        _ => [ny * (nx + 1) + edge, ny * (nx + 1) + edge + 1],
    }
}
fn element(h: [f64; 2], a: usize, b: usize) -> (f64, f64) {
    let mx = if (a & 1) == (b & 1) {
        1.0 / 3.0
    } else {
        1.0 / 6.0
    };
    let my = if (a >> 1) == (b >> 1) {
        1.0 / 3.0
    } else {
        1.0 / 6.0
    };
    let kx = if (a & 1) == (b & 1) { 1.0 } else { -1.0 };
    let ky = if (a >> 1) == (b >> 1) { 1.0 } else { -1.0 };
    (
        h[0] * h[1] * mx * my,
        h[1] / h[0] * kx * my + h[0] / h[1] * mx * ky,
    )
}
impl ThermalField {
    pub fn at(&self, point: [f64; 2]) -> Result<[C64; 3], Error> {
        if point
            .iter()
            .zip(self.dimensions)
            .any(|(x, l)| !x.is_finite() || *x < 0.0 || *x > l)
        {
            return Err(bad("rectangular thermal point lies outside its physical section").into());
        }
        let mut cell = [0; 2];
        let mut fraction = [0.0; 2];
        for a in 0..2 {
            let q = (point[a] / self.dimensions[a]) * self.cells[a] as f64;
            cell[a] = (q.floor() as usize).min(self.cells[a] - 1);
            fraction[a] = (q - cell[a] as f64).min(1.0);
        }
        Ok(self.in_cell(cell, fraction))
    }
    fn in_cell(&self, cell: [usize; 2], fraction: [f64; 2]) -> [C64; 3] {
        let h = [
            self.dimensions[0] / self.cells[0] as f64,
            self.dimensions[1] / self.cells[1] as f64,
        ];
        let ids = cell_nodes(self.cells[0], cell[0], cell[1]);
        let mut re = [Dd::ZERO; 3];
        let mut im = [Dd::ZERO; 3];
        for a in 0..4 {
            let x = if a & 1 == 0 {
                1.0 - fraction[0]
            } else {
                fraction[0]
            };
            let y = if a >> 1 == 0 {
                1.0 - fraction[1]
            } else {
                fraction[1]
            };
            let dx = if a & 1 == 0 { -1.0 / h[0] } else { 1.0 / h[0] };
            let dy = if a >> 1 == 0 { -1.0 / h[1] } else { 1.0 / h[1] };
            let shape = [x * y, dx * y, x * dy];
            for d in 0..3 {
                re[d] = re[d] + Dd::from_f64(shape[d]) * Dd::from_f64(self.stored[ids[a]].re);
                im[d] = im[d] + Dd::from_f64(shape[d]) * Dd::from_f64(self.stored[ids[a]].im);
            }
        }
        let mut result = core::array::from_fn(|i| C64::new(re[i].to_f64(), im[i].to_f64()));
        if self.complement {
            result[0] = C64::new((Dd::ONE - re[0]).to_f64(), -im[0].to_f64());
            result[1] = result[1].scale(-1.0);
            result[2] = result[2].scale(-1.0);
        }
        result
    }
    fn integrate(&mut self, cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
        let gauss = [0.5 - 0.5 / det::sqrt(3.0), 0.5 + 0.5 / det::sqrt(3.0)];
        let mut mean = [Dd::ZERO; 2];
        let mut gradient = Dd::ZERO;
        let mut field = Dd::ZERO;
        let h = [
            self.dimensions[0] / self.cells[0] as f64,
            self.dimensions[1] / self.cells[1] as f64,
        ];
        let weight = h[0] * h[1] * 0.25;
        for y in 0..self.cells[1] {
            poll(cancelled, y)?;
            for x in 0..self.cells[0] {
                for a in gauss {
                    for b in gauss {
                        let f = self.in_cell([x, y], [a, b]);
                        mean[0] = mean[0] + Dd::from_f64(weight) * Dd::from_f64(f[0].re);
                        mean[1] = mean[1] + Dd::from_f64(weight) * Dd::from_f64(f[0].im);
                        field = field + Dd::from_f64(weight) * Dd::from_f64(f[0].norm_sq());
                        gradient = gradient
                            + Dd::from_f64(weight)
                                * (Dd::from_f64(f[1].norm_sq()) + Dd::from_f64(f[2].norm_sq()));
                    }
                }
            }
        }
        let area = self.dimensions[0] * self.dimensions[1];
        self.mean = C64::new(
            (mean[0] / Dd::from_f64(area)).to_f64(),
            (mean[1] / Dd::from_f64(area)).to_f64(),
        );
        self.gradient_form = gradient.to_f64();
        self.field_form = field.to_f64();
        for side in 0..4 {
            let (count, length) = if side < 2 {
                (self.cells[1], h[1])
            } else {
                (self.cells[0], h[0])
            };
            let mut value = Dd::ZERO;
            for edge in 0..count {
                if edge % 32 == 0 {
                    poll(cancelled, edge)?;
                }
                for t in gauss {
                    let point = match side {
                        0 => [0.0, (edge as f64 + t) * h[1]],
                        1 => [self.dimensions[0], (edge as f64 + t) * h[1]],
                        2 => [(edge as f64 + t) * h[0], 0.0],
                        _ => [(edge as f64 + t) * h[0], self.dimensions[1]],
                    };
                    value = value
                        + Dd::from_f64(0.5 * length) * Dd::from_f64(self.at(point)?[0].norm_sq());
                }
            }
            self.wall_forms[side] = value.to_f64();
        }
        if !finite(self.mean)
            || [self.gradient_form, self.field_form]
                .into_iter()
                .chain(self.wall_forms)
                .any(|v| !v.is_finite() || v < 0.0)
        {
            return Err(bad("rectangular thermal physical integrals are unrepresentable").into());
        }
        Ok(())
    }
}
fn filled<T: Clone>(count: usize, value: T) -> Result<Vec<T>, Error> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(count)
        .map_err(|_| bad("rectangular thermal allocation failed"))?;
    result.resize(count, value);
    Ok(result)
}
#[cfg(test)]
mod tests;
