//! Conforming gas and nested solid Q1 cells, with the original full tensor K.
use super::super::{
    RectangularThermalLimits, RectangularThermalWallError as Error, bad, finite, poll,
};
use super::RectangularJacketBath as Bath;
use fs_math::{c64::C64, dd::Dd, det};
use fs_sparse::Coo;

#[derive(Clone, Copy, Debug)]
pub(super) struct Layer {
    pub thickness: [f64; 4],
    pub cells: [usize; 4],
    pub conductivity: [[f64; 2]; 2],
    pub capacity: f64,
}
#[derive(Clone, Copy, Debug, Default)]
pub(super) struct Forms {
    pub temperature_integral: C64,
    pub field: f64,
    pub conduction: f64,
}
#[derive(Debug)]
pub(super) struct JacketField {
    pub axes: [Vec<f64>; 2],
    regions: [Vec<usize>; 2],
    dimensions: [f64; 2],
    stored: Vec<C64>,
    lift: Option<C64>,
    pub forms: Vec<Forms>,
    pub bath_forms: [f64; 4],
    pub bath_flux: C64,
    pub mean_gas: C64,
    pub relative_residual: f64,
    pub relative_heat_defect: f64,
    pub relative_work_defect: f64,
    pub peak_payload_bytes: usize,
    pub factor_nonzeros: usize,
}

#[allow(clippy::too_many_arguments)]
pub(super) fn prepare(
    dimensions: [f64; 2],
    conductivity: f64,
    capacity: C64,
    omega: f64,
    layers: &[Layer],
    baths: [Bath; 4],
    limits: RectangularThermalLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<JacketField, Error> {
    poll(cancelled, 0)?;
    if layers.is_empty()
        || layers.len() > 8
        || !dimensions.into_iter().all(positive)
        || !positive(dimensions[0] * dimensions[1])
        || !positive(conductivity)
        || !positive(omega)
        || !finite(capacity)
        || capacity.re <= 0.0
        || capacity.im > 0.0
        || limits.cells.into_iter().any(|v| v < 2)
        || !positive(limits.maximum_cell_diffusion_phase)
        || limits.maximum_cell_diffusion_phase > 1.0
        || !positive(limits.maximum_relative_residual)
        || limits.maximum_relative_residual >= 1.0
    {
        return Err(bad("invalid original rectangular gas/jacket coefficients or limits").into());
    }
    let mut k = vec![[[conductivity, 0.0], [0.0, conductivity]]];
    let mut mass = vec![capacity];
    let mut diffusion = vec![capacity.abs() / conductivity];
    for (index, layer) in layers.iter().enumerate() {
        poll(cancelled, index)?;
        let a = layer.conductivity;
        let (lower, upper) = principal_bounds(a)?;
        if !layer.thickness.into_iter().all(positive)
            || layer.cells.into_iter().any(|v| v == 0 || v > 4096)
            || !positive(layer.capacity)
            || !positive(lower)
            || !positive(upper)
        {
            return Err(bad("positive resolved source-backed rectangular layers required").into());
        }
        k.push(a);
        mass.push(C64::from_re(layer.capacity));
        diffusion.push(layer.capacity / lower);
    }
    let conductance = [
        baths[0].conductance()?,
        baths[1].conductance()?,
        baths[2].conductance()?,
        baths[3].conductance()?,
    ];
    let nx = axis_cells(limits.cells[0], 0, layers)?;
    let ny = axis_cells(limits.cells[1], 1, layers)?;
    let nodes = nx
        .checked_add(1)
        .and_then(|x| ny.checked_add(1).and_then(|y| x.checked_mul(y)))
        .ok_or_else(|| bad("rectangular jacket node extent overflow"))?;
    if nodes > super::super::physics::MAX_NODES {
        return Err(bad("rectangular gas/jacket exceeds 16384 physical nodes").into());
    }
    let staged = nx
        .checked_mul(ny)
        .and_then(|v| v.checked_mul(64))
        .and_then(|v| (nx + ny).checked_mul(32).and_then(|b| v.checked_add(b)))
        .ok_or_else(|| bad("rectangular jacket staging extent overflow"))?;
    let fixed = staged
        .checked_mul(48)
        .and_then(|v| nodes.checked_mul(640).and_then(|n| v.checked_add(n)))
        .and_then(|v| v.checked_add(16384))
        .ok_or_else(|| bad("rectangular jacket payload extent overflow"))?;
    if fixed >= limits.maximum_payload_bytes {
        return Err(bad("rectangular jacket staging exceeds declared payload").into());
    }
    let (x, rx) = axis(dimensions[0], limits.cells[0], 0, layers)?;
    let (y, ry) = axis(dimensions[1], limits.cells[1], 1, layers)?;
    let z = C64::new(-capacity.im, capacity.re).scale(omega / conductivity);
    let lift =
        (z.abs() * dimensions[0].max(dimensions[1]).powi(2) >= 1.0).then(|| capacity.recip());
    let mut map = filled(nodes, usize::MAX)?;
    let mut free = 0usize;
    for j in 0..=ny {
        poll(cancelled, j)?;
        for i in 0..=nx {
            let prescribed = (i == 0 && baths[0] == Bath::Isothermal)
                || (i == nx && baths[1] == Bath::Isothermal)
                || (j == 0 && baths[2] == Bath::Isothermal)
                || (j == ny && baths[3] == Bath::Isothermal);
            if !prescribed {
                map[j * (nx + 1) + i] = free;
                free += 1;
            }
        }
    }
    let mut field = JacketField {
        axes: [x, y],
        regions: [rx, ry],
        dimensions,
        stored: filled(nodes, C64::ZERO)?,
        lift,
        forms: vec![Forms::default(); layers.len() + 1],
        bath_forms: [0.0; 4],
        bath_flux: C64::ZERO,
        mean_gas: C64::ZERO,
        relative_residual: 0.0,
        relative_heat_defect: 0.0,
        relative_work_defect: 0.0,
        peak_payload_bytes: fixed,
        factor_nonzeros: 0,
    };
    let mut coo = Coo::new(2 * free, 2 * free);
    let mut rhs = filled(2 * free, 0.0)?;
    for j in 0..ny {
        poll(cancelled, j)?;
        for i in 0..nx {
            let region = field.region([i, j]);
            let h = field.width([i, j]);
            let phase = h[0].max(h[1]) * det::sqrt(omega * diffusion[region]);
            if !phase.is_finite() || phase > limits.maximum_cell_diffusion_phase {
                return Err(
                    bad("original gas or solid cell does not resolve thermal diffusion").into(),
                );
            }
            let ids = cell_nodes(nx, i, j);
            let a = element(h, k[region], mass[region], omega);
            for row in 0..4 {
                if map[ids[row]] == usize::MAX {
                    continue;
                }
                let forcing = if let Some(lift) = field.lift {
                    if region == 0 {
                        C64::ZERO
                    } else {
                        let mut sum = Sum::default();
                        for column in 0..4 {
                            if field.gas_node(ids[column]) {
                                sum.add(a[row][column] * lift);
                            }
                        }
                        sum.value()
                    }
                } else if region == 0 {
                    C64::new(0.0, omega * h[0] * h[1] * 0.25)
                } else {
                    C64::ZERO
                };
                add_rhs(&mut rhs, map[ids[row]], forcing);
                for column in 0..4 {
                    add(&mut coo, &map, ids[row], ids[column], a[row][column]);
                }
            }
        }
    }
    for side in 0..4 {
        for edge in 0..field.edge_count(side) {
            if edge % 32 == 0 {
                poll(cancelled, edge)?;
            }
            let (ids, length) = field.edge(side, edge);
            for a in 0..2 {
                for b in 0..2 {
                    add(
                        &mut coo,
                        &map,
                        ids[a],
                        ids[b],
                        C64::from_re(
                            conductance[side] * length * if a == b { 1.0 / 3.0 } else { 1.0 / 6.0 },
                        ),
                    );
                }
            }
        }
    }
    poll(cancelled, staged)?;
    let original = coo.assemble();
    drop(coo);
    let solved =
        super::super::physics::solve_assembled(&original, &rhs, free, limits, fixed, cancelled)?;
    for (node, &index) in map.iter().enumerate() {
        if index != usize::MAX {
            field.stored[node] = C64::new(solved.answer[2 * index], solved.answer[2 * index + 1]);
        }
    }
    field.relative_residual = solved.relative_residual;
    field.peak_payload_bytes = solved.peak_payload_bytes;
    field.factor_nonzeros = solved.factor_nonzeros;
    field.integrate(&k, &mass, omega, conductance, &map, cancelled)?;
    if !field.relative_heat_defect.is_finite()
        || !field.relative_work_defect.is_finite()
        || field.relative_heat_defect > limits.maximum_relative_residual
        || field.relative_work_defect > limits.maximum_relative_residual
    {
        return Err(bad(
            "original gas/solid heat or acoustic availability balance exceeds admission",
        )
        .into());
    }
    poll(cancelled, nodes)?;
    Ok(field)
}

fn axis_cells(gas_cells: usize, axis: usize, layers: &[Layer]) -> Result<usize, Error> {
    let count = layers.iter().try_fold(gas_cells, |n, layer| {
        n.checked_add(layer.cells[2 * axis])
            .and_then(|v| v.checked_add(layer.cells[2 * axis + 1]))
            .ok_or_else(|| bad("rectangular layered axis extent overflow"))
    })?;
    if count >= super::super::physics::MAX_NODES {
        return Err(bad("rectangular layered axis exceeds its physical extent").into());
    }
    Ok(count)
}
fn axis(
    length: f64,
    gas_cells: usize,
    axis: usize,
    layers: &[Layer],
) -> Result<(Vec<f64>, Vec<usize>), Error> {
    let mut left = vec![0.0];
    let mut right = vec![length];
    let count = axis_cells(gas_cells, axis, layers)?;
    for layer in layers {
        left.push(left.last().copied().unwrap() - layer.thickness[2 * axis]);
        right.push(right.last().copied().unwrap() + layer.thickness[2 * axis + 1]);
    }
    if count >= super::super::physics::MAX_NODES
        || !left.iter().chain(&right).all(|v| v.is_finite())
    {
        return Err(bad("rectangular layered axis exceeds its physical extent").into());
    }
    let mut coordinates = Vec::with_capacity(count + 1);
    let mut regions = Vec::with_capacity(count);
    coordinates.push(*left.last().unwrap());
    let mut segment = |end: f64, cells: usize, region: usize| {
        let start = *coordinates.last().unwrap();
        for i in 1..=cells {
            coordinates.push(if i == cells {
                end
            } else {
                start + (end - start) * (i as f64 / cells as f64)
            });
            regions.push(region);
        }
    };
    for index in (0..layers.len()).rev() {
        segment(left[index], layers[index].cells[2 * axis], index + 1);
    }
    segment(length, gas_cells, 0);
    for (index, layer) in layers.iter().enumerate() {
        segment(right[index + 1], layer.cells[2 * axis + 1], index + 1);
    }
    if coordinates.windows(2).any(|w| !positive(w[1] - w[0])) {
        return Err(bad("original rectangular layer cell width is unrepresentable").into());
    }
    Ok((coordinates, regions))
}

impl JacketField {
    fn region(&self, cell: [usize; 2]) -> usize {
        self.regions[0][cell[0]].max(self.regions[1][cell[1]])
    }
    fn width(&self, cell: [usize; 2]) -> [f64; 2] {
        core::array::from_fn(|d| self.axes[d][cell[d] + 1] - self.axes[d][cell[d]])
    }
    fn gas_node(&self, node: usize) -> bool {
        let x = self.axes[0][node % self.axes[0].len()];
        let y = self.axes[1][node / self.axes[0].len()];
        x >= 0.0 && x <= self.dimensions[0] && y >= 0.0 && y <= self.dimensions[1]
    }
    fn node(&self, node: usize) -> C64 {
        self.lift.map_or(self.stored[node], |lift| {
            if self.gas_node(node) {
                lift - self.stored[node]
            } else {
                -self.stored[node]
            }
        })
    }
    pub fn at(&self, region: usize, point: [f64; 2]) -> Result<[C64; 3], Error> {
        if region >= self.forms.len() || point.into_iter().any(|x| !x.is_finite()) {
            return Err(bad("original gas or solid region and finite point required").into());
        }
        let mut candidates = [[0usize; 2]; 2];
        for d in 0..2 {
            let axis = &self.axes[d];
            if point[d] < axis[0] || point[d] > axis[axis.len() - 1] {
                return Err(
                    bad("rectangular conjugate point lies outside its source geometry").into(),
                );
            }
            let right = axis
                .partition_point(|x| *x <= point[d])
                .saturating_sub(1)
                .min(axis.len() - 2);
            candidates[d] = [
                right,
                if point[d] == axis[right] {
                    right.saturating_sub(1)
                } else {
                    right
                },
            ];
        }
        for &y in &candidates[1] {
            for &x in &candidates[0] {
                if self.region([x, y]) == region {
                    let h = self.width([x, y]);
                    let fraction = [
                        (point[0] - self.axes[0][x]) / h[0],
                        (point[1] - self.axes[1][y]) / h[1],
                    ];
                    return Ok(self.in_cell([x, y], fraction));
                }
            }
        }
        Err(bad("observation must select the original side of a gas/solid interface").into())
    }
    fn in_cell(&self, cell: [usize; 2], fraction: [f64; 2]) -> [C64; 3] {
        let h = self.width(cell);
        let ids = cell_nodes(self.axes[0].len() - 1, cell[0], cell[1]);
        let mut sums = [Sum::default(); 3];
        let mut lifts = [Sum::default(); 3];
        for a in 0..4 {
            let n = shape(a, fraction, h);
            for d in 0..3 {
                sums[d].add(self.stored[ids[a]].scale(n[d]));
                if let Some(lift) = self.lift {
                    if self.gas_node(ids[a]) {
                        lifts[d].add(lift.scale(n[d]));
                    }
                }
            }
        }
        if let Some(lift) = self.lift {
            if self.region(cell) == 0 {
                [lift - sums[0].value(), -sums[1].value(), -sums[2].value()]
            } else {
                core::array::from_fn(|d| lifts[d].value() - sums[d].value())
            }
        } else {
            sums.map(Sum::value)
        }
    }
    fn edge_count(&self, side: usize) -> usize {
        self.axes[usize::from(side < 2)].len() - 1
    }
    fn edge(&self, side: usize, edge: usize) -> ([usize; 2], f64) {
        let nx = self.axes[0].len() - 1;
        let ny = self.axes[1].len() - 1;
        let ids = match side {
            0 => [edge * (nx + 1), (edge + 1) * (nx + 1)],
            1 => [edge * (nx + 1) + nx, (edge + 1) * (nx + 1) + nx],
            2 => [edge, edge + 1],
            _ => [ny * (nx + 1) + edge, ny * (nx + 1) + edge + 1],
        };
        let axis = &self.axes[usize::from(side < 2)];
        (ids, axis[edge + 1] - axis[edge])
    }
    fn integrate(
        &mut self,
        k: &[[[f64; 2]; 2]],
        capacity: &[C64],
        omega: f64,
        bath: [f64; 4],
        map: &[usize],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        let gauss = [0.5 - 0.5 / det::sqrt(3.0), 0.5 + 0.5 / det::sqrt(3.0)];
        let mut integral = vec![Sum::default(); self.forms.len()];
        let mut norm = vec![Dd::ZERO; self.forms.len()];
        let mut heat = norm.clone();
        let mut reaction = Sum::default();
        let nx = self.axes[0].len() - 1;
        let ny = self.axes[1].len() - 1;
        for j in 0..ny {
            poll(cancelled, j)?;
            for i in 0..nx {
                let h = self.width([i, j]);
                let region = self.region([i, j]);
                let weight = h[0] * h[1] * 0.25;
                for x in gauss {
                    for y in gauss {
                        let f = self.in_cell([i, j], [x, y]);
                        integral[region].add(f[0].scale(weight));
                        norm[region] =
                            norm[region] + Dd::from_f64(weight) * Dd::from_f64(f[0].norm_sq());
                        heat[region] = heat[region]
                            + Dd::from_f64(weight)
                                * Dd::from_f64(quadratic(k[region], [f[1], f[2]])?);
                    }
                }
                let ids = cell_nodes(nx, i, j);
                let matrix = element(h, k[region], capacity[region], omega);
                for a in 0..4 {
                    if map[ids[a]] == usize::MAX {
                        for b in 0..4 {
                            reaction.add(matrix[a][b] * self.node(ids[b]));
                        }
                        if region == 0 {
                            reaction.add(C64::new(0.0, -omega * weight));
                        }
                    }
                }
            }
        }
        let mut flux = Sum::default();
        for side in 0..4 {
            let mut form = Dd::ZERO;
            for edge in 0..self.edge_count(side) {
                if edge % 32 == 0 {
                    poll(cancelled, edge)?;
                }
                let (ids, length) = self.edge(side, edge);
                let t = ids.map(|i| self.node(i));
                for s in gauss {
                    let value = t[0].scale(1.0 - s) + t[1].scale(s);
                    form = form + Dd::from_f64(0.5 * length) * Dd::from_f64(value.norm_sq());
                    flux.add(value.scale(bath[side] * length * 0.5));
                }
                for a in 0..2 {
                    if map[ids[a]] == usize::MAX {
                        for b in 0..2 {
                            reaction.add(t[b].scale(
                                bath[side] * length * if a == b { 1.0 / 3.0 } else { 1.0 / 6.0 },
                            ));
                        }
                    }
                }
            }
            self.bath_forms[side] = form.to_f64();
        }
        for region in 0..self.forms.len() {
            self.forms[region] = Forms {
                temperature_integral: integral[region].value(),
                field: norm[region].to_f64(),
                conduction: heat[region].to_f64(),
            };
        }
        self.bath_flux = flux.value() - reaction.value();
        let area = self.dimensions[0] * self.dimensions[1];
        self.mean_gas = self.forms[0].temperature_integral.scale(1.0 / area);
        let mut energy = Sum::default();
        for (form, c) in self.forms.iter().zip(capacity) {
            energy.add(form.temperature_integral * *c);
        }
        let storage_rate = C64::new(0.0, omega) * energy.value();
        let forcing = C64::new(0.0, omega * area);
        self.relative_heat_defect = relative(
            (storage_rate + self.bath_flux - forcing).abs(),
            storage_rate.abs() + self.bath_flux.abs() + forcing.abs(),
        );
        let work = self.forms.iter().map(|f| f.conduction).sum::<f64>()
            - omega * capacity[0].im * self.forms[0].field
            + bath
                .iter()
                .zip(self.bath_forms)
                .map(|(g, f)| g * f)
                .sum::<f64>();
        let input = omega * self.forms[0].temperature_integral.im;
        self.relative_work_defect = relative((work - input).abs(), work.abs() + input.abs());
        if self.forms.iter().any(|f| {
            !finite(f.temperature_integral)
                || !f.field.is_finite()
                || f.field < 0.0
                || !f.conduction.is_finite()
                || f.conduction < 0.0
        }) || !finite(self.bath_flux)
            || !finite(self.mean_gas)
            || self
                .bath_forms
                .into_iter()
                .any(|f| !f.is_finite() || f < 0.0)
        {
            return Err(
                bad("rectangular original physical thermal forms are unrepresentable").into(),
            );
        }
        Ok(())
    }
}

pub(super) fn quadratic(k: [[f64; 2]; 2], gradient: [C64; 2]) -> Result<f64, Error> {
    let schur = k[1][1] - k[0][1] * (k[0][1] / k[0][0]);
    if !positive(schur) {
        return Err(bad("original tensor conduction metric is unrepresentable").into());
    }
    let a = gradient[0] + gradient[1].scale(k[0][1] / k[0][0]);
    Ok(a.scale(det::sqrt(k[0][0])).norm_sq() + gradient[1].scale(det::sqrt(schur)).norm_sq())
}
fn principal_bounds(k: [[f64; 2]; 2]) -> Result<(f64, f64), Error> {
    if k.iter().flatten().any(|v| !v.is_finite())
        || k[0][1] != k[1][0]
        || !positive(k[0][0])
        || !positive(k[1][1])
    {
        return Err(bad("original finite symmetric positive conductivity tensor required").into());
    }
    let scale = k[0][0].max(k[1][1]).max(k[0][1].abs());
    let [a, b, d] = [k[0][0] / scale, k[0][1] / scale, k[1][1] / scale];
    let high = 0.5 * (a + d + det::hypot(a - d, 2.0 * b));
    let low = (a * d - b * b) / high;
    if !positive(low) {
        return Err(bad("original solid conductivity must be positive definite").into());
    }
    Ok((low * scale, high * scale))
}
fn cell_nodes(nx: usize, x: usize, y: usize) -> [usize; 4] {
    let a = y * (nx + 1) + x;
    [a, a + 1, a + nx + 1, a + nx + 2]
}
fn shape(a: usize, f: [f64; 2], h: [f64; 2]) -> [f64; 3] {
    let x = if a & 1 == 0 { 1.0 - f[0] } else { f[0] };
    let y = if a >> 1 == 0 { 1.0 - f[1] } else { f[1] };
    [
        x * y,
        if a & 1 == 0 { -y / h[0] } else { y / h[0] },
        if a >> 1 == 0 { -x / h[1] } else { x / h[1] },
    ]
}
fn element(h: [f64; 2], k: [[f64; 2]; 2], capacity: C64, omega: f64) -> [[C64; 4]; 4] {
    let gauss = [0.5 - 0.5 / det::sqrt(3.0), 0.5 + 0.5 / det::sqrt(3.0)];
    let mut out = [[C64::ZERO; 4]; 4];
    for a in 0..4 {
        for b in 0..=a {
            let mut conduction = Dd::ZERO;
            for x in gauss {
                for y in gauss {
                    let u = shape(a, [x, y], h);
                    let v = shape(b, [x, y], h);
                    conduction = conduction
                        + Dd::from_f64(h[0] * h[1] * 0.25)
                            * Dd::from_f64(
                                k[0][0] * u[1] * v[1]
                                    + k[1][1] * u[2] * v[2]
                                    + k[0][1] * (u[1] * v[2] + u[2] * v[1]),
                            );
                }
            }
            let mass = h[0]
                * h[1]
                * if a & 1 == b & 1 { 1.0 / 3.0 } else { 1.0 / 6.0 }
                * if a >> 1 == b >> 1 {
                    1.0 / 3.0
                } else {
                    1.0 / 6.0
                };
            let value = C64::from_re(conduction.to_f64())
                + C64::new(-capacity.im, capacity.re).scale(omega * mass);
            out[a][b] = value;
            out[b][a] = value;
        }
    }
    out
}
fn add(coo: &mut Coo, map: &[usize], a: usize, b: usize, value: C64) {
    let (a, b) = (map[a], map[b]);
    if a == usize::MAX || b == usize::MAX {
        return;
    }
    coo.push(2 * a, 2 * b, value.re);
    coo.push(2 * a, 2 * b + 1, -value.im);
    coo.push(2 * a + 1, 2 * b, -value.im);
    coo.push(2 * a + 1, 2 * b + 1, -value.re);
}
fn add_rhs(rhs: &mut [f64], row: usize, value: C64) {
    rhs[2 * row] += value.re;
    rhs[2 * row + 1] -= value.im;
}
#[derive(Clone, Copy, Default)]
struct Sum {
    re: Dd,
    im: Dd,
}
impl Sum {
    fn add(&mut self, value: C64) {
        self.re = self.re + Dd::from_f64(value.re);
        self.im = self.im + Dd::from_f64(value.im);
    }
    fn value(self) -> C64 {
        C64::new(self.re.to_f64(), self.im.to_f64())
    }
}
fn positive(v: f64) -> bool {
    v.is_finite() && v > 0.0
}
fn relative(error: f64, scale: f64) -> f64 {
    if scale == 0.0 {
        if error == 0.0 { 0.0 } else { f64::INFINITY }
    } else {
        error / scale
    }
}
fn filled<T: Clone>(length: usize, value: T) -> Result<Vec<T>, Error> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(length)
        .map_err(|_| bad("rectangular jacket allocation refused"))?;
    result.resize(length, value);
    Ok(result)
}

#[cfg(test)]
mod tests;
