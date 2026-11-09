//! Nonuniform thermal source cells inside the ORIGINAL optical layers.

use super::{
    AbsorbingLayer, AbsorbingSlab, Arc, ConductionError, GL8, GL16, MAX_SLAB_LAYERS,
    SlabDepthReport, SlabQuadrature, add, beam::optical_path, checkpoint, decayed_fraction,
    difference_product, exterior_difference_product, finite, invalid, layered, nonnegative,
    remaining_fraction,
};

/// Dense cell exchange is explicitly bounded; original optical layers are not split.
pub const MAX_SLAB_CELLS: usize = MAX_SLAB_LAYERS;

/// One thermal integration interval in an original material's local depth axis.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SlabCell {
    pub layer: usize,
    pub interval_m: [f64; 2],
}

/// Prepared radiative exchange for independently heated cells inside materials.
///
/// Cells partition the original layers completely. Each cell may have its own
/// uniform temperature; absorption/index stay those of its original material.
/// All reflection returns belong to ORIGINAL physical interfaces. No additional
/// interface, material, source, angular mixing or optical interpolation is made.
/// Within-layer exchange includes direct rays and every original reflected path.
///
/// Preparation costs O(Q*(L+C²)), including the original Snell critical bands;
/// storage is O(L²+C²), with L<=64, C<=64 and Q explicitly budgeted. Subsequent
/// applications cost O(C²), allocate nothing and perform no angular solves.
#[derive(Clone, Debug, PartialEq)]
pub struct SlabCellResponse {
    slab: AbsorbingSlab,
    cells: Arc<[SlabCell]>,
    face_coupling: Arc<[[f64; 2]]>,
    cell_coupling: Arc<[f64]>,
    report: SlabDepthReport,
}

impl SlabCellResponse {
    pub fn slab(&self) -> &AbsorbingSlab {
        &self.slab
    }
    pub fn cells(&self) -> &[SlabCell] {
        &self.cells
    }
    pub fn cell_count(&self) -> usize {
        self.cells.len()
    }
    /// Fractions of physical front/back irradiance absorbed in each cell.
    pub fn face_coupling(&self) -> &[[f64; 2]] {
        &self.face_coupling
    }
    /// Symmetric row-major thermal exchange, with a zero diagonal.
    pub fn cell_coupling(&self) -> &[f64] {
        &self.cell_coupling
    }
    pub const fn quadrature(&self) -> SlabDepthReport {
        self.report
    }

    /// Net heat INTO original integration cells [W/m² of common slab area].
    /// Exterior irradiation is physical flux; source exitances are vacuum
    /// blackbody powers in the same gray, finite-band or wavelength scope.
    /// Exterior n² and trapped modes are retained. Each internal pair supplies
    /// one equal/opposite heat; equal-temperature equilibrium cancels exactly.
    /// The destination changes only after the complete finite result is admitted.
    pub fn net_fluxes(
        &self,
        irradiation: [f64; 2],
        exitances: &[f64],
        out: &mut [f64],
    ) -> Result<(), ConductionError> {
        let n = self.cells.len();
        if exitances.len() != n
            || out.len() != n
            || irradiation.iter().chain(exitances).any(|x| !x.is_finite())
        {
            return Err(invalid(
                "cell heating needs finite fields on the prepared original-cell axis",
            ));
        }
        let mut heat = [0.0; MAX_SLAB_CELLS];
        let mut correction = [0.0; MAX_SLAB_CELLS];
        let indices = self.slab.exterior_refractive_indices();
        for i in 0..n {
            for side in 0..2 {
                add(
                    &mut heat[i],
                    &mut correction[i],
                    exterior_difference_product(
                        self.face_coupling[i][side],
                        irradiation[side],
                        exitances[i],
                        indices[side].powi(2),
                    )?,
                );
            }
            for j in i + 1..n {
                let into_i =
                    difference_product(self.cell_coupling[i * n + j], exitances[j], exitances[i])?;
                add(&mut heat[i], &mut correction[i], into_i);
                add(&mut heat[j], &mut correction[j], -into_i);
            }
        }
        for &value in &heat[..n] {
            finite(value)?;
        }
        out.copy_from_slice(&heat[..n]);
        Ok(())
    }
}

impl AbsorbingSlab {
    /// Prepare a complete ordered partition of the ORIGINAL physical layers.
    /// Cell intervals must be positive, adjacent and cover every layer exactly.
    /// A thermal temperature discontinuity creates no Fresnel interface.
    pub fn cell_response(
        &self,
        cells: Vec<SlabCell>,
        controls: SlabQuadrature,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SlabCellResponse, ConductionError> {
        checkpoint(&mut cancelled, 0)?;
        if !controls.is_valid() || !(1..=MAX_SLAB_CELLS).contains(&cells.len()) {
            return Err(invalid(
                "need 1..=64 original depth cells and bounded angular controls",
            ));
        }
        let (mut layer, mut depth) = (0, 0.0);
        for cell in &cells {
            if layer >= self.layer_count()
                || cell.layer != layer
                || cell.interval_m[0] != depth
                || cell.interval_m.iter().any(|x| !nonnegative(*x))
                || cell.interval_m[1] <= depth
                || cell.interval_m[1] > self.layers[layer].thickness_m()
            {
                return Err(invalid(
                    "thermal cells must partition every original layer without gaps or overlap",
                ));
            }
            depth = cell.interval_m[1];
            if depth == self.layers[layer].thickness_m() {
                layer += 1;
                depth = 0.0;
            }
        }
        if layer != self.layer_count() {
            return Err(invalid(
                "thermal cell partition does not cover all original material",
            ));
        }
        let n = cells.len();
        let channels = 2 * n + n * n;
        let mut total = vec![0.0; channels];
        let mut correction = vec![0.0; channels];
        let mut errors = vec![0.0; channels];
        let mut coarse = vec![0.0; channels];
        let mut fine = vec![0.0; channels];
        let mut rule_correction = vec![0.0; channels];
        let mut report = SlabDepthReport {
            controls,
            angular_evaluations: 0,
            accepted_panels: 0,
            maximum_estimated_coefficient_error: 0.0,
        };
        let exterior = self.exterior_refractive_indices();
        let mut edges = vec![0.0, exterior[0], exterior[1]];
        edges.extend(self.layers.iter().map(AbsorbingLayer::refractive_index));
        edges.sort_by(f64::total_cmp);
        edges.dedup();
        let bands: Vec<_> = edges.windows(2).map(|x| [x[0], x[1]]).collect();
        let mut pending = Vec::new();
        if self.layers.iter().any(|x| x.absorption_m_inv() > 0.0) {
            for band in (0..bands.len()).rev() {
                pending.push((0.0, 1.0, 0_usize, band));
            }
        }
        while let Some((lo, hi, depth, band)) = pending.pop() {
            checkpoint(&mut cancelled, report.angular_evaluations)?;
            if report.angular_evaluations + 24 > controls.maximum_evaluations {
                return Err(invalid("cell angular evaluation budget exhausted"));
            }
            for (nodes, output) in [(&GL8[..], &mut coarse), (&GL16[..], &mut fine)] {
                layered::rule(
                    &self.layers,
                    exterior,
                    bands[band],
                    lo,
                    hi,
                    nodes,
                    output,
                    &mut rule_correction,
                    &mut cancelled,
                    &mut report.angular_evaluations,
                    layered::Resolution::Cells(&cells),
                )?;
            }
            if coarse.iter().zip(&fine).all(|(&a, &b)| {
                (b - a).abs()
                    <= controls.absolute_coefficient_error * (hi - lo) / bands.len() as f64
                        + controls.relative_coefficient_error * b.abs()
            }) {
                for i in 0..channels {
                    add(&mut total[i], &mut correction[i], fine[i]);
                    errors[i] += (fine[i] - coarse[i]).abs();
                }
                report.accepted_panels += 1;
            } else {
                if depth >= controls.maximum_depth {
                    return Err(invalid("cell angular subdivision depth exhausted"));
                }
                let mid = 0.5 * (lo + hi);
                if mid <= lo || mid >= hi {
                    return Err(invalid("cell angular panel is unrepresentable"));
                }
                pending.push((mid, hi, depth + 1, band));
                pending.push((lo, mid, depth + 1, band));
            }
        }
        if total.iter().any(|&x| !nonnegative(x)) {
            return Err(invalid("unrepresentable cell radiative coupling"));
        }
        let face: Vec<_> = (0..n).map(|i| [total[2 * i], total[2 * i + 1]]).collect();
        let mut exchange = vec![0.0; n * n];
        for i in 0..n {
            for j in i + 1..n {
                exchange[i * n + j] = total[2 * n + i * n + j];
                exchange[j * n + i] = exchange[i * n + j];
            }
        }
        report.maximum_estimated_coefficient_error = errors.into_iter().fold(0.0_f64, f64::max);
        checkpoint(&mut cancelled, report.angular_evaluations)?;
        Ok(SlabCellResponse {
            slab: self.clone(),
            cells: cells.into(),
            face_coupling: face.into(),
            cell_coupling: exchange.into(),
            report,
        })
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn accumulate(
    layers: &[AbsorbingLayer],
    cells: &[SlabCell],
    mu: &[f64],
    t: &[f64],
    a: &[f64],
    r: &[f64],
    e: &[f64],
    exterior: [bool; 2],
    q: f64,
    exterior_q: [f64; 2],
    output: &mut [f64],
    correction: &mut [f64],
) -> Result<(), ConductionError> {
    let n = layers.len();
    let m = cells.len();
    let returns = layered::returns(t, a, r, e);
    let mut loss = [0.0; MAX_SLAB_CELLS];
    let mut passage = [0.0; MAX_SLAB_CELLS];
    let mut attenuation = [[0.0; 2]; MAX_SLAB_CELLS];
    let mut denominator = [0.0; MAX_SLAB_LAYERS];
    for i in 0..n {
        denominator[i] = a[i] * (1.0 + t[i])
            + t[i] * t[i] * (returns.left_deficit[i] + returns.left[i] * returns.right_deficit[i]);
    }
    for (i, cell) in cells.iter().enumerate() {
        let l = cell.layer;
        if mu[l] == 0.0 {
            continue;
        }
        let path = |distance| optical_path(layers[l].absorption_m_inv(), distance, mu[l]);
        let width = path(cell.interval_m[1] - cell.interval_m[0]);
        loss[i] = decayed_fraction(width);
        passage[i] = remaining_fraction(width);
        attenuation[i] = [
            remaining_fraction(path(cell.interval_m[0])),
            remaining_fraction(path(layers[l].thickness_m() - cell.interval_m[1])),
        ];
        if loss[i] > 0.0 && denominator[l] <= 0.0 {
            return Err(invalid("unrepresentable original-layer return denominator"));
        }
    }
    let mut face_fields = [[0.0; 2]; MAX_SLAB_LAYERS];
    if exterior[0] {
        let (_, _, mut field) = layered::joined(
            r[0],
            e[0],
            t[0],
            a[0],
            returns.right[0],
            returns.right_deficit[0],
        );
        for i in 0..n {
            face_fields[i][0] = field;
            if i + 1 < n {
                field *= t[i] * returns.cross_right[i];
            }
        }
    }
    if exterior[1] {
        let (_, _, mut field) = layered::joined(
            r[n],
            e[n],
            t[n - 1],
            a[n - 1],
            returns.left[n - 1],
            returns.left_deficit[n - 1],
        );
        for i in (0..n).rev() {
            face_fields[i][1] = field;
            if i > 0 {
                field *= t[i] * returns.cross_left[i];
            }
        }
    }
    for i in 0..m {
        if loss[i] == 0.0 {
            continue;
        }
        let l = cells[i].layer;
        let front = loss[i] * (attenuation[i][0] + t[l] * returns.right[l] * attenuation[i][1]);
        let back = loss[i] * (attenuation[i][1] + t[l] * returns.left[l] * attenuation[i][0]);
        for (side, capture) in [front, back].into_iter().enumerate() {
            add(
                &mut output[2 * i + side],
                &mut correction[2 * i + side],
                exterior_q[side] * face_fields[l][side] * capture,
            );
        }
        let emitted = (loss[i] / denominator[l])
            * (attenuation[i][1] + t[l] * returns.left[l] * attenuation[i][0]);
        let (mut target_layer, mut path) = (l, 1.0);
        for j in i + 1..m {
            let k = cells[j].layer;
            while target_layer < k {
                if target_layer > l {
                    path *= t[target_layer];
                }
                path *= returns.cross_right[target_layer];
                target_layer += 1;
            }
            if loss[j] == 0.0 {
                continue;
            }
            let exchange = if k == l {
                let gap = remaining_fraction(optical_path(
                    layers[l].absorption_m_inv(),
                    cells[j].interval_m[0] - cells[i].interval_m[1],
                    mu[l],
                ));
                // Factor the sum of direct, left, right and double-return
                // paths. Divide the source loss FIRST: weak trapped cavities
                // must not underflow a_i*a_j before their small denominator.
                (loss[i] / denominator[l])
                    * loss[j]
                    * gap
                    * (1.0 + returns.left[l] * passage[i] * attenuation[i][0] * attenuation[i][0])
                    * (1.0 + returns.right[l] * passage[j] * attenuation[j][1] * attenuation[j][1])
            } else {
                emitted
                    * path
                    * loss[j]
                    * (attenuation[j][0] + t[k] * returns.right[k] * attenuation[j][1])
            };
            let channel = 2 * m + i * m + j;
            add(&mut output[channel], &mut correction[channel], q * exchange);
        }
    }
    Ok(())
}
