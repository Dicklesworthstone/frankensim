//! Coherent film energy on explicitly authored physical volume quadrature.
//!
//! All samples at one depth share each original frequency reconstruction.
//! Fields interfere before positive volume and original material scales meet
//! their products. Repeated queries retain original coordinates, not bright
//! pair spectra or rounded point-energy densities.

use super::*;

const MAX_SAMPLES: usize = 65_536;

/// One physical volume weight at an ORIGINAL row-major transverse field sample
/// and an explicit depth in the selected film. Repeated points are permitted;
/// weights are never normalized and zero-weight entries still undergo admission.
/// This sampled rule does not certify spatial/depth quadrature accuracy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AngularSpectrumEnergyVolumeSample {
    pub depth_m: f64,
    pub point_index: usize,
    pub volume_m3: f64,
}

/// Bounded preparation of a reusable physical volume-energy observer.
#[derive(Clone, Copy, Debug)]
pub struct AngularSpectrumEnergyVolumeLimits {
    /// Authored samples including zero weights; hard ceiling 65,536.
    pub maximum_quadrature_points: usize,
    /// Distinct depths having at least one positive volume weight.
    pub maximum_depths: usize,
    /// Active depths times the original owner's point/source-line work.
    pub maximum_point_line_pairs: usize,
    /// Active depths times the original owner's field synthesis work.
    pub maximum_fft_work: u64,
    /// 3*Q_positive*sum_b F_b*(2+D+2L), including silent coordinates.
    pub maximum_coordinate_frequency_terms: usize,
    /// Original layer workspace plus retained coordinates/banks, active-index
    /// sorting and the largest temporary depth/frequency-point cache. Borrowed
    /// quadrature and allocator/call-stack overhead are excluded.
    pub maximum_workspace_bytes: usize,
}

/// Original electromagnetic and polarization energy in the authored volume [J].
/// For a change query these accounts are signed end-minus-start changes.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct AngularSpectrumExposureLayerVolumeEnergy {
    pub electromagnetic_j: Dd,
    pub polarization_restoring_j: Dd,
    pub polarization_kinetic_j: Dd,
    pub total_j: Dd,
}

/// Retained original-field energy on a fixed sampled volume in one film.
///
/// Queries allocate no heap and perform no material/field solve. The original
/// layer and complete authored quadrature remain borrowed. Sampled products may
/// alias even when the field grid is resolved; this owner does not replace the
/// continuous transverse Fourier measure or infer an enclosing surface.
pub struct PreparedAngularSpectrumExposureLayerVolumeEnergy<'p, 'q, 'e, 'b, 's, 'm> {
    layer: &'p PreparedAngularSpectrumExposureLayer<'e, 'b, 's, 'm>,
    quadrature: &'q [AngularSpectrumEnergyVolumeSample],
    banks: Vec<EnergyBank>,
    lines: Vec<EnergyLine>,
    volume_m3: Dd,
    depths: usize,
    coordinate_frequency_terms: usize,
    workspace_bytes: usize,
}

impl<'p, 'q, 'e, 'b, 's, 'm>
    PreparedAngularSpectrumExposureLayerVolumeEnergy<'p, 'q, 'e, 'b, 's, 'm>
{
    pub fn layer(&self) -> &'p PreparedAngularSpectrumExposureLayer<'e, 'b, 's, 'm> {
        self.layer
    }
    pub fn quadrature(&self) -> &'q [AngularSpectrumEnergyVolumeSample] {
        self.quadrature
    }
    pub fn represented_volume_m3_dd(&self) -> Dd {
        self.volume_m3
    }
    pub fn prepared_depth_count(&self) -> usize {
        self.depths
    }
    pub fn coordinate_frequency_terms(&self) -> usize {
        self.coordinate_frequency_terms
    }
    pub fn preparation_workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }

    /// Expected instantaneous storage, already in joules. Do not multiply the
    /// result by volume again. Original real/proper source laws are retained.
    pub fn evaluate_from_clock_dd(
        &self,
        time_s: Dd,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureLayerVolumeEnergy, AngularSpectrumError> {
        volume_accounts(&self.banks, &self.lines, time_s, None, &mut cancelled)
    }

    /// Signed physical storage change over the original positive interval.
    /// Analytic field increments precede weighted products and spatial sums;
    /// no subtraction of two bright endpoint energies is required.
    pub fn change_from_clock_dd(
        &self,
        start_s: Dd,
        duration_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureLayerVolumeEnergy, AngularSpectrumError> {
        poll(&mut cancelled)?;
        if !duration_s.is_finite() || duration_s <= 0.0 {
            return Err(bad(
                "film volume energy interval needs positive finite duration",
            ));
        }
        volume_accounts(
            &self.banks,
            &self.lines,
            start_s,
            Some(duration_s),
            &mut cancelled,
        )
    }
}

struct Admission {
    order: Vec<usize>,
    volume_m3: Dd,
    depths: usize,
    lines: usize,
    banks: usize,
    temporary_points: usize,
    coordinate_frequency_terms: usize,
    workspace_bytes: usize,
}

impl<'e, 'b, 's, 'm> PreparedAngularSpectrumExposureLayer<'e, 'b, 's, 'm> {
    /// Prepare original E/H and every authorized polarization coordinate on
    /// physical sampled volume quadrature. All samples at one depth share
    /// each source-bank frequency reconstruction, independent of point and
    /// reservoir counts. Positive volume scales are retained before powers.
    pub fn volume_energy<'p, 'q>(
        &'p self,
        quadrature: &'q [AngularSpectrumEnergyVolumeSample],
        limits: AngularSpectrumEnergyVolumeLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<
        PreparedAngularSpectrumExposureLayerVolumeEnergy<'p, 'q, 'e, 'b, 's, 'm>,
        AngularSpectrumError,
    > {
        let admission = admit(self, quadrature, limits, &mut cancelled)?;
        let mut points = allocated(admission.temporary_points)?;
        let mut lines = allocated(admission.lines)?;
        let mut banks = allocated(admission.banks)?;
        let mut first = 0;
        while first < admission.order.len() {
            poll(&mut cancelled)?;
            let depth_m = quadrature[admission.order[first]].depth_m;
            let mut end = first + 1;
            while end < admission.order.len() && quadrature[admission.order[end]].depth_m == depth_m
            {
                end += 1;
            }
            let samples = &admission.order[first..end];
            for (bank, proper) in ensemble::banks(self.mean, self.factors, self.statistics) {
                poll(&mut cancelled)?;
                if bank.frequency_count() == 0 {
                    continue;
                }
                points.clear();
                for group in 0..bank.frequency_count() {
                    poll(&mut cancelled)?;
                    let omega = bank.frequencies[group];
                    let cache = self
                        .frequencies
                        .binary_search_by(|v| v.omega.total_cmp(&omega))
                        .map_err(|_| bad("original film frequency cache is missing"))?;
                    let layer = &self.frequencies[cache].layer;
                    let fields =
                        bank.frequency_field_with(group, &mut cancelled, |line, cancelled| {
                            line.modal_layer_field(
                                layer.as_ref(),
                                self.layer_index,
                                depth_m,
                                cancelled,
                            )
                        })?;
                    for &index in samples {
                        poll(&mut cancelled)?;
                        points.push(PointField {
                            omega,
                            cache,
                            field: fields[quadrature[index].point_index],
                        });
                    }
                }
                for (offset, &index) in samples.iter().enumerate() {
                    retain_coordinates(
                        self,
                        points.iter().skip(offset).step_by(samples.len()),
                        proper,
                        ScaledDd::from_f64(quadrature[index].volume_m3),
                        &mut banks,
                        &mut lines,
                        &mut cancelled,
                    )?;
                }
            }
            first = end;
        }
        poll(&mut cancelled)?;
        Ok(PreparedAngularSpectrumExposureLayerVolumeEnergy {
            layer: self,
            quadrature,
            banks,
            lines,
            volume_m3: admission.volume_m3,
            depths: admission.depths,
            coordinate_frequency_terms: admission.coordinate_frequency_terms,
            workspace_bytes: admission.workspace_bytes,
        })
    }
}

fn admit(
    layer: &PreparedAngularSpectrumExposureLayer<'_, '_, '_, '_>,
    quadrature: &[AngularSpectrumEnergyVolumeSample],
    limits: AngularSpectrumEnergyVolumeLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Admission, AngularSpectrumError> {
    poll(cancelled)?;
    if quadrature.len() > MAX_SAMPLES || quadrature.len() > limits.maximum_quadrature_points {
        return Err(bad("film energy volume exceeds its quadrature budget"));
    }
    let mut volume_m3 = Dd::ZERO;
    let mut active = 0usize;
    for sample in quadrature {
        poll(cancelled)?;
        if !sample.depth_m.is_finite()
            || sample.depth_m < 0.0
            || sample.depth_m > layer.thickness_m
            || sample.point_index >= layer.anchor.modes.len()
            || !sample.volume_m3.is_finite()
            || sample.volume_m3 < 0.0
        {
            return Err(bad("film energy volume quadrature is invalid"));
        }
        if sample.volume_m3 != 0.0 {
            active += 1;
            volume_m3 = volume_m3 + Dd::from_f64(sample.volume_m3);
        }
    }
    if !volume_m3.hi.is_finite() || !volume_m3.lo.is_finite() {
        return Err(bad("film energy represented volume is unrepresentable"));
    }
    let index_bytes = active
        .checked_mul(core::mem::size_of::<usize>())
        .ok_or_else(overflow)?;
    // Bound the only preflight allocation before sorting active depths.
    // All field/cache allocation is admitted below, before reconstruction.
    if index_bytes
        .checked_add(layer.workspace_bytes)
        .ok_or_else(overflow)?
        > limits.maximum_workspace_bytes
    {
        return Err(bad("film energy volume exceeds its workspace budget"));
    }
    let mut order = allocated(active)?;
    for (i, sample) in quadrature.iter().enumerate() {
        poll(cancelled)?;
        if sample.volume_m3 != 0.0 {
            order.push(i);
        }
    }
    order.sort_unstable_by(|&a, &b| {
        quadrature[a]
            .depth_m
            .total_cmp(&quadrature[b].depth_m)
            .then_with(|| quadrature[a].point_index.cmp(&quadrature[b].point_index))
            .then_with(|| a.cmp(&b))
    });
    poll(cancelled)?;
    let mut depths = 0usize;
    let mut maximum_points = 0usize;
    let mut first = 0;
    while first < order.len() {
        poll(cancelled)?;
        let mut end = first + 1;
        while end < order.len()
            && quadrature[order[end]].depth_m == quadrature[order[first]].depth_m
        {
            end += 1;
        }
        depths += 1;
        maximum_points = maximum_points.max(end - first);
        first = end;
    }
    let original_lines = layer
        .point_line_pairs
        .checked_mul(depths)
        .ok_or_else(overflow)?;
    let fft_work = layer
        .field_fft_work
        .checked_mul(depths as u64)
        .ok_or_else(overflow)?;
    if depths > limits.maximum_depths
        || original_lines > limits.maximum_point_line_pairs
        || fft_work > limits.maximum_fft_work
    {
        return Err(bad("film energy volume exceeds its reconstruction budget"));
    }
    let coordinates = 2usize
        .checked_add(layer.law.pole_count())
        .and_then(|n| n.checked_add(layer.law.lorentz_poles.len()))
        .ok_or_else(overflow)?;
    let mut frequency_sum = 0usize;
    let mut source_banks = 0usize;
    let mut maximum_frequencies = 0usize;
    for (bank, _) in ensemble::banks(layer.mean, layer.factors, layer.statistics) {
        poll(cancelled)?;
        let f = bank.frequency_count();
        frequency_sum = frequency_sum.checked_add(f).ok_or_else(overflow)?;
        if f != 0 {
            source_banks = source_banks.checked_add(1).ok_or_else(overflow)?;
        }
        maximum_frequencies = maximum_frequencies.max(f);
    }
    let lines = active
        .checked_mul(coordinates)
        .and_then(|n| n.checked_mul(frequency_sum))
        .ok_or_else(overflow)?;
    let banks = active
        .checked_mul(coordinates)
        .and_then(|n| n.checked_mul(source_banks))
        .ok_or_else(overflow)?;
    let temporary_points = maximum_points
        .checked_mul(maximum_frequencies)
        .ok_or_else(overflow)?;
    let coordinate_frequency_terms = lines.checked_mul(3).ok_or_else(overflow)?;
    let workspace_bytes = lines
        .checked_mul(core::mem::size_of::<EnergyLine>())
        .and_then(|n| n.checked_add(banks.checked_mul(core::mem::size_of::<EnergyBank>())?))
        .and_then(|n| {
            n.checked_add(temporary_points.checked_mul(core::mem::size_of::<PointField>())?)
        })
        .and_then(|n| n.checked_add(index_bytes))
        .and_then(|n| n.checked_add(layer.workspace_bytes))
        .ok_or_else(overflow)?;
    if coordinate_frequency_terms > limits.maximum_coordinate_frequency_terms
        || workspace_bytes > limits.maximum_workspace_bytes
    {
        return Err(bad(
            "film energy volume exceeds its coordinate or workspace budget",
        ));
    }
    poll(cancelled)?;
    Ok(Admission {
        order,
        volume_m3,
        depths,
        lines,
        banks,
        temporary_points,
        coordinate_frequency_terms,
        workspace_bytes,
    })
}

fn volume_accounts(
    banks: &[EnergyBank],
    lines: &[EnergyLine],
    time: Dd,
    duration: Option<f64>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<AngularSpectrumExposureLayerVolumeEnergy, AngularSpectrumError> {
    let [
        electromagnetic_j,
        polarization_restoring_j,
        polarization_kinetic_j,
        total_j,
    ] = energy_accounts(banks, lines, time, duration, cancelled)?;
    Ok(AngularSpectrumExposureLayerVolumeEnergy {
        electromagnetic_j,
        polarization_restoring_j,
        polarization_kinetic_j,
        total_j,
    })
}

#[cfg(test)]
mod tests;
