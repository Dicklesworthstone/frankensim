//! Spatially integrated irreversible film heat with explicit physical support.
//!
//! Transverse integration acts on actual coherent internal fields, before
//! quadratic tones are emitted. The complete periodic cell can be integrated
//! continuously in the original Fourier basis without synthesizing samples.
//! Increasing depth/grid resolution does not multiply the retained spectral
//! term count. Material reservoirs and clocks remain those of the prepared layer.

use super::heat_spectrum::{nonnegative_heat, reservoir_range, spectrum_error};
use super::*;
use fs_fft::harmonic_window::{
    AngularHarmonicSpectrumTerm, HarmonicProductTone, HarmonicSpectrumTermLimits,
    MAX_PRODUCT_SPECTRUM_TERMS, PreparedHarmonicProductSpectrum,
};

mod rectangle;

/// One authored positive depth-quadrature weight, in metres.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AngularSpectrumHeatDepthSample {
    pub depth_m: f64,
    pub weight_m: f64,
}

/// Physical transverse support on the ORIGINAL row-major periodic grid.
#[derive(Clone, Copy, Debug)]
pub enum AngularSpectrumHeatAreaQuadrature<'q> {
    /// Every original FIELD SAMPLE contributes its full dx*dy area. Quadratic
    /// products can alias even when the represented electric field is resolved.
    PeriodicCell,
    /// Continuous integral over the full original periodic cell. Difference
    /// products pair equal signed Fourier modes; sum products pair opposite
    /// modes. A negative Nyquist mode has no positive partner. This avoids
    /// transverse product aliasing and requires no spatial synthesis FFTs.
    /// Depth integration still uses the explicitly authored quadrature below.
    ContinuousPeriodicCell,
    /// Exact continuous integral over a uniform axis-aligned rectangle in the
    /// original x/y coordinates. Both sum and difference spatial interference
    /// are retained, including between distinct Fourier modes. The rectangle
    /// observes the original periodic continuation and may cross cell edges;
    /// it does not impose a new field/material boundary at the patch edges.
    ContinuousRectangle(AngularSpectrumAperture),
    /// Per-cell areas [m²], each in [0,dx*dy]. Zero masks a cell; positive
    /// fractions describe partial cells without relocating field samples.
    CellAreasM2(&'q [f64]),
}

/// Caller-authored physical volume quadrature, retained without rescaling.
///
/// The positive weights must sum to depth_span_m[1]-depth_span_m[0] within
/// 32*f64::EPSILON relative allowance. Every sample must lie in that declared
/// span, itself inside the original film. Actual weight sums determine the
/// represented volume and the load; no normalization silently changes them.
#[derive(Clone, Copy, Debug)]
pub struct AngularSpectrumHeatVolumeQuadrature<'q> {
    pub depth_span_m: [f64; 2],
    pub depth_samples: &'q [AngularSpectrumHeatDepthSample],
    pub transverse: AngularSpectrumHeatAreaQuadrature<'q>,
}

/// Total query work, including repeated original field reconstruction.
#[derive(Clone, Copy, Debug)]
pub struct AngularSpectrumIntegratedHeatLimits {
    pub maximum_depth_samples: usize,
    /// Depth samples times the original owner's point/source-line work.
    pub maximum_point_line_pairs: usize,
    /// Repeated field synthesis; zero for either continuous support. Original
    /// beam analysis and layer preparation retain their own admission bounds.
    pub maximum_fft_work: u64,
    /// Depths * original grid points * sum F_r*(F_r+1)/2. Each visits three
    /// vector components, with conjugate and (for real banks) plain products.
    /// ContinuousPeriodicCell visits Fourier modes instead of field samples.
    /// ContinuousRectangle visits N² ordered mode pairs plus one center-field
    /// product at each frequency pair and depth; it charges N²+1 rather than N.
    pub maximum_spatial_point_frequency_pairs: usize,
    /// Selected reservoirs * sum F_r*(F_r+1)/2, cached once for all depths.
    pub maximum_reservoir_frequency_pairs: usize,
    /// Two tones per upper frequency pair for real banks, one for proper.
    pub maximum_generated_terms: usize,
    /// Conservative original complete workspace plus retained frequency
    /// planes, pair caches, raw contracted terms and shared tone preparation.
    /// Borrowed quadrature/source inputs and allocator overhead are excluded.
    pub maximum_workspace_bytes: usize,
}

/// Expected irreversible thermal load [W] from an explicit film volume.
///
/// This retains original material/field ownership, physical support and actual
/// quadrature weights. ContinuousPeriodicCell integrates the represented
/// transverse Fourier field exactly. ContinuousRectangle integrates its actual
/// finite patch exactly and retains cross-mode interference; sampled area
/// options can alias spatial products. Depth quadrature must resolve
/// standing waves and attenuation independently. No physical field/depth
/// convergence certificate or isolated finite-aperture model is inferred.
/// Reversible field/polarization storage and boundary momentum work are absent.
/// Do not also debit this same loss via local heat or integrated absorption.
pub struct AngularSpectrumExposureLayerHeatLoadSpectrum<'p, 'q, 'e, 'b, 's, 'm> {
    layer: &'p PreparedAngularSpectrumExposureLayer<'e, 'b, 's, 'm>,
    quadrature: AngularSpectrumHeatVolumeQuadrature<'q>,
    channel: AngularSpectrumExposureHeatChannel,
    spectrum: PreparedHarmonicProductSpectrum,
    area_m2: Dd,
    depth_weight_m: Dd,
    volume_m3: Dd,
    workspace_bytes: usize,
}

impl<'p, 'q, 'e, 'b, 's, 'm> AngularSpectrumExposureLayerHeatLoadSpectrum<'p, 'q, 'e, 'b, 's, 'm> {
    pub fn layer(&self) -> &'p PreparedAngularSpectrumExposureLayer<'e, 'b, 's, 'm> {
        self.layer
    }
    pub fn grid(&self) -> AngularSpectrumGrid {
        self.layer.anchor.grid()
    }
    pub fn quadrature(&self) -> AngularSpectrumHeatVolumeQuadrature<'q> {
        self.quadrature
    }
    pub fn channel(&self) -> AngularSpectrumExposureHeatChannel {
        self.channel
    }
    /// Physical full-cell/rectangle area or sum of authored cell areas [m²].
    pub fn represented_area_m2_dd(&self) -> Dd {
        self.area_m2
    }
    /// Actual sum of the authored depth weights [m], without normalization.
    pub fn represented_depth_m_dd(&self) -> Dd {
        self.depth_weight_m
    }
    pub fn represented_volume_m3_dd(&self) -> Dd {
        self.volume_m3
    }
    /// Already in W: downstream thermal owners must not multiply by volume.
    pub fn spectrum(&self) -> &PreparedHarmonicProductSpectrum {
        &self.spectrum
    }
    /// Move the already prepared watt spectrum into a downstream owner. The
    /// consumer must retain its own material, quadrature and source-clock binding;
    /// detaching coefficients neither refreshes fields nor changes their units.
    pub fn into_spectrum(self) -> PreparedHarmonicProductSpectrum {
        self.spectrum
    }
    pub fn preparation_workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }
    pub fn evaluate(&self, time_s: f64) -> Result<f64, AngularSpectrumError> {
        self.evaluate_cancellable(time_s, || false)
    }
    pub fn evaluate_cancellable(
        &self,
        time_s: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<f64, AngularSpectrumError> {
        nonnegative_heat(
            self.spectrum
                .evaluate_cancellable(time_s, cancelled)
                .map_err(spectrum_error)?,
        )
    }
    pub fn observe_window(&self, window: HarmonicWindow) -> Result<f64, AngularSpectrumError> {
        self.observe_window_cancellable(window, || false)
    }
    pub fn observe_window_cancellable(
        &self,
        window: HarmonicWindow,
        cancelled: impl FnMut() -> bool,
    ) -> Result<f64, AngularSpectrumError> {
        nonnegative_heat(
            self.spectrum
                .mean_cancellable(window, cancelled)
                .map_err(spectrum_error)?,
        )
    }
}

#[derive(Clone, Copy)]
struct VolumeAdmission {
    grid: AngularSpectrumGrid,
    cell_area_m2: Dd,
    area_m2: Dd,
    depth_weight_m: Dd,
    volume_m3: Dd,
    points: usize,
}

#[derive(Clone, Copy, Default)]
struct PairIntegral {
    difference: [ScaledDd; 2],
    sum: [ScaledDd; 2],
}

struct FrequencyPair {
    left: usize,
    right: usize,
    reservoir: PairIntegral,
    spatial: PairIntegral,
}

impl<'e, 'b, 's, 'm> PreparedAngularSpectrumExposureLayer<'e, 'b, 's, 'm> {
    /// Integrate actual internal heat over a declared subslab/full film and
    /// full/partial original cells or a continuous rectangle. Each bank-frequency
    /// plane is reconstructed once per depth, independent of the reservoir count.
    /// ContinuousPeriodicCell retains its modal vectors and integrates signed
    /// Fourier products directly; a continuous rectangle retains all mode pairs
    /// and their original spatial phases. Sampled supports synthesize fields.
    ///
    /// Reservoir transfer products are cached per frequency pair; spatial
    /// E_f dot conj(E_g) and E_f dot E_g are contracted once and then weighted
    /// by the SUM of the selected reservoir coefficients. Upper triangular
    /// frequency pairs combine the two symmetric off-diagonal contributions.
    /// Tone storage is O(sum F_r²), independent of quadrature sample counts.
    pub fn integrated_heat_spectrum<'p, 'q>(
        &'p self,
        quadrature: AngularSpectrumHeatVolumeQuadrature<'q>,
        channel: AngularSpectrumExposureHeatChannel,
        limits: AngularSpectrumIntegratedHeatLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<
        AngularSpectrumExposureLayerHeatLoadSpectrum<'p, 'q, 'e, 'b, 's, 'm>,
        AngularSpectrumError,
    > {
        poll(&mut cancelled)?;
        if quadrature.depth_samples.len() > limits.maximum_depth_samples {
            return Err(bad(
                "integrated heat depth quadrature exceeds its sample budget",
            ));
        }
        let volume = admit_volume(
            self.anchor.grid(),
            self.thickness_m,
            quadrature,
            &mut cancelled,
        )?;
        let (first, reservoirs) =
            reservoir_range(channel, self.law.conductivity, self.law.poles.len())?;
        let active = reservoirs != 0 && volume.area_m2.to_f64() > 0.0;
        let depths = if active {
            quadrature.depth_samples.len()
        } else {
            0
        };
        let mut pair_count = 0usize;
        let mut generated = 0usize;
        let mut maximum_pairs = 0usize;
        let mut maximum_frequencies = 0usize;
        if active {
            for (bank, proper) in ensemble::banks(self.mean, self.factors, self.statistics) {
                poll(&mut cancelled)?;
                let f = bank.frequency_count();
                let pairs = upper_pairs(f)?;
                pair_count = pair_count
                    .checked_add(pairs)
                    .ok_or_else(|| bad("integrated heat pair work overflows"))?;
                maximum_pairs = maximum_pairs.max(pairs);
                maximum_frequencies = maximum_frequencies.max(f);
                generated = pairs
                    .checked_mul(if proper { 1 } else { 2 })
                    .and_then(|v| generated.checked_add(v))
                    .ok_or_else(|| bad("integrated heat spectrum count overflows"))?;
            }
        }
        let original_line_work = self
            .point_line_pairs
            .checked_mul(depths)
            .ok_or_else(|| bad("integrated heat original field work overflows"))?;
        let continuous = matches!(
            quadrature.transverse,
            AngularSpectrumHeatAreaQuadrature::ContinuousPeriodicCell
                | AngularSpectrumHeatAreaQuadrature::ContinuousRectangle(_)
        );
        let rectangle = match quadrature.transverse {
            AngularSpectrumHeatAreaQuadrature::ContinuousRectangle(aperture)
                if active && pair_count != 0 =>
            {
                Some(aperture)
            }
            _ => None,
        };
        let original_fft_work = if continuous {
            0
        } else {
            u64::try_from(depths)
                .ok()
                .and_then(|d| self.field_fft_work.checked_mul(d))
                .ok_or_else(|| bad("integrated heat FFT work overflows"))?
        };
        let spatial_points = if rectangle.is_some() {
            volume
                .points
                .checked_mul(volume.points)
                .and_then(|v| v.checked_add(1))
                .ok_or_else(|| bad("rectangular heat mode-pair work overflows"))?
        } else {
            volume.points
        };
        let spatial_work = pair_count
            .checked_mul(spatial_points)
            .and_then(|v| v.checked_mul(depths))
            .ok_or_else(|| bad("integrated heat spatial work overflows"))?;
        let reservoir_work = pair_count
            .checked_mul(reservoirs)
            .ok_or_else(|| bad("integrated heat reservoir work overflows"))?;
        let raw_bytes = generated
            .checked_mul(core::mem::size_of::<AngularHarmonicSpectrumTerm>())
            .ok_or_else(|| bad("integrated heat spectrum storage overflows"))?;
        let tone_bytes = generated
            .checked_mul(core::mem::size_of::<HarmonicProductTone>())
            .ok_or_else(|| bad("integrated heat tone storage overflows"))?;
        let rectangle_bytes = if rectangle.is_some() {
            rectangle::PreparedRectangle::workspace_bytes(volume.grid)?
                .checked_add(
                    volume
                        .points
                        .checked_mul(core::mem::size_of::<[[Dd; 2]; 6]>())
                        .ok_or_else(|| bad("rectangular heat DD field workspace overflows"))?,
                )
                .ok_or_else(|| bad("rectangular heat workspace overflows"))?
        } else {
            0
        };
        let field_bytes = if rectangle.is_some() {
            core::mem::size_of::<rectangle::Electric>()
        } else {
            core::mem::size_of::<[C64; 6]>()
        };
        let field_metadata_bytes = if rectangle.is_some() {
            core::mem::size_of::<rectangle::Plane>()
        } else {
            core::mem::size_of::<Vec<[C64; 6]>>()
        };
        let workspace_bytes = maximum_frequencies
            .checked_mul(volume.points)
            .and_then(|v| v.checked_mul(field_bytes))
            .and_then(|v| {
                v.checked_add(
                    maximum_frequencies
                        .checked_mul(core::mem::size_of::<usize>() + field_metadata_bytes)?,
                )
            })
            .and_then(|v| {
                v.checked_add(maximum_pairs.checked_mul(core::mem::size_of::<FrequencyPair>())?)
            })
            .and_then(|v| v.checked_add(raw_bytes))
            .and_then(|v| v.checked_add(tone_bytes))
            .and_then(|v| v.checked_add(self.workspace_bytes))
            .and_then(|v| v.checked_add(rectangle_bytes))
            .ok_or_else(|| bad("integrated heat workspace overflows"))?;
        if original_line_work > limits.maximum_point_line_pairs
            || original_fft_work > limits.maximum_fft_work
            || spatial_work > limits.maximum_spatial_point_frequency_pairs
            || reservoir_work > limits.maximum_reservoir_frequency_pairs
            || generated > limits.maximum_generated_terms
            || generated > MAX_PRODUCT_SPECTRUM_TERMS
            || workspace_bytes > limits.maximum_workspace_bytes
        {
            return Err(bad(
                "integrated heat exceeds its original/spatial/spectral budget",
            ));
        }
        let rectangle = rectangle
            .map(|aperture| {
                rectangle::PreparedRectangle::prepare(
                    volume.grid,
                    aperture,
                    |i| self.anchor.modes[i].kx,
                    |i| self.anchor.modes[i * volume.grid.width].ky,
                    &mut cancelled,
                )
            })
            .transpose()?;
        let mut terms = Vec::new();
        terms
            .try_reserve_exact(generated)
            .map_err(|_| bad("integrated heat term allocation failed"))?;
        if active {
            for (bank, proper) in ensemble::banks(self.mean, self.factors, self.statistics) {
                poll(&mut cancelled)?;
                if bank.frequency_count() == 0 {
                    continue;
                }
                let mut mapping = Vec::new();
                mapping
                    .try_reserve_exact(bank.frequency_count())
                    .map_err(|_| bad("integrated heat frequency map allocation failed"))?;
                for &omega in &bank.frequencies {
                    poll(&mut cancelled)?;
                    mapping.push(
                        self.frequencies
                            .binary_search_by(|cache| cache.omega.total_cmp(&omega))
                            .map_err(|_| bad("original film frequency cache is missing"))?,
                    );
                }
                let reservoir = |i: usize, j: usize, cancelled: &mut _| {
                    reservoir_pair(
                        self.law.conductivity,
                        &self.frequencies[mapping[i]].poles,
                        &self.frequencies[mapping[j]].poles,
                        first,
                        reservoirs,
                        proper,
                        cancelled,
                    )
                };
                if let Some(rectangle) = &rectangle {
                    rectangle.contract_bank(
                        &bank.frequencies,
                        proper,
                        quadrature,
                        volume,
                        reservoir,
                        |group, depth, cancelled| {
                            let layer = &self.frequencies[mapping[group]].layer;
                            bank.frequency_modal_parts_with(group, cancelled, |line, cancelled| {
                                line.modal_layer_field(
                                    layer.as_ref(),
                                    self.layer_index,
                                    depth,
                                    cancelled,
                                )
                            })
                        },
                        &mut terms,
                        &mut cancelled,
                    )?;
                    continue;
                }
                contract_bank(
                    &bank.frequencies,
                    proper,
                    quadrature,
                    volume,
                    reservoir,
                    |group, depth, cancelled| {
                        let layer = &self.frequencies[mapping[group]].layer;
                        let mut field = bank.frequency_modal_field_with(
                            group,
                            cancelled,
                            |line, cancelled| {
                                line.modal_layer_field(
                                    layer.as_ref(),
                                    self.layer_index,
                                    depth,
                                    cancelled,
                                )
                            },
                        )?;
                        if !continuous {
                            bank.synthesize_frequency_field(group, &mut field, cancelled)?;
                        }
                        Ok(field)
                    },
                    &mut terms,
                    &mut cancelled,
                )?;
            }
        }
        let spectrum = PreparedHarmonicProductSpectrum::from_angular_terms_cancellable(
            terms,
            HarmonicSpectrumTermLimits {
                maximum_terms: generated,
                maximum_workspace_bytes: tone_bytes,
            },
            &mut cancelled,
        )
        .map_err(spectrum_error)?;
        poll(&mut cancelled)?;
        Ok(AngularSpectrumExposureLayerHeatLoadSpectrum {
            layer: self,
            quadrature,
            channel,
            spectrum,
            area_m2: volume.area_m2,
            depth_weight_m: volume.depth_weight_m,
            volume_m3: volume.volume_m3,
            workspace_bytes,
        })
    }
}

fn upper_pairs(f: usize) -> Result<usize, AngularSpectrumError> {
    f.checked_add(1)
        .and_then(|g| f.checked_mul(g))
        .map(|v| v / 2)
        .ok_or_else(|| bad("integrated heat frequency pair count overflows"))
}

fn admit_volume(
    grid: AngularSpectrumGrid,
    thickness: f64,
    quadrature: AngularSpectrumHeatVolumeQuadrature<'_>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<VolumeAdmission, AngularSpectrumError> {
    let [start, end] = quadrature.depth_span_m;
    if !start.is_finite()
        || !end.is_finite()
        || start < 0.0
        || end > thickness
        || end <= start
        || quadrature.depth_samples.is_empty()
    {
        return Err(bad(
            "integrated heat needs a nonempty depth span within the original film",
        ));
    }
    let mut depth_weight_m = Dd::ZERO;
    for sample in quadrature.depth_samples {
        poll(cancelled)?;
        if !sample.depth_m.is_finite()
            || sample.depth_m < start
            || sample.depth_m > end
            || !sample.weight_m.is_finite()
            || sample.weight_m <= 0.0
        {
            return Err(bad(
                "integrated heat depth nodes/weights are outside the physical span",
            ));
        }
        depth_weight_m = depth_weight_m + Dd::from_f64(sample.weight_m);
    }
    let span = Dd::from_f64(end) - Dd::from_f64(start);
    if !finite_dd(depth_weight_m)
        || (depth_weight_m - span).to_f64().abs() > 32.0 * f64::EPSILON * span.to_f64()
    {
        return Err(bad(
            "integrated heat depth weights must integrate the declared span",
        ));
    }
    let points = grid
        .width
        .checked_mul(grid.height)
        .ok_or_else(|| bad("integrated heat grid overflows"))?;
    let cell_area_m2 = Dd::from_f64(grid.spacing_m[0]) * Dd::from_f64(grid.spacing_m[1]);
    if !finite_dd(cell_area_m2) || cell_area_m2.to_f64() <= 0.0 {
        return Err(bad("integrated heat original cell area is unrepresentable"));
    }
    let area_m2 = match quadrature.transverse {
        AngularSpectrumHeatAreaQuadrature::PeriodicCell
        | AngularSpectrumHeatAreaQuadrature::ContinuousPeriodicCell => {
            cell_area_m2 * Dd::from_f64(points as f64)
        }
        AngularSpectrumHeatAreaQuadrature::ContinuousRectangle(aperture) => {
            let extent = aperture.extent_m();
            Dd::from_f64(extent[0]) * Dd::from_f64(extent[1])
        }
        AngularSpectrumHeatAreaQuadrature::CellAreasM2(areas) => {
            if areas.len() != points {
                return Err(bad(
                    "integrated heat cell areas must match the original grid",
                ));
            }
            let maximum_area = grid.spacing_m[0] * grid.spacing_m[1];
            let mut total = Dd::ZERO;
            for &area in areas {
                poll(cancelled)?;
                if !area.is_finite() || area < 0.0 || area > maximum_area {
                    return Err(bad(
                        "integrated heat cell areas must be zero or positive partial-cell areas",
                    ));
                }
                total = total + Dd::from_f64(area);
            }
            total
        }
    };
    let volume_m3 = area_m2 * depth_weight_m;
    if !finite_dd(area_m2)
        || !finite_dd(volume_m3)
        || (area_m2.to_f64() > 0.0 && volume_m3.to_f64() <= 0.0)
    {
        return Err(bad("integrated heat represented volume is unrepresentable"));
    }
    poll(cancelled)?;
    Ok(VolumeAdmission {
        grid,
        cell_area_m2,
        area_m2,
        depth_weight_m,
        volume_m3,
        points,
    })
}

fn area_weight(
    transverse: AngularSpectrumHeatAreaQuadrature<'_>,
    volume: VolumeAdmission,
    point: usize,
) -> Dd {
    match transverse {
        AngularSpectrumHeatAreaQuadrature::PeriodicCell => volume.cell_area_m2,
        AngularSpectrumHeatAreaQuadrature::ContinuousPeriodicCell => volume.area_m2,
        AngularSpectrumHeatAreaQuadrature::ContinuousRectangle(_) => volume.area_m2,
        AngularSpectrumHeatAreaQuadrature::CellAreasM2(areas) => Dd::from_f64(areas[point]),
    }
}

fn reservoir_pair(
    sigma: f64,
    a: &[PoleTransfer],
    b: &[PoleTransfer],
    first: usize,
    count: usize,
    proper: bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PairIntegral, AngularSpectrumError> {
    let mut result = PairIntegral::default();
    for reservoir in first..first + count {
        poll(cancelled)?;
        if reservoir == 0 {
            result.difference[0] = result.difference[0] + ScaledDd::from_f64(sigma);
            if !proper {
                result.sum[0] = result.sum[0] + ScaledDd::from_f64(sigma);
            }
        } else {
            let left = parts(a[reservoir - 1].heat);
            let right = parts(b[reservoir - 1].heat);
            add_parts(&mut result.difference, multiply_conjugate(left, right));
            if !proper {
                add_parts(&mut result.sum, multiply(left, right));
            }
        }
    }
    if !finite_parts(result.difference) || !finite_parts(result.sum) {
        return Err(bad("integrated heat reservoir pair is unrepresentable"));
    }
    Ok(result)
}

fn prepare_pairs<F: FnMut() -> bool>(
    frequencies: &[f64],
    mut reservoir: impl FnMut(usize, usize, &mut F) -> Result<PairIntegral, AngularSpectrumError>,
    cancelled: &mut F,
) -> Result<Vec<FrequencyPair>, AngularSpectrumError> {
    let mut pairs = Vec::new();
    pairs
        .try_reserve_exact(upper_pairs(frequencies.len())?)
        .map_err(|_| bad("integrated heat pair cache allocation failed"))?;
    for i in 0..frequencies.len() {
        for j in i..frequencies.len() {
            poll(cancelled)?;
            pairs.push(FrequencyPair {
                left: i,
                right: j,
                reservoir: reservoir(i, j, cancelled)?,
                spatial: PairIntegral::default(),
            });
        }
    }
    Ok(pairs)
}

#[allow(clippy::too_many_arguments)]
fn contract_bank<F: FnMut() -> bool>(
    frequencies: &[f64],
    proper: bool,
    quadrature: AngularSpectrumHeatVolumeQuadrature<'_>,
    volume: VolumeAdmission,
    reservoir: impl FnMut(usize, usize, &mut F) -> Result<PairIntegral, AngularSpectrumError>,
    mut field: impl FnMut(usize, f64, &mut F) -> Result<Vec<[C64; 6]>, AngularSpectrumError>,
    terms: &mut Vec<AngularHarmonicSpectrumTerm>,
    cancelled: &mut F,
) -> Result<(), AngularSpectrumError> {
    poll(cancelled)?;
    if matches!(
        quadrature.transverse,
        AngularSpectrumHeatAreaQuadrature::ContinuousRectangle(_)
    ) {
        return Err(bad(
            "rectangular heat requires its original modal-pair integral",
        ));
    }
    if frequencies.is_empty() {
        return Ok(());
    }
    let continuous = matches!(
        quadrature.transverse,
        AngularSpectrumHeatAreaQuadrature::ContinuousPeriodicCell
    );
    let mut pairs = prepare_pairs(frequencies, reservoir, cancelled)?;
    let mut fields = Vec::new();
    fields
        .try_reserve_exact(frequencies.len())
        .map_err(|_| bad("integrated heat field plane allocation failed"))?;
    for sample in quadrature.depth_samples {
        poll(cancelled)?;
        fields.clear();
        for group in 0..frequencies.len() {
            poll(cancelled)?;
            let plane = field(group, sample.depth_m, cancelled)?;
            if plane.len() != volume.points {
                return Err(bad(
                    "integrated heat field plane differs from the original grid",
                ));
            }
            fields.push(plane);
        }
        for pair in &mut pairs {
            poll(cancelled)?;
            let mut spatial = PairIntegral::default();
            for point in 0..volume.points {
                if point % 64 == 0 {
                    poll(cancelled)?;
                }
                let area = area_weight(quadrature.transverse, volume, point);
                if area.to_f64() == 0.0 {
                    continue;
                }
                // This is continuous Fourier orthogonality, not modular DFT
                // index addition. In particular an unpaired negative Nyquist
                // bin cannot create a spatially constant sum-frequency tone.
                let partner = if proper {
                    None
                } else if continuous {
                    opposite_mode(volume.grid, point)
                } else {
                    Some(point)
                };
                let dot = spatial_pair(&fields[pair.left], &fields[pair.right], point, partner);
                add_scaled(&mut spatial.difference, dot.difference, area);
                if !proper {
                    add_scaled(&mut spatial.sum, dot.sum, area);
                }
            }
            let weight = Dd::from_f64(sample.weight_m);
            add_scaled(&mut pair.spatial.difference, spatial.difference, weight);
            if !proper {
                add_scaled(&mut pair.spatial.sum, spatial.sum, weight);
            }
        }
    }
    emit_pairs(frequencies, proper, pairs, terms, cancelled)
}

fn emit_pairs<F: FnMut() -> bool>(
    frequencies: &[f64],
    proper: bool,
    pairs: Vec<FrequencyPair>,
    terms: &mut Vec<AngularHarmonicSpectrumTerm>,
    cancelled: &mut F,
) -> Result<(), AngularSpectrumError> {
    for pair in pairs {
        poll(cancelled)?;
        // Each off-diagonal pair combines its symmetric counterpart:
        // 1/2+1/2=1; a diagonal contributes the original 1/2.
        let scale = ScaledDd::from_f64(if pair.left == pair.right { 0.5 } else { 1.0 });
        let a = Dd::from_f64(frequencies[pair.left]);
        let b = Dd::from_f64(frequencies[pair.right]);
        let difference = multiply(pair.spatial.difference, pair.reservoir.difference);
        let negative_time = physical_parts([difference[0] * scale, -difference[1] * scale])?;
        terms.push(AngularHarmonicSpectrumTerm {
            angular_frequency_rad_s: a - b,
            peak_amplitude: negative_time,
        });
        if !proper {
            let sum = multiply(pair.spatial.sum, pair.reservoir.sum);
            let negative_time = physical_parts([sum[0] * scale, -sum[1] * scale])?;
            terms.push(AngularHarmonicSpectrumTerm {
                angular_frequency_rad_s: a + b,
                peak_amplitude: negative_time,
            });
        }
    }
    poll(cancelled)
}

// Sampled products use the same point for both terms. Continuous cell products
// use equal modes for differences and a genuine opposite for sums. None also
// represents the absence of pseudocovariance in an authored proper factor.
fn spatial_pair(
    left: &[[C64; 6]],
    right: &[[C64; 6]],
    point: usize,
    partner: Option<usize>,
) -> PairIntegral {
    let mut dot = PairIntegral::default();
    for axis in 0..3 {
        let a = parts(left[point][axis]);
        add_parts(
            &mut dot.difference,
            multiply_conjugate(a, parts(right[point][axis])),
        );
        if let Some(partner) = partner {
            add_parts(&mut dot.sum, multiply(a, parts(right[partner][axis])));
        }
    }
    dot
}

fn parts(a: C64) -> [ScaledDd; 2] {
    [ScaledDd::from_f64(a.re), ScaledDd::from_f64(a.im)]
}
fn multiply(a: [ScaledDd; 2], b: [ScaledDd; 2]) -> [ScaledDd; 2] {
    [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]
}
fn multiply_conjugate(a: [ScaledDd; 2], b: [ScaledDd; 2]) -> [ScaledDd; 2] {
    [a[0] * b[0] + a[1] * b[1], a[1] * b[0] - a[0] * b[1]]
}
fn add_parts(sum: &mut [ScaledDd; 2], value: [ScaledDd; 2]) {
    for c in 0..2 {
        sum[c] = sum[c] + value[c];
    }
}
fn add_scaled(sum: &mut [ScaledDd; 2], value: [ScaledDd; 2], weight: Dd) {
    let weight = ScaledDd::from_dd(weight);
    for c in 0..2 {
        sum[c] = sum[c] + value[c] * weight;
    }
}
fn finite_dd(value: Dd) -> bool {
    value.hi.is_finite() && value.lo.is_finite() && value.to_f64().is_finite()
}
fn finite_parts(value: [ScaledDd; 2]) -> bool {
    value.into_iter().all(ScaledDd::is_finite)
}

fn physical_parts(value: [ScaledDd; 2]) -> Result<[Dd; 2], AngularSpectrumError> {
    Ok([value[0].physical()?, value[1].physical()?])
}

// Keep spatial and reservoir scales separate from the DD significand until
// their physical product is assembled. Squaring a tiny E before a large sigma,
// or squaring a tiny pole response before a large E, must not erase finite heat.
// Every real/imaginary accumulator has its own exponent. This extends dynamic
// range, not DD precision, and does not certify arbitrary cancellation nulls.
#[derive(Clone, Copy, Default)]
struct ScaledDd {
    mantissa: Dd,
    exponent: i32,
}

impl ScaledDd {
    const ZERO: Self = Self {
        mantissa: Dd::ZERO,
        exponent: 0,
    };
    const INVALID: Self = Self {
        mantissa: Dd {
            hi: f64::NAN,
            lo: 0.0,
        },
        exponent: 0,
    };

    fn from_f64(value: f64) -> Self {
        Self::from_dd(Dd::from_f64(value))
    }

    fn from_dd(value: Dd) -> Self {
        Self::normalized(value, 0)
    }

    fn normalized(value: Dd, exponent: i32) -> Self {
        if !value.hi.is_finite() || !value.lo.is_finite() {
            return Self::INVALID;
        }
        let value = Dd::from_pair(value.hi, value.lo);
        if !value.hi.is_finite() || !value.lo.is_finite() {
            return Self::INVALID;
        }
        if value.hi == 0.0 {
            return Self::ZERO;
        }
        let bits = value.hi.to_bits() & 0x7fff_ffff_ffff_ffff;
        let encoded = ((bits >> 52) & 0x7ff) as i32;
        let shift = if encoded == 0 {
            63 - bits.leading_zeros() as i32 - 1074
        } else {
            encoded - 1023
        };
        let Some(exponent) = exponent.checked_add(shift) else {
            return Self::INVALID;
        };
        Self {
            mantissa: Dd::from_pair(
                scale_by_two(value.hi, -shift),
                scale_by_two(value.lo, -shift),
            ),
            exponent,
        }
    }

    fn is_finite(self) -> bool {
        self.mantissa.hi.is_finite() && self.mantissa.lo.is_finite()
    }

    fn is_zero(self) -> bool {
        self.mantissa.hi == 0.0 && self.mantissa.lo == 0.0
    }

    // Restore units only after spatial, reservoir and diagonal scales meet.
    // A nonzero final tone outside the physical DD range must refuse, rather
    // than publishing a silently dark spectrum.
    fn physical(self) -> Result<Dd, AngularSpectrumError> {
        if !self.is_finite() {
            return Err(bad("integrated heat scaled contraction is unrepresentable"));
        }
        let value = Dd::from_pair(
            scale_by_two(self.mantissa.hi, self.exponent),
            scale_by_two(self.mantissa.lo, self.exponent),
        );
        if !finite_dd(value) || (!self.is_zero() && value.hi == 0.0 && value.lo == 0.0) {
            return Err(bad("integrated heat final watt tone is unrepresentable"));
        }
        Ok(value)
    }
}

impl core::ops::Neg for ScaledDd {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            mantissa: -self.mantissa,
            exponent: self.exponent,
        }
    }
}

impl core::ops::Add for ScaledDd {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        if !self.is_finite() || !other.is_finite() {
            return Self::INVALID;
        }
        if self.is_zero() {
            return other;
        }
        if other.is_zero() {
            return self;
        }
        let (large, small) = if self.exponent >= other.exponent {
            (self, other)
        } else {
            (other, self)
        };
        let Some(shift) = small.exponent.checked_sub(large.exponent) else {
            return Self::INVALID;
        };
        let aligned = Dd::from_pair(
            scale_by_two(small.mantissa.hi, shift),
            scale_by_two(small.mantissa.lo, shift),
        );
        Self::normalized(large.mantissa + aligned, large.exponent)
    }
}

impl core::ops::Sub for ScaledDd {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        self + (-other)
    }
}

impl core::ops::Mul for ScaledDd {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        if !self.is_finite() || !other.is_finite() {
            return Self::INVALID;
        }
        if self.is_zero() || other.is_zero() {
            return Self::ZERO;
        }
        let Some(exponent) = self.exponent.checked_add(other.exponent) else {
            return Self::INVALID;
        };
        Self::normalized(self.mantissa * other.mantissa, exponent)
    }
}

// Exact binary powers, at most three bounded steps. No platform pow/log is
// used. Outside +/-2098, every nonzero finite f64 necessarily over/underflows.
// Such underflow is harmless when aligning a term far below DD precision;
// restoring a nonzero physical coefficient checks it explicitly above.
fn scale_by_two(mut value: f64, mut exponent: i32) -> f64 {
    if value == 0.0 {
        return value;
    }
    if exponent > 2098 {
        return if value.is_sign_negative() {
            f64::NEG_INFINITY
        } else {
            f64::INFINITY
        };
    }
    if exponent < -2098 {
        return 0.0_f64.copysign(value);
    }
    let up = f64::from_bits(2046_u64 << 52); // 2^1023
    let down = f64::from_bits(1_u64 << 52); // 2^-1022
    while exponent > 1023 {
        value *= up;
        exponent -= 1023;
    }
    while exponent < -1022 {
        value *= down;
        exponent += 1022;
    }
    value * f64::from_bits(((exponent + 1023) as u64) << 52)
}

#[cfg(test)]
mod tests;
