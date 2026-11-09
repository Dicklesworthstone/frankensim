//! Explicit field-preserving temporal integration of the original volume.
//!
//! Each depth reconstructs its original frequency planes once. A bounded
//! positive temporal factor acts before spatial products. Continuous supports
//! retain complex Fourier quadratures; they are never replaced by samples.

use super::*;
use fs_fft::harmonic_window::{
    MAX_UNIFORM_HARMONIC_BAND_MODES, MAX_UNIFORM_HARMONIC_BAND_TONES, UniformHarmonicBandRef,
    UniformHarmonicBandWorkspace, phase_rad_s_dd_from_clock_dd,
};
use rectangle::Electric;

/// Explicit reconstruction costs for a field-preserving volume observation.
/// All limits are checked before reconstruction. The old cached spectrum query
/// keeps its original cost; selecting this API opts into fresh field work.
#[derive(Clone, Copy, Debug)]
pub struct AngularSpectrumHeatWindowLimits {
    pub maximum_depth_samples: usize,
    pub maximum_point_line_pairs: usize,
    pub maximum_fft_work: u64,
    /// Three axes * depths * points * reservoirs * frequencies * temporal
    /// modes. Real continuous fields require two projections to retain both
    /// spatial quadratures; other admitted fields require one.
    pub maximum_temporal_projection_terms: usize,
    /// Three axes per temporal-mode spatial integral, plus reservoir pair
    /// preparation and the original spatial pair work for unsupported banks.
    /// Rectangles charge ordered N*N pairs, center sums and center products,
    /// plus the width²+height² axis weights prepared once per observation.
    pub maximum_spatial_products: usize,
    /// Original quadratic terms from banks outside the positive factor's
    /// temporal admission. Admission is independent for each source bank.
    pub maximum_fallback_terms: usize,
    /// Original owner workspace, rectangle caches, at most one depth's F
    /// original and M temporal electric planes, and fallback pair/tone storage.
    /// Borrowed inputs, allocator overhead and fixed bounded temporal factor
    /// and projection scratch are excluded. No depth*F*N history is retained.
    pub maximum_workspace_bytes: usize,
}

impl<'e, 'b, 's, 'm> PreparedAngularSpectrumExposureLayer<'e, 'b, 's, 'm> {
    /// Expected irreversible volume heat [W], projected before quadratic
    /// frequency contraction when the shared positive temporal factor admits
    /// a bank. Unsupported banks retain their original exact spectral law.
    /// The authored depth quadrature and sampled/continuous spatial support
    /// are unchanged. This explicitly repeats original field reconstruction.
    pub fn integrated_heat_window(
        &self,
        quadrature: AngularSpectrumHeatVolumeQuadrature<'_>,
        channel: AngularSpectrumExposureHeatChannel,
        window: HarmonicWindow,
        limits: AngularSpectrumHeatWindowLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64, AngularSpectrumError> {
        let local = HarmonicWindow::with_profile(0.0, window.duration_s(), window.profile())
            .map_err(|_| bad("integrated heat window is unrepresentable"))?;
        let result = self.field_heat_window(
            quadrature,
            channel,
            local,
            window.start_clock_dd(),
            Dd::ONE,
            limits,
            &mut cancelled,
        )?;
        nonnegative_heat(result.to_f64())
    }

    /// Irreversible volume heat [J] over a uniform accepted source interval.
    /// The DD epoch is retained independently of duration; a short accepted
    /// step is never rounded into a large absolute clock. Volume and duration
    /// scales meet before the positive result is restored to physical units.
    pub fn integrated_heat_interval_from_clock_dd(
        &self,
        quadrature: AngularSpectrumHeatVolumeQuadrature<'_>,
        channel: AngularSpectrumExposureHeatChannel,
        start_s: Dd,
        duration_s: f64,
        limits: AngularSpectrumHeatWindowLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Dd, AngularSpectrumError> {
        let window = HarmonicWindow::new(0.0, duration_s)
            .map_err(|_| bad("integrated heat interval needs positive finite duration"))?;
        self.field_heat_window(
            quadrature,
            channel,
            window,
            start_s,
            Dd::from_f64(duration_s),
            limits,
            &mut cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn field_heat_window(
        &self,
        quadrature: AngularSpectrumHeatVolumeQuadrature<'_>,
        channel: AngularSpectrumExposureHeatChannel,
        window: HarmonicWindow,
        epoch: Dd,
        output_scale: Dd,
        limits: AngularSpectrumHeatWindowLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Dd, AngularSpectrumError> {
        poll(cancelled)?;
        if !finite_dd(epoch) || quadrature.depth_samples.len() > limits.maximum_depth_samples {
            return Err(bad("integrated heat clock or depth budget is invalid"));
        }
        let volume = admit_volume(self.anchor.grid(), self.thickness_m, quadrature, cancelled)?;
        let (first, reservoirs) = reservoir_range(
            channel,
            self.law.conductivity,
            self.law.debye_poles.len(),
            self.law.lorentz_poles.len(),
        )?;
        let active = reservoirs != 0 && volume.area_m2.to_f64() > 0.0;
        let depths = if active {
            quadrature.depth_samples.len()
        } else {
            0
        };
        let continuous = is_continuous(quadrature.transverse);
        let shifted_cell = matches!(
            quadrature.transverse,
            AngularSpectrumHeatAreaQuadrature::ContinuousPeriodicCell
        ) && self.anchor.carrier_rad_m != [0.0; 2]
            && ensemble::banks(self.mean, self.factors, self.statistics)
                .any(|(bank, proper)| !proper && bank.frequency_count() != 0);
        let aperture = match quadrature.transverse {
            AngularSpectrumHeatAreaQuadrature::ContinuousRectangle(aperture) if active => {
                Some(aperture)
            }
            AngularSpectrumHeatAreaQuadrature::ContinuousPeriodicCell if active && shifted_cell => {
                Some(AngularSpectrumAperture::try_new(
                    [0.0; 2],
                    [
                        volume.grid.width as f64 * volume.grid.spacing_m[0],
                        volume.grid.height as f64 * volume.grid.spacing_m[1],
                    ],
                )?)
            }
            _ => None,
        };
        let mut work = Work::default();
        if active {
            for (bank, proper) in ensemble::banks(self.mean, self.factors, self.statistics) {
                poll(cancelled)?;
                work.add_bank(
                    bank.frequency_count(),
                    temporal_modes(window, &bank.frequencies, proper)?,
                    proper,
                    continuous,
                    aperture.is_some() && (!shifted_cell || !proper),
                    volume.points,
                    depths,
                    reservoirs,
                )?;
            }
        }
        let rectangle_bytes = if aperture.is_some() {
            work.spatial = sum(
                work.spatial,
                sum(
                    product(&[volume.grid.width, volume.grid.width])?,
                    product(&[volume.grid.height, volume.grid.height])?,
                )?,
            )?;
            rectangle::PreparedRectangle::workspace_bytes(volume.grid)?
        } else {
            0
        };
        work.admit(
            limits,
            self.point_line_pairs,
            self.field_fft_work,
            self.workspace_bytes,
            rectangle_bytes,
            depths,
            continuous,
        )?;
        if !active {
            return Ok(Dd::ZERO);
        }
        let rectangle = aperture
            .map(|aperture| {
                rectangle::PreparedRectangle::prepare(
                    volume.grid,
                    aperture,
                    self.anchor.carrier_rad_m,
                    |i| self.anchor.envelope_wavenumber(0, i),
                    |i| self.anchor.envelope_wavenumber(1, i),
                    cancelled,
                )
            })
            .transpose()?;
        let mut total = ScaledDd::ZERO;
        for (bank, proper) in ensemble::banks(self.mean, self.factors, self.statistics) {
            poll(cancelled)?;
            if bank.frequency_count() == 0 {
                continue;
            }
            let mut mapping = allocated(bank.frequency_count(), 0usize)?;
            for (group, &omega) in bank.frequencies.iter().enumerate() {
                poll(cancelled)?;
                mapping[group] = self
                    .frequencies
                    .binary_search_by(|cache| cache.omega.total_cmp(&omega))
                    .map_err(|_| bad("original film frequency cache is missing"))?;
            }
            let rectangle = rectangle.as_ref().filter(|_| !shifted_cell || !proper);
            if temporal_modes(window, &bank.frequencies, proper)?.is_some() {
                total = total
                    + positive_bank(
                        &bank.frequencies,
                        proper,
                        window,
                        epoch,
                        quadrature,
                        volume,
                        rectangle,
                        reservoirs,
                        |r| {
                            ScaledDd::from_f64(if first + r == 0 {
                                self.law.conductivity
                            } else {
                                1.0
                            })
                        },
                        |group, r| {
                            if first + r == 0 {
                                [ScaledDd::from_f64(1.0), ScaledDd::ZERO]
                            } else {
                                parts(self.frequencies[mapping[group]].poles[first + r - 1].heat)
                            }
                        },
                        |group, depth, cancelled| {
                            let layer = &self.frequencies[mapping[group]].layer;
                            if continuous {
                                let raw = bank.frequency_modal_parts_with(
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
                                if let Some(rectangle) = rectangle {
                                    return rectangle.centered_electric(raw, cancelled);
                                }
                                electric_dd(raw, cancelled)
                            } else {
                                let mut raw = bank.frequency_modal_field_with(
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
                                bank.synthesize_frequency_field(group, &mut raw, cancelled)?;
                                electric_f64(raw, cancelled)
                            }
                        },
                        cancelled,
                    )?;
            } else {
                let generated = product(&[
                    upper_pairs(bank.frequency_count())?,
                    if proper { 1 } else { 2 },
                ])?;
                let mut terms = Vec::new();
                terms
                    .try_reserve_exact(generated)
                    .map_err(|_| bad("integrated heat fallback allocation failed"))?;
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
                if let Some(rectangle) = rectangle {
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
                        cancelled,
                    )?;
                } else {
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
                        cancelled,
                    )?;
                }
                total = total + fallback_mean(terms, window, epoch, cancelled)?;
            }
        }
        poll(cancelled)?;
        nonnegative_scaled(total * ScaledDd::from_dd(output_scale))?.physical()
    }
}

impl AngularSpectrumExposureLayerHeatLoadSpectrum<'_, '_, '_, '_, '_, '_> {
    /// Explicitly reconstruct retained original fields using separate limits.
    /// `observe_window` continues to query only the cached spectrum.
    pub fn observe_window_from_fields(
        &self,
        window: HarmonicWindow,
        limits: AngularSpectrumHeatWindowLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<f64, AngularSpectrumError> {
        self.layer
            .integrated_heat_window(self.quadrature, self.channel, window, limits, cancelled)
    }
}

#[derive(Default)]
struct Work {
    temporal: usize,
    spatial: usize,
    fallback: usize,
    bytes: usize,
}

impl Work {
    #[allow(clippy::too_many_arguments)]
    fn add_bank(
        &mut self,
        f: usize,
        modes: Option<usize>,
        proper: bool,
        continuous: bool,
        rectangle: bool,
        points: usize,
        depths: usize,
        reservoirs: usize,
    ) -> Result<(), AngularSpectrumError> {
        if f == 0 {
            return Ok(());
        }
        let mut bytes = product(&[f, core::mem::size_of::<usize>()])?;
        bytes = sum(
            bytes,
            product(&[points, core::mem::size_of::<[[Dd; 2]; 6]>()])?,
        )?;
        if let Some(modes) = modes {
            self.temporal = sum(
                self.temporal,
                product(&[
                    3,
                    depths,
                    points,
                    reservoirs,
                    f,
                    modes,
                    if continuous && !proper { 2 } else { 1 },
                ])?,
            )?;
            let spatial = if rectangle {
                sum(product(&[points, points])?, sum(points, 1)?)?
            } else {
                points
            };
            self.spatial = sum(
                self.spatial,
                product(&[3, depths, reservoirs, modes, spatial])?,
            )?;
            // Original center rotations, source-clock rotations and filtered fields.
            self.spatial = sum(
                self.spatial,
                product(&[3, depths, points, f, sum(reservoirs, 2)?])?,
            )?;
            bytes = sum(
                bytes,
                product(&[sum(f, modes)?, points, core::mem::size_of::<Electric>()])?,
            )?;
            bytes = sum(bytes, product(&[f, core::mem::size_of::<Vec<Electric>>()])?)?;
        } else {
            let pairs = upper_pairs(f)?;
            let generated = product(&[pairs, if proper { 1 } else { 2 }])?;
            self.fallback = sum(self.fallback, generated)?;
            let spatial = if rectangle {
                sum(product(&[points, points])?, 1)?
            } else {
                points
            };
            self.spatial = sum(
                self.spatial,
                sum(
                    product(&[3, depths, pairs, spatial])?,
                    product(&[pairs, reservoirs])?,
                )?,
            )?;
            let field_bytes = if rectangle {
                core::mem::size_of::<Electric>()
            } else {
                core::mem::size_of::<[C64; 6]>()
            };
            bytes = sum(bytes, product(&[f, points, field_bytes])?)?;
            bytes = sum(
                bytes,
                product(&[f, core::mem::size_of::<rectangle::Plane>()])?,
            )?;
            bytes = sum(
                bytes,
                product(&[pairs, core::mem::size_of::<FrequencyPair>()])?,
            )?;
            bytes = sum(
                bytes,
                product(&[
                    generated,
                    sum(
                        core::mem::size_of::<AngularHarmonicSpectrumTerm>(),
                        core::mem::size_of::<HarmonicProductTone>(),
                    )?,
                ])?,
            )?;
        }
        self.bytes = self.bytes.max(bytes);
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn admit(
        &self,
        limits: AngularSpectrumHeatWindowLimits,
        lines: usize,
        fft: u64,
        original_bytes: usize,
        rectangle_bytes: usize,
        depths: usize,
        continuous: bool,
    ) -> Result<(), AngularSpectrumError> {
        let fft = if continuous {
            0
        } else {
            u64::try_from(depths)
                .ok()
                .and_then(|d| fft.checked_mul(d))
                .ok_or_else(|| bad("integrated heat window FFT work overflows"))?
        };
        if product(&[lines, depths])? > limits.maximum_point_line_pairs
            || fft > limits.maximum_fft_work
            || self.temporal > limits.maximum_temporal_projection_terms
            || self.spatial > limits.maximum_spatial_products
            || self.fallback > limits.maximum_fallback_terms
            || self.fallback > MAX_PRODUCT_SPECTRUM_TERMS
            || sum(sum(original_bytes, rectangle_bytes)?, self.bytes)?
                > limits.maximum_workspace_bytes
        {
            return Err(bad(
                "integrated heat window exceeds reconstruction/projection/workspace limits",
            ));
        }
        Ok(())
    }
}

fn sum(a: usize, b: usize) -> Result<usize, AngularSpectrumError> {
    a.checked_add(b)
        .ok_or_else(|| bad("integrated heat window admission overflows"))
}
fn product(values: &[usize]) -> Result<usize, AngularSpectrumError> {
    values.iter().try_fold(1usize, |a, &b| {
        a.checked_mul(b)
            .ok_or_else(|| bad("integrated heat window admission overflows"))
    })
}
fn is_continuous(area: AngularSpectrumHeatAreaQuadrature<'_>) -> bool {
    matches!(
        area,
        AngularSpectrumHeatAreaQuadrature::ContinuousPeriodicCell
            | AngularSpectrumHeatAreaQuadrature::ContinuousRectangle(_)
    )
}
fn temporal_modes(
    window: HarmonicWindow,
    frequencies: &[f64],
    proper: bool,
) -> Result<Option<usize>, AngularSpectrumError> {
    (if proper {
        window.uniform_band_mode_count_rad_s(frequencies)
    } else {
        window.uniform_real_band_mode_count_rad_s(frequencies)
    })
    .map_err(|_| bad("integrated heat temporal band is unrepresentable"))
}
fn allocated<T: Clone>(count: usize, value: T) -> Result<Vec<T>, AngularSpectrumError> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(count)
        .map_err(|_| bad("integrated heat window allocation failed"))?;
    result.resize(count, value);
    Ok(result)
}
fn electric_dd(
    raw: Vec<[[Dd; 2]; 6]>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<Electric>, AngularSpectrumError> {
    let mut result = allocated(raw.len(), [[ScaledDd::ZERO; 2]; 3])?;
    for (point, (out, field)) in result.iter_mut().zip(raw).enumerate() {
        if point % 32 == 0 {
            poll(cancelled)?;
        }
        *out = core::array::from_fn(|axis| field[axis].map(ScaledDd::from_dd));
    }
    Ok(result)
}
fn electric_f64(
    raw: Vec<[C64; 6]>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<Electric>, AngularSpectrumError> {
    let mut result = allocated(raw.len(), [[ScaledDd::ZERO; 2]; 3])?;
    for (point, (out, field)) in result.iter_mut().zip(raw).enumerate() {
        if point % 32 == 0 {
            poll(cancelled)?;
        }
        *out = core::array::from_fn(|axis| parts(field[axis]));
    }
    Ok(result)
}
fn nonnegative_scaled(value: ScaledDd) -> Result<ScaledDd, AngularSpectrumError> {
    if !value.is_finite() || value.mantissa.to_f64() < 0.0 {
        return Err(bad(
            "integrated heat temporal/spatial contraction is negative or unrepresentable",
        ));
    }
    Ok(value)
}

fn fallback_mean(
    mut terms: Vec<AngularHarmonicSpectrumTerm>,
    window: HarmonicWindow,
    epoch: Dd,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<ScaledDd, AngularSpectrumError> {
    for term in &mut terms {
        poll(cancelled)?;
        let rotation = phase_rad_s_dd_from_clock_dd(term.angular_frequency_rad_s, epoch)
            .map_err(|_| bad("integrated heat fallback source phase is unrepresentable"))?
            .map(ScaledDd::from_dd);
        term.peak_amplitude = physical_parts(multiply(
            term.peak_amplitude.map(ScaledDd::from_dd),
            rotation,
        ))?;
    }
    let count = terms.len();
    let spectrum = PreparedHarmonicProductSpectrum::from_angular_terms_cancellable(
        terms,
        HarmonicSpectrumTermLimits {
            maximum_terms: count,
            maximum_workspace_bytes: product(&[
                count,
                core::mem::size_of::<HarmonicProductTone>(),
            ])?,
        },
        &mut *cancelled,
    )
    .map_err(spectrum_error)?;
    let mut total = ScaledDd::from_dd(spectrum.dc_dd());
    for tone in spectrum.tones() {
        poll(cancelled)?;
        let value = window
            .product_kernel_dd(tone.angular_frequency_dd(), Dd::ZERO)
            .and_then(|kernel| {
                kernel.mean_real_product_parts_dd(tone.peak_amplitude_dd(), [Dd::ONE, Dd::ZERO])
            })
            .map_err(|_| bad("integrated heat fallback window is unrepresentable"))?;
        total = total + ScaledDd::from_dd(value);
    }
    poll(cancelled)?;
    nonnegative_scaled(total)
}

#[allow(clippy::too_many_arguments)]
fn positive_bank<F: FnMut() -> bool>(
    frequencies: &[f64],
    proper: bool,
    window: HarmonicWindow,
    epoch: Dd,
    quadrature: AngularSpectrumHeatVolumeQuadrature<'_>,
    volume: VolumeAdmission,
    rectangle: Option<&rectangle::PreparedRectangle>,
    reservoirs: usize,
    reservoir_scale: impl Fn(usize) -> ScaledDd,
    transfer: impl Fn(usize, usize) -> [ScaledDd; 2],
    mut field: impl FnMut(usize, f64, &mut F) -> Result<Vec<Electric>, AngularSpectrumError>,
    cancelled: &mut F,
) -> Result<ScaledDd, AngularSpectrumError> {
    poll(cancelled)?;
    let mut scratch = UniformHarmonicBandWorkspace::new();
    let mut negative = [0.0; MAX_UNIFORM_HARMONIC_BAND_TONES];
    if frequencies.len() > negative.len() {
        return Err(bad("integrated heat temporal band exceeds fixed scratch"));
    }
    for (out, &omega) in negative.iter_mut().zip(frequencies) {
        *out = -omega;
    }
    // Negative native clocks keep EM's original exp(+i omega t) fields in
    // their ORIGINAL complex spatial convention, including continuous real
    // quadratures. No complex conjugation of the Fourier envelope is hidden.
    let factor = if proper {
        window.prepare_uniform_band_rad_s_in(&negative[..frequencies.len()], &mut scratch)
    } else {
        window.prepare_uniform_real_band_rad_s_in(&negative[..frequencies.len()], &mut scratch)
    }
    .map_err(|_| bad("integrated heat temporal preparation is unrepresentable"))?
    .ok_or_else(|| bad("integrated heat positive factor was not admitted"))?;
    let mut rotations = [[ScaledDd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_TONES];
    let reference = if proper {
        Dd::from_f64(frequencies[0])
    } else {
        Dd::ZERO
    };
    for (rotation, &omega) in rotations.iter_mut().zip(frequencies) {
        poll(cancelled)?;
        *rotation = phase_rad_s_dd_from_clock_dd(reference - Dd::from_f64(omega), epoch)
            .map_err(|_| bad("integrated heat original source clock is unrepresentable"))?
            .map(ScaledDd::from_dd);
    }
    let mut fields = Vec::new();
    fields
        .try_reserve_exact(frequencies.len())
        .map_err(|_| bad("integrated heat frequency planes allocation failed"))?;
    let mut modes = allocated(
        product(&[volume.points, factor.mode_count()])?,
        [[ScaledDd::ZERO; 2]; 3],
    )?;
    let continuous = is_continuous(quadrature.transverse);
    let mut total = ScaledDd::ZERO;
    for sample in quadrature.depth_samples {
        poll(cancelled)?;
        fields.clear();
        for (group, &rotation) in rotations[..frequencies.len()].iter().enumerate() {
            let mut plane = field(group, sample.depth_m, cancelled)?;
            if plane.len() != volume.points {
                return Err(bad("integrated heat field differs from original grid"));
            }
            for (point, electric) in plane.iter_mut().enumerate() {
                if point % 32 == 0 {
                    poll(cancelled)?;
                }
                for component in electric {
                    *component = multiply(*component, rotation);
                }
            }
            fields.push(plane);
        }
        for reservoir in 0..reservoirs {
            poll(cancelled)?;
            let scale = reservoir_scale(reservoir);
            if scale.is_zero() {
                continue;
            }
            let mut filters = [[ScaledDd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_TONES];
            for (group, out) in filters[..frequencies.len()].iter_mut().enumerate() {
                *out = transfer(group, reservoir);
            }
            for point in 0..volume.points {
                if point % 16 == 0 {
                    poll(cancelled)?;
                }
                for axis in 0..3 {
                    let mut filtered = [[ScaledDd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_TONES];
                    let mut exponent = i32::MIN;
                    for group in 0..frequencies.len() {
                        filtered[group] = multiply(fields[group][point][axis], filters[group]);
                        if !finite_parts(filtered[group]) {
                            return Err(bad("integrated heat filtered field is unrepresentable"));
                        }
                        for component in filtered[group] {
                            if !component.is_zero() {
                                exponent = exponent.max(component.exponent);
                            }
                        }
                    }
                    if exponent == i32::MIN {
                        exponent = 0;
                    }
                    let mut source = [[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_TONES];
                    for group in 0..frequencies.len() {
                        source[group] = filtered[group].map(|v| {
                            Dd::from_pair(
                                scale_by_two(v.mantissa.hi, v.exponent - exponent),
                                scale_by_two(v.mantissa.lo, v.exponent - exponent),
                            )
                        });
                    }
                    let projected = project_spatial_parts(
                        &factor,
                        &mut source[..frequencies.len()],
                        proper,
                        continuous,
                    )?;
                    for mode in 0..factor.mode_count() {
                        modes[mode * volume.points + point][axis] =
                            projected[mode].map(|v| ScaledDd::normalized(v, exponent));
                    }
                }
            }
            for mode in 0..factor.mode_count() {
                poll(cancelled)?;
                let spatial = spatial_square(
                    &modes[mode * volume.points..(mode + 1) * volume.points],
                    proper,
                    quadrature.transverse,
                    volume,
                    rectangle,
                    cancelled,
                )?;
                let weight = scale
                    * ScaledDd::from_f64(sample.weight_m)
                    * ScaledDd::from_f64((2 * mode + 1) as f64)
                    * ScaledDd::from_f64(if proper { 0.5 } else { 1.0 });
                total = total + nonnegative_scaled(spatial)? * weight;
            }
        }
    }
    poll(cancelled)?;
    nonnegative_scaled(total)
}

fn project_spatial_parts(
    factor: &UniformHarmonicBandRef<'_>,
    source: &mut [[Dd; 2]],
    proper: bool,
    continuous: bool,
) -> Result<[[Dd; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES], AngularSpectrumError> {
    let mut projected = [[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES];
    factor
        .project_parts_dd(source, &mut projected[..factor.mode_count()])
        .map_err(|_| bad("integrated heat temporal projection is unrepresentable"))?;
    if !proper && continuous {
        // A real temporal factor projects one quadrature per parity. Spatial
        // Fourier coefficients need BOTH quadratures: apply it also to i*C.
        // For odd l the factor retains the imaginary component before the
        // common (-i)^l phase, so its reconstruction has the opposite ordering.
        for value in source.iter_mut() {
            *value = [-value[1], value[0]];
        }
        let mut second = [[Dd::ZERO; 2]; MAX_UNIFORM_HARMONIC_BAND_MODES];
        factor
            .project_parts_dd(source, &mut second[..factor.mode_count()])
            .map_err(|_| bad("integrated heat spatial quadrature projection is unrepresentable"))?;
        for mode in 0..factor.mode_count() {
            projected[mode] = if mode % 2 == 0 {
                [projected[mode][0], -second[mode][0]]
            } else {
                [second[mode][1], projected[mode][1]]
            };
        }
    }
    for (mode, value) in projected[..factor.mode_count()].iter_mut().enumerate() {
        let [re, im] = *value;
        *value = match mode % 4 {
            0 => [re, im],
            1 => [im, -re],
            2 => [-re, -im],
            _ => [-im, re],
        };
    }
    Ok(projected)
}

fn spatial_square(
    field: &[Electric],
    proper: bool,
    area: AngularSpectrumHeatAreaQuadrature<'_>,
    volume: VolumeAdmission,
    rectangle: Option<&rectangle::PreparedRectangle>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<ScaledDd, AngularSpectrumError> {
    if let Some(rectangle) = rectangle {
        return Ok(rectangle.square_centered(field, proper, cancelled)?
            * ScaledDd::from_dd(volume.area_m2));
    }
    let continuous = is_continuous(area);
    let half = ScaledDd::from_f64(0.5);
    let mut result = ScaledDd::ZERO;
    for point in 0..volume.points {
        if point % 32 == 0 {
            poll(cancelled)?;
        }
        let weight = ScaledDd::from_dd(area_weight(area, volume, point));
        if weight.is_zero() {
            continue;
        }
        let partner = if continuous && !proper {
            opposite_mode(volume.grid, point)
        } else {
            None
        };
        if partner.is_some_and(|other| other < point) {
            continue;
        }
        for axis in 0..3 {
            let value = field[point][axis];
            let square = if proper {
                multiply_conjugate(value, value)[0]
            } else if !continuous || partner == Some(point) {
                value[0] * value[0]
            } else if let Some(other) = partner {
                let real = value[0] + field[other][axis][0];
                let imaginary = value[1] - field[other][axis][1];
                half * (real * real + imaginary * imaginary)
            } else {
                half * multiply_conjugate(value, value)[0]
            };
            result = result + square * weight;
        }
    }
    poll(cancelled)?;
    Ok(result)
}

#[cfg(test)]
mod tests;
