//! Prepared finite-shutter interference between explicitly analyzed arms.
//!
//! Frequency grouping and real/proper source laws belong to the exposure
//! owners. This observer keeps their DD amplitudes through arm recombination,
//! then uses the original harmonic-window kernels for the real square.

use super::super::coherence::correlation::{
    ComplexDd, admit_analyzers, admit_arm_count, analyze, conjugate, lift, multiply, nonnegative,
    plus, rounded, synthesize_electric_parts,
};
use super::*;

/// One exterior observation and its dimensionless Hermitian analyzer `w† E`.
/// Positions use the original periodic envelope and Cartesian carrier. All
/// six original illumination sides are supported; no path or splitter gain
/// is inferred from the observation position.
#[derive(Clone, Copy, Debug)]
pub struct AngularSpectrumExposureInterferenceArm {
    pub side: AngularSpectrumIlluminationSide,
    pub distance_from_face_m: f64,
    pub position_m: [f64; 2],
    pub analyzer: [C64; 3],
}

/// Admission for preparation and every later reading. Mode work conservatively
/// counts A*N*sum(L_r); equal observation planes share their reconstructed
/// fields. Temporal work counts sum(F_r*(F_r+1)/2), including the mean bank.
/// Storage includes retained amplitudes/kernels, scratch, grouping indices and
/// one streamed original modal field. Borrowed source owners, arms, gains,
/// delays and allocator overhead are excluded. No FFT is performed here.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AngularSpectrumExposureInterferenceLimits {
    pub maximum_mode_evaluations: usize,
    pub maximum_frequency_pairs: usize,
    pub maximum_workspace_bytes: usize,
}

enum Kernels {
    Real(Vec<HarmonicProduct>),
    Proper(Vec<ProperHarmonicProduct>),
}

struct Bank {
    frequencies: Vec<f64>,
    // Frequency-major, then original arm order. Different banks never share
    // products, even when they contain exactly equal frequencies.
    amplitudes: Vec<ComplexDd>,
    kernels: Kernels,
}

/// Original spectral fields sampled once for a fixed physical shutter.
///
/// Each reading costs O(A*sum F_r + sum F_r²), with no allocations, FFTs or
/// Maxwell/material reconstruction. The mutable receiver reuses an admitted
/// frequency scratch vector. Cancellation leaves the preparation reusable.
/// Changing the source, observation, analyzer or shutter requires preparation.
pub struct PreparedAngularSpectrumExposureInterferometer {
    arm_count: usize,
    window: HarmonicWindow,
    banks: Vec<Bank>,
    scratch: Vec<ComplexDd>,
}

impl PreparedAngularSpectrumExposureInterferometer {
    pub fn arm_count(&self) -> usize {
        self.arm_count
    }
    pub fn window(&self) -> HarmonicWindow {
        self.window
    }

    /// Sum of distinct-frequency counts within the original source banks.
    /// Equal frequencies in different uncorrelated factors stay separate.
    pub fn frequency_group_count(&self) -> usize {
        self.banks.iter().map(|bank| bank.frequencies.len()).sum()
    }

    /// Finite-shutter electric mean square in V²/m² after arm recombination.
    ///
    /// At angular frequency omega, arm a contributes
    /// `gain[a] * exp(-i*omega*delay_s[a]) * w[a]† E[a,omega]`.
    /// Positive delay means the physical waveform is evaluated at t-delay;
    /// finite negative relative delays are also admitted. Gains are authored,
    /// frequency-independent complex amplitude multipliers, not intensities.
    /// Delays act on each actual frequency, never on an inferred carrier.
    ///
    /// The deterministic mean and real factors retain sum and difference
    /// tones. Proper factors retain their original difference terms only.
    /// Banks add second moments after their own coherent frequency/arm sums;
    /// frequency equality alone creates no cross-factor coherence. No PSD,
    /// optical-cycle averaging, detector watts or back-action is inferred.
    pub fn electric_mean_square_v2_m2(
        &mut self,
        gains: &[C64],
        delays_s: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<f64, AngularSpectrumError> {
        poll(&mut cancelled)?;
        if gains.len() != self.arm_count
            || delays_s.len() != self.arm_count
            || gains.iter().any(|&gain| !finite(gain))
            || delays_s.iter().any(|delay| !delay.is_finite())
        {
            return Err(bad(
                "exposure interference needs one finite gain and delay per arm",
            ));
        }
        let mut total = Dd::ZERO;
        for bank in &self.banks {
            poll(&mut cancelled)?;
            for (frequency, (&omega, amplitudes)) in bank
                .frequencies
                .iter()
                .zip(bank.amplitudes.chunks_exact(self.arm_count))
                .enumerate()
            {
                poll(&mut cancelled)?;
                let mut field = [Dd::ZERO; 2];
                for (arm, ((&amplitude, &gain), &delay)) in
                    amplitudes.iter().zip(gains).zip(delays_s).enumerate()
                {
                    if arm % 32 == 0 {
                        poll(&mut cancelled)?;
                    }
                    let phase = Dd::from_f64(omega) * Dd::from_f64(delay);
                    if !phase.hi.is_finite() || !phase.lo.is_finite() {
                        return Err(bad("interferometer delay phase is unrepresentable"));
                    }
                    let (sine, cosine) = det::sin_cos_dd(phase);
                    field = plus(
                        field,
                        multiply(lift(gain), multiply(amplitude, [cosine, -sine])),
                    );
                }
                // Source phasors use +i omega t; the original temporal owner
                // uses -i omega t. Conjugate AFTER gain/delay recombination.
                self.scratch[frequency] = conjugate(field);
            }
            let amplitudes = &self.scratch[..bank.frequencies.len()];
            let value = match &bank.kernels {
                Kernels::Real(kernels) => contract(
                    kernels,
                    amplitudes,
                    HarmonicProduct::mean_real_product_parts_dd,
                    &mut cancelled,
                )?,
                Kernels::Proper(kernels) => contract(
                    kernels,
                    amplitudes,
                    ProperHarmonicProduct::mean_real_product_parts_dd,
                    &mut cancelled,
                )?,
            };
            nonnegative(value)?;
            total = total + value;
        }
        let result = nonnegative(total)?;
        poll(&mut cancelled)?;
        Ok(result)
    }
}

// Symmetry halves both cached kernels and repeated quadratic work. Equal-
// frequency source/arm fields have already combined before any square.
fn contract<P: Copy>(
    kernels: &[P],
    amplitudes: &[ComplexDd],
    product: impl Fn(P, ComplexDd, ComplexDd) -> Result<Dd, fs_fft::harmonic_window::WindowError>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Dd, AngularSpectrumError> {
    let mut sum = Dd::ZERO;
    let mut pair = 0;
    for (i, &a) in amplitudes.iter().enumerate() {
        for (j, &b) in amplitudes.iter().enumerate().skip(i) {
            if pair % 32 == 0 {
                poll(cancelled)?;
            }
            let value = product(kernels[pair], a, b)
                .map_err(|_| bad("interferometer temporal product is unrepresentable"))?;
            sum = sum
                + if i == j {
                    value
                } else {
                    value * Dd::from_f64(2.0)
                };
            pair += 1;
        }
    }
    Ok(sum)
}

struct Sources<'e, 'b, 's, 'm> {
    anchor: &'b PlanarAngularSpectrum<'s, 'm>,
    mean: Option<&'e AngularSpectrumExposure<'b, 's, 'm>>,
    factors: &'e [AngularSpectrumExposure<'b, 's, 'm>],
    statistics: Option<&'e [AngularSpectrumFactorStatistics]>,
}

impl<'e, 'b, 's, 'm> Sources<'e, 'b, 's, 'm> {
    fn banks(&self) -> impl Iterator<Item = (&'e AngularSpectrumExposure<'b, 's, 'm>, bool)> {
        ensemble::banks(self.mean, self.factors, self.statistics)
    }

    fn requirements(
        &self,
        arms: usize,
    ) -> Result<AngularSpectrumExposureInterferenceLimits, AngularSpectrumError> {
        admit_arm_count(arms)?;
        let overflow = || bad("exposure interferometer work or storage overflow");
        let mut modes = 0usize;
        let mut pairs = 0usize;
        let mut bytes = core::mem::size_of::<PreparedAngularSpectrumExposureInterferometer>();
        let mut max_frequencies = 0;
        for (bank, proper) in self.banks() {
            let f = bank.frequency_count();
            max_frequencies = max_frequencies.max(f);
            modes = bank
                .line_count()
                .checked_mul(self.anchor.modes.len())
                .and_then(|count| count.checked_mul(arms))
                .and_then(|count| modes.checked_add(count))
                .ok_or_else(overflow)?;
            let count = f.checked_mul(f + 1).ok_or_else(overflow)? / 2;
            pairs = pairs.checked_add(count).ok_or_else(overflow)?;
            let kernel_size = if proper {
                core::mem::size_of::<ProperHarmonicProduct>()
            } else {
                core::mem::size_of::<HarmonicProduct>()
            };
            let amplitude_bytes = f
                .checked_mul(arms)
                .and_then(|count| count.checked_mul(core::mem::size_of::<ComplexDd>()))
                .ok_or_else(overflow)?;
            for extra in [
                core::mem::size_of::<Bank>(),
                amplitude_bytes,
                f.checked_mul(core::mem::size_of::<f64>())
                    .ok_or_else(overflow)?,
                count.checked_mul(kernel_size).ok_or_else(overflow)?,
            ] {
                bytes = bytes.checked_add(extra).ok_or_else(overflow)?;
            }
        }
        let modal_bytes = if modes == 0 {
            0
        } else {
            self.anchor
                .modes
                .len()
                .checked_mul(6 * (core::mem::size_of::<ComplexDd>() + core::mem::size_of::<C64>()))
                .ok_or_else(overflow)?
        };
        for extra in [
            modal_bytes,
            arms.checked_mul(core::mem::size_of::<usize>())
                .ok_or_else(overflow)?,
            max_frequencies
                .checked_mul(core::mem::size_of::<ComplexDd>())
                .ok_or_else(overflow)?,
        ] {
            bytes = bytes.checked_add(extra).ok_or_else(overflow)?;
        }
        Ok(AngularSpectrumExposureInterferenceLimits {
            maximum_mode_evaluations: modes,
            maximum_frequency_pairs: pairs,
            maximum_workspace_bytes: bytes,
        })
    }

    fn prepare(
        &self,
        arms: &[AngularSpectrumExposureInterferenceArm],
        window: HarmonicWindow,
        limits: AngularSpectrumExposureInterferenceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedAngularSpectrumExposureInterferometer, AngularSpectrumError> {
        poll(cancelled)?;
        let needed = self.requirements(arms.len())?;
        if needed.maximum_mode_evaluations > limits.maximum_mode_evaluations
            || needed.maximum_frequency_pairs > limits.maximum_frequency_pairs
            || needed.maximum_workspace_bytes > limits.maximum_workspace_bytes
        {
            return Err(bad(
                "exposure interferometer exceeds its work or storage budget",
            ));
        }
        for arm in arms {
            poll(cancelled)?;
            if !arm.distance_from_face_m.is_finite()
                || arm.distance_from_face_m < 0.0
                || arm.position_m.iter().any(|value| !value.is_finite())
            {
                return Err(bad(
                    "interference arms need finite coordinates and nonnegative depth",
                ));
            }
            admit_analyzers([arm.analyzer; 2])?;
        }
        // Stream one original modal plane per side/depth. Spatially separate
        // arms in the same plane share all original field reconstructions.
        let mut order = reserved(arms.len())?;
        order.extend(0..arms.len());
        order.sort_unstable_by(|&a, &b| {
            (arms[a].side as u8)
                .cmp(&(arms[b].side as u8))
                .then_with(|| {
                    arms[a]
                        .distance_from_face_m
                        .total_cmp(&arms[b].distance_from_face_m)
                })
                .then_with(|| a.cmp(&b))
        });
        let mut banks = reserved(self.factors.len() + usize::from(self.mean.is_some()))?;
        let mut max_frequencies = 0;
        for (source, proper) in self.banks() {
            let f = source.frequency_count();
            max_frequencies = max_frequencies.max(f);
            let mut frequencies = reserved(f)?;
            frequencies.extend_from_slice(&source.frequencies);
            let kernels = if proper {
                Kernels::Proper(prepare_kernels(
                    &frequencies,
                    |a, b| window.proper_product_kernel(a, b),
                    cancelled,
                )?)
            } else {
                Kernels::Real(prepare_kernels(
                    &frequencies,
                    |a, b| window.product_kernel(a, b),
                    cancelled,
                )?)
            };
            let mut amplitudes = reserved(f * arms.len())?;
            for index in 0..f * arms.len() {
                if index % 256 == 0 {
                    poll(cancelled)?;
                }
                amplitudes.push([Dd::ZERO; 2]);
            }
            for group in 0..f {
                let mut begin = 0;
                while begin < order.len() {
                    let plane = arms[order[begin]];
                    let parts = source.frequency_modal_parts(
                        group,
                        ExposureSide::TwoSided(plane.side),
                        plane.distance_from_face_m,
                        cancelled,
                    )?;
                    let mut end = begin;
                    while end < order.len() {
                        let index = order[end];
                        let arm = arms[index];
                        if arm.side != plane.side
                            || arm.distance_from_face_m != plane.distance_from_face_m
                        {
                            break;
                        }
                        let [field] = synthesize_electric_parts(
                            source.anchor,
                            [arm.position_m],
                            |mode, _| Ok([parts[mode][0], parts[mode][1], parts[mode][2]]),
                            cancelled,
                        )?;
                        let amplitude = analyze(field, arm.analyzer);
                        rounded(amplitude)?;
                        amplitudes[group * arms.len() + index] = amplitude;
                        end += 1;
                    }
                    begin = end;
                }
            }
            banks.push(Bank {
                frequencies,
                amplitudes,
                kernels,
            });
        }
        let mut scratch = reserved(max_frequencies)?;
        scratch.resize(max_frequencies, [Dd::ZERO; 2]);
        poll(cancelled)?;
        Ok(PreparedAngularSpectrumExposureInterferometer {
            arm_count: arms.len(),
            window,
            banks,
            scratch,
        })
    }
}

fn reserved<T>(count: usize) -> Result<Vec<T>, AngularSpectrumError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| bad("exposure interferometer allocation failed"))?;
    Ok(values)
}

fn prepare_kernels<P>(
    frequencies: &[f64],
    kernel: impl Fn(f64, f64) -> Result<P, fs_fft::harmonic_window::WindowError>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<P>, AngularSpectrumError> {
    let mut kernels = reserved(frequencies.len() * (frequencies.len() + 1) / 2)?;
    for (i, &a) in frequencies.iter().enumerate() {
        for &b in &frequencies[i..] {
            poll(cancelled)?;
            kernels.push(
                kernel(a, b)
                    .map_err(|_| bad("interferometer temporal kernel is unrepresentable"))?,
            );
        }
    }
    Ok(kernels)
}

impl AngularSpectrumExposure<'_, '_, '_> {
    /// Conservative complete bounds for 1..=4096 observation arms.
    pub fn interferometer_requirements(
        &self,
        arms: usize,
    ) -> Result<AngularSpectrumExposureInterferenceLimits, AngularSpectrumError> {
        Sources {
            anchor: self.anchor,
            mean: Some(self),
            factors: &[],
            statistics: None,
        }
        .requirements(arms)
    }

    /// Cache the original phase-related lines for a fixed finite shutter.
    pub fn prepare_interferometer(
        &self,
        arms: &[AngularSpectrumExposureInterferenceArm],
        window: HarmonicWindow,
        limits: AngularSpectrumExposureInterferenceLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedAngularSpectrumExposureInterferometer, AngularSpectrumError> {
        Sources {
            anchor: self.anchor,
            mean: Some(self),
            factors: &[],
            statistics: None,
        }
        .prepare(arms, window, limits, &mut cancelled)
    }
}

impl AngularSpectrumExposureEnsemble<'_, '_, '_, '_> {
    /// Conservative bounds over the mean and every original real/proper factor.
    pub fn interferometer_requirements(
        &self,
        arms: usize,
    ) -> Result<AngularSpectrumExposureInterferenceLimits, AngularSpectrumError> {
        Sources {
            anchor: self.anchor,
            mean: self.mean,
            factors: self.factors,
            statistics: self.statistics,
        }
        .requirements(arms)
    }

    /// Prepare all original second-moment banks without cross-factor products.
    pub fn prepare_interferometer(
        &self,
        arms: &[AngularSpectrumExposureInterferenceArm],
        window: HarmonicWindow,
        limits: AngularSpectrumExposureInterferenceLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedAngularSpectrumExposureInterferometer, AngularSpectrumError> {
        Sources {
            anchor: self.anchor,
            mean: self.mean,
            factors: self.factors,
            statistics: self.statistics,
        }
        .prepare(arms, window, limits, &mut cancelled)
    }
}
