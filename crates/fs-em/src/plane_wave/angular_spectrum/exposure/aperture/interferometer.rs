//! Reusable interference measurements through a physical rectangular aperture.
//! Original Maxwell fields, source statistics and shutter kernels remain the
//! owners of propagation, coherence and integration.

use super::super::super::coherence::correlation::admit_arm_count;
use super::gram::{GramWorkspace, clear};
use super::*;

mod temporal;

/// One parallel copy of the original field in the common observed half-space.
///
/// The copy samples `F(x + translation_x, y + translation_y, depth, t - delay)`.
/// All arms keep the same Cartesian orientation and observation medium.
/// Depth is measured outwards from the selected stack face. No splitter,
/// remapping of propagation directions, impedance conversion or throughput is
/// inferred; complex amplitude gains and relative delays are supplied per read.
#[derive(Clone, Copy, Debug)]
pub struct AngularSpectrumExposureApertureArm {
    pub distance_from_face_m: f64,
    pub translation_m: [f64; 2],
}

/// Bounds for preparation, retained fields and repeated detector work.
///
/// With A arms, N = Nx*Ny modes and F_r distinct frequencies in source bank r:
/// source storage scales as A*N*sum(F_r), while the separable Gram application
/// costs sum(N*F_r*(Nx+Ny+F_r)). It never builds an (N*F)² table. Logical
/// workspace includes retained kernels and DD fields, reusable work planes,
/// arm ordering and one original modal reconstruction. Borrowed source owners,
/// caller arrays and allocator overhead are excluded.
/// Window-specific requirements additionally reserve bounded positive temporal
/// projection where admitted, with work N*M*(F+Nx+Ny+1), M<=17, plus N*F
/// original-clock rotations for proper factors and A*(F+1) additional proper
/// arm carrier-phase/composition terms. The larger
/// work bound is retained per selected bank; tighter legacy limits keep the
/// original separable pair-kernel path.
#[derive(Clone, Copy, Debug)]
pub struct AngularSpectrumApertureInterferenceLimits {
    /// Conservative original mode/line reconstruction work: A*N*sum(L_r).
    pub maximum_mode_evaluations: usize,
    /// Six-component complex DD field entries: A*N*sum(F_r).
    pub maximum_source_field_entries: usize,
    /// Scalar spatial/time kernel applications: sum(N*F_r*(Nx+Ny+F_r)).
    /// Each application handles six field components and a fixed parity count.
    pub maximum_gram_products: usize,
    pub maximum_workspace_bytes: usize,
}

enum Temporal {
    Real(Vec<RealKernel>),
    Proper(Vec<ProperKernel>),
}
struct Bank {
    frequencies: Vec<f64>,
    // Frequency-major, arm-major, then original Cartesian modal order.
    fields: Vec<FieldParts>,
    temporal: Temporal,
    positive: Option<temporal::Positive>,
}

/// Fixed aperture/shutter observation of coherently recombined original fields.
///
/// Preparation retains the six original E/H components, including evanescent
/// and counterpropagating fields. Repeated reads allocate nothing and perform
/// no FFT, material lookup or Maxwell reconstruction. Mean and real factors
/// retain carrier sums and differences; proper factors retain differences.
/// Different banks always add moments independently. Cancellation leaves the
/// immutable field basis reusable; the next read overwrites all scratch.
/// Original fields remain DD phasors: arbitrary combined complex-gain, phase
/// and spatial nulls have no universal relative-error guarantee.
///
/// The returned normal power is signed Cartesian z-directed Poynting watts.
/// This ideal observation does not infer absorption, detector back-action,
/// passivity of authored gains, or a positive optical intensity from |E|².
pub struct PreparedAngularSpectrumApertureInterferometer {
    arm_count: usize,
    modes: usize,
    side: AngularSpectrumIlluminationSide,
    aperture: AngularSpectrumAperture,
    window: HarmonicWindow,
    banks: Vec<Bank>,
    spatial_real: Option<Spatial<RealKernel>>,
    spatial_proper: Option<Spatial<ProperKernel>>,
    gram: GramWorkspace,
    projected: Vec<FieldParts>,
}

impl PreparedAngularSpectrumApertureInterferometer {
    pub fn arm_count(&self) -> usize {
        self.arm_count
    }
    pub fn side(&self) -> AngularSpectrumIlluminationSide {
        self.side
    }
    pub fn aperture(&self) -> AngularSpectrumAperture {
        self.aperture
    }
    pub fn window(&self) -> HarmonicWindow {
        self.window
    }

    /// All 21 real E/H second moments, vector Poynting flux and normal watts.
    ///
    /// Arm a contributes `gain[a] * exp(-i*omega*delay_s[a])` at each actual
    /// angular frequency omega. Positive delay means t-delay; finite negative
    /// relative delays are also valid. Gains multiply E and H together before
    /// any products. Equal frequencies are coherent only within their original
    /// bank; the preparation introduces no cross-factor coherence.
    pub fn measure_moments(
        &mut self,
        gains: &[C64],
        delays_s: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureApertureMoments, AngularSpectrumError> {
        poll(&mut cancelled)?;
        if gains.len() != self.arm_count
            || delays_s.len() != self.arm_count
            || gains.iter().any(|&gain| !finite(gain))
            || delays_s.iter().any(|delay| !delay.is_finite())
        {
            return Err(bad(
                "aperture interference needs one finite gain and delay per arm",
            ));
        }
        let mut moments = [Dd::ZERO; 21];
        for bank in &self.banks {
            let f = bank.frequencies.len();
            if f == 0 {
                continue;
            }
            let fields = &mut self.gram.residual[..self.modes * f];
            combine(
                bank,
                self.arm_count,
                self.modes,
                gains,
                delays_s,
                fields,
                &mut cancelled,
            )?;
            if let Some(positive) = &bank.positive {
                positive.accumulate(
                    &mut self.gram,
                    self.spatial_real.as_ref(),
                    self.spatial_proper.as_ref(),
                    self.modes,
                    &mut self.projected,
                    &mut moments,
                    &mut cancelled,
                )?;
                continue;
            }
            match &bank.temporal {
                Temporal::Real(temporal) => {
                    self.gram.real(
                        self.spatial_real.as_ref().expect("prepared real bank"),
                        temporal,
                        f,
                        &mut cancelled,
                    )?;
                    // The real parity graph uses the input as a work plane.
                    combine(
                        bank,
                        self.arm_count,
                        self.modes,
                        gains,
                        delays_s,
                        &mut self.gram.residual[..self.modes * f],
                        &mut cancelled,
                    )?;
                }
                Temporal::Proper(temporal) => self.gram.proper(
                    self.spatial_proper.as_ref().expect("prepared proper bank"),
                    temporal,
                    f,
                    &mut cancelled,
                )?,
            }
            for (mode, (field, covector)) in self.gram.residual[..self.modes * f]
                .iter()
                .zip(&self.gram.covectors[..self.modes * f])
                .enumerate()
            {
                if mode % 32 == 0 {
                    poll(&mut cancelled)?;
                }
                let mut pair = 0;
                for i in 0..6 {
                    for j in i..6 {
                        moments[pair] = moments[pair]
                            + field[i][0] * covector[j][0]
                            + field[i][1] * covector[j][1];
                        pair += 1;
                    }
                }
            }
        }
        poll(&mut cancelled)?;
        polarization::Full(moments).finish(self.aperture)
    }
}

// Keep arm cancellation in DD before applying any spatial or temporal Gram
// operator. The same finite gain/delay acts on the original E and H phasors.
fn combine(
    bank: &Bank,
    arms: usize,
    modes: usize,
    gains: &[C64],
    delays: &[f64],
    output: &mut [FieldParts],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), AngularSpectrumError> {
    clear(output, cancelled)?;
    // A proper bank's common carrier phase cancels from every second moment.
    // Remove the first active arm's carrier delay before composing the small
    // frequency offsets. Independently rotating each bright optical phasor
    // otherwise perturbs a resolved dark envelope before it is projected.
    // Real banks keep their absolute phase. An unrepresentable relative
    // carrier phase retains the original independent phase path.
    let mut common = if bank.positive.is_some() && matches!(&bank.temporal, Temporal::Proper(_)) {
        gains
            .iter()
            .position(|gain| gain.re != 0.0 || gain.im != 0.0)
            .map(|reference| (bank.frequencies[0], delays[reference]))
    } else {
        None
    };
    if let Some((omega, reference_delay)) = common {
        for &delay in delays {
            poll(cancelled)?;
            let relative = Dd::from_f64(delay) - Dd::from_f64(reference_delay);
            if !finite_dd(relative) || !finite_dd(Dd::from_f64(omega) * relative) {
                common = None;
                break;
            }
        }
    }
    for arm in 0..arms {
        poll(cancelled)?;
        let carrier = if let Some((omega, reference_delay)) = common {
            Some(
                fs_fft::harmonic_window::phase_rad_s_dd_from_clock_dd(
                    Dd::from_f64(omega),
                    Dd::from_f64(delays[arm]) - Dd::from_f64(reference_delay),
                )
                .map_err(|_| bad("aperture arm carrier phase is unrepresentable"))?,
            )
        } else {
            None
        };
        for (frequency, &omega) in bank.frequencies.iter().enumerate() {
            poll(cancelled)?;
            let [cosine, negative_sine] =
                if let (Some((reference, _)), Some([cr, ci])) = (common, carrier) {
                    let [dr, di] = fs_fft::harmonic_window::phase_rad_s_dd_from_clock_dd(
                        Dd::from_f64(omega) - Dd::from_f64(reference),
                        Dd::from_f64(delays[arm]),
                    )
                    .map_err(|_| bad("aperture arm envelope phase is unrepresentable"))?;
                    [cr * dr - ci * di, cr * di + ci * dr]
                } else {
                    let phase = Dd::from_f64(omega) * Dd::from_f64(delays[arm]);
                    if !finite_dd(phase) {
                        return Err(bad("aperture arm delay phase is unrepresentable"));
                    }
                    let (sine, cosine) = fs_math::det::sin_cos_dd(phase);
                    [cosine, -sine]
                };
            let [gr, gi] = [gains[arm].re, gains[arm].im].map(Dd::from_f64);
            let (real, imaginary) = (
                gr * cosine - gi * negative_sine,
                gi * cosine + gr * negative_sine,
            );
            let begin = (frequency * arms + arm) * modes;
            for (mode, input) in bank.fields[begin..begin + modes].iter().enumerate() {
                if mode % 32 == 0 {
                    poll(cancelled)?;
                }
                for (sum, &[r, i]) in output[frequency * modes + mode].iter_mut().zip(input) {
                    sum[0] = sum[0] + real * r - imaginary * i;
                    sum[1] = sum[1] + real * i + imaginary * r;
                }
            }
        }
    }
    for tile in output.chunks(128) {
        poll(cancelled)?;
        if tile
            .iter()
            .flatten()
            .flatten()
            .any(|&part| !finite_dd(part))
        {
            return Err(bad("recombined aperture field is unrepresentable"));
        }
    }
    Ok(())
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
    ) -> Result<AngularSpectrumApertureInterferenceLimits, AngularSpectrumError> {
        admit_arm_count(arms)?;
        let overflow = || bad("aperture interferometer work or storage overflow");
        let (nx, ny) = (self.anchor.grid.width, self.anchor.grid.height);
        let n = self.anchor.modes.len();
        let mut mode_work = 0usize;
        let mut entries = 0usize;
        let mut products = 0usize;
        let mut bytes = core::mem::size_of::<PreparedAngularSpectrumApertureInterferometer>();
        let (mut max_f, mut max_real_f) = (0usize, 0usize);
        let (mut real, mut proper) = (false, false);
        for (bank, is_proper) in self.banks() {
            let f = bank.frequency_count();
            max_f = max_f.max(f);
            if is_proper {
                proper |= f > 0;
            } else {
                real |= f > 0;
                max_real_f = max_real_f.max(f);
            }
            mode_work = bank
                .line_count()
                .checked_mul(n)
                .and_then(|x| x.checked_mul(arms))
                .and_then(|x| mode_work.checked_add(x))
                .ok_or_else(overflow)?;
            let count = f
                .checked_mul(n)
                .and_then(|x| x.checked_mul(arms))
                .ok_or_else(overflow)?;
            entries = entries.checked_add(count).ok_or_else(overflow)?;
            let passes = nx
                .checked_add(ny)
                .and_then(|x| x.checked_add(f))
                .ok_or_else(overflow)?;
            products = f
                .checked_mul(n)
                .and_then(|x| x.checked_mul(passes))
                .and_then(|x| products.checked_add(x))
                .ok_or_else(overflow)?;
            let kernel_size = if is_proper {
                core::mem::size_of::<ProperKernel>()
            } else {
                core::mem::size_of::<RealKernel>()
            };
            for extra in [
                core::mem::size_of::<Bank>(),
                count
                    .checked_mul(core::mem::size_of::<FieldParts>())
                    .ok_or_else(overflow)?,
                f.checked_mul(core::mem::size_of::<f64>())
                    .ok_or_else(overflow)?,
                f.checked_mul(f)
                    .and_then(|x| x.checked_mul(kernel_size))
                    .ok_or_else(overflow)?,
            ] {
                bytes = bytes.checked_add(extra).ok_or_else(overflow)?;
            }
        }
        let spatial_pairs = nx
            .checked_mul(nx)
            .and_then(|x| ny.checked_mul(ny).and_then(|y| x.checked_add(y)))
            .ok_or_else(overflow)?;
        let scratch = max_f
            .checked_mul(3)
            .and_then(|x| x.checked_add(max_real_f))
            .and_then(|x| x.checked_mul(n))
            .and_then(|x| x.checked_mul(core::mem::size_of::<FieldParts>()))
            .ok_or_else(overflow)?;
        for extra in [
            spatial_pairs
                .checked_mul(
                    usize::from(real) * core::mem::size_of::<RealKernel>()
                        + usize::from(proper) * core::mem::size_of::<ProperKernel>(),
                )
                .ok_or_else(overflow)?,
            scratch,
            arms.checked_mul(core::mem::size_of::<usize>())
                .ok_or_else(overflow)?,
            if mode_work == 0 {
                0
            } else {
                n.checked_mul(core::mem::size_of::<FieldParts>() + 6 * core::mem::size_of::<C64>())
                    .ok_or_else(overflow)?
            },
        ] {
            bytes = bytes.checked_add(extra).ok_or_else(overflow)?;
        }
        Ok(AngularSpectrumApertureInterferenceLimits {
            maximum_mode_evaluations: mode_work,
            maximum_source_field_entries: entries,
            maximum_gram_products: products,
            maximum_workspace_bytes: bytes,
        })
    }

    fn prepare(
        &self,
        side: AngularSpectrumIlluminationSide,
        arms: &[AngularSpectrumExposureApertureArm],
        aperture: AngularSpectrumAperture,
        window: HarmonicWindow,
        limits: AngularSpectrumApertureInterferenceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedAngularSpectrumApertureInterferometer, AngularSpectrumError> {
        poll(cancelled)?;
        let needed = self.requirements(arms.len())?;
        if needed.maximum_mode_evaluations > limits.maximum_mode_evaluations
            || needed.maximum_source_field_entries > limits.maximum_source_field_entries
            || needed.maximum_gram_products > limits.maximum_gram_products
            || needed.maximum_workspace_bytes > limits.maximum_workspace_bytes
        {
            return Err(bad(
                "aperture interferometer exceeds its work or storage budget",
            ));
        }
        let mut projection_work = limits.maximum_gram_products - needed.maximum_gram_products;
        let mut projection_bytes = limits.maximum_workspace_bytes - needed.maximum_workspace_bytes;
        let mut projected_modes = 0;
        for arm in arms {
            poll(cancelled)?;
            if !arm.distance_from_face_m.is_finite()
                || arm.distance_from_face_m < 0.0
                || arm.translation_m.iter().any(|v| !v.is_finite())
            {
                return Err(bad(
                    "aperture arms need finite translations and nonnegative depth",
                ));
            }
        }
        let n = self.anchor.modes.len();
        let nx = self.anchor.grid.width;
        let mut order = reserved(arms.len())?;
        order.extend(0..arms.len());
        order.sort_unstable_by(|&a, &b| {
            arms[a]
                .distance_from_face_m
                .total_cmp(&arms[b].distance_from_face_m)
                .then_with(|| a.cmp(&b))
        });
        let mut banks = reserved(self.factors.len() + usize::from(self.mean.is_some()))?;
        let (mut real, mut proper) = (false, false);
        let (mut max_f, mut max_real_f) = (0usize, 0usize);
        for (source, is_proper) in self.banks() {
            let f = source.frequency_count();
            max_f = max_f.max(f);
            if is_proper {
                proper |= f > 0;
            } else {
                real |= f > 0;
                max_real_f = max_real_f.max(f);
            }
            let mut frequencies = reserved(f)?;
            frequencies.extend_from_slice(&source.frequencies);
            let temporal = if is_proper {
                Temporal::Proper(pair_kernels(f, window, |i| -frequencies[i], cancelled)?)
            } else {
                Temporal::Real(pair_kernels(f, window, |i| -frequencies[i], cancelled)?)
            };
            let positive = if let Some(plan) = temporal::Plan::new(
                self.anchor.grid(),
                window,
                &frequencies,
                is_proper,
                projected_modes,
            )? {
                let extra_work = plan
                    .extra_work
                    .checked_add(delay_projection_work(f, arms.len(), is_proper)?)
                    .ok_or_else(|| bad("aperture interferometer delay work overflow"))?;
                if extra_work <= projection_work && plan.extra_bytes <= projection_bytes {
                    let positive =
                        temporal::Positive::new(window, &frequencies, is_proper, cancelled)?;
                    projection_work -= extra_work;
                    projection_bytes -= plan.extra_bytes;
                    projected_modes = projected_modes.max(plan.modes);
                    Some(positive)
                } else {
                    None
                }
            } else {
                None
            };
            let mut fields = reserved(f * arms.len() * n)?;
            for i in 0..f * arms.len() * n {
                if i % 256 == 0 {
                    poll(cancelled)?;
                }
                fields.push([[Dd::ZERO; 2]; 6]);
            }
            for (frequency, &omega) in frequencies.iter().enumerate() {
                let mut begin = 0;
                while begin < order.len() {
                    let depth = arms[order[begin]].distance_from_face_m;
                    let parts = source.frequency_modal_parts(
                        frequency,
                        ExposureSide::TwoSided(side),
                        depth,
                        cancelled,
                    )?;
                    let mut end = begin;
                    while end < order.len() && arms[order[end]].distance_from_face_m == depth {
                        let arm = order[end];
                        let offset = (frequency * arms.len() + arm) * n;
                        let target = &mut fields[offset..offset + n];
                        for (destination, input) in target.chunks_mut(128).zip(parts.chunks(128)) {
                            poll(cancelled)?;
                            destination.copy_from_slice(input);
                        }
                        let shift = arms[arm].translation_m;
                        if is_proper {
                            translate(
                                target,
                                shift,
                                |i| transverse_wavenumber_dd(self.anchor, 0, i % nx),
                                |i| transverse_wavenumber_dd(self.anchor, 1, i / nx),
                                cancelled,
                            )?;
                        } else {
                            center_translated_fields_dd(
                                target,
                                omega,
                                |i| transverse_wavenumber_dd(self.anchor, 0, i % nx),
                                |i| transverse_wavenumber_dd(self.anchor, 1, i / nx),
                                aperture,
                                window,
                                shift,
                                cancelled,
                            )?;
                        }
                        end += 1;
                    }
                    begin = end;
                }
            }
            banks.push(Bank {
                frequencies,
                fields,
                temporal,
                positive,
            });
        }
        let spatial_real = if real {
            Some(Spatial::prepare(self.anchor, aperture, cancelled)?)
        } else {
            None
        };
        let spatial_proper = if proper {
            Some(Spatial::prepare(self.anchor, aperture, cancelled)?)
        } else {
            None
        };
        let gram = GramWorkspace::new(n, max_f, max_real_f, cancelled)?;
        let mut projected = reserved(n * projected_modes)?;
        for i in 0..n * projected_modes {
            if i % 128 == 0 {
                poll(cancelled)?;
            }
            projected.push([[Dd::ZERO; 2]; 6]);
        }
        poll(cancelled)?;
        Ok(PreparedAngularSpectrumApertureInterferometer {
            arm_count: arms.len(),
            modes: n,
            side,
            aperture,
            window,
            banks,
            spatial_real,
            spatial_proper,
            gram,
            projected,
        })
    }

    fn window_requirements(
        &self,
        arms: usize,
        window: HarmonicWindow,
    ) -> Result<AngularSpectrumApertureInterferenceLimits, AngularSpectrumError> {
        let mut limits = self.requirements(arms)?;
        let mut modes = 0;
        for (bank, proper) in self.banks() {
            if let Some(plan) =
                temporal::Plan::new(self.anchor.grid(), window, &bank.frequencies, proper, modes)?
            {
                limits.maximum_gram_products = limits
                    .maximum_gram_products
                    .checked_add(plan.extra_work)
                    .and_then(|work| {
                        work.checked_add(
                            delay_projection_work(bank.frequency_count(), arms, proper).ok()?,
                        )
                    })
                    .ok_or_else(|| bad("aperture projection work overflow"))?;
                limits.maximum_workspace_bytes = limits
                    .maximum_workspace_bytes
                    .checked_add(plan.extra_bytes)
                    .ok_or_else(|| bad("aperture projection workspace overflow"))?;
                modes = modes.max(plan.modes);
            }
        }
        Ok(limits)
    }
}

// A shared carrier rotation per arm and one complex phase composition per
// original frequency. Envelope rotations replace the original delay phases.
fn delay_projection_work(
    frequencies: usize,
    arms: usize,
    proper: bool,
) -> Result<usize, AngularSpectrumError> {
    if proper {
        frequencies
            .checked_add(1)
            .and_then(|count| count.checked_mul(arms))
            .ok_or_else(|| bad("aperture interferometer delay work overflow"))
    } else {
        Ok(0)
    }
}

fn reserved<T>(count: usize) -> Result<Vec<T>, AngularSpectrumError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| bad("aperture interferometer allocation failed"))?;
    Ok(values)
}
fn translate(
    fields: &mut [FieldParts],
    shift: [f64; 2],
    kx: impl Fn(usize) -> Dd,
    ky: impl Fn(usize) -> Dd,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), AngularSpectrumError> {
    for (i, components) in fields.iter_mut().enumerate() {
        if i % 32 == 0 {
            poll(cancelled)?;
        }
        let phase = kx(i) * Dd::from_f64(shift[0]) + ky(i) * Dd::from_f64(shift[1]);
        if !finite_dd(phase) {
            return Err(bad("aperture arm translation phase is unrepresentable"));
        }
        let (sine, cosine) = fs_math::det::sin_cos_dd(phase);
        for value in components {
            let [real, imaginary] = *value;
            *value = [
                real * cosine + imaginary * sine,
                imaginary * cosine - real * sine,
            ];
            if value.iter().any(|&part| !finite_dd(part)) {
                return Err(bad("translated aperture field is unrepresentable"));
            }
        }
    }
    Ok(())
}

impl AngularSpectrumExposure<'_, '_, '_> {
    /// Conservative bounds, before any reconstruction or allocation.
    pub fn aperture_interferometer_requirements(
        &self,
        arms: usize,
    ) -> Result<AngularSpectrumApertureInterferenceLimits, AngularSpectrumError> {
        Sources {
            anchor: self.anchor,
            mean: Some(self),
            factors: &[],
            statistics: None,
        }
        .requirements(arms)
    }

    /// Bounds for the actual shutter, reserving field-before-power projection
    /// for every supported narrow bank. No reconstruction or allocation occurs.
    pub fn aperture_interferometer_window_requirements(
        &self,
        arms: usize,
        window: HarmonicWindow,
    ) -> Result<AngularSpectrumApertureInterferenceLimits, AngularSpectrumError> {
        Sources {
            anchor: self.anchor,
            mean: Some(self),
            factors: &[],
            statistics: None,
        }
        .window_requirements(arms, window)
    }

    /// Prepare a shared-side finite detector from copies of this coherent bank.
    pub fn prepare_aperture_interferometer(
        &self,
        side: AngularSpectrumIlluminationSide,
        arms: &[AngularSpectrumExposureApertureArm],
        aperture: AngularSpectrumAperture,
        window: HarmonicWindow,
        limits: AngularSpectrumApertureInterferenceLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedAngularSpectrumApertureInterferometer, AngularSpectrumError> {
        Sources {
            anchor: self.anchor,
            mean: Some(self),
            factors: &[],
            statistics: None,
        }
        .prepare(side, arms, aperture, window, limits, &mut cancelled)
    }
}
impl AngularSpectrumExposureEnsemble<'_, '_, '_, '_> {
    /// Conservative bounds over the mean and all independent real/proper factors.
    pub fn aperture_interferometer_requirements(
        &self,
        arms: usize,
    ) -> Result<AngularSpectrumApertureInterferenceLimits, AngularSpectrumError> {
        Sources {
            anchor: self.anchor,
            mean: self.mean,
            factors: self.factors,
            statistics: self.statistics,
        }
        .requirements(arms)
    }

    /// Shutter-specific bounds for all independent original mean/factor banks.
    /// Positive temporal factors never create coherence between these banks.
    pub fn aperture_interferometer_window_requirements(
        &self,
        arms: usize,
        window: HarmonicWindow,
    ) -> Result<AngularSpectrumApertureInterferenceLimits, AngularSpectrumError> {
        Sources {
            anchor: self.anchor,
            mean: self.mean,
            factors: self.factors,
            statistics: self.statistics,
        }
        .window_requirements(arms, window)
    }

    /// Prepare independent source banks and coherently combine arms within each.
    pub fn prepare_aperture_interferometer(
        &self,
        side: AngularSpectrumIlluminationSide,
        arms: &[AngularSpectrumExposureApertureArm],
        aperture: AngularSpectrumAperture,
        window: HarmonicWindow,
        limits: AngularSpectrumApertureInterferenceLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedAngularSpectrumApertureInterferometer, AngularSpectrumError> {
        Sources {
            anchor: self.anchor,
            mean: self.mean,
            factors: self.factors,
            statistics: self.statistics,
        }
        .prepare(side, arms, aperture, window, limits, &mut cancelled)
    }
}
