//! Joint finite-area and finite-time observation of the original modal fields.
//!
//! Real fields use exp(+i omega t - i k.x). Separability of the rectangular
//! area and the authored shutter permits exact products of the shared harmonic
//! integrals. Deterministic means and real factors retain both carrier sum and
//! difference terms; proper-complex factors retain only differences. No spatial
//! sampling, optical-cycle average, detector back-action or source statistics
//! beyond the original exposure bank are introduced.

use super::*;

type FieldParts = [[Dd; 2]; 6];

/// Additional bounds for a joint aperture/shutter observation.
///
/// Pair work is N² times the sum of F² across the mean and factor banks, where
/// N is the modal grid size and each F counts that bank's distinct frequencies.
/// Reconstruction retains the original banks' admitted line-work bounds. Peak
/// logical workspace conservatively includes the largest original bank query,
/// retained DD fields, spatial caches and compensated result. Banks stream one at a time; no
/// dense source covariance, modal-pair table or new spatial FFT is allocated.
/// Borrowed owners, previous outputs and allocator overhead are excluded.
#[derive(Clone, Copy, Debug)]
pub struct AngularSpectrumExposureApertureLimits {
    pub maximum_mode_frequency_pairs: usize,
    pub maximum_workspace_bytes: usize,
}

/// Uniform area average and normalized shutter average of physical moments.
///
/// `mean` contains electric/magnetic mean squares and signed vector Poynting
/// flux in the units of [`AngularSpectrumExposurePoint`]. `normal_power_w` is
/// the mean Cartesian z-directed power through the rectangle, in watts, not
/// energy or intensity. The shutter may be uniform or explicitly shaped.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AngularSpectrumExposureApertureReading {
    pub mean: AngularSpectrumExposurePoint,
    pub normal_power_w: f64,
}

impl AngularSpectrumExposure<'_, '_, '_> {
    /// Integrate the continuous original Fourier fields over a finite detector
    /// and the authored shutter. Its rectangle may cross periodic-cell edges;
    /// the original periodic continuation remains the represented source.
    pub fn measure_aperture(
        &self,
        side: AngularSpectrumSide,
        distance_from_face_m: f64,
        aperture: AngularSpectrumAperture,
        window: HarmonicWindow,
        limits: AngularSpectrumExposureApertureLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureApertureReading, AngularSpectrumError> {
        measure(
            self.anchor,
            || core::iter::once((self, false)),
            ExposureSide::OneSided(side),
            distance_from_face_m,
            aperture,
            window,
            limits,
            &mut cancelled,
        )
    }

    /// Observe incoming, outgoing or total fields at either original face.
    /// The native two-sided modal reconstruction retains correlations between
    /// faces, polarizations, spatial modes and frequencies on one source clock.
    pub fn measure_illumination_aperture(
        &self,
        side: AngularSpectrumIlluminationSide,
        distance_from_face_m: f64,
        aperture: AngularSpectrumAperture,
        window: HarmonicWindow,
        limits: AngularSpectrumExposureApertureLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureApertureReading, AngularSpectrumError> {
        measure(
            self.anchor,
            || core::iter::once((self, false)),
            ExposureSide::TwoSided(side),
            distance_from_face_m,
            aperture,
            window,
            limits,
            &mut cancelled,
        )
    }
}

impl AngularSpectrumExposureEnsemble<'_, '_, '_, '_> {
    /// Expected finite-area, finite-shutter moments of the deterministic mean
    /// and authored real/proper factors. Their moments add without inventing
    /// cross-factor phase coherence or discarding within-factor interference.
    pub fn measure_aperture(
        &self,
        side: AngularSpectrumSide,
        distance_from_face_m: f64,
        aperture: AngularSpectrumAperture,
        window: HarmonicWindow,
        limits: AngularSpectrumExposureApertureLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureApertureReading, AngularSpectrumError> {
        measure(
            self.anchor,
            || self.banks(),
            ExposureSide::OneSided(side),
            distance_from_face_m,
            aperture,
            window,
            limits,
            &mut cancelled,
        )
    }

    /// Joint spatial/temporal observation of the original two-face ensemble.
    /// Right-face incoming radiation and local backflow retain negative z flux.
    pub fn measure_illumination_aperture(
        &self,
        side: AngularSpectrumIlluminationSide,
        distance_from_face_m: f64,
        aperture: AngularSpectrumAperture,
        window: HarmonicWindow,
        limits: AngularSpectrumExposureApertureLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureApertureReading, AngularSpectrumError> {
        measure(
            self.anchor,
            || self.banks(),
            ExposureSide::TwoSided(side),
            distance_from_face_m,
            aperture,
            window,
            limits,
            &mut cancelled,
        )
    }
}

#[allow(clippy::too_many_arguments)]
fn measure<'e, 'b: 'e, 's: 'b, 'm: 's, I>(
    anchor: &PlanarAngularSpectrum<'s, 'm>,
    banks: impl Fn() -> I,
    side: ExposureSide,
    distance: f64,
    aperture: AngularSpectrumAperture,
    window: HarmonicWindow,
    limits: AngularSpectrumExposureApertureLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<AngularSpectrumExposureApertureReading, AngularSpectrumError>
where
    I: Iterator<Item = (&'e AngularSpectrumExposure<'b, 's, 'm>, bool)>,
{
    poll(cancelled)?;
    if !distance.is_finite() || distance < 0.0 {
        return Err(bad(
            "aperture field distance must be finite and nonnegative",
        ));
    }
    let mut frequency_pairs = 0_usize;
    let mut bank_workspace = anchor.workspace_bytes();
    let (mut real, mut proper) = (false, false);
    for (bank, is_proper) in banks() {
        poll(cancelled)?;
        let f = bank.frequency_count();
        frequency_pairs = frequency_pairs
            .checked_add(
                f.checked_mul(f)
                    .ok_or_else(|| bad("aperture work overflow"))?,
            )
            .ok_or_else(|| bad("aperture work overflow"))?;
        // Retain grouped modal fields in DD through their joint center rotation.
        // The original query reserves their C64 counterpart, grouped DD scratch
        // and complete original reconstruction, which remain included below.
        let extra_fields = anchor
            .modes
            .len()
            .checked_mul(f)
            .and_then(|v| {
                v.checked_mul(core::mem::size_of::<FieldParts>() - core::mem::size_of::<[C64; 6]>())
            })
            .ok_or_else(|| bad("aperture DD field workspace overflow"))?;
        bank_workspace = bank_workspace.max(
            bank.workspace_bytes()
                .checked_add(extra_fields)
                .ok_or_else(|| bad("aperture DD field workspace overflow"))?,
        );
        real |= f != 0 && !is_proper;
        proper |= f != 0 && is_proper;
    }
    admit(
        anchor.grid(),
        frequency_pairs,
        bank_workspace,
        real,
        proper,
        limits,
    )?;
    // Both caches use the ACTUAL admitted modal wavevectors, not differences
    // derived from a rounded fundamental that could shift a distant fringe.
    let spatial_real = real
        .then(|| Spatial::<RealKernel>::prepare(anchor, aperture, cancelled))
        .transpose()?;
    let spatial_proper = proper
        .then(|| Spatial::<ProperKernel>::prepare(anchor, aperture, cancelled))
        .transpose()?;
    let mut sum = [Dd::ZERO; 9];
    for (bank, is_proper) in banks() {
        poll(cancelled)?;
        if bank.frequency_count() == 0 {
            continue;
        }
        if is_proper {
            accumulate_bank(
                bank,
                spatial_proper
                    .as_ref()
                    .ok_or_else(|| bad("missing proper aperture cache"))?,
                side,
                distance,
                window,
                &mut sum,
                cancelled,
            )?;
        } else {
            accumulate_bank(
                bank,
                spatial_real
                    .as_ref()
                    .ok_or_else(|| bad("missing real aperture cache"))?,
                side,
                distance,
                window,
                &mut sum,
                cancelled,
            )?;
        }
    }
    let values = sum.map(Dd::to_f64);
    let normal_power_w = (sum[8] * Dd::from_f64(aperture.area_m2())).to_f64();
    if values.iter().any(|v| !v.is_finite())
        || values[..6].iter().any(|v| *v < 0.0)
        || !normal_power_w.is_finite()
    {
        return Err(bad("aperture exposure moment is unrepresentable"));
    }
    poll(cancelled)?;
    Ok(AngularSpectrumExposureApertureReading {
        mean: AngularSpectrumExposurePoint {
            electric_mean_squared_v2_m2: [values[0], values[1], values[2]],
            magnetic_mean_squared_a2_m2: [values[3], values[4], values[5]],
            poynting_w_m2: [values[6], values[7], values[8]],
        },
        normal_power_w,
    })
}

fn admit(
    grid: AngularSpectrumGrid,
    frequency_pairs: usize,
    bank_workspace: usize,
    real: bool,
    proper: bool,
    limits: AngularSpectrumExposureApertureLimits,
) -> Result<(), AngularSpectrumError> {
    let n = grid
        .width
        .checked_mul(grid.height)
        .ok_or_else(|| bad("aperture grid overflow"))?;
    let pairs = n
        .checked_mul(n)
        .and_then(|n2| n2.checked_mul(frequency_pairs))
        .ok_or_else(|| bad("aperture work overflow"))?;
    let axis_pairs = grid
        .width
        .checked_mul(grid.width)
        .and_then(|x| {
            grid.height
                .checked_mul(grid.height)
                .and_then(|y| x.checked_add(y))
        })
        .ok_or_else(|| bad("aperture spatial cache overflow"))?;
    let kernel_size = usize::from(real) * core::mem::size_of::<RealKernel>()
        + usize::from(proper) * core::mem::size_of::<ProperKernel>();
    // The original exposure query reserves frequency fields, DD modal sums,
    // original line reconstruction and F² HarmonicProducts. Both temporal
    // kernel representations below fit in that reservation.
    let bytes = axis_pairs
        .checked_mul(kernel_size)
        .and_then(|v| v.checked_add(bank_workspace))
        .and_then(|v| {
            v.checked_add(
                core::mem::size_of::<Spatial<RealKernel>>()
                    + core::mem::size_of::<Spatial<ProperKernel>>()
                    + core::mem::size_of::<[Dd; 9]>()
                    + core::mem::size_of::<AngularSpectrumExposureApertureReading>(),
            )
        })
        .ok_or_else(|| bad("aperture workspace overflow"))?;
    if pairs > limits.maximum_mode_frequency_pairs || bytes > limits.maximum_workspace_bytes {
        return Err(bad(
            "aperture exposure exceeds its pair work or workspace budget",
        ));
    }
    Ok(())
}

trait Kernel: Copy {
    const CENTER_FIELDS: bool;
    fn prepare(window: HarmonicWindow, a: f64, b: f64) -> Result<Self, AngularSpectrumError>;
    fn compose(self, other: Self) -> Self;
    fn product(self, a: [Dd; 2], b: [Dd; 2]) -> Dd;
}

/// Centered <cos(a)cos(b)> and <sin(a)sin(b)>; centered mixed means vanish.
/// Modal amplitudes receive their joint center rotation once, before products.
/// Thus a translated dark fringe does not arise from cancellation between four
/// large uncentered Gram entries and lose the small positive window variance.
#[derive(Clone, Copy)]
struct RealKernel([Dd; 2]);
impl Kernel for RealKernel {
    const CENTER_FIELDS: bool = true;
    fn prepare(window: HarmonicWindow, a: f64, b: f64) -> Result<Self, AngularSpectrumError> {
        // Centered weights are translation-independent. Reuse the original
        // small-width/sinc implementation, including its subnormal handling.
        let kernel = HarmonicWindow::with_profile(0.0, window.duration_s(), window.profile())
            .map_err(|_| bad("aperture centered interval is unrepresentable"))?
            .product_kernel(a, b)
            .map_err(|_| bad("aperture real phase kernel is unrepresentable"))?;
        Ok(Self(kernel.centered_gram_dd()))
    }

    fn compose(self, other: Self) -> Self {
        let [c, s] = self.0;
        let [d, t] = other.0;
        Self([c * d + s * t, c * t + s * d])
    }

    fn product(self, a: [Dd; 2], b: [Dd; 2]) -> Dd {
        a[0] * b[0] * self.0[0] + a[1] * b[1] * self.0[1]
    }
}

/// Complex mean exp(-i (a-b) coordinate), with the physical factor 1/2 applied
/// once at contraction. Difference kernels preserve slow resolved beats even
/// when individual carrier phases would round differently at distant clocks.
#[derive(Clone, Copy)]
struct ProperKernel([Dd; 2]);
impl Kernel for ProperKernel {
    const CENTER_FIELDS: bool = false;
    fn prepare(window: HarmonicWindow, a: f64, b: f64) -> Result<Self, AngularSpectrumError> {
        let kernel = window
            .proper_product_kernel(a, b)
            .map_err(|_| bad("aperture proper phase kernel is unrepresentable"))?;
        let real = kernel
            .mean_real_product_dd(C64::ONE, C64::ONE)
            .map_err(|_| bad("aperture proper phase moment is unrepresentable"))?;
        let negative_imaginary = kernel
            .mean_real_product_dd(C64::new(0.0, 1.0), C64::ONE)
            .map_err(|_| bad("aperture proper phase moment is unrepresentable"))?;
        Ok(Self([
            real * Dd::from_f64(2.0),
            negative_imaginary * Dd::from_f64(-2.0),
        ]))
    }
    fn compose(self, other: Self) -> Self {
        let [a, b] = self.0;
        let [c, d] = other.0;
        Self([a * c - b * d, a * d + b * c])
    }
    fn product(self, [ar, ai]: [Dd; 2], [br, bi]: [Dd; 2]) -> Dd {
        ((ar * br + ai * bi) * self.0[0] - (ai * br - ar * bi) * self.0[1]) * Dd::from_f64(0.5)
    }
}

struct Spatial<K> {
    grid: AngularSpectrumGrid,
    aperture: AngularSpectrumAperture,
    x: Vec<K>,
    y: Vec<K>,
}
impl<K: Kernel> Spatial<K> {
    fn prepare(
        anchor: &PlanarAngularSpectrum<'_, '_>,
        aperture: AngularSpectrumAperture,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AngularSpectrumError> {
        let grid = anchor.grid();
        Self::from_axes(
            grid,
            aperture,
            |i| anchor.modes[i].kx,
            |i| anchor.modes[i * grid.width].ky,
            cancelled,
        )
    }

    fn from_axes(
        grid: AngularSpectrumGrid,
        aperture: AngularSpectrumAperture,
        kx: impl Fn(usize) -> f64,
        ky: impl Fn(usize) -> f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AngularSpectrumError> {
        let origin = aperture.origin_m();
        let extent = aperture.extent_m();
        let x_window = HarmonicWindow::new(origin[0], extent[0])
            .map_err(|_| bad("aperture x interval is unrepresentable"))?;
        let y_window = HarmonicWindow::new(origin[1], extent[1])
            .map_err(|_| bad("aperture y interval is unrepresentable"))?;
        Ok(Self {
            grid,
            aperture,
            x: pair_kernels(grid.width, x_window, kx, cancelled)?,
            y: pair_kernels(grid.height, y_window, ky, cancelled)?,
        })
    }

    fn pair(&self, i: usize, j: usize) -> K {
        self.x[(i % self.grid.width) * self.grid.width + j % self.grid.width]
            .compose(self.y[(i / self.grid.width) * self.grid.height + j / self.grid.width])
    }
}

fn pair_kernels<K: Kernel>(
    count: usize,
    window: HarmonicWindow,
    frequency: impl Fn(usize) -> f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<K>, AngularSpectrumError> {
    let pairs = count
        .checked_mul(count)
        .ok_or_else(|| bad("aperture kernel count overflow"))?;
    let mut kernels = Vec::new();
    kernels
        .try_reserve_exact(pairs)
        .map_err(|_| bad("aperture kernel allocation failed"))?;
    for i in 0..count {
        for j in 0..count {
            if j % 32 == 0 {
                poll(cancelled)?;
            }
            kernels.push(K::prepare(window, frequency(i), frequency(j))?);
        }
    }
    Ok(kernels)
}

#[allow(clippy::too_many_arguments)]
fn accumulate_bank<K: Kernel>(
    bank: &AngularSpectrumExposure<'_, '_, '_>,
    spatial: &Spatial<K>,
    side: ExposureSide,
    distance: f64,
    window: HarmonicWindow,
    sum: &mut [Dd; 9],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), AngularSpectrumError> {
    let f = bank.frequency_count();
    // Shared harmonic kernels use negative-time phasors; EM uses positive time.
    let temporal = pair_kernels(f, window, |i| -bank.frequencies[i], cancelled)?;
    let mut fields = Vec::new();
    fields
        .try_reserve_exact(f)
        .map_err(|_| bad("aperture frequency field allocation failed"))?;
    for group in 0..f {
        let mut field = bank.frequency_modal_parts(group, side, distance, cancelled)?;
        if K::CENTER_FIELDS {
            center_fields(
                &mut field,
                bank.frequencies[group],
                |i| bank.anchor.modes[i].kx,
                |i| bank.anchor.modes[i].ky,
                spatial.aperture,
                window,
                cancelled,
            )?;
        }
        fields.push(field);
    }
    integrate(spatial, &temporal, &fields, sum, cancelled)
}

fn integrate<K: Kernel>(
    spatial: &Spatial<K>,
    temporal: &[K],
    fields: &[Vec<FieldParts>],
    sum: &mut [Dd; 9],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), AngularSpectrumError> {
    let n = spatial.grid.width * spatial.grid.height;
    let f = fields.len();
    for i in 0..n {
        poll(cancelled)?;
        for j in 0..n {
            let space = spatial.pair(i, j);
            for (a, left) in fields.iter().enumerate() {
                for (b, right) in fields.iter().enumerate() {
                    if b % 16 == 0 {
                        poll(cancelled)?;
                    }
                    let kernel = space.compose(temporal[a * f + b]);
                    let (a, b) = (left[i], right[j]);
                    for c in 0..3 {
                        sum[c] = sum[c] + kernel.product(a[c], b[c]);
                        sum[c + 3] = sum[c + 3] + kernel.product(a[c + 3], b[c + 3]);
                        let (u, v) = ((c + 1) % 3, (c + 2) % 3);
                        sum[c + 6] = sum[c + 6] + kernel.product(a[u], b[v + 3])
                            - kernel.product(a[v], b[u + 3]);
                    }
                }
            }
        }
    }
    Ok(())
}

// Apply exp(+i omega*t_center - i k.x_center) once per grouped mode, in DD.
// Add all start phases BEFORE trigonometry so a translated spatiotemporal
// null cancels as a phase. Keep half-width phases separately: a very short
// detector interval must not disappear into a large source coordinate/clock.
#[allow(clippy::too_many_arguments)]
fn center_fields(
    field: &mut [FieldParts],
    omega: f64,
    kx: impl Fn(usize) -> f64,
    ky: impl Fn(usize) -> f64,
    aperture: AngularSpectrumAperture,
    window: HarmonicWindow,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), AngularSpectrumError> {
    let [x, y] = aperture.origin_m();
    let [width, height] = aperture.extent_m();
    for (i, components) in field.iter_mut().enumerate() {
        if i % 32 == 0 {
            poll(cancelled)?;
        }
        let (mut start, mut half_width) = (Dd::ZERO, Dd::ZERO);
        for (frequency, origin, extent) in [
            (kx(i), x, width),
            (ky(i), y, height),
            (-omega, window.start_s(), window.duration_s()),
        ] {
            let frequency_dd = Dd::from_f64(frequency);
            start = start + frequency_dd * Dd::from_f64(origin);
            // Halve the larger operand to retain representable products when
            // either the interval or frequency is subnormal.
            half_width = half_width
                + if frequency.abs() >= extent {
                    (frequency_dd * Dd::from_f64(0.5)) * Dd::from_f64(extent)
                } else {
                    frequency_dd * Dd::from_f64(0.5 * extent)
                };
        }
        if !finite_dd(start) || !finite_dd(half_width) {
            return Err(bad("joint aperture center phase is unrepresentable"));
        }
        let (s, c) = fs_math::det::sin_cos_dd(start);
        let (t, d) = fs_math::det::sin_cos_dd(half_width);
        let (cosine, sine) = (c * d - s * t, s * d + c * t);
        for value in components {
            let [real, imaginary] = *value;
            *value = [
                real * cosine + imaginary * sine,
                imaginary * cosine - real * sine,
            ];
            if value.iter().any(|&v| !finite_dd(v)) {
                return Err(bad("joint aperture centered field is unrepresentable"));
            }
        }
    }
    Ok(())
}

fn finite_dd(value: Dd) -> bool {
    value.hi.is_finite() && value.lo.is_finite()
}

#[cfg(test)]
mod tests;
