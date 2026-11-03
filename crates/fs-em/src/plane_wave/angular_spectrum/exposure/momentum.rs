//! Finite-window electromagnetic momentum from the original exterior fields.
//!
//! In one homogeneous, nondispersive lossless exterior the explicitly named
//! Minkowski density is g_M = D cross B = epsilon*mu*(E cross H).
//! In vacuum this is the electromagnetic momentum density S/c^2. A material
//! exterior does not acquire a mechanical momentum, striction or moving-body
//! law from this observable. The volume integral of its rate is the field
//! term needed alongside a closed Maxwell-stress surface integral; one plane
//! or its transverse quadrature is not a net force on a body.
//!
//! Reference for the vacuum conservation law: MIT 22.105, chapter 3,
//! https://ocw.mit.edu/courses/22-105-electromagnetic-interactions-fall-2005/225368a969b5f565d469525582fdb139_chap3.pdf.
//! Real fields use exp(+i omega t). Both the density and its TIME DERIVATIVE
//! are averaged with the original normalized shutter. Only a uniform shutter
//! makes the latter equal to an endpoint difference divided by duration.

use super::*;

/// Complete query storage, including original reconstruction/FFT workspace,
/// cached pair kernels, compensated accumulation and the returned plane.
/// Existing line/pair/FFT work admissions remain authoritative. A real pair
/// uses eighteen scalar temporal contractions; a proper pair uses at most
/// twelve. Each query contracts both N point samples and Q native spatial
/// quadratures, with Q <= 2*N, within the original O(N*F^2) pair-work bound.
/// No additional FFT or retained frequency bank is needed. Fixed modal
/// contraction scratch is included. Earlier outputs, borrowed owners and
/// allocator overhead are excluded.
#[derive(Clone, Copy, Debug)]
pub struct AngularSpectrumExposureMomentumLimits {
    pub maximum_workspace_bytes: usize,
}

/// Cartesian, signed window/ensemble means in the original physical frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AngularSpectrumExposureMomentumPoint {
    /// Mean D cross B [N s/m^3]. This equals S/c^2 in vacuum.
    pub minkowski_density_n_s_m3: [f64; 3],
    /// Mean d(D cross B)/dt [N/m^3], with the SAME normalized shutter.
    /// This is field-momentum storage rate, not a mechanical force density.
    pub minkowski_rate_n_m3: [f64; 3],
}

/// Original-grid momentum samples and exact continuous periodic-cell integrals.
/// Modal orthogonality integrates the represented fields before sampling
/// their quadratic products, retaining signed Nyquist modes without aliasing.
/// Integrating these plane integrals through depth yields N s and N. A caller
/// owns that volume quadrature and its convergence; no finite-aperture edge
/// model or interior material momentum is inferred.
#[derive(Debug)]
pub struct AngularSpectrumExposureMomentumField {
    grid: AngularSpectrumGrid,
    points: Vec<AngularSpectrumExposureMomentumPoint>,
    integrated: [Dd; 6],
}

impl AngularSpectrumExposureMomentumField {
    pub fn grid(&self) -> AngularSpectrumGrid {
        self.grid
    }

    pub fn points(&self) -> &[AngularSpectrumExposureMomentumPoint] {
        &self.points
    }

    /// Mean density integrated over the continuous periodic transverse cell.
    /// This can differ from summing the original point samples.
    pub fn momentum_per_depth_n_s_m(&self) -> [f64; 3] {
        std::array::from_fn(|c| self.integrated[c].to_f64())
    }

    /// Weighted momentum rate integrated over the continuous periodic cell.
    /// This can differ from summing the original point samples.
    pub fn momentum_rate_per_depth_n_m(&self) -> [f64; 3] {
        std::array::from_fn(|c| self.integrated[c + 3].to_f64())
    }
}

impl AngularSpectrumExposure<'_, '_, '_> {
    pub fn momentum_workspace_bytes(&self) -> Result<usize, AngularSpectrumError> {
        momentum_workspace(
            self.anchor.modes.len(),
            self.workspace_bytes(),
            self.frequency_count(),
        )
    }

    /// Observe original incoming/reflected/total/transmitted exterior fields.
    /// Both exteriors must retain the same source-identical lossless medium,
    /// as for the original Maxwell-stress owner. Complete Cartesian fields
    /// interfere before products; no cycle average or phase independence is
    /// assumed. No new film solve or constitutive sampling per frequency is
    /// performed. Work is O(L*N + F*N*log N + F^2*N).
    pub fn momentum_field(
        &self,
        side: AngularSpectrumSide,
        distance_from_face_m: f64,
        window: HarmonicWindow,
        limits: AngularSpectrumExposureMomentumLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureMomentumField, AngularSpectrumError> {
        momentum_field(
            self.anchor,
            core::iter::once((self, false)),
            self.momentum_workspace_bytes()?,
            ExposureSide::OneSided(side),
            distance_from_face_m,
            window,
            limits,
            &mut cancelled,
        )
    }

    /// Observe the original incoming, outgoing or TOTAL field at either face,
    /// retaining counterpropagating, polarization and spatial interference.
    pub fn illumination_momentum_field(
        &self,
        side: AngularSpectrumIlluminationSide,
        distance_from_face_m: f64,
        window: HarmonicWindow,
        limits: AngularSpectrumExposureMomentumLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureMomentumField, AngularSpectrumError> {
        momentum_field(
            self.anchor,
            core::iter::once((self, false)),
            self.momentum_workspace_bytes()?,
            ExposureSide::TwoSided(side),
            distance_from_face_m,
            window,
            limits,
            &mut cancelled,
        )
    }
}

impl AngularSpectrumExposureEnsemble<'_, '_, '_, '_> {
    pub fn momentum_workspace_bytes(&self) -> Result<usize, AngularSpectrumError> {
        let frequencies = self.banks().map(|(bank, _)| bank.frequency_count()).max().unwrap_or(0);
        momentum_workspace(self.anchor.modes.len(), self.workspace_bytes(), frequencies)
    }

    /// Mean plus authored real/proper covariance factors. Only the declared
    /// zero pseudocovariance of a proper factor removes sum-frequency terms.
    /// One bank's original frequency fields are retained at a time.
    pub fn momentum_field(
        &self,
        side: AngularSpectrumSide,
        distance_from_face_m: f64,
        window: HarmonicWindow,
        limits: AngularSpectrumExposureMomentumLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureMomentumField, AngularSpectrumError> {
        momentum_field(
            self.anchor,
            self.banks(),
            self.momentum_workspace_bytes()?,
            ExposureSide::OneSided(side),
            distance_from_face_m,
            window,
            limits,
            &mut cancelled,
        )
    }

    pub fn illumination_momentum_field(
        &self,
        side: AngularSpectrumIlluminationSide,
        distance_from_face_m: f64,
        window: HarmonicWindow,
        limits: AngularSpectrumExposureMomentumLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureMomentumField, AngularSpectrumError> {
        momentum_field(
            self.anchor,
            self.banks(),
            self.momentum_workspace_bytes()?,
            ExposureSide::TwoSided(side),
            distance_from_face_m,
            window,
            limits,
            &mut cancelled,
        )
    }
}

#[derive(Clone, Copy)]
struct MomentumKernel<P> {
    product: P,
    mean_frequency: Dd,
    half_difference: Dd,
}

fn momentum_workspace(
    points: usize,
    original: usize,
    frequencies: usize,
) -> Result<usize, AngularSpectrumError> {
    points
        .checked_mul(
            core::mem::size_of::<[Dd; 6]>()
                + core::mem::size_of::<AngularSpectrumExposureMomentumPoint>(),
        )
        .and_then(|v| v.checked_add(original))
        .and_then(|v| {
            v.checked_add(frequencies.checked_mul(frequencies)?
                .checked_mul(2 * core::mem::size_of::<Dd>())?)
        })
        // Retained and transient cosine/sine DD vectors, plus the cell,
        // per-mode and accumulated-bank six-component integrals.
        .and_then(|v| v.checked_add(
            4 * core::mem::size_of::<[[[Dd; 2]; 2]; 6]>()
                + 3 * core::mem::size_of::<[Dd; 6]>(),
        ))
        .ok_or_else(|| bad("exposure momentum workspace overflow"))
}

fn momentum_field<'e, 'b, 's, 'm>(
    anchor: &PlanarAngularSpectrum<'s, 'm>,
    banks: impl IntoIterator<Item = (&'e AngularSpectrumExposure<'b, 's, 'm>, bool)>,
    workspace: usize,
    side: ExposureSide,
    distance: f64,
    window: HarmonicWindow,
    limits: AngularSpectrumExposureMomentumLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<AngularSpectrumExposureMomentumField, AngularSpectrumError>
where
    'b: 'e,
    's: 'b,
    'm: 's,
{
    poll(cancelled)?;
    if workspace > limits.maximum_workspace_bytes {
        return Err(bad("exposure momentum exceeds its workspace budget"));
    }
    if !distance.is_finite() || distance < 0.0 {
        return Err(bad("momentum field distance must be finite and nonnegative"));
    }
    anchor.source().admit_momentum_exterior()?;
    let medium = anchor.source().incident_material();
    let epsilon = linear_electric_permittivity(medium).map_err(ReflectionError::from)?;
    let mu = linear_magnetic_permeability(medium).map_err(ReflectionError::from)?;
    let n = anchor.modes.len();
    let mut sums = Vec::new();
    sums.try_reserve_exact(n)
        .map_err(|_| bad("exposure momentum accumulation allocation failed"))?;
    for i in 0..n {
        if i % 64 == 0 { poll(cancelled)?; }
        sums.push([Dd::ZERO; 6]);
    }
    let mut integrated = [Dd::ZERO; 6];
    for (bank, proper) in banks {
        poll(cancelled)?;
        let cell = if proper {
            accumulate_bank(
                bank, side, distance, true,
                |a, b| window.proper_product_kernel(a, b),
                proper_product, ProperHarmonicProduct::mean_real_product_parts_dd,
                &mut sums, cancelled,
            )?
        } else {
            accumulate_bank(
                bank, side, distance, false,
                |a, b| window.product_kernel(a, b),
                product, HarmonicProduct::mean_real_product_parts_dd,
                &mut sums, cancelled,
            )?
        };
        for c in 0..6 { integrated[c] = integrated[c] + cell[c]; }
    }
    let mut points = Vec::new();
    points.try_reserve_exact(n)
        .map_err(|_| bad("exposure momentum result allocation failed"))?;
    for (i, sum) in sums.into_iter().enumerate() {
        if i % 64 == 0 { poll(cancelled)?; }
        let mut scaled = [Dd::ZERO; 6];
        for c in 0..6 {
            scaled[c] = scale_momentum(sum[c], epsilon, mu)?;
        }
        points.push(AngularSpectrumExposureMomentumPoint {
            minkowski_density_n_s_m3: std::array::from_fn(|c| scaled[c].to_f64()),
            minkowski_rate_n_m3: std::array::from_fn(|c| scaled[c + 3].to_f64()),
        });
    }
    for value in &mut integrated {
        *value = scale_momentum(*value, epsilon, mu)?;
    }
    let grid = anchor.grid();
    poll(cancelled)?;
    Ok(AngularSpectrumExposureMomentumField { grid, points, integrated })
}

// Avoid creating an underflowed epsilon*mu before multiplying a resolved
// field moment. Mixed coefficients are combined first; otherwise the larger
// coefficient scales the moment first. No gain clipping or zero repair.
fn scale_momentum(value: Dd, epsilon: f64, mu: f64) -> Result<Dd, AngularSpectrumError> {
    let result = if (epsilon >= 1.0) != (mu >= 1.0) {
        value * (Dd::from_f64(epsilon) * Dd::from_f64(mu))
    } else {
        (value * Dd::from_f64(epsilon.max(mu))) * Dd::from_f64(epsilon.min(mu))
    };
    if !result.to_f64().is_finite() || (value != Dd::ZERO && result == Dd::ZERO) {
        return Err(bad("exterior momentum scaling is unrepresentable"));
    }
    Ok(result)
}

fn accumulate_bank<P: Copy>(
    bank: &AngularSpectrumExposure<'_, '_, '_>,
    side: ExposureSide,
    distance: f64,
    proper: bool,
    kernel: impl Fn(f64, f64) -> Result<P, fs_fft::harmonic_window::WindowError>,
    product: impl Fn(P, C64, C64) -> Result<Dd, MomentError>,
    product_parts: impl Fn(P, [Dd; 2], [Dd; 2])
        -> Result<Dd, fs_fft::harmonic_window::WindowError>,
    sums: &mut [[Dd; 6]],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[Dd; 6], AngularSpectrumError> {
    let f = bank.frequency_count();
    if f == 0 { return Ok([Dd::ZERO; 6]); }
    let mut kernels = Vec::new();
    kernels.try_reserve_exact(f * f)
        .map_err(|_| bad("exposure momentum kernel allocation failed"))?;
    let half = Dd::from_f64(0.5);
    for &a in &bank.frequencies {
        for &b in &bank.frequencies {
            poll(cancelled)?;
            kernels.push(MomentumKernel {
                product: kernel(a, b)
                    .map_err(|_| bad("exposure momentum temporal kernel is unrepresentable"))?,
                mean_frequency: Dd::from_f64(a) * half + Dd::from_f64(b) * half,
                half_difference: (Dd::from_f64(a) - Dd::from_f64(b)) * half,
            });
        }
    }
    let mut fields = Vec::new();
    fields.try_reserve_exact(f)
        .map_err(|_| bad("exposure momentum field allocation failed"))?;
    for group in 0..f {
        fields.push(bank.frequency_modal_field(group, side, distance, cancelled)?);
    }
    let cell = integrate_cell_momentum(
        bank.anchor.grid(), &fields, &kernels, proper, product_parts, cancelled,
    )?;
    // The same retained modal bank becomes the original sample bank in place.
    for (group, field) in fields.iter_mut().enumerate() {
        bank.synthesize_frequency_field(group, field, cancelled)?;
    }
    for (point, sum) in sums.iter_mut().enumerate() {
        poll(cancelled)?;
        let value = integrate_momentum(
            f, &kernels, |i| fields[i][point], proper, rotate_i, &product, cancelled,
        ).map_err(momentum_error)?;
        for c in 0..6 { sum[c] = sum[c] + value[c]; }
    }
    Ok(cell)
}

fn momentum_error(error: MomentError) -> AngularSpectrumError {
    match error {
        MomentError::Cancelled => AngularSpectrumError::Cancelled,
        MomentError::Product(_) | MomentError::Unrepresentable =>
            bad("exposure momentum product is unrepresentable"),
    }
}

// Reuse the original signed-mode pairing and DD spatial quadratures. A
// sampled product aliases a Nyquist traveling wave's carrier into a false
// spatially uniform momentum/storage term; its continuous cell integral
// must instead retain the half-weighted sine AND cosine components.
// The spatial product is contracted in DD before carrier/beat differentiation,
// so stationary spatial sums cancel before multiplying a large frequency.
fn integrate_cell_momentum<P: Copy>(
    grid: AngularSpectrumGrid,
    fields: &[Vec<[C64; 6]>],
    kernels: &[MomentumKernel<P>],
    proper: bool,
    product: impl Fn(P, [Dd; 2], [Dd; 2])
        -> Result<Dd, fs_fft::harmonic_window::WindowError>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[Dd; 6], AngularSpectrumError> {
    let mut sum = [Dd::ZERO; 6];
    for mode in 0..grid.width * grid.height {
        if mode % 32 == 0 { poll(cancelled)?; }
        let opposite = opposite_mode(grid, mode);
        if opposite.is_some_and(|other| other < mode) { continue; }
        let constant = opposite == Some(mode);
        let partner = opposite.filter(|_| !constant);
        let quadratures = if constant { 1 } else { 2 };
        let weight = Dd::from_f64(if constant { 1.0 } else { 0.5 });
        let value = integrate_momentum(
            fields.len(), kernels,
            |i| {
                let parts = spatial_quadratures(
                    fields[i][mode], partner.map(|other| fields[i][other]),
                );
                core::array::from_fn(|c| [parts[0][c], parts[1][c]])
            },
            proper,
            // spatial_quadratures already conjugates amplitudes into the
            // shared NEGATIVE-time convention. Its derivative rotates -i,
            // unlike the original +i C64 point amplitudes below. No DD
            // modal sum/difference is rounded back to C64 here.
            |value: [[Dd; 2]; 2]| value.map(|part| [part[1], -part[0]]),
            |kernel, a: [[Dd; 2]; 2], b: [[Dd; 2]; 2]| {
                let mut value = Dd::ZERO;
                for q in 0..quadratures {
                    value = value + product(kernel, a[q], b[q])
                        .map_err(MomentError::Product)?;
                }
                Ok(value * weight)
            },
            cancelled,
        ).map_err(momentum_error)?;
        for c in 0..6 { sum[c] = sum[c] + value[c]; }
    }
    let area = Dd::from_f64(grid.width as f64 * grid.spacing_m[0])
        * Dd::from_f64(grid.height as f64 * grid.spacing_m[1]);
    let integrated = sum.map(|v| v * area);
    if integrated.iter().any(|v| !v.hi.is_finite() || !v.lo.is_finite()) {
        return Err(bad("integrated exposure momentum is unrepresentable"));
    }
    poll(cancelled)?;
    Ok(integrated)
}

fn rotate_i(value: C64) -> C64 {
    C64::new(-value.im, value.re)
}

fn integrate_momentum<P: Copy, A: Copy>(
    frequencies: usize,
    kernels: &[MomentumKernel<P>],
    field_at: impl Fn(usize) -> [A; 6],
    proper: bool,
    rotate: impl Fn(A) -> A,
    product: impl Fn(P, A, A) -> Result<Dd, MomentError>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[Dd; 6], MomentError> {
    let mut sum = [Dd::ZERO; 6];
    for i in 0..frequencies {
        let a = field_at(i);
        for j in 0..frequencies {
            if j % 64 == 0 && cancelled() { return Err(MomentError::Cancelled); }
            let b = field_at(j);
            let kernel = kernels[i * frequencies + j];
            for c in 0..3 {
                let (u, v) = ((c + 1) % 3, (c + 2) % 3);
                for (e, h, sign) in [(u, v + 3, 1.0), (v, u + 3, -1.0)] {
                    let value = product(kernel.product, a[e], b[h])?;
                    let rate = if proper && kernel.half_difference == Dd::ZERO {
                        // A proper monochromatic second moment is stationary,
                        // including polarized/coherent spatial interference.
                        Dd::ZERO
                    } else {
                        let left = product(kernel.product, rotate(a[e]), b[h])?;
                        if proper {
                            left * (kernel.half_difference * Dd::from_f64(2.0))
                        } else {
                            let right = product(kernel.product, a[e], rotate(b[h]))?;
                            // Separate carrier and beat terms BEFORE multiplying
                            // a large carrier. No rounded subtraction occurs.
                            (left + right) * kernel.mean_frequency
                                + (left - right) * kernel.half_difference
                        }
                    };
                    let sign = Dd::from_f64(sign);
                    sum[c] = sum[c] + sign * value;
                    sum[c + 3] = sum[c + 3] + sign * rate;
                }
            }
        }
    }
    if sum.iter().any(|v| !v.to_f64().is_finite()) {
        return Err(MomentError::Unrepresentable);
    }
    Ok(sum)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_fft::harmonic_window::HarmonicWindowProfile;

    fn evaluate(
        fields: &[[C64; 6]],
        frequencies: &[f64],
        window: HarmonicWindow,
        proper: bool,
    ) -> [f64; 6] {
        let half = Dd::from_f64(0.5);
        if proper {
            let kernels = frequencies.iter().flat_map(|&a| {
                frequencies.iter().map(move |&b| MomentumKernel {
                    product: window.proper_product_kernel(a, b).unwrap(),
                    mean_frequency: Dd::from_f64(a) * half + Dd::from_f64(b) * half,
                    half_difference: (Dd::from_f64(a) - Dd::from_f64(b)) * half,
                })
            }).collect::<Vec<_>>();
            integrate_momentum(frequencies.len(), &kernels, |i| fields[i],
                true, rotate_i, proper_product, &mut || false).unwrap().map(Dd::to_f64)
        } else {
            let kernels = frequencies.iter().flat_map(|&a| {
                frequencies.iter().map(move |&b| MomentumKernel {
                    product: window.product_kernel(a, b).unwrap(),
                    mean_frequency: Dd::from_f64(a) * half + Dd::from_f64(b) * half,
                    half_difference: (Dd::from_f64(a) - Dd::from_f64(b)) * half,
                })
            }).collect::<Vec<_>>();
            integrate_momentum(frequencies.len(), &kernels, |i| fields[i],
                false, rotate_i, product, &mut || false).unwrap().map(Dd::to_f64)
        }
    }

    #[test]
    fn g1_mean_and_weighted_rate_match_real_fields_and_authored_phase_ensemble() {
        // Arbitrary Cartesian phasors isolate the vector/temporal contraction.
        // The source-backed integration tests separately exercise Maxwell fields.
        let fields = [
            [C64::new(0.8, 0.2), C64::new(-0.2, 0.5), C64::new(0.3, -0.4),
             C64::new(-0.5, 0.1), C64::new(0.6, 0.2), C64::new(-0.1, 0.3)],
            [C64::new(-0.4, 0.7), C64::new(0.2, 0.3), C64::new(-0.6, 0.1),
             C64::new(0.3, -0.2), C64::new(0.4, -0.1), C64::new(0.5, 0.6)],
        ];
        let frequencies = [2.0, 2.3];
        for profile in [HarmonicWindowProfile::Uniform, HarmonicWindowProfile::Triangular] {
            let window = HarmonicWindow::with_profile(0.17, 0.31, profile).unwrap();
            for proper in [false, true] {
                let phases = if proper {
                    [C64::ONE, -C64::ONE, C64::new(0.0, 1.0), C64::new(0.0, -1.0)]
                } else { [C64::ONE; 4] };
                let mut reference = [0.0; 6];
                let steps = 2048;
                for phase in phases {
                    for step in 0..=steps {
                        let u = step as f64 / steps as f64;
                        let t = window.start_s() + window.duration_s() * u;
                        let mut value = [0.0; 6];
                        let mut rate = [0.0; 6];
                        for (field, omega) in fields.iter().zip(frequencies) {
                            let (s, c) = (omega * t).sin_cos();
                            for k in 0..6 {
                                let a = field[k] * phase;
                                value[k] += a.re * c - a.im * s;
                                rate[k] += -omega * (a.re * s + a.im * c);
                            }
                        }
                        let simpson = if step == 0 || step == steps { 1.0 }
                            else if step % 2 == 0 { 2.0 } else { 4.0 };
                        let shutter = match profile {
                            HarmonicWindowProfile::Uniform => 1.0,
                            HarmonicWindowProfile::Triangular => 2.0 * (1.0 - (2.0*u-1.0).abs()),
                        };
                        for k in 0..3 {
                            let (a, b) = ((k + 1) % 3, (k + 2) % 3);
                            reference[k] += simpson * shutter
                                * (value[a]*value[b+3] - value[b]*value[a+3]);
                            reference[k+3] += simpson * shutter
                                * (rate[a]*value[b+3] + value[a]*rate[b+3]
                                    - rate[b]*value[a+3] - value[b]*rate[a+3]);
                        }
                    }
                }
                let actual = evaluate(&fields, &frequencies, window, proper);
                for k in 0..6 {
                    let expected = reference[k] / (12.0 * steps as f64);
                    assert!((actual[k]-expected).abs() < 3e-12,
                        "component {k}: {} != {expected}", actual[k]);
                }
            }
        }
    }

    #[test]
    fn g1_proper_slow_beat_rate_has_no_large_common_carrier() {
        // Two co-propagating unit-admittance waves, amplitudes 1 and 2i.
        // Proper mean Sz = 2.5 - 2 sin(delta*t); its derivative is a slow beat.
        let fields = [
            [C64::ONE, C64::ZERO, C64::ZERO, C64::ZERO, C64::ONE, C64::ZERO],
            [C64::new(0.0, 2.0), C64::ZERO, C64::ZERO,
             C64::ZERO, C64::new(0.0, 2.0), C64::ZERO],
        ];
        let delta = 0.25;
        let window = HarmonicWindow::new(0.17, 0.12).unwrap();
        let half_width = delta * window.duration_s() * 0.5;
        let sinc = half_width.sin() / half_width;
        let center = delta * (window.start_s() + window.duration_s() * 0.5);
        let expected = [2.5 - 2.0*sinc*center.sin(), -2.0*delta*sinc*center.cos()];
        let mut previous = None;
        for carrier in [0.0, 1e9, 1e12, 2_f64.powi(49)] {
            let actual = evaluate(&fields, &[carrier, carrier + delta], window, true);
            assert!((actual[2] - expected[0]).abs() < 2e-15);
            assert!((actual[5] - expected[1]).abs() < 2e-15);
            if let Some(prior) = previous { assert_eq!(actual, prior); }
            previous = Some(actual);
        }
        let monochromatic = evaluate(&fields[..1], &[1e15], window, true);
        assert_eq!(&monochromatic[3..], &[0.0; 3]);
        assert!((scale_momentum(Dd::from_f64(1e200), 1e-200, 1e-200)
            .unwrap().to_f64() / 1e-200 - 1.0).abs() < 2e-15);
        assert!(scale_momentum(Dd::from_f64(1e300), 1e100, 1e100).is_err());
        let kernels = [MomentumKernel {
            product: window.product_kernel(2.0, 2.0).unwrap(),
            mean_frequency: Dd::from_f64(2.0), half_difference: Dd::ZERO,
        }];
        assert!(matches!(integrate_momentum(1, &kernels, |_| fields[0],
            false, rotate_i, product, &mut || true), Err(MomentError::Cancelled)));
    }
}
