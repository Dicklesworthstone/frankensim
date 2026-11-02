//! Finite-time interference of original frequency-resolved Maxwell fields.
//!
//! For authored spectral lines, the physical fields are
//! E(t)=Re(sum_j E_j exp(+i omega_j t)), and likewise H. The shared harmonic
//! window integrates their REAL products, retaining both sum and difference
//! frequencies. An ultrashort shutter therefore does not acquire an implicit
//! optical-cycle average. Source amplitudes carry their phase at t=0; they are
//! peak line amplitudes, not a power spectral density or frequency-bin weight.
//!
//! The rectangular time-average definition and real/complex field distinction
//! follow TU Delft, Interference and Coherence, equations 295 and 303–306:
//! <https://interactivetextbooks.tudelft.nl/tn2421/content/Chap5_Interference/InterferenceCoherence_2022_01Clean.html>.
//! No stationary-source, random-phase or slowly varying envelope assumption is
//! inferred. Each frequency uses its original material/film/diffraction solve.

use super::*;
use fs_fft::harmonic_window::{HarmonicProduct, HarmonicWindow};
use fs_math::dd::Dd;

/// Whole-line-bank bounds for a query. Pair work counts N*F² for N spatial
/// points and F distinct frequencies. Logical numerical storage includes one
/// original beam's complete workspace, all retained frequency fields, temporal
/// kernels, a compensated accumulation plane and the output. Other borrowed
/// beams, earlier outputs and allocator overhead are outside this bound.
#[derive(Clone, Copy, Debug)]
pub struct AngularSpectrumExposureLimits {
    pub maximum_lines: usize,
    pub maximum_point_line_pairs: usize,
    pub maximum_point_frequency_pairs: usize,
    pub maximum_workspace_bytes: usize,
    pub maximum_fft_work: u64,
}

/// Rectangular-exposure averages of the real physical fields. Electric and
/// magnetic squared components are mean squares in the named SI units.
/// Poynting flux is <E(t) cross H(t)> in W/m² and retains its physical sign;
/// no monochromatic complex/reactive flux is assigned to a finite exposure.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AngularSpectrumExposurePoint {
    pub electric_mean_squared_v2_m2: [f64; 3],
    pub magnetic_mean_squared_a2_m2: [f64; 3],
    pub poynting_w_m2: [f64; 3],
}

/// Exposure-averaged observables on the original periodic transverse cell.
#[derive(Debug)]
pub struct AngularSpectrumExposureField {
    grid: AngularSpectrumGrid,
    points: Vec<AngularSpectrumExposurePoint>,
    normal_power_w: f64,
}
impl AngularSpectrumExposureField {
    pub fn grid(&self) -> AngularSpectrumGrid {
        self.grid
    }
    pub fn points(&self) -> &[AngularSpectrumExposurePoint] {
        &self.points
    }
    /// Mean signed normal power over the original periodic cell, in watts.
    pub fn normal_power_w(&self) -> f64 {
        self.normal_power_w
    }
}

/// Phase-related spectral lines through one physical stack. Lines may differ
/// in frequency, polarization and spatial profile, but must retain the same
/// original incident/exit/film material owners, layer order and thicknesses.
/// Equal-frequency lines are summed as fields before forming any products.
pub struct AngularSpectrumExposure<'b, 's, 'm> {
    anchor: &'b PlanarAngularSpectrum<'s, 'm>,
    lines: &'b [PlanarAngularSpectrum<'s, 'm>],
    frequencies: Vec<f64>,
    source_indices: Vec<usize>,
    group_offsets: Vec<usize>,
    workspace_bytes: usize,
    field_fft_work: u64,
    point_frequency_pairs: usize,
}

impl<'b, 's, 'm> AngularSpectrumExposure<'b, 's, 'm> {
    /// The anchor fixes the physical stack/grid and adds no illumination unless
    /// it is listed as a line. An empty bank is darkness. Each supplied line
    /// includes its physical amplitude; no normalization or spectrum is inferred.
    /// O(L log L) preparation groups source indices by increasing frequency,
    /// preserving the original line order within every frequency.
    pub fn try_new(
        anchor: &'b PlanarAngularSpectrum<'s, 'm>,
        lines: &'b [PlanarAngularSpectrum<'s, 'm>],
        limits: AngularSpectrumExposureLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AngularSpectrumError> {
        poll(&mut cancelled)?;
        if lines.len() > 4096 || lines.len() > limits.maximum_lines {
            return Err(bad("exposure line count exceeds its budget"));
        }
        // Reserve complete grouping capacity before any allocation. Offsets
        // retain L+1 slots even when several lines share a frequency.
        let retained_group_bytes = lines
            .len()
            .checked_mul(core::mem::size_of::<f64>() + 2 * core::mem::size_of::<usize>())
            .and_then(|v| v.checked_add(core::mem::size_of::<usize>()))
            .ok_or_else(|| bad("exposure workspace overflow"))?;
        let minimum_workspace = retained_group_bytes
            .checked_add(anchor.workspace_bytes())
            .ok_or_else(|| bad("exposure workspace overflow"))?;
        if minimum_workspace > limits.maximum_workspace_bytes {
            return Err(bad("exposure exceeds its workspace budget"));
        }
        let mut frequencies = Vec::new();
        frequencies
            .try_reserve_exact(lines.len())
            .map_err(|_| bad("exposure frequency allocation failed"))?;
        let mut source_indices = Vec::new();
        source_indices
            .try_reserve_exact(lines.len())
            .map_err(|_| bad("exposure source index allocation failed"))?;
        let mut group_offsets = Vec::new();
        group_offsets
            .try_reserve_exact(lines.len() + 1)
            .map_err(|_| bad("exposure group offset allocation failed"))?;
        let mut original_workspace = anchor.workspace_bytes();
        for (index, line) in lines.iter().enumerate() {
            poll(&mut cancelled)?;
            if line.grid() != anchor.grid()
                || !line.source().same_sources_and_geometry(anchor.source())
            {
                return Err(bad(
                    "exposure lines require the same original materials, layers and grid",
                ));
            }
            source_indices.push(index);
            original_workspace = original_workspace.max(line.workspace_bytes());
        }
        // In-place sorting needs no extra index allocation. The secondary
        // index gives stable within-frequency accumulation without relying on
        // sorting stability. Admission caps this bounded sort at 4096 lines.
        source_indices.sort_unstable_by(|left, right| {
            lines[*left]
                .source()
                .angular_frequency_rad_s()
                .total_cmp(&lines[*right].source().angular_frequency_rad_s())
                .then_with(|| left.cmp(right))
        });
        poll(&mut cancelled)?;
        for (offset, &index) in source_indices.iter().enumerate() {
            if offset % 64 == 0 {
                poll(&mut cancelled)?;
            }
            let omega = lines[index].source().angular_frequency_rad_s();
            if frequencies.last().copied() != Some(omega) {
                frequencies.push(omega);
                group_offsets.push(offset);
            }
        }
        group_offsets.push(lines.len());
        let n = anchor.modes.len();
        let f = frequencies.len();
        let line_pairs = n
            .checked_mul(lines.len())
            .ok_or_else(|| bad("exposure line work overflow"))?;
        let field_fft_work = anchor
            .field_fft_work()
            .checked_mul(f as u64)
            .ok_or_else(|| bad("exposure FFT work overflow"))?;
        let pairs = f
            .checked_mul(f)
            .ok_or_else(|| bad("exposure pair overflow"))?;
        let point_frequency_pairs = n
            .checked_mul(pairs)
            .ok_or_else(|| bad("exposure work overflow"))?;
        let workspace_bytes = n
            .checked_mul(core::mem::size_of::<[C64; 6]>())
            .and_then(|v| v.checked_mul(f))
            .and_then(|v| v.checked_add(f.checked_mul(core::mem::size_of::<Vec<[C64; 6]>>())?))
            .and_then(|v| {
                v.checked_add(pairs.checked_mul(core::mem::size_of::<HarmonicProduct>())?)
            })
            .and_then(|v| {
                v.checked_add(n.checked_mul(
                    core::mem::size_of::<[Dd; 12]>()
                        + core::mem::size_of::<AngularSpectrumExposurePoint>(),
                )?)
            })
            .and_then(|v| v.checked_add(retained_group_bytes))
            .and_then(|v| v.checked_add(original_workspace))
            .ok_or_else(|| bad("exposure workspace overflow"))?;
        if line_pairs > limits.maximum_point_line_pairs
            || point_frequency_pairs > limits.maximum_point_frequency_pairs
            || field_fft_work > limits.maximum_fft_work
            || workspace_bytes > limits.maximum_workspace_bytes
        {
            return Err(bad("exposure exceeds its work or workspace budget"));
        }
        Ok(Self {
            anchor,
            lines,
            frequencies,
            source_indices,
            group_offsets,
            workspace_bytes,
            field_fft_work,
            point_frequency_pairs,
        })
    }

    pub fn line_count(&self) -> usize {
        self.lines.len()
    }
    pub fn frequency_count(&self) -> usize {
        self.frequencies.len()
    }
    pub fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }
    pub fn field_fft_work(&self) -> u64 {
        self.field_fft_work
    }
    pub fn point_frequency_pairs(&self) -> usize {
        self.point_frequency_pairs
    }

    /// Reconstruct each original line once, cache O(F²) temporal kernels, and
    /// integrate all vector cross terms. O(L*N + F*N*log N + F²*N) work; trigonometry
    /// for the temporal window is O(F²), independent of the grid size. There is
    /// no time stepping, time aliasing or repeated constitutive solve.
    pub fn field(
        &self,
        side: AngularSpectrumSide,
        distance_from_face_m: f64,
        window: HarmonicWindow,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureField, AngularSpectrumError> {
        poll(&mut cancelled)?;
        if !distance_from_face_m.is_finite() || distance_from_face_m < 0.0 {
            return Err(bad("field distance must be finite and nonnegative"));
        }
        let f = self.frequencies.len();
        let mut kernels = Vec::new();
        kernels
            .try_reserve_exact(f * f)
            .map_err(|_| bad("exposure kernel allocation failed"))?;
        for &a in &self.frequencies {
            for &b in &self.frequencies {
                poll(&mut cancelled)?;
                kernels.push(
                    window
                        .product_kernel(a, b)
                        .map_err(|_| bad("exposure temporal kernel is unrepresentable"))?,
                );
            }
        }
        let mut fields = Vec::new();
        fields
            .try_reserve_exact(f)
            .map_err(|_| bad("exposure field allocation failed"))?;
        for group in 0..f {
            fields.push(self.frequency_field(group, side, distance_from_face_m, &mut cancelled)?);
        }
        let n = self.anchor.modes.len();
        let mut points = Vec::new();
        points
            .try_reserve_exact(n)
            .map_err(|_| bad("exposure result allocation failed"))?;
        let mut total = Dd::ZERO;
        for point in 0..n {
            poll(&mut cancelled)?;
            let values = integrate_point(f, &kernels, |i| fields[i][point], &mut cancelled)
                .map_err(|error| match error {
                    MomentError::Cancelled => AngularSpectrumError::Cancelled,
                    MomentError::Product(_) => bad("exposure field product is unrepresentable"),
                    MomentError::Unrepresentable => {
                        bad("exposure second moment is unrepresentable")
                    }
                })?;
            add(&mut total, values[8]);
            points.push(AngularSpectrumExposurePoint {
                electric_mean_squared_v2_m2: [values[0], values[1], values[2]],
                magnetic_mean_squared_a2_m2: [values[3], values[4], values[5]],
                poynting_w_m2: [values[6], values[7], values[8]],
            });
        }
        let grid = self.anchor.grid();
        let normal_power_w = total.to_f64() * grid.spacing_m[0] * grid.spacing_m[1];
        if !normal_power_w.is_finite() {
            return Err(bad("integrated exposure power is unrepresentable"));
        }
        poll(&mut cancelled)?;
        Ok(AngularSpectrumExposureField {
            grid,
            points,
            normal_power_w,
        })
    }

    // Sum equal-frequency sources BEFORE products, preserving complete nulls
    // and small residual fields instead of subtracting large source powers.
    fn frequency_field(
        &self,
        group: usize,
        side: AngularSpectrumSide,
        distance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<[C64; 6]>, AngularSpectrumError> {
        let n = self.anchor.modes.len();
        let mut sums = Vec::new();
        sums.try_reserve_exact(n)
            .map_err(|_| bad("exposure sum allocation failed"))?;
        for i in 0..n {
            if i % 256 == 0 {
                poll(cancelled)?;
            }
            sums.push([Dd::ZERO; 12]);
        }
        let indices =
            &self.source_indices[self.group_offsets[group]..self.group_offsets[group + 1]];
        for &index in indices {
            poll(cancelled)?;
            let line = &self.lines[index];
            let (electric, magnetic) = line.modal_field(side, distance, &mut *cancelled)?;
            for (i, ((sum, e), h)) in sums.iter_mut().zip(&electric).zip(&magnetic).enumerate() {
                if i % 64 == 0 {
                    poll(cancelled)?;
                }
                for (c, z) in e.iter().chain(h).enumerate() {
                    add(&mut sum[2 * c], z.re);
                    add(&mut sum[2 * c + 1], z.im);
                }
            }
        }
        let mut electric = Vec::new();
        let mut magnetic = Vec::new();
        electric
            .try_reserve_exact(n)
            .map_err(|_| bad("exposure group allocation failed"))?;
        magnetic
            .try_reserve_exact(n)
            .map_err(|_| bad("exposure group allocation failed"))?;
        for (i, sum) in sums.into_iter().enumerate() {
            if i % 64 == 0 {
                poll(cancelled)?;
            }
            let value: [C64; 6] =
                std::array::from_fn(|c| C64::new(sum[2 * c].to_f64(), sum[2 * c + 1].to_f64()));
            if value.iter().any(|&z| !finite(z)) {
                return Err(bad("coherent exposure field is unrepresentable"));
            }
            electric.push([value[0], value[1], value[2]]);
            magnetic.push([value[3], value[4], value[5]]);
        }
        let source = &self.lines[indices[0]];
        let synthesized = source.synthesize(electric, magnetic, &mut *cancelled)?;
        let mut field = Vec::new();
        field
            .try_reserve_exact(n)
            .map_err(|_| bad("exposure group allocation failed"))?;
        for (i, (e, h)) in synthesized
            .electric_v_per_m
            .into_iter()
            .zip(synthesized.magnetic_a_per_m)
            .enumerate()
        {
            if i % 64 == 0 {
                poll(cancelled)?;
            }
            field.push([e[0], e[1], e[2], h[0], h[1], h[2]]);
        }
        Ok(field)
    }
}

// Pure real-vector contraction, shared by the original-source exposure path
// and the analytic Maxwell-wave regressions below.
#[derive(Debug)]
enum MomentError {
    Cancelled,
    Product(fs_fft::harmonic_window::WindowError),
    Unrepresentable,
}

fn integrate_point(
    frequencies: usize,
    kernels: &[HarmonicProduct],
    field_at: impl Fn(usize) -> [C64; 6],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[f64; 9], MomentError> {
    let mut sum = [Dd::ZERO; 9];
    for i in 0..frequencies {
        let a = field_at(i);
        for j in 0..frequencies {
            if j % 64 == 0 && cancelled() {
                return Err(MomentError::Cancelled);
            }
            let b = field_at(j);
            let kernel = kernels[i * frequencies + j];
            for c in 0..3 {
                sum[c] = sum[c] + product(kernel, a[c], b[c])?;
                sum[c + 3] = sum[c + 3] + product(kernel, a[c + 3], b[c + 3])?;
                let (u, v) = ((c + 1) % 3, (c + 2) % 3);
                sum[c + 6] = sum[c + 6] + product(kernel, a[u], b[v + 3])?;
                sum[c + 6] = sum[c + 6] - product(kernel, a[v], b[u + 3])?;
            }
        }
    }
    let values = sum.map(Dd::to_f64);
    if values.iter().any(|x| !x.is_finite()) || values[..6].iter().any(|x| *x < 0.0) {
        return Err(MomentError::Unrepresentable);
    }
    Ok(values)
}

fn add(sum: &mut Dd, value: f64) {
    *sum = *sum + Dd::from_f64(value);
}
fn product(kernel: HarmonicProduct, a: C64, b: C64) -> Result<Dd, MomentError> {
    // The shared window uses exp(-i omega t); EM uses exp(+i omega t).
    kernel
        .mean_real_product_dd(a.conj(), b.conj())
        .map_err(MomentError::Product)
}

#[cfg(test)]
mod tests {
    use super::super::modes::fields;
    use super::*;

    fn wave(omega: f64, backwards: bool, transverse: [C64; 2]) -> [C64; 6] {
        let (e, h) = fields(
            0.0,
            0.0,
            C64::from_re(if backwards { -omega } else { omega }),
            omega,
            transverse,
        )
        .unwrap();
        [e[0], e[1], e[2], h[0], h[1], h[2]]
    }
    fn evaluate(waves: &[[C64; 6]], frequencies: &[f64], start: f64, duration: f64) -> [f64; 9] {
        let window = HarmonicWindow::new(start, duration).unwrap();
        let kernels = frequencies
            .iter()
            .flat_map(|&a| {
                frequencies
                    .iter()
                    .map(move |&b| window.product_kernel(a, b).unwrap())
            })
            .collect::<Vec<_>>();
        integrate_point(frequencies.len(), &kernels, |i| waves[i], &mut || false).unwrap()
    }

    #[test]
    fn g1_real_maxwell_exposure_matches_direct_simpson_integration() {
        let frequencies = [4.0, 5.3];
        let waves = [
            wave(4.0, false, [C64::new(0.7, -0.3), C64::new(0.2, 0.5)]),
            wave(5.3, true, [C64::new(-0.4, 0.2), C64::new(0.1, -0.6)]),
        ];
        let (start, duration) = (0.17, 0.83);
        let actual = evaluate(&waves, &frequencies, start, duration);
        let mut expected = [0.0; 9];
        let n = 4096;
        for k in 0..=n {
            let t = start + duration * k as f64 / n as f64;
            let mut real = [0.0; 6];
            for (wave, omega) in waves.iter().zip(frequencies) {
                let (s, c) = (omega * t).sin_cos();
                for j in 0..6 {
                    real[j] += wave[j].re * c - wave[j].im * s;
                }
            }
            let weight = if k == 0 || k == n {
                1.0
            } else if k % 2 == 0 {
                2.0
            } else {
                4.0
            };
            for j in 0..3 {
                expected[j] += weight * real[j] * real[j];
                expected[j + 3] += weight * real[j + 3] * real[j + 3];
                let (u, v) = ((j + 1) % 3, (j + 2) % 3);
                expected[j + 6] += weight * (real[u] * real[v + 3] - real[v] * real[u + 3]);
            }
        }
        for (actual, expected) in actual
            .into_iter()
            .zip(expected.map(|v| v / (3.0 * n as f64)))
        {
            assert!((actual - expected).abs() < 3e-13, "{actual} != {expected}");
        }
    }

    #[test]
    fn g1_distinct_frequency_cancellation_retains_weak_beating_power() {
        let frequencies = [10.0, 10.0001];
        let waves = [
            wave(frequencies[0], false, [C64::from_re(1e3), C64::ZERO]),
            wave(frequencies[1], false, [C64::from_re(-1e3), C64::ZERO]),
        ];
        let actual = evaluate(&waves, &frequencies, 0.3, 0.7);
        let mut expected = 0.0;
        let n = 8192;
        for k in 0..=n {
            let t = 0.3 + 0.7 * k as f64 / n as f64;
            // Independent stable identity for cos(wa*t)-cos(wb*t).
            let e = 2e3
                * (((frequencies[0] + frequencies[1]) * 0.5) * t).sin()
                * (((frequencies[1] - frequencies[0]) * 0.5) * t).sin();
            let weight = if k == 0 || k == n {
                1.0
            } else if k % 2 == 0 {
                2.0
            } else {
                4.0
            };
            expected += weight * e * e;
        }
        expected /= 3.0 * n as f64;
        assert!(expected > 0.0);
        assert!(
            (actual[0] / expected - 1.0).abs() < 2e-9,
            "{} != {expected}",
            actual[0]
        );
        // The original Maxwell reconstruction obtains H via k cross E / omega;
        // its rounded multiplication/division need not equal E bit for bit.
        for component in [4, 8] {
            assert!((actual[component] / actual[0] - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn g1_short_exposures_keep_phase_polarization_and_reflected_flux_sign() {
        let sine = wave(10.0, false, [C64::new(0.0, 1.0), C64::ZERO]);
        let tiny = evaluate(&[sine], &[10.0], 0.0, 1e-10);
        assert!((tiny[0] / (1e-18 / 3.0) - 1.0).abs() < 1e-12);
        let reflected = [
            wave(10.0, true, [C64::ONE, C64::ZERO]),
            wave(11.0, true, [C64::ZERO, C64::ONE]),
        ];
        let long = evaluate(&reflected, &[10.0, 11.0], 0.0, core::f64::consts::TAU);
        assert!((long[8] + 1.0).abs() < 1e-14);
        let mut cancelled = || true;
        let kernel = HarmonicWindow::new(0.0, 1.0)
            .unwrap()
            .product_kernel(10.0, 10.0)
            .unwrap();
        assert!(matches!(
            integrate_point(1, &[kernel], |_| sine, &mut cancelled),
            Err(MomentError::Cancelled)
        ));
    }
}
