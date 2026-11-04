//! Exact rectangle averages of the original transverse Fourier fields.
//!
//! For exp(-i k.r), every E_i cross H_j* term is weighted by
//! exp(-i (k_i-k_j).center) sinc(dx*width/2) sinc(dy*height/2).
//! The shared harmonic-window owner supplies those integrals in double-double
//! precision, with metres and radians/metre replacing seconds and radians/second.
//! No stochastic source assumption is made: its difference-only kernel is
//! simply the same bilinear integral as the cycle-mean spatial field product.

use super::{AngularSpectrumError, AngularSpectrumGrid, bad, poll};
use fs_fft::harmonic_window::HarmonicWindow;
use fs_math::{c64::C64, dd::Dd};

/// Axis-aligned detector area in the original plane's physical x/y coordinates.
/// Origin and extent are retained separately so a small detector is not rounded
/// into a large coordinate. Its response is uniform across the rectangle.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AngularSpectrumAperture {
    origin_m: [f64; 2],
    extent_m: [f64; 2],
    area_m2: f64,
}

impl AngularSpectrumAperture {
    pub fn try_new(origin_m: [f64; 2], extent_m: [f64; 2]) -> Result<Self, AngularSpectrumError> {
        for axis in 0..2 {
            HarmonicWindow::new(origin_m[axis], extent_m[axis])
                .map_err(|_| bad("aperture needs finite origins and positive finite extents"))?;
        }
        let area_m2 = extent_m[0] * extent_m[1];
        if !area_m2.is_finite() || area_m2 <= 0.0 {
            return Err(bad("aperture area is unrepresentable"));
        }
        Ok(Self {
            origin_m,
            extent_m,
            area_m2,
        })
    }

    pub fn origin_m(self) -> [f64; 2] {
        self.origin_m
    }
    pub fn extent_m(self) -> [f64; 2] {
        self.extent_m
    }
    pub fn area_m2(self) -> f64 {
        self.area_m2
    }
}

/// Logical query bounds, checked before allocation or modal reconstruction.
/// For N=W*H original modes, `maximum_mode_pairs` bounds N*(W+H) axis
/// contractions of a three-component magnetic vector. Another N vector cross
/// products complete the query; W²+H² scalar kernel evaluations build the two
/// axis caches. Numerical buffer storage includes the original beam workspace,
/// both axis-pair caches, and N three-component complex double-double vectors.
/// Caller inputs, previously retained outputs and allocator overhead are excluded.
#[derive(Clone, Copy, Debug)]
pub struct AngularSpectrumApertureLimits {
    pub maximum_mode_pairs: usize,
    pub maximum_workspace_bytes: usize,
}

/// Cycle-mean electromagnetic transport through the declared rectangle.
/// Flux is the area average of the full vector E cross H* / 2, in W/m².
/// Normal power integrates its Cartesian z component, in watts; reflected
/// radiation and local energy backflow keep their negative signs.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AngularSpectrumApertureReading {
    pub mean_poynting_w_m2: [f64; 3],
    pub normal_power_w: f64,
}

pub(super) struct PreparedAperture {
    grid: AngularSpectrumGrid,
    aperture: AngularSpectrumAperture,
    x: Vec<[Dd; 2]>,
    y: Vec<[Dd; 2]>,
}

impl PreparedAperture {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn prepare(
        grid: AngularSpectrumGrid,
        aperture: AngularSpectrumAperture,
        beam_workspace_bytes: usize,
        limits: AngularSpectrumApertureLimits,
        kx: impl Fn(usize) -> f64,
        ky: impl Fn(usize) -> f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AngularSpectrumError> {
        poll(cancelled)?;
        let n = grid
            .width
            .checked_mul(grid.height)
            .ok_or_else(|| bad("aperture mode count overflow"))?;
        let pairs = grid
            .width
            .checked_add(grid.height)
            .and_then(|axes| n.checked_mul(axes))
            .ok_or_else(|| bad("aperture pair count overflow"))?;
        let axis_pairs = grid
            .width
            .checked_mul(grid.width)
            .and_then(|x| {
                grid.height
                    .checked_mul(grid.height)
                    .and_then(|y| x.checked_add(y))
            })
            .ok_or_else(|| bad("aperture axis cache extent overflow"))?;
        let workspace_bytes = axis_pairs
            .checked_mul(core::mem::size_of::<[Dd; 2]>())
            .and_then(|bytes| {
                n.checked_mul(core::mem::size_of::<[[Dd; 2]; 3]>())
                    .and_then(|scratch| bytes.checked_add(scratch))
            })
            .and_then(|bytes| bytes.checked_add(beam_workspace_bytes))
            .ok_or_else(|| bad("aperture workspace extent overflow"))?;
        if pairs > limits.maximum_mode_pairs || workspace_bytes > limits.maximum_workspace_bytes {
            return Err(bad("aperture exceeds its pair work or workspace budget"));
        }
        // Cache ACTUAL admitted axis-wavevector pairs. Deriving a difference
        // from a rounded fundamental instead would change resolved phases at
        // a distant detector. The Cartesian product needs no more trigonometry.
        let x = axis_weights(
            grid.width,
            aperture.origin_m[0],
            aperture.extent_m[0],
            kx,
            cancelled,
        )?;
        let y = axis_weights(
            grid.height,
            aperture.origin_m[1],
            aperture.extent_m[1],
            ky,
            cancelled,
        )?;
        poll(cancelled)?;
        Ok(Self {
            grid,
            aperture,
            x,
            y,
        })
    }

    #[allow(clippy::float_cmp)] // Exact binary scale/component round trips are admission checks.
    pub(super) fn measure(
        &self,
        electric: &[[C64; 3]],
        magnetic: &[[C64; 3]],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AngularSpectrumApertureReading, AngularSpectrumError> {
        poll(cancelled)?;
        let n = self.grid.width * self.grid.height;
        if electric.len() != n || magnetic.len() != n {
            return Err(bad("aperture fields must cover the original modal grid"));
        }
        let restore = magnetic_scale(magnetic, cancelled)?;
        let reduction = 1.0 / restore;
        // W(i,j)=Wx(ix,jx)*Wy(iy,jy), so apply the two exact axis
        // kernels to H* before contracting with E. Keep every complex part
        // in DD: rounding an intermediate field to C64 would erase weak
        // detector-null power. This local scratch also leaves a prepared
        // aperture reusable without mutation after cancellation or failure.
        let mut y_magnetic = Vec::new();
        y_magnetic
            .try_reserve_exact(n)
            .map_err(|_| bad("aperture field workspace allocation failed"))?;
        for row in 0..self.grid.height {
            for column in 0..self.grid.width {
                poll(cancelled)?;
                let mut field = [[Dd::ZERO; 2]; 3];
                for source_row in 0..self.grid.height {
                    if source_row % 32 == 0 {
                        poll(cancelled)?;
                    }
                    let weight = self.y[row * self.grid.height + source_row];
                    let h = magnetic[source_row * self.grid.width + column];
                    for (output, value) in field.iter_mut().zip(h) {
                        add_weighted(
                            output,
                            [
                                Dd::from_f64(value.re * reduction),
                                Dd::from_f64(-value.im * reduction),
                            ],
                            weight,
                        );
                    }
                }
                y_magnetic.push(field);
            }
        }
        let mut sum = [Dd::ZERO; 3];
        for (i, e) in electric.iter().enumerate() {
            poll(cancelled)?;
            let row = i / self.grid.width;
            let column = i % self.grid.width;
            let mut field = [[Dd::ZERO; 2]; 3];
            for source_column in 0..self.grid.width {
                if source_column % 32 == 0 {
                    poll(cancelled)?;
                }
                let weight = self.x[column * self.grid.width + source_column];
                let h = y_magnetic[row * self.grid.width + source_column];
                for (output, value) in field.iter_mut().zip(h) {
                    add_weighted(output, value, weight);
                }
            }
            for (axis, (a, b)) in [(1, 2), (2, 0), (0, 1)].into_iter().enumerate() {
                sum[axis] = sum[axis] + Dd::from_f64(e[a].re) * field[b][0]
                    - Dd::from_f64(e[a].im) * field[b][1]
                    - Dd::from_f64(e[b].re) * field[a][0]
                    + Dd::from_f64(e[b].im) * field[a][1];
            }
        }
        let half = Dd::from_f64(0.5);
        let mut mean = sum.map(|value| value * half);
        if restore != 1.0 {
            for (value, original) in mean.iter_mut().zip(sum) {
                let (hi, lo) = (value.hi * restore, value.lo * restore);
                if !hi.is_finite() || !lo.is_finite() {
                    return Err(bad("aperture Poynting mean is unrepresentable"));
                }
                // Refuse a lost DD component even if the rounded mean survived.
                if value.hi * 2.0 != original.hi || value.lo * 2.0 != original.lo {
                    return Err(bad("aperture scaled Poynting mean underflows"));
                }
                *value = Dd::from_pair(hi, lo);
            }
        }
        let mean_poynting_w_m2 = mean.map(Dd::to_f64);
        let normal_power_w = (mean[2] * Dd::from_f64(self.aperture.area_m2)).to_f64();
        if mean_poynting_w_m2
            .iter()
            .chain(&[normal_power_w])
            .any(|value| !value.is_finite())
        {
            return Err(bad("aperture Poynting transport is unrepresentable"));
        }
        poll(cancelled)?;
        Ok(AngularSpectrumApertureReading {
            mean_poynting_w_m2,
            normal_power_w,
        })
    }
}

/// |Wx|,|Wy| <= 1 bound a complex component by 4*N times its input maximum.
/// Keep another factor-four headroom for finite DD sums. Binary scaling is
/// conditional, exact on every admitted source component, and restored only
/// after the electric contraction; ordinary field arithmetic is unchanged.
#[allow(clippy::float_cmp)] // Exact binary scale/component round trips are admission checks.
fn magnetic_scale(
    magnetic: &[[C64; 3]],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, AngularSpectrumError> {
    let mut largest = 0.0_f64;
    for fields in magnetic.chunks(128) {
        poll(cancelled)?;
        for value in fields.iter().flatten().flat_map(|v| [v.re, v.im]) {
            if !value.is_finite() {
                return Err(bad("aperture magnetic field is unrepresentable"));
            }
            largest = largest.max(value.abs());
        }
    }
    let limit = f64::MAX / (16.0 * magnetic.len() as f64);
    let mut restore = 1.0;
    while largest > limit {
        largest *= 0.5;
        restore *= 2.0;
    }
    if restore != 1.0 {
        for fields in magnetic.chunks(128) {
            poll(cancelled)?;
            for value in fields.iter().flatten().flat_map(|v| [v.re, v.im]) {
                if (value / restore) * restore != value {
                    return Err(bad("aperture magnetic scaling loses a nonzero component"));
                }
            }
        }
    }
    Ok(restore)
}

fn axis_weights(
    count: usize,
    origin: f64,
    extent: f64,
    wavenumber: impl Fn(usize) -> f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<[Dd; 2]>, AngularSpectrumError> {
    let window = HarmonicWindow::new(origin, extent)
        .map_err(|_| bad("aperture coordinate is unrepresentable"))?;
    let mut weights = Vec::new();
    weights
        .try_reserve_exact(count * count)
        .map_err(|_| bad("aperture kernel allocation failed"))?;
    for i in 0..count {
        for j in 0..count {
            if j % 32 == 0 {
                poll(cancelled)?;
            }
            let kernel = window
                .proper_product_kernel(wavenumber(i), wavenumber(j))
                .map_err(|_| bad("aperture spatial phase is unrepresentable"))?;
            // Recover the complex integral without rounding its low part.
            // The proper kernel includes 1/2; the final E cross H* contraction
            // applies that physical cycle-mean normalization only once.
            let real = kernel
                .mean_real_product_dd(C64::ONE, C64::ONE)
                .map_err(|_| bad("aperture spatial kernel is unrepresentable"))?;
            let negative_imaginary = kernel
                .mean_real_product_dd(C64::new(0.0, 1.0), C64::ONE)
                .map_err(|_| bad("aperture spatial kernel is unrepresentable"))?;
            weights.push([
                real * Dd::from_f64(2.0),
                negative_imaginary * Dd::from_f64(-2.0),
            ]);
        }
    }
    Ok(weights)
}

fn add_weighted(output: &mut [Dd; 2], value: [Dd; 2], weight: [Dd; 2]) {
    output[0] = output[0] + value[0] * weight[0] - value[1] * weight[1];
    output[1] = output[1] + value[0] * weight[1] + value[1] * weight[0];
}

#[cfg(test)]
mod tests {
    use super::*;

    fn limits() -> AngularSpectrumApertureLimits {
        AngularSpectrumApertureLimits {
            maximum_mode_pairs: 4096,
            maximum_workspace_bytes: 1 << 20,
        }
    }

    fn close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "{actual:e} != {expected:e}, tolerance {tolerance:e}"
        );
    }

    #[test]
    fn g1_detector_width_and_position_resolve_two_wave_fringe_visibility() {
        let q = 0.75;
        let admittance = 0.8;
        let grid = AngularSpectrumGrid {
            width: 2,
            height: 1,
            spacing_m: [1.0; 2],
        };
        let e = [[C64::ZERO, C64::ONE, C64::ZERO]; 2];
        let h = [
            [C64::from_re(-admittance), C64::ZERO, C64::from_re(q)],
            [C64::from_re(-admittance), C64::ZERO, C64::from_re(-q)],
        ];
        for center in [-0.37, 0.0, 0.83] {
            for width in [0.013, 0.7, core::f64::consts::PI / q, 8.3] {
                let aperture =
                    AngularSpectrumAperture::try_new([center - width / 2.0, -2.0], [width, 3.0])
                        .unwrap();
                let prepared = PreparedAperture::prepare(
                    grid,
                    aperture,
                    0,
                    limits(),
                    |i| if i == 0 { q } else { -q },
                    |_| 0.0,
                    &mut || false,
                )
                .unwrap();
                let actual = prepared.measure(&e, &h, &mut || false).unwrap();
                // Ey = 2 cos(q*x), Hx = -Y Ey. The detector integrates
                // the intensity, not the mean amplitude squared.
                let expected =
                    admittance * (1.0 + (2.0 * q * center).cos() * (q * width).sin() / (q * width));
                close(actual.mean_poynting_w_m2[2], expected, 8e-15);
                close(actual.normal_power_w, expected * width * 3.0, 2e-13);
                close(actual.mean_poynting_w_m2[0], 0.0, 1e-15);
                close(actual.mean_poynting_w_m2[1], 0.0, 0.0);
            }
        }
    }

    #[test]
    fn g1_subwavelength_dark_fringe_retains_positive_aperture_power() {
        let grid = AngularSpectrumGrid {
            width: 2,
            height: 1,
            spacing_m: [1.0; 2],
        };
        let e = [
            [C64::ZERO, C64::ONE, C64::ZERO],
            [C64::ZERO, -C64::ONE, C64::ZERO],
        ];
        let h = [
            [-C64::ONE, C64::ZERO, C64::ONE],
            [C64::ONE, C64::ZERO, C64::ONE],
        ];
        for width in [1e-6, 1e-8, 1e-10] {
            let aperture =
                AngularSpectrumAperture::try_new([-width / 2.0, 0.0], [width, 2.0]).unwrap();
            let prepared = PreparedAperture::prepare(
                grid,
                aperture,
                0,
                limits(),
                |i| if i == 0 { 1.0 } else { -1.0 },
                |_| 0.0,
                &mut || false,
            )
            .unwrap();
            let actual = prepared.measure(&e, &h, &mut || false).unwrap();
            // <2 sin²(x)> = 1-sinc(width); evaluate the independent series
            // directly, without subtracting two rounded values near one.
            let expected = width * width / 6.0 * (1.0 - width * width / 20.0);
            assert!(actual.normal_power_w > 0.0);
            close(actual.mean_poynting_w_m2[2], expected, expected * 2e-10);
            close(
                actual.normal_power_w,
                expected * width * 2.0,
                expected * width * 4e-10,
            );
        }
    }

    #[test]
    fn g1_two_axis_dark_fringe_retains_positive_aperture_power() {
        let grid = AngularSpectrumGrid {
            width: 2,
            height: 2,
            spacing_m: [1.0; 2],
        };
        let signs = [C64::ONE, -C64::ONE, -C64::ONE, C64::ONE];
        let e = signs.map(|value| [C64::ZERO, value, C64::ZERO]);
        let h = signs.map(|value| [-value, C64::ZERO, C64::ZERO]);
        for width_x in [1e-4, 1e-8, 1e-10] {
            let width_y = 1.7 * width_x;
            let aperture = AngularSpectrumAperture::try_new(
                [-width_x / 2.0, -width_y / 2.0],
                [width_x, width_y],
            )
            .unwrap();
            let prepared = PreparedAperture::prepare(
                grid,
                aperture,
                0,
                limits(),
                |i| if i == 0 { 1.0 } else { -1.0 },
                |i| if i == 0 { 1.0 } else { -1.0 },
                &mut || false,
            )
            .unwrap();
            let actual = prepared.measure(&e, &h, &mut || false).unwrap();
            // Ey=-4 sin(x) sin(y), Hx=-Ey, so <Sz> is exactly
            // 2*(1-sinc(width_x))*(1-sinc(width_y)). Use independent
            // small-window series: both axis contractions must retain their
            // dark-fringe residual before these residuals are multiplied.
            let one_minus_sinc = |width: f64| {
                let square = width * width;
                square / 6.0 * (1.0 - square / 20.0 + square * square / 840.0)
            };
            let expected = 2.0 * one_minus_sinc(width_x) * one_minus_sinc(width_y);
            let expected_power = expected * aperture.area_m2();
            assert!(actual.mean_poynting_w_m2[2] > 0.0);
            assert!(actual.normal_power_w > 0.0);
            assert_eq!(actual.mean_poynting_w_m2[..2], [0.0; 2]);
            close(actual.mean_poynting_w_m2[2], expected, expected * 3e-10);
            close(
                actual.normal_power_w,
                expected_power,
                expected_power * 3e-10,
            );
        }
    }

    // Independent plane waves satisfy k.E=0 and H=k cross E/(omega*mu).
    // Complex kz covers a decaying evanescent component; negative real kz
    // covers reflected flux. This test never calls the production mode helper.
    fn wave(kx: f64, ky: f64, kz: C64, ex: C64, ey: C64) -> ([C64; 3], [C64; 3]) {
        let e = [ex, ey, -(ex.scale(kx) + ey.scale(ky)) / kz];
        let h = [
            e[2].scale(ky) - kz * ey,
            kz * ex - e[2].scale(kx),
            ey.scale(kx) - ex.scale(ky),
        ];
        (e, h)
    }

    #[test]
    fn g1_vector_and_evanescent_interference_match_direct_area_quadrature() {
        let grid = AngularSpectrumGrid {
            width: 2,
            height: 2,
            spacing_m: [1.0; 2],
        };
        let kx: [f64; 2] = [0.0, 1.3];
        let ky: [f64; 2] = [0.0, -0.4];
        let mut e = Vec::new();
        let mut h = Vec::new();
        for y in ky {
            for x in kx {
                let q2 = 1.0 - x * x - y * y;
                let kz = if q2 < 0.0 {
                    C64::new(0.0, -(-q2).sqrt())
                } else {
                    C64::from_re(-q2.sqrt())
                };
                let (electric, magnetic) =
                    wave(x, y, kz, C64::new(0.3 + x, -0.2), C64::new(-0.4, 0.7 + y));
                e.push(electric);
                h.push(magnetic);
            }
        }
        let aperture = AngularSpectrumAperture::try_new([-0.6, 0.3], [1.3, 0.9]).unwrap();
        let prepared = PreparedAperture::prepare(
            grid,
            aperture,
            0,
            limits(),
            |i| kx[i],
            |i| ky[i],
            &mut || false,
        )
        .unwrap();
        let actual = prepared.measure(&e, &h, &mut || false).unwrap();
        // Tensor-product Simpson integration of independently synthesized REAL
        // and imaginary vector fields, before their local Poynting product.
        let intervals = 160;
        let mut expected = [0.0; 3];
        for row in 0..=intervals {
            let y = aperture.origin_m[1] + aperture.extent_m[1] * row as f64 / intervals as f64;
            let wy = if row == 0 || row == intervals {
                1.0
            } else if row % 2 == 0 {
                2.0
            } else {
                4.0
            };
            for col in 0..=intervals {
                let x = aperture.origin_m[0] + aperture.extent_m[0] * col as f64 / intervals as f64;
                let wx = if col == 0 || col == intervals {
                    1.0
                } else if col % 2 == 0 {
                    2.0
                } else {
                    4.0
                };
                let mut electric = [C64::ZERO; 3];
                let mut magnetic = [C64::ZERO; 3];
                for i in 0..4 {
                    let phase = -(kx[i % 2] * x + ky[i / 2] * y);
                    let phase = C64::new(phase.cos(), phase.sin());
                    for c in 0..3 {
                        electric[c] = electric[c] + e[i][c] * phase;
                        magnetic[c] = magnetic[c] + h[i][c] * phase;
                    }
                }
                for (axis, (a, b)) in [(1, 2), (2, 0), (0, 1)].into_iter().enumerate() {
                    expected[axis] += wx
                        * wy
                        * 0.5
                        * (electric[a] * magnetic[b].conj() - electric[b] * magnetic[a].conj()).re;
                }
            }
        }
        for axis in 0..3 {
            expected[axis] /= 9.0 * (intervals * intervals) as f64;
            close(actual.mean_poynting_w_m2[axis], expected[axis], 3e-9);
        }
        assert!(actual.normal_power_w < 0.0);
        close(actual.normal_power_w, expected[2] * aperture.area_m2, 5e-9);
    }

    #[test]
    fn g3_whole_cell_orthogonality_and_detector_partition_preserve_flux() {
        let grid = AngularSpectrumGrid {
            width: 2,
            height: 2,
            spacing_m: [2.0, 3.0],
        };
        let kx = [0.0, -core::f64::consts::TAU / 4.0];
        let ky = [0.0, -core::f64::consts::TAU / 6.0];
        let e = [[C64::new(0.1, -0.2), C64::new(0.3, 0.1), C64::ONE]; 4];
        let h = [[C64::new(0.2, 0.3), C64::new(-0.4, 0.5), C64::ZERO]; 4];
        let measure = |origin, extent| {
            PreparedAperture::prepare(
                grid,
                AngularSpectrumAperture::try_new(origin, extent).unwrap(),
                0,
                limits(),
                |i| kx[i],
                |i| ky[i],
                &mut || false,
            )
            .unwrap()
            .measure(&e, &h, &mut || false)
            .unwrap()
        };
        let whole = measure([-0.7, 1.2], [4.0, 6.0]);
        for (axis, (a, b)) in [(1, 2), (2, 0), (0, 1)].into_iter().enumerate() {
            let diagonal: f64 = e
                .iter()
                .zip(h)
                .map(|(e, h)| 0.5 * (e[a] * h[b].conj() - e[b] * h[a].conj()).re)
                .sum();
            close(whole.mean_poynting_w_m2[axis], diagonal, 2e-14);
        }
        let left = measure([-0.7, 1.2], [1.1, 6.0]);
        let right = measure([0.4, 1.2], [2.9, 6.0]);
        close(
            left.normal_power_w + right.normal_power_w,
            whole.normal_power_w,
            2e-13,
        );
        let dark = PreparedAperture::prepare(
            grid,
            AngularSpectrumAperture::try_new([0.0; 2], [1.0; 2]).unwrap(),
            0,
            limits(),
            |i| kx[i],
            |i| ky[i],
            &mut || false,
        )
        .unwrap()
        .measure(&[[C64::ZERO; 3]; 4], &[[C64::ZERO; 3]; 4], &mut || false)
        .unwrap();
        assert_eq!(dark.mean_poynting_w_m2, [0.0; 3]);
        assert_eq!(dark.normal_power_w, 0.0);
    }

    #[test]
    fn g3_separable_transport_matches_all_pairs_at_exact_axis_work_limit() {
        let grid = AngularSpectrumGrid {
            width: 4,
            height: 8,
            spacing_m: [1.0; 2],
        };
        let n = grid.width * grid.height;
        let work = n * (grid.width + grid.height);
        assert!(work < n * n);
        let kx = [-0.8, 0.17, 1.3, -0.23];
        let ky = [0.0, 0.14, -0.23, 1.1, -0.75, 0.37, -0.41, 1.23];
        let e: Vec<[C64; 3]> = (0..n)
            .map(|i| {
                core::array::from_fn(|c| {
                    C64::new(0.13 * (i + c + 1) as f64, -0.07 * (i * c + 2) as f64)
                })
            })
            .collect();
        let h: Vec<[C64; 3]> = (0..n)
            .map(|i| {
                core::array::from_fn(|c| {
                    C64::new(-0.11 * (i + 2 * c + 1) as f64, 0.19 * (i + c + 3) as f64)
                })
            })
            .collect();
        let exact = AngularSpectrumApertureLimits {
            maximum_mode_pairs: work,
            maximum_workspace_bytes: (grid.width * grid.width + grid.height * grid.height)
                * core::mem::size_of::<[Dd; 2]>()
                + n * core::mem::size_of::<[[Dd; 2]; 3]>(),
        };
        for origin in [[-0.41, 0.27], [1e13 + 0.25, -3e12 + 0.5]] {
            let aperture = AngularSpectrumAperture::try_new(origin, [1.1, 0.9]).unwrap();
            // Refuse before evaluating any mode or allocating either cache.
            assert!(
                PreparedAperture::prepare(
                    grid,
                    aperture,
                    0,
                    AngularSpectrumApertureLimits {
                        maximum_mode_pairs: work - 1,
                        ..exact
                    },
                    |_| panic!("work refused after reading a wavevector"),
                    |_| panic!("work refused after reading a wavevector"),
                    &mut || false,
                )
                .is_err()
            );
            let prepared = PreparedAperture::prepare(
                grid,
                aperture,
                0,
                exact,
                |i| kx[i],
                |i| ky[i],
                &mut || false,
            )
            .unwrap();
            let actual = prepared.measure(&e, &h, &mut || false).unwrap();
            // Independent all-pairs reduction of the original E_i cross H_j*.
            // Nonuniform axes and distant windows expose phase/sign/index swaps.
            let mut expected = [Dd::ZERO; 3];
            for (i, e) in e.iter().enumerate() {
                for (j, h) in h.iter().enumerate() {
                    let x = prepared.x[(i % grid.width) * grid.width + j % grid.width];
                    let y = prepared.y[(i / grid.width) * grid.height + j / grid.width];
                    let weight = [x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0]];
                    let product = |e: C64, h: C64| {
                        let [er, ei, hr, hi] = [e.re, e.im, h.re, h.im].map(Dd::from_f64);
                        (er * hr + ei * hi) * weight[0] - (ei * hr - er * hi) * weight[1]
                    };
                    for (axis, (a, b)) in [(1, 2), (2, 0), (0, 1)].into_iter().enumerate() {
                        expected[axis] = expected[axis] + product(e[a], h[b]) - product(e[b], h[a]);
                    }
                }
            }
            for (actual, expected) in actual.mean_poynting_w_m2.into_iter().zip(expected) {
                let expected = (expected * Dd::from_f64(0.5)).to_f64();
                close(actual, expected, 1e-14 * expected.abs().max(1.0));
            }
            close(
                actual.normal_power_w,
                (expected[2] * Dd::from_f64(0.5) * Dd::from_f64(aperture.area_m2)).to_f64(),
                1e-14 * actual.normal_power_w.abs().max(1.0),
            );
            // Both axis passes can be cancelled without changing later results.
            for stop in [3, 2 * n + 4] {
                let mut polls = 0;
                assert!(matches!(
                    prepared.measure(&e, &h, &mut || {
                        polls += 1;
                        polls == stop
                    }),
                    Err(AngularSpectrumError::Cancelled)
                ));
                assert_eq!(prepared.measure(&e, &h, &mut || false).unwrap(), actual);
            }
        }
    }

    #[test]
    fn g1_scaled_magnetic_action_preserves_finite_all_pairs_transport() {
        let grid = AngularSpectrumGrid {
            width: 2,
            height: 2,
            spacing_m: [1.0; 2],
        };
        let kx = [1.0, -1.0];
        let ky = [0.7, -0.7];
        let aperture = AngularSpectrumAperture::try_new([-0.12, -0.125], [0.5, 0.25]).unwrap();
        let prepared = PreparedAperture::prepare(
            grid,
            aperture,
            0,
            limits(),
            |i| kx[i],
            |i| ky[i],
            &mut || false,
        )
        .unwrap();
        let e: [[C64; 3]; 4] = core::array::from_fn(|i| {
            [
                C64::new(0.4e-308, -0.2e-308),
                C64::new(1e-308, (i + 1) as f64 * 0.2e-308),
                C64::new(-0.3e-308, 0.1e-308),
            ]
        });
        let h = [[
            C64::new(-1e308, 0.2e308),
            C64::new(0.3e308, -0.1e308),
            C64::new(0.2e308, 0.15e308),
        ]; 4];
        // This first-axis Hx sum exceeds f64 even though every E_i H_j*
        // product and the actual vector flux remain finite.
        let half_phase = (ky[0] - ky[1]) * aperture.extent_m[1] * 0.5;
        assert!(!(h[0][0].re * (1.0 + half_phase.sin() / half_phase)).is_finite());
        let actual = prepared.measure(&e, &h, &mut || false).unwrap();
        // Independent finite all-pairs reference: integrate each exponential
        // directly and multiply E by H before applying its window weight.
        let window = |delta: f64, axis: usize| {
            let half_phase = delta * aperture.extent_m[axis] * 0.5;
            let sinc = if delta == 0.0 {
                1.0
            } else {
                half_phase.sin() / half_phase
            };
            let phase = -delta * (aperture.origin_m[axis] + aperture.extent_m[axis] * 0.5);
            C64::new(sinc * phase.cos(), sinc * phase.sin())
        };
        let mut expected = [0.0; 3];
        for (i, e) in e.iter().enumerate() {
            for (j, h) in h.iter().enumerate() {
                let weight = window(kx[i % 2] - kx[j % 2], 0) * window(ky[i / 2] - ky[j / 2], 1);
                for (axis, (a, b)) in [(1, 2), (2, 0), (0, 1)].into_iter().enumerate() {
                    expected[axis] += 0.5 * ((e[a] * h[b].conj() - e[b] * h[a].conj()) * weight).re;
                }
            }
        }
        for (actual, expected) in actual.mean_poynting_w_m2.into_iter().zip(expected) {
            close(actual, expected, expected.abs().max(1.0) * 2e-14);
        }
        close(
            actual.normal_power_w,
            expected[2] * aperture.area_m2(),
            2e-14,
        );
        // A representable integrated power cannot hide an overflowing mean.
        let large_e = [[C64::ZERO, C64::ONE, C64::ZERO]; 4];
        assert!(prepared.measure(&large_e, &h, &mut || false).is_err());
    }

    #[test]
    fn g0_scaled_magnetic_range_failures_do_not_hide_nonzero_components() {
        let grid = AngularSpectrumGrid {
            width: 1,
            height: 1,
            spacing_m: [1.0; 2],
        };
        let prepared = PreparedAperture::prepare(
            grid,
            AngularSpectrumAperture::try_new([0.0; 2], [1.0; 2]).unwrap(),
            0,
            limits(),
            |_| 0.0,
            |_| 0.0,
            &mut || false,
        )
        .unwrap();
        let e = [[C64::ONE, C64::ZERO, C64::ZERO]];
        for bits in [1, 17] {
            // The common factor is 16: one source vanishes, while the other
            // stays nonzero but loses a bit. Both must be explicitly refused.
            let tiny = f64::from_bits(bits);
            let h = [[C64::from_re(1e308), C64::from_re(tiny), C64::ZERO]];
            assert!(prepared.measure(&e, &h, &mut || false).is_err());
        }
        // Source reduction is now exact, but halving the scaled DD total
        // would erase a nonzero physical mean. Refuse instead of returning zero.
        let tiny = f64::from_bits(16);
        assert!(tiny * 0.5 > 0.0);
        let h = [[C64::from_re(1e308), C64::from_re(tiny), C64::ZERO]];
        assert!(prepared.measure(&e, &h, &mut || false).is_err());
    }

    #[test]
    fn g0_g4_area_work_storage_and_cancellation_refusals_preserve_inputs() {
        for (origin, extent) in [
            ([f64::NAN, 0.0], [1.0; 2]),
            ([0.0; 2], [0.0, 1.0]),
            ([0.0; 2], [-1.0, 1.0]),
            ([0.0; 2], [1e300, 1e300]),
            ([0.0; 2], [1e-300, 1e-300]),
        ] {
            assert!(AngularSpectrumAperture::try_new(origin, extent).is_err());
        }
        let grid = AngularSpectrumGrid {
            width: 2,
            height: 2,
            spacing_m: [1.0; 2],
        };
        let aperture = AngularSpectrumAperture::try_new([0.0; 2], [1.0; 2]).unwrap();
        let exact = AngularSpectrumApertureLimits {
            maximum_mode_pairs: 16,
            maximum_workspace_bytes: 123
                + 8 * core::mem::size_of::<[Dd; 2]>()
                + 4 * core::mem::size_of::<[[Dd; 2]; 3]>(),
        };
        let prepare = |limits, cancelled: &mut dyn FnMut() -> bool| {
            PreparedAperture::prepare(
                grid,
                aperture,
                123,
                limits,
                |i| i as f64,
                |i| -(i as f64),
                &mut || cancelled(),
            )
        };
        assert!(
            prepare(
                AngularSpectrumApertureLimits {
                    maximum_mode_pairs: 15,
                    ..exact
                },
                &mut || false
            )
            .is_err()
        );
        assert!(
            prepare(
                AngularSpectrumApertureLimits {
                    maximum_workspace_bytes: exact.maximum_workspace_bytes - 1,
                    ..exact
                },
                &mut || false
            )
            .is_err()
        );
        assert!(matches!(
            prepare(exact, &mut || true),
            Err(AngularSpectrumError::Cancelled)
        ));
        let prepared = prepare(exact, &mut || false).unwrap();
        let e = [[C64::ONE; 3]; 4];
        let h = [[C64::new(0.1, 0.2), C64::new(0.3, -0.4), C64::new(-0.5, 0.6)]; 4];
        let before = prepared.measure(&e, &h, &mut || false).unwrap();
        let mut polls = 0;
        assert!(matches!(
            prepared.measure(&e, &h, &mut || {
                polls += 1;
                polls == 4
            }),
            Err(AngularSpectrumError::Cancelled)
        ));
        assert_eq!(prepared.measure(&e, &h, &mut || false).unwrap(), before);
        let mut bad_e = e;
        bad_e[0][0] = C64::new(f64::NAN, 0.0);
        assert!(prepared.measure(&bad_e, &h, &mut || false).is_err());
        assert!(prepared.measure(&e[..3], &h, &mut || false).is_err());
    }
}
