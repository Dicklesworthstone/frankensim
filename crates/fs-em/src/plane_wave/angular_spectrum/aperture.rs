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
/// A query contracts N² ordered mode pairs; the two axis caches need at most
/// 2*N² additional bounded scalar kernel evaluations. Numerical storage includes
/// the original beam workspace and both axis-pair caches. Caller inputs,
/// previously retained outputs and allocator overhead are excluded.
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
        let pairs = n
            .checked_mul(n)
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

    pub(super) fn measure(
        &self,
        electric: &[[C64; 3]],
        magnetic: &[[C64; 3]],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AngularSpectrumApertureReading, AngularSpectrumError> {
        let n = self.grid.width * self.grid.height;
        if electric.len() != n || magnetic.len() != n {
            return Err(bad("aperture fields must cover the original modal grid"));
        }
        let mut sum = [Dd::ZERO; 3];
        for (i, e) in electric.iter().enumerate() {
            poll(cancelled)?;
            let row_i = i / self.grid.width;
            let column_i = i % self.grid.width;
            for (j, h) in magnetic.iter().enumerate() {
                if j % 32 == 0 {
                    poll(cancelled)?;
                }
                let x = self.x[column_i * self.grid.width + j % self.grid.width];
                let y = self.y[row_i * self.grid.height + j / self.grid.width];
                let weight = [x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0]];
                for (axis, (a, b)) in [(1, 2), (2, 0), (0, 1)].into_iter().enumerate() {
                    sum[axis] =
                        sum[axis] + product(e[a], h[b], weight) - product(e[b], h[a], weight);
                }
            }
        }
        let half = Dd::from_f64(0.5);
        let mean_poynting_w_m2 = sum.map(|value| (value * half).to_f64());
        let normal_power_w = (sum[2] * half * Dd::from_f64(self.aperture.area_m2)).to_f64();
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

fn product(e: C64, h: C64, weight: [Dd; 2]) -> Dd {
    let [er, ei, hr, hi] = [e.re, e.im, h.re, h.im].map(Dd::from_f64);
    let real = er * hr + ei * hi;
    let imaginary = ei * hr - er * hi;
    real * weight[0] - imaginary * weight[1]
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
            maximum_workspace_bytes: 123 + 8 * core::mem::size_of::<[Dd; 2]>(),
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
