use super::*;
use fs_fft::harmonic_window::HarmonicWindowProfile;
use fs_material::field_storage::dielectric::DebyePole;

const DEPTH: [AngularSpectrumHeatDepthSample; 3] = [
    AngularSpectrumHeatDepthSample {
        depth_m: 0.0,
        weight_m: 0.1,
    },
    AngularSpectrumHeatDepthSample {
        depth_m: 0.3,
        weight_m: 0.4,
    },
    AngularSpectrumHeatDepthSample {
        depth_m: 0.6,
        weight_m: 0.1,
    },
];
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "actual={a:e} expected={b:e} tolerance={tolerance:e}"
    );
}
fn grid(w: usize, h: usize) -> AngularSpectrumGrid {
    AngularSpectrumGrid {
        width: w,
        height: h,
        spacing_m: [1.0; 2],
    }
}
fn volume(aperture: AngularSpectrumAperture) -> AngularSpectrumHeatVolumeQuadrature<'static> {
    AngularSpectrumHeatVolumeQuadrature {
        depth_span_m: [0.0, 0.6],
        depth_samples: &DEPTH,
        transverse: AngularSpectrumHeatAreaQuadrature::ContinuousRectangle(aperture),
    }
}
fn spectrum(terms: Vec<AngularHarmonicSpectrumTerm>) -> PreparedHarmonicProductSpectrum {
    let n = terms.len();
    PreparedHarmonicProductSpectrum::from_angular_terms(
        terms,
        HarmonicSpectrumTermLimits {
            maximum_terms: n,
            maximum_workspace_bytes: n * core::mem::size_of::<HarmonicProductTone>(),
        },
    )
    .unwrap()
}
fn dd(values: [C64; 3]) -> [[Dd; 2]; 6] {
    core::array::from_fn(|i| {
        if i < 3 {
            [Dd::from_f64(values[i].re), Dd::from_f64(values[i].im)]
        } else {
            [Dd::ZERO; 2]
        }
    })
}
fn pole() -> DebyePole {
    DebyePole {
        permittivity_increment_f_m: 0.7,
        relaxation_time_s: 0.4,
        capacitance_f: 1.0,
        resistance_ohm: 0.4,
    }
}
fn field(bank: usize, frequency: usize, mode: usize, depth: f64) -> [C64; 3] {
    let (s, c) = (0.3 * mode as f64 + (0.4 + frequency as f64) * depth).sin_cos();
    core::array::from_fn(|axis| {
        let a = C64::new(
            0.2 * (axis + 1) as f64 / (mode + 1) as f64,
            0.1 * (mode as f64 - 1.2),
        );
        let b = C64::new(0.03 * (frequency + 1) as f64, -0.02 * (axis + 1) as f64);
        (a * C64::new(c, -s) + b * C64::new(c, s)).scale(1.0 / (bank + 1) as f64)
    })
}

fn prepared(aperture: AngularSpectrumAperture) -> (PreparedHarmonicProductSpectrum, usize, usize) {
    let omega = [2.0, 2.3];
    let kx = [0.0, -0.4];
    let ky = [0.0, -0.3];
    let g = grid(2, 2);
    let quadrature = volume(aperture);
    let admitted = admit_volume(g, 0.6, quadrature, &mut || false).unwrap();
    let rectangle =
        PreparedRectangle::prepare(g, aperture, |i| kx[i], |i| ky[i], &mut || false).unwrap();
    let poles = omega.map(|w| [PoleTransfer::new(pole(), w).unwrap()]);
    let mut terms = Vec::new();
    let mut calls = 0;
    for bank in 0..3 {
        let proper = bank == 2;
        rectangle
            .contract_bank(
                &omega,
                proper,
                quadrature,
                admitted,
                |i, j, c| reservoir_pair(0.3, &poles[i], &poles[j], 0, 2, proper, c),
                |f, z, _| {
                    calls += 1;
                    Ok((0..4).map(|m| dd(field(bank, f, m, z))).collect())
                },
                &mut terms,
                &mut || false,
            )
            .unwrap();
    }
    let count = terms.len();
    (spectrum(terms), calls, count)
}

// Independent spatial Simpson quadrature of the actual positive-time E and
// Debye current. The proper factor is averaged over {1,-1,i,-i}; its one
// coefficient remains shared across every frequency and spatial mode.
fn reference(aperture: AngularSpectrumAperture, time: f64) -> f64 {
    let steps = 16;
    let omega = [2.0, 2.3];
    let kx = [0.0, -0.4];
    let ky = [0.0, -0.3];
    let weight = |i| {
        if i == 0 || i == steps {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        }
    };
    let mut total = 0.0;
    for ix in 0..=steps {
        for iy in 0..=steps {
            let x = aperture.origin_m()[0] + aperture.extent_m()[0] * ix as f64 / steps as f64;
            let y = aperture.origin_m()[1] + aperture.extent_m()[1] * iy as f64 / steps as f64;
            for depth in DEPTH {
                for bank in 0..3 {
                    let phases = if bank == 2 {
                        vec![C64::ONE, -C64::ONE, C64::new(0.0, 1.0), C64::new(0.0, -1.0)]
                    } else {
                        vec![C64::ONE]
                    };
                    for z in &phases {
                        let mut electric = [0.0; 3];
                        let mut current = [0.0; 3];
                        for f in 0..2 {
                            for m in 0..4 {
                                let (s, c) =
                                    (omega[f] * time - kx[m % 2] * x - ky[m / 2] * y).sin_cos();
                                for axis in 0..3 {
                                    let a = field(bank, f, m, depth.depth_m)[axis] * *z;
                                    electric[axis] += (a * C64::new(c, s)).re;
                                    let p = a.scale(0.7) / C64::new(1.0, omega[f] * 0.4);
                                    current[axis] +=
                                        (C64::new(-p.im, p.re).scale(omega[f]) * C64::new(c, s)).re;
                                }
                            }
                        }
                        let heat = 0.3 * electric.iter().map(|v| v * v).sum::<f64>()
                            + (0.4 / 0.7) * current.iter().map(|v| v * v).sum::<f64>();
                        total +=
                            weight(ix) * weight(iy) * depth.weight_m * heat / phases.len() as f64;
                    }
                }
            }
        }
    }
    total * aperture.area_m2() / (3.0 * steps as f64).powi(2)
}

#[test]
fn rectangle_reservoir_heat_matches_independent_real_field_current_quadrature() {
    let area = AngularSpectrumAperture::try_new([0.21, -0.37], [0.7, 0.9]).unwrap();
    let (heat, calls, tones) = prepared(area);
    assert_eq!(calls, 3 * 2 * DEPTH.len());
    assert_eq!(tones, 15);
    for t in [-0.2, 0.17, 1.3] {
        close(heat.evaluate(t).unwrap(), reference(area, t), 2e-8);
    }
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(0.11, 0.37, profile).unwrap();
        let n = 64;
        let mut expected = 0.0;
        for i in 0..=n {
            let u = i as f64 / n as f64;
            let w = if i == 0 || i == n {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            let w = w * if profile == HarmonicWindowProfile::Triangular {
                2.0 * (1.0 - (2.0 * u - 1.0).abs())
            } else {
                1.0
            };
            expected += w * reference(area, window.start_s() + u * window.duration_s());
        }
        close(
            heat.mean(window).unwrap(),
            expected / (3.0 * n as f64),
            5e-8,
        );
    }
}

#[test]
fn continuous_rectangle_partitions_add_without_relocating_fields() {
    let whole = AngularSpectrumAperture::try_new([-0.4, 0.2], [1.3, 0.7]).unwrap();
    let left = AngularSpectrumAperture::try_new([-0.4, 0.2], [0.3, 0.7]).unwrap();
    let right = AngularSpectrumAperture::try_new([-0.1, 0.2], [1.0, 0.7]).unwrap();
    let (a, _, _) = prepared(whole);
    let (b, _, _) = prepared(left);
    let (c, _, _) = prepared(right);
    for t in [-0.7, 0.2, 1.1] {
        close(
            a.evaluate(t).unwrap(),
            b.evaluate(t).unwrap() + c.evaluate(t).unwrap(),
            2e-14,
        );
    }
}

#[test]
fn whole_rectangle_matches_continuous_periodic_cell_without_nyquist_aliasing() {
    let g = AngularSpectrumGrid {
        width: 4,
        height: 1,
        spacing_m: [core::f64::consts::FRAC_PI_2, 1.0],
    };
    let area =
        AngularSpectrumAperture::try_new([0.23, 0.0], [core::f64::consts::TAU, 1.0]).unwrap();
    let rect_volume = volume(area);
    let full_volume = AngularSpectrumHeatVolumeQuadrature {
        transverse: AngularSpectrumHeatAreaQuadrature::ContinuousPeriodicCell,
        ..rect_volume
    };
    let rectangle =
        PreparedRectangle::prepare(g, area, |i| [0.0, 1.0, -2.0, -1.0][i], |_| 0.0, &mut || {
            false
        })
        .unwrap();
    let field = |m| {
        if m == 2 {
            [C64::ZERO, C64::new(0.0, 1.0), C64::ZERO]
        } else {
            [C64::ZERO; 3]
        }
    };
    let mut a = Vec::new();
    let mut b = Vec::new();
    rectangle
        .contract_bank(
            &[2.0],
            false,
            rect_volume,
            admit_volume(g, 0.6, rect_volume, &mut || false).unwrap(),
            |_, _, c| reservoir_pair(1.0, &[], &[], 0, 1, false, c),
            |_, _, _| Ok((0..4).map(|m| dd(field(m))).collect()),
            &mut a,
            &mut || false,
        )
        .unwrap();
    super::super::contract_bank(
        &[2.0],
        false,
        full_volume,
        admit_volume(g, 0.6, full_volume, &mut || false).unwrap(),
        |_, _, c| reservoir_pair(1.0, &[], &[], 0, 1, false, c),
        |_, _, _| {
            Ok((0..4)
                .map(|m| {
                    let e = field(m);
                    [e[0], e[1], e[2], C64::ZERO, C64::ZERO, C64::ZERO]
                })
                .collect())
        },
        &mut b,
        &mut || false,
    )
    .unwrap();
    let (a, b) = (spectrum(a), spectrum(b));
    for t in [0.0, 0.31, 1.2] {
        close(a.evaluate(t).unwrap(), b.evaluate(t).unwrap(), 3e-15);
        close(a.evaluate(t).unwrap(), 0.3 * core::f64::consts::TAU, 3e-15);
    }
}

#[test]
fn translated_narrow_dark_patch_retains_positive_variation_and_grouped_residual() {
    for (width, residual) in [(1e-20, 0.0), (1e-24, 1e-20)] {
        let g = grid(2, 2);
        let area = AngularSpectrumAperture::try_new([1.0, 1.0], [width, width]).unwrap();
        let q = volume(area);
        let volume = admit_volume(g, 0.6, q, &mut || false).unwrap();
        let rectangle =
            PreparedRectangle::prepare(g, area, |i| [0.0, 2.0][i], |i| [0.0, -2.0][i], &mut || {
                false
            })
            .unwrap();
        for proper in [false, true] {
            let mut terms = Vec::new();
            rectangle
                .contract_bank(
                    &[2.0],
                    proper,
                    q,
                    volume,
                    |_, _, c| reservoir_pair(1.0, &[], &[], 0, 1, proper, c),
                    |_, _, _| {
                        let mut f = vec![[[Dd::ZERO; 2]; 6]; 4];
                        f[0][0][0] = Dd::from_pair(1.0, residual);
                        f[3][0][0] = -Dd::ONE;
                        Ok(f)
                    },
                    &mut terms,
                    &mut || false,
                )
                .unwrap();
            let actual = spectrum(terms).evaluate(0.31).unwrap();
            // E=(1+eps)-exp(-2i(x-y)); <|E|²>=eps²+(1+eps)*2L²/3+O(L⁴).
            // A real bank retains the carrier phase and its sum-frequency
            // term: leading variance is multiplied by sin²(omega*t).
            let variance = 2.0 * width * width / 3.0;
            let mean_square = if proper {
                0.5 * (residual * residual + (1.0 + residual) * variance)
            } else {
                let (s, c) = (2.0_f64 * 0.31).sin_cos();
                residual * residual * c * c + variance * s * s
            };
            let expected = 0.6 * width * width * mean_square;
            assert!(actual > 0.0);
            close(actual / expected, 1.0, 3e-13);
        }
    }
}

#[test]
fn rectangle_preserves_scaled_heat_and_refuses_cancelled_or_mismatched_fields() {
    let g = grid(1, 1);
    let area = AngularSpectrumAperture::try_new([0.0; 2], [1e-100, 1.0]).unwrap();
    let q = volume(area);
    let admitted = admit_volume(g, 0.6, q, &mut || false).unwrap();
    let rectangle = PreparedRectangle::prepare(g, area, |_| 0.0, |_| 0.0, &mut || false).unwrap();
    let mut terms = Vec::new();
    rectangle
        .contract_bank(
            &[2.0],
            true,
            q,
            admitted,
            |_, _, c| reservoir_pair(1e-200, &[], &[], 0, 1, true, c),
            |_, _, _| Ok(vec![dd([C64::from_re(1e200), C64::ZERO, C64::ZERO])]),
            &mut terms,
            &mut || false,
        )
        .unwrap();
    close(spectrum(terms).evaluate(0.0).unwrap() / 3e99, 1.0, 3e-14);
    assert!(matches!(
        PreparedRectangle::prepare(g, area, |_| 0.0, |_| 0.0, &mut || true),
        Err(AngularSpectrumError::Cancelled)
    ));
    assert!(PreparedRectangle::workspace_bytes(grid(usize::MAX, 2)).is_err());
    assert!(rectangle.electric(Vec::new(), &mut || false).is_err());
    let mut calls = 0;
    assert!(matches!(
        rectangle.contract_bank(
            &[2.0],
            false,
            q,
            admitted,
            |_, _, c| reservoir_pair(1.0, &[], &[], 0, 1, false, c),
            |_, _, _| Ok(vec![dd([C64::ONE, C64::ZERO, C64::ZERO])]),
            &mut Vec::new(),
            &mut || {
                calls += 1;
                calls == 7
            }
        ),
        Err(AngularSpectrumError::Cancelled)
    ));
}
