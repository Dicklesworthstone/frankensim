use super::*;

fn source(length: f64) -> Slab {
    Slab::new(
        4.0,
        [C64::new(60.0, 20.0), C64::new(30.0, 10.0)],
        20.0,
        length,
        [20.0, 40.0],
    )
    .unwrap()
}

fn close(a: f64, b: f64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}

#[test]
fn g1_positive_depth_integrals_match_constitutive_quadrature_and_partition() {
    let slab = source(2.0);
    let driven = slab
        .drive([C64::new(0.7, -0.2), C64::new(0.3, 0.1)], 0.7)
        .unwrap();
    let edges = [0.0, 0.03, 0.23, 0.7, 1.49, 2.0];
    let mut sum = [Dd::ZERO; 2];
    for pair in edges.windows(2) {
        let interval = [pair[0], pair[1]];
        let actual = slab
            .interval_heat(driven.inward, driven.parity, interval, 0.7)
            .unwrap();
        let count = 4096;
        let dx = (pair[1] - pair[0]) / f64::from(count);
        let mut quadrature = [0.0; 2];
        for i in 0..=count {
            let depth = if i == count {
                pair[1]
            } else {
                pair[0] + f64::from(i) * dx
            };
            let field = slab.field(driven.inward, driven.parity, depth).unwrap();
            // Independent constitutive work from the original strain field.
            let square = field.strain.abs().powi(2);
            let density = [
                0.5 * slab.omega * slab.moduli[0].im * square,
                (2.0 / 3.0) * slab.omega * slab.moduli[1].im * square,
            ];
            let weight = if i == 0 || i == count {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            for channel in 0..2 {
                quadrature[channel] += weight * density[channel];
            }
        }
        for channel in 0..2 {
            close(actual[channel], quadrature[channel] * dx * 0.7 / 3.0, 3e-11);
            sum[channel] = sum[channel] + Dd::from_f64(actual[channel]);
        }
    }
    for channel in 0..2 {
        close(sum[channel].to_f64(), driven.heat[channel], 3e-13);
    }
    assert_eq!(
        slab.interval_heat(driven.inward, driven.parity, [0.0, 2.0], 0.7)
            .unwrap(),
        driven.heat
    );
}

#[test]
fn g3_depth_partition_keeps_weak_thin_standing_wave_heat() {
    for length in [1e-3, 1e-10, 1e-80] {
        let slab = source(length);
        let driven = slab.drive([C64::ONE, -C64::ONE], 1.0).unwrap();
        let mut sums = [Dd::ZERO; 2];
        for i in 0..8 {
            let a = length * f64::from(i) / 8.0;
            let b = length * f64::from(i + 1) / 8.0;
            let q = slab
                .interval_heat(driven.inward, driven.parity, [a, b], 1.0)
                .unwrap();
            for channel in 0..2 {
                assert!(q[channel] > 0.0);
                sums[channel] = sums[channel] + Dd::from_f64(q[channel]);
            }
        }
        for channel in 0..2 {
            close(sums[channel].to_f64(), driven.heat[channel], 2e-12);
        }
    }
}

#[test]
fn g3_opaque_depth_heat_survives_square_underflow_without_cumulative_subtraction() {
    let mut slab = source(1.0);
    let attenuation = -slab.k.im;
    slab = source(1000.0 / attenuation);
    let driven = slab.drive([C64::ONE, C64::ZERO], 1e300).unwrap();
    assert_eq!(driven.inward[1], C64::ZERO);
    let depths = [400.0 / attenuation, 401.0 / attenuation];
    let heat = slab
        .interval_heat(driven.inward, driven.parity, depths, 1e300)
        .unwrap();
    for channel in 0..2 {
        let expected = det::exp(det::ln(driven.heat[channel]) - 800.0 + det::ln(-det::expm1(-2.0)));
        assert!(heat[channel] > 0.0);
        close(heat[channel], expected, 4e-12);
    }
}

#[test]
fn g0_interval_bounds_zero_measure_and_dark_fields_are_explicit() {
    let slab = source(2.0);
    let drive = slab.drive([C64::ONE, C64::ZERO], 1.0).unwrap();
    assert_eq!(
        slab.interval_heat(drive.inward, drive.parity, [1.0, 1.0], 1.0)
            .unwrap(),
        [0.0; 2]
    );
    for range in [[-1.0, 1.0], [1.0, 0.0], [0.0, 2.1], [0.0, f64::NAN]] {
        assert!(
            slab.interval_heat(drive.inward, drive.parity, range, 1.0)
                .is_err()
        );
    }
    assert!(
        slab.interval_heat(drive.inward, drive.parity, [0.0, 1.0], 0.0)
            .is_err()
    );
    let dark = slab.drive([C64::ZERO; 2], 1.0).unwrap();
    assert_eq!(
        slab.interval_heat(dark.inward, dark.parity, [0.3, 0.7], 1.0)
            .unwrap(),
        [0.0; 2]
    );
}
