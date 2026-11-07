//! Spatial accumulator regressions with isolated, already-solved wave payloads.
use super::*;
use fs_math::det;

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}

fn wave(k: C64, e: C64, h: C64) -> SlabWave {
    SlabWave {
        wave_kind: if k.im == 0.0 {
            UniaxialWaveKind::Propagating
        } else {
            UniaxialWaveKind::Evanescent
        },
        normal_wavenumber_rad_m: k,
        electric_v_m: [e, C64::ZERO, C64::ZERO],
        magnetic_a_m: [C64::ZERO, h, C64::ZERO],
    }
}

fn isolated(entrance: f64, thickness: f64) -> UniaxialStackResponse {
    // Start from the actual admitted Maxwell solve. Replace only its field
    // payloads to isolate spatial summation from the independently tested solve.
    let data = super::tests::fixture(1.0, 1.0, &[(1.0, 1.0, 1.0, [1.0, 0.0, 0.0])]);
    let mut result = solve(
        &data,
        Incidence::TangentialIndex(0.0),
        [C64::ONE, C64::ZERO],
        [1.0, 0.0],
        &mut || false,
    )
    .unwrap();
    let zero = wave(C64::ONE, C64::ZERO, C64::ZERO);
    result.layers[0] = UniaxialLayerField {
        entrance_z_m: entrance,
        exit_z_m: entrance + thickness,
        thickness_m: thickness,
        forward: [zero; 2],
        backward: [zero; 2],
    };
    result.thickness_m = entrance + thickness;
    result
}

fn phasor(phase: Dd) -> C64 {
    let (s, c) = det::sin_cos_dd(phase);
    C64::new(c.to_f64(), -s.to_f64())
}

#[test]
fn authored_exit_and_compensated_local_depth_preserve_layer_phase() {
    let mut result = isolated(1e16, 3.0);
    result.layers[0].forward[0] = wave(C64::from_re(0.4), C64::ONE, C64::ONE);
    result.layers[0].backward[0] = wave(C64::from_re(-0.4), C64::from_re(0.2), C64::from_re(-0.2));
    let layer = &result.layers[0];
    assert_ne!(layer.exit_z_m - layer.entrance_z_m, layer.thickness_m);
    let at = result
        .field_at(StackRegion::Layer(0), [0.0, 0.0, layer.exit_z_m])
        .unwrap();
    let forward = phasor(Dd::from_f64(0.4) * Dd::from_f64(3.0));
    close(
        (at.electric_v_m[0] - (forward + C64::from_re(0.2))).abs(),
        0.0,
        3e-15,
    );
    close(
        (at.magnetic_a_m[1] - (forward - C64::from_re(0.2))).abs(),
        0.0,
        3e-15,
    );
    close(at.poynting_w_m2[2], 0.48, 3e-15);
    assert_eq!(at, layer.face(true).unwrap());

    // The small entrance coordinate disappears in a rounded local f64 depth,
    // but its phase at large k remains an observable one-radian difference.
    let mut result = isolated(1e-20, 1.0);
    result.layers[0].forward[0] = wave(C64::from_re(1e20), C64::ONE, C64::ONE);
    let at = result
        .field_at(StackRegion::Layer(0), [0.0, 0.0, 0.5])
        .unwrap();
    let expected =
        phasor(Dd::from_f64(1e20) * Dd::from_f64(0.5) - Dd::from_f64(1e20) * Dd::from_f64(1e-20));
    close((at.electric_v_m[0] - expected).abs(), 0.0, 3e-15);
    let rounded = phasor(Dd::from_f64(1e20 * (0.5 - 1e-20)));
    assert!((at.electric_v_m[0] - rounded).abs() > 0.5);
}

#[test]
fn internal_near_degenerate_waves_keep_weak_interference_before_power() {
    let mut result = isolated(0.0, 1.0);
    let k = [1.0, f64::from_bits(1.0_f64.to_bits() + 1)];
    result.layers[0].forward = [
        wave(C64::from_re(k[0]), C64::ONE, C64::ONE),
        wave(C64::from_re(k[1]), -C64::ONE, -C64::ONE),
    ];
    let depth = 1e-4;
    let at = result
        .field_at(StackRegion::Layer(0), [0.0, 0.0, depth])
        .unwrap();
    let amplitude = (k[1] - k[0]) * depth;
    close(at.electric_v_m[0].abs() / amplitude, 1.0, 3e-14);
    close(at.magnetic_a_m[1].abs() / amplitude, 1.0, 3e-14);
    close(
        at.poynting_w_m2[2] / (0.5 * amplitude * amplitude),
        1.0,
        3e-14,
    );
}

#[test]
fn forward_and_backward_evanescent_tails_survive_underflowed_raw_fractions() {
    assert_eq!((-800.0_f64).exp(), 0.0);
    for backward in [false, true] {
        let mut result = isolated(10.0, 1000.0);
        let amplitude = 1e150_f64;
        let sign = if backward { 1.0 } else { -1.0 };
        let payload = wave(
            C64::new(0.0, sign),
            C64::from_re(amplitude),
            C64::new(0.0, sign * amplitude),
        );
        let depth = if backward {
            result.layers[0].backward[0] = payload;
            200.0
        } else {
            result.layers[0].forward[0] = payload;
            800.0
        };
        let at = result
            .field_at(StackRegion::Layer(0), [0.0, 0.0, 10.0 + depth])
            .unwrap();
        let expected = (amplitude.ln() - 800.0).exp();
        assert!(expected > 0.0);
        close(at.electric_v_m[0].re / expected, 1.0, 2e-13);
        close(at.magnetic_a_m[1].im / (sign * expected), 1.0, 2e-13);
        assert_eq!(at.poynting_w_m2, [0.0; 3]);
    }
}
