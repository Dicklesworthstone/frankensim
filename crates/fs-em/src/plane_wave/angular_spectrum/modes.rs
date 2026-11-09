//! Maxwell reconstruction for the parent's exp(+i omega t - i k.r) convention.
use fs_math::{c64::C64, dd::Dd, det};

pub(super) fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}

/// Positive incident cosine from the original transverse wave numbers.
/// A rounded hypot can turn a propagating mode into exact grazing before its
/// small longitudinal wave number is formed. Preserve the squared-distance
/// residual instead. Binary scaling keeps the products finite at every wave
/// number scale, and subtracting the larger transverse component first keeps
/// the cancellation in the compensated products. No light-cone tolerance or
/// clipping admits grazing/evanescent input.
pub(in crate::plane_wave) fn incident_cosine(k: f64, kx: f64, ky: f64) -> Option<f64> {
    let (a, b) = if kx.abs() >= ky.abs() {
        (kx.abs(), ky.abs())
    } else {
        (ky.abs(), kx.abs())
    };
    if !k.is_finite() || k <= 0.0 || !a.is_finite() || !b.is_finite() || a >= k {
        return None;
    }
    // The exponent field is an exact power of two. Subnormal k uses the
    // smallest normal divisor, leaving k² normal even for the smallest input.
    let scale = f64::from_bits(k.to_bits() & 0x7ff0_0000_0000_0000).max(f64::MIN_POSITIVE);
    let scaled_k = k / scale;
    let k = Dd::from_f64(scaled_k);
    let a = Dd::from_f64(a / scale);
    let b = Dd::from_f64(b / scale);
    let square = (k * k - a * a - b * b).to_f64();
    if square <= 0.0 {
        return None;
    }
    let cosine = det::sqrt(square) / scaled_k;
    (cosine.is_finite() && cosine > 0.0 && cosine <= 1.0).then_some(cosine)
}

/// Fixed tangential TM axis p=(kx,ky)/kt, and TE axis s=(-py,px).
/// At normal incidence the otherwise arbitrary basis is p=x, s=y.
pub(super) fn basis(kx: f64, ky: f64) -> [f64; 2] {
    let kt = det::hypot(kx, ky);
    if kt == 0.0 {
        [1.0, 0.0]
    } else {
        [kx / kt, ky / kt]
    }
}

pub(super) fn polarized(e: [C64; 2], p: [f64; 2]) -> [C64; 2] {
    [
        e[1].scale(p[0]) - e[0].scale(p[1]),
        e[0].scale(p[0]) + e[1].scale(p[1]),
    ]
}

pub(super) fn cartesian(te_tm: [C64; 2], p: [f64; 2]) -> [C64; 2] {
    [
        te_tm[1].scale(p[0]) - te_tm[0].scale(p[1]),
        te_tm[1].scale(p[1]) + te_tm[0].scale(p[0]),
    ]
}

/// The decaying branch is used for transmitted fields. Negative distance is
/// permitted only for lossless incident waves; growing evanescent waves refuse.
pub(in crate::plane_wave) fn propagation(kz: C64, distance: f64) -> Option<C64> {
    let (phase, decay) = propagation_phase(kz, distance)?;
    Some(phase.scale(det::exp(decay)))
}

fn propagation_phase(kz: C64, distance: f64) -> Option<(C64, f64)> {
    let decay = kz.im * distance;
    let phase = kz.re * distance;
    if decay.is_nan() || decay > 0.0 || !phase.is_finite() {
        return None;
    }
    // Keep the product residual until argument reduction. Optical paths can
    // accumulate millions of cycles even for metre-scale authored geometry.
    // This retains the supplied binary64 operands, not their physical uncertainty.
    let (sine, cosine) = det::sin_cos_dd(Dd::from_f64(kz.re) * Dd::from_f64(distance));
    let phase = C64::new(cosine.to_f64(), -sine.to_f64());
    finite(phase).then_some((phase, decay))
}

/// One phase/decay preparation shared by all six Maxwell field components.
/// A unit-amplitude transfer may underflow before its product with a physical
/// field does. Split only that late attenuation, as in fs_math::decay, while
/// retaining the original multiplication order on the ordinary path.
pub(super) struct FieldPropagation {
    transfer: C64,
    final_attenuation: f64,
}

impl FieldPropagation {
    pub(super) fn new(kz: C64, distance: f64) -> Option<Self> {
        let (phase, decay) = propagation_phase(kz, distance)?;
        let amplitude = det::exp(decay);
        let (amplitude, final_attenuation) = if amplitude < f64::MIN_POSITIVE {
            let half = det::exp(0.5 * decay);
            (half, half)
        } else {
            (amplitude, 1.0)
        };
        Some(Self {
            transfer: phase.scale(amplitude),
            final_attenuation,
        })
    }

    pub(super) fn apply(&self, value: C64) -> C64 {
        (value * self.transfer).scale(self.final_attenuation)
    }
}

/// Transversality determines Ez; Faraday determines the entire H. The caller
/// admits nonzero kz and positive omega*mu. No paraxial or impedance shortcut.
pub(super) fn fields(
    kx: f64,
    ky: f64,
    kz: C64,
    omega_mu: f64,
    transverse: [C64; 2],
) -> Option<([C64; 3], [C64; 3])> {
    let ez = -(transverse[0].scale(kx) + transverse[1].scale(ky)) / kz;
    let e = [transverse[0], transverse[1], ez];
    let h = [
        (e[2].scale(ky) - kz * e[1]).scale(1.0 / omega_mu),
        (kz * e[0] - e[2].scale(kx)).scale(1.0 / omega_mu),
        (e[1].scale(kx) - e[0].scale(ky)).scale(1.0 / omega_mu),
    ];
    e.iter().chain(&h).all(|&z| finite(z)).then_some((e, h))
}

pub(super) fn normal_flux(e: [C64; 3], h: [C64; 3]) -> f64 {
    0.5 * (e[0] * h[1].conj() - e[1] * h[0].conj()).re
}

/// Total TE/TM tangential E and power-conjugate H determine all six fields,
/// including a critical internal mode. This avoids dividing by kz entirely.
pub(super) fn fields_from_tangential(
    kx: f64,
    ky: f64,
    omega_epsilon: C64,
    omega_mu: f64,
    electric_te_tm: [C64; 2],
    magnetic_te_tm: [C64; 2],
) -> Option<([C64; 3], [C64; 3])> {
    let p = basis(kx, ky);
    let transverse = cartesian(electric_te_tm, p);
    let kt = det::hypot(kx, ky);
    let e = [
        transverse[0],
        transverse[1],
        -magnetic_te_tm[1].scale(kt) / omega_epsilon,
    ];
    let h = [
        -magnetic_te_tm[0].scale(p[0]) - magnetic_te_tm[1].scale(p[1]),
        -magnetic_te_tm[0].scale(p[1]) + magnetic_te_tm[1].scale(p[0]),
        electric_te_tm[0].scale(kt / omega_mu),
    ];
    e.iter().chain(&h).all(|&z| finite(z)).then_some((e, h))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: C64, b: C64) {
        assert!(
            (a - b).abs() < 4e-14 * (1.0 + a.abs() + b.abs()),
            "{a:?} != {b:?}"
        );
    }

    #[test]
    fn light_cone_keeps_adjacent_propagating_directions_and_exact_grazing() {
        // G1/G3: independent 90-digit evaluation of sqrt(1-x*x-y*y), using the
        // exact binary64 operands. The first hypot rounds to one; the second
        // rounded-hypot path gives an axial wave number about 29% too large.
        let cases = [
            (0.8, 0.6_f64.next_down(), 9.424_321_830_774_483e-9),
            (0.6, 0.8_f64.next_down(), 1.154_238_982_858_484e-8),
        ];
        for exponent in [-800, -400, 0, 400, 800] {
            let scale = 2.0_f64.powi(exponent);
            for (x, y, expected) in cases {
                let actual = incident_cosine(scale, x * scale, y * scale).unwrap();
                assert!(
                    (actual / expected - 1.0).abs() < 5e-16,
                    "{actual} != {expected}"
                );
                assert_eq!(incident_cosine(scale, -y * scale, x * scale), Some(actual));
            }
            assert!(incident_cosine(scale, 0.8 * scale, 0.6 * scale).is_none());
            assert!(incident_cosine(scale, scale, 0.0).is_none());
            // An exactly represented Pythagorean triple remains on the cone;
            // division rounding must not introduce a false positive cosine.
            assert!(incident_cosine(5.0 * scale, 3.0 * scale, 4.0 * scale).is_none());
            assert!(incident_cosine(5.0 * scale, 3.0 * scale, (4.0 * scale).next_up()).is_none());
            assert_eq!(incident_cosine(scale, 0.0, 0.0), Some(1.0));
        }
        assert_eq!(incident_cosine(f64::from_bits(1), 0.0, 0.0), Some(1.0));
    }

    #[test]
    fn arbitrary_polarization_satisfies_both_maxwell_curls_and_divergences() {
        let kx = 0.8;
        let ky = -0.6;
        for kz in [C64::from_re(2.0), C64::from_re(-2.0), C64::new(1.2, -0.4)] {
            let omega_mu = 3.0;
            let omega_epsilon = (kz * kz + C64::from_re(kx * kx + ky * ky)).scale(1.0 / omega_mu);
            let (e, h) = fields(
                kx,
                ky,
                kz,
                omega_mu,
                [C64::new(0.7, -0.3), C64::new(-0.2, 1.1)],
            )
            .unwrap();
            close(e[0].scale(kx) + e[1].scale(ky) + kz * e[2], C64::ZERO);
            close(h[0].scale(kx) + h[1].scale(ky) + kz * h[2], C64::ZERO);
            let curl_h = [
                h[2].scale(ky) - kz * h[1],
                kz * h[0] - h[2].scale(kx),
                h[1].scale(kx) - h[0].scale(ky),
            ];
            for c in 0..3 {
                close(curl_h[c], -omega_epsilon * e[c]);
            }
        }
    }

    #[test]
    fn basis_is_invertible_and_flux_uses_te_tm_admittances() {
        let tangential = [C64::new(1.0, -0.2), C64::new(0.4, 0.9)];
        for (kx, ky) in [(0.0, 0.0), (0.6, 0.8), (-0.8, 0.6)] {
            let p = basis(kx, ky);
            let amplitudes = polarized(tangential, p);
            let back = cartesian(amplitudes, p);
            close(back[0], tangential[0]);
            close(back[1], tangential[1]);
            for sign in [-1.0, 1.0] {
                let kz = sign * 2.0;
                let omega_mu = 3.0;
                let (e, h) = fields(kx, ky, C64::from_re(kz), omega_mu, tangential).unwrap();
                let expected = 0.5
                    * (amplitudes[0].norm_sq() * kz / omega_mu
                        + amplitudes[1].norm_sq() * (kx * kx + ky * ky + kz * kz)
                            / (omega_mu * kz));
                assert!((normal_flux(e, h) - expected).abs() < 2e-15);
            }
        }
    }

    #[test]
    fn evanescent_fields_decay_without_false_normal_power() {
        let kz = C64::new(0.0, -2.0);
        let (e, h) = fields(3.0, 0.0, kz, 1.0, [C64::ONE, C64::new(0.0, 0.7)]).unwrap();
        assert!(normal_flux(e, h).abs() < 1e-15);
        close(propagation(kz, 0.3).unwrap(), C64::from_re(det::exp(-0.6)));
        assert_eq!(propagation(kz, 1e300), Some(C64::ZERO));
        assert!(propagation(kz, -0.1).is_none());
        assert!(fields(3.0, 0.0, C64::ZERO, 1.0, [C64::ONE, C64::ZERO]).is_none());
        let forward = propagation(C64::from_re(2.0), 0.3).unwrap();
        close(
            propagation(C64::from_re(2.0), -0.3).unwrap(),
            forward.conj(),
        );
    }

    #[test]
    fn optical_phase_products_preserve_meter_and_kilometer_dark_ports() {
        // Independent 110-digit decimal reduction/Taylor values of exp(-i*k*d),
        // with each authored binary64 operand converted exactly. k=1e7 rad/m
        // corresponds to a 628 nm wavelength. Rounding k*d first loses
        // 3.58e-10 rad over this metre and 7.70e-7 rad over this kilometre.
        for (distance, expected) in [
            (
                1.000_000_000_000_000_2,
                C64::new(-0.907_270_387_115_543_2, -0.420_547_791_176_237_53),
            ),
            (
                1000.000_000_000_000_1,
                C64::new(0.873_120_176_906_475_4, 0.487_505_032_465_107_1),
            ),
        ] {
            let actual = propagation(C64::from_re(1e7), distance).unwrap();
            // A matched opposing reference must not acquire spurious dark-port
            // power from a prematurely rounded propagation phase.
            assert!(
                (actual - expected).norm_sq() < 1e-28,
                "distance={distance}: {actual:?} != {expected:?}"
            );
            close(
                propagation(C64::from_re(1e7), -distance).unwrap(),
                expected.conj(),
            );
            close(
                propagation(C64::new(1e7, -1e-4), distance).unwrap(),
                expected.scale(det::exp(-1e-4 * distance)),
            );
        }
    }

    #[test]
    fn physical_fields_retain_deep_evanescent_and_absorbing_tails() {
        // G1: independent high-precision exp(-x)*1e100. The first raw
        // exponential is subnormal; the second has already rounded to zero.
        for (exposure, surviving) in [
            (720.0, 2.032_230_802_424_293_2e-213),
            (800.0, 3.667_874_584_177_687e-248),
        ] {
            for phase_rate in [0.0, 0.125] {
                let transfer = FieldPropagation::new(C64::new(phase_rate, -1.0), exposure).unwrap();
                let actual = transfer.apply(C64::new(1e100, -3e100));
                let (sine, cosine) = (phase_rate * exposure).sin_cos();
                let expected = C64::new(surviving, -3.0 * surviving) * C64::new(cosine, -sine);
                assert!((actual.re / expected.re - 1.0).abs() < 3e-13);
                assert!((actual.im / expected.im - 1.0).abs() < 3e-13);
            }
        }
        assert_eq!(
            FieldPropagation::new(C64::new(0.0, -1.0), 744.1)
                .unwrap()
                .apply(C64::from_re(0.4))
                .re,
            f64::from_bits(1),
        );
        assert_eq!(
            FieldPropagation::new(C64::new(0.0, -1.0), 1e300)
                .unwrap()
                .apply(C64::new(1e100, -3e100)),
            C64::ZERO,
        );
        assert!(FieldPropagation::new(C64::new(0.0, -1.0), -1.0).is_none());
    }

    #[test]
    fn physical_propagation_keeps_ordinary_bits_and_depth_composition() {
        // G3: the normal-attenuation path retains the original multiplication
        // order. Deep decay composes without erasing a representable field.
        let input = C64::new(0.7, -0.4);
        for kz in [C64::ZERO, C64::from_re(3.0), C64::new(1.2, -0.3)] {
            for distance in [0.0, 0.3, 100.0] {
                assert_eq!(
                    FieldPropagation::new(kz, distance).unwrap().apply(input),
                    input * propagation(kz, distance).unwrap(),
                );
            }
        }
        let kz = C64::new(0.125, -1.0);
        let initial = C64::new(1e100, -3e100);
        let first = FieldPropagation::new(kz, 400.0).unwrap().apply(initial);
        let split = FieldPropagation::new(kz, 400.0).unwrap().apply(first);
        let whole = FieldPropagation::new(kz, 800.0).unwrap().apply(initial);
        assert!((split.re / whole.re - 1.0).abs() < 4e-14);
        assert!((split.im / whole.im - 1.0).abs() < 4e-14);
    }

    #[test]
    fn tangential_state_reconstruction_agrees_with_forward_and_backward_plane_waves() {
        for kz in [
            C64::from_re(2.0),
            C64::from_re(-2.0),
            C64::new(0.7, -0.3),
            C64::new(0.0, -1.2),
        ] {
            let (kx, ky, omega_mu) = (0.6, -0.8, 1.4);
            let omega_epsilon = (kz * kz + C64::ONE).scale(1.0 / omega_mu);
            let et = [C64::new(0.4, -0.2), C64::new(0.1, 0.8)];
            let p = basis(kx, ky);
            let a = polarized(et, p);
            let conjugate_h = [a[0] * kz.scale(1.0 / omega_mu), a[1] * omega_epsilon / kz];
            let (expected_e, expected_h) = fields(kx, ky, kz, omega_mu, et).unwrap();
            let (e, h) =
                fields_from_tangential(kx, ky, omega_epsilon, omega_mu, a, conjugate_h).unwrap();
            for c in 0..3 {
                close(e[c], expected_e[c]);
                close(h[c], expected_h[c]);
            }
        }
    }

    #[test]
    fn critical_and_absorbing_internal_states_satisfy_both_full_maxwell_curls() {
        let (kx, ky, omega_mu) = (0.6, -0.8, 1.4);
        for omega_epsilon in [C64::from_re(1.0 / omega_mu), C64::new(2.0, -0.7)] {
            let et = [C64::new(0.4, -0.2), C64::new(0.1, 0.8)];
            let ht = [C64::new(0.2, 0.3), C64::new(-0.5, 0.6)];
            let i = C64::new(0.0, 1.0);
            let de = [
                -i * ht[0].scale(omega_mu),
                -i * (C64::from_re(omega_mu) - C64::ONE / omega_epsilon) * ht[1],
            ];
            let dh = [
                -i * (omega_epsilon - C64::from_re(1.0 / omega_mu)) * et[0],
                -i * omega_epsilon * et[1],
            ];
            let (e, h) = fields_from_tangential(kx, ky, omega_epsilon, omega_mu, et, ht).unwrap();
            let (dz_e, dz_h) =
                fields_from_tangential(kx, ky, omega_epsilon, omega_mu, de, dh).unwrap();
            let curl = |f: [C64; 3], dz: [C64; 3]| {
                [
                    -i * f[2].scale(ky) - dz[1],
                    dz[0] + i * f[2].scale(kx),
                    -i * f[1].scale(kx) + i * f[0].scale(ky),
                ]
            };
            let curl_e = curl(e, dz_e);
            let curl_h = curl(h, dz_h);
            for c in 0..3 {
                close(curl_e[c], -i * h[c].scale(omega_mu));
                close(curl_h[c], i * omega_epsilon * e[c]);
            }
            close(e[0].scale(kx) + e[1].scale(ky), -i * dz_e[2]);
            close(h[0].scale(kx) + h[1].scale(ky), -i * dz_h[2]);
        }
    }
}
