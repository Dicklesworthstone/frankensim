//! Maxwell reconstruction for the parent's exp(+i omega t - i k.r) convention.
use fs_math::{c64::C64, det};

pub(super) fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
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
pub(super) fn propagation(kz: C64, distance: f64) -> Option<C64> {
    let decay = kz.im * distance;
    let phase = kz.re * distance;
    if decay.is_nan() || decay > 0.0 || !phase.is_finite() {
        return None;
    }
    let amplitude = det::exp(decay);
    let value = C64::new(det::cos(phase), -det::sin(phase)).scale(amplitude);
    finite(value).then_some(value)
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
}
