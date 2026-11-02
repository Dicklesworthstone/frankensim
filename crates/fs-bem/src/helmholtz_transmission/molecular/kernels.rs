//! Analytic continuation of the original factored transmission differences.
use super::*;

pub(super) fn difference(
    outside: DampedWaveNumber,
    inside: DampedWaveNumber,
    x: [f64; 3],
    nx: [f64; 3],
    y: [f64; 3],
    ny: [f64; 3],
) -> [C64; 4] {
    if outside.k == inside.k {
        return [C64::ZERO; 4];
    }
    let delta = sub(x, y);
    let r = norm(delta);
    let h = (inside.k - outside.k).scale(0.5 * r);
    let mid = outside.k.scale(r) + h;
    let h2 = h * h;
    let imaginary = C64::new(0.0, 1.0);
    // Weak contrasts retain the O(h^3) moment before cancellation. At large
    // contrast the endpoint polynomials multiply only decaying exponentials:
    // no hyperbolic overflow or cancellation between large midpoint squares.
    let (ds, f, g) = if h.abs() < 0.01 {
        let e = expik(mid.re).scale(det::exp(-mid.im));
        let sine = h
            * (C64::ONE
                + h2 * (C64::from_re(-1.0 / 6.0)
                    + h2 * (C64::from_re(1.0 / 120.0)
                        + h2 * (C64::from_re(-1.0 / 5040.0) + h2.scale(1.0 / 362880.0)))));
        let moment = h
            * h2
            * (C64::from_re(1.0 / 3.0)
                + h2 * (C64::from_re(-1.0 / 30.0)
                    + h2 * (C64::from_re(1.0 / 840.0) - h2.scale(1.0 / 45360.0))));
        let (sine, moment) = (e * sine, e * moment);
        (
            (imaginary * sine).scale(2.0),
            (mid * sine + imaginary * moment).scale(2.0),
            ((mid * moment).scale(-2.0)
                + imaginary * ((mid * mid + h2) * sine - moment.scale(2.0)))
            .scale(2.0),
        )
    } else {
        let eo = outside.outgoing(r);
        let ei = inside.outgoing(r);
        let (zo, zi) = (outside.k.scale(r), inside.k.scale(r));
        (
            ei - eo,
            ei * (C64::ONE - imaginary * zi) - eo * (C64::ONE - imaginary * zo),
            ei * (zi * zi + (imaginary * zi).scale(2.0) - C64::from_re(2.0))
                - eo * (zo * zo + (imaginary * zo).scale(2.0) - C64::from_re(2.0)),
        )
    };
    let xn = dot(nx, delta) / r;
    let yn = dot(ny, delta) / r;
    let scale = 1.0 / (4.0 * core::f64::consts::PI * r);
    [
        ds.scale(scale),
        f.scale(scale * yn / r),
        f.scale(-scale * xn / r),
        (f.scale(dot(nx, ny) - xn * yn) + g.scale(xn * yn)).scale(scale / (r * r)),
    ]
}

pub(super) fn self_contrast(
    outside: DampedWaveNumber,
    inside: DampedWaveNumber,
    area: f64,
) -> (C64, C64, C64) {
    let a = (area / core::f64::consts::PI).sqrt();
    let (zo, zi, dz) = (
        outside.k.scale(a),
        inside.k.scale(a),
        (inside.k - outside.k).scale(a),
    );
    let (mut factorial, mut po, mut pi, mut difference) = (1.0, C64::ONE, C64::ONE, C64::ZERO);
    let (mut phase, mut single, mut ds, mut dn) =
        (C64::ONE, C64::from_re(0.5 * a), C64::ZERO, C64::ZERO);
    for n in 1..=24 {
        difference = pi * dz + zo * difference;
        po = po * zo;
        pi = pi * zi;
        factorial *= f64::from(n);
        phase = phase * C64::new(0.0, 1.0);
        single = single + (phase * po).scale(0.5 * a / (factorial * f64::from(n + 1)));
        ds = ds + (phase * difference).scale(0.5 * a / (factorial * f64::from(n + 1)));
        if n >= 2 {
            dn = dn - (phase * difference).scale(1.0 / (2.0 * a * factorial));
        }
    }
    (single, ds, dn)
}
