//! Analytic volume integral of the original anisotropic dielectric work.
//!
//! E(z) = Eo exp(-i ko z) + Ee exp(-i ke z). After projection onto the
//! positive parallel/perpendicular loss channels, each scalar intensity has
//! two self terms and 2 Re(conj(Eo) Ee integral(exp(s z))), where
//! s = Im(ko)+Im(ke) + i(Re(ko)-Re(ke)). Integrating this expression retains
//! the mode interference that independent Beer--Lambert powers discard.
//! The spatial power balance follows the Poynting treatment in Passler,
//! Jeannin and Paarmann (2020), <https://arxiv.org/abs/2002.03832>.

use super::{AbsorbingCrystalResponse, Result, Vector};
use fs_math::{dd::Dd, decay, det};

type Complex = [Dd; 2];
const ZERO: Complex = [Dd::ZERO; 2];

pub(super) fn integrated_heat(
    response: &AbsorbingCrystalResponse,
    area: f64,
    start: f64,
    end: f64,
) -> Result<f64> {
    if response.loss_f_m == [0.0; 2] {
        return Ok(0.0);
    }
    let width = Dd::from_f64(end) - Dd::from_f64(start);
    let normals = response
        .transmitted
        .map(|mode| mode.normal_wavenumber_rad_m);
    let alpha = normals.map(|k| -k.im);
    let beat = Dd::from_f64(normals[0].re) - Dd::from_f64(normals[1].re);
    let mut electric = core::array::from_fn::<_, 2, _>(|m| {
        attenuated(
            response.transmitted[m].electric_v_m,
            alpha[m],
            Dd::from_f64(start),
        )
    });
    let field_scale = electric
        .iter()
        .flatten()
        .flatten()
        .fold(0.0_f64, |s, v| s.max(v.to_f64().abs()));
    if field_scale == 0.0 {
        return Ok(0.0);
    }
    let scale = Dd::from_f64(field_scale);
    for component in electric.iter_mut().flatten() {
        *component = component.map(|v| v / scale);
    }
    if electric.iter().all(|mode| mode.iter().any(|v| *v != ZERO)) {
        let phase = beat * Dd::from_f64(start);
        if !phase.hi.is_finite() || !phase.lo.is_finite() {
            return Err("crystal absorption relative phase is unrepresentable");
        }
        let (sine, cosine) = det::sin_cos_dd(phase);
        electric[1] = electric[1].map(|v| multiply(v, [cosine, sine]));
    }
    let axis = response.axis.map(Dd::from_f64);
    let norm = axis.iter().fold(Dd::ZERO, |s, &v| s + v * v);
    let projected = electric.map(|e| {
        let dot = [0, 1].map(|part| (0..3).fold(Dd::ZERO, |s, j| s + e[j][part] * axis[j]));
        core::array::from_fn::<_, 4, _>(|j| {
            if j == 3 {
                dot
            } else {
                [0, 1].map(|part| e[j][part] - dot[part] * axis[j] / norm)
            }
        })
    });
    let exposures = alpha.map(|a| product(a, width));
    let phase = beat * width;
    let short = exposures.iter().all(|&v| v < 0.125)
        && phase.hi.is_finite()
        && phase.to_f64().abs() < 0.125;
    let (integrals, length): ([Dd; 4], f64) = if alpha[0] == alpha[1] && beat == Dd::ZERO {
        let weight = decay_integral(alpha[0], width);
        (
            core::array::from_fn(|j| {
                let sum = [0, 1].map(|part| projected[0][j][part] + projected[1][j][part]);
                dot(sum, sum) * weight
            }),
            1.0,
        )
    } else if short {
        // Sum the field Taylor coefficients BEFORE squaring. At a coherent
        // node this retains O(width^3) heat instead of cancelling O(width)
        // diagonal and interference integrals. Relative carrier phase alone
        // matters, so an arbitrarily short optical period adds no work.
        let rates = [
            [-Dd::from_f64(alpha[0]) * width, Dd::ZERO],
            [-Dd::from_f64(alpha[1]) * width, phase],
        ];
        (
            core::array::from_fn(|j| short_mean([projected[0][j], projected[1][j]], rates)),
            width.to_f64(),
        )
    } else {
        let diagonal = alpha.map(|a| decay_integral(a, width));
        let cross = cross_integral(alpha, beat, width)?;
        (
            core::array::from_fn(|j| {
                let [a, b] = [projected[0][j], projected[1][j]];
                let overlap = [dot(a, b), a[0] * b[1] - a[1] * b[0]];
                dot(a, a) * diagonal[0]
                    + dot(b, b) * diagonal[1]
                    + (overlap[0] * cross[0] - overlap[1] * cross[1]) * Dd::from_f64(2.0)
            }),
            1.0,
        )
    };
    let mut power = Dd::ZERO;
    for (j, integral) in integrals.into_iter().enumerate() {
        let loss = if j == 3 {
            response.loss_f_m[1] / norm.to_f64()
        } else {
            response.loss_f_m[0]
        };
        if loss == 0.0 {
            continue;
        }
        let integral = integral.to_f64();
        if !integral.is_finite() || integral < 0.0 {
            return Err("coherent crystal absorption integral is unresolved");
        }
        power = power
            + Dd::from_f64(physical_power(
                integral,
                length,
                loss,
                response.omega,
                area,
                field_scale,
            ));
    }
    let power = power.to_f64();
    if !power.is_finite() || power < 0.0 {
        return Err("crystal slab absorbed power is unrepresentable");
    }
    Ok(power)
}

pub(in super::super) fn attenuated(electric: Vector, alpha: f64, start: Dd) -> [Complex; 3] {
    let exposure = product(alpha, start);
    let remaining = decay::remaining_fraction(exposure);
    let factor = if exposure < 0.5 {
        Dd::ONE + Dd::from_f64(det::expm1(-exposure))
    } else {
        Dd::from_f64(remaining)
    };
    electric.map(|e| {
        [e.re, e.im].map(|v| {
            if remaining >= f64::MIN_POSITIVE {
                Dd::from_f64(v) * factor
            } else {
                Dd::from_f64(decay::attenuate(v, exposure))
            }
        })
    })
}

pub(in super::super) fn product(rate: f64, width: Dd) -> f64 {
    let leading = rate * width.hi;
    if leading.is_infinite() {
        leading
    } else {
        (Dd::from_f64(rate) * width).to_f64()
    }
}

pub(in super::super) fn decay_integral(alpha: f64, width: Dd) -> Dd {
    if alpha == 0.0 {
        return width;
    }
    let x = product(alpha, width);
    if x < 0.125 {
        // Form the mean attenuation before the depth factor. A subnormal
        // opacity has no representable relative correction, but rounding it
        // before dividing by alpha could erase or perturb a finite interval.
        let mean = if x < f64::MIN_POSITIVE {
            1.0
        } else {
            (-0.5 * det::expm1(-2.0 * x)) / x
        };
        return width * Dd::from_f64(mean);
    }
    // Divide after expm1 and before doubling the rate: a finite skin depth
    // need not have a representable 2*alpha or alpha*width.
    Dd::from_f64(-0.5 * det::expm1(-2.0 * x)) / Dd::from_f64(alpha)
}

pub(in super::super) fn cross_integral(alpha: [f64; 2], beat: Dd, width: Dd) -> Result<Complex> {
    let reference = alpha[0].max(alpha[1]).max(beat.to_f64().abs());
    if reference == 0.0 {
        return Ok([width, Dd::ZERO]);
    }
    let scale = Dd::from_f64(reference);
    let denominator = [
        Dd::from_f64(alpha[0]) / scale + Dd::from_f64(alpha[1]) / scale,
        -beat / scale,
    ];
    let exposure = product(alpha[0], width) + product(alpha[1], width);
    let remaining = decay::remaining_fraction(exposure);
    let numerator = if remaining == 0.0 {
        [Dd::ONE, Dd::ZERO]
    } else {
        let phase = beat * width;
        if !phase.hi.is_finite() || !phase.lo.is_finite() {
            return Err("crystal absorption beat integral is unrepresentable");
        }
        let (sine, _) = det::sin_cos_dd(phase);
        let half_sine = det::sin_cos_dd(phase * Dd::from_f64(0.5)).0;
        let attenuation = Dd::from_f64(remaining);
        [
            Dd::from_f64(decay::decayed_fraction(exposure))
                + Dd::from_f64(2.0) * attenuation * half_sine * half_sine,
            -attenuation * sine,
        ]
    };
    let norm = dot(denominator, denominator);
    let quotient = multiply(numerator, [denominator[0], -denominator[1]]);
    Ok(quotient.map(|v| (v / norm) / scale))
}

pub(in super::super) fn short_mean<const N: usize>(
    mut modes: [Complex; N],
    rates: [Complex; N],
) -> Dd {
    // |rate| < sqrt(2)/8. Per unit modal amplitude, the first omitted field
    // coefficient is < 1e-27;
    // its effect is below binary64 roundoff even when the constant field
    // cancels, because all higher coefficients retain the same mode sum.
    const TERMS: usize = 17;
    let mut coefficients = [ZERO; TERMS];
    for (n, coefficient) in coefficients.iter_mut().enumerate() {
        *coefficient = [0, 1].map(|part| modes.iter().fold(Dd::ZERO, |sum, mode| sum + mode[part]));
        let divisor = Dd::from_f64((n + 1) as f64);
        for m in 0..N {
            modes[m] = multiply(modes[m], rates[m]).map(|v| v / divisor);
        }
    }
    let mut integral = Dd::ZERO;
    for (n, a) in coefficients.iter().copied().enumerate() {
        integral = integral + dot(a, a) / Dd::from_f64((2 * n + 1) as f64);
        for (m, b) in coefficients[..n].iter().copied().enumerate() {
            integral = integral + dot(a, b) * Dd::from_f64(2.0) / Dd::from_f64((n + m + 1) as f64);
        }
    }
    integral
}

pub(in super::super) fn physical_power(
    integral: f64,
    length: f64,
    loss: f64,
    omega: f64,
    area: f64,
    field_scale: f64,
) -> f64 {
    if integral == 0.0 {
        return 0.0;
    }
    let mut direct = 1.0;
    let mut normal = true;
    for factor in [
        integral,
        loss,
        0.5 * omega,
        length,
        area,
        field_scale,
        field_scale,
    ] {
        direct *= factor;
        normal &= direct.is_normal();
    }
    if normal {
        return direct;
    }
    // Preserve a representable physical result when a very thin interval,
    // weak original loss and large field/area exceed an intermediate range.
    let logarithm = det::ln(integral) + det::ln(loss) + det::ln(omega) - core::f64::consts::LN_2
        + det::ln(length)
        + det::ln(area)
        + 2.0 * det::ln(field_scale);
    det::exp(logarithm)
}

pub(in super::super) fn multiply(a: Complex, b: Complex) -> Complex {
    [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]
}

pub(in super::super) fn dot(a: Complex, b: Complex) -> Dd {
    a[0] * b[0] + a[1] * b[1]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn thin_coherent_node_retains_cubic_volume_heat() {
        // Integral_0^1 |1-exp(i delta u)|² du = 2(1-sin(delta)/delta).
        // Its small-delta limit delta²/3 is below the cancellation error of
        // subtracting separately rounded diagonal and cross integrals.
        for delta in [1e-4, 1e-8, 1e-16] {
            let actual = short_mean(
                [[Dd::ONE, Dd::ZERO], [-Dd::ONE, Dd::ZERO]],
                [ZERO, [Dd::ZERO, Dd::from_f64(delta)]],
            )
            .to_f64();
            let expected = delta * delta * (1.0 / 3.0 - delta * delta / 60.0);
            assert!((actual - expected).abs() < 3e-15 * expected);
        }
    }

    #[test]
    fn four_coherent_modes_retain_a_third_order_field_node() {
        // |1-3 exp(i d u)+3 exp(2 i d u)-exp(3 i d u)|²
        // = |1-exp(i d u)|^6, whose leading mean is d^6/7.
        for delta in [2.0_f64.powi(-12), 2.0_f64.powi(-20), 2.0_f64.powi(-30)] {
            let actual = short_mean(
                [1.0, -3.0, 3.0, -1.0].map(|v| [Dd::from_f64(v), Dd::ZERO]),
                [0.0, delta, 2.0 * delta, 3.0 * delta].map(|v| [Dd::ZERO, Dd::from_f64(v)]),
            )
            .to_f64();
            let d2 = delta * delta;
            let expected = d2 * d2 * d2 * (1.0 / 7.0 - d2 / 36.0 + 7.0 * d2 * d2 / 2640.0);
            assert!((actual / expected - 1.0).abs() < 3e-13);
        }
    }

    #[test]
    fn underflowed_opacity_and_intermediate_products_keep_representable_heat() {
        for width in [1e-120, 1e-200] {
            assert_eq!(decay_integral(1e-200, Dd::from_f64(width)).to_f64(), width);
        }
        let power = physical_power(1e-200, 1.0, 1e-120, 2.0, 1.0, 1e100);
        assert!((power / 1e-120 - 1.0).abs() < 3e-13);
        let power = physical_power(1e-200, 1e-200, 1e-200, 1e200, 1e200, 1e100);
        assert!((power / 0.5 - 1.0).abs() < 3e-13);
    }
}
