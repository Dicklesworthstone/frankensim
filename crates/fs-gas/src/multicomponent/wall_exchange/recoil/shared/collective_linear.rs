//! Exact common-linear-rate collective motion and its physical heat integrals.
//!
//! The gas mean/support slip decays at a*(1+m/M); gas deviations decay at a.
//! Positive Gram factors evaluate each cell's squared-slip integral without
//! subtracting modal energies. Work is integrated in the initial support frame.

use super::{SharedRecoilCell, SlidingSupport, Stage, Sum, WallExchangeError, poll, positive};
use fs_math::det;

pub(super) fn advance(
    cells: &[SharedRecoilCell],
    states: &mut [Stage],
    dt: f64,
    support: SlidingSupport,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Option<SlidingSupport>, WallExchangeError> {
    let mut rate = None;
    let mut total_mass = Sum::default();
    for cell in cells {
        poll(cancelled)?;
        if cell.linear_rate_per_s == 0.0 && cell.quadratic_per_m == 0.0 {
            continue;
        }
        if cell.quadratic_per_m != 0.0 || rate.is_some_and(|a| a != cell.linear_rate_per_s) {
            return Ok(None);
        }
        rate = Some(cell.linear_rate_per_s);
        total_mass.add(cell.density_kg_m3 * cell.volume_m3);
    }
    let Some(rate) = rate else { return Ok(None) };
    let mass = total_mass.value();
    let ratio = mass / support.mass_kg;
    let z = rate * dt;
    // Preserve the general owner for unrepresentable aggregate coordinates.
    // Positive infinite exponents are the resolved fully relaxed limit.
    if !positive(mass) || !positive(ratio) || z == 0.0 {
        return Ok(None);
    }
    let gas_weight = 1.0 / (1.0 + ratio);
    let support_weight = ratio / (1.0 + ratio);
    let mut mean = Sum::default();
    let mut laboratory_mean = Sum::default();
    for cell in cells {
        poll(cancelled)?;
        if cell.linear_rate_per_s == 0.0 {
            continue;
        }
        let weight = (cell.density_kg_m3 * cell.volume_m3) / mass;
        let velocity = cell.momentum_kg_m2_s / cell.density_kg_m3;
        mean.add(weight * (velocity - support.velocity_m_s));
        laboratory_mean.add(weight * velocity);
    }
    let mean = mean.value();
    let laboratory_mean = laboratory_mean.value();
    let support_scale = support_weight * mean;
    // Keep ratio independently: (1+ratio)-1 loses a weak, still resolved
    // coupling to a heavy support and erases initially-comoving cell motion.
    let contrast = product([ratio, rate, dt]);
    let fast = z + contrast;
    let slow_removed = -det::expm1(-z);
    let fast_removed = -det::expm1(-fast);
    let slow_factor = det::exp(-z);
    let fast_factor = det::exp(-fast);
    let support_delta = support_scale * fast_removed;
    let next_support = SlidingSupport {
        velocity_m_s: support.velocity_m_s + support_delta,
        ..support
    };
    let moments = Moments::new(z, contrast, fast, ratio, gas_weight)?;
    for (cell, state) in cells.iter().zip(states) {
        poll(cancelled)?;
        if cell.linear_rate_per_s == 0.0 {
            continue;
        }
        let velocity = cell.momentum_kg_m2_s / cell.density_kg_m3;
        let slip = velocity - support.velocity_m_s;
        let deviation = slip - mean;
        let delta = least_cancellation(
            -deviation * slow_removed,
            -product([gas_weight, mean, fast_removed]),
            -slip * slow_removed,
            product([
                support_scale,
                moments.feedback[0],
                moments.feedback[1],
                moments.feedback[2],
            ]),
        );
        let momentum = if slow_removed <= 0.5 {
            cell.momentum_kg_m2_s + cell.density_kg_m3 * delta
        } else if support_weight >= 0.5 {
            // Form the center from its closer gas endpoint. This retains a
            // small final velocity beside a large initial support slip.
            let center = laboratory_mean - gas_weight * mean;
            cell.density_kg_m3
                * compensated([
                    center,
                    product([gas_weight, mean, fast_factor]),
                    deviation * slow_factor,
                ])
        } else if slow_factor == 0.0 {
            // A finite input can have an overflowing decay exponent while
            // its small mass-contrast exponent remains finite. The decayed
            // slip is exactly zero; do not form an irrelevant 0*infinity.
            cell.density_kg_m3 * next_support.velocity_m_s
        } else {
            // The contrast is evaluated independently even if 1+ratio rounded
            // to one; otherwise a comoving cell loses its induced motion.
            let contrast_factor = if contrast <= 1.0 {
                fast * exposure_ratio(contrast)
            } else {
                ((1.0 + ratio) / ratio) * (-det::expm1(-contrast))
            };
            cell.density_kg_m3
                * compensated([
                    support.velocity_m_s,
                    support_delta,
                    slow_factor * compensated([slip, -support_scale * contrast_factor]),
                ])
        };
        let gas_mass = cell.density_kg_m3 * cell.volume_m3;
        let heat_offset = support_scale * moments.heat_scale;
        let relative = compensated([slip, -heat_offset * moments.heat_mean]);
        let root_mass = det::sqrt(gas_mass);
        let first = product([root_mass, moments.heat_root, relative]);
        let second = product([root_mass, moments.variance_root, heat_offset]);
        let heat = compensated([first * first, second * second]);
        let a = moments.work_slow;
        let b = moments.work_induced;
        // Choose equivalent modal/slip forms by their term magnitudes. The
        // slip form retains O(h^3) work in an initially comoving cell; the
        // modal form retains weak common-mode work for a very light support.
        let work = least_cancellation(
            product([-gas_mass, support_scale, deviation, a[0], a[1], a[2]]),
            product([
                -gas_mass,
                support_scale,
                mean,
                gas_weight,
                0.5,
                fast_removed,
                fast_removed,
            ]),
            product([-gas_mass, support_scale, slip, a[0], a[1], a[2]]),
            product([
                gas_mass,
                support_scale,
                support_scale,
                b[0],
                b[1],
                b[2],
                b[3],
            ]),
        );
        if ![momentum, heat, work].into_iter().all(f64::is_finite) || heat < 0.0 {
            return Err(WallExchangeError::Invalid(
                "collective linear recoil endpoint",
            ));
        }
        state.momentum = momentum;
        state.work.add(work);
        state.heat.add(heat);
    }
    if ![mean, laboratory_mean, next_support.velocity_m_s]
        .into_iter()
        .all(f64::is_finite)
    {
        return Err(WallExchangeError::Invalid("collective linear recoil mean"));
    }
    poll(cancelled)?;
    Ok(Some(next_support))
}

struct Moments {
    heat_root: f64,
    heat_mean: f64,
    variance_root: f64,
    heat_scale: f64,
    feedback: [f64; 3],
    work_slow: [f64; 3],
    work_induced: [f64; 4],
}

impl Moments {
    fn new(
        z: f64,
        contrast: f64,
        fast: f64,
        ratio: f64,
        gas_weight: f64,
    ) -> Result<Self, WallExchangeError> {
        let result = if fast < 0.5 {
            Self::short(z, contrast, fast)
        } else {
            let e = det::exp(-z);
            let fast_e = det::exp(-fast);
            let removed = -det::expm1(-z);
            let h = -0.5 * det::expm1(-2.0 * z);
            let pz = z + fast;
            let c = (gas_weight / (1.0 + gas_weight)) * (-det::expm1(-pz));
            let exposure = if contrast <= 1.0 {
                z * exposure_ratio(contrast)
            } else {
                (-det::expm1(-contrast)) / ratio
            };
            let weighted = |exponent: f64| {
                let factor = det::exp(-exponent);
                if factor == 0.0 {
                    0.0
                } else {
                    factor * exposure
                }
            };
            let f1 = if e == 0.0 { 0.0 } else { e * exposure };
            let n = compensated([removed * (1.0 + fast_e), -2.0 * f1]);
            let plus = c + f1;
            let qp = 1.0 / (1.0 + gas_weight);
            let feedback = removed - f1;
            // V/(support_weight^2) is factored before squaring the mean.
            // Its numerator is conditioned for fast>=1/2, including ratio~0.
            let variance_root = 0.5 * det::sqrt(qp) * det::sqrt(n / h) * det::sqrt(plus);
            Self {
                heat_root: det::sqrt(h),
                heat_mean: qp * ((h - weighted(2.0 * z)) / h),
                variance_root,
                heat_scale: 1.0,
                feedback: [1.0, 1.0, feedback],
                work_slow: [1.0, 1.0, removed - c],
                work_induced: [1.0, 1.0, 1.0, feedback - 0.5 * (c - weighted(pz))],
            }
        };
        if ![
            result.heat_root,
            result.heat_mean,
            result.variance_root,
            result.heat_scale,
            result.feedback[2],
            result.work_slow[2],
            result.work_induced[3],
        ]
        .into_iter()
        .all(positive)
        {
            return Err(WallExchangeError::Invalid(
                "collective linear moment resolution",
            ));
        }
        Ok(result)
    }

    fn short(z: f64, contrast: f64, fast: f64) -> Self {
        // Eight-point Gauss-Legendre evaluates the SAME exponential moments.
        // With fast<1/2 the total exponent is <=1; its scalar remainder
        // coefficient is 8!^4/(17*(16!)^3) < 1.7e-23. Positive weights and
        // centered variance retain the first nonzero induced heat at small h.
        const NODES: [(f64, f64); 8] = [
            (0.019855071751231884, 0.05061426814518813),
            (0.10166676129318664, 0.11119051722668724),
            (0.2372337950418355, 0.15685332293894365),
            (0.4082826787521751, 0.181341891689181),
            (0.5917173212478248, 0.181341891689181),
            (0.7627662049581645, 0.15685332293894365),
            (0.8983332387068134, 0.11119051722668724),
            (0.9801449282487681, 0.05061426814518813),
        ];
        let mut h = Sum::default();
        let mut mean = Sum::default();
        let mut feedback = Sum::default();
        let mut slow_work = Sum::default();
        let mut induced_work = Sum::default();
        let mut weights = [0.0; 8];
        let mut contrasts = [0.0; 8];
        for (i, (node, weight)) in NODES.into_iter().enumerate() {
            let e = det::exp(-z * node);
            let d = node * exposure_ratio(contrast * node);
            let f = node * exposure_ratio(fast * node);
            weights[i] = weight * e * e;
            contrasts[i] = d;
            h.add(weights[i]);
            mean.add(weights[i] * d);
            feedback.add(weight * e * d);
            slow_work.add(weight * e * f);
            induced_work.add(weight * e * f * d);
        }
        let h = h.value();
        let mean = mean.value() / h;
        let mut variance = Sum::default();
        for i in 0..8 {
            let centered = contrasts[i] - mean;
            variance.add(weights[i] * centered * centered);
        }
        let root_z = det::sqrt(z);
        Self {
            heat_root: root_z * det::sqrt(h),
            heat_mean: mean,
            variance_root: root_z * det::sqrt(variance.value()),
            heat_scale: fast,
            feedback: [z, fast, feedback.value()],
            work_slow: [z, fast, slow_work.value()],
            work_induced: [z, fast, fast, induced_work.value()],
        }
    }
}

fn exposure_ratio(x: f64) -> f64 {
    if x == 0.0 { 1.0 } else { -det::expm1(-x) / x }
}

fn compensated<const N: usize>(values: [f64; N]) -> f64 {
    let mut result = Sum::default();
    for value in values {
        result.add(value);
    }
    result.value()
}

fn least_cancellation(a: f64, b: f64, c: f64, d: f64) -> f64 {
    if a.abs().max(b.abs()) <= c.abs().max(d.abs()) {
        compensated([a, b])
    } else {
        compensated([c, d])
    }
}

// Most physical scales take only the ordinary product. Normalize binary
// exponents only if an intermediate over/underflow could hide a representable
// complete product (notably small-time work times a large gas mass).
fn product<const N: usize>(values: [f64; N]) -> f64 {
    let mut result = 1.0;
    let mut scaled = false;
    for value in values {
        result *= value;
        scaled |= !result.is_finite() || result.abs() < f64::MIN_POSITIVE;
    }
    if !scaled || values.iter().any(|x| !x.is_finite()) {
        return result;
    }
    if values.contains(&0.0) {
        return 0.0;
    }
    let mut mantissa = 1.0;
    let mut exponent = 0_i32;
    let mut negative = false;
    for value in values {
        negative ^= value.is_sign_negative();
        let mut magnitude = value.abs();
        let mut adjustment = 0;
        if magnitude < f64::MIN_POSITIVE {
            magnitude *= 18_014_398_509_481_984.0;
            adjustment = -54;
        }
        let bits = magnitude.to_bits();
        exponent += ((bits >> 52) & 0x7ff) as i32 - 1023 + adjustment;
        mantissa *= f64::from_bits((bits & 0x000f_ffff_ffff_ffff) | (1023_u64 << 52));
    }
    while mantissa >= 2.0 {
        mantissa *= 0.5;
        exponent += 1;
    }
    let value = if exponent > 1023 {
        f64::INFINITY
    } else if exponent < -1075 {
        0.0
    } else if exponent < -1022 {
        (mantissa * f64::from_bits(((exponent + 1022 + 1023) as u64) << 52)) * f64::MIN_POSITIVE
    } else {
        mantissa * f64::from_bits(((exponent + 1023) as u64) << 52)
    };
    if negative { -value } else { value }
}
