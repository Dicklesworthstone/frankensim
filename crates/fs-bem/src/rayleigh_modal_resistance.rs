//! Stable low-frequency evaluation of the SAME point-area resistance Gram.
//!
//! sinc(k|x-y|) = average_S2 exp(-ik d.x) exp(+ik d.y); see DLMF
//! https://dlmf.nist.gov/10.54.E2 at n=0. Projecting nearly equal sinc
//! entries loses quadrupole and higher radiation before the modal projection.
//! Instead sum each signed Cartesian moment first and integrate products of
//! the resulting outgoing amplitudes with positive angular weights.
use super::*;
use fs_math::dd::Dd;
use fs_math::eft::two_sum;

const DEGREE: usize = 16;
const ORDERS: usize = DEGREE + 1;
const MOMENTS: usize = (DEGREE + 1) * (DEGREE + 2) / 2;
const DIRECTIONS: usize = (DEGREE / 2 + 1) * (DEGREE + 1);

fn add(sum: &mut f64, correction: &mut f64, term: f64) {
    let (next, error) = two_sum(*sum, term);
    *sum = next;
    *correction += error;
}

fn moments(surface: &PreparedSurface) -> Result<Vec<f64>, HelmholtzError> {
    let n = surface.modes;
    let mut sum = reserved(MOMENTS * n, "Rayleigh signed source moments")?;
    sum.resize(MOMENTS * n, 0.0);
    let mut correction = reserved(sum.len(), "Rayleigh source-moment corrections")?;
    correction.resize(sum.len(), 0.0);
    let mut lo = surface.positions[0];
    let mut hi = lo;
    for point in &surface.positions {
        for axis in 0..2 {
            lo[axis] = lo[axis].min(point[axis]);
            hi[axis] = hi[axis].max(point[axis]);
        }
    }
    // Translate before taking powers: global coordinates must not set the
    // expansion scale. Every normalized point lies within the source diameter.
    let origin = [lo[0].midpoint(hi[0]), lo[1].midpoint(hi[1])];
    for (cell, point) in surface.positions.iter().enumerate() {
        let x = (point[0] - origin[0]) / surface.maximum_distance_m;
        let y = (point[1] - origin[1]) / surface.maximum_distance_m;
        let mut xp = [1.0; DEGREE + 1];
        let mut yp = xp;
        for order in 1..=DEGREE {
            xp[order] = xp[order - 1] * x;
            yp[order] = yp[order - 1] * y;
        }
        for order in 0..=DEGREE {
            for q in 0..=order {
                let offset = (order * (order + 1) / 2 + q) * n;
                let monomial = xp[order - q] * yp[q];
                for mode in 0..n {
                    add(
                        &mut sum[offset + mode],
                        &mut correction[offset + mode],
                        monomial * surface.weighted_shape[cell * n + mode],
                    );
                }
            }
        }
    }
    for (value, error) in sum.iter_mut().zip(correction) {
        *value += error;
        if !value.is_finite() {
            return Err(invalid(
                "Rayleigh source moment exceeds finite numerical range",
            ));
        }
    }
    Ok(sum)
}

struct Direction {
    x: f64,
    y: f64,
    weight: f64,
}

fn directions() -> &'static [Direction] {
    static RULE: OnceLock<Vec<Direction>> = OnceLock::new();
    RULE.get_or_init(|| {
        let (nodes, weights) = super::super::gauss_legendre(DEGREE + 1);
        let mut rule = Vec::with_capacity(DIRECTIONS);
        // Degree-16 amplitudes have degree-32 products. A 17-point polar
        // Gauss rule and 34 periodic azimuth points integrate these products
        // exactly. The planar source is even in polar cosine; Re(Fm conj Fn)
        // is also pi-periodic in azimuth. Retain one of each equal pair.
        for i in 0..=DEGREE / 2 {
            let transverse = det::sqrt(1.0 - nodes[i] * nodes[i]);
            let polar_weight = if i == DEGREE / 2 {
                0.5 * weights[i]
            } else {
                weights[i]
            };
            for j in 0..=DEGREE {
                let angle = PI * j as f64 / (DEGREE + 1) as f64;
                rule.push(Direction {
                    x: transverse * det::cos(angle),
                    y: transverse * det::sin(angle),
                    weight: polar_weight / (DEGREE + 1) as f64,
                });
            }
        }
        rule
    })
}

fn coefficients(surface: &PreparedSurface) -> Result<Vec<f64>, HelmholtzError> {
    let moments = moments(surface)?;
    let n = surface.modes;
    let mut coefficients = reserved(
        DIRECTIONS * ORDERS * n,
        "Rayleigh directional Taylor coefficients",
    )?;
    coefficients.resize(DIRECTIONS * ORDERS * n, 0.0);
    for (index, direction) in directions().iter().enumerate() {
        let mut xp = [1.0; ORDERS];
        let mut yp = xp;
        for order in 1..ORDERS {
            xp[order] = xp[order - 1] * direction.x / order as f64;
            yp[order] = yp[order - 1] * direction.y / order as f64;
        }
        for order in 0..ORDERS {
            let sign = if matches!(order % 4, 0 | 3) {
                1.0
            } else {
                -1.0
            };
            let output = (index * ORDERS + order) * n;
            let mut correction = [0.0; MAX_MODAL_RAYLEIGH_FIELDS];
            for q in 0..=order {
                let coefficient = sign * xp[order - q] * yp[q];
                let input = (order * (order + 1) / 2 + q) * n;
                for mode in 0..n {
                    add(
                        &mut coefficients[output + mode],
                        &mut correction[mode],
                        coefficient * moments[input + mode],
                    );
                }
            }
            for mode in 0..n {
                let value = &mut coefficients[output + mode];
                *value += correction[mode];
                if !value.is_finite() {
                    return Err(invalid(
                        "Rayleigh directional coefficient exceeds finite numerical range",
                    ));
                }
            }
        }
    }
    Ok(coefficients)
}

pub(super) fn matrix(
    surface: &PreparedSurface,
    k: f64,
    scale: f64,
) -> Result<Vec<f64>, HelmholtzError> {
    let cells = surface.count as u128;
    let fields = surface.modes as u128;
    // Conservatively include both geometry projections on every admission,
    // even when the cache exists. Ordinary reactance work also remains;
    // no field is dropped to fit. At the field cap the coefficient cache
    // occupies 153 * 17 * 128 doubles, about 2.54 MiB.
    admit_modal_work(
        cells * fields * (cells + fields)
            + MOMENTS as u128 * cells * fields
            + DIRECTIONS as u128 * fields * (MOMENTS as u128 + ORDERS as u128 + fields),
    )?;
    let coefficients = surface
        .resistance_coefficients
        .get_or_init(|| coefficients(surface))
        .as_ref()
        .map_err(Clone::clone)?;
    let n = surface.modes;
    let mut sum = reserved(n * n, "Rayleigh angular resistance Gram")?;
    sum.resize(n * n, 0.0);
    let mut correction = reserved(n * n, "Rayleigh angular Gram corrections")?;
    correction.resize(n * n, 0.0);
    let kd = k * surface.maximum_distance_m;
    let mut powers = [1.0; ORDERS];
    for order in 1..ORDERS {
        powers[order] = powers[order - 1] * kd;
    }
    // Geometry, direction cosines and factorials are already contracted.
    // Each frequency needs 17 terms per direction/field rather than 153.
    for (index, direction) in directions().iter().enumerate() {
        let (mut re, mut im, mut cr, mut ci) = (
            [0.0; MAX_MODAL_RAYLEIGH_FIELDS],
            [0.0; MAX_MODAL_RAYLEIGH_FIELDS],
            [0.0; MAX_MODAL_RAYLEIGH_FIELDS],
            [0.0; MAX_MODAL_RAYLEIGH_FIELDS],
        );
        for (order, power) in powers.iter().enumerate() {
            let offset = (index * ORDERS + order) * n;
            for mode in 0..n {
                let term = power * coefficients[offset + mode];
                if order % 2 == 0 {
                    add(&mut re[mode], &mut cr[mode], term);
                } else {
                    add(&mut im[mode], &mut ci[mode], term);
                }
            }
        }
        for mode in 0..n {
            re[mode] += cr[mode];
            im[mode] += ci[mode];
        }
        for row in 0..n {
            for column in row..n {
                let index = row * n + column;
                let value = (re[row] * re[column] + im[row] * im[column]) * direction.weight;
                add(&mut sum[index], &mut correction[index], value);
            }
        }
    }
    for row in 0..n {
        for column in row..n {
            let index = row * n + column;
            let value = scale * (sum[index] + correction[index]);
            if !value.is_finite() {
                return Err(invalid(
                    "Rayleigh angular resistance exceeds finite numerical range",
                ));
            }
            sum[index] = value;
            sum[column * n + row] = value;
        }
    }
    Ok(sum)
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), HelmholtzError> {
    if cancelled() {
        Err(HelmholtzError::Cancelled)
    } else {
        Ok(())
    }
}

fn squared_pair(even: [Dd; 2], odd: [Dd; 2]) -> Dd {
    // For arbitrary COMPLEX modal velocities, |F(d)|^2 alone is not
    // pi-periodic. Average the two opposite azimuths analytically:
    // (|E+iO|^2 + |E-iO|^2)/2 = |E|^2 + |O|^2.
    even[0] * even[0] + even[1] * even[1] + odd[0] * odd[0] + odd[1] * odd[1]
}

fn compact_power(
    surface: &PreparedSurface,
    kd: f64,
    velocity: &[C64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Dd, HelmholtzError> {
    let cells = surface.count as u128;
    let fields = surface.modes as u128;
    admit_modal_work(
        MOMENTS as u128 * cells * fields + DIRECTIONS as u128 * fields * (MOMENTS + ORDERS) as u128,
    )?;
    let coefficients = surface
        .resistance_coefficients
        .get_or_init(|| coefficients(surface))
        .as_ref()
        .map_err(Clone::clone)?;
    let mut powers = [Dd::ONE; ORDERS];
    for degree in 1..ORDERS {
        powers[degree] = powers[degree - 1] * Dd::from_f64(kd);
    }
    let mut result = Dd::ZERO;
    for (index, direction) in directions().iter().enumerate() {
        poll(cancelled)?;
        let (mut even, mut odd) = ([Dd::ZERO; 2], [Dd::ZERO; 2]);
        for (degree, &power) in powers.iter().enumerate() {
            let offset = (index * ORDERS + degree) * surface.modes;
            let mut amplitude = [Dd::ZERO; 2];
            for (mode, v) in velocity.iter().enumerate() {
                let coefficient = Dd::from_f64(coefficients[offset + mode]);
                amplitude[0] = amplitude[0] + coefficient * Dd::from_f64(v.re);
                amplitude[1] = amplitude[1] + coefficient * Dd::from_f64(v.im);
            }
            let parity = if degree % 2 == 0 { &mut even } else { &mut odd };
            for axis in 0..2 {
                parity[axis] = parity[axis] + amplitude[axis] * power;
            }
        }
        result = result + squared_pair(even, odd) * Dd::from_f64(direction.weight);
    }
    Ok(result)
}

// An n-point polar / 2n-point azimuth rule integrates the squared
// degree-(n-1) plane-wave polynomial exactly. For real phase |t|<=kd,
// |exp(it)-P_(n-1)(it)| <= kd^n/n!, from the integral Taylor remainder.
fn angular_order(kd: f64) -> Result<(usize, f64), HelmholtzError> {
    if kd > MAX_RAYLEIGH_POWER_POLAR_ORDER as f64 {
        return Err(invalid(
            "Rayleigh outgoing power requires a higher angular order",
        ));
    }
    let mut remainder = 1.0;
    for order in 1..=MAX_RAYLEIGH_POWER_POLAR_ORDER {
        remainder *= kd / order as f64;
        if order >= ORDERS && remainder <= 1.0e-18 {
            return Ok((order, remainder));
        }
    }
    Err(invalid(
        "Rayleigh outgoing power requires a higher angular order",
    ))
}

fn resolved_power(
    surface: &PreparedSurface,
    k: f64,
    velocity: &[C64],
    order: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Dd, HelmholtzError> {
    let directions = order.div_ceil(2) * order;
    admit_modal_work(surface.count as u128 * (surface.modes as u128 + 8 * directions as u128))?;
    let mut flow = reserved(surface.count, "coherent Rayleigh cell volume velocities")?;
    let mut lo = surface.positions[0];
    let mut hi = lo;
    for (cell, point) in surface.positions.iter().enumerate() {
        if cell % 64 == 0 {
            poll(cancelled)?;
        }
        let mut value = [Dd::ZERO; 2];
        for (mode, v) in velocity.iter().enumerate() {
            let weight = Dd::from_f64(surface.weighted_shape[cell * surface.modes + mode]);
            value[0] = value[0] + weight * Dd::from_f64(v.re);
            value[1] = value[1] + weight * Dd::from_f64(v.im);
        }
        flow.push(value);
        for axis in 0..2 {
            lo[axis] = lo[axis].min(point[axis]);
            hi[axis] = hi[axis].max(point[axis]);
        }
    }
    let origin = [lo[0].midpoint(hi[0]), lo[1].midpoint(hi[1])];
    let (nodes, weights) = super::super::gauss_legendre(order);
    let mut result = Dd::ZERO;
    for polar in 0..order.div_ceil(2) {
        let transverse = det::sqrt(1.0 - nodes[polar] * nodes[polar]);
        let weight = if order % 2 == 1 && polar == order / 2 {
            0.5 * weights[polar]
        } else {
            weights[polar]
        } / order as f64;
        for azimuth in 0..order {
            poll(cancelled)?;
            let angle = PI * azimuth as f64 / order as f64;
            let direction = [transverse * det::cos(angle), transverse * det::sin(angle)];
            let (mut even, mut odd) = ([Dd::ZERO; 2], [Dd::ZERO; 2]);
            for (point, value) in surface.positions.iter().zip(&flow) {
                let phase = -k
                    * (direction[0] * (point[0] - origin[0])
                        + direction[1] * (point[1] - origin[1]));
                let (sine, cosine) = det::sin_cos(phase);
                for axis in 0..2 {
                    even[axis] = even[axis] + value[axis] * Dd::from_f64(cosine);
                    odd[axis] = odd[axis] + value[axis] * Dd::from_f64(sine);
                }
            }
            result = result + squared_pair(even, odd) * Dd::from_f64(weight);
        }
    }
    Ok(result)
}

pub(super) fn coherent_power(
    surface: &PreparedSurface,
    omega: f64,
    medium: Medium,
    velocity: &[C64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<CoherentRadiatedPower, HelmholtzError> {
    poll(cancelled)?;
    let ppw = surface.equivalent_disc_ppw(omega, medium)?;
    if ppw < MIN_PANELS_PER_WAVELENGTH {
        return Err(HelmholtzError::TooCoarse {
            panels_per_wavelength: ppw,
        });
    }
    if velocity.len() != surface.modes
        || velocity
            .iter()
            .any(|v| !v.re.is_finite() || !v.im.is_finite())
    {
        return Err(invalid(
            "Rayleigh outgoing power requires every finite modal velocity",
        ));
    }
    let k = omega / medium.sound_speed;
    let scale = (0.5 * medium.density * omega / (2.0 * PI)) * k;
    let kd = k * surface.maximum_distance_m;
    if !scale.is_finite() || scale <= 0.0 || !kd.is_finite() {
        return Err(invalid(
            "Rayleigh outgoing power scale exceeds finite numerical range",
        ));
    }
    let mut magnitudes = [0.0; MAX_MODAL_RAYLEIGH_FIELDS];
    for (mode, v) in velocity.iter().enumerate() {
        magnitudes[mode] = det::hypot(v.re, v.im);
    }
    let mut absolute_amplitude = Dd::ZERO;
    for cell in 0..surface.count {
        if cell % 64 == 0 {
            poll(cancelled)?;
        }
        for (mode, &magnitude) in magnitudes[..surface.modes].iter().enumerate() {
            absolute_amplitude = absolute_amplitude
                + Dd::from_f64(surface.weighted_shape[cell * surface.modes + mode].abs())
                    * Dd::from_f64(magnitude);
        }
    }
    if !absolute_amplitude.to_f64().is_finite() {
        return Err(invalid(
            "Rayleigh source amplitude exceeds finite numerical range",
        ));
    }
    let (value, remainder, directions, quadrature_factor) = if surface.count == 1 {
        let mut value = [Dd::ZERO; 2];
        for (mode, v) in velocity.iter().enumerate() {
            let weight = Dd::from_f64(surface.weighted_shape[mode]);
            value[0] = value[0] + weight * Dd::from_f64(v.re);
            value[1] = value[1] + weight * Dd::from_f64(v.im);
        }
        (squared_pair(value, [Dd::ZERO; 2]), 0.0, 1, 1.0)
    } else if kd <= 0.5 {
        let mut remainder = 1.0;
        for order in 1..=ORDERS {
            remainder *= kd / order as f64;
        }
        (
            compact_power(surface, kd, velocity, cancelled)?,
            remainder,
            DIRECTIONS,
            1.0,
        )
    } else {
        let (order, remainder) = angular_order(kd)?;
        (
            resolved_power(surface, k, velocity, order, cancelled)?,
            remainder,
            order.div_ceil(2) * order,
            2.0,
        )
    };
    let power_w = (value * Dd::from_f64(scale)).to_f64();
    // Compact evaluation integrates the polynomial itself. Exact-exponential
    // quadrature has two remainder terms: the integral and the sampled rule.
    let angular_truncation_bound_w = (absolute_amplitude
        * absolute_amplitude
        * Dd::from_f64(remainder * (2.0 + remainder) * quadrature_factor)
        * Dd::from_f64(scale))
    .to_f64();
    if !power_w.is_finite() || power_w < 0.0 || !angular_truncation_bound_w.is_finite() {
        return Err(invalid(
            "Rayleigh outgoing power exceeds finite numerical range",
        ));
    }
    poll(cancelled)?;
    Ok(CoherentRadiatedPower {
        power_w,
        angular_truncation_bound_w,
        angular_directions: directions,
    })
}
