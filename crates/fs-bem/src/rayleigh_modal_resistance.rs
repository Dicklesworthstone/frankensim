//! Stable low-frequency evaluation of the SAME point-area resistance Gram.
//!
//! sinc(k|x-y|) = average_S2 exp(-ik d.x) exp(+ik d.y); see DLMF
//! https://dlmf.nist.gov/10.54.E2 at n=0. Projecting nearly equal sinc
//! entries loses quadrupole and higher radiation before the modal projection.
//! Instead sum each signed Cartesian moment first and integrate products of
//! the resulting outgoing amplitudes with positive angular weights.
use super::*;
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
