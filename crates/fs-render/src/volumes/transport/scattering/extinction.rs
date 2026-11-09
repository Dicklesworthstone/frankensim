//! Shared extinction arithmetic for host matter and independent absorbers.
//!
//! Host coefficients are m²/kg; additional absorption is absolute 1/m and can
//! occupy a cell with zero host density. The same linear extinction determines
//! flight distances, spectral mixture densities and opaque source integrals.
//! Inputs have already passed the original optical snapshot admission.

use fs_math::det;

pub(crate) fn log_total_linear_extinction(
    density: f64,
    host_extinction: f64,
    additional_absorption: f64,
) -> f64 {
    let host = if density == 0.0 || host_extinction == 0.0 {
        f64::NEG_INFINITY
    } else {
        det::ln(density) + det::ln(host_extinction)
    };
    if additional_absorption == 0.0 {
        return host;
    }
    let additional = det::ln(additional_absorption);
    let high = host.max(additional);
    let low = host.min(additional);
    high + det::ln(1.0 + det::exp(low - high))
}

pub(crate) fn event_distance_with_additional(
    optical_depth: f64,
    density: f64,
    host_extinction: f64,
    additional_absorption: f64,
) -> Option<f64> {
    if optical_depth == 0.0 {
        return Some(0.0);
    }
    let extinction = density * host_extinction + additional_absorption;
    let direct = optical_depth / extinction;
    let distance = if extinction.is_finite()
        && extinction >= f64::MIN_POSITIVE
        && direct.is_finite()
        && direct > 0.0
    {
        direct
    } else {
        det::exp(
            det::ln(optical_depth)
                - log_total_linear_extinction(density, host_extinction, additional_absorption),
        )
    };
    (distance.is_finite() && distance > 0.0).then_some(distance)
}

/// Original host-only mixture: the common positive density cancels exactly.
pub(crate) fn collision_weights<const N: usize>(
    depth: &[f64; N],
    extinction: &[f64; N],
    scattering: &[f64; N],
) -> Option<[f64; N]> {
    weights_from_logs(depth, &extinction.map(det::ln), &scattering.map(det::ln))
}

pub(crate) fn collision_weights_with_additional<const N: usize>(
    depth: &[f64; N],
    density: f64,
    host_extinction: &[f64; N],
    host_scattering: &[f64; N],
    additional_absorption: &[f64; N],
) -> Option<[f64; N]> {
    if density > 0.0 && additional_absorption.iter().all(|&value| value == 0.0) {
        return collision_weights(depth, host_extinction, host_scattering);
    }
    let log_extinction = std::array::from_fn(|lane| {
        log_total_linear_extinction(density, host_extinction[lane], additional_absorption[lane])
    });
    let log_density = det::ln(density);
    let log_scattering = host_scattering.map(|value| log_density + det::ln(value));
    weights_from_logs(depth, &log_extinction, &log_scattering)
}

fn weights_from_logs<const N: usize>(
    depth: &[f64; N],
    log_extinction: &[f64; N],
    log_scattering: &[f64; N],
) -> Option<[f64; N]> {
    let log_density = std::array::from_fn::<_, N, _>(|lane| log_extinction[lane] - depth[lane]);
    let peak = log_density.into_iter().fold(f64::NEG_INFINITY, f64::max);
    if !peak.is_finite() {
        return None;
    }
    let sum: f64 = log_density
        .iter()
        .map(|value| det::exp(*value - peak))
        .sum();
    let weights = std::array::from_fn(|lane| {
        if log_scattering[lane] == f64::NEG_INFINITY {
            0.0
        } else {
            det::exp(log_scattering[lane] - depth[lane] - peak) * (N as f64 / sum)
        }
    });
    weights
        .iter()
        .all(|weight| weight.is_finite())
        .then_some(weights)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn relative(actual: f64, expected: f64) {
        assert!(
            (actual / expected - 1.0).abs() < 4e-12,
            "{actual:e} != {expected:e}"
        );
    }

    #[test]
    fn spectral_absorption_changes_the_actual_flight_and_mixture_density() {
        let distance = event_distance_with_additional(0.35, 2.0, 0.5, 0.4).unwrap();
        relative(distance, 0.25);
        let beta = [1.4_f64, 1.6];
        let sigma = [0.6, 0.2];
        let depth = beta.map(|value| 0.25 * value);
        let weights =
            collision_weights_with_additional(&depth, 2.0, &[0.5, 0.8], &[0.3, 0.1], &[0.4, 0.0])
                .unwrap();
        let proposal = (beta[0] * (-depth[0]).exp() + beta[1] * (-depth[1]).exp()) / 2.0;
        for lane in 0..2 {
            relative(weights[lane], sigma[lane] * (-depth[lane]).exp() / proposal);
        }
        assert_eq!(
            collision_weights_with_additional(&depth, 2.0, &[0.5, 0.8], &[0.3, 0.1], &[0.0; 2]),
            collision_weights(&depth, &[0.5, 0.8], &[0.3, 0.1]),
        );
    }

    #[test]
    fn independent_absorber_in_zero_density_host_has_no_scattering_continuation() {
        relative(
            event_distance_with_additional(0.35, 0.0, 1.0, 0.7).unwrap(),
            0.5,
        );
        relative(log_total_linear_extinction(0.0, 1.0, 0.7).exp(), 0.7);
        assert_eq!(
            collision_weights_with_additional(&[0.35, 0.1], 0.0, &[1.0; 2], &[1.0; 2], &[0.7, 0.2]),
            Some([0.0; 2]),
        );
    }

    #[test]
    fn overflowing_linear_extinction_retains_finite_distance_weight_and_source_scale() {
        let log_beta = log_total_linear_extinction(1e155, 1e154, 1e308);
        relative((1e300_f64.ln() - log_beta).exp(), 1e-9 / 1.1);
        relative(
            event_distance_with_additional(0.35, 1e155, 1e154, 1e308).unwrap(),
            0.35e-309 / 1.1,
        );
        let weights =
            collision_weights_with_additional(&[0.35], 1e155, &[1e154], &[2e153], &[1e308])
                .unwrap();
        relative(weights[0], 2.0 / 11.0);
    }
}
