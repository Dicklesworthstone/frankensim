//! G1/G2 analytic piston qualification of the shared time-domain owner.
//! No recordings or other piston implementations enter these oracles.
use super::*;
use fs_material::gas::{GasSpec, GasState};
use std::f64::consts::{PI, TAU};
use std::time::Instant;

fn medium() -> RayleighMedium {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap();
    RayleighMedium {
        density: gas.density,
        sound_speed: gas.sound_speed,
    }
}

// Independent entire-function series: J1(x) = sum (-1)^m (x/2)^(2m+1)/(m!(m+1)!).
fn directivity(x: f64) -> f64 {
    if x == 0.0 {
        return 1.0;
    }
    let mut term = 0.5 * x;
    let mut sum = term;
    for m in 1..64 {
        term *= -0.25 * x * x / (m * (m + 1)) as f64;
        sum += term;
        if term.abs() < 1e-17 {
            break;
        }
    }
    2.0 * sum / x
}

fn config(position_m: [f64; 3], rings: usize, angles: usize, band: f64) -> CircularOutletReceiver {
    CircularOutletReceiver {
        position_m,
        radial_rings: rings,
        angular_points: angles,
        maximum_frequency_hz: band,
    }
}

// Least-squares quadrature fits BOTH sine and cosine, so no whole-cycle sampling
// assumption or half-sample/filter phase adjustment can bias the measured gain.
fn amplitude(mut mic: BaffledPressure, rate: u32, frequency: f64, flow: f64) -> f64 {
    let (a, b) = harmonics(&mut mic, rate, frequency, flow);
    det::sqrt(a * a + b * b)
}

fn harmonics(mic: &mut BaffledPressure, rate: u32, frequency: f64, flow: f64) -> (f64, f64) {
    let omega = TAU * frequency / f64::from(rate);
    let warmup = mic.delay_samples.1 + 256;
    let (mut ss, mut cc, mut sc, mut sp, mut cp) = (0.0, 0.0, 0.0, 0.0, 0.0);
    for n in 0..warmup + 2048 {
        let (s, c) = (det::sin(omega * n as f64), det::cos(omega * n as f64));
        let p = mic.step(&[flow * s]).unwrap();
        if n >= warmup {
            ss += s * s;
            cc += c * c;
            sc += s * c;
            sp += s * p;
            cp += c * p;
        }
    }
    let determinant = ss * cc - sc * sc;
    assert!(determinant > 0.0);
    let (a, b) = (
        (sp * cc - cp * sc) / determinant,
        (cp * ss - sp * sc) / determinant,
    );
    (a, b)
}

#[allow(clippy::too_many_arguments)]
fn log(case: &str, axis: &str, rung: f64, observed: f64, oracle: f64, tol: f64, start: Instant) {
    let err = (observed - oracle).abs();
    println!(
        "{{\"case\":\"{case}\",\"inputs\":\"analytic uniform baffled piston; SI; GasState dry-air USSA1976 293.15 K 101325 Pa\",\"axis\":\"{axis}\",\"rung\":{rung},\"observed\":{observed:e},\"oracle\":{oracle:e},\"err\":{err:e},\"tol\":{tol:e},\"verdict\":\"{}\",\"elapsed_s\":{}}}",
        if err <= tol { "pass" } else { "fail" },
        start.elapsed().as_secs_f64()
    );
    assert!(
        err <= tol,
        "{case} {axis}={rung}: {observed:e} versus {oracle:e}, error {err:e}"
    );
}

#[test]
fn baffled_pressure_on_axis_piston() {
    let start = Instant::now();
    let medium = medium();
    // 25 mm keeps z/a=2 inside the owner's >=5 cm receiver domain. The
    // dimensionless oracle is unchanged; 7 mm is used in the far-field tests.
    let (radius, rate, frequency, velocity) = (0.025, 192_000, 800.0, 0.01);
    let area = PI * radius * radius;
    for ratio in [2.0, 10.0, 100.0] {
        let z = radius * ratio;
        let delta = radius * radius / (det::sqrt(z * z + radius * radius) + z);
        let exact = 2.0
            * medium.density
            * medium.sound_speed
            * velocity
            * det::sin(PI * frequency * delta / medium.sound_speed).abs();
        let mic = BaffledPressure::circular_outlet(
            radius,
            rate,
            config([0.0, 0.0, z], 128, 128, frequency),
            medium,
        )
        .unwrap();
        let measured = amplitude(mic, rate, frequency, area * velocity);
        log(
            "baffled_pressure_on_axis_piston",
            "z_over_a",
            ratio,
            measured / exact,
            1.0,
            1e-3,
            start,
        );
    }
}

#[test]
fn baffled_pressure_far_field_directivity() {
    let start = Instant::now();
    let medium = medium();
    let (radius, rate, velocity) = (0.007, 192_000, 0.01);
    let distance = 2000.0 * radius;
    for ka in [0.5, 2.0, 5.0] {
        let frequency = ka * medium.sound_speed / (TAU * radius);
        let axial = medium.density * medium.sound_speed * ka * radius * velocity / (2.0 * distance);
        for degrees in [0.0, 30.0, 60.0, 85.0] {
            let theta = degrees * PI / 180.0;
            let mut mic = BaffledPressure::circular_outlet_multirate_differentiating(
                radius,
                rate,
                rate,
                config(
                    [distance * det::sin(theta), 0.0, distance * det::cos(theta)],
                    128,
                    512,
                    frequency,
                ),
                medium,
            )
            .unwrap();
            let phase = TAU
                * frequency
                * (distance / medium.sound_speed + mic.filter_delay_samples() / f64::from(rate));
            let (sine, cosine) =
                harmonics(&mut mic, rate, frequency, PI * radius * radius * velocity);
            let measured = (cosine * det::cos(phase) + sine * det::sin(phase)) / axial;
            let quadrature = (sine * det::cos(phase) - cosine * det::sin(phase)) / axial;
            let exact = directivity(ka * det::sin(theta));
            log(
                &format!("baffled_pressure_far_field_directivity_ka_{ka}"),
                "angle_deg",
                degrees,
                measured,
                exact,
                1e-3,
                start,
            );
            log(
                &format!("baffled_pressure_far_field_phase_ka_{ka}"),
                "angle_deg",
                degrees,
                quadrature,
                0.0,
                1e-3,
                start,
            );
        }
    }
}

#[test]
fn baffled_pressure_ring_angle_band_ladder() {
    let start = Instant::now();
    let medium = medium();
    let (radius, rate, frequency) = (0.025, 192_000, 1200.0);
    let z = 2.0 * radius;
    let delta = radius * radius / (det::sqrt(z * z + radius * radius) + z);
    let exact = 2.0
        * medium.density
        * medium.sound_speed
        * det::sin(PI * frequency * delta / medium.sound_speed).abs();
    let mut radial_errors = Vec::new();
    for rings in [8, 32, 128] {
        let mic = BaffledPressure::circular_outlet_multirate_differentiating(
            radius,
            rate,
            rate,
            config([0.0, 0.0, z], rings, 128, frequency),
            medium,
        )
        .unwrap();
        let observed = amplitude(mic, rate, frequency, PI * radius * radius) / exact;
        radial_errors.push((observed - 1.0).abs());
        log(
            "baffled_pressure_ring_angle_band_ladder",
            "radial_rings",
            rings as f64,
            observed,
            1.0,
            1e-3,
            start,
        );
    }
    assert!(
        radial_errors[0] > radial_errors[1] && radial_errors[1] > radial_errors[2],
        "radial axis must converge, errors {radial_errors:?}"
    );

    // Angle refinement uses an independent, dense angular integral at fixed
    // radial nodes. This isolates azimuthal error from the radial/time floors.
    // Work at ka=5, where 8 points do not yet resolve the angular interference.
    let radius = 0.007;
    let frequency = 5.0 * medium.sound_speed / (TAU * radius);
    let position = [0.07, 0.0, 0.05];
    let spatial = |angles: usize| {
        let mut re = 0.0;
        let mut im = 0.0;
        for ring in 0..128 {
            let r = radius * det::sqrt((ring as f64 + 0.5) / 128.0);
            for angle in 0..angles {
                let phi = TAU * (angle as f64 + 0.5) / angles as f64;
                let distance = det::sqrt(
                    (position[0] - r * det::cos(phi)).powi(2)
                        + (r * det::sin(phi)).powi(2)
                        + position[2] * position[2],
                );
                let phase = TAU * frequency * distance / medium.sound_speed;
                re += det::cos(phase) / distance / (128 * angles) as f64;
                im += det::sin(phase) / distance / (128 * angles) as f64;
            }
        }
        (re, im)
    };
    let reference = spatial(2048);
    let norm = det::sqrt(reference.0 * reference.0 + reference.1 * reference.1);
    let mut angular_errors = Vec::new();
    for angles in [8, 16, 32] {
        // circular_surface's conservative pi/4 cell admission is itself checked;
        // lower declared bands admit these meshes, then the test deliberately
        // probes high-frequency angular content without claiming it resolved.
        let declared = 0.4 * medium.sound_speed / (TAU * radius);
        let mic = BaffledPressure::circular_outlet_multirate_differentiating(
            radius,
            rate,
            rate,
            config(position, 128, angles, declared),
            medium,
        )
        .unwrap();
        let expected = spatial(angles);
        let actual = amplitude(mic, rate, frequency, 1.0);
        let oracle = medium.density
            * frequency
            * det::sqrt(expected.0 * expected.0 + expected.1 * expected.1);
        log(
            "baffled_pressure_ring_angle_band_ladder",
            "angular_points",
            angles as f64,
            actual / oracle,
            1.0,
            1e-3,
            start,
        );
        let error =
            det::sqrt((expected.0 - reference.0).powi(2) + (expected.1 - reference.1).powi(2))
                / norm;
        println!(
            "{{\"case\":\"baffled_pressure_angular_convergence\",\"inputs\":\"analytic Rayleigh angular integral at fixed 128 radial midpoints; reference 2048 azimuths; ka=5; a=0.007m; x=0.07m; z=0.05m\",\"axis\":\"angular_points\",\"rung\":{angles},\"observed\":{error:e},\"oracle\":0,\"err\":{error:e},\"tol\":1e-3,\"verdict\":\"{}\",\"elapsed_s\":{}}}",
            if error <= 1e-3 { "pass" } else { "fail" },
            start.elapsed().as_secs_f64()
        );
        angular_errors.push(error);
    }
    // By 32 points this analytic angular integrand is already at roundoff.
    assert!(
        angular_errors[0] > 1e-8
            && angular_errors[1] < angular_errors[0] * 0.01
            && angular_errors[2] < 1e-10,
        "angular axis must be live then converge: {angular_errors:?}"
    );

    // maximum_frequency_hz is admission metadata, NOT a tunable filter. Keep
    // all other inputs fixed and prove admitted kernels bit-identical. Its
    // liveness is the phase-span refusal at higher bands, never a fake audio gain.
    let mut previous = None;
    for band in [1000.0, 2000.0, 4000.0] {
        let mic = BaffledPressure::circular_outlet(
            0.007,
            48_000,
            config([0.0, 0.0, 0.2], 8, 32, band),
            medium,
        )
        .unwrap();
        if let Some(kernel) = previous {
            assert_eq!(mic.kernels, kernel);
        }
        previous = Some(mic.kernels);
        log(
            "baffled_pressure_ring_angle_band_ladder",
            "maximum_frequency_hz",
            band,
            0.0,
            0.0,
            0.0,
            start,
        );
    }
    assert!(
        BaffledPressure::circular_outlet(
            0.025,
            48_000,
            config([0.0, 0.0, 0.2], 1, 8, 4000.0),
            medium
        )
        .is_err()
    );
}

#[test]
fn circular_outlet_band_cap_refusal() {
    let medium = medium();
    // Small physical disk isolates the temporal cap from spatial phase admission.
    for (differentiating, fraction) in [(false, 0.1), (true, 0.45)] {
        for band in [fraction * 48_000.0, fraction * 48_000.0 + 1.0] {
            let receiver = config([0.0, 0.0, 0.2], 128, 512, band);
            let result = if differentiating {
                BaffledPressure::circular_outlet_multirate_differentiating(
                    0.001, 48_000, 48_000, receiver, medium,
                )
            } else {
                BaffledPressure::circular_outlet(0.001, 48_000, receiver, medium)
            };
            assert_eq!(result.is_ok(), band == fraction * 48_000.0);
        }
    }
}

#[test]
fn observer_bitwise_repeat() {
    let receiver = config([0.06, 0.0, 0.2], 16, 64, 2000.0);
    for differentiating in [false, true] {
        let build = || {
            if differentiating {
                BaffledPressure::circular_outlet_multirate_differentiating(
                    0.007,
                    48_000,
                    48_000,
                    receiver,
                    medium(),
                )
                .unwrap()
            } else {
                BaffledPressure::circular_outlet(0.007, 48_000, receiver, medium()).unwrap()
            }
        };
        let (mut a, mut b) = (build(), build());
        let mut peak = 0.0_f64;
        for n in 0..1024 {
            let q = 1e-5 * det::sin(n as f64 * 0.13);
            let p = a.step(&[q]).unwrap();
            assert_eq!(p.to_bits(), b.step(&[q]).unwrap().to_bits());
            peak = peak.max(p.abs());
        }
        assert!(peak > 0.0);
    }
}
