//! G0/G1: original finite-bulk secant inversion, work, compliance and branch admission.

use fs_material::hyper::discrete_gradient::bulk::{
    DiscreteBulkPoint, DiscreteBulkSettings, MAX_BULK_SECANT_ITERATIONS,
};
use fs_material::hyper::{Hyperelastic, HyperelasticModel, OgdenTerm, Volumetric};

fn settings(kappa: f64) -> DiscreteBulkSettings {
    DiscreteBulkSettings {
        minimum_volume_ratio: 0.2,
        maximum_volume_ratio: 2.6,
        maximum_iterations: 96,
        absolute_pressure_tolerance_pa: 2e-14 * kappa,
        relative_pressure_tolerance: 2e-14,
        absolute_work_tolerance_pa: 2e-14 * kappa,
        relative_work_tolerance: 2e-14,
    }
}

fn law(volume: Volumetric) -> Hyperelastic {
    let model = match volume {
        Volumetric::Quadratic { kappa } => HyperelasticModel::MooneyRivlin {
            c10: 4.0,
            c01: 2.0,
            kappa,
        },
        Volumetric::LogSquared { .. } => HyperelasticModel::Ogden {
            terms: vec![OgdenTerm {
                mu_pa: 6.0,
                alpha: 4.0,
            }],
            volumetric: volume,
        },
    };
    Hyperelastic::new(model, 4.0).unwrap()
}

fn close(actual: f64, expected: f64, allowance: f64) {
    assert!(
        actual.is_finite() && expected.is_finite() && (actual - expected).abs() <= allowance,
        "actual={actual:.17e}, expected={expected:.17e}, allowance={allowance:e}"
    );
}

// Independent scalar oracle uses platform logarithms and direct original U.
fn energy(volume: Volumetric, j: f64) -> f64 {
    let strain = match volume {
        Volumetric::Quadratic { .. } => j - 1.0,
        Volumetric::LogSquared { .. } => j.ln(),
    };
    0.5 * volume.kappa() * strain * strain
}

fn endpoint_pressure(volume: Volumetric, j: f64) -> f64 {
    volume.kappa()
        * match volume {
            Volumetric::Quadratic { .. } => j - 1.0,
            Volumetric::LogSquared { .. } => j.ln() / j,
        }
}

fn secant(volume: Volumetric, a: f64, b: f64) -> f64 {
    if (b - a).abs() <= 1e-6 * a {
        // Independent Simpson average of U', sufficient to fourth order in
        // this small interval, including the exact equality limit.
        (endpoint_pressure(volume, a)
            + 4.0 * endpoint_pressure(volume, 0.5 * a + 0.5 * b)
            + endpoint_pressure(volume, b))
            / 6.0
    } else {
        (energy(volume, b) - energy(volume, a)) / (b - a)
    }
}

fn integrate_pressure(volume: Volumetric, a: f64, b: f64, intervals: usize) -> f64 {
    let mut sum = endpoint_pressure(volume, a) + endpoint_pressure(volume, b);
    for i in 1..intervals {
        let j = a + (b - a) * (i as f64 / intervals as f64);
        sum += (if i % 2 == 0 { 2.0 } else { 4.0 }) * endpoint_pressure(volume, j);
    }
    sum / (3 * intervals) as f64
}

#[test]
fn g0_quadratic_original_storage_and_twice_endpoint_compliance_are_exact() {
    for kappa in [0.7, 60.0, 8e5, 3e12] {
        let volume = Volumetric::Quadratic { kappa };
        let material = law(volume);
        for a in [0.31, 1.0, 1.27, 2.1] {
            let point = DiscreteBulkPoint::try_new(&material, a, settings(kappa)).unwrap();
            for b in [0.25, 0.83, 1.0, 1.97, 2.55] {
                let pressure = kappa * (0.5 * (a - 1.0) + 0.5 * (b - 1.0));
                let response = point.response(pressure).unwrap();
                close(response.volume_ratio_after, b, 2e-15);
                close(response.energy_after_pa, energy(volume, b), 8e-15 * kappa);
                close(
                    response.energy_change_pa,
                    energy(volume, b) - energy(volume, a),
                    8e-15 * kappa,
                );
                close(
                    response.work_pa,
                    energy(volume, b) - energy(volume, a),
                    8e-15 * kappa,
                );
                close(
                    response.endpoint_pressure_pa,
                    endpoint_pressure(volume, b),
                    8e-15 * kappa,
                );
                close(response.compliance_per_pa, 2.0 / kappa, 2e-15 / kappa);
                assert_eq!(response.iterations, 0);
                assert_eq!(point.accepted_volume_ratio(), a);
            }
        }
    }
}

#[test]
fn g1_log_bulk_recovers_original_energy_and_independent_pressure_integral() {
    for kappa in [7.0, 9e5] {
        let volume = Volumetric::LogSquared { kappa };
        let material = law(volume);
        for (a, b) in [
            (0.31, 2.45),
            (2.4, 0.25),
            (0.73, 1.37),
            (1.73, 0.8),
            (1.0, 1.0),
        ] {
            let point = DiscreteBulkPoint::try_new(&material, a, settings(kappa)).unwrap();
            let pressure = secant(volume, a, b);
            let coarse = integrate_pressure(volume, a, b, 4096);
            let fine = integrate_pressure(volume, a, b, 8192);
            close(coarse, fine, 5e-11 * kappa);
            close(pressure, fine, 4e-12 * kappa);
            let response = point.response(pressure).unwrap();
            close(response.volume_ratio_after, b, 3e-12);
            close(response.energy_after_pa, energy(volume, b), 3e-13 * kappa);
            close(
                response.energy_change_pa,
                energy(volume, b) - energy(volume, a),
                3e-13 * kappa,
            );
            close(
                response.work_pa,
                energy(volume, b) - energy(volume, a),
                3e-13 * kappa,
            );
            close(
                response.endpoint_pressure_pa,
                endpoint_pressure(volume, b),
                3e-12 * kappa,
            );
            assert!(response.compliance_per_pa > 0.0);
            assert!(response.iterations <= settings(kappa).maximum_iterations);
            let reverse = point
                .rebind(response.volume_ratio_after)
                .unwrap()
                .response(pressure)
                .unwrap();
            close(reverse.volume_ratio_after, a, 3e-12);
            close(reverse.work_pa + response.work_pa, 0.0, 5e-13 * kappa);
            assert_eq!(point.accepted_volume_ratio(), a);
        }
    }
}

#[test]
fn g0_consistent_compliance_matches_pressure_finite_differences_and_small_pair_limit() {
    let mut worst = 0.0_f64;
    for volume in [
        Volumetric::Quadratic { kappa: 63.0 },
        Volumetric::LogSquared { kappa: 8e5 },
    ] {
        let kappa = volume.kappa();
        let material = law(volume);
        for (a, b) in [(0.42, 0.91), (1.16, 0.73), (0.83, 1.41), (2.2, 2.4)] {
            let point = DiscreteBulkPoint::try_new(&material, a, settings(kappa)).unwrap();
            let pressure = secant(volume, a, b);
            let response = point.response(pressure).unwrap();
            let h = 2e-6 * kappa;
            let difference = (point.response(pressure + h).unwrap().volume_ratio_after
                - point.response(pressure - h).unwrap().volume_ratio_after)
                / (2.0 * h);
            worst = worst.max((difference / response.compliance_per_pa - 1.0).abs());
            close(
                difference,
                response.compliance_per_pa,
                8e-8 * response.compliance_per_pa,
            );
        }
        for a in [0.83_f64, 1.0, 1.9, 2.5] {
            let point = DiscreteBulkPoint::try_new(&material, a, settings(kappa)).unwrap();
            let equal = point.response(point.endpoint_pressure_pa()).unwrap();
            let compliance = match volume {
                Volumetric::Quadratic { .. } => 2.0 / kappa,
                Volumetric::LogSquared { .. } => 2.0 * a * a / (kappa * (1.0 - a.ln())),
            };
            close(equal.volume_ratio_after, a, 5e-16);
            close(equal.compliance_per_pa, compliance, 8e-15 * compliance);
            for delta in [-1e-8, -1e-12, 1e-12, 1e-8] {
                let b = a + delta;
                let response = point.response(secant(volume, a, b)).unwrap();
                close(response.volume_ratio_after, b, 8e-13);
                let stable_change = match volume {
                    Volumetric::Quadratic { .. } => {
                        kappa * (0.5 * (a - 1.0) + 0.5 * (b - 1.0)) * (b - a)
                    }
                    Volumetric::LogSquared { .. } => {
                        let dl = ((b - a) / a).ln_1p();
                        kappa * (a.ln() + 0.5 * dl) * dl
                    }
                };
                close(response.energy_change_pa, stable_change, 3e-13 * kappa);
            }
        }
    }
    println!("original bulk compliance worst relative finite-difference discrepancy={worst:e}");
}

#[test]
fn g0_original_endpoint_pressure_initializes_compatible_bulk_and_zero_pressure_reverses_log_strain()
{
    for volume in [
        Volumetric::Quadratic { kappa: 60.0 },
        Volumetric::LogSquared { kappa: 700.0 },
    ] {
        let material = law(volume);
        for j in [0.25, 0.8, 1.0, 1.37, 2.5] {
            let pressure = endpoint_pressure(volume, j);
            let point = DiscreteBulkPoint::from_endpoint_pressure(
                &material,
                pressure,
                settings(volume.kappa()),
            )
            .unwrap();
            close(point.accepted_volume_ratio(), j, 3e-14);
            close(
                point.energy_before_pa(),
                energy(volume, j),
                2e-14 * volume.kappa(),
            );
            close(
                point.response(pressure).unwrap().volume_ratio_after,
                j,
                3e-13,
            );
        }
    }
    let volume = Volumetric::LogSquared { kappa: 71.0 };
    let material = law(volume);
    let point = DiscreteBulkPoint::try_new(&material, 0.5, settings(71.0)).unwrap();
    let response = point.response(0.0).unwrap();
    close(response.volume_ratio_after, 2.0, 2e-13);
    close(response.energy_change_pa, 0.0, 2e-14 * 71.0);
    assert_eq!(response.work_pa, 0.0);
    assert!(response.endpoint_pressure_pa > 0.0);
}

#[test]
fn g0_invalid_domains_exhausted_scalar_work_and_insufficient_work_allowances_refuse() {
    let volume = Volumetric::LogSquared { kappa: 60.0 };
    let material = law(volume);
    let controls = settings(60.0);
    let point = DiscreteBulkPoint::try_new(&material, 0.73, controls).unwrap();
    for pressure in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -1e6, 1e6] {
        assert!(point.response(pressure).is_err());
    }
    for pressure in [f64::NAN, 60.0, endpoint_pressure(volume, 0.15)] {
        assert!(DiscreteBulkPoint::from_endpoint_pressure(&material, pressure, controls).is_err());
    }
    assert!(point.rebind(0.1).is_err());
    assert!(point.rebind(2.7).is_err());
    assert_eq!(point.accepted_volume_ratio(), 0.73);
    for controls in [
        DiscreteBulkSettings {
            minimum_volume_ratio: 0.0,
            ..controls
        },
        DiscreteBulkSettings {
            maximum_volume_ratio: core::f64::consts::E,
            ..controls
        },
        DiscreteBulkSettings {
            maximum_volume_ratio: 0.1,
            ..controls
        },
        DiscreteBulkSettings {
            maximum_iterations: 0,
            ..controls
        },
        DiscreteBulkSettings {
            maximum_iterations: MAX_BULK_SECANT_ITERATIONS + 1,
            ..controls
        },
        DiscreteBulkSettings {
            absolute_work_tolerance_pa: f64::NAN,
            ..controls
        },
        DiscreteBulkSettings {
            relative_pressure_tolerance: -1.0,
            ..controls
        },
    ] {
        assert!(DiscreteBulkPoint::try_new(&material, 0.73, controls).is_err());
    }
    let pressure = secant(volume, 0.73, 1.57);
    let short = DiscreteBulkPoint::try_new(
        &material,
        0.73,
        DiscreteBulkSettings {
            maximum_iterations: 1,
            ..controls
        },
    )
    .unwrap();
    assert!(
        short
            .response(pressure)
            .unwrap_err()
            .to_string()
            .contains("exhausted")
    );
    let work_limited = DiscreteBulkPoint::try_new(
        &material,
        0.73,
        DiscreteBulkSettings {
            maximum_iterations: 1,
            absolute_pressure_tolerance_pa: 1e6,
            absolute_work_tolerance_pa: 0.0,
            relative_work_tolerance: 0.0,
            ..controls
        },
    )
    .unwrap();
    assert!(work_limited.response(pressure).is_err());
    let coupled = Hyperelastic::new(
        HyperelasticModel::NeoHookean {
            mu: 6.0,
            lambda: 60.0,
        },
        4.0,
    )
    .unwrap();
    assert!(DiscreteBulkPoint::try_new(&coupled, 1.0, controls).is_err());
    let mut invalid = law(Volumetric::Quadratic { kappa: 60.0 });
    invalid.stretch_limit = f64::NAN;
    assert!(DiscreteBulkPoint::try_new(&invalid, 1.0, controls).is_err());
}
