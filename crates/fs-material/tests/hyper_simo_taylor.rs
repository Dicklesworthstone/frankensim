//! G0/G1: collapse-resistant compression, its pressure dual and finite-step work.
//! Independent scalar derivatives follow U=K/4 (J^2-1-2 ln J), Hartmann &
//! Neff (2003), doi:10.1016/S0020-7683(03)00086-6, Table 4.

use fs_ad::{Dual, Dual64};
use fs_material::hyper::discrete_gradient::bulk::{DiscreteBulkPoint, DiscreteBulkSettings};
use fs_material::hyper::mixed::{
    MixedFiniteStrainPotential, VolumetricTreatment, pressure_potential,
};
use fs_material::hyper::{Hyperelastic, HyperelasticModel, OgdenTerm, Volumetric};

fn close(actual: f64, expected: f64, scale: f64, tolerance: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= tolerance * scale,
        "actual={actual:.17e}, expected={expected:.17e}, scale={scale:e}"
    );
}

fn pressure(kappa: f64, j: f64) -> f64 {
    0.5 * kappa * (j - 1.0) * (1.0 + 1.0 / j)
}

fn energy(kappa: f64, j: f64) -> f64 {
    0.25 * kappa * ((j - 1.0) * (j + 1.0) - 2.0 * j.ln())
}

#[test]
fn g1_bulk_energy_stress_and_curvature_survive_weak_and_finite_compression() {
    for kappa in [1.0, 8e5] {
        let law = Volumetric::SimoTaylor { kappa };
        for j in [
            1e-3_f64,
            0.125,
            0.874999999999,
            0.875,
            0.875000000001,
            1.0 - 1e-12,
            1.0,
            1.0 + 1e-12,
            1.124999999999,
            1.125,
            1.125000000001,
            3.0,
            40.0,
        ] {
            let variable = Dual::<Dual64<1>, 1>::variable(Dual64::variable(j, 0), 0);
            let response = law.energy(variable);
            let expected_pressure = pressure(kappa, j);
            let curvature = 0.5 * kappa * (1.0 + 1.0 / (j * j));
            close(
                response.eps[0].re,
                expected_pressure,
                expected_pressure.abs().max(kappa * 1e-14),
                2e-14,
            );
            close(response.eps[0].eps[0], curvature, curvature, 2e-14);
            assert!(response.re.re >= 0.0 && curvature > 0.0);
            let x = j - 1.0;
            let expected_energy = if x.abs() < 1e-5 {
                0.5 * kappa * x * x * (1.0 - x / 3.0 + x * x / 4.0)
            } else {
                energy(kappa, j)
            };
            close(
                response.re.re,
                expected_energy,
                expected_energy.abs().max(kappa * 1e-30),
                3e-14,
            );
        }
        assert!(law.energy(1e-9) > law.energy(1e-6));
        assert!(law.energy(1e-6) > law.energy(1e-3));
    }
}

#[test]
fn g0_pressure_inverse_has_one_positive_volume_through_compression_and_expansion() {
    let kappa = 8e5;
    let law = Volumetric::SimoTaylor { kappa };
    for ratio in [
        -1e150_f64, -1e8, -10.0, -0.02, -1e-18, 0.0, 1e-18, 0.02, 10.0, 1e8, 1e150,
    ] {
        let p = kappa * ratio;
        let result = pressure_potential(law, p).unwrap();
        assert!(result.volume_ratio > 0.0 && result.compliance_per_pa > 0.0);
        assert!(result.bulk_stored_pa >= 0.0 && result.complementary_pa >= 0.0);
        if ratio.abs() >= 1e-12 {
            close(pressure(kappa, result.volume_ratio), p, p.abs(), 5e-14);
        } else if ratio != 0.0 {
            // Pressure still carries sub-ulp reference strain even if J rounds to 1.
            close(result.volume_increment, ratio, ratio.abs(), 2e-14);
            close(result.bulk_stored_pa, 0.5 * p * ratio, p * ratio, 2e-14);
        }
        let j = result.volume_ratio;
        let expected_compliance = if j <= 1.0 {
            (2.0 * j / kappa) * (j / (1.0 + j * j))
        } else {
            (2.0 / kappa) / (1.0 + (1.0 / j).powi(2))
        };
        close(
            result.compliance_per_pa,
            expected_compliance,
            expected_compliance,
            5e-14,
        );
    }
    // p/kappa would overflow, but this compressed volume and compliance exist.
    let extreme = pressure_potential(Volumetric::SimoTaylor { kappa: 1e-300 }, -1e10).unwrap();
    close(extreme.volume_ratio, 5e-311, 5e-311, 1e-12);
    assert!(extreme.compliance_per_pa > 0.0);
}

#[test]
fn g1_hydrostatic_piola_and_mixed_condensation_retain_the_same_bulk_energy() {
    let kappa = 8e5;
    let volume = Volumetric::SimoTaylor { kappa };
    let models = [
        HyperelasticModel::Gent {
            mu: 3e4,
            locking_invariant: 20.0,
            volumetric: volume,
        },
        HyperelasticModel::Ogden {
            terms: vec![OgdenTerm {
                mu_pa: 3e4,
                alpha: 2.0,
            }],
            volumetric: volume,
        },
    ];
    for model in models {
        let law = Hyperelastic::new(model, 8.0).unwrap();
        assert_eq!(law.admissibility().polyconvex, Some(true));
        for stretch in [0.2_f64, 0.9, 1.0, 1.8] {
            let f = [stretch, 0.0, 0.0, 0.0, stretch, 0.0, 0.0, 0.0, stretch];
            let j = stretch.powi(3);
            let p = pressure(kappa, j);
            let native = law.piola(&f).unwrap();
            let tangent = law.tangent(&f).unwrap();
            let mixed = law
                .mixed_response(&f, p, VolumetricTreatment::OriginalMaterial)
                .unwrap();
            close(law.energy(&f), energy(kappa, j), kappa, 3e-13);
            close(mixed.state.pressure_residual, 0.0, j, 2e-13);
            for i in 0..9 {
                close(
                    native[i],
                    if i % 4 == 0 {
                        p * stretch * stretch
                    } else {
                        0.0
                    },
                    kappa.max(p.abs()),
                    3e-13,
                );
                close(
                    mixed.state.first_piola_pa[i],
                    native[i],
                    kappa.max(p.abs()),
                    3e-13,
                );
                for col in 0..9 {
                    let condensed = mixed.tangent_pa[i][col]
                        - mixed.pressure_coupling[i] * mixed.pressure_coupling[col]
                            / mixed.pressure_tangent_per_pa;
                    close(condensed, tangent[i][col], kappa.max(p.abs()), 2e-12);
                }
            }
        }
    }
}

#[test]
fn g1_finite_compression_work_is_reversible_across_the_full_positive_volume_branch() {
    let kappa = 8e5;
    let law = Hyperelastic::new(
        HyperelasticModel::Gent {
            mu: 3e4,
            locking_invariant: 20.0,
            volumetric: Volumetric::SimoTaylor { kappa },
        },
        8.0,
    )
    .unwrap();
    let settings = DiscreteBulkSettings {
        minimum_volume_ratio: 0.01,
        maximum_volume_ratio: 6.0,
        maximum_iterations: 96,
        absolute_pressure_tolerance_pa: 2e-14 * kappa,
        relative_pressure_tolerance: 2e-14,
        absolute_work_tolerance_pa: 2e-14 * kappa,
        relative_work_tolerance: 2e-14,
    };
    for (a, b) in [
        (1.0, 0.03),
        (0.03, 1.0),
        (0.94, 1.07),
        (1.0, 5.0),
        (4.0, 3.0),
    ] {
        let point = DiscreteBulkPoint::try_new(&law, a, settings).unwrap();
        let change = energy(kappa, b) - energy(kappa, a);
        let interval_pressure = change / (b - a);
        let response = point.response(interval_pressure).unwrap();
        close(response.volume_ratio_after, b, b, 3e-12);
        close(response.work_pa, change, kappa.max(change.abs()), 3e-12);
        close(
            response.energy_change_pa,
            change,
            kappa.max(change.abs()),
            3e-12,
        );
        close(
            response.endpoint_pressure_pa,
            pressure(kappa, b),
            kappa.max(interval_pressure.abs()),
            3e-12,
        );
        let reverse = point
            .rebind(b)
            .unwrap()
            .response(interval_pressure)
            .unwrap();
        close(reverse.volume_ratio_after, a, a, 3e-12);
        close(
            reverse.work_pa + response.work_pa,
            0.0,
            kappa.max(change.abs()),
            5e-12,
        );
        assert_eq!(point.accepted_volume_ratio(), a);
    }
}
