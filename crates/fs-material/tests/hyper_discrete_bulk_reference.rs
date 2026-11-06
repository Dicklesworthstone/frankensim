//! Held-reference finite bulk retains total volume and original physical work.

use fs_material::hyper::discrete_gradient::bulk::{DiscreteBulkPoint, DiscreteBulkSettings};
use fs_material::hyper::reference::ReferenceDeformation;
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
    Hyperelastic::new(
        match volume {
            Volumetric::Quadratic { kappa } => HyperelasticModel::MooneyRivlin {
                c10: 4.0,
                c01: 2.0,
                kappa,
            },
            Volumetric::LogSquared { .. } | Volumetric::SimoTaylor { .. } => {
                HyperelasticModel::Ogden {
                    terms: vec![OgdenTerm {
                        mu_pa: 6.0,
                        alpha: 4.0,
                    }],
                    volumetric: volume,
                }
            }
        },
        4.0,
    )
    .unwrap()
}

fn close(a: f64, b: f64, allowance: f64) {
    assert!(
        a.is_finite() && b.is_finite() && (a - b).abs() <= allowance,
        "{a:.17e} versus {b:.17e}, allowance={allowance:e}"
    );
}

// Independent endpoint oracle: platform logarithms, not the production secant.
fn energy(volume: Volumetric, je: f64) -> f64 {
    if matches!(volume, Volumetric::SimoTaylor { .. }) {
        return 0.25 * volume.kappa() * ((je - 1.0) * (je + 1.0) - 2.0 * je.ln());
    }
    let strain = match volume {
        Volumetric::Quadratic { .. } => je - 1.0,
        Volumetric::LogSquared { .. } => je.ln(),
        Volumetric::SimoTaylor { .. } => unreachable!(),
    };
    0.5 * volume.kappa() * strain * strain
}

fn pressure(volume: Volumetric, je: f64) -> f64 {
    volume.kappa()
        * match volume {
            Volumetric::Quadratic { .. } => je - 1.0,
            Volumetric::LogSquared { .. } => je.ln() / je,
            Volumetric::SimoTaylor { .. } => 0.5 * (je - 1.0) * (1.0 + 1.0 / je),
        }
}

#[test]
fn g0_reference_and_modulus_work_resolves_second_order_change_under_large_preload() {
    // At retained J=2, U0=k/2. Prescribe g1=1+h and k1=k(1+3h).
    // Exact rational reduction gives
    // U1-U0 = k h^2 (-5+3h)/(2(1+h)). The O(h) source terms cancel;
    // neither rounded absolute endpoint can retain this physical change.
    let kappa = 2.0_f64.powi(80);
    let original = law(Volumetric::Quadratic { kappa });
    let before = DiscreteBulkPoint::try_new(&original, 2.0, settings(kappa)).unwrap();
    for exponent in [35, 40, 45] {
        let h = 2.0_f64.powi(-exponent);
        let next_kappa = kappa * (1.0 + 3.0 * h);
        let changed = law(Volumetric::Quadratic { kappa: next_kappa });
        let reference =
            ReferenceDeformation::try_new([1.0 + h, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0])
                .unwrap();
        let after = DiscreteBulkPoint::try_new_with_reference(
            &changed,
            &reference,
            2.0,
            settings(next_kappa),
        )
        .unwrap();
        let expected = 0.5 * kappa * h * h * (-5.0 + 3.0 * h) / (1.0 + h);
        assert_eq!(after.energy_before_pa(), before.energy_before_pa());
        close(
            before.energy_increment_to(&after).unwrap(),
            expected,
            2e-12 * expected.abs(),
        );
        close(
            after.energy_increment_to(&before).unwrap(),
            -expected,
            2e-12 * expected.abs(),
        );
        assert_eq!(before.energy_increment_to(&before).unwrap(), 0.0);
        assert_eq!(
            before.accepted_volume_ratio(),
            after.accepted_volume_ratio()
        );
    }
}

#[test]
fn g0_all_bulk_sources_retain_small_natural_volume_work_at_fixed_elastic_strain() {
    let h = 2.0_f64.powi(-48);
    let reference =
        ReferenceDeformation::try_new([1.0 + h, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]).unwrap();
    for volume in [
        Volumetric::Quadratic { kappa: 1e24 },
        Volumetric::LogSquared { kappa: 1e24 },
        Volumetric::SimoTaylor { kappa: 1e24 },
    ] {
        let material = law(volume);
        let before = DiscreteBulkPoint::try_new(&material, 2.0, settings(volume.kappa())).unwrap();
        let after = DiscreteBulkPoint::try_new_with_reference(
            &material,
            &reference,
            2.0 * (1.0 + h),
            settings(volume.kappa()),
        )
        .unwrap();
        // The elastic strain is unchanged, so only the relaxed volume scales
        // the independent original stored-energy density.
        let expected = h * energy(volume, 2.0);
        close(
            before.energy_increment_to(&after).unwrap(),
            expected,
            2e-14 * expected.abs(),
        );
        close(
            after.energy_increment_to(&before).unwrap(),
            -expected,
            2e-14 * expected.abs(),
        );
    }
}

#[test]
fn g1_reference_bulk_snapshot_changes_keep_each_original_constitutive_family() {
    let before_reference =
        ReferenceDeformation::try_new([1.2, 0.1, 0.0, 0.0, 1.1, 0.0, 0.0, 0.0, 0.9]).unwrap();
    let after_reference =
        ReferenceDeformation::try_new([0.9, 0.0, 0.0, 0.0, 1.3, 0.2, 0.0, 0.0, 1.1]).unwrap();
    let families = [
        Volumetric::Quadratic { kappa: 9e5 },
        Volumetric::LogSquared { kappa: 7e5 },
        Volumetric::SimoTaylor { kappa: 6e5 },
    ];
    for before_volume in families {
        for after_volume in families {
            let before_law = law(before_volume);
            let after_law = law(after_volume);
            let a = DiscreteBulkPoint::try_new_with_reference(
                &before_law,
                &before_reference,
                1.4,
                settings(before_volume.kappa()),
            )
            .unwrap();
            let b = DiscreteBulkPoint::try_new_with_reference(
                &after_law,
                &after_reference,
                0.8,
                settings(after_volume.kappa()),
            )
            .unwrap();
            let expected = after_reference.volume_ratio()
                * energy(after_volume, 0.8 / after_reference.volume_ratio())
                - before_reference.volume_ratio()
                    * energy(before_volume, 1.4 / before_reference.volume_ratio());
            close(a.energy_increment_to(&b).unwrap(), expected, 5e-10);
            close(b.energy_increment_to(&a).unwrap(), -expected, 5e-10);
        }
    }
}

#[test]
fn g0_opposing_reference_and_modulus_scales_do_not_overflow_finite_storage_change() {
    let large =
        ReferenceDeformation::try_new([1e200, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]).unwrap();
    let small =
        ReferenceDeformation::try_new([1e-200, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]).unwrap();
    for family in [0, 1, 2] {
        let volume = |kappa| match family {
            0 => Volumetric::Quadratic { kappa },
            1 => Volumetric::LogSquared { kappa },
            _ => Volumetric::SimoTaylor { kappa },
        };
        let soft = law(volume(1e-200));
        let stiff = law(volume(1e200));
        let before =
            DiscreteBulkPoint::try_new_with_reference(&soft, &large, 2e200, settings(1e-200))
                .unwrap();
        let after =
            DiscreteBulkPoint::try_new_with_reference(&stiff, &small, 2e-200, settings(1e200))
                .unwrap();
        close(before.energy_increment_to(&after).unwrap(), 0.0, 3e-16);
        close(after.energy_increment_to(&before).unwrap(), 0.0, 3e-16);
    }
}

#[test]
fn reference_bulk_secant_stores_total_j_and_original_energy_with_physical_pressure() {
    let g = ReferenceDeformation::try_new([1.7, 0.24, -0.17, 0.08, 1.3, 0.19, 0.14, -0.11, 1.5])
        .unwrap();
    let jg = g.volume_ratio();
    for volume in [
        Volumetric::Quadratic { kappa: 8e5 },
        Volumetric::LogSquared { kappa: 8e5 },
        Volumetric::SimoTaylor { kappa: 8e5 },
    ] {
        let material = law(volume);
        for (a, b) in [(0.41, 1.79), (1.63, 0.77), (0.93, 1.21)] {
            let accepted_total = jg * a;
            let point = DiscreteBulkPoint::try_new_with_reference(
                &material,
                &g,
                accepted_total,
                settings(volume.kappa()),
            )
            .unwrap();
            let p = (energy(volume, b) - energy(volume, a)) / (b - a);
            let response = point.response(p).unwrap();
            close(response.volume_ratio_after, jg * b, 3e-11);
            close(point.energy_before_pa(), jg * energy(volume, a), 2e-8);
            close(response.energy_after_pa, jg * energy(volume, b), 3e-7);
            close(
                response.energy_change_pa,
                jg * (energy(volume, b) - energy(volume, a)),
                3e-7,
            );
            close(
                response.work_pa,
                p * (response.volume_ratio_after - accepted_total),
                0.0,
            );
            close(response.work_pa, response.energy_change_pa, 3e-7);
            close(response.endpoint_pressure_pa, pressure(volume, b), 3e-7);
            let native = DiscreteBulkPoint::try_new(
                &material,
                accepted_total / jg,
                settings(volume.kappa()),
            )
            .unwrap()
            .response(p)
            .unwrap();
            close(
                response.compliance_per_pa,
                jg * native.compliance_per_pa,
                3e-16,
            );
            let dp = 1e-5 * volume.kappa();
            let plus = point.response(p + dp).unwrap().volume_ratio_after;
            let minus = point.response(p - dp).unwrap().volume_ratio_after;
            close(
                response.compliance_per_pa,
                (plus - minus) / (2.0 * dp),
                2e-8 * jg / volume.kappa(),
            );
            let allowance = settings(volume.kappa()).absolute_work_tolerance_pa
                + settings(volume.kappa()).relative_work_tolerance
                    * point
                        .energy_before_pa()
                        .abs()
                        .max(response.energy_after_pa.abs())
                        .max(response.work_pa.abs());
            assert!(response.work_defect_pa.abs() <= allowance);
            assert_eq!(
                point.accepted_volume_ratio().to_bits(),
                accepted_total.to_bits()
            );
        }
    }
}

#[test]
fn reference_bulk_initialization_rebind_and_elastic_convex_brackets_keep_the_reference() {
    let g = ReferenceDeformation::try_new([5.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 5.0]).unwrap();
    let jg = g.volume_ratio();
    for volume in [
        Volumetric::Quadratic { kappa: 60.0 },
        Volumetric::LogSquared { kappa: 60.0 },
    ] {
        let material = law(volume);
        let branch = DiscreteBulkSettings {
            minimum_volume_ratio: 0.8,
            maximum_volume_ratio: 1.2,
            ..settings(volume.kappa())
        };
        let point = DiscreteBulkPoint::from_endpoint_pressure_with_reference(
            &material,
            &g,
            pressure(volume, 1.07),
            branch,
        )
        .unwrap();
        close(point.accepted_volume_ratio(), jg * 1.07, 3e-12);
        close(point.endpoint_pressure_pa(), pressure(volume, 1.07), 3e-13);
        let retained = point.rebind(jg * 0.94).unwrap();
        assert_eq!(retained.accepted_volume_ratio(), jg * 0.94);
        close(
            retained.energy_before_pa(),
            jg * energy(volume, 0.94),
            3e-12,
        );
        close(
            retained.endpoint_pressure_pa(),
            pressure(volume, 0.94),
            3e-13,
        );
        assert!(
            retained.rebind(1.0).is_err(),
            "a total-volume value must not be read as elastic Je"
        );
        assert!(retained.rebind(jg * 0.79).is_err());
        assert!(retained.rebind(jg * 1.21).is_err());
        let natural = point.rebind(jg).unwrap().response(0.0).unwrap();
        assert_eq!(natural.volume_ratio_after, jg);
        assert_eq!(natural.energy_after_pa, 0.0);
        assert_eq!(natural.work_pa, 0.0);
        close(point.accepted_volume_ratio(), jg * 1.07, 3e-12);
    }
    let logarithmic = law(Volumetric::LogSquared { kappa: 60.0 });
    assert!(
        DiscreteBulkPoint::try_new_with_reference(
            &logarithmic,
            &g,
            jg,
            DiscreteBulkSettings {
                minimum_volume_ratio: 0.8,
                maximum_volume_ratio: 3.0,
                ..settings(60.0)
            },
        )
        .is_err(),
        "reference scaling cannot admit the nonconvex elastic log branch"
    );
}

#[test]
fn identity_bulk_is_exact_and_unrepresentable_original_branch_is_refused_atomically() {
    let identity = ReferenceDeformation::identity();
    for volume in [
        Volumetric::Quadratic { kappa: 8e5 },
        Volumetric::LogSquared { kappa: 8e5 },
    ] {
        let material = law(volume);
        let ordinary =
            DiscreteBulkPoint::try_new(&material, 1.07, settings(volume.kappa())).unwrap();
        let referenced = DiscreteBulkPoint::try_new_with_reference(
            &material,
            &identity,
            1.07,
            settings(volume.kappa()),
        )
        .unwrap();
        for p in [-0.08 * volume.kappa(), 0.0, 0.11 * volume.kappa()] {
            assert_eq!(
                ordinary.response(p).unwrap(),
                referenced.response(p).unwrap()
            );
            assert_eq!(
                ordinary.rebind(0.91).unwrap().response(p).unwrap(),
                referenced.rebind(0.91).unwrap().response(p).unwrap()
            );
        }
        let initial_p = pressure(volume, 1.07);
        let a = DiscreteBulkPoint::from_endpoint_pressure(
            &material,
            initial_p,
            settings(volume.kappa()),
        )
        .unwrap();
        let b = DiscreteBulkPoint::from_endpoint_pressure_with_reference(
            &material,
            &identity,
            initial_p,
            settings(volume.kappa()),
        )
        .unwrap();
        assert_eq!(
            a.accepted_volume_ratio().to_bits(),
            b.accepted_volume_ratio().to_bits()
        );
        assert_eq!(
            a.energy_before_pa().to_bits(),
            b.energy_before_pa().to_bits()
        );
    }
    let huge =
        ReferenceDeformation::try_new([1e100, 0.0, 0.0, 0.0, 1e100, 0.0, 0.0, 0.0, 1e100]).unwrap();
    let material = law(Volumetric::Quadratic { kappa: 60.0 });
    assert!(
        DiscreteBulkPoint::try_new_with_reference(
            &material,
            &huge,
            huge.volume_ratio(),
            DiscreteBulkSettings {
                maximum_volume_ratio: 1e10,
                ..settings(60.0)
            },
        )
        .is_err()
    );
}
