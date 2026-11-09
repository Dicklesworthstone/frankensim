use super::*;
use crate::elastic::acoustic::viscoelastic::ViscoelasticSlabSource;
use crate::visco::GeneralizedMaxwell;
use fs_blake3::ContentHash;
use fs_math::det;

fn laws() -> [GeneralizedMaxwell; 4] {
    [
        GeneralizedMaxwell {
            e_inf: 60.0,
            terms: vec![(20.0, 0.05)],
        },
        GeneralizedMaxwell {
            e_inf: 30.0,
            terms: vec![(15.0, 0.02)],
        },
        GeneralizedMaxwell {
            e_inf: 130.0,
            terms: vec![(15.0, 0.02)],
        },
        GeneralizedMaxwell {
            e_inf: 45.0,
            terms: vec![(27.0, 0.03)],
        },
    ]
}

fn source<'a>(
    bulk: &'a GeneralizedMaxwell,
    shear: &'a GeneralizedMaxwell,
    density: f64,
    identity: u8,
) -> ViscoelasticSlabSource<'a> {
    ViscoelasticSlabSource {
        bulk,
        shear,
        density_kg_m3: density,
        frequency_band_rad_s: (10.0, 100.0),
        maximum_absolute_strain: 0.1,
        source_identity: ContentHash([identity; 32]),
    }
}

fn layers(laws: &[GeneralizedMaxwell; 4]) -> [ViscoelasticLayer<'_>; 2] {
    [
        ViscoelasticLayer {
            source: source(&laws[0], &laws[1], 4.0, 3),
            thickness_m: 0.2,
        },
        ViscoelasticLayer {
            source: source(&laws[2], &laws[3], 7.0, 8),
            thickness_m: 0.3,
        },
    ]
}

fn budget(count: usize) -> ViscoelasticLayerBudget {
    ViscoelasticLayerBudget {
        maximum_layers: count,
        maximum_total_branches: 2 * count,
    }
}

fn close(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() <= 3.0e-12 * expected.abs().max(1.0e-300),
        "actual={actual:?}, expected={expected:?}"
    );
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 5.0e-12 * expected.abs().max(1.0e-300),
        "actual={actual:e}, expected={expected:e}"
    );
}

fn phase(distance: f64) -> C64 {
    let (sine, cosine) = det::sin_cos(4.0 * distance);
    C64::new(cosine, -sine)
}

#[test]
fn g1_supported_original_layers_match_independent_boundary_and_positive_work_reference() {
    let laws = laws();
    let layers = layers(&laws);
    let fluid = FluidAcousticState::try_new(2.0, 10.0).unwrap();
    let backing = ViscoelasticBacking::Impedance(C64::new(13.0, 7.0));
    let prepared = SupportedViscoelasticLayers::try_new(
        &layers,
        &fluid,
        backing,
        20.0,
        budget(2),
        &mut || false,
    )
    .unwrap();
    assert_eq!(prepared.layer_count(), 2);
    assert_eq!(prepared.incident_state(), fluid);
    assert_eq!(prepared.backing(), backing);
    assert_eq!(prepared.angular_frequency_rad_s(), 20.0);
    assert!(core::ptr::eq(
        prepared.incident_source(),
        &fluid as &dyn FluidAcoustics
    ));
    for (i, original) in layers.iter().enumerate() {
        let retained = prepared.layer(i).unwrap();
        assert!(core::ptr::eq(retained.source.bulk, original.source.bulk));
        assert!(core::ptr::eq(retained.source.shear, original.source.shear));
        assert_eq!(retained.source.source_identity, original.source.source_identity);
        let expected = [original.source.bulk, original.source.shear].map(|law| {
            let (storage, loss) = law.modulus(20.0);
            C64::new(storage, loss)
        });
        assert_eq!(prepared.layer_moduli_pa(i).unwrap(), expected);
    }

    let driven = prepared.drive(C64::ONE, 1.0, &mut || false).unwrap();
    assert!(core::ptr::eq(driven.prepared(), &prepared));
    assert_eq!(driven.area_m2(), 1.0);

    // Independent 80-digit transfer solve and quadrature of the ORIGINAL bulk
    // and shear strain work; the two heat channels are not flux differences.
    let reflected = C64::new(0.21661124142492808, -0.4698365313672166);
    close(prepared.pressure_reflection_coefficient(), reflected);
    close(driven.reflected_pressure_pa(), reflected);
    close(
        driven.support_compressive_traction_pa(),
        C64::new(0.008713485141651073, -0.5957053709938446),
    );
    close(
        driven.support_normal_velocity_m_s(),
        C64::new(-0.01860854261520848, -0.03580350559133733),
    );
    close(
        driven.support_compressive_traction_pa(),
        C64::new(13.0, 7.0) * driven.support_normal_velocity_m_s(),
    );
    let heat = [
        [0.0033753594392632683, 0.0023278340960436333],
        [0.0004967436606334002, 0.0015252952402978522],
    ];
    for (i, expected) in heat.into_iter().enumerate() {
        let local = driven.layer_absorption(i, &mut || false).unwrap();
        assert_eq!(local.source_identity, layers[i].source.source_identity);
        assert!(local.maximum_absolute_strain_bound <= 0.1);
        near(local.bulk_absorbed_power_w, expected[0]);
        near(local.shear_absorbed_power_w, expected[1]);
        let depth = layers[i].thickness_m;
        let left = driven
            .absorption_between(i, [0.0, 0.4 * depth], &mut || false)
            .unwrap();
        let right = driven
            .absorption_between(i, [0.4 * depth, depth], &mut || false)
            .unwrap();
        near(
            left.bulk_absorbed_power_w + right.bulk_absorbed_power_w,
            expected[0],
        );
        near(
            left.shear_absorbed_power_w + right.shear_absorbed_power_w,
            expected[1],
        );
    }
    let power = driven.power();
    near(power.incident_power_w, 0.025);
    near(power.reflected_power_w, 0.006691669902970649);
    near(power.backing_absorbed_power_w, 0.010583097660791197);
    near(
        power.bulk_absorbed_power_w,
        0.0033753594392632683 + 0.0004967436606334002,
    );
    near(
        power.shear_absorbed_power_w,
        0.0023278340960436333 + 0.0015252952402978522,
    );
    assert!(power.balance_residual_w.abs() < 2.0e-14 * power.incident_power_w);

    let front = driven.field(0, 0.0, &mut || false).unwrap();
    close(front.compressive_traction_pa, C64::ONE + reflected);
    close(front.normal_velocity_m_s, (C64::ONE - reflected).scale(0.05));
    let left_join = driven.field(0, 0.2, &mut || false).unwrap();
    let right_join = driven.field(1, 0.0, &mut || false).unwrap();
    close(left_join.compressive_traction_pa, right_join.compressive_traction_pa);
    close(left_join.normal_velocity_m_s, right_join.normal_velocity_m_s);
    let rear = driven.field(1, 0.3, &mut || false).unwrap();
    assert_eq!(
        rear.compressive_traction_pa,
        driven.support_compressive_traction_pa()
    );
    assert_eq!(rear.normal_velocity_m_s, driven.support_normal_velocity_m_s());
}

#[test]
fn g1_empty_zero_thickness_and_lossless_matched_coatings_keep_original_phase() {
    let bulk = GeneralizedMaxwell {
        e_inf: 60.0,
        terms: vec![],
    };
    let shear = GeneralizedMaxwell {
        e_inf: 30.0,
        terms: vec![],
    };
    // M=100 Pa, rho=4 kg/m3: Z=20 Pa s/m and k=4 rad/m at omega=20.
    let original = source(&bulk, &shear, 4.0, 4);
    let fluid = FluidAcousticState::try_new(2.0, 10.0).unwrap();
    let thicknesses: &[&[f64]] = &[&[], &[0.0], &[0.3], &[0.12, 0.18]];
    for backing in [
        ViscoelasticBacking::TractionFree,
        ViscoelasticBacking::Clamped,
        ViscoelasticBacking::Impedance(C64::new(13.0, 7.0)),
        ViscoelasticBacking::Impedance(C64::new(0.0, 7.0)),
        ViscoelasticBacking::Impedance(C64::from_re(20.0)),
    ] {
        let (reflection, rear_pressure, rear_velocity) = match backing {
            ViscoelasticBacking::TractionFree => (
                -C64::ONE,
                C64::ZERO,
                C64::from_re(0.1),
            ),
            ViscoelasticBacking::Clamped => (
                C64::ONE,
                C64::from_re(2.0),
                C64::ZERO,
            ),
            ViscoelasticBacking::Impedance(z) => {
                let denominator = z + C64::from_re(20.0);
                (
                    (z - C64::from_re(20.0)) / denominator,
                    z.scale(2.0) / denominator,
                    C64::from_re(2.0) / denominator,
                )
            }
        };
        for lengths in thicknesses {
            let layers: Vec<_> = lengths
                .iter()
                .map(|&thickness_m| ViscoelasticLayer {
                    source: original,
                    thickness_m,
                })
                .collect();
            let distance: f64 = lengths.iter().sum();
            let prepared = SupportedViscoelasticLayers::try_new(
                &layers,
                &fluid,
                backing,
                20.0,
                budget(layers.len()),
                &mut || false,
            )
            .unwrap();
            let driven = prepared.drive(C64::ONE, 1.0, &mut || false).unwrap();
            let expected_reflection = reflection * phase(2.0 * distance);
            if expected_reflection == C64::ZERO {
                assert!(driven.reflected_pressure_pa().abs() < 2.0e-14);
            } else {
                close(driven.reflected_pressure_pa(), expected_reflection);
            }
            close(
                driven.support_compressive_traction_pa(),
                rear_pressure * phase(distance),
            );
            close(
                driven.support_normal_velocity_m_s(),
                rear_velocity * phase(distance),
            );
            let power = driven.power();
            assert_eq!(power.bulk_absorbed_power_w, 0.0);
            assert_eq!(power.shear_absorbed_power_w, 0.0);
            let expected_backing = match backing {
                ViscoelasticBacking::Impedance(z) => {
                    0.5 * z.re * rear_velocity.abs().powi(2)
                }
                _ => 0.0,
            };
            near(power.backing_absorbed_power_w, expected_backing);
            assert!(power.balance_residual_w.abs() < 2.0e-14 * power.incident_power_w);
            if distance == 0.0 {
                assert!(driven.field(0, 0.0, &mut || false).is_err());
            }
        }
    }
}

#[test]
fn g1_free_clamped_and_very_stiff_backing_keep_actual_rear_motion_and_work() {
    let laws = laws();
    let layers = layers(&laws);
    let fluid = FluidAcousticState::try_new(2.0, 10.0).unwrap();
    for (backing, reflected, velocity, intrinsic, support_work) in [
        (
            ViscoelasticBacking::TractionFree,
            C64::new(0.4904617257586931, -0.6046444356642918),
            C64::new(-0.023297401014702475, -0.04563001812531362),
            0.009846310049650366,
            0.0,
        ),
        (
            ViscoelasticBacking::Clamped,
            C64::new(-0.7806572538207887, 0.06407651098154609),
            C64::ZERO,
            0.009661711319935414,
            0.0,
        ),
        (
            ViscoelasticBacking::Impedance(C64::new(13.0e30, 7.0e30)),
            C64::new(-0.7806572538207887, 0.06407651098154609),
            C64::new(-9.205606637296277e-32, -1.6606316675136072e-31),
            0.009661711319935414,
            2.3433341559909816e-31,
        ),
    ] {
        let prepared = SupportedViscoelasticLayers::try_new(
            &layers,
            &fluid,
            backing,
            20.0,
            budget(2),
            &mut || false,
        )
        .unwrap();
        let driven = prepared.drive(C64::ONE, 1.0, &mut || false).unwrap();
        close(driven.reflected_pressure_pa(), reflected);
        close(driven.support_normal_velocity_m_s(), velocity);
        let power = driven.power();
        near(
            power.bulk_absorbed_power_w + power.shear_absorbed_power_w,
            intrinsic,
        );
        near(power.backing_absorbed_power_w, support_work);
        let rear = driven.field(1, 0.3, &mut || false).unwrap();
        assert_eq!(rear.normal_velocity_m_s, driven.support_normal_velocity_m_s());
        match backing {
            ViscoelasticBacking::TractionFree => {
                assert_eq!(driven.support_compressive_traction_pa(), C64::ZERO);
            }
            ViscoelasticBacking::Clamped => {
                close(
                    driven.support_compressive_traction_pa(),
                    C64::new(-0.034286695588991044, -2.803213632378429),
                );
            }
            ViscoelasticBacking::Impedance(z) => {
                assert!(power.backing_absorbed_power_w > 0.0);
                close(
                    driven.support_compressive_traction_pa(),
                    z * driven.support_normal_velocity_m_s(),
                );
            }
        }
    }
}

#[test]
fn g1_weak_intrinsic_heat_and_opaque_shielding_are_not_flux_deficits() {
    let bulk = GeneralizedMaxwell {
        e_inf: 60.0,
        terms: vec![(1.0e-14, 0.05)],
    };
    let shear = GeneralizedMaxwell {
        e_inf: 30.0,
        terms: vec![(1.0e-14, 0.02)],
    };
    let thin = [ViscoelasticLayer {
        source: source(&bulk, &shear, 4.0, 9),
        thickness_m: 1.0e-8,
    }];
    let fluid = FluidAcousticState::try_new(2.0, 10.0).unwrap();
    let prepared = SupportedViscoelasticLayers::try_new(
        &thin,
        &fluid,
        ViscoelasticBacking::Clamped,
        20.0,
        budget(1),
        &mut || false,
    )
    .unwrap();
    let driven = prepared.drive(C64::ONE, 1.0, &mut || false).unwrap();
    let power = driven.power();
    // For k*d=4e-8, |p|^2=4+O((k*d)^2) and |strain|^2=4/100^2.
    let bulk_reference = 20.0 * bulk.modulus(20.0).1 * 4.0 * 1.0e-8 / (2.0 * 10000.0);
    let shear_reference =
        (2.0 / 3.0) * 20.0 * shear.modulus(20.0).1 * 4.0 * 1.0e-8 / 10000.0;
    near(power.bulk_absorbed_power_w, bulk_reference);
    near(power.shear_absorbed_power_w, shear_reference);
    for heat in [power.bulk_absorbed_power_w, power.shear_absorbed_power_w] {
        assert!(heat > 0.0 && heat < f64::EPSILON * power.incident_power_w);
    }

    let laws = laws();
    let mut opaque = layers(&laws);
    opaque[0].thickness_m = 1.0e5;
    let prepared = SupportedViscoelasticLayers::try_new(
        &opaque,
        &fluid,
        ViscoelasticBacking::Impedance(C64::new(13.0, 7.0)),
        20.0,
        budget(2),
        &mut || false,
    )
    .unwrap();
    let driven = prepared.drive(C64::ONE, 1.0, &mut || false).unwrap();
    let moduli = prepared.layer_moduli_pa(0).unwrap();
    let z = (moduli[0] + moduli[1].scale(4.0 / 3.0)).sqrt().scale(2.0);
    close(
        driven.reflected_pressure_pa(),
        (z - C64::from_re(20.0)) / (z + C64::from_re(20.0)),
    );
    assert_eq!(driven.support_compressive_traction_pa(), C64::ZERO);
    assert_eq!(driven.support_normal_velocity_m_s(), C64::ZERO);
    assert_eq!(driven.power().backing_absorbed_power_w, 0.0);
    let hidden = driven.layer_absorption(1, &mut || false).unwrap();
    assert_eq!(hidden.bulk_absorbed_power_w, 0.0);
    assert_eq!(hidden.shear_absorbed_power_w, 0.0);
}

#[test]
fn g3_g4_supported_source_domains_scaling_and_cancellation_are_binding() {
    let laws = laws();
    let layers = layers(&laws);
    let fluid = FluidAcousticState::try_new(2.0, 10.0).unwrap();
    let backing = ViscoelasticBacking::Impedance(C64::new(13.0, 7.0));
    let prepared = SupportedViscoelasticLayers::try_new(
        &layers,
        &fluid,
        backing,
        20.0,
        budget(2),
        &mut || false,
    )
    .unwrap();
    let unit = prepared.drive(C64::ONE, 1.0, &mut || false).unwrap();
    let amplitude = C64::new(0.0, 2.0);
    let scaled = prepared.drive(amplitude, 0.3, &mut || false).unwrap();
    close(scaled.reflected_pressure_pa(), amplitude * unit.reflected_pressure_pa());
    close(
        scaled.support_normal_velocity_m_s(),
        amplitude * unit.support_normal_velocity_m_s(),
    );
    close(
        scaled.support_compressive_traction_pa(),
        amplitude * unit.support_compressive_traction_pa(),
    );
    let a = unit.power();
    let b = scaled.power();
    for (base, actual) in [
        (a.incident_power_w, b.incident_power_w),
        (a.reflected_power_w, b.reflected_power_w),
        (a.bulk_absorbed_power_w, b.bulk_absorbed_power_w),
        (a.shear_absorbed_power_w, b.shear_absorbed_power_w),
        (a.backing_absorbed_power_w, b.backing_absorbed_power_w),
    ] {
        near(actual, 1.2 * base);
    }
    assert!(matches!(
        prepared.drive(C64::from_re(1.0e5), 1.0, &mut || false),
        Err(ViscoelasticSlabError::OutsideSourceSupport)
    ));
    assert!(matches!(
        prepared.drive(C64::ONE, 1.0, &mut || true),
        Err(ViscoelasticSlabError::Cancelled)
    ));
    assert!(matches!(
        unit.field(0, 0.1, &mut || true),
        Err(ViscoelasticSlabError::Cancelled)
    ));
    assert!(unit.field(0, 0.21, &mut || false).is_err());
    assert!(unit.absorption_between(0, [0.1, 0.21], &mut || false).is_err());
    assert!(unit.layer_absorption(2, &mut || false).is_err());
    assert!(matches!(
        SupportedViscoelasticLayers::try_new(
            &layers,
            &fluid,
            backing,
            20.0,
            budget(2),
            &mut || true,
        ),
        Err(ViscoelasticSlabError::Cancelled)
    ));
    assert!(matches!(
        SupportedViscoelasticLayers::try_new(
            &layers,
            &fluid,
            backing,
            20.0,
            budget(1),
            &mut || false,
        ),
        Err(ViscoelasticSlabError::Budget)
    ));
    for invalid in [
        C64::new(-1.0, 0.0),
        C64::new(f64::NAN, 0.0),
        C64::new(0.0, f64::INFINITY),
    ] {
        assert!(
            SupportedViscoelasticLayers::try_new(
                &layers,
                &fluid,
                ViscoelasticBacking::Impedance(invalid),
                20.0,
                budget(2),
                &mut || false,
            )
            .is_err()
        );
    }
    let mut zero = layers[0];
    zero.thickness_m = 0.0;
    assert!(matches!(
        SupportedViscoelasticLayers::try_new(
            &[zero],
            &fluid,
            backing,
            5.0,
            budget(1),
            &mut || false,
        ),
        Err(ViscoelasticSlabError::OutsideSourceSupport)
    ));
}
