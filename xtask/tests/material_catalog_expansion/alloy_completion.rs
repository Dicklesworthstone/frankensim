//! Source-specific alloy stiffness and expansion reach their matching consumers.
use super::*;

fn reference() -> QueryPoint {
    flags(&[
        ("engineering-reference-model", 1.0),
        ("source-specimens-matched", 0.0),
    ])
}

#[test]
fn a36_and_1018_references_drive_axial_and_slab_models_without_invented_temperature_support() {
    use fs_conduction::transient::{MaterialHeatCapacity, MaterialHeatCapacitySource};
    use fs_matdb::{MaterialStateId, NormalizedMaterialCardPack};
    use fs_material::axial_reference::{
        AxialLoadingMode, AxialReferenceRequest, ResolvedAxialModulusReference,
    };
    use fs_material::state_point::{
        MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    };

    // Independent conversions from the cited supplier's printed imperial facts.
    let psi = 0.45359237 * 9.80665 / 0.0254_f64.powi(2);
    let modulus = 29_000_000.0 * psi;
    let cp = 0.11 * 1055.05585262 / (0.45359237 * (5.0 / 9.0));
    for (
        slug,
        process_axis,
        pounds_per_inch3,
        nu,
        shear_ksi,
        k_imperial,
        uts,
        ys,
        elongation,
        area,
    ) in [
        (
            "steel-a36-onlinemetals-reference",
            "source-hot-rolled",
            0.282,
            0.26,
            11_500.0,
            348.0,
            58.0,
            47.7,
            0.22,
            0.50,
        ),
        (
            "steel-1018-onlinemetals-reference",
            "source-cold-finished",
            0.284,
            0.29,
            11_300.0,
            360.0,
            63.8,
            53.7,
            0.15,
            0.40,
        ),
    ] {
        let pack = compile(slug, 10);
        let rho = pounds_per_inch3 * 0.45359237 / 0.0254_f64.powi(3);
        let conductivity =
            k_imperial * 1055.05585262 * 0.0254 / (3600.0 * 0.3048_f64.powi(2) * (5.0 / 9.0));
        let at = flags(&[
            ("engineering-reference-model", 1.0),
            ("source-specimens-matched", 0.0),
            ("source-test-temperature-known", 0.0),
            ("source-test-method-known", 0.0),
            (process_axis, 1.0),
        ]);
        for (name, dims, value) in [
            ("density", RHO, rho),
            ("young_modulus", PA, modulus),
            ("poisson_ratio", Dims::NONE, nu),
            ("shear-modulus", PA, shear_ksi * 1000.0 * psi),
            ("specific-heat-capacity", CP, cp),
            ("thermal-conductivity", K, conductivity),
            ("tensile-strength", PA, uts * 1000.0 * psi),
            ("tensile-yield-strength", PA, ys * 1000.0 * psi),
            ("elongation-at-break", Dims::NONE, elongation),
            ("reduction-of-area", Dims::NONE, area),
        ] {
            let answer = sample(&pack, name, dims, &at, value);
            assert_eq!(
                answer.evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
            for bad in [
                QueryPoint::new(),
                at.clone().with("engineering-reference-model", 0.0).unwrap(),
                at.clone().with("source-specimens-matched", 1.0).unwrap(),
                at.clone()
                    .with("source-test-temperature-known", 1.0)
                    .unwrap(),
                at.clone().with("source-test-method-known", 1.0).unwrap(),
                at.clone().with(process_axis, 0.0).unwrap(),
            ] {
                refuses(&pack, name, dims, &bad);
            }
        }
        let material = NormalizedMaterialCardPack::new(
            MaterialStateId {
                chemistry: slug.into(),
                phase: "solid".into(),
                process: "unmatched supplier-guide engineering reference; test temperature unknown"
                    .into(),
                revision: 0,
            },
            pack.clone(),
        )
        .unwrap();
        let requirement = |name, dims| {
            ScalarPropertyRequirement::try_with_key(
                &key(name, dims),
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap()
        };
        let axial = ResolvedAxialModulusReference::resolve(
            material.card(),
            &at,
            &AxialReferenceRequest {
                density: requirement("density", RHO),
                modulus: requirement("young_modulus", PA),
                mode: AxialLoadingMode::Tension,
                maximum_absolute_strain: 0.0002,
                crosshead_speed_axis: None,
            },
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let mut coupon = axial.coupon(0.2, 2e-4).unwrap();
        let response = coupon.advance_displacement(2e-5, 1.0).unwrap();
        assert!((response.mass_kg / (rho * 4e-5) - 1.0).abs() < 1e-13);
        assert!((response.reference_force_n / (modulus * 2e-8) - 1.0).abs() < 1e-13);
        assert!((response.model_loading_work_j / (modulus * 2e-13) - 1.0).abs() < 1e-13);
        assert!(coupon.advance_displacement(5e-5, 1.0).is_err());
        let thermal = fs_conduction::ThermalResistanceTerm::slab_from_card(
            "10 mm supplier-reference steel plate",
            0.01,
            0.02,
            material.card(),
            &key("thermal-conductivity", K),
            &at,
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
        assert!((thermal.value_k_per_w() / (0.5 / conductivity) - 1.0).abs() < 1e-13);

        // An unknown test temperature must not become a validity-wide heat law.
        let density = key("density", RHO);
        let heat_capacity = key("specific-heat-capacity", CP);
        assert!(
            MaterialHeatCapacity::from_card_at_query_points(
                material.card(),
                MaterialHeatCapacitySource {
                    density: &density,
                    specific_heat: &heat_capacity,
                    temperature_axis: "temperature",
                    fixed_source_selectors: &[
                        "engineering-reference-model",
                        "source-specimens-matched",
                        "source-test-temperature-known",
                        "source-test-method-known",
                        process_axis,
                    ],
                },
                &temperature(at, "temperature", 293.15),
                [293.15, 294.15],
                MaterialPropertySelection::SingleClaimOnly,
            )
            .is_err()
        );
    }
}

#[test]
fn inconel_625_dynamic_curves_keep_heat_treatments_and_method_separate() {
    let pack = compile("inconel-625-special-metals-2013", 10);
    let condition = |kelvin, solution| {
        temperature(
            flags(&[
                ("engineering-reference-model", 1.0),
                ("source-specimens-matched", 0.0),
                ("source-solution-treated", solution),
                ("source-dynamic-method", 1.0),
                ("source-frequency-known", 0.0),
            ])
            .with_quantity(
                "source-rod-diameter",
                QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0])),
                0.01905,
            )
            .unwrap(),
            "temperature",
            kelvin,
        )
    };
    for (t, treatment, e, g, nu) in [
        (294.15, 0.0, 207.5e9, 81.4e9, 0.278),
        (294.15, 1.0, 204.8e9, 78.0e9, 0.312),
        (330.15, 0.0, 205.8e9, 80.7e9, 0.279),
        (330.15, 1.0, 202.7e9, 77.25e9, 0.3115),
        (1144.15, 0.0, 147.5e9, 55.2e9, 0.336),
        (1144.15, 1.0, 148.2e9, 57.2e9, 0.289),
    ] {
        let at = condition(t, treatment);
        for (name, dims, value) in [
            ("dynamic-young-modulus", PA, e),
            ("dynamic-shear-modulus", PA, g),
            ("dynamic-poisson-ratio", Dims::NONE, nu),
        ] {
            let answer = sample(&pack, name, dims, &at, value);
            assert_eq!(
                answer.evidence.numerical.kind,
                fs_evidence::NumericalKind::Estimate
            );
        }
        // A dynamic measurement must not silently become a static law.
        refuses(&pack, "young_modulus", PA, &at);
        refuses(&pack, "tensile-modulus", PA, &at);
    }
    let base = condition(294.15, 0.0);
    for bad in [
        QueryPoint::new(),
        condition(294.14, 0.0),
        condition(1144.16, 0.0),
        condition(330.15, 0.5),
        base.clone().with("source-dynamic-method", 0.0).unwrap(),
        base.clone().with("source-frequency-known", 1.0).unwrap(),
        base.clone().with("source-rod-diameter", 0.02).unwrap(),
        base.clone().with("source-specimens-matched", 1.0).unwrap(),
    ] {
        refuses(&pack, "dynamic-young-modulus", PA, &bad);
    }
    // Table 4's unspecified anneal cannot supply Table 3's thermal state.
    refuses(&pack, "thermal-conductivity", K, &base);
}

#[test]
fn monel_reference_drives_axial_force_mass_and_work_without_isotropic_fusion() {
    use fs_matdb::{MaterialStateId, NormalizedMaterialCardPack};
    use fs_material::axial_reference::{
        AxialLoadingMode, AxialReferenceRequest, ResolvedAxialModulusReference,
    };
    use fs_material::state_point::{
        MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    };

    let pack = compile("monel-400-special-metals-2005", 8);
    let at = temperature(
        flags(&[
            ("engineering-reference-model", 1.0),
            ("source-specimens-matched", 0.0),
            ("source-elastic-temperature-known", 0.0),
            ("source-elastic-method-known", 0.0),
        ]),
        "temperature",
        294.15,
    );
    let psi = 0.45359237 * 9.80665 / (0.0254 * 0.0254);
    let e = 26e6 * psi;
    sample(&pack, "tensile-modulus", PA, &at, e);
    sample(&pack, "compression-modulus", PA, &at, e);
    sample(&pack, "torsional-modulus", PA, &at, 9.5e6 * psi);
    sample(&pack, "poisson-ratio", Dims::NONE, &at, 0.32);
    assert!((e / (2.0 * 9.5e6 * psi) - 1.0 - 0.32).abs() > 0.04);
    let material = NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: "MONEL 400 UNS N04400".into(),
            phase: "solid".into(),
            process: "SMC-053 nominal engineering reference; source specimens unmatched".into(),
            revision: 0,
        },
        pack.clone(),
    )
    .unwrap();
    for (mode, name, sign) in [
        (AxialLoadingMode::Tension, "tensile-modulus", 1.0),
        (AxialLoadingMode::Compression, "compression-modulus", -1.0),
    ] {
        let requirement = |name, dims| {
            ScalarPropertyRequirement::try_with_key(
                &key(name, dims),
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap()
        };
        let request = AxialReferenceRequest {
            density: requirement("density", RHO),
            modulus: requirement(name, PA),
            mode,
            maximum_absolute_strain: 0.001,
            crosshead_speed_axis: None,
        };
        let resolved = ResolvedAxialModulusReference::resolve(
            material.card(),
            &at,
            &request,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let mut coupon = resolved.coupon(0.1, 1e-4).unwrap();
        let response = coupon.advance_displacement(sign * 1e-5, 1.0).unwrap();
        assert!((response.mass_kg - 0.088).abs() < 1e-14);
        assert!((response.reference_force_n - sign * e * 1e-8).abs() < 1e-9);
        assert!((response.model_loading_work_j - e * 5e-14).abs() < 1e-14);
    }
    for bad in [
        temperature(reference(), "temperature", 294.15),
        at.clone().with("temperature", 300.0).unwrap(),
        at.clone()
            .with("source-elastic-temperature-known", 1.0)
            .unwrap(),
        at.clone().with("source-elastic-method-known", 1.0).unwrap(),
    ] {
        refuses(&pack, "tensile-modulus", PA, &bad);
    }
}

#[test]
fn timetal_interval_means_drive_endpoint_expansion_with_distinct_reference_lengths() {
    use fs_matdb::{MaterialStateId, NormalizedMaterialCardPack};
    use fs_material::state_point::{
        MaterialPropertySelection, resolve_mean_thermal_expansion_strain,
    };

    let pack = compile("timetal-6-4-timet-2000", 7);
    let material = NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: "TIMETAL 6-4 ASTM Grade 5".into(),
            phase: "solid".into(),
            process: "TMC-0150 nominal interval means; direction unspecified".into(),
            revision: 0,
        },
        pack.clone(),
    )
    .unwrap();
    let coefficient = key(
        "mean-linear-thermal-expansion-coefficient",
        Dims([0, 0, 0, -1, 0, 0]),
    );
    let interval = |start, end| {
        temperature(
            temperature(
                reference()
                    .with_quantity(
                        "source-expansion-direction-known",
                        QuantitySpec::dimensional(Dims::NONE),
                        0.0,
                    )
                    .unwrap(),
                "source-reference-temperature",
                start,
            ),
            "source-range-end-temperature",
            end,
        )
    };
    for (start, end, alpha, strain) in [
        (273.15, 373.15, 9e-6, 0.0009),
        (293.15, 698.15, 9.4e-6, 0.003807),
        (293.15, 923.15, 9.7e-6, 0.006111),
    ] {
        let at = interval(start, end);
        sample(&pack, coefficient.name(), coefficient.dims(), &at, alpha);
        let response = resolve_mean_thermal_expansion_strain(
            material.card(),
            &coefficient,
            &at,
            "source-reference-temperature",
            "source-range-end-temperature",
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        assert!((response.engineering_strain() - strain).abs() < 1e-14);
        // A two-metre reference segment lengthens by 1.8, 7.614 or 12.222 mm.
        assert!((2.0 * response.engineering_strain() - 2.0 * strain).abs() < 1e-14);
    }
    for bad in [
        interval(273.15, 698.15),
        interval(293.15, 800.0),
        interval(293.15, 923.16),
        interval(293.15, 923.15)
            .with("source-expansion-direction-known", 1.0)
            .unwrap(),
    ] {
        refuses(&pack, coefficient.name(), coefficient.dims(), &bad);
    }
    refuses(
        &pack,
        "linear-thermal-expansion-coefficient",
        coefficient.dims(),
        &interval(293.15, 923.15),
    );
}
