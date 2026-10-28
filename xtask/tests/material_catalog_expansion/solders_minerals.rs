//! Solder and rock references retain specimen, method, phase and state meaning.
use super::*;

use fs_matdb::{MaterialStateId, NormalizedMaterialCardPack};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_mean_thermal_expansion_strain,
};

const SPECIFIC_ENERGY: Dims = Dims([2, 0, -2, 0, 0, 0]);
const MASS: Dims = Dims([0, 1, 0, 0, 0, 0]);
const TEMPERATURE_RATE: Dims = Dims([0, 0, -1, 1, 0, 0]);
const EXPANSION: Dims = Dims([0, 0, 0, -1, 0, 0]);

fn replace_flag(point: QueryPoint, axis: &str, value: f64) -> QueryPoint {
    let quantity = QuantitySpec::dimensional(Dims::NONE);
    assert_eq!(point.axis_quantities().get(axis), Some(&quantity));
    point.with_quantity(axis, quantity, value).unwrap()
}

fn reference(alloy_axis: &str) -> QueryPoint {
    flags(&[
        ("engineering-reference-model", 1.0),
        ("source-specimens-matched", 0.0),
        ("source-test-method-known", 0.0),
        (alloy_axis, 1.0),
    ])
}

fn unknown_temperature(alloy_axis: &str) -> QueryPoint {
    reference(alloy_axis)
        .with_quantity(
            "source-test-temperature-known",
            QuantitySpec::dimensional(Dims::NONE),
            0.0,
        )
        .unwrap()
}

fn card(pack: &NormalizedPack) -> NormalizedMaterialCardPack {
    NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: pack.pack_id().into(),
            phase: "source-specific solder reference".into(),
            process: "retain original supplier or DSC protocol; no qualified joint law".into(),
            revision: 0,
        },
        pack.clone(),
    )
    .unwrap()
}

fn absolute_key(name: &str) -> PropertyKey {
    PropertyKey::with_quantity(
        name,
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        )),
    )
}

fn phase_endpoint(pack: &NormalizedPack, name: &str, at: &QueryPoint, expected: f64) {
    let answer = pack
        .claims()
        .query_typed(&absolute_key(name), at, SelectionPolicy::SingleClaimOnly)
        .unwrap();
    assert!((answer.evidence.value.value - expected).abs() < 1e-12);
    assert_eq!(answer.evidence.value.dims, TEMP);
    assert_eq!(
        answer.evidence.value.uncertainty,
        UncertaintyModel::Unstated
    );
    assert!(answer.receipt.observation_backed);
    assert!(!answer.receipt.source_hashes.is_empty());
    pack.claims().verify_receipt(&answer.receipt).unwrap();
    // A temperature dimension alone must not erase the absolute-temperature meaning.
    refuses(pack, name, TEMP, at);
    assert!(
        pack.claims()
            .query_typed(
                &absolute_key(name),
                &QueryPoint::new(),
                SelectionPolicy::SingleClaimOnly,
            )
            .is_err()
    );
}

#[test]
fn sac305_manufacturer_values_reach_bounded_axial_and_thermal_references() {
    use fs_material::axial_reference::{
        AxialLoadingMode, AxialReferenceRequest, ResolvedAxialModulusReference,
    };

    let pack = compile("solder-sac305-indium-2024", 10);
    let at = unknown_temperature("source-alloy-sac305");
    phase_endpoint(&pack, "solidus-temperature", &at, 217.0 + 273.15);
    phase_endpoint(&pack, "liquidus-temperature", &at, 220.0 + 273.15);
    for (name, dims, value) in [
        ("density", RHO, 7.35 * 1000.0),
        ("thermal-conductivity", K, 61.0),
        (
            "electrical-resistivity",
            Dims([3, 1, -3, 0, -2, 0]),
            13.2e-8,
        ),
        ("young_modulus", PA, 6.9e9),
        ("poisson_ratio", Dims::NONE, 0.35),
        ("tensile-strength", PA, 39.58e6),
        ("tensile-yield-strength", PA, 22e6),
        ("elongation-at-break", Dims::NONE, 0.49),
    ] {
        let answer = sample(&pack, name, dims, &at, value);
        assert_eq!(
            answer.evidence.value.uncertainty,
            UncertaintyModel::Unstated
        );
        for bad in [
            QueryPoint::new(),
            unknown_temperature("source-alloy-sn63pb37"),
            replace_flag(at.clone(), "source-alloy-sac305", 0.0),
            replace_flag(at.clone(), "source-specimens-matched", 1.0),
            replace_flag(at.clone(), "source-test-temperature-known", 1.0),
            replace_flag(at.clone(), "source-test-method-known", 1.0),
        ] {
            refuses(&pack, name, dims, &bad);
        }
    }
    let material = card(&pack);
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
            maximum_absolute_strain: 1e-4,
            crosshead_speed_axis: None,
        },
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let (length, area, displacement) = (0.02, 1e-6, 1e-7);
    let mut coupon = axial.coupon(length, area).unwrap();
    let response = coupon.advance_displacement(displacement, 1.0).unwrap();
    assert!((response.mass_kg / (7350.0 * length * area) - 1.0).abs() < 1e-13);
    assert!(
        (response.reference_force_n / (6.9e9 * area * displacement / length) - 1.0).abs() < 1e-13
    );
    assert!(
        (response.model_loading_work_j
            / (0.5 * 6.9e9 * area * displacement * displacement / length)
            - 1.0)
            .abs()
            < 1e-13
    );
    assert!(coupon.advance_displacement(3e-6, 1.0).is_err());
    let slab = fs_conduction::ThermalResistanceTerm::slab_from_card(
        "SAC305 nominal reference slab",
        0.001,
        1e-4,
        material.card(),
        &key("thermal-conductivity", K),
        &at,
        SelectionPolicy::SingleClaimOnly,
    )
    .unwrap();
    assert!((slab.value_k_per_w() - 0.001 / (61.0 * 1e-4)).abs() < 1e-14);
    refuses(&pack, "latent-heat-fusion", SPECIFIC_ENERGY, &at);
    refuses(&pack, "creep-strain-rate", Dims([0, 0, -1, 0, 0, 0]), &at);
}

#[test]
fn kester_mean_expansion_has_exact_endpoints_and_no_invented_heat_capacity_path() {
    use fs_conduction::transient::{MaterialHeatCapacity, MaterialHeatCapacitySource};

    let pack = compile("solder-sac305-kester-2023", 5);
    let unknown = unknown_temperature("source-alloy-sac305");
    phase_endpoint(&pack, "melting-range-lower-temperature", &unknown, 490.15);
    phase_endpoint(&pack, "melting-range-upper-temperature", &unknown, 492.15);
    sample(&pack, "density", RHO, &unknown, 7370.0);
    sample(
        &pack,
        "specific-heat-capacity",
        CP,
        &unknown,
        0.232 * 1000.0,
    );
    assert!(
        pack.claims()
            .query_typed(
                &absolute_key("liquidus-temperature"),
                &unknown,
                SelectionPolicy::SingleClaimOnly,
            )
            .is_err()
    );
    let at = temperature(
        temperature(
            reference("source-alloy-sac305"),
            "source-reference-temperature",
            293.15,
        ),
        "source-range-end-temperature",
        373.15,
    );
    let coefficient = key("mean-linear-thermal-expansion-coefficient", EXPANSION);
    sample(&pack, coefficient.name(), coefficient.dims(), &at, 21.9e-6);
    let material = card(&pack);
    let strain = resolve_mean_thermal_expansion_strain(
        material.card(),
        &coefficient,
        &at,
        "source-reference-temperature",
        "source-range-end-temperature",
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert!((strain.engineering_strain() - 0.001752).abs() < 1e-14);
    for bad in [
        temperature(at.clone(), "source-reference-temperature", 273.15),
        temperature(at.clone(), "source-range-end-temperature", 333.15),
        temperature(at.clone(), "source-range-end-temperature", 373.16),
        replace_flag(at.clone(), "source-alloy-sac305", 0.0),
    ] {
        refuses(&pack, coefficient.name(), coefficient.dims(), &bad);
    }
    refuses(
        &pack,
        "linear-thermal-expansion-coefficient",
        EXPANSION,
        &at,
    );
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
                    "source-alloy-sac305",
                ],
            },
            &temperature(unknown, "temperature", 293.15),
            [293.15, 294.15],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .is_err()
    );
}

#[test]
fn sn63pb37_transport_uses_the_85c_column_and_retains_eutectic_endpoints() {
    let pack = compile("solder-sn63pb37-indium-ribbon-reference", 5);
    let unknown = unknown_temperature("source-alloy-sn63pb37");
    phase_endpoint(&pack, "solidus-temperature", &unknown, 183.0 + 273.15);
    phase_endpoint(&pack, "liquidus-temperature", &unknown, 183.0 + 273.15);
    sample(&pack, "density", RHO, &unknown, 8.39 * 1000.0);
    let psi = 0.45359237 * 9.80665 / 0.0254_f64.powi(2);
    sample(&pack, "tensile-strength", PA, &unknown, 8720.0 * psi);
    let at = temperature(reference("source-alloy-sn63pb37"), "temperature", 358.15);
    sample(&pack, "thermal-conductivity", K, &at, 0.48 * 100.0);
    let material = card(&pack);
    let slab = fs_conduction::ThermalResistanceTerm::slab_from_card(
        "Sn63Pb37 reference conductivity at 85 degC",
        0.001,
        1e-4,
        material.card(),
        &key("thermal-conductivity", K),
        &at,
        SelectionPolicy::SingleClaimOnly,
    )
    .unwrap();
    assert!((slab.value_k_per_w() - 0.001 / (48.0 * 1e-4)).abs() < 1e-14);
    for bad in [
        reference("source-alloy-sn63pb37"),
        temperature(at.clone(), "temperature", 293.15),
        temperature(at.clone(), "temperature", 358.16),
        temperature(reference("source-alloy-sac305"), "temperature", 358.15),
        replace_flag(at.clone(), "source-alloy-sn63pb37", 0.0),
    ] {
        refuses(&pack, "thermal-conductivity", K, &bad);
    }
    refuses(&pack, "young_modulus", PA, &unknown);
    refuses(&pack, "latent-heat-fusion", SPECIFIC_ENERGY, &unknown);
}

fn dsc_point(alloy_axis: &str, mass_kg: f64, ramp_k_s: f64) -> QueryPoint {
    temperature(
        flags(&[
            ("source-dsc-131", 1.0),
            ("source-nitrogen-atmosphere", 1.0),
            (alloy_axis, 1.0),
            ("source-aging-history-known", 0.0),
            ("source-pressure-known", 0.0),
        ]),
        "source-maximum-program-temperature",
        523.15,
    )
    .with_quantity(
        "source-specimen-mass",
        QuantitySpec::dimensional(MASS),
        mass_kg,
    )
    .unwrap()
    .with_quantity(
        "temperature-ramp-rate",
        QuantitySpec::dimensional(TEMPERATURE_RATE),
        ramp_k_s,
    )
    .unwrap()
}

#[test]
fn dsc_latent_energies_preserve_alloy_mass_and_heating_cooling_paths() {
    for (slug, alloy_axis, heating_j_g, cooling_j_g) in [
        (
            "solder-sac305-fros-dsc-2022",
            "source-alloy-sac305",
            [64.8, 65.1, 65.8],
            [62.4, 62.6, 62.9],
        ),
        (
            "solder-sn63pb37-fros-dsc-2022",
            "source-alloy-sn63pb37",
            [44.1, 42.6, 42.6],
            [39.7, 37.9, 38.9],
        ),
    ] {
        let pack = compile(slug, 2);
        for (index, mass_kg) in [7e-6, 14e-6, 21e-6].into_iter().enumerate() {
            let heating = sample(
                &pack,
                "latent-heat-fusion",
                SPECIFIC_ENERGY,
                &dsc_point(alloy_axis, mass_kg, 3.0 / 60.0),
                heating_j_g[index] * 1000.0,
            );
            let cooling = sample(
                &pack,
                "latent-heat-solidification",
                SPECIFIC_ENERGY,
                &dsc_point(alloy_axis, mass_kg, -3.0 / 60.0),
                cooling_j_g[index] * 1000.0,
            );
            for answer in [&heating, &cooling] {
                assert_eq!(
                    answer.evidence.value.uncertainty,
                    UncertaintyModel::Unstated
                );
            }
            // Integrated event energies, not a reversible h(T) or a heat-release rate.
            let absorbed_j = mass_kg * heating.evidence.value.value;
            let released_j = mass_kg * cooling.evidence.value.value;
            assert!(absorbed_j > released_j);
            assert!((absorbed_j - mass_kg * heating_j_g[index] * 1000.0).abs() < 1e-14);
            assert!((released_j - mass_kg * cooling_j_g[index] * 1000.0).abs() < 1e-14);
        }
        for (property, rate) in [
            ("latent-heat-fusion", 0.05),
            ("latent-heat-solidification", -0.05),
        ] {
            let at = dsc_point(alloy_axis, 14e-6, rate);
            for bad in [
                QueryPoint::new(),
                dsc_point(alloy_axis, 6e-6, rate),
                dsc_point(alloy_axis, 10.5e-6, rate),
                dsc_point(alloy_axis, 22e-6, rate),
                dsc_point(alloy_axis, 14e-6, -rate),
                dsc_point(alloy_axis, 14e-6, 2.0 * rate),
                replace_flag(at.clone(), alloy_axis, 0.0),
                replace_flag(at.clone(), "source-nitrogen-atmosphere", 0.0),
                replace_flag(at.clone(), "source-pressure-known", 1.0),
                replace_flag(at.clone(), "source-aging-history-known", 1.0),
                temperature(at.clone(), "source-maximum-program-temperature", 523.16),
            ] {
                refuses(&pack, property, SPECIFIC_ENERGY, &bad);
            }
            // Missing sensible branches and datum cannot become a phase trajectory.
            refuses(&pack, "specific-heat-capacity", CP, &at);
            refuses(&pack, "density", RHO, &at);
            refuses(&pack, "specific-enthalpy", SPECIFIC_ENERGY, &at);
        }
    }
}

#[test]
fn fusible_alloy_references_keep_exact_composition_and_conductivity_temperature() {
    let psi = 0.45359237 * 9.80665 / 0.0254_f64.powi(2);
    for (slug, alloy_axis, melting_c, density_g_cc, tensile_psi, conductivity, test_c) in [
        (
            "solder-sn42bi58-indium-ribbon-reference",
            "source-alloy-sn42bi58",
            138.0,
            8.56,
            11_700.0,
            16.0,
            20.0,
        ),
        (
            "solder-in52sn48-indium-ribbon-reference",
            "source-alloy-in52sn48",
            118.0,
            7.29,
            3_100.0,
            41.0,
            85.0,
        ),
        (
            "solder-sn96-5ag3-5-indium-ribbon-reference",
            "source-alloy-sn96-5ag3-5",
            221.0,
            7.36,
            5_690.0,
            33.0,
            85.0,
        ),
    ] {
        let pack = compile(slug, 5);
        let unknown = unknown_temperature(alloy_axis);
        phase_endpoint(&pack, "solidus-temperature", &unknown, melting_c + 273.15);
        phase_endpoint(&pack, "liquidus-temperature", &unknown, melting_c + 273.15);
        for (name, dims, expected) in [
            ("density", RHO, density_g_cc * 1000.0),
            ("tensile-strength", PA, tensile_psi * psi),
        ] {
            let answer = sample(&pack, name, dims, &unknown, expected);
            assert_eq!(
                answer.evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
            refuses(
                &pack,
                name,
                dims,
                &unknown_temperature("source-alloy-sac305"),
            );
            refuses(
                &pack,
                name,
                dims,
                &replace_flag(unknown.clone(), "source-test-temperature-known", 1.0),
            );
        }
        let at = temperature(reference(alloy_axis), "temperature", test_c + 273.15);
        sample(&pack, "thermal-conductivity", K, &at, conductivity);
        let material = card(&pack);
        let slab = fs_conduction::ThermalResistanceTerm::slab_from_card(
            "source-alloy conductivity reference slab",
            0.001,
            1e-4,
            material.card(),
            &key("thermal-conductivity", K),
            &at,
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
        assert!((slab.value_k_per_w() - 0.001 / (conductivity * 1e-4)).abs() < 1e-14);
        for bad in [
            QueryPoint::new(),
            reference(alloy_axis),
            temperature(at.clone(), "temperature", test_c + 273.16),
            temperature(at.clone(), "temperature", 313.15),
            replace_flag(at.clone(), alloy_axis, 0.0),
            replace_flag(at.clone(), "source-test-method-known", 1.0),
            replace_flag(at.clone(), "source-specimens-matched", 1.0),
        ] {
            refuses(&pack, "thermal-conductivity", K, &bad);
        }
        refuses(&pack, "latent-heat-fusion", SPECIFIC_ENERGY, &at);
        refuses(&pack, "young_modulus", PA, &unknown);
    }
}

#[test]
fn jining_granite_baseline_retains_unresolved_method_and_thermal_history() {
    let pack = compile("granite-jining-chen-yin-2024", 5);
    let at = flags(&[
        ("engineering-reference-model", 1.0),
        ("source-granite-jining", 1.0),
        ("source-before-thermal-disturbance", 1.0),
        ("source-specimens-matched", 0.0),
        ("source-test-temperature-known", 0.0),
        ("source-test-method-known", 0.0),
        ("source-moisture-content-known", 0.0),
    ]);
    for (name, dims, value) in [
        ("density", RHO, 2.6 * 1000.0),
        ("p-wave-speed", SPEED, 3054.85),
        ("elastic-modulus-method-unspecified", PA, 21.73e9),
        ("poisson-ratio-method-unspecified", Dims::NONE, 0.24),
        ("compressive-strength-method-unspecified", PA, 156.27e6),
    ] {
        let answer = sample(&pack, name, dims, &at, value);
        assert_eq!(
            answer.evidence.value.uncertainty,
            UncertaintyModel::Unstated
        );
        for bad in [
            QueryPoint::new(),
            replace_flag(at.clone(), "source-granite-jining", 0.0),
            replace_flag(at.clone(), "source-before-thermal-disturbance", 0.0),
            replace_flag(at.clone(), "source-specimens-matched", 1.0),
            replace_flag(at.clone(), "source-test-temperature-known", 1.0),
            replace_flag(at.clone(), "source-test-method-known", 1.0),
            replace_flag(at.clone(), "source-moisture-content-known", 1.0),
        ] {
            refuses(&pack, name, dims, &bad);
        }
    }

    // The source's Brazilian-test heating groups are not Table 1 measurement
    // temperatures. Its unspecified E cannot silently supply a static or
    // dynamic modulus, and wave speed cannot supply missing elastic constants.
    for (name, dims) in [
        ("young_modulus", PA),
        ("young-modulus-static", PA),
        ("young-modulus-dynamic", PA),
        ("shear_modulus", PA),
        ("bulk-modulus", PA),
        ("poisson_ratio", Dims::NONE),
        ("uniaxial-compressive-strength", PA),
        ("tensile-strength", PA),
        ("thermal-conductivity", K),
        ("specific-heat-capacity", CP),
    ] {
        refuses(&pack, name, dims, &at);
    }
}
