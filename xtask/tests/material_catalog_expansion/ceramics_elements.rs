//! Source transport, exact applicability and reduced physical use of new
//! elemental, stainless-steel, mineral-crystal and engineering-ceramic packs.
//! These check authored reference models, not measured coupon validation.

use super::*;

const RESISTIVITY: Dims = Dims([3, 1, -3, 0, -2, 0]);
const INVERSE_TEMPERATURE: Dims = Dims([0, 0, 0, -1, 0, 0]);
const FREQUENCY: Dims = Dims([0, 0, -1, 0, 0, 0]);

fn frequency(point: QueryPoint, hz: f64) -> QueryPoint {
    point
        .with_quantity("frequency", QuantitySpec::dimensional(FREQUENCY), hz)
        .unwrap()
}

fn selector(point: QueryPoint, name: &str, value: f64) -> QueryPoint {
    point
        .with_quantity(name, QuantitySpec::dimensional(Dims::NONE), value)
        .unwrap()
}

fn reference_card(pack: &NormalizedPack) -> fs_matdb::NormalizedMaterialCardPack {
    fs_matdb::NormalizedMaterialCardPack::new(
        fs_matdb::MaterialStateId {
            chemistry: pack.pack_id().into(),
            phase: "solid".into(),
            process: "original source grade, temperature and explicit unknown-condition flags"
                .into(),
            revision: 0,
        },
        pack.clone(),
    )
    .unwrap()
}

/// Only a frozen conductivity at the exact accepted source point is used.
/// This does not widen a singleton fact into a temperature-span law.
fn slab_reference(pack: &NormalizedPack, at: &QueryPoint, conductivity: f64) {
    let normalized = reference_card(pack);
    let slab = fs_conduction::ThermalResistanceTerm::slab_from_card(
        "20 mm source-reference slab",
        0.02,
        0.004,
        normalized.card(),
        &key("thermal-conductivity", K),
        at,
        SelectionPolicy::SingleClaimOnly,
    )
    .unwrap();
    assert!((slab.value_k_per_w() / (5.0 / conductivity) - 1.0).abs() < 1e-13);
}

/// A source Young modulus drives the existing axial owner at a caller-authored
/// small strain. No tensile-test protocol, failure limit or isotropy is inferred.
fn axial_reference(pack: &NormalizedPack, at: &QueryPoint, density: f64, modulus: f64) {
    use fs_material::axial_reference::{
        AxialLoadingMode, AxialReferenceRequest, ResolvedAxialModulusReference,
    };
    use fs_material::state_point::{
        MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    };
    let normalized = reference_card(pack);
    let request = AxialReferenceRequest {
        density: ScalarPropertyRequirement::try_new(
            "density",
            RHO,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
        modulus: ScalarPropertyRequirement::try_new(
            "young-modulus",
            PA,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
        mode: AxialLoadingMode::Tension,
        maximum_absolute_strain: 0.0002,
        crosshead_speed_axis: None,
    };
    let material = ResolvedAxialModulusReference::resolve(
        normalized.card(),
        at,
        &request,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let mut coupon = material.coupon(0.1, 1e-4).unwrap();
    let response = coupon.advance_displacement(1e-5, 1.0).unwrap();
    let force = modulus * 1e-4 * (1e-5 / 0.1);
    assert!((response.reference_force_n / force - 1.0).abs() < 1e-13);
    assert!((response.mass_kg / (density * 1e-5) - 1.0).abs() < 1e-13);
    assert!((response.model_loading_work_j / (0.5 * force * 1e-5) - 1.0).abs() < 1e-13);
    for property in material.source().properties() {
        normalized
            .card()
            .claims()
            .verify_receipt(&property.answer().receipt)
            .unwrap();
    }
    assert!(coupon.advance_displacement(3e-5, 1.0).is_err());
}

fn mean_expansion_reference(
    pack: &NormalizedPack,
    at: &QueryPoint,
    name: &str,
    reference: f64,
    endpoint: f64,
    alpha: f64,
) {
    use fs_material::state_point::{
        MaterialPropertySelection, resolve_mean_thermal_expansion_strain,
    };
    let normalized = reference_card(pack);
    let result = resolve_mean_thermal_expansion_strain(
        normalized.card(),
        &key(name, INVERSE_TEMPERATURE),
        at,
        "source-reference-temperature",
        "source-range-end-temperature",
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert!((result.engineering_strain() - alpha * (endpoint - reference)).abs() < 1e-14);
    assert_eq!(result.reference_temperature_k(), reference);
    assert_eq!(result.endpoint_temperature_k(), endpoint);
}

#[test]
fn g0_g3_elemental_bismuth_and_vanadium_preserve_units_and_distinct_temperatures() {
    for (slug, rho, cp, k, resistivity, resistivity_temperature) in [
        (
            "bismuth-fine-metals-reference",
            9747.0,
            123.8464,
            7.92,
            1.068e-6,
            273.15,
        ),
        (
            "vanadium-espi-reference",
            6110.0,
            485.344,
            30.7,
            2.54e-7,
            293.15,
        ),
    ] {
        let pack = compile(slug, 4);
        let base = flags(&[("source-process-known", 0.0)]);
        let unknown_density = selector(base.clone(), "source-density-temperature-known", 0.0);
        sample(&pack, "density", RHO, &unknown_density, rho);
        refuses(&pack, "density", RHO, &base);
        let conductivity = temperature(base.clone(), "temperature", 298.2);
        sample(&pack, "thermal-conductivity", K, &conductivity, k);
        slab_reference(&pack, &conductivity, k);
        let heat_capacity = temperature(
            selector(base.clone(), "thermochemical-calorie-conversion", 1.0),
            "temperature",
            298.15,
        );
        sample(&pack, "specific-heat-capacity", CP, &heat_capacity, cp);
        refuses(&pack, "specific-heat-capacity", CP, &conductivity);
        refuses(&pack, "thermal-conductivity", K, &heat_capacity);
        sample(
            &pack,
            "electrical-resistivity",
            RESISTIVITY,
            &temperature(base.clone(), "temperature", resistivity_temperature),
            resistivity,
        );
        refuses(
            &pack,
            "specific-heat-capacity",
            CP,
            &selector(heat_capacity, "thermochemical-calorie-conversion", 0.0),
        );
        refuses(&pack, "young-modulus", PA, &base);
    }
}

#[test]
fn g0_g3_430_and_316l_keep_grade_and_interval_conditions_in_physical_use() {
    for (slug, grade, rho, k, cp, resistivity, e, e_mid, alpha) in [
        (
            "stainless-14016-thyssenkrupp-reference",
            "source-grade-stainless-14016",
            7700.0,
            25.0,
            460.0,
            6e-7,
            220e9,
            215e9,
            11e-6,
        ),
        (
            "stainless-14404-thyssenkrupp-reference",
            "source-grade-stainless-14404",
            8000.0,
            15.0,
            500.0,
            7.5e-7,
            200e9,
            193e9,
            16.5e-6,
        ),
    ] {
        let pack = compile(slug, 6);
        let base = flags(&[
            (grade, 1.0),
            ("source-heat-treatment-known", 0.0),
            ("source-pressure-known", 0.0),
        ]);
        let at = temperature(base.clone(), "temperature", 293.15);
        for (name, dims, value) in [
            ("density", RHO, rho),
            ("specific-heat-capacity", CP, cp),
            ("thermal-conductivity", K, k),
            ("electrical-resistivity", RESISTIVITY, resistivity),
            ("young-modulus", PA, e),
        ] {
            sample(&pack, name, dims, &at, value);
            refuses(&pack, name, dims, &selector(at.clone(), grade, 0.0));
        }
        axial_reference(&pack, &at, rho, e);
        slab_reference(&pack, &at, k);
        sample(
            &pack,
            "young-modulus",
            PA,
            // The midpoint of the 20 and 200 degC source knots is 110 degC.
            &temperature(base.clone(), "temperature", 383.15),
            e_mid,
        );
        refuses(
            &pack,
            "density",
            RHO,
            &temperature(base.clone(), "temperature", 383.15),
        );
        refuses(
            &pack,
            "young-modulus",
            PA,
            &temperature(base.clone(), "temperature", 800.0),
        );
        let mean = temperature(
            temperature(base.clone(), "source-reference-temperature", 293.15),
            "source-range-end-temperature",
            473.15,
        );
        let name = "mean-linear-expansion-coefficient-from-20c";
        sample(&pack, name, INVERSE_TEMPERATURE, &mean, alpha);
        mean_expansion_reference(&pack, &mean, name, 293.15, 473.15, alpha);
        refuses(
            &pack,
            name,
            INVERSE_TEMPERATURE,
            &temperature(mean.clone(), "source-range-end-temperature", 423.15),
        );
        refuses(
            &pack,
            name,
            INVERSE_TEMPERATURE,
            &temperature(mean, "source-reference-temperature", 298.15),
        );
        refuses(&pack, "tensile-strength", PA, &at);
    }
}

#[test]
fn g0_g3_alumina_grades_keep_20c_mechanics_100c_heat_capacity_and_dielectric_frequency() {
    for (slug, grade, rho, e, nu, flexure, compression, k, eps, loss, resistivity) in [
        (
            "alumina-ad96-coorstek-reference",
            "source-grade-alumina-ad96",
            3720.0,
            303e9,
            0.21,
            343e6,
            2068e6,
            24.7,
            9.0,
            0.0002,
            4e7,
        ),
        (
            "alumina-ad995-coorstek-reference",
            "source-grade-alumina-ad995",
            3900.0,
            370e9,
            0.22,
            300e6,
            2600e6,
            30.0,
            9.7,
            0.0001,
            2e8,
        ),
    ] {
        let pack = compile(slug, 11);
        let base = flags(&[(grade, 1.0), ("source-process-known", 0.0)]);
        let at = temperature(base.clone(), "temperature", 293.15);
        for (name, dims, value) in [
            ("young-modulus", PA, e),
            ("poisson-ratio", Dims::NONE, nu),
            ("four-point-flexural-strength", PA, flexure),
            ("compressive-strength", PA, compression),
            ("thermal-conductivity", K, k),
        ] {
            sample(&pack, name, dims, &at, value);
        }
        let density_at = selector(at.clone(), "source-density-temperature-known", 0.0);
        sample(&pack, "density", RHO, &density_at, rho);
        axial_reference(&pack, &density_at, rho, e);
        slab_reference(&pack, &at, k);
        refuses(&pack, "density", RHO, &at);
        refuses(&pack, "specific-heat-capacity", CP, &at);
        let hot = temperature(base.clone(), "temperature", 373.15);
        sample(&pack, "specific-heat-capacity", CP, &hot, 880.0);
        refuses(&pack, "thermal-conductivity", K, &hot);
        let dielectric = frequency(temperature(base.clone(), "temperature", 298.15), 1e6);
        sample(&pack, "relative-permittivity", Dims::NONE, &dielectric, eps);
        sample(
            &pack,
            "dielectric-loss-tangent",
            Dims::NONE,
            &dielectric,
            loss,
        );
        refuses(
            &pack,
            "relative-permittivity",
            Dims::NONE,
            &frequency(at.clone(), 1e6),
        );
        refuses(
            &pack,
            "relative-permittivity",
            Dims::NONE,
            &frequency(temperature(base.clone(), "temperature", 298.15), 1000.0),
        );
        sample(
            &pack,
            "electrical-resistivity",
            RESISTIVITY,
            &temperature(base.clone(), "temperature", 773.15),
            resistivity,
        );
        for t in [298.15, 1000.0] {
            refuses(
                &pack,
                "electrical-resistivity",
                RESISTIVITY,
                &temperature(base.clone(), "temperature", t),
            );
        }
        let mean = temperature(
            temperature(base.clone(), "source-reference-temperature", 298.15),
            "source-range-end-temperature",
            1273.15,
        );
        let name = "mean-linear-expansion-coefficient-from-25c";
        sample(&pack, name, INVERSE_TEMPERATURE, &mean, 8.2e-6);
        mean_expansion_reference(&pack, &mean, name, 298.15, 1273.15, 8.2e-6);
        refuses(
            &pack,
            name,
            INVERSE_TEMPERATURE,
            &temperature(mean, "source-range-end-temperature", 873.15),
        );
        refuses(&pack, "tensile-strength", PA, &at);
    }
}

#[test]
fn g0_g3_aluminum_nitride_exposes_high_conductivity_without_inventing_joint_thermal_support() {
    let pack = compile("aluminum-nitride-180-coorstek-reference", 8);
    let base = flags(&[
        ("source-grade-aluminum-nitride-180", 1.0),
        ("source-process-known", 0.0),
    ]);
    let at = temperature(base.clone(), "temperature", 293.15);
    for (name, dims, value) in [
        ("young-modulus", PA, 310e9),
        ("four-point-flexural-strength", PA, 300e6),
        ("compressive-strength", PA, 2000e6),
        ("thermal-conductivity", K, 180.0),
    ] {
        sample(&pack, name, dims, &at, value);
    }
    let density_at = selector(at.clone(), "source-density-temperature-known", 0.0);
    sample(&pack, "density", RHO, &density_at, 3300.0);
    axial_reference(&pack, &density_at, 3300.0, 310e9);
    slab_reference(&pack, &at, 180.0);
    sample(
        &pack,
        "specific-heat-capacity",
        CP,
        &temperature(base.clone(), "temperature", 373.15),
        738.0,
    );
    refuses(&pack, "specific-heat-capacity", CP, &at);
    refuses(&pack, "poisson-ratio", Dims::NONE, &at);
    refuses(
        &pack,
        "dielectric-loss-tangent",
        Dims::NONE,
        &frequency(temperature(base.clone(), "temperature", 298.15), 1e6),
    );
    sample(
        &pack,
        "relative-permittivity",
        Dims::NONE,
        &frequency(temperature(base.clone(), "temperature", 298.15), 1e6),
        8.6,
    );
    let mean = temperature(
        temperature(base, "source-reference-temperature", 298.15),
        "source-range-end-temperature",
        1273.15,
    );
    sample(
        &pack,
        "mean-linear-expansion-coefficient-from-25c",
        INVERSE_TEMPERATURE,
        &mean,
        5.6e-6,
    );
    mean_expansion_reference(
        &pack,
        &mean,
        "mean-linear-expansion-coefficient-from-25c",
        298.15,
        1273.15,
        5.6e-6,
    );
}

#[test]
fn g0_g3_mineral_crystals_separate_electrical_source_facts_from_thermal_models() {
    for (slug, rho, cp, k, k_temperature, eps, dielectric_temperature) in [
        (
            "halite-nacl-crystran-reference",
            2170.0,
            854.0,
            6.5,
            273.0,
            5.9,
            None,
        ),
        (
            "sylvite-kcl-crystran-reference",
            1990.0,
            690.0,
            6.53,
            322.0,
            4.64,
            Some(300.0),
        ),
        (
            "periclase-mgo-crystran-reference",
            3580.0,
            877.0,
            42.0,
            273.0,
            9.65,
            None,
        ),
    ] {
        let pack = compile(slug, 4);
        let base = flags(&[("source-optical-crystal", 1.0)]);
        let unknown = selector(base.clone(), "source-test-temperature-known", 0.0);
        let conductivity = if slug == "sylvite-kcl-crystran-reference" {
            sample(&pack, "density", RHO, &unknown, rho);
            sample(&pack, "specific-heat-capacity", CP, &unknown, cp);
            refuses(&pack, "density", RHO, &base);
            let at = temperature(base.clone(), "temperature", k_temperature);
            refuses(
                &pack,
                "thermal-conductivity",
                K,
                &temperature(base.clone(), "temperature", k_temperature + 0.15),
            );
            refuses(
                &pack,
                "thermal-conductivity",
                K,
                &selector(at.clone(), "source-optical-crystal", 0.0),
            );
            at
        } else {
            // Preserve the existing explicitly authored frozen thermal model;
            // its one-kelvin interval does not qualify the new electrical fact.
            let at = temperature(
                flags(&[
                    ("engineering-reference-model", 1.0),
                    ("source-temperatures-matched", 0.0),
                ]),
                "temperature",
                273.5,
            );
            sample(&pack, "density", RHO, &at, rho);
            sample(&pack, "specific-heat-capacity", CP, &at, cp);
            refuses(
                &pack,
                "relative-permittivity",
                Dims::NONE,
                &frequency(at.clone(), 1e6),
            );
            refuses(
                &pack,
                "thermal-conductivity",
                K,
                &temperature(at.clone(), "temperature", 274.01),
            );
            at
        };
        sample(&pack, "thermal-conductivity", K, &conductivity, k);
        slab_reference(&pack, &conductivity, k);
        let dielectric = frequency(
            match dielectric_temperature {
                Some(t) => temperature(base.clone(), "temperature", t),
                None => unknown.clone(),
            },
            1e6,
        );
        sample(&pack, "relative-permittivity", Dims::NONE, &dielectric, eps);
        refuses(
            &pack,
            "relative-permittivity",
            Dims::NONE,
            &frequency(dielectric, 1000.0),
        );
        refuses(&pack, "young-modulus", PA, &unknown);
        refuses(&pack, "poisson-ratio", Dims::NONE, &unknown);
    }
}

#[test]
fn g0_g3_hexoloy_preserves_source_edition_room_temperature_and_flexure_geometry() {
    let pack = compile("silicon-carbide-hexoloy-sa-2012", 10);
    let rt = flags(&[
        ("source-grade-hexoloy-sa", 1.0),
        ("source-test-temperature-known", 0.0),
    ]);
    for (name, dims, value) in [
        ("density", RHO, 3100.0),
        ("young-modulus", PA, 430e9),
        ("poisson-ratio", Dims::NONE, 0.14),
        ("four-point-flexural-strength", PA, 380e6),
        ("three-point-flexural-strength", PA, 550e6),
        ("compressive-strength", PA, 3900e6),
        ("specific-heat-capacity", CP, 670.0),
        ("thermal-conductivity", K, 125.6),
    ] {
        sample(&pack, name, dims, &rt, value);
        refuses(&pack, name, dims, &QueryPoint::new());
    }
    axial_reference(&pack, &rt, 3100.0, 430e9);
    slab_reference(&pack, &rt, 125.6);
    let known = flags(&[
        ("source-grade-hexoloy-sa", 1.0),
        ("source-test-temperature-known", 1.0),
    ]);
    for (t, k) in [(473.15, 102.6), (573.15, 90.05), (673.15, 77.5)] {
        let at = temperature(known.clone(), "temperature", t);
        sample(&pack, "thermal-conductivity", K, &at, k);
        slab_reference(&pack, &at, k);
        refuses(&pack, "specific-heat-capacity", CP, &at);
    }
    for t in [293.15, 473.14, 673.16] {
        refuses(
            &pack,
            "thermal-conductivity",
            K,
            &temperature(known.clone(), "temperature", t),
        );
    }
    sample(
        &pack,
        "mean-linear-expansion-coefficient-from-room-temperature",
        INVERSE_TEMPERATURE,
        &temperature(rt.clone(), "source-range-end-temperature", 973.15),
        4.02e-6,
    );
    refuses(&pack, "tensile-strength", PA, &rt);
    refuses(&pack, "electrical-resistivity", RESISTIVITY, &rt);
}

#[test]
fn g0_g3_macor_preserves_interval_means_frequency_knots_and_minimum_average_strength() {
    let pack = compile("macor-corning-2013", 11);
    let base = flags(&[("source-grade-macor", 1.0)]);
    let at = temperature(base.clone(), "temperature", 298.15);
    for (name, dims, value) in [
        ("young-modulus", PA, 66.9e9),
        ("shear-modulus", PA, 25.5e9),
        ("flexural-strength-minimum-specified-average", PA, 94e6),
        ("thermal-conductivity", K, 1.46),
        ("specific-heat-capacity", CP, 790.0),
        ("electrical-resistivity", RESISTIVITY, 1e15),
    ] {
        sample(&pack, name, dims, &at, value);
        refuses(
            &pack,
            name,
            dims,
            &temperature(base.clone(), "temperature", 293.15),
        );
    }
    let unknown = selector(at.clone(), "source-test-temperature-known", 0.0);
    sample(&pack, "density", RHO, &unknown, 2520.0);
    sample(&pack, "poisson-ratio", Dims::NONE, &unknown, 0.29);
    axial_reference(&pack, &unknown, 2520.0, 66.9e9);
    slab_reference(&pack, &at, 1.46);
    refuses(&pack, "density", RHO, &at);
    refuses(&pack, "flexural-strength", PA, &at);
    refuses(&pack, "tensile-strength", PA, &at);
    for (hz, eps, loss) in [(1000.0, 6.01, 0.0040), (8.5e9, 5.64, 0.0025)] {
        let electrical = frequency(at.clone(), hz);
        sample(&pack, "relative-permittivity", Dims::NONE, &electrical, eps);
        sample(
            &pack,
            "dielectric-loss-tangent",
            Dims::NONE,
            &electrical,
            loss,
        );
    }
    refuses(
        &pack,
        "relative-permittivity",
        Dims::NONE,
        &frequency(at.clone(), 1e6),
    );
    for (end, alpha) in [(573.15, 9e-6), (873.15, 11.2e-6), (1073.15, 12.3e-6)] {
        let mean = temperature(
            temperature(base.clone(), "source-reference-temperature", 298.15),
            "source-range-end-temperature",
            end,
        );
        let name = "mean-linear-expansion-coefficient-from-25c";
        sample(&pack, name, INVERSE_TEMPERATURE, &mean, alpha);
        mean_expansion_reference(&pack, &mean, name, 298.15, end, alpha);
        refuses(
            &pack,
            name,
            INVERSE_TEMPERATURE,
            &temperature(mean, "source-range-end-temperature", end - 1.0),
        );
    }
}
