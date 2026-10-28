//! Source conditions and native consumers for the engineering-material additions.
use super::*;

const RESISTIVITY: Dims = Dims([3, 1, -3, 0, -2, 0]);
const EXPANSION: Dims = Dims([0, 0, 0, -1, 0, 0]);

#[test]
fn g0_g3_precious_metals_keep_product_and_temper_conditions() {
    for (slug, count, selector, rho, cp, k) in [
        (
            "indium-wire-indium-corporation-reference",
            6,
            "source-product-indium-wire-4n",
            7310.0,
            233.0,
            86.0,
        ),
        (
            "palladium-deringer-ney-reference",
            9,
            "source-product-ney-palladium",
            12020.0,
            250.0,
            70.0,
        ),
        (
            "platinum-heraeus-reference",
            11,
            "source-product-heraeus-platinum",
            21450.0,
            133.0,
            74.0,
        ),
    ] {
        let pack = compile(slug, count);
        let at = flags(&[
            (selector, 1.0),
            ("source-temperature-known", 0.0),
            ("source-pressure-known", 0.0),
        ]);
        for (property, dims, value) in [
            ("density", RHO, rho),
            ("specific-heat-capacity", CP, cp),
            ("thermal-conductivity", K, k),
        ] {
            sample(&pack, property, dims, &at, value);
            for wrong in [
                QueryPoint::new(),
                at.clone().with(selector, 0.0).unwrap(),
                at.clone().with("source-temperature-known", 1.0).unwrap(),
                at.clone().with("source-pressure-known", 1.0).unwrap(),
            ] {
                refuses(&pack, property, dims, &wrong);
            }
        }
        if slug.starts_with("palladium") {
            sample(&pack, "electrical-resistivity", RESISTIVITY, &at, 1.08e-7);
            let wire = at
                .clone()
                .with_quantity(
                    "annealed-fine-wire",
                    QuantitySpec::dimensional(Dims::NONE),
                    1.0,
                )
                .unwrap();
            sample(&pack, "yield-strength", PA, &wire, 270e6);
            sample(&pack, "ultimate-tensile-strength", PA, &wire, 280e6);
            sample(&pack, "elongation-at-break", Dims::NONE, &wire, 0.25);
            refuses(&pack, "yield-strength", PA, &at);
            refuses(
                &pack,
                "yield-strength",
                PA,
                &wire.with("annealed-fine-wire", 0.0).unwrap(),
            );
        }
        if slug.starts_with("platinum") {
            // Work hardening changes strength and ductility, not a generic
            // material label that may silently select either source state.
            for (hardened, yield_pa, ultimate_pa, elongation) in
                [(0.0, 70e6, 150e6, 0.40), (1.0, 290e6, 330e6, 0.03)]
            {
                let temper = at
                    .clone()
                    .with_quantity(
                        "work-hardened",
                        QuantitySpec::dimensional(Dims::NONE),
                        hardened,
                    )
                    .unwrap();
                sample(&pack, "yield-strength", PA, &temper, yield_pa);
                sample(&pack, "ultimate-tensile-strength", PA, &temper, ultimate_pa);
                sample(
                    &pack,
                    "elongation-at-break",
                    Dims::NONE,
                    &temper,
                    elongation,
                );
            }
            refuses(&pack, "yield-strength", PA, &at);
            refuses(
                &pack,
                "yield-strength",
                PA,
                &at.with_quantity("work-hardened", QuantitySpec::dimensional(Dims::NONE), 0.5)
                    .unwrap(),
            );
        }
    }
}

#[test]
fn g0_g3_iridium_point_observations_do_not_become_an_invented_curve() {
    let pack = compile("iridium-espi-reference", 7);
    let reference = flags(&[
        ("source-reference-espi-iridium", 1.0),
        ("source-pressure-known", 0.0),
    ]);
    for (kelvin, resistivity) in [(273.15, 4.71e-8), (293.15, 5.3e-8)] {
        let at = temperature(reference.clone(), "temperature", kelvin);
        sample(
            &pack,
            "electrical-resistivity",
            RESISTIVITY,
            &at,
            resistivity,
        );
    }
    for kelvin in [273.14, 283.15, 293.16] {
        refuses(
            &pack,
            "electrical-resistivity",
            RESISTIVITY,
            &temperature(reference.clone(), "temperature", kelvin),
        );
    }
    let at = temperature(reference.clone(), "temperature", 293.15);
    sample(&pack, "density", RHO, &at, 22650.0);
    sample(
        &pack,
        "linear-thermal-expansion-coefficient",
        EXPANSION,
        &at,
        6.8e-6,
    );
    refuses(&pack, "density", RHO, &reference);
    refuses(
        &pack,
        "density",
        RHO,
        &at.clone().with("temperature", 293.16).unwrap(),
    );
    let unknown = reference
        .clone()
        .with_quantity(
            "source-temperature-known",
            QuantitySpec::dimensional(Dims::NONE),
            0.0,
        )
        .unwrap();
    sample(&pack, "specific-heat-capacity", CP, &unknown, 130.0);
    refuses(&pack, "specific-heat-capacity", CP, &reference);
    // Only conductivity is reported for the source's 0--100 C interval.
    for kelvin in [273.15, 323.15, 373.15] {
        sample(
            &pack,
            "thermal-conductivity",
            K,
            &temperature(reference.clone(), "temperature", kelvin),
            147.0,
        );
    }
    for kelvin in [273.14, 373.16] {
        refuses(
            &pack,
            "thermal-conductivity",
            K,
            &temperature(reference.clone(), "temperature", kelvin),
        );
    }
}

#[test]
fn g0_g3_iridium_conductivity_drives_a_bounded_native_slab() {
    use fs_conduction::{SlabConductivityQuery, ThermalResistanceTerm};
    use fs_matdb::{MaterialStateId, NormalizedMaterialCardPack};

    let pack = compile("iridium-espi-reference", 7);
    let normalized = NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: pack.pack_id().into(),
            phase: "solid iridium".into(),
            process: "ESPI general conductivity reference over 0--100 C".into(),
            revision: 0,
        },
        pack,
    )
    .unwrap();
    let card = normalized.card();
    let conductivity = key("thermal-conductivity", K);
    let reference = flags(&[
        ("source-reference-espi-iridium", 1.0),
        ("source-pressure-known", 0.0),
    ]);
    let lower = temperature(reference.clone(), "temperature", 273.15);
    let upper = temperature(reference, "temperature", 373.15);
    let term = ThermalResistanceTerm::slab_over_temperature_span_from_card(
        "iridium reference wall",
        0.002,
        0.01,
        card,
        SlabConductivityQuery {
            property: &conductivity,
            temperature_axis: "temperature",
            lower: &lower,
            upper: &upper,
            selection: SelectionPolicy::SingleClaimOnly,
        },
    )
    .unwrap();
    assert!((term.value_k_per_w() - 0.002 / (0.01 * 147.0)).abs() < 1e-14);
    assert!(
        ThermalResistanceTerm::slab_over_temperature_span_from_card(
            "unsupported hotter iridium wall",
            0.002,
            0.01,
            card,
            SlabConductivityQuery {
                property: &conductivity,
                temperature_axis: "temperature",
                lower: &lower,
                upper: &upper.with("temperature", 373.16).unwrap(),
                selection: SelectionPolicy::SingleClaimOnly,
            },
        )
        .is_err()
    );
}

#[test]
fn g0_g3_concrete_fresh_density_and_hardened_strength_remain_separate() {
    let pack = compile("quikrete-5000-concrete-1007-2022", 5);
    let fresh = flags(&[
        ("source-product-quikrete-5000-1007", 1.0),
        ("source-temperature-known", 0.0),
        ("source-pressure-known", 0.0),
        ("fresh-mixture", 1.0),
    ]);
    sample(&pack, "density", RHO, &fresh, 2242.5);
    let hardened = fresh.clone().with("fresh-mixture", 0.0).unwrap();
    refuses(&pack, "density", RHO, &hardened);
    let at_age = |days: f64| {
        hardened
            .clone()
            .with_quantity(
                "reported-age",
                QuantitySpec::dimensional(Dims([0, 0, 1, 0, 0, 0])),
                days * 86400.0,
            )
            .unwrap()
    };
    for (days, mpa) in [(1.0, 10.3), (3.0, 17.2), (7.0, 24.1), (28.0, 34.4)] {
        let at = at_age(days);
        sample(&pack, "compression-strength", PA, &at, mpa * 1e6);
        refuses(&pack, "density", RHO, &at);
        refuses(
            &pack,
            "compression-strength",
            PA,
            &at.with("fresh-mixture", 1.0).unwrap(),
        );
    }
    for days in [0.0, 2.0, 14.0, 29.0] {
        refuses(&pack, "compression-strength", PA, &at_age(days));
    }
    refuses(&pack, "compression-strength", PA, &hardened);
    refuses(&pack, "young-modulus", PA, &at_age(28.0));
}

#[test]
fn g0_g3_repair_cement_preserves_age_origins_flow_and_test_mode() {
    let pack = compile("cts-cement-all-ds024-2026", 11);
    let reference = temperature(
        flags(&[
            ("source-product-cement-all-ds024-2026", 1.0),
            ("flow-consistency", 102.0),
            ("source-pressure-known", 0.0),
        ]),
        "temperature",
        294.15,
    );
    for (property, seconds) in [
        ("initial-setting-time", 900.0),
        ("final-setting-time", 2100.0),
    ] {
        sample(
            &pack,
            property,
            Dims([0, 0, 1, 0, 0, 0]),
            &reference,
            seconds,
        );
    }
    let after_final_set = reference
        .clone()
        .with_quantity(
            "source-age-origin-known",
            QuantitySpec::dimensional(Dims::NONE),
            1.0,
        )
        .unwrap()
        .with_quantity(
            "age-after-final-set",
            QuantitySpec::dimensional(Dims([0, 0, 1, 0, 0, 0])),
            3600.0,
        )
        .unwrap();
    sample(&pack, "compression-strength", PA, &after_final_set, 20.7e6);
    refuses(
        &pack,
        "compression-strength",
        PA,
        &after_final_set
            .clone()
            .with("age-after-final-set", 3599.0)
            .unwrap(),
    );
    let at_age = |seconds| {
        reference
            .clone()
            .with_quantity(
                "source-age-origin-known",
                QuantitySpec::dimensional(Dims::NONE),
                0.0,
            )
            .unwrap()
            .with_quantity(
                "reported-age",
                QuantitySpec::dimensional(Dims([0, 0, 1, 0, 0, 0])),
                seconds,
            )
            .unwrap()
    };
    for (seconds, mpa) in [
        (10800.0, 34.5),
        (86400.0, 41.4),
        (604800.0, 48.3),
        (2419200.0, 62.1),
    ] {
        sample(
            &pack,
            "compression-strength",
            PA,
            &at_age(seconds),
            mpa * 1e6,
        );
    }
    // One hour since an unspecified origin is not one hour after final set.
    refuses(&pack, "compression-strength", PA, &at_age(3600.0));
    refuses(&pack, "compression-strength", PA, &at_age(14.0 * 86400.0));
    for (days, flexural, splitting) in [(7.0, 4.14e6, 4.82e6), (28.0, 5.51e6, 6.06e6)] {
        let at = at_age(days * 86400.0);
        sample(&pack, "flexural-strength", PA, &at, flexural);
        sample(&pack, "splitting-tensile-strength", PA, &at, splitting);
        refuses(&pack, "ultimate-tensile-strength", PA, &at);
        for bad in [
            at.clone().with("temperature", 294.16).unwrap(),
            at.clone().with("flow-consistency", 103.0).unwrap(),
            at.clone().with("source-age-origin-known", 1.0).unwrap(),
        ] {
            refuses(&pack, "compression-strength", PA, &bad);
        }
    }
}
