use super::*;

const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);
const EXPANSION: Dims = Dims([0, 0, 0, -1, 0, 0]);
const ELECTRIC_FIELD: Dims = Dims([1, 1, -3, 0, -1, 0]);

fn dim(point: QueryPoint, name: &str, dims: Dims, value: f64) -> QueryPoint {
    point
        .with_quantity(name, QuantitySpec::dimensional(dims), value)
        .unwrap()
}

#[test]
fn cured_epoxies_preserve_cure_history_and_hardener_selection() {
    let henkel = compile("henkel-loctite-ea9466-seven-day-film-2025", 5);
    let film = |cure_axis: &str, age, thickness, test_temperature_known| {
        let point = temperature(
            flags(&[
                ("source-grade-ea9466", 1.0),
                ("resin-to-hardener-mass-ratio", 2.0),
                ("source-moisture-condition-known", 0.0),
                ("source-test-temperature-known", test_temperature_known),
                ("source-test-rate-known", 0.0),
                ("source-thermal-scan-conditions-known", 0.0),
            ]),
            cure_axis,
            295.15,
        );
        dim(
            dim(point, "cure-age", TIME, age),
            "film-thickness",
            LENGTH,
            thickness,
        )
    };
    let point = film("cure-temperature", 604_800.0, 0.0012, 0.0);
    for (name, dims, expected) in [
        ("tensile-modulus", PA, 1.718e9),
        ("tensile-strength", PA, 32.0e6),
        ("tensile-elongation", Dims::NONE, 0.03),
        ("dielectric-breakdown-strength", ELECTRIC_FIELD, 30.0e6),
    ] {
        sample(&henkel, name, dims, &point, expected);
    }
    let transition = henkel
        .claims()
        .query_typed(
            &PropertyKey::with_quantity(
                "glass-transition-temperature",
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::AbsoluteTemperature,
                    ValueForm::Static,
                )),
            ),
            &point,
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
    assert_eq!(transition.evidence.value.value, 335.15);
    henkel.claims().verify_receipt(&transition.receipt).unwrap();
    // A cure temperature cannot silently become a measured test temperature.
    for bad in [
        film("temperature", 604_800.0, 0.0012, 0.0),
        film("cure-temperature", 604_800.0, 0.0012, 1.0),
        film("cure-temperature", 86_400.0, 0.0012, 0.0),
        film("cure-temperature", 604_800.0, 0.0024, 0.0),
    ] {
        refuses(&henkel, "tensile-modulus", PA, &bad);
    }

    let filled = compile("smooth-on-epoxacast655-103-seven-day-2024", 9);
    let epoxy = |hardener, parts, kelvin, age| {
        dim(
            temperature(
                flags(&[
                    ("source-grade-epoxacast655", 1.0),
                    ("hardener-103", hardener),
                    ("hardener-parts-per-100-resin-mass", parts),
                    ("source-processing-history-known", 0.0),
                    ("source-moisture-condition-known", 0.0),
                    ("source-test-rate-known", 0.0),
                    ("source-thermal-test-method-known", 0.0),
                ]),
                "temperature",
                kelvin,
            ),
            "cure-age",
            TIME,
            age,
        )
    };
    let point = epoxy(1.0, 12.0, 296.15, 604_800.0);
    let psi = 6_894.757_293_168;
    for (name, dims, expected) in [
        ("tensile-modulus", PA, 2_040_000.0 * psi),
        ("tensile-strength", PA, 4_810.0 * psi),
        ("tensile-elongation", Dims::NONE, 0.0054),
        ("flexural-modulus", PA, 1_403_000.0 * psi),
        ("flexural-strength", PA, 7_660.0 * psi),
        ("compression-modulus", PA, 125_500.0 * psi),
        ("compressive-strength", PA, 12_500.0 * psi),
        ("linear-molding-shrinkage", Dims::NONE, 0.0006),
        ("thermal-conductivity", K, 0.675),
    ] {
        sample(&filled, name, dims, &point, expected);
    }
    for bad in [
        epoxy(0.0, 12.0, 296.15, 604_800.0),
        epoxy(1.0, 10.0, 296.15, 604_800.0),
        epoxy(1.0, 12.0, 297.15, 604_800.0),
        epoxy(1.0, 12.0, 296.15, 86_400.0),
    ] {
        refuses(&filled, "tensile-modulus", PA, &bad);
    }
    // Mixed-liquid specific gravity does not supply cured density.
    refuses(&filled, "density", RHO, &point);
}

#[test]
fn ceramic_adhesive_and_melamine_glass_keep_units_and_coupon_directions() {
    let ceramic = compile("cotronics-resbond940-zircon-reference", 5);
    let adhesive = |btu_it, btu_known| {
        flags(&[
            ("source-grade-resbond940", 1.0),
            ("source-test-conditions-known", 0.0),
            ("source-cure-history-known", 0.0),
            ("source-moisture-condition-known", 0.0),
            ("source-expansion-interval-known", 0.0),
            ("btu-it-conversion", btu_it),
            ("source-btu-convention-known", btu_known),
        ])
    };
    let point = adhesive(1.0, 0.0);
    // Btu_IT, inch, square foot, hour and Fahrenheit differences are explicit.
    let conductivity = 8.0 * 1_055.055_852_62 * 0.0254 * 1.8 / (3_600.0 * 0.092_903_04);
    for (name, dims, expected) in [
        ("thermal-conductivity", K, conductivity),
        ("thermal-expansion-coefficient", EXPANSION, 4.5e-6 * 1.8),
        ("compressive-strength", PA, 4_000.0 * 6_894.757_293_168),
        ("flexural-strength", PA, 1_800.0 * 6_894.757_293_168),
        (
            "dielectric-breakdown-strength",
            ELECTRIC_FIELD,
            125.0 / 0.000_025_4,
        ),
    ] {
        sample(&ceramic, name, dims, &point, expected);
    }
    refuses(&ceramic, "thermal-conductivity", K, &adhesive(0.0, 0.0));
    refuses(&ceramic, "thermal-conductivity", K, &adhesive(1.0, 1.0));

    let laminate = compile("norplex-np509-melamine-glass-coupons-2022", 11);
    let coupon = |thickness| {
        dim(
            flags(&[
                ("source-grade-np509", 1.0),
                ("source-test-temperature-known", 0.0),
                ("np509-condition-a", 1.0),
                ("source-water-reference-temperature-known", 0.0),
                ("source-expansion-interval-known", 0.0),
                ("source-expansion-condition-known", 0.0),
            ]),
            "nominal-thickness",
            LENGTH,
            thickness,
        )
    };
    for (name, dims, thickness, expected) in [
        ("flexural-modulus-longitudinal", PA, 0.001_587_5, 13.8e9),
        ("flexural-modulus-transverse", PA, 0.001_587_5, 11.7e9),
        ("flexural-strength-longitudinal", PA, 0.001_587_5, 486.1e6),
        ("flexural-strength-transverse", PA, 0.001_587_5, 435.7e6),
        ("shear-strength-perpendicular", PA, 0.001_587_5, 124.1e6),
        ("tensile-strength-longitudinal", PA, 0.003_175, 344.7e6),
        ("tensile-strength-transverse", PA, 0.003_175, 241.3e6),
        ("compressive-strength-flatwise", PA, 0.0127, 586.1e6),
        ("specific-gravity", Dims::NONE, 0.0127, 1.9),
        (
            "thermal-expansion-coefficient-source-x",
            EXPANSION,
            0.003_175,
            15e-6,
        ),
        (
            "thermal-expansion-coefficient-source-y",
            EXPANSION,
            0.003_175,
            18e-6,
        ),
    ] {
        sample(&laminate, name, dims, &coupon(thickness), expected);
    }
    refuses(
        &laminate,
        "flexural-modulus-longitudinal",
        PA,
        &coupon(0.003_175),
    );
    refuses(
        &laminate,
        "tensile-strength-longitudinal",
        PA,
        &coupon(0.001_587_5),
    );
    // Neither isotropic elasticity nor an X/Y-to-LW/CW map is supplied.
    refuses(&laminate, "young-modulus", PA, &coupon(0.001_587_5));
    refuses(
        &laminate,
        "thermal-expansion-coefficient-longitudinal",
        EXPANSION,
        &coupon(0.003_175),
    );
}
