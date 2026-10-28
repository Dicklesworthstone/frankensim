use super::*;

const INVERSE_K: Dims = Dims([0, 0, 0, -1, 0, 0]);
const RESISTIVITY: Dims = Dims([3, 1, -3, 0, -2, 0]);

fn with_flags(point: QueryPoint, values: &[(&str, f64)]) -> QueryPoint {
    values.iter().fold(point, |p, (name, value)| {
        p.with_quantity(*name, QuantitySpec::dimensional(Dims::NONE), *value)
            .unwrap()
    })
}

fn zirconium_point(kelvin: f64) -> QueryPoint {
    temperature(
        flags(&[
            ("engineering-reference-model", 1.0),
            ("source-specimens-matched", 0.0),
        ]),
        "temperature",
        kelvin,
    )
}

fn ceramic_point(slug: &str, kelvin: f64) -> QueryPoint {
    let point = temperature(
        flags(&[
            ("engineering-reference-model", 1.0),
            ("source-temperatures-matched", 0.0),
        ]),
        "temperature",
        kelvin,
    );
    if slug == "silicon-carbide-cvd-morgan-2013" {
        with_flags(point, &[("extremely-low-resistivity-grade", 0.0)])
    } else {
        point
    }
}

#[test]
fn g0_g3_zirconium_temperature_curves_keep_processing_and_direction_average() {
    for (slug, rho, cp, k, room_e, warm_e, ultimate, proof) in [
        (
            "zirconium-702-ati-2012",
            6510.0,
            284.7,
            22.0,
            99.2e9,
            90.0e9,
            333.85e6,
            231.65e6,
        ),
        (
            "zirconium-705-ati-2012",
            6640.0,
            280.5,
            17.1,
            94.7e9,
            88.9e9,
            441.8e6,
            331.5e6,
        ),
    ] {
        let pack = compile(slug, 9);
        sample(&pack, "density", RHO, &zirconium_point(293.15), rho);
        sample(
            &pack,
            "specific-heat-capacity",
            CP,
            &zirconium_point(323.15),
            cp,
        );
        sample(&pack, "thermal-conductivity", K, &zirconium_point(500.0), k);
        refuses(&pack, "density", RHO, &zirconium_point(500.0));
        refuses(
            &pack,
            "specific-heat-capacity",
            CP,
            &zirconium_point(373.16),
        );
        refuses(&pack, "thermal-conductivity", K, &zirconium_point(800.01));

        let elastic = |kelvin, known| {
            with_flags(
                zirconium_point(kelvin),
                &[
                    ("source-elastic-temperature-known", known),
                    ("source-elastic-direction-known", 0.0),
                    ("source-elastic-temper-known", 0.0),
                ],
            )
        };
        sample(&pack, "young-modulus", PA, &elastic(293.15, 0.0), room_e);
        // Independent arithmetic between the printed 93 and 149 C rows.
        sample(&pack, "young-modulus", PA, &elastic(394.15, 1.0), warm_e);
        refuses(&pack, "young-modulus", PA, &elastic(293.15, 1.0));
        refuses(&pack, "young-modulus", PA, &elastic(394.15, 0.0));
        refuses(&pack, "young-modulus", PA, &elastic(644.16, 1.0));
        refuses(
            &pack,
            "young-modulus",
            PA,
            &elastic(394.15, 1.0)
                .with("source-elastic-direction-known", 1.0)
                .unwrap(),
        );

        let point = with_flags(
            zirconium_point(394.15),
            &[
                ("cold-worked-and-annealed", 1.0),
                ("longitudinal-transverse-average", 1.0),
                ("source-test-rate-known", 0.0),
            ],
        );
        sample(&pack, "ultimate-tensile-strength", PA, &point, ultimate);
        sample(&pack, "proof-stress-0.2-percent", PA, &point, proof);
        for wrong in [
            point.clone().with("cold-worked-and-annealed", 0.0).unwrap(),
            point
                .clone()
                .with("longitudinal-transverse-average", 0.0)
                .unwrap(),
            point.clone().with("source-test-rate-known", 1.0).unwrap(),
            point.clone().with("temperature", 644.16).unwrap(),
        ] {
            refuses(&pack, "ultimate-tensile-strength", PA, &wrong);
        }
        // Neither ASTM minima nor a directional yield law was manufactured.
        refuses(&pack, "yield-strength", PA, &point);
    }
}

#[test]
fn g0_g3_engineering_ceramics_keep_frozen_thermal_intervals_and_elastic_assumptions() {
    for (slug, count, start, rho, cp, k, elastic) in [
        (
            "macor-corning-2012",
            11,
            298.15,
            2520.0,
            790.0,
            1.46,
            Some((66.9e9, 0.29)),
        ),
        (
            "silicon-nitride-sn241o-kyocera-reference",
            7,
            293.15,
            3200.0,
            660.0,
            54.0,
            Some((290e9, 0.28)),
        ),
        (
            "silicon-nitride-sn240o-kyocera-reference",
            7,
            293.15,
            3300.0,
            650.0,
            27.0,
            Some((300e9, 0.28)),
        ),
        (
            "silicon-carbide-cvd-morgan-2013",
            6,
            293.15,
            3210.0,
            660.0,
            250.0,
            None,
        ),
    ] {
        let pack = compile(slug, count);
        let point = ceramic_point(slug, start + 0.5);
        sample(&pack, "density", RHO, &point, rho);
        sample(&pack, "specific-heat-capacity", CP, &point, cp);
        sample(&pack, "thermal-conductivity", K, &point, k);
        for wrong in [
            point
                .clone()
                .with("engineering-reference-model", 0.0)
                .unwrap(),
            point
                .clone()
                .with("source-temperatures-matched", 1.0)
                .unwrap(),
            ceramic_point(slug, start + 1.01),
        ] {
            refuses(&pack, "thermal-conductivity", K, &wrong);
        }
        let original = ceramic_point(slug, start);
        refuses(&pack, "young-modulus", PA, &original);
        if let Some((young, nu)) = elastic {
            let isotropic = with_flags(original, &[("isotropic-engineering-model", 1.0)]);
            sample(&pack, "young-modulus", PA, &isotropic, young);
            sample(&pack, "poisson-ratio", Dims::NONE, &isotropic, nu);
            refuses(
                &pack,
                "young-modulus",
                PA,
                &isotropic.with("temperature", start + 0.1).unwrap(),
            );
        }
    }
}

#[test]
fn g0_g3_ceramic_interval_expansion_keeps_both_endpoints() {
    let mean = "mean-linear-thermal-expansion-coefficient";
    for (slug, count, reference, end, expected, unreported) in [
        ("macor-corning-2012", 11, 298.15, 573.15, 9.0e-6, 673.15),
        (
            "silicon-nitride-sn241o-kyocera-reference",
            7,
            313.15,
            673.15,
            2.9e-6,
            873.15,
        ),
        (
            "silicon-nitride-sn240o-kyocera-reference",
            7,
            313.15,
            673.15,
            2.8e-6,
            873.15,
        ),
        (
            "silicon-carbide-cvd-morgan-2013",
            6,
            233.15,
            1223.15,
            4.5e-6,
            1100.0,
        ),
    ] {
        let pack = compile(slug, count);
        let base = ceramic_point(slug, end);
        refuses(&pack, mean, INVERSE_K, &base);
        let point = temperature(base.clone(), "reference-temperature", reference);
        sample(&pack, mean, INVERSE_K, &point, expected);
        refuses(
            &pack,
            mean,
            INVERSE_K,
            &point.clone().with("temperature", unreported).unwrap(),
        );
        refuses(
            &pack,
            mean,
            INVERSE_K,
            &temperature(base, "reference-temperature", reference + 1.0),
        );
        refuses(
            &pack,
            "linear-thermal-expansion-coefficient",
            INVERSE_K,
            &point,
        );
    }
    let macor = compile("macor-corning-2012", 11);
    let cold = temperature(
        ceramic_point("macor-corning-2012", 298.15),
        "reference-temperature",
        173.15,
    );
    sample(&macor, mean, INVERSE_K, &cold, 8.1e-6);
}

#[test]
fn g0_g3_ceramic_strength_modes_and_censored_resistivity_remain_distinct() {
    let macor = compile("macor-corning-2012", 11);
    let point = ceramic_point("macor-corning-2012", 298.15);
    sample(
        &macor,
        "flexural-strength-minimum-specified-average",
        PA,
        &point,
        94e6,
    );
    refuses(&macor, "flexural-strength", PA, &point);
    let compression = with_flags(point, &[("polished-strength-enhancement", 0.0)]);
    sample(&macor, "compressive-strength", PA, &compression, 345e6);
    refuses(
        &macor,
        "compressive-strength",
        PA,
        &compression
            .with("polished-strength-enhancement", 1.0)
            .unwrap(),
    );

    for (slug, strength) in [
        ("silicon-nitride-sn241o-kyocera-reference", 790e6),
        ("silicon-nitride-sn240o-kyocera-reference", 1020e6),
    ] {
        let pack = compile(slug, 7);
        let point = ceramic_point(slug, 293.15);
        sample(&pack, "flexural-strength-three-point", PA, &point, strength);
        refuses(&pack, "ultimate-tensile-strength", PA, &point);
    }

    let slug = "silicon-carbide-cvd-morgan-2013";
    let sic = compile(slug, 6);
    for (t, strength) in [(293.15, 370e6), (1573.15, 560e6)] {
        sample(
            &sic,
            "flexural-strength",
            PA,
            &ceramic_point(slug, t),
            strength,
        );
    }
    refuses(&sic, "flexural-strength", PA, &ceramic_point(slug, 800.0));
    let point = ceramic_point(slug, 293.15);
    sample(
        &sic,
        "electrical-resistivity-exclusive-lower-limit",
        RESISTIVITY,
        &point,
        0.01,
    );
    refuses(&sic, "electrical-resistivity", RESISTIVITY, &point);
    refuses(
        &sic,
        "thermal-conductivity",
        K,
        &point.with("extremely-low-resistivity-grade", 1.0).unwrap(),
    );
}
