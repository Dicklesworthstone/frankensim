use super::*;

fn polymer_point(grade: &str, value: f64) -> QueryPoint {
    flags(&[
        (grade, value),
        ("source-temperature-known", 0.0),
        ("source-moisture-condition-known", 0.0),
        ("source-pressure-known", 0.0),
        ("source-orientation-known", 0.0),
    ])
}

fn speed(point: QueryPoint, metres_per_second: f64) -> QueryPoint {
    point
        .with_quantity(
            "test-speed",
            QuantitySpec::dimensional(SPEED),
            metres_per_second,
        )
        .unwrap()
}

#[test]
fn engineering_polymers_bind_grade_rate_and_discrete_compression() {
    for (slug, count, density, modulus, conductivity, compression) in [
        (
            "ensinger-tecapeek-gf30-natural",
            6,
            1530.0,
            6.3e9,
            Some(0.35),
            52e6,
        ),
        (
            "ensinger-tecatron-gf40-natural",
            9,
            1630.0,
            6.5e9,
            Some(0.35),
            41e6,
        ),
        (
            "ensinger-tecatron-sx-natural",
            9,
            1360.0,
            4.0e9,
            Some(0.25),
            38e6,
        ),
        (
            "ensinger-tecapai-cm-xp403-green",
            6,
            1410.0,
            3.6e9,
            None,
            32e6,
        ),
    ] {
        let pack = compile(slug, count);
        let grade = format!("source-grade-{}", slug.strip_prefix("ensinger-").unwrap());
        let base = polymer_point(&grade, 1.0);
        sample(&pack, "density", RHO, &base, density);
        refuses(&pack, "density", RHO, &polymer_point(&grade, 0.0));
        sample(
            &pack,
            "tensile-modulus",
            PA,
            &speed(base.clone(), 1.0 / 60000.0),
            modulus,
        );
        refuses(&pack, "tensile-modulus", PA, &base);
        refuses(
            &pack,
            "tensile-modulus",
            PA,
            &speed(base.clone(), 5.0 / 60000.0),
        );
        if let Some(k) = conductivity {
            sample(&pack, "thermal-conductivity", K, &base, k);
            sample(&pack, "specific-heat-capacity", CP, &base, 1000.0);
        } else {
            refuses(&pack, "thermal-conductivity", K, &base);
            refuses(&pack, "specific-heat-capacity", CP, &base);
            sample(
                &pack,
                "tensile-strength-at-break",
                PA,
                &speed(base.clone(), 5.0 / 60000.0),
                122e6,
            );
            refuses(&pack, "tensile-strength", PA, &base);
        }
        let compression_point = if conductivity.is_some() {
            speed(base.clone(), 5.0 / 60000.0)
        } else {
            base.clone()
                .with_quantity(
                    "source-test-rate-known",
                    QuantitySpec::dimensional(Dims::NONE),
                    0.0,
                )
                .unwrap()
        };
        let at_strain = |strain| {
            compression_point
                .clone()
                .with_quantity(
                    "compressive-strain",
                    QuantitySpec::dimensional(Dims::NONE),
                    strain,
                )
                .unwrap()
        };
        sample(
            &pack,
            "compressive-stress",
            PA,
            &at_strain(0.02),
            compression,
        );
        refuses(&pack, "compressive-stress", PA, &at_strain(0.03));
    }
}

fn foam_point(grade: &str, kelvin: f64) -> QueryPoint {
    temperature(
        flags(&[
            (grade, 1.0),
            ("source-moisture-condition-known", 0.0),
            ("source-pressure-known", 0.0),
        ]),
        "temperature",
        kelvin,
    )
}

fn flag(point: QueryPoint, name: &str, value: f64) -> QueryPoint {
    point
        .with_quantity(name, QuantitySpec::dimensional(Dims::NONE), value)
        .unwrap()
}

#[test]
fn pvc_foam_grades_keep_nominal_minimum_direction_and_temperature_distinct() {
    let keys = [
        "compression-strength",
        "compression-modulus",
        "tensile-strength",
        "tensile-modulus",
        "shear-strength",
        "shear-modulus",
    ];
    for (density, k, nominal, minimum) in [
        (
            60,
            0.029,
            [0.9, 70.0, 1.8, 75.0, 0.76, 20.0],
            [0.7, 60.0, 1.5, 57.0, 0.63, 16.0],
        ),
        (
            80,
            0.031,
            [1.4, 90.0, 2.5, 95.0, 1.15, 27.0],
            [1.15, 80.0, 2.2, 85.0, 0.95, 23.0],
        ),
        (
            130,
            0.036,
            [3.0, 170.0, 4.8, 175.0, 2.2, 50.0],
            [2.4, 145.0, 3.5, 135.0, 1.9, 40.0],
        ),
        (
            200,
            0.044,
            [5.4, 310.0, 7.1, 250.0, 3.5, 73.0],
            [4.5, 265.0, 6.3, 210.0, 3.2, 65.0],
        ),
    ] {
        let pack = compile(&format!("diab-divinycell-h{density}-2021"), 14);
        let grade = format!("source-grade-divinycell-h{density}-2021");
        let base = foam_point(&grade, 296.15);
        sample(
            &pack,
            "density",
            RHO,
            &flag(base.clone(), "source-is-nominal-value", 1.0),
            f64::from(density),
        );
        let thermal = |kelvin| flag(foam_point(&grade, kelvin), "source-orientation-known", 0.0);
        sample(&pack, "thermal-conductivity", K, &thermal(293.15), k);
        refuses(&pack, "thermal-conductivity", K, &thermal(296.15));
        refuses(&pack, "specific-heat-capacity", CP, &thermal(293.15));

        for (index, property) in keys.iter().enumerate() {
            let mechanical = flag(base.clone(), "source-test-rate-known", 0.0);
            let oriented = if property.starts_with("shear-") {
                flag(mechanical.clone(), "source-orientation-known", 0.0)
            } else {
                flag(mechanical.clone(), "normal-to-sheet", 1.0)
            };
            let average = flag(oriented.clone(), "source-is-nominal-value", 1.0);
            let lower = flag(oriented, "source-is-nominal-value", 0.0);
            let minimum_property = format!("{property}-minimum");
            sample(&pack, property, PA, &average, nominal[index] * 1e6);
            sample(&pack, &minimum_property, PA, &lower, minimum[index] * 1e6);
            refuses(&pack, &minimum_property, PA, &average);
            refuses(
                &pack,
                property,
                PA,
                &flag(mechanical, "source-is-nominal-value", 1.0),
            );
        }
    }
}
