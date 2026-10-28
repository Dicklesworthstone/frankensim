use super::*;

const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);
const TEAR: Dims = Dims([0, 1, -2, 0, 0, 0]);

fn unknown_foam_conditions(grade: &str) -> QueryPoint {
    flags(&[
        (grade, 1.0),
        ("source-temperature-known", 0.0),
        ("source-moisture-condition-known", 0.0),
        ("source-pressure-known", 0.0),
        ("source-orientation-known", 0.0),
        ("source-test-rate-known", 0.0),
        ("source-specimen-geometry-known", 0.0),
    ])
}

#[test]
fn g0_g3_cushioning_foam_compression_and_recovery_keep_their_history() {
    for (grade, count, rho, stress25, stress50, tensile, elongation, tear, hold, set30, set24, k) in [
        (
            "plastazote-ld33",
            8,
            33.0,
            66000.0,
            134000.0,
            416000.0,
            1.48,
            1888.0,
            0.25,
            0.09,
            0.04,
            Some(0.039),
        ),
        (
            "plastazote-ld45",
            8,
            45.0,
            91000.0,
            164000.0,
            537000.0,
            1.55,
            2418.0,
            0.25,
            0.07,
            0.03,
            Some(0.041),
        ),
        (
            "evazote-ev30",
            7,
            30.0,
            39000.0,
            100000.0,
            520000.0,
            2.05,
            1675.0,
            0.5,
            0.13,
            0.04,
            None,
        ),
        (
            "evazote-ev50",
            7,
            50.0,
            50000.0,
            116000.0,
            749000.0,
            2.43,
            2880.0,
            0.5,
            0.09,
            0.02,
            None,
        ),
    ] {
        let pack = compile(&format!("zotefoams-{grade}-2026"), count);
        let selector = format!("source-grade-{grade}-2026");
        let base = unknown_foam_conditions(&selector);
        let density = base
            .clone()
            .with_quantity(
                "density-includes-skins",
                QuantitySpec::dimensional(Dims::NONE),
                1.0,
            )
            .unwrap();
        sample(&pack, "density", RHO, &density, rho);
        refuses(&pack, "density", RHO, &base);
        refuses(
            &pack,
            "density",
            RHO,
            &density.clone().with("density-includes-skins", 0.0).unwrap(),
        );
        for (strain, expected) in [(0.25, stress25), (0.5, stress50)] {
            let p = base
                .clone()
                .with_quantity(
                    "compressive-strain",
                    QuantitySpec::dimensional(Dims::NONE),
                    strain,
                )
                .unwrap();
            let answer = sample(&pack, "compressive-stress", PA, &p, expected);
            assert_eq!(
                answer.evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
        }
        for strain in [0.0, 0.249, 0.375, 0.501] {
            let p = base
                .clone()
                .with_quantity(
                    "compressive-strain",
                    QuantitySpec::dimensional(Dims::NONE),
                    strain,
                )
                .unwrap();
            refuses(&pack, "compressive-stress", PA, &p);
        }
        let p = base
            .clone()
            .with_quantity(
                "compressive-strain",
                QuantitySpec::dimensional(Dims::NONE),
                0.25,
            )
            .unwrap();
        for bad in [
            p.clone().with(&selector, 0.0).unwrap(),
            p.clone().with("source-test-rate-known", 1.0).unwrap(),
            p.clone()
                .with("source-specimen-geometry-known", 1.0)
                .unwrap(),
            p.clone().with("source-temperature-known", 1.0).unwrap(),
        ] {
            refuses(&pack, "compressive-stress", PA, &bad);
        }
        sample(
            &pack,
            "cell-diameter",
            Dims([1, 0, 0, 0, 0, 0]),
            &base,
            0.0004,
        );
        sample(&pack, "tensile-strength", PA, &base, tensile);
        sample(&pack, "tensile-elongation", Dims::NONE, &base, elongation);
        sample(&pack, "tear-strength", TEAR, &base, tear);
        refuses(&pack, "compression-modulus", PA, &base);

        let history = temperature(base.clone(), "compression-temperature", 296.15)
            .with_quantity(
                "compression-duration",
                QuantitySpec::dimensional(TIME),
                79200.0,
            )
            .unwrap()
            .with_quantity(
                "compressive-strain",
                QuantitySpec::dimensional(Dims::NONE),
                hold,
            )
            .unwrap()
            .with_quantity(
                "source-recovery-temperature-known",
                QuantitySpec::dimensional(Dims::NONE),
                0.0,
            )
            .unwrap();
        for (seconds, set) in [(1800.0, set30), (86400.0, set24)] {
            let p = history
                .clone()
                .with_quantity("recovery-time", QuantitySpec::dimensional(TIME), seconds)
                .unwrap();
            sample(&pack, "compression-set", Dims::NONE, &p, set);
        }
        let p = history
            .with_quantity("recovery-time", QuantitySpec::dimensional(TIME), 86400.0)
            .unwrap();
        for bad in [
            p.clone().with("recovery-time", 3600.0).unwrap(),
            p.clone().with("compression-duration", 3600.0).unwrap(),
            p.clone().with("compressive-strain", 0.75 - hold).unwrap(),
            p.clone().with("compression-temperature", 293.15).unwrap(),
            p.clone()
                .with("source-recovery-temperature-known", 1.0)
                .unwrap(),
        ] {
            refuses(&pack, "compression-set", Dims::NONE, &bad);
        }
        let at = temperature(base.clone(), "mean-temperature", 283.15);
        if let Some(conductivity) = k {
            sample(&pack, "thermal-conductivity", K, &at, conductivity);
            refuses(
                &pack,
                "thermal-conductivity",
                K,
                &at.clone().with("mean-temperature", 293.15).unwrap(),
            );
        } else {
            refuses(&pack, "thermal-conductivity", K, &at);
        }
        refuses(&pack, "specific-heat-capacity", CP, &at);
    }
}

#[test]
fn g0_g3_poron_preserves_density_grade_speed_and_reported_limits() {
    let pack = compile("rogers-poron-4701-30-320-2023", 11);
    let base = unknown_foam_conditions("source-grade-poron-4701-30-320-2023")
        .with_quantity(
            "source-thickness-band-1p57-to-3p18-mm",
            QuantitySpec::dimensional(Dims::NONE),
            1.0,
        )
        .unwrap()
        .with_quantity(
            "source-is-typical-value",
            QuantitySpec::dimensional(Dims::NONE),
            1.0,
        )
        .unwrap();
    sample(&pack, "density", RHO, &base, 320.0);
    sample(&pack, "thermal-conductivity", K, &base, 0.076);
    sample(&pack, "tensile-strength", PA, &base, 346000.0);
    sample(&pack, "tensile-elongation", Dims::NONE, &base, 1.55);
    sample(&pack, "tear-strength", TEAR, &base, 1200.0);
    let limits = base.clone().with("source-is-typical-value", 0.0).unwrap();
    sample(&pack, "tensile-strength-minimum", PA, &limits, 207000.0);
    sample(
        &pack,
        "tensile-elongation-minimum",
        Dims::NONE,
        &limits,
        1.0,
    );
    sample(&pack, "tear-strength-minimum", TEAR, &limits, 500.0);
    refuses(&pack, "tensile-strength", PA, &limits);
    refuses(&pack, "tensile-strength-minimum", PA, &base);

    let p = base
        .clone()
        .with_quantity(
            "compressive-strain",
            QuantitySpec::dimensional(Dims::NONE),
            0.25,
        )
        .unwrap()
        .with_quantity("test-speed", QuantitySpec::dimensional(SPEED), 0.000085)
        .unwrap();
    sample(&pack, "compressive-stress", PA, &p, 35000.0);
    let limits = p.clone().with("source-is-typical-value", 0.0).unwrap();
    sample(
        &pack,
        "compressive-stress-reported-lower-limit",
        PA,
        &limits,
        21000.0,
    );
    sample(
        &pack,
        "compressive-stress-reported-upper-limit",
        PA,
        &limits,
        55000.0,
    );
    for bad in [
        p.clone().with("compressive-strain", 0.5).unwrap(),
        p.clone().with("test-speed", 0.0001).unwrap(),
        p.clone()
            .with("source-thickness-band-1p57-to-3p18-mm", 0.0)
            .unwrap(),
        p.clone()
            .with("source-grade-poron-4701-30-320-2023", 0.0)
            .unwrap(),
    ] {
        refuses(&pack, "compressive-stress", PA, &bad);
    }
    refuses(&pack, "compressive-stress", PA, &base);
    refuses(&pack, "compression-modulus", PA, &p);
    refuses(&pack, "specific-heat-capacity", CP, &base);
}

#[test]
fn g0_g3_bisco_silicone_keeps_specifications_and_compression_history_distinct() {
    let pack = compile("rogers-bisco-ht800-2026", 10);
    let base = unknown_foam_conditions("source-grade-bisco-ht800-2026")
        .with_quantity(
            "source-deflection-known",
            QuantitySpec::dimensional(Dims::NONE),
            0.0,
        )
        .unwrap()
        .with_quantity(
            "source-is-typical-value",
            QuantitySpec::dimensional(Dims::NONE),
            1.0,
        )
        .unwrap();
    sample(&pack, "density", RHO, &base, 352.0);
    sample(&pack, "compression-force-deflection", PA, &base, 67000.0);
    sample(&pack, "thermal-conductivity", K, &base, 0.076);
    let spec = base.clone().with("source-is-typical-value", 0.0).unwrap();
    sample(&pack, "density-specification-minimum", RHO, &spec, 300.0);
    sample(&pack, "density-specification-maximum", RHO, &spec, 473.0);
    sample(
        &pack,
        "compression-force-deflection-specification-minimum",
        PA,
        &spec,
        41000.0,
    );
    sample(
        &pack,
        "compression-force-deflection-specification-maximum",
        PA,
        &spec,
        97000.0,
    );
    refuses(&pack, "compression-force-deflection", PA, &spec);
    refuses(
        &pack,
        "compression-force-deflection",
        PA,
        &base.clone().with("source-deflection-known", 1.0).unwrap(),
    );
    refuses(&pack, "compressive-stress", PA, &base);
    let history = temperature(base.clone(), "compression-temperature", 373.15)
        .with_quantity(
            "compression-duration",
            QuantitySpec::dimensional(TIME),
            79200.0,
        )
        .unwrap()
        .with_quantity(
            "compressive-strain",
            QuantitySpec::dimensional(Dims::NONE),
            0.5,
        )
        .unwrap()
        .with_quantity(
            "source-recovery-time-known",
            QuantitySpec::dimensional(Dims::NONE),
            0.0,
        )
        .unwrap()
        .with_quantity(
            "source-recovery-temperature-known",
            QuantitySpec::dimensional(Dims::NONE),
            0.0,
        )
        .unwrap();
    sample(&pack, "compression-set", Dims::NONE, &history, 0.024);
    sample(
        &pack,
        "compression-set-exclusive-upper-limit",
        Dims::NONE,
        &history
            .clone()
            .with("source-is-typical-value", 0.0)
            .unwrap(),
        0.05,
    );
    for bad in [
        history
            .clone()
            .with("compression-temperature", 296.15)
            .unwrap(),
        history.clone().with("compressive-strain", 0.25).unwrap(),
        history
            .clone()
            .with("source-recovery-time-known", 1.0)
            .unwrap(),
    ] {
        refuses(&pack, "compression-set", Dims::NONE, &bad);
    }
    let electrical = base
        .clone()
        .with_quantity(
            "source-electrical-test-voltage-known",
            QuantitySpec::dimensional(Dims::NONE),
            0.0,
        )
        .unwrap();
    sample(
        &pack,
        "electrical-resistivity",
        Dims([3, 1, -3, 0, -2, 0]),
        &electrical,
        1e12,
    );
    refuses(&pack, "specific-heat-capacity", CP, &base);
}
