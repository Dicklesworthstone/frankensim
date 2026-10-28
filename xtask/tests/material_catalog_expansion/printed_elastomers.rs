use super::*;

const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);
const TEAR: Dims = Dims([0, 1, -2, 0, 0, 0]);

fn dim(point: QueryPoint, name: &str, dims: Dims, value: f64) -> QueryPoint {
    point
        .with_quantity(name, QuantitySpec::dimensional(dims), value)
        .unwrap()
}

fn prusa_point(grade: &str, nozzle: f64, bed: f64, orientation: f64) -> QueryPoint {
    let point = flags(&[
        (grade, 1.0),
        ("source-print-profile-mk3-slic3r140", 1.0),
        ("print-orientation", orientation),
        ("infill-fraction", 1.0),
        ("source-test-temperature-known", 0.0),
        ("source-moisture-condition-known", 0.0),
        ("source-test-rate-known", 0.0),
        ("source-postprocessing-known", 0.0),
    ]);
    let point = dim(point, "layer-height", LENGTH, 0.0002);
    let point = dim(point, "infill-print-speed", SPEED, 0.2);
    temperature(
        temperature(point, "nozzle-temperature", nozzle),
        "bed-temperature",
        bed,
    )
}

#[test]
fn printed_pla_petg_keep_build_profile_orientation_and_test_meaning() {
    let properties = [
        ("tensile-modulus", PA),
        ("tensile-yield-strength", PA),
        ("tensile-elongation-at-yield", Dims::NONE),
        ("flexural-modulus", PA),
        ("flexural-strength", PA),
    ];
    for (resin, nozzle, bed, expected) in [
        (
            "pla",
            488.15,
            333.15,
            [
                [2.3e9, 51e6, 0.029, 3.1e9, 83e6],
                [2.4e9, 59e6, 0.032, 3.2e9, 99e6],
            ],
        ),
        (
            "petg",
            523.15,
            353.15,
            [
                [1.5e9, 47e6, 0.051, 1.7e9, 66e6],
                [1.6e9, 50e6, 0.051, 1.6e9, 70e6],
            ],
        ),
    ] {
        let pack = compile(&format!("prusament-{resin}-mk3-printed-2022"), 10);
        let grade = format!("source-grade-prusament-{resin}");
        for (orientation, values) in expected.into_iter().enumerate() {
            let point = prusa_point(&grade, nozzle, bed, orientation as f64);
            for ((property, dims), value) in properties.iter().zip(values) {
                sample(&pack, property, *dims, &point, value);
            }
            // Feedstock density and flexural stiffness cannot complete a
            // three-dimensional elastic model of these printed specimens.
            for (property, dims) in [
                ("density", RHO),
                ("young-modulus", PA),
                ("compression-modulus", PA),
                ("tensile-strength", PA),
                ("poisson-ratio", Dims::NONE),
                ("thermal-conductivity", K),
            ] {
                refuses(&pack, property, dims, &point);
            }
            for bad in [
                dim(point.clone(), "layer-height", LENGTH, 0.0003),
                dim(point.clone(), "infill-fraction", Dims::NONE, 0.2),
                dim(point.clone(), "infill-print-speed", SPEED, 0.1),
                dim(point.clone(), &grade, Dims::NONE, 0.0),
                dim(point.clone(), "print-orientation", Dims::NONE, 2.0),
                dim(point.clone(), "print-orientation", Dims::NONE, 0.5),
                dim(
                    point.clone(),
                    "source-test-temperature-known",
                    Dims::NONE,
                    1.0,
                ),
                dim(point.clone(), "source-test-rate-known", Dims::NONE, 1.0),
                temperature(point.clone(), "nozzle-temperature", nozzle + 5.0),
                temperature(point.clone(), "bed-temperature", bed + 5.0),
            ] {
                refuses(&pack, "tensile-modulus", PA, &bad);
            }
        }
        refuses(&pack, "tensile-modulus", PA, &QueryPoint::new());
    }
}

fn tpu_point(orientation: f64) -> QueryPoint {
    let point = flags(&[
        ("source-grade-ultrafuse-tpu95a", 1.0),
        ("source-print-profile-mk3-tpu95a-2022", 1.0),
        ("print-orientation", orientation),
        ("conditioning-relative-humidity", 0.5),
        ("source-print-infill-known", 0.0),
        ("source-print-layer-height-known", 0.0),
        ("source-test-temperature-known", 0.0),
        ("source-test-rate-known", 0.0),
        ("source-filament-conditioning-known", 0.0),
    ]);
    let point = dim(point, "tempering-duration", TIME, 72_000.0);
    let point = dim(point, "conditioning-duration", TIME, 259_200.0);
    temperature(
        temperature(point, "tempering-temperature", 373.15),
        "conditioning-temperature",
        296.15,
    )
}

#[test]
fn printed_tpu_preserves_tempering_rates_and_discrete_large_strain_data() {
    let pack = compile("ultrafuse-tpu95a-mk3-printed-2022", 12);
    for (orientation, modulus, strength, elongation, stresses) in [
        (0.0, 48.4e6, 44.2e6, 6.61, [8.3e6, 10.5e6]),
        (2.0, 46.7e6, 12.2e6, 1.92, [7.9e6, 9.9e6]),
    ] {
        let point = tpu_point(orientation);
        let slow = dim(point.clone(), "test-speed", SPEED, 1.0 / 60_000.0);
        let fast = dim(point.clone(), "test-speed", SPEED, 200.0 / 60_000.0);
        sample(&pack, "tensile-modulus", PA, &slow, modulus);
        sample(&pack, "tensile-strength-at-break", PA, &fast, strength);
        sample(
            &pack,
            "tensile-elongation-at-break",
            Dims::NONE,
            &fast,
            elongation,
        );
        for (strain, stress) in [0.5, 1.0].into_iter().zip(stresses) {
            sample(
                &pack,
                "tensile-stress",
                PA,
                &dim(fast.clone(), "tensile-strain", Dims::NONE, strain),
                stress,
            );
        }
        for strain in [0.0, 0.75, 2.0, 4.0] {
            refuses(
                &pack,
                "tensile-stress",
                PA,
                &dim(fast.clone(), "tensile-strain", Dims::NONE, strain),
            );
        }
        refuses(&pack, "tensile-modulus", PA, &fast);
        refuses(&pack, "tensile-strength-at-break", PA, &slow);
        refuses(&pack, "density", RHO, &point);
        sample(&pack, "filament-density", RHO, &point, 1149.0);
        for bad in [
            temperature(slow.clone(), "tempering-temperature", 296.15),
            dim(slow.clone(), "tempering-duration", TIME, 0.0),
            dim(slow.clone(), "conditioning-duration", TIME, 0.0),
            dim(
                slow.clone(),
                "conditioning-relative-humidity",
                Dims::NONE,
                0.2,
            ),
            dim(slow.clone(), "source-print-infill-known", Dims::NONE, 1.0),
            dim(
                slow.clone(),
                "source-test-temperature-known",
                Dims::NONE,
                1.0,
            ),
            dim(
                slow.clone(),
                "source-print-profile-mk3-tpu95a-2022",
                Dims::NONE,
                0.0,
            ),
        ] {
            refuses(&pack, "tensile-modulus", PA, &bad);
        }
    }
    let strain_three = dim(
        dim(tpu_point(0.0), "test-speed", SPEED, 200.0 / 60_000.0),
        "tensile-strain",
        Dims::NONE,
        3.0,
    );
    sample(&pack, "tensile-stress", PA, &strain_three, 20.3e6);
    refuses(
        &pack,
        "tensile-stress",
        PA,
        &dim(strain_three, "print-orientation", Dims::NONE, 2.0),
    );
    for (orientation, tear) in [(0.0, 90_000.0), (1.0, 8_000.0), (2.0, 14_000.0)] {
        sample(
            &pack,
            "tear-strength-iso34-1a",
            TEAR,
            &tpu_point(orientation),
            tear,
        );
    }
    refuses(
        &pack,
        "tensile-modulus",
        PA,
        &dim(tpu_point(1.0), "test-speed", SPEED, 1.0 / 60_000.0),
    );
}

#[test]
fn named_epdm_nbr_preserve_tensile_stress_and_compression_set_history() {
    for (slug, grade, density, strength, elongation, stress, sets) in [
        (
            "parker-e3609-epdm-2023",
            "source-grade-parker-e3609",
            1080.0,
            15.1e6,
            1.79,
            6.4e6,
            [(398.15, 0.09), (423.15, 0.28)],
        ),
        (
            "parker-n0674-nbr-2024",
            "source-grade-parker-n0674",
            1220.0,
            16.2e6,
            3.31,
            5.4e6,
            [(373.15, 0.14), (398.15, 0.34)],
        ),
    ] {
        let pack = compile(slug, 6);
        let point = flags(&[
            (grade, 1.0),
            ("source-specimen-din53504-s2", 1.0),
            ("source-method-iso815-1-b-a", 1.0),
            ("source-test-temperature-known", 0.0),
            ("source-test-rate-known", 0.0),
            ("source-moisture-condition-known", 0.0),
            ("source-cure-history-known", 0.0),
            ("source-recovery-history-known", 0.0),
            ("imposed-compression-strain", 0.25),
        ]);
        let answer = sample(&pack, "density", RHO, &point, density);
        assert_eq!(
            answer.evidence.value.uncertainty,
            UncertaintyModel::Unstated
        );
        sample(&pack, "tensile-strength", PA, &point, strength);
        sample(
            &pack,
            "tensile-elongation-at-break",
            Dims::NONE,
            &point,
            elongation,
        );
        sample(
            &pack,
            "tensile-stress",
            PA,
            &dim(point.clone(), "tensile-strain", Dims::NONE, 1.0),
            stress,
        );
        for strain in [0.0, 0.5, 1.5] {
            refuses(
                &pack,
                "tensile-stress",
                PA,
                &dim(point.clone(), "tensile-strain", Dims::NONE, strain),
            );
        }
        for (dwell_temperature, set) in sets {
            let dwell = temperature(
                dim(point.clone(), "dwell-duration", TIME, 252_000.0),
                "dwell-temperature",
                dwell_temperature,
            );
            sample(&pack, "compression-set", Dims::NONE, &dwell, set);
            for bad in [
                dim(dwell.clone(), "dwell-duration", TIME, 3600.0),
                dim(dwell.clone(), "imposed-compression-strain", Dims::NONE, 0.5),
                dim(dwell.clone(), "source-method-iso815-1-b-a", Dims::NONE, 0.0),
                dim(
                    dwell.clone(),
                    "source-recovery-history-known",
                    Dims::NONE,
                    1.0,
                ),
                temperature(dwell.clone(), "dwell-temperature", dwell_temperature + 1.0),
            ] {
                refuses(&pack, "compression-set", Dims::NONE, &bad);
            }
        }
        for (name, dims) in [
            ("young-modulus", PA),
            ("tensile-modulus", PA),
            ("poisson-ratio", Dims::NONE),
            ("thermal-conductivity", K),
        ] {
            refuses(&pack, name, dims, &point);
        }
        for bad in [
            dim(point.clone(), grade, Dims::NONE, 0.0),
            dim(point.clone(), "source-cure-history-known", Dims::NONE, 1.0),
            dim(
                point.clone(),
                "source-test-temperature-known",
                Dims::NONE,
                1.0,
            ),
            dim(
                point.clone(),
                "source-specimen-din53504-s2",
                Dims::NONE,
                0.0,
            ),
        ] {
            refuses(&pack, "tensile-strength", PA, &bad);
        }
    }
}

#[test]
fn cast_rubbers_keep_seven_day_cure_and_distinct_tear_geometries() {
    const PSI: f64 = 6894.757293168;
    const PLI: f64 = 4.4482216152605 / 0.0254;
    for (slug, grade, strength_psi, stress_psi, elongation, die, tear_pli) in [
        (
            "smooth-on-ecoflex0030-seven-day-2024",
            "source-grade-ecoflex0030",
            200.0,
            10.0,
            9.0,
            "b",
            38.0,
        ),
        (
            "smooth-on-dragon-skin20-seven-day-2025",
            "source-grade-dragon-skin20",
            550.0,
            49.0,
            6.2,
            "b",
            120.0,
        ),
        (
            "smooth-on-vytaflex40-seven-day-2024",
            "source-grade-vytaflex40",
            522.0,
            100.0,
            6.6,
            "c",
            82.0,
        ),
    ] {
        let pack = compile(slug, 4);
        let point = flags(&[
            (grade, 1.0),
            ("part-a-to-b-volume-ratio", 1.0),
            ("source-moisture-condition-known", 0.0),
            ("source-test-rate-known", 0.0),
            ("source-processing-history-known", 0.0),
        ]);
        let point = temperature(
            dim(point, "cure-age", TIME, 604_800.0),
            "temperature",
            296.15,
        );
        sample(&pack, "tensile-strength", PA, &point, strength_psi * PSI);
        sample(
            &pack,
            "tensile-elongation-at-break",
            Dims::NONE,
            &point,
            elongation,
        );
        sample(
            &pack,
            &format!("tear-strength-die-{die}"),
            TEAR,
            &point,
            tear_pli * PLI,
        );
        sample(
            &pack,
            "tensile-stress",
            PA,
            &dim(point.clone(), "tensile-strain", Dims::NONE, 1.0),
            stress_psi * PSI,
        );
        for strain in [0.0, 0.5, 2.0] {
            refuses(
                &pack,
                "tensile-stress",
                PA,
                &dim(point.clone(), "tensile-strain", Dims::NONE, strain),
            );
        }
        let other_die = if die == "b" { "c" } else { "b" };
        refuses(
            &pack,
            &format!("tear-strength-die-{other_die}"),
            TEAR,
            &point,
        );
        for (name, dims) in [
            ("density", RHO),
            ("young-modulus", PA),
            ("tensile-modulus", PA),
            ("poisson-ratio", Dims::NONE),
        ] {
            refuses(&pack, name, dims, &point);
        }
        for bad in [
            dim(point.clone(), "cure-age", TIME, 14_400.0),
            dim(point.clone(), "cure-age", TIME, 57_600.0),
            dim(point.clone(), "part-a-to-b-volume-ratio", Dims::NONE, 2.0),
            dim(point.clone(), grade, Dims::NONE, 0.0),
            dim(point.clone(), "source-test-rate-known", Dims::NONE, 1.0),
            temperature(point.clone(), "temperature", 298.15),
        ] {
            refuses(&pack, "tensile-strength", PA, &bad);
        }
    }
}

fn poron_point(density: f64, nominal: f64) -> QueryPoint {
    let point = flags(&[
        ("source-grade-poron4701-40", 1.0),
        ("source-is-nominal-value", nominal),
        ("compressive-strain", 0.25),
        ("source-test-temperature-known", 0.0),
        ("source-moisture-condition-known", 0.0),
        ("source-test-rate-known", 0.0),
        ("source-specimen-thickness-known", 0.0),
        ("source-orientation-known", 0.0),
    ]);
    dim(
        dim(point, "source-nominal-density", RHO, density),
        "test-speed",
        SPEED,
        0.000085,
    )
}

#[test]
fn flexible_foam_preserves_density_grade_limits_and_compression_test_speed() {
    for (slug, density, compression, minimum, maximum, tensile, tear, tear_minimum, conductivity) in [
        (
            "rogers-poron-4701-40-15-2024",
            240.0,
            41_000.0,
            27_000.0,
            76_000.0,
            276_000.0,
            1600.0,
            500.0,
            0.065,
        ),
        (
            "rogers-poron-4701-40-20-2024",
            320.0,
            76_000.0,
            48_000.0,
            90_000.0,
            518_000.0,
            2100.0,
            900.0,
            0.080,
        ),
        (
            "rogers-poron-4701-40-30-2024",
            480.0,
            173_000.0,
            104_000.0,
            276_000.0,
            829_000.0,
            3000.0,
            2100.0,
            0.127,
        ),
    ] {
        let pack = compile(slug, 9);
        let nominal = poron_point(density, 1.0);
        let limits = poron_point(density, 0.0);
        let answer = sample(&pack, "density", RHO, &nominal, density);
        assert_eq!(
            answer.evidence.value.uncertainty,
            UncertaintyModel::Unstated
        );
        sample(&pack, "compressive-stress", PA, &nominal, compression);
        sample(&pack, "compressive-stress-minimum", PA, &limits, minimum);
        sample(&pack, "compressive-stress-maximum", PA, &limits, maximum);
        sample(&pack, "tensile-strength-minimum", PA, &limits, tensile);
        sample(
            &pack,
            "tensile-elongation-minimum",
            Dims::NONE,
            &limits,
            1.0,
        );
        sample(&pack, "tear-strength-die-c", TEAR, &nominal, tear);
        sample(
            &pack,
            "tear-strength-die-c-minimum",
            TEAR,
            &limits,
            tear_minimum,
        );
        sample(&pack, "thermal-conductivity", K, &nominal, conductivity);
        refuses(&pack, "compressive-stress", PA, &limits);
        refuses(&pack, "compressive-stress-minimum", PA, &nominal);
        refuses(&pack, "compressive-stress-maximum", PA, &nominal);
        refuses(&pack, "tensile-strength-minimum", PA, &nominal);
        for bad in [
            dim(
                nominal.clone(),
                "source-nominal-density",
                RHO,
                density + 1.0,
            ),
            dim(
                nominal.clone(),
                "source-grade-poron4701-40",
                Dims::NONE,
                0.0,
            ),
            dim(nominal.clone(), "compressive-strain", Dims::NONE, 0.5),
            dim(nominal.clone(), "test-speed", SPEED, 0.00017),
            dim(
                nominal.clone(),
                "source-test-temperature-known",
                Dims::NONE,
                1.0,
            ),
            dim(
                nominal.clone(),
                "source-specimen-thickness-known",
                Dims::NONE,
                1.0,
            ),
            dim(nominal.clone(), "source-orientation-known", Dims::NONE, 1.0),
        ] {
            refuses(&pack, "compressive-stress", PA, &bad);
        }
        for (name, dims) in [
            ("tensile-strength", PA),
            ("compression-modulus", PA),
            ("young-modulus", PA),
            ("specific-heat", CP),
        ] {
            refuses(&pack, name, dims, &nominal);
        }
        refuses(&pack, "thermal-conductivity", K, &QueryPoint::new());
        refuses(
            &pack,
            "thermal-conductivity",
            K,
            &dim(nominal, "source-test-temperature-known", Dims::NONE, 1.0),
        );
    }
}
