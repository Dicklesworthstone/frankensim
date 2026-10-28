//! Construction-specific PCB electrical facts and cured adhesive references.
//! Source transport and bounded model checks do not certify an assembly.

use super::*;

const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);
const FREQUENCY: Dims = Dims([0, 0, -1, 0, 0, 0]);
const EXPANSION: Dims = Dims([0, 0, 0, -1, 0, 0]);
const RESISTIVITY: Dims = Dims([3, 1, -3, 0, -2, 0]);

fn dim(point: QueryPoint, name: &str, dims: Dims, value: f64) -> QueryPoint {
    point
        .with_quantity(name, QuantitySpec::dimensional(dims), value)
        .unwrap()
}

fn frequency(point: QueryPoint, hz: f64) -> QueryPoint {
    dim(point, "frequency", FREQUENCY, hz)
}

#[test]
fn fr408hr_core_queries_preserve_construction_and_discrete_frequency_support() {
    let pack = compile("isola-fr408hr-core-dielectric-2026", 8);
    for (style, plies, resin, thickness, dk, df) in [
        (
            106.0,
            1.0,
            0.72,
            0.000_051,
            [3.37, 3.36, 3.34, 3.32, 3.30, 3.30],
            [0.0075, 0.0089, 0.0096, 0.0101, 0.0107, 0.0107],
        ),
        (
            2116.0,
            1.0,
            0.50,
            0.000_114,
            [3.82, 3.81, 3.80, 3.79, 3.78, 3.76],
            [0.0068, 0.0075, 0.0083, 0.0086, 0.0091, 0.0089],
        ),
        (
            2116.0,
            1.0,
            0.54,
            0.000_127,
            [3.72, 3.71, 3.70, 3.69, 3.67, 3.67],
            [0.0070, 0.0077, 0.0086, 0.0089, 0.0095, 0.0093],
        ),
        (
            2116.0,
            2.0,
            0.54,
            0.000_254,
            [3.72, 3.71, 3.70, 3.69, 3.67, 3.67],
            [0.0070, 0.0077, 0.0086, 0.0089, 0.0095, 0.0093],
        ),
    ] {
        let construction = dim(
            flags(&[
                ("source-grade-fr408hr", 1.0),
                ("source-core", 1.0),
                ("source-glass-style", style),
                ("source-core-ply-count", plies),
                ("source-resin-content", resin),
                ("source-test-temperature-known", 0.0),
                ("source-test-direction-known", 0.0),
                ("source-copper-condition-known", 0.0),
                ("source-moisture-condition-known", 0.0),
                ("source-cure-history-known", 0.0),
            ]),
            "source-core-thickness",
            LENGTH,
            thickness,
        );
        for (index, hz) in [1e8, 5e8, 1e9, 2e9, 5e9, 1e10].into_iter().enumerate() {
            let at = frequency(construction.clone(), hz);
            sample(&pack, "relative-permittivity", Dims::NONE, &at, dk[index]);
            sample(&pack, "dielectric-loss-tangent", Dims::NONE, &at, df[index]);
        }
        let at = frequency(construction.clone(), 1e9);
        for wrong in [
            frequency(construction.clone(), 0.0),
            frequency(construction.clone(), 3e9),
            frequency(construction.clone(), 1.1e10),
            dim(at.clone(), "source-core", Dims::NONE, 0.0),
            dim(at.clone(), "source-core-ply-count", Dims::NONE, 3.0),
            dim(at.clone(), "source-resin-content", Dims::NONE, 0.60),
            dim(
                at.clone(),
                "source-core-thickness",
                LENGTH,
                thickness * 1.01,
            ),
            dim(at.clone(), "source-copper-condition-known", Dims::NONE, 1.0),
            construction,
            QueryPoint::new(),
        ] {
            refuses(&pack, "relative-permittivity", Dims::NONE, &wrong);
            refuses(&pack, "dielectric-loss-tangent", Dims::NONE, &wrong);
        }
        // Electrical construction data do not create a generic PCB thermal law.
        refuses(&pack, "thermal-conductivity", K, &at);
        refuses(&pack, "density", RHO, &at);
    }
}

fn gray_2216() -> QueryPoint {
    flags(&[
        ("source-grade-2216-gray", 1.0),
        ("base-to-accelerator-volume-ratio", 2.0 / 3.0),
        ("source-test-temperature-known", 0.0),
        ("source-cure-history-known", 0.0),
        ("source-moisture-condition-known", 0.0),
        ("source-test-direction-known", 0.0),
    ])
}

#[test]
fn cured_2216_thermal_and_dielectric_references_preserve_their_own_conditions() {
    let pack = compile("3m-2216-gray-cured-2024", 6);
    let bulk = gray_2216();
    sample(&pack, "thermal-conductivity", K, &bulk, 0.3943);
    let electrical = temperature(frequency(bulk.clone(), 1000.0), "temperature", 296.15);
    sample(
        &pack,
        "relative-permittivity",
        Dims::NONE,
        &electrical,
        5.51,
    );
    sample(
        &pack,
        "dielectric-loss-tangent",
        Dims::NONE,
        &electrical,
        0.112,
    );
    sample(
        &pack,
        "electrical-resistivity",
        RESISTIVITY,
        &electrical,
        1.9e10,
    );
    for wrong in [
        frequency(electrical.clone(), 2000.0),
        temperature(electrical.clone(), "temperature", 297.15),
        dim(
            electrical.clone(),
            "base-to-accelerator-volume-ratio",
            Dims::NONE,
            1.0,
        ),
        bulk.clone(),
    ] {
        refuses(&pack, "relative-permittivity", Dims::NONE, &wrong);
    }
    let interval = temperature(
        temperature(bulk.clone(), "source-reference-temperature", 313.15),
        "source-range-end-temperature",
        353.15,
    );
    sample(
        &pack,
        "mean-linear-thermal-expansion-coefficient",
        EXPANSION,
        &interval,
        134e-6,
    );
    refuses(
        &pack,
        "linear-thermal-expansion-coefficient",
        EXPANSION,
        &interval,
    );
    refuses(
        &pack,
        "mean-linear-thermal-expansion-coefficient",
        EXPANSION,
        &temperature(interval, "source-range-end-temperature", 363.15),
    );
    for property in ["density", "specific-heat-capacity"] {
        refuses(
            &pack,
            property,
            if property == "density" { RHO } else { CP },
            &bulk,
        );
    }
    let normalized = fs_matdb::NormalizedMaterialCardPack::new(
        fs_matdb::MaterialStateId {
            chemistry: pack.pack_id().into(),
            phase: "cured epoxy".into(),
            process: "2216 Gray source reference; thermal test and cure conditions unstated".into(),
            revision: 0,
        },
        pack,
    )
    .unwrap();
    // Actual conduction owner, with a frozen source conductivity and caller
    // geometry. Unknown source temperature is not promoted into a k(T) law.
    let slab = fs_conduction::ThermalResistanceTerm::slab_from_card(
        "1 mm cured adhesive reference layer",
        0.001,
        0.01,
        normalized.card(),
        &key("thermal-conductivity", K),
        &bulk,
        SelectionPolicy::SingleClaimOnly,
    )
    .unwrap();
    assert!((slab.value_k_per_w() - 0.001 / (0.3943 * 0.01)).abs() < 1e-14);
    assert!(
        fs_conduction::ThermalResistanceTerm::slab_from_card(
            "different formulation",
            0.001,
            0.01,
            normalized.card(),
            &key("thermal-conductivity", K),
            &dim(bulk, "source-grade-2216-gray", Dims::NONE, 0.0),
            SelectionPolicy::SingleClaimOnly,
        )
        .is_err()
    );
}

#[test]
fn adhesive_joint_modulus_never_becomes_bulk_elasticity_or_a_cure_temperature_curve() {
    let pack = compile("3m-2216-gray-cured-2024", 6);
    let joint = temperature(
        dim(
            dim(
                dim(gray_2216(), "source-etched-aluminum-joint", Dims::NONE, 1.0),
                "source-bondline-thickness-known",
                Dims::NONE,
                0.0,
            ),
            "source-dwell-duration",
            TIME,
            7200.0,
        ),
        "source-dwell-temperature",
        339.15,
    )
    .with_quantity(
        "source-environmental-pressure",
        QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
        14000.0,
    )
    .unwrap();
    for (kelvin, mpa) in [
        (173.15, 2745.0),
        (213.15, 2199.0),
        (233.15, 1947.0),
        (273.15, 1500.0),
        (296.15, 342.0),
    ] {
        let at = temperature(joint.clone(), "temperature", kelvin);
        sample(&pack, "bondline-shear-modulus", PA, &at, mpa * 1e6);
        refuses(&pack, "shear-modulus", PA, &at);
        refuses(&pack, "young-modulus", PA, &at);
    }
    let at = temperature(joint.clone(), "temperature", 296.15);
    for wrong in [
        temperature(joint.clone(), "temperature", 250.0),
        temperature(joint.clone(), "temperature", 339.15),
        dim(at.clone(), "source-etched-aluminum-joint", Dims::NONE, 0.0),
        dim(at.clone(), "source-dwell-duration", TIME, 3600.0),
        at.clone()
            .with_quantity(
                "source-environmental-pressure",
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::Pressure,
                    ValueForm::Static,
                )),
                0.0,
            )
            .unwrap(),
        temperature(gray_2216(), "temperature", 296.15),
        joint,
    ] {
        refuses(&pack, "bondline-shear-modulus", PA, &wrong);
    }
}

fn gfrp_coupon(family: u8, wet: bool) -> QueryPoint {
    let (plies, fiber, minimum, maximum) = match family {
        1 => (12.0, 0.4742, 0.000_259_08, 0.000_271_78),
        2 => (12.0, 0.4628, 0.000_254, 0.000_266_7),
        3 => (8.0, 0.4615, 0.000_254, 0.000_271_78),
        _ => panic!("only the three sourced coupon families exist"),
    };
    let p = temperature(
        dim(
            dim(
                flags(&[
                    ("source-mtm45-1-6781-35rc", 1.0),
                    ("source-cure-cycle-mh", 1.0),
                    ("source-measured-not-normalized", 1.0),
                    ("source-coupon-family", f64::from(family)),
                    ("source-ply-count", plies),
                    ("source-fiber-volume", fiber),
                    ("source-dry-condition", if wet { 0.0 } else { 1.0 }),
                    ("source-wet-equilibrium", if wet { 1.0 } else { 0.0 }),
                    ("source-density-test-temperature-known", 0.0),
                ]),
                "source-cured-ply-thickness-min",
                LENGTH,
                minimum,
            ),
            "source-cured-ply-thickness-max",
            LENGTH,
            maximum,
        ),
        "source-conditioning-temperature",
        344.261_111_111_111_06,
    );
    if wet {
        dim(p, "source-conditioning-relative-humidity", Dims::NONE, 0.85)
    } else {
        p
    }
}

#[test]
fn s2_glass_epoxy_preserves_direction_layup_measured_basis_and_dry_wet_states() {
    let pack = compile("gfrp-mtm45-1-6781-s2-ncamp-2013", 17);
    let psi = 6894.757_293_168;
    for (family, property, scale, values) in [
        (1, "tensile-modulus-warp", 1e6, [4.21, 4.14, 3.80, 3.68]),
        (
            1,
            "tensile-strength-warp",
            1e3,
            [90.05, 79.92, 56.11, 56.18],
        ),
        (2, "tensile-modulus-fill", 1e6, [4.06, 3.99, 3.69, 3.53]),
        (
            2,
            "tensile-strength-fill",
            1e3,
            [88.17, 78.92, 56.07, 54.10],
        ),
        (3, "shear-modulus-in-plane", 1e6, [0.71, 0.55, 0.34, 0.27]),
        (
            3,
            "shear-stress-at-5-percent-strain",
            1e3,
            [12.29, 9.16, 5.64, 4.57],
        ),
        (
            3,
            "shear-stress-at-0.2-percent-offset",
            1e3,
            [7.25, 5.45, 3.38, 2.66],
        ),
    ] {
        for (index, kelvin) in [
            219.261_111_111_111_1,
            297.038_888_888_888_9,
            366.483_333_333_333_3,
            394.261_111_111_111_06,
        ]
        .into_iter()
        .enumerate()
        {
            let wet = index >= 2;
            let at = temperature(gfrp_coupon(family, wet), "temperature", kelvin);
            sample(&pack, property, PA, &at, values[index] * psi * scale);
            // A hot-wet value cannot stand in for the unmeasured hot-dry
            // condition; nor can another cure or normalized basis be selected.
            for wrong in [
                temperature(gfrp_coupon(family, !wet), "temperature", kelvin),
                dim(at.clone(), "source-cure-cycle-mh", Dims::NONE, 0.0),
                dim(
                    at.clone(),
                    "source-measured-not-normalized",
                    Dims::NONE,
                    0.0,
                ),
                dim(at.clone(), "source-ply-count", Dims::NONE, 16.0),
                dim(at.clone(), "source-fiber-volume", Dims::NONE, 0.60),
                dim(at.clone(), "source-coupon-family", Dims::NONE, 4.0),
                dim(
                    at.clone(),
                    "source-cured-ply-thickness-min",
                    LENGTH,
                    0.000_256_54,
                ),
                temperature(at.clone(), "source-conditioning-temperature", kelvin),
                temperature(gfrp_coupon(family, wet), "temperature", kelvin + 0.01),
                gfrp_coupon(family, wet),
                QueryPoint::new(),
            ] {
                refuses(&pack, property, PA, &wrong);
            }
            if wet {
                refuses(
                    &pack,
                    property,
                    PA,
                    &dim(
                        at,
                        "source-conditioning-relative-humidity",
                        Dims::NONE,
                        0.50,
                    ),
                );
            }
        }
    }
    for family in 1..=3 {
        let p = gfrp_coupon(family, false);
        sample(&pack, "source-group-density", RHO, &p, 1790.0);
        refuses(
            &pack,
            "source-group-density",
            RHO,
            &dim(
                p.clone(),
                "source-density-test-temperature-known",
                Dims::NONE,
                1.0,
            ),
        );
        // Group physical tests have no reported temperature/moisture match to
        // individual mechanical tests, so generic density and tensor completion
        // remain unavailable instead of manufacturing a common specimen state.
        for (property, dims) in [
            ("density", RHO),
            ("young-modulus", PA),
            ("shear-modulus", PA),
            ("shear-strength", PA),
            ("thermal-conductivity", K),
        ] {
            refuses(&pack, property, dims, &p);
        }
    }
}
