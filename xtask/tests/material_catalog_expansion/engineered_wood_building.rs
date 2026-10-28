//! Exact engineered-wood constructions and bounded building-material uses.
//! Declared grade statistics remain distinct from measured coupon responses.

use super::*;

trait DimensionlessSelector {
    fn with_selector(self, name: &str, value: f64) -> Result<QueryPoint, fs_matdb::MatDbError>;
}

impl DimensionlessSelector for QueryPoint {
    fn with_selector(self, name: &str, value: f64) -> Result<QueryPoint, fs_matdb::MatDbError> {
        self.with_quantity(name, QuantitySpec::dimensional(Dims::NONE), value)
    }
}

fn thickness(point: QueryPoint, metres: f64) -> QueryPoint {
    point
        .with_quantity(
            "nominal-thickness",
            QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0])),
            metres,
        )
        .unwrap()
}

fn reference_card(pack: &NormalizedPack) -> fs_matdb::NormalizedMaterialCardPack {
    fs_matdb::NormalizedMaterialCardPack::new(
        fs_matdb::MaterialStateId {
            chemistry: pack.pack_id().into(),
            phase: "solid".into(),
            process: "exact source product, construction and property-specific conditions".into(),
            revision: 0,
        },
        pack.clone(),
    )
    .unwrap()
}

#[test]
fn kerto_qp_exact_cross_veneer_layups_preserve_direction_and_statistics() {
    for (slug, thick, longitudinal, axial_mean, axial_fifth, flat_mean, flat_fifth) in [
        (
            "metsa-kerto-qp-39mm-2023",
            0.039,
            11.0,
            11700e6,
            9800e6,
            11300e6,
            9500e6,
        ),
        (
            "metsa-kerto-qp-75mm-2023",
            0.075,
            23.0,
            12300e6,
            10300e6,
            11400e6,
            9600e6,
        ),
    ] {
        let pack = compile(slug, 10);
        let base = thickness(
            flags(&[
                ("source-product-kerto-qp-2023", 1.0),
                ("source-longitudinal-veneers", longitudinal),
                ("source-cross-veneers", 2.0),
                ("source-temperature-known", 0.0),
                ("source-moisture-condition-known", 0.0),
            ]),
            thick,
        );
        let parallel = base
            .clone()
            .with_selector("load-parallel-to-face-veneer", 1.0)
            .unwrap();
        sample(&pack, "density", RHO, &base, 510.0);
        sample(&pack, "density-fifth-percentile", RHO, &base, 480.0);
        let e = sample(&pack, "young-modulus", PA, &parallel, axial_mean);
        assert_eq!(e.evidence.value.uncertainty, UncertaintyModel::Unstated);
        for (property, expected) in [
            ("longitudinal-modulus-fifth-percentile", axial_fifth),
            ("bending-modulus-flatwise-parallel-mean", flat_mean),
            (
                "bending-modulus-flatwise-parallel-fifth-percentile",
                flat_fifth,
            ),
            ("shear-modulus-edgewise-mean", 600e6),
            ("shear-modulus-flatwise-parallel-mean", 120e6),
            ("shear-modulus-edgewise-fifth-percentile", 400e6),
            ("shear-modulus-flatwise-parallel-fifth-percentile", 100e6),
        ] {
            sample(&pack, property, PA, &base, expected);
        }
        refuses(&pack, "young-modulus", PA, &base);
        refuses(
            &pack,
            "young-modulus",
            PA,
            &base
                .with_selector("load-parallel-to-face-veneer", 0.0)
                .unwrap(),
        );
        for bad in [
            thickness(parallel.clone(), thick + 0.001),
            parallel
                .clone()
                .with_selector("source-longitudinal-veneers", longitudinal + 1.0)
                .unwrap(),
            parallel
                .clone()
                .with_selector("source-cross-veneers", 0.0)
                .unwrap(),
            parallel
                .clone()
                .with_selector("source-moisture-condition-known", 1.0)
                .unwrap(),
        ] {
            refuses(&pack, "young-modulus", PA, &bad);
            refuses(&pack, "density", RHO, &bad);
        }
        refuses(&pack, "poisson-ratio", Dims::NONE, &parallel);
        refuses(&pack, "thermal-conductivity", K, &parallel);

        use fs_material::axial_reference::{
            AxialLoadingMode, AxialReferenceRequest, ResolvedAxialModulusReference,
        };
        use fs_material::state_point::{
            MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
        };
        let normalized = reference_card(&pack);
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
        let resolved = ResolvedAxialModulusReference::resolve(
            normalized.card(),
            &parallel,
            &request,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        // Authored small-strain reference using the declared mean Et,0.
        // This is neither a Eurocode capacity check nor a measured elastic range.
        let area = thick * 0.1;
        let mut coupon = resolved.coupon(1.0, area).unwrap();
        let response = coupon.advance_displacement(0.0001, 1.0).unwrap();
        let force = axial_mean * area * 0.0001;
        assert!((response.reference_force_n / force - 1.0).abs() < 1e-13);
        assert!((response.mass_kg / (510.0 * area) - 1.0).abs() < 1e-13);
        assert!((response.model_loading_work_j / (0.5 * force * 0.0001) - 1.0).abs() < 1e-13);
        assert!(coupon.advance_displacement(0.0003, 1.0).is_err());
        for property in resolved.source().properties() {
            normalized
                .card()
                .claims()
                .verify_receipt(&property.answer().receipt)
                .unwrap();
        }
    }
}

fn frozen_slab_reference(pack: &NormalizedPack, at: &QueryPoint, thick: f64, k: f64) {
    let normalized = reference_card(pack);
    let slab = fs_conduction::ThermalResistanceTerm::slab_from_card(
        "exact nominal-thickness source-reference slab, unit area",
        thick,
        1.0,
        normalized.card(),
        &key("thermal-conductivity", K),
        at,
        SelectionPolicy::SingleClaimOnly,
    )
    .unwrap();
    assert!((slab.value_k_per_w() / (thick / k) - 1.0).abs() < 1e-13);
}

#[test]
fn p5_particleboards_keep_wet_tests_and_thermal_density_basis_separate() {
    for (slug, thick, bending, strength, dry_bond, boil_bond) in [
        (
            "egger-p5-particleboard-18mm-2026",
            0.018,
            2400e6,
            16e6,
            0.45e6,
            0.14e6,
        ),
        (
            "egger-p5-particleboard-22mm-2026",
            0.022,
            2150e6,
            14e6,
            0.40e6,
            0.12e6,
        ),
    ] {
        let pack = compile(slug, 9);
        let base = thickness(
            flags(&[
                ("source-product-egger-p5-rec225-2026", 1.0),
                ("source-temperature-known", 0.0),
                ("source-moisture-condition-known", 0.0),
            ]),
            thick,
        );
        let unconditioned = base
            .clone()
            .with_selector("after-en1087-1-boil", 0.0)
            .unwrap();
        for (property, expected) in [
            ("bending-modulus-declared-board-mean", bending),
            ("bending-strength-declared-board-mean", strength),
            ("internal-bond-strength-declared-board-mean", dry_bond),
            ("surface-soundness-declared-board-mean", 0.95e6),
        ] {
            sample(&pack, property, PA, &unconditioned, expected);
            refuses(&pack, property, PA, &base);
            refuses(
                &pack,
                property,
                PA,
                &thickness(unconditioned.clone(), thick + 0.001),
            );
        }
        let boiled = base
            .clone()
            .with_selector("after-en1087-1-boil", 1.0)
            .unwrap();
        sample(
            &pack,
            "internal-bond-strength-declared-board-mean",
            PA,
            &boiled,
            boil_bond,
        );
        refuses(&pack, "bending-modulus-declared-board-mean", PA, &boiled);
        let swelling = base
            .clone()
            .with_selector("en317-immersion", 1.0)
            .unwrap()
            .with_quantity(
                "immersion-duration",
                QuantitySpec::dimensional(Dims([0, 0, 1, 0, 0, 0])),
                86400.0,
            )
            .unwrap();
        sample(
            &pack,
            "thickness-swelling-24h-declared-board-mean",
            Dims::NONE,
            &swelling,
            0.10,
        );
        refuses(
            &pack,
            "thickness-swelling-24h-declared-board-mean",
            Dims::NONE,
            &swelling
                .with_quantity(
                    "immersion-duration",
                    QuantitySpec::dimensional(Dims([0, 0, 1, 0, 0, 0])),
                    3600.0,
                )
                .unwrap(),
        );
        let thermal = base
            .clone()
            .with_quantity(
                "source-bulk-density-reference",
                QuantitySpec::dimensional(RHO),
                600.0,
            )
            .unwrap();
        sample(&pack, "thermal-conductivity", K, &thermal, 0.12);
        frozen_slab_reference(&pack, &thermal, thick, 0.12);
        refuses(&pack, "thermal-conductivity", K, &base);
        refuses(
            &pack,
            "thermal-conductivity",
            K,
            &thermal
                .clone()
                .with_quantity(
                    "source-bulk-density-reference",
                    QuantitySpec::dimensional(RHO),
                    650.0,
                )
                .unwrap(),
        );
        for (wet, resistance) in [(0.0, 50.0), (1.0, 15.0)] {
            sample(
                &pack,
                "water-vapour-diffusion-resistance-factor",
                Dims::NONE,
                &thermal
                    .clone()
                    .with_selector("source-vapour-condition-wet", wet)
                    .unwrap(),
                resistance,
            );
        }
        refuses(
            &pack,
            "water-vapour-diffusion-resistance-factor",
            Dims::NONE,
            &thermal,
        );
        // The 600 kg/m3 basis must not silently become product density.
        refuses(&pack, "density", RHO, &thermal);
        refuses(&pack, "specific-heat-capacity", CP, &thermal);
        refuses(&pack, "young-modulus", PA, &unconditioned);
        refuses(&pack, "poisson-ratio", Dims::NONE, &unconditioned);
    }
}

#[test]
fn glulam_building_element_requires_joint_direction_for_thermal_use() {
    let pack = compile("mm-hbe-gl24h-100mm-2024", 4);
    let base = thickness(
        flags(&[
            ("source-product-mm-hbe-gl24h-2024", 1.0),
            ("source-adhesive-muf", 1.0),
            ("source-temperature-known", 0.0),
            ("source-moisture-test-condition-known", 0.0),
        ]),
        0.1,
    );
    sample(&pack, "density", RHO, &base, 430.0);
    sample(
        &pack,
        "moisture-movement-transverse-average",
        Dims::NONE,
        &base,
        0.24,
    );
    sample(
        &pack,
        "moisture-movement-longitudinal-average",
        Dims::NONE,
        &base,
        0.01,
    );
    let parallel = base
        .clone()
        .with_selector("heat-flow-parallel-to-bonded-joints", 1.0)
        .unwrap();
    sample(&pack, "thermal-conductivity", K, &parallel, 0.13);
    frozen_slab_reference(&pack, &parallel, 0.1, 0.13);
    refuses(&pack, "thermal-conductivity", K, &base);
    refuses(
        &pack,
        "thermal-conductivity",
        K,
        &base
            .clone()
            .with_selector("heat-flow-parallel-to-bonded-joints", 0.0)
            .unwrap(),
    );
    for bad in [
        thickness(parallel.clone(), 0.12),
        parallel
            .clone()
            .with_selector("source-adhesive-muf", 0.0)
            .unwrap(),
        parallel
            .clone()
            .with_selector("source-moisture-test-condition-known", 1.0)
            .unwrap(),
    ] {
        refuses(&pack, "thermal-conductivity", K, &bad);
        refuses(&pack, "density", RHO, &bad);
    }
    refuses(&pack, "specific-heat-capacity", CP, &parallel);
    refuses(&pack, "young-modulus", PA, &parallel);
    // Published average movements have neither a strain datum nor a tested
    // moisture path. They do not satisfy the integrated hygroelastic contract.
    for property in fs_material::state_point::ORTHOTROPIC_HYGROSCOPIC_EXPANSION_PROPERTIES {
        refuses(&pack, property, Dims::NONE, &parallel);
    }
}

#[test]
fn aircrete_dry_and_service_condition_conductivities_are_discrete() {
    let pack = compile("hh-celcon-standard-100mm-bba2026", 12);
    let base = thickness(flags(&[("source-grade-celcon-standard-bba2026", 1.0)]), 0.1);
    let dry = base
        .clone()
        .with_selector("source-dry-reference", 1.0)
        .unwrap()
        .with_selector("source-density-temperature-known", 0.0)
        .unwrap();
    for (property, rho) in [
        ("density", 600.0),
        ("dry-density-range-lower-bound", 550.0),
        ("dry-density-range-upper-bound", 650.0),
    ] {
        sample(&pack, property, RHO, &dry, rho);
        refuses(
            &pack,
            property,
            RHO,
            &dry.clone()
                .with_selector("source-dry-reference", 0.0)
                .unwrap(),
        );
    }
    let dry_thermal = temperature(
        dry.clone()
            .with_selector("source-thermal-condition", 0.0)
            .unwrap()
            .with_selector("source-conductivity-temperature-known", 1.0)
            .unwrap(),
        "temperature",
        283.15,
    );
    sample(&pack, "thermal-conductivity", K, &dry_thermal, 0.13);
    frozen_slab_reference(&pack, &dry_thermal, 0.1, 0.13);
    refuses(
        &pack,
        "thermal-conductivity",
        K,
        &temperature(dry_thermal.clone(), "temperature", 293.15),
    );
    refuses(&pack, "thermal-conductivity", K, &base);
    for (condition, k) in [(1.0, 0.15), (2.0, 0.16), (3.0, 0.24)] {
        let at = base
            .clone()
            .with_selector("source-thermal-condition", condition)
            .unwrap()
            .with_selector("source-dry-reference", 0.0)
            .unwrap()
            .with_selector("source-conductivity-temperature-known", 0.0)
            .unwrap();
        sample(&pack, "thermal-conductivity", K, &at, k);
        frozen_slab_reference(&pack, &at, 0.1, k);
        refuses(
            &pack,
            "thermal-conductivity",
            K,
            &at.clone()
                .with_selector("source-thermal-condition", condition + 0.5)
                .unwrap(),
        );
        refuses(
            &pack,
            "thermal-conductivity",
            K,
            &at.clone()
                .with_selector("source-dry-reference", 1.0)
                .unwrap(),
        );
        // No dry nominal density is silently reused as an in-service density.
        refuses(
            &pack,
            "density",
            RHO,
            &at.with_selector("source-density-temperature-known", 0.0)
                .unwrap(),
        );
    }
    let unknown_tests = base
        .clone()
        .with_selector("source-test-temperature-known", 0.0)
        .unwrap()
        .with_selector("source-test-moisture-known", 0.0)
        .unwrap();
    let intended = unknown_tests
        .clone()
        .with_selector("source-intended-loading-orientation", 1.0)
        .unwrap();
    for (property, requirement) in [
        ("compressive-strength-mean-requirement", 3.6e6),
        ("compressive-strength-individual-minimum-requirement", 2.9e6),
    ] {
        sample(&pack, property, PA, &intended, requirement);
        refuses(&pack, property, PA, &unknown_tests);
        refuses(
            &pack,
            property,
            PA,
            &intended
                .clone()
                .with_selector("source-intended-loading-orientation", 0.0)
                .unwrap(),
        );
    }
    for (into, factor) in [(1.0, 5.0), (0.0, 10.0)] {
        sample(
            &pack,
            "water-vapour-diffusion-resistance-factor",
            Dims::NONE,
            &unknown_tests
                .clone()
                .with_selector("vapour-diffusion-into-block", into)
                .unwrap(),
            factor,
        );
    }
    refuses(
        &pack,
        "water-vapour-diffusion-resistance-factor",
        Dims::NONE,
        &unknown_tests,
    );
    sample(
        &pack,
        "moisture-movement-maximum",
        Dims::NONE,
        &unknown_tests,
        0.0004,
    );
    refuses(&pack, "compression-modulus", PA, &intended);
    refuses(&pack, "specific-heat-capacity", CP, &dry_thermal);
}

#[test]
fn cellular_glass_keeps_declared_bounds_separate_from_sample_statistics() {
    let pack = compile("foamglas-t4-plus-100mm-2025", 13);
    let base = thickness(
        flags(&[
            ("source-product-foamglas-t4-plus-2025", 1.0),
            ("source-temperature-known", 0.0),
            ("source-moisture-condition-known", 0.0),
        ]),
        0.1,
    );
    sample(&pack, "density", RHO, &base, 110.0);
    sample(
        &pack,
        "density-declared-relative-tolerance",
        Dims::NONE,
        &base,
        0.15,
    );
    let thickness_tolerance = sample(
        &pack,
        "nominal-thickness-declared-tolerance",
        Dims([1, 0, 0, 0, 0, 0]),
        &base,
        0.002,
    );
    sample(
        &pack,
        "thermal-resistance-area-specific-declared",
        Dims([0, -1, 3, 1, 0, 0]),
        &base,
        2.40,
    );
    let bound_key = PropertyKey::with_quantity(
        "thermal-conductivity-declared-upper-bound",
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::ThermalConductivity,
            ValueForm::Static,
        )),
    );
    let upper = pack
        .claims()
        .query_typed(&bound_key, &base, SelectionPolicy::SingleClaimOnly)
        .unwrap();
    assert_eq!(upper.evidence.value.value, 0.041);
    assert_eq!(upper.evidence.value.dims, K);
    assert!(upper.receipt.observation_backed);
    assert!(!upper.receipt.source_hashes.is_empty());
    pack.claims().verify_receipt(&upper.receipt).unwrap();
    for bad in [
        thickness(base.clone(), 0.12),
        base.clone()
            .with_selector("source-temperature-known", 1.0)
            .unwrap(),
    ] {
        assert!(
            pack.claims()
                .query_typed(&bound_key, &bad, SelectionPolicy::SingleClaimOnly)
                .is_err()
        );
    }
    // Explicit lower-bound calculation: L >= nominal minus its declared
    // tolerance, k <= its declared limit, and caller-fixed unit area.
    // The scalar slab owner receives a declared engineering input, never a
    // fabricated nominal-conductivity material receipt or a transient law.
    let minimum_thickness = 0.1 - thickness_tolerance.evidence.value.value;
    let lower_resistance = fs_conduction::ThermalResistanceTerm::slab(
        "cellular-glass material-only resistance lower bound, unit area",
        minimum_thickness,
        upper.evidence.value.value,
        1.0,
        fs_conduction::ResistanceUncertainty::Unstated,
        format!(
            "Authored monotonic bound L_min/(k_max*A), using declared thickness tolerance and conductivity upper limit from {}; frozen source conditions, no joints or finishes, not actual nominal conductivity",
            pack.content_hash()
        ),
    )
    .unwrap();
    assert!((lower_resistance.value_k_per_w() / (0.098 / 0.041) - 1.0).abs() < 1e-13);
    assert!(matches!(
        lower_resistance.origin(),
        fs_conduction::ResistanceOrigin::Declared { .. }
    ));

    for (property, limit) in [
        ("compressive-strength-declared-minimum", 600e3),
        ("bending-strength-declared-minimum", 450e3),
        ("tensile-strength-declared-minimum", 150e3),
    ] {
        sample(&pack, property, PA, &base, limit);
    }
    let statistics = base
        .clone()
        .with_selector("source-sample-size", 50.0)
        .unwrap();
    let mean = sample(
        &pack,
        "compressive-stress-sample-mean",
        PA,
        &statistics,
        750e3,
    );
    assert_eq!(mean.evidence.value.uncertainty, UncertaintyModel::Unstated);
    sample(
        &pack,
        "compressive-stress-sample-standard-deviation",
        PA,
        &statistics,
        55e3,
    );
    refuses(&pack, "compressive-stress-sample-mean", PA, &base);
    let characteristic = statistics
        .clone()
        .with_selector("source-characteristic-fractile", 0.05)
        .unwrap()
        .with_selector("source-one-sided-confidence", 0.75)
        .unwrap();
    sample(
        &pack,
        "compressive-stress-characteristic",
        PA,
        &characteristic,
        633e3,
    );
    refuses(&pack, "compressive-stress-characteristic", PA, &statistics);
    for bad in [
        characteristic
            .clone()
            .with_selector("source-characteristic-fractile", 0.50)
            .unwrap(),
        characteristic
            .clone()
            .with_selector("source-one-sided-confidence", 0.95)
            .unwrap(),
        characteristic
            .with_selector("source-sample-size", 25.0)
            .unwrap(),
    ] {
        refuses(&pack, "compressive-stress-characteristic", PA, &bad);
    }
    sample(&pack, "specific-heat-capacity", CP, &base, 1000.0);
    let expansion = base
        .clone()
        .with_selector("source-expansion-reference-temperature-known", 0.0)
        .unwrap()
        .with_selector("source-expansion-interval-known", 0.0)
        .unwrap();
    sample(
        &pack,
        "thermal-expansion-coefficient-reference",
        Dims([0, 0, 0, -1, 0, 0]),
        &expansion,
        9e-6,
    );
    refuses(
        &pack,
        "thermal-expansion-coefficient-reference",
        Dims([0, 0, 0, -1, 0, 0]),
        &base,
    );
    refuses(&pack, "thermal-conductivity", K, &base);
    refuses(&pack, "compression-modulus", PA, &base);
    refuses(&pack, "young-modulus", PA, &base);
}
