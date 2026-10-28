//! G0/G1/G3 source-to-pack/query and reduced thermal-consumer checks.
//! These establish source transport and bounded model evaluation, not experimental validation.
#![deny(unsafe_code)]

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use fs_matdb::{
    MaterialAnswer, NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

const RHO: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const PA: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const CP: Dims = Dims([2, 0, -2, -1, 0, 0]);
const K: Dims = Dims([1, 1, -3, -1, 0, 0]);
const TEMP: Dims = Dims([0, 0, 0, 1, 0, 0]);
const SPEED: Dims = Dims([1, 0, -1, 0, 0, 0]);

fn workspace() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .to_path_buf()
}

fn compile(slug: &str, count: usize) -> NormalizedPack {
    static NEXT: AtomicU64 = AtomicU64::new(0);
    let time = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "frankensim-material-expansion-{}-{time}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir(&dir).unwrap();
    let output = dir.join("material.fsmatpk");
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(workspace().join(format!("data/matdb/seed-v1/{slug}/manifest.tsv")))
        .arg("--out")
        .arg(&output)
        .env("CARGO_WORKSPACE_DIR", workspace())
        .current_dir(workspace())
        .output()
        .expect("run the actual material compiler");
    assert!(
        run.status.success(),
        "{slug}: {}",
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(output).unwrap();
    let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes, "{slug}: canonical round trip");
    assert_eq!(pack.pack_id(), slug);
    assert_eq!(pack.claims().claims_ordered().count(), count, "{slug}");
    pack
}

fn flags(values: &[(&str, f64)]) -> QueryPoint {
    values.iter().fold(QueryPoint::new(), |p, (name, value)| {
        p.with_quantity(*name, QuantitySpec::dimensional(Dims::NONE), *value)
            .unwrap()
    })
}

fn temperature(point: QueryPoint, name: &str, kelvin: f64) -> QueryPoint {
    point
        .with_quantity(
            name,
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            kelvin,
        )
        .unwrap()
}

fn key(name: &str, dims: Dims) -> PropertyKey {
    let quantity = if matches!(
        name,
        "thermal-conductivity" | "thermal-conductivity-transverse"
    ) {
        assert_eq!(dims, K);
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::ThermalConductivity,
            ValueForm::Static,
        ))
    } else {
        QuantitySpec::dimensional(dims)
    };
    PropertyKey::with_quantity(name, quantity)
}

fn sample(
    pack: &NormalizedPack,
    name: &str,
    dims: Dims,
    point: &QueryPoint,
    expected: f64,
) -> MaterialAnswer {
    let answer = pack
        .claims()
        .query_typed(&key(name, dims), point, SelectionPolicy::SingleClaimOnly)
        .unwrap_or_else(|e| panic!("{} {name}: {e:?}", pack.pack_id()));
    let actual = answer.evidence.value.value;
    assert!(
        (actual - expected).abs() <= 1e-11 * expected.abs().max(1.0),
        "{} {name}: {actual} != {expected}",
        pack.pack_id()
    );
    assert_eq!(answer.evidence.value.dims, dims);
    assert!(answer.receipt.observation_backed);
    assert!(!answer.receipt.source_hashes.is_empty());
    assert!(!answer.receipt.axis_quantities.is_empty());
    pack.claims().verify_receipt(&answer.receipt).unwrap();
    answer
}

fn refuses(pack: &NormalizedPack, name: &str, dims: Dims, point: &QueryPoint) {
    assert!(
        pack.claims()
            .query_typed(&key(name, dims), point, SelectionPolicy::SingleClaimOnly)
            .is_err(),
        "{} {name} must reject {point:?}",
        pack.pack_id()
    );
}

#[test]
fn g0_g3_clear_wood_density_stiffness_and_moisture_are_preserved() {
    let mut final_temperatures = Vec::new();
    for (species, count, sg, bending_mpa, rupture_kpa) in [
        ("pine-loblolly", 25, 0.51, 12300.0, 88000.0),
        ("pine-longleaf", 25, 0.59, 13700.0, 100000.0),
        ("pine-slash", 23, 0.59, 13700.0, 112000.0),
        ("pine-ponderosa", 25, 0.40, 8900.0, 65000.0),
        ("oak-northern-red", 16, 0.63, 12500.0, 99000.0),
        ("oak-white", 16, 0.68, 12300.0, 105000.0),
        ("cedar-northern-white", 25, 0.31, 5500.0, 45000.0),
        ("fir-subalpine", 24, 0.32, 8900.0, 59000.0),
    ] {
        let pack = compile(&format!("{species}-fpl-gtr282"), count);
        let p = flags(&[
            ("moisture-content", 0.12),
            ("clear-specimen", 1.0),
            ("source-test-temperature-known", 0.0),
        ]);
        sample(&pack, "density", RHO, &p, sg * 1000.0 * 1.12);
        sample(
            &pack,
            "modulus_of_elasticity_bending",
            PA,
            &p,
            bending_mpa * 1e6,
        );
        sample(
            &pack,
            "young_modulus_longitudinal",
            PA,
            &p,
            bending_mpa * 1.1e6,
        );
        sample(&pack, "bending_mor", PA, &p, rupture_kpa * 1000.0);
        for bad in [
            p.clone().with("moisture-content", 0.20).unwrap(),
            p.clone().with("clear-specimen", 0.0).unwrap(),
            QueryPoint::new(),
        ] {
            refuses(&pack, "density", RHO, &bad);
        }
        refuses(&pack, "young_modulus", PA, &p); // no isotropic substitution
        if species.starts_with("oak-") {
            refuses(&pack, "er_over_el", Dims::NONE, &p);
        }

        // FPL Chapter 4 engineering correlations complete the same species'
        // thermal inputs. They are estimates at fixed moisture, not added
        // specimen measurements or permission to treat wood as isotropic.
        let thermal = p
            .clone()
            .with_quantity(
                "engineering-reference-model",
                QuantitySpec::dimensional(Dims::NONE),
                1.0,
            )
            .unwrap()
            .with_quantity(
                "source-pressure-known",
                QuantitySpec::dimensional(Dims::NONE),
                0.0,
            )
            .unwrap();
        for t in [280.0, 297.15, 350.0, 420.0] {
            let cp = 1000.0
                * ((0.1031 + 0.003867 * t + 4.18 * 0.12) / 1.12
                    + 12.0 * (-0.06191 + 0.000236 * t - 0.000133 * 12.0));
            sample(
                &pack,
                "specific-heat-capacity",
                CP,
                &temperature(thermal.clone(), "temperature", t),
                cp,
            );
        }
        let at = temperature(thermal, "temperature", 297.15);
        sample(
            &pack,
            "thermal-conductivity-transverse",
            K,
            &at,
            sg * (0.1941 + 0.004064 * 12.0) + 0.01864,
        );
        // No isotropic/longitudinal coefficient or conductivity-temperature
        // extrapolation is manufactured to complete a solver requirement set.
        refuses(&pack, "thermal-conductivity", K, &at);
        for bad in [
            at.clone().with("temperature", 297.16).unwrap(),
            at.clone().with("engineering-reference-model", 0.0).unwrap(),
            at.clone().with("moisture-content", 0.10).unwrap(),
            at.clone().with("source-pressure-known", 1.0).unwrap(),
        ] {
            refuses(&pack, "thermal-conductivity-transverse", K, &bad);
        }
        for t in [279.99, 420.01] {
            refuses(
                &pack,
                "specific-heat-capacity",
                CP,
                &at.clone().with("temperature", t).unwrap(),
            );
        }
        final_temperatures.push(catalog_wood_sensible_heat(&pack, &at));
    }
    // Equal reference volumes and deposited heat: the denser white oak heats
    // less than northern white cedar. No wood conductivity is extrapolated.
    assert!(final_temperatures[6] > final_temperatures[5] + 15.0);
}

/// Fixed-mass sensible heating uses the new Cp(T) in the production enthalpy
/// law. Moisture is fixed and the initial mass uses the species' reference
/// density; this does not infer thermal expansion, drying or internal gradients.
fn catalog_wood_sensible_heat(pack: &NormalizedPack, point: &QueryPoint) -> f64 {
    use fs_matdb::{MaterialStateId, NormalizedMaterialCardPack};
    use fs_material::phase::{
        HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase,
    };

    let normalized = NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: pack.pack_id().into(),
            phase: "solid clear wood at fixed 12 percent dry-basis moisture".into(),
            process: "FPL engineering thermal correlation; fixed reference mass".into(),
            revision: 0,
        },
        pack.clone(),
    )
    .unwrap();
    let card = normalized.card();
    let lower = point.clone().with("temperature", 280.0).unwrap();
    let upper = point.clone().with("temperature", 420.0).unwrap();
    let resolved = ResolvedHeatCapacityCurve::try_new(
        card,
        HeatCapacitySourceQuery {
            heat_capacity: &key("specific-heat-capacity", CP),
            density: &key("density", RHO),
            temperature_axis: "temperature",
            lower: &lower,
            upper: &upper,
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        0.01,
    )
    .unwrap();
    for receipt in resolved.receipts() {
        card.claims().verify_receipt(receipt).unwrap();
    }
    let curve = resolved.curve();
    let mut state = curve.state_at_specific_enthalpy(0.0).unwrap();
    let mass_kg = state.bulk_density_kg_m3() * 0.0001;
    // Independent integration/inversion of FPL 4-18 through 4-20 at MC=12%.
    let cp0 = 1000.0
        * ((0.1031 + 0.003867 * 280.0 + 4.18 * 0.12) / 1.12
            + 12.0 * (-0.06191 + 0.000236 * 280.0 - 0.000133 * 12.0));
    let slope = 1000.0 * (0.003867 / 1.12 + 12.0 * 0.000236);
    for step in 1..=20 {
        state = curve
            .advance_specific_energy(state, 125.0 / mass_kg)
            .unwrap();
        let heat = f64::from(step) * 125.0 / mass_kg;
        let expected_t = 280.0 + 2.0 * heat / (cp0 + (cp0 * cp0 + 2.0 * slope * heat).sqrt());
        assert!((state.temperature_k() - expected_t).abs() < 0.01 / cp0 + 1e-8);
        assert!((mass_kg * state.specific_enthalpy_j_kg() - f64::from(step) * 125.0).abs() < 1e-9);
        assert_eq!(state.material_card_identity(), card.content_hash());
        assert_eq!(state.phase(), SolidLiquidPhase::Solid);
    }
    // Using a single initial Cp would give a different thermal response.
    assert!((280.0 + 2500.0 / mass_kg / cp0 - state.temperature_k()).abs() > 0.5);
    assert!(curve.advance_specific_energy(state, 1e6).is_err());
    eprintln!(
        "{}: 2500 J into fixed reference volume = {:.6} K",
        pack.pack_id(),
        state.temperature_k()
    );
    state.temperature_k()
}

#[test]
fn g0_g3_elements_minerals_and_alloys_keep_reference_model_guards() {
    let mut cooling_times = Vec::new();
    for (slug, count, t, rho, cp, k) in [
        (
            "diamond-crystran-thermal-reference",
            3,
            273.0,
            3510.0,
            502.0,
            2600.0,
        ),
        (
            "fluorite-caf2-crystran-reference",
            3,
            293.0,
            3180.0,
            854.0,
            9.71,
        ),
        (
            "fused-silica-crystran-reference",
            5,
            298.15,
            2203.0,
            703.0,
            1.38,
        ),
        (
            "germanium-crystran-thermal-reference",
            3,
            293.0,
            5330.0,
            310.0,
            58.62,
        ),
        (
            "halite-nacl-crystran-reference",
            4,
            273.0,
            2170.0,
            854.0,
            6.5,
        ),
        (
            "periclase-mgo-crystran-reference",
            4,
            273.0,
            3580.0,
            877.0,
            42.0,
        ),
        (
            "aluminum-5005-austral-wright-reference",
            7,
            293.15,
            2700.0,
            900.0,
            201.0,
        ),
        (
            "aluminum-5052-austral-wright-reference",
            7,
            293.15,
            2680.0,
            900.0,
            138.0,
        ),
        (
            "bronze-c93200-cda-reference",
            4,
            293.15,
            8912.929316685405,
            376.812,
            58.15268479007874,
        ),
        (
            "bronze-c95400-cda-reference",
            4,
            293.15,
            7445.89436704464,
            418.68,
            58.67190518999016,
        ),
        (
            "cupronickel-c70600-cda-reference",
            6,
            293.15,
            8940.609221395608,
            376.812,
            44.99910132565617,
        ),
        (
            "magnesium-espi-reference",
            4,
            298.15,
            1738.0,
            1016.712,
            156.0,
        ),
        (
            "molybdenum-attl-reference",
            5,
            293.15,
            10220.0,
            240.0,
            140.0,
        ),
        ("tantalum-efinea-reference", 4, 293.15, 16650.0, 140.0, 57.5),
    ] {
        let pack = compile(slug, count);
        let p = temperature(
            flags(&[
                ("engineering-reference-model", 1.0),
                ("source-temperatures-matched", 0.0),
            ]),
            "temperature",
            t,
        );
        sample(&pack, "density", RHO, &p, rho);
        sample(&pack, "specific-heat-capacity", CP, &p, cp);
        sample(&pack, "thermal-conductivity", K, &p, k);
        cooling_times.push(catalog_slab_cooling(&pack, &p, t, rho, cp, k));
        for bad in [
            p.clone().with("engineering-reference-model", 0.0).unwrap(),
            p.clone().with("source-temperatures-matched", 1.0).unwrap(),
            p.clone().with("temperature", t - 0.1).unwrap(),
            QueryPoint::new(),
        ] {
            refuses(&pack, "density", RHO, &bad);
        }
        if count == 3 {
            refuses(&pack, "young-modulus", PA, &p);
        }
        if slug == "aluminum-5052-austral-wright-reference" {
            sample(&pack, "young-modulus", PA, &p, 69.3e9);
            sample(
                &pack,
                "electrical-resistivity",
                Dims([3, 1, -3, 0, -2, 0]),
                &p,
                5e-8,
            );
        }
    }
    // A real material change reaches the physical screening calculation.
    // The thermal-reference source assumptions remain explicit for both.
    assert!(cooling_times[2] > 1000.0 * cooling_times[0]); // fused silica / diamond
}

/// Compile-derived material data feed the production constant-capacity binding
/// and the existing first-mode slab screen. The source bundles explicitly
/// declare a frozen 1 K reference interval; this is not a broad heating law.
fn catalog_slab_cooling(
    pack: &NormalizedPack,
    point: &QueryPoint,
    low_k: f64,
    rho: f64,
    cp: f64,
    k: f64,
) -> f64 {
    use fs_conduction::reduced_slab::{MaterialFirstModeSlabCooling, MaterialFirstModeSlabInput};
    use fs_conduction::transient::{MaterialHeatCapacity, MaterialHeatCapacitySource};
    use fs_matdb::{MaterialStateId, NormalizedMaterialCardPack};
    use fs_material::state_point::MaterialPropertySelection;

    let normalized = NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: pack.pack_id().into(),
            phase: "solid".into(),
            process: "explicit frozen engineering reference; source temperatures unmatched".into(),
            revision: 0,
        },
        pack.clone(),
    )
    .unwrap();
    let card = normalized.card();
    let rho_key = key("density", RHO);
    let cp_key = key("specific-heat-capacity", CP);
    let k_key = key("thermal-conductivity", K);
    let source = MaterialHeatCapacitySource {
        density: &rho_key,
        specific_heat: &cp_key,
        temperature_axis: "temperature",
        fixed_source_selectors: &["engineering-reference-model", "source-temperatures-matched"],
    };
    let capacity = MaterialHeatCapacity::from_card_at_query_points(
        card,
        source,
        point,
        [low_k, low_k + 1.0],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap_or_else(|error| panic!("{}: {error}", pack.pack_id()));
    assert!((capacity.value_j_per_m3_k() / (rho * cp) - 1.0).abs() < 1e-12);
    for state in [
        capacity.reference(),
        &capacity.endpoints()[0],
        &capacity.endpoints()[1],
    ] {
        assert_eq!(state.card_identity(), card.content_hash());
        for property in state.properties() {
            card.claims()
                .verify_receipt(&property.answer().receipt)
                .unwrap();
        }
    }
    let input = MaterialFirstModeSlabInput {
        thickness_m: 0.001,
        initial_temperature_k: low_k + 0.75,
        boundary_temperature_k: low_k + 0.25,
        threshold_temperature_k: low_k + 0.5,
    };
    let screen = MaterialFirstModeSlabCooling::from_card(
        &capacity,
        card,
        &k_key,
        SelectionPolicy::SingleClaimOnly,
        input,
    )
    .unwrap();
    let result = screen.step();
    assert_eq!(screen.source(), &capacity);
    for receipt in screen.conductivity_receipts() {
        card.claims().verify_receipt(receipt).unwrap();
    }
    let expected_tau = 1e-6 * rho * cp / (std::f64::consts::PI.powi(2) * k);
    assert!((result.time_constant_s / expected_tau - 1.0).abs() < 1e-12);
    assert!((result.time_to_threshold_s / (expected_tau * 2.0_f64.ln()) - 1.0).abs() < 1e-12);
    let crossing = screen.temperature_at_s(result.time_to_threshold_s).unwrap();
    assert!((crossing - input.threshold_temperature_k).abs() < 1e-10);
    let thick = MaterialFirstModeSlabCooling::from_card(
        &capacity,
        card,
        &k_key,
        SelectionPolicy::SingleClaimOnly,
        MaterialFirstModeSlabInput {
            thickness_m: 0.002,
            ..input
        },
    )
    .unwrap()
    .step();
    assert!((thick.time_to_threshold_s / result.time_to_threshold_s - 4.0).abs() < 1e-12);
    assert!(
        MaterialHeatCapacity::from_card_at_query_points(
            card,
            source,
            point,
            [low_k, low_k + 1.001],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .is_err()
    );
    assert!(
        MaterialHeatCapacity::from_card_at_query_points(
            card,
            MaterialHeatCapacitySource {
                fixed_source_selectors: &[],
                ..source
            },
            point,
            [low_k, low_k + 1.0],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .is_err()
    );
    eprintln!(
        "{}: first-mode 1 mm slab half-excess cooling = {:.9} s",
        pack.pack_id(),
        result.time_to_threshold_s
    );
    result.time_to_threshold_s
}

fn plastic_point(grade: &str) -> QueryPoint {
    let mut p = flags(&[
        (grade, 1.0),
        ("source-temperature-known", 0.0),
        ("source-moisture-condition-known", 0.0),
        ("source-pressure-known", 0.0),
    ]);
    if grade.contains("tecamid") {
        p = p
            .with_quantity(
                "immediately-after-machining",
                QuantitySpec::dimensional(Dims::NONE),
                1.0,
            )
            .unwrap();
    }
    p
}

#[test]
fn g0_g3_plastics_keep_grade_test_speed_and_discrete_compression() {
    for (grade, rho, cp, k, modulus) in [
        ("tecaform-ah-natural", 1410.0, 1400.0, 0.39, 2.8e9),
        ("tecamid-6-natural", 1140.0, 1600.0, 0.37, 3.3e9),
        ("tecapet-white", 1360.0, 1200.0, 0.31, 3.1e9),
        ("tecaflon-pvdf-natural", 1780.0, 1300.0, 0.25, 2.2e9),
        ("tecapeek-natural", 1310.0, 1100.0, 0.27, 4.2e9),
        ("tecapei-natural", 1280.0, 1200.0, 0.21, 3.2e9),
        ("tecason-s-natural", 1240.0, 1200.0, 0.21, 2.7e9),
        ("tecamid-66-natural", 1150.0, 1500.0, 0.36, 3.5e9),
    ] {
        let count = if matches!(grade, "tecapeek-natural" | "tecapei-natural") {
            11
        } else {
            10
        };
        let pack = compile(&format!("ensinger-{grade}"), count);
        let grade_axis = format!("source-grade-{grade}");
        let p = plastic_point(&grade_axis);
        sample(&pack, "density", RHO, &p, rho);
        sample(&pack, "specific-heat-capacity", CP, &p, cp);
        sample(&pack, "thermal-conductivity", K, &p, k);
        let tension = p
            .clone()
            .with_quantity(
                "test-speed",
                QuantitySpec::dimensional(SPEED),
                1.6666666666666667e-5,
            )
            .unwrap();
        sample(&pack, "tensile-modulus", PA, &tension, modulus);
        refuses(&pack, "tensile-modulus", PA, &p);
        refuses(
            &pack,
            "density",
            RHO,
            &p.clone().with(&grade_axis, 0.0).unwrap(),
        );
        if grade == "tecamid-66-natural" {
            let missing_machining = flags(&[
                (&grade_axis, 1.0),
                ("source-temperature-known", 0.0),
                ("source-moisture-condition-known", 0.0),
                ("source-pressure-known", 0.0),
            ]);
            // All eight original claims share the producer's condition, not
            // only the newly transcribed thermal-expansion intervals.
            for (name, dims, value, speed, strain) in [
                ("density", RHO, 1150.0, None, None),
                ("thermal-conductivity", K, 0.36, None, None),
                ("specific-heat-capacity", CP, 1500.0, None, None),
                (
                    "tensile-modulus",
                    PA,
                    3.5e9,
                    Some(1.6666666666666667e-5),
                    None,
                ),
                (
                    "tensile-strength",
                    PA,
                    85e6,
                    Some(8.333333333333334e-4),
                    None,
                ),
                (
                    "flexural-modulus",
                    PA,
                    3.1e9,
                    Some(3.3333333333333335e-5),
                    None,
                ),
                (
                    "compression-modulus",
                    PA,
                    2.7e9,
                    Some(8.333333333333333e-5),
                    None,
                ),
                (
                    "compressive-stress",
                    PA,
                    20e6,
                    Some(8.333333333333333e-5),
                    Some(0.01),
                ),
            ] {
                let test_conditions = |mut query: QueryPoint| {
                    if let Some(speed) = speed {
                        query = query
                            .with_quantity("test-speed", QuantitySpec::dimensional(SPEED), speed)
                            .unwrap();
                    }
                    if let Some(strain) = strain {
                        query = query
                            .with_quantity(
                                "compressive-strain",
                                QuantitySpec::dimensional(Dims::NONE),
                                strain,
                            )
                            .unwrap();
                    }
                    query
                };
                sample(&pack, name, dims, &test_conditions(p.clone()), value);
                for wrong in [
                    missing_machining.clone(),
                    p.clone().with("immediately-after-machining", 0.0).unwrap(),
                ] {
                    refuses(&pack, name, dims, &test_conditions(wrong));
                }
            }
        }
        if grade == "tecaform-ah-natural" {
            let compression = p
                .with_quantity(
                    "test-speed",
                    QuantitySpec::dimensional(SPEED),
                    8.333333333333333e-5,
                )
                .unwrap();
            for (strain, stress) in [(0.01, 20e6), (0.02, 35e6), (0.05, 68e6)] {
                let at = compression
                    .clone()
                    .with_quantity(
                        "compressive-strain",
                        QuantitySpec::dimensional(Dims::NONE),
                        strain,
                    )
                    .unwrap();
                sample(&pack, "compressive-stress", PA, &at, stress);
            }
            let between = compression
                .with_quantity(
                    "compressive-strain",
                    QuantitySpec::dimensional(Dims::NONE),
                    0.03,
                )
                .unwrap();
            refuses(&pack, "compressive-stress", PA, &between);
        }
    }
}

#[test]
fn g0_g3_four_polymer_longitudinal_intervals_drive_source_bound_length() {
    for (grade, expected) in [
        ("tecaform-ah-natural", [(13e-5, 0.00481), (14e-5, 0.01078)]),
        ("tecamid-6-natural", [(12e-5, 0.00444), (13e-5, 0.01001)]),
        ("tecapet-white", [(8e-5, 0.00296), (10e-5, 0.00770)]),
        (
            "tecaflon-pvdf-natural",
            [(16e-5, 0.00592), (18e-5, 0.01386)],
        ),
    ] {
        let pack = compile(&format!("ensinger-{grade}"), 10);
        let at = plastic_point(&format!("source-grade-{grade}"));
        catalog_longitudinal_expansion(&pack, &at, expected);
    }
}

#[test]
fn g0_g3_additional_polymer_intervals_keep_distinct_reference_lengths() {
    for (grade, expected, upper_interval) in [
        (
            "tecapeek-natural",
            [(5e-5, 0.00185), (5e-5, 0.00385)],
            Some((7e-5, 0.0035)),
        ),
        (
            "tecapei-natural",
            [(5e-5, 0.00185), (5e-5, 0.00385)],
            Some((6e-5, 0.003)),
        ),
        (
            "tecason-s-natural",
            [(6e-5, 0.00222), (6e-5, 0.00462)],
            None,
        ),
        (
            "tecamid-66-natural",
            [(11e-5, 0.00407), (12e-5, 0.00924)],
            None,
        ),
    ] {
        let pack = compile(
            &format!("ensinger-{grade}"),
            if upper_interval.is_some() { 11 } else { 10 },
        );
        let at = plastic_point(&format!("source-grade-{grade}"));
        catalog_longitudinal_expansion(&pack, &at, expected);
        if let Some((alpha, strain)) = upper_interval {
            catalog_longitudinal_intervals(
                &pack,
                &at,
                "mean-longitudinal-thermal-expansion-coefficient-from-100c",
                373.15,
                &[(423.15, alpha, strain)],
            );
            // The separate 100 C reference observation cannot extend the
            // 23 C reference law to 150 C or supply a stitched strain law.
            let unsupported = temperature(
                temperature(
                    at.clone().with("source-temperature-known", 1.0).unwrap(),
                    "source-reference-temperature",
                    296.15,
                ),
                "source-range-end-temperature",
                423.15,
            );
            refuses(
                &pack,
                "mean-longitudinal-thermal-expansion-coefficient-from-23c",
                Dims([0, 0, 0, -1, 0, 0]),
                &unsupported,
            );
        }
    }
}

/// Longitudinal interval means reach the original engineering-strain owner.
/// This supplies no transverse strain or continuous instantaneous-alpha law.
fn catalog_longitudinal_expansion(
    pack: &NormalizedPack,
    original: &QueryPoint,
    expected: [(f64, f64); 2],
) {
    catalog_longitudinal_intervals(
        pack,
        original,
        "mean-longitudinal-thermal-expansion-coefficient-from-23c",
        296.15,
        &[
            (333.15, expected[0].0, expected[0].1),
            (373.15, expected[1].0, expected[1].1),
        ],
    );
}

fn catalog_longitudinal_intervals(
    pack: &NormalizedPack,
    original: &QueryPoint,
    coefficient_name: &str,
    reference: f64,
    expected: &[(f64, f64, f64)],
) {
    use fs_matdb::{MaterialStateId, NormalizedMaterialCardPack};
    use fs_material::state_point::{
        MaterialPropertySelection, resolve_mean_thermal_expansion_strain,
    };
    let normalized = NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: pack.pack_id().into(),
            phase: "solid polymer stock shape".into(),
            process: "producer longitudinal interval means; original grade and conditioning".into(),
            revision: 0,
        },
        pack.clone(),
    )
    .unwrap();
    let card = normalized.card();
    let coefficient = key(coefficient_name, Dims([0, 0, 0, -1, 0, 0]));
    for &(endpoint, alpha, strain) in expected {
        let at = temperature(
            temperature(
                original
                    .clone()
                    .with("source-temperature-known", 1.0)
                    .unwrap(),
                "source-reference-temperature",
                reference,
            ),
            "source-range-end-temperature",
            endpoint,
        );
        let resolve = |p: &QueryPoint| {
            resolve_mean_thermal_expansion_strain(
                card,
                &coefficient,
                p,
                "source-reference-temperature",
                "source-range-end-temperature",
                MaterialPropertySelection::SingleClaimOnly,
            )
        };
        sample(pack, coefficient.name(), coefficient.dims(), &at, alpha);
        let result = resolve(&at).unwrap();
        assert!((result.engineering_strain() - strain).abs() < 1e-14);
        assert_eq!(result.reference_temperature_k(), reference);
        assert_eq!(result.endpoint_temperature_k(), endpoint);
        assert_eq!(result.resolved().card_identity(), card.content_hash());
        for property in result.resolved().properties() {
            card.claims()
                .verify_receipt(&property.answer().receipt)
                .unwrap();
        }
        // Free endpoint length of a longitudinal coupon measuring 250 mm at
        // THIS interval's reference temperature, including the 100 C cases.
        let length_m = 0.25 * (1.0 + result.engineering_strain());
        assert!((length_m - 0.25 - 0.25 * strain).abs() < 1e-15);
        for wrong in [
            at.clone()
                .with("source-range-end-temperature", endpoint - 1.0)
                .unwrap(),
            at.clone()
                .with("source-reference-temperature", 293.15)
                .unwrap(),
            at.clone()
                .with(
                    "source-reference-temperature",
                    if reference == 296.15 { 373.15 } else { 296.15 },
                )
                .unwrap(),
            at.clone().with("source-temperature-known", 0.0).unwrap(),
            at.clone().with("source-pressure-known", 1.0).unwrap(),
        ] {
            assert!(resolve(&wrong).is_err());
        }
        if pack.pack_id() == "ensinger-tecamid-66-natural" {
            let missing = temperature(
                temperature(
                    flags(&[
                        ("source-grade-tecamid-66-natural", 1.0),
                        ("source-temperature-known", 1.0),
                        ("source-moisture-condition-known", 0.0),
                        ("source-pressure-known", 0.0),
                    ]),
                    "source-reference-temperature",
                    reference,
                ),
                "source-range-end-temperature",
                endpoint,
            );
            assert!(resolve(&missing).is_err());
            assert!(
                resolve(&at.clone().with("immediately-after-machining", 0.0).unwrap()).is_err()
            );
        }
        refuses(
            pack,
            "linear_thermal_expansion_coefficient",
            coefficient.dims(),
            &at,
        );
        eprintln!(
            "{}: longitudinal interval {reference}–{endpoint} K gives strain {strain}, length {length_m} m",
            pack.pack_id()
        );
    }
}

#[test]
fn g0_g3_stainless_reference_properties_drive_only_the_supported_slab_state() {
    use fs_conduction::{SlabConductivityQuery, ThermalResistanceTerm};
    use fs_matdb::{MaterialStateId, NormalizedMaterialCardPack};
    let pack = compile("stainless-14401-thyssenkrupp-2017", 10);
    let at = temperature(
        flags(&[
            ("source-heat-treatment-known", 0.0),
            ("source-pressure-known", 0.0),
        ]),
        "temperature",
        293.15,
    );
    for (name, dims, value) in [
        ("density", RHO, 8000.0),
        ("specific-heat-capacity", CP, 500.0),
        ("thermal-conductivity", K, 15.0),
        ("electrical-resistivity", Dims([3, 1, -3, 0, -2, 0]), 7.5e-7),
    ] {
        sample(&pack, name, dims, &at, value);
        for wrong in [
            at.clone().with("temperature", 293.16).unwrap(),
            at.clone().with("source-heat-treatment-known", 1.0).unwrap(),
            at.clone().with("source-pressure-known", 1.0).unwrap(),
        ] {
            refuses(&pack, name, dims, &wrong);
        }
    }
    let normalized = NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: pack.pack_id().into(),
            phase: "solid".into(),
            process: "producer grade reference at exactly 20 C; heat treatment unknown".into(),
            revision: 0,
        },
        pack,
    )
    .unwrap();
    let conductivity = key("thermal-conductivity", K);
    let card = normalized.card();
    let term = ThermalResistanceTerm::slab_from_card(
        "reference steel wall",
        0.02,
        0.004,
        card,
        &conductivity,
        &at,
        SelectionPolicy::SingleClaimOnly,
    )
    .unwrap();
    assert!((term.value_k_per_w() - 1.0 / 3.0).abs() < 1e-14);
    // The neighboring E(T) curve does not widen these singleton thermal facts.
    assert!(
        ThermalResistanceTerm::slab_over_temperature_span_from_card(
            "unsupported warm wall",
            0.02,
            0.004,
            card,
            SlabConductivityQuery {
                property: &conductivity,
                temperature_axis: "temperature",
                lower: &at,
                upper: &at.clone().with("temperature", 294.15).unwrap(),
                selection: SelectionPolicy::SingleClaimOnly
            }
        )
        .is_err()
    );
}

#[test]
fn g0_g3_foam_grade_and_cured_resin_conditions_survive_compilation() {
    for (grade, rho, modulus) in [("51", 52.0, 43e6), ("110", 110.0, 120e6)] {
        let pack = compile(&format!("evonik-rohacell-{grade}-igf-2022"), 7);
        let grade_axis = format!("source-grade-rohacell-{grade}-igf-2022");
        let p = flags(&[
            (&grade_axis, 1.0),
            ("source-temperature-known", 0.0),
            ("source-moisture-condition-known", 0.0),
            ("source-test-rate-known", 0.0),
            ("full-sheet-at-least-10mm", 1.0),
        ]);
        sample(&pack, "density", RHO, &p, rho);
        sample(&pack, "compression-modulus", PA, &p, modulus);
        refuses(
            &pack,
            "density",
            RHO,
            &p.clone().with("full-sheet-at-least-10mm", 0.0).unwrap(),
        );
    }
    for (slug, grade, modulus) in [
        (
            "west-system-105-205-neat-cured",
            "source-grade-west-105-205",
            2813060975.612544,
        ),
        (
            "west-system-gflex-650-neat-cured",
            "source-grade-west-gflex-650",
            992845050.216192,
        ),
    ] {
        let pack = compile(slug, 7);
        let p = flags(&[
            (grade, 1.0),
            ("room-temperature-cure", 1.0),
            ("neat-unreinforced-resin", 1.0),
            ("source-temperature-known", 0.0),
            ("source-test-rate-known", 0.0),
        ])
        .with_quantity(
            "cure-age",
            QuantitySpec::dimensional(Dims([0, 0, 1, 0, 0, 0])),
            1209600.0,
        )
        .unwrap();
        sample(&pack, "tensile-modulus", PA, &p, modulus);
        refuses(
            &pack,
            "tensile-modulus",
            PA,
            &p.clone().with("cure-age", 604800.0).unwrap(),
        );
        refuses(&pack, "density", RHO, &p); // specific gravity alone has no declared water reference
    }
}

#[test]
fn g0_g3_insulation_curves_and_density_temperature_surface_keep_typed_axes() {
    let rock = compile("rockwool-prorox-sl960-en-2020", 3);
    let p = temperature(
        flags(&[
            ("source-grade-prorox-sl960-en-2020", 1.0),
            ("source-temperature-known", 1.0),
            ("source-moisture-condition-known", 0.0),
            ("source-pressure-known", 0.0),
        ]),
        "temperature",
        348.15,
    );
    sample(&rock, "thermal-conductivity", K, &p, 0.0425);
    refuses(
        &rock,
        "thermal-conductivity",
        K,
        &p.clone().with("temperature", 300.0).unwrap(),
    );
    refuses(
        &rock,
        "thermal-conductivity",
        K,
        &p.with_quantity("temperature", QuantitySpec::dimensional(TEMP), 348.15)
            .unwrap(),
    );

    let eps = compile("eps-nist-srm1453-2013", 1);
    let p = temperature(
        flags(&[
            ("source-lot-srm1453", 1.0),
            ("source-conditioning-protocol-satisfied", 1.0),
            ("unstacked-unsliced-uncontaminated-board", 1.0),
        ]),
        "mean-temperature",
        297.0,
    )
    .with_quantity("bulk-density", QuantitySpec::dimensional(RHO), 41.5)
    .unwrap()
    .with_quantity(
        "temperature-difference",
        QuantitySpec::dimensional(TEMP),
        20.0,
    )
    .unwrap();
    let answer = sample(
        &eps,
        "thermal-conductivity",
        K,
        &p,
        0.00111 - 0.0000424 * 41.5 + 0.000115 * 297.0,
    );
    assert!(matches!(answer.evidence.value.uncertainty,
        UncertaintyModel::RelativeHalfWidth { fraction, confidence } if fraction == 0.015 && confidence == 0.95));
    for bad in [
        p.clone().with("bulk-density", 36.0).unwrap(),
        p.clone().with("mean-temperature", 280.0).unwrap(),
        p.clone().with("temperature-difference", 10.0).unwrap(),
        p.with("source-lot-srm1453", 0.0).unwrap(),
    ] {
        refuses(&eps, "thermal-conductivity", K, &bad);
    }
    let repeated = compile("eps-nist-srm1453-2013", 1);
    assert_eq!(
        eps.to_bytes(),
        repeated.to_bytes(),
        "source surface compilation is deterministic"
    );
}

#[test]
fn g0_g3_xps_keeps_declared_conductivity_thickness_gaps_and_grade_limits() {
    let cases = [
        (
            "30",
            9,
            30.0,
            12e6,
            300e3,
            &[
                (0.03, 0.033),
                (0.04, 0.032),
                (0.05, 0.032),
                (0.06, 0.033),
                (0.07, 0.035),
                (0.12, 0.035),
                (0.14, 0.036),
                (0.16, 0.036),
            ] as &[(f64, f64)],
            [0.035, 0.13, 0.20],
        ),
        (
            "70",
            7,
            39.0,
            25e6,
            700e3,
            &[
                (0.05, 0.033),
                (0.06, 0.033),
                (0.08, 0.035),
                (0.14, 0.035),
                (0.16, 0.036),
            ] as &[(f64, f64)],
            [0.07, 0.15, 0.20],
        ),
    ];
    for (grade, count, rho, modulus, minimum_strength, thicknesses, gaps) in cases {
        let slug = format!("austrotherm-xps-top{grade}-sf-2026");
        let pack = compile(&slug, count);
        let grade_axis = format!("source-grade-{slug}");
        let p = flags(&[
            (&grade_axis, 1.0),
            ("source-temperature-known", 0.0),
            ("source-moisture-condition-known", 0.0),
            ("source-pressure-known", 0.0),
        ]);
        sample(&pack, "density", RHO, &p, rho);
        sample(&pack, "closed-cell-fraction-minimum", Dims::NONE, &p, 0.95);
        let mechanical = p
            .clone()
            .with_quantity(
                "source-test-rate-known",
                QuantitySpec::dimensional(Dims::NONE),
                0.0,
            )
            .unwrap();
        sample(&pack, "compression-modulus", PA, &mechanical, modulus);
        let strength = mechanical
            .with_quantity(
                "source-cs10y-class",
                QuantitySpec::dimensional(Dims::NONE),
                1.0,
            )
            .unwrap();
        sample(
            &pack,
            "compression-strength-minimum",
            PA,
            &strength,
            minimum_strength,
        );

        let declared = p
            .clone()
            .with_quantity(
                "source-is-declared-value",
                QuantitySpec::dimensional(Dims::NONE),
                1.0,
            )
            .unwrap();
        let at_thickness = |metres| {
            declared
                .clone()
                .with_quantity(
                    "nominal-thickness",
                    QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0])),
                    metres,
                )
                .unwrap()
        };
        for &(metres, expected) in thicknesses {
            sample(
                &pack,
                "thermal-conductivity",
                K,
                &at_thickness(metres),
                expected,
            );
        }
        for unsupported in gaps {
            refuses(&pack, "thermal-conductivity", K, &at_thickness(unsupported));
        }
        refuses(
            &pack,
            "thermal-conductivity",
            K,
            &at_thickness(0.05)
                .with("source-is-declared-value", 0.0)
                .unwrap(),
        );
        refuses(&pack, "density", RHO, &p.with(&grade_axis, 0.0).unwrap());
    }
}

#[test]
fn g0_g3_peek_model_curves_preserve_cubic_accuracy_and_explicit_model_domain() {
    let pack = compile("peek-nasa-thermic-plate-linearized-model", 3);
    let at = |kelvin| {
        temperature(
            flags(&[
                ("source-pressure-atmospheric", 1.0),
                ("source-constant-density-model", 1.0),
                ("source-test-temperature-known", 0.0),
            ]),
            "temperature",
            kelvin,
        )
    };
    let lower = at(300.0);
    let upper = at(525.0);
    let conductivity = key("thermal-conductivity", K);
    let heat_capacity = key("specific-heat-capacity", CP);
    let selection = fs_matdb::ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
    let k_line = pack
        .claims()
        .query_line_typed(&conductivity, "temperature", &lower, &upper, selection)
        .unwrap();
    let cp_line = pack
        .claims()
        .query_line_typed(&heat_capacity, "temperature", &lower, &upper, selection)
        .unwrap();
    assert_eq!(k_line.knots.len(), 46);
    assert_eq!(cp_line.knots.len(), 46);
    let cubic_k =
        |t: f64| -0.040607 + 0.0018749 * t - 0.0000044753 * t * t + 0.0000000039027 * t * t * t;
    let cubic_cp = |t: f64| 336.06 + 3.3066 * t - 0.0061332 * t * t + 0.000010477 * t * t * t;
    for index in 0..46 {
        let t = 300.0 + 5.0 * index as f64;
        assert_eq!(k_line.knots[index].0, t);
        assert_eq!(cp_line.knots[index].0, t);
        sample(&pack, "thermal-conductivity", K, &at(t), cubic_k(t));
        sample(&pack, "specific-heat-capacity", CP, &at(t), cubic_cp(t));
        sample(&pack, "density", RHO, &at(t), 1264.0);
        if index < 45 {
            let midpoint = t + 2.5;
            for (name, dims, polynomial, bound) in [
                (
                    "thermal-conductivity",
                    K,
                    cubic_k(midpoint),
                    0.000010446578125,
                ),
                (
                    "specific-heat-capacity",
                    CP,
                    cubic_cp(midpoint),
                    0.06480046875,
                ),
            ] {
                let answer = pack
                    .claims()
                    .query_typed(
                        &key(name, dims),
                        &at(midpoint),
                        SelectionPolicy::SingleClaimOnly,
                    )
                    .unwrap();
                assert!((answer.evidence.value.value - polynomial).abs() <= bound);
                assert!(matches!(
                    answer.evidence.value.uncertainty,
                    UncertaintyModel::Unstated,
                ));
                pack.claims().verify_receipt(&answer.receipt).unwrap();
            }
        }
    }
    let integrated_cp: f64 = cp_line
        .knots
        .windows(2)
        .map(|pair| {
            (pair[1].0 - pair[0].0)
                * 0.5
                * (pair[0].1.evidence.value.value + pair[1].1.evidence.value.value)
        })
        .sum();
    // Independent analytic integral of the reported cubic from 300 to 525 K.
    let cubic_integral = 319_641.290_332_031_25;
    assert!((integrated_cp - cubic_integral).abs() <= 14.580_105_468_75);
    assert!((integrated_cp - 319_647.695_414_062_5).abs() < 1e-7);
    for (name, dims) in [
        ("thermal-conductivity", K),
        ("specific-heat-capacity", CP),
        ("density", RHO),
    ] {
        refuses(&pack, name, dims, &at(299.99));
        refuses(&pack, name, dims, &at(525.01));
        refuses(
            &pack,
            name,
            dims,
            &lower
                .clone()
                .with("source-constant-density-model", 0.0)
                .unwrap(),
        );
        refuses(
            &pack,
            name,
            dims,
            &temperature(QueryPoint::new(), "temperature", 400.0),
        );
    }
}

#[path = "material_catalog_expansion/engineering.rs"]
mod engineering;

#[test]
fn g0_g3_airex_moduli_require_method_direction_and_preserve_minimums() {
    for (grade, iso_mpa, astm_mpa, iso_min_mpa, astm_min_mpa, strength_mpa, strength_min_mpa) in [
        ("55", 69.0, 58.0, 55.0, 45.0, 0.90, 0.75),
        ("130", 170.0, 145.0, 145.0, 123.0, 3.0, 2.6),
    ] {
        let slug = format!("airex-c70-{grade}-2022");
        let pack = compile(&slug, 18);
        let grade_axis = format!("source-grade-{slug}");
        let compression = flags(&[
            (&grade_axis, 1.0),
            ("source-temperature-known", 0.0),
            ("source-moisture-condition-known", 0.0),
            ("source-pressure-known", 0.0),
            ("source-test-rate-known", 0.0),
            ("block-center-specimen", 1.0),
            ("source-specimen-at-least-40mm", 1.0),
            ("loading-in-sheet-plane", 0.0),
            ("loading-normal-to-sheet", 1.0),
        ]);
        // The strength row lists ISO 844 / ASTM C365 together, but the
        // modulus rows report different results under those methods.
        sample(
            &pack,
            "compression-strength",
            PA,
            &compression,
            strength_mpa * 1e6,
        );
        sample(
            &pack,
            "compression-strength-minimum",
            PA,
            &compression,
            strength_min_mpa * 1e6,
        );
        for name in ["compression-modulus", "compression-modulus-minimum"] {
            refuses(&pack, name, PA, &compression);
        }
        for (method, average, minimum) in
            [(0.0, iso_mpa, iso_min_mpa), (1.0, astm_mpa, astm_min_mpa)]
        {
            let point = compression
                .clone()
                .with_quantity(
                    "compression-modulus-method-astm-c365",
                    QuantitySpec::dimensional(Dims::NONE),
                    method,
                )
                .unwrap();
            // DNV-GL-defined minimums remain distinct source properties,
            // not uncertainty bounds on the approximate nominal averages.
            sample(&pack, "compression-modulus", PA, &point, average * 1e6);
            sample(
                &pack,
                "compression-modulus-minimum",
                PA,
                &point,
                minimum * 1e6,
            );
            let in_plane = point
                .clone()
                .with("loading-in-sheet-plane", 1.0)
                .unwrap()
                .with("loading-normal-to-sheet", 0.0)
                .unwrap();
            let unknown_method = point
                .with("compression-modulus-method-astm-c365", 2.0)
                .unwrap();
            for name in ["compression-modulus", "compression-modulus-minimum"] {
                refuses(&pack, name, PA, &in_plane);
                refuses(&pack, name, PA, &unknown_method);
            }
        }
    }
}

#[test]
fn g0_g3_ppsu_compression_points_do_not_create_an_interpolated_curve() {
    let pack = compile("ensinger-tecason-p-black-radel-5500", 10);
    let base = temperature(
        flags(&[
            ("source-grade-tecason-p-black-radel-5500", 1.0),
            ("source-moisture-condition-known", 0.0),
            ("source-pressure-known", 0.0),
            ("source-measurement-orientation-known", 0.0),
            ("source-test-rate-known", 0.0),
        ]),
        "temperature",
        295.92777777777775, // 73 F; the mechanical table's test condition
    );
    let at = |strain| {
        base.clone()
            .with_quantity(
                "compressive-strain",
                QuantitySpec::dimensional(Dims::NONE),
                strain,
            )
            .unwrap()
    };
    // Independently convert the source's 1,990 and 14,350 psi values.
    // The 10% row retains its producer-secondary observation in the pack.
    let psi_to_pa = 0.45359237 * 9.80665 / (0.0254 * 0.0254);
    sample(
        &pack,
        "compressive-stress",
        PA,
        &at(0.01),
        1990.0 * psi_to_pa,
    );
    sample(
        &pack,
        "compressive-stress",
        PA,
        &at(0.10),
        14350.0 * psi_to_pa,
    );
    refuses(&pack, "compressive-stress", PA, &base);
    for unsupported_strain in [0.0, 0.02, 0.05, 0.11] {
        refuses(&pack, "compressive-stress", PA, &at(unsupported_strain));
    }
    refuses(
        &pack,
        "compressive-stress",
        PA,
        &at(0.01).with("temperature", 296.15).unwrap(),
    );
    assert!(
        pack.claims()
            .query_line_typed(
                &key("compressive-stress", PA),
                "compressive-strain",
                &at(0.01),
                &at(0.10),
                fs_matdb::ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
            )
            .is_err(),
        "two isolated test points must not become a compression constitutive curve"
    );
}

#[path = "material_catalog_expansion/alloys_minerals.rs"]
mod alloys_minerals;

#[path = "material_catalog_expansion/polymers_foams.rs"]
mod polymers_foams;

#[path = "material_catalog_expansion/adhesives.rs"]
mod adhesives;

#[path = "material_catalog_expansion/cushioning.rs"]
mod cushioning;

#[path = "material_catalog_expansion/thermal_engineering.rs"]
mod thermal_engineering;

#[path = "material_catalog_expansion/compression_reference.rs"]
mod compression_reference;

#[path = "material_catalog_expansion/alloy_completion.rs"]
mod alloy_completion;

#[path = "material_catalog_expansion/printed_elastomers.rs"]
mod printed_elastomers;

#[path = "material_catalog_expansion/ceramics_elements.rs"]
mod ceramics_elements;
