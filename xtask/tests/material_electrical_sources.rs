//! Actual catalog sources drive the existing Ohmic law with their exact keys.
//! These checks establish source transport and calculation, not experimental validation.
#![deny(unsafe_code)]

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use fs_matdb::{
    MaterialStateId, NormalizedMaterialCardPack, NormalizedPack, PropertyKey, QueryPoint,
};
use fs_material::conductor::{
    ELECTRICAL_RESISTIVITY_DIMS, ElectricalPropertyBasis, OhmicDrive, ResolvedConductorState,
    planar::PlanarConductorMaterial, resolve_conductor_state,
    resolve_conductor_state_with_property, resolve_uniform_conductor_with_property,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

const STEEL: &str = "stainless-14401-thyssenkrupp-2017";
const STEEL_FLAGS: &[(&str, f64)] = &[
    ("source-heat-treatment-known", 0.0),
    ("source-pressure-known", 0.0),
];
const MODEL_FLAGS: &[(&str, f64)] = &[
    ("engineering-reference-model", 1.0),
    ("source-temperatures-matched", 0.0),
];

fn workspace() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .to_path_buf()
}

fn compile(slug: &str, count: usize) -> NormalizedMaterialCardPack {
    static NEXT: AtomicU64 = AtomicU64::new(0);
    let time = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "frankensim-electrical-source-{}-{time}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed),
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
    NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: slug.into(),
            phase: "solid".into(),
            process: "source-qualified reference-state electrical calculation".into(),
            revision: 0,
        },
        pack,
    )
    .unwrap()
}

fn temperature(point: QueryPoint, kelvin: f64) -> QueryPoint {
    point
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            kelvin,
        )
        .unwrap()
}

fn point(kelvin: f64, flags: &[(&str, f64)]) -> QueryPoint {
    temperature(
        flags
            .iter()
            .fold(QueryPoint::new(), |point, (name, value)| {
                point
                    .with_quantity(*name, QuantitySpec::dimensional(Dims::NONE), *value)
                    .unwrap()
            }),
        kelvin,
    )
}

fn electrical_key() -> PropertyKey {
    PropertyKey::new("electrical-resistivity", ELECTRICAL_RESISTIVITY_DIMS)
}

fn conductivity_key() -> PropertyKey {
    PropertyKey::with_quantity(
        "thermal-conductivity",
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::ThermalConductivity,
            ValueForm::Static,
        )),
    )
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-12 * expected.abs(),
        "{actual} != {expected}"
    );
}

#[test]
fn six_catalog_resistivities_drive_current_heat_and_exact_source_refresh() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let key = electrical_key();
    for (slug, count, kelvin, rho, resistance) in [
        (STEEL, 10, 293.15, 7.5e-7, 1.5),
        (
            "aluminum-5005-austral-wright-reference",
            7,
            293.15,
            3.4e-8,
            0.068,
        ),
        (
            "aluminum-5052-austral-wright-reference",
            7,
            293.15,
            5.0e-8,
            0.1,
        ),
        ("magnesium-espi-reference", 4, 298.15, 4.45e-8, 0.089),
        ("molybdenum-attl-reference", 5, 293.15, 5.3e-8, 0.106),
        ("tantalum-efinea-reference", 4, 293.15, 1.25e-7, 0.25),
    ] {
        let normalized = compile(slug, count);
        let card = normalized.card();
        let flags = if slug == STEEL {
            STEEL_FLAGS
        } else {
            MODEL_FLAGS
        };
        let query = point(kelvin, flags);
        let bulk = resolve_conductor_state_with_property(
            card,
            &query,
            MaterialPropertySelection::SingleClaimOnly,
            basis,
            &key,
        )
        .unwrap_or_else(|error| panic!("{slug}: {error:?}"));
        assert_eq!(bulk.source_key(), &key);
        assert_eq!(bulk.source_basis(), basis);
        assert_eq!(bulk.material().query(), query);
        assert_eq!(bulk.material().card_identity(), card.content_hash());
        assert_eq!(
            bulk.source_property(),
            bulk.material().property_by_key(&key).unwrap()
        );
        close(bulk.resistivity_ohm_m(), rho);
        let receipt = bulk.source_property().answer().receipt.clone();
        assert!(receipt.observation_backed);
        assert!(!receipt.source_hashes.is_empty());
        card.claims().verify_receipt(&receipt).unwrap();
        assert_eq!(
            bulk.source_property_names().collect::<Vec<_>>(),
            [key.name()]
        );
        assert_eq!(
            bulk.source_properties().collect::<Vec<_>>(),
            [bulk.source_property()]
        );
        // Applicability-only source flags qualify this point calculation;
        // they do not satisfy the separate automatic coupled-source policy.
        assert!(!bulk.has_admitted_static_sources());
        close(bulk.conductivity().matrix_s_m()[0][0], rho.recip());

        let wire = bulk.clone().into_uniform_conductor(2.0, 1e-6).unwrap();
        assert_eq!(wire.source_key(), &key);
        assert_eq!(wire.material(), bulk.material());
        assert_eq!(wire.source_property().answer().receipt, receipt);
        close(wire.resistance_ohm(), resistance);
        close(wire.joule_power_w(3.0).unwrap(), 9.0 * resistance);
        close(wire.joule_power_w(-3.0).unwrap(), 9.0 * resistance);
        close(
            wire.current_a_for_drive(OhmicDrive::Voltage(3.0 * resistance))
                .unwrap(),
            3.0,
        );
        assert_eq!(
            wire,
            resolve_uniform_conductor_with_property(
                card,
                &query,
                MaterialPropertySelection::SingleClaimOnly,
                basis,
                &key,
                2.0,
                1e-6,
            )
            .unwrap()
        );

        // Pin the original source name; no underscore alias or reciprocal claim is inserted.
        let pinned = resolve_conductor_state_with_property(
            card,
            &query,
            MaterialPropertySelection::PinnedByProperty(vec![(
                key.name().into(),
                receipt.selected,
            )]),
            basis,
            &key,
        )
        .unwrap();
        assert_eq!(
            pinned.source_property().answer().receipt.selected,
            receipt.selected
        );
        assert!(
            resolve_conductor_state(card, &query, MaterialPropertySelection::SingleClaimOnly)
                .is_err()
        );
        assert!(
            resolve_conductor_state_with_property(
                card,
                &query,
                MaterialPropertySelection::PinnedByProperty(vec![(
                    "electrical_resistivity".into(),
                    receipt.selected
                )]),
                basis,
                &key,
            )
            .is_err()
        );
        assert!(
            resolve_conductor_state_with_property(
                card,
                &query,
                MaterialPropertySelection::SingleClaimOnly,
                ElectricalPropertyBasis::Conductivity,
                &key,
            )
            .is_err()
        );
        assert!(
            resolve_conductor_state_with_property(
                card,
                &query,
                MaterialPropertySelection::SingleClaimOnly,
                basis,
                &conductivity_key(),
            )
            .is_err()
        );

        for invalid in [
            point(kelvin, &[]),
            query.clone().with(flags[0].0, 1.0 - flags[0].1).unwrap(),
            temperature(query.clone(), kelvin - 1.0),
            temperature(query.clone(), kelvin + 2.0),
        ] {
            assert!(
                resolve_conductor_state_with_property(
                    card,
                    &invalid,
                    MaterialPropertySelection::SingleClaimOnly,
                    basis,
                    &key,
                )
                .is_err(),
                "{slug}: unsupported query admitted"
            );
        }

        // Only the authored reference models admit a changed temperature. The
        // steel source remains exactly at 20 C and cannot acquire a wider law.
        let refresh_query = temperature(
            query.clone(),
            kelvin + if slug == STEEL { 0.0 } else { 0.5 },
        );
        let refreshed = bulk
            .material()
            .at_query_point(card, &refresh_query)
            .unwrap();
        let planar = PlanarConductorMaterial::with_material(&bulk, refreshed.clone()).unwrap();
        let refreshed_wire = wire.with_material(refreshed.clone()).unwrap();
        assert_eq!(planar.source_key(), &key);
        assert_eq!(refreshed_wire.source_key(), &key);
        assert!(planar.has_same_electrical_source(&bulk));
        assert_eq!(planar.material(), &refreshed);
        assert_eq!(refreshed_wire.material(), &refreshed);
        assert_eq!(
            planar.source_property().answer().receipt.selected,
            receipt.selected
        );
        assert_eq!(
            refreshed_wire.source_property().answer().receipt.selected,
            receipt.selected
        );
        assert_eq!(
            planar.source_property().answer().receipt.source_hashes,
            receipt.source_hashes
        );
        assert_eq!(refreshed_wire.length_m(), wire.length_m());
        assert_eq!(refreshed_wire.area_m2(), wire.area_m2());
        close(refreshed_wire.resistance_ohm(), resistance);
        assert!(
            bulk.material()
                .at_query_point(card, &temperature(query, kelvin + 2.0))
                .is_err()
        );
    }
}

#[test]
fn stainless_conductor_retains_the_original_four_property_thermal_bundle() {
    let normalized = compile(STEEL, 10);
    let card = normalized.card();
    let query = point(293.15, STEEL_FLAGS);
    let electrical = electrical_key();
    let keys = [
        (
            PropertyKey::new("density", Dims([-3, 1, 0, 0, 0, 0])),
            8000.0,
        ),
        (
            PropertyKey::new("specific-heat-capacity", Dims([2, 0, -2, -1, 0, 0])),
            500.0,
        ),
        (conductivity_key(), 15.0),
        (electrical.clone(), 7.5e-7),
    ];
    let requirements = keys
        .iter()
        .map(|(key, _)| {
            ScalarPropertyRequirement::try_with_key(key, ScalarAdmissibility::StrictlyPositive)
                .unwrap()
        })
        .collect::<Vec<_>>();
    let bundle = resolve_material_state_point(
        card,
        &query,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let bulk = ResolvedConductorState::try_from_material_with_property(
        bundle.clone(),
        ElectricalPropertyBasis::Resistivity,
        &electrical,
    )
    .unwrap();
    let wire = bulk.clone().into_uniform_conductor(2.0, 1e-6).unwrap();
    assert_eq!(bulk.material(), &bundle);
    assert_eq!(wire.material(), &bundle);
    assert_eq!(wire.material().properties().len(), 4);
    for (key, expected) in &keys {
        let original = bundle.property_by_key(key).unwrap();
        close(original.value_si(), *expected);
        assert_eq!(wire.material().property_by_key(key), Some(original));
        assert_eq!(bulk.material().property_by_key(key), Some(original));
        card.claims()
            .verify_receipt(&original.answer().receipt)
            .unwrap();
    }
    close(wire.resistance_ohm(), 1.5);
    close(wire.joule_power_w(3.0).unwrap(), 13.5);
    assert!(
        bundle
            .at_query_point(card, &temperature(query, 294.15))
            .is_err()
    );
    assert_eq!(wire.material(), &bundle); // a refused refresh never mutates accepted source state
}
