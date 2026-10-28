//! Real compiler -> durable store/reopen -> typed discovery -> layered heat.
//! CFAST entries are nominal model inputs, not measured specimens or allowables.
#![deny(unsafe_code)]

use fs_conduction::layered::{LayeredWall, PlanarThermalLayer};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSelection, MaterialStateId, NormalizedMaterialCardPack, NormalizedPack, PropertyKey,
    QueryPoint, SelectionPolicy,
};
use fs_matdb_store::{
    CatalogPack, DiscoveryDomain, DiscoveryRequest, DiscoveryStatus, DiscoveryTarget, MaterialStore,
};
use fs_material::graph::LawRegistry;
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::path::{Path, PathBuf};
use std::process::Command;

fn scratch() -> PathBuf {
    let time = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let path = std::env::temp_dir().join(format!(
        "frankensim-catalog-wall-{}-{time}",
        std::process::id()
    ));
    std::fs::create_dir(&path).unwrap();
    path
}

fn compile(dir: &Path, slug: &str) -> NormalizedMaterialCardPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let output = dir.join(format!("{slug}.fsmatpk"));
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(root.join(format!("data/matdb/seed-v1/{slug}/manifest.tsv")))
        .arg("--out")
        .arg(&output)
        .current_dir(root)
        .env("CARGO_WORKSPACE_DIR", root)
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{slug}: {}",
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(output).unwrap();
    let pack = NormalizedPack::from_bytes(&bytes).unwrap();
    let verified = NormalizedPack::from_bytes_verified(pack.content_hash(), &bytes).unwrap();
    assert_eq!(verified, pack);
    NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: slug.into(),
            phase: "solid".into(),
            process: "CFAST nominal layer-input reference; source test conditions unspecified"
                .into(),
            revision: 0,
        },
        pack,
    )
    .unwrap()
}

fn point(kelvin: f64, condition: &str, moisture_known: f64) -> QueryPoint {
    let mut point = QueryPoint::new();
    for (axis, value) in [
        ("source-test-temperature-known", 0.0),
        ("source-moisture-condition-known", moisture_known),
        ("source-measurement-method-known", 0.0),
        (condition, 0.0),
    ] {
        point = point
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
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

fn with_cx<T>(run: impl FnOnce(&Cx<'_>) -> T) -> T {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        run(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0,
                kernel_id: 51,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 1e-10 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn g1_g3_catalog_wall_survives_store_discovery_and_material_substitution() {
    let directory = scratch();
    let specifications = [
        (
            "gypsum-board-5-8-cfast",
            0.016,
            0.16,
            "source-product-and-formulation-known",
        ),
        (
            "glass-fiber-insulation-cfast",
            0.088,
            0.04,
            "source-product-and-install-density-known",
        ),
        (
            "brick-common-cfast-sp1041",
            0.076,
            0.72,
            "source-material-processing-known",
        ),
    ];
    let original = specifications.map(|(slug, _, _, _)| compile(&directory, slug));
    let database = directory.join("wall.sqlite");
    {
        let store = MaterialStore::open(database.to_str().unwrap()).unwrap();
        let bundle = original
            .iter()
            .cloned()
            .map(CatalogPack::MaterialCard)
            .collect::<Vec<_>>();
        store.ingest_bundle(&bundle).unwrap();
        store.seal_corpus().unwrap();
    }
    let store = MaterialStore::open(database.to_str().unwrap()).unwrap();
    let key = PropertyKey::with_quantity(
        "thermal-conductivity",
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::ThermalConductivity,
            ValueForm::Static,
        )),
    );
    let mut layers = Vec::new();
    for (original, (slug, thickness, k, condition)) in original.iter().zip(specifications) {
        let points = [point(273.15, condition, 0.0), point(293.15, condition, 0.0)];
        let request = DiscoveryRequest {
            target: DiscoveryTarget::Materials,
            properties: vec![key.clone()],
            models: Vec::new(),
            domain: DiscoveryDomain::Envelope {
                lower: points[0].clone(),
                upper: points[1].clone(),
            },
            selection: SelectionPolicy::SingleClaimOnly,
        };
        let discovered = store.discover(&request, &LawRegistry::new()).unwrap();
        let candidate = discovered
            .candidates
            .iter()
            .find(|c| c.pack.pack_id == original.pack_id())
            .unwrap();
        assert_eq!(candidate.status, DiscoveryStatus::Complete);
        assert_eq!(candidate.pack.content_hash, original.content_hash());
        let CatalogPack::MaterialCard(loaded) =
            store.load_catalog_pack(original.pack_id()).unwrap()
        else {
            panic!("stored card changed family")
        };
        assert_eq!(&loaded, original);
        let layer = PlanarThermalLayer::from_card(
            slug,
            thickness,
            loaded.card(),
            &key,
            "temperature",
            &points,
            ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        )
        .unwrap();
        // Source values are independently transcribed literals from the CFAST
        // entries. No density or heat capacity is needed for this steady solve.
        for p in &points {
            let answer = loaded
                .card()
                .claims()
                .query_typed(&key, p, SelectionPolicy::SingleClaimOnly)
                .unwrap();
            close(answer.evidence.value.value, k);
            loaded
                .card()
                .claims()
                .verify_receipt(&answer.receipt)
                .unwrap();
        }
        let incompatible = [point(273.15, condition, 1.0), point(293.15, condition, 1.0)];
        assert!(
            PlanarThermalLayer::from_card(
                slug,
                thickness,
                loaded.card(),
                &key,
                "temperature",
                &incompatible,
                ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
            )
            .is_err(),
            "unknown source moisture must not become known"
        );
        layers.push(layer);
    }
    let wall = LayeredWall::new(layers.clone(), vec![0.01, 0.002], [0.13, 0.04]).unwrap();
    let solved = with_cx(|cx| wall.solve(cx, [293.15, 273.15], 10.0)).unwrap();
    let expected_resistance =
        0.13 + 0.016 / 0.16 + 0.01 + 0.088 / 0.04 + 0.002 + 0.076 / 0.72 + 0.04;
    let expected_flux = 20.0 / expected_resistance;
    // Report field names follow the production LayeredWallSolution API.
    close(solved.heat_flux_w_per_m2, expected_flux);
    close(solved.heat_rate_w, 10.0 * expected_flux);
    let again = with_cx(|cx| wall.solve(cx, [293.15, 273.15], 10.0)).unwrap();
    assert_eq!(solved, again, "same-ISA deterministic replay");

    // Halve only the insulation geometry; retain the exact same source table.
    let thinner = PlanarThermalLayer::new(
        "half-thickness-insulation",
        0.044,
        layers[1].conductivity().clone(),
    )
    .unwrap();
    layers[1] = thinner;
    let changed = LayeredWall::new(layers, vec![0.01, 0.002], [0.13, 0.04]).unwrap();
    let changed = with_cx(|cx| changed.solve(cx, [293.15, 273.15], 10.0)).unwrap();
    close(
        changed.heat_flux_w_per_m2,
        20.0 / (expected_resistance - 1.1),
    );
    assert!(changed.heat_flux_w_per_m2 > solved.heat_flux_w_per_m2);
    println!(
        "CFAST wall: R''={expected_resistance} m2 K/W; q''={expected_flux} W/m2; \
         heat={} W; database={}",
        solved.heat_rate_w,
        database.display()
    );
}
