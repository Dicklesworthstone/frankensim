//! Source-to-geometry and density checks for FPL's bounded moisture model.
//! These exercise the published engineering approximation, not experimental
//! validation, a moisture transport solver or a longitudinal shrinkage law.
#![deny(unsafe_code)]

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use fs_matdb::{
    MaterialStateId, NormalizedMaterialCardPack, NormalizedPack, PropertyKey, QueryPoint,
    UncertaintyModel,
};
use fs_material::state_point::MaterialPropertySelection;
use fs_material::wood_moisture::{ResolvedWoodMoisture, WoodMoistureSource};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec};

const RHO: Dims = Dims([-3, 1, 0, 0, 0, 0]);

struct Species {
    slug: &'static str,
    claims: usize,
    specific_gravity: f64,
    shrinkage_percent: [f64; 3],
}

const fn species(
    slug: &'static str,
    claims: usize,
    specific_gravity: f64,
    shrinkage_percent: [f64; 3],
) -> Species {
    Species {
        slug,
        claims,
        specific_gravity,
        shrinkage_percent,
    }
}

const SPECIES: [Species; 33] = [
    species("pine-loblolly", 25, 0.51, [4.8, 7.4, 12.3]),
    species("pine-longleaf", 25, 0.59, [5.1, 7.5, 12.2]),
    species("pine-slash", 23, 0.59, [5.4, 7.6, 12.1]),
    species("pine-ponderosa", 25, 0.40, [3.9, 6.2, 9.7]),
    species("oak-northern-red", 16, 0.63, [4.0, 8.6, 13.7]),
    species("oak-white", 16, 0.68, [5.6, 10.5, 16.3]),
    species("cedar-northern-white", 25, 0.31, [2.2, 4.9, 7.2]),
    species("fir-subalpine", 24, 0.32, [2.6, 7.4, 9.4]),
    species("douglas-fir-coast", 21, 0.48, [4.8, 7.6, 12.4]),
    species("ash-white", 20, 0.60, [4.9, 7.8, 13.3]),
    species("baldcypress", 19, 0.46, [3.8, 6.2, 10.5]),
    species("basswood-american", 20, 0.37, [6.6, 9.3, 15.8]),
    species("birch-yellow", 21, 0.62, [7.3, 9.5, 16.8]),
    species("cedar-western-red", 19, 0.32, [2.4, 5.0, 6.8]),
    species("cherry-black", 20, 0.50, [3.7, 7.1, 11.5]),
    species("cottonwood-eastern", 20, 0.40, [3.9, 9.2, 13.9]),
    species("hemlock-western", 19, 0.45, [4.2, 7.8, 12.4]),
    species("larch-western", 19, 0.52, [4.5, 9.1, 14.0]),
    species("maple-red", 20, 0.54, [4.0, 8.2, 12.6]),
    species("maple-sugar", 22, 0.63, [4.8, 9.9, 14.7]),
    species("pine-sugar", 21, 0.36, [2.9, 5.6, 7.9]),
    species("redwood-old-growth", 19, 0.40, [2.6, 4.4, 6.8]),
    species("spruce-engelmann", 23, 0.35, [3.8, 7.1, 11.0]),
    species("spruce-sitka", 23, 0.40, [4.3, 7.5, 11.5]),
    species("stick-hickory-shagbark", 24, 0.72, [7.0, 10.5, 16.7]),
    species("sweetgum", 21, 0.52, [5.3, 10.2, 15.8]),
    species("walnut-black", 21, 0.55, [5.5, 7.8, 12.8]),
    species("yellow-poplar", 21, 0.42, [4.6, 8.2, 12.7]),
    // These five specific gravities use GREEN volume. Their source packs
    // explicitly derive the 12%-moisture density with the volumetric model.
    species("balsa", 20, 0.16, [3.0, 7.6, 10.8]),
    species("mahogany-african", 20, 0.42, [2.5, 4.5, 8.8]),
    species("mahogany-honduras", 20, 0.45, [3.0, 4.1, 7.8]),
    species("rosewood-brazilian", 9, 0.80, [2.9, 4.6, 8.5]),
    species("rosewood-indian", 9, 0.75, [2.7, 5.8, 8.5]),
];

impl Species {
    fn moisture_axis(&self) -> &'static str {
        if matches!(
            self.slug,
            "pine-loblolly"
                | "pine-longleaf"
                | "pine-slash"
                | "pine-ponderosa"
                | "oak-northern-red"
                | "oak-white"
                | "cedar-northern-white"
                | "fir-subalpine"
        ) {
            "moisture-content"
        } else {
            "moisture_content"
        }
    }

    fn reference_density(&self) -> f64 {
        let reference_volume_over_sg_volume = if matches!(
            self.slug,
            "balsa"
                | "mahogany-african"
                | "mahogany-honduras"
                | "rosewood-brazilian"
                | "rosewood-indian"
        ) {
            1.0 - self.shrinkage_percent[2] / 100.0 * (1.0 - 0.12 / 0.30)
        } else {
            1.0
        };
        self.specific_gravity * 1120.0 / reference_volume_over_sg_volume
    }
}

fn workspace() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .to_path_buf()
}

fn compile(species: &Species) -> NormalizedMaterialCardPack {
    static NEXT: AtomicU64 = AtomicU64::new(0);
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "frankensim-wood-moisture-{}-{nonce}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir(&dir).unwrap();
    let slug = format!("{}-fpl-gtr282", species.slug);
    let output = dir.join("material.fsmatpk");
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(workspace().join(format!("data/matdb/seed-v1/{slug}/manifest.tsv")))
        .arg("--out")
        .arg(&output)
        .env("CARGO_WORKSPACE_DIR", workspace())
        .current_dir(workspace())
        .output()
        .expect("run the actual source compiler");
    assert!(
        run.status.success(),
        "{slug}: compiler status {}; {}",
        run.status,
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(output).unwrap();
    let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
    assert_eq!(pack.pack_id(), slug);
    assert_eq!(pack.claims().claims_ordered().count(), species.claims);
    NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: slug,
            phase: "clear wood; moisture on dry-mass basis".into(),
            process: "FPL green-to-ovendry shrinkage with explicit fiber-saturation approximation"
                .into(),
            revision: 0,
        },
        pack,
    )
    .unwrap()
}

fn keys() -> [PropertyKey; 5] {
    [
        PropertyKey::new("density", RHO),
        PropertyKey::new("radial_shrinkage_green_to_ovendry", Dims::NONE),
        PropertyKey::new("tangential_shrinkage_green_to_ovendry", Dims::NONE),
        PropertyKey::new("volumetric_shrinkage_green_to_ovendry", Dims::NONE),
        PropertyKey::new("fiber_saturation_moisture_content", Dims::NONE),
    ]
}

fn source<'a>(
    keys: &'a [PropertyKey; 5],
    axis: &'a str,
    reference: &'a QueryPoint,
    shrinkage: &'a QueryPoint,
) -> WoodMoistureSource<'a> {
    WoodMoistureSource {
        density: &keys[0],
        radial_shrinkage: &keys[1],
        tangential_shrinkage: &keys[2],
        volumetric_shrinkage: &keys[3],
        fiber_saturation: &keys[4],
        moisture_axis: axis,
        reference_point: reference,
        shrinkage_point: shrinkage,
    }
}

fn point(axis: &str, moisture: f64) -> QueryPoint {
    let legacy = axis == "moisture_content";
    [
        (axis, moisture),
        (
            if legacy {
                "clear_specimen"
            } else {
                "clear-specimen"
            },
            1.0,
        ),
        (
            if legacy {
                "engineering_reference_model"
            } else {
                "engineering-reference-model"
            },
            1.0,
        ),
        (
            if legacy {
                "source_test_temperature_known"
            } else {
                "source-test-temperature-known"
            },
            0.0,
        ),
    ]
    .into_iter()
    .fold(QueryPoint::new(), |p, (name, value)| {
        p.with_quantity(name, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap()
    })
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-12 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn g0_g1_all_33_source_species_drive_anisotropic_geometry_and_conserve_dry_mass() {
    let keys = keys();
    for species in &SPECIES {
        let normalized = compile(species);
        let card = normalized.card();
        let axis = species.moisture_axis();
        let mut reference = point(axis, 0.12);
        if species.slug == "stick-hickory-shagbark" {
            reference = reference
                .with_quantity(
                    "T",
                    QuantitySpec::semantic(SemanticType::new(
                        QuantityKind::AbsoluteTemperature,
                        ValueForm::Static,
                    )),
                    293.15,
                )
                .unwrap();
        }
        let law = ResolvedWoodMoisture::from_card(
            card,
            source(&keys, axis, &reference, &reference),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap_or_else(|e| panic!("{}: {e:?}", species.slug));
        close(law.fiber_saturation(), 0.30);
        assert_eq!(law.moisture_span(), [0.0, 0.30]);
        let expected_source = [
            species.reference_density(),
            species.shrinkage_percent[0] / 100.0,
            species.shrinkage_percent[1] / 100.0,
            species.shrinkage_percent[2] / 100.0,
            0.30,
        ];
        for (i, key) in keys.iter().enumerate() {
            let owner = if i == 0 {
                law.density_source()
            } else {
                law.shrinkage_source()
            };
            assert_eq!(owner.card_identity(), card.content_hash());
            let property = owner.property_by_key(key).unwrap();
            close(property.value_si(), expected_source[i]);
            assert_eq!(
                property.answer().evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
            let receipt = &property.answer().receipt;
            assert!(receipt.observation_backed);
            assert!(!receipt.source_hashes.is_empty());
            card.claims().verify_receipt(receipt).unwrap();
        }

        let at_reference = law.at_moisture_content(0.12).unwrap();
        close(at_reference.density_kg_m3, expected_source[0]);
        assert_eq!(at_reference.radial_length_factor, 1.0);
        assert_eq!(at_reference.tangential_length_factor, 1.0);
        assert_eq!(at_reference.volume_factor, 1.0);
        let dry_mass_kg = expected_source[0] / 1.12 * 0.002;
        for m in [0.0, 0.08, 0.20, 0.30] {
            let state = law.at_moisture_content(m).unwrap();
            assert_eq!(state.material_card_identity, card.content_hash());
            assert_eq!(state.moisture_content, m);
            // FPL Eq.4-9 independently applied to each tabulated R/T/V value.
            // No longitudinal shrinkage is inferred and V is never R times T.
            let factors = species.shrinkage_percent.map(|percent| {
                let shrinkage = percent / 100.0;
                (1.0 - shrinkage * (1.0 - m / 0.30)) / (1.0 - shrinkage * (1.0 - 0.12 / 0.30))
            });
            close(state.radial_length_factor, factors[0]);
            close(state.tangential_length_factor, factors[1]);
            close(state.volume_factor, factors[2]);
            close(
                state.density_kg_m3,
                expected_source[0] / 1.12 * (1.0 + m) / factors[2],
            );
            let volume_m3 = 0.002 * state.volume_factor;
            let mass_kg = state.density_kg_m3 * volume_m3;
            close(mass_kg / (1.0 + m), dry_mass_kg);
            close(mass_kg - dry_mass_kg, m * dry_mass_kg);
            assert!(
                (state.volume_factor - state.radial_length_factor * state.tangential_length_factor)
                    .abs()
                    > 1e-5,
                "{}: independent volume shrinkage was replaced",
                species.slug
            );
        }
        let dry = law.at_moisture_content(0.08).unwrap();
        let wet = law.at_moisture_content(0.20).unwrap();
        assert!(dry.tangential_length_factor < dry.radial_length_factor);
        assert!(wet.tangential_length_factor > wet.radial_length_factor);
        assert!(dry.volume_factor < 1.0 && wet.volume_factor > 1.0);
        assert!(dry.density_kg_m3 < at_reference.density_kg_m3);
        assert!(wet.density_kg_m3 > at_reference.density_kg_m3);
        eprintln!(
            "{}: rho8={:.6}, rho12={:.6}, rho20={:.6} kg/m3; R20={:.9}, T20={:.9}, V20={:.9}",
            species.slug,
            dry.density_kg_m3,
            at_reference.density_kg_m3,
            wet.density_kg_m3,
            wet.radial_length_factor,
            wet.tangential_length_factor,
            wet.volume_factor,
        );

        for invalid_moisture in [-0.001, 0.300001, f64::NAN, f64::INFINITY] {
            assert!(law.at_moisture_content(invalid_moisture).is_err());
        }
        let wrong_reference = reference.clone().with(axis, 0.08).unwrap();
        assert!(
            ResolvedWoodMoisture::from_card(
                card,
                source(&keys, axis, &wrong_reference, &reference),
                MaterialPropertySelection::SingleClaimOnly,
            )
            .is_err()
        );
        // Turning off the named engineering approximation cannot silently
        // turn the model into a measured material response.
        let approximation_axis = if axis == "moisture_content" {
            "engineering_reference_model"
        } else {
            "engineering-reference-model"
        };
        let wrong_context = reference.clone().with(approximation_axis, 0.0).unwrap();
        assert!(
            ResolvedWoodMoisture::from_card(
                card,
                source(&keys, axis, &reference, &wrong_context),
                MaterialPropertySelection::SingleClaimOnly,
            )
            .is_err()
        );
    }
}
