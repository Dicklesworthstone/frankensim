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
use fs_qty::{Dims, QuantitySpec};

const RHO: Dims = Dims([-3, 1, 0, 0, 0, 0]);

struct Species {
    slug: &'static str,
    claims: usize,
    specific_gravity: f64,
    shrinkage_percent: [f64; 3],
}

const SPECIES: [Species; 9] = [
    Species {
        slug: "pine-loblolly",
        claims: 25,
        specific_gravity: 0.51,
        shrinkage_percent: [4.8, 7.4, 12.3],
    },
    Species {
        slug: "pine-longleaf",
        claims: 25,
        specific_gravity: 0.59,
        shrinkage_percent: [5.1, 7.5, 12.2],
    },
    Species {
        slug: "pine-slash",
        claims: 23,
        specific_gravity: 0.59,
        shrinkage_percent: [5.4, 7.6, 12.1],
    },
    Species {
        slug: "pine-ponderosa",
        claims: 25,
        specific_gravity: 0.40,
        shrinkage_percent: [3.9, 6.2, 9.7],
    },
    Species {
        slug: "oak-northern-red",
        claims: 16,
        specific_gravity: 0.63,
        shrinkage_percent: [4.0, 8.6, 13.7],
    },
    Species {
        slug: "oak-white",
        claims: 16,
        specific_gravity: 0.68,
        shrinkage_percent: [5.6, 10.5, 16.3],
    },
    Species {
        slug: "cedar-northern-white",
        claims: 25,
        specific_gravity: 0.31,
        shrinkage_percent: [2.2, 4.9, 7.2],
    },
    Species {
        slug: "fir-subalpine",
        claims: 24,
        specific_gravity: 0.32,
        shrinkage_percent: [2.6, 7.4, 9.4],
    },
    Species {
        slug: "douglas-fir-coast",
        claims: 21,
        specific_gravity: 0.48,
        shrinkage_percent: [4.8, 7.6, 12.4],
    },
];

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
        "{slug}: {}",
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
fn g0_g1_nine_source_species_drive_anisotropic_geometry_and_conserve_dry_mass() {
    let keys = keys();
    for species in &SPECIES {
        let normalized = compile(species);
        let card = normalized.card();
        let axis = if species.slug == "douglas-fir-coast" {
            "moisture_content"
        } else {
            "moisture-content"
        };
        let reference = point(axis, 0.12);
        let law = ResolvedWoodMoisture::from_card(
            card,
            source(&keys, axis, &reference, &reference),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap_or_else(|e| panic!("{}: {e:?}", species.slug));
        close(law.fiber_saturation(), 0.30);
        let expected_source = [
            species.specific_gravity * 1120.0,
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
        let dry_mass_kg = species.specific_gravity * 1000.0 * 0.002;
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
                species.specific_gravity * 1000.0 * (1.0 + m) / factors[2],
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
            "{}: rho8={:.6},rho12={:.6},rho20={:.6} kg/m3; R20={:.9},T20={:.9},V20={:.9}",
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
        let wrong_reference = point(axis, 0.08);
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
