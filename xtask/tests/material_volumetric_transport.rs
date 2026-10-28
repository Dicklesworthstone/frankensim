//! Source-bound thermal transport needs volumetric storage, not invented density.
#![deny(unsafe_code)]

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    NormalizedMaterialCardPack, NormalizedPack, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::thermal::transport::{
    ThermalTransportError,
    volumetric::{
        ResolvedVolumetricThermalTransport, VOLUMETRIC_HEAT_CAPACITY_DIMS as C,
        VolumetricThermalSourceQuery,
    },
};
use fs_qty::{
    Dims, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::{path::Path, process::Command};

const K: Dims = Dims([1, 1, -3, -1, 0, 0]);
const CP: Dims = Dims([2, 0, -2, -1, 0, 0]);
const POLICY: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);

fn semantic(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn keys() -> [PropertyKey; 2] {
    [
        PropertyKey::new("volumetric-heat-capacity", C),
        PropertyKey::with_quantity(
            "thermal-conductivity",
            semantic(QuantityKind::ThermalConductivity),
        ),
    ]
}
fn point(t: f64, moisture: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "temperature",
            semantic(QuantityKind::AbsoluteTemperature),
            t,
        )
        .unwrap()
        .with_quantity("moisture", QuantitySpec::dimensional(Dims::NONE), moisture)
        .unwrap()
}
fn resolve(
    card: &MaterialCard,
    point: &QueryPoint,
    keys: &[PropertyKey; 2],
    selections: [ClaimSelection; 2],
) -> Result<ResolvedVolumetricThermalTransport, ThermalTransportError> {
    ResolvedVolumetricThermalTransport::try_new(
        card,
        VolumetricThermalSourceQuery {
            volumetric_heat_capacity: &keys[0],
            conductivity: &keys[1],
            temperature_axis: "temperature",
            point,
            selections,
        },
    )
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-13 * expected.abs(),
        "{actual:e} != {expected:e}"
    );
}
fn identity(revision: u32) -> MaterialStateId {
    MaterialStateId {
        chemistry: "source-bound volumetric transport".into(),
        phase: "solid".into(),
        process: "explicit source conditions".into(),
        revision,
    }
}
fn compiled(slug: &str) -> NormalizedMaterialCardPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let path = std::env::temp_dir().join(format!(
        "volumetric-source-{}-{nonce}.fsmatpk",
        std::process::id()
    ));
    let result = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(root.join(format!("data/matdb/seed-v1/{slug}/manifest.tsv")))
        .arg("--out")
        .arg(&path)
        .env("CARGO_WORKSPACE_DIR", root)
        .current_dir(root)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{slug}: {}",
        String::from_utf8_lossy(&result.stderr)
    );
    let bytes = std::fs::read(path).unwrap();
    let pack = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(pack.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes);
    NormalizedMaterialCardPack::new(identity(0), pack).unwrap()
}

#[test]
fn g1_g3_real_tps_packs_supply_transport_without_mass_specific_inputs() {
    let keys = keys();
    for (slug, k, capacity, reported_alpha) in [
        ("loctite-aa-h4500-thermtest-23c", 0.200, 1.717e6, 0.117e-6),
        (
            "prima-solder-eg8020-thermtest-23c",
            1.905,
            1.030e6,
            1.849e-6,
        ),
        (
            "polypropylene-mtps-thermtest-26c-reference",
            0.280,
            1.536e6,
            1.83e-7,
        ),
    ] {
        let pack = compiled(slug);
        let card = pack.card();
        let source = card
            .claims()
            .claims_ordered()
            .find(|(_, claim)| {
                claim.key == keys[0]
                    && claim
                        .validity
                        .bounds()
                        .get("source-statistic-is-mean")
                        .is_none_or(|&(lo, hi)| lo == 1.0 && hi == 1.0)
            })
            .unwrap()
            .1;
        // Select the exact published singleton context, including every descriptor.
        let p = source
            .validity
            .bounds()
            .iter()
            .fold(QueryPoint::new(), |p, (axis, &(lo, hi))| {
                assert_eq!(lo, hi);
                p.with_quantity(
                    axis,
                    *source.validity.axis_quantities().get(axis).unwrap(),
                    lo,
                )
                .unwrap()
            });
        let model = resolve(card, &p, &keys, [POLICY; 2]).unwrap();
        assert_eq!(model.material_card_identity(), card.content_hash());
        assert_eq!(model.query_point(), &p);
        close(model.volumetric_heat_capacity_j_m3_k(), capacity);
        close(model.conductivity_w_m_k(), k);
        close(model.diffusivity_m2_s(), k / capacity);
        close(model.effusivity_j_m2_k_sqrt_s(), (k * capacity).sqrt());
        for (key, answer) in keys.iter().zip(model.sources()) {
            let original = card
                .claims()
                .query_typed(key, &p, SelectionPolicy::SingleClaimOnly)
                .unwrap();
            assert_eq!(answer, &original);
            assert_eq!(
                answer.evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
            assert!(answer.receipt.observation_backed);
            card.claims().verify_receipt(&answer.receipt).unwrap();
        }
        let alpha_key = PropertyKey::new("thermal-diffusivity", Dims([2, 0, -1, 0, 0, 0]));
        let reported = card
            .claims()
            .query_typed(&alpha_key, &p, SelectionPolicy::SingleClaimOnly)
            .unwrap();
        close(reported.evidence.value.value, reported_alpha);
        assert_ne!(model.diffusivity_m2_s(), reported.evidence.value.value);
        assert!((model.diffusivity_m2_s() / reported_alpha - 1.0).abs() < 0.005);
        let again = model.at_query_point(card, &p).unwrap();
        close(again.diffusivity_m2_s(), model.diffusivity_m2_s());
        assert_eq!(
            card.claims()
                .query_typed(&alpha_key, &p, SelectionPolicy::SingleClaimOnly)
                .unwrap(),
            reported
        );
        for (name, dims) in [
            ("density", Dims([-3, 1, 0, 0, 0, 0])),
            ("specific-heat-capacity", CP),
        ] {
            assert!(
                card.claims()
                    .query_typed(
                        &PropertyKey::new(name, dims),
                        &p,
                        SelectionPolicy::SingleClaimOnly
                    )
                    .is_err()
            );
        }
        let unsupported = p
            .clone()
            .with("temperature", p.axes()["temperature"] + 0.01)
            .unwrap();
        assert!(model.at_query_point(card, &unsupported).is_err());
    }
}

fn claim(index: usize, value: PropertyValue, high: f64) -> PropertyClaim {
    PropertyClaim {
        key: keys()[index].clone(),
        value,
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                "temperature",
                semantic(QuantityKind::AbsoluteTemperature),
                300.0,
                high,
            )
            .with_quantity("moisture", QuantitySpec::dimensional(Dims::NONE), 0.0, 1.0),
        interpolation: InterpolationPolicy::LinearInside,
        uncertainty: UncertaintyModel::RelativeHalfWidth {
            fraction: 0.02,
            confidence: 0.95,
        },
        observations: Vec::new(),
        provenance: Provenance {
            source: "authored numerical fixture; not measured data".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}
fn curve(index: usize, low: f64, high: f64) -> PropertyClaim {
    claim(
        index,
        PropertyValue::Curve {
            abscissa: "temperature".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(300.0, low), (500.0, high)],
            dims: [C, K][index],
        },
        500.0,
    )
}
fn card(claims: impl IntoIterator<Item = PropertyClaim>) -> MaterialCard {
    let mut set = ClaimSet::new();
    for claim in claims {
        set.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(identity(0), set, Vec::new()).unwrap()
}

#[test]
fn g3_multiaxis_refresh_keeps_exact_pins_and_rejects_switching_or_changed_schema() {
    let conductivity = claim(
        1,
        PropertyValue::Surface2D {
            axes: ["temperature".into(), "moisture".into()],
            axis_dims: [Temperature::DIMS, Dims::NONE],
            samples: vec![
                ([300.0, 0.0], 1.0),
                ([400.0, 0.0], 3.0),
                ([300.0, 1.0], 2.0),
                ([400.0, 1.0], 4.0),
            ],
            triangles: vec![[0, 1, 2], [1, 3, 2]],
            dims: K,
        },
        400.0,
    );
    let mut alternate = curve(1, 50.0, 60.0);
    alternate.validity = alternate.validity.with_quantity(
        "temperature",
        semantic(QuantityKind::AbsoluteTemperature),
        400.0,
        500.0,
    );
    alternate.value = PropertyValue::Curve {
        abscissa: "temperature".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(400.0, 50.0), (500.0, 60.0)],
        dims: K,
    };
    let card = card([curve(0, 2e6, 6e6), conductivity, alternate]);
    let keys = keys();
    let start = point(325.0, 0.25);
    let selected = keys.each_ref().map(|key| {
        card.claims()
            .query_typed(key, &start, SelectionPolicy::SingleClaimOnly)
            .unwrap()
            .receipt
            .selected
    });
    let model = resolve(&card, &start, &keys, selected.map(ClaimSelection::Pinned)).unwrap();
    let before = model.clone();
    let next = point(375.0, 0.75);
    let updated = model.at_query_point(&card, &next).unwrap();
    close(updated.volumetric_heat_capacity_j_m3_k(), 3.5e6);
    close(updated.conductivity_w_m_k(), 3.25);
    close(updated.diffusivity_m2_s(), 3.25 / 3.5e6);
    assert_ne!(updated.identity(), model.identity());
    for (answer, id) in updated.sources().iter().zip(selected) {
        assert_eq!(answer.receipt.selected, id);
        assert_eq!(
            answer.evidence.value.uncertainty,
            UncertaintyModel::RelativeHalfWidth {
                fraction: 0.02,
                confidence: 0.95
            }
        );
        card.claims().verify_receipt(&answer.receipt).unwrap();
    }
    // Ordinary selection could use the alternate source at450 K. Refresh must not.
    let switched = point(450.0, 0.75);
    assert!(
        card.claims()
            .query_typed(&keys[1], &switched, SelectionPolicy::SingleClaimOnly)
            .is_ok()
    );
    assert!(model.at_query_point(&card, &switched).is_err());
    for invalid in [
        point(600.0, 0.75),
        point(375.0, 1.01),
        next.clone().with("extra", 1.0).unwrap(),
        next.clone()
            .with_quantity(
                "temperature",
                QuantitySpec::dimensional(Temperature::DIMS),
                375.0,
            )
            .unwrap(),
    ] {
        assert!(model.at_query_point(&card, &invalid).is_err());
    }
    let changed_id = MaterialStateId {
        process: "different source material".into(),
        ..identity(0)
    };
    let changed = MaterialCard::assemble(changed_id, card.claims().clone(), Vec::new()).unwrap();
    assert!(model.at_query_point(&changed, &next).is_err());
    assert_eq!(model, before);
}

#[test]
fn g0_quantity_static_support_and_representability_are_admitted_before_publication() {
    let keys = keys();
    let p = point(350.0, 0.5);
    let supported = card([curve(0, 2e6, 2e6), curve(1, 1.0, 1.0)]);
    let mut mass_specific = keys.clone();
    mass_specific[0] = PropertyKey::new("specific-heat-capacity", CP);
    assert!(matches!(
        resolve(&supported, &p, &mass_specific, [POLICY; 2]),
        Err(ThermalTransportError::Invalid(_))
    ));
    let dynamic = p
        .clone()
        .with_quantity(
            "test-frequency",
            QuantitySpec::dimensional(Dims([0, 0, -1, 0, 0, 0])),
            1000.0,
        )
        .unwrap();
    let difference = p
        .clone()
        .with_quantity(
            "temperature",
            semantic(QuantityKind::TemperatureDifference),
            350.0,
        )
        .unwrap();
    for bad in [dynamic, difference, point(299.0, 0.5)] {
        assert!(resolve(&supported, &bad, &keys, [POLICY; 2]).is_err());
    }
    let mut unbounded = curve(0, 2e6, 2e6);
    unbounded.value = PropertyValue::Scalar {
        value: 2e6,
        dims: C,
    };
    unbounded.interpolation = InterpolationPolicy::ConstantWithinValidity;
    unbounded.validity = ValidityDomain::unconstrained();
    assert!(
        resolve(
            &card([unbounded, curve(1, 1.0, 1.0)]),
            &p,
            &keys,
            [POLICY; 2]
        )
        .is_err()
    );
    for (capacity, k) in [
        (0.0, 1.0),
        (-1.0, 1.0),
        (1.0, 0.0),
        (1.0, -1.0),
        (1e308, 1e-300),
        (1e-300, 1e308),
    ] {
        let bad = card([curve(0, capacity, capacity), curve(1, k, k)]);
        assert!(resolve(&bad, &p, &keys, [POLICY; 2]).is_err());
    }
    let large = card([curve(0, 1e200, 1e200), curve(1, 1e200, 1e200)]);
    let large = resolve(&large, &p, &keys, [POLICY; 2]).unwrap();
    close(large.diffusivity_m2_s(), 1.0);
    close(large.effusivity_j_m2_k_sqrt_s(), 1e200);
}
