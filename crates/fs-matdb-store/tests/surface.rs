//! Discovery respects source patches and delegates every value to fs-matdb.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, NormalizedPack, ObservationDataset, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_matdb_store::{CatalogPack, MaterialStore};
use fs_qty::{Dims, QuantitySpec};

fn pack(name: &str, interpolation: InterpolationPolicy) -> NormalizedPack {
    let source_artifact = fs_blake3::hash_bytes(
        b"synthetic disconnected patches\nx,y,response\n0,0,10\n1,0,20\n0,1,10\n3,0,40\n4,0,50\n3,1,40\ntriangles\n0,1,2\n3,4,5\n",
    );
    let provenance = Provenance {
        source: "manufactured disconnected source patches; not measured material data".into(),
        license: "test-only".into(),
        artifact: Some(source_artifact),
    };
    let mut claims = ClaimSet::new();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "two independent dimensionless source patches".into(),
            method: "authored sample table with explicit triangle connectivity".into(),
            artifact: source_artifact,
            caveats: "deterministic synthetic fixture; no experimental uncertainty claim".into(),
            provenance: provenance.clone(),
        })
        .unwrap();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("surface_response", Dims::NONE),
            value: PropertyValue::Surface2D {
                axes: ["x".into(), "y".into()],
                axis_dims: [Dims::NONE; 2],
                samples: vec![
                    ([0.0, 0.0], 10.0),
                    ([1.0, 0.0], 20.0),
                    ([0.0, 1.0], 10.0),
                    ([3.0, 0.0], 40.0),
                    ([4.0, 0.0], 50.0),
                    ([3.0, 1.0], 40.0),
                ],
                triangles: vec![[0, 1, 2], [3, 4, 5]],
                dims: Dims::NONE,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity("x", QuantitySpec::dimensional(Dims::NONE), 0.0, 4.0)
                .with_quantity("y", QuantitySpec::dimensional(Dims::NONE), 0.0, 1.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation,
            observations: vec![observation],
            provenance,
        })
        .unwrap();
    NormalizedPack::new(
        name,
        "store-surface-test",
        source_artifact,
        "synthetic redistribution permitted for tests",
        claims,
        vec![],
        vec![],
    )
    .unwrap()
}

#[test]
fn g0_surface_discovery_preserves_disconnected_ranges_and_canonical_evaluation() {
    let linear = pack("surface-linear", InterpolationPolicy::LinearInside);
    let exact = pack("surface-exact", InterpolationPolicy::TabulatedOnly);
    let store = MaterialStore::open(":memory:").unwrap();
    store
        .ingest_bundle(&[
            CatalogPack::Properties(linear.clone()),
            CatalogPack::Properties(exact),
        ])
        .unwrap();
    let seal = store.seal_corpus().unwrap();
    for name in ["surface-linear", "surface-exact"] {
        store.verify_index(name).unwrap();
    }
    let inside = store
        .materials_with("surface_response", Some((15.0, 15.0)))
        .unwrap();
    assert_eq!(inside.len(), 1);
    assert_eq!(inside[0].pack_id, "surface-linear");
    assert_eq!(inside[0].kind, "surface2d");
    assert_eq!(inside[0].scalar_value, None);
    assert_eq!(
        store
            .materials_with("surface_response", Some((10.0, 10.0)))
            .unwrap()
            .len(),
        2
    );
    // Overall extrema are [10,50], but neither declared patch attains [25,35].
    assert!(
        store
            .materials_with("surface_response", Some((25.0, 35.0)))
            .unwrap()
            .is_empty()
    );
    let point = QueryPoint::new()
        .with_quantity("x", QuantitySpec::dimensional(Dims::NONE), 0.5)
        .unwrap()
        .with_quantity("y", QuantitySpec::dimensional(Dims::NONE), 0.25)
        .unwrap();
    let key = PropertyKey::new("surface_response", Dims::NONE);
    let direct = linear
        .claims()
        .query_typed(&key, &point, SelectionPolicy::SingleClaimOnly)
        .unwrap();
    let stored = store
        .evaluate_typed(
            "surface-linear",
            &key,
            &point,
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
    assert_eq!(stored, direct);
    assert_eq!(stored.evidence.value.value, 15.0);
    assert!(
        store
            .evaluate_typed(
                "surface-exact",
                &key,
                &point,
                SelectionPolicy::SingleClaimOnly
            )
            .is_err()
    );
    let gap = point.with("x", 2.0).unwrap();
    assert!(
        store
            .evaluate_typed(
                "surface-linear",
                &key,
                &gap,
                SelectionPolicy::SingleClaimOnly
            )
            .is_err()
    );
    assert_eq!(store.seal_corpus().unwrap(), seal);
}
