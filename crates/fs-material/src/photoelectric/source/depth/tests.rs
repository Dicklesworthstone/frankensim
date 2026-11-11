use super::*;
use crate::photoelectric::ELEMENTARY_CHARGE_C;
use crate::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyValue, Provenance,
    UncertaintyModel,
};
use fs_qty::{
    QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

fn temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", temperature(), t)
        .unwrap()
}
fn claim(key: PropertyKey, cold: f64, hot: f64) -> PropertyClaim {
    let dims = key.quantity().dims();
    PropertyClaim {
        key,
        value: PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: fs_qty::Temperature::DIMS,
            knots: vec![(300.0, cold), (400.0, hot)],
            dims,
        },
        validity: ValidityDomain::unconstrained().with_quantity("T", temperature(), 300.0, 400.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: Vec::new(),
        provenance: Provenance {
            source: "synthetic conditional escape fixture; not measured data".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}
fn fixture() -> (
    MaterialCard,
    ResolvedMaterialStatePoint,
    PropertyKey,
    PropertyKey,
    Vec<PhotoelectricDepthBin>,
) {
    let work = PropertyKey::new("surface-removal-work", Dims([2, 1, -2, 0, 0, 0]));
    let thickness = PropertyKey::new("calibrated-film-thickness", fs_qty::Length::DIMS);
    let mut claims = vec![
        claim(
            work.clone(),
            2.0 * ELEMENTARY_CHARGE_C,
            2.2 * ELEMENTARY_CHARGE_C,
        ),
        claim(thickness.clone(), 12e-9, 12e-9),
        claim(
            PropertyKey::new("density", fs_qty::Density::DIMS),
            3000.0,
            2900.0,
        ),
    ];
    let mut bins = Vec::new();
    for (b, (end_fraction, probability)) in [(0.25, 0.8), (1.0, 0.2)].into_iter().enumerate() {
        let mut knots = Vec::new();
        for (k, wavelength_nm) in [200.0, 400.0].into_iter().enumerate() {
            let yield_key = PropertyKey::new(format!("conditional-escape-{b}-{k}"), Dims::NONE);
            let energy_key = PropertyKey::new(format!("mean-excess-{b}-{k}"), Dims::NONE);
            let q = probability / (k + 1) as f64;
            claims.push(claim(yield_key.clone(), q, q * 0.5));
            claims.push(claim(energy_key.clone(), 0.5, 0.25));
            knots.push(PhotoelectricKnot {
                wavelength_nm,
                absorbed_quantum_yield_property: yield_key,
                mean_excess_fraction_property: energy_key,
            });
        }
        bins.push(PhotoelectricDepthBin {
            end_fraction,
            knots,
        });
    }
    let requirements: Vec<_> = claims
        .iter()
        .map(|c| {
            ScalarPropertyRequirement::try_with_key(&c.key, ScalarAdmissibility::Finite).unwrap()
        })
        .collect();
    let mut set = ClaimSet::new();
    for claim in claims {
        set.insert_claim(claim).unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic cathode".into(),
            phase: "solid".into(),
            process: "named surface".into(),
            revision: 0,
        },
        set,
        Vec::new(),
    )
    .unwrap();
    let resolved = resolve_material_state_point(
        &card,
        &point(300.0),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    (card, resolved, work, thickness, bins)
}

/// Every f64 a refresh must reproduce: the held state coordinates, each
/// resolved property value, the bound work function and thickness, and each
/// bin's end fraction and authored wavelengths.
fn physical_bits(material: &DepthPhotoelectricMaterial) -> Vec<u64> {
    let resolved = material.resolved();
    let mut values: Vec<f64> = Vec::new();
    values.extend(resolved.query_point().iter().map(|(_, value)| *value));
    values.extend(resolved.properties().iter().map(|p| p.value_si()));
    values.extend([material.work_function_j(), material.thickness_m()]);
    for bin in material.bins() {
        values.push(bin.end_fraction);
        values.extend(bin.knots.iter().map(|knot| knot.wavelength_nm));
    }
    values.into_iter().map(f64::to_bits).collect()
}

#[test]
fn g1_depth_conditional_projection_and_original_source_refresh() {
    let (card, resolved, work, thickness, bins) = fixture();
    let material = DepthPhotoelectricMaterial::try_from_resolved(
        resolved.clone(),
        work,
        thickness,
        &bins,
        200.0,
    )
    .unwrap();
    assert_eq!(material.resolved(), &resolved);
    assert_eq!(material.bin_bounds_fraction(0).unwrap(), [0.0, 0.25]);
    assert_eq!(material.bin_bounds_fraction(1).unwrap(), [0.25, 1.0]);
    assert_eq!(material.thickness_m(), 12e-9);
    let near = material.at_wavelength(0, 300.0).unwrap();
    let far = material.at_wavelength(1, 300.0).unwrap();
    assert!((near.absorbed_quantum_yield() - 0.6).abs() < 1e-15);
    assert!((far.absorbed_quantum_yield() - 0.15).abs() < 1e-15);
    assert_eq!(near.mean_excess_fraction(), 0.5);
    // Equal absorbed watts at different depths do not imply equal emission.
    let near_rate = near.evaluate(300e-9, 1.0).unwrap().emitted_electron_rate_s;
    let far_rate = far.evaluate(300e-9, 1.0).unwrap().emitted_electron_rate_s;
    assert!((near_rate / far_rate - 4.0).abs() < 1e-14);
    let hot = material.at_query_point(&card, &point(400.0)).unwrap();
    assert_eq!(
        hot.resolved(),
        &resolved.at_query_point(&card, &point(400.0)).unwrap()
    );
    assert_ne!(hot.identity(), material.identity());
    assert_eq!(
        hot.at_wavelength(0, 200.0)
            .unwrap()
            .absorbed_quantum_yield(),
        0.4
    );
    assert_eq!(
        hot.at_wavelength(1, 200.0).unwrap().mean_excess_fraction(),
        0.25
    );
    assert_eq!(
        hot.resolved().property("density").unwrap().value_si(),
        2900.0
    );
    // A refresh pins each selected claim; `material` was resolved with
    // SingleClaimOnly. The library states the rule in
    // src/liquid/density/tests.rs:285-287: "A refresh records pinned-claim
    // selection; an independent resolution records SingleClaimOnly. Their
    // physical tuples, not policy receipts, must agree exactly." So the restored
    // law equals the original refreshed through the same pinned path at its own
    // state, and every physical value equals the single-claim original bit for
    // bit; only the receipt policies, their provenance and the identities that
    // bind them differ.
    let restored = hot.at_query_point(&card, &point(300.0)).unwrap();
    assert_eq!(
        restored,
        material.at_query_point(&card, &point(300.0)).unwrap()
    );
    assert_eq!(restored.resolved().query(), material.resolved().query());
    assert_eq!(physical_bits(&restored), physical_bits(&material));
    assert!(material.at_query_point(&card, &point(401.0)).is_err());
    for (a, b) in resolved
        .properties()
        .iter()
        .zip(hot.resolved().properties())
    {
        assert_eq!(a.answer().receipt.selected, b.answer().receipt.selected);
        card.claims().verify_receipt(&b.answer().receipt).unwrap();
    }
}

#[test]
fn g0_depth_coverage_spectral_support_and_role_checks() {
    let (_, resolved, work, thickness, bins) = fixture();
    let bind = |bins: &[PhotoelectricDepthBin], gap| {
        DepthPhotoelectricMaterial::try_from_resolved(
            resolved.clone(),
            work.clone(),
            thickness.clone(),
            bins,
            gap,
        )
    };
    let material = bind(&bins, 200.0).unwrap();
    assert!(bind(&[], 200.0).is_err());
    for end in [0.0, -0.1, 1.1, f64::NAN, f64::INFINITY] {
        let mut bad = bins.clone();
        bad[0].end_fraction = end;
        assert!(bind(&bad, 200.0).is_err());
    }
    for end in [0.25, 0.9] {
        let mut bad = bins.clone();
        bad[1].end_fraction = end;
        assert!(bind(&bad, 200.0).is_err());
    }
    for nm in [199.0, 401.0, f64::NAN] {
        assert!(material.at_wavelength(0, nm).is_err());
    }
    assert!(material.at_wavelength(2, 300.0).is_err());
    assert!(material.bin_bounds_fraction(2).is_err());
    assert!(bind(&bins, 199.0).is_err());
    let mut bad = bins.clone();
    bad[1].knots[0].wavelength_nm = 201.0;
    assert!(bind(&bad, 200.0).is_err());
    let mut bad = bins.clone();
    bad[1].knots[0].mean_excess_fraction_property =
        bins[0].knots[0].absorbed_quantum_yield_property.clone();
    assert!(bind(&bad, 200.0).is_err());
    let mut bad = bins.clone();
    bad[1].knots[0].absorbed_quantum_yield_property =
        PropertyKey::new("density", fs_qty::Density::DIMS);
    assert!(bind(&bad, 200.0).is_err());
    // Explicitly repeated conditional source roles admit harmless refinement.
    let mut split = bins.clone();
    split.insert(
        0,
        PhotoelectricDepthBin {
            end_fraction: 0.125,
            knots: bins[0].knots.clone(),
        },
    );
    let refined = bind(&split, 200.0).unwrap();
    assert_eq!(
        refined.at_wavelength(0, 300.0).unwrap(),
        refined.at_wavelength(1, 300.0).unwrap()
    );
    assert_ne!(refined.identity(), material.identity());
}
