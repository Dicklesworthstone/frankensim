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
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

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

fn claim(key: PropertyKey, values: [f64; 2]) -> PropertyClaim {
    let dims = key.quantity().dims();
    PropertyClaim {
        key,
        value: PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: fs_qty::Temperature::DIMS,
            knots: vec![(300.0, values[0]), (400.0, values[1])],
            dims,
        },
        validity: ValidityDomain::unconstrained().with_quantity("T", temperature(), 300.0, 400.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: Vec::new(),
        provenance: Provenance {
            source: "synthetic prepared-film equation fixture; not measured photocathode data"
                .into(),
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
    Vec<PhotoelectricKnot>,
) {
    let work = PropertyKey::new("surface-removal-work", Dims([2, 1, -2, 0, 0, 0]));
    let thickness = PropertyKey::new("calibrated-film-thickness", fs_qty::Length::DIMS);
    let mut claims = vec![
        claim(
            work.clone(),
            [4.0 * ELEMENTARY_CHARGE_C, 4.2 * ELEMENTARY_CHARGE_C],
        ),
        claim(thickness.clone(), [12e-9; 2]),
        claim(
            PropertyKey::new("density", fs_qty::Density::DIMS),
            [3000.0, 2900.0],
        ),
    ];
    let mut knots = Vec::new();
    for (index, wavelength_nm) in [200.0, 300.0, 400.0].into_iter().enumerate() {
        let yield_key = PropertyKey::new(format!("absorbed-yield-{index}"), Dims::NONE);
        let energy_key = PropertyKey::new(format!("mean-excess-fraction-{index}"), Dims::NONE);
        claims.push(claim(
            yield_key.clone(),
            [0.2 + index as f64 * 0.1, 0.4 + index as f64 * 0.1],
        ));
        claims.push(claim(energy_key.clone(), [0.25, 0.5]));
        knots.push(PhotoelectricKnot {
            wavelength_nm,
            absorbed_quantum_yield_property: yield_key,
            mean_excess_fraction_property: energy_key,
        });
    }
    let requirements: Vec<_> = claims
        .iter()
        .map(|claim| {
            ScalarPropertyRequirement::try_with_key(&claim.key, ScalarAdmissibility::Finite)
                .unwrap()
        })
        .collect();
    let mut set = ClaimSet::new();
    for claim in claims {
        set.insert_claim(claim).unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic photocathode".into(),
            phase: "solid".into(),
            process: "named clean surface preparation".into(),
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
    (card, resolved, work, thickness, knots)
}

#[test]
fn g1_original_yield_energy_and_complete_state_survive_spectral_projection() {
    let (_, resolved, work, thickness, knots) = fixture();
    let model =
        PhotoelectricMaterial::try_from_resolved(resolved.clone(), work, thickness, &knots, 100.0)
            .unwrap();
    assert_eq!(model.resolved(), &resolved);
    assert_eq!(model.thickness_m(), 12e-9);
    assert_eq!(model.knots(), knots);
    assert_eq!(model.wavelength_bounds_nm(), [200.0, 400.0]);
    let law = model.at_wavelength(250.0).unwrap();
    assert_eq!(law.absorbed_quantum_yield(), 0.25);
    assert_eq!(law.mean_excess_fraction(), 0.25);
    let emission = law.evaluate(250e-9, 1.0).unwrap();
    assert!(emission.emitted_electron_rate_s > 0.0);
    let below = model
        .at_wavelength(400.0)
        .unwrap()
        .evaluate(400e-9, 100.0)
        .unwrap();
    assert_eq!(below.emitted_electron_rate_s, 0.0);
    assert_eq!(below.heat_power_w, 100.0);
    for nm in [199.0, 401.0, f64::NAN, f64::INFINITY] {
        assert!(model.at_wavelength(nm).is_err());
    }
}

/// Every f64 a refresh must reproduce: the held state coordinates, each
/// resolved property value, the bound work function and thickness, and each
/// authored wavelength.
fn physical_bits(material: &PhotoelectricMaterial) -> Vec<u64> {
    let resolved = material.resolved();
    let mut values: Vec<f64> = Vec::new();
    values.extend(resolved.query_point().iter().map(|(_, value)| *value));
    values.extend(resolved.properties().iter().map(|p| p.value_si()));
    values.extend([material.work_function_j(), material.thickness_m()]);
    values.extend(material.knots().iter().map(|knot| knot.wavelength_nm));
    values.into_iter().map(f64::to_bits).collect()
}

#[test]
fn g3_refresh_reuses_original_claims_and_refuses_unsupported_conditions() {
    let (card, resolved, work, thickness, knots) = fixture();
    let model =
        PhotoelectricMaterial::try_from_resolved(resolved.clone(), work, thickness, &knots, 100.0)
            .unwrap();
    let hot = model.at_query_point(&card, &point(400.0)).unwrap();
    assert_eq!(
        hot.resolved(),
        &resolved.at_query_point(&card, &point(400.0)).unwrap()
    );
    assert_ne!(hot.identity(), model.identity());
    assert_eq!(hot.work_function_j(), 4.2 * ELEMENTARY_CHARGE_C);
    assert_eq!(
        hot.at_wavelength(200.0).unwrap().absorbed_quantum_yield(),
        0.4
    );
    assert_eq!(
        hot.at_wavelength(200.0).unwrap().mean_excess_fraction(),
        0.5
    );
    assert_eq!(
        hot.resolved().property("density").unwrap().value_si(),
        2900.0
    );
    // A refresh pins each selected claim; `model` was resolved with
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
        model.at_query_point(&card, &point(300.0)).unwrap()
    );
    assert_eq!(restored.resolved().query(), model.resolved().query());
    assert_eq!(physical_bits(&restored), physical_bits(&model));
    assert!(model.at_query_point(&card, &point(401.0)).is_err());
    for (a, b) in model
        .resolved()
        .properties()
        .iter()
        .zip(hot.resolved().properties())
    {
        assert_eq!(a.answer().receipt.selected, b.answer().receipt.selected);
        card.claims().verify_receipt(&b.answer().receipt).unwrap();
    }
}

#[test]
fn g0_spectral_roles_units_gaps_and_missing_sources_refuse() {
    let (_, resolved, work, thickness, knots) = fixture();
    let bind = |knots: &[PhotoelectricKnot], gap| {
        PhotoelectricMaterial::try_from_resolved(
            resolved.clone(),
            work.clone(),
            thickness.clone(),
            knots,
            gap,
        )
    };
    for count in 0..2 {
        assert!(bind(&knots[..count], 100.0).is_err());
    }
    for gap in [0.0, 99.0, f64::NAN, f64::INFINITY] {
        assert!(bind(&knots, gap).is_err());
    }
    let mut duplicate = knots.clone();
    duplicate[1].absorbed_quantum_yield_property =
        duplicate[0].absorbed_quantum_yield_property.clone();
    assert!(bind(&duplicate, 100.0).is_err());
    let mut wrong = knots.clone();
    wrong[0].mean_excess_fraction_property = PropertyKey::new("missing-energy-source", Dims::NONE);
    assert!(bind(&wrong, 100.0).is_err());
    wrong[0].mean_excess_fraction_property = PropertyKey::new("density", fs_qty::Density::DIMS);
    assert!(bind(&wrong, 100.0).is_err());
    for wavelength in [0.0, 200.0, f64::NAN, f64::INFINITY] {
        wrong = knots.clone();
        wrong[1].wavelength_nm = wavelength;
        assert!(bind(&wrong, 100.0).is_err());
    }
    let per_mole = PropertyKey::new(work.name(), Dims([2, 1, -2, 0, 0, -1]));
    assert!(
        PhotoelectricMaterial::try_from_resolved(resolved, per_mole, thickness, &knots, 100.0)
            .is_err()
    );
}
