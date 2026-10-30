//! Original strain sources change the native conductor geometry and exact tangent.
use super::*;
use fs_matdb::{ClaimSelection, SelectionPolicy};
use fs_material::thermal::transport::expansion::{
    FreeExpansionGeometryQuery, FreeExpansionStrainSources, PreparedFreeExpansion,
};
use fs_qty::Dims;

fn with_strains(
    card: MaterialCard,
    original: ResolvedConductorState,
    strains: &[Vec<(f64, f64)>],
) -> (MaterialCard, ResolvedConductorState, Vec<PropertyKey>) {
    let mut claims = card.claims().clone();
    let keys: Vec<_> = strains
        .iter()
        .enumerate()
        .map(|(i, knots)| {
            let key = PropertyKey::new(format!("engineering_length_change_{i}"), Dims::NONE);
            claims
                .insert_claim(PropertyClaim {
                    key: key.clone(),
                    value: PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots: knots.clone(),
                        dims: Dims::NONE,
                    },
                    interpolation: InterpolationPolicy::LinearInside,
                    validity: ValidityDomain::unconstrained().with_quantity(
                        "T",
                        absolute(),
                        300.0,
                        400.0,
                    ),
                    uncertainty: UncertaintyModel::Unstated,
                    observations: vec![],
                    provenance: Provenance {
                        source: "manufactured engineering length-change table".into(),
                        license: "test-only".into(),
                        artifact: None,
                    },
                })
                .unwrap();
            key
        })
        .collect();
    let card = MaterialCard::assemble(card.id().clone(), claims, vec![]).unwrap();
    let requirements: Vec<_> = original
        .material()
        .properties()
        .iter()
        .map(|property| property.requirement().clone())
        .collect();
    let material = resolve_material_state_point(
        &card,
        &original.material().query(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let original = ResolvedConductorState::try_from_material_with_property(
        material,
        original.source_basis(),
        original.source_key(),
    )
    .unwrap();
    (card, original, keys)
}

fn expansion(
    card: &MaterialCard,
    reference: &QueryPoint,
    keys: &[PropertyKey],
    range: [f64; 2],
) -> PreparedFreeExpansion {
    let selection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
    let strains = if keys.len() == 1 {
        FreeExpansionStrainSources::Isotropic {
            strain: &keys[0],
            selection,
        }
    } else {
        FreeExpansionStrainSources::FixedPrincipal {
            frame: fs_blake3::hash_bytes(b"manufactured orthogonal specimen axes"),
            strains: [&keys[0], &keys[1], &keys[2]],
            selections: [selection; 3],
        }
    };
    PreparedFreeExpansion::try_new(
        card,
        FreeExpansionGeometryQuery {
            reference,
            strains,
            temperature_axis: "T",
            lower: &reference.clone().with("T", range[0]).unwrap(),
            upper: &reference.clone().with("T", range[1]).unwrap(),
        },
    )
    .unwrap()
}

#[test]
fn native_and_graph_geometry_share_the_same_pre_resolved_temperature_bundle() {
    let (card, original) = curve(
        ElectricalPropertyBasis::Conductivity,
        vec![(300.0, 0.5), (400.0, 0.25)],
    );
    let (card, original, keys) = with_strains(card, original, &[vec![(300.0, 0.1), (400.0, 0.3)]]);
    let path = expansion(&card, &original.material().query(), &keys, [300.0, 400.0]);
    let law = MaterialOhmicNode::try_new(&card, original, "T", [300.0, 400.0], 2.0, 0.5)
        .unwrap()
        .with_free_expansion(path, 0)
        .unwrap();
    let query = law
        .source()
        .original()
        .material()
        .query()
        .with("T", 375.0)
        .unwrap();
    let state = law.source().at_query_point(&query).unwrap();
    let bundle = state.material().clone();
    let native = law.resolve_from_state(state).unwrap();
    assert_eq!(native.material(), &bundle);
    assert_eq!(native, law.resolve_at_temperature_k(375.0).unwrap());
    near(native.length_m(), 2.0 * 1.25 / 1.1);
    near(native.area_m2(), 0.5 * (1.25_f64 / 1.1).powi(2));
    // A canonical source refresh may resolve outside this consumer's declared
    // temperature slice. Reusing its material receipt must not bypass that slice.
    let narrow = MaterialOhmicNode::try_new(
        &card,
        law.source().original().clone(),
        "T",
        [300.0, 350.0],
        2.0,
        0.5,
    )
    .unwrap();
    assert!(
        narrow
            .resolve_from_state(law.source().at_query_point(&query).unwrap())
            .is_err()
    );
    let (foreign, state) = curve(
        ElectricalPropertyBasis::Conductivity,
        vec![(300.0, 0.5), (400.0, 0.25)],
    );
    assert_ne!(foreign.content_hash(), card.content_hash());
    assert!(law.resolve_from_state(state).is_err());
}

#[test]
fn original_strain_tangent_keeps_fitted_datum_and_one_sided_source_knots() {
    let (card, original) = curve(
        ElectricalPropertyBasis::Resistivity,
        vec![(300.0, 2.0), (400.0, 3.0)],
    );
    let (card, original, keys) = with_strains(
        card,
        original,
        &[
            vec![(300.0, 0.2), (350.0, 0.5), (400.0, 1.1)],
            vec![(300.0, 0.1), (350.0, 0.2), (400.0, 0.3)],
            vec![(300.0, -0.1), (350.0, 0.0), (400.0, 0.3)],
        ],
    );
    let path = expansion(&card, &original.material().query(), &keys, [300.0, 400.0]);
    let original_path = path.clone();
    assert_eq!(
        path.at_temperature(300.0).unwrap().principal_stretches(),
        [1.0; 3]
    );
    for (temperature, strains, slopes) in [
        (300.0, [0.2, 0.1, -0.1], [0.006, 0.002, 0.002]),
        (340.0, [0.44, 0.18, -0.02], [0.006, 0.002, 0.002]),
        (350.0, [0.5, 0.2, 0.0], [0.012, 0.002, 0.006]),
        (400.0, [1.1, 0.3, 0.3], [0.012, 0.002, 0.006]),
    ] {
        let reference = [1.2, 1.1, 0.9];
        let stretches: [f64; 3] = std::array::from_fn(|i| (1.0 + strains[i]) / reference[i]);
        let derivatives: [f64; 3] = std::array::from_fn(|i| slopes[i] / reference[i]);
        let tangent = path.at_temperature_with_tangent(temperature).unwrap();
        assert_eq!(
            tangent.geometry(),
            path.at_temperature(temperature).unwrap()
        );
        assert_eq!(tangent.geometry().path_identity(), path.identity());
        for i in 0..3 {
            near(tangent.geometry().principal_stretches()[i], stretches[i]);
            near(
                tangent.principal_stretch_derivatives_per_k()[i],
                derivatives[i],
            );
        }
        near(
            tangent.relative_volume_derivative_per_k(),
            derivatives[0] * stretches[1] * stretches[2]
                + stretches[0] * derivatives[1] * stretches[2]
                + stretches[0] * stretches[1] * derivatives[2],
        );
    }
    assert_eq!(path, original_path);
    assert!(path.at_temperature_with_tangent(299.0).is_err());
    assert!(path.at_temperature_with_tangent(401.0).is_err());
    assert!(path.at_temperature_with_tangent(f64::NAN).is_err());
}

#[test]
fn principal_geometry_and_each_electrical_basis_share_the_native_chain_rule() {
    for basis in [
        ElectricalPropertyBasis::Resistivity,
        ElectricalPropertyBasis::Conductivity,
    ] {
        let (card, original) = curve(basis, vec![(300.0, 2.0), (400.0, 3.0)]);
        let (card, original, keys) = with_strains(
            card,
            original,
            &[
                vec![(300.0, 0.2), (400.0, 0.68)],
                vec![(300.0, -0.1), (400.0, 0.44)],
                vec![(300.0, 0.1), (400.0, -0.01)],
            ],
        );
        let path = expansion(&card, &original.material().query(), &keys, [300.0, 400.0]);
        for axis in 0..3 {
            let law = node(&card, original.clone())
                .with_free_expansion(path.clone(), axis)
                .unwrap();
            let native = law.resolve_at_temperature_k(350.0).unwrap();
            let stretches = [1.2, 1.3, 0.95];
            let slopes = [0.004, 0.006, -0.001];
            let transverse = [(axis + 1) % 3, (axis + 2) % 3];
            let length = 2.0 * stretches[axis];
            let area = 0.5 * stretches[transverse[0]] * stretches[transverse[1]];
            let (rho, material_relative) = match basis {
                ElectricalPropertyBasis::Resistivity => (2.5, 0.01 / 2.5),
                ElectricalPropertyBasis::Conductivity => (1.0 / 2.5, -0.01 / 2.5),
            };
            let resistance = rho * length / area;
            let derivative = resistance
                * (material_relative + slopes[axis] / stretches[axis]
                    - slopes[transverse[0]] / stretches[transverse[0]]
                    - slopes[transverse[1]] / stretches[transverse[1]]);
            near(native.length_m(), length);
            near(native.area_m2(), area);
            near(native.resistance_ohm(), resistance);
            assert_eq!(native.source_basis(), basis);
            assert_eq!(native.source_key(), original.source_key());
            assert_eq!(
                native.source_property().answer().receipt.selected,
                original.source_property().answer().receipt.selected
            );
            // Electrical geometry cannot rewrite an independent density claim.
            near(
                native.material().property("density").unwrap().value_si(),
                8500.0,
            );
            assert_eq!(law.free_expansion().unwrap().0, &path);
            let tangent = law.tangent(&[], &[50.0, -3.0]).unwrap();
            for (actual, expected) in tangent.into_iter().zip([
                -3.0 * derivative,
                resistance,
                9.0 * derivative,
                -6.0 * resistance,
            ]) {
                near(actual, expected);
            }
            check_consistent_tangent("expanding wire", &law, &[], &[47.0, -3.0], 1e-7).unwrap();
            let mut graph = ConstitutiveGraph::new();
            graph.add_node("wire", Box::new(law)).unwrap();
            let value = graph
                .execute(
                    graph.state_schema().version(),
                    &[],
                    &[
                        (("wire".into(), "temperature_change".into()), 50.0),
                        (("wire".into(), "current".into()), -3.0),
                    ]
                    .into(),
                )
                .unwrap();
            near(
                value.outputs[&("wire".into(), "voltage".into())],
                -3.0 * resistance,
            );
            near(
                value.outputs[&("wire".into(), "heat".into())],
                9.0 * resistance,
            );
            near(value.total_dissipation, 9.0 * resistance);
        }
    }
}

#[test]
fn true_material_geometry_cancellation_and_constant_geometry_remain_admissible() {
    for final_strain in [0.0, 0.5] {
        let (card, original) = curve(
            ElectricalPropertyBasis::Resistivity,
            vec![(300.0, 2.0), (400.0, 3.0)],
        );
        let (card, original, keys) =
            with_strains(card, original, &[vec![(300.0, 0.0), (400.0, final_strain)]]);
        let fixed = node(&card, original.clone());
        let law = fixed
            .clone()
            .with_free_expansion(
                expansion(&card, &original.material().query(), &keys, [300.0, 400.0]),
                0,
            )
            .unwrap();
        for change in [0.0, 25.0, 50.0, 100.0] {
            let values = law.evaluate(&[], &[change, 2.0]).unwrap();
            let tangent = law.tangent(&[], &[change, 2.0]).unwrap();
            if final_strain == 0.0 {
                assert_eq!(
                    values.outputs,
                    fixed.evaluate(&[], &[change, 2.0]).unwrap().outputs
                );
                assert_eq!(tangent, fixed.tangent(&[], &[change, 2.0]).unwrap());
            } else {
                // rho/rho_ref and isotropic stretch rise together: R stays 8 ohm.
                near(values.outputs[0], 16.0);
                near(values.outputs[1], 32.0);
                near(tangent[0], 0.0);
                near(tangent[2], 0.0);
            }
        }
    }
}

#[test]
fn expansion_cannot_cross_material_reference_range_or_fixed_context_ownership() {
    let (card, original) = surface(false);
    let strains = [vec![(300.0, 0.0), (400.0, 0.1)]];
    let (card, original, keys) = with_strains(card, original, &strains);
    let reference = original.material().query();
    let path = expansion(&card, &reference, &keys, [300.0, 400.0]);
    let original_path = path.clone();
    let law = node(&card, original.clone());
    assert!(law.clone().with_free_expansion(path.clone(), 3).is_err());
    let bound = law.clone().with_free_expansion(path.clone(), 1).unwrap();
    assert!(bound.clone().with_free_expansion(path.clone(), 1).is_err());
    for (context, range) in [
        (reference.clone().with("T", 325.0).unwrap(), [300.0, 400.0]),
        (
            reference.clone().with("pressure", 1.8e5).unwrap(),
            [300.0, 400.0],
        ),
        (reference.clone(), [300.0, 390.0]),
    ] {
        let changed = expansion(&card, &context, &keys, range);
        assert!(law.clone().with_free_expansion(changed, 0).is_err());
    }
    let (foreign, foreign_state) = surface(false);
    let (foreign, foreign_state, foreign_keys) =
        with_strains(foreign, foreign_state, &[vec![(300.0, 0.0), (400.0, 0.2)]]);
    let foreign_path = expansion(
        &foreign,
        &foreign_state.material().query(),
        &foreign_keys,
        [300.0, 400.0],
    );
    assert!(law.with_free_expansion(foreign_path, 0).is_err());
    assert!(bound.evaluate(&[], &[101.0, 1.0]).is_err());
    assert_eq!(bound.free_expansion().unwrap().0, &original_path);
    assert_eq!(original.material().query(), reference);
}

#[test]
fn overlapping_strain_requirement_keeps_its_original_selected_claim() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let (card, original) = curve(basis, vec![(300.0, 2.0), (400.0, 3.0)]);
    let (card, original, keys) = with_strains(card, original, &[vec![(300.0, 0.0), (400.0, 0.1)]]);
    let (first, strain_claim) = card.claims_for(keys[0].name())[0];
    let mut competitor = strain_claim.clone();
    competitor.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(300.0, 0.0), (400.0, 0.2)],
        dims: Dims::NONE,
    };
    let mut claims = card.claims().clone();
    let second = claims.insert_claim(competitor).unwrap();
    let card = MaterialCard::assemble(card.id().clone(), claims, vec![]).unwrap();
    let mut requirements: Vec<_> = original
        .material()
        .properties()
        .iter()
        .map(|property| property.requirement().clone())
        .collect();
    requirements.push(
        ScalarPropertyRequirement::try_with_key(&keys[0], ScalarAdmissibility::Finite).unwrap(),
    );
    let mut pins: Vec<_> = original
        .material()
        .properties()
        .iter()
        .map(|property| {
            (
                property.requirement().name().to_owned(),
                property.answer().receipt.selected,
            )
        })
        .collect();
    pins.push((keys[0].name().to_owned(), first));
    let reference = original.material().query();
    let material = resolve_material_state_point(
        &card,
        &reference,
        &requirements,
        MaterialPropertySelection::PinnedByProperty(pins),
    )
    .unwrap();
    let original = ResolvedConductorState::try_from_material_with_property(
        material,
        basis,
        original.source_key(),
    )
    .unwrap();
    let law = node(&card, original);
    for selected in [first, second] {
        let path = PreparedFreeExpansion::try_new(
            &card,
            FreeExpansionGeometryQuery {
                reference: &reference,
                strains: FreeExpansionStrainSources::Isotropic {
                    strain: &keys[0],
                    selection: ClaimSelection::Pinned(selected),
                },
                temperature_axis: "T",
                lower: &reference,
                upper: &reference.clone().with("T", 400.0).unwrap(),
            },
        )
        .unwrap();
        let result = law.clone().with_free_expansion(path, 0);
        assert_eq!(result.is_ok(), selected == first);
        if let Ok(admitted) = result {
            let resolved = admitted.resolve_at_temperature_k(350.0).unwrap();
            assert_eq!(
                resolved
                    .material()
                    .property_by_key(&keys[0])
                    .unwrap()
                    .answer()
                    .receipt
                    .selected,
                first
            );
            near(resolved.length_m(), 2.1);
        }
    }
}
