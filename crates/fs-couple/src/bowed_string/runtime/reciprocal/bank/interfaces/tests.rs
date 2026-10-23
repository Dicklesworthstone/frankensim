use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_qty::Density;

fn material(t: f64) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("density", Density::DIMS),
            value: PropertyValue::Scalar {
                value: 700.0,
                dims: Density::DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", 280.0, 310.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            provenance: Provenance {
                source: "synthetic footprint state".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
            observations: vec![],
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "same material name".into(),
            phase: "solid".into(),
            process: "fixture".into(),
            // MaterialCard::assemble starts a lineage at revision 0.
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_material_state_point(
        &card,
        &QueryPoint::new().with("T", t).unwrap(),
        &[ScalarPropertyRequirement::try_new(
            "density",
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

#[test]
fn physical_footprints_bind_incident_complete_states_and_include_shell_moments() {
    let a = material(293.15);
    let b = material(300.0);
    assert_eq!(a.material(), b.material());
    assert_ne!(a.identity(), b.identity());
    let triangles = [[0, 1, 2], [1, 3, 2]];
    for dofs in [1, 6] {
        // A shell's rotational-only moment arm still belongs to its real node.
        let offset = dofs - 1;
        let mut row = vec![0.0; 4 * dofs];
        row[offset] = 0.02;
        let mut selected = None;
        footprint(
            4,
            &triangles,
            &row,
            dofs,
            |e| Some(if e == 0 { &a } else { &b }),
            &mut selected,
            None,
        )
        .unwrap();
        assert_eq!(selected.unwrap().identity(), a.identity());
        row.fill(0.0);
        row[3 * dofs + offset] = 0.02;
        selected = None;
        footprint(
            4,
            &triangles,
            &row,
            dofs,
            |e| Some(if e == 0 { &a } else { &b }),
            &mut selected,
            None,
        )
        .unwrap();
        assert_eq!(selected.unwrap().identity(), b.identity());
        row.fill(0.0);
        row[dofs + offset] = 0.02; // Actual boundary node touches both states.
        selected = None;
        assert!(
            footprint(
                4,
                &triangles,
                &row,
                dofs,
                |e| Some(if e == 0 { &a } else { &b }),
                &mut selected,
                None
            )
            .is_err()
        );
        selected = None;
        footprint(4, &triangles, &row, dofs, |_| Some(&a), &mut selected, None).unwrap();
        assert_eq!(selected.unwrap().identity(), a.identity());
    }
    // Assembly parts use their ORIGINAL rows. An inactive part cannot select
    // its state, and later active parts cannot erase an earlier physical state.
    let row = [0.0; 24];
    let mut selected = None;
    footprint(4, &triangles, &row, 6, |_| Some(&a), &mut selected, None).unwrap();
    assert!(selected.is_none());
    let mut row = row;
    row[2] = 1.0;
    footprint(4, &triangles, &row, 6, |_| Some(&b), &mut selected, None).unwrap();
    assert_eq!(selected.unwrap().identity(), b.identity());
    assert!(footprint(4, &triangles, &row, 6, |_| Some(&a), &mut selected, None).is_err());
}

#[test]
fn absolute_temperature_uses_source_semantics_and_refuses_conflicting_coordinates() {
    let absolute = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let mut quantities = BTreeMap::new();
    quantities.insert("surface_kelvin".into(), absolute);
    assert_eq!(
        temperature(&[("surface_kelvin".into(), 293.15)], &quantities).unwrap(),
        Some(293.15)
    );
    assert!(
        temperature(
            &[("surface_kelvin".into(), 293.15), ("T".into(), 293.15)],
            &quantities
        )
        .is_err()
    );
    assert!(temperature(&[("T".into(), 0.0)], &BTreeMap::new()).is_err());
    quantities.insert("T".into(), QuantitySpec::dimensional(Density::DIMS));
    assert!(temperature(&[("T".into(), 293.15)], &quantities).is_err());
}
