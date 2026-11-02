//! Independent Faraday/Ohm reduction; manufactured source data only.
use super::*;
use fs_material::{
    conductor::{ElectricalPropertyBasis, OhmicDrive, resolve_conductor_state_with_basis},
    field_storage::hysteresis::lamination::{LaminatedCoreEddyCurrent, ThinLaminationGeometry},
};

fn conducting_card(rho: f64) -> Arc<MaterialCard> {
    let original = card(|_| {});
    let mut claims = original.claims().clone();
    let basis = ElectricalPropertyBasis::Resistivity;
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(basis.property_name(), basis.dimensions()),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(250.0, 0.95 * rho), (600.0, 1.3 * rho)],
                dims: basis.dimensions(),
            },
            validity: ValidityDomain::unconstrained().with_quantity(
                "T",
                QuantitySpec::dimensional(Temperature::DIMS),
                250.0,
                600.0,
            ),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured lamination Ohmic source".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    Arc::new(MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap())
}

fn loss(card: Arc<MaterialCard>, turns: f64, thickness: f64) -> LaminatedCoreEddyCurrent {
    let core = UniformHystereticCore::from_source(resolve(card.clone()), turns, 0.5, 2.0).unwrap();
    let conductor = resolve_conductor_state_with_basis(
        &card,
        &point(300.0),
        MaterialPropertySelection::SingleClaimOnly,
        ElectricalPropertyBasis::Resistivity,
    )
    .unwrap();
    LaminatedCoreEddyCurrent::new(
        &core,
        conductor,
        ThinLaminationGeometry {
            thickness_m: thickness,
            width_m: 0.2,
            maximum_thickness_ratio: 0.1,
        },
    )
    .unwrap()
}

#[test]
fn native_shunt_matches_integrated_faraday_heat_and_geometry_substitution() {
    let rho = 2e-5;
    let card = conducting_card(rho);
    let a = loss(card.clone(), 3.0, 0.01);
    let b = loss(card.clone(), -3.0, 0.02);
    let c = loss(conducting_card(2.0 * rho), 3.0, 0.01);
    for (model, sigma) in [(&a, 1.0 / rho), (&b, 1.0 / rho), (&c, 0.5 / rho)] {
        let conductor = model.equivalent_conductor();
        let d = model.geometry().thickness_m;
        let b_dot: f64 = 0.7;
        let voltage = model.turns() * model.magnetic_area_m2() * b_dot;
        let current = conductor
            .current_a_for_drive(OhmicDrive::Voltage(voltage))
            .unwrap();
        // Integral y² dy from -d/2 to d/2, times sigma*width*length.
        let expected = model.volume_m3() * sigma * d.powi(2) * b_dot.powi(2) / 12.0;
        near(conductor.joule_power_w(current).unwrap(), expected);
        near(
            conductor.length_m() * conductor.area_m2(),
            model.volume_m3(),
        );
        near(model.magnetic_area_m2(), 0.5);
        near(model.path_length_m(), 2.0);
    }
    near(
        b.equivalent_conductor().resistance_ohm(),
        a.equivalent_conductor().resistance_ohm() / 4.0,
    );
    near(
        c.equivalent_conductor().resistance_ohm(),
        a.equivalent_conductor().resistance_ohm() * 2.0,
    );
    let original = a.equivalent_conductor();
    let snapshot = original
        .material()
        .at_query_point(&card, &point(350.0))
        .unwrap();
    let refreshed = a
        .with_equivalent_conductor(&original.with_material(snapshot).unwrap())
        .unwrap();
    near(
        refreshed.equivalent_conductor().resistance_ohm(),
        original.resistance_ohm() * 1.05,
    );
    assert_eq!(
        refreshed
            .equivalent_conductor()
            .source_property()
            .answer()
            .receipt
            .selected,
        original.source_property().answer().receipt.selected
    );
    near(original.material().query().axes()["T"], 300.0);
    assert!(
        a.with_equivalent_conductor(c.equivalent_conductor())
            .is_err()
    );
    assert!(
        a.with_equivalent_conductor(b.equivalent_conductor())
            .is_err()
    );
}

#[test]
fn lamination_requires_its_original_core_and_admitted_physical_sheet_geometry() {
    let card = conducting_card(2e-5);
    let core = UniformHystereticCore::from_source(resolve(card.clone()), 3.0, 0.5, 2.0).unwrap();
    let conductor = resolve_conductor_state_with_basis(
        &card,
        &point(300.0),
        MaterialPropertySelection::SingleClaimOnly,
        ElectricalPropertyBasis::Resistivity,
    )
    .unwrap();
    for geometry in [
        ThinLaminationGeometry {
            thickness_m: 0.0,
            width_m: 0.2,
            maximum_thickness_ratio: 0.1,
        },
        ThinLaminationGeometry {
            thickness_m: 0.03,
            width_m: 0.2,
            maximum_thickness_ratio: 0.1,
        },
        ThinLaminationGeometry {
            thickness_m: 0.01,
            width_m: 0.2,
            maximum_thickness_ratio: 1.0,
        },
    ] {
        assert!(LaminatedCoreEddyCurrent::new(&core, conductor.clone(), geometry).is_err());
    }
    let original = loss(card.clone(), 3.0, 0.01);
    assert!(original.matches_core(&core));
    let other = UniformHystereticCore::from_source(resolve(card), 6.0, 0.5, 2.0).unwrap();
    assert!(!original.matches_core(&other));
    let foreign = conducting_card(4e-5);
    let different = resolve_conductor_state_with_basis(
        &foreign,
        &point(300.0),
        MaterialPropertySelection::SingleClaimOnly,
        ElectricalPropertyBasis::Resistivity,
    )
    .unwrap();
    assert!(LaminatedCoreEddyCurrent::new(&core, different, original.geometry()).is_err());
}
