//! Fixed environmental data and both live body temperatures on the same contact.
//! Synthetic constitutive laws isolate wiring and conservation, not qualification.
use super::*;
use fs_conduction::transient::deposition::contact::{
    ContactCell, MaterialContactLink, apply_material_contacts,
};
use fs_conduction::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RD, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RP,
};
use fs_evidence::ValidityDomain;
use fs_material::electrical_interface::{
    ELECTRICAL_CONTACT_RESISTANCE_DIMS as ED, ELECTRICAL_CONTACT_RESISTANCE_PROPERTY as EP,
};
use fs_qty::semantic::{CompositionBasis, QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec, Temperature};

fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

fn claim(name: &str, dims: Dims, value: PropertyValue, validity: ValidityDomain) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value,
        validity,
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic environmental contact response".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}

fn interface(cards: &[MaterialCard; 2], claims: ClaimSet) -> InterfaceSystemCard {
    InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: cards[0].id().clone(),
            texture_frame: "first finish".into(),
        },
        SurfaceSpec {
            material: cards[1].id().clone(),
            texture_frame: "second finish".into(),
        },
        SystemContext {
            medium: "authored contact".into(),
            third_body: None,
            environment: "explicit fixed conditions".into(),
            history: "fixed".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn validity() -> ValidityDomain {
    ValidityDomain::unconstrained()
        .with("Ta", 200.0, 500.0)
        .with("Tb", 200.0, 500.0)
}

fn pressure_point(pressure_pa: f64) -> QueryPoint {
    source::point()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), pressure_pa)
        .unwrap()
}

#[test]
fn pressure_and_composition_data_drive_analytic_current_and_both_body_heats() {
    let (a, ca) = material_with_conductivity(1.0, 1.0, 0.0, 0.0);
    let (b, cb) = material_with_conductivity(3.0, 2.0, 0.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let context = point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(
        &cards,
        [curves[0].curve(), curves[1].curve()],
        &context,
        &key,
    );
    let resistors = resistors();
    let links = links();
    let bases = bases();
    let bulk = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &links,
    };
    let thermal = initial(&curves, 999.0, None);
    for (axis, kind, lo, hi) in [
        ("pressure", QuantityKind::Pressure, 100_000.0, 200_000.0),
        (
            "oxide_fraction",
            QuantityKind::Composition(CompositionBasis::MassFraction),
            0.0,
            1.0,
        ),
    ] {
        let q = quantity(kind);
        let mut claims = ClaimSet::new();
        let pin = claims
            .insert_claim(claim(
                EP,
                ED,
                PropertyValue::Curve {
                    abscissa: axis.into(),
                    abscissa_dims: q.dims(),
                    knots: vec![(lo, 2.0), (hi, 1.0)],
                    dims: ED,
                },
                validity().with_quantity(axis, q, lo, hi),
            ))
            .unwrap();
        let card = interface(&cards, claims);
        let mut charges = Vec::new();
        for (coordinate, contact_r) in [(lo, 2.0), (hi, 1.0)] {
            let interface_point = source::point().with_quantity(axis, q, coordinate).unwrap();
            let bindings = [contact(&cards, &card, &interface_point, 0.25)];
            let network = ElectricalContactNetwork {
                bulk: &bulk,
                contacts: &bindings,
                electrical_bases: &bases,
            };
            let electrical = circuit(&thermal, &cards, &curves)
                .rebind_contacts_candidate(
                    network.resolve(&thermal, || false).unwrap().contacts,
                    &[],
                )
                .unwrap();
            let next = network
                .step_candidate(
                    &thermal,
                    0.25,
                    &electrical,
                    &[],
                    config(),
                    &[0.0; 2],
                    || false,
                )
                .unwrap();
            // Independent midpoint discharge of C=1 F through 1+3+R_contact ohm.
            let total_r = 4.0 + contact_r;
            let alpha = 0.25 / (2.0 * total_r);
            let q1 = 10.0 * (1.0 - alpha) / (1.0 + alpha);
            let current = (10.0 + q1) / (2.0 * total_r);
            let unit_heat = 0.25 * current * current;
            let loss = &next.report().electrical.contacts[0];
            near(next.circuit().charges_c()[0], q1);
            near(loss.current_a, current);
            near(loss.heat_j, unit_heat * contact_r);
            near(
                next.thermal().temperature_k(CellId(30)).unwrap(),
                300.0 + unit_heat * (1.0 + 0.25 * contact_r),
            );
            near(
                next.thermal().temperature_k(CellId(10)).unwrap(),
                300.0 + unit_heat * (3.0 + 0.75 * contact_r) / 2.0,
            );
            near(
                next.thermal().energy_j() - thermal.energy_j()
                    + next.report().electrical.stored_after_j
                    - 50.0,
                0.0,
            );
            near(next.report().deposited_joule_j, unit_heat * total_r);
            let resolved = loss.contact.source().resolved();
            assert_eq!(
                resolved.property(EP).unwrap().answer().receipt.selected,
                pin
            );
            assert_eq!(resolved.card_identity(), card.content_hash());
            assert_eq!(resolved.query().axes()[axis], coordinate);
            assert_eq!(resolved.query().axis_quantities()[axis], q);
            charges.push(q1);
        }
        assert!(
            charges[1] < charges[0],
            "only external source data changed the discharge"
        );
    }
}

fn surface(axis: &str, dims: Dims, law: impl Fn(f64, f64) -> f64, hole: bool) -> PropertyValue {
    PropertyValue::Surface2D {
        axes: [axis.into(), "pressure".into()],
        axis_dims: [Temperature::DIMS, quantity(QuantityKind::Pressure).dims()],
        samples: [
            [200.0, 100_000.0],
            [500.0, 100_000.0],
            [200.0, 200_000.0],
            [500.0, 200_000.0],
        ]
        .map(|p| (p, law(p[0], p[1])))
        .to_vec(),
        triangles: if hole {
            vec![[0, 1, 2]]
        } else {
            vec![[0, 1, 2], [1, 3, 2]]
        },
        dims,
    }
}

fn thermal_contact<'a>(
    cards: &'a [MaterialCard; 2],
    curves: [&'a EquilibriumEnthalpyPhaseCurve; 2],
    card: &'a InterfaceSystemCard,
    context: &'a QueryPoint,
    pin: ClaimId,
) -> MaterialContactLink<'a> {
    MaterialContactLink {
        link: 7,
        sides: [0, 1].map(|i| ContactCell {
            cell: DESTINATIONS[i].cell,
            material_id: MaterialId(i as u64),
            curve: curves[i],
        }),
        point: ThermalInterfacePoint {
            interface: card,
            materials: [&cards[0], &cards[1]],
            context,
            temperature_axes: ["Ta", "Tb"],
            liquid_fraction_axes: [Some("melt_a"), Some("melt_b")],
            mean_temperature_axis: None,
        },
        selection: ClaimSelection::Pinned(pin),
    }
}

#[test]
fn one_interface_pressure_changes_both_live_transport_surfaces_and_closed_history() {
    let (a, ca) = material_with_conductivity(1.0, 1.0, 0.005, 1.0);
    let (b, cb) = material_with_conductivity(3.0, 2.0, 0.0, 2.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let context = point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(
        &cards,
        [curves[0].curve(), curves[1].curve()],
        &context,
        &key,
    );
    let resistors = resistors();
    let links = links();
    let bases = bases();
    let bulk = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &links,
    };
    let electrical_r =
        |t: f64, p: f64| 2.0 + 0.01 * (t - 300.0) - 0.5 * (p - 100_000.0) / 100_000.0;
    let thermal_r = |t: f64, p: f64| 1.0 + 0.002 * (t - 300.0) - 0.5 * (p - 100_000.0) / 100_000.0;
    let domain = validity().with_quantity(
        "pressure",
        quantity(QuantityKind::Pressure),
        100_000.0,
        200_000.0,
    );
    let mut claims = ClaimSet::new();
    let epin = claims
        .insert_claim(claim(
            EP,
            ED,
            surface("Ta", ED, electrical_r, false),
            domain.clone(),
        ))
        .unwrap();
    let tpin = claims
        .insert_claim(claim(RP, RD, surface("Tb", RD, thermal_r, false), domain))
        .unwrap();
    let card = interface(&cards, claims);
    let run = |pressure_pa| {
        let interface_point = pressure_point(pressure_pa);
        let bindings = [contact(&cards, &card, &interface_point, 0.25)];
        let thermal_bindings = [thermal_contact(
            &cards,
            [curves[0].curve(), curves[1].curve()],
            &card,
            &interface_point,
            tpin,
        )];
        let network = ElectricalContactNetwork {
            bulk: &bulk,
            contacts: &bindings,
            electrical_bases: &bases,
        };
        let mut thermal = initial(&curves, 999.0, None);
        let prepared = network.prepare(&thermal, || false).unwrap();
        let mut electrical = circuit(&thermal, &cards, &curves)
            .rebind_contacts_candidate(prepared.resolve(&thermal, || false).unwrap().contacts, &[])
            .unwrap();
        let total = thermal.energy_j() + 50.0;
        for i in 1..=20 {
            let target = f64::from(i) * 0.1;
            let dt = target - thermal.time_s();
            let next = contacts::step_candidate(
                &thermal,
                target,
                input(&electrical),
                &routes(0.25),
                config(),
                &[0.0; 2],
                |state| {
                    let mut coefficients = prepared
                        .resolve(state, || false)
                        .map_err(|e| e.to_string())?;
                    let samples = apply_material_contacts(
                        state,
                        &thermal_bindings,
                        &mut coefficients.bulk.conductances_w_k,
                        || false,
                    )
                    .map_err(|e| e.to_string())?;
                    Ok::<_, String>(contacts::ContactThermalCoefficients {
                        bulk: CircuitThermalCoefficients {
                            resistors: coefficients.bulk.resistors,
                            sources_w: coefficients.bulk.sources_w,
                            conductances_w_k: coefficients.bulk.conductances_w_k,
                            source: (coefficients.bulk.source, samples),
                        },
                        contacts: coefficients.contacts,
                    })
                },
                || false,
            )
            .unwrap();
            let report = next.report();
            let contact = &report.source.material.1[0];
            let ta = next.thermal().temperature_k(CellId(30)).unwrap();
            let tb = next.thermal().temperature_k(CellId(10)).unwrap();
            let loss = &report.electrical.contacts[0];
            near(loss.contact.resistance_ohm(), electrical_r(ta, pressure_pa));
            near(
                contact.interface.value_m2_k_per_w(),
                thermal_r(tb, pressure_pa),
            );
            near(
                contact.conductance_w_k,
                1.0 / (0.75 + thermal_r(tb, pressure_pa)),
            );
            assert_eq!(contact.interface.receipt().selected, tpin);
            assert_eq!(
                loss.contact
                    .source()
                    .resolved()
                    .property(EP)
                    .unwrap()
                    .answer()
                    .receipt
                    .selected,
                epin
            );
            let exchange = dt * contact.conductance_w_k * (tb - ta);
            let da = next.thermal().cells()[&CellId(30)].specific_h_j_kg
                - thermal.cells()[&CellId(30)].specific_h_j_kg;
            let db = next.thermal().cells()[&CellId(10)].specific_h_j_kg
                - thermal.cells()[&CellId(10)].specific_h_j_kg;
            near(
                da,
                report.electrical.resistors[0].heat_j + 0.25 * loss.heat_j + exchange,
            );
            near(
                db,
                report.electrical.resistors[1].heat_j + 0.75 * loss.heat_j - exchange,
            );
            near(
                next.thermal().energy_j() + report.electrical.stored_after_j,
                total,
            );
            near(report.deposited_joule_j, report.electrical.heat_j);
            let retained = loss.contact.source().resolved().query();
            near(retained.axes()["Ta"], ta);
            near(retained.axes()["Tb"], tb);
            assert_eq!(retained.axes()["pressure"], pressure_pa);
            let (t, c, _) = next.into_parts();
            thermal = t;
            electrical = c;
        }
        let ta = thermal.temperature_k(CellId(30)).unwrap();
        let tb = thermal.temperature_k(CellId(10)).unwrap();
        assert!(
            (ta - tb).abs() > 0.01,
            "the two body states remain independent"
        );
        (ta, tb, electrical.charges_c()[0])
    };
    let low = run(100_000.0);
    let high = run(200_000.0);
    assert!((high.0 - low.0).abs() > 0.01 && (high.1 - low.1).abs() > 0.01);
    assert!(
        high.2 < low.2,
        "the same authored source responds to the pressure data"
    );
}

#[test]
fn unsupported_response_axes_and_surface_holes_leave_physical_state_and_coefficients_unchanged() {
    let (a, ca) = material(1.0, 1.0, 0.0);
    let (b, cb) = material(3.0, 2.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let context = point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(
        &cards,
        [curves[0].curve(), curves[1].curve()],
        &context,
        &key,
    );
    let resistors = resistors();
    let links = links();
    let bases = bases();
    let bulk = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &links,
    };
    let thermal = initial(&curves, 999.0, None);
    let saved = thermal.cells().clone();
    let mut conductances = bulk.resolve(&thermal, || false).unwrap().conductances_w_k;
    let original = conductances.clone();
    // Same dimensions do not turn an evolving mechanical stress into fixed pressure.
    for q in [
        quantity(QuantityKind::Stress),
        QuantitySpec::dimensional(quantity(QuantityKind::Pressure).dims()),
    ] {
        let mut claims = ClaimSet::new();
        let value = |dims| PropertyValue::Curve {
            abscissa: "pressure".into(),
            abscissa_dims: q.dims(),
            knots: vec![(100_000.0, 2.0), (200_000.0, 1.0)],
            dims,
        };
        let domain = validity().with_quantity("pressure", q, 100_000.0, 200_000.0);
        claims
            .insert_claim(claim(EP, ED, value(ED), domain.clone()))
            .unwrap();
        let tpin = claims
            .insert_claim(claim(RP, RD, value(RD), domain))
            .unwrap();
        let card = interface(&cards, claims);
        let query = source::point()
            .with_quantity("pressure", q, 150_000.0)
            .unwrap();
        let bindings = [contact(&cards, &card, &query, 0.25)];
        let network = ElectricalContactNetwork {
            bulk: &bulk,
            contacts: &bindings,
            electrical_bases: &bases,
        };
        assert!(network.resolve(&thermal, || false).is_err());
        let thermal_bindings = [thermal_contact(
            &cards,
            [curves[0].curve(), curves[1].curve()],
            &card,
            &query,
            tpin,
        )];
        assert!(
            apply_material_contacts(&thermal, &thermal_bindings, &mut conductances, || false)
                .is_err()
        );
        assert_eq!(conductances, original);
    }
    let mut claims = ClaimSet::new();
    let tpin = claims
        .insert_claim(claim(
            RP,
            RD,
            surface("Tb", RD, |_, _| 1.0, true),
            validity().with_quantity(
                "pressure",
                quantity(QuantityKind::Pressure),
                100_000.0,
                200_000.0,
            ),
        ))
        .unwrap();
    let card = interface(&cards, claims);
    let query = pressure_point(200_000.0);
    let thermal_bindings = [thermal_contact(
        &cards,
        [curves[0].curve(), curves[1].curve()],
        &card,
        &query,
        tpin,
    )];
    // (300 K, 200 kPa) is inside the bounding box, outside the declared triangle.
    assert!(
        apply_material_contacts(&thermal, &thermal_bindings, &mut conductances, || false).is_err()
    );
    assert_eq!(conductances, original);
    assert_eq!(thermal.cells(), &saved);
}
