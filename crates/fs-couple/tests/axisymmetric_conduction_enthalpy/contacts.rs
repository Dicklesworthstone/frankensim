//! Swept contact work, authored heat fractions and original source feedback.
//! Synthetic coefficients below are independent manufactured-model fixtures.
use super::*;
use fs_couple::circuit_enthalpy::conduction::contacts::ConductionContactThermal;
use fs_em::conduction::ConductionContact;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, PropertyClaim, PropertyValue, Provenance,
    SurfaceSpec, SystemContext, UncertaintyModel,
};
use fs_material::{
    electrical_interface::{
        ELECTRICAL_CONTACT_RESISTANCE_DIMS, ELECTRICAL_CONTACT_RESISTANCE_PROPERTY,
        resolve_electrical_interface_state,
    },
    thermal::interface::ThermalInterfacePoint,
};
use fs_qty::{
    QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

fn absolute_temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn pressure_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static))
}

pub(super) fn interface(
    cards: &[MaterialCard; 2],
    r0: f64,
    beta: f64,
    upper: f64,
) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                ELECTRICAL_CONTACT_RESISTANCE_PROPERTY,
                ELECTRICAL_CONTACT_RESISTANCE_DIMS,
            ),
            value: PropertyValue::Curve {
                abscissa: "Tm".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![
                    (200.0, r0 * (1.0 - 100.0 * beta)),
                    (upper, r0 * (1.0 + (upper - 300.0) * beta)),
                ],
                dims: ELECTRICAL_CONTACT_RESISTANCE_DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity("Tm", absolute_temperature(), 200.0, upper)
                .with_quantity("pressure", pressure_quantity(), 1e5, 1e5),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured electrothermal contact, no experimental claim".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: cards[0].id().clone(),
            texture_frame: "first surface".into(),
        },
        SurfaceSpec {
            material: cards[1].id().clone(),
            texture_frame: "second surface".into(),
        },
        SystemContext {
            medium: "synthetic".into(),
            third_body: None,
            environment: "fixed pressure".into(),
            history: "fixed manufactured interface".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}

pub(super) fn interface_context(t: f64) -> QueryPoint {
    let mut point = QueryPoint::new();
    for axis in ["Ta", "Tb", "Tm"] {
        point = point
            .with_quantity(axis, absolute_temperature(), t)
            .unwrap();
    }
    point
        .with_quantity("pressure", pressure_quantity(), 1e5)
        .unwrap()
}

pub(super) fn binding<'a>(
    cards: &'a [MaterialCard; 2],
    card: &'a InterfaceSystemCard,
    context: &'a QueryPoint,
    fraction: f64,
) -> ConductionContactThermal<'a> {
    ConductionContactThermal {
        contact: 0,
        point: ThermalInterfacePoint {
            interface: card,
            materials: [&cards[0], &cards[1]],
            context,
            temperature_axes: ["Ta", "Tb"],
            liquid_fraction_axes: [None, None],
            mean_temperature_axis: Some("Tm"),
        },
        heat_fraction_to_first: fraction,
    }
}

fn contact_mesh(angle: f64) -> TriComplex2 {
    TriComplex2::from_indexed_triangles(
        "two cylinder layers with a matching swept disk contact",
        vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [1.0, 2.0],
            [0.0, 2.0],
        ],
        vec![[0, 1, 2], [0, 2, 3], [4, 5, 6], [4, 6, 7]],
        Metric2::axisymmetric(angle).unwrap(),
    )
    .unwrap()
}
fn specimen(
    cards: &[MaterialCard; 2],
    interface: &InterfaceSystemCard,
    angle: f64,
) -> AxisymmetricConduction {
    let sources = [0, 0, 1, 1]
        .map(|i| {
            resolve_conductor_state_with_basis(
                &cards[i],
                &point(300.0),
                MaterialPropertySelection::SingleClaimOnly,
                ElectricalPropertyBasis::Resistivity,
            )
            .unwrap()
        })
        .to_vec();
    let source = resolve_electrical_interface_state(
        interface,
        &interface_context(300.0),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    AxisymmetricConduction::new_with_contacts(
        contact_mesh(angle),
        sources,
        vec![ConductionContact {
            first_vertices: [2, 3],
            second_vertices: [5, 4],
            source,
        }],
        limits(),
        &mut || false,
    )
    .unwrap()
}

#[test]
fn swept_interface_feedback_closes_both_body_balances_for_current_and_voltage() {
    let alpha = 0.002;
    let beta = 0.003;
    let fraction = 0.25;
    let dt = 0.5;
    let j2 = 4.0;
    for angle in [0.7, 2.0, core::f64::consts::TAU] {
        let area = angle / 2.0;
        for r0 in [0.5, 1.5] {
            let (a, ca) = material(1.0, 1.0, alpha);
            let (b, cb) = material(3.0, 2.0, 0.0);
            let cards = [a, b];
            let curves = [ca, cb];
            let thermal = initial(&curves, angle);
            let interface = interface(&cards, r0, beta, 1000.0);
            let original = specimen(&cards, &interface, angle);
            let context = point(300.0);
            let surface_context = interface_context(300.0);
            let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
            let cells = cells(&cards, &curves, &context, &key);
            let layout = network(&cells);
            let contacts = [binding(&cards, &interface, &surface_context, fraction)];
            // Divide both finite-body equations by actual swept area A.
            // C1=A, C2=2*A, G=A, bulk P=j^2*A*rho; contact P=j^2*A*R''.
            let slope = r0 * beta / 2.0;
            let a11 = 1.0 + dt - dt * j2 * (alpha + fraction * slope);
            let a12 = -dt - dt * j2 * fraction * slope;
            let a21 = -dt - dt * j2 * (1.0 - fraction) * slope;
            let a22 = 2.0 + dt - dt * j2 * (1.0 - fraction) * slope;
            let rhs1 = dt * j2 * (1.0 + fraction * r0);
            let rhs2 = dt * j2 * (3.0 + (1.0 - fraction) * r0);
            let det = a11 * a22 - a12 * a21;
            let d1 = (rhs1 * a22 - a12 * rhs2) / det;
            let d2 = (a11 * rhs2 - a21 * rhs1) / det;
            let ra = 1.0 + alpha * d1;
            let rc = r0 * (1.0 + beta * (d1 + d2) / 2.0);
            let v = 2.0 * (ra + rc + 3.0);
            let voltages = [(0, 0.0), (1, 0.0), (6, v), (7, v)];
            for drive in [
                ConductionDrive::Current {
                    positive: &[6, 7],
                    negative: &[0, 1],
                    current_a: 2.0 * area,
                },
                ConductionDrive::Voltage(&voltages),
            ] {
                let next = layout
                    .step_axisymmetric_contact_candidate(
                        &thermal,
                        dt,
                        &original,
                        drive,
                        &contacts,
                        limits(),
                        field_controls(),
                        config(),
                        &[0.0; 2],
                        || false,
                    )
                    .unwrap();
                near(
                    next.thermal().temperature_k(CellId(30)).unwrap(),
                    300.0 + d1,
                );
                near(
                    next.thermal().temperature_k(CellId(10)).unwrap(),
                    300.0 + d2,
                );
                let step = &next.report().electrical;
                let heat = &step.contacts[0];
                near(step.voltage_v[6], v);
                near(heat.response.current_a, -2.0 * area);
                for jump in heat.response.voltage_jump_v {
                    near(jump, -2.0 * rc);
                }
                for density in heat.response.current_density_a_m2 {
                    near(density, -2.0);
                }
                near(heat.response.power_w, j2 * area * rc);
                assert_eq!(heat.cells, [CellId(30), CellId(10)]);
                near(heat.heat_fraction_to_first, fraction);
                for (triangle, power) in step.cell_power_w.iter().enumerate() {
                    let rho = if triangle < 2 { ra } else { 3.0 };
                    near(
                        *power,
                        j2 * rho * original.mesh().face_measure(triangle).unwrap(),
                    );
                    near(step.current_density_a_m2[triangle][0], 0.0);
                    near(step.current_density_a_m2[triangle][1], -2.0);
                }
                near(step.cell_power_w.iter().sum(), j2 * area * (ra + 3.0));
                near(step.heat_j, dt * j2 * area * (ra + rc + 3.0));
                near(step.supplied_work_j, step.heat_j);
                near(next.report().deposited_joule_j, step.heat_j);
                near(
                    next.thermal().energy_j() - thermal.energy_j(),
                    area * (d1 + 2.0 * d2),
                );
                assert!(
                    next.report().total_energy_residual_j.abs() < config().total_energy_tolerance_j
                );
                let refreshed = &next.problem().contacts()[0].source;
                near(refreshed.resistance_ohm_m2(), rc);
                for (axis, expected) in [
                    ("Ta", 300.0 + d1),
                    ("Tb", 300.0 + d2),
                    ("Tm", 300.0 + (d1 + d2) / 2.0),
                    ("pressure", 1e5),
                ] {
                    near(
                        refreshed.resolved().state_coordinate(axis).unwrap(),
                        expected,
                    );
                }
                assert_eq!(
                    refreshed.resolved().card_identity(),
                    original.contacts()[0].source.resolved().card_identity()
                );
                assert_eq!(
                    refreshed.resolved().properties()[0]
                        .answer()
                        .receipt
                        .selected,
                    original.contacts()[0].source.resolved().properties()[0]
                        .answer()
                        .receipt
                        .selected
                );
                for (triangle, material) in next.problem().conductors().iter().enumerate() {
                    near(
                        material.material().query().axes()["T"],
                        300.0 + if triangle < 2 { d1 } else { d2 },
                    );
                }
                let continuing = layout
                    .step_axisymmetric_contact_candidate(
                        next.thermal(),
                        2.0 * dt,
                        next.problem(),
                        drive,
                        &contacts,
                        limits(),
                        field_controls(),
                        config(),
                        &[0.0; 2],
                        || false,
                    )
                    .unwrap();
                assert!(
                    continuing.problem().contacts()[0]
                        .source
                        .resistance_ohm_m2()
                        > rc
                );
            }
            assert_eq!(thermal.time_s(), 0.0);
            near(original.contacts()[0].source.resistance_ohm_m2(), r0);
            assert!(
                layout
                    .step_axisymmetric_candidate(
                        &thermal,
                        dt,
                        &original,
                        ConductionDrive::Current {
                            positive: &[6, 7],
                            negative: &[0, 1],
                            current_a: 2.0 * area
                        },
                        limits(),
                        field_controls(),
                        config(),
                        &[0.0; 2],
                        || false,
                    )
                    .is_err()
            );
        }
    }
}

#[test]
fn contact_source_partition_support_and_late_cancellation_keep_inputs_reusable() {
    let angle = 1.3;
    let area = angle / 2.0;
    let (a, ca) = material(1.0, 1.0, 0.0);
    let (b, cb) = material(3.0, 2.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let thermal = initial(&curves, angle);
    let interface = interface(&cards, 0.5, 0.002, 1000.0);
    let original = specimen(&cards, &interface, angle);
    let context = point(300.0);
    let surface_context = interface_context(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &context, &key);
    let layout = network(&cells);
    let contact = binding(&cards, &interface, &surface_context, 0.25);
    let run = |problem: &AxisymmetricConduction,
               bindings: &[ConductionContactThermal<'_>],
               cancel: &mut dyn FnMut() -> bool| {
        layout.step_axisymmetric_contact_candidate(
            &thermal,
            0.5,
            problem,
            ConductionDrive::Current {
                positive: &[6, 7],
                negative: &[0, 1],
                current_a: 2.0 * area,
            },
            bindings,
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            cancel,
        )
    };
    let energy = thermal.energy_j();
    let source = original.contacts()[0].source.clone();
    let mut count = 0;
    let baseline = run(&original, &[contact], &mut || {
        count += 1;
        false
    })
    .unwrap();
    assert!(count > 100);
    let mut polls = 0;
    assert!(matches!(
        run(&original, &[contact], &mut || {
            polls += 1;
            polls >= count - 1
        }),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    assert!(run(&original, &[], &mut || false).is_err());
    let mut invalid_fraction = contact;
    invalid_fraction.heat_fraction_to_first = 1.1;
    assert!(run(&original, &[invalid_fraction], &mut || false).is_err());
    let pressure_context = surface_context.clone().with("pressure", 100001.0).unwrap();
    let mut wrong_pressure = contact;
    wrong_pressure.point.context = &pressure_context;
    assert!(run(&original, &[wrong_pressure], &mut || false).is_err());
    let mut reversed_materials = contact;
    reversed_materials.point.materials = [&cards[1], &cards[0]];
    assert!(run(&original, &[reversed_materials], &mut || false).is_err());
    let stale_context = surface_context.clone().with("Ta", 300.005).unwrap();
    let stale = AxisymmetricConduction::new_with_contacts(
        original.mesh().clone(),
        original.conductors().to_vec(),
        vec![ConductionContact {
            source: source.at_query_point(&interface, &stale_context).unwrap(),
            ..original.contacts()[0].clone()
        }],
        limits(),
        &mut || false,
    )
    .unwrap();
    near(
        stale.contacts()[0].source.resistance_ohm_m2(),
        source.resistance_ohm_m2(),
    );
    assert!(run(&stale, &[contact], &mut || false).is_err());
    let narrow = self::interface(&cards, 0.5, 0.002, 300.01);
    let unsupported = specimen(&cards, &narrow, angle);
    let narrow_binding = binding(&cards, &narrow, &surface_context, 0.25);
    assert!(run(&unsupported, &[narrow_binding], &mut || false).is_err());
    assert_eq!(thermal.energy_j().to_bits(), energy.to_bits());
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(original.contacts()[0].source, source);
    let retry = run(&original, &[contact], &mut || false).unwrap();
    assert_eq!(
        retry.report().electrical.voltage_v,
        baseline.report().electrical.voltage_v
    );
    assert_eq!(
        retry.report().electrical.contacts,
        baseline.report().electrical.contacts
    );
    assert_eq!(
        retry.thermal().energy_j().to_bits(),
        baseline.thermal().energy_j().to_bits()
    );
}

#[test]
fn swept_contact_refreshes_both_phase_coordinates_on_a_latent_plateau() {
    use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
    use fs_qty::Dims;
    use fs_thermal::enthalpy::{
        CellGeometry, DepositedCell, DepositionBudget, DepositionSimulation, Event, StepConfig,
        ThermalLaw, TimedEvent,
    };
    let scalar = QuantitySpec::dimensional(Dims::NONE);
    let bulk_point = QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 310.0)
        .unwrap()
        .with_quantity("f", scalar, 0.2)
        .unwrap();
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                "electrical_resistivity",
                ElectricalPropertyBasis::Resistivity.dimensions(),
            ),
            value: PropertyValue::Curve {
                abscissa: "f".into(),
                abscissa_dims: Dims::NONE,
                knots: vec![(0.0, 1.0), (1.0, 1.0)],
                dims: ElectricalPropertyBasis::Resistivity.dimensions(),
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    300.0,
                    350.0,
                )
                .with_quantity("f", scalar, 0.0, 1.0),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic constant-density phase conductor".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        fs_matdb::MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "confined solid-liquid".into(),
            process: "authored test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let cards = [card.clone(), card];
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(
        cards[0].content_hash(),
        [
            (0.0, 300.0, 0.0),
            (20.0, 310.0, 0.0),
            (120.0, 310.0, 1.0),
            (200.0, 350.0, 1.0),
        ]
        .map(|(h, t, f)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: f,
            bulk_density_kg_m3: 1.0,
        })
        .to_vec(),
    )
    .unwrap();
    let events = [CellId(30), CellId(10)]
        .into_iter()
        .enumerate()
        .map(|(i, id)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id,
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MaterialId(i as u64),
                mass_kg: 1.0,
                specific_h_j_kg: 40.0,
                geometry: CellGeometry {
                    start_m: [0.0, 0.0, i as f64],
                    end_m: [0.0, 0.0, i as f64 + 1.0],
                    volume_m3: 1.0, // radius=1, angle=2, axial length=1.
                    surface_area_m2: 4.0,
                },
                boundary: common::boundary(),
            }),
        })
        .collect();
    let thermal = DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), &curve as &dyn ThermalLaw),
            (MaterialId(1), &curve as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 2,
            max_links: 0,
            max_events: 2,
            max_substeps: 100,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.1,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                ELECTRICAL_CONTACT_RESISTANCE_PROPERTY,
                ELECTRICAL_CONTACT_RESISTANCE_DIMS,
            ),
            value: PropertyValue::Curve {
                abscissa: "fa".into(),
                abscissa_dims: Dims::NONE,
                knots: vec![(0.0, 1.0), (1.0, 3.0)],
                dims: ELECTRICAL_CONTACT_RESISTANCE_DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with("Ta", 300.0, 350.0)
                .with("Tb", 300.0, 350.0)
                .with_quantity("fa", scalar, 0.0, 1.0)
                .with_quantity("fb", scalar, 0.0, 1.0),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic phase-dependent interface".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let interface = InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: cards[0].id().clone(),
            texture_frame: "A".into(),
        },
        SurfaceSpec {
            material: cards[1].id().clone(),
            texture_frame: "B".into(),
        },
        SystemContext {
            medium: "confined material interface".into(),
            third_body: None,
            environment: "fixed".into(),
            history: "fixed".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    let ip = interface_context(310.0)
        .with_quantity("fa", scalar, 0.2)
        .unwrap()
        .with_quantity("fb", scalar, 0.2)
        .unwrap();
    let conductor = resolve_conductor_state_with_basis(
        &cards[0],
        &bulk_point,
        MaterialPropertySelection::SingleClaimOnly,
        ElectricalPropertyBasis::Resistivity,
    )
    .unwrap()
    .with_liquid_mass_fraction_axis("f")
    .unwrap();
    let contacts = vec![ConductionContact {
        first_vertices: [2, 3],
        second_vertices: [5, 4],
        source: resolve_electrical_interface_state(
            &interface,
            &ip,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap(),
    }];
    // The original fixed-external API must still refuse unowned phase response.
    assert!(
        AxisymmetricConduction::new_with_contacts(
            contact_mesh(2.0),
            vec![conductor.clone(); 4],
            contacts.clone(),
            limits(),
            &mut || false
        )
        .is_err()
    );
    let problem = AxisymmetricConduction::new_with_phase_contacts(
        contact_mesh(2.0),
        vec![conductor; 4],
        contacts,
        vec![[Some("fa".into()), Some("fb".into())]],
        limits(),
        &mut || false,
    )
    .unwrap();
    let cells = [CellId(30), CellId(10)]
        .into_iter()
        .enumerate()
        .map(|(i, cell)| MaterialCell {
            cell,
            material_id: MaterialId(i as u64),
            curve: &curve,
            point: ThermalMaterialPoint {
                material: &cards[i],
                context: &bulk_point,
                temperature_axis: "T",
                liquid_fraction_axis: Some("f"),
            },
            conductivity: None,
        })
        .collect::<Vec<_>>();
    let network = ConductionThermalNetwork {
        cells: &cells,
        links: &[],
        regions: &REGIONS,
        volume_tolerance_m3: 1e-12,
    };
    let mut contact_binding = binding(&cards, &interface, &ip, 0.25);
    contact_binding.point.liquid_fraction_axes = [Some("fa"), Some("fb")];
    let next = network
        .step_axisymmetric_contact_candidate(
            &thermal,
            0.5,
            &problem,
            ConductionDrive::Current {
                positive: &[6, 7],
                negative: &[0, 1],
                current_a: 2.0,
            },
            &[contact_binding],
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    // Independent latent balance: 100*(fa-.2)=2*(1+.25*(1+2*fa)).
    let fa = 22.5 / 99.0;
    let fb = 0.2 + (3.5 + 3.0 * fa) / 100.0;
    let resistance = 1.0 + 2.0 * fa;
    near(next.thermal().temperature_k(CellId(30)).unwrap(), 310.0);
    near(next.thermal().temperature_k(CellId(10)).unwrap(), 310.0);
    let source = next.problem().contacts()[0].source.resolved();
    near(source.state_coordinate("fa").unwrap(), fa);
    near(source.state_coordinate("fb").unwrap(), fb);
    near(
        next.report().electrical.contacts[0].response.power_w,
        4.0 * resistance,
    );
    near(
        next.report().electrical.voltage_v[6],
        2.0 * (2.0 + resistance),
    );
    near(
        next.thermal().energy_j() - thermal.energy_j(),
        100.0 * (fa + fb - 0.4),
    );
    near(next.report().deposited_joule_j, 2.0 * (2.0 + resistance));
    let mut missing = contact_binding;
    missing.point.liquid_fraction_axes = [Some("fa"), None];
    assert!(
        network
            .step_axisymmetric_contact_candidate(
                &thermal,
                0.5,
                &problem,
                ConductionDrive::Current {
                    positive: &[6, 7],
                    negative: &[0, 1],
                    current_a: 2.0
                },
                &[missing],
                limits(),
                field_controls(),
                config(),
                &[0.0; 2],
                || false
            )
            .is_err()
    );

    assert_eq!(thermal.time_s(), 0.0);
    near(thermal.energy_j(), 80.0);
}
