//! G1/G3/G4: spatial contact heat, source refresh and native circuit feedback.
//! Manufactured laws test implementation; they are not measured material data.

use super::*;
use fs_couple::circuit_enthalpy::conduction::{
    circuit::{ConductionCircuit, ConductionCircuitTerminal},
    contacts::ConductionContactThermal,
};
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
use fs_phs::circuit::{Branch, CircuitGraph};
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

pub(super) fn interface(cards: &[MaterialCard; 2], r0: f64, beta: f64, upper: f64) -> InterfaceSystemCard {
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

pub(super) fn contact_mesh() -> TriComplex2 {
    TriComplex2::from_indexed_triangles(
        "two original bodies with matching finite-resistance traces",
        vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            [1.0, 0.0],
            [2.0, 0.0],
            [2.0, 1.0],
            [1.0, 1.0],
        ],
        vec![[0, 1, 2], [0, 2, 3], [4, 5, 6], [4, 6, 7]],
        Metric2::planar(1.0).unwrap(),
    )
    .unwrap()
}

fn specimen(cards: &[MaterialCard; 2], interface: &InterfaceSystemCard) -> PlanarConduction {
    let mesh = contact_mesh();
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
    PlanarConduction::new_with_contacts(
        mesh,
        sources,
        vec![ConductionContact {
            first_vertices: [1, 2],
            second_vertices: [4, 7],
            source,
        }],
        limits(),
        &mut || false,
    )
    .unwrap()
}

#[test]
fn g1_two_sided_contact_heat_and_temperature_feedback_match_independent_balance() {
    let alpha = 0.002;
    let beta = 0.003;
    let fraction = 0.25;
    let dt = 0.5;
    let i2 = 4.0;
    let mut work = Vec::new();
    for r0 in [0.5, 1.5] {
        let (a, ca) = material(1.0, 1.0, alpha);
        let (b, cb) = material(3.0, 2.0, 0.0);
        let cards = [a, b];
        let curves = [ca, cb];
        let thermal = initial(&curves, 999.0, None);
        let interface = interface(&cards, r0, beta, 1000.0);
        let original = specimen(&cards, &interface);
        let context = point(300.0);
        let surface_context = interface_context(300.0);
        // Re-querying pins the same physical claim; its policy receipt changes.
        // Admission must check source/state identity without requiring that
        // an initial SingleClaimOnly receipt already be a pinned receipt.
        let pinned = original.contacts()[0]
            .source
            .at_query_point(&interface, &surface_context)
            .unwrap();
        assert_ne!(
            pinned.resolved().identity(),
            original.contacts()[0].source.resolved().identity()
        );
        assert_eq!(
            pinned.resolved().properties()[0].answer().receipt.selected,
            original.contacts()[0].source.resolved().properties()[0]
                .answer()
                .receipt
                .selected
        );
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = cells(&cards, &curves, &context, &key);
        let layout = network(&cells);
        let contacts = [binding(&cards, &interface, &surface_context, fraction)];
        let drive = ConductionDrive::Current {
            positive: &[5, 6],
            negative: &[0, 3],
            current_a: 2.0,
        };
        let next = layout
            .step_contact_candidate(
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
        // Independent two-body backward-Euler equations: each native body has
        // R=rho*length/area, the interface has R=R''/area, G=1 W/K.
        let half_slope = r0 * beta / 2.0;
        let a11 = 1.0 + dt - dt * i2 * (alpha + fraction * half_slope);
        let a12 = -dt - dt * i2 * fraction * half_slope;
        let a21 = -dt - dt * i2 * (1.0 - fraction) * half_slope;
        let a22 = 2.0 + dt - dt * i2 * (1.0 - fraction) * half_slope;
        let rhs1 = dt * i2 * (1.0 + fraction * r0);
        let rhs2 = dt * i2 * (3.0 + (1.0 - fraction) * r0);
        let det = a11 * a22 - a12 * a21;
        let d1 = (rhs1 * a22 - a12 * rhs2) / det;
        let d2 = (a11 * rhs2 - a21 * rhs1) / det;
        near(
            next.thermal().temperature_k(CellId(30)).unwrap(),
            300.0 + d1,
        );
        near(
            next.thermal().temperature_k(CellId(10)).unwrap(),
            300.0 + d2,
        );
        let ra = 1.0 + alpha * d1;
        let rc = r0 * (1.0 + beta * (d1 + d2) / 2.0);
        let step = &next.report().electrical;
        near(step.voltage_v[5], 2.0 * (ra + rc + 3.0));
        near(step.contacts[0].response.current_a, -2.0);
        near(step.contacts[0].response.power_w, i2 * rc);
        assert_eq!(step.contacts[0].cells, [CellId(30), CellId(10)]);
        near(step.contacts[0].heat_fraction_to_first, fraction);
        near(step.cell_power_w.iter().sum(), i2 * (ra + 3.0));
        near(step.heat_j, dt * i2 * (ra + rc + 3.0));
        near(step.supplied_work_j, step.heat_j);
        near(next.report().deposited_joule_j, step.heat_j);
        near(
            next.thermal().energy_j() - thermal.energy_j(),
            d1 + 2.0 * d2,
        );
        near(next.problem().contacts()[0].source.resistance_ohm_m2(), rc);
        let voltages = [
            (0, 0.0),
            (3, 0.0),
            (5, 2.0 * (ra + rc + 3.0)),
            (6, 2.0 * (ra + rc + 3.0)),
        ];
        let voltage_driven = layout
            .step_contact_candidate(
                &thermal,
                dt,
                &original,
                ConductionDrive::Voltage(&voltages),
                &contacts,
                limits(),
                field_controls(),
                config(),
                &[0.0; 2],
                || false,
            )
            .unwrap();
        near(
            voltage_driven.thermal().temperature_k(CellId(30)).unwrap(),
            300.0 + d1,
        );
        near(
            voltage_driven.thermal().temperature_k(CellId(10)).unwrap(),
            300.0 + d2,
        );
        near(
            voltage_driven.report().electrical.contacts[0]
                .response
                .power_w,
            i2 * rc,
        );
        assert_eq!(
            next.problem().contacts()[0]
                .source
                .resolved()
                .card_identity(),
            original.contacts()[0].source.resolved().card_identity()
        );
        let continuing = layout
            .step_contact_candidate(
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
        assert!(continuing.report().electrical.voltage_v[5] > step.voltage_v[5]);
        // Ordinary APIs may not silently throw away a surface heat account.
        assert!(
            layout
                .step_driven_candidate(
                    &thermal,
                    dt,
                    &original,
                    drive,
                    limits(),
                    field_controls(),
                    config(),
                    &[0.0; 2],
                    || false
                )
                .is_err()
        );
        work.push(step.heat_j);
    }
    assert!(work[1] > work[0]);
}

#[test]
fn g1_contact_field_drives_native_capacitor_discharge_and_exact_heat_split() {
    let (a, ca) = material(1.0, 1.0, 0.0);
    let (b, cb) = material(3.0, 2.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let thermal = initial(&curves, 999.0, None);
    let interface = interface(&cards, 0.5, 0.0, 1000.0);
    let context = point(300.0);
    let surface_context = interface_context(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &context, &key);
    let layout = network(&cells);
    let contacts = [binding(&cards, &interface, &surface_context, 0.25)];
    let circuit = ConductionCircuit::new(
        specimen(&cards, &interface),
        CircuitGraph {
            node_count: 2,
            branches: vec![
                (1, 0, Branch::Resistor { ohms: 999.0 }),
                (1, 0, Branch::Capacitor { farads: 1.0 }),
            ],
            transformers: vec![],
        },
        ConductionCircuitTerminal {
            branch: 0,
            positive: vec![5, 6],
            negative: vec![0, 3],
        },
        0.0,
        vec![],
        vec![9.0],
        &[],
        field_controls(),
        || false,
    )
    .unwrap();
    let dt = 0.5;
    let next = layout
        .step_contact_circuit_candidate(
            &thermal,
            dt,
            &circuit,
            &[],
            &[],
            &contacts,
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    // Independent midpoint RC solution; C=1, total resistance=4.5 ohms.
    let q1 = 9.0 * (1.0 - dt / 9.0) / (1.0 + dt / 9.0);
    let current = 0.5 * (9.0 + q1) / 4.5;
    let surface_heat = dt * current * current * 0.5;
    near(next.circuit().charges_c()[0], q1);
    let step = &next.report().electrical.spatial;
    near(step.contacts[0].response.power_w * dt, surface_heat);
    near(step.heat_j, 0.5 * (81.0 - q1 * q1));
    near(next.report().deposited_joule_j, step.heat_j);
    near(next.thermal().energy_j() - thermal.energy_j(), step.heat_j);
    let p1 = current * current * (1.0 + 0.25 * 0.5);
    let p2 = current * current * (3.0 + 0.75 * 0.5);
    let det = (1.0 + dt) * (2.0 + dt) - dt * dt;
    near(
        next.thermal().temperature_k(CellId(30)).unwrap(),
        300.0 + dt * ((2.0 + dt) * p1 + dt * p2) / det,
    );
    near(
        next.thermal().temperature_k(CellId(10)).unwrap(),
        300.0 + dt * (dt * p1 + (1.0 + dt) * p2) / det,
    );
    assert!(next.report().total_energy_residual_j.abs() < config().total_energy_tolerance_j);
}

#[test]
fn g4_contact_source_partition_support_and_cancellation_refuse_transactionally() {
    let (a, ca) = material(1.0, 1.0, 0.0);
    let (b, cb) = material(3.0, 2.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let thermal = initial(&curves, 999.0, None);
    let interface = interface(&cards, 0.5, 0.002, 300.01);
    let original = specimen(&cards, &interface);
    let context = point(300.0);
    let surface_context = interface_context(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &context, &key);
    let layout = network(&cells);
    let contacts = [binding(&cards, &interface, &surface_context, 0.25)];
    let drive = ConductionDrive::Current {
        positive: &[5, 6],
        negative: &[0, 3],
        current_a: 2.0,
    };
    let energy = thermal.energy_j();
    let source = original.contacts()[0].source.clone();
    assert!(
        layout
            .step_contact_candidate(
                &thermal,
                0.5,
                &original,
                drive,
                &contacts,
                limits(),
                field_controls(),
                config(),
                &[0.0; 2],
                || false
            )
            .is_err()
    );
    let mut wrong = contacts;
    wrong[0].heat_fraction_to_first = 1.1;
    assert!(
        layout
            .step_contact_candidate(
                &thermal,
                0.5,
                &original,
                drive,
                &wrong,
                limits(),
                field_controls(),
                config(),
                &[0.0; 2],
                || false
            )
            .is_err()
    );
    let mut cancelled = 0;
    assert!(matches!(
        layout.step_contact_candidate(
            &thermal,
            0.5,
            &original,
            drive,
            &contacts,
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || {
                cancelled += 1;
                cancelled > 10
            }
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    // Equal resistance is insufficient: the original surface temperatures
    // still have to match the actual receiving bodies exactly.
    let stale_context = surface_context.clone().with("Ta", 300.005).unwrap();
    let stale = PlanarConduction::new_with_contacts(
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
    assert!(
        layout
            .step_contact_candidate(
                &thermal,
                0.5,
                &stale,
                drive,
                &contacts,
                limits(),
                field_controls(),
                config(),
                &[0.0; 2],
                || false,
            )
            .is_err()
    );
    assert_eq!(thermal.energy_j().to_bits(), energy.to_bits());
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(original.contacts()[0].source, source);
    // A valid low-energy retry proves that refusal did not consume the state.
    layout
        .step_contact_candidate(
            &thermal,
            0.5,
            &original,
            ConductionDrive::Current {
                positive: &[5, 6],
                negative: &[0, 3],
                current_a: 0.001,
            },
            &contacts,
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
}

#[test]
fn g1_contact_phase_fractions_follow_latent_enthalpy_at_unchanged_temperature() {
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
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 1.0, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: 6.0,
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
        first_vertices: [1, 2],
        second_vertices: [4, 7],
        source: resolve_electrical_interface_state(
            &interface,
            &ip,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap(),
    }];
    // The original fixed-external API must still refuse unowned phase response.
    assert!(
        PlanarConduction::new_with_contacts(
            contact_mesh(),
            vec![conductor.clone(); 4],
            contacts.clone(),
            limits(),
            &mut || false
        )
        .is_err()
    );
    let problem = PlanarConduction::new_with_phase_contacts(
        contact_mesh(),
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
        .step_contact_candidate(
            &thermal,
            0.5,
            &problem,
            ConductionDrive::Current {
                positive: &[5, 6],
                negative: &[0, 3],
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
        next.report().electrical.voltage_v[5],
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
            .step_contact_candidate(
                &thermal,
                0.5,
                &problem,
                ConductionDrive::Current {
                    positive: &[5, 6],
                    negative: &[0, 3],
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

    // The original circuit path must retain exactly the same phase ownership.
    // Independent scalar balance for a 1 F capacitor starting at 9 C:
    // I_mid=9/(R_total+dt/2), R_total=3+2*fa, dt=0.5 s.
    // Only the surface-A fraction enters its latent-energy balance.
    let residual = |f: f64| {
        let current = 9.0 / (3.25 + 2.0 * f);
        100.0 * (f - 0.2) - 0.5 * current * current * (1.25 + 0.5 * f)
    };
    let (mut low, mut high) = (0.2, 0.5);
    assert!(residual(low) < 0.0 && residual(high) > 0.0);
    for _ in 0..64 {
        let middle = f64::midpoint(low, high);
        if residual(middle) < 0.0 {
            low = middle;
        } else {
            high = middle;
        }
    }
    let fa = f64::midpoint(low, high);
    let resistance = 1.0 + 2.0 * fa;
    let current = 9.0 / (2.25 + resistance);
    let charge = 9.0 - 0.5 * current;
    let fb = 0.2 + 0.5 * current * current * (1.0 + 0.75 * resistance) / 100.0;
    let circuit = ConductionCircuit::new(
        problem,
        CircuitGraph {
            node_count: 2,
            branches: vec![
                (1, 0, Branch::Resistor { ohms: 999.0 }),
                (1, 0, Branch::Capacitor { farads: 1.0 }),
            ],
            transformers: vec![],
        },
        ConductionCircuitTerminal {
            branch: 0,
            positive: vec![5, 6],
            negative: vec![0, 3],
        },
        0.0,
        vec![],
        vec![9.0],
        &[],
        field_controls(),
        || false,
    )
    .unwrap();
    let discharged = network
        .step_contact_circuit_candidate(
            &thermal,
            0.5,
            &circuit,
            &[],
            &[],
            &[contact_binding],
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    near(discharged.circuit().charges_c()[0], charge);
    near(
        discharged.thermal().temperature_k(CellId(30)).unwrap(),
        310.0,
    );
    near(
        discharged.thermal().temperature_k(CellId(10)).unwrap(),
        310.0,
    );
    let source = discharged.circuit().problem().contacts()[0]
        .source
        .resolved();
    near(source.state_coordinate("fa").unwrap(), fa);
    near(source.state_coordinate("fb").unwrap(), fb);
    near(
        discharged.report().electrical.spatial.contacts[0]
            .response
            .power_w,
        current * current * resistance,
    );
    near(
        discharged.report().deposited_joule_j,
        0.5 * (81.0 - charge * charge),
    );
    near(
        discharged.thermal().energy_j() - thermal.energy_j(),
        100.0 * (fa + fb - 0.4),
    );
}
