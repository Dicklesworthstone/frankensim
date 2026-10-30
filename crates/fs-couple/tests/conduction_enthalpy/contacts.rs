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

fn interface(cards: &[MaterialCard; 2], r0: f64, beta: f64, upper: f64) -> InterfaceSystemCard {
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

fn interface_context(t: f64) -> QueryPoint {
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

fn binding<'a>(
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

fn specimen(cards: &[MaterialCard; 2], interface: &InterfaceSystemCard) -> PlanarConduction {
    let mesh = TriComplex2::from_indexed_triangles(
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
    .unwrap();
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
