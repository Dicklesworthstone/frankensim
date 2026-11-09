//! Native current-dependent contact work heats both original finite bodies.
use super::*;
use crate::response::{boundary_response, controls};
use fs_couple::circuit_enthalpy::conduction::spatial::response::{
    SpatialConductionResponseEnthalpyCandidate, step_spatial_conduction_response_candidate,
};
use fs_em::{
    conduction::response::{
        contact::ContactResponseControls,
        terminals::{ConductionTerminalControl, ConductionTerminalDrive, ConductionTerminalLimits},
    },
    spatial::conduction::response::SpatialConductionResponse,
};
use fs_material::{
    conductor::response::ScalarConductorResponse,
    electrical_interface::response::{
        ElectricalInterfaceCoordinate, ScalarElectricalInterfaceResponse,
    },
};

fn nonlinear_interface(cards: &[MaterialCard; 2]) -> InterfaceSystemCard {
    let original = interface(cards, 0.5, 2000.0);
    let mut claims = ClaimSet::new();
    for (_, original) in original.claims().claims_ordered() {
        if original.key.name() != ELECTRICAL_RP {
            claims.insert_claim(original.clone()).unwrap();
        }
    }
    let current = ElectricalInterfaceCoordinate::CurrentDensityMagnitude;
    let temperature = static_quantity(QuantityKind::AbsoluteTemperature);
    claims
        .insert_claim(claim(
            PropertyKey::new(ELECTRICAL_RP, ELECTRICAL_RD),
            PropertyValue::Surface2D {
                axes: ["Tm".into(), "contact_J".into()],
                axis_dims: [Temperature::DIMS, current.quantity().dims()],
                samples: [[250.0, 0.0], [2000.0, 0.0], [250.0, 20.0], [2000.0, 20.0]]
                    .into_iter()
                    .map(|p| (p, 0.5 + 0.002 * (p[0] - 300.0) + 0.1 * p[1]))
                    .collect(),
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims: ELECTRICAL_RD,
            },
            ValidityDomain::unconstrained()
                .with_quantity("Ta", temperature, 250.0, 2000.0)
                .with_quantity("Tb", temperature, 250.0, 2000.0)
                .with_quantity("Tm", temperature, 250.0, 2000.0)
                .with_quantity("contact_J", current.quantity(), 0.0, 20.0),
        ))
        .unwrap();
    interface_card(cards, claims)
}

fn nonlinear_problem(
    cards: &[MaterialCard; 2],
    card: &InterfaceSystemCard,
    context: &QueryPoint,
) -> SpatialConductionResponse {
    let sources = (0..12)
        .map(|tet| {
            let card = &cards[tet / 6];
            let conductor = resolve_conductor_state_with_basis(
                card,
                &point(300.0),
                MaterialPropertySelection::SingleClaimOnly,
                ElectricalPropertyBasis::Resistivity,
            )
            .unwrap();
            ScalarConductorResponse::new(card, conductor, &[]).unwrap()
        })
        .collect();
    let source = resolve_electrical_interface_state(
        card,
        context,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let law = ScalarElectricalInterfaceResponse::new(
        card,
        source.clone(),
        &[(
            "contact_J",
            ElectricalInterfaceCoordinate::CurrentDensityMagnitude,
        )],
    )
    .unwrap();
    let (positions, tets) = mesh();
    SpatialConductionResponse::new_with_contact_response(
        positions,
        tets,
        sources,
        contacts(source),
        vec![law.clone(), law],
        ContactResponseControls {
            absolute_resistance_tolerance_ohm_m2: 1e-11,
            relative_resistance_tolerance: 1e-11,
            voltage_tolerance_v: 1e-10,
            current_density_tolerance_a_m2: 1e-10,
        },
        limits(),
        &mut || false,
    )
    .unwrap()
    .with_relaxation(0.5)
    .unwrap()
}

// Independent two-cell BE balance, J=2, dt=.5, C=[1,2],
// rho_a=1+.002*theta_a; rho_b=3; R''=.5+.1*J+.001*(theta_a+theta_b).
fn reference(old: [f64; 2], g: f64) -> ([f64; 2], f64) {
    let (p, f, rc0, slope) = (2.0, 0.25, 0.7, 0.001);
    let a = 1.0 + 0.5 * g - p * (0.002 + f * slope);
    let b = -0.5 * g - p * f * slope;
    let c = -0.5 * g - p * (1.0 - f) * slope;
    let d = 2.0 + 0.5 * g - p * (1.0 - f) * slope;
    let x = old[0] + p * (1.0 + f * rc0);
    let y = 2.0 * old[1] + p * (3.0 + (1.0 - f) * rc0);
    let theta = [
        (x * d - b * y) / (a * d - b * c),
        (a * y - c * x) / (a * d - b * c),
    ];
    (theta, rc0 + slope * (theta[0] + theta[1]))
}

fn assert_response<S>(
    next: &SpatialConductionResponseEnthalpyCandidate<'_, S>,
    thermal: &DepositionSimulation<'_>,
    old: [f64; 2],
    g: f64,
) -> [f64; 2] {
    let (theta, resistance) = reference(old, g);
    near(
        next.thermal().temperature_k(CellId(30)).unwrap(),
        300.0 + theta[0],
    );
    near(
        next.thermal().temperature_k(CellId(10)).unwrap(),
        300.0 + theta[1],
    );
    let step = &next.report().electrical;
    near(
        step.field.field.voltage_v[15],
        2.0 * (4.0 + 0.002 * theta[0] + resistance),
    );
    near(
        step.field.field.cell_power_w.iter().sum(),
        4.0 * (4.0 + 0.002 * theta[0]),
    );
    assert_eq!(step.contacts.len(), 2);
    assert_eq!(step.field.evaluated_contacts.len(), 2);
    for (contact, sample) in step.contacts.iter().zip(&step.field.evaluated_contacts) {
        assert_eq!(contact.cells, [CellId(30), CellId(10)]);
        near(contact.heat_fraction_to_first, 0.25);
        near(contact.response.current_a, -1.0);
        near(contact.response.power_w, 2.0 * resistance);
        for source in &sample.evaluated {
            near(source.resistance_ohm_m2(), resistance);
            near(
                source.resolved().state_coordinate("Ta").unwrap(),
                300.0 + theta[0],
            );
            near(
                source.resolved().state_coordinate("Tb").unwrap(),
                300.0 + theta[1],
            );
            near(
                source.resolved().state_coordinate("Tm").unwrap(),
                300.0 + (theta[0] + theta[1]) / 2.0,
            );
            near(
                source.resolved().state_coordinate("contact_J").unwrap(),
                2.0,
            );
        }
    }
    near(step.heat_j, 2.0 * (4.0 + 0.002 * theta[0] + resistance));
    near(step.supplied_work_j, step.heat_j);
    near(next.report().deposited_joule_j, step.heat_j);
    near(next.thermal().energy_j() - thermal.energy_j(), step.heat_j);
    near(next.report().total_energy_residual_j, 0.0);
    assert_eq!(next.problem().relaxation(), 0.5);
    assert!(next.problem().contact_response_sources().is_some());
    theta
}

#[test]
fn current_dependent_contact_heats_both_bodies_and_complete_sourced_heat_feedback() {
    let (a, ca) = material(1.0, 1.0, 0.002);
    let (b, cb) = material(3.0, 2.0, 0.0);
    let (cards, curves) = ([a, b], [ca, cb]);
    let card = nonlinear_interface(&cards);
    let ip = interface_point(300.0, 1e5)
        .with_quantity(
            "Tm",
            static_quantity(QuantityKind::AbsoluteTemperature),
            300.0,
        )
        .unwrap()
        .with_quantity(
            "contact_J",
            ElectricalInterfaceCoordinate::CurrentDensityMagnitude.quantity(),
            0.0,
        )
        .unwrap();
    let mut surface = surface_bindings(&cards, &card, &ip);
    for contact in &mut surface {
        contact.point.mean_temperature_axis = Some("Tm");
    }
    let context = point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = bindings(&cards, &curves, &context, &key);
    let regions = regions();
    let layout = SpatialConductionThermalNetwork {
        cells: &cells,
        links: &LINKS,
        regions: &regions,
        volume_tolerance_m3: 1e-12,
    };
    let original = nonlinear_problem(&cards, &card, &ip);
    let initial = initial(&curves, 999.0, None);
    let drive = boundary_response(&original, 2.0, None);
    let mut polls = 0;
    let first = layout
        .step_response_contact_candidate(
            &initial,
            0.5,
            &original,
            &drive,
            &surface,
            limits(),
            field_controls(),
            controls(),
            config(),
            &[0.0; 2],
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let theta = assert_response(&first, &initial, [0.0; 2], 1.0);
    let electrodes = [0.0, 2.0].map(|x| {
        original
            .exterior_faces()
            .filter(|&face| {
                original.mesh().faces[face].map(|v| original.positions_m()[v as usize][0]) == [x; 3]
            })
            .collect::<Vec<_>>()
    });
    let terminals = electrodes.each_ref().map(Vec::as_slice);
    let direct = layout
        .step_response_contact_candidate(
            &initial,
            0.5,
            &original,
            ConductionTerminalDrive {
                terminals: &terminals,
                control: ConductionTerminalControl::Currents(&[-2.0, 2.0]),
                limits: ConductionTerminalLimits {
                    maximum_terminals: 2,
                    maximum_workspace_values: 20_000,
                },
            },
            &surface,
            limits(),
            field_controls(),
            controls(),
            config(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    assert_response(&direct, &initial, [0.0; 2], 1.0);
    let actual = direct.report().electrical.field.terminals.as_ref().unwrap();
    near(actual.current_a[0], -2.0);
    near(actual.current_a[1], 2.0);
    near(
        0.5 * actual.supplied_power_w,
        direct.report().electrical.heat_j,
    );
    assert_eq!(
        direct.report().electrical.field.field_solves,
        direct.report().electrical.field.iterations,
    );
    let continued = layout
        .step_response_contact_candidate(
            first.thermal(),
            1.0,
            first.problem(),
            &drive,
            &surface,
            limits(),
            field_controls(),
            controls(),
            config(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    assert_response(&continued, first.thermal(), theta, 1.0);
    let mut stop = 0;
    assert!(matches!(
        layout.step_response_contact_candidate(
            &initial,
            0.5,
            &original,
            &drive,
            &surface,
            limits(),
            field_controls(),
            controls(),
            config(),
            &[0.0; 2],
            || {
                stop += 1;
                stop >= polls - 1
            },
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    assert_eq!(initial.time_s(), 0.0);
    near(original.contacts()[0].source.resistance_ohm_m2(), 0.5);
    assert_eq!(ip.axes()["contact_J"], 0.0);

    let bulk = MaterialNetwork {
        cells: &cells,
        resistors: &[],
        links: &LINKS,
    };
    let electrical = ElectricalContactNetwork {
        bulk: &bulk,
        contacts: &[],
        electrical_bases: &[],
    };
    let thermal_contacts = [MaterialContactLink {
        link: 7,
        point: surface[0].point,
        sides: [0, 1].map(|i| ContactCell {
            cell: cells[i].cell,
            material_id: cells[i].material_id,
            curve: curves[i].curve(),
        }),
        selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
    }];
    let heat = MaterialHeatNetwork {
        electrical: &electrical,
        thermal_contacts: &thermal_contacts,
        radiation: None,
    };
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 73,
                kernel_id: 51,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        let prepared = heat.prepare(&cx, &initial).unwrap();
        let next = step_spatial_conduction_response_candidate(
            &prepared,
            &cx,
            &initial,
            0.5,
            &original,
            &drive,
            &surface,
            &regions,
            1e-12,
            limits(),
            field_controls(),
            controls(),
            config(),
            &[0.0; 2],
        )
        .unwrap();
        assert_response(&next, &initial, [0.0; 2], 1.0 / 3.0);
        near(
            next.report().source.thermal_contacts[0]
                .interface
                .value_m2_k_per_w(),
            2.0,
        );

        // A differently named voltage applicability axis is still electrical.
        // The heat source has no owner for that work and cannot retain a V=0
        // slice after the original contact voltage changes.
        let voltage = ElectricalInterfaceCoordinate::VoltageJump.quantity();
        let mut claims = ClaimSet::new();
        for (_, original) in card.claims().claims_ordered() {
            let mut claim = original.clone();
            if claim.key.name() == THERMAL_RP {
                claim.validity = claim.validity.with_quantity("heat_voltage", voltage, 0.0, 20.0);
            }
            claims.insert_claim(claim).unwrap();
        }
        let bounded_card = interface_card(&cards, claims);
        let bounded_context = ip.clone().with_quantity("heat_voltage", voltage, 0.0).unwrap();
        let bounded_contacts = [MaterialContactLink {
            link: 7,
            point: ThermalInterfacePoint {
                interface: &bounded_card, context: &bounded_context, ..surface[0].point
            },
            sides: [0, 1].map(|i| ContactCell {
                cell: cells[i].cell, material_id: cells[i].material_id, curve: curves[i].curve(),
            }),
            selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        }];
        let bounded_heat = MaterialHeatNetwork {
            electrical: &electrical, thermal_contacts: &bounded_contacts, radiation: None,
        };
        let bounded = bounded_heat.prepare(&cx, &initial).unwrap();
        assert!(matches!(step_spatial_conduction_response_candidate(
            &bounded, &cx, &initial, 0.5, &original, &drive, &surface, &regions,
            1e-12, limits(), field_controls(), controls(), config(), &[0.0; 2],
        ), Err(CircuitEnthalpyError::Source(
            fs_couple::circuit_enthalpy::conduction::heat::ConductionHeatError::Heat(
                fs_couple::circuit_enthalpy::contacts::heat::HeatNetworkError::Invalid(
                    "thermal contacts and emissivity require field-independent response and applicability"
                )
            )
        ))));
    });
}
