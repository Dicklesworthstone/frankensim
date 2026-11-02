//! Matching original tetrahedra, two source-side states and their actual heat.
//! The laws below are manufactured integration fixtures, not measured data.
use super::*;
use fs_conduction::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as THERMAL_RD,
    AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as THERMAL_RP,
    transient::deposition::contact::{ContactCell, MaterialContactLink},
};
use fs_couple::circuit_enthalpy::{
    conduction::contacts::ConductionContactThermal,
    contacts::{heat::MaterialHeatNetwork, material::ElectricalContactNetwork},
    material::MaterialNetwork,
};
use fs_exec::{CancelGate, Cx};
use fs_matdb::{ClaimSelection, InterfaceSystemCard, SurfaceSpec, SystemContext};
use fs_material::{
    electrical_interface::{
        ELECTRICAL_CONTACT_RESISTANCE_DIMS as ELECTRICAL_RD,
        ELECTRICAL_CONTACT_RESISTANCE_PROPERTY as ELECTRICAL_RP, ResolvedElectricalInterface,
        resolve_electrical_interface_state,
    },
    thermal::interface::ThermalInterfacePoint,
};
use fs_qty::{
    Pressure,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

fn static_quantity(kind: QuantityKind) -> QuantitySpec {
    // A source's named environmental role needs its physical quantity kind;
    // SI dimensions alone cannot distinguish pressure from evolving stress.
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

fn claim(key: PropertyKey, value: PropertyValue, validity: ValidityDomain) -> PropertyClaim {
    PropertyClaim {
        key,
        value,
        validity,
        interpolation: InterpolationPolicy::LinearInside,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "manufactured two-sided 3D material contact; not empirical data".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}

fn interface_card(cards: &[MaterialCard; 2], claims: ClaimSet) -> InterfaceSystemCard {
    InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: cards[0].id().clone(),
            texture_frame: "first original side".into(),
        },
        SurfaceSpec {
            material: cards[1].id().clone(),
            texture_frame: "second original side".into(),
        },
        SystemContext {
            medium: "explicit confined interface".into(),
            third_body: None,
            environment: "fixed authored conditions".into(),
            history: "fixed preparation".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn interface(cards: &[MaterialCard; 2], r0: f64, upper: f64) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    let temperature = static_quantity(QuantityKind::AbsoluteTemperature);
    claims
        .insert_claim(claim(
            PropertyKey::new(ELECTRICAL_RP, ELECTRICAL_RD),
            PropertyValue::Surface2D {
                axes: ["Ta".into(), "Tb".into()],
                axis_dims: [Temperature::DIMS; 2],
                samples: [
                    [250.0, 250.0],
                    [upper, 250.0],
                    [250.0, upper],
                    [upper, upper],
                ]
                .map(|p| {
                    (
                        p,
                        r0 * (1.0 + 0.003 * (p[0] - 300.0) + 0.007 * (p[1] - 300.0)),
                    )
                })
                .to_vec(),
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims: ELECTRICAL_RD,
            },
            ValidityDomain::unconstrained()
                .with_quantity("Ta", temperature, 250.0, upper)
                .with_quantity("Tb", temperature, 250.0, upper),
        ))
        .unwrap();
    claims
        .insert_claim(claim(
            PropertyKey::new(THERMAL_RP, THERMAL_RD),
            PropertyValue::Curve {
                abscissa: "pressure".into(),
                abscissa_dims: Pressure::DIMS,
                knots: vec![(1e5, 2.0), (2e5, 0.2)],
                dims: THERMAL_RD,
            },
            ValidityDomain::unconstrained().with_quantity(
                "pressure",
                static_quantity(QuantityKind::Pressure),
                1e5,
                2e5,
            ),
        ))
        .unwrap();
    interface_card(cards, claims)
}

fn interface_point(temperature: f64, pressure: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "Ta",
            static_quantity(QuantityKind::AbsoluteTemperature),
            temperature,
        )
        .unwrap()
        .with_quantity(
            "Tb",
            static_quantity(QuantityKind::AbsoluteTemperature),
            temperature,
        )
        .unwrap()
        .with_quantity(
            "pressure",
            static_quantity(QuantityKind::Pressure),
            pressure,
        )
        .unwrap()
}

fn surface_bindings<'a>(
    cards: &'a [MaterialCard; 2],
    interface: &'a InterfaceSystemCard,
    context: &'a QueryPoint,
) -> [ConductionContactThermal<'a>; 2] {
    [0, 1].map(|contact| ConductionContactThermal {
        contact,
        point: ThermalInterfacePoint {
            interface,
            materials: [&cards[0], &cards[1]],
            context,
            temperature_axes: ["Ta", "Tb"],
            liquid_fraction_axes: [None, None],
            mean_temperature_axis: None,
        },
        heat_fraction_to_first: 0.25,
    })
}

fn contacts(source: ResolvedElectricalInterface) -> Vec<SpatialConductionContact> {
    vec![
        SpatialConductionContact {
            first_vertices: [4, 5, 7],
            second_vertices: [8, 9, 11],
            source: source.clone(),
        },
        SpatialConductionContact {
            first_vertices: [4, 6, 7],
            second_vertices: [8, 10, 11],
            source,
        },
    ]
}

fn mesh() -> (Vec<[f64; 3]>, Vec<[u32; 4]>) {
    let (mut positions, mut tets) = slab(&[0.0, 1.0]);
    let (other, other_tets) = slab(&[1.0, 2.0]);
    positions.extend(other);
    tets.extend(other_tets.into_iter().map(|t| t.map(|v| v + 8)));
    (positions, tets)
}

fn specimen(
    cards: &[MaterialCard; 2],
    interface: &InterfaceSystemCard,
    context: &QueryPoint,
) -> SpatialConduction {
    let materials = cards.each_ref().map(|card| {
        resolve_conductor_state_with_basis(
            card,
            &point(300.0),
            MaterialPropertySelection::SingleClaimOnly,
            ElectricalPropertyBasis::Resistivity,
        )
        .unwrap()
    });
    let source = resolve_electrical_interface_state(
        interface,
        context,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let (positions, tets) = mesh();
    SpatialConduction::new_with_contacts(
        positions,
        tets,
        (0..12).map(|i| materials[i / 6].clone()).collect(),
        contacts(source),
        limits(),
        &mut || false,
    )
    .unwrap()
}

fn regions() -> Vec<SpatialConductionRegion> {
    // Deliberately reverse registration; routes retain original tetrahedra.
    (0..12)
        .rev()
        .map(|tetrahedron| SpatialConductionRegion {
            tetrahedron,
            cell: if tetrahedron < 6 {
                CellId(30)
            } else {
                CellId(10)
            },
        })
        .collect()
}

fn expected(r0: f64, conductance: f64) -> [f64; 2] {
    // Two independent backward-Euler volume balances at I=2 A, dt=.5 s,
    // Cp*m=(1,2), rho_A=1+.002*dTa, rho_B=3,
    // R''=r0*(1+.003*dTa+.007*dTb), authored interface split=(.25,.75).
    let a = 1.0 + 0.5 * conductance - 2.0 * (0.002 + 0.25 * r0 * 0.003);
    let b = -0.5 * conductance - 2.0 * 0.25 * r0 * 0.007;
    let c = -0.5 * conductance - 2.0 * 0.75 * r0 * 0.003;
    let d = 2.0 + 0.5 * conductance - 2.0 * 0.75 * r0 * 0.007;
    let r = 2.0 * (1.0 + 0.25 * r0);
    let s = 2.0 * (3.0 + 0.75 * r0);
    [
        (r * d - b * s) / (a * d - b * c),
        (a * s - c * r) / (a * d - b * c),
    ]
}

fn assert_balance<S>(
    next: &SpatialConductionEnthalpyCandidate<'_, ResolvedConductorState, S>,
    thermal: &DepositionSimulation<'_>,
    r0: f64,
    conductance: f64,
) {
    let [a, b] = expected(r0, conductance);
    let resistance = r0 * (1.0 + 0.003 * a + 0.007 * b);
    near(next.thermal().temperature_k(CellId(30)).unwrap(), 300.0 + a);
    near(next.thermal().temperature_k(CellId(10)).unwrap(), 300.0 + b);
    let step = &next.report().electrical;
    near(
        step.field.voltage_v[15],
        2.0 * (1.0 + 0.002 * a + resistance + 3.0),
    );
    near(
        step.field.cell_power_w.iter().sum(),
        4.0 * (4.0 + 0.002 * a),
    );
    assert_eq!(step.contacts.len(), 2);
    for (contact, source) in step.contacts.iter().zip(next.problem().contacts()) {
        assert_eq!(contact.cells, [CellId(30), CellId(10)]);
        near(contact.heat_fraction_to_first, 0.25);
        near(contact.response.current_a, -1.0);
        near(contact.response.power_w, 2.0 * resistance);
        near(source.source.resistance_ohm_m2(), resistance);
        near(
            source.source.resolved().state_coordinate("Ta").unwrap(),
            300.0 + a,
        );
        near(
            source.source.resolved().state_coordinate("Tb").unwrap(),
            300.0 + b,
        );
    }
    near(step.heat_j, 2.0 * (4.0 + 0.002 * a + resistance));
    near(step.supplied_work_j, step.heat_j);
    near(next.report().deposited_joule_j, step.heat_j);
    near(next.thermal().energy_j() - thermal.energy_j(), a + 2.0 * b);
    near(next.report().total_energy_residual_j, 0.0);
}

#[test]
fn two_source_side_temperatures_and_pressure_close_both_native_heat_entry_points() {
    for (r0, pressure) in [(0.5, 1e5), (0.5, 2e5), (1.5, 1e5)] {
        let (a, ca) = material(1.0, 1.0, 0.002);
        let (b, cb) = material(3.0, 2.0, 0.0);
        let cards = [a, b];
        let curves = [ca, cb];
        let context = point(300.0);
        let interface = interface(&cards, r0, 1000.0);
        let ip = interface_point(300.0, pressure);
        let surface = surface_bindings(&cards, &interface, &ip);
        let field = specimen(&cards, &interface, &ip);
        let thermal = initial(&curves, 999.0, None);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = bindings(&cards, &curves, &context, &key);
        let regions = regions();
        let layout = SpatialConductionThermalNetwork {
            cells: &cells,
            links: &LINKS,
            regions: &regions,
            volume_tolerance_m3: 1e-12,
        };
        let drive = boundary(&field, 2.0, 2.0, None);
        let [positive, negative] = electrodes(&field, 2.0);
        let current_drive = SpatialConductionDrive::total_current(&positive, &negative, 2.0);
        let bare = layout
            .step_contact_candidate(
                &thermal,
                0.5,
                &field,
                &drive,
                &surface,
                limits(),
                field_controls(),
                config(),
                &[0.0; 2],
                || false,
            )
            .unwrap();
        assert_balance(&bare, &thermal, r0, 1.0);
        let current_bare = layout
            .step_contact_candidate(
                &thermal,
                0.5,
                &field,
                current_drive,
                &surface,
                limits(),
                field_controls(),
                config(),
                &[0.0; 2],
                || false,
            )
            .unwrap();
        assert_balance(&current_bare, &thermal, r0, 1.0);
        assert!(
            layout
                .step_candidate(
                    &thermal,
                    0.5,
                    &field,
                    &drive,
                    limits(),
                    field_controls(),
                    config(),
                    &[0.0; 2],
                    || false
                )
                .is_err()
        );
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
            let prepared = heat.prepare(&cx, &thermal).unwrap();
            let next = prepared
                .step_spatial_conduction_contact_candidate(
                    &cx,
                    &thermal,
                    0.5,
                    &field,
                    current_drive,
                    &surface,
                    &regions,
                    1e-12,
                    limits(),
                    field_controls(),
                    config(),
                    &[0.0; 2],
                )
                .unwrap();
            let thermal_resistance = 2.0 - 1.8 * (pressure - 1e5) / 1e5;
            assert_balance(&next, &thermal, r0, 1.0 / (1.0 + thermal_resistance));
            near(
                next.report().source.thermal_contacts[0]
                    .interface
                    .value_m2_k_per_w(),
                thermal_resistance,
            );
            let voltage = next.report().electrical.field.voltage_v[15];
            let voltage_driven = prepared
                .step_spatial_conduction_contact_candidate(
                    &cx,
                    &thermal,
                    0.5,
                    &field,
                    &boundary(&field, 2.0, 0.0, Some(voltage)),
                    &surface,
                    &regions,
                    1e-12,
                    limits(),
                    field_controls(),
                    config(),
                    &[0.0; 2],
                )
                .unwrap();
            assert_balance(
                &voltage_driven,
                &thermal,
                r0,
                1.0 / (1.0 + thermal_resistance),
            );
            let continued = prepared
                .step_spatial_conduction_contact_candidate(
                    &cx,
                    next.thermal(),
                    1.0,
                    next.problem(),
                    current_drive,
                    &surface,
                    &regions,
                    1e-12,
                    limits(),
                    field_controls(),
                    config(),
                    &[0.0; 2],
                )
                .unwrap();
            assert!(
                continued.problem().contacts()[0].source.resistance_ohm_m2()
                    > next.problem().contacts()[0].source.resistance_ohm_m2()
            );
            assert_eq!(
                continued.problem().contacts()[0]
                    .source
                    .resolved()
                    .card_identity(),
                field.contacts()[0].source.resolved().card_identity()
            );
        });
        assert_eq!(thermal.time_s(), 0.0);
    }
}

#[test]
fn both_latent_fractions_refresh_contacts_without_a_temperature_change() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let (card, curve, keys) = phase_material(basis, 0.0);
    let cards = [card.clone(), card];
    let scalar = QuantitySpec::dimensional(Dims::NONE);
    let context = point(400.0).with_quantity("f", scalar, 0.2).unwrap();
    let requirements =
        [
            ScalarPropertyRequirement::try_with_key(
                &keys[0],
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap(),
        ];
    let snapshot = resolve_material_state_point(
        &cards[0],
        &context,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let conductor =
        ResolvedConductorState::try_from_material_with_property(snapshot, basis, &keys[0])
            .unwrap()
            .with_liquid_mass_fraction_axis("f")
            .unwrap();
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(claim(
            PropertyKey::new(ELECTRICAL_RP, ELECTRICAL_RD),
            PropertyValue::Surface2D {
                axes: ["fa".into(), "fb".into()],
                axis_dims: [Dims::NONE; 2],
                samples: [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]
                    .map(|p| (p, 1.0 + 2.0 * p[0] + 3.0 * p[1]))
                    .to_vec(),
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims: ELECTRICAL_RD,
            },
            ValidityDomain::unconstrained()
                .with_quantity("fa", scalar, 0.0, 1.0)
                .with_quantity("fb", scalar, 0.0, 1.0),
        ))
        .unwrap();
    let interface = interface_card(&cards, claims);
    let ip = interface_point(400.0, 1e5)
        .with_quantity("fa", scalar, 0.2)
        .unwrap()
        .with_quantity("fb", scalar, 0.2)
        .unwrap();
    let source = resolve_electrical_interface_state(
        &interface,
        &ip,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let (positions, tets) = mesh();
    let field = SpatialConduction::new_with_phase_contacts(
        positions,
        tets,
        vec![conductor; 12],
        contacts(source),
        vec![[Some("fa".into()), Some("fb".into())]; 2],
        limits(),
        &mut || false,
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
                specific_h_j_kg: 20.0,
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
    let cells = [CellId(30), CellId(10)]
        .into_iter()
        .enumerate()
        .map(|(i, cell)| MaterialCell {
            cell,
            material_id: MaterialId(i as u64),
            curve: &curve,
            point: ThermalMaterialPoint {
                material: &cards[i],
                context: &context,
                temperature_axis: "T",
                liquid_fraction_axis: Some("f"),
            },
            conductivity: None,
        })
        .collect::<Vec<_>>();
    let regions = regions();
    let layout = SpatialConductionThermalNetwork {
        cells: &cells,
        links: &[],
        regions: &regions,
        volume_tolerance_m3: 1e-12,
    };
    let mut surfaces = surface_bindings(&cards, &interface, &ip);
    for surface in &mut surfaces {
        surface.point.liquid_fraction_axes = [Some("fa"), Some("fb")];
    }
    let [positive, negative] = electrodes(&field, 2.0);
    let drive = SpatialConductionDrive::total_current(&positive, &negative, 2.0);
    let next = layout
        .step_contact_candidate(
            &thermal,
            0.5,
            &field,
            drive,
            &surfaces,
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    let trace_next = layout
        .step_contact_candidate(
            &thermal,
            0.5,
            &field,
            &boundary(&field, 2.0, 2.0, None),
            &surfaces,
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    near(trace_next.thermal().energy_j(), next.thermal().energy_j());
    for (trace, total) in trace_next
        .report()
        .electrical
        .field
        .voltage_v
        .iter()
        .zip(&next.report().electrical.field.voltage_v)
    {
        near(*trace, *total);
    }
    // Independent latent balances: 99 fa - 1.5 fb = 24.5,
    // -3 fa + 95.5 fb = 25.5. Both temperatures stay at 400 K.
    let det = 99.0 * 95.5 - 4.5;
    let fa = (24.5 * 95.5 + 1.5 * 25.5) / det;
    let fb = (99.0 * 25.5 + 3.0 * 24.5) / det;
    let resistance = 1.0 + 2.0 * fa + 3.0 * fb;
    for cell in [CellId(30), CellId(10)] {
        near(next.thermal().temperature_k(cell).unwrap(), 400.0);
    }
    for contact in next.problem().contacts() {
        near(
            contact.source.resolved().state_coordinate("fa").unwrap(),
            fa,
        );
        near(
            contact.source.resolved().state_coordinate("fb").unwrap(),
            fb,
        );
        near(contact.source.resistance_ohm_m2(), resistance);
    }
    near(
        next.report().electrical.field.voltage_v[15],
        2.0 * (4.0 + resistance),
    );
    near(next.report().deposited_joule_j, 2.0 * (4.0 + resistance));
    near(
        next.thermal().energy_j() - thermal.energy_j(),
        100.0 * (fa + fb - 0.4),
    );
    surfaces[0].point.liquid_fraction_axes = [Some("fa"), None];
    assert!(
        layout
            .step_contact_candidate(
                &thermal,
                0.5,
                &field,
                drive,
                &surfaces,
                limits(),
                field_controls(),
                config(),
                &[0.0; 2],
                || false
            )
            .is_err()
    );
    near(thermal.energy_j(), 40.0);
}

#[test]
fn source_support_and_late_cancellation_preserve_both_accepted_histories() {
    for upper in [300.05, 1000.0] {
        let (a, ca) = material(1.0, 1.0, 0.002);
        let (b, cb) = material(3.0, 2.0, 0.0);
        let cards = [a, b];
        let curves = [ca, cb];
        let context = point(300.0);
        let interface = interface(&cards, 0.5, upper);
        let ip = interface_point(300.0, 1e5);
        let surface = surface_bindings(&cards, &interface, &ip);
        let field = specimen(&cards, &interface, &ip);
        let thermal = initial(&curves, 999.0, None);
        let saved = thermal.cells().clone();
        let original_contacts = field.contacts().to_vec();
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = bindings(&cards, &curves, &context, &key);
        let regions = regions();
        let layout = SpatialConductionThermalNetwork {
            cells: &cells,
            links: &LINKS,
            regions: &regions,
            volume_tolerance_m3: 1e-12,
        };
        let drive = boundary(&field, 2.0, 2.0, None);
        let mut polls = 0;
        let result = layout.step_contact_candidate(
            &thermal,
            0.5,
            &field,
            &drive,
            &surface,
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || {
                polls += 1;
                false
            },
        );
        if upper < 301.0 {
            assert!(result.is_err());
        } else {
            let next = result.unwrap();
            assert_balance(&next, &thermal, 0.5, 1.0);
            let mut calls = 0;
            assert!(matches!(
                layout.step_contact_candidate(
                    &thermal,
                    0.5,
                    &field,
                    &drive,
                    &surface,
                    limits(),
                    field_controls(),
                    config(),
                    &[0.0; 2],
                    || {
                        calls += 1;
                        calls == polls
                    }
                ),
                Err(CircuitEnthalpyError::Cancelled)
            ));
            let mut wrong = surface;
            wrong[1].heat_fraction_to_first = f64::NAN;
            assert!(
                layout
                    .step_contact_candidate(
                        &thermal,
                        0.5,
                        &field,
                        &drive,
                        &wrong,
                        limits(),
                        field_controls(),
                        config(),
                        &[0.0; 2],
                        || false
                    )
                    .is_err()
            );
        }
        assert_eq!(thermal.cells(), &saved);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(field.contacts(), original_contacts);
    }
}

#[test]
fn capacitor_discharge_refreshes_two_sided_contact_sources_and_deposits_surface_work_once() {
    use fs_couple::circuit_enthalpy::conduction::spatial::circuit::SpatialConductionCircuit;
    use fs_phs::circuit::{Branch, CircuitGraph};
    let r0 = 0.75;
    for upper in [1000.0, 300.0001] {
        let (a, ca) = material(1.0, 1.0, 0.002);
        let (b, cb) = material(3.0, 2.0, 0.0);
        let cards = [a, b];
        let curves = [ca, cb];
        let context = point(300.0);
        let interface = interface(&cards, r0, upper);
        let ip = interface_point(300.0, 1e5);
        let surface = surface_bindings(&cards, &interface, &ip);
        let field = specimen(&cards, &interface, &ip);
        let [positive, negative] = electrodes(&field, 2.0);
        let thermal = initial(&curves, 999.0, None);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = bindings(&cards, &curves, &context, &key);
        let regions = regions();
        let layout = SpatialConductionThermalNetwork {
            cells: &cells,
            links: &LINKS,
            regions: &regions,
            volume_tolerance_m3: 1e-12,
        };
        let response = SpatialConductanceControls {
            maximum_terminals: 2,
            maximum_field_solves: 1,
            maximum_response_values: 1000,
            maximum_dense_work: 1000,
            reference_voltage_v: 1.0,
            coefficient_tolerance_s: 1e-10,
        };
        let circuit = SpatialConductionCircuit::new(
            field,
            CircuitGraph {
                node_count: 2,
                branches: vec![(1, 0, Branch::Capacitor { farads: 1.0 })],
                transformers: vec![],
            },
            vec![negative, positive],
            vec![0, 1],
            0.0,
            vec![],
            vec![4.0],
            &[],
            response,
            field_controls(),
            || false,
        )
        .unwrap();
        let next = layout.step_circuit_candidate(
            &thermal,
            0.5,
            &circuit,
            &[],
            &surface,
            response,
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || false,
        );
        if upper < 301.0 {
            assert!(next.is_err());
            near(circuit.charges_c()[0], 4.0);
            assert_eq!(thermal.time_s(), 0.0);
            near(
                circuit.problem().contacts()[0]
                    .source
                    .resolved()
                    .state_coordinate("Ta")
                    .unwrap(),
                300.0,
            );
            continue;
        }
        let next = next.unwrap();
        // Independent BE elimination for a proposed midpoint current, then
        // capacitor closure I*(R_total(Ta,Tb)+dt/(2C))=q_old.
        let endpoint = |current: f64| {
            let s = 0.5 * current * current;
            let a = 1.5 - s * (0.002 + 0.25 * r0 * 0.003);
            let b = -0.5 - s * 0.25 * r0 * 0.007;
            let c = -0.5 - s * 0.75 * r0 * 0.003;
            let d = 2.5 - s * 0.75 * r0 * 0.007;
            let ra = s * (1.0 + 0.25 * r0);
            let rb = s * (3.0 + 0.75 * r0);
            let theta = [
                (ra * d - b * rb) / (a * d - b * c),
                (a * rb - c * ra) / (a * d - b * c),
            ];
            let contact = r0 * (1.0 + 0.003 * theta[0] + 0.007 * theta[1]);
            let resistance = 4.0 + 0.002 * theta[0] + contact;
            (
                theta,
                contact,
                resistance,
                current * (resistance + 0.25) - 4.0,
            )
        };
        let (mut lo, mut hi) = (0.0, 4.0 / (4.0 + r0 + 0.25));
        assert!(endpoint(lo).3 < 0.0 && endpoint(hi).3 > 0.0);
        for _ in 0..90 {
            let mid = 0.5 * (lo + hi);
            if endpoint(mid).3 < 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let current = 0.5 * (lo + hi);
        let (theta, contact, resistance, _) = endpoint(current);
        let final_charge = 4.0 - 0.5 * current;
        near(next.circuit().charges_c()[0], final_charge);
        near(
            next.thermal().temperature_k(CellId(30)).unwrap(),
            300.0 + theta[0],
        );
        near(
            next.thermal().temperature_k(CellId(10)).unwrap(),
            300.0 + theta[1],
        );
        let step = &next.report().electrical;
        near(step.terminal_current_a[1], current);
        near(step.terminal_voltage_v[1], current * resistance);
        near(step.spatial.heat_j, 0.5 * current * current * resistance);
        near(step.supplied_work_j, 0.0);
        near(step.native.supplied, 0.0);
        near(
            step.stored_change_j,
            0.5 * (final_charge * final_charge - 16.0),
        );
        near(step.native.dissipated, step.spatial.heat_j);
        near(step.energy_residual_j, 0.0);
        near(next.report().total_energy_residual_j, 0.0);
        near(
            next.thermal().energy_j() - thermal.energy_j(),
            step.spatial.heat_j,
        );
        assert_eq!(step.spatial.contacts.len(), 2);
        for (heat, source) in step
            .spatial
            .contacts
            .iter()
            .zip(next.circuit().problem().contacts())
        {
            near(heat.heat_fraction_to_first, 0.25);
            near(heat.response.current_a, -0.5 * current);
            near(heat.response.power_w, 0.5 * current * current * contact);
            near(source.source.resistance_ohm_m2(), contact);
            near(
                source.source.resolved().state_coordinate("Ta").unwrap(),
                300.0 + theta[0],
            );
            near(
                source.source.resolved().state_coordinate("Tb").unwrap(),
                300.0 + theta[1],
            );
        }
        assert_eq!(thermal.time_s(), 0.0);
        near(circuit.charges_c()[0], 4.0);
    }
}

#[path = "contacts/material_resistors.rs"]
mod material_resistors;
