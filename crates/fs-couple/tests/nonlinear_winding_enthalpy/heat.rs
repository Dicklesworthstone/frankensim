//! Independent two-temperature circuit/field/radiation/contact balance.
use super::*;
#[path = "../../../fs-conduction/tests/radiation_material/common.rs"]
mod radiative;
use fs_conduction::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RD, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RP,
    EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY, SurfaceEmissivity,
    transient::deposition::{contact::*, radiation::*},
};
use fs_couple::circuit_enthalpy::{
    contacts::{heat::*, material::ElectricalContactNetwork},
    material::{CellResistor, ConductionLink, MaterialNetwork},
};
use fs_exec::{CancelGate, Cx};
use fs_matdb::*;
use fs_material::{
    state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS},
    thermal::interface::ThermalInterfacePoint,
};
use fs_qty::Temperature;

const SINK: CellId = CellId(8);
const SINK_MATERIAL: MaterialId = MaterialId(10);
const SLOPE: f64 = 0.004;
const GRAPH_ALPHA: f64 = 0.002;
const DESTINATIONS: [ResistorCell; 1] = [ResistorCell {
    branch: 2,
    cell: SINK,
}];
const LINKS: [ConductionLink; 1] = [ConductionLink {
    link: 17,
    distance_a_m: 0.5,
    distance_b_m: 0.5,
}];

fn claim(name: &str, dims: fs_qty::Dims, value: PropertyValue, upper: f64) -> PropertyClaim {
    let interpolation = if matches!(value, PropertyValue::Curve { .. }) {
        InterpolationPolicy::LinearInside
    } else {
        InterpolationPolicy::ConstantWithinValidity
    };
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value,
        validity: fs_evidence::ValidityDomain::unconstrained().with("T", 200.0, upper),
        uncertainty: UncertaintyModel::Unstated,
        interpolation,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic magnetic/contact/radiation fixture".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}

fn cards(e: f64, upper: f64) -> [MaterialCard; 2] {
    [
        card(SLOPE),
        linear::card(0.5, GRAPH_ALPHA, 1.0, ElectricalPropertyBasis::Resistivity),
    ]
    .into_iter()
    .enumerate()
    .map(|(i, base)| {
        let mut claims = base.claims().clone();
        claims
            .insert_claim(claim(
                "thermal_conductivity",
                THERMAL_CONDUCTIVITY_DIMS,
                PropertyValue::Scalar {
                    value: 1.0,
                    dims: THERMAL_CONDUCTIVITY_DIMS,
                },
                500.0,
            ))
            .unwrap();
        let epsilon = if i == 0 { e } else { 0.6 };
        claims
            .insert_claim(claim(
                SURFACE_EMISSIVITY_PROPERTY,
                EMISSIVITY_DIMS,
                PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: [200.0, upper]
                        .map(|t| (t, epsilon * (1.0 + 0.0001 * (t - 300.0))))
                        .to_vec(),
                    dims: EMISSIVITY_DIMS,
                },
                upper,
            ))
            .unwrap();
        MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap()
    })
    .collect::<Vec<_>>()
    .try_into()
    .unwrap()
}

fn interface(cards: &[MaterialCard; 2], resistance: f64) -> InterfaceSystemCard {
    let mut value = claim(
        RP,
        RD,
        PropertyValue::Curve {
            abscissa: "aT".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(200.0, 0.9 * resistance), (500.0, 1.2 * resistance)],
            dims: RD,
        },
        500.0,
    );
    value.validity = fs_evidence::ValidityDomain::unconstrained()
        .with("aT", 200.0, 500.0)
        .with("bT", 200.0, 500.0);
    let mut claims = ClaimSet::new();
    claims.insert_claim(value).unwrap();
    InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: cards[0].id().clone(),
            texture_frame: "a".into(),
        },
        SurfaceSpec {
            material: cards[1].id().clone(),
            texture_frame: "b".into(),
        },
        SystemContext {
            medium: "fixed contact".into(),
            third_body: None,
            environment: "authored".into(),
            history: "fixed".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn with_heat<R>(
    e: f64,
    resistance: f64,
    upper: f64,
    run: impl FnOnce(
        &Cx<'_>,
        &CancelGate,
        &MaterialHeatNetwork<'_>,
        DepositionSimulation<'_>,
        &PlanarNonlinearMagnetostatics,
        &[MagneticThermalLaw],
    ) -> R,
) -> R {
    let cards = cards(e, upper);
    let curves = cards.each_ref().map(caloric);
    let context = linear::absolute_point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = [0, 1].map(|i| MaterialCell {
        cell: if i == 0 { CELL } else { SINK },
        material_id: if i == 0 { MATERIAL } else { SINK_MATERIAL },
        curve: curves[i].curve(),
        point: ThermalMaterialPoint {
            material: &cards[i],
            context: &context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: Some((&key, MaterialPropertySelection::SingleClaimOnly)),
    });
    let resistors = [CellResistor {
        destination: DESTINATIONS[0],
        length_m: 1.0,
        area_m2: VOLUME,
        selection: MaterialPropertySelection::SingleClaimOnly,
    }];
    let bulk = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &LINKS,
    };
    let electrical = ElectricalContactNetwork {
        bulk: &bulk,
        contacts: &[],
        electrical_bases: &[(2, ElectricalPropertyBasis::Resistivity)],
    };
    let interface = interface(&cards, resistance);
    let interface_context = QueryPoint::new()
        .with("aT", 300.0)
        .unwrap()
        .with("bT", 300.0)
        .unwrap();
    let contacts = [MaterialContactLink {
        link: 17,
        point: ThermalInterfacePoint {
            interface: &interface,
            materials: [&cards[0], &cards[1]],
            context: &interface_context,
            temperature_axes: ["aT", "bT"],
            liquid_fraction_axes: [None, None],
            mean_temperature_axis: None,
        },
        sides: [0, 1].map(|i| ContactCell {
            cell: cells[i].cell,
            material_id: cells[i].material_id,
            curve: curves[i].curve(),
        }),
        selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
    }];
    let surfaces = [0, 1].map(|i| RadiatingCell {
        surface: i,
        cell: cells[i].cell,
        material_id: cells[i].material_id,
        curve: curves[i].curve(),
        material: cells[i].point,
        selection: MaterialPropertySelection::SingleClaimOnly,
    });
    let enclosure = radiative::enclosure([0, 1].map(|i| {
        SurfaceEmissivity::from_thermal_state(
            "original magnetic assembly surface",
            cells[i].point,
            curves[i].curve().state_at_specific_enthalpy(400.0).unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap()
        .0
    }));
    let radiation = MaterialRadiation {
        enclosure: &enclosure,
        cells: &surfaces,
        power_tolerance_w: 1e-8,
        radiosity_tolerance_w_m2: 1e-8,
    };
    let assembly = MaterialHeatNetwork {
        electrical: &electrical,
        thermal_contacts: &contacts,
        radiation: Some(&radiation),
    };
    let mut events: Vec<_> = [0, 1]
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: cells[i].cell,
                segment_id: i as u64,
                body_id: i as u64,
                material_id: cells[i].material_id,
                mass_kg: VOLUME,
                specific_h_j_kg: 400.0,
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 1.0, 0.0, 0.0],
                    volume_m3: VOLUME,
                    surface_area_m2: 1.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .into();
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 17,
            a: CELL,
            b: SINK,
            area_m2: 1.0,
            conductance_w_k: 999.0,
        }),
    });
    let thermal = DepositionSimulation::new(
        events,
        &[
            (MATERIAL, curves[0].curve() as &dyn ThermalLaw),
            (SINK_MATERIAL, curves[1].curve() as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: DT,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let field = problem(&cards[0]);
    let bindings = laws(&cards[0], &field);
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 91,
                kernel_id: 14,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        run(&cx, &gate, &assembly, thermal, &field, &bindings)
    })
}

fn config() -> CircuitEnthalpySolve {
    CircuitEnthalpySolve {
        maximum_cells: 2,
        ..controls()
    }
}

// Independent manufactured P1 modes: B=mu(Tc)*I, lambda=2*depth*B.
// A scalar sink heat balance and total physical energy determine both T's.
// No native constitutive, radiation, electrical or nonlinear solver is called.
fn reference(t0: [f64; 2], i0: f64, e: f64, contact: f64) -> ([f64; 2], f64) {
    let b0 = mu(t0[0], SLOPE) * i0;
    let residual = |tc: f64, ts: f64| {
        let resistance = [
            1.0 + RHO_SLOPE * (tc - 300.0),
            1.0 + GRAPH_ALPHA * (ts - 300.0),
        ];
        let current = (2.0 * DEPTH * b0 + END_H * i0 + DT * VOLTAGE)
            / (2.0 * DEPTH * mu(tc, SLOPE) + END_H + DT * resistance.iter().sum::<f64>());
        let b = mu(tc, SLOPE) * current;
        let g = 1.0 / (1.0 + contact * (1.0 + 0.001 * (tc - 300.0)));
        let radiated = 5.670_374_419e-8 * (tc.powi(4) - ts.powi(4))
            / (1.0 / (e * (1.0 + 0.0001 * (tc - 300.0)))
                + 1.0 / (0.6 * (1.0 + 0.0001 * (ts - 300.0)))
                - 1.0);
        let transfer = g * (tc - ts) + radiated;
        let numerical = VOLUME * (b - b0).powi(2) / (2.0 * mu(tc, SLOPE))
            + 0.5 * END_H * (current - i0).powi(2);
        (
            HEAT_CAPACITY * (ts - t0[1]) - DT * (resistance[1] * current * current + transfer),
            HEAT_CAPACITY * (tc - t0[0]) + HEAT_CAPACITY * (ts - t0[1]) + field_u(tc, b, SLOPE)
                - field_u(t0[0], b0, SLOPE)
                + 0.5 * END_H * (current * current - i0 * i0)
                + numerical
                - DT * VOLTAGE * current,
            current,
        )
    };
    let at = |tc| {
        let (mut lo, mut hi) = (280.0, 400.0);
        assert!(residual(tc, lo).0 < 0.0 && residual(tc, hi).0 > 0.0);
        for _ in 0..80 {
            let mid = f64::midpoint(lo, hi);
            if residual(tc, mid).0 > 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        let ts = f64::midpoint(lo, hi);
        (ts, residual(tc, ts))
    };
    let (mut lo, mut hi) = (280.0, 360.0);
    assert!(at(lo).1.1 < 0.0 && at(hi).1.1 > 0.0);
    for _ in 0..80 {
        let mid = f64::midpoint(lo, hi);
        if at(mid).1.1 > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let tc = f64::midpoint(lo, hi);
    let (ts, r) = at(tc);
    ([tc, ts], r.2)
}

#[test]
fn material_contact_radiation_and_magnetic_entropy_close_same_voltage_driven_endpoint() {
    let mut outputs = Vec::new();
    for (e, contact) in [(0.2, 0.2), (0.8, 0.2), (0.2, 2.0)] {
        outputs.push(with_heat(
            e,
            contact,
            500.0,
            |cx, _, assembly, mut thermal, field, laws| {
                let mut circuit = circuit(field, 0.0);
                let prepared = assembly.prepare(cx, &thermal).unwrap();
                for n in 1..=2 {
                    let expected = reference(
                        [
                            thermal.temperature_k(CELL).unwrap(),
                            thermal.temperature_k(SINK).unwrap(),
                        ],
                        circuit.currents_a()[0],
                        e,
                        contact,
                    );
                    let next = prepared
                        .step_nonlinear_winding_candidate(
                            cx,
                            &thermal,
                            n as f64 * DT,
                            &circuit,
                            &[VOLTAGE],
                            &routes(),
                            &DESTINATIONS,
                            1e-12,
                            Some((laws, linear::limits())),
                            config(),
                            field::controls(),
                            &[0.0; 2],
                        )
                        .unwrap();
                    for (cell, t) in [CELL, SINK].into_iter().zip(expected.0) {
                        near(next.thermal().temperature_k(cell).unwrap(), t);
                    }
                    near(next.circuit().currents_a()[0], expected.1);
                    assert_eq!(next.report().source.network.thermal_contacts.len(), 1);
                    assert!(next.report().source.network.radiation.is_some());
                    assert_eq!(next.circuit().graph_resistor_materials().len(), 1);
                    assert!(next.caloric_balance().unwrap().external_heat_j.abs() < 1e-8);
                    assert!(
                        next.caloric_balance()
                            .unwrap()
                            .total_energy_residual_j
                            .abs()
                            < config().total_energy_tolerance_j
                    );
                    let parts = next.into_parts();
                    thermal = parts.0;
                    circuit = parts.1;
                }
                [
                    thermal.temperature_k(CELL).unwrap(),
                    thermal.temperature_k(SINK).unwrap(),
                    circuit.currents_a()[0],
                ]
            },
        ));
    }
    for changed in &outputs[1..] {
        for i in 0..3 {
            assert!((changed[i] - outputs[0][i]).abs() > 1e-5);
        }
    }
}

#[test]
fn empty_extensions_preserve_original_magnetic_step_and_cancellation_preserves_history() {
    with_heat(
        0.2,
        0.2,
        500.0,
        |cx, gate, assembly, thermal, field, laws| {
            let circuit = circuit(field, 0.0);
            let bare = MaterialHeatNetwork {
                electrical: assembly.electrical,
                thermal_contacts: &[],
                radiation: None,
            };
            let prepared = bare.prepare(cx, &thermal).unwrap();
            let next = prepared
                .step_nonlinear_winding_candidate(
                    cx,
                    &thermal,
                    DT,
                    &circuit,
                    &[VOLTAGE],
                    &routes(),
                    &DESTINATIONS,
                    1e-12,
                    Some((laws, linear::limits())),
                    config(),
                    field::controls(),
                    &[0.0; 2],
                )
                .unwrap();
            let bulk = assembly.electrical.bulk;
            let original = WindingThermalNetwork {
                cells: bulk.cells,
                links: bulk.links,
                magnetic_regions: &routes(),
                graph_resistors: &DESTINATIONS,
                volume_tolerance_m3: 1e-12,
            }
            .step_nonlinear_thermal_candidate(
                &thermal,
                DT,
                &circuit,
                &[VOLTAGE],
                laws,
                Some((bulk.resistors, assembly.electrical.electrical_bases)),
                config(),
                field::controls(),
                linear::limits(),
                &[0.0; 2],
                || false,
            )
            .unwrap();
            assert_eq!(next.thermal().cells(), original.thermal().cells());
            assert_eq!(next.circuit().currents_a(), original.circuit().currents_a());
            assert_eq!(next.caloric_balance(), original.caloric_balance());
            gate.request();
            assert!(matches!(
                prepared.step_nonlinear_winding_candidate(
                    cx,
                    &thermal,
                    DT,
                    &circuit,
                    &[VOLTAGE],
                    &routes(),
                    &DESTINATIONS,
                    1e-12,
                    Some((laws, linear::limits())),
                    config(),
                    field::controls(),
                    &[0.0; 2]
                ),
                Err(CircuitEnthalpyError::Cancelled)
            ));
            assert_eq!(thermal.time_s(), 0.0);
            assert_eq!(circuit.time_s(), 0.0);
            assert_eq!(circuit.currents_a(), &[0.0]);
        },
    );
}

#[test]
fn missing_endpoint_emissivity_refuses_without_advancing_any_native_history() {
    with_heat(0.2, 0.2, 300.01, |cx, _, assembly, thermal, field, laws| {
        let circuit = circuit(field, 0.0);
        let prepared = assembly.prepare(cx, &thermal).unwrap();
        let before = thermal.cells().clone();
        assert!(
            prepared
                .step_nonlinear_winding_candidate(
                    cx,
                    &thermal,
                    DT,
                    &circuit,
                    &[VOLTAGE],
                    &routes(),
                    &DESTINATIONS,
                    1e-12,
                    Some((laws, linear::limits())),
                    config(),
                    field::controls(),
                    &[0.0; 2]
                )
                .is_err()
        );
        assert_eq!(thermal.cells(), &before);
        assert_eq!(circuit.time_s(), 0.0);
        assert_eq!(circuit.field().time_s(), 0.0);
        assert_eq!(circuit.currents_a(), &[0.0]);
    });
}
