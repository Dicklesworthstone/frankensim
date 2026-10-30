//! G1/G3/G4: complete heat transfer and nonlinear material field storage.
//! Manufactured sources verify coupled equations, not empirical material data.
use super::*;
use crate::{
    REGIONS,
    heat::{interface as thermal_interface, radiative},
};
use fs_conduction::{
    EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY, SurfaceEmissivity,
    transient::deposition::{contact::*, radiation::*},
};
use fs_couple::circuit_enthalpy::{
    conduction::contacts::ConductionContactThermal, contacts::heat::PreparedHeatNetwork,
    material::ConductionLink,
};
use fs_matdb::ClaimSelection;
use fs_material::thermal::{
    conduction::{
        LayerConductionConfig,
        network::{LayerLinkBinding, PreparedLayerConduction},
        source::ResolvedConductionLayer,
    },
    interface::ThermalInterfacePoint,
};

const FIELD_SLOPE: f64 = 0.001;
const K_SLOPE: f64 = 0.002;
const DT: f64 = 0.1;
const LINKS: [ConductionLink; 2] = [
    ConductionLink {
        link: 7,
        distance_a_m: 0.5,
        distance_b_m: 0.5,
    },
    ConductionLink {
        link: 8,
        distance_a_m: 0.5,
        distance_b_m: 0.5,
    },
];

fn layer_controls() -> LayerConductionConfig {
    LayerConductionConfig {
        max_iterations: 80,
        absolute_heat_tolerance_w: 1e-11,
        relative_heat_tolerance: 1e-12,
    }
}

fn conductivity_claim(scale: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
        value: PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(200.0, 0.8 * scale), (2000.0, 4.4 * scale)],
            dims: THERMAL_CONDUCTIVITY_DIMS,
        },
        validity: ValidityDomain::unconstrained().with("T", 200.0, 2000.0),
        interpolation: InterpolationPolicy::LinearInside,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "manufactured affine conductivity".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}

// Every changed card is rebound to its original caloric source. The electrical
// and field-storage laws themselves remain the parent's original source data.
fn with_transport(
    original: MaterialCard,
    scale: f64,
    surface: Option<(f64, f64, bool)>,
) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (_, claim) in original.claims().claims_ordered() {
        if claim.key.name() != "thermal_conductivity" {
            claims.insert_claim(claim.clone()).unwrap();
        }
    }
    claims.insert_claim(conductivity_claim(scale)).unwrap();
    if let Some((emissivity, upper, field_dependent)) = surface {
        let mut validity = ValidityDomain::unconstrained().with("T", 200.0, upper);
        if field_dependent {
            validity = validity.with_quantity(
                "operating",
                ConductorMagnitude::ElectricField.quantity(),
                0.0,
                20.0,
            );
        }
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
                value: PropertyValue::Scalar {
                    value: emissivity,
                    dims: EMISSIVITY_DIMS,
                },
                validity,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured gray surface".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap()
}

fn caloric(
    card: &MaterialCard,
    lower: &QueryPoint,
    upper: &QueryPoint,
) -> ResolvedHeatCapacityCurve {
    ResolvedHeatCapacityCurve::try_new(
        card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower,
            upper,
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-9,
    )
    .unwrap()
}

#[allow(clippy::too_many_arguments)]
fn with_heat<R>(
    magnetic: bool,
    emissivity: f64,
    contact_resistance: f64,
    conductivity_scale: f64,
    emissivity_upper: f64,
    field_dependent_surface: bool,
    run: impl FnOnce(
        &Cx<'_>,
        &CancelGate,
        &PreparedHeatNetwork<'_, '_>,
        DepositionSimulation<'_>,
        ConductionResponseCircuit,
        &[ConductionContactThermal<'_>],
        &[FieldThermalBinding],
    ) -> R,
) -> R {
    let coordinate = ConductorMagnitude::ElectricField;
    let cards = [0, 1].map(|i| {
        with_transport(
            material((i + 1) as f64, CP[i], ALPHA[i], BETA[i], coordinate, false).0,
            conductivity_scale,
            Some((
                if i == 0 { emissivity } else { 0.6 },
                if i == 0 { emissivity_upper } else { 2000.0 },
                i == 0 && field_dependent_surface,
            )),
        )
    });
    let curves = cards
        .each_ref()
        .map(|card| caloric(card, &point(200.0, coordinate), &point(2000.0, coordinate)));
    let field_card = with_transport(
        field_card(magnetic, FIELD_SLOPE).0,
        conductivity_scale,
        None,
    );
    let field_curve = caloric(&field_card, &field_point(250.0), &field_point(500.0));
    let electrical_interface = interface(&cards, 0.5, CONTACT_SLOPE, 1000.0);
    let context = point(300.0, coordinate);
    let contact_context = interface_context(300.0);
    let field_context = field_point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let mut cells = cells(&cards, &curves, &context, &key);
    cells.push(MaterialCell {
        cell: FIELD_CELL,
        material_id: MaterialId(2),
        curve: field_curve.curve(),
        point: ThermalMaterialPoint {
            material: &field_card,
            context: &field_context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: Some((&key, MaterialPropertySelection::SingleClaimOnly)),
    });
    let resistors = [resistor()];
    let bulk = MaterialNetwork {
        cells: &cells,
        resistors: &resistors,
        links: &LINKS,
    };
    let electrical = ElectricalContactNetwork {
        bulk: &bulk,
        contacts: &[],
        electrical_bases: &[],
    };
    let contacts = [binding(
        &cards,
        &electrical_interface,
        &contact_context,
        FRACTION,
    )];
    let thermal = history_with_links(
        &curves,
        &field_curve,
        1.0,
        &[ThermalLink {
            id: 8,
            a: CellId(30),
            b: FIELD_CELL,
            area_m2: 1.0,
            conductance_w_k: 999.0,
        }],
    );
    let (circuit, field) = initial_circuit(
        specimen(&cards, &electrical_interface, coordinate),
        &bulk,
        &thermal,
        &field_card,
        magnetic,
    );
    let fields = [field];
    let thermal_interface = thermal_interface(&cards, contact_resistance);
    let thermal_context = QueryPoint::new()
        .with("aT", 300.0)
        .unwrap()
        .with("bT", 300.0)
        .unwrap()
        .with("pressure", 100000.0)
        .unwrap();
    let thermal_contacts = [MaterialContactLink {
        link: 7,
        point: ThermalInterfacePoint {
            interface: &thermal_interface,
            materials: [&cards[0], &cards[1]],
            context: &thermal_context,
            temperature_axes: ["aT", "bT"],
            liquid_fraction_axes: [None, None],
            mean_temperature_axis: None,
        },
        sides: [0, 1].map(|i| ContactCell {
            cell: cells[i].cell,
            material_id: cells[i].material_id,
            curve: cells[i].curve,
        }),
        selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
    }];
    let surfaces = [0, 1].map(|i| RadiatingCell {
        surface: i,
        cell: cells[i].cell,
        material_id: cells[i].material_id,
        curve: cells[i].curve,
        material: cells[i].point,
        selection: MaterialPropertySelection::SingleClaimOnly,
    });
    let emissivities = [0, 1].map(|i| {
        SurfaceEmissivity::from_thermal_state(
            "original source surface",
            cells[i].point,
            curves[i]
                .curve()
                .state_at_specific_enthalpy(100.0 * CP[i])
                .unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap()
        .0
    });
    let enclosure = radiative::enclosure(emissivities);
    let radiation = MaterialRadiation {
        enclosure: &enclosure,
        cells: &surfaces,
        power_tolerance_w: 1e-9,
        radiosity_tolerance_w_m2: 1e-9,
    };
    let network = MaterialHeatNetwork {
        electrical: &electrical,
        thermal_contacts: &thermal_contacts,
        radiation: Some(&radiation),
    };
    let layers = [0, 2].map(|i| {
        ResolvedConductionLayer::new(
            cells[i].curve,
            cells[i].point,
            &key,
            ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
            SolidLiquidPhase::Solid,
            [250.0, 500.0],
            || Ok(()),
        )
        .unwrap()
    });
    let integrated = PreparedLayerConduction::new(
        &thermal,
        vec![LayerLinkBinding {
            link: 8,
            layers: [&layers[0], &layers[1]],
            distances_m: [0.5; 2],
            contact_area_resistance_m2_k_w: 0.0,
        }],
        1,
        || Ok(()),
    )
    .unwrap();
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 73,
                kernel_id: 57,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        let prepared = network
            .prepare(&cx, &thermal)
            .unwrap()
            .with_integrated_conduction(&cx, &thermal, &integrated, layer_controls())
            .unwrap();
        run(&cx, &gate, &prepared, thermal, circuit, &contacts, &fields)
    })
}

#[test]
fn source_storage_and_every_heat_path_close_independent_continuing_balances() {
    for magnetic in [false, true] {
        let mut outcomes = Vec::new();
        for (emissivity, contact_resistance, conductivity_scale) in [
            (0.2, 0.2, 1.0),
            (0.8, 0.2, 1.0),
            (0.2, 2.0, 1.0),
            (0.2, 0.2, 4.0),
        ] {
            let outcome = with_heat(
                magnetic,
                emissivity,
                contact_resistance,
                conductivity_scale,
                2000.0,
                false,
                |cx, _, prepared, mut thermal, mut circuit, contacts, fields| {
                    let (mut old, mut z0, mut final_current) = ([300.0; 3], 3.0, 0.0);
                    let initial = thermal.energy_j() + internal(FIELD_SLOPE, old[2], z0);
                    for step in 1..=3 {
                        let next = prepared
                            .step_conduction_response_field_circuit_candidate(
                                cx,
                                &thermal,
                                step as f64 * DT,
                                &circuit,
                                &[],
                                contacts,
                                fields,
                                &REGIONS,
                                1e-12,
                                limits(),
                                field_controls(),
                                controls(),
                                config(),
                                &[0.0; 3],
                            )
                            .unwrap();
                        let t = [CellId(30), CellId(10), FIELD_CELL]
                            .map(|id| next.thermal().temperature_k(id).unwrap());
                        let theta = t.map(|value| value - 300.0);
                        let z = if magnetic {
                            next.circuit().fluxes_wb()[0]
                        } else {
                            next.circuit().charges_c()[0]
                        };
                        let report = next.report();
                        let record = &report.coupled.electrical.circuit;
                        let heat = &report.coupled.source.network;
                        let field = &report.coupled.source.fields[0];
                        let j = record.spatial_current_a;
                        final_current = j;
                        // Independent exact algebra for the manufactured source
                        // laws, not the solver's fitted operator or heat maps.
                        let e = [0, 1].map(|i| {
                            j * ((i + 1) as f64 + ALPHA[i] * theta[i]) / (1.0 - BETA[i] * j)
                        });
                        let rc = 0.5 * (1.0 + CONTACT_SLOPE * (theta[0] + theta[1]) / 2.0);
                        let v = e[0] + e[1] + j * rc;
                        let load = LOAD_R + LOAD_SLOPE * theta[2];
                        let current = j + v / load;
                        let coefficient = 1.0 + FIELD_SLOPE * theta[2];
                        near(record.spatial_voltage_v, v);
                        near(z, z0 - DT * if magnetic { v } else { current });
                        near(
                            0.5 * (z0 + z),
                            coefficient * if magnetic { current } else { v },
                        );
                        near(
                            next.circuit().material_resistors()[0]
                                .conductor
                                .resistance_ohm(),
                            load,
                        );
                        near(field.material.coefficient(), coefficient);
                        near(record.spatial.contacts[0].response.power_w, j * j * rc);
                        near(record.graph_resistors[0].heat_j, DT * v * v / load);

                        let k = theta.map(|x| conductivity_scale * (1.0 + K_SLOPE * x));
                        let rt = contact_resistance * (1.0 + 0.001 * theta[0]);
                        let g01 = 1.0 / (0.5 / k[0] + 0.5 / k[1] + rt);
                        let radiation = 5.670_374_419e-8 * (t[0].powi(4) - t[1].powi(4))
                            / (1.0 / emissivity + 1.0 / 0.6 - 1.0);
                        let q01 = g01 * (t[0] - t[1]) + radiation;
                        // Equal affine source laws in the two half-layers:
                        // integral k(T)dT / deltaT is their arithmetic mean.
                        let g02 =
                            conductivity_scale * (1.0 + K_SLOPE * (theta[0] + theta[2]) / 2.0);
                        let q02 = g02 * (t[0] - t[2]);
                        near(
                            CP[0] * (t[0] - old[0]),
                            DT * (j * e[0] + FRACTION * j * j * rc - q01 - q02),
                        );
                        near(
                            CP[1] * (t[1] - old[1]),
                            DT * (j * e[1] + (1.0 - FRACTION) * j * j * rc + q01),
                        );

                        let delta_u =
                            internal(FIELD_SLOPE, t[2], z) - internal(FIELD_SLOPE, old[2], z0);
                        // Circuit work into the finite field cell is -v*I.
                        // Its physical U and caloric energy include the load,
                        // actual layer flux and only the authored 1 W source.
                        near(
                            FIELD_CP * (t[2] - old[2]) + delta_u,
                            DT * (-v * current + v * v / load + q02 + 1.0),
                        );
                        near(
                            FIELD_CP * (t[2] - old[2]),
                            DT * (v * v / load + q02 + 1.0) + field.reversible_heat_j,
                        );
                        near(
                            field.state[0].internal_energy_j,
                            internal(FIELD_SLOPE, old[2], z0),
                        );
                        near(
                            field.state[1].internal_energy_j,
                            internal(FIELD_SLOPE, t[2], z),
                        );
                        near(
                            next.thermal().energy_j() + internal(FIELD_SLOPE, t[2], z),
                            initial + step as f64 * DT,
                        );
                        near(report.external_heat_j, DT);
                        near(report.total_energy_residual_j, 0.0);
                        near(
                            report.coupled.deposited_joule_j,
                            DT * (j * v + v * v / load),
                        );
                        assert_eq!(heat.thermal_contacts.len(), 1);
                        near(heat.thermal_contacts[0].interface.value_m2_k_per_w(), rt);
                        assert_eq!(heat.integrated_conduction.len(), 1);
                        let layer = &heat.integrated_conduction[0];
                        assert_eq!(layer.link, 8);
                        near(layer.response.conductance_w_k, g02);
                        near(layer.applied_power_a_to_b_w, q02);
                        assert!(heat.radiation.is_some());
                        near(heat.radiation_net_power_w(), 0.0);
                        for surface in heat.radiation.as_ref().unwrap().surfaces() {
                            let cell = heat
                                .material
                                .cells
                                .iter()
                                .find(|c| c.cell == surface.cell)
                                .unwrap();
                            assert_eq!(cell.phase, surface.phase);
                        }
                        assert!(report.coupled.electrical.field_solves > 1);
                        assert_eq!(next.circuit().storage().len(), 1);
                        assert_eq!(circuit.time_s(), (step - 1) as f64 * DT);
                        old = t;
                        z0 = z;
                        (thermal, circuit, _) = next.into_parts();
                    }
                    [old[0], old[1], old[2], final_current]
                },
            );
            outcomes.push(outcome);
        }
        for changed in &outcomes[1..] {
            for (actual, original) in changed.iter().zip(outcomes[0]) {
                assert!((actual - original).abs() > 1e-9);
            }
        }
    }
}

#[test]
fn full_heat_field_candidate_retains_sources_and_refuses_missing_work_atomically() {
    with_heat(
        false,
        0.2,
        0.2,
        1.0,
        2000.0,
        false,
        |cx, gate, prepared, thermal, circuit, contacts, fields| {
            let before = (
                thermal.energy_j(),
                circuit.graph().clone(),
                circuit.conductors().to_vec(),
                circuit.contacts().to_vec(),
                circuit.storage().to_vec(),
            );
            let mut exhausted = config();
            exhausted.maximum_circuit_trials = 1;
            assert!(matches!(
                prepared.step_conduction_response_field_circuit_candidate(
                    cx,
                    &thermal,
                    DT,
                    &circuit,
                    &[],
                    contacts,
                    fields,
                    &REGIONS,
                    1e-12,
                    limits(),
                    field_controls(),
                    controls(),
                    exhausted,
                    &[0.0; 3],
                ),
                Err(CircuitEnthalpyError::CircuitBudget)
            ));
            for (contacts, fields) in [(&[][..], fields), (contacts, &[][..])] {
                assert!(
                    prepared
                        .step_conduction_response_field_circuit_candidate(
                            cx,
                            &thermal,
                            DT,
                            &circuit,
                            &[],
                            contacts,
                            fields,
                            &REGIONS,
                            1e-12,
                            limits(),
                            field_controls(),
                            controls(),
                            config(),
                            &[0.0; 3],
                        )
                        .is_err()
                );
            }
            // The previous ideal-storage heat route must not drop a graph source.
            assert!(
                prepared
                    .step_conduction_response_contact_circuit_candidate(
                        cx,
                        &thermal,
                        DT,
                        &circuit,
                        &[],
                        &[resistor().destination],
                        contacts,
                        &REGIONS,
                        1e-12,
                        limits(),
                        field_controls(),
                        controls(),
                        config(),
                        &[0.0; 3],
                    )
                    .is_err()
            );
            let accepted = prepared
                .step_conduction_response_field_circuit_candidate(
                    cx,
                    &thermal,
                    DT,
                    &circuit,
                    &[],
                    contacts,
                    fields,
                    &REGIONS,
                    1e-12,
                    limits(),
                    field_controls(),
                    controls(),
                    config(),
                    &[0.0; 3],
                )
                .unwrap();
            assert!(accepted.thermal().temperature_k(FIELD_CELL).unwrap() > 300.0);
            gate.request();
            assert!(matches!(
                prepared.step_conduction_response_field_circuit_candidate(
                    cx,
                    &thermal,
                    DT,
                    &circuit,
                    &[],
                    contacts,
                    fields,
                    &REGIONS,
                    1e-12,
                    limits(),
                    field_controls(),
                    controls(),
                    config(),
                    &[0.0; 3],
                ),
                Err(CircuitEnthalpyError::Cancelled)
            ));
            assert_eq!(thermal.energy_j(), before.0);
            assert_eq!(circuit.graph(), &before.1);
            assert_eq!(circuit.conductors(), before.2);
            assert_eq!(circuit.contacts(), before.3);
            assert_eq!(circuit.storage(), before.4);
            assert_eq!(circuit.charges_c(), &[3.0]);
            assert_eq!(circuit.time_s(), 0.0);
            assert_eq!(thermal.time_s(), 0.0);
        },
    );
    for (upper, dependent) in [(300.01, false), (2000.0, true)] {
        with_heat(
            false,
            0.2,
            0.2,
            1.0,
            upper,
            dependent,
            |cx, _, prepared, thermal, circuit, contacts, fields| {
                assert!(
                    prepared
                        .step_conduction_response_field_circuit_candidate(
                            cx,
                            &thermal,
                            DT,
                            &circuit,
                            &[],
                            contacts,
                            fields,
                            &REGIONS,
                            1e-12,
                            limits(),
                            field_controls(),
                            controls(),
                            config(),
                            &[0.0; 3],
                        )
                        .is_err()
                );
                assert_eq!(thermal.time_s(), 0.0);
                assert_eq!(circuit.time_s(), 0.0);
                assert_eq!(circuit.charges_c(), &[3.0]);
                assert_eq!(circuit.storage().len(), 1);
            },
        );
    }
}
