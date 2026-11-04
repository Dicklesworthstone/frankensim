//! Two actual walls, native circuit loss and radiation share one gas-gap link.
use super::*;
use fs_airflow::graph::thermal::medium::convection::deposition::source::{
    enclosure::{
        SourceEnclosedGasLayer, evaluate_source_enclosed_gas_layers,
        heat::SourceEnclosedGasHeatAssembly,
    },
    heat::composition::SourceConvectionHeatLayer,
};
use fs_conduction::transient::deposition::{
    contact::{ContactCell, MaterialContactLink},
    radiation::{MaterialRadiation, RadiatingCell},
};
use fs_conduction::{
    EMISSIVITY_DIMS, GrayDiffuseEnclosure, RadiationSurface, SURFACE_EMISSIVITY_PROPERTY,
    SurfaceEmissivity, ViewFactorMatrix,
};
use fs_couple::circuit_enthalpy::{
    CircuitEnthalpyInput, contacts,
    contacts::heat::{HeatAssembly, ThermalLinkMetric, observe_heat_assembly},
    material::ConductionLink,
};
use fs_exec::{CancelGate, Cx};
use fs_matdb::{ClaimSelection, InterfaceSystemCard, SurfaceSpec, SystemContext};
use fs_material::thermal::interface::ThermalInterfacePoint;
use fs_material::thermochemical_gas::ThermochemicalGasModel;

#[path = "../../../finite_gas_convection/source.rs"]
mod gas_source;

const UPPER: CellId = CellId(72);
const WALLS: [CellId; 2] = [BODY, UPPER];
const GAP_LINK: u64 = 23;
const GAP: f64 = 0.018;
const SPANS: [f64; 2] = [0.2, 0.3];
const GAP_AREA: f64 = SPANS[0] * SPANS[1];
const EMISSIVITY: f64 = 0.6;

fn gap() -> SourceEnclosedGasLayer {
    SourceEnclosedGasLayer {
        link: GAP_LINK,
        lower_cell: BODY,
        gap: Length::new(GAP),
        plate_spans: SPANS.map(Length::new),
        gravity_m_s2: 9.81,
        inclination_degrees: 0.0,
        maximum_relative_density_departure: 0.1,
    }
}

fn wall_material() -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let (original, _) = body_material();
    let mut claims = original.claims().clone();
    for (name, dims, value) in [
        (SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS, EMISSIVITY),
        (
            THERMAL_CONDUCTIVITY_PROPERTY,
            THERMAL_CONDUCTIVITY_DIMS,
            1.0,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained()
                    .with_quantity("T", temperature_quantity(), 280.0, 320.0)
                    .with("lot", 7.0, 7.0),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured isothermal gap-wall fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new(
                SPECIFIC_HEAT_CAPACITY_PROPERTY,
                SPECIFIC_HEAT_CAPACITY_DIMS,
            ),
            density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            temperature_axis: "T",
            lower: &body_point(280.0),
            upper: &body_point(320.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-6,
    )
    .unwrap();
    (card, curve)
}

fn with_case<T>(
    electrical: bool,
    radiating: bool,
    run: impl FnOnce(
        &Cx<'_>,
        &CancelGate,
        &MaterialHeatNetwork<'_>,
        DepositionSimulation<'_>,
        &ThermochemicalGasModel,
    ) -> T,
) -> T {
    let (card, curve) = wall_material();
    let contexts = [body_point(310.0), body_point(290.0)];
    let conductivity = PropertyKey::new(THERMAL_CONDUCTIVITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS);
    let cells = [0, 1].map(|i| MaterialCell {
        cell: WALLS[i],
        material_id: MaterialId(0),
        curve: curve.curve(),
        point: ThermalMaterialPoint {
            material: &card,
            context: &contexts[i],
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: Some((&conductivity, MaterialPropertySelection::SingleClaimOnly)),
    });
    let resistors = [CellResistor {
        destination: ResistorCell {
            branch: 1,
            cell: BODY,
        },
        selection: MaterialPropertySelection::SingleClaimOnly,
        length_m: 0.1,
        area_m2: 1e-5,
    }];
    let bulk = MaterialNetwork {
        cells: &cells,
        resistors: if electrical { &resistors } else { &[] },
        links: &[],
    };
    let electrical = ElectricalContactNetwork {
        bulk: &bulk,
        contacts: &[],
        electrical_bases: if electrical {
            &[(1, ElectricalPropertyBasis::Resistivity)]
        } else {
            &[]
        },
    };
    let mut events = [310.0, 290.0]
        .into_iter()
        .enumerate()
        .map(|(i, t)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: WALLS[i],
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MaterialId(0),
                mass_kg: MASS,
                specific_h_j_kg: CP * (t - 280.0),
                geometry: CellGeometry {
                    start_m: [0.0, 0.0, i as f64 * GAP],
                    end_m: [SPANS[0], 0.0, i as f64 * GAP],
                    volume_m3: 1e-6,
                    surface_area_m2: 2.0 * GAP_AREA,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    source_w: 0.0,
                    reservoir: None,
                },
            }),
        })
        .collect::<Vec<_>>();
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: GAP_LINK,
            a: UPPER,
            b: BODY,
            area_m2: GAP_AREA,
            conductance_w_k: 0.0,
        }),
    });
    let thermal = DepositionSimulation::new(
        events,
        &[(MaterialId(0), curve.curve() as &dyn ThermalLaw)],
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
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let (complex, positions) =
        fs_conduction::fixtures::box_grid([1, 1, 1], [SPANS[0], SPANS[1], GAP]);
    let mesh = fs_conduction::mesh::ConductionMesh::new(complex, positions).unwrap();
    let surfaces = [0, 1].map(|i| {
        let emissivity = SurfaceEmissivity::from_thermal_state(
            "gap wall",
            cells[i].point,
            curve
                .curve()
                .state_at_specific_enthalpy(thermal.cells()[&WALLS[i]].specific_h_j_kg)
                .unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap()
        .0;
        RadiationSurface::new(
            &mesh,
            format!("gap-wall-{i}"),
            |face| fs_conduction::fixtures::on_box_face(face.centroid[2], i as f64 * GAP),
            emissivity,
        )
        .unwrap()
    });
    // Analytic parallel-plate radiation limit; not a ray model of finite edges.
    let enclosure = GrayDiffuseEnclosure::new(
        surfaces.to_vec(),
        ViewFactorMatrix::infinite_parallel_plates(GAP_AREA).unwrap(),
    )
    .unwrap();
    let radiating_cells = [0, 1].map(|i| RadiatingCell {
        surface: i,
        cell: WALLS[i],
        material_id: MaterialId(0),
        curve: curve.curve(),
        material: cells[i].point,
        selection: MaterialPropertySelection::SingleClaimOnly,
    });
    let radiation = MaterialRadiation {
        enclosure: &enclosure,
        cells: &radiating_cells,
        power_tolerance_w: 1e-9,
        radiosity_tolerance_w_m2: 1e-9,
    };
    let heat = MaterialHeatNetwork {
        electrical: &electrical,
        thermal_contacts: &[],
        radiation: radiating.then_some(&radiation),
    };
    // Original source packs, deliberately manufactured: Pr ~0.7195 exercises
    // the explicit air-policy API without claiming empirical fluid validation.
    let gas = gas_source::model(1.0);
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 229,
                kernel_id: 29,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        run(&cx, &gate, &heat, thermal, &gas)
    })
}

// Independent gas-layer algebra and Stefan-Boltzmann parallel-plate balance.
// The gas properties still come from the original manufactured material law.
fn gap_conductance(
    gas: &ThermochemicalGasModel,
    pressure: f64,
    lower: f64,
    upper: f64,
    inclination_degrees: f64,
) -> f64 {
    let sample = gas
        .evaluate(&gas_source::point(0.5 * lower + 0.5 * upper, pressure))
        .unwrap();
    let s = sample.state();
    let ra = 9.81 * (lower - upper) / s.temperature
        * GAP.powi(3)
        * s.density.powi(2)
        * s.specific_heat_cp
        / (s.dynamic_viscosity * s.thermal_conductivity);
    let angle = inclination_degrees.to_radians();
    let normal_ra = ra * angle.cos();
    let nu = if inclination_degrees == 90.0 {
        let magnitude = ra.abs();
        let tall = if magnitude <= 1e4 {
            1.0 + 1.75967e-10 * magnitude.powf(2.2984755)
        } else if magnitude <= 5e4 {
            0.028154 * magnitude.powf(0.4134)
        } else {
            0.0673838 * magnitude.cbrt()
        };
        tall.max(0.242 * (magnitude / (SPANS[1] / GAP)).powf(0.272))
    } else if normal_ra <= 1708.0 {
        1.0
    } else {
        1.0 + 1.44
            * (1.0 - 1708.0 / normal_ra)
            * (1.0 - 1708.0 * (1.8 * angle).sin().powf(1.6) / normal_ra)
            + ((normal_ra / 5830.0).cbrt() - 1.0).max(0.0)
    };
    GAP_AREA * s.thermal_conductivity / GAP * nu
}

fn reference(gas: &ThermochemicalGasModel, pressure: f64, inclination_degrees: f64) -> [f64; 3] {
    let evaluate = |lower: f64| {
        let resistance = 1.0 + 0.02 * (lower - 300.0);
        let ratio = DT / (2.0 * CAPACITANCE * resistance);
        let charge = INITIAL_CHARGE * (1.0 - ratio) / (1.0 + ratio);
        let joule = (INITIAL_CHARGE.powi(2) - charge.powi(2)) / (2.0 * CAPACITANCE);
        let upper = 290.0 + joule / (MASS * CP) - (lower - 310.0);
        let gas_power =
            gap_conductance(gas, pressure, lower, upper, inclination_degrees) * (lower - upper);
        let radiation = 5.670_374_419e-8 * GAP_AREA * (lower.powi(4) - upper.powi(4))
            / (2.0 / EMISSIVITY - 1.0);
        (
            MASS * CP * (lower - 310.0) - joule + DT * (gas_power + radiation),
            [lower, upper, charge],
        )
    };
    let (mut lower, mut upper) = (300.0, 319.0);
    assert!(evaluate(lower).0 < 0.0 && evaluate(upper).0 > 0.0);
    for _ in 0..80 {
        let midpoint = 0.5 * lower + 0.5 * upper;
        if evaluate(midpoint).0 > 0.0 {
            upper = midpoint;
        } else {
            lower = midpoint;
        }
    }
    let (residual, state) = evaluate(0.5 * lower + 0.5 * upper);
    assert!(residual.abs() < 1e-10);
    state
}

#[test]
fn gas_gap_changes_actual_circuit_and_wall_storage_with_closed_radiation() {
    for inclination_degrees in [0.0, 45.0, 90.0] {
        with_case(true, true, |cx, _, heat, thermal, gas| {
            let saved = thermal.cells().clone();
            let layers = [SourceEnclosedGasLayer {
                inclination_degrees,
                ..gap()
            }];
            let mut results = Vec::new();
            for pressure in [50000.0, 100000.0] {
                let resolve =
                    |_: &SourceEnclosedGasLayer, _: &DepositionSimulation<'_>, t: Temperature| {
                        gas.evaluate(&gas_source::point(t.value(), pressure))
                    };
                let initial =
                    evaluate_source_enclosed_gas_layers(&thermal, &layers, 1, resolve, || false)
                        .unwrap();
                assert!(heat.prepare(cx, &thermal).is_err());
                let prepared = heat
                    .prepare_with_link_conductances(cx, &thermal, initial.conductances_w_k())
                    .unwrap();
                let assembly = SourceEnclosedGasHeatAssembly {
                    heat: &prepared,
                    layers: &layers,
                    maximum_layers: 1,
                    resolve_gas: resolve,
                };
                let coefficients = assembly.resolve_heat(cx, &thermal).unwrap();
                let circuit = MaterialCircuit::new(
                    CircuitGraph {
                        node_count: 2,
                        branches: vec![
                            (
                                1,
                                0,
                                Branch::Capacitor {
                                    farads: CAPACITANCE,
                                },
                            ),
                            (1, 0, Branch::Resistor { ohms: 0.0 }),
                        ],
                        transformers: vec![],
                    },
                    coefficients.bulk.resistors,
                    vec![],
                    vec![INITIAL_CHARGE],
                    &[],
                )
                .unwrap();
                let destinations = prepared.resistor_destinations();
                let routes = prepared.heat_destinations(&thermal).unwrap();
                let mut solve = controls();
                solve.maximum_cells = 2;
                solve.newton.linear_restart = 2;
                // The existing contact/circuit driver consumes this complete
                // HeatAssembly resolver; no gap-specific evolution exists.
                let next = contacts::step_candidate(
                    &thermal,
                    DT,
                    CircuitEnthalpyInput {
                        circuit: &circuit,
                        ports: &[],
                        destinations: &destinations,
                    },
                    &routes,
                    solve,
                    &[0.0; 2],
                    |state| assembly.resolve_heat(cx, state),
                    || false,
                )
                .unwrap();
                let actual = [
                    next.thermal().temperature_k(BODY).unwrap(),
                    next.thermal().temperature_k(UPPER).unwrap(),
                    next.circuit().charges_c()[0],
                ];
                for (a, b) in actual
                    .into_iter()
                    .zip(reference(gas, pressure, inclination_degrees))
                {
                    near(a, b, 1e-8);
                }
                let sample = &next.report().source.material;
                let gap = &sample.gas_layers.samples()[0];
                assert_eq!(gap.inclination_degrees, inclination_degrees);
                near(gap.lower_temperature_k, actual[0], 1e-12);
                near(gap.upper_temperature_k, actual[1], 1e-12);
                near(
                    gap.property_temperature_k,
                    0.5 * actual[0] + 0.5 * actual[1],
                    1e-12,
                );
                near(gap.material.state().pressure, pressure, 1e-8);
                assert_eq!(
                    gap.material.source_packs(),
                    initial.samples()[0].material.source_packs()
                );
                near(
                    gap.conductance_w_k,
                    gap_conductance(gas, pressure, actual[0], actual[1], inclination_degrees),
                    1e-12,
                );
                assert!(gap.heat_into_upper_w > 0.0);
                let radiation = sample.heat.radiation.as_ref().unwrap();
                assert!(radiation.sources_w()[&BODY] < 0.0 && radiation.sources_w()[&UPPER] > 0.0);
                let external = next.report().non_electrical_heat_j - DT * radiation.net_heat_w();
                near(external, 0.0, 1e-10);
                near(
                    next.thermal().energy_j() - thermal.energy_j()
                        + (actual[2].powi(2) - INITIAL_CHARGE.powi(2)) / (2.0 * CAPACITANCE),
                    0.0,
                    1e-8,
                );
                assert_eq!(thermal.cells(), &saved);
                assert_eq!(circuit.charges_c(), &[INITIAL_CHARGE]);
                assert_eq!(next.thermal().links()[&GAP_LINK].conductance_w_k, 0.0);
                results.push(actual);
            }
            assert!((results[0][0] - results[1][0]).abs() > 1e-3);
            assert!((results[0][2] - results[1][2]).abs() > 1e-6);
        });
    }
}

#[test]
fn outside_reservoir_composes_with_gap_radiation_and_native_explicit_bounds() {
    with_case(false, true, |cx, _, heat, thermal, gas| {
        let layers = [gap()];
        let resolve = |_: &SourceEnclosedGasLayer, _: &DepositionSimulation<'_>, t: Temperature| {
            gas.evaluate(&gas_source::point(t.value(), 100000.0))
        };
        let initial =
            evaluate_source_enclosed_gas_layers(&thermal, &layers, 1, resolve, || false).unwrap();
        let prepared = heat
            .prepare_with_link_conductances(cx, &thermal, initial.conductances_w_k())
            .unwrap();
        let gaps = SourceEnclosedGasHeatAssembly {
            heat: &prepared,
            layers: &layers,
            maximum_layers: 1,
            resolve_gas: resolve,
        };
        let reservoirs = [SourceConvectiveCell {
            cell: BODY,
            reference_temperature: Temperature::new(300.0),
            characteristic_length: Length::new(0.1),
            property_temperature: ConvectionPropertyTemperature::FilmMean,
            law: SourceCellConvectionLaw::Forced {
                correlation: CorrelationId::FlatPlateLaminarAverage,
                speed_m_s: 1.0,
                settings: ForcedConvectionSettings::default(),
            },
        }];
        let assembly = SourceConvectionHeatLayer {
            inner: &gaps,
            cells: &reservoirs,
            maximum_cells: 1,
            resolve_fluid: |_: &SourceConvectiveCell,
                            _: &DepositionSimulation<'_>,
                            t: Temperature| {
                gas.evaluate(&gas_source::point(t.value(), 100000.0))
            },
        };
        let result = assembly.resolve_heat(cx, &thermal).unwrap();
        let sample = &result.bulk.source;
        let bath = &sample.convection.convection().samples()[0];
        near(bath.exposed_area_m2, GAP_AREA, 1e-15);
        near(
            result.bulk.conductances_w_k[&GAP_LINK],
            initial.conductances_w_k()[&GAP_LINK],
            1e-15,
        );
        near(
            result.bulk.sources_w.values().sum::<f64>(),
            bath.heat_into_cell_w + sample.heat.heat.radiation_net_power_w(),
            1e-12,
        );
        let radiation = sample
            .heat
            .heat
            .explicit_source_conductance_bounds_w_k(&mut || Ok(()))
            .unwrap();
        let bounds = assembly
            .explicit_source_conductance_bounds_w_k(sample, &mut || Ok(()))
            .unwrap();
        let reservoir_g = bath.evaluation.as_ref().unwrap().coefficient.value.value() * GAP_AREA;
        near(bounds[&BODY], radiation[&BODY] + reservoir_g, 1e-12);
        near(bounds[&UPPER], radiation[&UPPER], 1e-12);
        let observer = observe_heat_assembly(&assembly, cx, (), 1e-8).unwrap();
        let (next, accepted, audit) = thermal
            .advance_observed(
                &observer,
                0.5,
                StepConfig {
                    max_dt_s: 0.1,
                    energy_tolerance_j: 1e-8,
                },
                || Ok(()),
            )
            .unwrap();
        let endpoint = accepted.first.sample().unwrap();
        let gap = &endpoint.heat.gas_layers.samples()[0];
        near(
            gap.lower_temperature_k,
            next.temperature_k(BODY).unwrap(),
            1e-12,
        );
        near(
            gap.upper_temperature_k,
            next.temperature_k(UPPER).unwrap(),
            1e-12,
        );
        assert!(audit.external_heat_j - accepted.first.closed_radiation_net_j() < 0.0);
        near(
            next.energy_j() - thermal.energy_j(),
            audit.external_heat_j,
            1e-8,
        );
        assert_eq!(thermal.time_s(), 0.0);
        assert!(observer.first.sample().is_none());
    });
}

#[test]
fn complete_gap_owner_refuses_missing_competing_and_untransformed_paths() {
    with_case(false, false, |cx, _, heat, thermal, gas| {
        let layers = [gap()];
        let resolve = |_: &SourceEnclosedGasLayer, _: &DepositionSimulation<'_>, t: Temperature| {
            gas.evaluate(&gas_source::point(t.value(), 100000.0))
        };
        let initial =
            evaluate_source_enclosed_gas_layers(&thermal, &layers, 1, resolve, || false).unwrap();
        let prepared = heat
            .prepare_with_link_conductances(cx, &thermal, initial.conductances_w_k())
            .unwrap();
        let assembly = SourceEnclosedGasHeatAssembly {
            heat: &prepared,
            layers: &layers,
            maximum_layers: 1,
            resolve_gas: resolve,
        };
        assert!(prepared.resolve(cx, &thermal).is_err());
        assert!(
            SourceEnclosedGasHeatAssembly {
                layers: &[],
                ..assembly
            }
            .resolve_heat(cx, &thermal)
            .is_err()
        );
        assert!(
            SourceEnclosedGasHeatAssembly {
                maximum_layers: 0,
                ..assembly
            }
            .resolve_heat(cx, &thermal)
            .is_err()
        );
        assert!(
            assembly
                .resolve_heat_with_contact_metrics(
                    cx,
                    &thermal,
                    &[],
                    &[ThermalLinkMetric {
                        link: GAP_LINK,
                        area_m2: GAP_AREA,
                        distances_m: [GAP / 2.0; 2],
                    }]
                )
                .is_err()
        );
        let fourier = [ConductionLink {
            link: GAP_LINK,
            distance_a_m: GAP / 2.0,
            distance_b_m: GAP / 2.0,
        }];
        let bulk = MaterialNetwork {
            cells: heat.electrical.bulk.cells,
            resistors: &[],
            links: &fourier,
        };
        let electrical = ElectricalContactNetwork {
            bulk: &bulk,
            contacts: &[],
            electrical_bases: &[],
        };
        let competing = MaterialHeatNetwork {
            electrical: &electrical,
            thermal_contacts: &[],
            radiation: None,
        };
        assert!(
            competing
                .prepare_with_link_conductances(cx, &thermal, initial.conductances_w_k())
                .is_err()
        );
        let id = heat.electrical.bulk.cells[0].point.material.id();
        let interface = InterfaceSystemCard::assemble(
            SurfaceSpec {
                material: id.clone(),
                texture_frame: "lower".into(),
            },
            SurfaceSpec {
                material: id.clone(),
                texture_frame: "upper".into(),
            },
            SystemContext {
                medium: "competing gap contact".into(),
                third_body: None,
                environment: "manufactured".into(),
                history: "fixed".into(),
            },
            ClaimSet::new(),
            vec![],
        )
        .unwrap();
        let point = QueryPoint::new();
        let cells = heat.electrical.bulk.cells;
        let contact = [MaterialContactLink {
            link: GAP_LINK,
            point: ThermalInterfacePoint {
                interface: &interface,
                materials: [cells[0].point.material, cells[1].point.material],
                context: &point,
                temperature_axes: ["lower_T", "upper_T"],
                liquid_fraction_axes: [None; 2],
                mean_temperature_axis: None,
            },
            sides: [1, 0].map(|i| ContactCell {
                cell: cells[i].cell,
                material_id: cells[i].material_id,
                curve: cells[i].curve,
            }),
            selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        }];
        let competing = MaterialHeatNetwork {
            electrical: heat.electrical,
            thermal_contacts: &contact,
            radiation: None,
        };
        assert!(
            competing
                .prepare_with_link_conductances(cx, &thermal, initial.conductances_w_k())
                .is_err()
        );
    });
}

#[test]
fn native_gap_observer_rolls_back_late_source_mismatch_and_cancellation() {
    with_case(false, false, |cx, gate, heat, thermal, gas| {
        let layers = [gap()];
        let initial = evaluate_source_enclosed_gas_layers(
            &thermal,
            &layers,
            1,
            |_, _, t| gas.evaluate(&gas_source::point(t.value(), 100000.0)),
            || false,
        )
        .unwrap();
        let prepared = heat
            .prepare_with_link_conductances(cx, &thermal, initial.conductances_w_k())
            .unwrap();
        let saved = thermal.cells().clone();
        for cancel in [false, true] {
            let calls = std::cell::Cell::new(0);
            let assembly = SourceEnclosedGasHeatAssembly {
                heat: &prepared,
                layers: &layers,
                maximum_layers: 1,
                resolve_gas: |_: &SourceEnclosedGasLayer,
                              state: &DepositionSimulation<'_>,
                              t: Temperature| {
                    calls.set(calls.get() + 1);
                    if cancel && state.time_s() > 0.0 {
                        gate.request();
                    }
                    // A fresh original source at the wrong requested T must
                    // fail at the endpoint, after one private thermal step.
                    let offset = if !cancel && state.time_s() > 0.0 {
                        1.0
                    } else {
                        0.0
                    };
                    gas.evaluate(&gas_source::point(t.value() + offset, 100000.0))
                },
            };
            let observer = observe_heat_assembly(&assembly, cx, (), 1e-8).unwrap();
            assert!(
                thermal
                    .advance_observed(
                        &observer,
                        0.2,
                        StepConfig {
                            max_dt_s: 0.1,
                            energy_tolerance_j: 1e-8
                        },
                        || Ok(())
                    )
                    .is_err()
            );
            assert!(calls.get() >= 3);
            assert_eq!(thermal.cells(), &saved);
            assert_eq!(thermal.time_s(), 0.0);
            assert!(observer.first.sample().is_none());
        }
    });
}
