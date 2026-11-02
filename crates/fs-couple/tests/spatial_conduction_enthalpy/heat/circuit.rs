//! G1/G3/G4: full three-electrode C/L storage and complete original heat transport.
use super::*;
use fs_couple::circuit_enthalpy::{
    ResistorCell,
    conduction::{
        heat::{ConductionHeatError, step_spatial_conduction_circuit_candidate},
        spatial::circuit::SpatialConductionCircuit,
    },
    material::CellResistor,
};
use fs_phs::circuit::{Branch, CircuitGraph};

fn response() -> SpatialConductanceControls {
    SpatialConductanceControls {
        maximum_terminals: 3,
        maximum_field_solves: 2,
        maximum_response_values: 1000,
        maximum_dense_work: 1000,
        reference_voltage_v: 0.5,
        coefficient_tolerance_s: 1e-10,
    }
}

fn initial_cells<'a>(cells: &[MaterialCell<'a>]) -> DepositionSimulation<'a> {
    let mut events = cells
        .iter()
        .enumerate()
        .map(|(i, cell)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: cell.cell,
                segment_id: i as u64,
                body_id: i as u64,
                material_id: cell.material_id,
                mass_kg: 1.0,
                specific_h_j_kg: 100.0 * (i + 1) as f64,
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 1.0, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: 6.0,
                },
                boundary: Boundary {
                    reservoir: (i == 2).then_some(Reservoir {
                        area_m2: 1.0,
                        conductance_w_k: 0.4,
                        temperature_k: 310.0,
                    }),
                    ..common::boundary()
                },
            }),
        })
        .collect::<Vec<_>>();
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 7,
            a: CellId(30),
            b: CellId(10),
            area_m2: 1.0,
            conductance_w_k: 999.0,
        }),
    });
    let laws = cells
        .iter()
        .map(|cell| (cell.material_id, cell.curve as &dyn ThermalLaw))
        .collect::<Vec<_>>();
    DepositionSimulation::new(
        events,
        &laws,
        DepositionBudget {
            max_cells: 8,
            max_links: 8,
            max_events: 16,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.1,
            energy_tolerance_j: 1e-7,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

fn with_circuit_heat<R>(
    emissivity: f64,
    pressure: f64,
    upper: f64,
    sign: f64,
    run: impl FnOnce(
        &Cx<'_>,
        &CancelGate,
        &MaterialHeatNetwork<'_>,
        DepositionSimulation<'_>,
        SpatialConductionCircuit,
        &[SpatialConductionRegion],
    ) -> R,
) -> R {
    with_heat(
        emissivity,
        pressure,
        upper,
        |cx, gate, original, _, original_field, _| {
            let (r, cr) = material(2.0, 3.0, 0.003);
            let context = point(300.0);
            let mut cells = original
                .electrical
                .bulk
                .cells
                .iter()
                .map(|cell| MaterialCell {
                    cell: cell.cell,
                    material_id: cell.material_id,
                    curve: cell.curve,
                    point: cell.point,
                    conductivity: cell.conductivity,
                })
                .collect::<Vec<_>>();
            cells.push(MaterialCell {
                cell: CellId(20),
                material_id: MaterialId(2),
                curve: cr.curve(),
                point: ThermalMaterialPoint {
                    material: &r,
                    context: &context,
                    temperature_axis: "T",
                    liquid_fraction_axis: None,
                },
                conductivity: None,
            });
            let resistors = [CellResistor {
                destination: ResistorCell {
                    branch: 2,
                    cell: CellId(20),
                },
                selection: MaterialPropertySelection::SingleClaimOnly,
                length_m: 1.0,
                area_m2: 1.0,
            }];
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
            let heat = MaterialHeatNetwork {
                electrical: &electrical,
                thermal_contacts: original.thermal_contacts,
                radiation: original.radiation,
            };
            let thermal = initial_cells(&cells);
            let sources = bulk.resolve(&thermal, || false).unwrap().resistors;
            // Two original spatial volumes meet only through the middle electrode;
            // the thermal link and enclosure retain their independent physical laws.
            let (mut positions, mut tets) = slab(&[0.0, 0.5, 1.0]);
            let (other, other_tets) = slab(&[1.0, 1.5, 2.0]);
            let offset = positions.len() as u32;
            positions.extend(other);
            tets.extend(other_tets.into_iter().map(|tet| tet.map(|v| v + offset)));
            let field = SpatialConduction::new(
                positions,
                tets,
                (0..24)
                    .map(|i| original_field.conductors()[6 * (i / 12)].clone())
                    .collect(),
                limits(),
                &mut || false,
            )
            .unwrap();
            let faces = [0.0, 1.0, 2.0].map(|x| {
                field
                    .exterior_faces()
                    .filter(|&face| {
                        field.mesh().faces[face].map(|v| field.positions_m()[v as usize][0])
                            == [x; 3]
                    })
                    .collect::<Vec<_>>()
            });
            let regions = (0..24)
                .map(|tetrahedron| SpatialConductionRegion {
                    tetrahedron,
                    cell: cells[tetrahedron / 12].cell,
                })
                .collect::<Vec<_>>();
            let circuit = SpatialConductionCircuit::new_with_material_resistors(
                field,
                CircuitGraph {
                    node_count: 4,
                    branches: vec![
                        (3, 0, Branch::Capacitor { farads: 1.0 }),
                        (2, 0, Branch::Inductor { henries: 2.0 }),
                        (3, 1, Branch::Resistor { ohms: 999.0 }),
                        (3, 0, Branch::CurrentSource { port: 0 }),
                    ],
                    transformers: vec![],
                },
                faces.to_vec(),
                vec![0, 1, 2],
                sources,
                0.0,
                vec![sign * 0.5],
                vec![sign * 4.0],
                &[sign],
                response(),
                field_controls(),
                || false,
            )
            .unwrap();
            run(cx, gate, &heat, thermal, circuit, &regions)
        },
    )
}

struct Expected {
    theta: [f64; 3],
    resistance: [f64; 3],
    joule: [f64; 3],
    voltage: f64,
    inductor_current: f64,
    branch_current: f64,
    capacitor_voltage: f64,
    bath: f64,
    residual: f64,
}

fn bisect(mut lo: f64, mut hi: f64, residual: impl Fn(f64) -> f64) -> f64 {
    assert!(residual(lo) <= 0.0 && residual(hi) >= 0.0);
    for _ in 0..90 {
        let mid = 0.5 * (lo + hi);
        if residual(mid) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

// Independent scalar physical elimination, never invoking a native field,
// circuit, thermal coefficient resolver or the production nonlinear solver.
// Midpoint C/L fixes I_L=(phi+h V)/(2+3h). Two-body total heat eliminates T_B;
// the A balance uses Fourier/contact transfer and exact two-grey-surface
// radiosity. The separate resistor's BE balance includes its 310 K bath.
fn expected(q: f64, phi: f64, old: [f64; 3], dt: f64, emissivity: f64, pressure: f64) -> Expected {
    let contact_r = 2.0 - 1.8 * (pressure - 100000.0) / 100000.0;
    let endpoint = |v: f64| {
        let h = 0.5 * dt;
        let il = (phi + h * v) / (2.0 + 3.0 * h);
        let qb = dt * 3.0 * il * il;
        let balance = |ta: f64| {
            let ra = 1.0 + 0.002 * ta;
            let qa = dt * v * v / ra;
            let tb = (old[0] + 2.0 * old[1] + qa + qb - ta) / 2.0;
            let ea = emissivity * (1.0 + 0.0001 * ta);
            let eb = 0.6 * (1.0 + 0.0001 * tb);
            let transfer = (ta - tb) / (1.0 + contact_r)
                + 5.670_374_419e-8 * ((300.0 + ta).powi(4) - (300.0 + tb).powi(4))
                    / (1.0 / ea + 1.0 / eb - 1.0);
            (ra, qa, tb, ta - old[0] - qa + dt * transfer)
        };
        let ta = bisect(0.0, old[0] + 2.0 * old[1] + dt * v * v + qb, |ta| {
            balance(ta).3
        });
        let (ra, qa, tb, _) = balance(ta);
        let ir = v / ra + il;
        let tr = (3.0 * old[2] + dt * (2.0 * ir * ir + 0.4 * 10.0))
            / (3.0 + dt * 0.4 - dt * 0.006 * ir * ir);
        let rr = 2.0 * (1.0 + 0.003 * tr);
        let vc = v + rr * ir;
        Expected {
            theta: [ta, tb, tr],
            resistance: [ra, 3.0, rr],
            joule: [qa, qb, dt * rr * ir * ir],
            voltage: v,
            inductor_current: il,
            branch_current: ir,
            capacitor_voltage: vc,
            bath: dt * 0.4 * (10.0 - tr),
            residual: vc + h * ir - q - h,
        }
    };
    endpoint(bisect(0.0, q + 0.5 * dt, |v| endpoint(v).residual))
}

#[test]
fn complete_heat_drives_three_electrode_storage_and_material_resistor_over_two_steps() {
    let mut substitutions = Vec::new();
    for (emissivity, pressure) in [(0.2, 100000.0), (0.8, 100000.0), (0.2, 200000.0)] {
        for sign in [1.0, -1.0] {
            let result = with_circuit_heat(
                emissivity,
                pressure,
                2000.0,
                sign,
                |cx, _, heat, mut thermal, mut circuit, regions| {
                    let prepared = heat.prepare(cx, &thermal).unwrap();
                    let (mut q, mut phi, mut old) = (4.0, 0.5, [0.0; 3]);
                    let mut final_current = 0.0;
                    for interval in 1..=2 {
                        let dt = 0.5;
                        let oracle = expected(q, phi, old, dt, emissivity, pressure);
                        let target = interval as f64 * dt;
                        let next = if interval == 1 {
                            step_spatial_conduction_circuit_candidate(
                                &prepared,
                                cx,
                                &thermal,
                                target,
                                &circuit,
                                &[sign],
                                &[],
                                regions,
                                1e-12,
                                response(),
                                limits(),
                                field_controls(),
                                config(),
                                &[0.0; 3],
                            )
                            .unwrap()
                        } else {
                            prepared
                                .step_spatial_conduction_circuit_candidate(
                                    cx,
                                    &thermal,
                                    target,
                                    &circuit,
                                    &[sign],
                                    &[],
                                    regions,
                                    1e-12,
                                    response(),
                                    limits(),
                                    field_controls(),
                                    config(),
                                    &[0.0; 3],
                                )
                                .unwrap()
                        };
                        let report = next.report();
                        let step = &report.electrical;
                        for (id, theta) in [CellId(30), CellId(10), CellId(20)]
                            .into_iter()
                            .zip(oracle.theta)
                        {
                            near(next.thermal().temperature_k(id).unwrap(), 300.0 + theta);
                        }
                        near(step.terminal_voltage_v[1], sign * oracle.voltage);
                        near(
                            step.terminal_voltage_v[2],
                            sign * (oracle.voltage - 3.0 * oracle.inductor_current),
                        );
                        near(
                            step.terminal_current_a[0],
                            -sign * oracle.voltage / oracle.resistance[0],
                        );
                        near(step.terminal_current_a[1], sign * oracle.branch_current);
                        near(step.terminal_current_a[2], -sign * oracle.inductor_current);
                        for i in 0..2 {
                            near(
                                step.spatial.field.cell_power_w[i * 12..(i + 1) * 12]
                                    .iter()
                                    .sum::<f64>()
                                    * dt,
                                oracle.joule[i],
                            );
                        }
                        assert_eq!(step.resistors.len(), 1);
                        near(step.resistors[0].heat_j, oracle.joule[2]);
                        near(step.resistors[0].current_a, sign * oracle.branch_current);
                        near(
                            step.resistors[0].material.conductor.resistance_ohm(),
                            oracle.resistance[2],
                        );
                        near(
                            step.resistors[0]
                                .material
                                .conductor
                                .material()
                                .state_coordinate("T")
                                .unwrap(),
                            300.0 + oracle.theta[2],
                        );
                        near(
                            report.source.thermal_contacts[0]
                                .interface
                                .value_m2_k_per_w(),
                            2.0 - 1.8 * (pressure - 100000.0) / 100000.0,
                        );
                        near(step.spatial.heat_j, oracle.joule[0] + oracle.joule[1]);
                        near(step.heat_j, oracle.joule.iter().sum());
                        near(step.native.dissipated, step.heat_j);
                        near(report.deposited_joule_j, step.heat_j);
                        near(
                            report.non_electrical_heat_j
                                - dt * report.source.radiation_net_power_w(),
                            oracle.bath,
                        );
                        let next_q = 2.0 * oracle.capacitor_voltage - q;
                        let next_phi = 4.0 * oracle.inductor_current - phi;
                        near(next.circuit().charges_c()[0], sign * next_q);
                        near(next.circuit().fluxes_wb()[0], sign * next_phi);
                        let stored = 0.5 * (next_q * next_q - q * q)
                            + 0.25 * (next_phi * next_phi - phi * phi);
                        near(step.stored_change_j, stored);
                        near(step.supplied_work_j, dt * oracle.capacitor_voltage);
                        near(
                            thermal.energy_j() + step.heat_j + oracle.bath,
                            next.thermal().energy_j(),
                        );
                        near(
                            stored + next.thermal().energy_j()
                                - thermal.energy_j()
                                - step.supplied_work_j
                                - oracle.bath,
                            0.0,
                        );
                        near(
                            report.total_energy_residual_j
                                + dt * report.source.radiation_net_power_w(),
                            0.0,
                        );
                        assert_eq!(report.circuit_trials % 3, 0);
                        near(circuit.charges_c()[0], sign * q);
                        near(circuit.fluxes_wb()[0], sign * phi);
                        assert_eq!(thermal.time_s(), target - dt);
                        q = next_q;
                        phi = next_phi;
                        old = oracle.theta;
                        final_current = sign * step.terminal_current_a[1];
                        (thermal, circuit, _) = next.into_parts();
                    }
                    [old[0], old[1], old[2], final_current]
                },
            );
            if sign > 0.0 {
                substitutions.push(result);
            } else {
                for i in 0..4 {
                    near(result[i], substitutions.last().unwrap()[i]);
                }
            }
        }
    }
    for changed in &substitutions[1..] {
        for i in [0, 1, 3] {
            assert!(
                (changed[i] - substitutions[0][i]).abs() > 1e-8,
                "the original heat-source substitution must change body temperatures and circuit current"
            );
        }
    }
}

#[test]
fn unsupported_complete_heat_and_late_refusal_publish_no_circuit_or_caloric_history() {
    for upper in [2000.0, 300.0001] {
        with_circuit_heat(
            0.2,
            100000.0,
            upper,
            1.0,
            |cx, gate, heat, thermal, circuit, regions| {
                let prepared = heat.prepare(cx, &thermal).unwrap();
                let energy = thermal.energy_j();
                let sources = circuit.problem().conductors().to_vec();
                let run = |controls| {
                    step_spatial_conduction_circuit_candidate(
                        &prepared,
                        cx,
                        &thermal,
                        0.5,
                        &circuit,
                        &[1.0],
                        &[],
                        regions,
                        1e-12,
                        response(),
                        limits(),
                        field_controls(),
                        controls,
                        &[0.0; 3],
                    )
                };
                if upper < 301.0 {
                    assert!(
                        run(config()).is_err(),
                        "unsupported hot emissivity must not retain its accepted 300 K slice"
                    );
                } else {
                    let mut short = config();
                    short.maximum_circuit_trials = 2;
                    assert!(matches!(
                        run(short),
                        Err(CircuitEnthalpyError::CircuitBudget)
                    ));
                    let wrong = [CellResistor {
                        destination: ResistorCell {
                            branch: 2,
                            cell: CellId(20),
                        },
                        selection: MaterialPropertySelection::SingleClaimOnly,
                        length_m: 0.5,
                        area_m2: 2.0,
                    }];
                    let bulk = MaterialNetwork {
                        cells: heat.electrical.bulk.cells,
                        resistors: &wrong,
                        links: heat.electrical.bulk.links,
                    };
                    let electrical = ElectricalContactNetwork {
                        bulk: &bulk,
                        contacts: &[],
                        electrical_bases: &[],
                    };
                    let changed = MaterialHeatNetwork {
                        electrical: &electrical,
                        thermal_contacts: heat.thermal_contacts,
                        radiation: heat.radiation,
                    };
                    let changed = changed.prepare(cx, &thermal).unwrap();
                    assert!(
                        matches!(
                            step_spatial_conduction_circuit_candidate(
                                &changed,
                                cx,
                                &thermal,
                                0.5,
                                &circuit,
                                &[1.0],
                                &[],
                                regions,
                                1e-12,
                                response(),
                                limits(),
                                field_controls(),
                                config(),
                                &[0.0; 3],
                            ),
                            Err(CircuitEnthalpyError::Source(
                                ConductionHeatError::Conduction(_)
                            ))
                        ),
                        "same-volume branch geometry cannot replace the accepted original conductor"
                    );
                }
                gate.request();
                assert!(matches!(
                    run(config()),
                    Err(CircuitEnthalpyError::Cancelled)
                ));
                assert_eq!(thermal.time_s(), 0.0);
                near(thermal.energy_j(), energy);
                near(circuit.charges_c()[0], 4.0);
                near(circuit.fluxes_wb()[0], 0.5);
                assert_eq!(circuit.problem().conductors(), sources);
                near(
                    circuit.material_resistors()[0]
                        .conductor
                        .material()
                        .state_coordinate("T")
                        .unwrap(),
                    300.0,
                );
            },
        );
    }
}
