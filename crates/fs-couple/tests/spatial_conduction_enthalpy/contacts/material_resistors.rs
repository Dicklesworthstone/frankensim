//! Material branch heat, full field heat and contact heat share one endpoint.
use super::*;
use fs_couple::circuit_enthalpy::{
    ResistorCell, conduction::spatial::circuit::SpatialConductionCircuit, material::CellResistor,
};
use fs_phs::circuit::{Branch, CircuitGraph};

fn initial_three(curves: &[ResolvedHeatCapacityCurve; 3]) -> DepositionSimulation<'_> {
    let ids = [CellId(30), CellId(10), CellId(20)];
    let schedule = ids
        .into_iter()
        .enumerate()
        .map(|(i, id)| {
            let knots = curves[i].curve().knots();
            let cp = (knots.last().unwrap().specific_enthalpy_j_kg
                - knots[0].specific_enthalpy_j_kg)
                / 1800.0;
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id,
                    segment_id: i as u64,
                    body_id: i as u64,
                    material_id: MaterialId(i as u64),
                    mass_kg: 1.0,
                    specific_h_j_kg: 100.0 * cp,
                    geometry: CellGeometry {
                        start_m: [i as f64, 0.0, 0.0],
                        end_m: [i as f64 + 1.0, 0.0, 0.0],
                        volume_m3: 1.0,
                        surface_area_m2: 6.0,
                    },
                    boundary: common::boundary(),
                }),
            }
        })
        .collect();
    let laws = curves
        .iter()
        .enumerate()
        .map(|(i, curve)| (MaterialId(i as u64), curve.curve() as &dyn ThermalLaw))
        .collect::<Vec<_>>();
    DepositionSimulation::new(
        schedule,
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

fn three_cells<'a>(
    cards: [&'a MaterialCard; 3],
    curves: &'a [ResolvedHeatCapacityCurve; 3],
    contexts: [&'a QueryPoint; 3],
) -> Vec<MaterialCell<'a>> {
    [CellId(30), CellId(10), CellId(20)]
        .into_iter()
        .enumerate()
        .map(|(i, cell)| MaterialCell {
            cell,
            material_id: MaterialId(i as u64),
            curve: curves[i].curve(),
            point: ThermalMaterialPoint {
                material: cards[i],
                context: contexts[i],
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: None,
        })
        .collect()
}

fn branch_binding(branch: usize, cell: CellId) -> CellResistor {
    CellResistor {
        destination: ResistorCell { branch, cell },
        selection: MaterialPropertySelection::SingleClaimOnly,
        length_m: 1.0,
        area_m2: 1.0,
    }
}

fn response(terminals: usize) -> SpatialConductanceControls {
    SpatialConductanceControls {
        maximum_terminals: terminals,
        maximum_field_solves: terminals - 1,
        maximum_response_values: 1000,
        maximum_dense_work: 1000,
        reference_voltage_v: 0.5,
        coefficient_tolerance_s: 1e-10,
    }
}

#[test]
fn material_resistor_and_three_electrode_field_refresh_under_one_signed_rlc_history() {
    for sign in [1.0, -1.0] {
        let (a, ca) = material(1.0, 1.0, 0.002);
        let (b, cb) = material(3.0, 2.0, 0.0);
        let (r, cr) = material(2.0, 3.0, 0.003);
        let cards = [a, b];
        let curves = [ca, cb, cr];
        let context = point(300.0);
        let cells = three_cells([&cards[0], &cards[1], &r], &curves, [&context; 3]);
        let resistor_cells = [branch_binding(2, CellId(20))];
        let material_network = MaterialNetwork {
            cells: &cells,
            resistors: &resistor_cells,
            links: &[],
        };
        let mut thermal = initial_three(&curves);
        let sources = material_network
            .resolve(&thermal, || false)
            .unwrap()
            .resistors;
        let conductors = cards.each_ref().map(|card| {
            resolve_conductor_state_with_basis(
                card,
                &context,
                MaterialPropertySelection::SingleClaimOnly,
                ElectricalPropertyBasis::Resistivity,
            )
            .unwrap()
        });
        let (mut positions, mut tets) = slab(&[0.0, 0.5, 1.0]);
        let (other, other_tets) = slab(&[1.0, 1.5, 2.0]);
        let offset = positions.len() as u32;
        positions.extend(other);
        tets.extend(other_tets.into_iter().map(|tet| tet.map(|v| v + offset)));
        let field = SpatialConduction::new(
            positions,
            tets,
            (0..24).map(|i| conductors[i / 12].clone()).collect(),
            limits(),
            &mut || false,
        )
        .unwrap();
        let faces = [0.0, 1.0, 2.0].map(|x| {
            field
                .exterior_faces()
                .filter(|&f| {
                    field.mesh().faces[f].map(|v| field.positions_m()[v as usize][0]) == [x; 3]
                })
                .collect::<Vec<_>>()
        });
        let regions = (0..24)
            .map(|tetrahedron| SpatialConductionRegion {
                tetrahedron,
                cell: cells[tetrahedron / 12].cell,
            })
            .collect::<Vec<_>>();
        let layout = SpatialConductionThermalNetwork {
            cells: &cells,
            links: &[],
            regions: &regions,
            volume_tolerance_m3: 1e-12,
        };
        // The source and C feed a separate material resistor before the middle
        // field electrode. The right electrode retains the original L history.
        let graph = CircuitGraph {
            node_count: 4,
            branches: vec![
                (3, 0, Branch::Capacitor { farads: 1.0 }),
                (2, 0, Branch::Inductor { henries: 2.0 }),
                (3, 1, Branch::Resistor { ohms: 999.0 }),
                (3, 0, Branch::CurrentSource { port: 0 }),
            ],
            transformers: vec![],
        };
        let ports = [sign];
        let (mut charge, mut flux) = (4.0, 0.5);
        let mut circuit = SpatialConductionCircuit::new_with_material_resistors(
            field,
            graph,
            faces.to_vec(),
            vec![0, 1, 2],
            sources,
            0.0,
            vec![sign * flux],
            vec![sign * charge],
            &ports,
            response(3),
            field_controls(),
            || false,
        )
        .unwrap();
        let mut theta = [0.0; 3];
        let dt = 0.5;
        for interval in 1..=2 {
            // Independent scalar elimination: BE heat changes each original R;
            // midpoint C/L determine current. These three volumes have no heat link.
            let endpoint = |v: f64| {
                let h = 0.5 * dt;
                let base = 1.0 + 0.002 * theta[0];
                let rhs = dt * v * v;
                let da = 2.0 * rhs / (base + (base * base + 0.008 * rhs).sqrt());
                let ta = theta[0] + da;
                let ra = 1.0 + 0.002 * ta;
                let il = (flux + h * v) / (2.0 + 3.0 * h);
                let ir = v / ra + il;
                let gain = dt * 2.0 * ir * ir / 3.0;
                let tr = (theta[2] + gain) / (1.0 - 0.003 * gain);
                let rr = 2.0 * (1.0 + 0.003 * tr);
                let vc = v + rr * ir;
                (ta, tr, ra, il, ir, rr, vc, vc + h * ir - charge - h)
            };
            let (mut lo, mut hi) = (0.0, charge + 0.5 * dt);
            assert!(endpoint(lo).7 < 0.0 && endpoint(hi).7 > 0.0);
            for _ in 0..90 {
                let mid = 0.5 * (lo + hi);
                if endpoint(mid).7 < 0.0 {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            let v = 0.5 * (lo + hi);
            let (ta, tr, ra, il, ir, rr, vc, _) = endpoint(v);
            let heat = [dt * v * v / ra, dt * 3.0 * il * il, dt * rr * ir * ir];
            let tb = theta[1] + heat[1] / 2.0;
            let final_charge = 2.0 * vc - charge;
            let final_flux = 4.0 * il - flux;
            let mut polls = 0;
            let next = layout
                .step_material_circuit_candidate(
                    &thermal,
                    interval as f64 * dt,
                    &circuit,
                    &ports,
                    &resistor_cells,
                    &[],
                    response(3),
                    limits(),
                    field_controls(),
                    config(),
                    &[0.0; 3],
                    || {
                        polls += 1;
                        false
                    },
                )
                .unwrap();
            for (cell, t) in [(CellId(30), ta), (CellId(10), tb), (CellId(20), tr)] {
                near(next.thermal().temperature_k(cell).unwrap(), 300.0 + t);
            }
            let step = &next.report().electrical;
            near(next.circuit().charges_c()[0], sign * final_charge);
            near(next.circuit().fluxes_wb()[0], sign * final_flux);
            near(step.terminal_voltage_v[1], sign * v);
            near(step.terminal_voltage_v[2], sign * (v - 3.0 * il));
            near(step.terminal_current_a[0], -sign * v / ra);
            near(step.terminal_current_a[1], sign * ir);
            near(step.terminal_current_a[2], -sign * il);
            near(
                step.spatial.field.cell_power_w[..12].iter().sum::<f64>() * dt,
                heat[0],
            );
            near(
                step.spatial.field.cell_power_w[12..].iter().sum::<f64>() * dt,
                heat[1],
            );
            assert_eq!(step.resistors.len(), 1);
            near(step.resistors[0].current_a, sign * ir);
            near(step.resistors[0].heat_j, heat[2]);
            near(step.resistors[0].material.conductor.resistance_ohm(), rr);
            near(
                step.resistors[0]
                    .material
                    .conductor
                    .material()
                    .state_coordinate("T")
                    .unwrap(),
                300.0 + tr,
            );
            near(step.spatial.heat_j, heat[0] + heat[1]);
            near(step.heat_j, heat.iter().sum());
            near(step.native.dissipated, step.heat_j);
            near(step.supplied_work_j, dt * vc);
            let stored = 0.5 * (final_charge * final_charge - charge * charge)
                + 0.25 * (final_flux * final_flux - flux * flux);
            near(step.stored_change_j, stored);
            near(stored + step.heat_j, step.supplied_work_j);
            near(next.thermal().energy_j() - thermal.energy_j(), step.heat_j);
            near(next.report().total_energy_residual_j, 0.0);
            assert_eq!(next.report().circuit_trials % 3, 0);
            if interval == 1 {
                // Missing sources, branch changes and overlapping whole-cell
                // receivers refuse before either accepted history is published.
                for bindings in [
                    vec![],
                    vec![branch_binding(2, CellId(30))],
                    vec![branch_binding(1, CellId(20))],
                ] {
                    assert!(
                        layout
                            .step_material_circuit_candidate(
                                &thermal,
                                dt,
                                &circuit,
                                &ports,
                                &bindings,
                                &[],
                                response(3),
                                limits(),
                                field_controls(),
                                config(),
                                &[0.0; 3],
                                || false
                            )
                            .is_err()
                    );
                }
                let mut short = config();
                short.maximum_circuit_trials = 2;
                assert!(matches!(
                    layout.step_material_circuit_candidate(
                        &thermal,
                        dt,
                        &circuit,
                        &ports,
                        &resistor_cells,
                        &[],
                        response(3),
                        limits(),
                        field_controls(),
                        short,
                        &[0.0; 3],
                        || false
                    ),
                    Err(CircuitEnthalpyError::CircuitBudget)
                ));
                let mut calls = 0;
                assert!(matches!(
                    layout.step_material_circuit_candidate(
                        &thermal,
                        dt,
                        &circuit,
                        &ports,
                        &resistor_cells,
                        &[],
                        response(3),
                        limits(),
                        field_controls(),
                        config(),
                        &[0.0; 3],
                        || {
                            calls += 1;
                            calls == polls
                        }
                    ),
                    Err(CircuitEnthalpyError::Cancelled)
                ));
                near(circuit.charges_c()[0], sign * charge);
                near(circuit.fluxes_wb()[0], sign * flux);
                near(thermal.energy_j(), 600.0);
                assert_eq!(thermal.time_s(), 0.0);
            }
            theta = [ta, tb, tr];
            charge = final_charge;
            flux = final_flux;
            (thermal, circuit, _) = next.into_parts();
        }
    }
}

#[test]
fn separate_conductivity_branch_bulk_and_two_sided_contact_heat_close_capacitor_work() {
    let (r0, dt) = (0.75, 0.5);
    for upper in [1000.0, 300.0001] {
        let (a, ca) = material(1.0, 1.0, 0.002);
        let (b, cb) = material(3.0, 2.0, 0.0);
        // Deliberately use the other original electrical basis and retain a
        // nonthermal source coordinate; no reciprocal source card is invented.
        let (r, cr) = common::material_with_electrical_basis(
            ElectricalPropertyBasis::Conductivity,
            0.5,
            3.0,
            0.004,
            1.0,
        );
        let cards = [a, b];
        let curves = [ca, cb, cr];
        let context = point(300.0);
        let resistor_context = common::point(300.0);
        let cells = three_cells(
            [&cards[0], &cards[1], &r],
            &curves,
            [&context, &context, &resistor_context],
        );
        let resistors = [branch_binding(1, CellId(20))];
        let material_network = MaterialNetwork {
            cells: &cells,
            resistors: &resistors,
            links: &[],
        };
        let thermal = initial_three(&curves);
        let sources = material_network
            .resolve_with_bases(
                &thermal,
                &[(1, ElectricalPropertyBasis::Conductivity)],
                || false,
            )
            .unwrap()
            .resistors;
        let interface = interface(&cards, r0, upper);
        let ip = interface_point(300.0, 1e5);
        let surface = surface_bindings(&cards, &interface, &ip);
        let field = specimen(&cards, &interface, &ip);
        let [positive, negative] = electrodes(&field, 2.0);
        let regions = regions();
        let layout = SpatialConductionThermalNetwork {
            cells: &cells,
            links: &[],
            regions: &regions,
            volume_tolerance_m3: 1e-12,
        };
        let circuit = SpatialConductionCircuit::new_with_material_resistors(
            field,
            CircuitGraph {
                node_count: 3,
                branches: vec![
                    (2, 0, Branch::Capacitor { farads: 1.0 }),
                    (2, 1, Branch::Resistor { ohms: 999.0 }),
                ],
                transformers: vec![],
            },
            vec![negative, positive],
            vec![0, 1],
            sources,
            0.0,
            vec![],
            vec![4.0],
            &[],
            response(2),
            field_controls(),
            || false,
        )
        .unwrap();
        let next = layout.step_material_circuit_candidate(
            &thermal,
            dt,
            &circuit,
            &[],
            &resistors,
            &surface,
            response(2),
            limits(),
            field_controls(),
            config(),
            &[0.0; 3],
            || false,
        );
        if upper < 301.0 {
            assert!(next.is_err());
            assert_eq!(thermal.time_s(), 0.0);
            near(thermal.energy_j(), 600.0);
            near(circuit.charges_c()[0], 4.0);
            near(
                circuit.material_resistors()[0]
                    .conductor
                    .material()
                    .state_coordinate("T")
                    .unwrap(),
                300.0,
            );
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
        let unsupported_context = point(300.0).with("lot", 8.0).unwrap();
        let unsupported_cells = three_cells(
            [&cards[0], &cards[1], &r],
            &curves,
            [&context, &context, &unsupported_context],
        );
        let unsupported = SpatialConductionThermalNetwork {
            cells: &unsupported_cells,
            links: &[],
            regions: &regions,
            volume_tolerance_m3: 1e-12,
        };
        assert!(
            unsupported
                .step_material_circuit_candidate(
                    &thermal,
                    dt,
                    &circuit,
                    &[],
                    &resistors,
                    &surface,
                    response(2),
                    limits(),
                    field_controls(),
                    config(),
                    &[0.0; 3],
                    || false,
                )
                .is_err()
        );
        near(circuit.charges_c()[0], 4.0);
        near(thermal.energy_j(), 600.0);
        let next = next.unwrap();
        let endpoint = |current: f64| {
            let s = dt * current * current;
            let a = 1.0 - s * (0.002 + 0.25 * r0 * 0.003);
            let b = -s * 0.25 * r0 * 0.007;
            let c = -s * 0.75 * r0 * 0.003;
            let d = 2.0 - s * 0.75 * r0 * 0.007;
            let ra = s * (1.0 + 0.25 * r0);
            let rb = s * (3.0 + 0.75 * r0);
            let ta = (ra * d - b * rb) / (a * d - b * c);
            let tb = (a * rb - c * ra) / (a * d - b * c);
            let rhs = 2.0 * s / 3.0;
            let tr = 2.0 * rhs / (1.0 + (1.0 + 0.016 * rhs).sqrt());
            let bulk_a = 1.0 + 0.002 * ta;
            let contact = r0 * (1.0 + 0.003 * ta + 0.007 * tb);
            let branch = 2.0 / (1.0 + 0.004 * tr);
            (
                ta,
                tb,
                tr,
                bulk_a,
                contact,
                branch,
                current * (bulk_a + 3.0 + contact + branch + 0.5 * dt) - 4.0,
            )
        };
        let (mut lo, mut hi) = (0.0, 1.0);
        assert!(endpoint(lo).6 < 0.0 && endpoint(hi).6 > 0.0);
        for _ in 0..90 {
            let mid = 0.5 * (lo + hi);
            if endpoint(mid).6 < 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let current = 0.5 * (lo + hi);
        let (ta, tb, tr, ra, rc, rr, _) = endpoint(current);
        let s = dt * current * current;
        let step = &next.report().electrical;
        for (cell, theta) in [(CellId(30), ta), (CellId(10), tb), (CellId(20), tr)] {
            near(next.thermal().temperature_k(cell).unwrap(), 300.0 + theta);
        }
        near(
            step.spatial.field.cell_power_w[..6].iter().sum::<f64>() * dt,
            s * ra,
        );
        near(
            step.spatial.field.cell_power_w[6..].iter().sum::<f64>() * dt,
            s * 3.0,
        );
        near(
            step.spatial
                .field
                .contacts
                .iter()
                .map(|c| c.power_w)
                .sum::<f64>()
                * dt,
            s * rc,
        );
        near(step.resistors[0].heat_j, s * rr);
        near(step.spatial.heat_j, s * (ra + 3.0 + rc));
        near(step.heat_j, s * (ra + 3.0 + rc + rr));
        near(next.thermal().energy_j() - thermal.energy_j(), step.heat_j);
        let final_charge = 4.0 - dt * current;
        near(next.circuit().charges_c()[0], final_charge);
        near(
            step.stored_change_j,
            0.5 * (final_charge * final_charge - 16.0),
        );
        near(step.stored_change_j + step.heat_j, 0.0);
        near(step.supplied_work_j, 0.0);
        near(step.native.dissipated, step.heat_j);
        near(next.report().total_energy_residual_j, 0.0);
        let conductor = &next.circuit().material_resistors()[0].conductor;
        assert_eq!(
            conductor.source_basis(),
            ElectricalPropertyBasis::Conductivity
        );
        near(conductor.material().state_coordinate("lot").unwrap(), 7.0);
        near(
            conductor.material().state_coordinate("T").unwrap(),
            300.0 + tr,
        );
        near(conductor.resistance_ohm(), rr);
        for contact in next.circuit().problem().contacts() {
            near(
                contact.source.resolved().state_coordinate("Ta").unwrap(),
                300.0 + ta,
            );
            near(
                contact.source.resolved().state_coordinate("Tb").unwrap(),
                300.0 + tb,
            );
        }
    }
}
