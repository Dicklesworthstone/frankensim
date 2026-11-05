//! Mixed source controls close the original material/contact/enthalpy endpoint.
use super::*;
use std::collections::BTreeSet;

#[test]
fn mixed_voltage_current_and_floating_electrodes_refresh_material_and_contact_heat() {
    // Two unit-volume receivers: C_A=1, C_B=2 J/K, G=1 W/K, dt=.5 s.
    // The manufactured endpoint has R_A=5 ohm (bulk plus internal contact),
    // R_B=3 ohm, I_A=2 A and I_B=1 or 2 A. The middle ideal electrode joins
    // otherwise separate volumes; I_middle=2-I_B, including a floating case.
    // Independent finite-cell balances give these exact temperature rises:
    // 1.5*dTa-.5*dTb=10, -.5*dTa+2.5*dTb=1.5*I_B².
    for right_current in [1.0_f64, 2.0] {
        let theta = [
            (25.0 + 0.75 * right_current.powi(2)) / 3.5,
            (5.0 + 2.25 * right_current.powi(2)) / 3.5,
        ];
        let bulk_r = 1.0 + 0.002 * theta[0];
        let contact_r = 5.0 - bulk_r;
        let (a, ca) = material(1.0, 1.0, 0.002);
        let (b, cb) = material(3.0, 2.0, 0.0);
        let cards = [a, b];
        let curves = [ca, cb];
        let context = point(300.0);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = bindings(&cards, &curves, &context, &key);
        let thermal = initial(&curves, 999.0, None);
        // Both sides of this internal seam belong to receiver A. Its actual
        // endpoint temperature drives BOTH source axes, not a guessed average.
        let r0 = contact_r / (1.0 + 0.01 * theta[0]);
        let interface = interface(&[cards[0].clone(), cards[0].clone()], r0, 1000.0);
        let ip = interface_point(300.0, 1e5);
        let surface = [0, 1].map(|contact| ConductionContactThermal {
            contact,
            point: ThermalInterfacePoint {
                interface: &interface,
                materials: [&cards[0], &cards[0]],
                context: &ip,
                temperature_axes: ["Ta", "Tb"],
                liquid_fraction_axes: [None, None],
                mean_temperature_axis: None,
            },
            heat_fraction_to_first: 0.25,
        });
        let conductors = cards.each_ref().map(|card| {
            resolve_conductor_state_with_basis(
                card,
                &context,
                MaterialPropertySelection::SingleClaimOnly,
                ElectricalPropertyBasis::Resistivity,
            )
            .unwrap()
        });
        let source = resolve_electrical_interface_state(
            &interface,
            &ip,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let (mut positions, mut tets) = slab(&[0.0, 0.5]);
        for xs in [&[0.5, 1.0][..], &[1.0, 1.5, 2.0][..]] {
            let (part, elements) = slab(xs);
            let offset = positions.len() as u32;
            positions.extend(part);
            tets.extend(elements.into_iter().map(|t| t.map(|v| v + offset)));
        }
        let problem = SpatialConduction::new_with_contacts(
            positions,
            tets,
            (0..24).map(|i| conductors[i / 12].clone()).collect(),
            contacts(source),
            limits(),
            &mut || false,
        )
        .unwrap();
        let regions = (0..24)
            .map(|tetrahedron| SpatialConductionRegion {
                tetrahedron,
                cell: cells[tetrahedron / 12].cell,
            })
            .collect::<Vec<_>>();
        let network = SpatialConductionThermalNetwork {
            cells: &cells,
            links: &LINKS,
            regions: &regions,
            volume_tolerance_m3: 1e-12,
        };
        let faces = [0.0, 1.0, 2.0].map(|x| {
            problem
                .exterior_faces()
                .filter(|&face| {
                    problem.mesh().faces[face].map(|v| problem.positions_m()[v as usize][0])
                        == [x; 3]
                })
                .collect::<Vec<_>>()
        });
        let terminals = faces.each_ref().map(Vec::as_slice);
        let response = SpatialConductanceControls {
            maximum_terminals: 3,
            maximum_field_solves: 2,
            maximum_response_values: 1000,
            maximum_dense_work: 1000,
            reference_voltage_v: 0.5,
            coefficient_tolerance_s: 1e-10,
        };
        let mixed = SpatialMixedConductanceControls {
            maximum_current_terminals: 2,
            maximum_response_values: 100,
            maximum_dense_work: 100,
        };
        let power = 20.0 + 3.0 * right_current.powi(2);
        let end_voltage = 10.0 + 3.0 * right_current;
        assert!((1.0 + r0 - 5.0).abs() > 0.1);
        for sign in [-1.0, 1.0] {
            for bias in [0.0, 1e10] {
                let currents = [sign * (2.0 - right_current), -2.0 * sign];
                let voltages = [bias, bias + sign * end_voltage];
                // Also prescribe original terminal zero's current, in reversed
                // caller order. The sole voltage electrode is then terminal 2.
                for current_terminals in [&[1][..], &[1, 0][..]] {
                    let drive = SpatialConductionDrive::MixedTerminals {
                        terminals: &terminals,
                        current_terminals,
                        voltages_v: if current_terminals.len() == 1 {
                            &voltages
                        } else {
                            &voltages[1..]
                        },
                        currents_a: &currents[..current_terminals.len()],
                        response,
                        mixed,
                    };
                    let mut polls = 0;
                    let next = network
                        .step_contact_candidate(
                            &thermal,
                            0.5,
                            &problem,
                            drive,
                            &surface,
                            limits(),
                            field_controls(),
                            config(),
                            &[0.0; 2],
                            || {
                                polls += 1;
                                false
                            },
                        )
                        .unwrap();
                    for (cell, rise) in [CellId(30), CellId(10)].into_iter().zip(theta) {
                        near(next.thermal().temperature_k(cell).unwrap(), 300.0 + rise);
                    }
                    let report = next.report();
                    let field = &report.electrical.field;
                    for ((terminal, voltage), current) in terminals
                        .iter()
                        .zip([0.0, 10.0 * sign, end_voltage * sign])
                        .zip([-2.0 * sign, currents[0], right_current * sign])
                    {
                        let vertices = terminal
                            .iter()
                            .flat_map(|&face| problem.mesh().faces[face])
                            .collect::<BTreeSet<_>>();
                        near(
                            vertices
                                .iter()
                                .map(|&v| field.nodal_boundary_current_a[v as usize].unwrap())
                                .sum(),
                            current,
                        );
                        for vertex in vertices {
                            near(field.voltage_v[vertex as usize], voltage);
                        }
                    }
                    near(field.current_density_a_m2[0][0], -2.0 * sign);
                    near(field.current_density_a_m2[12][0], -right_current * sign);
                    near(field.cell_power_w[..12].iter().sum(), 4.0 * bulk_r);
                    near(
                        field.cell_power_w[12..].iter().sum(),
                        3.0 * right_current.powi(2),
                    );
                    near(field.power_w, power);
                    near(report.electrical.heat_j, 0.5 * power);
                    near(report.electrical.supplied_work_j, 0.5 * power);
                    near(report.deposited_joule_j, 0.5 * power);
                    near(
                        next.thermal().energy_j() - thermal.energy_j(),
                        theta[0] + 2.0 * theta[1],
                    );
                    near(report.total_energy_residual_j, 0.0);
                    assert!(report.circuit_trials >= 2);
                    assert_eq!(report.circuit_trials % 2, 0);
                    near(
                        next.problem().conductors()[0].material().query().axes()["T"],
                        300.0 + theta[0],
                    );
                    for (heat, contact) in report
                        .electrical
                        .contacts
                        .iter()
                        .zip(next.problem().contacts())
                    {
                        assert_eq!(heat.cells, [CellId(30); 2]);
                        near(heat.heat_fraction_to_first, 0.25);
                        near(heat.response.current_a, -sign);
                        near(heat.response.power_w, 2.0 * contact_r);
                        near(contact.source.resistance_ohm_m2(), contact_r);
                        for axis in ["Ta", "Tb"] {
                            near(
                                contact.source.resolved().state_coordinate(axis).unwrap(),
                                300.0 + theta[0],
                            );
                        }
                    }
                    if right_current == 1.0
                        && sign == 1.0
                        && bias == 0.0
                        && current_terminals.len() == 1
                    {
                        let mut bounded = config();
                        bounded.maximum_circuit_trials = 1;
                        assert!(matches!(
                            network.step_contact_candidate(
                                &thermal,
                                0.5,
                                &problem,
                                drive,
                                &surface,
                                limits(),
                                field_controls(),
                                bounded,
                                &[0.0; 2],
                                || false,
                            ),
                            Err(CircuitEnthalpyError::CircuitBudget)
                        ));
                        let mut cancelled_polls = 0;
                        assert!(matches!(
                            network.step_contact_candidate(
                                &thermal,
                                0.5,
                                &problem,
                                drive,
                                &surface,
                                limits(),
                                field_controls(),
                                config(),
                                &[0.0; 2],
                                || {
                                    cancelled_polls += 1;
                                    cancelled_polls == polls
                                },
                            ),
                            Err(CircuitEnthalpyError::Cancelled)
                        ));
                    }
                }
            }
        }
        // Every trial and late cancellation leaves the accepted source/state.
        assert_eq!(thermal.time_s(), 0.0);
        near(thermal.temperature_k(CellId(30)).unwrap(), 300.0);
        near(
            problem.conductors()[0].material().query().axes()["T"],
            300.0,
        );
        near(problem.contacts()[0].source.resistance_ohm_m2(), r0);
    }
}
