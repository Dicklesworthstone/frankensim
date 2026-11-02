//! A total-current electrode lets the local current redistribute with R(T).
use super::*;
use std::collections::BTreeSet;

#[test]
fn multiterminal_currents_close_source_feedback_and_charge_each_native_basis_solve() {
    let (a, ca) = material(1.0, 1.0, 0.002);
    let (b, cb) = material(3.0, 2.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let context = point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = bindings(&cards, &curves, &context, &key);
    let thermal = initial(&curves, 999.0, None);
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
    let (second_positions, second_tets) = slab(&[1.0, 1.5, 2.0]);
    let offset = positions.len() as u32;
    positions.extend(second_positions);
    tets.extend(second_tets.into_iter().map(|t| t.map(|v| v + offset)));
    let problem = SpatialConduction::new(
        positions,
        tets,
        (0..24).map(|i| conductors[i / 12].clone()).collect(),
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
    // The middle equipotential electrode joins the two separate volume meshes.
    // A prescribed middle current changes the two bars' currents independently.
    let faces = [0.0, 1.0, 2.0].map(|x| {
        problem
            .exterior_faces()
            .filter(|&f| {
                problem.mesh().faces[f].map(|v| problem.positions_m()[v as usize][0]) == [x; 3]
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
    for right_current in [1.0, 2.0] {
        // Independent local BE equations at dt=.5:
        // (1.5-.004)*dTa-.5*dTb=2, -.5*dTa+2.5*dTb=1.5*I_b².
        let aa = 1.5 - 0.004;
        let bb = 2.5;
        let determinant = aa * bb - 0.25;
        let rhs_b = 1.5 * right_current * right_current;
        let theta = [
            (2.0 * bb + 0.5 * rhs_b) / determinant,
            (aa * rhs_b + 1.0) / determinant,
        ];
        let rho_a = 1.0 + 0.002 * theta[0];
        let endpoint_voltage = [0.0, 2.0 * rho_a, 2.0 * rho_a + 3.0 * right_current];
        let heat = 0.5 * (4.0 * rho_a + 3.0 * right_current * right_current);
        for sign in [1.0, -1.0] {
            let currents = [
                -2.0 * sign,
                (2.0 - right_current) * sign,
                right_current * sign,
            ];
            let drive = SpatialConductionDrive::TerminalCurrents {
                terminals: &terminals,
                currents_a: &currents,
                response,
            };
            let next = network
                .step_candidate(
                    &thermal,
                    0.5,
                    &problem,
                    drive,
                    limits(),
                    field_controls(),
                    config(),
                    &[0.0; 2],
                    || false,
                )
                .unwrap();
            near(
                next.thermal().temperature_k(CellId(30)).unwrap(),
                300.0 + theta[0],
            );
            near(
                next.thermal().temperature_k(CellId(10)).unwrap(),
                300.0 + theta[1],
            );
            let report = next.report();
            for (i, terminal) in terminals.iter().enumerate() {
                near(
                    terminal_current(&problem, &report.electrical.field, terminal),
                    currents[i],
                );
                for &face in *terminal {
                    for vertex in problem.mesh().faces[face] {
                        near(
                            report.electrical.field.voltage_v[vertex as usize],
                            sign * endpoint_voltage[i],
                        );
                    }
                }
            }
            for (i, density) in report
                .electrical
                .field
                .current_density_a_m2
                .iter()
                .enumerate()
            {
                near(density[0], -sign * if i < 12 { 2.0 } else { right_current });
                near(density[1], 0.0);
                near(density[2], 0.0);
                near(
                    report.electrical.field.cell_power_w[i],
                    if i < 12 {
                        4.0 * rho_a / 12.0
                    } else {
                        3.0 * right_current * right_current / 12.0
                    },
                );
            }
            near(
                next.problem().conductors()[0].material().query().axes()["T"],
                300.0 + theta[0],
            );
            near(report.electrical.heat_j, heat);
            near(report.electrical.supplied_work_j, heat);
            near(
                next.thermal().energy_j() - thermal.energy_j(),
                theta[0] + 2.0 * theta[1],
            );
            near(report.total_energy_residual_j, 0.0);
            assert!(report.circuit_trials >= 2);
            assert_eq!(report.circuit_trials % 2, 0);
            let mut bounded = config();
            bounded.maximum_circuit_trials = 1;
            assert!(matches!(
                network.step_candidate(
                    &thermal,
                    0.5,
                    &problem,
                    drive,
                    limits(),
                    field_controls(),
                    bounded,
                    &[0.0; 2],
                    || false,
                ),
                Err(CircuitEnthalpyError::CircuitBudget)
            ));
        }
    }
    assert_eq!(thermal.time_s(), 0.0);
    for conductor in problem.conductors() {
        near(conductor.material().query().axes()["T"], 300.0);
    }
}

fn terminal_current(
    problem: &SpatialConduction,
    field: &SpatialConductionField,
    faces: &[usize],
) -> f64 {
    faces
        .iter()
        .flat_map(|&face| problem.mesh().faces[face])
        .collect::<BTreeSet<_>>()
        .into_iter()
        .map(|vertex| field.nodal_boundary_current_a[vertex as usize].unwrap())
        .sum()
}

#[test]
fn total_current_redistributes_with_temperature_and_closes_both_parallel_body_balances() {
    // Two electrically isolated parallel bars share the same two electrodes
    // and one authored thermal link G=1 W/K. C=[1,2] J/K, dt=.5 s.
    // Prescribe an independently constructed exact endpoint dT=[10,3] K:
    // P=[(10+.5*7)/.5, (2*3-.5*7)/.5]=[27,5] W.
    // V=3 V therefore requires R=[1/3,9/5] ohm and I=32/3 A.
    // Only the first resistance changes with T; its initial value differs.
    let alpha = 0.004;
    let initial_rho = (1.0 / 3.0) / (1.0 + alpha * 10.0);
    let (a, ca) = material(initial_rho, 1.0, alpha);
    let (b, cb) = material(9.0 / 5.0, 2.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let context = point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = bindings(&cards, &curves, &context, &key);
    let thermal = initial(&curves, 999.0, None);
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
    let (second_positions, second_tets) = slab(&[0.0, 0.5, 1.0]);
    let offset = positions.len() as u32;
    positions.extend(
        second_positions
            .into_iter()
            .map(|[x, y, z]| [x, y + 1.0, z]),
    );
    tets.extend(second_tets.into_iter().map(|tet| tet.map(|v| v + offset)));
    let problem = SpatialConduction::new(
        positions,
        tets,
        (0..24).map(|i| conductors[i / 12].clone()).collect(),
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
    let [positive, negative] = electrodes(&problem, 1.0);
    let current = 32.0 / 3.0;
    let cold_voltage = current / (1.0 / initial_rho + 5.0 / 9.0);
    assert!((cold_voltage - 3.0).abs() > 0.01);
    let before_energy = thermal.energy_j();
    for sign in [1.0, -1.0] {
        for reference_voltage_v in [0.25, 4.0] {
            let next = network
                .step_candidate(
                    &thermal,
                    0.5,
                    &problem,
                    SpatialConductionDrive::TotalCurrent {
                        positive_faces: &positive,
                        negative_faces: &negative,
                        current_a: sign * current,
                        reference_voltage_v,
                    },
                    limits(),
                    field_controls(),
                    config(),
                    &[0.0; 2],
                    || false,
                )
                .unwrap();
            near(next.thermal().temperature_k(CellId(30)).unwrap(), 310.0);
            near(next.thermal().temperature_k(CellId(10)).unwrap(), 303.0);
            let report = next.report();
            let field = &report.electrical.field;
            near(terminal_current(&problem, field, &positive), sign * current);
            near(
                terminal_current(&problem, field, &negative),
                -sign * current,
            );
            for (position, &voltage) in problem.positions_m().iter().zip(&field.voltage_v) {
                near(voltage, sign * 3.0 * position[0]);
            }
            for (i, density) in field.current_density_a_m2.iter().enumerate() {
                let power = [27.0, 5.0][i / 12];
                near(density[0], -sign * power / 3.0);
                near(density[1], 0.0);
                near(density[2], 0.0);
                near(field.cell_power_w[i], power / 12.0);
            }
            near(
                next.problem().conductors()[0].material().query().axes()["T"],
                310.0,
            );
            near(
                next.problem().conductors()[12].material().query().axes()["T"],
                303.0,
            );
            near(field.power_w, 32.0);
            near(field.supplied_power_w, 32.0);
            near(report.electrical.heat_j, 16.0);
            near(report.electrical.supplied_work_j, 16.0);
            near(report.deposited_joule_j, 16.0);
            near(next.thermal().energy_j() - before_energy, 16.0);
            near(report.total_energy_residual_j, 0.0);
        }
    }

    let zero = network
        .step_candidate(
            &thermal,
            0.5,
            &problem,
            SpatialConductionDrive::total_current(&positive, &negative, 0.0),
            limits(),
            field_controls(),
            config(),
            &[0.0; 2],
            || false,
        )
        .unwrap();
    near(zero.thermal().energy_j(), before_energy);
    near(zero.report().electrical.heat_j, 0.0);
    near(zero.report().electrical.supplied_work_j, 0.0);
    assert!(
        zero.report()
            .electrical
            .field
            .voltage_v
            .iter()
            .all(|&v| v == 0.0)
    );

    for drive in [
        SpatialConductionDrive::total_current(&positive, &positive, current),
        SpatialConductionDrive::total_current(&positive, &[], 0.0),
        SpatialConductionDrive::total_current(&positive, &negative, f64::NAN),
        SpatialConductionDrive::TotalCurrent {
            positive_faces: &positive,
            negative_faces: &negative,
            current_a: current,
            reference_voltage_v: 0.0,
        },
    ] {
        assert!(
            network
                .step_candidate(
                    &thermal,
                    0.5,
                    &problem,
                    drive,
                    limits(),
                    field_controls(),
                    config(),
                    &[0.0; 2],
                    || false,
                )
                .is_err()
        );
    }
    near(thermal.energy_j(), before_energy);
    assert_eq!(thermal.time_s(), 0.0);
    for source in problem.conductors() {
        near(source.material().query().axes()["T"], 300.0);
    }
}
