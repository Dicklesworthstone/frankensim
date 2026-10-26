//! Real shared-network transport with original source convection and hydraulics.
use super::*;
use fs_material::gas_field::steady::network::{GasOutlet, SteadyGasNetwork};
use fs_thermal::enthalpy::CellId;
use fs_thermal::enthalpy::advection::MassFlow;

fn recirculating(variable: bool) -> SteadyGasNetwork {
    let mut network = problem(variable).to_network().unwrap();
    network.flows = vec![
        MassFlow {
            id: 7,
            a: CellId(0),
            b: CellId(1),
            kg_s: 3.0 * FLOW,
        },
        MassFlow {
            id: 11,
            a: CellId(1),
            b: CellId(0),
            kg_s: 2.0 * FLOW,
        },
    ];
    network.inlets = vec![
        GasInjection {
            cell: 0,
            mass_flow_kg_s: FLOW,
            state: primitive(400.0, 1.0),
        },
        GasInjection {
            cell: 1,
            mass_flow_kg_s: FLOW,
            state: primitive(360.0, 0.0),
        },
    ];
    network.outlets = vec![GasOutlet {
        cell: 1,
        mass_flow_kg_s: 2.0 * FLOW,
    }];
    network
}

#[test]
fn recirculating_materials_close_species_wall_heat_and_original_enthalpy() {
    with_cx(|cx, _| {
        let p = recirculating(true);
        let saved = p.clone();
        let solved = solve_material_gas_network(cx, &p, &laws(), controls(), feedback()).unwrap();
        assert!(
            solved.iterations > 2,
            "actual temperature-dependent conductivity feeds back"
        );
        assert_eq!(solved.thermal.mass_flow_kg_s(), &[3.0 * FLOW, 4.0 * FLOW]);
        near(
            solved.thermal.species().cells()[&CellId(0)][0],
            2.0 / 3.0,
            2e-13,
        );
        near(solved.thermal.species().cells()[&CellId(1)][0], 0.5, 2e-13);
        let mut h = [0.0; 2];
        for (index, cell) in solved.cells.iter().enumerate() {
            let gas = cell.gas.state();
            h[index] = cell.gas.specific_enthalpy_j_kg();
            assert_eq!(&cell.gas, &solved.thermal.states()[index]);
            assert_eq!(cell.gas.source_packs(), p.source.source_packs());
            assert_eq!(cell.hydraulic.transport, cell.transport);
            // Independent fully-developed circular-duct relation, using the
            // conductivity of THIS final source-resolved mixture and state.
            let expected_g = 3.66 * gas.thermal_conductivity / DIAMETER
                * p.cells[index].geometry.surface_area_m2;
            near(cell.conductance_w_k, expected_g, 1e-12);
            near(
                cell.wall_heat_w,
                expected_g * (p.cells[index].wall_temperature_k - gas.temperature),
                1e-11,
            );
            let mass = solved.thermal.mass_flow_kg_s()[index];
            near(cell.velocity_m_s, mass / (gas.density * area()), 1e-12);
            near(
                cell.hydraulic.pressure_drop_pa,
                32.0 * gas.dynamic_viscosity * LENGTH * cell.velocity_m_s / DIAMETER.powi(2),
                1e-11,
            );
        }
        let inlet_h = |index: usize| {
            let supply = &p.inlets[index];
            p.source
                .with_composition(supply.state.composition.clone())
                .unwrap()
                .evaluate(&source::point(
                    supply.state.temperature_k,
                    p.static_pressure_pa,
                ))
                .unwrap()
                .specific_enthalpy_j_kg()
        };
        let q = [solved.cells[0].wall_heat_w, solved.cells[1].wall_heat_w];
        // Two independently reconstructed cell balances; their circulation
        // terms cancel from the boundary account exactly once.
        near(
            FLOW * inlet_h(0) + 2.0 * FLOW * h[1] - 3.0 * FLOW * h[0] + q[0],
            0.0,
            2e-7,
        );
        near(
            FLOW * inlet_h(1) + 3.0 * FLOW * h[0] - 4.0 * FLOW * h[1] + q[1],
            0.0,
            2e-7,
        );
        near(
            FLOW * (inlet_h(0) + inlet_h(1)) + q[0] + q[1] - 2.0 * FLOW * h[1],
            0.0,
            4e-7,
        );
        assert!(solved.maximum_power_residual_bound_w <= controls().steady_power_tolerance_w);
        assert!(solved.substeps <= feedback().maximum_substeps);
        assert_eq!(p, saved);
    });
}

#[test]
fn data_only_species_and_pressure_changes_reach_recirculating_heat_and_flow() {
    with_cx(|cx, _| {
        let mut p = recirculating(true);
        let mixed = solve_material_gas_network(cx, &p, &laws(), controls(), feedback()).unwrap();
        p.inlets[1].state.composition = primitive(360.0, 1.0).composition;
        let replaced = solve_material_gas_network(cx, &p, &laws(), controls(), feedback()).unwrap();
        for index in 0..2 {
            assert_eq!(
                replaced.thermal.species().cells()[&CellId(index as u64)],
                vec![1.0, 0.0]
            );
            let a = &mixed.cells[index];
            let b = &replaced.cells[index];
            assert!((a.gas.state().temperature - b.gas.state().temperature).abs() > 1.0);
            assert_ne!(a.conductance_w_k, b.conductance_w_k);
            assert_ne!(a.gas.state().sound_speed, b.gas.state().sound_speed);
            assert_ne!(a.hydraulic.pressure_drop_pa, b.hydraulic.pressure_drop_pa);
        }
        p.static_pressure_pa *= 2.0;
        let high = solve_material_gas_network(cx, &p, &laws(), controls(), feedback()).unwrap();
        for index in 0..2 {
            near(
                high.cells[index].gas.state().temperature,
                replaced.cells[index].gas.state().temperature,
                2e-5,
            );
            near(
                high.cells[index].gas.state().density / replaced.cells[index].gas.state().density,
                2.0,
                2e-6,
            );
            near(
                high.cells[index].hydraulic.pressure_drop_pa,
                0.5 * replaced.cells[index].hydraulic.pressure_drop_pa,
                1e-5,
            );
        }
    });
}

#[test]
fn branched_sections_use_local_states_and_axial_api_keeps_the_same_feedback() {
    with_cx(|cx, _| {
        let axial = problem(true);
        let legacy =
            solve_material_gas_transport(cx, &axial, &laws(), controls(), feedback()).unwrap();
        let mut p = axial.to_network().unwrap();
        let general = solve_material_gas_network(cx, &p, &laws(), controls(), feedback()).unwrap();
        for (a, b) in legacy.cells.iter().zip(&general.cells) {
            near(a.gas.state().temperature, b.gas.state().temperature, 1e-9);
            near(a.conductance_w_k, b.conductance_w_k, 1e-12);
            near(
                a.hydraulic.pressure_drop_pa,
                b.hydraulic.pressure_drop_pa,
                1e-9,
            );
        }
        p.cells.push(p.cells[1].clone());
        p.cells[2].wall_temperature_k = 340.0;
        p.flows = vec![
            MassFlow {
                id: 2,
                a: CellId(0),
                b: CellId(1),
                kg_s: FLOW,
            },
            MassFlow {
                id: 3,
                a: CellId(0),
                b: CellId(2),
                kg_s: FLOW,
            },
        ];
        p.inlets[0].mass_flow_kg_s = 2.0 * FLOW;
        p.outlets = vec![
            GasOutlet {
                cell: 1,
                mass_flow_kg_s: FLOW,
            },
            GasOutlet {
                cell: 2,
                mass_flow_kg_s: FLOW,
            },
        ];
        let branched =
            solve_material_gas_network(cx, &p, &[laws()[0]; 3], controls(), feedback()).unwrap();
        assert_eq!(branched.thermal.mass_flow_kg_s(), &[2.0 * FLOW, FLOW, FLOW]);
        assert!(
            branched.cells[2].gas.state().temperature > branched.cells[1].gas.state().temperature
        );
        assert_ne!(
            branched.cells[2].conductance_w_k,
            branched.cells[1].conductance_w_k
        );
        near(
            branched.total_pressure_drop_pa,
            branched
                .cells
                .iter()
                .map(|c| c.hydraulic.pressure_drop_pa)
                .sum(),
            1e-12,
        );
        assert!(branched.maximum_power_residual_bound_w <= controls().steady_power_tolerance_w);
    });
}

#[test]
fn network_owner_geometry_source_work_and_cancellation_fail_without_mutation() {
    with_cx(|cx, gate| {
        let p = recirculating(true);
        let saved = p.clone();
        let mut invalid = p.clone();
        invalid.cells[0].wall_conductance_w_k = 1.0;
        assert!(matches!(
            solve_material_gas_network(cx, &invalid, &laws(), controls(), feedback()),
            Err(MaterialSteadyGasError::Invalid(_))
        ));
        invalid = p.clone();
        invalid.cells[0].flow_area_m2 *= 2.0;
        assert!(matches!(
            solve_material_gas_network(cx, &invalid, &laws(), controls(), feedback()),
            Err(MaterialSteadyGasError::Invalid(_))
        ));
        assert!(matches!(
            solve_material_gas_network(
                cx,
                &p,
                &laws(),
                controls(),
                GasTransportFeedback {
                    maximum_iterations: 1,
                    ..feedback()
                }
            ),
            Err(MaterialSteadyGasError::NotConverged { .. })
        ));
        invalid = p.clone();
        invalid.inlets[0].state.temperature_k = 1001.0;
        assert!(matches!(
            solve_material_gas_network(cx, &invalid, &laws(), controls(), feedback()),
            Err(MaterialSteadyGasError::Thermal(_))
        ));
        assert_eq!(p, saved);
        gate.request();
        assert!(matches!(
            solve_material_gas_network(cx, &p, &laws(), controls(), feedback()),
            Err(MaterialSteadyGasError::Cancelled)
        ));
        assert_eq!(p, saved);
    });
}
