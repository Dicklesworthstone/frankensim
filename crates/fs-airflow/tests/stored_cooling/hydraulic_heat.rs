//! Pump losses heat the same native spatial solid and finite reservoir histories.
use super::*;

fn heated_controls() -> StoredCoolingControls {
    let mut selected = controls();
    selected
        .pump
        .network
        .work
        .maximum_neglected_hydraulic_power_w = 0.0;
    selected.pump.operating_point.pressure_bracket_pa = (10e6, 30e6);
    selected.reservoir.maximum_neglected_hydraulic_energy_j = 0.0;
    selected
}

fn high_head_pump() -> PumpHeadCurve {
    PumpHeadCurve::try_new(
        vec![
            PumpHeadPoint {
                flow: VolumetricFlowRate::new(0.0),
                head: Length::new(2000.0),
            },
            PumpHeadPoint {
                flow: VolumetricFlowRate::new(0.001),
                head: Length::new(2000.0),
            },
        ],
        10.0,
        (1e-6, 3e-6),
        SourceProvenance::new(
            "synthetic high-head incompressible energy regression, not equipment data",
            "hydraulic-storage-test",
        ),
    )
    .unwrap()
}

fn high_resistance() -> [LiquidResistanceSection<'static>; 1] {
    [LiquidResistanceSection {
        region: "body",
        resistance: MaterialFlowResistance::try_new(
            1e16,
            0.0,
            Area::new(core::f64::consts::PI * 0.01 * 0.01 / 4.0),
            Length::new(0.01),
            (0.0, 2300.0),
        )
        .unwrap(),
    }]
}

#[test]
fn hydraulic_heat_advances_native_solid_and_reservoir_once_across_steps() {
    fixture(1.0, |cx, _, ordinary, _, state, mesh, material| {
        let resistance = high_resistance();
        let original = ordinary.branches()[0];
        let heated = ParallelLiquidNetwork::try_new(
            cx,
            &[ParallelLiquidBranch {
                resistances: &resistance,
                pressure_budget: PressureDrivenBudget {
                    maximum_neglected_hydraulic_power_w: 0.0,
                    ..original.pressure_budget
                },
                ..original
            }],
        )
        .unwrap()
        .with_hydraulic_dissipation();
        let pump = high_head_pump();
        let source = ScalarField::uniform("no body heater", 0.0).unwrap();
        let initial_solid = state.solid_specific_enthalpy_j_kg().to_vec();
        let initial_tank = state.reservoir().state().specific_enthalpy_j_kg();
        let mut walls = vec![300.0];
        let mut previous_tank_temperature = 300.0;
        for end in [1.0, 2.0] {
            let result = heated
                .advance_stored_solid_with_pump(
                    cx,
                    state,
                    end,
                    0.0,
                    &walls,
                    &pump,
                    StoredSolidProblem {
                        mesh,
                        material,
                        source: &source,
                        element_materials: None,
                        interfaces: None,
                    },
                    heated_controls(),
                    |_, rows| boundary(mesh, rows),
                )
                .unwrap();
            let liquid = &result.liquid;
            let flow = &liquid.circuit.state.fluid;
            let solid = &liquid.circuit.state.solid.state;
            // dp = K*mu*mdot/rho selects 0.001 kg/s, hence dp*mdot/rho = 20 W.
            assert!((flow.total_mass_flow_kg_per_s - 0.001).abs() < 1e-10);
            assert!((liquid.deposited_hydraulic_energy_j - 20.0).abs() < 1e-6);
            assert!((liquid.circuit.supply.hydraulic_power_w - 20.0).abs() < 1e-6);
            assert_eq!(liquid.neglected_hydraulic_energy_j, 0.0);
            assert_eq!(result.solid_external_heat_j, 0.0);
            assert_eq!(liquid.external_heat_j, 0.0);
            assert!(solid.stored_energy_change_j > 0.0);
            assert!(flow.total_heat_rate_w < 0.0); // Heated fluid gives heat TO the body.
            assert!(liquid.returned_heat_j > 0.0);
            assert!((liquid.returned_heat_j - flow.total_enthalpy_gain_w).abs() < 1e-8);
            assert!(
                (solid.stored_energy_change_j + liquid.stored_enthalpy_change_j - 20.0).abs()
                    < 2e-6
            );
            assert!(
                result.interface_energy_residual_j.abs()
                    <= heated_controls().joint_energy_tolerance_j
            );
            assert!(
                result.combined_energy_residual_j.abs()
                    <= heated_controls().joint_energy_tolerance_j
            );
            assert!(state.reservoir().state().temperature_k() > previous_tank_temperature);
            assert_eq!(state.reservoir().state(), liquid.final_state);
            assert_eq!(state.reservoir().time_s(), end);
            assert_eq!(
                state.solid_specific_enthalpy_j_kg(),
                solid.specific_enthalpy_j_kg.as_slice()
            );
            assert!(liquid.total_solid_solves <= heated_controls().pump.max_total_solid_solves);
            assert!(
                liquid.total_flow_evaluations
                    <= heated_controls().pump.network.work.max_flow_evaluations
            );
            for branch in &flow.branches {
                assert_eq!(branch.solution.thermal.regions[0].inlet, liquid.final_state);
                assert!(std::ptr::eq(
                    branch.solution.thermal.source,
                    state.reservoir().source()
                ));
            }
            previous_tank_temperature = state.reservoir().state().temperature_k();
            walls = liquid
                .circuit
                .state
                .solid
                .regions
                .iter()
                .map(|r| r.mean_wall_temperature_k)
                .collect();
        }
        let solid_j: f64 = state
            .storage()
            .reference_nodal_masses_kg()
            .iter()
            .zip(state.solid_specific_enthalpy_j_kg())
            .zip(initial_solid)
            .map(|((&mass, &h), old)| mass * (h - old))
            .sum();
        let tank_j = state.reservoir().mass_kg()
            * (state.reservoir().state().specific_enthalpy_j_kg() - initial_tank);
        assert!((solid_j + tank_j - 40.0).abs() < 4e-6);
    });
}

#[test]
fn selected_hydraulic_model_keeps_joint_storage_transaction_and_work_gates() {
    fixture(1.0, |cx, cancel, ordinary, _, state, mesh, material| {
        let resistance = high_resistance();
        let original = ordinary.branches()[0];
        let heated = ParallelLiquidNetwork::try_new(
            cx,
            &[ParallelLiquidBranch {
                resistances: &resistance,
                pressure_budget: PressureDrivenBudget {
                    maximum_neglected_hydraulic_power_w: 0.0,
                    ..original.pressure_budget
                },
                ..original
            }],
        )
        .unwrap()
        .with_hydraulic_dissipation();
        let pump = high_head_pump();
        let source = ScalarField::uniform("no body heater", 0.0).unwrap();
        let problem = StoredSolidProblem {
            mesh,
            material,
            source: &source,
            element_materials: None,
            interfaces: None,
        };
        let old_tank = state.reservoir().clone();
        let old_solid = state.solid_specific_enthalpy_j_kg().to_vec();
        let called = std::cell::Cell::new(false);
        let mut invalid = heated_controls();
        invalid.reservoir.maximum_neglected_hydraulic_energy_j = 1.0;
        assert!(
            heated
                .advance_stored_solid_with_pump(
                    cx,
                    state,
                    1.0,
                    0.0,
                    &[300.0],
                    &pump,
                    problem,
                    invalid,
                    |_, rows| {
                        called.set(true);
                        boundary(mesh, rows)
                    }
                )
                .is_err()
        );
        assert!(!called.get());
        assert_eq!(state.reservoir(), &old_tank);
        assert_eq!(state.solid_specific_enthalpy_j_kg(), old_solid.as_slice());
        let mut limited = heated_controls();
        limited.pump.max_total_solid_solves = 1;
        assert!(
            heated
                .advance_stored_solid_with_pump(
                    cx,
                    state,
                    1.0,
                    0.0,
                    &[300.0],
                    &pump,
                    problem,
                    limited,
                    |_, rows| boundary(mesh, rows)
                )
                .is_err()
        );
        assert_eq!(state.reservoir(), &old_tank);
        assert_eq!(state.solid_specific_enthalpy_j_kg(), old_solid.as_slice());
        assert!(
            heated
                .advance_stored_solid_with_pump(
                    cx,
                    state,
                    1.0,
                    0.0,
                    &[300.0],
                    &pump,
                    problem,
                    heated_controls(),
                    |_, rows| {
                        cancel.request();
                        boundary(mesh, rows)
                    }
                )
                .is_err()
        );
        assert_eq!(state.reservoir(), &old_tank);
        assert_eq!(state.solid_specific_enthalpy_j_kg(), old_solid.as_slice());
    });
}
