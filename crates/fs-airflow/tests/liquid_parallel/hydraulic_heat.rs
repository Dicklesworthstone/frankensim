//! Pressure loss participates in branch transport, mixing and material feedback.
use super::*;

fn heated_work() -> ParallelPressureBudget {
    ParallelPressureBudget {
        maximum_neglected_hydraulic_power_w: 0.0,
        ..budget()
    }
}

fn with_heated_network<R>(
    variable_mu: bool,
    variable_cp: bool,
    run: impl FnOnce(&Cx<'_>, &ParallelLiquidNetwork<'_, '_, '_>) -> R,
) -> R {
    with_source(variable_mu, variable_cp, |cx, law, source| {
        let a = path(cx, law, source, "a-wall");
        let b = path(cx, law, source, "b-wall");
        let ar = loss("a-wall", 1e6);
        let br = loss("b-wall", 2e6);
        let branch = PressureDrivenBudget {
            mass_flow_bracket_kg_per_s: (0.0002, 0.003),
            maximum_neglected_hydraulic_power_w: 0.0,
            ..branch_budget()
        };
        let network = ParallelLiquidNetwork::try_new(
            cx,
            &[
                ParallelLiquidBranch {
                    name: "a",
                    path: &a,
                    resistances: &ar,
                    pressure_budget: branch,
                },
                ParallelLiquidBranch {
                    name: "b",
                    path: &b,
                    resistances: &br,
                    pressure_budget: branch,
                },
            ],
        )
        .unwrap()
        .with_hydraulic_dissipation();
        run(cx, &network)
    })
}

#[test]
fn heated_parallel_outlets_match_independent_wall_and_pressure_energy() {
    with_heated_network(false, false, |cx, network| {
        assert!(network.deposits_hydraulic_dissipation());
        let walls = [300.0, 330.0];
        let solved = network
            .solve_pressure_driven(cx, &walls, Pressure::new(20e6), local(), heated_work())
            .unwrap();
        let mut expected_wall = 0.0;
        let mut expected_gain = 0.0;
        for ((branch, mass), wall) in solved.branches.iter().zip([0.001, 0.0005]).zip(walls) {
            // Constant Cp solution of mdot*Cp*dT/dx = UA*(Tw-T) + dp*Q.
            let capacity = mass * 2000.0;
            let hydraulic = 20e6 * mass / 1000.0;
            let conductance = 183.0 * 0.003;
            let outlet = 300.0
                + (wall - 300.0 + hydraulic / conductance) * (-(-conductance / capacity).exp_m1());
            let gain = capacity * (outlet - 300.0);
            near(branch.solution.mass_flow_kg_per_s, mass);
            near(
                branch.solution.hydraulics.total_dissipated_power_w,
                hydraulic,
            );
            near(
                branch.solution.thermal.regions[0].outlet.temperature_k(),
                outlet,
            );
            near(branch.solution.thermal.total_heat_rate_w, gain - hydraulic);
            expected_wall += gain - hydraulic;
            expected_gain += gain;
        }
        assert!(solved.branches[0].solution.thermal.total_heat_rate_w < 0.0);
        assert!(
            solved.branches[0].solution.thermal.regions[0]
                .outlet
                .temperature_k()
                > 300.0
        );
        near(solved.deposited_hydraulic_power_w, 30.0);
        near(solved.total_dissipated_power_w, 30.0);
        near(solved.total_heat_rate_w, expected_wall);
        near(solved.total_enthalpy_gain_w, expected_gain);
        near(
            solved.mixed_outlet.temperature_k(),
            300.0 + expected_gain / (0.0015 * 2000.0),
        );
        assert_eq!(
            solved.flow_evaluations,
            solved
                .branches
                .iter()
                .map(|b| b.solution.flow_evaluations)
                .sum::<usize>()
        );
        assert!(
            network
                .solve_pressure_driven(
                    cx,
                    &walls,
                    Pressure::new(20e6),
                    local(),
                    ParallelPressureBudget {
                        max_flow_evaluations: solved.flow_evaluations - 1,
                        ..heated_work()
                    }
                )
                .is_err()
        );
    });
}

#[test]
fn heated_material_feedback_and_mixing_keep_original_source_endpoints() {
    let response = |variable_mu| {
        with_heated_network(variable_mu, true, |cx, network| {
            let solved = network
                .solve_pressure_driven(
                    cx,
                    &[330.0, 290.0],
                    Pressure::new(20e6),
                    local(),
                    heated_work(),
                )
                .unwrap();
            let inlet_h = network.branches()[0].path.inlet().specific_enthalpy_j_kg();
            let mut endpoint_gain = 0.0;
            let mut mean_temperature = 0.0;
            for branch in &solved.branches {
                let thermal = &branch.solution.thermal;
                assert!(std::ptr::eq(
                    thermal.source,
                    solved.branches[0].solution.thermal.source
                ));
                let row = &thermal.regions[0];
                let mass = branch.solution.mass_flow_kg_per_s;
                endpoint_gain += mass * (row.outlet.specific_enthalpy_j_kg() - inlet_h);
                mean_temperature += mass * row.outlet.temperature_k();
                let bulk = 0.5 * (row.inlet.temperature_k() + row.outlet.temperature_k());
                let mu = if variable_mu {
                    0.005 - 1e-5 * bulk
                } else {
                    0.002
                };
                near(
                    branch.solution.hydraulics.regions[0]
                        .material
                        .dynamic_viscosity_pa_s(),
                    mu,
                );
                near(row.property_temperature_k, bulk);
                near(branch.solution.hydraulics.total_pressure_drop_pa, 20e6);
            }
            near(solved.total_enthalpy_gain_w, endpoint_gain);
            near(
                solved.total_enthalpy_gain_w,
                solved.total_heat_rate_w + solved.deposited_hydraulic_power_w,
            );
            near(
                solved.mixed_outlet.specific_enthalpy_j_kg(),
                inlet_h + endpoint_gain / solved.total_mass_flow_kg_per_s,
            );
            assert!(
                (solved.mixed_outlet.temperature_k()
                    - mean_temperature / solved.total_mass_flow_kg_per_s)
                    .abs()
                    > 1e-4
            );
            (
                solved.total_mass_flow_kg_per_s,
                solved.mixed_outlet.temperature_k(),
            )
        })
    };
    let fixed = response(false);
    let variable = response(true);
    assert!(variable.0 > fixed.0);
    assert!((variable.1 - fixed.1).abs() > 1e-3);
}

#[test]
fn heating_is_explicit_and_omission_allowances_cannot_mask_it() {
    with_network(false, false, 1.0, |cx, ordinary| {
        assert!(!ordinary.deposits_hydraulic_dissipation());
        let cold = ordinary
            .solve_pressure_driven(cx, &[300.0, 300.0], Pressure::new(50.0), local(), budget())
            .unwrap();
        assert_eq!(cold.deposited_hydraulic_power_w, 0.0);
        assert_eq!(cold.total_enthalpy_gain_w, 0.0);
        assert!(cold.total_dissipated_power_w > 0.0);
        let selected = ParallelLiquidNetwork::try_new(cx, ordinary.branches())
            .unwrap()
            .with_hydraulic_dissipation();
        assert!(matches!(
            selected.solve_pressure_driven(
                cx,
                &[300.0, 300.0],
                Pressure::new(50.0),
                local(),
                heated_work(),
            ),
            Err(ParallelPressureError::Invalid(_))
        )); // Original nonzero branch allowance.
    });
    with_heated_network(false, false, |cx, network| {
        assert!(matches!(
            network.solve_pressure_driven(
                cx,
                &[300.0, 300.0],
                Pressure::new(20e6),
                local(),
                budget(),
            ),
            Err(ParallelPressureError::Invalid(_))
        )); // Nonzero network allowance.
    });
}
