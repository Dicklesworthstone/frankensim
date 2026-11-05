//! Resolvable heat loads must survive a large common temperature baseline.
//! The slab references are exact nodal P1 solutions on the structured mesh;
//! these are algebraic/roundoff checks, not continuum error estimates.

mod support;

use fs_conduction::bc::{ThermalBc, ThermalBoundary, ThermalBoundaryBuilder};
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::solve::{
    ConductionProblem, ConductionSolver, InitialGuess, Nonlinearity, SolveConfig,
    element_heat_flux, solve,
};
use fs_conduction::{
    ConductionError, ConductionMesh, ConductivityModel, ConductivityTable, ScalarField,
};
use support::{with_cancelled_cx, with_cx};

fn slab() -> ConductionMesh {
    let (complex, positions) = box_grid([4, 2, 2], [1.0, 1.0, 1.0]);
    ConductionMesh::new(complex, positions).unwrap()
}

fn walls(mesh: &ConductionMesh, temperature: f64) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .region(
            "walls",
            |face| on_box_face(face.centroid[0], 0.0) || on_box_face(face.centroid[0], 1.0),
            ThermalBc::dirichlet(temperature).unwrap(),
        )
        .unwrap()
        .adiabatic_remainder()
        .finish()
        .unwrap()
}

fn drivers() -> [Nonlinearity; 2] {
    [
        Nonlinearity::FixedPoint {
            relaxation: 1.0,
            max_backtracks: 8,
        },
        Nonlinearity::default(),
    ]
}

#[test]
fn weak_internal_heating_is_resolved_at_both_temperature_baselines() {
    let mesh = slab();
    let conductivity = 10.0;
    let generation = 1e-6;
    let material = ConductivityModel::isotropic_declared(conductivity).unwrap();
    let source = ScalarField::Uniform(generation);
    // Only 12.5 nK at the midpoint: readily representable, but smaller than
    // the old residual tolerance scaled by the absolute-temperature lift.
    let peak_rise = generation / (8.0 * conductivity);
    for baseline in [300.0, 900.0] {
        let boundary = walls(&mesh, baseline);
        for nonlinearity in drivers() {
            let problem = ConductionProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material,
                element_materials: None,
                source: &source,
            };
            let mut config = SolveConfig {
                nonlinearity,
                initial: InitialGuess::Uniform(baseline),
                ..SolveConfig::default()
            };
            // The linear tolerance applies to the small heat correction. An
            // absolute-temperature Picard solve cannot deliver this accuracy
            // by subtracting two large, loosely solved temperature fields.
            config.linear.tolerance = 1e-5;
            let solution = with_cx(|cx| solve(cx, problem, config).unwrap());
            assert!(
                solution.report.iterations > 0,
                "the heater cannot disappear"
            );
            for (&temperature, position) in solution.temperature.iter().zip(mesh.positions()) {
                let x = position[0];
                let expected_rise = generation * x * (1.0 - x) / (2.0 * conductivity);
                assert!(
                    ((temperature - baseline) - expected_rise).abs() < 3e-11,
                    "baseline={baseline}, x={x}, rise={}, expected={expected_rise}",
                    temperature - baseline
                );
            }
            let measured_rise = solution
                .temperature
                .iter()
                .copied()
                .fold(baseline, f64::max)
                - baseline;
            assert!(measured_rise > 0.99 * peak_rise);
            assert!(solution.report.final_residual <= solution.report.residual_threshold);
            assert!(
                solution.report.residual_threshold < 0.1 * solution.report.residual_history[0],
                "the explicit floating-point allowance must still resolve this heat load"
            );
            assert!(
                (solution.report.energy.dirichlet_in_w + generation).abs() < 2e-9,
                "the walls must remove the independently specified generated heat"
            );
        }
    }
}

#[test]
fn weak_neumann_heat_reaches_a_robin_reservoir_without_dirichlet_rows() {
    let mesh = slab();
    let conductivity = 10.0;
    let htc = 5.0;
    let inflow = 1e-7;
    let ambient = 600.0;
    let material = ConductivityModel::isotropic_declared(conductivity).unwrap();
    let source = ScalarField::Uniform(0.0);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "heater",
            |face| on_box_face(face.centroid[0], 0.0),
            ThermalBc::neumann(-inflow).unwrap(),
        )
        .unwrap()
        .region(
            "reservoir",
            |face| on_box_face(face.centroid[0], 1.0),
            ThermalBc::robin(htc, ambient).unwrap(),
        )
        .unwrap()
        .adiabatic_remainder()
        .finish()
        .unwrap();
    for nonlinearity in drivers() {
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let config = SolveConfig {
            nonlinearity,
            initial: InitialGuess::Uniform(ambient),
            ..SolveConfig::default()
        };
        let solution = with_cx(|cx| solve(cx, problem, config).unwrap());
        assert!(solution.report.iterations > 0);
        for (&temperature, position) in solution.temperature.iter().zip(mesh.positions()) {
            let expected_rise = inflow * ((1.0 - position[0]) / conductivity + 1.0 / htc);
            assert!(((temperature - ambient) - expected_rise).abs() < 3e-11);
        }
        assert!(solution.report.final_residual <= solution.report.residual_threshold);
        assert!((solution.report.energy.neumann_out_w + inflow).abs() < 1e-20);
        assert!(
            (solution.report.energy.robin_out_w - inflow).abs() < 1e-10,
            "the Robin reservoir must receive the weak prescribed heat"
        );
    }
}

#[test]
fn weak_robin_reservoir_cannot_inflate_the_reference_heat_scale() {
    let mesh = slab();
    let conductivity = 10.0;
    let material = ConductivityModel::isotropic_declared(conductivity).unwrap();
    let source = ScalarField::Uniform(0.0);
    let cold = 300.0;
    let hot = 600.0;
    let weak_htc = 1e-9;
    // With a prescribed cold wall the weak bath produces a 30 nK rise.
    // Also exercise the Robin-only case: its strong reservoir must dominate
    // the reference through conductance, not share equal weight with the bath.
    for cold_htc in [None, Some(10.0)] {
        let cold_rule = match cold_htc {
            Some(htc) => ThermalBc::robin(htc, cold).unwrap(),
            None => ThermalBc::dirichlet(cold).unwrap(),
        };
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .region("cold", |face| on_box_face(face.centroid[0], 0.0), cold_rule)
            .unwrap()
            .region(
                "weak hot bath",
                |face| on_box_face(face.centroid[0], 1.0),
                ThermalBc::robin(weak_htc, hot).unwrap(),
            )
            .unwrap()
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let cold_resistance = cold_htc.map_or(0.0, |htc| 1.0 / htc);
        let heat_in = (hot - cold) / (1.0 / weak_htc + 1.0 / conductivity + cold_resistance);
        for nonlinearity in drivers() {
            let solution = with_cx(|cx| {
                solve(
                    cx,
                    ConductionProblem {
                        mesh: &mesh,
                        boundary: &boundary,
                        material: &material,
                        element_materials: None,
                        source: &source,
                    },
                    SolveConfig {
                        nonlinearity,
                        initial: InitialGuess::Uniform(cold),
                        ..SolveConfig::default()
                    },
                )
                .unwrap()
            });
            assert!(
                solution.report.iterations > 0,
                "weak Robin heat must not disappear; cold h={cold_htc:?}"
            );
            for (&temperature, position) in solution.temperature.iter().zip(mesh.positions()) {
                let expected_rise = heat_in * (cold_resistance + position[0] / conductivity);
                assert!(
                    ((temperature - cold) - expected_rise).abs() < 3e-11,
                    "cold h={cold_htc:?}, x={}, rise={}, expected={expected_rise}",
                    position[0],
                    temperature - cold
                );
            }
            assert!(solution.report.final_residual <= solution.report.residual_threshold);
            assert!(
                solution.report.residual_threshold < 0.1 * solution.report.residual_history[0],
                "a weak reservoir must not inflate the relative convergence target"
            );
            let incoming = solution
                .report
                .robin_fluxes
                .iter()
                .find(|flux| flux.region == "weak hot bath")
                .unwrap();
            assert!((incoming.heat_rate_w + heat_in).abs() < 1e-10);
        }
    }
}

#[test]
fn skew_anisotropic_isothermal_equilibrium_does_not_chase_roundoff() {
    let (complex, mut positions) = box_grid([3, 2, 2], [1.0, 1.0, 1.0]);
    for point in &mut positions {
        *point = [
            point[0] + 0.17 * point[1],
            point[1] + 0.23 * point[2],
            0.7 * point[2],
        ];
    }
    let mesh = ConductionMesh::new(complex, positions).unwrap();
    let material =
        ConductivityModel::constant_tensor([[10.0, 1.0, 0.5], [1.0, 4.0, 0.25], [0.5, 0.25, 2.0]])
            .unwrap();
    let source = ScalarField::Uniform(0.0);
    let temperature = 315.125;
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "reservoir",
            |_| true,
            ThermalBc::robin(8.0, temperature).unwrap(),
        )
        .unwrap()
        .finish()
        .unwrap();
    for nonlinearity in drivers() {
        let solution = with_cx(|cx| {
            solve(
                cx,
                ConductionProblem {
                    mesh: &mesh,
                    boundary: &boundary,
                    material: &material,
                    element_materials: None,
                    source: &source,
                },
                SolveConfig {
                    nonlinearity,
                    initial: InitialGuess::Uniform(temperature),
                    ..SolveConfig::default()
                },
            )
            .unwrap()
        });
        assert_eq!(solution.report.iterations, 0);
        assert!(
            solution
                .temperature
                .iter()
                .all(|&value| value == temperature)
        );
        assert!(
            solution.report.residual_threshold > 1e-30,
            "the reported threshold must expose the numerical floor"
        );
        assert!(solution.report.final_residual <= solution.report.residual_threshold);
        assert_eq!(solution.report.energy.robin_out_w, 0.0);
    }
}

#[test]
fn a_small_damped_step_cannot_publish_unbalanced_heat_as_converged() {
    let mesh = slab();
    let boundary = walls(&mesh, 300.0);
    let material = ConductivityModel::isotropic_declared(10.0).unwrap();
    let source = ScalarField::Uniform(100.0);
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &material,
        element_materials: None,
        source: &source,
    };
    let mut config = SolveConfig {
        nonlinearity: Nonlinearity::FixedPoint {
            relaxation: 1e-6,
            max_backtracks: 8,
        },
        ..SolveConfig::default()
    };
    config.stop.step_atol = 1.0;
    with_cx(|cx| {
        let mut solver = ConductionSolver::new(problem, config).unwrap();
        let first = solver.step(cx).unwrap();
        assert!(!first.converged);
        assert!(first.step_norm < 1.0);
        let before = solver.state().clone();
        match solver.step(cx).unwrap_err() {
            ConductionError::NotConverged {
                iterations,
                residual,
                threshold,
            } => {
                assert_eq!(iterations, 1);
                assert!(residual > 0.99 * first.residual);
                assert!(residual > threshold);
            }
            other => panic!("expected unresolved heat, got {other}"),
        }
        assert_eq!(solver.state(), &before);
    });
}

#[test]
fn fully_prescribed_anisotropic_field_retains_flux_and_every_heat_balance_term() {
    let (complex, positions) = box_grid([1, 1, 1], [1.0, 1.0, 1.0]);
    let mesh = ConductionMesh::new(complex, positions).unwrap();
    let prescribed: Vec<_> = mesh
        .positions()
        .iter()
        .map(|p| 300.0 + 20.0 * p[0] - 10.0 * p[1] + 5.0 * p[2])
        .collect();
    // Both x faces together prescribe every vertex of this coarse mesh.
    // The other faces can still exchange heat with external reservoirs.
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "temperature walls",
            |face| on_box_face(face.centroid[0], 0.0) || on_box_face(face.centroid[0], 1.0),
            ThermalBc::Dirichlet {
                temperature: ScalarField::Nodal(prescribed.clone()),
            },
        )
        .unwrap()
        .region(
            "heater",
            |face| on_box_face(face.centroid[1], 0.0),
            ThermalBc::neumann(-4.0).unwrap(),
        )
        .unwrap()
        .region(
            "coolant",
            |face| on_box_face(face.centroid[1], 1.0),
            ThermalBc::robin(5.0, 295.0).unwrap(),
        )
        .unwrap()
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let material =
        ConductivityModel::constant_tensor([[10.0, 1.0, 0.5], [1.0, 4.0, 0.25], [0.5, 0.25, 2.0]])
            .unwrap();
    let source = ScalarField::Uniform(120.0);
    for nonlinearity in drivers() {
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let config = SolveConfig {
            nonlinearity,
            initial: InitialGuess::Free(Vec::new()),
            ..SolveConfig::default()
        };
        assert_eq!(
            ConductionSolver::initial_temperature(problem, &config).unwrap(),
            prescribed
        );
        let solution = with_cx(|cx| solve(cx, problem, config).unwrap());
        assert_eq!(solution.temperature, prescribed);
        assert_eq!(solution.report.free_dofs, 0);
        assert_eq!(solution.report.iterations, 0);
        assert!(solution.report.linear.is_empty());
        assert_eq!(solution.report.final_residual, 0.0);
        let energy = &solution.report.energy;
        assert!((energy.source_w - 120.0).abs() < 1e-12);
        assert!((energy.neumann_out_w + 4.0).abs() < 1e-12);
        assert!((energy.robin_out_w - 37.5).abs() < 1e-12);
        assert!((energy.dirichlet_in_w + 86.5).abs() < 1e-10);
        assert!(energy.closure_w.abs() < 1e-10);
        for flux in element_heat_flux(&mesh, &material, &solution.temperature).unwrap() {
            for (actual, expected) in flux.into_iter().zip([-192.5, 18.75, -17.5]) {
                assert!((actual - expected).abs() < 1e-12);
            }
        }
    }
}

#[test]
fn fully_prescribed_field_uses_actual_nonlinear_conductivity_and_keeps_admission_checks() {
    let (complex, positions) = box_grid([1, 1, 1], [1.0, 1.0, 1.0]);
    let mesh = ConductionMesh::new(complex, positions).unwrap();
    let prescribed: Vec<_> = mesh
        .positions()
        .iter()
        .map(|p| 300.0 + 20.0 * p[0])
        .collect();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .remainder(
            "temperature walls",
            ThermalBc::Dirichlet {
                temperature: ScalarField::Nodal(prescribed.clone()),
            },
        )
        .unwrap()
        .finish()
        .unwrap();
    let material = ConductivityModel::isotropic(
        ConductivityTable::declared_curve(vec![(300.0, 2.0), (320.0, 2.2)]).unwrap(),
    );
    let source = ScalarField::Nodal(mesh.positions().iter().map(|p| 7.0 + 2.0 * p[0]).collect());
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &material,
        element_materials: None,
        source: &source,
    };
    let solution = with_cx(|cx| solve(cx, problem, SolveConfig::default()).unwrap());
    assert_eq!(solution.temperature, prescribed);
    assert!((solution.report.energy.source_w - 8.0).abs() < 1e-12);
    assert!((solution.report.energy.dirichlet_in_w + 8.0).abs() < 1e-11);
    let fluxes = element_heat_flux(&mesh, &material, &solution.temperature).unwrap();
    for (tet, flux) in mesh.complex().tets.iter().zip(fluxes) {
        let mean_x = tet
            .iter()
            .map(|&v| mesh.positions()[v as usize][0])
            .sum::<f64>()
            / 4.0;
        assert!((flux[0] + 20.0 * (2.0 + 0.2 * mean_x)).abs() < 1e-12);
        assert!(flux[1].abs() < 1e-12 && flux[2].abs() < 1e-12);
    }
    with_cancelled_cx(|cx| {
        assert!(matches!(
            solve(cx, problem, SolveConfig::default()),
            Err(ConductionError::Cancelled { .. })
        ));
    });
    let limited_material = ConductivityModel::isotropic(
        ConductivityTable::declared_curve(vec![(300.0, 2.0), (310.0, 2.1)]).unwrap(),
    );
    with_cx(|cx| {
        assert!(matches!(
            solve(
                cx,
                ConductionProblem {
                    material: &limited_material,
                    ..problem
                },
                SolveConfig::default()
            ),
            Err(ConductionError::OutsideTemperatureSpan { .. })
        ));
    });
}
