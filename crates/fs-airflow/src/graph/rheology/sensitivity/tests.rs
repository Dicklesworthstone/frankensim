//! Independent Poiseuille algebra and perturbed ORIGINAL network solves.
//! Finite differences are test-only; production uses retained constitutive
//! derivatives and one grounded linear solve per direction or objective.

use super::*;
use crate::graph::rheology::CarreauPipeBranch;
use crate::graph::tests::with_cx;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_material::fluid_transport::rheology::{
    CarreauYasuda, ThermalCarreauYasuda,
    pipe::{CarreauPipeControls, LaminarCarreauPipe},
};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};

fn law(mode: ViscosityShiftMode, newtonian: bool) -> ThermalCarreauYasuda {
    ThermalCarreauYasuda::new(
        CarreauYasuda::new(
            if newtonian { 3.0 } else { 12.0 },
            if newtonian { 3.0 } else { 0.2 },
            0.6,
            0.45,
            2.0,
            [0.0, 1e7],
        )
        .unwrap(),
        TemperatureShift::new(
            TemperatureShiftModel::Arrhenius {
                activation_energy_j_mol: 12_000.0,
            },
            300.0,
            [280.0, 360.0],
        )
        .unwrap(),
        mode,
    )
}
fn pipe(
    law: ThermalCarreauYasuda,
    temperature: f64,
    radius: f64,
    length: f64,
) -> LaminarCarreauPipe {
    LaminarCarreauPipe::new(law, temperature, radius, length, 1e4).unwrap()
}
fn branch(name: &str, from: usize, to: usize, pipe: LaminarCarreauPipe) -> CarreauPipeBranch {
    CarreauPipeBranch {
        name: name.into(),
        from,
        to,
        pipe,
    }
}
fn boundary(node: usize, pressure: f64) -> FixedPressure {
    FixedPressure {
        node,
        pressure: Pressure::new(pressure),
    }
}
fn linear() -> CheckedCgConfig {
    CheckedCgConfig {
        tolerance: 2e-12,
        max_iterations: 128,
        max_corrections: 2,
    }
}
fn primal() -> CarreauGraphControls {
    CarreauGraphControls {
        absolute_flow_tolerance: VolumetricFlowRate::new(2e-16),
        relative_flow_tolerance: 1e-11,
        pipe: CarreauPipeControls {
            relative_constitutive_tolerance: 1e-12,
            relative_quadrature_tolerance: 1e-11,
            ..CarreauPipeControls::default()
        },
        ..CarreauGraphControls::default()
    }
}
fn near(a: f64, b: f64, relative: f64, absolute: f64) {
    assert!(
        a.is_finite()
            && b.is_finite()
            && (a - b).abs() <= absolute + relative * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}
fn bridge(mode: ViscosityShiftMode) -> CarreauPipeGraph {
    let source = law(mode, false);
    CarreauPipeGraph::new(
        4,
        vec![
            branch("feed", 0, 1, pipe(source, 295.0, 0.010, 0.5)),
            branch("left", 1, 3, pipe(source, 310.0, 0.013, 1.0)),
            branch("cross", 2, 1, pipe(source, 335.0, 0.009, 0.75)),
            branch("right-feed", 0, 2, pipe(source, 315.0, 0.012, 2.0)),
            branch("right", 2, 3, pipe(source, 305.0, 0.014, 0.8)),
        ],
    )
    .unwrap()
}
fn ports() -> (Vec<FixedPressure>, Vec<FlowInjection>) {
    (
        vec![boundary(0, 3000.0), boundary(3, 50.0)],
        vec![
            FlowInjection {
                node: 1,
                flow: VolumetricFlowRate::new(2e-6),
            },
            FlowInjection {
                node: 2,
                flow: VolumetricFlowRate::new(-1e-6),
            },
            FlowInjection {
                node: 1,
                flow: VolumetricFlowRate::new(0.5e-6),
            },
            FlowInjection {
                node: 0,
                flow: VolumetricFlowRate::new(-3e-7),
            },
        ],
    )
}
fn perturbed(
    graph: &CarreauPipeGraph,
    boundaries: &[FixedPressure],
    injections: &[FlowInjection],
    d: CarreauGraphDirection<'_>,
    h: f64,
) -> CarreauGraphSolution {
    let branches = graph
        .branches()
        .iter()
        .enumerate()
        .map(|(i, b)| {
            branch(
                &b.name,
                b.from,
                b.to,
                LaminarCarreauPipe::new(
                    *b.pipe.law(),
                    b.pipe.temperature_k() + h * d.temperatures_k[i],
                    b.pipe.radius_m() * (h * d.log_radii[i]).exp(),
                    b.pipe.length_m() * (h * d.log_lengths[i]).exp(),
                    b.pipe.maximum_laminar_wall_rate_per_s(),
                )
                .unwrap(),
            )
        })
        .collect();
    let changed = CarreauPipeGraph::new(graph.node_count(), branches).unwrap();
    let boundaries: Vec<_> = boundaries
        .iter()
        .zip(d.boundary_pressures)
        .map(|(b, v)| boundary(b.node, b.pressure.value() + h * v.value()))
        .collect();
    let injections: Vec<_> = injections
        .iter()
        .zip(d.injections)
        .map(|(p, v)| FlowInjection {
            node: p.node,
            flow: VolumetricFlowRate::new(p.flow.value() + h * v.value()),
        })
        .collect();
    with_cx(|cx| changed.solve_with_injections(&boundaries, &injections, primal(), cx)).unwrap()
}

#[test]
fn g1_original_geometry_and_both_temperature_shifts_match_resolved_network_differences() {
    let r = [0.3, -0.1, 0.2, 0.15, -0.05];
    let l = [-0.1, 0.25, 0.0, -0.2, 0.1];
    let t = [2.0, -1.0, 3.0, -0.5, 1.5];
    let bp = [Pressure::new(27.0), Pressure::new(-13.0)];
    let ip = [1e-6, -0.7e-6, 0.2e-6, -0.8e-6].map(VolumetricFlowRate::new);
    for mode in [
        ViscosityShiftMode::Vertical,
        ViscosityShiftMode::VerticalAndRate,
    ] {
        let graph = bridge(mode);
        let (boundaries, injections) = ports();
        let owner = with_cx(|cx| {
            graph.solve_with_injections_linearized(&boundaries, &injections, primal(), linear(), cx)
        })
        .unwrap();
        for family in 0..4 {
            // Three isolated physical directions plus one simultaneous
            // geometry/thermal/pressure/injection perturbation.
            let direction = CarreauGraphDirection {
                log_radii: if family == 0 || family == 3 {
                    &r
                } else {
                    &[0.0; 5]
                },
                log_lengths: if family == 1 || family == 3 {
                    &l
                } else {
                    &[0.0; 5]
                },
                temperatures_k: if family == 2 || family == 3 {
                    &t
                } else {
                    &[0.0; 5]
                },
                boundary_pressures: if family == 3 {
                    &bp
                } else {
                    &[Pressure::new(0.0); 2]
                },
                injections: if family == 3 {
                    &ip
                } else {
                    &[VolumetricFlowRate::new(0.0); 4]
                },
            };
            let tangent = with_cx(|cx| owner.tangent(direction, cx)).unwrap();
            let h = 1e-4;
            let plus = perturbed(&graph, &boundaries, &injections, direction, h);
            let minus = perturbed(&graph, &boundaries, &injections, direction, -h);
            for i in 0..graph.node_count() {
                near(
                    tangent.pressures[i].value(),
                    (plus.pressures[i].value() - minus.pressures[i].value()) / (2.0 * h),
                    3e-5,
                    2e-6,
                );
            }
            for i in 0..graph.branches().len() {
                near(
                    tangent.flows[i].value(),
                    (plus.branches[i].profile.flow().volumetric_flow_m3_s
                        - minus.branches[i].profile.flow().volumetric_flow_m3_s)
                        / (2.0 * h),
                    3e-5,
                    2e-12,
                );
            }
            for i in 0..boundaries.len() {
                near(
                    tangent.reservoir_outflows[i].value(),
                    (plus.reservoir_outflows[i].value() - minus.reservoir_outflows[i].value())
                        / (2.0 * h),
                    3e-5,
                    3e-12,
                );
            }
            assert!(tangent.max_node_imbalance.value() < 2e-16);
            assert!(tangent.linear_report.converged_euclidean());
        }
        for (edge, original) in owner.solution().branches.iter().zip(graph.branches()) {
            assert_eq!(edge.profile.pipe(), &original.pipe);
        }
    }
}

#[test]
fn g1_g3_adjoint_dot_product_includes_repeated_and_fixed_pressure_injections() {
    let r = [0.3, -0.1, 0.2, 0.15, -0.05];
    let l = [-0.1, 0.25, 0.0, -0.2, 0.1];
    let t = [2.0, -1.0, 3.0, -0.5, 1.5];
    let bp = [Pressure::new(27.0), Pressure::new(-13.0)];
    let ip = [1e-6, -0.7e-6, 0.2e-6, -0.8e-6].map(VolumetricFlowRate::new);
    let wp = [3e-4, -4e-4, 2e-4, 1.5e-4].map(PressureWeight::new);
    let wq = [1e4, -2e4, 3e4, -4e4, 5e4].map(FlowWeight::new);
    for mode in [
        ViscosityShiftMode::Vertical,
        ViscosityShiftMode::VerticalAndRate,
    ] {
        let graph = bridge(mode);
        let (boundaries, injections) = ports();
        let owner = with_cx(|cx| {
            graph.solve_with_injections_linearized(&boundaries, &injections, primal(), linear(), cx)
        })
        .unwrap();
        let tangent = with_cx(|cx| {
            owner.tangent(
                CarreauGraphDirection {
                    log_radii: &r,
                    log_lengths: &l,
                    temperatures_k: &t,
                    boundary_pressures: &bp,
                    injections: &ip,
                },
                cx,
            )
        })
        .unwrap();
        let gradient = with_cx(|cx| owner.adjoint(&wp, &wq, cx)).unwrap();
        let forward = wp
            .iter()
            .zip(&tangent.pressures)
            .map(|(w, p)| w.value() * p.value())
            .sum::<f64>()
            + wq.iter()
                .zip(&tangent.flows)
                .map(|(w, q)| w.value() * q.value())
                .sum::<f64>();
        let reverse = gradient
            .log_radii
            .iter()
            .zip(r)
            .map(|(g, d)| g * d)
            .sum::<f64>()
            + gradient
                .log_lengths
                .iter()
                .zip(l)
                .map(|(g, d)| g * d)
                .sum::<f64>()
            + gradient
                .temperatures_per_k
                .iter()
                .zip(t)
                .map(|(g, d)| g * d)
                .sum::<f64>()
            + gradient
                .boundary_pressures
                .iter()
                .zip(bp)
                .map(|(g, d)| g.value() * d.value())
                .sum::<f64>()
            + gradient
                .injections
                .iter()
                .zip(ip)
                .map(|(g, d)| g.value() * d.value())
                .sum::<f64>();
        near(forward, reverse, 1e-10, 1e-12);
        assert_eq!(gradient.injections[0], gradient.injections[2]);
        assert_eq!(gradient.injections[3].value(), 0.0);
        assert!(gradient.linear_report.converged_euclidean());
        let pressure_direction = CarreauGraphDirection {
            log_radii: &[0.0; 5],
            log_lengths: &[0.0; 5],
            temperatures_k: &[0.0; 5],
            boundary_pressures: &bp,
            injections: &[VolumetricFlowRate::new(0.0); 4],
        };
        let ordinary = with_cx(|cx| owner.tangent(pressure_direction, cx)).unwrap();
        let shifted = with_cx(|cx| {
            owner.tangent(
                CarreauGraphDirection {
                    boundary_pressures: &[Pressure::new(1e12 + 27.0), Pressure::new(1e12 - 13.0)],
                    ..pressure_direction
                },
                cx,
            )
        })
        .unwrap();
        assert_eq!(ordinary.flows, shifted.flows);
        assert_eq!(ordinary.reservoir_outflows, shifted.reservoir_outflows);
        for (a, b) in shifted.pressures.iter().zip(&ordinary.pressures) {
            near(a.value() - 1e12, b.value(), 0.0, 2e-4);
        }
        let fixed_only = with_cx(|cx| {
            owner.tangent(
                CarreauGraphDirection {
                    log_radii: &[0.0; 5],
                    log_lengths: &[0.0; 5],
                    temperatures_k: &[0.0; 5],
                    boundary_pressures: &[Pressure::new(0.0); 2],
                    injections: &[
                        VolumetricFlowRate::new(0.0),
                        VolumetricFlowRate::new(0.0),
                        VolumetricFlowRate::new(0.0),
                        VolumetricFlowRate::new(1.0),
                    ],
                },
                cx,
            )
        })
        .unwrap();
        assert!(fixed_only.pressures.iter().all(|p| p.value() == 0.0));
        assert!(fixed_only.flows.iter().all(|q| q.value() == 0.0));
        assert_eq!(fixed_only.reservoir_outflows[0].value(), -1.0);
    }
}

#[test]
fn g1_newtonian_limit_matches_independent_poiseuille_geometry_and_thermal_derivatives() {
    let source = law(ViscosityShiftMode::VerticalAndRate, true);
    let temperatures = [300.0, 330.0];
    let radii = [0.01_f64, 0.02];
    let lengths = [0.5, 0.7];
    let graph = CarreauPipeGraph::new(
        3,
        vec![
            branch(
                "a",
                0,
                1,
                pipe(source, temperatures[0], radii[0], lengths[0]),
            ),
            branch(
                "b",
                1,
                2,
                pipe(source, temperatures[1], radii[1], lengths[1]),
            ),
        ],
    )
    .unwrap();
    let boundaries = [boundary(0, 800.0), boundary(2, 50.0)];
    let inputs = [FlowInjection {
        node: 1,
        flow: VolumetricFlowRate::new(3e-6),
    }];
    let owner = with_cx(|cx| {
        graph.solve_with_injections_linearized(&boundaries, &inputs, primal(), linear(), cx)
    })
    .unwrap();
    let r = [0.3, -0.2];
    let l = [-0.1, 0.15];
    let t = [2.0, -1.0];
    let tangent = with_cx(|cx| {
        owner.tangent(
            CarreauGraphDirection {
                log_radii: &r,
                log_lengths: &l,
                temperatures_k: &t,
                boundary_pressures: &[Pressure::new(27.0), Pressure::new(-13.0)],
                injections: &[VolumetricFlowRate::new(0.8e-6)],
            },
            cx,
        )
    })
    .unwrap();
    // This oracle uses only Poiseuille, the declared Arrhenius algebra, and a
    // one-unknown Kirchhoff equation, not material or native derivative APIs.
    let activation = 12_000.0 / 8.314_462_618_153_24;
    let g: [f64; 2] = core::array::from_fn(|i| {
        let mu = 3.0 * (activation * (1.0 / temperatures[i] - 1.0 / 300.0)).exp();
        core::f64::consts::PI * radii[i].powi(4) / (8.0 * lengths[i] * mu)
    });
    let p = (g[0] * 800.0 + g[1] * 50.0 + 3e-6) / (g[0] + g[1]);
    let dg: [f64; 2] = core::array::from_fn(|i| {
        g[i] * (4.0 * r[i] - l[i] + activation / temperatures[i].powi(2) * t[i])
    });
    let dp = (g[0] * 27.0 - g[1] * 13.0 + 0.8e-6 + dg[0] * (800.0 - p) + dg[1] * (50.0 - p))
        / (g[0] + g[1]);
    near(tangent.pressures[1].value(), dp, 2e-12, 2e-12);
    near(
        tangent.flows[0].value(),
        dg[0] * (800.0 - p) + g[0] * (27.0 - dp),
        2e-12,
        1e-20,
    );
    near(
        tangent.flows[1].value(),
        dg[1] * (p - 50.0) + g[1] * (dp + 13.0),
        2e-12,
        1e-20,
    );
}

#[test]
fn g3_rest_and_separate_component_gauges_keep_the_original_finite_conductance() {
    let source = law(ViscosityShiftMode::Vertical, false);
    let p = pipe(source, 300.0, 0.01, 0.5);
    let graph = CarreauPipeGraph::new(
        6,
        vec![
            branch("a", 0, 1, p),
            branch("b", 1, 2, p),
            branch("c", 3, 4, p),
            branch("d", 4, 5, p),
        ],
    )
    .unwrap();
    let boundaries = [
        boundary(0, 100_000.0),
        boundary(2, 100_000.0),
        boundary(3, -7.0),
        boundary(5, -7.0),
    ];
    let owner = with_cx(|cx| graph.solve_linearized(&boundaries, primal(), linear(), cx)).unwrap();
    let tangent = with_cx(|cx| {
        owner.tangent(
            CarreauGraphDirection {
                log_radii: &[1.0; 4],
                log_lengths: &[1.0; 4],
                temperatures_k: &[1.0; 4],
                boundary_pressures: &[
                    Pressure::new(1.0),
                    Pressure::new(0.0),
                    Pressure::new(-2.0),
                    Pressure::new(-2.0),
                ],
                injections: &[],
            },
            cx,
        )
    })
    .unwrap();
    let g = core::f64::consts::PI * 0.01_f64.powi(4) / (8.0 * 12.0 * 0.5);
    for (actual, expected) in tangent
        .pressures
        .iter()
        .zip([1.0, 0.5, 0.0, -2.0, -2.0, -2.0])
    {
        near(actual.value(), expected, 2e-12, 1e-14);
    }
    for (actual, expected) in tangent.flows.iter().zip([0.5 * g, 0.5 * g, 0.0, 0.0]) {
        near(actual.value(), expected, 2e-12, 1e-22);
    }
    let gradient = with_cx(|cx| {
        owner.adjoint(
            &[PressureWeight::new(0.0); 6],
            &[FlowWeight::new(1.0); 4],
            cx,
        )
    })
    .unwrap();
    assert!(
        gradient
            .log_radii
            .iter()
            .chain(&gradient.log_lengths)
            .chain(&gradient.temperatures_per_k)
            .all(|x| *x == 0.0)
    );
}

#[test]
fn g0_source_domain_boundaries_do_not_become_two_sided_derivative_authority() {
    let source = law(ViscosityShiftMode::Vertical, false);
    let at_endpoint = CarreauPipeGraph::new(
        2,
        vec![branch("endpoint", 0, 1, pipe(source, 280.0, 0.01, 0.5))],
    )
    .unwrap();
    let boundaries = [boundary(0, 0.0), boundary(1, 0.0)];
    assert!(matches!(
        with_cx(|cx| at_endpoint.solve_linearized(&boundaries, primal(), linear(), cx)),
        Err(CarreauGraphError::Graph(GraphError::InvalidBranch {
            branch: 0,
            ..
        }))
    ));
    let limited = LaminarCarreauPipe::new(source, 300.0, 1.0, 0.5, 1.0).unwrap();
    let ceiling_pressure = limited
        .profile_at_wall_rate(1.0, primal().pipe)
        .unwrap()
        .flow()
        .pressure_drop_pa;
    let at_rate = CarreauPipeGraph::new(2, vec![branch("rate", 0, 1, limited)]).unwrap();
    assert!(matches!(
        with_cx(|cx| at_rate.solve_linearized(
            &[boundary(0, ceiling_pressure), boundary(1, 0.0)],
            primal(),
            linear(),
            cx
        )),
        Err(CarreauGraphError::Graph(GraphError::InvalidBranch {
            branch: 0,
            ..
        }))
    ));
    for mode in [
        ViscosityShiftMode::Vertical,
        ViscosityShiftMode::VerticalAndRate,
    ] {
        let p = LaminarCarreauPipe::new(law(mode, false), 300.0, 0.01, 0.5, 1e7).unwrap();
        let graph = CarreauPipeGraph::new(2, vec![branch("envelope", 0, 1, p)]).unwrap();
        let result = with_cx(|cx| graph.solve_linearized(&boundaries, primal(), linear(), cx));
        // Only the rate-shift convention moves the full prepared envelope.
        assert_eq!(result.is_ok(), mode == ViscosityShiftMode::Vertical);
    }
}

#[test]
fn g4_bad_directions_work_limits_and_cancellation_preserve_the_accepted_owner() {
    let graph = bridge(ViscosityShiftMode::Vertical);
    let (boundaries, inputs) = ports();
    let owner = with_cx(|cx| {
        graph.solve_with_injections_linearized(&boundaries, &inputs, primal(), linear(), cx)
    })
    .unwrap();
    let before_pressures = owner.solution().pressures.clone();
    let before_branches = owner.solution().branches.clone();
    let d = CarreauGraphDirection {
        log_radii: &[0.2; 5],
        log_lengths: &[-0.1; 5],
        temperatures_k: &[1.0; 5],
        boundary_pressures: &[Pressure::new(2.0), Pressure::new(-1.0)],
        injections: &[VolumetricFlowRate::new(1e-6); 4],
    };
    assert!(
        with_cx(|cx| owner.tangent(
            CarreauGraphDirection {
                log_radii: &[],
                ..d
            },
            cx
        ))
        .is_err()
    );
    assert!(
        with_cx(|cx| owner.tangent(
            CarreauGraphDirection {
                temperatures_k: &[f64::NAN; 5],
                ..d
            },
            cx
        ))
        .is_err()
    );
    assert!(with_cx(|cx| owner.adjoint(&[], &[], cx)).is_err());
    assert!(
        with_cx(|cx| owner.adjoint(
            &[PressureWeight::new(f64::NAN); 4],
            &[FlowWeight::new(0.0); 5],
            cx
        ))
        .is_err()
    );
    assert!(matches!(
        with_cx(|cx| graph.solve_linearized(
            &boundaries,
            primal(),
            CheckedCgConfig {
                max_iterations: 0,
                ..linear()
            },
            cx
        )),
        Err(CarreauGraphError::Graph(GraphError::InvalidInput(_)))
    ));
    let short = with_cx(|cx| {
        graph.solve_with_injections_linearized(
            &boundaries,
            &inputs,
            primal(),
            CheckedCgConfig {
                max_iterations: 1,
                max_corrections: 0,
                tolerance: 1e-14,
            },
            cx,
        )
    })
    .unwrap();
    assert!(matches!(
        with_cx(|cx| short.tangent(d, cx)),
        Err(CarreauGraphError::LinearDidNotConverge(_))
    ));
    let gate = CancelGate::new();
    gate.request();
    let cancelled = ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 17,
                kernel_id: 711,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        owner.tangent(d, &cx)
    });
    assert!(matches!(
        cancelled,
        Err(CarreauGraphError::Graph(GraphError::Interrupted))
    ));
    assert_eq!(owner.solution().pressures, before_pressures);
    assert_eq!(owner.solution().branches, before_branches);
    let first = with_cx(|cx| owner.tangent(d, cx)).unwrap();
    let retry = with_cx(|cx| owner.tangent(d, cx)).unwrap();
    assert_eq!(first.pressures, retry.pressures);
    assert_eq!(first.flows, retry.flows);
}
