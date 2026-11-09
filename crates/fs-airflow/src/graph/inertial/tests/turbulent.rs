//! Original Colebrook drag on the native, energy-audited momentum trajectory.

use super::*;
use crate::graph::TurbulentPipeLoss;
use crate::graph::thermal::medium::series::hydraulic::{MaterialFlowResistance, ResistanceError};
use fs_qty::Length;

const RHO: f64 = 1000.0;
const DIAMETER: f64 = 0.01;
const LENGTH: f64 = 1.0;
const LINEAR: f64 = 3e7;

fn area() -> f64 {
    core::f64::consts::FRAC_PI_4 * DIAMETER * DIAMETER
}

fn pipe(mu: f64, relative_roughness: f64) -> TurbulentPipeLoss {
    TurbulentPipeLoss::new(
        "original-pipe",
        MaterialFlowResistance::circular_turbulent(
            Length::new(LENGTH),
            Length::new(DIAMETER),
            Length::new(relative_roughness * DIAMETER),
            (4000.0, 1e6),
        )
        .unwrap(),
        FluidTransportState::try_new(RHO, 4000.0, 0.6, mu).unwrap(),
        RHO,
        0.0,
        SourceProvenance::new("Colebrook (1939)", "native-inertial-pipe-fixture"),
        ToleranceBasis::EngineeringAllowance,
    )
    .unwrap()
}

fn network(mu: f64, roughness: f64, series: bool) -> LossGraph {
    let mut branches = vec![PassiveBranch::Turbulent {
        from: 0,
        to: 1,
        loss: pipe(mu, roughness),
    }];
    if series {
        branches.push(PassiveBranch::Linear {
            from: 1,
            to: 2,
            loss: LinearLossElement::new(
                "original-linear",
                LinearLossResistance::new(LINEAR),
                0.0,
                SourceProvenance::new("authored linear fixture", "native-inertial-series"),
                ToleranceBasis::Analytic,
            )
            .unwrap(),
        });
    }
    LossGraph::new_mixed(if series { 3 } else { 2 }, branches).unwrap()
}

fn endpoints(graph: &LossGraph) -> [FixedPressure; 2] {
    [boundary(0, 0.0), boundary(graph.node_count() - 1, 0.0)]
}

fn state(cx: &Cx<'_>, graph: &LossGraph, q: f64) -> InertialGraph {
    let count = graph.branches().len();
    InertialGraph::try_new(
        cx,
        graph,
        Density::new(RHO),
        &vec![
            InertialBranch {
                length_m: LENGTH,
                flow_area_m2: area()
            };
            count
        ],
        &endpoints(graph),
        &vec![q; count],
        0.0,
        16,
    )
    .unwrap()
}

// Formula-disjoint reference: bisection in inverse sqrt(f), platform log10,
// then Darcy-Weisbach. Production uses safeguarded Newton and det::ln.
fn pressure(q: f64, mu: f64, roughness: f64) -> f64 {
    let reynolds = RHO * q.abs() * DIAMETER / (mu * area());
    let (mut low, mut high) = (1.0, 32.0);
    for _ in 0..96 {
        let x = 0.5 * low + 0.5 * high;
        if x + 2.0 * (roughness / 3.7 + 2.51 * x / reynolds).log10() > 0.0 {
            high = x;
        } else {
            low = x;
        }
    }
    let x = 0.5 * low + 0.5 * high;
    (LENGTH / DIAMETER) * RHO * (q / area()) * (q.abs() / area()) / (2.0 * x * x)
}

fn midpoint(q0: f64, dt: f64, drive: f64, mu: f64, roughness: f64, series: bool) -> f64 {
    let sign = q0.signum();
    let inertia = RHO * LENGTH / area() * if series { 2.0 } else { 1.0 };
    let linear = if series { LINEAR } else { 0.0 };
    let residual = |q| {
        2.0 * inertia * (q - q0.abs()) / dt + pressure(q, mu, roughness) + linear * q - sign * drive
    };
    let (mut low, mut high) = (0.5 * q0.abs(), 1.5 * q0.abs());
    assert!(residual(low) < 0.0 && residual(high) > 0.0);
    for _ in 0..96 {
        let q = 0.5 * low + 0.5 * high;
        if residual(q) > 0.0 {
            high = q;
        } else {
            low = q;
        }
    }
    sign * (0.5 * low + 0.5 * high)
}

fn propose(
    cx: &Cx<'_>,
    owner: &InertialGraph,
    graph: &LossGraph,
    average: f64,
    drive: f64,
    dt: f64,
) -> Result<InertialGraphTrial, InertialGraphError> {
    let mut drives = vec![0.0; graph.branches().len()];
    drives[0] = drive;
    owner.propose(
        cx,
        graph,
        &endpoints(graph),
        &drives,
        &vec![average; drives.len()],
        owner.accepted_time_s() + dt,
        controls(),
    )
}

#[test]
fn turbulent_midpoint_closes_original_loss_junction_impulse_and_work_in_both_directions() {
    with_cx(|cx| {
        for roughness in [0.0, 0.005] {
            let graph = network(0.001, roughness, true);
            for sign in [-1.0, 1.0] {
                let q0 = sign * 2e-4;
                let initial = state(cx, &graph, q0);
                let dt = 0.02;
                let drive = sign * 500.0;
                let average = midpoint(q0, dt, drive, 0.001, roughness, true);
                let (next, report) = propose(cx, &initial, &graph, average, drive, dt)
                    .unwrap()
                    .accept(cx, 1e-8, 1e-10)
                    .unwrap();
                for q in &report.average_volume_flows_m3_s {
                    near(*q, average, 2e-17);
                }
                for q in next.volume_flows_m3_s() {
                    near(q, 2.0 * average - q0, 2e-17);
                }
                near(
                    report.dissipative_pressure_drops_pa[0],
                    pressure(average, 0.001, roughness),
                    2e-9,
                );
                near(
                    report.dissipative_pressure_drops_pa[1],
                    LINEAR * average,
                    2e-9,
                );
                let loss = dt * average * (pressure(average, 0.001, roughness) + LINEAR * average);
                near(report.hydraulic_loss_j, loss, 2e-12);
                near(
                    next.kinetic_energy_j() - initial.kinetic_energy_j() + loss,
                    dt * drive * average,
                    2e-12,
                );
                near(report.node_outflows_m3_s[1], 0.0, 1e-17);
                near(
                    report.average_pressures_pa[1],
                    LINEAR * average + (RHO * LENGTH / area()) * (2.0 * average - 2.0 * q0) / dt,
                    2e-8,
                );
                assert!(report.hydraulic_loss_j > 0.0);
                assert_eq!(initial.accepted_time_s(), 0.0);
            }
        }
    });
}

#[test]
fn turbulent_native_decay_is_second_order_and_preserves_the_original_steady_limit() {
    with_cx(|cx| {
        let (mu, roughness, q0, duration) = (0.001, 0.001, 2e-4, 0.1);
        let graph = network(mu, roughness, false);
        let initial = state(cx, &graph, q0);
        let equilibrium = pressure(q0, mu, roughness);
        let (steady, _) = propose(cx, &initial, &graph, q0, equilibrium, 0.01)
            .unwrap()
            .accept(cx, 1e-8, 1e-10)
            .unwrap();
        near(steady.volume_flows_m3_s()[0], q0, 2e-17);

        // Independent continuous scalar momentum reference, classical RK4.
        let rate = |q| -pressure(q, mu, roughness) / (RHO * LENGTH / area());
        let h = duration / 1024.0;
        let mut reference = q0;
        for _ in 0..1024 {
            let k1 = rate(reference);
            let k2 = rate(reference + 0.5 * h * k1);
            let k3 = rate(reference + 0.5 * h * k2);
            let k4 = rate(reference + h * k3);
            reference += h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
        }
        let mut errors = Vec::new();
        for steps in [8, 16, 32] {
            let dt = duration / steps as f64;
            let mut current = initial.clone();
            for _ in 0..steps {
                let q = current.volume_flows_m3_s()[0];
                let average = midpoint(q, dt, 0.0, mu, roughness, false);
                let (next, report) = propose(cx, &current, &graph, average, 0.0, dt)
                    .unwrap()
                    .accept(cx, 1e-8, 1e-10)
                    .unwrap();
                assert!(next.kinetic_energy_j() < current.kinetic_energy_j());
                near(
                    current.kinetic_energy_j() - next.kinetic_energy_j(),
                    report.hydraulic_loss_j,
                    2e-12,
                );
                current = next;
            }
            errors.push((current.volume_flows_m3_s()[0] - reference).abs());
        }
        for pair in errors.windows(2) {
            let ratio = pair[0] / pair[1];
            assert!(
                (3.8..4.2).contains(&ratio),
                "midpoint decay errors {errors:?}"
            );
        }
    });
}

#[test]
fn refreshed_material_viscosity_changes_drag_without_changing_the_native_trajectory_on_recheck() {
    with_cx(|cx| {
        let cold = network(0.001, 0.001, false);
        let hot = network(0.0005, 0.001, false);
        let q0 = 2e-4;
        let initial = state(cx, &cold, q0);
        let dt = 0.02;
        let cold_average = midpoint(q0, dt, 0.0, 0.001, 0.001, false);
        let hot_average = midpoint(q0, dt, 0.0, 0.0005, 0.001, false);
        assert!(hot_average > cold_average);
        let cold_trial = propose(cx, &initial, &cold, cold_average, 0.0, dt).unwrap();
        let original_impulses = cold_trial.report().momentum.x.clone();
        let rechecked = cold_trial
            .recheck(cx, &hot, &endpoints(&hot), &[0.0], controls())
            .unwrap();
        assert_eq!(rechecked.report().momentum.x, original_impulses);
        assert!(matches!(
            rechecked.accept(cx, 1e-8, 1e-10),
            Err(InertialGraphError::MomentumResidual { .. })
        ));
        let (next, report) = propose(cx, &initial, &hot, hot_average, 0.0, dt)
            .unwrap()
            .accept(cx, 1e-8, 1e-10)
            .unwrap();
        near(next.volume_flows_m3_s()[0], 2.0 * hot_average - q0, 2e-17);
        near(
            report.dissipative_pressure_drops_pa[0],
            pressure(hot_average, 0.0005, 0.001),
            2e-9,
        );
        let changed_geometry = network(0.0005, 0.002, false);
        assert!(matches!(
            propose(cx, &initial, &changed_geometry, hot_average, 0.0, dt),
            Err(InertialGraphError::Graph(GraphError::InvalidBranch {
                branch: 0,
                ..
            }))
        ));
        assert_eq!(initial.accepted_time_s(), 0.0);
        near(initial.volume_flows_m3_s()[0], q0, 2e-17);
    });
}

#[test]
fn turbulent_support_and_transition_refusals_leave_the_accepted_state_intact() {
    with_cx(|cx| {
        let graph = network(0.001, 0.001, false);
        let geometry = [InertialBranch {
            length_m: LENGTH,
            flow_area_m2: area(),
        }];
        for q in [0.0, 1e-5, 0.01] {
            assert!(matches!(
                InertialGraph::try_new(
                    cx,
                    &graph,
                    Density::new(RHO),
                    &geometry,
                    &endpoints(&graph),
                    &[q],
                    0.0,
                    16
                ),
                Err(InertialGraphError::Graph(GraphError::Resistance {
                    branch: 0,
                    source: ResistanceError::Reynolds { .. }
                }))
            ));
        }
        assert!(matches!(
            InertialGraph::try_new(
                cx,
                &graph,
                Density::new(RHO + 1.0),
                &geometry,
                &endpoints(&graph),
                &[2e-4],
                0.0,
                16
            ),
            Err(InertialGraphError::Invalid(_))
        ));
        for (q0, q1, crosses_zero) in [(6e-5, 2e-5, false), (4e-4, -2e-4, true)] {
            let initial = state(cx, &graph, q0);
            let before = initial.volume_flows_m3_s();
            let average = 0.5 * q0 + 0.5 * q1;
            let dt = 0.001;
            // The candidate and actual midpoint are supported; only the
            // proposed endpoint or complete sign interval is inadmissible.
            pipe(0.001, 0.001)
                .evaluate(VolumetricFlowRate::new(average))
                .unwrap();
            let drive = (RHO * LENGTH / area()) * (q1 - q0) / dt + pressure(average, 0.001, 0.001);
            let error = propose(cx, &initial, &graph, average, drive, dt)
                .err()
                .unwrap();
            if crosses_zero {
                assert!(matches!(
                    error,
                    InertialGraphError::Graph(GraphError::InvalidBranch { branch: 0, .. })
                ));
            } else {
                assert!(matches!(
                    error,
                    InertialGraphError::Graph(GraphError::Resistance {
                        branch: 0,
                        source: ResistanceError::Reynolds { .. }
                    })
                ));
            }
            assert_eq!(initial.volume_flows_m3_s(), before);
            assert_eq!(initial.accepted_time_s(), 0.0);
        }
    });
}
