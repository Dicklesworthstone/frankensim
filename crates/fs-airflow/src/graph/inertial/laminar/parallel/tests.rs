//! G1/G3/G4: physical flow splitting, branch heat and transactional continuation.

use super::*;
use crate::graph::inertial::laminar::{LaminarPipeSpec, tests::startup};
use crate::graph::tests::with_cx;
use fs_material::fluid_transport::FluidTransportState;
use fs_phs::ZkDiffusiveBasis;
use fs_qty::Length;

fn fluid(mu: f64) -> FluidTransportState {
    FluidTransportState::try_new(1000.0, 4200.0, 0.6, mu).unwrap()
}
fn pipe(
    cx: &Cx<'_>,
    radius: f64,
    mu: f64,
    maximum_reynolds: f64,
) -> LaminarPipe<FluidTransportState> {
    LaminarPipe::try_new(
        cx,
        fluid(mu),
        LaminarPipeSpec {
            length: Length::new(0.4),
            radius: Length::new(radius),
            maximum_reynolds,
            band_rad_s: [0.01, 100.0],
            maximum_relative_wall_error: 1e-4,
        },
        &ZkDiffusiveBasis::new(16).unwrap(),
    )
    .unwrap()
}
fn bundle(cx: &Cx<'_>) -> LaminarParallel<FluidTransportState> {
    LaminarParallel::try_new(
        cx,
        vec![
            pipe(cx, 0.001, 0.001, 1000.0),
            pipe(cx, 0.002, 0.001, 1000.0),
        ],
    )
    .unwrap()
}
fn near(actual: f64, expected: f64, scale: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * scale,
        "{actual:e} != {expected:e}; scale {scale:e}"
    );
}
fn unchanged(actual: &LaminarParallelState, expected: &LaminarParallelState) {
    assert_eq!(actual.accepted_time_s(), expected.accepted_time_s());
    for (actual, expected) in actual.branches().iter().zip(expected.branches()) {
        assert_eq!(actual.coordinates(), expected.coordinates());
        assert_eq!(actual.accepted_time_s(), expected.accepted_time_s());
    }
}
fn balance(report: &LaminarParallelStep, branches: &[LaminarPipeStep], pressure: f64) {
    let scale = report
        .native
        .supplied
        .abs()
        .max(report.native.dissipated)
        .max(report.initial_mean_kinetic_energy_j + report.initial_profile_memory_energy_j);
    near(report.native.balance_residual(), 0.0, scale, 3e-10);
    near(
        report.native.supplied,
        pressure * report.displaced_volume_m3,
        scale,
        3e-12,
    );
    near(
        branches.iter().map(|r| r.native.dissipated).sum(),
        report.native.dissipated,
        scale,
        1e-14,
    );
    assert!(branches.iter().all(|r| r.native.dissipated >= 0.0));
    assert!(report.native.solver_certificate.is_none());
}

#[test]
fn parallel_startup_and_reversal_match_each_radial_stokes_solution() {
    with_cx(|cx| {
        let bundle = bundle(cx);
        // Equal length and viscosity, radii 1:2: steady Poiseuille flows are
        // 1:16, whereas the initial inertial accelerations are 1:4.
        let pi = core::f64::consts::PI;
        let steady = [0.001_f64, 0.002].map(|r| 100.0 * pi * r.powi(4) / (8.0 * 0.001 * 0.4));
        for dt in [0.04, 0.2, 1.0] {
            let mut state = bundle.steady_state(Pressure::new(0.0), 0.0).unwrap();
            let mut prepared = bundle.prepare(cx, dt).unwrap();
            let (history, mobility) = prepared.average_flow_relation(&state).unwrap();
            assert_eq!(history, 0.0);
            let report = prepared
                .advance(cx, &mut state, [Pressure::new(100.0), Pressure::new(0.0)])
                .unwrap();
            for (index, tau) in [1.0, 4.0].into_iter().enumerate() {
                let (endpoint, average, energy) = startup(dt / tau);
                let branch = &prepared.branch_steps()[index];
                let plug = 0.5
                    * bundle.branches()[index].mean_inertance_pa_s2_m3()
                    * steady[index].powi(2);
                near(
                    branch.final_volume_flow_m3_s,
                    steady[index] * endpoint,
                    steady[index],
                    3e-9,
                );
                near(
                    branch.average_volume_flow_m3_s,
                    steady[index] * average,
                    steady[index],
                    3e-9,
                );
                near(
                    state.branches()[index].kinetic_energy_j(),
                    plug * energy,
                    plug,
                    3e-9,
                );
                let expected_heat = 100.0 * dt * steady[index] * average - plug * energy;
                near(
                    branch.native.dissipated,
                    expected_heat,
                    100.0 * dt * steady[index],
                    3e-9,
                );
            }
            near(
                report.average_volume_flow_m3_s,
                100.0 * mobility,
                steady.iter().sum(),
                3e-12,
            );
            near(
                report.final_volume_flow_m3_s,
                state.volume_flow_m3_s(),
                steady.iter().sum(),
                1e-14,
            );
            balance(&report, prepared.branch_steps(), 100.0);
            if dt == 0.04 {
                let ratio =
                    state.branches()[1].volume_flow_m3_s() / state.branches()[0].volume_flow_m3_s();
                assert!(
                    ratio > 4.0 && ratio < 8.0,
                    "startup must not use the steady split: {ratio}"
                );
            }
        }
        let mut reverse = bundle.steady_state(Pressure::new(100.0), 0.0).unwrap();
        let mut prepared = bundle.prepare(cx, 0.2).unwrap();
        let report = prepared
            .advance(
                cx,
                &mut reverse,
                [Pressure::new(-100.0), Pressure::new(0.0)],
            )
            .unwrap();
        for (index, tau) in [1.0, 4.0].into_iter().enumerate() {
            near(
                reverse.branches()[index].volume_flow_m3_s(),
                steady[index] * (1.0 - 2.0 * startup(0.2 / tau).0),
                steady[index],
                5e-9,
            );
        }
        assert!(reverse.branches()[0].volume_flow_m3_s() < 0.0);
        assert!(reverse.branches()[1].volume_flow_m3_s() > 0.0);
        balance(&report, prepared.branch_steps(), -100.0);
    });
}

#[test]
fn identical_parallel_branches_preserve_the_original_pipe_and_zero_drive_heat() {
    with_cx(|cx| {
        let original = pipe(cx, 0.001, 0.001, 1000.0);
        let bundle = LaminarParallel::try_new(
            cx,
            vec![
                pipe(cx, 0.001, 0.001, 1000.0),
                pipe(cx, 0.001, 0.001, 1000.0),
            ],
        )
        .unwrap();
        let mut state = bundle.steady_state(Pressure::new(0.0), 2.0).unwrap();
        let mut single = original
            .steady_state(VolumetricFlowRate::new(0.0), 2.0)
            .unwrap();
        let mut prepared = bundle.prepare(cx, 0.025).unwrap();
        let mut one = original.prepare(cx, 0.025).unwrap();
        for pressure in [100.0, 100.0, -50.0, 0.0] {
            let ports = [Pressure::new(1e9 + pressure), Pressure::new(1e9)];
            let expected = one.advance(cx, &mut single, ports).unwrap();
            let actual = prepared.advance(cx, &mut state, ports).unwrap();
            for branch in state.branches() {
                assert_eq!(branch.coordinates(), single.coordinates());
            }
            near(
                actual.average_volume_flow_m3_s,
                2.0 * expected.average_volume_flow_m3_s,
                1e-7,
                1e-14,
            );
            near(
                state.kinetic_energy_j(),
                2.0 * single.kinetic_energy_j(),
                single.kinetic_energy_j(),
                1e-14,
            );
            balance(&actual, prepared.branch_steps(), pressure);
            if pressure == 0.0 {
                assert_eq!(actual.native.supplied, 0.0);
                assert!(actual.native.delta_h < 0.0 && actual.native.dissipated > 0.0);
            }
        }
        for omega in [0.01, 1.0, 10.0, 100.0] {
            let expected = original.impedance_at(omega).unwrap().scale(0.5);
            let actual = bundle.impedance_at(omega).unwrap();
            near(actual.re, expected.re, expected.abs(), 1e-14);
            near(actual.im, expected.im, expected.abs(), 1e-14);
        }
        assert!(bundle.impedance_at(101.0).is_err());
    });
}

#[test]
fn total_mean_flow_pressure_preserves_opposing_branches_and_real_work() {
    with_cx(|cx| {
        let bundle = bundle(cx);
        let mut state = bundle.steady_state(Pressure::new(20.0), 0.0).unwrap();
        let mut prepared = bundle.prepare(cx, 0.2).unwrap();
        for target in [0.0, 4e-7, -2e-7] {
            let pressure = prepared
                .pressure_for_average_flow(&state, VolumetricFlowRate::new(target))
                .unwrap();
            let report = prepared
                .advance(cx, &mut state, [pressure, Pressure::new(0.0)])
                .unwrap();
            near(report.average_volume_flow_m3_s, target, 1e-6, 3e-12);
            balance(&report, prepared.branch_steps(), pressure.value());
            if target == 0.0 {
                let branch = prepared.branch_steps();
                assert!(branch[0].displaced_volume_m3 * branch[1].displaced_volume_m3 < 0.0);
                assert!(report.native.delta_h < 0.0);
                assert!(report.native.supplied.abs() < 1e-10 * report.native.dissipated);
            }
        }
    });
}

#[test]
fn viscosity_continuation_retains_original_nonsteady_profiles_and_branch_heat() {
    with_cx(|cx| {
        let original = bundle(cx);
        let mut state = original.steady_state(Pressure::new(0.0), 3.0).unwrap();
        original
            .prepare(cx, 0.05)
            .unwrap()
            .advance(cx, &mut state, [Pressure::new(100.0), Pressure::new(0.0)])
            .unwrap();
        let before = state.clone();
        let updated = original
            .with_sources(cx, vec![fluid(0.004), fluid(0.002)])
            .unwrap();
        unchanged(&state, &before);
        let independent = original
            .branches()
            .iter()
            .zip([0.004, 0.002])
            .map(|(p, mu)| p.with_source(cx, fluid(mu)).unwrap())
            .collect::<Vec<_>>();
        let mut singles = before.branches().to_vec();
        let mut single_steps = independent
            .iter()
            .map(|p| p.prepare(cx, 0.03).unwrap())
            .collect::<Vec<_>>();
        let mut prepared = updated.prepare(cx, 0.03).unwrap();
        for pressure in [-50.0, 0.0, 40.0] {
            let ports = [Pressure::new(pressure), Pressure::new(0.0)];
            let report = prepared.advance(cx, &mut state, ports).unwrap();
            for index in 0..2 {
                let single = single_steps[index]
                    .advance(cx, &mut singles[index], ports)
                    .unwrap();
                assert_eq!(
                    state.branches()[index].coordinates(),
                    singles[index].coordinates()
                );
                assert_eq!(prepared.branch_steps()[index].native, single.native);
            }
            balance(&report, prepared.branch_steps(), pressure);
        }
        let assembled = updated
            .state_from_branches(state.branches().to_vec())
            .unwrap();
        unchanged(&assembled, &state);
        assert!(updated.with_sources(cx, vec![fluid(0.003)]).is_err());
        let changed_density = FluidTransportState::try_new(1000.1, 4200.0, 0.6, 0.002).unwrap();
        assert!(
            updated
                .with_sources(cx, vec![fluid(0.003), changed_density])
                .is_err()
        );
    });
}

#[test]
fn late_branch_refusal_preserves_all_histories_and_last_successful_reports() {
    with_cx(|cx| {
        let bundle = LaminarParallel::try_new(
            cx,
            vec![pipe(cx, 0.001, 0.004, 1000.0), pipe(cx, 0.002, 0.001, 1.0)],
        )
        .unwrap();
        let mut state = bundle.steady_state(Pressure::new(0.0), 4.0).unwrap();
        let mut prepared = bundle.prepare(cx, 0.05).unwrap();
        assert!(prepared.branch_steps().is_empty());
        prepared
            .advance(cx, &mut state, [Pressure::new(0.1), Pressure::new(0.0)])
            .unwrap();
        let before = state.clone();
        let reports = prepared.branch_steps().to_vec();
        let gate = fs_exec::CancelGate::new();
        gate.request();
        let cancelled = Cx::new(&gate, cx.arena(), cx.stream_key(), cx.budget(), cx.mode());
        assert!(matches!(
            prepared.advance(
                &cancelled,
                &mut state,
                [Pressure::new(0.1), Pressure::new(0.0)]
            ),
            Err(LaminarPipeError::Execution(
                crate::graph::GraphError::Interrupted
            ))
        ));
        unchanged(&state, &before);
        for ports in [
            [Pressure::new(100.0), Pressure::new(0.0)],
            [Pressure::new(f64::NAN), Pressure::new(0.0)],
        ] {
            assert!(prepared.advance(cx, &mut state, ports).is_err());
            unchanged(&state, &before);
            for (actual, expected) in prepared.branch_steps().iter().zip(&reports) {
                assert_eq!(actual.native, expected.native);
                assert_eq!(
                    actual.final_volume_flow_m3_s,
                    expected.final_volume_flow_m3_s
                );
            }
        }
        // A failed second branch cannot contaminate the next valid trial.
        let mut expected = before.clone();
        let ports = [Pressure::new(-0.1), Pressure::new(0.0)];
        bundle
            .prepare(cx, 0.05)
            .unwrap()
            .advance(cx, &mut expected, ports)
            .unwrap();
        prepared.advance(cx, &mut state, ports).unwrap();
        unchanged(&state, &expected);
        let foreign = self::bundle(cx);
        let mut foreign_state = foreign.steady_state(Pressure::new(0.0), 4.0).unwrap();
        assert!(prepared.advance(cx, &mut foreign_state, ports).is_err());
        assert!(
            bundle
                .state_from_branches(foreign_state.branches().to_vec())
                .is_err()
        );
        let mut wrong_clock = state.branches().to_vec();
        wrong_clock[1].accepted_time_s += 1.0;
        assert!(bundle.state_from_branches(wrong_clock).is_err());
        assert!(LaminarParallel::<FluidTransportState>::try_new(cx, vec![]).is_err());
    });
}
