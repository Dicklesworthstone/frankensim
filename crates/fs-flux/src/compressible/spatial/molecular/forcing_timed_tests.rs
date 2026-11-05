//! Complete source-driven SSP/RKL splits with replayable physical bath times.

use super::*;
use crate::compressible::spatial::{ForcingWorkspace3d, ReservoirFaceState};

const OPEN_AXIS: AxisBoundaries = AxisBoundaries::Ends {
    lower: BoundaryCondition::Transmissive,
    upper: BoundaryCondition::Transmissive,
};
const OPEN: Boundaries3d = Boundaries3d {
    x: OPEN_AXIS,
    y: OPEN_AXIS,
    z: OPEN_AXIS,
};

fn bath(time: f64, side: BoundarySide) -> ReservoirFaceState {
    ReservoirFaceState {
        density: 1.05 + 0.03 * time,
        velocity: [0.13 + 0.2 * time, 0.01 * side.sign(), -0.02 + 0.05 * time],
        pressure_pa: 2.02 + 0.1 * (12.0 * time).sin(),
    }
}

fn profiled() -> Euler3d {
    let mut flow = field(
        gas(0.0),
        CartesianGrid::new(1, 1, 1, 1.0, 1.0, 1.0).unwrap(),
        OPEN,
        |_| [1.0, 0.1, -0.02, 0.03, 2.0],
    );
    let bytes = flow
        .required_reservoir_profile_workspace_bytes(&BoundarySide::ALL)
        .unwrap();
    flow.prepare_reservoir_profiles(
        &BoundarySide::ALL,
        bytes,
        |side, _, _| Ok(bath(0.0, side)),
        || false,
    )
    .unwrap();
    flow
}

fn make(conductivity: f64) -> (MolecularEuler3d<ConstantTransport>, ForcingWorkspace3d) {
    let flow = profiled();
    let source = law(flow.gas(), 0.02, 0.003, conductivity);
    let flow = MolecularEuler3d::new(flow, source, || false)
        .unwrap()
        .with_reservoir_exchange();
    let work = flow.flow().forcing_workspace(40).unwrap();
    (flow, work)
}

fn cap() -> ForceControl3d {
    ForceControl3d { maximum_mach: 10.0 }
}

#[allow(clippy::too_many_arguments)]
fn forced(
    flow: &mut MolecularEuler3d<ConstantTransport>,
    stabilized: bool,
    duration: f64,
    force: BodyForce3d<'_>,
    cap: ForceControl3d,
    work: &mut ForcingWorkspace3d,
    budget: usize,
    sample: impl FnMut(f64, BoundarySide, usize, [f64; 3]) -> Result<ReservoirFaceState, Error>,
    cancelled: impl FnMut() -> bool,
) -> Result<MolecularStepReport3d, MolecularError> {
    if stabilized {
        flow.step_stabilized_with_forcing_and_reservoir_profiles(
            Time::new(duration),
            MolecularStepControl::default(),
            MolecularStabilization::default(),
            force,
            cap,
            work,
            budget,
            sample,
            cancelled,
        )
        .map(|r| r.transport)
    } else {
        flow.step_with_forcing_and_reservoir_profiles(
            Time::new(duration),
            MolecularStepControl::default(),
            force,
            cap,
            work,
            budget,
            sample,
            cancelled,
        )
    }
}

#[allow(clippy::too_many_arguments)]
fn molecular_half(
    flow: &mut Euler3d,
    stabilized: bool,
    conductivity: f64,
    start: f64,
    duration: f64,
    times: &mut Vec<f64>,
) -> (MolecularReport3d, [Vec<Totals>; 6], usize) {
    let mut op =
        MolecularTransport3d::new(flow, law(flow.gas(), 0.02, 0.003, conductivity), || false)
            .unwrap()
            .with_reservoir_exchange();
    let sample = |time, side, _, _| {
        if side == BoundarySide::XMin {
            times.push(time);
        }
        Ok(bath(time, side))
    };
    let (report, stages) = if stabilized {
        let report = op
            .advance_stabilized_with_reservoir_profiles(
                flow,
                Time::new(start),
                Time::new(duration),
                MolecularControl::default(),
                MolecularStabilization::default(),
                10000,
                sample,
                || false,
            )
            .unwrap();
        (report.transport, report.maximum_stages_used)
    } else {
        (
            op.advance_with_reservoir_profiles(
                flow,
                Time::new(start),
                Time::new(duration),
                MolecularControl::default(),
                10000,
                sample,
                || false,
            )
            .unwrap(),
            2,
        )
    };
    (
        report,
        BoundarySide::ALL.map(|side| op.boundary_impulses(side).to_vec()),
        stages,
    )
}

#[test]
fn g1_forced_timed_ssp_rkl_matches_original_phases_and_separate_source_work() {
    let force = BodyForce3d::UniformForceDensity([0.6, -0.2, 0.1]);
    for stabilized in [false, true] {
        let conductivity = 20.0;
        let (mut flow, mut work) = make(conductivity);
        let mut times = Vec::new();
        let report = forced(
            &mut flow,
            stabilized,
            0.04,
            force,
            cap(),
            &mut work,
            10000,
            |time, side, _, _| {
                if side == BoundarySide::XMin {
                    times.push(time);
                }
                Ok(bath(time, side))
            },
            || false,
        )
        .unwrap();
        assert_eq!(report.coupling_retries, 0);
        assert_eq!(report.elapsed.value(), 0.04);
        let mut reference = profiled();
        let mut expected_times = vec![0.0]; // Complete-split exterior admission.
        let (first, first_loads, first_stages) = molecular_half(
            &mut reference,
            stabilized,
            conductivity,
            0.0,
            0.02,
            &mut expected_times,
        );
        let before_kicks = reference.states()[0];
        let mut reference_work = reference.forcing_workspace(40).unwrap();
        let central = reference
            .step_with_forcing_and_reservoir_profiles(
                report.elapsed,
                StepControl::default(),
                force,
                cap(),
                &mut reference_work,
                12,
                |time, side, _, _| {
                    if side == BoundarySide::XMin {
                        expected_times.push(time);
                    }
                    Ok(bath(time, side))
                },
                || false,
            )
            .unwrap();
        let central_loads =
            BoundarySide::ALL.map(|side| reference.boundary_impulses(side).to_vec());
        let after_kicks = reference.states()[0];
        let (second, second_loads, second_stages) = molecular_half(
            &mut reference,
            stabilized,
            conductivity,
            0.02,
            0.02,
            &mut expected_times,
        );
        if stabilized {
            assert!(first_stages > 2 && second_stages > 2);
        }
        assert_eq!(report.first, first);
        assert_eq!(report.advection, central);
        assert_eq!(report.second, second);
        assert_eq!(flow.flow().states(), reference.states());
        assert_eq!(times, expected_times);
        assert_eq!(times.last(), Some(&report.elapsed.value()));
        assert_eq!(report.source_change, central.source_change);
        let half_impulse = [0.6, -0.2, 0.1].map(|f| 0.02 * f);
        let expected_work: f64 = (0..3)
            .map(|a| {
                half_impulse[a] * (before_kicks.momentum[a] + 0.5 * half_impulse[a])
                    / before_kicks.density
                    + half_impulse[a] * (after_kicks.momentum[a] - 0.5 * half_impulse[a])
                        / after_kicks.density
            })
            .sum();
        close(report.source_change.work, expected_work, 3e-13);
        for a in 0..3 {
            close(
                report.source_change.impulse[a],
                2.0 * half_impulse[a],
                3e-13,
            );
        }
        assert!(report.source_change.work > 0.0);
        assert!(report.first.boundary_change.energy.abs() > 0.0);
        assert!(report.second.boundary_change.energy.abs() > 0.0);
        for (i, side) in BoundarySide::ALL.into_iter().enumerate() {
            let actual = flow.boundary_impulses(side)[0].array();
            let first = first_loads[i][0].array();
            let central = central_loads[i][0].array();
            let second = second_loads[i][0].array();
            for d in 0..5 {
                close(actual[d], first[d] + central[d] + second[d], 3e-13);
            }
            assert_eq!(flow.flow().reservoir_face(side, 0), Some(bath(0.04, side)));
        }
        for residual in report.conservation_residual.array() {
            close(residual, 0.0, 3e-12);
        }
    }
}

#[test]
fn g0_zero_forced_timed_ssp_rkl_preserves_complete_original_bits() {
    for stabilized in [false, true] {
        let (mut flow, mut work) = make(0.3);
        let (mut reference, _) = make(0.3);
        for _ in 0..2 {
            let expected = if stabilized {
                reference
                    .step_stabilized_with_reservoir_profiles(
                        Time::new(0.01),
                        MolecularStepControl::default(),
                        MolecularStabilization::default(),
                        1000,
                        |time, side, _, _| Ok(bath(time, side)),
                        || false,
                    )
                    .unwrap()
                    .transport
            } else {
                reference
                    .step_with_reservoir_profiles(
                        Time::new(0.01),
                        MolecularStepControl::default(),
                        1000,
                        |time, side, _, _| Ok(bath(time, side)),
                        || false,
                    )
                    .unwrap()
            };
            let actual = forced(
                &mut flow,
                stabilized,
                0.01,
                BodyForce3d::UniformAcceleration([0.0; 3]),
                cap(),
                &mut work,
                1000,
                |time, side, _, _| Ok(bath(time, side)),
                || false,
            )
            .unwrap();
            assert_eq!(actual, expected);
            assert_eq!(flow.flow().states(), reference.flow().states());
            assert_eq!(flow.flow().time(), reference.flow().time());
            assert_eq!(flow.flow().last_step(), reference.flow().last_step());
            for side in BoundarySide::ALL {
                assert_eq!(
                    flow.boundary_impulses(side),
                    reference.boundary_impulses(side)
                );
                assert_eq!(
                    flow.flow().reservoir_face(side, 0),
                    reference.flow().reservoir_face(side, 0)
                );
            }
        }
    }
}

#[test]
fn g2_forced_timed_whole_split_retry_resamples_both_halves_without_budget_refund() {
    for stabilized in [false, true] {
        let (mut flow, mut work) = make(0.3);
        let force = BodyForce3d::UniformAcceleration([80.0, 0.0, 0.0]);
        let cap = ForceControl3d { maximum_mach: 0.25 };
        let mut times = Vec::new();
        let mut calls = 0;
        let report = forced(
            &mut flow,
            stabilized,
            0.06,
            force,
            cap,
            &mut work,
            10000,
            |time, side, _, _| {
                calls += 1;
                if side == BoundarySide::XMin {
                    times.push(time);
                }
                Ok(bath(time, side))
            },
            || false,
        )
        .unwrap();
        assert!(report.coupling_retries > 0);
        assert!(times.windows(2).any(|t| t[1] < t[0]));
        assert_eq!(times.last(), Some(&report.elapsed.value()));
        assert_eq!(report.first.elapsed.value() * 2.0, report.elapsed.value());
        assert_eq!(report.second.elapsed, report.first.elapsed);
        let (mut direct, mut direct_work) = make(0.3);
        let short = forced(
            &mut direct,
            stabilized,
            report.elapsed.value(),
            force,
            cap,
            &mut direct_work,
            10000,
            |time, side, _, _| Ok(bath(time, side)),
            || false,
        )
        .unwrap();
        assert_eq!(short.coupling_retries, 0);
        assert_eq!(short.first, report.first);
        assert_eq!(short.advection, report.advection);
        assert_eq!(short.second, report.second);
        assert_eq!(flow.flow().states(), direct.flow().states());
        assert_eq!(short.source_change, report.source_change);
        for side in BoundarySide::ALL {
            assert_eq!(flow.boundary_impulses(side), direct.boundary_impulses(side));
        }
        let (mut refused, mut refused_work) = make(0.3);
        let original = refused.flow().states().to_vec();
        let mut used = 0;
        assert_eq!(
            forced(
                &mut refused,
                stabilized,
                0.06,
                force,
                cap,
                &mut refused_work,
                calls - 6,
                |time, side, _, _| {
                    used += 1;
                    Ok(bath(time, side))
                },
                || false,
            ),
            Err(MolecularError::Flow(Error::Invalid(
                "timed reservoir profile face evaluation budget exhausted"
            )))
        );
        assert_eq!(used, calls - 6);
        assert_eq!(refused.flow().states(), original);
        assert_eq!(refused.flow().time().value(), 0.0);
        assert!(refused.last_step().is_none());
        for side in BoundarySide::ALL {
            assert_eq!(
                refused.flow().reservoir_face(side, 0),
                Some(bath(0.0, side))
            );
        }
    }
}

#[test]
fn g4_forced_timed_second_half_refusals_preserve_all_accepted_data_and_held_recovery() {
    let force = BodyForce3d::UniformForceDensity([0.6, -0.2, 0.1]);
    for stabilized in [false, true] {
        let accepted = || {
            let (mut flow, mut work) = make(0.3);
            forced(
                &mut flow,
                stabilized,
                0.01,
                force,
                cap(),
                &mut work,
                1000,
                |time, side, _, _| Ok(bath(time, side)),
                || false,
            )
            .unwrap();
            (flow, work)
        };
        let (mut probe, mut probe_work) = accepted();
        let mut total_polls = 0;
        let mut total_calls = 0;
        forced(
            &mut probe,
            stabilized,
            0.01,
            force,
            cap(),
            &mut probe_work,
            1000,
            |time, side, _, _| {
                total_calls += 1;
                Ok(bath(time, side))
            },
            || {
                total_polls += 1;
                false
            },
        )
        .unwrap();
        for failure in 0..3 {
            let (mut flow, mut work) = accepted();
            let (mut reference, mut reference_work) = accepted();
            let original = flow.flow().states().to_vec();
            let previous = flow.last_step().copied();
            let euler_report = flow.flow().last_step().copied();
            let loads = BoundarySide::ALL.map(|side| flow.boundary_impulses(side).to_vec());
            let euler_loads =
                BoundarySide::ALL.map(|side| flow.flow().boundary_impulses(side).to_vec());
            let mut calls = 0;
            let mut polls = 0;
            let result = forced(
                &mut flow,
                stabilized,
                0.01,
                force,
                cap(),
                &mut work,
                if failure == 0 { total_calls - 6 } else { 1000 },
                |time, side, _, _| {
                    calls += 1;
                    if failure == 1 && calls == total_calls - 4 {
                        Err(Error::Invalid("late forced molecular trace unavailable"))
                    } else {
                        Ok(bath(time, side))
                    }
                },
                || {
                    polls += 1;
                    failure == 2 && polls == total_polls
                },
            );
            assert_eq!(
                result,
                Err(match failure {
                    0 => MolecularError::Flow(Error::Invalid(
                        "timed reservoir profile face evaluation budget exhausted"
                    )),
                    1 => MolecularError::Flow(Error::Invalid(
                        "late forced molecular trace unavailable"
                    )),
                    _ => MolecularError::Cancelled,
                })
            );
            assert_eq!(flow.flow().states(), original);
            assert_eq!(flow.flow().time().value(), 0.01);
            assert_eq!(flow.last_step().copied(), previous);
            assert_eq!(flow.flow().last_step().copied(), euler_report);
            for (i, side) in BoundarySide::ALL.into_iter().enumerate() {
                assert_eq!(flow.boundary_impulses(side), loads[i]);
                assert_eq!(flow.flow().boundary_impulses(side), euler_loads[i]);
                assert_eq!(flow.flow().reservoir_face(side, 0), Some(bath(0.01, side)));
            }
            let held = |flow: &mut MolecularEuler3d<ConstantTransport>,
                        work: &mut ForcingWorkspace3d| {
                if stabilized {
                    flow.step_stabilized_with_forcing(
                        Time::new(0.01),
                        MolecularStepControl::default(),
                        MolecularStabilization::default(),
                        force,
                        cap(),
                        work,
                        || false,
                    )
                    .unwrap()
                    .transport
                } else {
                    flow.step_with_forcing(
                        Time::new(0.01),
                        MolecularStepControl::default(),
                        force,
                        cap(),
                        work,
                        || false,
                    )
                    .unwrap()
                }
            };
            assert_eq!(
                held(&mut flow, &mut work),
                held(&mut reference, &mut reference_work)
            );
            assert_eq!(flow.flow().states(), reference.flow().states());
            for side in BoundarySide::ALL {
                assert_eq!(
                    flow.boundary_impulses(side),
                    reference.boundary_impulses(side)
                );
            }
        }
    }
}
