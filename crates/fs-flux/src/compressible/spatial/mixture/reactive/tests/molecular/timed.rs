//! Time-resolved chemistry, species/energy exchange and complete-trial rollback.

use super::*;
use crate::compressible::spatial::mixture::molecular::{
    DissipationError, InertialWall3d, MolecularWall3d, WallInertiaControl3d, WallThermalCondition,
};
use crate::compressible::spatial::mixture::{Error, ReservoirFaceState};
use fs_qty::Mass;

const SIDE: BoundarySide = BoundarySide::XMin;

fn prescribed(time: f64, face: usize, densities: &mut [f64]) -> ReservoirFaceState {
    let a = 0.65 + 0.02 * face as f64 + 2.0 * time;
    densities.copy_from_slice(&[a, 1.0 - a]);
    ReservoirFaceState {
        velocity: [0.15 + time, -0.02 + 0.5 * time, 0.03],
        pressure_pa: 1.05 + 3.0 * time,
    }
}

fn owner(heat: f64, rate: f64) -> ReactingEuler3d {
    let (gas, source) = mechanism(heat, rate, 40.0);
    let mut flow = field(
        gas,
        CartesianGrid::new(2, 2, 2, 0.2, 0.25, 0.3).unwrap(),
        Boundaries3d {
            x: AxisBoundaries::Ends {
                lower: BoundaryCondition::Transmissive,
                upper: BoundaryCondition::Transmissive,
            },
            y: WALLS,
            z: AxisBoundaries::Periodic,
        },
        |x| {
            (
                [0.8 + 0.01 * x[2], 0.2 - 0.01 * x[2]],
                [0.12, 0.0, 0.03],
                1.0,
            )
        },
    );
    let bytes = flow
        .required_reservoir_profile_workspace_bytes(&[SIDE])
        .unwrap();
    flow.prepare_reservoir_profiles(
        &[SIDE],
        bytes,
        |_, face, _, densities| Ok(prescribed(0.0, face, densities)),
        || false,
    )
    .unwrap();
    coupled(flow, source)
}

fn profiles(flow: &ReactingEuler3d) -> Vec<[f64; 6]> {
    (0..flow.flow().grid().face_count(SIDE))
        .map(|face| {
            let q = flow.flow().reservoir_face(SIDE, face).unwrap();
            [
                q.partial_densities[0],
                q.partial_densities[1],
                q.velocity[0],
                q.velocity[1],
                q.velocity[2],
                q.pressure_pa,
            ]
        })
        .collect()
}

#[derive(Debug, PartialEq)]
struct Published {
    states: Vec<f64>,
    time: Time,
    report: Option<ReactiveStepReport3d>,
    transport: Option<StepReport>,
    profiles: Vec<[f64; 6]>,
    boundary: Vec<[f64; 6]>,
    molecular: Option<MolecularReport>,
    molecular_boundary: Vec<[f64; 6]>,
    walls: [MolecularWall3d; 6],
}

fn published(flow: &ReactingEuler3d, op: &MolecularTransport3d<ConstantTransport>) -> Published {
    Published {
        states: flow.flow().states().to_vec(),
        time: flow.flow().time(),
        report: flow.last_step().copied(),
        transport: flow.flow().last_step().copied(),
        profiles: profiles(flow),
        boundary: loads(flow),
        molecular: op.last_step().copied(),
        molecular_boundary: molecular_loads(op),
        walls: BoundarySide::ALL.map(|side| op.wall(side)),
    }
}

#[allow(clippy::too_many_arguments)]
fn timed_split(
    flow: &mut ReactingEuler3d,
    op: &mut MolecularTransport3d<ConstantTransport>,
    stabilized: bool,
    dt: f64,
    maximum_faces: usize,
    sample: impl FnMut(
        f64,
        BoundarySide,
        usize,
        [f64; 3],
        &mut [f64],
    ) -> Result<ReservoirFaceState, Error>,
    cancelled: impl FnMut() -> bool,
) -> Result<ReactiveStepReport3d, ReactiveError> {
    let bytes = if stabilized {
        flow.required_stabilized_workspace_bytes(op).unwrap()
    } else {
        budget(flow, op)
    };
    if stabilized {
        flow.step_with_stabilized_dissipation_and_reservoir_profiles(
            op,
            Time::new(dt),
            ReactiveStepControl::default(),
            StabilizedDissipationControl::default(),
            bytes,
            maximum_faces,
            sample,
            cancelled,
        )
        .map(|report| report.coupled)
    } else {
        flow.step_with_dissipation_and_reservoir_profiles(
            op,
            Time::new(dt),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            bytes,
            maximum_faces,
            sample,
            cancelled,
        )
    }
}

fn warmup(flow: &mut ReactingEuler3d, op: &mut MolecularTransport3d<ConstantTransport>) {
    let bytes = budget(flow, op);
    let report = flow
        .step_with_dissipation(
            op,
            Time::new(0.0001),
            ReactiveStepControl::default(),
            DissipationControl::default(),
            bytes,
            || false,
        )
        .unwrap();
    balance(&report);
}

#[test]
fn g4_timed_reacting_constant_history_matches_zero_rate_and_original_forced_chemistry() {
    for rate in [0.0, 8.0] {
        for forced in [false, true] {
            let mut held = owner(5.0, rate);
            let mut timed = owner(5.0, rate);
            let duration = Time::new(0.0004);
            let control = ReactiveStepControl::default();
            let force = BodyForce3d::UniformAcceleration([0.3, -0.1, 0.2]);
            let cap = ForceControl3d { maximum_mach: 5.0 };
            let mut held_work = held.flow().forcing_workspace(usize::MAX).unwrap();
            let mut timed_work = timed.flow().forcing_workspace(usize::MAX).unwrap();
            let expected = if forced {
                held.step_with_forcing(duration, control, force, cap, &mut held_work, || false)
                    .unwrap()
            } else {
                held.step(duration, control, || false).unwrap()
            };
            let sample =
                |_, _, face, _, densities: &mut [f64]| Ok(prescribed(0.0, face, densities));
            let actual = if forced {
                timed
                    .step_with_forcing_and_reservoir_profiles(
                        duration,
                        control,
                        force,
                        cap,
                        &mut timed_work,
                        1000,
                        sample,
                        || false,
                    )
                    .unwrap()
            } else {
                timed
                    .step_with_reservoir_profiles(duration, control, 1000, sample, || false)
                    .unwrap()
            };
            balance(&actual);
            assert_eq!(actual, expected);
            assert_eq!(timed.flow().states(), held.flow().states());
            assert_eq!(profiles(&timed), profiles(&held));
            assert_eq!(loads(&timed), loads(&held));
            if rate == 0.0 {
                assert_eq!(actual.chemistry.species_change_kg, [0.0; MAX_SPECIES]);
                assert_eq!(actual.chemistry.referenced_energy_change_j, 0.0);
            } else {
                assert!(actual.chemistry.species_change_kg[0] < 0.0);
                assert!(actual.chemistry.referenced_energy_change_j > 0.0);
            }
        }
    }
}

#[test]
fn g4_timed_reacting_constant_history_matches_original_ssp_rkl_and_force_splits() {
    for stabilized in [false, true] {
        for forced in [false, true] {
            let mut held = owner(5.0, 8.0);
            let mut timed = owner(5.0, 8.0);
            let mut held_op = operator(held.flow(), 20.0);
            let mut timed_op = operator(timed.flow(), 20.0);
            let mut held_work = held.flow().forcing_workspace(usize::MAX).unwrap();
            let mut timed_work = timed.flow().forcing_workspace(usize::MAX).unwrap();
            let bytes = if stabilized {
                held.required_stabilized_workspace_bytes(&held_op).unwrap()
            } else {
                budget(&held, &held_op)
            } + if forced {
                held_work.workspace_bytes()
            } else {
                0
            };
            let duration = Time::new(0.002);
            let control = ReactiveStepControl::default();
            let force = BodyForce3d::UniformForceDensity([0.3, -0.1, 0.2]);
            let cap = ForceControl3d { maximum_mach: 5.0 };
            let sample =
                |_, _, face, _, densities: &mut [f64]| Ok(prescribed(0.0, face, densities));
            let (expected, actual) = match (stabilized, forced) {
                (false, false) => (
                    held.step_with_dissipation(
                        &mut held_op,
                        duration,
                        control,
                        DissipationControl::default(),
                        bytes,
                        || false,
                    )
                    .unwrap(),
                    timed
                        .step_with_dissipation_and_reservoir_profiles(
                            &mut timed_op,
                            duration,
                            control,
                            DissipationControl::default(),
                            bytes,
                            10000,
                            sample,
                            || false,
                        )
                        .unwrap(),
                ),
                (false, true) => (
                    held.step_with_dissipation_and_forcing(
                        &mut held_op,
                        duration,
                        control,
                        DissipationControl::default(),
                        bytes,
                        force,
                        cap,
                        &mut held_work,
                        || false,
                    )
                    .unwrap(),
                    timed
                        .step_with_dissipation_forcing_and_reservoir_profiles(
                            &mut timed_op,
                            duration,
                            control,
                            DissipationControl::default(),
                            bytes,
                            force,
                            cap,
                            &mut timed_work,
                            10000,
                            sample,
                            || false,
                        )
                        .unwrap(),
                ),
                (true, false) => {
                    let expected = held
                        .step_with_stabilized_dissipation(
                            &mut held_op,
                            duration,
                            control,
                            StabilizedDissipationControl::default(),
                            bytes,
                            || false,
                        )
                        .unwrap();
                    let actual = timed
                        .step_with_stabilized_dissipation_and_reservoir_profiles(
                            &mut timed_op,
                            duration,
                            control,
                            StabilizedDissipationControl::default(),
                            bytes,
                            10000,
                            sample,
                            || false,
                        )
                        .unwrap();
                    assert_eq!(actual, expected);
                    assert!(actual.maximum_stages_used > 2);
                    (expected.coupled, actual.coupled)
                }
                (true, true) => {
                    let expected = held
                        .step_with_stabilized_dissipation_and_forcing(
                            &mut held_op,
                            duration,
                            control,
                            StabilizedDissipationControl::default(),
                            bytes,
                            force,
                            cap,
                            &mut held_work,
                            || false,
                        )
                        .unwrap();
                    let actual = timed
                        .step_with_stabilized_dissipation_forcing_and_reservoir_profiles(
                            &mut timed_op,
                            duration,
                            control,
                            StabilizedDissipationControl::default(),
                            bytes,
                            force,
                            cap,
                            &mut timed_work,
                            10000,
                            sample,
                            || false,
                        )
                        .unwrap();
                    assert_eq!(actual, expected);
                    assert!(actual.maximum_stages_used > 2);
                    (expected.coupled, actual.coupled)
                }
            };
            balance(&actual);
            assert_eq!(actual, expected);
            assert_eq!(published(&timed, &timed_op), published(&held, &held_op));
        }
    }
}

#[test]
fn g4_timed_reacting_species_inlet_preserves_chemical_formation_and_body_energy_balances() {
    let mut flow = owner(5.0, 8.0);
    let mut op = operator(flow.flow(), 20.0);
    let mut force_work = flow.flow().forcing_workspace(usize::MAX).unwrap();
    let bytes =
        flow.required_stabilized_workspace_bytes(&op).unwrap() + force_work.workspace_bytes();
    let report = flow
        .step_with_stabilized_dissipation_forcing_and_reservoir_profiles(
            &mut op,
            Time::new(0.002),
            ReactiveStepControl::default(),
            StabilizedDissipationControl::default(),
            bytes,
            BodyForce3d::UniformAcceleration([0.3, -0.1, 0.2]),
            ForceControl3d { maximum_mach: 5.0 },
            &mut force_work,
            10000,
            |time, _, face, _, densities| Ok(prescribed(time, face, densities)),
            || false,
        )
        .unwrap();
    let physical = report.coupled;
    balance(&physical);
    assert_eq!(physical.elapsed.value(), 0.002);
    assert_eq!(physical.chemistry.cell_solves, 2 * flow.flow().cell_count());
    assert!(physical.chemistry.species_change_kg[0] < 0.0);
    assert!(physical.chemistry.referenced_energy_change_j > 0.0);
    assert!(physical.source_change.work != 0.0);
    assert!(physical.boundary_change.energy.abs() > 1e-8);
    assert!(report.maximum_stages_used > 2);
    for k in 0..2 {
        close(
            physical.after.species_masses[k] - physical.before.species_masses[k],
            physical.boundary_change.species_masses[k] + physical.chemistry.species_change_kg[k],
            2e-11,
        );
    }
    for face in 0..flow.flow().grid().face_count(SIDE) {
        let mut densities = [0.0; 2];
        let endpoint = prescribed(flow.flow().time().value(), face, &mut densities);
        let actual = flow.flow().reservoir_face(SIDE, face).unwrap();
        assert_eq!(actual.partial_densities, densities);
        assert_eq!(actual.velocity, endpoint.velocity);
        assert_eq!(actual.pressure_pa, endpoint.pressure_pa);
    }
    assert_eq!(op.last_step(), Some(&physical.dissipation.unwrap().second));
    assert_eq!(flow.last_step(), Some(&physical));
}

#[test]
fn g4_timed_reacting_heat_release_replays_whole_trials_and_counts_all_source_calls() {
    let mut flow = owner(40.0, 500.0);
    let mut times = Vec::new();
    let mut calls = 0usize;
    let actual = flow
        .step_with_reservoir_profiles(
            Time::new(0.2),
            ReactiveStepControl::default(),
            10000,
            |time, _, face, _, densities| {
                calls += 1;
                if face == 0 {
                    times.push(time);
                }
                Ok(prescribed(time, face, densities))
            },
            || false,
        )
        .unwrap();
    balance(&actual);
    assert!(actual.coupling_retries > 0);
    assert!(times.iter().any(|time| *time > actual.elapsed.value()));
    assert!(
        times
            .windows(2)
            .any(|times| times[0] > 0.0 && times[1] == 0.0)
    );
    let mut replay = owner(40.0, 500.0);
    let expected = replay
        .step_with_reservoir_profiles(
            actual.elapsed,
            ReactiveStepControl::default(),
            10000,
            |time, _, face, _, densities| Ok(prescribed(time, face, densities)),
            || false,
        )
        .unwrap();
    assert_eq!(expected.coupling_retries, 0);
    assert_eq!(flow.flow().states(), replay.flow().states());
    assert_eq!(actual.transport, expected.transport);
    assert_eq!(actual.chemistry, expected.chemistry);
    assert_eq!(actual.boundary_change, expected.boundary_change);
    assert_eq!(profiles(&flow), profiles(&replay));

    let mut refused = owner(40.0, 500.0);
    let op = operator(refused.flow(), 0.02);
    let before = published(&refused, &op);
    let mut refused_calls = 0usize;
    let error = refused
        .step_with_reservoir_profiles(
            Time::new(0.2),
            ReactiveStepControl::default(),
            calls - 1,
            |time, _, face, _, densities| {
                refused_calls += 1;
                Ok(prescribed(time, face, densities))
            },
            || false,
        )
        .unwrap_err();
    assert!(!matches!(error, ReactiveError::RetryLimit { .. }));
    assert_eq!(refused_calls, calls - 1);
    assert_eq!(published(&refused, &op), before);
}

#[test]
fn g4_timed_reacting_split_samples_physical_stage_times_and_exact_rounded_endpoint() {
    for stabilized in [false, true] {
        let mut flow = owner(5.0, 8.0);
        let mut op = operator(flow.flow(), 0.02);
        let start = 6e-6;
        let end: f64 = 7e-6;
        let dt = end - start;
        flow.accepted.time = start;
        flow.trial.time = start;
        assert!((start + 0.5 * dt) + 0.5 * dt > end);
        let mut times = Vec::new();
        let report = timed_split(
            &mut flow,
            &mut op,
            stabilized,
            dt,
            10000,
            |time, _, face, _, densities| {
                if time < start || time > end {
                    return Err(Error::Invalid(
                        "history queried outside its physical interval",
                    ));
                }
                if face == 0 {
                    times.push(time);
                }
                let mut boundary = prescribed(time, face, densities);
                boundary.velocity[2] = time;
                Ok(boundary)
            },
            || false,
        )
        .unwrap();
        let offsets: &[f64] = if stabilized {
            &[0.0, 0.0, 1.0 / 6.0, 0.5, 0.0, 1.0, 0.5, 2.0 / 3.0, 1.0]
        } else {
            &[0.0, 0.0, 0.5, 0.0, 1.0, 0.5, 1.0]
        };
        assert_eq!(times.len(), offsets.len());
        for (time, offset) in times.iter().zip(offsets) {
            close(*time, start + offset * dt, 1e-15);
        }
        assert_eq!(flow.flow().time().value().to_bits(), end.to_bits());
        assert_eq!(times.last().unwrap().to_bits(), end.to_bits());
        assert_eq!(report.elapsed.value().to_bits(), dt.to_bits());
        for face in 0..flow.flow().grid().face_count(SIDE) {
            assert_eq!(
                flow.flow().reservoir_face(SIDE, face).unwrap().velocity[2].to_bits(),
                end.to_bits()
            );
        }
        balance(&report);
    }
}

#[test]
fn g4_timed_reacting_late_provider_admission_and_work_failures_preserve_held_recovery() {
    for stabilized in [false, true] {
        for failure in 0..3 {
            let mut flow = owner(5.0, 8.0);
            let mut op = operator(flow.flow(), 0.02);
            warmup(&mut flow, &mut op);
            let before = published(&flow, &op);
            let mut calls = 0usize;
            let cap = if failure == 0 { 24 } else { 10000 };
            let error = timed_split(
                &mut flow,
                &mut op,
                stabilized,
                0.001,
                cap,
                |time, _, face, _, densities| {
                    calls += 1;
                    let boundary = prescribed(time, face, densities);
                    if calls == 25 {
                        if failure == 1 {
                            return Err(Error::Invalid("prescribed chemical inflow refused"));
                        }
                        if failure == 2 {
                            densities[0] = -1.0;
                        }
                    }
                    Ok(boundary)
                },
                || false,
            )
            .unwrap_err();
            assert!(!matches!(error, ReactiveError::RetryLimit { .. }));
            assert_eq!(calls, if failure == 0 { 24 } else { 25 });
            assert_eq!(published(&flow, &op), before);

            let mut reference = owner(5.0, 8.0);
            let mut reference_op = operator(reference.flow(), 0.02);
            warmup(&mut reference, &mut reference_op);
            warmup(&mut reference, &mut reference_op);
            // An abandoned timed trial has privately published intermediate
            // profiles; the ordinary next call must restore accepted history.
            warmup(&mut flow, &mut op);
            assert_eq!(published(&flow, &op), published(&reference, &reference_op));
        }
    }
}

#[test]
fn g4_timed_reacting_final_cancel_is_atomic_and_finite_wall_inertia_refuses_before_sampling() {
    for stabilized in [false, true] {
        let mut reference = owner(5.0, 8.0);
        let mut reference_op = operator(reference.flow(), 0.02);
        warmup(&mut reference, &mut reference_op);
        let mut total_polls = 0usize;
        timed_split(
            &mut reference,
            &mut reference_op,
            stabilized,
            0.001,
            10000,
            |time, _, face, _, densities| Ok(prescribed(time, face, densities)),
            || {
                total_polls += 1;
                false
            },
        )
        .unwrap();
        let mut flow = owner(5.0, 8.0);
        let mut op = operator(flow.flow(), 0.02);
        warmup(&mut flow, &mut op);
        let before = published(&flow, &op);
        let mut polls = 0usize;
        let error = timed_split(
            &mut flow,
            &mut op,
            stabilized,
            0.001,
            10000,
            |time, _, face, _, densities| Ok(prescribed(time, face, densities)),
            || {
                polls += 1;
                polls == total_polls
            },
        )
        .unwrap_err();
        assert!(matches!(error, ReactiveError::Cancelled));
        assert_eq!(polls, total_polls);
        assert_eq!(published(&flow, &op), before);

        let mut held = owner(5.0, 8.0);
        let mut held_op = operator(held.flow(), 0.02);
        warmup(&mut held, &mut held_op);
        warmup(&mut held, &mut held_op);
        warmup(&mut flow, &mut op);
        assert_eq!(published(&flow, &op), published(&held, &held_op));

        let mut finite = owner(5.0, 8.0);
        let ordinary = operator(finite.flow(), 0.02);
        let source =
            SpatialTransport::new(finite.flow().gas().clone(), ordinary.law().clone(), || {
                false
            })
            .unwrap();
        let bytes = MolecularTransport3d::required_workspace_bytes(finite.flow(), &source).unwrap()
            + WallInertiaControl3d::solver_workspace_bytes();
        let mut finite_op = MolecularTransport3d::new(finite.flow(), source, bytes, || false)
            .unwrap()
            .with_wall(
                BoundarySide::YMax,
                MolecularWall3d::NoSlip {
                    velocity: [0.1, 0.0, 0.0],
                    thermal: WallThermalCondition::Adiabatic,
                },
            )
            .unwrap()
            .with_wall_inertia(
                &[InertialWall3d {
                    side: BoundarySide::YMax,
                    mass: Mass::new(0.5),
                }],
                WallInertiaControl3d::default(),
            )
            .unwrap();
        let before = published(&finite, &finite_op);
        let mut calls = 0usize;
        let error = timed_split(
            &mut finite,
            &mut finite_op,
            stabilized,
            0.001,
            10000,
            |time, _, face, _, densities| {
                calls += 1;
                Ok(prescribed(time, face, densities))
            },
            || false,
        )
        .unwrap_err();
        assert_eq!(
            error,
            ReactiveError::Dissipation(DissipationError::Invalid(
                "time-resolved molecular reservoirs require prescribed walls without finite inertia",
            ))
        );
        assert_eq!(calls, 0);
        assert_eq!(published(&finite, &finite_op), before);
    }
}
