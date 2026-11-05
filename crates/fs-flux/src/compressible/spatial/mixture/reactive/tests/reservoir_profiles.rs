use super::*;
use crate::compressible::spatial::mixture::ReservoirFaceState;

#[test]
fn g4_reacting_profile_copy_update_and_last_poll_are_one_transaction() {
    let (gas, mechanism) = mechanism(0.0, 2.0, 0.0);
    let grid = CartesianGrid::new(2, 2, 2, 0.1, 0.1, 0.1).unwrap();
    let side = BoundarySide::XMin;
    let mut input = field(
        gas,
        grid,
        Boundaries3d {
            x: AxisBoundaries::Ends {
                lower: BoundaryCondition::Transmissive,
                upper: BoundaryCondition::Transmissive,
            },
            y: AxisBoundaries::Periodic,
            z: AxisBoundaries::Periodic,
        },
        |_| ([0.5, 0.5], [0.1, 0.0, 0.0], 1.0),
    );
    let profile_bytes = input
        .required_reservoir_profile_workspace_bytes(&[side])
        .unwrap();
    input
        .prepare_reservoir_profiles(
            &[side],
            profile_bytes,
            |_, face, _, densities| {
                let first = 0.3 + 0.1 * face as f64;
                densities.copy_from_slice(&[first, 1.0 - first]);
                Ok(ReservoirFaceState {
                    velocity: [0.1, 0.0, 0.0],
                    pressure_pa: 1.0,
                })
            },
            || false,
        )
        .unwrap();
    let mut reacting = coupled(input, mechanism);
    let bytes = reacting.workspace_bytes();
    assert_eq!(
        reacting.accepted.reservoir_profile_workspace_bytes(),
        profile_bytes
    );
    assert_eq!(
        reacting.trial.reservoir_profile_workspace_bytes(),
        profile_bytes
    );
    let mut calls = 0usize;
    let mut total_polls = 0usize;
    reacting
        .update_reservoir_profiles(
            |_, face, _, densities| {
                calls += 1;
                let first = 0.5 + 0.05 * face as f64;
                densities.copy_from_slice(&[first, 1.0 - first]);
                Ok(ReservoirFaceState {
                    velocity: [0.2, 0.0, 0.0],
                    pressure_pa: 1.1,
                })
            },
            || {
                total_polls += 1;
                false
            },
        )
        .unwrap();
    assert_eq!(calls, grid.face_count(side));
    let mut polls = 0usize;
    let result = reacting.update_reservoir_profiles(
        |_, _, _, densities| {
            densities.copy_from_slice(&[0.9, 0.1]);
            Ok(ReservoirFaceState {
                velocity: [0.4, 0.0, 0.0],
                pressure_pa: 1.5,
            })
        },
        || {
            polls += 1;
            polls == total_polls
        },
    );
    assert!(matches!(result, Err(ReactiveError::Cancelled)));
    for face in 0..grid.face_count(side) {
        let actual = reacting.accepted.reservoir_face(side, face).unwrap();
        assert_eq!(actual, reacting.trial.reservoir_face(side, face).unwrap());
        assert_eq!(actual.partial_densities[0], 0.5 + 0.05 * face as f64);
        assert_eq!(actual.velocity, [0.2, 0.0, 0.0]);
        assert_eq!(actual.pressure_pa, 1.1);
    }
    assert_eq!(reacting.flow().time().value(), 0.0);
    assert!(reacting.last_step().is_none());
    assert_eq!(reacting.workspace_bytes(), bytes);
    let report = reacting
        .step(Time::new(0.001), ReactiveStepControl::default(), || false)
        .unwrap();
    balance(&report);
    assert!(report.chemistry.cell_solves > 0);
    assert!(report.boundary_change.mass() > 0.0);
    for face in 0..grid.face_count(side) {
        assert_eq!(
            reacting.accepted.reservoir_face(side, face),
            reacting.trial.reservoir_face(side, face)
        );
    }
}
