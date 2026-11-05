use super::*;
use crate::compressible::spatial::ReservoirFaceState;

fn exterior(face: usize, revision: usize) -> ReservoirFaceState {
    ReservoirFaceState {
        density: 1.0,
        velocity: [0.02 * revision as f64, 0.01 * face as f64, -0.01],
        pressure_pa: 1.2 + 0.03 * face as f64 + 0.02 * revision as f64,
    }
}

#[test]
fn g1_g4_profiles_survive_molecular_owner_swaps_and_cancelled_control_updates() {
    let source = gas(0.4);
    let grid = CartesianGrid::new(3, 2, 2, 0.1, 0.1, 0.1).unwrap();
    let boundary = Boundaries3d {
        x: AxisBoundaries::Ends {
            lower: BoundaryCondition::Transmissive,
            upper: BoundaryCondition::Reflecting,
        },
        y: AxisBoundaries::Periodic,
        z: AxisBoundaries::Periodic,
    };
    let side = BoundarySide::XMin;
    let make = || {
        field(source.clone(), grid, boundary, |_| {
            [1.0, 0.0, 0.0, 0.0, 1.0]
        })
    };
    for prepare_before_composition in [false, true] {
        let mut reference = make();
        let required = reference
            .required_reservoir_profile_workspace_bytes(&[side])
            .unwrap();
        reference
            .prepare_reservoir_profiles(
                &[side],
                required,
                |_, face, _| Ok(exterior(face, 0)),
                || false,
            )
            .unwrap();
        let mut initial = make();
        if prepare_before_composition {
            initial
                .prepare_reservoir_profiles(
                    &[side],
                    required,
                    |_, face, _| Ok(exterior(face, 0)),
                    || false,
                )
                .unwrap();
        }
        let mut coupled =
            MolecularEuler3d::new(initial, law(&source, 0.0, 0.0, 0.0), || false).unwrap();
        if !prepare_before_composition {
            assert_eq!(
                coupled
                    .required_reservoir_profile_workspace_bytes(&[side])
                    .unwrap(),
                2 * required
            );
            assert!(
                coupled
                    .prepare_reservoir_profiles(
                        &[side],
                        2 * required - 1,
                        |_, _, _| unreachable!("combined budget refuses before sampling"),
                        || false
                    )
                    .is_err()
            );
            coupled
                .prepare_reservoir_profiles(
                    &[side],
                    2 * required,
                    |_, face, _| Ok(exterior(face, 0)),
                    || false,
                )
                .unwrap();
        }
        assert_eq!(coupled.reservoir_profile_workspace_bytes(), 2 * required);
        for revision in 1..=3 {
            reference
                .update_reservoir_profiles(|_, face, _| Ok(exterior(face, revision)), || false)
                .unwrap();
            coupled
                .update_reservoir_profiles(|_, face, _| Ok(exterior(face, revision)), || false)
                .unwrap();
            // Two successive complete steps exercise both ownership swaps,
            // including a step with no intervening boundary update.
            for _ in 0..2 {
                let expected = reference
                    .step(Time::new(0.001), StepControl::default(), || false)
                    .unwrap();
                let actual = coupled
                    .step(Time::new(0.001), MolecularStepControl::default(), || false)
                    .unwrap();
                assert_eq!(actual.advection, expected);
                assert_eq!(coupled.flow().states(), reference.states());
                assert_eq!(coupled.flow().time(), reference.time());
                for face in 0..grid.face_count(side) {
                    assert_eq!(
                        coupled.flow().reservoir_face(side, face),
                        Some(exterior(face, revision))
                    );
                }
                close(actual.conservation_residual.energy, 0.0, 2e-13);
            }
        }
        let cells = coupled.flow().states().to_vec();
        let time = coupled.flow().time();
        let report = coupled.last_step().copied();
        let impulses = coupled.boundary_impulses(side).to_vec();
        let stop = std::cell::Cell::new(false);
        assert!(
            coupled
                .update_reservoir_profiles(
                    |_, face, _| {
                        stop.set(face + 1 == grid.face_count(side));
                        Ok(exterior(face, 4))
                    },
                    || stop.get()
                )
                .is_err()
        );
        assert_eq!(coupled.flow().states(), cells);
        assert_eq!(coupled.flow().time(), time);
        assert_eq!(coupled.last_step().copied(), report);
        assert_eq!(coupled.boundary_impulses(side), impulses);
        for face in 0..grid.face_count(side) {
            assert_eq!(
                coupled.flow().reservoir_face(side, face),
                Some(exterior(face, 3))
            );
        }
    }
}
