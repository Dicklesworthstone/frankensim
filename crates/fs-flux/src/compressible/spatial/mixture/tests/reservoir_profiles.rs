use super::*;

fn open_axis(axis: usize) -> Boundaries3d {
    let mut axes = [AxisBoundaries::Periodic; 3];
    axes[axis] = OPEN;
    Boundaries3d {
        x: axes[0],
        y: axes[1],
        z: axes[2],
    }
}

fn input(
    densities: &mut [f64],
    first: f64,
    velocity: [f64; 3],
    pressure_pa: f64,
) -> Result<ReservoirFaceState, Error> {
    densities.copy_from_slice(&[first, 1.0 - first, 0.0]);
    Ok(ReservoirFaceState {
        velocity,
        pressure_pa,
    })
}

fn snapshot(flow: &MixtureEuler3d) -> Vec<f64> {
    let mut rows = Vec::new();
    for side in BoundarySide::ALL {
        for face in 0..flow.grid().face_count(side) {
            if let Some(p) = flow.reservoir_face(side, face) {
                rows.extend_from_slice(p.partial_densities);
                rows.extend_from_slice(&p.velocity);
                rows.push(p.pressure_pa);
            }
        }
    }
    rows
}

#[test]
fn g0_mixture_face_profiles_preserve_uniform_native_results_and_canonical_order() {
    let gas = standard(0.2);
    let grid = CartesianGrid::new(2, 3, 2, 0.1, 0.2, 0.3).unwrap();
    let topology = Boundaries3d {
        x: OPEN,
        y: OPEN,
        z: OPEN,
    };
    let densities = [0.7, 1.0 - 0.7, 0.0];
    let velocity = [0.15, -0.1, 0.08];
    for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
        let initial = field(gas.clone(), grid, topology, |_| (densities, velocity, 2.0));
        let inputs = BoundarySide::ALL.map(|side| ReservoirBoundary {
            side,
            partial_densities: &densities,
            velocity,
            pressure_pa: 2.0,
        });
        let bytes =
            MixtureEuler3d::required_workspace_bytes_with_reservoirs(&gas, &grid, 6).unwrap();
        let mut uniform = MixtureEuler3d::new_with_reservoirs(
            gas.clone(),
            grid,
            initial.states().to_vec(),
            topology,
            &inputs,
            bytes,
            || false,
        )
        .unwrap();
        let mut profile = initial;
        let base = profile.workspace_bytes();
        let sides = [
            BoundarySide::ZMax,
            BoundarySide::XMin,
            BoundarySide::YMax,
            BoundarySide::XMax,
            BoundarySide::ZMin,
            BoundarySide::YMin,
        ];
        let bytes = profile
            .required_reservoir_profile_workspace_bytes(&sides)
            .unwrap();
        let mut visited = Vec::new();
        profile
            .prepare_reservoir_profiles(
                &sides,
                bytes,
                |side, face, position, densities| {
                    assert_eq!(position, grid.face_center(side, face).unwrap());
                    assert!(densities.iter().all(|v| v.is_nan()));
                    visited.push((side, face));
                    input(densities, 0.7, velocity, 2.0)
                },
                || false,
            )
            .unwrap();
        let expected: Vec<_> = BoundarySide::ALL
            .into_iter()
            .flat_map(|side| (0..grid.face_count(side)).map(move |face| (side, face)))
            .collect();
        assert_eq!(visited, expected);
        assert_eq!(profile.prepared_reservoir_profile_sides(), [true; 6]);
        assert_eq!(
            profile.workspace_bytes(),
            base + profile.reservoir_profile_workspace_bytes()
        );
        for side in BoundarySide::ALL {
            assert!(profile.has_reservoir_profile(side));
            assert!(profile.reservoir(side).is_none());
            for face in 0..grid.face_count(side) {
                assert_eq!(
                    profile.reservoir_face(side, face),
                    uniform.reservoir_face(side, face)
                );
            }
            assert!(
                profile
                    .reservoir_face(side, grid.face_count(side))
                    .is_none()
            );
        }
        for _ in 0..3 {
            let controls = StepControl {
                reconstruction,
                ..StepControl::default()
            };
            let a = profile.step(Time::new(0.001), controls, || false).unwrap();
            let b = uniform.step(Time::new(0.001), controls, || false).unwrap();
            assert_eq!(profile.states(), uniform.states());
            assert_eq!(a, b);
            balance(a);
        }
    }
}

#[test]
fn g1_each_profile_face_carries_its_own_species_shear_and_enthalpy_on_every_side() {
    let grid = CartesianGrid::new(3, 3, 3, 0.1, 0.1, 0.1).unwrap();
    for side in BoundarySide::ALL {
        for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
            let axis = side.axis();
            let [a, b] = tangents(axis);
            let mut velocity = [0.0; 3];
            velocity[axis] = -side.sign() * 4.0;
            let mut flow = field(standard(0.2), grid, open_axis(axis), |_| {
                ([1.0, 0.0, 0.0], velocity, 1.0)
            });
            let bytes = flow
                .required_reservoir_profile_workspace_bytes(&[side])
                .unwrap();
            flow.prepare_reservoir_profiles(
                &[side],
                bytes,
                |actual, face, position, densities| {
                    assert_eq!(actual, side);
                    assert_eq!(position, grid.face_center(side, face).unwrap());
                    let mut u = velocity;
                    u[a] = 0.1 * face as f64;
                    u[b] = -0.05 * face as f64;
                    input(densities, 0.1 + 0.08 * face as f64, u, 1.0)
                },
                || false,
            )
            .unwrap();
            let mut expected = Vec::new();
            for face in 0..grid.face_count(side) {
                let exterior = flow.reservoir_face(side, face).unwrap();
                let p = flow.reservoirs.primitive_face(side, face).unwrap();
                let r = exterior.partial_densities[0] + 0.5 * exterior.partial_densities[1];
                let temperature = 1.0 / r;
                let cp = 3.5 * exterior.partial_densities[0]
                    + 0.5 * exterior.partial_densities[1] * (3.5 + 0.2 * temperature);
                close(p.temperature, temperature, 3e-13);
                close(p.sound_speed.powi(2), cp / (cp - r), 3e-13);
                let mut q = [0.0; 7];
                state_into(
                    flow.gas(),
                    exterior.partial_densities,
                    exterior.velocity,
                    exterior.pressure_pa,
                    &mut q,
                )
                .unwrap();
                let normal = exterior.velocity[axis];
                // Supersonic incoming flow makes the physical Euler flux equal
                // to the prescribed exterior state, independently of HLLC details.
                expected.push((
                    std::array::from_fn::<_, 3, _>(|s| side.sign() * q[s] * normal),
                    std::array::from_fn::<_, 3, _>(|d| {
                        side.sign() * (q[d + 3] * normal + if d == axis { 1.0 } else { 0.0 })
                    }),
                    side.sign() * (q[6] + 1.0) * normal,
                ));
            }
            let report = flow
                .step(
                    Time::new(1e-4),
                    StepControl {
                        reconstruction,
                        ..StepControl::default()
                    },
                    || false,
                )
                .unwrap();
            balance(report);
            let scale = report.elapsed.value() * grid.face_area(side);
            for (load, (mass, momentum, energy)) in flow.boundary_impulses(side).zip(expected) {
                for d in 0..3 {
                    close(load.species_masses[d], scale * mass[d], 3e-12);
                    close(load.momentum[d], scale * momentum[d], 3e-12);
                }
                close(load.energy, scale * energy, 3e-11);
                assert_eq!(load.species_masses[2], 0.0);
            }
            assert!(report.boundary_change.species_masses[1] > 0.0);
        }
    }
}

#[test]
fn g4_profile_budget_admission_and_last_poll_replacement_are_atomic() {
    let grid = CartesianGrid::new(2, 3, 2, 0.1, 0.2, 0.3).unwrap();
    let topology = Boundaries3d {
        x: OPEN,
        y: OPEN,
        z: OPEN,
    };
    let make = |boundaries| {
        field(standard(0.0), grid, boundaries, |_| {
            ([1.0, 0.0, 0.0], [0.0; 3], 1.0)
        })
    };
    let sides = [BoundarySide::YMax, BoundarySide::XMin];
    let mut flow = make(topology);
    let bytes = flow
        .required_reservoir_profile_workspace_bytes(&sides)
        .unwrap();
    let initial_bytes = flow.workspace_bytes();
    let mut calls = 0usize;
    assert!(
        flow.prepare_reservoir_profiles(
            &sides,
            bytes - 1,
            |_, _, _, densities| {
                calls += 1;
                input(densities, 0.3, [0.0; 3], 1.0)
            },
            || false
        )
        .is_err()
    );
    assert_eq!(calls, 0);
    assert_eq!(flow.workspace_bytes(), initial_bytes);
    assert_eq!(flow.prepared_reservoir_profile_sides(), [false; 6]);
    assert!(
        flow.required_reservoir_profile_workspace_bytes(&[])
            .is_err()
    );
    assert!(
        flow.required_reservoir_profile_workspace_bytes(&[sides[0], sides[0]])
            .is_err()
    );
    for axis in [
        AxisBoundaries::Periodic,
        WALLS,
        AxisBoundaries::Ends {
            lower: BoundaryCondition::PressureOutlet { pressure_pa: 1.0 },
            upper: BoundaryCondition::Transmissive,
        },
    ] {
        let other = make(Boundaries3d {
            x: axis,
            ..topology
        });
        assert!(
            other
                .required_reservoir_profile_workspace_bytes(&[BoundarySide::XMin])
                .is_err()
        );
    }
    let uniform = ReservoirBoundary {
        side: BoundarySide::XMin,
        partial_densities: &[1.0, 0.0, 0.0],
        velocity: [0.0; 3],
        pressure_pa: 1.0,
    };
    let required =
        MixtureEuler3d::required_workspace_bytes_with_reservoirs(flow.gas(), &grid, 1).unwrap();
    let other = MixtureEuler3d::new_with_reservoirs(
        flow.gas().clone(),
        grid,
        flow.states().to_vec(),
        topology,
        &[uniform],
        required,
        || false,
    )
    .unwrap();
    assert!(
        other
            .required_reservoir_profile_workspace_bytes(&[uniform.side])
            .is_err()
    );
    flow.prepare_reservoir_profiles(
        &sides,
        bytes,
        |_, _, _, densities| input(densities, 0.3, [0.0; 3], 1.0),
        || false,
    )
    .unwrap();
    let previous = snapshot(&flow);
    let workspace = flow.workspace_bytes();
    let final_side = BoundarySide::YMax;
    let final_face = grid.face_count(final_side) - 1;
    for missing_species in [false, true] {
        let result = flow.update_reservoir_profiles(
            |side, face, _, densities| {
                let value = input(densities, 0.6, [0.01, 0.02, -0.01], 1.0)?;
                if side == final_side && face == final_face {
                    densities[2] = if missing_species { f64::NAN } else { -0.1 };
                }
                Ok(value)
            },
            || false,
        );
        assert!(result.is_err());
        assert_eq!(snapshot(&flow), previous);
        assert_eq!(flow.workspace_bytes(), workspace);
    }
    let mut total_polls = 0usize;
    flow.update_reservoir_profiles(
        |_, _, _, densities| input(densities, 0.3, [0.0; 3], 1.0),
        || {
            total_polls += 1;
            false
        },
    )
    .unwrap();
    let cells = flow.states().to_vec();
    let mut polls = 0usize;
    let result = flow.update_reservoir_profiles(
        |_, _, _, densities| input(densities, 0.6, [0.01, 0.02, -0.01], 1.0),
        || {
            polls += 1;
            polls == total_polls
        },
    );
    assert!(matches!(result, Err(Error::Cancelled)));
    assert_eq!(snapshot(&flow), previous);
    assert_eq!(flow.states(), cells);
    assert_eq!(flow.time().value(), 0.0);
    assert!(flow.last_step().is_none());
    flow.update_reservoir_profiles(
        |_, _, _, densities| input(densities, 0.6, [0.01, 0.02, -0.01], 1.0),
        || false,
    )
    .unwrap();
    assert_ne!(snapshot(&flow), previous);
    assert_eq!(flow.workspace_bytes(), workspace);
    balance(
        flow.step(Time::new(0.001), StepControl::default(), || false)
            .unwrap(),
    );
}

#[test]
fn g1_one_hot_face_bounds_the_initial_native_three_axis_interval() {
    let grid = CartesianGrid::new(2, 2, 2, 0.1, 0.1, 0.1).unwrap();
    let side = BoundarySide::XMin;
    let mut flow = field(standard(0.0), grid, open_axis(0), |_| {
        ([1.0, 0.0, 0.0], [0.0; 3], 0.02)
    });
    let bytes = flow
        .required_reservoir_profile_workspace_bytes(&[side])
        .unwrap();
    flow.prepare_reservoir_profiles(
        &[side],
        bytes,
        |_, face, _, densities| {
            input(
                densities,
                0.5,
                [0.0; 3],
                if face == 3 { 32.0 } else { 0.02 },
            )
        },
        || false,
    )
    .unwrap();
    let controls = StepControl {
        cfl: 0.2,
        maximum_halvings: 0,
        reconstruction: Reconstruction::FirstOrder,
    };
    let report = flow.step(Time::new(1.0), controls, || false).unwrap();
    let bound = controls.cfl * 0.1 / (3.0 * (1.4_f64 * 32.0).sqrt());
    assert!(report.elapsed.value() <= bound * (1.0 + 2e-14));
    assert_eq!(report.halvings, 0);
    balance(report);
    flow.update_reservoir_profiles(
        |_, _, _, densities| input(densities, 0.5, [0.0; 3], 0.02),
        || false,
    )
    .unwrap();
    close(
        flow.reservoirs.maximum_rate(&grid),
        3.0 * (1.4_f64 * 0.02).sqrt() / 0.1,
        3e-13,
    );
}
