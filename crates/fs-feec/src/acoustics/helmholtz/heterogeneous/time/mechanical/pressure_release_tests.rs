use super::*;
use crate::acoustics::helmholtz::HelmholtzFe;
use crate::acoustics::helmholtz::heterogeneous::time::tests::{MEDIUM, brick, p2, plane, uniform};
use crate::acoustics::helmholtz::pressure_release::EmptyFreePressure;
use fs_math::det;

fn restricted(prepared: &PreparedMechanicalTime<'_, '_, '_>, field: &[f64]) -> Vec<f64> {
    prepared
        .wall_dofs(0)
        .unwrap()
        .iter()
        .map(|&dof| field[dof])
        .collect()
}

fn assert_same_bits(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len());
    for (&a, &b) in actual.iter().zip(expected) {
        assert_eq!(a.to_bits(), b.to_bits());
    }
}

#[test]
fn pressure_release_wall_spring_balance_and_accelerating_fluid_have_exact_work() {
    let pin = WallPin {
        surface_density: 0.4,
        stiffness_per_area: 1.7,
        resistance: 0.3,
    };
    for height in [0.1, 0.3] {
        let lengths = [1.0, height, 0.2];
        let area = lengths[1] * lengths[2];
        let (complex, positions) = brick(2, lengths);
        let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
        let model = uniform(&geometry);
        let inlet = plane(&geometry, 0.0);
        let outlet = plane(&geometry, 1.0);
        let h = 0.02;
        let prepared = model
            .prepare_time_with_walls_and_pressure_release(
                h,
                &[],
                &[MechanicalBoundary {
                    faces: inlet.clone(),
                    wall: pin,
                    relaxation: Vec::new(),
                }],
                &outlet,
                MechanicalTimeLimits::default(),
                &mut || false,
            )
            .unwrap();
        for pressure_scale in [0.08, 8e-10] {
            // Exact continuum solution for L=1, including the original wall:
            // p=P(1-x), psi=t P(1-x), v_x=t P/rho, w=P/K, w_dot=0.
            // The imposed inlet velocity is t P/rho; the release outlet has
            // zero pressure work although its particle velocity is nonzero.
            let pressure = p2(&geometry, |point| pressure_scale * (1.0 - point[0]));
            let zero = vec![0.0; model.ndof()];
            let displacement = restricted(
                &prepared,
                &p2(&geometry, |_| pressure_scale / pin.stiffness_per_area),
            );
            let velocity = vec![0.0; displacement.len()];
            let mut state = prepared
                .state_with_walls(&zero, &pressure, &[displacement], &[velocity])
                .unwrap();
            let initial = prepared.energy(&state).unwrap().total_j();
            let compression = area * pressure_scale.powi(2)
                / (6.0 * MEDIUM.density_kg_m3 * MEDIUM.sound_speed_m_s.powi(2));
            let spring = area * pressure_scale.powi(2) / (2.0 * pin.stiffness_per_area);
            let mut total_work = 0.0;
            let mut total_heat = 0.0;
            for index in 1..=8 {
                let time = index as f64 * h;
                let midpoint_time = time - 0.5 * h;
                let drive = inlet
                    .iter()
                    .map(|&face| {
                        (
                            face,
                            [midpoint_time * pressure_scale / MEDIUM.density_kg_m3; 3],
                        )
                    })
                    .collect::<Vec<_>>();
                let next = prepared.step(&state, &drive, &mut || false).unwrap();
                let energy = prepared.energy(&next.state).unwrap();
                let kinetic = area * (time * pressure_scale).powi(2) / (2.0 * MEDIUM.density_kg_m3);
                for (actual, expected) in [
                    (energy.acoustic.kinetic_j, kinetic),
                    (energy.acoustic.compression_j, compression),
                    (energy.wall_elastic_j, spring),
                ] {
                    assert!(
                        (actual / expected - 1.0).abs() < 3e-10,
                        "actual={actual:e}, expected={expected:e}"
                    );
                }
                assert!(energy.wall_kinetic_j < 3e-10 * initial);
                assert_eq!(energy.wall_relaxation_j, 0.0);
                for (dof, &expected) in pressure.iter().enumerate() {
                    assert!(
                        (next.state.pressure_pa()[dof] - expected).abs() < 3e-10 * pressure_scale
                    );
                    assert!(
                        (next.state.pressure_potential_pa_s()[dof] - time * expected).abs()
                            < 3e-10 * time * pressure_scale
                    );
                    if prepared.condensed.pressure_release.constrained[dof] {
                        assert_eq!(next.state.pressure_pa()[dof], 0.0);
                        assert_eq!(next.state.pressure_potential_pa_s()[dof], 0.0);
                    }
                }
                let input =
                    h * area * pressure_scale.powi(2) * midpoint_time / MEDIUM.density_kg_m3;
                assert!((next.work.input_work_j / input - 1.0).abs() < 3e-10);
                assert!((next.work.acoustic_delta_energy_j / input - 1.0).abs() < 3e-10);
                assert!(next.work.wall_work_j.abs() < 3e-10 * initial);
                assert!(next.work.wall_dissipation_j < 3e-10 * initial);
                assert_eq!(next.work.boundary_dissipation_j, 0.0);
                total_work += next.work.input_work_j;
                total_heat += next.work.wall_dissipation_j;
                assert!(
                    (energy.total_j() + total_heat - initial - total_work).abs() < 3e-10 * initial
                );
                state = next.state;
            }
        }
    }
}

#[test]
fn pressure_release_adjacent_wall_retains_full_maxwell_memory_and_passive_pulse() {
    let lengths = [1.0, 0.2, 0.3];
    let area = lengths[0] * lengths[2];
    let (complex, positions) = brick(2, lengths);
    let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let wall_faces = geometry
        .boundary_faces()
        .into_iter()
        .filter(|face| face.iter().all(|&v| positions[v as usize][1] == 0.0))
        .collect::<Vec<_>>();
    let k = 2.5;
    let h = 0.005;
    let prepared = model
        .prepare_time_with_walls_and_pressure_release(
            h,
            &[],
            &[MechanicalBoundary {
                faces: wall_faces,
                wall: WallPin {
                    surface_density: 0.6,
                    stiffness_per_area: 3.0,
                    resistance: 0.0,
                },
                relaxation: vec![AcousticWallRelaxation {
                    stiffness_per_area: k,
                    relaxation_time_s: 0.04,
                }],
            }],
            &plane(&geometry, 1.0),
            MechanicalTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let amplitude = 0.01;
    let memory = restricted(
        &prepared,
        &p2(&geometry, |point| {
            det::sqrt(k) * amplitude * (1.0 + point[0]).powi(2)
        }),
    );
    let wall_dofs = prepared.wall_dofs(0).unwrap();
    assert!(wall_dofs.iter().enumerate().any(|(local, &original)| {
        prepared.condensed.pressure_release.constrained[original] && memory[local] != 0.0
    }));
    assert!(
        wall_dofs
            .iter()
            .enumerate()
            .any(|(local, &original)| { original >= positions.len() && memory[local] != 0.0 })
    );
    let zero = vec![0.0; model.ndof()];
    let wall_zero = vec![0.0; wall_dofs.len()];
    let mut state = prepared
        .state_with_wall_memory(
            &zero,
            &zero,
            &[wall_zero.clone()],
            &[wall_zero],
            &[vec![memory.clone()]],
        )
        .unwrap();
    // Integral_0^1 (1+x)^4 dx=31/5. This includes the original P2 edge
    // functions and wall coordinates on the neighboring fluid-release edge.
    let expected = 0.5 * k * area * amplitude * amplitude * 31.0 / 5.0;
    let initial = prepared.energy(&state).unwrap();
    assert!((initial.wall_relaxation_j / expected - 1.0).abs() < 3e-12);
    assert_eq!(initial.acoustic.total_j(), 0.0);
    let mut heat = 0.0;
    let mut exchange = 0.0;
    for _ in 0..12 {
        let next = prepared.step(&state, &[], &mut || false).unwrap();
        assert_eq!(next.work.input_work_j, 0.0);
        assert_eq!(next.work.boundary_dissipation_j, 0.0);
        assert_eq!(
            next.work.wall_dissipation_j,
            next.work.wall_relaxation_dissipation_j
        );
        assert!(next.work.wall_relaxation_dissipation_j > 0.0);
        assert!(next.work.acoustic_balance_residual_j().abs() < 3e-10 * expected);
        assert!(next.work.wall_balance_residual_j().abs() < 3e-10 * expected);
        for (dof, &fixed) in prepared
            .condensed
            .pressure_release
            .constrained
            .iter()
            .enumerate()
        {
            if fixed {
                assert_eq!(next.state.pressure_pa()[dof], 0.0);
                assert_eq!(next.state.pressure_potential_pa_s()[dof], 0.0);
            }
        }
        heat += next.work.wall_dissipation_j;
        exchange += next.work.wall_work_j;
        state = next.state;
    }
    let energy = prepared.energy(&state).unwrap();
    assert!(energy.acoustic.total_j() > 0.0);
    assert!(energy.wall_relaxation_j < expected);
    assert!(exchange < 0.0);
    assert!(
        state
            .wall_normal_velocity_m_s(0)
            .unwrap()
            .iter()
            .any(|v| v.abs() > 1e-8)
    );
    assert_ne!(
        state.wall_relaxation_memory(0, 0).unwrap(),
        memory.as_slice()
    );
    assert!((energy.total_j() + heat - expected).abs() < 3e-10 * expected);
}

#[test]
fn pressure_release_wall_admission_is_transactional_and_empty_faces_preserve_bits() {
    let (complex, positions) = brick(1, [1.0, 0.2, 0.3]);
    let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let release = plane(&geometry, 1.0);
    let boundary = MechanicalBoundary {
        faces: plane(&geometry, 0.0),
        wall: WallPin {
            surface_density: 0.4,
            stiffness_per_area: 3.0,
            resistance: 0.1,
        },
        relaxation: vec![AcousticWallRelaxation {
            stiffness_per_area: 2.0,
            relaxation_time_s: 0.1,
        }],
    };
    let free = geometry
        .pressure_release_constraints(&release, EmptyFreePressure::Refuse, &mut || false)
        .unwrap()
        .free_dofs
        .len();
    let limits = MechanicalTimeLimits {
        acoustic: HeterogeneousTimeLimits {
            max_dofs: free,
            ..HeterogeneousTimeLimits::default()
        },
        ..MechanicalTimeLimits::default()
    };
    let prepared = model
        .prepare_time_with_walls_and_pressure_release(
            0.01,
            &[],
            &[boundary.clone()],
            &release,
            limits,
            &mut || false,
        )
        .unwrap();
    let count = prepared.wall_dofs(0).unwrap().len();
    for bounded in [
        MechanicalTimeLimits {
            acoustic: HeterogeneousTimeLimits {
                max_dofs: free - 1,
                ..limits.acoustic
            },
            ..limits
        },
        MechanicalTimeLimits {
            max_wall_dofs: count - 1,
            ..limits
        },
        MechanicalTimeLimits {
            max_wall_relaxation_dofs: count - 1,
            ..limits
        },
    ] {
        assert!(
            model
                .prepare_time_with_walls_and_pressure_release(
                    0.01,
                    &[],
                    &[boundary.clone()],
                    &release,
                    bounded,
                    &mut || false,
                )
                .is_err()
        );
    }
    let mut overlap = boundary.clone();
    overlap.faces = release
        .iter()
        .map(|face| [face[2], face[0], face[1]])
        .collect();
    assert!(
        model
            .prepare_time_with_walls_and_pressure_release(
                0.01,
                &[],
                &[overlap],
                &release,
                limits,
                &mut || false,
            )
            .is_err()
    );
    assert!(
        model
            .prepare_time_with_walls_and_pressure_release(
                0.01,
                &[ResistiveBoundary {
                    faces: release.clone(),
                    admittance_m_pa_s: 0.0
                }],
                &[boundary.clone()],
                &release,
                limits,
                &mut || false,
            )
            .is_err()
    );
    let zero = vec![0.0; model.ndof()];
    let mut incompatible = zero.clone();
    let fixed_edge = prepared
        .condensed
        .pressure_release
        .constrained
        .iter()
        .enumerate()
        .find(|(dof, fixed)| *dof >= positions.len() && **fixed)
        .unwrap()
        .0;
    incompatible[fixed_edge] = f64::MIN_POSITIVE;
    assert!(prepared.state(&incompatible, &zero).is_err());
    assert!(prepared.state(&zero, &incompatible).is_err());
    let wall_zero = vec![vec![0.0; count]];
    let memories = vec![vec![restricted(
        &prepared,
        &p2(&geometry, |point| 0.01 * (1.0 + point[1]).powi(2)),
    )]];
    let state = prepared
        .state_with_wall_memory(&zero, &zero, &wall_zero, &wall_zero, &memories)
        .unwrap();
    let saved = state.clone();
    let forbidden = release
        .iter()
        .map(|&face| (face, [0.0; 3]))
        .collect::<Vec<_>>();
    assert!(prepared.step(&state, &forbidden, &mut || false).is_err());
    let drive = boundary
        .faces
        .iter()
        .map(|&face| (face, [0.01; 3]))
        .collect::<Vec<_>>();
    let mut polls = 0;
    assert!(matches!(
        prepared.step(&state, &drive, &mut || {
            polls += 1;
            polls == 20
        }),
        Err(HeterogeneousError::Cancelled)
    ));
    assert_eq!(polls, 20);
    assert_same_bits(state.pressure_pa(), saved.pressure_pa());
    assert_same_bits(
        state.pressure_potential_pa_s(),
        saved.pressure_potential_pa_s(),
    );
    assert_same_bits(
        state.wall_displacement_m(0).unwrap(),
        saved.wall_displacement_m(0).unwrap(),
    );
    assert_same_bits(
        state.wall_normal_velocity_m_s(0).unwrap(),
        saved.wall_normal_velocity_m_s(0).unwrap(),
    );
    assert_same_bits(
        state.wall_relaxation_memory(0, 0).unwrap(),
        saved.wall_relaxation_memory(0, 0).unwrap(),
    );
    prepared.step(&state, &drive, &mut || false).unwrap();

    let legacy = model
        .prepare_time_with_walls(
            0.01,
            &[],
            &[boundary.clone()],
            MechanicalTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let explicit = model
        .prepare_time_with_walls_and_pressure_release(
            0.01,
            &[],
            &[boundary],
            &[],
            MechanicalTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let old = legacy
        .state_with_wall_memory(&zero, &zero, &wall_zero, &wall_zero, &memories)
        .unwrap();
    let new = explicit
        .state_with_wall_memory(&zero, &zero, &wall_zero, &wall_zero, &memories)
        .unwrap();
    assert!(prepared.step(&old, &drive, &mut || false).is_err());
    let old = legacy.step(&old, &drive, &mut || false).unwrap();
    let new = explicit.step(&new, &drive, &mut || false).unwrap();
    assert_same_bits(new.state.pressure_pa(), old.state.pressure_pa());
    assert_same_bits(
        new.state.pressure_potential_pa_s(),
        old.state.pressure_potential_pa_s(),
    );
    assert_same_bits(
        new.state.wall_displacement_m(0).unwrap(),
        old.state.wall_displacement_m(0).unwrap(),
    );
    assert_same_bits(
        new.state.wall_normal_velocity_m_s(0).unwrap(),
        old.state.wall_normal_velocity_m_s(0).unwrap(),
    );
    assert_same_bits(
        new.state.wall_relaxation_memory(0, 0).unwrap(),
        old.state.wall_relaxation_memory(0, 0).unwrap(),
    );
    let work_values = |work: MechanicalTimeWork| {
        [
            work.initial_energy_j,
            work.final_energy_j,
            work.delta_energy_j,
            work.acoustic_delta_energy_j,
            work.wall_delta_energy_j,
            work.input_work_j,
            work.boundary_dissipation_j,
            work.wall_dissipation_j,
            work.wall_relaxation_dissipation_j,
            work.wall_work_j,
            work.backward_error,
        ]
    };
    assert_same_bits(&work_values(new.work), &work_values(old.work));
}
