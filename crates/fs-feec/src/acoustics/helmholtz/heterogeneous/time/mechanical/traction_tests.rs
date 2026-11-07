//! G1/G3 applied mechanical excitation of the original coupled acoustic field.

use super::super::boundary_source_tests::{quadratic_profile, quadratic_work};
use super::traction::WallTractionProfile;
use super::*;
use crate::acoustics::helmholtz::HelmholtzFe;
use crate::acoustics::helmholtz::heterogeneous::time::tests::{MEDIUM, brick, p2, plane, uniform};
use fs_rep_mesh::TetComplex;

fn profile(
    owner: &PreparedMechanicalTime<'_, '_, '_>,
    coefficients: &[f64],
) -> WallTractionProfile {
    WallTractionProfile {
        wall: 0,
        coefficients: owner
            .wall_dofs(0)
            .unwrap()
            .iter()
            .map(|&global| (global, coefficients[global]))
            .collect(),
    }
}

fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() < 3e-9 * scale,
        "actual={actual:e}, expected={expected:e}, scale={scale:e}"
    );
}

#[test]
fn continuum_translation_separates_external_traction_from_acoustic_pressure_work() {
    // Exact continuum solution: u=v, p=p0, psi=-rho*v*x+p0*t,
    // w=w0+v*t. Applied traction r*v+K*w-p0 supplies the wall force that
    // is not supplied by fluid pressure. Midpoint reproduces this polynomial.
    let pin = WallPin {
        surface_density: 0.4,
        stiffness_per_area: 1.7,
        resistance: 0.3,
    };
    for height in [0.1, 0.3] {
        let (mesh, positions) = brick(2, [1.0, height, 0.2]);
        let area = height * 0.2;
        let geometry = HelmholtzFe::new(&mesh, &positions, 2, MEDIUM).unwrap();
        let model = uniform(&geometry);
        let h = 0.0125;
        let owner = model
            .prepare_time_with_walls(
                h,
                &[],
                &[MechanicalBoundary {
                    faces: plane(&geometry, 1.0),
                    wall: pin,
                    relaxation: vec![],
                }],
                MechanicalTimeLimits::default(),
                &mut || false,
            )
            .unwrap();
        let drive = owner
            .prepare_wall_traction(&profile(&owner, &p2(&geometry, |_| 1.0)), &mut || false)
            .unwrap();
        for amplitude in [1.0, 1e-8] {
            let v = 0.02 * amplitude;
            let w0 = 0.03 * amplitude;
            let p0 = 0.09 * amplitude;
            let restrict = |field: Vec<f64>| {
                owner
                    .wall_dofs(0)
                    .unwrap()
                    .iter()
                    .map(|&i| field[i])
                    .collect::<Vec<_>>()
            };
            let mut state = owner
                .state_with_walls(
                    &p2(&geometry, |x| -MEDIUM.density_kg_m3 * v * x[0]),
                    &p2(&geometry, |_| p0),
                    &[restrict(p2(&geometry, |_| w0))],
                    &[restrict(p2(&geometry, |_| v))],
                )
                .unwrap();
            let source = plane(&geometry, 0.0)
                .into_iter()
                .map(|face| (face, [v; 3]))
                .collect::<Vec<_>>();
            let initial = owner.energy(&state).unwrap().total_j();
            let mut work = 0.0;
            let mut heat = 0.0;
            for step in 1..=8 {
                let time = step as f64 * h;
                let midpoint_w = w0 + v * (time - h / 2.0);
                let traction = pin.resistance * v + pin.stiffness_per_area * midpoint_w - p0;
                let next = owner
                    .step_with_wall_tractions(
                        &state,
                        &source,
                        &[],
                        &[],
                        &[(&drive, traction)],
                        &mut || false,
                    )
                    .unwrap();
                for (&actual, expected) in
                    next.state.pressure_pa().iter().zip(p2(&geometry, |_| p0))
                {
                    close(actual, expected, amplitude);
                }
                for (&actual, expected) in next
                    .state
                    .wall_displacement_m(0)
                    .unwrap()
                    .iter()
                    .zip(restrict(p2(&geometry, |_| w0 + v * time)))
                {
                    close(actual, expected, amplitude);
                }
                let scale = h * area * amplitude * amplitude;
                close(next.work.acoustic_input_work_j, h * area * p0 * v, scale);
                close(next.work.wall_work_j, h * area * p0 * v, scale);
                close(
                    next.work.wall_traction_work_j,
                    h * area * traction * v,
                    scale,
                );
                close(
                    next.work.wall_dissipation_j,
                    h * area * pin.resistance * v * v,
                    scale,
                );
                assert!(next.work.wall_traction_work_j < 0.0);
                assert!(next.work.acoustic_input_work_j > 0.0);
                close(next.work.acoustic_delta_energy_j, 0.0, scale);
                for residual in [
                    next.work.balance_residual_j(),
                    next.work.acoustic_balance_residual_j(),
                    next.work.wall_balance_residual_j(),
                ] {
                    close(residual, 0.0, scale);
                }
                work += next.work.input_work_j;
                heat += next.work.wall_dissipation_j;
                close(
                    owner.energy(&next.state).unwrap().total_j() + heat - initial,
                    work,
                    initial,
                );
                state = next.state;
            }
        }
    }
}

#[test]
fn nonuniform_traction_drives_mechanical_maxwell_and_spatial_absorptive_states() {
    let points = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &points, 1, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let residue = model
        .spatial_resistive_boundary(
            &[[0, 1, 3]],
            4,
            &mut |_, _| Ok::<_, crate::acoustics::helmholtz::AcousticFeError>(0.7),
            &mut || false,
        )
        .unwrap();
    let law = model
        .spatial_relaxing_boundary(None, &[(&residue, 0.03)], &mut || false)
        .unwrap();
    let h = 0.001;
    let owner = model
        .prepare_time_with_walls_and_relaxing_boundaries(
            h,
            &[],
            &[],
            &[MechanicalBoundary {
                faces: vec![[0, 1, 2]],
                wall: WallPin {
                    surface_density: 0.3,
                    stiffness_per_area: 7.0,
                    resistance: 0.2,
                },
                relaxation: vec![AcousticWallRelaxation {
                    stiffness_per_area: 2.0,
                    relaxation_time_s: 0.04,
                }],
            }],
            &[&law],
            &[],
            MechanicalTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    // This shape has ZERO net force on the triangle. Its spatial motion and
    // sound must survive; a patch-average piston would produce silence.
    let drive = owner
        .prepare_wall_traction(
            &WallTractionProfile {
                wall: 0,
                coefficients: vec![(0, 1.0), (1, -1.0)],
            },
            &mut || false,
        )
        .unwrap();
    let mut state = owner.zero_state();
    let mut input = 0.0;
    let mut heat = 0.0;
    for _ in 0..16 {
        let next = owner
            .step_with_wall_tractions(&state, &[], &[], &[], &[(&drive, 2.0)], &mut || false)
            .unwrap();
        let velocity = state
            .wall_normal_velocity_m_s(0)
            .unwrap()
            .iter()
            .zip(next.state.wall_normal_velocity_m_s(0).unwrap())
            .map(|(&a, &b)| f64::midpoint(a, b))
            .collect::<Vec<_>>();
        // Independent P1 triangle moments W_ii=1/12, W_ij=1/24.
        let force = [2.0, -2.0, 0.0];
        let mut expected = 0.0;
        for i in 0..3 {
            for j in 0..3 {
                expected += h * force[i] * velocity[j] / if i == j { 12.0 } else { 24.0 };
            }
        }
        close(
            next.work.wall_traction_work_j,
            expected,
            expected.abs().max(1e-10),
        );
        assert_eq!(next.work.acoustic_input_work_j, 0.0);
        assert_eq!(next.work.input_work_j, next.work.wall_traction_work_j);
        for residual in [
            next.work.balance_residual_j(),
            next.work.acoustic_balance_residual_j(),
            next.work.wall_balance_residual_j(),
            next.work.admittance_balance_residual_j(),
        ] {
            close(residual, 0.0, next.work.input_work_j.abs().max(1e-10));
        }
        input += next.work.input_work_j;
        heat += next.work.wall_dissipation_j + next.work.admittance_dissipation_j;
        state = next.state;
    }
    let energy = owner.energy(&state).unwrap();
    assert!(
        energy.acoustic.total_j() > 0.0
            && energy.wall_relaxation_j > 0.0
            && energy.admittance_memory_j > 0.0
    );
    assert!(heat > 0.0);
    close(energy.total_j() + heat, input, input);
}

#[test]
fn quadratic_traction_retains_weak_coherent_force_and_exact_surface_work() {
    let lengths = [1.0, 0.2, 0.3];
    let (mesh, points) = brick(2, lengths);
    let geometry = HelmholtzFe::new(&mesh, &points, 2, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let faces = plane(&geometry, 1.0);
    let h = 0.005;
    let owner = model
        .prepare_time_with_walls(
            h,
            &[],
            &[MechanicalBoundary {
                faces: faces.clone(),
                wall: WallPin {
                    surface_density: 0.4,
                    stiffness_per_area: 2.0,
                    resistance: 0.1,
                },
                relaxation: vec![],
            }],
            MechanicalTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let shape = |x: [f64; 3]| {
        let y = x[1] / lengths[1];
        4.0 * y * (1.0 - y)
    };
    let coefficients = p2(&geometry, shape);
    assert!(coefficients[..points.len()].iter().all(|&x| x == 0.0));
    let drive = owner
        .prepare_wall_traction(&profile(&owner, &coefficients), &mut || false)
        .unwrap();
    let zero = owner.zero_state();
    let single = owner
        .step_with_wall_tractions(&zero, &[], &[], &[], &[(&drive, 1.0)], &mut || false)
        .unwrap();
    let coherent = owner
        .step_with_wall_tractions(
            &zero,
            &[],
            &[],
            &[],
            &[(&drive, 1e16), (&drive, 1.0), (&drive, -1e16)],
            &mut || false,
        )
        .unwrap();
    assert_eq!(single.state.pressure_pa(), coherent.state.pressure_pa());
    assert_eq!(
        single.state.wall_normal_velocity_m_s(0),
        coherent.state.wall_normal_velocity_m_s(0)
    );
    assert_eq!(
        single.work.wall_traction_work_j,
        coherent.work.wall_traction_work_j
    );
    assert!(single.work.wall_traction_work_j > 0.0);
    let before = vec![0.0; geometry.ndof()];
    let mut after = before.clone();
    for (&global, &velocity) in owner
        .wall_dofs(0)
        .unwrap()
        .iter()
        .zip(single.state.wall_normal_velocity_m_s(0).unwrap())
    {
        after[global] = velocity;
    }
    // The helper integrates the P2 polynomial product with an independent
    // exact nodal triangle mass; here that product is traction times velocity.
    let expected = quadratic_work(
        &geometry,
        &quadratic_profile(&geometry, &faces, shape),
        &before,
        &after,
        1.0,
        h,
    );
    close(single.work.wall_traction_work_j, expected, expected);
}

#[test]
fn traction_preserves_release_shared_wall_coordinates_and_refuses_wrong_ownership() {
    let points = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &points, 1, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let wall = MechanicalBoundary {
        faces: vec![[0, 1, 2]],
        wall: WallPin {
            surface_density: 0.3,
            stiffness_per_area: 7.0,
            resistance: 0.2,
        },
        relaxation: vec![],
    };
    let prepare = || {
        model
            .prepare_time_with_walls_and_pressure_release(
                0.001,
                &[],
                &[wall.clone()],
                &[[0, 1, 3]],
                MechanicalTimeLimits::default(),
                &mut || false,
            )
            .unwrap()
    };
    let owner = prepare();
    let other = prepare();
    let shape = WallTractionProfile {
        wall: 0,
        coefficients: vec![(0, 1.0)],
    };
    let drive = owner.prepare_wall_traction(&shape, &mut || false).unwrap();
    let wrong = other.prepare_wall_traction(&shape, &mut || false).unwrap();
    let zero = owner.zero_state();
    let next = owner
        .step_with_wall_tractions(&zero, &[], &[], &[], &[(&drive, 1.0)], &mut || false)
        .unwrap();
    for dof in [0, 1, 3] {
        assert_eq!(next.state.pressure_pa()[dof], 0.0);
    }
    assert!(next.state.wall_normal_velocity_m_s(0).unwrap()[0] > 0.0);
    assert_ne!(next.state.pressure_pa()[2], 0.0);
    for invalid in [
        WallTractionProfile {
            wall: 1,
            coefficients: vec![],
        },
        WallTractionProfile {
            wall: 0,
            coefficients: vec![(3, 0.0)],
        },
        WallTractionProfile {
            wall: 0,
            coefficients: vec![(0, 1.0), (0, 2.0)],
        },
        WallTractionProfile {
            wall: 0,
            coefficients: vec![(0, f64::NAN)],
        },
    ] {
        assert!(
            owner
                .prepare_wall_traction(&invalid, &mut || false)
                .is_err()
        );
    }
    assert!(
        owner
            .step_with_wall_tractions(&zero, &[], &[], &[], &[(&wrong, 1.0)], &mut || false)
            .is_err()
    );
    assert!(
        owner
            .step_with_wall_tractions(&zero, &[], &[], &[], &[(&drive, f64::NAN)], &mut || false)
            .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        owner.step_with_wall_tractions(&zero, &[], &[], &[], &[(&drive, 1.0)], &mut || {
            polls += 1;
            polls == 5
        }),
        Err(HeterogeneousError::Cancelled)
    ));
    assert!(zero.pressure_pa().iter().all(|&p| p == 0.0));
    let replay = owner
        .step_with_wall_tractions(&zero, &[], &[], &[], &[(&drive, 1.0)], &mut || false)
        .unwrap();
    assert_eq!(replay.state.pressure_pa(), next.state.pressure_pa());
    let legacy = owner.step(&zero, &[], &mut || false).unwrap();
    let empty = owner
        .step_with_wall_tractions(&zero, &[], &[], &[], &[], &mut || false)
        .unwrap();
    assert_eq!(legacy.state.pressure_pa(), empty.state.pressure_pa());
    assert_eq!(empty.work.wall_traction_work_j, 0.0);
}
