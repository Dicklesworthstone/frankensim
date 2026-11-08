//! Original material histories on the same resolved fluid/structure trajectory.

use super::*;
use crate::acoustics::thermoviscous::moving_wall::relaxation::MovingWallRelaxationLimits;

fn memory_limits() -> MovingWallRelaxationLimits {
    MovingWallRelaxationLimits {
        maximum_channels: 16,
        maximum_workspace_bytes: 128 * 1024,
        maximum_projection_terms: 4096,
    }
}

fn tetrahedron() -> (TetComplex, Vec<[f64; 3]>, f64) {
    let size = [0.002, 0.001, 0.001];
    (
        TetComplex::from_tets(4, vec![[0, 1, 2, 3]]),
        vec![
            [0.0, 0.0, 0.0],
            [size[0], 0.0, 0.0],
            [0.0, size[1], 0.0],
            [0.0, 0.0, size[2]],
        ],
        size.iter().product::<f64>() / 6.0,
    )
}

fn rigid_translations<'c>(
    mesh: &'c TetComplex,
    positions: &[[f64; 3]],
    volume: f64,
    modes: usize,
) -> (ThermoviscousFe<'c>, f64) {
    // Manufactured zero-expansion liquid. On one tetrahedron every P2 velocity
    // DOF is on the boundary. Uniform vector translations therefore constrain
    // the complete velocity field exactly, with no free Cartesian velocity row.
    // Each orthogonal trace carries added fluid mass rho*V*amplitude^2 = 1;
    // original structural modal mass is 1, giving the independent exact M=2.
    let medium = ThermoviscousMedium {
        density_kg_m3: 1000.0,
        pressure_pa: 101325.0,
        temperature_k: 300.0,
        cp_j_kg_k: 4000.0,
        dynamic_viscosity_pa_s: 0.001,
        thermal_conductivity_w_m_k: 0.6,
        bulk_viscosity_pa_s: 0.002,
    };
    let thermo = ThermoviscousThermodynamics {
        isothermal_compressibility_pa_inv: 1.0 / 2e9,
        isobaric_expansion_per_k: 0.0,
    };
    let geometry = HelmholtzFe::new(
        mesh,
        positions,
        2,
        Medium {
            density_kg_m3: medium.density_kg_m3,
            sound_speed_m_s: thermo.sound_speed_m_s(medium).unwrap(),
        },
    )
    .unwrap();
    let faces = geometry.boundary_faces();
    let amplitude = 1.0 / det::sqrt(medium.density_kg_m3 * volume);
    let mut shapes = vec![vec![[0.0; 3]; SimplexSpace::new(mesh, 2).ndof]; modes];
    for (axis, shape) in shapes.iter_mut().enumerate() {
        for vector in &mut shape[..mesh.vertex_count] {
            vector[axis] = amplitude;
        }
    }
    let structure = ModalAcousticTimeModel::try_new_with_free_coordinates(
        48_000,
        vec![
            ModalAcousticMode {
                angular_frequency_rad_s: 0.0,
                damping_ratio: 0.0,
                pressure_per_modal_velocity: C64::ZERO,
            };
            modes
        ],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    let model = ThermoviscousFe::new_with_moving_wall(
        mesh,
        positions,
        medium,
        thermo,
        &ThermoviscousBoundary {
            no_slip_faces: faces.clone(),
            isothermal_faces: vec![],
            ports: vec![],
        },
        MovingWall {
            faces,
            velocity_shapes: shapes,
            structure,
        },
        limits(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(model.velocity_count, 0);
    (model, amplitude)
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-8 * expected.abs().max(1e-30),
        "{actual:e} != {expected:e}",
    );
}

#[test]
fn g1_original_preload_and_dark_history_release_the_exact_energy_into_motion_and_heat() {
    let (mesh, positions, volume) = tetrahedron();
    let (model, _) = rigid_translations(&mesh, &positions, volume, 1);
    let (c, tau, dark_tau) = (40.0, 0.05, 0.03);
    let model = model
        .with_moving_wall_relaxation_channels(
            vec![vec![c], vec![0.0]],
            vec![tau, dark_tau],
            memory_limits(),
            &mut || false,
        )
        .unwrap();
    let dt = 0.017;
    let prepared = model
        .prepare_time(dt, ThermoviscousTimeLimits::default(), &mut || false)
        .unwrap();
    let mut initial = prepared.zero_state();
    assert_eq!(initial.wall_relaxation_coordinates_sqrt_j, vec![0.0, 0.0]);
    let (preload, dark_preload) = (7e-7, -2e-7);
    initial.wall_relaxation_coordinates_sqrt_j = vec![preload, dark_preload];
    let snapshot = initial.clone();
    let step = prepared.step(&initial, &[], &[0.0], &mut || false).unwrap();
    assert_eq!(initial, snapshot);

    // Independent closed solution of M*v'=-c*x, x'=c*v-x/tau,
    // q'=v with M=2 and q(0)=v(0)=0. The dark channel only relaxes.
    let mass = 2.0;
    let a = 0.5 / tau;
    let frequency = det::sqrt(c * c / mass - a * a);
    let decay = (-a * dt).exp();
    let sine = (frequency * dt).sin();
    let cosine = (frequency * dt).cos();
    let velocity = -c * preload / (mass * frequency) * decay * sine;
    let displacement = -preload / c * (1.0 - decay * (cosine + a / frequency * sine));
    let memory = preload * decay * (cosine - a / frequency * sine);
    let dark = dark_preload * (-dt / dark_tau).exp();
    close(step.state.wall_velocity_m_sqrt_kg_per_s[0], velocity);
    close(step.state.wall_displacement_m_sqrt_kg[0], displacement);
    close(step.state.wall_relaxation_coordinates_sqrt_j[0], memory);
    close(step.state.wall_relaxation_coordinates_sqrt_j[1], dark);
    let initial_energy = 0.5 * (preload * preload + dark_preload * dark_preload);
    let final_energy = 0.5 * (mass * velocity * velocity + memory * memory + dark * dark);
    close(prepared.stored_energy_j(&initial).unwrap(), initial_energy);
    close(step.work.initial_energy_j, initial_energy);
    close(step.work.final_energy_j, final_energy);
    close(
        step.work.structural_relaxation_heat_j,
        initial_energy - final_energy,
    );
    assert!(step.work.structural_relaxation_heat_j > 0.0);
    assert!(step.work.supplied_j().abs() < 1e-12 * initial_energy);
    assert!(
        (step.work.dissipated_j() - step.work.structural_relaxation_heat_j).abs()
            < 2e-8 * initial_energy
    );

    // The entire original history survives two accepted intervals. Heat and
    // endpoint agree with the same exact full-state trajectory at one interval.
    let half = model
        .prepare_time(0.5 * dt, ThermoviscousTimeLimits::default(), &mut || false)
        .unwrap();
    let first = half.step(&initial, &[], &[0.0], &mut || false).unwrap();
    let second = half.step(&first.state, &[], &[0.0], &mut || false).unwrap();
    close(second.state.wall_velocity_m_sqrt_kg_per_s[0], velocity);
    close(second.state.wall_displacement_m_sqrt_kg[0], displacement);
    close(second.state.wall_relaxation_coordinates_sqrt_j[0], memory);
    close(second.state.wall_relaxation_coordinates_sqrt_j[1], dark);
    close(
        first.work.structural_relaxation_heat_j + second.work.structural_relaxation_heat_j,
        step.work.structural_relaxation_heat_j,
    );
}

#[test]
fn g1_signed_cross_mode_memory_preserves_undamped_common_motion_and_opposite_reaction() {
    let (mesh, positions, volume) = tetrahedron();
    let (model, amplitude) = rigid_translations(&mesh, &positions, volume, 2);
    let model = model
        .with_moving_wall_relaxation_channels(
            vec![vec![30.0, -30.0]],
            vec![0.05],
            memory_limits(),
            &mut || false,
        )
        .unwrap();
    let dt = 0.01;
    let prepared = model
        .prepare_time(dt, ThermoviscousTimeLimits::default(), &mut || false)
        .unwrap();
    let mut common = prepared.zero_state();
    let common_velocity = 2.5e-8;
    common.wall_velocity_m_sqrt_kg_per_s = vec![common_velocity, common_velocity];
    for vector in &mut common.velocity_m_s[..mesh.vertex_count] {
        *vector = [
            common_velocity * amplitude,
            common_velocity * amplitude,
            0.0,
        ];
    }
    let step = prepared
        .step(&common, &[], &[0.0, 0.0], &mut || false)
        .unwrap();
    for &velocity in &step.state.wall_velocity_m_sqrt_kg_per_s {
        close(velocity, common_velocity);
    }
    for &displacement in &step.state.wall_displacement_m_sqrt_kg {
        close(displacement, common_velocity * dt);
    }
    assert!(step.state.wall_relaxation_coordinates_sqrt_j[0].abs() < 1e-12);
    assert!(step.work.structural_relaxation_heat_j.abs() < 1e-12 * step.work.initial_energy_j);
    assert!(step.work.dissipated_j().abs() < 2e-8 * step.work.initial_energy_j);

    let mut preload = prepared.zero_state();
    preload.wall_relaxation_coordinates_sqrt_j[0] = 3e-7;
    let reaction = prepared
        .step(&preload, &[], &[0.0, 0.0], &mut || false)
        .unwrap();
    let velocity = &reaction.state.wall_velocity_m_sqrt_kg_per_s;
    assert!(velocity[0] < 0.0 && velocity[1] > 0.0);
    assert!((velocity[0] + velocity[1]).abs() < 2e-8 * velocity[0].abs());
    assert!(reaction.work.structural_relaxation_heat_j > 0.0);
    assert!(reaction.work.balance_residual_j().abs() < 2e-8 * reaction.work.initial_energy_j);
}

#[test]
fn g1_heat_source_trajectory_retains_the_same_material_history_and_heat() {
    let (mesh, positions, volume) = tetrahedron();
    let (model, _) = rigid_translations(&mesh, &positions, volume, 1);
    let model = model
        .with_moving_wall_relaxation_channels(
            vec![vec![40.0]],
            vec![0.05],
            memory_limits(),
            &mut || false,
        )
        .unwrap();
    let source = model
        .volume_heat_source(4, &mut |_, _| Ok(1.0 / volume), &mut || false)
        .unwrap();
    let dt = 0.01;
    let ordinary = model
        .prepare_time(dt, ThermoviscousTimeLimits::default(), &mut || false)
        .unwrap();
    let heated = model
        .prepare_time_with_heat_sources(
            dt,
            model.medium().cp_j_kg_k,
            &[],
            None,
            None,
            &[&source],
            ThermoviscousTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let mut initial = ordinary.zero_state();
    initial.wall_relaxation_coordinates_sqrt_j[0] = 7e-7;
    let expected = ordinary.step(&initial, &[], &[0.0], &mut || false).unwrap();
    let actual = heated
        .step(&initial, &[], &[0.0], &[0.0], &mut || false)
        .unwrap()
        .field;
    close(
        actual.state.wall_relaxation_coordinates_sqrt_j[0],
        expected.state.wall_relaxation_coordinates_sqrt_j[0],
    );
    close(
        actual.state.wall_velocity_m_sqrt_kg_per_s[0],
        expected.state.wall_velocity_m_sqrt_kg_per_s[0],
    );
    close(
        actual.work.structural_relaxation_heat_j,
        expected.work.structural_relaxation_heat_j,
    );
    close(actual.work.final_energy_j, expected.work.final_energy_j);
    assert_eq!(actual.work.heat_source_work_j, 0.0);
}

#[test]
fn g0_complete_history_capacity_finiteness_and_post_propagation_cancellation_refuse() {
    let (mesh, positions, volume) = tetrahedron();
    let (plain, _) = rigid_translations(&mesh, &positions, volume, 1);
    let without_memory = plain
        .prepare_time(0.01, ThermoviscousTimeLimits::default(), &mut || false)
        .unwrap()
        .system()
        .state_dim();
    let model = plain
        .with_moving_wall_relaxation_channels(
            vec![vec![40.0]],
            vec![0.05],
            memory_limits(),
            &mut || false,
        )
        .unwrap();
    assert!(
        model
            .prepare_time(
                0.01,
                ThermoviscousTimeLimits {
                    max_states: without_memory,
                    ..ThermoviscousTimeLimits::default()
                },
                &mut || false
            )
            .is_err()
    );
    let prepared = model
        .prepare_time(0.01, ThermoviscousTimeLimits::default(), &mut || false)
        .unwrap();
    assert_eq!(prepared.system().state_dim(), without_memory + 1);
    let mut initial = prepared.zero_state();
    initial.wall_relaxation_coordinates_sqrt_j[0] = 7e-7;
    let original = initial.clone();
    let mut calls = 0;
    assert!(matches!(
        prepared.step(&initial, &[], &[0.0], &mut || {
            calls += 1;
            calls == 2
        }),
        Err(ThermoviscousTimeError::Fluid(ThermoviscousError::Cancelled))
    ));
    assert_eq!(initial, original);
    assert_eq!(calls, 2);
    let mut missing = initial.clone();
    missing.wall_relaxation_coordinates_sqrt_j.clear();
    assert!(prepared.step(&missing, &[], &[0.0], &mut || false).is_err());
    let mut nonfinite = initial.clone();
    nonfinite.wall_relaxation_coordinates_sqrt_j[0] = f64::INFINITY;
    assert!(prepared.stored_energy_j(&nonfinite).is_err());
}

#[test]
fn g1_empty_memory_attachment_preserves_the_original_full_trajectory() {
    let (mesh, positions, volume) = tetrahedron();
    let (plain, _) = rigid_translations(&mesh, &positions, volume, 1);
    let (empty, _) = rigid_translations(&mesh, &positions, volume, 1);
    let empty = empty
        .with_moving_wall_relaxation_channels(vec![], vec![], memory_limits(), &mut || false)
        .unwrap();
    let a = plain
        .prepare_time(0.01, ThermoviscousTimeLimits::default(), &mut || false)
        .unwrap();
    let b = empty
        .prepare_time(0.01, ThermoviscousTimeLimits::default(), &mut || false)
        .unwrap();
    let state = a.zero_state();
    let first = a.step(&state, &[], &[1e-6], &mut || false).unwrap();
    let second = b.step(&state, &[], &[1e-6], &mut || false).unwrap();
    assert_eq!(first.state, second.state);
    assert_eq!(first.work.initial_energy_j, second.work.initial_energy_j);
    assert_eq!(first.work.final_energy_j, second.work.final_energy_j);
    assert_eq!(first.work.structural_work_j, second.work.structural_work_j);
    assert_eq!(first.work.structural_relaxation_heat_j, 0.0);
    assert_eq!(second.work.structural_relaxation_heat_j, 0.0);
}
