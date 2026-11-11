use super::*;
use crate::acoustics::thermoviscous::moving_wall::relaxation::MovingWallRelaxationLimits;

const STIFFNESS: f64 = 2000.0;
const TAU: f64 = 0.02;
const SIZE: [f64; 3] = [0.002; 3];

fn memory_limits() -> MovingWallRelaxationLimits {
    MovingWallRelaxationLimits {
        maximum_channels: 16,
        maximum_workspace_bytes: 4096,
        maximum_projection_terms: 1024,
    }
}

fn two_coordinates<'c>(mesh: &'c TetComplex, pos: &[[f64; 3]]) -> ThermoviscousFe<'c> {
    let (medium, thermo) = liquid();
    let geometry = HelmholtzFe::new(
        mesh,
        pos,
        2,
        Medium {
            density_kg_m3: medium.density_kg_m3,
            sound_speed_m_s: thermo.sound_speed_m_s(medium).unwrap(),
        },
    )
    .unwrap();
    let faces = geometry.boundary_faces();
    let amplitude = 1.0 / det::sqrt(medium.density_kg_m3 * SIZE.iter().product::<f64>());
    let mut shapes = vec![vec![[0.0; 3]; SimplexSpace::new(mesh, 2).ndof]; 2];
    for shape in &mut shapes {
        for face in &faces {
            for &v in face {
                shape[v as usize][0] = amplitude;
            }
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
            2
        ],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    ThermoviscousFe::new_with_moving_wall(
        mesh,
        pos,
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
    .unwrap()
}

fn relative_rows() -> Vec<Vec<f64>> {
    let c = det::sqrt(STIFFNESS);
    vec![vec![c, -c]]
}

fn close(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() < 2e-8 * expected.abs().max(1e-25),
        "actual={actual:?}, expected={expected:?}"
    );
}

#[test]
fn g1_original_memory_has_correct_reactive_impedance_storage_and_heat() {
    let (mesh, pos) = box_mesh([1, 1, 1], SIZE);
    let model = two_coordinates(&mesh, &pos)
        .with_moving_wall_relaxation_channels(
            relative_rows(),
            vec![TAU],
            memory_limits(),
            &mut || false,
        )
        .unwrap();
    // Differential motion has identically zero fluid trace. Its closed-form
    // response isolates the complete off-diagonal material impedance without
    // deleting the actual resolved fluid or its common-mode inertial loading.
    for frequency in [1.0, 30.0, 1000.0] {
        let omega = std::f64::consts::TAU * frequency;
        let gain = C64::ONE / C64::new(1.0 / TAU, omega);
        let velocity = C64::ONE / (C64::new(0.0, omega) + gain.scale(2.0 * STIFFNESS));
        let history = gain * velocity.scale(2.0 * det::sqrt(STIFFNESS));
        let prepared = model.prepare(frequency, &mut || false).unwrap();
        for phase in [C64::ONE, C64::new(0.0, 1.0)] {
            let field = prepared
                .response_with_wall_forces(&[], &[phase, -phase], &mut || false)
                .unwrap();
            close(
                field.moving_wall.modal_velocity_m_sqrt_kg_per_s[0],
                velocity * phase,
            );
            close(
                field.moving_wall.modal_velocity_m_sqrt_kg_per_s[1],
                -velocity * phase,
            );
            close(
                field.moving_wall.relaxation_coordinates_sqrt_j[0],
                history * phase,
            );
            let heat = history.norm_sq() / (2.0 * TAU);
            let stored = history.norm_sq() / 4.0;
            assert!((field.power.structural_relaxation_w / heat - 1.0).abs() < 2e-8);
            assert!((field.moving_wall.relaxation_mean_energy_j / stored - 1.0).abs() < 2e-8);
            assert!((field.power.structural_relaxation_w / velocity.re - 1.0).abs() < 2e-8);
            assert_eq!(field.power.structural_damping_w, 0.0);
            assert!(
                field
                    .moving_wall
                    .fluid_force_n_per_sqrt_kg
                    .iter()
                    .all(|f| f.abs() < 2e-8)
            );
            assert_work(&field);
        }
    }
}

#[test]
fn g1_relative_memory_preserves_common_fluid_inertia_and_its_nullspace() {
    let (mesh, pos) = box_mesh([1, 1, 1], SIZE);
    let original = two_coordinates(&mesh, &pos);
    let memory = two_coordinates(&mesh, &pos)
        .with_moving_wall_relaxation_channels(
            relative_rows(),
            vec![TAU],
            memory_limits(),
            &mut || false,
        )
        .unwrap();
    let a = original
        .prepare(30.0, &mut || false)
        .unwrap()
        .response_with_wall_forces(&[], &[C64::ONE; 2], &mut || false)
        .unwrap();
    let b = memory
        .prepare(30.0, &mut || false)
        .unwrap()
        .response_with_wall_forces(&[], &[C64::ONE; 2], &mut || false)
        .unwrap();
    for i in 0..2 {
        close(
            a.moving_wall.modal_velocity_m_sqrt_kg_per_s[i],
            b.moving_wall.modal_velocity_m_sqrt_kg_per_s[i],
        );
        close(
            a.moving_wall.fluid_force_n_per_sqrt_kg[i],
            b.moving_wall.fluid_force_n_per_sqrt_kg[i],
        );
    }
    // With no structural damping and a memory channel that is relative only, a common-mode force F moves both
    // coordinates with the purely reactive velocity u = F / (i omega M), M the real total modal mass, so its real input
    // power Re(conj(F) u) / 2 is zero in exact arithmetic and the fluid moves as a rigid body. What the solver reports
    // is roundoff: on the first run of this test input 4.6e-21 W next to viscous plus bulk dissipation of 7.4e-21 W.
    // `assert_work` scales the balance by the larger of those two, which is itself that roundoff, and cannot certify it.
    // Here the balance is scaled by the apparent power sum(|u| / 2), as `viscosity.rs:108-109` does for its common-mode
    // case, and the viscous and bulk power are pinned to the memory-free model `a`.
    let apparent = b
        .moving_wall
        .modal_velocity_m_sqrt_kg_per_s
        .iter()
        .map(|v| 0.5 * v.abs())
        .sum::<f64>();
    let consistent = |m: &ThermoviscousField| {
        let p = m.power;
        p.viscous_w >= 0.0
            && p.bulk_w >= 0.0
            && p.balance_residual_w().abs() < 2e-8 * apparent
            && (p.viscous_w - a.power.viscous_w).abs() < 2e-8 * a.power.viscous_w
            && (p.bulk_w - a.power.bulk_w).abs() < 2e-8 * a.power.bulk_w
            && p.structural_relaxation_w < 1e-20
            && m.moving_wall.relaxation_mean_energy_j < 1e-20
    };
    assert!(consistent(&b), "{:?}", b.power);
    // Planted negative: a channel that couples the common mode (rows [c, c]) stores and dissipates energy in it, and the
    // same checks refuse it.
    let c = det::sqrt(STIFFNESS);
    let coupled = two_coordinates(&mesh, &pos)
        .with_moving_wall_relaxation_channels(
            vec![vec![c, c]],
            vec![TAU],
            memory_limits(),
            &mut || false,
        )
        .unwrap()
        .prepare(30.0, &mut || false)
        .unwrap()
        .response_with_wall_forces(&[], &[C64::ONE; 2], &mut || false)
        .unwrap();
    assert!(!consistent(&coupled), "{:?}", coupled.power);
}

#[test]
fn g0_original_memory_requires_complete_finite_bounded_channels_and_one_owner() {
    let (mesh, pos) = box_mesh([1, 1, 1], SIZE);
    for (rows, times) in [
        (vec![vec![1.0]], vec![TAU]),
        (vec![vec![1.0, f64::NAN]], vec![TAU]),
        (vec![vec![1e200; 2]], vec![TAU]),
        (vec![vec![1e-200; 2]], vec![TAU]),
        (relative_rows(), vec![0.0]),
        (relative_rows(), vec![f64::INFINITY]),
        (relative_rows(), vec![]),
    ] {
        assert!(
            two_coordinates(&mesh, &pos)
                .with_moving_wall_relaxation_channels(rows, times, memory_limits(), &mut || false)
                .is_err()
        );
    }
    for limits in [
        MovingWallRelaxationLimits {
            maximum_channels: 0,
            ..memory_limits()
        },
        MovingWallRelaxationLimits {
            maximum_workspace_bytes: 0,
            ..memory_limits()
        },
        MovingWallRelaxationLimits {
            maximum_projection_terms: 3,
            ..memory_limits()
        },
    ] {
        assert!(matches!(
            two_coordinates(&mesh, &pos).with_moving_wall_relaxation_channels(
                relative_rows(),
                vec![TAU],
                limits,
                &mut || false
            ),
            Err(ThermoviscousError::Budget { .. })
        ));
    }
    let zero = two_coordinates(&mesh, &pos)
        .with_moving_wall_relaxation_channels(
            vec![],
            vec![],
            MovingWallRelaxationLimits {
                maximum_channels: 0,
                maximum_workspace_bytes: 0,
                maximum_projection_terms: 0,
            },
            &mut || false,
        )
        .unwrap();
    assert_eq!(
        zero.moving_wall_relaxation_channels_per_s(),
        Some([].as_slice())
    );
    assert!(
        zero.with_moving_wall_relaxation_channels(vec![], vec![], memory_limits(), &mut || false)
            .is_err()
    );
    assert!(matches!(
        two_coordinates(&mesh, &pos).with_moving_wall_relaxation_channels(
            relative_rows(),
            vec![TAU],
            memory_limits(),
            &mut || true
        ),
        Err(ThermoviscousError::Cancelled)
    ));
}
