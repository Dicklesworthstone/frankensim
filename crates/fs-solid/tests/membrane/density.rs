//! G1/G3 reference-density placement, exact inertia, and native membrane modes.
use super::*;
use fs_solid::membrane::dynamics::MembraneMotion;
use fs_solid::membrane::inertia::MembraneMassImage;
use fs_solid::membrane::vibration::VibrationSettings;

fn on_mesh(mesh: Mesh2) -> MembranePatch {
    let count = mesh.elems.len();
    MembranePatch::try_new(
        mesh,
        vec![sheet(8000.0, 8000.0, 4000.0, 0.2); count],
        vec![0.0; count],
        domain(),
        budget(),
        || false,
    )
    .unwrap()
}

fn triangle() -> MembranePatch {
    on_mesh(Mesh2 {
        nodes: vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
        elems: vec![vec![0, 1, 2]],
        patches: vec![],
    })
}

fn kinetic(p: &MembranePatch, image: MembraneMassImage, velocity: &[[f64; 3]]) -> f64 {
    let momentum = p.reference_mass_apply(image, velocity, || false).unwrap();
    0.5 * velocity
        .iter()
        .flatten()
        .zip(momentum.iter().flatten())
        .map(|(v, p)| v * p)
        .sum::<f64>()
}

#[test]
fn graded_areal_mass_has_exact_cubic_inertia_and_reference_inventory() {
    let original = triangle();
    // rho = 2 + 3x + 5y. Independent simplex integration gives the matrix
    // below, total mass 7/3, first moments (19/24, 7/8).
    let p = original
        .clone()
        .with_reference_areal_density(vec![[2.0, 5.0, 7.0]], || false)
        .unwrap();
    let exact = [
        [3.0 / 10.0, 7.0 / 40.0, 23.0 / 120.0],
        [7.0 / 40.0, 2.0 / 5.0, 13.0 / 60.0],
        [23.0 / 120.0, 13.0 / 60.0, 7.0 / 15.0],
    ];
    for column in 0..3 {
        let mut input = [[0.0; 3]; 3];
        input[column][2] = 1.0;
        let result = p
            .reference_mass_apply(MembraneMassImage::Consistent, &input, || false)
            .unwrap();
        for row in 0..3 {
            close(result[row][2], exact[row][column], 2e-16, 2e-15);
            assert_eq!(result[row][0], 0.0);
            assert_eq!(result[row][1], 0.0);
        }
    }
    close(p.reference_mass_kg(), 7.0 / 3.0, 1e-15, 2e-15);
    for (actual, expected) in p
        .nodal_mass_kg()
        .iter()
        .zip([2.0 / 3.0, 19.0 / 24.0, 7.0 / 8.0])
    {
        close(*actual, expected, 1e-15, 2e-15);
    }
    let gravity = p.body_load_n([0.0, 0.0, -9.0]).unwrap();
    close(gravity.iter().map(|f| f[2]).sum(), -21.0, 1e-14, 1e-14);
    close(gravity[1][2], -9.0 * 19.0 / 24.0, 1e-14, 1e-14);
    close(gravity[2][2], -9.0 * 7.0 / 8.0, 1e-14, 1e-14);
    let consistent_gravity = p
        .reference_mass_apply(
            MembraneMassImage::Consistent,
            &[[0.0, 0.0, -9.0]; 3],
            || false,
        )
        .unwrap();
    for (a, b) in gravity
        .iter()
        .flatten()
        .zip(consistent_gravity.iter().flatten())
    {
        close(*a, *b, 1e-14, 1e-14);
    }

    let u = [[0.0; 3], [0.05, 0.0, 0.0], [0.0, 0.03, 0.0]];
    let current = p.current_corner_areal_density_kg_m2(&u, || false).unwrap();
    let jacobian = 1.05 * 1.03;
    for (rho, rho0) in current[0].iter().zip([2.0, 5.0, 7.0]) {
        close(*rho * jacobian, rho0, 2e-15, 2e-15);
    }
    close(
        current[0].iter().sum::<f64>() * (0.5 * jacobian / 3.0),
        7.0 / 3.0,
        1e-15,
        1e-14,
    );
    assert_eq!(p.reference_corner_areal_density_kg_m2(), &[[2.0, 5.0, 7.0]]);
    assert_eq!(p.sheets(), original.sheets());
    assert_eq!(
        p.evaluate(&u, || false).unwrap().strain_energy_j,
        original.evaluate(&u, || false).unwrap().strain_energy_j
    );
}

#[test]
fn equal_inventory_density_placement_changes_actual_supported_modes() {
    let base = patch(8000.0, 8000.0, 4000.0, 0.2, 0.0);
    let central = base
        .clone()
        .with_reference_areal_density(vec![[0.4, 0.1, 0.1]; 4], || false)
        .unwrap();
    let peripheral = base
        .with_reference_areal_density(vec![[0.1, 0.25, 0.25]; 4], || false)
        .unwrap();
    close(
        central.reference_mass_kg(),
        peripheral.reference_mass_kg(),
        1e-18,
        1e-14,
    );
    let fixed = prescribed(1.02);
    let zero = [[0.0; 3]; 5];
    let stiffness = 2.0 * 8000.0 * (1.02_f64.powi(2) - 1.0);
    let modal = VibrationSettings {
        max_free_dofs: 8,
        relative_residual_tolerance: 1e-8,
    };
    for (image, masses) in [
        (MembraneMassImage::Lumped, [0.02 / 12.0, 0.014 / 12.0]),
        (MembraneMassImage::Consistent, [0.028 / 30.0, 0.016 / 30.0]),
    ] {
        let mut eigenvalues = Vec::new();
        for (p, mass) in [&central, &peripheral].into_iter().zip(masses) {
            let modes = p
                .solve_vibration_with_inertia(
                    &zero,
                    0.0,
                    &fixed,
                    &zero,
                    settings(),
                    modal,
                    image,
                    || false,
                )
                .unwrap();
            assert_eq!(modes.inertia_image(), image);
            assert_eq!(modes.free_dofs(), &[2]);
            close(modes.modes()[0].lambda, stiffness / mass, 1e-5, 1e-11);
            close(
                modes.modes()[0].phi[0].abs(),
                1.0 / mass.sqrt(),
                1e-12,
                1e-12,
            );
            let velocity = modes.reconstruct(&[0.003]).unwrap();
            close(
                kinetic(p, image, &velocity),
                0.5 * 0.003_f64.powi(2),
                1e-18,
                1e-12,
            );
            eigenvalues.push(modes.modes()[0].lambda);
        }
        assert!(eigenvalues[0] < eigenvalues[1]);
    }
}

#[test]
fn consistent_native_modes_retain_inter_node_inertia_and_pressure_work() {
    let p = triangle()
        .with_reference_areal_density(vec![[2.0, 5.0, 7.0]], || false)
        .unwrap();
    let stretch = 1.02;
    let fixed = [
        (0, 0.0),
        (1, 0.0),
        (2, 0.0),
        (3, stretch - 1.0),
        (4, 0.0),
        (6, 0.0),
        (7, stretch - 1.0),
    ];
    let modes = p
        .solve_vibration_with_inertia(
            &[[0.0; 3]; 3],
            0.0,
            &fixed,
            &[[0.0; 3]; 3],
            settings(),
            VibrationSettings {
                max_free_dofs: 8,
                relative_residual_tolerance: 1e-8,
            },
            MembraneMassImage::Consistent,
            || false,
        )
        .unwrap();
    assert_eq!(modes.free_dofs(), &[5, 8]);
    let k = 0.5 * 8000.0 * 0.5 * (stretch * stretch - 1.0);
    let center = 13.0 / 30.0;
    let radius = ((1.0_f64 / 30.0).powi(2) + (13.0_f64 / 60.0).powi(2)).sqrt();
    for (mode, mass) in modes.modes().iter().zip([center + radius, center - radius]) {
        close(mode.lambda, k / mass, 1e-9, 1e-11);
    }
    let qdot = [0.02, -0.03];
    let velocity = modes.reconstruct(&qdot).unwrap();
    close(
        kinetic(&p, MembraneMassImage::Consistent, &velocity),
        0.5 * (0.02_f64.powi(2) + 0.03_f64.powi(2)),
        1e-16,
        1e-12,
    );
    let force: Vec<_> = modes
        .pressure_load_per_pa()
        .iter()
        .map(|f| f.map(|v| 2.5 * v))
        .collect();
    let physical: f64 = force
        .iter()
        .flatten()
        .zip(velocity.iter().flatten())
        .map(|(f, v)| f * v)
        .sum();
    let modal: f64 = modes
        .project_nodal_force(&force)
        .unwrap()
        .iter()
        .zip(qdot)
        .map(|(f, v)| f * v)
        .sum();
    close(physical, modal, 1e-16, 1e-12);
}

#[test]
fn affine_density_and_velocity_integrals_survive_partition_and_permutation() {
    let coarse_mesh = triangle().mesh().clone();
    let refined_mesh = Mesh2 {
        nodes: vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [0.5, 0.0],
            [0.5, 0.5],
            [0.0, 0.5],
        ],
        elems: vec![vec![0, 3, 5], vec![3, 1, 4], vec![5, 4, 2], vec![3, 4, 5]],
        patches: vec![],
    };
    let mut reference = None;
    for mut mesh in [coarse_mesh, refined_mesh] {
        for rotate_corners in [false, true] {
            if rotate_corners {
                for element in &mut mesh.elems {
                    element.rotate_left(1);
                }
            }
            let density = mesh
                .elems
                .iter()
                .map(|t| {
                    std::array::from_fn(|a| {
                        let x = mesh.nodes[t[a]];
                        2.0 + 3.0 * x[0] + 5.0 * x[1]
                    })
                })
                .collect();
            let velocity: Vec<_> = mesh
                .nodes
                .iter()
                .map(|x| [1.0 + x[0], 2.0 - x[1], 0.5 + 2.0 * x[0] - 3.0 * x[1]])
                .collect();
            let p = on_mesh(mesh.clone())
                .with_reference_areal_density(density, || false)
                .unwrap();
            let moments = [
                p.reference_mass_kg(),
                p.nodal_mass_kg()
                    .iter()
                    .zip(&mesh.nodes)
                    .map(|(m, x)| m * x[0])
                    .sum(),
                p.nodal_mass_kg()
                    .iter()
                    .zip(&mesh.nodes)
                    .map(|(m, x)| m * x[1])
                    .sum(),
                kinetic(&p, MembraneMassImage::Consistent, &velocity),
            ];
            if let Some(first) = reference {
                for (a, b) in moments.into_iter().zip(first) {
                    close(a, b, 1e-14, 1e-14);
                }
            } else {
                reference = Some(moments);
            }
        }
    }
}

#[test]
fn uniform_corner_density_keeps_original_lumped_modes_and_motion_bits() {
    let original = patch(8000.0, 8000.0, 4000.0, 0.2, 0.0);
    let field = original
        .clone()
        .with_reference_areal_density(vec![[0.2; 3]; 4], || false)
        .unwrap();
    assert_eq!(original.nodal_mass_kg(), field.nodal_mass_kg());
    let zero = [[0.0; 3]; 5];
    let fixed = prescribed(1.02);
    let modal = VibrationSettings {
        max_free_dofs: 8,
        relative_residual_tolerance: 1e-8,
    };
    let a = original
        .solve_vibration(&zero, &fixed, &zero, settings(), modal, || false)
        .unwrap();
    let b = field
        .solve_vibration(&zero, &fixed, &zero, settings(), modal, || false)
        .unwrap();
    assert_eq!(a.modes()[0].lambda.to_bits(), b.modes()[0].lambda.to_bits());
    assert_eq!(a.modes()[0].phi, b.modes()[0].phi);
    let mut u = zero;
    for &(dof, value) in &fixed {
        u[dof / 3][dof % 3] = value;
    }
    u[0][2] = 0.003;
    let mut a = MembraneMotion::new(&original, &fixed, &u, &zero, 0.0, || false).unwrap();
    let mut b = MembraneMotion::new(&field, &fixed, &u, &zero, 0.0, || false).unwrap();
    for _ in 0..8 {
        assert_eq!(
            a.step(1e-5, &zero, || false).unwrap(),
            b.step(1e-5, &zero, || false).unwrap()
        );
        assert_eq!(a.displacement_m(), b.displacement_m());
        assert_eq!(a.velocity_m_s(), b.velocity_m_s());
    }
}

#[test]
fn bad_density_and_cancelled_observations_leave_the_reference_owner_intact() {
    let original = triangle();
    for density in [
        vec![],
        vec![[0.0, 1.0, 2.0]],
        vec![[-1.0; 3]],
        vec![[f64::NAN; 3]],
        vec![[f64::INFINITY; 3]],
        vec![[1.0; 3]; 2],
    ] {
        assert!(
            original
                .clone()
                .with_reference_areal_density(density, || false)
                .is_err()
        );
    }
    assert!(matches!(
        original
            .clone()
            .with_reference_areal_density(vec![[2.0, 5.0, 7.0]], || true),
        Err(MembraneError::Cancelled)
    ));
    let mut polls = 0;
    assert!(matches!(
        original
            .clone()
            .with_reference_areal_density(vec![[2.0, 5.0, 7.0]], || {
                polls += 1;
                polls == 3
            }),
        Err(MembraneError::Cancelled)
    ));
    let p = original
        .clone()
        .with_reference_areal_density(vec![[2.0, 5.0, 7.0]], || false)
        .unwrap();
    assert!(
        p.reference_mass_apply(MembraneMassImage::Consistent, &[[0.0; 3]; 2], || false)
            .is_err()
    );
    assert!(matches!(
        p.reference_mass_apply(MembraneMassImage::Consistent, &[[1.0; 3]; 3], || true),
        Err(MembraneError::Cancelled)
    ));
    let collapsed = [[0.0; 3], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]];
    assert!(matches!(
        p.current_corner_areal_density_kg_m2(&collapsed, || false),
        Err(MembraneError::State { .. })
    ));
    assert_eq!(original.reference_corner_areal_density_kg_m2(), &[[0.2; 3]]);
    assert_eq!(p.reference_corner_areal_density_kg_m2(), &[[2.0, 5.0, 7.0]]);
}
