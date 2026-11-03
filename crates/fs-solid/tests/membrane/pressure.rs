use super::*;
use fs_solid::membrane::vibration::VibrationSettings;

#[test]
fn follower_pressure_tracks_current_area_and_rigid_rotation() {
    let p = patch(8000.0, 8000.0, 4000.0, 0.2, 0.0);
    let stretch = 1.07;
    let pressure = 1200.0;
    let (s, c) = 0.63_f64.sin_cos();
    let rotate = |v: [f64; 3]| [v[0], c * v[1] - s * v[2], s * v[1] + c * v[2]];
    let u: Vec<_> = p
        .mesh()
        .nodes
        .iter()
        .map(|x| [(stretch - 1.0) * x[0], (stretch - 1.0) * x[1], 0.0])
        .collect();
    let rotated: Vec<_> = p
        .mesh()
        .nodes
        .iter()
        .zip(&u)
        .map(|(x, u)| {
            let y = rotate([x[0] + u[0], x[1] + u[1], u[2]]);
            [y[0] - x[0] + 0.7, y[1] - x[1] - 0.3, y[2] + 0.2]
        })
        .collect();
    let force = p.pressure_load_n(&u, pressure, || false).unwrap();
    let changed = p.pressure_load_n(&rotated, pressure, || false).unwrap();
    close(
        force.iter().map(|v| v[2]).sum(),
        pressure * 0.02 * stretch * stretch,
        1e-12,
        1e-13,
    );
    for (original, actual) in force.iter().zip(&changed) {
        for (expected, actual) in rotate(*original).iter().zip(actual) {
            close(*expected, *actual, 1e-12, 1e-12);
        }
    }
    assert!(changed.iter().any(|v| v[1].abs() > 1.0));
}

#[test]
fn follower_pressure_tangent_matches_every_force_derivative() {
    let p = patch(8000.0, 2000.0, 500.0, 0.2, 0.37);
    let u = displaced();
    let pressure = -1400.0;
    let h = 1e-7;
    for dof in 0..15 {
        let mut lo = u.clone();
        let mut hi = u.clone();
        lo[dof / 3][dof % 3] -= h;
        hi[dof / 3][dof % 3] += h;
        let a = p.pressure_load_n(&lo, pressure, || false).unwrap();
        let b = p.pressure_load_n(&hi, pressure, || false).unwrap();
        let mut direction = vec![[0.0; 3]; 5];
        direction[dof / 3][dof % 3] = 1.0;
        let tangent = p
            .pressure_tangent_apply(&u, &direction, pressure, || false)
            .unwrap();
        for node in 0..5 {
            for axis in 0..3 {
                close(
                    (b[node][axis] - a[node][axis]) / (2.0 * h),
                    tangent[node][axis],
                    2e-8,
                    2e-8,
                );
            }
        }
    }
    let rigid = p
        .pressure_tangent_apply(&u, &vec![[0.3, -0.2, 0.1]; 5], pressure, || false)
        .unwrap();
    assert!(rigid.iter().flatten().all(|x| *x == 0.0));
}

#[test]
fn pressure_work_integrates_changing_area_exactly() {
    let p = patch(8000.0, 8000.0, 4000.0, 0.2, 0.0);
    let zero = vec![[0.0; 3]; 5];
    let dilation = 0.2_f64;
    let rise = 0.01;
    let pressure = 2000.0;
    let u: Vec<_> = p
        .mesh()
        .nodes
        .iter()
        .map(|x| [dilation * x[0], dilation * x[1], rise])
        .collect();
    let actual = p.pressure_work_j(&zero, &u, pressure, || false).unwrap();
    // A(t)=A0(1+d*t)^2 while every point rises at constant normal speed.
    let exact = pressure * 0.02 * rise * (1.0 + dilation + dilation * dilation / 3.0);
    close(actual, exact, 1e-14, 1e-13);
    let midpoint_work = pressure * 0.02 * rise * (1.0 + 0.5 * dilation).powi(2);
    assert!((actual - midpoint_work).abs() > 1e-4);
    close(
        p.pressure_work_j(&u, &zero, pressure, || false).unwrap(),
        -exact,
        1e-14,
        1e-13,
    );
    assert_eq!(p.pressure_work_j(&u, &u, pressure, || false).unwrap(), 0.0);
}

#[test]
fn pressure_loaded_equilibrium_matches_exact_pyramid_and_reactions() {
    let p = patch(8000.0, 8000.0, 4000.0, 0.2, 0.0);
    let h = 0.006_f64;
    let stretch = 1.02_f64;
    let side = 0.1_f64;
    let area = 2.0 * (side * stretch).powi(2);
    let restoring =
        2.0 * 8000.0 * (stretch.powi(2) - 1.0) * h + 4.0 * 8000.0 * h.powi(3) / side.powi(2);
    let pressure = restoring / (area / 3.0);
    let zero = vec![[0.0; 3]; 5];
    for sign in [-1.0, 1.0] {
        let result = p
            .solve_with_pressure(
                &zero,
                sign * pressure,
                &prescribed(stretch),
                &zero,
                settings(),
                || false,
            )
            .unwrap();
        close(result.displacement_m[0][2], sign * h, 1e-10, 1e-10);
        close(result.residual_n[0][2], 0.0, 1e-9, 0.0);
        close(
            result.residual_n.iter().skip(1).map(|r| r[2]).sum(),
            -sign * pressure * area,
            1e-9,
            1e-10,
        );
        close(
            result.pressure_work_j,
            pressure * area * h / 3.0,
            1e-12,
            1e-10,
        );
        close(
            result.potential_energy_j,
            result.strain_energy_j - result.pressure_work_j,
            1e-14,
            1e-13,
        );
        assert_eq!(result.pressure_pa, sign * pressure);
    }
    let plain = p
        .solve(&zero, &prescribed(stretch), &zero, settings(), || false)
        .unwrap();
    let unloaded = p
        .solve_with_pressure(&zero, 0.0, &prescribed(stretch), &zero, settings(), || {
            false
        })
        .unwrap();
    assert_eq!(plain.displacement_m, unloaded.displacement_m);
    assert_eq!(plain.newton, unloaded.newton);
}

#[test]
fn pressure_loaded_modes_include_the_follower_tangent() {
    // Two interior nodes give nonzero pressure stiffness on the free space;
    // a one-centre pyramid alone cannot expose its accidental omission.
    let mut nodes = Vec::new();
    for j in 0..3 {
        for i in 0..4 {
            nodes.push([0.1 * i as f64, 0.1 * j as f64]);
        }
    }
    let mut elems = Vec::new();
    for j in 0..2 {
        for i in 0..3 {
            let a = 4 * j + i;
            elems.extend([vec![a, a + 1, a + 5], vec![a, a + 5, a + 4]]);
        }
    }
    let mut fixed = Vec::new();
    for (node, x) in nodes.iter().enumerate() {
        if node != 5 && node != 6 {
            fixed.extend([
                (3 * node, 0.02 * x[0]),
                (3 * node + 1, 0.02 * x[1]),
                (3 * node + 2, 0.0),
            ]);
        }
    }
    let p = MembranePatch::try_new(
        Mesh2 {
            nodes,
            elems,
            patches: vec![],
        },
        vec![sheet(8000.0, 8000.0, 4000.0, 0.2); 12],
        vec![0.0; 12],
        domain(),
        budget(),
        || false,
    )
    .unwrap();
    let zero = vec![[0.0; 3]; 12];
    let pressure = 200.0;
    let result = p
        .solve_vibration_with_pressure(
            &zero,
            pressure,
            &fixed,
            &zero,
            settings(),
            VibrationSettings {
                max_free_dofs: 6,
                relative_residual_tolerance: 1e-8,
            },
            || false,
        )
        .unwrap();
    let u = &result.equilibrium().displacement_m;
    let mut largest_pressure_action = 0.0_f64;
    for mode in result.modes() {
        let mut direction = zero.clone();
        for (&dof, &phi) in result.free_dofs().iter().zip(&mode.phi) {
            direction[dof / 3][dof % 3] = phi;
        }
        let material = p.tangent_apply(u, &direction, || false).unwrap();
        let load = p
            .pressure_tangent_apply(u, &direction, pressure, || false)
            .unwrap();
        for (&dof, &phi) in result.free_dofs().iter().zip(&mode.phi) {
            largest_pressure_action = largest_pressure_action.max(load[dof / 3][dof % 3].abs());
            close(
                material[dof / 3][dof % 3] - load[dof / 3][dof % 3],
                mode.lambda * p.nodal_mass_kg()[dof / 3] * phi,
                1e-6,
                1e-8,
            );
        }
    }
    assert!(largest_pressure_action > 1.0);
    assert_eq!(result.equilibrium().pressure_pa, pressure);
    assert_eq!(
        result.pressure_load_per_pa(),
        p.pressure_load_n(u, 1.0, || false).unwrap()
    );
}

#[test]
fn pressure_volume_admission_validates_the_public_coordinate_map() {
    let p = patch(8000.0, 8000.0, 4000.0, 0.2, 0.0);
    let mut fixed = vec![None; 15];
    for (dof, value) in prescribed(1.02) {
        fixed[dof] = Some(value);
    }
    p.admit_pressure_perimeter(&fixed, || false).unwrap();
    for invalid in [&fixed[..14], &vec![Some(0.0); 16][..]] {
        assert!(matches!(
            p.admit_pressure_perimeter(invalid, || false),
            Err(MembraneError::Invalid(_))
        ));
    }
    let mut nonfinite = fixed.clone();
    nonfinite[0] = Some(f64::NAN);
    assert!(matches!(
        p.admit_pressure_perimeter(&nonfinite, || false),
        Err(MembraneError::Invalid(_))
    ));
    fixed[5] = None;
    assert!(matches!(
        p.admit_pressure_perimeter(&fixed, || false),
        Err(MembraneError::Invalid(_))
    ));
    assert!(matches!(
        p.admit_pressure_perimeter(&fixed, || true),
        Err(MembraneError::Cancelled)
    ));
}

#[test]
fn pressure_solve_refuses_free_perimeter_and_cancelled_work() {
    let p = patch(8000.0, 8000.0, 4000.0, 0.2, 0.0);
    let zero = vec![[0.0; 3]; 5];
    let mut open = prescribed(1.02);
    open.retain(|(dof, _)| *dof != 5);
    assert!(matches!(
        p.solve_with_pressure(&zero, 100.0, &open, &zero, settings(), || false),
        Err(MembraneError::Invalid(_))
    ));
    assert!(
        p.solve_with_pressure(
            &zero,
            f64::NAN,
            &prescribed(1.02),
            &zero,
            settings(),
            || false
        )
        .is_err()
    );
    assert!(matches!(
        p.pressure_load_n(&zero, 100.0, || true),
        Err(MembraneError::Cancelled)
    ));
    let collapsed: Vec<_> = p.mesh().nodes.iter().map(|x| [-x[0], -x[1], 0.0]).collect();
    assert!(matches!(
        p.pressure_load_n(&collapsed, 100.0, || false),
        Err(MembraneError::State { .. })
    ));
    assert!(matches!(
        p.pressure_tangent_apply(&collapsed, &zero, 100.0, || false),
        Err(MembraneError::State { .. })
    ));
    let mut polls = 0;
    assert!(matches!(
        p.solve_with_pressure(&zero, 100.0, &prescribed(1.02), &zero, settings(), || {
            polls += 1;
            polls == 45
        }),
        Err(MembraneError::Cancelled)
    ));
    assert!(
        p.solve_with_pressure(&zero, 100.0, &prescribed(1.02), &zero, settings(), || false)
            .is_ok()
    );
}
