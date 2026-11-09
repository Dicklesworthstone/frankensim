use super::*;
use fs_solid::membrane::{MembranePressure, vibration::VibrationSettings};

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

/// G1/G3: exact P1 traction and its first spatial moment, including signed
/// pressure and a rigidly rotated/translated current surface.
#[test]
fn material_pressure_profile_preserves_force_and_moment() {
    let mesh = Mesh2 {
        nodes: vec![[0.0, 0.0], [2.0, 0.0], [0.0, 1.0]],
        elems: vec![vec![0, 1, 2]],
        patches: vec![],
    };
    let p = MembranePatch::try_new(
        mesh,
        vec![sheet(8000.0, 8000.0, 4000.0, 0.2)],
        vec![0.0],
        domain(),
        budget(),
        || false,
    )
    .unwrap();
    let samples = [[3.0, 9.0, -3.0]];
    let pressure = MembranePressure::TriangleLinear(&samples);
    let zero = [[0.0; 3]; 3];
    let force = p.pressure_load_n(&zero, pressure, || false).unwrap();
    // Integrating p=3+3*x-6*y over x/2+y<=1 gives Fz=3,
    // integral y*p=1/2 and integral x*p=3. Uniform averaging loses this moment.
    for (actual, expected) in force.iter().zip([1.0, 1.5, 0.5]) {
        assert_eq!(*actual, [0.0, 0.0, expected]);
    }
    let cross = |a: [f64; 3], b: [f64; 3]| {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    };
    let (sin, cos) = 0.61_f64.sin_cos();
    let rotate = |x: [f64; 3]| [x[0], cos * x[1] - sin * x[2], sin * x[1] + cos * x[2]];
    let translation = [0.7, -0.4, 0.2];
    let positions: Vec<_> = p
        .mesh()
        .nodes
        .iter()
        .map(|x| {
            let rotated = rotate([x[0], x[1], 0.0]);
            std::array::from_fn::<_, 3, _>(|i| rotated[i] + translation[i])
        })
        .collect();
    let displacement: Vec<_> = positions
        .iter()
        .zip(&p.mesh().nodes)
        .map(|(x, r)| [x[0] - r[0], x[1] - r[1], x[2]])
        .collect();
    let rotated = p
        .pressure_load_n(&displacement, pressure, || false)
        .unwrap();
    let mut moment = [0.0; 3];
    let mut resultant = [0.0; 3];
    for ((x, actual), original) in positions.iter().zip(&rotated).zip(&force) {
        let expected = rotate(*original);
        let torque = cross(*x, *actual);
        for i in 0..3 {
            close(actual[i], expected[i], 1e-14, 1e-13);
            resultant[i] += actual[i];
            moment[i] += torque[i];
        }
    }
    let expected_force = rotate([0.0, 0.0, 3.0]);
    let expected_moment = rotate([0.5, -3.0, 0.0]);
    let translated_moment = cross(translation, expected_force);
    for i in 0..3 {
        close(resultant[i], expected_force[i], 1e-14, 1e-13);
        close(
            moment[i],
            expected_moment[i] + translated_moment[i],
            1e-14,
            1e-13,
        );
    }
    let constant = [[7.0; 3]];
    let field = MembranePressure::TriangleLinear(&constant);
    assert_eq!(
        p.pressure_load_n(&displacement, field, || false).unwrap(),
        p.pressure_load_n(&displacement, 7.0, || false).unwrap()
    );
    assert_eq!(
        p.pressure_tangent_apply(&displacement, &zero, field, || false)
            .unwrap(),
        p.pressure_tangent_apply(&displacement, &zero, 7.0, || false)
            .unwrap()
    );
    assert_eq!(
        p.pressure_work_j(&zero, &displacement, field, || false)
            .unwrap(),
        p.pressure_work_j(&zero, &displacement, 7.0, || false)
            .unwrap()
    );
    for bad in [vec![], vec![[1.0; 3]; 2], vec![[1.0, f64::NAN, 2.0]]] {
        let invalid = MembranePressure::TriangleLinear(&bad);
        assert!(p.pressure_load_n(&zero, invalid, || false).is_err());
        assert!(
            p.pressure_tangent_apply(&zero, &zero, invalid, || false)
                .is_err()
        );
        assert!(p.pressure_work_j(&zero, &zero, invalid, || false).is_err());
    }
    assert!(matches!(
        p.pressure_load_n(&zero, pressure, || true),
        Err(MembraneError::Cancelled)
    ));
}

/// G1/G3: exact follower derivative and drift quadrature for distinct original
/// facet traces. Nonuniform material loading is not a conservative cavity law.
#[test]
fn material_pressure_profile_tangent_and_path_work_match_independent_quadrature() {
    let p = patch(8000.0, 2000.0, 500.0, 0.2, 0.37);
    let samples = [
        [3.0, 9.0, -3.0],
        [-2.0, 6.0, 1.0],
        [9.0, 0.0, 4.0],
        [3.0, -7.0, 2.0],
    ];
    let pressure = MembranePressure::TriangleLinear(&samples);
    let before = displaced();
    let direction = [
        [0.003, -0.002, 0.004],
        [0.001, 0.0, 0.002],
        [0.0, -0.001, 0.003],
        [-0.002, 0.001, 0.0],
        [0.001, 0.002, -0.001],
    ];
    let at = |s: f64| {
        before
            .iter()
            .zip(&direction)
            .map(|(u, d)| std::array::from_fn::<_, 3, _>(|i| u[i] + s * d[i]))
            .collect::<Vec<_>>()
    };
    let after = at(1.0);
    let h = 1e-5;
    let lo = p.pressure_load_n(&at(-h), pressure, || false).unwrap();
    let hi = p.pressure_load_n(&at(h), pressure, || false).unwrap();
    let tangent = p
        .pressure_tangent_apply(&before, &direction, pressure, || false)
        .unwrap();
    for ((lo, hi), exact) in lo
        .iter()
        .flatten()
        .zip(hi.iter().flatten())
        .zip(tangent.iter().flatten())
    {
        close((hi - lo) / (2.0 * h), *exact, 1e-12, 1e-8);
    }
    // Two-point Gauss is exact for the quadratic current-area force along the
    // straight path. It also resolves unequal nodal pressure, unlike p_mean*Q.
    let mut quadrature = 0.0;
    for s in [0.5 - (1.0_f64 / 12.0).sqrt(), 0.5 + (1.0_f64 / 12.0).sqrt()] {
        let force = p.pressure_load_n(&at(s), pressure, || false).unwrap();
        quadrature += 0.5
            * force
                .iter()
                .flatten()
                .zip(direction.iter().flatten())
                .map(|(f, d)| f * d)
                .sum::<f64>();
    }
    let work = p
        .pressure_work_j(&before, &after, pressure, || false)
        .unwrap();
    close(work, quadrature, 1e-17, 1e-13);
    close(
        p.pressure_work_j(&after, &before, pressure, || false)
            .unwrap(),
        -work,
        1e-17,
        1e-13,
    );
    let force = p.pressure_load_n(&before, pressure, || false).unwrap();
    let power: f64 = force
        .iter()
        .flatten()
        .zip(direction.iter().flatten())
        .map(|(f, d)| f * d)
        .sum();
    let differential = p
        .pressure_work_j(&at(-h), &at(h), pressure, || false)
        .unwrap()
        / (2.0 * h);
    close(differential, power, 1e-12, 1e-8);

    // The outer rim stays fixed through this closed loop of the center node.
    let localized = [[0.0, 6.0, 0.0], [0.0; 3], [0.0; 3], [0.0; 3]];
    let mut corners = vec![vec![[0.0; 3]; 5]; 5];
    corners[1][0][0] = 0.01;
    corners[2][0] = [0.01, 0.0, 0.002];
    corners[3][0][2] = 0.002;
    let cycle: f64 = corners
        .windows(2)
        .map(|pair| {
            p.pressure_work_j(
                &pair[0],
                &pair[1],
                MembranePressure::TriangleLinear(&localized),
                || false,
            )
            .unwrap()
        })
        .sum();
    assert!(
        cycle.abs() > 1e-7,
        "nonuniform fixed-rim work cannot be treated as p*dV"
    );
}

fn curved_pressure_segment() -> (MembranePatch, Vec<[f64; 3]>, Vec<[f64; 3]>) {
    let mut nodes = Vec::new();
    let mut elems = Vec::new();
    for y in 0..4 {
        for x in 0..4 {
            nodes.push([0.1 * x as f64, 0.1 * y as f64]);
        }
    }
    for y in 0..3 {
        for x in 0..3 {
            let a = 4 * y + x;
            elems.push(vec![a, a + 1, a + 5]);
            elems.push(vec![a, a + 5, a + 4]);
        }
    }
    let count = elems.len();
    let p = MembranePatch::try_new(
        Mesh2 {
            nodes,
            elems,
            patches: vec![],
        },
        vec![sheet(8000.0, 8000.0, 4000.0, 0.2); count],
        vec![0.0; count],
        domain(),
        budget(),
        || false,
    )
    .unwrap();
    let mut before = vec![[0.0; 3]; 16];
    let mut after = before.clone();
    // The outer perimeter stays fixed, while an interior triangle has three
    // moving vertices. This gives a genuinely cubic swept-volume polynomial.
    for (node, height, increment) in [
        (5, 0.004, [0.0012, -0.0007, 0.0021]),
        (6, 0.007, [-0.0008, 0.0014, -0.0011]),
        (9, 0.008, [0.0006, 0.0009, -0.0017]),
        (10, 0.006, [-0.0015, -0.0004, 0.0008]),
    ] {
        before[node][2] = height;
        after[node] = std::array::from_fn(|i| before[node][i] + increment[i]);
    }
    (p, before, after)
}

/// G1/G3: the exact path-mean traction exchanges the SAME cubic swept volume
/// as its gas port. Midpoint geometry alone demonstrably loses pressure work.
#[test]
fn mean_pressure_traction_closes_curved_cavity_work_without_a_correction() {
    let (p, before, after) = curved_pressure_segment();
    let pressure = 130000.0;
    let delta: Vec<_> = before
        .iter()
        .zip(&after)
        .map(|(a, b)| std::array::from_fn::<_, 3, _>(|i| b[i] - a[i]))
        .collect();
    let mean = p
        .pressure_segment_load_n(&before, &after, pressure, || false)
        .unwrap();
    let work: f64 = mean
        .iter()
        .flatten()
        .zip(delta.iter().flatten())
        .map(|(f, d)| f * d)
        .sum();
    let exact = p
        .pressure_work_j(&before, &after, pressure, || false)
        .unwrap();
    close(work, exact, 1e-13, 1e-12);
    let midpoint: Vec<_> = before
        .iter()
        .zip(&delta)
        .map(|(a, d)| std::array::from_fn::<_, 3, _>(|i| a[i] + 0.5 * d[i]))
        .collect();
    let middle = p.pressure_load_n(&midpoint, pressure, || false).unwrap();
    let midpoint_work: f64 = middle
        .iter()
        .flatten()
        .zip(delta.iter().flatten())
        .map(|(f, d)| f * d)
        .sum();
    assert!(
        (midpoint_work - exact).abs() > 1e-8,
        "fixture must expose midpoint pressure-work loss: {midpoint_work:e}, {exact:e}"
    );
    let reverse = p
        .pressure_segment_load_n(&after, &before, pressure, || false)
        .unwrap();
    let mut quadrature = vec![[0.0; 3]; before.len()];
    for t in [0.5 - (1.0_f64 / 12.0).sqrt(), 0.5 + (1.0_f64 / 12.0).sqrt()] {
        let u: Vec<_> = before
            .iter()
            .zip(&delta)
            .map(|(a, d)| std::array::from_fn::<_, 3, _>(|i| a[i] + t * d[i]))
            .collect();
        let force = p.pressure_load_n(&u, pressure, || false).unwrap();
        for (sum, f) in quadrature.iter_mut().flatten().zip(force.iter().flatten()) {
            *sum += 0.5 * f;
        }
    }
    for ((mean, reverse), quadrature) in mean
        .iter()
        .flatten()
        .zip(reverse.iter().flatten())
        .zip(quadrature.iter().flatten())
    {
        close(*mean, *reverse, 1e-11, 1e-13);
        close(*mean, *quadrature, 1e-11, 1e-13);
    }
    assert_eq!(
        p.pressure_segment_load_n(&before, &before, pressure, || false)
            .unwrap(),
        p.pressure_load_n(&before, pressure, || false).unwrap()
    );
}

/// G1/G3/G4: endpoint differentiation retains its one-half coincident-endpoint
/// limit, signed material profiles, rigid translations and admission behavior.
#[test]
fn mean_pressure_endpoint_tangent_matches_force_differences_and_midpoint_chain_rule() {
    let p = patch(8000.0, 2000.0, 500.0, 0.2, 0.37);
    let before = displaced();
    let change = [
        [0.003, -0.002, 0.004],
        [0.001, 0.0, 0.002],
        [0.0, -0.001, 0.003],
        [-0.002, 0.001, 0.0],
        [0.001, 0.002, -0.001],
    ];
    let after: Vec<_> = before
        .iter()
        .zip(change)
        .map(|(a, d)| std::array::from_fn::<_, 3, _>(|i| a[i] + d[i]))
        .collect();
    let direction = [
        [0.2, -0.1, 0.4],
        [-0.1, 0.3, -0.2],
        [0.1, 0.2, 0.3],
        [0.3, -0.2, -0.1],
        [-0.2, 0.1, 0.2],
    ];
    let samples = [
        [3.0, 9.0, -3.0],
        [-2.0, 6.0, 1.0],
        [9.0, 0.0, 4.0],
        [3.0, -7.0, 2.0],
    ];
    for pressure in [
        MembranePressure::Uniform(-1400.0),
        MembranePressure::TriangleLinear(&samples),
    ] {
        let h = 1e-7;
        let at = |scale: f64| {
            after
                .iter()
                .zip(&direction)
                .map(|(u, d)| std::array::from_fn::<_, 3, _>(|i| u[i] + scale * d[i]))
                .collect::<Vec<_>>()
        };
        let lo = p
            .pressure_segment_load_n(&before, &at(-h), pressure, || false)
            .unwrap();
        let hi = p
            .pressure_segment_load_n(&before, &at(h), pressure, || false)
            .unwrap();
        let tangent = p
            .pressure_segment_tangent_apply(&before, &after, &direction, pressure, || false)
            .unwrap();
        for ((a, b), t) in lo
            .iter()
            .flatten()
            .zip(hi.iter().flatten())
            .zip(tangent.iter().flatten())
        {
            close((b - a) / (2.0 * h), *t, 1e-8, 2e-8);
        }
        let coincident = p
            .pressure_segment_tangent_apply(&before, &before, &direction, pressure, || false)
            .unwrap();
        let instantaneous = p
            .pressure_tangent_apply(&before, &direction, pressure, || false)
            .unwrap();
        for (a, b) in coincident
            .iter()
            .flatten()
            .zip(instantaneous.iter().flatten())
        {
            close(2.0 * a, *b, 1e-13, 1e-13);
        }
        let rigid = p
            .pressure_segment_tangent_apply(
                &before,
                &after,
                &[[0.2, -0.1, 0.3]; 5],
                pressure,
                || false,
            )
            .unwrap();
        assert!(rigid.iter().flatten().all(|v| *v == 0.0));
        let mean = p
            .pressure_segment_load_n(&before, &after, pressure, || false)
            .unwrap();
        let work: f64 = mean
            .iter()
            .flatten()
            .zip(change.iter().flatten())
            .map(|(f, d)| f * d)
            .sum();
        close(
            work,
            p.pressure_work_j(&before, &after, pressure, || false)
                .unwrap(),
            1e-14,
            1e-12,
        );
    }
    assert!(matches!(
        p.pressure_segment_load_n(&before, &after, 1.0, || true),
        Err(MembraneError::Cancelled)
    ));
    assert!(
        p.pressure_segment_tangent_apply(&before, &after, &[[f64::NAN; 3]; 5], 1.0, || false)
            .is_err()
    );
    assert!(
        p.pressure_segment_load_n(
            &before,
            &after,
            MembranePressure::TriangleLinear(&[]),
            || false
        )
        .is_err()
    );
}
