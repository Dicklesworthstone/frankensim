//! Independent accuracy and virtual-work checks for the distributed functional.
use super::*;
use fs_em::spatial::vector_potential::{
    SpatialVectorMagnetostatics, VectorMagneticBoundary, VectorMagneticControls,
    VectorPotentialFace,
};

fn weights(positions: &[[f64; 3]], smooth: bool) -> Vec<f64> {
    positions
        .iter()
        .map(|p| {
            let radius = p[0].abs().max(p[1].abs()).max(p[2].abs());
            let t = (2.0 - radius).clamp(0.0, 1.0);
            if smooth { t * t * (3.0 - 2.0 * t) } else { t }
        })
        .collect()
}
fn load(field: &fs_em::spatial::SpatialMagneticField<'_>, smooth: bool) -> Wrench {
    let w = weights(field.problem().positions_m(), smooth);
    let shell = field
        .problem()
        .prepare_eggshell(&w, force_limits(), &mut || false)
        .unwrap();
    assert!(shell.plateau_cell_count() > 0 && shell.collar_cell_count() > 0);
    assert!(shell.collar_volume_m3() > 0.0);
    assert_eq!(shell.nodal_weights(), &w);
    field
        .eggshell_wrench(&shell, &Se3::identity(), [0.0; 3], &mut || false)
        .unwrap()
}
fn applied(problem: &SpatialMagnetostatics) -> SpatialMagneticBoundary {
    potentials(problem, |p| {
        cube_magnet_potential(p) - 0.4 * p[0] - 0.7 * p[0] * p[2]
    })
}

#[test]
fn distributed_cuboid_force_and_torque_improve_raw_accuracy_and_refine_without_extrapolation() {
    let expected = Wrench::new(Vec3::new(0.0, 0.4, 0.0), Vec3::new(0.7, 0.0, 0.0));
    let mut previous = None;
    for n in [8, 16] {
        let problem = magnet(n);
        let field = problem
            .solve(&applied(&problem), controls(), &mut || false)
            .unwrap();
        let surface = field
            .maxwell_wrench(
                &region(&problem, 1.0),
                &Se3::identity(),
                [0.0; 3],
                force_limits(),
                &mut || false,
            )
            .unwrap();
        let raw_surface_error = vector_error(surface.force, expected.force)
            + vector_error(surface.torque, expected.torque);
        let mut errors = [0.0; 2];
        for (index, smooth) in [false, true].into_iter().enumerate() {
            let actual = load(&field, smooth);
            errors[index] = vector_error(actual.force, expected.force)
                + vector_error(actual.torque, expected.torque);
            eprintln!(
                "eggshell n={n} smooth={smooth} force={:?} torque={:?} error={} surface_error={raw_surface_error}",
                actual.force, actual.torque, errors[index]
            );
            assert!(actual.force.x > 0.0 && actual.torque.y > 0.0);
            // Both are actual raw production outputs: no Richardson estimate
            // or field blending is used to achieve this accuracy comparison.
            assert!(errors[index] < 0.25 * raw_surface_error);
        }
        if let Some(old) = previous {
            let old: [f64; 2] = old;
            for index in 0..2 {
                assert!(
                    errors[index] < 0.6 * old[index],
                    "collar profile {index} must refine: {old:?} -> {errors:?}"
                );
                assert!(
                    errors[index] < 0.06,
                    "independent raw force/torque accuracy: {errors:?}"
                );
            }
        }
        previous = Some(errors);
    }
}

fn signed_volume(p: &[[f64; 3]], t: [u32; 4]) -> f64 {
    let a = p[t[0] as usize];
    let b = std::array::from_fn::<_, 3, _>(|i| p[t[1] as usize][i] - a[i]);
    let c = std::array::from_fn::<_, 3, _>(|i| p[t[2] as usize][i] - a[i]);
    let d = std::array::from_fn::<_, 3, _>(|i| p[t[3] as usize][i] - a[i]);
    (b[0] * (c[1] * d[2] - c[2] * d[1]) - b[1] * (c[0] * d[2] - c[2] * d[0])
        + b[2] * (c[0] * d[1] - c[1] * d[0]))
        / 6.0
}
fn coenergy(field: &fs_em::spatial::SpatialMagneticField<'_>) -> f64 {
    let p = field.problem();
    field.stored_energy_j()
        + p.materials()
            .iter()
            .enumerate()
            .map(|(cell, material)| {
                if let SpatialMagneticMaterial::Recoil(magnet) = material {
                    let h = field.field_strength_a_m()[cell];
                    let br = magnet.remanence_t();
                    signed_volume(p.positions_m(), p.mesh().tets[cell]).abs()
                        * (br[0] * h[0] + br[1] * h[1] + br[2] * h[2])
                } else {
                    0.0
                }
            })
            .sum::<f64>()
}

#[test]
fn scalar_translation_matches_resolved_coenergy_shape_derivative_at_fixed_boundary_potential() {
    let problem = magnet(8);
    let boundary = applied(&problem);
    let field = problem.solve(&boundary, controls(), &mut || false).unwrap();
    let actual = load(&field, false);
    let w = weights(problem.positions_m(), false);
    let shifted = |amount: f64| {
        let p: Vec<_> = problem
            .positions_m()
            .iter()
            .zip(&w)
            .map(|(p, w)| [p[0] + amount * w, p[1], p[2]])
            .collect();
        let changed = SpatialMagnetostatics::new(
            p,
            problem.mesh().tets.clone(),
            problem.materials().to_vec(),
            limits(),
            &mut || false,
        )
        .unwrap();
        // Boundary coordinates and potentials are held exactly fixed because
        // the authored virtual displacement vanishes on every exterior node.
        let field = changed.solve(&boundary, controls(), &mut || false).unwrap();
        coenergy(&field)
    };
    let mut old = None;
    for h in [2e-3, 1e-3] {
        let derivative = (shifted(h) - shifted(-h)) / (2.0 * h);
        let error = (derivative - actual.force.x).abs();
        eprintln!(
            "eggshell coenergy h={h} derivative={derivative} force={} defect={error}",
            actual.force.x
        );
        assert!(error < 2e-5);
        if let Some(previous) = old {
            assert!(error < 0.45 * previous);
        }
        old = Some(error);
    }
}

#[test]
fn prepared_shell_reuses_geometry_and_respects_frames_origins_and_cancellation() {
    let problem = magnet(8);
    let field = problem
        .solve(&applied(&problem), controls(), &mut || false)
        .unwrap();
    let w = weights(problem.positions_m(), true);
    let shell = problem
        .prepare_eggshell(&w, force_limits(), &mut || false)
        .unwrap();
    let original = field
        .eggshell_wrench(&shell, &Se3::identity(), [0.0; 3], &mut || false)
        .unwrap();
    let repeat = field
        .eggshell_wrench(&shell, &Se3::identity(), [0.0; 3], &mut || false)
        .unwrap();
    assert_eq!(original, repeat);
    let pose = Se3::from_parts(
        So3::exp(So3Tangent {
            angular: Vec3::new(0.3, -0.2, 0.4),
        })
        .unwrap(),
        Vec3::new(4.0, -3.0, 2.0),
    )
    .unwrap();
    let origin = Vec3::new(0.4, -0.7, 0.2);
    let transformed = original.transform_by(&pose);
    let shift = origin.cross(transformed.force);
    let expected = Wrench::new(
        Vec3::new(
            transformed.torque.x - shift.x,
            transformed.torque.y - shift.y,
            transformed.torque.z - shift.z,
        ),
        transformed.force,
    );
    let actual = field
        .eggshell_wrench(&shell, &pose, [origin.x, origin.y, origin.z], &mut || false)
        .unwrap();
    assert_wrench(actual, expected, 2e-11);
    let other = magnet(8);
    let other_field = other
        .solve(&applied(&other), controls(), &mut || false)
        .unwrap();
    assert!(matches!(
        other_field.eggshell_wrench(&shell, &pose, [0.0; 3], &mut || false),
        Err(SpatialMagneticError::Invalid(_))
    ));
    let mut polls = 0;
    assert!(matches!(
        field.eggshell_wrench(&shell, &pose, [0.0; 3], &mut || {
            polls += 1;
            polls > shell.collar_cell_count() / 2
        }),
        Err(SpatialMagneticError::Cancelled)
    ));
    let after = field
        .eggshell_wrench(&shell, &Se3::identity(), [0.0; 3], &mut || false)
        .unwrap();
    assert_eq!(original, after);
}

#[test]
fn complete_closed_weights_original_source_identity_and_fractional_sources_are_admitted() {
    let problem = linear(4, 1.0);
    let w = weights(problem.positions_m(), false);
    let field = problem
        .solve(
            &potentials(&problem, |p| -2.0 * p[0] + p[1] - 3.0 * p[2]),
            controls(),
            &mut || false,
        )
        .unwrap();
    assert_wrench(load(&field, false), Wrench::default(), 1e-9);
    let mut bad = w.clone();
    bad[0] = 0.01;
    assert!(
        problem
            .prepare_eggshell(&bad, force_limits(), &mut || false)
            .is_err()
    );
    for value in [f64::NAN, -0.01, 1.01] {
        let mut bad = w.clone();
        bad[31] = value;
        assert!(
            problem
                .prepare_eggshell(&bad, force_limits(), &mut || false)
                .is_err()
        );
    }
    assert!(
        problem
            .prepare_eggshell(&w[..w.len() - 1], force_limits(), &mut || false)
            .is_err()
    );
    assert!(
        problem
            .prepare_eggshell(&vec![0.0; w.len()], force_limits(), &mut || false)
            .is_err()
    );
    let short = SpatialMaxwellForceLimits {
        max_cells: problem.mesh().tets.len() - 1,
        max_faces: 100000,
    };
    assert!(problem.prepare_eggshell(&w, short, &mut || false).is_err());
    assert!(matches!(
        problem.prepare_eggshell(&w, force_limits(), &mut || true),
        Err(SpatialMagneticError::Cancelled)
    ));

    let card = common::source(
        PropertyValue::Scalar {
            value: 1.0,
            dims: MAGNETIC_PERMEABILITY_DIMS,
        },
        ValidityDomain::unconstrained(),
    );
    let point = QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
        .unwrap();
    let distinct = common::resolve(&card, &point);
    let mut sources = problem.materials().to_vec();
    let cell = problem
        .mesh()
        .tets
        .iter()
        .position(|t| {
            t.iter().any(|&v| w[v as usize] == 0.0) && t.iter().any(|&v| w[v as usize] > 0.0)
        })
        .unwrap();
    sources[cell] = SpatialMagneticMaterial::Linear(distinct);
    let changed = SpatialMagnetostatics::new(
        problem.positions_m().to_vec(),
        problem.mesh().tets.clone(),
        sources,
        limits(),
        &mut || false,
    )
    .unwrap();
    let field = changed
        .solve(&potentials(&changed, |p| -p[0]), controls(), &mut || false)
        .unwrap();
    let shell = changed
        .prepare_eggshell(&w, force_limits(), &mut || false)
        .unwrap();
    assert!(matches!(
        field.eggshell_wrench(&shell, &Se3::identity(), [0.0; 3], &mut || false),
        Err(SpatialMagneticError::Invalid(_))
    ));

    let magnet = magnet(8);
    let mut partial = weights(magnet.positions_m(), false);
    for (p, w) in magnet.positions_m().iter().zip(&mut partial) {
        if p.iter().all(|x| x.abs() <= 0.5) {
            *w = 0.5;
        }
    }
    let shell = magnet
        .prepare_eggshell(&partial, force_limits(), &mut || false)
        .unwrap();
    let field = magnet
        .solve(&applied(&magnet), controls(), &mut || false)
        .unwrap();
    assert!(matches!(
        field.eggshell_wrench(&shell, &Se3::identity(), [0.0; 3], &mut || false),
        Err(SpatialMagneticError::Invalid(_))
    ));
}

#[test]
fn saturated_core_is_enclosed_by_the_same_source_admitted_volume_functional() {
    let q = [-2.0, -1.0, -0.5, 0.5, 1.0, 2.0];
    let (p, t) = grid(&q, &q, &q);
    let law = saturation_law();
    let air = material(1.0);
    let materials: Vec<_> = t
        .iter()
        .map(|&tet| {
            if center(&p, tet).iter().all(|x| x.abs() < 0.5) {
                SpatialNonlinearMaterial::Reversible(law.clone())
            } else {
                SpatialNonlinearMaterial::Linear(air.clone())
            }
        })
        .collect();
    let problem =
        SpatialNonlinearMagnetostatics::new(p, t, materials, limits(), &mut || false).unwrap();
    let initial: Vec<_> = problem
        .positions_m()
        .iter()
        .map(|p| -p.iter().sum::<f64>())
        .collect();
    let boundary = SpatialMagneticBoundary {
        potential: problem
            .exterior_faces()
            .into_iter()
            .map(|face| MagneticPotentialFace {
                face,
                potential_a: problem.mesh().faces[face].map(|v| initial[v as usize]),
            })
            .collect(),
        flux: vec![],
    };
    let settings = SpatialNonlinearControls {
        field: controls(),
        relative_residual_tolerance: 1e-10,
        maximum_newton_iterations: 35,
        maximum_backtracks: 35,
    };
    let field = problem
        .solve(&boundary, &initial, settings, &mut || false)
        .unwrap();
    assert!(
        problem
            .materials()
            .iter()
            .enumerate()
            .filter(|(_, m)| matches!(m, SpatialNonlinearMaterial::Reversible(_)))
            .all(|(i, _)| {
                let b = field.flux_density_t()[i];
                b[0].hypot(b[1]).hypot(b[2]) > 1.0
            })
    );
    let w = weights(problem.positions_m(), false);
    let shell = problem
        .prepare_eggshell(&w, force_limits(), &mut || false)
        .unwrap();
    let force = field
        .eggshell_wrench(&shell, &Se3::identity(), [0.0; 3], &mut || false)
        .unwrap();
    assert_wrench(force, Wrench::default(), 3e-8);
    let mut bad = w.clone();
    let inner = problem
        .positions_m()
        .iter()
        .position(|p| p.iter().all(|x| x.abs() == 0.5))
        .unwrap();
    bad[inner] = 0.9;
    let shell = problem
        .prepare_eggshell(&bad, force_limits(), &mut || false)
        .unwrap();
    assert!(matches!(
        field.eggshell_wrench(&shell, &Se3::identity(), [0.0; 3], &mut || false),
        Err(SpatialMagneticError::Invalid(_))
    ));
}

fn vector_controls() -> VectorMagneticControls {
    VectorMagneticControls {
        field: fs_solver::CheckedCgConfig {
            tolerance: 1e-11,
            max_iterations: 2500,
            max_corrections: 2,
        },
        gauge: fs_solver::CheckedCgConfig {
            tolerance: 1e-12,
            max_iterations: 2500,
            max_corrections: 2,
        },
        maximum_gauge_product_terms: 2_000_000,
        current_tolerance_a: 2e-8,
        flux_tolerance_wb: 2e-9,
        energy_tolerance_j: 3e-8,
    }
}
fn vector_trace(
    problem: &SpatialVectorMagnetostatics,
    potential: impl Fn([f64; 3]) -> [f64; 3],
) -> VectorMagneticBoundary {
    let integrals: Vec<_> = problem
        .mesh()
        .edges
        .iter()
        .map(|edge| {
            let a = problem.positions_m()[edge[0] as usize];
            let b = problem.positions_m()[edge[1] as usize];
            let m = std::array::from_fn(|i| f64::midpoint(a[i], b[i]));
            let (aa, am, ab) = (potential(a), potential(m), potential(b));
            (0..3)
                .map(|i| (b[i] - a[i]) * (aa[i] + 4.0 * am[i] + ab[i]) / 6.0)
                .sum::<f64>()
        })
        .collect();
    VectorMagneticBoundary {
        potential: problem
            .exterior_faces()
            .into_iter()
            .map(|face| VectorPotentialFace {
                face,
                edge_integrals_wb: problem.face_edges(face).unwrap().map(|e| integrals[e]),
            })
            .collect(),
    }
}
#[test]
fn vector_field_reuses_eggshell_and_refuses_current_through_its_transition() {
    use fs_em::spatial::vector_potential::nonlinear::{
        NonlinearVectorControls, SpatialNonlinearVectorMagnetostatics,
    };
    let base = linear(4, 1.0);
    let problem = SpatialVectorMagnetostatics::new(
        base.positions_m().to_vec(),
        base.mesh().tets.clone(),
        base.materials().to_vec(),
        limits(),
        &mut || false,
    )
    .unwrap();
    let nonlinear = SpatialNonlinearVectorMagnetostatics::new(
        base.positions_m().to_vec(),
        base.mesh().tets.clone(),
        base.materials()
            .iter()
            .map(|m| match m {
                SpatialMagneticMaterial::Linear(m) => SpatialNonlinearMaterial::Linear(m.clone()),
                SpatialMagneticMaterial::Recoil(_) => unreachable!(),
            })
            .collect(),
        limits(),
        &mut || false,
    )
    .unwrap();
    let w = weights(problem.positions_m(), false);
    let shell = problem
        .prepare_eggshell(&w, force_limits(), &mut || false)
        .unwrap();
    let nonlinear_shell = nonlinear
        .prepare_eggshell(&w, force_limits(), &mut || false)
        .unwrap();
    let settings = NonlinearVectorControls {
        field: vector_controls(),
        relative_residual_tolerance: 1e-10,
        maximum_newton_iterations: 35,
        maximum_backtracks: 35,
    };
    for current in [false, true] {
        let boundary = if current {
            vector_trace(&problem, |p| {
                [0.0, 0.0, -0.25 * (p[0] * p[0] + p[1] * p[1])]
            })
        } else {
            vector_trace(&problem, |p| [-0.5 * p[1], 0.5 * p[0], 0.0])
        };
        let currents = vec![[0.0, 0.0, if current { 1.0 } else { 0.0 }]; problem.mesh().tets.len()];
        let field = problem
            .solve(&currents, &boundary, vector_controls(), &mut || false)
            .unwrap();
        let nonlinear_field = nonlinear
            .solve(
                &currents,
                &boundary,
                field.potential_wb(),
                settings,
                &mut || false,
            )
            .unwrap();
        let force = field.eggshell_wrench(&shell, &Se3::identity(), [0.0; 3], &mut || false);
        let nonlinear_force = nonlinear_field.eggshell_wrench(
            &nonlinear_shell,
            &Se3::identity(),
            [0.0; 3],
            &mut || false,
        );
        if current {
            assert!(matches!(force, Err(SpatialMagneticError::Invalid(_))));
            assert!(matches!(
                nonlinear_force,
                Err(SpatialMagneticError::Invalid(_))
            ));
        } else {
            assert_wrench(force.unwrap(), Wrench::default(), 2e-8);
            assert_wrench(nonlinear_force.unwrap(), Wrench::default(), 2e-8);
        }
    }
}
