//! G1/G3 prescribed eigenstrain through the canonical plane CutFEM owner.

use fs_cutfem::elastic::eigenstrain::Q1Eigenstrain;
use fs_cutfem::{
    Circle, CutElasticity, CutElasticityOperator, CutFemError, CutSdf, CutStabilizationScaling,
    ElasticPlane, HalfPlane, Quadtree,
};
use fs_material::{IsotropicElastic, Voigt};
use std::ops::ControlFlow;

fn proceed() -> ControlFlow<()> {
    ControlFlow::Continue(())
}
fn zero(_: f64, _: f64) -> [f64; 2] {
    [0.0; 2]
}
fn problem<'a>(
    grid: &'a Quadtree,
    sdf: &'a dyn CutSdf,
    material: &'a IsotropicElastic,
) -> CutElasticity<'a> {
    CutElasticity {
        grid,
        sdf,
        material,
        nitsche_beta: 80.0,
        ghost_gamma: 0.1,
        stabilization_scaling: CutStabilizationScaling::LongitudinalModulus,
        quad_depth: 4,
        clamp: None,
        boundary_traction: None,
        traction_free_interface: false,
        solver_tol: 1e-11,
        solver_max_iters: 20_000,
    }
}
fn assemble(
    cut: &CutElasticity<'_>,
    plane: ElasticPlane,
    body: &dyn Fn(f64, f64) -> [f64; 2],
    boundary: &dyn Fn(f64, f64) -> [f64; 2],
) -> CutElasticityOperator {
    match plane {
        ElasticPlane::Strain => cut.assemble(body, boundary).unwrap(),
        ElasticPlane::Stress => cut.assemble_plane_stress(body, boundary).unwrap(),
    }
}
fn continued<T>(value: ControlFlow<(), T>) -> T {
    match value {
        ControlFlow::Continue(value) => value,
        ControlFlow::Break(()) => panic!("unexpected cancellation"),
    }
}
fn close(actual: f64, expected: f64, scale: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * scale.max(actual.abs()).max(expected.abs()),
        "actual={actual:e}, expected={expected:e}"
    );
}
fn coefficients(
    grid: &Quadtree,
    operator: &CutElasticityOperator,
    exact: impl Fn(f64, f64) -> [f64; 2],
) -> Vec<f64> {
    let mut result = vec![0.0; operator.dof_count()];
    for (&node, &id) in operator.node_ids() {
        let [x, y] = grid.node_pos(node);
        let value = exact(x, y);
        result[2 * id..2 * id + 2].copy_from_slice(&value);
    }
    result
}
fn assert_equilibrium(operator: &CutElasticityOperator, exact: &[f64], force_scale: f64) {
    for (applied, rhs) in operator.apply_vec(exact).iter().zip(operator.rhs()) {
        close(*applied, *rhs, force_scale, 3e-12);
    }
}

#[test]
fn g1_isotropic_expansion_retains_free_thickness_and_plane_strain_stress() {
    let expansion = 0.001;
    for graded in [false, true] {
        let mut grid = Quadtree::with_room(2, 3);
        if graded {
            grid.refine_where(3, &|lo, _| lo[0] >= 0.5);
        }
        let sdf = Circle {
            center: [0.5, 0.5],
            radius: 0.37,
        };
        for (plane, nu) in [
            (ElasticPlane::Strain, 0.3),
            (ElasticPlane::Stress, 0.3),
            (ElasticPlane::Stress, 0.499_999_9),
        ] {
            let material = IsotropicElastic::new(7.3e9, nu, 0.02).unwrap();
            let cut = problem(&grid, &sdf, &material);
            // A laterally free plane-strain layer must expand by (1+nu)*e:
            // sigma_xx=sigma_yy=0 while sigma_zz=-E*e remains constrained.
            let planar = if plane == ElasticPlane::Strain {
                (1.0 + nu) * expansion
            } else {
                expansion
            };
            let exact = |x: f64, y: f64| [planar * (x - 0.5), planar * (y - 0.5)];
            let original = assemble(&cut, plane, &zero, &exact);
            let matrix = original.matrix().clone();
            let operator = continued(
                original
                    .with_eigenstrain(
                        0.1,
                        |_, _| Q1Eigenstrain::try_isotropic([expansion; 4]),
                        proceed,
                    )
                    .unwrap(),
            );
            assert_eq!(*operator.matrix(), matrix);
            let values = coefficients(&grid, &operator, exact);
            assert_equilibrium(&operator, &values, material.youngs * expansion);
            let solution = continued(
                operator
                    .solve_controlled(1e-11, 20_000, 23, |_| proceed())
                    .unwrap(),
            );
            for (&actual, &expected) in solution.coefficients().iter().zip(&values) {
                close(actual, expected, expansion, 3e-9);
            }
            for &cell in operator.active_cells() {
                let (lo, hi) = grid.rect(cell);
                let point = [f64::midpoint(lo[0], hi[0]), f64::midpoint(lo[1], hi[1])];
                let response = operator
                    .strain_stress(&grid, solution.coefficients(), cell, point)
                    .unwrap();
                for i in [0, 1, 3, 4, 5] {
                    close(response.stress[i], 0.0, material.youngs * expansion, 3e-8);
                }
                let expected_zz = if plane == ElasticPlane::Strain {
                    -material.youngs * expansion
                } else {
                    0.0
                };
                close(
                    response.stress[2],
                    expected_zz,
                    material.youngs * expansion,
                    3e-8,
                );
                close(
                    response.strain[2],
                    if plane == ElasticPlane::Strain {
                        0.0
                    } else {
                        expansion
                    },
                    expansion,
                    3e-8,
                );
                let energy = if plane == ElasticPlane::Strain {
                    0.5 * material.youngs * expansion * expansion
                } else {
                    0.0
                };
                close(
                    response.energy_density,
                    energy,
                    material.youngs * expansion * expansion,
                    3e-8,
                );
                assert_eq!(
                    response.free_strain,
                    [expansion, expansion, expansion, 0.0, 0.0, 0.0]
                );
            }
        }
    }
}

#[test]
fn g1_restrained_expansion_has_correct_nitsche_sign_and_compressive_storage() {
    let grid = Quadtree::uniform(2);
    let sdf = HalfPlane {
        normal: [1.0, 0.0],
        offset: 0.83,
    };
    let material = IsotropicElastic::new(210e9, 0.3, 0.02).unwrap();
    let expansion = 0.001;
    for plane in [ElasticPlane::Strain, ElasticPlane::Stress] {
        let mut cut = problem(&grid, &sdf, &material);
        cut.clamp = Some(&|_, _| true);
        let original = assemble(&cut, plane, &zero, &zero);
        let operator = continued(
            original
                .with_eigenstrain(
                    0.1,
                    |_, _| Q1Eigenstrain::try_isotropic([expansion; 4]),
                    proceed,
                )
                .unwrap(),
        );
        let values = vec![0.0; operator.dof_count()];
        // The volume term and negative initial-traction term cancel for the
        // boundary-connected free test functions. A missing/wrong surface term
        // would produce O(E*expansion) spurious displacement at the cut face.
        assert_equilibrium(&operator, &values, material.youngs * expansion);
        for &cell in operator.active_cells() {
            let (lo, _) = grid.rect(cell);
            let response = operator.strain_stress(&grid, &values, cell, lo).unwrap();
            let denominator = if plane == ElasticPlane::Strain {
                1.0 - 2.0 * material.poisson
            } else {
                1.0 - material.poisson
            };
            let pressure = -material.youngs * expansion / denominator;
            close(response.stress[0], pressure, pressure.abs(), 2e-14);
            close(response.stress[1], pressure, pressure.abs(), 2e-14);
            close(
                response.stress[2],
                if plane == ElasticPlane::Strain {
                    pressure
                } else {
                    0.0
                },
                pressure.abs(),
                2e-14,
            );
            let expected_energy = if plane == ElasticPlane::Strain {
                -1.5 * pressure * expansion
            } else {
                -pressure * expansion
            };
            close(
                response.energy_density,
                expected_energy,
                expected_energy,
                2e-14,
            );
            let expected_thickness = if plane == ElasticPlane::Stress {
                expansion * (1.0 + material.poisson) / (1.0 - material.poisson)
            } else {
                0.0
            };
            close(response.strain[2], expected_thickness, expansion, 2e-14);
        }
    }
}

#[test]
fn g1_compatible_spatial_eigenstrain_uses_native_hanging_node_and_cut_operators() {
    let mut grid = Quadtree::with_room(2, 3);
    grid.refine_where(3, &|lo, _| lo[0] >= 0.5);
    let sdf = Circle {
        center: [0.5, 0.5],
        radius: 0.39,
    };
    let material = IsotropicElastic::new(12e6, 0.25, 0.04).unwrap();
    let (a, b) = (0.009, -0.004);
    let exact = |x: f64, y: f64| [a * x * y, b * x * y];
    for plane in [ElasticPlane::Strain, ElasticPlane::Stress] {
        let cut = problem(&grid, &sdf, &material);
        let original = assemble(&cut, plane, &zero, &exact);
        let operator = continued(
            original
                .with_eigenstrain(
                    0.1,
                    |_, corners| {
                        Q1Eigenstrain::try_new(corners.map(|[x, y]| {
                            [
                                a * y,
                                b * x,
                                if plane == ElasticPlane::Stress {
                                    0.002 * x * y
                                } else {
                                    0.0
                                },
                                0.5 * (a * x + b * y),
                                if plane == ElasticPlane::Stress {
                                    0.001 * x
                                } else {
                                    0.0
                                },
                                if plane == ElasticPlane::Stress {
                                    -0.001 * y
                                } else {
                                    0.0
                                },
                            ]
                        }))
                    },
                    proceed,
                )
                .unwrap(),
        );
        let values = coefficients(&grid, &operator, exact);
        assert_equilibrium(&operator, &values, material.youngs * a);
        let solution = continued(
            operator
                .solve_controlled(1e-11, 20_000, 7, |_| proceed())
                .unwrap(),
        );
        for (actual, expected) in solution.coefficients().iter().zip(&values) {
            close(*actual, *expected, a, 3e-8);
        }
        for &cell in operator.active_cells() {
            let (lo, hi) = grid.rect(cell);
            let point = [0.3 * lo[0] + 0.7 * hi[0], 0.6 * lo[1] + 0.4 * hi[1]];
            let response = operator.strain_stress(&grid, &values, cell, point).unwrap();
            for stress in response.stress {
                close(stress, 0.0, material.youngs * a, 2e-13);
            }
            close(response.energy_density, 0.0, material.youngs * a * a, 2e-25);
            if plane == ElasticPlane::Stress {
                close(response.strain[2], 0.002 * point[0] * point[1], a, 2e-14);
                close(response.strain[4], 0.001 * point[0], a, 2e-14);
                close(response.strain[5], -0.001 * point[1], a, 2e-14);
            }
        }
    }
}

#[test]
fn g1_nonuniform_thermal_stress_balances_the_original_body_force() {
    let grid = Quadtree::uniform(2);
    let sdf = HalfPlane {
        normal: [1.0, 0.0],
        offset: 0.83,
    };
    let material = IsotropicElastic::new(1.1e8, 0.2, 0.03).unwrap();
    let slope = 0.002;
    for plane in [ElasticPlane::Strain, ElasticPlane::Stress] {
        let mut cut = problem(&grid, &sdf, &material);
        cut.clamp = Some(&|_, _| true);
        let denominator = if plane == ElasticPlane::Strain {
            1.0 - 2.0 * material.poisson
        } else {
            1.0 - material.poisson
        };
        let thermal_modulus = material.youngs / denominator;
        let body = |_: f64, _: f64| [thermal_modulus * slope, 0.0];
        let operator = continued(
            assemble(&cut, plane, &body, &zero)
                .with_eigenstrain(
                    0.1,
                    |_, corners| Q1Eigenstrain::try_isotropic(corners.map(|[x, _]| slope * x)),
                    proceed,
                )
                .unwrap(),
        );
        let values = vec![0.0; operator.dof_count()];
        assert_equilibrium(&operator, &values, thermal_modulus * slope);
        for &cell in operator.active_cells() {
            let (lo, hi) = grid.rect(cell);
            let point = [f64::midpoint(lo[0], hi[0]), f64::midpoint(lo[1], hi[1])];
            let response = operator.strain_stress(&grid, &values, cell, point).unwrap();
            close(
                response.stress[0],
                -thermal_modulus * slope * point[0],
                thermal_modulus * slope,
                2e-14,
            );
            close(
                response.stress[1],
                response.stress[0],
                thermal_modulus * slope,
                2e-14,
            );
        }
    }
}

#[test]
fn g1_natural_cut_surface_releases_directional_initial_strain() {
    let mut grid = Quadtree::with_room(2, 3);
    grid.refine_where(3, &|lo, _| lo[1] >= 0.5);
    let sdf = HalfPlane {
        normal: [1.0, 0.0],
        offset: 0.83,
    };
    let material = IsotropicElastic::new(8e8, 0.3, 0.02).unwrap();
    let free: Voigt = [0.003, 0.0, 0.0, 0.0, 0.0, 0.0];
    let mut cut = problem(&grid, &sdf, &material);
    cut.traction_free_interface = true;
    cut.clamp = Some(&|x, _| x == 0.0);
    for plane in [ElasticPlane::Strain, ElasticPlane::Stress] {
        let operator = continued(
            assemble(&cut, plane, &zero, &zero)
                .with_eigenstrain(0.1, |_, _| Q1Eigenstrain::try_uniform(free), proceed)
                .unwrap(),
        );
        let exact = coefficients(&grid, &operator, |x, _| [free[0] * x, 0.0]);
        assert_equilibrium(&operator, &exact, material.youngs * free[0]);
        let solved = continued(
            operator
                .solve_controlled(1e-11, 20_000, 11, |_| proceed())
                .unwrap(),
        );
        for (value, expected) in solved.coefficients().iter().zip(exact) {
            close(*value, expected, free[0], 5e-8);
        }
    }
}

#[test]
fn g0_eigenstrain_is_retained_once_and_failures_publish_no_partial_field() {
    let grid = Quadtree::uniform(2);
    let sdf = HalfPlane {
        normal: [1.0, 0.0],
        offset: 0.83,
    };
    let material = IsotropicElastic::new(10.0, 0.3, 0.1).unwrap();
    let original = assemble(
        &problem(&grid, &sdf, &material),
        ElasticPlane::Stress,
        &zero,
        &zero,
    );
    let mut calls = 0;
    let operator = continued(
        original
            .clone()
            .with_eigenstrain(
                0.1,
                |_, _| {
                    calls += 1;
                    Q1Eigenstrain::try_uniform([0.0; 6])
                },
                proceed,
            )
            .unwrap(),
    );
    assert_eq!(calls, operator.active_cells().len());
    assert_eq!(*operator.matrix(), *original.matrix());
    assert_eq!(operator.rhs(), original.rhs());
    assert!(
        operator
            .clone()
            .with_eigenstrain(0.1, |_, _| Q1Eigenstrain::try_uniform([0.0; 6]), proceed)
            .is_err()
    );
    assert!(Q1Eigenstrain::try_uniform([f64::NAN; 6]).is_err());
    let mut polls = 0;
    let stopped = original
        .clone()
        .with_eigenstrain(
            0.1,
            |_, _| Q1Eigenstrain::try_uniform([0.0; 6]),
            || {
                polls += 1;
                if polls == 4 {
                    ControlFlow::Break("stop")
                } else {
                    ControlFlow::Continue(())
                }
            },
        )
        .unwrap();
    assert!(matches!(stopped, ControlFlow::Break("stop")));
    assert!(
        original
            .active_cells()
            .iter()
            .all(|&cell| original.eigenstrain_cell(cell).is_none())
    );
    let refused = original.clone().with_eigenstrain(
        0.1,
        |_, _| {
            Err(CutFemError::InvalidElasticityInput {
                what: "source temperature left its measured domain".into(),
            })
        },
        proceed,
    );
    assert!(
        matches!(refused, Err(CutFemError::InvalidElasticityInput { what }) if what.contains("source temperature"))
    );
}

#[test]
fn g0_both_total_and_elastic_vertex_envelopes_remain_binding() {
    let grid = Quadtree::uniform(1);
    let sdf = HalfPlane {
        normal: [1.0, 0.0],
        offset: 0.83,
    };
    let material = IsotropicElastic::new(10.0, 0.3, 0.01).unwrap();
    let cut = problem(&grid, &sdf, &material);
    for plane in [ElasticPlane::Strain, ElasticPlane::Stress] {
        let operator = continued(
            assemble(&cut, plane, &zero, &zero)
                .with_eigenstrain(
                    0.01,
                    |_, _| {
                        let mut corners = [[0.0; 6]; 4];
                        corners[0][2] = 0.02;
                        Q1Eigenstrain::try_new(corners)
                    },
                    proceed,
                )
                .unwrap(),
        );
        let values = vec![0.0; operator.dof_count()];
        let cell = operator.active_cells()[0];
        let (lo, hi) = grid.rect(cell);
        let midpoint = [f64::midpoint(lo[0], hi[0]), f64::midpoint(lo[1], hi[1])];
        assert!(
            operator
                .strain_stress(&grid, &values, cell, midpoint)
                .is_ok()
        );
        assert!(
            operator
                .admit_strain_domain(&grid, &values, proceed)
                .is_err()
        );
        // In plane strain the elastic transverse strain is excessive. In plane
        // stress it relaxes exactly, but TOTAL thickness strain still exceeds
        // the original small-strain kinematic domain at the same corner.
        assert!(operator.strain_stress(&grid, &values, cell, lo).is_err());
    }
}

#[test]
fn g1_free_expansion_uses_a_separate_kinematic_budget_without_relaxing_the_elastic_domain() {
    let grid = Quadtree::uniform(2);
    let sdf = Circle {
        center: [0.5, 0.5],
        radius: 0.37,
    };
    let material = IsotropicElastic::new(7.3e9, 0.3, 0.001).unwrap();
    let expansion = 0.01;
    let exact = |x: f64, y: f64| [expansion * (x - 0.5), expansion * (y - 0.5)];
    let cut = problem(&grid, &sdf, &material);
    let original = assemble(&cut, ElasticPlane::Stress, &zero, &exact);
    let free = |_, _| Q1Eigenstrain::try_isotropic([expansion; 4]);
    let operator = continued(
        original
            .clone()
            .with_eigenstrain(0.04, free, proceed)
            .unwrap(),
    );
    let solution = continued(
        operator
            .solve_controlled(1e-11, 20_000, 13, |_| proceed())
            .unwrap(),
    );
    let cell = operator.active_cells()[0];
    let (lo, hi) = grid.rect(cell);
    let point = [f64::midpoint(lo[0], hi[0]), f64::midpoint(lo[1], hi[1])];
    let response = operator
        .strain_stress(&grid, solution.coefficients(), cell, point)
        .unwrap();
    assert!(response.strain_magnitude > material.strain_limit);
    assert!(response.elastic_strain_magnitude < 1e-8 * material.strain_limit);
    close(
        response.strain_magnitude,
        3.0_f64.sqrt() * expansion,
        expansion,
        3e-8,
    );
    assert_eq!(operator.maximum_kinematic_strain(), 0.04);

    let limited = continued(
        original
            .clone()
            .with_eigenstrain(0.01, free, proceed)
            .unwrap(),
    );
    assert!(
        limited
            .strain_stress(&grid, solution.coefficients(), cell, point)
            .is_err()
    );
    assert!(
        limited
            .solve_controlled(1e-11, 20_000, 13, |_| proceed())
            .is_err()
    );

    // Blocking the same thermal expansion produces elastic strain far beyond
    // its original 0.001 calibration, even though total strain is small. The
    // separate kinematic budget must not weaken that constitutive refusal.
    let constrained = vec![0.0; operator.dof_count()];
    assert!(
        operator
            .strain_stress(&grid, &constrained, cell, point)
            .is_err()
    );
    assert!(
        operator
            .admit_strain_domain(&grid, &constrained, proceed)
            .is_err()
    );
    for invalid in [0.0, -0.1, 1.0, f64::NAN, f64::INFINITY] {
        assert!(
            original
                .clone()
                .with_eigenstrain(invalid, free, proceed)
                .is_err()
        );
    }
}
