//! G0/G1/G3 physical strain recovery on the original CutFEM trial space.

use fs_cutfem::elastic::ElasticPlane;
use fs_cutfem::{
    BoundaryTraction, CutElasticity, CutElasticityOperator, CutStabilizationScaling,
    DesignBoxEdge, EdgeBand, HalfPlane, Quadtree,
};
use fs_material::IsotropicElastic;
use std::ops::ControlFlow;

fn problem<'a>(
    grid: &'a Quadtree,
    sdf: &'a HalfPlane,
    material: &'a IsotropicElastic,
) -> CutElasticity<'a> {
    CutElasticity {
        grid,
        sdf,
        material,
        nitsche_beta: 20.0,
        ghost_gamma: 0.5,
        stabilization_scaling: CutStabilizationScaling::MuScaled,
        quad_depth: 2,
        clamp: None,
        boundary_traction: None,
        traction_free_interface: true,
        solver_tol: 1e-12,
        solver_max_iters: 10_000,
    }
}

fn coefficients(
    grid: &Quadtree,
    op: &CutElasticityOperator,
    field: impl Fn([f64; 2]) -> [f64; 2],
) -> Vec<f64> {
    let mut values = vec![0.0; op.dof_count()];
    for (&node, &id) in op.node_ids() {
        let u = field(grid.node_pos(node));
        values[2 * id..2 * id + 2].copy_from_slice(&u);
    }
    values
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 3e-12 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn affine_full_tensors_match_independent_plane_laws_and_assembled_work() {
    for graded in [false, true] {
        let mut grid = Quadtree::with_room(1, 3);
        if graded {
            grid.refine_where(3, &|lo, _| lo[0] < 0.5);
        }
        let sdf = HalfPlane {
            normal: [1.0, 0.0],
            offset: 2.0,
        };
        for nu in [-0.4, 0.3, 0.499_999_999_999] {
            let material = IsotropicElastic::new(12.0, nu, 0.2).unwrap();
            let cut = problem(&grid, &sdf, &material);
            for plane in [ElasticPlane::Strain, ElasticPlane::Stress] {
                if plane == ElasticPlane::Strain && nu > 0.4 {
                    continue;
                }
                let op = match plane {
                    ElasticPlane::Strain => cut.assemble(&|_, _| [0.0; 2], &|_, _| [0.0; 2]),
                    ElasticPlane::Stress => {
                        cut.assemble_plane_stress(&|_, _| [0.0; 2], &|_, _| [0.0; 2])
                    }
                }
                .unwrap();
                let u = coefficients(&grid, &op, |[x, y]| {
                    [0.02 * x + 0.012 * y, 0.004 * x - 0.006 * y]
                });
                let mu = material.youngs / (2.0 * (1.0 + nu));
                let lambda = match plane {
                    ElasticPlane::Strain => material.youngs * nu / ((1.0 + nu) * (1.0 - 2.0 * nu)),
                    ElasticPlane::Stress => material.youngs * nu / (1.0 - nu * nu),
                };
                let zz = if plane == ElasticPlane::Stress {
                    -nu / (1.0 - nu) * 0.014
                } else {
                    0.0
                };
                let expected_strain = [0.02, -0.006, zz, 0.008, 0.0, 0.0];
                let expected_stress = [
                    2.0 * mu * 0.02 + lambda * 0.014,
                    -2.0 * mu * 0.006 + lambda * 0.014,
                    if plane == ElasticPlane::Strain {
                        lambda * 0.014
                    } else {
                        0.0
                    },
                    mu * 0.016,
                    0.0,
                    0.0,
                ];
                let expected_energy = 0.5
                    * (expected_strain[0] * expected_stress[0]
                        + expected_strain[1] * expected_stress[1]
                        + 2.0 * expected_strain[3] * expected_stress[3]);
                let expected_norm =
                    (0.02_f64.powi(2) + 0.006_f64.powi(2) + zz * zz + 2.0 * 0.008_f64.powi(2))
                        .sqrt();
                for &cell in op.active_cells() {
                    let (lo, hi) = grid.rect(cell);
                    let sample = op
                        .strain_stress(
                            &grid,
                            &u,
                            cell,
                            [0.37 * lo[0] + 0.63 * hi[0], 0.61 * lo[1] + 0.39 * hi[1]],
                        )
                        .unwrap();
                    for i in 0..6 {
                        close(sample.strain[i], expected_strain[i]);
                        close(sample.stress[i], expected_stress[i]);
                    }
                    close(sample.energy_density, expected_energy);
                    close(sample.strain_magnitude, expected_norm);
                }
                let ControlFlow::Continue(bound) = op
                    .admit_strain_domain(&grid, &u, || ControlFlow::<()>::Continue(()))
                    .unwrap()
                else {
                    panic!("cancelled")
                };
                close(bound, expected_norm);
                let ku = op.apply_vec(&u);
                let energy: f64 = 0.5 * u.iter().zip(ku).map(|(a, b)| a * b).sum::<f64>();
                close(energy, expected_energy);
            }
        }
    }
}

#[test]
fn rigid_translation_and_infinitesimal_rotation_have_no_strain() {
    let grid = Quadtree::uniform(2);
    let sdf = HalfPlane {
        normal: [1.0, 0.0],
        offset: 2.0,
    };
    let material = IsotropicElastic::new(1e9, 0.3, 0.01).unwrap();
    let op = problem(&grid, &sdf, &material)
        .assemble(&|_, _| [0.0; 2], &|_, _| [0.0; 2])
        .unwrap();
    for rotation in [0.0, 1.0 / 32.0] {
        let u = coefficients(&grid, &op, |[x, y]| {
            [1_048_576.0 - rotation * y, -2_097_152.0 + rotation * x]
        });
        for &cell in op.active_cells() {
            let (lo, hi) = grid.rect(cell);
            let sample = op
                .strain_stress(
                    &grid,
                    &u,
                    cell,
                    [
                        lo[0] + 0.237 * (hi[0] - lo[0]),
                        lo[1] + 0.731 * (hi[1] - lo[1]),
                    ],
                )
                .unwrap();
            assert!(sample.strain_magnitude < 1e-16);
            assert!(sample.energy_density < 1e-22);
        }
    }
}

#[test]
fn whole_cell_gate_finds_a_strain_peak_missed_by_gauss_sites() {
    let grid = Quadtree::uniform(0);
    let sdf = HalfPlane {
        normal: [1.0, 0.0],
        offset: 2.0,
    };
    let material = IsotropicElastic::new(10.0, 0.0, 0.11).unwrap();
    let op = problem(&grid, &sdf, &material)
        .assemble(&|_, _| [0.0; 2], &|_, _| [0.0; 2])
        .unwrap();
    let u = coefficients(&grid, &op, |[x, y]| [0.1 * x * y, 0.0]);
    let a = 0.5 / 3.0_f64.sqrt();
    for x in [0.5 - a, 0.5 + a] {
        for y in [0.5 - a, 0.5 + a] {
            assert!(op.strain_stress(&grid, &u, (0, 0, 0), [x, y]).is_ok());
        }
    }
    assert!(
        op.admit_strain_domain(&grid, &u, || ControlFlow::<()>::Continue(()))
            .is_err()
    );
    assert!(op.strain_stress(&grid, &u, (0, 0, 0), [1.0, 1.0]).is_err());
}

#[test]
fn full_tensor_domain_includes_relaxed_thickness_and_both_shear_entries() {
    let grid = Quadtree::uniform(0);
    let sdf = HalfPlane {
        normal: [1.0, 0.0],
        offset: 2.0,
    };
    let material = IsotropicElastic::new(10.0, 0.49, 0.08).unwrap();
    let op = problem(&grid, &sdf, &material)
        .assemble_plane_stress(&|_, _| [0.0; 2], &|_, _| [0.0; 2])
        .unwrap();
    // In-plane norm = .0707 < .08, but thickness contraction makes the full
    // material norm .1193. Checking only the visible components would pass.
    let u = coefficients(&grid, &op, |[x, y]| [0.05 * x, 0.05 * y]);
    assert!(op.strain_stress(&grid, &u, (0, 0, 0), [0.5; 2]).is_err());
    // Each tensor shear component is .06 < .08, but its two symmetric
    // entries give sqrt(2)*.06 > .08. Engineering shear is .12.
    let u = coefficients(&grid, &op, |[_, y]| [0.12 * y, 0.0]);
    assert!(op.strain_stress(&grid, &u, (0, 0, 0), [0.5; 2]).is_err());
}

#[test]
fn solved_and_controlled_fields_recover_the_original_material_response() {
    let grid = Quadtree::uniform(2);
    let sdf = HalfPlane {
        normal: [1.0, 0.0],
        offset: 0.83,
    };
    let material = IsotropicElastic::new(12.0, 0.3, 0.2).unwrap();
    let mut cut = problem(&grid, &sdf, &material);
    cut.traction_free_interface = false;
    cut.solver_tol = 1e-13;
    // Uniaxial tension: the embedded right edge prescribes the affine field,
    // the left edge carries its matching traction, and top/bottom are free.
    // Prescribing an arbitrary affine field only on the embedded edge would
    // not make that field the solution on the remaining traction-free edges.
    let extension = 0.01;
    let exact = |x: f64, y: f64| [extension * x, -material.poisson * extension * y];
    let traction = |_: f64, _: f64| [-material.youngs * extension, 0.0];
    let boundary = BoundaryTraction::EdgeBand {
        support: EdgeBand::new(DesignBoxEdge::Left, 0.0, 1.0).unwrap(),
        value: &traction,
    };
    let solution = cut
        .solve_plane_stress_with_boundary_traction(&|_, _| [0.0; 2], &exact, boundary)
        .unwrap();
    let op = cut
        .assemble_plane_stress_with_boundary_traction(&|_, _| [0.0; 2], &exact, boundary)
        .unwrap();
    let ControlFlow::Continue(controlled) = op
        .solve_controlled(cut.solver_tol, cut.solver_max_iters, 7, |_| {
            ControlFlow::<()>::Continue(())
        })
        .unwrap()
    else {
        panic!("cancelled")
    };
    let expected_strain = [
        extension,
        -material.poisson * extension,
        -material.poisson * extension,
        0.0,
        0.0,
        0.0,
    ];
    let expected_stress = [material.youngs * extension, 0.0, 0.0, 0.0, 0.0, 0.0];
    for &cell in solution.active_cells() {
        let (lo, hi) = grid.rect(cell);
        let point = [
            f64::midpoint(lo[0], hi[0].min(sdf.offset)),
            f64::midpoint(lo[1], hi[1]),
        ];
        let observed = solution.strain_stress(&grid, cell, point).unwrap();
        let other = op
            .strain_stress(&grid, controlled.coefficients(), cell, point)
            .unwrap();
        for i in 0..6 {
            close(observed.strain[i], other.strain[i]);
            close(observed.strain[i], expected_strain[i]);
            close(observed.stress[i], other.stress[i]);
            close(observed.stress[i], expected_stress[i]);
        }
        close(observed.energy_density, other.energy_density);
        close(observed.energy_density, 0.5 * material.youngs * extension * extension);
        close(observed.stress[2], 0.0);
    }
    let wrong_lattice = Quadtree::uniform(3);
    assert!(
        solution
            .strain_stress(&wrong_lattice, solution.active_cells()[0], [0.0; 2])
            .is_err()
    );
}

#[test]
fn observation_refuses_bad_inputs_and_cancelled_scan_returns_no_partial_bound() {
    let grid = Quadtree::uniform(1);
    let sdf = HalfPlane {
        normal: [1.0, 0.0],
        offset: 2.0,
    };
    let material = IsotropicElastic::new(10.0, 0.3, 0.1).unwrap();
    let op = problem(&grid, &sdf, &material)
        .assemble(&|_, _| [0.0; 2], &|_, _| [0.0; 2])
        .unwrap();
    let mut u = vec![0.0; op.dof_count()];
    assert!(
        op.strain_stress(&grid, &u[..u.len() - 1], (1, 0, 0), [0.2; 2])
            .is_err()
    );
    assert!(op.strain_stress(&grid, &u, (2, 0, 0), [0.2; 2]).is_err());
    assert!(
        op.strain_stress(&grid, &u, (1, 0, 0), [f64::NAN, 0.2])
            .is_err()
    );
    assert!(op.strain_stress(&grid, &u, (1, 0, 0), [0.6, 0.2]).is_err());
    let mut calls = 0;
    let stopped = op
        .admit_strain_domain(&grid, &u, || {
            calls += 1;
            if calls == 2 {
                ControlFlow::Break("stop")
            } else {
                ControlFlow::Continue(())
            }
        })
        .unwrap();
    assert_eq!(stopped, ControlFlow::Break("stop"));
    u[0] = f64::NAN;
    assert!(op.strain_stress(&grid, &u, (1, 0, 0), [0.2; 2]).is_err());
}
