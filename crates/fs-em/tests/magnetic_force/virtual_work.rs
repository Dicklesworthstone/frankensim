//! Native shape derivatives, held electrical controls and a nonzero rotor torque.

use super::{AMBIENT_MU, DEPTH, boundary, common, config, limits, material, specimen};
use fs_em::planar::{MagneticError, PlanarField, PlanarMagnetostatics};
use fs_rep_mesh::TriComplex2;

fn displaced(
    original: &PlanarMagnetostatics,
    displacement: &[[f64; 2]],
    step: f64,
) -> PlanarMagnetostatics {
    let mesh = original.mesh();
    let positions = mesh
        .vertices()
        .iter()
        .zip(displacement)
        .map(|(x, u)| std::array::from_fn(|axis| step.mul_add(u[axis], x[axis])))
        .collect();
    PlanarMagnetostatics::new(
        TriComplex2::from_indexed_triangles(
            "displaced original magnetic problem",
            positions,
            mesh.faces().to_vec(),
            mesh.metric(),
        )
        .unwrap(),
        original.materials().to_vec(),
        limits(),
        &mut || false,
    )
    .unwrap()
}

fn potential(field: &PlanarField<'_>) -> f64 {
    let mesh = field.problem().mesh();
    let currents = field.impressed_current_density_a_m2().unwrap();
    let source_pairing: f64 = mesh
        .faces()
        .iter()
        .enumerate()
        .map(|(cell, face)| {
            currents[cell]
                * mesh.face_measure(cell).unwrap()
                * face
                    .iter()
                    .map(|&v| field.potential_wb_m()[v as usize])
                    .sum::<f64>()
                / 3.0
        })
        .sum();
    field.energy_j() - source_pairing
}

#[test]
fn fixed_boundary_flux_work_has_analytic_sign_scaling_and_rigid_motion_limits() {
    let mu = 2e-6;
    for (flux, depth) in [(1e-3, 0.2), (-1e-3, 0.2), (2e-3, 0.2), (1e-3, 0.6)] {
        let problem = PlanarMagnetostatics::new(
            common::layered(0.5, 1.0, depth),
            vec![material(mu); 4],
            limits(),
            &mut || false,
        )
        .unwrap();
        let field = problem
            .solve(&[0.0; 4], &common::ends(flux), config(), &mut || false)
            .unwrap();
        let energy = 0.5 * depth * flux * flux / mu;
        for (motion, expected) in [
            (
                problem
                    .mesh()
                    .vertices()
                    .iter()
                    .map(|x| [x[0], 0.0])
                    .collect::<Vec<_>>(),
                energy,
            ),
            (
                problem
                    .mesh()
                    .vertices()
                    .iter()
                    .map(|x| [0.0, x[1]])
                    .collect(),
                -energy,
            ),
            (
                problem
                    .mesh()
                    .vertices()
                    .iter()
                    .map(|x| [-x[1], x[0]])
                    .collect(),
                0.0,
            ),
            (problem.mesh().vertices().to_vec(), 0.0),
        ] {
            let work = field.virtual_work(&motion, &mut || false).unwrap();
            assert!((work - expected).abs() < 1e-10 * energy);
        }
        assert_eq!(
            field
                .virtual_work(&[[2_f64.powi(50), -2_f64.powi(50)]; 6], &mut || false)
                .unwrap(),
            0.0
        );
    }
}

#[test]
fn rotor_torque_matches_resolved_fixed_current_potential_derivative() {
    let (problem, _, body) = specimen(16, AMBIENT_MU, false, DEPTH);
    let mesh = problem.mesh();
    let mut currents = vec![0.0; mesh.faces().len()];
    for cell in body {
        let center_y = mesh.faces()[cell]
            .iter()
            .map(|&v| mesh.vertices()[v as usize][1])
            .sum::<f64>()
            / 3.0;
        currents[cell] = 1000.0 * center_y.signum();
    }
    let bc = boundary(&problem, |[x, _]| 1e-3 * x);
    let field = problem
        .solve(&currents, &bc, config(), &mut || false)
        .unwrap();
    // Equal outgoing/return currents in the upper/lower body halves give a
    // pure couple in B_y=-1e-3: zero net force, tau_z=-integral(|y| dV).
    // Unit rotation about the center. The body rotates rigidly;
    // only the current-free exterior deforms, and the outer boundary is fixed.
    // The independent continuum moment is -0.03125*depth N m.
    let motion: Vec<_> = mesh
        .vertices()
        .iter()
        .map(|&[x, y]| {
            let radius = x.abs().max(y.abs());
            let weight = ((0.75 - radius) / 0.5).clamp(0.0, 1.0);
            [-weight * y, weight * x]
        })
        .collect();
    let work = field.virtual_work(&motion, &mut || false).unwrap();
    let epsilon = 1e-4;
    let values: Vec<_> = [-epsilon, epsilon]
        .into_iter()
        .map(|step| {
            let moved = displaced(&problem, &motion, step);
            // Preserve integrated current, not current density, on EACH cell.
            let moved_current: Vec<_> = currents
                .iter()
                .enumerate()
                .map(|(cell, &j)| {
                    j * mesh.face_measure(cell).unwrap() / moved.mesh().face_measure(cell).unwrap()
                })
                .collect();
            let solution = moved
                .solve(&moved_current, &bc, config(), &mut || false)
                .unwrap();
            potential(&solution)
        })
        .collect();
    let difference = -(values[1] - values[0]) / (2.0 * epsilon);
    assert!(
        (work - difference).abs() < 2e-5 * work.abs(),
        "{work} != {difference}"
    );
    let exact = -0.03125 * DEPTH;
    assert!(
        (work - exact).abs() < 0.005 * exact.abs(),
        "{work} != {exact}"
    );
}

#[test]
fn malformed_displacements_and_late_cancellation_leave_the_field_unchanged() {
    let (problem, _, _) = specimen(8, AMBIENT_MU, false, DEPTH);
    let field = problem
        .solve(
            &vec![0.0; problem.mesh().faces().len()],
            &boundary(&problem, |[x, _]| 1e-3 * x),
            config(),
            &mut || false,
        )
        .unwrap();
    let motion = problem.mesh().vertices().to_vec();
    let original = field.potential_wb_m().to_vec();
    assert!(field.virtual_work(&[], &mut || false).is_err());
    let mut invalid = motion.clone();
    invalid[0][0] = f64::NAN;
    assert!(field.virtual_work(&invalid, &mut || false).is_err());
    let mut calls = 0;
    let reference = field
        .virtual_work(&motion, &mut || {
            calls += 1;
            false
        })
        .unwrap();
    let mut repeated = 0;
    assert!(matches!(
        field.virtual_work(&motion, &mut || {
            repeated += 1;
            repeated == calls
        }),
        Err(MagneticError::Cancelled)
    ));
    assert_eq!(field.potential_wb_m(), original);
    assert_eq!(
        field.virtual_work(&motion, &mut || false).unwrap(),
        reference
    );
}
