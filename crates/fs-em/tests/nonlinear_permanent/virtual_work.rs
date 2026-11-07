//! Permanent-magnet torque from original fields and independently rotated solves.

use super::{common, induction, recoil};
use fs_em::planar::{
    MagneticError, PlanarMagnetostatics,
    nonlinear::{IsotropicMagneticMaterial, PlanarNonlinearMagnetostatics},
    permanent::{PlanarMagneticMaterial, PlanarRecoilMaterial},
};
use fs_rep_mesh::{Metric2, TriComplex2};

const DEPTH: f64 = 0.2;
const REMANENCE: f64 = 0.6;
const BACKGROUND: f64 = 1.2;

fn mesh(n: usize) -> TriComplex2 {
    let vertices = (0..=n)
        .flat_map(|y| {
            (0..=n).map(move |x| {
                [
                    -1.0 + 2.0 * x as f64 / n as f64,
                    -1.0 + 2.0 * y as f64 / n as f64,
                ]
            })
        })
        .collect();
    let faces = (0..n)
        .flat_map(|y| {
            (0..n).flat_map(move |x| {
                let a = (y * (n + 1) + x) as u32;
                let d = a + n as u32 + 1;
                [[a, a + 1, d + 1], [a, d + 1, d]]
            })
        })
        .collect();
    TriComplex2::from_indexed_triangles(
        "rigid recoil rotor",
        vertices,
        faces,
        Metric2::planar(DEPTH).unwrap(),
    )
    .unwrap()
}

fn center(mesh: &TriComplex2, face: &[u32; 3]) -> [f64; 2] {
    [0, 1].map(|axis| {
        face.iter()
            .map(|&v| mesh.vertices()[v as usize][axis])
            .sum::<f64>()
            / 3.0
    })
}

fn materials(mesh: &TriComplex2, nonlinear: bool) -> Vec<IsotropicMagneticMaterial> {
    let (source, _) = recoil(1.0, REMANENCE, [-10.0, 10.0], 10.0, 18.0);
    let magnet = PlanarRecoilMaterial::new(source, [1.0, 0.0]).unwrap();
    let exterior = if nonlinear {
        IsotropicMagneticMaterial::Reversible(induction::law(&induction::source(
            &[(0.0, 0.0), (1.0, 1.0), (4.0, 13.0)],
            18.0,
        )))
    } else {
        IsotropicMagneticMaterial::Linear(common::material(1.0))
    };
    mesh.faces()
        .iter()
        .map(|face| {
            let [x, y] = center(mesh, face);
            if x.abs().max(y.abs()) < 0.25 {
                IsotropicMagneticMaterial::Recoil(magnet.clone())
            } else {
                exterior.clone()
            }
        })
        .collect()
}

fn linear(mesh: TriComplex2) -> PlanarMagnetostatics {
    let models = materials(&mesh, false)
        .into_iter()
        .map(|law| match law {
            IsotropicMagneticMaterial::Linear(law) => PlanarMagneticMaterial::Linear(law),
            IsotropicMagneticMaterial::Recoil(law) => PlanarMagneticMaterial::Recoil(law),
            IsotropicMagneticMaterial::Reversible(_) => unreachable!(),
        })
        .collect();
    PlanarMagnetostatics::new_with_material_models(mesh, models, common::limits(), &mut || false)
        .unwrap()
}

fn boundary(mesh: &TriComplex2) -> Vec<(usize, f64)> {
    mesh.boundary_trace()
        .unwrap()
        .vertices
        .into_iter()
        .map(|v| (v as usize, BACKGROUND * mesh.vertices()[v as usize][0]))
        .collect()
}

fn weight([x, y]: [f64; 2]) -> f64 {
    ((0.75 - x.abs().max(y.abs())) / 0.5).clamp(0.0, 1.0)
}

fn motion(mesh: &TriComplex2) -> Vec<[f64; 2]> {
    mesh.vertices()
        .iter()
        .map(|&x| {
            let w = weight(x);
            [-w * x[1], w * x[0]]
        })
        .collect()
}

fn rotate([x, y]: [f64; 2], angle: f64) -> [f64; 2] {
    let (s, c) = angle.sin_cos();
    [c * x - s * y, s * x + c * y]
}

// Exact finite rotation in every magnet; a smooth-in-angle P1 extension in
// the current-free collar. Outer stator cells and the boundary remain fixed.
fn rotated_mesh(mesh: &TriComplex2, angle: f64) -> TriComplex2 {
    let vertices = mesh
        .vertices()
        .iter()
        .map(|&x| {
            let rotated = rotate(x, angle);
            [0, 1].map(|axis| weight(x).mul_add(rotated[axis] - x[axis], x[axis]))
        })
        .collect();
    TriComplex2::from_indexed_triangles(
        "rotated original recoil mesh",
        vertices,
        mesh.faces().to_vec(),
        mesh.metric(),
    )
    .unwrap()
}

fn rotated_magnet(magnet: &PlanarRecoilMaterial, angle: f64) -> PlanarRecoilMaterial {
    PlanarRecoilMaterial::new(magnet.material().clone(), rotate(magnet.direction(), angle)).unwrap()
}

fn source_pairing(mesh: &TriComplex2, potential: &[f64], currents: &[f64]) -> f64 {
    mesh.faces()
        .iter()
        .enumerate()
        .map(|(cell, face)| {
            currents[cell]
                * mesh.face_measure(cell).unwrap()
                * face.iter().map(|&v| potential[v as usize]).sum::<f64>()
                / 3.0
        })
        .sum()
}

#[test]
fn rigid_permanent_magnet_torque_matches_finite_rotations_and_requires_carried_remanence() {
    let problem = linear(mesh(16));
    let currents = vec![0.0; problem.mesh().faces().len()];
    let bc = boundary(problem.mesh());
    let field = problem
        .solve(&currents, &bc, common::config(), &mut || false)
        .unwrap();
    let work = field
        .rigid_magnet_virtual_work(&motion(problem.mesh()), &mut || false)
        .unwrap();
    // Independent continuum dipole moment: m=volume*Br/mu along +x, B along
    // -y. Equal magnet and exterior permeability isolates alignment torque.
    let exact = -0.25 * DEPTH * REMANENCE * BACKGROUND;
    assert!(
        (work - exact).abs() < 0.005 * exact.abs(),
        "{work} != {exact}"
    );
    let epsilon = 1e-4;
    for carried in [true, false] {
        let values: Vec<_> = [-epsilon, epsilon]
            .into_iter()
            .map(|angle| {
                let models = problem
                    .materials()
                    .iter()
                    .zip(problem.recoil_materials())
                    .map(|(law, magnet)| match magnet {
                        Some(magnet) => PlanarMagneticMaterial::Recoil(rotated_magnet(
                            magnet,
                            if carried { angle } else { 0.0 },
                        )),
                        None => PlanarMagneticMaterial::Linear(law.clone()),
                    })
                    .collect();
                let moved = PlanarMagnetostatics::new_with_material_models(
                    rotated_mesh(problem.mesh(), angle),
                    models,
                    common::limits(),
                    &mut || false,
                )
                .unwrap();
                let solved = moved
                    .solve(&currents, &bc, common::config(), &mut || false)
                    .unwrap();
                solved.energy_j() - source_pairing(moved.mesh(), solved.potential_wb_m(), &currents)
            })
            .collect();
        let difference = -(values[1] - values[0]) / (2.0 * epsilon);
        if carried {
            assert!(
                (work - difference).abs() < 2e-6 * work.abs(),
                "{work} != {difference}"
            );
        } else {
            // Rotating geometry but freezing the magnet's direction is a
            // different physical control and misses the alignment torque.
            assert!((work - difference).abs() > 0.5 * work.abs());
        }
    }
}

#[test]
fn rigid_magnet_in_nonlinear_iron_matches_the_original_mixed_control_potential() {
    let mesh = mesh(8);
    let models = materials(&mesh, true);
    let problem =
        PlanarNonlinearMagnetostatics::new(mesh, models, common::limits(), &mut || false).unwrap();
    let mesh = problem.mesh();
    // Fixed stator current lives entirely outside the moving collar.
    let currents: Vec<_> = mesh
        .faces()
        .iter()
        .map(|face| {
            let [x, y] = center(mesh, face);
            if x.abs().max(y.abs()) > 0.75 {
                0.4 * y.signum()
            } else {
                0.0
            }
        })
        .collect();
    let initial: Vec<_> = mesh.vertices().iter().map(|x| BACKGROUND * x[0]).collect();
    let bc = boundary(mesh);
    let controls = induction::controls().nonlinear;
    let field = problem
        .solve(&currents, &bc, &initial, controls, &mut || false)
        .unwrap();
    assert!(
        problem
            .materials()
            .iter()
            .zip(field.magnetic_flux_density_t())
            .any(|(law, b)| {
                matches!(law, IsotropicMagneticMaterial::Reversible(_)) && b[0].hypot(b[1]) > 1.0
            })
    );
    let work = field
        .rigid_magnet_virtual_work(&motion(mesh), &mut || false)
        .unwrap();
    let epsilon = 1e-4;
    let values: Vec<_> = [-epsilon, epsilon]
        .into_iter()
        .map(|angle| {
            let models = problem
                .materials()
                .iter()
                .map(|law| match law {
                    IsotropicMagneticMaterial::Recoil(magnet) => {
                        IsotropicMagneticMaterial::Recoil(rotated_magnet(magnet, angle))
                    }
                    other => other.clone(),
                })
                .collect();
            let moved = PlanarNonlinearMagnetostatics::new(
                rotated_mesh(mesh, angle),
                models,
                common::limits(),
                &mut || false,
            )
            .unwrap();
            let moved_current: Vec<_> = currents
                .iter()
                .enumerate()
                .map(|(cell, j)| {
                    j * mesh.face_measure(cell).unwrap() / moved.mesh().face_measure(cell).unwrap()
                })
                .collect();
            let solved = moved
                .solve(
                    &moved_current,
                    &bc,
                    field.potential_wb_m(),
                    controls,
                    &mut || false,
                )
                .unwrap();
            solved.energy_j()
                - source_pairing(moved.mesh(), solved.potential_wb_m(), &moved_current)
        })
        .collect();
    let difference = -(values[1] - values[0]) / (2.0 * epsilon);
    assert!(
        (work - difference).abs() < 2e-6 * work.abs(),
        "{work} != {difference}"
    );
}

#[test]
fn rigid_recoil_motion_is_objective_and_straining_magnets_refuse_atomically() {
    let problem = linear(mesh(8));
    let field = problem
        .solve(
            &vec![0.0; problem.mesh().faces().len()],
            &boundary(problem.mesh()),
            common::config(),
            &mut || false,
        )
        .unwrap();
    let n = problem.mesh().vertices().len();
    let original = field.potential_wb_m().to_vec();
    let rotation: Vec<_> = problem
        .mesh()
        .vertices()
        .iter()
        .map(|x| [-x[1], x[0]])
        .collect();
    assert!(
        field
            .rigid_magnet_virtual_work(&rotation, &mut || false)
            .unwrap()
            .abs()
            < 1e-11
    );
    assert_eq!(
        field
            .rigid_magnet_virtual_work(&vec![[2_f64.powi(40), -2_f64.powi(40)]; n], &mut || false)
            .unwrap(),
        0.0
    );
    assert!(field.rigid_magnet_virtual_work(&[], &mut || false).is_err());
    for shear in [false, true] {
        let strained: Vec<_> = problem
            .mesh()
            .vertices()
            .iter()
            .map(|x| if shear { [x[1], 0.0] } else { *x })
            .collect();
        assert!(
            field
                .rigid_magnet_virtual_work(&strained, &mut || false)
                .is_err()
        );
    }
    let motion = motion(problem.mesh());
    let mut calls = 0;
    let reference = field
        .rigid_magnet_virtual_work(&motion, &mut || {
            calls += 1;
            false
        })
        .unwrap();
    let mut repeated = 0;
    assert!(matches!(
        field.rigid_magnet_virtual_work(&motion, &mut || {
            repeated += 1;
            repeated == calls
        }),
        Err(MagneticError::Cancelled)
    ));
    assert_eq!(field.potential_wb_m(), original);
    assert_eq!(
        field
            .rigid_magnet_virtual_work(&motion, &mut || false)
            .unwrap(),
        reference
    );
}
