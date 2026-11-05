//! G1/G3: directional sheet transport against analytic flux and decay modes.
use super::*;

const AXES: [[f64; 3]; 2] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];

fn rotate(v: [f64; 3]) -> [f64; 3] {
    [
        (v[0] + 2.0 * v[1] - 2.0 * v[2]) / 3.0,
        (2.0 * v[0] + v[1] + 2.0 * v[2]) / 3.0,
        (2.0 * v[0] - 2.0 * v[1] - v[2]) / 3.0,
    ]
}

fn tetrahedron(rotated: bool) -> (ConductionMesh, usize) {
    let mut points = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    if rotated {
        points.iter_mut().for_each(|p| *p = rotate(*p));
    }
    let mesh = ConductionMesh::new(TetComplex::from_tets(4, vec![[0, 1, 2, 3]]), points).unwrap();
    let face = mesh
        .boundary()
        .iter()
        .position(|f| !f.vertices.contains(&3))
        .unwrap();
    (mesh, face)
}

#[test]
fn rotated_p1_tensors_retain_exact_flux_energy_and_constant_nullspace() {
    with_cx(|cx| {
        for rotated in [false, true] {
            let (mesh, face) = tetrahedron(rotated);
            let axes = if rotated { AXES.map(rotate) } else { AXES };
            // Mean tensor [[8,2],[2,4]]. T = 300 + 3x - 2y + 5z.
            // The z gradient does not enter the intrinsic sheet law.
            let nodal = [1.0, 2.0, 3.0].map(|s| [[4.0 * s, s], [s, 2.0 * s]]);
            let film = SurfaceHeatConduction::nodal_tensor(face, axes, nodal).unwrap();
            assert_eq!(film.nodal_conductance_w_per_k(), None);
            assert_eq!(film.tangent_axes(), Some(axes));
            assert_eq!(film.nodal_tensor_conductance_w_per_k(), Some(nodal));
            let temperature = [300.0, 303.0, 298.0, 305.0];
            let flux = film.mean_heat_flux_w_per_m(&mesh, &temperature).unwrap();
            let expected = if rotated {
                rotate([-20.0, 2.0, 0.0])
            } else {
                [-20.0, 2.0, 0.0]
            };
            for i in 0..3 {
                close(flux[i], expected[i], 2e-13);
            }
            assert_eq!(
                film.mean_heat_flux_w_per_m(&mesh, &[1e12; 4]).unwrap(),
                [0.0; 3]
            );
            let base = ThermalBoundaryBuilder::new(&mesh)
                .adiabatic_remainder()
                .finish()
                .unwrap();
            let coated = base
                .clone()
                .with_surface_heat_conduction(cx, &mesh, &[film])
                .unwrap();
            let material = ConductivityModel::isotropic_declared(2.0).unwrap();
            let source = ScalarField::Uniform(0.0);
            let residual = |boundary, values: &[f64]| {
                let system =
                    assemble_operator(cx, &mesh, boundary, &material, &source, values).unwrap();
                full_residual(&system, values)
            };
            let bare = residual(&base, &temperature);
            let total = residual(&coated, &temperature);
            let transfer: Vec<_> = total.iter().zip(&bare).map(|(a, b)| a - b).collect();
            close(transfer.iter().sum(), 0.0, 2e-12);
            close(transfer[3], 0.0, 1e-13);
            // A * grad(T)^T K grad(T) = 1/2 * (3*20 + (-2)*(-2)).
            let dissipation: f64 = transfer
                .iter()
                .zip(temperature)
                .map(|(r, t)| r * (t - 300.0))
                .sum();
            close(dissipation, 32.0, 2e-11);
            for value in residual(&coated, &[300.0; 4]) {
                close(value, 0.0, 2e-12);
            }
        }
    });
}

#[test]
fn tensor_and_scalar_isotropic_limits_share_the_original_operator() {
    with_cx(|cx| {
        let (mesh, face) = tetrahedron(false);
        let values = [2.0, 3.0, 7.0];
        let scalar = SurfaceHeatConduction::nodal(face, values).unwrap();
        let tensor =
            SurfaceHeatConduction::nodal_tensor(face, AXES, values.map(|k| [[k, 0.0], [0.0, k]]))
                .unwrap();
        let base = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let a = base
            .clone()
            .with_surface_heat_conduction(cx, &mesh, &[scalar])
            .unwrap();
        let b = base
            .with_surface_heat_conduction(cx, &mesh, &[tensor])
            .unwrap();
        let material = ConductivityModel::isotropic_declared(2.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        for vertex in 0..4 {
            let mut temperature = vec![300.0; 4];
            temperature[vertex] += 1.0;
            let system_a =
                assemble_operator(cx, &mesh, &a, &material, &source, &temperature).unwrap();
            let system_b =
                assemble_operator(cx, &mesh, &b, &material, &source, &temperature).unwrap();
            for (ra, rb) in full_residual(&system_a, &temperature)
                .iter()
                .zip(full_residual(&system_b, &temperature))
            {
                close(*ra, rb, 2e-12);
            }
        }
    });
}

#[test]
fn anisotropic_sheet_and_original_storage_follow_the_analytic_decay_mode() {
    with_cx(|cx| {
        let (mesh, face) = tetrahedron(false);
        let film =
            SurfaceHeatConduction::uniform_tensor(face, AXES, [[5.0, 3.0], [3.0, 5.0]]).unwrap();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap()
            .with_surface_heat_conduction(cx, &mesh, &[film])
            .unwrap();
        // V=1/6, lumped capacity C_i=(24*V)/4=1 J/K. This mode has
        // bulk eigenvalue k*V=1 W/K and sheet eigenvalue (5-3)/2=1 W/K.
        let material = ConductivityModel::isotropic_declared(6.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let engine =
            BackwardEuler::uniform(cx, &mesh, VolumetricHeatCapacity::declared(24.0).unwrap())
                .unwrap();
        let old = [300.0, 310.0, 290.0, 300.0];
        let step = engine
            .advance(cx, problem, None, &old, 0.5, config())
            .unwrap();
        for (actual, expected) in step.temperature.iter().zip([300.0, 305.0, 295.0, 300.0]) {
            close(*actual, expected, 2e-9);
        }
        close(step.stored_energy_change_j, 0.0, 1e-8);
        close(step.energy_residual_j, 0.0, 1e-8);
        assert_eq!(old, [300.0, 310.0, 290.0, 300.0]);
    });
}

#[test]
fn rank_one_sheets_admit_insulated_directions_and_refuse_nonphysical_laws() {
    with_cx(|cx| {
        let (mesh, face) = tetrahedron(false);
        let rank_one =
            SurfaceHeatConduction::uniform_tensor(face, AXES, [[4.0, 2.0], [2.0, 1.0]]).unwrap();
        assert_eq!(
            rank_one
                .mean_heat_flux_w_per_m(&mesh, &[300.0, 301.0, 298.0, 300.0])
                .unwrap(),
            [0.0; 3]
        );
        let flux = rank_one
            .mean_heat_flux_w_per_m(&mesh, &[300.0, 303.0, 298.0, 300.0])
            .unwrap();
        assert_eq!(flux, [-8.0, -4.0, 0.0]);
        for tensor in [
            [[1.0, 2.0], [2.0, 1.0]],
            [[-1.0, 0.0], [0.0, 2.0]],
            [[1.0, 0.5], [0.0, 2.0]],
            [[f64::NAN, 0.0], [0.0, 2.0]],
        ] {
            assert!(SurfaceHeatConduction::uniform_tensor(face, AXES, tensor).is_err());
        }
        assert!(
            SurfaceHeatConduction::uniform_tensor(face, [AXES[0]; 2], [[1.0, 0.0], [0.0, 1.0]])
                .is_err()
        );
        let base = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let wrong_plane = SurfaceHeatConduction::uniform_tensor(
            face,
            [AXES[0], [0.0, 0.0, 1.0]],
            [[1.0, 0.0], [0.0, 1.0]],
        )
        .unwrap();
        assert!(
            base.clone()
                .with_surface_heat_conduction(cx, &mesh, &[wrong_plane])
                .is_err()
        );
        assert!(
            base.clone()
                .with_surface_heat_conduction(cx, &mesh, &[rank_one, rank_one])
                .is_err()
        );
        let zero = SurfaceHeatConduction::uniform_tensor(face, AXES, [[0.0; 2]; 2]).unwrap();
        assert!(
            base.with_surface_heat_conduction(cx, &mesh, &[zero])
                .unwrap()
                .surface_heat_conduction()
                .is_empty()
        );
    });
}
