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

#[test]
fn bulk_zero_normal_flux_reduction_matches_relaxed_fourier_flux_and_original_assembly() {
    with_cx(|cx| {
        let bulk = [[6.0, 2.0, 2.0], [2.0, 5.0, 1.0], [2.0, 1.0, 2.0]];
        // For grad_t=[3,-2], q_n=0 requires grad_n=-2. Independent 3-D
        // Fourier multiplication gives q_bulk=[-10,6,0], hence t*q_t=[-2.5,1.5].
        let relaxed_gradient = [3.0, -2.0, -2.0];
        let bulk_flux = bulk.map(|row| {
            -row.iter()
                .zip(relaxed_gradient)
                .map(|(k, g)| k * g)
                .sum::<f64>()
        });
        assert_eq!(bulk_flux, [-10.0, 6.0, 0.0]);
        for rotated in [false, true] {
            let (mesh, face) = tetrahedron(rotated);
            let axes = if rotated { AXES.map(rotate) } else { AXES };
            let world_bulk = if rotated {
                let columns = [AXES[0], AXES[1], [0.0, 0.0, 1.0]].map(rotate);
                let mut result = [[0.0; 3]; 3];
                for i in 0..3 {
                    for j in i..3 {
                        let mut entry = 0.0;
                        for a in 0..3 {
                            for b in 0..3 {
                                entry += columns[a][i] * bulk[a][b] * columns[b][j];
                            }
                        }
                        result[i][j] = entry;
                        result[j][i] = entry;
                    }
                }
                result
            } else {
                bulk
            };
            let reduced = SurfaceHeatConduction::uniform_from_bulk_zero_normal_flux(
                face, axes, world_bulk, 0.25,
            )
            .unwrap();
            let expected_tensor = [[1.0, 0.25], [0.25, 1.125]];
            for tensor in reduced.nodal_tensor_conductance_w_per_k().unwrap() {
                for i in 0..2 {
                    for j in 0..2 {
                        close(tensor[i][j], expected_tensor[i][j], 3e-14);
                    }
                }
            }
            let temperature = [300.0, 303.0, 298.0, 305.0];
            let expected_flux = bulk_flux.map(|q| 0.25 * q);
            let expected_flux = if rotated {
                rotate(expected_flux)
            } else {
                expected_flux
            };
            let flux = reduced.mean_heat_flux_w_per_m(&mesh, &temperature).unwrap();
            for i in 0..3 {
                close(flux[i], expected_flux[i], 2e-13);
            }
            let boundary = ThermalBoundaryBuilder::new(&mesh)
                .adiabatic_remainder()
                .finish()
                .unwrap();
            let reduced_boundary = boundary
                .clone()
                .with_surface_heat_conduction(cx, &mesh, &[reduced])
                .unwrap();
            let analytic =
                SurfaceHeatConduction::uniform_tensor(face, axes, expected_tensor).unwrap();
            let analytic_boundary = boundary
                .with_surface_heat_conduction(cx, &mesh, &[analytic])
                .unwrap();
            let material = ConductivityModel::isotropic_declared(2.0).unwrap();
            let source = ScalarField::Uniform(0.0);
            let residual = |boundary| {
                let system =
                    assemble_operator(cx, &mesh, boundary, &material, &source, &temperature)
                        .unwrap();
                full_residual(&system, &temperature)
            };
            for (actual, expected) in residual(&reduced_boundary)
                .iter()
                .zip(residual(&analytic_boundary))
            {
                close(*actual, expected, 2e-12);
            }
        }
    });
}

#[test]
fn bulk_sheet_reduction_recovers_isotropy_and_normal_reversal() {
    let (mesh, face) = tetrahedron(false);
    let bulk = [[6.0, 2.0, 2.0], [2.0, 5.0, 1.0], [2.0, 1.0, 2.0]];
    let temperature = [300.0, 303.0, 298.0, 305.0];
    for axes in [AXES, [AXES[0], [0.0, -1.0, 0.0]]] {
        let film =
            SurfaceHeatConduction::uniform_from_bulk_zero_normal_flux(face, axes, bulk, 0.25)
                .unwrap();
        for (actual, expected) in film
            .mean_heat_flux_w_per_m(&mesh, &temperature)
            .unwrap()
            .into_iter()
            .zip([-2.5, 1.5, 0.0])
        {
            close(actual, expected, 2e-13);
        }
        let isotropic = SurfaceHeatConduction::uniform_from_bulk_zero_normal_flux(
            face,
            axes,
            [[6.0, 0.0, 0.0], [0.0, 6.0, 0.0], [0.0, 0.0, 6.0]],
            0.25,
        )
        .unwrap();
        assert_eq!(
            isotropic.nodal_tensor_conductance_w_per_k(),
            Some([[[1.5, 0.0], [0.0, 1.5]]; 3])
        );
    }
}

#[test]
fn bulk_sheet_reduction_retains_rank_one_and_refuses_invalid_constitutive_constraints() {
    let (mesh, face) = tetrahedron(false);
    let bulk = [[1.0, 0.0, 2.0], [0.0, 2.0, 0.0], [2.0, 0.0, 4.0]];
    let film =
        SurfaceHeatConduction::uniform_from_bulk_zero_normal_flux(face, AXES, bulk, 0.25).unwrap();
    assert_eq!(
        film.nodal_tensor_conductance_w_per_k(),
        Some([[[0.0, 0.0], [0.0, 0.5]]; 3])
    );
    assert_eq!(
        film.mean_heat_flux_w_per_m(&mesh, &[300.0, 301.0, 300.0, 300.0])
            .unwrap(),
        [0.0; 3]
    );
    for invalid_bulk in [
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]],
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]],
        [[1.0, 0.0, 2.0], [0.0, 1.0, 0.0], [2.0, 0.0, 1.0]],
        [[1.0, 0.5, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        [[f64::NAN, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        [
            [f64::from_bits(1), 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1e308],
        ],
    ] {
        assert!(
            SurfaceHeatConduction::uniform_from_bulk_zero_normal_flux(
                face,
                AXES,
                invalid_bulk,
                0.25
            )
            .is_err()
        );
    }
    for thickness in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(
            SurfaceHeatConduction::uniform_from_bulk_zero_normal_flux(face, AXES, bulk, thickness)
                .is_err()
        );
    }
    assert!(
        SurfaceHeatConduction::uniform_from_bulk_zero_normal_flux(face, [AXES[0]; 2], bulk, 0.25)
            .is_err()
    );
}

#[test]
fn bulk_sheet_reduction_avoids_an_unrepresentable_common_conductance_scale() {
    let bulk = [[1e300, 0.0, 0.0], [0.0, 2e300, 0.0], [0.0, 0.0, 1e308]];
    // max(|K|)*t overflows, although both physical sheet entries are finite.
    let film =
        SurfaceHeatConduction::uniform_from_bulk_zero_normal_flux(0, AXES, bulk, 10.0).unwrap();
    let sheet = film.nodal_tensor_conductance_w_per_k().unwrap()[0];
    close(sheet[0][0] / 1e301, 1.0, 3e-15);
    close(sheet[1][1] / 2e301, 1.0, 3e-15);
    assert_eq!(sheet[0][1], 0.0);
    assert!(
        SurfaceHeatConduction::uniform_from_bulk_zero_normal_flux(
            0,
            AXES,
            [[1e308, 0.0, 0.0], [0.0, 1e308, 0.0], [0.0, 0.0, 1e308]],
            10.0,
        )
        .is_err()
    );
}
