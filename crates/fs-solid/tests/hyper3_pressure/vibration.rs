//! G1/G3: native modal pressure stiffness against scalar volume integrals.

use super::*;
use fs_solid::hyper3::vibration::{VibrationError, VibrationSettings};

fn settings() -> VibrationSettings {
    VibrationSettings {
        maximum_free_dofs: 128,
        relative_residual_tolerance: 1e-9,
    }
}

fn problem<'a>(
    material: &'a Hyperelastic,
    fixed: &'a [(usize, f64)],
    dead: &'a [(usize, f64)],
) -> HyperTetProblem<'a> {
    HyperTetProblem {
        nodes_m: &NODES,
        tetrahedra: &TETS,
        material: HyperelasticLawField::Uniform(material),
        prescribed_m: fixed,
        nodal_forces_n: dead,
        body_force_n_m3: [0.0; 3],
        budget: TetAssemblyBudget::standard(),
        settings: HyperTetSettings::default(),
    }
}

fn material() -> Hyperelastic {
    Hyperelastic::new(
        HyperelasticModel::NeoHookean {
            mu: 20.0,
            lambda: 30.0,
        },
        4.0,
    )
    .unwrap()
}

/// The sole free coordinate is the unit tet's height s. Differentiate its
/// original scalar neo-Hookean energy V0 W(diag(1,1,s)) independently.
fn elastic_force_and_stiffness(s: f64) -> (f64, f64) {
    (
        (20.0 * s + (30.0 * s.ln() - 20.0) / s) / 6.0,
        (20.0 + (20.0 + 30.0 * (1.0 - s.ln())) / (s * s)) / 6.0,
    )
}

fn vertical_supports() -> Vec<(usize, f64)> {
    (0..12)
        .filter(|&dof| dof != 11)
        .map(|dof| (dof, 0.0))
        .collect()
}

fn window_plan(lambda: f64) -> fs_modal::WindowPlan {
    fs_modal::WindowPlan {
        windows: vec![(0.5 * lambda, 1.5 * lambda)],
        threads: 1,
        maximum_windows: 1,
        maximum_bytes: 32 * 1024 * 1024,
        maximum_relative_halfwidth: 1e-8,
        coordinates: None,
    }
}

fn window_options() -> fs_modal::SliceOptions {
    fs_modal::SliceOptions {
        max_lanczos: 64,
        ritz_tol: 1e-12,
        ldlt: fs_sparse::LdltOptions {
            max_factor_bytes: Some(16 * 1024 * 1024),
            ..Default::default()
        },
        ..Default::default()
    }
}

#[test]
fn affine_field_modes_keep_world_gradient_and_original_graded_inertia() {
    with_cx(|_, cx| {
        let material = material();
        let fixed = vertical_supports();
        let surface =
            ClosedTetPressureSurface::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx)
                .unwrap();
        let (p0, gradient) = (9.0, -3.5);
        let field = TetPressure::Affine {
            origin_m: [0.0; 3],
            pressure_pa: p0,
            gradient_pa_m: [0.0, 0.0, gradient],
        };
        let density = [3.0, 4.0, 5.0, 9.0];
        for s in [1.0, 1.25] {
            // Integral p dV = p0*s/6 + gradient*s^2/24. Its Hessian is
            // gradient/12, even though a uniform follower pressure has zero
            // stiffness in this particular height coordinate.
            let (internal, material_k) = elastic_force_and_stiffness(s);
            let dead = [(11, internal + p0 / 6.0 + gradient * s / 12.0)];
            let p = problem(&material, &fixed, &dead);
            let expected_k = material_k + gradient / 12.0;
            for graded in [false, true] {
                let dense = if graded {
                    p.solve_vibration_with_nodal_density_and_pressure(
                        &surface,
                        field,
                        &density,
                        settings(),
                        cx,
                    )
                } else {
                    p.solve_vibration_with_pressure(&surface, field, &[6.0], settings(), cx)
                }
                .unwrap();
                let expected_mass = if graded {
                    (density.iter().sum::<f64>() + 2.0 * density[3]) / 360.0
                } else {
                    6.0 / 60.0
                };
                let expected_lambda = expected_k / expected_mass;
                near(dense.equilibrium().displacement_m[3][2], s - 1.0, 2e-9);
                assert_eq!(dense.free_dofs(), &[11]);
                near(dense.stiffness_n_m()[0], expected_k, 2e-8);
                near(dense.mass_kg()[0], expected_mass, 2e-14);
                near(dense.modes()[0].lambda, expected_lambda, 4e-7);
                assert!(dense.held_pressure().is_none());
                let (bound_surface, bound_field) = dense.held_pressure_field().unwrap();
                assert!(core::ptr::eq(bound_surface, &surface));
                assert!(matches!(
                    bound_field,
                    TetPressure::Affine {
                        origin_m: [0.0, 0.0, 0.0],
                        pressure_pa: 9.0,
                        gradient_pa_m: [0.0, 0.0, -3.5],
                    }
                ));
                let plan = window_plan(expected_lambda);
                let options = window_options();
                let sparse = if graded {
                    p.solve_vibration_windows_with_nodal_density_and_pressure(
                        &surface, field, &density, &plan, &options, cx,
                    )
                } else {
                    p.solve_vibration_windows_with_pressure(
                        &surface,
                        field,
                        &[6.0],
                        &plan,
                        &options,
                        cx,
                    )
                }
                .unwrap();
                assert_eq!(sparse.equilibrium(), dense.equilibrium());
                assert_eq!(sparse.stiffness_n_m().to_dense(), dense.stiffness_n_m());
                assert_eq!(sparse.mass_kg().to_dense(), dense.mass_kg());
                assert_eq!(sparse.spectrum().reports[0].expected, 1);
                near(
                    sparse.spectrum().reports[0].modes[0].lambda,
                    expected_lambda,
                    4e-7,
                );
                assert!(sparse.held_pressure().is_none());
                assert!(matches!(
                    sparse.held_pressure_field().unwrap().1,
                    TetPressure::Affine {
                        gradient_pa_m: [0.0, 0.0, -3.5],
                        ..
                    }
                ));
            }
        }
    });
}

#[test]
fn clipped_pressure_modes_use_the_actual_wet_volume_potential() {
    with_cx(|_, cx| {
        let material = material();
        let fixed = vertical_supports();
        let surface =
            ClosedTetPressureSurface::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx)
                .unwrap();
        let s = 1.2;
        let rate = 7.3;
        let (internal, material_k) = elastic_force_and_stiffness(s);
        for h in [-0.2, 0.5, 1.8] {
            // Integrate (h-z)_+ * (1-z/s)^2/2 dz directly over the current
            // tet, then differentiate in s. This never samples face tractions.
            let (load_gradient, load_hessian) = if h <= 0.0 {
                (0.0, 0.0)
            } else if h >= s {
                (rate * (h / 6.0 - s / 12.0), -rate / 12.0)
            } else {
                (
                    rate * (h.powi(3) / (6.0 * s * s) - h.powi(4) / (12.0 * s.powi(3))),
                    rate * (-h.powi(3) / (3.0 * s.powi(3)) + h.powi(4) / (4.0 * s.powi(4))),
                )
            };
            let field = TetPressure::FreeSurface {
                origin_m: [0.0; 3],
                pressure_pa: rate * h,
                gradient_pa_m: [0.0, 0.0, -rate],
            };
            let dead = [(11, internal + load_gradient)];
            let p = problem(&material, &fixed, &dead);
            let dense = p
                .solve_vibration_with_pressure(&surface, field, &[6.0], settings(), cx)
                .unwrap();
            near(dense.equilibrium().displacement_m[3][2], s - 1.0, 2e-9);
            let expected_k = material_k + load_hessian;
            near(dense.stiffness_n_m()[0], expected_k, 2e-8);
            near(dense.modes()[0].lambda, 10.0 * expected_k, 4e-7);
            let sparse = p
                .solve_vibration_windows_with_pressure(
                    &surface,
                    field,
                    &[6.0],
                    &window_plan(10.0 * expected_k),
                    &window_options(),
                    cx,
                )
                .unwrap();
            assert_eq!(sparse.stiffness_n_m().to_dense(), dense.stiffness_n_m());
            assert!(matches!(
                sparse.held_pressure_field().unwrap().1,
                TetPressure::FreeSurface { .. }
            ));
        }
    });
}

#[test]
fn immersed_body_translation_has_archimedes_waterplane_restoring_stiffness() {
    with_cx(|_, cx| {
        let material = material();
        let surface =
            ClosedTetPressureSurface::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx)
                .unwrap();
        let rate = 7.3;
        for h in [0.35_f64, 1.4] {
            let field = TetPressure::FreeSurface {
                origin_m: [0.0; 3],
                pressure_pa: rate * h,
                gradient_pa_m: [0.0, 0.0, -rate],
            };
            // Balance the existing load owner to isolate its equilibrium
            // tangent at the original shape. The restoring oracle below is
            // the independently known waterplane area, not that owner's Hessian.
            let dead: Vec<_> = surface
                .force_n(&[0.0; 12], field, cx)
                .unwrap()
                .into_iter()
                .enumerate()
                .map(|(dof, value)| (dof, -value))
                .collect();
            let p = problem(&material, &[], &dead);
            let dense = p
                .solve_vibration_with_pressure(&surface, field, &[6.0], settings(), cx)
                .unwrap();
            assert_eq!(dense.equilibrium().displacement_m, vec![[0.0; 3]; 4]);
            assert_eq!(dense.modes().len(), 12);
            let mut rigid_vertical_k = 0.0;
            for a in 0..4 {
                for b in 0..4 {
                    rigid_vertical_k += dense.stiffness_n_m()[12 * (3 * a + 2) + 3 * b + 2];
                }
            }
            let area = if h < 1.0 {
                (1.0 - h).powi(2) / 2.0
            } else {
                0.0
            };
            near(rigid_vertical_k, rate * area, 2e-12);
        }
    });
}

#[test]
fn conservative_pressure_modes_keep_uniform_limits_and_reject_face_patches() {
    with_cx(|gate, cx| {
        let material = material();
        let fixed = vertical_supports();
        let surface =
            ClosedTetPressureSurface::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx)
                .unwrap();
        for q in [0.0, 1.0] {
            let dead = [(11, q / 6.0)];
            let p = problem(&material, &fixed, &dead);
            let uniform = p
                .solve_vibration_with_pressure(&surface, q, &[6.0], settings(), cx)
                .unwrap();
            assert_eq!(uniform.held_pressure().unwrap().1, q);
            for field in [
                TetPressure::Affine {
                    origin_m: [2.0, -3.0, 7.0],
                    pressure_pa: q,
                    gradient_pa_m: [0.0; 3],
                },
                TetPressure::FreeSurface {
                    origin_m: [2.0, -3.0, 7.0],
                    pressure_pa: q,
                    gradient_pa_m: [0.0; 3],
                },
            ] {
                let spatial = p
                    .solve_vibration_with_pressure(&surface, field, &[6.0], settings(), cx)
                    .unwrap();
                assert_eq!(spatial.equilibrium(), uniform.equilibrium());
                assert_eq!(spatial.stiffness_n_m(), uniform.stiffness_n_m());
                assert_eq!(spatial.mass_kg(), uniform.mass_kg());
                assert_eq!(
                    spatial.modes()[0].lambda.to_bits(),
                    uniform.modes()[0].lambda.to_bits()
                );
            }
        }
        let p = problem(&material, &fixed, &[]);
        let patches = [1.0; 4];
        for field in [
            TetPressure::PerFace(&patches),
            TetPressure::Affine {
                origin_m: [f64::NAN, 0.0, 0.0],
                pressure_pa: 1.0,
                gradient_pa_m: [0.0; 3],
            },
            TetPressure::FreeSurface {
                origin_m: [0.0; 3],
                pressure_pa: 1.0,
                gradient_pa_m: [0.0, 0.0, f64::INFINITY],
            },
        ] {
            assert!(matches!(
                p.solve_vibration_with_pressure(&surface, field, &[6.0], settings(), cx),
                Err(VibrationError::Solid(HyperTetError::InvalidInput { .. }))
            ));
            assert!(matches!(
                p.solve_vibration_windows_with_pressure(
                    &surface,
                    field,
                    &[6.0],
                    &window_plan(1.0),
                    &window_options(),
                    cx
                ),
                Err(VibrationError::Solid(HyperTetError::InvalidInput { .. }))
            ));
        }
        gate.request();
        let field = TetPressure::Affine {
            origin_m: [0.0; 3],
            pressure_pa: 1.0,
            gradient_pa_m: [0.0, 0.0, -1.0],
        };
        assert!(matches!(
            p.solve_vibration_with_pressure(&surface, field, &[6.0], settings(), cx),
            Err(VibrationError::Solid(HyperTetError::Cancelled))
        ));
        assert!(matches!(
            p.solve_vibration_windows_with_pressure(
                &surface,
                field,
                &[6.0],
                &window_plan(1.0),
                &window_options(),
                cx
            ),
            Err(VibrationError::Solid(HyperTetError::Cancelled))
        ));
    });
}
