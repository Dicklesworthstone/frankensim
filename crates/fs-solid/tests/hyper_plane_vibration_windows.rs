//! G1/G3: loaded sparse sheet spectra on the original physical pencil.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::hyper::reference::{ReferenceDeformation, ReferencedHyperelastic};
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_modal::{ModalError, SliceOptions, WindowPlan};
use fs_solid::hyper2d::HyperelasticMaterials;
use fs_solid::hyper2d::vibration::{
    PlaneReferenceDensity, PlaneVibrationError, PlaneVibrationFormulation, PlaneVibrationSettings,
    VibrationSettings,
};
use fs_solid::{HyperProblem, Mesh2, NewtonSettings, Patch, SolidError};

const MASS_BYTES: usize = 4 * 1024 * 1024;

fn with_cx(run: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 73,
                    kernel_id: 61,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            &gate,
        )
    });
}
fn problem<'a>(mesh: &'a Mesh2, material: &'a dyn HyperelasticMaterials) -> HyperProblem<'a> {
    HyperProblem {
        mesh,
        material,
        dirichlet: vec![],
        traction: vec![],
        settings: NewtonSettings {
            tol: 2e-9,
            max_iters: 40,
            load_steps: 1,
            max_backtracks: 30,
        },
    }
}
fn dense_settings() -> PlaneVibrationSettings {
    PlaneVibrationSettings {
        spectrum: VibrationSettings {
            maximum_free_dofs: 256,
            relative_residual_tolerance: 1e-8,
        },
        maximum_mass_assembly_bytes: MASS_BYTES,
    }
}
fn plan(windows: Vec<(f64, f64)>) -> WindowPlan {
    WindowPlan {
        windows,
        threads: 2,
        maximum_windows: 4,
        maximum_bytes: 32 * 1024 * 1024,
        maximum_relative_halfwidth: 1e-8,
        coordinates: None,
    }
}
fn options() -> SliceOptions {
    SliceOptions {
        max_lanczos: 256,
        ritz_tol: 1e-11,
        ldlt: fs_sparse::LdltOptions {
            max_factor_bytes: Some(2 * 1024 * 1024),
            ..Default::default()
        },
        ..Default::default()
    }
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 3e-8 * b.abs().max(1.0), "{a:e} != {b:e}");
}

#[test]
fn loaded_strip_exceeds_dense_cap_and_matches_independent_discrete_modes() {
    with_cx(|cx, _| {
        let nx = 67;
        let (width, height, mu, rho, h, stretch) = (1.0, 0.7, 2.0, 3.0, 0.02, 1.3_f64);
        let mesh = Mesh2::quads(width, height, nx, 2);
        let law =
            Hyperelastic::new(HyperelasticModel::NeoHookean { mu, lambda: 4000.0 }, 4.0).unwrap();
        let boundary = |x: f64, y: f64| [(stretch - 1.0) * x, (stretch - 1.0) * y];
        let mut p = problem(&mesh, &law);
        p.dirichlet = [Patch::Left, Patch::Right, Patch::Bottom, Patch::Top]
            .into_iter()
            .map(|patch| (patch, &boundary as _))
            .collect();
        let initial: Vec<_> = mesh.nodes.iter().map(|&[x, y]| boundary(x, y)).collect();
        let density = vec![rho; mesh.elems.len()];
        let form =
            PlaneVibrationFormulation::IncompressiblePlaneStress(PlaneStressSettings::default());
        assert!(matches!(
            p.solve_vibration_from(
                &initial,
                form,
                h,
                PlaneReferenceDensity::PerElement(&density),
                dense_settings(),
                cx,
            ),
            Err(PlaneVibrationError::Solid(SolidError::InvalidInput { .. }))
        ));

        // One free interior row: opposite-cell mixed derivatives cancel.
        // Sx/Mx has exact sine eigenvalues 6(1-cos(theta))/(hx²(2+cos(theta))).
        // Sy/My=3/hy². Full transverse relaxation gives A=mu+3mu/s^6.
        let (hx, hy) = (width / nx as f64, height / 2.0);
        let a = mu + 3.0 * mu / stretch.powi(6);
        let mut expected = Vec::new();
        for j in 1..nx {
            let cosine = (core::f64::consts::PI * j as f64 / nx as f64).cos();
            let x = 6.0 * (1.0 - cosine) / (hx * hx * (2.0 + cosine));
            expected.extend([
                (a * x + mu * 3.0 / (hy * hy)) / rho,
                (mu * x + a * 3.0 / (hy * hy)) / rho,
            ]);
        }
        expected.sort_by(f64::total_cmp);
        let split = 0.5 * (expected[1] + expected[2]);
        let windows = plan(vec![
            (0.5 * expected[0], split),
            (split, 0.5 * (expected[3] + expected[4])),
        ]);
        let sparse = p
            .solve_vibration_windows_from(
                &initial,
                form,
                h,
                PlaneReferenceDensity::PerElement(&density),
                MASS_BYTES,
                &windows,
                &options(),
                cx,
            )
            .unwrap();
        assert_eq!(sparse.free_dofs().len(), 132);
        assert_eq!(sparse.stiffness_n_m().nrows(), 132);
        assert!(sparse.stiffness_n_m().nnz() < 10 * 132);
        near(sparse.reference_mass_kg(), rho * h * width * height);
        for (u, expected) in sparse.displacement_m().iter().zip(&initial) {
            near(u[0], expected[0]);
            near(u[1], expected[1]);
        }
        let report = sparse.spectrum();
        assert_eq!(report.receipt.pencil_symbolic_analyses, 1);
        assert_eq!(report.receipt.endpoint_factorizations, 3);
        for (index, window) in report.reports.iter().enumerate() {
            assert_eq!(
                (window.below_low, window.below_high, window.expected),
                (2 * index, 2 * index + 2, 2)
            );
            for (mode, &value) in window.modes.iter().zip(&expected[2 * index..2 * index + 2]) {
                near(mode.lambda, value);
                assert!(mode.residual <= windows.maximum_relative_halfwidth * mode.lambda.abs());
                let mut mphi = vec![0.0; 132];
                sparse.mass_kg().spmv(&mode.phi, &mut mphi);
                near(mode.phi.iter().zip(mphi).map(|(a, b)| a * b).sum(), 1.0);
            }
        }
    });
}

#[test]
fn sparse_and_dense_keep_identical_loaded_references_and_graded_nodal_inertia() {
    with_cx(|cx, _| {
        let law = Hyperelastic::new(
            HyperelasticModel::NeoHookean {
                mu: 2.0,
                lambda: 40.0,
            },
            4.0,
        )
        .unwrap();
        let reference = ReferencedHyperelastic::new(
            &law,
            ReferenceDeformation::try_new([1.02, 0.04, 0.0, 0.0, 0.97, 0.0, 0.0, 0.0, 1.08])
                .unwrap(),
        );
        let boundary = |x: f64, y: f64| [0.15 * x, 0.15 * y];
        for mesh in [
            Mesh2::triangles(1.0, 1.0, 2, 2),
            Mesh2::quads(1.0, 1.0, 2, 2),
        ] {
            let density: Vec<_> = mesh.nodes.iter().map(|&[x, y]| 2.0 + x + 2.0 * y).collect();
            let local = PlaneStressSettings::default();
            for form in [
                PlaneVibrationFormulation::PlaneStress(local),
                PlaneVibrationFormulation::IncompressiblePlaneStress(local),
                PlaneVibrationFormulation::MeanDilatation,
            ] {
                let mean = matches!(form, PlaneVibrationFormulation::MeanDilatation);
                if mean && mesh.elems[0].len() == 3 {
                    continue;
                }
                let material: &dyn HyperelasticMaterials = if mean { &law } else { &reference };
                let mut p = problem(&mesh, material);
                p.dirichlet = [Patch::Left, Patch::Right, Patch::Bottom, Patch::Top]
                    .into_iter()
                    .map(|patch| (patch, &boundary as _))
                    .collect();
                let input = PlaneReferenceDensity::Nodal(&density);
                let dense = p
                    .solve_vibration(form, 0.02, input, dense_settings(), cx)
                    .unwrap();
                let windows = plan(vec![(
                    0.5 * dense.modes()[0].lambda,
                    1.5 * dense.modes()[1].lambda,
                )]);
                let sparse = p
                    .solve_vibration_windows(
                        form,
                        0.02,
                        input,
                        MASS_BYTES,
                        &windows,
                        &options(),
                        cx,
                    )
                    .unwrap();
                assert!(core::ptr::eq(sparse.problem(), &p));
                assert_eq!(sparse.free_dofs(), dense.free_dofs());
                assert_eq!(sparse.displacement_m(), dense.displacement_m());
                assert_eq!(sparse.reactions_n(), dense.reactions_n());
                assert!(
                    matches!(sparse.reference_density(), PlaneReferenceDensity::Nodal(values) if values == density.as_slice())
                );
                for i in 0..2 {
                    for j in 0..2 {
                        assert_eq!(
                            sparse.stiffness_n_m().get(i, j),
                            dense.stiffness_n_m()[2 * i + j]
                        );
                        assert_eq!(sparse.mass_kg().get(i, j), dense.mass_kg()[2 * i + j]);
                    }
                }
                assert_eq!(sparse.spectrum().reports[0].expected, 2);
                for (a, b) in sparse.spectrum().reports[0].modes.iter().zip(dense.modes()) {
                    near(a.lambda, b.lambda);
                }
                near(sparse.reference_mass_kg(), 3.5 * 0.02);
                let loads = sparse
                    .reference_acceleration_load_n(&vec![[2.0, -1.0]; mesh.nodes.len()], cx)
                    .unwrap();
                near(
                    loads.iter().step_by(2).sum(),
                    2.0 * sparse.reference_mass_kg(),
                );
                near(
                    loads.iter().skip(1).step_by(2).sum(),
                    -sparse.reference_mass_kg(),
                );
            }
        }
    });
}

#[test]
fn unstable_windows_remain_signed_and_native_budgets_refuse_without_partial_modes() {
    with_cx(|cx, gate| {
        let mesh = Mesh2::quads(1.0, 1.0, 2, 2);
        let (mu, lambda, rho, s) = (2.0, 100.0, 3.0, 2.0_f64);
        let law = Hyperelastic::new(HyperelasticModel::NeoHookean { mu, lambda }, 4.0).unwrap();
        let boundary = |x: f64, y: f64| [(s - 1.0) * x, (s - 1.0) * y];
        let mut p = problem(&mesh, &law);
        p.dirichlet = [Patch::Left, Patch::Right, Patch::Bottom, Patch::Top]
            .into_iter()
            .map(|patch| (patch, &boundary as _))
            .collect();
        let initial: Vec<_> = mesh.nodes.iter().map(|&[x, y]| boundary(x, y)).collect();
        let density = vec![rho; mesh.elems.len()];
        let input = PlaneReferenceDensity::PerElement(&density);
        let form = PlaneVibrationFormulation::MeanDilatation;
        let expected = (mu * (24.0 + 12.0 / s.powi(2))
            + 9.0 * lambda * (1.0 - 2.0 * s.ln()) / s.powi(2))
            / rho;
        let windows = plan(vec![(1.2 * expected, 0.8 * expected), (1.0, 2.0)]);
        let native = options();
        let sparse = p
            .solve_vibration_windows_from(
                &initial, form, 0.02, input, MASS_BYTES, &windows, &native, cx,
            )
            .unwrap();
        assert_eq!(sparse.spectrum().reports[0].expected, 2);
        assert_eq!(sparse.spectrum().reports[1].expected, 0);
        for mode in &sparse.spectrum().reports[0].modes {
            near(mode.lambda, expected);
            assert!(mode.interval.1 < 0.0);
        }
        let mut short = native;
        short.max_lanczos = 1;
        assert!(matches!(
            p.solve_vibration_windows_from(
                &initial, form, 0.02, input, MASS_BYTES, &windows, &short, cx,
            ),
            Err(PlaneVibrationError::Modal(
                ModalError::WindowUnresolved { .. }
            ))
        ));
        assert!(matches!(
            p.solve_vibration_windows_from(&initial, form, 0.02, input, 0, &windows, &native, cx,),
            Err(PlaneVibrationError::Solid(SolidError::InvalidInput { .. }))
        ));
        let mut insufficient = windows.clone();
        insufficient.maximum_windows = 1;
        assert!(matches!(
            p.solve_vibration_windows_from(
                &initial,
                form,
                0.02,
                input,
                MASS_BYTES,
                &insufficient,
                &native,
                cx,
            ),
            Err(PlaneVibrationError::Modal(ModalError::WindowBudget { .. }))
        ));
        gate.request();
        assert!(matches!(
            p.solve_vibration_windows_from(
                &initial, form, 0.02, input, MASS_BYTES, &windows, &native, cx,
            ),
            Err(PlaneVibrationError::Solid(SolidError::Cancelled))
        ));
    });
}
