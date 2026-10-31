//! G1 open-system heat/reflection balances and a full nonlinear solid solve.
//! Prescribed view factors below are algebra fixtures, not geometric claims.
use super::*;
use fs_conduction::radiation::{OpenDiffuseEnclosure, solve_with_open_gray_diffuse_enclosure};

fn evidence() -> ViewFactorEvidence {
    ViewFactorEvidence::Analytic {
        geometry: "authored-open-exchange-algebra-fixture".into(),
    }
}

fn near(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-30),
        "{actual:e} versus {expected:e}"
    );
}

fn surface(mesh: &ConductionMesh, name: &str, x: f64, epsilon: f64) -> RadiationSurface {
    RadiationSurface::new(
        mesh,
        name,
        |f| on_box_face(f.centroid[0], x),
        emissivity(name, epsilon, name),
    )
    .unwrap()
}

fn unit_mesh() -> ConductionMesh {
    let (complex, positions) = box_grid([3, 2, 2], [1.0; 3]);
    ConductionMesh::new(complex, positions).unwrap()
}

#[test]
fn g1_partial_views_keep_multiple_reflection_and_independent_ambient_heat() {
    let mesh = two_disconnected_slabs();
    with_cx(|cx| {
        let model = OpenDiffuseEnclosure::new(
            cx,
            vec![
                surface(&mesh, "hot", 0.5, 0.4),
                surface(&mesh, "cold", 1.5, 0.7),
            ],
            vec![vec![0.0, 0.25], vec![0.25, 0.0]],
            evidence(),
            ViewFactorTolerance::default(),
        )
        .unwrap();
        let actual = model.solve(cx, &[400.0, 300.0], 280.0).unwrap();
        // Independent two-equation elimination in unshifted absolute radiosity.
        let e = [0.4, 0.7];
        let eb = STEFAN_BOLTZMANN_W_M2_K4 * 280.0_f64.powi(4);
        let emission = [400.0_f64, 300.0].map(|t| STEFAN_BOLTZMANN_W_M2_K4 * t.powi(4));
        let rhs = [
            e[0] * emission[0] + (1.0 - e[0]) * 0.75 * eb,
            e[1] * emission[1] + (1.0 - e[1]) * 0.75 * eb,
        ];
        let j0 = (rhs[0] + (1.0 - e[0]) * 0.25 * rhs[1])
            / (1.0 - (1.0 - e[0]) * (1.0 - e[1]) * 0.25 * 0.25);
        let j1 = rhs[1] + (1.0 - e[1]) * 0.25 * j0;
        let heat = [
            e[0] * (emission[0] - 0.25 * j1 - 0.75 * eb),
            e[1] * (emission[1] - 0.25 * j0 - 0.75 * eb),
        ];
        near(actual.surface_heat_w()[0], heat[0], 2e-14);
        near(actual.surface_heat_w()[1], heat[1], 2e-14);
        near(actual.ambient_heat_w(), 0.75 * (2.0 * eb - j0 - j1), 2e-14);
        assert!(actual.surface_heat_w().iter().sum::<f64>() > 0.0);
        assert!(actual.ambient_heat_w() < 0.0);
        assert!(actual.exchange.enclosure_energy_closure_w.abs() < 1e-12);
        assert_eq!(model.view_factors().factors()[0], vec![0.0, 0.25, 0.75]);
        assert_eq!(actual, model.solve(cx, &[400.0, 300.0], 280.0).unwrap());
    });
}

#[test]
fn g1_single_patch_cavity_equilibrium_and_adjacent_temperature_exchange() {
    let mesh = unit_mesh();
    with_cx(|cx| {
        let model = OpenDiffuseEnclosure::new(
            cx,
            vec![surface(&mesh, "wall", 1.0, 0.2)],
            vec![vec![0.75]],
            evidence(),
            ViewFactorTolerance::default(),
        )
        .unwrap();
        let equilibrium = model.solve(cx, &[350.0], 350.0).unwrap();
        assert_eq!(equilibrium.exchange.net_outward_heat_w, vec![0.0, 0.0]);
        let hot = f64::from_bits(350.0_f64.to_bits() + 1);
        let actual = model.solve(cx, &[hot], 350.0).unwrap();
        // First derivative is accurate to O(delta_T/T), far below this band.
        let expected = 0.125 * 4.0 * STEFAN_BOLTZMANN_W_M2_K4 * 350.0_f64.powi(3) * (hot - 350.0);
        near(actual.surface_heat_w()[0], expected, 5e-14);
        near(actual.ambient_heat_w(), -expected, 5e-14);
        assert!(actual.surface_heat_w()[0] > 0.0);
        // A cold space-like bath does not query the warm surface card domain.
        let cold_bath = model.solve(cx, &[350.0], 3.0).unwrap();
        near(
            cold_bath.surface_heat_w()[0],
            0.125 * STEFAN_BOLTZMANN_W_M2_K4 * (350.0_f64.powi(4) - 3.0_f64.powi(4)),
            2e-14,
        );
        assert!(model.solve(cx, &[600.0], 350.0).is_err());
    });
}

#[test]
fn g1_open_radiation_drives_the_real_conduction_solver_with_separate_external_heat() {
    let mesh = unit_mesh();
    let material = ConductivityModel::isotropic_declared(10.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    for h in [0.0, 5.0] {
        let condition = if h == 0.0 {
            ThermalBc::neumann(0.0).unwrap()
        } else {
            ThermalBc::robin(h, 300.0).unwrap()
        };
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .region(
                "hot",
                |f| on_box_face(f.centroid[0], 0.0),
                ThermalBc::dirichlet(400.0).unwrap(),
            )
            .unwrap()
            .region("wall", |f| on_box_face(f.centroid[0], 1.0), condition)
            .unwrap()
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &material,
            element_materials: None,
            source: &source,
        };
        with_cx(|cx| {
            let model = OpenDiffuseEnclosure::new(
                cx,
                vec![surface(&mesh, "wall", 1.0, 0.8)],
                vec![vec![0.25]],
                evidence(),
                ViewFactorTolerance::default(),
            )
            .unwrap();
            let solution = solve_with_open_gray_diffuse_enclosure(
                cx,
                problem,
                None,
                &model,
                280.0,
                config(),
                CoupledRadiationConfig::default(),
            )
            .unwrap();
            // Scalar conduction/convection/radiation balance, independent of FEM,
            // shifted loads, LU and the shared fixed-point iteration.
            let effective = 0.8 * 0.75 / (1.0 - 0.2 * 0.25);
            let mut low = 280.0_f64;
            let mut high = 400.0_f64;
            for _ in 0..80 {
                let t = 0.5 * (low + high);
                let q = effective * STEFAN_BOLTZMANN_W_M2_K4 * (t.powi(4) - 280.0_f64.powi(4));
                if 10.0 * (400.0 - t) - h * (t - 300.0) - q > 0.0 {
                    low = t;
                } else {
                    high = t;
                }
            }
            let expected = 0.5 * (low + high);
            let mean = model.surfaces()[0]
                .mean_temperature(&mesh, &solution.conduction.temperature)
                .unwrap();
            assert!((mean - expected).abs() < 2e-6, "{mean} versus {expected}");
            assert!((solution.convective_out_w - h * (mean - 300.0)).abs() < 1e-8);
            assert!((solution.radiative_out_w - solution.ambient_absorbed_w).abs() < 1e-10);
            let energy = &solution.conduction.report.energy;
            assert!(
                (energy.source_w + energy.dirichlet_in_w
                    - solution.prescribed_neumann_out_w
                    - solution.convective_out_w
                    - solution.radiative_out_w)
                    .abs()
                    < 2e-5
            );
            assert!(
                solution.radiation.equilibrium.residual_w
                    <= solution.radiation.equilibrium.residual_threshold_w
            );
            assert_eq!(solution.convective_robin_fluxes.len(), usize::from(h > 0.0));
            let exhausted = solve_with_open_gray_diffuse_enclosure(
                cx,
                problem,
                None,
                &model,
                280.0,
                config(),
                CoupledRadiationConfig {
                    max_iterations: 1,
                    ..CoupledRadiationConfig::default()
                },
            );
            assert!(exhausted.is_err());
        });
    }
}

#[test]
fn g0_g4_open_exchange_refuses_overfull_nonreciprocal_or_unexposed_systems() {
    let mesh = two_disconnected_slabs();
    for matrix in [
        vec![vec![0.0, 1.0], vec![1.0, 0.0]],
        vec![vec![0.0, 0.5], vec![0.25, 0.0]],
        vec![vec![0.6, 0.5], vec![0.5, 0.0]],
        vec![vec![0.0, f64::NAN], vec![0.0, 0.0]],
    ] {
        with_cx(|cx| {
            assert!(
                OpenDiffuseEnclosure::new(
                    cx,
                    vec![
                        surface(&mesh, "hot", 0.5, 0.8),
                        surface(&mesh, "cold", 1.5, 0.8)
                    ],
                    matrix,
                    evidence(),
                    ViewFactorTolerance::default()
                )
                .is_err()
            );
        });
    }
    let cancelled = with_cancelled_cx(|cx| {
        OpenDiffuseEnclosure::new(
            cx,
            vec![surface(&mesh, "hot", 0.5, 0.8)],
            vec![vec![0.0]],
            evidence(),
            ViewFactorTolerance::default(),
        )
    });
    assert!(matches!(cancelled, Err(ConductionError::Cancelled { .. })));
    with_cx(|cx| {
        let model = OpenDiffuseEnclosure::new(
            cx,
            vec![surface(&mesh, "hot", 0.5, 0.8)],
            vec![vec![0.0]],
            evidence(),
            ViewFactorTolerance::default(),
        )
        .unwrap();
        assert!(model.solve(cx, &[], 300.0).is_err());
        for bath in [0.0, -1.0, f64::NAN] {
            assert!(model.solve(cx, &[350.0], bath).is_err());
        }
        with_cancelled_cx(|cancelled| {
            assert!(matches!(
                model.solve(cancelled, &[350.0], 300.0),
                Err(ConductionError::Cancelled { .. })
            ));
        });
    });
}
