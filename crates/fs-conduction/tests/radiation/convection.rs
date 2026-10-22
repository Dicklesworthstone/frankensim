//! G1 physical heat balance and G0/G4 admission for coincident channels.
use super::*;
use fs_conduction::solve_with_convective_gray_diffuse_enclosure;

fn boundary(mesh: &ConductionMesh, hot: ThermalBc) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .region(
            "hot-reservoir",
            |f| on_box_face(f.centroid[0], 0.0),
            ThermalBc::dirichlet(400.0).unwrap(),
        )
        .unwrap()
        .region(
            "cold-reservoir",
            |f| on_box_face(f.centroid[0], 2.0),
            ThermalBc::dirichlet(300.0).unwrap(),
        )
        .unwrap()
        .region("hot-facing", |f| on_box_face(f.centroid[0], 0.5), hot)
        .unwrap()
        .region(
            "cold-facing",
            |f| on_box_face(f.centroid[0], 1.5),
            ThermalBc::robin(8.0, 320.0).unwrap(),
        )
        .unwrap()
        .adiabatic_remainder()
        .finish()
        .unwrap()
}

// Independent two-resistance thermal circuit. Each plane is one square metre;
// k/L=20 W/K, h_hot=5 and h_cold=8. The view-factor-one limiting model remains
// an authored algebra fixture, not a claim about finite separated squares.
fn reference() -> (f64, f64, f64) {
    let mut bracket = [0.0, 1000.0];
    for _ in 0..80 {
        let heat = 0.5 * (bracket[0] + bracket[1]);
        let hot: f64 = (20.0 * 400.0 + 5.0 * 300.0 - heat) / 25.0;
        let cold: f64 = (20.0 * 300.0 + 8.0 * 320.0 + heat) / 28.0;
        let radiative = STEFAN_BOLTZMANN_W_M2_K4 * (hot.powi(4) - cold.powi(4)) / 1.5;
        bracket[usize::from(heat > radiative)] = heat;
    }
    let heat = 0.5 * (bracket[0] + bracket[1]);
    ((9500.0 - heat) / 25.0, (8560.0 + heat) / 28.0, heat)
}

#[test]
fn g1_steady_convection_and_enclosure_exchange_close_independent_heat_balance() {
    let mesh = two_disconnected_slabs();
    let enclosure = two_slab_enclosure(&mesh);
    let boundary = boundary(&mesh, ThermalBc::robin(5.0, 300.0).unwrap());
    let original = boundary.clone();
    let material = ConductivityModel::isotropic_declared(10.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &material,
        element_materials: None,
        source: &source,
    };
    let solve = || {
        with_cx(|cx| {
            solve_with_convective_gray_diffuse_enclosure(
                cx,
                problem,
                None,
                &enclosure,
                config(),
                CoupledRadiationConfig::default(),
            )
            .unwrap()
        })
    };
    let solution = solve();
    assert_eq!(
        solution,
        solve(),
        "same admitted numerical profile replays exactly"
    );
    assert_eq!(
        boundary, original,
        "caller-owned physical boundary remains unchanged"
    );
    let (hot, cold, heat) = reference();
    let means = enclosure
        .surfaces()
        .iter()
        .map(|surface| {
            surface
                .mean_temperature(&mesh, &solution.conduction.temperature)
                .unwrap()
        })
        .collect::<Vec<_>>();
    assert!((means[0] - hot).abs() < 2e-6, "hot: {} vs {hot}", means[0]);
    assert!(
        (means[1] - cold).abs() < 2e-6,
        "cold: {} vs {cold}",
        means[1]
    );
    let convection = &solution.convective_robin_fluxes;
    assert_eq!(convection.len(), 2);
    assert_eq!(convection[0].region, "hot-facing");
    assert_eq!(convection[1].region, "cold-facing");
    assert!((convection[0].heat_rate_w - 5.0 * (hot - 300.0)).abs() < 1e-5);
    assert!((convection[1].heat_rate_w - 8.0 * (cold - 320.0)).abs() < 1e-5);
    assert_eq!(convection[0].mean_reference_temperature_k, 300.0);
    assert_eq!(convection[1].mean_reference_temperature_k, 320.0);
    let fresh = &solution.radiation.equilibrium.nonlinear_radiosity;
    assert!((fresh.net_outward_heat_w[0] - heat).abs() < 2e-5);
    assert!((fresh.net_outward_heat_w[1] + heat).abs() < 2e-5);
    assert!(
        (20.0 * (400.0 - means[0]) - convection[0].heat_rate_w - fresh.net_outward_heat_w[0]).abs()
            < 4e-5
    );
    assert!(
        (20.0 * (300.0 - means[1]) - convection[1].heat_rate_w - fresh.net_outward_heat_w[1]).abs()
            < 4e-5
    );
    for (combined, original) in solution
        .conduction
        .report
        .robin_fluxes
        .iter()
        .zip(convection)
    {
        let i = usize::from(combined.region == "cold-facing");
        assert!(
            (combined.heat_rate_w
                - original.heat_rate_w
                - solution.radiation.radiosity.net_outward_heat_w[i])
                .abs()
                < 1e-10,
            "radiation must not be counted as air heating"
        );
    }
    let energy = &solution.conduction.report.energy;
    assert!(
        (energy.source_w + energy.dirichlet_in_w
            - energy.neumann_out_w
            - solution.convective_out_w
            - fresh.enclosure_energy_closure_w)
            .abs()
            < 4e-5
    );
    assert!(
        solution.radiation.equilibrium.residual_w
            <= solution.radiation.equilibrium.residual_threshold_w
    );
    let reconstructed = with_cx(|cx| {
        fs_conduction::solve(
            cx,
            ConductionProblem {
                boundary: &solution.combined_boundary,
                ..problem
            },
            config(),
        )
        .unwrap()
    });
    assert!(
        reconstructed
            .temperature
            .iter()
            .zip(&solution.conduction.temperature)
            .all(|(a, b)| (a - b).abs() < 2e-6)
    );
}

#[test]
fn g0_enclosure_convection_requires_exact_trace_coverage() {
    let mesh = two_disconnected_slabs();
    let hot = RadiationSurface::new(
        &mesh,
        "hot-facing",
        |f| on_box_face(f.centroid[0], 0.5) && f.centroid[1] < 0.5,
        emissivity("hot-facing", 0.8, "hot"),
    )
    .unwrap();
    let cold = RadiationSurface::new(
        &mesh,
        "cold-facing",
        |f| on_box_face(f.centroid[0], 1.5) && f.centroid[1] < 0.5,
        emissivity("cold-facing", 0.8, "cold"),
    )
    .unwrap();
    let factors = ViewFactorMatrix::infinite_parallel_plates(hot.area_m2()).unwrap();
    let enclosure = GrayDiffuseEnclosure::new(vec![hot, cold], factors).unwrap();
    let boundary = boundary(&mesh, ThermalBc::robin(5.0, 300.0).unwrap());
    let material = ConductivityModel::isotropic_declared(10.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    let error = with_cx(|cx| {
        solve_with_convective_gray_diffuse_enclosure(
            cx,
            ConductionProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material,
                element_materials: None,
                source: &source,
            },
            None,
            &enclosure,
            config(),
            CoupledRadiationConfig::default(),
        )
    })
    .unwrap_err();
    assert!(matches!(error, ConductionError::Radiation { what, .. }
        if what.contains("exactly cover")));
}

#[test]
fn g0_g4_enclosure_convection_refuses_missing_law_cancellation_and_exhaustion() {
    let mesh = two_disconnected_slabs();
    let enclosure = two_slab_enclosure(&mesh);
    let material = ConductivityModel::isotropic_declared(10.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    for condition in [
        ThermalBc::adiabatic(),
        ThermalBc::Robin {
            htc: ScalarField::Nodal(vec![5.0; mesh.vertex_count()]),
            t_ref: ScalarField::Uniform(300.0),
        },
    ] {
        let boundary = boundary(&mesh, condition);
        let error = with_cx(|cx| {
            solve_with_convective_gray_diffuse_enclosure(
                cx,
                ConductionProblem {
                    mesh: &mesh,
                    boundary: &boundary,
                    material: &material,
                    element_materials: None,
                    source: &source,
                },
                None,
                &enclosure,
                config(),
                CoupledRadiationConfig::default(),
            )
        })
        .unwrap_err();
        assert!(matches!(error, ConductionError::Radiation { what, .. }
            if what.contains("uniform Robin")));
    }
    let boundary = boundary(&mesh, ThermalBc::robin(5.0, 300.0).unwrap());
    let original = boundary.clone();
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &material,
        element_materials: None,
        source: &source,
    };
    let cancelled = with_cancelled_cx(|cx| {
        solve_with_convective_gray_diffuse_enclosure(
            cx,
            problem,
            None,
            &enclosure,
            config(),
            CoupledRadiationConfig::default(),
        )
    });
    assert!(matches!(cancelled, Err(ConductionError::Cancelled { .. })));
    let exhausted = with_cx(|cx| {
        solve_with_convective_gray_diffuse_enclosure(
            cx,
            problem,
            None,
            &enclosure,
            config(),
            CoupledRadiationConfig {
                max_iterations: 1,
                ..CoupledRadiationConfig::default()
            },
        )
    });
    assert!(matches!(exhausted, Err(ConductionError::Radiation { .. })));
    assert_eq!(boundary, original);
}
