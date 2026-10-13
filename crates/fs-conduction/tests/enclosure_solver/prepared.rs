//! Compare the reusable operator with the unchanged one-shot elimination.
use super::*;

#[test]
fn prepared_reports_match_fresh_solves_for_independent_temperature_histories() {
    let (_mesh, left, right) = opposite_surfaces();
    let matrix = ViewFactorMatrix::infinite_parallel_plates(left.area_m2()).unwrap();
    let enclosure = GrayDiffuseEnclosure::new(vec![left, right], matrix).unwrap();
    let original = enclosure.clone();
    with_cx(|cx| {
        let prepared = enclosure.prepare_radiosity(cx).unwrap();
        // Reversed, isothermal, and repeated states exercise no hidden history.
        for temperatures in [[400.0, 300.0], [270.0, 450.0], [350.0, 350.0],
            [499.0, 251.0], [400.0, 300.0]] {
            assert_eq!(prepared.solve(cx, &temperatures).unwrap(),
                enclosure.solve(cx, &temperatures).unwrap());
        }
        assert!(prepared.solve(cx, &[300.0]).is_err());
        assert!(prepared.solve(cx, &[300.0, f64::NAN]).is_err());
        assert!(prepared.solve(cx, &[300.0, 501.0]).is_err());
        assert_eq!(enclosure, original);
    });
}

#[test]
fn pivoted_unequal_area_and_black_surface_operators_keep_exact_report_bits() {
    let (complex, positions) = box_grid([2, 2, 2], [1.0, 1.0, 1.0]);
    let mesh = ConductionMesh::new(complex, positions).unwrap();
    for eps in [[0.01, 0.6, 0.1], [1.0, 1.0, 1.0], [1.0, 0.02, 0.9]] {
        let first = RadiationSurface::new(&mesh, "first",
            |f| on_box_face(f.centroid[0], 0.0), emissivity("first", eps[0], "test")).unwrap();
        let second = RadiationSurface::new(&mesh, "second",
            |f| on_box_face(f.centroid[0], 1.0), emissivity("second", eps[1], "test")).unwrap();
        let third = RadiationSurface::new(&mesh, "half",
            |f| on_box_face(f.centroid[1], 0.0) && f.centroid[0] < 0.5,
            emissivity("half", eps[2], "test")).unwrap();
        let factors = vec![vec![0.9, 0.0, 0.1], vec![0.0, 0.7, 0.3], vec![0.2, 0.6, 0.2]];
        let matrix = ViewFactorMatrix::admit(
            vec![first.area_m2(), second.area_m2(), third.area_m2()], factors,
            ViewFactorEvidence::Analytic { geometry: "synthetic-reciprocal-algebra-only".into() },
            ViewFactorTolerance::default(),
        ).unwrap();
        let enclosure = GrayDiffuseEnclosure::new(vec![first, second, third], matrix).unwrap();
        with_cx(|cx| {
            let prepared = enclosure.prepare_radiosity(cx).unwrap();
            for temperatures in [[500.0, 350.0, 250.0], [290.0, 410.0, 320.0], [300.0; 3]] {
                assert_eq!(prepared.solve(cx, &temperatures).unwrap(),
                    enclosure.solve(cx, &temperatures).unwrap());
            }
        });
    }
}

#[test]
fn preparation_and_reuse_keep_cancellation_and_singular_limit_refusals() {
    let (_mesh, left, right) = opposite_surfaces();
    let matrix = ViewFactorMatrix::infinite_parallel_plates(left.area_m2()).unwrap();
    let enclosure = GrayDiffuseEnclosure::new(vec![left, right], matrix).unwrap();
    assert!(matches!(with_cancelled_cx(|cx| enclosure.prepare_radiosity(cx)),
        Err(ConductionError::Cancelled { .. })));
    let prepared = with_cx(|cx| enclosure.prepare_radiosity(cx)).unwrap();
    assert!(matches!(with_cancelled_cx(|cx| prepared.solve(cx, &[400.0, 300.0])),
        Err(ConductionError::Cancelled { .. })));
    // A cancelled application cannot poison the immutable reusable operator.
    with_cx(|cx| assert_eq!(prepared.solve(cx, &[400.0, 300.0]).unwrap(),
        enclosure.solve(cx, &[400.0, 300.0]).unwrap()));

    let (_mesh, left, right) = opposite_surfaces_with_emissivities(
        emissivity("left", f64::MIN_POSITIVE, "singular-limit"),
        emissivity("right", f64::MIN_POSITIVE, "singular-limit"),
    );
    let matrix = ViewFactorMatrix::infinite_parallel_plates(left.area_m2()).unwrap();
    let singular = GrayDiffuseEnclosure::new(vec![left, right], matrix).unwrap();
    with_cx(|cx| {
        assert!(matches!(singular.prepare_radiosity(cx), Err(ConductionError::Radiation { .. })));
        assert!(matches!(singular.solve(cx, &[400.0, 300.0]), Err(ConductionError::Radiation { .. })));
    });
}

/// The same physical acceptance, with one fresh radiosity factorization
/// per iteration. This comparison intentionally uses the existing public
/// producers, not a surrogate thermal model or a copied dense factorization.
fn fresh_factor_coupling(
    cx: &fs_exec::Cx<'_>, problem: ConductionProblem<'_>, enclosure: &GrayDiffuseEnclosure,
    cfg: SolveConfig, outer: CoupledRadiationConfig,
) -> fs_conduction::CoupledRadiationSolution {
    let means = |temperature: &[f64]| enclosure.surfaces().iter()
        .map(|surface| surface.mean_temperature(problem.mesh, temperature).unwrap())
        .collect::<Vec<_>>();
    let mut conduction = fs_conduction::solve(cx, problem, cfg.clone()).unwrap();
    let mut driving = means(&conduction.temperature);
    let mut history = Vec::new();
    let mut equilibrium_checks = 0usize;
    for iteration in 0..outer.max_iterations {
        let radiosity = enclosure.solve(cx, &driving).unwrap();
        let rows = enclosure.surfaces().iter().zip(&radiosity.net_outward_flux_w_m2)
            .map(|(surface, &q)| (surface.name().to_owned(), surface.face_slots().to_vec(), q))
            .collect::<Vec<_>>();
        let boundary = problem.boundary.with_uniform_outward_flux_overlays(problem.mesh, &rows).unwrap();
        let dofs = fs_conduction::DofMap::new(&boundary, problem.mesh.vertex_count()).unwrap();
        let mut next_config = cfg.clone();
        next_config.initial = InitialGuess::Free(dofs.gather(&conduction.temperature));
        conduction = fs_conduction::solve(cx, ConductionProblem { boundary: &boundary, ..problem },
            next_config).unwrap();
        assert!(conduction.report.final_residual <= conduction.report.residual_threshold);
        let actual = means(&conduction.temperature);
        let update = actual.iter().zip(&driving).map(|(a,b)| (a-b).abs()).fold(0.0f64, f64::max);
        let scale = actual.iter().copied().map(f64::abs).fold(0.0f64, f64::max);
        let threshold = outer.surface_temperature_rtol.mul_add(scale, outer.surface_temperature_atol_k);
        history.push(update);
        if update <= threshold {
            let equilibrium = super::equilibrium::actual_equilibrium(
                cx, problem, enclosure, &conduction, &radiosity);
            equilibrium_checks += 1;
            if equilibrium.residual_w <= equilibrium.residual_threshold_w {
                return fs_conduction::CoupledRadiationSolution {
                    conduction,
                    radiation: fs_conduction::CoupledRadiationReport { iterations: iteration+1,
                        surface_update_history_k: history, final_threshold_k: threshold, radiosity,
                        equilibrium, equilibrium_checks },
                };
            }
        }
        for (old, new) in driving.iter_mut().zip(actual) {
            *old = outer.relaxation.mul_add(new, (1.0-outer.relaxation)* *old);
        }
    }
    panic!("fresh-factor reference exhausted its original iteration cap");
}

#[test]
fn the_actual_enclosure_driver_preserves_the_complete_fresh_factor_trajectory() {
    let mesh = two_disconnected_slabs();
    let boundary = two_slab_boundary(&mesh);
    let enclosure = two_slab_enclosure(&mesh);
    let fallback = ConductivityModel::isotropic_declared(30.0).unwrap();
    let materials = fs_conduction::ElementMaterials::new(
        fs_conduction::MaterialTable::new([
            (fs_conduction::MaterialId(0), ConductivityModel::isotropic_declared(4.0).unwrap()),
            (fs_conduction::MaterialId(1), ConductivityModel::isotropic_declared(16.0).unwrap()),
        ]).unwrap(),
        mesh.complex().tets.iter().map(|tet| fs_conduction::MaterialId(
            if mesh.positions()[tet[0] as usize][0] < 1.0 { 0 } else { 1 })).collect(),
    ).unwrap();
    let source = ScalarField::Uniform(0.0);
    let problem = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &fallback,
        element_materials: Some(&materials), source: &source };
    let outer = CoupledRadiationConfig { surface_temperature_rtol: 0.0,
        surface_temperature_atol_k: 1e-8, relaxation: 0.25, max_iterations: 256 };
    with_cx(|cx| {
        let fresh = fresh_factor_coupling(cx, problem, &enclosure, config(), outer);
        let reused = solve_with_gray_diffuse_enclosure(cx, problem, None, &enclosure, config(), outer).unwrap();
        assert!(reused.radiation.iterations > 1);
        assert_eq!(reused, fresh);
    });
}
